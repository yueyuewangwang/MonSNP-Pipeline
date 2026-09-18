#!/usr/bin/env python3
"""Build the study-compatible per-window PLINK LD tag set table.

PLINK ``--show-tags`` reports LD tag sets. This script preserves the
two-column table structure and records the exact generation parameters in its
summary.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


CANONICAL_ID = re.compile(r"^([^:]+):(\d+):([ACGT]+):([ACGT]+)$")


@dataclass(frozen=True)
class Window:
    index: int
    chrom: str
    start0: int
    end0: int
    window_id: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", required=True, help="PLINK BED/BIM/FAM prefix")
    parser.add_argument("--windows", required=True, help="BED3/BED4 window file")
    parser.add_argument("--plink", default="plink", help="PLINK 1.9 executable")
    parser.add_argument("--output", required=True, help="aggregated two-column LD tag set table")
    parser.add_argument("--summary", required=True, help="JSON run summary")
    parser.add_argument("--maf", type=float, default=0.05)
    parser.add_argument("--tag-kb", type=int, default=100)
    parser.add_argument("--tag-r2", type=float, default=0.8)
    parser.add_argument("--jobs", type=int, default=1, help="concurrent PLINK window processes")
    parser.add_argument("--memory-mb", type=int, default=2048, help="memory limit per PLINK process")
    parser.add_argument("--chr-set", type=int, default=31)
    parser.add_argument(
        "--window-boundary",
        choices=("study_compatible", "bed"),
        default="study_compatible",
        help=("study_compatible uses --from-bp START0 --to-bp END0 as in the "
              "original analysis; bed uses START0+1 through END0"),
    )
    return parser.parse_args()


def read_windows(path: str) -> list[Window]:
    windows: list[Window] = []
    seen: set[str] = set()
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"window row has fewer than three columns at {path}:{line_no}")
            try:
                start0, end0 = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"invalid window coordinates at {path}:{line_no}") from exc
            if start0 < 0 or end0 <= start0:
                raise ValueError(f"invalid BED interval at {path}:{line_no}")
            window_id = fields[3].strip() if len(fields) >= 4 and fields[3].strip() else f"{fields[0]}-{start0}"
            if window_id in seen:
                raise ValueError(f"duplicate window ID {window_id!r}")
            seen.add(window_id)
            windows.append(Window(len(windows), fields[0], start0, end0, window_id))
    if not windows:
        raise ValueError(f"no windows found in {path}")
    return windows


def validate_bfile(prefix: str) -> dict[str, int]:
    for suffix in (".bed", ".bim", ".fam"):
        path = prefix + suffix
        if not os.path.isfile(path):
            raise FileNotFoundError(f"missing PLINK binary input: {path}")

    variants = invalid_ids = allele_mismatches = 0
    examples: list[str] = []
    with open(prefix + ".bim", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split()
            if len(fields) < 6:
                raise ValueError(f"invalid BIM row at {prefix}.bim:{line_no}")
            variants += 1
            match = CANONICAL_ID.match(fields[1])
            if match is None or match.group(1) != fields[0] or int(match.group(2)) != int(fields[3]):
                invalid_ids += 1
                if len(examples) < 3:
                    examples.append(fields[1])
                continue
            # PLINK BED convention after --set-all-var-ids @:#:$r:$a stores ALT in A1
            # and REF in A2.  Exact orientation is required by downstream VEP matching.
            if match.group(3) != fields[5].upper() or match.group(4) != fields[4].upper():
                allele_mismatches += 1
                if len(examples) < 3:
                    examples.append(fields[1])
    if not variants:
        raise ValueError(f"no variants found in {prefix}.bim")
    if invalid_ids or allele_mismatches:
        raise ValueError(
            "BIM variants must use CHROM:POS:REF:ALT IDs with BIM A1=ALT and A2=REF; "
            f"invalid_ids={invalid_ids}, allele_mismatches={allele_mismatches}, examples={examples}"
        )
    with open(prefix + ".fam", encoding="utf-8", errors="replace") as handle:
        samples = sum(1 for line in handle if line.strip())
    if not samples:
        raise ValueError(f"no samples found in {prefix}.fam")
    return {"variants": variants, "samples": samples}


def parse_tags_list(path: str) -> list[str]:
    tag_sets: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter=" ", skipinitialspace=True)
        for fields in reader:
            fields = [value for value in fields if value]
            if not fields or fields[0] == "SNP":
                continue
            if len(fields) < 8:
                raise ValueError(f"unexpected PLINK --show-tags row in {path}: {fields}")
            if CANONICAL_ID.match(fields[0]) is None:
                raise ValueError(f"noncanonical target SNP ID in {path}: {fields[0]}")
            tags = fields[-1]
            if tags != "NONE":
                for tag in tags.split("|"):
                    if CANONICAL_ID.match(tag) is None:
                        raise ValueError(f"noncanonical tag SNP ID in {path}: {tag}")
            tag_sets.append(f"{fields[0]}|{tags}")
    return tag_sets


def main() -> None:
    cfg = parse_args()
    if not (0 <= cfg.maf <= 0.5):
        raise ValueError("--maf must be between 0 and 0.5")
    if not (0 <= cfg.tag_r2 <= 1):
        raise ValueError("--tag-r2 must be between 0 and 1")
    if cfg.jobs < 1 or cfg.memory_mb < 128:
        raise ValueError("--jobs must be positive and --memory-mb must be at least 128")
    executable = shutil.which(cfg.plink) if os.sep not in cfg.plink else cfg.plink
    if not executable or not os.path.isfile(executable):
        raise FileNotFoundError(f"PLINK 1.9 executable not found: {cfg.plink}")
    version_run = subprocess.run(
        [executable, "--version"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    version_text = version_run.stdout.strip()
    plink_version = version_text.splitlines()[0] if version_text else "unknown"

    bfile_stats = validate_bfile(cfg.bfile)
    windows = read_windows(cfg.windows)
    output = Path(cfg.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    empty_windows = 0
    ld_tag_sets_total = 0
    with tempfile.TemporaryDirectory(prefix="plink_tag_sets_", dir=str(output.parent)) as temp_dir:
        def run_window(window: Window):
            start1 = max(1, window.start0 if cfg.window_boundary == "study_compatible" else window.start0 + 1)
            prefix = os.path.join(temp_dir, f"window_{window.index:07d}")
            command = [
                executable, "--bfile", cfg.bfile, "--chr", window.chrom,
                "--from-bp", str(start1), "--to-bp", str(window.end0),
                "--maf", str(cfg.maf), "--show-tags", "all",
                "--tag-kb", str(cfg.tag_kb), "--tag-r2", str(cfg.tag_r2),
                "--memory", str(cfg.memory_mb), "--out", prefix,
                "--allow-extra-chr", "--chr-set", str(cfg.chr_set),
            ]
            process = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, check=False)
            tags_path = prefix + ".tags.list"
            if process.returncode != 0:
                if "All variants excluded" in process.stdout or "No variants remaining" in process.stdout:
                    return window, [], True
                tail = "\n".join(process.stdout.splitlines()[-20:])
                raise RuntimeError(f"PLINK failed for {window.window_id}:\n{tail}")
            tag_sets = parse_tags_list(tags_path) if os.path.isfile(tags_path) else []
            return window, tag_sets, not tag_sets

        with ThreadPoolExecutor(max_workers=cfg.jobs) as executor, output.open("w", encoding="utf-8") as handle:
            for window, tag_sets, empty in executor.map(run_window, windows):
                handle.write(f"{window.window_id}\t{','.join(tag_sets)}\n")
                empty_windows += int(empty)
                ld_tag_sets_total += len(tag_sets)

    summary = {
        "method": "PLINK 1.9 --show-tags LD tag sets",
        "bfile_prefix": os.path.abspath(cfg.bfile),
        "plink_executable": os.path.abspath(executable),
        "plink_version": plink_version,
        "windows": len(windows),
        "windows_without_ld_tag_sets": empty_windows,
        "ld_tag_sets_total": ld_tag_sets_total,
        "bfile_variants": bfile_stats["variants"],
        "bfile_samples": bfile_stats["samples"],
        "parameters": {
            "maf": cfg.maf,
            "tag_kb": cfg.tag_kb,
            "tag_r2": cfg.tag_r2,
            "window_boundary": cfg.window_boundary,
            "jobs": cfg.jobs,
            "memory_mb_per_job": cfg.memory_mb,
            "chr_set": cfg.chr_set,
        },
    }
    Path(cfg.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(cfg.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
