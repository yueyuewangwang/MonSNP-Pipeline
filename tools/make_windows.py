#!/usr/bin/env python3
"""Index a reference FASTA and generate BED design windows."""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="reference FASTA")
    parser.add_argument("--output", required=True, help="output BED4")
    parser.add_argument("--fai", default="", help="FAI path; defaults to <reference>.fai")
    parser.add_argument("--window-size", type=int, default=40000)
    parser.add_argument("--step", type=int, default=0,
                        help="window step; defaults to window size")
    parser.add_argument("--min-contig-length", type=int, default=1)
    parser.add_argument("--min-window-size", type=int, default=1,
                        help="omit a final partial window shorter than this")
    parser.add_argument("--include-regex", default="", help="keep matching contig names")
    parser.add_argument("--exclude-regex", default="", help="drop matching contig names")
    parser.add_argument("--exclude-contig", action="append", default=[])
    parser.add_argument("--summary", default="", help="JSON summary path")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def build_fai(reference: Path, fai: Path) -> list[tuple[str, int, int, int, int]]:
    entries: list[tuple[str, int, int, int, int]] = []
    seen: set[str] = set()
    with reference.open("rb") as handle:
        name = None
        length = offset = line_bases = line_width = 0
        previous_full_width: tuple[int, int] | None = None
        saw_short_line = False
        while True:
            line_offset = handle.tell()
            raw = handle.readline()
            if not raw:
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                break
            if raw.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                header = raw[1:].strip().split(None, 1)
                if not header:
                    raise ValueError(f"empty FASTA header at byte {line_offset}")
                name = header[0].decode("utf-8")
                if name in seen:
                    raise ValueError(f"duplicate FASTA sequence name: {name}")
                seen.add(name)
                length = line_bases = line_width = 0
                offset = handle.tell()
                previous_full_width = None
                saw_short_line = False
                continue
            if name is None:
                raise ValueError("FASTA sequence data appears before the first header")
            sequence = raw.rstrip(b"\r\n")
            if not sequence:
                raise ValueError(f"blank FASTA sequence line in {name}")
            if re.search(rb"[^A-Za-z*.-]", sequence):
                raise ValueError(f"invalid FASTA character in {name}")
            bases, width = len(sequence), len(raw)
            if not line_bases:
                line_bases, line_width = bases, width
                previous_full_width = (bases, width)
            elif saw_short_line:
                raise ValueError(f"sequence line follows a short final line in {name}")
            elif (bases, width) != previous_full_width:
                if bases > line_bases or width > line_width:
                    raise ValueError(f"inconsistent FASTA wrapping in {name}")
                saw_short_line = True
            length += bases
    if not entries:
        raise ValueError(f"no FASTA records found: {reference}")
    fai.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=fai.parent, delete=False,
                                     encoding="utf-8") as out:
        temp_name = out.name
        for row in entries:
            out.write("\t".join(map(str, row)) + "\n")
    os.replace(temp_name, fai)
    return entries


def read_fai(path: Path) -> list[tuple[str, int, int, int, int]]:
    entries = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"FAI line {line_no} has fewer than five columns")
            try:
                row = (fields[0], *map(int, fields[1:5]))
            except ValueError as exc:
                raise ValueError(f"invalid FAI numeric value at line {line_no}") from exc
            if (row[0] in seen or row[1] < 1 or row[2] < 0
                    or row[3] < 1 or row[4] < row[3]):
                raise ValueError(f"invalid or duplicate FAI sequence at line {line_no}: {row[0]}")
            seen.add(row[0])
            entries.append(row)
    if not entries:
        raise ValueError(f"no sequences in FAI: {path}")
    return entries


def main() -> None:
    cfg = parse_args()
    step = cfg.step or cfg.window_size
    if min(cfg.window_size, step, cfg.min_contig_length, cfg.min_window_size) < 1:
        raise ValueError("window, step, and minimum lengths must be positive")
    reference = Path(cfg.reference).expanduser().resolve()
    output = Path(cfg.output).expanduser().resolve()
    fai = Path(cfg.fai).expanduser().resolve() if cfg.fai else Path(str(reference) + ".fai")
    summary_path = (Path(cfg.summary).expanduser().resolve() if cfg.summary
                    else Path(str(output) + ".summary.json"))
    if not reference.is_file():
        raise FileNotFoundError(f"missing reference FASTA: {reference}")
    for path in (output, summary_path):
        if path.exists() and not cfg.force:
            raise FileExistsError(f"output exists; use --force to replace it: {path}")

    entries = read_fai(fai) if fai.exists() else build_fai(reference, fai)
    include = re.compile(cfg.include_regex) if cfg.include_regex else None
    exclude = re.compile(cfg.exclude_regex) if cfg.exclude_regex else None
    excluded_names = set(cfg.exclude_contig)
    kept_contigs = skipped_contigs = windows = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=output.parent, delete=False,
                                     encoding="utf-8") as out:
        temp_name = out.name
        for chrom, length, _offset, _bases, _width in entries:
            if (length < cfg.min_contig_length or chrom in excluded_names
                    or (include and include.search(chrom) is None)
                    or (exclude and exclude.search(chrom) is not None)):
                skipped_contigs += 1
                continue
            kept_contigs += 1
            for start0 in range(0, length, step):
                end0 = min(start0 + cfg.window_size, length)
                if end0 - start0 < cfg.min_window_size:
                    continue
                out.write(f"{chrom}\t{start0}\t{end0}\t{chrom}:{start0}-{end0}\n")
                windows += 1
    if not windows:
        Path(temp_name).unlink(missing_ok=True)
        raise ValueError("contig filters and minimum sizes produced zero windows")
    os.replace(temp_name, output)
    summary = {
        "status": "PASS", "reference": str(reference), "fai": str(fai),
        "output": str(output), "window_size": cfg.window_size, "step": step,
        "min_contig_length": cfg.min_contig_length,
        "min_window_size": cfg.min_window_size, "contigs_total": len(entries),
        "contigs_kept": kept_contigs, "contigs_skipped": skipped_contigs,
        "windows": windows, "coordinate_system": "0-based half-open BED",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
