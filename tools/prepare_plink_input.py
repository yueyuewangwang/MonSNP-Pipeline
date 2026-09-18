#!/usr/bin/env python3
"""Normalize a PLINK BED/BIM/FAM set for the LD tag set workflow."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


CANONICAL_ID = re.compile(r"^([^:]+):(\d+):([ACGT]):([ACGT])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use PLINK2 and a reference FASTA to retain biallelic A/C/G/T SNPs, "
            "set REF/ALT orientation, and assign CHROM:POS:REF:ALT IDs."
        )
    )
    parser.add_argument("--bfile", required=True, help="input PLINK prefix")
    parser.add_argument("--reference", required=True, help="matching reference FASTA")
    parser.add_argument("--out", required=True, help="output PLINK prefix")
    parser.add_argument("--plink2", default="plink2", help="PLINK2 executable")
    parser.add_argument("--chr-set", type=int, default=None,
                        help="autosome count passed to PLINK2, e.g. 22 or 31")
    parser.add_argument("--chrom-map", default="",
                        help="optional two-column OLD_CHROM NEW_CHROM mapping")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--force", action="store_true", help="replace existing output files")
    return parser.parse_args()


def reject_blank_records(path: Path, label: str) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(
                    f"{label} contains an empty record at line {line_no}: {path}. "
                    "Remove blank lines before running PLINK2."
                )
            count += 1
    if not count:
        raise ValueError(f"{label} has no records: {path}")
    return count


def read_chrom_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    targets: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(
                    f"chromosome map line {line_no} must contain OLD_CHROM NEW_CHROM"
                )
            old, new = fields
            if old in mapping:
                raise ValueError(f"duplicate source chromosome in map: {old}")
            if new in targets:
                raise ValueError(f"multiple source chromosomes map to the same target: {new}")
            mapping[old] = new
            targets.add(new)
    if not mapping:
        raise ValueError(f"chromosome map has no entries: {path}")
    return mapping


def write_remapped_source(source: Path, destination: Path,
                          mapping: dict[str, str]) -> tuple[Path, int]:
    for suffix in (".bed", ".fam"):
        shutil.copyfile(Path(str(source) + suffix), Path(str(destination) + suffix))
    remapped = 0
    observed = set()
    with Path(str(source) + ".bim").open(encoding="utf-8", errors="replace") as source_bim, \
            Path(str(destination) + ".bim").open("w", encoding="utf-8") as output_bim:
        for line_no, line in enumerate(source_bim, 1):
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"BIM line {line_no} does not have six columns")
            observed.add(fields[0])
            replacement = mapping.get(fields[0], fields[0])
            remapped += replacement != fields[0]
            fields[0] = replacement
            output_bim.write("\t".join(fields) + "\n")
    unused = sorted(set(mapping) - observed)
    if unused:
        raise ValueError(f"chromosome map source names absent from BIM: {unused[:10]}")
    return destination, remapped


def validate_output(prefix: Path) -> dict[str, int]:
    bed = Path(str(prefix) + ".bed")
    bim = Path(str(prefix) + ".bim")
    fam = Path(str(prefix) + ".fam")
    if bed.read_bytes()[:3] != b"\x6c\x1b\x01":
        raise ValueError(f"output is not a variant-major PLINK BED: {bed}")

    samples = reject_blank_records(fam, "FAM")
    variants = 0
    seen_sites: set[tuple[str, int]] = set()
    with bim.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f"BIM contains an empty record at line {line_no}: {bim}")
            fields = line.split()
            if len(fields) != 6:
                raise ValueError(f"BIM line {line_no} does not have six columns: {bim}")
            match = CANONICAL_ID.fullmatch(fields[1])
            if match is None:
                raise ValueError(f"noncanonical output variant ID at {bim}:{line_no}: {fields[1]}")
            chrom, pos = match.group(1), int(match.group(2))
            ref, alt = match.group(3), match.group(4)
            if chrom != fields[0] or pos != int(fields[3]):
                raise ValueError(f"ID coordinate does not match BIM columns at {bim}:{line_no}")
            if fields[4].upper() != alt or fields[5].upper() != ref:
                raise ValueError(f"expected BIM A1=ALT and A2=REF at {bim}:{line_no}")
            site = (chrom, pos)
            if site in seen_sites:
                raise ValueError(f"multiple variants remain at {chrom}:{pos}")
            seen_sites.add(site)
            variants += 1
    if not variants:
        raise ValueError("no biallelic A/C/G/T SNPs remain after normalization")

    expected_size = 3 + ((samples + 3) // 4) * variants
    actual_size = bed.stat().st_size
    if actual_size != expected_size:
        raise ValueError(
            f"unexpected BED size: observed {actual_size}, expected {expected_size} "
            f"for {samples} samples and {variants} variants"
        )
    return {"samples": samples, "variants": variants}


def run(command: list[str], label: str) -> str:
    process = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    if process.returncode:
        tail = "\n".join(process.stdout.splitlines()[-40:])
        raise RuntimeError(f"PLINK2 {label} failed (exit {process.returncode}):\n{tail}")
    return process.stdout


def main() -> None:
    cfg = parse_args()
    if cfg.threads < 1:
        raise ValueError("--threads must be positive")
    if cfg.chr_set is not None and cfg.chr_set < 1:
        raise ValueError("--chr-set must be positive")

    source = Path(cfg.bfile).expanduser().resolve()
    reference = Path(cfg.reference).expanduser().resolve()
    destination = Path(cfg.out).expanduser().resolve()
    chrom_map_path = (Path(cfg.chrom_map).expanduser().resolve()
                      if cfg.chrom_map else None)
    if source == destination:
        raise ValueError("--out must differ from --bfile; never normalize in place")
    for suffix in (".bed", ".bim", ".fam"):
        path = Path(str(source) + suffix)
        if not path.is_file():
            raise FileNotFoundError(f"missing PLINK input: {path}")
    if not reference.is_file():
        raise FileNotFoundError(f"missing reference FASTA: {reference}")
    if chrom_map_path and not chrom_map_path.is_file():
        raise FileNotFoundError(f"missing chromosome map: {chrom_map_path}")

    source_variants = reject_blank_records(Path(str(source) + ".bim"), "BIM")
    source_samples = reject_blank_records(Path(str(source) + ".fam"), "FAM")
    chrom_map = read_chrom_map(chrom_map_path) if chrom_map_path else {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    final_paths = [Path(str(destination) + suffix)
                   for suffix in (".bed", ".bim", ".fam", ".log", ".prepare.json")]
    existing = [str(path) for path in final_paths if path.exists()]
    if existing and not cfg.force:
        raise FileExistsError("output exists; use --force to replace it: " + ", ".join(existing))

    executable = shutil.which(cfg.plink2) if os.sep not in cfg.plink2 else cfg.plink2
    if not executable or not Path(executable).is_file():
        raise FileNotFoundError(f"PLINK2 executable not found: {cfg.plink2}")

    with tempfile.TemporaryDirectory(prefix="prepare_plink_", dir=destination.parent) as temp_dir:
        temp_dir_path = Path(temp_dir)
        plink_source = source
        remapped_variants = 0
        if chrom_map:
            plink_source, remapped_variants = write_remapped_source(
                source, temp_dir_path / "chrom_remapped", chrom_map
            )

        pgen_prefix = temp_dir_path / "normalized_pgen"
        temp_prefix = temp_dir_path / "normalized"
        normalize_command = [
            str(executable), "--bfile", str(plink_source), "--fa", str(reference),
            "--ref-from-fa", "force", "--snps-only", "just-acgt",
            "--max-alleles", "2", "--set-all-var-ids", "@:#:$r:$a",
            "--new-id-max-allele-len", "1000", "--sort-vars", "natural",
            "--make-pgen", "--allow-extra-chr", "--threads", str(cfg.threads),
            "--out", str(pgen_prefix),
        ]
        if cfg.chr_set is not None:
            insertion = normalize_command.index("--allow-extra-chr")
            normalize_command[insertion:insertion] = ["--chr-set", str(cfg.chr_set)]
        run(normalize_command, "reference normalization")

        bed_command = [
            str(executable), "--pfile", str(pgen_prefix), "--make-bed",
            "--allow-extra-chr", "--threads", str(cfg.threads),
            "--out", str(temp_prefix),
        ]
        if cfg.chr_set is not None:
            insertion = bed_command.index("--allow-extra-chr")
            bed_command[insertion:insertion] = ["--chr-set", str(cfg.chr_set)]
        run(bed_command, "BED conversion")
        stats = validate_output(temp_prefix)

        combined_log = Path(str(temp_prefix) + ".log")
        combined_log.write_text(
            "### Reference normalization and sorting ###\n"
            + Path(str(pgen_prefix) + ".log").read_text(encoding="utf-8", errors="replace")
            + "\n### BED conversion ###\n"
            + combined_log.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )
        version_run = subprocess.run(
            [str(executable), "--version"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False,
        )
        summary = {
            "status": "PASS",
            "method": "PLINK2 reference-orientation and canonical-ID normalization",
            "input_prefix": str(source), "reference": str(reference),
            "output_prefix": str(destination), "source_samples": source_samples,
            "source_variants": source_variants, "output_samples": stats["samples"],
            "output_biallelic_acgt_snps": stats["variants"],
            "filtered_variants": source_variants - stats["variants"],
            "chromosome_map": str(chrom_map_path) if chrom_map_path else None,
            "chromosome_map_entries": len(chrom_map),
            "variants_with_renamed_chromosome": remapped_variants,
            "plink2": str(Path(executable).resolve()),
            "plink2_version": version_run.stdout.strip().splitlines()[0]
            if version_run.stdout.strip() else "unknown",
            "chr_set": cfg.chr_set, "threads": cfg.threads,
            "commands": [normalize_command, bed_command],
        }
        Path(str(temp_prefix) + ".prepare.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        for suffix in (".bed", ".bim", ".fam", ".log", ".prepare.json"):
            source_path = Path(str(temp_prefix) + suffix)
            if not source_path.exists():
                raise FileNotFoundError(f"PLINK2 did not create expected output: {source_path}")
            os.replace(source_path, Path(str(destination) + suffix))

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
