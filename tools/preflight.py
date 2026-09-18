#!/usr/bin/env python3
"""Validate workflow inputs before launching a study-scale run."""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

try:
    import yaml
except ImportError as exc:  # pragma: no cover - exercised by user environments
    raise SystemExit("PyYAML is required; activate the workflow Conda environment") from exc


CANONICAL_ID = re.compile(r"^([^:]+):(\d+):([ACGT]):([ACGT])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="workflow YAML configuration")
    parser.add_argument(
        "--base-dir", default=".",
        help="directory used to resolve relative paths, normally the repository root",
    )
    parser.add_argument("--report", default="", help="write a JSON report")
    parser.add_argument("--reference-sample", type=int, default=1000,
                        help="number of BIM variants checked against FASTA")
    return parser.parse_args()


class Checks:
    def __init__(self) -> None:
        self.error_count = 0
        self.warning_count = 0
        self.message_limit = 100
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.metrics: dict[str, object] = {}

    def error(self, message: str) -> None:
        self.error_count += 1
        if len(self.errors) < self.message_limit:
            self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warning_count += 1
        if len(self.warnings) < self.message_limit:
            self.warnings.append(message)


def resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def read_fai(path: Path, checks: Checks) -> dict[str, tuple[int, int, int, int]]:
    index: dict[str, tuple[int, int, int, int]] = {}
    if not path.is_file():
        checks.error(f"missing FASTA index: {path}; run tools/make_windows.py")
        return index
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                checks.error(f"FAI line {line_no} has fewer than five columns: {path}")
                continue
            try:
                values = tuple(map(int, fields[1:5]))
            except ValueError:
                checks.error(f"FAI line {line_no} has a nonnumeric field: {path}")
                continue
            length, offset, line_bases, line_width = values
            if fields[0] in index:
                checks.error(f"duplicate sequence in FAI: {fields[0]}")
            elif (length < 1 or offset < 0 or line_bases < 1
                    or line_width < line_bases):
                checks.error(f"invalid FAI values for {fields[0]}")
            else:
                index[fields[0]] = values
    checks.metrics["reference_sequences"] = len(index)
    return index


class FastaAccess:
    def __init__(self, fasta: Path, index: dict[str, tuple[int, int, int, int]]):
        self.handle = fasta.open("rb")
        self.index = index

    def base(self, chrom: str, pos1: int) -> str:
        length, offset, line_bases, line_width = self.index[chrom]
        if pos1 < 1 or pos1 > length:
            return ""
        pos0 = pos1 - 1
        line_no, column = divmod(pos0, line_bases)
        self.handle.seek(offset + line_no * line_width + column)
        return self.handle.read(1).decode("ascii", errors="replace").upper()

    def close(self) -> None:
        self.handle.close()


def scan_fam(path: Path, checks: Checks) -> int:
    if not path.is_file():
        checks.error(f"missing FAM: {path}")
        return 0
    samples = 0
    identifiers = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                checks.error(f"empty FAM record at line {line_no}: {path}")
                continue
            fields = line.split()
            if len(fields) < 6:
                checks.error(f"FAM line {line_no} has fewer than six columns")
                continue
            key = (fields[0], fields[1])
            if key in identifiers:
                checks.error(f"duplicate FID/IID in FAM at line {line_no}: {key}")
            identifiers.add(key)
            samples += 1
    if not samples:
        checks.error(f"FAM has no samples: {path}")
    checks.metrics["samples"] = samples
    return samples


def scan_bim(path: Path, checks: Checks, sample_size: int):
    if not path.is_file():
        checks.error(f"missing BIM: {path}")
        return 0, set(), []
    variants = 0
    chromosomes = set()
    reference_sample: list[tuple[str, int, str, str]] = []
    rng = random.Random(1729)
    previous_site = None
    nonadjacent_order = 0
    last_position_by_chrom: dict[str, int] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                checks.error(f"empty BIM record at line {line_no}: {path}")
                continue
            fields = line.split()
            if len(fields) != 6:
                checks.error(f"BIM line {line_no} does not have six columns")
                continue
            try:
                pos = int(fields[3])
            except ValueError:
                checks.error(f"BIM line {line_no} has a nonnumeric position")
                continue
            match = CANONICAL_ID.fullmatch(fields[1])
            if match is None:
                checks.error(
                    f"BIM ID is not CHROM:POS:REF:ALT at line {line_no}: {fields[1]}; "
                    "run tools/prepare_plink_input.py"
                )
                continue
            chrom, id_pos, ref, alt = match.group(1), int(match.group(2)), match.group(3), match.group(4)
            if chrom != fields[0] or id_pos != pos:
                checks.error(f"BIM ID coordinate mismatch at line {line_no}: {fields[1]}")
            if fields[4].upper() != alt or fields[5].upper() != ref:
                checks.error(f"BIM must have A1=ALT and A2=REF at line {line_no}: {fields[1]}")
            site = (chrom, pos)
            if site == previous_site:
                checks.error(f"multiple records at BIM coordinate {chrom}:{pos}")
            previous_site = site
            if chrom in last_position_by_chrom and pos < last_position_by_chrom[chrom]:
                nonadjacent_order += 1
            last_position_by_chrom[chrom] = pos
            variants += 1
            chromosomes.add(chrom)
            item = (chrom, pos, ref, fields[1])
            if len(reference_sample) < sample_size:
                reference_sample.append(item)
            elif sample_size:
                replacement = rng.randrange(variants)
                if replacement < sample_size:
                    reference_sample[replacement] = item
    if not variants:
        checks.error(f"BIM has no valid variants: {path}")
    if nonadjacent_order:
        checks.warn(
            f"BIM positions decrease {nonadjacent_order} times within chromosomes; "
            "sorting is recommended so duplicate coordinates remain detectable"
        )
    checks.metrics["variants"] = variants
    checks.metrics["variant_chromosomes"] = sorted(chromosomes)
    return variants, chromosomes, reference_sample


def check_bed(path: Path, variants: int, samples: int, checks: Checks) -> None:
    if not path.is_file():
        checks.error(f"missing BED: {path}")
        return
    with path.open("rb") as handle:
        magic = handle.read(3)
    if magic != b"\x6c\x1b\x01":
        checks.error(f"not a variant-major PLINK BED (bad magic bytes): {path}")
        return
    expected = 3 + math.ceil(samples / 4) * variants
    observed = path.stat().st_size
    if variants and samples and observed != expected:
        checks.error(f"BED size is {observed} bytes; expected {expected} from BIM/FAM counts")
    checks.metrics["bed_size_bytes"] = observed


def check_reference_sample(fasta: Path, index, rows, checks: Checks) -> None:
    if not fasta.is_file():
        checks.error(f"missing reference FASTA: {fasta}")
        return
    if not index:
        return
    mismatches = missing = 0
    examples = []
    reader = FastaAccess(fasta, index)
    try:
        for chrom, pos, expected, identifier in rows:
            if chrom not in index:
                missing += 1
                if len(examples) < 5:
                    examples.append(f"{identifier}: chromosome absent from FASTA")
                continue
            observed = reader.base(chrom, pos)
            if observed != expected:
                mismatches += 1
                if len(examples) < 5:
                    examples.append(f"{identifier}: FASTA={observed}, BIM_REF={expected}")
    finally:
        reader.close()
    if missing:
        checks.error(f"{missing} sampled BIM variants use chromosomes absent from FASTA")
    if mismatches:
        checks.error(f"{mismatches} sampled BIM REF alleles disagree with FASTA: {examples}")
    checks.metrics["reference_variants_sampled"] = len(rows)
    checks.metrics["reference_sample_mismatches"] = mismatches


def scan_bed_intervals(path: Path, label: str, index, checks: Checks,
                       require_ids: bool = False):
    count = 0
    chromosomes = set()
    identifiers = set()
    if not path.is_file():
        checks.error(f"missing {label}: {path}")
        return count, chromosomes
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                checks.error(f"{label} line {line_no} has fewer than three tab-separated columns")
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                checks.error(f"{label} line {line_no} has nonnumeric coordinates")
                continue
            if start < 0 or end <= start:
                checks.error(f"invalid 0-based half-open interval at {label} line {line_no}")
                continue
            chrom = fields[0]
            if chrom not in index:
                checks.error(f"{label} chromosome absent from FAI at line {line_no}: {chrom}")
            elif end > index[chrom][0]:
                checks.error(f"{label} interval exceeds {chrom} length at line {line_no}")
            if require_ids:
                identifier = fields[3].strip() if len(fields) >= 4 else f"{chrom}-{start}"
                if identifier in identifiers:
                    checks.error(f"duplicate window ID at {label} line {line_no}: {identifier}")
                identifiers.add(identifier)
            chromosomes.add(chrom)
            count += 1
    if not count:
        checks.error(f"{label} has no intervals: {path}")
    return count, chromosomes


def main() -> None:
    cfg = parse_args()
    if cfg.reference_sample < 0:
        raise ValueError("--reference-sample cannot be negative")
    base = Path(cfg.base_dir).expanduser().resolve()
    config_path = resolve(base, cfg.config)
    if not config_path.is_file():
        raise FileNotFoundError(f"missing config: {config_path}")
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    checks = Checks()
    inputs = config.get("inputs") or {}
    scoring = config.get("scoring") or {}
    parameters = config.get("parameters") or {}
    genotype_source = config.get("genotype_source")

    for key in ("reference", "reference_fai", "windows"):
        if not inputs.get(key):
            checks.error(f"config is missing inputs.{key}")
    if checks.error_count:
        reference = fai = windows = base / "__missing__"
        index = {}
    else:
        reference = resolve(base, inputs["reference"])
        fai = resolve(base, inputs["reference_fai"])
        windows = resolve(base, inputs["windows"])
        index = read_fai(fai, checks)

    variant_chromosomes = set()
    if genotype_source == "plink_bfile":
        if not inputs.get("bfile"):
            checks.error("config is missing inputs.bfile")
        else:
            prefix = resolve(base, inputs["bfile"])
            samples = scan_fam(Path(str(prefix) + ".fam"), checks)
            variants, variant_chromosomes, sample = scan_bim(
                Path(str(prefix) + ".bim"), checks, cfg.reference_sample
            )
            check_bed(Path(str(prefix) + ".bed"), variants, samples, checks)
            check_reference_sample(reference, index, sample, checks)
            checks.metrics["bfile_prefix"] = str(prefix)
    elif genotype_source == "precomputed":
        for key in ("afreq", "ld_tag_sets"):
            if not inputs.get(key) or not resolve(base, inputs[key]).is_file():
                checks.error(f"missing precomputed inputs.{key}")
        checks.warn("reference-allele checks are limited for precomputed genotype input")
    else:
        checks.error("genotype_source must be plink_bfile or precomputed")

    window_chromosomes = set()
    if index:
        window_count, window_chromosomes = scan_bed_intervals(
            windows, "design windows", index, checks, require_ids=True
        )
        checks.metrics["design_windows"] = window_count
        if variant_chromosomes:
            absent = sorted(window_chromosomes - variant_chromosomes)
            if absent:
                checks.warn(f"window chromosomes without BIM variants: {absent[:10]}")

    requested_chromosomes = config.get("target_chromosomes", "all")
    excluded_chromosomes = set(parameters.get("exclude_chromosomes") or [])
    if requested_chromosomes == "all":
        resolved_chromosomes = sorted(window_chromosomes - excluded_chromosomes)
    elif isinstance(requested_chromosomes, list) and requested_chromosomes:
        resolved_chromosomes = requested_chromosomes
        if len(resolved_chromosomes) != len(set(resolved_chromosomes)):
            checks.error("target_chromosomes contains duplicate values")
        missing = sorted(set(resolved_chromosomes) - window_chromosomes)
        if missing:
            checks.error(f"target_chromosomes absent from design windows: {missing[:10]}")
        blocked = sorted(set(resolved_chromosomes) & excluded_chromosomes)
        if blocked:
            checks.error(f"target_chromosomes also excluded by parameters: {blocked[:10]}")
    else:
        resolved_chromosomes = []
        checks.error('target_chromosomes must be "all" or a non-empty YAML list')
    if not resolved_chromosomes and window_chromosomes:
        checks.error("no target chromosomes remain after chromosome filtering")
    checks.metrics["target_chromosomes"] = resolved_chromosomes

    profile = scoring.get("profile")
    if profile not in {"study", "custom"}:
        checks.error("scoring.profile must be study or custom")
    aggregation = scoring.get("region_aggregation", "sum")
    if aggregation not in {"sum", "max", "first"}:
        checks.error("scoring.region_aggregation must be sum, max, or first")
    if not isinstance(scoring.get("include_maf", True), bool):
        checks.error("scoring.include_maf must be true or false")
    if not isinstance(scoring.get("include_vep", False), bool):
        checks.error("scoring.include_vep must be true or false")
    if not config.get("weights") or not resolve(base, config["weights"]).is_file():
        checks.error("missing scoring weight table configured by weights")
    if profile == "study":
        for key in ("vep", "conserved_200m", "conserved_100v", "qtl",
                    "literature", "expression", "local_adaptation",
                    "pur_local_adaptation", "domestication"):
            if not inputs.get(key) or not resolve(base, inputs[key]).is_file():
                checks.error(f"study profile requires inputs.{key}")

    tracks = scoring.get("weighted_regions") or []
    if not isinstance(tracks, list):
        checks.error("scoring.weighted_regions must be a YAML list")
        tracks = []
    names = set()
    track_metrics = []
    for position, track in enumerate(tracks, 1):
        if not isinstance(track, dict):
            checks.error(f"weighted_regions entry {position} is not a mapping")
            continue
        name = str(track.get("name", "")).strip()
        if not name:
            checks.error(f"weighted_regions entry {position} has no name")
        elif name in names:
            checks.error(f"duplicate weighted-region name: {name}")
        names.add(name)
        try:
            weight = float(track["weight"])
            if not math.isfinite(weight):
                raise ValueError("weight must be finite")
        except (KeyError, TypeError, ValueError):
            checks.error(f"weighted region {name or position} has a nonnumeric weight")
            weight = None
        if not track.get("path"):
            checks.error(f"weighted region {name or position} has no path")
            continue
        path = resolve(base, track["path"])
        interval_count, chromosomes = scan_bed_intervals(
            path, f"weighted region {name or position}", index, checks
        )
        if variant_chromosomes:
            absent = sorted(chromosomes - variant_chromosomes)
            if absent:
                checks.warn(f"weighted region {name} has chromosomes without BIM variants: {absent[:10]}")
        track_metrics.append({"name": name, "path": str(path), "weight": weight,
                              "intervals": interval_count})
    checks.metrics["weighted_regions"] = track_metrics
    if scoring.get("profile") == "custom" and not tracks:
        checks.warn("custom scoring has no weighted_regions; ranking will use only enabled MAF/VEP")
    if scoring.get("include_vep"):
        if not inputs.get("vep") or not resolve(base, inputs["vep"]).is_file():
            checks.error("scoring.include_vep is true but inputs.vep is missing")

    for key in ("plink_executable", "plink2_executable"):
        executable = str(parameters.get(key, "plink2" if key.startswith("plink2") else "plink"))
        located = shutil.which(executable) if os.sep not in executable else str(resolve(base, executable))
        if not located or not Path(located).is_file() or not os.access(located, os.X_OK):
            checks.warn(
                f"{key} is not currently executable ({executable}); Snakemake can provide it "
                "through --software-deployment-method conda"
            )
    if config.get("run_probe_design"):
        checks.warn("full probe design builds Bowtie2 and Jellyfish indexes and may require substantial RAM/disk")

    status = "FAIL" if checks.error_count else ("PASS_WITH_WARNINGS" if checks.warning_count else "PASS")
    report = {
        "status": status, "config": str(config_path), "base_dir": str(base),
        "error_count": checks.error_count, "warning_count": checks.warning_count,
        "errors": checks.errors, "warnings": checks.warnings,
        "metrics": checks.metrics,
    }
    report_path = resolve(base, cfg.report) if cfg.report else None
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Preflight: {status}")
    for message in checks.errors:
        print(f"ERROR: {message}")
    for message in checks.warnings:
        print(f"WARNING: {message}")
    print(json.dumps(checks.metrics, indent=2))
    if checks.error_count:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
