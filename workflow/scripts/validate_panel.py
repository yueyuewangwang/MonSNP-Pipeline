#!/usr/bin/env python3
"""Validate final probe geometry, sequence, reference match, and uniqueness."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from prepare_oligominer_candidates import FastaRandomAccess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    cfg = parse_args()
    failures = Counter()
    rows = 0
    windows = set()
    sites = set()
    reader = FastaRandomAccess(cfg.reference, cfg.fai)
    try:
        with open(cfg.panel, encoding="utf-8", newline="") as handle:
            table = csv.DictReader(handle, delimiter="\t")
            required = {
                "window_id", "selected_snp_id", "chrom", "pos1",
                "probe_chrom", "probe_start", "probe_end", "probe_sequence",
                "probe_tier",
            }
            missing = required - set(table.fieldnames or [])
            if missing:
                raise ValueError(f"panel table missing columns: {sorted(missing)}")
            for row in table:
                rows += 1
                if row["window_id"] in windows:
                    failures["duplicate_window"] += 1
                windows.add(row["window_id"])
                if row["selected_snp_id"] in sites:
                    failures["duplicate_site"] += 1
                sites.add(row["selected_snp_id"])
                try:
                    pos1 = int(row["pos1"])
                    start1 = int(row["probe_start"])
                    end1 = int(row["probe_end"])
                except ValueError:
                    failures["invalid_coordinate"] += 1
                    continue
                sequence = row["probe_sequence"].upper()
                if row["probe_tier"] != "strict_oligominer":
                    failures["non_strict_probe"] += 1
                if row["chrom"] != row["probe_chrom"]:
                    failures["chromosome_mismatch"] += 1
                if not start1 <= pos1 <= end1:
                    failures["probe_does_not_cover_snp"] += 1
                if end1 - start1 + 1 != 120 or len(sequence) != 120:
                    failures["probe_length_not_120"] += 1
                if set(sequence) - set("ACGT"):
                    failures["non_acgt_sequence"] += 1
                try:
                    reference_sequence = reader.fetch(row["probe_chrom"], start1 - 1, end1)
                except KeyError:
                    failures["reference_chromosome_missing"] += 1
                    continue
                if sequence != reference_sequence:
                    failures["reference_sequence_mismatch"] += 1
    finally:
        reader.close()

    result = {
        "status": "PASS" if not failures else "FAIL",
        "panel_rows": rows,
        "unique_windows": len(windows),
        "unique_sites": len(sites),
        "failure_counts": dict(sorted(failures.items())),
    }
    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    Path(cfg.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
