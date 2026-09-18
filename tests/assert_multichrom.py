#!/usr/bin/env python3
"""Confirm that two chromosomes were independently partitioned and gathered."""

import csv
import json
import sys
from pathlib import Path


root = Path(sys.argv[1])
for chrom in ("ECA1", "ECA2"):
    assert (root / f"00_inputs/by_chrom/{chrom}/windows.bed").is_file()
    assert (root / f"00_plink/by_chrom/{chrom}/variants.afreq").is_file()
    assert (root / f"00_plink/by_chrom/{chrom}/window_ld_tag_sets.tsv").is_file()

summary = json.loads((root / "00_plink/window_ld_tag_sets.summary.json").read_text())
assert summary["target_chromosomes"] == ["ECA1", "ECA2"]
assert summary["windows"] == 2

with (root / "02_scoring/selection_ld_tag_set_candidates.tsv").open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
assert {row["chrom"] for row in rows} == {"ECA1", "ECA2"}
print("Multi-chromosome scatter/gather test passed: ECA1 and ECA2")
