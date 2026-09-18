#!/usr/bin/env python3
"""Assertions for the dependency-free core smoke test."""

import csv
import json
import sys
from pathlib import Path


root = Path(sys.argv[1])
annotation = json.loads((root / "01_annotation/vep_scores.summary.json").read_text())
assert annotation["variant_count"] == 6
assert annotation["variants_with_exact_vep_match"] == 6

with (root / "02_scoring/snp_filter_important_ld_tag_set_snp.tsv").open() as handle:
    winners = {row["window_id"]: row["best_snp_id"]
               for row in csv.DictReader(handle, delimiter="\t")}
assert winners == {"ECA1-0": "ECA1:200:T:C", "ECA1-400": "ECA1:450:C:A"}

with (root / "02_scoring/selection_ld_tag_set_candidates.tsv").open() as handle:
    candidates = list(csv.DictReader(handle, delimiter="\t"))
assert len(candidates) == 6
assert {row["window_id"] for row in candidates} == {"ECA1-0", "ECA1-400"}

fasta_records = sum(1 for line in (root / "03_probe_candidates/oligominer_candidate_flanks.fa").read_text().splitlines()
                    if line.startswith(">"))
assert fasta_records == 6

report = (root / "06_final_panel/summary_report.txt").read_text()
assert "final_strict_probes\t2" in report
assert "unresolved_windows\t0" in report
print("Core smoke test passed: 6 VEP matches, 2 window winners, 6 flanks, summary report")
