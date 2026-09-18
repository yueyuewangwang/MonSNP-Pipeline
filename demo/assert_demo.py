#!/usr/bin/env python3
"""Confirm that the bundled custom-data demo produced the expected results."""

import csv
import json
from pathlib import Path


root = Path(__file__).resolve().parent
work = root / "work"
prepared = json.loads((work / "prepared/demo.prepare.json").read_text())
assert prepared["status"] == "PASS"
assert prepared["output_biallelic_acgt_snps"] == 6
assert prepared["variants_with_renamed_chromosome"] == 6

preflight = json.loads((work / "reports/preflight.json").read_text())
assert preflight["status"] == "PASS"
assert preflight["error_count"] == 0

score_path = work / "results/02_scoring/snp_filter_snp_scores.tsv"
with score_path.open() as handle:
    scores = {row["snp_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
assert len(scores) == 6
assert float(scores["ECA1:100:T:G"]["selection_bonus"]) == 25
assert float(scores["ECA1:150:C:T"]["selection_bonus"]) == 37
assert float(scores["ECA1:450:C:A"]["selection_bonus"]) == 12
assert float(scores["ECA1:600:T:C"]["selection_bonus"]) == 0
print("Demo assertions passed: 6 SNPs, chromosome map, preflight, and +37 overlap")
