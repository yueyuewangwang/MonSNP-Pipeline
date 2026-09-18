#!/usr/bin/env python3
"""Assertions for the PLINK-entry and custom-region workflow test."""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from score_candidates import aggregate_weighted


root = Path(sys.argv[1])
with (root / "00_plink/window_ld_tag_sets.tsv").open() as handle:
    rows = [line.rstrip("\n").split("\t", 1) for line in handle]
assert [row[0] for row in rows] == ["ECA1-0", "ECA1-400"]
assert all(row[1] for row in rows)
tag_summary = json.loads(
    (root / "00_plink/window_ld_tag_sets.summary.json").read_text()
)
assert tag_summary["plink_version"] == "PLINK v1.90b7.7 mock"

matches = [("first", 5.0, []), ("second", 20.0, [])]
assert aggregate_weighted(matches, "max") == ("custom:second", 20.0)
assert aggregate_weighted(matches, "first") == ("custom:first", 5.0)


with (root / "02_scoring/snp_filter_snp_scores.tsv").open() as handle:
    scores = {row["snp_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
assert float(scores["ECA1:100:T:G"]["selection_bonus"]) == 5
assert float(scores["ECA1:150:C:T"]["selection_bonus"]) == 25
assert float(scores["ECA1:200:T:C"]["selection_bonus"]) == 20
assert scores["ECA1:150:C:T"]["selection_category"] == "custom:altitude+custom:priority"

with (root / "02_scoring/snp_filter_important_ld_tag_set_snp.tsv").open() as handle:
    winners = {row["window_id"]: row["best_snp_id"]
               for row in csv.DictReader(handle, delimiter="\t")}
assert winners["ECA1-0"] == "ECA1:150:C:T"

summary = json.loads(
    (root / "02_scoring/score_filter_chip_summary.json").read_text()
)
assert summary["scoring_profile"] == "custom"
assert summary["vep_scoring_enabled"] is False
assert summary["custom_region_policy"]["aggregation"] == "sum"
assert [item["name"] for item in summary["custom_region_policy"]["tracks"]] == [
    "altitude", "priority"
]
print("PLINK/custom smoke test passed")
