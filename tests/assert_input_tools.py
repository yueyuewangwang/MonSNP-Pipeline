#!/usr/bin/env python3
"""Assertions for PLINK preparation, window generation, and preflight tests."""

import json
import sys
from pathlib import Path


root = Path(sys.argv[1])
prepared = root / "prepared/tiny"
with Path(str(prepared) + ".bim").open() as handle:
    rows = [line.split() for line in handle if line.strip()]
assert len(rows) == 6
assert rows[0][1] == "ECA1:100:T:G"
assert rows[0][4:] == ["G", "T"]
summary = json.loads(Path(str(prepared) + ".prepare.json").read_text())
assert summary["status"] == "PASS"
assert summary["output_biallelic_acgt_snps"] == 6

windows = [line.split("\t") for line in (root / "windows.bed").read_text().splitlines()]
assert [(row[1], row[2]) for row in windows] == [
    ("0", "400"), ("400", "800"), ("800", "1000")
]
assert (root / "generated_reference.fa.fai").is_file()

preflight = json.loads((root / "preflight.json").read_text())
assert preflight["status"] in {"PASS", "PASS_WITH_WARNINGS"}
assert preflight["errors"] == []
assert preflight["metrics"]["samples"] == 4
assert preflight["metrics"]["variants"] == 6
assert preflight["metrics"]["reference_sample_mismatches"] == 0
print("Input-tool tests passed: normalized PLINK, 3 windows, clean preflight")
