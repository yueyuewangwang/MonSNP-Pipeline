#!/usr/bin/env python3
"""Verify archived manuscript files against the frozen checksum table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checksums", required=True)
    parser.add_argument("--bundle-manifest", required=True,
                        help="TSV with item and relative_path columns")
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def table(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    cfg = parse_args()
    expected = {row["item"]: row for row in table(cfg.checksums)}
    mapping = {row["item"]: row["relative_path"] for row in table(cfg.bundle_manifest)}
    bundle = Path(cfg.bundle).resolve()
    results = []
    failures = 0
    for item, record in expected.items():
        relative = mapping.get(item, "")
        candidate = (bundle / relative).resolve() if relative else None
        if candidate is None or bundle not in candidate.parents or not candidate.is_file():
            results.append({"item": item, "status": "MISSING", "path": relative})
            failures += 1
            continue
        observed_size = candidate.stat().st_size
        observed_hash = sha256(candidate)
        size_ok = not record["size_bytes"] or observed_size == int(record["size_bytes"])
        hash_ok = observed_hash == record["sha256"]
        status = "PASS" if size_ok and hash_ok else "FAIL"
        failures += status != "PASS"
        results.append({
            "item": item, "status": status, "path": relative,
            "observed_size": observed_size, "observed_sha256": observed_hash,
        })
    report = {"status": "PASS" if failures == 0 else "FAIL",
              "checked": len(results), "failures": failures, "files": results}
    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    Path(cfg.output).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
