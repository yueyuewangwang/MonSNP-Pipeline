#!/usr/bin/env python
"""Merge target-covering OligoMiner results into the candidate table.

This file intentionally remains compatible with Python 2.7 because the
archived OligoMiner release used by the study requires Python 2.7.
"""

from __future__ import print_function

import argparse
import csv
import os


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--filtered-bed", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    cfg = parse_args()
    with open(cfg.candidates, "r") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fields = list(reader.fieldnames or [])
        required = {"candidate_id", "probe_status", "pos1", "flank_start0"}
        missing = required - set(fields)
        if missing:
            raise ValueError("candidate table missing columns: %s" % sorted(missing))
        rows = list(reader)

    target_local1 = {}
    for row in rows:
        try:
            target_local1[row["candidate_id"]] = (
                int(row["target_offset0"]) + 1 if row.get("target_offset0", "")
                else int(row["pos1"]) - int(row["flank_start0"])
            )
        except (KeyError, ValueError):
            target_local1[row["candidate_id"]] = None

    passing = {}
    with open(cfg.filtered_bed, "r") as handle:
        for line in handle:
            values = line.rstrip("\n\r").split("\t")
            if len(values) < 4:
                continue
            candidate_id = values[0]
            if candidate_id not in target_local1:
                continue
            try:
                start1, end1 = int(values[1]), int(values[2])
            except ValueError:
                continue
            target = target_local1[candidate_id]
            if target is not None and start1 <= target <= end1:
                passing[candidate_id] = passing.get(candidate_id, 0) + 1

    output_fields = fields + ["oligominer_pass", "oligominer_overlap_count"]
    status_counts = {}
    for row in rows:
        candidate_id = row["candidate_id"]
        if row["probe_status"] == "reference_allele_mismatch":
            status = "reference_allele_mismatch"
        elif passing.get(candidate_id, 0):
            status = "pass_oligominer"
        else:
            status = "fail_oligominer"
        row["probe_status"] = status
        row["oligominer_pass"] = "yes" if status == "pass_oligominer" else "no"
        row["oligominer_overlap_count"] = str(passing.get(candidate_id, 0))
        status_counts[status] = status_counts.get(status, 0) + 1

    output_dir = os.path.dirname(os.path.abspath(cfg.output))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    with open(cfg.output, "w") as output:
        writer = csv.DictWriter(output, fieldnames=output_fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with open(os.path.splitext(cfg.output)[0] + ".summary.tsv", "w") as output:
        output.write("status\tcount\n")
        for status in sorted(status_counts):
            output.write("%s\t%s\n" % (status, status_counts[status]))


if __name__ == "__main__":
    main()
