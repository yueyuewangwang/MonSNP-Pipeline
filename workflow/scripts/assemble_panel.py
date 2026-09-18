#!/usr/bin/env python3
"""Select the highest-ranked strict probe per window and write a final panel."""

from __future__ import annotations

import argparse
import csv
import json
from collections import OrderedDict
from pathlib import Path


PROBE_FIELDS = [
    "probe_chrom", "probe_start", "probe_end", "probe_sequence", "probe_tm",
    "center_distance_bp", "probe_tier",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--strict-designs", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--unresolved", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--summary-report", required=True)
    return parser.parse_args()


def read_table(path: str):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def canonical_site(value: str) -> str:
    fields = value.split(":")
    if len(fields) < 2:
        return ""
    try:
        return f"{fields[0]}:{int(fields[1])}"
    except ValueError:
        return ""


def write_table(path: str, fields, rows) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    cfg = parse_args()
    candidate_fields, candidates = read_table(cfg.candidates)
    required_candidates = {
        "window_id", "candidate_id", "candidate_rank", "snp_id", "chrom",
        "pos1", "ref", "alt", "base_score", "selection_bonus",
        "total_score", "selection_category",
    }
    missing = required_candidates - set(candidate_fields)
    if missing:
        raise ValueError(f"candidate table missing columns: {sorted(missing)}")

    design_fields, designs = read_table(cfg.strict_designs)
    required_designs = {"candidate_id", "probe_tier", *PROBE_FIELDS[:-1]}
    missing = required_designs - set(design_fields)
    if missing:
        raise ValueError(f"strict design table missing columns: {sorted(missing)}")
    designs_by_id = {}
    for row in designs:
        if row["probe_tier"] != "strict_oligominer":
            continue
        candidate_id = row["candidate_id"]
        if candidate_id in designs_by_id:
            raise ValueError(f"duplicate strict design: {candidate_id}")
        designs_by_id[candidate_id] = row

    by_window = OrderedDict()
    candidate_ids = set()
    for row in candidates:
        candidate_id = row["candidate_id"]
        if candidate_id in candidate_ids:
            raise ValueError(f"duplicate candidate ID: {candidate_id}")
        candidate_ids.add(candidate_id)
        by_window.setdefault(row["window_id"], []).append(row)
    unknown_designs = set(designs_by_id) - candidate_ids
    if unknown_designs:
        raise ValueError(f"strict designs absent from candidate table: {len(unknown_designs)}")

    panel_rows = []
    unresolved_rows = []
    used_sites = set()
    duplicate_site_skips = 0
    for window_id, rows in by_window.items():
        rows.sort(key=lambda row: (int(row["candidate_rank"]), row["snp_id"]))
        strict_count = 0
        chosen = None
        chosen_design = None
        for row in rows:
            design = designs_by_id.get(row["candidate_id"])
            if design is None:
                continue
            strict_count += 1
            site = canonical_site(row["snp_id"])
            if site in used_sites:
                duplicate_site_skips += 1
                continue
            chosen, chosen_design = row, design
            used_sites.add(site)
            break
        if chosen is None:
            reason = ("no_strict_probe" if strict_count == 0
                      else "strict_sites_already_selected_in_other_windows")
            unresolved_rows.append({
                "window_id": window_id,
                "candidate_count": str(len(rows)),
                "strict_candidate_count": str(strict_count),
                "reason": reason,
            })
            continue
        output = dict(chosen)
        output.update({field: chosen_design.get(field, "") for field in PROBE_FIELDS})
        output["selected_snp_id"] = canonical_site(chosen["snp_id"])
        output["selection_method"] = "strict_oligominer"
        panel_rows.append(output)

    panel_fields = candidate_fields + [
        "selected_snp_id", "selection_method",
        *[field for field in PROBE_FIELDS if field not in candidate_fields],
    ]
    unresolved_fields = ["window_id", "candidate_count", "strict_candidate_count", "reason"]
    write_table(cfg.output, panel_fields, panel_rows)
    write_table(cfg.unresolved, unresolved_fields, unresolved_rows)

    Path(cfg.bed).parent.mkdir(parents=True, exist_ok=True)
    with open(cfg.bed, "w", encoding="utf-8") as handle:
        for row in panel_rows:
            handle.write(
                f"{row['chrom']}\t{int(row['pos1']) - 1}\t{row['pos1']}\t"
                f"{row['selected_snp_id']}\t{row['window_id']}\n"
            )

    summary = {
        "windows_with_candidates": len(by_window),
        "selected_strict_probes": len(panel_rows),
        "unresolved_windows": len(unresolved_rows),
        "duplicate_cross_window_site_skips": duplicate_site_skips,
        "selected_windows_unique": len({row["window_id"] for row in panel_rows}),
        "selected_sites_unique": len(used_sites),
        "policy": "lowest candidate_rank with a strict OligoMiner probe; sites remain globally unique",
    }
    Path(cfg.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(cfg.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    report_lines = [
        "MonSNP-Pipeline summary report",
        "stage\tcount",
        f"ranked_candidate_variants\t{len(candidates)}",
        f"windows_with_ranked_candidates\t{len(by_window)}",
        f"candidates_with_strict_oligominer_designs\t{len(designs_by_id)}",
        f"final_strict_probes\t{len(panel_rows)}",
        f"unresolved_windows\t{len(unresolved_rows)}",
        f"duplicate_cross_window_site_skips\t{duplicate_site_skips}",
    ]
    Path(cfg.summary_report).write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
