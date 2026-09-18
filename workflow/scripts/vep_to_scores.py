#!/usr/bin/env python3
"""Convert VEP consequences to one maximum functional score per variant."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


SNP_RE = re.compile(r"^([^:]+):(\d+):([^:]+):([^:]+)$")
WEIGHT_KEYS = {
    "missense_variant": "vep_missense",
    "splice_site_variant": "vep_splice",
    "UTR_variant": "vep_utr",
    "upstream_downstream_variant": "vep_gene_end",
    "synonymous_variant": "vep_synonymous",
    "intronic_variant": "vep_intronic",
    "intergenic_variant": "vep_intergenic",
    "unknown_other": "vep_unknown",
}


def open_text(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode.replace("t", ""), encoding="utf-8", errors="replace")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants", required=True, help="PLINK2 .afreq file")
    parser.add_argument("--vep", required=True, help="VEP tabular output, optionally gzip-compressed")
    parser.add_argument("--weights", required=True, help="key/value scoring table")
    parser.add_argument("--output", required=True, help="output score track (.tsv or .tsv.gz)")
    parser.add_argument("--summary", required=True, help="JSON run summary")
    return parser.parse_args()


def read_weights(path: str) -> dict[str, float]:
    values: dict[str, float] = {}
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2 or fields[0].strip().lower() == "key":
                continue
            try:
                values[fields[0].strip()] = float(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid weight at {path}:{line_no}") from exc
    missing = set(WEIGHT_KEYS.values()) - set(values)
    if missing:
        raise ValueError(f"weight table is missing VEP keys: {sorted(missing)}")
    return {category: values[key] for category, key in WEIGHT_KEYS.items()}


def parse_snp_id(value: str):
    match = SNP_RE.match(value)
    if not match:
        return None
    return match.group(1), int(match.group(2)), match.group(3).upper(), match.group(4).upper()


def read_variants(path: str):
    variants: dict[str, tuple[str, int, str, str]] = {}
    by_alt: dict[tuple[str, int, str], list[str]] = defaultdict(list)
    occupied_sites: dict[tuple[str, int], str] = {}
    invalid = duplicate = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 4:
                invalid += 1
                continue
            parsed = parse_snp_id(fields[1])
            if parsed is None:
                invalid += 1
                continue
            chrom, pos1, ref, alt = parsed
            if "," in alt:
                invalid += 1
                continue
            if fields[0] != chrom or fields[2].upper() != ref or fields[3].upper() != alt:
                raise ValueError(f"AFREQ columns disagree with ID: {fields[1]}")
            if fields[1] in variants:
                duplicate += 1
                continue
            site = (chrom, pos1)
            if site in occupied_sites:
                raise ValueError(
                    f"multiple ALT records at {chrom}:{pos1}; normalize to one "
                    "biallelic record per coordinate before running this workflow"
                )
            occupied_sites[site] = fields[1]
            variants[fields[1]] = parsed
            by_alt[(chrom, pos1, alt)].append(fields[1])
    if not variants:
        raise ValueError(f"no canonical biallelic variants found in {path}")
    return variants, by_alt, invalid, duplicate


def classify_term(term: str) -> str:
    if term == "missense_variant":
        return "missense_variant"
    if "splice" in term:
        return "splice_site_variant"
    if term in {"5_prime_UTR_variant", "3_prime_UTR_variant"}:
        return "UTR_variant"
    if term in {"upstream_gene_variant", "downstream_gene_variant"}:
        return "upstream_downstream_variant"
    if term == "synonymous_variant":
        return "synonymous_variant"
    if term == "intron_variant":
        return "intronic_variant"
    if term == "intergenic_variant":
        return "intergenic_variant"
    return "unknown_other"


def classify_consequence(text: str, weights: dict[str, float]):
    categories = {classify_term(term) for term in text.split("&") if term and term != "-"}
    if not categories:
        categories = {"unknown_other"}
    maximum = max(weights[category] for category in categories)
    return maximum, {category for category in categories if weights[category] == maximum}


def uploaded_to_snp_id(uploaded: str) -> str:
    parsed = parse_snp_id(uploaded)
    if parsed is not None:
        return uploaded
    parts = uploaded.rsplit("_", 2)
    if len(parts) != 3 or "/" not in parts[2]:
        return ""
    ref, alt = parts[2].split("/", 1)
    try:
        pos1 = int(parts[1])
    except ValueError:
        return ""
    return f"{parts[0]}:{pos1}:{ref.upper()}:{alt.upper()}"


def location_key(location: str, allele: str):
    match = re.match(r"^([^:]+):(\d+)(?:-\d+)?$", location)
    if not match:
        return None
    return match.group(1), int(match.group(2)), allele.upper()


def stream_vep(path: str, variants, by_alt, weights):
    best = {snp_id: weights["unknown_other"] for snp_id in variants}
    categories = {snp_id: {"unknown_other"} for snp_id in variants}
    accepted = Counter()
    counts = Counter()
    term_counts = Counter()
    with open_text(path) as handle:
        for line in handle:
            counts["lines_scanned"] += 1
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 7:
                counts["malformed_rows"] += 1
                continue
            snp_id = uploaded_to_snp_id(fields[0])
            candidate_ids = [snp_id] if snp_id in variants else []
            if not candidate_ids:
                key = location_key(fields[1], fields[2])
                candidate_ids = by_alt.get(key, []) if key else []
            if not candidate_ids:
                counts["rows_outside_variant_universe"] += 1
                continue
            score, row_categories = classify_consequence(fields[6], weights)
            for term in fields[6].split("&"):
                if term and term != "-":
                    term_counts[term] += 1
            for candidate_id in candidate_ids:
                expected_alt = variants[candidate_id][3]
                if fields[2].upper() != expected_alt:
                    counts["allele_mismatch_rows"] += 1
                    continue
                if not accepted[candidate_id] or score > best[candidate_id]:
                    best[candidate_id] = score
                    categories[candidate_id] = set(row_categories)
                elif score == best[candidate_id]:
                    categories[candidate_id].update(row_categories)
                accepted[candidate_id] += 1
                counts["accepted_variant_rows"] += 1
    return best, categories, accepted, counts, term_counts


def main() -> None:
    cfg = parse_args()
    weights = read_weights(cfg.weights)
    variants, by_alt, invalid, duplicate = read_variants(cfg.variants)
    best, categories, accepted, counts, term_counts = stream_vep(
        cfg.vep, variants, by_alt, weights)

    Path(cfg.output).parent.mkdir(parents=True, exist_ok=True)
    with open_text(cfg.output, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#chrom", "pos1", "vep_score", "vep_category", "accepted_vep_rows", "snp_id"])
        for snp_id, (chrom, pos1, _ref, _alt) in variants.items():
            writer.writerow([
                chrom, pos1, f"{best[snp_id]:.8f}",
                ",".join(sorted(categories[snp_id])), accepted[snp_id], snp_id,
            ])

    summary = {
        "variant_count": len(variants),
        "variants_with_exact_vep_match": sum(value > 0 for value in accepted.values()),
        "variants_without_vep_match": sum(value == 0 for value in accepted.values()),
        "invalid_afreq_rows": invalid,
        "duplicate_afreq_ids": duplicate,
        "vep_scan_counts": dict(sorted(counts.items())),
        "vep_consequence_term_counts": dict(sorted(term_counts.items())),
        "unannotated_policy": "assign vep_unknown weight",
    }
    Path(cfg.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(cfg.summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
