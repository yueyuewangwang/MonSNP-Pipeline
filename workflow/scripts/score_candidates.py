#!/usr/bin/env python3
"""Score variants and rank probe candidates within window LD tag sets.

The PLINK2 ``.afreq`` IDs define the variant universe. ``study`` mode
reproduces the original annotation and selection-region policy. ``custom`` mode
optionally combines MAF and VEP with user-supplied weighted BED tracks. Large
LD tag set inputs are streamed one line at a time.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


SNP_RE = re.compile(r"^([^:]+):(\d+)(?::([^:]+):([^:]+))?$")
COMPONENTS = ("maf_score", "annotation_score", "con200m_score", "con100v_score",
              "qtl_score", "literature_score", "expression_breadth",
              "expression_intensity")

MAF_WEIGHTS = {
    "maf_lower_bound", "maf_upper_bound", "maf_multiplier",
    "maf_above_upper",
}
STUDY_WEIGHTS = MAF_WEIGHTS | {
    "conserved_200m", "conserved_100v", "qtl",
    "literature", "selection_pur", "selection_local",
    "selection_domestication",
}


def read_weights(path, required):
    """Read ``key<TAB>value`` weights, allowing comments and a header."""
    weights = {}
    with text_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 2 or fields[0].strip().lower() == "key":
                continue
            try:
                weights[fields[0].strip()] = float(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid weight at {path}:{line_no}") from exc
    missing = set(required) - set(weights)
    if missing:
        raise ValueError(f"weight table is missing keys: {sorted(missing)}")
    return weights


def text_open(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_snp(value):
    m = SNP_RE.match(value.strip())
    return None if not m else (m.group(1), int(m.group(2)))


def site_key(value):
    p = parse_snp(value)
    return "" if p is None else f"{p[0]}:{p[1]}"


def read_windows(path):
    rows, by_id = [], {}
    with text_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 3:
                f = line.rstrip("\n\r").split()
            if len(f) < 3:
                continue
            try:
                chrom, start, end = f[0], int(f[1]), int(f[2])
            except ValueError:
                raise ValueError(f"invalid window at {path}:{line_no}")
            wid = f[3].strip() if len(f) >= 4 and f[3].strip() else f"{chrom}-{start}"
            if wid in by_id:
                raise ValueError(f"duplicate window ID {wid!r} in {path}")
            row = {"window_id": wid, "chrom": chrom, "start0": start, "end0": end}
            rows.append(row)
            by_id[wid] = row
    if not rows:
        raise ValueError(f"no windows found in {path}")
    return rows, by_id


def read_regions(path, prefix):
    out = []
    with text_open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 3:
                f = line.rstrip("\n\r").split()
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                raise ValueError(f"invalid interval at {path}:{line_no}")
            if end <= start:
                continue
            out.append((f[0], start, end,
                        f[3] if len(f) >= 4 and f[3] else f"{prefix}_{len(out)+1:04d}"))
    return out


def interval_index(regions):
    by_chrom = defaultdict(list)
    for chrom, start, end, name in regions:
        by_chrom[chrom].append((start, end, name))
    starts, prefix_end = {}, {}
    for chrom, values in by_chrom.items():
        values.sort(key=lambda x: (x[0], x[1], x[2]))
        merged = []
        for start, end, name in values:
            if merged and start <= merged[-1][1]:
                old_start, old_end, old_name = merged[-1]
                merged[-1] = (old_start, max(old_end, end), old_name + "," + name)
            else:
                merged.append((start, end, name))
        by_chrom[chrom] = merged
        starts[chrom] = [x[0] for x in merged]
        prefix_end[chrom] = [x[1] for x in merged]
    return dict(by_chrom), starts, prefix_end


def read_weighted_regions(path):
    """Read NAME, PATH, WEIGHT tracks from a tab-separated manifest."""
    if not path:
        return []
    manifest = Path(path).resolve()
    tracks, names = [], set()
    with text_open(str(manifest)) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if fields[0].strip().lower() == "name":
                continue
            if len(fields) != 3:
                raise ValueError(
                    f"weighted-region manifest needs name, path, weight at {path}:{line_no}"
                )
            name, region_path, raw_weight = (value.strip() for value in fields)
            if not name or name in names:
                raise ValueError(f"empty or duplicate weighted-region name {name!r}")
            resolved = Path(region_path)
            if not resolved.is_absolute():
                resolved = manifest.parent / resolved
            if not resolved.is_file():
                raise FileNotFoundError(f"weighted-region BED not found: {resolved}")
            try:
                weight = float(raw_weight)
            except ValueError as exc:
                raise ValueError(f"invalid region weight at {path}:{line_no}") from exc
            regions = read_regions(str(resolved), name)
            index, starts, prefix_end = interval_index(regions)
            tracks.append({"name": name, "path": str(resolved), "weight": weight,
                           "regions": regions, "index": index, "starts": starts,
                           "prefix_end": prefix_end})
            names.add(name)
    if not tracks:
        raise ValueError(f"no weighted-region tracks found in {path}")
    return tracks


def hits(chrom, pos1, index, starts, prefix_end):
    values = index.get(chrom, ())
    if not values:
        return []
    i = bisect.bisect_right(starts[chrom], pos1) - 1
    result = []
    while i >= 0:
        start, end, name = values[i]
        if end < pos1:
            break
        if start < pos1 <= end:
            result.append((start, end, name))
        i -= 1
    return result


def maf_score(alt_freq, weights):
    maf = min(alt_freq, 1.0 - alt_freq)
    lower = weights["maf_lower_bound"]
    upper = weights["maf_upper_bound"]
    score = (weights["maf_multiplier"] * maf if lower <= maf <= upper
             else (weights["maf_above_upper"] if maf > upper else 0.0))
    return maf, score


def read_maf(path, weights, include_maf=True):
    records, sites, invalid = {}, set(), 0
    with text_open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 6:
                invalid += 1
                continue
            try:
                chrom, pos = parse_snp(f[1])
                alt = float(f[5])
            except (TypeError, ValueError):
                invalid += 1
                continue
            maf, score = maf_score(alt, weights)
            site = f"{chrom}:{pos}"
            if site in sites:
                raise ValueError(
                    f"multiple ALT records at {site}; normalize to one "
                    "biallelic record per coordinate before scoring"
                )
            if not include_maf:
                score = 0.0
            records[f[1]] = (chrom, pos, maf, score, f[2], f[3])
            sites.add(site)
    if not records:
        raise ValueError(f"no valid SNPs in {path}")
    return records, sites, invalid


def add_track(path, site_scores, sites, component, value_fn, counts):
    if not os.path.isfile(path):
        counts[f"{component}_file_missing"] += 1
        return
    index = COMPONENTS.index(component)
    with text_open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 2:
                continue
            try:
                site = f"{f[0]}:{int(f[1])}"
                if site not in sites:
                    continue
                value = float(value_fn(f))
            except (IndexError, ValueError, TypeError):
                counts[f"{component}_invalid"] += 1
                continue
            site_scores.setdefault(site, [0.0] * len(COMPONENTS))[index] += value
            counts[f"{component}_hits"] += 1


def add_expression(path, site_scores, sites, component, last_index, counts):
    if not os.path.isfile(path):
        counts[f"{component}_file_missing"] += 1
        return
    index = COMPONENTS.index(component)
    with text_open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            if len(f) < 3:
                continue
            try:
                site = f"{f[0]}:{int(f[2])}"
                if site not in sites:
                    continue
                value = float(f[last_index])
            except (IndexError, ValueError):
                counts[f"{component}_invalid"] += 1
                continue
            site_scores.setdefault(site, [0.0] * len(COMPONENTS))[index] += value
            counts[f"{component}_hits"] += 1


def existing_sites(path):
    result = set()
    if not path or not os.path.isfile(path):
        return result
    with text_open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n\r").split("\t")
            for value in ([f[3]] if len(f) >= 4 else []) + [f[0]]:
                k = site_key(value)
                if k:
                    result.add(k)
                    break
            else:
                if len(f) >= 3:
                    try:
                        result.add(f"{f[0]}:{int(f[1])+1}")
                    except ValueError:
                        pass
    return result


def score_ld_tag_sets(text, records, site_scores, selection_hits, top_n):
    scored_tag_sets, missing = [], 0
    for tag_set_index, tag_set in enumerate(text.split(","), 1):
        by_site = {}
        for raw in tag_set.split("|"):
            raw = raw.strip()
            if not raw or raw == "NONE":
                continue
            rec = records.get(raw)
            if rec is None:
                missing += 1
                continue
            chrom, pos, _maf, maf_component, ref, alt = rec
            site = f"{chrom}:{pos}"
            comp = site_scores.get(site)
            base = maf_component + (sum(comp) if comp else 0.0)
            sel = selection_hits.get(site)
            bonus = sel[1] if sel else 0
            val = {"token": raw, "site": site, "chrom": chrom, "pos1": pos,
                   "ref": ref, "alt": alt, "base": base, "bonus": bonus,
                   "total": base + bonus, "category": sel[0] if sel else "none"}
            old = by_site.get(site)
            if old is None or val["total"] > old["total"]:
                by_site[site] = val
        snps = sorted(by_site.values(), key=lambda x: (-x["total"], -x["base"], x["pos1"], x["token"]))
        if snps:
            scored_tag_sets.append((sum(x["total"] for x in snps), len(snps), tag_set_index, snps))
    scored_tag_sets.sort(key=lambda x: (-x[0], -x[1], x[2]))
    candidates = []
    for tag_set_rank, (tag_set_score, _count, tag_set_index, snps) in enumerate(scored_tag_sets, 1):
        for snp_rank, value in enumerate(snps, 1):
            if len(candidates) >= top_n:
                break
            candidates.append((tag_set_rank, tag_set_index, tag_set_score, snp_rank, value))
        if len(candidates) >= top_n:
            break
    metrics = {"ld_tag_set_count": len(text.split(",")) if text else 0,
               "valid_ld_tag_set_count": len(scored_tag_sets),
               "scored_snp_count": sum(len(x[3]) for x in scored_tag_sets),
               "missing_score_token_count": missing}
    return (candidates[0] if candidates else None), candidates, metrics


def write_rows(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def aggregate_weighted(matches, method):
    """Return one category and bonus from ordered (name, weight, hits) items."""
    if not matches:
        return "none", 0.0
    if method == "first":
        selected = matches[:1]
    elif method == "max":
        selected = [max(matches, key=lambda item: item[1])]
    else:
        selected = matches
    category = "+".join(f"custom:{item[0]}" for item in selected)
    return category, sum(item[1] for item in selected)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("windows", "ld_tag_sets", "maf", "weights", "outdir"):
        flag = "--" + name.replace("_", "-")
        p.add_argument(flag, required=True)
    for name in ("annotation", "con200m", "con100v", "qtl", "literature",
                 "expression", "local_core", "pur_core", "domestication_core"):
        flag = "--" + name.replace("_", "-")
        p.add_argument(flag, default="")
    p.add_argument("--profile", choices=("study", "custom"), default="study")
    p.add_argument("--weighted-regions", default="",
                   help="TSV manifest with name, path, and weight columns")
    p.add_argument("--region-aggregation", choices=("sum", "max", "first"),
                   default="sum")
    p.add_argument("--disable-maf-score", action="store_true")
    p.add_argument("--top-per-window", type=int, default=50)
    p.add_argument("--current-sites", default="")
    p.add_argument("--exclude-chrom", action="append", default=[])
    a = p.parse_args()

    for name in ("windows", "ld_tag_sets", "maf", "weights"):
        if not os.path.isfile(getattr(a, name)):
            p.error(f"input file not found: {getattr(a, name)}")
    if a.profile == "study":
        required = ("annotation", "con200m", "con100v", "qtl", "literature",
                    "expression", "local_core", "pur_core", "domestication_core")
        for name in required:
            if not os.path.isfile(getattr(a, name)):
                p.error(f"study profile requires --{name.replace('_', '-')}")
    elif a.annotation and not os.path.isfile(a.annotation):
        p.error(f"input file not found: {a.annotation}")
    if a.weighted_regions and not os.path.isfile(a.weighted_regions):
        p.error(f"input file not found: {a.weighted_regions}")
    if a.top_per_window < 1:
        p.error("--top-per-window must be positive")
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    required_weights = STUDY_WEIGHTS if a.profile == "study" else MAF_WEIGHTS
    weights = read_weights(a.weights, required_weights)
    windows, by_id = read_windows(a.windows)
    custom_tracks = read_weighted_regions(a.weighted_regions)
    if a.profile == "study":
        local = read_regions(a.local_core, "local")
        pur = read_regions(a.pur_core, "pur")
        dom = read_regions(a.domestication_core, "domestication")
    else:
        local, pur, dom = [], [], []
    li, ls, lp = interval_index(local)
    pi, ps, pp = interval_index(pur)
    di, ds, dp = interval_index(dom)

    affected, affected_rows = set(), []
    for row in windows:
        def ov(regs):
            return [r for r in regs if r[0] == row["chrom"] and r[1] < row["end0"] and row["start0"] < r[2]]
        lh, ph, dh = ov(local), ov(pur), ov(dom)
        if not lh and not ph and not dh:
            continue
        affected.add(row["window_id"])
        category, bonus = (
            ("local_adaptation_PUR", weights["selection_pur"]) if ph else
            (("local_adaptation", weights["selection_local"]) if lh else
             ("domestication", weights["selection_domestication"]))
        )
        region_names = ([f"local:{x[3]}" for x in lh] +
                        [f"PUR:{x[3]}" for x in ph] +
                        [f"domestication:{x[3]}" for x in dh])
        affected_rows.append([row["window_id"], row["chrom"], row["start0"],
                              row["end0"], category, bonus, ",".join(region_names)])
    for row in windows:
        matches, region_names = [], []
        for track in custom_tracks:
            overlap = [r for r in track["regions"]
                       if r[0] == row["chrom"] and
                       r[1] < row["end0"] and row["start0"] < r[2]]
            if overlap:
                matches.append((track["name"], track["weight"], overlap))
                region_names.extend(f"{track['name']}:{value[3]}" for value in overlap)
        if matches:
            affected.add(row["window_id"])
            category, bonus = aggregate_weighted(matches, a.region_aggregation)
            affected_rows.append([row["window_id"], row["chrom"], row["start0"],
                                  row["end0"], category, bonus,
                                  ",".join(region_names)])
    affected_rows.sort(key=lambda x: (x[1], int(x[2]), x[4]))
    write_rows(out / "selection_affected_windows.tsv",
               ["window_id", "chrom", "start0", "end0", "region_categories",
                "region_bonus", "overlap_regions"], affected_rows)

    records, sites, maf_invalid = read_maf(
        a.maf, weights, include_maf=not a.disable_maf_score
    )
    site_scores, track_counts = {}, Counter()
    if a.annotation:
        add_track(a.annotation, site_scores, sites, "annotation_score",
                  lambda f: f[2], track_counts)
    if a.profile == "study":
        add_track(a.con200m, site_scores, sites, "con200m_score",
                  lambda f: weights["conserved_200m"], track_counts)
        add_track(a.con100v, site_scores, sites, "con100v_score",
                  lambda f: weights["conserved_100v"], track_counts)
        add_track(a.qtl, site_scores, sites, "qtl_score",
                  lambda f: weights["qtl"], track_counts)
        add_track(a.literature, site_scores, sites, "literature_score",
                  lambda f: weights["literature"], track_counts)
        add_expression(a.expression, site_scores, sites, "expression_breadth",
                       -2, track_counts)
        add_expression(a.expression, site_scores, sites, "expression_intensity",
                       -1, track_counts)

    selection_hits = {}
    for chrom, pos in ((v[0], v[1]) for v in records.values()):
        site = f"{chrom}:{pos}"
        if site in selection_hits:
            continue
        if a.profile == "study":
            ph = hits(chrom, pos, pi, ps, pp)
            lh = hits(chrom, pos, li, ls, lp)
            dh = hits(chrom, pos, di, ds, dp)
            if ph:
                selection_hits[site] = (
                    "local_adaptation_PUR", weights["selection_pur"]
                )
            elif lh:
                selection_hits[site] = (
                    "local_adaptation", weights["selection_local"]
                )
            elif dh:
                selection_hits[site] = (
                    "domestication", weights["selection_domestication"]
                )
        matches = []
        for track in custom_tracks:
            overlap = hits(chrom, pos, track["index"], track["starts"],
                           track["prefix_end"])
            if overlap:
                matches.append((track["name"], track["weight"], overlap))
                track_counts[f"custom_region_{track['name']}_hit_sites"] += 1
        if matches:
            category, bonus = aggregate_weighted(matches, a.region_aggregation)
            previous = selection_hits.get(site)
            if previous:
                selection_hits[site] = (
                    previous[0] + "+" + category, previous[1] + bonus
                )
            else:
                selection_hits[site] = (category, bonus)

    # Full score table is streamed to keep output memory bounded.
    score_header = ["snp_id", "chrom", "pos1", "ref", "alt", "maf", "base_score", "selection_bonus",
                    "total_score", "selection_category", *COMPONENTS, "has_heavy_track"]
    with open(out / "snp_filter_snp_scores.tsv", "w", encoding="utf-8", newline="") as fh, \
         open(out / "snp_filter_snp_heavy.txt", "w", encoding="utf-8") as compact:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(score_header)
        for sid in sorted(records):
            chrom, pos, maf, maf_component, ref, alt = records[sid]
            site = f"{chrom}:{pos}"
            comps = [0.0] * len(COMPONENTS)
            comps[0] = maf_component
            if site in site_scores:
                for i, value in enumerate(site_scores[site]):
                    comps[i] += value
            base = sum(comps)
            sel = selection_hits.get(site)
            bonus = sel[1] if sel else 0
            total = base + bonus
            writer.writerow([sid, chrom, pos, ref, alt, f"{maf:.8f}", f"{base:.8f}", bonus, f"{total:.8f}",
                             sel[0] if sel else "none", *[f"{x:.8f}" for x in comps],
                             "yes" if site in site_scores else "no"])
            compact.write(f"{sid}\t{total:.8f}\n")

    old_sites = existing_sites(a.current_sites)
    important_rows, location_all, location_probe, status_rows = [], [], [], []
    candidate_rows = []
    candidate_sites = set()
    seen, tag_set_rows_total = set(), 0
    malformed = unknown = duplicate = 0
    with text_open(a.ld_tag_sets) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            tag_set_rows_total += 1
            f = line.rstrip("\n\r").split("\t", 1)
            if len(f) < 2:
                malformed += 1
                continue
            wid, text = f
            if wid not in by_id:
                unknown += 1
                continue
            if wid in seen:
                duplicate += 1
                continue
            seen.add(wid)
            best, candidates, metrics = score_ld_tag_sets(text, records, site_scores, selection_hits, a.top_per_window)
            if best:
                tag_set_rank, tag_set_index, tag_set_score, snp_rank, val = best
                important_rows.append([wid, val["token"], f"{val['total']:.8f}", f"{tag_set_score:.8f}", tag_set_index, val["category"]])
                location_all.append([val["chrom"], val["pos1"] - 1, val["pos1"]])
                if val["chrom"] not in set(a.exclude_chrom):
                    location_probe.append([val["chrom"], val["pos1"] - 1, val["pos1"], f"{wid}|{val['token']}", f"{val['total']:.8f}"])
                status, top_site, top_cat = "best_candidate", val["site"], val["category"]
                for rank, (tag_set_rank, tag_set_index, tag_set_score, s_rank, v) in enumerate(candidates, 1):
                    candidate_rows.append([wid, v["chrom"], by_id[wid]["start0"], by_id[wid]["end0"], rank,
                                           tag_set_rank, tag_set_index, f"{tag_set_score:.8f}", s_rank, v["token"], v["ref"], v["alt"],
                                           v["pos1"], f"{v['base']:.8f}", v["bonus"], f"{v['total']:.8f}",
                                           v["category"], "yes" if v["site"] in old_sites else "no",
                                           f"cand{len(candidate_rows)+1:08d}"])
                    candidate_sites.add(v["site"])
            else:
                status, top_site, top_cat = "no_valid_snp_score", "", ""
            status_rows.append([wid, status, metrics["ld_tag_set_count"], metrics["valid_ld_tag_set_count"], top_site,
                                top_cat, metrics["scored_snp_count"], metrics["missing_score_token_count"]])
    for row in windows:
        if row["window_id"] not in seen:
            status_rows.append([row["window_id"], "missing_ld_tag_set_row", 0, 0, "", "", 0, 0])

    write_rows(out / "snp_filter_important_ld_tag_set_snp.tsv",
               ["window_id", "best_snp_id", "best_snp_total_score", "best_ld_tag_set_score", "ld_tag_set_index", "selection_category"], important_rows)
    with open(out / "snp_filter_important_ld_tag_set_snp.txt", "w", encoding="utf-8") as fh:
        for row in important_rows:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(out / "bigest_snp_location.all.txt", "w", encoding="utf-8") as fh:
        for row in location_all:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(out / "bigest_snp_location.txt", "w", encoding="utf-8") as fh, \
         open(out / "bigest_snp_location.bed", "w", encoding="utf-8") as bed:
        for row in location_probe:
            fh.write("\t".join(map(str, row[:3])) + "\n")
            bed.write("\t".join(map(str, row)) + "\n")

    candidate_header = ["window_id", "chrom", "window_start0", "window_end0", "candidate_rank", "ld_tag_set_rank",
                        "ld_tag_set_index", "ld_tag_set_score", "snp_rank", "snp_id", "ref", "alt", "pos1", "base_score",
                        "selection_bonus", "total_score", "selection_category", "in_existing_panel", "candidate_id"]
    write_rows(out / "selection_ld_tag_set_candidates.tsv", candidate_header, candidate_rows)
    write_rows(out / "snp_filter_window_status.tsv",
               ["window_id", "status", "ld_tag_set_count", "valid_ld_tag_set_count", "top_snp_id", "top_category",
                "scored_snp_count", "missing_score_token_count"], status_rows)
    with open(out / "selection_candidates.bed", "w", encoding="utf-8") as fh:
        for row in candidate_rows:
            fh.write(f"{row[1]}\t{int(row[12])-1}\t{row[12]}\t{row[18]}\n")

    summary = {
        "windows_total": len(windows), "windows_affected_by_selection": len(affected),
        "windows_with_ld_tag_set_rows": len(seen), "windows_with_best_candidate": len(important_rows),
        "windows_without_valid_candidate": len(windows) - len(important_rows),
        "new_population_snp_ids": len(records), "new_population_sites": len(sites),
        "selection_hit_sites": len(selection_hits), "selection_hit_site_categories": dict(Counter(v[0] for v in selection_hits.values())),
        "candidate_rows_all_windows": len(candidate_rows), "candidate_sites_unique": len(candidate_sites),
        "ld_tag_set_rows_total": tag_set_rows_total, "malformed_ld_tag_set_rows": malformed, "unknown_ld_tag_set_windows": unknown,
        "duplicate_ld_tag_set_rows": duplicate, "invalid_maf_lines": maf_invalid, "heavy_track_hits": dict(track_counts),
        "existing_panel_sites_loaded": len(old_sites), "excluded_chromosomes_from_probe_bed": sorted(set(a.exclude_chrom)),
        "scoring_profile": a.profile,
        "maf_scoring_enabled": not a.disable_maf_score,
        "vep_scoring_enabled": bool(a.annotation),
        "study_selection_policy": (
            {"pur_local_adaptation": weights["selection_pur"],
             "local_adaptation_core": weights["selection_local"],
             "domestication_core": weights["selection_domestication"],
             "priority": "PUR local-adaptation > local-adaptation > domestication"}
            if a.profile == "study" else None
        ),
        "custom_region_policy": {
            "aggregation": a.region_aggregation,
            "manifest": os.path.abspath(a.weighted_regions)
            if a.weighted_regions else None,
            "tracks": [
                {"name": track["name"], "path": track["path"],
                 "weight": track["weight"], "regions": len(track["regions"])}
                for track in custom_tracks
            ],
        },
        "weight_policy": {
            "source": os.path.abspath(a.weights),
            "annotation": "VEP score track column 3" if a.annotation else "disabled",
            "expression": ("last two columns summed"
                           if a.profile == "study" else "disabled"),
        },
    }
    with open(out / "score_filter_chip_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
