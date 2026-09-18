#!/usr/bin/env python3
"""Select one target-covering probe per site, with explicit rescue tiers."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


TIERS = (
    "strict_oligominer",
    "unique_only",
    "blockparse_only",
    "reference_centered_unfiltered",
)


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    site: str
    chrom: str
    pos1: int
    flank_start0: int
    target_local1: int
    flank_sequence: str = ""


@dataclass(frozen=True)
class Probe:
    candidate_id: str
    rel_start1: int
    rel_end1: int
    sequence: str
    tm: float | None
    tier: str

    def rank(self, target_local1: int) -> tuple[float, float, int, int, str]:
        midpoint = (self.rel_start1 + self.rel_end1) / 2.0
        tm_rank = -self.tm if self.tm is not None else 0.0
        return (abs(midpoint - target_local1), tm_rank,
                self.rel_start1, self.rel_end1, self.sequence)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True,
                        help="oligominer_candidate_flanks.with_status.tsv")
    parser.add_argument("--filtered-bed", required=True,
                        help="kmerFilter BED (strict tier)")
    parser.add_argument("--unique-bed", default="",
                        help="outputClean unique BED (optional rescue tier)")
    parser.add_argument("--blockparse-fastq", default="",
                        help="raw blockParse FASTQ (optional rescue tier)")
    parser.add_argument("--flanks-fasta", default="",
                        help="241-bp flank FASTA (optional final fallback)")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def canonical_site(value: str) -> str:
    fields = value.strip().split(":")
    if len(fields) < 2:
        return ""
    return f"{fields[0]}:{int(fields[1])}"


def read_fasta(path: str) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name:
                parts.append(line)
    if name:
        sequences[name] = "".join(parts).upper()
    return sequences


def read_candidates(path: str, flanks_fasta: str) -> tuple[list[dict[str, str]], dict[str, Candidate]]:
    flanks = read_fasta(flanks_fasta) if flanks_fasta else {}
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        required = {"candidate_id", "snp_id", "chrom", "pos1", "flank_start0"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"candidate table missing columns: {sorted(missing)}")
    metadata: dict[str, Candidate] = {}
    for row in rows:
        candidate_id = row["candidate_id"]
        pos1 = int(row["pos1"])
        flank_start0 = int(row["flank_start0"])
        target_local1 = (int(row["target_offset0"]) + 1
                         if row.get("target_offset0", "")
                         else pos1 - flank_start0)
        metadata[candidate_id] = Candidate(
            candidate_id=candidate_id,
            site=canonical_site(row["snp_id"]),
            chrom=row["chrom"],
            pos1=pos1,
            flank_start0=flank_start0,
            target_local1=target_local1,
            flank_sequence=flanks.get(candidate_id, ""),
        )
    return rows, metadata


def retain(best: dict[str, Probe], candidate: Candidate, probe: Probe) -> None:
    if not (probe.rel_start1 <= candidate.target_local1 <= probe.rel_end1):
        return
    if len(probe.sequence) != probe.rel_end1 - probe.rel_start1 + 1:
        return
    previous = best.get(candidate.candidate_id)
    if previous is None or probe.rank(candidate.target_local1) < previous.rank(candidate.target_local1):
        best[candidate.candidate_id] = probe


def read_bed(path: str, tier: str, metadata: dict[str, Candidate],
             excluded: set[str]) -> dict[str, Probe]:
    best: dict[str, Probe] = {}
    if not path:
        return best
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 4:
                continue
            candidate_id = fields[0]
            if candidate_id not in metadata or candidate_id in excluded:
                continue
            try:
                start1, end1 = int(fields[1]), int(fields[2])
                tm = float(fields[4]) if len(fields) >= 5 and fields[4] else None
            except ValueError:
                continue
            probe = Probe(candidate_id, start1, end1, fields[3].upper(), tm, tier)
            retain(best, metadata[candidate_id], probe)
    return best


def fastq_records(path: str) -> Iterable[tuple[str, str]]:
    with open(path, encoding="utf-8", errors="replace") as handle:
        while True:
            header = handle.readline().rstrip("\n\r")
            if not header:
                break
            sequence = handle.readline().rstrip("\n\r").upper()
            plus = handle.readline()
            quality = handle.readline()
            if header.startswith("@") and plus and quality:
                yield header[1:], sequence


def read_fastq(path: str, metadata: dict[str, Candidate],
               excluded: set[str]) -> dict[str, Probe]:
    best: dict[str, Probe] = {}
    if not path:
        return best
    for header, sequence in fastq_records(path):
        try:
            candidate_id, coordinates = header.split(":", 1)
            start1, end1 = (int(value) for value in coordinates.split("-", 1))
        except ValueError:
            continue
        if candidate_id not in metadata or candidate_id in excluded:
            continue
        probe = Probe(candidate_id, start1, end1, sequence, None, "blockparse_only")
        retain(best, metadata[candidate_id], probe)
    return best


def centered_probe(candidate: Candidate) -> Probe | None:
    sequence = candidate.flank_sequence
    if len(sequence) < 120:
        return None
    start0 = max(0, min(candidate.target_local1 - 1 - 59, len(sequence) - 120))
    probe_sequence = sequence[start0:start0 + 120]
    probe = Probe(candidate.candidate_id, start0 + 1, start0 + 120,
                  probe_sequence, None, "reference_centered_unfiltered")
    if probe.rel_start1 <= candidate.target_local1 <= probe.rel_end1:
        return probe
    return None


def select_by_tier(cfg: argparse.Namespace, metadata: dict[str, Candidate]) -> dict[str, Probe]:
    chosen = read_bed(cfg.filtered_bed, "strict_oligominer", metadata, set())
    unique = read_bed(cfg.unique_bed, "unique_only", metadata, set(chosen))
    chosen.update(unique)
    raw = read_fastq(cfg.blockparse_fastq, metadata, set(chosen))
    chosen.update(raw)
    if cfg.flanks_fasta:
        for candidate_id, candidate in metadata.items():
            if candidate_id in chosen:
                continue
            probe = centered_probe(candidate)
            if probe is not None:
                chosen[candidate_id] = probe
    return chosen


def main() -> None:
    cfg = parse_args()
    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows, metadata = read_candidates(cfg.candidates, cfg.flanks_fasta)
    chosen = select_by_tier(cfg, metadata)

    fields = [
        "snp_id", "candidate_id", "probe_status", "probe_tier", "probe_chrom",
        "probe_start", "probe_end", "probe_sequence", "probe_tm",
        "center_distance_bp", "strict_oligominer_pass",
    ]
    outputs: list[dict[str, str]] = []
    tier_counts = {tier: 0 for tier in TIERS}
    unresolved = 0
    for row in rows:
        candidate_id = row["candidate_id"]
        candidate = metadata[candidate_id]
        probe = chosen.get(candidate_id)
        if probe is None:
            unresolved += 1
            outputs.append({
                "snp_id": candidate.site, "candidate_id": candidate_id,
                "probe_status": "unresolved", "probe_tier": "unresolved",
                "probe_chrom": candidate.chrom, "probe_start": "", "probe_end": "",
                "probe_sequence": "", "probe_tm": "", "center_distance_bp": "",
                "strict_oligominer_pass": "no",
            })
            continue
        tier_counts[probe.tier] += 1
        start = candidate.flank_start0 + probe.rel_start1
        end = candidate.flank_start0 + probe.rel_end1
        center_distance = abs((start + end) / 2.0 - candidate.pos1)
        outputs.append({
            "snp_id": candidate.site, "candidate_id": candidate_id,
            "probe_status": "designed", "probe_tier": probe.tier,
            "probe_chrom": candidate.chrom, "probe_start": str(start),
            "probe_end": str(end), "probe_sequence": probe.sequence,
            "probe_tm": "" if probe.tm is None else f"{probe.tm:.2f}",
            "center_distance_bp": f"{center_distance:.1f}",
            "strict_oligominer_pass": "yes" if probe.tier == "strict_oligominer" else "no",
        })

    all_path = outdir / "missing_probe_designs.all.tsv"
    with all_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(outputs)
    strict_path = outdir / "missing_probe_designs.strict.tsv"
    with strict_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row in outputs if row["probe_tier"] == "strict_oligominer")
    rescue_path = outdir / "missing_probe_designs.rescue.tsv"
    with rescue_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row in outputs
                         if row["probe_tier"] not in ("strict_oligominer", "unresolved"))

    summary_path = outdir / "missing_probe_designs.summary.tsv"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"candidate_sites\t{len(rows)}\n")
        for tier in TIERS:
            handle.write(f"selected_{tier}\t{tier_counts[tier]}\n")
        handle.write(f"strict_pass_sites\t{tier_counts['strict_oligominer']}\n")
        handle.write(f"rescue_only_sites\t{sum(tier_counts[tier] for tier in TIERS[1:])}\n")
        handle.write(f"unresolved_sites\t{unresolved}\n")

    print(f"candidate_sites\t{len(rows)}")
    for tier in TIERS:
        print(f"selected_{tier}\t{tier_counts[tier]}")
    print(f"unresolved_sites\t{unresolved}")


if __name__ == "__main__":
    main()
