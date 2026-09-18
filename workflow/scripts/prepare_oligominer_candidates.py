#!/usr/bin/env python3
"""Prepare allele-aware flank FASTA/BED input for an OligoMiner prescreen."""

from __future__ import annotations

import argparse
import csv
import os


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidates", required=True, help="selection_ld_tag_set_candidates.tsv")
    p.add_argument("--reference", required=True, help="reference FASTA")
    p.add_argument("--fai", default="", help="FASTA index; defaults to <reference>.fai")
    p.add_argument("--outdir", required=True)
    p.add_argument("--chromosome", default="", help="optional chromosome shard")
    p.add_argument("--max-rank", type=int, default=50)
    p.add_argument("--flank-left", type=int, default=121)
    p.add_argument("--flank-right", type=int, default=120)
    return p.parse_args()


class FastaRandomAccess:
    def __init__(self, fasta: str, fai: str):
        self.handle = open(fasta, "rb")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with open(fai, encoding="utf-8") as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                self.index[fields[0]] = tuple(map(int, fields[1:5]))

    def fetch(self, chrom: str, start0: int, end0: int) -> str:
        if chrom not in self.index:
            raise KeyError(f"reference chromosome not in FAI: {chrom}")
        length, offset, line_bases, line_width = self.index[chrom]
        start0 = max(0, start0)
        end0 = min(length, end0)
        if end0 <= start0:
            return ""
        out = bytearray()
        cursor = start0
        while cursor < end0:
            line_no, col = divmod(cursor, line_bases)
            take = min(end0 - cursor, line_bases - col)
            self.handle.seek(offset + line_no * line_width + col)
            out.extend(self.handle.read(take))
            cursor += take
        return bytes(out).decode("ascii").upper()

    def close(self) -> None:
        self.handle.close()


def main() -> None:
    cfg = args()
    if cfg.max_rank < 1 or cfg.flank_left < 1 or cfg.flank_right < 1:
        raise ValueError("rank and flank lengths must be positive")
    fai = cfg.fai or cfg.reference + ".fai"
    os.makedirs(cfg.outdir, exist_ok=True)
    reader = FastaRandomAccess(cfg.reference, fai)
    bed_path = os.path.join(cfg.outdir, "oligominer_candidate_flanks.bed")
    site_bed_path = os.path.join(cfg.outdir, "oligominer_candidate_sites.bed")
    table_path = os.path.join(cfg.outdir, "oligominer_candidate_flanks.tsv")
    fasta_path = os.path.join(cfg.outdir, "oligominer_candidate_flanks.fa")
    selected = 0
    written = 0
    mismatches = 0
    records = []
    try:
        with open(cfg.candidates, encoding="utf-8", newline="") as source:
            reader_csv = csv.DictReader(source, delimiter="\t")
            required = {"candidate_id", "window_id", "snp_id", "chrom", "pos1", "ref", "alt", "candidate_rank"}
            missing = required - set(reader_csv.fieldnames or [])
            if missing:
                raise ValueError(f"candidate table missing columns: {sorted(missing)}")
            for row in reader_csv:
                if cfg.chromosome and row["chrom"] != cfg.chromosome:
                    continue
                rank = int(row["candidate_rank"])
                if rank > cfg.max_rank:
                    continue
                selected += 1
                chrom = row["chrom"]
                pos1 = int(row["pos1"])
                if chrom not in reader.index:
                    row["probe_status"] = "missing_reference_chromosome"
                    records.append(row)
                    continue
                chrom_len = reader.index[chrom][0]
                start0 = max(0, pos1 - cfg.flank_left)
                end0 = min(chrom_len, pos1 + cfg.flank_right)
                sequence = reader.fetch(chrom, start0, end0)
                target_offset = pos1 - 1 - start0
                reference_base = sequence[target_offset] if 0 <= target_offset < len(sequence) else ""
                expected_ref = row["ref"].upper()
                if expected_ref and reference_base != expected_ref:
                    row["probe_status"] = "reference_allele_mismatch"
                    mismatches += 1
                else:
                    row["probe_status"] = "ready_for_oligominer"
                    written += 1
                row["flank_chrom"] = chrom
                row["flank_start0"] = str(start0)
                row["flank_end0"] = str(end0)
                row["target_offset0"] = str(target_offset)
                row["reference_base"] = reference_base
                row["flank_length"] = str(len(sequence))
                records.append((row, sequence))
    finally:
        reader.close()

    with open(bed_path, "w", encoding="utf-8") as bed, open(site_bed_path, "w", encoding="utf-8") as site_bed, open(fasta_path, "w", encoding="utf-8") as fasta:
        for item in records:
            row, sequence = item if isinstance(item, tuple) else (item, "")
            if not sequence or row.get("probe_status") == "reference_allele_mismatch":
                continue
            candidate_id = row["candidate_id"]
            bed.write(f"{row['flank_chrom']}\t{row['flank_start0']}\t{row['flank_end0']}\t{candidate_id}\n")
            site_bed.write(f"{row['flank_chrom']}\t{int(row['pos1']) - 1}\t{row['pos1']}\t{candidate_id}\n")
            fasta.write(f">{candidate_id}\n{sequence}\n")

    # Keep the original candidate columns first and append reproducibility fields.
    original_fields = []
    with open(cfg.candidates, encoding="utf-8", newline="") as source:
        original_fields = list(csv.DictReader(source, delimiter="\t").fieldnames or [])
    extras = [
        "probe_status", "flank_chrom", "flank_start0", "flank_end0",
        "target_offset0", "reference_base", "flank_length",
    ]
    fields = original_fields + [x for x in extras if x not in original_fields]
    with open(table_path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for item in records:
            row = item[0] if isinstance(item, tuple) else item
            writer.writerow(row)

    with open(os.path.join(cfg.outdir, "prepare_summary.tsv"), "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"selected_candidates\t{selected}\n")
        out.write(f"reference_checked_candidates\t{written + mismatches}\n")
        out.write(f"reference_allele_mismatches\t{mismatches}\n")
        out.write(f"fasta_records_written\t{written}\n")
        out.write(f"flank_left\t{cfg.flank_left}\n")
        out.write(f"flank_right\t{cfg.flank_right}\n")


if __name__ == "__main__":
    main()
