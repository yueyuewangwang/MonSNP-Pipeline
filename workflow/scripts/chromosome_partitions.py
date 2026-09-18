#!/usr/bin/env python3
"""Create deterministic chromosome shards for windows, AFREQ, and LD tag sets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("windows", "afreq"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--input", required=True)
        sub.add_argument("--chromosome", required=True)
        sub.add_argument("--output", required=True)
    ld = subparsers.add_parser("ld-tag-sets")
    ld.add_argument("--input", required=True)
    ld.add_argument("--windows", required=True)
    ld.add_argument("--chromosome", required=True)
    ld.add_argument("--output", required=True)
    ld.add_argument("--summary", required=True)
    return parser.parse_args()


def output_handle(path: str):
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination.open("w", encoding="utf-8")


def slice_windows(cfg: argparse.Namespace) -> None:
    records = 0
    with open(cfg.input, encoding="utf-8", errors="replace") as source, output_handle(cfg.output) as output:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            if line.rstrip("\n\r").split("\t", 1)[0] == cfg.chromosome:
                output.write(line)
                records += 1
    if not records:
        raise ValueError(f"no design windows found for target chromosome {cfg.chromosome}")


def slice_afreq(cfg: argparse.Namespace) -> None:
    records = 0
    header = ""
    with open(cfg.input, encoding="utf-8", errors="replace") as source, output_handle(cfg.output) as output:
        for line in source:
            if not line.strip():
                continue
            if line.startswith("#"):
                if not header:
                    header = line
                    output.write(line)
                continue
            fields = line.rstrip("\n\r").split()
            if fields and fields[0] == cfg.chromosome:
                output.write(line)
                records += 1
    if not header:
        raise ValueError(f"AFREQ input has no #CHROM header: {cfg.input}")
    if not records:
        raise ValueError(f"no AFREQ variants found for target chromosome {cfg.chromosome}")


def slice_ld_tag_sets(cfg: argparse.Namespace) -> None:
    window_ids = set()
    with open(cfg.windows, encoding="utf-8", errors="replace") as windows:
        for line in windows:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            window_ids.add(fields[3] if len(fields) >= 4 and fields[3] else f"{fields[0]}-{fields[1]}")
    found = set()
    tag_sets = 0
    empty = 0
    with open(cfg.input, encoding="utf-8", errors="replace") as source, output_handle(cfg.output) as output:
        for line_no, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t", 1)
            if len(fields) != 2:
                raise ValueError(f"LD tag-set line {line_no} does not contain two columns")
            if fields[0] not in window_ids:
                continue
            if fields[0] in found:
                raise ValueError(f"duplicate LD tag-set window ID: {fields[0]}")
            found.add(fields[0])
            output.write(line)
            values = [value for value in fields[1].split(",") if value]
            tag_sets += len(values)
            empty += int(not values)
    missing = window_ids - found
    if missing:
        raise ValueError(
            f"precomputed LD tag sets omit {len(missing)} windows on {cfg.chromosome}; "
            f"examples: {sorted(missing)[:3]}"
        )
    summary = {
        "method": "precomputed LD tag sets",
        "chromosome": cfg.chromosome,
        "windows": len(window_ids),
        "windows_without_ld_tag_sets": empty,
        "ld_tag_sets_total": tag_sets,
    }
    destination = Path(cfg.summary)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    cfg = parse_args()
    if cfg.command == "windows":
        slice_windows(cfg)
    elif cfg.command == "afreq":
        slice_afreq(cfg)
    else:
        slice_ld_tag_sets(cfg)


if __name__ == "__main__":
    main()
