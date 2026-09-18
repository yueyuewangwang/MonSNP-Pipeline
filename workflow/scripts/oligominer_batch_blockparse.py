#!/usr/bin/env python
"""Run the Python-2 OligoMiner block parser on every FASTA record.

The upstream blockParse.py intentionally supports one FASTA record per call;
this driver keeps one Python process but invokes the same library routine once
per candidate flank, then concatenates the resulting FASTQ records.
"""

from __future__ import print_function

import argparse
import contextlib
import imp
import os
import shutil
import sys
import sys
import tempfile

from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt


@contextlib.contextmanager
def redirect_stdout(stream):
    previous = sys.stdout
    sys.stdout = stream
    try:
        yield
    finally:
        sys.stdout = previous


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--oligominer", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    module_path = os.path.join(args.oligominer, "blockParse.py")
    block_parse = imp.load_source("oligominer_blockParse", module_path)
    temp_dir = tempfile.mkdtemp(prefix="oligominer_blockparse_")
    count = 0
    records = 0
    try:
        with open(args.out, "w") as combined:
            for record in SeqIO.parse(args.fasta, "fasta"):
                records += 1
                fasta_path = os.path.join(temp_dir, "record.fa")
                out_stem = os.path.join(temp_dir, "record")
                with open(fasta_path, "w") as handle:
                    handle.write(">%s\n%s\n" % (record.id, str(record.seq).upper()))
                # blockParse prints a short per-record summary unconditionally;
                # suppress it so the batch log remains useful.
                with open(os.devnull, "w") as sink:
                    with redirect_stdout(sink):
                        block_parse.runSequenceCrawler(
                            fasta_path, 120, 120, 20, 80, mt.DNA_NN3,
                            40, 200, "AAAAAA,TTTTTT,CCCCCC,GGGGGG",
                            390, 50, 0, 25, 25, None, False, True,
                            False, False, False, False, out_stem)
                fastq_path = out_stem + ".fastq"
                if os.path.exists(fastq_path):
                    with open(fastq_path) as source:
                        data = source.read()
                    if data:
                        combined.write(data)
                        if not data.endswith("\n"):
                            combined.write("\n")
                        # blockParse joins records without a trailing newline;
                        # count the final line as well before dividing FASTQ
                        # records into four lines.
                        count += (data.count("\n") + 1) // 4
                if records % 500 == 0:
                    print("processed %d FASTA records; %d candidate probes" % (records, count))
        with open(args.out + ".complete", "w") as marker:
            marker.write("records\t%d\nprobes\t%d\n" % (records, count))
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
    print("batch blockParse complete: %d records, %d candidate probes" % (records, count))


if __name__ == "__main__":
    main()
