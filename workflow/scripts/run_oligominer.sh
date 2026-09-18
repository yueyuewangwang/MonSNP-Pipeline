#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: run_oligominer.sh --input-dir DIR --reference FASTA --outdir DIR \
  --oligominer DIR [--python PYTHON] [--threads N] [--index PREFIX] \
  [--jellyfish DB]
EOF
    exit 2
}

INPUT_DIR=""
REFERENCE=""
OUTDIR=""
OLIGOMINER=""
PYTHON_BIN="python"
THREADS=4
INDEX_PREFIX=""
JELLYFISH_DB=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir) INPUT_DIR=$2; shift 2 ;;
        --reference) REFERENCE=$2; shift 2 ;;
        --outdir) OUTDIR=$2; shift 2 ;;
        --oligominer) OLIGOMINER=$2; shift 2 ;;
        --python) PYTHON_BIN=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --index) INDEX_PREFIX=$2; shift 2 ;;
        --jellyfish) JELLYFISH_DB=$2; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -n "$INPUT_DIR" && -n "$REFERENCE" && -n "$OUTDIR" && -n "$OLIGOMINER" ]] || usage
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
FASTA="$INPUT_DIR/oligominer_candidate_flanks.fa"
CANDIDATES="$INPUT_DIR/oligominer_candidate_flanks.tsv"
BLOCKPARSE="$OLIGOMINER/blockParse.py"
OUTPUTCLEAN="$OLIGOMINER/outputClean.py"
KMERFILTER="$OLIGOMINER/kmerFilter.py"

for file in "$FASTA" "$CANDIDATES" "$BLOCKPARSE" "$OUTPUTCLEAN" "$KMERFILTER"; do
    [[ -f "$file" ]] || { echo "Missing file: $file" >&2; exit 1; }
done
for tool in "$PYTHON_BIN" bowtie2 bowtie2-build jellyfish; do
    command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

finish_empty_screen() {
    : > "$OUTDIR/blockparse.flt.kmer.bed"
    "$PYTHON_BIN" "$SCRIPT_DIR/merge_oligominer_results.py" \
        --candidates "$CANDIDATES" \
        --filtered-bed "$OUTDIR/blockparse.flt.kmer.bed" \
        --output "$OUTDIR/selection_probe_status.tsv"
}

if [[ ! -s "$FASTA" ]]; then
    echo "No candidate flanks in this chromosome shard; writing empty screening results"
    finish_empty_screen
    exit 0
fi

BLOCK_PREFIX="$OUTDIR/blockparse"
"$PYTHON_BIN" "$SCRIPT_DIR/oligominer_batch_blockparse.py" \
    --fasta "$FASTA" --oligominer "$OLIGOMINER" --out "$BLOCK_PREFIX.fastq"
if [[ ! -s "$BLOCK_PREFIX.fastq" ]]; then
    echo "OligoMiner generated no candidate probes; marking the shard unresolved"
    finish_empty_screen
    exit 0
fi

if [[ -z "$INDEX_PREFIX" ]]; then INDEX_PREFIX="$OUTDIR/reference_index"; fi
if [[ ! -f "${INDEX_PREFIX}.1.bt2" && ! -f "${INDEX_PREFIX}.1.bt2l" ]]; then
    bowtie2-build --threads "$THREADS" "$REFERENCE" "$INDEX_PREFIX"
fi
bowtie2 -x "$INDEX_PREFIX" -U "$BLOCK_PREFIX.fastq" --no-hd -k 2 \
    --very-sensitive-local -p "$THREADS" -S "$OUTDIR/blockparse.sam"

"$PYTHON_BIN" "$OUTPUTCLEAN" -u -f "$OUTDIR/blockparse.sam" \
    -o "$OUTDIR/blockparse.uniq"
if [[ ! -s "$OUTDIR/blockparse.uniq.bed" ]]; then
    echo "No uniquely aligned probes; marking the shard unresolved"
    finish_empty_screen
    exit 0
fi

if [[ -z "$JELLYFISH_DB" ]]; then JELLYFISH_DB="$OUTDIR/reference.25mer.jf"; fi
if [[ ! -f "$JELLYFISH_DB" && ! -f "${JELLYFISH_DB}_0" ]]; then
    jellyfish count -s 3G -m 25 -o "$JELLYFISH_DB" --out-counter-len 1 \
        -L 2 -t "$THREADS" "$REFERENCE"
fi
if [[ -f "${JELLYFISH_DB}_0" ]]; then
    JF_QUERY="${JELLYFISH_DB}_0"
else
    JF_QUERY="$JELLYFISH_DB"
fi
"$PYTHON_BIN" "$KMERFILTER" -f "$OUTDIR/blockparse.uniq.bed" -m 25 \
    -j "$JF_QUERY" -k 3 -o "$OUTDIR/blockparse.flt.kmer"
if [[ ! -s "$OUTDIR/blockparse.flt.kmer.bed" ]]; then
    echo "No probes passed k-mer filtering; marking the shard unresolved"
    finish_empty_screen
    exit 0
fi

"$PYTHON_BIN" "$SCRIPT_DIR/merge_oligominer_results.py" \
    --candidates "$CANDIDATES" \
    --filtered-bed "$OUTDIR/blockparse.flt.kmer.bed" \
    --output "$OUTDIR/selection_probe_status.tsv"
