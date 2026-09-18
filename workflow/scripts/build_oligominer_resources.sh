#!/usr/bin/env bash
set -euo pipefail

REFERENCE=""
INDEX_PREFIX=""
JELLYFISH_DB=""
THREADS=4
MARKER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference) REFERENCE=$2; shift 2 ;;
        --index) INDEX_PREFIX=$2; shift 2 ;;
        --jellyfish) JELLYFISH_DB=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --marker) MARKER=$2; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -n "$REFERENCE" && -n "$INDEX_PREFIX" && -n "$JELLYFISH_DB" && -n "$MARKER" ]] || {
    echo "Required options: --reference --index --jellyfish --marker" >&2
    exit 2
}

for tool in bowtie2-build jellyfish; do
    command -v "$tool" >/dev/null || { echo "Missing executable: $tool" >&2; exit 1; }
done

mkdir -p "$(dirname "$INDEX_PREFIX")" "$(dirname "$MARKER")"
if [[ ! -f "${INDEX_PREFIX}.1.bt2" && ! -f "${INDEX_PREFIX}.1.bt2l" ]]; then
    bowtie2-build --threads "$THREADS" "$REFERENCE" "$INDEX_PREFIX"
fi
if [[ ! -f "$JELLYFISH_DB" && ! -f "${JELLYFISH_DB}_0" ]]; then
    jellyfish count -s 3G -m 25 -o "$JELLYFISH_DB" --out-counter-len 1 \
        -L 2 -t "$THREADS" "$REFERENCE"
fi

{
    printf 'reference\t%s\n' "$REFERENCE"
    printf 'bowtie2_index\t%s\n' "$INDEX_PREFIX"
    printf 'jellyfish_25mer\t%s\n' "$JELLYFISH_DB"
} > "$MARKER"
