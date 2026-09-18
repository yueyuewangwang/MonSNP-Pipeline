#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

for tool in python3 plink plink2 snakemake; do
    command -v "$tool" >/dev/null || {
        echo "Missing executable: $tool" >&2
        echo "Create and activate environment.yaml before running the demo." >&2
        exit 1
    }
done

mkdir -p demo/work/raw demo/work/prepared demo/work/reports

plink --file demo/input/demo \
    --make-bed --allow-extra-chr --chr-set 31 \
    --out demo/work/raw/demo

python3 tools/prepare_plink_input.py \
    --bfile demo/work/raw/demo \
    --reference demo/input/reference.fa \
    --chrom-map demo/input/chromosome_map.tsv \
    --out demo/work/prepared/demo \
    --chr-set 31 --threads 1 --force

python3 tools/make_windows.py \
    --reference demo/input/reference.fa \
    --output demo/work/windows.bed \
    --window-size 400 --step 400 --force

python3 tools/preflight.py \
    --config demo/config.yaml --base-dir . \
    --report demo/work/reports/preflight.json

snakemake -s workflow/Snakefile \
    --configfile demo/config.yaml --cores 1 \
    --rerun-incomplete \
    --runtime-source-cache-path "$ROOT/demo/work/.snakemake-source-cache" \
    run_scoring

python3 demo/assert_demo.py
echo
echo "Demo completed. Main candidate table:"
echo "  demo/work/results/02_scoring/selection_ld_tag_set_candidates.tsv"
