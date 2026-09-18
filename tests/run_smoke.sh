#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TEST_TMP=$(mktemp -d "${TMPDIR:-/tmp}/monsnp-pipeline-test.XXXXXX")
trap 'rm -rf -- "$TEST_TMP"' EXIT

python3 "$ROOT/workflow/scripts/vep_to_scores.py" \
    --variants "$ROOT/tests/data/variants.afreq" \
    --vep "$ROOT/tests/data/vep.tsv" \
    --weights "$ROOT/resources/scoring_weights.tsv" \
    --output "$TEST_TMP/01_annotation/vep_scores.tsv.gz" \
    --summary "$TEST_TMP/01_annotation/vep_scores.summary.json"

python3 "$ROOT/workflow/scripts/score_candidates.py" \
    --windows "$ROOT/tests/data/windows.bed" \
    --ld-tag-sets "$ROOT/tests/data/ld_tag_sets.tsv" \
    --maf "$ROOT/tests/data/variants.afreq" \
    --annotation "$TEST_TMP/01_annotation/vep_scores.tsv.gz" \
    --con200m "$ROOT/tests/data/conserved_200m.tsv" \
    --con100v "$ROOT/tests/data/conserved_100v.tsv" \
    --qtl "$ROOT/tests/data/qtl.tsv" \
    --literature "$ROOT/tests/data/literature.tsv" \
    --expression "$ROOT/tests/data/expression.tsv" \
    --local-core "$ROOT/tests/data/local_adaptation.bed" \
    --pur-core "$ROOT/tests/data/pur_local_adaptation.bed" \
    --domestication-core "$ROOT/tests/data/domestication.bed" \
    --weights "$ROOT/resources/scoring_weights.tsv" \
    --outdir "$TEST_TMP/02_scoring" --top-per-window 3

python3 "$ROOT/workflow/scripts/prepare_oligominer_candidates.py" \
    --candidates "$TEST_TMP/02_scoring/selection_ld_tag_set_candidates.tsv" \
    --reference "$ROOT/tests/data/reference.fa" \
    --fai "$ROOT/tests/data/reference.fa.fai" \
    --outdir "$TEST_TMP/03_probe_candidates" --max-rank 3

python3 "$ROOT/workflow/scripts/assemble_panel.py" \
    --candidates "$TEST_TMP/02_scoring/selection_ld_tag_set_candidates.tsv" \
    --strict-designs "$ROOT/tests/data/strict_designs.tsv" \
    --output "$TEST_TMP/06_final_panel/final_strict_probe_panel.tsv" \
    --bed "$TEST_TMP/06_final_panel/final_strict_probe_panel.bed" \
    --unresolved "$TEST_TMP/06_final_panel/unresolved_windows.tsv" \
    --summary "$TEST_TMP/06_final_panel/summary.json" \
    --summary-report "$TEST_TMP/06_final_panel/summary_report.txt"

python3 "$ROOT/tests/assert_smoke.py" "$TEST_TMP"

CUSTOM_ROOT="$TEST_TMP/plink_custom"
"$ROOT/tests/mock_bin/plink2" \
    --bfile "$ROOT/tests/data/plink/tiny" \
    --freq --out "$CUSTOM_ROOT/00_plink/variants"

python3 "$ROOT/workflow/scripts/build_plink_tag_sets.py" \
    --bfile "$ROOT/tests/data/plink/tiny" \
    --windows "$ROOT/tests/data/windows.bed" \
    --plink "$ROOT/tests/mock_bin/plink" \
    --output "$CUSTOM_ROOT/00_plink/window_ld_tag_sets.tsv" \
    --summary "$CUSTOM_ROOT/00_plink/window_ld_tag_sets.summary.json" \
    --jobs 1 --memory-mb 128 --window-boundary bed

python3 "$ROOT/workflow/scripts/score_candidates.py" \
    --windows "$ROOT/tests/data/windows.bed" \
    --ld-tag-sets "$CUSTOM_ROOT/00_plink/window_ld_tag_sets.tsv" \
    --maf "$CUSTOM_ROOT/00_plink/variants.afreq" \
    --weights "$ROOT/resources/scoring_weights.tsv" \
    --outdir "$CUSTOM_ROOT/02_scoring" \
    --top-per-window 3 --profile custom \
    --weighted-regions "$ROOT/tests/data/custom_regions.tsv" \
    --region-aggregation sum

python3 "$ROOT/tests/assert_plink_custom.py" "$CUSTOM_ROOT"

TOOLS_ROOT="$TEST_TMP/input_tools"
mkdir -p "$TOOLS_ROOT"
python3 "$ROOT/tools/prepare_plink_input.py" \
    --bfile "$ROOT/tests/data/plink/tiny" \
    --reference "$ROOT/tests/data/reference.fa" \
    --out "$TOOLS_ROOT/prepared/tiny" \
    --plink2 "$ROOT/tests/mock_bin/plink2" --chr-set 31
python3 "$ROOT/tools/make_windows.py" \
    --reference "$ROOT/tests/data/reference.fa" \
    --fai "$TOOLS_ROOT/generated_reference.fa.fai" \
    --output "$TOOLS_ROOT/windows.bed" --window-size 400 --step 400
python3 "$ROOT/tools/preflight.py" \
    --config "$ROOT/tests/config_plink_custom.yaml" --base-dir "$ROOT" \
    --report "$TOOLS_ROOT/preflight.json"
python3 "$ROOT/tests/assert_input_tools.py" "$TOOLS_ROOT"
