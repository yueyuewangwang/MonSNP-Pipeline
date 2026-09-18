# Configuration

Run from the repository root and pass a project YAML with `--configfile`.
Relative paths are resolved from the repository root.

## Genotype source

The root-level `tools/` directory contains independent pre-processing
utilities used before the main workflow. For a new population, first run
`tools/prepare_plink_input.py` and
`tools/make_windows.py`, then validate the project with `tools/preflight.py`.
See [`docs/custom_data_quickstart.md`](../docs/custom_data_quickstart.md).

`genotype_source: plink_bfile` starts from `inputs.bfile`, the prefix shared by
`PREFIX.bed`, `PREFIX.bim`, and `PREFIX.fam`. The workflow creates AFREQ and
per-window PLINK 1.9 LD tag sets. This is the recommended entry point.

`genotype_source: precomputed` skips those jobs and requires:

| Key | Format |
|---|---|
| `inputs.afreq` | Biallelic PLINK2 `.afreq` |
| `inputs.ld_tag_sets` | Headerless two-column window LD tag set TSV |


## Scoring profiles

`scoring.profile: study` reproduces the original scoring policy and requires
VEP, two conservation tracks, QTL, literature, expression, local-adaptation,
PUR-local-adaptation, and domestication inputs.

`scoring.profile: custom` requires none of those study tracks. Configure:

| Key | Meaning |
|---|---|
| `include_maf` | Apply the MAF function from `resources/scoring_weights.tsv` |
| `include_vep` | Add VEP scores; requires `inputs.vep` |
| `weighted_regions` | List of `{name, path, weight}` BED tracks |
| `region_aggregation` | `sum`, `max`, or `first` for multi-file overlap |

See [`custom_regions.example.yaml`](custom_regions.example.yaml) and
[`docs/custom_scoring.md`](../docs/custom_scoring.md). The study profile can
also add `weighted_regions`; doing so is a sensitivity analysis, not an exact
published-score reproduction.

## Common inputs

| Key | Format | Purpose |
|---|---|---|
| `reference` | FASTA | Assembly used for coordinates and probes |
| `reference_fai` | samtools FAI | Index for the exact FASTA |
| `windows` | BED3/BED4 | Design windows |
| `current_sites` | optional site table | Existing panel membership |

All coordinate-bearing files must use the same assembly and chromosome names.
BIM IDs must be `CHROM:POS:REF:ALT`, with BIM A1=ALT and A2=REF.

## Chromosome scope

Use `target_chromosomes: all` to resolve every chromosome present in
`inputs.windows`, excluding `parameters.exclude_chromosomes`. For a real-data
trial on one chromosome, use an explicit list:

```yaml
target_chromosomes:
  - ECA1
```

LD calculations and OligoMiner screening are scattered across this list and
gathered into the canonical output paths. Explicit chromosome names must occur
in the design-window file and must not also be excluded.

## LD parameters

The original analysis settings are exposed under `parameters`:

```yaml
ld_maf: 0.05
ld_tag_kb: 100
ld_tag_r2: 0.8
ld_jobs: 4
ld_memory_mb_per_job: 2048
ld_window_boundary: study_compatible
```

`study_compatible` preserves the original `--from-bp START0 --to-bp END0`
calls. Use `bed` for strict 0-based, half-open BED conversion (`START0+1`
through `END0`). Total potential PLINK memory is approximately
`ld_jobs * ld_memory_mb_per_job`.

## Workflow stages

The public stage targets are `run_annotation` (frequencies, target windows,
and optional VEP), `run_scoring` (LD tag sets and ranked candidates), and
`run_full` (OligoMiner plus the final panel). `run_probe_design: false` makes
the default target stop after scoring; `true` makes it complete the full
workflow. Explicit stage targets override that default choice.

`candidates_per_window` controls both retained ranking depth and probe-screen
depth. The published alternative screen used 10.
