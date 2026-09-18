# Custom Weighted-Region Scoring

Use the `custom` profile when the study-specific conservation, QTL,
literature, expression, adaptation, and domestication tracks are not
appropriate for a new project.

## Minimal configuration

```yaml
genotype_source: plink_bfile

inputs:
  bfile: data/snp_filter
  reference: data/assembly.fa
  reference_fai: data/assembly.fa.fai
  windows: data/design_windows.bed

scoring:
  profile: custom
  include_maf: true
  include_vep: false
  region_aggregation: sum
  weighted_regions:
    - name: altitude
      path: data/altitude_regions.bed
      weight: 25
    - name: immune
      path: data/immune_regions.bed
      weight: 12
```

Each `path` is a BED3/BED4 file using 0-based, half-open coordinates. A SNP
at 1-based position `POS` matches BED interval `[START0, END0)` when
`START0 < POS <= END0`. Names must be unique. Weights may be positive, zero,
or negative and are recorded in the run summary.

When one SNP matches more than one file, `region_aggregation` controls the
bonus:

| Value | Behavior |
|---|---|
| `sum` | Add every matching file's weight. |
| `max` | Use the numerically largest matching weight. |
| `first` | Use the first matching entry in YAML order. |

`include_maf: true` applies the MAF function in
`resources/scoring_weights.tsv`. Set it to `false` for coordinate weights
alone. To add VEP, set `include_vep: true` and provide `inputs.vep`; VEP
category weights still come from the same weight table.

The complete runnable template is
[`config/custom_regions.example.yaml`](../config/custom_regions.example.yaml).

## Study profile plus additional regions

The `study` profile always applies the original MAF, VEP, conservation, QTL,
literature, expression, and mutually exclusive selection policy. It may also
contain `weighted_regions`; their aggregated bonus is added to the original
selection bonus. This is useful for a sensitivity analysis, but it no longer
represents an exact reproduction of the published scoring policy.

## Provenance output

`score_filter_chip_summary.json` records the profile, MAF and VEP switches,
aggregation method, absolute path, weight, and interval count for every
custom track. `snp_filter_snp_scores.tsv` records the resulting category and
bonus for every SNP.
