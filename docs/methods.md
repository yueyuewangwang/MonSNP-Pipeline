# Scoring and Probe Selection Methods

## PLINK binary input and LD tag sets

With `genotype_source: plink_bfile`, PLINK2 computes allele frequencies from
`PREFIX.bed/.bim/.fam`. PLINK 1.9 then runs `--show-tags all` per design
window. The study-compatible defaults are MAF 0.05, tag radius 100 kb, r²
0.8, and historical START0/END0 boundary handling. The output is represented as LD tag sets. Each run records the executable,
input counts, parameters, empty-window count, and LD tag set count.

`target_chromosomes` defines the run scope. The workflow creates one frequency,
window, and LD shard per target chromosome; independent LD rules can run in
parallel before deterministic gather rules restore the standard AFREQ and LD
tag-set paths.

## Variant Universe and MAF

The biallelic PLINK2 `.afreq` records define the eligible variant universe.
For ALT frequency `p`, MAF is `min(p, 1-p)`. The default MAF score is:

```text
0                         MAF < 0.05
10 * MAF                  0.05 <= MAF <= 0.30
3                         MAF > 0.30
```

## Functional Annotation

VEP rows are matched to the expected chromosome, 1-based position, REF, and
ALT. For a variant with multiple transcripts or consequence terms, the
highest configured category weight is retained. Tied maximum categories are
reported together. The study weights were missense 9, any splice 9, UTR 8,
upstream/downstream 8, synonymous 4, intron 3, intergenic 2, and unknown or
unmatched 1.

## Composite SNP Score

The base score is the sum of:

```text
MAF + VEP + conserved_200m + conserved_100v + QTL
    + literature + expression_breadth + expression_intensity
```

Each position-track occurrence adds its configured score. Expression values
are already transformed upstream and are added as numeric values.

Selection bonuses are applied once per SNP with strict priority: PUR-local
adaptation 50, other local adaptation 30, then domestication 15. A SNP never
receives more than one selection bonus.

## Ranking Within LD Tag Sets

Within each window, duplicate coordinates inside one LD tag set are collapsed
to their best allele-specific score. The LD tag set score is the sum of its
unique SNP total scores. LD tag sets are ordered by:

1. descending LD tag set score;
2. descending number of scored SNPs;
3. original LD tag set order.

Within an LD tag set, SNPs are ordered by descending total score, descending
base score, ascending position, then SNP ID. Candidates are emitted from the
highest-ranked LD tag set first until `candidates_per_window` is reached.

## Custom weighted regions

The `custom` profile can rank SNPs using MAF, optional VEP, and any number of
weighted BED files without requiring the study-specific annotation tracks. For a
SNP matching multiple files, weights are summed, reduced to the maximum, or
resolved by YAML order according to `region_aggregation`. The `study` profile
retains the original formula and may optionally add custom-region bonuses.
See [custom scoring](custom_scoring.md).

## Probe Screening

The workflow extracts 121 bp to the left and 120 bp to the right of each
1-based SNP position, producing a 241-bp reference flank with the SNP at
0-based offset 120 when sequence bounds permit.

The archived OligoMiner release runs `blockParse` with fixed 120-nt probes,
20-80% GC, Tm 40-200 degrees C, and exclusion of six-base homopolymers.
Candidate probes are aligned to the same reference with Bowtie2 using
`--very-sensitive-local -k 2`; only unique alignments are retained. Jellyfish
then rejects probes whose 25-mers exceed the configured occurrence threshold
of 3. A passing probe must also span the target SNP.

Candidate flanks and OligoMiner jobs are partitioned by target chromosome.
Bowtie2 and Jellyfish reference assets are built once and shared read-only by
all chromosome jobs. The filtered BED and status tables are then gathered
before global per-window probe selection, preserving cross-window uniqueness.

The final assembler walks candidate ranks in order and chooses the first
strict probe for each window. Target coordinates are kept globally unique;
if an already selected target recurs in another window, the next strict
candidate is considered.

## Validation

Every final record must have a unique window and target site, a 120-nt A/C/G/T
probe, matching target/probe chromosomes, target coverage by the probe, and a
probe sequence identical to the forward reference interval. These are
computational checks, not assay-manufacturing or wet-lab validation.

The workflow exposes three target rules: `run_annotation`, `run_scoring`, and
`run_full`. Their dependency boundaries are designed for inspecting outputs
before either expensive stage rather than for changing scientific parameters.
