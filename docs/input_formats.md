# Input Formats and Coordinate Conventions

All inputs must refer to the same reference assembly and use identical
chromosome labels. Mixing `1` with `ECA1`, or mixing EquCab3 coordinates with
the Mongolian horse assembly, is an error.

## PLINK binary genotype input

PLINK `.bed` is binary genotype data; it is not the tab-delimited BED interval
format used for windows and region tracks.

Source files may use rsIDs and an arbitrary allele order. Run
`tools/prepare_plink_input.py` with the matching reference before configuring
the workflow. It retains biallelic A/C/G/T SNPs and writes a new canonical
BED/BIM/FAM set; the source is never edited in place.

With `genotype_source: plink_bfile`, `inputs.bfile` is the prefix shared by
`PREFIX.bed`, `PREFIX.bim`, and `PREFIX.fam`. The BIM must use one
biallelic record per coordinate and canonical `CHROM:POS:REF:ALT` IDs. The
ID chromosome and 1-based position must match BIM columns 1 and 4; BIM A1
(column 5) must equal ALT and A2 (column 6) must equal REF. The workflow
validates these conditions before launching window jobs.

PLINK2 generates the AFREQ table. PLINK 1.9 applies `--show-tags all` in each
window with configurable MAF, tag distance, and r² thresholds. Empty windows
are retained as empty records rather than silently dropped.

## Windows

Tab-separated BED3 or BED4:

```text
ECA1    0    40000    ECA1-0
```

Coordinates are 0-based and half-open. If column 4 is absent, the workflow
creates the identifier `CHROM-START0`. IDs must be unique. The study used
40-kb windows. `ld_window_boundary: study_compatible` preserves its historical
`--from-bp START0 --to-bp END0` call; `bed` converts the left boundary to
`START0+1`. `60K` in the legacy filename describes the intended marker
count, not window width.

## Window LD tag sets

One tab-separated record per window, with no header:

```text
ECA1-0    ECA1:100:T:G|ECA1:150:C:T,ECA1:200:T:C|NONE
```

- Column 1 is exactly the window ID.
- Column 2 contains comma-separated PLINK LD tag sets.
- Each LD tag set starts with its target SNP; LD tag SNPs follow, separated by `|`.
- Every SNP token uses `CHROM:POS:REF:ALT`, with 1-based `POS`.
- `NONE` is ignored.

These records are PLINK 1.9 LD tag sets. The precomputed-input configuration
key is `inputs.ld_tag_sets`. With `plink_bfile`, the workflow creates
this table directly and records all LD parameters in a JSON summary.

## Variant Frequencies

PLINK2 `.afreq` output with at least these columns:

```text
#CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT
ECA1 ECA1:100:T:G T G Y 0.40 200
```

Requirements:

- SNP IDs must be `CHROM:POS:REF:ALT`.
- Records must be biallelic and `ALT_FREQS` must contain one numeric value.
- Only one ALT record may remain at each chromosome/position.
- Split multiallelic variants and normalize alleles before creating `.afreq`.

## VEP Output

Standard tab-separated VEP output is accepted, plain text or gzip-compressed.
The first seven columns must be:

```text
Uploaded_variation Location Allele Gene Feature Feature_type Consequence
```

The preferred `Uploaded_variation` is `CHROM_POS_REF/ALT`; canonical
`CHROM:POS:REF:ALT` is also accepted. A row contributes only when the
coordinate and VEP `Allele` match the expected ALT. Multiple transcripts and
multiple consequence terms are reduced to the maximum configured weight.
Variants without an exact match receive `vep_unknown` rather than zero.

## Position Annotation Tracks

The conservation, QTL, and literature inputs are position tables, despite
some historical filenames ending in `.bed`:

```text
CHROM    POS1    OPTIONAL_VALUE
```

`POS1` is 1-based. Each matching row adds that track's configured weight.
Duplicate rows therefore count more than once and should only be retained when
they represent intentionally distinct evidence records.

## Expression Annotation

The expected layout is:

```text
SNP_CHROM SNP_START0 SNP_END1 FEATURE_CHROM FEATURE_START0 FEATURE_END0 FEATURE_ID BREADTH INTENSITY
```

`SNP_END1` is used as the 1-based SNP coordinate. The last two columns are
numeric expression breadth and intensity scores and are added directly. A SNP
overlapping multiple features receives the sum across rows.

## Selection Intervals

Local adaptation, PUR-local adaptation, and domestication inputs are BED3 or
BED4, using 0-based half-open intervals. For overlapping categories the bonus
is mutually exclusive:

```text
PUR-local adaptation > other local adaptation > domestication
```

## Custom weighted intervals

Each `scoring.weighted_regions` entry points to BED3/BED4 intervals and assigns
one numeric weight to the entire file. Matching uses standard 0-based,
half-open BED semantics. See [custom scoring](custom_scoring.md) for overlap
aggregation and optional MAF/VEP components.

## Reference FASTA

Provide the exact FASTA used for all coordinates and its `samtools faidx`
index. The workflow checks the reference base against the REF allele while
preparing 241-bp flanks. Final probe coordinates are 1-based and inclusive;
final panel BED output is 0-based and half-open.
