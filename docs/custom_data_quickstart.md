# Quick Start With Your Own Population

This path is for users who have a PLINK binary genotype set, the exact
reference assembly used for its coordinates, and one or more regions of
interest.

## 1. Required files

| Input | Required format |
|---|---|
| Genotypes | One prefix shared by `PREFIX.bed`, `PREFIX.bim`, and `PREFIX.fam` |
| Reference | Uncompressed FASTA; sequence names must match PLINK chromosome names |
| Weighted regions | BED3/BED4, 0-based, half-open; one numeric weight per file |
| Design windows | BED3/BED4, 0-based, half-open; generate these with `make_windows.py` |

PLINK `.bed` is a binary genotype file and is unrelated to the text BED3/BED4
interval format used for windows and weighted regions.

The source BIM may contain rsIDs and either allele order. The preparation tool
uses the FASTA with PLINK2, retains biallelic A/C/G/T SNPs, establishes REF and
ALT, sorts variants, and writes IDs as `CHROM:POS:REF:ALT`. It writes a new
fileset and never edits the source in place.

The FASTA, BIM, windows, and all region BEDs must use the same assembly,
chromosome labels, and coordinates.

## 2. Prepare the PLINK files

```bash
python tools/prepare_plink_input.py \
  --bfile data/my_population \
  --reference data/assembly.fa \
  --out data/prepared/my_population \
  --chr-set 31 --threads 8
```

When PLINK and FASTA chromosome names differ, provide a headerless two-column
`OLD_CHROM NEW_CHROM` file with `--chrom-map`; unlisted names are unchanged.

Set `--chr-set` to the autosome count expected by PLINK for the species (for
example, 22 for human and 31 for horse). The output is a new BED/BIM/FAM set
plus `my_population.prepare.json`. Indels, multiallelic sites, and non-ACGT
alleles are intentionally excluded from this SNP-array workflow.

## 3. Generate design windows

```bash
python tools/make_windows.py \
  --reference data/assembly.fa \
  --output data/design_windows.bed \
  --window-size 40000 --step 40000 \
  --exclude-regex '^(MT|chrM)$'
```

If `assembly.fa.fai` is absent, the tool creates it. By default it keeps the
short final window of each contig. Use `--min-contig-length`,
`--min-window-size`, `--include-regex`, `--exclude-regex`, or repeated
`--exclude-contig` options when needed.

Design windows determine where candidate SNPs are ranked. Weighted regions
only add scores; they do not replace the design-window file.

## 4. Supply weighted regions

Each BED file may have three or four columns:

```text
chr1    99999    120000    altitude_region_1
chr1    450000   470000    altitude_region_2
```

The first base is 0-based and included; the end is excluded. Blank lines and
lines beginning with `#` are ignored. Configure one weight for the entire
file:

```yaml
scoring:
  profile: custom
  include_maf: true
  include_vep: false
  region_aggregation: sum
  weighted_regions:
    - name: altitude
      path: data/altitude.bed
      weight: 25
    - name: immune
      path: data/immune.bed
      weight: 12
```

Weights may be positive, zero, or negative. For a SNP matching multiple files,
`sum` adds all weights, `max` keeps the numerically largest, and `first` uses
the first matching YAML entry.

## 5. Select chromosomes, run preflight, and execute

Copy `config/custom_regions.example.yaml`, replace all paths and species
parameters, then run:

```yaml
# Recommended first real-data trial
target_chromosomes:
  - ECA1
```

Then validate and run through candidate scoring:

```bash
python tools/preflight.py \
  --config config/my_population.yaml --base-dir . \
  --report preflight.json

snakemake -s workflow/Snakefile \
  --configfile config/my_population.yaml \
  --cores 8 --software-deployment-method conda \
  run_scoring
```

Preflight checks the PLINK binary structure and counts, canonical variant IDs,
sampled REF alleles against FASTA, FAI bounds, BED semantics, chromosome-name
agreement, configured weights, and required files. Errors stop the run;
warnings identify conditions that require review.

Inspect `results/02_scoring/selection_ld_tag_set_candidates.tsv`. When ready,
set `target_chromosomes: all`, allocate suitable cores, and invoke `run_full`
to build shared Bowtie2/Jellyfish indexes and run chromosome-scattered
OligoMiner screening. See the root README for all three checkpoint commands.

## Demo

Run `bash demo/run_demo.sh` for a six-SNP example. See
[`demo/README.md`](../demo/README.md) for its inputs and expected behavior.
