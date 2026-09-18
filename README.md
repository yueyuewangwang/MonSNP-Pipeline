# MonSNP-Pipeline: A Reproducible Workflow for SNP Array Design

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)

[中文说明](README.zh-CN.md)

MonSNP-Pipeline turns WGS/PLINK data into ranked SNP candidates and a strict
probe panel for SNP-array design. You provide genotypes, a matching reference
assembly, design windows, and the annotations or regions you want to use. The
workflow handles allele frequencies, LD tag sets, variant scoring, candidate
ranking, and 120-nt probe screening.

Use it to design an array for your own population or species, or to reproduce
the Mongolian horse 58K design with the separate manuscript data bundle. A
six-SNP toy dataset is included so you can check the installation before using
study-scale data.

## Highlights

- **Genotypes to probe panel:** start from PLINK BED/BIM/FAM and automate
  candidate selection and probe design.
- **Configurable scoring:** combine MAF, VEP, conservation, QTL, expression,
  selection signals, or your own weighted regions.
- **Start with one chromosome:** validate real inputs on a small scope before
  switching to the whole genome.
- **Review results between stages:** stop before or after LD calculation and
  probe screening.
- **Reproducible execution:** manage jobs and software versions with Snakemake
  and Conda.

## What You Get

After a complete run, the main deliverables are in
`<results_dir>/06_final_panel/`. Here, `<results_dir>` is the `results_dir`
value in your configuration (the main template uses `results`):

| File | Purpose |
|---|---|
| `final_strict_probe_panel.tsv` | Final probes with variant coordinates, sequences, ranks, and scores |
| `final_strict_probe_panel.bed` | Coordinates of the selected variants |
| `summary_report.txt` | Counts retained at the main filtering stages |
| `validation.json` | Machine-readable sequence, coordinate, and uniqueness checks |

When you stop after scoring, inspect
`<results_dir>/02_scoring/selection_ld_tag_set_candidates.tsv` before starting probe
design.

## Check the Installation First

Requirements: Linux x86-64, Conda/Miniforge, and Git LFS.

```bash
git lfs install
git clone https://github.com/YourName/MonSNP-Pipeline.git
cd MonSNP-Pipeline
conda env create -f environment.yaml
conda activate monsnp-pipeline
bash demo/run_demo.sh
```

The installation check is complete when the command prints `Demo assertions
passed`. Demo files are written to `demo/work/` and the demo intentionally stops
before the expensive OligoMiner screen. See [`demo/README.md`](demo/README.md)
for the small input set and expected behavior.

## Run Your Own Data

### 1. Prepare the Inputs

You need at least:

| Input | Description |
|---|---|
| `PREFIX.bed/.bim/.fam` | PLINK binary genotypes sharing one prefix |
| `assembly.fa` | Reference assembly matching the genotype coordinates |
| `assembly.fa.fai` | FASTA index; the preparation tools can create it |
| `design_windows.bed` | 0-based, half-open design windows |

The study profile additionally uses VEP, conservation, QTL, expression, and
selection tracks. The custom profile lets you use your own weighted regions
without the manuscript-specific annotation tracks.

### 2. Normalize PLINK and Create Windows

The root-level `tools/` directory contains independent pre-processing utilities.
They run before Snakemake and do not modify source inputs in place.

```bash
python tools/prepare_plink_input.py \
  --bfile data/raw \
  --reference data/assembly.fa \
  --out data/prepared/population \
  --chr-set 31 \
  --threads 8

python tools/make_windows.py \
  --reference data/assembly.fa \
  --output data/design_windows.bed \
  --window-size 40000 \
  --step 40000
```

Use the chromosome count expected by PLINK for your species (`31` is typical
for horse; human data usually use `22`). See
[`docs/input_formats.md`](docs/input_formats.md) for format and coordinate rules.

### 3. Copy a Configuration Template

For a new population with custom regions, start here:

```bash
cp config/custom_regions.example.yaml config/my_run.yaml
```

Edit `inputs.bfile`, `inputs.reference`, `inputs.reference_fai`, `inputs.windows`,
and the paths and weights under `scoring.weighted_regions`. Then validate before
launching Snakemake:

```bash
python tools/preflight.py \
  --config config/my_run.yaml \
  --base-dir . \
  --report preflight.json
```

Continue when preflight reports `PASS` or `PASS_WITH_WARNINGS`. The configuration
reference is in [`config/README.md`](config/README.md); scoring options are in
[`docs/custom_scoring.md`](docs/custom_scoring.md).

The custom template uses `results_dir: results_custom`. Change it to
`results_dir: results` if you want the output paths shown in the main template.

## Start Small, Then Scale Up

For a first real-data run, select one chromosome in `config/my_run.yaml`:

```yaml
target_chromosomes:
  - ECA1
```

After checking the candidate table, change it to:

```yaml
target_chromosomes: all
```

### Three Run Stages

| Stage | Command target | Produces | Use it when |
|---|---|---|---|
| A: annotation preparation | `run_annotation` | Windows, frequencies, and optional VEP scores | Checking inputs for the first time |
| B: candidate ranking | `run_scoring` | LD tag sets, scores, and ranked candidates | Reviewing results before probe design |
| C: complete design | `run_full` | OligoMiner results and the final panel | Ready for the full computation |

```bash
# Stage A: quick input and annotation check
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 8 \
  --software-deployment-method conda run_annotation

# Stage B: LD calculation and candidate ranking
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 16 \
  --software-deployment-method conda run_scoring

# Stage C: full probe design and panel assembly
snakemake -s workflow/Snakefile \
  --configfile config/my_run.yaml --cores 32 \
  --software-deployment-method conda run_full
```

> **Time and resource warning:** LD tag-set calculation and OligoMiner probe
> screening are the two most expensive parts of the workflow. They can require
> substantial CPU time, memory, and disk space on real WGS data. Test one
> chromosome first, inspect the intermediate results, and run `--dry-run` before
> a whole-genome job. Set `--cores` to the resources available on your host or
> scheduler.

## Reproduce the Manuscript Result

Download the complete manuscript input bundle from the paper data archive,
place it under `data/manuscript/`, and update `config/config.yaml` as described
in [`docs/manuscript_reproduction.md`](docs/manuscript_reproduction.md). Then run:

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml --cores 32 \
  --software-deployment-method conda run_full
```

Compare the final panel size and SHA-256 values with
`resources/manuscript_checksums.tsv`. Replace the placeholder DOI, archive
accession, repository URL, and citation metadata before public release.

## More Documentation

- [Input formats and coordinate conventions](docs/input_formats.md)
- [Quick start with custom data](docs/custom_data_quickstart.md)
- [Custom scoring](docs/custom_scoring.md)
- [Methods and algorithms](docs/methods.md)
- [Version history and reproducibility](docs/manuscript_reproduction.md)
- [Data availability](docs/data_availability.md)
- [Release checklist](docs/release_checklist.md)

## Citation and License

If you use MonSNP-Pipeline, please cite:

> 

Update `CITATION.cff` before release. MonSNP-Pipeline is distributed under the
[MIT License](LICENSE); the vendored OligoMiner files retain their original
license as documented in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
