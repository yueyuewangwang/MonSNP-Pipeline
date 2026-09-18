# GitHub and Journal Release Checklist

## Required Before Public Release

- Replace every `TO_BE_ADDED` entry in `docs/data_availability.md`.
- Update `CITATION.cff` with all authors, affiliations, repository URL,
  article title, journal, DOI, and version.
- Confirm the MIT license with all authors and the host institution.
- Deposit the manuscript reproduction bundle and verify every checksum.
- Confirm that no nonsynthetic `.fam`, phenotype, sample ID, access token, absolute server
  path, Conda environment, SAM file, or restricted dataset is tracked.
- Run `bash tests/run_smoke.sh` from a clean clone.
- Run the Snakemake dry run and core test configuration from a clean clone.
- Run `git lfs ls-files` and confirm that all BED/BIM/FAM/VCF.GZ and other
  configured large formats are LFS objects rather than regular Git blobs.
- Tag the submitted version, for example `v1.0.0`.
- Create a GitHub release from that tag and archive it with Zenodo for a DOI.
- Cite the immutable release tag/DOI in the manuscript, not only the moving
  default branch.

## Suggested Commands

```bash
git init
git lfs install
# For files already committed before .gitattributes was updated:
git lfs migrate import --include="*.bed,*.bim,*.fam,*.vcf.gz,*.bcf,*.pgen,*.bam,*.cram,*.jf"
git add .
git status --short
git grep -nE '/data[0-9]+/|BEGIN (RSA|OPENSSH) PRIVATE KEY|password[[:space:]]*='
bash tests/run_smoke.sh
snakemake -s workflow/Snakefile --configfile tests/config.yaml --cores 1 --dry-run
snakemake -s workflow/Snakefile --configfile tests/config.yaml --cores 1 run_full --dry-run
git lfs ls-files
```

Inspect `git status` and the grep output manually before committing. The grep
may legitimately match documentation that explains prohibited absolute paths;
it must not match executable defaults or credentials.

## Suggested GitHub Metadata

- Description: `MonSNP-Pipeline: reproducible LD-aware SNP prioritization and probe design for SNP arrays`
- Topics: `bioinformatics`, `snakemake`, `snp-array`, `horse-genomics`,
  `probe-design`, `reproducible-research`
- Protect the default branch and require the smoke-test workflow.
- Enable Zenodo only after the repository metadata and author list are final.
