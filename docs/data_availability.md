# Data Availability

## Repository Policy

This GitHub repository contains source code, environment definitions,
documentation, scoring parameters, and a small synthetic test. It must not
contain individual-level horse genotype or phenotype files, multi-gigabyte
VEP/LD tag set tables, Bowtie2 indexes, SAM files, Conda environments, or other
regenerable intermediates.

## Study Data Release

The following fields are release blockers and must be completed before the
repository is made public:

| Item | Public location | Identifier | Status |
|---|---|---|---|
| Reference assembly FASTA and FAI | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| Population variant data or authorized derived AFREQ | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| Frozen window LD tag set table | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| VEP 113 output | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| Conservation/QTL/literature/expression tracks | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| Selection-signal BED files | `TO_BE_ADDED` | `TO_BE_ADDED` | pending |
| Manuscript reproduction bundle | `TO_BE_ADDED` | DOI/accession pending | pending |

Recommended split:

- GitHub: code and small test data;
- Zenodo or an institutional repository: derived inputs, final panel, run
  manifests, checksums, and release archive;
- a controlled-access genetics repository: individual-level genotypes when
  consent or policy prevents open distribution.

If genotype data are controlled, the derived non-individual LD tag set and
frequency files needed for the workflow should be made openly available when
ethically and legally permitted. Otherwise document the access committee,
eligibility, application URL, and expected response time.

## Integrity

The manuscript release bundle should retain the exact filenames and SHA-256
values in `resources/manuscript_checksums.tsv`. After depositing, add the
archive DOI and a download command here, then test that command from a clean
machine or CI runner.
