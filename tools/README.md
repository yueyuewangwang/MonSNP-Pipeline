# Pre-processing Utilities

These standalone input-preparation and validation tools run before the main
Snakemake workflow. Run them from the repository root after activating
`environment.yaml`; none edits source input in place.

| Tool | Purpose | Main outputs |
|---|---|---|
| `prepare_plink_input.py` | Orient alleles against FASTA, retain biallelic A/C/G/T SNPs, assign canonical IDs | `.bed`, `.bim`, `.fam`, `.log`, `.prepare.json` |
| `make_windows.py` | Create a FASTA index when absent and generate BED design windows | `.fai`, BED4, summary JSON |
| `preflight.py` | Validate the complete configured input set before Snakemake | terminal report and optional JSON |

All three tools implement `--help`. Detailed formats and complete commands are
in [`docs/custom_data_quickstart.md`](../docs/custom_data_quickstart.md) and
[`docs/custom_data_quickstart.zh-CN.md`](../docs/custom_data_quickstart.zh-CN.md).
