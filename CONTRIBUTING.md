# Contributing

Please open an issue before changing scoring policy, coordinate conventions,
or the definition of a strict probe. Bug fixes should include a minimal test
that fails before the change and passes after it.

Run these checks before submitting a pull request:

```bash
bash tests/run_smoke.sh
snakemake -s workflow/Snakefile --configfile tests/config_multichrom.yaml --cores 2 run_scoring
python3 tests/assert_multichrom.py tests/work/multichrom_results
python3 -m py_compile workflow/scripts/*.py tools/*.py demo/assert_demo.py
snakemake -s workflow/Snakefile --configfile tests/config.yaml --cores 1 --dry-run
snakemake -s workflow/Snakefile --configfile tests/config.yaml --cores 1 run_full --dry-run
```

Do not commit individual-level data, restricted accessions, machine-specific
absolute paths, environments, or large intermediate files.
