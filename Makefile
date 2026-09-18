SNAKEMAKE ?= snakemake
CONFIG ?= config/config.yaml
CORES ?= 8

.PHONY: test preflight dry-run annotation scoring full run demo

test:
	bash tests/run_smoke.sh

preflight:
	python tools/preflight.py --config $(CONFIG) --base-dir . --report preflight.json

dry-run:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) --cores $(CORES) --dry-run

annotation:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) --cores $(CORES) run_annotation --software-deployment-method conda

scoring:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) --cores $(CORES) run_scoring --software-deployment-method conda

full:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) --cores $(CORES) run_full --software-deployment-method conda

demo:
	bash demo/run_demo.sh

run:
	$(SNAKEMAKE) -s workflow/Snakefile --configfile $(CONFIG) --cores $(CORES) --software-deployment-method conda
