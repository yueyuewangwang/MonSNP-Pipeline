# Small Custom-Data Demo

This demo starts with six SNPs in four individuals, a 1-kb reference, a
`1`-to-`ECA1` chromosome-name map, and two weighted BED tracks. It
exercises PLINK conversion, REF/ALT normalization, window generation, input
preflight, LD tag-set generation, and custom scoring.

The readable PED/MAP files are converted to a real source BED/BIM/FAM set at
the start of the script; the public workflow itself then starts from that PLINK
binary prefix.

From the repository root:

```bash
conda env create -f environment.yaml
conda activate monsnp-pipeline
bash demo/run_demo.sh
```

The two tracks assign `+25` to the altitude interval and `+12` to the immune
interval. Because `region_aggregation: sum`, SNPs in both receive `+37`.
The script finishes by checking these values and printing `Demo assertions passed`.

The demo stops before OligoMiner probe screening. Delete `demo/work/` to reset
all generated files.
