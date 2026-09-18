# Changelog

## 0.3.0 - 2026-09-18

- Renamed the public project to MonSNP-Pipeline and rewrote the bilingual entry-point documentation.
- Added `target_chromosomes` with chromosome-scattered LD and OligoMiner jobs and deterministic gather outputs.
- Added stable `run_annotation`, `run_scoring`, and `run_full` stage targets.
- Added a final plain-text summary report and Git LFS tracking for large bioinformatics formats.
- Reframed testing fixtures and documentation with public-facing bioinformatics terminology.

## 0.2.0 - 2026-09-17

- Added PLINK normalization, automatic window generation, and comprehensive input preflight tools.
- Added a runnable six-SNP custom-region demo and bilingual own-data quick starts.
- Added direct PLINK BED/BIM/FAM input with PLINK2 allele frequencies.
- Added study-compatible per-window PLINK 1.9 LD tag set generation.
- Added study and custom scoring profiles with weighted BED tracks.
- Added sum, max, and first overlap aggregation with full JSON provenance.
- Renamed configuration keys, CLI arguments, outputs, table columns, and summaries to use LD tag set terminology consistently.

## 0.1.0 - 2026-09-16

- Created a portable Snakemake workflow from the study analysis scripts.
- Replaced absolute server paths and hard-coded scoring weights with config files.
- Added direct VEP-to-score conversion with exact allele matching.
- Added de novo all-window candidate ranking, strict probe assembly, and validation.
- Preserved OligoMiner commit `fcb124f` and its Python 2.7 runtime definition.
- Added synthetic smoke data, CI, bilingual usage documentation, and manuscript checksums.
