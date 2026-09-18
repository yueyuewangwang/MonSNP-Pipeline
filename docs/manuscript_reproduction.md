# Version History and Reproducibility

## Why There Are Two Paths

The public de novo workflow applies VEP weights from the start and is the
recommended reusable implementation. During the study, an earlier functional
annotation score was first used, then replaced by VEP 113 while all other
score components were frozen. That correction changed the best SNP in 12,620
of 61,175 comparable windows. Only those windows were redesigned, while valid
strict probes outside them were retained.

The historical path was:

```text
legacy full SNP scores and window winners
  -> replace only annotation_score with exact VEP score
  -> verify zero mismatches when reproducing legacy winners
  -> identify 12,620 changed windows
  -> rerank up to 10 candidates in each changed window
  -> validate and reuse existing strict probes where possible
  -> run OligoMiner/Bowtie2/Jellyfish by rank for unresolved candidates
  -> assemble changed-window outcomes with unchanged old-panel records
  -> independently validate reference, BIM alleles, probes, windows, and sites
```

## Frozen Result

The recommended manuscript panel is
`revised_strict_probe_panel.tsv`, with SHA-256
`5779bdd23c862d9cb35765bad8b43b602a5d51f4ddf529c3bad0e880c693a8bc`.
It contains 58,063 unique strict probe records:

- 46,005 unchanged-window records retained from the previous strict panel;
- 12,058 resolved VEP-changed windows;
- 10,329 changed windows resolved at VEP rank 1;
- 1,729 changed windows resolved at ranks 2-10;
- 562 changed windows unresolved after rank 10.

The corresponding BED SHA-256 is
`d46466fc05369edb5eb908004ab662632248c837efaa34a7e0a6d4a6d01745e2`.
The independent validation JSON SHA-256 is
`a3f0ec1d0f3228e88abb7279f5f6640fd886985837b11b1cf26a3f59e2fd225f`.

## Required Archive Contents

The public reproduction bundle must contain the large inputs and staged
outputs listed in `resources/manuscript_checksums.tsv`. These files are not
appropriate for GitHub. Preserve the VEP comparison outputs, per-rank strict
probe screens, final assembly tables, input manifest, and checksums so a
reviewer can inspect every decision without rerunning all alignments.

For a fresh scientifically equivalent run, use `workflow/Snakefile` with the
same reference, LD tag sets, annotation tracks, VEP file, and weights. A fully
fresh probe run may select a different equally passing 120-nt sequence for a
target when multiple probes tie; the frozen archive remains the authority for
byte-level reproduction of the manuscript panel.
