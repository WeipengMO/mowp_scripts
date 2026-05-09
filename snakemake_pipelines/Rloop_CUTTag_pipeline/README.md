# R-loop CUT&Tag Snakemake pipeline

This pipeline calls candidate R-loop peaks from untreated S9.6/HBD CUT&Tag libraries, then filters high-confidence R-loop peaks by enrichment over IgG/no-primary controls and depletion after RNase H treatment.

Run:

```bash
snakemake --use-conda -np
snakemake --use-conda -j 40
```
