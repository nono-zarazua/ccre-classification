# Chromosome-matched GRCh38 EVS source files

Each chromosome directory contains the selected enhancer and silencer loci and
their reference-genome sequences:

```text
chromN_GRCh38_enhancer.bed
chromN_GRCh38_enhancer.fasta
chromN_GRCh38_silencer.bed
chromN_GRCh38_silencer.fasta
```

All cleaned silencers are retained. Enhancers are sampled with seed 42 to equal
the silencer count independently on every chromosome. Enhancer and silencer BED
intervals are mutually non-overlapping, and BED/FASTA counts and sequence lengths
are validated by the builder. ChrY is excluded.

Exact per-chromosome counts are recorded in `chromosome_counts.csv` and
`../evs_manifest.tsv`. Chr1 is validation, chr2 is held-out test, and chr3–22/X
form the cumulative chromosome-learning training series.
