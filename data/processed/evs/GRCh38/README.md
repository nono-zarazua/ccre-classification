# GRCh38 enhancer-versus-silencer data

This directory contains the balanced enhancer-versus-silencer (EVS) benchmark.
Enhancers are the positive class and silencers are the negative class; no genomic
or shuffled background sequences are used in this task.

## Source cleaning

`src/get_whole_body_regulatory_sets.sh` deduplicates EVN and SVN cCRE IDs and
removes every silencer-overlapping interval from EVN before chromosome splitting.
The clean source collections contain 1,522,586 enhancers and 196,083 silencers;
chrY is excluded from the benchmark.

## Balancing and splits

`src/make_GRCh38_evs_training.py` retains all 196,083 silencers on chr1–22/X and,
with seed 42, samples exactly the same number of enhancers separately within each
chromosome. This prevents chromosome composition from identifying the class.

| Split | Chromosomes | Enhancers | Silencers | Total observations |
|---|---|---:|---:|---:|
| Validation | chr1 | 28,920 | 28,920 | 57,840 |
| Test | chr2 | 14,291 | 14,291 | 28,582 |
| Training | chr3–22, chrX | 152,872 | 152,872 | 305,744 |
| **Complete benchmark** | chr1–22, chrX | **196,083** | **196,083** | **392,166** |

## Layout

- `chroms/chromN/`: selected enhancer and silencer BED/FASTA files for each chromosome.
- `real/chromosomes/chromN/`: chromosome training shards using the common
  `train_pos`/`train_neg` interface expected by model scripts.
- `real/validation_*`: chr1 evaluation data.
- `real/test_*`: held-out chr2 data.
- `real/train_100_*`: complete chr3–22/X training pool.
- `evs_manifest.tsv`: per-chromosome inputs, outputs, counts, and split assignment.
- `evs_summary.json`: machine-readable design and totals.

The directory name `real/` is retained for compatibility with existing training
and reporting scripts. In EVS it denotes real biological enhancer and silencer
sequences, not regulatory loci versus generated genomic backgrounds.
