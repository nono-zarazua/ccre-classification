# GRCh38 silencer-versus-negative data

This directory contains the cleaned silencer-positive source collection and the final GRCh38 silencer-versus-negative (SVN) benchmark inputs.

## Positive-source cleaning

`src/get_whole_body_regulatory_sets.sh` consolidates seven silencer sources (excluding `REST-Enhancers.bed`), verifies repeated cCRE IDs, retains one record per ID, preserves distinct partially overlapping cCREs, removes the resulting SVN loci from SCREEN ELS, and splits both tasks by canonical chromosome.

| Dataset | Source rows | Duplicates removed | SVN loci removed from EVN | Clean loci |
|---|---:|---:|---:|---:|
| EVN | 1,718,669 | 0 | 196,083 | **1,522,586** |
| SVN | 202,493 | 6,410 | — | **196,083** |

The 6,410 repeated rows are cCREs reported by multiple silencer studies, not conflicting loci. SVN has 196,083 unique loci on chr1–22/X and an empty chrY placeholder.

## Joint benchmark

All SVN loci and a chromosome-proportional 196,083-locus EVN sample are combined per chromosome before candidate generation. A positive is retained only with an acceptable matched negative.

| Background | Validation (chr1) | Test (chr2) | Training (chr3–22/X) | Total per class |
|---|---:|---:|---:|---:|
| Real genomic | 28,623 | 14,129 | 151,189 | **193,941** |
| Synthetic Antinoise `mix0` | 28,920 | 14,291 | 152,870 | **196,081** |

Every final split is balanced 1:1.

## Layout and leakage controls

- `GRCh38_Silencer_ELS.bed`: deduplicated six-column records.
- `GRCh38_Silencer_ELS.trimmed.bed`: corresponding four-column loci.
- `chroms/`: complete cleaned chromosome source pools.
- `real/` and `synthetic/`: final paired benchmark data and chr3–22/X training shards.

Chr1 is validation and chr2 is held-out test. EVN and SVN positives are mutually exclusive. Real negatives overlapping any SCREEN cCRE are rejected and exact selected loci are not reused. Synthetic candidates matching a foreground or reverse complement are rejected, and selected identities are globally unique. Audit tables are under `data/processed/consolidation/GRCh38/antinoise*_selection/`.
