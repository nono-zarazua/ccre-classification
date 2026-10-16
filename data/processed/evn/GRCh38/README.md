# GRCh38 enhancer-versus-negative data

This directory contains the cleaned enhancer-positive source collection and the final GRCh38 enhancer-versus-negative (EVN) benchmark inputs.

## Positive-source cleaning

`src/get_whole_body_regulatory_sets.sh` prepares EVN and SVN together: it deduplicates the silencer sources by cCRE ID, preserves distinct partially overlapping cCREs, removes every SCREEN ELS record overlapping SVN by at least one base, and splits the cleaned records by canonical chromosome.

| Dataset | Source rows | Duplicates removed | SVN loci removed from EVN | Clean loci |
|---|---:|---:|---:|---:|
| EVN | 1,718,669 | 0 | 196,083 | **1,522,586** |
| SVN | 202,493 | 6,410 | — | **196,083** |

EVN and SVN have no genomic intersection after cleaning. EVN contains 1,520,749 benchmark-eligible loci on chr1–22/X and 1,837 excluded chrY loci.

## Joint benchmark

With seed 42, `src/make_GRCh38_positive_consolidation.py` selects 196,083 EVN loci while preserving the original chromosome proportions with Hamilton allocation. These are combined with all 196,083 SVN loci chromosome by chromosome before real and synthetic candidate generation. A positive is retained only with an acceptable matched negative, so final real and synthetic positive membership differs slightly.

| Background | Validation (chr1) | Test (chr2) | Training (chr3–22/X) | Total per class |
|---|---:|---:|---:|---:|
| Real genomic | 15,889 | 16,630 | 162,069 | **194,588** |
| Synthetic Antinoise `mix0` | 16,024 | 16,757 | 163,301 | **196,082** |

Every final split is balanced 1:1.

## Layout and leakage controls

- `GRCh38_ELS.trimmed.bed`: 1,522,586 cleaned EVN loci.
- `GRCh38_ELS.silencer_overlaps.bed`: ELS records excluded due to SVN overlap.
- `chroms/`: complete cleaned chromosome source pools.
- `real/` and `synthetic/`: final paired benchmark data and chr3–22/X training shards.
- `full-benchmark/`: separate full-EVN LS-GKM supervisor handoff.

Chr1 is validation, chr2 is held-out test, and chr3–22/X are cumulatively added during training. Real negatives overlapping any SCREEN cCRE are rejected and exact selected loci are not reused. Synthetic candidates identical to any foreground or reverse complement are rejected, and selected synthetic identities are globally unique. Audit tables are under `data/processed/consolidation/GRCh38/antinoise*_selection/`.
