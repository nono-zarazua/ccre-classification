# HepG2 enhancer fold composition

This legacy HepG2 benchmark groups chromosome-held loci into five approximately balanced folds. Each fold contains positive and 1x negative BED/FASTA pools. Missing negatives are the positive count minus the available negative count.

| Fold | Chromosomes | Positive ELS loci | Available negatives | Missing negatives |
|---|---|---:|---:|---:|
| 1 | chr1, chr20, chr14, chr22 | 9,744 | 9,727 | 17 |
| 2 | chr2, chr12, chr9, chr13 | 9,583 | 9,573 | 10 |
| 3 | chr17, chr11, chr7, chr15, chrX | 9,995 | 9,975 | 20 |
| 4 | chr6, chr19, chr5, chr4, chr21 | 9,874 | 9,859 | 15 |
| 5 | chr3, chr16, chr10, chr8, chr18 | 10,229 | 10,214 | 15 |
| **Total** | chr1–22/X | **49,425** | **49,348** | **77** |

This fold-based HepG2 workflow is independent of the current GRCh38 whole-body EVN/SVN chromosome-learning benchmark.
