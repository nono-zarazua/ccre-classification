# Clean GRCh38 EVN source pools by chromosome

These are the complete cleaned enhancer-positive source pools. They are inputs to the matched-size benchmark, not final training sets. EVN is downsampled with seed 42 to 196,083 loci while preserving chromosome proportions; final positives are retained only when a valid paired background is available.

| chromosome | clean source | benchmark foreground | real pairs | missing real | synthetic pairs | missing synthetic |
|---|---:|---:|---:|---:|---:|---:|
| chr1 | 124,276 | 16,024 | 15,889 | 135 | 16,024 | 0 |
| chr2 | 129,963 | 16,757 | 16,630 | 127 | 16,757 | 0 |
| chr3 | 104,013 | 13,411 | 13,307 | 104 | 13,411 | 0 |
| chr4 | 75,397 | 9,722 | 9,637 | 85 | 9,722 | 0 |
| chr5 | 86,971 | 11,214 | 11,157 | 57 | 11,214 | 0 |
| chr6 | 85,010 | 10,961 | 10,907 | 54 | 10,961 | 0 |
| chr7 | 82,690 | 10,662 | 10,566 | 96 | 10,662 | 0 |
| chr8 | 75,704 | 9,761 | 9,677 | 84 | 9,761 | 0 |
| chr9 | 63,734 | 8,218 | 8,157 | 61 | 8,218 | 0 |
| chr10 | 77,579 | 10,003 | 9,929 | 74 | 10,003 | 0 |
| chr11 | 73,477 | 9,474 | 9,420 | 54 | 9,474 | 0 |
| chr12 | 76,235 | 9,830 | 9,760 | 70 | 9,830 | 0 |
| chr13 | 46,418 | 5,985 | 5,933 | 52 | 5,985 | 0 |
| chr14 | 50,617 | 6,526 | 6,470 | 56 | 6,526 | 0 |
| chr15 | 46,358 | 5,977 | 5,924 | 53 | 5,977 | 0 |
| chr16 | 48,093 | 6,201 | 6,153 | 48 | 6,201 | 0 |
| chr17 | 62,278 | 8,030 | 7,955 | 75 | 8,029 | 1 |
| chr18 | 36,726 | 4,735 | 4,693 | 42 | 4,735 | 0 |
| chr19 | 44,617 | 5,753 | 5,710 | 43 | 5,753 | 0 |
| chr20 | 41,635 | 5,368 | 5,325 | 43 | 5,368 | 0 |
| chr21 | 19,179 | 2,473 | 2,462 | 11 | 2,473 | 0 |
| chr22 | 29,670 | 3,826 | 3,792 | 34 | 3,826 | 0 |
| chrX | 40,109 | 5,172 | 5,135 | 37 | 5,172 | 0 |
| **chr1–22/X** | **1,520,749** | **196,083** | **194,588** | **1,495** | **196,082** | **1** |
| chrY | 1,837 | excluded | — | — | — | — |

Source BEDs are `chrom<chrom>/chrom<chrom>_GRCh38_ELS.bed`. Quotas are recorded in `data/processed/consolidation/GRCh38/manifest.tsv`. Final pairs are under `../real/` and `../synthetic/`. Real negatives have BED and FASTA records; synthetic negatives have no genomic locus and therefore only FASTA records.
