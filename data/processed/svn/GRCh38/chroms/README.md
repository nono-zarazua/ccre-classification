# Clean GRCh38 SVN source pools by chromosome

These are the complete deduplicated silencer-positive source pools. All 196,083 loci enter the joint foreground pools; final positives are retained only when a valid paired background is available.

| chromosome | clean source and benchmark foreground | real pairs | missing real | synthetic pairs | missing synthetic |
|---|---:|---:|---:|---:|---:|
| chr1 | 28,920 | 28,623 | 297 | 28,920 | 0 |
| chr2 | 14,291 | 14,129 | 162 | 14,291 | 0 |
| chr3 | 9,583 | 9,473 | 110 | 9,583 | 0 |
| chr4 | 10,835 | 10,728 | 107 | 10,835 | 0 |
| chr5 | 10,498 | 10,350 | 148 | 10,498 | 0 |
| chr6 | 16,221 | 16,018 | 203 | 16,221 | 0 |
| chr7 | 8,265 | 8,188 | 77 | 8,265 | 0 |
| chr8 | 7,324 | 7,250 | 74 | 7,324 | 0 |
| chr9 | 10,466 | 10,351 | 115 | 10,466 | 0 |
| chr10 | 7,817 | 7,722 | 95 | 7,817 | 0 |
| chr11 | 12,053 | 11,904 | 149 | 12,053 | 0 |
| chr12 | 7,555 | 7,464 | 91 | 7,555 | 0 |
| chr13 | 1,427 | 1,414 | 13 | 1,427 | 0 |
| chr14 | 3,992 | 3,955 | 37 | 3,992 | 0 |
| chr15 | 8,639 | 8,559 | 80 | 8,639 | 0 |
| chr16 | 8,233 | 8,153 | 80 | 8,233 | 0 |
| chr17 | 7,648 | 7,573 | 75 | 7,647 | 1 |
| chr18 | 3,541 | 3,510 | 31 | 3,541 | 0 |
| chr19 | 2,979 | 2,950 | 29 | 2,979 | 0 |
| chr20 | 6,927 | 6,859 | 68 | 6,927 | 0 |
| chr21 | 2,632 | 2,596 | 36 | 2,631 | 1 |
| chr22 | 4,828 | 4,784 | 44 | 4,828 | 0 |
| chrX | 1,409 | 1,388 | 21 | 1,409 | 0 |
| **Total** | **196,083** | **193,941** | **2,142** | **196,081** | **2** |

Source BEDs are `chrom<chrom>/chrom<chrom>_GRCh38_Silencer_ELS.bed`; chrY is an empty placeholder. Final pairs are under `../real/` and `../synthetic/`. Real negatives have BED and FASTA records; synthetic negatives have no genomic coordinates and therefore only FASTA records.
