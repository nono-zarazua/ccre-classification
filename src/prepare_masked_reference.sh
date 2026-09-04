#!/bin/bash

set -euo pipefail

mkdir -p data/raw data/processed data/reference

wget -O data/raw/GRCh38-cCREs.bed https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
wget -O data/reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip -f data/reference/hg38.fa.gz

mkdir -p data/reference/hg38-PLAIN-canonical-merged
mkdir -p data/reference/hg38-PLAIN-canonical-masked
mkdir -p data/reference/hg38-PLAIN-canonical
mkdir -p data/reference/hg38-chroms

../antinoise/bin/linux/fasta_muliplefiles.exe \
    data/reference/hg38.fa \
    data/reference/hg38-chroms/ \
    0

for chr in {1..22} X Y; do
    mv data/reference/hg38-chroms/${chr}.fa data/reference/hg38-PLAIN-canonical/chr${chr}.fa
done

../antinoise/bin/linux/fasta_to_plain0.exe \
    data/reference/hg38-PLAIN-canonical/ \
    hg38

# 1. Sort all cCRE annotations
../antinoise/bin/linux/bed_sort.exe \
    data/raw/GRCh38-cCREs.bed \
    data/processed/GRCh38-cCREs_sorted.bed \
    hg38

# 2. Split the sorted BED into chromosome-specific BED files
../antinoise/bin/linux/bed_chr_separation.exe \
    data/processed/GRCh38-cCREs_sorted.bed \
    data/reference/hg38-PLAIN-canonical-merged/GRCh38-cCREs \
    hg38

# 3. Merge overlapping regions
../antinoise/bin/linux/area_self_overlap.exe \
    data/reference/hg38-PLAIN-canonical-merged/GRCh38-cCREs \
    data/reference/hg38-PLAIN-canonical-merged/GRCh38-cCREs_merged \
    hg38 \
    data/reference/hg38-PLAIN-canonical-merged/overlap_stats.txt


../antinoise/bin/linux/bed_chr_mask.exe \
    data/reference/hg38-PLAIN-canonical/ \
    data/reference/hg38-PLAIN-canonical-merged/ \
    GRCh38-cCREs_merged \
    chr \
    chr \
    -1 \
    hg38

for chr in {1..22} X Y; do
    mv data/reference/hg38-PLAIN-canonical-merged/chr${chr}.plain data/reference/hg38-PLAIN-canonical-masked/
done