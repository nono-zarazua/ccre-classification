#!/bin/bash

wget -qO data/raw/hepg2_cCREs.bed https://downloads.wenglab.org/Registry-V4/ENCFF546MZK_ENCFF732PJK_ENCFF795ONN_ENCFF357NFO.bed

awk -F'\t' '$10 == "dELS" || $10 == "pELS"' data/raw/hepg2_cCREs.bed > data/processed/hepg2_ELS.bed

awk 'BEGIN { FS="\t"} { print $10 }' data/processed/hepg2_ELS.bed | sort | uniq -c | sort -rn

echo "Keeping only first 4 columns."

awk 'BEGIN { OFS="\t"} { print $1,$2,$3,$4 }' data/processed/hepg2_ELS.bed > data/processed/hepg2_ELS.trimmed.bed
