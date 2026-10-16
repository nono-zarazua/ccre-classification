#!/bin/bash

set -euo pipefail

mkdir -p data/raw data/processed data/reference/antinoise

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

# These two antinoise programs return 1 after successful completion. Run them
# without allowing set -e to abort, then validate their complete output sets.
run_return_one() {
    set +e
    "$@"
    status=$?
    set -e
    [[ "${status}" -eq 1 ]] || fail "$1 returned unexpected status ${status}"
}

if [[ ! -s data/raw/GRCh38-cCREs.bed ]]; then
    wget -O data/raw/GRCh38-cCREs.bed https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
fi
if [[ ! -s data/reference/hg38.fa ]]; then
    wget -O data/reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
    gunzip -f data/reference/hg38.fa.gz
fi

canonical=data/reference/antinoise/hg38-PLAIN-canonical
merged=data/reference/antinoise/hg38-PLAIN-canonical-merged
masked=data/reference/antinoise/hg38-PLAIN-canonical-masked
chroms=data/reference/antinoise/hg38-chroms

mkdir -p "${merged}" "${masked}" "${canonical}" "${chroms}"

run_return_one ../antinoise/bin/linux/fasta_muliplefiles.exe \
    data/reference/hg38.fa \
    "${chroms}/" \
    0

for chr in {1..22} X Y; do
    [[ -s "${chroms}/${chr}.fa" ]] || fail "splitter did not create ${chroms}/${chr}.fa"
    mv "${chroms}/${chr}.fa" "${canonical}/chr${chr}.fa"
done

run_return_one ../antinoise/bin/linux/fasta_to_plain0.exe \
    "${canonical}/" \
    hg38

for chr in {1..22} X Y; do
    [[ -s "${canonical}/chr${chr}.plain" ]] || fail "converter did not create ${canonical}/chr${chr}.plain"
done

# 1. Sort all cCRE annotations.
../antinoise/bin/linux/bed_sort.exe \
    data/raw/GRCh38-cCREs.bed \
    data/processed/GRCh38-cCREs_sorted.bed \
    hg38
[[ -s data/processed/GRCh38-cCREs_sorted.bed ]] || fail "sorted cCRE BED was not created"

# 2. Split the sorted BED into chromosome-specific BED files.
../antinoise/bin/linux/bed_chr_separation.exe \
    data/processed/GRCh38-cCREs_sorted.bed \
    "${merged}/GRCh38-cCREs" \
    hg38

# 3. Merge overlapping regions.
../antinoise/bin/linux/area_self_overlap.exe \
    "${merged}/GRCh38-cCREs" \
    "${merged}/GRCh38-cCREs_merged" \
    hg38 \
    "${merged}/overlap_stats.txt"

../antinoise/bin/linux/bed_chr_mask.exe \
    "${canonical}/" \
    "${merged}/" \
    GRCh38-cCREs_merged \
    chr \
    chr \
    -1 \
    hg38

for chr in {1..22} X Y; do
    chromosome_dir="${masked}/chr${chr}"
    chromosome_file="${chromosome_dir}/chr${chr}.plain"
    mkdir -p "${chromosome_dir}"
    mv "${merged}/chr${chr}.plain" "${chromosome_file}"
    rm -f "${masked}/chr${chr}.plain"
    [[ -s "${chromosome_file}" ]] || fail "missing masked chr${chr}.plain"
    plain_count=$(find "${chromosome_dir}" -maxdepth 1 -type f -name "*.plain" | wc -l)
    [[ "${plain_count}" -eq 1 ]] || fail "chr${chr} reference directory must contain exactly one PLAIN file"
done
