#!/bin/bash

set -euo pipefail

RAW_BED="data/raw/GRCh38_ELS.bed"
PROCESSED_DIR="data/processed/evn/GRCh38"
TRIMMED_BED="${PROCESSED_DIR}/GRCh38_ELS.trimmed.bed"
CHROMS_DIR="${PROCESSED_DIR}/chroms"
DOWNLOAD_URL="https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed"

mkdir -p "data/raw" "${PROCESSED_DIR}" "${CHROMS_DIR}"

download_tmp=$(mktemp "${RAW_BED}.download.XXXXXX")
trap 'rm -f "${download_tmp}"' EXIT
wget -qO "${download_tmp}" "${DOWNLOAD_URL}"
mv "${download_tmp}" "${RAW_BED}"
trap - EXIT

awk 'BEGIN { FS="\t" } { print $1 }' "${RAW_BED}" | sort | uniq -c | sort -rn

awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4 }'     "${RAW_BED}" > "${TRIMMED_BED}"

chromosomes=({1..22} X Y)
for chromosome in "${chromosomes[@]}"; do
    chrom_dir="${CHROMS_DIR}/chrom${chromosome}"
    mkdir -p "${chrom_dir}"
    : > "${chrom_dir}/chrom${chromosome}_GRCh38_ELS.bed"
done

awk -v chroms_dir="${CHROMS_DIR}" '
    BEGIN { OFS="\t" }
    {
        chromosome = $1
        sub(/^chr/, "", chromosome)
        output = chroms_dir "/chrom" chromosome "/chrom" chromosome "_GRCh38_ELS.bed"
        print $0 > output
    }
' "${TRIMMED_BED}"

trimmed_count=$(wc -l < "${TRIMMED_BED}")
chrom_count=$(find "${CHROMS_DIR}" -name 'chrom*_GRCh38_ELS.bed' -type f -exec cat {} + | wc -l)

if [[ "${trimmed_count}" -ne "${chrom_count}" ]]; then
    echo "ERROR: chromosome BEDs contain ${chrom_count} records; expected ${trimmed_count}" >&2
    exit 1
fi

echo "Wrote ${trimmed_count} ELS records across ${#chromosomes[@]} chromosome BED files"
