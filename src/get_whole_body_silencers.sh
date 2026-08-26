#!/bin/bash

set -euo pipefail

PROCESSED_DIR="data/processed/svn/GRCh38"
FILTERED_BED="${PROCESSED_DIR}/GRCh38_Silencer_ELS.bed"
TRIMMED_BED="${PROCESSED_DIR}/GRCh38_Silencer_ELS.trimmed.bed"
CHROMS_DIR="${PROCESSED_DIR}/chroms"

input_beds=(
    "data/raw/Cai-Fullwood-2021.Silencer-cCREs.bed"
    "data/raw/Huan-Ovcharenko-2019.Silencer-cCREs.bed"
    "data/raw/Jayavelu-Hawkins-2020.Silencer-cCREs.bed"
    "data/raw/Pang-Snyder-2020.Silencer-cCREs.bed"
    "data/raw/REST-Silencers.bed"
    "data/raw/STARR-Silencers.Robust.bed"
    "data/raw/STARR-Silencers.Stringent.bed"
)

mkdir -p "${PROCESSED_DIR}" "${CHROMS_DIR}"

for input_bed in "${input_beds[@]}"; do
    if [[ ! -s "${input_bed}" ]]; then
        echo "ERROR: required silencer BED is missing or empty: ${input_bed}" >&2
        exit 1
    fi
    if ! awk -F '\t' 'NF != 6 { exit 1 }' "${input_bed}"; then
        echo "ERROR: expected exactly six columns in ${input_bed}" >&2
        exit 1
    fi
done

filtered_tmp=$(mktemp "${FILTERED_BED}.tmp.XXXXXX")
trimmed_tmp=$(mktemp "${TRIMMED_BED}.tmp.XXXXXX")
trap 'rm -f "${filtered_tmp}" "${trimmed_tmp}"' EXIT

# REST-Enhancers.bed is deliberately absent from input_beds. Keep every pELS
# and dELS row from the silencer files; do not deduplicate across studies.
awk -F '\t' '$6 == "pELS" || $6 == "dELS"' "${input_beds[@]}" > "${filtered_tmp}"

# genNullSeqs uses the BED interval and identifier. Match the established
# enhancer input by retaining only chromosome, start, end, and cCRE ID.
awk 'BEGIN { FS=OFS="\t" } { print $1, $2, $3, $4 }' \
    "${filtered_tmp}" > "${trimmed_tmp}"

filtered_count=$(wc -l < "${filtered_tmp}")
trimmed_count=$(wc -l < "${trimmed_tmp}")

if [[ "${filtered_count}" -eq 0 ]]; then
    echo "ERROR: no pELS/dELS silencer records passed the filter" >&2
    exit 1
fi
if [[ "${filtered_count}" -ne "${trimmed_count}" ]]; then
    echo "ERROR: filtered and trimmed record counts differ" >&2
    exit 1
fi
if ! awk -F '\t' 'NF != 4 { exit 1 }' "${trimmed_tmp}"; then
    echo "ERROR: trimmed output does not contain exactly four columns" >&2
    exit 1
fi

mv "${filtered_tmp}" "${FILTERED_BED}"
mv "${trimmed_tmp}" "${TRIMMED_BED}"
trap - EXIT
# Split the trimmed four-column BED into one file per canonical chromosome.
# Empty chromosome files are created so the output layout is always complete.
chromosomes=({1..22} X Y)
for chromosome in "${chromosomes[@]}"; do
    chrom_dir="${CHROMS_DIR}/chrom${chromosome}"
    mkdir -p "${chrom_dir}"
    : > "${chrom_dir}/chrom${chromosome}_GRCh38_Silencer_ELS.bed"
done

awk -v chroms_dir="${CHROMS_DIR}" '
    BEGIN { OFS="\t" }
    {
        chromosome = $1
        sub(/^chr/, "", chromosome)
        output = chroms_dir "/chrom" chromosome "/chrom" chromosome "_GRCh38_Silencer_ELS.bed"
        print $0 > output
    }
' "${TRIMMED_BED}"
chrom_count=$(find "${CHROMS_DIR}" -name 'chrom*_GRCh38_Silencer_ELS.bed' -type f -exec cat {} + | wc -l)

if [[ "${trimmed_count}" -ne "${chrom_count}" ]]; then
    echo "ERROR: chromosome BEDs contain ${chrom_count} records; expected ${trimmed_count}" >&2
    exit 1
fi

echo "Filtered silencer pELS/dELS records by input:"
for input_bed in "${input_beds[@]}"; do
    count=$(awk -F '\t' '$6 == "pELS" || $6 == "dELS" { n++ } END { print n + 0 }' "${input_bed}")
    printf '%8d  %s\n' "${count}" "${input_bed}"
done
echo "Wrote ${filtered_count} six-column records to ${FILTERED_BED}"
echo "Wrote ${trimmed_count} four-column records to ${TRIMMED_BED}"
echo "Wrote ${chrom_count} records across ${#chromosomes[@]} chromosome BED files"
