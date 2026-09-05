#!/bin/bash

# Build mutually exclusive GRCh38 enhancer (EVN) and silencer (SVN) positive sets.
set -euo pipefail
export LC_ALL=C

ELS_URL="https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed"
SILENCER_URL="https://downloads.wenglab.org/Human_Silencers.tar.gz"
BEDTOOLS_BIN="${BEDTOOLS_BIN:-/home/zarazuanav/.conda/envs/bstool/bin/bedtools}"
RAW_DIR="data/raw"
EVN_DIR="data/processed/evn/GRCh38"
SVN_DIR="data/processed/svn/GRCh38"
EVN_CHROMS_DIR="${EVN_DIR}/chroms"
SVN_CHROMS_DIR="${SVN_DIR}/chroms"

silencer_files=(
  "Cai-Fullwood-2021.Silencer-cCREs.bed"
  "Huan-Ovcharenko-2019.Silencer-cCREs.bed"
  "Jayavelu-Hawkins-2020.Silencer-cCREs.bed"
  "Pang-Snyder-2020.Silencer-cCREs.bed"
  "REST-Silencers.bed"
  "STARR-Silencers.Robust.bed"
  "STARR-Silencers.Stringent.bed"
)
archive_files=("${silencer_files[@]}" "REST-Enhancers.bed")
chromosomes=({1..22} X Y)

[[ -x "${BEDTOOLS_BIN}" ]] || {
  echo "ERROR: bedtools not found: ${BEDTOOLS_BIN}; set BEDTOOLS_BIN if installed elsewhere" >&2
  exit 1
}
mkdir -p "${RAW_DIR}" "${EVN_DIR}" "${SVN_DIR}" "${EVN_CHROMS_DIR}" "${SVN_CHROMS_DIR}"
work_dir=$(mktemp -d "data/.prepare_GRCh38_regulatory_sets.XXXXXX")
trap 'rm -rf "${work_dir}"' EXIT
mkdir -p "${work_dir}/raw" "${work_dir}/evn_chroms" "${work_dir}/svn_chroms"

# Download to temporary paths so an interrupted download cannot replace valid data.
wget -qO "${work_dir}/raw/GRCh38_ELS.bed" "${ELS_URL}"
wget -qO "${work_dir}/raw/Human_Silencers.tar.gz" "${SILENCER_URL}"
tar -xzf "${work_dir}/raw/Human_Silencers.tar.gz" -C "${work_dir}/raw"

els_raw="${work_dir}/raw/GRCh38_ELS.bed"
[[ -s "${els_raw}" ]] || { echo "ERROR: downloaded SCREEN ELS BED is empty" >&2; exit 1; }
awk -F '\t' 'NF != 6 { exit 1 }' "${els_raw}" || {
  echo "ERROR: expected six columns in SCREEN ELS BED" >&2; exit 1;
}

silencer_inputs=()
for filename in "${silencer_files[@]}"; do
  input="${work_dir}/raw/${filename}"
  [[ -s "${input}" ]] || { echo "ERROR: missing archive member: ${filename}" >&2; exit 1; }
  awk -F '\t' 'NF != 6 { exit 1 }' "${input}" || {
    echo "ERROR: expected six columns in ${filename}" >&2; exit 1;
  }
  silencer_inputs+=("${input}")
done

# REST-Enhancers is intentionally excluded. The same cCRE can occur in multiple
# silencer studies. Verify repeated IDs, retain one record per cCRE ID, and sort
# unique loci without merging distinct overlapping cCREs.
svn_all="${work_dir}/GRCh38_Silencer_ELS.all_records.bed"
svn_unique_unsorted="${work_dir}/GRCh38_Silencer_ELS.unique.unsorted.bed"
svn_six="${work_dir}/GRCh38_Silencer_ELS.bed"
svn_four="${work_dir}/GRCh38_Silencer_ELS.trimmed.bed"
awk -F '\t' '$6 == "pELS" || $6 == "dELS"' "${silencer_inputs[@]}" > "${svn_all}"
awk -F '\t' '
  { current=$1 FS $2 FS $3 FS $4 FS $5 FS $6 }
  ($4 in first) && first[$4] != current {
    printf "ERROR: cCRE ID %s has conflicting silencer records\n", $4 > "/dev/stderr"
    conflicts++
    next
  }
  { if (!($4 in first)) { first[$4]=current; print } }
  END { if (conflicts) exit 1 }
' "${svn_all}" > "${svn_unique_unsorted}"
"${BEDTOOLS_BIN}" sort -i "${svn_unique_unsorted}" > "${svn_six}"
awk 'BEGIN { FS=OFS="\t" } { print $1,$2,$3,$4 }' "${svn_six}" > "${svn_four}"

els_unique_unsorted="${work_dir}/GRCh38_ELS.unique.unsorted.bed"
els_unique_six="${work_dir}/GRCh38_ELS.unique.bed"
els_all="${work_dir}/GRCh38_ELS.all.trimmed.bed"
evn_clean="${work_dir}/GRCh38_ELS.trimmed.bed"
evn_removed="${work_dir}/GRCh38_ELS.silencer_overlaps.bed"
awk -F '\t' '
  { current=$1 FS $2 FS $3 FS $4 FS $5 FS $6 }
  ($4 in first) && first[$4] != current {
    printf "ERROR: cCRE ID %s has conflicting enhancer records\n", $4 > "/dev/stderr"
    conflicts++
    next
  }
  { if (!($4 in first)) { first[$4]=current; print } }
  END { if (conflicts) exit 1 }
' "${els_raw}" > "${els_unique_unsorted}"
"${BEDTOOLS_BIN}" sort -i "${els_unique_unsorted}" > "${els_unique_six}"
awk 'BEGIN { FS=OFS="\t" } { print $1,$2,$3,$4 }' "${els_unique_six}" > "${els_all}"
"${BEDTOOLS_BIN}" intersect -v -a "${els_all}" -b "${svn_four}" > "${evn_clean}"
"${BEDTOOLS_BIN}" intersect -u -a "${els_all}" -b "${svn_four}" > "${evn_removed}"

enhancer_source_count=$(wc -l < "${els_raw}")
source_count=$(wc -l < "${els_all}")
enhancer_duplicate_count=$((enhancer_source_count - source_count))
clean_count=$(wc -l < "${evn_clean}")
removed_count=$(wc -l < "${evn_removed}")
silencer_source_count=$(wc -l < "${svn_all}")
silencer_count=$(wc -l < "${svn_four}")
silencer_duplicate_count=$((silencer_source_count - silencer_count))
[[ "${silencer_count}" -gt 0 ]] || { echo "ERROR: no silencer pELS/dELS records found" >&2; exit 1; }
[[ $((clean_count + removed_count)) -eq "${source_count}" ]] || {
  echo "ERROR: clean + removed enhancer counts do not equal source count" >&2; exit 1;
}
remaining_overlap=$("${BEDTOOLS_BIN}" intersect -u -a "${evn_clean}" -b "${svn_four}" | wc -l)
[[ "${remaining_overlap}" -eq 0 ]] || {
  echo "ERROR: ${remaining_overlap} retained enhancers still overlap silencers" >&2; exit 1;
}

# Build complete and predictable chromosome layouts in the temporary workspace.
for chromosome in "${chromosomes[@]}"; do
  mkdir -p "${work_dir}/evn_chroms/chrom${chromosome}" "${work_dir}/svn_chroms/chrom${chromosome}"
  : > "${work_dir}/evn_chroms/chrom${chromosome}/chrom${chromosome}_GRCh38_ELS.bed"
  : > "${work_dir}/svn_chroms/chrom${chromosome}/chrom${chromosome}_GRCh38_Silencer_ELS.bed"
done
awk -v base="${work_dir}/evn_chroms" 'BEGIN{OFS="\t"}{c=$1;sub(/^chr/,"",c);print > base "/chrom" c "/chrom" c "_GRCh38_ELS.bed"}' "${evn_clean}"
awk -v base="${work_dir}/svn_chroms" 'BEGIN{OFS="\t"}{c=$1;sub(/^chr/,"",c);print > base "/chrom" c "/chrom" c "_GRCh38_Silencer_ELS.bed"}' "${svn_four}"

evn_chrom_count=$(find "${work_dir}/evn_chroms" -name 'chrom*_GRCh38_ELS.bed' -type f -exec cat {} + | wc -l)
svn_chrom_count=$(find "${work_dir}/svn_chroms" -name 'chrom*_GRCh38_Silencer_ELS.bed' -type f -exec cat {} + | wc -l)
[[ "${evn_chrom_count}" -eq "${clean_count}" ]] || { echo "ERROR: EVN chromosome count mismatch" >&2; exit 1; }
[[ "${svn_chrom_count}" -eq "${silencer_count}" ]] || { echo "ERROR: SVN chromosome count mismatch" >&2; exit 1; }

# Install only after every validation has passed.
mv "${els_raw}" "${RAW_DIR}/GRCh38_ELS.bed"
mv "${work_dir}/raw/Human_Silencers.tar.gz" "${RAW_DIR}/Human_Silencers.tar.gz"
for filename in "${archive_files[@]}"; do mv "${work_dir}/raw/${filename}" "${RAW_DIR}/${filename}"; done
mv "${svn_six}" "${SVN_DIR}/GRCh38_Silencer_ELS.bed"
mv "${svn_four}" "${SVN_DIR}/GRCh38_Silencer_ELS.trimmed.bed"
mv "${evn_clean}" "${EVN_DIR}/GRCh38_ELS.trimmed.bed"
mv "${evn_removed}" "${EVN_DIR}/GRCh38_ELS.silencer_overlaps.bed"
for chromosome in "${chromosomes[@]}"; do
  mkdir -p "${EVN_CHROMS_DIR}/chrom${chromosome}" "${SVN_CHROMS_DIR}/chrom${chromosome}"
  mv "${work_dir}/evn_chroms/chrom${chromosome}/chrom${chromosome}_GRCh38_ELS.bed" "${EVN_CHROMS_DIR}/chrom${chromosome}/"
  mv "${work_dir}/svn_chroms/chrom${chromosome}/chrom${chromosome}_GRCh38_Silencer_ELS.bed" "${SVN_CHROMS_DIR}/chrom${chromosome}/"
done

printf 'Silencer pELS/dELS source rows:  %d\n' "${silencer_source_count}"
printf 'Unique silencer cCRE loci:       %d\n' "${silencer_count}"
printf 'Removed duplicate silencer rows: %d\n' "${silencer_duplicate_count}"
printf 'Complete SCREEN ELS source rows: %d\n' "${enhancer_source_count}"
printf 'Unique SCREEN ELS loci:          %d\n' "${source_count}"
printf 'Removed duplicate ELS rows:      %d\n' "${enhancer_duplicate_count}"
printf 'Excluded overlap ELS records:    %d\n' "${removed_count}"
printf 'Retained clean enhancer records: %d\n' "${clean_count}"
printf 'Residual EVN/SVN overlaps:       %d\n' "${remaining_overlap}"
