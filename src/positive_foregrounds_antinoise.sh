#!/bin/bash

# Extract the selected EVN+SVN foreground sequences from the unmasked hg38
# reference. The masked reference is reserved for background sampling.
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${repo_root}"

extractor="../antinoise/bin/linux/longext_many.exe"
reference_dir="data/reference/antinoise/hg38-PLAIN-canonical/"
chromosomes_dir="data/processed/evn/GRCh38/chroms"

[[ -x "${extractor}" ]] || { echo "ERROR: missing executable: ${extractor}" >&2; exit 1; }

for chr in {1..22} X; do
    directory="${chromosomes_dir}/chrom${chr}"
    reference="${reference_dir}chr${chr}.plain"
    input_bed="${directory}/chrom${chr}_GRCh38_ELS.bed"
    output_fasta="${directory}/chrom${chr}_GRCh38_ELS.fa"
    error_log="${directory}/longext_many.log"

    [[ -s "${reference}" ]] || { echo "ERROR: missing reference: ${reference}" >&2; exit 1; }
    [[ -s "${input_bed}" ]] || { echo "ERROR: missing BED: ${input_bed}" >&2; exit 1; }

    temporary_fasta=$(mktemp "${directory}/.chrom${chr}_GRC38_ELS.fa.XXXXXX")
    temporary_log=$(mktemp "${directory}/.longext_many.log.XXXXXX")

    # longext_many returns 1 after successful completion, so validate its files
    # rather than interpreting status 1 as an automatic failure.
    set +e
    "${extractor}" \
        "${reference_dir}" \
        "${input_bed}" \
        "${temporary_fasta}" \
        0 \
        0 \
        1000 \
        hg38 \
        "${temporary_log}"
    status=$?
    set -e

    if [[ "${status}" -ne 1 ]]; then
        rm -f "${temporary_fasta}" "${temporary_log}"
        echo "ERROR: longext_many returned unexpected status ${status} for chr${chr}" >&2
        exit 1
    fi
    if [[ -s "${temporary_log}" ]]; then
        cat "${temporary_log}" >&2
        rm -f "${temporary_fasta}" "${temporary_log}"
        echo "ERROR: longext_many reported extraction problems for chr${chr}" >&2
        exit 1
    fi
    if [[ ! -s "${temporary_fasta}" ]]; then
        rm -f "${temporary_fasta}" "${temporary_log}"
        echo "ERROR: longext_many created an empty FASTA for chr${chr}" >&2
        exit 1
    fi

    bed_records=$(wc -l < "${input_bed}")
    fasta_records=$(grep -c '^>' "${temporary_fasta}" || true)
    if [[ "${bed_records}" -ne "${fasta_records}" ]]; then
        rm -f "${temporary_fasta}" "${temporary_log}"
        echo "ERROR: chr${chr} BED/FASTA mismatch: ${bed_records} vs ${fasta_records}" >&2
        exit 1
    fi

    mv "${temporary_fasta}" "${output_fasta}"
    mv "${temporary_log}" "${error_log}"
    printf 'chr%s: %d foreground sequences\n' "${chr}" "${fasta_records}"
done
