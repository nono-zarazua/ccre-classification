#!/bin/bash

set -euo pipefail

consolidation="data/processed/consolidation/GRCh38/chromosomes"
executable="../antinoise/bin/linux/mix0.exe"

[[ -x "${executable}" ]] || {
    echo "ERROR: executable not found: ${executable}" >&2
    exit 1
}

for chromosome in {1..22} X; do
    chromosome_dir="${consolidation}/chrom${chromosome}"
    input="${chromosome_dir}/consolidation.fa"
    output="${chromosome_dir}/synthetic_antinoise_candidates_2x.fa"
    temporary="${output}.tmp.$$"

    [[ -s "${input}" ]] || {
        echo "ERROR: missing or empty foreground FASTA: ${input}" >&2
        exit 1
    }
    rm -f "${temporary}"

    echo "Processing chromosome ${chromosome}"

    # mix0.exe returns 1 even after successful completion. Keep strict mode for
    # the rest of the script, capture this command's status explicitly, and
    # determine success from both the known return code and validated output.
    if "${executable}" "${input}" "${temporary}" 2; then
        status=0
    else
        status=$?
    fi
    if [[ "${status}" -ne 1 ]]; then
        echo "ERROR: mix0.exe returned unexpected status ${status} for chr${chromosome}" >&2
        rm -f "${temporary}"
        exit 1
    fi
    [[ -s "${temporary}" ]] || {
        echo "ERROR: mix0.exe did not create a non-empty FASTA for chr${chromosome}" >&2
        rm -f "${temporary}"
        exit 1
    }

    foreground_count=$(grep -c '^>' "${input}")
    synthetic_count=$(grep -c '^>' "${temporary}")
    expected_count=$((2 * foreground_count))
    if [[ "${synthetic_count}" -ne "${expected_count}" ]]; then
        echo "ERROR: chr${chromosome}: expected ${expected_count} synthetic records, observed ${synthetic_count}" >&2
        rm -f "${temporary}"
        exit 1
    fi

    mv "${temporary}" "${output}"
    echo "Generated ${synthetic_count} candidates from ${foreground_count} foregrounds"
done
