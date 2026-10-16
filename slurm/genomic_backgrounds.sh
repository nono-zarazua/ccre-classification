#!/bin/bash
## Written by Nono Zarazua
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=genomic_backgrounds
#SBATCH --comment=
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/genomic_backgrounds.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/genomic_backgrounds.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
# antinoise is single-threaded. Previous jobs used less than 100 MB MaxRSS;
# 2 GB retains a large safety margin.
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-23
################################################################################

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
start_time="$(date -u +%s)"

echo "################################################################################"
echo "Used SLURM script:"
echo
cat "$0"
echo
echo "################################################################################"
echo "Output:"
echo "################################################################################"

task_id="${SLURM_ARRAY_TASK_ID}"
if (( task_id < 1 || task_id > 23 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 23" >&2
    exit 1
fi
if (( task_id <= 22 )); then
    chromosome="${task_id}"
else
    chromosome="X"
fi

reference_dir="data/reference/antinoise/hg38-PLAIN-canonical-masked/chr${chromosome}/"
reference_file="${reference_dir}chr${chromosome}.plain"
chromosome_dir="data/processed/evn/GRCh38/chroms/chrom${chromosome}"
foreground_fasta="${chromosome_dir}/chrom${chromosome}_GRCh38_ELS.fa"
output_prefix="${chromosome_dir}/antinoise_candidates_1x"
executable="bin/background_genome_mono_split.exe"

[[ -x "${executable}" ]] || { echo "ERROR: executable not found: ${executable}" >&2; exit 1; }
[[ -d "${reference_dir}" ]] || { echo "ERROR: reference directory not found: ${reference_dir}" >&2; exit 1; }
[[ -s "${reference_file}" ]] || { echo "ERROR: masked chromosome reference missing: ${reference_file}" >&2; exit 1; }
[[ -s "${foreground_fasta}" ]] || { echo "ERROR: foreground FASTA missing or empty: ${foreground_fasta}" >&2; exit 1; }

mapfile -t plain_files < <(find "${reference_dir}" -maxdepth 1 -type f -name '*.plain' -print)
if [[ "${#plain_files[@]}" -ne 1 || "${plain_files[0]}" != "${reference_file}" ]]; then
    echo "ERROR: ${reference_dir} must contain only chr${chromosome}.plain" >&2
    exit 1
fi

foreground_count=$(grep -c '^>' "${foreground_fasta}")
[[ "${foreground_count}" -gt 0 ]] || { echo "ERROR: no foreground records for chr${chromosome}" >&2; exit 1; }

echo "Array task ID: ${task_id}"
echo "Foreground collection: consolidated EVN+SVN"
echo "Chromosome: ${chromosome}"
echo "Reference directory: ${reference_dir}"
echo "Foreground FASTA: ${foreground_fasta}"
echo "Foreground records: ${foreground_count}"
echo "Output prefix: ${output_prefix}"
echo "Requested candidates per foreground: 1"

"${executable}" \
    "${reference_dir}" \
    "${foreground_fasta}" \
    "${output_prefix}" \
    1 \
    0.01 \
    10000 \
    hg38 \
    0.99 \
    0 \
    "${output_prefix}.outm" \
    "${output_prefix}.outd" \
    "${output_prefix}.outm_one" \
    "${output_prefix}.outd_one" \
    "${output_prefix}.outlog" \
    "${chromosome}"

for suffix in fa bed outm outd outm_one outd_one outlog; do
    output="${output_prefix}.${suffix}"
    [[ -s "${output}" ]] || {
        echo "ERROR: antinoise did not create a non-empty output: ${output}" >&2
        exit 1
    }
done

fasta_count=$(grep -c '^>' "${output_prefix}.fa")
bed_count=$(wc -l < "${output_prefix}.bed")
maximum_count=$((2 * foreground_count))

[[ "${fasta_count}" -eq "${bed_count}" ]] || {
    echo "ERROR: candidate BED/FASTA mismatch: ${bed_count} vs ${fasta_count}" >&2
    exit 1
}
(( fasta_count <= maximum_count )) || {
    echo "ERROR: generated ${fasta_count} candidates; maximum is ${maximum_count}" >&2
    exit 1
}
awk -v expected="chr${chromosome}" 'BEGIN{FS="\t"} $1 != expected {bad++} END{exit bad ? 1 : 0}' \
    "${output_prefix}.bed" || {
    echo "ERROR: candidate BED contains a chromosome other than chr${chromosome}" >&2
    exit 1
}

echo "Generated candidates: ${fasta_count}"
echo "Maximum requested candidates: ${maximum_count}"
echo "Candidate yield: $(awk -v n="${fasta_count}" -v d="${maximum_count}" 'BEGIN { printf "%.4f", n/d }')"
echo "Final antinoise status:"
cat "${output_prefix}.outlog"

end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"
echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
