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
# antinoise is single-threaded. Completed large-chromosome jobs used <100 MB
# MaxRSS, so 1G retains a greater than 10x observed-memory safety margin.
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-46
################################################################################
### Safer Bash behavior
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

# Tasks 1-23: EVN chromosomes 1-22, X.
# Tasks 24-46: SVN chromosomes 1-22, X.
if (( task_id >= 1 && task_id <= 23 )); then
    dataset="evn"
    chromosome_index="${task_id}"
elif (( task_id >= 24 && task_id <= 46 )); then
    dataset="svn"
    chromosome_index="$((task_id - 23))"
else
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 46" >&2
    exit 1
fi

if (( chromosome_index <= 22 )); then
    chromosome="${chromosome_index}"
elif (( chromosome_index == 23 )); then
    chromosome="X"
else
    echo "ERROR: invalid task mapping for SLURM_ARRAY_TASK_ID=${task_id}" >&2
    exit 1
fi

reference_dir="data/reference/hg38-PLAIN-canonical-masked/"
if [[ "${dataset}" == "evn" ]]; then
    dataset_name="ELS"
else
    dataset_name="Silencer_ELS"
fi

chromosome_dir="data/processed/${dataset}/GRCh38/chroms/chrom${chromosome}"
foreground_fasta="${chromosome_dir}/chrom${chromosome}_GRCh38_${dataset_name}.fa"
output_prefix="${chromosome_dir}/real_neg1x_chrom${chromosome}_GRCh38_${dataset_name}"
executable="../antinoise/bin/linux/background_genome_mono.exe"

[[ -x "${executable}" ]] || { echo "ERROR: executable not found: ${executable}" >&2; exit 1; }
[[ -d "${reference_dir}" ]] || { echo "ERROR: reference directory not found: ${reference_dir}" >&2; exit 1; }
[[ -s "${reference_dir}chr${chromosome}.plain" ]] || { echo "ERROR: masked chromosome reference missing: ${reference_dir}chr${chromosome}.plain" >&2; exit 1; }
[[ -s "${foreground_fasta}" ]] || { echo "ERROR: foreground FASTA missing or empty: ${foreground_fasta}" >&2; exit 1; }

echo "Array task ID: ${task_id}"
echo "Dataset: ${dataset}"
echo "Chromosome: ${chromosome}"
echo "Foreground FASTA: ${foreground_fasta}"
echo "Output prefix: ${output_prefix}"
echo "Backgrounds per foreground sequence: 1"

"${executable}" \
    "${reference_dir}" \
    "${foreground_fasta}" \
    "${output_prefix}" \
    1 \
    0.01 \
    10000 \
    hg38 \
    0.99 \
    "${output_prefix}.outm" \
    "${output_prefix}.outd" \
    "${output_prefix}.outm_one" \
    "${output_prefix}.outd_one" \
    "${output_prefix}.outlog"

[[ -s "${output_prefix}.fa" ]] || {
    echo "ERROR: antinoise did not create a non-empty FASTA: ${output_prefix}.fa" >&2
    exit 1
}
[[ -s "${output_prefix}.bed" ]] || {
    echo "ERROR: antinoise did not create a non-empty BED: ${output_prefix}.bed" >&2
    exit 1
}

foreground_count=$(grep -c '^>' "${foreground_fasta}")
background_count=$(grep -c '^>' "${output_prefix}.fa")
echo "Foreground records: ${foreground_count}"
echo "Generated real-background records: ${background_count}"
if (( background_count > foreground_count )); then
    echo "ERROR: generated more than one background per foreground sequence" >&2
    exit 1
fi

end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
