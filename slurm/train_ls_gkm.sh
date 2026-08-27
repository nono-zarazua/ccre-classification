#!/bin/bash
## Written by Nono Zarazua
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=train_ls_gkm
#SBATCH --comment=
### changedir / output
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_ls_gkm.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_ls_gkm.%A_%a.err

#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#
### Array / Node / CPU / Memory Settings
#SBATCH --array=1-10%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=05-00:00:00
# --qos=long
################################################################################
### Safer Bash behavior
################################################################################
set -euo pipefail

################################################################################
### Job Info
################################################################################
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
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

################################################################################
### Load modules / source environment
################################################################################
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm

################################################################################
### Map array task to one training dataset
################################################################################
echo "Conda environment: ${CONDA_PREFIX}"

task_id="${SLURM_ARRAY_TASK_ID}"
if (( task_id < 1 || task_id > 10 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 10" >&2
    exit 1
fi

percentages=(10 20 40 60 80 100)
if (( task_id <= 6 )); then
    percentage="${percentages[task_id - 1]}"
    dataset_name="learning_curve_${percentage}"
    train_dir="data/processed/svn/GRCh38/learning_curves"
    train_pos="${train_dir}/train_${percentage}_pos.fa"
    train_neg="${train_dir}/train_${percentage}_neg.fa"
    model_dir="models/svn/GRCh38/ls-gkm/learning_curves/${percentage}"
    model_prefix="${model_dir}/learning_curve_${percentage}"
    echo "Dataset type: fixed learning curve"
    echo "Training fraction: ${percentage}%"
else
    outer="$((task_id - 6))"
    dataset_name="outer${outer}"
    train_dir="data/processed/svn/GRCh38/training/outer${outer}"
    train_pos="${train_dir}/train_pos.fa"
    train_neg="${train_dir}/train_neg.fa"
    model_dir="models/svn/GRCh38/ls-gkm/outer${outer}"
    model_prefix="${model_dir}/outer${outer}"
    echo "Dataset type: chromosome-based outer"
    echo "Outer set: ${outer}"
fi

model_file="${model_prefix}.model.txt"
for input_fasta in "${train_pos}" "${train_neg}"; do
    if [[ ! -s "${input_fasta}" ]]; then
        echo "ERROR: required training FASTA is missing or empty: ${input_fasta}" >&2
        exit 1
    fi
done

if [[ -s "${model_file}" ]]; then
    echo "SKIP: non-empty model already exists: ${model_file}"
    exit 0
fi
if [[ -e "${model_file}" ]]; then
    echo "ERROR: existing model file is empty: ${model_file}" >&2
    exit 1
fi

mkdir -p "${model_dir}"

echo "Dataset: ${dataset_name}"
echo "Training positive FASTA: ${train_pos}"
echo "Training negative FASTA: ${train_neg}"
echo "Training positive sequences: $(grep -c '^>' "${train_pos}")"
echo "Training negative sequences: $(grep -c '^>' "${train_neg}")"
echo "Model output: ${model_file}"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmtrain \
    -T "${SLURM_CPUS_PER_TASK}" \
    -m 160000 \
    "${train_pos}" \
    "${train_neg}" \
    "${model_prefix}"

if [[ ! -s "${model_file}" ]]; then
    echo "ERROR: gkmtrain exited successfully but did not create ${model_file}" >&2
    exit 1
fi
################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
