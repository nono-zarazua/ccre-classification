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
### Node / CPU / Memory Settings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
# --qos=long
#SBATCH --array=1-24%6
################################################################################
### Safer Bash behavior
################################################################################
set -euo pipefail
################################################################################
### Job Info
################################################################################
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
###############################################################################
### Load modules / source environment
################################################################################
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm
################################################################################
### Running code
################################################################################
echo "Conda environment: ${CONDA_PREFIX}"

# Map 24 tasks to four outers x six nested training percentages.
percentages=(100 80 60 40 20 10)
task_index="$((SLURM_ARRAY_TASK_ID - 1))"
outer="$((task_index / ${#percentages[@]} + 1))"
percentage_index="$((task_index % ${#percentages[@]}))"
percentage="${percentages[${percentage_index}]}"

if (( outer < 1 || outer > 4 )); then
    echo "ERROR: derived outer number is invalid: ${outer}" >&2
    exit 1
fi

TRAIN_DIR="data/processed/evn/GRCh38/learning_curves/outer${outer}"
TRAIN_POS="${TRAIN_DIR}/train_${percentage}_pos.fa"
TRAIN_NEG="${TRAIN_DIR}/train_${percentage}_neg.fa"

# Percentage-specific directories prevent the six models for an outer from
# overwriting one another while retaining the requested outerN.model.txt name.
MODEL_DIR="models/evn/GRCh38/ls-gkm/outer${outer}/${percentage}"
MODEL_PREFIX="${MODEL_DIR}/outer${outer}"

for input_fasta in "${TRAIN_POS}" "${TRAIN_NEG}"; do
    if [[ ! -s "${input_fasta}" ]]; then
        echo "ERROR: required training FASTA is missing or empty: ${input_fasta}" >&2
        exit 1
    fi
done

mkdir -p "${MODEL_DIR}"

echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Outer: ${outer}"
echo "Training fraction: ${percentage}%"
echo "Training positive FASTA: ${TRAIN_POS}"
echo "Training negative FASTA: ${TRAIN_NEG}"
echo "Training positive sequences: $(grep -c '^>' "${TRAIN_POS}")"
echo "Training negative sequences: $(grep -c '^>' "${TRAIN_NEG}")"
echo "Model output: ${MODEL_PREFIX}.model.txt"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmtrain \
    -T "${SLURM_CPUS_PER_TASK}" \
    -m 64000 \
    "${TRAIN_POS}" \
    "${TRAIN_NEG}" \
    "${MODEL_PREFIX}"

################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
