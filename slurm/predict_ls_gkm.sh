#!/bin/bash
## Written by Nono Zarazua
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=predict_ls_gkm
#SBATCH --comment=
### changedir / output
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_ls_gkm.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_ls_gkm.%A_%a.err

#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#
### Node / CPU / Memory Settings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
# --qos=long
#SBATCH --array=1-5
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
set -euo pipefail

echo "Conda environment: ${CONDA_PREFIX}"

TEST_FOLD="${SLURM_ARRAY_TASK_ID}"

FOLD_DIR="data/processed/folds"
MODEL="models/ls-gkm/outer${TEST_FOLD}.model.txt"
OUTPUT_DIR="predictions/ls-gkm/outer${TEST_FOLD}"

TEST_POS="${FOLD_DIR}/fold${TEST_FOLD}/fold${TEST_FOLD}.fa"
TEST_NEG="${FOLD_DIR}/fold${TEST_FOLD}/neg1x_fold${TEST_FOLD}.fa"

OUT_POS="${OUTPUT_DIR}/pos_scores.txt"
OUT_NEG="${OUTPUT_DIR}/neg_scores.txt"

mkdir -p "predictions/ls-gkm/outer${TEST_FOLD}"

echo "Positive predictions outer${TEST_FOLD}"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
	-T "${SLURM_CPUS_PER_TASK}" \
	"${TEST_POS}" \
	"${MODEL}" \
	"${OUT_POS}"


echo "Negative predictions outer${TEST_FOLD}"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
	-T "${SLURM_CPUS_PER_TASK}" \
	"${TEST_NEG}" \
	"${MODEL}" \
	"${OUT_NEG}"### Running code
################################################################################
set -euo pipefail

echo "Conda environment: ${CONDA_PREFIX}"

TEST_FOLD="${SLURM_ARRAY_TASK_ID}"

FOLD_DIR="data/processed/folds"
MODEL="models/ls-gkm/outer${TEST_FOLD}.model.txt"
OUTPUT_DIR="predictions/ls-gkm/outer${TEST_FOLD}"

TEST_POS="${FOLD_DIR}/fold${TEST_FOLD}/fold${TEST_FOLD}.fa"
TEST_NEG="${FOLD_DIR}/fold${TEST_FOLD}/neg1x_fold${TEST_FOLD}.fa"

OUT_POS="${OUTPUT_DIR}/pos_scores.txt"
OUT_NEG="${OUTPUT_DIR}/neg_scores.txt"

mkdir -p "${OUTPUT_DIR}"

for FILE in "${TEST_POS}" "${TEST_NEG}" "${MODEL}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: missing or empty file: ${FILE}" >&2
        exit 1
    fi
done

echo "Positive predictions outer${TEST_FOLD}"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
        -T "${SLURM_CPUS_PER_TASK}" \
        "${TEST_POS}" \
        "${MODEL}" \
        "${OUT_POS}"


echo "Negative predictions outer${TEST_FOLD}"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
        -T "${SLURM_CPUS_PER_TASK}" \
        "${TEST_NEG}" \
        "${MODEL}" \
        "${OUT_NEG}"

################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
