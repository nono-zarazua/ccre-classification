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

bash src/make_outer_files.sh "${TEST_FOLD}"

TRAIN_DIR="data/processed/ls-gkm/training/outer${TEST_FOLD}"
MODEL_DIR="models/ls-gkm"

TRAIN_POS="${TRAIN_DIR}/train_pos.fa"
TRAIN_NEG="${TRAIN_DIR}/train_neg.fa"

mkdir -p "${MODEL_DIR}"

echo "Outer test fold: ${TEST_FOLD}"
echo "Training positive sequences: $(grep -c '^>' "${TRAIN_POS}")"
echo "Training negative sequences: $(grep -c '^>' "${TRAIN_NEG}")"

/work/zarazuanav/workspace/repos/lsgkm/bin/gkmtrain \
	-T "${SLURM_CPUS_PER_TASK}" \
	-m 16000 \
	"${TRAIN_POS}" \
	"${TRAIN_NEG}" \
	"${MODEL_DIR}/outer${TEST_FOLD}"
################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
