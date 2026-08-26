#!/bin/bash
## Train frozen-encoder AlphaGenome enhancer classifiers.
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=train_alphagenome
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_alphagenome.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_alphagenome.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-24%1

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
start_time="$(date -u +%s)"

module load lang/Miniforge3/24.1.2-0
source activate alphagenome

# Map tasks to four outer folds x six independently trained nested fractions.
# Tasks 1-6 are outer1, 7-12 outer2, 13-18 outer3, and 19-24 outer4.
percentages=(10 20 40 60 80 100)
task_index="$((SLURM_ARRAY_TASK_ID - 1))"
outer="$((task_index / ${#percentages[@]} + 1))"
percentage_index="$((task_index % ${#percentages[@]}))"
percentage="${percentages[${percentage_index}]}"

if (( outer < 1 || outer > 4 )); then
    echo "ERROR: derived outer number is invalid: ${outer}" >&2
    exit 1
fi

echo "Conda environment: ${CONDA_PREFIX}"
echo "Python: $(command -v python)"
echo "Outer fold: ${outer}"
echo "Training percentage: ${percentage}%"

python alphagenome_trial/train_enhancer_classifier.py \
    --outer "${outer}" \
    --percentage "${percentage}"

end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"
echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
