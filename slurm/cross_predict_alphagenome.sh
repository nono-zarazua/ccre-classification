#!/usr/bin/env bash
#SBATCH --job-name=cross_ag
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=slurmlogs/cross_alphagenome.%A_%a.out
#SBATCH --error=slurmlogs/cross_alphagenome.%A_%a.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodelist=n-hpc-gz6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=1-12%4

set -euo pipefail
module load lang/Miniforge3/24.1.2-0
set +u
source activate alphagenome
set -u

fractions=(10 20 40 60 80 100)
task_id=${SLURM_ARRAY_TASK_ID}
if (( task_id < 1 || task_id > 12 )); then
    echo "ERROR: task must be 1-12" >&2
    exit 1
fi
if (( task_id <= 6 )); then
    model_task=evn
    evaluation_task=svn
    fraction=${fractions[task_id-1]}
else
    model_task=svn
    evaluation_task=evn
    fraction=${fractions[task_id-7]}
fi

prefix="models/${model_task}/GRCh38/alphagenome/learning_curves/${fraction}/learning_curve_${fraction}"
checkpoint="${prefix}.best_head.pt"
output_dir="predictions/cross_task/GRCh38/alphagenome/${model_task}_model_on_${evaluation_task}/${fraction}"
mkdir -p "${output_dir}"

for split in validation test; do
    output="${output_dir}/${split}.tsv"
    if [[ -s "${output}" ]]; then
        echo "SKIP: ${output}"
        continue
    fi
    [[ ! -e "${output}" ]] || { echo "ERROR: empty/partial output: ${output}" >&2; exit 1; }
    "${CONDA_PREFIX}/bin/python" src/predict_alphagenome.py \
        "${checkpoint}" \
        "data/processed/${evaluation_task}/GRCh38/training/${split}_pos.fa" \
        "data/processed/${evaluation_task}/GRCh38/training/${split}_neg.fa" \
        "${output}" \
        --weights weights/fold_0_weights.safetensors \
        --batch-size 32 \
        --num-workers "${SLURM_CPUS_PER_TASK}"
done
