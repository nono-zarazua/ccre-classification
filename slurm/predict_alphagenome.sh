#!/usr/bin/env bash
## Evaluate every trained AlphaGenome checkpoint on the shared train_100 pool.
## Select DATASET=evn|svn at submission.
#SBATCH --job-name=predict_alphagenome
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_alphagenome.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_alphagenome.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-21%4
set -euo pipefail
echo "Job ${SLURM_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}, node $(hostname), start $(date)"
start_time="$(date -u +%s)"
module load lang/Miniforge3/24.1.2-0
set +u; source activate alphagenome; set -u
python_bin="${CONDA_PREFIX}/bin/python"
[[ -x "${python_bin}" ]] || { echo "ERROR: missing AlphaGenome Python" >&2; exit 1; }
"${python_bin}" -c 'import numpy, torch, sklearn, alphagenome_pytorch; assert torch.cuda.is_available(), "CUDA unavailable"; print(torch.cuda.get_device_name(0))'
dataset="${DATASET:-evn}"
type="${TYPE:-real}"
[[ "${dataset}" == evn || "${dataset}" == svn ]] || { echo "ERROR: DATASET must be evn or svn" >&2; exit 1; }
[[ "${type}" == real || "${type}" == synthetic ]] || { echo "ERROR: TYPE must be real or synthetic" >&2; exit 1; }
chromosomes=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 1 && task_id <= ${#chromosomes[@]} )) || { echo "ERROR: task must be 1-${#chromosomes[@]}" >&2; exit 1; }
endpoint="${chromosomes[task_id-1]}"; name="through_chr${endpoint}"
root="data/processed/${dataset}/GRCh38/${type}"
model_root="models/${dataset}/GRCh38/alphagenome/${type}/chromosome_learning_curves/${name}"
checkpoint="${model_root}/${name}.best_head.pt"
metrics="${model_root}/${name}.train_100_metrics.json"
weights="weights/fold_0_weights.safetensors"
for file in "${checkpoint}" "${weights}" "${root}/train_100_pos.fa" "${root}/train_100_neg.fa"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty ${file}" >&2; exit 1; }
done
[[ ! -s "${metrics}" ]] || { echo "SKIP: ${metrics}"; exit 0; }
[[ ! -e "${metrics}" ]] || { echo "ERROR: empty metrics exists: ${metrics}" >&2; exit 1; }
"${python_bin}" src/predict_alphagenome.py "${checkpoint}" "${root}/train_100_pos.fa" "${root}/train_100_neg.fa" --metrics-json "${metrics}" --weights "${weights}" --batch-size 32 --num-workers 0
[[ -s "${metrics}" ]] || { echo "ERROR: metrics missing: ${metrics}" >&2; exit 1; }
elapsed="$(($(date -u +%s)-start_time))"
echo "PASS ${dataset}/${type}/${name}; elapsed ${elapsed}s"
