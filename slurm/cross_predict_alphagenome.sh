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
#SBATCH --array=1-42%4
set -euo pipefail
module load lang/Miniforge3/24.1.2-0
set +u; source activate alphagenome; set -u
py="$CONDA_PREFIX/bin/python"; "$py" -c 'import torch; assert torch.cuda.is_available()'
chroms=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
id="$SLURM_ARRAY_TASK_ID"; ((id>=1 && id<=42)) || exit 1
if ((id<=21)); then mt=evn; et=svn; ci=$id; else mt=svn; et=evn; ci=$((id-21)); fi
name="through_chr${chroms[ci-1]}"
checkpoint="models/$mt/GRCh38/alphagenome/chromosome_learning_curves/$name/$name.best_head.pt"
root="data/processed/$et/GRCh38/training"
out="predictions/cross_task/GRCh38/alphagenome/${mt}_model_on_$et/chromosome_learning_curves/$name"
[[ -s "$checkpoint" ]] || { echo "ERROR: missing $checkpoint" >&2; exit 1; }; mkdir -p "$out"
for split in train_100 validation test; do
  output="$out/${split}_metrics.json"; [[ ! -s "$output" ]] || continue; [[ ! -e "$output" ]] || exit 1
  "$py" src/predict_alphagenome.py "$checkpoint" "$root/${split}_pos.fa" "$root/${split}_neg.fa" --metrics-json "$output" --weights weights/fold_0_weights.safetensors --batch-size 32 --num-workers "${AG_NUM_WORKERS:-0}"
done
