#!/usr/bin/env bash
#SBATCH --job-name=cross_lsgkm
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=slurmlogs/cross_ls_gkm.%A_%a.out
#SBATCH --error=slurmlogs/cross_ls_gkm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=10:00:00
#SBATCH --array=1-42%10
set -euo pipefail
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm
chroms=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
id="$SLURM_ARRAY_TASK_ID"; ((id>=1 && id<=42)) || exit 1
if ((id<=21)); then mt=evn; et=svn; ci=$id; else mt=svn; et=evn; ci=$((id-21)); fi
name="through_chr${chroms[ci-1]}"
model="models/$mt/GRCh38/ls-gkm/chromosome_learning_curves/$name/$name.model.txt"
root="data/processed/$et/GRCh38/training"
out="predictions/cross_task/GRCh38/ls-gkm/${mt}_model_on_$et/chromosome_learning_curves/$name"
[[ -s "$model" ]] || { echo "ERROR: missing $model" >&2; exit 1; }; mkdir -p "$out"
for split in train_100 validation test; do for cls in pos neg; do
  input="$root/${split}_$cls.fa"; output="$out/${split}_${cls}_scores.txt"
  [[ -s "$input" ]] || exit 1; [[ ! -s "$output" ]] || continue; [[ ! -e "$output" ]] || exit 1
  /work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict -T "$SLURM_CPUS_PER_TASK" "$input" "$model" "$output"
  [[ -s "$output" ]] || exit 1
done; done
