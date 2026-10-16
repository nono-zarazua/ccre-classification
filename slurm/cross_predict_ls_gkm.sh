#!/usr/bin/env bash
## Cross-task evaluation of cumulative chromosome LS-GKM models.
##
## Submit once per background type:
##   sbatch --export=ALL,TYPE=real slurm/cross_predict_ls_gkm.sh
##
## Array tasks  1-21: EVN models evaluated on SVN data.
## Array tasks 22-42: SVN models evaluated on EVN data.
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

echo "Job ${SLURM_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}, node $(hostname), start $(date)"
start_time="$(date -u +%s)"

module load lang/Miniforge3/24.1.2-0
set +u
source activate ls-gkm
set -u

background_type="${TYPE:-real}"
if [[ "${background_type}" != "real" ]]; then
  echo "ERROR: cross-task prediction is currently restricted to TYPE=real; received: ${background_type}" >&2
  exit 1
fi

chromosomes=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
models_per_task="${#chromosomes[@]}"
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 1 && task_id <= 2 * models_per_task )) || {
  echo "ERROR: SLURM_ARRAY_TASK_ID must be 1-$((2 * models_per_task)); received: ${task_id}" >&2
  exit 1
}

if (( task_id <= models_per_task )); then
  model_task="evn"
  evaluation_task="svn"
  checkpoint_index="${task_id}"
else
  model_task="svn"
  evaluation_task="evn"
  checkpoint_index="$((task_id - models_per_task))"
fi

endpoint="${chromosomes[checkpoint_index-1]}"
checkpoint="through_chr${endpoint}"
model="models/${model_task}/GRCh38/ls-gkm/${background_type}/chromosome_learning_curves/${checkpoint}/${checkpoint}.model.txt"
evaluation_root="data/processed/${evaluation_task}/GRCh38/${background_type}"
output_dir="predictions/cross_task/GRCh38/ls-gkm/${background_type}/${model_task}_model_on_${evaluation_task}/chromosome_learning_curves/${checkpoint}"
gkmpredict="/work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict"

[[ -x "${gkmpredict}" ]] || { echo "ERROR: executable not found: ${gkmpredict}" >&2; exit 1; }
[[ -s "${model}" ]] || { echo "ERROR: model missing or empty: ${model}" >&2; exit 1; }

for split in train_100 validation test; do
  for class_name in pos neg; do
    input="${evaluation_root}/${split}_${class_name}.fa"
    [[ -s "${input}" ]] || { echo "ERROR: evaluation FASTA missing or empty: ${input}" >&2; exit 1; }
  done
done

mkdir -p "${output_dir}"
echo "Background: ${background_type}"
echo "Direction: ${model_task} model -> ${evaluation_task} data"
echo "Checkpoint: ${checkpoint}"

for split in train_100 validation test; do
  for class_name in pos neg; do
    input="${evaluation_root}/${split}_${class_name}.fa"
    output="${output_dir}/${split}_${class_name}_scores.txt"

    if [[ -s "${output}" ]]; then
      echo "SKIP: ${output}"
      continue
    fi
    [[ ! -e "${output}" ]] || { echo "ERROR: empty output already exists: ${output}" >&2; exit 1; }

    "${gkmpredict}" -T "${SLURM_CPUS_PER_TASK}" "${input}" "${model}" "${output}"
    [[ -s "${output}" ]] || { echo "ERROR: prediction output missing or empty: ${output}" >&2; exit 1; }
  done
done

elapsed="$(( $(date -u +%s) - start_time ))"
echo "PASS ${model_task}->${evaluation_task}/${background_type}/${checkpoint}; elapsed ${elapsed}s"
