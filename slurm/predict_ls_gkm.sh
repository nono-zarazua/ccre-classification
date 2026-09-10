#!/usr/bin/env bash
## Evaluate every cumulative chromosome LS-GKM model on train_100, validation, and test.
#SBATCH --job-name=predict_ls_gkm
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_ls_gkm.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/predict_ls_gkm.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=10:00:00
#SBATCH --array=1-21%10

set -euo pipefail
start_time="$(date -u +%s)"
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm

dataset="${DATASET:-evn}"
[[ "${dataset}" == evn || "${dataset}" == svn ]] || { echo "ERROR: DATASET must be evn or svn" >&2; exit 1; }
negative_kind="${NEGATIVE_KIND:-synthetic}"
[[ "${negative_kind}" == synthetic || "${negative_kind}" == real ]] || { echo "ERROR: NEGATIVE_KIND must be synthetic or real" >&2; exit 1; }

chromosomes=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 1 && task_id <= ${#chromosomes[@]} )) || { echo "ERROR: task must be 1-${#chromosomes[@]}" >&2; exit 1; }
endpoint="${chromosomes[task_id-1]}"; name="through_chr${endpoint}"
root="data/processed/${dataset}/GRCh38/training/${negative_kind}"
model="models/${dataset}/GRCh38/${negative_kind}/ls-gkm/chromosome_learning_curves/${name}/${name}.model.txt"
output_dir="predictions/${dataset}/GRCh38/${negative_kind}/ls-gkm/chromosome_learning_curves/${name}"

for file in "${model}" "${root}/train_100_pos.fa" "${root}/train_100_neg.fa" "${root}/validation_pos.fa" "${root}/validation_neg.fa" "${root}/test_pos.fa" "${root}/test_neg.fa"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty ${file}" >&2; exit 1; }
done
mkdir -p "${output_dir}"
for split in train_100 validation test; do
  for class_name in pos neg; do
    input="${root}/${split}_${class_name}.fa"
    output="${output_dir}/${split}_${class_name}_scores.txt"
    [[ ! -s "${output}" ]] || { echo "SKIP: ${output}"; continue; }
    [[ ! -e "${output}" ]] || { echo "ERROR: empty output exists: ${output}" >&2; exit 1; }
    /work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict -T "${SLURM_CPUS_PER_TASK}" "${input}" "${model}" "${output}"
    [[ -s "${output}" ]] || { echo "ERROR: prediction missing: ${output}" >&2; exit 1; }
  done
done
elapsed="$(($(date -u +%s)-start_time))"
echo "PASS ${dataset}/${negative_kind}/${name}; elapsed ${elapsed}s"
