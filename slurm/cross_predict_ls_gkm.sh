#!/usr/bin/env bash
#SBATCH --job-name=cross_lsgkm
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=slurmlogs/cross_ls_gkm.%A_%a.out
#SBATCH --error=slurmlogs/cross_ls_gkm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-12%6

set -euo pipefail
module load lang/Miniforge3/24.1.2-0
set +u
source activate ls-gkm
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

model="models/${model_task}/GRCh38/ls-gkm/learning_curves/${fraction}/learning_curve_${fraction}.model.txt"
output_dir="predictions/cross_task/GRCh38/ls-gkm/${model_task}_model_on_${evaluation_task}/${fraction}"
[[ -s "${model}" ]] || { echo "ERROR: missing model: ${model}" >&2; exit 1; }
mkdir -p "${output_dir}"

for split in validation test; do
    for class_name in pos neg; do
        input="data/processed/${evaluation_task}/GRCh38/training/${split}_${class_name}.fa"
        output="${output_dir}/${split}_${class_name}_scores.txt"
        [[ -s "${input}" ]] || { echo "ERROR: missing input: ${input}" >&2; exit 1; }
        if [[ -s "${output}" ]]; then
            echo "SKIP: ${output}"
            continue
        fi
        [[ ! -e "${output}" ]] || { echo "ERROR: empty/partial output: ${output}" >&2; exit 1; }
        /work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
            -T "${SLURM_CPUS_PER_TASK}" "${input}" "${model}" "${output}"
        [[ -s "${output}" ]] || { echo "ERROR: no output: ${output}" >&2; exit 1; }
    done
done
