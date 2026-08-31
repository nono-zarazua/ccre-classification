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
#SBATCH --array=7-10
################################################################################
### Safer Bash behavior
################################################################################
set -euo pipefail
################################################################################
### Job Info
################################################################################
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
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

task_id="${SLURM_ARRAY_TASK_ID}"
if (( task_id < 1 || task_id > 10 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 10" >&2
    exit 1
fi

percentages=(10 20 40 60 80 100)
if (( task_id <= 6 )); then
    percentage="${percentages[task_id - 1]}"
    dataset_name="learning_curve_${percentage}"
    validation_pos="data/processed/svn/GRCh38/training/validation_pos.fa"
    validation_neg="data/processed/svn/GRCh38/training/validation_neg.fa"
    test_pos="data/processed/svn/GRCh38/training/test_pos.fa"
    test_neg="data/processed/svn/GRCh38/training/test_neg.fa"
    model="models/svn/GRCh38/ls-gkm/learning_curves/${percentage}/learning_curve_${percentage}.model.txt"
    output_dir="predictions/svn/GRCh38/ls-gkm/learning_curves/${percentage}"
    echo "Dataset type: fixed learning curve"
    echo "Predicting fraction: ${percentage}%"
else
    outer="$((task_id - 6))"
    dataset_name="outer${outer}"
    validation_pos="data/processed/svn/GRCh38/training/outer${outer}/validation_pos.fa"
    validation_neg="data/processed/svn/GRCh38/training/outer${outer}/validation_neg.fa"
    test_pos="data/processed/svn/GRCh38/folds/fold5/fold5.fa"
    test_neg="data/processed/svn/GRCh38/folds/fold5/neg1x_fold5.fa"
    model="models/svn/GRCh38/ls-gkm/outer${outer}/outer${outer}.model.txt"
    output_dir="predictions/svn/GRCh38/ls-gkm/outer${outer}"
    echo "Dataset type: chromosome-based outer"
    echo "Outer set: ${outer}"
fi

for input_file in \
    "${model}" \
    "${validation_pos}" \
    "${validation_neg}" \
    "${test_pos}" \
    "${test_neg}"; do
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: missing or empty input: ${input_file}" >&2
        exit 1
    fi
done

mkdir -p "${output_dir}"

for split in validation test; do
    for class_name in pos neg; do
        input_variable="${split}_${class_name}"
        input_fasta="${!input_variable}"
        output_file="${output_dir}/${split}_${class_name}_scores.txt"

        if [[ -s "${output_file}" ]]; then
            echo "SKIP: non-empty prediction already exists: ${output_file}"
            continue
        fi
        if [[ -e "${output_file}" ]]; then
            echo "ERROR: existing prediction file is empty: ${output_file}" >&2
            exit 1
        fi

        echo "Predicting ${dataset_name}: ${split} ${class_name}"
        echo "Input: ${input_fasta}"
        echo "Model: ${model}"
        echo "Output: ${output_file}"

        /work/zarazuanav/workspace/repos/lsgkm/bin/gkmpredict \
            -T "${SLURM_CPUS_PER_TASK}" \
            "${input_fasta}" \
            "${model}" \
            "${output_file}"

        if [[ ! -s "${output_file}" ]]; then
            echo "ERROR: gkmpredict did not create ${output_file}" >&2
            exit 1
        fi
    done
done

################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
