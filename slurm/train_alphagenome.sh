#!/usr/bin/env bash
## Train frozen-encoder AlphaGenome binary classifiers for svn learning curves
## and chromosome-composition stability experiments.
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
#SBATCH --array=1-10%4

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
start_time="$(date -u +%s)"

module load lang/Miniforge3/24.1.2-0
# Conda activation hooks inspect optional variables that may be unset.
set +u
source activate alphagenome
set -u

python_bin="${CONDA_PREFIX}/bin/python"
if [[ ! -x "${python_bin}" ]]; then
    echo "ERROR: AlphaGenome Python is not executable: ${python_bin}" >&2
    exit 1
fi
echo "Conda environment: ${CONDA_PREFIX}"
echo "Python: ${python_bin}"
"${python_bin}" -c 'import numpy, torch, sklearn, alphagenome_pytorch; print("AlphaGenome environment import preflight: PASS")'

################################################################################
### Map one array task to one independent experiment
################################################################################
task_id="${SLURM_ARRAY_TASK_ID}"
if (( task_id < 1 || task_id > 10 )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 10" >&2
    exit 1
fi

percentages=(10 20 40 60 80 100)
if (( task_id <= 6 )); then
    percentage="${percentages[task_id - 1]}"
    dataset_name="learning_curve_${percentage}"

    train_pos="data/processed/svn/GRCh38/learning_curves/train_${percentage}_pos.fa"
    train_neg="data/processed/svn/GRCh38/learning_curves/train_${percentage}_neg.fa"
    validation_pos="data/processed/svn/GRCh38/training/validation_pos.fa"
    validation_neg="data/processed/svn/GRCh38/training/validation_neg.fa"
    test_pos="data/processed/svn/GRCh38/training/test_pos.fa"
    test_neg="data/processed/svn/GRCh38/training/test_neg.fa"

    output_dir="models/svn/GRCh38/alphagenome/learning_curves/${percentage}"
    output_prefix="${output_dir}/${dataset_name}"

    echo "Experiment: fixed learning curve"
    echo "Training fraction: ${percentage}%"
else
    outer="$((task_id - 6))"
    dataset_name="outer${outer}"

    train_pos="data/processed/svn/GRCh38/training/outer${outer}/train_pos.fa"
    train_neg="data/processed/svn/GRCh38/training/outer${outer}/train_neg.fa"
    validation_pos="data/processed/svn/GRCh38/training/outer${outer}/validation_pos.fa"
    validation_neg="data/processed/svn/GRCh38/training/outer${outer}/validation_neg.fa"
    test_pos="data/processed/svn/GRCh38/folds/fold5/fold5.fa"
    test_neg="data/processed/svn/GRCh38/folds/fold5/neg1x_fold5.fa"

    output_dir="models/svn/GRCh38/alphagenome/outer${outer}"
    output_prefix="${output_dir}/${dataset_name}"

    echo "Experiment: chromosome-composition stability"
    echo "Outer set: ${outer}"
fi

################################################################################
### Validate inputs and protect completed/partial runs
################################################################################
weights="weights/fold_0_weights.safetensors"
required_inputs=(
    "${train_pos}"
    "${train_neg}"
    "${validation_pos}"
    "${validation_neg}"
    "${test_pos}"
    "${test_neg}"
    "${weights}"
)
for input_file in "${required_inputs[@]}"; do
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: required input is missing or empty: ${input_file}" >&2
        exit 1
    fi
done

artifact_suffixes=(
    .best_head.pt
    .config.json
    .history.tsv
    .validation_metrics.json
    .test_metrics.json
)
existing_artifacts=0
complete_artifacts=0
for suffix in "${artifact_suffixes[@]}"; do
    artifact="${output_prefix}${suffix}"
    [[ -e "${artifact}" ]] && existing_artifacts="$((existing_artifacts + 1))"
    [[ -s "${artifact}" ]] && complete_artifacts="$((complete_artifacts + 1))"
done

if (( complete_artifacts == ${#artifact_suffixes[@]} )); then
    echo "SKIP: all output artifacts already exist and are non-empty for ${dataset_name}"
    exit 0
fi
if (( existing_artifacts > 0 )); then
    echo "ERROR: partial AlphaGenome output set exists for ${dataset_name}." >&2
    echo "Remove or archive the partial files before resubmitting this task." >&2
    exit 1
fi

mkdir -p "${output_dir}"

echo "Dataset: ${dataset_name}"
echo "Training positive FASTA: ${train_pos}"
echo "Training negative FASTA: ${train_neg}"
echo "Validation positive FASTA: ${validation_pos}"
echo "Validation negative FASTA: ${validation_neg}"
echo "Held-out test positive FASTA: ${test_pos}"
echo "Held-out test negative FASTA: ${test_neg}"
echo "Training positive sequences: $(grep -c '^>' "${train_pos}")"
echo "Training negative sequences: $(grep -c '^>' "${train_neg}")"
echo "Validation positive sequences: $(grep -c '^>' "${validation_pos}")"
echo "Validation negative sequences: $(grep -c '^>' "${validation_neg}")"
echo "Test positive sequences: $(grep -c '^>' "${test_pos}")"
echo "Test negative sequences: $(grep -c '^>' "${test_neg}")"
echo "Output prefix: ${output_prefix}"

################################################################################
### Frozen production configuration
################################################################################
# Every task starts a new Python process and therefore a newly initialized head.
# Validation AUPRC selects the best epoch; test FASTAs are loaded only afterward.
"${python_bin}" src/train_alphagenome.py \
    "${train_pos}" \
    "${train_neg}" \
    "${validation_pos}" \
    "${validation_neg}" \
    "${test_pos}" \
    "${test_neg}" \
    "${output_prefix}" \
    --weights "${weights}" \
    --sequence-length 256 \
    --hidden-size 1024 \
    --dropout 0.1 \
    --batch-size 32 \
    --num-workers "${SLURM_CPUS_PER_TASK}" \
    --epochs 10 \
    --learning-rate 0.001 \
    --weight-decay 0.0 \
    --patience 5 \
    --seed 42 \
    --threshold 0.5 \
    --max-shift 15 \
    --rc-prob 0.5 \
    --shift-prob 0.5

for suffix in "${artifact_suffixes[@]}"; do
    artifact="${output_prefix}${suffix}"
    if [[ ! -s "${artifact}" ]]; then
        echo "ERROR: trainer did not produce a non-empty artifact: ${artifact}" >&2
        exit 1
    fi
done

end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"
echo "Training PASS: ${dataset_name}"
echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
