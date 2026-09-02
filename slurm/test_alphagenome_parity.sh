#!/usr/bin/env bash
## GPU parity smoke test for src/train_alphagenome.py.
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=test_ag_parity
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/test_alphagenome_parity.%j.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/test_alphagenome_parity.%j.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"

module load lang/Miniforge3/24.1.2-0
# Conda activation hooks may inspect optional variables that are unset.
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

train_records_per_class=256
evaluation_records_per_class=128
parity_tmp="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/alphagenome-parity.XXXXXX")"
echo "Temporary FASTA directory: ${parity_tmp}"

# FASTA records may span multiple lines. Stop only when the header following the
# requested final record is reached, thereby preserving every selected record.
take_fasta_prefix() {
    local input_fasta="$1"
    local output_fasta="$2"
    local records="$3"
    awk -v limit="${records}" '
        /^>/ {
            seen++
            if (seen > limit) exit
        }
        { print }
    ' "${input_fasta}" > "${output_fasta}"

    local observed
    observed="$(grep -c '^>' "${output_fasta}")"
    if [[ "${observed}" -ne "${records}" ]]; then
        echo "ERROR: expected ${records} records in ${output_fasta}, found ${observed}" >&2
        exit 1
    fi
}

take_fasta_prefix \
    data/processed/evn/GRCh38/learning_curves/train_10_pos.fa \
    "${parity_tmp}/train_pos.fa" "${train_records_per_class}"
take_fasta_prefix \
    data/processed/evn/GRCh38/learning_curves/train_10_neg.fa \
    "${parity_tmp}/train_neg.fa" "${train_records_per_class}"
take_fasta_prefix \
    data/processed/evn/GRCh38/training/validation_pos.fa \
    "${parity_tmp}/validation_pos.fa" "${evaluation_records_per_class}"
take_fasta_prefix \
    data/processed/evn/GRCh38/training/validation_neg.fa \
    "${parity_tmp}/validation_neg.fa" "${evaluation_records_per_class}"
take_fasta_prefix \
    data/processed/evn/GRCh38/training/test_pos.fa \
    "${parity_tmp}/test_pos.fa" "${evaluation_records_per_class}"
take_fasta_prefix \
    data/processed/evn/GRCh38/training/test_neg.fa \
    "${parity_tmp}/test_neg.fa" "${evaluation_records_per_class}"

output_prefix="models/evn/GRCh38/alphagenome/parity/parity_smoke"

"${python_bin}" src/train_alphagenome.py \
    "${parity_tmp}/train_pos.fa" \
    "${parity_tmp}/train_neg.fa" \
    "${parity_tmp}/validation_pos.fa" \
    "${parity_tmp}/validation_neg.fa" \
    "${parity_tmp}/test_pos.fa" \
    "${parity_tmp}/test_neg.fa" \
    "${output_prefix}" \
    --epochs 3 \
    --patience 2 \
    --batch-size 32 \
    --num-workers "${SLURM_CPUS_PER_TASK}" \
    --seed 42 \
    --overwrite

required_suffixes=(
    .best_head.pt
    .config.json
    .history.tsv
    .validation_metrics.json
    .test_metrics.json
)
for suffix in "${required_suffixes[@]}"; do
    artifact="${output_prefix}${suffix}"
    if [[ ! -s "${artifact}" ]]; then
        echo "ERROR: missing or empty parity artifact: ${artifact}" >&2
        exit 1
    fi
done

echo "Parity test PASS"
echo "Validation metrics:"
cat "${output_prefix}.validation_metrics.json"
echo "Test metrics:"
cat "${output_prefix}.test_metrics.json"
echo "End time: $(date)"
