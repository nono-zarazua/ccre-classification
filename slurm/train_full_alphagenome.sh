#!/usr/bin/env bash
## Train the full-EVN AlphaGenome classifier on the common full benchmark.
#SBATCH --job-name=train_full_ag
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_full_alphagenome.%j.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_full_alphagenome.%j.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodelist=n-hpc-gz6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=05-00:00:00

set -euo pipefail
echo "Job ${SLURM_JOB_ID}, node $(hostname), start $(date)"

module load lang/Miniforge3/24.1.2-0
set +u
source activate alphagenome
set -u
python_bin="${CONDA_PREFIX}/bin/python"
[[ -x "${python_bin}" ]] || { echo "ERROR: AlphaGenome Python is missing" >&2; exit 1; }
"${python_bin}" -c 'import numpy, sklearn, torch, alphagenome_pytorch; assert torch.cuda.is_available(), "CUDA unavailable"; print(torch.__version__, torch.cuda.get_device_name(0))'


benchmark="data/processed/evn/GRCh38/full-benchmark"
train_pos="${benchmark}/positive/train_pos.fa"
validation_pos="${benchmark}/positive/validation_pos.fa"
test_pos="${benchmark}/positive/test_pos.fa"
train_neg="${benchmark}/train_neg.fa"
validation_neg="${benchmark}/validation_neg.fa"
test_neg="${benchmark}/test_neg.fa"
weights="weights/fold_0_weights.safetensors"
out_dir="models/evn/GRCh38/alphagenome/full_benchmark"
prefix="${out_dir}/full_evn"

for file in "${train_pos}" "${train_neg}" "${validation_pos}" \
            "${validation_neg}" "${test_pos}" "${test_neg}" "${weights}"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty ${file}" >&2; exit 1; }
done

suffixes=(.best_head.pt .config.json .history.tsv .validation_metrics.json .test_metrics.json)
existing=0
complete=0
for suffix in "${suffixes[@]}"; do
  [[ -e "${prefix}${suffix}" ]] && ((existing+=1))
  [[ -s "${prefix}${suffix}" ]] && ((complete+=1))
done
if (( complete == ${#suffixes[@]} )); then
  echo "Training outputs already complete; evaluating train_100 only"
elif (( existing != 0 )); then
  echo "ERROR: partial training outputs exist" >&2
  exit 1
else
  mkdir -p "${out_dir}"
  echo "Full EVN benchmark; training chr3-chr22,chrX; validation chr1; test chr2"
  "${python_bin}" src/train_alphagenome.py \
    "${train_pos}" "${train_neg}" "${validation_pos}" "${validation_neg}" \
    "${test_pos}" "${test_neg}" "${prefix}" \
    --weights "${weights}" --sequence-length 256 --hidden-size 1024 --dropout 0.1 \
    --batch-size 32 --num-workers 0 --epochs 10 \
    --learning-rate 0.001 --weight-decay 0 --patience 5 --seed 42 \
    --threshold 0.5 --max-shift 15 --rc-prob 0.5 --shift-prob 0.5 \
    --dataset-name evn --checkpoint-name "full_evn" \
    --training-chromosomes chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX
fi

for suffix in "${suffixes[@]}"; do
  [[ -s "${prefix}${suffix}" ]] || { echo "ERROR: missing ${prefix}${suffix}" >&2; exit 1; }
done

train_metrics="${prefix}.train_100_metrics.json"
if [[ -s "${train_metrics}" ]]; then
  echo "SKIP: train_100 metrics already exist: ${train_metrics}"
elif [[ -e "${train_metrics}" ]]; then
  echo "ERROR: empty train_100 metrics file exists: ${train_metrics}" >&2
  exit 1
else
  "${python_bin}" src/predict_alphagenome.py \
    "${prefix}.best_head.pt" "${train_pos}" "${train_neg}" \
    --metrics-json "${train_metrics}" --weights "${weights}" \
    --batch-size 32 --num-workers 0
fi
[[ -s "${train_metrics}" ]] || { echo "ERROR: missing ${train_metrics}" >&2; exit 1; }
echo "PASS: full EVN, including train_100 evaluation with genomic backgrounds; end $(date)"
