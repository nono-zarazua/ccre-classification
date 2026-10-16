#!/usr/bin/env bash
## Independent EVN/SVN/EVS frozen-encoder AlphaGenome cumulative chromosome learning curves.
## Select DATASET=evn|svn|evs at submission.
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
if [[ "${dataset}" != "evn" && "${dataset}" != "svn" && "${dataset}" != "evs" ]]; then
  echo "ERROR: DATASET must be evn, svn, or evs; received: ${dataset}" >&2
  exit 1
fi
if [[ "${type}" != "real" && "${type}" != "synthetic" ]]; then
  echo "ERROR: TYPE must be real or synthetic; received: ${type}" >&2
  exit 1
fi
if [[ "${dataset}" == "evs" && "${type}" != "real" ]]; then
  echo "ERROR: EVS is available only under TYPE=real" >&2
  exit 1
fi


task_id="${SLURM_ARRAY_TASK_ID}"
chromosomes=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
(( task_id >= 1 && task_id <= ${#chromosomes[@]} )) || { echo "ERROR: task must be 1-${#chromosomes[@]}" >&2; exit 1; }
endpoint="${chromosomes[task_id-1]}"; name="through_chr${endpoint}"
root="data/processed/${dataset}/GRCh38/${type}"
train_100_pos="${root}/train_100_pos.fa"; train_100_neg="${root}/train_100_neg.fa"
validation_pos="${root}/validation_pos.fa"; validation_neg="${root}/validation_neg.fa"
test_pos="${root}/test_pos.fa"; test_neg="${root}/test_neg.fa"
weights="weights/fold_0_weights.safetensors"
out_dir="models/${dataset}/GRCh38/alphagenome/${type}/chromosome_learning_curves/${name}"
prefix="${out_dir}/${name}"
for f in "${train_100_pos}" "${train_100_neg}" "${validation_pos}" "${validation_neg}" "${test_pos}" "${test_neg}" "${weights}"; do
  [[ -s "${f}" ]] || { echo "ERROR: missing/empty ${f}" >&2; exit 1; }
done
suffixes=(.best_head.pt .config.json .history.tsv .validation_metrics.json .test_metrics.json .train_100_metrics.json)
existing=0; complete=0
for s in "${suffixes[@]}"; do [[ -e "${prefix}${s}" ]] && ((existing+=1)); [[ -s "${prefix}${s}" ]] && ((complete+=1)); done
(( complete != ${#suffixes[@]} )) || { echo "SKIP: complete outputs exist for ${name}"; exit 0; }
(( existing == 0 )) || { echo "ERROR: partial outputs exist for ${name}" >&2; exit 1; }

work_dir="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/train_ag.${SLURM_JOB_ID}.${task_id}.XXXXXX")"
trap 'rm -rf -- "${work_dir}"' EXIT
pos="${work_dir}/train_pos.fa"; neg="${work_dir}/train_neg.fa"; : >"${pos}"; : >"${neg}"
included=()
for ((i=0;i<task_id;i++)); do
  chrom="${chromosomes[i]}"; shard="${root}/chromosomes/chrom${chrom}"
  [[ -s "${shard}/train_pos.fa" && -s "${shard}/train_neg.fa" ]] || { echo "ERROR: missing chr${chrom} shard" >&2; exit 1; }
  cat "${shard}/train_pos.fa" >>"${pos}"; cat "${shard}/train_neg.fa" >>"${neg}"; included+=("chr${chrom}")
done
n_pos="$(grep -c '^>' "${pos}")"; n_neg="$(grep -c '^>' "${neg}")"
(( n_pos > 0 && n_pos == n_neg )) || { echo "ERROR: unbalanced cumulative input" >&2; exit 1; }
mkdir -p "${out_dir}"
echo "Task dataset: ${dataset}; background: ${type}; checkpoint: ${name}; chromosomes: ${included[*]}; observations: $((n_pos+n_neg))"
echo "Validation chr1; test chr2; output ${prefix}"
"${python_bin}" src/train_alphagenome.py "${pos}" "${neg}" "${validation_pos}" "${validation_neg}" "${test_pos}" "${test_neg}" "${prefix}" \
  --weights "${weights}" --sequence-length 256 --hidden-size 1024 --dropout 0.1 \
  --batch-size 32 --num-workers "${SLURM_CPUS_PER_TASK}" --epochs 10 \
  --dataset-name "${dataset}" --background-type "${type}" --checkpoint-name "${name}" \
  --training-chromosomes "$(IFS=,; echo "${included[*]}")" \
  --learning-rate 0.001 --weight-decay 0.0 --patience 5 --seed 42 \
  --threshold 0.5 --max-shift 15 --rc-prob 0.5 --shift-prob 0.5
"${python_bin}" src/predict_alphagenome.py "${prefix}.best_head.pt" "${train_100_pos}" "${train_100_neg}" \
  --metrics-json "${prefix}.train_100_metrics.json" --weights "${weights}" \
  --batch-size 32 --num-workers "${SLURM_CPUS_PER_TASK}"
for s in "${suffixes[@]}"; do [[ -s "${prefix}${s}" ]] || { echo "ERROR: missing ${prefix}${s}" >&2; exit 1; }; done
elapsed="$(($(date -u +%s)-start_time))"
echo "Training PASS: ${name}; elapsed $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
