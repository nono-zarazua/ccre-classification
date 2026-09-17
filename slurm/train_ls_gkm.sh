#!/usr/bin/env bash
## Independent cumulative EVN/SVN chromosome learning curves.
## Select DATASET=evn|svn at submission.
#SBATCH --job-name=train_ls_gkm
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_ls_gkm.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/train_ls_gkm.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-13%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
set -euo pipefail
echo "Job ${SLURM_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}, node $(hostname), start $(date)"
start_time="$(date -u +%s)"
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm

dataset="${DATASET:-evn}"
if [[ "${dataset}" != "evn" && "${dataset}" != "svn" ]]; then
  echo "ERROR: DATASET must be evn or svn; received: ${dataset}" >&2
  exit 1
fi


task_id="${SLURM_ARRAY_TASK_ID}"
chromosomes=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
(( task_id >= 1 && task_id <= ${#chromosomes[@]} )) || { echo "ERROR: task must be 1-${#chromosomes[@]}" >&2; exit 1; }
endpoint="${chromosomes[task_id-1]}"
name="through_chr${endpoint}"
root="data/processed/${dataset}/GRCh38/training"
model_dir="models/${dataset}/GRCh38/ls-gkm/chromosome_learning_curves/${name}"
prefix="${model_dir}/${name}"
model="${prefix}.model.txt"

# EVN crosses 75% at task 14; SVN crosses it at task 13.
if [[ "${dataset}" == "evn" ]]; then high_start=14; else high_start=13; fi
if (( task_id < high_start )); then gkm_mem=64000; required_mb=100000; else gkm_mem=160000; required_mb=200000; fi
allocated_mb="${SLURM_MEM_PER_NODE:-0}"
if (( allocated_mb < required_mb )); then
  echo "ERROR: task ${task_id} needs ${required_mb} MB; got ${allocated_mb} MB." >&2
  echo "High group: sbatch --export=ALL,DATASET=${dataset} --array=${high_start}-21%5 --mem=200G slurm/train_ls_gkm.sh" >&2
  exit 1
fi
[[ ! -s "${model}" ]] || { echo "SKIP: ${model} exists"; exit 0; }
[[ ! -e "${model}" ]] || { echo "ERROR: empty model exists: ${model}" >&2; exit 1; }

work_dir="$(mktemp -d "${SLURM_TMPDIR:-/tmp}/train_lsgkm.${SLURM_JOB_ID}.${task_id}.XXXXXX")"
trap 'rm -rf -- "${work_dir}"' EXIT
pos="${work_dir}/train_pos.fa"; neg="${work_dir}/train_neg.fa"
: >"${pos}"; : >"${neg}"
included=()
for ((i=0;i<task_id;i++)); do
  chrom="${chromosomes[i]}"; shard="${root}/chromosomes/chrom${chrom}"
  [[ -s "${shard}/train_pos.fa" && -s "${shard}/train_neg.fa" ]] || { echo "ERROR: missing chr${chrom} shard" >&2; exit 1; }
  cat "${shard}/train_pos.fa" >>"${pos}"; cat "${shard}/train_neg.fa" >>"${neg}"
  included+=("chr${chrom}")
done
n_pos="$(grep -c '^>' "${pos}")"; n_neg="$(grep -c '^>' "${neg}")"
(( n_pos > 0 && n_pos == n_neg )) || { echo "ERROR: unbalanced cumulative input" >&2; exit 1; }
mkdir -p "${model_dir}"
echo "Task dataset: ${dataset}; checkpoint: ${name}; chromosomes: ${included[*]}"
echo "Observations: $((n_pos+n_neg)); per class: ${n_pos}; gkm -m ${gkm_mem}; Slurm ${allocated_mb} MB"
/work/zarazuanav/workspace/repos/lsgkm/bin/gkmtrain -T "${SLURM_CPUS_PER_TASK}" -m "${gkm_mem}" "${pos}" "${neg}" "${prefix}"
[[ -s "${model}" ]] || { echo "ERROR: model was not created" >&2; exit 1; }
elapsed="$(($(date -u +%s)-start_time))"
echo "End $(date); elapsed $((elapsed/3600))h $(((elapsed%3600)/60))m $((elapsed%60))s"
