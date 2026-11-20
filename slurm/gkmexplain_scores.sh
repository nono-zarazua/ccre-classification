#!/usr/bin/env bash
## Generate observed, hypothetical, and shuffled-control GkmExplain scores.
#SBATCH --job-name=gkmexplain
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/gkmexplain.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/gkmexplain.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-17%2
set -euo pipefail

module load lang/Miniforge3/24.1.2-0
set +u
source activate ls-gkm
set -u

tasks=(evn svn evs)
classes=(pos neg)
job_types=(observed hypothetical null)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 0 && task_id < 18 )) || { echo "ERROR: task must be 0-17" >&2; exit 1; }

dataset="${tasks[task_id / 6]}"
within=$((task_id % 6))
class_name="${classes[within / 3]}"
job_type="${job_types[within % 3]}"
root="data/processed/${dataset}/GRCh38/real"
model="models/${dataset}/GRCh38/ls-gkm/real/chromosome_learning_curves/through_chrX/through_chrX.model.txt"
output_dir="predictions/${dataset}/GRCh38/ls-gkm/real/gkmexplain/through_chrX"
input="${root}/test_${class_name}.fa"
mode=0
output="${output_dir}/test_${class_name}_mode0.txt"

if [[ "${job_type}" == hypothetical ]]; then
  mode=1
  output="${output_dir}/test_${class_name}_mode1.txt"
elif [[ "${job_type}" == null ]]; then
  null_fasta="${output_dir}/test_${class_name}_dinucleotide_shuffled.fa"
  "${CONDA_PREFIX}/bin/python" src/gkmexplain_analysis.py make-null \
    --input-fasta "${input}" --output-fasta "${null_fasta}" \
    --n-sequences 1000 --seed $((42 + task_id))
  input="${null_fasta}"
  output="${output_dir}/test_${class_name}_null_mode0.txt"
fi

bash src/run_gkmexplain_parallel.sh \
  "${input}" "${model}" "${output}" "${mode}" \
  "${SLURM_CPUS_PER_TASK}" \
  repos/lsgkm/bin/gkmexplain
