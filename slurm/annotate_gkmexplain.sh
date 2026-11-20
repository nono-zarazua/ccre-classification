#!/usr/bin/env bash
## Annotate existing GkmExplain TF-MoDISco patterns; does not rerun discovery.
#SBATCH --job-name=annotate_gkmexplain
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/annotate_gkmexplain.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/annotate_gkmexplain.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-2%2
set -euo pipefail

module load lang/Miniforge3/24.1.2-0
set +u
source activate ls-gkm
set -u

tasks=(evn svn evs)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 0 && task_id < ${#tasks[@]} )) || {
  echo "ERROR: SLURM_ARRAY_TASK_ID must be 0-2" >&2
  exit 1
}
task="${tasks[task_id]}"
discovery="results/${task}/GRCh38/ls-gkm/real/gkmexplain/through_chrX"
output="results/${task}/GRCh38/ls-gkm/real/gkmexplain_hocomoco14/through_chrX"

"${CONDA_PREFIX}/bin/python" src/annotate_gkmexplain_patterns.py \
  --task "${task}" \
  --positive-modisco "${discovery}/positive_modisco.h5" \
  --negative-modisco "${discovery}/negative_modisco.h5" \
  --motif-meme data/motifs/hocomoco14/H14CORE-CLUSTERED.meme \
  --motif-metadata data/motifs/hocomoco14/cluster_list.tsv \
  --output-dir "${output}" \
  --n-jobs "${SLURM_CPUS_PER_TASK}" \
  --fdr 0.05 \
  --overwrite
