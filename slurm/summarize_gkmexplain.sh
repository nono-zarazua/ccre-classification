#!/usr/bin/env bash
## Summarize final-model GkmExplain contributions with TF-MoDISco and k-mer ranks.
#SBATCH --job-name=summarize_gkmexplain
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/summarize_gkmexplain.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/summarize_gkmexplain.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2%1
set -euo pipefail

module load lang/Miniforge3/24.1.2-0
set +u
source activate ls-gkm
set -u

tasks=(evn svn evs)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 0 && task_id < 3 )) || { echo "ERROR: task must be 0-2" >&2; exit 1; }
dataset="${tasks[task_id]}"

"${CONDA_PREFIX}/bin/python" src/gkmexplain_analysis.py summarize \
  --task "${dataset}" \
  --model "models/${dataset}/GRCh38/ls-gkm/real/chromosome_learning_curves/through_chrX/through_chrX.model.txt" \
  --positive-fasta "data/processed/${dataset}/GRCh38/real/test_pos.fa" \
  --negative-fasta "data/processed/${dataset}/GRCh38/real/test_neg.fa" \
  --raw-dir "predictions/${dataset}/GRCh38/ls-gkm/real/gkmexplain/through_chrX" \
  --motif-meme "data/motifs/JASPAR2026_CORE_non-redundant_pfms_homo_sapiens_meme.txt" \
  --output-dir "results/${dataset}/GRCh38/ls-gkm/real/gkmexplain/through_chrX" \
  --n-jobs "${SLURM_CPUS_PER_TASK}" \
  --score-shards 8 \
  --overwrite
