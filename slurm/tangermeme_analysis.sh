#!/usr/bin/env bash
## Clustered motif-importance analysis for final real-background EVN/SVN/EVS models.
#SBATCH --job-name=tangermeme
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/tangermeme.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/tangermeme.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100_40gb:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2%1
set -euo pipefail

echo "Job ${SLURM_JOB_ID}, task ${SLURM_ARRAY_TASK_ID}, node $(hostname), start $(date)"
module load lang/Miniforge3/24.1.2-0
set +u
source activate alphagenome
set -u

python_bin="${CONDA_PREFIX}/bin/python"
"${python_bin}" - <<'PY'
import torch
import alphagenome_pytorch
import tangermeme
assert torch.cuda.is_available(), "CUDA unavailable"
print("PyTorch:", torch.__version__)
print("AlphaGenome:", alphagenome_pytorch.__version__)
print("Tangermeme:", tangermeme.__version__)
print("GPU:", torch.cuda.get_device_name(0))
PY

tasks=(evn svn evs)
task_id="${SLURM_ARRAY_TASK_ID}"
(( task_id >= 0 && task_id < ${#tasks[@]} )) || {
  echo "ERROR: SLURM_ARRAY_TASK_ID must be 0-2" >&2
  exit 1
}
task="${tasks[task_id]}"
seqlet_selection="${SEQLET_SELECTION:-recursive}"
case "${seqlet_selection}" in
  fixed|recursive|modisco) ;;
  *) echo "ERROR: SEQLET_SELECTION must be fixed, recursive, or modisco" >&2; exit 1 ;;
esac

checkpoint="models/${task}/GRCh38/alphagenome/real/chromosome_learning_curves/through_chrX/through_chrX.best_head.pt"
positive="data/processed/${task}/GRCh38/real/test_pos.fa"
negative="data/processed/${task}/GRCh38/real/test_neg.fa"
weights="weights/fold_0_weights.safetensors"
motifs="data/motifs/hocomoco14/H14CORE-CLUSTERED.meme"
motif_metadata="data/motifs/hocomoco14/cluster_list.tsv"
if [[ "${seqlet_selection}" == "fixed" ]]; then
  output="results/${task}/GRCh38/alphagenome/real/tangermeme_clustered/through_chrX"
else
  output="results/${task}/GRCh38/alphagenome/real/tangermeme_clustered_${seqlet_selection}/through_chrX"
fi

for input in "${checkpoint}" "${positive}" "${negative}" "${weights}" "${motifs}" "${motif_metadata}"; do
  [[ -s "${input}" ]] || { echo "ERROR: missing/empty ${input}" >&2; exit 1; }
done

"${python_bin}" src/tangermeme_analysis.py \
  --task "${task}" \
  --checkpoint "${checkpoint}" \
  --weights "${weights}" \
  --positive-fasta "${positive}" \
  --negative-fasta "${negative}" \
  --motif-meme "${motifs}" \
  --motif-metadata "${motif_metadata}" \
  --output-dir "${output}" \
  --all-sequences \
  --seed 42 \
  --method auto \
  --batch-size 1 \
  --n-shuffles 10 \
  --deeplift-tolerance 0.001 \
  --audit-sequences 4 \
  --audit-shuffles 5 \
  --tomtom-pvalue 0.05 \
  --n-jobs "${SLURM_CPUS_PER_TASK}" \
  --seqlet-selection "${seqlet_selection}" \
  --candidate-window-size 21 \
  --candidate-fraction 0.02 \
  --max-candidates-per-sign 3 \
  --max-candidate-sites 30000 \
  --controls-per-candidate 4 \
  --min-known-fraction 1.0 \
  --ablation-shuffles 20 \
  --ablation-site-batch 64 \
  --region-fdr 0.05 \
  --control-quantile 0.95 \
  --motif-fdr 0.05 \
  --bootstrap-replicates 1000 \
  --n-ism-sequences 40 \
  --ism-batch-size 32 \
  --n-logo-sequences 8 \
  --logo-flank 50 \
  --overwrite

echo "Tangermeme analysis PASS: ${task}; end $(date)"
