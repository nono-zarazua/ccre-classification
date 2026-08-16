#!/bin/bash
## Written by Nono Zarazua
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=negs
#SBATCH --comment=
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/negs.%A_%a.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classification/slurmlogs/negs.%A_%a.err
#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
#SBATCH --array=1-24
################################################################################
### Safer Bash behavior
################################################################################
set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
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

module load lang/Miniforge3/24.1.2-0
source activate ls-gkm

echo "Conda environment: ${CONDA_PREFIX}"
echo "Rscript: $(command -v Rscript)"
ls -lh src/build_negs.R

task_id="${SLURM_ARRAY_TASK_ID}"
if (( task_id >= 1 && task_id <= 22 )); then
    i="${task_id}"
elif (( task_id == 23 )); then
    i="X"
elif (( task_id == 24 )); then
    i="Y"
else
    echo "ERROR: SLURM_ARRAY_TASK_ID must be between 1 and 24" >&2
    exit 1
fi

input_bed="data/processed/evn/GRCh38/chroms/chrom${i}/chrom${i}_GRCh38_ELS.bed"
negative_bed="data/processed/evn/GRCh38/chroms/chrom${i}/neg1x_chrom${i}_GRCh38_ELS.bed"

echo "Array task ID: ${task_id}"
echo "Chromosome: ${i}"
echo "Input BED: ${input_bed}"
echo "Expected negative BED: ${negative_bed}"

Rscript --vanilla src/build_negs.R "${input_bed}"

end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
