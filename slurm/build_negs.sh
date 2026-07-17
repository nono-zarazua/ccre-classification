#!/bin/bash
## Written by Nono Zarazua
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name=negs
#SBATCH --comment=
### changedir / output
#SBATCH --chdir=/work/zarazuanav/workspace/repos/ccre-classification
#SBATCH --output=/work/zarazuanav/workspace/repos/ccre-classificatino/slurmlogs/negs.%j.out
#SBATCH --error=/work/zarazuanav/workspace/repos/ccre-classificatino/slurmlogs/negs.%j.err

#SBATCH --mail-user=zarazuanav@uni-potsdam.de
#SBATCH --mail-type=BEGIN,END,FAIL
#
### Node / CPU / Memory Settings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=05-00:00:00
# --qos=long
################################################################################
### Safer Bash behavior
################################################################################
set -euo pipefail
################################################################################
### Job Info
################################################################################
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
###############################################################################
### Load modules / source environment
################################################################################
module load lang/Miniforge3/24.1.2-0
source activate ls-gkm
################################################################################
### Running code
################################################################################
echo "Conda environment: ${CONDA_PREFIX}"
echo "Rscript: $(command -v Rscript)"
ls -lh build_negs.R
Rscript --vanilla src/build_negs.R
################################################################################
### Runtime
################################################################################
end_time="$(date -u +%s)"
elapsed="$((end_time - start_time))"

echo "End time: $(date)"
echo "Elapsed time: $((elapsed / 3600))h $(((elapsed % 3600) / 60))m $((elapsed % 60))s"
