#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

background_type="${1:-}"
case "${background_type}" in
  real) ;;
  *) echo "ERROR: cross-task reports are currently restricted to real backgrounds" >&2; exit 2 ;;
esac

metrics="results/cross_task/GRCh38/${background_type}/alphagenome_native_vs_cross_metrics.tsv"
outdir="reports/rendered/cross-task/${background_type}"

python3 src/prepare_cross_task_metrics.py \
  --background-type "${background_type}" \
  --models alphagenome \
  --output "${metrics}"

mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/cross_task_chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],negative_kind=a[3]),output_file="alphagenome_native_vs_cross_GRCh38_chromosome_learning.html",output_dir=a[2],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${outdir}" "${background_type}"
