#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

background_type="${1:-}"
case "${background_type}" in
  real|synthetic) ;;
  *) echo "ERROR: background type must be real or synthetic" >&2; exit 2 ;;
esac

metrics="results/comparison/GRCh38/${background_type}/alphagenome_native_metrics.tsv"
outdir="reports/rendered/comparison/${background_type}"
python3 src/prepare_cross_task_metrics.py \
  --background-type "${background_type}" \
  --models alphagenome \
  --native-only \
  --output "${metrics}"

mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/native_comparison_chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],negative_kind=a[3],full_model_prefix=""),output_file="alphagenome_native_evn_vs_svn_GRCh38_chromosome_learning.html",output_dir=a[2],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${outdir}" "${background_type}"
