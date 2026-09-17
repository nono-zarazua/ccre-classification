#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"
metrics="results/cross_task/GRCh38/native_vs_cross_metrics.tsv"
full_prefix="models/evn/GRCh38/alphagenome/full_benchmark/full_evn"
outdir="reports/rendered/comparison"
[[ -s "${metrics}" ]] || { echo "ERROR: missing/empty metrics table: ${metrics}" >&2; exit 1; }
mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/native_comparison_chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],negative_kind="genomic",full_model_prefix=a[3]),output_file="native_evn_vs_svn_GRCh38_chromosome_learning.html",output_dir=a[2],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${outdir}" "${full_prefix}"
