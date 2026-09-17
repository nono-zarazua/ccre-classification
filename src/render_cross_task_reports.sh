#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"
metrics="results/cross_task/GRCh38/native_vs_cross_metrics.tsv"
outdir="reports/rendered/comparison"
python3 src/prepare_cross_task_metrics.py --output "${metrics}"
mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/cross_task_chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],negative_kind="genomic"),output_file="native_vs_cross_GRCh38_chromosome_learning.html",output_dir=a[2],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${outdir}"
