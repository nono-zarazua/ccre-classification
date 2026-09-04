#!/usr/bin/env bash
set -euo pipefail

Rscript -e 'rmarkdown::render("reports/cross_task_method_GRCh38_results.Rmd", params=list(metrics="results/cross_task/GRCh38/native_vs_cross_metrics.tsv",model="alphagenome",feature_dir="results/cross_task/GRCh38/alphagenome_features"), output_file="alphagenome_GRCh38_native_vs_cross.html",output_dir="reports/rendered/comparison",knit_root_dir=getwd())'
Rscript -e 'rmarkdown::render("reports/cross_task_method_GRCh38_results.Rmd", params=list(metrics="results/cross_task/GRCh38/native_vs_cross_metrics.tsv",model="lsgkm",feature_dir=NULL), output_file="lsgkm_GRCh38_native_vs_cross.html",output_dir="reports/rendered/comparison",knit_root_dir=getwd())'
Rscript -e 'rmarkdown::render("reports/cross_task_GRCh38_results.Rmd", output_file="cross_task_GRCh38_results.html",output_dir="reports/rendered/comparison",knit_root_dir=getwd())'
