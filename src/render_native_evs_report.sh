#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

alphagenome_metrics="results/evs/GRCh38/alphagenome/real/chromosome_learning_curves/metrics.tsv"
lsgkm_metrics="results/evs/GRCh38/ls-gkm/real/chromosome_learning_curves/metrics.tsv"
outdir="reports/rendered/evs"

mkdir -p "${outdir}"

Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/native_evs_chromosome_learning_GRCh38_results.Rmd",params=list(alphagenome_metrics_path=a[1],lsgkm_metrics_path=a[2]),output_file="native_evs_GRCh38_chromosome_learning.html",output_dir=a[3],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${alphagenome_metrics}" "${lsgkm_metrics}" "${outdir}"
