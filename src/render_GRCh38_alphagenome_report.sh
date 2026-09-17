#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"
analysis="${1:-evn}"
case "${analysis}" in
  evn) analysis_name="enhancer-versus-negative"; positive_class="Enhancer" ;;
  svn) analysis_name="silencer-versus-negative"; positive_class="Silencer" ;;
  *) echo "ERROR: analysis must be evn or svn" >&2; exit 2 ;;
esac
metrics="results/${analysis}/GRCh38/alphagenome/chromosome_learning_curves/metrics.tsv"
full_prefix=""
if [[ "${analysis}" == evn ]]; then
  full_prefix="models/evn/GRCh38/alphagenome/full_benchmark/full_evn"
fi
outdir="reports/rendered/${analysis}"
mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],analysis_name=a[2],positive_class=a[3],model_name="AlphaGenome",negative_kind="genomic",full_model_prefix=a[5]),output_file="alphagenome_GRCh38_chromosome_learning.html",output_dir=a[4],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${analysis_name}" "${positive_class}" "${outdir}" "${full_prefix}"
