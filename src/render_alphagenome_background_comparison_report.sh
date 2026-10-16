#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

analysis="${1:-}"
case "${analysis}" in
  evn) analysis_name="enhancer-versus-negative"; positive_class="Enhancer" ;;
  svn) analysis_name="silencer-versus-negative"; positive_class="Silencer" ;;
  *) echo "ERROR: analysis must be evn or svn" >&2; exit 2 ;;
esac

real_metrics="results/${analysis}/GRCh38/alphagenome/real/chromosome_learning_curves/metrics.tsv"
synthetic_metrics="results/${analysis}/GRCh38/alphagenome/synthetic/chromosome_learning_curves/metrics.tsv"
outdir="reports/rendered/${analysis}/real_vs_synthetic"

for file in "${real_metrics}" "${synthetic_metrics}"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty metrics table: ${file}" >&2; exit 1; }
done

mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/alphagenome_background_comparison_GRCh38_results.Rmd",params=list(real_metrics_path=a[1],synthetic_metrics_path=a[2],analysis_name=a[3],positive_class=a[4]),output_file="alphagenome_real_vs_synthetic_GRCh38_chromosome_learning.html",output_dir=a[5],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${real_metrics}" "${synthetic_metrics}" "${analysis_name}" "${positive_class}" "${outdir}"
