#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

analysis="${1:-}"
case "${analysis}" in
  evn) analysis_name="enhancer-versus-negative"; positive_class="Enhancer" ;;
  svn) analysis_name="silencer-versus-negative"; positive_class="Silencer" ;;
  *) echo "ERROR: analysis must be evn or svn" >&2; exit 2 ;;
esac

alphagenome_real_metrics="results/${analysis}/GRCh38/alphagenome/real/chromosome_learning_curves/metrics.tsv"
alphagenome_synthetic_metrics="results/${analysis}/GRCh38/alphagenome/synthetic/chromosome_learning_curves/metrics.tsv"
lsgkm_real_metrics="results/${analysis}/GRCh38/ls-gkm/real/chromosome_learning_curves/metrics.tsv"
lsgkm_synthetic_metrics="results/${analysis}/GRCh38/ls-gkm/synthetic/chromosome_learning_curves/metrics.tsv"
outdir="reports/rendered/${analysis}/real_vs_synthetic"

for file in "${alphagenome_real_metrics}" "${alphagenome_synthetic_metrics}" \
            "${lsgkm_real_metrics}" "${lsgkm_synthetic_metrics}"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty metrics table: ${file}" >&2; exit 1; }
done

mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/alphagenome_background_comparison_GRCh38_results.Rmd",params=list(alphagenome_real_metrics_path=a[1],alphagenome_synthetic_metrics_path=a[2],lsgkm_real_metrics_path=a[3],lsgkm_synthetic_metrics_path=a[4],analysis_name=a[5],positive_class=a[6]),output_file="real_vs_synthetic_GRCh38_chromosome_learning.html",output_dir=a[7],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${alphagenome_real_metrics}" "${alphagenome_synthetic_metrics}" "${lsgkm_real_metrics}" "${lsgkm_synthetic_metrics}" "${analysis_name}" "${positive_class}" "${outdir}"
