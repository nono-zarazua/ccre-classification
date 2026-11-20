#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

dataset="${1:-}"
background_type="${2:-}"

case "${dataset}" in
  evn) analysis_name="Enhancer versus negative"; positive_class="Enhancer" ;;
  svn) analysis_name="Silencer versus negative"; positive_class="Silencer" ;;
  *) echo "ERROR: dataset must be evn or svn" >&2; exit 2 ;;
esac
case "${background_type}" in
  real|synthetic) ;;
  *) echo "ERROR: background type must be real or synthetic" >&2; exit 2 ;;
esac

alphagenome_metrics="results/${dataset}/GRCh38/alphagenome/${background_type}/chromosome_learning_curves/metrics.tsv"
lsgkm_metrics="results/${dataset}/GRCh38/ls-gkm/${background_type}/chromosome_learning_curves/metrics.tsv"
outdir="reports/rendered/${dataset}/${background_type}"

for file in "${alphagenome_metrics}" "${lsgkm_metrics}"; do
  [[ -s "${file}" ]] || { echo "ERROR: missing/empty metrics table: ${file}" >&2; exit 1; }
done

mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/native_task_chromosome_learning_GRCh38_results.Rmd",params=list(alphagenome_metrics_path=a[1],lsgkm_metrics_path=a[2],analysis_name=a[3],positive_class=a[4],negative_kind=a[5]),output_file="native_GRCh38_chromosome_learning.html",output_dir=a[6],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${alphagenome_metrics}" "${lsgkm_metrics}" "${analysis_name}" "${positive_class}" "${background_type}" "${outdir}"
