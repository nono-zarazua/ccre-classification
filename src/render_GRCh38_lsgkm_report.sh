#!/usr/bin/env bash
set -euo pipefail
export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"
analysis="${1:-evn}"
background_type="${2:-}"
case "${analysis}" in
  evn) analysis_name="enhancer-versus-negative"; positive_class="Enhancer"; negative_class="Background sequence" ;;
  svn) analysis_name="silencer-versus-negative"; positive_class="Silencer"; negative_class="Background sequence" ;;
  evs) analysis_name="enhancer-versus-silencer"; positive_class="Enhancer"; negative_class="Silencer" ;;
  *) echo "ERROR: analysis must be evn, svn, or evs" >&2; exit 2 ;;
esac
case "${background_type}" in
  real|synthetic) ;;
  *) echo "ERROR: background type must be real or synthetic" >&2; exit 2 ;;
esac
[[ "${analysis}" != evs || "${background_type}" == real ]] || { echo "ERROR: EVS is real only" >&2; exit 2; }
metrics="results/${analysis}/GRCh38/ls-gkm/${background_type}/chromosome_learning_curves/metrics.tsv"
if [[ "${analysis}" == evs ]]; then outdir="reports/rendered/evs"; else outdir="reports/rendered/${analysis}/${background_type}"; fi
mkdir -p "${outdir}"
Rscript -e 'a<-commandArgs(TRUE); rmarkdown::render("reports/chromosome_learning_GRCh38_results.Rmd",params=list(metrics_path=a[1],analysis_name=a[2],positive_class=a[3],negative_class=a[6],model_name="LS-GKM",negative_kind=a[5],full_model_prefix=""),output_file="lsgkm_GRCh38_chromosome_learning.html",output_dir=a[4],knit_root_dir=getwd(),envir=new.env(parent=globalenv()))' "${metrics}" "${analysis_name}" "${positive_class}" "${outdir}" "${background_type}" "${negative_class}"
