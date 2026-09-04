#!/usr/bin/env bash
set -euo pipefail

analysis="${1:-evn}"
case "${analysis}" in
    evn)
        export LSGKM_RESULTS_DIR="results/evn/GRCh38/ls-gkm"
        export LSGKM_TRAINING_DIR="data/processed/evn/GRCh38/learning_curves"
        export LSGKM_ANALYSIS_NAME="enhancer-versus-negative"
        export LSGKM_POSITIVE_CLASS="Enhancer"
        export LSGKM_OUTPUT_DIR="reports/rendered/evn"
        export LSGKM_OUTPUT_FILE="lsgkm_GRCh38_results.html"
        ;;
    svn)
        export LSGKM_RESULTS_DIR="results/svn/GRCh38/ls-gkm"
        export LSGKM_TRAINING_DIR="data/processed/svn/GRCh38/learning_curves"
        export LSGKM_ANALYSIS_NAME="silencer-versus-negative"
        export LSGKM_POSITIVE_CLASS="Silencer"
        export LSGKM_OUTPUT_DIR="reports/rendered/svn"
        export LSGKM_OUTPUT_FILE="lsgkm_GRCh38_results_svn.html"
        ;;
    *)
        echo "ERROR: analysis must be 'evn' or 'svn', received: ${analysis}" >&2
        exit 2
        ;;
esac

Rscript -e 'rmarkdown::render(
  input = "reports/lsgkm_GRCh38_results.Rmd",
  params = list(
    results_dir = Sys.getenv("LSGKM_RESULTS_DIR"),
    training_dir = Sys.getenv("LSGKM_TRAINING_DIR"),
    analysis_name = Sys.getenv("LSGKM_ANALYSIS_NAME"),
    positive_class = Sys.getenv("LSGKM_POSITIVE_CLASS")
  ),
  output_file = Sys.getenv("LSGKM_OUTPUT_FILE"),
  output_dir = Sys.getenv("LSGKM_OUTPUT_DIR"),
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv())
)'
