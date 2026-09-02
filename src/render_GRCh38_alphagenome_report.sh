#!/usr/bin/env bash
set -euo pipefail

analysis="${1:-evn}"
case "${analysis}" in
    evn)
        export AG_RESULTS_DIR="results/evn/GRCh38/alphagenome"
        export AG_ANALYSIS_NAME="enhancer-versus-negative"
        export AG_POSITIVE_CLASS="Enhancer"
        export AG_OUTPUT_DIR="reports/rendered/evn"
        ;;
    svn)
        export AG_RESULTS_DIR="results/svn/GRCh38/alphagenome"
        export AG_ANALYSIS_NAME="silencer-versus-negative"
        export AG_POSITIVE_CLASS="Silencer"
        export AG_OUTPUT_DIR="reports/rendered/svn"
        ;;
    *)
        echo "ERROR: analysis must be 'evn' or 'svn', received: ${analysis}" >&2
        exit 2
        ;;
esac

Rscript -e 'rmarkdown::render(
  input = "reports/alphagenome_GRCh38_results.Rmd",
  params = list(
    results_dir = Sys.getenv("AG_RESULTS_DIR"),
    analysis_name = Sys.getenv("AG_ANALYSIS_NAME"),
    positive_class = Sys.getenv("AG_POSITIVE_CLASS")
  ),
  output_file = "alphagenome_GRCh38_results.html",
  output_dir = Sys.getenv("AG_OUTPUT_DIR"),
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv())
)'
