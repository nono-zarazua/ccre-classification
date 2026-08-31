#!/bin/bash
set -euo pipefail

Rscript -e 'rmarkdown::render(
  input = "reports/lsgkm_GRCh38_results.Rmd",
  params = list(
    results_dir = "results/svn/GRCh38/ls-gkm",
    training_dir = "data/processed/svn/GRCh38/learning_curves"
  ),
  output_file = "lsgkm_GRCh38_results_svn.html",
  output_dir = "reports/rendered/svn",
  knit_root_dir = getwd()
)'
