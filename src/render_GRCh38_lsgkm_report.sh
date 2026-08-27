#!/bin/bash
set -euo pipefail

Rscript -e 'rmarkdown::render(
  input = "reports/lsgkm_GRCh38_results.Rmd",
  params = list(
    results_dir = "results/evn/GRCh38/ls-gkm",
    training_dir = "data/processed/evn/GRCh38/learning_curves"
  ),
  output_file = "lsgkm_GRCh38_results.html",
  output_dir = "reports/rendered/evn",
  knit_root_dir = getwd()
)'
