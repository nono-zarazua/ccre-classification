#!/usr/bin/env bash
# Render the dedicated EVN-versus-SVN LS-GKM learning-curve comparison.
set -euo pipefail

Rscript -e 'rmarkdown::render(
  input = "reports/lsgkm_GRCh38_evn_vs_svn.Rmd",
  params = list(
    evn_results_dir = "results/evn/GRCh38/ls-gkm",
    svn_results_dir = "results/svn/GRCh38/ls-gkm",
    evn_training_dir = "data/processed/evn/GRCh38/learning_curves",
    svn_training_dir = "data/processed/svn/GRCh38/learning_curves"
  ),
  output_file = "lsgkm_GRCh38_evn_vs_svn.html",
  output_dir = "reports/rendered/comparison",
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv())
)'
