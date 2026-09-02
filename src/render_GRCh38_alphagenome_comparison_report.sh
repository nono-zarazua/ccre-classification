#!/usr/bin/env bash
set -euo pipefail

Rscript -e 'rmarkdown::render(
  input = "reports/alphagenome_GRCh38_evn_vs_svn.Rmd",
  params = list(
    evn_results_dir = "results/evn/GRCh38/alphagenome",
    svn_results_dir = "results/svn/GRCh38/alphagenome"
  ),
  output_file = "alphagenome_GRCh38_evn_vs_svn.html",
  output_dir = "reports/rendered/comparison",
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv())
)'
