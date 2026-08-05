#!/bin/bash

Rscript -e 'rmarkdown::render(
  input = "reports/lsgkm_results.Rmd",
  params = list(
    results_dir = "results/evn/ls-gkm"
  ),
  output_file = "lsgkm_results.html",
  output_dir = "reports/rendered/evn",
  knit_root_dir = getwd()
)'
