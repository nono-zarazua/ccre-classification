#!/usr/bin/env bash
set -euo pipefail
Rscript -e 'rmarkdown::render("reports/cross_task_GRCh38_results.Rmd", output_file="cross_task_GRCh38_results.html", output_dir="reports/rendered/comparison", knit_root_dir=getwd())'
