#!/usr/bin/env bash
set -euo pipefail

export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

outdir="reports/rendered/feature-importance"
mkdir -p "${outdir}"

Rscript -e 'rmarkdown::render(
  "reports/gkmexplain_GRCh38_results.Rmd",
  output_file="lsgkm_gkmexplain_GRCh38.html",
  output_dir="reports/rendered/feature-importance",
  knit_root_dir=getwd(),
  envir=new.env(parent=globalenv())
)'
