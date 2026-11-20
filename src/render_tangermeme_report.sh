#!/usr/bin/env bash
set -euo pipefail

export PATH="/home/zarazuanav/.conda/envs/ls-gkm/bin:${PATH}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

outdir="reports/rendered/feature-importance"
mkdir -p "${outdir}"

motif_logo_dir="$(mktemp -d /tmp/ccre-tangermeme-motif-logos.XXXXXX)"
trap 'rm -rf "${motif_logo_dir}"' EXIT

branch_complete() {
  local branch="$1" task file
  local required=(metadata.json audit.json candidate_region_tests.tsv
                  motif_cluster_summary.tsv motif_annotation_summary.tsv
                  ism_validation.tsv)
  for task in evn svn evs; do
    for file in "${required[@]}"; do
      [[ -s "results/${task}/GRCh38/alphagenome/real/${branch}/through_chrX/${file}" ]] || return 1
    done
  done
}

if branch_complete tangermeme_clustered; then
  Rscript -e 'rmarkdown::render(
  "reports/tangermeme_GRCh38_results.Rmd",
  output_file="alphagenome_tangermeme_GRCh38.html",
  output_dir="reports/rendered/feature-importance",
  params=list(
    result_branch="tangermeme_clustered",
    analysis_label="Fixed-window seqlet selection"
  ),
  knit_root_dir=getwd(),
  envir=new.env(parent=globalenv())
)'
else
  echo "Skipping incomplete legacy fixed-window Tangermeme results" >&2
fi

/home/zarazuanav/.conda/envs/alphagenome/bin/python \
  src/render_tangermeme_motif_logos.py \
  --result-branch tangermeme_clustered_recursive \
  --motif-meme data/motifs/hocomoco14/H14CORE-CLUSTERED.meme \
  --output-dir "${motif_logo_dir}"

if branch_complete tangermeme_clustered_modisco; then
  Rscript -e 'rmarkdown::render(
  "reports/tangermeme_GRCh38_results.Rmd",
  output_file="alphagenome_tangermeme_tfmodisco_GRCh38.html",
  output_dir="reports/rendered/feature-importance",
  params=list(
    result_branch="tangermeme_clustered_modisco",
    analysis_label="TF-MoDISco seqlet selection"
  ),
  knit_root_dir=getwd(),
  envir=new.env(parent=globalenv())
)'
else
  echo "Skipping incomplete legacy TF-MoDISco Tangermeme results" >&2
fi

MOTIF_LOGO_MANIFEST="${motif_logo_dir}/motif_logo_manifest.tsv" Rscript -e 'rmarkdown::render(
  "reports/tangermeme_GRCh38_results.Rmd",
  output_file="alphagenome_tangermeme_recursive_GRCh38.html",
  output_dir="reports/rendered/feature-importance",
  params=list(
    result_branch="tangermeme_clustered_recursive",
    analysis_label="Recursive variable-length seqlet selection",
    include_logos=TRUE,
    motif_logo_manifest=Sys.getenv("MOTIF_LOGO_MANIFEST")
  ),
  knit_root_dir=getwd(),
  envir=new.env(parent=globalenv())
)'
