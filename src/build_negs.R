#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(gkmSVM)
	library(BSgenome.Hsapiens.UCSC.hg38.masked)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
	stop(
	     "Usage: Rscript src/build_negs.R <fold_number>\n",
	     "Example: Rscript src/build_negs.R 1"
	)
}

fold <- suppressWarnings(as.integer(args[1]))

if (is.na(fold) || fold < 1 || fold > 5) {
	stop("fold_number must be an integer from 1 to 5")
}

fold_dir <- file.path("data", "processed", "ls-gkm", "folds", paste0("fold", fold))

input_bed <- file.path(fold_dir, paste0("fold", fold, ".bed"))
pos_fa    <- file.path(fold_dir, paste0("fold", fold, ".fa"))
neg_bed   <- file.path(fold_dir, paste0("neg1x_fold", fold, ".bed"))
neg_fa    <- file.path(fold_dir, paste0("neg1x_fold", fold, ".fa"))

stopifnot(file.exists(input_bed))

message("Starting negative-sequence generation")
message("Input BED: ", input_bed)
message("Start time: ", Sys.time())

genNullSeqs(
  input_bed,
  nMaxTrials = 1000,
  xfold = 1,
  genome = BSgenome.Hsapiens.UCSC.hg38.masked,
  outputPosFastaFN = pos_fa,
  outputBedFN = neg_bed,
  outputNegFastaFN = neg_fa
)

message("Finished at: ", Sys.time())
