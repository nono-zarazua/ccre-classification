#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gkmSVM)
  library(BSgenome.Hsapiens.UCSC.hg38.masked)
})

input_bed <- "data/processed/hepg2_ELS.trimmed.bed"
pos_fa    <- "data/processed/hepg2_ELS.trimmed.fa"
neg_bed   <- "data/processed/negatives/neg1x_hepg2_ELS.trimmed.bed"
neg_fa    <- "data/processed/negatives/neg1x_hepg2_ELS.trimmed.fa"

stopifnot(file.exists(input_bed))

message("Starting negative-sequence generation")
message("Input BED: ", input_bed)
message("Start time: ", Sys.time())

genNullSeqs(
  input_bed,
  nMaxTrials = 10,
  xfold = 1,
  genome = BSgenome.Hsapiens.UCSC.hg38.masked,
  outputPosFastaFN = pos_fa,
  outputBedFN = neg_bed,
  outputNegFastaFN = neg_fa
)

message("Finished at: ", Sys.time())
