#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(gkmSVM)
    library(BSgenome.Hsapiens.UCSC.hg38.masked)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop(
        "Usage: Rscript src/build_negs.R <input_bed>\n",
        "Example: Rscript src/build_negs.R ",
        "data/processed/evn/GRCh38/chroms/chrom1/chrom1_GRCh38_ELS.bed"
    )
}

input_bed <- args[1]

if (!file.exists(input_bed)) {
    stop("Input BED does not exist: ", input_bed)
}

input_dir <- dirname(input_bed)
input_stem <- tools::file_path_sans_ext(basename(input_bed))

pos_fa <- file.path(input_dir, paste0(input_stem, ".fa"))
neg_prefix <- file.path(input_dir, paste0("neg1x_", input_stem))
neg_bed <- paste0(neg_prefix, ".bed")
neg_fa <- paste0(neg_prefix, ".fa")

message("Starting negative-sequence generation")
message("Input BED: ", input_bed)
message("Positive FASTA: ", pos_fa)
message("Negative BED: ", neg_bed)
message("Negative FASTA: ", neg_fa)
message("Start time: ", Sys.time())

genNullSeqs(
    input_bed,
    nMaxTrials = 500,
    xfold = 1,
    genome = BSgenome.Hsapiens.UCSC.hg38.masked,
    outputPosFastaFN = pos_fa,
    outputBedFN = neg_bed,
    outputNegFastaFN = neg_fa
)

message("Finished at: ", Sys.time())
