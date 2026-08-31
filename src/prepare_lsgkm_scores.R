#!/home/zarazuanav/.conda/envs/ls-gkm/bin/Rscript

library(tidyr)
library(dplyr)

read_scores <- function(file_path, model, task, dataset, split) {
	df <- read.table(
		file_path,
		header = FALSE,
		sep = "\t"
	)

	df <- df |>
		extract(
			V1,
			into = c("chrom", "start", "end", "class", "index"),
			regex = "^(.*)_([0-9]+)_([0-9]+)_(pos|neg)_([0-9]+)$",
			convert = TRUE
		) |>
		rename(score = V2) |>
		mutate(
			label = case_when(
				class == "pos" ~ 1L,
				class == "neg" ~ 0L,
				TRUE ~ NA_integer_
			),
			split = split,
			model = model,
			task = task,
			dataset = dataset,
			seq_length = end - start
		)

	if (anyNA(df$label)) {
		stop("Unexpected or malformed sequence identifier in ", file_path)
	}

	return(df)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
	stop(
		"Usage: Rscript src/prepare_lsgkm_scores.R <prediction_dir>\n",
		"Examples:\n",
		"  Rscript src/prepare_lsgkm_scores.R predictions/svn/GRCh38/ls-gkm/outer1\n",
		"  Rscript src/prepare_lsgkm_scores.R predictions/svn/GRCh38/ls-gkm/learning_curves/10"
	)
}

dir_path <- normalizePath(args[1], mustWork = TRUE)
run_name <- basename(dir_path)
parent_name <- basename(dirname(dir_path))

if (parent_name == "learning_curves") {
	dataset <- paste0("learning_curve_", run_name)
	output_subdir <- "learning_curves"
} else if (grepl("^outer[1-9][0-9]*$", run_name)) {
	dataset <- run_name
	output_subdir <- ""
} else {
	stop("Unrecognized prediction directory layout: ", dir_path)
}

model <- "lsgkm"
task <- "GRCh38"
results_dir <- file.path("results", "svn", task, "ls-gkm", output_subdir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(results_dir, paste0(run_name, ".tsv"))

score_specs <- expand.grid(
	split = c("validation", "test"),
	class = c("pos", "neg"),
	stringsAsFactors = FALSE
)

score_tables <- lapply(seq_len(nrow(score_specs)), function(i) {
	split <- score_specs$split[i]
	class <- score_specs$class[i]
	file_path <- file.path(dir_path, paste0(split, "_", class, "_scores.txt"))

	if (!file.exists(file_path)) {
		stop("Missing required score file: ", file_path)
	}

	df <- read_scores(
		file_path = file_path,
		model = model,
		task = task,
		dataset = dataset,
		split = split
	)

	if (any(df$class != class)) {
		stop("Sequence class does not match file name in ", file_path)
	}

	return(df)
})

df_combined <- bind_rows(score_tables)

write.table(
	df_combined,
	file = out_file,
	sep = "\t",
	quote = FALSE,
	col.names = TRUE,
	row.names = FALSE
)

cat("Wrote ", nrow(df_combined), " rows to ", out_file, "\n", sep = "")
