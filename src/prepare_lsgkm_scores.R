#!/home/zarazuanav/.conda/envs/ls-gkm/bin/Rscript

library(tidyr)
library(dplyr)

read_scores <- function(file_path, model, task, outer) {
	df <- read.table(file_path,
			 header = FALSE,
			 sep = "\t")

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
		model = model,
		task = task,
		outer = outer,
		seq_length = end - start
	    )
	if (anyNA(df$label)){
		stop("Unexpected value found in the class columns.")
	}

	return(df)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
	stop(
	     "Usage: Rscript src/prepare_lsgkm_scores.R <outer_dir>\n",
	     "Example: Rscript src/prepare_lsgkm_scores.R predictions/evn/ls-gkm/outer1/"
	)
}

dir_path <- normalizePath(
			  args[1],
			  mustWork = FALSE
			  )

pos_file <- file.path(dir_path, "pos_scores.txt")
neg_file <- file.path(dir_path, "neg_scores.txt")

if (!file.exists(pos_file) || !file.exists(neg_file)){
	cat("Error: the directory must contain pos_scores.txt and neg_scores.txt.\n", file = stderr()
	)
	quit(status = 1)
}


outer <- basename(dir_path)
model <- basename(dirname(dir_path))
task <- basename(dirname(dirname(dir_path)))

results_dir <- "results/evn/ls-gkm/"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(results_dir, paste0(outer,".tsv")) 

df_pos <- read_scores(file_path = pos_file, model = model, task = task, outer = outer)
df_neg <- read_scores(file_path = neg_file, model = model, task = task, outer = outer)

df_combined <- bind_rows(df_pos, df_neg)

write.table(df_combined,
	    file=out_file,
	    sep="\t",quote=FALSE,
	    col.names=TRUE,
	    row.names=FALSE
)
