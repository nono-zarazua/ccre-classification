chromosome_order <- c(paste0("chr", 1:22), "chrX")

chromosome_axis <- function(data) {
  data |>
    dplyr::distinct(checkpoint_index, endpoint_chromosome) |>
    dplyr::arrange(checkpoint_index)
}

chromosome_x_scale <- function(data) {
  axis <- chromosome_axis(data)
  ggplot2::scale_x_continuous(
    breaks=axis$checkpoint_index,
    labels=sub("^chr", "", axis$endpoint_chromosome),
    expand=ggplot2::expansion(mult=c(.02, .03))
  )
}

chromosome_x_theme <- function() {
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
}

chromosome_composition <- function(data) {
  required <- c("Group", "Split", "checkpoint_index", "endpoint_chromosome",
                "n_training_total", "n_evaluation_total")
  missing <- setdiff(required, names(data))
  if (length(missing))
    stop("Missing chromosome-composition columns: ", paste(missing, collapse=", "))

  training <- data |>
    dplyr::distinct(Group, checkpoint_index, endpoint_chromosome, n_training_total) |>
    dplyr::group_by(Group) |>
    dplyr::arrange(checkpoint_index, .by_group=TRUE) |>
    dplyr::mutate(
      Chromosome=endpoint_chromosome,
      Observations=n_training_total - dplyr::lag(n_training_total, default=0),
      Percentage=Observations / max(n_training_total),
      Set="Training"
    ) |>
    dplyr::ungroup() |>
    dplyr::select(Group, Chromosome, Observations, Percentage, Set)

  evaluation <- data |>
    dplyr::mutate(
      Split=as.character(Split),
      Chromosome=dplyr::case_when(
        Split %in% c("validation", "Validation") ~ "chr1",
        Split %in% c("test", "Test") ~ "chr2",
        TRUE ~ NA_character_
      ),
      Set=dplyr::case_when(
        Chromosome == "chr1" ~ "Validation",
        Chromosome == "chr2" ~ "Test",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(Chromosome)) |>
    dplyr::distinct(Group, Chromosome, n_evaluation_total, Set) |>
    dplyr::transmute(Group, Chromosome, Observations=n_evaluation_total,
                     Percentage=NA_real_, Set)

  dplyr::bind_rows(evaluation, training) |>
    dplyr::mutate(
      Chromosome=factor(Chromosome, levels=rev(chromosome_order)),
      Label=dplyr::if_else(
        Set == "Training",
        paste0(scales::comma(Observations), " (", scales::percent(Percentage, accuracy=.1), ")"),
        scales::comma(Observations)
      )
    ) |>
    dplyr::arrange(Group, Chromosome)
}

chromosome_composition_plot <- function(data, title="Observations per chromosome",
                                        palette=NULL) {
  p <- ggplot2::ggplot(data, ggplot2::aes(Chromosome, Observations, fill=Group)) +
    ggplot2::geom_col(width=.72) +
    ggplot2::geom_text(ggplot2::aes(label=Label), hjust=-.08, size=3) +
    ggplot2::facet_wrap(~Group, ncol=1, scales="free_x") +
    ggplot2::coord_flip(clip="off") +
    ggplot2::scale_y_continuous(labels=scales::comma,
                                expand=ggplot2::expansion(mult=c(0, .28))) +
    ggplot2::labs(
      title=title,
      subtitle="Parentheses show each chr3–chrX contribution to the corresponding training pool",
      x=NULL, y="Observations", fill=NULL
    ) +
    ggplot2::guides(fill="none") +
    ggplot2::theme(plot.margin=ggplot2::margin(5.5, 80, 5.5, 5.5))
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_manual(values=palette)
  p
}
