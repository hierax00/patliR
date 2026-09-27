#' @include AllGenerics.R internal.R plot-helpers.R network_condition_compare.R
NULL

#' Plot `network_condition_compare()` results: which extract has the best effect
#'
#' @description
#' One point per `(condition, disease)`, `aggregate`'s z-score (`mean_z` or
#' `median_z`) on the x-axis, conditions sorted with the most negative
#' (closest to the disease module, averaged over the disease(s) drawn) at
#' the top -- the same reading as [plot_proximity()]'s `view = "z"`, one row
#' per condition instead of per compound. A range bar shows one SD (`sd_z`)
#' when at least two compounds contributed a finite z-score;
#' `n_compounds_present` is annotated beside each point.
#'
#' @section More than one disease:
#' Pass `disease` as a vector (or leave it `NULL` when
#' `network_condition_compare()` was run for more than one disease) to draw
#' every condition's z-score for **each** disease side by side in the same
#' panel, colour = disease -- this is how to answer "which extract has the
#' strongest effect, across diseases" in one figure, instead of one figure
#' per disease. Conditions are still a single ranking (by the mean of the
#' plotted diseases' z, so the ordering does not depend on which disease
#' happens to be listed first).
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character disease id, a vector of several, or `NULL`
#'   (default) for every `disease_id` present in `patliRResults(proj,
#'   "network_condition_compare")`.
#' @param aggregate `"mean"` (default) or `"median"` -- must match a
#'   [network_condition_compare()] call already stored in `proj`.
#'
#' @return A `ggplot`. `attr(., "table")` holds the exact rows plotted
#'   (`patliRResults(proj, "network_condition_compare")`, filtered to
#'   `disease` and sorted). If `save = TRUE` (default), also writes a PNG
#'   and logs it to `patliRResults(proj, "condition_compare_plot_log")`.
#'
#' @seealso [network_condition_compare()], [plot_proximity()], [plot_condition_disease_flow()]
#' @export
plot_condition_compare <- function(proj, disease = NULL, aggregate = c("mean", "median"),
                                   save = TRUE, out_dir = NULL, width = 7, height = NULL, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  aggregate <- match.arg(aggregate)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_condition_compare}.")
  }

  cmp <- patliRResults(proj, "network_condition_compare")
  if (is.null(cmp) || nrow(cmp) == 0) {
    cli::cli_abort(c("No {.val network_condition_compare} entry in {.arg proj}.",
                     "i" = "Run {.fn network_condition_compare} first."))
  }
  if (is.null(disease)) disease <- unique(cmp$disease_id)
  dat <- cmp[cmp$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat) == 0) cli::cli_abort("No rows for {.val {disease}} in {.val network_condition_compare}.")
  multi <- length(unique(dat$disease_id)) > 1

  z_col <- paste0(aggregate, "_z")
  overall <- stats::aggregate(dat[[z_col]], by = list(condition = dat$condition), FUN = mean, na.rm = TRUE)
  ord <- overall$condition[order(overall$x)]                     # most negative (best) first
  dat$condition <- factor(dat$condition, levels = rev(ord))       # ... drawn at the top
  dat$ymin <- dat[[z_col]] - ifelse(is.na(dat$sd_z), 0, dat$sd_z)
  dat$ymax <- dat[[z_col]] + ifelse(is.na(dat$sd_z), 0, dat$sd_z)
  ## disease_label is display text only -- two distinct disease_ids can
  ## share a label (e.g. two imported custom modules named the same thing),
  ## and grouping/colouring by the label rather than the id would silently
  ## fuse them into one dodge position and one legend entry.
  dat$disease_label <- vapply(dat$disease_id, function(d) .plot_disease_label(proj, d, with_id = FALSE), character(1))
  disease_label <- paste(unique(dat$disease_label), collapse = " / ")

  if (is.null(height)) height <- max(3, 0.5 * length(ord) + 1.5)
  dodge <- if (multi) ggplot2::position_dodge(width = 0.5) else ggplot2::position_identity()

  p <- ggplot2::ggplot(dat, ggplot2::aes(y = .data$condition, x = .data[[z_col]],
                                         colour = .data$disease_id, group = .data$disease_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$ymin, xmax = .data$ymax), height = 0.15,
                            position = dodge, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_compounds_present), position = dodge) +
    ggplot2::scale_size_continuous(name = "Compounds\npresent", range = c(2, 6))
  if (multi) {
    disease_id_labels <- stats::setNames(dat$disease_label, dat$disease_id)
    disease_id_labels <- disease_id_labels[!duplicated(names(disease_id_labels))]
    p <- p + ggplot2::scale_colour_brewer(palette = "Set1", name = "Disease", labels = disease_id_labels)
  } else {
    p <- p + ggplot2::scale_colour_manual(values = "#2980b9", guide = "none") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("n=%d", .data$n_compounds_present)),
                        hjust = -0.35, size = 3, colour = "grey30")
  }
  p <- p +
    ggplot2::labs(
      title = .plot_wrap(paste0("Extracts ranked by mean proximity z -- ", disease_label), .plot_wrap_width(width, 12)),
      subtitle = .plot_wrap(paste0(
        "More negative ", aggregate, " z = that condition's compounds sit, on average, closer ",
        "to the disease module than the degree-matched null (no significance test at the condition level; ",
        "error bars = 1 SD across its compounds, not an uncertainty of the mean).",
        if (multi) " Conditions are ranked by their mean z across the disease(s) shown." else ""
      ), .plot_wrap_width(width, 8)),
      x = paste0(aggregate, " proximity z-score (per compound, this condition)"), y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                  plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
                  plot.title.position = "plot")

  result <- .plot_finish(
    proj, p,
    name = "condition_compare_plot_log",
    filename = paste0("condition_compare_", paste(sort(unique(dat$disease_id)), collapse = "+"), "_", aggregate, ".png"),
    log_row = data.frame(disease_id = paste(sort(unique(dat$disease_id)), collapse = "+"), aggregate = aggregate,
                         n_conditions = length(ord), path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("disease_id", "aggregate"),
    engine = "static", save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- dat
  result
}
