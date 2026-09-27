#' @include AllGenerics.R internal.R plot-helpers.R network_condition_compare.R
NULL

## Uses ggalluvial (Suggests), same rationale as plot_bowtie.R/plot_gpcr_flow.R.

#' Alluvial flow: condition (stage) -> disease, by proximity effect
#'
#' @description
#' Every condition in `patliRResults(proj, "network_condition_compare")`
#' flowing into every disease it was compared against, flow thickness = how
#' strong that condition's effect on that disease is (`-mean_z`/`-median_z`
#' -- proximity z is negative when the extract sits closer to the disease
#' module than the null, so its sign is flipped here for "thicker = a
#' stronger effect"). Reads like a small metro map of "which disease does
#' each phenological stage/extract line up with best" -- the [plot_gpcr_flow()]
#' idea applied one level up, at the whole-extract level instead of
#' compound-to-receptor. Every condition and disease already in
#' [network_condition_compare()]'s table is drawn (no `top_n`): this is
#' meant to be read as one complete picture across all of them, not a
#' shortlist.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character vector of disease ids, or `NULL` (default) for
#'   every one in `patliRResults(proj, "network_condition_compare")`.
#' @param aggregate `"mean"` (default) or `"median"` -- must match a
#'   [network_condition_compare()] call already stored in `proj`.
#'
#' @return A `ggplot` (`ggalluvial` has no `ggiraph` equivalent). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "condition_disease_flow_plot_log")`.
#'   `attr(., "table")` holds the exact `(condition, disease_id, z, weight)`
#'   rows drawn.
#'
#' @seealso [network_condition_compare()], [plot_condition_compare()], [plot_gpcr_flow()]
#' @export
plot_condition_disease_flow <- function(proj, disease = NULL, aggregate = c("mean", "median"),
                                        save = TRUE, out_dir = NULL, width = 8, height = NULL, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  aggregate <- match.arg(aggregate)
  .plot_require(extra = "ggalluvial")

  cmp <- patliRResults(proj, "network_condition_compare")
  if (is.null(cmp) || nrow(cmp) == 0) {
    cli::cli_abort(c("No {.val network_condition_compare} entry in {.arg proj}.",
                     "i" = "Run {.fn network_condition_compare} first."))
  }
  if (is.null(disease)) disease <- unique(cmp$disease_id)
  dat <- cmp[cmp$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat) == 0) cli::cli_abort("No rows for {.val {disease}} in {.val network_condition_compare}.")
  if (length(unique(dat$disease_id)) < 2) {
    cli::cli_abort(c("{.fn plot_condition_disease_flow} needs at least two diseases to draw a flow.",
                     "i" = "Only {.val {unique(dat$disease_id)}} is available; use {.fn plot_condition_compare} for a single disease."))
  }

  z_col <- paste0(aggregate, "_z")
  ## disease_label is display text only -- two distinct disease_ids sharing
  ## a label would otherwise merge into a single alluvial stratum (axis2
  ## groups by whatever value it is given), silently combining two different
  ## diseases' flows into one node.
  dat$disease_label <- vapply(dat$disease_id, function(d) .plot_disease_label(proj, d, with_id = FALSE), character(1))
  ## flow thickness = strength of the effect (more negative z -> thicker); a z at or
  ## above 0 (no signal, or weaker than the null) still gets a thin, floored flow
  ## rather than vanishing or breaking the alluvial's positive-weight requirement
  dat$weight <- pmax(-dat[[z_col]], 0.05)
  ## TRUE for the disease each condition lines up with most strongly (its most negative z)
  dat$best_for_condition <- ave(dat[[z_col]], dat$condition, FUN = function(z) z == min(z)) == 1

  ## conditions ordered by their strongest (most negative) z across the diseases drawn
  cond_rank <- stats::aggregate(dat[[z_col]], by = list(condition = dat$condition), FUN = min)
  cond_order <- cond_rank$condition[order(cond_rank$x)]
  dat$condition <- factor(dat$condition, levels = cond_order)

  if (is.null(height)) height <- max(4, 0.4 * length(cond_order) + 2)
  disease_ids_present <- sort(unique(dat$disease_id))
  disease_colors <- .plot_contrast_palette(length(disease_ids_present))
  names(disease_colors) <- disease_ids_present
  disease_id_labels <- stats::setNames(dat$disease_label, dat$disease_id)
  disease_id_labels <- disease_id_labels[!duplicated(names(disease_id_labels))]
  ## the stratum text is drawn for BOTH axes (condition names as well as
  ## disease ids) -- a condition maps to itself, a disease_id to its label
  stratum_labels <- c(stats::setNames(as.character(cond_order), as.character(cond_order)), disease_id_labels)

  p <- ggplot2::ggplot(dat, ggplot2::aes(axis1 = .data$condition, axis2 = .data$disease_id, y = .data$weight)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$disease_id), width = 1 / 5, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1 / 5, fill = "grey92", colour = "grey40") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(
                         .plot_truncate(unname(stratum_labels[as.character(stratum)]), 22)
                       )), size = 2.9) +
    ggplot2::scale_x_discrete(limits = c("condition", "disease"), expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::scale_fill_manual(values = disease_colors, name = "Disease", labels = disease_id_labels) +
    ggplot2::labs(
      title = "Which disease does each extract line up with best?",
      subtitle = .plot_wrap(paste0(
        "Flow thickness = strength of the ", aggregate, " proximity effect (-", aggregate, "_z; thicker = closer ",
        "to that disease's module than the degree-matched null). Every condition and disease compared is drawn."
      ), .plot_wrap_width(width, 8)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_blank()
    )

  result <- .plot_finish(
    proj, p,
    name = "condition_disease_flow_plot_log",
    filename = paste0("condition_disease_flow_", paste(sort(unique(dat$disease_id)), collapse = "+"), "_", aggregate, ".png"),
    log_row = data.frame(disease_id = paste(sort(unique(dat$disease_id)), collapse = "+"), aggregate = aggregate,
                         n_conditions = length(cond_order), path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("disease_id", "aggregate"),
    engine = NULL, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- dat[, c("condition", "disease_id", z_col, "weight", "best_for_condition")]
  result
}
