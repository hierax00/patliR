#' @include AllGenerics.R internal.R plot-helpers.R network_condition_compare.R
NULL

#' Trend line: how the extract's total effect changes across an ordered
#' sequence of conditions
#'
#' @description
#' [plot_condition_compare()] ranks conditions best-to-worst; this instead
#' walks them **in a caller-given sequence** (a phenological order, a
#' treatment order, ...) and draws `aggregate`'s z-score as a line, so a
#' *trend* across that sequence is visible -- does the effect grow, shrink,
#' or hold steady from one condition to the next -- optionally split into
#' several lines (`group`, e.g. greenhouse vs. wild) and/or one panel per
#' disease. `n_compounds_present` is drawn as a second, lighter line on its
#' own axis alongside it, so a change in the effect can be read against
#' whether it tracks a change in how many compounds were even present.
#'
#' @inheritParams network_build
#' @param conditions Character vector, **the sequence's left-to-right
#'   order** -- not resolved or sorted for you. Must already have a row in
#'   `patliRResults(proj, "network_condition_compare")` for `disease`
#'   (run [network_condition_compare()] first; this function does not
#'   compute anything).
#' @param disease Character vector of one or more disease ids, or `NULL`
#'   (default) for every one in `patliRResults(proj,
#'   "network_condition_compare")`. More than one draws a panel per disease.
#' @param group `NULL` (default), or a named character vector mapping each
#'   condition in `conditions` to a group (e.g. `c("EVEG-I" = "Invernadero",
#'   "EVEG-S" = "Silvestre")`) -- drawn as separate coloured lines sharing
#'   the same x sequence, for comparing e.g. two origins across the same
#'   developmental stages.
#' @param condition_label `NULL` (default), or a named character vector
#'   mapping each condition in `conditions` to a **display label on x**,
#'   independent of the actual join key -- e.g. `c("EVEG-I" = "Vegetativa",
#'   "EVEG-S" = "Vegetativa", ...)` so two conditions from different
#'   `group`s that represent the same developmental stage share one x
#'   position (and their lines are visibly aligned) instead of each getting
#'   its own tick along an interleaved sequence. Conditions that map to the
#'   same label must be adjacent (or interleaved consistently) in
#'   `conditions`; x order follows each label's first occurrence.
#' @param aggregate `"mean"` (default) or `"median"` -- must match a
#'   [network_condition_compare()] call already stored in `proj`.
#' @inheritParams plot_save_params
#'
#' @return A `ggplot`. `attr(., "table")` holds the exact rows plotted. If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "condition_trend_plot_log")`.
#'
#' @seealso [network_condition_compare()], [plot_condition_compare()], [plot_condition_flow()]
#' @export
plot_condition_trend <- function(proj, conditions, disease = NULL, group = NULL, condition_label = NULL,
                                 aggregate = c("mean", "median"),
                                 save = TRUE, out_dir = NULL, width = 8, height = 5, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  aggregate <- match.arg(aggregate)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_condition_trend}.")
  }
  if (!is.character(conditions) || length(conditions) < 2 || anyNA(conditions)) {
    cli::cli_abort("{.arg conditions} must be a character vector of at least 2 condition names, in the order to draw them.")
  }
  if (!is.null(group)) {
    if (is.null(names(group)) || !all(conditions %in% names(group))) {
      cli::cli_abort("{.arg group} must be a character vector named by every value in {.arg conditions}.")
    }
  }
  if (!is.null(condition_label)) {
    if (is.null(names(condition_label)) || !all(conditions %in% names(condition_label))) {
      cli::cli_abort("{.arg condition_label} must be a character vector named by every value in {.arg conditions}.")
    }
  }

  cmp <- patliRResults(proj, "network_condition_compare")
  if (is.null(cmp) || nrow(cmp) == 0) {
    cli::cli_abort(c("No {.val network_condition_compare} entry in {.arg proj}.",
                     "i" = "Run {.fn network_condition_compare} first."))
  }
  if (is.null(disease)) disease <- unique(cmp$disease_id)
  dat <- cmp[cmp$condition %in% conditions & cmp$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat) == 0) cli::cli_abort("No rows for {.arg conditions} x {.val {disease}} in {.val network_condition_compare}.")
  expected <- length(conditions) * length(disease)
  if (nrow(dat) < expected) {
    cli::cli_warn("Only {nrow(dat)} of {expected} (condition, disease) combinations were found; the missing ones are left out of the plot.")
  }

  z_col <- paste0(aggregate, "_z")
  x_lab <- if (is.null(condition_label)) conditions else unname(condition_label[conditions])
  dat$x <- factor(if (is.null(condition_label)) dat$condition else unname(condition_label[as.character(dat$condition)]),
                  levels = unique(x_lab))
  dat$condition <- factor(dat$condition, levels = conditions)
  dat$disease_label <- vapply(dat$disease_id, function(d) .plot_disease_label(proj, d, with_id = FALSE), character(1))
  dat$group <- if (is.null(group)) "" else unname(group[as.character(dat$condition)])
  multi_disease <- length(unique(dat$disease_id)) > 1
  multi_group <- !is.null(group) && length(unique(dat$group)) > 1

  ## n_compounds_present rescaled onto the z-score's own range, so it can share the panel
  ## as a second, visually distinct line without a second y-axis's usual misreading risk
  rng_z <- range(dat[[z_col]], na.rm = TRUE)
  rng_n <- range(dat$n_compounds_present, na.rm = TRUE)
  scale_n <- function(n) {
    if (diff(rng_n) == 0) return(rep(mean(rng_z), length(n)))
    rng_z[1] + (n - rng_n[1]) / diff(rng_n) * diff(rng_z)
  }
  dat$n_scaled <- scale_n(dat$n_compounds_present)

  colour_var <- if (multi_group) "group" else if (multi_disease) "disease_label" else NULL
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data[[z_col]]))
  if (!is.null(colour_var)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = .data[[colour_var]], group = .data[[colour_var]]), linewidth = 0.9) +
      ggplot2::geom_point(ggplot2::aes(colour = .data[[colour_var]]), size = 2.2) +
      ggplot2::geom_line(ggplot2::aes(y = .data$n_scaled, colour = .data[[colour_var]], group = .data[[colour_var]]),
                         linetype = "dotted", linewidth = 0.6, alpha = 0.7) +
      ggplot2::scale_colour_manual(values = .plot_contrast_palette(length(unique(dat[[colour_var]]))) |>
                                     stats::setNames(sort(unique(dat[[colour_var]]))),
                                   name = if (multi_group) "Group" else "Disease")
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = 1), colour = "#2980b9", linewidth = 0.9) +
      ggplot2::geom_point(colour = "#2980b9", size = 2.2) +
      ggplot2::geom_line(ggplot2::aes(y = .data$n_scaled, group = 1), colour = "#2980b9",
                         linetype = "dotted", linewidth = 0.6, alpha = 0.7)
  }
  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::scale_y_continuous(
      name = paste0(aggregate, " proximity z-score (solid)"),
      sec.axis = ggplot2::sec_axis(
        function(z) rng_n[1] + (z - rng_z[1]) / diff(rng_z) * diff(rng_n),
        name = "Compounds present (dotted)"
      )
    ) +
    ggplot2::labs(
      title = "How the extract's total effect changes across conditions",
      subtitle = .plot_wrap(paste0(
        "Solid line = ", aggregate, " proximity z-score (more negative = stronger effect); ",
        "dotted line = number of compounds present in that condition, on its own secondary axis."
      ), .plot_wrap_width(width, 8)),
      x = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                  plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
                  plot.title.position = "plot")
  if (multi_disease) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$disease_label))

  result <- .plot_finish(
    proj, p,
    name = "condition_trend_plot_log",
    filename = paste0("condition_trend_", paste(conditions, collapse = "-"), "_",
                      paste(sort(unique(dat$disease_id)), collapse = "+"), "_", aggregate, ".png"),
    log_row = data.frame(conditions = paste(conditions, collapse = ">"),
                         disease_id = paste(sort(unique(dat$disease_id)), collapse = "+"),
                         aggregate = aggregate, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("conditions", "disease_id", "aggregate"),
    engine = "static", save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- dat[, c("condition", "disease_id", "group", z_col, "n_compounds_present")]
  result
}
