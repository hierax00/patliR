#' @include AllGenerics.R internal.R network_proximity.R plot_network_layers.R
NULL

## plot_proximity() -- a real constraint discovered while designing this,
## not a bug: network_proximity() only persists the *summary* of its null
## model (d_random_mean, d_random_sd, z_score), not the n_random raw
## resampled distances themselves (see network_proximity.R's @return --
## storing 1000 raw distances per compound/condition in a CSV would bloat
## results/ substantially for a number nobody downstream actually needs
## once the z-score is computed). So this cannot be a literal permutation
## histogram with d_observed overlaid, however common that is in the
## literature -- it is a z-score forest/lollipop plot instead, which is
## the other standard way this kind of result gets shown, and is an
## honest reflection of what is actually stored.

#' Lollipop plot of `network_proximity()` z-scores per compound
#'
#' @description
#' One point per (`condition`, `compound_id`) row of
#' `patliRResults(proj, "network_proximity")`, `z_score` on the y-axis,
#' compounds sorted by `z_score`, with reference lines at `z = 0` and the
#' usual \eqn{z = \pm 1.96} two-sided 95% bands. Points are colored by
#' whether they cross that threshold -- a very negative `z_score` means
#' the compound's targets sit significantly *closer* to the disease module
#' than the degree-matched null model, i.e. topologically meaningful
#' proximity, which is what [network_proximity()] is actually testing for.
#'
#' @section Not a permutation histogram:
#' `network_proximity()` only stores the null model's summary
#' (`d_random_mean`, `d_random_sd`) and the resulting `z_score`, not the
#' `n_random` raw resampled distances themselves -- so there is nothing to
#' draw a literal permutation histogram from. This plot uses what is
#' actually persisted.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character EFO ID, or `NULL` (default) for every disease
#'   present in `patliRResults(proj, "network_proximity")` (faceted).
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "proximity_plot_log")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' abundance <- read.csv(
#'   system.file("extdata", "input_abundance_matrix.csv", package = "patliR"),
#'   check.names = FALSE
#' )
#' proj <- prep_binarize(proj, abundance)
#' proj <- targets_import_batch(
#'   proj,
#'   system.file("extdata", "import_targets", package = "patliR"),
#'   platform = "superpred"
#' )
#' proj <- network_build(proj)
#' proj <- targets_disease_filter(proj, disease = "type 2 diabetes mellitus")
#' proj <- network_proximity(
#'   proj, condition = "FLO-ET",
#'   disease = unique(patliRResults(proj, "targets_disease")$disease_id)[1]
#' )
#' plot_proximity(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_proximity <- function(proj, condition = NULL, disease = NULL,
                            engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                            width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_proximity}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  prox_all <- patliRResults(proj, "network_proximity")
  if (is.null(prox_all) || nrow(prox_all) == 0) {
    cli::cli_abort(c("No {.val network_proximity} entry in {.arg proj}.", "i" = "Run {.fn network_proximity} first."))
  }
  dat <- prox_all[prox_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) dat <- dat[dat$disease_id %in% disease, , drop = FALSE]
  dat <- dat[!is.na(dat$z_score), , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_proximity} rows with a non-NA z_score for the requested condition(s)/disease.")
  }
  dat$label <- .network_layers_labels(proj, conditions, data.frame(name = dat$compound_id, layer = "compound", stringsAsFactors = FALSE))
  dat <- dat[order(dat$z_score), , drop = FALSE]
  dat$label <- factor(dat$label, levels = unique(dat$label))
  dat$significant <- abs(dat$z_score) >= 1.96
  dat$panel <- paste(dat$condition, dat$disease_id, sep = " / ")
  dat$tooltip <- sprintf("%s\nz = %.2f (d_obs = %.2f, null = %.2f +/- %.2f)", dat$label, dat$z_score, dat$d_observed, dat$d_random_mean, dat$d_random_sd)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$label, y = .data$z_score, colour = .data$significant)) +
    ggplot2::geom_hline(yintercept = c(-1.96, 0, 1.96), linetype = c("22", "solid", "22"), colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$z_score), linewidth = 0.4, alpha = 0.6)
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$label), size = 2.5)
  } else {
    ggplot2::geom_point(size = 2.5)
  }
  p <- p +
    ggplot2::coord_flip() +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey40"), name = "|z| >= 1.96") +
    ggplot2::labs(
      title = paste0("Network proximity z-scores -- ", scope_label),
      subtitle = "More negative = targets significantly closer to the disease module than the degree-matched null model",
      x = NULL, y = "z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"))
  if (length(unique(dat$panel)) > 1) p <- p + ggplot2::facet_wrap(~panel, scales = "free_y")

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("proximity_", scope_label, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
    log_row <- data.frame(condition = scope_label, path = path, stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "proximity_plot_log", log_row, "condition")
    patliRResults(proj, "proximity_plot_log") <- log_df
    .write_results_csv(proj, "proximity_plot_log", log_df)
  }

  result <- if (engine == "static") p else ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  if (save) attr(result, "proj") <- proj
  result
}
