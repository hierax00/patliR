#' @include AllGenerics.R internal.R network_synergy.R plot-helpers.R
NULL

#' Volcano-style scatter of `network_synergy()` compound pairs
#'
#' @description
#' One point per compound pair from `patliRResults(proj, "network_synergy")`:
#' `complementarity` (targets that differ between the pair -- 0 = same
#' targets, 1 = fully distinct) on the x-axis, `joint_closeness` (how close
#' *both* compounds sit to the disease module) on the y-axis, point size/
#' colour by `synergy_score`. High `synergy_score` needs both distinct
#' targets *and* both close to the disease module -- not just one of the
#' two -- so the interesting pairs are the ones in the upper-right, not
#' just whichever axis is largest alone (see [network_synergy()]'s own
#' roxygen for why `synergy_score` is built that way). The `top_n` highest-
#' scoring pairs are labeled by compound name.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character EFO ID, or `NULL` (default) for every disease
#'   present (faceted).
#' @param top_n Integer, default `10`. Number of highest-`synergy_score`
#'   pairs to label.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "synergy_plot_log")`.
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
#' proj <- network_synergy(proj, condition = "FLO-ET")
#' plot_synergy(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_synergy <- function(proj, condition = NULL, disease = NULL, top_n = 10,
                          engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                          width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 0)
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  syn_all <- patliRResults(proj, "network_synergy")
  if (is.null(syn_all) || nrow(syn_all) == 0) {
    cli::cli_abort(c("No {.val network_synergy} entry in {.arg proj}.", "i" = "Run {.fn network_synergy} first."))
  }
  dat <- syn_all[syn_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) dat <- dat[dat$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_synergy} rows for the requested condition(s)/disease.")
  }
  ## `synergy_score` is `NA` for every pair that is not Complementary
  ## Exposure (P2) -- or, in `separation = "network"` mode, whenever the
  ## proximity BH gate was arithmetically unreachable. An all-NA score
  ## column would otherwise make every point vanish with a bare ggplot
  ## warning, so draw a self-explanatory empty panel instead.
  n_scored <- sum(!is.na(dat$synergy_score))
  if (n_scored == 0) {
    cheng <- if ("cheng_class" %in% names(dat)) {
      tb <- table(dat$cheng_class, useNA = "ifany")
      paste0(names(tb), " x", as.integer(tb), collapse = ", ")
    } else NA_character_
    msg <- paste0(
      "No scored compound pairs for ", scope_label, ".\n\n",
      "synergy_score is only set for Complementary-Exposure (P2) pairs",
      if (!is.na(cheng)) paste0(".\nCheng classes present: ", cheng) else ".",
      "\n\nIf network_synergy() warned that the BH-adjusted proximity p\n",
      "cannot clear alpha, rerun network_proximity() with a larger n_random."
    )
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3.4,
                        colour = "grey25", lineheight = 1.1) +
      ggplot2::labs(title = paste0("Compound pair synergy -- ", scope_label)) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
    return(.plot_finish(
      proj, p,
      name = "synergy_plot_log",
      filename = paste0("synergy_", scope_label, ".png"),
      log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
      key_cols = "condition",
      engine = engine, save = save, out_dir = out_dir,
      width = width, height = height, dpi = dpi
    ))
  }
  if (n_scored < nrow(dat)) {
    cli::cli_inform(c(
      "i" = "{nrow(dat) - n_scored} of {nrow(dat)} pair{?s} {?is/are} unscored ({.field synergy_score} {.val NA}) and {?is/are} not drawn -- only Complementary-Exposure (P2) pairs get a score."
    ))
  }

  label_a <- .plot_label_nodes(proj, conditions, dat$compound_a, "compound")
  label_b <- .plot_label_nodes(proj, conditions, dat$compound_b, "compound")
  dat$pair_label <- paste(label_a, "+", label_b)
  dat$panel <- paste(dat$condition, dat$disease_id, sep = " / ")
  dat$tooltip <- sprintf("%s\ncomplementarity %.2f, joint closeness %.2f, synergy %.2f", dat$pair_label, dat$complementarity, dat$joint_closeness, dat$synergy_score)

  top <- utils::head(dat[order(-dat$synergy_score), , drop = FALSE], top_n)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$complementarity, y = .data$joint_closeness, size = .data$synergy_score, colour = .data$synergy_score))
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$pair_label), alpha = 0.75)
  } else {
    ggplot2::geom_point(alpha = 0.75)
  }
  if (nrow(top) > 0) {
    p <- p + ggplot2::geom_text(data = top, ggplot2::aes(label = .data$pair_label), size = 2.6, colour = "grey15", vjust = -1, show.legend = FALSE)
  }
  p <- p +
    ggplot2::scale_colour_gradient(low = "grey70", high = "#c0392b", name = "Synergy\nscore") +
    ggplot2::scale_size(range = c(1.5, 6), guide = "none") +
    ggplot2::labs(
      title = paste0("Compound pair synergy -- ", scope_label),
      subtitle = "Top-right = distinct targets (complementarity) AND both close to the disease module (joint_closeness); top pairs labeled",
      x = "Complementarity (target Jaccard distance)", y = "Joint closeness to disease module"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"))
  if (length(unique(dat$panel)) > 1) p <- p + ggplot2::facet_wrap(~panel)

  .plot_finish(
    proj, p,
    name = "synergy_plot_log",
    filename = paste0("synergy_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}
