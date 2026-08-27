#' @include AllGenerics.R internal.R network_module_robustness.R
NULL

#' Percolation curve(s) from `network_module_robustness()`
#'
#' @description
#' Line plot of `largest_component_fraction` vs. `n_removed` (the raw
#' percolation curve behind each module's `r_index`, from
#' `network_robustness_curve`), one line per module, faceted by
#' `condition`. `r_index` (the curve's area, from `network_module_robustness`)
#' is annotated on each panel -- a module that stays near 1.0 for most of
#' the curve (robust to node removal) has `r_index` close to 1; one that
#' collapses after a few removals (fragile, likely hub-dependent) has
#' `r_index` close to 0 (Schneider et al. 2011's R-index, cited in
#' [network_module_robustness()]'s own roxygen).
#'
#' @inheritParams network_build
#' @param module_id Character vector of module IDs to restrict to, or
#'   `NULL` (default) for every module in scope.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "robustness_plot_log")`.
#'
#' @examples
#' \donttest{
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
#' proj <- network_module_robustness(proj, condition = "FLO-ET")
#' plot_robustness(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_robustness <- function(proj, condition = NULL, module_id = NULL,
                             engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                             width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_robustness}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  curve_all <- patliRResults(proj, "network_robustness_curve")
  summary_all <- patliRResults(proj, "network_module_robustness")
  if (is.null(curve_all) || nrow(curve_all) == 0) {
    cli::cli_abort(c("No {.val network_robustness_curve} entry in {.arg proj}.", "i" = "Run {.fn network_module_robustness} first."))
  }
  curve <- curve_all[curve_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(module_id)) curve <- curve[curve$module_id %in% module_id, , drop = FALSE]
  if (nrow(curve) == 0) {
    cli::cli_abort("No {.val network_robustness_curve} rows for condition(s) {.val {conditions}}{if (!is.null(module_id)) paste0(' / module_id ', toString(module_id)) else ''}.")
  }
  curve$panel <- paste(curve$condition, curve$module_id, sep = " / ")

  ann <- NULL
  if (!is.null(summary_all) && nrow(summary_all) > 0) {
    ann <- summary_all[summary_all$condition %in% conditions & (is.null(module_id) | summary_all$module_id %in% module_id), , drop = FALSE]
    ann$panel <- paste(ann$condition, ann$module_id, sep = " / ")
    ann$label <- sprintf("R = %.3f", ann$r_index)
  }

  p <- ggplot2::ggplot(curve, ggplot2::aes(x = .data$n_removed, y = .data$largest_component_fraction))
  p <- if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p + ggiraph::geom_line_interactive(ggplot2::aes(tooltip = .data$panel, data_id = .data$panel), colour = "#2980b9", linewidth = 0.8)
  } else {
    p + ggplot2::geom_line(colour = "#2980b9", linewidth = 0.8)
  }
  p <- p + ggplot2::geom_point(size = 1, colour = "#2980b9") + ggplot2::facet_wrap(~panel, scales = "free_x")
  if (!is.null(ann) && nrow(ann) > 0) {
    p <- p + ggplot2::geom_text(
      data = ann, ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      hjust = 1.1, vjust = 1.5, size = 3.2, colour = "grey30", inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::labs(
      title = paste0("Module percolation robustness -- ", scope_label),
      subtitle = "Fraction of the module remaining in the largest component as nodes are removed one at a time; R = area under the curve (Schneider et al. 2011)",
      x = "Nodes removed", y = "Largest component fraction"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"))

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("robustness_", scope_label, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
    log_row <- data.frame(condition = scope_label, path = path, stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "robustness_plot_log", log_row, "condition")
    patliRResults(proj, "robustness_plot_log") <- log_df
    .write_results_csv(proj, "robustness_plot_log", log_df)
  }

  result <- if (engine == "static") p else ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  if (save) attr(result, "proj") <- proj
  result
}
