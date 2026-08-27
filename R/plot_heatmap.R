#' @include AllGenerics.R internal.R network_build.R plot_network_layers.R
NULL

## Uses pheatmap::pheatmap() for real hierarchical clustering on both axes.
## Contract differs from the rest of plot_* on purpose: pheatmap is
## grid/base graphics, not ggplot2, so there is no `engine` argument
## (ggiraph has nothing to wrap), and `save = TRUE` writes PDF (via
## pheatmap's own `filename` argument) rather than PNG -- a clustered
## heatmap with printed cell values reads better as vector output.

#' Heatmap of compound-target association or compound-condition presence
#' (real `pheatmap()`)
#'
#' @description
#' `what = "compound_target"` (default): compound x target matrix, filled
#' by `network_edges$weight` (the import probability), for one or more
#' conditions pooled. `what = "compound_condition"`: compound x condition
#' matrix straight from [binarizedMatrix()] (presence/absence, 0/1) --
#' which compounds were actually detected in which extract. Both are
#' rendered with `pheatmap::pheatmap()` -- hierarchical clustering on both
#' axes (dropped automatically, with a warning, on whichever axis has
#' fewer than 2 rows/columns -- `hclust()` cannot cluster a single item),
#' cell values printed on top of the color scale.
#'
#' @inheritParams network_build
#' @param what `"compound_target"` (default) or `"compound_condition"`.
#' @param save Logical, default `TRUE`. If `TRUE`, writes a PDF via
#'   `pheatmap::pheatmap(..., filename = ...)` and logs the path to
#'   `patliRResults(proj, "heatmap_plot_log")`. If `FALSE`, draws to the
#'   current graphics device instead (same as calling `pheatmap()`
#'   directly without `filename`).
#' @param out_dir Directory for the PDF, default `file.path(projectDir(proj), "plots")`.
#' @param width,height `NULL` (default) to size the plot to the matrix
#'   (`max(8, n_cols * 0.9 + 2)` / `max(6, n_rows * 0.5 + 2)`), or an
#'   explicit number of inches.
#'
#' @return The `pheatmap` object returned by `pheatmap::pheatmap()`
#'   (`list(tree_row, tree_col, kmeans, gtable)`) -- re-drawable with
#'   `grid::grid.draw(result$gtable)`. If `save = TRUE`, also carries the
#'   updated `proj` as `attr(result, "proj")` (same pattern as every other
#'   `plot_*` function).
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
#' plot_heatmap(proj, what = "compound_condition", save = FALSE)
#' plot_heatmap(proj, condition = "FLO-ET", what = "compound_target", save = FALSE)
#' }
#'
#' @export
plot_heatmap <- function(proj, condition = NULL, what = c("compound_target", "compound_condition"),
                          save = TRUE, out_dir = NULL, width = NULL, height = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  what <- match.arg(what)
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    cli::cli_abort("The {.pkg pheatmap} package is required for {.fn plot_heatmap}.")
  }

  if (what == "compound_condition") {
    bin <- binarizedMatrix(proj)
    if (nrow(bin) == 0) {
      cli::cli_abort(c("No binarized abundance matrix in {.arg proj}.", "i" = "Run {.fn prep_binarize} first."))
    }
    conds <- setdiff(names(bin), "compound_id")
    mat <- as.matrix(bin[, conds, drop = FALSE])
    rownames(mat) <- .network_layers_labels(proj, conds, data.frame(name = bin$compound_id, layer = "compound", stringsAsFactors = FALSE))
    storage.mode(mat) <- "double"
    scope_label <- "ALL"
    title <- "Compound x Condition Presence"
  } else {
    conditions <- .network_resolve_conditions(proj, condition)
    scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")
    edges_all <- patliRResults(proj, "network_edges")
    long <- edges_all[edges_all$condition %in% conditions, c("compound_id", "uniprot_id", "weight")]
    long <- stats::aggregate(weight ~ compound_id + uniprot_id, long, max) # pooled scope: keep the strongest hit per pair
    if (nrow(long) == 0) {
      cli::cli_abort("Nothing to plot for {.arg what} = {.val {what}} in this scope.")
    }
    long$compound_label <- .network_layers_labels(proj, conditions, data.frame(name = long$compound_id, layer = "compound", stringsAsFactors = FALSE))
    long$target_label <- .network_layers_labels(proj, conditions, data.frame(name = long$uniprot_id, layer = "target", stringsAsFactors = FALSE))
    mat <- as.matrix(stats::xtabs(weight ~ compound_label + target_label, data = long))
    title <- paste0("Compound-Target Binding Probability -- ", scope_label)
  }

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    cli::cli_abort("Nothing to plot for {.arg what} = {.val {what}} in this scope.")
  }

  if (is.null(width)) width <- max(8, ncol(mat) * 0.9 + 2)
  if (is.null(height)) height <- max(6, nrow(mat) * 0.5 + 2)

  filename <- NA
  path <- NULL
  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("heatmap_", what, "_", scope_label, ".pdf"))
    filename <- path
  }

  result <- pheatmap::pheatmap(
    mat,
    color = grDevices::colorRampPalette(c("white", "#f39c12", "#c0392b"))(50),
    display_numbers = TRUE, number_format = "%.2f", fontsize_number = 9,
    cluster_rows = nrow(mat) > 1, cluster_cols = ncol(mat) > 1,
    main = title, angle_col = 45, fontsize_row = 10, fontsize_col = 10,
    filename = filename, width = width, height = height
  )

  if (save) {
    log_row <- data.frame(condition = scope_label, what = what, path = path, stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "heatmap_plot_log", log_row, c("condition", "what"))
    patliRResults(proj, "heatmap_plot_log") <- log_df
    .write_results_csv(proj, "heatmap_plot_log", log_df)
    attr(result, "proj") <- proj
  }
  result
}
