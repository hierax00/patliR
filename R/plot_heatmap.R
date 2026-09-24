#' @include AllGenerics.R internal.R network_build.R plot-helpers.R
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
#' @param top_n_targets `NULL` (default, every target) or an integer: for
#'   `what = "compound_target"`, keep only the `top_n_targets` targets hit
#'   by the most compounds (ties by summed weight). With a few hundred
#'   targets the full matrix is unreadable; a `cli_inform` suggests this
#'   above 60 targets. Ignored for `what = "compound_condition"`.
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
                          top_n_targets = NULL,
                          save = TRUE, out_dir = NULL, width = NULL, height = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  what <- match.arg(what)
  stopifnot(is.null(top_n_targets) ||
              (is.numeric(top_n_targets) && length(top_n_targets) == 1 && !is.na(top_n_targets) && top_n_targets >= 1))
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
    rownames(mat) <- .plot_unique_labels(proj, conds, bin$compound_id, "compound")
    storage.mode(mat) <- "double"
    scope_label <- "ALL"
    title <- "Compound x Condition Presence"
  } else {
    conditions <- .network_resolve_conditions(proj, condition)
    scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")
    edges_all <- patliRResults(proj, "network_edges")
    long <- edges_all[edges_all$condition %in% conditions, c("compound_id", "uniprot_id", "weight")]
    if (nrow(long) == 0) {
      cli::cli_abort("Nothing to plot for {.arg what} = {.val {what}} in this scope.")
    }
    n_targets_all <- length(unique(long$uniprot_id))
    if (!is.null(top_n_targets)) {
      long <- long[long$uniprot_id %in% .plot_heatmap_top_targets(long, top_n_targets), , drop = FALSE]
    } else if (n_targets_all > 60) {
      cli::cli_inform(c(
        "i" = "{.fn plot_heatmap}: {n_targets_all} target columns -- the PDF will be {round(n_targets_all * 0.9 + 2)} in wide and the cells unreadable; pass {.arg top_n_targets} (e.g. {.code top_n_targets = 40}) to keep the targets shared by the most compounds."
      ))
    }
    mat <- .plot_edge_matrix(proj, conditions, long)
    title <- paste0(
      "Compound-Target Binding Probability -- ", scope_label,
      if (ncol(mat) < n_targets_all) paste0(" (top ", ncol(mat), " of ", n_targets_all, " targets)")
    )
  }

  ## top_n_targets joined the log key later -- older rows are full matrices
  top_n_used <- if (what == "compound_target" && !is.null(top_n_targets)) as.numeric(top_n_targets) else NA_real_
  proj <- .plot_log_backfill(proj, "heatmap_plot_log", "top_n_targets", NA_real_)
  .plot_pheatmap_render(
    proj, mat, title, save, out_dir, width, height,
    log_name = "heatmap_plot_log",
    log_row = data.frame(condition = scope_label, what = what, top_n_targets = top_n_used,
                         path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "what", "top_n_targets"),
    filename = paste0("heatmap_", what, "_", scope_label, if (!is.na(top_n_used)) paste0("_top", top_n_used), ".pdf")
  )
}

#' Targets shared by the most compounds, for `plot_heatmap(top_n_targets =)`
#'
#' @description
#' Ranked by number of distinct compounds hitting the target, then by the
#' summed `weight` (`NA` counted as 0), then by `uniprot_id`.
#' @return Character vector of at most `top_n` `uniprot_id`s.
#' @keywords internal
.plot_heatmap_top_targets <- function(long, top_n) {
  ids <- sort(unique(long$uniprot_id))
  breadth <- vapply(ids, function(t) length(unique(long$compound_id[long$uniprot_id == t])), integer(1))
  weight <- vapply(ids, function(t) sum(long$weight[long$uniprot_id == t], na.rm = TRUE), numeric(1))
  utils::head(ids[order(-breadth, -weight, ids)], top_n)
}

#' Build by IDs, preserving unknown weights and labelling only after pooling
#' @keywords internal
.plot_edge_matrix <- function(proj, conditions, edges) {
  long <- stats::aggregate(edges["weight"], edges[c("compound_id", "uniprot_id")],
                           function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  compound_ids <- sort(unique(long$compound_id))
  target_ids <- sort(unique(long$uniprot_id))
  mat <- matrix(0, length(compound_ids), length(target_ids), dimnames = list(compound_ids, target_ids))
  mat[cbind(match(long$compound_id, compound_ids), match(long$uniprot_id, target_ids))] <- long$weight
  rownames(mat) <- .plot_unique_labels(proj, conditions, compound_ids, "compound")
  colnames(mat) <- .plot_unique_labels(proj, conditions, target_ids, "target")
  mat
}

#' Render a compound x target/condition matrix with `pheatmap()`, save/log
#' it, and attach the updated project
#'
#' @description
#' The rendering tail shared by [plot_heatmap()] and [plot_rank()]'s
#' `view = "heatmap"` -- both need "hierarchical clustering on both axes,
#' cell values printed, save as PDF, log the path" and differ only in how
#' `mat`/`title`/`filename` are built. Requires `pheatmap` (checked by the
#' caller, not here, so the caller's error message can name itself).
#'
#' @param mat Numeric matrix to plot.
#' @param log_row One-row `data.frame` for the log table; a `path` column
#'   (if present) is overwritten with the resolved path.
#' @param key_cols Passed to [.network_upsert()] for the log upsert.
#' @return The `pheatmap` object, with `attr(., "proj")` set when
#'   `save = TRUE`.
#' @keywords internal
.plot_pheatmap_render <- function(proj, mat, title, save, out_dir, width, height,
                                   log_name, log_row, key_cols, filename) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    cli::cli_abort("Nothing to plot -- the matrix has zero rows or columns.")
  }
  if (is.null(width)) width <- max(8, ncol(mat) * 0.9 + 2)
  if (is.null(height)) height <- max(6, nrow(mat) * 0.5 + 2)

  path <- NULL
  pheatmap_filename <- NA
  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, filename)
    pheatmap_filename <- path
  }

  finite <- mat[is.finite(mat)]
  limits <- if (length(finite)) range(finite) else c(0, 1)
  if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
  clusterable <- function(x) nrow(x) > 1 && all(is.finite(stats::dist(x)))
  result <- pheatmap::pheatmap(
    mat,
    color = grDevices::colorRampPalette(c("white", "#f39c12", "#c0392b"))(50),
    breaks = seq(limits[1], limits[2], length.out = 51), na_col = "grey80",
    display_numbers = TRUE, number_format = "%.2f", fontsize_number = 9,
    cluster_rows = clusterable(mat), cluster_cols = clusterable(t(mat)),
    main = title, angle_col = 45, fontsize_row = 10, fontsize_col = 10,
    filename = pheatmap_filename, width = width, height = height
  )

  if (save) {
    log_row$path <- path
    log_df <- .network_upsert(proj, log_name, log_row, key_cols)
    patliRResults(proj, log_name) <- log_df
    .write_results_csv(proj, log_name, log_df)
    attr(result, "proj") <- proj
  }
  result
}
