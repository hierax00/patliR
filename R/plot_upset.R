#' @include AllGenerics.R internal.R network_build.R
NULL

## plot_upset() -- deliberately NOT UpSetR/ComplexUpset (both would be new,
## fairly heavy Suggests: UpSetR pulls in its own plotting stack, ComplexUpset
## pulls in ComplexHeatmap's dependency tree). An UpSet plot is really just
## a bar chart (intersection sizes) stacked on a dot-matrix (which sets are
## in each intersection) sharing an x-axis -- entirely reasonable to
## hand-roll in plain ggplot2 the same way plot_admet_radar()'s spiderweb
## and plot_network_layers()'s hairball were. The one new Suggests this
## needs is `patchwork`, purely to stack the two ggplot2 panels with an
## aligned x-axis -- patchwork's own dependency footprint is just ggplot2
## + grid, about as light as a "combine two ggplots" package can be.

#' UpSet-style plot: target set overlaps across conditions
#'
#' @description
#' Which targets are shared across which conditions (extracts), and how
#' many. Builds a `condition x target` presence matrix from
#' `network_edges`, groups targets by their exact combination of
#' conditions ("intersection"), and shows the `top_n` largest
#' intersections as a bar chart (size) stacked on a dot-matrix (which
#' conditions make up that intersection) -- the standard UpSet layout
#' (Lex et al. 2014), hand-rolled in `ggplot2` + `patchwork` rather than
#' pulling in `UpSetR`/`ComplexUpset`.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param top_n Integer, default `15`. Only the `top_n` largest
#'   intersections are shown -- with more than a handful of conditions the
#'   number of possible combinations grows fast.
#'
#' @return A `patchwork` object (bar chart on top of the dot-matrix,
#'   `print()`-able and `ggsave()`-able like a single `ggplot`). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "upset_plot_log")`.
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
#' plot_upset(proj, save = FALSE)
#' }
#'
#' @export
plot_upset <- function(proj, condition = NULL, top_n = 15,
                        save = TRUE, out_dir = NULL, width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_upset}.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_upset} needs the {.pkg patchwork} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"patchwork\")} -- pure CRAN, ggplot2 + grid only."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)
  if (length(conditions) < 2) {
    cli::cli_abort("{.fn plot_upset} needs at least 2 conditions in scope to compare; got {.val {conditions}}.")
  }
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  edges_all <- patliRResults(proj, "network_edges")
  dat <- unique(edges_all[edges_all$condition %in% conditions, c("condition", "uniprot_id")])
  membership <- lapply(split(dat$condition, dat$uniprot_id), function(x) conditions %in% x)
  combo <- vapply(membership, function(x) paste(as.integer(x), collapse = ""), character(1))
  members_by_combo <- membership[match(unique(combo), combo)]
  names(members_by_combo) <- unique(combo)
  sizes <- as.data.frame(table(combo), stringsAsFactors = FALSE)
  names(sizes) <- c("combo", "n_targets")
  sizes <- sizes[order(-sizes$n_targets), , drop = FALSE]
  sizes <- utils::head(sizes, top_n)
  sizes$combo <- factor(sizes$combo, levels = sizes$combo)

  dot <- do.call(rbind, lapply(as.character(sizes$combo), function(cb) {
    data.frame(combo = cb, condition = conditions, in_set = members_by_combo[[cb]], stringsAsFactors = FALSE)
  }))
  dot$combo <- factor(dot$combo, levels = levels(sizes$combo))

  p_bar <- ggplot2::ggplot(sizes, ggplot2::aes(x = .data$combo, y = .data$n_targets)) +
    ggplot2::geom_col(fill = "#2980b9") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_targets), vjust = -0.4, size = 3) +
    ggplot2::labs(title = paste0("Target set intersections -- ", scope_label), x = NULL, y = "Targets") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )

  p_dot <- ggplot2::ggplot(dot, ggplot2::aes(x = .data$combo, y = .data$condition)) +
    ggplot2::geom_line(data = dot[dot$in_set, , drop = FALSE], ggplot2::aes(group = .data$combo), colour = "grey30", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$in_set), size = 3) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey85"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())

  p <- patchwork::wrap_plots(p_bar, p_dot, ncol = 1, heights = c(2, 1))

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("upset_", scope_label, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
    log_row <- data.frame(condition = scope_label, path = path, n_intersections = nrow(sizes), stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "upset_plot_log", log_row, "condition")
    patliRResults(proj, "upset_plot_log") <- log_df
    .write_results_csv(proj, "upset_plot_log", log_df)
  }

  if (save) attr(p, "proj") <- proj
  p
}
