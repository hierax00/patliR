#' @include AllGenerics.R internal.R plot_network_layers.R network_degeneracy.R
NULL

## Draws network_degeneracy()'s compound-compound convergence links on top
## of the compound-target graph. Reuses plot_network_layers()'s layout
## pipeline rather than recomputing one -- same base network, one extra
## edge layer.

#' Compound-target network with `network_degeneracy()` convergence links
#' overlaid
#'
#' @description
#' Draws the same compound-target layout [plot_network_layers()] does
#' (reuses its layout code directly, so node positions match exactly), then
#' adds one extra curved edge per compound pair from
#' [network_degeneracy()] whose `degeneracy_score >= min_degeneracy` --
#' visualizing directly on the network which compounds converge on the same
#' targets/pathways, rather than reading `degeneracy_score` off a table
#' with no structural context.
#' Degeneracy edges are drawn as a slight arc (not a straight line) purely
#' so they are visually distinguishable from the underlying compound-target
#' edges even when a pair of compounds happens to sit close together.
#'
#' @inheritParams plot_network_layers
#' @param min_degeneracy Single number in `[0, 1]`, default `0.3`. Only
#'   compound pairs with `degeneracy_score >= min_degeneracy` get a link
#'   drawn -- `network_degeneracy()` reports every pair, including
#'   near-zero ones, which would otherwise clutter the plot with
#'   effectively meaningless links.
#'
#' @return A `ggplot` object (`engine = "static"`) or a `girafe` htmlwidget
#'   (`engine = "ggiraph"`). Same `save`/logging behavior as
#'   [plot_network_layers()], under `patliRResults(proj,
#'   "network_degeneracy_plot_log")`.
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
#' proj <- network_enrich(proj, condition = "FLO-ET", db = "go")
#' proj <- network_degeneracy(proj, condition = "FLO-ET")
#' plot_network_degeneracy(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_network_degeneracy <- function(proj, condition = NULL, engine = c("static", "ggiraph"),
                                     layout = c("fr", "kk", "drl"), top_hub_n = 15, seed = 1,
                                     min_degeneracy = 0.3, save = TRUE, out_dir = NULL,
                                     width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(min_degeneracy), length(min_degeneracy) == 1, min_degeneracy >= 0, min_degeneracy <= 1)
  engine <- match.arg(engine)
  layout <- match.arg(layout)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_network_degeneracy}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  degeneracy_all <- patliRResults(proj, "network_degeneracy")
  if (is.null(degeneracy_all) || nrow(degeneracy_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_degeneracy} entry in {.arg proj}.",
      "i" = "Run {.fn network_degeneracy} first."
    ))
  }
  deg <- degeneracy_all[degeneracy_all$condition %in% conditions & degeneracy_all$degeneracy_score >= min_degeneracy, , drop = FALSE]

  g <- .network_layered_graph_multi(proj, conditions)
  g <- .network_layers_filter(proj, g, conditions, layers = c("compound", "target"), max_pathways = 0)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort("Condition(s) {.val {conditions}} have zero compound-target nodes -- nothing to plot.")
  }

  pd <- .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed)
  deg_edges <- .network_degeneracy_edges(pd$nodes, deg)
  p <- .network_degeneracy_ggplot(pd, deg_edges, engine = engine, title_suffix = scope_label)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    p_static <- if (engine == "static") p else .network_degeneracy_ggplot(pd, deg_edges, engine = "static", title_suffix = scope_label)
    path <- file.path(out_dir, paste0("network_degeneracy_", scope_label, ".png"))
    ggplot2::ggsave(path, p_static, width = width, height = height, dpi = dpi)
    log_row <- data.frame(
      condition = scope_label, path = path, min_degeneracy = min_degeneracy,
      n_links = nrow(deg_edges), stringsAsFactors = FALSE
    )
    log_df <- .network_upsert(proj, "network_degeneracy_plot_log", log_row, "condition")
    patliRResults(proj, "network_degeneracy_plot_log") <- log_df
    .write_results_csv(proj, "network_degeneracy_plot_log", log_df)
  }

  result <- if (engine == "static") {
    p
  } else {
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  }
  if (save) attr(result, "proj") <- proj
  result
}

#' Resolve `network_degeneracy()` compound pairs to node coordinates from
#' an already-laid-out `plot_network_layers()`-style node table
#'
#' @return `data.frame(x, y, xend, yend, degeneracy_score)` -- pairs where
#'   either compound is not present in `nodes` (e.g. it had zero targets in
#'   this scope) are silently dropped, not an error, since
#'   `network_degeneracy()` covers every pair in the condition regardless
#'   of whether both ended up in this particular graph.
#' @keywords internal
.network_degeneracy_edges <- function(nodes, deg) {
  if (nrow(deg) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                       degeneracy_score = numeric(0), stringsAsFactors = FALSE))
  }
  idx <- stats::setNames(seq_len(nrow(nodes)), nodes$name)
  ia <- idx[deg$compound_a]
  ib <- idx[deg$compound_b]
  ok <- !is.na(ia) & !is.na(ib)
  data.frame(
    x = nodes$x[ia[ok]], y = nodes$y[ia[ok]],
    xend = nodes$x[ib[ok]], yend = nodes$y[ib[ok]],
    degeneracy_score = deg$degeneracy_score[ok],
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.network_degeneracy_ggplot <- function(pd, deg_edges, engine, title_suffix) {
  p <- .network_layers_ggplot(pd, engine = engine, title_suffix = title_suffix)

  if (nrow(deg_edges) > 0) {
    ## curvature purely so a degeneracy link is visually distinct from the
    ## underlying straight compound-target edges even when the two
    ## compounds happen to be near each other in the layout.
    p <- p + ggplot2::geom_curve(
      data = deg_edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                    linewidth = .data$degeneracy_score),
      colour = "#8e44ad", alpha = 0.6, curvature = 0.25
    ) +
      ggplot2::scale_linewidth(range = c(0.3, 1.8), name = "Degeneracy\nscore")
  }
  p + ggplot2::labs(
    title = paste0("Compound-target network + degeneracy links -- ", title_suffix),
    subtitle = "Purple curves = compound pairs converging on shared targets/pathways (network_degeneracy())"
  )
}
