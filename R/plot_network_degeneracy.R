#' @include AllGenerics.R internal.R plot-helpers.R plot_network_layers.R network_degeneracy.R
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
#' adds one extra curved edge per compound pair from [network_degeneracy()]
#' selected by `filter` -- visualizing directly on the network which
#' compounds converge on the same targets/pathways, rather than reading
#' `degeneracy_score`/`p_adjusted` off a table with no structural context.
#' Degeneracy edges are drawn as a slight arc (not a straight line) purely
#' so they are visually distinguishable from the underlying compound-target
#' edges even when a pair of compounds happens to sit close together, and
#' coloured by `z_score` -- how far the pair's functional similarity sits
#' above the permutation null, the bounded, significance-bearing quantity
#' (unlike `degeneracy_score`, which is a convenience scalarisation with no
#' fixed scale, see [network_degeneracy()]).
#'
#' @section `filter`: which pairs get a link drawn:
#' `filter = "p_adjusted"` (default) draws every pair with a non-`NA`
#' `p_adjusted < alpha` -- the significance test [network_degeneracy()]'s
#' default `annotation = "direct"` mode actually computes.
#' `filter = "score"` reproduces the historical behaviour: every pair with
#' `degeneracy_score >= min_degeneracy`, an arbitrary threshold on an
#' unbounded, non-significance-tested scalar. `p_adjusted` is `NA` for
#' every row when `network_degeneracy()` was run with
#' `annotation %in% c("enriched", "jaccard")` (piece 14 -- no permutation
#' null exists for those modes); if every pair in scope has `p_adjusted =
#' NA`, `filter = "p_adjusted"` (including the default) transparently
#' falls back to `filter = "score"` and `cli_inform`s why.
#'
#' @inheritParams plot_network_layers
#' @param filter `"p_adjusted"` (default) or `"score"` -- see the section
#'   above.
#' @param min_degeneracy Single number in `[0, 1]`, default `0.3`. Only
#'   used when `filter` resolves to `"score"`: compound pairs with
#'   `degeneracy_score >= min_degeneracy` get a link drawn --
#'   `network_degeneracy()` reports every pair, including near-zero ones,
#'   which would otherwise clutter the plot with effectively meaningless
#'   links.
#' @param alpha Single number in `(0, 1]`, default `0.05`. Only used when
#'   `filter` resolves to `"p_adjusted"`: the significance cutoff on
#'   `p_adjusted`.
#'
#' @return A `ggplot` object (`engine = "static"`) or a `girafe` htmlwidget
#'   (`engine = "ggiraph"`). Same `save`/logging behavior as
#'   [plot_network_layers()], under `patliRResults(proj,
#'   "network_degeneracy_plot_log")` -- keyed (and filenamed) by
#'   `(condition, layout, colour_by, filter)`, same reasoning as
#'   [plot_network_layers()]'s own log.
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
#' # "direct" mode (default) has its own permutation null -- filter = "p_adjusted" works
#' proj <- network_degeneracy(proj, condition = "FLO-ET", seed = 1)
#' plot_network_degeneracy(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_network_degeneracy <- function(proj, condition = NULL, engine = c("static", "ggiraph"),
                                     layout = c("fr", "kk", "drl", "bipartite"),
                                     colour_by = c("layer", "module", "node_type"),
                                     top_hub_n = 15, seed = 1,
                                     filter = c("p_adjusted", "score"),
                                     min_degeneracy = 0.3, alpha = 0.05, save = TRUE, out_dir = NULL,
                                     width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(min_degeneracy), length(min_degeneracy) == 1, min_degeneracy >= 0, min_degeneracy <= 1)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha <= 1)
  engine <- match.arg(engine)
  layout <- match.arg(layout)
  colour_by <- match.arg(colour_by)
  filter <- match.arg(filter)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  degeneracy_all <- patliRResults(proj, "network_degeneracy")
  if (is.null(degeneracy_all) || nrow(degeneracy_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_degeneracy} entry in {.arg proj}.",
      "i" = "Run {.fn network_degeneracy} first."
    ))
  }

  used_filter <- filter
  if (filter == "p_adjusted") {
    in_scope <- degeneracy_all[degeneracy_all$condition %in% conditions, , drop = FALSE]
    has_p_adjusted <- "p_adjusted" %in% names(degeneracy_all) && any(!is.na(in_scope$p_adjusted))
    if (!has_p_adjusted) {
      cli::cli_inform(c(
        "i" = "{.fn plot_network_degeneracy}: {.field p_adjusted} is {.val NA} for every pair in scope (likely {.code annotation != \"direct\"} was used in {.fn network_degeneracy}, which has no permutation null for {.val enriched}/{.val jaccard}) -- falling back to {.code filter = \"score\"} ({.arg min_degeneracy} = {min_degeneracy})."
      ))
      used_filter <- "score"
    }
  }

  ## filter = "score": kept byte-for-byte identical to the pre-revision
  ## expression (including its NA-row subsetting quirk when
  ## degeneracy_score is NA) -- this is the regression-safe path.
  deg <- if (used_filter == "p_adjusted") {
    degeneracy_all[degeneracy_all$condition %in% conditions &
      !is.na(degeneracy_all$p_adjusted) & degeneracy_all$p_adjusted < alpha, , drop = FALSE]
  } else {
    degeneracy_all[degeneracy_all$condition %in% conditions & degeneracy_all$degeneracy_score >= min_degeneracy, , drop = FALSE]
  }

  g <- .network_layered_graph_multi(proj, conditions)
  g <- .network_layers_filter(proj, g, conditions, layers = c("compound", "target"), max_pathways = 0)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort("Condition(s) {.val {conditions}} have zero compound-target nodes -- nothing to plot.")
  }

  pd <- .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed, colour_by = colour_by)
  deg_edges <- .network_degeneracy_edges(pd$nodes, deg)
  p <- .network_degeneracy_ggplot(pd, deg_edges, engine = engine, title_suffix = scope_label)

  .plot_finish(
    proj, p,
    name = "network_degeneracy_plot_log",
    filename = paste0("network_degeneracy_", scope_label, "_", layout, "_", colour_by, "_", used_filter, ".png"),
    log_row = data.frame(
      condition = scope_label, layout = layout, colour_by = colour_by, filter = used_filter,
      path = NA_character_, min_degeneracy = min_degeneracy, alpha = alpha,
      n_links = nrow(deg_edges), stringsAsFactors = FALSE
    ),
    key_cols = c("condition", "layout", "colour_by", "filter"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi,
    static = if (engine == "static") p else .network_degeneracy_ggplot(pd, deg_edges, engine = "static", title_suffix = scope_label)
  )
}

#' Resolve `network_degeneracy()` compound pairs to node coordinates from
#' an already-laid-out `plot_network_layers()`-style node table
#'
#' @return `data.frame(x, y, xend, yend, degeneracy_score, z_score)` --
#'   pairs where either compound is not present in `nodes` (e.g. it had
#'   zero targets in this scope) are silently dropped, not an error, since
#'   `network_degeneracy()` covers every pair in the condition regardless
#'   of whether both ended up in this particular graph. `z_score` is
#'   `NA_real_` for every row when `deg` has no `z_score` column (a `deg`
#'   built by hand, e.g. in tests, or a call site that only ever produces
#'   `degeneracy_score`).
#' @keywords internal
.network_degeneracy_edges <- function(nodes, deg) {
  if (nrow(deg) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                       degeneracy_score = numeric(0), z_score = numeric(0), stringsAsFactors = FALSE))
  }
  idx <- stats::setNames(seq_len(nrow(nodes)), nodes$name)
  ia <- idx[deg$compound_a]
  ib <- idx[deg$compound_b]
  ok <- !is.na(ia) & !is.na(ib)
  z <- if ("z_score" %in% names(deg)) deg$z_score[ok] else rep(NA_real_, sum(ok))
  data.frame(
    x = nodes$x[ia[ok]], y = nodes$y[ia[ok]],
    xend = nodes$x[ib[ok]], yend = nodes$y[ib[ok]],
    degeneracy_score = deg$degeneracy_score[ok],
    z_score = z,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.network_degeneracy_ggplot <- function(pd, deg_edges, engine, title_suffix) {
  p <- .network_layers_ggplot(pd, engine = engine, title_suffix = title_suffix)

  if (nrow(deg_edges) > 0) {
    ## curvature purely so a degeneracy link is visually distinct from the
    ## underlying straight compound-target edges even when the two
    ## compounds happen to be near each other in the layout. Coloured by
    ## z_score (bounded, significance-bearing) rather than sized by the
    ## unbounded degeneracy_score -- `colour` is free here because node
    ## grouping uses `fill`, not `colour` (.network_layers_ggplot()).
    p <- p + ggplot2::geom_curve(
      data = deg_edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                    colour = .data$z_score),
      linewidth = 0.9, alpha = 0.8, curvature = 0.25
    ) +
      ggplot2::scale_colour_gradient(low = "#d2b4de", high = "#4a235a", name = "Degeneracy\nz-score", na.value = "grey60")
  }
  p + ggplot2::labs(
    title = paste0("Compound-target network + degeneracy links -- ", title_suffix),
    subtitle = "Curves = compound pairs converging on shared targets/pathways (network_degeneracy()), coloured by z-score"
  )
}
