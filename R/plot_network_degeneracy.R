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
#' @section Views and layouts:
#' `view = "network"` (default) draws the links over the compound-target
#' network. With a few hundred targets the force-directed layouts (`"fr"`,
#' `"kk"`, `"drl"`) bury the compounds in the middle of the target mesh;
#' `layout = "circle"` drops the targets and puts the compounds on a
#' circle (ordered so mutually degenerate compounds sit side by side), which
#' makes every link readable as a chord. `view = "heatmap"` is a compound x
#' compound tile map of every scored pair's `z_score` (the pairs passing
#' `filter` marked `*`), the densest summary when many pairs pass.
#' Links are coloured by `z_score` and their width follows
#' `degeneracy_score`.
#'
#' @section A condition with no `network_degeneracy()` rows is an error:
#' [network_degeneracy()] is run per condition. Plotting a condition it was
#' never run on aborts (naming the conditions it *was* run on) rather than
#' silently drawing the bare network; when rows exist but none passes
#' `filter`, the figure is drawn with a warning.
#'
#' @inheritParams plot_network_layers
#' @param layout As in [plot_network_layers()], plus `"circle"` --
#'   compounds only, on a circle; see the section above.
#' @param view `"network"` (default) or `"heatmap"` -- see the section
#'   above. `layout`, `colour_by`, `top_hub_n`, `seed` and `engine` only
#'   apply to the network view.
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
#'   `(condition, layout, colour_by, filter, view)`, same reasoning as
#'   [plot_network_layers()]'s own log (the heatmap is saved as
#'   `network_degeneracy_heatmap_<condition>_<filter>.png`).
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
#' # compounds on a circle, targets hidden -- every link visible
#' plot_network_degeneracy(proj, condition = "FLO-ET", layout = "circle", save = FALSE)
#' # compound x compound z-score heatmap
#' plot_network_degeneracy(proj, condition = "FLO-ET", view = "heatmap", save = FALSE)
#' }
#'
#' @export
plot_network_degeneracy <- function(proj, condition = NULL, engine = c("static", "ggiraph"),
                                     layout = c("fr", "kk", "drl", "bipartite", "circle"),
                                     colour_by = c("layer", "module", "node_type"),
                                     top_hub_n = 15, seed = 1,
                                     filter = c("p_adjusted", "score"),
                                     min_degeneracy = 0.3, alpha = 0.05,
                                     view = c("network", "heatmap"), save = TRUE, out_dir = NULL,
                                     width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(min_degeneracy), length(min_degeneracy) == 1, min_degeneracy >= 0, min_degeneracy <= 1)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha <= 1)
  engine <- match.arg(engine)
  layout <- match.arg(layout)
  colour_by <- match.arg(colour_by)
  filter <- match.arg(filter)
  view <- match.arg(view)
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
  ## network_degeneracy() is run per condition; a figure for a condition it
  ## was never run on would silently draw the bare network with no link at
  ## all -- indistinguishable from "no degenerate pair". Refuse instead.
  in_scope <- degeneracy_all[degeneracy_all$condition %in% conditions, , drop = FALSE]
  if (nrow(in_scope) == 0) {
    available <- sort(unique(degeneracy_all$condition))
    cli::cli_abort(c(
      "No {.val network_degeneracy} rows for condition(s) {.val {conditions}} -- there is nothing to draw.",
      "i" = "{.fn network_degeneracy} has been run for: {.val {available}}.",
      "i" = "Pass one of those as {.arg condition}, or run {.code network_degeneracy(proj, condition = ...)} first."
    ))
  }

  used_filter <- filter
  if (filter == "p_adjusted") {
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
  n_pass <- sum(!is.na(deg$compound_a))
  if (n_pass == 0) {
    cli::cli_warn(c(
      "!" = "None of the {nrow(in_scope)} compound pair{?s} in scope passes {.code filter = \"{used_filter}\"} ({if (used_filter == 'p_adjusted') paste0('p_adjusted < ', alpha) else paste0('degeneracy_score >= ', min_degeneracy)}) -- no degeneracy link is drawn.",
      "i" = "Relax {.arg alpha}/{.arg min_degeneracy}, or use {.code view = \"heatmap\"} to see every pair's z-score."
    ))
  }

  ## `view` joined the log key after the first release -- older logs are
  ## all network views.
  proj <- .plot_log_backfill(proj, "network_degeneracy_plot_log", "view", "network")

  if (view == "heatmap") {
    p <- .network_degeneracy_heatmap(proj, in_scope, deg, conditions, scope_label, used_filter, fig_width = width)
    return(.plot_finish(
      proj, p,
      name = "network_degeneracy_plot_log",
      filename = paste0("network_degeneracy_heatmap_", scope_label, "_", used_filter, ".png"),
      log_row = data.frame(
        condition = scope_label, layout = NA_character_, colour_by = NA_character_, filter = used_filter,
        view = "heatmap", path = NA_character_, min_degeneracy = min_degeneracy, alpha = alpha,
        n_links = n_pass, stringsAsFactors = FALSE
      ),
      key_cols = c("condition", "layout", "colour_by", "filter", "view"),
      engine = NULL, save = save, out_dir = out_dir,
      width = width, height = height, dpi = dpi
    ))
  }

  g <- .network_layered_graph_multi(proj, conditions)
  g <- .network_layers_filter(proj, g, conditions, layers = c("compound", "target"), max_pathways = 0)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort("Condition(s) {.val {conditions}} have zero compound-target nodes -- nothing to plot.")
  }

  pd <- if (layout == "circle") {
    .network_degeneracy_circle_data(proj, g, conditions, in_scope, colour_by = colour_by, label = top_hub_n > 0)
  } else {
    .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed, colour_by = colour_by)
  }
  deg_edges <- .network_degeneracy_edges(pd$nodes, deg)
  p <- .network_degeneracy_ggplot(pd, deg_edges, engine = engine, title_suffix = scope_label, fig_width = width)

  .plot_finish(
    proj, p,
    name = "network_degeneracy_plot_log",
    filename = paste0("network_degeneracy_", scope_label, "_", layout, "_", colour_by, "_", used_filter, ".png"),
    log_row = data.frame(
      condition = scope_label, layout = layout, colour_by = colour_by, filter = used_filter,
      view = "network", path = NA_character_, min_degeneracy = min_degeneracy, alpha = alpha,
      n_links = nrow(deg_edges), stringsAsFactors = FALSE
    ),
    key_cols = c("condition", "layout", "colour_by", "filter", "view"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi,
    static = if (engine == "static") p else .network_degeneracy_ggplot(pd, deg_edges, engine = "static", title_suffix = scope_label, fig_width = width)
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

#' Order compounds so that mutually degenerate ones sit next to each other
#'
#' @description
#' Average-linkage clustering of the pairwise `z_score` matrix
#' (`degeneracy_score` when no `z_score` exists), distance = max - value,
#' missing pairs treated as the least similar. Used for the heatmap's axes
#' and the circle layout's order, so a block of convergent compounds shows
#' up as a block / as short chords. Fewer than 3 compounds: alphabetical.
#' @param in_scope `network_degeneracy` rows in scope.
#' @param ids Compound IDs to order.
#' @return `ids`, reordered.
#' @keywords internal
.network_degeneracy_order <- function(in_scope, ids) {
  ids <- sort(unique(ids))
  if (length(ids) < 3) return(ids)
  value_col <- if ("z_score" %in% names(in_scope) && any(is.finite(in_scope$z_score))) "z_score" else "degeneracy_score"
  d <- in_scope[!is.na(in_scope$compound_a) & in_scope$compound_a %in% ids & in_scope$compound_b %in% ids, , drop = FALSE]
  v <- d[[value_col]]
  if (length(v) == 0 || !any(is.finite(v))) return(ids)
  m <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  ## pooled conditions: a pair present in several keeps its strongest value
  for (i in seq_len(nrow(d))) {
    if (!is.finite(v[i])) next
    a <- d$compound_a[i]
    b <- d$compound_b[i]
    m[a, b] <- m[b, a] <- max(m[a, b], v[i], na.rm = TRUE)
  }
  top <- max(m, na.rm = TRUE)
  low <- min(m, na.rm = TRUE)
  dist_m <- top - m
  dist_m[is.na(dist_m)] <- top - low + 1
  diag(dist_m) <- 0
  hc <- stats::hclust(stats::as.dist(dist_m), method = "average")
  ids[hc$order]
}

#' `layout = "circle"`: compounds only, evenly spaced on a circle
#'
#' @description
#' The force-directed compound-target layout buries the compounds (and the
#' degeneracy links between them) in the middle of a few hundred target
#' nodes. This view drops the targets altogether: each compound is a node on
#' a unit circle, ordered by `.network_degeneracy_order()`, sized by its
#' number of targets, and the degeneracy links become chords. Same return
#' shape as `.network_layers_plot_data()`, with no compound-target edges;
#' `nodes$label_x`/`label_y`/`label_hjust` place each label radially just
#' outside the circle.
#' @keywords internal
.network_degeneracy_circle_data <- function(proj, g, conditions, in_scope, colour_by, label = TRUE) {
  is_compound <- igraph::V(g)$layer == "compound"
  ids <- .network_degeneracy_order(in_scope, igraph::V(g)$name[is_compound])
  degree <- igraph::degree(g, mode = "all")[match(ids, igraph::V(g)$name)]
  theta <- pi / 2 - 2 * pi * (seq_along(ids) - 1) / length(ids) # 12 o'clock, clockwise
  nodes <- data.frame(
    name = ids, layer = "compound", degree = as.integer(degree),
    x = cos(theta), y = sin(theta), stringsAsFactors = FALSE
  )
  nodes$label <- .plot_label_nodes(proj, conditions, nodes$name, nodes$layer)
  nodes$is_hub <- FALSE # labels are drawn radially by .network_degeneracy_ggplot()
  nodes$show_label <- label
  nodes$label_x <- 1.1 * nodes$x
  nodes$label_y <- 1.1 * nodes$y
  nodes$label_hjust <- ifelse(abs(nodes$x) < 0.05, 0.5, ifelse(nodes$x > 0, 0, 1))
  nodes$label_vjust <- ifelse(abs(nodes$x) < 0.05, ifelse(nodes$y > 0, 0, 1), 0.5)
  colour_info <- .network_layers_colour_info(proj, conditions, nodes, colour_by)
  nodes$colour_group <- colour_info$group
  edges <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                      edge_kind = character(0), is_derived = logical(0), stringsAsFactors = FALSE)
  list(nodes = nodes, edges = edges, palette = colour_info$palette,
       legend_name = colour_info$legend_name, circle = TRUE,
       xlim = c(-2, 2), ylim = c(-1.25, 1.25)) # room for the radial labels
}

#' @keywords internal
.network_degeneracy_ggplot <- function(pd, deg_edges, engine, title_suffix, fig_width = 9) {
  circle <- isTRUE(pd$circle)
  ## the base compound-target edges are dimmed under the links so the
  ## links -- the point of this figure -- are what the eye lands on
  p <- .network_layers_ggplot(pd, engine = engine, title_suffix = title_suffix, fig_width = fig_width,
                              edge_alpha = if (nrow(deg_edges) > 0) 0.5 else 1)

  if (nrow(deg_edges) > 0) {
    ## curvature purely so a degeneracy link is visually distinct from the
    ## underlying straight compound-target edges even when the two
    ## compounds happen to be near each other in the layout. Coloured by
    ## z_score (bounded, significance-bearing), width by degeneracy_score
    ## -- `colour` is free here because node grouping uses `fill`, not
    ## `colour` (.network_layers_ggplot()). A saturated plasma ramp, not a
    ## pale one: the links must stand out against the grey edge mesh.
    deg_edges <- deg_edges[order(deg_edges$z_score, na.last = FALSE), , drop = FALSE] # strongest on top
    p <- p + ggplot2::geom_curve(
      data = deg_edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                    colour = .data$z_score, linewidth = .data$degeneracy_score),
      alpha = 0.85, curvature = if (circle) 0.15 else 0.25, lineend = "round"
    ) +
      ggplot2::scale_colour_viridis_c(option = "plasma", begin = 0.05, end = 0.85, direction = -1,
                                      name = "Degeneracy\nz-score", na.value = "grey45") +
      ggplot2::scale_linewidth(range = c(0.4, 1.8), name = "Degeneracy\nscore")
  }

  if (circle && any(pd$nodes$show_label)) {
    lab <- pd$nodes
    lab$short_label <- .plot_truncate(lab$label)
    p <- p + ggplot2::geom_text(
      data = lab,
      ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data$short_label,
                   hjust = .data$label_hjust, vjust = .data$label_vjust),
      size = 3, colour = "grey15", show.legend = FALSE
    )
  }

  subtitle <- if (circle) {
    "Compounds on a circle (ordered so mutually degenerate compounds are adjacent; node size = number of targets, targets not drawn); chords = compound pairs converging on shared targets/pathways (network_degeneracy()), colour = z-score, width = degeneracy score"
  } else {
    "Curves = compound pairs converging on shared targets/pathways (network_degeneracy()), colour = z-score, width = degeneracy score"
  }
  p + ggplot2::labs(
    title = paste0(if (circle) "Compound degeneracy links -- " else "Compound-target network + degeneracy links -- ", title_suffix),
    subtitle = .plot_wrap(subtitle, .plot_wrap_width(fig_width, 8))
  )
}

#' `view = "heatmap"`: compound x compound tile map of the degeneracy z-score
#'
#' @description
#' Every pair `network_degeneracy()` scored in scope, not only those
#' passing the filter -- the filter's pairs are marked with `*`. Axes are
#' ordered by `.network_degeneracy_order()`; the diagonal is blank. Falls
#' back to `degeneracy_score` as the fill when no `z_score` exists (the
#' `enriched`/`jaccard` annotations). One facet per condition when several
#' are pooled.
#' @keywords internal
.network_degeneracy_heatmap <- function(proj, in_scope, deg, conditions, scope_label, used_filter, fig_width = 9) {
  in_scope <- in_scope[!is.na(in_scope$compound_a) & !is.na(in_scope$compound_b), , drop = FALSE]
  use_z <- "z_score" %in% names(in_scope) && any(is.finite(in_scope$z_score))
  value_col <- if (use_z) "z_score" else "degeneracy_score"
  ids <- .network_degeneracy_order(in_scope, c(in_scope$compound_a, in_scope$compound_b))
  labels <- .plot_truncate(.plot_unique_labels(proj, conditions, ids, "compound"))

  pass_key <- paste(deg$condition, deg$compound_a, deg$compound_b, sep = "\r")
  tiles <- data.frame(
    condition = rep(in_scope$condition, 2),
    row_id = c(in_scope$compound_a, in_scope$compound_b),
    col_id = c(in_scope$compound_b, in_scope$compound_a),
    value = rep(in_scope[[value_col]], 2),
    pass = rep(paste(in_scope$condition, in_scope$compound_a, in_scope$compound_b, sep = "\r") %in% pass_key, 2),
    stringsAsFactors = FALSE
  )
  tiles$row_id <- factor(tiles$row_id, levels = rev(ids))
  tiles$col_id <- factor(tiles$col_id, levels = ids)
  label_of <- stats::setNames(labels, ids)

  fill_scale <- if (use_z) {
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                                  na.value = "grey90", name = "Degeneracy\nz-score")
  } else {
    ggplot2::scale_fill_gradient(low = "white", high = "#b2182b", na.value = "grey90", name = "Degeneracy\nscore")
  }
  subtitle <- paste0(
    "Cell = ", if (use_z) "z-score of the pair's functional similarity against network_degeneracy()'s permutation null" else "degeneracy_score (no z-score for this annotation mode)",
    "; * = pair passes filter = \"", used_filter, "\"; axes ordered so mutually degenerate compounds are adjacent"
  )

  p <- ggplot2::ggplot(tiles, ggplot2::aes(x = .data$col_id, y = .data$row_id, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(data = tiles[tiles$pass, , drop = FALSE], label = "*", size = 4, colour = "grey10", vjust = 0.75) +
    fill_scale +
    ggplot2::scale_x_discrete(labels = label_of, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = label_of, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("Compound degeneracy -- ", scope_label),
      subtitle = .plot_wrap(subtitle, .plot_wrap_width(fig_width, 8)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7.5),
      axis.text.y = ggplot2::element_text(size = 7.5),
      panel.grid = ggplot2::element_blank()
    )
  if (length(unique(tiles$condition)) > 1) p <- p + ggplot2::facet_wrap(~condition)
  p
}
