#' @include AllGenerics.R internal.R network_build.R
NULL

#' Node centrality (degree, betweenness, hub score) per condition
#'
#' @description
#' Computes classic node-level centrality measures over each condition's
#' compound-target graph from [network_build()], for **every node** --
#' compounds and targets alike (`node_type` distinguishes them). This is
#' deliberate: a target's raw degree here is exactly what
#' [network_hub_penalty()] needs as input ("how many compounds touch this
#' target"), and a compound's degree/betweenness is a first, structural
#' read of how central it is to the condition's target profile.
#'
#' @section Unweighted on purpose:
#' [network_build()]'s edges carry a `weight` column, but it is the
#' **import probability** of a compound-target prediction (higher = more
#' confident), not a distance/cost. `igraph::betweenness()` and
#' `igraph::hits_scores()` both default to *using* an edge attribute
#' literally named `weight` if present, treating it as something to
#' minimize along a path -- silently applying that default here would
#' invert the intended meaning (weak/uncertain predictions would look like
#' the "shortest, cheapest" paths). All three measures are computed
#' unweighted (`weights = NA`) to avoid that mismatch. If a future version
#' wants probability-weighted centrality, it needs to explicitly transform
#' `weight` into a proper distance (e.g. `1 - probability`) first, not pass
#' the raw column through.
#'
#' @section Hub score is eigenvector centrality here:
#' The graph is **undirected**, so Kleinberg's hub and authority scores
#' coincide and both equal the principal eigenvector of the adjacency
#' matrix -- i.e. `hub_score` is eigenvector centrality, not a directional
#' "hub vs. authority" distinction (that only exists on a directed graph).
#' It is kept under the name `hub_score` for continuity with earlier
#' versions. Implementation note: `igraph::hub_score()` is deprecated as of
#' igraph 2.0.3; the code uses `igraph::hits_scores(g, weights = NA)$hub`,
#' numerically the same.
#'
#' On a **disconnected** graph the principal eigenvector of the whole
#' adjacency matrix concentrates on the component with the largest spectral
#' radius, so every node in every other component gets `hub_score` ~ 0
#' regardless of how central it is *within its own component*. When
#' `normalize = TRUE` an additional `hub_score_component` column is
#' computed by running the HITS eigenvector **per connected component** and
#' rescaling each to `max = 1` inside that component; `component_id`
#' records which component each node is in. Use `hub_score` for a
#' whole-graph ranking, `hub_score_component` to compare nodes across
#' components of different size.
#'
#' @section Bipartite-aware normalisation (`normalize = TRUE`, the default):
#' This is a two-mode (bipartite) graph -- a compound's degree ranges over
#' `[0, |T|]`, a target's over `[0, |C|]`, and `|T|` is usually much larger,
#' so a joint ranking of the raw `degree` (or `betweenness`) column mostly
#' reflects which mode a node is in, not its importance. `normalize = TRUE`
#' adds mode-aware `_norm` siblings that all land in `[0, 1]` and *are*
#' comparable across `node_type`, following Borgatti & Everett (1997):
#'
#' * `degree_norm` = `degree(v) / |T|` for a compound, `degree(v) / |C|`
#'   for a target -- the opposite mode's size is the structural maximum.
#' * `betweenness_norm` = `betweenness(v) / B_max`, with `B_max` the
#'   mode-specific maximum betweenness of a bipartite graph of these mode
#'   sizes (Borgatti & Everett 1997; cross-checked against NetworkX's
#'   `bipartite.betweenness_centrality`). `igraph::betweenness()` on an
#'   undirected graph already counts each unordered pair once, so no
#'   factor-of-two correction is applied to the raw score.
#'
#' The raw `degree` / `betweenness` / `hub_score` columns are always kept.
#' Two per-condition constants, `n_compounds` and `n_targets`, are added so
#' the normalisation is self-documenting in the CSV.
#'
#' `normalize = FALSE` reproduces the pre-normalisation output exactly
#' (`condition`, `node_id`, `node_type`, plus one column per requested
#' measure) -- an escape hatch for code that does `names(result)` checks.
#'
#' @section Comparable across `node_type`, not across conditions:
#' `|C|` and `|T|` are per-condition constants, so a `_norm` value is
#' comparable between a compound and a target *of the same condition* but
#' **not** between two conditions of different size -- exactly as the
#' R-index of [network_module_robustness()] is not comparable across
#' modules of different size. A disconnected graph makes `betweenness_norm`
#' a slight *under*-estimate (cross-component pairs contribute no paths,
#' while `B_max` assumes connectivity).
#'
#' @inheritParams network_build
#' @param measures Character vector, any of `"degree"`, `"betweenness"`,
#'   `"hub_score"` (default: all three).
#' @param normalize Logical (default `TRUE`). When `TRUE`, adds the
#'   bipartite-aware `_norm` columns for whichever `measures` were
#'   requested, plus `hub_score_component` / `component_id` (if
#'   `"hub_score"` was requested) and the `n_compounds` / `n_targets`
#'   constants. When `FALSE`, the output is exactly the raw columns.
#'
#' @return The updated `proj`, with a `network_centrality` entry in
#'   [patliRResults()] (columns `condition`, `node_id`, `node_type`
#'   (`"compound"`/`"target"`), one column per requested measure, and -- when
#'   `normalize = TRUE` -- `degree_norm` / `betweenness_norm` /
#'   `hub_score_component` / `component_id` for the requested measures plus
#'   `n_compounds` / `n_targets`), also written to
#'   `results/network_centrality.csv`.
#'
#' @references
#' Borgatti SP, Everett MG (1997). "Network analysis of 2-mode data."
#' *Social Networks* 19(3):243-269.
#' \doi{10.1016/S0378-8733(96)00301-2}
#'
#' Freeman LC (1978). "Centrality in social networks: conceptual
#' clarification." *Social Networks* 1(3):215-239.
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
#' proj <- network_centrality(proj)
#' patliRResults(proj, "network_centrality")
#' }
#'
#' @export
network_centrality <- function(proj, condition = NULL,
                                measures = c("degree", "betweenness", "hub_score"),
                                normalize = TRUE) {
  stopifnot(is(proj, "PatliRProject"))
  measures <- match.arg(measures, c("degree", "betweenness", "hub_score"), several.ok = TRUE)
  stopifnot(is.logical(normalize), length(normalize) == 1L, !is.na(normalize))
  conditions <- .network_resolve_conditions(proj, condition)

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    g <- .network_graph(proj, cond)
    rows[[cond]] <- .network_centrality_one(g, cond, measures, normalize)
  }

  result <- do.call(rbind, rows)
  if (is.null(result) || nrow(result) == 0) result <- .empty_network_centrality_row()
  rownames(result) <- NULL
  result <- .network_upsert(
    proj, "network_centrality", result, "condition",
    touched_keys = data.frame(condition = conditions, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "network_centrality") <- result
  .write_results_csv(proj, "network_centrality", result)
  proj <- .log_append(
    proj, step = "network_centrality", id = NA_character_,
    message = paste0("condition '", conditions, "': centrality computed for ",
                      vapply(conditions, function(cond) sum(result$condition == cond), integer(1)),
                      " nodes", if (normalize) " (bipartite-normalised)" else "")
  )
  .write_log_csv(proj)
  proj
}

#' One condition's centrality table
#' @keywords internal
.network_centrality_one <- function(g, cond, measures, normalize) {
  node_ids  <- igraph::V(g)$name
  is_target <- as.logical(igraph::V(g)$type)
  node_type <- ifelse(is_target, "target", "compound")
  n_comp <- sum(!is_target)
  n_targ <- sum(is_target)

  df <- data.frame(condition = cond, node_id = node_ids, node_type = node_type,
                   stringsAsFactors = FALSE)

  deg <- igraph::degree(g)
  btw <- if ("betweenness" %in% measures) igraph::betweenness(g, weights = NA) else NULL

  if ("degree" %in% measures) df$degree <- deg
  if ("betweenness" %in% measures) df$betweenness <- btw
  if ("hub_score" %in% measures) df$hub_score <- igraph::hits_scores(g, weights = NA)$hub

  if (!normalize) return(df)

  ## --- Borgatti & Everett (1997) bipartite normalisation --------------
  if ("degree" %in% measures) {
    dn <- rep(NA_real_, length(node_ids))
    if (n_targ > 0) dn[!is_target] <- deg[!is_target] / n_targ   # compound / |T|
    if (n_comp > 0) dn[is_target]  <- deg[is_target]  / n_comp   # target   / |C|
    df$degree_norm <- dn
  }

  if ("betweenness" %in% measures) {
    bn     <- rep(NA_real_, length(node_ids))
    bmax_u <- .network_bipartite_bmax(n_comp, n_targ)   # compound mode
    bmax_v <- .network_bipartite_bmax(n_targ, n_comp)   # target mode
    if (isTRUE(is.finite(bmax_u)) && bmax_u > 0) bn[!is_target] <- btw[!is_target] / bmax_u
    if (isTRUE(is.finite(bmax_v)) && bmax_v > 0) bn[is_target]  <- btw[is_target]  / bmax_v
    df$betweenness_norm <- bn
    zero_u <- isTRUE(is.finite(bmax_u)) && bmax_u == 0 && n_comp > 0
    zero_v <- isTRUE(is.finite(bmax_v)) && bmax_v == 0 && n_targ > 0
    if (zero_u || zero_v) {
      cli::cli_inform(c(
        "i" = "Condition {.val {cond}}: bipartite {.field B_max} is 0 for a mode of size 1; {.field betweenness_norm} is {.val NA} for that mode."
      ))
    }
  }

  if ("hub_score" %in% measures) {
    hc <- .network_hub_score_component(g)
    df$hub_score_component <- hc$hub
    df$component_id        <- hc$membership
  }

  df$n_compounds <- n_comp
  df$n_targets   <- n_targ
  df
}

#' Borgatti & Everett (1997) maximum betweenness for one mode of a
#' bipartite graph
#'
#' @description
#' `B_max` for the mode with `n_mode` nodes, given the opposite mode has
#' `n_other` nodes. Cross-checked 2026-08-28 against
#' `networkx.algorithms.bipartite.centrality.betweenness_centrality`
#' (which attributes the maxima to Borgatti & Everett). NetworkX halves its
#' raw scores before dividing by these maxima *because* its base routine
#' counts ordered pairs; `igraph::betweenness()` on an undirected graph
#' already counts each unordered pair once, so the `0.5` here is the only
#' factor and the raw score is passed through unmodified.
#'
#' @param n_mode Size of the mode being normalised.
#' @param n_other Size of the opposite mode.
#' @return A single numeric `B_max`, or `NA_real_` if either mode is empty.
#' @keywords internal
.network_bipartite_bmax <- function(n_mode, n_other) {
  if (n_mode < 1L || n_other < 1L) return(NA_real_)
  s <- (n_mode - 1L) %/% n_other
  t <- (n_mode - 1L) %%  n_other
  0.5 * (
    n_other^2 * (s + 1)^2 +
      n_other * (s + 1) * (2 * t - s - 1) -
      t * (2 * s - t + 3)
  )
}

#' HITS hub score recomputed per connected component
#'
#' @description
#' Fixes the disconnected-graph pathology of a whole-graph principal
#' eigenvector: runs `igraph::hits_scores()` on each connected component
#' and rescales the hub vector to `max = 1` within that component. A
#' singleton (isolated node) component gets `hub = 1` for its one node.
#'
#' @param g An `igraph` object.
#' @return `list(membership = <integer>, hub = <double>)`, both in `V(g)`
#'   order.
#' @keywords internal
.network_hub_score_component <- function(g) {
  nv <- igraph::vcount(g)
  if (nv == 0L) return(list(membership = integer(0), hub = numeric(0)))
  comp <- igraph::components(g)
  hub  <- numeric(nv)
  for (k in seq_len(comp$no)) {
    idx <- which(comp$membership == k)
    sub <- igraph::induced_subgraph(g, idx)
    h   <- igraph::hits_scores(sub, weights = NA)$hub
    mx  <- suppressWarnings(max(h))
    hub[idx] <- if (isTRUE(is.finite(mx)) && mx > 0) h / mx else rep(1, length(idx))
  }
  list(membership = as.integer(comp$membership), hub = hub)
}

#' Maximal (all-measures, normalised) empty schema for `network_centrality`
#' @keywords internal
.empty_network_centrality_row <- function() {
  data.frame(
    condition = character(0), node_id = character(0), node_type = character(0),
    degree = double(0), betweenness = double(0), hub_score = double(0),
    degree_norm = double(0), betweenness_norm = double(0),
    hub_score_component = double(0), component_id = integer(0),
    n_compounds = integer(0), n_targets = integer(0),
    stringsAsFactors = FALSE
  )
}
