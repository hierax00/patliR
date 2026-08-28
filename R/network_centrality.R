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
#' @section Degree and betweenness are comparable only within `node_type`:
#' This is a bipartite graph -- a compound's degree ranges over `[0, |T|]`,
#' a target's over `[0, |C|]`, and `|T|` is usually much larger, so a joint
#' ranking of the raw `degree` (or `betweenness`) column mostly reflects
#' which mode a node is in, not its importance. Compare within
#' `node_type`. A bipartite-aware normalisation (Borgatti & Everett 1997)
#' is planned (see `ROADMAP.md`).
#'
#' @inheritParams network_build
#' @param measures Character vector, any of `"degree"`, `"betweenness"`,
#'   `"hub_score"` (default: all three).
#'
#' @return The updated `proj`, with a `network_centrality` entry in
#'   [patliRResults()] (columns `condition`, `node_id`, `node_type`
#'   (`"compound"`/`"target"`), plus one column per requested measure), also
#'   written to `results/network_centrality.csv`.
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
                                measures = c("degree", "betweenness", "hub_score")) {
  stopifnot(is(proj, "PatliRProject"))
  measures <- match.arg(measures, c("degree", "betweenness", "hub_score"), several.ok = TRUE)
  conditions <- .network_resolve_conditions(proj, condition)

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    g <- .network_graph(proj, cond)
    node_ids <- igraph::V(g)$name
    node_type <- ifelse(igraph::V(g)$type, "target", "compound")

    df <- data.frame(condition = cond, node_id = node_ids, node_type = node_type, stringsAsFactors = FALSE)
    if ("degree" %in% measures) df$degree <- igraph::degree(g)
    if ("betweenness" %in% measures) df$betweenness <- igraph::betweenness(g, weights = NA)
    if ("hub_score" %in% measures) df$hub_score <- igraph::hits_scores(g, weights = NA)$hub
    rows[[cond]] <- df
  }

  result <- do.call(rbind, rows)
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
                      vapply(conditions, function(cond) sum(result$condition == cond), integer(1)), " nodes")
  )
  .write_log_csv(proj)
  proj
}
