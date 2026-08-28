#' @include AllGenerics.R internal.R network_build.R network_layers.R
NULL

## network_motifs() -- feed-forward / feedback triple detection on the
## directed layered graph (.network_layered_graph(), network_layers.R). A
## feed-forward loop is a transitive triple A->B->C with A->C also present;
## enumerated as a vectorised join over the edge list (every 2-path, then a
## hashed closing-edge test), not a per-node pair loop -- see
## .network_find_3node_motifs(). A researcher wants the actual
## compound/target/pathway triples, not igraph::motifs()'s isomorphism-
## class counts, so the enumeration lists node identities directly.

#' Feed-forward and feedback loop motifs on the directed layered graph
#'
#' @description
#' Lists the condition's directed compound-target-pathway-disease graph
#' (see `.network_layered_graph()`, internal) for two 3-node patterns
#' (the shapes described by Milo et al. 2002, *Science* 298, 824-827):
#'
#' - **Feed-forward loop**: A -> B, B -> C, A -> C (a transitive triple).
#'   In this graph that is a compound acting on a pathway/disease both
#'   directly (the "derived" edges built by `.network_layered_graph()`) and
#'   indirectly through one of its targets -- e.g. compound X hits pathway
#'   P because target T (which X also hits) participates in P.
#' - **Feedback loop**: A -> B, B -> C, C -> A (a directed 3-cycle).
#'
#' @section This is an occurrence list, not a validated motif set:
#' A "motif" in the Milo et al. sense is a subgraph that occurs
#' significantly more often than in degree-preserving randomised networks.
#' `network_motifs()` does **not** run that null model -- and it could not
#' meaningfully, because most feed-forward loops here are *guaranteed by
#' construction*: `.network_layered_graph()` builds the closing
#' `compound -> pathway` / `compound -> disease` edge as exactly the
#' transitive shortcut of the two-hop path, so every such path is
#' automatically a feed-forward triple. Read the output as "which triples
#' realise this shape", useful for tracing a compound's routes to a
#' pathway, not as evidence that the biology favours feed-forward control.
#'
#' @section Feedback loops are structurally impossible in this graph:
#' Every edge in `.network_layered_graph()` points strictly downstream
#' (compound -> target -> pathway/disease), so a directed 3-cycle cannot
#' exist. The search still runs (cheaply, on the same enumeration) so a
#' future cyclic substrate would be picked up, but on the current graph it
#' always returns zero -- expected, not a bug.
#'
#' @inheritParams network_build
#' @param n_cores Deprecated and ignored. The enumeration is now a
#'   vectorised join over the edge list (all two-paths at once, then a
#'   hashed closing-edge test), fast enough that a `parallel` cluster would
#'   cost more than it saves. Kept only so old calls do not error.
#'
#' @return The updated `proj`, with a `network_motifs` entry in
#'   [patliRResults()] (columns `condition`, `motif_type`
#'   (`"feed_forward"`/`"feedback"`), `node_a`, `node_b`, `node_c`,
#'   `layer_a`, `layer_b`, `layer_c`), also written to
#'   `results/network_motifs.csv`. Conditions with zero motifs of either
#'   type contribute no rows (not a row full of `NA`s) -- see
#'   `projectLog(proj)` for a summary count either way.
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
#' proj <- network_motifs(proj, condition = "FLO-ET")
#' patliRResults(proj, "network_motifs")
#' }
#'
#' @param pathway_db See `.network_target_pathway_edges()`, internal --
#'   restrict the pathway layer to specific `network_enrich()` `db` value(s)
#'   (e.g. `"kegg"` alone, far coarser than GO and so a much smaller
#'   `compound -> pathway` fan-out per compound), or `NULL` (default) for
#'   every `db` combined.
#' @export
network_motifs <- function(proj, condition = NULL, n_cores = 1L, pathway_db = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  if (!missing(n_cores) && !identical(as.integer(n_cores), 1L)) {
    cli::cli_warn("{.arg n_cores} is deprecated and ignored -- the motif enumeration is vectorised now.")
  }
  conditions <- .network_resolve_conditions(proj, condition)

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    g <- .network_layered_graph(proj, cond, pathway_db = pathway_db)
    found <- .network_find_3node_motifs(g)
    rows[[cond]] <- if (nrow(found) == 0) {
      .empty_network_motifs_row()
    } else {
      cbind(condition = cond, found, stringsAsFactors = FALSE)
    }
    proj <- .log_append(
      proj, step = "network_motifs", id = NA_character_,
      message = paste0(
        "condition '", cond, "': ", sum(rows[[cond]]$motif_type == "feed_forward"), " feed-forward, ",
        sum(rows[[cond]]$motif_type == "feedback"), " feedback loop(s) found"
      )
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_motifs", result, "condition")

  patliRResults(proj, "network_motifs") <- result
  .write_results_csv(proj, "network_motifs", result)
  .write_log_csv(proj)
  proj
}

#' Directed 3-node pattern enumeration over the layered graph, as a
#' vectorised join over the edge list
#'
#' @description
#' A feed-forward loop is a transitive triple `A -> B -> C` with `A -> C`
#' present; a feedback loop is `A -> B -> C` with `C -> A` present. Both are
#' found by materialising every two-path once (`merge()` joining each
#' edge's head to the next edge's tail) and testing the closing edge with a
#' single hashed `match()` against the edge set -- work bounded by the
#' number of two-paths, done inside compiled `merge`/`match` rather than an
#' R loop over out-neighbour pairs. This replaced a per-node
#' `choose(out_degree, 2)` search that took minutes on a hub-heavy graph.
#' The edge key joins `from`/`to` with an ASCII `0x01` byte, which cannot
#' occur in a node name, so `"A"` + `"BC"` and `"AB"` + `"C"` never collide.
#' @return `data.frame(motif_type, node_a, node_b, node_c, layer_a, layer_b, layer_c)`.
#' @keywords internal
.network_find_3node_motifs <- function(g) {
  el <- igraph::as_data_frame(g, what = "edges")[, c("from", "to")]
  layer_of <- stats::setNames(igraph::V(g)$layer, igraph::V(g)$name)
  empty <- .empty_network_motifs_row()[, -1]
  if (nrow(el) == 0) return(empty)

  sep <- "\x01"
  edge_key <- paste(el$from, el$to, sep = sep)

  ## every two-path A -> B -> C
  tp <- merge(
    data.frame(a = el$from, b = el$to, stringsAsFactors = FALSE),
    data.frame(b = el$from, c = el$to, stringsAsFactors = FALSE),
    by = "b", sort = FALSE
  )
  tp <- tp[tp$a != tp$c, c("a", "b", "c"), drop = FALSE]
  if (nrow(tp) == 0) return(empty)

  is_ffl <- !is.na(match(paste(tp$a, tp$c, sep = sep), edge_key))  # A -> C closes
  is_fbl <- !is.na(match(paste(tp$c, tp$a, sep = sep), edge_key))  # C -> A closes

  build <- function(sub, type) {
    if (nrow(sub) == 0) return(NULL)
    data.frame(
      motif_type = type,
      node_a = sub$a, node_b = sub$b, node_c = sub$c,
      layer_a = unname(layer_of[sub$a]),
      layer_b = unname(layer_of[sub$b]),
      layer_c = unname(layer_of[sub$c]),
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(build(tp[is_ffl, , drop = FALSE], "feed_forward"),
               build(tp[is_fbl, , drop = FALSE], "feedback"))
  if (is.null(out) || nrow(out) == 0) return(empty)
  rownames(out) <- NULL
  out
}

#' @keywords internal
.empty_network_motifs_row <- function() {
  data.frame(condition = character(0), motif_type = character(0), node_a = character(0),
             node_b = character(0), node_c = character(0), layer_a = character(0),
             layer_b = character(0), layer_c = character(0), stringsAsFactors = FALSE)
}
