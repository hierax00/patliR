#' @include AllGenerics.R internal.R network_build.R network_layers.R
NULL

## network_motifs() -- feed-forward/feedback loop detection on the directed
## layered graph (.network_layered_graph(), network_layers.R). Enumerated
## directly (brute-force over connected triples) rather than via
## igraph::motifs()'s isomorphism-class counts: our graphs are small
## (tens of nodes), and a researcher wants to know *which* compound/target/
## pathway/disease triple forms a motif, not just a count per abstract
## isomorphism class -- igraph::motifs() would need translating its class
## IDs back to actual node triples anyway, with real risk of getting that
## numbering wrong (another "verify, don't assume" case, avoided here by
## just checking the three edges directly).

#' Feed-forward and feedback loop motifs on the directed layered graph
#'
#' @description
#' Searches the condition's directed compound-target-pathway-disease graph
#' (see `.network_layered_graph()`, internal) for two classic 3-node
#' directed network motifs (Milo et al. 2002, *Science* 298, 824-827):
#'
#' - **Feed-forward loop**: A -> B, A -> C, B -> C. In this graph, the only
#'   way this occurs is a compound acting on a pathway/disease both
#'   directly (the "derived" edges built by `.network_layered_graph()`) and
#'   indirectly through one of its targets -- e.g. compound X hits pathway
#'   P both because target T (which X also hits) participates in P.
#' - **Feedback loop**: A -> B, B -> C, C -> A (a directed 3-cycle).
#'
#' @section Feedback loops are structurally impossible in this graph, by construction:
#' Every edge in `.network_layered_graph()` points strictly "downstream"
#' (compound -> target -> {pathway, disease}, plus same-direction derived
#' shortcuts) -- there is no edge type that ever points back toward a
#' compound or from a pathway/disease back to a target. A directed 3-cycle
#' therefore cannot exist in the current layered graph, so
#' `network_motifs()` finding zero feedback loops is an expected, correct
#' result, not evidence of a bug. Feedback loops of the kind the original
#' design envisioned would need a genuinely cyclic substrate (e.g. a
#' protein-protein interaction network, which does have feedback) --
#' `network_bowtie()` runs into the same structural issue, see its own
#' documentation.
#'
#' @inheritParams network_build
#' @param n_cores Integer, default `1L` (sequential, the original
#'   behaviour). The feed-forward-loop search is embarrassingly parallel
#'   over nodes (see `.network_find_3node_motifs()`, internal) -- real
#'   layered graphs can have a handful of hub nodes (compounds with a large
#'   derived `compound -> pathway` fan-out, see [network_enrich()]'s
#'   `simplify_go` for reducing that fan-out at the source) whose
#'   `choose(degree, 2)` pair count dominates total runtime (confirmed:
#'   ~10 min on a real 2838-node/48575-edge Chilcuague network, 2026-08-25).
#'   `n_cores > 1` splits the node loop across a `parallel::makePSOCKcluster()`
#'   cluster (works identically on Windows/macOS/Linux, unlike
#'   fork-based `mclapply()`) -- worth it mainly when a handful of nodes have
#'   very high degree, since those nodes' `O(d^2)` cost dominates and a
#'   PSOCK worker gets one whole such node's pair-enumeration to itself,
#'   not because the *number* of nodes is large. `parallel::detectCores()`
#'   is a reasonable starting point; do not request more workers than the
#'   number of *high-degree* nodes, since idle workers with only
#'   low-degree nodes left add cluster-startup overhead for no real gain.
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
  stopifnot(is.numeric(n_cores), length(n_cores) == 1, n_cores >= 1, n_cores == as.integer(n_cores))
  conditions <- .network_resolve_conditions(proj, condition)

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    g <- .network_layered_graph(proj, cond, pathway_db = pathway_db)
    found <- .network_find_3node_motifs(g, n_cores = as.integer(n_cores))
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

#' Directed 3-node motif search over the layered graph, using precomputed
#' adjacency lists and a hashed edge-existence index instead of repeated
#' full-edge-list scans
#'
#' @description
#' The original implementation called `has_edge(a, b)` (a full `O(E)` scan
#' of the edge list) inside doubly- and triply-nested loops over every
#' node's out-neighbors -- `O(n * d^2 * E)` overall. That is fine for
#' `patliR`'s bundled example (tens of nodes), but real conditions with a
#' `network_enrich()` pathway layer can push some nodes' out-degree well
#' past that assumption, and the brute-force scan cost then dominates
#' (confirmed: ~221s for this one context in a real run). This version
#' precomputes each node's out-neighbors once
#' (`split()`) and answers `has_edge()` via a hashed lookup (an environment
#' used as a hash map, keyed `from`/`to` joined with an ASCII 0x01 control
#' byte -- a character that cannot appear in a node name, so this cannot
#' collide the way plain string concatenation could for names like
#' `"AB"` vs `"A"`+`"B"`). Output (row content and order) is unchanged from the
#' original brute-force version.
#'
#' @section Feed-forward-loop search is the part that gets parallelized:
#' Per-node cost here is `choose(degree(a), 2)` (every pair of `a`'s
#' out-neighbors) -- for most nodes this is tiny, but a handful of hub
#' nodes (e.g. a compound with a large derived `compound -> pathway`
#' fan-out) can have `degree` in the thousands, and `choose(2500, 2)` ~=
#' 3.1 million pairs dominates the whole search all by itself (confirmed
#' against real Chilcuague data, 2026-08-25: ~10 min total, driven by
#' exactly this -- see [network_enrich()]'s `simplify_go` for reducing
#' that fan-out at the source instead). The per-node computation
#' (`ff_worker`, defined below) is embarrassingly parallel -- nodes never
#' share mutable state -- so `n_cores > 1` splits the node loop across a
#' `parallel::makePSOCKcluster()` cluster (works identically on
#' Windows/macOS/Linux, unlike fork-based `mclapply()`). `parLapply()`
#' serializes the worker closure (which captures `out_neighbors`/
#' `has_edge`/`layer_of` by lexical scoping) to each worker automatically,
#' so no manual `clusterExport()` bookkeeping is needed. The feedback-loop
#' search is left sequential -- it is structurally guaranteed empty for
#' this graph (see the exported function's docs) and is cheap in practice
#' (most out-neighbors two hops downstream are terminal pathway/disease
#' nodes with no further out-neighbors), so parallelizing it would only
#' add cluster overhead for no real gain.
#' @return `data.frame(motif_type, node_a, node_b, node_c, layer_a, layer_b, layer_c)`.
#' @keywords internal
.network_find_3node_motifs <- function(g, n_cores = 1L) {
  el <- igraph::as_data_frame(g, what = "edges")[, c("from", "to")]
  layer_of <- stats::setNames(igraph::V(g)$layer, igraph::V(g)$name)
  nodes <- igraph::V(g)$name

  out_neighbors <- split(el$to, el$from)
  edge_index <- new.env(hash = TRUE, parent = emptyenv())
  edge_keys <- paste(el$from, el$to, sep = "")
  for (k in edge_keys) assign(k, TRUE, envir = edge_index)
  has_edge <- function(a, b) exists(paste(a, b, sep = ""), envir = edge_index, inherits = FALSE)

  ff_worker <- function(a) {
    out_a <- out_neighbors[[a]]
    if (is.null(out_a) || length(out_a) < 2) return(NULL)
    pairs <- utils::combn(unique(out_a), 2, simplify = FALSE)
    node_rows <- list()
    for (p in pairs) {
      b <- p[1]; c <- p[2]
      if (has_edge(b, c)) {
        node_rows[[length(node_rows) + 1]] <- data.frame(
          motif_type = "feed_forward", node_a = a, node_b = b, node_c = c,
          layer_a = layer_of[[a]], layer_b = layer_of[[b]], layer_c = layer_of[[c]],
          stringsAsFactors = FALSE
        )
      }
      if (has_edge(c, b)) {
        node_rows[[length(node_rows) + 1]] <- data.frame(
          motif_type = "feed_forward", node_a = a, node_b = c, node_c = b,
          layer_a = layer_of[[a]], layer_b = layer_of[[c]], layer_c = layer_of[[b]],
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(node_rows) == 0) return(NULL)
    do.call(rbind, node_rows)
  }

  ff_per_node <- if (n_cores > 1L && length(nodes) > 1L) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, nodes, ff_worker)
  } else {
    lapply(nodes, ff_worker)
  }
  ff_per_node <- ff_per_node[!vapply(ff_per_node, is.null, logical(1))]
  ff_rows <- if (length(ff_per_node) == 0) list() else list(do.call(rbind, ff_per_node))

  fb_rows <- list()

  for (a in nodes) {
    out_a <- out_neighbors[[a]]
    if (is.null(out_a)) next
    for (b in out_a) {
      out_b <- out_neighbors[[b]]
      if (is.null(out_b)) next
      for (cc in out_b) {
        if (cc != a && has_edge(cc, a)) {
          fb_rows[[length(fb_rows) + 1]] <- data.frame(
            motif_type = "feedback", node_a = a, node_b = b, node_c = cc,
            layer_a = layer_of[[a]], layer_b = layer_of[[b]], layer_c = layer_of[[cc]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  all_rows <- c(ff_rows, fb_rows)
  if (length(all_rows) == 0) return(.empty_network_motifs_row()[, -1])
  do.call(rbind, all_rows)
}

#' @keywords internal
.empty_network_motifs_row <- function() {
  data.frame(condition = character(0), motif_type = character(0), node_a = character(0),
             node_b = character(0), node_c = character(0), layer_a = character(0),
             layer_b = character(0), layer_c = character(0), stringsAsFactors = FALSE)
}
