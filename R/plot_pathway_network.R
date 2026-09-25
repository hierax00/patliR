#' @include AllGenerics.R internal.R network_build.R network_module_robustness.R plot-helpers.R plot_network_layers.R
NULL

## Two figures for the pathway level of the analysis, which until now only
## existed as tables (network_enrichment, network_kegg_topology) and as
## KEGG's own pre-drawn diagrams (network_pathview()):
##
##  * plot_pathway_network() -- an "enrichment map": the most significant
##    enriched terms as nodes, linked when their target gene sets overlap,
##    clusters of related terms outlined. Answers "which pathways does the
##    extract hit, and how do they connect".
##  * plot_kegg_topology() -- the directed gene-gene relations KEGG itself
##    draws inside a pathway (activation, inhibition, phosphorylation, ...),
##    parsed by network_kegg_topology(), with the extract's targets
##    highlighted.

# ---------------------------------------------------------------------------
# plot_pathway_network()
# ---------------------------------------------------------------------------

#' Enrichment map: enriched pathways as a network of overlapping gene sets
#'
#' @description
#' Draws the `top_n` most significant terms of one database in
#' [network_enrich()]'s results as a network: one node per term, an edge
#' between two terms when the extract's target genes annotated to them
#' overlap (Jaccard or overlap coefficient at or above `min_similarity`).
#' Node area is the number of the extract's targets in the term (`Count`),
#' node colour its significance (\eqn{-\log_{10}} adjusted p), edge width
#' the overlap. Groups of related terms are found by community detection
#' (Louvain) and outlined, each with a cluster name taken from its most
#' central term, so redundant terms (for instance a family of cancer
#' pathways that all share the same PI3K/MAPK genes) read as one theme.
#'
#' Optionally, the extract's compounds can be added as a second node type
#' (`show_compounds = TRUE`), and the terms whose genes are enriched in a
#' disease gene set can be ringed (`disease`).
#'
#' @section Gene-set overlap:
#' For two terms with target-gene sets \eqn{A} and \eqn{B} (the `geneID`
#' column of `network_enrichment`, i.e. only the extract's targets in each
#' term, not the full pathway):
#' \deqn{J(A,B) = |A \cap B| / |A \cup B| \quad \textrm{(Jaccard)}}
#' \deqn{O(A,B) = |A \cap B| / \min(|A|, |B|) \quad \textrm{(overlap coefficient)}}
#' The overlap coefficient reaches 1 when one set is contained in the other
#' (useful for GO, where child terms nest inside parents); Jaccard penalises
#' size differences. Merico et al. (2010) use a 0.5 overlap-coefficient
#' cutoff, or a combination of both measures, to draw an edge.
#'
#' @section Disease-gene enrichment (`disease`):
#' For each drawn term, a one-sided hypergeometric test asks whether its
#' target genes contain more genes of the disease set than expected. The
#' universe is the extract's own targets in `condition` (mapped to Entrez,
#' `N` genes, `K` of them disease genes); a term with `n` targets, `k` of
#' them disease genes, gets \eqn{p = P(X \ge k)},
#' \eqn{X \sim \mathrm{Hypergeometric}(N, K, n)}. p-values are adjusted
#' (Benjamini-Hochberg) across the drawn terms; a term with adjusted p
#' below `disease_fdr` is drawn with a black ring. Every term with at least
#' one disease gene also gets a dark inner dot whose area is the share
#' \eqn{k / n} of the node's area, so the disease load stays visible when
#' no term reaches significance. Because the universe is
#' the extract's targets, not the genome, the test asks "within what this
#' extract hits, is this pathway especially disease-related", which does not
#' re-test the pathway's enrichment itself.
#'
#' @section Compounds (`show_compounds`):
#' A compound is linked to a term when at least `min_compound_hits` of its
#' predicted targets (`network_edges`, UniProt mapped to Entrez) are among
#' the term's genes. The `top_n_compounds` compounds with the most (term,
#' gene) hits over the drawn terms are placed, evenly spaced, on an outer
#' ring, in the circular order of the centre of the terms they hit. Each
#' compound keeps only its `max_links_per_compound` strongest links (largest
#' share of the term's genes hit): most compounds hit most terms a little,
#' and drawing every link hides which terms a compound hits hardest.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param condition Character scalar, a condition with rows in
#'   `network_enrichment` for `db`. `NULL` (default) works only when exactly
#'   one condition has enrichment results for `db`.
#' @param db One of `"kegg"` (default), `"reactome"`, `"go"`.
#' @param top_n Integer, default `30`: number of most significant terms
#'   (lowest adjusted p) drawn.
#' @param p_cutoff Numeric in (0, 1], default `0.05`: terms with adjusted p
#'   above this are never drawn.
#' @param similarity `"jaccard"` (default) or `"overlap"`, see the section
#'   above.
#' @param min_similarity Minimum overlap for an edge. `NULL` (default)
#'   means `0.3` for Jaccard and `0.5` for the overlap coefficient.
#' @param layout `"kk"` (Kamada-Kawai, default: spreads nodes evenly, keeps clusters apart)
#'   or `"fr"` (Fruchterman-Reingold), both with edge weights = overlap and a fixed `seed`.
#' @param clusters Logical, default `TRUE`: outline Louvain clusters of
#'   terms (clusters of a single term are not outlined).
#' @param show_compounds Logical, default `FALSE`: add the extract's
#'   compounds as a second node type (see the section above).
#' @param top_n_compounds Integer, default `8`.
#' @param min_compound_hits Integer, default `2`.
#' @param max_links_per_compound Integer, default `5`.
#' @param disease `NULL` (default) or a single `disease_id` present in
#'   `patliRResults(proj, "disease_genes")` (see [disease_genes_fetch()],
#'   [disease_genes_import()]).
#' @param disease_fdr Numeric in (0, 1], default `0.05`.
#' @param label_width Integer, default `24`: characters per line of the
#'   wrapped term labels (labels are also cut at three lines).
#' @param seed Integer seed for the layout and the community detection;
#'   the caller's random state is restored afterwards.
#' @param engine `"static"` (default) or `"ggiraph"` (hover tooltips with
#'   the full term name, genes and statistics).
#'
#' @return A `ggplot` (or `girafe`) object with attributes `terms` (one row
#'   per drawn term: statistics, cluster, degree and, with `disease`, the
#'   hypergeometric result), `edges` (term pairs, `similarity`,
#'   `n_shared`) and, with `show_compounds`, `compound_links`. With
#'   `save = TRUE` the PNG is written and logged in
#'   `patliRResults(proj, "pathway_network_plot_log")` (keyed by
#'   `condition`, `db`, `disease`, `compounds`), and the updated project is
#'   in `attr(., "proj")`.
#'
#' @references
#' Merico, D., Isserlin, R., Stueker, O., Emili, A. & Bader, G. D. (2010),
#' "Enrichment Map: a network-based method for gene-set enrichment
#' visualization and interpretation", *PLoS ONE* 5(11), e13984.
#' \doi{10.1371/journal.pone.0013984}
#'
#' Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. (2008),
#' "Fast unfolding of communities in large networks", *Journal of
#' Statistical Mechanics* P10008. \doi{10.1088/1742-5468/2008/10/P10008}
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... build the project, run network_build() and
#' # network_enrich(condition = "FLO-ET", db = "kegg") ...
#' plot_pathway_network(proj, condition = "FLO-ET", db = "kegg", save = FALSE)
#' # rings on the pathways enriched in hypertension genes
#' proj <- disease_genes_fetch(proj, "hypertension")
#' plot_pathway_network(proj, condition = "FLO-ET", disease = "MONDO_0005044")
#' }
#'
#' @seealso [network_enrich()], [plot_enrichment()], [plot_kegg_topology()]
#' @export
plot_pathway_network <- function(proj, condition = NULL, db = c("kegg", "reactome", "go"), top_n = 30,
                                 p_cutoff = 0.05, similarity = c("jaccard", "overlap"), min_similarity = NULL,
                                 layout = c("kk", "fr"), clusters = TRUE,
                                 show_compounds = FALSE, top_n_compounds = 8, min_compound_hits = 2,
                                 max_links_per_compound = 5, disease = NULL, disease_fdr = 0.05, label_width = 24, seed = 1,
                                 engine = c("static", "ggiraph"),
                                 save = TRUE, out_dir = NULL, width = 11, height = 8.5, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  db <- match.arg(db)
  similarity <- match.arg(similarity)
  layout <- match.arg(layout)
  engine <- match.arg(engine)
  .pathway_check_count(top_n, "top_n", min = 1)
  .pathway_check_prob(p_cutoff, "p_cutoff")
  .pathway_check_prob(disease_fdr, "disease_fdr")
  if (is.null(min_similarity)) min_similarity <- if (similarity == "jaccard") 0.3 else 0.5
  if (!is.numeric(min_similarity) || length(min_similarity) != 1L || is.na(min_similarity) ||
      min_similarity < 0 || min_similarity > 1) {
    cli::cli_abort("{.arg min_similarity} must be a single number in [0, 1], or NULL.")
  }
  .pathway_check_flag(clusters, "clusters")
  .pathway_check_flag(show_compounds, "show_compounds")
  .pathway_check_count(top_n_compounds, "top_n_compounds", min = 1)
  .pathway_check_count(min_compound_hits, "min_compound_hits", min = 1)
  .pathway_check_count(max_links_per_compound, "max_links_per_compound", min = 1)
  .pathway_check_count(label_width, "label_width", min = 8)
  .pathway_check_count(seed, "seed", min = -Inf)
  if (!is.null(disease) && (!is.character(disease) || length(disease) != 1L || is.na(disease) || !nzchar(disease))) {
    cli::cli_abort("{.arg disease} must be NULL or a single {.field disease_id} string.")
  }
  engine <- .plot_require(engine)

  enr_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enr_all) || nrow(enr_all) == 0) {
    cli::cli_abort(c("No {.val network_enrichment} entry in {.arg proj}.", "i" = "Run {.fn network_enrich} first."))
  }
  enr_db <- enr_all[enr_all$db == db, , drop = FALSE]
  if (nrow(enr_db) == 0) {
    cli::cli_abort(c(
      "No {.code db = \"{db}\"} rows in {.val network_enrichment}.",
      "i" = "Available: {.val {unique(enr_all$db)}}. Run {.fn network_enrich} with {.code db = \"{db}\"} first."
    ))
  }
  cond <- .pathway_single_condition(unique(enr_db$condition), condition, "plot_pathway_network",
                                    paste0("network_enrichment (db = \"", db, "\")"))

  terms <- .pathway_net_select_terms(enr_db, cond, top_n, p_cutoff)
  if (nrow(terms) == 0) {
    cli::cli_abort("No {.val {db}} term of condition {.val {cond}} has adjusted p <= {p_cutoff}.")
  }
  sets <- .pathway_net_gene_sets(terms$geneID, terms$ID)
  sim <- .pathway_net_overlap_matrix(sets, similarity)
  edges <- .pathway_net_edges(sim, sets, min_similarity)

  restore_rng <- .with_seed(seed)
  on.exit(restore_rng(), add = TRUE)
  membership <- if (clusters) .pathway_net_clusters(terms$ID, edges) else stats::setNames(rep(NA_integer_, nrow(terms)), terms$ID)
  aspect <- max(0.6, (width - 3.4) / max(1, height - 1.8))
  coords <- .pathway_net_layout(terms$ID, edges, layout, aspect = aspect)
  restore_rng()

  terms$x <- coords[terms$ID, "x"]
  terms$y <- coords[terms$ID, "y"]
  terms$cluster <- unname(membership[terms$ID])
  terms$degree <- vapply(terms$ID, function(id) sum(edges$from == id | edges$to == id), integer(1))
  terms$neg_log10_padj <- -log10(pmax(terms$p.adjust, .Machine$double.xmin))
  cluster_info <- .pathway_net_cluster_names(terms, edges)

  ## Entrez mapping of the extract's targets: needed for the disease test's
  ## universe and for the compound links.
  map <- NULL
  if (show_compounds || !is.null(disease)) {
    map <- .pathway_target_entrez(proj, cond)
  }

  dis <- NULL
  if (!is.null(disease)) {
    dis <- .pathway_net_disease(proj, disease, sets, map, disease_fdr)
    terms <- merge(terms, dis$table, by = "ID", sort = FALSE)
    terms <- terms[order(terms$p.adjust), , drop = FALSE]
  }

  links <- NULL
  cmp_nodes <- NULL
  if (show_compounds) {
    links <- .pathway_net_compound_links(map, sets, min_hits = min_compound_hits, top_n = top_n_compounds,
                                         max_links = max_links_per_compound)
    if (nrow(links) == 0) {
      cli::cli_warn("No compound hits at least {min_compound_hits} gene{?s} of a drawn term; no compound drawn.")
    } else {
      cmp_nodes <- .pathway_net_compound_ring(links, terms)
      cmp_nodes$label <- .plot_truncate(.plot_unique_labels(proj, cond, cmp_nodes$compound_id, "compound"), 26)
    }
  }

  p <- .pathway_net_ggplot(
    terms = terms, edges = edges, cluster_info = cluster_info, links = links, cmp_nodes = cmp_nodes,
    dis = dis, db = db, cond = cond, similarity = similarity, min_similarity = min_similarity,
    label_width = label_width, engine = engine, width = width, seed = seed
  )

  attr(p, "terms") <- .pathway_net_term_table(terms, cluster_info)
  attr(p, "edges") <- edges
  if (!is.null(links)) attr(p, "compound_links") <- links

  disease_tag <- if (is.null(disease)) "none" else disease
  filename <- paste0(
    "pathway_network_", cond, "_", db,
    if (!is.null(disease)) paste0("_", gsub("[^A-Za-z0-9]+", "-", disease)) else "",
    if (show_compounds) "_compounds" else "", ".png"
  )
  log_row <- data.frame(
    condition = cond, db = db, disease = disease_tag, compounds = show_compounds,
    n_terms = nrow(terms), n_edges = nrow(edges),
    n_clusters = nrow(cluster_info), similarity = similarity, min_similarity = min_similarity,
    n_disease_enriched = if (is.null(dis)) NA_integer_ else sum(terms$disease_enriched),
    path = NA_character_, stringsAsFactors = FALSE
  )
  .plot_finish(
    proj, p,
    name = "pathway_network_plot_log", filename = filename, log_row = log_row,
    key_cols = c("condition", "db", "disease", "compounds"),
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Argument checks shared by the two pathway plots
#' @keywords internal
.pathway_check_count <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min || x != round(x)) {
    cli::cli_abort("{.arg {arg}} must be a single whole number{if (is.finite(min)) paste0(' >= ', min) else ''}.")
  }
  invisible(TRUE)
}

#' @rdname dot-pathway_check_count
#' @keywords internal
.pathway_check_prob <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x > 1) {
    cli::cli_abort("{.arg {arg}} must be a single number in (0, 1].")
  }
  invisible(TRUE)
}

#' @rdname dot-pathway_check_count
#' @keywords internal
.pathway_check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) cli::cli_abort("{.arg {arg}} must be TRUE or FALSE.")
  invisible(TRUE)
}

#' Resolve a single condition among those that have results
#'
#' @description
#' Pathway-level figures are per condition (the same term has different
#' statistics in each), so `NULL` is only accepted when exactly one
#' condition has results -- the rule [network_pathview()] and
#' [plot_gochord()] use.
#' @param available Conditions that have rows in the table.
#' @param what Name of the table, for the error message.
#' @return A single condition name.
#' @keywords internal
.pathway_single_condition <- function(available, condition, fn, what) {
  available <- sort(unique(available))
  if (is.null(condition)) {
    if (length(available) != 1L) {
      cli::cli_abort(c(
        "{.fn {fn}} draws one condition at a time; {.arg condition} is required.",
        "i" = "{what} has results for {.val {available}}."
      ))
    }
    return(available)
  }
  if (!is.character(condition) || length(condition) != 1L || is.na(condition)) {
    cli::cli_abort("{.arg condition} must be a single condition name, or NULL.")
  }
  if (!condition %in% available) {
    cli::cli_abort(c(
      "No {what} rows for condition {.val {condition}}.",
      "i" = "Available: {.val {available}}."
    ))
  }
  condition
}

#' The `top_n` most significant terms of one condition (one row per ID)
#' @return `enr` rows of `cond`, `p.adjust <= p_cutoff`, lowest `p.adjust`
#'   first, duplicated IDs dropped, at most `top_n` rows.
#' @keywords internal
.pathway_net_select_terms <- function(enr, cond, top_n, p_cutoff) {
  enr <- enr[enr$condition == cond & !is.na(enr$p.adjust) & enr$p.adjust <= p_cutoff, , drop = FALSE]
  enr <- enr[!is.na(enr$geneID) & nzchar(enr$geneID), , drop = FALSE]
  enr <- enr[order(enr$p.adjust, enr$pvalue, -enr$Count), , drop = FALSE]
  enr <- enr[!duplicated(enr$ID), , drop = FALSE]
  enr <- utils::head(enr, top_n)
  rownames(enr) <- NULL
  enr
}

#' Split `/`-separated `geneID` strings into a named list of gene sets
#' @keywords internal
.pathway_net_gene_sets <- function(gene_id, ids) {
  sets <- lapply(as.character(gene_id), function(s) {
    g <- trimws(strsplit(s, "/", fixed = TRUE)[[1]])
    unique(g[nzchar(g) & !is.na(g)])
  })
  stats::setNames(sets, ids)
}

#' Pairwise overlap between gene sets
#'
#' @param sets Named list of character vectors.
#' @param method `"jaccard"` (\eqn{|A \cap B| / |A \cup B|}) or
#'   `"overlap"` (\eqn{|A \cap B| / \min(|A|, |B|)}).
#' @return Symmetric numeric matrix, `dimnames = names(sets)`, diagonal 1
#'   (0 for an empty set).
#' @keywords internal
.pathway_net_overlap_matrix <- function(sets, method = c("jaccard", "overlap")) {
  method <- match.arg(method)
  ids <- names(sets)
  n <- length(sets)
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  if (n == 0) return(m)
  genes <- sort(unique(unlist(sets, use.names = FALSE)))
  inc <- matrix(0L, n, length(genes))
  for (i in seq_len(n)) inc[i, match(sets[[i]], genes)] <- 1L
  shared <- inc %*% t(inc)
  size <- lengths(sets)
  denom <- if (method == "jaccard") {
    outer(size, size, `+`) - shared
  } else {
    outer(size, size, pmin)
  }
  m[] <- ifelse(denom > 0, shared / denom, 0)
  dimnames(m) <- list(ids, ids)
  m
}

#' Term-term edges at or above a similarity threshold
#' @return `data.frame(from, to, similarity, n_shared)`, one row per
#'   unordered pair (`from` listed before `to` in `names(sets)`).
#' @keywords internal
.pathway_net_edges <- function(sim, sets, min_similarity) {
  empty <- data.frame(from = character(0), to = character(0), similarity = numeric(0),
                      n_shared = integer(0), stringsAsFactors = FALSE)
  if (nrow(sim) < 2) return(empty)
  idx <- which(upper.tri(sim) & sim >= min_similarity & sim > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) return(empty)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  ids <- rownames(sim)
  out <- data.frame(from = ids[idx[, 1]], to = ids[idx[, 2]], similarity = sim[idx],
                    stringsAsFactors = FALSE)
  out$n_shared <- mapply(function(a, b) length(intersect(sets[[a]], sets[[b]])), out$from, out$to, USE.NAMES = FALSE)
  out
}

#' Louvain communities of the term graph
#' @return Named integer vector (names = `ids`): cluster number, `1` = the
#'   largest cluster (ties: the one holding the most significant term, i.e.
#'   the earliest in `ids`); `NA` for terms alone in their cluster.
#' @keywords internal
.pathway_net_clusters <- function(ids, edges) {
  out <- stats::setNames(rep(NA_integer_, length(ids)), ids)
  if (nrow(edges) == 0) return(out)
  g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                                     vertices = data.frame(name = ids, stringsAsFactors = FALSE))
  igraph::E(g)$weight <- edges$similarity
  memb <- igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
  memb <- stats::setNames(as.integer(memb), igraph::V(g)$name)[ids]
  sizes <- table(memb)
  multi <- as.integer(names(sizes)[sizes >= 2])
  if (length(multi) == 0) return(out)
  first_pos <- vapply(multi, function(k) min(which(memb == k)), integer(1))
  ord <- multi[order(-as.integer(sizes[as.character(multi)]), first_pos)]
  out[memb %in% multi] <- match(memb[memb %in% multi], ord)
  out
}

#' Layout of the term graph, scaled to a `aspect`:1 box
#'
#' @description
#' Each connected component is laid out on its own (FR or KK, weights =
#' overlap) and components are packed by [igraph::layout_components()]
#' (terms with no edge land around the main components). The result is
#' rescaled to `x` in `[0, 10 * aspect]` and `y` in `[0, 10]`.
#' @return Matrix with columns `x`, `y`, rownames `ids`.
#' @keywords internal
.pathway_net_layout <- function(ids, edges, layout = "fr", aspect = 1.3) {
  g <- igraph::graph_from_data_frame(
    if (nrow(edges)) edges[, c("from", "to")] else data.frame(from = character(0), to = character(0)),
    directed = FALSE, vertices = data.frame(name = ids, stringsAsFactors = FALSE)
  )
  if (nrow(edges)) igraph::E(g)$weight <- edges$similarity
  n <- length(ids)
  xy <- if (n == 1) {
    matrix(c(0, 0), 1)
  } else if (layout == "kk") {
    igraph::layout_components(g, layout = function(h) {
      if (igraph::vcount(h) < 2) return(matrix(0, igraph::vcount(h), 2))
      w <- igraph::E(h)$weight
      igraph::layout_with_kk(h, weights = if (length(w)) 1 / (w + 0.05) else NULL)
    })
  } else {
    igraph::layout_components(g, layout = function(h) {
      if (igraph::vcount(h) < 2) return(matrix(0, igraph::vcount(h), 2))
      igraph::layout_with_fr(h, niter = 2000, weights = igraph::E(h)$weight)
    })
  }
  xy <- .pathway_rescale_xy(xy, 10 * aspect, 10)
  rownames(xy) <- igraph::V(g)$name
  colnames(xy) <- c("x", "y")
  xy[ids, , drop = FALSE]
}

#' Rescale each column of a two-column coordinate matrix to `[0, w]` and
#' `[0, h]` (a constant column goes to the middle)
#' @keywords internal
.pathway_rescale_xy <- function(xy, w, h) {
  xy <- as.matrix(xy)
  resc <- function(v, to) {
    r <- range(v)
    if (!is.finite(diff(r)) || diff(r) == 0) return(rep(to / 2, length(v)))
    (v - r[1]) / diff(r) * to
  }
  cbind(resc(xy[, 1], w), resc(xy[, 2], h))
}

#' Name each term cluster after its most central term
#' @return `data.frame(cluster, id, n_terms, representative, label)` --
#'   `representative` is the Description of the term with the largest summed
#'   within-cluster overlap (ties: lowest adjusted p); `label` the legend
#'   text, naming the two most central terms.
#' @keywords internal
.pathway_net_cluster_names <- function(terms, edges) {
  ks <- sort(unique(stats::na.omit(terms$cluster)))
  if (length(ks) == 0) {
    return(data.frame(cluster = integer(0), id = character(0), n_terms = integer(0),
                      representative = character(0), label = character(0), stringsAsFactors = FALSE))
  }
  rows <- lapply(ks, function(k) {
    members <- terms$ID[!is.na(terms$cluster) & terms$cluster == k]
    e <- edges[edges$from %in% members & edges$to %in% members, , drop = FALSE]
    strength <- vapply(members, function(id) sum(e$similarity[e$from == id | e$to == id]), numeric(1))
    padj <- terms$p.adjust[match(members, terms$ID)]
    top <- members[order(-strength, padj)]
    rep_desc <- terms$Description[match(top[1], terms$ID)]
    shown <- sub(" (signaling )?pathway$", "", terms$Description[match(utils::head(top, 2), terms$ID)])
    data.frame(cluster = k, id = paste0("C", k), n_terms = length(members), representative = rep_desc,
               label = paste0("C", k, "  ", paste(.plot_truncate(shown, 26), collapse = " / "),
                              " (", length(members), " terms)"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Map a condition's targets (UniProt) to Entrez gene IDs
#' @return `data.frame(compound_id, uniprot_id, ENTREZID)`.
#' @keywords internal
.pathway_target_entrez <- function(proj, cond) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.arg show_compounds}/{.arg disease} need {.pkg clusterProfiler} and {.pkg org.Hs.eg.db} (Bioconductor)",
      "i" = "They map the extract's UniProt targets to the Entrez IDs of {.field network_enrichment$geneID}."
    ))
  }
  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c("No {.val network_edges} entry in {.arg proj}.", "i" = "Run {.fn network_build} first."))
  }
  ct <- unique(edges_all[edges_all$condition == cond & !is.na(edges_all$uniprot_id) & nzchar(edges_all$uniprot_id),
                         c("compound_id", "uniprot_id"), drop = FALSE])
  if (nrow(ct) == 0) cli::cli_abort("No compound-target edges for condition {.val {cond}}.")
  map <- .pathway_bitr(unique(ct$uniprot_id), "UNIPROT", "ENTREZID")
  if (nrow(map) == 0) cli::cli_abort("None of condition {.val {cond}}'s targets mapped to an Entrez ID.")
  out <- merge(ct, map, by.x = "uniprot_id", by.y = "UNIPROT")
  unique(out[, c("compound_id", "uniprot_id", "ENTREZID")])
}

#' Quiet `clusterProfiler::bitr()` against `org.Hs.eg.db`
#' @return `data.frame` with columns named `from` and `to`; zero rows on
#'   failure or when the packages are missing.
#' @keywords internal
.pathway_bitr <- function(ids, from, to) {
  empty <- stats::setNames(data.frame(character(0), character(0), stringsAsFactors = FALSE), c(from, to))
  ids <- unique(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) == 0 || !requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(empty)
  }
  out <- tryCatch(
    suppressMessages(suppressWarnings(
      clusterProfiler::bitr(ids, fromType = from, toType = to, OrgDb = "org.Hs.eg.db", drop = TRUE)
    )),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0) return(empty)
  out[, c(from, to)]
}

#' Hypergeometric enrichment of each term's genes in a disease gene set
#'
#' @param term_sets Named list of Entrez gene sets (one per term).
#' @param disease_genes Character vector of Entrez IDs of the disease set.
#' @param universe Character vector of Entrez IDs (the extract's targets);
#'   term genes outside it are added to it.
#' @return `data.frame(ID, disease_k, term_n, disease_K, universe_N,
#'   disease_expected, disease_fold, disease_p, disease_padj,
#'   disease_genes)`; `disease_padj` is Benjamini-Hochberg across the terms.
#' @keywords internal
.pathway_net_disease_test <- function(term_sets, disease_genes, universe) {
  universe <- unique(c(universe, unlist(term_sets, use.names = FALSE)))
  dis <- intersect(unique(disease_genes), universe)
  N <- length(universe)
  K <- length(dis)
  n <- lengths(term_sets)
  k <- vapply(term_sets, function(s) length(intersect(s, dis)), integer(1))
  p <- stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
  p[k == 0] <- 1
  expected <- n * K / max(N, 1)
  data.frame(
    ID = names(term_sets), disease_k = k, term_n = n, disease_K = K, universe_N = N,
    disease_expected = expected, disease_fold = ifelse(expected > 0, k / expected, NA_real_),
    disease_p = p, disease_padj = stats::p.adjust(p, method = "BH"),
    disease_genes = vapply(term_sets, function(s) paste(intersect(s, dis), collapse = "/"), character(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

#' Disease layer of [plot_pathway_network()]: gene set lookup + test
#' @return `list(table, disease_id, disease_name, N, K)`, `table` as
#'   [.pathway_net_disease_test()] plus logical `disease_enriched`.
#' @keywords internal
.pathway_net_disease <- function(proj, disease, sets, map, fdr) {
  dg <- patliRResults(proj, "disease_genes")
  if (is.null(dg) || nrow(dg) == 0) {
    cli::cli_abort(c("No {.val disease_genes} entry in {.arg proj}.",
                     "i" = "Run {.fn disease_genes_fetch} or {.fn disease_genes_import} first."))
  }
  if (!disease %in% dg$disease_id) {
    cli::cli_abort(c("{.arg disease} {.val {disease}} is not a {.field disease_id} of {.val disease_genes}.",
                     "i" = "Available: {.val {unique(dg$disease_id)}}."))
  }
  dg <- dg[dg$disease_id == disease, , drop = FALSE]
  dmap <- .pathway_bitr(unique(dg$uniprot_id), "UNIPROT", "ENTREZID")
  dis_entrez <- unique(dmap$ENTREZID)
  if ("gene_symbol" %in% names(dg)) {
    smap <- .pathway_bitr(unique(dg$gene_symbol[!dg$uniprot_id %in% dmap$UNIPROT]), "SYMBOL", "ENTREZID")
    dis_entrez <- unique(c(dis_entrez, smap$ENTREZID))
  }
  if (length(dis_entrez) == 0) cli::cli_abort("No gene of disease {.val {disease}} mapped to an Entrez ID.")
  tab <- .pathway_net_disease_test(sets, dis_entrez, unique(map$ENTREZID))
  tab$disease_enriched <- tab$disease_padj < fdr & tab$disease_k > 0
  name <- dg$disease_name[!is.na(dg$disease_name) & nzchar(dg$disease_name)][1]
  list(table = tab, disease_id = disease, disease_name = if (is.na(name)) disease else name,
       N = tab$universe_N[1], K = tab$disease_K[1], fdr = fdr)
}

#' Compound-term links: how many of a term's genes each compound hits
#'
#' @description
#' A link needs at least `min_hits` of the term's genes hit by the compound.
#' Compounds are ranked by the total number of (term, gene) hits over the
#' drawn terms (ties: more linked terms, then compound id) and the `top_n`
#' kept; each keeps only its `max_links` strongest links, by the share of
#' the term's genes it hits (`share = n_hit / term_n`; ties: more hits).
#' Without that cap a promiscuous compound links to nearly every term and
#' the figure becomes a hairball.
#' @param map `data.frame(compound_id, ENTREZID)`.
#' @return `data.frame(compound_id, ID, n_hit, term_n, share)`.
#' @keywords internal
.pathway_net_compound_links <- function(map, sets, min_hits = 2, top_n = 8, max_links = 5) {
  empty <- data.frame(compound_id = character(0), ID = character(0), n_hit = integer(0), term_n = integer(0),
                      share = numeric(0), stringsAsFactors = FALSE)
  map <- unique(map[, c("compound_id", "ENTREZID")])
  rows <- lapply(names(sets), function(id) {
    hit <- map[map$ENTREZID %in% sets[[id]], , drop = FALSE]
    if (nrow(hit) == 0) return(NULL)
    tab <- table(hit$compound_id)
    data.frame(compound_id = names(tab), ID = id, n_hit = as.integer(tab), term_n = length(sets[[id]]),
               stringsAsFactors = FALSE)
  })
  links <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(links)) return(empty)
  links <- links[links$n_hit >= min_hits, , drop = FALSE]
  if (nrow(links) == 0) return(empty)
  links$share <- links$n_hit / links$term_n
  n_terms <- tapply(links$ID, links$compound_id, length)
  n_hits <- tapply(links$n_hit, links$compound_id, sum)
  cmp <- names(n_hits)
  keep <- utils::head(cmp[order(-n_hits, -n_terms[cmp], cmp)], top_n)
  links <- links[links$compound_id %in% keep, , drop = FALSE]
  links <- do.call(rbind, lapply(keep, function(c) {
    l <- links[links$compound_id == c, , drop = FALSE]
    utils::head(l[order(-l$share, -l$n_hit, l$ID), , drop = FALSE], max_links)
  }))
  rownames(links) <- NULL
  links
}

#' Place the compounds on an ellipse around the term network
#'
#' @description
#' Each compound goes to the angle of the (hit-weighted) centre of the
#' terms it is linked to; angles closer than a minimum gap are spread
#' apart so compounds do not overlap.
#' @return `data.frame(compound_id, n_terms, x, y, angle)`.
#' @keywords internal
.pathway_net_compound_ring <- function(links, terms) {
  cx <- mean(range(terms$x))
  cy <- mean(range(terms$y))
  rx <- diff(range(terms$x)) / 2 + 2.4
  ry <- diff(range(terms$y)) / 2 + 2.0
  cmp <- unique(links$compound_id)
  ang <- vapply(cmp, function(c) {
    l <- links[links$compound_id == c, , drop = FALSE]
    i <- match(l$ID, terms$ID)
    mx <- stats::weighted.mean(terms$x[i], l$n_hit) - cx
    my <- stats::weighted.mean(terms$y[i], l$n_hit) - cy
    if (abs(mx) < 1e-9 && abs(my) < 1e-9) 0 else atan2(my / ry, mx / rx)
  }, numeric(1))
  ## evenly spaced slots, assigned in the circular order of the preferred
  ## angles (rotated to best match them) -- no two compounds collide
  ang <- .pathway_slot_angles(ang)
  data.frame(compound_id = cmp, n_terms = as.integer(table(links$compound_id)[cmp]),
             x = cx + rx * cos(ang), y = cy + ry * sin(ang), angle = ang, stringsAsFactors = FALSE)
}

#' Evenly spaced angular slots in the circular order of preferred angles
#'
#' @description
#' `n` angles `2 * pi / n` apart; the i-th smallest preferred angle gets
#' the i-th slot, and the whole set is rotated by the circular mean of the
#' differences so slots sit as close as possible to the preferred angles.
#' @return Numeric vector, same order as `ang`.
#' @keywords internal
.pathway_slot_angles <- function(ang) {
  n <- length(ang)
  if (n < 2) return(ang)
  o <- order(ang)
  step <- 2 * pi / n
  base <- (seq_len(n) - 1) * step
  d <- ang[o] - base
  offset <- atan2(mean(sin(d)), mean(cos(d)))
  out <- numeric(n)
  out[o] <- offset + base
  out
}

#' Clean per-term table returned as `attr(., "terms")`
#' @keywords internal
.pathway_net_term_table <- function(terms, cluster_info) {
  keep <- intersect(c("ID", "Description", "Count", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "neg_log10_padj",
                      "cluster", "degree", "x", "y", "disease_k", "term_n", "disease_expected", "disease_fold",
                      "disease_p", "disease_padj", "disease_enriched", "disease_genes", "geneID"), names(terms))
  out <- terms[, keep, drop = FALSE]
  out$cluster_name <- cluster_info$representative[match(out$cluster, cluster_info$cluster)]
  rownames(out) <- NULL
  out
}

#' Build the enrichment-map ggplot
#' @keywords internal
.pathway_net_ggplot <- function(terms, edges, cluster_info, links, cmp_nodes, dis, db, cond, similarity,
                                min_similarity, label_width, engine, width, seed) {
  pos <- stats::setNames(seq_len(nrow(terms)), terms$ID)
  terms$label <- vapply(terms$Description, function(d) {
    w <- strwrap(d, width = label_width)
    if (length(w) > 3) w <- c(w[1:2], .plot_truncate(paste(w[-(1:2)], collapse = " "), label_width))
    paste(w, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
  terms$tooltip <- sprintf("%s (%s)\n%d targets, adj. p = %.2g%s", terms$Description, terms$ID, terms$Count,
                           terms$p.adjust,
                           if (!is.null(dis)) sprintf("\n%d/%d %s genes, FDR = %.2g", terms$disease_k, terms$term_n,
                                                      dis$disease_name, terms$disease_padj) else "")
  terms$name <- terms$ID

  p <- ggplot2::ggplot()

  ## cluster outlines
  if (nrow(cluster_info) > 0) {
    in_cl <- !is.na(terms$cluster)
    hulls <- .network_view_hulls(terms$x[in_cl], terms$y[in_cl], paste0("C", terms$cluster[in_cl]), pad = 0.55)
    pal <- .network_view_qual_palette(cluster_info$id)[cluster_info$id]
    hulls$group <- factor(hulls$group, levels = cluster_info$id)
    p <- p +
      ggplot2::geom_polygon(data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$group),
                            alpha = 0.13, colour = NA) +
      ggplot2::geom_polygon(data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                            fill = NA, colour = unname(pal[as.character(hulls$group)]), linewidth = 0.45,
                            linetype = "22", show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = pal, breaks = cluster_info$id, labels = cluster_info$label,
                                 name = "Term clusters (Louvain)",
                                 guide = ggplot2::guide_legend(order = 4, override.aes = list(alpha = 0.35)))
  }

  ## compound spokes (under everything else)
  if (!is.null(cmp_nodes) && nrow(cmp_nodes) > 0) {
    ln <- links[links$compound_id %in% cmp_nodes$compound_id, , drop = FALSE]
    ln$x <- cmp_nodes$x[match(ln$compound_id, cmp_nodes$compound_id)]
    ln$y <- cmp_nodes$y[match(ln$compound_id, cmp_nodes$compound_id)]
    ln$xend <- terms$x[pos[ln$ID]]
    ln$yend <- terms$y[pos[ln$ID]]
    p <- p + ggplot2::geom_segment(data = ln, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                                   colour = "#0072B2", alpha = 0.22, linewidth = 0.3)
  }

  ## term-term overlap edges
  if (nrow(edges) > 0) {
    ed <- edges
    ed$x <- terms$x[pos[ed$from]]
    ed$y <- terms$y[pos[ed$from]]
    ed$xend <- terms$x[pos[ed$to]]
    ed$yend <- terms$y[pos[ed$to]]
    ed <- ed[order(ed$similarity), , drop = FALSE]
    p <- p +
      ggplot2::geom_segment(data = ed, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                                                    linewidth = .data$similarity),
                            colour = "grey55", alpha = 0.45, lineend = "round") +
      ggplot2::scale_linewidth(range = c(0.25, 2.4), limits = c(min_similarity, 1),
                               name = paste0("Gene overlap\n(", if (similarity == "jaccard") "Jaccard" else "overlap coef.", ")"),
                               guide = ggplot2::guide_legend(order = 3))
  }

  ## term nodes: white-bordered underlay + coloured disc
  p <- p +
    ggplot2::geom_point(data = terms, ggplot2::aes(x = .data$x, y = .data$y, size = .data$Count),
                        shape = 21, fill = "white", colour = "grey20", stroke = 0.5) +
    .network_view_points(engine, terms, list(size = "Count", colour = "neg_log10_padj"), shape = 19) +
    ggplot2::scale_size_area(max_size = 11.5, name = "Targets in term",
                             guide = ggplot2::guide_legend(order = 2, override.aes = list(colour = "grey45"))) +
    ggplot2::scale_colour_gradientn(
      colours = c("#FDE0A1", "#FDB462", "#F46D43", "#D53E4F", "#9E0142"),
      name = expression(-log[10] ~ "adj. p"),
      guide = ggplot2::guide_colourbar(order = 1, barheight = ggplot2::unit(2.6, "cm"))
    )

  ## disease: inner dot whose area is the disease share of the term's
  ## targets (same area scale as the node), ring when FDR-significant.
  ## Disease layers, compounds share the shape legend.
  shape_values <- c()
  if (!is.null(dis) && any(terms$disease_k > 0)) {
    inner <- terms[terms$disease_k > 0, , drop = FALSE]
    inner$inner_size <- inner$Count * inner$disease_k / inner$term_n
    inner$mark <- "share"
    p <- p + ggplot2::geom_point(data = inner, ggplot2::aes(x = .data$x, y = .data$y, size = .data$inner_size, shape = .data$mark),
                                 colour = "#1D3557", alpha = 0.9)
    shape_values <- c(shape_values, share = 19)
  }
  if (!is.null(dis) && any(terms$disease_enriched)) {
    ring <- terms[terms$disease_enriched, , drop = FALSE]
    ring$mark <- "ring"
    p <- p + ggplot2::geom_point(data = ring, ggplot2::aes(x = .data$x, y = .data$y, size = .data$Count, shape = .data$mark),
                                 fill = NA, colour = "grey5", stroke = 1.5)
    shape_values <- c(shape_values, ring = 21)
  }
  if (!is.null(cmp_nodes) && nrow(cmp_nodes) > 0) {
    cmp_nodes$mark <- "compound"
    cmp_nodes$name <- cmp_nodes$compound_id
    cmp_nodes$tooltip <- paste0(cmp_nodes$label, "\nlinked to ", cmp_nodes$n_terms, " term(s)")
    p <- p + .network_view_points(engine, cmp_nodes, list(shape = "mark"), size = 4.2, colour = "#0072B2")
    shape_values <- c(shape_values, compound = 18)
  }
  if (length(shape_values)) {
    dname <- if (!is.null(dis)) .plot_truncate(sub(" \\([^()]*\\)$", "", dis$disease_name), 30) else ""
    shape_labels <- c(
      share = paste0("Inner dot: share of the term's\ntargets that are ", dname, " genes"),
      ring = paste0("Ring: enriched in ", dname, " genes\n(hypergeometric, FDR < ", if (!is.null(dis)) dis$fdr else "", ")"),
      compound = paste0("Compound (top ", if (!is.null(cmp_nodes)) nrow(cmp_nodes) else 0, ")")
    )[names(shape_values)]
    p <- p + ggplot2::scale_shape_manual(
      values = shape_values, breaks = names(shape_values), labels = shape_labels, name = NULL,
      guide = ggplot2::guide_legend(order = 0, override.aes = list(
        size = c(share = 3, ring = 4.5, compound = 4)[names(shape_values)], fill = NA, stroke = 1.3,
        colour = c(share = "#1D3557", ring = "grey5", compound = "#0072B2")[names(shape_values)]))
    )
  }

  ## labels
  terms$fontface <- if (!is.null(dis)) ifelse(terms$disease_enriched, "bold", "plain") else "plain"
  ## ggrepel keeps labels off the discs when it knows their drawn size
  ## (scale_size_area(): size = max_size * sqrt(Count / max(Count)))
  terms$point_size <- 11.5 * sqrt(terms$Count / max(terms$Count)) + 1
  label_aes <- if (.plot_use_repel()) {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, fontface = .data$fontface, point.size = .data$point_size)
  } else {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, fontface = .data$fontface)
  }
  p <- p + .plot_text_layer(
    data = terms, mapping = label_aes,
    size = 2.45, colour = "grey10", lineheight = 0.9, inherit.aes = FALSE,
    repel_args = list(box.padding = 0.35, point.padding = 0.25, min.segment.length = 0.3, segment.colour = "grey45",
                      segment.size = 0.25, bg.colour = "white", bg.r = 0.12, max.overlaps = Inf, seed = seed,
                      force = 3, max.time = 3, max.iter = 20000),
    text_args = list(vjust = -1.2)
  )
  if (nrow(cluster_info) > 0) {
    tags <- do.call(rbind, lapply(seq_len(nrow(cluster_info)), function(i) {
      k <- cluster_info$cluster[i]
      mem <- terms[!is.na(terms$cluster) & terms$cluster == k, , drop = FALSE]
      data.frame(id = cluster_info$id[i], x = mean(range(mem$x)), y = max(mem$y) + 0.9, stringsAsFactors = FALSE)
    }))
    tag_cols <- .network_view_qual_palette(cluster_info$id)[tags$id]
    p <- p + ggplot2::geom_label(data = tags, ggplot2::aes(x = .data$x, y = .data$y, label = .data$id),
                                 size = 3.1, fontface = "bold", colour = unname(tag_cols), fill = "white",
                                 label.padding = ggplot2::unit(0.12, "lines"), linewidth = 0.3,
                                 inherit.aes = FALSE)
  }
  if (!is.null(cmp_nodes) && nrow(cmp_nodes) > 0) {
    cmp_nodes$hjust <- ifelse(abs(cos(cmp_nodes$angle)) < 0.3, 0.5, ifelse(cos(cmp_nodes$angle) > 0, 0, 1))
    cmp_nodes$lx <- cmp_nodes$x + 0.4 * cos(cmp_nodes$angle)
    cmp_nodes$ly <- cmp_nodes$y + 0.45 * sin(cmp_nodes$angle)
    p <- p + ggplot2::geom_text(data = cmp_nodes, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$label, hjust = .data$hjust),
                                size = 2.5, colour = "#0B4F7C", fontface = "italic", inherit.aes = FALSE)
  }

  db_name <- c(kegg = "KEGG", reactome = "Reactome", go = "GO")[[db]]
  sub <- sprintf(paste0(
    "The %d most significant %s terms of %s (adjusted p <= cutoff). Node area = number of the extract's targets in the term; ",
    "colour = -log10 adjusted p; an edge joins two terms whose target genes overlap (%s >= %.2f). ",
    "%d term-term edge%s; %d cluster%s of related terms outlined (Louvain)."
  ), nrow(terms), db_name, cond, if (similarity == "jaccard") "Jaccard" else "overlap coefficient", min_similarity,
  nrow(edges), if (nrow(edges) == 1) "" else "s", nrow(cluster_info), if (nrow(cluster_info) == 1) "" else "s")
  if (!is.null(dis)) {
    sub <- paste0(sub, sprintf(
      " Dark inner dot: share of the term's targets that are %s genes; ring: that share is higher than chance (hypergeometric, universe = the extract's %d mapped targets, %d of them disease genes; BH FDR < %s): %d of %d terms.",
      dis$disease_name, dis$N, dis$K, dis$fdr, sum(terms$disease_enriched), nrow(terms)))
  }
  if (!is.null(cmp_nodes) && nrow(cmp_nodes) > 0) {
    sub <- paste0(sub, sprintf(paste0(
      " Blue diamonds: the %d compounds with the most target hits in these terms, each linked to the (at most %d) terms ",
      "of which it hits the largest share of genes (>= %d genes)."),
      nrow(cmp_nodes), max(table(links$compound_id)), min(links$n_hit)))
  }
  xr <- range(c(terms$x, if (!is.null(cmp_nodes)) cmp_nodes$x))
  yr <- range(c(terms$y, if (!is.null(cmp_nodes)) cmp_nodes$y))
  padx <- if (!is.null(cmp_nodes) && nrow(cmp_nodes) > 0) 3.4 else 1.2

  p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = padx)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.9, 1.4))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = .plot_wrap(paste0(db_name, " pathway network -- ", cond,
                                if (!is.null(dis)) paste0(" -- ", dis$disease_name) else ""), .plot_wrap_width(width, 12)),
      subtitle = .plot_wrap(sub, .plot_wrap_width(width, 8))
    ) +
    .network_view_theme("right") +
    ggplot2::theme(legend.spacing.y = ggplot2::unit(0.2, "cm"), legend.key.spacing.y = ggplot2::unit(0.12, "cm"))
}

# ---------------------------------------------------------------------------
# plot_kegg_topology()
# ---------------------------------------------------------------------------

#' KEGG pathway topology: the directed gene relations of a pathway, with
#' the extract's targets highlighted
#'
#' @description
#' Draws, for one or more KEGG pathways, the gene-gene relations parsed from
#' KEGG's KGML by [network_kegg_topology()] (`PPrel` protein-protein,
#' `GErel` gene-expression, `ECrel` enzyme-enzyme): one box per gene, one
#' directed edge per relation. Boxes of genes the extract's compounds are
#' predicted to hit (`network_edges`) are filled by the highest probability
#' among those compounds and drawn larger the more compounds hit them;
#' other genes stay white. Edges follow KEGG's notation: colour and end
#' marker give the effect (arrow = activation / expression, flat bar =
#' inhibition / repression, no end = binding), line type the mechanism
#' (dashed = phosphorylation, dotted = indirect, ...).
#'
#' When `pathway_id` is `NULL`, the `top_n_pathways` pathways with the most
#' target genes of the extract among their relations are drawn, one panel
#' each.
#'
#' @section Paralog boxes:
#' A KGML box usually bundles several paralogs (e.g. `MAPK1`, `MAPK3`) that
#' take part in exactly the same relations; [network_kegg_topology()]
#' expands such a box into one row per gene pair. With `collapse_paralogs =
#' TRUE` (default) genes with identical incoming and outgoing relations are
#' merged back into one box, labelled `MAPK1/3` (hit genes listed first) --
#' close to how KEGG draws the pathway, and much less cluttered.
#'
#' @section Layout:
#' `"sugiyama"` puts upstream genes on the left and downstream genes on the
#' right, in layers (Sugiyama et al. 1981; [igraph::layout_with_sugiyama()]),
#' which reads like a signalling cascade. `"auto"` (default) uses it for
#' every connected component whose relations are close to acyclic (a
#' feedback arc set of at most 20% of its edges), and Fruchterman-Reingold
#' otherwise. Components are packed side by side, largest first.
#'
#' @section Full pathway vs. extract-only topology:
#' [network_kegg_topology()] keeps, by default (`restrict_to_network =
#' TRUE`), only relations whose two genes are both targets of the extract --
#' every box is then filled. Run it with `restrict_to_network = FALSE` to
#' store the whole pathway, and the figure shows the extract's targets
#' within the full cascade (unfilled boxes = genes it does not hit).
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param condition Character scalar with rows in `network_kegg_topology`.
#'   `NULL` (default) works when exactly one condition has topology.
#' @param pathway_id `NULL` (default, see above) or a character vector of
#'   KEGG pathway IDs present in `network_kegg_topology` for `condition`.
#' @param top_n_pathways Integer, default `6`. Ignored when `pathway_id` is
#'   given.
#' @param layout `"auto"` (default), `"sugiyama"`, `"fr"` or `"kk"`.
#' @param collapse_paralogs Logical, default `TRUE` (see the section above).
#' @param label_max_chars Integer, default `18`: longer box labels are cut.
#' @param ncol `NULL` (default, automatic) or the number of panel columns.
#' @param seed Integer seed for the force-directed layouts.
#' @param width,height Figure size in inches; `NULL` (default) picks a size
#'   from the number of panels.
#' @param engine `"static"` (default) or `"ggiraph"`.
#'
#' @return A `ggplot` (or `girafe`) object with attributes
#'   `pathway_summary` (one row per drawn pathway: `pathway_id`,
#'   `pathway_title`, `n_genes`, `n_hit`, `share_hit`, `n_relations`,
#'   `n_compounds`, `n_boxes`, `layout`), `nodes` and `edges`. With
#'   `save = TRUE` the PNG is written and logged in
#'   `patliRResults(proj, "kegg_topology_plot_log")` (keyed by `condition`,
#'   `pathways`), and the updated project is in `attr(., "proj")`.
#'
#' @references
#' Kanehisa, M. & Goto, S. (2000), "KEGG: Kyoto Encyclopedia of Genes and
#' Genomes", *Nucleic Acids Research* 28(1), 27-30.
#' \doi{10.1093/nar/28.1.27}
#'
#' Sugiyama, K., Tagawa, S. & Toda, M. (1981), "Methods for visual
#' understanding of hierarchical system structures", *IEEE Transactions on
#' Systems, Man, and Cybernetics* 11(2), 109-125.
#' \doi{10.1109/TSMC.1981.4308636}
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... network_build(), network_enrich(db = "kegg"), then the full topology:
#' proj <- network_kegg_topology(proj, condition = "FLO-ET",
#'                               pathway_ids = c("hsa04010", "hsa04020"),
#'                               restrict_to_network = FALSE) # needs internet
#' plot_kegg_topology(proj, condition = "FLO-ET", pathway_id = "hsa04010")
#' plot_kegg_topology(proj, condition = "FLO-ET") # the most-hit pathways
#' }
#'
#' @seealso [network_kegg_topology()], [network_pathview()],
#'   [plot_pathway_network()]
#' @export
plot_kegg_topology <- function(proj, condition = NULL, pathway_id = NULL, top_n_pathways = 6,
                               layout = c("auto", "sugiyama", "fr", "kk"), collapse_paralogs = TRUE,
                               label_max_chars = 18, ncol = NULL, seed = 1,
                               engine = c("static", "ggiraph"),
                               save = TRUE, out_dir = NULL, width = NULL, height = NULL, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  layout <- match.arg(layout)
  engine <- match.arg(engine)
  .pathway_check_count(top_n_pathways, "top_n_pathways", min = 1)
  .pathway_check_flag(collapse_paralogs, "collapse_paralogs")
  .pathway_check_count(label_max_chars, "label_max_chars", min = 6)
  .pathway_check_count(seed, "seed", min = -Inf)
  if (!is.null(ncol)) .pathway_check_count(ncol, "ncol", min = 1)
  if (!is.null(pathway_id) && (!is.character(pathway_id) || length(pathway_id) == 0 || anyNA(pathway_id))) {
    cli::cli_abort("{.arg pathway_id} must be NULL or a character vector of KEGG pathway IDs.")
  }
  for (a in c("width", "height")) {
    v <- get(a)
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || is.na(v) || v <= 0)) {
      cli::cli_abort("{.arg {a}} must be NULL or a single positive number.")
    }
  }
  engine <- .plot_require(engine)

  topo_all <- patliRResults(proj, "network_kegg_topology")
  if (is.null(topo_all) || nrow(topo_all) == 0) {
    cli::cli_abort(c("No {.val network_kegg_topology} entry in {.arg proj}.",
                     "i" = "Run {.fn network_kegg_topology} first."))
  }
  cond <- .pathway_single_condition(unique(topo_all$condition), condition, "plot_kegg_topology", "network_kegg_topology")
  topo <- topo_all[topo_all$condition == cond & !is.na(topo_all$from_kegg) & !is.na(topo_all$to_kegg), , drop = FALSE]

  edges_all <- patliRResults(proj, "network_edges")
  ct <- if (is.null(edges_all)) NULL else edges_all[edges_all$condition == cond, , drop = FALSE]
  genes <- .kegg_topology_gene_stats(topo, ct)
  summary_all <- .kegg_topology_pathway_summary(topo, genes)
  pids <- .kegg_topology_select(summary_all, pathway_id, top_n_pathways, cond)
  summ <- summary_all[match(pids, summary_all$pathway_id), , drop = FALSE]

  genes$symbol <- .kegg_topology_symbols(proj, cond, genes)

  restore_rng <- .with_seed(seed)
  on.exit(restore_rng(), add = TRUE)
  prepared <- lapply(pids, function(pid) {
    .kegg_topology_prepare(topo[topo$pathway_id == pid, , drop = FALSE], genes, collapse_paralogs = collapse_paralogs,
                           layout = layout, label_max_chars = label_max_chars)
  })
  restore_rng()

  ## figure / panel size: fixed when given, otherwise grown with the
  ## layout's extent (abstract units: one layer ~ 1 in wide, one stacked
  ## box ~ 0.34 in tall) within sensible bounds
  n_panels <- length(pids)
  if (is.null(ncol)) ncol <- if (n_panels <= 1) 1 else if (n_panels <= 4) 2 else 3
  ncol <- min(ncol, n_panels)
  nrow_p <- ceiling(n_panels / ncol)
  ext <- vapply(prepared, function(b) c(diff(range(b$xy[, 1])), diff(range(b$xy[, 2]))), numeric(2))
  pad_w <- max(vapply(prepared, function(b) 2 * max(b$nodes$hw), numeric(1))) + 0.2
  pad_h <- max(vapply(prepared, function(b) 2 * max(b$nodes$hh), numeric(1))) + 0.2
  want_w <- min(max(ext[1, ]) * 0.42 + pad_w, if (n_panels == 1) 15 else 6.5)
  want_h <- min(max(ext[2, ]) * 0.3 + pad_h, if (n_panels == 1) 13 else 5)
  strip_h <- if (n_panels > 1) 0.5 else 0
  if (is.null(width)) width <- max(if (n_panels == 1) 8.5 else 5, ncol * max(want_w, 3.8) + 3.1)
  if (is.null(height)) height <- max(5.2, nrow_p * (max(want_h, 3.2) + strip_h) + 2.2)
  panel_w <- max(2, (width - 3.1) / ncol - 0.15)
  panel_h <- max(1.6, (height - 2.2) / nrow_p - strip_h - 0.1)

  built <- lapply(prepared, .kegg_topology_finish_panel, panel_w = panel_w, panel_h = panel_h)
  nodes <- do.call(rbind, lapply(built, `[[`, "nodes"))
  edges <- do.call(rbind, lapply(built, `[[`, "edges"))
  paths <- do.call(rbind, lapply(built, `[[`, "paths"))
  tees <- do.call(rbind, lapply(built, `[[`, "tees"))
  summ$n_boxes <- vapply(built, function(b) nrow(b$nodes), integer(1))
  summ$layout <- vapply(built, `[[`, character(1), "layout")
  rownames(summ) <- NULL

  restricted <- "restrict_to_network" %in% names(topo) && nrow(topo) > 0 &&
    isTRUE(all(as.logical(topo$restrict_to_network[topo$pathway_id %in% pids])))
  p <- .kegg_topology_ggplot(nodes, edges, paths, tees, summ, cond, ncol, panel_w, panel_h, engine, width,
                             restricted = restricted)
  attr(p, "pathway_summary") <- summ
  attr(p, "nodes") <- nodes
  attr(p, "edges") <- edges

  pathways_key <- paste(pids, collapse = "+")
  file_tag <- if (length(pids) <= 3) paste(pids, collapse = "+") else paste0("top", length(pids), "_", substr(rlang::hash(pids), 1, 8))
  log_row <- data.frame(
    condition = cond, pathways = pathways_key, n_pathways = length(pids),
    n_genes = sum(summ$n_genes), n_hit = sum(summ$n_hit), layout = layout,
    collapse_paralogs = collapse_paralogs, path = NA_character_, stringsAsFactors = FALSE
  )
  .plot_finish(
    proj, p,
    name = "kegg_topology_plot_log", filename = paste0("kegg_topology_", cond, "_", file_tag, ".png"),
    log_row = log_row, key_cols = c("condition", "pathways"),
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Per-KEGG-gene target statistics
#' @param topo `network_kegg_topology` rows of one condition.
#' @param ct `network_edges` rows of the same condition (or `NULL`).
#' @return `data.frame(kegg, uniprots, is_target, max_weight, n_compounds,
#'   compounds)` -- one row per KEGG gene of `topo`, sorted; `uniprots` and
#'   `compounds` are `;`-joined.
#' @keywords internal
.kegg_topology_gene_stats <- function(topo, ct) {
  gu <- unique(rbind(
    data.frame(kegg = topo$from_kegg, uniprot = topo$from_uniprot, stringsAsFactors = FALSE),
    data.frame(kegg = topo$to_kegg, uniprot = topo$to_uniprot, stringsAsFactors = FALSE)
  ))
  gu <- gu[!is.na(gu$kegg), , drop = FALSE]
  genes <- sort(unique(gu$kegg))
  out <- data.frame(kegg = genes, uniprots = "", is_target = FALSE, max_weight = NA_real_,
                    n_compounds = 0L, compounds = "", stringsAsFactors = FALSE)
  if (length(genes) == 0) return(out)
  ups <- tapply(gu$uniprot, gu$kegg, function(u) paste(sort(unique(stats::na.omit(u))), collapse = ";"))
  out$uniprots <- unname(ups[genes])
  if (is.null(ct) || nrow(ct) == 0) return(out)
  hit <- merge(gu[!is.na(gu$uniprot), , drop = FALSE],
               ct[, c("compound_id", "uniprot_id", "weight")], by.x = "uniprot", by.y = "uniprot_id")
  if (nrow(hit) == 0) return(out)
  i <- match(sort(unique(hit$kegg)), genes)
  w <- tapply(hit$weight, hit$kegg, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  cmp <- tapply(hit$compound_id, hit$kegg, function(v) sort(unique(v)), simplify = FALSE)
  out$is_target[i] <- TRUE
  out$max_weight[i] <- unname(w[genes[i]])
  out$n_compounds[i] <- unname(lengths(cmp[genes[i]]))
  out$compounds[i] <- vapply(cmp[genes[i]], paste, character(1), collapse = ";")
  out
}

#' Pathway-level summary of the topology table
#' @return `data.frame(pathway_id, pathway_title, n_genes, n_hit, share_hit,
#'   n_relations, n_compounds)`, most-hit pathways first (ties: larger
#'   share, then more relations, then ID).
#' @keywords internal
.kegg_topology_pathway_summary <- function(topo, genes) {
  pids <- unique(topo$pathway_id)
  by_pid <- split(seq_len(nrow(topo)), topo$pathway_id)
  rows <- lapply(pids, function(pid) {
    t <- topo[by_pid[[pid]], , drop = FALSE]
    g <- unique(c(t$from_kegg, t$to_kegg))
    gi <- genes[match(g, genes$kegg), , drop = FALSE]
    pairs <- unique(paste(t$from_kegg[t$from_kegg != t$to_kegg], t$to_kegg[t$from_kegg != t$to_kegg]))
    cmp <- unique(unlist(strsplit(gi$compounds[gi$is_target], ";", fixed = TRUE)))
    data.frame(pathway_id = pid, pathway_title = t$pathway_title[1], n_genes = length(g),
               n_hit = sum(gi$is_target), share_hit = sum(gi$is_target) / length(g),
               n_relations = length(pairs), n_compounds = length(cmp), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame(pathway_id = character(0), pathway_title = character(0), n_genes = integer(0),
                                      n_hit = integer(0), share_hit = numeric(0), n_relations = integer(0),
                                      n_compounds = integer(0), stringsAsFactors = FALSE))
  out <- out[order(-out$n_hit, -out$share_hit, -out$n_relations, out$pathway_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Which pathways to draw
#' @return Character vector of pathway IDs.
#' @keywords internal
.kegg_topology_select <- function(summary, pathway_id, top_n, cond) {
  if (!is.null(pathway_id)) {
    pathway_id <- unique(pathway_id)
    unknown <- setdiff(pathway_id, summary$pathway_id)
    if (length(unknown) > 0) {
      hint <- sprintf("network_kegg_topology(proj, condition = \"%s\", pathway_ids = \"%s\", restrict_to_network = FALSE)",
                      cond, unknown[1])
      cli::cli_abort(c(
        "No relation in {.val network_kegg_topology} for condition {.val {cond}} and pathway{?s} {.val {unknown}}.",
        "i" = "With the default {.code restrict_to_network = TRUE}, {.fn network_kegg_topology} keeps only relations between two of the extract's targets; a pathway can end up empty.",
        "i" = "Fetch the whole pathway with {.code {hint}}."
      ))
    }
    return(pathway_id)
  }
  hit <- summary[summary$n_hit > 0, , drop = FALSE]
  if (nrow(hit) == 0) hit <- summary
  utils::head(hit$pathway_id, top_n)
}

#' Map a relation subtype (KGML `<subtype name>`) to an effect and a
#' mechanism
#'
#' @description
#' `effect` (edge colour + end marker): `activation`, `expression` (arrow);
#' `inhibition`, `repression` (flat bar); `indirect effect`, `modification`
#' (a phosphorylation/... with no stated sign; arrow); `binding/association`
#' and `other` (dissociation, state change, missing interaction, compound,
#' none; no end marker).
#' `mode` (line type): `phosphorylation`, `dephosphorylation`,
#' `other modification` (glycosylation, ubiquitination, methylation),
#' `via compound`, `indirect`, `direct`.
#' @return `data.frame(subtype, effect, mode)`, one row per input.
#' @keywords internal
.kegg_topology_subtype_map <- function(subtype) {
  s <- tolower(trimws(as.character(subtype)))
  s[is.na(s)] <- ""
  effect <- ifelse(s == "activation", "activation",
            ifelse(s == "expression", "expression",
            ifelse(s == "inhibition", "inhibition",
            ifelse(s == "repression", "repression",
            ifelse(s == "indirect effect", "indirect effect",
            ifelse(s %in% c("phosphorylation", "dephosphorylation", "glycosylation", "ubiquitination", "methylation"), "modification",
            ifelse(s %in% c("binding/association", "binding", "association"), "binding/association", "other")))))))
  mode <- ifelse(s == "phosphorylation", "phosphorylation",
          ifelse(s == "dephosphorylation", "dephosphorylation",
          ifelse(s %in% c("glycosylation", "ubiquitination", "methylation"), "other modification",
          ifelse(s == "compound", "via compound",
          ifelse(s == "indirect effect", "indirect", "direct")))))
  data.frame(subtype = subtype, effect = effect, mode = mode, stringsAsFactors = FALSE)
}

#' @rdname dot-kegg_topology_subtype_map
#' @keywords internal
.kegg_topology_effect_levels <- function() {
  c("activation", "expression", "inhibition", "repression", "indirect effect", "modification",
    "binding/association", "other")
}

#' @rdname dot-kegg_topology_subtype_map
#' @keywords internal
.kegg_topology_mode_levels <- function() {
  c("phosphorylation", "dephosphorylation", "other modification", "via compound", "indirect", "direct")
}

#' One edge per directed gene pair, with the dominant effect and mechanism
#'
#' @description
#' A gene pair can carry several subtypes (e.g. `activation` +
#' `phosphorylation`). The drawn effect is the first present in
#' `inhibition > repression > activation > expression > indirect effect >
#' modification > binding/association > other`, the mechanism the first in
#' `phosphorylation > dephosphorylation > other modification > via compound
#' > indirect > direct`; all subtypes are kept in `subtypes`, and the
#' pathways the pair appears in in `pathways`. Self-relations are dropped.
#' @return `data.frame(from, to, effect, mode, subtypes, relation_types,
#'   pathways)`, in order of first appearance.
#' @keywords internal
.kegg_topology_collapse_relations <- function(topo) {
  t <- topo[!is.na(topo$from_kegg) & !is.na(topo$to_kegg) & topo$from_kegg != topo$to_kegg, , drop = FALSE]
  if (nrow(t) == 0) {
    return(data.frame(from = character(0), to = character(0), effect = character(0), mode = character(0),
                      subtypes = character(0), relation_types = character(0), pathways = character(0),
                      stringsAsFactors = FALSE))
  }
  sm <- .kegg_topology_subtype_map(t$relation_subtype)
  eff_lev <- c("inhibition", "repression", "activation", "expression", "indirect effect", "modification",
               "binding/association", "other")
  mode_lev <- .kegg_topology_mode_levels()
  key <- factor(paste(t$from_kegg, t$to_kegg, sep = "\r"), levels = unique(paste(t$from_kegg, t$to_kegg, sep = "\r")))
  first <- !duplicated(key)
  eff <- tapply(match(sm$effect, eff_lev), key, min)
  mod <- tapply(match(sm$mode, mode_lev), key, min)
  collapse <- function(v) tapply(v, key, function(x) {
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    if (length(x)) paste(x, collapse = "; ") else "none"
  })
  out <- data.frame(
    from = t$from_kegg[first], to = t$to_kegg[first],
    effect = eff_lev[eff], mode = mode_lev[mod],
    subtypes = unname(collapse(t$relation_subtype)),
    relation_types = unname(collapse(t$relation_type)),
    pathways = if ("pathway_id" %in% names(t)) unname(collapse(t$pathway_id)) else NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Merge genes that take part in exactly the same relations (paralog boxes)
#' @param edges Output of [.kegg_topology_collapse_relations()].
#' @return Named character vector gene -> box id (the box id is the first
#'   member gene in sorted order).
#' @keywords internal
.kegg_topology_paralog_groups <- function(edges) {
  genes <- sort(unique(c(edges$from, edges$to)))
  if (length(genes) == 0) return(stats::setNames(character(0), character(0)))
  sig <- vapply(genes, function(g) {
    o <- edges[edges$from == g, , drop = FALSE]
    i <- edges[edges$to == g, , drop = FALSE]
    paste(c(sort(paste0(">", o$to, "|", o$effect, "|", o$mode)), sort(paste0("<", i$from, "|", i$effect, "|", i$mode))),
          collapse = ",")
  }, character(1))
  first <- tapply(genes, sig, function(x) sort(x)[1])
  stats::setNames(unname(first[sig]), genes)
}

#' Compact label of a paralog box: common prefix + suffixes, hit genes first
#' ("MAPK1/3", "PIK3CA/CB/CD +1")
#' @keywords internal
.kegg_topology_box_label <- function(symbols, hit = rep(FALSE, length(symbols)), max_members = 4, max_chars = 18) {
  o <- order(!hit, nchar(symbols), symbols)
  s <- unique(symbols[o])
  if (length(s) == 1) return(.plot_truncate(s, max_chars))
  shown <- utils::head(s, max_members)
  extra <- length(s) - length(shown)
  pre <- shown[1]
  for (x in shown[-1]) {
    n <- 0
    while (n < min(nchar(pre), nchar(x)) && substr(pre, n + 1, n + 1) == substr(x, n + 1, n + 1)) n <- n + 1
    pre <- substr(pre, 1, n)
  }
  ## keep the prefix at a letter-digit boundary-ish length and never empty the first symbol
  if (nchar(pre) >= 3 && all(nchar(shown) > nchar(pre))) {
    lab <- paste0(shown[1], "/", paste(substring(shown[-1], nchar(pre) + 1), collapse = "/"))
  } else {
    lab <- paste(shown, collapse = "/")
  }
  if (nchar(lab) > max_chars) {
    lab <- shown[1]
    extra <- length(s) - 1
  }
  if (extra > 0) lab <- paste0(lab, " +", extra)
  lab
}

#' Gene symbols for KEGG gene IDs (Entrez for `hsa:`), falling back to the
#' UniProt label, then to the bare ID
#' @keywords internal
.kegg_topology_symbols <- function(proj, cond, genes) {
  if (nrow(genes) == 0) return(character(0))
  entrez <- sub("^[a-z]+:", "", genes$kegg)
  sym <- rep(NA_character_, nrow(genes))
  human <- grepl("^hsa:", genes$kegg)
  if (any(human)) {
    map <- .pathway_bitr(entrez[human], "ENTREZID", "SYMBOL")
    if (nrow(map) > 0) sym[human] <- map$SYMBOL[match(entrez[human], map$ENTREZID)]
  }
  miss <- is.na(sym)
  if (any(miss)) {
    first_up <- sub(";.*$", "", genes$uniprots[miss])
    has_up <- nzchar(first_up)
    if (any(has_up)) {
      lab <- .plot_label_nodes(proj, cond, first_up[has_up], "target")
      lab[lab == first_up[has_up]] <- NA_character_
      sym[which(miss)[has_up]] <- lab
    }
  }
  ifelse(is.na(sym), entrez, sym)
}

#' Box geometry (inches) for a label drawn by `geom_label()` at `size` mm
#' @return `list(hw, hh)` half-width and half-height in inches.
#' @keywords internal
.kegg_topology_box_size <- function(label, size) {
  fs <- size * ggplot2::.pt / 72
  pad <- 0.25 * fs * 1.2
  list(hw = nchar(label) * 0.56 * fs / 2 + pad, hh = fs * 0.62 + pad)
}

#' Boxes, box-level edges and abstract layout of one pathway panel
#' @return `list(pathway_id, nodes, e, xy, layout)`; `xy` in abstract
#'   layout units (see [.kegg_topology_layout()]).
#' @keywords internal
.kegg_topology_prepare <- function(t, genes, collapse_paralogs, layout, label_max_chars) {
  pid <- t$pathway_id[1]
  rel <- .kegg_topology_collapse_relations(t)
  gene_ids <- sort(unique(c(t$from_kegg, t$to_kegg)))
  group <- if (collapse_paralogs && nrow(rel) > 0) .kegg_topology_paralog_groups(rel) else stats::setNames(gene_ids, gene_ids)
  missing <- setdiff(gene_ids, names(group))
  group <- c(group, stats::setNames(missing, missing))

  gi <- genes[match(names(group), genes$kegg), , drop = FALSE]
  gi$box <- unname(group)
  boxes <- unique(gi$box)
  nodes <- do.call(rbind, lapply(boxes, function(b) {
    m <- gi[gi$box == b, , drop = FALSE]
    w <- suppressWarnings(max(m$max_weight, na.rm = TRUE))
    data.frame(
      pathway_id = pid, name = paste0(pid, "::", b), box = b, members = paste(m$kegg, collapse = ";"),
      symbols = paste(m$symbol[order(!m$is_target, m$symbol)], collapse = "/"),
      n_members = nrow(m), n_hit_members = sum(m$is_target), is_target = any(m$is_target),
      max_weight = if (is.finite(w)) w else NA_real_,
      n_compounds = max(m$n_compounds),
      label = .kegg_topology_box_label(m$symbol, m$is_target, max_chars = label_max_chars),
      stringsAsFactors = FALSE
    )
  }))

  e <- rel
  e$from_box <- unname(group[e$from])
  e$to_box <- unname(group[e$to])
  e <- e[e$from_box != e$to_box, , drop = FALSE]
  e <- e[!duplicated(paste(e$from_box, e$to_box, e$effect, e$mode)), , drop = FALSE]
  ## keep one edge per box pair (the dominant effect) so parallel edges do not stack
  eff_rank <- stats::setNames(c(3, 4, 1, 2, 5, 6, 7, 8), .kegg_topology_effect_levels())
  e <- e[order(paste(e$from_box, e$to_box), eff_rank[e$effect]), , drop = FALSE]
  e <- e[!duplicated(paste(e$from_box, e$to_box)), , drop = FALSE]

  ## label size (mm): hit genes grow with the number of compounds hitting them
  nodes$size <- .kegg_topology_label_size(nodes$n_compounds, nodes$is_target)
  bs <- .kegg_topology_box_size(nodes$label, nodes$size)
  nodes$hw <- bs$hw
  nodes$hh <- bs$hh

  lay <- .kegg_topology_layout(nodes$box, e[, c("from_box", "to_box")], layout)
  list(pathway_id = pid, nodes = nodes, e = e, xy = lay$xy, bends = lay$bends, layout = lay$used)
}

#' Fit a prepared panel into a `panel_w` x `panel_h` inch panel and compute
#' the edge geometry
#' @param b Output of [.kegg_topology_prepare()].
#' @return `list(nodes, edges, paths, tees, layout)`, each table with a
#'   `pathway_id` column.
#' @keywords internal
.kegg_topology_finish_panel <- function(b, panel_w, panel_h) {
  nodes <- b$nodes
  pad_x <- max(nodes$hw) + 0.05
  pad_y <- max(nodes$hh) + 0.05
  fit <- .kegg_topology_fit(b$xy, panel_w, panel_h, pad_x, pad_y)
  nodes$x <- fit$xy[match(nodes$box, rownames(fit$xy)), 1]
  nodes$y <- fit$xy[match(nodes$box, rownames(fit$xy)), 2]
  bends <- b$bends
  if (!is.null(bends) && nrow(bends)) {
    m <- fit$tf(as.matrix(bends[, c("x", "y")]))
    bends$x <- m[, 1]
    bends$y <- m[, 2]
  }
  eg <- .kegg_topology_edge_geometry(b$e, nodes, bends)
  for (nm in c("edges", "paths", "tees")) eg[[nm]]$pathway_id <- rep(b$pathway_id, nrow(eg[[nm]]))
  list(nodes = nodes[, setdiff(names(nodes), c("hw", "hh"))], edges = eg$edges, paths = eg$paths, tees = eg$tees,
       layout = b$layout)
}

#' Label size (mm) of a box from its number of hitting compounds
#' @keywords internal
.kegg_topology_label_size <- function(n_compounds, is_target) {
  is_target <- rep_len(is_target, length(n_compounds))
  ifelse(is_target, 2.45 + 1.0 * pmin(1, sqrt(pmax(n_compounds - 1, 0)) / sqrt(24)), 2.25)
}

#' Layout of one pathway's box graph in abstract units
#'
#' @description
#' Components are laid out separately (Sugiyama left-to-right, or
#' FR/KK) and packed in rows, largest first. For a layered component, an
#' edge spanning several layers is routed through the dummy vertices of the
#' Sugiyama layout (returned as `bends`), which keeps long edges out of the
#' boxes they would otherwise cross.
#' @param e Two-column data frame of box-level edges (from, to).
#' @return `list(xy = matrix (rownames = boxes), bends = data.frame(from,
#'   to, k, x, y), used = "sugiyama"/"fr"/"kk"/"mixed"/"none")`.
#' @keywords internal
.kegg_topology_layout <- function(boxes, e, layout = "auto") {
  el <- data.frame(from = as.character(e[[1]]), to = as.character(e[[2]]), stringsAsFactors = FALSE)
  g <- igraph::graph_from_data_frame(el, directed = TRUE,
                                     vertices = data.frame(name = boxes, stringsAsFactors = FALSE))
  comp <- igraph::components(g, mode = "weak")
  parts <- list()
  part_bends <- list()
  used <- character(0)
  for (k in order(-comp$csize, seq_along(comp$csize))) {
    vs <- igraph::V(g)$name[comp$membership == k]
    h <- igraph::induced_subgraph(g, vs)
    vs <- igraph::V(h)$name
    n <- igraph::vcount(h)
    bends <- NULL
    if (n == 1) {
      xy <- matrix(0, 1, 2)
      u <- "single"
    } else {
      use_sug <- layout == "sugiyama" || (layout == "auto" && .kegg_topology_is_layered(h))
      if (use_sug) {
        s <- igraph::layout_with_sugiyama(h, hgap = 1, vgap = 1, maxiter = 200)
        L <- rbind(s$layout, s$layout.dummy)
        ## igraph puts the first layer (sources) at the largest y: map
        ## layers to x so the cascade reads left (upstream) to right
        XY <- cbind(-L[, 2] * 2.6, -L[, 1])
        xy <- XY[seq_len(n), , drop = FALSE]
        bends <- .kegg_topology_sugiyama_bends(h, s$extd_graph, XY)
        u <- "sugiyama"
      } else if (layout == "kk") {
        xy <- .kegg_topology_unit_scale(igraph::layout_with_kk(igraph::as_undirected(h, mode = "collapse")), n)
        u <- "kk"
      } else {
        xy <- .kegg_topology_unit_scale(igraph::layout_with_fr(igraph::as_undirected(h, mode = "collapse"), niter = 1500), n)
        u <- "fr"
      }
    }
    rownames(xy) <- vs
    parts[[length(parts) + 1]] <- xy
    part_bends[[length(part_bends) + 1]] <- bends
    used <- c(used, u)
  }
  packed <- .kegg_topology_pack(parts)
  bends <- do.call(rbind, lapply(seq_along(part_bends), function(i) {
    b <- part_bends[[i]]
    if (is.null(b) || nrow(b) == 0) return(NULL)
    b$x <- b$x + packed$shift[i, 1]
    b$y <- b$y + packed$shift[i, 2]
    b
  }))
  if (is.null(bends)) {
    bends <- data.frame(from = character(0), to = character(0), k = integer(0), x = numeric(0), y = numeric(0),
                        stringsAsFactors = FALSE)
  }
  used <- setdiff(unique(used), "single")
  list(xy = packed$xy, bends = bends,
       used = if (length(used) == 0) "none" else if (length(used) == 1) used else "mixed")
}

#' Interior points (dummy vertices) of each edge of a Sugiyama layout
#' @param h The laid-out graph (vertex names = boxes).
#' @param ext `extd_graph` of [igraph::layout_with_sugiyama()]; its edges
#'   carry the original edge id in `orig`.
#' @param XY Coordinates of all extended-graph vertices (original first).
#' @return `data.frame(from, to, k, x, y)` -- `k` orders the bends of one
#'   edge from `from` to `to`.
#' @keywords internal
.kegg_topology_sugiyama_bends <- function(h, ext, XY) {
  n <- igraph::vcount(h)
  seg <- igraph::as_edgelist(ext, names = FALSE)
  orig <- igraph::E(ext)$orig
  ends <- igraph::ends(h, igraph::E(h), names = FALSE)
  nm <- igraph::V(h)$name
  out <- list()
  for (eid in unique(orig[duplicated(orig)])) {
    s <- seg[orig == eid, , drop = FALSE]
    cur <- ends[eid, 1]
    used <- rep(FALSE, nrow(s))
    pts <- integer(0)
    for (step in seq_len(nrow(s))) {
      j <- which(!used & (s[, 1] == cur | s[, 2] == cur))[1]
      if (is.na(j)) break
      used[j] <- TRUE
      nxt <- if (s[j, 1] == cur) s[j, 2] else s[j, 1]
      if (nxt <= n) break
      pts <- c(pts, nxt)
      cur <- nxt
    }
    if (length(pts)) {
      out[[length(out) + 1]] <- data.frame(from = nm[ends[eid, 1]], to = nm[ends[eid, 2]], k = seq_along(pts),
                                           x = XY[pts, 1], y = XY[pts, 2], stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

#' Is a directed component close enough to acyclic for a layered layout?
#' (feedback arc set of at most 20% of its edges)
#' @keywords internal
.kegg_topology_is_layered <- function(h) {
  m <- igraph::ecount(h)
  if (m == 0) return(FALSE)
  fas <- tryCatch(length(igraph::feedback_arc_set(h, algo = "approx_eades")), error = function(e) m)
  fas <= 0.2 * m
}

#' Rescale a force-directed layout so its extent is about sqrt(n) units,
#' stretched horizontally like the layered layout (labels are wide)
#' @keywords internal
.kegg_topology_unit_scale <- function(xy, n) {
  xy <- sweep(xy, 2, colMeans(xy))
  ext <- max(apply(xy, 2, function(v) diff(range(v))), 1e-9)
  xy <- xy / ext * (1.4 * sqrt(n))
  xy[, 1] <- xy[, 1] * 1.8
  xy
}

#' Pack component layouts in rows (largest first), in abstract units
#' @return `list(xy, shift)` -- `shift` is a two-column matrix, the
#'   translation applied to each part.
#' @keywords internal
.kegg_topology_pack <- function(parts, gap_x = 2.6, gap_y = 1.2) {
  if (length(parts) == 1) return(list(xy = parts[[1]], shift = matrix(0, 1, 2)))
  ext <- lapply(parts, function(p) c(w = diff(range(p[, 1])), h = diff(range(p[, 2]))))
  total_area <- sum(vapply(ext, function(e) (e[["w"]] + gap_x) * (e[["h"]] + gap_y), numeric(1)))
  row_w <- max(max(vapply(ext, `[[`, numeric(1), "w")), sqrt(total_area * 2.2))
  out <- list()
  shift <- matrix(0, length(parts), 2)
  x0 <- 0
  y0 <- 0
  row_h <- 0
  for (i in seq_along(parts)) {
    p <- parts[[i]]
    w <- ext[[i]][["w"]]
    h <- ext[[i]][["h"]]
    if (x0 > 0 && x0 + w > row_w) {
      x0 <- 0
      y0 <- y0 - row_h - gap_y
      row_h <- 0
    }
    shift[i, ] <- c(x0 - min(p[, 1]), y0 - max(p[, 2]))
    q <- p
    q[, 1] <- p[, 1] + shift[i, 1]
    q[, 2] <- p[, 2] + shift[i, 2]
    out[[i]] <- q
    x0 <- x0 + w + gap_x
    row_h <- max(row_h, h)
  }
  list(xy = do.call(rbind, out), shift = shift)
}

#' Fit abstract coordinates into a `w` x `h` inch panel (independent x/y
#' scaling, capped so a small graph is not blown up)
#' @return `list(xy, tf)` -- `tf(m)` applies the same transform to any
#'   two-column matrix (the edge bends).
#' @keywords internal
.kegg_topology_fit <- function(xy, w, h, pad_x, pad_y) {
  mk <- function(v, to, pad, cap) {
    r <- range(v)
    avail <- max(to - 2 * pad, 0.1)
    if (diff(r) == 0) return(function(z) rep(to / 2, length(z)))
    s <- min(avail / diff(r), cap)
    m <- mean(r)
    function(z) (z - m) * s + to / 2
  }
  fx <- mk(xy[, 1], w, pad_x, cap = 0.75)
  fy <- mk(xy[, 2], h, pad_y, cap = 0.55)
  tf <- function(m) cbind(fx(m[, 1]), fy(m[, 2]))
  out <- tf(xy)
  rownames(out) <- rownames(xy)
  list(xy = out, tf = tf)
}

#' Edge paths clipped at the label boxes, tee bars for inhibitions,
#' parallel offsets for straight reciprocal pairs
#'
#' @param e Box-level edges (`from_box`, `to_box`, `effect`, `mode`,
#'   `subtypes`).
#' @param nodes Boxes with `box`, `x`, `y`, `hw`, `hh` (inches).
#' @param bends `data.frame(from, to, k, x, y)` of interior points (inches).
#' @return `list(edges, paths, tees)`: `edges` one row per edge (`edge_id`,
#'   `from`, `to`, `effect`, `mode`, `subtypes`, `head`); `paths` the drawn
#'   polyline points (`edge_id`, `ord`, `x`, `y` + the edge columns);
#'   `tees` one bar per inhibition-type edge (`x`, `y`, `xend`, `yend`).
#' @keywords internal
.kegg_topology_edge_geometry <- function(e, nodes, bends = NULL) {
  empty_paths <- data.frame(edge_id = integer(0), ord = integer(0), x = numeric(0), y = numeric(0),
                            effect = character(0), mode = character(0), head = character(0), stringsAsFactors = FALSE)
  if (nrow(e) == 0) {
    return(list(
      edges = data.frame(edge_id = integer(0), from = character(0), to = character(0), effect = character(0),
                         mode = character(0), subtypes = character(0), head = character(0), stringsAsFactors = FALSE),
      paths = empty_paths,
      tees = data.frame(edge_id = integer(0), x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
                        effect = character(0), stringsAsFactors = FALSE)
    ))
  }
  edges <- data.frame(edge_id = seq_len(nrow(e)), from = e$from_box, to = e$to_box, effect = e$effect, mode = e$mode,
                      subtypes = e$subtypes, stringsAsFactors = FALSE)
  edges$head <- ifelse(edges$effect %in% c("inhibition", "repression"), "tee",
                       ifelse(edges$effect %in% c("binding/association", "other"), "none", "arrow"))
  pair <- paste(edges$from, edges$to)
  rev_pair <- paste(edges$to, edges$from)
  exit <- function(hw, hh, ux, uy) min(if (abs(ux) > 1e-9) hw / abs(ux) else Inf, if (abs(uy) > 1e-9) hh / abs(uy) else Inf)
  gap <- 0.02
  tee_half <- 0.05
  paths <- vector("list", nrow(edges))
  tees <- vector("list", nrow(edges))
  for (i in seq_len(nrow(edges))) {
    a <- match(edges$from[i], nodes$box)
    b <- match(edges$to[i], nodes$box)
    mid <- if (!is.null(bends) && nrow(bends)) bends[bends$from == edges$from[i] & bends$to == edges$to[i], , drop = FALSE] else NULL
    mid <- if (!is.null(mid) && nrow(mid)) mid[order(mid$k), c("x", "y"), drop = FALSE] else NULL
    px <- c(nodes$x[a], if (!is.null(mid)) mid$x, nodes$x[b])
    py <- c(nodes$y[a], if (!is.null(mid)) mid$y, nodes$y[b])
    m <- length(px)
    ## start: leave the source box along the first segment
    dx <- px[2] - px[1]; dy <- py[2] - py[1]; l1 <- sqrt(dx^2 + dy^2)
    ## end: enter the target box along the last segment
    ex <- px[m] - px[m - 1]; ey <- py[m] - py[m - 1]; l2 <- sqrt(ex^2 + ey^2)
    if (l1 < 1e-9 || l2 < 1e-9) next
    u1 <- c(dx, dy) / l1
    u2 <- c(ex, ey) / l2
    s0 <- exit(nodes$hw[a], nodes$hh[a], u1[1], u1[2]) + gap
    s1 <- exit(nodes$hw[b], nodes$hh[b], u2[1], u2[2]) + gap
    if (m == 2 && l1 <= s0 + s1 + 0.02) next
    px[1] <- px[1] + u1[1] * min(s0, 0.9 * l1)
    py[1] <- py[1] + u1[2] * min(s0, 0.9 * l1)
    px[m] <- px[m] - u2[1] * min(s1, 0.9 * l2)
    py[m] <- py[m] - u2[2] * min(s1, 0.9 * l2)
    if (m == 2 && rev_pair[i] %in% pair) {
      off <- 0.035
      px <- px + u1[2] * off
      py <- py - u1[1] * off
    }
    paths[[i]] <- data.frame(edge_id = i, ord = seq_len(m), x = px, y = py, effect = edges$effect[i],
                             mode = edges$mode[i], head = edges$head[i], stringsAsFactors = FALSE)
    if (edges$head[i] == "tee") {
      tees[[i]] <- data.frame(edge_id = i, x = px[m] - u2[2] * tee_half, y = py[m] + u2[1] * tee_half,
                              xend = px[m] + u2[2] * tee_half, yend = py[m] - u2[1] * tee_half,
                              effect = edges$effect[i], stringsAsFactors = FALSE)
    }
  }
  paths <- do.call(rbind, paths)
  tees <- do.call(rbind, tees)
  list(
    edges = edges,
    paths = if (is.null(paths)) empty_paths else paths,
    tees = if (is.null(tees)) data.frame(edge_id = integer(0), x = numeric(0), y = numeric(0), xend = numeric(0),
                                         yend = numeric(0), effect = character(0), stringsAsFactors = FALSE) else tees
  )
}

#' Build the topology ggplot (one facet per pathway)
#' @keywords internal
.kegg_topology_ggplot <- function(nodes, edges, paths, tees, summ, cond, ncol, panel_w, panel_h, engine, width,
                                  restricted = FALSE) {
  strip <- stats::setNames(sprintf("%s (%s)\n%d of %d genes hit by the extract (%.0f%%)", summ$pathway_title, summ$pathway_id,
                                   summ$n_hit, summ$n_genes, 100 * summ$share_hit), summ$pathway_id)
  lev <- unname(strip[summ$pathway_id])
  nodes$panel <- factor(strip[nodes$pathway_id], levels = lev)
  if (nrow(paths)) paths$panel <- factor(strip[paths$pathway_id], levels = lev)
  if (nrow(tees)) tees$panel <- factor(strip[tees$pathway_id], levels = lev)

  eff_cols <- c(activation = "#009E73", expression = "#0072B2", inhibition = "#D55E00", repression = "#CC79A7",
                `indirect effect` = "#56B4E9", modification = "#E69F00", `binding/association` = "grey35", other = "grey60")
  eff_labs <- c(activation = "activation  -->", expression = "expression (GErel)  -->", inhibition = "inhibition  --|",
                repression = "repression (GErel)  --|", `indirect effect` = "indirect effect  ..>",
                modification = "phosphorylation etc., sign not stated  -->", `binding/association` = "binding / association  ---",
                other = "other (state change, via compound...)")
  mode_lty <- c(phosphorylation = "42", dephosphorylation = "4212", `other modification` = "longdash",
                `via compound` = "dotdash", indirect = "11", direct = "solid")
  mode_labs <- c(phosphorylation = "phosphorylation (+p)", dephosphorylation = "dephosphorylation (-p)",
                 `other modification` = "ubiquitination / glycosylation / methylation",
                 `via compound` = "via a compound", indirect = "indirect", direct = "direct / unspecified")
  eff_present <- intersect(names(eff_cols), unique(edges$effect))
  mode_present <- intersect(names(mode_lty), unique(edges$mode))

  p <- ggplot2::ggplot()
  if (nrow(paths)) {
    path_aes <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id, colour = .data$effect, linetype = .data$mode)
    paths <- paths[order(paths$pathway_id, paths$edge_id, paths$ord), , drop = FALSE]
    paths$edge_id <- paste(paths$pathway_id, paths$edge_id)
    arr <- paths[paths$head == "arrow", , drop = FALSE]
    rest <- paths[paths$head != "arrow", , drop = FALSE]
    if (nrow(rest)) p <- p + ggplot2::geom_path(data = rest, path_aes, linewidth = 0.45, alpha = 0.9, lineend = "butt")
    if (nrow(arr)) {
      p <- p + ggplot2::geom_path(data = arr, path_aes, linewidth = 0.45, alpha = 0.95, linejoin = "round",
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.065, "inches"), type = "closed", angle = 22))
    }
    if (nrow(tees)) {
      p <- p + ggplot2::geom_segment(data = tees, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                                                               colour = .data$effect), linewidth = 0.9, show.legend = FALSE)
    }
    p <- p +
      ggplot2::scale_colour_manual(values = eff_cols, breaks = eff_present, labels = eff_labs[eff_present],
                                   name = "Relation (KEGG)", drop = TRUE,
                                   guide = ggplot2::guide_legend(order = 3, override.aes = list(linewidth = 0.9, linetype = "solid"))) +
      ggplot2::scale_linetype_manual(values = mode_lty, breaks = mode_present, labels = mode_labs[mode_present],
                                     name = "Mechanism", drop = TRUE,
                                     guide = ggplot2::guide_legend(order = 4, override.aes = list(linewidth = 0.6, colour = "grey25")))
  }

  hits <- nodes[nodes$is_target, , drop = FALSE]
  non_hits <- nodes[!nodes$is_target, , drop = FALSE]
  nodes$tooltip <- sprintf("%s\n%s%s", nodes$symbols, ifelse(nodes$is_target, sprintf("hit by %d compound(s), max p = %.2f",
                                                                                        nodes$n_compounds, nodes$max_weight),
                                                                "not a predicted target"), "")
  if (nrow(non_hits)) {
    p <- p + ggplot2::geom_label(data = non_hits, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                                 size = .kegg_topology_label_size(0, FALSE), fill = "white", colour = "grey40", linewidth = 0.25,
                                 label.padding = ggplot2::unit(0.25, "lines"), label.r = ggplot2::unit(0.08, "lines"),
                                 show.legend = FALSE)
  }
  if (nrow(hits)) {
    dark <- .kegg_topology_fill_is_dark(hits$max_weight)
    legend_done <- FALSE
    for (is_dark in c(FALSE, TRUE)) {
      d <- hits[dark == is_dark, , drop = FALSE]
      if (nrow(d) == 0) next
      p <- p + ggplot2::geom_label(data = d, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                                          fill = .data$max_weight, size = .data$size),
                                   show.legend = c(fill = !legend_done, size = !legend_done, colour = FALSE,
                                                   linetype = FALSE, linewidth = FALSE),
                                   border.colour = "grey10", text.colour = if (is_dark) "white" else "grey5",
                                   fontface = "bold", linewidth = 0.4, label.padding = ggplot2::unit(0.25, "lines"),
                                   label.r = ggplot2::unit(0.08, "lines"))
      legend_done <- TRUE
    }
    size_br <- .kegg_topology_size_breaks(hits$n_compounds)
    p <- p +
      ggplot2::scale_fill_gradientn(colours = .kegg_topology_fill_colours(), limits = c(0.5, 1), oob = .pathway_squish,
                                    name = "Max. probability\n(extract's compounds)",
                                    guide = ggplot2::guide_colourbar(order = 1, barheight = ggplot2::unit(2.4, "cm"))) +
      ggplot2::scale_size_identity(guide = if (length(size_br$counts) > 1) "legend" else "none",
                                   breaks = size_br$sizes, labels = size_br$counts,
                                   name = "Compounds\nhitting the gene")
    if (length(size_br$counts) > 1) {
      p <- p + ggplot2::guides(size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "grey85", label = "GENE")))
    }
  }
  if (length(unique(nodes$panel)) > 1) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = ncol)
  }
  single <- length(unique(nodes$panel)) == 1
  title <- if (single) paste0(summ$pathway_title[1], " (", summ$pathway_id[1], ") -- KEGG topology, ", cond)
           else paste0("KEGG pathway topology -- ", cond)
  sub <- paste0(
    if (single) sprintf("%d of %d genes are predicted targets of the extract (%.0f%%); %d relations drawn. ",
                        summ$n_hit[1], summ$n_genes[1], 100 * summ$share_hit[1], summ$n_relations[1]) else
      sprintf("The %d pathways with the most extract targets among their KEGG relations. ", nrow(summ)),
    "Boxes = genes (paralogs with identical relations merged, e.g. MAPK1/3); filled boxes are targets of the extract, ",
    "fill = highest predicted probability, box size = number of compounds hitting it",
    if (any(!nodes$is_target)) "; white boxes are not targeted. " else ". ",
    "Edges = KGML relations (PPrel/GErel/ECrel), colour and end = effect, line type = mechanism",
    if (restricted) paste0(
      ". Only relations between two targets of the extract are stored (network_kegg_topology(restrict_to_network = TRUE)); ",
      "rerun it with restrict_to_network = FALSE to draw the whole pathway."
    ) else "; the layout runs from upstream (left) to downstream (right) where the relations allow it."
  )
  p +
    ggplot2::scale_x_continuous(limits = c(0, panel_w), expand = ggplot2::expansion(add = 0.02)) +
    ggplot2::scale_y_continuous(limits = c(0, panel_h), expand = ggplot2::expansion(add = 0.02)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = .plot_wrap(title, .plot_wrap_width(width, 12)), subtitle = .plot_wrap(sub, .plot_wrap_width(width, 8))) +
    .network_view_theme("right") +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0, margin = ggplot2::margin(2, 2, 4, 2)),
      panel.border = if (single) ggplot2::element_blank() else ggplot2::element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
      panel.spacing = ggplot2::unit(0.15, "in"),
      legend.key.width = ggplot2::unit(0.9, "cm")
    )
}

#' Sequential fill of the target boxes (light yellow -> dark red)
#' @keywords internal
.kegg_topology_fill_colours <- function() c("#FFF3C4", "#FDD08A", "#F99B5B", "#E4572E", "#A11D21")

#' Is the box fill for these probabilities dark enough for white text?
#' @keywords internal
.kegg_topology_fill_is_dark <- function(w) {
  ramp <- grDevices::colorRamp(.kegg_topology_fill_colours())
  v <- pmin(pmax((w - 0.5) / 0.5, 0), 1)
  v[is.na(v)] <- 0
  rgb <- ramp(v)
  (0.299 * rgb[, 1] + 0.587 * rgb[, 2] + 0.114 * rgb[, 3]) / 255 < 0.5
}

#' Squish out-of-bounds values into the limits (avoids a `scales` import)
#' @keywords internal
.pathway_squish <- function(x, range = c(0, 1), only.finite = TRUE) {
  finite <- if (only.finite) is.finite(x) else TRUE
  x[finite & x < range[1]] <- range[1]
  x[finite & x > range[2]] <- range[2]
  x
}

#' Legend breaks for the box-size (compound count) legend
#' @keywords internal
.kegg_topology_size_breaks <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(list(counts = integer(0), sizes = numeric(0)))
  br <- unique(round(stats::quantile(counts, c(0, 0.5, 1), names = FALSE)))
  br <- unique(pmax(1L, as.integer(br)))
  sizes <- .kegg_topology_label_size(br, TRUE)
  keep <- !duplicated(round(sizes, 6))
  list(counts = br[keep], sizes = sizes[keep])
}

# ---------------------------------------------------------------------------
# kegg_routes() / kegg_routes_table() / plot_kegg_routes()
# ---------------------------------------------------------------------------

#' Routes along KEGG relations from the extract's targets to disease genes
#'
#' @description
#' "How do the extract's targets reach the disease?" -- a route planner
#' over the directed gene-gene relations KEGG draws inside its pathways
#' (as stored by [network_kegg_topology()]). All pathways of `condition`
#' are merged into one directed graph (a gene pair seen in several
#' pathways becomes one edge that remembers them all). Starting points
#' ("sources") are the extract's predicted targets; end points
#' ("destinations") are the genes of a disease gene set (`disease`) or an
#' explicit list (`to`). For every reachable (source, destination) pair the
#' shortest directed route is found (breadth-first search); the best pairs
#' are kept, each with up to `n_alternatives` equally short alternatives.
#'
#' @section Relations, direction and sign:
#' Each KGML relation keeps its subtype (activation, inhibition,
#' phosphorylation, expression, ...). A hop is `+1` for activation /
#' expression, `-1` for inhibition / repression and unknown otherwise
#' (phosphorylation, indirect effect, binding, ...: KEGG does not state
#' whether they activate). The route sign is the product of its hops:
#' `"+"` (the source, when active, promotes the destination), `"-"`
#' (represses it), or `"?"` when any hop is unsigned. The sign describes the
#' relation chain, not what a compound does to its target (which is
#' generally unknown). Binding/association relations have no direction and
#' can be travelled both ways (`undirected_binding = TRUE`).
#'
#' @section Ranking:
#' Only shortest routes are considered, of `min_hops` to `max_len` hops.
#' Each (source, destination) pair gets
#' \deqn{\mathrm{score} = \frac{w_s \, \log_2(1 + c_s) \, a_d}{\mathrm{hops}}}
#' with \eqn{w_s} the highest predicted probability among the compounds
#' hitting the source, \eqn{c_s} the number of those compounds and
#' \eqn{a_d} the destination's Open Targets association score (1 when not
#' available, e.g. with `to`).
#'
#' * `rank_by = "destination"` (default) answers "how does the extract
#'   reach the genes that matter most for the disease": destinations are
#'   taken by decreasing association score (ties: best pair score), each
#'   with its best-scoring route.
#' * `rank_by = "pair"` takes pairs by decreasing score -- the shortest,
#'   best-supported routes overall (typically one-hop neighbours).
#'
#' In both, a source is used at most `max_per_gene` times (and, with
#' `"pair"`, a destination too), so the result is not dominated by one hub.
#' Among equally short alternatives, those passing through more of the
#' extract's other targets come first. Sources that are themselves
#' destinations need no route; they are returned in
#' `attr(., "direct_hits")`.
#'
#' @inheritParams network_build
#' @param condition Character scalar with rows in `network_kegg_topology`;
#'   `NULL` (default) when only one condition has topology.
#' @param disease `NULL` or a `disease_id` of
#'   `patliRResults(proj, "disease_genes")`. One of `disease` / `to` is
#'   required.
#' @param to `NULL` or a character vector of destination genes (UniProt
#'   accessions, KEGG gene IDs such as `"hsa:5594"`, or gene symbols).
#'   When both are given, destinations are the disease genes in `to`.
#' @param from `NULL` (default, every target of the extract in the
#'   topology) or a character vector restricting the sources (same ID types
#'   as `to`).
#' @param compounds `NULL` (default) or `compounds(proj)$id` values: only
#'   targets of these compounds are sources.
#' @param max_len Integer, default `6`: longest route (hops) considered.
#' @param top_n Integer, default `10`: number of (source, destination) pairs
#'   returned.
#' @param n_alternatives Integer, default `3`: equally short routes kept
#'   per pair.
#' @param max_per_gene Integer, default `2`.
#' @param min_hops Integer, default `1`: shortest route length kept (e.g.
#'   `2` to look only at mechanisms through at least one intermediate gene).
#' @param rank_by `"destination"` (default) or `"pair"`, see the section
#'   above.
#' @param undirected_binding Logical, default `TRUE`.
#'
#' @return A `data.frame`, one row per route, best first: `route_id`
#'   (`"R1"`, `"R1b"`, ...), `pair_rank`, `alternative`, `source`,
#'   `source_symbol`, `compounds`, `n_compounds`, `source_max_weight`,
#'   `destination`, `destination_symbol`, `destination_score`, `n_hops`,
#'   `path` (symbols with the hop effects, e.g. `"EGFR --> GRB2 --- SOS1"`),
#'   `path_kegg`, `sign`, `pathways` (KEGG IDs the hops belong to),
#'   `pathway_titles`, `score`. Attributes: `hops` (one row per hop: route,
#'   hop number, genes, effect, sign, pathways of the hop), `direct_hits`
#'   (targets that are destinations themselves), `graph_size`.
#'
#' @examples
#' \dontrun{
#' proj <- network_kegg_topology(proj, condition = "FLO-ET", restrict_to_network = FALSE)
#' proj <- disease_genes_fetch(proj, "hypertension")
#' r <- kegg_routes(proj, condition = "FLO-ET", disease = "MONDO_0005044")
#' kegg_routes_table(r)
#' }
#' @seealso [plot_kegg_routes()], [network_kegg_topology()], [network_proximity()]
#' @export
kegg_routes <- function(proj, condition = NULL, disease = NULL, to = NULL, from = NULL, compounds = NULL,
                        max_len = 6, top_n = 10, n_alternatives = 3, max_per_gene = 2,
                        min_hops = 1, rank_by = c("destination", "pair"), undirected_binding = TRUE) {
  stopifnot(is(proj, "PatliRProject"))
  rank_by <- match.arg(rank_by)
  .pathway_check_count(min_hops, "min_hops", min = 1)
  if (min_hops > max_len) cli::cli_abort("{.arg min_hops} must not exceed {.arg max_len}.")
  .pathway_check_count(max_len, "max_len", min = 1)
  .pathway_check_count(top_n, "top_n", min = 1)
  .pathway_check_count(n_alternatives, "n_alternatives", min = 1)
  .pathway_check_count(max_per_gene, "max_per_gene", min = 1)
  .pathway_check_flag(undirected_binding, "undirected_binding")
  for (a in c("to", "from", "compounds")) {
    v <- get(a)
    if (!is.null(v) && (!is.character(v) || length(v) == 0 || anyNA(v))) {
      cli::cli_abort("{.arg {a}} must be NULL or a character vector without NA.")
    }
  }
  if (!is.null(disease) && (!is.character(disease) || length(disease) != 1L || is.na(disease) || !nzchar(disease))) {
    cli::cli_abort("{.arg disease} must be NULL or a single {.field disease_id} string.")
  }
  if (is.null(disease) && is.null(to)) cli::cli_abort("Give the destinations: {.arg disease} and/or {.arg to}.")

  topo_all <- patliRResults(proj, "network_kegg_topology")
  if (is.null(topo_all) || nrow(topo_all) == 0) {
    cli::cli_abort(c("No {.val network_kegg_topology} entry in {.arg proj}.",
                     "i" = "Run {.fn network_kegg_topology} first (ideally with {.code restrict_to_network = FALSE})."))
  }
  cond <- .pathway_single_condition(unique(topo_all$condition), condition, "kegg_routes", "network_kegg_topology")
  topo <- topo_all[topo_all$condition == cond, , drop = FALSE]
  edges_all <- patliRResults(proj, "network_edges")
  ct <- if (is.null(edges_all)) NULL else edges_all[edges_all$condition == cond, , drop = FALSE]
  if (!is.null(compounds)) {
    unknown <- setdiff(compounds, ct$compound_id)
    if (length(unknown)) cli::cli_warn("Compound{?s} {.val {unknown}} {?has/have} no target in condition {.val {cond}}.")
    ct <- ct[ct$compound_id %in% compounds, , drop = FALSE]
  }
  genes <- .kegg_topology_gene_stats(topo, ct)
  genes$symbol <- .kegg_topology_symbols(proj, cond, genes)

  ## destinations
  dest <- NULL
  dest_score <- stats::setNames(rep(1, nrow(genes)), genes$kegg)
  if (!is.null(disease)) {
    dg <- patliRResults(proj, "disease_genes")
    if (is.null(dg) || nrow(dg) == 0 || !disease %in% dg$disease_id) {
      cli::cli_abort(c("{.arg disease} {.val {disease}} is not in {.val disease_genes}.",
                       "i" = "Available: {.val {unique(dg$disease_id)}}; see {.fn disease_genes_fetch}."))
    }
    dg <- dg[dg$disease_id == disease, , drop = FALSE]
    dest <- .kegg_routes_resolve(genes, c(dg$uniprot_id, stats::na.omit(dg$gene_symbol)), warn = FALSE)
    sc <- .kegg_routes_gene_score(genes, dg)
    dest_score[] <- NA_real_
    dest_score[names(sc)] <- sc
  }
  if (!is.null(to)) {
    to_k <- .kegg_routes_resolve(genes, to, warn = TRUE, arg = "to")
    dest <- if (is.null(dest)) to_k else intersect(dest, to_k)
  }
  if (length(dest) == 0) {
    cli::cli_abort("None of the destination genes appears in condition {.val {cond}}'s KEGG topology.")
  }

  ## sources
  src <- genes$kegg[genes$is_target]
  if (!is.null(from)) src <- intersect(src, .kegg_routes_resolve(genes, from, warn = TRUE, arg = "from"))
  if (length(src) == 0) {
    cli::cli_abort("No target of the extract (in scope) appears in condition {.val {cond}}'s KEGG topology.")
  }
  direct <- intersect(src, dest)

  rel <- .kegg_topology_collapse_relations(topo)
  g_edges <- .kegg_routes_edges(rel, undirected_binding)
  g <- igraph::graph_from_data_frame(g_edges[, c("from", "to")], directed = TRUE,
                                     vertices = data.frame(name = genes$kegg, stringsAsFactors = FALSE))
  pairs <- .kegg_routes_pairs(g, src, dest, max_len, genes, dest_score, min_hops = min_hops)
  if (rank_by == "destination" && nrow(pairs)) {
    a <- dest_score[pairs$destination]
    a[is.na(a)] <- 0
    pairs <- pairs[order(-a, -pairs$score, pairs$n_hops, pairs$source), , drop = FALSE]
  }
  titles <- unique(topo[, c("pathway_id", "pathway_title")])
  title_of <- stats::setNames(titles$pathway_title, titles$pathway_id)
  out <- .kegg_routes_expand(g, g_edges, pairs, genes, title_of, top_n, n_alternatives, max_per_gene,
                             max_per_destination = if (rank_by == "destination") 1L else max_per_gene)
  if (nrow(out) && !is.null(disease)) out$destination_score <- unname(dest_score[out$destination])
  attr(out, "direct_hits") <- genes[genes$kegg %in% direct, c("kegg", "symbol", "n_compounds", "max_weight"), drop = FALSE]
  attr(out, "graph_size") <- c(genes = igraph::vcount(g), relations = nrow(rel), sources = length(src),
                               destinations = length(dest))
  attr(out, "condition") <- cond
  attr(out, "disease") <- disease
  out
}

#' Resolve gene identifiers (KEGG gene ID, UniProt, symbol) to KEGG gene
#' IDs of the topology
#' @keywords internal
.kegg_routes_resolve <- function(genes, ids, warn = TRUE, arg = "to") {
  ids <- unique(trimws(ids[!is.na(ids) & nzchar(ids)]))
  up <- strsplit(genes$uniprots, ";", fixed = TRUE)
  up_tab <- data.frame(kegg = rep(genes$kegg, lengths(up)), uniprot = unlist(up), stringsAsFactors = FALSE)
  hit <- lapply(ids, function(x) {
    unique(c(genes$kegg[genes$kegg == x], up_tab$kegg[up_tab$uniprot == x],
             genes$kegg[!is.na(genes$symbol) & toupper(genes$symbol) == toupper(x)]))
  })
  miss <- ids[lengths(hit) == 0]
  if (warn && length(miss)) {
    cli::cli_warn("{length(miss)} {.arg {arg}} gene{?s} not in the KEGG topology: {.val {utils::head(miss, 10)}}.")
  }
  unique(unlist(hit))
}

#' Highest disease association score per KEGG gene
#' @keywords internal
.kegg_routes_gene_score <- function(genes, dg) {
  up <- strsplit(genes$uniprots, ";", fixed = TRUE)
  up_tab <- data.frame(kegg = rep(genes$kegg, lengths(up)), uniprot = unlist(up), stringsAsFactors = FALSE)
  m <- merge(up_tab, dg[, c("uniprot_id", "association_score")], by.x = "uniprot", by.y = "uniprot_id")
  m <- m[!is.na(m$association_score), , drop = FALSE]
  if (nrow(m) == 0) return(numeric(0))
  s <- tapply(m$association_score, m$kegg, max)
  stats::setNames(as.numeric(s), names(s))
}

#' Routable edges: collapsed relations with a sign; binding added both ways
#' @return `data.frame(from, to, effect, sign, pathways, reversed)`.
#' @keywords internal
.kegg_routes_edges <- function(rel, undirected_binding = TRUE) {
  sign <- ifelse(rel$effect %in% c("activation", "expression"), 1L,
                 ifelse(rel$effect %in% c("inhibition", "repression"), -1L, NA_integer_))
  e <- data.frame(from = rel$from, to = rel$to, effect = rel$effect, sign = sign, pathways = rel$pathways,
                  reversed = FALSE, stringsAsFactors = FALSE)
  if (undirected_binding) {
    b <- e[e$effect == "binding/association", , drop = FALSE]
    if (nrow(b)) {
      b <- data.frame(from = b$to, to = b$from, effect = b$effect, sign = b$sign, pathways = b$pathways,
                      reversed = TRUE, stringsAsFactors = FALSE)
      b <- b[!paste(b$from, b$to) %in% paste(e$from, e$to), , drop = FALSE]
      e <- rbind(e, b)
    }
  }
  rownames(e) <- NULL
  e
}

#' Reachable (source, destination) pairs with their route length and score
#' @return `data.frame(source, destination, n_hops, score)`, best first.
#' @keywords internal
.kegg_routes_pairs <- function(g, src, dest, max_len, genes, dest_score, min_hops = 1) {
  src_only <- setdiff(src, dest)
  empty <- data.frame(source = character(0), destination = character(0), n_hops = integer(0), score = numeric(0),
                      stringsAsFactors = FALSE)
  if (length(src_only) == 0) return(empty)
  d <- igraph::distances(g, v = src_only, to = dest, mode = "out")
  idx <- which(is.finite(d) & d >= min_hops & d <= max_len, arr.ind = TRUE)
  if (nrow(idx) == 0) return(empty)
  out <- data.frame(source = rownames(d)[idx[, 1]], destination = colnames(d)[idx[, 2]], n_hops = as.integer(d[idx]),
                    stringsAsFactors = FALSE)
  gi <- match(out$source, genes$kegg)
  w <- genes$max_weight[gi]
  w[is.na(w)] <- 0.5
  a <- dest_score[out$destination]
  a[is.na(a)] <- 1
  out$score <- w * log2(1 + genes$n_compounds[gi]) * a / out$n_hops
  out <- out[order(-out$score, out$n_hops, out$source, out$destination), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Pick the pairs (per-gene caps) and expand them into routes and hops
#' @keywords internal
.kegg_routes_expand <- function(g, g_edges, pairs, genes, title_of, top_n, n_alternatives, max_per_gene,
                                max_per_destination = max_per_gene) {
  empty <- .kegg_routes_empty()
  if (nrow(pairs) == 0) {
    attr(empty, "hops") <- .kegg_routes_empty_hops()
    return(empty)
  }
  n_src <- integer(0)
  n_dst <- integer(0)
  chosen <- integer(0)
  for (i in seq_len(nrow(pairs))) {
    s <- pairs$source[i]
    d <- pairs$destination[i]
    if (max(0L, n_src[s], na.rm = TRUE) >= max_per_gene || max(0L, n_dst[d], na.rm = TRUE) >= max_per_destination) next
    chosen <- c(chosen, i)
    n_src[s] <- max(0L, n_src[s], na.rm = TRUE) + 1L
    n_dst[d] <- max(0L, n_dst[d], na.rm = TRUE) + 1L
    if (length(chosen) >= top_n) break
  }
  pairs <- pairs[chosen, , drop = FALSE]
  ekey <- paste(g_edges$from, g_edges$to, sep = "\r")
  sym <- stats::setNames(genes$symbol, genes$kegg)
  vnames <- igraph::V(g)$name
  arrow_of <- c(activation = "-->", expression = "-->", inhibition = "--|", repression = "--|",
                `indirect effect` = "..>", modification = "-p>", `binding/association` = "---", other = "-?-")
  routes <- list()
  hops <- list()
  for (k in seq_len(nrow(pairs))) {
    sp <- igraph::all_shortest_paths(g, from = pairs$source[k], to = pairs$destination[k], mode = "out")$vpaths
    paths <- lapply(sp, function(v) vnames[as.integer(v)])
    support <- vapply(paths, function(p) sum(genes$is_target[match(p[-c(1, length(p))], genes$kegg)]), numeric(1))
    ord <- order(-support, vapply(paths, paste, character(1), collapse = " "))
    paths <- utils::head(paths[ord], n_alternatives)
    for (a in seq_along(paths)) {
      p <- paths[[a]]
      rid <- paste0("R", k, if (a > 1) letters[a] else "")
      e <- g_edges[match(paste(p[-length(p)], p[-1], sep = "\r"), ekey), , drop = FALSE]
      hop_pw <- e$pathways
      pw <- sort(unique(unlist(strsplit(hop_pw, "; ", fixed = TRUE))))
      sgn <- if (anyNA(e$sign)) "?" else if (prod(e$sign) > 0) "+" else "-"
      txt <- sym[p[1]]
      for (h in seq_len(nrow(e))) txt <- paste(txt, arrow_of[[e$effect[h]]], sym[p[h + 1]])
      gi <- match(p[1], genes$kegg)
      routes[[length(routes) + 1]] <- data.frame(
        route_id = rid, pair_rank = k, alternative = a, source = p[1], source_symbol = unname(sym[p[1]]),
        compounds = genes$compounds[gi], n_compounds = genes$n_compounds[gi], source_max_weight = genes$max_weight[gi],
        destination = p[length(p)], destination_symbol = unname(sym[p[length(p)]]),
        destination_score = NA_real_, n_hops = length(p) - 1L, path = unname(txt), path_kegg = paste(p, collapse = ">"),
        sign = sgn, pathways = paste(pw, collapse = "; "),
        pathway_titles = paste(unname(ifelse(is.na(title_of[pw]), pw, title_of[pw])), collapse = "; "),
        score = pairs$score[k], stringsAsFactors = FALSE
      )
      hops[[length(hops) + 1]] <- data.frame(
        route_id = rid, hop = seq_len(nrow(e)), from = p[-length(p)], to = p[-1],
        from_symbol = unname(sym[p[-length(p)]]), to_symbol = unname(sym[p[-1]]), effect = e$effect, sign = e$sign,
        reversed = e$reversed, pathways = hop_pw, stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, routes)
  rownames(out) <- NULL
  attr(out, "hops") <- do.call(rbind, hops)
  out
}

#' @keywords internal
.kegg_routes_empty <- function() {
  data.frame(route_id = character(0), pair_rank = integer(0), alternative = integer(0), source = character(0),
             source_symbol = character(0), compounds = character(0), n_compounds = integer(0),
             source_max_weight = numeric(0), destination = character(0), destination_symbol = character(0),
             destination_score = numeric(0), n_hops = integer(0), path = character(0), path_kegg = character(0),
             sign = character(0), pathways = character(0), pathway_titles = character(0), score = numeric(0),
             stringsAsFactors = FALSE)
}

#' @keywords internal
.kegg_routes_empty_hops <- function() {
  data.frame(route_id = character(0), hop = integer(0), from = character(0), to = character(0),
             from_symbol = character(0), to_symbol = character(0), effect = character(0), sign = integer(0),
             reversed = logical(0), pathways = character(0), stringsAsFactors = FALSE)
}

#' Compact table of routes: route -> hops -> pathways
#'
#' @description
#' One line per route of [kegg_routes()], for a report or the console:
#' the stops with KEGG's arrow notation between them (`-->` activation or
#' expression, `--|` inhibition or repression, `-p>` phosphorylation or
#' another modification of unstated sign, `..>` indirect effect, `---`
#' binding), the net sign and the pathways each hop comes from.
#'
#' @param routes Output of [kegg_routes()].
#' @param alternatives Logical, default `FALSE`: include the alternative
#'   routes (`R1b`, ...).
#' @return A `data.frame(route, hops, sign, route_path, compounds,
#'   pathways_per_hop)`.
#' @examples
#' \dontrun{
#' kegg_routes_table(kegg_routes(proj, condition = "FLO-ET", disease = "MONDO_0005044"))
#' }
#' @export
kegg_routes_table <- function(routes, alternatives = FALSE) {
  if (!is.data.frame(routes) || !all(c("route_id", "path", "n_hops") %in% names(routes))) {
    cli::cli_abort("{.arg routes} must be the output of {.fn kegg_routes}.")
  }
  .pathway_check_flag(alternatives, "alternatives")
  r <- if (alternatives) routes else routes[routes$alternative == 1, , drop = FALSE]
  hops <- attr(routes, "hops")
  per_hop <- vapply(r$route_id, function(id) {
    h <- hops[hops$route_id == id, , drop = FALSE]
    if (is.null(h) || nrow(h) == 0) return("")
    paste(sprintf("%s>%s [%s]", h$from_symbol, h$to_symbol, h$pathways), collapse = " | ")
  }, character(1))
  data.frame(route = r$route_id, hops = r$n_hops, sign = r$sign, route_path = r$path,
             compounds = r$n_compounds, pathways_per_hop = unname(per_hop), stringsAsFactors = FALSE,
             row.names = NULL)
}

#' Metro map of the routes from the extract's targets to disease genes
#'
#' @description
#' Draws the best routes of [kegg_routes()] like a transit map: each route
#' is a coloured line, genes are stations. Source stations (the extract's
#' targets) are filled and labelled with the compounds that hit them;
#' destination stations (disease genes) carry a red ring; a station used by
#' several routes is an interchange (larger, thick outline). Stations are
#' placed in columns from the sources (left) to the destinations (right)
#' (a layered, Sugiyama-type layout) and lines run horizontally, diagonally
#' at 45 degrees or vertically between them. Other KEGG relations among the
#' drawn genes are shown as faint grey links, and every hop carries its
#' effect mark (arrow = activation/expression, bar = inhibition/repression).
#'
#' @inheritParams kegg_routes
#' @inheritParams plot_save_params
#' @param routes `NULL` (default: computed with [kegg_routes()] from the
#'   other arguments) or an existing [kegg_routes()] result.
#' @param top_n_routes Integer, default `6`: routes drawn (best pairs,
#'   first alternative each unless `alternatives = TRUE`).
#' @param alternatives Logical, default `FALSE`.
#' @param background Logical, default `TRUE`: faint grey KEGG relations
#'   among the drawn stations.
#' @param max_compound_labels Integer, default `2`: compound names printed
#'   under a source station (the rest as "+n").
#' @param engine `"static"` (default) or `"ggiraph"`.
#' @param ... Passed to [kegg_routes()] when `routes` is `NULL`.
#'
#' @return A `ggplot`/`girafe` object with attribute `routes` (the drawn
#'   routes). With `save = TRUE` the PNG is written and logged in
#'   `patliRResults(proj, "kegg_routes_plot_log")` (keyed by `condition`,
#'   `target_set`), and the updated project is in `attr(., "proj")`.
#'
#' @examples
#' \dontrun{
#' plot_kegg_routes(proj, condition = "FLO-ET", disease = "MONDO_0005044")
#' }
#' @seealso [kegg_routes()], [kegg_routes_table()], [plot_kegg_topology()]
#' @export
plot_kegg_routes <- function(proj, condition = NULL, disease = NULL, to = NULL, routes = NULL, top_n_routes = 6,
                             alternatives = FALSE, background = TRUE, max_compound_labels = 2,
                             engine = c("static", "ggiraph"),
                             save = TRUE, out_dir = NULL, width = NULL, height = NULL, dpi = 150, ...) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  .pathway_check_count(top_n_routes, "top_n_routes", min = 1)
  .pathway_check_flag(alternatives, "alternatives")
  .pathway_check_flag(background, "background")
  .pathway_check_count(max_compound_labels, "max_compound_labels", min = 0)
  engine <- .plot_require(engine)
  if (is.null(routes)) {
    routes <- kegg_routes(proj, condition = condition, disease = disease, to = to,
                          top_n = max(top_n_routes, 1), ...)
  } else if (!is.data.frame(routes) || is.null(attr(routes, "hops"))) {
    cli::cli_abort("{.arg routes} must be the output of {.fn kegg_routes}.")
  }
  cond <- attr(routes, "condition") %||% condition
  disease <- attr(routes, "disease") %||% disease
  r <- if (alternatives) routes else routes[routes$alternative == 1, , drop = FALSE]
  r <- utils::head(r, if (alternatives) nrow(r) else top_n_routes)
  if (nrow(r) == 0) {
    cli::cli_abort(c("No route from the extract's targets to the destination genes within the length limit.",
                     "i" = "Increase {.arg max_len}, or use the full topology ({.code network_kegg_topology(restrict_to_network = FALSE)})."))
  }
  hops <- attr(routes, "hops")
  hops <- hops[hops$route_id %in% r$route_id, , drop = FALSE]

  topo_all <- patliRResults(proj, "network_kegg_topology")
  topo <- topo_all[topo_all$condition == cond, , drop = FALSE]
  lay <- .kegg_routes_layout(r, hops)
  st <- lay$stations
  cmp <- compounds(proj)
  cmp_name <- stats::setNames(ifelse(is.na(cmp$name) | !nzchar(cmp$name), cmp$id, cmp$name), cmp$id)

  pal <- .kegg_routes_palette(r$route_id)
  lines <- .kegg_routes_lines(r, hops, st)
  p <- ggplot2::ggplot()
  if (background) {
    bg <- .kegg_routes_background(topo, st, hops)
    if (nrow(bg)) {
      p <- p + ggplot2::geom_segment(data = bg, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                                     colour = "grey80", linewidth = 0.35, linetype = "22")
    }
  }
  p <- p +
    ggplot2::geom_path(data = lines$paths, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, colour = .data$route_id),
                       linewidth = 2.2, lineend = "round", linejoin = "round") +
    ggplot2::scale_colour_manual(values = pal, breaks = r$route_id, labels = .kegg_routes_legend_labels(r),
                                 name = "Routes (hops; pathways)",
                                 guide = ggplot2::guide_legend(ncol = 1, override.aes = list(linewidth = 2.2)))
  if (nrow(lines$marks)) {
    tee <- lines$marks[lines$marks$type == "tee", , drop = FALSE]
    arr <- lines$marks[lines$marks$type == "arrow", , drop = FALSE]
    if (nrow(arr)) {
      p <- p + ggplot2::geom_segment(data = arr, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                                     colour = "grey15", linewidth = 0.5,
                                     arrow = ggplot2::arrow(length = ggplot2::unit(0.07, "inches"), type = "closed", angle = 25))
    }
    if (nrow(tee)) {
      p <- p + ggplot2::geom_segment(data = tee, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                                     colour = "grey15", linewidth = 1.1)
    }
  }
  st$tooltip <- st$symbol
  st$name <- st$kegg
  inter <- st[st$n_routes > 1 & st$role == "via", , drop = FALSE]
  plain <- st[st$n_routes <= 1 & st$role == "via", , drop = FALSE]
  srcs <- st[st$role == "source", , drop = FALSE]
  dsts <- st[st$role == "destination", , drop = FALSE]
  if (nrow(plain)) p <- p + ggplot2::geom_point(data = plain, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = 3.4,
                                                fill = "white", colour = "grey15", stroke = 1)
  if (nrow(inter)) p <- p + ggplot2::geom_point(data = inter, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = 5,
                                                fill = "white", colour = "grey5", stroke = 1.8)
  if (nrow(dsts)) {
    p <- p +
      ggplot2::geom_point(data = dsts, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = 6.6,
                          fill = "white", colour = "#C0392B", stroke = 1.8) +
      ggplot2::geom_point(data = dsts, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = 3.2,
                          fill = "white", colour = "grey10", stroke = 1)
  }
  if (nrow(srcs)) p <- p + ggplot2::geom_point(data = srcs, ggplot2::aes(x = .data$x, y = .data$y), shape = 21, size = 5.4,
                                               fill = "grey10", colour = "grey10", stroke = 1)
  ## station names above, compound names under the sources
  st$label_y <- st$y + 0.42
  p <- p + ggplot2::geom_label(data = st, ggplot2::aes(x = .data$x, y = .data$label_y, label = .data$symbol),
                               size = 3, fontface = "bold", fill = "white", border.colour = NA,
                               label.padding = ggplot2::unit(0.08, "lines"), colour = "grey10")
  if (nrow(srcs) && max_compound_labels > 0) {
    srcs$cmp_label <- vapply(srcs$compounds, function(s) {
      ids <- strsplit(s, ";", fixed = TRUE)[[1]]
      nm <- .plot_truncate(unname(ifelse(is.na(cmp_name[ids]), ids, cmp_name[ids])), 24)
      shown <- utils::head(nm, max_compound_labels)
      paste0(paste(shown, collapse = "\n"), if (length(nm) > length(shown)) paste0("\n+", length(nm) - length(shown), " more") else "")
    }, character(1))
    p <- p + ggplot2::geom_text(data = srcs, ggplot2::aes(x = .data$x, y = .data$y - 0.38, label = .data$cmp_label),
                                size = 2.3, fontface = "italic", colour = "#0B4F7C", vjust = 1, lineheight = 0.85)
  }
  ## legend for the station types (dummy, keyed by shape)
  key <- data.frame(x = NA_real_, y = NA_real_, kind = factor(c("source", "via", "interchange", "destination"),
                                                              levels = c("source", "via", "interchange", "destination")))
  p <- p + ggplot2::geom_point(data = key, ggplot2::aes(x = .data$x, y = .data$y, shape = .data$kind), na.rm = TRUE) +
    ggplot2::scale_shape_manual(
      values = c(source = 21, via = 21, interchange = 21, destination = 21), drop = FALSE, name = "Stations",
      labels = c(source = "extract target (compounds below)", via = "intermediate gene",
                 interchange = "interchange (several routes)",
                 destination = if (!is.null(disease)) "disease gene" else "destination gene"),
      guide = ggplot2::guide_legend(order = 2, override.aes = list(
        size = c(4.4, 3, 4.2, 5), fill = c("grey10", "white", "white", "white"),
        colour = c("grey10", "grey15", "grey5", "#C0392B"), stroke = c(1, 1, 1.8, 1.8)))
    )

  gs <- attr(routes, "graph_size")
  dname <- .kegg_routes_disease_name(proj, disease)
  n_direct <- nrow(attr(routes, "direct_hits") %||% data.frame())
  sub <- paste0(
    "Routes along KEGG relations (all pathways of ", cond, " merged: ", gs[["genes"]], " genes, ", gs[["relations"]],
    " relations) from the extract's targets (left) to ", if (!is.null(dname)) paste0(dname, " genes") else "the destination genes",
    " (right). Each coloured line is one shortest route; arrowheads mark activation/expression hops, bars inhibition/repression, ",
    "plain joins are unsigned (phosphorylation, binding, indirect). Ranked by source support (compounds, probability), ",
    "destination association and length. Grey dashed: other KEGG relations among these genes.",
    if (n_direct > 0) sprintf(" %d target%s already disease gene%s (no route needed).", n_direct,
                              if (n_direct == 1) " is" else "s are", if (n_direct == 1) "" else "s") else ""
  )
  xr <- range(c(st$x, lines$paths$x))
  yr <- range(c(st$y, lines$paths$y))
  ## figure size from the map's own aspect (coord_equal): u inches per map
  ## unit, shrunk when the map is large; + legend column and title block
  data_w <- diff(xr) + 2.4
  data_h <- diff(yr) + 2.0
  u <- min(0.62, 11 / data_w, 9 / data_h)
  if (is.null(width)) width <- max(8, data_w * u + 5)
  if (is.null(height)) height <- max(5.5, data_h * u + 2.6, 0.62 * nrow(r) + 3.4)
  title <- paste0("Routes from the extract's targets to ", if (!is.null(dname)) dname else "the destination genes", " -- ", cond)
  p <- p +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(1.2, 1.2))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(1.2, 0.8))) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(title = .plot_wrap(title, .plot_wrap_width(width, 12)), subtitle = .plot_wrap(sub, .plot_wrap_width(width, 8))) +
    .network_view_theme("right") +
    ggplot2::theme(legend.key.width = ggplot2::unit(0.9, "cm"),
                   legend.key.spacing.y = ggplot2::unit(0.18, "cm"),
                   legend.text = ggplot2::element_text(size = 7.5, lineheight = 1))
  attr(p, "routes") <- r

  target_set <- if (!is.null(disease)) disease else paste0("custom_", substr(rlang::hash(sort(unique(r$destination))), 1, 8))
  log_row <- data.frame(condition = cond, target_set = target_set, n_routes = nrow(r),
                        max_hops = max(r$n_hops), path = NA_character_, stringsAsFactors = FALSE)
  .plot_finish(
    proj, p, name = "kegg_routes_plot_log",
    filename = paste0("kegg_routes_", cond, "_", gsub("[^A-Za-z0-9]+", "-", target_set), ".png"),
    log_row = log_row, key_cols = c("condition", "target_set"),
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Disease display name from `disease_genes`
#' @keywords internal
.kegg_routes_disease_name <- function(proj, disease) {
  if (is.null(disease)) return(NULL)
  dg <- patliRResults(proj, "disease_genes")
  nm <- if (!is.null(dg)) dg$disease_name[dg$disease_id == disease & !is.na(dg$disease_name)][1] else NA
  if (is.na(nm) || !nzchar(nm)) disease else nm
}

#' One colour per route (Okabe-Ito, then `hcl.colors`)
#' @keywords internal
.kegg_routes_palette <- function(ids) {
  oi <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#999933", "#882255")
  cols <- if (length(ids) <= length(oi)) oi[seq_along(ids)] else grDevices::hcl.colors(length(ids), "Dark 3")
  stats::setNames(cols, ids)
}

#' Legend text of the routes: "R1  EGFR > ... > NOS3 (3 hops, +; MAPK, Calcium)"
#' @keywords internal
.kegg_routes_legend_labels <- function(r) {
  short_pw <- vapply(strsplit(r$pathway_titles, "; ", fixed = TRUE), function(t) {
    t <- sub(" signaling pathway$", "", t)
    t <- .plot_truncate(t, 22)
    paste0(paste(utils::head(t, 2), collapse = ", "), if (length(t) > 2) paste0(" +", length(t) - 2) else "")
  }, character(1))
  path <- gsub(" (-->|--\\||\\.\\.>|-p>|---|-\\?-) ", " \\1 ", r$path)
  .plot_wrap(sprintf("%s  %s\n(%d hop%s, sign %s; %s)", r$route_id, path, r$n_hops, ifelse(r$n_hops == 1, "", "s"),
                     r$sign, short_pw), 46)
}

#' Station positions of the metro map
#'
#' @description
#' The union of the drawn routes is layered: every station at its longest
#' distance from a route start along the route edges (so sources that no
#' route passes through are in column 0), destinations that no route
#' continues from pulled to the last column. Rows inside a column come
#' from [igraph::layout_with_sugiyama()] (crossing reduction), so the map
#' reads left to right. Columns are 3 units apart, rows 1.6.
#' @return `list(stations = data.frame(kegg, symbol, role, n_routes,
#'   compounds, x, y))`.
#' @keywords internal
.kegg_routes_layout <- function(r, hops) {
  ed <- unique(hops[, c("from", "to")])
  nodes <- unique(c(r$source, hops$from, hops$to, r$destination))
  g <- igraph::graph_from_data_frame(ed, directed = TRUE, vertices = data.frame(name = nodes, stringsAsFactors = FALSE))
  if (!igraph::is_dag(g)) {
    fas <- igraph::feedback_arc_set(g, algo = "approx_eades")
    g <- igraph::delete_edges(g, fas)
  }
  topo_order <- igraph::topo_sort(g, mode = "out")
  layer <- stats::setNames(rep(0L, length(nodes)), nodes)
  vn <- igraph::V(g)$name
  for (v in vn[as.integer(topo_order)]) {
    succ <- vn[as.integer(igraph::neighbors(g, v, mode = "out"))]
    for (s in succ) layer[s] <- max(layer[s], layer[v] + 1L)
  }
  ## destinations that no route continues from go to the last column; a
  ## gene that is both a source/destination and a stop of another route
  ## keeps its place along the line
  dest <- unique(r$destination)
  sinks <- dest[igraph::degree(g, dest, mode = "out") == 0]
  layer[sinks] <- max(layer)
  s <- igraph::layout_with_sugiyama(g, layers = unname(layer[vn]), hgap = 1, vgap = 1, maxiter = 200)
  pos <- s$layout[seq_along(vn), , drop = FALSE]
  sym <- stats::setNames(c(r$source_symbol, r$destination_symbol, hops$from_symbol, hops$to_symbol),
                         c(r$source, r$destination, hops$from, hops$to))
  n_routes <- vapply(vn, function(v) length(unique(hops$route_id[hops$from == v | hops$to == v])), integer(1))
  role <- ifelse(vn %in% r$source, "source", ifelse(vn %in% dest, "destination", "via"))
  cmp <- stats::setNames(r$compounds, r$source)
  st <- data.frame(kegg = vn, symbol = unname(sym[vn]), role = role, n_routes = n_routes,
                   compounds = ifelse(role == "source", unname(cmp[vn]), ""),
                   x = layer[vn] * 3, y = -pos[, 1] * 1.6, stringsAsFactors = FALSE)
  st$y <- st$y - min(st$y)
  list(stations = st)
}

#' Octilinear polyline between two stations: horizontal, 45-degree
#' diagonal, horizontal (or diagonal, vertical, diagonal when the rise is
#' larger than the run)
#' @return Two-column matrix of points.
#' @keywords internal
.kegg_routes_octilinear <- function(x1, y1, x2, y2) {
  dx <- x2 - x1
  dy <- y2 - y1
  if (abs(dy) < 1e-9 || abs(dx) < 1e-9) return(rbind(c(x1, y1), c(x2, y2)))
  sx <- sign(dx)
  sy <- sign(dy)
  if (abs(dy) <= abs(dx)) {
    xm <- (x1 + x2) / 2
    rbind(c(x1, y1), c(xm - sx * abs(dy) / 2, y1), c(xm + sx * abs(dy) / 2, y2), c(x2, y2))
  } else {
    h <- abs(dx) / 2
    rbind(c(x1, y1), c(x1 + sx * h, y1 + sy * h), c(x1 + sx * h, y2 - sy * h), c(x2, y2))
  }
}

#' Route polylines (parallel offsets on shared hops) and hop effect marks
#' @return `list(paths = data.frame(route_id, group, x, y), marks =
#'   data.frame(type, x, y, xend, yend))`.
#' @keywords internal
.kegg_routes_lines <- function(r, hops, st) {
  pos <- stats::setNames(seq_len(nrow(st)), st$kegg)
  hk <- paste(pmin(hops$from, hops$to), pmax(hops$from, hops$to))
  share <- tapply(hops$route_id, hk, function(v) unique(v))
  off_step <- 0.2
  paths <- list()
  marks <- list()
  for (i in seq_len(nrow(hops))) {
    a <- pos[hops$from[i]]
    b <- pos[hops$to[i]]
    pts <- .kegg_routes_octilinear(st$x[a], st$y[a], st$x[b], st$y[b])
    users <- share[[hk[i]]]
    k <- match(hops$route_id[i], users)
    off <- (k - (length(users) + 1) / 2) * off_step
    pts[, 2] <- pts[, 2] + off
    paths[[i]] <- data.frame(route_id = hops$route_id[i], group = paste(hops$route_id[i], hops$hop[i]),
                             x = pts[, 1], y = pts[, 2], stringsAsFactors = FALSE)
    ## effect mark at 70% of the last segment of the hop
    m <- nrow(pts)
    p0 <- pts[m - 1, ]
    p1 <- pts[m, ]
    L <- sqrt(sum((p1 - p0)^2))
    if (L < 1e-9 || is.na(hops$sign[i])) next
    u <- (p1 - p0) / L
    tip <- p1 - u * 0.55
    if (hops$sign[i] > 0) {
      marks[[length(marks) + 1]] <- data.frame(type = "arrow", x = tip[1] - u[1] * 0.3, y = tip[2] - u[2] * 0.3,
                                               xend = tip[1], yend = tip[2])
    } else {
      marks[[length(marks) + 1]] <- data.frame(type = "tee", x = tip[1] - u[2] * 0.2, y = tip[2] + u[1] * 0.2,
                                               xend = tip[1] + u[2] * 0.2, yend = tip[2] - u[1] * 0.2)
    }
  }
  list(paths = do.call(rbind, paths),
       marks = if (length(marks)) do.call(rbind, marks) else
         data.frame(type = character(0), x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0)))
}

#' Other KEGG relations among the drawn stations (straight, for the
#' faint background)
#' @keywords internal
.kegg_routes_background <- function(topo, st, hops) {
  t <- topo[topo$from_kegg %in% st$kegg & topo$to_kegg %in% st$kegg & topo$from_kegg != topo$to_kegg, , drop = FALSE]
  if (nrow(t) == 0) return(data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0)))
  k <- unique(data.frame(a = pmin(t$from_kegg, t$to_kegg), b = pmax(t$from_kegg, t$to_kegg), stringsAsFactors = FALSE))
  used <- paste(pmin(hops$from, hops$to), pmax(hops$from, hops$to))
  k <- k[!paste(k$a, k$b) %in% used, , drop = FALSE]
  i <- match(k$a, st$kegg)
  j <- match(k$b, st$kegg)
  data.frame(x = st$x[i], y = st$y[i], xend = st$x[j], yend = st$y[j])
}
