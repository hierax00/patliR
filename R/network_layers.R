#' @include AllGenerics.R internal.R network_build.R
NULL

## Shared internal builder: the directed compound -> target -> {pathway,
## disease} layered graph that gives network_motifs() (and, once feasible,
## network_bowtie()) real directed-graph semantics instead of forcing a
## bipartite-undirected reinterpretation. Confirmed with Uriel (2026-07-23,
## chat) before building this: rather than fabricate a "proxy" motif/bowtie
## definition on the undirected compound-target graph, build a genuinely
## directed structure from data patliR already computes:
##
##   compound --acts on--> target --participates in--> pathway
##                              \--associated with-----> disease
##
## Two "derived" edge kinds (compound_pathway_derived, compound_disease_
## derived) are added too: if a compound reaches a pathway/disease only via
## one of its targets, that transitive relationship is made an explicit
## edge. This is not fabricating a new relationship -- it is exactly what
## "this compound affects this pathway" already means once you have
## compound->target->pathway -- but it is what turns a compound->target->
## pathway chain into an actual feed-forward-loop-eligible triangle
## (compound->target, compound->pathway, target->pathway), which is the
## real point of building this layer: network_motifs() needs triangles to
## find, not just a 3-hop chain.

#' Build the directed compound-target-pathway-disease layered graph for one
#' condition
#'
#' @description
#' Combines three already-computed tables -- [network_build()]'s
#' `network_edges` (compound -> target), [network_enrich()]'s
#' `network_enrichment` (target -> pathway, reconstructed from the
#' `geneID` column of significant pathways for this condition), and
#' [targets_disease_filter()]'s `targets_disease` (target -> disease) --
#' into one directed graph. Only `network_edges` is required; the pathway
#' and/or disease layers are silently omitted (not an error) if
#' [network_enrich()]/[targets_disease_filter()] have not been run, since
#' this graph is meant to reflect whatever is actually available, not force
#' every consumer to have completed the whole pipeline.
#'
#' @return An `igraph` directed graph. Every vertex has a `layer` attribute
#'   (`"compound"`, `"target"`, `"pathway"`, or `"disease"`); every edge has
#'   an `edge_kind` attribute (`"compound_target"`, `"target_pathway"`,
#'   `"target_disease"`, `"compound_pathway_derived"`,
#'   `"compound_disease_derived"`).
#' @param pathway_db See `.network_target_pathway_edges()`, internal.
#' @keywords internal
.network_layered_graph <- function(proj, condition, pathway_db = NULL) {
  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || !condition %in% edges_all$condition) {
    cli::cli_abort(c(
      "No network built for condition {.val {condition}}.",
      "i" = "Run {.fn network_build} first."
    ))
  }

  ct <- unique(edges_all[edges_all$condition == condition, c("compound_id", "uniprot_id")])
  edges <- data.frame(
    from = ct$compound_id, to = ct$uniprot_id,
    from_layer = "compound", to_layer = "target", edge_kind = "compound_target",
    stringsAsFactors = FALSE
  )

  target_pathway <- .network_target_pathway_edges(proj, condition, unique(ct$uniprot_id), pathway_db = pathway_db)
  if (!is.null(target_pathway) && nrow(target_pathway) > 0) {
    edges <- rbind(edges, data.frame(
      from = target_pathway$uniprot_id, to = target_pathway$pathway_id,
      from_layer = "target", to_layer = "pathway", edge_kind = "target_pathway",
      stringsAsFactors = FALSE
    ))
    cp <- unique(merge(ct, target_pathway, by = "uniprot_id")[, c("compound_id", "pathway_id")])
    if (nrow(cp) > 0) {
      edges <- rbind(edges, data.frame(
        from = cp$compound_id, to = cp$pathway_id,
        from_layer = "compound", to_layer = "pathway", edge_kind = "compound_pathway_derived",
        stringsAsFactors = FALSE
      ))
    }
  }

  disease_all <- patliRResults(proj, "targets_disease")
  if (!is.null(disease_all) && nrow(disease_all) > 0) {
    td <- unique(disease_all[disease_all$target_id %in% ct$uniprot_id, c("target_id", "disease_id")])
    if (nrow(td) > 0) {
      edges <- rbind(edges, data.frame(
        from = td$target_id, to = td$disease_id,
        from_layer = "target", to_layer = "disease", edge_kind = "target_disease",
        stringsAsFactors = FALSE
      ))
      cd <- unique(merge(ct, td, by.x = "uniprot_id", by.y = "target_id")[, c("compound_id", "disease_id")])
      if (nrow(cd) > 0) {
        edges <- rbind(edges, data.frame(
          from = cd$compound_id, to = cd$disease_id,
          from_layer = "compound", to_layer = "disease", edge_kind = "compound_disease_derived",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  vertices <- unique(rbind(
    data.frame(name = edges$from, layer = edges$from_layer, stringsAsFactors = FALSE),
    data.frame(name = edges$to, layer = edges$to_layer, stringsAsFactors = FALSE)
  ))
  igraph::graph_from_data_frame(d = edges[, c("from", "to", "edge_kind")], directed = TRUE, vertices = vertices)
}

#' Reconstruct target(UniProt) -> pathway edges from a condition's
#' significant `network_enrichment` rows
#'
#' @description
#' `network_enrich()` stores `geneID` as a `/`-separated list of *Entrez*
#' IDs per significant pathway (`clusterProfiler`'s own convention -- see
#' `network_enrich()`), so this maps the condition's target UniProt IDs to
#' Entrez the same way `network_enrich()`/`network_degeneracy()` do, then
#' checks which of a pathway's member Entrez IDs are among them.
#'
#' @param pathway_db Character vector of `network_enrichment$db` values to
#'   restrict to (any of `"go"`, `"reactome"`, `"kegg"`), or `NULL`
#'   (default) for every `db` [network_enrich()] has been run with for this
#'   condition. `network_enrich()` results accumulate across separate calls
#'   with different `db` values (see its own docs) -- without this filter,
#'   running `db = "go"` and `db = "kegg"` for the same condition silently
#'   mixes both into the pathway layer. Worth restricting to `"kegg"` alone
#'   in practice: KEGG's gene sets are far coarser/fewer than GO's, so a
#'   `compound -> pathway` layer built only from KEGG has a much smaller
#'   derived-edge fan-out per compound than one that includes GO -- a
#'   direct, no-tuning alternative to [network_enrich()]'s own `simplify_go`
#'   for keeping [network_motifs()]/[network_degeneracy()]/
#'   [plot_network_layers()] from turning into an unreadable hairball on a
#'   real, densely-annotated target set (confirmed against real Chilcuague
#'   data, 2026-08-25).
#' @return `data.frame(uniprot_id, pathway_id)`, or `NULL` if
#'   `network_enrich()` was never run for this condition (with a matching
#'   `db`, if `pathway_db` was given), or the mapping packages are not
#'   installed (same optional `Suggests` as [network_enrich()] --
#'   silently omitting the pathway layer here is consistent with this being
#'   a best-effort enrichment of the graph, not a hard requirement).
#' @keywords internal
.network_target_pathway_edges <- function(proj, condition, uniprot_ids, pathway_db = NULL) {
  enrichment_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enrichment_all) || length(uniprot_ids) == 0) return(NULL)
  enr <- enrichment_all[enrichment_all$condition == condition, , drop = FALSE]
  if (!is.null(pathway_db)) enr <- enr[enr$db %in% pathway_db, , drop = FALSE]
  if (nrow(enr) == 0) return(NULL)
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(NULL)
  }

  map <- tryCatch(
    clusterProfiler::bitr(uniprot_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE),
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) return(NULL)

  rows <- lapply(seq_len(nrow(enr)), function(i) {
    entrez_in_pathway <- strsplit(enr$geneID[i], "/")[[1]]
    hit_uniprot <- unique(map$UNIPROT[map$ENTREZID %in% entrez_in_pathway])
    if (length(hit_uniprot) == 0) return(NULL)
    data.frame(uniprot_id = hit_uniprot, pathway_id = enr$ID[i], stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}
