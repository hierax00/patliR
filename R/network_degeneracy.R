#' @include AllGenerics.R internal.R network_build.R network_layers.R
NULL

## "Degeneracy" in the systems-biology sense (Edelman & Gally, 2001, PNAS
## 98(24), 13763-8): two compounds are degenerate when they act on largely
## DIFFERENT targets but those targets converge on largely the SAME
## enriched pathways -- redundant routes to the same effect. A pairwise
## compound property. Unlike the rest of network_*, this has a hard
## dependency on network_enrich()'s output, not just network_build().

#' Pairwise functional degeneracy between compounds (pathway-based)
#'
#' @description
#' For each condition, computes, for every pair of compounds present in that
#' condition's network:
#'
#' - `target_jaccard`: Jaccard similarity of the two compounds' target sets
#'   (structural overlap -- do they act on the same proteins?).
#' - `pathway_jaccard`: Jaccard similarity of the two compounds' *enriched*
#'   pathway sets (functional overlap -- reconstructed from
#'   [network_enrich()]'s `geneID` column via the same UniProt -> Entrez
#'   mapping [network_enrich()] itself uses, see `.network_target_pathway_edges()`,
#'   internal).
#' - `degeneracy_score = pathway_jaccard * (1 - target_jaccard)`: high only
#'   when two compounds are functionally similar (`pathway_jaccard` near 1)
#'   *despite* being structurally distinct (`target_jaccard` near 0) --
#'   different elements, same function, in the spirit of degeneracy as
#'   Edelman & Gally (2001, *PNAS* 98(24), 13763-8) describe it
#'   qualitatively (their quantitative measure, Tononi et al. 1999, is a
#'   different, information-theoretic construction). A pair that hits the
#'   same targets is redundant, not degenerate; a pair that shares neither
#'   targets nor pathways is unrelated. Only the different-targets /
#'   same-pathways quadrant scores high.
#'
#'   The multiplicative form is a convenience scalarisation, not a derived
#'   result -- and in sparse target data `target_jaccard` is ~0 for almost
#'   every pair, so `degeneracy_score` collapses to roughly
#'   `pathway_jaccard`. Read `target_jaccard` and `pathway_jaccard`
#'   together (the columns are both in the output) rather than leaning on
#'   the scalar alone; `pathway_jaccard` is also an unweighted set overlap
#'   over GO terms, so it does not account for term specificity (see
#'   `ROADMAP.md` -- a `GOSemSim` semantic-similarity version is planned).
#'
#' @section Requires [network_enrich()] to have already run for the condition:
#' Unlike [network_centrality()] or [network_hub_penalty()] (which only need
#' [network_build()]), this function needs the pathway layer, so it errors
#' clearly (naming the exact missing prerequisite) rather than silently
#' returning empty/`NA` results if [network_enrich()] has not been run, or if
#' its optional package dependencies (`clusterProfiler`, `org.Hs.eg.db`) are
#' not installed.
#'
#' @section Compound pairs with no pathway data are omitted, not `NA`-filled:
#' If neither compound in a pair has any target participating in a
#' significant pathway (`pathway_jaccard` would be `0/0`), that pair
#' contributes no row -- consistent with the rest of `patliR`, "nothing to
#' compare" is not the same as "zero functional similarity". See
#' `projectLog(proj)` for how many pairs were skipped this way.
#'
#' @inheritParams network_build
#' @param pathway_db See `.network_target_pathway_edges()`, internal --
#'   restrict `pathway_jaccard`'s pathway sets to specific
#'   `network_enrich()` `db` value(s) (e.g. `"kegg"` alone), or `NULL`
#'   (default) for every `db` combined. Worth restricting when a real
#'   target set makes GO's pathway_jaccard uninformatively close to 1 for
#'   nearly every pair -- GO's hierarchical redundancy means most compounds
#'   end up sharing at least some very generic term, which flattens
#'   `degeneracy_score` toward "everything looks degenerate".
#'   `network_enrich()`'s own
#'   `simplify_go` addresses the same root cause upstream (fewer, less
#'   redundant GO terms in the first place); this argument is the
#'   complementary downstream lever, for when you want KEGG's coarser
#'   pathway sets specifically rather than a de-duplicated but still
#'   GO-scale set of terms.
#'
#' @return The updated `proj`, with a `network_degeneracy` entry in
#'   [patliRResults()] (columns `condition`, `compound_a`, `compound_b`,
#'   `n_targets_a`, `n_targets_b`, `n_pathways_a`, `n_pathways_b`,
#'   `target_jaccard`, `pathway_jaccard`, `degeneracy_score`), also written
#'   to `results/network_degeneracy.csv`.
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
#' proj <- network_enrich(proj, condition = "FLO-ET", db = "go")
#' proj <- network_degeneracy(proj, condition = "FLO-ET")
#' patliRResults(proj, "network_degeneracy")
#' }
#'
#' @export
network_degeneracy <- function(proj, condition = NULL, pathway_db = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  conditions <- .network_resolve_conditions(proj, condition)

  enrichment_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enrichment_all) || nrow(enrichment_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_enrichment} entry in {.arg proj}.",
      "i" = "{.fn network_degeneracy} needs the pathway layer -- run {.fn network_enrich} first."
    ))
  }
  missing_enrich <- setdiff(conditions, unique(enrichment_all$condition))
  if (length(missing_enrich) > 0) {
    cli::cli_abort(c(
      "Condition(s) {.val {missing_enrich}} have no {.fn network_enrich} results.",
      "i" = "Run {.fn network_enrich} for {.val {missing_enrich}} first."
    ))
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_degeneracy} needs {.pkg clusterProfiler} and {.pkg org.Hs.eg.db} (same UniProt -> Entrez mapping as {.fn network_enrich}), not installed.",
      "i" = "Install them with {.code BiocManager::install(c(\"clusterProfiler\", \"org.Hs.eg.db\"))}."
    ))
  }

  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions
  n_skipped_total <- 0L

  for (cond in conditions) {
    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    all_compounds <- unique(ct$compound_id)

    if (length(all_compounds) < 2) {
      rows[[cond]] <- .empty_network_degeneracy_row()
      proj <- .log_append(
        proj, step = "network_degeneracy", id = NA_character_,
        message = paste0("condition '", cond, "': fewer than 2 compounds in the network, nothing to compare")
      )
      next
    }

    target_pathway <- .network_target_pathway_edges(proj, cond, unique(ct$uniprot_id), pathway_db = pathway_db)
    cp <- if (!is.null(target_pathway) && nrow(target_pathway) > 0) {
      unique(merge(ct, target_pathway, by = "uniprot_id")[, c("compound_id", "pathway_id")])
    } else {
      data.frame(compound_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE)
    }

    target_sets <- split(ct$uniprot_id, ct$compound_id)
    pathway_sets <- split(cp$pathway_id, cp$compound_id)

    pairs <- utils::combn(sort(all_compounds), 2, simplify = FALSE)
    pair_rows <- vector("list", length(pairs))
    n_skipped <- 0L

    for (i in seq_along(pairs)) {
      a <- pairs[[i]][1]; b <- pairs[[i]][2]
      ta <- target_sets[[a]]; tb <- target_sets[[b]]
      pa <- pathway_sets[[a]]; pb <- pathway_sets[[b]]
      if (is.null(pa)) pa <- character(0)
      if (is.null(pb)) pb <- character(0)

      pathway_union <- length(union(pa, pb))
      if (pathway_union == 0) {
        n_skipped <- n_skipped + 1L
        pair_rows[[i]] <- NULL
        next
      }
      pathway_jaccard <- length(intersect(pa, pb)) / pathway_union
      target_jaccard <- length(intersect(ta, tb)) / length(union(ta, tb))
      degeneracy_score <- pathway_jaccard * (1 - target_jaccard)

      pair_rows[[i]] <- data.frame(
        condition = cond, compound_a = a, compound_b = b,
        n_targets_a = length(ta), n_targets_b = length(tb),
        n_pathways_a = length(pa), n_pathways_b = length(pb),
        target_jaccard = target_jaccard, pathway_jaccard = pathway_jaccard,
        degeneracy_score = degeneracy_score, stringsAsFactors = FALSE
      )
    }

    pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
    rows[[cond]] <- if (length(pair_rows) == 0) .empty_network_degeneracy_row() else do.call(rbind, pair_rows)
    n_skipped_total <- n_skipped_total + n_skipped

    proj <- .log_append(
      proj, step = "network_degeneracy", id = NA_character_,
      message = paste0(
        "condition '", cond, "': ", nrow(rows[[cond]]), " compound pair(s) scored, ",
        n_skipped, " pair(s) skipped (neither compound has a target in a significant pathway)"
      )
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_degeneracy", result, "condition")

  patliRResults(proj, "network_degeneracy") <- result
  .write_results_csv(proj, "network_degeneracy", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.empty_network_degeneracy_row <- function() {
  data.frame(
    condition = character(0), compound_a = character(0), compound_b = character(0),
    n_targets_a = integer(0), n_targets_b = integer(0),
    n_pathways_a = integer(0), n_pathways_b = integer(0),
    target_jaccard = double(0), pathway_jaccard = double(0), degeneracy_score = double(0),
    stringsAsFactors = FALSE
  )
}
