#' @include AllGenerics.R internal.R network_build.R
NULL

## Deliberately narrow: a reporting/export filter over network_edges, not a
## new graph substrate. Wiring the filtered result into the rest of
## network_* (centrality, robustness, ...) would need a new scope axis
## alongside "condition" in .network_resolve_conditions()/.network_graph().
## See ROADMAP.md.

#' Filter a project's compound-target network to a proteome of interest
#'
#' @description
#' Restricts [patliRResults()]`(proj, "network_edges")` to rows whose
#' `uniprot_id` is in `proteome` -- e.g. the proteome of a specific cell
#' line, tissue, or organism/species of interest (FOPER 2026 proposal,
#' month 5) -- so only compounds and targets relevant to that biological
#' context remain visible. Never modifies `network_edges` itself (still
#' the full network every other `network_*` function operates on); this
#' produces a separate, inspectable/exportable filtered view alongside it.
#'
#' @inheritParams network_build
#' @param proteome Character vector of UniProt IDs defining the proteome
#'   of interest.
#'
#' @return The updated `proj`, with a `network_filtered_edges` entry in
#'   [patliRResults()] (same columns as `network_edges`: `condition`,
#'   `compound_id`, `uniprot_id`, `weight`, `disease_association_score`),
#'   also written to `results/network_filtered_edges.csv`. How many of
#'   each condition's edges were kept is always logged, even when the
#'   answer is zero.
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
#' proj <- network_filter_proteome(proj, proteome = c("P03372", "P08183"))
#' patliRResults(proj, "network_filtered_edges")
#' }
#'
#' @export
network_filter_proteome <- function(proj, proteome, condition = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(proteome), length(proteome) > 0)

  conditions <- .network_resolve_conditions(proj, condition)
  edges_all <- patliRResults(proj, "network_edges")
  edges <- edges_all[edges_all$condition %in% conditions, , drop = FALSE]

  new_filtered <- edges[edges$uniprot_id %in% proteome, , drop = FALSE]
  rownames(new_filtered) <- NULL

  n_kept <- vapply(conditions, function(cond) sum(new_filtered$condition == cond), integer(1))
  n_total <- vapply(conditions, function(cond) sum(edges$condition == cond), integer(1))
  proj <- .log_append(
    proj, step = "network_filter_proteome", id = NA_character_,
    message = paste0(
      "condition '", conditions, "': ", n_kept, " of ", n_total,
      " network_edges rows kept (uniprot_id %in% proteome, n_proteome=", length(proteome), ")"
    )
  )

  result <- .network_upsert(proj, "network_filtered_edges", new_filtered, "condition")
  patliRResults(proj, "network_filtered_edges") <- result
  .write_results_csv(proj, "network_filtered_edges", result)
  .write_log_csv(proj)
  proj
}
