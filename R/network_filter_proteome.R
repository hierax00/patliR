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
#' @param proteome_label `NULL` (default) auto-derives a stable label from a
#'   cheap hash of `sort(unique(proteome))` (so the same set of accessions,
#'   in any order, always gets the same label); or a single non-empty
#'   character string to name this proteome yourself (e.g. `"liver"`,
#'   `"mitochondrion"`). Without *some* label, `network_filtered_edges.csv`
#'   is uninterpretable on its own (which proteome produced which rows) --
#'   against the CSV-is-truth contract. Part of the upsert key alongside
#'   `condition`, so filtering the same condition against two different
#'   proteomes keeps both results instead of one overwriting the other.
#'
#' @return The updated `proj`, with a `network_filtered_edges` entry in
#'   [patliRResults()] (same columns as `network_edges` -- `condition`,
#'   `compound_id`, `uniprot_id`, `weight` -- plus `proteome_label` and
#'   `n_proteome` (`length(proteome)`, the size of the proteome that
#'   produced this row set)), also written to
#'   `results/network_filtered_edges.csv`. How many of each condition's
#'   edges were kept is always logged, even when the answer is zero. A
#'   `proteome` that matches nothing in `network_edges` AND does not look
#'   like UniProt accessions (see the UniProt accession regex in the
#'   source) gets a `cli_warn` -- the most common cause is passing gene
#'   symbols by mistake.
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
network_filter_proteome <- function(proj, proteome, condition = NULL, proteome_label = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(proteome), length(proteome) > 0)
  if (is.null(proteome_label)) {
    proteome_label <- .network_proteome_auto_label(proteome)
  } else {
    stopifnot(is.character(proteome_label), length(proteome_label) == 1, nzchar(proteome_label))
  }

  conditions <- .network_resolve_conditions(proj, condition)
  edges_all <- patliRResults(proj, "network_edges")
  edges <- edges_all[edges_all$condition %in% conditions, , drop = FALSE]

  ## Sanity check: gene symbols (or any other non-UniProt identifier) as
  ## input give a silent zero-row result otherwise, with no clue why. Only
  ## warn when BOTH signals point the same way -- nothing matched at all,
  ## and the input doesn't even look like UniProt accessions -- so a
  ## deliberately narrow, correctly-spelled proteome that just happens to
  ## share nothing with this condition's targets does not get a spurious
  ## warning.
  uniprot_shaped <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$", proteome)
  if (length(intersect(proteome, edges$uniprot_id)) == 0 && mean(uniprot_shaped) < 0.5) {
    cli::cli_warn(c(
      "None of {.arg proteome} matched any {.field uniprot_id} in {.val network_edges}, and fewer than half of its {length(proteome)} entries look like UniProt accessions.",
      "i" = "{.fn network_filter_proteome} expects UniProt accessions (e.g. {.val P08183}), not gene symbols -- double check {.arg proteome}."
    ))
  }

  new_filtered <- edges[edges$uniprot_id %in% proteome, , drop = FALSE]
  rownames(new_filtered) <- NULL
  ## rep(..., nrow(new_filtered)), not a bare scalar assignment: a
  ## zero-row new_filtered (nothing matched) rejects `df$col <- scalar`
  ## with "replacement has 1 row, data has 0" -- exactly the zero-match
  ## case this function must handle cleanly (spec 1.11).
  new_filtered$proteome_label <- rep(proteome_label, nrow(new_filtered))
  new_filtered$n_proteome <- rep(length(proteome), nrow(new_filtered))

  n_kept <- vapply(conditions, function(cond) sum(new_filtered$condition == cond), integer(1))
  n_total <- vapply(conditions, function(cond) sum(edges$condition == cond), integer(1))
  proj <- .log_append(
    proj, step = "network_filter_proteome", id = NA_character_,
    message = paste0(
      "condition '", conditions, "', proteome_label '", proteome_label, "': ", n_kept, " of ", n_total,
      " network_edges rows kept (uniprot_id %in% proteome, n_proteome=", length(proteome), ")"
    )
  )

  result <- .network_upsert(
    proj, "network_filtered_edges", new_filtered, c("condition", "proteome_label"),
    touched_keys = data.frame(condition = conditions, proteome_label = proteome_label, stringsAsFactors = FALSE)
  )
  patliRResults(proj, "network_filtered_edges") <- result
  .write_results_csv(proj, "network_filtered_edges", result)
  .write_log_csv(proj)
  proj
}

#' Stable auto-label for a proteome vector, from a cheap hash of its sorted,
#' deduplicated content -- so provenance always exists even when the caller
#' does not name one, and the same set of accessions (any order) round-trips
#' to the same label
#' @return `character(1)`.
#' @keywords internal
.network_proteome_auto_label <- function(proteome) {
  x <- paste(sort(unique(proteome)), collapse = "\x01")
  codes <- utf8ToInt(enc2utf8(x))
  if (length(codes) == 0) return("unnamed")
  h <- sum(as.numeric(codes) * seq_along(codes)) %% (2^31 - 1)
  paste0("proteome_", sprintf("%08x", as.integer(h)))
}
