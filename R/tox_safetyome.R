#' @include AllGenerics.R internal.R
NULL

#' Systemic safety-liability screen for predicted targets (Safetyome core panel) -- local, no network
#'
#' @description
#' Flags each predicted protein target (from [targets_import()]/
#' [targets_import_batch()]) against the "core panel" of genes identified
#' by Liu, X. et al. (2026), "Safetyome and specialized panels for over
#' 3,000 phenotypes: a systematic and translational approach using human
#' genetics and pharmacology", *Toxicological Sciences* 209(3), kfag021,
#' \doi{10.1093/toxsci/kfag021}.
#'
#' Unlike [tox_local()] (structural SMARTS alerts on the *compound*), this
#' is annotation at the *target* level: given a predicted `uniprot_id`,
#' mark whether it is one of the ~500 genes the paper's own prioritization
#' pipeline ranked as carrying the strongest cross-phenotype systemic
#' safety signal (tissue-specificity + evolutionary conservation +
#' human-genetics evidence density, across 22 MedDRA System Organ
#' Classes), and for which organ system(s).
#'
#' @section This never produces a pass/fail verdict either:
#' Same reasoning as [tox_local()] -- a target in the core panel is not
#' automatically a reason to drop a compound; it is a flag worth reading
#' in context (see `tox_report()`'s disclaimer, which applies here too).
#'
#' @section Coverage -- the 500-gene core panel, not the full ~11,000-gene catalog:
#' `inst/extdata/safetyome_core_panel.csv` is transcribed directly from
#' Supplementary Table 4 of Liu et al. (2026) ("Core gene (top 500) panels
#' selected by different metrics") -- 500 unique genes (507 rows; 7 genes
#' are listed twice under two different organ systems, kept as separate
#' rows rather than collapsed). The paper's full prioritized catalog
#' (Supplementary Table 2, ~11,300 genes across every score tier) is
#' **not** bundled -- it is the entire genome-scale ranking, not a curated
#' alert set, and at that size it stops being a targeted screen and starts
#' being "does this gene exist" for most real target lists. If broader
#' coverage turns out to matter in practice, extending to Supplementary
#' Table 2 is a config change, not a redesign.
#'
#' @section UniProt -> gene symbol mapping:
#' The Safetyome panel is keyed by HGNC gene symbol; `targets_imported` is
#' keyed by UniProt accession. Mapped via
#' `clusterProfiler::bitr(..., fromType = "UNIPROT", toType = "SYMBOL",
#' OrgDb = "org.Hs.eg.db")` -- the same mapping mechanism already used by
#' [network_enrich()] for UniProt -> Entrez, applied here to UniProt ->
#' SYMBOL instead, so no new ID-mapping dependency is introduced. UniProt
#' IDs with no symbol in `org.Hs.eg.db` are excluded and logged rather
#' than silently dropped.
#'
#' @param proj A [`PatliRProject-class`] object with a non-empty
#'   `targets_imported` entry in [patliRResults()] (run [targets_import()]
#'   or [targets_import_batch()] first).
#' @param compound_ids Character vector of `compounds(proj)$id` to restrict
#'   to, or `NULL` (default) for every compound with imported targets.
#'
#' @return The updated `proj`, with a `tox_safetyome` entry in
#'   [patliRResults()] (columns `compound_id`, `uniprot_id`,
#'   `gene_symbol`, `in_core_panel`, `organ_system`,
#'   `safetyome_scaled_score`, `tau_score_percentile`,
#'   `conservation_score_percentile`, `safetyome_scaled_score_percentile`,
#'   `median_tau_conservation_score`, `median_all_scores`, `source`), also
#'   written to `results/tox_safetyome.csv`. Every predicted target gets a
#'   row (`in_core_panel = FALSE` and score columns `NA` when it is not in
#'   the panel) -- absence of a flag is always visible in the table, never
#'   just an omitted row.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- targets_import_batch(
#'   proj,
#'   system.file("extdata", "import_targets", package = "patliR"),
#'   platform = "superpred"
#' )
#' proj <- tox_safetyome(proj)
#' patliRResults(proj, "tox_safetyome")
#' }
#'
#' @export
tox_safetyome <- function(proj, compound_ids = NULL) {
  stopifnot(is(proj, "PatliRProject"))

  tgt <- patliRResults(proj, "targets_imported")
  if (is.null(tgt) || nrow(tgt) == 0) {
    cli::cli_abort(c(
      "No {.val targets_imported} found in {.arg proj}.",
      "i" = "Run {.fn targets_import} or {.fn targets_import_batch} first."
    ))
  }
  if (!is.null(compound_ids)) {
    tgt <- tgt[tgt$compound_id %in% compound_ids, , drop = FALSE]
  }
  if (nrow(tgt) == 0) {
    cli::cli_abort("No rows in {.val targets_imported} match the requested {.arg compound_ids}.")
  }

  needed <- c("clusterProfiler", "org.Hs.eg.db")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "{.pkg {missing}} required for {.fn tox_safetyome} (UniProt -> gene symbol mapping).",
      "i" = "See {.file TESTING_GUIDE.Rmd}, section 0, for the {.fn BiocManager::install} chunk."
    ))
  }

  panel <- .load_safetyome_core_panel()

  uids <- unique(tgt$uniprot_id)
  map <- clusterProfiler::bitr(
    uids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = TRUE
  )

  unmapped <- setdiff(uids, map$UNIPROT)
  if (length(unmapped) > 0) {
    proj <- .log_append(
      proj, step = "tox_safetyome", id = unmapped,
      message = "tox_safetyome_unmapped: no gene symbol found in org.Hs.eg.db for this UniProt ID; excluded from the safetyome match"
    )
  }

  tgt$gene_symbol <- map$SYMBOL[match(tgt$uniprot_id, map$UNIPROT)]
  tgt <- tgt[!is.na(tgt$gene_symbol), , drop = FALSE]

  hits <- merge(
    tgt[, c("compound_id", "uniprot_id", "gene_symbol")],
    panel,
    by.x = "gene_symbol", by.y = "gene", all.x = TRUE, sort = FALSE
  )
  hits$in_core_panel <- !is.na(hits$safetyome_scaled_score)
  hits$source <- "safetyome_2026"

  out <- hits[, c(
    "compound_id", "uniprot_id", "gene_symbol", "in_core_panel", "organ_system",
    "safetyome_scaled_score", "tau_score_percentile", "conservation_score_percentile",
    "safetyome_scaled_score_percentile", "median_tau_conservation_score", "median_all_scores",
    "source"
  )]

  out <- .network_upsert(proj, "tox_safetyome", out, "compound_id")
  patliRResults(proj, "tox_safetyome") <- out
  .write_results_csv(proj, "tox_safetyome", out)
  .write_log_csv(proj)
  proj
}

#' Load the bundled Safetyome core panel (500 genes, Liu et al. 2026)
#'
#' See the "Coverage" section of [tox_safetyome()] for provenance.
#' @return `data.frame(gene, organ_system, evidence_sources,
#'   safetyome_scaled_score, tau_score_percentile,
#'   conservation_score_percentile, safetyome_scaled_score_percentile,
#'   median_tau_conservation_score, median_all_scores)`.
#' @keywords internal
.load_safetyome_core_panel <- function() {
  path <- system.file("extdata", "safetyome_core_panel.csv", package = "patliR")
  if (!nzchar(path) || !file.exists(path)) {
    cli::cli_abort(c(
      "{.val inst/extdata/safetyome_core_panel.csv} not found in this {.pkg patliR} install.",
      "i" = "Re-install the package, or see {.file TESTING_GUIDE.Rmd} for provenance/regeneration."
    ))
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}
