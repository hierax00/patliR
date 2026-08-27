#' @include AllGenerics.R internal.R
NULL

#' Import toxicity properties exported from an external platform
#'
#' @description
#' Imports a CSV exported from a third-party toxicity platform (SwissADME,
#' ADMETlab, pkCSM, ...) -- hERG liability, hepatotoxicity (DILI), Ames
#' mutagenicity, and similar model-derived predictions that are not
#' calculable locally in R (see [tox_local()] for what *is* local: PAINS/
#' Brenk structural alerts). Same import philosophy as [adme_import()]:
#' `patliR` never scrapes these platforms, you download the CSV yourself,
#' and this function reconciles it against [compounds()].
#'
#' @inheritParams adme_import
#' @param mapping_file Path to the `<out_file>_map.csv` written by
#'   [tox_export_smiles()], or `NULL` (default). When given, rows are
#'   matched to `compound_id` **by position** against that file (the
#'   platform returns its output in input order), not by SMILES -- this is
#'   the reliable way to reconcile an ADMETlab/pkCSM export, since a
#'   third-party toolkit's canonical SMILES is not guaranteed to match
#'   patliR's (see `DESIGN.md`). Takes precedence over the `platform`
#'   matching rule.
#'
#' @details
#' Built-in presets (confirmed against the column headers the bundled
#' example file actually uses):
#' \itemize{
#'   \item `platform = "admetlab"`: `smiles` (canonicalized the same way
#'     [prep_compounds()] does, matched against
#'     `compounds(proj)$canonical_smiles` unless `mapping_file` is given),
#'     `Ames`, `DILI`, `hERG` (all mapped to lowercase property names).
#'   \item `platform = "swissadme"`: no built-in preset -- the bundled
#'     `import_adme_swissadme.csv` example only carries ADME properties, not
#'     toxicity ones; pass your own `column_map` if your SwissADME export
#'     includes toxicity-relevant columns.
#' }
#'
#' @return The updated `proj`, with a `tox_imported` entry in
#'   [patliRResults()] in long format (columns `compound_id`, `property`,
#'   `value`, `source`, `import_date`), also written to
#'   `results/tox_imported.csv`. Rows that could not be matched to a known
#'   compound are logged (`"tox_import_unmatched"`) and excluded.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- tox_import(
#'   proj,
#'   system.file("extdata", "import_tox_admetlab.csv", package = "patliR"),
#'   platform = "admetlab"
#' )
#' patliRResults(proj, "tox_imported")
#'
#' @export
tox_import <- function(proj, path, platform = c("admetlab", "swissadme", "other"),
                        column_map = NULL, mapping_file = NULL) {
  stopifnot(is(proj, "PatliRProject"), file.exists(path))
  platform <- match.arg(platform)

  raw <- .read_csv_safe(path)
  cmp <- compounds(proj)

  matched_id <- if (!is.null(mapping_file)) {
    .match_by_export_mapping(nrow(raw), mapping_file)
  } else {
    .tox_import_match_compounds(raw, cmp, platform)$compound_id
  }
  unmatched <- is.na(matched_id)

  if (any(unmatched)) {
    proj <- .log_append(
      proj, step = "tox_import", id = NA_character_,
      message = paste0("tox_import_unmatched: row ", which(unmatched), " of '", basename(path), "' did not match any compound")
    )
  }

  props <- .tox_import_apply_column_map(raw, platform, column_map)
  long <- .to_long_import(matched_id[!unmatched], props[!unmatched, , drop = FALSE], source = platform)

  existing <- patliRResults(proj, "tox_imported")
  patliRResults(proj, "tox_imported") <- if (is.null(existing)) long else rbind(existing, long)

  .write_results_csv(proj, "tox_imported", patliRResults(proj, "tox_imported"))
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.tox_import_match_compounds <- function(raw, cmp, platform) {
  if (platform == "admetlab" || ("smiles" %in% names(raw) && nrow(cmp) > 0)) {
    keys <- .compute_structure_key(as.character(raw[["smiles"]]))
    hit <- match(keys, cmp$canonical_smiles)
    return(data.frame(compound_id = cmp$id[hit], stringsAsFactors = FALSE))
  }
  if (platform == "swissadme" || ("Molecule" %in% names(raw) && nrow(cmp) > 0)) {
    key <- as.character(raw[["Molecule"]])
    hit <- match(key, cmp$pubchem_id)
    return(data.frame(compound_id = cmp$id[hit], stringsAsFactors = FALSE))
  }
  cli::cli_abort(c(
    "Could not determine how to match rows of {.arg path} to {.fn compounds}.",
    "i" = "Expected a {.val smiles} or {.val Molecule} (PubChem CID) column, or run {.fn prep_compounds} first."
  ))
}

#' @keywords internal
.tox_import_apply_column_map <- function(raw, platform, column_map) {
  preset <- switch(platform,
    admetlab = c(Ames = "ames", DILI = "dili", hERG = "herg"),
    swissadme = character(0), # no bundled toxicity example for swissadme; user must pass column_map
    character(0)
  )
  map <- if (is.null(column_map)) preset else column_map
  if (length(map) == 0) {
    cli::cli_abort("No {.arg column_map} available for platform {.val {platform}}; provide one explicitly.")
  }
  available <- intersect(names(map), names(raw))
  out <- raw[, available, drop = FALSE]
  names(out) <- map[available]
  out
}

#' Export compounds as a plain SMILES list, ready to paste into a toxicity
#' platform
#'
#' @description
#' The toxicity-side counterpart of [adme_export_smiles()]. Platforms like
#' ADMETlab or pkCSM take their input as a pasted list of SMILES, one per
#' line, and return their output in that same order with no compound id of
#' their own. This writes that list plus a companion mapping file so the
#' export can be reconciled back to `compound_id` by row position --
#' pass the mapping file to [tox_import()]`(mapping_file = ...)`, which is
#' more reliable than matching on SMILES (a third-party toolkit's canonical
#' SMILES need not match patliR's; see `DESIGN.md`).
#'
#' Round trip: `tox_export_smiles()` -> paste into the platform -> download
#' its export -> `tox_import(..., mapping_file = "<out_file>_map.csv")`.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound currently in [compounds()].
#' @param out_file Character scalar, base path to write to (e.g.
#'   `"chilcuague_tox_smiles"`), or `NULL` (default) to only return the
#'   result. When given, writes `<out_file>.txt` (the plain SMILES list)
#'   and `<out_file>_map.csv` (columns `row_order`, `compound_id`, `name`,
#'   `smiles`).
#'
#' @return Invisibly, `list(smiles_text, mapping)` -- see
#'   [adme_export_smiles()], which shares its implementation.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' export <- tox_export_smiles(proj)
#' cat(export$smiles_text)
#' export$mapping
#'
#' @export
tox_export_smiles <- function(proj, compound_ids = NULL, out_file = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  .export_smiles_list(cmp, out_file)
}
