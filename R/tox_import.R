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
#'
#' @details
#' Built-in presets (confirmed against the column headers the bundled
#' example file actually uses):
#' \itemize{
#'   \item `platform = "admetlab"`: `smiles` (canonicalized the same way
#'     [prep_compounds()] does, matched against
#'     `compounds(proj)$canonical_smiles`), `Ames`, `DILI`, `hERG` (all
#'     mapped to lowercase property names).
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
                        column_map = NULL) {
  stopifnot(is(proj, "PatliRProject"), file.exists(path))
  platform <- match.arg(platform)

  raw <- .read_csv_safe(path)
  cmp <- compounds(proj)

  match_result <- .tox_import_match_compounds(raw, cmp, platform)
  matched_id <- match_result$compound_id
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
