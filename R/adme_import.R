#' @include AllGenerics.R internal.R
NULL

#' Import ADME properties exported from an external platform
#'
#' @description
#' Imports a CSV exported from a third-party ADME platform (SwissADME,
#' ADMETlab, pkCSM, ...). `patliR` never scrapes these platforms -- you run
#' them yourself in your browser, download the CSV they export, and this
#' function reconciles it against [compounds()].
#'
#' Different platforms key their export by different identifiers: the
#' bundled example for SwissADME uses `Molecule` (a PubChem CID), while the
#' ADMETlab example uses `smiles` directly. `adme_import()` reconciles
#' either against [compounds()] via a canonical-SMILES structure key
#' (recomputing it from `smiles` when needed -- see the note in
#' `.check_structures()`, internal, for why this is not an InChIKey), so it
#' never assumes a shared identifier across platforms.
#'
#' @inheritParams compounds
#' @param path Path to the CSV file to import.
#' @param platform `"swissadme"`, `"admetlab"`, or `"other"`. Selects a
#'   built-in `column_map` (see Details) so you do not have to write one by
#'   hand for a supported platform's default export headers.
#' @param column_map Named character vector mapping *your* file's column
#'   names to patliR's property names, e.g.
#'   `c(MW = "mw", Consensus.Log.Po.w = "logp")`. Required when
#'   `platform = "other"`; ignored (a preset is used instead) for
#'   `"swissadme"`/`"admetlab"` unless you explicitly pass one, in which
#'   case yours takes precedence.
#'
#' @details
#' Built-in presets (confirmed against the column headers real exports
#' actually use):
#' \itemize{
#'   \item `platform = "swissadme"`: `Molecule` (PubChem CID, used to match
#'     `compounds(proj)$pubchem_id`), `MW`, `Consensus.Log.Po.w` (-> `logp`),
#'     `HBD`, `HBA`, `TPSA`, `GI.absorption`, `BBB.permeant`,
#'     `Pgp.substrate`, `Bioavailability.Score`.
#'   \item `platform = "admetlab"`: `smiles` (canonicalized the same way
#'     [prep_compounds()] does, and matched against
#'     `compounds(proj)$canonical_smiles` -- no shared numeric id with
#'     SwissADME). Note: the bundled `import_tox_admetlab.csv` example is
#'     toxicity data (see `tox_import()`, not yet implemented), not ADME --
#'     there is no ready-to-run ADMETlab *ADME* example yet.
#' }
#'
#' @return The updated `proj`, with an `adme_imported` entry in
#'   [patliRResults()] in long format (columns `compound_id`, `property`,
#'   `value`, `source`, `import_date`), also written to
#'   `results/adme_imported.csv`. Rows that could not be matched to a known
#'   compound are logged (`"adme_import_unmatched"`) and excluded.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_import(
#'   proj,
#'   system.file("extdata", "import_adme_swissadme.csv", package = "patliR"),
#'   platform = "swissadme"
#' )
#' patliRResults(proj, "adme_imported")
#'
#' @export
adme_import <- function(proj, path, platform = c("swissadme", "admetlab", "other"),
                         column_map = NULL) {
  stopifnot(is(proj, "PatliRProject"), file.exists(path))
  platform <- match.arg(platform)

  raw <- .read_csv_safe(path)
  cmp <- compounds(proj)

  match_result <- .import_match_compounds(raw, cmp, platform)
  matched_id <- match_result$compound_id
  unmatched <- is.na(matched_id)

  if (any(unmatched)) {
    proj <- .log_append(
      proj, step = "adme_import", id = NA_character_,
      message = paste0("adme_import_unmatched: row ", which(unmatched), " of '", basename(path), "' did not match any compound")
    )
  }

  props <- .import_apply_column_map(raw, platform, column_map)
  long <- .to_long_import(matched_id[!unmatched], props[!unmatched, , drop = FALSE], source = platform)

  existing <- patliRResults(proj, "adme_imported")
  patliRResults(proj, "adme_imported") <- if (is.null(existing)) long else rbind(existing, long)

  .write_results_csv(proj, "adme_imported", patliRResults(proj, "adme_imported"))
  .write_log_csv(proj)
  proj
}

#' Export compounds as a plain SMILES list, ready to paste into SwissADME
#' (or similar)
#'
#' @description
#' `adme_local()` covers Ro5/Veber/Ghose/Egan/Oprea and BOILED-Egg from
#' `rcdk` alone, no network needed -- but a platform like SwissADME
#' computes properties `patliR` does not (e.g. its own consensus LogP,
#' synthetic accessibility, P-gp substrate/CYP inhibition predictions) that
#' can only be brought in via [adme_import()]. Those platforms take their
#' input as a pasted list of SMILES, one per line, in their own web form --
#' `patliR` does not scrape them (same policy as [adme_import()]'s own
#' documentation) -- so this writes exactly that: one canonical SMILES per
#' line, in a fixed, deterministic order, plus a companion mapping file so
#' the platform's output (which is only ever ordered by input position, no
#' compound id of its own) can be matched back to `compound_id`/`name`
#' afterwards. Round-trip: `adme_export_smiles()` -> paste into the
#' platform -> download its export -> [adme_import()].
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound currently in [compounds()].
#' @param out_file Character scalar, base path to write to (e.g.
#'   `"chilcuague_smiles"`), or `NULL` (default) to only return the result
#'   without writing anything. When given, writes `<out_file>.txt` (the
#'   plain SMILES list, ready to paste) and `<out_file>_map.csv` (columns
#'   `row_order`, `compound_id`, `name`, `smiles` -- `row_order` matches
#'   the line number in the `.txt` file, i.e. the row order the platform's
#'   own output will come back in).
#'
#' @return Invisibly, a `list(smiles_text, mapping)`: `smiles_text` the
#'   full newline-joined SMILES list as a single string, `mapping` the
#'   `data.frame` described above (also written to `<out_file>_map.csv`
#'   when `out_file` is given).
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' export <- adme_export_smiles(proj)
#' cat(export$smiles_text)
#' export$mapping
#'
#' @export
adme_export_smiles <- function(proj, compound_ids = NULL, out_file = NULL) {
  stopifnot(is(proj, "PatliRProject"))

  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  if (nrow(cmp) == 0) {
    cli::cli_abort("No matching compounds in {.arg proj}; run {.fn prep_compounds} first.")
  }

  smiles <- ifelse(!is.na(cmp$canonical_smiles) & nzchar(cmp$canonical_smiles), cmp$canonical_smiles, cmp$smiles)
  missing <- is.na(smiles) | !nzchar(smiles)
  if (any(missing)) {
    cli::cli_warn("{sum(missing)} compound(s) have no SMILES at all and are excluded from the export ({.val {cmp$id[missing]}}).")
  }

  mapping <- data.frame(
    row_order = seq_len(sum(!missing)),
    compound_id = cmp$id[!missing],
    name = cmp$name[!missing],
    smiles = smiles[!missing],
    stringsAsFactors = FALSE
  )
  smiles_text <- paste(mapping$smiles, collapse = "\n")

  if (!is.null(out_file)) {
    writeLines(mapping$smiles, paste0(out_file, ".txt"))
    utils::write.csv(mapping, paste0(out_file, "_map.csv"), row.names = FALSE)
  }

  invisible(list(smiles_text = smiles_text, mapping = mapping))
}

#' @keywords internal
.import_match_compounds <- function(raw, cmp, platform) {
  if (platform == "swissadme" || ("Molecule" %in% names(raw) && nrow(cmp) > 0)) {
    key <- as.character(raw[["Molecule"]])
    hit <- match(key, cmp$pubchem_id)
    return(data.frame(compound_id = cmp$id[hit], stringsAsFactors = FALSE))
  }
  if (platform == "admetlab" || ("smiles" %in% names(raw) && nrow(cmp) > 0)) {
    ## Canonicalize with the same CDK flavor prep_compounds() used, so the
    ## join is against patliR's own canonical_smiles key (see the note in
    ## .check_structures(), internal, for why this isn't an InChIKey).
    keys <- .compute_structure_key(as.character(raw[["smiles"]]))
    hit <- match(keys, cmp$canonical_smiles)
    return(data.frame(compound_id = cmp$id[hit], stringsAsFactors = FALSE))
  }
  cli::cli_abort(c(
    "Could not determine how to match rows of {.arg path} to {.fn compounds}.",
    "i" = "Expected a {.val Molecule} (PubChem CID) or {.val smiles} column, or run {.fn prep_compounds} first."
  ))
}

#' @keywords internal
.import_apply_column_map <- function(raw, platform, column_map) {
  preset <- switch(platform,
    swissadme = c(
      MW = "mw", `Consensus.Log.Po.w` = "logp", HBD = "hbd", HBA = "hba",
      TPSA = "tpsa", `GI.absorption` = "gi_absorption", `BBB.permeant` = "bbb_permeant",
      `Pgp.substrate` = "pgp_substrate", `Bioavailability.Score` = "bioavailability_score"
    ),
    admetlab = character(0), # admetlab example only carries toxicity columns; see tox_import()
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

#' @keywords internal
.to_long_import <- function(compound_id, props, source) {
  if (nrow(props) == 0 || ncol(props) == 0) {
    return(data.frame(compound_id = character(), property = character(), value = character(),
                       source = character(), import_date = as.Date(character()), stringsAsFactors = FALSE))
  }
  do.call(rbind, lapply(names(props), function(p) {
    data.frame(
      compound_id = compound_id, property = p, value = as.character(props[[p]]),
      source = source, import_date = Sys.Date(), stringsAsFactors = FALSE
    )
  }))
}
