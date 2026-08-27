#' @include AllGenerics.R internal.R
NULL

#' Import predicted targets exported from an external platform
#'
#' @description
#' Imports a CSV exported from a target-prediction platform
#' (SwissTargetPrediction, SuperPred, ...). `patliR` does not run its own
#' target-prediction model (see `ROADMAP.md`) -- you run the platform
#' yourself, download its export, and this function reconciles it against
#' [compounds()].
#'
#' Most of these platforms only accept one compound per run, so their
#' natural export convention is one file per compound. The default
#' (`id_from = "filename"`) expects files named exactly `Targets<pubchem_id>.csv`
#' (no separator, e.g. `Targets5280443.csv`).
#'
#' @inheritParams compounds
#' @param path Path to the CSV file to import.
#' @param platform `"swisstargetprediction"`, `"superpred"`, or `"other"`.
#'   `"superpred"` has a built-in default for `target_col`/`probability_col`/
#'   `confidence_col` (confirmed against a real SuperPred export: `"UniProt
#'   ID"`, `"Probability"`, `"Model accuracy"`) which is used whenever you
#'   don't pass those arguments explicitly. `"swisstargetprediction"` and
#'   `"other"` have no default -- `target_col` and `probability_col` are
#'   always required for them.
#' @param target_col Name of the column holding the UniProt ID in `path`.
#'   Always required unless `platform = "superpred"` supplies the default.
#'   Matched case/whitespace/decoration-insensitively against `path`'s
#'   actual columns (e.g. `"*Probability"` still matches `"Probability"`),
#'   so a cosmetic export difference between runs of the same platform
#'   does not abort the import.
#' @param probability_col Name of the column holding the prediction
#'   probability/score. Always required unless `platform = "superpred"`
#'   supplies the default. Parsed tolerant of either convention real
#'   platforms use -- a plain `0..1` fraction, or a `"96.55%"` percentage
#'   string -- and always stored as a `0..1` fraction in the result (see
#'   `.parse_percent_column()`, internal).
#' @param confidence_col Name of an optional model-confidence/accuracy
#'   column, or `NULL` if the platform does not provide one. Same
#'   percentage/fraction handling as `probability_col`.
#' @param id_from `"filename"` (default): extract the compound's PubChem CID
#'   from the file name (`Targets<pubchem_id>.csv`) and match against
#'   `compounds(proj)$pubchem_id`. `"column"`: match using a column inside
#'   the file itself, named by `compound_col`.
#' @param compound_col Only used when `id_from = "column"`: name of the
#'   column in `path` holding the PubChem CID to match against
#'   `compounds(proj)$pubchem_id`.
#'
#' @return The updated `proj`, with a `targets_imported` entry in
#'   [patliRResults()] (columns `compound_id`, `uniprot_id`, `probability`,
#'   `confidence`, `source`, `import_date`), also written to
#'   `results/targets_imported.csv`. Rows that could not be matched to a
#'   known compound are logged (`"targets_import_unmatched"`) and excluded.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- targets_import(
#'   proj,
#'   system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
#'   platform = "superpred"
#' )
#' patliRResults(proj, "targets_imported")
#' }
#'
#' @export
targets_import <- function(proj, path, platform = c("swisstargetprediction", "superpred", "other"),
                            target_col = NULL, probability_col = NULL, confidence_col = NULL,
                            id_from = c("filename", "column"), compound_col = NULL) {
  stopifnot(is(proj, "PatliRProject"), file.exists(path))
  platform <- match.arg(platform)
  id_from <- match.arg(id_from)

  if (platform == "superpred") {
    if (is.null(target_col)) target_col <- "UniProt ID"
    if (is.null(probability_col)) probability_col <- "Probability"
    if (is.null(confidence_col)) confidence_col <- "Model accuracy"
  }
  if (is.null(target_col) || is.null(probability_col)) {
    cli::cli_abort(c(
      "{.arg target_col} and {.arg probability_col} are required.",
      "i" = "Only {.val superpred} has a built-in default; for {.val swisstargetprediction}/{.val other} pass them explicitly."
    ))
  }

  raw <- .read_csv_safe(path)

  ## Real platform exports are not always consistent about a column's exact
  ## spelling (seen on real SuperPred data: "Probability" vs.
  ## "*Probability") -- match tolerant of case/whitespace/decoration before
  ## giving up, so a cosmetic export difference doesn't abort the whole
  ## import (see .match_column_flexible(), internal).
  target_col_resolved <- .match_column_flexible(names(raw), target_col)
  probability_col_resolved <- .match_column_flexible(names(raw), probability_col)
  if (is.na(target_col_resolved)) {
    cli::cli_abort("Column {.val {target_col}} not found in {.arg path} (even allowing for case/whitespace/decoration differences). Available columns: {.val {names(raw)}}.")
  }
  if (is.na(probability_col_resolved)) {
    cli::cli_abort("Column {.val {probability_col}} not found in {.arg path} (even allowing for case/whitespace/decoration differences). Available columns: {.val {names(raw)}}.")
  }
  confidence_col_resolved <- if (!is.null(confidence_col)) .match_column_flexible(names(raw), confidence_col) else NA_character_

  cmp <- compounds(proj)
  ## .targets_match_compound() already returns one compound_id per row of
  ## `raw` (rep()'d internally for id_from = "filename", matched row-by-row
  ## for id_from = "column") -- it must NOT be rep()'d again here. An
  ## earlier version of this line did `rep(compound_id, n)` on top of that
  ## already-full-length vector; since data.frame() silently recycles a
  ## column whose length evenly divides the frame's row count, that bug
  ## never errored -- it just produced n^2 rows (every real row repeated n
  ## times) for every single import, undetected because no existing test
  ## asserted an exact row count.
  compound_id <- .targets_match_compound(path, raw, cmp, id_from, compound_col)

  n <- nrow(raw)
  long <- data.frame(
    compound_id = compound_id,
    uniprot_id = as.character(raw[[target_col_resolved]]),
    ## Real exports mix conventions for this column: a plain 0..1 fraction,
    ## or a "96.55%" string (which as.numeric() would silently turn into
    ## NA because of the '%') -- .parse_percent_column() (internal) handles
    ## both and always returns a 0..1 fraction.
    probability = .parse_percent_column(raw[[probability_col_resolved]]),
    confidence = if (!is.na(confidence_col_resolved)) {
      .parse_percent_column(raw[[confidence_col_resolved]])
    } else {
      NA_real_
    },
    source = platform,
    import_date = Sys.Date(),
    stringsAsFactors = FALSE
  )

  unmatched <- is.na(long$compound_id)
  if (any(unmatched)) {
    proj <- .log_append(
      proj, step = "targets_import", id = NA_character_,
      message = paste0("targets_import_unmatched: '", basename(path), "' did not match any compound (id_from = '", id_from, "')")
    )
  }
  long <- long[!unmatched, , drop = FALSE]

  existing <- patliRResults(proj, "targets_imported")
  patliRResults(proj, "targets_imported") <- if (is.null(existing)) long else rbind(existing, long)

  .write_results_csv(proj, "targets_imported", patliRResults(proj, "targets_imported"))
  .write_log_csv(proj)
  proj
}

#' Import a whole folder of per-compound target-prediction exports
#'
#' @description
#' Convenience wrapper around [targets_import()] for the common case: a
#' folder full of `Targets<pubchem_id>.csv` files (one per compound,
#' downloaded one at a time from a platform like SuperPred). Reads every
#' file matching that naming pattern, calls [targets_import()] on each, and
#' returns everything unified into `proj`'s single `targets_imported` table.
#'
#' @section One bad file never loses the rest of the batch:
#' Real exports are not always uniform -- a platform can change its export
#' columns between compounds/sessions (e.g. a file with only `"PDB
#' Visualization"`/`"TTD ID"` and no `"Probability"` column at all, seen
#' against real Chilcuague data). [targets_import()] itself still aborts
#' loudly on a single direct call (a bad file is a real problem to
#' surface immediately when you only import one), but here each file goes
#' through its own `tryCatch()` -- same per-item pattern already used by
#' [network_pathview()] for per-pathway rendering failures -- so one file
#' failing (missing/renamed columns, unparseable content, ...) is logged
#' and excluded, never crashes the whole folder and loses every file
#' already imported before it.
#'
#' @inheritParams targets_import
#' @param dir Directory containing `Targets<pubchem_id>.csv` files.
#'
#' @return The updated `proj` (same `targets_imported` slot as
#'   [targets_import()], accumulated across every file in `dir` that
#'   imported successfully), plus a `targets_import_batch_log` entry in
#'   [patliRResults()] (also written to
#'   `results/targets_import_batch_log.csv`) with **one row per file found
#'   in `dir`** -- columns `path`, `pubchem_id`, `compound_id`,
#'   `compound_name` (the last two resolved from [compounds()] whenever
#'   `pubchem_id` matches, so an excluded file is identifiable by name, not
#'   just by CID), `ok`, `reason` (`"imported"`, `"skipped_naming"` -- file
#'   name does not match `^Targets\\d+\\.csv$`, or `"import_failed"` --
#'   matched the naming pattern but [targets_import()] errored on it), and
#'   `message` (`NA` when `ok`, the real error text otherwise). Every file
#'   is logged (`"targets_import_batch_skipped"`/
#'   `"targets_import_batch_excluded"`), not silently ignored.
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
#' patliRResults(proj, "targets_imported")
#' patliRResults(proj, "targets_import_batch_log")
#' }
#'
#' @export
targets_import_batch <- function(proj, dir, platform = c("swisstargetprediction", "superpred", "other"),
                                  target_col = NULL, probability_col = NULL, confidence_col = NULL) {
  stopifnot(is(proj, "PatliRProject"), dir.exists(dir))
  platform <- match.arg(platform)

  all_files <- list.files(dir, full.names = TRUE)
  file_names <- basename(all_files)
  is_valid <- grepl("^Targets\\d+\\.csv$", file_names)

  log_rows <- list()

  skipped <- file_names[!is_valid]
  if (length(skipped) > 0) {
    proj <- .log_append(
      proj, step = "targets_import_batch", id = NA_character_,
      message = paste0("targets_import_batch_skipped: '", skipped, "' does not match the 'Targets<pubchem_id>.csv' naming pattern")
    )
    log_rows[["skipped"]] <- data.frame(
      path = skipped, pubchem_id = NA_character_, compound_id = NA_character_,
      compound_name = NA_character_, ok = FALSE, reason = "skipped_naming",
      message = "file name does not match 'Targets<pubchem_id>.csv'", stringsAsFactors = FALSE
    )
  }

  valid_files <- all_files[is_valid]
  imported_rows <- vector("list", length(valid_files))
  for (i in seq_along(valid_files)) {
    f <- valid_files[i]
    cmp <- compounds(proj) # re-read every iteration: earlier files in this same batch may have added compounds
    pubchem_id <- .targets_pubchem_id_from_filename(f)
    hit <- match(pubchem_id, cmp$pubchem_id)
    compound_id <- if (!is.na(hit)) cmp$id[hit] else NA_character_
    compound_name <- if (!is.na(hit)) cmp$name[hit] else NA_character_

    res <- tryCatch(
      targets_import(
        proj, f, platform = platform,
        target_col = target_col, probability_col = probability_col, confidence_col = confidence_col,
        id_from = "filename"
      ),
      error = function(e) e
    )
    ok <- !inherits(res, "error")
    if (ok) {
      proj <- res
    } else {
      proj <- .log_append(
        proj, step = "targets_import_batch", id = pubchem_id,
        message = paste0(
          "targets_import_batch_excluded: '", basename(f), "' (pubchem_id=", pubchem_id,
          if (!is.na(compound_name)) paste0(", name='", compound_name, "'") else " -- pubchem_id not found in compounds(proj)",
          ") could not be imported -- ", conditionMessage(res)
        )
      )
    }

    imported_rows[[i]] <- data.frame(
      path = basename(f), pubchem_id = pubchem_id, compound_id = compound_id,
      compound_name = compound_name, ok = ok, reason = if (ok) "imported" else "import_failed",
      message = if (ok) NA_character_ else conditionMessage(res),
      stringsAsFactors = FALSE
    )
  }
  log_rows[["imported"]] <- do.call(rbind, imported_rows)

  batch_log <- do.call(rbind, log_rows)
  if (!is.null(batch_log) && nrow(batch_log) > 0) {
    rownames(batch_log) <- NULL
    patliRResults(proj, "targets_import_batch_log") <- batch_log
    .write_results_csv(proj, "targets_import_batch_log", batch_log)
  }
  .write_log_csv(proj)

  proj
}

#' Extract the PubChem CID out of a `Targets<pubchem_id>.csv` file name
#' @return Character scalar, or `NA_character_` if the file name has no
#'   digits at all.
#' @keywords internal
.targets_pubchem_id_from_filename <- function(path) {
  m <- regmatches(basename(path), regexpr("\\d+", basename(path)))
  if (length(m) == 0 || !nzchar(m)) NA_character_ else m
}

#' @keywords internal
.targets_match_compound <- function(path, raw, cmp, id_from, compound_col) {
  if (id_from == "filename") {
    m <- .targets_pubchem_id_from_filename(path)
    if (is.na(m)) {
      cli::cli_abort(c(
        "Could not extract a PubChem CID from the file name {.file {basename(path)}}.",
        "i" = "Expected the pattern {.val Targets<pubchem_id>.csv}, e.g. {.val Targets5280443.csv}."
      ))
    }
    hit <- match(m, cmp$pubchem_id)
    return(rep(cmp$id[hit], nrow(raw)))
  }
  ## id_from == "column"
  if (is.null(compound_col)) {
    cli::cli_abort("{.arg compound_col} is required when {.code id_from = \"column\"}.")
  }
  if (!compound_col %in% names(raw)) {
    cli::cli_abort("Column {.val {compound_col}} not found in {.arg path}.")
  }
  hit <- match(as.character(raw[[compound_col]]), cmp$pubchem_id)
  cmp$id[hit]
}
