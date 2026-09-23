#' @include AllGenerics.R internal.R
NULL

#' Import and validate compound identity
#'
#' @description
#' `prep_compounds()` is the main entry point of the `patliR` pipeline for
#' **Scenario B**: a curated table of compounds (one row per compound,
#' already identified by name/CAS/PubChem CID and, ideally, SMILES). For raw
#' GC-MS abundance matrices (**Scenario A**), start with [prep_binarize()]
#' instead.
#'
#' Every row is validated as a real chemical structure (via `rcdk`), given a
#' canonical SMILES (CDK's own canonicalization -- see the note in
#' `.check_structures()`, internal, for why this is not an InChIKey),
#' checked for duplicates, and appended to
#' [compounds()]. Rows that fail any of these checks are **never silently
#' dropped**: they are excluded from `compounds()` but always recorded in
#' [projectLog()] with the reason (`"invalid_structure"`,
#' `"duplicate"`, or `"smiles_not_resolved"`).
#'
#' @param proj A [`PatliRProject-class`] object, typically from
#'   [patliR_project()].
#' @param data A `data.frame` with one row per compound.
#' @param identifier `"pubchem"` (default) if `data` is keyed by PubChem CID
#'   (column `id_col`, default `"PubChemCID"`); `"smiles"` if `data` is keyed
#'   directly by a SMILES column (default `"SMILES"`). In both cases, if a
#'   `SMILES` column is present it is always used directly (no network
#'   lookup needed) -- `identifier` only changes what the *default* `id_col`
#'   is. `on_missing_smiles` (below) applies the same way regardless of
#'   `identifier`: any row with a `PubChemCID` (whether or not that is what
#'   you set `identifier` to) but no SMILES is eligible for the PubChem
#'   fetch when `on_missing_smiles = "fetch"`.
#' @param id_col Character scalar or `NULL`. Name of the identifier column
#'   in `data`. Defaults to `"PubChemCID"` when `identifier = "pubchem"` and
#'   to `"SMILES"` when `identifier = "smiles"`.
#' @param name_col Character scalar or `NULL`. Name of the compound-name
#'   column in `data`. Defaults to `"Name"` if present, else `NA` names are
#'   used.
#' @param dedup Logical. If `TRUE` (default), duplicate compounds (matched
#'   by canonical SMILES, against both `data` itself and any compounds
#'   already in `proj`) are removed, keeping the first occurrence, and
#'   logged as `"duplicate"`.
#' @param on_missing_smiles What to do with rows that have an identifier but
#'   no SMILES. `"abort"` (default): stop with an informative error, since
#'   silently guessing structure is unsafe. `"fetch"`: try to resolve SMILES
#'   from PubChem via the \pkg{webchem} package (Suggests), following the
#'   `abort`/`warn_and_cache` policy of `.fetch_external()` (internal)
#'   with `fetch_mode`. `"drop"`: exclude these rows and log them as
#'   `"smiles_not_resolved"`.
#' @param fetch_mode `"warn_and_cache"` (default) or `"abort"`; only used
#'   when `on_missing_smiles = "fetch"`. See `.fetch_external()` (internal).
#'
#' @return The updated `proj`, with [compounds()] extended and
#'   `01_compounds.csv` (re)written inside [projectDir()].
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' compounds(proj)
#' projectLog(proj)
#'
#' @seealso [prep_compound()] for a single compound, [prep_binarize()] for
#'   Scenario A (raw abundance matrices).
#' @export
prep_compounds <- function(proj, data,
                            identifier = c("pubchem", "smiles"),
                            id_col = NULL,
                            name_col = NULL,
                            dedup = TRUE,
                            on_missing_smiles = c("abort", "fetch", "drop"),
                            fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"), is.data.frame(data), nrow(data) > 0)
  identifier <- match.arg(identifier)
  on_missing_smiles <- match.arg(on_missing_smiles)
  fetch_mode <- match.arg(fetch_mode)

  if (is.null(id_col)) {
    id_col <- if (identifier == "pubchem") "PubChemCID" else "SMILES"
  }
  if (!id_col %in% names(data)) {
    cli::cli_abort("Column {.val {id_col}} (the {.arg id_col}) was not found in {.arg data}.")
  }
  if (is.null(name_col)) {
    name_col <- if ("Name" %in% names(data)) "Name" else NA_character_
  }

  n <- nrow(data)
  pubchem_id <- if ("PubChemCID" %in% names(data)) as.character(data[["PubChemCID"]]) else rep(NA_character_, n)
  name       <- if (!is.na(name_col) && name_col %in% names(data)) as.character(data[[name_col]]) else rep(NA_character_, n)
  smiles     <- if ("SMILES" %in% names(data)) as.character(data[["SMILES"]]) else rep(NA_character_, n)
  source_tag <- rep("input", n)

  ## Resolve missing SMILES. Only rows that also have a PubChemCID can be
  ## fetched (CID -> SMILES via PubChem) -- this used to be gated on
  ## `identifier == "pubchem"`, which meant a Scenario B table keyed by
  ## SMILES (identifier = "smiles") that *also* carries a PubChemCID column
  ## for some rows (e.g. real_data/compound_list_Chilcuague.csv, populated
  ## by a CAS -> CID resolution script run before prep_compounds()) never
  ## got its missing SMILES fetched at all -- rows were silently dropped
  ## ("smiles_not_resolved") even though the CID needed to fetch them was
  ## right there, which then cascaded into targets_import_batch() reporting
  ## every one of those compounds' Targets<cid>.csv as unmatched. What
  ## `identifier` actually changes is only the default `id_col` (see above)
  ## and, implicitly, whether `pubchem_id` was even filled from `data` in
  ## the first place -- SMILES resolution itself should not care which one
  ## the caller nominally chose.
  missing_smiles <- is.na(smiles) | !nzchar(smiles)
  has_pubchem_id <- !is.na(pubchem_id) & nzchar(pubchem_id)
  fetchable <- missing_smiles & has_pubchem_id
  if (any(missing_smiles)) {
    if (on_missing_smiles == "abort" && any(missing_smiles & !has_pubchem_id)) {
      cli::cli_abort(c(
        "{sum(missing_smiles & !has_pubchem_id)} row(s) have no SMILES and no PubChemCID to fetch one from, and {.arg on_missing_smiles} is {.val abort}.",
        "i" = "Add a SMILES column to {.arg data}, or set {.arg on_missing_smiles} to {.val fetch} or {.val drop}."
      ))
    }
    if (on_missing_smiles == "abort" && any(fetchable)) {
      cli::cli_abort(c(
        "{sum(fetchable)} row(s) have a PubChemCID but no SMILES, and {.arg on_missing_smiles} is {.val abort}.",
        "i" = "Set {.arg on_missing_smiles} to {.val fetch} to resolve them from PubChem, or {.val drop} to exclude them."
      ))
    }
    if (on_missing_smiles == "fetch" && any(fetchable)) {
      fetched <- .fetch_smiles_from_pubchem(pubchem_id[fetchable], cacheDir(proj), fetch_mode)
      smiles[fetchable] <- fetched
      ## `fetched` is only as long as `sum(fetchable)`; build a full-length
      ## logical mask before indexing `source_tag` so this can never
      ## silently misalign via recycling.
      fetched_ok <- fetchable
      fetched_ok[fetchable] <- !is.na(fetched)
      source_tag[fetched_ok] <- "pubchem_fetch"
    }
    ## Rows with no PubChemCID either (nothing to fetch from) and rows left
    ## over after on_missing_smiles == "drop": leave as NA, handled by the
    ## validity check below (logged as "smiles_not_resolved").
  }

  has_smiles <- !is.na(smiles) & nzchar(smiles)
  valid <- rep(FALSE, n)
  canonical_smiles <- rep(NA_character_, n)
  reason <- rep(NA_character_, n)
  reason[!has_smiles] <- "smiles_not_resolved"

  if (any(has_smiles)) {
    checked <- .check_structures(smiles[has_smiles])
    valid[has_smiles] <- checked$valid
    canonical_smiles[has_smiles] <- checked$canonical_smiles
    ## Real rcdk/rJava error message when available, not a generic label --
    ## if rcdk/Java itself is misconfigured, this is what tells you (rather
    ## than every compound silently vanishing with no clue why).
    reason[has_smiles] <- ifelse(checked$valid, NA_character_,
                                  paste0("invalid_structure: ", checked$error))
  }

  new_rows <- data.frame(
    pubchem_id = pubchem_id, name = name, smiles = smiles,
    canonical_smiles = canonical_smiles, source = source_tag,
    stringsAsFactors = FALSE
  )

  ## Log and drop rows that never resolved to a valid structure.
  bad <- !valid
  if (any(bad)) {
    proj <- .log_append(
      proj, step = "prep_compounds",
      id = ifelse(is.na(new_rows$pubchem_id[bad]), new_rows$name[bad], new_rows$pubchem_id[bad]),
      message = paste0(reason[bad], " (row ", which(bad), " of input)")
    )
  }
  new_rows <- new_rows[valid, , drop = FALSE]

  existing <- compounds(proj)

  if (dedup) {
    dup_within <- duplicated(new_rows$canonical_smiles)
    dup_existing <- new_rows$canonical_smiles %in% existing$canonical_smiles
    dup <- dup_within | dup_existing
    if (any(dup)) {
      proj <- .log_append(
        proj, step = "prep_compounds",
        id = ifelse(is.na(new_rows$pubchem_id[dup]), new_rows$name[dup], new_rows$pubchem_id[dup]),
        message = "duplicate (matched by canonical SMILES)"
      )
    }
    new_rows <- new_rows[!dup, , drop = FALSE]
  }

  if (nrow(new_rows) > 0) {
    new_rows$id <- .next_compound_ids(existing$id, nrow(new_rows))
    new_rows <- new_rows[, c("id", "pubchem_id", "name", "smiles", "canonical_smiles", "source")]
    compounds(proj) <- rbind(existing, new_rows)
  }

  .write_step_csv(proj, "01_compounds.csv", compounds(proj))
  .write_log_csv(proj)
  proj
}

#' Fetch missing SMILES from PubChem by CID, preferring the modern
#' stereo-bearing response key
#'
#' @description
#' PubChem renamed its PUG-REST SMILES properties (2025): the old
#' `CanonicalSMILES` now comes back as `ConnectivitySMILES` and the
#' stereo-bearing `IsomericSMILES` as the unqualified `SMILES`. The old
#' names are still accepted on the *request* (HTTP 200) but the response
#' key uses the new name, so reading `$CanonicalSMILES` silently returns
#' `NULL`. We request the new names and read whichever column is present,
#' preferring the stereo-bearing one, and fall back to the legacy keys so a
#' future revert cannot silently break this again.
#' @keywords internal
.fetch_smiles_from_pubchem <- function(pubchem_ids, cache_dir, fetch_mode) {
  pick_smiles <- function(res) {
    for (col in c("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES")) {
      if (col %in% names(res)) {
        v <- as.character(res[[col]][[1]])
        if (!is.na(v) && nzchar(v)) return(v)
      }
    }
    NA_character_
  }
  vapply(pubchem_ids, function(cid) {
    if (is.na(cid) || !nzchar(cid)) return(NA_character_)
    fetched <- .fetch_external(
      fetch_fun = function() {
        if (!requireNamespace("webchem", quietly = TRUE)) {
          stop("the 'webchem' package is required to fetch SMILES from PubChem; install it or provide SMILES directly.")
        }
        res <- webchem::pc_prop(as.character(cid),
                                 properties = c("SMILES", "ConnectivitySMILES"))
        s <- pick_smiles(res)
        if (is.na(s)) stop("PubChem returned no SMILES for CID ", cid)
        s
      },
      cache_dir = cache_dir,
      ## key bumped (was pubchem_smiles_) so caches poisoned by the
      ## pre-rename NULL bug are not read back.
      cache_key = paste0("pubchem_cid_smiles_", cid),
      mode = fetch_mode
    )
    if (is.null(fetched)) NA_character_ else as.character(fetched)
  }, character(1), USE.NAMES = FALSE)
}

#' Import and validate a single compound
#'
#' Thin wrapper around [prep_compounds()] for the common case of adding one
#' compound at a time (e.g. interactively, or one row read from a form).
#'
#' @inheritParams prep_compounds
#' @param data_row A single-row `data.frame` (or a named `list`/vector
#'   coercible to one), with the same columns [prep_compounds()] expects.
#'
#' @return The updated `proj`. See [prep_compounds()].
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' one <- read.csv(
#'   system.file("extdata", "input_single_compound.csv", package = "patliR")
#' )
#' proj <- prep_compound(proj, one, identifier = "pubchem")
#' compounds(proj)
#'
#' @export
prep_compound <- function(proj, data_row,
                           identifier = c("pubchem", "smiles"),
                           id_col = NULL,
                           name_col = NULL,
                           dedup = TRUE,
                           on_missing_smiles = c("abort", "fetch", "drop"),
                           fetch_mode = c("warn_and_cache", "abort")) {
  if (!is.data.frame(data_row)) data_row <- as.data.frame(data_row, stringsAsFactors = FALSE)
  if (nrow(data_row) != 1) {
    cli::cli_abort("{.arg data_row} must have exactly one row (got {nrow(data_row)}); use {.fn prep_compounds} for more than one.")
  }
  prep_compounds(
    proj, data_row,
    identifier = identifier, id_col = id_col, name_col = name_col,
    dedup = dedup, on_missing_smiles = on_missing_smiles, fetch_mode = fetch_mode
  )
}

#' @keywords internal
.write_log_csv <- function(proj) {
  path <- file.path(projectDir(proj), "patliR_log.csv")
  .atomic_write_csv(projectLog(proj), path)
  invisible(path)
}
