## Internal helpers, not exported. Naming convention: prefixed with '.'.
## Kept in one file since they are short, single-purpose, and reused across
## several exported functions (prep_*, refdb_*, and future families).

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

#' Point-in-polygon test (standard ray-casting / even-odd rule)
#'
#' Vectorized over `(x, y)`; `poly_x`/`poly_y` describe a single closed (or
#' unclosed -- both work) polygon. Used by [plot_boiled_egg()] and
#' `adme_local()`'s BOILED-Egg classification to test whether a compound's
#' (TPSA, WLogP-proxy) point falls inside the published GIA/BBB ellipses.
#'
#' @return Logical vector, same length as `x`/`y`.
#' @keywords internal
.point_in_polygon <- function(x, y, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    intersects <- (poly_y[i] > y) != (poly_y[j] > y)
    slope <- (poly_x[j] - poly_x[i]) * (y - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i]
    cross <- intersects & (x < slope)
    inside <- xor(inside, ifelse(is.na(cross), FALSE, cross))
    j <- i
  }
  inside
}

#' Load the digitized BOILED-Egg GIA (white) and BBB (yolk) ellipse
#' boundaries bundled in `inst/extdata`
#'
#' @description
#' These coordinates are the published model's own numeric parameters
#' from the supporting information of Daina, A. & Zoete, V. (2016), "A
#' BOILED-Egg To Predict Gastrointestinal Absorption and Brain
#' Penetration of Small Molecules", *ChemMedChem* 11, 1117-1121 -- copied
#' directly from the hard-coded coordinate lists in
#' <https://github.com/bfmilne/PyBOILEDegg>'s source (Milne, B.F., 2021,
#' Zenodo \doi{10.5281/zenodo.4725530}), whose own comment states the
#' same original source. Not produced by running that program (it only
#' ever outputs a GI/BBB classification, never boundary coordinates) --
#' this is a direct transcription of published numeric data, not a
#' derivative of any code Milne wrote. Used as-is (TPSA on the x-axis,
#' WLogP on the y-axis), not re-derived by us.
#' @return A `list(gia = data.frame(tpsa, wlogp), bbb = data.frame(tpsa, wlogp))`.
#' @keywords internal
.load_boiled_egg_polygons <- function() {
  list(
    gia = utils::read.csv(system.file("extdata", "boiled_egg_gia.csv", package = "patliR"), stringsAsFactors = FALSE),
    bbb = utils::read.csv(system.file("extdata", "boiled_egg_bbb.csv", package = "patliR"), stringsAsFactors = FALSE)
  )
}

#' Append rows to a project's event log and return the updated project
#'
#' @param proj A [`PatliRProject-class`] object.
#' @param step Character scalar naming the calling function, e.g.
#'   `"prep_compounds"`.
#' @param id Character vector of affected compound/entity ids, recycled
#'   against `message` if needed. Use `NA_character_` for project-level
#'   events that are not about one specific compound.
#' @param message Character vector of human-readable event descriptions.
#' @return `proj`, with the new rows appended to its `log` slot.
#' @keywords internal
.log_append <- function(proj, step, id = NA_character_, message) {
  n <- max(length(id), length(message))
  new_rows <- data.frame(
    step = rep_len(step, n),
    id = rep_len(id, n),
    message = rep_len(message, n),
    timestamp = rep(Sys.time(), n),
    stringsAsFactors = FALSE
  )
  projectLog(proj) <- rbind(projectLog(proj), new_rows)
  proj
}

#' Write a CSV without ever leaving a partially-written file at `path`
#'
#' @description
#' `write.csv(data, path, ...)` writes directly to the final path; a crash,
#' power loss, or kill signal partway through leaves `path` truncated --
#' for the three callers below, that file *is* DESIGN.md's "durable,
#' portable source of truth", so a partial write there is silent data
#' loss on the next [patliR_load()], not just a rerun-and-move-on
#' inconvenience. Write to a temp file in the *same* directory (so the
#' rename below stays on one filesystem/volume, a precondition for it
#' being atomic at all) and rename it over the destination -- readers only
#' ever see the old complete file or the new complete file, never a
#' half-written one. `file.rename()` replacing an existing destination is
#' reliable on the platforms `patliR` targets (confirmed on Windows, where
#' older guidance sometimes assumed otherwise); the copy+remove fallback
#' below exists only for the unlikely case it is not, on some filesystem
#' this was never tested against.
#' @keywords internal
.atomic_write_csv <- function(data, path, row.names = FALSE) {
  dir <- dirname(path)
  tmp <- tempfile(pattern = ".patliR_tmp_", tmpdir = dir, fileext = ".csv")
  utils::write.csv(data, tmp, row.names = row.names)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
  }
  invisible(path)
}

#' Write one of the numbered "source of truth" CSV files for a project
#' @keywords internal
.write_step_csv <- function(proj, filename, data) {
  path <- file.path(projectDir(proj), filename)
  .atomic_write_csv(data, path)
  invisible(path)
}

#' Write one entry of the results bag to `results/<name>.csv`
#' @keywords internal
.write_results_csv <- function(proj, name, data) {
  results_dir <- file.path(projectDir(proj), "results")
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  .atomic_write_csv(data, file.path(results_dir, paste0(name, ".csv")))
  invisible(NULL)
}

#' Escape literal curly braces so arbitrary text is safe to embed as a cli
#' bullet
#'
#' @description
#' `cli::cli_warn()`/`cli::cli_abort()` treat every bullet string as a
#' glue template -- any `{...}` inside it is parsed as an R expression to
#' evaluate, not shown literally. `.fetch_external()` (below) embeds
#' `conditionMessage()` from whatever error the external call raised --
#' text `patliR` does not control, e.g. a raw JSON fragment echoed back by
#' a failed API parse. A message containing `{`/`}` (real case: NPClassifier
#' returning `"class_results": [...`) makes cli try to `parse()` that
#' fragment as R code and throw a confusing meta-error that fully masks
#' the real one. Doubling braces (`{` -> `{{`, `}` -> `}}`) is glue's own
#' escape convention for a literal brace -- this makes the substitution a
#' single, safe pass: the escaped text is inserted verbatim, never
#' re-parsed for further `{...}` expressions.
#' @return Character scalar, brace-escaped.
#' @keywords internal
.cli_escape <- function(x) {
  gsub("\\}", "}}", gsub("\\{", "{{", x))
}

#' Call an external resource with a documented, shared failure policy
#'
#' @description
#' Every function in `patliR` that talks to an external API or web resource
#' (PubChem, ChEMBL, COCONUT, Open Targets, ...) goes through this helper,
#' so that failure handling is consistent and never silently swallowed.
#'
#' @param fetch_fun A zero-argument function (typically a closure) that
#'   performs the actual request and returns the parsed result.
#' @param cache_dir Directory to store/read the `.rds` cache for this
#'   specific request.
#' @param cache_key Character scalar identifying this request (used to name
#'   the cache file, and in log/error messages).
#' @param mode `"abort"`: no fallback, fail loudly if the request fails --
#'   use this when there is no sensible way to continue without the data.
#'   `"warn_and_cache"`: use the local cache if one exists, otherwise warn
#'   and return `NULL` -- use this for optional enrichment steps that should
#'   never break the pipeline.
#' @return The fetched result, the cached result, or `NULL` (only reachable
#'   with `mode = "warn_and_cache"`).
#' @keywords internal
.fetch_external <- function(fetch_fun, cache_dir, cache_key,
                             mode = c("abort", "warn_and_cache")) {
  mode <- match.arg(mode)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_path <- file.path(cache_dir, paste0(cache_key, ".rds"))

  result <- tryCatch(fetch_fun(), error = function(e) e)

  if (!inherits(result, "error")) {
    saveRDS(result, cache_path)
    return(result)
  }

  if (mode == "abort") {
    cli::cli_abort(c(
      "Could not reach the external resource for {.val {cache_key}}.",
      "x" = .cli_escape(conditionMessage(result))
    ))
  }

  if (file.exists(cache_path)) {
    cli::cli_warn(c(
      "Could not reach the external resource for {.val {cache_key}}; using the local cache instead.",
      "x" = .cli_escape(conditionMessage(result))
    ))
    return(readRDS(cache_path))
  }

  cli::cli_warn(c(
    "Could not reach the external resource for {.val {cache_key}}, and no local cache exists.",
    "x" = .cli_escape(conditionMessage(result)),
    "i" = "Continuing without this data; re-run later once the resource is reachable."
  ))
  NULL
}

#' Parse SMILES strings into rcdk molecule objects, one at a time
#' @return A list the same length as `smiles`, each element either an rcdk
#'   `IAtomContainer` (jobjRef) or `NULL` if parsing failed.
#' @keywords internal
.parse_smiles_safe <- function(smiles) {
  lapply(smiles, function(s) {
    if (is.na(s) || !nzchar(s)) return(NULL)
    tryCatch(rcdk::parse.smiles(s)[[1]], error = function(e) NULL)
  })
}

#' Which SMILES strings parse into a chemically sane (non-empty) structure
#' @return Logical vector, same length as `smiles`.
#' @keywords internal
.validate_smiles <- function(smiles) {
  .check_structures(smiles)$valid
}

#' Compute a canonical-SMILES structure key for a vector of SMILES strings
#' @return Character vector, same length as `smiles`; `NA_character_` where
#'   parsing/canonicalization failed.
#' @keywords internal
.compute_structure_key <- function(smiles) {
  .check_structures(smiles)$canonical_smiles
}

#' Parse, validate, and canonicalize a vector of SMILES in one pass
#'
#' @description
#' The structure-identity key is a CDK canonical SMILES (`Canonical` +
#' `UseAromaticSymbols` + `Stereo` flavor), not an InChIKey -- `rcdk` on
#' CRAN has no InChIKey function. Deterministic within CDK, but **not**
#' guaranteed to match the canonical SMILES another toolkit (RDKit,
#' OpenBabel, PubChem, ChEMBL) produces for the same molecule;
#' cross-platform joins (`adme_import()`, `refdb_build()`) account for
#' this. See `DESIGN.md`.
#'
#' @section Stereochemistry is part of the identity key:
#' The `Stereo` flavor bit is required, not optional -- without it, CDK's
#' SMILES writer drops chirality (`@`/`@@`) and double-bond geometry
#' (`/`/`\\`) markers, so two true enantiomers (or E/Z isomers) parse to
#' the *same* canonical SMILES and [prep_compounds()]'s
#' `duplicated(canonical_smiles)` dedup (`dedup = TRUE`, the default)
#' silently drops one of them as if it were the same compound. Confirmed
#' directly with CDK: `(R)`-2-chlorobutane and `(S)`-2-chlorobutane both
#' canonicalize to `"ClC(C)CC"` without `Stereo`; with it, to
#' `"CC[C@H](C)Cl"` / `"CC[C@@@@H](C)Cl"`. This key still does not
#' distinguish isotopes, salts, or protonation states -- those are a
#' separate, still-open question of what "same molecule" should mean for
#' deduplication (see `DESIGN.md`/`ROADMAP.md`); this fix only restores
#' the two stereo descriptors CDK is capable of writing.
#'
#' Also keeps the *real* rcdk/rJava error message instead of collapsing
#' every failure into a generic "invalid" flag -- silently swallowing the
#' underlying Java exception is exactly what made the InChIKey bug above
#' impossible to diagnose from the log alone.
#'
#' @return A `data.frame` with columns `valid` (logical), `canonical_smiles`
#'   (character, `NA` if invalid), and `error` (character, `NA` if valid --
#'   otherwise the actual condition message from rcdk/rJava).
#' @keywords internal
.check_structures <- function(smiles) {
  flavor <- rcdk::smiles.flavors(c("Canonical", "UseAromaticSymbols", "Stereo"))
  rows <- lapply(smiles, function(s) {
    if (is.na(s) || !nzchar(s)) {
      return(data.frame(valid = FALSE, canonical_smiles = NA_character_, error = "empty SMILES", stringsAsFactors = FALSE))
    }
    mol <- tryCatch(rcdk::parse.smiles(s)[[1]], error = function(e) e)
    if (inherits(mol, "error")) {
      return(data.frame(valid = FALSE, canonical_smiles = NA_character_,
                         error = paste("parse.smiles() failed:", conditionMessage(mol)), stringsAsFactors = FALSE))
    }
    if (is.null(mol)) {
      return(data.frame(valid = FALSE, canonical_smiles = NA_character_,
                         error = "parse.smiles() returned NULL (unparseable SMILES)", stringsAsFactors = FALSE))
    }
    n_atoms <- tryCatch(rcdk::get.atom.count(mol), error = function(e) e)
    if (inherits(n_atoms, "error") || is.null(n_atoms) || n_atoms == 0) {
      msg <- if (inherits(n_atoms, "error")) paste("get.atom.count() failed:", conditionMessage(n_atoms)) else "structure has 0 atoms"
      return(data.frame(valid = FALSE, canonical_smiles = NA_character_, error = msg, stringsAsFactors = FALSE))
    }
    key <- tryCatch(rcdk::get.smiles(mol, flavor), error = function(e) e)
    if (inherits(key, "error")) {
      return(data.frame(valid = FALSE, canonical_smiles = NA_character_,
                         error = paste("get.smiles() failed:", conditionMessage(key)), stringsAsFactors = FALSE))
    }
    data.frame(valid = TRUE, canonical_smiles = key, error = NA_character_, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Read a CSV from an external/user-provided file, without the cosmetic
#' "incomplete final line" warning drowning out real problems
#'
#' @description
#' Real exports from third-party platforms (SuperPred, SwissADME, ...)
#' routinely lack a trailing newline on the last line -- harmless (R still
#' reads every row correctly), but [utils::read.csv()] raises a warning for
#' it every time, indistinguishable at a glance from a warning that
#' actually matters. This wraps [utils::read.csv()] and muffles only that
#' specific warning, letting every other warning (malformed quoting,
#' embedded nulls, ...) through unchanged.
#' @return A `data.frame`, as [utils::read.csv()] with
#'   `stringsAsFactors = FALSE, check.names = FALSE`.
#' @keywords internal
.read_csv_safe <- function(path) {
  withCallingHandlers(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    warning = function(w) {
      if (grepl("incomplete final line", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Find a column by name, tolerant of case, surrounding whitespace, and a
#' leading/trailing decoration (`*`, `.`, `_`) real platform exports add
#'
#' @description
#' Real exports from prediction/ADME platforms are not always consistent
#' about a column's exact spelling across runs (seen on real SuperPred
#' data: `"Probability"` vs. `"*Probability"`). Tries an exact match first
#' (so a real column named e.g. `"probability_2"` is never mistaken for
#' `"Probability"`), then falls back to comparing both names with
#' everything but letters/digits stripped and lower-cased.
#' @param available Character vector of column names actually present.
#' @param target The column name being looked for.
#' @return The matching element of `available`, or `NA_character_` if none
#'   matched even loosely.
#' @keywords internal
.match_column_flexible <- function(available, target) {
  if (target %in% available) return(target)
  norm <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
  hit <- which(norm(available) == norm(target))
  if (length(hit) > 0) return(available[hit[1]])
  NA_character_
}

#' Parse a numeric column that may be a plain fraction (`0.92`) or a
#' percentage (`"96.55%"`, `"85"`), always returning a `0..1` fraction
#'
#' @description
#' Real target-prediction exports mix conventions: some give probability
#' as a `0..1` fraction, others as a `"96.55%"` string (fails
#' [as.numeric()] outright because of the `%`), others as a bare `85`
#' meaning 85%. Strips `%` and any other non-numeric decoration first,
#' then rescales to `0..1` if the column looks like it was on a `0..100`
#' scale (either it had a literal `%` somewhere, or its max parsed value
#' is `> 1`) -- decided once per column, not per row, so one column never
#' ends up a mix of both scales.
#' @param x Character or numeric vector.
#' @return Numeric vector in `0..1` (`NA` where unparseable), same length
#'   as `x`.
#' @keywords internal
.parse_percent_column <- function(x) {
  x_chr <- as.character(x)
  had_percent <- grepl("%", x_chr, fixed = TRUE)
  cleaned <- gsub("[^0-9eE+.\\-]", "", x_chr)
  val <- suppressWarnings(as.numeric(cleaned))
  is_percent_scale <- any(had_percent, na.rm = TRUE) || (any(!is.na(val)) && max(val, na.rm = TRUE) > 1)
  if (is_percent_scale) val / 100 else val
}

#' Export a compound set as a plain SMILES list plus a row-order bridge CSV
#'
#' Shared body of [adme_export_smiles()] and [tox_export_smiles()]: one
#' canonical SMILES per line in a fixed order, plus a `data.frame` mapping
#' each line back to `compound_id`/`name` (external platforms only ever
#' order their output by input position, with no id of their own).
#'
#' @param cmp A `compounds(proj)` data frame, already subset to the
#'   compounds to export.
#' @param out_file Base path, or `NULL` to only return the result.
#' @return `list(smiles_text, mapping)`; `mapping` has columns `row_order`,
#'   `compound_id`, `name`, `smiles`.
#' @keywords internal
.export_smiles_list <- function(cmp, out_file = NULL) {
  if (nrow(cmp) == 0) {
    cli::cli_abort("No matching compounds in {.arg proj}; run {.fn prep_compounds} first.")
  }
  smiles <- ifelse(!is.na(cmp$canonical_smiles) & nzchar(cmp$canonical_smiles),
                   cmp$canonical_smiles, cmp$smiles)
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
  if (!is.null(out_file)) {
    writeLines(mapping$smiles, paste0(out_file, ".txt"))
    utils::write.csv(mapping, paste0(out_file, "_map.csv"), row.names = FALSE)
  }
  invisible(list(smiles_text = paste(mapping$smiles, collapse = "\n"), mapping = mapping))
}

#' Match imported-platform rows to compound ids via an export bridge CSV
#'
#' The reliable half of the [adme_export_smiles()]/[tox_export_smiles()]
#' round trip: rows come back from the platform in the same order they were
#' pasted in, so row `i` of the export is `mapping_file`'s `row_order == i`.
#' No SMILES canonicalization involved.
#'
#' @param n_rows Number of data rows in the platform export.
#' @param mapping_file Path to the `<out_file>_map.csv` written by the
#'   matching `*_export_smiles()` call.
#' @return Character vector of length `n_rows`; `NA` for any row past the
#'   end of the mapping.
#' @keywords internal
.match_by_export_mapping <- function(n_rows, mapping_file) {
  if (!file.exists(mapping_file)) {
    cli::cli_abort("{.arg mapping_file} {.path {mapping_file}} does not exist.")
  }
  map <- utils::read.csv(mapping_file, stringsAsFactors = FALSE)
  if (!all(c("row_order", "compound_id") %in% names(map))) {
    cli::cli_abort("{.arg mapping_file} must have {.val row_order} and {.val compound_id} columns (written by {.fn adme_export_smiles}/{.fn tox_export_smiles}).")
  }
  if (n_rows != nrow(map)) {
    cli::cli_warn(paste(
      "The export has {n_rows} row(s) but {.arg mapping_file} has {nrow(map)};",
      "matching by position -- any extra rows on either side stay unmatched."
    ))
  }
  map$compound_id[match(seq_len(n_rows), map$row_order)]
}

#' Next block of sequential internal compound ids, continuing from existing
#' @param existing_ids Character vector of already-used ids, e.g. `"C0007"`.
#' @param n How many new ids to generate.
#' @return Character vector of length `n`, e.g. `c("C0008", "C0009")`.
#' @keywords internal
.next_compound_ids <- function(existing_ids, n) {
  used_numbers <- suppressWarnings(as.integer(sub("^C", "", existing_ids)))
  start <- if (length(used_numbers) == 0 || all(is.na(used_numbers))) {
    1L
  } else {
    max(used_numbers, na.rm = TRUE) + 1L
  }
  sprintf("C%04d", seq.int(start, length.out = n))
}
