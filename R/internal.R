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
#' These coordinates are digitized from the supporting information of
#' Daina, A. & Zoete, V. (2016), "A BOILED-Egg To Predict Gastrointestinal
#' Absorption and Brain Penetration of Small Molecules", *ChemMedChem*
#' 11, 1117-1121 -- taken from the GPL-3 implementation
#' <https://github.com/bfmilne/PyBOILEDegg> (Milne, B.F., 2021, Zenodo
#' \doi{10.5281/zenodo.4725530}), which states the same provenance. Used
#' as-is (TPSA on the x-axis, WLogP on the y-axis), not re-derived by us.
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

#' Write one of the numbered "source of truth" CSV files for a project
#' @keywords internal
.write_step_csv <- function(proj, filename, data) {
  path <- file.path(projectDir(proj), filename)
  utils::write.csv(data, path, row.names = FALSE)
  invisible(path)
}

#' Write one entry of the results bag to `results/<name>.csv`
#' @keywords internal
.write_results_csv <- function(proj, name, data) {
  results_dir <- file.path(projectDir(proj), "results")
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  utils::write.csv(data, file.path(results_dir, paste0(name, ".csv")), row.names = FALSE)
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
#' This used to compute an InChIKey via `rcdk::get.inchi.key()` -- which
#' turned out not to exist in `rcdk` at all (that function lives only in
#' the GitHub-only, non-CRAN `CDK-R/rinchi` package, so depending on it
#' would have blocked a Bioconductor submission outright; see the
#' `patliR_manual.md` decision log). We use a canonical SMILES (CDK's own
#' `Canonical` + `UseAromaticSymbols` flavor) as the structure-identity key
#' instead -- deterministic within `patliR`/CDK, but **not** guaranteed to
#' match the canonical SMILES another toolkit (RDKit, OpenBabel, PubChem,
#' ChEMBL) would generate for the same molecule. Cross-platform joins
#' (`adme_import()`, `refdb_build()`) account for this -- see their source.
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
  flavor <- rcdk::smiles.flavors(c("Canonical", "UseAromaticSymbols"))
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
