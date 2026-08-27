#' @include AllGenerics.R internal.R
NULL

#' Average replicates and binarize a raw abundance matrix (presence/absence)
#'
#' @description
#' `prep_binarize()` is the main entry point of the `patliR` pipeline for
#' **Scenario A**: a raw abundance matrix (e.g. GC-MS peak areas) with
#' compounds in rows and one column per replicate-condition combination. For
#' an already-curated compound table (**Scenario B**), use
#' [prep_compounds()] instead.
#'
#' Condition columns must follow the pattern `R<n>-<CONDITION>` (replicate
#' first, e.g. `R1-LEA-ET`, `R2-FLO-AQ`; the condition itself may contain
#' hyphens). Replicates are averaged per condition, and a per-condition
#' first-quartile (Q1) threshold is applied to the *averaged* values to call
#' presence/absence: a compound is called present (`1`) in a condition if
#' its average abundance is greater than 0 and at or above that condition's
#' Q1 across all compounds; otherwise absent (`0`). Columns that are already
#' binary (only `0`/`1` values) are detected and passed through unchanged
#' (no Q1 thresholding applied), and this is recorded in the log.
#'
#' @inheritParams compounds
#' @param data A `data.frame` with one identifier column (`id_col`) and one
#'   or more `R<n>-<CONDITION>` columns.
#' @param id_col Character scalar. Name of the compound identifier column in
#'   `data` (default `"Name"`). If [compounds()] already has entries, each
#'   value is matched against `compounds(proj)$name` (case- and
#'   whitespace-insensitive); unmatched rows are still included (using the
#'   raw identifier from `data`) but logged as `"compound_not_found"`, so
#'   this function never blocks on [prep_compounds()] not having run yet.
#' @param average_replicates Logical. If `TRUE` (default), replicate columns
#'   sharing the same condition are averaged. If `FALSE`, every condition
#'   must already have exactly one column (no `R<n>` replication) or this
#'   errors.
#' @param q Numeric scalar in `(0, 1)`. Quantile used as the per-condition
#'   presence threshold (default `0.25`, i.e. Q1).
#'
#' @return The updated `proj`, with [matrixRaw()] and [binarizedMatrix()]
#'   set (never overwriting one with the other) and `02_matrix_raw.csv` /
#'   `03_binarized.csv` written inside [projectDir()].
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' abundance <- read.csv(
#'   system.file("extdata", "input_abundance_matrix.csv", package = "patliR"),
#'   check.names = FALSE
#' )
#' proj <- prep_binarize(proj, abundance)
#' matrixRaw(proj)
#' binarizedMatrix(proj)
#'
#' @seealso [prep_compounds()] for Scenario B (already-curated compound
#'   tables).
#' @export
prep_binarize <- function(proj, data, id_col = "Name",
                           average_replicates = TRUE, q = 0.25) {
  stopifnot(is(proj, "PatliRProject"), is.data.frame(data), nrow(data) > 0)
  stopifnot(is.numeric(q), length(q) == 1, q > 0, q < 1)
  if (!id_col %in% names(data)) {
    cli::cli_abort("Column {.val {id_col}} (the {.arg id_col}) was not found in {.arg data}.")
  }

  ids <- as.character(data[[id_col]])
  value_cols <- setdiff(names(data), id_col)

  parsed <- .parse_replicate_columns(value_cols)
  conditions <- unique(parsed$condition)

  ## Match ids against proj@compounds$name, best-effort (never blocking).
  existing <- compounds(proj)
  matched_id <- ids
  if (nrow(existing) > 0) {
    norm <- function(x) trimws(tolower(x))
    hit <- match(norm(ids), norm(existing$name))
    found <- !is.na(hit)
    matched_id[found] <- existing$id[hit[found]]
    if (any(!found)) {
      proj <- .log_append(
        proj, step = "prep_binarize", id = ids[!found],
        message = "compound_not_found (no match in compounds(proj)$name; using raw identifier)"
      )
    }
  }

  matrix_raw <- data.frame(compound_id = matched_id, stringsAsFactors = FALSE)
  binarized  <- data.frame(compound_id = matched_id, stringsAsFactors = FALSE)
  log_rows <- list()

  for (cond in conditions) {
    cols <- value_cols[parsed$condition == cond]
    vals <- as.matrix(data[, cols, drop = FALSE])
    storage.mode(vals) <- "double"

    already_binary <- all(vals %in% c(0, 1, NA))

    if (already_binary) {
      if (length(cols) > 1) {
        if (!average_replicates) {
          cli::cli_abort(c(
            "Condition {.val {cond}} has {length(cols)} replicate columns but {.arg average_replicates} is {.val FALSE}.",
            "i" = "Set {.arg average_replicates = TRUE}, or provide exactly one column per condition."
          ))
        }
        avg <- rowMeans(vals, na.rm = TRUE)
      } else {
        avg <- vals[, 1]
      }
      ## Majority rule across already-binary replicates (>= half present),
      ## same "average then threshold" shape as the non-binary branch below
      ## -- a single replicate column used to silently swallow the others
      ## here before this fix.
      bin <- ifelse(avg >= 0.5, 1L, 0L)
      used_q1 <- NA_real_
    } else {
      if (length(cols) > 1) {
        if (!average_replicates) {
          cli::cli_abort(c(
            "Condition {.val {cond}} has {length(cols)} replicate columns but {.arg average_replicates} is {.val FALSE}.",
            "i" = "Set {.arg average_replicates = TRUE}, or provide exactly one column per condition."
          ))
        }
        avg <- rowMeans(vals, na.rm = TRUE)
      } else {
        avg <- vals[, 1]
      }
      nonzero <- avg[avg > 0]
      used_q1 <- if (length(nonzero) > 0) stats::quantile(nonzero, probs = q, na.rm = TRUE, names = FALSE) else 0
      bin <- ifelse(avg > 0 & avg >= used_q1, 1L, 0L)
    }

    matrix_raw[[cond]] <- avg
    binarized[[cond]]  <- bin
    log_rows[[cond]] <- data.frame(
      step = "prep_binarize", id = NA_character_,
      message = paste0(
        "condition '", cond, "': already_binary=", already_binary,
        ", n_replicates=", length(cols),
        if (!already_binary) paste0(", q", q * 100, "_used=", signif(used_q1, 4)) else ""
      ),
      timestamp = Sys.time(), stringsAsFactors = FALSE
    )
  }

  matrixRaw(proj)      <- matrix_raw
  binarizedMatrix(proj) <- binarized
  projectLog(proj) <- rbind(projectLog(proj), do.call(rbind, log_rows))

  .write_step_csv(proj, "02_matrix_raw.csv", matrix_raw)
  .write_step_csv(proj, "03_binarized.csv", binarized)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.parse_replicate_columns <- function(cols) {
  pattern <- "^R([0-9]+)-(.+)$"
  ok <- grepl(pattern, cols)
  if (!all(ok)) {
    cli::cli_abort(c(
      "These column names do not match the expected {.code R<n>-<CONDITION>} pattern:",
      "x" = "{.val {cols[!ok]}}",
      "i" = "Rename replicate/condition columns (e.g. {.val R1-LEA-ET}) or exclude non-value columns from {.arg data}."
    ))
  }
  data.frame(
    column = cols,
    replicate = sub(pattern, "\\1", cols),
    condition = sub(pattern, "\\2", cols),
    stringsAsFactors = FALSE
  )
}
