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
#' @param min_replicates `NULL` (default: the Q1 rule above) or a single
#'   integer: a compound is called present in a condition when it is detected
#'   (abundance > 0) in **at least** `min_replicates` of that condition's
#'   replicate columns, and absent otherwise -- e.g. `2` with three
#'   replicates means "absent if 2 or 3 replicates are zero, present if 2 or
#'   3 have a value". When set, `q` is ignored (no abundance threshold is
#'   applied) and the rule is recorded in the log.
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
                           average_replicates = TRUE, q = 0.25, min_replicates = NULL) {
  stopifnot(is(proj, "PatliRProject"), is.data.frame(data), nrow(data) > 0)
  stopifnot(is.numeric(q), length(q) == 1, q > 0, q < 1)
  if (!is.null(min_replicates)) {
    stopifnot(is.numeric(min_replicates), length(min_replicates) == 1, is.finite(min_replicates),
              min_replicates >= 1, min_replicates == round(min_replicates))
  }
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
    if (!is.null(min_replicates) && min_replicates > length(cols)) {
      cli::cli_abort("{.arg min_replicates} = {min_replicates} but condition {.val {cond}} has only {length(cols)} replicate column{?s}.")
    }

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
      if (!is.null(min_replicates)) bin <- ifelse(rowSums(vals > 0, na.rm = TRUE) >= min_replicates, 1L, 0L)
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
      used_q1 <- if (any(!is.na(avg))) stats::quantile(avg, probs = q, na.rm = TRUE, names = FALSE) else 0
      bin <- ifelse(avg > 0 & avg >= used_q1, 1L, 0L)
      if (!is.null(min_replicates)) {
        ## replicate-consistency rule: present iff at least `min_replicates`
        ## replicates are detected (> 0); the quartile threshold is not used
        bin <- ifelse(rowSums(vals > 0, na.rm = TRUE) >= min_replicates, 1L, 0L)
        used_q1 <- NA_real_
      }
    }

    matrix_raw[[cond]] <- avg
    binarized[[cond]]  <- bin
    log_rows[[cond]] <- data.frame(
      step = "prep_binarize", id = NA_character_,
      message = paste0(
        "condition '", cond, "': already_binary=", already_binary,
        ", n_replicates=", length(cols),
        if (!is.null(min_replicates)) paste0(", rule=detected_in_at_least_", min_replicates, "_replicates")
        else if (!already_binary) paste0(", q", q * 100, "_used=", signif(used_q1, 4)) else ""
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

#' Treat a compound list as a single experimental condition (no abundance matrix)
#'
#' @description
#' The `network_*` family keys everything off [binarizedMatrix()] -- one
#' network per condition column. When you have a plain compound list rather
#' than a GC-MS abundance matrix (Scenario B, [prep_compounds()]), there is
#' no matrix and no conditions, so `network_build()` has nothing to fork on.
#' `prep_as_condition()` fills that gap: it marks every compound (or a
#' subset) as present in one named condition, so the rest of the pipeline
#' -- `network_*`, `plot_*` -- runs on the list as if it were a single
#' extract. This is the entry point for a list-only, exploratory flow.
#'
#' @inheritParams compounds
#' @param condition Character scalar, the condition name to create (default
#'   `"all"`).
#' @param compound_ids Character vector of `compounds(proj)$id` to mark
#'   present, or `NULL` (default) for every compound currently in
#'   [compounds()].
#'
#' @return The updated `proj`, with [binarizedMatrix()] (and a matching
#'   [matrixRaw()] of all `1`s) carrying a single `condition` column, and
#'   `02_matrix_raw.csv` / `03_binarized.csv` written. Re-running with a
#'   different `condition` adds a column rather than replacing the matrix.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- prep_as_condition(proj, condition = "my_extract")
#' binarizedMatrix(proj)
#'
#' @seealso [prep_binarize()] for a real replicate-abundance matrix.
#' @export
prep_as_condition <- function(proj, condition = "all", compound_ids = NULL) {
  stopifnot(is(proj, "PatliRProject"), is.character(condition), length(condition) == 1, nzchar(condition))
  cmp <- compounds(proj)
  if (nrow(cmp) == 0) {
    cli::cli_abort("No compounds in {.arg proj}; run {.fn prep_compounds} first.")
  }
  ids <- cmp$id
  if (!is.null(compound_ids)) {
    unknown <- setdiff(compound_ids, ids)
    if (length(unknown) > 0) {
      cli::cli_warn("{length(unknown)} {.arg compound_ids} not in {.fn compounds} and ignored: {.val {unknown}}.")
    }
    ids <- intersect(ids, compound_ids)
  }
  if (length(ids) == 0) cli::cli_abort("No compounds left to mark present.")

  present <- as.integer(cmp$id %in% ids)

  raw <- matrixRaw(proj)
  bin <- binarizedMatrix(proj)
  if (nrow(raw) == 0) raw <- data.frame(compound_id = cmp$id, stringsAsFactors = FALSE)
  if (nrow(bin) == 0) bin <- data.frame(compound_id = cmp$id, stringsAsFactors = FALSE)
  if (condition %in% names(bin)) {
    cli::cli_warn("Condition {.val {condition}} already exists; overwriting it.")
  }
  missing_ids <- setdiff(cmp$id, raw$compound_id)
  if (length(missing_ids) > 0) {
    extra <- raw[rep(NA_integer_, length(missing_ids)), , drop = FALSE]
    extra$compound_id <- missing_ids
    raw <- rbind(raw, extra)
  }
  missing_ids <- setdiff(cmp$id, bin$compound_id)
  if (length(missing_ids) > 0) {
    extra <- bin[rep(NA_integer_, length(missing_ids)), , drop = FALSE]
    extra$compound_id <- missing_ids
    bin <- rbind(bin, extra)
  }
  raw[[condition]] <- as.integer(raw$compound_id %in% ids)
  bin[[condition]] <- as.integer(bin$compound_id %in% ids)

  ## set both slots together -- the S4 validity check requires matrix_raw
  ## and binarized to carry the same condition columns at all times.
  proj@matrix_raw <- raw
  proj@binarized <- bin
  methods::validObject(proj)
  proj <- .log_append(
    proj, step = "prep_as_condition", id = NA_character_,
    message = paste0("condition '", condition, "': ", sum(present), " of ", nrow(cmp), " compounds marked present (list-as-condition, no abundance matrix)")
  )
  .write_step_csv(proj, "02_matrix_raw.csv", raw)
  .write_step_csv(proj, "03_binarized.csv", bin)
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
