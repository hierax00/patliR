#' The patliR project object
#'
#' @description
#' `PatliRProject` is the single state-carrying object of the `patliR`
#' pipeline. It is deliberately split into two kinds of slots:
#'
#' * **Structural slots** (`compounds`, `matrix_raw`, `binarized`) hold the
#'   handful of tables that almost every family of functions depends on, and
#'   whose shape is stable across the whole pipeline. These have a fixed
#'   `data.frame` type and are checked by [methods::validObject()].
#' * **The `results` slot** is a named `list` "bag" that holds every other,
#'   more heterogeneous result the pipeline produces (ADME tables, target
#'   predictions, one network per condition, bias audits, rankings, ...).
#'   New families of functions (or future versions of `patliR`) can add new
#'   named entries to `results` without ever requiring a change to the class
#'   definition -- the same pattern
#'   `SummarizedExperiment::SummarizedExperiment` uses for its `metadata()`
#'   slot, applied here to the whole pipeline rather than to one assay.
#'
#' None of this replaces the CSV files on disk, which remain the durable
#' source of truth (see [patliR_project()] and [patliR_load()]); the S4
#' object is the in-memory, self-validating, pipeable mirror of those files
#' for the duration of one R session.
#'
#' @slot project_dir Character scalar. Path to the project directory where
#'   every pipeline step writes its numbered CSV output (e.g.
#'   `01_compounds.csv`).
#' @slot cache_dir Character scalar. Path to a directory for derived,
#'   always-regenerable caches (e.g. the `.rds` built by
#'   [refdb_rebuild_cache()]). Never required for the analysis to continue.
#' @slot compounds A `data.frame` with one row per compound and columns
#'   `id`, `pubchem_id`, `name`, `smiles`, `canonical_smiles`, `source`.
#'   `canonical_smiles` (CDK's own canonical form, via [prep_compounds()])
#'   is `patliR`'s structure-identity key for deduplication and joins --
#'   see the note in `.check_structures()` (internal) for why this is not
#'   an InChIKey.
#' @slot matrix_raw A `data.frame` with one row per compound and one column
#'   per condition, holding the (averaged-over-replicates) raw abundance
#'   values. Never binarized, kept so that later steps can still use
#'   relative concentration.
#' @slot binarized A `data.frame` with the same shape as `matrix_raw`, with
#'   0/1 presence/absence calls.
#' @slot results A named `list` holding every other pipeline result (see
#'   Description). Access it with [patliRResults()].
#' @slot log A `data.frame` with columns `step`, `id`, `message`, `timestamp`,
#'   accumulating one row per notable event (dropped compound, fetch
#'   failure, threshold used, ...) across every function that has touched
#'   this project. This is the backbone of the fully-logged report produced
#'   by `report_generate()`.
#' @slot version Character scalar. The `patliR` package version that created
#'   or last updated this object, for provenance.
#'
#' @seealso [patliR_project()] to create a new project, [patliR_load()] to
#'   reconstruct one from disk.
#' @name PatliRProject-class
#' @rdname PatliRProject-class
#' @exportClass PatliRProject
setClass(
  "PatliRProject",
  slots = c(
    project_dir = "character",
    cache_dir   = "character",
    compounds   = "data.frame",
    matrix_raw  = "data.frame",
    binarized   = "data.frame",
    results     = "list",
    log         = "data.frame",
    version     = "character"
  ),
  prototype = list(
    project_dir = NA_character_,
    cache_dir   = NA_character_,
    compounds   = data.frame(
      id = character(), pubchem_id = character(), name = character(),
      smiles = character(), canonical_smiles = character(), source = character(),
      stringsAsFactors = FALSE
    ),
    matrix_raw  = data.frame(),
    binarized   = data.frame(),
    results     = list(),
    log         = data.frame(
      step = character(), id = character(), message = character(),
      timestamp = as.POSIXct(character()), stringsAsFactors = FALSE
    ),
    version     = NA_character_
  )
)

#' @name PatliRProject-class
#' @rdname PatliRProject-class
setValidity("PatliRProject", function(object) {
  errors <- character()

  req_compound_cols <- c("id", "pubchem_id", "name", "smiles", "canonical_smiles", "source")
  if (nrow(object@compounds) > 0 && !all(req_compound_cols %in% names(object@compounds))) {
    errors <- c(errors, paste0(
      "'compounds' must contain columns: ", paste(req_compound_cols, collapse = ", "), "."
    ))
  }

  if (nrow(object@matrix_raw) > 0 && nrow(object@binarized) > 0) {
    if (!identical(sort(colnames(object@matrix_raw)), sort(colnames(object@binarized)))) {
      errors <- c(errors, "'matrix_raw' and 'binarized' must have the same condition columns.")
    }
  }

  req_log_cols <- c("step", "id", "message", "timestamp")
  if (!all(req_log_cols %in% names(object@log))) {
    errors <- c(errors, paste0(
      "'log' must contain columns: ", paste(req_log_cols, collapse = ", "), "."
    ))
  }

  if (length(errors) == 0) TRUE else errors
})
