#' @include AllClasses.R
NULL

#' Create or reload a patliR project
#'
#' @description
#' `patliR_project()` creates a brand new, empty project rooted at
#' `project_dir`. `patliR_load()` reconstructs a project from the CSV files
#' a previous session left in `project_dir` -- this is what makes the
#' pipeline resumable without depending on the R session (or even the
#' `patliR` version) that produced those files staying alive. See
#' [`PatliRProject-class`] for what "reconstructed" means in detail.
#'
#' @param project_dir Character scalar. Path to the project directory.
#'   Created if it does not exist yet.
#' @param cache_dir Character scalar or `NULL`. Path to a directory for
#'   derived caches. Defaults to `file.path(project_dir, ".patliR_cache")`.
#'   Created if it does not exist yet. Never required for the analysis to
#'   continue -- see [refdb_rebuild_cache()].
#'
#' @return A new, empty [`PatliRProject-class`] object.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' proj
#'
#' @export
patliR_project <- function(project_dir, cache_dir = NULL) {
  stopifnot(is.character(project_dir), length(project_dir) == 1, nzchar(project_dir))

  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE)
  }
  if (is.null(cache_dir)) {
    cache_dir <- file.path(project_dir, ".patliR_cache")
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  methods::new(
    "PatliRProject",
    project_dir = normalizePath(project_dir, mustWork = TRUE),
    cache_dir   = normalizePath(cache_dir, mustWork = TRUE),
    version     = .patliR_version()
  )
}

#' @rdname patliR_project
#'
#' @details
#' `patliR_load()` looks for the following files inside `project_dir` and
#' reloads whichever of them are present (all are optional -- a project can
#' be reloaded at any point of the pipeline):
#' \itemize{
#'   \item `01_compounds.csv` -> [compounds()]
#'   \item `02_matrix_raw.csv` -> [matrixRaw()]
#'   \item `03_binarized.csv` -> [binarizedMatrix()]
#'   \item `patliR_log.csv` -> [projectLog()]
#'   \item any `<name>.csv` inside `<project_dir>/results/` -> that entry of
#'     [patliRResults()]
#' }
#'
#' @export
patliR_load <- function(project_dir, cache_dir = NULL) {
  stopifnot(is.character(project_dir), length(project_dir) == 1, dir.exists(project_dir))
  proj <- patliR_project(project_dir, cache_dir = cache_dir)

  compounds_path <- file.path(project_dir, "01_compounds.csv")
  if (file.exists(compounds_path)) {
    compounds(proj) <- utils::read.csv(compounds_path, stringsAsFactors = FALSE)
  }

  matrix_raw_path <- file.path(project_dir, "02_matrix_raw.csv")
  if (file.exists(matrix_raw_path)) {
    matrixRaw(proj) <- utils::read.csv(matrix_raw_path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  binarized_path <- file.path(project_dir, "03_binarized.csv")
  if (file.exists(binarized_path)) {
    binarizedMatrix(proj) <- utils::read.csv(binarized_path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  log_path <- file.path(project_dir, "patliR_log.csv")
  if (file.exists(log_path)) {
    log_df <- utils::read.csv(log_path, stringsAsFactors = FALSE)
    log_df$timestamp <- as.POSIXct(log_df$timestamp, tz = "UTC")
    projectLog(proj) <- log_df
  }

  results_dir <- file.path(project_dir, "results")
  if (dir.exists(results_dir)) {
    result_files <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)
    for (f in result_files) {
      entry_name <- tools::file_path_sans_ext(basename(f))
      patliRResults(proj, entry_name) <- utils::read.csv(f, stringsAsFactors = FALSE)
    }
  }

  proj
}

.patliR_version <- function() {
  tryCatch(as.character(utils::packageVersion("patliR")), error = function(e) "0.0.0.dev")
}
