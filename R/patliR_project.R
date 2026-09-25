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
    compounds(proj) <- utils::read.csv(
      compounds_path, stringsAsFactors = FALSE, colClasses = .patliR_compounds_colclasses
    )
  }

  matrix_raw_path <- file.path(project_dir, "02_matrix_raw.csv")
  if (file.exists(matrix_raw_path)) {
    matrixRaw(proj) <- utils::read.csv(
      matrix_raw_path, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = .patliR_compound_id_colclass
    )
  }

  binarized_path <- file.path(project_dir, "03_binarized.csv")
  if (file.exists(binarized_path)) {
    binarizedMatrix(proj) <- utils::read.csv(
      binarized_path, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = .patliR_compound_id_colclass
    )
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
      ## read.csv() silently ignores a colClasses name that is one of
      ## several column names present in the file, but warns
      ## "not all columns named in 'colClasses' exist" when NONE of them
      ## are (e.g. disease_genes.csv, which has neither string_version nor
      ## condition) -- harmless, but noisy on every patliR_load() of such a
      ## table. Pre-filtering to the header's actual names avoids it
      ## without changing which columns get pinned.
      header <- names(utils::read.csv(f, nrows = 0, check.names = FALSE))
      colclasses <- .patliR_results_colclasses[names(.patliR_results_colclasses) %in% header]
      patliRResults(proj, entry_name) <- utils::read.csv(
        f, stringsAsFactors = FALSE, colClasses = colclasses
      )
    }
  }

  proj
}

#' Columns of any `results/*.csv` that are character *by contract* and must
#' not be re-typed by `utils::read.csv()`'s value inference
#'
#' @description
#' Most character results columns hold non-numeric-looking values
#' (`"C0001"`, `"P12345"`, `"EFO_0000537"`, `"P2"`, `"leiden"`) and survive
#' a `read.csv()` round-trip unchanged. A few do not: `string_version`
#' holds `"12.0"` / `"11.5"`, which `read.csv()` reads back as a numeric,
#' and `condition` is a user-chosen label that could be all digits. When
#' such a column is later merged by [network_build()]'s `.network_upsert()`
#' against a freshly-built (still character) column of the same name, the
#' type mismatch aborts the upsert. Pinning them here fixes it at the read
#' boundary; `read.csv()` silently ignores names not present in a given
#' file, so one shared vector covers every results CSV.
#' @keywords internal
.patliR_results_colclasses <- c(
  string_version = "character",
  version        = "character",   # network_bowtie_summary: STRING release "11.0"
  condition      = "character"
)

#' Identifier columns of `01_compounds.csv` that must survive
#' [patliR_load()] as `character`, never re-typed by `read.csv()`'s value
#' inference
#'
#' @description
#' All six are validated as required by [`PatliRProject-class`] and every
#' one of them can look purely numeric for a given project and so get
#' silently re-typed on reload if not pinned: `pubchem_id` is usually all
#' digits (becomes `integer`, dropping any leading zero a source database
#' used), and `name`/`smiles`/`canonical_smiles`/`source` are free text
#' that *happens* to contain no compound named e.g. `"TRUE"`/`"2e5"` only
#' by chance for a given dataset -- `read.csv()` infers a column's type
#' from the values actually present, not from any contract on what the
#' column is supposed to hold. `id` (`"C0001"`-style, see
#' `.next_compound_ids()`) never collides with a purely-numeric parse
#' because of its letter prefix, but is included here too so this list is
#' the full, load-bearing schema rather than "whatever happened not to
#' need it yet".
#' @keywords internal
.patliR_compounds_colclasses <- c(
  id = "character", pubchem_id = "character", name = "character",
  smiles = "character", canonical_smiles = "character", source = "character"
)

#' `compound_id` in `02_matrix_raw.csv`/`03_binarized.csv` must survive
#' [patliR_load()] as `character` for the same reason as
#' `.patliR_compounds_colclasses()` -- these two tables otherwise have
#' one column per experimental *condition* (user-chosen names,
#' `check.names = FALSE` on purpose), so only the one column name common
#' to every project can be pinned generically here.
#' @keywords internal
.patliR_compound_id_colclass <- c(compound_id = "character")

.patliR_version <- function() {
  tryCatch(as.character(utils::packageVersion("patliR")), error = function(e) "0.0.0.dev")
}
