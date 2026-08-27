#' @include AllClasses.R
NULL

#' Display a patliR project
#'
#' Prints a short, human-readable summary instead of dumping every slot --
#' full detail is always available via the accessors ([compounds()],
#' [matrixRaw()], [projectLog()], [patliRResults()], ...).
#'
#' @param object A [`PatliRProject-class`] object.
#' @return `object`, invisibly. Called for its side effect (printing).
#' @name show,PatliRProject-method
#' @aliases show
#' @importFrom methods show
setMethod("show", "PatliRProject", function(object) {
  ## Deliberately cat()/writeLines(), not cli::cli_*() -- show() methods
  ## display an object (real stdout output an R user or a captured script
  ## can rely on), while cli's alert/bullet functions route through R's
  ## condition system like message() does (stderr-like), which is right
  ## for diagnostics but wrong for "print my object".
  n_compounds  <- nrow(object@compounds)
  n_conditions <- if (nrow(object@matrix_raw) > 0) ncol(object@matrix_raw) else 0L
  n_results    <- length(object@results)
  n_log        <- nrow(object@log)
  results_detail <- if (n_results > 0) paste0(" (", paste(names(object@results), collapse = ", "), ")") else ""

  writeLines(c(
    "<PatliRProject>",
    paste0("  Project directory: ", object@project_dir),
    paste0("  Cache directory: ", object@cache_dir),
    paste0("  Compounds: ", n_compounds),
    paste0("  Conditions (abundance matrix): ", n_conditions),
    paste0("  Additional results computed: ", n_results, results_detail),
    paste0("  Log entries: ", n_log),
    paste0("  Created/updated with patliR ", object@version)
  ))

  invisible(object)
})
