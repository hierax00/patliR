#' Launch the patliR Shiny wizard
#'
#' @description
#' `launch_app()` starts an optional graphical wizard over the same
#' `prep_*`/`refdb_*`/`adme_*`/`tox_*`/`targets_*`/
#' `network_*`/`plot_*` functions already exported by `patliR` -- the GUI
#' orchestrates step by step, it never reimplements any analysis logic.
#' `library(patliR)` continues to work standalone without ever calling this
#' function; it is purely an optional convenience layer on top of the core
#' package, requiring the \pkg{shiny} package (Suggests).
#'
#' @param ... Passed on to [shiny::runApp()] (e.g. `launch.browser`, `port`).
#' @return Called for its side effect (blocks the R session while the app is
#'   running, same as any [shiny::runApp()] call). Invisibly returns `NULL`.
#'
#' @export
launch_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn launch_app} needs the {.pkg shiny} package.",
      "i" = "{.code install.packages(\"shiny\")}"
    ))
  }
  app_dir <- system.file("shiny", package = "patliR")
  if (!nzchar(app_dir)) {
    cli::cli_abort("Could not find the bundled Shiny app inside the installed {.pkg patliR} package.")
  }
  shiny::runApp(app_dir, ...)
  invisible(NULL)
}
