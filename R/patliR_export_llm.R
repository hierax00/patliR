#' @include AllGenerics.R internal.R patliR_project.R
NULL

## Loops over whatever patliRResults() currently holds rather than a fixed
## per-family summary, so it stays correct as new result types are added.
## A raw data dump, not a curated narrative (that is report_generate()'s
## job -- see ROADMAP.md).

#' Export a project's full state as plain text, for an LLM to read
#'
#' @description
#' Serializes [compounds()], every entry of [patliRResults()], and the
#' full [projectLog()] into one plain-text document (CSV-formatted
#' tables under Markdown-style headings) -- meant to be pasted into or
#' uploaded to an LLM chat for analysis/discussion of a project's
#' results, without needing R itself. Read-only: nothing is written to
#' `proj`, only optionally to `out_file`.
#'
#' @inheritParams compounds
#' @param out_file Character scalar path to write the text to, or `NULL`
#'   (default) to only return it.
#' @param max_rows Integer, default `200`. Tables longer than this are
#'   truncated (noted explicitly, never silently) so the export stays a
#'   reasonable size for an LLM context window.
#'
#' @return Invisibly, the full export as a single character string.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' txt <- patliR_export_llm(proj)
#' cat(substr(txt, 1, 200))
#' }
#'
#' @export
patliR_export_llm <- function(proj, out_file = NULL, max_rows = 200) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(max_rows), length(max_rows) == 1, max_rows >= 1)

  .df_to_text <- function(df) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return("(empty)")
    truncated <- nrow(df) > max_rows
    if (truncated) df <- utils::head(df, max_rows)
    out_lines <- NULL
    con <- textConnection("out_lines", "w", local = TRUE)
    utils::write.csv(df, con, row.names = FALSE)
    close(con)
    txt <- paste(out_lines, collapse = "\n")
    if (truncated) txt <- paste0(txt, "\n... (truncated at ", max_rows, " of ", nrow(df), "+ rows)")
    txt
  }

  results <- patliRResults(proj)

  parts <- c(
    "# patliR project export (for LLM analysis)",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Project directory: ", projectDir(proj)),
    paste0("patliR version: ", proj@version %||% "unknown"),
    "",
    "## Compounds",
    .df_to_text(compounds(proj)),
    ""
  )

  for (name in names(results)) {
    parts <- c(parts, paste0("## Result: ", name), .df_to_text(results[[name]]), "")
  }

  parts <- c(parts, "## Full project log (projectLog)", .df_to_text(projectLog(proj)), "")

  text <- paste(parts, collapse = "\n")
  if (!is.null(out_file)) writeLines(text, out_file)
  invisible(text)
}
