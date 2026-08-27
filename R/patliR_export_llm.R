#' @include AllGenerics.R internal.R patliR_project.R
NULL

## patliR_export_llm() -- Uriel asked directly in chat for "un export del
## objeto que permita subírselo a un LLM". Deliberately generic over
## patliRResults() (loops over whatever is actually present) rather than
## hardcoding a bespoke summary per family -- keeps this function correct
## automatically as new network_*/bias_*/etc. results get added, at the
## cost of a plainer, less curated narrative than a bespoke report would
## give. That curated narrative is exactly what the still-unimplemented
## report_generate() and "narrador" module (both in patliR_manual.md,
## "Planeado / no implementado") are for -- this function is a raw,
## complete data dump for an LLM to read and reason over itself, not a
## pre-written summary.

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
