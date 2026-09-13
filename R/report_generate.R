#' @include AllGenerics.R internal.R network_build.R adme_filter.R bias_audit.R network_module_robustness.R rank_candidates.R
NULL

## v1, deliberately plain: a self-contained HTML file built with base R
## string concatenation, no rmarkdown/pandoc. Two reasons, not one: (1) this
## dev environment has no pandoc (confirmed via
## rmarkdown::pandoc_available() before writing this), so an Rmd-based
## report could not even be rendered/tested here; (2) DESIGN.md's own
## "minimal Suggests" philosophy already prefers hand-rolling a simple
## rendering over a heavier dependency when the geometry (here: a handful
## of HTML tables) does not need one. Every section auto-includes when its
## source table exists in patliRResults(proj) and shows a one-line
## placeholder otherwise -- the same "no silent behaviour change, always
## visible" convention rank_candidates() already established, not a hard
## prerequisite abort.

#' Generate a self-contained HTML report per condition
#'
#' @description
#' One plain HTML file per condition, covering compounds present, ADME
#' filtering, the database-bias audit, modular robustness, the
#' [rank_candidates()] ranking, and the full [projectLog()] -- each section
#' auto-included when its source table exists in [patliRResults()] and
#' replaced by a one-line "not run" placeholder otherwise, so a project that
#' skipped a step still gets a report, just a shorter one. Deliberately
#' plain HTML (base R string building, no `rmarkdown`/`pandoc`, no CSS
#' framework) -- see the file header for why.
#'
#' @inheritParams network_build
#' @param out_dir Directory to write the report(s) to. `NULL` (default) ->
#'   `file.path(projectDir(proj), "reports")`.
#' @param top_n How many rows of the bias-audit "promiscuous" table and the
#'   [rank_candidates()] ranking to show. Default `15`.
#'
#' @return The updated `proj`, with a `report_log` entry in
#'   [patliRResults()] (columns `condition`, `path`), also written to
#'   `results/report_log.csv`. Each condition's HTML file is written to
#'   `out_dir/report_<condition>.html`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... build the project, run whichever of adme_local/adme_filter,
#' # bias_audit/bias_reweight, network_module_robustness, rank_candidates ...
#' proj <- report_generate(proj)
#' patliRResults(proj, "report_log")
#' }
#'
#' @export
report_generate <- function(proj, condition = NULL, out_dir = NULL, top_n = 15) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)

  conditions <- .network_resolve_conditions(proj, condition)
  if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "reports")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  paths <- character(length(conditions))
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    html <- .report_render_condition(proj, cond, top_n)
    path <- file.path(out_dir, paste0("report_", cond, ".html"))
    writeLines(html, path)
    paths[i] <- path
    proj <- .log_append(
      proj, step = "report_generate", id = NA_character_,
      message = paste0("condition '", cond, "': report written to '", path, "'")
    )
  }

  log_row <- data.frame(condition = conditions, path = paths, stringsAsFactors = FALSE)
  log_df <- .network_upsert(proj, "report_log", log_row, "condition")
  patliRResults(proj, "report_log") <- log_df
  .write_results_csv(proj, "report_log", log_df)
  .write_log_csv(proj)
  proj
}

#' Build one condition's report as a single HTML string
#' @keywords internal
.report_render_condition <- function(proj, cond, top_n) {
  cmp <- compounds(proj)
  edges <- patliRResults(proj, "network_edges")
  present_ids <- if (!is.null(edges)) sort(unique(edges$compound_id[edges$condition == cond])) else character(0)
  cmp_cond <- cmp[cmp$id %in% present_ids, , drop = FALSE]

  header <- sprintf(
    "<h1>patliR report -- condition '%s'</h1>\n<p>Generated %s by patliR %s. Project: <code>%s</code>.</p>",
    .report_html_escape(cond), .report_html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    as.character(utils::packageVersion("patliR")), .report_html_escape(projectDir(proj))
  )

  compounds_section <- .report_section(
    sprintf("Compounds present (%d)", nrow(cmp_cond)),
    .report_html_table(cmp_cond[, intersect(c("id", "name", "pubchem_id", "smiles"), names(cmp_cond)), drop = FALSE])
  )

  adme <- patliRResults(proj, "adme_filtered")
  adme_body <- if (is.null(adme) || nrow(adme) == 0) {
    .report_not_run("adme_local()/adme_filter()")
  } else {
    adme_sub <- adme[adme$compound_id %in% present_ids, , drop = FALSE]
    if (nrow(adme_sub) == 0) {
      .report_not_run("adme_filter() for this condition's compounds")
    } else {
      agg <- stats::aggregate(pass ~ rule, adme_sub, function(x) sprintf("%d/%d pass", sum(x), length(x)))
      names(agg) <- c("rule", "pass_fraction")
      .report_html_table(agg)
    }
  }
  adme_section <- .report_section("ADME filtering", adme_body)

  reweighted <- patliRResults(proj, "bias_reweighted")
  bias_body <- if (is.null(reweighted) || nrow(reweighted) == 0) {
    .report_not_run("bias_audit()/bias_reweight()")
  } else {
    rep <- bias_report(proj)
    paste0("<p>", .report_html_escape(rep$note), "</p>", .report_html_table(rep$summary, max_rows = top_n))
  }
  bias_section <- .report_section("Database-bias audit", bias_body)

  rob <- patliRResults(proj, "network_module_robustness")
  rob_body <- if (is.null(rob) || nrow(rob) == 0) {
    .report_not_run("network_module_robustness()")
  } else {
    rob_cond <- rob[rob$condition == cond, , drop = FALSE]
    if (nrow(rob_cond) == 0) {
      .report_not_run("network_module_robustness() for this condition")
    } else {
      cols <- intersect(c("module_id", "module_type", "n_nodes", "r_index", "r_index_random", "clustering"), names(rob_cond))
      .report_html_table(rob_cond[, cols, drop = FALSE])
    }
  }
  robustness_section <- .report_section("Modular robustness", rob_body)

  rc <- patliRResults(proj, "rank_candidates")
  rank_body <- if (is.null(rc) || nrow(rc) == 0) {
    .report_not_run("rank_candidates()")
  } else {
    rc_cond <- rc[rc$condition == cond, , drop = FALSE]
    if (nrow(rc_cond) == 0) {
      .report_not_run("rank_candidates() for this condition")
    } else {
      rc_cond <- rc_cond[order(rc_cond$rra_rank), , drop = FALSE]
      .report_html_table(rc_cond, max_rows = top_n)
    }
  }
  rank_section <- .report_section(sprintf("Candidate ranking (top %d by rra_rank)", top_n), rank_body)

  log_section <- .report_section("Full decision / parameter / seed log", .report_html_table(projectLog(proj), max_rows = 1000))

  paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>patliR report -- ",
    .report_html_escape(cond), "</title>",
    "<style>body{font-family:sans-serif;margin:2em;max-width:1100px;}",
    "table{border-collapse:collapse;margin-bottom:1.5em;width:100%;}",
    "th,td{padding:4px 8px;border:1px solid #ccc;font-size:0.85em;text-align:left;}",
    "th{background:#f0f0f0;}</style></head><body>",
    header, compounds_section, adme_section, bias_section, robustness_section, rank_section, log_section,
    "</body></html>"
  )
}

#' Wrap a report body in an `<h2>` section heading
#' @keywords internal
.report_section <- function(title, body_html) {
  sprintf("<h2>%s</h2>\n%s\n", .report_html_escape(title), body_html)
}

#' A one-line "this step has not been run" placeholder for a report section
#' @keywords internal
.report_not_run <- function(step_label) {
  sprintf("<p><em>%s has not been run (or produced nothing for this condition).</em></p>", .report_html_escape(step_label))
}

#' Escape `&`/`<`/`>` for safe embedding in the report's HTML (order matters
#' -- `&` first)
#' @keywords internal
.report_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  ifelse(is.na(x), "", x)
}

#' Render a `data.frame` as a plain HTML `<table>`
#' @param max_rows `NULL` (default, every row) or a cap (`utils::head()`).
#' @keywords internal
.report_html_table <- function(df, max_rows = NULL) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return("<p><em>No rows.</em></p>")
  if (!is.null(max_rows) && nrow(df) > max_rows) df <- utils::head(df, max_rows)

  fmt_cell <- function(v) {
    if (is.numeric(v) && !is.na(v)) formatC(v, digits = 4, format = "g") else .report_html_escape(v)
  }
  header <- paste0("<tr>", paste0("<th>", .report_html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body_rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(df[i, , drop = FALSE], fmt_cell, character(1))
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1))
  paste0("<table>", header, paste(body_rows, collapse = ""), "</table>")
}
