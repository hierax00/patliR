## report_generate() is deliberately plain HTML (base R, no rmarkdown/
## pandoc -- see R/report_generate.R's file header). These tests check the
## graceful-degradation contract (each section appears or shows a
## placeholder depending on what has been run) and the file/log mechanics,
## not HTML rendering fidelity.

.report_test_setup <- function() {
  .network_stats_test_setup()
}

test_that("report_generate() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(report_generate(proj), "network_edges")
})

test_that("report_generate() writes one HTML file per condition and logs it", {
  proj <- .report_test_setup()
  conditions <- unique(patliRResults(proj, "network_edges")$condition)
  proj <- report_generate(proj)
  log_df <- patliRResults(proj, "report_log")

  expect_setequal(log_df$condition, conditions)
  for (p in log_df$path) {
    expect_true(file.exists(p))
    expect_match(p, "\\.html$")
  }
})

test_that("report_generate() shows a placeholder for every section when nothing else has run", {
  proj <- .report_test_setup()
  proj <- report_generate(proj, condition = "FLO-ET")
  html <- readLines(patliRResults(proj, "report_log")$path[1])
  txt <- paste(html, collapse = "\n")

  expect_match(txt, "Compounds present")
  expect_match(txt, "adme_local\\(\\)/adme_filter\\(\\) has not been run")
  expect_match(txt, "bias_audit\\(\\)/bias_reweight\\(\\) has not been run")
  expect_match(txt, "network_module_robustness\\(\\) has not been run")
  expect_match(txt, "rank_candidates\\(\\) has not been run")
  expect_match(txt, "Full decision / parameter / seed log")
})

test_that("report_generate() includes the ADME table once adme_filter() has run", {
  proj <- .report_test_setup()
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  proj <- report_generate(proj, condition = "FLO-ET")
  txt <- paste(readLines(patliRResults(proj, "report_log")$path[1]), collapse = "\n")

  expect_match(txt, "pass_fraction")
  expect_false(grepl("adme_local\\(\\)/adme_filter\\(\\) has not been run", txt))
})

test_that("report_generate() includes the rank_candidates table, capped at top_n, once it has run", {
  testthat::skip_if_not_installed("RobustRankAggreg")
  proj <- .report_test_setup()
  proj <- network_centrality(proj)
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  proj <- rank_candidates(proj, condition = "FLO-ET")
  n_compounds <- nrow(patliRResults(proj, "rank_candidates"))
  testthat::skip_if(n_compounds < 2, "need at least 2 compounds to test the top_n cap")

  proj <- report_generate(proj, condition = "FLO-ET", top_n = 1)
  txt <- paste(readLines(patliRResults(proj, "report_log")$path[1]), collapse = "\n")
  expect_match(txt, ">rra_rank<")

  ## the ranking section's table must show exactly top_n = 1 data row (1
  ## header <tr> + 1 body <tr>), not every compound in the project.
  section <- sub(".*(<h2>Candidate ranking.*?)<h2>Full decision.*", "\\1", txt)
  expect_equal(lengths(regmatches(section, gregexpr("<tr>", section))), 2L)
})

test_that("report_generate() escapes HTML-special characters in compound names", {
  proj <- .report_test_setup()
  cmp <- compounds(proj)
  cmp$name[1] <- "Foo <script>alert(1)</script> & Bar"
  compounds(proj) <- cmp

  proj <- report_generate(proj, condition = "FLO-ET")
  txt <- paste(readLines(patliRResults(proj, "report_log")$path[1]), collapse = "\n")

  expect_false(grepl("<script>alert", txt, fixed = TRUE))
  expect_match(txt, "&lt;script&gt;", fixed = TRUE)
  expect_match(txt, "&amp; Bar", fixed = TRUE)
})

test_that("report_generate() rebuilding one condition does not touch the others' report_log rows", {
  proj <- .report_test_setup()
  conditions <- unique(patliRResults(proj, "network_edges")$condition)
  testthat::skip_if(length(conditions) < 2, "fixture needs at least 2 conditions")

  proj <- report_generate(proj)
  before <- patliRResults(proj, "report_log")
  proj <- report_generate(proj, condition = conditions[1])
  after <- patliRResults(proj, "report_log")

  expect_setequal(after$condition, before$condition)
})

test_that(".report_html_table() renders an empty/degenerate data.frame without error", {
  expect_match(patliR:::.report_html_table(NULL), "No rows")
  expect_match(patliR:::.report_html_table(data.frame()), "No rows")
  expect_match(patliR:::.report_html_table(data.frame(a = character(0))), "No rows")
})
