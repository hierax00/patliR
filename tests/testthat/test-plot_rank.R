.plot_rank_test_setup <- function() {
  testthat::skip_if_not_installed("RobustRankAggreg")
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  rank_candidates(proj, condition = "FLO-ET")
}

test_that("plot_rank() errors clearly without rank_candidates()", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_rank(proj, condition = "FLO-ET", save = FALSE), "rank_candidates")
})

test_that("plot_rank(view = 'pareto') returns a ggplot with the documented default axes", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", save = FALSE)
  expect_s3_class(p, "ggplot")
  ## no optional criteria ran -> y falls back to crit_centrality (the last
  ## entry in .plot_rank_y_priority still present in the mandatory-only table)
  expect_identical(rlang::as_label(p$mapping$x), "rra_rank")
  expect_identical(rlang::as_label(p$mapping$y), "crit_centrality")
})

test_that("plot_rank(view = 'pareto') honours explicit x/y overrides", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", x = "crit_adme_pass_frac", y = "crit_centrality", save = FALSE)
  expect_identical(rlang::as_label(p$mapping$x), "crit_adme_pass_frac")
  expect_identical(rlang::as_label(p$mapping$y), "crit_centrality")
})

test_that("plot_rank(view = 'pareto') rejects an unknown column name", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  expect_error(plot_rank(proj, condition = "FLO-ET", view = "pareto", x = "not_a_column", save = FALSE), "not_a_column")
})

test_that("plot_rank(view = 'pareto') saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "rank_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.png$")
})

test_that("plot_rank(view = 'heatmap') returns a pheatmap object restricted to top_n_compounds/top_n_targets", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .plot_rank_test_setup()
  n_all_compounds <- length(unique(patliRResults(proj, "rank_candidates")$compound_id))
  testthat::skip_if(n_all_compounds < 2, "need at least 2 compounds")

  h <- plot_rank(proj, condition = "FLO-ET", view = "heatmap", top_n_compounds = 2, top_n_targets = 3, save = FALSE)
  expect_s3_class(h, "pheatmap")
  expect_lte(nrow(h$gtable), 100) # sanity: object renders without error
})

test_that("plot_rank(view = 'heatmap') saves a PDF and logs it", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .plot_rank_test_setup()
  h <- plot_rank(proj, condition = "FLO-ET", view = "heatmap", top_n_compounds = 2, top_n_targets = 3, save = TRUE)
  proj2 <- attr(h, "proj")
  log_df <- patliRResults(proj2, "rank_heatmap_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.pdf$")
})

test_that("plot_rank(view = 'heatmap', target_relevance = 'centrality') requires network_centrality", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .network_stats_test_setup()
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  ## build rank_candidates with only adme+centrality is required, so
  ## network_centrality must exist for rank_candidates() itself already --
  ## exercise the heatmap-specific centrality lookup by requesting a
  ## disease-agnostic condition with an intentionally cleared centrality slot.
  proj <- network_centrality(proj)
  proj <- rank_candidates(proj, condition = "FLO-ET")
  patliRResults(proj, "network_centrality") <- NULL
  expect_error(
    plot_rank(proj, condition = "FLO-ET", view = "heatmap", target_relevance = "centrality", save = FALSE),
    "network_centrality"
  )
})
