test_that("plot_centrality() returns a ggplot and errors clearly without network_centrality()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_centrality(proj, condition = "FLO-ET", save = FALSE), "network_centrality")

  proj <- network_centrality(proj)
  p <- plot_centrality(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_centrality() errors clearly on a measure not present in network_centrality()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, measures = "degree")
  expect_error(plot_centrality(proj, condition = "FLO-ET", measure = "betweenness", save = FALSE), "betweenness")
})

test_that("plot_centrality() saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  p <- plot_centrality(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "centrality_plot_log")
  expect_true(file.exists(log_df$path[1]))
})
