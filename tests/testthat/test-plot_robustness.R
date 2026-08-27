test_that("plot_robustness() returns a ggplot and errors clearly without network_module_robustness()", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dbscan")

  proj <- .network_stats_test_setup()
  expect_error(plot_robustness(proj, condition = "FLO-ET", save = FALSE), "network_module_robustness")

  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 1)
  p <- plot_robustness(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_robustness() saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("dbscan")

  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 1)
  p <- plot_robustness(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "robustness_plot_log")
  expect_true(file.exists(log_df$path[1]))
})
