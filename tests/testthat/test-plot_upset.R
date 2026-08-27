test_that("plot_upset() returns a patchwork object and needs >= 2 conditions", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  proj <- .network_stats_test_setup()
  p <- plot_upset(proj, save = FALSE) # fixture builds every condition -> pooled default
  expect_s3_class(p, "patchwork")

  expect_error(plot_upset(proj, condition = "FLO-ET", save = FALSE), "at least 2")
})

test_that("plot_upset() saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")

  proj <- .network_stats_test_setup()
  p <- plot_upset(proj, save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "upset_plot_log")
  expect_true(file.exists(log_df$path[1]))
})
