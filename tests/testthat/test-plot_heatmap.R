test_that("plot_heatmap() returns a pheatmap object for both 'what' modes", {
  testthat::skip_if_not_installed("pheatmap")

  proj <- .network_stats_test_setup()
  p1 <- plot_heatmap(proj, what = "compound_condition", save = FALSE)
  expect_s3_class(p1, "pheatmap")

  p2 <- plot_heatmap(proj, condition = "FLO-ET", what = "compound_target", save = FALSE)
  expect_s3_class(p2, "pheatmap")
})

test_that("plot_heatmap() saves a PDF and logs it", {
  testthat::skip_if_not_installed("pheatmap")

  proj <- .network_stats_test_setup()
  p <- plot_heatmap(proj, what = "compound_condition", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "heatmap_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.pdf$")
})
