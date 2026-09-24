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

test_that("plot_upset() preserves structured membership when condition names contain pipes", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = c("A|B", "A", "B", "A|B"), compound_id = "c1",
    uniprot_id = c("t1", "t2", "t2", "t3"), weight = 1
  )
  p <- plot_upset(proj, save = FALSE)
  bars <- p[[1]]$data
  dots <- p[[2]]$data
  expect_equal(nrow(bars), 2L)
  expect_equal(sort(bars$n_targets), c(1L, 2L))
  largest <- as.character(bars$combo[which.max(bars$n_targets)])
  expect_equal(dots$condition[dots$combo == largest & dots$in_set], "A|B")
  expect_setequal(dots$condition[dots$combo != largest & dots$in_set], c("A", "B"))
})
