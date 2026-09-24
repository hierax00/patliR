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
test_that("legacy robustness summaries retain each module's own R annotation", {
  skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "t1")
  patliRResults(proj, "network_robustness_curve") <- data.frame(
    condition = "A", module_id = rep(c("m1", "m2"), each = 2),
    n_removed = rep(0:1, 2), largest_component_fraction = c(1, 0, 1, 0)
  )
  patliRResults(proj, "network_module_robustness") <- data.frame(
    condition = "A", module_id = c("m1", "m2"), r_index = c(0.1, 0.3)
  )
  p <- plot_robustness(proj, condition = "A", save = FALSE)
  ann <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]$data
  expect_equal(ann$label, c("R = 0.100", "R = 0.300"))
})
