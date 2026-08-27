test_that("plot_structure2d() returns a patchwork object and errors clearly without prep_structure2d()", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("png")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(plot_structure2d(proj, save = FALSE), "structure2d_log")

  proj <- prep_structure2d(proj)
  p <- plot_structure2d(proj, save = FALSE)
  expect_s3_class(p, "patchwork")
})

test_that("plot_structure2d() saves a grid PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("png")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_structure2d(proj)
  p <- plot_structure2d(proj, save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "structure2d_grid_log")
  expect_true(file.exists(log_df$path[1]))
})
