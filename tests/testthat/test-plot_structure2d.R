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

test_that(".structure2d_height() grows with the number of grid rows, never below the minimum", {
  expect_equal(patliR:::.structure2d_height(4, 4, 8), 8)
  expect_equal(patliR:::.structure2d_height(108, 4, 8), 27 * 2.3)
  expect_equal(patliR:::.structure2d_height(9, 1, 8), 9 * 2.3)
})

test_that("plot_structure2d() saves a taller figure when height is not given", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("png")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_structure2d(proj)
  n <- sum(patliRResults(proj, "structure2d_log")$generated_ok)
  testthat::skip_if(n < 4, "need at least 4 depictions")
  p <- plot_structure2d(proj, ncol = 1, dpi = 20, save = TRUE)
  img <- png::readPNG(patliRResults(attr(p, "proj"), "structure2d_grid_log")$path)
  expect_equal(dim(img)[1], round(patliR:::.structure2d_height(n, 1, 8) * 20), tolerance = 2)
  p2 <- plot_structure2d(proj, ncol = 1, height = 5, dpi = 20, save = TRUE)
  img2 <- png::readPNG(patliRResults(attr(p2, "proj"), "structure2d_grid_log")$path)
  expect_equal(dim(img2)[1], 100)
})
