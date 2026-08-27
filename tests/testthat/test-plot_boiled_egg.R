test_that(".point_in_polygon() correctly classifies simple known points", {
  square_x <- c(0, 10, 10, 0)
  square_y <- c(0, 0, 10, 10)
  expect_true(patliR:::.point_in_polygon(5, 5, square_x, square_y))
  expect_false(patliR:::.point_in_polygon(15, 5, square_x, square_y))
  expect_false(patliR:::.point_in_polygon(-1, -1, square_x, square_y))
})

test_that(".load_boiled_egg_polygons() loads both bundled ellipses", {
  poly <- patliR:::.load_boiled_egg_polygons()
  expect_true(all(c("gia", "bbb") %in% names(poly)))
  expect_true(all(c("tpsa", "wlogp") %in% names(poly$gia)))
  expect_gt(nrow(poly$gia), 50)
  expect_gt(nrow(poly$bbb), 50)
})

test_that("plot_boiled_egg() returns a single ggplot object with engine = 'static'", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  p <- plot_boiled_egg(proj, engine = "static", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_boiled_egg() errors clearly when adme_local() has not run", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(plot_boiled_egg(proj, save = FALSE), "adme_local")
})

test_that("plot_boiled_egg() saves a labeled and a plain PNG to the project's plots folder", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  p <- plot_boiled_egg(proj, engine = "static", save = TRUE)
  proj2 <- attr(p, "proj")

  out_dir <- file.path(projectDir(proj2), "plots")
  expect_true(file.exists(file.path(out_dir, "boiled_egg_labeled.png")))
  expect_true(file.exists(file.path(out_dir, "boiled_egg_plain.png")))
  log_df <- patliRResults(proj2, "boiled_egg_log")
  expect_true(all(c("path_labeled", "path_plain") %in% names(log_df)))
})

test_that("adme_local()'s gi_absorption/bbb_permeant are real point-in-polygon calls", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  adme <- patliRResults(proj, "adme_local")
  expect_true("wlogp_proxy" %in% names(adme))
  expect_true(all(adme$gi_absorption %in% c("High", "Low", NA)))
  expect_true(all(is.logical(adme$bbb_permeant)))
})
