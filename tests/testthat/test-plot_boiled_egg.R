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

test_that(".boiled_egg_limits() keeps the published frame and widens it to every point", {
  lim <- patliR:::.boiled_egg_limits(c(0, 50), c(1, 2))
  expect_equal(lim$x, c(-20, 220))
  expect_equal(lim$y, c(-3, 8))
  lim <- patliR:::.boiled_egg_limits(c(0, 300, NA), c(-5, 17.3, Inf))
  expect_true(lim$x[1] == -20 && lim$x[2] > 300)
  expect_true(lim$y[1] < -5 && lim$y[2] > 17.3)
})

test_that("plot_boiled_egg() shows points above WLogP 8 and labels by name on request", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)
  adme <- patliRResults(proj, "adme_local")
  adme$wlogp_proxy[1] <- 12
  patliRResults(proj, "adme_local") <- adme
  p <- plot_boiled_egg(proj, engine = "static", save = FALSE)
  expect_gt(p$coordinates$limits$y[2], 12)
  expect_false(is.null(p$labels$caption))

  p_named <- plot_boiled_egg(proj, engine = "static", save = TRUE, label = "name")
  expect_true(file.exists(patliRResults(attr(p_named, "proj"), "boiled_egg_log")$path_labeled))
  expect_error(plot_boiled_egg(proj, engine = "static", save = FALSE, label = "bogus"))
})
