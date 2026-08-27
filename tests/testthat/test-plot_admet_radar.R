test_that("plot_admet_radar() returns one named ggplot per compound with engine = 'static'", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  radars <- plot_admet_radar(proj, engine = "static", save = FALSE)
  expect_type(radars, "list")
  expect_equal(length(radars), nrow(compounds(proj)))
  expect_setequal(names(radars), compounds(proj)$id)
  expect_s3_class(radars[[1]], "ggplot")
})

test_that("plot_admet_radar() errors clearly when adme_local() has not run", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(plot_admet_radar(proj, save = FALSE), "adme_local")
})

test_that("plot_admet_radar() can subset to specific compound_ids", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  one_id <- compounds(proj)$id[1]
  radars <- plot_admet_radar(proj, compound_ids = one_id, engine = "static", save = FALSE)
  expect_equal(length(radars), 1)
  expect_equal(names(radars), one_id)
  expect_s3_class(radars[[one_id]], "ggplot")
})

test_that("plot_admet_radar() saves a labeled and a plain PNG per compound", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  one_id <- compounds(proj)$id[1]
  radars <- plot_admet_radar(proj, compound_ids = one_id, engine = "static", save = TRUE)
  proj2 <- attr(radars, "proj")

  out_dir <- file.path(projectDir(proj2), "plots", "admet_radar")
  expect_true(file.exists(file.path(out_dir, paste0(one_id, "_labeled.png"))))
  expect_true(file.exists(file.path(out_dir, paste0(one_id, "_plain.png"))))

  log_df <- patliRResults(proj2, "admet_radar_log")
  expect_true(all(c("compound_id", "path_labeled", "path_plain") %in% names(log_df)))
})

test_that("plot_admet_radar()'s subtitle cites the SwissADME radar source, not just 'published'", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  one_id <- compounds(proj)$id[1]
  radars <- plot_admet_radar(proj, compound_ids = one_id, engine = "static", save = FALSE)
  subtitle <- radars[[one_id]]$labels$subtitle
  expect_true(grepl("Daina", subtitle))
  expect_true(grepl("2017", subtitle))
})

test_that(".smiles_approx_descriptors() gives sane fractions for a simple aromatic SMILES", {
  ## Apigenin: mix of aromatic rings and one aliphatic carbonyl carbon.
  d <- patliR:::.smiles_approx_descriptors(
    "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12", 20L
  )
  expect_true(d$fraction_csp3_approx >= 0 && d$fraction_csp3_approx <= 1)
  expect_true(d$aromatic_proportion_approx >= 0 && d$aromatic_proportion_approx <= 1)
})
