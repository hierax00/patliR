test_that("plot_adme_upset() counts every all-failing compound in the empty intersection", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  proj <- .test_project()
  dat <- expand.grid(compound_id = c("c1", "c2", "c3"), rule = c("ro5", "veber"),
                     stringsAsFactors = FALSE)
  dat$pass <- FALSE
  patliRResults(proj, "adme_filtered") <- dat

  p <- plot_adme_upset(proj, save = FALSE)
  expect_s3_class(p, "patchwork")
  expect_equal(as.character(p[[1]]$data$combo), "(none)")
  expect_equal(p[[1]]$data$n_compounds, 3L)
  expect_false(any(p[[2]]$data$in_set))
  expect_no_error(ggplot2::ggplot_build(p[[1]]))
})
