test_that("plot_gochord() returns a ggplot and errors clearly without network_enrich()", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GOplot")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  expect_error(plot_gochord(proj, condition = "FLO-ET", save = FALSE), "network_enrichment")

  proj <- tryCatch(network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE), error = function(e) NULL)
  testthat::skip_if(is.null(proj), "network_enrich() could not run in this environment")
  enr <- patliRResults(proj, "network_enrichment")
  testthat::skip_if(is.null(enr) || nrow(enr) == 0, "no significant terms in this fixture")

  p <- plot_gochord(proj, condition = "FLO-ET", top_n_terms = 5, save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_gochord() rejects pooling across conditions", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GOplot")

  proj <- .network_stats_test_setup()
  expect_error(plot_gochord(proj, save = FALSE), "single")
})
