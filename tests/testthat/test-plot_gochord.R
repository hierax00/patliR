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

test_that("plot_gochord() rejects explicit empty and multiple conditions", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GOplot")
  proj <- .network_stats_test_setup()
  conditions <- unique(patliRResults(proj, "network_edges")$condition)
  expect_error(plot_gochord(proj, condition = conditions, save = FALSE), "single")
  expect_error(plot_gochord(proj, condition = character(), save = FALSE), "single")
})

test_that("plot_gochord() saves a separate file for each enrichment database", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GOplot")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "t1", weight = 1)
  patliRResults(proj, "network_enrichment") <- data.frame(condition = "A", db = c("go", "reactome"), Description = "term", p.adjust = 0.01, geneID = "1/2")
  testthat::local_mocked_bindings(GOChord = function(...) ggplot2::ggplot(), .package = "GOplot")
  testthat::local_mocked_bindings(.network_gochord_gene_labels = function(ids) stats::setNames(ids, ids), .package = "patliR")
  for (db in c("go", "reactome")) {
    p <- plot_gochord(proj, condition = "A", db = db, width = 2, height = 2)
    proj <- attr(p, "proj")
  }
  log <- patliRResults(proj, "gochord_plot_log")
  expect_equal(nrow(log), 2L)
  expect_equal(length(unique(log$path)), 2L)
  expect_setequal(basename(log$path), c("gochord_A_go.png", "gochord_A_reactome.png"))
  expect_true(all(file.exists(log$path)))
})

test_that("plot_gochord() keeps only the top_n_genes most shared genes and wraps long term names", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("GOplot")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  entrez <- as.character(c(7157, 1956, 3569, 4790, 5594, 5595, 207, 2064, 673, 5290, 1017, 3479))
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = "FLO-ET", db = "go", ID = c("GO:1", "GO:2"),
    Description = c("a very long enriched biological process description that must be wrapped in the legend", "short term"),
    p.adjust = c(0.001, 0.01),
    geneID = c(paste(entrez, collapse = "/"), paste(entrez[1:4], collapse = "/")),
    stringsAsFactors = FALSE
  )
  seen <- NULL
  testthat::local_mocked_bindings(
    GOChord = function(data, ...) { seen <<- data; ggplot2::ggplot() },
    .package = "GOplot"
  )
  plot_gochord(proj, condition = "FLO-ET", top_n_genes = 5, save = FALSE)
  expect_equal(nrow(seen), 5)
  ## the 4 genes shared by both terms are always kept
  expect_true(all(rowSums(seen)[order(-rowSums(seen))][1:4] == 2))
  expect_true(any(grepl("\n", colnames(seen), fixed = TRUE)))

  plot_gochord(proj, condition = "FLO-ET", top_n_genes = Inf, save = FALSE)
  expect_equal(nrow(seen), length(entrez))
})
