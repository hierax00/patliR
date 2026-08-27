.network_enrich_test_setup <- function() {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )
  proj <- network_build(proj)
  proj
}

test_that("network_enrich() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  ## no network_build() -- but skip straight past the dependency check so
  ## this test is meaningful even on a machine without clusterProfiler/org.Hs.eg.db.
  ## db = "go" on purpose (not the "reactome" default): this test isn't about
  ## reactome specifically, and ReactomePA may not be installed even when
  ## clusterProfiler/org.Hs.eg.db are -- that dependency error would otherwise
  ## fire before ever reaching the network_edges check this test targets.
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  expect_error(network_enrich(proj, db = "go"), "network_edges")
})

test_that("network_enrich() errors on a condition network_build() never built", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, condition = "NOT-A-REAL-CONDITION", db = "go"), "not built")
})

test_that("network_enrich() validates pvalueCutoff/qvalueCutoff", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, pvalueCutoff = 0), "pvalueCutoff")
  expect_error(network_enrich(proj, pvalueCutoff = 1.5), "pvalueCutoff")
  expect_error(network_enrich(proj, qvalueCutoff = -0.1), "qvalueCutoff")
})

test_that("network_enrich() runs end-to-end for db = 'go' when dependencies are installed", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  ## simplify_go defaults to TRUE -- this is the one test in the suite that
  ## deliberately pays the GOSemSim::godata() build cost, so the default
  ## path itself stays covered (see test-network_degeneracy.R/
  ## test-network_motifs.R/test-plot_*.R for why they opt out with
  ## simplify_go = FALSE instead of each rebuilding it too).
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go")

  result <- patliRResults(proj, "network_enrichment")
  expect_true(all(c("condition", "db", "ID", "Description", "GeneRatio", "BgRatio",
                     "pvalue", "p.adjust", "qvalue", "geneID", "Count") %in% names(result)))
  if (nrow(result) > 0) {
    expect_true(all(result$condition == "FLO-ET"))
    expect_true(all(result$db == "go"))
    ## simplify_go ran (not skipped/failed silently) -- either a reduction
    ## message or, for this tiny example, possibly nothing to simplify.
    log_msg <- projectLog(proj)$message[projectLog(proj)$step == "network_enrich"]
    expect_false(any(grepl("network_enrich_simplify_failed", log_msg)))
  }
})

test_that("network_enrich() with simplify_go = FALSE skips GOSemSim entirely", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)

  result <- patliRResults(proj, "network_enrichment")
  expect_true(all(c("condition", "db", "ID", "Description") %in% names(result)))
  log_msg <- projectLog(proj)$message[projectLog(proj)$step == "network_enrich"]
  expect_false(any(grepl("simplify_go", log_msg)))
})

test_that("network_enrich() with db = 'reactome' requires ReactomePA specifically", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if(requireNamespace("ReactomePA", quietly = TRUE), "ReactomePA is installed -- nothing to test here")

  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, db = "reactome"), "ReactomePA")
})
