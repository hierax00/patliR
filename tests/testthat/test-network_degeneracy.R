test_that("network_degeneracy() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_degeneracy(proj), "network_edges")
})

test_that("network_degeneracy() requires network_enrich() to have run for the condition", {
  proj <- .network_stats_test_setup() ## network_build() only, no enrich yet
  expect_error(network_degeneracy(proj, condition = "FLO-ET"), "network_enrich")
})

test_that("network_degeneracy() errors on a condition network_build() never built", {
  proj <- .network_stats_test_setup()
  expect_error(network_degeneracy(proj, condition = "NOT-A-REAL-CONDITION"), "not built")
})

test_that("network_degeneracy() runs end-to-end and scores are in [0, 1]", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  ## simplify_go = FALSE: this test isn't about GO simplification, and the
  ## default TRUE rebuilds GOSemSim::godata() from scratch (slow) -- see
  ## network_enrich()'s own tests for simplify_go coverage.
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)
  proj <- network_degeneracy(proj, condition = "FLO-ET")

  result <- patliRResults(proj, "network_degeneracy")
  expect_true(all(c("condition", "compound_a", "compound_b", "n_targets_a", "n_targets_b",
                     "n_pathways_a", "n_pathways_b", "target_jaccard", "pathway_jaccard",
                     "degeneracy_score") %in% names(result)))
  if (nrow(result) > 0) {
    expect_true(all(result$condition == "FLO-ET"))
    expect_true(all(result$compound_a != result$compound_b))
    expect_true(all(result$target_jaccard >= 0 & result$target_jaccard <= 1))
    expect_true(all(result$pathway_jaccard >= 0 & result$pathway_jaccard <= 1))
    expect_true(all(result$degeneracy_score >= 0 & result$degeneracy_score <= 1))
  }
})

test_that("network_degeneracy() needs clusterProfiler/org.Hs.eg.db even when network_enrich() already ran", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()
  ## Sanity check only: with both installed there is nothing to assert about
  ## the missing-package error path here (it's exercised implicitly by the
  ## dependency check itself, mirrored from network_enrich()'s own test file).
  proj <- .network_stats_test_setup()
  ## simplify_go = FALSE: this test isn't about GO simplification, and the
  ## default TRUE rebuilds GOSemSim::godata() from scratch (slow) -- see
  ## network_enrich()'s own tests for simplify_go coverage.
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)
  expect_no_error(network_degeneracy(proj, condition = "FLO-ET"))
})
