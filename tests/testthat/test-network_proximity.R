test_that("network_proximity() requires network_build() to have run first", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_proximity(proj, disease = "EFO_0000000"), "network_edges")
})

test_that("network_proximity() requires targets_disease_filter() to have run", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  expect_error(network_proximity(proj, condition = "FLO-ET", disease = "EFO_0000000"), "targets_disease")
})

test_that("network_proximity() errors on a condition network_build() never built", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  expect_error(
    network_proximity(proj, condition = "NOT-A-REAL-CONDITION", disease = "EFO_0000000"),
    "not built"
  )
})

test_that("network_proximity() errors clearly on an unknown disease_id", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  ## Fabricate a minimal targets_disease table directly (bypassing a live
  ## Open Targets call) so this test is self-contained and offline.
  ct <- unique(patliRResults(proj, "network_edges")[, c("compound_id", "uniprot_id")])
  fake_disease <- data.frame(
    compound_id = ct$compound_id[1], target_id = ct$uniprot_id[1],
    disease_id = "SOME_REAL_DISEASE", association_score = 0.5,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease") <- fake_disease

  expect_error(
    network_proximity(proj, condition = "FLO-ET", disease = "NOT-A-REAL-DISEASE"),
    "targets_disease"
  )
})

test_that(".network_degree_bins()/.network_resample_degree_matched() preserve length and vocabulary", {
  testthat::skip_if_not_installed("STRINGdb")
  set.seed(1)
  degree_all <- stats::setNames(sample(1:50, 30, replace = TRUE), paste0("n", 1:30))
  bins <- patliR:::.network_degree_bins(degree_all)
  resampled <- patliR:::.network_resample_degree_matched(names(degree_all)[1:5], degree_all, bins)
  expect_length(resampled, 5)
  expect_true(all(resampled %in% names(degree_all)))
})

test_that(".network_closest_distance() tolerates duplicate source/target ids", {
  ## Real bug, found by Uriel running network_proximity() against a live
  ## STRINGdb download (2026-07-23): igraph::distances()'s underlying C
  ## routine errors outright ("Target vertex list must not have any
  ## duplicates") if `to=` has a repeated vertex -- and
  ## .network_resample_degree_matched() legitimately produces duplicates
  ## (independent sampling *with replacement* per node). No STRINGdb/
  ## internet needed for this one -- a tiny synthetic igraph reproduces it.
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"), stringsAsFactors = FALSE),
    directed = FALSE
  )
  expect_no_error(
    d <- patliR:::.network_closest_distance(g, c("a", "a", "b"), c("d", "d"))
  )
  expect_true(is.numeric(d) && length(d) == 1)
})

test_that("network_proximity() end-to-end is not exercised automatically -- needs a real STRING download", {
  testthat::skip_if_not_installed("STRINGdb")
  skip_on_cran()
  testthat::skip(
    "network_proximity() end-to-end needs a real STRINGdb flat-file download (tens-hundreds of MB) -- verify manually against real data, not in the automated suite."
  )
})
