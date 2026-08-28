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

test_that(".network_degree_bins() builds >=min_per_bin bins from consecutive degrees (Guney rule)", {
  set.seed(1)
  degree_all <- stats::setNames(
    pmax(1L, rpois(2000, 4) + rbinom(2000, 40, 0.1)),
    paste0("n", 1:2000)
  )
  bins <- patliR:::.network_degree_bins(degree_all, min_per_bin = 100)

  sizes <- lengths(bins)
  expect_true(all(utils::head(sizes, -1) >= 100))          # all but last >= 100
  deg_ranges <- lapply(bins, function(ix) range(degree_all[ix]))
  for (i in seq_len(length(bins) - 1)) {                   # contiguous in degree
    expect_lte(deg_ranges[[i]][2], deg_ranges[[i + 1]][1])
  }
})

test_that(".network_resample_degree_matched() returns an equal-size set of distinct, degree-matched nodes", {
  set.seed(2)
  degree_all <- stats::setNames(sample(1:60, 400, replace = TRUE), paste0("n", 1:400))
  bins <- patliR:::.network_degree_bins(degree_all, min_per_bin = 50)

  input <- names(degree_all)[1:12]
  out <- patliR:::.network_resample_degree_matched(input, degree_all, bins)
  expect_length(out, length(input))
  expect_false(anyDuplicated(out) > 0)
  expect_true(all(out %in% names(degree_all)))

  bin_of <- integer(length(degree_all))
  for (b in seq_along(bins)) bin_of[bins[[b]]] <- b
  names(bin_of) <- names(degree_all)
  expect_equal(unname(bin_of[out]), unname(bin_of[input]))
})

test_that(".network_closest_distance() tolerates duplicate source/target ids", {
  ## igraph::distances()'s underlying C routine errors outright ("Target
  ## vertex list must not have any duplicates") if `to=` has a repeated
  ## vertex. The observed source/target sets can still carry a repeat, so
  ## .network_closest_distance() dedups defensively. A tiny synthetic
  ## igraph reproduces the constraint, no STRINGdb needed.
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
