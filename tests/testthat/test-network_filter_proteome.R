test_that("network_filter_proteome() keeps only edges whose target is in the given proteome", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  keep_target <- edges$uniprot_id[1]

  proj <- network_filter_proteome(proj, proteome = keep_target)
  filtered <- patliRResults(proj, "network_filtered_edges")

  expect_true(all(c("condition", "compound_id", "uniprot_id", "weight") %in% names(filtered)))
  expect_true(all(filtered$uniprot_id == keep_target))
  expect_true(nrow(filtered) == sum(edges$uniprot_id == keep_target))
})

test_that("network_filter_proteome() returns zero rows (not an error) when nothing matches", {
  proj <- .network_stats_test_setup()
  proj <- network_filter_proteome(proj, proteome = "NOT-A-REAL-UNIPROT-ID")
  filtered <- patliRResults(proj, "network_filtered_edges")
  expect_equal(nrow(filtered), 0)
})

test_that("network_filter_proteome() never modifies network_edges itself", {
  proj <- .network_stats_test_setup()
  before <- patliRResults(proj, "network_edges")
  proj <- network_filter_proteome(proj, proteome = before$uniprot_id[1])
  expect_identical(patliRResults(proj, "network_edges"), before)
})
