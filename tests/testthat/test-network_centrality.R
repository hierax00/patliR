test_that("network_centrality() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_centrality(proj), "network_edges")
})

test_that("network_centrality() errors on a condition network_build() never built", {
  proj <- .network_stats_test_setup()
  expect_error(network_centrality(proj, condition = "NOT-A-REAL-CONDITION"), "not built")
})

test_that("network_centrality() computes degree/betweenness/hub_score for every node", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, condition = "FLO-ET")

  result <- patliRResults(proj, "network_centrality")
  expect_true(all(c("condition", "node_id", "node_type", "degree", "betweenness", "hub_score") %in% names(result)))
  expect_setequal(unique(result$node_type), c("compound", "target"))

  g <- patliR:::.network_graph(proj, "FLO-ET")
  expect_equal(nrow(result), igraph::vcount(g))
})

test_that("network_centrality() respects the measures argument", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, condition = "FLO-ET", measures = "degree")

  result <- patliRResults(proj, "network_centrality")
  expect_true("degree" %in% names(result))
  expect_false("betweenness" %in% names(result))
  expect_false("hub_score" %in% names(result))
})

test_that("network_centrality() rebuilding one condition does not touch the others", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  before <- patliRResults(proj, "network_centrality")

  proj <- network_centrality(proj, condition = "FLO-ET")
  after <- patliRResults(proj, "network_centrality")

  expect_setequal(unique(after$condition), unique(before$condition))
})
