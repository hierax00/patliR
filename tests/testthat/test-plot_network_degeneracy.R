test_that("plot_network_degeneracy() returns a ggplot and errors clearly without network_degeneracy()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_network_degeneracy(proj, condition = "FLO-ET", save = FALSE), "network_degeneracy")

  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    n_targets_a = 3L, n_targets_b = 3L, n_pathways_a = 2L, n_pathways_b = 2L,
    target_jaccard = 0.5, functional_similarity = 0.5, degeneracy_score = 0.5,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg

  p <- plot_network_degeneracy(proj, condition = "FLO-ET", min_degeneracy = 0.1, save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that(".network_degeneracy_edges() drops pairs where a compound is missing from the node table", {
  nodes <- data.frame(name = c("C1", "C2"), x = c(0, 1), y = c(0, 1), stringsAsFactors = FALSE)
  deg <- data.frame(compound_a = c("C1", "C1"), compound_b = c("C2", "C3"), degeneracy_score = c(0.5, 0.9), stringsAsFactors = FALSE)
  edges <- patliR:::.network_degeneracy_edges(nodes, deg)
  expect_equal(nrow(edges), 1)
  expect_equal(edges$degeneracy_score, 0.5)
})
