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

test_that(".network_degeneracy_edges() carries z_score through when present, NA_real_ when absent", {
  nodes <- data.frame(name = c("C1", "C2"), x = c(0, 1), y = c(0, 1), stringsAsFactors = FALSE)
  deg <- data.frame(compound_a = "C1", compound_b = "C2", degeneracy_score = 0.5, z_score = 2.1, stringsAsFactors = FALSE)
  edges <- patliR:::.network_degeneracy_edges(nodes, deg)
  expect_equal(edges$z_score, 2.1)

  deg_no_z <- data.frame(compound_a = "C1", compound_b = "C2", degeneracy_score = 0.5, stringsAsFactors = FALSE)
  edges_no_z <- patliR:::.network_degeneracy_edges(nodes, deg_no_z)
  expect_true(is.na(edges_no_z$z_score))
})

test_that("plot_network_degeneracy() default filter = 'p_adjusted' falls back to 'score' and informs when p_adjusted is all-NA", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  ## annotation = "jaccard"-shaped row: p_adjusted/z_score all NA, as piece 14 leaves them
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    n_targets_a = 3L, n_targets_b = 3L, n_pathways_a = 2L, n_pathways_b = 2L,
    target_jaccard = 0.5, functional_similarity = 0.5,
    z_score = NA_real_, p_adjusted = NA_real_, degeneracy_score = 0.5,
    annotation = "jaccard", stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg

  expect_message(
    p <- plot_network_degeneracy(proj, condition = "FLO-ET", min_degeneracy = 0.1, save = FALSE),
    "falling back"
  )
  expect_s3_class(p, "ggplot")
})

test_that("plot_network_degeneracy() filter = 'score' reproduces the historical min_degeneracy cutoff exactly (regression)", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    n_targets_a = 3L, n_targets_b = 3L, n_pathways_a = 2L, n_pathways_b = 2L,
    target_jaccard = 0.5, functional_similarity = 0.5,
    z_score = 2.5, p_adjusted = 0.5, degeneracy_score = 0.5,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg

  ## min_degeneracy = 0.6 excludes the one pair (score 0.5) -> zero links, still a valid plot
  p_excluded <- plot_network_degeneracy(proj, condition = "FLO-ET", filter = "score", min_degeneracy = 0.6, save = FALSE)
  expect_s3_class(p_excluded, "ggplot")
  ## min_degeneracy = 0.1 includes it -- same cutoff arithmetic as the pre-revision code
  p_included <- plot_network_degeneracy(proj, condition = "FLO-ET", filter = "score", min_degeneracy = 0.1, save = FALSE)
  expect_s3_class(p_included, "ggplot")
})

test_that("plot_network_degeneracy() filter = 'p_adjusted' selects only significant pairs", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    n_targets_a = 3L, n_targets_b = 3L, n_pathways_a = 2L, n_pathways_b = 2L,
    target_jaccard = 0.5, functional_similarity = 0.9,
    z_score = 3.2, p_adjusted = 0.01, degeneracy_score = 0.9,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg

  p <- plot_network_degeneracy(proj, condition = "FLO-ET", filter = "p_adjusted", alpha = 0.05, save = FALSE)
  expect_s3_class(p, "ggplot")

  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  g <- patliR:::.network_layers_filter(proj, g, "FLO-ET", layers = c("compound", "target"), max_pathways = 0)
  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 15, seed = 1)
  deg_scope <- fake_deg[fake_deg$p_adjusted < 0.05, , drop = FALSE]
  edges <- patliR:::.network_degeneracy_edges(pd$nodes, deg_scope)
  expect_equal(nrow(edges), 1)
  expect_equal(edges$z_score, 3.2)
})

test_that("plot_network_degeneracy() inherits colour_by = 'module' from plot_network_layers()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    z_score = 3.2, p_adjusted = 0.01, degeneracy_score = 0.9,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)

  p <- plot_network_degeneracy(proj, condition = "FLO-ET", colour_by = "module", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_network_degeneracy() supports layout = 'bipartite' (its scope is always two-layer)", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_deg <- data.frame(
    condition = "FLO-ET", compound_a = ct[1], compound_b = ct[2],
    z_score = 3.2, p_adjusted = 0.01, degeneracy_score = 0.9,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- fake_deg
  p <- plot_network_degeneracy(proj, condition = "FLO-ET", layout = "bipartite", save = FALSE)
  expect_s3_class(p, "ggplot")
})
