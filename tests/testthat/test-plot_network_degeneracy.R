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

  ## min_degeneracy = 0.6 excludes the one pair (score 0.5) -> zero links,
  ## still a valid plot, but no longer a silent one
  expect_warning(
    p_excluded <- plot_network_degeneracy(proj, condition = "FLO-ET", filter = "score", min_degeneracy = 0.6, save = FALSE),
    "no degeneracy link"
  )
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

.degeneracy_fixture <- function(proj) {
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 3, "need at least 3 compounds in this fixture's FLO-ET condition")
  pairs <- utils::combn(ct[1:3], 2)
  patliRResults(proj, "network_degeneracy") <- data.frame(
    condition = "FLO-ET", compound_a = pairs[1, ], compound_b = pairs[2, ],
    z_score = c(4, 1, -0.5), p_adjusted = c(0.01, 0.3, 0.9), degeneracy_score = c(0.8, 0.4, 0.1),
    stringsAsFactors = FALSE
  )
  proj
}

test_that("plot_network_degeneracy() aborts, naming the available conditions, when the requested one has no rows", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .degeneracy_fixture(.network_stats_test_setup())
  deg <- patliRResults(proj, "network_degeneracy")
  deg$condition <- "OTHER"
  patliRResults(proj, "network_degeneracy") <- deg
  expect_error(plot_network_degeneracy(proj, condition = "FLO-ET", save = FALSE), "OTHER")
})

test_that("plot_network_degeneracy() warns when no pair passes the filter", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .degeneracy_fixture(.network_stats_test_setup())
  expect_warning(
    p <- plot_network_degeneracy(proj, condition = "FLO-ET", alpha = 0.001, save = FALSE),
    "no degeneracy link"
  )
  expect_s3_class(p, "ggplot")
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomCurve"), logical(1))))
})

test_that("degeneracy links are drawn only for passing pairs, coloured by z and sized by score", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .degeneracy_fixture(.network_stats_test_setup())
  p <- plot_network_degeneracy(proj, condition = "FLO-ET", save = FALSE)
  curves <- Filter(function(l) inherits(l$geom, "GeomCurve"), p$layers)[[1]]
  expect_equal(nrow(curves$data), 1)
  expect_equal(curves$data$z_score, 4)
  expect_true(all(c("colour", "linewidth") %in% names(curves$mapping)))
})

test_that("layout = 'circle' draws only the compounds, on the unit circle, with every link", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .degeneracy_fixture(.network_stats_test_setup())
  p <- plot_network_degeneracy(proj, condition = "FLO-ET", layout = "circle", filter = "score",
                               min_degeneracy = 0.3, save = FALSE)
  pts <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_true(all(pts$layer == "compound"))
  expect_equal(sqrt(pts$x^2 + pts$y^2), rep(1, nrow(pts)))
  curves <- Filter(function(l) inherits(l$geom, "GeomCurve"), p$layers)[[1]]
  expect_equal(nrow(curves$data), 2)
})

test_that("view = 'heatmap' tiles every scored pair symmetrically and marks the passing ones", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .degeneracy_fixture(.network_stats_test_setup())
  p <- plot_network_degeneracy(proj, condition = "FLO-ET", view = "heatmap", save = TRUE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomTile"), logical(1))))
  expect_equal(nrow(p$data), 6)
  stars <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]$data
  expect_equal(nrow(stars), 2)
  expect_true(all(stars$value == 4))
  log_df <- patliRResults(attr(p, "proj"), "network_degeneracy_plot_log")
  expect_identical(log_df$view, "heatmap")
  expect_true(file.exists(log_df$path))
})

test_that(".network_degeneracy_order() puts the most degenerate pair next to each other", {
  deg <- data.frame(condition = "A", compound_a = c("a", "a", "b", "a"), compound_b = c("d", "b", "c", "c"),
                    z_score = c(9, 0, 0, 0), degeneracy_score = 0.1, stringsAsFactors = FALSE)
  ord <- patliR:::.network_degeneracy_order(deg, c("a", "b", "c", "d"))
  expect_setequal(ord, c("a", "b", "c", "d"))
  expect_equal(abs(diff(match(c("a", "d"), ord))), 1)
})
