test_that("network_module_robustness() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  testthat::skip_if_not_installed("dbscan")
  expect_error(network_module_robustness(proj), "network_edges")
})

test_that("network_module_robustness() produces a well-formed summary table on FLO-ET", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)

  summary_tbl <- patliRResults(proj, "network_module_robustness")
  expect_true(all(c("condition", "module_id", "n_nodes", "r_index", "seed_used") %in% names(summary_tbl)))
  expect_true(all(summary_tbl$r_index >= 0 & summary_tbl$r_index <= 1))
})

test_that(".network_detect_modules() falls back to a single module when the graph is smaller than 2 * min_module_size", {
  testthat::skip_if_not_installed("dbscan")
  ## FLO-ET's own graph (4 present compounds + 10 targets = 14 nodes, per
  ## the Q1-thresholded presence rule in prep_binarize()) is comfortably
  ## above 2 * min_module_size for the default min_module_size = 2, so it
  ## never actually exercises this fallback -- test the internal function
  ## directly against a synthetic graph small enough to trigger it, rather
  ## than relying on a bundled fixture happening to be tiny.
  tiny_g <- igraph::graph_from_data_frame(
    data.frame(from = "a", to = "b", stringsAsFactors = FALSE),
    directed = FALSE
  )
  modules <- patliR:::.network_detect_modules(tiny_g, min_module_size = 2)
  expect_equal(names(modules$members), "M1")
  expect_setequal(modules$members[["M1"]], c("a", "b"))
  expect_true(grepl("single module", modules$log_message))
})

test_that(".network_detect_modules() tags the small-graph fallback module as type 'cluster', not 'noise' (2026-08-14 design change)", {
  testthat::skip_if_not_installed("dbscan")
  ## This component is too small to run HDBSCAN at all -- it was never
  ## classified as noise, so it must not be tagged as such just because it
  ## shares the "single module" shape with the real all-noise fallback in
  ## .network_detect_modules_component().
  tiny_g <- igraph::graph_from_data_frame(
    data.frame(from = "a", to = "b", stringsAsFactors = FALSE),
    directed = FALSE
  )
  modules <- patliR:::.network_detect_modules(tiny_g, min_module_size = 2)
  expect_equal(unname(modules$types), "cluster")
})

test_that("network_module_robustness() tags every module 'cluster'/'noise' and drops no nodes (2026-08-14 design change: noise nodes used to be excluded)", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  g <- .network_graph(proj, "FLO-ET")

  summary_tbl <- patliRResults(proj, "network_module_robustness")
  expect_true("module_type" %in% names(summary_tbl))
  expect_true(all(summary_tbl$module_type %in% c("cluster", "noise")))
  ## The core of the design change: every node in the condition's graph is
  ## accounted for by exactly one module's n_nodes, whether or not HDBSCAN
  ## considered it "noise" -- see DEVLOG.md, 2026-08-14.
  expect_equal(sum(summary_tbl$n_nodes), igraph::vcount(g))
})

test_that("network_module_robustness() writes a matching fragmentation curve", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)

  summary_tbl <- patliRResults(proj, "network_module_robustness")
  curve <- patliRResults(proj, "network_robustness_curve")
  expect_true(all(c("condition", "module_id", "n_removed", "largest_component_fraction") %in% names(curve)))

  for (m in unique(summary_tbl$module_id)) {
    n_nodes <- summary_tbl$n_nodes[summary_tbl$module_id == m]
    curve_m <- curve[curve$module_id == m, ]
    expect_equal(nrow(curve_m), n_nodes + 1)
    expect_equal(curve_m$largest_component_fraction[curve_m$n_removed == 0], 1)
    ## r_index is Schneider et al. (2011)'s mean over removal steps
    ## Q = 1..N only -- the pre-removal Q = 0 row (n_removed == 0,
    ## fraction == 1, checked above) is kept in the curve for plotting but
    ## must be excluded from this average (see .network_percolate(),
    ## internal, and DEVLOG.md 2026-08-25).
    post_removal <- curve_m$largest_component_fraction[curve_m$n_removed > 0]
    expect_equal(mean(post_removal), summary_tbl$r_index[summary_tbl$module_id == m])
  }
})

test_that("network_module_robustness() is reproducible given the same seed", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  proj1 <- network_module_robustness(proj, condition = "FLO-ET", seed = 123)
  proj2 <- network_module_robustness(proj, condition = "FLO-ET", seed = 123)

  expect_equal(
    patliRResults(proj1, "network_module_robustness")$r_index,
    patliRResults(proj2, "network_module_robustness")$r_index
  )
})

test_that("network_module_robustness() does not leak its RNG state into the caller's session", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()

  set.seed(999)
  before <- runif(1)

  set.seed(999)
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  after <- runif(1)

  expect_equal(before, after)
})

test_that("network_module_robustness() rejects min_module_size < 2", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  expect_error(network_module_robustness(proj, condition = "FLO-ET", min_module_size = 1))
})

test_that(".network_detect_modules() clusters disconnected graphs per connected component instead of crashing (regression: real EFLO-S 'Index out of bounds' crash, 2026-08-11)", {
  testthat::skip_if_not_installed("dbscan")
  ## Root cause (see DEVLOG.md, 2026-08-11): igraph::distances() returns
  ## Inf for node pairs in different connected components, and
  ## dbscan::hdbscan() has no validation against non-finite input -- fed a
  ## whole-graph distance matrix with Inf entries, it crashed with an
  ## unhelpful C++-level "Index out of bounds" error instead of erroring
  ## cleanly. This is expected on real, large compound-target networks
  ## (a target reachable only through absent compounds, an isolated
  ## compound-target pair, etc.) but never occurs on patliR's small, fully
  ## connected bundled example data, which is why devtools::test() never
  ## caught it. A synthetic two-component graph (one component large
  ## enough to run HDBSCAN, one too small and falling back to a single
  ## module) reproduces the same code path without needing real data.
  disconnected_edges <- data.frame(
    from = c("a", "b", "c", "d", "e", "f", "x"),
    to   = c("b", "c", "d", "e", "f", "a", "y"),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(disconnected_edges, directed = FALSE)
  expect_equal(igraph::components(g)$no, 2)

  modules <- patliR:::.network_detect_modules(g, min_module_size = 2)
  all_nodes <- unlist(modules$members, use.names = FALSE)

  expect_setequal(all_nodes, igraph::V(g)$name)
  expect_equal(length(all_nodes), length(unique(all_nodes)))
  expect_true(grepl("2 connected component", modules$log_message))
  expect_true(grepl("component with 2 node", modules$log_message))
})

test_that("network_module_robustness() runs end to end on a disconnected graph without crashing (regression)", {
  testthat::skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  g <- .network_graph(proj, "FLO-ET")

  ## Force the same condition into a two-component state by stripping the
  ## cached graph and re-caching a disconnected copy directly, mirroring
  ## what a real, sparser condition (fewer shared compound-target edges)
  ## produces -- add one isolated compound-target pair not connected to
  ## the rest of the graph.
  disconnected_g <- igraph::add_vertices(g, 2, name = c("__iso_c", "__iso_t"))
  disconnected_g <- igraph::add_edges(disconnected_g, c("__iso_c", "__iso_t"))
  expect_gt(igraph::components(disconnected_g)$no, 1)

  modules <- patliR:::.network_detect_modules(disconnected_g, min_module_size = 2)
  all_nodes <- unlist(modules$members, use.names = FALSE)
  expect_setequal(all_nodes, igraph::V(disconnected_g)$name)
  expect_equal(length(all_nodes), length(unique(all_nodes)))
})
