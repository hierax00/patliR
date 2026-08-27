test_that("plot_network_layers() returns a single ggplot object with engine = 'static', per condition", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p <- plot_network_layers(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_network_layers() defaults to condition = NULL, pooling every built condition", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p_all <- plot_network_layers(proj, save = FALSE)
  expect_s3_class(p_all, "ggplot")

  g_all <- patliR:::.network_layered_graph_multi(proj, patliR:::.network_resolve_conditions(proj, NULL))
  g_one <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  # pooling across every condition should never have fewer nodes than any single condition
  expect_gte(igraph::vcount(g_all), igraph::vcount(g_one))
})

test_that("plot_network_layers() errors clearly on an unbuilt condition", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_network_layers(proj, condition = "not_a_real_condition", save = FALSE), "not built")
})

test_that("plot_network_layers() saves a PNG and logs it to network_layers_plot_log, keyed by scope", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p <- plot_network_layers(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")

  out_dir <- file.path(projectDir(proj2), "plots")
  expect_true(file.exists(file.path(out_dir, "network_layers_FLO-ET.png")))
  log_df <- patliRResults(proj2, "network_layers_plot_log")
  expect_true(all(c("condition", "path", "n_nodes", "n_edges") %in% names(log_df)))
  expect_true("FLO-ET" %in% log_df$condition)

  p_all <- plot_network_layers(proj2, save = TRUE)
  proj3 <- attr(p_all, "proj")
  expect_true(file.exists(file.path(out_dir, "network_layers_ALL.png")))
  log_df3 <- patliRResults(proj3, "network_layers_plot_log")
  expect_true("ALL" %in% log_df3$condition)
})

test_that(".network_layered_graph_multi() pools compound-target edges across conditions without duplicating shared nodes", {
  proj <- .network_stats_test_setup()
  edges_all <- patliRResults(proj, "network_edges")
  built <- unique(edges_all$condition)
  testthat::skip_if(length(built) < 2, "fixture only builds one condition, nothing to pool")

  g_multi <- patliR:::.network_layered_graph_multi(proj, built)
  ct_pool <- unique(edges_all[edges_all$condition %in% built, c("compound_id", "uniprot_id")])
  expect_equal(igraph::vcount(g_multi), length(unique(c(ct_pool$compound_id, ct_pool$uniprot_id))))
})

test_that("plot_network_layers() excludes pathway and disease layers by default", {
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, patliR:::.network_resolve_conditions(proj, NULL))
  testthat::skip_if_not(any(c("pathway", "disease") %in% igraph::V(g)$layer), "no pathway/disease layer in this fixture")

  g_filtered <- patliR:::.network_layers_filter(
    proj, g, patliR:::.network_resolve_conditions(proj, NULL),
    layers = eval(formals(plot_network_layers)$layers), max_pathways = eval(formals(plot_network_layers)$max_pathways)
  )
  expect_false(any(c("pathway", "disease") %in% igraph::V(g_filtered)$layer))
})

test_that(".network_layers_filter() caps pathway nodes to max_pathways when the pathway layer is requested", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  proj <- tryCatch(network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE), error = function(e) NULL)
  testthat::skip_if(is.null(proj), "network_enrich() could not run in this environment")

  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  n_pathways <- sum(igraph::V(g)$layer == "pathway")
  testthat::skip_if(n_pathways <= 5, "not enough pathway nodes in this fixture to exercise the cap")

  g_capped <- patliR:::.network_layers_filter(proj, g, "FLO-ET", layers = c("compound", "target", "pathway"), max_pathways = 5)
  expect_lte(sum(igraph::V(g_capped)$layer == "pathway"), 5)
})

test_that(".network_layers_plot_data() places every graph vertex and marks the top-degree nodes as hubs", {
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")

  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 2, seed = 1)
  expect_equal(nrow(pd$nodes), igraph::vcount(g))
  expect_equal(nrow(pd$edges), igraph::ecount(g))
  expect_true(sum(pd$nodes$is_hub) <= 2)
  expect_true(all(c("x", "y", "label", "degree") %in% names(pd$nodes)))
})

test_that(".network_layers_plot_data() is reproducible for a fixed seed", {
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")

  pd1 <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 5, seed = 42)
  pd2 <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 5, seed = 42)
  expect_equal(pd1$nodes$x, pd2$nodes$x)
  expect_equal(pd1$nodes$y, pd2$nodes$y)
})
