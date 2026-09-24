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

test_that("plot_network_layers() saves a PNG and logs it to network_layers_plot_log, keyed by (condition, layout, colour_by)", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p <- plot_network_layers(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")

  out_dir <- file.path(projectDir(proj2), "plots")
  expect_true(file.exists(file.path(out_dir, "network_layers_FLO-ET_fr_layer.png")))
  log_df <- patliRResults(proj2, "network_layers_plot_log")
  expect_true(all(c("condition", "layout", "colour_by", "path", "n_nodes", "n_edges") %in% names(log_df)))
  expect_true("FLO-ET" %in% log_df$condition)

  p_all <- plot_network_layers(proj2, save = TRUE)
  proj3 <- attr(p_all, "proj")
  expect_true(file.exists(file.path(out_dir, "network_layers_ALL_fr_layer.png")))
  log_df3 <- patliRResults(proj3, "network_layers_plot_log")
  expect_true("ALL" %in% log_df3$condition)
})

test_that("plot_network_layers() keeps distinct PNGs/log rows per (layout, colour_by) combination", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p1 <- plot_network_layers(proj, condition = "FLO-ET", layout = "fr", colour_by = "layer", save = TRUE)
  proj2 <- attr(p1, "proj")
  p2 <- plot_network_layers(proj2, condition = "FLO-ET", layout = "bipartite", colour_by = "node_type", save = TRUE)
  proj3 <- attr(p2, "proj")

  out_dir <- file.path(projectDir(proj3), "plots")
  expect_true(file.exists(file.path(out_dir, "network_layers_FLO-ET_fr_layer.png")))
  expect_true(file.exists(file.path(out_dir, "network_layers_FLO-ET_bipartite_node_type.png")))

  log_df <- patliRResults(proj3, "network_layers_plot_log")
  expect_equal(nrow(log_df[log_df$condition == "FLO-ET", ]), 2)
})

test_that("plot_network_layers() layout = 'bipartite' returns a ggplot on a two-layer scope", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p <- plot_network_layers(proj, condition = "FLO-ET", layout = "bipartite", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that(".network_layers_bipartite_coords() places compounds at x = 0, targets at x = 1, ordered by degree", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("c1", "c1", "c2"), to = c("t1", "t2", "t1"), stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = data.frame(name = c("c1", "c2", "t1", "t2"), layer = c("compound", "compound", "target", "target"), stringsAsFactors = FALSE)
  )
  coords <- patliR:::.network_layers_bipartite_coords(g)
  vlayer <- igraph::V(g)$layer
  expect_true(all(coords[vlayer == "compound", 1] == 0))
  expect_true(all(coords[vlayer == "target", 1] == 1))
  ## c1 has degree 2 (t1, t2), c2 has degree 1 -- distinct y within the column
  expect_equal(length(unique(coords[vlayer == "compound", 2])), 2)
})

test_that(".network_layers_bipartite_coords() aborts when the graph has more than two layers", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("c1", "t1"), to = c("t1", "p1"), stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = data.frame(name = c("c1", "t1", "p1"), layer = c("compound", "target", "pathway"), stringsAsFactors = FALSE)
  )
  expect_error(patliR:::.network_layers_bipartite_coords(g), "two-layer")
})

test_that("plot_network_layers() colour_by = 'node_type' groups nodes into compound/target", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  p <- plot_network_layers(proj, condition = "FLO-ET", colour_by = "node_type", save = FALSE)
  expect_s3_class(p, "ggplot")

  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 5, seed = 1, colour_by = "node_type")
  ## default layers = compound/target only -> node_type grouping equals layer exactly
  expect_equal(pd$nodes$colour_group, pd$nodes$layer)
})

test_that("plot_network_layers() colour_by = 'module' errors clearly without network_module_robustness()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(
    plot_network_layers(proj, condition = "FLO-ET", colour_by = "module", save = FALSE),
    "network_module_membership"
  )
})

test_that("plot_network_layers() colour_by = 'module' joins network_module_membership by node_id", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  p <- plot_network_layers(proj, condition = "FLO-ET", colour_by = "module", save = FALSE)
  expect_s3_class(p, "ggplot")

  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 5, seed = 1, colour_by = "module")
  memb <- patliRResults(proj, "network_module_membership")
  lookup <- stats::setNames(memb$module_id, memb$node_id)
  expect_equal(pd$nodes$colour_group, unname(lookup[pd$nodes$name]))
  ## every node in this compound-target graph was clustered -> no "not clustered" fallback here
  expect_false(any(pd$nodes$colour_group == "not clustered"))
})

test_that(".network_layers_module_colour_info() falls back to 'not clustered' for nodes absent from the membership table", {
  proj <- .network_stats_test_setup()
  fake_memb <- data.frame(
    condition = "FLO-ET", node_id = "A", node_type = "compound",
    module_id = "M1", module_type = "cluster", stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_module_membership") <- fake_memb
  info <- patliR:::.network_layers_module_colour_info(proj, "FLO-ET", c("A", "B"))
  expect_equal(info$group, c("M1", "not clustered"))
  expect_equal(unname(info$palette["not clustered"]), "grey70")
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
  ## is_hub is "degree >= the top_hub_n'th highest degree", so a tie at that
  ## cutoff is kept, not arbitrarily broken -- top_hub_n is a floor on how
  ## many nodes get marked, not a hard cap. On this fixture (after
  ## prep_binarize()'s Q1-threshold fix enlarged it) there is a 5-way tie
  ## at the cutoff degree, so is_hub is TRUE for 6 nodes with top_hub_n = 2.
  expect_equal(sum(pd$nodes$is_hub), 6)
  expect_true(all(pd$nodes$degree[pd$nodes$is_hub] >= sort(pd$nodes$degree, decreasing = TRUE)[2]))
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

test_that("pooled module colours distinguish condition-local module IDs", {
  proj <- .test_project()
  patliRResults(proj, "network_module_membership") <- data.frame(
    condition = c("A", "B", "B"), node_id = c("c1", "c2", "c1"), module_id = "M1"
  )
  info <- patliR:::.network_layers_module_colour_info(proj, c("A", "B"), c("c1", "c2", "missing"))
  expect_equal(info$group, c("A:M1", "B:M1", "not clustered"))
  expect_false(identical(unname(info$palette["A:M1"]), unname(info$palette["B:M1"])))
  single <- patliR:::.network_layers_module_colour_info(proj, "B", c("c1", "c2"))
  expect_equal(single$group, c("M1", "M1"))
})

test_that("zero hubs disables labels and invalid hub counts or empty scopes error clearly", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", "bipartite", 0, 1)
  expect_identical(pd$nodes$is_hub, rep(FALSE, nrow(pd$nodes)))
  expect_no_error(plot_network_layers(proj, condition = "FLO-ET", top_hub_n = 0, save = FALSE))
  for (n in list(-1, 0.5, NA_real_, Inf, numeric(), c(1, 2), "1")) {
    expect_error(patliR:::.network_layers_plot_data(proj, g, "FLO-ET", "bipartite", n, 1), "non-negative integer")
  }
  expect_error(plot_network_layers(proj, condition = character(), save = FALSE), "at least one condition")
  expect_error(patliR:::.network_layered_graph_multi(proj, character()), "at least one condition")
})

test_that("min_target_degree keeps only shared targets and drops the compounds left isolated", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  g <- patliR:::.network_layers_filter(proj, g, "FLO-ET", layers = c("compound", "target"), max_pathways = 0)
  expect_identical(patliR:::.network_layers_subset(proj, g, "FLO-ET"), g)

  g2 <- patliR:::.network_layers_subset(proj, g, "FLO-ET", min_target_degree = 2)
  ct <- igraph::as_data_frame(g, what = "edges")
  n_hit <- tapply(ct$from, ct$to, function(x) length(unique(x)))
  kept_targets <- igraph::V(g2)$name[igraph::V(g2)$layer == "target"]
  expect_setequal(kept_targets, names(n_hit)[n_hit >= 2])
  expect_true(all(igraph::degree(g2) > 0))
})

test_that("compound_ids draws only that compound's sub-network, in its own file and log row", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  one <- unique(edges$compound_id[edges$condition == "FLO-ET"])[1]
  p <- plot_network_layers(proj, condition = "FLO-ET", compound_ids = one, save = TRUE)
  pts <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_identical(unique(pts$name[pts$layer == "compound"]), one)
  expect_true(all(pts$name[pts$layer == "target"] %in% edges$uniprot_id[edges$compound_id == one]))

  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  tag <- patliR:::.network_layers_subset_tag(compound_ids = one)
  expect_true(startsWith(tag, "cmp-"))
  expect_identical(log_df$subset, tag)
  expect_true(file.exists(file.path(projectDir(proj), "plots", paste0("network_layers_FLO-ET_fr_layer_", tag, ".png"))))
  expect_match(p$labels$title, "1 selected compound")
})

test_that("module_id restricts to one module and needs network_module_robustness()", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_network_layers(proj, condition = "FLO-ET", module_id = "M1", save = FALSE), "network_module_robustness")
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  mem <- patliRResults(proj, "network_module_membership")
  m <- stats::na.omit(mem$module_id[mem$condition == "FLO-ET"])[1]
  p <- plot_network_layers(proj, condition = "FLO-ET", module_id = m, colour_by = "module", save = FALSE)
  pts <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_true(all(pts$name %in% mem$node_id[mem$condition == "FLO-ET" & mem$module_id %in% m]))
  expect_error(plot_network_layers(proj, condition = "FLO-ET", module_id = "no_such_module", save = FALSE), "No node")
})

test_that("an old network_layers_plot_log without `subset` is back-filled and the default row replaced", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_layers_plot_log") <- data.frame(
    condition = "FLO-ET", layout = "fr", colour_by = "layer", path = "old.png", n_nodes = 1, n_edges = 1
  )
  p <- plot_network_layers(proj, condition = "FLO-ET", save = TRUE)
  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  expect_equal(nrow(log_df), 1)
  expect_identical(log_df$subset, "")
  expect_false(identical(log_df$path, "old.png"))
})

test_that("hub labels are truncated for drawing but the tooltip keeps the full name", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  pd <- patliR:::.network_layers_plot_data(proj, g, "FLO-ET", layout = "fr", top_hub_n = 50, seed = 1)
  pd$nodes$label[1] <- strrep("Long-compound-name-", 4)
  pd$nodes$is_hub[1] <- TRUE
  p <- patliR:::.network_layers_ggplot(pd, engine = "static", title_suffix = "X")
  lab <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p$layers)[[1]]$data
  expect_true(all(nchar(lab$short_label) <= 28))
  expect_identical(lab$label[lab$name == pd$nodes$name[1]], pd$nodes$label[1])
})
