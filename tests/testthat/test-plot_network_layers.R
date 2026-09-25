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

## ---- simplified views (view = "compound" / "bipartite" / "module" / "pathway") ----

.view_ct <- function() {
  ## t1: A,B,C  t2: A,B  t3: A only  t4: C only  t5: B,C
  data.frame(
    compound_id = c("A", "B", "C", "A", "B", "A", "C", "B", "C"),
    uniprot_id  = c("t1", "t1", "t1", "t2", "t2", "t3", "t4", "t5", "t5"),
    weight      = c(0.9, 0.8, 0.7, 0.6, 0.6, 0.5, 0.55, 0.95, 0.65),
    stringsAsFactors = FALSE
  )
}

test_that(".network_view_collapse() keeps targets shared by >= k compounds and counts the rest in badges", {
  col <- patliR:::.network_view_collapse(.view_ct(), min_shared = 2)
  expect_identical(col$k, 2L)
  expect_setequal(col$targets$uniprot_id[col$targets$drawn], c("t1", "t2", "t5"))
  expect_equal(nrow(col$edges), 7)
  expect_true(all(col$edges$uniprot_id %in% c("t1", "t2", "t5")))
  b <- col$badges[match(c("A", "B", "C"), col$badges$compound_id), ]
  expect_equal(b$n_targets, c(3, 3, 3))
  expect_equal(b$n_unique, c(1, 0, 1))
  expect_equal(b$n_drawn, c(2, 3, 2))
  expect_equal(b$n_hidden, b$n_targets - b$n_drawn)
  ## targets sorted by sharing, then mean probability
  expect_identical(col$targets$uniprot_id[1], "t1")
  expect_false(is.unsorted(rev(col$targets$n_compounds)))
  expect_identical(col$targets$uniprot_id[2:3], c("t5", "t2")) # t5 mean 0.80 > t2 mean 0.60
})

test_that(".network_view_collapse() picks the smallest k >= 2 within max_targets automatically", {
  ct <- .view_ct()
  expect_identical(patliR:::.network_view_collapse(ct, max_targets = 5)$k, 2L)
  expect_identical(patliR:::.network_view_collapse(ct, max_targets = 1)$k, 3L)
  col <- patliR:::.network_view_collapse(ct, max_targets = 1)
  expect_identical(col$targets$uniprot_id[col$targets$drawn], "t1")
  ## budget unreachable: the largest possible k is used
  ct2 <- rbind(ct, data.frame(compound_id = c("A", "B", "C"), uniprot_id = "t6", weight = 0.7))
  expect_identical(patliR:::.network_view_collapse(ct2, max_targets = 1)$k, 3L)
  ## one compound: nothing is "shared", every target is drawn
  one <- ct[ct$compound_id == "A", ]
  expect_identical(patliR:::.network_view_collapse(one)$k, 1L)
  expect_true(all(patliR:::.network_view_collapse(one)$targets$drawn))
})

test_that(".network_view_compound_order() places compounds with similar target sets next to each other", {
  ct <- data.frame(
    compound_id = c("A", "A", "B", "B", "C", "C", "D", "D", "E"),
    uniprot_id = c("t1", "t2", "t1", "t2", "t3", "t4", "t3", "t4", "t5"),
    stringsAsFactors = FALSE
  )
  o <- patliR:::.network_view_compound_order(ct)
  expect_setequal(o, c("A", "B", "C", "D", "E"))
  expect_equal(abs(diff(match(c("A", "B"), o))), 1)
  expect_equal(abs(diff(match(c("C", "D"), o))), 1)
  expect_identical(patliR:::.network_view_compound_order(ct[ct$compound_id %in% c("A", "E"), ]), c("A", "E"))
})

test_that(".network_view_compound_layout() puts compounds on the ring and more shared targets closer to the centre", {
  col <- patliR:::.network_view_collapse(.view_ct(), min_shared = 2)
  lay <- patliR:::.network_view_compound_layout(col, c("A", "B", "C"))
  expect_equal(sqrt(lay$compounds$x^2 + lay$compounds$y^2), rep(1, 3))
  r <- stats::setNames(sqrt(lay$targets$x^2 + lay$targets$y^2), lay$targets$uniprot_id)
  expect_lt(r[["t1"]], r[["t2"]])
  expect_equal(r[["t2"]], r[["t5"]])
  expect_true(all(r < 1))
  ## same shell, distinct positions
  xy <- lay$targets[lay$targets$uniprot_id %in% c("t2", "t5"), c("x", "y")]
  expect_gt(sqrt(sum((xy[1, ] - xy[2, ])^2)), 0.1)
})

test_that(".network_view_bipartite_layout() groups targets by module (largest first) and sorts them by sharing", {
  col <- patliR:::.network_view_collapse(.view_ct(), min_shared = 1)
  module <- c(t1 = "M2", t2 = "M1", t3 = "M1", t4 = "M1", t5 = NA)
  lay <- patliR:::.network_view_bipartite_layout(col$edges, col$targets, c("A", "B", "C"), module = module)
  tg <- lay$targets
  expect_identical(tg$module, c("M1", "M1", "M1", "M2", "not clustered"))
  expect_identical(tg$uniprot_id[tg$module == "M1"][1], "t2") # hit by 2 compounds, before t3/t4
  expect_false(is.unsorted(rev(tg$y)))                        # top to bottom
  expect_true(all(tg$x == 1) && all(lay$compounds$x == 0))
  ## compounds ordered by the mean height of their targets
  expect_false(is.unsorted(rev(lay$compounds$bary)))
  ## compound modules group the compounds in the target module order
  lay2 <- patliR:::.network_view_bipartite_layout(col$edges, col$targets, c("A", "B", "C"), module = module,
                                                  compound_module = c(A = "M2", B = "M1", C = "M1"))
  expect_identical(lay2$compounds$module, c("M1", "M1", "M2"))
  ## no modules: sorted by sharing only
  lay3 <- patliR:::.network_view_bipartite_layout(col$edges, col$targets, c("A", "B", "C"))
  expect_identical(lay3$targets$uniprot_id[1], "t1")
})

test_that(".network_layers_resolve_view() simplifies only large plain compound-target graphs", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = "c1", to = paste0("t", 1:10)), directed = TRUE,
    vertices = data.frame(name = c("c1", paste0("t", 1:10)), layer = c("compound", rep("target", 10)))
  )
  expect_identical(patliR:::.network_layers_resolve_view("auto", g, auto_max_nodes = 300)$view, "network")
  auto <- patliR:::.network_layers_resolve_view("auto", g, auto_max_nodes = 5)
  expect_identical(auto$view, "compound")
  expect_match(auto$note, "11 nodes")
  expect_identical(patliR:::.network_layers_resolve_view("auto", g, layout_given = TRUE, auto_max_nodes = 5)$view, "network")
  expect_identical(patliR:::.network_layers_resolve_view("bipartite", g, auto_max_nodes = 5), list(view = "bipartite", note = NULL))
  igraph::V(g)$layer[2] <- "pathway"
  expect_identical(patliR:::.network_layers_resolve_view("auto", g, auto_max_nodes = 5)$view, "network")
})

test_that("plot_network_layers() validates the view arguments", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "hairball", save = FALSE), "should be one of")
  for (bad in list(0, 1.5, "2", c(2, 3), NA_real_)) {
    expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "compound", min_shared = bad, save = FALSE), "min_shared")
  }
  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "compound", max_targets = 0, save = FALSE), "max_targets")
  expect_error(plot_network_layers(proj, condition = "FLO-ET", auto_max_nodes = -1, save = FALSE), "auto_max_nodes")
  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "compound", top_hub_n = -1, save = FALSE), "non-negative integer")
})

test_that("view = \"compound\" draws only shared targets, logs its own row, and badge counts match the edges", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  e <- unique(edges[edges$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  n_hit <- table(e$uniprot_id)

  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  g <- patliR:::.network_layers_filter(proj, g, "FLO-ET", layers = c("compound", "target"), max_pathways = 0)
  ct <- patliR:::.network_layers_ct_table(proj, g, "FLO-ET")
  expect_equal(nrow(ct), nrow(e))
  w <- edges[edges$condition == "FLO-ET", ]
  expect_equal(ct$weight, w$weight[match(paste(ct$compound_id, ct$uniprot_id), paste(w$compound_id, w$uniprot_id))])

  built <- patliR:::.network_view_compound(proj, ct, "FLO-ET", min_shared = 2)
  expect_s3_class(built$plot, "ggplot")
  expect_setequal(built$data$targets$uniprot_id, names(n_hit)[n_hit >= 2])
  expect_equal(built$n_edges, sum(e$uniprot_id %in% names(n_hit)[n_hit >= 2]))
  b <- built$data$compounds
  expect_equal(b$n_targets, as.integer(table(e$compound_id)[b$name]))
  single <- names(n_hit)[n_hit == 1]
  expect_equal(b$n_unique, vapply(b$name, function(id) sum(e$uniprot_id[e$compound_id == id] %in% single), integer(1), USE.NAMES = FALSE))
  expect_identical(built$view_tag, "k2")

  p <- plot_network_layers(proj, condition = "FLO-ET", view = "compound", min_shared = 2, save = TRUE)
  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  row <- log_df[log_df$view == "compound", ]
  expect_equal(nrow(row), 1)
  expect_identical(row$subset, "k2")
  expect_true(is.na(row$layout))
  expect_true(file.exists(file.path(projectDir(proj), "plots", "network_layers_FLO-ET_compound_k2.png")))
})

test_that("view = \"auto\" keeps the network figure for small graphs and switches (with a subtitle note) above auto_max_nodes", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  p_small <- plot_network_layers(proj, condition = "FLO-ET", save = FALSE)
  expect_match(p_small$labels$title, "network")
  expect_false(grepl("Simplified view", p_small$labels$subtitle))
  p_big <- plot_network_layers(proj, condition = "FLO-ET", auto_max_nodes = 2, save = FALSE)
  expect_match(p_big$labels$title, "Compound-target summary")
  expect_match(gsub("\n", " ", p_big$labels$subtitle), "Simplified view chosen automatically")
  ## an explicit layout keeps the requested force layout
  p_kk <- plot_network_layers(proj, condition = "FLO-ET", layout = "kk", auto_max_nodes = 2, save = FALSE)
  expect_false(grepl("summary", p_kk$labels$title))
})

test_that("view = \"bipartite\" and view = \"module\" draw every compound-target edge in scope", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  e <- unique(edges[edges$condition == "FLO-ET", c("compound_id", "uniprot_id")])

  p <- plot_network_layers(proj, condition = "FLO-ET", view = "bipartite", save = TRUE)
  expect_s3_class(p, "ggplot")
  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  expect_equal(log_df$n_edges[log_df$view == "bipartite"], nrow(e))
  expect_true(file.exists(file.path(projectDir(proj), "plots", "network_layers_FLO-ET_bipartite.png")))

  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "module", save = FALSE), "network_module_robustness")
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  g <- patliR:::.network_layered_graph_multi(proj, "FLO-ET")
  ct <- patliR:::.network_layers_ct_table(proj, g, "FLO-ET")
  built <- patliR:::.network_view_module(proj, ct, "FLO-ET")
  expect_s3_class(built$plot, "ggplot")
  expect_equal(built$n_edges, nrow(e))
  expect_equal(built$n_nodes, length(unique(c(e$compound_id, e$uniprot_id))))
  mem <- patliRResults(proj, "network_module_membership")
  mem <- mem[mem$condition == "FLO-ET", ]
  nd <- built$data$nodes
  expected <- mem$module_id[match(nd$name, mem$node_id)]
  expect_equal(nd$module, ifelse(is.na(expected), "not clustered", expected))
  ## bipartite with modules: each module is one contiguous block of targets
  b <- patliR:::.network_view_bipartite(proj, ct, "FLO-ET")
  blocks <- rle(b$data$targets$module)$values
  expect_equal(length(blocks), length(unique(blocks)))
})

test_that(".network_view_pathway_cells() counts each compound's targets per term", {
  ct <- .view_ct()
  cells <- patliR:::.network_view_pathway_cells(ct, list(P1 = c("t1", "t2"), P2 = "t4", P3 = "zzz"))
  expect_setequal(unique(cells$ID), c("P1", "P2"))
  p1 <- cells[cells$ID == "P1", ]
  expect_equal(p1$n_hit[match(c("A", "B", "C"), p1$compound_id)], c(2, 2, 1))
  expect_true(all(p1$n_term_targets == 2))
  expect_equal(p1$frac[p1$compound_id == "C"], 0.5)
  expect_equal(nrow(patliR:::.network_view_pathway_cells(ct, list(P3 = "zzz"))), 0)
})

test_that("view = \"pathway\" needs network_enrich() results", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "pathway", save = FALSE), "network_enrich")
})

test_that(".network_view_hulls() returns one padded polygon per group", {
  h <- patliR:::.network_view_hulls(c(0, 1, 5), c(0, 0, 5), c("a", "a", "b"), pad = 0.2)
  expect_setequal(unique(h$group), c("a", "b"))
  ## a single-node group still gets a polygon around its point
  hb <- h[h$group == "b", ]
  expect_gte(nrow(hb), 8)
  expect_equal(max(sqrt((hb$x - 5)^2 + (hb$y - 5)^2)), 0.2, tolerance = 1e-8)
})

test_that("an old network_layers_plot_log without `view` is back-filled as the network view", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_layers_plot_log") <- data.frame(
    condition = "FLO-ET", layout = "fr", colour_by = "layer", subset = "", path = "old.png", n_nodes = 1, n_edges = 1
  )
  p <- plot_network_layers(proj, condition = "FLO-ET", save = TRUE)
  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  expect_equal(nrow(log_df), 1)
  expect_identical(log_df$view, "network")
  expect_false(identical(log_df$path, "old.png"))
})

test_that("view = \"pathway\" keeps the max_pathways most significant terms and orders the dots by clustering", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  up <- unique(edges$uniprot_id[edges$condition == "FLO-ET"])
  map <- suppressMessages(suppressWarnings(tryCatch(
    clusterProfiler::bitr(up, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE),
    error = function(e) NULL
  )))
  testthat::skip_if(is.null(map) || length(unique(map$ENTREZID)) < 4, "fixture targets do not map to Entrez")
  ez <- unique(map$ENTREZID)
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = "FLO-ET", db = c("kegg", "kegg", "kegg", "go"),
    ID = c("hsa1", "hsa2", "hsa3", "GO:1"), Description = c("Term one", "Term two", "Term three", "GO term"),
    p.adjust = c(0.01, 0.001, 0.05, 1e-6),
    geneID = c(paste(ez[1:2], collapse = "/"), paste(ez, collapse = "/"), ez[3], ez[1]),
    stringsAsFactors = FALSE
  )
  p <- suppressMessages(plot_network_layers(proj, condition = "FLO-ET", view = "pathway", pathway_db = "kegg",
                                            max_pathways = 2, save = TRUE))
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(as.character(p$data$ID)), c("hsa1", "hsa2"))
  expect_true(all(p$data$frac > 0 & p$data$frac <= 1))
  log_df <- patliRResults(attr(p, "proj"), "network_layers_plot_log")
  expect_identical(log_df$subset[log_df$view == "pathway"], "top2_kegg")
  expect_true(file.exists(file.path(projectDir(proj), "plots", "network_layers_FLO-ET_pathway_top2_kegg.png")))
  expect_error(plot_network_layers(proj, condition = "FLO-ET", view = "pathway", pathway_db = "reactome", save = FALSE), "pathway_db")
})

test_that("module colours are colour-blind safe and stable across views and subsets", {
  full <- patliR:::.network_view_qual_palette(c("M10", "M2", "M1", "not clustered"))
  sub <- patliR:::.network_view_qual_palette(c("M2", "M1"))
  expect_identical(unname(full[c("M1", "M2")]), unname(sub[c("M1", "M2")]))
  expect_identical(unname(full["M1"]), "#0072B2")
  expect_identical(unname(full["not clustered"]), "grey70")
  ## the palette covers every module in scope, not just the nodes asked about
  proj <- .test_project()
  patliRResults(proj, "network_module_membership") <- data.frame(
    condition = "A", node_id = c("c1", "c2", "c3"), module_id = c("M1", "M2", "M3")
  )
  info <- patliR:::.network_layers_module_colour_info(proj, "A", "c3")
  expect_true(all(c("M1", "M2", "M3") %in% names(info$palette)))
  expect_identical(unname(info$palette["M3"]), "#009E73")
})
