## Two K_{3,3} blocks. `bridge = TRUE` joins them by a single edge (the
## spec's graph); `bridge = FALSE` leaves them disjoint. Compounds = *c*,
## targets = *t*.
.two_k33 <- function(bridge = TRUE) {
  block <- function(pre) {
    expand.grid(
      c = paste0(pre, c("c1", "c2", "c3")),
      t = paste0(pre, c("t1", "t2", "t3")),
      stringsAsFactors = FALSE
    )
  }
  el <- rbind(block("A"), block("B"))
  g <- igraph::graph_from_data_frame(el, directed = FALSE)
  if (bridge) g <- igraph::add_edges(g, c("At1", "Bc1"))
  igraph::V(g)$type <- grepl("t[0-9]$", igraph::V(g)$name)
  igraph::E(g)$weight <- 0.9
  g
}
.two_k33_joined <- function() .two_k33(bridge = TRUE)

.block_of <- function(node_name) substr(node_name, 1, 1)
.blocks_separated <- function(members) {
  all(vapply(members, function(m) length(unique(vapply(m, .block_of, character(1)))) == 1L, logical(1)))
}

test_that("network_module_robustness() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_module_robustness(proj), "network_edges")
})

test_that("network_module_robustness() produces a well-formed summary table on FLO-ET (default clustering = leiden)", {
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)

  summary_tbl <- patliRResults(proj, "network_module_robustness")
  expect_true(all(c(
    "condition", "module_id", "module_type", "n_nodes", "r_index", "r_index_random",
    "n_random", "clustering", "modularity", "resolution", "method_detail", "seed_used"
  ) %in% names(summary_tbl)))
  expect_true(all(summary_tbl$clustering == "leiden"))
  expect_true(all(summary_tbl$r_index >= 0 & summary_tbl$r_index <= 0.5))
  ## attack defaults to "targeted" -> no random baseline
  expect_true(all(is.na(summary_tbl$r_index_random)))
})

test_that("network_module_robustness() writes network_module_membership, one row per node", {
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", seed = 42)
  g <- .network_graph(proj, "FLO-ET")

  memb <- patliRResults(proj, "network_module_membership")
  expect_setequal(names(memb), c("condition", "node_id", "node_type", "module_id", "module_type"))
  expect_setequal(memb$node_id, igraph::V(g)$name)
  expect_equal(nrow(memb), igraph::vcount(g))
  expect_true(all(memb$node_type %in% c("compound", "target")))
  ## membership module_id must join back to the summary module_id
  summ <- patliRResults(proj, "network_module_robustness")
  expect_true(all(memb$module_id %in% summ$module_id))
})

for (cl in c("leiden", "bipartite", "hdbscan")) {
  test_that(sprintf("clustering = '%s': two disjoint K_{3,3} blocks land in different modules", cl), {
    if (cl == "bipartite") skip_if_not_installed("bipartite")
    if (cl == "hdbscan") skip_if_not_installed("dbscan")
    g <- .two_k33(bridge = FALSE)
    mods <- patliR:::.network_detect_modules(
      g, clustering = cl, min_module_size = 2, resolution = 1,
      n_iterations = 5L, min_component_size = 3L, seed = 1L
    )
    all_nodes <- unlist(mods$members, use.names = FALSE)
    expect_setequal(all_nodes, igraph::V(g)$name)
    expect_equal(length(all_nodes), length(unique(all_nodes)))
    expect_true(.blocks_separated(mods$members))
    expect_gte(length(mods$members), 2L)
  })
}

test_that("clustering = 'leiden'/'bipartite': a single bridge edge does not merge two K_{3,3} blocks", {
  ## hdbscan on integer graph distances degenerates to single linkage and
  ## cannot resist the bridge here -- documented in the roxygen; only the
  ## two modularity backends are asserted.
  g <- .two_k33(bridge = TRUE)
  ml <- patliR:::.network_detect_modules(g, clustering = "leiden", seed = 1L)
  expect_true(.blocks_separated(ml$members))
  expect_gte(length(ml$members), 2L)

  skip_if_not_installed("bipartite")
  mb <- patliR:::.network_detect_modules(g, clustering = "bipartite", seed = 1L)
  expect_true(.blocks_separated(mb$members))
  expect_gte(length(mb$members), 2L)
})

for (cl in c("leiden", "bipartite", "hdbscan")) {
  for (atk in c("targeted", "random", "both")) {
    test_that(sprintf("sum(n_nodes) == vcount(g) and r_index <= 0.5 -- clustering='%s' attack='%s'", cl, atk), {
      if (cl == "bipartite") skip_if_not_installed("bipartite")
      if (cl == "hdbscan") skip_if_not_installed("dbscan")
      proj <- .network_stats_test_setup()
      proj <- network_module_robustness(
        proj, condition = "FLO-ET", clustering = cl, attack = atk,
        n_random = 8L, seed = 42
      )
      g <- .network_graph(proj, "FLO-ET")
      summ <- patliRResults(proj, "network_module_robustness")
      expect_equal(sum(summ$n_nodes), igraph::vcount(g))
      expect_true(all(summ$r_index <= 0.5, na.rm = TRUE))
      expect_true(all(summ$r_index_random <= 0.5, na.rm = TRUE))
      if (atk == "targeted") expect_true(all(is.na(summ$r_index_random)))
      if (atk == "random") expect_true(all(is.na(summ$r_index)))
      if (atk == "both") expect_true(all(!is.na(summ$r_index) & !is.na(summ$r_index_random)))
    })
  }
}

test_that("Leiden with objective_function = 'modularity' escapes the CPM singleton trap", {
  ## two K4 cliques joined by one edge
  g <- igraph::make_full_graph(4) + igraph::make_full_graph(4)
  g <- igraph::set_vertex_attr(g, "name", value = paste0("v", seq_len(8)))
  g <- igraph::add_edges(g, c("v4", "v5"))
  mods <- patliR:::.network_detect_modules_leiden(g, resolution = 1, n_iterations = 5L)
  expect_lt(length(mods$members), igraph::vcount(g))   # NOT all singletons
  memb <- integer(8)
  for (k in seq_along(mods$members)) memb[match(mods$members[[k]], paste0("v", 1:8))] <- k
  expect_gt(igraph::modularity(g, memb, weights = NA), 0.3)
})

test_that("module2constraints() round-trip survives a zero-sum row/column (empty() bug)", {
  skip_if_not_installed("bipartite")
  ## 2 clean blocks + one isolated compound (all-zero row) + isolated target (all-zero col)
  mat <- rbind(
    Ac1 = c(At1 = 1, At2 = 1, Bt1 = 0, Bt2 = 0, Zt = 0),
    Ac2 = c(1, 1, 0, 0, 0),
    Bc1 = c(0, 0, 1, 1, 0),
    Bc2 = c(0, 0, 1, 1, 0),
    Zc  = c(0, 0, 0, 0, 0)
  )
  colnames(mat) <- c("At1", "At2", "Bt1", "Bt2", "Zt")
  bm <- patliR:::.network_bipartite_membership(mat, seed = 7L)
  m <- bm$membership
  expect_equal(length(m), nrow(mat) + ncol(mat))
  ## isolated nodes -> NA
  expect_true(is.na(m[["Zc"]]))
  expect_true(is.na(m[["Zt"]]))
  ## the two real blocks recovered correctly
  expect_equal(unname(m[["Ac1"]]), unname(m[["Ac2"]]))
  expect_equal(unname(m[["Bc1"]]), unname(m[["Bc2"]]))
  expect_false(unname(m[["Ac1"]]) == unname(m[["Bc1"]]))
  expect_equal(unname(m[["Ac1"]]), unname(m[["At1"]]))
  expect_equal(unname(m[["Bc1"]]), unname(m[["Bt1"]]))
})

test_that("same seed -> identical membership, for Leiden and bipartite", {
  g <- .two_k33_joined()
  l1 <- patliR:::.network_detect_modules(g, clustering = "leiden", seed = 99L)
  l2 <- patliR:::.network_detect_modules(g, clustering = "leiden", seed = 99L)
  expect_identical(l1$members, l2$members)

  skip_if_not_installed("bipartite")
  b1 <- patliR:::.network_detect_modules(g, clustering = "bipartite", seed = 99L)
  b2 <- patliR:::.network_detect_modules(g, clustering = "bipartite", seed = 99L)
  expect_identical(b1$members, b2$members)
})

test_that("network_module_robustness() does not leak RNG state -- leiden, bipartite, random attack", {
  proj <- .network_stats_test_setup()
  for (cl in c("leiden", "bipartite")) {
    if (cl == "bipartite" && !requireNamespace("bipartite", quietly = TRUE)) next
    set.seed(999); before <- runif(1)
    set.seed(999)
    network_module_robustness(proj, condition = "FLO-ET", clustering = cl,
                              attack = "both", n_random = 5L, seed = 42)
    after <- runif(1)
    expect_equal(before, after)
  }
})

test_that("clustering = 'hdbscan' reproduces the pre-0.2.0 output given the same seed", {
  skip_if_not_installed("dbscan")
  proj <- .network_stats_test_setup()
  p1 <- network_module_robustness(proj, condition = "FLO-ET", clustering = "hdbscan", seed = 42)
  p2 <- network_module_robustness(proj, condition = "FLO-ET", clustering = "hdbscan", seed = 42)
  s1 <- patliRResults(p1, "network_module_robustness")
  s2 <- patliRResults(p2, "network_module_robustness")
  expect_equal(s1$r_index, s2$r_index)
  expect_equal(s1$module_id, s2$module_id)
  expect_equal(s1$n_nodes, s2$n_nodes)
  ## frozen against patliR 2a27dde (pre-0.2.0): FLO-ET / seed 42 / hdbscan
  expect_equal(s1$module_id, c("M1", "M2", "M3"))
  expect_equal(s1$module_type, c("cluster", "cluster", "noise"))
  expect_equal(s1$n_nodes, c(4L, 2L, 8L))
  expect_equal(round(s1$r_index, 6), c(0.187500, 0.250000, 0.109375))
  ## the historical module_id / n_nodes scheme is unchanged: modules cover
  ## every node, exactly once
  g <- .network_graph(p1, "FLO-ET")
  expect_equal(sum(s1$n_nodes), igraph::vcount(g))
  expect_true(all(is.na(s1$modularity)))
})

test_that("attack = 'random' curve sits above the targeted curve on average (Schneider contrast)", {
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(
    proj, condition = "FLO-ET", clustering = "leiden", attack = "both",
    n_random = 20L, seed = 42
  )
  summ <- patliRResults(proj, "network_module_robustness")
  ## random failure fragments a module slower than a hub-targeted attack
  expect_gt(mean(summ$r_index_random), mean(summ$r_index))
})

test_that(".network_is_bipartite(): TRUE on a compound-target graph, FALSE within-mode, no throw on named non-bipartite", {
  proj <- .network_stats_test_setup()
  g <- .network_graph(proj, "FLO-ET")
  expect_true(patliR:::.network_is_bipartite(g))

  ## add a compound-compound edge -> no longer bipartite
  cmp <- igraph::V(g)$name[!as.logical(igraph::V(g)$type)][1:2]
  g_bad <- igraph::add_edges(g, cmp)
  expect_false(patliR:::.network_is_bipartite(g_bad))

  ## a named, clearly non-bipartite graph (triangle) must return, not throw
  tri <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")), directed = FALSE
  )
  expect_silent(res <- patliR:::.network_is_bipartite(tri))
  expect_false(res)
})

test_that("clustering = 'bipartite' aborts on a non-bipartite graph before touching the bipartite package", {
  skip_if_not_installed("bipartite")
  proj <- .network_stats_test_setup()
  g <- .network_graph(proj, "FLO-ET")
  cmp <- igraph::V(g)$name[!as.logical(igraph::V(g)$type)][1:2]
  g_bad <- igraph::add_edges(g, cmp)
  .save_cache_graph(g_bad, proj, "FLO-ET")
  expect_error(
    network_module_robustness(proj, condition = "FLO-ET", clustering = "bipartite", seed = 1),
    "two-mode"
  )
})

test_that("clustering = 'bipartite' runs end-to-end on a star (1xk) component (B1)", {
  skip_if_not_installed("bipartite")
  proj <- .network_stats_test_setup()
  g <- .network_graph(proj, "FLO-ET")
  ## add an isolated compound with two targets nothing else touches: a 1x2
  ## star component. bipartite::computeModules() errors on a 1-row/1-col
  ## biadjacency matrix, so this used to abort the whole run.
  g2 <- igraph::add_vertices(
    g, 3, attr = list(name = c("Cstar", "TstarA", "TstarB"), type = c(FALSE, TRUE, TRUE))
  )
  g2 <- igraph::add_edges(g2, c("Cstar", "TstarA", "Cstar", "TstarB"))
  .save_cache_graph(g2, proj, "FLO-ET")

  for (cl in c("leiden", "bipartite", "hdbscan")) {
    if (cl == "hdbscan") skip_if_not_installed("dbscan")
    expect_no_error(
      p <- network_module_robustness(proj, condition = "FLO-ET", clustering = cl, seed = 1)
    )
    s <- patliRResults(p, "network_module_robustness")
    expect_equal(sum(s$n_nodes), igraph::vcount(g2))
  }
})

test_that("clustering = 'leiden' with an NA edge weight does not collapse to all-singletons", {
  ## behaviour is correct via delete_edge_attr; this pins it (Phase G item 18)
  g <- .two_k33(bridge = TRUE)
  igraph::E(g)$weight[1] <- NA
  mods <- patliR:::.network_detect_modules(g, clustering = "leiden", seed = 1L)
  expect_lt(length(mods$members), igraph::vcount(g))
  expect_true(.blocks_separated(mods$members))
})

test_that("a near-.Machine$integer.max seed does not overflow the random sub-seed (S2)", {
  proj <- .network_stats_test_setup()
  expect_no_error(
    network_module_robustness(
      proj, condition = "FLO-ET", clustering = "leiden", attack = "both",
      n_random = 4L, seed = .Machine$integer.max - 3L
    )
  )
})

test_that("random robustness does not reuse replicate seeds across modules above 1000 draws", {
  proj <- .network_stats_test_setup()
  .save_cache_graph(.two_k33(bridge = FALSE), proj, "FLO-ET")
  seeds <- integer(0)
  with_seed <- patliR:::.with_seed
  local_mocked_bindings(
    .with_seed = function(seed) {
      seeds <<- c(seeds, seed)
      with_seed(seed)
    },
    .package = "patliR"
  )
  proj <- network_module_robustness(
    proj, condition = "FLO-ET", attack = "random", min_component_size = 7L,
    n_random = 1001L, seed = .Machine$integer.max - 3L
  )
  result <- patliRResults(proj, "network_module_robustness")
  expect_equal(nrow(result), 2L)
  ## One clustering seed, then 1001 independent draws for each module.
  expect_length(seeds, 2003L)
  expect_equal(anyDuplicated(seeds[-1]), 0L)
  expect_true(all(is.finite(result$r_index_random)))
})

test_that("network_module_robustness() rejects n_random < 2 for a random attack (S3)", {
  proj <- .network_stats_test_setup()
  expect_error(
    network_module_robustness(proj, condition = "FLO-ET", attack = "random", n_random = 1L),
    "n_random"
  )
  expect_error(
    network_module_robustness(proj, condition = "FLO-ET", attack = "both", n_random = 1L),
    "n_random"
  )
  ## still fine for the default targeted attack
  expect_no_error(
    network_module_robustness(proj, condition = "FLO-ET", n_random = 1L, seed = 1)
  )
})

test_that(".network_detect_modules() falls back to a single module below min_component_size", {
  tiny_g <- igraph::graph_from_data_frame(
    data.frame(from = "a", to = "b", stringsAsFactors = FALSE), directed = FALSE
  )
  modules <- patliR:::.network_detect_modules(tiny_g, clustering = "hdbscan", min_module_size = 2)
  expect_equal(names(modules$members), "M1")
  expect_setequal(modules$members[["M1"]], c("a", "b"))
  expect_true(grepl("single module", modules$log_message))
  expect_equal(unname(modules$types), "cluster")
})

test_that(".network_detect_modules() clusters disconnected graphs per connected component", {
  skip_if_not_installed("dbscan")
  disconnected_edges <- data.frame(
    from = c("a", "b", "c", "d", "e", "f", "x"),
    to   = c("b", "c", "d", "e", "f", "a", "y"),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(disconnected_edges, directed = FALSE)
  expect_equal(igraph::components(g)$no, 2)

  modules <- patliR:::.network_detect_modules(g, clustering = "hdbscan", min_module_size = 2)
  all_nodes <- unlist(modules$members, use.names = FALSE)
  expect_setequal(all_nodes, igraph::V(g)$name)
  expect_equal(length(all_nodes), length(unique(all_nodes)))
  expect_true(grepl("2 connected component", modules$log_message))
  expect_true(grepl("component with 2 node", modules$log_message))
})

test_that("network_module_robustness() writes a matching fragmentation curve (targeted)", {
  proj <- .network_stats_test_setup()
  proj <- network_module_robustness(proj, condition = "FLO-ET", clustering = "leiden", seed = 42)

  summary_tbl <- patliRResults(proj, "network_module_robustness")
  curve <- patliRResults(proj, "network_robustness_curve")
  expect_true(all(c(
    "condition", "module_id", "n_removed", "largest_component_fraction",
    "removal_strategy", "largest_component_sd", "n_replicates"
  ) %in% names(curve)))

  for (m in unique(summary_tbl$module_id)) {
    n_nodes <- summary_tbl$n_nodes[summary_tbl$module_id == m]
    curve_m <- curve[curve$module_id == m & curve$removal_strategy == "targeted", ]
    expect_equal(nrow(curve_m), n_nodes + 1)
    expect_equal(curve_m$largest_component_fraction[curve_m$n_removed == 0], 1)
    post_removal <- curve_m$largest_component_fraction[curve_m$n_removed > 0]
    expect_equal(mean(post_removal), summary_tbl$r_index[summary_tbl$module_id == m])
  }
})

test_that("network_module_robustness() is reproducible given the same seed", {
  proj <- .network_stats_test_setup()
  proj1 <- network_module_robustness(proj, condition = "FLO-ET", seed = 123)
  proj2 <- network_module_robustness(proj, condition = "FLO-ET", seed = 123)
  expect_equal(
    patliRResults(proj1, "network_module_robustness")$r_index,
    patliRResults(proj2, "network_module_robustness")$r_index
  )
})

test_that("network_module_robustness() rejects min_module_size < 2", {
  proj <- .network_stats_test_setup()
  expect_error(network_module_robustness(proj, condition = "FLO-ET", min_module_size = 1))
})

test_that("a legacy network_module_robustness.csv (pre-0.2.0 schema) is backfilled, not aborted", {
  proj <- .network_stats_test_setup()
  ## simulate an old-schema results bag + CSV
  old <- data.frame(
    condition = "FLO-ET", module_id = "M1", module_type = "cluster",
    n_nodes = 14L, r_index = 0.2, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_module_robustness") <- old
  dir.create(file.path(projectDir(proj), "results"), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(old, file.path(projectDir(proj), "results", "network_module_robustness.csv"), row.names = FALSE)

  expect_no_error(
    proj <- network_module_robustness(proj, condition = "LEA-ET", clustering = "leiden", seed = 1)
  )
  s <- patliRResults(proj, "network_module_robustness")
  expect_true(all(c("clustering", "modularity", "resolution", "method_detail", "r_index_random") %in% names(s)))
  ## the untouched legacy FLO-ET row is kept, backfilled NA on the new cols
  flo <- s[s$condition == "FLO-ET", ]
  expect_equal(nrow(flo), 1L)
  expect_true(is.na(flo$clustering))
})

test_that("a mixed project (one < 2-node condition + one normal) does not fail the rbind, all clusterings", {
  for (cl in c("leiden", "bipartite", "hdbscan")) {
    if (cl == "bipartite" && !requireNamespace("bipartite", quietly = TRUE)) next
    if (cl == "hdbscan" && !requireNamespace("dbscan", quietly = TRUE)) next
    proj <- .network_stats_test_setup()
    ## force one condition's cached graph to a single node
    one_node <- igraph::make_empty_graph(n = 1, directed = FALSE)
    one_node <- igraph::set_vertex_attr(one_node, "name", value = "solo")
    one_node <- igraph::set_vertex_attr(one_node, "type", value = FALSE)
    .save_cache_graph(one_node, proj, "LEA-AQ")
    expect_no_error(
      proj <- network_module_robustness(proj, clustering = cl, seed = 1)
    )
    s <- patliRResults(proj, "network_module_robustness")
    expect_true("FLO-ET" %in% s$condition)
  }
})
