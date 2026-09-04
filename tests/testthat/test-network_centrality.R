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

## ---- bipartite normalisation (spec 2.3) -----------------------------

## build a bipartite igraph the way network_build() does (type = FALSE for
## compounds, TRUE for targets)
.bip_graph <- function(comp_ids, targ_ids, edges) {
  v <- data.frame(
    name = c(comp_ids, targ_ids),
    type = c(rep(FALSE, length(comp_ids)), rep(TRUE, length(targ_ids))),
    stringsAsFactors = FALSE
  )
  igraph::graph_from_data_frame(
    as.data.frame(edges, stringsAsFactors = FALSE), directed = FALSE, vertices = v
  )
}

test_that("degree_norm and betweenness_norm stay in [0, 1] on random bipartite graphs (factor-of-2 guard)", {
  set.seed(42)
  for (i in seq_len(25)) {
    n <- sample(1:15, 1); m <- sample(1:15, 1)
    comp <- paste0("C", seq_len(n)); targ <- paste0("T", seq_len(m))
    grid <- expand.grid(c = comp, t = targ, stringsAsFactors = FALSE)
    grid <- grid[stats::runif(nrow(grid)) < 0.4, , drop = FALSE]
    if (nrow(grid) == 0) grid <- data.frame(c = comp[1], t = targ[1], stringsAsFactors = FALSE)
    g <- .bip_graph(comp, targ, grid)

    df <- patliR:::.network_centrality_one(g, "x", c("degree", "betweenness", "hub_score"), TRUE)
    dn <- df$degree_norm[!is.na(df$degree_norm)]
    bn <- df$betweenness_norm[!is.na(df$betweenness_norm)]
    expect_true(all(dn >= 0 & dn <= 1 + 1e-9), info = sprintf("n=%d m=%d degree_norm range", n, m))
    expect_true(all(bn >= 0 & bn <= 1 + 1e-9), info = sprintf("n=%d m=%d betweenness_norm range", n, m))
  }
})

test_that("complete bipartite K_{n,m}: every degree_norm == 1, betweenness_norm equal within a mode", {
  n <- 4; m <- 6
  comp <- paste0("C", seq_len(n)); targ <- paste0("T", seq_len(m))
  grid <- expand.grid(c = comp, t = targ, stringsAsFactors = FALSE)
  g <- .bip_graph(comp, targ, grid)

  df <- patliR:::.network_centrality_one(g, "x", c("degree", "betweenness"), TRUE)
  expect_equal(df$degree_norm, rep(1, n + m))
  bn_c <- df$betweenness_norm[df$node_type == "compound"]
  bn_t <- df$betweenness_norm[df$node_type == "target"]
  expect_equal(length(unique(round(bn_c, 8))), 1L)
  expect_equal(length(unique(round(bn_t, 8))), 1L)
})

test_that("bipartite star: degree_norm(compound) == degree_norm(target) == 1, betweenness_norm(compound) == 1", {
  m <- 5
  comp <- "C1"; targ <- paste0("T", seq_len(m))
  g <- .bip_graph(comp, targ, data.frame(c = comp, t = targ, stringsAsFactors = FALSE))

  df <- patliR:::.network_centrality_one(g, "x", c("degree", "betweenness"), TRUE)
  expect_equal(df$degree_norm[df$node_type == "compound"], 1)
  expect_equal(df$degree_norm[df$node_type == "target"], rep(1, m))
  ## compound betweenness = C(5,2) = 10; B_max(U) for n=1,m=5 = 5*4/2 = 10
  expect_equal(df$betweenness_norm[df$node_type == "compound"], 1)
})

test_that("hand-computed B_max for n = 3, m = 5 matches the formula", {
  ## B_max(U) (compound mode, n_mode = 3, n_other = 5):
  ##   s = (3-1) %/% 5 = 0 ; t = (3-1) %% 5 = 2
  ##   0.5 * (5^2 * 1^2 + 5 * 1 * (2*2 - 0 - 1) - 2 * (2*0 - 2 + 3))
  ##   = 0.5 * (25 + 15 - 2) = 0.5 * 38 = 19
  expect_equal(patliR:::.network_bipartite_bmax(3, 5), 19)
  ## B_max(V) (target mode, n_mode = 5, n_other = 3):
  ##   p = (5-1) %/% 3 = 1 ; r = (5-1) %% 3 = 1
  ##   0.5 * (3^2 * 2^2 + 3 * 2 * (2*1 - 1 - 1) - 1 * (2*1 - 1 + 3))
  ##   = 0.5 * (36 + 0 - 4) = 0.5 * 32 = 16
  expect_equal(patliR:::.network_bipartite_bmax(5, 3), 16)
})

test_that("B_max is NA for an empty opposite mode; degree_norm is NA (not NaN/0) when a mode is empty", {
  expect_true(is.na(patliR:::.network_bipartite_bmax(3, 0)))
  expect_true(is.na(patliR:::.network_bipartite_bmax(0, 3)))

  ## a lone compound, no targets, no edges -> m == 0
  g <- igraph::graph_from_data_frame(
    data.frame(from = character(0), to = character(0)),
    directed = FALSE,
    vertices = data.frame(name = "C1", type = FALSE, stringsAsFactors = FALSE)
  )
  df <- patliR:::.network_centrality_one(g, "x", c("degree", "betweenness"), TRUE)
  expect_true(is.na(df$degree_norm))
  expect_true(is.na(df$betweenness_norm))
  expect_false(is.nan(df$degree_norm))
})

test_that("hub_score_component: max 1 within each edge-bearing component, NA for an isolated node", {
  ## three disjoint compound-target edges + one isolated compound = 4 components
  g <- .bip_graph(
    c("C1", "C2", "C3", "C4"), c("T1", "T2", "T3"),
    data.frame(c = c("C1", "C2", "C3"), t = c("T1", "T2", "T3"), stringsAsFactors = FALSE)
  )
  df <- patliR:::.network_centrality_one(g, "x", "hub_score", TRUE)
  expect_equal(length(unique(df$component_id)), 4L)
  ## isolated compound C4 -> NA (eigenvector centrality is undefined with no edges)
  expect_true(is.na(df$hub_score_component[df$node_id == "C4"]))
  non_iso <- df[df$node_id != "C4", ]
  comp_max <- tapply(non_iso$hub_score_component, non_iso$component_id, max)
  expect_true(all(abs(comp_max - 1) < 1e-9))
})

test_that("hub_score / hub_score_component are deterministic across RNG seeds (B1)", {
  ## igraph::hits_scores()$hub on a bipartite graph drew an RNG-seeded
  ## arbitrary vector from a degenerate eigenspace, so two runs of the same
  ## pipeline under different seeds gave hub scores differing by the full
  ## 0..1 range. eigen_centrality() (Perron vector) is unique per component.
  proj <- .network_stats_test_setup()
  set.seed(1L)
  r1 <- patliRResults(network_centrality(proj, condition = "FLO-ET"), "network_centrality")
  set.seed(20260903L)
  r2 <- patliRResults(network_centrality(proj, condition = "FLO-ET"), "network_centrality")
  ## ARPACK contributes ~1e-15 float noise, never the 0<->1 swing of the bug
  expect_equal(r1$hub_score, r2$hub_score, tolerance = 1e-8)
  expect_equal(r1$hub_score_component, r2$hub_score_component, tolerance = 1e-8)
  expect_identical(round(r1$hub_score, 8), round(r2$hub_score, 8))
})

test_that("hub_score_component is strictly positive for both modes of a connected bipartite graph (B1 non-degeneracy)", {
  n <- 3; m <- 4
  comp <- paste0("C", seq_len(n)); targ <- paste0("T", seq_len(m))
  g <- .bip_graph(comp, targ, expand.grid(c = comp, t = targ, stringsAsFactors = FALSE))
  df <- patliR:::.network_centrality_one(g, "x", "hub_score", TRUE)
  ## the old HITS call put one entire mode at exactly 0 on some seeds
  expect_true(all(df$hub_score_component > 0))
})

test_that("network_centrality(normalize = FALSE) keeps the raw values and leaves every normalisation column NA", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, condition = "FLO-ET", normalize = FALSE)
  result <- patliRResults(proj, "network_centrality")
  ## schema does not depend on `normalize`
  expect_true(all(c(
    "condition", "node_id", "node_type", "degree", "betweenness", "hub_score",
    "degree_norm", "betweenness_norm", "hub_score_component", "component_id",
    "n_compounds", "n_targets"
  ) %in% names(result)))
  ## raw columns unchanged
  expect_true(all(is.finite(result$degree)))
  ## normalisation columns all NA
  for (col in c("degree_norm", "betweenness_norm", "hub_score_component",
                "component_id", "n_compounds", "n_targets")) {
    expect_true(all(is.na(result[[col]])), info = col)
  }
})

test_that("network_centrality(normalize = FALSE) upserts cleanly onto a normalised table (S3)", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)                       # normalize = TRUE, all conditions
  expect_no_error(
    proj <- network_centrality(proj, condition = "FLO-ET", normalize = FALSE)
  )
  result <- patliRResults(proj, "network_centrality")
  flo <- result[result$condition == "FLO-ET", ]
  expect_true(all(is.na(flo$degree_norm)))
  ## the other conditions keep their normalised values
  other <- result[result$condition != "FLO-ET", ]
  expect_true(any(!is.na(other$degree_norm)))
})

test_that("n_compounds + n_targets == vcount(g) for every condition", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  result <- patliRResults(proj, "network_centrality")
  for (cond in unique(result$condition)) {
    g <- patliR:::.network_graph(proj, cond)
    row <- result[result$condition == cond, ][1, ]
    expect_equal(row$n_compounds + row$n_targets, igraph::vcount(g))
  }
})

test_that("the empty-row constructor declares the normalisation columns (zero-row schema invariant)", {
  emp <- patliR:::.empty_network_centrality_row()
  expect_true(all(c(
    "degree_norm", "betweenness_norm", "hub_score_component",
    "component_id", "n_compounds", "n_targets"
  ) %in% names(emp)))
  expect_equal(nrow(emp), 0L)
})

test_that("network_centrality() rebuilding one condition does not touch the others", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  before <- patliRResults(proj, "network_centrality")

  proj <- network_centrality(proj, condition = "FLO-ET")
  after <- patliRResults(proj, "network_centrality")

  expect_setequal(unique(after$condition), unique(before$condition))
})
