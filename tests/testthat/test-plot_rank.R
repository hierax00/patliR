.plot_rank_test_setup <- function() {
  testthat::skip_if_not_installed("RobustRankAggreg")
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  rank_candidates(proj, condition = "FLO-ET")
}

test_that("plot_rank() chooses a finite default axis within the current scope", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = c("A", "B"), compound_id = "c1", uniprot_id = "t1", weight = 1)
  rc <- data.frame(condition = c("A", "B"), compound_id = "c1", rra_rank = 1L,
                   pareto_front = 1L, crit_synergy_best = c(NA, 1), crit_proximity_z = c(Inf, 2),
                   crit_centrality = c(0.5, 0.8))
  patliRResults(proj, "rank_candidates") <- rc
  p <- plot_rank(proj, condition = "A", save = FALSE)
  expect_identical(rlang::as_label(p$mapping$y), "crit_centrality")
  rc$crit_centrality[1] <- NaN
  patliRResults(proj, "rank_candidates") <- rc
  expect_error(plot_rank(proj, condition = "A", save = FALSE), "No usable numeric column")
})

test_that("rank heatmaps save distinct files for each target relevance", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "t1", weight = 1)
  patliRResults(proj, "rank_candidates") <- data.frame(condition = "A", compound_id = "c1", rra_rank = 1L)
  patliRResults(proj, "network_centrality") <- data.frame(condition = "A", node_id = "t1", node_type = "target", degree_norm = 1)
  for (relevance in c("breadth", "weight", "centrality")) {
    p <- plot_rank(proj, condition = "A", view = "heatmap", target_relevance = relevance)
    proj <- attr(p, "proj")
  }
  log <- patliRResults(proj, "rank_heatmap_plot_log")
  expect_equal(nrow(log), 3L)
  expect_equal(length(unique(log$path)), 3L)
  expect_true(all(file.exists(log$path)))
  expect_true(all(vapply(seq_len(nrow(log)), function(i) grepl(log$target_relevance[i], basename(log$path[i]), fixed = TRUE), logical(1))))
})

test_that("rank heatmaps retain IDs and NA weights and select only connected centrality targets", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = c("c1", "c2"), uniprot_id = c("t1", "t2"), weight = NA_real_)
  patliRResults(proj, "rank_candidates") <- data.frame(condition = "A", compound_id = c("c1", "c2"), rra_rank = 1:2)
  cent <- data.frame(condition = "A", node_id = c("outside", "t1", "t2"), node_type = "target", degree_norm = c(1, NA, NA))
  patliRResults(proj, "network_centrality") <- cent
  testthat::local_mocked_bindings(
    .plot_label_nodes = function(proj, conditions, ids, type) rep("same", length(ids)),
    .plot_pheatmap_render = function(proj, mat, ...) mat,
    .package = "patliR"
  )
  for (relevance in c("breadth", "weight", "centrality")) {
    mat <- plot_rank(proj, condition = "A", view = "heatmap", target_relevance = relevance, save = FALSE)
    expect_equal(dim(mat), c(2L, 2L))
    expect_equal(colnames(mat), c("same (t1)", "same (t2)"))
    expect_true(all(is.na(diag(mat))))
  }
  mat <- plot_rank(proj, condition = "A", view = "heatmap", target_relevance = "centrality", top_n_targets = 1, save = FALSE)
  expect_equal(ncol(mat), 1L)
  cent$degree_norm <- NULL
  patliRResults(proj, "network_centrality") <- cent
  expect_no_error(plot_rank(proj, condition = "A", view = "heatmap", target_relevance = "centrality", save = FALSE))
  cent <- cent[cent$node_id == "outside", , drop = FALSE]
  patliRResults(proj, "network_centrality") <- cent
  expect_error(plot_rank(proj, condition = "A", view = "heatmap", target_relevance = "centrality", save = FALSE), "No centrality rows")
  rc <- patliRResults(proj, "rank_candidates")
  rc$compound_id <- "absent"
  patliRResults(proj, "rank_candidates") <- rc
  expect_error(plot_rank(proj, condition = "A", view = "heatmap", save = FALSE), "No.*network_edges")
})

test_that("plot_rank() errors clearly without rank_candidates()", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_rank(proj, condition = "FLO-ET", save = FALSE), "rank_candidates")
})

test_that("plot_rank(view = 'pareto') returns a ggplot with the documented default axes", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", save = FALSE)
  expect_s3_class(p, "ggplot")
  ## no optional criteria ran -> y falls back to crit_centrality (the last
  ## entry in .plot_rank_y_priority still present in the mandatory-only table)
  expect_identical(rlang::as_label(p$mapping$x), "rra_rank")
  expect_identical(rlang::as_label(p$mapping$y), "crit_centrality")
})

test_that("plot_rank(view = 'pareto') honours explicit x/y overrides", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", x = "crit_adme_pass_frac", y = "crit_centrality", save = FALSE)
  expect_identical(rlang::as_label(p$mapping$x), "crit_adme_pass_frac")
  expect_identical(rlang::as_label(p$mapping$y), "crit_centrality")
})

test_that("plot_rank(view = 'pareto') rejects an unknown column name", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  expect_error(plot_rank(proj, condition = "FLO-ET", view = "pareto", x = "not_a_column", save = FALSE), "not_a_column")
})

test_that("plot_rank(view = 'pareto') saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .plot_rank_test_setup()
  p <- plot_rank(proj, condition = "FLO-ET", view = "pareto", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "rank_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.png$")
})

test_that("plot_rank(view = 'heatmap') returns a pheatmap object restricted to top_n_compounds/top_n_targets", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .plot_rank_test_setup()
  n_all_compounds <- length(unique(patliRResults(proj, "rank_candidates")$compound_id))
  testthat::skip_if(n_all_compounds < 2, "need at least 2 compounds")

  h <- plot_rank(proj, condition = "FLO-ET", view = "heatmap", top_n_compounds = 2, top_n_targets = 3, save = FALSE)
  expect_s3_class(h, "pheatmap")
  expect_lte(nrow(h$gtable), 100) # sanity: object renders without error
})

test_that("plot_rank(view = 'heatmap') saves a PDF and logs it", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .plot_rank_test_setup()
  h <- plot_rank(proj, condition = "FLO-ET", view = "heatmap", top_n_compounds = 2, top_n_targets = 3, save = TRUE)
  proj2 <- attr(h, "proj")
  log_df <- patliRResults(proj2, "rank_heatmap_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.pdf$")
})

test_that("plot_rank(view = 'heatmap', target_relevance = 'centrality') requires network_centrality", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .network_stats_test_setup()
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  ## build rank_candidates with only adme+centrality is required, so
  ## network_centrality must exist for rank_candidates() itself already --
  ## exercise the heatmap-specific centrality lookup by requesting a
  ## disease-agnostic condition with an intentionally cleared centrality slot.
  proj <- network_centrality(proj)
  proj <- rank_candidates(proj, condition = "FLO-ET")
  patliRResults(proj, "network_centrality") <- NULL
  expect_error(
    plot_rank(proj, condition = "FLO-ET", view = "heatmap", target_relevance = "centrality", save = FALSE),
    "network_centrality"
  )
})

test_that("plot_rank(view = 'pareto') draws shortened front-1 labels with x headroom", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  long_id <- "Naphthalene, 1,2-dihydro-1,1,6-trimethyl-"
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = c(long_id, "c2"), uniprot_id = "t1", weight = 1)
  patliRResults(proj, "rank_candidates") <- data.frame(
    condition = "A", compound_id = c(long_id, "c2"), rra_rank = 1:2, pareto_front = c(1L, 2L),
    crit_centrality = c(0.9, 0.1), stringsAsFactors = FALSE
  )
  p <- plot_rank(proj, condition = "A", save = FALSE)
  lab <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p$layers)[[1]]$data
  expect_equal(nrow(lab), 1)
  expect_true(nchar(lab$short_label) <= 28)
  expect_identical(lab$compound_label, long_id)
  expect_equal(p$scales$get_scales("x")$expand, ggplot2::expansion(mult = 0.12))
})
