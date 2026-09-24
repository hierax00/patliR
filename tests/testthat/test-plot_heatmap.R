test_that("plot_heatmap() returns a pheatmap object for both 'what' modes", {
  testthat::skip_if_not_installed("pheatmap")

  proj <- .network_stats_test_setup()
  p1 <- plot_heatmap(proj, what = "compound_condition", save = FALSE)
  expect_s3_class(p1, "pheatmap")

  p2 <- plot_heatmap(proj, condition = "FLO-ET", what = "compound_target", save = FALSE)
  expect_s3_class(p2, "pheatmap")
})

test_that("plot_heatmap() saves a PDF and logs it", {
  testthat::skip_if_not_installed("pheatmap")

  proj <- .network_stats_test_setup()
  p <- plot_heatmap(proj, what = "compound_condition", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "heatmap_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.pdf$")
})

test_that("plot_heatmap() keeps duplicate labels and unknown edge weights separate", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .test_project()
  edges <- data.frame(condition = "A", compound_id = c("c1", "c1", "c2", "c2"),
                      uniprot_id = c("t1", "t1", "t1", "t2"), weight = c(0.2, 0.8, NA, 0.4))
  patliRResults(proj, "network_edges") <- edges
  testthat::local_mocked_bindings(
    .plot_label_nodes = function(proj, conditions, ids, type) rep("same", length(ids)),
    .plot_pheatmap_render = function(proj, mat, ...) mat,
    .package = "patliR"
  )
  mat <- plot_heatmap(proj, condition = "A", save = FALSE)
  expect_equal(dim(mat), c(2L, 2L))
  expect_equal(rownames(mat), c("same (c1)", "same (c2)"))
  expect_equal(colnames(mat), c("same (t1)", "same (t2)"))
  expect_equal(unname(mat[1, ]), c(0.8, 0))
  expect_true(is.na(mat[2, 1]))
  expect_equal(mat[2, 2], 0.4)
  expect_error(plot_heatmap(proj, condition = character(), save = FALSE), "Nothing to plot")
})

test_that("plot_heatmap() renders an entirely unknown matrix without clustering errors", {
  testthat::skip_if_not_installed("pheatmap")
  proj <- .test_project()
  edges <- expand.grid(compound_id = c("c1", "c2"), uniprot_id = c("t1", "t2"), stringsAsFactors = FALSE)
  edges$condition <- "A"
  edges$weight <- NA_real_
  patliRResults(proj, "network_edges") <- edges
  mat <- patliR:::.plot_edge_matrix(proj, "A", edges)
  expect_true(all(is.na(mat)))
  p <- plot_heatmap(proj, condition = "A", save = TRUE)
  expect_s3_class(p, "pheatmap")
  expect_identical(p$tree_row, NA)
  expect_identical(p$tree_col, NA)
  expect_true(file.exists(patliRResults(attr(p, "proj"), "heatmap_plot_log")$path))
})
