test_that("plot_centrality() returns a ggplot and errors clearly without network_centrality()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_centrality(proj, condition = "FLO-ET", save = FALSE), "network_centrality")

  proj <- network_centrality(proj)
  p <- plot_centrality(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_centrality() errors clearly on a measure not present in network_centrality()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, measures = "degree")
  expect_error(plot_centrality(proj, condition = "FLO-ET", measure = "betweenness", save = FALSE), "betweenness")
})

test_that("plot_centrality() defaults to degree_norm when node_type = 'both', informing why", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  expect_message(
    p <- plot_centrality(proj, condition = "FLO-ET", save = FALSE),
    "degree_norm"
  )
  expect_true("degree_norm" %in% names(p$data))
  expect_match(as.character(p$labels$y), "\\|T\\|")
})

test_that("plot_centrality() falls back to raw degree with a warning when degree_norm has no values (normalize = FALSE)", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, normalize = FALSE)
  expect_warning(
    p <- plot_centrality(proj, condition = "FLO-ET", save = FALSE),
    "normalize = FALSE"
  )
  ## the degree_norm column exists but is all NA, so the plot uses raw degree
  expect_true(all(is.na(p$data$degree_norm)))
  expect_identical(as.character(p$labels$y), "degree")
})

test_that("plot_centrality() saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  p <- plot_centrality(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "centrality_plot_log")
  expect_true(file.exists(log_df$path[1]))
})

test_that("centrality scope keys separate pooled scopes from literal condition names", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  conditions <- c("A", "B", "A+B", "ALL")
  patliRResults(proj, "network_edges") <- data.frame(condition = conditions, compound_id = "c1", uniprot_id = "t1", weight = 1)
  patliRResults(proj, "network_centrality") <- data.frame(condition = conditions, node_id = "c1", node_type = "compound", degree = 1)
  scope <- patliR:::.plot_scope
  expect_identical(scope(proj, "A")$scope_label, "A")
  expect_identical(scope(proj, "A+B")$scope_label, "A+B")
  expect_identical(scope(proj, NULL)$scope_label, "ALL")
  expect_identical(scope(proj, c("A", "B"))$scope_label, scope(proj, c("B", "A"))$scope_label)
  expect_false(identical(scope(proj, c("A", "B"))$scope_label, scope(proj, "A+B")$scope_label))
  expect_false(identical(scope(proj, "ALL")$scope_label, scope(proj, NULL)$scope_label))
  for (condition in list("A+B", c("A", "B"), "ALL", NULL)) {
    p <- plot_centrality(proj, condition = condition, measure = "degree", width = 3, height = 3)
    proj <- attr(p, "proj")
  }
  log <- patliRResults(proj, "centrality_plot_log")
  expect_equal(nrow(log), 4L)
  expect_equal(length(unique(log$path)), 4L)
  expect_true(all(file.exists(log$path)))
  p <- plot_centrality(proj, condition = c("B", "A"), measure = "degree", width = 3, height = 3)
  expect_equal(nrow(patliRResults(attr(p, "proj"), "centrality_plot_log")), 4L)
})
