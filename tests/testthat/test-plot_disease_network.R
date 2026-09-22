## targets_disease_profile()'s own live-network calls are not mocked here
## (see test-targets_disease.R) -- plot_disease_network() only ever READS an
## already-populated targets_disease_profile table, so these tests build
## that table by hand, fully offline, and focus on the plotting logic
## itself (in particular the convex-hull branch, which needs >= 3 targets
## sharing one disease to engage -- see .chemical_space_ggplot()'s identical
## chull()-skip-below-3-points convention this reuses).

.disease_network_fake_profile <- function(compound_id, target_ids, condition_targets,
                                           disease_id = "D1", disease_name = "Test Disease") {
  data.frame(
    compound_id = compound_id, target_id = target_ids,
    disease_id = disease_id, disease_name = disease_name,
    association_score = seq(0.9, 0.5, length.out = length(target_ids))[seq_along(target_ids)],
    evidence = NA_character_, rank = 1L, mode = "explore",
    stringsAsFactors = FALSE
  )
}

test_that("plot_disease_network() errors clearly without targets_disease_profile", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_disease_network(proj, condition = "FLO-ET", save = FALSE), "targets_disease_profile")
})

test_that("plot_disease_network() draws a convex hull when a disease has >= 3 targets, and skips it otherwise", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  targets <- unique(flo$uniprot_id)
  testthat::skip_if(length(targets) < 3, "fixture needs at least 3 FLO-ET targets")
  cmp <- flo$compound_id[1]

  ## disease D1 gets 3 targets (hull should be drawn), disease D2 gets 1 (no hull)
  profile <- rbind(
    .disease_network_fake_profile(cmp, targets[1:3], disease_id = "D1", disease_name = "Hull Disease"),
    .disease_network_fake_profile(cmp, targets[1], disease_id = "D2", disease_name = "Single-target Disease")
  )
  patliRResults(proj, "targets_disease_profile") <- profile

  p <- plot_disease_network(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)
  ## first layer is the hull polygon (conditionally added) or the edge
  ## segments if no hull was drawn -- check a GeomPolygon layer exists with
  ## exactly one group (D1), not two.
  poly_layers <- vapply(built$plot$layers, function(l) inherits(l$geom, "GeomPolygon"), logical(1))
  expect_true(any(poly_layers))
  hull_data <- built$data[[which(poly_layers)[1]]]
  expect_equal(length(unique(hull_data$group)), 1)
})

test_that("plot_disease_network() disease nodes are drawn larger than compound/target nodes", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  cmp <- flo$compound_id[1]
  target <- flo$uniprot_id[1]

  patliRResults(proj, "targets_disease_profile") <- .disease_network_fake_profile(cmp, target)

  p <- plot_disease_network(proj, condition = "FLO-ET", save = FALSE)
  built <- ggplot2::ggplot_build(p)
  point_layers <- vapply(built$plot$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
  pts <- built$data[[which(point_layers)[1]]]
  ## disease point_size (9) must exceed both compound (3.2) and target (2.6)
  expect_equal(max(pts$size), 9)
  expect_true(min(pts$size) < 9)
})

test_that("plot_disease_network() disease/max_rank filters restrict the plotted table", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  cmp <- flo$compound_id[1]
  target <- flo$uniprot_id[1]

  profile <- rbind(
    data.frame(compound_id = cmp, target_id = target, disease_id = "D1", disease_name = "Keep",
               association_score = 0.9, evidence = NA_character_, rank = 1L, mode = "explore", stringsAsFactors = FALSE),
    data.frame(compound_id = cmp, target_id = target, disease_id = "D2", disease_name = "Drop",
               association_score = 0.5, evidence = NA_character_, rank = 2L, mode = "explore", stringsAsFactors = FALSE)
  )
  patliRResults(proj, "targets_disease_profile") <- profile

  p1 <- plot_disease_network(proj, condition = "FLO-ET", max_rank = 1, save = FALSE)
  expect_s3_class(p1, "ggplot")

  p2 <- plot_disease_network(proj, condition = "FLO-ET", disease = "D1", save = FALSE)
  expect_s3_class(p2, "ggplot")

  expect_error(
    plot_disease_network(proj, condition = "FLO-ET", disease = "not_a_real_disease", save = FALSE),
    "targets_disease_profile"
  )
})

test_that("plot_disease_network() saves a PNG and logs it", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  patliRResults(proj, "targets_disease_profile") <- .disease_network_fake_profile(flo$compound_id[1], flo$uniprot_id[1])

  p <- plot_disease_network(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "disease_network_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_match(log_df$path[1], "\\.png$")
})
