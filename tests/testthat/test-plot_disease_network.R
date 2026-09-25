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
  ## disease D1 gets 3 targets (hull should be drawn), disease D2 gets 1 (no hull)
  profile <- rbind(
    .disease_network_fake_profile(flo$compound_id[match(targets[1:3], flo$uniprot_id)],
                                  targets[1:3], disease_id = "D1", disease_name = "Hull Disease"),
    .disease_network_fake_profile(flo$compound_id[match(targets[1], flo$uniprot_id)],
                                  targets[1], disease_id = "D2", disease_name = "Single-target Disease")
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
test_that("disease networks only include profiles matching retained compound-target pairs", {
  skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = c("A", "A", "A", "B"), compound_id = c("c1", "c2", "c2", "c1"),
    uniprot_id = c("t1", "t2", "t1", "t3"), weight = 1
  )
  profile <- data.frame(compound_id = c("c1", "c1", "c1"),
                        target_id = c("t1", "t2", "t3"), disease_id = c("d1", "d2", "d3"),
                        disease_name = c("Keep", "Wrong pair", "Other condition"),
                        association_score = 0.9, rank = 1L)
  patliRResults(proj, "targets_disease_profile") <- profile
  p <- plot_disease_network(proj, condition = "A", save = FALSE)
  nodes <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_setequal(nodes$id, c("c1", "t1", "d1"))
  edges <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)[[1]]$data
  expect_equal(nrow(edges), 2L)
  expect_no_error(ggplot2::ggplot_build(p))

  patliRResults(proj, "targets_disease_profile") <- profile[-1, ]
  expect_error(plot_disease_network(proj, condition = "A", save = FALSE),
               "No .*targets_disease_profile.* rows")
})

test_that("plot_disease_network() labels only the top_n_labels diseases with the most targets", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  targets <- unique(flo$uniprot_id)
  testthat::skip_if(length(targets) < 3, "fixture needs at least 3 FLO-ET targets")
  cmp <- flo$compound_id[match(targets[1:3], flo$uniprot_id)]
  profile <- rbind(
    .disease_network_fake_profile(cmp, targets[1:3], disease_id = "D1", disease_name = "Three-target disease"),
    .disease_network_fake_profile(cmp[1:2], targets[1:2], disease_id = "D2", disease_name = "Two-target disease"),
    .disease_network_fake_profile(cmp[1], targets[1], disease_id = "D3", disease_name = "One-target disease")
  )
  patliRResults(proj, "targets_disease_profile") <- profile
  p <- plot_disease_network(proj, condition = "FLO-ET", top_n_labels = 2, save = FALSE)
  lab <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p$layers)[[1]]$data
  expect_setequal(lab$id[lab$layer == "disease"], c("D1", "D2"))
  p0 <- plot_disease_network(proj, condition = "FLO-ET", top_n_labels = 0, top_n_targets = 0,
                             label_compounds = FALSE, save = FALSE)
  lab0 <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p0$layers)
  expect_length(lab0, 0)
})

test_that(".disease_network_top_targets() ranks by compounds hitting the target, then association score", {
  ct <- data.frame(compound_id = c("c1", "c2", "c3", "c1", "c2", "c1", "c3"),
                   uniprot_id = c("T1", "T1", "T1", "T2", "T2", "T3", "T4"), stringsAsFactors = FALSE)
  td <- data.frame(target_id = c("T1", "T2", "T3", "T4"), association_score = c(0.1, 0.2, 0.3, 0.9),
                   stringsAsFactors = FALSE)
  expect_identical(patliR:::.disease_network_top_targets(ct, td, 10), c("T1", "T2", "T4", "T3"))
  expect_identical(patliR:::.disease_network_top_targets(ct, td, 2), c("T1", "T2"))
  expect_identical(patliR:::.disease_network_top_targets(ct, td, 0), character(0))
})

test_that("plot_disease_network() labels diseases, the top targets and the compounds", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", , drop = FALSE]
  targets <- unique(flo$uniprot_id)
  testthat::skip_if(length(targets) < 3, "fixture needs at least 3 FLO-ET targets")
  cmp <- flo$compound_id[match(targets[1:3], flo$uniprot_id)]
  patliRResults(proj, "targets_disease_profile") <- .disease_network_fake_profile(cmp, targets[1:3], disease_id = "D1",
                                                                                 disease_name = "Named disease")
  p <- plot_disease_network(proj, condition = "FLO-ET", top_n_targets = 2, save = FALSE)
  lab <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")), p$layers)[[1]]$data
  expect_equal(sum(lab$layer == "target"), 2)
  expect_setequal(lab$id[lab$layer == "compound"], unique(cmp))
  expect_identical(lab$short_label[lab$layer == "disease"], "Named disease")
  expect_match(p$labels$subtitle, "Named disease (D1)", fixed = TRUE)

  expect_error(plot_disease_network(proj, condition = "FLO-ET", top_n_targets = -1, save = FALSE), "top_n_targets")
  expect_error(plot_disease_network(proj, condition = "FLO-ET", label_compounds = NA, save = FALSE), "label_compounds")
})
