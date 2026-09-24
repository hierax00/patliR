## Fabricate a network_synergy() row directly -- these tests do not need
## STRINGdb or a real network_synergy() run, only the post-Phase-3-piece-12
## schema (network_synergy.R's @return / .empty_network_synergy_row()).
.plot_synergy_fake_row <- function(a, b, cheng_class, condition = "FLO-ET",
                                    disease_id = "D1", synergy_score = NA_real_) {
  sep <- if (is.na(cheng_class)) NA else cheng_class %in% c("P2", "P4", "P6")
  n_prox <- if (is.na(cheng_class)) NA else switch(
    cheng_class, P1 = 2L, P2 = 2L, P3 = 1L, P4 = 1L, P5 = 0L, P6 = 0L
  )
  proximal_a <- if (is.na(cheng_class)) NA else n_prox >= 1L
  proximal_b <- if (is.na(cheng_class)) NA else n_prox >= 1L
  za <- if (isTRUE(proximal_a)) -2 else 1
  zb <- if (isTRUE(proximal_b)) -1.5 else 0.8
  data.frame(
    condition = condition, disease_id = disease_id, compound_a = a, compound_b = b,
    n_targets_a = 3L, n_targets_b = 3L, n_targets_a_mapped = 3L, n_targets_b_mapped = 3L,
    target_jaccard = 0.2, complementarity = 0.8, singleton_a = FALSE, singleton_b = FALSE,
    d_aa = 1, d_bb = 1, d_ab = 2,
    s_ab = if (is.na(cheng_class)) NA_real_ else ifelse(sep, 1.5, -0.5),
    separated = sep, z_score_a = za, z_score_b = zb, p_adjusted_a = 0.01, p_adjusted_b = 0.01,
    proximal_a = proximal_a, proximal_b = proximal_b,
    both_proximal = isTRUE(proximal_a) && isTRUE(proximal_b),
    cheng_class = cheng_class,
    complementary_exposure = if (is.na(cheng_class)) NA else cheng_class == "P2",
    joint_closeness = 1.5, synergy_score = synergy_score,
    separation_method = "network", alpha = 0.05, pairs_mode = "all",
    species = 9606, string_version = "12.0", disease_gene_source = "disease_genes",
    stringsAsFactors = FALSE
  )
}

test_that("plot_synergy() errors clearly without network_synergy()", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_synergy(proj, condition = "FLO-ET", save = FALSE), "network_synergy")
})

test_that("plot_synergy() draws the Cheng quadrant over all P1..P6 + NA classes, and the P2 annotation count matches a hand count", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  skip_if(length(ct) < 4, "need at least 4 FLO-ET compounds for 6 distinct pairs")

  pairs <- utils::combn(ct[1:4], 2, simplify = FALSE)
  classes <- c("P1", "P2", "P3", "P4", "P5", "P6")
  fake <- do.call(rbind, lapply(seq_along(pairs), function(i) {
    cc <- classes[i]
    score <- if (cc == "P2") 0.5 else NA_real_
    .plot_synergy_fake_row(pairs[[i]][1], pairs[[i]][2], cc, synergy_score = score)
  }))
  ## a second, unclassified pair -- separated/proximal/cheng_class all NA
  ## (e.g. separation = "jaccard" for this one row, or an unmapped pair)
  fake <- rbind(fake, .plot_synergy_fake_row(ct[1], ct[2], NA_character_, disease_id = "D2"))
  ## separated TRUE/FALSE/NA all present across `fake` by construction
  expect_true(all(c(TRUE, FALSE) %in% fake$separated))
  expect_true(anyNA(fake$separated))

  patliRResults(proj, "network_synergy") <- fake

  p <- expect_no_error(plot_synergy(proj, condition = "FLO-ET", save = FALSE))
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)
  ann_labels <- unlist(lapply(built$data, function(d) {
    if ("label" %in% names(d)) d$label[grepl("Complementary Exposure|P2 pair", d$label)] else NULL
  }))
  expect_true(any(grepl("Complementary Exposure", ann_labels)))

  hand_count_d1 <- sum(fake$cheng_class == "P2" & fake$disease_id == "D1", na.rm = TRUE)
  hand_count_d2 <- sum(fake$cheng_class == "P2" & fake$disease_id == "D2", na.rm = TRUE)
  expect_true(any(grepl(paste0("n = ", hand_count_d1, "\\b"), ann_labels)))
  expect_true(any(grepl(paste0("n = ", hand_count_d2, "\\b"), ann_labels)))
  expect_equal(hand_count_d1, 1L)
  expect_equal(hand_count_d2, 0L)
})

test_that("plot_synergy() draws a self-explanatory empty panel when no row has a Cheng classification", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  skip_if(length(ct) < 2, "need at least 2 compounds")

  fake_jaccard <- .plot_synergy_fake_row(ct[1], ct[2], NA_character_)
  fake_jaccard$separation_method <- "jaccard"
  patliRResults(proj, "network_synergy") <- fake_jaccard

  p <- expect_no_error(plot_synergy(proj, condition = "FLO-ET", save = FALSE))
  expect_s3_class(p, "ggplot")
  ## the empty panel is built with annotate("text", ...), not geom_point
  expect_true(length(p$layers) >= 1)
})

test_that("plot_synergy() handles a legacy table missing the Menche/Cheng columns entirely (defensive back-fill)", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  skip_if(length(ct) < 2, "need at least 2 compounds")

  legacy <- data.frame(
    condition = "FLO-ET", disease_id = "D1", compound_a = ct[1], compound_b = ct[2],
    target_jaccard = 0.2, complementarity = 0.8, z_score_a = -1.5, z_score_b = -1.2,
    joint_closeness = 1.35, synergy_score = 1.08, pairs_mode = "all", stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_synergy") <- legacy
  ## no cheng_class column at all -> n_class == 0 -> empty panel, no error
  p <- expect_no_error(plot_synergy(proj, condition = "FLO-ET", save = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("plot_synergy() saves both the main scatter and the s_AB companion PNG, logging both", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  skip_if(length(ct) < 2, "need at least 2 compounds")

  patliRResults(proj, "network_synergy") <- .plot_synergy_fake_row(ct[1], ct[2], "P2", synergy_score = 0.5)
  p <- plot_synergy(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")

  main_log <- patliRResults(proj2, "synergy_plot_log")
  sab_log <- patliRResults(proj2, "synergy_sab_plot_log")
  expect_true(file.exists(main_log$path[1]))
  expect_true(file.exists(sab_log$path[1]))
  expect_false(identical(main_log$path[1], sab_log$path[1]))
})

test_that("plot_synergy() top_n default is 5 (spec breaking-change from 10)", {
  expect_equal(formals(plot_synergy)$top_n, 5)
})
test_that("unclassified synergy pairs still contribute to the separation histogram", {
  skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "t1")
  patliRResults(proj, "network_synergy") <- data.frame(
    condition = "A", disease_id = "d1", cheng_class = NA_character_, s_ab = c(-0.5, 1, NA)
  )
  testthat::local_mocked_bindings(
    .synergy_finish_both = function(proj, p_main, p_sab, ...) list(main = p_main, sab = p_sab),
    .package = "patliR"
  )
  plots <- plot_synergy(proj, condition = "A", save = FALSE)
  expect_equal(plots$sab$data$s_ab, c(-0.5, 1))
  expect_equal(sum(ggplot2::ggplot_build(plots$sab)$data[[1]]$count), 2)
  expect_true(inherits(plots$main$layers[[1]]$geom, "GeomText"))
})

test_that("plot_synergy() draws the P2 quadrant name and its count as separate single-line labels", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "t1")
  patliRResults(proj, "network_synergy") <- data.frame(
    condition = "A", disease_id = "d1", compound_a = c("c1", "c2"), compound_b = c("c2", "c3"),
    z_score_a = c(-3, -1), z_score_b = c(-2, -4), s_ab = c(0.5, -0.2), separated = c(TRUE, FALSE),
    cheng_class = c("P2", "P1"), synergy_score = c(1, NA), stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .synergy_finish_both = function(proj, p_main, p_sab, ...) list(main = p_main, sab = p_sab),
    .package = "patliR"
  )
  plots <- plot_synergy(proj, condition = "A", save = FALSE)
  txt <- Filter(function(l) inherits(l$geom, "GeomText"), plots$main$layers)
  labels <- unlist(lapply(txt, function(l) c(l$data$quadrant_label, l$data$count_label)))
  expect_true(any(grepl("Complementary Exposure (P2)", labels, fixed = TRUE)))
  expect_true(any(grepl("n = 1 P2 pair", labels, fixed = TRUE)))
  expect_false(any(grepl("\n", labels, fixed = TRUE)))
})
