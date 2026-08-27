test_that("plot_synergy() returns a ggplot and errors clearly without network_synergy()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_synergy(proj, condition = "FLO-ET", save = FALSE), "network_synergy")

  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", "compound_id"])
  testthat::skip_if(length(ct) < 2, "need at least 2 compounds in this fixture's FLO-ET condition")
  fake_syn <- data.frame(
    condition = "FLO-ET", disease_id = "SOME_DISEASE", compound_a = ct[1], compound_b = ct[2],
    target_jaccard = 0.2, complementarity = 0.8, z_score_a = -1.5, z_score_b = -1.2,
    joint_closeness = 1.35, synergy_score = 1.08, pairs_mode = "all", stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_synergy") <- fake_syn

  p <- plot_synergy(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})
