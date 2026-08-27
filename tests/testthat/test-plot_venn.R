test_that("plot_venn() returns a ggplot and errors clearly without targets_disease()", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggVennDiagram")

  proj <- .network_stats_test_setup()
  expect_error(plot_venn(proj, condition = "FLO-ET", save = FALSE), "targets_disease")

  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  fake_disease <- data.frame(
    compound_id = ct$compound_id, target_id = ct$uniprot_id,
    disease_id = "SOME_DISEASE", association_score = 0.5,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease") <- fake_disease

  p <- plot_venn(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_venn() requires an explicit disease when more than one is present", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggVennDiagram")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  fake_disease <- data.frame(
    compound_id = ct$compound_id, target_id = ct$uniprot_id,
    disease_id = rep(c("DISEASE_A", "DISEASE_B"), length.out = nrow(ct)), association_score = 0.5,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease") <- fake_disease

  expect_error(plot_venn(proj, condition = "FLO-ET", save = FALSE), "disease")
  p <- plot_venn(proj, condition = "FLO-ET", disease = "DISEASE_A", save = FALSE)
  expect_s3_class(p, "ggplot")
})
