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

test_that(".venn_shared_targets() ranks shared targets by the number of compounds hitting them", {
  edges <- data.frame(compound_id = c("c1", "c2", "c3", "c1", "c2", "c1", "c1"),
                      uniprot_id = c("T1", "T1", "T1", "T2", "T2", "T3", "T3"), stringsAsFactors = FALSE)
  out <- patliR:::.venn_shared_targets(edges, c("T3", "T2", "T1", "T9"))
  expect_identical(out$target_id, c("T1", "T2", "T3", "T9"))
  expect_identical(out$n_compounds, c(3L, 2L, 1L, 0L))
  expect_equal(nrow(patliR:::.venn_shared_targets(edges, character(0))), 0L)
})

test_that("plot_venn() names the top shared targets and the disease, and validates label_top", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggVennDiagram")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  patliRResults(proj, "targets_disease") <- data.frame(
    compound_id = ct$compound_id, target_id = ct$uniprot_id, disease_id = "MONDO_X",
    association_score = 0.5, stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(disease_id = "MONDO_X", disease_name = "test disorder",
                                                     uniprot_id = "P1", stringsAsFactors = FALSE)
  expect_error(plot_venn(proj, condition = "FLO-ET", label_top = -1, save = FALSE), "label_top")
  expect_error(plot_venn(proj, condition = "FLO-ET", label_top = NA, save = FALSE), "label_top")

  p <- plot_venn(proj, condition = "FLO-ET", label_top = 3, save = FALSE)
  shared <- attr(p, "shared_targets")
  expect_setequal(shared$target_id, unique(ct$uniprot_id))
  expect_false(is.unsorted(rev(shared$n_compounds)))
  skip_if_not_installed("patchwork")
  expect_s3_class(p, "patchwork")
  expect_match(p$patches$annotation$title, "test disorder (MONDO_X)", fixed = TRUE)

  bare <- plot_venn(proj, condition = "FLO-ET", label_top = 0, save = FALSE)
  expect_false(inherits(bare, "patchwork"))
  expect_match(bare$labels$title, "test disorder", fixed = TRUE)
})
