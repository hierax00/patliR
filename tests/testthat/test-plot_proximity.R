test_that("plot_proximity() returns a ggplot and errors clearly without network_proximity()", {
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_proximity(proj, condition = "FLO-ET", save = FALSE), "network_proximity")

  ## Fabricate network_proximity results directly -- bypassing the real
  ## STRINGdb-dependent computation -- same pattern test-network_synergy.R
  ## already uses for the same reason (offline, self-contained test).
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  compounds <- unique(ct$compound_id)
  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = compounds, disease_id = "SOME_DISEASE",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, d_observed = 2,
    d_random_mean = 3, d_random_sd = 1, z_score = seq(-2, 2, length.out = length(compounds)),
    n_random = 100L, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  p <- plot_proximity(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})
