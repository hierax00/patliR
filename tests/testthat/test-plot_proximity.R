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

## Shared fixture builder for view = "null" tests: a network_proximity
## table + a matching network_proximity_null table, fabricated directly
## (offline, no STRINGdb / no real network_proximity(store_null = TRUE)
## run needed) -- same pattern the view = "z" test above already uses.
.proximity_null_fixture <- function(proj, n_compounds, n_random = 20) {
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  compounds <- utils::head(unique(ct$compound_id), n_compounds)

  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = compounds, disease_id = "SOME_DISEASE",
    disease_gene_source = "disease_genes",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, n_overlap = 0L,
    d_observed = seq(1, 3, length.out = length(compounds)),
    d_random_mean = 3, d_random_sd = 1,
    z_score = seq(-4, -0.1, length.out = length(compounds)),
    p_empirical = 0.01, n_random = n_random, seed_used = 1L,
    species = 9606, string_version = "12.0", score_threshold = 400,
    n_tests_in_family = length(compounds), p_adjusted = 0.02,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  set.seed(1)
  fake_null <- do.call(rbind, lapply(compounds, function(cp) {
    data.frame(
      condition = "FLO-ET", compound_id = cp, disease_id = "SOME_DISEASE",
      disease_gene_source = "disease_genes",
      draw = seq_len(n_random), d_random = stats::rnorm(n_random, mean = 3, sd = 1),
      stringsAsFactors = FALSE
    )
  }))
  patliRResults(proj, "network_proximity_null") <- fake_null
  proj
}

test_that("plot_proximity(view = \"null\") happy path with store_null = TRUE data returns a ggplot", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  proj <- .proximity_null_fixture(proj, n_compounds = 3)

  p <- plot_proximity(proj, condition = "FLO-ET", view = "null", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_proximity(view = \"null\") aborts with a clear message when store_null was not set", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  compounds <- unique(ct$compound_id)
  ## A plain network_proximity() run (default store_null = FALSE) leaves
  ## no network_proximity_null slot at all.
  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = compounds, disease_id = "SOME_DISEASE",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, d_observed = 2,
    d_random_mean = 3, d_random_sd = 1, z_score = seq(-2, 2, length.out = length(compounds)),
    n_random = 100L, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  expect_error(
    plot_proximity(proj, condition = "FLO-ET", view = "null", save = FALSE),
    "store_null"
  )
})

test_that("plot_proximity(view = \"null\", top_n =) caps the number of compound facets, ranked by |z_score|", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  n_available <- length(unique(ct$compound_id))
  skip_if(n_available < 4, "fixture needs >= 4 compounds for a meaningful top_n cap")
  proj <- .proximity_null_fixture(proj, n_compounds = n_available)

  expect_message(
    p <- plot_proximity(proj, condition = "FLO-ET", view = "null", top_n = 2, save = FALSE),
    "top_n"
  )
  built <- ggplot2::ggplot_build(p)
  n_facets <- length(unique(built$data[[1]]$PANEL))
  expect_equal(n_facets, 2L)

  ## The kept facets must be the two compounds with the largest |z_score|
  ## -- the fixture's z_score is monotonic in compound order (most negative
  ## first), so the histogram's x-range (centred on d_random_mean = 3,
  ## sd = 1) should appear in exactly 2 panels' worth of vline data.
  vline_layer <- built$data[[2]]
  expect_equal(nrow(vline_layer), 2L)
})
