test_that("network_synergy() requires network_proximity() to have run first", {
  proj <- .network_stats_test_setup()
  expect_error(network_synergy(proj, condition = "FLO-ET", disease = "EFO_0000000"), "network_proximity")
})

test_that("network_synergy() errors clearly when network_proximity() ran for a different condition/disease", {
  proj <- .network_stats_test_setup()
  ## Fabricate network_proximity results directly -- bypassing the real
  ## STRINGdb-dependent computation -- so this test is self-contained and
  ## offline, same pattern as test-network_proximity.R's fake targets_disease.
  ct <- unique(patliRResults(proj, "network_edges")[, c("compound_id", "uniprot_id")])
  compounds <- unique(ct$compound_id)
  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = compounds, disease_id = "SOME_DISEASE",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, d_observed = 2,
    d_random_mean = 3, d_random_sd = 1, z_score = seq(-2, 2, length.out = length(compounds)),
    n_random = 100L, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  expect_error(
    network_synergy(proj, condition = "FLO-ET", disease = "NOT-THE-SAME-DISEASE"),
    "network_proximity"
  )
})

test_that("network_synergy() with pairs = 'rank_top' and pairs = 'all' both run end-to-end", {
  proj <- .network_stats_test_setup()
  ## Filter to FLO-ET specifically -- network_edges spans every condition
  ## network_build() built, and network_synergy() only pairs compounds
  ## within the requested condition, so the expected pair count below must
  ## be based on FLO-ET's own compound set, not every compound in proj.
  edges_flo_et <- patliRResults(proj, "network_edges")
  edges_flo_et <- edges_flo_et[edges_flo_et$condition == "FLO-ET", , drop = FALSE]
  ct <- unique(edges_flo_et[, c("compound_id", "uniprot_id")])
  compounds <- unique(ct$compound_id)
  skip_if(length(compounds) < 2, "needs at least 2 compounds in the FLO-ET condition")

  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = compounds, disease_id = "SOME_DISEASE",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, d_observed = 2,
    d_random_mean = 3, d_random_sd = 1, z_score = seq(-2, 2, length.out = length(compounds)),
    n_random = 100L, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  proj_top <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", pairs = "rank_top", top_n = 2)
  result_top <- patliRResults(proj_top, "network_synergy")
  expect_true(all(c("condition", "disease_id", "compound_a", "compound_b", "target_jaccard",
                     "complementarity", "z_score_a", "z_score_b", "joint_closeness",
                     "synergy_score", "pairs_mode") %in% names(result_top)))
  expect_true(all(result_top$pairs_mode == "rank_top"))
  expect_true(all(result_top$compound_a != result_top$compound_b))

  proj_all <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", pairs = "all")
  result_all <- patliRResults(proj_all, "network_synergy")
  expect_equal(nrow(result_all), choose(length(compounds), 2))
  expect_true(all(result_all$complementarity >= 0 & result_all$complementarity <= 1))

  ## Cheng conjunction: synergy_score is scored iff BOTH z_score < 0, and
  ## is never a negative-signed "far from disease" artefact.
  expect_true("both_proximal" %in% names(result_all))
  scored <- !is.na(result_all$synergy_score)
  expect_equal(scored, result_all$both_proximal)
  expect_equal(result_all$both_proximal,
               result_all$z_score_a < 0 & result_all$z_score_b < 0)
  expect_true(all(result_all$synergy_score[scored] >= 0))
})

test_that("network_synergy(pairs = 'all') does not error when a compound is absent from network_proximity", {
  proj <- .network_stats_test_setup()
  edges_flo_et <- patliRResults(proj, "network_edges")
  edges_flo_et <- edges_flo_et[edges_flo_et$condition == "FLO-ET", , drop = FALSE]
  compounds <- unique(edges_flo_et$compound_id)
  skip_if(length(compounds) < 3, "needs at least 3 compounds in FLO-ET")

  ## Fabricate proximity results for every compound EXCEPT the last one --
  ## network_proximity() legitimately skips compounds with no
  ## STRING-mappable target, so `pairs = "all"` must yield NA for such a
  ## compound's pairs, not `subscript out of bounds` (spec headline #3).
  scored_compounds <- compounds[-length(compounds)]
  fake_prox <- data.frame(
    condition = "FLO-ET", compound_id = scored_compounds, disease_id = "SOME_DISEASE",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, d_observed = 2,
    d_random_mean = 3, d_random_sd = 1,
    z_score = seq(-2, 1, length.out = length(scored_compounds)),
    n_random = 100L, seed_used = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- fake_prox

  expect_no_error(
    proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", pairs = "all")
  )
  res <- patliRResults(proj, "network_synergy")
  missing_compound <- compounds[length(compounds)]
  missing_rows <- res[res$compound_a == missing_compound | res$compound_b == missing_compound, , drop = FALSE]
  expect_true(nrow(missing_rows) > 0)
  ## the absent compound's own z is NA on every pair it appears in
  za_is_missing <- ifelse(missing_rows$compound_a == missing_compound, missing_rows$z_score_a, missing_rows$z_score_b)
  expect_true(all(is.na(za_is_missing)))
  expect_true(all(!missing_rows$both_proximal))
  expect_true(all(is.na(missing_rows$synergy_score)))
})
