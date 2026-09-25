## rank_candidates() combines already-tested source tables, so these tests
## focus on the composition logic itself (criteria detection, roll-up,
## direction normalisation, RRA/Pareto correctness) rather than re-testing
## e.g. network_centrality()'s own formulas. Optional-criteria tables are
## built directly via patliRResults(proj, "<name>") <- <fake df>, same
## convention test-network_synergy.R uses for its network_proximity fixture
## -- avoids needing STRINGdb/disease_genes_fetch()/GOSemSim for a test that
## only needs a couple of specific column values.

testthat::skip_if_not_installed("RobustRankAggreg")

.rank_test_setup <- function() {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  proj <- adme_local(proj)
  adme_filter(proj)
}

.rank_flo_compounds <- function(proj) {
  e <- patliRResults(proj, "network_edges")
  sort(unique(e$compound_id[e$condition == "FLO-ET"]))
}

test_that("rank_candidates() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(rank_candidates(proj), "network_edges")
})

test_that("rank_candidates() requires adme_filtered", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj)
  expect_error(rank_candidates(proj, condition = "FLO-ET"), "adme_filtered")
})

test_that("rank_candidates() requires network_centrality", {
  proj <- .network_stats_test_setup()
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  expect_error(rank_candidates(proj, condition = "FLO-ET"), "network_centrality")
})

test_that("rank_candidates() aborts when network_centrality was run with normalize = FALSE", {
  proj <- .network_stats_test_setup()
  proj <- network_centrality(proj, normalize = FALSE)
  proj <- adme_local(proj)
  proj <- adme_filter(proj)
  expect_error(rank_candidates(proj, condition = "FLO-ET"), "normalis")
})

test_that("rank_candidates() aborts if criteria= drops a mandatory criterion", {
  proj <- .rank_test_setup()
  expect_error(rank_candidates(proj, condition = "FLO-ET", criteria = "centrality"), "mandatory")
})

test_that("rank_candidates() with only mandatory criteria produces a well-formed table", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  proj <- rank_candidates(proj, condition = "FLO-ET")
  rc <- patliRResults(proj, "rank_candidates")

  expect_setequal(rc$compound_id, cps)
  expect_true(all(c(
    "crit_adme_pass_frac", "crit_centrality", "rank_adme_pass_frac", "rank_centrality",
    "n_p2_partners", "n_criteria_used", "rra_score", "rra_rank", "pareto_front"
  ) %in% names(rc)))
  expect_false(any(c("crit_hub_penalty", "crit_proximity_z", "crit_synergy_best", "crit_module_r_index") %in% names(rc)))
  expect_true(all(rc$n_criteria_used == 2L))
  expect_equal(rc$rank_adme_pass_frac, rank(-rc$crit_adme_pass_frac, ties.method = "average"))
  expect_equal(rc$rank_centrality, rank(-rc$crit_centrality, ties.method = "average"))
  expect_true(all(rc$pareto_front >= 1L))
  expect_true(all(is.na(rc$n_p2_partners)))
})

test_that("rank_candidates() logs which optional criteria were skipped when their tables are absent", {
  proj <- .rank_test_setup()
  proj <- rank_candidates(proj, condition = "FLO-ET")
  log <- projectLog(proj)
  msg <- log$message[log$step == "rank_candidates"][1]
  expect_match(msg, "skipped \\(no data\\)")
  expect_match(msg, "hub_penalty")
  expect_match(msg, "proximity")
  expect_match(msg, "synergy")
  expect_match(msg, "module_robustness")
})

test_that("rank_candidates() auto-detects an optional criterion once its table exists", {
  proj <- .rank_test_setup()
  proj <- network_hub_penalty(proj, condition = "FLO-ET")
  proj <- rank_candidates(proj, condition = "FLO-ET")
  rc <- patliRResults(proj, "rank_candidates")
  expect_true("crit_hub_penalty" %in% names(rc))
  expect_true(all(rc$n_criteria_used == 3L))
  log <- projectLog(proj)
  msg <- log$message[log$step == "rank_candidates"][1]
  expect_match(msg, "hub_penalty")
  skipped_part <- sub(".*skipped \\(no data\\): ", "", msg)
  expect_false(grepl("hub_penalty", skipped_part))
})

test_that("rank_candidates() roll_up = 'max' never scores a compound below roll_up = 'weighted_mean'", {
  proj <- .rank_test_setup()
  proj_wm <- rank_candidates(proj, condition = "FLO-ET", roll_up = "weighted_mean")
  proj_mx <- rank_candidates(proj, condition = "FLO-ET", roll_up = "max")
  rc_wm <- patliRResults(proj_wm, "rank_candidates")
  rc_mx <- patliRResults(proj_mx, "rank_candidates")
  ord <- match(rc_wm$compound_id, rc_mx$compound_id)
  expect_true(all(rc_mx$crit_centrality[ord] >= rc_wm$crit_centrality - 1e-9))
})

test_that("rank_candidates() output shows the raw z_score for proximity, direction-flipped only for ranking", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  ## An ascending, all-distinct sequence -- length tracks cps so this
  ## doesn't assume a fixed fixture size; preserves c(-2, -1, 0, 1) when
  ## cps has exactly 4 elements.
  prox <- data.frame(condition = "FLO-ET", compound_id = cps, disease_id = "D1",
                      z_score = seq(-2, by = 1, length.out = length(cps)), stringsAsFactors = FALSE)
  patliRResults(proj, "network_proximity") <- prox
  proj <- rank_candidates(proj, condition = "FLO-ET", disease = "D1")
  rc <- patliRResults(proj, "rank_candidates")
  expect_setequal(rc$crit_proximity_z, prox$z_score)
  ## most negative z_score (best proximity) must rank 1st on that criterion
  best <- rc$compound_id[rc$rank_proximity_z == 1]
  expect_identical(best, prox$compound_id[which.min(prox$z_score)])
})

test_that("rank_candidates() aborts when network_proximity has multiple diseases and none is specified", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  prox <- rbind(
    data.frame(condition = "FLO-ET", compound_id = cps, disease_id = "D1", z_score = -1, stringsAsFactors = FALSE),
    data.frame(condition = "FLO-ET", compound_id = cps, disease_id = "D2", z_score = 1, stringsAsFactors = FALSE)
  )
  patliRResults(proj, "network_proximity") <- prox
  expect_error(rank_candidates(proj, condition = "FLO-ET"), "multiple diseases")
  expect_no_error(rank_candidates(proj, condition = "FLO-ET", disease = "D1"))
})

test_that("rank_candidates() synergy criterion takes the best (max) partner score and counts P2 partnerships", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  testthat::skip_if(length(cps) < 3, "fixture needs at least 3 compounds")
  a <- cps[1]; b <- cps[2]; d <- cps[3]

  prox <- data.frame(condition = "FLO-ET", compound_id = cps, disease_id = "D1", z_score = -1, stringsAsFactors = FALSE)
  syn <- data.frame(
    condition = "FLO-ET", disease_id = "D1",
    compound_a = c(a, a), compound_b = c(b, d),
    synergy_score = c(0.2, 0.9), cheng_class = c("P2", "P2"),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- prox
  patliRResults(proj, "network_synergy") <- syn
  proj <- rank_candidates(proj, condition = "FLO-ET", disease = "D1", criteria = c("adme", "centrality", "synergy"))
  rc <- patliRResults(proj, "rank_candidates")

  expect_equal(rc$crit_synergy_best[rc$compound_id == a], 0.9)
  expect_equal(rc$n_p2_partners[rc$compound_id == a], 2L)
  expect_equal(rc$crit_synergy_best[rc$compound_id == b], 0.2)
  expect_true(all(is.na(rc$crit_synergy_best[!rc$compound_id %in% c(a, b, d)])))
})

test_that("rank_candidates() module_robustness joins r_index via network_module_membership", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  testthat::skip_if(length(cps) < 4, "fixture needs at least 4 compounds for this module split")
  ## Deliberately a controlled 4-compound/2-module scenario -- restricted
  ## to the first 4 of cps regardless of how many the fixture actually
  ## has, so this test doesn't depend on the ambient fixture's size.
  membership <- data.frame(
    condition = "FLO-ET", node_id = cps[1:4], node_type = "compound",
    module_id = c("m1", "m1", "m2", "m2"), module_type = "cluster",
    stringsAsFactors = FALSE
  )
  robustness <- data.frame(condition = "FLO-ET", module_id = c("m1", "m2"), r_index = c(0.8, 0.2), stringsAsFactors = FALSE)
  patliRResults(proj, "network_module_membership") <- membership
  patliRResults(proj, "network_module_robustness") <- robustness
  proj <- rank_candidates(proj, condition = "FLO-ET", criteria = c("adme", "centrality", "module_robustness"))
  rc <- patliRResults(proj, "rank_candidates")

  expect_equal(rc$crit_module_r_index[rc$compound_id == cps[1]], 0.8)
  expect_equal(rc$crit_module_r_index[rc$compound_id == cps[3]], 0.2)
})

test_that("rank_candidates() rewards a compound that dominates on one criterion with Pareto front 1", {
  proj <- .rank_test_setup()
  cps <- .rank_flo_compounds(proj)
  best <- cps[1]
  prox <- data.frame(condition = "FLO-ET", compound_id = cps, disease_id = "D1",
                      z_score = ifelse(cps == best, -5, 1), stringsAsFactors = FALSE)
  patliRResults(proj, "network_proximity") <- prox
  proj <- rank_candidates(proj, condition = "FLO-ET", disease = "D1")
  rc <- patliRResults(proj, "rank_candidates")
  expect_equal(rc$pareto_front[rc$compound_id == best], 1L)
})

test_that("rank_candidates() rebuilding one condition does not touch the others", {
  proj <- .rank_test_setup()
  conditions <- unique(patliRResults(proj, "network_edges")$condition)
  testthat::skip_if(length(conditions) < 2, "fixture needs at least 2 conditions")
  proj <- rank_candidates(proj)
  before <- patliRResults(proj, "rank_candidates")

  proj <- rank_candidates(proj, condition = conditions[1])
  after <- patliRResults(proj, "rank_candidates")

  expect_setequal(unique(after$condition), unique(before$condition))
})

test_that("rank_candidates() export = 'sdf'/'smi' writes the top-N structures", {
  proj <- .rank_test_setup()
  proj <- rank_candidates(proj, condition = "FLO-ET", top_n = 2, export = "sdf")
  sdf_path <- file.path(projectDir(proj), "results", "rank_candidates_top2_FLO-ET.sdf")
  expect_true(file.exists(sdf_path))
  expect_equal(sum(readLines(sdf_path) == "$$$$"), 2)

  proj <- rank_candidates(proj, condition = "FLO-ET", top_n = 2, export = "smi")
  smi_path <- file.path(projectDir(proj), "results", "rank_candidates_top2_FLO-ET.smi")
  expect_true(file.exists(smi_path))
  expect_equal(length(readLines(smi_path)), 2)
})

test_that(".rank_pareto_front() identifies non-dominated tiers correctly", {
  mat <- matrix(c(
    3, 1,
    3, 2,
    1, 3,
    1, 1
  ), ncol = 2, byrow = TRUE)
  front <- patliR:::.rank_pareto_front(mat)
  expect_equal(front, c(2L, 1L, 1L, 3L))
})

test_that("rank_candidates() ties are invariant to compound renaming and input order", {
  testthat::skip_if_not_installed("RobustRankAggreg")
  run_fixture <- function(ids, row_order = seq_along(ids), proximity = c(-2, -2, NA, NA)) {
    proj <- .test_project()
    edges <- data.frame(condition = "test", compound_id = ids,
                        uniprot_id = paste0("T", seq_along(ids)), weight = 1)
    adme <- data.frame(compound_id = ids, pass = c(TRUE, TRUE, TRUE, FALSE))
    cent <- data.frame(condition = "test", node_type = "target",
                       node_id = edges$uniprot_id, degree_norm = c(1, 1, 1, 0))
    prox <- data.frame(condition = "test", compound_id = ids,
                       disease_id = "D1", z_score = proximity)
    patliRResults(proj, "network_edges") <- edges[row_order, ]
    patliRResults(proj, "adme_filtered") <- adme[row_order, ]
    patliRResults(proj, "network_centrality") <- cent[row_order, ]
    patliRResults(proj, "network_proximity") <- prox[row_order, ]
    proj <- rank_candidates(proj, condition = "test",
                            criteria = c("adme", "centrality", "proximity"))
    out <- patliRResults(proj, "rank_candidates")
    out <- out[match(ids, out$compound_id), ]
    rownames(out) <- NULL
    out
  }
  original <- run_fixture(c("a", "b", "c", "d"))
  reordered <- run_fixture(c("a", "b", "c", "d"), c(4, 2, 1, 3))
  renamed <- run_fixture(c("z", "x", "y", "w"), c(3, 1, 4, 2))
  columns <- setdiff(names(original), "compound_id")
  expect_equal(reordered[columns], original[columns])
  expect_equal(renamed[columns], original[columns])
  expect_equal(renamed$rra_rank[1:3], original$rra_rank[1:3])
  expect_equal(original$rank_adme_pass_frac, c(2, 2, 2, 4))
  expect_equal(original$rank_centrality, c(2, 2, 2, 4))
  expect_equal(original$rank_proximity_z, c(1.5, 1.5, 4, 4))
  expect_equal(original$n_criteria_used, c(3, 3, 2, 2))
  expect_equal(original$rra_score[1], original$rra_score[2])
  expect_equal(original$rra_rank[1], original$rra_rank[2])

  all_missing <- run_fixture(c("a", "b", "c", "d"), proximity = rep(NA_real_, 4))
  expect_equal(all_missing$rank_proximity_z, rep(4, 4))
  expect_equal(all_missing$n_criteria_used, rep(2, 4))
  expect_true(all(is.finite(all_missing$rra_score)))
})

test_that(".rank_rra_rank() breaks saturated RRA ties by mean per-criterion rank, keeping min ties", {
  score <- c(0.1, 1, 1, 1, 1)
  mean_rank <- c(0.5, 0.9, 0.3, 0.6, 0.6)
  expect_equal(patliR:::.rank_rra_rank(score, mean_rank), c(1L, 5L, 2L, 3L, 3L))
  ## without ties on the score it is the plain rank
  expect_equal(patliR:::.rank_rra_rank(c(0.3, 0.1, 0.2), c(1, 1, 1)), c(3L, 1L, 2L))
})
