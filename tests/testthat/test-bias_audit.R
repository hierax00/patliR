## Fabricate a minimal reference_bioactivity table directly (bypassing a
## live PubChem/ChEMBL call, same offline pattern as
## test-network_proximity.R's fake_disease) -- a hub target T_HUB hit by
## every one of 10 compounds (degree 10), and 5 more targets with a
## non-tied, spread-out degree distribution (1..5) so stats::mad() does not
## collapse to zero (ties at the mode zero out MAD -- a real edge case
## worth keeping in mind, exercised separately below).
.bias_test_bioactivity <- function() {
  hub_rows <- data.frame(
    compound_id = sprintf("C%02d", 1:10), target_chembl_id = "T_HUB", stringsAsFactors = FALSE
  )
  spread_rows <- do.call(rbind, lapply(1:5, function(n_hit) {
    data.frame(
      compound_id = sprintf("C%02d", seq_len(n_hit)),
      target_chembl_id = paste0("T_", n_hit),
      stringsAsFactors = FALSE
    )
  }))
  rbind(hub_rows, spread_rows)
}

.bias_test_project <- function() {
  proj <- .test_project()
  patliRResults(proj, "reference_bioactivity") <- .bias_test_bioactivity()
  proj
}

test_that("bias_audit() requires reference_bioactivity to exist", {
  proj <- .test_project()
  expect_error(bias_audit(proj), "reference_bioactivity")
})

test_that("bias_audit() rejects a non-NULL categories with an informative error", {
  proj <- .bias_test_project()
  expect_error(bias_audit(proj, categories = "mesh_disease"), "not implemented yet")
})

test_that("bias_audit() rejects check_homogeneity = FALSE with categories = NULL (nothing to do)", {
  proj <- .bias_test_project()
  expect_error(bias_audit(proj, check_homogeneity = FALSE), "Nothing to audit")
})

test_that("bias_audit() flags the hub target as 'promiscuo' and leaves the spread targets 'normal'", {
  proj <- .bias_test_project()
  proj <- bias_audit(proj, mad_threshold = 2.5)
  homog <- patliRResults(proj, "bias_homogeneity")

  expect_true(all(c("id", "tipo", "frecuencia_global_refdb", "mad_score", "categoria") %in% names(homog)))
  expect_setequal(unique(homog$tipo), c("compound", "target"))

  hub_row <- homog[homog$tipo == "target" & homog$id == "T_HUB", ]
  expect_equal(nrow(hub_row), 1)
  expect_equal(hub_row$frecuencia_global_refdb, 10L)
  expect_identical(hub_row$categoria, "promiscuo")

  spread_targets <- homog[homog$tipo == "target" & homog$id != "T_HUB", ]
  expect_true(all(spread_targets$categoria == "normal"))
})

test_that("bias_audit() handles a fully-tied degree distribution (MAD = 0) without crashing", {
  proj <- .test_project()
  ## Every compound hits exactly one, distinct target -- degree 1 for
  ## everyone on both sides, so stats::mad() is 0 and mad_score must come
  ## out NA (never Inf/NaN from a division by zero).
  patliRResults(proj, "reference_bioactivity") <- data.frame(
    compound_id = sprintf("C%02d", 1:5), target_chembl_id = sprintf("T%02d", 1:5),
    stringsAsFactors = FALSE
  )
  proj <- bias_audit(proj)
  homog <- patliRResults(proj, "bias_homogeneity")
  expect_true(all(is.na(homog$mad_score)))
  expect_true(all(homog$categoria == "normal"))
})

test_that("bias_reweight() requires bias_audit() to have run first", {
  proj <- .bias_test_project()
  expect_error(bias_reweight(proj), "bias_homogeneity")
})

test_that("bias_reweight() matches the documented log-ratio formula for 'promiscuo' entries, and leaves 'normal' entries unadjusted", {
  proj <- .bias_test_project()
  proj <- bias_audit(proj, mad_threshold = 2.5)
  proj <- bias_reweight(proj)
  result <- patliRResults(proj, "bias_reweighted")

  expect_true(all(c("id", "tipo", "score_crudo", "score_ajustado", "categoria") %in% names(result)))

  normal_rows <- result[result$categoria == "normal", ]
  expect_equal(normal_rows$score_ajustado, normal_rows$score_crudo)

  hub_row <- result[result$tipo == "target" & result$id == "T_HUB", ]
  n_compounds_total <- sum(result$tipo == "compound")
  expected <- hub_row$score_crudo * log(n_compounds_total / hub_row$score_crudo)
  expect_equal(hub_row$score_ajustado, expected)
  expect_true(hub_row$score_ajustado < hub_row$score_crudo)
})

test_that("bias_report() warns and returns an empty (but present) summary when nothing has run yet", {
  proj <- .test_project()
  expect_warning(report <- bias_report(proj), "bias_reweighted")
  expect_equal(nrow(report$summary), 0)
  expect_true(all(c("id", "tipo", "score_crudo", "score_ajustado", "categoria") %in% names(report$summary)))
  expect_true(nzchar(report$note))
})

test_that("bias_report() returns only 'promiscuo' entries, ordered by score_ajustado ascending, and never mutates proj", {
  proj <- .bias_test_project()
  proj <- bias_audit(proj, mad_threshold = 2.5)
  proj <- bias_reweight(proj)
  before <- patliRResults(proj)

  report <- bias_report(proj)

  expect_true(all(report$summary$categoria == "promiscuo"))
  expect_true("T_HUB" %in% report$summary$id)
  expect_true(!is.unsorted(report$summary$score_ajustado))
  expect_identical(patliRResults(proj), before)
})
