test_that("network_hub_penalty() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_hub_penalty(proj), "network_edges")
})

test_that("network_hub_penalty() does not require network_centrality() to have run", {
  proj <- .network_stats_test_setup()
  expect_false("network_centrality" %in% names(patliRResults(proj)))
  expect_no_error(network_hub_penalty(proj, condition = "FLO-ET"))
})

test_that("network_hub_penalty() matches the documented formula and bounds", {
  proj <- .network_stats_test_setup()
  proj <- network_hub_penalty(proj, condition = "FLO-ET")
  result <- patliRResults(proj, "network_hub_penalty")

  expect_true(all(c("condition", "uniprot_id", "degree_raw", "n_compounds_total", "score_adjusted") %in% names(result)))
  expect_true(all(result$degree_raw <= result$n_compounds_total))

  expected <- result$degree_raw * log(result$n_compounds_total / result$degree_raw)
  expect_equal(result$score_adjusted, expected)

  ## a target touching every present compound is maximally promiscuous -> score 0
  fully_promiscuous <- result[result$degree_raw == result$n_compounds_total, , drop = FALSE]
  if (nrow(fully_promiscuous) > 0) {
    expect_true(all(abs(fully_promiscuous$score_adjusted) < 1e-9))
  }
})

test_that("network_hub_penalty() rebuilding one condition does not touch the others", {
  proj <- .network_stats_test_setup()
  proj <- network_hub_penalty(proj)
  before <- patliRResults(proj, "network_hub_penalty")

  proj <- network_hub_penalty(proj, condition = "FLO-ET")
  after <- patliRResults(proj, "network_hub_penalty")

  expect_setequal(unique(after$condition), unique(before$condition))
})
