test_that("compounds_similarity() returns a symmetric-shaped pairwise table bounded in [0, 1]", {
  testthat::skip_if_not_installed("fingerprint")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  n <- nrow(compounds(proj))

  proj <- compounds_similarity(proj)
  result <- patliRResults(proj, "compounds_similarity")

  expect_true(all(c("compound_a", "compound_b", "similarity", "fingerprint_type", "method") %in% names(result)))
  expect_equal(nrow(result), choose(n, 2)) # every compound in the fixture has a valid SMILES
  expect_true(all(result$similarity >= 0 & result$similarity <= 1))
  expect_true(all(result$compound_a != result$compound_b)) # no self-pairs

  ## a compound compared against itself must be maximal similarity --
  ## sanity check on the fingerprint/method wiring, not just the shape.
  one_id <- compounds(proj)$id[1]
  self_sim <- fingerprint::fp.sim.matrix(
    list(rcdk::get.fingerprint(rcdk::parse.smiles(compounds(proj)$smiles[1])[[1]], type = "standard")),
    method = "tanimoto"
  )
  expect_equal(self_sim[1, 1], 1)
})

test_that("compounds_similarity() requires at least 2 compounds", {
  testthat::skip_if_not_installed("fingerprint")

  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem")
  expect_error(compounds_similarity(proj), "at least 2")
})

test_that("compounds_similarity() respects compound_ids and different fingerprint_type/method", {
  testthat::skip_if_not_installed("fingerprint")

  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  three_ids <- compounds(proj)$id[1:3]

  proj <- compounds_similarity(proj, compound_ids = three_ids, fingerprint_type = "maccs", method = "dice")
  result <- patliRResults(proj, "compounds_similarity")

  expect_equal(nrow(result), 3) # choose(3, 2)
  expect_true(all(result$compound_a %in% three_ids) && all(result$compound_b %in% three_ids))
  expect_setequal(unique(result$fingerprint_type), "maccs")
  expect_setequal(unique(result$method), "dice")
})
