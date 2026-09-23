test_that("compounds_classify() does not reuse a different structure's cached family", {
  fail_lookup <- FALSE
  testthat::local_mocked_bindings(
    .npclassifier_lookup = function(smiles) {
      if (fail_lookup) stop("NPClassifier unavailable")
      list(pathway = "original family", superclass = "original superclass",
           class = "original class", isglycoside = FALSE)
    },
    .package = "patliR"
  )
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_single_compound(), identifier = "pubchem")
  proj <- compounds_classify(proj)
  expect_equal(patliRResults(proj, "compounds_classified")$pathway, "original family")

  fail_lookup <- TRUE
  expect_warning(cached <- compounds_classify(proj))
  expect_equal(patliRResults(cached, "compounds_classified")$pathway, "original family")

  ## Keep the compound ID but replace its structure. An offline retry
  ## must not attach the old structure's cached classification to it.
  cmp <- compounds(proj)
  cmp$canonical_smiles <- "CCO"
  cmp$smiles <- "CCO"
  compounds(proj) <- cmp
  fail_lookup <- TRUE
  expect_warning(proj <- compounds_classify(proj))
  expect_true(is.na(patliRResults(proj, "compounds_classified")$pathway))
})
