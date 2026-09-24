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

test_that(".npclassifier_parse_body() accepts valid JSON regardless of the Content-Type header", {
  skip_if_not_installed("httr2")
  json <- '{"class_results": ["Hydrocarbons"], "superclass_results": ["Fatty acyls"], "pathway_results": ["Fatty acids"], "isglycoside": false}'
  as_type <- function(type, body) httr2::response(status_code = 200L, headers = list(`Content-Type` = type), body = charToRaw(body))

  for (type in c("application/json", "text/html; charset=utf-8")) {
    parsed <- patliR:::.npclassifier_parse_body(as_type(type, json))
    expect_equal(parsed$pathway_results[[1]], "Fatty acids")
    expect_false(isTRUE(parsed$isglycoside))
  }
  expect_null(patliR:::.npclassifier_parse_body(as_type("text/html", "<html>rate limited</html>")))
  expect_null(patliR:::.npclassifier_parse_body(as_type("application/json", '{"error": "x"}')))
})
