test_that("refdb_build() rejects 'coconut' as a source with an informative error", {
  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem")
  expect_error(
    refdb_build(proj, sources = "coconut"),
    "not implemented yet"
  )
})

test_that("refdb_build() warns and returns proj unchanged when there are no compounds", {
  proj <- .test_project()
  expect_warning(out <- refdb_build(proj, sources = "pubchem"), "No matching compounds")
  expect_equal(nrow(compounds(out)), 0)
})

test_that("refdb_update() requires at least one compound_id", {
  proj <- .test_project()
  expect_error(refdb_update(proj, character(0)))
})

test_that("refdb_rebuild_cache() warns when there is nothing to rebuild from", {
  proj <- .test_project()
  expect_warning(refdb_rebuild_cache(proj), "No reference database CSVs found")
})

test_that(".chembl_url_encode_smiles() leaves '=', '(', ')' untouched", {
  ## Fully percent-encoding these (utils::URLencode(smiles, reserved = TRUE))
  ## broke the live ChEMBL similarity-search endpoint with HTTP 500,
  ## confirmed against the real API -- see R/refdb.R.
  expect_equal(patliR:::.chembl_url_encode_smiles("C1=CC=CC=C1"), "C1=CC=CC=C1")
  expect_equal(patliR:::.chembl_url_encode_smiles("CC(=O)Oc1ccccc1C(=O)O"), "CC(=O)Oc1ccccc1C(=O)O") # aspirin
})

test_that(".chembl_url_encode_smiles() escapes '#', '/', and '%'", {
  expect_equal(patliR:::.chembl_url_encode_smiles("C#N"), "C%23N")
  expect_equal(patliR:::.chembl_url_encode_smiles("C1CC1/C=C/C"), "C1CC1%2FC=C%2FC")
  expect_equal(patliR:::.chembl_url_encode_smiles("C%10"), "C%2510") # ring-closure SMILES like %10
})

test_that(".chembl_url_encode_smiles() escapes a literal '%' first, before '#'/'/' introduce their own '%' escapes", {
  ## If '%' were not escaped first, the '%23'/'%2F' this function itself
  ## produces would get double-escaped into '%2523'/'%252F'.
  expect_equal(patliR:::.chembl_url_encode_smiles("%#"), "%25%23")
})

## Live-network tests (actual PubChem/ChEMBL calls) are intentionally left
## out of the automated suite -- see TESTING_GUIDE.Rmd for how we will
## exercise refdb_build()/refdb_update() together against the real APIs.
## Note (2026-08-11): for real natural-product compounds (real_data/), a
## "no ChEMBL molecule found for this SMILES (flexmatch)" warning per
## compound is EXPECTED, not a bug -- .chembl_lookup() only ever does a
## SMILES flexmatch (no PubChemCID branch, see R/refdb.R), and most
## natural-product secondary metabolites genuinely have no ChEMBL entry.
## warn_and_cache (the default fetch_mode) is designed exactly for this:
## log it, keep going, do not fail the pipeline.
