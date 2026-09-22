test_that("patliR_project() creates an empty project with the given directories", {
  proj <- .test_project()
  expect_s4_class(proj, "PatliRProject")
  expect_true(dir.exists(projectDir(proj)))
  expect_true(dir.exists(cacheDir(proj)))
})

test_that("patliR_load() round-trips compounds()'s text columns without read.csv() re-typing them (regression: pubchem_id/name/etc used to lose type on reload)", {
  proj <- .test_project()
  ## A PubChemCID with a leading zero (an all-digit string, so
  ## read.csv()'s value inference reads it back as integer -- dropping the
  ## zero -- unless the column is pinned to character) and a compound
  ## literally named "TRUE" (read back as logical without pinning). Both
  ## survive prep_compounds() unchanged (it does as.character() on its
  ## inputs) and are written to 01_compounds.csv as-is; the risk is purely
  ## at the patliR_load() read boundary.
  tricky <- data.frame(
    Name = c("TRUE", "Normal Name"),
    PubChemCID = c("000123", "5280443"),
    SMILES = c("CCO", "c1ccccc1"),
    stringsAsFactors = FALSE
  )
  proj <- prep_compounds(proj, tricky, identifier = "pubchem")
  expect_true(file.exists(file.path(projectDir(proj), "01_compounds.csv")))

  reloaded <- patliR_load(projectDir(proj))
  cmp <- compounds(reloaded)

  expect_type(cmp$pubchem_id, "character")
  expect_type(cmp$name, "character")
  expect_true("000123" %in% cmp$pubchem_id) ## leading zero preserved
  expect_true("TRUE" %in% cmp$name) ## not coerced to logical
})

test_that("patliR_load() round-trips matrixRaw()/binarizedMatrix()'s compound_id column as character", {
  proj <- .test_project()
  tricky <- data.frame(
    Name = "Some Compound", PubChemCID = "000123", SMILES = "CCO",
    stringsAsFactors = FALSE
  )
  proj <- prep_compounds(proj, tricky, identifier = "pubchem")
  ## compound_id itself is always "C%04d"-style (never purely numeric, see
  ## .next_compound_ids()), so this asserts the colClasses pin is present
  ## and harmless, not that it fixes an observable symptom on this column
  ## specifically.
  abundance <- data.frame(
    Name = "Some Compound", `R1-COND` = 1,
    check.names = FALSE
  )
  proj <- prep_binarize(proj, abundance, id_col = "Name")
  expect_true(file.exists(file.path(projectDir(proj), "03_binarized.csv")))

  reloaded <- patliR_load(projectDir(proj))
  expect_type(binarizedMatrix(reloaded)$compound_id, "character")
  expect_type(matrixRaw(reloaded)$compound_id, "character")
})
