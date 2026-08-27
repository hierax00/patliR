test_that("prep_binarize() parses R<n>-<CONDITION> columns and averages replicates", {
  proj <- .test_project()
  proj <- prep_binarize(proj, .test_abundance_matrix())

  mr <- matrixRaw(proj)
  bin <- binarizedMatrix(proj)

  expected_conditions <- c("LEA-ET", "LEA-AQ", "FLO-ET", "FLO-AQ", "STE-ET")
  expect_true(all(expected_conditions %in% names(mr)))
  expect_true(all(expected_conditions %in% names(bin)))
  expect_equal(nrow(mr), 8) # 8 compounds in the example matrix

  ## Apigenin (row 1) is 0 across all LEA-ET replicates.
  expect_equal(mr$`LEA-ET`[mr$compound_id == "Apigenin"], 0)
  expect_equal(bin$`LEA-ET`[bin$compound_id == "Apigenin"], 0)

  ## Every binarized value must be 0 or 1.
  bin_values <- unlist(bin[expected_conditions])
  expect_true(all(bin_values %in% c(0, 1)))
})

test_that("prep_binarize() errors on columns that do not match the naming convention", {
  proj <- .test_project()
  bad <- data.frame(Name = "X", `LEA_ET_R1` = 10, check.names = FALSE)
  expect_error(prep_binarize(proj, bad), "R<n>-<CONDITION>", fixed = TRUE)
})

test_that("prep_binarize() logs unmatched compound names when compounds(proj) is populated", {
  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem") # only Apigenin
  proj <- prep_binarize(proj, .test_abundance_matrix())

  log <- projectLog(proj)
  expect_true(any(grepl("compound_not_found", log$message)))
  ## Apigenin itself should have matched and use the internal id.
  expect_true("C0001" %in% matrixRaw(proj)$compound_id)
})

test_that("an already-binary column is passed through without Q1 thresholding", {
  proj <- .test_project()
  data <- data.frame(
    Name = c("A", "B", "C"),
    `R1-COND` = c(0, 1, 1),
    check.names = FALSE
  )
  proj <- prep_binarize(proj, data)
  expect_equal(binarizedMatrix(proj)$COND, c(0, 1, 1))
  log <- projectLog(proj)
  expect_true(any(grepl("already_binary=TRUE", log$message)))
})
