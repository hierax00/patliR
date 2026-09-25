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

test_that("prep_as_condition() turns a plain compound list into one network-ready condition", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  n <- nrow(compounds(proj))

  proj <- prep_as_condition(proj, condition = "my_extract")
  bin <- binarizedMatrix(proj)
  expect_true("my_extract" %in% names(bin))
  expect_equal(sum(bin$my_extract), n)
  expect_equal(nrow(bin), n)

  ## a second call adds a column, does not wipe the first
  sub <- compounds(proj)$id[1:2]
  proj <- prep_as_condition(proj, condition = "subset", compound_ids = sub)
  bin <- binarizedMatrix(proj)
  expect_true(all(c("my_extract", "subset") %in% names(bin)))
  expect_equal(sum(bin$subset), 2)
})

test_that("prep_as_condition() errors before prep_compounds()", {
  proj <- .test_project()
  expect_error(prep_as_condition(proj), "run .*prep_compounds")
})

test_that("prep_binarize() includes zero abundances in the condition quantile", {
  data <- data.frame(Name = letters[1:5], `R1-X` = c(0, 0, 2, 10, NA), check.names = FALSE)
  proj <- prep_binarize(.test_project(), data)
  expect_equal(binarizedMatrix(proj)$X, c(0L, 0L, 1L, 1L, NA_integer_))
  data$`R1-X` <- c(0, 2, 4, 6, NA)
  proj <- prep_binarize(.test_project(), data, q = 0.5)
  expect_equal(binarizedMatrix(proj)$X, c(0L, 0L, 1L, 1L, NA_integer_))
})

test_that("prep_as_condition() aligns presence by ID and includes missing compounds", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  cmp <- compounds(proj)
  data <- data.frame(Name = rev(cmp$name[-1]), `R1-X` = 1, check.names = FALSE)
  proj <- prep_binarize(proj, data)
  proj <- prep_as_condition(proj, "subset", cmp$id[1:2])
  for (result in list(matrixRaw(proj), binarizedMatrix(proj))) {
    expect_setequal(result$compound_id, cmp$id)
    expect_equal(result$subset, as.integer(result$compound_id %in% cmp$id[1:2]))
    expect_true(is.na(result$X[result$compound_id == cmp$id[1]]))
    expect_true(all(result$X[result$compound_id != cmp$id[1]] == 1))
  }
})

test_that("prep_binarize(min_replicates = 2) calls presence by replicate consistency, not by the quartile threshold", {
  proj <- .test_project()
  dat <- data.frame(
    Name = c("all3", "two", "one", "none", "tiny2", "hi1", "hi2"),
    `R1-A` = c(5, 4, 3, 0, 0.001, 10, 12), `R2-A` = c(6, 0, 0, 0, 0.002, 10, 12), `R3-A` = c(7, 2, 0, 0, 0, 10, 12),
    check.names = FALSE
  )
  ## default Q1 rule: the very low-abundance compound fails the quartile threshold
  q1 <- binarizedMatrix(prep_binarize(proj, dat))
  expect_equal(q1$A[q1$compound_id == "tiny2"], 0L)

  b <- binarizedMatrix(prep_binarize(proj, dat, min_replicates = 2))
  got <- setNames(b$A, b$compound_id)
  expect_equal(got[["all3"]], 1L)   # 3 of 3 detected
  expect_equal(got[["two"]], 1L)    # 2 of 3 detected
  expect_equal(got[["one"]], 0L)    # 2 of 3 are zero
  expect_equal(got[["none"]], 0L)   # 3 of 3 are zero
  expect_equal(got[["tiny2"]], 1L)  # low abundance but detected in 2 replicates: no abundance threshold

  ## min_replicates = 3 requires every replicate
  b3 <- binarizedMatrix(prep_binarize(proj, dat, min_replicates = 3))
  expect_equal(setNames(b3$A, b3$compound_id)[["two"]], 0L)

  ## the rule is logged
  log <- projectLog(prep_binarize(proj, dat, min_replicates = 2))
  expect_true(any(grepl("detected_in_at_least_2_replicates", log$message)))
})

test_that("prep_binarize() validates min_replicates", {
  proj <- .test_project()
  dat <- data.frame(Name = c("a", "b"), `R1-A` = c(1, 0), `R2-A` = c(1, 0), check.names = FALSE)
  expect_error(prep_binarize(proj, dat, min_replicates = 3), "only 2 replicate")
  expect_error(prep_binarize(proj, dat, min_replicates = 0))
  expect_error(prep_binarize(proj, dat, min_replicates = 1.5))
})

test_that("prep_binarize() can be re-run on a project that already has an extra condition", {
  proj <- .test_project()
  dat <- data.frame(Name = c("a", "b", "c"), `R1-A` = c(5, 0, 3), `R2-A` = c(6, 0, 2), check.names = FALSE)
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_as_condition(proj, condition = "EXTRA")
  expect_true("EXTRA" %in% names(binarizedMatrix(proj)))

  ## used to abort with "matrix_raw and binarized must have the same condition columns"
  proj2 <- prep_binarize(proj, dat, min_replicates = 2)
  expect_setequal(names(binarizedMatrix(proj2))[-1], "A")
  expect_setequal(names(matrixRaw(proj2))[-1], "A")
})
