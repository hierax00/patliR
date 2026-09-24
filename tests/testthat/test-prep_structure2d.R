test_that("prep_structure2d() rejects the unimplemented engine before doing work", {
  proj <- .test_project()
  path <- file.path(projectDir(proj), "depictions")
  expect_error(prep_structure2d(proj, engine = "chemminer", out_dir = path), "not implemented")
  expect_false(dir.exists(path))
})

test_that("prep_structure2d() preserves each underlying failure in both logs", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  cmp <- compounds(proj)[1:3, ]
  compounds(proj) <- cmp
  testthat::local_mocked_bindings(
    .depict_2d = function(smiles, path, engine) {
      if (smiles == cmp$smiles[1]) stop("could not parse test SMILES")
      if (smiles == cmp$smiles[2]) stop("test depiction renderer failed")
      invisible(path)
    }, .package = "patliR"
  )
  result <- prep_structure2d(proj)
  log <- patliRResults(result, "structure2d_log")
  expect_identical(log$generated_ok, c(FALSE, FALSE, TRUE))
  expect_identical(log$failure_reason, c("could not parse test SMILES", "test depiction renderer failed", NA_character_))
  failures <- projectLog(result)
  failures <- failures[failures$step == "prep_structure2d", ]
  expect_equal(failures$id, cmp$id[1:2])
  expect_match(failures$message[1], "could not parse test SMILES")
  expect_match(failures$message[2], "test depiction renderer failed")
  disk <- utils::read.csv(file.path(projectDir(result), "results", "structure2d_log.csv"))
  expect_identical(disk$failure_reason, log$failure_reason)
})
