test_that("patliR_export_llm() produces a text export covering compounds, results, and the log, and mutates nothing", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- tox_local(proj, alert_sets = "pains")
  before <- patliRResults(proj)
  before_log <- projectLog(proj)

  txt <- patliR_export_llm(proj)

  expect_type(txt, "character")
  expect_length(txt, 1)
  expect_match(txt, "# patliR project export", fixed = TRUE)
  expect_match(txt, "## Compounds", fixed = TRUE)
  expect_match(txt, "## Result: tox_local", fixed = TRUE)
  expect_match(txt, "## Full project log", fixed = TRUE)
  expect_match(txt, compounds(proj)$name[1], fixed = TRUE)

  expect_identical(patliRResults(proj), before)
  expect_identical(projectLog(proj), before_log)
})

test_that("patliR_export_llm() writes to out_file when given one", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  out_path <- tempfile(fileext = ".txt")
  txt <- patliR_export_llm(proj, out_file = out_path)

  expect_true(file.exists(out_path))
  expect_equal(paste(readLines(out_path), collapse = "\n"), txt)
})

test_that("patliR_export_llm() truncates long tables instead of dumping everything", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  txt <- patliR_export_llm(proj, max_rows = 2)
  expect_match(txt, "truncated at 2")
})
