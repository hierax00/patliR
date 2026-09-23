test_that("targets_import() imports a single SuperPred-style file via id_from = 'filename'", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  imported <- patliRResults(proj, "targets_imported")
  expect_true(all(c("compound_id", "uniprot_id", "probability", "confidence", "source", "import_date") %in% names(imported)))
  expect_true(all(!is.na(imported$compound_id)))
  ## Targets5280443.csv -> PubChem CID 5280443 -> Apigenin
  apigenin_id <- compounds(proj)$id[compounds(proj)$pubchem_id == "5280443"]
  expect_true(all(imported$compound_id == apigenin_id))
  expect_true("P03372" %in% imported$uniprot_id) # Estrogen receptor, first row of the fixture
})

test_that("targets_import() requires target_col/probability_col for non-superpred platforms", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  expect_error(
    targets_import(
      proj,
      system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
      platform = "swisstargetprediction"
    ),
    "target_col"
  )
})

test_that("targets_import() errors clearly when the filename does not match the expected pattern", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  bad_path <- file.path(tempdir(), "not_a_targets_file.csv")
  file.copy(
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    bad_path, overwrite = TRUE
  )

  expect_error(
    targets_import(proj, bad_path, platform = "superpred", id_from = "filename"),
    "PubChem CID"
  )
})

test_that("targets_import_batch() reads every 'Targets<cid>.csv' in a folder and skips the rest", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )

  imported <- patliRResults(proj, "targets_imported")
  expect_true(nrow(imported) > 0)
  expect_setequal(unique(imported$source), "superpred")
  ## every compound_id present should be a real id from compounds(proj)
  expect_true(all(imported$compound_id %in% compounds(proj)$id))
})

test_that("targets_import() resolves a decorated column name and parses a percent-formatted probability", {
  ## Regression coverage: targets_import() now resolves target_col/
  ## probability_col via .match_column_flexible() and parses probability/
  ## confidence via .parse_percent_column() -- neither was exercised by a
  ## real fixture file before (real SuperPred exports have used both
  ## "Probability" and "*Probability", and either a 0..1 fraction or a
  ## "96.55%" string).
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  tmp_dir <- file.path(tempdir(), "patliR_test_percent_col")
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir)
  path <- file.path(tmp_dir, "Targets5280443.csv")
  writeLines(
    c(
      "UniProt ID,*Probability,*Model accuracy",
      "P03372,96.55%,0.8",
      "P04150,50%,80%"
    ),
    path
  )

  proj <- targets_import(
    proj, path, platform = "other",
    target_col = "UniProt ID", probability_col = "Probability", confidence_col = "Model accuracy"
  )
  imported <- patliRResults(proj, "targets_imported")

  expect_equal(imported$probability, c(0.9655, 0.5))
  ## Same per-column all-or-nothing rescaling for confidence: one "80%"
  ## value forces the whole column (including the already-fractional 0.8)
  ## to be divided by 100.
  expect_equal(imported$confidence, c(0.008, 0.8))

  unlink(tmp_dir, recursive = TRUE)
})

test_that("targets_import_batch() excludes a file with the wrong columns instead of losing the whole batch", {
  ## Regression: a real SuperPred export for Chilcuague (Targets17100.csv)
  ## had no "Probability" column at all (only "Target Name", "ChEMBL-ID",
  ## "UniProt ID", "PDB Visualization", "TTD ID", "Min Activity", "Assay
  ## type") -- targets_import() correctly aborts on it, but that used to
  ## crash targets_import_batch() entirely, losing every file already
  ## imported before it in the same call.
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  ## Caffeic acid, CID 689043 -- picked distinct from the good fixture
  ## (Apigenin, CID 5280443) so the two are unambiguous in assertions.
  caffeic_id <- compounds(proj)$id[compounds(proj)$pubchem_id == "689043"]

  batch_dir <- file.path(tempdir(), "patliR_test_targets_batch")
  if (dir.exists(batch_dir)) unlink(batch_dir, recursive = TRUE)
  dir.create(batch_dir)
  file.copy(
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    file.path(batch_dir, "Targets5280443.csv")
  )
  writeLines(
    c(
      "Target Name,ChEMBL-ID,UniProt ID,PDB Visualization,TTD ID,Min Activity,Assay type",
      "Prelamin-A,CHEMBL12,P02545,6JLB,Not Available,3.5 nm,Potency"
    ),
    file.path(batch_dir, "Targets689043.csv")
  )

  proj <- targets_import_batch(proj, batch_dir, platform = "superpred")

  ## The good file's compound is still imported -- the bad file did not
  ## crash the batch and lose it.
  imported <- patliRResults(proj, "targets_imported")
  expect_true(nrow(imported) > 0)
  expect_false(caffeic_id %in% imported$compound_id)

  batch_log <- patliRResults(proj, "targets_import_batch_log")
  expect_true(all(c("path", "pubchem_id", "compound_id", "compound_name", "ok", "reason", "message") %in% names(batch_log)))
  bad_row <- batch_log[batch_log$path == "Targets689043.csv", ]
  expect_equal(nrow(bad_row), 1)
  expect_false(bad_row$ok)
  expect_equal(bad_row$reason, "import_failed")
  expect_equal(bad_row$pubchem_id, "689043")
  expect_equal(bad_row$compound_id, caffeic_id)
  expect_equal(bad_row$compound_name, "Caffeic acid")
  expect_match(bad_row$message, "Probability")

  good_row <- batch_log[batch_log$path == "Targets5280443.csv", ]
  expect_equal(nrow(good_row), 1)
  expect_true(good_row$ok)
  expect_equal(good_row$reason, "imported")

  ## The exclusion is visible in the project log by pubchem_id AND name,
  ## not just a raw id.
  log_msg <- projectLog(proj)$message[projectLog(proj)$step == "targets_import_batch"]
  expect_true(any(grepl("689043", log_msg) & grepl("Caffeic acid", log_msg)))

  unlink(batch_dir, recursive = TRUE)
})

test_that("targets_import() accepts header-only exports", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  batch_dir <- tempfile("patliR_empty_targets_")
  dir.create(batch_dir)
  path <- file.path(batch_dir, "Targets5280443.csv")
  writeLines("UniProt ID,Probability", path)
  on.exit(unlink(batch_dir, recursive = TRUE))
  proj <- targets_import(proj, path, platform = "superpred")
  result <- patliRResults(proj, "targets_imported")
  expect_equal(nrow(result), 0L)
  expect_s3_class(result$import_date, "Date")
  expect_type(result$confidence, "double")
  proj <- targets_import_batch(proj, dirname(path), platform = "superpred")
  log <- patliRResults(proj, "targets_import_batch_log")
  expect_true(log$ok[log$path == basename(path)])
})
