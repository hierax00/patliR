test_that("tox_local() flags PAINS alerts and never adds a hard_cutoff-style pass/fail", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- tox_local(proj, alert_sets = "pains")

  result <- patliRResults(proj, "tox_local")
  expect_true(all(c("compound_id", "alert_set", "alert_name", "smarts", "matched") %in% names(result)))
  expect_true(all(result$alert_set == "pains"))
  expect_true(is.logical(result$matched) || all(is.na(result$matched)))
  expect_false("hard_cutoff" %in% names(formals(tox_local)))
})

test_that("tox_local() re-running on a subset of compounds does not wipe out the rest of the table", {
  ## Same regression class as adme_local() -- see DEVLOG.md. tox_local()
  ## has multiple rows per compound (one per alert); the upsert must
  ## replace all of a re-run compound's rows, not just one, while leaving
  ## every other compound's rows untouched.
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- tox_local(proj, alert_sets = "pains")
  n_compounds_before <- length(unique(patliRResults(proj, "tox_local")$compound_id))
  one_id <- compounds(proj)$id[1]

  proj <- tox_local(proj, compound_ids = one_id, alert_sets = "pains")
  result <- patliRResults(proj, "tox_local")

  expect_equal(length(unique(result$compound_id)), n_compounds_before)
})

test_that("tox_local() flags Brenk alerts now that the full 105-alert set is bundled", {
  ## Was: "tox_local() errors clearly on 'brenk' since it is not bundled
  ## yet" -- stale since 2026-08-11 (see DEVLOG.md), when the complete
  ## Brenk set (inst/extdata/brenk_smarts.csv, 105 alerts, cross-validated
  ## against RDKit's own compiled FilterCatalogs.BRENK) got bundled and
  ## `alert_sets = c("pains", "brenk")` became the default. Rewritten to
  ## assert the new, real behaviour instead of the old placeholder error.
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- tox_local(proj, alert_sets = "brenk")

  result <- patliRResults(proj, "tox_local")
  expect_true(all(c("compound_id", "alert_set", "alert_name", "smarts", "matched") %in% names(result)))
  expect_true(all(result$alert_set == "brenk"))
  expect_true(is.logical(result$matched) || all(is.na(result$matched)))
})

test_that(".load_brenk_smarts() loads the complete, well-formed 105-alert Brenk set", {
  brenk <- patliR:::.load_brenk_smarts()
  expect_true(all(c("smarts", "name") %in% names(brenk)))
  expect_equal(nrow(brenk), 105)
  expect_true(all(nzchar(brenk$smarts)))
  expect_equal(length(unique(brenk$smarts)), nrow(brenk))
})

test_that("tox_local() can subset to specific compounds", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  one_id <- compounds(proj)$id[1]
  proj <- tox_local(proj, compound_ids = one_id, alert_sets = "pains")

  result <- patliRResults(proj, "tox_local")
  expect_true(all(result$compound_id == one_id))
})

test_that(".load_pains_smarts() loads a non-trivial, well-formed SMARTS table", {
  pains <- patliR:::.load_pains_smarts()
  expect_true(all(c("smarts", "name", "frequency") %in% names(pains)))
  expect_gt(nrow(pains), 50)
  expect_true(all(nzchar(pains$smarts)))
  ## `name` alone is NOT guaranteed unique in the upstream WEHI PAINS list --
  ## confirmed on the real bundled data: "amino_acridine_A" appears twice,
  ## for two genuinely different SMARTS, distinguished only by `frequency`
  ## (their original regId strings were "amino_acridine_A(1)" and
  ## "amino_acridine_A(46)"; our loader keeps the base name and the
  ## frequency in separate columns, so the name repeats). `smarts` is the
  ## real structural identity and is always unique; (name, frequency)
  ## together are also always unique. See tox_local.R for the same note.
  expect_equal(length(unique(pains$smarts)), nrow(pains))
  expect_equal(nrow(unique(pains[, c("name", "frequency")])), nrow(pains))
})

test_that("tox_import() imports ADMETlab toxicity properties and reconciles via canonical SMILES", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- tox_import(
    proj,
    system.file("extdata", "import_tox_admetlab.csv", package = "patliR"),
    platform = "admetlab"
  )

  imported <- patliRResults(proj, "tox_imported")
  expect_true(all(c("compound_id", "property", "value", "source", "import_date") %in% names(imported)))
  expect_true(all(c("ames", "dili", "herg") %in% imported$property))
  expect_true(all(!is.na(imported$compound_id)))
})

test_that("tox_import() requires an explicit column_map for platform = 'swissadme' (no bundled example)", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(
    tox_import(proj, system.file("extdata", "import_adme_swissadme.csv", package = "patliR"), platform = "swissadme"),
    "column_map"
  )
})

test_that("tox_report() always returns the fixed disclaimer note, regardless of results", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  ## No tox_* results at all yet -- should warn, not error, and still carry the note.
  report_empty <- suppressWarnings(tox_report(proj))
  expect_true(nzchar(report_empty$note))
  expect_match(report_empty$note, "pharmacological")

  proj <- tox_local(proj, alert_sets = "pains")
  proj <- tox_import(
    proj,
    system.file("extdata", "import_tox_admetlab.csv", package = "patliR"),
    platform = "admetlab"
  )
  report <- tox_report(proj)
  expect_true(all(c("compound_id", "n_pains_alerts", "n_imported_properties") %in% names(report$summary)))
  expect_identical(report$note, report_empty$note)
})
