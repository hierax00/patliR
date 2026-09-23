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
  ## tox_local() has multiple rows per compound (one per alert); the upsert must
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
  ## yet" -- stale since 2026-08-11, when the complete
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

test_that("tox_export_smiles() mirrors adme_export_smiles(): SMILES list + row_order bridge", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  cmp <- compounds(proj)

  export <- tox_export_smiles(proj)
  expect_equal(export$smiles_text, paste(cmp$canonical_smiles, collapse = "\n"))
  expect_true(all(c("row_order", "compound_id", "name", "smiles") %in% names(export$mapping)))
  expect_equal(export$mapping$row_order, seq_len(nrow(cmp)))
  expect_equal(export$mapping$compound_id, cmp$id)
})

test_that("tox_import(mapping_file=) matches by row position, not by SMILES", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  map_base <- file.path(projectDir(proj), "tox_bridge")
  tox_export_smiles(proj, out_file = map_base)
  map_path <- paste0(map_base, "_map.csv")
  expect_true(file.exists(map_path))

  ## The bundled ADMETlab example is keyed by its own `smiles` column; feed
  ## it back through the bridge instead. Row order in the example file is
  ## the same order prep_compounds() ingested, so row i -> row_order i.
  proj <- tox_import(
    proj,
    system.file("extdata", "import_tox_admetlab.csv", package = "patliR"),
    platform = "admetlab",
    mapping_file = map_path
  )
  imported <- patliRResults(proj, "tox_imported")
  expect_true(all(!is.na(imported$compound_id)))
  expect_true(all(imported$compound_id %in% compounds(proj)$id))
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

  ## the summary is persisted as plain CSV, like every other family's output
  csv_path <- file.path(projectDir(proj), "results", "tox_report.csv")
  expect_true(file.exists(csv_path))
  expect_equal(nrow(utils::read.csv(csv_path)), nrow(report$summary))
})

test_that("tox_report() returns a typed empty summary for empty result tables", {
  proj <- .test_project()
  patliRResults(proj, "tox_local") <- data.frame(
    compound_id = character(), alert_set = character(),
    alert_name = character(), matched = logical()
  )
  report <- suppressWarnings(tox_report(proj))
  expect_s3_class(report$summary, "data.frame")
  expect_equal(nrow(report$summary), 0L)
  expect_type(report$summary$n_pains_alerts, "integer")
  expect_false(file.exists(file.path(projectDir(proj), "results", "tox_report.csv")))
})

test_that("tox_local() preserves unknown results when SMARTS matching fails", {
  skip_if_not_installed("rcdk")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  local_mocked_bindings(
    .load_alert_smarts = function(set) data.frame(name = "invalid", smarts = "[INVALID"),
    .package = "patliR"
  )
  proj <- tox_local(proj, alert_sets = "pains")
  expect_true(all(is.na(patliRResults(proj, "tox_local")$matched)))
})

test_that("tox_safetyome() removes stale rows even when mapping returns no rows", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  proj <- .test_project()
  panel <- patliR:::.load_safetyome_core_panel()
  local_mocked_bindings(
    bitr = function(geneID, ...) data.frame(
      UNIPROT = geneID[geneID == "mapped"],
      SYMBOL = rep(panel$gene[1], sum(geneID == "mapped"))
    ),
    .package = "clusterProfiler"
  )
  patliRResults(proj, "targets_imported") <- data.frame(
    compound_id = c("A", "B"), uniprot_id = "mapped"
  )
  proj <- tox_safetyome(proj)
  before <- patliRResults(proj, "tox_safetyome")
  patliRResults(proj, "targets_imported") <- data.frame(
    compound_id = c("A", "B"), uniprot_id = c("unmapped", "mapped")
  )
  proj <- tox_safetyome(proj, compound_ids = "A")
  expect_equal(patliRResults(proj, "tox_safetyome"), before[before$compound_id == "B", ])
  patliRResults(proj, "targets_imported") <- data.frame(compound_id = "A", uniprot_id = "unmapped")
  proj <- tox_safetyome(proj)
  expect_equal(nrow(patliRResults(proj, "tox_safetyome")), 0L)
  expect_true(any(grepl("tox_safetyome_unmapped", projectLog(proj)$message)))
})
