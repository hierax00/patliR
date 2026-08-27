test_that("adme_local() computes descriptors and rule pass/fail for all compounds", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  adme <- patliRResults(proj, "adme_local")
  expect_equal(nrow(adme), 8)
  expect_true(all(c("mw", "logp", "hbd", "hba", "tpsa", "ro5_pass", "veber_pass", "ghose_pass") %in% names(adme)))
  expect_true(file.exists(file.path(projectDir(proj), "results", "adme_local.csv")))
})

test_that("adme_local() re-running on a subset of compounds does not wipe out the rest of the table", {
  ## Regression: patliRResults(proj, "adme_local") <- out used to be a
  ## bare overwrite -- calling adme_local() again for just one compound
  ## after already running it for the whole project silently collapsed
  ## the table to that one compound. Fixed via the same .network_upsert()
  ## pattern the network_* family already used. See DEVLOG.md.
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)
  one_id <- compounds(proj)$id[1]

  proj <- adme_local(proj, compound_ids = one_id)
  adme <- patliRResults(proj, "adme_local")

  expect_equal(nrow(adme), 8) # still every compound, not just re-run
  expect_equal(length(unique(adme$compound_id)), 8)
})

test_that("adme_local() warns and no-ops when there are no compounds", {
  proj <- .test_project()
  expect_warning(out <- adme_local(proj), "No matching compounds")
  expect_null(patliRResults(out, "adme_local"))
})

test_that("adme_import() reconciles a SwissADME export via PubChem CID", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_import(
    proj,
    system.file("extdata", "import_adme_swissadme.csv", package = "patliR"),
    platform = "swissadme"
  )
  imported <- patliRResults(proj, "adme_imported")
  expect_true(nrow(imported) > 0)
  expect_true(all(c("compound_id", "property", "value", "source", "import_date") %in% names(imported)))
  expect_true(all(!is.na(imported$compound_id)))
})

test_that("adme_filter() marks pass/fail without removing anything by default", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)
  n_before <- nrow(compounds(proj))

  proj <- adme_filter(proj, rules = c("ro5", "veber"))
  filtered <- patliRResults(proj, "adme_filtered")
  expect_true(all(c("compound_id", "rule", "pass") %in% names(filtered)))
  expect_equal(nrow(compounds(proj)), n_before) # nothing removed
})

test_that("adme_filter() with hard_cutoff and ask=FALSE removes failing compounds and logs it", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  ## Force every compound to fail 'ro5' so the cutoff path is exercised
  ## deterministically, regardless of the real descriptor values.
  adme <- patliRResults(proj, "adme_local")
  adme$ro5_pass <- FALSE
  patliRResults(proj, "adme_local") <- adme

  proj <- adme_filter(proj, rules = "ro5", hard_cutoff = TRUE, ask = FALSE)
  expect_equal(nrow(compounds(proj)), 0)
  expect_true(any(grepl("hard_cutoff removal", projectLog(proj)$message)))
})

test_that("adme_export_smiles() returns a newline-joined SMILES list and a matching mapping", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  export <- adme_export_smiles(proj)
  cmp <- compounds(proj)

  expect_true(all(c("row_order", "compound_id", "name", "smiles") %in% names(export$mapping)))
  expect_equal(nrow(export$mapping), nrow(cmp))
  expect_equal(export$mapping$row_order, seq_len(nrow(cmp)))
  expect_equal(export$mapping$compound_id, cmp$id)
  expect_equal(export$mapping$smiles, cmp$canonical_smiles)
  expect_equal(strsplit(export$smiles_text, "\n")[[1]], cmp$canonical_smiles)
})

test_that("adme_export_smiles() writes '<out_file>.txt' and '<out_file>_map.csv' when out_file is given", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  out_base <- file.path(tempdir(), "patliR_test_adme_export_smiles")
  txt_path <- paste0(out_base, ".txt")
  map_path <- paste0(out_base, "_map.csv")
  on.exit(unlink(c(txt_path, map_path)))

  export <- adme_export_smiles(proj, out_file = out_base)

  expect_true(file.exists(txt_path))
  expect_true(file.exists(map_path))
  expect_equal(readLines(txt_path), export$mapping$smiles)
  map_csv <- utils::read.csv(map_path, stringsAsFactors = FALSE)
  expect_equal(nrow(map_csv), nrow(export$mapping))
  expect_true(all(c("row_order", "compound_id", "name", "smiles") %in% names(map_csv)))
})

test_that("adme_export_smiles() supports compound_ids subsetting", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  one_id <- compounds(proj)$id[1]

  export <- adme_export_smiles(proj, compound_ids = one_id)
  expect_equal(nrow(export$mapping), 1)
  expect_equal(export$mapping$compound_id, one_id)
})

test_that("adme_export_smiles() errors when no compound matches compound_ids", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(adme_export_smiles(proj, compound_ids = "not_a_real_id"), "No matching compounds")
})

test_that("adme_export_smiles() warns and excludes compounds with no SMILES at all", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  cmp <- compounds(proj)
  cmp$smiles[1] <- NA_character_
  cmp$canonical_smiles[1] <- NA_character_
  compounds(proj) <- cmp

  expect_warning(export <- adme_export_smiles(proj), "no SMILES at all")
  expect_equal(nrow(export$mapping), nrow(cmp) - 1)
  expect_false(cmp$id[1] %in% export$mapping$compound_id)
})

test_that("adme_filter() errors on source = 'imported' (not implemented yet)", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)
  expect_error(adme_filter(proj, source = "imported"), "not implemented yet")
})
