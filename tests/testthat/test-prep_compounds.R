test_that("prep_compounds() imports a curated compound list", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  cmp <- compounds(proj)
  expect_equal(nrow(cmp), 8) # all 8 example compounds have valid SMILES
  expect_true(all(c("id", "pubchem_id", "name", "smiles", "canonical_smiles", "source") %in% names(cmp)))
  expect_true(all(!is.na(cmp$canonical_smiles)))
  expect_true(file.exists(file.path(projectDir(proj), "01_compounds.csv")))
})

test_that("prep_compound() adds exactly one compound", {
  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem")
  expect_equal(nrow(compounds(proj)), 1)
  expect_equal(compounds(proj)$name, "Apigenin")
})

test_that("prep_compound() errors on more than one row", {
  proj <- .test_project()
  expect_error(prep_compound(proj, .test_compound_list()), "exactly one row")
})

test_that("canonical_smiles is deterministic across two equivalent SMILES for the same molecule", {
  proj <- .test_project()
  ## Same molecule (ethanol), two differently-written (non-canonical) SMILES.
  two_forms <- data.frame(
    Name = c("Ethanol A", "Ethanol B"), CAS = NA,
    PubChemCID = c("702", "702b"),
    SMILES = c("CCO", "OCC"), stringsAsFactors = FALSE
  )
  proj <- prep_compounds(proj, two_forms, identifier = "pubchem", dedup = TRUE)
  ## Both should canonicalize to the same structure key and collapse to 1 row.
  expect_equal(nrow(compounds(proj)), 1)
})

test_that("prep_compounds() keeps two enantiomers as distinct compounds under dedup = TRUE (regression: missing Stereo flavor used to collapse them)", {
  proj <- .test_project()
  ## (R)- and (S)-2-chlorobutane: same connectivity, genuinely different
  ## molecules. Before .check_structures() included the CDK `Stereo`
  ## flavor bit, both canonicalized to the same achiral SMILES and the
  ## second was silently dropped here as a "duplicate" of the first.
  enantiomers <- data.frame(
    Name = c("(R)-2-chlorobutane", "(S)-2-chlorobutane"), CAS = NA,
    PubChemCID = c("999901", "999902"),
    SMILES = c("CC[C@H](C)Cl", "CC[C@@H](C)Cl"), stringsAsFactors = FALSE
  )
  proj <- prep_compounds(proj, enantiomers, identifier = "pubchem", dedup = TRUE)

  cmp <- compounds(proj)
  expect_equal(nrow(cmp), 2)
  expect_false(identical(cmp$canonical_smiles[1], cmp$canonical_smiles[2]))
})

test_that("duplicate compounds (same canonical SMILES) are logged and dropped", {
  proj <- .test_project()
  one <- .test_single_compound()
  proj <- prep_compound(proj, one, identifier = "pubchem")
  proj <- prep_compound(proj, one, identifier = "pubchem")

  expect_equal(nrow(compounds(proj)), 1)
  log <- projectLog(proj)
  expect_true(any(grepl("duplicate", log$message)))
})

test_that("an invalid SMILES is logged and excluded, not silently dropped", {
  proj <- .test_project()
  bad <- data.frame(
    Name = "not_a_molecule", CAS = NA, PubChemCID = "999999999",
    SMILES = "this-is-not-smiles(((", stringsAsFactors = FALSE
  )
  proj <- prep_compounds(proj, bad, identifier = "pubchem")
  expect_equal(nrow(compounds(proj)), 0)
  log <- projectLog(proj)
  expect_true(any(grepl("invalid_structure", log$message)))
})

test_that("prep_compounds() requires the identifier column to exist", {
  proj <- .test_project()
  bad_data <- data.frame(SMILES = "CCO", stringsAsFactors = FALSE)
  expect_error(
    prep_compounds(proj, bad_data, identifier = "pubchem"),
    "PubChemCID"
  )
})

test_that("prep_compounds(id_col = ) actually reads the named column (codex audit regression)", {
  ## Before the fix, id_col was validated for presence and then ignored --
  ## pubchem_id/smiles were always read from the literal "PubChemCID"/
  ## "SMILES" column names, so a table using a nonstandard column name for
  ## the identifier silently produced all-NA identity for every row.
  proj <- .test_project()
  data_cid <- data.frame(Name = "Ethanol", CID = "702", SMILES = "CCO", stringsAsFactors = FALSE)
  proj <- prep_compounds(proj, data_cid, identifier = "pubchem", id_col = "CID")
  cmp <- compounds(proj)
  expect_equal(nrow(cmp), 1)
  expect_equal(cmp$pubchem_id, "702")

  proj2 <- .test_project()
  data_struct <- data.frame(Name = "Ethanol", PubChemCID = "702", structure = "CCO", stringsAsFactors = FALSE)
  proj2 <- prep_compounds(proj2, data_struct, identifier = "smiles", id_col = "structure")
  cmp2 <- compounds(proj2)
  expect_equal(nrow(cmp2), 1)
  expect_false(is.na(cmp2$canonical_smiles))
  ## the conventional PubChemCID column is still picked up as the secondary field
  expect_equal(cmp2$pubchem_id, "702")
})

test_that("prep_compounds() fetches a missing SMILES via PubChemCID even when identifier = 'smiles'", {
  ## Regression: the PubChem fetch used to only trigger when
  ## identifier == "pubchem". A Scenario B table keyed by SMILES
  ## (identifier = "smiles") that also carries a PubChemCID for some rows
  ## (real case: real_data/compound_list_Chilcuague.csv) never got its
  ## missing SMILES fetched -- rows were silently dropped as
  ## "smiles_not_resolved" even though the CID needed to fetch them was
  ## right there. Fetchability should only depend on having a PubChemCID,
  ## not on what `identifier` was set to. No network call here: mock
  ## .fetch_smiles_from_pubchem() to prove the fetch path is reached at all.
  proj <- .test_project()
  data <- data.frame(
    Name = "Ethanol", CAS = NA, PubChemCID = "702", SMILES = NA_character_,
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .fetch_smiles_from_pubchem = function(pubchem_ids, cache_dir, fetch_mode) {
      expect_equal(pubchem_ids, "702")
      "CCO"
    },
    .package = "patliR"
  )
  proj <- prep_compounds(proj, data, identifier = "smiles", on_missing_smiles = "fetch")
  cmp <- compounds(proj)
  expect_equal(nrow(cmp), 1)
  expect_equal(cmp$source, "pubchem_fetch")
  expect_false(is.na(cmp$canonical_smiles))
})

test_that("prep_compounds() never reissues a removed compound's id (codex audit regression)", {
  ## Before the fix, .next_compound_ids() derived the next id purely from
  ## compounds(proj)$id -- so removing a compound (e.g. adme_filter()'s
  ## hard_cutoff, simulated here directly) and then adding a new one reissued
  ## the removed compound's id, silently attaching any leftover downstream
  ## result rows (adme_local, targets_imported, ...) still keyed by that old
  ## id to the unrelated new molecule.
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list()[1:2, ], identifier = "pubchem")
  expect_equal(compounds(proj)$id, c("C0001", "C0002"))

  ## simulate a hard_cutoff-style removal of C0002
  compounds(proj) <- compounds(proj)[compounds(proj)$id == "C0001", , drop = FALSE]
  expect_equal(compounds(proj)$id, "C0001")

  proj <- prep_compounds(proj, .test_compound_list()[3, , drop = FALSE], identifier = "pubchem")
  new_id <- compounds(proj)$id[compounds(proj)$name == .test_compound_list()$Name[3]]
  expect_equal(new_id, "C0003") # not C0002

  ## the high-water mark survives a save/reload round trip
  proj2 <- patliR_load(projectDir(proj))
  proj2 <- prep_compounds(proj2, .test_compound_list()[4, , drop = FALSE], identifier = "pubchem")
  newer_id <- compounds(proj2)$id[compounds(proj2)$name == .test_compound_list()$Name[4]]
  expect_equal(newer_id, "C0004")
})
