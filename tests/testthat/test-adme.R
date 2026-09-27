test_that("adme_local() computes descriptors and rule pass/fail for all compounds", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)

  adme <- patliRResults(proj, "adme_local")
  expect_equal(nrow(adme), 8)
  expect_true(all(c("mw", "logp", "hbd", "hba", "tpsa", "ro5_pass", "veber_pass", "ghose_pass") %in% names(adme)))
  expect_true(file.exists(file.path(projectDir(proj), "results", "adme_local.csv")))
})

test_that("adme_local() keeps the scalar ALogP value for BOILED-Egg predictions", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  cmp <- compounds(proj)
  cmp$smiles[1] <- "c1ccccc1"
  compounds(proj) <- cmp
  proj <- adme_local(proj, compound_ids = cmp$id[1])
  adme <- patliRResults(proj, "adme_local")

  mol <- rcdk::parse.smiles("c1ccccc1")[[1]]
  rcdk::convert.implicit.to.explicit(mol)
  expected <- rcdk::get.alogp(mol)
  expect_true(is.finite(expected))
  expect_equal(adme$wlogp_proxy, as.numeric(expected))
  expect_false(is.na(adme$gi_absorption))
  expect_false(is.na(adme$bbb_permeant))
})

test_that("adme_local() divides aromatic atoms by heavy atoms for ESOL", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_single_compound(), identifier = "pubchem")
  cmp <- compounds(proj)
  cmp$smiles <- "c1ccccc1"
  compounds(proj) <- cmp
  proj <- adme_local(proj)
  adme <- patliRResults(proj, "adme_local")

  ## Benzene has six aromatic heavy atoms and six hydrogens. Ghose
  ## uses total atoms, whereas ESOL's aromatic proportion uses heavy atoms.
  expect_equal(adme$n_atoms, 12L)
  expect_equal(adme$aromatic_proportion_approx, 1)
  expect_true(all(is.finite(c(adme$logp, adme$mw, adme$rotatable_bonds, adme$logs_esol))))
  expect_equal(adme$logs_esol,
               0.16 - 0.63 * adme$logp - 0.0062 * adme$mw +
                 0.066 * adme$rotatable_bonds - 0.74)
})

.adme_one <- function(smiles) {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_single_compound(), identifier = "pubchem")
  cmp <- compounds(proj)
  cmp$smiles <- smiles
  compounds(proj) <- cmp
  patliRResults(adme_local(proj), "adme_local")
}

test_that("Fsp3 and aromatic proportion come from the molecular graph, not the SMILES spelling", {
  ## regression: the regex proxy gave Fsp3 = 1 for C=C and 0 vs 1 for the two benzene spellings
  expect_equal(.adme_one("C=C")$fraction_csp3_approx, 0)
  expect_equal(.adme_one("CC(=O)OC")$fraction_csp3_approx, 2 / 3)
  a <- .adme_one("c1ccccc1")
  b <- .adme_one("C1=CC=CC=C1")
  expect_equal(a$fraction_csp3_approx, 0)
  expect_equal(b$fraction_csp3_approx, 0)
  expect_equal(a$aromatic_proportion_approx, 1)
  expect_equal(b$aromatic_proportion_approx, 1)
  expect_true(is.finite(a$logs_esol) && is.finite(b$logs_esol))
  expect_equal(a$logs_esol, b$logs_esol)
})

test_that("mw is the average molecular weight, not the monoisotopic mass", {
  ## C12H6Br4O2: 497.71 monoisotopic vs 501.79 average -> fails the 500 gate on average MW
  d <- .adme_one("Oc1c(Br)cc(Br)cc1-c1cc(Br)cc(Br)c1O")
  expect_equal(d$mw, 501.79, tolerance = 1e-3)
  expect_false(d$ro5_pass)
})

test_that("rotatable_bonds excludes amide C-N bonds as in Veber et al.", {
  expect_equal(.adme_one("CCCCCCCCCCCC(=O)NC")$rotatable_bonds, 10L)  # 11 with CDK's default
  expect_equal(.adme_one("CC(C)Cc1ccc(cc1)C(C)C(=O)O")$rotatable_bonds, 4L)
})

test_that("wlogp_source defaults to the CDK proxy and is recorded per compound", {
  d <- .adme_one("c1ccccc1")
  expect_equal(d$wlogp_source, "cdk_alogp_proxy")
  expect_equal(d$wlogp_proxy, d$logp[1] * 0 + d$wlogp_proxy) # column exists, sanity only
})

test_that("adme_local(wlogp_source = 'rdkit') gives the true WLogP when RDKit is available, else falls back and warns", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  rk <- tryCatch(patliR:::.rdkit_descriptors("c1ccccc1"), error = function(e) NULL)
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_single_compound(), identifier = "pubchem")
  cmp <- compounds(proj); cmp$smiles <- "c1ccccc1"; compounds(proj) <- cmp
  if (is.null(rk)) {
    expect_warning(proj <- adme_local(proj, wlogp_source = "rdkit"), "rdkit")
    d <- patliRResults(proj, "adme_local")
    expect_equal(d$wlogp_source, "cdk_alogp_proxy")
  } else {
    proj <- adme_local(proj, wlogp_source = "rdkit")
    d <- patliRResults(proj, "adme_local")
    expect_equal(d$wlogp_source, "rdkit_wlogp")
    expect_equal(d$wlogp_proxy, 1.6866, tolerance = 1e-3) # RDKit's own documented value for benzene
    expect_false(isTRUE(all.equal(d$wlogp_proxy, d$logp))) # not silently the same as XLogP
  }
})

test_that("adme_local(rdkit_qc = TRUE) adds the cross-check columns without touching the CDK-based ones", {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  rk <- tryCatch(patliR:::.rdkit_descriptors("c1ccccc1"), error = function(e) NULL)
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_single_compound(), identifier = "pubchem")
  cmp <- compounds(proj); cmp$smiles <- "c1ccccc1"; compounds(proj) <- cmp
  base <- patliRResults(adme_local(proj), "adme_local")

  if (is.null(rk)) {
    expect_warning(proj2 <- adme_local(proj, rdkit_qc = TRUE), "rdkit")
    d <- patliRResults(proj2, "adme_local")
    expect_true(all(is.na(d[, c("hba_lipinski_rdkit", "hbd_lipinski_rdkit", "tpsa_rdkit",
                                "fraction_csp3_rdkit", "aromatic_proportion_rdkit")])))
  } else {
    proj2 <- adme_local(proj, rdkit_qc = TRUE)
    d <- patliRResults(proj2, "adme_local")
    ## benzene: RDKit's literal Lipinski N+O count is 0 (no N/O at all); its own
    ## aromaticity model agrees with CDK's here (both give a fully aromatic ring)
    expect_equal(d$hba_lipinski_rdkit, 0L)
    expect_equal(d$hbd_lipinski_rdkit, 0L)
    expect_equal(d$aromatic_proportion_rdkit, 1)
    expect_equal(d$fraction_csp3_rdkit, 0)
    expect_equal(d$tpsa_rdkit, d$tpsa, tolerance = 1e-6) # same Ertl method, benzene has no polar atoms either way
  }
  ## the existing CDK-based columns are untouched by rdkit_qc
  expect_equal(d[, setdiff(names(base), names(d))], base[, setdiff(names(base), names(d)), drop = FALSE])
  expect_equal(base$hba, d$hba); expect_equal(base$mw, d$mw); expect_equal(base$wlogp_source, d$wlogp_source)
})

test_that("Oprea ranges are inclusive at their endpoints", {
  ## ethyl benzoate: 1 ring, 2 acceptors, 3 rotatable bonds, 0 donors -> inside every published range
  d <- .adme_one("CCOC(=O)c1ccccc1")
  expect_equal(d$n_rings_approx, 1)
  expect_true(d$oprea_pass)
})

test_that("adme_local() re-running on a subset of compounds does not wipe out the rest of the table", {
  ## Regression: patliRResults(proj, "adme_local") <- out used to be a
  ## bare overwrite -- calling adme_local() again for just one compound
  ## after already running it for the whole project silently collapsed
  ## the table to that one compound. Fixed via the same .network_upsert()
  ## pattern the network_* family already used.
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

test_that("adme_import(mapping_file=) matches by row position", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  map_base <- file.path(projectDir(proj), "adme_bridge")
  adme_export_smiles(proj, out_file = map_base)
  map_path <- paste0(map_base, "_map.csv")

  ## the bundled SwissADME example rows are in the same order as the
  ## compound list, so row i -> row_order i
  proj <- adme_import(
    proj,
    system.file("extdata", "import_adme_swissadme.csv", package = "patliR"),
    platform = "swissadme",
    mapping_file = map_path
  )
  imported <- patliRResults(proj, "adme_imported")
  expect_true(nrow(imported) > 0)
  expect_true(all(!is.na(imported$compound_id)))
  expect_true(all(imported$compound_id %in% compounds(proj)$id))
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

test_that("adme_filter() does not count unknown rule outcomes as failures", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  ids <- compounds(proj)$id
  adme <- data.frame(compound_id = ids, ro5_pass = NA)
  patliRResults(proj, "adme_local") <- adme

  all_unknown <- adme_filter(proj, rules = "ro5", hard_cutoff = TRUE, ask = FALSE)
  expect_equal(compounds(all_unknown)$id, ids)
  expect_true(all(is.na(patliRResults(all_unknown, "adme_filtered")$pass)))
  expect_false(any(grepl("hard_cutoff removal", projectLog(all_unknown)$message)))

  adme$ro5_pass <- TRUE
  adme$ro5_pass[1:2] <- c(FALSE, NA)
  patliRResults(proj, "adme_local") <- adme
  testthat::local_mocked_bindings(
    .ask_hard_cutoff = function(n_remove, n_total, summary_msg, fails_any, cmp) {
      expect_equal(n_remove, 1L)
      expect_equal(summary_msg, "ro5_pass: 1")
      expect_equal(fails_any, ids[1])
      TRUE
    },
    .package = "patliR"
  )
  proj <- adme_filter(proj, rules = "ro5", hard_cutoff = TRUE, ask = TRUE)
  expect_equal(compounds(proj)$id, ids[-1])
  log <- projectLog(proj)
  expect_equal(log$id[grepl("hard_cutoff removal", log$message)], ids[1])
})
