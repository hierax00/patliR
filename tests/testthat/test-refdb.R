test_that("refdb_build() rejects 'coconut' as a source with an informative error", {
  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem")
  expect_error(
    refdb_build(proj, sources = "coconut"),
    "not implemented yet"
  )
})

test_that("refdb_build()'s default sources do not include the unimplemented 'coconut'", {
  ## a bare refdb_build(proj) used to abort because match.arg(several.ok)
  ## returned the whole default vector, which contained "coconut".
  expect_false("coconut" %in% eval(formals(refdb_build)$sources))
})

test_that(".refdb_merge() replaces a compound's rows and drops key-duplicates (idempotent)", {
  old <- data.frame(
    compound_id = c("C0001", "C0001", "C0002"),
    source = c("pubchem", "chembl", "pubchem"),
    external_id = c("111", "CHEMBL1", "222"),
    name = c("a", "a", "b"), stringsAsFactors = FALSE
  )
  ## re-run for C0001 only, with a changed external_id
  new <- data.frame(
    compound_id = c("C0001", "C0001"),
    source = c("pubchem", "chembl"),
    external_id = c("111", "CHEMBL_NEW"),
    name = c("a", "a"), stringsAsFactors = FALSE
  )
  out <- patliR:::.refdb_merge(old, new, touched_ids = "C0001",
                               key = c("compound_id", "source"))
  ## C0002 untouched, C0001 replaced (not duplicated), 3 rows total
  expect_equal(nrow(out), 3)
  expect_equal(out$external_id[out$compound_id == "C0001" & out$source == "chembl"], "CHEMBL_NEW")
  expect_true("C0002" %in% out$compound_id)

  ## running the SAME merge twice is a no-op
  out2 <- patliR:::.refdb_merge(out, new, touched_ids = "C0001", key = c("compound_id", "source"))
  expect_equal(out2[order(out2$compound_id, out2$source), ],
               out[order(out$compound_id, out$source), ], ignore_attr = TRUE)
})

test_that(".refdb_merge() tolerates an old table written by an earlier schema", {
  old <- data.frame(compound_id = "C0002", source = "chembl", external_id = "CHEMBL9",
                    name = "keep", stringsAsFactors = FALSE)
  new <- patliR:::.empty_reference_compounds()
  new[1, ] <- list("C0001", "chembl", "CHEMBL1", "new", "AAA-BBB-N", "inchikey_exact", "CHEMBL1", "t")
  out <- patliR:::.refdb_merge(old, new, touched_ids = "C0001", key = c("compound_id", "source"))
  expect_setequal(out$compound_id, c("C0001", "C0002"))
  expect_true(all(c("inchikey", "match_type", "parent_chembl_id", "fetched_at") %in% names(out)))
  expect_true(is.na(out$inchikey[out$compound_id == "C0002"]))
})

test_that("refdb_build() warns and returns proj unchanged when there are no compounds", {
  proj <- .test_project()
  expect_warning(out <- refdb_build(proj, sources = "pubchem"), "No matching compounds")
  expect_equal(nrow(compounds(out)), 0)
})

test_that("refdb_update() requires at least one compound_id", {
  proj <- .test_project()
  expect_error(refdb_update(proj, character(0)))
})

test_that("refdb_rebuild_cache() warns when there is nothing to rebuild from", {
  proj <- .test_project()
  expect_warning(refdb_rebuild_cache(proj), "No reference database CSVs found")
})

test_that(".chembl_pick() prefers own-parent, then pref_name, then lowest id", {
  mols <- data.frame(
    query_inchikey = NA_character_,
    molecule_chembl_id = c("CHEMBL999", "CHEMBL10", "CHEMBL5"),
    pref_name = c(NA, "NAMED", "NAMED"),
    parent_chembl_id = c("CHEMBL999", "CHEMBL10", "CHEMBL7"), # row 3 is a child (salt)
    molecule_type = "Small molecule", structure_type = "MOL",
    stringsAsFactors = FALSE
  )
  ## row 3 loses (not own parent) despite lowest id; row 2 wins over row 1 (has pref_name)
  expect_equal(patliR:::.chembl_pick(mols)$molecule_chembl_id, "CHEMBL10")

  ## tie on own-parent + pref_name -> lowest numeric id
  mols2 <- data.frame(
    query_inchikey = NA_character_,
    molecule_chembl_id = c("CHEMBL10", "CHEMBL5"),
    pref_name = c("A", "B"),
    parent_chembl_id = c("CHEMBL10", "CHEMBL5"),
    molecule_type = "Small molecule", structure_type = "MOL",
    stringsAsFactors = FALSE
  )
  expect_equal(patliR:::.chembl_pick(mols2)$molecule_chembl_id, "CHEMBL5")
  expect_null(patliR:::.chembl_pick(patliR:::.empty_chembl_molecule_df()))
})

test_that(".chembl_resolve() walks CID -> InChIKey -> exact ChEMBL match", {
  testthat::local_mocked_bindings(
    .refdb_get_json = function(url, query = NULL) {
      if (grepl("pubchem", url)) {
        return(list(PropertyTable = list(Properties = list(
          list(CID = 8181L, Title = "Methyl Palmitate",
               InChIKey = "FLIACVVOZYBSBS-UHFFFAOYSA-N",
               ConnectivitySMILES = "CCCCCCCCCCCCCCCC(=O)OC")
        ))))
      }
      ## ChEMBL exact __in
      expect_true(grepl("FLIACVVOZYBSBS-UHFFFAOYSA-N", query[["molecule_structures__standard_inchi_key__in"]], fixed = TRUE))
      list(
        molecules = list(list(
          molecule_chembl_id = "CHEMBL335125",
          pref_name = "HEXADECANOIC ACID METHYL ESTER",
          molecule_hierarchy = list(parent_chembl_id = "CHEMBL335125"),
          molecule_type = "Small molecule", structure_type = "MOL",
          molecule_structures = list(standard_inchi_key = "FLIACVVOZYBSBS-UHFFFAOYSA-N")
        )),
        page_meta = list(`next` = NULL)
      )
    },
    .package = "patliR"
  )
  out <- patliR:::.chembl_resolve(smiles = "CCCCCCCCCCCCCCCC(=O)OC", cids = "8181")
  expect_equal(out$chembl_id, "CHEMBL335125")
  expect_equal(out$match_type, "inchikey_exact")
  expect_equal(out$inchikey, "FLIACVVOZYBSBS-UHFFFAOYSA-N")
})

test_that(".chembl_resolve() falls back to skeleton then flexmatch, flagging the weaker match", {
  testthat::local_mocked_bindings(
    .refdb_get_json = function(url, query = NULL) {
      if (grepl("pubchem", url)) {
        return(list(PropertyTable = list(Properties = list(
          list(CID = 1L, Title = "X", InChIKey = "ABCDEFGHIJKLMN-QRSTUVWXYZ-N",
               ConnectivitySMILES = "CCO")
        ))))
      }
      if (!is.null(query[["molecule_structures__standard_inchi_key__in"]])) {
        return(list(molecules = list(), page_meta = list(`next` = NULL))) # no exact hit
      }
      if (!is.null(query[["molecule_structures__standard_inchi_key__startswith"]])) {
        return(list(molecules = list(list(
          molecule_chembl_id = "CHEMBL777", pref_name = "SKEL",
          molecule_hierarchy = list(parent_chembl_id = "CHEMBL777"),
          molecule_structures = list(standard_inchi_key = "ABCDEFGHIJKLMN-DIFFERENT-N")
        )), page_meta = list(`next` = NULL)))
      }
      list(molecules = list(), page_meta = list(`next` = NULL))
    },
    .package = "patliR"
  )
  out <- patliR:::.chembl_resolve(smiles = "CCO", cids = "1")
  expect_equal(out$chembl_id, "CHEMBL777")
  expect_equal(out$match_type, "inchikey_skeleton")
})

test_that(".chembl_bioactivity() paginates, drops non-target rows, keeps relation/validity columns", {
  page1 <- list(
    activities = list(
      list(target_chembl_id = "CHEMBL218", target_pref_name = "Cannabinoid receptor 1",
           target_organism = "Homo sapiens", standard_type = "IC50", standard_relation = ">",
           standard_value = "100000", standard_units = "nM", pchembl_value = NULL,
           assay_chembl_id = "CHEMBL1", assay_type = "B", data_validity_comment = NULL),
      list(target_chembl_id = NA, target_pref_name = "No relevant target",
           standard_type = "LogP", standard_relation = "=", standard_value = "3.1",
           standard_units = NULL, assay_chembl_id = "CHEMBL2", assay_type = "A")
    ),
    page_meta = list(`next` = "/chembl/api/data/activity.json?limit=1000&offset=1000&molecule_chembl_id=CHEMBL335125")
  )
  page2 <- list(
    activities = list(
      list(target_chembl_id = "CHEMBL253", target_pref_name = "Cannabinoid receptor 2",
           target_organism = "Homo sapiens", standard_type = "Ki", standard_relation = "=",
           standard_value = "50", standard_units = "nM", pchembl_value = "7.3",
           assay_chembl_id = "CHEMBL9", assay_type = "B", data_validity_comment = "Outside typical range")
    ),
    page_meta = list(`next` = NULL)
  )
  calls <- 0
  testthat::local_mocked_bindings(
    .refdb_get_json = function(url, query = NULL) {
      calls <<- calls + 1
      if (calls == 1) page1 else page2
    },
    .package = "patliR"
  )
  out <- patliR:::.chembl_bioactivity("CHEMBL335125")
  expect_equal(calls, 2)                       # followed page_meta$next
  expect_equal(nrow(out), 2)                   # non-target "No relevant target" row dropped
  expect_true(all(c("standard_relation", "pchembl_value", "target_organism",
                    "assay_type", "data_validity_comment") %in% names(out)))
  expect_equal(out$standard_relation[out$target_chembl_id == "CHEMBL218"], ">")
  expect_equal(out$data_validity_comment[out$target_chembl_id == "CHEMBL253"], "Outside typical range")
})

test_that("refdb_build() writes both CSVs, the new schema, and a refdb_build log row", {
  proj <- .test_project()
  proj <- prep_compound(proj, .test_single_compound(), identifier = "pubchem")

  testthat::local_mocked_bindings(
    .pubchem_lookup = function(cid = NA_character_, smiles = NA_character_) {
      list(cid = "999", name = "Demo", inchikey = "AAAAAAAAAAAAAA-BBBBBBBBBB-N",
           connectivity_smiles = "CCO")
    },
    .chembl_resolve = function(smiles, cids) {
      data.frame(row = seq_along(smiles), inchikey = "AAAAAAAAAAAAAA-BBBBBBBBBB-N",
                 chembl_id = "CHEMBL999", parent_chembl_id = "CHEMBL999",
                 pref_name = "Demurol", match_type = "inchikey_exact",
                 stringsAsFactors = FALSE)
    },
    .chembl_bioactivity = function(chembl_id) {
      data.frame(target_chembl_id = "CHEMBL218", target_name = "CB1",
                 standard_type = "IC50", standard_value = 10, standard_units = "nM",
                 assay_chembl_id = "CHEMBL1", standard_relation = "=",
                 pchembl_value = 8, target_organism = "Homo sapiens",
                 assay_type = "B", data_validity_comment = NA_character_,
                 stringsAsFactors = FALSE)
    },
    .package = "patliR"
  )

  proj <- refdb_build(proj, sources = c("pubchem", "chembl"))
  rc <- patliRResults(proj, "reference_compounds")
  rb <- patliRResults(proj, "reference_bioactivity")
  expect_true(all(c("inchikey", "match_type", "parent_chembl_id", "fetched_at") %in% names(rc)))
  expect_true("inchikey_exact" %in% rc$match_type)
  expect_true(all(c("standard_relation", "pchembl_value", "target_organism") %in% names(rb)))
  expect_true(file.exists(file.path(projectDir(proj), "results", "reference_compounds.csv")))
  expect_true(file.exists(file.path(projectDir(proj), "results", "reference_bioactivity.csv")))
  log_df <- projectLog(proj)
  expect_true(any(log_df$step == "refdb_build"))
})

## Live-network tests (actual PubChem/ChEMBL calls) are intentionally left
## out of the automated suite -- see TESTING_GUIDE.Rmd for how we will
## exercise refdb_build()/refdb_update() together against the real APIs.
## The identity chain (.chembl_resolve) is a deterministic InChIKey/SMILES
## match, not a similarity search: for a real natural-product compound with
## no ChEMBL entry the expected outcome is chembl_id = NA and no
## bioactivity rows, logged by warn_and_cache, never a wrong molecule.
