test_that("targets_disease_filter() requires a non-empty disease argument", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(targets_disease_filter(proj), "disease")
  expect_error(targets_disease_filter(proj, disease = NULL), "disease")
  expect_error(targets_disease_filter(proj, disease = ""), "disease")
  expect_error(targets_disease_filter(proj, disease = c("a", "b")), "disease")
})

test_that("targets_disease_filter() requires targets_imported to exist first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  expect_error(
    targets_disease_filter(proj, disease = "EFO_0000537"),
    "targets_imported"
  )
})

test_that("targets_disease_filter() rejects an unsupported source", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(
    targets_disease_filter(proj, disease = "EFO_0000537", source = "genecards"),
    "arg"
  )
})

test_that("targets_disease_filter() validates min_score is a single number in [0, 1]", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(targets_disease_filter(proj, disease = "EFO_0000537", min_score = 1.5), "min_score")
  expect_error(targets_disease_filter(proj, disease = "EFO_0000537", min_score = -0.1), "min_score")
  expect_error(targets_disease_filter(proj, disease = "EFO_0000537", min_score = c(0.1, 0.2)), "min_score")
})

## Live-network tests (actual Open Targets GraphQL calls) are intentionally
## left out of the automated suite, same convention as the rest of this
## file -- targets_disease_profile() was exercised against the real API by
## hand (both modes) before being committed; see NEWS.md.

test_that("targets_disease_profile() requires targets_imported to exist first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")

  expect_error(targets_disease_profile(proj), "targets_imported")
})

test_that("targets_disease_profile() validates top_n_diseases", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(targets_disease_profile(proj, top_n_diseases = 0), "top_n_diseases")
  expect_error(targets_disease_profile(proj, top_n_diseases = c(1, 2)), "top_n_diseases")
})

test_that("targets_disease_profile() rejects a non-single-string disease argument", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(targets_disease_profile(proj, disease = ""), "disease")
  expect_error(targets_disease_profile(proj, disease = c("a", "b")), "disease")
})

test_that("targets_disease_profile() validates min_score is a single number in [0, 1]", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import(
    proj,
    system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
    platform = "superpred"
  )

  expect_error(targets_disease_profile(proj, min_score = 1.5), "min_score")
  expect_error(targets_disease_profile(proj, min_score = -0.1), "min_score")
})

test_that("targets_disease_profile() writes an empty result when no associations survive", {
  top <- data.frame(
    disease_id = "EFO_TEST", disease_name = "test disease",
    association_score = 0.4, evidence = "genetic_association=0.4", rank = 1L
  )
  testthat::local_mocked_bindings(
    .open_targets_map_id = function(...) "ENSG1",
    .open_targets_resolve_disease = function(...) "EFO_TEST",
    .open_targets_target_top_diseases = function(...) top,
    .package = "patliR"
  )
  for (disease in list(NULL, "EFO_TEST")) {
    proj <- .test_project()
    patliRResults(proj, "targets_imported") <- data.frame(
      compound_id = c("c1", "c2"), uniprot_id = c("P1", "P1")
    )
    proj <- targets_disease_profile(proj, disease = disease, min_score = 0.4)
    before <- patliRResults(proj, "targets_disease_profile")
    expect_equal(nrow(before), 2L)
    expect_true(all(before$association_score == 0.4))

    ## Removing the final association must replace the previous result,
    ## with a zero-row table that retains the complete output schema.
    proj <- targets_disease_profile(proj, disease = disease, min_score = 0.5)
    out <- patliRResults(proj, "targets_disease_profile")
    expect_equal(nrow(out), 0L)
    expect_named(out, names(before))
    saved <- utils::read.csv(file.path(projectDir(proj), "results", "targets_disease_profile.csv"))
    expect_equal(nrow(saved), 0L)
    expect_named(saved, names(before))
  }

  top <- NULL
  for (disease in list(NULL, "EFO_TEST")) {
    proj <- .test_project()
    patliRResults(proj, "targets_imported") <- data.frame(compound_id = "c1", uniprot_id = "P1")
    proj <- targets_disease_profile(proj, disease = disease)
    expect_equal(nrow(patliRResults(proj, "targets_disease_profile")), 0L)
    expect_true(any(grepl("no_association", projectLog(proj)$message)))
  }
})
