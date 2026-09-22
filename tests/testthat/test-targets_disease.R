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
