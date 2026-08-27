## Shared fixtures for tests. testthat automatically sources every
## helper-*.R file in this directory before running tests.

.test_single_compound <- function() {
  utils::read.csv(system.file("extdata", "input_single_compound.csv", package = "patliR"))
}

.test_compound_list <- function() {
  utils::read.csv(system.file("extdata", "input_compound_list.csv", package = "patliR"))
}

.test_abundance_matrix <- function() {
  utils::read.csv(
    system.file("extdata", "input_abundance_matrix.csv", package = "patliR"),
    check.names = FALSE
  )
}

.test_project <- function() {
  patliR_project(tempfile("patliR_test_"))
}

## Shared setup for the network_* family: compounds + binarized matrix +
## imported targets + network_build() already run. Used across
## test-network_centrality.R, test-network_hub_penalty.R, and later
## network_* test files.
.network_stats_test_setup <- function() {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )
  network_build(proj)
}
