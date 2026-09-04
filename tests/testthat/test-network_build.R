.network_test_setup <- function() {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )
  proj
}

test_that("network_build() requires binarizedMatrix() to be populated first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )
  expect_error(network_build(proj), "binarized")
})

test_that("network_build() requires targets_imported to exist first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  expect_error(network_build(proj), "targets_imported")
})

test_that("network_build() rejects unimplemented target_source values with a clear error", {
  proj <- .network_test_setup()
  expect_error(network_build(proj, target_source = "consensus"), "not implemented yet")
  expect_error(network_build(proj, target_source = "bipartite"), "not implemented yet")
})

test_that("network_build() errors on an unknown condition name", {
  proj <- .network_test_setup()
  expect_error(network_build(proj, condition = "NOT-A-REAL-CONDITION"), "Unknown condition")
})

test_that("network_build() builds one network_edges row set per condition, restricted to present compounds", {
  proj <- .network_test_setup()
  proj <- network_build(proj)

  edges <- patliRResults(proj, "network_edges")
  expect_true(all(c("condition", "compound_id", "uniprot_id", "weight", "disease_association_score") %in% names(edges)))

  bin <- binarizedMatrix(proj)
  all_conditions <- setdiff(names(bin), "compound_id")
  expect_setequal(unique(edges$condition), all_conditions)

  ## Apigenin (CID 5280443) is only present (1) in FLO-ET/FLO-AQ in the fixture matrix
  apigenin_id <- compounds(proj)$id[compounds(proj)$pubchem_id == "5280443"]
  cond_with_apigenin <- edges$condition[edges$compound_id == apigenin_id]
  expect_true(all(cond_with_apigenin %in% c("FLO-ET", "FLO-AQ")))
})

test_that("network_build() can rebuild a single condition without touching the others", {
  proj <- .network_test_setup()
  proj <- network_build(proj)
  before <- patliRResults(proj, "network_edges")

  proj <- network_build(proj, condition = "FLO-ET")
  after <- patliRResults(proj, "network_edges")

  expect_setequal(unique(after$condition), unique(before$condition))

  row_key <- function(df) {
    df <- df[df$condition != "FLO-ET", , drop = FALSE]
    sort(apply(df, 1, paste, collapse = "|"))
  }
  expect_identical(row_key(after), row_key(before))
})

test_that("network_build() min_score drops low-probability edges", {
  proj <- .network_test_setup()
  proj_all <- network_build(proj, min_score = NULL)
  proj_strict <- network_build(proj, min_score = 0.99)

  edges_all <- patliRResults(proj_all, "network_edges")
  edges_strict <- patliRResults(proj_strict, "network_edges")
  expect_true(nrow(edges_strict) <= nrow(edges_all))
  expect_true(all(edges_strict$weight >= 0.99))
})

test_that("network_build() rejects a malformed min_score", {
  proj <- .network_test_setup()
  expect_error(network_build(proj, min_score = 1.5), "min_score")
  expect_error(network_build(proj, min_score = -1), "min_score")
})

test_that(".network_graph() returns an igraph object with compound and target node types", {
  proj <- .network_test_setup()
  proj <- network_build(proj)

  g <- patliR:::.network_graph(proj, "FLO-ET")
  expect_true(igraph::is_igraph(g))
  expect_true("type" %in% igraph::vertex_attr_names(g))
})

test_that(".network_graph() rebuilds from network_edges when the cache file is missing", {
  proj <- .network_test_setup()
  proj <- network_build(proj)

  cache_path <- patliR:::.network_cache_path(proj, "FLO-ET")
  expect_true(file.exists(cache_path))
  file.remove(cache_path)
  expect_false(file.exists(cache_path))

  g <- patliR:::.network_graph(proj, "FLO-ET")
  expect_true(igraph::is_igraph(g))
  expect_true(file.exists(cache_path)) # repopulated
})

test_that(".network_graph() errors clearly for a condition that was never built", {
  proj <- .network_test_setup()
  proj <- network_build(proj, condition = "FLO-ET")
  expect_error(patliR:::.network_graph(proj, "LEA-ET"), "network_build")
})

test_that(".network_upsert(touched_keys=) drops a recomputed slot even when it now yields zero rows", {
  proj <- .test_project()
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = c("X", "X", "Y"), value = c(1, 2, 3), stringsAsFactors = FALSE
  )

  ## a rerun for condition X that produced nothing must remove X's old
  ## rows and leave Y's intact -- headline bug #2 / spec 2.5
  out <- patliR:::.network_upsert(
    proj, "demo_slot",
    new_rows = data.frame(condition = character(0), value = double(0), stringsAsFactors = FALSE),
    key_cols = "condition",
    touched_keys = data.frame(condition = "X", stringsAsFactors = FALSE)
  )
  expect_equal(nrow(out), 1L)
  expect_equal(out$condition, "Y")
})

test_that(".network_upsert() without touched_keys keeps the historical derive-from-new_rows behaviour", {
  proj <- .test_project()
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = c("X", "Y"), value = c(1, 3), stringsAsFactors = FALSE
  )

  out <- patliR:::.network_upsert(
    proj, "demo_slot",
    data.frame(condition = "X", value = 9, stringsAsFactors = FALSE), "condition"
  )
  expect_equal(out$value[out$condition == "X"], 9)
  expect_equal(out$value[out$condition == "Y"], 3)

  ## zero-row new_rows with no touched_keys -> previous table returned unchanged
  out0 <- patliR:::.network_upsert(
    proj, "demo_slot",
    data.frame(condition = character(0), value = double(0), stringsAsFactors = FALSE), "condition"
  )
  expect_equal(nrow(out0), 2L)
})

test_that(".network_upsert() fails loudly, not with R's generic rbind message, on column removal", {
  proj <- .test_project()
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = "X", value = 1, stringsAsFactors = FALSE
  )
  ## new_rows has VALUE but not value -> `value` was removed -> abort
  expect_error(
    patliR:::.network_upsert(
      proj, "demo_slot",
      data.frame(condition = "Y", VALUE = 2, stringsAsFactors = FALSE), "condition"
    ),
    "column removed|not reconcilable"
  )
})

test_that(".network_upsert() aborts on a type conflict on a shared column", {
  proj <- .test_project()
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = "X", value = 1L, stringsAsFactors = FALSE
  )
  expect_error(
    patliR:::.network_upsert(
      proj, "demo_slot",
      data.frame(condition = "Y", value = "two", stringsAsFactors = FALSE), "condition"
    ),
    "[Tt]ype changed|not reconcilable"
  )
})

test_that(".network_upsert() still aborts on numeric-looking character vs number (the loosening is reverted)", {
  proj <- .test_project()
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = "X", value = c("12"), stringsAsFactors = FALSE
  )
  expect_error(
    patliR:::.network_upsert(
      proj, "demo_slot",
      data.frame(condition = "Y", value = 12, stringsAsFactors = FALSE), "condition"
    ),
    "[Tt]ype changed|not reconcilable"
  )
})

test_that("results CSV round-trip keeps by-contract character columns character", {
  reg <- patliR:::.patliR_results_colclasses
  expect_identical(unname(reg[["string_version"]]), "character")

  proj <- .test_project()
  dir.create(file.path(projectDir(proj), "results"), showWarnings = FALSE)
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = "2024", string_version = "12.0", score_threshold = 400,
    stringsAsFactors = FALSE
  )
  # write it the way the package does, then reload
  utils::write.csv(
    patliRResults(proj, "demo_slot"),
    file.path(projectDir(proj), "results", "demo_slot.csv"), row.names = FALSE
  )
  reloaded <- patliR_load(projectDir(proj))
  got <- patliRResults(reloaded, "demo_slot")
  expect_type(got$string_version, "character")
  expect_identical(got$string_version, "12.0")
  expect_type(got$condition, "character")
  expect_identical(got$condition, "2024")
  expect_true(is.numeric(got$score_threshold))  # genuine numbers stay numeric
})

test_that(".network_upsert() back-fills a purely-added column as NA on pre-existing rows (S1: legacy CSV migration)", {
  proj <- .test_project()
  ## a pre-Phase-1 table with no `disease_gene_source` / `n_overlap`
  patliRResults(proj, "demo_slot") <- data.frame(
    condition = c("X", "Y"), z = c(-1.2, 0.3), stringsAsFactors = FALSE
  )
  new_rows <- data.frame(
    condition = "X", z = -2.0,
    disease_gene_source = "disease_genes", n_overlap = 0L,
    stringsAsFactors = FALSE
  )
  expect_message(
    out <- patliR:::.network_upsert(proj, "demo_slot", new_rows, "condition",
                                    touched_keys = data.frame(condition = "X")),
    "back-filled"
  )
  expect_setequal(names(out), names(new_rows))
  ## Y's old row survives, with NA in the new columns; X is the recomputed row
  y <- out[out$condition == "Y", ]
  expect_equal(nrow(y), 1L)
  expect_true(is.na(y$disease_gene_source) && is.na(y$n_overlap))
  expect_equal(out$z[out$condition == "X"], -2.0)
  expect_type(out$n_overlap, "integer")
})
