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
  expect_true(all(c("condition", "compound_id", "uniprot_id", "weight") %in% names(edges)))

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

test_that("network_build() no longer emits disease_association_score (breaking schema change, spec 1.1 [MUST] #1)", {
  proj <- .network_test_setup()
  proj <- network_build(proj)
  edges <- patliRResults(proj, "network_edges")
  expect_false("disease_association_score" %in% names(edges))
  expect_true(all(c("condition", "compound_id", "uniprot_id", "weight") %in% names(edges)))
})

test_that(".cache_key_slug() no longer collides on punctuation-only differences (e.g. 'FLO-ET' vs 'FLO_ET')", {
  s1 <- patliR:::.cache_key_slug("FLO-ET")
  s2 <- patliR:::.cache_key_slug("FLO_ET")
  expect_false(identical(s1, s2))
  ## same string in -> same slug out (deterministic, needed for the cache to
  ## actually be reused across calls)
  expect_identical(patliR:::.cache_key_slug("FLO-ET"), s1)
})

test_that(".network_build_igraph() aborts with a clear message when a compound_id collides with a uniprot_id (bipartite invariant)", {
  ## "X" is both a compound_id (edge 2) and a uniprot_id (edge 1's target)
  ## -- it is typed as a compound (compound_ids includes "X"), so edge 1
  ## (compound_id="C1", uniprot_id="X") becomes an intra-mode
  ## compound-compound edge.
  edges <- data.frame(
    compound_id = c("C1", "X"), uniprot_id = c("X", "P12345"), weight = c(0.9, 0.5),
    stringsAsFactors = FALSE
  )
  expect_error(patliR:::.network_build_igraph(c("C1", "X"), edges), "not bipartite")
  expect_error(patliR:::.network_build_igraph(c("C1", "X"), edges), "X")
})

test_that(".network_graph() rebuilds when the cache is present but its row count disagrees with the live network_edges (hand-edited CSV removed a row)", {
  proj <- .network_test_setup()
  proj <- network_build(proj, condition = "FLO-ET")

  g1 <- patliR:::.network_graph(proj, "FLO-ET")
  n_edges_before <- igraph::ecount(g1)
  testthat::skip_if(n_edges_before == 0, "no edges in the FLO-ET fixture")

  ## Simulate a hand-edit of results/network_edges.csv removing one row for
  ## FLO-ET -- .network_graph() reads patliRResults() first (falling back
  ## to the CSV only when that slot is NULL), so setting it directly here
  ## exercises the same "live rows disagree with the cache" path.
  edges <- patliRResults(proj, "network_edges")
  drop_idx <- which(edges$condition == "FLO-ET")[1]
  edited <- edges[-drop_idx, , drop = FALSE]
  patliRResults(proj, "network_edges") <- edited

  g2 <- patliR:::.network_graph(proj, "FLO-ET")
  expect_lt(igraph::ecount(g2), n_edges_before) ## rebuilt from the edited edge set, not served stale
})

test_that(".network_graph() rebuilds when the cache's edge CONTENT disagrees with network_edges even at the same row count (checksum catches it)", {
  proj <- .network_test_setup()
  proj <- network_build(proj, condition = "FLO-ET")

  g1 <- patliR:::.network_graph(proj, "FLO-ET")
  checksum_before <- attr(g1, "edges_checksum")

  edges <- patliRResults(proj, "network_edges")
  idx <- which(edges$condition == "FLO-ET")[1]
  edges$weight[idx] <- edges$weight[idx] + 0.0001234 ## same nrow, different content
  patliRResults(proj, "network_edges") <- edges

  g2 <- patliR:::.network_graph(proj, "FLO-ET")
  expect_equal(igraph::ecount(g2), igraph::ecount(g1)) ## same edge count
  expect_false(identical(attr(g2, "edges_checksum"), checksum_before)) ## but rebuilt, not served stale
})

test_that(".network_graph() rebuilds when binarizedMatrix()'s presence changed for the condition even though network_edges itself is untouched (regression: only edges_nrow/edges_checksum were stamped, never presence)", {
  proj <- .network_test_setup()
  proj <- network_build(proj, condition = "FLO-ET")

  g1 <- patliR:::.network_graph(proj, "FLO-ET")
  edges_before <- patliRResults(proj, "network_edges")

  ## Flip one currently-absent compound to present for FLO-ET, WITHOUT
  ## re-running network_build() -- network_edges stays byte-identical
  ## (same nrow, same checksum), but .network_build_igraph() adds every
  ## present compound as a vertex regardless of whether it has any edges,
  ## so the graph itself must change: a new isolated compound vertex.
  bm <- binarizedMatrix(proj)
  absent_idx <- which(bm[["FLO-ET"]] == 0)
  testthat::skip_if(length(absent_idx) == 0, "no currently-absent compound in this fixture for FLO-ET")
  flipped_id <- bm$compound_id[absent_idx[1]]
  bm[["FLO-ET"]][absent_idx[1]] <- 1
  binarizedMatrix(proj) <- bm

  edges_after <- patliRResults(proj, "network_edges")
  expect_equal(edges_after[edges_after$condition == "FLO-ET", ], edges_before[edges_before$condition == "FLO-ET", ]) ## network_edges genuinely untouched

  g2 <- patliR:::.network_graph(proj, "FLO-ET")
  expect_gt(igraph::vcount(g2), igraph::vcount(g1)) ## gained the newly-present, edge-less vertex
  expect_true(flipped_id %in% igraph::V(g2)$name)
  expect_false(flipped_id %in% igraph::V(g1)$name)
})

test_that(".network_edges_checksum() distinguishes edge sets with identical row count and total string length (regression: length-sum was not a real hash)", {
  ## Two single-edge tables with the SAME total character count (same nrow,
  ## same combined nchar across compound_id/uniprot_id/weight) but
  ## different content -- a length-based "checksum" collides on these;
  ## a real content hash must not.
  edges_a <- data.frame(compound_id = "C0001", uniprot_id = "P12345",
                         weight = 0.1, stringsAsFactors = FALSE)
  edges_b <- data.frame(compound_id = "C0001", uniprot_id = "P67890",
                         weight = 0.9, stringsAsFactors = FALSE)

  expect_false(identical(
    patliR:::.network_edges_checksum(edges_a),
    patliR:::.network_edges_checksum(edges_b)
  ))
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

test_that("a STRING release column named `version` survives a results CSV round-trip as character", {
  expect_identical(unname(patliR:::.patliR_results_colclasses[["version"]]), "character")
  proj <- .test_project()
  dir.create(file.path(projectDir(proj), "results"), showWarnings = FALSE)
  utils::write.csv(
    data.frame(species = 9606, version = "11.0", n_nodes = 10L, stringsAsFactors = FALSE),
    file.path(projectDir(proj), "results", "network_bowtie_summary.csv"), row.names = FALSE
  )
  got <- patliRResults(patliR_load(projectDir(proj)), "network_bowtie_summary")
  expect_identical(got$version, "11.0")
})
