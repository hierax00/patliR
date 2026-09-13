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

## A stand-in for the object `.network_stringdb()` returns, so
## network_proximity() can be exercised offline (the real path needs a
## tens-to-hundreds-of-MB STRING flat-file download). `uniprot_ids` are the
## accessions the test will ask about; each maps to a distinct STRING node
## `s<i>`, and the graph is a connected star+path over those nodes so
## degree binning and BFS distances are well defined.
.fake_string_db <- function(uniprot_ids) {
  ids <- unique(as.character(uniprot_ids))
  string_ids <- stats::setNames(paste0("s", seq_along(ids)), ids)
  n <- length(string_ids)
  el <- rbind(
    cbind(string_ids[rep(1L, n - 1L)], string_ids[-1L]),  # star from node 1
    cbind(string_ids[-n], string_ids[-1L])                # + a path for degree variety
  )
  g <- igraph::simplify(igraph::graph_from_edgelist(matrix(as.character(el), ncol = 2), directed = FALSE))
  list(
    get_graph = function() g,
    map = function(my_data_frame, my_data_frame_id_col_name,
                   removeUnmappedRows = FALSE, quiet = TRUE) {
      key <- as.character(my_data_frame[[my_data_frame_id_col_name]])
      my_data_frame$STRING_id <- unname(string_ids[key])
      if (isTRUE(removeUnmappedRows)) {
        my_data_frame <- my_data_frame[!is.na(my_data_frame$STRING_id), , drop = FALSE]
      }
      my_data_frame
    }
  )
}

## Save a hand-modified igraph object directly to a condition's cache slot
## so a test can exercise network_module_robustness()/etc. on a synthetic
## graph without a full network_build() rebuild. .network_graph() (since
## the network_build.R cache-staleness fix) only trusts a cached graph
## whose edges_nrow/edges_checksum attrs match the project's CURRENT
## network_edges for that condition -- an unstamped saveRDS() is treated as
## stale and silently rebuilt from network_edges, discarding the
## hand-edited graph entirely. Stamping against the *unchanged* edges
## table (rather than deriving from `g` itself) keeps the cache "fresh"
## from that check's point of view while still letting the test's `g`
## differ structurally from what network_build() would produce.
.save_cache_graph <- function(g, proj, condition) {
  edges_all <- patliRResults(proj, "network_edges")
  edges <- edges_all[edges_all$condition == condition, , drop = FALSE]
  attr(g, "edges_nrow") <- nrow(edges)
  attr(g, "edges_checksum") <- patliR:::.network_edges_checksum(edges)
  saveRDS(g, patliR:::.network_cache_path(proj, condition))
  invisible(g)
}
