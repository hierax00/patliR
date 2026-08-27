test_that("network_motifs() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_motifs(proj), "network_edges")
})

test_that("network_motifs() finds zero motifs with only the compound-target layer", {
  proj <- .network_stats_test_setup() ## network_build() only, no enrich/disease yet
  proj <- network_motifs(proj, condition = "FLO-ET")

  result <- patliRResults(proj, "network_motifs")
  expect_true(all(c("condition", "motif_type", "node_a", "node_b", "node_c",
                     "layer_a", "layer_b", "layer_c") %in% names(result)))
  expect_equal(nrow(result), 0) ## compound->target only, no B->C edges possible between two targets
})

test_that("network_motifs() never finds feedback loops -- structurally impossible by construction", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  ## simplify_go = FALSE: not testing GO simplification here, avoids a slow
  ## GOSemSim::godata() rebuild -- see network_enrich()'s own tests.
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)
  proj <- network_motifs(proj, condition = "FLO-ET")

  result <- patliRResults(proj, "network_motifs")
  expect_equal(sum(result$motif_type == "feedback"), 0)
})

test_that(".network_layered_graph() includes the pathway layer once network_enrich() has run", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  ## simplify_go = FALSE: not testing GO simplification here, avoids a slow
  ## GOSemSim::godata() rebuild -- see network_enrich()'s own tests.
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)

  g <- patliR:::.network_layered_graph(proj, "FLO-ET")
  expect_true(igraph::is_igraph(g))
  expect_true(igraph::is_directed(g))
  ## whether the pathway layer actually has members depends on whether any
  ## GO term came out significant for this tiny example gene set -- just
  ## check the layer vocabulary is a subset of what's expected, not that
  ## "pathway" is necessarily present.
  expect_true(all(igraph::V(g)$layer %in% c("compound", "target", "pathway", "disease")))
})

test_that("network_motifs() errors on a condition network_build() never built", {
  proj <- .network_stats_test_setup()
  expect_error(network_motifs(proj, condition = "NOT-A-REAL-CONDITION"), "not built")
})

test_that(".network_layered_graph()'s pathway_db restricts the layer to matching network_enrichment rows", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  ## Fabricate network_enrichment rows directly (no real enrichment call
  ## needed) with two dbs for the same condition/target, so the test
  ## doesn't pay for a live clusterProfiler enrichment call to check the
  ## filter itself.
  uniprot_ids <- unique(patliRResults(proj, "network_edges")$uniprot_id)
  testthat::skip_if(length(uniprot_ids) == 0, "no targets in the FLO-ET fixture")
  one_target <- uniprot_ids[1]

  enr <- data.frame(
    condition = "FLO-ET", db = c("go", "kegg"),
    ID = c("GO:0000001", "hsa00001"), Description = c("go term", "kegg pathway"),
    GeneRatio = "1/1", BgRatio = "1/1", pvalue = 0.01, p.adjust = 0.01, qvalue = 0.01,
    geneID = "1", Count = 1L, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_enrichment") <- enr

  ## Mock the UniProt -> Entrez mapping .network_target_pathway_edges() does
  ## internally, so this test exercises the pathway_db filter itself rather
  ## than depending on whether a toy fixture's UniProt ID happens to have a
  ## real org.Hs.eg.db mapping.
  testthat::local_mocked_bindings(
    bitr = function(x, ...) data.frame(UNIPROT = one_target, ENTREZID = "1", stringsAsFactors = FALSE),
    .package = "clusterProfiler"
  )

  g_all <- patliR:::.network_layered_graph(proj, "FLO-ET")
  g_kegg <- patliR:::.network_layered_graph(proj, "FLO-ET", pathway_db = "kegg")

  pathway_nodes_all <- igraph::V(g_all)$name[igraph::V(g_all)$layer == "pathway"]
  pathway_nodes_kegg <- igraph::V(g_kegg)$name[igraph::V(g_kegg)$layer == "pathway"]

  expect_setequal(pathway_nodes_all, c("GO:0000001", "hsa00001"))
  expect_setequal(pathway_nodes_kegg, "hsa00001")
})

test_that(".network_find_3node_motifs() gives identical results with n_cores = 1 and n_cores > 1", {
  testthat::skip_if_not_installed("parallel")
  ## Small synthetic directed graph with genuine feed-forward-loop
  ## triangles (A -> B, A -> C, B -> C), built directly rather than via a
  ## real project -- exercises the n_cores > 1 PSOCK-cluster code path
  ## without needing network_enrich()/a real pathway layer.
  edges <- data.frame(
    from = c("A1", "A1", "B1", "A2", "A2", "B2", "C1"),
    to   = c("B1", "C1", "C1", "B2", "C2", "C2", "A1")
  )
  vertices <- data.frame(
    name = c("A1", "B1", "C1", "A2", "B2", "C2"),
    layer = c("compound", "target", "pathway", "compound", "target", "pathway"),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = vertices)

  seq_res <- patliR:::.network_find_3node_motifs(g, n_cores = 1L)
  par_res <- patliR:::.network_find_3node_motifs(g, n_cores = 2L)

  seq_sorted <- seq_res[do.call(order, seq_res), ]
  par_sorted <- par_res[do.call(order, par_res), ]
  rownames(seq_sorted) <- NULL
  rownames(par_sorted) <- NULL

  expect_true(nrow(seq_res) > 0) ## sanity: the fixture graph actually has motifs to find
  expect_equal(seq_sorted, par_sorted)
})
