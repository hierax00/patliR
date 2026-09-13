test_that(".network_bowtie_classify() finds a real core/in/out split on a synthetic directed graph", {
  ## a -> b -> c -> a (a 3-cycle, the core) ; d -> a (feeds into the core,
  ## in_component) ; c -> e (fed by the core, out_component) ; f (isolated,
  ## other).
  edges <- data.frame(
    from = c("a", "b", "c", "d", "c"),
    to   = c("b", "c", "a", "a", "e"),
    stringsAsFactors = FALSE
  )
  ## "f" has no edges at all -- included only via the explicit vertices
  ## list, to check an isolated node is classified as "other" without error.
  g <- igraph::graph_from_data_frame(
    edges, directed = TRUE,
    vertices = data.frame(name = c("a", "b", "c", "d", "e", "f"), stringsAsFactors = FALSE)
  )

  bowtie <- patliR:::.network_bowtie_classify(g)
  expect_setequal(names(bowtie)[bowtie == "core"], c("a", "b", "c"))
  expect_equal(unname(bowtie["d"]), "in_component")
  expect_equal(unname(bowtie["e"]), "out_component")
  expect_equal(unname(bowtie["f"]), "other")
})

test_that(".network_bowtie_classify() handles a graph with no cycles (no core) gracefully", {
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b"), to = c("b", "c"), stringsAsFactors = FALSE),
    directed = TRUE
  )
  bowtie <- patliR:::.network_bowtie_classify(g)
  ## every node is its own trivial SCC of size 1 here; whichever single
  ## node `which.max` picks becomes a (degenerate, single-node) "core" --
  ## the important thing is this does not error and returns a full,
  ## well-formed classification for every node.
  expect_length(bowtie, 3)
  expect_true(all(bowtie %in% c("core", "in_component", "out_component", "other")))
})

test_that("network_bowtie() requires STRINGdb to be installed", {
  testthat::skip_if(requireNamespace("STRINGdb", quietly = TRUE), "STRINGdb is installed -- nothing to test here")
  proj <- .test_project()
  expect_error(network_bowtie(proj), "STRINGdb")
})

test_that("network_bowtie() requires network_build() to have run first", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_bowtie(proj), "network_edges")
})

test_that("network_bowtie() end-to-end is not exercised automatically -- needs a real STRINGdb actions download", {
  testthat::skip_if_not_installed("STRINGdb")
  skip_on_cran()
  testthat::skip(
    "network_bowtie() end-to-end needs a real STRING 'actions' flat-file download -- verify manually against real data, not in the automated suite."
  )
})

## Spec 1.10: mocking .network_stringdb_actions_raw()/.network_stringdb_aliases_map()
## makes the rest of network_bowtie() exercisable fully offline -- no real
## STRING download, no STRINGdb object at all.

test_that(".network_stringdb_actions_graph(score_threshold=) filters low-confidence directed actions before the graph is built", {
  ## n1->n2->n3->n1 all score 900 (a 3-cycle); n3->n4 and n4->n1 score 100
  ## -- including those two low-confidence edges merges n4 into the same
  ## cycle (n1->n2->n3->n4->n1).
  raw_actions <- data.frame(
    item_id_a = c("n1", "n2", "n3", "n3", "n4"),
    item_id_b = c("n2", "n3", "n1", "n4", "n1"),
    is_directional = "t", a_is_acting = "t",
    score = c(900, 900, 900, 100, 100),
    stringsAsFactors = FALSE
  )
  proj <- .test_project()
  testthat::local_mocked_bindings(
    .network_stringdb_actions_raw = function(...) raw_actions, .package = "patliR"
  )

  g_hi <- patliR:::.network_stringdb_actions_graph(proj, 9606, "11.0", score_threshold = 400)$graph
  g_lo <- patliR:::.network_stringdb_actions_graph(proj, 9606, "11.0", score_threshold = 50)$graph
  g_none <- patliR:::.network_stringdb_actions_graph(proj, 9606, "11.0")$graph ## score_threshold = NULL default

  expect_equal(igraph::vcount(g_hi), 3L) ## only the 3 high-confidence nodes
  expect_equal(igraph::vcount(g_lo), 4L) ## n4 pulled in at the lower threshold
  expect_equal(igraph::ecount(g_none), 5L) ## unfiltered -- plot_target_chord()'s expectation, unchanged
})

test_that("network_bowtie(actions_score_threshold=) changes the core size on a fixture with mixed-confidence actions edges", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  raw_actions <- data.frame(
    item_id_a = c("n1", "n2", "n3", "n3", "n4"),
    item_id_b = c("n2", "n3", "n1", "n4", "n1"),
    is_directional = "t", a_is_acting = "t",
    score = c(900, 900, 900, 100, 100),
    stringsAsFactors = FALSE
  )
  edges <- patliRResults(proj, "network_edges")
  uniprot_ids <- unique(edges$uniprot_id[edges$condition == "FLO-ET"])
  testthat::skip_if(length(uniprot_ids) == 0, "no targets in the FLO-ET fixture")
  fake_aliases <- data.frame(
    string_protein_id = rep(c("n1", "n2", "n3", "n4"), length.out = length(uniprot_ids)),
    alias = uniprot_ids, stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    .network_stringdb_actions_raw = function(...) raw_actions,
    .network_stringdb_aliases_map = function(...) fake_aliases,
    .package = "patliR"
  )

  proj_hi <- suppressWarnings(network_bowtie(proj, condition = "FLO-ET", actions_score_threshold = 400))
  proj_lo <- suppressWarnings(network_bowtie(proj, condition = "FLO-ET", actions_score_threshold = 50))

  core_hi <- patliRResults(proj_hi, "network_bowtie_summary")$n_core
  core_lo <- patliRResults(proj_lo, "network_bowtie_summary")$n_core
  expect_equal(core_hi, 3L)
  expect_equal(core_lo, 4L)
  expect_gt(core_lo, core_hi)
})

test_that("network_bowtie() warns on a degenerate (tiny) core", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  raw_actions <- data.frame(
    item_id_a = "n1", item_id_b = "n2", is_directional = "t", a_is_acting = "t", score = 900,
    stringsAsFactors = FALSE
  )
  edges <- patliRResults(proj, "network_edges")
  uniprot_ids <- unique(edges$uniprot_id[edges$condition == "FLO-ET"])
  testthat::skip_if(length(uniprot_ids) == 0, "no targets in the FLO-ET fixture")
  fake_aliases <- data.frame(
    string_protein_id = rep(c("n1", "n2"), length.out = length(uniprot_ids)),
    alias = uniprot_ids, stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    .network_stringdb_actions_raw = function(...) raw_actions,
    .network_stringdb_aliases_map = function(...) fake_aliases,
    .package = "patliR"
  )

  ## n1<->n2 is not even a cycle (a single directed edge), so the core is
  ## degenerate (a trivial 1-node SCC) -- well under both the 50-node floor
  ## and the 1% floor.
  expect_warning(
    network_bowtie(proj, condition = "FLO-ET", actions_score_threshold = 400),
    "degenerate|noise|core"
  )
})

test_that("network_bowtie() no longer instantiates a second STRINGdb object just for UniProt -> STRING_id mapping", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  raw_actions <- data.frame(
    item_id_a = "n1", item_id_b = "n2", is_directional = "t", a_is_acting = "t", score = 900,
    stringsAsFactors = FALSE
  )
  edges <- patliRResults(proj, "network_edges")
  uniprot_ids <- unique(edges$uniprot_id[edges$condition == "FLO-ET"])
  testthat::skip_if(length(uniprot_ids) == 0, "no targets in the FLO-ET fixture")
  fake_aliases <- data.frame(
    string_protein_id = rep(c("n1", "n2"), length.out = length(uniprot_ids)),
    alias = uniprot_ids, stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    .network_stringdb_actions_raw = function(...) raw_actions,
    .network_stringdb_aliases_map = function(...) fake_aliases,
    .network_stringdb = function(...) stop("network_bowtie() should not instantiate a STRINGdb object at all any more"),
    .package = "patliR"
  )

  expect_no_error(suppressWarnings(network_bowtie(proj, condition = "FLO-ET")))
})
