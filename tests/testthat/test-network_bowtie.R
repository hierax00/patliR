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
