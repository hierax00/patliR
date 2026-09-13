## plot_target_chord() (network-design-spec.md Sec 3.9) reuses
## network_bowtie()'s STRING "actions" machinery, which itself needs a real
## flat-file download on first use -- every test here mocks
## .network_stringdb() / .network_stringdb_actions_graph() the same way
## test-network_proximity.R mocks .network_stringdb() (.fake_string_db(),
## from helper-data.R), so the whole suite stays offline.

## ---- .target_chord_nodes(): pure, no proj/mocking needed -------------

test_that(".target_chord_nodes() orders by ascending degree (ties alphabetical) and caps labels by descending degree (same tie-break)", {
  targets <- c("A", "B", "C", "D", "E")
  edge_df <- data.frame(
    uniprot_a = c("A", "A", "B"),
    uniprot_b = c("B", "C", "C"),
    stringsAsFactors = FALSE
  )
  ## degrees: A=2, B=2, C=2, D=0, E=0
  nodes <- patliR:::.target_chord_nodes(targets, edge_df, top_n_labels = 2, labels = targets)

  expect_equal(nrow(nodes), 5)
  expect_equal(nodes$x, 1:5)
  expect_equal(nodes$degree, c(0, 0, 2, 2, 2))
  ## D/E tie at degree 0 -> alphabetical -> D before E, at the low-x end.
  expect_equal(nodes$uniprot_id[1:2], c("D", "E"))
  ## A/B/C tie at degree 2 -> alphabetical, at the high-x end.
  expect_equal(nodes$uniprot_id[3:5], c("A", "B", "C"))

  ## top_n_labels = 2: A/B/C tie at the highest degree -> alphabetical
  ## tie-break keeps A and B; C, D, E get no label (unlabeled tick).
  expect_equal(sum(!is.na(nodes$label)), 2)
  expect_setequal(nodes$label[!is.na(nodes$label)], c("A", "B"))
})

test_that(".target_chord_nodes(): top_n_labels = 0 labels nothing; top_n_labels >= n labels everything", {
  targets <- c("A", "B", "C")
  edge_df <- data.frame(uniprot_a = "A", uniprot_b = "B", stringsAsFactors = FALSE)

  none <- patliR:::.target_chord_nodes(targets, edge_df, top_n_labels = 0, labels = targets)
  expect_true(all(is.na(none$label)))

  all_labeled <- patliR:::.target_chord_nodes(targets, edge_df, top_n_labels = 99, labels = targets)
  expect_true(all(!is.na(all_labeled$label)))
})

test_that(".target_chord_nodes() is deterministic across repeated calls with identical input", {
  targets <- c("T5", "T1", "T3", "T2", "T4")
  edge_df <- data.frame(uniprot_a = c("T1", "T2"), uniprot_b = c("T2", "T3"), stringsAsFactors = FALSE)
  n1 <- patliR:::.target_chord_nodes(targets, edge_df, top_n_labels = 2, labels = targets)
  n2 <- patliR:::.target_chord_nodes(targets, edge_df, top_n_labels = 2, labels = targets)
  expect_identical(n1, n2)
})

## ---- .target_chord_edges(): pure, no proj/mocking needed --------------

test_that(".target_chord_edges() builds undirected, deduplicated target-target pairs and drops self-loops/unmapped nodes", {
  uni_to_string <- c(A = "s1", B = "s2", C = "s3", D = NA_character_, E = "s5")
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = c("s1", "s2", "s1"),
      to   = c("s2", "s1", "s1"), # s2->s1 reciprocal of s1->s2; s1->s1 self-loop
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )
  ed <- patliR:::.target_chord_edges(uni_to_string, g)
  ## Only A ("s1") and B ("s2") are actually vertices of g (C's "s3", D's
  ## missing STRING_id, and E's "s5" are all absent from g); s1<->s2
  ## reciprocal collapses to one row, and the s1->s1 self-loop (both ends
  ## resolve to "A") is dropped.
  expect_equal(nrow(ed), 1)
  expect_equal(ed$uniprot_a, "A")
  expect_equal(ed$uniprot_b, "B")
  expect_setequal(names(ed), c("uniprot_a", "uniprot_b", "string_a", "string_b"))
})

test_that(".target_chord_edges() returns zero rows (right columns) when fewer than 2 targets map onto the actions graph", {
  uni_to_string <- c(A = "s1", B = NA_character_)
  g <- igraph::graph_from_data_frame(data.frame(from = "s1", to = "s9", stringsAsFactors = FALSE), directed = TRUE)
  ed <- patliR:::.target_chord_edges(uni_to_string, g)
  expect_equal(nrow(ed), 0)
  expect_setequal(names(ed), c("uniprot_a", "uniprot_b", "string_a", "string_b"))
})

test_that(".target_chord_edges() returns zero rows when >= 2 targets map onto the actions graph but none share an edge", {
  uni_to_string <- c(A = "s1", B = "s2", C = "s3")
  g <- igraph::graph_from_data_frame(
    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE),
    directed = TRUE, vertices = data.frame(name = c("s1", "s2", "s3"), stringsAsFactors = FALSE)
  )
  ed <- patliR:::.target_chord_edges(uni_to_string, g)
  expect_equal(nrow(ed), 0)
})

## ---- .target_chord_actions_scores(): file-based, no proj mocking needed

test_that(".target_chord_actions_scores() parses the cached raw actions file, filters by is_directional/a_is_acting, and keeps the max score per undirected pair", {
  proj <- .test_project()
  dir <- file.path(cacheDir(proj), "stringdb")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  raw <- data.frame(
    item_id_a = c("s1", "s2", "s1", "s3"),
    item_id_b = c("s2", "s1", "s2", "s4"),
    mode = c("activation", "activation", "inhibition", "binding"),
    action = c("activation", "activation", "inhibition", ""),
    is_directional = c("t", "t", "t", "f"),
    a_is_acting = c("t", "t", "t", "f"),
    score = c(300, 900, 700, 999),
    stringsAsFactors = FALSE
  )
  raw_gz <- file.path(dir, "9606.protein.actions.v11.0.txt.gz")
  gz <- gzfile(raw_gz, "w")
  utils::write.table(raw, gz, sep = "\t", row.names = FALSE, quote = FALSE)
  close(gz)

  sc <- patliR:::.target_chord_actions_scores(proj, 9606, "11.0")
  expect_false(is.null(sc))
  ## s1/s2: two directional rows (s1->s2 score 300, s1->s2 (dup mode) score
  ## 700, s2->s1 score 900) -- all kept, undirected-normalised to the same
  ## key, the strongest (900) wins.
  row <- sc[sc$string_a == "s1" & sc$string_b == "s2", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$score, 900)
  ## s3/s4: is_directional = "f" -> excluded entirely, no row at all.
  expect_equal(nrow(sc[sc$string_a == "s3" & sc$string_b == "s4", ]), 0)
})

test_that(".target_chord_actions_scores() returns NULL when the raw actions file is not on disk", {
  proj <- .test_project()
  expect_null(patliR:::.target_chord_actions_scores(proj, 9606, "11.0"))
})

## ---- .network_target_chord_bezier(): pure, no proj/mocking needed -----

test_that(".network_target_chord_bezier() draws a baseline-y=0 arc whose peak height scales linearly with node span", {
  edge_df <- data.frame(edge_id = 1:2, x0 = c(1, 1), x1 = c(2, 10), stringsAsFactors = FALSE)
  arcs <- patliR:::.network_target_chord_bezier(edge_df, n_points = 21, height_scale = 0.5)
  expect_equal(nrow(arcs), 42)
  by_edge <- split(arcs, arcs$edge_id)
  expect_equal(by_edge[["1"]]$y[1], 0)
  expect_equal(by_edge[["1"]]$y[21], 0, tolerance = 1e-9)
  expect_equal(by_edge[["2"]]$y[1], 0)
  expect_equal(by_edge[["2"]]$y[21], 0, tolerance = 1e-9)
  ## edge 2 spans 9x edge 1's span -> its peak arc height is 9x as tall.
  expect_equal(max(by_edge[["2"]]$y) / max(by_edge[["1"]]$y), 9, tolerance = 1e-6)
})

test_that(".network_target_chord_bezier() returns zero rows for a zero-row edge_df", {
  edge_df <- data.frame(edge_id = integer(0), x0 = numeric(0), x1 = numeric(0), stringsAsFactors = FALSE)
  arcs <- patliR:::.network_target_chord_bezier(edge_df)
  expect_equal(nrow(arcs), 0)
})

## ---- plot_target_chord(): end-to-end, STRINGdb + actions graph mocked -

test_that("plot_target_chord() requires STRINGdb to be installed", {
  testthat::skip_if(requireNamespace("STRINGdb", quietly = TRUE), "STRINGdb is installed -- nothing to test here")
  proj <- .test_project()
  expect_error(plot_target_chord(proj), "STRINGdb")
})

test_that("plot_target_chord() requires network_build() to have run first", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(plot_target_chord(proj), "network_edges")
})

test_that("plot_target_chord() aborts clearly when a condition has fewer than 2 targets", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  ## Fabricate a single-target condition directly (every real fixture
  ## condition has >= 10 targets) -- same "assign patliRResults() directly"
  ## pattern test-plot_bowtie.R already uses for network_bowtie.
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "SOLO", compound_id = "C1", uniprot_id = "U1", stringsAsFactors = FALSE
  )
  expect_error(plot_target_chord(proj, condition = "SOLO", save = FALSE), "2 target")
})

test_that("plot_target_chord() aborts clearly (not opaquely) when the actions graph has zero edges among the plotted targets", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  targets <- sort(unique(edges$uniprot_id[edges$condition == cond]))
  fake <- .fake_string_db(targets)
  string_ids <- paste0("s", seq_along(targets)) # same convention .fake_string_db() uses internally
  g_no_edges <- igraph::graph_from_data_frame(
    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE),
    directed = TRUE, vertices = data.frame(name = string_ids, stringsAsFactors = FALSE)
  )
  testthat::local_mocked_bindings(
    .network_stringdb = function(...) fake,
    .network_stringdb_actions_graph = function(...) list(graph = g_no_edges),
    .package = "patliR"
  )

  expect_error(plot_target_chord(proj, condition = cond, save = FALSE), "actions")
})

test_that("plot_target_chord() happy path (STRINGdb + actions graph mocked): no error, returns a ggplot, save = FALSE leaves proj untouched", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  targets <- sort(unique(edges$uniprot_id[edges$condition == cond]))
  expect_gte(length(targets), 8) # fixture-dependent sanity check

  fake <- .fake_string_db(targets)
  ## Hub-and-spoke-ish directed graph on the mocked "s1".."s10" STRING ids
  ## (same numbering .fake_string_db() assigns, in sorted-target order):
  ## s1 hub (degree 5), s2/s3 degree 2, s4-s8 degree 1, s9/s10 isolated.
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = c("s1", "s1", "s1", "s1", "s1", "s2", "s7"),
      to   = c("s2", "s3", "s4", "s5", "s6", "s3", "s8"),
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )
  testthat::local_mocked_bindings(
    .network_stringdb = function(...) fake,
    .network_stringdb_actions_graph = function(...) list(graph = g),
    .package = "patliR"
  )

  ## The raw actions flat file was never downloaded (only the mocked graph
  ## exists) -- score is unavailable, exercising network-design-spec.md's
  ## own documented fallback ("otherwise a constant colour").
  p <- NULL
  expect_message(p <- plot_target_chord(proj, condition = cond, top_n_labels = 3, save = FALSE), "score")
  expect_s3_class(p, "ggplot")
  expect_null(attr(p, "proj"))
})

test_that("plot_target_chord() saves a PNG and logs it to target_chord_plot_log", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  targets <- sort(unique(edges$uniprot_id[edges$condition == cond]))
  fake <- .fake_string_db(targets)
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("s1", "s1"), to = c("s2", "s3"), stringsAsFactors = FALSE),
    directed = TRUE
  )
  testthat::local_mocked_bindings(
    .network_stringdb = function(...) fake,
    .network_stringdb_actions_graph = function(...) list(graph = g),
    .package = "patliR"
  )

  p <- suppressMessages(plot_target_chord(proj, condition = cond, save = TRUE))
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "target_chord_plot_log")
  expect_equal(nrow(log_df), 1)
  expect_true(file.exists(log_df$path[1]))
  expect_equal(log_df$n_targets[1], length(targets))
  expect_false(log_df$score_available[1])
})

test_that("plot_target_chord() node ordering/labels are deterministic across repeated calls with identical (mocked) input", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  targets <- sort(unique(edges$uniprot_id[edges$condition == cond]))
  fake <- .fake_string_db(targets)
  g <- igraph::graph_from_data_frame(
    data.frame(
      from = c("s1", "s1", "s2"),
      to   = c("s2", "s3", "s3"),
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )
  testthat::local_mocked_bindings(
    .network_stringdb = function(...) fake,
    .network_stringdb_actions_graph = function(...) list(graph = g),
    .package = "patliR"
  )

  p1 <- suppressMessages(plot_target_chord(proj, condition = cond, top_n_labels = 4, save = FALSE))
  p2 <- suppressMessages(plot_target_chord(proj, condition = cond, top_n_labels = 4, save = FALSE))
  ## Both calls reconstruct the node table from scratch (no cached
  ## randomness anywhere in this plot) -- the underlying node/arc layer
  ## data must match exactly.
  expect_identical(p1$layers[[2]]$data, p2$layers[[2]]$data)
})
