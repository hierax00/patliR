## ---------------------------------------------------------------------------
## Fixtures
## ---------------------------------------------------------------------------

## A network_proximity() result frame carrying the full post-piece-12
## schema, so network_synergy() never has to fall back on a missing column.
.synergy_fake_prox <- function(compounds, condition = "FLO-ET", disease = "SOME_DISEASE",
                                z = NULL, source = "disease_genes", n_random = 1000L,
                                n_tests_in_family = NULL, p_adjusted = NULL,
                                species = 9606, string_version = "12.0",
                                score_threshold = 400) {
  n <- length(compounds)
  if (is.null(z)) z <- seq(-3, 2, length.out = n)
  if (is.null(p_adjusted)) p_adjusted <- ifelse(z < 0, 0.001, 0.5)
  if (is.null(n_tests_in_family)) n_tests_in_family <- n
  data.frame(
    condition = condition, compound_id = compounds, disease_id = disease,
    disease_gene_source = source,
    n_targets_mapped = 3L, n_disease_genes_mapped = 5L, n_overlap = 0L,
    d_observed = 2, d_random_mean = 3, d_random_sd = 1,
    z_score = z, p_empirical = 0.001, n_random = as.integer(n_random),
    seed_used = 1L,
    species = species, string_version = string_version, score_threshold = score_threshold,
    n_tests_in_family = as.integer(n_tests_in_family), p_adjusted = p_adjusted,
    stringsAsFactors = FALSE
  )
}

## Stand-ins for .network_stringdb() ($map only) and .network_string_lcc()
## ($graph only -- synergy uses nothing else). `singleton_targets` are
## UniProt accessions the map should refuse, so a compound can be forced to
## exactly one mapped target. `graph` overrides the default connected path.
.synergy_mock_string <- function(uniprot_ids, singleton_targets = NULL, graph = NULL) {
  ids <- unique(as.character(uniprot_ids))
  string_ids <- stats::setNames(paste0("s", seq_along(ids)), ids)
  g <- if (!is.null(graph)) {
    graph
  } else {
    n <- length(ids)
    el <- cbind(string_ids[-n], string_ids[-1])
    igraph::simplify(igraph::graph_from_edgelist(matrix(as.character(el), ncol = 2), directed = FALSE))
  }
  map_fn <- function(my_data_frame, my_data_frame_id_col_name,
                     removeUnmappedRows = FALSE, quiet = TRUE) {
    key <- as.character(my_data_frame[[my_data_frame_id_col_name]])
    sid <- unname(string_ids[key])
    if (!is.null(singleton_targets)) sid[key %in% singleton_targets] <- NA_character_
    my_data_frame$STRING_id <- sid
    if (isTRUE(removeUnmappedRows)) my_data_frame <- my_data_frame[!is.na(my_data_frame$STRING_id), , drop = FALSE]
    my_data_frame
  }
  list(
    stringdb = list(map = map_fn, get_graph = function() g),
    lcc = list(graph = g)
  )
}

.synergy_use_mock <- function(mock) {
  testthat::local_mocked_bindings(
    .network_stringdb    = function(...) mock$stringdb,
    .network_string_lcc  = function(...) mock$lcc,
    .package = "patliR",
    .env = parent.frame()
  )
}

.synergy_flo_compounds <- function(proj) {
  e <- patliRResults(proj, "network_edges")
  unique(e$compound_id[e$condition == "FLO-ET"])
}

## Append a fabricated single-condition network to an existing edges table.
.synergy_add_condition <- function(proj, condition, compound_targets) {
  e <- patliRResults(proj, "network_edges")
  pairs <- do.call(rbind, lapply(names(compound_targets), function(cp) {
    data.frame(compound_id = cp, uniprot_id = compound_targets[[cp]], stringsAsFactors = FALSE)
  }))
  extra <- e[rep(1, nrow(pairs)), , drop = FALSE]
  extra$condition <- condition
  extra$compound_id <- pairs$compound_id
  extra$uniprot_id <- pairs$uniprot_id
  patliRResults(proj, "network_edges") <- rbind(e, extra)
  proj
}

## ---------------------------------------------------------------------------
## Distance / separation helpers -- spec 2.1 tests 1-5 (+ pre-review 12-SE2)
## ---------------------------------------------------------------------------

test_that(".network_set_distance() singleton convention and D path", {
  g <- igraph::graph_from_literal(1 - 2, 1 - 3, 2 - 3, 4 - 5, 4 - 6, 5 - 6, 3 - 4)
  A <- c("1", "2", "3")

  expect_identical(patliR:::.network_set_distance(g, "1"), 0)          # |ids| < 2 => 0
  expect_identical(patliR:::.network_set_distance(g, character(0)), 0) # empty => 0
  expect_equal(patliR:::.network_set_distance(g, A), 1)               # all adjacent

  D <- igraph::distances(g, v = igraph::V(g)$name, to = igraph::V(g)$name, weights = NA)
  expect_equal(patliR:::.network_set_distance(g, A, D), patliR:::.network_set_distance(g, A))
})

test_that("s_AA == -.network_set_distance(g, x) exactly (pre-review 12-SE2)", {
  g <- igraph::graph_from_literal(1 - 2, 1 - 3, 2 - 3, 4 - 5, 4 - 6, 5 - 6, 3 - 4)
  for (x in list(c("1", "2", "3"), c("1", "2"), "1", c("4", "5", "6"))) {
    expect_equal(patliR:::.network_separation(g, x, x), -patliR:::.network_set_distance(g, x))
    expect_lte(patliR:::.network_separation(g, x, x), 0)
  }
})

test_that("s_AB on two 3-cliques joined by one bridge edge matches hand arithmetic", {
  ## G: cliques {1,2,3} and {4,5,6}, plus the single bridge edge 3-4.
  ##   <d_AA> = mean over A of nearest-other-A = (1+1+1)/3 = 1  ; <d_BB> = 1
  ##   A->B nearest: 1->4 = 2, 2->4 = 2, 3->4 = 1            => sum 5
  ##   B->A nearest: 4->3 = 1, 5->3 = 2, 6->3 = 2            => sum 5
  ##   <d_AB> = (5 + 5) / (3 + 3) = 10/6
  ##   s_AB   = 10/6 - (1 + 1)/2 = 2/3
  g <- igraph::graph_from_literal(1 - 2, 1 - 3, 2 - 3, 4 - 5, 4 - 6, 5 - 6, 3 - 4)
  A <- c("1", "2", "3"); B <- c("4", "5", "6")
  expect_equal(patliR:::.network_between_set_distance(g, A, B), 10 / 6, tolerance = 1e-12)
  expect_equal(patliR:::.network_separation(g, A, B), 2 / 3, tolerance = 1e-12)
})

test_that("s_AB == s_BA (symmetry)", {
  g <- igraph::graph_from_literal(1 - 2, 1 - 3, 2 - 3, 4 - 5, 4 - 6, 5 - 6, 3 - 4)
  A <- c("1", "2", "3"); B <- c("4", "5", "6")
  expect_equal(patliR:::.network_separation(g, A, B), patliR:::.network_separation(g, B, A))
  expect_equal(
    patliR:::.network_between_set_distance(g, A, B),
    patliR:::.network_between_set_distance(g, B, A)
  )
})

test_that("two disjoint components => d_AB is NA, not Inf; separation is NA", {
  g <- igraph::graph_from_literal(1 - 2, 2 - 3, 4 - 5, 5 - 6)
  expect_true(is.na(patliR:::.network_between_set_distance(g, c("1", "2"), c("4", "5"))))
  expect_true(is.na(patliR:::.network_separation(g, c("1", "2", "3"), c("4", "5", "6"))))
})

test_that(".synergy_cheng_class(): P2 iff separated & proximal_a & proximal_b, over all 8 combos", {
  for (sep in c(TRUE, FALSE)) for (pa in c(TRUE, FALSE)) for (pb in c(TRUE, FALSE)) {
    cc <- patliR:::.synergy_cheng_class(sep, pa, pb)
    expect_equal(identical(cc, "P2"), sep && pa && pb)
    expect_false(is.na(cc))
  }
  ## P1..P6 spot checks
  expect_identical(patliR:::.synergy_cheng_class(FALSE, TRUE,  TRUE),  "P1")
  expect_identical(patliR:::.synergy_cheng_class(TRUE,  TRUE,  FALSE), "P4")
  expect_identical(patliR:::.synergy_cheng_class(FALSE, FALSE, FALSE), "P5")
  expect_identical(patliR:::.synergy_cheng_class(TRUE,  FALSE, FALSE), "P6")
  ## any NA input => NA class
  expect_true(is.na(patliR:::.synergy_cheng_class(NA, TRUE, TRUE)))
  expect_true(is.na(patliR:::.synergy_cheng_class(TRUE, NA, FALSE)))
})

test_that(".synergy_proximal(): three-valued, with a per-row p_adjusted fallback", {
  expect_true(is.na(patliR:::.synergy_proximal(NA_real_, 0.01, 0.05)))
  expect_true(patliR:::.synergy_proximal(-2, 0.01, 0.05))
  expect_false(patliR:::.synergy_proximal(-2, 0.20, 0.05))
  expect_false(patliR:::.synergy_proximal(1, 0.01, 0.05))
  ## p_adjusted NA => fall back to the z-sign alone
  expect_true(patliR:::.synergy_proximal(-2, NA_real_, 0.05))
  expect_false(patliR:::.synergy_proximal(0.5, NA_real_, 0.05))
})

## ---------------------------------------------------------------------------
## network_synergy() -- error paths / regressions
## ---------------------------------------------------------------------------

test_that("network_synergy() requires network_proximity() to have run first", {
  proj <- .network_stats_test_setup()
  expect_error(
    network_synergy(proj, condition = "FLO-ET", disease = "EFO_0000000", separation = "jaccard"),
    "network_proximity"
  )
})

test_that("network_synergy() errors when network_proximity() ran for a different disease", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(cps, disease = "SOME_DISEASE")
  expect_error(
    network_synergy(proj, condition = "FLO-ET", disease = "NOT-THE-SAME", separation = "jaccard"),
    "network_proximity"
  )
})

test_that("network_synergy(pairs = 'all') does not error when a compound is absent from network_proximity (S1.9 regression)", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  skip_if(length(cps) < 3, "needs >= 3 FLO-ET compounds")
  scored <- cps[-length(cps)]
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(scored, z = seq(-2, 1, length.out = length(scored)))

  expect_no_error(
    proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE",
                            separation = "jaccard", pairs = "all")
  )
  res <- patliRResults(proj, "network_synergy")
  miss <- cps[length(cps)]
  mr <- res[res$compound_a == miss | res$compound_b == miss, , drop = FALSE]
  expect_gt(nrow(mr), 0)
  za_missing <- ifelse(mr$compound_a == miss, mr$z_score_a, mr$z_score_b)
  expect_true(all(is.na(za_missing)))
  expect_true(all(!mr$both_proximal))
  expect_true(all(is.na(mr$synergy_score)))
})

## ---------------------------------------------------------------------------
## Schema
## ---------------------------------------------------------------------------

test_that("network_synergy() output column set + order is stable across zero-row and populated conditions", {
  empty <- patliR:::.empty_network_synergy_row()
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)

  ## populated
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(cps)
  proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", separation = "jaccard", pairs = "all")
  res_pop <- patliRResults(proj, "network_synergy")
  expect_gt(nrow(res_pop), 0)
  expect_identical(names(res_pop), names(empty))

  ## zero-row: proximity for a single compound => < 2 candidates
  proj2 <- .network_stats_test_setup()
  patliRResults(proj2, "network_proximity") <- .synergy_fake_prox(cps[1])
  proj2 <- network_synergy(proj2, condition = "FLO-ET", disease = "SOME_DISEASE",
                           separation = "jaccard", pairs = "rank_top", top_n = 5)
  res_zero <- patliRResults(proj2, "network_synergy")
  expect_equal(nrow(res_zero), 0)
  expect_identical(names(res_zero), names(empty))
})

test_that(".network_upsert() migrates a legacy 12-column network_synergy.csv without aborting", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  legacy <- data.frame(
    condition = "OLD-COND", disease_id = "OLD_D", compound_a = "X", compound_b = "Y",
    target_jaccard = 0.1, complementarity = 0.9, z_score_a = -1, z_score_b = -1,
    both_proximal = TRUE, joint_closeness = 1, synergy_score = 0.9, pairs_mode = "all",
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_synergy") <- legacy
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(cps)

  expect_no_error(
    proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", separation = "jaccard", pairs = "all")
  )
  res <- patliRResults(proj, "network_synergy")
  expect_true("OLD_D" %in% res$disease_id)                       # legacy row retained
  expect_true(all(c("cheng_class", "s_ab", "singleton_a", "separation_method") %in% names(res)))
  expect_true(all(is.na(res$cheng_class[res$disease_id == "OLD_D"])))  # back-filled NA
})

## ---------------------------------------------------------------------------
## disease_gene_source gating (pre-review 12-R2)
## ---------------------------------------------------------------------------

test_that("a mixed disease_gene_source proximity table warns and does not pick the wrong z", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  prox_dg <- .synergy_fake_prox(cps, source = "disease_genes",   z = rep(-2, length(cps)))
  prox_td <- .synergy_fake_prox(cps, source = "targets_disease", z = rep(3,  length(cps)))
  ## targets_disease rows first: a buggy (no-subset) setNames() would resolve z_of[a]
  ## to the +3 rows; the correct path subsets to disease_genes and sees -2.
  patliRResults(proj, "network_proximity") <- rbind(prox_td, prox_dg)

  expect_warning(
    proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", separation = "jaccard", pairs = "all"),
    "multiple"
  )
  res <- patliRResults(proj, "network_synergy")
  expect_true(all(res$z_score_a < 0))
  expect_true(all(res$z_score_b < 0))
  expect_true(all(res$disease_gene_source == "disease_genes"))
})

test_that("disease_gene_source = 'targets_disease' warns about inherited circularity", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(cps, source = "targets_disease")
  expect_warning(
    network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE",
                    separation = "jaccard", pairs = "all", disease_gene_source = "targets_disease"),
    "circularity"
  )
})

## ---------------------------------------------------------------------------
## unreachable-alpha guard (pre-review 12-U3)
## ---------------------------------------------------------------------------

test_that("the unreachable-alpha warning fires at n_random = 100 with 30 compounds", {
  proj <- .network_stats_test_setup()
  cps <- sprintf("C%02d", 1:30)
  ct <- stats::setNames(lapply(seq_along(cps), function(i) paste0("U", c(i, i + 1))), cps)
  proj <- .synergy_add_condition(proj, "SYN", ct)
  prox <- .synergy_fake_prox(cps, condition = "SYN", n_random = 100L, n_tests_in_family = 30L)
  patliRResults(proj, "network_proximity") <- prox

  expect_warning(
    network_synergy(proj, condition = "SYN", disease = "SOME_DISEASE", separation = "jaccard", pairs = "all"),
    "cannot clear"
  )
})

## ---------------------------------------------------------------------------
## separation = "jaccard" regression (spec test 6)
## ---------------------------------------------------------------------------

test_that("separation = 'jaccard' reproduces the historical complementarity + synergy_score exactly", {
  proj <- .network_stats_test_setup()
  cps <- .synergy_flo_compounds(proj)
  z <- stats::setNames(seq(-2, 1, length.out = length(cps)), cps)
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(cps, z = unname(z))
  proj <- network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", separation = "jaccard", pairs = "all")
  res <- patliRResults(proj, "network_synergy")

  ## recompute the pre-piece-12 quantities directly
  e <- patliRResults(proj, "network_edges")
  ct <- unique(e[e$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  tset <- lapply(split(ct$uniprot_id, ct$compound_id), unique)
  for (i in seq_len(nrow(res))) {
    a <- res$compound_a[i]; b <- res$compound_b[i]
    ta <- tset[[a]]; tb <- tset[[b]]
    comp <- 1 - length(intersect(ta, tb)) / length(union(ta, tb))
    expect_equal(res$complementarity[i], comp)
    za <- unname(z[a]); zb <- unname(z[b])
    bp <- za < 0 && zb < 0
    expect_equal(res$both_proximal[i], bp)
    expected_syn <- if (bp) comp * (-max(za, zb)) else NA_real_
    expect_equal(res$synergy_score[i], expected_syn)
  }
  expect_true(all(is.na(res$s_ab)))
  expect_true(all(is.na(res$cheng_class)))
  expect_true(all(res$separation_method == "jaccard"))
})

## ---------------------------------------------------------------------------
## separation = "network" -- end to end (STRINGdb mocked)
## ---------------------------------------------------------------------------

test_that("network mode: hand-checkable s_AB, P2 class, and singleton gating", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  ## A: one target ; B: three targets ; path graph s1-s2-s3-s4
  proj <- .synergy_add_condition(proj, "SYN", list(A = "PA1", B = c("PB1", "PB2", "PB3")))
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(c("A", "B"), condition = "SYN", z = c(-2, -2))
  mock <- .synergy_mock_string(c("PA1", "PB1", "PB2", "PB3"))
  .synergy_use_mock(mock)

  proj <- network_synergy(proj, condition = "SYN", disease = "SOME_DISEASE", separation = "network", pairs = "all")
  res <- patliRResults(proj, "network_synergy")
  row <- res[res$compound_a == "A" & res$compound_b == "B", ]
  expect_equal(nrow(row), 1L)

  ## d_aa(A) = 0 (singleton) ; d_bb(B): mins among s2,s3,s4 = (1,1,1)/3 = 1
  ## d_ab: A={s1} B={s2,s3,s4}; row-min 1 ; col-mins 1,2,3 => (1+6)/(1+3) = 1.75
  ## s_ab = 1.75 - (0 + 1)/2 = 1.25
  expect_equal(row$s_ab, 1.25)
  expect_true(row$separated)
  expect_true(row$singleton_a)
  expect_false(row$singleton_b)
  expect_equal(row$n_targets_a_mapped, 1L)
  expect_equal(row$n_targets_b_mapped, 3L)
  expect_identical(row$cheng_class, "P2")           # separated & both proximal
  expect_true(row$complementary_exposure)
  expect_true(is.na(row$synergy_score))             # singleton => excluded from the scalar
  expect_true(all(res$species == 9606))
  expect_true(all(res$string_version == "12.0"))
})

test_that("network mode: zero mapped targets => s_ab NA, cheng_class NA, never dropped", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  proj <- .synergy_add_condition(proj, "SYN", list(A = c("PA1", "PA2"), B = c("PB1", "PB2")))
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(c("A", "B"), condition = "SYN", z = c(-2, -2))
  ## the map refuses every one of A's accessions
  mock <- .synergy_mock_string(c("PA1", "PA2", "PB1", "PB2"), singleton_targets = c("PA1", "PA2"))
  .synergy_use_mock(mock)

  proj <- network_synergy(proj, condition = "SYN", disease = "SOME_DISEASE", separation = "network", pairs = "all")
  row <- patliRResults(proj, "network_synergy")
  expect_equal(nrow(row), 1L)
  expect_equal(row$n_targets_a_mapped, 0L)
  expect_true(is.na(row$s_ab))
  expect_true(is.na(row$cheng_class))
  expect_true(is.na(row$synergy_score))
})

test_that("network mode: a disconnected interactome LCC => s_ab NA and a logged line", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  proj <- .synergy_add_condition(proj, "SYN", list(A = c("PA1", "PA2"), B = c("PB1", "PB2")))
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(c("A", "B"), condition = "SYN", z = c(-2, -2))
  g <- igraph::graph_from_literal(s1 - s2, s3 - s4)     # two components
  mock <- .synergy_mock_string(c("PA1", "PA2", "PB1", "PB2"), graph = g)
  .synergy_use_mock(mock)

  proj <- network_synergy(proj, condition = "SYN", disease = "SOME_DISEASE", separation = "network", pairs = "all")
  row <- patliRResults(proj, "network_synergy")
  expect_true(is.na(row$s_ab))
  expect_true(any(grepl("disconnected", projectLog(proj)$message)))
})

test_that("network mode: aborts on a species / version / score_threshold mismatch with the proximity rows", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  proj <- .synergy_add_condition(proj, "SYN", list(A = c("PA1", "PA2"), B = c("PB1", "PB2")))
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(
    c("A", "B"), condition = "SYN", z = c(-2, -2), score_threshold = 700
  )
  mock <- .synergy_mock_string(c("PA1", "PA2", "PB1", "PB2"))
  .synergy_use_mock(mock)
  expect_error(
    network_synergy(proj, condition = "SYN", disease = "SOME_DISEASE",
                    separation = "network", pairs = "all", score_threshold = 400),
    "different STRING interactome"
  )
})

test_that("network_synergy(separation = 'network') aborts cleanly when STRINGdb is absent", {
  skip_if(requireNamespace("STRINGdb", quietly = TRUE), "STRINGdb is installed")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_proximity") <- .synergy_fake_prox(.synergy_flo_compounds(proj))
  expect_error(
    network_synergy(proj, condition = "FLO-ET", disease = "SOME_DISEASE", separation = "network"),
    "STRINGdb"
  )
})
