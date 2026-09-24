test_that("network_proximity() requires network_build() to have run first", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_proximity(proj, disease = "EFO_0000000"), "network_edges")
})

test_that("network_proximity() default path requires a disease_genes slot", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  expect_error(network_proximity(proj, condition = "FLO-ET", disease = "EFO_0000000"), "disease_genes")
})

test_that("network_proximity(disease_genes = \"targets_disease\") requires targets_disease_filter() to have run", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  expect_error(
    network_proximity(proj, condition = "FLO-ET", disease = "EFO_0000000",
                      disease_genes = "targets_disease"),
    "targets_disease"
  )
})

test_that("network_proximity() errors on a condition network_build() never built", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  expect_error(
    network_proximity(proj, condition = "NOT-A-REAL-CONDITION", disease = "EFO_0000000"),
    "not built"
  )
})

test_that("network_proximity() errors clearly on an unknown disease_id", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  ## Fabricate a minimal disease_genes table directly (bypassing a live
  ## Open Targets call) so this test is self-contained and offline.
  proj <- disease_genes_import(
    proj,
    data.frame(uniprot_id = unique(patliRResults(proj, "network_edges")$uniprot_id)[1]),
    disease_id = "SOME_REAL_DISEASE", disease_name = "x", source = "synthetic"
  )

  expect_error(
    network_proximity(proj, condition = "FLO-ET", disease = "NOT-A-REAL-DISEASE"),
    "disease_genes"
  )
})

test_that(".network_value_bins() builds >=min_per_bin bins from consecutive degrees (Guney rule)", {
  set.seed(1)
  degree_all <- stats::setNames(
    pmax(1L, rpois(2000, 4) + rbinom(2000, 40, 0.1)),
    paste0("n", 1:2000)
  )
  bins <- patliR:::.network_value_bins(degree_all, min_per_bin = 100)

  sizes <- lengths(bins)
  expect_true(all(utils::head(sizes, -1) >= 100))          # all but last >= 100
  deg_ranges <- lapply(bins, function(ix) range(degree_all[ix]))
  for (i in seq_len(length(bins) - 1)) {                   # contiguous in degree
    expect_lte(deg_ranges[[i]][2], deg_ranges[[i + 1]][1])
  }
})

test_that(".network_resample_matched() returns an equal-size set of distinct, degree-matched nodes", {
  set.seed(2)
  degree_all <- stats::setNames(sample(1:60, 400, replace = TRUE), paste0("n", 1:400))
  bins <- patliR:::.network_value_bins(degree_all, min_per_bin = 50)

  node_names <- names(degree_all)
  bin_of_node <- integer(length(degree_all))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b

  input <- node_names[1:12]
  out <- patliR:::.network_resample_matched(input, node_names, bins, bin_of_node)
  expect_length(out, length(input))
  expect_false(anyDuplicated(out) > 0)
  expect_true(all(out %in% node_names))

  bin_of <- bin_of_node
  names(bin_of) <- node_names
  expect_equal(unname(bin_of[out]), unname(bin_of[input]))
})

test_that("network_proximity schema carries p_adjusted whether or not any row is produced (zero-row rerun stability)", {
  ## Spec 1.8: p_adjusted used to be added only when nrow(result) > 0 and
  ## the empty-row constructor did not declare it, so a zero-row rerun
  ## produced a 12-col frame against a 13-col table and .network_upsert()'s
  ## rbind failed. Both sides must now be 13 cols, same names, same order.
  empty <- patliR:::.empty_network_proximity_row()
  expect_true("p_adjusted" %in% names(empty))

  populated <- data.frame(
    condition = "X", compound_id = "C1", disease_id = "D",
    disease_gene_source = "disease_genes",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, n_overlap = 0L,
    d_observed = 1, d_random_mean = 2, d_random_sd = 1,
    z_score = -1, p_empirical = 0.1, n_random = 10L, seed_used = 1L,
    species = 9606, string_version = "12.0", score_threshold = 400,
    n_tests_in_family = NA_integer_, p_adjusted = NA_real_,   # how network_proximity() builds the row
    stringsAsFactors = FALSE
  )
  expect_identical(names(empty), names(populated))
})

test_that(".network_value_bins() count-per-degree pass matches a brute-force recount", {
  set.seed(7)
  degree_all <- stats::setNames(sample(1:40, 600, replace = TRUE), paste0("n", 1:600))
  bins <- patliR:::.network_value_bins(degree_all, min_per_bin = 100)
  ## every node accounted for exactly once, bins non-empty
  expect_equal(sort(unlist(bins, use.names = FALSE)), seq_along(degree_all))
  expect_true(all(lengths(bins) > 0))
  expect_true(all(utils::head(lengths(bins), -1) >= 100))
})

test_that(".network_closest_distance() tolerates duplicate source/target ids", {
  ## igraph::distances()'s underlying C routine errors outright ("Target
  ## vertex list must not have any duplicates") if `to=` has a repeated
  ## vertex. The observed source/target sets can still carry a repeat, so
  ## .network_closest_distance() dedups defensively. A tiny synthetic
  ## igraph reproduces the constraint, no STRINGdb needed.
  g <- igraph::graph_from_data_frame(
    data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"), stringsAsFactors = FALSE),
    directed = FALSE
  )
  expect_no_error(
    d <- patliR:::.network_closest_distance(g, c("a", "a", "b"), c("d", "d"))
  )
  expect_true(is.numeric(d) && length(d) == 1)
})

test_that("network_proximity(): a disease_genes set disjoint from targets_imported gives n_overlap == 0 (STRINGdb mocked)", {
  ## Spec 2.0 test (a): with an independent disease gene set, S and T do
  ## not overlap by construction, so n_overlap is 0.
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])
  disease_uni <- paste0("DIS", seq_len(6))            # disjoint synthetic accessions

  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_TEST", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_TEST",
                            n_random = 12, seed = 1)
  res <- patliRResults(proj, "network_proximity")
  expect_gt(nrow(res), 0)
  expect_true(all(res$n_overlap == 0))
  expect_true(all(res$disease_gene_source == "disease_genes"))
})

test_that("network_proximity(): d_observed == 0 exactly when S is a subset of T, surfaced in n_overlap (STRINGdb mocked)", {
  ## Spec 2.0 test (c).
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])

  ## disease module = every target in the condition => S subset of T for
  ## every compound.
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = cmp_uni),
    disease_id = "D_ALL", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(cmp_uni)
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_ALL",
                            n_random = 12, seed = 1)
  res <- patliRResults(proj, "network_proximity")
  expect_gt(nrow(res), 0)
  expect_true(all(res$d_observed == 0))
  expect_equal(res$n_overlap, res$n_targets_mapped)
})

test_that("network_proximity(): partial overlap gives n_overlap > 0 AND d_observed > 0 (spec 2.0 (c) converse, STRINGdb mocked)", {
  ## Spec 2.0 test (c), the untested direction: S not a subset of T but
  ## S ∩ T != empty must give n_overlap > 0 AND d_observed > 0 -- the
  ## assertion that would catch a future regression reintroducing a d == 0
  ## special case.
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  ct <- unique(edges[edges$condition == cond, c("compound_id", "uniprot_id")])
  by_c <- lapply(split(ct$uniprot_id, ct$compound_id), unique)
  multi <- names(by_c)[lengths(by_c) >= 2][1]
  skip_if(is.na(multi), "fixture has no compound with >= 2 distinct targets")
  shared <- by_c[[multi]][1]

  ## disease module = 1 of `multi`'s targets + 5 disjoint synthetic ones
  disease_uni <- c(shared, paste0("DIS", seq_len(5)))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_PARTIAL", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(unique(ct$uniprot_id), disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_PARTIAL",
                            n_random = 12, seed = 1)
  res <- patliRResults(proj, "network_proximity")
  row <- res[res$compound_id == multi, ]
  expect_equal(nrow(row), 1L)
  expect_gt(row$n_overlap, 0)                       # S ∩ T != empty
  expect_lt(row$n_overlap, row$n_targets_mapped)    # S not a subset of T
  expect_gt(row$d_observed, 0)                      # ... so d_observed > 0
})

test_that("network_proximity(disease_genes = \"targets_disease\") warns 'circular' and stamps the source column (STRINGdb mocked)", {
  ## Spec 2.0 test (b), plus the gap: the emitted rows must carry
  ## disease_gene_source == "targets_disease".
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  ct <- unique(edges[edges$condition == cond, c("compound_id", "uniprot_id")])
  patliRResults(proj, "targets_disease") <- data.frame(
    compound_id = ct$compound_id, target_id = ct$uniprot_id,
    disease_id = "TD_DIS", association_score = 0.5, stringsAsFactors = FALSE
  )
  fake <- .fake_string_db(unique(ct$uniprot_id))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  expect_warning(
    proj <- network_proximity(proj, condition = cond, disease = "TD_DIS",
                              disease_genes = "targets_disease", n_random = 12, seed = 1),
    "circular"
  )
  res <- patliRResults(proj, "network_proximity")
  expect_gt(nrow(res), 0)
  expect_true(all(res$disease_gene_source == "targets_disease"))
})

test_that("network_proximity() z_score is unchanged after the .network_string_lcc() extraction (behaviour-identical refactor)", {
  ## Piece 12 factored the LCC-restriction + degree-binning block out of
  ## network_proximity() into .network_string_lcc(). The graph, the bins
  ## and therefore every z_score must be byte-identical on the fixture.
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])
  disease_uni <- paste0("DIS", seq_len(6))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_PIN", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 1)
  res <- patliRResults(proj, "network_proximity")
  res <- res[order(res$compound_id), ]
  ## Frozen values updated when prep_binarize()'s Q1 threshold was fixed
  ## to include zero-abundance compounds in the quantile (its documented
  ## "Q1 across all compounds" contract) -- more compounds now correctly
  ## clear the presence threshold, so the FLO-ET network (and therefore
  ## this fixture's z_scores) is different.
  expect_equal(res$compound_id, c("C0001", "C0002", "C0004", "C0006", "C0007", "C0008"))
  expect_equal(
    round(res$z_score, 6),
    c(5.996314, 4.445597, 2.636711, 4.089164, 4.402743, 3.13524),
    tolerance = 1e-5
  )
  expect_true(all(res$n_tests_in_family == nrow(res)))
  expect_true(all(res$species == 9606 & res$string_version == "12.0" & res$score_threshold == 400))
})

test_that("network_proximity(store_null = FALSE) reproduces byte-identical main-table output to store_null = TRUE (regression)", {
  ## Phase 4 3.7: store_null must not change how d_random_mean/d_random_sd/
  ## z_score/p_empirical/p_adjusted are computed -- it only captures an
  ## intermediate value the loop already had. Same fixture/seed as the
  ## pinned LCC-extraction test above; the two calls must give identical
  ## main-table rows.
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])
  disease_uni <- paste0("DIS", seq_len(6))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_PIN", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj_false <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 1, store_null = FALSE)
  res_false <- patliRResults(proj_false, "network_proximity")
  res_false <- res_false[order(res_false$compound_id), ]

  proj_true <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 1, store_null = TRUE)
  res_true <- patliRResults(proj_true, "network_proximity")
  res_true <- res_true[order(res_true$compound_id), ]

  expect_identical(res_false, res_true)
  expect_null(patliRResults(proj_false, "network_proximity_null"))
})

test_that("network_proximity(store_null = TRUE) writes network_proximity_null with the right shape and schema", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])
  disease_uni <- paste0("DIS", seq_len(6))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_PIN", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 1, store_null = TRUE)
  main <- patliRResults(proj, "network_proximity")
  null <- patliRResults(proj, "network_proximity_null")

  expect_identical(
    names(null),
    c("condition", "compound_id", "disease_id", "disease_gene_source", "draw", "d_random")
  )
  ## one row per (compound, disease, draw) -- n_random rows per compound
  ## that received a main-table row (spec 3.7).
  expect_equal(nrow(null), 12L * nrow(main))
  expect_true(setequal(unique(null$compound_id), main$compound_id))
  for (cp in main$compound_id) {
    expect_equal(sort(null$draw[null$compound_id == cp]), 1:12)
  }
  expect_true(is.integer(null$draw))
  expect_true(is.double(null$d_random))
})

test_that("network_proximity(store_null = FALSE) does not wipe a previously stored network_proximity_null (opt-in, not opt-out)", {
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  cmp_uni <- unique(edges$uniprot_id[edges$condition == cond])
  disease_uni <- paste0("DIS", seq_len(6))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_PIN", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  proj <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 1, store_null = TRUE)
  null_before <- patliRResults(proj, "network_proximity_null")
  expect_gt(nrow(null_before), 0)

  ## A second call with store_null = FALSE must leave the stored draws
  ## untouched -- not wipe them, not re-key them, not touch them at all.
  proj <- network_proximity(proj, condition = cond, disease = "D_PIN", n_random = 12, seed = 2, store_null = FALSE)
  null_after <- patliRResults(proj, "network_proximity_null")
  expect_identical(null_before, null_after)
})

test_that(".network_upsert() back-fills the new proximity provenance columns on a legacy table (no abort)", {
  proj <- .network_stats_test_setup()
  legacy <- data.frame(
    condition = "X", compound_id = "C1", disease_id = "D",
    disease_gene_source = "disease_genes",
    n_targets_mapped = 1L, n_disease_genes_mapped = 1L, n_overlap = 0L,
    d_observed = 1, d_random_mean = 2, d_random_sd = 1,
    z_score = -1, p_empirical = 0.1, n_random = 10L, seed_used = 1L, p_adjusted = 0.1,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- legacy

  wide <- patliR:::.empty_network_proximity_row()
  wide[1, ] <- NA
  wide$condition <- "Y"; wide$compound_id <- "C2"; wide$disease_id <- "D2"
  wide$z_score <- -1.5

  expect_no_error(
    merged <- patliR:::.network_upsert(
      proj, "network_proximity", wide, c("condition", "disease_id", "compound_id")
    )
  )
  expect_true(all(c("species", "string_version", "score_threshold", "n_tests_in_family") %in% names(merged)))
  expect_true(all(c("X", "Y") %in% merged$condition))
})

test_that("network_proximity() drops stale scores after a compound loses all edges", {
  skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  ids <- unique(edges$uniprot_id)
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = ids),
    disease_id = "D_RERUN", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(ids)
  local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")
  proj <- network_proximity(proj, condition = c("FLO-ET", "LEA-ET"),
                            disease = "D_RERUN", n_random = 12, seed = 1)
  before <- patliRResults(proj, "network_proximity")
  removed <- before$compound_id[before$condition == "FLO-ET"][1]
  other_disease <- before[before$condition == "FLO-ET", , drop = FALSE]
  other_disease$disease_id <- "D_OTHER"
  patliRResults(proj, "network_proximity") <- rbind(before, other_disease)
  patliRResults(proj, "network_edges") <- edges[
    !(edges$condition == "FLO-ET" & edges$compound_id == removed), , drop = FALSE]

  proj <- network_proximity(proj, condition = "FLO-ET", disease = "D_RERUN",
                            n_random = 12, seed = 1)
  after <- patliRResults(proj, "network_proximity")
  current <- after[after$condition == "FLO-ET" & after$disease_id == "D_RERUN", ]
  expect_false(removed %in% current$compound_id)
  expect_equal(nrow(current), sum(before$condition == "FLO-ET") - 1L)
  expect_equal(after$z_score[after$condition == "LEA-ET"],
               before$z_score[before$condition == "LEA-ET"])
  expect_equal(after$compound_id[after$disease_id == "D_OTHER"], other_disease$compound_id)
})

## ---- Reference (pre-batching) null computation, kept as the oracle -------
## Verbatim copies of the per-draw implementation network_proximity() used
## before the null distances were batched into one multi-source BFS per
## draw (.network_closest_distance_null()): one character-based
## degree-matched draw and one igraph::distances() call per (compound,
## draw). Slow, but independent of the code under test.
.ref_resample_matched <- function(string_ids, node_names, bins, bin_of_node) {
  idx <- match(string_ids, node_names)
  by_bin <- split(seq_along(idx), bin_of_node[idx])
  out <- character(length(string_ids))
  for (bs in names(by_bin)) {
    slots <- by_bin[[bs]]
    pool <- bins[[as.integer(bs)]]
    k <- length(slots)
    picks <- if (k <= length(pool)) {
      pool[sample.int(length(pool), k)]
    } else {
      pool[sample.int(length(pool), k, replace = TRUE)]
    }
    out[slots] <- node_names[picks]
  }
  out
}

.ref_closest_distance <- function(g, source_ids, target_ids) {
  d <- igraph::distances(g, v = unique(source_ids), to = unique(target_ids), weights = NA)
  d[!is.finite(d)] <- NA_real_
  row_min <- suppressWarnings(apply(d, 1, min, na.rm = TRUE))
  mean(row_min[is.finite(row_min)])
}

## The null exactly as the old loop drew and scored it: all T' draws
## first, then each compound's S' draws in compound order.
.ref_proximity_null <- function(g, source_sets, target_string, n_random, seed, bins, bin_of_node) {
  node_names <- igraph::V(g)$name
  restore <- patliR:::.with_seed(seed)
  on.exit(restore(), add = TRUE)
  t_rand <- lapply(seq_len(n_random), function(j)
    .ref_resample_matched(target_string, node_names, bins, bin_of_node))
  lapply(source_sets, function(s) vapply(seq_len(n_random), function(j) {
    .ref_closest_distance(g, .ref_resample_matched(s, node_names, bins, bin_of_node), t_rand[[j]])
  }, numeric(1)))
}

test_that(".network_closest_distance_null() equals the per-draw .network_closest_distance() (chunked, Inf, duplicates)", {
  ## Synthetic scale-free graph plus a detached 3-node component, so some
  ## draws have sources with no finite path to T' (dropped, as before) and
  ## some have none at all (NaN, as before); a duplicated source node in a
  ## draw exercises the unique(). A tiny max_cells forces many chunks.
  set.seed(11)
  g <- igraph::sample_pa(300, m = 2, directed = FALSE)
  g <- igraph::add_vertices(g, 3)
  g <- igraph::add_edges(g, c(301, 302, 302, 303))
  igraph::V(g)$name <- paste0("n", seq_len(igraph::vcount(g)))
  nm <- igraph::V(g)$name
  n_random <- 23
  t_rand <- lapply(seq_len(n_random), function(j) sample.int(300, 15))
  s_rand <- lapply(1:4, function(k) lapply(seq_len(n_random), function(j) sample.int(303, 6)))
  s_rand[[1]][[1]] <- c(301L, 302L)              # entirely off T''s component
  s_rand[[2]][[3]] <- c(5L, 5L, 303L, 7L)        # duplicate + one unreachable
  t_rand[[4]] <- c(t_rand[[4]], t_rand[[4]][1])  # duplicate target

  got <- patliR:::.network_closest_distance_null(g, s_rand, t_rand, max_cells = 50)
  want <- lapply(s_rand, function(draws) vapply(seq_len(n_random), function(j)
    .ref_closest_distance(g, nm[draws[[j]]], nm[t_rand[[j]]]), numeric(1)))

  expect_identical(got, want)
  expect_true(is.nan(got[[1]][1]))
  expect_true(is.finite(got[[2]][3]))
  ## the chunk size does not change anything
  expect_identical(patliR:::.network_closest_distance_null(g, s_rand, t_rand), want)
})

test_that(".network_resample_matched_idx() consumes the RNG exactly like .network_resample_matched()", {
  set.seed(3)
  g <- igraph::sample_pa(400, m = 2, directed = FALSE)
  igraph::V(g)$name <- paste0("n", seq_len(400))
  deg <- igraph::degree(g)
  bins <- patliR:::.network_value_bins(deg, min_per_bin = 30)
  bon <- integer(length(deg))
  for (b in seq_along(bins)) bon[bins[[b]]] <- b
  ids <- names(deg)[c(1:5, 50, 399)]
  set.seed(99); a <- replicate(20, patliR:::.network_resample_matched(ids, names(deg), bins, bon), simplify = FALSE)
  set.seed(99); b <- replicate(20, patliR:::.network_resample_matched_idx(match(ids, names(deg)), bins, bon), simplify = FALSE)
  set.seed(99); r <- replicate(20, .ref_resample_matched(ids, names(deg), bins, bon), simplify = FALSE)
  expect_identical(a, r)
  expect_identical(lapply(b, function(i) names(deg)[i]), r)
})

test_that("network_proximity() batched null reproduces the per-draw reference implementation exactly (STRINGdb mocked)", {
  ## Regression for the null-distance batching: same seed => same random
  ## draws => identical d_random per draw, and therefore identical
  ## d_random_mean / d_random_sd / z_score / p_empirical, versus the old
  ## one-distances()-per-draw loop reimplemented above as the oracle.
  testthat::skip_if_not_installed("STRINGdb")
  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  ct <- unique(edges[edges$condition == cond, c("compound_id", "uniprot_id")])
  cmp_uni <- unique(ct$uniprot_id)
  disease_uni <- c(cmp_uni[1:2], paste0("DIS", seq_len(8)))
  proj <- disease_genes_import(
    proj, data.frame(uniprot_id = disease_uni),
    disease_id = "D_REF", disease_name = "synthetic", source = "synthetic"
  )
  fake <- .fake_string_db(c(cmp_uni, disease_uni))
  testthat::local_mocked_bindings(.network_stringdb = function(...) fake, .package = "patliR")

  n_random <- 40
  proj <- network_proximity(proj, condition = cond, disease = "D_REF",
                            n_random = n_random, seed = 20240, store_null = TRUE)
  res <- patliRResults(proj, "network_proximity")
  nul <- patliRResults(proj, "network_proximity_null")

  ## oracle, on the same LCC / bins / mapping network_proximity() used
  lcc <- patliR:::.network_string_lcc(proj, 9606, "12.0", 400)
  g <- lcc$graph
  to_string <- function(u) {
    s <- unique(stats::na.omit(fake$map(data.frame(uniprot_id = u), "uniprot_id")$STRING_id))
    s[s %in% igraph::V(g)$name]
  }
  compounds <- unique(ct$compound_id)
  source_sets <- lapply(split(ct$uniprot_id, ct$compound_id)[compounds], to_string)
  source_sets <- source_sets[lengths(source_sets) > 0]
  target_string <- to_string(disease_uni)
  ref <- .ref_proximity_null(g, source_sets, target_string, n_random, 20240,
                             lcc$bins, lcc$bin_of_node)

  expect_setequal(res$compound_id, names(ref))
  for (cp in names(ref)) {
    got <- nul[nul$compound_id == cp, ]
    got <- got$d_random[order(got$draw)]
    want_raw <- ref[[cp]]
    expect_identical(got, ifelse(is.finite(want_raw), want_raw, NA_real_))

    want <- want_raw[is.finite(want_raw)]
    row <- res[res$compound_id == cp, ]
    d_obs <- .ref_closest_distance(g, source_sets[[cp]], target_string)
    expect_identical(row$d_observed, d_obs)
    want_sd <- if (length(want) > 1) stats::sd(want) else NA_real_
    expect_identical(row$d_random_mean, mean(want))
    expect_identical(row$d_random_sd, want_sd)
    expect_identical(
      row$z_score,
      if (!is.na(want_sd) && want_sd > 0) (d_obs - mean(want)) / want_sd else NA_real_
    )
    expect_identical(row$p_empirical, (1 + sum(want <= d_obs)) / (length(want) + 1))
    expect_identical(row$n_random, length(want))
  }
})

test_that("network_proximity() end-to-end is not exercised automatically -- needs a real STRING download", {
  testthat::skip_if_not_installed("STRINGdb")
  skip_on_cran()
  testthat::skip(
    "network_proximity() end-to-end needs a real STRINGdb flat-file download (tens-hundreds of MB) -- verify manually against real data, not in the automated suite."
  )
})
