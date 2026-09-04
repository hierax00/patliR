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
  expect_equal(res$compound_id, c("C0001", "C0002", "C0006", "C0007"))
  expect_equal(
    round(res$z_score, 6),
    c(4.88244, 3.943095, 7.60343, 3.691396),
    tolerance = 1e-5
  )
  expect_true(all(res$n_tests_in_family == nrow(res)))
  expect_true(all(res$species == 9606 & res$string_version == "12.0" & res$score_threshold == 400))
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

test_that("network_proximity() end-to-end is not exercised automatically -- needs a real STRING download", {
  testthat::skip_if_not_installed("STRINGdb")
  skip_on_cran()
  testthat::skip(
    "network_proximity() end-to-end needs a real STRINGdb flat-file download (tens-hundreds of MB) -- verify manually against real data, not in the automated suite."
  )
})
