## ---------------------------------------------------------------------------
## Piece 14: network_degeneracy() with GO semantic similarity
## ---------------------------------------------------------------------------

## Probe for the local GOSemSim/Rcpp IC-path breakage (pre-review 14-E1):
## infoContentMethod_cpp fails to resolve its Rcpp entry point in this
## install, so Resnik/Lin/Rel/Jiang cannot be exercised here.
.ic_path_broken <- function() {
  if (!requireNamespace("GOSemSim", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(TRUE)
  }
  res <- tryCatch({
    sd <- GOSemSim::godata(annoDb = "org.Hs.eg.db", keytype = "ENTREZID",
                           ont = "BP", computeIC = TRUE)
    GOSemSim::clusterSim(c("207", "208"), c("2064", "1956"),
                         semData = sd, measure = "Resnik", drop = "IEA", combine = "BMA")
    FALSE
  }, error = function(e) TRUE)
  isTRUE(res)
}

## The shared test fixture's project universe is only ~10 genes, so the
## `annotation = "direct"` null legitimately warns every run (one
## annotation-count bin; n_pairs / (n_random + 1) above the BH-resolution
## threshold). Those two warnings are asserted once, in their own test
## below, and muffled everywhere else so they do not drown the suite.
.ndeg <- function(...) {
  withCallingHandlers(
    network_degeneracy(...),
    warning = function(w) {
      if (grepl("annotation-count bin|BH-adjusted p cannot resolve", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("network_degeneracy() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_degeneracy(proj), "network_edges")
})

test_that("network_degeneracy() errors on a condition network_build() never built", {
  proj <- .network_stats_test_setup()
  expect_error(network_degeneracy(proj, condition = "NOT-A-REAL-CONDITION"), "not built")
})

test_that("annotation = 'enriched'/'jaccard' still require network_enrich(); 'direct' does not (14-R3)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup() ## network_build() only, no enrich
  expect_error(network_degeneracy(proj, condition = "FLO-ET", annotation = "enriched"),
               "network_enrich")
  expect_error(network_degeneracy(proj, condition = "FLO-ET", annotation = "jaccard"),
               "network_enrich")
  ## default (direct) must NOT abort on the missing enrichment slot
  expect_no_error(
    .ndeg(proj, condition = "FLO-ET", n_random = 8, seed = 1)
  )
})

test_that("network_degeneracy(annotation = 'direct') warns on a degenerate pool and an unresolvable BH floor (14-U1, 14-SE8)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  msgs <- character(0)
  withCallingHandlers(
    network_degeneracy(proj, condition = "FLO-ET", n_random = 8, seed = 1),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  ## fixture project universe is ~10 genes -> one annotation-count bin (14-U1)
  expect_true(any(grepl("annotation-count bin", msgs)))
  ## 6 pairs vs n_random = 8 -> n_pairs / (n_random + 1) > 0.05 (14-SE8)
  expect_true(any(grepl("BH-adjusted p cannot resolve", msgs)))
})

test_that("network_degeneracy() n_random is rejected below 1; n_random = 1 pins sd/z = NA, p in {0.5, 1} (14-U8)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  expect_error(network_degeneracy(proj, condition = "FLO-ET", n_random = 0), "n_random")

  proj <- .ndeg(proj, condition = "FLO-ET", n_random = 1, seed = 42)
  res <- patliRResults(proj, "network_degeneracy")
  res <- res[!is.na(res$functional_similarity), ]
  skip_if(nrow(res) == 0, "no scorable pairs in the fixture")
  expect_true(all(is.na(res$sim_random_sd)))
  expect_true(all(is.na(res$z_score)))
  expect_true(all(res$p_empirical %in% c(0.5, 1)))
  expect_true(any(res$p_empirical %in% c(0.5, 1)))
})

test_that("network_degeneracy(annotation = 'direct') runs, functional_similarity in [0, 1] under Wang, schema is stable", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  proj <- .ndeg(proj, condition = "FLO-ET", n_random = 12, seed = 1)
  res <- patliRResults(proj, "network_degeneracy")

  expect_true(all(c(
    "condition", "compound_a", "compound_b", "n_targets_a", "n_targets_b",
    "n_genes_a_mapped", "n_genes_b_mapped", "n_pathways_a", "n_pathways_b",
    "target_jaccard", "functional_similarity", "sim_random_mean", "sim_random_sd",
    "z_score", "p_empirical", "p_adjusted", "degeneracy_score",
    "annotation", "ont", "measure", "combine", "drop", "universe",
    "n_random", "seed_used"
  ) %in% names(res)))
  expect_identical(names(res), names(patliR:::.empty_network_degeneracy_row()))
  expect_false("pathway_jaccard" %in% names(res))

  scor <- res[!is.na(res$functional_similarity), ]
  skip_if(nrow(scor) == 0, "no scorable pairs in the fixture")
  ## Wang/BMA is bounded [0, 1] (unlike Resnik) -- see 14-E1 / spec test 2.
  expect_true(all(scor$functional_similarity >= 0 & scor$functional_similarity <= 1))
  expect_true(all(scor$annotation == "direct"))
  expect_true(all(scor$ont == "BP" & scor$measure == "Wang" & scor$combine == "BMA"))
  expect_true(all(scor$drop == "IEA" & scor$universe == "project"))
  ## n_pathways_* are NA in direct mode but the columns are present (14-U5)
  expect_true(all(is.na(res$n_pathways_a)))
  ## degeneracy_score identity
  ok <- !is.na(scor$degeneracy_score)
  expect_equal(scor$degeneracy_score[ok],
               scor$functional_similarity[ok] * (1 - scor$target_jaccard[ok]))
})

test_that("network_degeneracy(annotation = 'direct'): identical target sets give functional_similarity == 1, degeneracy_score == 0 (spec test 3)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  flo <- edges[edges$condition == "FLO-ET", ]
  donor <- flo[flo$compound_id == "C0001", ]
  clone <- donor; clone$compound_id <- "CDUP"
  patliRResults(proj, "network_edges") <- rbind(edges, clone)

  proj <- .ndeg(proj, condition = "FLO-ET", n_random = 6, seed = 1)
  res <- patliRResults(proj, "network_degeneracy")
  row <- res[(res$compound_a == "C0001" & res$compound_b == "CDUP") |
             (res$compound_a == "CDUP" & res$compound_b == "C0001"), ]
  expect_equal(nrow(row), 1L)
  ## quantised to 0.001 by combineScores; it is exactly 1 so this passes
  ## regardless of tolerance -- the spec's "to 1e-8" was meaningless (14-SE6).
  expect_equal(row$functional_similarity, 1)
  expect_equal(row$target_jaccard, 1)
  expect_equal(row$degeneracy_score, 0)
})

test_that("network_degeneracy() fast BMA scorer reproduces GOSemSim::clusterSim() exactly on real gene-set pairs (14-SE5/SE6/SE7)", {
  testthat::skip_if_not_installed("GOSemSim")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  sd <- GOSemSim::godata(annoDb = "org.Hs.eg.db", keytype = "ENTREZID",
                         ont = "BP", computeIC = FALSE)
  ga <- sd@geneAnno
  ga <- ga[!ga$EVIDENCE %in% "IEA", ]
  g2go <- lapply(split(as.character(ga$GO), as.character(ga$ENTREZID)),
                 function(x) unique(x[!is.na(x)]))

  fast <- function(a, b) {
    tva <- unlist(g2go[a], use.names = FALSE)
    tvb <- unlist(g2go[b], use.names = FALSE)
    terms <- sort(unique(c(tva, tvb)))
    M <- GOSemSim::termSim(terms, terms, semData = sd, method = "Wang")
    patliR:::.network_go_bma(tva, tvb, M, "BMA")
  }

  pairs <- list(
    list(a = c("2099", "2064", "1956", "207"), b = c("2064", "1956", "5290")),
    list(a = c("367", "2099"),                 b = c("2064", "1956", "207")),
    list(a = c("2064", "1956", "207", "5290"), b = c("2064", "1956", "207", "5290"))
  )
  for (p in pairs) {
    ref <- GOSemSim::clusterSim(p$a, p$b, semData = sd, measure = "Wang",
                                drop = "IEA", combine = "BMA")
    expect_equal(fast(p$a, p$b), ref)
  }
})

test_that("network_degeneracy(annotation = 'direct'): disjoint unrelated gene sets score below a same-process control (spec test 4)", {
  testthat::skip_if_not_installed("GOSemSim")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  sd <- GOSemSim::godata(annoDb = "org.Hs.eg.db", keytype = "ENTREZID",
                         ont = "BP", computeIC = FALSE)
  ga <- sd@geneAnno
  ga <- ga[!ga$EVIDENCE %in% "IEA", ]
  g2go <- lapply(split(as.character(ga$GO), as.character(ga$ENTREZID)),
                 function(x) unique(x[!is.na(x)]))
  score <- function(a, b) {
    a <- a[a %in% names(g2go)]; b <- b[b %in% names(g2go)]
    tva <- unlist(g2go[a], use.names = FALSE); tvb <- unlist(g2go[b], use.names = FALSE)
    terms <- sort(unique(c(tva, tvb)))
    M <- GOSemSim::termSim(terms, terms, semData = sd, method = "Wang")
    patliR:::.network_go_bma(tva, tvb, M, "BMA")
  }
  ribosomal <- c("6122", "6124", "6125", "6128", "6129", "6132")   # RPL*
  olfactory <- c("26648", "26649", "26658", "26664", "8347")       # OR*
  kinases_a <- c("207", "208", "10000")                            # AKT1/2/3
  kinases_b <- c("5290", "5291", "5293")                           # PIK3C*
  unrelated <- score(ribosomal, olfactory)
  same_proc <- score(kinases_a, kinases_b)
  expect_lt(unrelated, same_proc)
})

test_that("network_degeneracy(annotation = 'direct') is seed-reproducible (spec tests 7, 8)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  seed_before <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  proj1 <- .ndeg(.network_stats_test_setup(), condition = "FLO-ET", n_random = 16, seed = 7)
  seed_after <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  expect_identical(seed_before, seed_after)             # RNG isolation

  proj2 <- .ndeg(.network_stats_test_setup(), condition = "FLO-ET", n_random = 16, seed = 7)
  r1 <- patliRResults(proj1, "network_degeneracy")
  r2 <- patliRResults(proj2, "network_degeneracy")
  r1 <- r1[order(r1$compound_a, r1$compound_b), ]
  r2 <- r2[order(r2$compound_a, r2$compound_b), ]
  expect_equal(r1$z_score, r2$z_score)
  expect_equal(r1$functional_similarity, r2$functional_similarity)
})

test_that("network_degeneracy(annotation = 'jaccard') reproduces the enriched-pathway Jaccard (regression, spec test 1)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  ## universe = "genome" (explicit opt-out of the new default, 20-U1): this
  ## test is about network_degeneracy()'s jaccard math reproducing the
  ## pathway-overlap computation, not about network_enrich()'s background
  ## choice -- universe = "project" restricts the background enough that
  ## this tiny fixture can come back with zero significant GO terms, which
  ## would make this a test of the enrichment universe instead.
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE, universe = "genome")
  proj <- network_degeneracy(proj, condition = "FLO-ET", annotation = "jaccard")
  res <- patliRResults(proj, "network_degeneracy")
  skip_if(nrow(res) == 0, "no enriched pathways for this fixture")

  ## independently recompute pathway Jaccard the way the pathway path does
  ct <- unique(patliRResults(proj, "network_edges"))
  ct <- ct[ct$condition == "FLO-ET", c("compound_id", "uniprot_id")]
  tpe <- patliR:::.network_target_pathway_edges(proj, "FLO-ET", unique(ct$uniprot_id))
  cp <- unique(merge(ct, tpe, by = "uniprot_id")[, c("compound_id", "pathway_id")])
  psets <- split(cp$pathway_id, cp$compound_id)
  for (i in seq_len(nrow(res))) {
    pa <- psets[[res$compound_a[i]]]; pb <- psets[[res$compound_b[i]]]
    if (is.null(pa)) pa <- character(0); if (is.null(pb)) pb <- character(0)
    expect_equal(res$functional_similarity[i],
                 length(intersect(pa, pb)) / length(union(pa, pb)))
  }
  expect_true(all(res$annotation == "jaccard"))
  expect_true(all(is.na(res$z_score) & is.na(res$p_adjusted)))
  ok <- !is.na(res$degeneracy_score)
  expect_equal(res$degeneracy_score[ok],
               res$functional_similarity[ok] * (1 - res$target_jaccard[ok]))
})

test_that("network_degeneracy() IC measure path (Resnik) -- skipped when the local Rcpp entry point is broken (14-E1)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()
  skip_if(.ic_path_broken(), "GOSemSim infoContentMethod_cpp Rcpp entry point broken in this install (14-E1)")

  proj <- .network_stats_test_setup()
  proj <- .ndeg(proj, condition = "FLO-ET", measure = "Resnik",
                n_random = 6, seed = 1)
  res <- patliRResults(proj, "network_degeneracy")
  expect_true(all(res$measure[!is.na(res$measure)] == "Resnik"))
  ## Resnik is UNBOUNDED -- do not assert into [0, 1] (spec test 2 note).
  expect_true(all(res$functional_similarity[!is.na(res$functional_similarity)] >= 0))
})

test_that(".network_value_bins() on GO annotation counts satisfies the contiguity + >= min_per_bin contract (spec test 5)", {
  set.seed(11)
  ## a heavy-tailed annotation-count-like vector
  counts <- stats::setNames(pmax(1L, rpois(1500, 6) + rnbinom(1500, size = 1, mu = 4)),
                            paste0("g", 1:1500))
  bins <- patliR:::.network_value_bins(counts, min_per_bin = 100)
  expect_true(all(utils::head(lengths(bins), -1) >= 100))
  rng <- lapply(bins, function(ix) range(counts[ix]))
  for (i in seq_len(length(bins) - 1)) expect_lte(rng[[i]][2], rng[[i + 1]][1])
  expect_equal(sort(unlist(bins, use.names = FALSE)), seq_along(counts))
})

test_that(".network_resample_matched() rejects an id outside the pool instead of emitting '' (14-U2)", {
  set.seed(3)
  vals <- stats::setNames(sample(1:40, 300, replace = TRUE), paste0("n", 1:300))
  bins <- patliR:::.network_value_bins(vals, min_per_bin = 60)
  bon <- integer(length(vals))
  for (b in seq_along(bins)) bon[bins[[b]]] <- b

  input <- names(vals)[1:10]
  out <- patliR:::.network_resample_matched(input, names(vals), bins, bon)
  expect_length(out, 10)
  expect_false(anyDuplicated(out) > 0)
  bn <- stats::setNames(bon, names(vals))
  expect_equal(unname(bn[out]), unname(bn[input]))          # bin multiset preserved

  expect_error(
    patliR:::.network_resample_matched(c(input, "NOT-IN-POOL"), names(vals), bins, bon),
    "string_ids"
  )
})

test_that("network_degeneracy(): a legacy 10-column network_degeneracy.csv survives a rerun (pathway_jaccard -> functional_similarity, 14-U3)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  legacy <- data.frame(
    condition = "LEGACY-COND", compound_a = "X1", compound_b = "X2",
    n_targets_a = 3L, n_targets_b = 4L, n_pathways_a = 2L, n_pathways_b = 5L,
    target_jaccard = 0.1, pathway_jaccard = 0.7, degeneracy_score = 0.63,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_degeneracy") <- legacy

  expect_no_error(
    proj <- .ndeg(proj, condition = "FLO-ET", n_random = 6, seed = 1)
  )
  res <- patliRResults(proj, "network_degeneracy")
  expect_false("pathway_jaccard" %in% names(res))
  expect_true("functional_similarity" %in% names(res))
  keep <- res[res$condition == "LEGACY-COND", ]
  expect_equal(nrow(keep), 1L)
  expect_equal(keep$functional_similarity, 0.7)   # value carried through under the new name
  expect_true(is.na(keep$annotation))             # per-call constant back-filled NA
})

## ---------------------------------------------------------------------------
## piece-14 verification: universe = "genome" pool capping (blocking fix)
## ---------------------------------------------------------------------------

test_that(".network_cap_pool_keep_observed() never evicts an observed gene", {
  pool <- paste0("G", seq_len(2000))
  observed <- c("G5", "G1999", "G2000", "G1")   # scattered, unlikely to survive a naive sample()
  set.seed(1)
  res <- patliR:::.network_cap_pool_keep_observed(pool, observed, cap = 1500L)
  expect_true(res$capped)
  expect_length(res$pool, 1500L)
  expect_true(all(observed %in% res$pool))
})

test_that(".network_cap_pool_keep_observed() keeps the whole observed set even if it alone exceeds cap", {
  pool <- paste0("G", seq_len(10))
  res <- patliR:::.network_cap_pool_keep_observed(pool, observed = pool, cap = 3L)
  expect_true(res$capped)
  expect_setequal(res$pool, pool)   # can't shrink below the observed set
})

test_that(".network_cap_pool_keep_observed() is a no-op under the cap", {
  pool <- paste0("G", seq_len(5))
  res <- patliR:::.network_cap_pool_keep_observed(pool, observed = "G1", cap = 1500L)
  expect_false(res$capped)
  expect_identical(res$pool, pool)
})

test_that(".network_cap_pool_keep_observed() tolerates observed genes absent from the pool", {
  pool <- paste0("G", seq_len(2000))
  res <- patliR:::.network_cap_pool_keep_observed(pool, observed = c("G1", "NOT_IN_POOL"), cap = 1500L)
  expect_true(res$capped)
  expect_true("G1" %in% res$pool)
  expect_false("NOT_IN_POOL" %in% res$pool)
  expect_length(res$pool, 1500L)
})

test_that("network_degeneracy(annotation = 'enriched') warns instead of silently scoring all-NA when the enrichment has no GO terms", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_not_installed("GOSemSim")
  skip_on_cran()

  proj <- .network_stats_test_setup()
  ## fabricate an enrichment result with only non-GO (e.g. KEGG-style) IDs
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = "FLO-ET", db = "kegg", ID = "hsa00010", Description = "Glycolysis",
    pvalue = 0.01, p.adjust = 0.02, geneID = "P12345/P23456",
    stringsAsFactors = FALSE
  )
  ## give the pathway-edge lookup something to merge against
  local_mocked_bindings(
    .network_target_pathway_edges = function(proj, cond, uniprot_ids, pathway_db = NULL) {
      data.frame(uniprot_id = uniprot_ids, pathway_id = "hsa00010", stringsAsFactors = FALSE)
    }
  )
  expect_warning(
    proj2 <- network_degeneracy(proj, condition = "FLO-ET", annotation = "enriched"),
    "no GO-prefixed"
  )
  res <- patliRResults(proj2, "network_degeneracy")
  if (nrow(res) > 0) expect_true(all(is.na(res$functional_similarity)))
})

## Frozen pre-optimization cores: keep bodies verbatim as regression oracles.

.ndeg_old_degeneracy_direct <- function(proj, cond, ct, all_compounds, sem_data,
                                        g2go, ann_counts, measure, combine, universe,
                                        n_random, project_pool_entrez, edges_all) {
  ## compound -> distinct Entrez (many-to-many UniProt -> Entrez)
  uni_all <- unique(ct$uniprot_id)
  u2e <- .network_uniprot_to_entrez_map(uni_all)
  compound_entrez <- lapply(split(ct$uniprot_id, ct$compound_id), function(u) {
    unique(stats::na.omit(unlist(u2e[unique(u)], use.names = FALSE)))
  })

  ## resampling pool
  pool <- switch(universe,
    project   = project_pool_entrez,
    condition = unique(stats::na.omit(unlist(u2e, use.names = FALSE))),
    genome    = names(g2go)
  )
  pool <- unique(as.character(pool))
  pool <- pool[pool %in% names(g2go)]        # annotated genes only (14-U2)

  ## Cap the pool for the term x term matrix's sake, but never at the cost
  ## of evicting a compound's own (annotated) targets: .network_resample_matched()
  ## looks up each observed gene's own annotation-count bin to draw its
  ## replacement, and the pairwise loop below intersects each compound's
  ## target set with `pool` before scoring, so a gene missing from `pool`
  ## silently drops out of every pair it's in.
  cp_res <- .network_cap_pool_keep_observed(
    pool, unlist(compound_entrez, use.names = FALSE), cap = 1500L
  )
  pool <- cp_res$pool
  pool_capped <- cp_res$capped

  if (length(pool) < 2) {
    proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
      message = paste0("condition '", cond, "': resampling pool (universe = '", universe,
                       "') has < 2 annotated genes; no rows"))
    return(list(rows = .empty_network_degeneracy_row(), proj = proj))
  }

  ## annotation-count bins over the pool (shared helper with proximity)
  pool_counts <- ann_counts[pool]
  bins <- .network_value_bins(pool_counts, min_per_bin = 100)
  bin_of_node <- integer(length(pool))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b
  if (length(bins) < 3) {
    cli::cli_warn(c(
      "!" = "{.fn network_degeneracy}: condition {.val {cond}} resampling pool (universe = {.val {universe}}, {length(pool)} gene{?s}) yields only {length(bins)} annotation-count bin{?s}.",
      "i" = "The permutation null is close to uniform sampling; a larger {.arg universe} gives a sharper null."
    ))
  }
  if (pool_capped) {
    cli::cli_warn(c(
      "!" = "{.fn network_degeneracy}: condition {.val {cond}} universe = {.val {universe}} pool capped at {length(pool)} gene{?s} (every observed target kept; the rest randomly subsampled) for the term-similarity matrix."
    ))
  }

  ## term x term similarity matrix over the pool's term union, ONCE
  terms_pool <- sort(unique(unlist(g2go[pool], use.names = FALSE)))
  sim_mat <- GOSemSim::termSim(terms_pool, terms_pool, semData = sem_data, method = measure)
  if (is.null(dim(sim_mat))) {
    sim_mat <- matrix(sim_mat, nrow = length(terms_pool), ncol = length(terms_pool),
                      dimnames = list(terms_pool, terms_pool))
  }

  ## genes that lost all annotation after the drop filter (log once)
  lost <- setdiff(unique(unlist(compound_entrez, use.names = FALSE)), names(g2go))
  if (length(lost) > 0) {
    proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
      message = paste0("condition '", cond, "': ", length(lost),
                       " Entrez gene(s) had no GO annotation after the drop filter and were excluded from scoring"))
  }

  pairs <- utils::combn(all_compounds, 2, simplify = FALSE)
  pair_rows <- vector("list", length(pairs))
  n_skipped <- 0L

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1]; b <- pairs[[i]][2]
    ea_all <- compound_entrez[[a]]; eb_all <- compound_entrez[[b]]
    if (is.null(ea_all)) ea_all <- character(0)
    if (is.null(eb_all)) eb_all <- character(0)
    ## scored sets: annotated AND in the resampling pool (14-U2)
    ea <- intersect(ea_all, pool)
    eb <- intersect(eb_all, pool)

    if (length(ea) == 0 || length(eb) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }

    tv_a <- unlist(g2go[ea], use.names = FALSE)
    tv_b <- unlist(g2go[eb], use.names = FALSE)
    sim_obs <- .network_go_bma(tv_a, tv_b, sim_mat, combine)

    sim_rand <- vapply(seq_len(n_random), function(j) {
      ra <- .network_resample_matched(ea, pool, bins, bin_of_node)
      rb <- .network_resample_matched(eb, pool, bins, bin_of_node)
      .network_go_bma(unlist(g2go[ra], use.names = FALSE),
                      unlist(g2go[rb], use.names = FALSE), sim_mat, combine, check = FALSE)
    }, numeric(1))
    sim_rand <- sim_rand[is.finite(sim_rand)]

    sim_mean <- if (length(sim_rand) > 0) mean(sim_rand) else NA_real_
    sim_sd   <- if (length(sim_rand) > 1) stats::sd(sim_rand) else NA_real_
    z_score  <- if (!is.na(sim_sd) && sim_sd > 0 && !is.na(sim_obs)) {
      (sim_obs - sim_mean) / sim_sd
    } else NA_real_
    p_emp <- if (length(sim_rand) > 0 && !is.na(sim_obs)) {
      (1 + sum(sim_rand >= sim_obs)) / (length(sim_rand) + 1)
    } else NA_real_

    tj <- length(intersect(ct$uniprot_id[ct$compound_id == a],
                           ct$uniprot_id[ct$compound_id == b])) /
          length(union(ct$uniprot_id[ct$compound_id == a],
                       ct$uniprot_id[ct$compound_id == b]))
    deg <- if (is.na(sim_obs)) NA_real_ else sim_obs * (1 - tj)

    pair_rows[[i]] <- data.frame(
      condition = cond, compound_a = a, compound_b = b,
      n_targets_a = length(ct$uniprot_id[ct$compound_id == a]),
      n_targets_b = length(ct$uniprot_id[ct$compound_id == b]),
      n_genes_a_mapped = length(ea_all), n_genes_b_mapped = length(eb_all),
      n_pathways_a = NA_integer_, n_pathways_b = NA_integer_,
      target_jaccard = tj, functional_similarity = sim_obs,
      sim_random_mean = sim_mean, sim_random_sd = sim_sd,
      z_score = z_score, p_empirical = p_emp, p_adjusted = NA_real_,
      degeneracy_score = deg,
      stringsAsFactors = FALSE
    )
  }

  pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
  out <- if (length(pair_rows) == 0) {
    .network_degeneracy_direct_empty()
  } else {
    do.call(rbind, pair_rows)
  }
  proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
    message = paste0("condition '", cond, "' (annotation = 'direct'): ", nrow(out),
                     " pair(s) scored, ", n_skipped,
                     " pair(s) skipped (a compound has no Entrez-mapped, annotated target), seed-driven null n_random = ", n_random))
  list(rows = out, proj = proj)
}

.ndeg_old_degeneracy_pathway <- function(proj, cond, ct, all_compounds, annotation,
                                         sem_data, measure, combine, pathway_db) {
  target_pathway <- .network_target_pathway_edges(proj, cond, unique(ct$uniprot_id), pathway_db = pathway_db)
  cp <- if (!is.null(target_pathway) && nrow(target_pathway) > 0) {
    unique(merge(ct, target_pathway, by = "uniprot_id")[, c("compound_id", "pathway_id")])
  } else {
    data.frame(compound_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE)
  }
  target_sets <- split(ct$uniprot_id, ct$compound_id)
  pathway_sets <- split(cp$pathway_id, cp$compound_id)

  if (annotation == "enriched" && nrow(cp) > 0 && !any(grepl("^GO:", cp$pathway_id))) {
    cli::cli_warn(c(
      "!" = "{.fn network_degeneracy}: condition {.val {cond}} (annotation = \"enriched\") has no GO-prefixed enriched pathway IDs.",
      "i" = "{.field functional_similarity} needs GO terms ({.fn GOSemSim::mgoSim}); every pair in this condition will score {.val NA}. Re-run {.fn network_enrich} with a GO database, or use {.code annotation = \"direct\"}."
    ))
  }

  pairs <- utils::combn(all_compounds, 2, simplify = FALSE)
  pair_rows <- vector("list", length(pairs))
  n_skipped <- 0L

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1]; b <- pairs[[i]][2]
    ta <- target_sets[[a]]; tb <- target_sets[[b]]
    pa <- pathway_sets[[a]]; pb <- pathway_sets[[b]]
    if (is.null(pa)) pa <- character(0)
    if (is.null(pb)) pb <- character(0)

    if (length(union(pa, pb)) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }
    target_jaccard <- length(intersect(ta, tb)) / length(union(ta, tb))

    fs <- if (annotation == "jaccard") {
      length(intersect(pa, pb)) / length(union(pa, pb))
    } else {
      ## enriched: GO semantic similarity over the enriched-term sets
      go_a <- unique(pa[grepl("^GO:", pa)])
      go_b <- unique(pb[grepl("^GO:", pb)])
      if (length(go_a) == 0 || length(go_b) == 0) NA_real_
      else GOSemSim::mgoSim(go_a, go_b, semData = sem_data, measure = measure, combine = combine)
    }
    deg <- if (is.na(fs)) NA_real_ else fs * (1 - target_jaccard)

    pair_rows[[i]] <- data.frame(
      condition = cond, compound_a = a, compound_b = b,
      n_targets_a = length(ta), n_targets_b = length(tb),
      n_genes_a_mapped = NA_integer_, n_genes_b_mapped = NA_integer_,
      n_pathways_a = length(pa), n_pathways_b = length(pb),
      target_jaccard = target_jaccard, functional_similarity = fs,
      sim_random_mean = NA_real_, sim_random_sd = NA_real_,
      z_score = NA_real_, p_empirical = NA_real_, p_adjusted = NA_real_,
      degeneracy_score = deg,
      stringsAsFactors = FALSE
    )
  }

  pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
  out <- if (length(pair_rows) == 0) .network_degeneracy_direct_empty() else do.call(rbind, pair_rows)
  proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
    message = paste0("condition '", cond, "' (annotation = '", annotation, "'): ", nrow(out),
                     " pair(s) scored, ", n_skipped,
                     " pair(s) skipped (neither compound has a target in a significant pathway)"))
  list(rows = out, proj = proj, n_skipped = n_skipped)
}

.ndeg_old_go_bma <- function(terms_a, terms_b, sim_mat, combine = "BMA", check = TRUE) {
  if (check) {
    terms_a <- terms_a[!is.na(terms_a) & terms_a %in% rownames(sim_mat)]
    terms_b <- terms_b[!is.na(terms_b) & terms_b %in% colnames(sim_mat)]
  }
  if (length(terms_a) == 0 || length(terms_b) == 0) return(NA_real_)
  S <- sim_mat[terms_a, terms_b, drop = FALSE]
  if (!sum(!is.na(S))) return(NA_real_)

  .cs <- function(M) {
    if (is.vector(M) || nrow(M) == 1 || ncol(M) == 1) {
      if (combine == "avg") return(round(mean(M, na.rm = TRUE), 3))
      return(round(max(M, na.rm = TRUE), 3))
    }
    na <- is.na(M)
    rna <- rowSums(na) == ncol(M)
    if (any(rna)) { M <- M[!rna, , drop = FALSE]; na <- na[!rna, , drop = FALSE] }
    cna <- colSums(na) == nrow(M)
    if (any(cna)) M <- M[, !cna, drop = FALSE]
    if (is.vector(M) || nrow(M) == 1 || ncol(M) == 1) {
      if (combine == "avg") return(round(mean(M, na.rm = TRUE), 3))
      return(round(max(M, na.rm = TRUE), 3))
    }
    result <- switch(combine,
      avg   = mean(M, na.rm = TRUE),
      max   = max(M, na.rm = TRUE),
      rcmax = max(mean(apply(M, 1, max, na.rm = TRUE)),
                  mean(apply(M, 2, max, na.rm = TRUE))),
      ## BMA (rcmax.avg)
      sum(apply(M, 1, max, na.rm = TRUE), apply(M, 2, max, na.rm = TRUE)) / sum(dim(M))
    )
    round(result, 3)
  }
  .cs(S)
}

## Synthetic semantic scores include duplicate annotations, missing similarities,
## a wholly unscorable term, and an empty gene annotation set.
.ndeg_synthetic <- function() {
  terms <- paste0('GO:', seq_len(7))
  sim <- outer(seq_len(7), seq_len(7), function(a, b) 1 / (1 + abs(a - b)))
  dimnames(sim) <- list(terms, terms)
  sim[2, 4] <- sim[4, 2] <- NA_real_
  sim[7, ] <- sim[, 7] <- NA_real_
  g2go <- list(g1 = terms[c(1, 2)], g2 = terms[c(2, 3, 4)],
               g3 = terms[5], g4 = terms[c(4, 6)], g5 = terms[7],
               g6 = character(), g7 = terms[c(1, 4, 5, 6)], g8 = terms[2])
  mapping <- list(u1 = c('g1', 'g2'), u2 = 'g3', u3 = c('g2', 'g4'),
                  u4 = 'g5', u5 = 'unannotated', u6 = character(), u7 = 'g6')
  ct <- data.frame(compound_id = c('A', 'A', 'B', 'B', 'C', 'D', 'E', 'F', 'G'),
                   uniprot_id = c('u1', 'u2', 'u2', 'u3', 'u4', 'u5', 'u6', 'u7', 'u1'))
  list(sim = sim, mapping = mapping,
       args = list(proj = list(), cond = 'synthetic', ct = ct,
                   all_compounds = sort(unique(ct$compound_id)), sem_data = NULL,
                   g2go = g2go, ann_counts = lengths(g2go), measure = 'Wang',
                   combine = 'BMA', universe = 'project', n_random = 70,
                   project_pool_entrez = names(g2go), edges_all = ct))
}

## Give the frozen functions the package helpers but the frozen GO reducer.
.ndeg_oracle <- function(fun) {
  env <- new.env(parent = asNamespace('patliR'))
  env$.network_go_bma <- .ndeg_old_go_bma
  environment(fun) <- env
  fun
}

test_that('optimized direct core is identical to the frozen core, including RNG', {
  skip_if_not_installed('GOSemSim')
  fixture <- .ndeg_synthetic()
  calls <- 0L
  local_mocked_bindings(
    termSim = function(t1, t2, semData, method) {
      calls <<- calls + 1L
      fixture$sim[t1, t2, drop = FALSE]
    }, .package = 'GOSemSim'
  )
  local_mocked_bindings(
    .network_uniprot_to_entrez_map = function(u) fixture$mapping[u],
    .log_append = function(proj, ...) proj
  )
  old <- .ndeg_oracle(.ndeg_old_degeneracy_direct)
  for (combine in c('BMA', 'max', 'avg', 'rcmax')) {
    for (n_random in c(1, 70)) {
      args <- fixture$args
      args$combine <- combine
      args$n_random <- n_random
      set.seed(492)
      expected <- suppressWarnings(do.call(old, args))
      expected_rng <- .Random.seed
      set.seed(492)
      calls <- 0L
      actual <- suppressWarnings(do.call(patliR:::.network_degeneracy_direct, args))
      expect_identical(actual, expected)
      expect_identical(.Random.seed, expected_rng)
      expect_identical(calls, 1L)
      expect_identical(p.adjust(actual$rows$p_empirical, 'BH'),
                       p.adjust(expected$rows$p_empirical, 'BH'))
    }
  }
  ## Capping consumes random numbers before the pair loop; preserve that too.
  extra <- paste0('extra', seq_len(1500))
  fixture$args$g2go[extra] <- rep(list('GO:1'), length(extra))
  fixture$args$ann_counts <- lengths(fixture$args$g2go)
  fixture$args$project_pool_entrez <- names(fixture$args$g2go)
  fixture$args$n_random <- 3
  set.seed(492)
  expected <- suppressWarnings(do.call(old, fixture$args))
  expected_rng <- .Random.seed
  set.seed(492)
  actual <- suppressWarnings(do.call(patliR:::.network_degeneracy_direct, fixture$args))
  expect_identical(actual, expected)
  expect_identical(.Random.seed, expected_rng)
})

test_that('compiled bin draws preserve multi-bin sampling and replacement RNG', {
  bins <- list(c(1L, 4L, 7L), c(2L, 5L, 8L), c(3L, 6L))
  bin_of_node <- c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L)
  for (idx in list(c(8L, 1L, 3L, 4L), c(3L, 6L, 3L))) {
    draw <- patliR:::.network_degeneracy_sampler(idx, bins, bin_of_node)
    set.seed(932)
    expected <- replicate(20, patliR:::.network_resample_matched_idx(idx, bins, bin_of_node),
                          simplify = FALSE)
    expected_rng <- .Random.seed
    set.seed(932)
    expect_identical(replicate(20, draw(), simplify = FALSE), expected)
    expect_identical(.Random.seed, expected_rng)
  }
})

test_that('batched term maxima preserve the pooled-term score and NA edge cases', {
  fixture <- .ndeg_synthetic()
  sim <- fixture$sim
  ## Deliberately asymmetric too: forward/reverse caches must remain distinct.
  sim[3, 1] <- 0.3145
  gene_terms <- lapply(fixture$args$g2go, match, table = rownames(sim))
  maxima <- patliR:::.network_go_gene_maxima(gene_terms, sim)
  sets <- list(1L, 2L, 3L, 5L, 6L, c(1L, 2L), c(1L, 1L, 7L), c(3L, 5L))
  pair <- expand.grid(a = seq_along(sets), b = seq_along(sets))
  a <- sets[pair$a]; b <- sets[pair$b]
  ## Cross the 64-draw batch boundary.
  a <- c(a, a); b <- c(b, b)
  for (combine in c('BMA', 'max', 'avg', 'rcmax')) {
    expected <- vapply(seq_along(a), function(i) {
      .ndeg_old_go_bma(unlist(gene_terms[a[[i]]], use.names = FALSE),
                       unlist(gene_terms[b[[i]]], use.names = FALSE),
                       sim, combine, check = FALSE)
    }, numeric(1))
    expect_identical(patliR:::.network_go_gene_scores(a, b, gene_terms, sim, maxima, combine),
                     expected)
  }
  ## Singleton term universe and a wholly empty batch.
  one <- matrix(0.1235, 1, 1)
  gt <- list(1L)
  expect_identical(patliR:::.network_go_gene_scores(list(1L), list(1L), gt, one,
                     patliR:::.network_go_gene_maxima(gt, one), 'BMA'), round(0.1235, 3))
  expect_identical(patliR:::.network_go_gene_scores(list(6L), list(6L), gene_terms,
                     sim, maxima, 'BMA'), NA_real_)
})

test_that('enriched core reuses one term matrix and preserves mgoSim results', {
  skip_if_not_installed('GOSemSim')
  fixture <- .ndeg_synthetic()
  calls <- 0L
  local_mocked_bindings(
    termSim = function(t1, t2, semData, method) {
      calls <<- calls + 1L
      fixture$sim[t1, t2, drop = FALSE]
    }, .package = 'GOSemSim'
  )
  local_mocked_bindings(
    .network_target_pathway_edges = function(...) {
      data.frame(uniprot_id = c('u1', 'u1', 'u2', 'u3', 'u4', 'u5'),
                 pathway_id = c('GO:1', 'GO:2', 'GO:3', 'GO:4', 'GO:7', 'hsa00010'))
    },
    .log_append = function(proj, ...) proj
  )
  args <- fixture$args[c('proj', 'cond', 'ct', 'all_compounds', 'sem_data', 'measure', 'combine')]
  args["pathway_db"] <- list(NULL)
  old <- .ndeg_oracle(.ndeg_old_degeneracy_pathway)
  for (annotation in c('enriched', 'jaccard')) {
    args$annotation <- annotation
    for (combine in c('BMA', 'max', 'avg', 'rcmax')) {
      args$combine <- combine
      calls <- 0L
      expected <- do.call(old, args)
      old_calls <- calls
      calls <- 0L
      actual <- do.call(patliR:::.network_degeneracy_pathway, args)
      expect_identical(actual, expected)
      expect_identical(calls, if (annotation == 'enriched') 1L else 0L)
      if (annotation == 'enriched') expect_gt(old_calls, calls)
    }
  }
})
