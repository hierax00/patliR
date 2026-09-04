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
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)
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
