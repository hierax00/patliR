.network_enrich_test_setup <- function() {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  proj <- targets_import_batch(
    proj,
    system.file("extdata", "import_targets", package = "patliR"),
    platform = "superpred"
  )
  proj <- network_build(proj)
  proj
}

test_that("network_enrich() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- prep_binarize(proj, .test_abundance_matrix())
  ## no network_build() -- but skip straight past the dependency check so
  ## this test is meaningful even on a machine without clusterProfiler/org.Hs.eg.db.
  ## db = "go" on purpose (not the "reactome" default): this test isn't about
  ## reactome specifically, and ReactomePA may not be installed even when
  ## clusterProfiler/org.Hs.eg.db are -- that dependency error would otherwise
  ## fire before ever reaching the network_edges check this test targets.
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  expect_error(network_enrich(proj, db = "go"), "network_edges")
})

test_that("network_enrich() errors on a condition network_build() never built", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, condition = "NOT-A-REAL-CONDITION", db = "go"), "not built")
})

test_that("network_enrich() validates pvalueCutoff/qvalueCutoff", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, pvalueCutoff = 0), "pvalueCutoff")
  expect_error(network_enrich(proj, pvalueCutoff = 1.5), "pvalueCutoff")
  expect_error(network_enrich(proj, qvalueCutoff = -0.1), "qvalueCutoff")
})

test_that("network_enrich() runs end-to-end for db = 'go' when dependencies are installed", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  ## simplify_go defaults to TRUE -- this is the one test in the suite that
  ## deliberately pays the GOSemSim::godata() build cost, so the default
  ## path itself stays covered (see test-network_degeneracy.R/
  ## test-network_motifs.R/test-plot_*.R for why they opt out with
  ## simplify_go = FALSE instead of each rebuilding it too).
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go")

  result <- patliRResults(proj, "network_enrichment")
  expect_true(all(c("condition", "db", "ID", "Description", "GeneRatio", "BgRatio",
                     "pvalue", "p.adjust", "qvalue", "geneID", "Count",
                     "ONTOLOGY", "universe") %in% names(result)))
  if (nrow(result) > 0) {
    expect_true(all(result$condition == "FLO-ET"))
    expect_true(all(result$db == "go"))
    ## universe defaults to "project" and is recorded per row
    expect_true(all(result$universe == "project"))
    ## simplify_go ran (not skipped/failed silently) -- either a reduction
    ## message or, for this tiny example, possibly nothing to simplify.
    log_msg <- projectLog(proj)$message[projectLog(proj)$step == "network_enrich"]
    expect_false(any(grepl("network_enrich_simplify_failed", log_msg)))
  }
})

test_that("network_enrich() with simplify_go = FALSE skips GOSemSim entirely", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  proj <- network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE)

  result <- patliRResults(proj, "network_enrichment")
  expect_true(all(c("condition", "db", "ID", "Description") %in% names(result)))
  log_msg <- projectLog(proj)$message[projectLog(proj)$step == "network_enrich"]
  expect_false(any(grepl("simplify_go", log_msg)))
})

test_that("network_enrich() with db = 'reactome' requires ReactomePA specifically", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if(requireNamespace("ReactomePA", quietly = TRUE), "ReactomePA is installed -- nothing to test here")

  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, db = "reactome"), "ReactomePA")
})

test_that("network_enrich() validates pAdjustMethod", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .network_enrich_test_setup()
  expect_error(network_enrich(proj, pAdjustMethod = "not-a-method"), "should be one of")
})

test_that(".network_enrich_result_df() fills NA when qvalue/ONTOLOGY are absent from clusterProfiler's output", {
  ## No Bioconductor dependency -- clusterProfiler omits `qvalue` whenever
  ## the qvalue package's pi0 estimation fails, routine for a small gene
  ## set (audit finding); this exercises the defensive column-by-column
  ## build directly rather than relying on that failure mode reproducing.
  df_no_qvalue <- data.frame(
    ID = c("GO:1", "GO:2"), Description = c("term a", "term b"),
    GeneRatio = c("2/5", "1/5"), BgRatio = c("20/100", "10/100"),
    pvalue = c(0.01, 0.2), p.adjust = c(0.02, 0.3),
    geneID = c("1/2", "3"), Count = c(2L, 1L),
    stringsAsFactors = FALSE
  )
  out <- patliR:::.network_enrich_result_df(df_no_qvalue, cond = "COND1", db = "go", universe = "project")
  expect_true(all(is.na(out$qvalue)))
  expect_true(all(is.na(out$ONTOLOGY)))
  expect_equal(out$condition, rep("COND1", 2))
  expect_equal(out$db, rep("go", 2))
  expect_equal(out$universe, rep("project", 2))
  expect_equal(out$ID, df_no_qvalue$ID)
  expect_equal(out$Count, df_no_qvalue$Count)

  ## and it passes real qvalue/ONTOLOGY values through unchanged when present
  df_with_both <- df_no_qvalue
  df_with_both$qvalue <- c(0.03, 0.4)
  df_with_both$ONTOLOGY <- c("BP", "BP")
  out2 <- patliR:::.network_enrich_result_df(df_with_both, cond = "COND2", db = "go", universe = "genome")
  expect_equal(out2$qvalue, df_with_both$qvalue)
  expect_equal(out2$ONTOLOGY, df_with_both$ONTOLOGY)
  expect_equal(out2$universe, rep("genome", 2))
})

test_that("network_enrich() resolves universe = 'project' to a non-NULL Entrez background and threads it to .network_enrich_run(); universe = 'genome' passes NULL", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()

  captured <- list()
  fake_run <- function(db, entrez_ids, ont, pvalueCutoff, qvalueCutoff, pAdjustMethod, universe_entrez = NULL) {
    ## `captured[[i]] <<- NULL` would DELETE that position instead of
    ## storing it (list[[i]] <- NULL removes an element) -- wrap in list()
    ## and concatenate so a NULL universe_entrez (universe = "genome") is
    ## actually recorded as a captured NULL, not silently dropped.
    captured <<- c(captured, list(universe_entrez))
    NULL
  }
  testthat::local_mocked_bindings(.network_enrich_run = fake_run)

  network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE, universe = "project")
  network_enrich(proj, condition = "FLO-ET", db = "go", simplify_go = FALSE, universe = "genome")

  expect_length(captured, 2)
  expect_true(is.character(captured[[1]]))
  expect_gt(length(captured[[1]]), 0)
  expect_null(captured[[2]])
})

test_that("network_enrich() universe = 'project' background is smaller than universe = 'genome' (BgRatio denominator, real clusterProfiler call)", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  ## pvalueCutoff/qvalueCutoff relaxed to 1 so terms are kept regardless of
  ## significance -- the comparison is about the background size recorded
  ## in BgRatio, not about which terms pass a cutoff (which could
  ## legitimately differ in either direction between the two universes).
  proj_bg <- network_enrich(
    proj, condition = "FLO-ET", db = "go", ont = "MF", simplify_go = FALSE,
    universe = "project", pvalueCutoff = 1, qvalueCutoff = 1
  )
  proj_gen <- network_enrich(
    proj, condition = "FLO-ET", db = "go", ont = "MF", simplify_go = FALSE,
    universe = "genome", pvalueCutoff = 1, qvalueCutoff = 1
  )

  res_bg  <- patliRResults(proj_bg,  "network_enrichment")
  res_gen <- patliRResults(proj_gen, "network_enrichment")
  skip_if(nrow(res_bg) == 0 || nrow(res_gen) == 0,
          "no GO MF terms survive minGSSize for this tiny fixture at one of the two universes")

  expect_true(all(res_bg$universe == "project"))
  expect_true(all(res_gen$universe == "genome"))

  bg_denom  <- as.numeric(sub(".*/", "", res_bg$BgRatio[1]))
  gen_denom <- as.numeric(sub(".*/", "", res_gen$BgRatio[1]))
  expect_lt(bg_denom, gen_denom)
})

test_that("network_enrich() pAdjustMethod is passed through and changes p.adjust ('BH' vs 'bonferroni')", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  proj_bh <- network_enrich(
    proj, condition = "FLO-ET", db = "go", ont = "MF", simplify_go = FALSE,
    universe = "genome", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"
  )
  proj_bonf <- network_enrich(
    proj, condition = "FLO-ET", db = "go", ont = "MF", simplify_go = FALSE,
    universe = "genome", pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "bonferroni"
  )

  res_bh   <- patliRResults(proj_bh,   "network_enrichment")
  res_bonf <- patliRResults(proj_bonf, "network_enrichment")
  skip_if(nrow(res_bh) == 0, "no GO MF terms for this tiny fixture")

  res_bh   <- res_bh[order(res_bh$ID), ]
  res_bonf <- res_bonf[order(res_bonf$ID), ]
  expect_identical(res_bh$ID, res_bonf$ID)
  expect_false(isTRUE(all.equal(res_bh$p.adjust, res_bonf$p.adjust)))
  ## bonferroni is at least as conservative as BH, term by term
  expect_true(all(res_bonf$p.adjust >= res_bh$p.adjust - 1e-12))
})

test_that("network_enrich() migrates a legacy network_enrichment table (no universe/ONTOLOGY column) via .network_upsert()", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  skip_on_cran()

  proj <- .network_enrich_test_setup()
  legacy <- data.frame(
    condition = "LEGACY-COND", db = "go", ID = "GO:0000001", Description = "legacy term",
    GeneRatio = "3/10", BgRatio = "30/18000", pvalue = 0.01, p.adjust = 0.02, qvalue = 0.03,
    geneID = "1/2/3", Count = 3L,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_enrichment") <- legacy

  expect_no_error(
    proj <- network_enrich(
      proj, condition = "FLO-ET", db = "go", simplify_go = FALSE,
      pvalueCutoff = 1, qvalueCutoff = 1
    )
  )
  res <- patliRResults(proj, "network_enrichment")
  expect_true(all(c("universe", "ONTOLOGY") %in% names(res)))

  legacy_row <- res[res$condition == "LEGACY-COND", ]
  expect_equal(nrow(legacy_row), 1L)
  expect_true(is.na(legacy_row$universe))
  expect_true(is.na(legacy_row$ONTOLOGY))

  new_rows <- res[res$condition == "FLO-ET", ]
  skip_if(nrow(new_rows) == 0, "no GO BP terms for this tiny fixture")
  expect_true(all(new_rows$universe == "project"))
})
