## Unit tests for .network_pathview_gene_vector() (spec 1.12 [SHOULD] #2:
## aggregate()'s implicit na.action = na.omit silently dropping NA-weight
## edges) do not need pathview/clusterProfiler/STRINGdb at all -- pure data
## frame in, list(proj, gene_vector) out.

test_that(".network_pathview_gene_vector() drops NA-weight edges explicitly and logs which Entrez IDs lost one", {
  proj <- .test_project()
  edges_cond <- data.frame(
    ENTREZID = c("1", "1", "2"), weight = c(0.5, NA_real_, 0.8), stringsAsFactors = FALSE
  )
  out <- patliR:::.network_pathview_gene_vector(proj, "COND", edges_cond, "max_weight")

  expect_equal(unname(out$gene_vector["1"]), 0.5) ## the surviving (non-NA) row for ENTREZID "1"
  expect_equal(unname(out$gene_vector["2"]), 0.8)

  msgs <- projectLog(out$proj)$message
  expect_true(any(grepl("network_pathview_na_weight_dropped", msgs)))
})

test_that(".network_pathview_gene_vector() aborts with a clear message when every weight is NA (would otherwise hand pathview() an empty vector)", {
  proj <- .test_project()
  edges_cond <- data.frame(
    ENTREZID = c("1", "2"), weight = c(NA_real_, NA_real_), stringsAsFactors = FALSE
  )
  expect_error(
    patliR:::.network_pathview_gene_vector(proj, "COND", edges_cond, "mean_weight"),
    "NA"
  )
})

test_that(".network_pathview_gene_vector() n_compounds branch ignores weight entirely, NA or not", {
  proj <- .test_project()
  edges_cond <- data.frame(
    ENTREZID = c("1", "1", "2"), compound_id = c("C1", "C2", "C1"),
    weight = c(NA_real_, NA_real_, NA_real_), stringsAsFactors = FALSE
  )
  out <- patliR:::.network_pathview_gene_vector(proj, "COND", edges_cond, "n_compounds")
  expect_equal(unname(out$gene_vector["1"]), 2L)
  expect_equal(unname(out$gene_vector["2"]), 1L)
})

## End-to-end network_pathview() tests need pathview/clusterProfiler/
## org.Hs.eg.db, but never touch the network (pathview::pathview() and
## clusterProfiler::bitr() are mocked) or the filesystem outside a temp dir.

test_that("network_pathview(): ok reflects whether pathview() actually wrote the PNG, not just the absence of an R error", {
  testthat::skip_if_not_installed("pathview")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  uniprot_ids <- unique(edges$uniprot_id[edges$condition == cond])
  testthat::skip_if(length(uniprot_ids) < 1, "need >= 1 target in the FLO-ET fixture")

  entrez_map <- data.frame(
    UNIPROT = uniprot_ids, ENTREZID = paste0("E", seq_along(uniprot_ids)), stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = cond, db = "kegg", ID = c("hsaOK", "hsaFAIL"),
    Description = c("ok pathway", "fail pathway"),
    GeneRatio = "1/1", BgRatio = "1/1", pvalue = 0.01, p.adjust = c(0.01, 0.02), qvalue = 0.01,
    geneID = "1", Count = 1L, stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(bitr = function(x, ...) entrez_map, .package = "clusterProfiler")
  ## Real pathview(): frequently warns and returns without writing anything
  ## (no mappable nodes / a KGML-download that returns an HTML error page)
  ## -- no R `error` condition either way. Mimic that: "hsaOK" writes the
  ## PNG pathview() would have written (into the current wd -- the caller
  ## setwd()s into out_dir first); "hsaFAIL" returns cleanly without
  ## writing anything.
  testthat::local_mocked_bindings(
    pathview = function(gene.data, pathway.id, out.suffix, ...) {
      if (identical(pathway.id, "hsaOK")) file.create(paste0(pathway.id, ".", out.suffix, ".png"))
      invisible(NULL)
    },
    .package = "pathview"
  )

  out_dir <- file.path(tempdir(), paste0("pathview_test_ok_", as.integer(Sys.time())))
  proj <- network_pathview(proj, condition = cond, pathway_id = c("hsaOK", "hsaFAIL"), out_dir = out_dir)
  log <- patliRResults(proj, "kegg_pathview_log")

  expect_true(log$ok[log$pathway_id == "hsaOK"])
  expect_true(file.exists(log$path[log$pathway_id == "hsaOK"]))

  expect_false(log$ok[log$pathway_id == "hsaFAIL"])
  expect_true(is.na(log$path[log$pathway_id == "hsaFAIL"]))
  expect_true(any(grepl("network_pathview_render_failed", projectLog(proj)$message)))
})

test_that("network_pathview() warns when a gene score exceeds 1 before pathview()'s [0, 1] colour clamp", {
  testthat::skip_if_not_installed("pathview")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  cond <- "FLO-ET"
  edges <- patliRResults(proj, "network_edges")
  uniprot_ids <- unique(edges$uniprot_id[edges$condition == cond])
  testthat::skip_if(length(uniprot_ids) < 1, "need >= 1 target in the FLO-ET fixture")

  ## Simulate a platform that reported a raw 0-100 percentage into `weight`
  ## instead of targets_import()'s [0, 1] convention.
  idx <- which(edges$condition == cond)[1]
  edges$weight[idx] <- 96.5
  patliRResults(proj, "network_edges") <- edges

  entrez_map <- data.frame(
    UNIPROT = uniprot_ids, ENTREZID = paste0("E", seq_along(uniprot_ids)), stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = cond, db = "kegg", ID = "hsaOK", Description = "ok pathway",
    GeneRatio = "1/1", BgRatio = "1/1", pvalue = 0.01, p.adjust = 0.01, qvalue = 0.01,
    geneID = "1", Count = 1L, stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(bitr = function(x, ...) entrez_map, .package = "clusterProfiler")
  testthat::local_mocked_bindings(
    pathview = function(gene.data, pathway.id, out.suffix, ...) {
      file.create(paste0(pathway.id, ".", out.suffix, ".png"))
      invisible(NULL)
    },
    .package = "pathview"
  )

  out_dir <- file.path(tempdir(), paste0("pathview_test_clamp_", as.integer(Sys.time())))
  expect_warning(
    network_pathview(proj, condition = cond, pathway_id = "hsaOK", out_dir = out_dir),
    "exceed 1"
  )
})

test_that("network_pathview() resolves relative directories and isolates fresh renders by condition", {
  testthat::skip_if_not_installed("pathview")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  conds <- unique(edges$condition)[1:2]
  testthat::skip_if(anyNA(conds), "need two conditions")
  uniprot_ids <- unique(edges$uniprot_id)
  entrez_map <- data.frame(
    UNIPROT = uniprot_ids, ENTREZID = paste0("E", seq_along(uniprot_ids)), stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_enrichment") <- data.frame(
    condition = conds, db = "kegg", ID = "hsaOK", Description = "ok pathway",
    GeneRatio = "1/1", BgRatio = "1/1", pvalue = 0.01, p.adjust = 0.01, qvalue = 0.01,
    geneID = "1", Count = 1L, stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(bitr = function(x, ...) entrez_map, .package = "clusterProfiler")
  render_text <- "first condition"
  write_png <- TRUE
  root <- tempfile("pathview_relative_")
  dir.create(root)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root)
  root_wd <- getwd()
  testthat::local_mocked_bindings(
    pathview = function(gene.data, pathway.id, out.suffix, kegg.dir, ...) {
      expect_identical(kegg.dir, normalizePath(file.path(root, "cache"), winslash = "/"))
      if (write_png) writeLines(render_text, paste0(pathway.id, ".", out.suffix, ".png"))
      invisible(NULL)
    },
    .package = "pathview"
  )

  proj <- network_pathview(proj, condition = conds[1], out_dir = "plots", kegg_dir = "cache")
  first <- patliRResults(proj, "kegg_pathview_log")
  expect_true(first$ok)
  expect_true(file.exists(first$path))
  expect_identical(getwd(), root_wd)

  render_text <- "second condition"
  proj <- network_pathview(proj, condition = conds[2], out_dir = "plots", kegg_dir = "cache")
  log <- patliRResults(proj, "kegg_pathview_log")
  expect_equal(length(unique(log$path)), 2L)
  expect_identical(readLines(first$path), "first condition")
  expect_true(all(log$ok))

  write_png <- FALSE
  proj <- network_pathview(proj, condition = conds[2], out_dir = "plots", kegg_dir = "cache")
  log <- patliRResults(proj, "kegg_pathview_log")
  expect_false(log$ok[log$condition == conds[2]])
  expect_true(is.na(log$path[log$condition == conds[2]]))
  expect_identical(readLines(first$path), "first condition")
})
