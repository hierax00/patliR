test_that("disease_genes_import() round-trips a small curated frame into the disease_genes slot", {
  ## Spec 2.0 test (d).
  proj <- .test_project()
  curated <- data.frame(
    uniprot_id = c("P37840", "P05067", "P10636"),
    gene_symbol = c("SNCA", "APP", "MAPT"),
    stringsAsFactors = FALSE
  )
  proj <- disease_genes_import(
    proj, curated,
    disease_id = "MONDO_0005180", disease_name = "Parkinson disease",
    source = "curated_demo"
  )
  dg <- patliRResults(proj, "disease_genes")

  expect_setequal(names(dg), names(patliR:::.empty_disease_genes_row()))
  expect_setequal(dg$uniprot_id, curated$uniprot_id)
  expect_true(all(dg$disease_id == "MONDO_0005180"))
  expect_true(all(dg$disease_name == "Parkinson disease"))
  expect_true(all(dg$source == "curated_demo"))
  expect_true(all(is.na(dg$association_score)))
  expect_true(file.exists(file.path(projectDir(proj), "results", "disease_genes.csv")))

  ## survives a reload from CSV
  reloaded <- patliR_load(projectDir(proj))
  expect_setequal(patliRResults(reloaded, "disease_genes")$uniprot_id, curated$uniprot_id)
})

test_that("disease_genes_import() reads a curated CSV path and keeps an association score column", {
  proj <- .test_project()
  f <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(uniprot = c("P12345", "Q99999"), score = c(0.8, 0.2)),
    f, row.names = FALSE
  )
  proj <- disease_genes_import(proj, f, disease_id = "D1", source = "GWAS Catalog")
  dg <- patliRResults(proj, "disease_genes")
  expect_setequal(dg$uniprot_id, c("P12345", "Q99999"))
  expect_equal(dg$association_score[dg$uniprot_id == "P12345"], 0.8)
})

test_that("disease_genes_import() rejects a symbol-only table with a clear message", {
  proj <- .test_project()
  expect_error(
    disease_genes_import(proj, data.frame(symbol = c("SNCA", "APP")), disease_id = "D1"),
    "UniProt"
  )
})

test_that("disease_genes_import() requires a non-empty disease_id", {
  proj <- .test_project()
  expect_error(
    disease_genes_import(proj, data.frame(uniprot_id = "P12345")),
    "disease_id"
  )
  expect_error(
    disease_genes_import(proj, data.frame(uniprot_id = "P12345"), disease_id = ""),
    "disease_id"
  )
})

test_that("disease_genes_import() upsert replaces only the given disease's rows", {
  proj <- .test_project()
  proj <- disease_genes_import(proj, data.frame(uniprot_id = c("P1", "P2")),
                               disease_id = "DA", source = "s")
  proj <- disease_genes_import(proj, data.frame(uniprot_id = c("P3")),
                               disease_id = "DB", source = "s")
  proj <- disease_genes_import(proj, data.frame(uniprot_id = c("P9")),
                               disease_id = "DA", source = "s")   # re-import DA
  dg <- patliRResults(proj, "disease_genes")
  expect_setequal(dg$uniprot_id[dg$disease_id == "DA"], "P9")
  expect_setequal(dg$uniprot_id[dg$disease_id == "DB"], "P3")
})

test_that("disease_genes_fetch() maps Swiss-Prot proteinIds and upserts by disease_id (GraphQL mocked)", {
  ## Spec 2.0 test (e): the fetch's GraphQL call is fully mocked -- no live
  ## network in the suite.
  fake_graphql <- function(query_string, variables = list()) {
    if (grepl("DiseaseById", query_string)) {
      return(list(disease = list(id = variables$id)))
    }
    if (grepl("DiseaseTargets", query_string)) {
      if (identical(variables$index, 0L) || identical(variables$index, 0)) {
        return(list(disease = list(
          id = "EFO_TEST", name = "testitis",
          associatedTargets = list(
            count = 3,
            rows = list(
              list(score = 0.90, target = list(
                id = "ENSG1", approvedSymbol = "AAA",
                proteinIds = list(
                  list(id = "P00001", source = "uniprot_swissprot"),
                  list(id = "ZZZ", source = "uniprot_trembl")
                ))),
              list(score = 0.40, target = list(
                id = "ENSG2", approvedSymbol = "BBB",
                proteinIds = list(
                  list(id = "P00002", source = "uniprot_swissprot"),
                  list(id = "P00003", source = "uniprot_swissprot")
                ))),
              list(score = 0.05, target = list(
                id = "ENSG3", approvedSymbol = "CCC",
                proteinIds = list(list(id = "ENSP9", source = "ensembl_PRO"))))
            )
          )
        )))
      }
      return(list(disease = list(id = "EFO_TEST", name = "testitis",
                                 associatedTargets = list(count = 3, rows = list()))))
    }
    stop("unexpected query in mock: ", query_string)
  }
  testthat::local_mocked_bindings(.open_targets_graphql = fake_graphql, .package = "patliR")

  proj <- .test_project()
  proj <- disease_genes_fetch(proj, disease = "EFO_TEST")
  dg <- patliRResults(proj, "disease_genes")

  expect_setequal(names(dg), names(patliR:::.empty_disease_genes_row()))
  ## ENSG2's two Swiss-Prot accessions both kept; ENSG3 (no swissprot) dropped
  expect_setequal(dg$uniprot_id, c("P00001", "P00002", "P00003"))
  expect_false("CCC" %in% dg$gene_symbol)
  expect_true(all(dg$disease_id == "EFO_TEST"))
  expect_true(all(dg$source == "open_targets"))
  expect_equal(unique(dg$association_score[dg$gene_symbol == "BBB"]), 0.40)

  ## min_score is applied and logged
  proj2 <- disease_genes_fetch(proj, disease = "EFO_TEST", min_score = 0.5)
  dg2 <- patliRResults(proj2, "disease_genes")
  expect_setequal(dg2$uniprot_id, "P00001")
  expect_true(any(grepl("disease_genes_below_min_score", projectLog(proj2)$message)))
})

test_that("disease_genes_fetch() re-fetch replaces only that disease's rows (GraphQL mocked)", {
  fake_graphql <- function(query_string, variables = list()) {
    if (grepl("DiseaseById", query_string)) return(list(disease = list(id = variables$id)))
    list(disease = list(
      id = variables$efoId, name = "d",
      associatedTargets = if (identical(as.integer(variables$index), 0L)) list(
        count = 1,
        rows = list(list(score = 0.7, target = list(
          id = "ENSGX", approvedSymbol = "XXX",
          proteinIds = list(list(id = "P77777", source = "uniprot_swissprot")))))
      ) else list(count = 1, rows = list())
    ))
  }
  testthat::local_mocked_bindings(.open_targets_graphql = fake_graphql, .package = "patliR")

  proj <- .test_project()
  proj <- disease_genes_import(proj, data.frame(uniprot_id = "P00000"),
                               disease_id = "KEEP_ME", source = "curated")
  proj <- disease_genes_fetch(proj, disease = "EFO_AAA")
  proj <- disease_genes_fetch(proj, disease = "EFO_AAA")   # idempotent re-fetch
  dg <- patliRResults(proj, "disease_genes")
  expect_true("P00000" %in% dg$uniprot_id)                                  # curated disease untouched
  expect_equal(sum(dg$disease_id == "EFO_AAA"), 1L)                         # not duplicated
})

test_that(".open_targets_disease_targets() paginates a partial final page, each row once (GraphQL mocked)", {
  ## Test gap: the fetch mocks only ever return "page 0 -> rows, page >= 1
  ## -> empty", never a partial final page (the real shape: T2D page 19
  ## returns 407 of 500). count 7 / page_size 5 exercises the
  ## `n_seen >= count` arithmetic and the "no wasted extra request" path.
  all_targets <- lapply(1:7, function(k) list(
    score = k / 10,
    target = list(
      id = paste0("ENSG", k), approvedSymbol = paste0("G", k),
      proteinIds = list(list(id = sprintf("P%05d", k), source = "uniprot_swissprot"))
    )
  ))
  seen <- new.env(parent = emptyenv()); seen$idx <- integer(0)
  fake_graphql <- function(query_string, variables = list()) {
    i <- as.integer(variables$index); size <- as.integer(variables$size)
    seen$idx <- c(seen$idx, i)
    start <- i * size + 1L
    rows <- if (start > 7L) list() else all_targets[start:min(start + size - 1L, 7L)]
    list(disease = list(id = variables$efoId, name = "d",
                        associatedTargets = list(count = 7L, rows = rows)))
  }
  testthat::local_mocked_bindings(.open_targets_graphql = fake_graphql, .package = "patliR")

  out <- patliR:::.open_targets_disease_targets("EFO_PAGED", page_size = 5L)
  expect_equal(nrow(out$rows), 7L)
  expect_setequal(out$rows$uniprot_id, sprintf("P%05d", 1:7))
  expect_false(anyDuplicated(out$rows$uniprot_id) > 0)
  expect_equal(seen$idx, c(0L, 1L))   # exactly two requests: full page + partial page
})

test_that("disease_genes_import(map_symbols = TRUE) maps a symbol-only table via org.Hs.eg.db", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  proj <- .test_project()
  proj <- disease_genes_import(
    proj, data.frame(symbol = c("SNCA", "APP", "MAPT")),
    disease_id = "MONDO_0005180", source = "curated", map_symbols = TRUE
  )
  dg <- patliRResults(proj, "disease_genes")
  expect_gt(nrow(dg), 0)
  expect_true(all(grepl("^[A-Z0-9]+$", dg$uniprot_id)))
  expect_true(all(dg$gene_symbol %in% c("SNCA", "APP", "MAPT")))
  expect_true(all(grepl("org.Hs.eg.db", dg$source)))
  expect_true(any(grepl("disease_genes_import_map_symbols", projectLog(proj)$message)))
})

test_that("disease_genes_import(map_symbols = FALSE) still rejects a symbol-only table", {
  proj <- .test_project()
  expect_error(
    disease_genes_import(proj, data.frame(symbol = c("SNCA", "APP")), disease_id = "D1"),
    "UniProt"
  )
})

test_that("disease_genes_fetch() validates its arguments", {
  proj <- .test_project()
  expect_error(disease_genes_fetch(proj), "disease")
  expect_error(disease_genes_fetch(proj, disease = ""), "disease")
  expect_error(disease_genes_fetch(proj, disease = "x", min_score = 2), "min_score")
  expect_error(disease_genes_fetch(proj, disease = "x", source = "genecards"), "arg")
})

test_that("disease_genes_fetch() requires a known score to meet min_score", {
  fake_graphql <- function(query_string, variables = list()) {
    if (grepl("DiseaseById", query_string)) return(list(disease = list(id = variables$id)))
    scores <- c(NA_real_, 0.4, 0.39)
    list(disease = list(
      id = variables$efoId, name = "test disease",
      associatedTargets = list(count = 3L, rows = lapply(seq_along(scores), function(i) {
        list(score = scores[i], target = list(
          id = paste0("ENSG", i), approvedSymbol = paste0("G", i),
          proteinIds = list(list(id = paste0("P", i), source = "uniprot_swissprot"))
        ))
      }))
    ))
  }
  testthat::local_mocked_bindings(.open_targets_graphql = fake_graphql, .package = "patliR")
  proj <- .test_project()
  proj <- disease_genes_import(proj, data.frame(uniprot_id = "P_OTHER"), disease_id = "OTHER")
  proj <- disease_genes_fetch(proj, disease = "EFO_TEST")
  dg <- patliRResults(proj, "disease_genes")
  expect_equal(dg$uniprot_id[dg$disease_id == "EFO_TEST"], "P2")
  expect_equal(dg$association_score[dg$disease_id == "EFO_TEST"], 0.4)

  proj <- disease_genes_fetch(proj, disease = "EFO_TEST", min_score = NULL)
  dg <- patliRResults(proj, "disease_genes")
  expect_setequal(dg$uniprot_id[dg$disease_id == "EFO_TEST"], c("P1", "P2", "P3"))
  expect_true(is.na(dg$association_score[dg$uniprot_id == "P1"]))

  proj <- disease_genes_fetch(proj, disease = "EFO_TEST", min_score = 0.5)
  dg <- patliRResults(proj, "disease_genes")
  expect_equal(dg$uniprot_id, "P_OTHER")
  expect_true(any(grepl("dropped 3", projectLog(proj)$message)))
})

test_that("disease_genes_import() keeps scores aligned through symbol multi-mapping", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::local_mocked_bindings(
    bitr = function(...) data.frame(
      SYMBOL = c("AAA", "AAA", "BBB"), UNIPROT = c("P1", "P2", "P3")
    ),
    .package = "clusterProfiler"
  )
  proj <- .test_project()
  proj <- disease_genes_import(
    proj, data.frame(symbol = c("BBB", " AAA ", "AAA", "UNMAPPED"),
                     score = c(0.2, 0.8, 0.6, 0.9)),
    disease_id = "D1", map_symbols = TRUE
  )
  dg <- patliRResults(proj, "disease_genes")
  expect_equal(nrow(dg), 5L)
  expect_setequal(dg$association_score[dg$uniprot_id == "P1"], c(0.8, 0.6))
  expect_setequal(dg$association_score[dg$uniprot_id == "P2"], c(0.8, 0.6))
  expect_equal(dg$association_score[dg$uniprot_id == "P3"], 0.2)
  expect_false("UNMAPPED" %in% dg$gene_symbol)
})
