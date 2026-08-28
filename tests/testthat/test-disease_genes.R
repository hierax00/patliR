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

test_that("disease_genes_fetch() validates its arguments", {
  proj <- .test_project()
  expect_error(disease_genes_fetch(proj), "disease")
  expect_error(disease_genes_fetch(proj, disease = ""), "disease")
  expect_error(disease_genes_fetch(proj, disease = "x", min_score = 2), "min_score")
  expect_error(disease_genes_fetch(proj, disease = "x", source = "genecards"), "arg")
})
