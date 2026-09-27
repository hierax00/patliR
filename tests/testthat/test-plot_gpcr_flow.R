## plot_gpcr_flow() is exercised fully offline: .ppi_gpcr_go() (org.Hs.eg.db)
## is mocked, exactly as test-plot_ppi_network.R does; the G-protein coupling
## table itself is the REAL bundled CSV, keyed by real UniProt accessions, so
## the merge/drop logic is tested against the actual shipped data.

.gflow_fake_gpcr <- function() {
  ## ADRB1 (Gs), CNR2 (Gi/o), S1PR2 (Gq/11), GPR55 (G12/13), HPGD (known
  ## false positive -- g_protein is NA in the real CSV) and one accession
  ## absent from the CSV entirely (also dropped, but for a different reason)
  data.frame(
    uniprot_id = c("P08588", "P34972", "O95136", "Q9Y2T6", "P15428", "X99999"),
    symbol = c("ADRB1", "CNR2", "S1PR2", "GPR55", "HPGD", "MADEUP"),
    evidence = "IDA", stringsAsFactors = FALSE
  )
}

.gflow_test_project <- function() {
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "A",
    compound_id = c("c1", "c1", "c1", "c2", "c2", "c3", "c3"),
    uniprot_id = c("P08588", "P34972", "O95136", "P08588", "Q9Y2T6", "P15428", "X99999"),
    weight = 0.6, stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "test disease",
    uniprot_id = "P08588", gene_symbol = "ADRB1", association_score = 0.9,
    stringsAsFactors = FALSE
  )
  proj
}

test_that("plot_gpcr_flow() drops non-GPCR and uncoupled targets, keeps the rest with their real coupling", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  proj <- .gflow_test_project()

  p <- plot_gpcr_flow(proj, condition = "A", save = FALSE)
  tab <- attr(p, "table")
  ## P15428 (HPGD, false positive) and X99999 (not in the coupling table) are dropped
  expect_setequal(tab$uniprot_id, c("P08588", "P34972", "O95136", "Q9Y2T6"))
  expect_equal(unique(tab$g_protein[tab$uniprot_id == "P08588"]), "Gs")
  expect_equal(unique(tab$g_protein[tab$uniprot_id == "P34972"]), "Gi/o")
  expect_equal(unique(tab$g_protein[tab$uniprot_id == "O95136"]), "Gq/11")
  expect_equal(unique(tab$g_protein[tab$uniprot_id == "Q9Y2T6"]), "G12/13")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_gpcr_flow(disease = ) marks the receptor that is also a disease gene", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  proj <- .gflow_test_project()

  p <- plot_gpcr_flow(proj, condition = "A", disease = "D1", save = FALSE)
  tab <- attr(p, "table")
  expect_true(all(tab$is_disease_gene[tab$uniprot_id == "P08588"]))
  expect_false(any(tab$is_disease_gene[tab$uniprot_id != "P08588"]))
})

test_that("plot_gpcr_flow(top_n_receptors = ) keeps the most-hit receptors", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  proj <- .gflow_test_project()

  p <- plot_gpcr_flow(proj, condition = "A", top_n_receptors = 1, save = FALSE)
  tab <- attr(p, "table")
  expect_equal(unique(tab$uniprot_id), "P08588") # hit by c1 AND c2 -- the only receptor with 2 compounds
})

test_that("plot_gpcr_flow(show_vascular_effect = TRUE) buckets uncovered receptors instead of guessing", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggalluvial")
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  proj <- .gflow_test_project()

  p <- plot_gpcr_flow(proj, condition = "A", show_vascular_effect = TRUE, save = FALSE)
  tab <- attr(p, "table")
  expect_equal(tab$vascular_effect[tab$uniprot_id == "P34972"], "not established / other tissue") # CNR2: NA in the CSV
  expect_false(any(is.na(tab$vascular_effect[tab$uniprot_id == "P08588"])))
})

test_that("plot_gpcr_flow() validates and errors clearly when nothing is a GPCR", {
  skip_if_not_installed("ggplot2")
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  proj <- .gflow_test_project()
  expect_error(plot_gpcr_flow(proj, condition = "A", top_n_receptors = 0), "top_n_receptors")

  proj2 <- .test_project()
  patliRResults(proj2, "network_edges") <- data.frame(condition = "A", compound_id = "c1", uniprot_id = "P99999", weight = 0.5, stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() .gflow_fake_gpcr(), .package = "patliR")
  expect_error(plot_gpcr_flow(proj2, condition = "A", save = FALSE), "classified as a GPCR")
})
