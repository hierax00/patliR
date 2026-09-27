.cflow_fixture <- function() {
  proj <- .test_project()
  bin <- data.frame(
    compound_id = c("c1", "c2", "c3", "c4"),
    A = c(1, 1, 0, 0), B = c(1, 0, 1, 0), C = c(0, 1, 1, 1),
    stringsAsFactors = FALSE
  )
  binarizedMatrix(proj) <- bin
  patliRResults(proj, "compounds_classified") <- data.frame(
    compound_id = c("c1", "c2", "c3", "c4"),
    pathway = c("Terpenoids", "Terpenoids", "Fatty acids", "Fatty acids"),
    superclass = c("Sesquiterpenoids", "Monoterpenoids", "Fatty acyls", "Fatty esters"),
    class = c("A1", "A2", "B1", "B2"),
    stringsAsFactors = FALSE
  )
  proj
}

test_that("plot_condition_flow() tracks compound presence across the given condition sequence", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .cflow_fixture()
  p <- plot_condition_flow(proj, conditions = c("A", "B", "C"), group_by = "pathway", save = FALSE)
  tab <- attr(p, "table")
  expect_setequal(tab$condition[tab$compound_id == "c1"], c("A", "B"))    # c1: in A and B, not C
  expect_setequal(tab$condition[tab$compound_id == "c3"], c("B", "C"))    # c3: in B and C, not A
  expect_equal(levels(tab$condition), c("A", "B", "C"))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_condition_flow() pools rare groups into 'Other' beyond top_n_groups", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .cflow_fixture()
  p <- plot_condition_flow(proj, conditions = c("A", "B", "C"), group_by = "class", top_n_groups = 1, save = FALSE)
  tab <- attr(p, "table")
  expect_true("Other" %in% tab$group)
  expect_equal(length(unique(tab$group)), 2L) # the most frequent class + Other
})

test_that("plot_condition_flow() validates and errors clearly on missing prerequisites", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .cflow_fixture()
  expect_error(plot_condition_flow(proj, conditions = "A", save = FALSE), "at least 2")
  expect_error(plot_condition_flow(proj, conditions = c("A", "nope"), save = FALSE), "not found")

  proj2 <- .test_project()
  binarizedMatrix(proj2) <- data.frame(compound_id = "c1", A = 1, B = 0, stringsAsFactors = FALSE)
  expect_error(plot_condition_flow(proj2, conditions = c("A", "B"), save = FALSE), "compounds_classified")
})

test_that("plot_condition_trend() aligns groups on a shared x via condition_label, and computes the secondary axis", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_condition_compare") <- data.frame(
    condition = c("A1", "A2", "B1", "B2"), disease_id = "D1",
    n_compounds_present = c(10, 12, 8, 14), n_compounds_scored = c(10, 12, 8, 14), n_compounds_finite = c(10, 12, 8, 14),
    mean_z = c(-3, -2, -1, 0), median_z = c(-3, -2, -1, 0), sd_z = 1, best_z = c(-3, -2, -1, 0),
    n_disease_genes_hit = 5, n_disease_genes_total = 10, disease_gene_coverage = 0.5,
    stringsAsFactors = FALSE
  )
  grp <- c(A1 = "G1", B1 = "G1", A2 = "G2", B2 = "G2")
  lbl <- c(A1 = "Stage1", B1 = "Stage1", A2 = "Stage2", B2 = "Stage2")

  p <- plot_condition_trend(proj, conditions = c("A1", "B1", "A2", "B2"), group = grp, condition_label = lbl, save = FALSE)
  tab <- attr(p, "table")
  expect_equal(nrow(tab), 4L)
  expect_no_error(ggplot2::ggplot_build(p))

  expect_error(plot_condition_trend(proj, conditions = "A1", save = FALSE), "at least 2")
  expect_error(plot_condition_trend(proj, conditions = c("A1", "B1"), group = c(A1 = "G1"), save = FALSE), "group")
  expect_warning(plot_condition_trend(proj, conditions = c("A1", "nope"), save = FALSE), "Only 1 of 2")
  expect_error(plot_condition_trend(proj, conditions = c("nope1", "nope2"), save = FALSE), "No rows")
})

test_that("plot_condition_trend() warns (not errors) when some requested (condition, disease) rows are missing", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .test_project()
  patliRResults(proj, "network_condition_compare") <- data.frame(
    condition = c("A1", "A2"), disease_id = "D1",
    n_compounds_present = 10, n_compounds_scored = 10, n_compounds_finite = 10,
    mean_z = -1, median_z = -1, sd_z = 1, best_z = -1,
    n_disease_genes_hit = 5, n_disease_genes_total = 10, disease_gene_coverage = 0.5,
    stringsAsFactors = FALSE
  )
  expect_warning(plot_condition_trend(proj, conditions = c("A1", "A2", "A3"), save = FALSE), "Only 2 of 3")
})
