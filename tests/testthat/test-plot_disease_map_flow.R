## Same offline fixture as test-plot_disease_map.R's .dmap_project()/.dmap_ids
## (testthat edition 3 isolates each test file's environment, so it is not
## visible across files and is duplicated here rather than shared): network_edges
## conditions "A" (C1/C2/C3, targets T1-T4) and "B" (C1, target T5 only),
## targets_disease_profile/disease_genes/targets_disease as the association sources.

.dmf_ids <- c(T1 = "P35354", T2 = "P23219", T3 = "P08253", T4 = "P00533", T5 = "P04637")

.dmf_project <- function() {
  id <- .dmf_ids
  proj <- patliR_project(tempfile("patliR_dmf_"))
  patliRResults(proj, "network_edges") <- data.frame(
    condition = c("A", "A", "A", "A", "A", "A", "B"),
    compound_id = c("C1", "C1", "C1", "C2", "C2", "C3", "C1"),
    uniprot_id = unname(id[c("T1", "T2", "T3", "T1", "T2", "T4", "T5")]),
    weight = c(0.9, 0.8, 0.7, 0.6, 0.95, 0.9, 0.9),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease_profile") <- data.frame(
    compound_id = "C1",
    target_id = unname(id[c("T1", "T1", "T2", "T2", "T3", "T3", "T5")]),
    disease_id = c("D1", "D2", "D1", "D3", "D4", "D2", "D1"),
    disease_name = c("breast cancer", "Alzheimer disease", "breast cancer", "asthma",
                     "smoking initiation", "Alzheimer disease", "breast cancer"),
    association_score = c(0.8, 0.5, 0.6, 0.4, 0.05, 0.7, 0.9),
    evidence = NA_character_, rank = 1L, mode = "explore", stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = c("F1", "F1", "GO_0006954"),
    disease_name = c("hypertensive disorder", "hypertensive disorder", "inflammatory response (GO:0006954)"),
    uniprot_id = unname(id[c("T1", "T3", "T2")]), gene_symbol = NA_character_,
    association_score = c(0.7, 0.5, NA), source = "test", fetched_at = NA_character_,
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease") <- data.frame(
    compound_id = "C1", target_id = unname(id[c("T1", "T2")]), disease_id = "F1",
    association_score = c(0.2, 0.01), evidence = NA_character_, stringsAsFactors = FALSE
  )
  proj
}

.dmf_plot <- function(...) suppressWarnings(suppressMessages(plot_disease_map_flow(...)))

test_that("plot_disease_map_flow() assigns each target its single most-associated area and drops targets with none", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  id <- .dmf_ids
  p <- .dmf_plot(proj, conditions = c("A", "B"), min_score = 0.1, save = FALSE)
  tab <- attr(p, "table")

  ## T4 (P00533) has no disease association in any source -> silently excluded, not an error
  expect_false(id[["T4"]] %in% tab$target_id)

  a <- tab[tab$condition == "A", ]
  expect_setequal(a$target_id, id[c("T1", "T2", "T3")])
  ## T1: Neoplasm (D1, 0.8) beats Nervous (D2, 0.5) and Cardiovascular (F1, 0.7) -- all n_diseases = 1, highest score wins
  expect_equal(a$disease_class[a$target_id == id[["T1"]]], "Neoplasm")
  ## T2: Neoplasm (D1, 0.6) beats Respiratory (D3, 0.4) and the unscored GO row -- same tie-break
  expect_equal(a$disease_class[a$target_id == id[["T2"]]], "Neoplasm")
  ## T3: D4 (smoking initiation, 0.05) dropped by min_score = 0.1, leaving only D2 (Alzheimer) -> Nervous system
  expect_equal(a$disease_class[a$target_id == id[["T3"]]], "Nervous system & psychiatric")

  b <- tab[tab$condition == "B", ]
  expect_equal(b$target_id, id[["T5"]])
  expect_equal(b$disease_class, "Neoplasm")

  expect_equal(levels(factor(tab$condition, levels = c("A", "B"))), c("A", "B"))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_disease_map_flow() writes a PNG and logs it when save = TRUE (the default)", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  p <- suppressWarnings(suppressMessages(plot_disease_map_flow(proj, conditions = c("A", "B"), min_score = 0.1)))
  proj2 <- attr(p, "proj")
  expect_true(file.exists(file.path(projectDir(proj2), "results", "disease_map_flow_plot_log.csv")))
  log <- patliRResults(proj2, "disease_map_flow_plot_log")
  expect_equal(log$conditions, "A>B")
  expect_true(file.exists(file.path(projectDir(proj2), "plots", "disease_map_flow_A-B.png")))
})

test_that("plot_disease_map_flow() lets the same target's ribbon persist across conditions where it is present in both", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  id <- .dmf_ids
  ne <- patliRResults(proj, "network_edges")
  ne <- rbind(ne, data.frame(condition = "B", compound_id = "C1", uniprot_id = id[["T1"]], weight = 0.9))
  patliRResults(proj, "network_edges") <- ne

  p <- .dmf_plot(proj, conditions = c("A", "B"), min_score = 0.1, save = FALSE)
  tab <- attr(p, "table")
  expect_true(id[["T1"]] %in% tab$target_id[tab$condition == "A"])
  expect_true(id[["T1"]] %in% tab$target_id[tab$condition == "B"])
  expect_equal(tab$disease_class[tab$target_id == id[["T1"]] & tab$condition == "A"],
              tab$disease_class[tab$target_id == id[["T1"]] & tab$condition == "B"])
})

test_that("plot_disease_map_flow() validates conditions", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  expect_error(.dmf_plot(proj, conditions = "A", save = FALSE), "at least 2 condition")
  expect_error(.dmf_plot(proj, conditions = c("A", "A"), save = FALSE), "duplicate")
  expect_error(.dmf_plot(proj, conditions = c("A", "NOPE"), save = FALSE), "not built")
})

test_that("plot_disease_map_flow() requires network_edges and validates min_score/sources", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- patliR_project(tempfile("patliR_dmf_"))
  expect_error(.dmf_plot(proj, conditions = c("A", "B"), save = FALSE), "network_edges")

  proj2 <- .dmf_project()
  expect_error(.dmf_plot(proj2, conditions = c("A", "B"), min_score = 1.5, save = FALSE), "min_score")
  expect_error(.dmf_plot(proj2, conditions = c("A", "B"), sources = "nope", save = FALSE), "sources")
})

test_that("plot_disease_map_flow() errors clearly when no target in scope has any disease association", {
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  patliRResults(proj, "targets_disease_profile") <- NULL
  patliRResults(proj, "disease_genes") <- NULL
  patliRResults(proj, "targets_disease") <- NULL
  expect_error(.dmf_plot(proj, conditions = c("A", "B"), save = FALSE), "No target-disease association")
})

test_that("plot_disease_map_flow() handles an association table with NO scores at all (regression)", {
  ## Formula aggregate() defaults to na.omit, dropping every row before
  ## grouping -- if EVERY association_score in scope is NA (a real case: a
  ## purely GO-sourced disease/target association, which carries no score),
  ## the tie-break aggregate() used to error outright ("no rows to
  ## aggregate") rather than reach its own all(is.na(x)) -Inf branch.
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .dmf_project()
  prof <- patliRResults(proj, "targets_disease_profile"); prof$association_score <- NA_real_
  patliRResults(proj, "targets_disease_profile") <- prof
  dg <- patliRResults(proj, "disease_genes"); dg$association_score <- NA_real_
  patliRResults(proj, "disease_genes") <- dg
  td <- patliRResults(proj, "targets_disease"); td$association_score <- NA_real_
  patliRResults(proj, "targets_disease") <- td

  p <- .dmf_plot(proj, conditions = c("A", "B"), save = FALSE)
  tab <- attr(p, "table")
  expect_gt(nrow(tab), 0)
})
