test_that("plot_bowtie() returns a ggplot and errors clearly without network_bowtie()", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")

  proj <- .network_stats_test_setup()
  expect_error(plot_bowtie(proj, condition = "FLO-ET", save = FALSE), "network_bowtie")

  ## Fabricate network_bowtie results directly -- bypassing the real
  ## STRINGdb-dependent computation, same offline pattern as
  ## test-network_synergy.R.
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  set.seed(1)
  fake_bowtie <- data.frame(
    condition = "FLO-ET", compound_id = ct$compound_id, uniprot_id = ct$uniprot_id,
    string_id = paste0("9606.", ct$uniprot_id),
    bowtie_component = sample(c("core", "in_component", "out_component", "other"), nrow(ct), replace = TRUE),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_bowtie") <- fake_bowtie

  p <- plot_bowtie(proj, condition = "FLO-ET", save = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_bowtie() gives every network_bowtie() component level -- including not_in_action_network -- a named colour", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")

  proj <- .network_stats_test_setup()
  ct <- unique(patliRResults(proj, "network_edges")[patliRResults(proj, "network_edges")$condition == "FLO-ET", c("compound_id", "uniprot_id")])
  all_levels <- c("core", "in_component", "out_component", "not_in_action_network", "other", "unmapped")
  fake_bowtie <- data.frame(
    condition = "FLO-ET", compound_id = ct$compound_id, uniprot_id = ct$uniprot_id,
    string_id = paste0("9606.", ct$uniprot_id),
    bowtie_component = rep_len(all_levels, nrow(ct)),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_bowtie") <- fake_bowtie

  p <- plot_bowtie(proj, condition = "FLO-ET", save = FALSE)
  ## before the fix the sixth level had no entry in component_colors, so
  ## ggplot_build() warned ("Removed ... unknown fill") and emitted an NA
  ## fill for what is often the largest stratum -- assert neither happens.
  built <- expect_no_warning(ggplot2::ggplot_build(p))
  expect_false(anyNA(built$data[[1]]$fill))
})

test_that("plot_bowtie() requires ggalluvial with a clear message", {
  testthat::skip_if(requireNamespace("ggalluvial", quietly = TRUE), "ggalluvial is installed -- nothing to test here")
  testthat::skip_if_not_installed("ggplot2")

  proj <- .network_stats_test_setup()
  expect_error(plot_bowtie(proj, condition = "FLO-ET", save = FALSE), "ggalluvial")
})
