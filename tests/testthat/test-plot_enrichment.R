.plot_enrichment_fake <- function(condition = "FLO-ET", n_go = 5, n_reactome = 5) {
  mk <- function(db, n, offset) {
    if (n == 0) {
      return(data.frame(
        condition = character(0), db = character(0), ID = character(0),
        Description = character(0), GeneRatio = character(0), BgRatio = character(0),
        pvalue = double(0), p.adjust = double(0), geneID = character(0), Count = integer(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      condition = condition, db = db,
      ID = paste0(db, "_T", seq_len(n)),
      Description = paste0(db, " term ", seq_len(n)),
      GeneRatio = paste0(seq_len(n) + offset, "/20"),
      BgRatio = "100/2000",
      pvalue = seq(0.001, 0.04, length.out = n),
      p.adjust = seq(0.001, 0.04, length.out = n),
      geneID = "G1/G2",
      Count = seq_len(n) + offset,
      stringsAsFactors = FALSE
    )
  }
  rbind(mk("go", n_go, 0), mk("reactome", n_reactome, 2))
}

test_that("plot_enrichment() errors clearly without network_enrichment()", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  expect_error(plot_enrichment(proj, condition = "FLO-ET", save = FALSE), "network_enrichment")
})

test_that("plot_enrichment() draws a bubble plot, facets by db, and parses GeneRatio correctly", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  fake <- .plot_enrichment_fake()
  ## qvalue is entirely absent from this table -- clusterProfiler drops it
  ## when qvalue's pi0 estimation fails; this plot must not require it.
  expect_false("qvalue" %in% names(fake))
  patliRResults(proj, "network_enrichment") <- fake

  p <- expect_no_error(plot_enrichment(proj, condition = "FLO-ET", save = FALSE))
  expect_s3_class(p, "ggplot")

  ## faceting happened: two distinct db values -> two panels
  built <- ggplot2::ggplot_build(p)
  expect_gt(length(unique(built$data[[1]]$PANEL)), 1)

  ## GeneRatio "3/20" -> 0.15
  expect_equal(patliR:::.enrichment_parse_ratio("3/20"), 0.15)
  expect_equal(patliR:::.enrichment_parse_ratio(c("1/10", "5/20", NA, "bogus")),
               c(0.1, 0.25, NA_real_, NA_real_))

  expect_true(all(p$data$gene_ratio == patliR:::.enrichment_parse_ratio(p$data$GeneRatio)))
})

test_that("plot_enrichment() tolerates a table missing qvalue entirely and rows with malformed GeneRatio", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  fake <- .plot_enrichment_fake(n_go = 3, n_reactome = 0)
  fake$GeneRatio[1] <- "not-a-ratio"
  patliRResults(proj, "network_enrichment") <- fake

  expect_warning(
    p <- plot_enrichment(proj, condition = "FLO-ET", save = FALSE),
    "GeneRatio"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2L) ## the malformed row was dropped
})

test_that("plot_enrichment() restricts to the requested db and errors on an unknown one", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_enrichment") <- .plot_enrichment_fake()

  p <- plot_enrichment(proj, condition = "FLO-ET", db = "go", save = FALSE)
  expect_true(all(p$data$db == "go"))

  expect_error(
    plot_enrichment(proj, condition = "FLO-ET", db = "not_a_db", save = FALSE),
    "network_enrichment"
  )
})

test_that("plot_enrichment() top_n restricts to the most-significant terms per (condition, db)", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_enrichment") <- .plot_enrichment_fake(n_go = 10, n_reactome = 0)

  p <- plot_enrichment(proj, condition = "FLO-ET", db = "go", top_n = 3, save = FALSE)
  expect_equal(nrow(p$data), 3L)
  expect_equal(sort(p$data$p.adjust), sort(head(sort(.plot_enrichment_fake(n_go = 10, n_reactome = 0)$p.adjust), 3)))
})

test_that("plot_enrichment() saves a PNG and logs it, keyed by (condition, db)", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_enrichment") <- .plot_enrichment_fake()

  p <- plot_enrichment(proj, condition = "FLO-ET", save = TRUE)
  proj2 <- attr(p, "proj")
  log_df <- patliRResults(proj2, "enrichment_plot_log")
  expect_true(file.exists(log_df$path[1]))
  expect_true(all(c("condition", "db", "path") %in% names(log_df)))
})

test_that("plot_enrichment() returns a girafe widget with engine = 'ggiraph'", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggiraph")
  proj <- .network_stats_test_setup()
  patliRResults(proj, "network_enrichment") <- .plot_enrichment_fake()
  p <- plot_enrichment(proj, condition = "FLO-ET", save = FALSE, engine = "ggiraph")
  expect_s3_class(p, "girafe")
})
