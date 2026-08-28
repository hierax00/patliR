test_that("network_filter_proteome() keeps only edges whose target is in the given proteome", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  keep_target <- edges$uniprot_id[1]

  proj <- network_filter_proteome(proj, proteome = keep_target)
  filtered <- patliRResults(proj, "network_filtered_edges")

  expect_true(all(c("condition", "compound_id", "uniprot_id", "weight") %in% names(filtered)))
  expect_true(all(filtered$uniprot_id == keep_target))
  expect_true(nrow(filtered) == sum(edges$uniprot_id == keep_target))
})

test_that("network_filter_proteome() returns zero rows (not an error) when nothing matches", {
  proj <- .network_stats_test_setup()
  proj <- network_filter_proteome(proj, proteome = "NOT-A-REAL-UNIPROT-ID")
  filtered <- patliRResults(proj, "network_filtered_edges")
  expect_equal(nrow(filtered), 0)
})

test_that("network_filter_proteome() zero-row rerun drops the recomputed condition's stale rows, keeps the others", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  conds <- unique(edges$condition)
  skip_if(length(conds) < 2, "needs >= 2 built conditions")
  cA <- conds[1]
  cB <- conds[2]
  tgtA <- edges$uniprot_id[edges$condition == cA][1]
  tgtB <- edges$uniprot_id[edges$condition == cB][1]

  proj <- network_filter_proteome(proj, proteome = c(tgtA, tgtB))
  f1 <- patliRResults(proj, "network_filtered_edges")
  expect_true(cA %in% f1$condition)
  expect_true(cB %in% f1$condition)

  ## rerun for condition cA only against a proteome that matches nothing
  proj <- network_filter_proteome(proj, proteome = "NOT-A-REAL-UNIPROT-ID", condition = cA)
  f2 <- patliRResults(proj, "network_filtered_edges")
  expect_false(cA %in% f2$condition)   # stale cA rows gone (was the bug)
  expect_true(cB %in% f2$condition)    # cB untouched
})

test_that("network_filter_proteome() never modifies network_edges itself", {
  proj <- .network_stats_test_setup()
  before <- patliRResults(proj, "network_edges")
  proj <- network_filter_proteome(proj, proteome = before$uniprot_id[1])
  expect_identical(patliRResults(proj, "network_edges"), before)
})
