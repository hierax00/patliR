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
  ## "NOT-A-REAL-UNIPROT-ID" matches nothing AND does not look like a
  ## UniProt accession, so this also exercises the new sanity-check warning
  ## (spec 1.11 [SHOULD] #2) -- expected here, asserted more directly below.
  expect_warning(
    proj <- network_filter_proteome(proj, proteome = "NOT-A-REAL-UNIPROT-ID"),
    "UniProt"
  )
  filtered <- patliRResults(proj, "network_filtered_edges")
  expect_equal(nrow(filtered), 0)
})

test_that("network_filter_proteome() warns when the input looks like gene symbols, not UniProt accessions", {
  proj <- .network_stats_test_setup()
  ## Neither "TP53" nor "EGFR" matches the 6-character UniProt accession
  ## shape, and neither will be present in the fixture's uniprot_id column.
  expect_warning(
    network_filter_proteome(proj, proteome = c("TP53", "EGFR")),
    "gene symbols"
  )
})

test_that("network_filter_proteome() does not warn when the input is UniProt-shaped even if it matches nothing", {
  proj <- .network_stats_test_setup()
  ## "P99999" is UniProt-accession-shaped (matches the regex) even though it
  ## is not present in this fixture -- more than half of a length-1 vector
  ## being accession-shaped means no warning should fire.
  expect_no_warning(network_filter_proteome(proj, proteome = "P99999"))
})

test_that("network_filter_proteome() records proteome_label and n_proteome, and keys the upsert on (condition, proteome_label)", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  tgt <- edges$uniprot_id[1]

  proj <- network_filter_proteome(proj, proteome = tgt, proteome_label = "my_proteome")
  f1 <- patliRResults(proj, "network_filtered_edges")
  expect_true(all(c("proteome_label", "n_proteome") %in% names(f1)))
  expect_true(all(f1$proteome_label == "my_proteome"))
  expect_true(all(f1$n_proteome == 1L))

  ## a second, differently-labelled filter over the SAME condition must be
  ## kept alongside the first, not overwrite it -- proteome_label is part
  ## of the upsert key.
  proj <- network_filter_proteome(proj, proteome = tgt, proteome_label = "other_proteome")
  f2 <- patliRResults(proj, "network_filtered_edges")
  expect_setequal(unique(f2$proteome_label), c("my_proteome", "other_proteome"))
  expect_equal(sum(f2$proteome_label == "my_proteome"), nrow(f1))
})

test_that("network_filter_proteome() auto-derives a stable proteome_label when none is given", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  tgt <- edges$uniprot_id[1]

  proj1 <- network_filter_proteome(proj, proteome = tgt)
  proj2 <- network_filter_proteome(proj, proteome = tgt)
  lab1 <- unique(patliRResults(proj1, "network_filtered_edges")$proteome_label)
  lab2 <- unique(patliRResults(proj2, "network_filtered_edges")$proteome_label)
  expect_equal(lab1, lab2) ## same proteome content -> same auto label
  expect_false(is.na(lab1) || !nzchar(lab1))
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

  ## Same proteome_label both calls: this models re-running the SAME named
  ## filter with an updated proteome list, which is what headline bug #2
  ## (a zero-row rerun silently keeping the previous run's rows) is about.
  ## A DIFFERENT proteome_label would legitimately coexist instead (see the
  ## "keys the upsert on (condition, proteome_label)" test above) -- that is
  ## the new, intended behaviour, not the bug this test guards against.
  proj <- network_filter_proteome(proj, proteome = c(tgtA, tgtB), proteome_label = "shared_label")
  f1 <- patliRResults(proj, "network_filtered_edges")
  expect_true(cA %in% f1$condition)
  expect_true(cB %in% f1$condition)

  ## rerun for condition cA only, same label, against a proteome that now
  ## matches nothing
  expect_warning(
    proj <- network_filter_proteome(
      proj, proteome = "NOT-A-REAL-UNIPROT-ID", condition = cA, proteome_label = "shared_label"
    ),
    "UniProt"
  )
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

test_that("network_filter_proteome() keeps distinct proteomes with colliding weighted character sums", {
  proj <- .network_stats_test_setup()
  edges <- patliRResults(proj, "network_edges")
  cond <- edges$condition[1]
  edges <- edges[rep(1L, 2L), , drop = FALSE]
  edges$uniprot_id <- c("P10040", "P10500")
  patliRResults(proj, "network_edges") <- edges

  proj <- network_filter_proteome(proj, "P10040", condition = cond)
  proj <- network_filter_proteome(proj, "P10500", condition = cond)
  filtered <- patliRResults(proj, "network_filtered_edges")
  expect_setequal(filtered$uniprot_id, c("P10040", "P10500"))
  expect_equal(length(unique(filtered$proteome_label)), 2L)
  expect_identical(
    patliR:::.network_proteome_auto_label(c("P10040", "P10500", "P10040")),
    patliR:::.network_proteome_auto_label(c("P10500", "P10040"))
  )
})
