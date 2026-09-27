.ccmp_fixture <- function() {
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = c("A", "A", "A", "A", "B", "B"),
    compound_id = c("c1", "c1", "c2", "c2", "c3", "c3"),
    uniprot_id = c("t1", "t2", "t1", "t3", "t1", "t4"),
    weight = c(0.9, 0.8, 0.9, 0.7, 0.9, 0.6),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "disease one",
    uniprot_id = c("t1", "t2"), gene_symbol = c("T1", "T2"),
    association_score = c(0.9, 0.5), source = "open_targets", fetched_at = Sys.time(),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- data.frame(
    condition = c("A", "A", "B"), compound_id = c("c1", "c2", "c3"), disease_id = "D1",
    disease_gene_source = "disease_genes", n_targets_mapped = 1, n_disease_genes_mapped = 2,
    n_overlap = 1, d_observed = 1, d_random_mean = 2, d_random_sd = 1,
    z_score = c(-2, -1, 0.5), p_empirical = 0.1, n_random = 100, seed_used = 1,
    species = 9606, string_version = "12.0", score_threshold = 400,
    n_tests_in_family = NA_integer_, p_adjusted = NA_real_, stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_centrality") <- data.frame(
    condition = c("A", "A", "B"), node_id = c("t1", "t2", "t1"), node_type = "target",
    degree = c(2, 1, 1), hub_score = c(0.8, 0.4, 0.5), stringsAsFactors = FALSE
  )
  proj
}

test_that("network_condition_compare() aggregates per-compound z to one row per condition, sorted most-negative first", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  cmp <- patliRResults(proj, "network_condition_compare")

  expect_equal(cmp$condition, c("A", "B")) # A (mean -1.5) before B (mean 0.5)
  expect_equal(cmp$mean_z, c(-1.5, 0.5))
  expect_equal(cmp$median_z, c(-1.5, 0.5))
  expect_equal(cmp$n_compounds_present, c(2L, 1L))
  expect_equal(cmp$n_disease_genes_total, c(2L, 2L))
  expect_equal(cmp$n_disease_genes_hit, c(2L, 1L)) # A: t1+t2 both hit; B: only t1
  expect_equal(cmp$disease_gene_coverage, c(1, 0.5))
  expect_true(file.exists(file.path(projectDir(proj), "results", "network_condition_compare.csv")))
})

test_that("network_condition_compare() ranks the disease-gene targets driving each condition's signal", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  tg <- patliRResults(proj, "network_condition_targets")

  a <- tg[tg$condition == "A", ]
  expect_equal(a$uniprot_id, c("t1", "t2")) # t1 (2 compounds) ranked above t2 (1 compound)
  expect_equal(a$n_compounds, c(2L, 1L))
  expect_false("t3" %in% a$uniprot_id) # not a disease gene -> excluded

  b <- tg[tg$condition == "B", ]
  expect_equal(b$uniprot_id, "t1")
  expect_false("t4" %in% b$uniprot_id) # not a disease gene
  expect_equal(b$degree, 1)
})

test_that("network_condition_compare() aborts when conditions have inconsistent proximity provenance", {
  proj <- .ccmp_fixture()
  prox <- patliRResults(proj, "network_proximity")
  prox$score_threshold[prox$condition == "B"] <- 700 # B was computed at a different STRING threshold than A
  patliRResults(proj, "network_proximity") <- prox
  expect_error(network_condition_compare(proj, disease = "D1"), "score_threshold")
})

test_that("network_condition_compare() re-run wipes a target that no longer qualifies (lost its only edge)", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  expect_true("t2" %in% patliRResults(proj, "network_condition_targets")$uniprot_id[patliRResults(proj, "network_condition_targets")$condition == "A"])

  ## c1 (the only compound hitting t2) loses that edge
  ne <- patliRResults(proj, "network_edges")
  patliRResults(proj, "network_edges") <- ne[!(ne$condition == "A" & ne$uniprot_id == "t2"), ]
  proj <- network_condition_compare(proj, disease = "D1")
  tg <- patliRResults(proj, "network_condition_targets")
  expect_false("t2" %in% tg$uniprot_id[tg$condition == "A"]) # not a stale leftover row
})

test_that("network_condition_compare() aborts when a later call's proximity provenance differs from a condition already stored for the same disease (regression)", {
  ## The pre-existing check above only guards conditions passed to the SAME
  ## call -- a single-condition call is trivially self-consistent, so a
  ## second call comparing just "B" (recomputed at a different STRING
  ## threshold) after "A" was already stored used to pass its own check and
  ## get stored right alongside "A", silently combinable by
  ## plot_condition_compare()/plot_condition_trend() reading the whole table.
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, conditions = "A", disease = "D1")
  prox <- patliRResults(proj, "network_proximity")
  prox$score_threshold[prox$condition == "B"] <- 700
  patliRResults(proj, "network_proximity") <- prox
  expect_error(network_condition_compare(proj, conditions = "B", disease = "D1"), "different")
})

test_that("network_condition_compare() lets a later call replace the SAME condition under new provenance (not a false positive)", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, conditions = "A", disease = "D1")
  prox <- patliRResults(proj, "network_proximity")
  prox$score_threshold[prox$condition == "A"] <- 700
  patliRResults(proj, "network_proximity") <- prox
  ## re-running "A" itself under new provenance is a legitimate change of
  ## mind about A, not a cross-condition mismatch
  expect_no_error(network_condition_compare(proj, conditions = "A", disease = "D1"))
})

test_that("network_condition_compare() stamps proximity provenance onto every output row", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  cmp <- patliRResults(proj, "network_condition_compare")
  expect_true(all(c("disease_gene_source", "score_threshold", "species", "string_version") %in% names(cmp)))
  expect_true(all(cmp$score_threshold == 400))
})

test_that(".condition_compare_disease_genes() collapses a duplicated accession to its highest association_score", {
  proj <- .ccmp_fixture()
  dg <- patliRResults(proj, "disease_genes")
  dg <- rbind(dg, data.frame(disease_id = "D1", disease_name = "disease one", uniprot_id = "t1",
                             gene_symbol = "T1", association_score = 0.3, source = "open_targets",
                             fetched_at = Sys.time(), stringsAsFactors = FALSE))
  patliRResults(proj, "disease_genes") <- dg
  out <- patliR:::.condition_compare_disease_genes(proj, "D1")
  expect_equal(nrow(out), 2L) # not 3 -- the duplicate t1 row collapsed
  expect_equal(out$association_score[out$uniprot_id == "t1"], 0.9) # kept the higher score, not the last row
})

test_that("network_condition_compare() honours a conditions= subset and aggregate = 'median'", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, conditions = "A", disease = "D1", aggregate = "median")
  cmp <- patliRResults(proj, "network_condition_compare")
  expect_equal(cmp$condition, "A")
})

test_that("network_condition_compare() validates its arguments and prerequisites", {
  proj <- .test_project()
  expect_error(network_condition_compare(proj, disease = "D1"), "network_edges")

  proj2 <- .ccmp_fixture()
  patliRResults(proj2, "disease_genes") <- NULL
  expect_error(network_condition_compare(proj2, disease = "D1"), "No rows for")
})

test_that("plot_condition_compare() draws the stored table, sorted, with a required disease", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  p <- plot_condition_compare(proj, disease = "D1", save = FALSE)
  tab <- attr(p, "table")
  expect_equal(as.character(tab$condition), c("A", "B"))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(plot_condition_compare(proj, disease = "nope"), "No rows")
})

test_that("plot_condition_compare() draws every disease side by side and ranks conditions by their mean z", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  ## a second disease where B (not A) is the stronger condition
  patliRResults(proj, "disease_genes") <- rbind(
    patliRResults(proj, "disease_genes"),
    data.frame(disease_id = "D2", disease_name = "disease two", uniprot_id = c("t1", "t4"),
              gene_symbol = c("T1", "T4"), association_score = c(0.5, 0.5), source = "open_targets",
              fetched_at = Sys.time(), stringsAsFactors = FALSE)
  )
  patliRResults(proj, "network_proximity") <- rbind(
    patliRResults(proj, "network_proximity"),
    data.frame(condition = c("A", "B"), compound_id = c("c1", "c3"), disease_id = "D2",
              disease_gene_source = "disease_genes", n_targets_mapped = 1, n_disease_genes_mapped = 2,
              n_overlap = 1, d_observed = 1, d_random_mean = 2, d_random_sd = 1,
              z_score = c(0.2, -3), p_empirical = 0.1, n_random = 100, seed_used = 1,
              species = 9606, string_version = "12.0", score_threshold = 400,
              n_tests_in_family = NA_integer_, p_adjusted = NA_real_, stringsAsFactors = FALSE)
  )
  proj <- network_condition_compare(proj, disease = "D2")

  p <- plot_condition_compare(proj, save = FALSE) # disease = NULL -> both
  tab <- attr(p, "table")
  expect_setequal(tab$disease_id, c("D1", "D2"))
  ## overall mean across D1 (-1.5)/D2 (0.2) for A = -0.65; D1 (0.5)/D2 (-3) for B = -1.25 -> B ranks first
  expect_equal(levels(tab$condition)[nlevels(tab$condition)], "B")
  expect_no_error(ggplot2::ggplot_build(p))
})

.ccmp_two_disease_fixture <- function() {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  patliRResults(proj, "disease_genes") <- rbind(
    patliRResults(proj, "disease_genes"),
    data.frame(disease_id = "D2", disease_name = "disease two", uniprot_id = c("t1", "t4"),
              gene_symbol = c("T1", "T4"), association_score = c(0.5, 0.5), source = "open_targets",
              fetched_at = Sys.time(), stringsAsFactors = FALSE)
  )
  patliRResults(proj, "network_proximity") <- rbind(
    patliRResults(proj, "network_proximity"),
    data.frame(condition = c("A", "B"), compound_id = c("c1", "c3"), disease_id = "D2",
              disease_gene_source = "disease_genes", n_targets_mapped = 1, n_disease_genes_mapped = 2,
              n_overlap = 1, d_observed = 1, d_random_mean = 2, d_random_sd = 1,
              z_score = c(0.2, -3), p_empirical = 0.1, n_random = 100, seed_used = 1,
              species = 9606, string_version = "12.0", score_threshold = 400,
              n_tests_in_family = NA_integer_, p_adjusted = NA_real_, stringsAsFactors = FALSE)
  )
  network_condition_compare(proj, disease = "D2")
}

test_that("plot_condition_disease_flow() draws every condition and disease, weight = -z (floored), best marked", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .ccmp_two_disease_fixture()

  p <- plot_condition_disease_flow(proj, save = FALSE)
  tab <- attr(p, "table")
  expect_equal(nrow(tab), 4L) # 2 conditions x 2 diseases, none dropped ("con todas")
  ## A: D1 z=-1.5 (best, weight 1.5), D2 z=0.2 (weight floored at 0.05, not negative/zero)
  expect_equal(tab$weight[tab$condition == "A" & tab$disease_id == "D1"], 1.5)
  expect_equal(tab$weight[tab$condition == "A" & tab$disease_id == "D2"], 0.05)
  expect_true(tab$best_for_condition[tab$condition == "A" & tab$disease_id == "D1"])
  expect_false(tab$best_for_condition[tab$condition == "A" & tab$disease_id == "D2"])
  ## B: D2 (z=-3) is its best match, not D1 (z=0.5)
  expect_true(tab$best_for_condition[tab$condition == "B" & tab$disease_id == "D2"])
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("network_condition_compare(permutation_test = TRUE) is reproducible with a fixed seed and BH-adjusts across conditions", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1", permutation_test = TRUE, n_perm = 500, seed = 1)
  cmp <- patliRResults(proj, "network_condition_compare")

  expect_true(all(c("perm_null_mean", "perm_null_sd", "perm_p_value", "perm_p_adjusted",
                    "perm_n_draws", "perm_seed_used") %in% names(cmp)))
  expect_equal(cmp$perm_n_draws, c(500L, 500L))
  expect_equal(cmp$perm_seed_used, c(1L, 1L))
  expect_equal(cmp$perm_p_adjusted, stats::p.adjust(cmp$perm_p_value, "BH"))
  ## same seed -> identical draws -> identical p-values
  proj2 <- .ccmp_fixture()
  proj2 <- network_condition_compare(proj2, disease = "D1", permutation_test = TRUE, n_perm = 500, seed = 1)
  expect_equal(patliRResults(proj2, "network_condition_compare")$perm_p_value, cmp$perm_p_value)
})

test_that("network_condition_compare(permutation_test = TRUE) does not perturb the caller's RNG state", {
  proj <- .ccmp_fixture()
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  proj <- network_condition_compare(proj, disease = "D1", permutation_test = TRUE, n_perm = 200)
  after <- runif(1)
  expect_equal(before, after)
})

test_that("network_condition_compare(permutation_test = TRUE) handles a singleton pooled z correctly (regression)", {
  ## base R's sample(x, n) treats a length-1 NUMERIC x as "sample from 1:x"
  ## (?sample's documented surprise), not "draw n values from the vector x".
  ## With a single-compound pool (zc of length 1), the old
  ## `sample(zc, length(ids))` therefore drew from 1:zc[[1]] instead of
  ## returning zc[[1]] -- for a z of 5.5 this produced a fabricated null
  ## distribution (integers 1-5) that never reached 5.5, giving a spuriously
  ## small two-sided p-value (~2/(n_perm+1)) for what has only one possible
  ## subset and is therefore not a real permutation test at all (true p = 1).
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "A", compound_id = "c1", uniprot_id = "t1", weight = 0.9, stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "d", uniprot_id = "t1", gene_symbol = "T1",
    association_score = 0.9, source = "test", fetched_at = Sys.time(), stringsAsFactors = FALSE
  )
  patliRResults(proj, "network_proximity") <- data.frame(
    condition = "A", compound_id = "c1", disease_id = "D1",
    disease_gene_source = "disease_genes", n_targets_mapped = 1, n_disease_genes_mapped = 1,
    n_overlap = 1, d_observed = 1, d_random_mean = 2, d_random_sd = 1,
    z_score = 5.5, p_empirical = 0.1, n_random = 100, seed_used = 1,
    species = 9606, string_version = "12.0", score_threshold = 400,
    n_tests_in_family = NA_integer_, p_adjusted = NA_real_, stringsAsFactors = FALSE
  )
  proj <- network_condition_compare(proj, disease = "D1", permutation_test = TRUE, n_perm = 500, seed = 1)
  cmp <- patliRResults(proj, "network_condition_compare")
  expect_equal(cmp$perm_null_mean, 5.5)
  expect_equal(cmp$perm_null_sd, 0)
  expect_equal(cmp$perm_p_value, 1)
})

test_that("network_condition_compare() rerun with permutation_test = FALSE after TRUE does not abort (regression)", {
  ## .network_upsert() used to abort on column removal even when every
  ## existing row for the touched keys was about to be replaced -- so
  ## turning permutation_test off again after a TRUE run crashed instead of
  ## just dropping the perm_* columns for this rerun.
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1", permutation_test = TRUE, n_perm = 200)
  expect_true("perm_p_value" %in% names(patliRResults(proj, "network_condition_compare")))

  proj <- network_condition_compare(proj, disease = "D1", permutation_test = FALSE)
  cmp <- patliRResults(proj, "network_condition_compare")
  expect_false("perm_p_value" %in% names(cmp))
  expect_equal(cmp$condition, c("A", "B"))
})

test_that("network_condition_compare() rejects a non-positive n_perm", {
  proj <- .ccmp_fixture()
  expect_error(network_condition_compare(proj, disease = "D1", permutation_test = TRUE, n_perm = 0),
              "n_perm")
})

test_that("network_condition_compare(permutation_test = TRUE) leaves perm_* absent (not error) when default", {
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  cmp <- patliRResults(proj, "network_condition_compare")
  expect_false(any(grepl("^perm_", names(cmp))))
})

test_that("plot_condition_disease_flow() requires at least two diseases", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .ccmp_fixture()
  proj <- network_condition_compare(proj, disease = "D1")
  expect_error(plot_condition_disease_flow(proj, save = FALSE), "at least two diseases")
})

test_that("plot_condition_compare()/plot_condition_disease_flow() keep two diseases with the same display label distinct (regression)", {
  ## Grouping/colouring/faceting by disease_label (display text) instead of
  ## disease_id used to fuse two different diseases that happen to share a
  ## name -- one dodge position / colour in plot_condition_compare(), one
  ## merged alluvial stratum in plot_condition_disease_flow().
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggalluvial")
  proj <- .ccmp_two_disease_fixture()
  dg <- patliRResults(proj, "disease_genes")
  dg$disease_name[dg$disease_id %in% c("D1", "D2")] <- "same label"
  patliRResults(proj, "disease_genes") <- dg

  p <- plot_condition_compare(proj, save = FALSE)
  tab <- attr(p, "table")
  expect_setequal(tab$disease_id, c("D1", "D2"))
  expect_equal(unique(tab$disease_label), "same label") # the label really did collide
  built <- ggplot2::ggplot_build(p)
  point_layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  expect_length(unique(built$data[[point_layer]]$colour), 2L) # still two distinct colours

  pf <- plot_condition_disease_flow(proj, save = FALSE)
  tabf <- attr(pf, "table")
  expect_setequal(tabf$disease_id, c("D1", "D2"))
  builtf <- ggplot2::ggplot_build(pf)
  alluvium_layer <- which(vapply(pf$layers, function(l) inherits(l$geom, "GeomAlluvium"), logical(1)))
  expect_length(unique(builtf$data[[alluvium_layer]]$fill), 2L) # still two distinct strata
})
