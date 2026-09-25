## plot_ppi_network() is exercised fully offline: the two STRING readers
## (.network_stringdb_aliases_map() and .ppi_string_links()) and the GO-based
## GPCR table (.ppi_gpcr_go()) are mocked, so no STRING download and no
## org.Hs.eg.db query happens here.
##
## Synthetic set-up (condition "A"):
##   targets P00001..P00008 hit by compounds c1..c3 (P00002 hit by all three)
##   STRING ids s1..s7 for P00001..P00007; P00008 is not in STRING
##   links >= 700: s1-s2 (900), s2-s3 (800), s1-s3 (750), s4-s5 (950)
##   link below 700: s1-s7 (650)  -> P00007 has no partner at 700
##   P00006 (s6) only links to a non-target (s9)
##   disease "D1" genes: P00001, P00004, Q00009 (s9, not a target), Q00010 (unmapped)
##   GPCR (mocked GO table): P00002

.ppi_test_project <- function() {
  proj <- .test_project()
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "A",
    compound_id = c("c1", "c1", "c1", "c2", "c2", "c2", "c3", "c3", "c3", "c3", "c1"),
    uniprot_id = c("P00001", "P00002", "P00003", "P00002", "P00004", "P00005",
                   "P00002", "P00006", "P00007", "P00008", "P00005"),
    weight = 0.5, stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "test disease",
    uniprot_id = c("P00001", "P00004", "Q00009", "Q00010"),
    gene_symbol = c("GENE1", "GENE4", "GENE9", "GENE10"),
    association_score = c(0.9, 0.8, 0.7, 0.6), stringsAsFactors = FALSE
  )
  proj
}

.ppi_fake_aliases <- function() {
  data.frame(
    string_protein_id = c(paste0("s", 1:7), "s9", paste0("s", 10:99)),
    alias = c(sprintf("P%05d", 1:7), "Q00009", paste0("X", 10:99)),
    stringsAsFactors = FALSE
  )
}

.ppi_fake_links <- function() {
  data.frame(
    protein1 = c("s1", "s2", "s1", "s4", "s1", "s6", "s1", "s4"),
    protein2 = c("s2", "s3", "s3", "s5", "s7", "s9", "s9", "s9"),
    combined_score = c(900L, 800L, 750L, 950L, 650L, 900L, 850L, 820L),
    stringsAsFactors = FALSE
  )
}

.ppi_fake_gpcr <- function() {
  data.frame(uniprot_id = "P00002", symbol = "GPCRX", evidence = "IDA", stringsAsFactors = FALSE)
}

## Mocks the three data sources in the calling test's scope.
.ppi_local_mocks <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(
    .network_stringdb_aliases_map = function(...) .ppi_fake_aliases(),
    .ppi_string_links = function(proj, species, version, score_threshold) {
      l <- .ppi_fake_links()
      l[l$combined_score >= score_threshold, , drop = FALSE]
    },
    .ppi_gpcr_go = function() .ppi_fake_gpcr(),
    .package = "patliR", .env = env
  )
}

test_that("plot_ppi_network() draws the targets with a STRING partner and flags disease genes", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()

  p <- plot_ppi_network(proj, condition = "A", disease = "D1", save = FALSE)
  expect_s3_class(p, "ggplot")
  d <- attr(p, "data")
  nodes <- d$nodes
  expect_setequal(nodes$uniprot_id, c("P00001", "P00002", "P00003", "P00004", "P00005"))
  expect_true(all(nodes$degree >= 1))
  expect_equal(sort(nodes$uniprot_id[nodes$is_disease]), c("P00001", "P00004"))
  ## disease genes are ranked first
  expect_equal(nodes$uniprot_id[1:2], c("P00001", "P00004"))
  ## the 650 link (s1-s7) is below the default threshold of 700
  expect_equal(nrow(d$ppi_edges), 4L)
  expect_false(any(d$ppi_edges$uniprot_b == "P00007"))

  s <- attr(p, "summary")
  expect_equal(s$n_targets, 8L)
  expect_equal(s$n_targets_mapped, 7L)
  expect_equal(s$n_no_partner, 3L) # P00006, P00007, P00008
  expect_equal(s$n_targets_drawn, 5L)
  expect_equal(s$n_disease_hit, 2L)
  expect_equal(s$n_disease_genes_mapped, 3L)
  expect_equal(s$n_universe, 98L)
  expect_equal(s$p_hypergeom, stats::phyper(1, 3, 98 - 3, 7, lower.tail = FALSE))
  ## compound view: every compound hitting a drawn target is placed
  expect_setequal(d$compounds$compound_id, c("c1", "c2", "c3"))
  expect_true(all(d$compound_edges$uniprot_id %in% nodes$uniprot_id))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("score_threshold changes which interactions (and targets) are drawn", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  p <- plot_ppi_network(proj, condition = "A", score_threshold = 600, highlight = "none", save = FALSE)
  expect_true("P00007" %in% attr(p, "data")$nodes$uniprot_id)
  expect_equal(attr(p, "summary")$n_no_partner, 2L)
  expect_error(plot_ppi_network(proj, condition = "A", score_threshold = 990, save = FALSE),
               "No STRING interaction")
})

test_that("GPCRs and their interactions are highlighted, and the membership table is returned", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()

  p <- plot_ppi_network(proj, condition = "A", disease = "D1", save = FALSE)
  nodes <- attr(p, "data")$nodes
  expect_equal(nodes$uniprot_id[nodes$is_gpcr], "P00002")
  e <- attr(p, "data")$ppi_edges
  expect_equal(sum(e$gpcr_edge), 2L) # P00001-P00002, P00002-P00003
  s <- attr(p, "summary")
  expect_equal(s$highlight, "GPCR")
  expect_equal(s$n_gpcr_edges_drawn, 2L)
  expect_equal(s$n_gpcr_targets, 1L)
  g <- attr(p, "gpcr")
  expect_equal(nrow(g), 8L)
  expect_equal(g$uniprot_id[g$is_highlighted], "P00002")
  expect_true(g$drawn[g$uniprot_id == "P00002"])
  ## GPCR edges are drawn in their own (thicker) segment layer
  seg <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  widths <- vapply(seg, function(l) as.numeric(l$aes_params$linewidth %||% NA), numeric(1))
  expect_true(any(widths > 1))
  ## user set on top of / instead of the GO classification
  p2 <- plot_ppi_network(proj, condition = "A", highlight = "none", highlight_set = "P00005", save = FALSE)
  expect_equal(attr(p2, "data")$nodes$uniprot_id[attr(p2, "data")$nodes$is_gpcr], "P00005")
  expect_equal(attr(p2, "summary")$highlight, "Highlighted")
  p3 <- plot_ppi_network(proj, condition = "A", highlight = "none", save = FALSE)
  expect_equal(attr(p3, "summary")$n_gpcr_edges_drawn, 0L)
  expect_equal(attr(p3, "summary")$highlight, "none")
})

test_that(".ppi_highlight_table() matches by accession, then by symbol, plus a user set", {
  testthat::local_mocked_bindings(.ppi_gpcr_go = function() {
    data.frame(uniprot_id = c("A1", NA), symbol = c("R1", "R2"), evidence = c("IDA", NA),
               stringsAsFactors = FALSE)
  }, .package = "patliR")
  h <- patliR:::.ppi_highlight_table(c("A1", "B2", "C3"), c("x", "R2", "y"), "gpcr", "C3")
  expect_equal(h$is_highlighted, c(TRUE, TRUE, TRUE))
  expect_equal(h$source, c("GO:0004930", "GO:0004930", "highlight_set"))
  expect_equal(h$evidence[2], "symbol match")
  h0 <- patliR:::.ppi_highlight_table(c("A1", "B2"), c("x", "R2"), "none")
  expect_false(any(h0$is_highlighted))
})

test_that("max_targets keeps only targets with a partner among the selected ones", {
  targets <- data.frame(
    uniprot_id = c("a", "b", "c", "d", "e"), n_compounds = c(5, 4, 3, 2, 1),
    is_disease = c(TRUE, TRUE, FALSE, FALSE, FALSE), ppi_degree_all = c(1, 1, 2, 1, 0),
    stringsAsFactors = FALSE
  )
  ppi <- data.frame(uniprot_a = c("a", "b"), uniprot_b = c("c", "d"), stringsAsFactors = FALSE)
  s <- patliR:::.ppi_select_targets(targets, ppi, max_targets = 2)
  ## a and b are the top two but not linked to each other: they are replaced
  ## by c and d, which are not linked either -> nothing connected fits in 2
  expect_true(length(s$drawn) <= 2)
  sub <- ppi[ppi$uniprot_a %in% s$drawn & ppi$uniprot_b %in% s$drawn, ]
  expect_true(all(s$drawn %in% c(sub$uniprot_a, sub$uniprot_b)))
  expect_equal(s$n_no_partner, 1L)
  s4 <- patliR:::.ppi_select_targets(targets, ppi, max_targets = 4)
  expect_equal(s4$drawn, c("a", "b", "c", "d"))
  s_all <- patliR:::.ppi_select_targets(targets, ppi, max_targets = Inf)
  expect_false("e" %in% s_all$drawn)
})

test_that("disease_partners adds non-target disease genes linked to >= 2 drawn targets", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", disease_partners = 3,
                        disease_node = TRUE, view = "ppi", save = FALSE)
  nodes <- attr(p, "data")$nodes
  expect_equal(nodes$uniprot_id[nodes$kind == "disease_partner"], "Q00009")
  expect_equal(attr(p, "summary")$n_disease_partners, 1L)
  expect_null(attr(p, "data")$compounds)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("labels: top_n_labels limits the labelled targets, GPCRs are always named", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", top_n_labels = 1, save = FALSE)
  txt <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")) && "text" %in% names(l$data), p$layers)[[1]]$data
  labelled <- txt$uniprot_id[nzchar(txt$text)]
  expect_setequal(labelled, c("P00001", "P00002"))
  p0 <- plot_ppi_network(proj, condition = "A", top_n_labels = 0, highlight = "none", save = FALSE)
  txt0 <- Filter(function(l) inherits(l$geom, c("GeomText", "GeomTextRepel")) && "text" %in% names(l$data), p0$layers)[[1]]$data
  expect_false(any(nzchar(txt0$text)))
})

test_that("colour_by = 'module' uses network_module_membership", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  expect_error(plot_ppi_network(proj, condition = "A", colour_by = "module", save = FALSE),
               "network_module_membership")
  patliRResults(proj, "network_module_membership") <- data.frame(
    condition = "A", node_id = c("P00001", "P00002", "P00004"), node_type = "target",
    module_id = c("M1", "M1", "M2"), module_type = "cluster", stringsAsFactors = FALSE
  )
  p <- plot_ppi_network(proj, condition = "A", colour_by = "module", save = FALSE)
  nodes <- Filter(function(l) inherits(l$geom, "GeomPoint") && "fill_group" %in% names(l$data), p$layers)[[1]]$data
  expect_setequal(unique(nodes$fill_group), c("M1", "M2", "not clustered"))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot_ppi_network() validates its arguments", {
  skip_if_not_installed("ggplot2")
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  expect_error(plot_ppi_network(proj, condition = "A", max_targets = 1, save = FALSE), "max_targets")
  expect_error(plot_ppi_network(proj, condition = "A", score_threshold = 2000, save = FALSE), "score_threshold")
  expect_error(plot_ppi_network(proj, condition = "A", top_n_labels = -1, save = FALSE), "top_n_labels")
  expect_error(plot_ppi_network(proj, condition = "A", disease_partners = NA, save = FALSE), "disease_partners")
  expect_error(plot_ppi_network(proj, condition = "A", disease_node = NA, save = FALSE), "disease_node")
  expect_error(plot_ppi_network(proj, condition = "A", highlight_set = 1, save = FALSE), "highlight_set")
  expect_error(plot_ppi_network(proj, condition = "A", disease = c("D1", "D2"), save = FALSE), "disease")
  expect_error(plot_ppi_network(proj, condition = "A", disease = "nope", save = FALSE), "Available disease IDs")
  expect_error(plot_ppi_network(proj, condition = "B", save = FALSE), "not built")
  expect_error(plot_ppi_network(proj, condition = "A", view = "other", save = FALSE))
  proj2 <- proj
  patliRResults(proj2, "disease_genes") <- NULL
  expect_error(plot_ppi_network(proj2, condition = "A", disease = "D1", save = FALSE), "disease_genes")
})

test_that("save = TRUE writes a PNG and upserts one log row per (condition, disease, view, colour_by)", {
  skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  out <- withr::local_tempdir()
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", out_dir = out, width = 6, height = 5, dpi = 40)
  proj <- attr(p, "proj")
  log1 <- patliRResults(proj, "ppi_network_plot_log")
  expect_equal(nrow(log1), 1L)
  expect_true(file.exists(log1$path))
  expect_match(basename(log1$path), "^ppi_network_A_D1_compound_ppi[.]png$")
  expect_equal(log1$n_targets_drawn, 5L)
  expect_true(all(c("condition", "disease_id", "view", "colour_by", "p_hypergeom", "n_gpcr_edges_drawn") %in% names(log1)))
  ## same key again -> replaced, not appended
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", out_dir = out, width = 6, height = 5, dpi = 40)
  proj <- attr(p, "proj")
  expect_equal(nrow(patliRResults(proj, "ppi_network_plot_log")), 1L)
  ## other view -> second row
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", view = "ppi", out_dir = out, width = 6, height = 5, dpi = 40)
  proj <- attr(p, "proj")
  expect_equal(nrow(patliRResults(proj, "ppi_network_plot_log")), 2L)
  ## a non-default overlay is its own variant (file + row)
  p <- plot_ppi_network(proj, condition = "A", disease = "D1", view = "ppi", disease_partners = 2,
                        out_dir = out, width = 6, height = 5, dpi = 40)
  log3 <- patliRResults(attr(p, "proj"), "ppi_network_plot_log")
  expect_equal(nrow(log3), 3L)
  expect_setequal(log3$variant, c("default", "partners1"))
  expect_true(file.exists(file.path(out, "ppi_network_A_D1_ppi_partners1.png")))
})

test_that("engine = 'ggiraph' returns a girafe widget with the summary attached", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggiraph")
  withr::local_options(patliR.repel = FALSE)
  .ppi_local_mocks()
  proj <- .ppi_test_project()
  w <- plot_ppi_network(proj, condition = "A", engine = "ggiraph", save = FALSE)
  expect_s3_class(w, "girafe")
  expect_equal(attr(w, "summary")$n_targets_drawn, 5L)
})

test_that(".ppi_string_network() maps accessions, shares a STRING id's edges and deduplicates pairs", {
  testthat::local_mocked_bindings(
    .network_stringdb_aliases_map = function(...) {
      data.frame(string_protein_id = c("s1", "s1", "s2", "s3"), alias = c("A", "A2", "B", "C"),
                 stringsAsFactors = FALSE)
    },
    .ppi_string_links = function(...) {
      data.frame(protein1 = c("s1", "s2"), protein2 = c("s2", "s3"), combined_score = c(900L, 400L),
                 stringsAsFactors = FALSE)
    },
    .package = "patliR"
  )
  proj <- .test_project()
  res <- patliR:::.ppi_string_network(proj, c("A", "A2", "B", "C", "Z"), 9606, "12.0", 700)
  expect_equal(unname(res$map[c("A", "A2", "Z")]), c("s1", "s1", NA))
  expect_equal(res$n_universe, 3L)
  expect_equal(nrow(res$edges), 2L) # A-B and A2-B; the 400 link is filtered
  expect_setequal(paste(res$edges$uniprot_a, res$edges$uniprot_b), c("A B", "A2 B"))
})

test_that(".ppi_string_links() reads the flat file once, caches it, and reuses a lower-threshold cache", {
  proj <- .test_project()
  dir <- file.path(cacheDir(proj), "stringdb")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  gz <- file.path(dir, "9606.protein.links.v12.0.txt.gz")
  con <- gzfile(gz, "w")
  writeLines(c("protein1 protein2 combined_score",
               "s1 s2 900", "s2 s1 900", "s1 s3 500", "s3 s1 500", "s2 s3 750", "s3 s2 750"), con)
  close(con)
  l7 <- patliR:::.ppi_string_links(proj, 9606, "12.0", 700)
  expect_equal(nrow(l7), 2L)
  expect_true(all(l7$protein1 < l7$protein2))
  expect_true(file.exists(file.path(dir, "ppi_links_9606_12.0_700.rds")))
  unlink(gz)
  ## no flat file any more: a stricter threshold is served from the 700 cache
  l8 <- patliR:::.ppi_string_links(proj, 9606, "12.0", 800)
  expect_equal(nrow(l8), 1L)
})

test_that(".ppi_hypergeom() and .ppi_declutter() behave", {
  h <- patliR:::.ppi_hypergeom(c("a", "b", "c"), c("b", "c", "d", "e"), 100)
  expect_equal(h$k, 2L)
  expect_equal(h$p, stats::phyper(1, 4, 96, 3, lower.tail = FALSE))
  expect_true(is.na(patliR:::.ppi_hypergeom(character(0), "a", 100)$p))
  xy <- matrix(c(0, 0, 0.01, 0, 0.5, 0.5), ncol = 2, byrow = TRUE)
  out <- patliR:::.ppi_declutter(xy, min_dist = 0.1)
  d <- as.matrix(stats::dist(out))
  expect_true(min(d[upper.tri(d)]) >= 0.1 - 1e-6)
  expect_true(all(sqrt(rowSums(out^2)) <= 1 + 1e-9))
})
