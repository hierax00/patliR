## plot_pathway_network(), plot_kegg_topology(), kegg_routes() and
## plot_kegg_routes() only READ result tables (network_enrichment,
## network_kegg_topology, network_edges, disease_genes), so every test
## here builds those tables by hand on an empty project -- fully offline.
## Figures are checked through their data (layers, attributes, log rows),
## never through pixels. KEGG gene IDs in the synthetic topology are
## made-up ("hsa:9000xx"), so symbol mapping falls back to the bare IDs.

.pw_enrichment <- function() {
  ## Entrez IDs of real genes (TP53 7157, EGFR 1956, TNF 7124, IL6 3569,
  ## IL1B 3553, PTGS2 5743, AKT1 207, MAPK1 5594, JUN 3725, FOS 2353)
  data.frame(
    condition = c(rep("A", 6), "B"),
    db = c("kegg", "kegg", "kegg", "kegg", "kegg", "kegg", "kegg"),
    ID = c("hsa1", "hsa2", "hsa3", "hsa4", "hsa5", "hsa6", "hsa1"),
    Description = c("Term one", "Term two", "Term three", "Term four", "Term five", "Not significant", "Term one"),
    GeneRatio = "1/10", BgRatio = "1/100",
    pvalue = c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.1, 1e-3),
    p.adjust = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.2, 1e-2),
    qvalue = 0.01,
    geneID = c("7157/1956/7124/3569", "7157/1956/7124/3553", "7157/1956/7124/5743",
               "207/5594/3725", "207/5594/2353", "7157/207", "7157/1956"),
    Count = c(4L, 4L, 4L, 3L, 3L, 2L, 2L),
    stringsAsFactors = FALSE
  )
}

.pw_edges <- function() {
  data.frame(
    condition = "A",
    compound_id = c("C1", "C1", "C1", "C2", "C2", "C2", "C3"),
    uniprot_id = c("P04637", "P00533", "P01375", "P31749", "P28482", "P05412", "P05231"),
    weight = c(0.9, 0.8, 0.7, 0.95, 0.6, 0.55, 0.5),
    stringsAsFactors = FALSE
  )
}

.pw_project <- function() {
  proj <- .test_project()
  patliRResults(proj, "network_enrichment") <- .pw_enrichment()
  patliRResults(proj, "network_edges") <- .pw_edges()
  proj
}

## --- enrichment map helpers -------------------------------------------------

test_that(".pathway_net_overlap_matrix() computes Jaccard and overlap coefficients", {
  sets <- list(a = c("1", "2", "3", "4"), b = c("3", "4", "5"), c = c("9"))
  j <- patliR:::.pathway_net_overlap_matrix(sets, "jaccard")
  o <- patliR:::.pathway_net_overlap_matrix(sets, "overlap")
  expect_equal(dim(j), c(3L, 3L))
  expect_equal(unname(diag(j)), c(1, 1, 1))
  expect_equal(j["a", "b"], 2 / 5)
  expect_equal(o["a", "b"], 2 / 3)
  expect_equal(j["a", "c"], 0)
  expect_true(isSymmetric(j))
})

test_that(".pathway_net_edges() keeps only pairs at or above the threshold, once each", {
  sets <- list(a = c("1", "2", "3", "4"), b = c("3", "4", "5"), c = c("4", "5"))
  m <- patliR:::.pathway_net_overlap_matrix(sets, "jaccard")
  e <- patliR:::.pathway_net_edges(m, sets, 0.4)
  expect_setequal(paste(e$from, e$to), c("a b", "b c"))
  expect_equal(e$n_shared[e$from == "a" & e$to == "b"], 2L)
  expect_equal(nrow(patliR:::.pathway_net_edges(m, sets, 0.99)), 0L)
  expect_named(patliR:::.pathway_net_edges(m[1, 1, drop = FALSE], sets[1], 0.1),
               c("from", "to", "similarity", "n_shared"))
})

test_that(".pathway_net_select_terms() filters condition and p cutoff, orders and caps", {
  enr <- .pw_enrichment()
  t <- patliR:::.pathway_net_select_terms(enr, "A", top_n = 3, p_cutoff = 0.05)
  expect_equal(t$ID, c("hsa1", "hsa2", "hsa3"))
  t_all <- patliR:::.pathway_net_select_terms(enr, "A", top_n = 50, p_cutoff = 0.05)
  expect_false("hsa6" %in% t_all$ID)
  expect_equal(nrow(t_all), 5L)
  dup <- rbind(enr, enr[1, ])
  expect_equal(sum(patliR:::.pathway_net_select_terms(dup, "A", 50, 0.05)$ID == "hsa1"), 1L)
})

test_that(".pathway_net_clusters() separates disconnected groups and leaves singletons NA", {
  ids <- c("a", "b", "c", "d", "e", "f")
  edges <- data.frame(from = c("a", "a", "b", "d"), to = c("b", "c", "c", "e"),
                      similarity = c(0.9, 0.8, 0.7, 0.9), stringsAsFactors = FALSE)
  set.seed(1)
  cl <- patliR:::.pathway_net_clusters(ids, edges)
  expect_equal(unname(cl[c("a", "b", "c")]), c(1L, 1L, 1L))
  expect_equal(unname(cl[c("d", "e")]), c(2L, 2L))
  expect_true(is.na(cl[["f"]]))
})

test_that(".pathway_net_disease_test() matches the hypergeometric tail", {
  sets <- list(t1 = c("1", "2", "3", "4"), t2 = c("5", "6"))
  universe <- as.character(1:20)
  res <- patliR:::.pathway_net_disease_test(sets, disease_genes = c("1", "2", "3", "7", "99"), universe = universe)
  expect_equal(res$disease_k, c(3L, 0L))
  expect_equal(res$disease_K[1], 4L) # "99" is outside the universe
  expect_equal(res$disease_p[1], stats::phyper(2, 4, 16, 4, lower.tail = FALSE))
  expect_equal(res$disease_p[2], 1)
  expect_equal(res$disease_padj, stats::p.adjust(res$disease_p, "BH"))
})

test_that(".pathway_net_compound_links() applies the hit minimum, top_n and the per-compound cap", {
  map <- data.frame(compound_id = c("c1", "c1", "c1", "c2", "c2", "c3"),
                    ENTREZID = c("1", "2", "3", "1", "2", "9"), stringsAsFactors = FALSE)
  sets <- list(t1 = c("1", "2", "3"), t2 = c("1", "2"), t3 = c("2", "3", "4", "5"))
  l <- patliR:::.pathway_net_compound_links(map, sets, min_hits = 2, top_n = 1, max_links = 2)
  expect_equal(unique(l$compound_id), "c1")
  expect_equal(nrow(l), 2L)
  expect_equal(l$ID, c("t1", "t2")) # share 1 each, before t3 (0.5)
  l2 <- patliR:::.pathway_net_compound_links(map, sets, min_hits = 2, top_n = 5, max_links = 5)
  expect_false("c3" %in% l2$compound_id)
  expect_equal(nrow(patliR:::.pathway_net_compound_links(map, sets, min_hits = 10)), 0L)
})

test_that(".pathway_slot_angles() spaces angles evenly in their circular order", {
  a <- patliR:::.pathway_slot_angles(c(0.1, 0.12, 0.11, 3))
  d <- diff(sort(a %% (2 * pi)))
  expect_equal(d, rep(pi / 2, 3), tolerance = 1e-9)
  expect_equal(order(c(0.1, 0.12, 0.11, 3)), order(a))
})

## --- plot_pathway_network() --------------------------------------------------

test_that("plot_pathway_network() validates its arguments", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .pw_project()
  expect_error(plot_pathway_network(proj, save = FALSE), "one condition at a time")
  expect_error(plot_pathway_network(proj, condition = "Z", save = FALSE), "Available")
  expect_error(plot_pathway_network(proj, condition = "A", db = "go", save = FALSE), "No .*go")
  expect_error(plot_pathway_network(proj, condition = "A", top_n = 0, save = FALSE), "top_n")
  expect_error(plot_pathway_network(proj, condition = "A", min_similarity = 2, save = FALSE), "min_similarity")
  expect_error(plot_pathway_network(proj, condition = "A", p_cutoff = 0, save = FALSE), "p_cutoff")
  expect_error(plot_pathway_network(proj, condition = "A", disease = c("x", "y"), save = FALSE), "disease")
  expect_error(plot_pathway_network(proj, condition = "A", p_cutoff = 1e-9, save = FALSE), "adjusted p")
  expect_error(plot_pathway_network(.test_project(), save = FALSE), "network_enrichment")
})

test_that("plot_pathway_network() draws the selected terms, overlap edges and clusters", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  proj <- .pw_project()
  p <- plot_pathway_network(proj, condition = "A", top_n = 10, min_similarity = 0.3, save = FALSE)
  expect_s3_class(p, "ggplot")
  terms <- attr(p, "terms")
  expect_setequal(terms$ID, c("hsa1", "hsa2", "hsa3", "hsa4", "hsa5"))
  e <- attr(p, "edges")
  expect_true(all(e$similarity >= 0.3))
  expect_setequal(paste(e$from, e$to), c("hsa1 hsa2", "hsa1 hsa3", "hsa2 hsa3", "hsa4 hsa5"))
  expect_equal(length(unique(stats::na.omit(terms$cluster))), 2L)
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot_pathway_network(proj, condition = "A", top_n = 2, save = FALSE)
  expect_equal(attr(p2, "terms")$ID, c("hsa1", "hsa2"))
  p3 <- plot_pathway_network(proj, condition = "A", clusters = FALSE, save = FALSE)
  expect_true(all(is.na(attr(p3, "terms")$cluster)))
})

test_that("plot_pathway_network() saves a PNG and upserts its log row", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(patliR.repel = FALSE)
  proj <- .pw_project()
  p <- plot_pathway_network(proj, condition = "A", width = 6, height = 5, dpi = 50)
  proj2 <- attr(p, "proj")
  log <- patliRResults(proj2, "pathway_network_plot_log")
  expect_equal(nrow(log), 1L)
  expect_equal(log$condition, "A")
  expect_equal(log$db, "kegg")
  expect_equal(log$disease, "none")
  expect_false(log$compounds)
  expect_true(file.exists(log$path))
  p <- plot_pathway_network(proj2, condition = "A", width = 6, height = 5, dpi = 50)
  expect_equal(nrow(patliRResults(attr(p, "proj"), "pathway_network_plot_log")), 1L)
})

test_that("plot_pathway_network() adds compounds and the disease layer", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  withr::local_options(patliR.repel = FALSE)
  proj <- .pw_project()
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "Test disease", uniprot_id = c("P04637", "P00533", "P01375"),
    gene_symbol = c("TP53", "EGFR", "TNF"), association_score = 0.8, source = "manual", fetched_at = "x",
    stringsAsFactors = FALSE
  )
  expect_error(plot_pathway_network(proj, condition = "A", disease = "nope", save = FALSE), "Available")
  p <- plot_pathway_network(proj, condition = "A", disease = "D1", show_compounds = TRUE, save = FALSE)
  terms <- attr(p, "terms")
  expect_true(all(c("disease_k", "disease_padj", "disease_enriched") %in% names(terms)))
  expect_equal(terms$disease_k[terms$ID == "hsa1"], 3L)
  expect_equal(terms$disease_k[terms$ID == "hsa4"], 0L)
  links <- attr(p, "compound_links")
  expect_true(all(links$n_hit >= 2))
  expect_true("C1" %in% links$compound_id)
  expect_no_error(ggplot2::ggplot_build(p))
})

## --- KEGG topology helpers ----------------------------------------------------

test_that(".kegg_topology_subtype_map() maps KGML subtypes to effect and mechanism", {
  m <- patliR:::.kegg_topology_subtype_map(c("activation", "inhibition", "expression", "repression", "phosphorylation",
                                             "dephosphorylation", "binding/association", "indirect effect",
                                             "compound", "ubiquitination", NA))
  expect_equal(m$effect, c("activation", "inhibition", "expression", "repression", "modification", "modification",
                           "binding/association", "indirect effect", "other", "modification", "other"))
  expect_equal(m$mode, c("direct", "direct", "direct", "direct", "phosphorylation", "dephosphorylation", "direct",
                         "indirect", "via compound", "other modification", "direct"))
})

.pw_topo_rows <- function(from, to, subtype, pid = "hsa09999", cond = "A") {
  data.frame(condition = cond, pathway_id = pid, pathway_title = paste("Pathway", pid),
             from_kegg = from, to_kegg = to, from_uniprot = sub("hsa:", "U", from), to_uniprot = sub("hsa:", "U", to),
             relation_type = "PPrel", relation_subtype = subtype, relation_value = NA_character_,
             restrict_to_network = FALSE, stringsAsFactors = FALSE)
}

test_that(".kegg_topology_collapse_relations() keeps one edge per pair with the dominant effect", {
  t <- .pw_topo_rows(c("hsa:1", "hsa:1", "hsa:2", "hsa:3"), c("hsa:2", "hsa:2", "hsa:3", "hsa:3"),
                     c("phosphorylation", "activation", "inhibition", "activation"))
  r <- patliR:::.kegg_topology_collapse_relations(t)
  expect_equal(nrow(r), 2L) # the self-relation hsa:3 -> hsa:3 is dropped
  expect_equal(r$effect[r$from == "hsa:1"], "activation")
  expect_equal(r$mode[r$from == "hsa:1"], "phosphorylation")
  expect_equal(r$subtypes[r$from == "hsa:1"], "activation; phosphorylation")
  expect_equal(r$pathways[1], "hsa09999")
})

test_that(".kegg_topology_paralog_groups() merges genes with identical relations", {
  e <- data.frame(from = c("a1", "a2", "b", "c"), to = c("b", "b", "d", "d"),
                  effect = "activation", mode = "direct", stringsAsFactors = FALSE)
  g <- patliR:::.kegg_topology_paralog_groups(e)
  expect_equal(g[["a1"]], g[["a2"]])
  expect_false(g[["b"]] == g[["c"]])
  expect_equal(patliR:::.kegg_topology_box_label(c("MAPK3", "MAPK1"), c(FALSE, TRUE)), "MAPK1/3")
  expect_equal(patliR:::.kegg_topology_box_label(c("AKT1", "AKT2", "AKT3", "PIK3CA", "PIK3CB"), max_chars = 12), "AKT1 +4")
})

test_that(".kegg_topology_edge_geometry() clips at the boxes and marks inhibitions with a bar", {
  nodes <- data.frame(box = c("a", "b", "c"), x = c(0, 2, 2), y = c(0, 0, 2), hw = 0.3, hh = 0.1)
  e <- data.frame(from_box = c("a", "a"), to_box = c("b", "c"), effect = c("activation", "inhibition"),
                  mode = "direct", subtypes = "x", stringsAsFactors = FALSE)
  g <- patliR:::.kegg_topology_edge_geometry(e, nodes)
  expect_equal(g$edges$head, c("arrow", "tee"))
  p1 <- g$paths[g$paths$edge_id == 1, ]
  expect_equal(p1$x[1], 0.32, tolerance = 1e-9)
  expect_equal(p1$x[2], 2 - 0.32, tolerance = 1e-9)
  expect_equal(nrow(g$tees), 1L)
  expect_equal(g$tees$edge_id, 2L)
})

.pw_topology_project <- function() {
  proj <- .test_project()
  t <- rbind(
    .pw_topo_rows(c("hsa:900001", "hsa:900002", "hsa:900003"), c("hsa:900002", "hsa:900003", "hsa:900004"),
                  c("activation", "inhibition", "phosphorylation"), pid = "hsa09001"),
    .pw_topo_rows(c("hsa:900005", "hsa:900006"), c("hsa:900006", "hsa:900007"),
                  c("activation", "expression"), pid = "hsa09002")
  )
  patliRResults(proj, "network_kegg_topology") <- t
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "A", compound_id = c("C1", "C2", "C1", "C3"),
    uniprot_id = c("U900001", "U900001", "U900002", "U900005"), weight = c(0.9, 0.7, 0.6, 0.8),
    stringsAsFactors = FALSE
  )
  proj
}

test_that("plot_kegg_topology() picks the most-hit pathways and summarises them", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .pw_topology_project()
  p <- plot_kegg_topology(proj, top_n_pathways = 1, save = FALSE)
  s <- attr(p, "pathway_summary")
  expect_equal(s$pathway_id, "hsa09001")
  expect_equal(s$n_genes, 4L)
  expect_equal(s$n_hit, 2L)
  expect_equal(s$share_hit, 0.5)
  expect_equal(s$n_compounds, 2L)
  nodes <- attr(p, "nodes")
  expect_equal(nodes$max_weight[nodes$box == "hsa:900001"], 0.9)
  expect_equal(nodes$n_compounds[nodes$box == "hsa:900001"], 2L)
  expect_false(nodes$is_target[nodes$box == "hsa:900004"])
  expect_no_error(ggplot2::ggplot_build(p))
  p2 <- plot_kegg_topology(proj, save = FALSE)
  expect_equal(attr(p2, "pathway_summary")$pathway_id, c("hsa09001", "hsa09002"))
  expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("plot_kegg_topology() validates arguments and explains a missing pathway", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .pw_topology_project()
  expect_error(plot_kegg_topology(proj, pathway_id = "hsa04020", save = FALSE), "restrict_to_network = FALSE")
  expect_error(plot_kegg_topology(proj, top_n_pathways = 0, save = FALSE), "top_n_pathways")
  expect_error(plot_kegg_topology(proj, pathway_id = 1, save = FALSE), "pathway_id")
  expect_error(plot_kegg_topology(proj, width = -1, save = FALSE), "width")
  expect_error(plot_kegg_topology(proj, condition = "Z", save = FALSE), "Available")
  expect_error(plot_kegg_topology(.test_project(), save = FALSE), "network_kegg_topology")
})

test_that("plot_kegg_topology() saves a PNG and logs it keyed by pathways", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .pw_topology_project()
  p <- plot_kegg_topology(proj, pathway_id = "hsa09002", dpi = 50)
  log <- patliRResults(attr(p, "proj"), "kegg_topology_plot_log")
  expect_equal(log$pathways, "hsa09002")
  expect_equal(log$condition, "A")
  expect_true(file.exists(log$path))
})

## --- kegg_routes() / plot_kegg_routes() ------------------------------------------

## sources: g1 (U900001), g2 (U900002), g9 (U900009, also a disease gene)
## g1 --> g3 --| g5 and g1 --> g4 --| g5  (two equally short routes, net "-")
## g2 --- g6 --> g7                       (binding hop: sign "?")
## g8 --> g1                              (g8 unreachable from the targets)
.pw_routes_project <- function() {
  proj <- .test_project()
  k <- function(i) paste0("hsa:9000", sprintf("%02d", i))
  t <- rbind(
    .pw_topo_rows(c(k(1), k(3), k(1), k(4)), c(k(3), k(5), k(4), k(5)),
                  c("activation", "inhibition", "activation", "inhibition"), pid = "hsa09001"),
    .pw_topo_rows(c(k(2), k(6), k(8)), c(k(6), k(7), k(1)),
                  c("binding/association", "activation", "activation"), pid = "hsa09002"),
    .pw_topo_rows(k(9), k(3), "activation", pid = "hsa09002")
  )
  patliRResults(proj, "network_kegg_topology") <- t
  patliRResults(proj, "network_edges") <- data.frame(
    condition = "A", compound_id = c("C1", "C2", "C3", "C1"),
    uniprot_id = c("U900001", "U900001", "U900002", "U900009"), weight = c(0.9, 0.8, 0.7, 0.6),
    stringsAsFactors = FALSE
  )
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = "D1", disease_name = "Test disease", uniprot_id = c("U900005", "U900007", "U900008", "U900009"),
    gene_symbol = NA_character_, association_score = c(0.9, 0.5, 0.95, 0.3), source = "manual", fetched_at = "x",
    stringsAsFactors = FALSE
  )
  proj
}

test_that("kegg_routes() finds shortest routes, alternatives and propagates signs", {
  proj <- .pw_routes_project()
  r <- kegg_routes(proj, disease = "D1")
  expect_true(all(c("route_id", "n_hops", "path", "sign", "pathways") %in% names(r)))
  to5 <- r[r$destination == "hsa:900005", , drop = FALSE]
  expect_equal(nrow(to5), 2L) # two equally short alternatives
  expect_equal(to5$route_id, c("R1", "R1b"))
  expect_true(all(to5$n_hops == 2L))
  expect_true(all(to5$sign == "-")) # activation x inhibition
  to7 <- r[r$destination == "hsa:900007", , drop = FALSE]
  expect_equal(to7$source, "hsa:900002")
  expect_equal(to7$sign, "?") # binding hop has no sign
  expect_false("hsa:900008" %in% r$destination) # only reachable against the edge direction
  expect_equal(attr(r, "direct_hits")$kegg, "hsa:900009")
  hops <- attr(r, "hops")
  expect_equal(hops$pathways[hops$route_id == "R1"], c("hsa09001", "hsa09001"))
  ## destination ranking: the most associated reachable gene first
  expect_equal(r$destination[1], "hsa:900005")
})

test_that("kegg_routes() honours n_alternatives, min_hops, to/from and undirected_binding", {
  proj <- .pw_routes_project()
  r1 <- kegg_routes(proj, disease = "D1", n_alternatives = 1)
  expect_equal(sum(r1$destination == "hsa:900005"), 1L)
  r2 <- kegg_routes(proj, disease = "D1", min_hops = 3)
  expect_equal(nrow(r2), 0L)
  r3 <- kegg_routes(proj, to = "hsa:900007")
  expect_equal(unique(r3$destination), "hsa:900007")
  r4 <- kegg_routes(proj, disease = "D1", from = "U900001")
  expect_true(all(r4$source == "hsa:900001"))
  ## binding/association is travelled against its stored direction only
  ## with undirected_binding = TRUE
  rel <- data.frame(from = c("x", "y"), to = c("y", "z"), effect = c("binding/association", "activation"),
                    pathways = "p", stringsAsFactors = FALSE)
  e_on <- patliR:::.kegg_routes_edges(rel, undirected_binding = TRUE)
  e_off <- patliR:::.kegg_routes_edges(rel, undirected_binding = FALSE)
  expect_true("y x" %in% paste(e_on$from, e_on$to))
  expect_false("y x" %in% paste(e_off$from, e_off$to))
  expect_false("z y" %in% paste(e_on$from, e_on$to)) # directed relations are never reversed
  expect_equal(e_on$sign[e_on$from == "y" & e_on$to == "z"], 1L)
  expect_true(is.na(e_on$sign[e_on$from == "x" & e_on$to == "y"]))
  expect_true(e_on$reversed[e_on$from == "y" & e_on$to == "x"])
  rp <- kegg_routes(proj, disease = "D1", rank_by = "pair")
  expect_equal(nrow(rp[rp$alternative == 1, ]), 2L)
})

test_that("kegg_routes() returns an empty table when nothing is reachable, and validates arguments", {
  proj <- .pw_routes_project()
  r <- kegg_routes(proj, to = "hsa:900008")
  expect_equal(nrow(r), 0L)
  expect_named(r, names(patliR:::.kegg_routes_empty()))
  expect_error(plot_kegg_routes(proj, routes = r, save = FALSE), "No route")
  expect_error(kegg_routes(proj), "destinations")
  expect_error(kegg_routes(proj, disease = "nope"), "disease_genes")
  expect_error(kegg_routes(proj, disease = "D1", max_len = 0), "max_len")
  expect_error(kegg_routes(proj, disease = "D1", min_hops = 5, max_len = 2), "min_hops")
  expect_error(kegg_routes(proj, disease = "D1", rank_by = "x"))
  expect_error(kegg_routes(proj, to = NA_character_), "to")
  expect_warning(kegg_routes(proj, to = c("hsa:900005", "GHOST")), "not in the KEGG topology")
  expect_error(kegg_routes(.test_project(), disease = "D1"), "network_kegg_topology")
})

test_that("kegg_routes_table() gives one compact line per route", {
  proj <- .pw_routes_project()
  r <- kegg_routes(proj, disease = "D1")
  tab <- kegg_routes_table(r)
  expect_equal(nrow(tab), sum(r$alternative == 1))
  expect_match(tab$route_path[tab$route == "R1"], "-->.*--\\|")
  expect_match(tab$pathways_per_hop[tab$route == "R1"], "hsa09001")
  expect_equal(nrow(kegg_routes_table(r, alternatives = TRUE)), nrow(r))
  expect_error(kegg_routes_table(data.frame(x = 1)), "kegg_routes")
})

test_that(".kegg_routes_octilinear() uses horizontal, 45-degree and vertical segments only", {
  pts <- patliR:::.kegg_routes_octilinear(0, 0, 4, 1)
  d <- diff(pts)
  ok <- abs(d[, 2]) < 1e-9 | abs(d[, 1]) < 1e-9 | abs(abs(d[, 1]) - abs(d[, 2])) < 1e-9
  expect_true(all(ok))
  pts2 <- patliR:::.kegg_routes_octilinear(0, 0, 1, 3)
  d2 <- diff(pts2)
  expect_true(all(abs(d2[, 2]) < 1e-9 | abs(d2[, 1]) < 1e-9 | abs(abs(d2[, 1]) - abs(d2[, 2])) < 1e-9))
  expect_equal(pts2[nrow(pts2), ], c(1, 3))
})

test_that("plot_kegg_routes() draws the routes as lines and logs the PNG", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .pw_routes_project()
  p <- plot_kegg_routes(proj, disease = "D1", dpi = 50)
  expect_s3_class(p, "ggplot")
  drawn <- attr(p, "routes")
  expect_true(all(drawn$alternative == 1))
  paths <- Filter(function(l) inherits(l$geom, "GeomPath"), p$layers)[[1]]$data
  expect_setequal(unique(paths$route_id), drawn$route_id)
  log <- patliRResults(attr(p, "proj"), "kegg_routes_plot_log")
  expect_equal(log$target_set, "D1")
  expect_true(file.exists(log$path))
  expect_error(plot_kegg_routes(proj, routes = data.frame(x = 1), save = FALSE), "kegg_routes")
  expect_error(plot_kegg_routes(proj, disease = "D1", top_n_routes = 0, save = FALSE), "top_n_routes")
})
