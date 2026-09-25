## plot_disease_map() only READS already-populated results tables, so the
## fixture below writes network_edges / targets_disease_profile /
## disease_genes / targets_disease by hand, fully offline. Tests check the
## data behind the figure (classification, counts, top-N selection,
## isolated compounds, log), never pixels.

.dmap_ids <- c(T1 = "P35354", T2 = "P23219", T3 = "P08253", T4 = "P00533", T5 = "P04637")

.dmap_project <- function() {
  id <- .dmap_ids
  proj <- patliR_project(tempfile("patliR_dmap_"))
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

.dmap_plot <- function(...) suppressWarnings(suppressMessages(plot_disease_map(...)))

test_that("keyword classification follows list order and falls back to Other", {
  cl <- disease_map_classes()
  nm <- c("non-small cell lung carcinoma", "viral pneumonia", "hypertensive disorder",
          "inflammatory response", "non-autoimmune hemolytic anemia", "autoimmune hemolytic anemia",
          "hypercholanemia, familial", "Alzheimer disease", "asthma", "ocular hypertension",
          "osteoarthritis, knee", "rheumatoid arthritis", "chronic kidney disease", "psoriasis",
          "smoking initiation", NA, "")
  got <- .disease_map_classify(nm, cl)
  expect_equal(got, c("Neoplasm", "Infectious", "Cardiovascular", "Immune & inflammatory", "Hematologic",
                      "Immune & inflammatory", "Metabolic & endocrine", "Nervous system & psychiatric",
                      "Respiratory", "Eye & ear", "Musculoskeletal", "Immune & inflammatory",
                      "Renal & urogenital", "Skin", "Other", "Other", "Other"))
  ## case-insensitive, custom rules replace the defaults
  custom <- list(Lungs = c("lung", "asthma"), Brain = "alzheimer")
  expect_equal(.disease_map_classify(c("ASTHMA", "Alzheimer disease", "gout"), custom),
               c("Lungs", "Brain", "Other"))
  expect_true(all(c("Neoplasm", "Cardiovascular", "Skin") %in% names(cl)))
  pal <- .disease_map_palette(c(names(cl), "Other"))
  expect_equal(length(pal), length(cl) + 1L)
  expect_false(anyNA(pal))
})

test_that("associations are pooled, filtered by min_score and deduplicated", {
  proj <- .dmap_project()
  id <- .dmap_ids
  a <- .disease_map_associations(proj, unname(id[c("T1", "T2", "T3", "T4")]), min_score = 0.1)
  key <- paste(a$target_id, a$disease_id)
  expect_false(anyDuplicated(key) > 0)
  ## T3-D4 (0.05) and T2-F1 (0.01) dropped; unscored GO row kept
  expect_false(paste(id[["T3"]], "D4") %in% key)
  expect_false(paste(id[["T2"]], "F1") %in% key)
  go <- a[a$disease_id == "GO_0006954", ]
  expect_equal(nrow(go), 1)
  expect_true(is.na(go$association_score))
  expect_equal(go$disease_name, "inflammatory response")
  ## T1-F1 in disease_genes (0.7) and targets_disease (0.2): one row, max score, both sources
  t1f1 <- a[a$target_id == id[["T1"]] & a$disease_id == "F1", ]
  expect_equal(t1f1$association_score, 0.7)
  expect_equal(t1f1$source, "disease_genes+targets_disease")
  expect_equal(t1f1$disease_name, "hypertensive disorder")
  ## T5 is not requested, T4 has no association
  expect_false(any(a$target_id %in% id[c("T4", "T5")]))
  ## sources restricts the tables read
  only_prof <- .disease_map_associations(proj, unname(id), sources = "profile", min_score = 0)
  expect_setequal(unique(only_prof$source), "profile")
})

test_that("disease and area summaries count distinct targets and rank diseases", {
  proj <- .dmap_project()
  id <- .dmap_ids
  a <- .disease_map_associations(proj, unname(id[c("T1", "T2", "T3", "T4")]), min_score = 0.1)
  cl <- disease_map_classes()
  a$disease_class <- .disease_map_classify(a$disease_name, cl)
  a$focus <- .disease_map_focus(proj, a, NULL)
  s <- .disease_map_summaries(a, names(cl))
  d <- s$diseases
  expect_equal(d$disease_id, c("D1", "D2", "F1", "D3", "GO_0006954"))
  expect_equal(d$n_targets, c(2L, 2L, 2L, 1L, 1L))
  expect_equal(d$mean_score[1], 0.7)
  expect_equal(d$focus, c(FALSE, FALSE, TRUE, FALSE, TRUE))
  k <- stats::setNames(s$classes$n_targets, s$classes$disease_class)
  expect_equal(unname(k[c("Neoplasm", "Nervous system & psychiatric", "Cardiovascular",
                          "Respiratory", "Immune & inflammatory", "Skin")]), c(2L, 2L, 2L, 1L, 1L, 0L))
  expect_equal(utils::tail(s$classes$disease_class, 1), "Other")
  ## explicit focus / no focus
  expect_equal(sum(.disease_map_focus(proj, a, "D3")), 1)
  expect_equal(sum(.disease_map_focus(proj, a, character(0))), 0)
})

test_that("class view: drawn targets, isolated compound, attributes", {
  proj <- .dmap_project()
  id <- .dmap_ids
  p <- .dmap_plot(proj, condition = "A", save = FALSE)
  expect_s3_class(p, "ggplot")
  expect_null(attr(p, "proj"))
  cls <- attr(p, "classification")
  expect_true(all(c("target_id", "target_label", "disease_id", "disease_name", "disease_class",
                    "association_score", "source", "focus") %in% names(cls)))
  expect_equal(nrow(cls), 8)
  expect_false(id[["T4"]] %in% cls$target_id)
  expect_equal(attr(p, "disease_summary")$disease_id[1], "D1")
  expect_equal(attr(p, "class_summary")$disease_class[nrow(attr(p, "class_summary"))], "Other")

  rn <- .disease_map_right_nodes(cls, list(diseases = attr(p, "disease_summary"), classes = attr(p, "class_summary")),
                                 view = "class")
  ## five areas with targets + the two focus diseases, each right after its area
  expect_equal(sum(rn$nodes$type == "class"), 5)
  expect_equal(rn$nodes$node_id[match("F1", rn$nodes$node_id) - 1], "Cardiovascular")

  ct <- .disease_map_ct(proj, "A")
  built <- .disease_map_view_network(proj, ct, cls, list(diseases = attr(p, "disease_summary"), classes = attr(p, "class_summary")),
                                     "A", view = "class")
  cmp <- built$data$compounds
  ## C3's only target (T4) has no disease association -> isolated, listed last, no edge
  expect_equal(cmp$name[nrow(cmp)], "C3")
  expect_true(cmp$isolated[cmp$name == "C3"])
  expect_false(any(cmp$isolated[cmp$name != "C3"]))
  expect_setequal(built$data$targets$uniprot_id, unname(id[c("T1", "T2", "T3")]))
  expect_equal(built$n_compounds, 3)
})

test_that("target selection: max_targets / min_shared", {
  ct <- data.frame(compound_id = c("C1", "C2", "C3", "C1", "C2", "C1"),
                   uniprot_id = c("X", "X", "X", "Y", "Y", "Z"), weight = 0.8, stringsAsFactors = FALSE)
  links <- data.frame(target_id = c("X", "Y", "Z", "Z"), node_id = c("a", "a", "a", "b"), stringsAsFactors = FALSE)
  expect_equal(.disease_map_select_targets(ct, links, max_targets = 2)$uniprot_id, c("X", "Y"))
  expect_equal(.disease_map_select_targets(ct, links, min_shared = 2)$uniprot_id, c("X", "Y"))
  expect_equal(.disease_map_select_targets(ct, links, min_shared = 1)$uniprot_id, c("X", "Y", "Z"))
  ## a target without link is never a candidate
  expect_equal(.disease_map_select_targets(ct, links[links$target_id == "Z", ], max_targets = 5)$uniprot_id, "Z")
})

test_that("disease view keeps the top N diseases plus every focus disease", {
  proj <- .dmap_project()
  p <- .dmap_plot(proj, condition = "A", view = "disease", top_n_diseases = 1, save = FALSE)
  sums <- list(diseases = attr(p, "disease_summary"), classes = attr(p, "class_summary"))
  rn <- .disease_map_right_nodes(attr(p, "classification"), sums, view = "disease", top_n = 1)
  expect_setequal(rn$nodes$node_id, c("D1", "F1", "GO_0006954"))
  expect_equal(sum(rn$nodes$type == "focus"), 2)
  rn2 <- .disease_map_right_nodes(attr(p, "classification"), sums, view = "disease", top_n = 2)
  expect_setequal(rn2$nodes$node_id, c("D1", "D2", "F1", "GO_0006954"))
  ## no focus at all
  p0 <- .dmap_plot(proj, condition = "A", view = "disease", focus = character(0), save = FALSE)
  expect_false(any(attr(p0, "classification")$focus))
})

test_that("heatmap counts compounds' targets per area and focus disease", {
  proj <- .dmap_project()
  p <- .dmap_plot(proj, condition = "A", view = "heatmap", save = FALSE)
  expect_s3_class(p, "ggplot")
  cls <- attr(p, "classification")
  sums <- list(diseases = attr(p, "disease_summary"), classes = attr(p, "class_summary"))
  ct <- .disease_map_ct(proj, "A")
  hm <- .disease_map_heatmap_counts(ct, cls, sums, rows = "compound", row_ids = c("C1", "C2", "C3"))
  n <- function(r, c) hm$n[hm$row_id == r & hm$col_id == c]
  expect_equal(n("C1", "Neoplasm"), 2)           # T1, T2
  expect_equal(n("C2", "Nervous system & psychiatric"), 1)  # T1 only
  expect_equal(n("C1", "F1"), 2)                 # T1, T3 in the hypertension set
  expect_equal(n("C2", "GO_0006954"), 1)
  expect_true(all(hm$n[hm$row_id == "C3"] == 0)) # isolated compound: a row of zeros
  expect_setequal(unique(hm$block), c("Disease area", "Focus disease"))
  pt <- .dmap_plot(proj, condition = "A", view = "heatmap", heatmap_rows = "target", save = FALSE)
  expect_s3_class(pt, "ggplot")
})

test_that("argument validation", {
  proj <- .dmap_project()
  expect_error(plot_disease_map(proj, condition = "A", view = "nope", save = FALSE))
  expect_error(plot_disease_map(proj, condition = "A", top_n_diseases = 0, save = FALSE), "top_n_diseases")
  expect_error(plot_disease_map(proj, condition = "A", max_targets = 2.5, save = FALSE), "max_targets")
  expect_error(plot_disease_map(proj, condition = "A", min_shared = -1, save = FALSE), "min_shared")
  expect_error(plot_disease_map(proj, condition = "A", top_hub_n = -1, save = FALSE), "top_hub_n")
  expect_error(plot_disease_map(proj, condition = "A", min_score = 2, save = FALSE), "min_score")
  expect_error(plot_disease_map(proj, condition = "A", sources = "omim", save = FALSE), "sources")
  expect_error(plot_disease_map(proj, condition = "A", focus = 1, save = FALSE), "focus")
  expect_error(plot_disease_map(proj, condition = "A", disease_classes = list("x"), save = FALSE), "disease_classes")
  expect_error(plot_disease_map(proj, condition = "A", disease_classes = list(Bad = "(unclosed"), save = FALSE),
               "regular expression")
  expect_error(plot_disease_map(proj, condition = "Z", save = FALSE), "not built")
  empty <- proj
  patliRResults(empty, "targets_disease_profile") <- NULL
  patliRResults(empty, "disease_genes") <- NULL
  patliRResults(empty, "targets_disease") <- NULL
  expect_error(plot_disease_map(empty, condition = "A", save = FALSE), "targets_disease_profile")
})

test_that("save = TRUE writes the PNG and upserts one log row per (condition, view, subset)", {
  proj <- .dmap_project()
  out <- tempfile("dmap_out_")
  p <- .dmap_plot(proj, condition = "A", out_dir = out, width = 8, height = 6, dpi = 50)
  proj2 <- attr(p, "proj")
  log <- patliRResults(proj2, "disease_map_plot_log")
  expect_equal(nrow(log), 1)
  expect_equal(log$condition, "A")
  expect_equal(log$view, "class")
  expect_equal(log$subset, "t30")
  expect_true(file.exists(log$path))
  expect_equal(basename(log$path), "disease_map_A_class_t30.png")
  expect_equal(log$n_compounds, 3)
  expect_equal(log$n_targets, 3)
  ## same call again replaces its row; another view adds one
  p <- .dmap_plot(proj2, condition = "A", out_dir = out, width = 8, height = 6, dpi = 50)
  p <- .dmap_plot(attr(p, "proj"), condition = "A", view = "heatmap", out_dir = out, width = 8, height = 6, dpi = 50)
  log <- patliRResults(attr(p, "proj"), "disease_map_plot_log")
  expect_equal(nrow(log), 2)
  expect_setequal(log$view, c("class", "heatmap"))
})
