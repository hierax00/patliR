## Shared fixture: a project with adme_local() results and a synthetic
## family classification (so color_by = "family" works without network).
.cs_project <- function() {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  proj <- adme_local(proj)
  ids <- compounds(proj)$id
  patliRResults(proj, "compounds_classified") <- data.frame(
    compound_id = ids,
    pathway = rep(c("Terpenoids", "Fatty acids", "Shikimates and Phenylpropanoids"),
                  length.out = length(ids)),
    superclass = "sc", class = "cl", isglycoside = FALSE,
    source = "test", fetch_date = Sys.Date(), stringsAsFactors = FALSE
  )
  proj
}

test_that("chemical-space projections drop constant descriptors before scaling", {
  mat <- cbind(mw = c(100, 200, 150, 300), logp = c(1, 3, 2, 1), hbd = 0)
  coords <- patliR:::.chemical_space_coords(mat, "pca", NULL)
  expected <- stats::prcomp(mat[, 1:2], center = TRUE, scale. = TRUE)
  expect_equal(unname(coords[, 1:2]), unname(expected$x))
  expect_equal(rownames(attr(coords, "pca")$rotation), c("mw", "logp"))
  for (method in c("pca", "umap")) {
    expect_error(patliR:::.chemical_space_coords(mat, method, NULL, dims = 3),
                 "at least 3 non-constant descriptors")
    expect_error(patliR:::.chemical_space_coords(mat[, c(1, 3)], method, NULL),
                 "at least 2 non-constant descriptors")
  }
})

test_that("UMAP receives finite scaled descriptors when a column is constant", {
  skip_if_not_installed("umap")
  mat <- cbind(mw = seq_len(20), logp = sin(seq_len(20)), hbd = 0)
  testthat::local_mocked_bindings(
    umap = function(d, config) {
      expect_equal(ncol(d), 2L)
      expect_true(all(is.finite(d)))
      list(layout = d)
    }, .package = "umap"
  )
  expect_equal(dim(patliR:::.chemical_space_coords(mat, "umap", 1)), c(20L, 2L))
})

test_that("chemical space renders unknown categories and retains exclusions without saving", {
  skip_if_not_installed("ggplot2")
  proj <- .test_project()
  adme <- data.frame(compound_id = paste0("c", 1:5), mw = c(100, 200, 150, 300, NA),
                     logp = c(1, 3, 2, 1, 4), hbd = 0, ro5_pass = NA)
  patliRResults(proj, "adme_local") <- adme
  before <- list.files(projectDir(proj), recursive = TRUE, full.names = TRUE)
  hashes <- tools::md5sum(before)
  p <- plot_chemical_space(proj, color_by = "ro5_pass", engine = "static", save = FALSE)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  points <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)[[1]]$data
  expect_equal(points$compound_id, paste0("c", 1:4))
  expect_true(all(is.na(points$color_value)))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))))
  log <- projectLog(attr(p, "proj"))
  expect_equal(log$id[grepl("plot_chemical_space_excluded_na", log$message)], "c5")
  expect_identical(list.files(projectDir(proj), recursive = TRUE, full.names = TRUE), before)
  expect_identical(tools::md5sum(before), hashes)
})

test_that("plot_chemical_space() 2D static returns a ggplot and logs PC scores", {
  skip_if_not_installed("ggplot2")
  proj <- .cs_project()
  p <- plot_chemical_space(proj, color_by = "family", engine = "static", save = TRUE)
  expect_s3_class(p, "ggplot")

  log_df <- patliRResults(attr(p, "proj"), "chemical_space_log")
  expect_true(all(c("compound_id", "dim1", "dim2", "color_value") %in% names(log_df)))
  expect_false("dim3" %in% names(log_df))
})

test_that("plot_chemical_space() dims = 3 static is a 3-panel patchwork with a dim3", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  proj <- .cs_project()
  p <- plot_chemical_space(proj, dims = 3, color_by = "family", engine = "static", save = TRUE)
  expect_s3_class(p, "patchwork")
  expect_true("dim3" %in% names(patliRResults(attr(p, "proj"), "chemical_space_log")))
})

test_that("plot_chemical_space() dims = 3 plotly is a plotly htmlwidget", {
  skip_if_not_installed("plotly")
  proj <- .cs_project()
  p <- plot_chemical_space(proj, dims = 3, color_by = "family", engine = "plotly", save = FALSE)
  expect_s3_class(p, "plotly")
})

test_that("plot_chemical_space() rejects bad dims and warns on plotly + dims = 2", {
  skip_if_not_installed("ggplot2")
  proj <- .cs_project()
  expect_error(plot_chemical_space(proj, dims = 4, save = FALSE))
  skip_if_not_installed("plotly")
  expect_warning(
    plot_chemical_space(proj, dims = 2, engine = "plotly", color_by = "family", save = FALSE),
    "only does something"
  )
})

test_that("plot_chemical_space() errors clearly without adme_local()", {
  skip_if_not_installed("ggplot2")
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(plot_chemical_space(proj, save = FALSE), "adme_local")
})

test_that("plot_chemical_space() keeps a compound-subset plot from overwriting the whole-project one", {
  skip_if_not_installed("ggplot2")
  proj <- .cs_project()
  sub_ids <- utils::head(compounds(proj)$id, 4)
  p_all <- plot_chemical_space(proj, color_by = "family", engine = "static", save = TRUE)
  proj2 <- attr(p_all, "proj")
  p_sub <- plot_chemical_space(proj2, compound_ids = sub_ids,
                               color_by = "family", engine = "static", save = TRUE)
  proj3 <- attr(p_sub, "proj")

  plots <- list.files(file.path(projectDir(proj3), "plots"), pattern = "[.]png$")
  expect_true("chemical_space_2d_pca_family.png" %in% plots)
  expect_equal(sum(grepl("^chemical_space_2d_pca_family_subset_", plots)), 1)

  log_df <- patliRResults(proj3, "chemical_space_log")
  expect_setequal(unique(log_df$scope), c("", unique(log_df$scope[nzchar(log_df$scope)])))
  expect_equal(length(unique(log_df$scope)), 2)
  ## re-plotting the same scope replaces its rows instead of duplicating them
  p_sub2 <- plot_chemical_space(proj3, compound_ids = sub_ids,
                                color_by = "family", engine = "static", save = TRUE)
  expect_equal(nrow(patliRResults(attr(p_sub2, "proj"), "chemical_space_log")), nrow(log_df))
})

test_that(".chemical_space_scope() is empty for the whole project and distinct per condition/subset", {
  expect_identical(patliR:::.chemical_space_scope(NULL, NULL), "")
  expect_identical(patliR:::.chemical_space_scope("EVEG-I", NULL), "EVEG-I")
  expect_false(identical(patliR:::.chemical_space_scope(NULL, c("a", "b")),
                         patliR:::.chemical_space_scope(NULL, c("a", "c"))))
  expect_identical(patliR:::.chemical_space_scope(NULL, c("b", "a")),
                   patliR:::.chemical_space_scope(NULL, c("a", "b")))
})

test_that("plot_chemical_space() draws a single family legend with matching hull and point colours", {
  testthat::skip_if_not_installed("ggplot2")
  proj <- .cs_project()
  p <- plot_chemical_space(proj, color_by = "family", engine = "static", save = FALSE)
  fill <- p$scales$get_scales("fill")
  colour <- p$scales$get_scales("colour")
  expect_identical(fill$guide, "none")
  expect_identical(fill$palette(3), colour$palette(3))
  built <- ggplot2::ggplot_build(p)
  expect_no_error(ggplot2::ggplot_gtable(built))
})

test_that(".chemical_space_label_key() numbers compounds by their compounds() row and falls back for missing CIDs", {
  cmp <- data.frame(id = c("C0001", "C0002", "C0003", "C0004"),
                    pubchem_id = c("7269", NA, "8089", ""),
                    name = c("Benzene, 1,2,4,5-tetramethyl-", "Isopropyl myristate", NA, "Squalane"),
                    stringsAsFactors = FALSE)
  adme <- data.frame(compound_id = c("C0003", "C0001", "C0004", "C9999"),
                     dim1 = c(0, 1, 2, 10), dim2 = c(0, 1, 0, 10),
                     color_value = c("Terpenoids", "Fatty acids", "Terpenoids", "Unclassified"),
                     stringsAsFactors = FALSE)
  key <- patliR:::.chemical_space_label_key(adme, cmp, "index", color_by = "family")
  expect_identical(key$compound_id, c("C0001", "C0003", "C0004", "C9999"))
  expect_identical(key$index, c(1L, 3L, 4L, 5L))
  expect_identical(key$plot_label, c("1", "3", "4", "5"))
  expect_true(all(c("index", "compound_id", "pubchem_id", "name", "family", "dim1", "dim2", "labeled") %in% names(key)))
  expect_true(all(key$labeled))

  pub <- patliR:::.chemical_space_label_key(adme, cmp, "pubchem")
  expect_identical(pub$plot_label, c("7269", "8089", "C0004", "C9999"))
  nm <- patliR:::.chemical_space_label_key(adme, cmp, "name")
  expect_identical(nm$plot_label[nm$compound_id == "C0003"], "C0003")
  expect_true(all(nchar(nm$plot_label) <= 18))
  expect_identical(patliR:::.chemical_space_label_key(adme, cmp, "id")$plot_label, key$compound_id)
})

test_that("label_top keeps the most isolated compounds", {
  cmp <- data.frame(id = paste0("c", 1:6), pubchem_id = NA, name = NA, stringsAsFactors = FALSE)
  adme <- data.frame(compound_id = paste0("c", 1:6),
                     dim1 = c(0, 0.1, 0.2, 0.1, 5, 0), dim2 = c(0, 0.1, 0, 0.2, 5, 9),
                     color_value = "a", stringsAsFactors = FALSE)
  iso <- patliR:::.chemical_space_isolation(adme[c("dim1", "dim2")])
  expect_setequal(order(-iso)[1:2], c(5L, 6L))
  key <- patliR:::.chemical_space_label_key(adme, cmp, "id", label_top = 2)
  expect_setequal(key$compound_id[key$labeled], c("c5", "c6"))
  expect_equal(sum(patliR:::.chemical_space_label_key(adme, cmp, "id", label_top = 50)$labeled), 6)
})

test_that(".chemical_space_file_tag() adds _labeled_<kind> only for labeled figures", {
  expect_identical(patliR:::.chemical_space_file_tag(2, "pca", "family"), "chemical_space_2d_pca_family")
  expect_identical(patliR:::.chemical_space_file_tag(2, "pca", "family", "EFLO-S", "index"),
                   "chemical_space_2d_pca_family_EFLO-S_labeled_index")
  expect_identical(patliR:::.chemical_space_file_tag(3, "umap", "ro5_pass", "", "pubchem"),
                   "chemical_space_3d_umap_ro5_pass_labeled_pubchem")
})

test_that("plot_chemical_space() validates label / label_top", {
  skip_if_not_installed("ggplot2")
  proj <- .cs_project()
  expect_error(plot_chemical_space(proj, label = "smiles", engine = "static", save = FALSE))
  expect_error(plot_chemical_space(proj, label = "id", label_top = 0, engine = "static", save = FALSE), "label_top")
  expect_error(plot_chemical_space(proj, label = "id", label_top = c(1, 2), engine = "static", save = FALSE), "label_top")
})

test_that("a labeled chemical space is saved beside the unlabeled one, with its own log rows and key CSV", {
  skip_if_not_installed("ggplot2")
  old <- options(patliR.repel = FALSE)
  on.exit(options(old), add = TRUE)
  proj <- .cs_project()
  p0 <- plot_chemical_space(proj, color_by = "family", engine = "static", save = TRUE)
  expect_null(attr(p0, "label_key"))
  p1 <- plot_chemical_space(attr(p0, "proj"), color_by = "family", engine = "static",
                            label = "index", label_top = 2, save = TRUE)
  key <- attr(p1, "label_key")
  expect_s3_class(key, "data.frame")
  expect_equal(sum(key$labeled), 2)
  expect_true(file.exists(attr(key, "path")))
  expect_identical(basename(attr(key, "path")), "chemical_space_2d_pca_family_labeled_index_key.csv")
  expect_equal(nrow(utils::read.csv(attr(key, "path"))), nrow(key))

  ## the labeled figure draws one text per selected compound (plus the family labels)
  txt <- Filter(function(l) inherits(l$geom, "GeomText"), p1$layers)
  expect_true(any(vapply(txt, function(l) is.data.frame(l$data) && nrow(l$data) == 2 &&
                           "plot_label" %in% names(l$data), logical(1))))

  proj2 <- attr(p1, "proj")
  plots <- list.files(file.path(projectDir(proj2), "plots"), pattern = "[.]png$")
  expect_true(all(c("chemical_space_2d_pca_family.png", "chemical_space_2d_pca_family_labeled_index.png") %in% plots))
  log_df <- patliRResults(proj2, "chemical_space_log")
  expect_setequal(unique(log_df$label), c("none", "index"))
  expect_equal(sum(log_df$label == "index"), sum(log_df$label == "none"))
  ## re-running the labeled figure replaces its own rows only
  p2 <- plot_chemical_space(proj2, color_by = "family", engine = "static", label = "index", save = TRUE)
  expect_equal(nrow(patliRResults(attr(p2, "proj"), "chemical_space_log")), nrow(log_df))
  ## a name-labeled figure carries the key but writes no CSV
  p3 <- plot_chemical_space(proj2, color_by = "family", engine = "static", label = "name", save = TRUE)
  expect_null(attr(attr(p3, "label_key"), "path"))
})

test_that("labels also work on the static 3-panel matrix", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  old <- options(patliR.repel = FALSE)
  on.exit(options(old), add = TRUE)
  proj <- .cs_project()
  p <- plot_chemical_space(proj, dims = 3, color_by = "family", engine = "static", label = "id", save = FALSE)
  expect_s3_class(p, "patchwork")
  expect_true("dim3" %in% names(attr(p, "label_key")))
  expect_no_error(ggplot2::ggplot_build(p[[1]]))
})
