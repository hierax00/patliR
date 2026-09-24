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
