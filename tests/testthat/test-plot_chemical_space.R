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
