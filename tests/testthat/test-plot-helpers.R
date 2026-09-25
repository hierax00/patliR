test_that(".plot_wrap() wraps to the width, keeps explicit line breaks and passes NULL/NA through", {
  long <- paste(rep("word", 40), collapse = " ")
  wrapped <- patliR:::.plot_wrap(long, 30)
  expect_true(all(nchar(strsplit(wrapped, "\n", fixed = TRUE)[[1]]) <= 30))
  expect_identical(gsub("\n", " ", wrapped, fixed = TRUE), long)

  two <- patliR:::.plot_wrap("first line\nsecond line", 100)
  expect_identical(two, "first line\nsecond line")
  expect_null(patliR:::.plot_wrap(NULL))
  expect_true(is.na(patliR:::.plot_wrap(NA_character_)))
  expect_length(patliR:::.plot_wrap(c("a b", "c d"), 10), 2)
})

test_that(".plot_wrap() breaks a space-free systematic name at hyphens", {
  name <- "2-(4-Ethenyl-4-methyl-3-(prop-1-en-2-yl)cyclohexyl)propan-2-ol"
  pieces <- strsplit(patliR:::.plot_wrap(name, 25), "\n", fixed = TRUE)[[1]]
  expect_gt(length(pieces), 1)
  expect_true(all(nchar(pieces) <= 25))
  expect_identical(paste(pieces, collapse = ""), name)
})

test_that(".plot_wrap_width() grows with the figure and shrinks with the font", {
  expect_gt(patliR:::.plot_wrap_width(12, 8), patliR:::.plot_wrap_width(8, 8))
  expect_gt(patliR:::.plot_wrap_width(8, 7.5), patliR:::.plot_wrap_width(8, 12))
  expect_identical(patliR:::.plot_wrap_width(0.1, 8), 20L)
})

test_that(".plot_truncate() shortens only long labels, to max_chars ending in '...'", {
  x <- c("Globulol", "Naphthalene, 1,2-dihydro-1,1,6-trimethyl-", NA)
  out <- patliR:::.plot_truncate(x, 28)
  expect_identical(out[1], "Globulol")
  expect_identical(nchar(out[2]), 28L)
  expect_true(endsWith(out[2], "..."))
  expect_true(startsWith(x[2], sub("[.]{3}$", "", out[2])))
  expect_true(is.na(out[3]))
})

test_that(".plot_text_layer() repels only when ggrepel is available and allowed", {
  testthat::skip_if_not_installed("ggplot2")
  d <- data.frame(x = 1:2, y = 1:2, lab = c("a", "b"))
  old <- options(patliR.repel = FALSE)
  on.exit(options(old), add = TRUE)
  plain <- patliR:::.plot_text_layer(
    data = d, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
    repel_args = list(max.overlaps = Inf), text_args = list(vjust = -1)
  )
  expect_true(inherits(plain$geom, "GeomText"))
  expect_identical(plain$aes_params$vjust, -1)

  testthat::skip_if_not_installed("ggrepel")
  options(patliR.repel = TRUE)
  rep <- patliR:::.plot_text_layer(
    data = d, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
    repel_args = list(max.overlaps = Inf), text_args = list(vjust = -1)
  )
  expect_true(inherits(rep$geom, "GeomTextRepel"))
  expect_null(rep$aes_params$vjust)
})

test_that(".plot_finish() saves theme_void() figures on an opaque white background", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("png")
  proj <- .test_project()
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point() + ggplot2::theme_void()
  out <- patliR:::.plot_finish(
    proj, p, name = "test_plot_log", filename = "void.png",
    log_row = data.frame(condition = "A", path = NA_character_), key_cols = "condition",
    save = TRUE, width = 2, height = 2, dpi = 30
  )
  path <- patliRResults(attr(out, "proj"), "test_plot_log")$path
  img <- png::readPNG(path)
  if (dim(img)[3] == 4) expect_true(all(img[, , 4] == 1))
  expect_equal(img[1, 1, 1:3], c(1, 1, 1))
})

test_that(".plot_log_backfill() adds a missing key column to an old log only", {
  proj <- .test_project()
  expect_identical(patliR:::.plot_log_backfill(proj, "x_log", "subset", ""), proj)
  patliRResults(proj, "x_log") <- data.frame(condition = "A", path = "p")
  proj <- patliR:::.plot_log_backfill(proj, "x_log", "subset", "")
  expect_identical(patliRResults(proj, "x_log")$subset, "")
  patliRResults(proj, "x_log") <- data.frame(condition = "A", subset = "keep")
  proj <- patliR:::.plot_log_backfill(proj, "x_log", "subset", "")
  expect_identical(patliRResults(proj, "x_log")$subset, "keep")
})

.lab_delta_e <- function(cols) {
  lab <- grDevices::convertColor(t(grDevices::col2rgb(cols) / 255), from = "sRGB", to = "Lab")
  sqrt(rowSums(diff(lab)^2))
}

test_that(".plot_contrast_palette() gives n distinct colours with contrasting neighbours", {
  for (n in c(1, 5, 10, 12, 13, 30)) {
    cols <- patliR:::.plot_contrast_palette(n)
    expect_length(cols, n)
    expect_false(anyDuplicated(cols) > 0)
  }
  cols <- patliR:::.plot_contrast_palette(12)
  ## every neighbour pair far apart in CIELAB, unlike the rainbow ramp
  expect_gt(min(.lab_delta_e(cols)), 60)
  expect_gt(min(.lab_delta_e(patliR:::.plot_contrast_palette(10))),
            min(.lab_delta_e(grDevices::rainbow(10))))
  expect_gt(min(.lab_delta_e(patliR:::.plot_contrast_palette(30))), 20)
  expect_identical(patliR:::.plot_contrast_palette(0), character(0))
  expect_identical(patliR:::.plot_contrast_palette(7, "default"), grDevices::rainbow(7))
  expect_error(patliR:::.plot_contrast_palette(3, "neon"))
})

test_that(".plot_contrast_palette('grey') is greys alternating dark / light", {
  cols <- patliR:::.plot_contrast_palette(10, "grey")
  rgb <- grDevices::col2rgb(cols)
  expect_true(all(rgb[1, ] == rgb[2, ] & rgb[2, ] == rgb[3, ]))
  expect_false(anyDuplicated(cols) > 0)
  steps <- diff(rgb[1, ])
  expect_true(all(sign(steps[-1]) == -sign(steps[-length(steps)])))
  expect_gt(min(abs(steps)), 60)
  expect_length(patliR:::.plot_contrast_palette(1, "grey"), 1)
})

test_that(".plot_disease_names() / .plot_disease_label() resolve names from the disease tables", {
  proj <- .test_project()
  expect_identical(patliR:::.plot_disease_label(proj, "D1"), "D1")
  patliRResults(proj, "disease_genes") <- data.frame(
    disease_id = c("MONDO_1", "GO_0006954"), disease_name = c("hypertensive disorder", "inflammatory response (GO:0006954)"),
    uniprot_id = "P1", stringsAsFactors = FALSE
  )
  patliRResults(proj, "targets_disease_profile") <- data.frame(
    disease_id = c("MONDO_1", "EFO_2"), disease_name = c("ignored", "asthma"), stringsAsFactors = FALSE
  )
  nm <- patliR:::.plot_disease_names(proj, c("MONDO_1", "EFO_2", "X"))
  expect_identical(unname(nm), c("hypertensive disorder", "asthma", NA))
  expect_identical(
    patliR:::.plot_disease_label(proj, c("MONDO_1", "GO_0006954", "X")),
    c("hypertensive disorder (MONDO_1)", "inflammatory response (GO:0006954)", "X")
  )
  expect_identical(patliR:::.plot_disease_label(proj, "MONDO_1", with_id = FALSE), "hypertensive disorder")
})
