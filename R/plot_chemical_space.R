#' @include AllGenerics.R internal.R adme_local.R
NULL

## plot_chemical_space() -- reduces adme_local()'s physicochemical
## descriptors to 2D (PCA by default, base R stats::prcomp(), no new
## Suggests; UMAP optional via the `umap` package, new Suggests) and colors
## points by chemical family (compounds_classify()'s NPClassifier pathway)
## or by any drug-likeness rule/route flag. Loosely inspired by the classic
## "chemical space" figures that color-code whole molecule categories
## (e.g. Reymond & Awale, 2012) -- reconstructed for patliR's own domain
## (natural product extracts), not copied: our axes are physicochemical
## descriptors of *our* compounds, not a universal chemical-space model
## trained across DNA/peptides/graphenes/etc. Convex-hull outlines per
## family group (base R grDevices::chull(), no new dependency) stand in
## for that figure's colored "wedges". Design discussed with Uriel in chat
## (2026-08-11); "overlay_targets" (superimposing target/protein shapes on
## the same space, an idea in patliR_manual.md's "Planeado / no
## implementado" section) is deliberately NOT attempted here -- proteins
## do not have logP/TPSA in the same sense
## compounds do, so a joint embedding is a separate, harder design question
## left for a future round rather than guessed at now.

#' Chemical space plot: 2D projection of compounds, colored by family
#'
#' @description
#' Projects every compound's physicochemical descriptors (from
#' [adme_local()]) to 2D -- PCA by default, UMAP optional -- and colors
#' points by chemical family ([compounds_classify()]'s NPClassifier
#' pathway) or by any drug-likeness rule/route flag [adme_local()]
#' computed. Meant to answer "where do my extract's compounds sit, and do
#' compounds from the same chemical family cluster together" at a glance --
#' a companion, not a replacement, to the compound-by-compound detail in
#' [plot_admet_radar()].
#'
#' @section Method:
#' `method = "pca"` (default) standardizes (`scale = TRUE`) every numeric
#' descriptor and runs [stats::prcomp()] -- deterministic, no new
#' dependency. `method = "umap"` needs the `umap` package (CRAN, not
#' installed by default) -- non-linear, can separate clusters PCA
#' compresses together, but the axes have no direct physicochemical
#' meaning and results can shift slightly between runs unless `seed` is
#' set.
#'
#' @section Descriptors used:
#' `mw`, `logp`, `hbd`, `hba`, `tpsa`, `rotatable_bonds`, `amr`,
#' `fraction_csp3_approx`, `aromatic_proportion_approx`, `n_rings_approx`
#' -- every numeric column [adme_local()] computes except the derived
#' pass/fail rule flags and the BOILED-Egg call (those are outcomes, not
#' independent axes). Compounds with `NA` in any of these are excluded and
#' logged (`"plot_chemical_space_excluded_na"`), the same principle
#' [network_proximity()] uses for compounds STRINGdb cannot map.
#'
#' @inheritParams compounds
#' @param condition Character scalar, a single condition column of
#'   [binarizedMatrix()], or `NULL` (default) for every compound with
#'   [adme_local()] results regardless of condition.
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default). Combined with `condition` if both given (intersection).
#' @param method `"pca"` (default) or `"umap"`.
#' @param color_by `"family"` (default, needs [compounds_classify()] to
#'   have run -- colors by NPClassifier pathway), or any of `"ro5_pass"`,
#'   `"veber_pass"`, `"ghose_pass"`, `"egan_pass"`, `"oprea_pass"`,
#'   `"route_oral"`, `"gi_absorption"`, `"bbb_permeant"` (any column
#'   [adme_local()] produced -- validated against the actual columns
#'   present, not a fixed list, so a route you did not request in
#'   [adme_local()] raises a clear error rather than silently plotting
#'   nothing).
#' @param show_hulls Logical, default `TRUE`. Draw a convex-hull outline
#'   (`grDevices::chull()`, base R) around each `color_by` group with 3+
#'   points. Only really legible for a categorical `color_by` (i.e.
#'   `"family"`); ignored with a warning for continuous/logical `color_by`.
#' @param seed Integer or `NULL` (default). Only used by `method = "umap"`
#'   (UMAP's neighbor search has a stochastic component); restores the
#'   pre-call RNG state afterwards, same principle as `network_proximity()`
#'   / `network_module_robustness()`'s `seed` argument.
#' @param engine `"ggiraph"` (default, if installed): interactive plot with
#'   a tooltip per point. `"static"`: plain `ggplot2`.
#' @param save Logical, default `TRUE`. If `TRUE`, also writes a PNG.
#' @param out_dir Directory to write the PNG to (only used if
#'   `save = TRUE`). Defaults to `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()].
#'
#' @return A `girafe` htmlwidget (`engine = "ggiraph"`) or a `ggplot`
#'   object (`engine = "static"`). If `save = TRUE` (the default), the PNG
#'   path and PC/UMAP scores per compound are also recorded in
#'   `patliRResults(proj, "chemical_space_log")` /
#'   `results/chemical_space_log.csv`, retrievable via
#'   `attr(result, "proj")`.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' plot_chemical_space(proj, color_by = "ro5_pass", engine = "static", save = FALSE)
#' }
#'
#' @export
plot_chemical_space <- function(proj, condition = NULL, compound_ids = NULL,
                                 method = c("pca", "umap"),
                                 color_by = "family",
                                 show_hulls = TRUE, seed = NULL,
                                 engine = c("ggiraph", "static"),
                                 save = TRUE, out_dir = NULL,
                                 width = 7, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  method <- match.arg(method)
  engine <- match.arg(engine)
  stopifnot(is.character(color_by), length(color_by) == 1)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_chemical_space}.")
  }
  if (method == "umap" && !requireNamespace("umap", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_chemical_space} needs the {.pkg umap} package (CRAN) for {.code method = \"umap\"}.",
      "i" = "{.code install.packages(\"umap\")}, or use {.code method = \"pca\"} (base R, no extra install)."
    ))
  }

  adme <- patliRResults(proj, "adme_local")
  if (is.null(adme) || nrow(adme) == 0) {
    cli::cli_abort("No {.val adme_local} results found; run {.fn adme_local} first.")
  }

  scope_ids <- adme$compound_id
  if (!is.null(compound_ids)) scope_ids <- intersect(scope_ids, compound_ids)
  if (!is.null(condition)) {
    stopifnot(is.character(condition), length(condition) == 1)
    bin <- binarizedMatrix(proj)
    if (nrow(bin) == 0 || !condition %in% names(bin)) {
      cli::cli_abort("Condition {.val {condition}} not found in {.fn binarizedMatrix}; run {.fn prep_binarize} first, or omit {.arg condition}.")
    }
    present <- bin$compound_id[bin[[condition]] == 1]
    scope_ids <- intersect(scope_ids, present)
  }
  adme <- adme[adme$compound_id %in% scope_ids, , drop = FALSE]
  if (nrow(adme) < 3) {
    cli::cli_abort("Need at least 3 compounds with {.val adme_local} results in scope to compute {.arg method}; got {nrow(adme)}.")
  }

  desc_cols <- c("mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "amr",
                  "fraction_csp3_approx", "aromatic_proportion_approx", "n_rings_approx")
  desc_cols <- intersect(desc_cols, names(adme))
  mat_full <- as.matrix(adme[, desc_cols, drop = FALSE])
  complete <- stats::complete.cases(mat_full)
  if (any(!complete)) {
    proj <- .log_append(
      proj, step = "plot_chemical_space", id = adme$compound_id[!complete],
      message = "plot_chemical_space_excluded_na: missing one or more physicochemical descriptors; excluded from the projection"
    )
  }
  adme <- adme[complete, , drop = FALSE]
  mat <- mat_full[complete, , drop = FALSE]
  if (nrow(adme) < 3) {
    cli::cli_abort("Fewer than 3 compounds have complete descriptors after excluding {.val NA}s; cannot compute {.arg method}.")
  }

  coords <- .chemical_space_coords(mat, method, seed)
  adme$dim1 <- coords[, 1]
  adme$dim2 <- coords[, 2]

  color_info <- .chemical_space_color_values(proj, adme, color_by)
  adme$color_value <- color_info$value
  is_categorical <- color_info$categorical

  cmp <- compounds(proj)
  name_lookup <- stats::setNames(cmp$name, cmp$id)
  adme$compound_label <- paste0(
    adme$compound_id,
    ifelse(is.na(name_lookup[adme$compound_id]), "", paste0(" (", name_lookup[adme$compound_id], ")"))
  )

  draw_hulls <- show_hulls && is_categorical
  if (show_hulls && !is_categorical) {
    cli::cli_warn("{.arg show_hulls} ignored: {.arg color_by = {color_by}} is not categorical.")
  }

  axis_labels <- if (method == "pca") {
    pca <- attr(coords, "pca")
    ve <- summary(pca)$importance["Proportion of Variance", 1:2] * 100
    c(sprintf("PC1 (%.1f%% var.)", ve[1]), sprintf("PC2 (%.1f%% var.)", ve[2]))
  } else {
    c("UMAP1", "UMAP2")
  }

  p <- .chemical_space_ggplot(adme, axis_labels, color_by, is_categorical, draw_hulls, engine, method)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    p_static <- .chemical_space_ggplot(adme, axis_labels, color_by, is_categorical, draw_hulls, "static", method)
    path <- file.path(out_dir, paste0("chemical_space_", method, "_", color_by, ".png"))
    ggplot2::ggsave(path, p_static, width = width, height = height, dpi = dpi)
    log_df <- adme[, c("compound_id", "dim1", "dim2", "color_value")]
    log_df$method <- method
    log_df$color_by <- color_by
    log_df$path <- path
    patliRResults(proj, "chemical_space_log") <- log_df
    .write_results_csv(proj, "chemical_space_log", log_df)
  }
  .write_log_csv(proj)

  result <- if (engine == "static") {
    p
  } else if (requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  } else {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    p
  }
  if (save) attr(result, "proj") <- proj
  result
}

#' @keywords internal
.chemical_space_coords <- function(mat, method, seed) {
  if (method == "pca") {
    pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
    out <- pca$x[, 1:2, drop = FALSE]
    attr(out, "pca") <- pca
    return(out)
  }
  ## method == "umap"
  scaled <- scale(mat)
  old_seed <- if (!is.null(seed) && exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  if (!is.null(seed)) {
    set.seed(seed)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }
  fit <- umap::umap(scaled)
  fit$layout[, 1:2, drop = FALSE]
}

#' @keywords internal
.chemical_space_color_values <- function(proj, adme, color_by) {
  if (color_by == "family") {
    classified <- patliRResults(proj, "compounds_classified")
    if (is.null(classified) || nrow(classified) == 0) {
      cli::cli_abort(c(
        "No {.val compounds_classified} results found; run {.fn compounds_classify} first.",
        "i" = "Or pass a different {.arg color_by}, e.g. {.val ro5_pass}."
      ))
    }
    lookup <- stats::setNames(classified$pathway, classified$compound_id)
    val <- lookup[adme$compound_id]
    val[is.na(val)] <- "Unclassified"
    return(list(value = as.character(val), categorical = TRUE))
  }
  if (!color_by %in% names(adme)) {
    cli::cli_abort("{.arg color_by = {color_by}} is not a column of {.val adme_local} results (and is not {.val \"family\"}).")
  }
  val <- adme[[color_by]]
  categorical <- is.logical(val) || is.character(val)
  list(value = if (categorical) as.character(val) else val, categorical = categorical)
}

#' Build the chemical-space ggplot: origin-crossing arrow axes, filled
#' family "blobs" (convex hulls), direct on-plot family labels near each
#' group's centroid
#'
#' @description
#' Visual language modeled on the classic "chemical space" figures Uriel
#' referenced in chat (2026-08-11) -- e.g. the GDB/chemical-universe style
#' plots that draw two arrows crossing at a shared origin instead of a
#' boxed axis, and label each colored region directly on the plot instead
#' of relying only on a legend. Reconstructed for patliR's own axes (PCA/
#' UMAP of *our* compounds' physicochemical descriptors), not copied --
#' see the header comment of this file for why a literal reproduction
#' would not even make sense here (our "categories" are chemical families
#' within one extract, not whole molecule classes like DNA/peptides/
#' graphenes).
#' @keywords internal
.chemical_space_ggplot <- function(adme, axis_labels, color_by, is_categorical, draw_hulls, engine, method) {
  tooltip_txt <- sprintf("%s\n%s: %s", adme$compound_label, color_by, adme$color_value)
  origin <- c(mean(range(adme$dim1)), mean(range(adme$dim2)))
  pad <- 0.12 * c(diff(range(adme$dim1)), diff(range(adme$dim2)))
  xlim <- range(adme$dim1) + c(-1, 1) * pad[1]
  ylim <- range(adme$dim2) + c(-1, 1) * pad[2]

  p <- ggplot2::ggplot()

  ## Origin-crossing arrow axes instead of a boxed axis -- theme_void()
  ## below removes ggplot's own axis lines/ticks, these two geom_segment
  ## calls stand in for them the way the reference figure draws its axes.
  p <- p +
    ggplot2::annotate("segment", x = xlim[1], xend = xlim[2], y = origin[2], yend = origin[2],
                       arrow = grid::arrow(length = grid::unit(0.2, "cm")), colour = "grey20", linewidth = 0.5) +
    ggplot2::annotate("segment", x = origin[1], xend = origin[1], y = ylim[1], yend = ylim[2],
                       arrow = grid::arrow(length = grid::unit(0.2, "cm")), colour = "grey20", linewidth = 0.5)

  hulls <- NULL
  if (draw_hulls) {
    hulls <- do.call(rbind, lapply(split(adme, adme$color_value), function(d) {
      if (nrow(d) < 3) return(NULL)
      idx <- grDevices::chull(d$dim1, d$dim2)
      d[idx, , drop = FALSE]
    }))
    if (!is.null(hulls) && nrow(hulls) > 0) {
      p <- p + ggplot2::geom_polygon(
        data = hulls,
        ggplot2::aes(x = .data$dim1, y = .data$dim2, fill = .data$color_value, group = .data$color_value),
        alpha = 0.28, colour = NA
      )
    }
  }

  if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p <- p + ggiraph::geom_point_interactive(
      data = adme,
      ggplot2::aes(x = .data$dim1, y = .data$dim2, colour = .data$color_value, tooltip = tooltip_txt),
      size = 3.2, alpha = 0.75
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = adme, ggplot2::aes(x = .data$dim1, y = .data$dim2, colour = .data$color_value),
      size = 3.2, alpha = 0.75
    )
  }

  ## Direct on-plot family labels near each group's centroid, in addition
  ## to the legend -- this is the part of the reference figure that makes
  ## it readable at a glance without cross-referencing a legend.
  if (is_categorical) {
    centroids <- stats::aggregate(cbind(dim1, dim2) ~ color_value, adme, mean)
    p <- p + ggplot2::geom_text(
      data = centroids,
      ggplot2::aes(x = .data$dim1, y = .data$dim2, label = .data$color_value, colour = .data$color_value),
      fontface = "bold", size = 3.4, show.legend = FALSE
    )
  }

  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::labs(
      colour = color_by, fill = color_by,
      title = paste0("Chemical space (", toupper(method), "), coloured by ", color_by),
      subtitle = paste0(axis_labels[1], "  |  ", axis_labels[2],
                         if (!is_categorical) "  -- continuous colour_by: hull outlines/labels not drawn" else "")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      legend.position = "right"
    )
}
