#' @include AllGenerics.R internal.R adme_local.R
NULL

## Reduces adme_local()'s descriptors to 2D or 3D (PCA via base
## stats::prcomp(); UMAP optional via the `umap` package) and colors points
## by chemical family or any drug-likeness/route flag. Convex-hull "halos"
## per family: chull() in 2D; a translucent 3D alpha-hull mesh (plotly) or
## a 3-panel PC-pair matrix (static) in 3D. Overlaying target/protein
## shapes on the same space is not attempted -- proteins have no logP/TPSA,
## so a joint embedding is a separate design question (see ROADMAP.md).

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
#' Constant descriptors are dropped before scaling; at least `dims`
#' non-constant descriptors must remain.
#'
#' @inheritParams compounds
#' @param condition Character scalar, a single condition column of
#'   [binarizedMatrix()], or `NULL` (default) for every compound with
#'   [adme_local()] results regardless of condition.
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default). Combined with `condition` if both given (intersection).
#' @param dims `2` (default) or `3`. `3` projects onto three components and
#'   renders either a rotatable 3D scatter with a translucent hull mesh
#'   ("halo") per family (`engine = "plotly"`, needs the `plotly` package)
#'   or a 3-panel matrix of PC pairs with 2D hulls (`engine = "static"`).
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
#' @param engine `"ggiraph"` (default, if installed): interactive 2D plot
#'   with a tooltip per point. `"static"`: plain `ggplot2`. `"plotly"`:
#'   rotatable 3D (only meaningful with `dims = 3`; needs the `plotly`
#'   package).
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
#'   `attr(result, "proj")`. This project attribute is also returned when
#'   `save = FALSE`, retaining descriptor-exclusion logs without writing files.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' plot_chemical_space(proj, color_by = "ro5_pass", engine = "static", save = FALSE)
#' # three axes, one hull per family, rotatable:
#' plot_chemical_space(proj, dims = 3, color_by = "family", engine = "plotly", save = FALSE)
#' }
#'
#' @export
plot_chemical_space <- function(proj, condition = NULL, compound_ids = NULL,
                                 dims = 2,
                                 method = c("pca", "umap"),
                                 color_by = "family",
                                 show_hulls = TRUE, seed = NULL,
                                 engine = c("ggiraph", "static", "plotly"),
                                 save = TRUE, out_dir = NULL,
                                 width = 7, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  method <- match.arg(method)
  engine <- match.arg(engine)
  stopifnot(is.character(color_by), length(color_by) == 1)
  stopifnot(dims %in% c(2, 3))
  dims <- as.integer(dims)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_chemical_space}.")
  }
  if (method == "umap" && !requireNamespace("umap", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_chemical_space} needs the {.pkg umap} package (CRAN) for {.code method = \"umap\"}.",
      "i" = "{.code install.packages(\"umap\")}, or use {.code method = \"pca\"} (base R, no extra install)."
    ))
  }
  if (engine == "plotly" && dims != 3) {
    cli::cli_warn("{.code engine = \"plotly\"} only does something for {.code dims = 3}; using {.val static}.")
    engine <- "static"
  }
  if (engine == "plotly" && !requireNamespace("plotly", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.code engine = \"plotly\"} needs the {.pkg plotly} package (CRAN).",
      "i" = "{.code install.packages(\"plotly\")}, or use {.code engine = \"static\"} for a 3-panel PC matrix."
    ))
  }
  if (dims == 3 && engine == "ggiraph") engine <- "static"

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

  coords <- .chemical_space_coords(mat, method, seed, dims)
  adme$dim1 <- coords[, 1]
  adme$dim2 <- coords[, 2]
  if (dims == 3) adme$dim3 <- coords[, 3]

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
    ve <- summary(pca)$importance["Proportion of Variance", seq_len(dims)] * 100
    sprintf("PC%d (%.1f%% var.)", seq_len(dims), ve)
  } else {
    paste0("UMAP", seq_len(dims))
  }

  result <- if (dims == 3 && engine == "plotly") {
    .chemical_space_plotly3d(adme, axis_labels, color_by, is_categorical, draw_hulls, method)
  } else if (dims == 3) {
    .chemical_space_static3d(adme, axis_labels, color_by, is_categorical, draw_hulls, method)
  } else {
    p <- .chemical_space_ggplot(adme, axis_labels, color_by, is_categorical, draw_hulls, engine, method)
    if (engine == "static") {
      p
    } else if (requireNamespace("ggiraph", quietly = TRUE)) {
      ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
    } else {
      cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
      p
    }
  }

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ## scope goes into the file name and the log key: a per-condition or
    ## compound-subset plot must not overwrite the whole-project one.
    scope <- .chemical_space_scope(condition, compound_ids)
    tag <- paste0("chemical_space_", dims, "d_", method, "_", color_by,
                  if (nzchar(scope)) paste0("_", scope) else "")
    if (dims == 3 && engine == "plotly") {
      ## plotly widget cannot go through ggsave(); save an html file if we can.
      if (requireNamespace("htmlwidgets", quietly = TRUE)) {
        path <- file.path(out_dir, paste0(tag, ".html"))
        ## selfcontained needs pandoc; fall back to a sidecar _files/ dir.
        htmlwidgets::saveWidget(result, path, selfcontained = nzchar(Sys.which("pandoc")))
      } else {
        path <- NA_character_
        cli::cli_warn("Install {.pkg htmlwidgets} to save the 3D widget to disk; returning it un-saved.")
      }
    } else {
      p_static <- if (dims == 3) result else
        .chemical_space_ggplot(adme, axis_labels, color_by, is_categorical, draw_hulls, "static", method)
      path <- file.path(out_dir, paste0(tag, ".png"))
      ggplot2::ggsave(path, p_static, width = width, height = if (dims == 3) max(height, 5) else height, dpi = dpi, bg = "white")
    }
    log_cols <- c("compound_id", "dim1", "dim2", if (dims == 3) "dim3", "color_value")
    log_df <- adme[, log_cols]
    log_df$method <- method
    log_df$color_by <- color_by
    log_df$dims <- dims
    log_df$scope <- scope
    log_df$path <- path
    old_log <- patliRResults(proj, "chemical_space_log")
    if (!is.null(old_log) && nrow(old_log) > 0) {
      ## legacy logs (no dims/scope columns) are treated as whole-project scope
      if (!"scope" %in% names(old_log)) old_log$scope <- ""
      if (!"dims" %in% names(old_log)) old_log$dims <- ifelse("dim3" %in% names(old_log), 3, 2)
      same <- old_log$scope == scope & old_log$method == method &
        old_log$color_by == color_by & old_log$dims == dims
      old_log <- old_log[!same, , drop = FALSE]
      for (nm in setdiff(names(old_log), names(log_df))) log_df[[nm]] <- NA
      for (nm in setdiff(names(log_df), names(old_log))) old_log[[nm]] <- NA
      log_df <- rbind(old_log[, names(log_df), drop = FALSE], log_df)
    }
    patliRResults(proj, "chemical_space_log") <- log_df
    .write_results_csv(proj, "chemical_space_log", log_df)
  }
  if (save) .write_log_csv(proj)

  attr(result, "proj") <- proj
  result
}

#' Scope suffix for file names / log keys: `""` for the whole project,
#' otherwise the (sanitised) condition and/or a hash of `compound_ids`.
#' @keywords internal
.chemical_space_scope <- function(condition, compound_ids) {
  parts <- character(0)
  if (!is.null(condition)) parts <- c(parts, gsub("[^A-Za-z0-9_.-]+", "_", condition))
  if (!is.null(compound_ids)) parts <- c(parts, paste0("subset_", substr(rlang::hash(sort(unique(compound_ids))), 1, 8)))
  paste(parts, collapse = "_")
}

#' @keywords internal
.chemical_space_coords <- function(mat, method, seed, dims = 2) {
  varying <- vapply(seq_len(ncol(mat)), function(j) {
    x <- mat[, j]
    variance <- stats::var(x)
    is.finite(variance) && variance > 0
  }, logical(1))
  mat <- mat[, varying, drop = FALSE]
  if (ncol(mat) < dims || nrow(mat) < dims) {
    cli::cli_abort("Need at least {dims} non-constant descriptors and {dims} compounds for a {dims}D projection after dropping zero-variance descriptors.")
  }
  if (method == "pca") {
    pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
    out <- pca$x[, seq_len(min(dims, ncol(pca$x))), drop = FALSE]
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
  cfg <- umap::umap.defaults
  cfg$n_components <- dims
  fit <- umap::umap(scaled, config = cfg)
  fit$layout[, seq_len(dims), drop = FALSE]
}

#' 3D chemical space as a rotatable plotly scatter with a translucent
#' alpha-hull mesh ("halo") per family
#' @keywords internal
.chemical_space_plotly3d <- function(adme, axis_labels, color_by, is_categorical, draw_hulls, method) {
  pal <- .chemical_space_palette(sort(unique(adme$color_value)))
  p <- plotly::plot_ly()
  if (draw_hulls) {
    for (grp in names(pal)) {
      d <- adme[adme$color_value == grp, , drop = FALSE]
      if (nrow(d) < 4) next  # a 3D hull needs >= 4 non-coplanar points
      p <- plotly::add_trace(
        p, type = "mesh3d", x = d$dim1, y = d$dim2, z = d$dim3,
        alphahull = 8, opacity = 0.15, facecolor = rep(pal[[grp]], 1),
        color = I(pal[[grp]]), hoverinfo = "skip", showlegend = FALSE, name = grp
      )
    }
  }
  p <- plotly::add_markers(
    p, data = adme, x = ~dim1, y = ~dim2, z = ~dim3,
    color = ~color_value, colors = unlist(pal),
    text = ~compound_label, hoverinfo = "text",
    marker = list(size = 4, opacity = 0.9)
  )
  plotly::layout(
    p,
    title = paste0("Chemical space (", toupper(method), ", 3D) -- ", color_by),
    scene = list(
      xaxis = list(title = axis_labels[1]),
      yaxis = list(title = axis_labels[2]),
      zaxis = list(title = axis_labels[3])
    )
  )
}

#' 3D chemical space, static: a 3-panel matrix of PC pairs, each with 2D
#' family hulls -- every compound shown against all three axes at once
#' @keywords internal
.chemical_space_static3d <- function(adme, axis_labels, color_by, is_categorical, draw_hulls, method) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "Static 3D (the 3-panel PC matrix) needs the {.pkg patchwork} package.",
      "i" = "{.code install.packages(\"patchwork\")}, or use {.code engine = \"plotly\"}."
    ))
  }
  pairs <- list(c(1, 2), c(1, 3), c(2, 3))
  panel <- function(ij) {
    d <- adme
    d$px <- d[[paste0("dim", ij[1])]]
    d$py <- d[[paste0("dim", ij[2])]]
    g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$px, y = .data$py))
    if (draw_hulls) {
      hull <- do.call(rbind, lapply(split(d, d$color_value), function(x) {
        if (nrow(x) < 3) return(NULL)
        x[grDevices::chull(x$px, x$py), , drop = FALSE]
      }))
      if (!is.null(hull) && nrow(hull) > 0) {
        g <- g + ggplot2::geom_polygon(
          data = hull,
          ggplot2::aes(fill = .data$color_value, group = .data$color_value),
          alpha = 0.22, colour = NA
        )
      }
    }
    g +
      ggplot2::geom_point(ggplot2::aes(colour = .data$color_value), size = 2.4, alpha = 0.8) +
      ggplot2::labs(x = axis_labels[ij[1]], y = axis_labels[ij[2]], colour = color_by) +
      ggplot2::theme_minimal(base_size = 9)
  }
  panels <- lapply(pairs, panel)
  ## same single-legend fix as .chemical_space_ggplot()
  if (is_categorical) {
    pal <- .chemical_space_palette(sort(unique(stats::na.omit(adme$color_value))))
    if (length(pal) > 0) {
      panels <- lapply(panels, function(g) g +
        ggplot2::scale_colour_manual(values = pal, name = color_by, na.value = "grey60") +
        ggplot2::scale_fill_manual(values = pal, guide = "none", na.value = "grey60"))
    }
  }
  patchwork::wrap_plots(panels, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0("Chemical space (", toupper(method), ", 3 axes) -- ", color_by),
      subtitle = "every compound against all three principal components; filled areas are per-family hulls"
    ) &
    ggplot2::theme(legend.position = "right")
}

#' @keywords internal
.chemical_space_palette <- function(groups) {
  base <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
            "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")
  stats::setNames(rep(base, length.out = length(groups)), groups)
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
#' Visual language modeled on the classic "chemical space" figures (e.g.
#' the GDB/chemical-universe style plots) that draw two arrows crossing at
#' a shared origin instead of a boxed axis, and label each colored region
#' directly on the plot instead of relying only on a legend. Reconstructed
#' for patliR's own axes (PCA/UMAP of our compounds' physicochemical
#' descriptors), not copied --
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
  if (is_categorical && any(!is.na(adme$color_value))) {
    ## Repelled (when ggrepel is installed) so two families whose centroids
    ## nearly coincide do not print on top of each other; a white halo
    ## keeps them readable over the hulls and points.
    centroids <- stats::aggregate(cbind(dim1, dim2) ~ color_value, adme, mean)
    p <- p + .plot_text_layer(
      data = centroids,
      mapping = ggplot2::aes(x = .data$dim1, y = .data$dim2, label = .data$color_value, colour = .data$color_value),
      fontface = "bold", size = 3.4, show.legend = FALSE,
      repel_args = list(bg.colour = "white", bg.r = 0.15, box.padding = 0.5, min.segment.length = Inf,
                        max.overlaps = Inf, seed = 1),
      text_args = list()
    )
  }

  ## One legend, not two: hulls (fill) and points (colour) share one named
  ## palette, and only the colour legend is drawn. Left to ggplot2, the two
  ## scales get different hue ramps (hulls exist only for families with >= 3
  ## compounds, so the fill scale has fewer levels) and two "family" legends
  ## whose colours do not even match.
  pal <- if (is_categorical) .chemical_space_palette(sort(unique(stats::na.omit(adme$color_value)))) else NULL
  if (length(pal) > 0) {
    p <- p +
      ggplot2::scale_colour_manual(values = pal, name = color_by, na.value = "grey60") +
      ggplot2::scale_fill_manual(values = pal, guide = "none", na.value = "grey60") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  }

  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off") +
    ggplot2::labs(
      colour = color_by,
      title = paste0("Chemical space (", toupper(method), "), coloured by ", color_by),
      subtitle = paste0(axis_labels[1], "  |  ", axis_labels[2],
                         if (!is_categorical) "  -- continuous colour_by: hull outlines/labels not drawn" else "")
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right"
    )
}
