#' @include AllGenerics.R internal.R
NULL

## Axis definitions for the bioavailability radar (plot_admet_radar()).
## `zone_min`/`zone_max` (the shaded "pink zone") are the published optimal
## physicochemical-space boundaries from Daina, A. & Zoete, V. (2017),
## "A BOILED-Egg To Predict Gastrointestinal Absorption and Brain
## Penetration of Small Molecules", and the companion SwissADME paper
## (Daina et al. 2017, Sci. Rep. 7:42717) that introduced this exact radar.
## `axis_min`/`axis_max` (the outer scale of each spoke) are *our own*
## choice for display purposes, wide enough to comfortably show the pink
## zone plus typical natural-product values -- they are NOT claimed to
## match SwissADME's own (undocumented) axis limits pixel-for-pixel.
##
## Axes are listed in the order they are drawn clockwise starting at 12
## o'clock: LIPO, SIZE, POLAR, INSOLU, INSATU, FLEX.
.radar_axes <- data.frame(
  axis = c("LIPO", "SIZE", "POLAR", "INSOLU", "INSATU", "FLEX"),
  label = c(
    "LIPO\n(XLOGP3)", "SIZE\n(MW)", "POLAR\n(TPSA)",
    "INSOLU\n(Log S)", "INSATU\n(Fsp3, approx.)", "FLEX\n(rotatable bonds)"
  ),
  property = c("logp", "mw", "tpsa", "logs_esol", "fraction_csp3_approx", "rotatable_bonds"),
  axis_min = c(-2, 0, 0, -10, 0, 0),
  axis_max = c(7, 700, 200, 2, 1, 20),
  zone_min = c(-0.7, 150, 20, -6, 0.25, 0),
  zone_max = c(5.0, 500, 130, 0, 1, 9),
  stringsAsFactors = FALSE
)

#' Bioavailability radar ("spiderweb") plot, SwissADME-style
#'
#' @description
#' Draws a six-axis radar plot -- lipophilicity (LIPO), molecular size
#' (SIZE), polarity (POLAR), solubility (INSOLU), saturation (INSATU), and
#' flexibility (FLEX) -- with the published "pink zone" of optimal oral
#' drug-likeness shaded behind the compound's own polygon, the same idea as
#' SwissADME's Bioavailability Radar. Needs [adme_local()] to have been run
#' first (that is where every one of these six properties is computed).
#'
#' **One plot per compound is returned** (a named list, keyed by
#' `compound_id`) rather than a single faceted plot -- with 6+ compounds a
#' facet grid makes every axis label unreadably small and cramped; a full
#' -size hexagon per compound is legible and closer to how SwissADME itself
#' presents one molecule at a time.
#'
#' @section Geometry note:
#' This is built with plain Cartesian trigonometry (`x = r*sin(theta)`,
#' `y = r*cos(theta)`) plus `ggplot2::coord_fixed()`, not
#' `ggplot2::coord_polar()`. `coord_polar()` interpolates *along the arc*
#' between two data points, so a hexagon drawn with it comes out with
#' visibly curved edges and can look rotated/off-center depending on the
#' panel's x-range; drawing the hexagon's vertices directly in Cartesian
#' space gives straight edges and full control over where axis 1 sits (12
#' o'clock here, going clockwise).
#'
#' @section On the approximations involved (please read):
#' The pink zone boundaries are the published ones (see `.radar_axes` in
#' the package source for the exact reference); the outer axis scale is
#' our own display choice, not a literature value. INSOLU (`logs_esol`) and
#' INSATU (`fraction_csp3_approx`) both inherit the approximations
#' documented in [adme_local()] -- in particular, INSATU is a rough
#' aliphatic-vs-aromatic-carbon proxy from the SMILES string, not a true
#' Fsp3 from hybridization perception. Treat this plot as a good first
#' read, not a certified reproduction of SwissADME's own figure.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound with [adme_local()] results.
#' @param engine `"ggiraph"` (default, if installed): interactive plot with
#'   a tooltip per vertex showing the raw value. `"static"`: plain
#'   `ggplot2`, no extra dependency beyond `ggplot2` itself.
#' @param save Logical, default `TRUE`. If `TRUE`, also writes two PNG files
#'   per compound to `out_dir` -- one with axis labels
#'   (`<compound_id>_labeled.png`) and one without
#'   (`<compound_id>_plain.png`), the latter meant for tiling into reports
#'   or side-by-side comparisons where six repeated axis labels are just
#'   noise.
#' @param out_dir Directory to write the PNG files to (only used if
#'   `save = TRUE`). Defaults to
#'   `file.path(projectDir(proj), "plots", "admet_radar")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] for the saved PNGs.
#'
#' @return A named `list` (one element per compound, named by
#'   `compound_id`) of either `girafe` htmlwidgets (`engine = "ggiraph"`)
#'   or `ggplot` objects (`engine = "static"`). Print/plot each element
#'   individually, e.g. `radars[["C0001"]]`, or `for (p in radars) print(p)`
#'   to see them all. If `save = TRUE` (the default), the PNG paths are also
#'   recorded in `patliRResults(proj, "admet_radar_log")` -- but note this
#'   function itself only returns the plot list, not the updated `proj`;
#'   grab the log from `attr(radars, "proj")` if you need the updated
#'   project object without re-running [adme_local()].
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' radars <- plot_admet_radar(proj, engine = "static", save = FALSE)
#' radars[["C0001"]]
#' }
#'
#' @export
plot_admet_radar <- function(proj, compound_ids = NULL, engine = c("ggiraph", "static"),
                              save = TRUE, out_dir = NULL,
                              width = 6, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_admet_radar}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }

  adme <- patliRResults(proj, "adme_local")
  if (is.null(adme)) {
    cli::cli_abort("No {.val adme_local} results found; run {.fn adme_local} first.")
  }
  if (!is.null(compound_ids)) adme <- adme[adme$compound_id %in% compound_ids, , drop = FALSE]
  if (nrow(adme) == 0) {
    cli::cli_abort("No matching compounds with {.val adme_local} results.")
  }

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots", "admet_radar")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  }

  cmp <- compounds(proj)
  name_lookup <- stats::setNames(cmp$name, cmp$id)
  n_axes <- nrow(.radar_axes)
  angles <- 2 * pi * (seq_len(n_axes) - 1) / n_axes # axis 1 at 12 o'clock, clockwise
  zone_xy <- .radar_zone_xy(angles)

  save_log <- vector("list", nrow(adme))

  plots <- stats::setNames(lapply(seq_len(nrow(adme)), function(i) {
    id <- adme$compound_id[i]
    label <- paste0(
      id,
      if (!is.na(name_lookup[id])) paste0(" (", name_lookup[id], ")") else ""
    )
    vertices <- .radar_compound_xy(adme[i, , drop = FALSE], angles)

    if (save) {
      p_labeled <- .radar_ggplot_one(vertices, zone_xy, label, engine = "static", show_labels = TRUE, fig_width = width)
      p_plain <- .radar_ggplot_one(vertices, zone_xy, label, engine = "static", show_labels = FALSE, fig_width = width)
      path_labeled <- file.path(out_dir, paste0(id, "_labeled.png"))
      path_plain <- file.path(out_dir, paste0(id, "_plain.png"))
      ggplot2::ggsave(path_labeled, p_labeled, width = width, height = height, dpi = dpi, bg = "white")
      ggplot2::ggsave(path_plain, p_plain, width = width, height = height, dpi = dpi, bg = "white")
      save_log[[i]] <<- data.frame(
        compound_id = id, path_labeled = path_labeled, path_plain = path_plain,
        stringsAsFactors = FALSE
      )
    }

    p <- .radar_ggplot_one(vertices, zone_xy, label, engine, show_labels = TRUE, fig_width = width)
    if (engine == "static") return(p)
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)),
                     width_svg = 5, height_svg = 5)
  }), adme$compound_id)

  if (save) {
    log_df <- do.call(rbind, save_log)
    patliRResults(proj, "admet_radar_log") <- log_df
    .write_results_csv(proj, "admet_radar_log", log_df)
    attr(plots, "proj") <- proj
  }

  plots
}

#' @keywords internal
.radar_zone_xy <- function(angles) {
  norm <- function(x) (x - .radar_axes$axis_min) / (.radar_axes$axis_max - .radar_axes$axis_min)
  r_min <- pmin(pmax(norm(.radar_axes$zone_min), 0), 1)
  r_max <- pmin(pmax(norm(.radar_axes$zone_max), 0), 1)
  list(
    inner = .close_xy(data.frame(x = r_min * sin(angles), y = r_min * cos(angles))),
    outer = .close_xy(data.frame(x = r_max * sin(angles), y = r_max * cos(angles)))
  )
}

#' @keywords internal
.radar_compound_xy <- function(adme_row, angles) {
  raw <- vapply(.radar_axes$property, function(p) adme_row[[p]], numeric(1))
  r <- pmin(pmax((raw - .radar_axes$axis_min) / (.radar_axes$axis_max - .radar_axes$axis_min), 0), 1)
  df <- data.frame(
    axis = .radar_axes$axis, label = .radar_axes$label,
    angle = angles, r = r, raw_value = raw,
    x = r * sin(angles), y = r * cos(angles),
    stringsAsFactors = FALSE
  )
  .close_xy(df)
}

#' Repeat the first row at the end so a polygon/path closes on itself
#' @keywords internal
.close_xy <- function(df) rbind(df, df[1, , drop = FALSE])

#' @keywords internal
.radar_ggplot_one <- function(vertices, zone_xy, label, engine, show_labels = TRUE, fig_width = 6) {
  n_axes <- nrow(.radar_axes)
  outer_ring <- .close_xy(data.frame(x = sin(2 * pi * (seq_len(n_axes) - 1) / n_axes),
                                      y = cos(2 * pi * (seq_len(n_axes) - 1) / n_axes)))
  spokes <- do.call(rbind, lapply(seq_len(n_axes), function(i) {
    data.frame(spoke = i, x = c(0, outer_ring$x[i]), y = c(0, outer_ring$y[i]))
  }))
  label_pos <- data.frame(
    x = 1.35 * sin(2 * pi * (seq_len(n_axes) - 1) / n_axes),
    y = 1.35 * cos(2 * pi * (seq_len(n_axes) - 1) / n_axes),
    label = .radar_axes$label,
    hjust = ifelse(abs(sin(2 * pi * (seq_len(n_axes) - 1) / n_axes)) < 0.1, 0.5,
                    ifelse(sin(2 * pi * (seq_len(n_axes) - 1) / n_axes) > 0, 0, 1)
    )
  )

  tooltip_txt <- sprintf("%s: %.2f", vertices$label, vertices$raw_value)

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = outer_ring, ggplot2::aes(x = .data$x, y = .data$y),
                          fill = NA, colour = "grey80", linewidth = 0.4) +
    ggplot2::geom_line(data = spokes, ggplot2::aes(x = .data$x, y = .data$y, group = .data$spoke),
                       colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_polygon(data = zone_xy$outer, ggplot2::aes(x = .data$x, y = .data$y),
                          fill = "#f4c2c2", alpha = 0.6, colour = NA) +
    ggplot2::geom_polygon(data = zone_xy$inner, ggplot2::aes(x = .data$x, y = .data$y),
                          fill = "white", colour = NA) +
    ggplot2::geom_polygon(data = vertices, ggplot2::aes(x = .data$x, y = .data$y),
                          fill = NA, colour = "#3366cc", linewidth = 1)

  if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p <- p + ggiraph::geom_point_interactive(
      data = vertices, ggplot2::aes(x = .data$x, y = .data$y, tooltip = tooltip_txt),
      colour = "#3366cc", size = 2
    )
  } else {
    p <- p + ggplot2::geom_point(data = vertices, ggplot2::aes(x = .data$x, y = .data$y),
                                  colour = "#3366cc", size = 2)
  }

  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = label_pos, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, hjust = .data$hjust),
      size = 3, lineheight = 0.85
    )
  }

  p +
    ggplot2::coord_fixed(clip = "off", xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
    ggplot2::labs(
      ## long GC-MS compound names would otherwise run past both edges of
      ## the centred title
      title = .plot_wrap(paste0("Bioavailability radar -- ", label), .plot_wrap_width(fig_width, 12, margin = 0.6)),
      subtitle = .plot_wrap("Pink zone = optimal range (Daina, Michielin & Zoete, 2017, Sci. Rep. 7:42717)",
                            .plot_wrap_width(fig_width, 9, margin = 0.6))
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, colour = "grey40"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 20, 10, 20)
    )
}
