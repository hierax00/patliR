#' @include AllGenerics.R internal.R
NULL

#' BOILED-Egg plot: predicted GI absorption and BBB penetration
#'
#' @description
#' Scatter of every compound's TPSA (x) against a WLogP proxy (y), with the
#' two published BOILED-Egg regions drawn to scale: the white "egg white"
#' ellipse (high probability of passive gastrointestinal absorption) and
#' the yellow "yolk" ellipse (high probability of BBB penetration).
#' Compounds falling outside both are flagged neither. Needs [adme_local()]
#' to have been run first.
#'
#' @section Provenance of the ellipse boundaries:
#' The two ellipses are the published model's own numeric parameters,
#' from the supporting information of Daina, A. & Zoete, V. (2016), "A
#' BOILED-Egg To Predict Gastrointestinal Absorption and Brain
#' Penetration of Small Molecules", *ChemMedChem* 11, 1117-1121,
#' \doi{10.1002/cmdc.201600182} -- copied directly from the hard-coded
#' coordinate lists in the GPL-3 reference implementation
#' <https://github.com/bfmilne/PyBOILEDegg> (Milne, B.F., 2021,
#' \doi{10.5281/zenodo.4725530}), whose own comment states the same
#' original source. Not produced by running that program -- it only ever
#' outputs a classification, never boundary coordinates. These are the
#' real published boundaries, not an approximation.
#'
#' @section The one real caveat -- WLogP:
#' The original model's y-axis is **WLogP** (Wildman & Crippen's
#' atom-contribution logP, RDKit's implementation). `rcdk`/CDK has no
#' identical descriptor, so this plot (and `adme_local()`'s
#' `gi_absorption`/`bbb_permeant` columns) use CDK's ALogP
#' (`rcdk::get.alogp()`) as the closest available proxy -- a different,
#' but methodologically related, atom-contribution method. A compound near
#' an ellipse boundary could land on a different side of the line than it
#' would in the real SwissADME tool.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound with [adme_local()] results.
#' @param engine `"ggiraph"` (default, if installed): interactive plot with
#'   a tooltip per point. `"static"`: plain `ggplot2`.
#' @param save Logical, default `TRUE`. If `TRUE`, also writes two static
#'   PNGs to `out_dir` -- `boiled_egg_labeled.png` (with each point labeled
#'   by `compound_id` or name, see `label`) and `boiled_egg_plain.png` (same scatter, no
#'   point labels -- useful once you have many compounds and the labels
#'   start overlapping).
#' @param out_dir Directory to write the PNG files to (only used if
#'   `save = TRUE`). Defaults to `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] for the saved PNG.
#' @param label What the labeled PNG (`boiled_egg_labeled.png`) prints
#'   next to each point: `"id"` (default, `compound_id`), `"name"` (the
#'   compound name, shortened to 24 characters; `compound_id` when there is
#'   none) or `"none"`. Labels are repelled when the optional `ggrepel`
#'   package is installed, and ones that still cannot be placed without
#'   overlapping are left out (the plain PNG and the `ggiraph` tooltips
#'   cover every compound).
#'
#' @section Axis range:
#' The plot opens on the published frame (TPSA -20 to 220, WLogP -3 to 8)
#' and widens it to include every compound -- lipophilic volatiles often
#' exceed WLogP 8 -- noting the widening in a caption. The two ellipses are
#' drawn from their published coordinates regardless.
#'
#' @return A `girafe` htmlwidget (`engine = "ggiraph"`) or a `ggplot`
#'   object (`engine = "static"`) -- one plot for all requested compounds
#'   together (this is a comparative overview, unlike [plot_admet_radar()]
#'   which is one compound at a time). If `save = TRUE` (the default), the
#'   PNG path is also recorded in `patliRResults(proj, "boiled_egg_log")`,
#'   retrievable via `attr(result, "proj")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' plot_boiled_egg(proj, engine = "static", save = FALSE)
#' }
#'
#' @export
plot_boiled_egg <- function(proj, compound_ids = NULL, engine = c("ggiraph", "static"),
                             save = TRUE, out_dir = NULL,
                             width = 7, height = 5, dpi = 150, label = c("id", "name", "none")) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  label <- match.arg(label)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_boiled_egg}.")
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

  cmp <- compounds(proj)
  name_lookup <- stats::setNames(cmp$name, cmp$id)
  adme$compound_label <- paste0(
    adme$compound_id,
    ifelse(is.na(name_lookup[adme$compound_id]), "", paste0(" (", name_lookup[adme$compound_id], ")"))
  )
  adme$point_label <- switch(label,
    id = adme$compound_id,
    name = ifelse(is.na(name_lookup[adme$compound_id]) | name_lookup[adme$compound_id] == "",
                  adme$compound_id, .plot_truncate(name_lookup[adme$compound_id], 24)),
    none = NA_character_
  )
  adme$call <- ifelse(
    !is.na(adme$bbb_permeant) & adme$bbb_permeant, "BBB permeant",
    ifelse(!is.na(adme$gi_absorption) & adme$gi_absorption == "High", "GI absorption (high)", "Neither")
  )

  poly <- .load_boiled_egg_polygons()

  p <- .boiled_egg_ggplot(adme, poly, engine = engine, show_labels = FALSE)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    p_labeled <- .boiled_egg_ggplot(adme, poly, engine = "static", show_labels = label != "none")
    p_plain <- .boiled_egg_ggplot(adme, poly, engine = "static", show_labels = FALSE)
    path_labeled <- file.path(out_dir, "boiled_egg_labeled.png")
    path_plain <- file.path(out_dir, "boiled_egg_plain.png")
    ggplot2::ggsave(path_labeled, p_labeled, width = width, height = height, dpi = dpi, bg = "white")
    ggplot2::ggsave(path_plain, p_plain, width = width, height = height, dpi = dpi, bg = "white")
    log_df <- data.frame(
      path_labeled = path_labeled, path_plain = path_plain, n_compounds = nrow(adme),
      stringsAsFactors = FALSE
    )
    patliRResults(proj, "boiled_egg_log") <- log_df
    .write_results_csv(proj, "boiled_egg_log", log_df)
  }

  result <- if (engine == "static") {
    p
  } else {
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  }
  if (save) attr(result, "proj") <- proj
  result
}

#' Axis limits for the BOILED-Egg plot
#'
#' @description
#' The published figure frame (TPSA -20..220, WLogP -3..8), widened -- with
#' a 4% margin on the widened side -- to take in every finite point:
#' lipophilic GC-MS volatiles routinely have a WLogP proxy above 8 and used
#' to be cut off. The ellipses are drawn from their published coordinates
#' either way, so they keep their shape and position.
#' @return `list(x = c(min, max), y = c(min, max))`.
#' @keywords internal
.boiled_egg_limits <- function(tpsa, wlogp) {
  widen <- function(default, v) {
    v <- v[is.finite(v)]
    lim <- range(c(default, v))
    pad <- 0.04 * diff(lim)
    c(if (lim[1] < default[1]) lim[1] - pad else lim[1], if (lim[2] > default[2]) lim[2] + pad else lim[2])
  }
  list(x = widen(c(-20, 220), tpsa), y = widen(c(-3, 8), wlogp))
}

#' @keywords internal
.boiled_egg_ggplot <- function(adme, poly, engine, show_labels = FALSE) {
  lim <- .boiled_egg_limits(adme$tpsa, adme$wlogp_proxy)
  if (is.null(adme$point_label)) adme$point_label <- adme$compound_id
  tooltip_txt <- sprintf(
    "%s\nTPSA: %.1f, WLogP (proxy): %.2f\n%s",
    adme$compound_label, adme$tpsa, adme$wlogp_proxy, adme$call
  )

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = lim$x[1], xmax = lim$x[2], ymin = lim$y[1], ymax = lim$y[2], fill = "grey96") +
    ggplot2::geom_polygon(data = poly$gia, ggplot2::aes(x = .data$tpsa, y = .data$wlogp),
                          fill = "white", colour = "black", linewidth = 0.6) +
    ggplot2::geom_polygon(data = poly$bbb, ggplot2::aes(x = .data$tpsa, y = .data$wlogp),
                          fill = "#f5d61e", colour = "black", alpha = 0.85, linewidth = 0.6)

  if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p <- p + ggiraph::geom_point_interactive(
      data = adme,
      ggplot2::aes(x = .data$tpsa, y = .data$wlogp_proxy, colour = .data$call, tooltip = tooltip_txt),
      size = 2.5
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = adme, ggplot2::aes(x = .data$tpsa, y = .data$wlogp_proxy, colour = .data$call),
      size = 2.5
    )
  }

  if (show_labels) {
    ## Repelled when the optional ggrepel is installed; labels that still
    ## cannot be placed without overlapping (dozens of compounds stacked at
    ## TPSA 0) are left out rather than printed on top of each other --
    ## the plain PNG and the ggiraph tooltips cover every point.
    p <- p + .plot_text_layer(
      data = adme[!is.na(adme$point_label), , drop = FALSE],
      mapping = ggplot2::aes(x = .data$tpsa, y = .data$wlogp_proxy, label = .data$point_label),
      size = 2.6, colour = "grey20", show.legend = FALSE,
      repel_args = list(box.padding = 0.2, point.padding = 0.1, min.segment.length = 0.3,
                        segment.colour = "grey60", max.overlaps = 15, seed = 1),
      text_args = list(vjust = -0.9)
    )
  }

  p +
    ggplot2::scale_colour_manual(values = c(
      "BBB permeant" = "#8e44ad", "GI absorption (high)" = "#2980b9", "Neither" = "grey30"
    )) +
    ggplot2::coord_cartesian(xlim = lim$x, ylim = lim$y) +
    ggplot2::labs(
      x = expression("TPSA (" * ring(A)^2 * ")"), y = "WLogP (CDK ALogP proxy)",
      colour = NULL,
      title = "BOILED-Egg: predicted GI absorption and BBB penetration",
      subtitle = "White = high GI absorption; yolk = high BBB penetration (Daina & Zoete, 2016)",
      caption = if (any(adme$wlogp_proxy > 8 | adme$tpsa > 220 | adme$wlogp_proxy < -3 | adme$tpsa < -20, na.rm = TRUE))
        "Axes widened beyond the published frame (TPSA -20..220, WLogP -3..8) to show every compound." else NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      legend.position = "bottom"
    )
}
