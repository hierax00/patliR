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
#'   by its `compound_id`) and `boiled_egg_plain.png` (same scatter, no
#'   point labels -- useful once you have many compounds and the labels
#'   start overlapping).
#' @param out_dir Directory to write the PNG files to (only used if
#'   `save = TRUE`). Defaults to `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] for the saved PNG.
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
                             width = 7, height = 5, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)

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
  adme$call <- ifelse(
    !is.na(adme$bbb_permeant) & adme$bbb_permeant, "BBB permeant",
    ifelse(!is.na(adme$gi_absorption) & adme$gi_absorption == "High", "GI absorption (high)", "Neither")
  )

  poly <- .load_boiled_egg_polygons()

  p <- .boiled_egg_ggplot(adme, poly, engine = engine, show_labels = FALSE)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    p_labeled <- .boiled_egg_ggplot(adme, poly, engine = "static", show_labels = TRUE)
    p_plain <- .boiled_egg_ggplot(adme, poly, engine = "static", show_labels = FALSE)
    path_labeled <- file.path(out_dir, "boiled_egg_labeled.png")
    path_plain <- file.path(out_dir, "boiled_egg_plain.png")
    ggplot2::ggsave(path_labeled, p_labeled, width = width, height = height, dpi = dpi)
    ggplot2::ggsave(path_plain, p_plain, width = width, height = height, dpi = dpi)
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

#' @keywords internal
.boiled_egg_ggplot <- function(adme, poly, engine, show_labels = FALSE) {
  tooltip_txt <- sprintf(
    "%s\nTPSA: %.1f, WLogP (proxy): %.2f\n%s",
    adme$compound_label, adme$tpsa, adme$wlogp_proxy, adme$call
  )

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = -20, xmax = 220, ymin = -3, ymax = 8, fill = "grey96") +
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
    ## Plain geom_text, no ggrepel: patliR keeps Suggests minimal, so labels
    ## can overlap with many compounds close together -- use the "plain"
    ## variant for that case instead of fighting label placement here.
    p <- p + ggplot2::geom_text(
      data = adme,
      ggplot2::aes(x = .data$tpsa, y = .data$wlogp_proxy, label = .data$compound_id),
      size = 2.8, vjust = -0.9, colour = "grey20", show.legend = FALSE
    )
  }

  p +
    ggplot2::scale_colour_manual(values = c(
      "BBB permeant" = "#8e44ad", "GI absorption (high)" = "#2980b9", "Neither" = "grey30"
    )) +
    ggplot2::coord_cartesian(xlim = c(-20, 220), ylim = c(-3, 8)) +
    ggplot2::labs(
      x = expression("TPSA (" * ring(A)^2 * ")"), y = "WLogP (CDK ALogP proxy)",
      colour = NULL,
      title = "BOILED-Egg: predicted GI absorption and BBB penetration",
      subtitle = "White = high GI absorption; yolk = high BBB penetration (Daina & Zoete, 2016)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey40"),
      legend.position = "bottom"
    )
}
