#' @include AllGenerics.R internal.R network_bowtie.R plot-helpers.R
NULL

## Uses ggalluvial (Suggests): curved flow ribbons between stacked strata
## are not reasonably hand-rollable with plain ggplot2 geoms, and
## ggalluvial is itself a ggplot2 extension (no separate plotting system).

#' Alluvial plot: which compounds' targets fall into which bowtie component
#'
#' @description
#' Two-axis alluvial of `patliRResults(proj, "network_bowtie")`: left
#' stratum is `compound_id` (labeled by name), right stratum is
#' `bowtie_component` (`"core"`/`"in_component"`/`"out_component"`/
#' `"not_in_action_network"`/`"other"`/`"unmapped"`), flow thickness is the
#' number of that compound's targets landing in that component for
#' `condition`. Answers
#' "for this extract, which compounds actually reach the core of the
#' directed action network, versus only its periphery (or don't map at
#' all)?" at a glance, instead of reading `network_bowtie()`'s row-per-
#' target table directly.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param top_n_compounds Integer or `NULL` (default `NULL`, all). If set,
#'   restricts to the `top_n_compounds` with the most targets in scope --
#'   with many compounds the alluvial gets visually crowded.
#'
#' @return A `ggplot` object (`ggalluvial` has no `ggiraph` equivalent, so
#'   this function has no `engine` argument, unlike the rest of the
#'   family). If `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "bowtie_plot_log")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' abundance <- read.csv(
#'   system.file("extdata", "input_abundance_matrix.csv", package = "patliR"),
#'   check.names = FALSE
#' )
#' proj <- prep_binarize(proj, abundance)
#' proj <- targets_import_batch(
#'   proj,
#'   system.file("extdata", "import_targets", package = "patliR"),
#'   platform = "superpred"
#' )
#' proj <- network_build(proj)
#' proj <- network_bowtie(proj, condition = "FLO-ET")
#' plot_bowtie(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_bowtie <- function(proj, condition = NULL, top_n_compounds = NULL,
                         save = TRUE, out_dir = NULL, width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  .plot_require(extra = "ggalluvial")
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  bowtie_all <- patliRResults(proj, "network_bowtie")
  if (is.null(bowtie_all) || nrow(bowtie_all) == 0) {
    cli::cli_abort(c("No {.val network_bowtie} entry in {.arg proj}.", "i" = "Run {.fn network_bowtie} first."))
  }
  dat <- bowtie_all[bowtie_all$condition %in% conditions, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_bowtie} rows for condition(s) {.val {conditions}}.")
  }
  flow <- as.data.frame(table(compound_id = dat$compound_id, bowtie_component = dat$bowtie_component), stringsAsFactors = FALSE)
  flow <- flow[flow$Freq > 0, , drop = FALSE]
  names(flow)[names(flow) == "Freq"] <- "n_targets"

  if (!is.null(top_n_compounds)) {
    totals <- stats::aggregate(n_targets ~ compound_id, flow, sum)
    keep <- utils::head(totals$compound_id[order(-totals$n_targets)], top_n_compounds)
    flow <- flow[flow$compound_id %in% keep, , drop = FALSE]
  }
  flow$compound_label <- .plot_unique_labels(proj, conditions, flow$compound_id, "compound")
  label_lookup <- stats::setNames(flow$compound_label, flow$compound_id)
  stratum_label <- function(x) {
    label <- unname(label_lookup[as.character(x)])
    ifelse(is.na(label), as.character(x), label)
  }

  component_colors <- c(
    core = "#c0392b", in_component = "#e67e22", out_component = "#2980b9",
    not_in_action_network = "#8e44ad", other = "grey60", unmapped = "grey85"
  )

  ## Stratum labels sit OUTSIDE the boxes -- compound names (shortened) to
  ## the left of axis 1, component names to the right of axis 2 -- instead
  ## of centred on a box a quarter of an axis unit wide, which long GC-MS
  ## names overflowed and ran off the figure. The x expansion is sized to
  ## the longest label on each side (.bowtie_x_expansion()).
  side_label <- function(stratum, x, side) {
    ifelse(x == side, .plot_truncate(stratum_label(stratum)), "")
  }
  left_labels <- .plot_truncate(unique(flow$compound_label))
  right_labels <- unique(as.character(flow$bowtie_component))
  x_expand <- .bowtie_x_expansion(max(nchar(left_labels)), max(nchar(right_labels)), width)
  label_offset <- 1 / 8 + 0.03 # half the stratum width plus a gap
  label_layer <- function(side) {
    .plot_text_layer(
      stat = ggalluvial::StatStratum,
      mapping = ggplot2::aes(label = ggplot2::after_stat(side_label(.data$stratum, .data$x, side))),
      size = 2.8, hjust = if (side == 1) 1 else 0,
      nudge_x = if (side == 1) -label_offset else label_offset,
      repel_args = list(direction = "y", min.segment.length = Inf, box.padding = 0.05,
                        point.padding = 0, force = 0.5, max.overlaps = Inf, seed = 1),
      text_args = list()
    )
  }

  p <- ggplot2::ggplot(
    flow,
    ggplot2::aes(axis1 = .data$compound_id, axis2 = .data$bowtie_component, y = .data$n_targets)
  ) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$bowtie_component), width = 1 / 4, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1 / 4, fill = "grey92", colour = "grey40") +
    label_layer(1) + label_layer(2) +
    ggplot2::scale_x_discrete(limits = c("compound", "bowtie_component"),
                              expand = ggplot2::expansion(add = x_expand)) +
    ggplot2::scale_fill_manual(values = component_colors, name = "Bowtie\ncomponent") +
    ggplot2::labs(
      title = paste0("Compound -> bowtie component -- ", scope_label),
      subtitle = .plot_wrap(
        "Flow thickness = number of that compound's targets landing in each component of the STRING directed-action network",
        .plot_wrap_width(width, 8)
      ),
      x = NULL, y = "Number of targets"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_blank()
    )

  .plot_finish(
    proj, p,
    name = "bowtie_plot_log",
    filename = paste0("bowtie_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = NULL, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}

#' Discrete x expansion that leaves room for the outside stratum labels
#'
#' @description
#' `plot_bowtie()` has two axes at x = 1, 2 and draws its labels outside
#' the strata (left of axis 1, right of axis 2), starting `offset` axis
#' units from each axis. With `e_l`/`e_r` the added expansion, the panel
#' spans `S = 1 + e_l + e_r` axis units over `panel_in` inches, and a
#' label `t` inches long fits when `(e - offset) * panel_in / S >= t`.
#' Writing `f = t / panel_in` for each side, the smallest such expansion
#' solves `e_l = offset + f_l * S`, `e_r = offset + f_r * S`, i.e.
#' `S = (1 + 2 * offset) / (1 - f_l - f_r)`. The label fractions are
#' capped at 0.8 together so the flows always keep some width.
#'
#' @param left_chars,right_chars Longest label (characters) on each side.
#' @param fig_width Figure width in inches.
#' @param char_in Average glyph width in inches (size-2.8 text ~0.061).
#' @param reserve_in Figure width not available to the panel (legend,
#'   y-axis title, margins).
#' @param offset Axis units between an axis and the start of its labels.
#' @return Numeric `c(left, right)` for `ggplot2::expansion(add = )`.
#' @keywords internal
.bowtie_x_expansion <- function(left_chars, right_chars, fig_width, char_in = 0.061,
                                reserve_in = 2.3, offset = 1 / 8 + 0.03) {
  panel_in <- max(fig_width - reserve_in, 1)
  f <- c(left_chars, right_chars) * char_in / panel_in
  if (sum(f) > 0.8) f <- f * 0.8 / sum(f)
  s <- (1 + 2 * offset) / (1 - sum(f))
  offset + f * s
}
