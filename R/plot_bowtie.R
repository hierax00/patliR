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

  p <- ggplot2::ggplot(
    flow,
    ggplot2::aes(axis1 = .data$compound_id, axis2 = .data$bowtie_component, y = .data$n_targets)
  ) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$bowtie_component), width = 1 / 4, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1 / 4, fill = "grey92", colour = "grey40") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, ggplot2::aes(label = stratum_label(ggplot2::after_stat(.data$stratum))), size = 2.8) +
    ggplot2::scale_x_discrete(limits = c("compound", "bowtie_component"), expand = c(0.1, 0.1)) +
    ggplot2::scale_fill_manual(values = component_colors, name = "Bowtie\ncomponent") +
    ggplot2::labs(
      title = paste0("Compound -> bowtie component -- ", scope_label),
      subtitle = "Flow thickness = number of that compound's targets landing in each component of the STRING directed-action network",
      x = NULL, y = "Number of targets"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
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
