#' @include AllGenerics.R internal.R network_bowtie.R plot_network_layers.R
NULL

## plot_bowtie() -- the one plot_* in this batch (2026-07-23) that needs a
## genuinely new Suggests dependency: ggalluvial. A real alluvial/Sankey
## needs curved flow ribbons between stacked strata, which is not
## reasonably hand-rollable with plain ggplot2 geoms the way the radar
## (plot_admet_radar(), manual trig) or the hairball (plot_network_layers(),
## manual igraph layout + geom_segment) were -- unlike circlize/GOplot/
## UpSetR (skipped in this same batch in favor of hand-rolled ggplot2
## equivalents), ggalluvial IS itself a ggplot2 extension (geom_alluvium/
## geom_stratum, no separate plotting system, no heavy transitive deps),
## so it fits the package's existing "ggplot2-first" Suggests philosophy
## rather than fighting it.

#' Alluvial plot: which compounds' targets fall into which bowtie component
#'
#' @description
#' Two-axis alluvial of `patliRResults(proj, "network_bowtie")`: left
#' stratum is `compound_id` (labeled by name), right stratum is
#' `bowtie_component` (`"core"`/`"in_component"`/`"out_component"`/
#' `"other"`/`"unmapped"`), flow thickness is the number of that
#' compound's targets landing in that component for `condition`. Answers
#' "for this extract, which compounds actually reach the core of the
#' directed action network, versus only its periphery (or don't map at
#' all)?" at a glance, instead of reading `network_bowtie()`'s row-per-
#' target table directly.
#'
#' @inheritParams network_build
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
#' \donttest{
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
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_bowtie}.")
  }
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_bowtie} needs the {.pkg ggalluvial} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"ggalluvial\")} -- pure CRAN, ggplot2-based, no heavy transitive deps."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  bowtie_all <- patliRResults(proj, "network_bowtie")
  if (is.null(bowtie_all) || nrow(bowtie_all) == 0) {
    cli::cli_abort(c("No {.val network_bowtie} entry in {.arg proj}.", "i" = "Run {.fn network_bowtie} first."))
  }
  dat <- bowtie_all[bowtie_all$condition %in% conditions, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_bowtie} rows for condition(s) {.val {conditions}}.")
  }
  dat$compound_label <- .network_layers_labels(proj, conditions, data.frame(name = dat$compound_id, layer = "compound", stringsAsFactors = FALSE))

  flow <- as.data.frame(table(compound_label = dat$compound_label, bowtie_component = dat$bowtie_component), stringsAsFactors = FALSE)
  flow <- flow[flow$Freq > 0, , drop = FALSE]
  names(flow)[names(flow) == "Freq"] <- "n_targets"

  if (!is.null(top_n_compounds)) {
    totals <- stats::aggregate(n_targets ~ compound_label, flow, sum)
    keep <- utils::head(totals$compound_label[order(-totals$n_targets)], top_n_compounds)
    flow <- flow[flow$compound_label %in% keep, , drop = FALSE]
  }

  component_colors <- c(
    core = "#c0392b", in_component = "#e67e22", out_component = "#2980b9",
    other = "grey60", unmapped = "grey85"
  )

  p <- ggplot2::ggplot(
    flow,
    ggplot2::aes(axis1 = .data$compound_label, axis2 = .data$bowtie_component, y = .data$n_targets)
  ) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$bowtie_component), width = 1 / 4, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1 / 4, fill = "grey92", colour = "grey40") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, ggplot2::aes(label = ggplot2::after_stat(.data$stratum)), size = 2.8) +
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

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("bowtie_", scope_label, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
    log_row <- data.frame(condition = scope_label, path = path, stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "bowtie_plot_log", log_row, "condition")
    patliRResults(proj, "bowtie_plot_log") <- log_df
    .write_results_csv(proj, "bowtie_plot_log", log_df)
  }

  if (save) attr(p, "proj") <- proj
  p
}
