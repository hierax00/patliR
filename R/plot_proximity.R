#' @include AllGenerics.R internal.R network_proximity.R plot-helpers.R
NULL

## plot_proximity() -- view = "z" (default) is a z-score forest/lollipop
## plot built only from what network_proximity() persists by default
## (d_random_mean, d_random_sd, z_score), not the n_random raw resampled
## distances themselves -- storing 1000 raw distances per compound/
## condition in a CSV by default would bloat results/ substantially for a
## number nobody downstream needs once the z-score is computed. view =
## "null" is the literal permutation-histogram figure the literature draws
## (Guney et al. 2016 Fig. 1; Menche et al. 2015 SI) -- it needs the raw
## draws, so it reads the opt-in network_proximity(store_null = TRUE)
## output instead, and aborts with a clear message if that was never run.

#' Lollipop plot, or per-compound null-distribution histograms, of
#' `network_proximity()` results
#'
#' @description
#' `view = "z"` (default): one point per (`condition`, `compound_id`) row
#' of `patliRResults(proj, "network_proximity")`, `z_score` on the y-axis,
#' compounds sorted by `z_score`, with reference lines at `z = 0` and the
#' usual \eqn{z = \pm 1.96} two-sided 95% bands. Points are colored by
#' whether they cross that threshold -- a very negative `z_score` means
#' the compound's targets sit significantly *closer* to the disease module
#' than the degree-matched null model, i.e. topologically meaningful
#' proximity, which is what [network_proximity()] is actually testing for.
#' Every compound is named on the axis, and each significant one also
#' carries its `z` value beside the point. Diseases are shown by name
#' (from [disease_genes_fetch()] / [targets_disease_profile()]) with the ID
#' in brackets, in the subtitle (single disease) or the panel strips.
#'
#' `view = "null"`: the literal permutation-histogram figure (Guney et al.
#' 2016 Fig. 1; Menche et al. 2015 SI) -- one facet per `(condition,
#' disease_id, compound_id)` in scope, a histogram of that compound's raw
#' `d_random` null draws from `patliRResults(proj,
#' "network_proximity_null")`, a vertical red line at its `d_observed`
#' (joined from `network_proximity`), and `z_score`/`p_adjusted` (or
#' `p_empirical` when no `p_adjusted` is available) annotated per panel.
#'
#' @section `view = "z"` -- not a permutation histogram:
#' `network_proximity()` only stores the null model's summary
#' (`d_random_mean`, `d_random_sd`) and the resulting `z_score`, not the
#' `n_random` raw resampled distances themselves, *unless* it was run with
#' `store_null = TRUE` -- so by default there is nothing to draw a literal
#' permutation histogram from. `view = "z"` uses what is actually persisted
#' by every run; `view = "null"` uses the opt-in `store_null = TRUE`
#' output (see below).
#'
#' @section `view = "null"` needs `network_proximity(store_null = TRUE)`:
#' This is the one place in the package where a plot's prerequisite is an
#' *optional* analysis output, not a required one. `view = "null"` reads
#' `patliRResults(proj, "network_proximity_null")`; if that slot does not
#' exist, or has no rows for the condition(s)/disease/compounds in scope,
#' it aborts naming `network_proximity(store_null = TRUE)` rather than
#' silently falling back to `view = "z"`. `top_n` (default `12`) caps how
#' many compound facets are drawn, ranked by `|z_score|` descending -- a
#' grid of 30+ histograms is illegible, and a message reports how many
#' compounds were omitted when the scope has more than `top_n`.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character EFO ID, or `NULL` (default) for every disease
#'   present in `patliRResults(proj, "network_proximity")` (faceted).
#' @param view `"z"` (default) -- the z-score lollipop, unchanged from
#'   earlier `patliR` versions. `"null"` -- the per-compound null-draw
#'   histogram figure; needs `network_proximity(store_null = TRUE)` (see
#'   Details).
#' @param top_n Integer, default `12`. `view = "null"` only: maximum
#'   number of compound facets to draw, ranked by `|z_score|` descending.
#'   Ignored for `view = "z"`.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it -- `view = "z"` to
#'   `patliRResults(proj, "proximity_plot_log")` (filename
#'   `proximity_<scope>.png`); `view = "null"` to its own slot,
#'   `patliRResults(proj, "proximity_null_plot_log")` (filename
#'   `proximity_null_<scope>.png`), so the two views never overwrite each
#'   other's log row or file.
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
#' proj <- targets_disease_filter(proj, disease = "type 2 diabetes mellitus")
#' proj <- network_proximity(
#'   proj, condition = "FLO-ET",
#'   disease = unique(patliRResults(proj, "targets_disease")$disease_id)[1]
#' )
#' plot_proximity(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @references
#' Guney, Menche, Vidal & Barabasi (2016), *Nat Commun* 7:10331, Fig. 1,
#' \doi{10.1038/ncomms10331}. Menche et al. (2015), *Science*
#' 347(6224):1257601, Supplementary Information.
#'
#' @export
plot_proximity <- function(proj, condition = NULL, disease = NULL,
                            view = c("z", "null"), top_n = 12,
                            engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                            width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  view <- match.arg(view)
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  if (view == "null") {
    return(.plot_proximity_null(
      proj, conditions, scope_label, disease, top_n,
      engine, save, out_dir, width, height, dpi
    ))
  }

  prox_all <- patliRResults(proj, "network_proximity")
  if (is.null(prox_all) || nrow(prox_all) == 0) {
    cli::cli_abort(c("No {.val network_proximity} entry in {.arg proj}.", "i" = "Run {.fn network_proximity} first."))
  }
  dat <- prox_all[prox_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) dat <- dat[dat$disease_id %in% disease, , drop = FALSE]
  dat <- dat[!is.na(dat$z_score), , drop = FALSE]

  ## network_proximity() records which disease-gene set each row was built
  ## from. If a results table mixes the independent "disease_genes" set with
  ## the legacy circular "targets_disease" set, the circular rows would
  ## render alongside -- and typically as the most significant points, which
  ## is exactly what the circularity produces. Drop them, loudly.
  if ("disease_gene_source" %in% names(dat) &&
      all(c("disease_genes", "targets_disease") %in% dat$disease_gene_source)) {
    n_drop <- sum(dat$disease_gene_source == "targets_disease")
    cli::cli_warn(c(
      "!" = "{.val network_proximity} mixes both disease-gene sources; dropping {n_drop} circular {.val targets_disease} row{?s} from the plot.",
      "i" = "Re-run {.fn network_proximity} with the default {.code disease_genes = \"disease_genes\"} to replace them."
    ))
    dat <- dat[dat$disease_gene_source != "targets_disease", , drop = FALSE]
  }
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_proximity} rows with a non-NA z_score for the requested condition(s)/disease.")
  }
  dat$label <- .plot_unique_labels(proj, conditions, dat$compound_id, "compound")
  dat <- dat[order(dat$z_score), , drop = FALSE]
  dat$compound_id <- factor(dat$compound_id, levels = unique(dat$compound_id))
  label_lookup <- stats::setNames(dat$label, dat$compound_id)
  dat$significant <- abs(dat$z_score) >= 1.96
  ## disease name, not only its ID, wherever the disease is named
  dat$disease_label <- .plot_disease_label(proj, dat$disease_id)
  dat$panel <- paste(dat$condition, dat$disease_label, sep = " / ")
  dat$z_text <- sprintf("%.2f", dat$z_score)
  dat$tooltip <- sprintf("%s\nz = %.2f (d_obs = %.2f, null = %.2f +/- %.2f)", dat$label, dat$z_score, dat$d_observed, dat$d_random_mean, dat$d_random_sd)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$compound_id, y = .data$z_score, colour = .data$significant)) +
    ## reference lines as data, not a length-3 linetype vector: a vector
    ## aesthetic is recycled per panel and fails once facet_wrap() gives more
    ## than one (condition, disease) panel
    ggplot2::geom_hline(
      data = data.frame(yintercept = c(-1.96, 0, 1.96), lt = c("22", "solid", "22")),
      ggplot2::aes(yintercept = .data$yintercept, linetype = .data$lt), colour = "grey60", inherit.aes = FALSE
    ) +
    ggplot2::scale_linetype_identity() +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$compound_id, y = 0, yend = .data$z_score), linewidth = 0.4, alpha = 0.6)
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$compound_id), size = 2.5)
  } else {
    ggplot2::geom_point(size = 2.5)
  }
  ## every significant compound is named on the axis; its z value is
  ## printed beside the point, on the outer side of the lollipop
  sig <- dat[dat$significant, , drop = FALSE]
  if (nrow(sig) > 0) {
    sig$hjust <- ifelse(sig$z_score < 0, 1.35, -0.35)
    p <- p + ggplot2::geom_text(
      data = sig, ggplot2::aes(x = .data$compound_id, y = .data$z_score, label = .data$z_text, hjust = .data$hjust),
      size = 2.5, colour = "#922b21", inherit.aes = FALSE, show.legend = FALSE
    )
    ## room for the printed value beyond the most extreme points
    span <- diff(range(c(dat$z_score, -1.96, 1.96)))
    p <- p + ggplot2::expand_limits(y = c(min(dat$z_score) - 0.2 * span, max(dat$z_score, 1.96) + 0.2 * span * any(sig$z_score > 0)))
  }
  diseases <- unique(dat$disease_label)
  p <- p +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(labels = label_lookup) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "grey40"), name = "|z| >= 1.96") +
    ggplot2::labs(
      title = paste0("Network proximity z-scores -- ", scope_label),
      subtitle = .plot_wrap(
        paste0(
          if (length(diseases) == 1) paste0("Disease: ", diseases, "\n") else "",
          "More negative = targets significantly closer to the disease module than the degree-matched null model; ",
          "numbers = z of the significant compounds"
        ),
        .plot_wrap_width(width, 8)
      ),
      x = NULL, y = "z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot"
    )
  ## Several (condition, disease) panels side by side sharing the compound
  ## axis: facet_wrap(scales = "free_y") repeated the full compound labels
  ## per panel and squeezed each panel to a sliver with clipped strips.
  if (length(unique(dat$panel)) > 1) {
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(.data$panel), labeller = ggplot2::label_wrap_gen(25))
  }

  .plot_finish(
    proj, p,
    name = "proximity_plot_log",
    filename = paste0("proximity_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}

#' `plot_proximity(view = "null")` -- one histogram facet per compound of
#' its raw `network_proximity_null` draws, `d_observed` overlaid
#'
#' @description
#' The literal permutation-histogram figure (Guney et al. 2016 Fig. 1;
#' Menche et al. 2015 SI), built from the opt-in
#' `network_proximity(store_null = TRUE)` output -- see
#' [plot_proximity()]'s `@section`s for why this cannot be drawn from the
#' default `network_proximity()` output.
#' @keywords internal
.plot_proximity_null <- function(proj, conditions, scope_label, disease, top_n,
                                  engine, save, out_dir, width, height, dpi) {
  null_all <- patliRResults(proj, "network_proximity_null")
  if (is.null(null_all) || nrow(null_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_proximity_null} entry in {.arg proj}.",
      "i" = "{.fn plot_proximity}'s {.code view = \"null\"} needs the raw null draws, which {.fn network_proximity} only stores when asked to.",
      "i" = "Re-run {.fn network_proximity} with {.code store_null = TRUE} for this condition/disease first."
    ))
  }
  dat_null <- null_all[null_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) dat_null <- dat_null[dat_null$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat_null) == 0) {
    cli::cli_abort(c(
      "No {.val network_proximity_null} rows for the requested condition(s)/disease.",
      "i" = "Re-run {.fn network_proximity} with {.code store_null = TRUE} for this scope."
    ))
  }

  prox_all <- patliRResults(proj, "network_proximity")
  if (is.null(prox_all) || nrow(prox_all) == 0) {
    cli::cli_abort(c("No {.val network_proximity} entry in {.arg proj}.", "i" = "Run {.fn network_proximity} first."))
  }
  main <- prox_all[prox_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) main <- main[main$disease_id %in% disease, , drop = FALSE]
  main <- main[!is.na(main$z_score), , drop = FALSE]

  ## Only compounds present in BOTH tables can be drawn -- d_observed comes
  ## from network_proximity, the histogram from network_proximity_null.
  key_null <- paste(dat_null$condition, dat_null$disease_id, dat_null$compound_id, sep = "\r")
  key_main <- paste(main$condition, main$disease_id, main$compound_id, sep = "\r")
  main_scope <- main[key_main %in% key_null, , drop = FALSE]
  if (nrow(main_scope) == 0) {
    cli::cli_abort(c(
      "No compound in scope has both a {.val network_proximity} row and stored {.val network_proximity_null} draws.",
      "i" = "Re-run {.fn network_proximity} with {.code store_null = TRUE} for this condition/disease."
    ))
  }

  ## top_n cap, ranked by |z_score| descending -- a grid of 30+ histograms
  ## is illegible (spec 3.7).
  main_scope <- main_scope[order(-abs(main_scope$z_score)), , drop = FALSE]
  n_in_scope <- nrow(main_scope)
  if (n_in_scope > top_n) {
    n_dropped <- n_in_scope - top_n
    cli::cli_inform(c(
      "i" = "{n_dropped} of {n_in_scope} compound{?s} omitted from the panel grid ({.code top_n = {top_n}}, ranked by {.code |z_score|})."
    ))
    main_scope <- utils::head(main_scope, top_n)
  }

  keep_key <- paste(main_scope$condition, main_scope$disease_id, main_scope$compound_id, sep = "\r")
  dat_null <- dat_null[paste(dat_null$condition, dat_null$disease_id, dat_null$compound_id, sep = "\r") %in% keep_key, , drop = FALSE]

  main_scope$label <- .plot_unique_labels(proj, conditions, main_scope$compound_id, "compound")
  main_scope$facet <- as.character(seq_len(nrow(main_scope)))
  dat_null$facet <- main_scope$facet[match(
    paste(dat_null$condition, dat_null$disease_id, dat_null$compound_id, sep = "\r"), keep_key
  )]
  main_scope$disease_label <- .plot_disease_label(proj, main_scope$disease_id, with_id = FALSE)
  diseases <- unique(.plot_disease_label(proj, main_scope$disease_id))
  facet_labels <- stats::setNames(
    paste0(.plot_truncate(main_scope$label, 32), "\n(", main_scope$condition, " / ",
           .plot_truncate(main_scope$disease_label, 28), ")"), main_scope$facet
  )

  ## p_adjusted may be all-NA (or absent, on a legacy/fabricated table) --
  ## fall back to p_empirical, same defensive pattern as plot_synergy().
  p_col <- if ("p_adjusted" %in% names(main_scope) && any(!is.na(main_scope$p_adjusted))) "p_adjusted" else "p_empirical"
  main_scope$annot <- sprintf("z = %.2f\np%s = %.3g", main_scope$z_score, if (p_col == "p_adjusted") "_adj" else "_emp", main_scope[[p_col]])

  p <- ggplot2::ggplot(dat_null[is.finite(dat_null$d_random), , drop = FALSE], ggplot2::aes(x = .data$d_random)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2980b9", colour = "white", alpha = 0.85) +
    ggplot2::geom_vline(
      data = main_scope, ggplot2::aes(xintercept = .data$d_observed),
      colour = "#c0392b", linewidth = 0.7, inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = main_scope, ggplot2::aes(x = Inf, y = Inf, label = .data$annot),
      hjust = 1.05, vjust = 1.2, size = 2.6, colour = "grey20", inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~facet, scales = "free", labeller = ggplot2::as_labeller(facet_labels)) +
    ggplot2::labs(
      title = paste0("Network proximity null distributions -- ", scope_label),
      subtitle = .plot_wrap(
        paste0(
          if (length(diseases) == 1) paste0("Disease: ", diseases, "\n") else "",
          "Degree-preserving null resamples (network_proximity(store_null = TRUE)); red line = d_observed (Guney et al. 2016 Fig. 1)"
        ),
        .plot_wrap_width(width, 8)
      ),
      x = "d_random (null draw)", y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      strip.text = ggplot2::element_text(size = 7.5)
    )

  .plot_finish(
    proj, p,
    name = "proximity_null_plot_log",
    filename = paste0("proximity_null_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}
