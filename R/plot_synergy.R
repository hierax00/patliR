#' @include AllGenerics.R internal.R network_synergy.R plot-helpers.R
NULL

#' Cheng et al. (2019) Complementary-Exposure quadrant for `network_synergy()`
#'
#' @description
#' One panel per `(condition, disease_id)`. Each point is a compound pair
#' from `patliRResults(proj, "network_synergy")`, plotted on the
#' `z_score_a` / `z_score_b` plane of Cheng, Kovacs & Barabasi (2019)
#' Fig. 2 -- the quantitative version of their six schematic drug-pair
#' topologies. Every visual channel maps onto one column of
#' [network_synergy()]'s output:
#'
#' - **x / y** -- `z_score_a` / `z_score_b`, but re-ordered *for this plot
#'   only* so the more disease-proximal compound of the pair (the smaller,
#'   more negative [network_proximity()] `z_score`) always lands on the
#'   x-axis. This makes the panel upper-triangular (every point has
#'   `x <= y`) and roughly halves the visual noise from an otherwise
#'   arbitrary `a`/`b` labeling. It does **not** change which compound is
#'   `compound_a` vs. `compound_b` in `patliRResults(proj,
#'   "network_synergy")` -- only how this one plot's two axes are filled
#'   in; rows where either `z_score` is `NA` are left unswapped.
#' - **reference lines and shading** -- `x = 0` / `y = 0` mark
#'   [network_proximity()]'s proximity threshold; the lower-left quadrant
#'   (both compounds individually proximal by raw `z < 0`) is shaded
#'   faintly, since that is the region every `P1`/`P2` pair lives in.
#' - **colour** -- `separated` (`s_AB >= 0`, Menche et al. (2015)'s sign
#'   convention), a two-level **manual** scale, deliberately not a
#'   gradient -- this is the categorical variable Cheng et al.'s whole
#'   result turns on: whether the two compounds' target modules sit in
#'   separate neighbourhoods of the interactome. `NA` (`separation =
#'   "jaccard"`, an unmapped/singleton target set, or a disconnected pair)
#'   is a distinct grey, drawn, never dropped.
#' - **shape** -- `cheng_class` (`"P1"`..`"P6"`, plus `NA`), a named shape
#'   scale. A second facet is deliberately not used for this -- the panel
#'   grid already carries `(condition, disease_id)`.
#' - **size** -- `abs(s_ab)` magnitude (how strongly separated or
#'   overlapping); `NA` maps to the smallest size class rather than
#'   dropping the point (a continuous `ggplot2` scale would otherwise omit
#'   `NA` rows with a bare warning).
#' - **in-plot annotation** -- the `P2` (Complementary Exposure) region
#'   is labeled with its literal name and the count of `P2` pairs in that
#'   panel (`cheng_class == "P2"`, i.e. `complementary_exposure`), since
#'   that is Cheng et al.'s only class shown to correlate with therapeutic
#'   efficacy.
#' - **text labels** -- the `top_n` `P2` pairs with the highest
#'   `synergy_score` (descending, `NA` last -- `synergy_score` is `NA` for
#'   singleton-flagged `P2` pairs even though they are still classified
#'   `P2`) are labeled by compound name, via [.plot_label_nodes()].
#'
#' @section Requires `separation = "network"` results to show anything but grey:
#' `cheng_class` / `s_ab` / `separated` are only non-`NA` when
#' [network_synergy()] was run with `separation = "network"` (the default
#' since patliR 0.2.0) and the pair's target sets mapped onto the STRING
#' LCC and were not singleton-flagged into a disconnected component. A
#' table built entirely with `separation = "jaccard"`, or a legacy table
#' that predates these columns, still draws the `z_score_a` x
#' `z_score_b` scatter (colour grey, shape `NA`), but with no `P2` region
#' to speak of; if **no** row anywhere in scope carries a `cheng_class`,
#' this function draws a self-explanatory empty panel instead of a
#' confusing all-grey one.
#'
#' @section Companion panel -- `s_AB` distribution:
#' Also produces (and, when `save = TRUE`, separately saves and logs) a
#' second figure: a histogram of `s_ab` across every pair in scope, with a
#' vertical line at `s_AB = 0` -- the direct analogue of Cheng et al.'s
#' Fig. 1 separation distribution, and the thing that tells the reader
#' whether *any* pair in this extract is topologically separated at all.
#' Logged to `patliRResults(proj, "synergy_sab_plot_log")` and saved as
#' `synergy_sab_<scope>.png` alongside the main `synergy_<scope>.png`.
#' Only the main scatter is returned by this function; call
#' `attr(plot_synergy(...), "proj")` and look up `synergy_sab_plot_log`
#' for the companion file's path.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character EFO ID, or `NULL` (default) for every disease
#'   present (faceted).
#' @param top_n Integer, default `5`. Number of highest-`synergy_score`
#'   `P2` (Complementary Exposure) pairs to label by compound name, per
#'   panel. **Breaking change from patliR <= 0.1.x**: the old default was
#'   `10`, chosen for the previous complementarity/joint_closeness scatter;
#'   `5` keeps the new, much more selective `P2`-only label set legible.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget (the main scatter). If
#'   `save = TRUE` (default), also writes two PNGs (`synergy_<scope>.png`
#'   and `synergy_sab_<scope>.png`) and logs them to
#'   `patliRResults(proj, "synergy_plot_log")` /
#'   `patliRResults(proj, "synergy_sab_plot_log")` respectively.
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
#' proj <- disease_genes_fetch(proj, disease = "type 2 diabetes mellitus")
#' disease_id <- unique(patliRResults(proj, "disease_genes")$disease_id)[1]
#' proj <- network_proximity(proj, condition = "FLO-ET", disease = disease_id)
#' proj <- network_synergy(proj, condition = "FLO-ET", disease = disease_id)
#' plot_synergy(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @references
#' Cheng, Kovacs & Barabasi (2019), *Nat Commun* 10:1197, Fig. 2,
#' \doi{10.1038/s41467-019-09186-x}. Menche et al. (2015), *Science*
#' 347(6224):1257601, \doi{10.1126/science.1257601}.
#'
#' @export
plot_synergy <- function(proj, condition = NULL, disease = NULL, top_n = 5,
                          engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                          width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 0)
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  syn_all <- patliRResults(proj, "network_synergy")
  if (is.null(syn_all) || nrow(syn_all) == 0) {
    cli::cli_abort(c("No {.val network_synergy} entry in {.arg proj}.", "i" = "Run {.fn network_synergy} first."))
  }
  dat <- syn_all[syn_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(disease)) dat <- dat[dat$disease_id %in% disease, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_synergy} rows for the requested condition(s)/disease.")
  }

  ## Defensive column back-fill: a legacy/fabricated network_synergy table
  ## (pre-Menche s_AB, or hand-built in a test) may not carry these columns
  ## at all -- .network_upsert() back-fills NA for a real rerun, but a
  ## table assigned directly to patliRResults() bypasses that. Treat an
  ## absent column exactly like an all-NA one.
  for (col in c("s_ab", "separated", "cheng_class", "complementary_exposure",
                "synergy_score", "p_adjusted_a", "p_adjusted_b")) {
    if (is.null(dat[[col]])) {
      dat[[col]] <- if (col == "cheng_class") NA_character_ else if (col == "separated" || col == "complementary_exposure") NA else NA_real_
    }
  }

  dat$panel <- paste(dat$condition, dat$disease_id, sep = " / ")

  ## Companion panel: s_AB distribution across every pair in scope.
  n_sab <- sum(!is.na(dat$s_ab))
  p_sab <- if (n_sab == 0) {
    .synergy_empty_panel(
      paste0("No s_AB values available for ", scope_label, " (separation = \"jaccard\", or no pair mapped onto the STRING LCC)."),
      paste0("s_AB distribution -- ", scope_label)
    )
  } else {
    ps <- ggplot2::ggplot(dat[!is.na(dat$s_ab), , drop = FALSE], ggplot2::aes(x = .data$s_ab)) +
      ggplot2::geom_histogram(bins = 30, fill = "#2980b9", colour = "white", alpha = 0.85) +
      ggplot2::geom_vline(xintercept = 0, colour = "#c0392b", linetype = "solid", linewidth = 0.6) +
      ggplot2::labs(
        title = paste0("s_AB distribution -- ", scope_label),
        subtitle = "Menche et al. (2015) network separation across every scored pair; s_AB >= 0 (right of the red line) = topologically separated",
        x = "s_AB", y = "Pair count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40"))
    if (length(unique(dat$panel)) > 1) ps <- ps + ggplot2::facet_wrap(~panel)
    ps
  }

  ## This plot needs a Cheng classification, not a synergy_score -- a pair
  ## can carry a real cheng_class (P1..P6) without ever getting a
  ## synergy_score (that column is gated on complementary_exposure AND not
  ## singleton-flagged). So the empty-state trigger is "no row anywhere in
  ## scope has a non-NA cheng_class", not "no row is scored" -- the latter
  ## is the common, uninteresting case (most pairs are not P2) and would
  ## make this guard fire on almost every real run.
  n_class <- sum(!is.na(dat$cheng_class))
  if (n_class == 0) {
    n_sab <- sum(!is.na(dat$s_ab))
    msg <- paste0(
      "No compound pair in ", scope_label, " has a computed Cheng classification (cheng_class).\n\n",
      "This plot needs network_synergy(separation = \"network\") results (s_ab / separated /\n",
      "cheng_class); separation = \"jaccard\" leaves those columns NA, as does a pair whose\n",
      "target sets did not map onto the STRING LCC.\n\n",
      "s_ab non-NA count: ", n_sab, " of ", nrow(dat), "."
    )
    p_main <- .synergy_empty_panel(msg, paste0("Compound pair Cheng classification -- ", scope_label))
    return(.synergy_finish_both(proj, p_main, p_sab, scope_label, engine, save, out_dir, width, height, dpi))
  }
  if (n_class < nrow(dat)) {
    cli::cli_inform(c(
      "i" = "{nrow(dat) - n_class} of {nrow(dat)} pair{?s} {?has/have} no Cheng classification (cheng_class {.val NA}) -- drawn in grey with no shape assigned."
    ))
  }

  ## Row-wise relabel for THIS PLOT ONLY: put the more proximal compound
  ## (smaller/more negative z_score) on the x-axis, so the panel is
  ## upper-triangular. compound_a/compound_b in the underlying data (and in
  ## `dat` itself) are untouched -- only the *_plot columns feed the
  ## aes(x=, y=).
  dat$za_plot <- dat$z_score_a
  dat$zb_plot <- dat$z_score_b
  swap <- !is.na(dat$z_score_a) & !is.na(dat$z_score_b) & dat$z_score_a > dat$z_score_b
  if (any(swap)) {
    dat$za_plot[swap] <- dat$z_score_b[swap]
    dat$zb_plot[swap] <- dat$z_score_a[swap]
  }

  label_a <- .plot_label_nodes(proj, conditions, dat$compound_a, "compound")
  label_b <- .plot_label_nodes(proj, conditions, dat$compound_b, "compound")
  dat$pair_label <- paste(label_a, "+", label_b)
  dat$tooltip <- sprintf(
    "%s\nz_a=%.2f  z_b=%.2f\ns_ab=%s (%s)\nclass=%s",
    dat$pair_label, dat$za_plot, dat$zb_plot,
    ifelse(is.na(dat$s_ab), "NA", sprintf("%.2f", dat$s_ab)),
    ifelse(is.na(dat$separated), "NA", ifelse(dat$separated, "separated", "overlapping")),
    ifelse(is.na(dat$cheng_class), "NA", dat$cheng_class)
  )

  ## colour: separated, a two-level MANUAL scale (not a gradient) -- see
  ## roxygen. NA (jaccard mode / unmapped / disconnected) is a distinct grey.
  dat$separated_f <- factor(dat$separated, levels = c(TRUE, FALSE))
  ## shape: cheng_class, P1..P6 named + NA.
  dat$cheng_f <- factor(dat$cheng_class, levels = paste0("P", 1:6))
  ## size: |s_ab|, NA -> the smallest class (never dropped).
  dat$s_ab_mag <- abs(dat$s_ab)
  dat$size_val <- ifelse(is.na(dat$s_ab_mag), 0, dat$s_ab_mag)

  ## P2 = Complementary Exposure: separated & proximal_a & proximal_b,
  ## already computed by network_synergy() as cheng_class == "P2" /
  ## complementary_exposure. Per-panel count for the in-plot annotation.
  dat$is_p2 <- dat$cheng_class == "P2" & !is.na(dat$cheng_class)
  is_p2 <- dat$is_p2
  ann_p2 <- stats::aggregate(is_p2 ~ panel, data = dat, FUN = sum)
  names(ann_p2)[2] <- "n_p2"
  ann_p2$label <- sprintf("Complementary Exposure (P2)\nn = %d", ann_p2$n_p2)

  ## top_n P2 pairs by synergy_score (descending, NA last), per panel.
  p2_rows <- dat[is_p2, , drop = FALSE]
  top <- if (nrow(p2_rows) > 0 && top_n > 0) {
    do.call(rbind, lapply(split(p2_rows, p2_rows$panel), function(d) {
      d <- d[order(d$synergy_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
      utils::head(d, top_n)
    }))
  } else {
    p2_rows[FALSE, , drop = FALSE]
  }

  shape_values <- c(P1 = 16, P2 = 17, P3 = 15, P4 = 18, P5 = 3, P6 = 4)
  sep_colours <- c(`TRUE` = "#2980b9", `FALSE` = "#c0392b")

  p <- ggplot2::ggplot(dat, ggplot2::aes(
    x = .data$za_plot, y = .data$zb_plot,
    colour = .data$separated_f, shape = .data$cheng_f, size = .data$size_val
  ))
  p <- p +
    ggplot2::geom_rect(
      data = data.frame(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0),
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
      inherit.aes = FALSE, fill = "grey40", alpha = 0.08
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = "22") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "22")
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$pair_label), alpha = 0.75)
  } else {
    ggplot2::geom_point(alpha = 0.75)
  }
  if (nrow(top) > 0) {
    p <- p + ggplot2::geom_text(
      data = top, ggplot2::aes(x = .data$za_plot, y = .data$zb_plot, label = .data$pair_label),
      size = 2.6, colour = "grey15", vjust = -1, inherit.aes = FALSE, show.legend = FALSE
    )
  }
  p <- p + ggplot2::geom_text(
    data = ann_p2, ggplot2::aes(x = -Inf, y = -Inf, label = .data$label),
    hjust = -0.05, vjust = -0.6, size = 3, colour = "grey25", inherit.aes = FALSE
  )
  p <- p +
    ggplot2::scale_colour_manual(
      values = sep_colours, na.value = "grey70", drop = FALSE,
      breaks = c("TRUE", "FALSE"),
      labels = c("separated (s_AB >= 0)", "overlapping (s_AB < 0)"),
      name = "Topological separation"
    ) +
    ggplot2::scale_shape_manual(values = shape_values, na.value = 8, drop = FALSE, name = "Cheng class") +
    ggplot2::scale_size(range = c(1.5, 6), name = "|s_AB|\n(NA -> smallest)") +
    ggplot2::labs(
      title = paste0("Compound pair Cheng classification -- ", scope_label),
      subtitle = paste0(
        "x = z_score (more proximal compound), y = z_score (less proximal); shaded region = both individually proximal (z < 0);\n",
        "colour = separated (s_AB >= 0, Menche et al. 2015); shape = Cheng class P1-P6 (Cheng, Kovacs & Barabasi 2019); size = |s_AB|"
      ),
      x = "z_score (more proximal compound of the pair)",
      y = "z_score (less proximal compound of the pair)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40"))
  if (length(unique(dat$panel)) > 1) p <- p + ggplot2::facet_wrap(~panel)

  .synergy_finish_both(proj, p, p_sab, scope_label, engine, save, out_dir, width, height, dpi)
}

#' A self-explanatory empty panel for `plot_synergy()`'s two figures
#' @keywords internal
.synergy_empty_panel <- function(msg, title) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 3.4, colour = "grey25", lineheight = 1.1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
}

#' Save/log both `plot_synergy()` figures and chain the returned `proj`
#'
#' @description
#' `.plot_finish()` is written for one figure; `plot_synergy()` produces
#' two (main scatter + `s_AB` companion). Runs it twice, feeding the `proj`
#' that the first call attached (when `save = TRUE`) into the second, so
#' both log rows land in the same project -- then re-attaches the final
#' `proj` to the main result, which is the only one returned to the caller.
#' @keywords internal
.synergy_finish_both <- function(proj, p_main, p_sab, scope_label, engine, save, out_dir, width, height, dpi) {
  result_main <- .plot_finish(
    proj, p_main,
    name = "synergy_plot_log",
    filename = paste0("synergy_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  proj_next <- if (save) attr(result_main, "proj") else proj
  result_sab <- .plot_finish(
    proj_next, p_sab,
    name = "synergy_sab_plot_log",
    filename = paste0("synergy_sab_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = "static", save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  if (save) attr(result_main, "proj") <- attr(result_sab, "proj")
  result_main
}
