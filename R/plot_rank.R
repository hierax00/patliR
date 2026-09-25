#' @include AllGenerics.R internal.R rank_candidates.R plot-helpers.R plot_heatmap.R
NULL

## Two views over rank_candidates()'s output, same view= convention as
## plot_proximity(view = "z"/"null"). "pareto" is new geometry (plain
## ggplot2, no dependency): a from-scratch non-dominated sort
## (.rank_pareto_front(), computed once inside rank_candidates() itself,
## over every crit_* column that call actually used -- not just the two
## plotted axes). "heatmap" reuses .plot_pheatmap_render() (extracted from
## plot_heatmap() in this same commit) after subsetting to the top-ranked
## compounds x most-relevant targets -- no new geometry, no new dependency.

#' Visualise `rank_candidates()`'s output: a Pareto scatter, or a top
#' compounds x top targets heatmap
#'
#' @description
#' `view = "pareto"` (default): a 2D scatter over `rank_candidates()`'s
#' output, coloured by `pareto_front` (the non-dominated tiering
#' [rank_candidates()] already computed over every criterion that call
#' actually used -- not recomputed here, and not limited to the two plotted
#' axes). Front-1 (non-dominated) compounds are labelled by name. `x`
#' defaults to `"rra_rank"` (the single combined answer); `y` defaults to
#' whichever optional criterion is present first in this priority order:
#' `crit_synergy_best`, `crit_proximity_z`, `crit_module_r_index`,
#' `crit_hub_penalty` (falling back to `crit_centrality` if none of those
#' ran) -- the criteria most likely to actually differentiate compounds,
#' since `x`/`y` are only ever
#' an illustrative 2D placement of a decision that used every criterion, not
#' the decision itself. Both are freely overridable to any numeric column
#' of `patliRResults(proj, "rank_candidates")` (any `crit_*`, `rank_*`,
#' `rra_score`, `rra_rank`, or `n_criteria_used`).
#'
#' `view = "heatmap"`: compound x target matrix restricted to the top
#' `top_n_compounds` compounds (by `rra_rank`) and the `top_n_targets` most
#' relevant targets among that restricted set (`target_relevance`), filled
#' by [network_build()]'s edge weight -- delegates to the same
#' [pheatmap::pheatmap()] rendering [plot_heatmap()] uses (via the shared
#' `.plot_pheatmap_render()` helper), just on a subset. Answers "which
#' targets do the best compounds actually converge on", which the full,
#' every-compound x every-target matrix in `plot_heatmap(what =
#' "compound_target")` is too dense to show.
#'
#' @section `target_relevance` (heatmap view only):
#' \describe{
#'   \item{`"breadth"` (default)}{Number of distinct top-`top_n_compounds`
#'     compounds hitting each target -- the targets several best compounds
#'     converge on, the mechanistically interesting overlap.}
#'   \item{`"weight"`}{Sum of [network_build()] edge weight (import
#'     probability) across the top compounds, per target.}
#'   \item{`"centrality"`}{The target's own composite centrality (same
#'     `rowMeans()` of available `*_norm` columns [rank_candidates()]'s
#'     `"centrality"` criterion uses), independent of which top compounds
#'     hit it. Requires `network_centrality`.}
#' }
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param view `"pareto"` (default) or `"heatmap"`.
#' @param x,y Column name(s) from `rank_candidates`'s output for
#'   `view = "pareto"`'s axes. `NULL` (default) uses `x = "rra_rank"` and
#'   the priority default for `y` described above.
#' @param top_n_compounds Integer, default `15` -- `view = "heatmap"` only.
#' @param top_n_targets Integer, default `20` -- `view = "heatmap"` only.
#' @param target_relevance `"breadth"` (default), `"weight"`, or
#'   `"centrality"` -- `view = "heatmap"` only, see the section above.
#' @param label `view = "pareto"` only: `"front"` (default) names every
#'   Pareto-front-1 (non-dominated) compound; `"all"` names every compound
#'   (front 1 in dark bold, the rest in grey); `"none"` draws no names.
#'   Names are truncated and repelled (`ggrepel`), never dropped for
#'   overlap.
#' @param engine `"static"` (default) or `"ggiraph"` -- `view = "pareto"`
#'   only (`pheatmap` has no interactive engine, same contract difference
#'   [plot_heatmap()] documents).
#'
#' @return `view = "pareto"`: a `ggplot`/`girafe` object. `view = "heatmap"`:
#'   the `pheatmap` object (`list(tree_row, tree_col, kmeans, gtable)`). Both
#'   carry the updated `proj` as `attr(., "proj")` when `save = TRUE`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... build the project, run network_build/network_centrality/adme_local/
#' # adme_filter/rank_candidates() ...
#' plot_rank(proj, condition = "FLO-ET", view = "pareto", save = FALSE)
#' plot_rank(proj, condition = "FLO-ET", view = "heatmap", save = FALSE)
#' }
#'
#' @seealso [rank_candidates()], [plot_heatmap()]
#' @export
plot_rank <- function(proj, condition = NULL, view = c("pareto", "heatmap"),
                       x = NULL, y = NULL,
                       top_n_compounds = 15, top_n_targets = 20,
                       target_relevance = c("breadth", "weight", "centrality"),
                       label = c("front", "all", "none"),
                       engine = c("static", "ggiraph"),
                       save = TRUE, out_dir = NULL, width = NULL, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  view <- match.arg(view)
  target_relevance <- match.arg(target_relevance)
  label <- match.arg(label)
  engine <- match.arg(engine)
  stopifnot(is.numeric(top_n_compounds), length(top_n_compounds) == 1, top_n_compounds >= 1)
  stopifnot(is.numeric(top_n_targets), length(top_n_targets) == 1, top_n_targets >= 1)

  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  rc_all <- patliRResults(proj, "rank_candidates")
  if (is.null(rc_all) || nrow(rc_all) == 0) {
    cli::cli_abort(c("No {.val rank_candidates} entry in {.arg proj}.", "i" = "Run {.fn rank_candidates} first."))
  }
  dat <- rc_all[rc_all$condition %in% conditions, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val rank_candidates} rows for the requested condition(s).")
  }

  if (view == "pareto") {
    engine <- .plot_require(engine)
    .plot_rank_pareto(proj, dat, conditions, scope_label, x, y, engine, save, out_dir, width, height, dpi, label)
  } else {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      cli::cli_abort("The {.pkg pheatmap} package is required for {.fn plot_rank}({.code view = \"heatmap\"}).")
    }
    .plot_rank_heatmap(proj, dat, conditions, scope_label, top_n_compounds, top_n_targets, target_relevance, save, out_dir, width, height)
  }
}

#' @keywords internal
.plot_rank_y_priority <- c("crit_synergy_best", "crit_proximity_z", "crit_module_r_index", "crit_hub_penalty", "crit_centrality", "crit_adme_pass_frac")

#' @keywords internal
.plot_rank_pareto <- function(proj, dat, conditions, scope_label, x, y, engine, save, out_dir, width, height, dpi,
                              label = "front") {
  if (is.null(x)) x <- "rra_rank"
  if (is.null(y)) {
    present <- .plot_rank_y_priority[.plot_rank_y_priority %in% names(dat)]
    present <- present[vapply(dat[present], function(x) is.numeric(x) && any(is.finite(x)), logical(1))]
    if (length(present) == 0) {
      cli::cli_abort("No usable numeric column found to default {.arg y} to; pass {.arg y} explicitly.")
    }
    y <- present[1]
  }
  if (!x %in% names(dat)) cli::cli_abort("{.arg x} = {.val {x}} is not a column of {.val rank_candidates}.")
  if (!y %in% names(dat)) cli::cli_abort("{.arg y} = {.val {y}} is not a column of {.val rank_candidates}.")

  dat$panel <- dat$condition
  dat$front_f <- factor(ifelse(dat$pareto_front == 1L, "1 (non-dominated)", "2+"), levels = c("1 (non-dominated)", "2+"))
  labels <- .plot_label_nodes(proj, conditions, dat$compound_id, "compound")
  dat$compound_label <- labels
  dat$tooltip <- sprintf("%s\n%s = %s\n%s = %s\nrra_rank = %d, front = %d", dat$compound_label, x, signif(dat[[x]], 3), y, signif(dat[[y]], 3), dat$rra_rank, dat$pareto_front)

  front1 <- .plot_rank_label_rows(dat, label)
  if (is.null(width)) width <- 8

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data$front_f))
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$compound_label), size = 2.5, alpha = 0.85)
  } else {
    ggplot2::geom_point(size = 2.5, alpha = 0.85)
  }
  ## Front-1 labels: shortened for drawing (tooltip keeps the full name) and
  ## repelled when ggrepel is installed -- a label on a point at the right
  ## edge (the worst rra_rank) otherwise runs off the panel.
  if (nrow(front1) > 0) {
    front1$short_label <- .plot_truncate(front1$compound_label)
    is_f1 <- front1$pareto_front == 1L & !is.na(front1$pareto_front)
    p <- p + .plot_text_layer(
      data = front1,
      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], label = .data$short_label),
      size = 2.6, inherit.aes = FALSE,
      colour = ifelse(is_f1, "grey10", "grey45"),
      fontface = ifelse(is_f1 & label == "all", "bold", "plain"),
      repel_args = list(box.padding = 0.3, min.segment.length = 0.2, segment.colour = "grey60",
                        bg.colour = "white", bg.r = 0.1, max.overlaps = Inf, seed = 1),
      text_args = list(vjust = -1)
    )
  }
  ## headroom on both sides of a numeric x so edge labels have room either way
  if (is.numeric(dat[[x]])) p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.12))
  p <- p +
    ggplot2::scale_colour_manual(values = c(`1 (non-dominated)` = "#c0392b", `2+` = "#7f8c8d"), name = "Pareto front") +
    ggplot2::labs(
      title = paste0("Candidate ranking -- ", scope_label),
      subtitle = .plot_wrap(
        paste0("x = ", x, ", y = ", y, "; colour = Pareto front (non-dominated over every criterion rank_candidates() used, not just x/y)"),
        .plot_wrap_width(width, 7.5)
      ),
      x = x, y = y
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40"),
      plot.title.position = "plot"
    )
  if (length(unique(dat$panel)) > 1) p <- p + ggplot2::facet_wrap(~panel)
  ## a non-default label set gets its own file and log row; logs written
  ## before `label` existed hold the default ("front") figure
  if (save) proj <- .plot_log_backfill(proj, "rank_plot_log", "label", "front")
  .plot_finish(
    proj, p,
    name = "rank_plot_log",
    filename = paste0("rank_pareto_", scope_label, if (label != "front") paste0("_labels_", label) else "", ".png"),
    log_row = data.frame(condition = scope_label, view = "pareto", label = label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "view", "label"),
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Compounds to name on the Pareto view
#' @description `"front"`: every front-1 (non-dominated) compound; `"all"`:
#'   every compound; `"none"`: nothing.
#' @return Subset of `dat`.
#' @keywords internal
.plot_rank_label_rows <- function(dat, label = c("front", "all", "none")) {
  label <- match.arg(label)
  switch(label,
    front = dat[!is.na(dat$pareto_front) & dat$pareto_front == 1L, , drop = FALSE],
    all = dat,
    none = dat[FALSE, , drop = FALSE]
  )
}

#' @keywords internal
.plot_rank_heatmap <- function(proj, dat, conditions, scope_label, top_n_compounds, top_n_targets, target_relevance, save, out_dir, width, height) {
  top_compounds <- unique(unlist(lapply(split(dat, dat$condition), function(d) {
    d <- d[order(d$rra_rank), , drop = FALSE]
    utils::head(d$compound_id, top_n_compounds)
  })))

  edges_all <- patliRResults(proj, "network_edges")
  edges_top <- edges_all[edges_all$condition %in% conditions & edges_all$compound_id %in% top_compounds, c("compound_id", "uniprot_id", "weight"), drop = FALSE]
  if (nrow(edges_top) == 0) {
    cli::cli_abort("No {.val network_edges} rows for the top-ranked compounds in this scope.")
  }

  relevance <- if (target_relevance == "breadth") {
    stats::aggregate(compound_id ~ uniprot_id, edges_top, function(x) length(unique(x)))
  } else if (target_relevance == "weight") {
    stats::aggregate(edges_top["weight"], edges_top["uniprot_id"],
                     function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE))
  } else {
    cent <- patliRResults(proj, "network_centrality")
    if (is.null(cent) || nrow(cent) == 0) {
      cli::cli_abort(c(
        "{.arg target_relevance} = \"centrality\" needs {.val network_centrality}.",
        "i" = "Run {.fn network_centrality} first, or use {.code target_relevance = \"breadth\"}/{.code \"weight\"}."
      ))
    }
    csub <- cent[cent$condition %in% conditions & cent$node_type == "target" &
                   cent$node_id %in% edges_top$uniprot_id, , drop = FALSE]
    if (nrow(csub) == 0) cli::cli_abort("No centrality rows for targets hit by the top compounds in this scope.")
    norm_cols <- intersect(c("degree_norm", "betweenness_norm", "hub_score_component"), names(csub))
    csub$composite <- if (length(norm_cols) > 0) rowMeans(as.matrix(csub[, norm_cols, drop = FALSE]), na.rm = TRUE) else NA_real_
    stats::aggregate(csub["composite"], csub["node_id"],
                     function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  }
  names(relevance) <- c("uniprot_id", "score")
  relevance <- relevance[order(relevance$score, decreasing = TRUE), , drop = FALSE]
  top_targets <- utils::head(relevance$uniprot_id, top_n_targets)

  long <- edges_top[edges_top$uniprot_id %in% top_targets, , drop = FALSE]
  if (nrow(long) == 0) {
    cli::cli_abort("Nothing to plot: no edges between the top compounds and the top targets in this scope.")
  }
  mat <- .plot_edge_matrix(proj, conditions, long)

  title <- paste0(
    "Top ", length(top_compounds), " Compounds x Top ", ncol(mat), " Targets (", target_relevance, ") -- ", scope_label
  )
  .plot_pheatmap_render(
    proj, mat, title, save, out_dir, width, height,
    log_name = "rank_heatmap_plot_log",
    log_row = data.frame(condition = scope_label, target_relevance = target_relevance, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "target_relevance"),
    filename = paste0("rank_heatmap_", scope_label, "_", target_relevance, ".pdf")
  )
}
