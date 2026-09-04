#' @include AllGenerics.R internal.R network_centrality.R plot-helpers.R
NULL

#' Bar plot of the top-`top_n` nodes by a `network_centrality()` measure
#'
#' @description
#' Horizontal bar chart of the `top_n` highest-scoring nodes (labeled with
#' their real name/symbol, same lookup as [plot_network_layers()]) for one
#' `measure` (`"degree"`, `"betweenness"`, or `"hub_score"`), colored by
#' `node_type` (`"compound"`/`"target"`), faceted by `condition` when more
#' than one is in scope.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param measure Single measure to plot: `"degree"`, `"betweenness"`,
#'   `"hub_score"`, or the bipartite-normalised `"degree_norm"` /
#'   `"betweenness_norm"` -- must already be a column of
#'   `patliRResults(proj, "network_centrality")`, i.e. requested (and, for
#'   the `_norm` variants, `normalize = TRUE`) in the [network_centrality()]
#'   call that produced it. When `measure` is not supplied **and**
#'   `node_type = "both"`, it defaults to `"degree_norm"` (raw `degree` is
#'   not comparable between compounds and targets; Borgatti & Everett
#'   1997), falling back to `"degree"` with a warning if no normalised
#'   `degree_norm` value is available for the conditions in scope (a result
#'   from before patliR 0.2.0, or one produced with `normalize = FALSE`).
#'   Otherwise the default is `"degree"`.
#' @param node_type `"both"` (default), `"compound"`, or `"target"` --
#'   restrict to one side of the bipartite graph.
#' @param top_n Integer, default `20`. Only the `top_n` highest-`measure`
#'   nodes (within the `condition`/`node_type` scope) are drawn -- with
#'   hundreds of targets, plotting all of them produces an unreadable wall
#'   of bars.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "centrality_plot_log")`.
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
#' proj <- network_centrality(proj)
#' plot_centrality(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_centrality <- function(proj, condition = NULL,
                             measure = c("degree", "betweenness", "hub_score",
                                         "degree_norm", "betweenness_norm"),
                             node_type = c("both", "compound", "target"), top_n = 20,
                             engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                             width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  measure_supplied <- !missing(measure)
  measure <- match.arg(measure)
  node_type <- match.arg(node_type)
  engine <- match.arg(engine)
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  cent <- patliRResults(proj, "network_centrality")
  if (is.null(cent) || nrow(cent) == 0) {
    cli::cli_abort(c("No {.val network_centrality} entry in {.arg proj}.", "i" = "Run {.fn network_centrality} first."))
  }
  if (!measure_supplied && node_type == "both") {
    in_scope <- cent[cent$condition %in% conditions, , drop = FALSE]
    if ("degree_norm" %in% names(cent) && any(!is.na(in_scope$degree_norm))) {
      measure <- "degree_norm"
      cli::cli_inform(c(
        "i" = "Defaulting {.arg measure} to {.val degree_norm}: with {.arg node_type} = {.val both}, raw {.val degree} is not comparable between compounds and targets (Borgatti & Everett 1997)."
      ))
    } else {
      measure <- "degree"
      cli::cli_warn(c(
        "!" = "{.arg node_type} is {.val both} but no bipartite-normalised {.val degree_norm} is available (this result predates patliR 0.2.0, or was produced with {.code normalize = FALSE}); using raw {.val degree}, which is not comparable across node types.",
        "i" = "Re-run {.fn network_centrality} with {.code normalize = TRUE}."
      ))
    }
  }
  if (!measure %in% names(cent)) {
    cli::cli_abort("Column {.val {measure}} not found in {.val network_centrality} -- was it requested in the {.arg measures} (and {.code normalize = TRUE}) of {.fn network_centrality}?")
  }
  dat <- cent[cent$condition %in% conditions, , drop = FALSE]
  if (node_type != "both") dat <- dat[dat$node_type == node_type, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_centrality} rows for condition(s) {.val {conditions}} / node_type {.val {node_type}}.")
  }
  dat$label <- .plot_label_nodes(proj, conditions, dat$node_id, dat$node_type)
  dat <- dat[order(dat$condition, -dat[[measure]]), , drop = FALSE]
  dat <- do.call(rbind, lapply(split(dat, dat$condition), utils::head, top_n))
  ## Two different nodes can resolve to the same gene symbol; without this
  ## they collapse into a single bar under factor(). Append the accession
  ## on collision (spec 2.3 / 3.3).
  ids_per_label <- tapply(dat$node_id, dat$label, function(x) length(unique(x)))
  collide <- names(ids_per_label)[ids_per_label > 1L]
  if (length(collide)) {
    hit <- dat$label %in% collide
    dat$label[hit] <- paste0(dat$label[hit], " (", dat$node_id[hit], ")")
  }
  dat$label <- factor(dat$label, levels = unique(dat$label[order(dat[[measure]])]))
  dat$tooltip <- sprintf("%s (%s)\n%s: %.3g", dat$label, dat$node_type, measure, dat[[measure]])

  y_lab <- .plot_centrality_ylab(measure, node_type, dat)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$label, y = .data[[measure]], fill = .data$node_type))
  p <- if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p + ggiraph::geom_col_interactive(ggplot2::aes(tooltip = .data$tooltip))
  } else {
    p + ggplot2::geom_col()
  }
  p <- p +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(compound = "#e67e22", target = "#2980b9"), name = "Node type") +
    ggplot2::labs(
      title = paste0("Top ", top_n, " nodes by ", measure, " -- ", scope_label),
      x = NULL, y = y_lab
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))
  if (length(conditions) > 1) p <- p + ggplot2::facet_wrap(~condition, scales = "free_y")

  .plot_finish(
    proj, p,
    name = "centrality_plot_log",
    filename = paste0("centrality_", measure, "_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, measure = measure, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "measure"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}

#' Axis label for `plot_centrality()`, stating the bipartite normaliser
#'
#' @description
#' For `"degree_norm"` / `"betweenness_norm"` the value axis carries a
#' subtitle naming what the raw score was divided by (Borgatti & Everett
#' 1997); raw measures keep the bare measure name.
#'
#' @param measure The resolved measure string.
#' @param node_type `"both"` / `"compound"` / `"target"`.
#' @param dat The plotting frame (needs `n_compounds` / `n_targets` for the
#'   `_norm` variants).
#' @return A single string for `ggplot2::labs(y = )`.
#' @keywords internal
.plot_centrality_ylab <- function(measure, node_type, dat) {
  collapse_num <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = "/")
  if (measure == "degree_norm") {
    nt <- if ("n_targets"   %in% names(dat)) collapse_num(dat$n_targets)   else "|T|"
    nc <- if ("n_compounds" %in% names(dat)) collapse_num(dat$n_compounds) else "|C|"
    sub <- switch(
      node_type,
      compound = paste0("/ |T| = ", nt),
      target   = paste0("/ |C| = ", nc),
      both     = paste0("compound / |T| = ", nt, ", target / |C| = ", nc)
    )
    return(paste0("degree_norm  (", sub, ")"))
  }
  if (measure == "betweenness_norm") {
    return("betweenness_norm  (/ B_max, Borgatti & Everett 1997)")
  }
  measure
}
