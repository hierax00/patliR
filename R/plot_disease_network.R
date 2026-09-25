#' @include AllGenerics.R internal.R network_build.R targets_disease.R plot-helpers.R plot_network_layers.R
NULL

## Compound-target-disease network, disease nodes drawn as large "blocks"
## with a convex-hull halo around the targets that belong to them --
## answers "which disease neighbourhood(s) does this compound's target set
## actually sit in", which plot_network_layers()'s generic disease layer
## (same point size as everything else, raw unresolved EFO ID as the label)
## does not make legible. Reuses the exact convex-hull pattern
## plot_chemical_space() already established (grDevices::chull() + a
## ggplot2::geom_polygon(alpha=..., colour=NA) fill, skipped below 3 points)
## rather than inventing a new one, and .network_layers_layer_colors() for
## palette consistency with the rest of the network plot family.

#' Compound-target-disease network, disease nodes as large blocks with a
#' convex-hull halo over their targets
#'
#' @description
#' Reads [targets_disease_profile()]'s output (not [targets_disease_filter()]
#' -- that table has no resolved disease name to label with) and draws a
#' three-layer network: compounds and targets at the usual small point size
#' (same palette as [plot_network_layers()]), disease nodes drawn
#' substantially larger, each with a translucent convex-hull polygon behind
#' it enclosing the positions of every target associated with that disease --
#' so a compound's target set visibly clusters toward the disease
#' neighbourhood(s) it is actually implicated in, rather than every disease
#' looking like just another same-size node in the layout.
#'
#' Only targets that have at least one disease association in
#' [targets_disease_profile()]'s output are drawn -- this is deliberately a
#' focused view of the compound-target-disease structure, not a full replay
#' of every predicted target ([plot_network_layers()] already does that).
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease `NULL` (default, every disease in scope) or a character
#'   vector of `disease_id`s to restrict to.
#' @param max_rank `NULL` (default, every row) or an integer -- keep only
#'   `targets_disease_profile()` rows with `rank <= max_rank` (its explore
#'   mode can carry several diseases per target; this trims the plot to each
#'   target's top few without re-fetching).
#' @param compound_ids Character vector of `compounds(proj)$id` to restrict
#'   to, or `NULL` (default) for every compound present in scope.
#' @param top_n_labels Integer, default `25`: only this many disease nodes
#'   -- those with the most associated targets -- get a text label (full
#'   names stay in the `ggiraph` tooltips). `0` draws no label.
#' @param top_n_targets Integer, default `15`: the targets hit by the most
#'   compounds in the figure (ties broken by their highest disease
#'   association score) get a text label (gene symbol when `clusterProfiler`
#'   + `org.Hs.eg.db` are installed, UniProt ID otherwise). `0` labels no
#'   target.
#' @param label_compounds Logical, default `TRUE`: label every compound node
#'   with its (truncated) name. Set `FALSE` for a very large compound set.
#' @param engine `"static"` (default) or `"ggiraph"`.
#'
#' @return A `ggplot`/`girafe` object, with `attr(., "proj")` set when
#'   `save = TRUE` (also writes a PNG, logged to
#'   `patliRResults(proj, "disease_network_plot_log")`).
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... build the project, run network_build() and targets_disease_profile() ...
#' plot_disease_network(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @seealso [targets_disease_profile()], [plot_network_layers()]
#' @export
plot_disease_network <- function(proj, condition = NULL, disease = NULL, max_rank = NULL,
                                  compound_ids = NULL, top_n_labels = 25,
                                  top_n_targets = 15, label_compounds = TRUE,
                                  engine = c("static", "ggiraph"),
                                  save = TRUE, out_dir = NULL, width = 9, height = 7, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n_labels), length(top_n_labels) == 1, !is.na(top_n_labels), top_n_labels >= 0)
  if (!(is.numeric(top_n_targets) && length(top_n_targets) == 1 && !is.na(top_n_targets) && top_n_targets >= 0)) {
    cli::cli_abort("{.arg top_n_targets} must be a single number >= 0.")
  }
  if (!(is.logical(label_compounds) && length(label_compounds) == 1 && !is.na(label_compounds))) {
    cli::cli_abort("{.arg label_compounds} must be {.code TRUE} or {.code FALSE}.")
  }
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  profile <- patliRResults(proj, "targets_disease_profile")
  if (is.null(profile) || nrow(profile) == 0) {
    cli::cli_abort(c(
      "No {.val targets_disease_profile} entry in {.arg proj}.",
      "i" = "Run {.fn targets_disease_profile} first."
    ))
  }
  if (!is.null(disease)) profile <- profile[profile$disease_id %in% disease, , drop = FALSE]
  if (!is.null(max_rank)) profile <- profile[profile$rank <= max_rank, , drop = FALSE]

  edges_all <- patliRResults(proj, "network_edges")
  edges <- edges_all[edges_all$condition %in% conditions, c("compound_id", "uniprot_id", "weight"), drop = FALSE]
  if (!is.null(compound_ids)) edges <- edges[edges$compound_id %in% compound_ids, , drop = FALSE]

  profile <- merge(profile, unique(edges[c("compound_id", "uniprot_id")]),
                   by.x = c("compound_id", "target_id"),
                   by.y = c("compound_id", "uniprot_id"))
  if (nrow(profile) == 0) {
    cli::cli_abort("No {.val targets_disease_profile} rows for the requested compounds/condition(s)/disease/max_rank.")
  }

  ## restrict compound-target edges to targets that actually carry disease
  ## data -- this plot is deliberately the disease-relevant subgraph, not
  ## every predicted target (plot_network_layers() already covers that).
  ct_edges <- unique(merge(edges, unique(profile[c("compound_id", "target_id")]),
                           by.x = c("compound_id", "uniprot_id"),
                           by.y = c("compound_id", "target_id")))
  td_edges <- unique(profile[, c("target_id", "disease_id", "disease_name", "association_score")])

  compound_nodes <- sort(unique(ct_edges$compound_id))
  target_nodes <- sort(unique(c(ct_edges$uniprot_id, td_edges$target_id)))
  disease_nodes <- sort(unique(td_edges$disease_id))
  disease_name_of <- stats::setNames(td_edges$disease_name, td_edges$disease_id)[disease_nodes]

  node_ids <- c(compound_nodes, target_nodes, disease_nodes)
  node_layer <- c(
    rep("compound", length(compound_nodes)), rep("target", length(target_nodes)),
    rep("disease", length(disease_nodes))
  )
  node_label <- c(
    .plot_label_nodes(proj, conditions, compound_nodes, "compound"),
    .plot_label_nodes(proj, conditions, target_nodes, "target"),
    unname(ifelse(is.na(disease_name_of) | !nzchar(disease_name_of), disease_nodes, disease_name_of))
  )

  g <- igraph::graph_from_data_frame(
    d = rbind(
      data.frame(from = ct_edges$compound_id, to = ct_edges$uniprot_id, stringsAsFactors = FALSE),
      data.frame(from = td_edges$target_id, to = td_edges$disease_id, stringsAsFactors = FALSE)
    ),
    vertices = data.frame(name = node_ids, layer = node_layer, label = node_label, stringsAsFactors = FALSE),
    directed = FALSE
  )

  coords <- igraph::layout_with_fr(g, niter = 2000)
  nodes <- data.frame(
    id = igraph::V(g)$name, layer = igraph::V(g)$layer, label = igraph::V(g)$label,
    x = coords[, 1], y = coords[, 2], stringsAsFactors = FALSE
  )
  pos_of <- stats::setNames(seq_len(nrow(nodes)), nodes$id)

  palette <- .network_layers_layer_colors()
  nodes$point_size <- ifelse(nodes$layer == "disease", 9, ifelse(nodes$layer == "target", 2.6, 3.2))

  ## disease hull: convex hull over the positions of every target connected
  ## to that disease -- same grDevices::chull() pattern as
  ## .chemical_space_ggplot()'s per-family halos, skipped below 3 points.
  hulls <- do.call(rbind, lapply(disease_nodes, function(d) {
    tids <- unique(td_edges$target_id[td_edges$disease_id == d])
    idx <- pos_of[tids]
    idx <- idx[!is.na(idx)]
    if (length(idx) < 3) return(NULL)
    pts <- nodes[idx, c("x", "y"), drop = FALSE]
    hull_idx <- grDevices::chull(pts$x, pts$y)
    data.frame(disease_id = d, x = pts$x[hull_idx], y = pts$y[hull_idx], stringsAsFactors = FALSE)
  }))

  edge_df <- as.data.frame(igraph::as_edgelist(g), stringsAsFactors = FALSE)
  names(edge_df) <- c("from", "to")
  edge_df$x <- nodes$x[pos_of[edge_df$from]]
  edge_df$y <- nodes$y[pos_of[edge_df$from]]
  edge_df$xend <- nodes$x[pos_of[edge_df$to]]
  edge_df$yend <- nodes$y[pos_of[edge_df$to]]

  nodes$tooltip <- sprintf("%s\n%s", nodes$label, nodes$layer)

  p <- ggplot2::ggplot()
  if (!is.null(hulls) && nrow(hulls) > 0) {
    p <- p + ggplot2::geom_polygon(
      data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$disease_id),
      fill = palette[["disease"]], alpha = 0.14, colour = NA
    )
  }
  p <- p + ggplot2::geom_segment(
    data = edge_df, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
    colour = "grey70", linewidth = 0.3, alpha = 0.6
  )
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$layer, size = .data$point_size,
                   tooltip = .data$tooltip, data_id = .data$id)
    )
  } else {
    ggplot2::geom_point(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$layer, size = .data$point_size))
  }
  ## drawn label shortened, full disease name kept in the tooltip
  ## Only the `top_n_labels` diseases with the most associated targets are
  ## named -- with hundreds of disease nodes, labelling all of them prints an
  ## unreadable block of text over the whole network.
  disease_labels <- nodes[nodes$layer == "disease", , drop = FALSE]
  n_targets_of <- table(td_edges$disease_id)
  disease_labels$n_targets <- as.integer(n_targets_of[disease_labels$id])
  disease_labels <- utils::head(disease_labels[order(-disease_labels$n_targets, disease_labels$label), , drop = FALSE], top_n_labels)
  disease_labels$short_label <- .plot_truncate(disease_labels$label)
  if (length(disease_nodes) > 60) {
    cli::cli_inform(c(
      "i" = "{.fn plot_disease_network}: {length(disease_nodes)} disease nodes -- only the {min(top_n_labels, length(disease_nodes))} with the most targets are labeled; narrow the figure with {.arg disease}, {.arg max_rank} (e.g. {.code max_rank = 1}) or {.arg compound_ids}."
    ))
  }
  ## targets: the top_n_targets hit by the most compounds (ties: highest
  ## disease association score); compounds: all of them, names shortened
  top_targets <- .disease_network_top_targets(ct_edges, td_edges, top_n_targets)
  target_labels <- nodes[nodes$layer == "target" & nodes$id %in% top_targets, , drop = FALSE]
  target_labels$short_label <- .plot_truncate(target_labels$label, 14)
  compound_labels <- if (label_compounds) nodes[nodes$layer == "compound", , drop = FALSE] else nodes[0, , drop = FALSE]
  compound_labels$short_label <- .plot_truncate(compound_labels$label, 20)
  label_df <- rbind(
    disease_labels[, c("id", "layer", "x", "y", "short_label")],
    target_labels[, c("id", "layer", "x", "y", "short_label")],
    compound_labels[, c("id", "layer", "x", "y", "short_label")]
  )
  if (nrow(label_df) > 0) {
    ## one repelled layer for every label, so labels of different node types
    ## avoid each other too; text colour = a darker shade of the node colour
    text_col <- c(disease = "#7b241c", target = "#1b4f72", compound = "#8a4b0f")
    p <- p + .plot_text_layer(
      data = label_df, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$short_label),
      size = unname(c(disease = 3.1, target = 2.5, compound = 2.4)[label_df$layer]),
      fontface = unname(c(disease = "bold", target = "bold.italic", compound = "plain")[label_df$layer]),
      colour = unname(text_col[label_df$layer]), inherit.aes = FALSE, show.legend = FALSE,
      ## point.size = the drawn node size, so labels clear the large disease blocks
      repel_args = list(point.size = unname(c(disease = 9, target = 2.6, compound = 3.2)[label_df$layer]),
                        box.padding = 0.25, point.padding = 0.15, min.segment.length = 0.2,
                        segment.colour = "grey55", segment.size = 0.25,
                        bg.colour = "white", bg.r = 0.12, max.overlaps = Inf, seed = 1),
      text_args = list(vjust = -1.2)
    )
  }
  n_dis <- length(unique(td_edges$disease_id))
  subtitle <- paste0(
    "disease nodes enlarged, with a convex-hull halo over their associated targets",
    if (n_dis == 1) paste0(" -- ", .plot_disease_label(proj, disease_nodes)) else "",
    "\nlabels: ",
    paste(c(
      if (nrow(disease_labels) > 0) if (nrow(disease_labels) < n_dis) paste0(nrow(disease_labels), " of ", n_dis, " diseases (most targets)") else "diseases",
      if (nrow(target_labels) > 0) paste0("top ", nrow(target_labels), " targets by number of compounds (italic)"),
      if (nrow(compound_labels) > 0) "compounds"
    ), collapse = ", ")
  )
  p <- p +
    ggplot2::scale_colour_manual(values = palette, name = "Layer") +
    ggplot2::scale_size_identity() +
    ggplot2::labs(
      title = paste0("Compound-target-disease network -- ", scope_label),
      subtitle = .plot_wrap(subtitle, .plot_wrap_width(width, 7.5))
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40"),
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right"
    )

  .plot_finish(
    proj, p,
    name = "disease_network_plot_log",
    filename = paste0("disease_network_", scope_label, ".png"),
    log_row = data.frame(condition = scope_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = "condition",
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Targets to label in `plot_disease_network()`
#' @description Ranked by the number of distinct compounds hitting the
#'   target in the drawn subgraph, then by its highest disease association
#'   score, then by ID.
#' @param ct_edges Compound-target edges (`compound_id`, `uniprot_id`).
#' @param td_edges Target-disease edges (`target_id`, `association_score`).
#' @param n Number of targets to return.
#' @return Character vector of at most `n` target IDs, best first.
#' @keywords internal
.disease_network_top_targets <- function(ct_edges, td_edges, n) {
  if (n <= 0 || nrow(ct_edges) == 0) return(character(0))
  e <- unique(ct_edges[, c("compound_id", "uniprot_id"), drop = FALSE])
  n_cmp <- table(e$uniprot_id)
  ids <- names(n_cmp)
  best <- vapply(ids, function(t) {
    s <- td_edges$association_score[td_edges$target_id == t]
    if (length(s) == 0 || all(is.na(s))) -Inf else max(s, na.rm = TRUE)
  }, numeric(1))
  ord <- order(-as.integer(n_cmp), -best, ids)
  utils::head(ids[ord], n)
}
