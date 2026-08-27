#' @include AllGenerics.R internal.R network_layers.R
NULL

## plot_network_layers() -- first of the pending "interactive/static plot_*
## for network_*" items (patliR_manual.md, section 6 wrap-up, 2026-07-23).
## Confirmed with Uriel (chat, 2026-07-23) after comparing ggiraph (already
## used by plot_boiled_egg()/plot_admet_radar(), zero new Suggests) against
## visNetwork/networkD3 (real in-browser force layout, but a new
## dependency and no direct ggsave()-equivalent for a publication PNG):
## start with the tripartite/layered network (.network_layered_graph(),
## network_layers.R) as a STATIC figure, ggplot2 only -- the reference
## figures Uriel shared (all from published network-pharmacology papers)
## are static, high-resolution, hairball-style plots with categorical node
## colors and only hub nodes labeled, which is exactly what ggsave() from
## plain ggplot2 already produces well. `visNetwork`-style live/interactive
## exploration is parked as a separate future item, not built here.
##
## No `ggraph` (not in Suggests, and would be the "normal" way to plot an
## igraph object in ggplot2) -- the layout is computed directly with
## `igraph::layout_with_*()` and turned into plain node/edge data.frames,
## then drawn with `geom_segment()`/`geom_point()`/`geom_text()`. Same
## "keep Suggests minimal" principle as plot_boiled_egg()'s no-ggrepel
## note.
##
## Two real corrections after Uriel's first live renders (2026-07-23),
## both confirmed in chat:
## 1. "Dandelion" starburst, not hairball -- network_enrich() easily
##    returns hundreds of significant GO/Reactome terms, each a degree-1
##    leaf hanging off its target, drowning the actual compound-target
##    structure. Fix: `layers` defaults to c("compound", "target") only;
##    `"pathway"`/`"disease"` are opt-in, and `max_pathways` caps the
##    former to the most significant terms when requested.
## 2. Uriel shared Yildirim et al. 2007 ("Drug-target network", Nat
##    Biotechnol 25, 1119-1126, doi:10.1038/nbt1338) as the real target
##    look -- confirmed by reading it: that figure is the FULL FDA-approved
##    drug-target network (890 drugs x 394 targets, every drug at once),
##    NOT filtered by disease at all; disease-gene proximity is a *separate*
##    analysis there (interactome shortest-path), matching what
##    network_proximity() already does in this package, not something drawn
##    as graph nodes. Confirms leaving "disease" out of the default layers
##    is correct, not just a pathway-specific fix. It also means the
##    clustering-by-category look only emerges at the *global* scale (many
##    compounds sharing targets) -- a single small condition's bipartite
##    graph structurally cannot look like that. Uriel confirmed (chat,
##    2026-07-23) both a per-condition view (existing) and a global,
##    all-conditions-pooled view should be supported; `condition = NULL`
##    is now the global/pooled mode, matching the rest of the `network_*`
##    family's own NULL-means-"every built condition" convention (see
##    `.network_resolve_conditions()`) -- except here it pools them into
##    one combined graph instead of looping per condition.

#' Static network plot of the compound-target(-pathway)(-disease) layered
#' graph, for one condition or pooled across every built condition
#'
#' @description
#' Force-directed layout of `.network_layered_graph()`'s underlying data
#' (internal, `network_layers.R`; built from `network_edges`,
#' `network_enrichment`, and `targets_disease`): every compound, its
#' targets, and (opt-in, see `layers`) whichever pathways/diseases those
#' targets reach, colored by layer, node size scaled by degree, and only
#' the highest-degree ("hub") nodes labeled -- the same visual language as
#' Yildirim et al. 2007's drug-target network figure (\doi{10.1038/nbt1338},
#' the reference Uriel shared, 2026-07-23): one big connected "hairball"
#' plus many small disconnected components, categorical node colors with a
#' fixed legend, hubs called out by name.
#'
#' @section `condition = NULL` pools every built condition into one graph:
#' Yildirim et al.'s clustering-by-category look is a property of the
#' *global* drug-target network (890 drugs at once) -- a single condition's
#' handful of compounds cannot show that kind of structure no matter how
#' the layout is tuned, there just are not enough shared targets yet.
#' `condition = NULL` (default) pools every condition [network_build()] has
#' built into one combined graph (same compound/target/pathway/disease
#' node, if it appears in more than one condition, becomes one node with
#' edges from all of them -- not duplicated); pass one or more condition
#' names to restrict to those specifically, e.g. `condition = "FLO-ET"`
#' for the original per-condition view.
#'
#' @section Layers are opt-in beyond compound/target, on purpose:
#' Default `layers = c("compound", "target")`, **not** `"pathway"` or
#' `"disease"`. Two independent reasons, both confirmed against real
#' renders (2026-07-23): `network_enrich()` can return hundreds of
#' significant GO/Reactome terms, each attaching to its target as a
#' degree-1 leaf, which turns the plot into an unreadable dandelion (see
#' `max_pathways` for how to include pathways without that happening); and
#' Yildirim et al.'s own reference figure does not draw disease genes as
#' graph nodes at all -- they measure drug-target-to-disease-gene distance
#' in the interactome as a *separate* analysis (exactly what
#' [network_proximity()] already does in this package), not as part of
#' this structural network. Both layers stay available for anyone who
#' wants them anyway -- just not the default.
#'
#' @section Labels are looked up, not raw IDs, when possible:
#' Compound nodes are labeled with [compounds()]`$name` (falling back to
#' `compound_id` if missing). Target nodes are labeled with their gene
#' symbol via `clusterProfiler::bitr()` + `org.Hs.eg.db` (same optional
#' `UNIPROT -> ...` mapping [network_enrich()]/`.network_target_pathway_edges()`
#' already use elsewhere in this family; falls back to the raw UniProt ID
#' if either package is not installed or a given ID does not map).
#' Pathway nodes are labeled with `network_enrichment$Description` when
#' available (falls back to the raw pathway ID). Disease nodes are labeled
#' with the raw EFO ID -- patliR does not currently store a resolved
#' disease display name anywhere, so there is nothing better to show yet.
#'
#' @section Derived edges are drawn lighter and dashed:
#' `edge_kind %in% c("compound_pathway_derived", "compound_disease_derived")`
#' (see `.network_layered_graph()`) are transitive relationships inferred
#' through a target, not a directly recorded interaction -- drawn dashed
#' and more transparent so they read as "consequence of", not "same
#' evidentiary weight as", the direct edges (`compound_target`,
#' `target_pathway`, `target_disease`).
#'
#' @param proj A `PatliRProject` object, with [network_build()] already run.
#' @param condition `NULL` (default -- pool every condition
#'   [network_build()] has built into one combined graph, see the section
#'   above) or a character vector of one or more condition names to
#'   restrict to.
#' @param engine `"static"` (default, plain `ggplot2`) or `"ggiraph"`
#'   (adds a tooltip per node with its id/layer/degree; falls back to
#'   `"static"` with a warning if `ggiraph` is not installed).
#' @param layout One of `"fr"` (default, Fruchterman-Reingold -- the
#'   classic "hairball" force layout, matching the reference figure),
#'   `"kk"` (Kamada-Kawai), or `"drl"`; matched against
#'   `igraph::layout_with_<layout>()`.
#' @param top_hub_n Integer, default `15`. Only the `top_hub_n`
#'   highest-degree nodes get a text label (matches the reference figure:
#'   a handful of named hubs, not every node).
#' @param seed Single integer or `NULL` (default `1`). Force-directed
#'   layouts are stochastic; a fixed seed makes the figure reproducible
#'   across reruns of the same graph. Pass `NULL` to not fix it.
#' @param layers Which layers to draw, subset of `c("compound", "target",
#'   "pathway", "disease")`. Default `c("compound", "target")` -- see the
#'   section above for why `"pathway"`/`"disease"` are opt-in.
#' @param max_pathways Integer, default `30`. Only used when `"pathway"`
#'   is in `layers`: keeps the `max_pathways` most significant pathway
#'   nodes (lowest `network_enrichment$p.adjust`, pooled across whichever
#'   conditions are in scope) and drops the rest, along with any
#'   `target_pathway`/`compound_pathway_derived` edge that only reached a
#'   dropped pathway. Set to `Inf` to keep every significant pathway (not
#'   recommended -- see the section above for why).
#' @param pathway_db See `.network_target_pathway_edges()`, internal --
#'   restrict the pathway layer to specific `network_enrich()` `db`
#'   value(s) (e.g. `"kegg"` alone) before `max_pathways` even applies, or
#'   `NULL` (default) for every `db` combined.
#' @param save Logical, default `TRUE`. If `TRUE`, also writes a PNG to
#'   `out_dir`.
#' @param out_dir Directory to write the PNG to (only used if
#'   `save = TRUE`). Defaults to `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] for the saved PNG.
#'
#' @return A `ggplot` object (`engine = "static"`) or a `girafe` htmlwidget
#'   (`engine = "ggiraph"`). If `save = TRUE` (the default), the PNG path
#'   is also recorded in `patliRResults(proj, "network_layers_plot_log")`
#'   (one row per distinct `condition` scope, `"ALL"` for the pooled
#'   default -- upserted, see `.network_upsert()`), retrievable via
#'   `attr(result, "proj")`.
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
#' plot_network_layers(proj, save = FALSE) # pooled across every condition
#' plot_network_layers(proj, condition = "FLO-ET", save = FALSE) # one condition
#'
#' ## To also draw the pathway layer (opt-in -- see the layers section
#' ## above for why), run network_enrich() first and pass it explicitly,
#' ## keeping max_pathways capped:
#' # proj <- network_enrich(proj, condition = "FLO-ET", db = "go")
#' # plot_network_layers(
#' #   proj, condition = "FLO-ET",
#' #   layers = c("compound", "target", "pathway"), max_pathways = 30,
#' #   save = FALSE
#' # )
#' }
#'
#' @export
plot_network_layers <- function(proj, condition = NULL, engine = c("static", "ggiraph"),
                                 layout = c("fr", "kk", "drl"), top_hub_n = 15,
                                 seed = 1, layers = c("compound", "target"),
                                 max_pathways = 30, pathway_db = NULL, save = TRUE, out_dir = NULL,
                                 width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  layout <- match.arg(layout)
  unknown_layers <- setdiff(layers, c("compound", "target", "pathway", "disease"))
  if (length(unknown_layers) > 0) {
    cli::cli_abort("Unknown {.arg layers} value(s) {.val {unknown_layers}}; must be a subset of {.val {c('compound', 'target', 'pathway', 'disease')}}.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_network_layers}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  conditions <- .network_resolve_conditions(proj, condition) # NULL -> every built condition
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  g <- .network_layered_graph_multi(proj, conditions, pathway_db = pathway_db)
  g <- .network_layers_filter(proj, g, conditions, layers = layers, max_pathways = max_pathways)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort("Condition(s) {.val {conditions}} have zero nodes for {.arg layers} = {.val {layers}} -- nothing to plot.")
  }

  pd <- .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed)
  p <- .network_layers_ggplot(pd, engine = engine, title_suffix = scope_label)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    p_static <- if (engine == "static") p else .network_layers_ggplot(pd, engine = "static", title_suffix = scope_label)
    path <- file.path(out_dir, paste0("network_layers_", scope_label, ".png"))
    ggplot2::ggsave(path, p_static, width = width, height = height, dpi = dpi)
    log_row <- data.frame(
      condition = scope_label, path = path,
      n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g),
      stringsAsFactors = FALSE
    )
    log_df <- .network_upsert(proj, "network_layers_plot_log", log_row, "condition")
    patliRResults(proj, "network_layers_plot_log") <- log_df
    .write_results_csv(proj, "network_layers_plot_log", log_df)
  }

  result <- if (engine == "static") {
    p
  } else {
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  }
  if (save) attr(result, "proj") <- proj
  result
}

#' Build the directed compound-target-pathway-disease layered graph pooled
#' across one or more conditions
#'
#' @description
#' Same edge-construction logic as `.network_layered_graph()`
#' (`network_layers.R`), generalized to pool `ct`/pathway/disease rows
#' across every condition in `conditions` before building a single graph,
#' instead of that function's one-condition-at-a-time contract (kept
#' as-is there since [network_motifs()] genuinely needs per-condition
#' motif counts -- this is a deliberate sibling, not a replacement).
#' A node/edge that appears in more than one condition is not duplicated:
#' `unique()` on the pooled tables collapses it to one.
#'
#' @param pathway_db See `.network_target_pathway_edges()`, internal.
#' @return An `igraph` directed graph, same vertex/edge attribute shape as
#'   `.network_layered_graph()`.
#' @keywords internal
.network_layered_graph_multi <- function(proj, conditions, pathway_db = NULL) {
  edges_all <- patliRResults(proj, "network_edges")
  ct <- unique(edges_all[edges_all$condition %in% conditions, c("compound_id", "uniprot_id")])
  edges <- data.frame(
    from = ct$compound_id, to = ct$uniprot_id,
    from_layer = "compound", to_layer = "target", edge_kind = "compound_target",
    stringsAsFactors = FALSE
  )

  tp_list <- lapply(conditions, function(cond) .network_target_pathway_edges(proj, cond, unique(ct$uniprot_id), pathway_db = pathway_db))
  tp_list <- tp_list[!vapply(tp_list, is.null, logical(1))]
  target_pathway <- if (length(tp_list) > 0) unique(do.call(rbind, tp_list)) else NULL
  if (!is.null(target_pathway) && nrow(target_pathway) > 0) {
    edges <- rbind(edges, data.frame(
      from = target_pathway$uniprot_id, to = target_pathway$pathway_id,
      from_layer = "target", to_layer = "pathway", edge_kind = "target_pathway",
      stringsAsFactors = FALSE
    ))
    cp <- unique(merge(ct, target_pathway, by = "uniprot_id")[, c("compound_id", "pathway_id")])
    if (nrow(cp) > 0) {
      edges <- rbind(edges, data.frame(
        from = cp$compound_id, to = cp$pathway_id,
        from_layer = "compound", to_layer = "pathway", edge_kind = "compound_pathway_derived",
        stringsAsFactors = FALSE
      ))
    }
  }

  disease_all <- patliRResults(proj, "targets_disease")
  if (!is.null(disease_all) && nrow(disease_all) > 0) {
    td <- unique(disease_all[disease_all$target_id %in% ct$uniprot_id, c("target_id", "disease_id")])
    if (nrow(td) > 0) {
      edges <- rbind(edges, data.frame(
        from = td$target_id, to = td$disease_id,
        from_layer = "target", to_layer = "disease", edge_kind = "target_disease",
        stringsAsFactors = FALSE
      ))
      cd <- unique(merge(ct, td, by.x = "uniprot_id", by.y = "target_id")[, c("compound_id", "disease_id")])
      if (nrow(cd) > 0) {
        edges <- rbind(edges, data.frame(
          from = cd$compound_id, to = cd$disease_id,
          from_layer = "compound", to_layer = "disease", edge_kind = "compound_disease_derived",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  vertices <- unique(rbind(
    data.frame(name = edges$from, layer = edges$from_layer, stringsAsFactors = FALSE),
    data.frame(name = edges$to, layer = edges$to_layer, stringsAsFactors = FALSE)
  ))
  igraph::graph_from_data_frame(d = edges[, c("from", "to", "edge_kind")], directed = TRUE, vertices = vertices)
}

#' Restrict a layered graph to the requested `layers`, capping pathway
#' nodes to the most significant `max_pathways` when included
#'
#' @return An `igraph` induced subgraph (vertex/edge attributes preserved).
#' @keywords internal
.network_layers_filter <- function(proj, g, conditions, layers, max_pathways) {
  vlayer <- igraph::V(g)$layer
  keep <- vlayer %in% layers

  if ("pathway" %in% layers && is.finite(max_pathways)) {
    pathway_names <- igraph::V(g)$name[vlayer == "pathway"]
    if (length(pathway_names) > max_pathways) {
      enrichment_all <- patliRResults(proj, "network_enrichment")
      rank_col <- if (!is.null(enrichment_all) && "p.adjust" %in% names(enrichment_all)) "p.adjust" else NULL
      keep_pathways <- pathway_names
      if (!is.null(rank_col)) {
        enr <- enrichment_all[enrichment_all$condition %in% conditions & enrichment_all$ID %in% pathway_names, , drop = FALSE]
        enr <- enr[order(enr[[rank_col]]), , drop = FALSE]
        keep_pathways <- utils::head(unique(enr$ID), max_pathways)
      } else {
        keep_pathways <- utils::head(pathway_names, max_pathways)
      }
      keep <- keep & !(vlayer == "pathway" & !(igraph::V(g)$name %in% keep_pathways))
    }
  }

  igraph::induced_subgraph(g, igraph::V(g)[keep])
}

#' Compute a force-directed layout + labels/colors/sizes for
#' `plot_network_layers()`
#'
#' @return `list(nodes = data.frame(name, layer, label, x, y, degree,
#'   is_hub), edges = data.frame(x, y, xend, yend, edge_kind, is_derived))`.
#' @keywords internal
.network_layers_plot_data <- function(proj, g, conditions, layout, top_hub_n, seed) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  ## niter well above igraph's default (500) -- with the pathway layer
  ## capped (or excluded) this is a few hundred nodes at most, cheap to
  ## run longer, and a half-settled FR layout is exactly what produced
  ## the "everything crushed into the middle, degree-1 leaves splayed out
  ## in a starburst" look on Uriel's first render (2026-07-23) instead of
  ## the loosely clustered look of the reference figure.
  coords <- switch(layout,
    fr = igraph::layout_with_fr(g, niter = 2000),
    kk = igraph::layout_with_kk(g, maxiter = 2000),
    drl = igraph::layout_with_drl(g)
  )

  vnames <- igraph::V(g)$name
  vlayer <- igraph::V(g)$layer
  degree <- igraph::degree(g, mode = "all")

  nodes <- data.frame(
    name = vnames, layer = vlayer, degree = degree,
    x = coords[, 1], y = coords[, 2], stringsAsFactors = FALSE
  )
  nodes$label <- .network_layers_labels(proj, conditions, nodes)
  hub_cut <- if (nrow(nodes) <= top_hub_n) -Inf else sort(nodes$degree, decreasing = TRUE)[top_hub_n]
  nodes$is_hub <- nodes$degree >= hub_cut & nodes$degree > 0

  ed <- igraph::as_data_frame(g, what = "edges")
  idx <- stats::setNames(seq_len(nrow(nodes)), nodes$name)
  edges <- data.frame(
    x = nodes$x[idx[ed$from]], y = nodes$y[idx[ed$from]],
    xend = nodes$x[idx[ed$to]], yend = nodes$y[idx[ed$to]],
    edge_kind = ed$edge_kind, stringsAsFactors = FALSE
  )
  edges$is_derived <- edges$edge_kind %in% c("compound_pathway_derived", "compound_disease_derived")

  list(nodes = nodes, edges = edges)
}

#' Look up display labels for every node of a (pooled) layered graph
#' @return Character vector, same length/order as `nodes$name`.
#' @keywords internal
.network_layers_labels <- function(proj, conditions, nodes) {
  label <- nodes$name

  is_compound <- nodes$layer == "compound"
  if (any(is_compound)) {
    cmp <- compounds(proj)
    lookup <- stats::setNames(cmp$name, cmp$id)
    resolved <- lookup[nodes$name[is_compound]]
    label[is_compound] <- ifelse(is.na(resolved) | resolved == "", nodes$name[is_compound], resolved)
  }

  is_target <- nodes$layer == "target"
  if (any(is_target) && requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    map <- tryCatch(
      clusterProfiler::bitr(nodes$name[is_target], fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = TRUE),
      error = function(e) NULL
    )
    if (!is.null(map) && nrow(map) > 0) {
      lookup <- stats::setNames(map$SYMBOL, map$UNIPROT)
      resolved <- lookup[nodes$name[is_target]]
      label[is_target] <- ifelse(is.na(resolved), nodes$name[is_target], resolved)
    }
  }

  is_pathway <- nodes$layer == "pathway"
  if (any(is_pathway)) {
    enrichment_all <- patliRResults(proj, "network_enrichment")
    if (!is.null(enrichment_all)) {
      enr <- enrichment_all[enrichment_all$condition %in% conditions, , drop = FALSE]
      lookup <- stats::setNames(enr$Description, enr$ID)
      resolved <- lookup[nodes$name[is_pathway]]
      label[is_pathway] <- ifelse(is.na(resolved) | resolved == "", nodes$name[is_pathway], resolved)
    }
  }

  label
}

#' @keywords internal
.network_layers_layer_colors <- function() {
  c(compound = "#e67e22", target = "#2980b9", pathway = "#27ae60", disease = "#c0392b")
}

#' @keywords internal
.network_layers_ggplot <- function(pd, engine, title_suffix) {
  nodes <- pd$nodes
  edges <- pd$edges
  ## sqrt-scaled point size: hub nodes read as visibly bigger without a
  ## couple of very high-degree nodes (compounds especially, if they map
  ## to many targets) blowing every other node down to invisible.
  nodes$point_size <- 1.2 + 2.6 * sqrt(nodes$degree / max(nodes$degree, 1))
  nodes$tooltip <- sprintf("%s\n%s, degree %d", nodes$label, nodes$layer, nodes$degree)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                    linetype = .data$is_derived, alpha = .data$is_derived),
      colour = "grey75", linewidth = 0.15
    ) +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.2), guide = "none")

  if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p <- p + ggiraph::geom_point_interactive(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$layer, size = .data$point_size,
                    tooltip = .data$tooltip, data_id = .data$name)
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$layer, size = .data$point_size)
    )
  }

  hubs <- nodes[nodes$is_hub, , drop = FALSE]
  if (nrow(hubs) > 0) {
    p <- p + ggplot2::geom_text(
      data = hubs,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 3, fontface = "bold", colour = "grey15", vjust = -1, show.legend = FALSE
    )
  }

  layer_title <- paste(unique(nodes$layer[order(match(nodes$layer, c("compound", "target", "pathway", "disease")))]), collapse = "-")
  p +
    ggplot2::scale_colour_manual(values = .network_layers_layer_colors(), name = "Layer") +
    ggplot2::scale_size_identity(guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(layer_title, " network -- ", title_suffix),
      subtitle = "Node size = degree; labeled nodes are the top-degree hubs; dashed edges are derived (transitive) relationships",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      legend.position = "right"
    )
}
