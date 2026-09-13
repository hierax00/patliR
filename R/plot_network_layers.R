#' @include AllGenerics.R internal.R network_layers.R plot-helpers.R
NULL

## Static ggplot2 figure (no ggraph): the layout is computed with
## igraph::layout_with_*() and drawn as plain node/edge data.frames, same
## "minimal Suggests" principle as the rest of plot_*.
##
## `layers` defaults to c("compound", "target") only -- network_enrich()
## returns hundreds of GO/Reactome terms, each a degree-1 leaf that drowns
## the compound-target structure, so "pathway"/"disease" are opt-in and
## `max_pathways` caps them. The category-clustering look only emerges when
## many compounds share targets, so `condition = NULL` pools every built
## condition into one graph (a single small condition cannot look like
## that).

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
#' Yildirim et al. 2007's drug-target network figure (\doi{10.1038/nbt1338}):
#' one big connected "hairball" plus many small disconnected components,
#' categorical node colors with a fixed legend, hubs called out by name.
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
#' `"disease"`. Two independent reasons: `network_enrich()` can return
#' hundreds of significant GO/Reactome terms, each attaching to its target
#' as a degree-1 leaf, which turns the plot into an unreadable dandelion
#' (see `max_pathways` for how to include pathways without that happening);
#' and Yildirim et al.'s own reference figure does not draw disease genes as
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
#'   `"kk"` (Kamada-Kawai), `"drl"` (matched against
#'   `igraph::layout_with_<layout>()`), or `"bipartite"` -- a deterministic
#'   two-column layout (compounds on the left, targets on the right, each
#'   column ordered by degree so fewer edges cross; the standard
#'   affiliation-network rendering, Borgatti & Everett 1997). `"bipartite"`
#'   only works when the scope's graph is genuinely two-layer
#'   (compound-target only, the default `layers`) -- `cli_abort`s naming
#'   the extra layer(s) otherwise.
#' @param colour_by One of `"layer"` (default -- the historical behaviour:
#'   compound/target/pathway/disease, `.network_layers_layer_colors()`),
#'   `"module"` (joins [network_module_robustness()]'s per-node
#'   `network_module_membership` table by `(condition, node_id)`, coloured
#'   with a qualitative palette (`grDevices::hcl.colors(n, "Dark 3")`) and a
#'   grey `"not clustered"` level for any node absent from that table --
#'   e.g. `network_module_robustness()` was never run; `cli_abort`s if the
#'   table does not exist at all, naming the function to run), or
#'   `"node_type"` (compound vs. target, two colours; any opted-in
#'   pathway/disease node groups as `"other"`).
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
#'   (one row per distinct `(condition, layout, colour_by)` combination --
#'   `"ALL"` for the pooled default `condition` -- upserted, see
#'   `.network_upsert()`; the filename encodes all three too, e.g.
#'   `network_layers_ALL_fr_layer.png`, so different `layout`/`colour_by`
#'   views never overwrite each other's log row or PNG), retrievable via
#'   `attr(result, "proj")`.
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
#' plot_network_layers(proj, save = FALSE) # pooled across every condition
#' plot_network_layers(proj, condition = "FLO-ET", save = FALSE) # one condition
#'
#' ## Two-column bipartite layout (compound-target only):
#' plot_network_layers(proj, condition = "FLO-ET", layout = "bipartite", save = FALSE)
#'
#' ## Colour by module (needs network_module_robustness() first):
#' # proj <- network_module_robustness(proj, condition = "FLO-ET")
#' # plot_network_layers(proj, condition = "FLO-ET", colour_by = "module", save = FALSE)
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
                                 layout = c("fr", "kk", "drl", "bipartite"),
                                 colour_by = c("layer", "module", "node_type"), top_hub_n = 15,
                                 seed = 1, layers = c("compound", "target"),
                                 max_pathways = 30, pathway_db = NULL, save = TRUE, out_dir = NULL,
                                 width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  layout <- match.arg(layout)
  colour_by <- match.arg(colour_by)
  unknown_layers <- setdiff(layers, c("compound", "target", "pathway", "disease"))
  if (length(unknown_layers) > 0) {
    cli::cli_abort("Unknown {.arg layers} value(s) {.val {unknown_layers}}; must be a subset of {.val {c('compound', 'target', 'pathway', 'disease')}}.")
  }
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition) # condition NULL -> every built condition
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  g <- .network_layered_graph_multi(proj, conditions, pathway_db = pathway_db)
  g <- .network_layers_filter(proj, g, conditions, layers = layers, max_pathways = max_pathways)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort("Condition(s) {.val {conditions}} have zero nodes for {.arg layers} = {.val {layers}} -- nothing to plot.")
  }

  pd <- .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed, colour_by = colour_by)
  p <- .network_layers_ggplot(pd, engine = engine, title_suffix = scope_label)

  .plot_finish(
    proj, p,
    name = "network_layers_plot_log",
    filename = paste0("network_layers_", scope_label, "_", layout, "_", colour_by, ".png"),
    log_row = data.frame(
      condition = scope_label, layout = layout, colour_by = colour_by, path = NA_character_,
      n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g),
      stringsAsFactors = FALSE
    ),
    key_cols = c("condition", "layout", "colour_by"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi,
    static = if (engine == "static") p else .network_layers_ggplot(pd, engine = "static", title_suffix = scope_label)
  )
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

#' Compute a layout + labels/colors/sizes for `plot_network_layers()`
#'
#' @param colour_by `"layer"` (default, historical behaviour), `"module"`,
#'   or `"node_type"` -- see [plot_network_layers()]'s own `colour_by` docs.
#' @return `list(nodes = data.frame(name, layer, label, x, y, degree,
#'   is_hub, colour_group), edges = data.frame(x, y, xend, yend, edge_kind,
#'   is_derived), palette = named colour vector, legend_name =
#'   character(1))`.
#' @keywords internal
.network_layers_plot_data <- function(proj, g, conditions, layout, top_hub_n, seed, colour_by = "layer") {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  ## niter well above igraph's default (500) -- with the pathway layer
  ## capped (or excluded) this is a few hundred nodes at most, cheap to
  ## run longer, and a half-settled FR layout is exactly what produced
  ## the "everything crushed into the middle, degree-1 leaves splayed out
  ## in a starburst" look instead of
  ## the loosely clustered look of the reference figure.
  coords <- switch(layout,
    fr = igraph::layout_with_fr(g, niter = 2000),
    kk = igraph::layout_with_kk(g, maxiter = 2000),
    drl = igraph::layout_with_drl(g),
    bipartite = .network_layers_bipartite_coords(g)
  )

  vnames <- igraph::V(g)$name
  vlayer <- igraph::V(g)$layer
  degree <- igraph::degree(g, mode = "all")

  nodes <- data.frame(
    name = vnames, layer = vlayer, degree = degree,
    x = coords[, 1], y = coords[, 2], stringsAsFactors = FALSE
  )
  nodes$label <- .plot_label_nodes(proj, conditions, nodes$name, nodes$layer)
  hub_cut <- if (nrow(nodes) <= top_hub_n) -Inf else sort(nodes$degree, decreasing = TRUE)[top_hub_n]
  nodes$is_hub <- nodes$degree >= hub_cut & nodes$degree > 0

  colour_info <- .network_layers_colour_info(proj, conditions, nodes, colour_by)
  nodes$colour_group <- colour_info$group

  ed <- igraph::as_data_frame(g, what = "edges")
  idx <- stats::setNames(seq_len(nrow(nodes)), nodes$name)
  edges <- data.frame(
    x = nodes$x[idx[ed$from]], y = nodes$y[idx[ed$from]],
    xend = nodes$x[idx[ed$to]], yend = nodes$y[idx[ed$to]],
    edge_kind = ed$edge_kind, stringsAsFactors = FALSE
  )
  edges$is_derived <- edges$edge_kind %in% c("compound_pathway_derived", "compound_disease_derived")

  list(nodes = nodes, edges = edges, palette = colour_info$palette, legend_name = colour_info$legend_name)
}

#' Deterministic two-column bipartite layout for `layout = "bipartite"`:
#' compounds on the left, targets on the right, each column ordered by
#' degree so fewer edges cross
#'
#' @description
#' Only valid on a genuinely two-layer (compound-target) graph -- aborts,
#' naming any extra layer, otherwise (a pathway/disease node has no
#' left/right meaning in a two-column affiliation-network drawing).
#' `igraph::layout_as_bipartite()` is called to exercise its `types=`
#' contract, but the coordinates actually used are this function's own:
#' its crossing-minimizing layout lays the two modes out as top/bottom
#' *rows*, which reads top-to-bottom, not left-to-right -- the axes are
#' swapped here, and each resulting column is then ordered by degree
#' (rather than kept in the algorithm's own within-row order), which is
#' what keeps a small (tens-to-hundreds-of-nodes) bipartite figure legible.
#' @return A 2-column numeric matrix, `x` (0 = compound, 1 = target) and
#'   `y` (degree rank within the column), same row order as `igraph::V(g)`.
#' @keywords internal
.network_layers_bipartite_coords <- function(g) {
  vlayer <- igraph::V(g)$layer
  extra_layers <- setdiff(unique(vlayer), c("compound", "target"))
  if (length(extra_layers) > 0) {
    cli::cli_abort(c(
      "{.arg layout = \"bipartite\"} only supports a two-layer compound-target graph.",
      "i" = "This scope also includes {.val {extra_layers}} node(s) -- restrict {.arg layers} to {.val {c('compound', 'target')}} (the default) to use the bipartite layout."
    ))
  }
  is_target <- vlayer == "target"
  deg <- igraph::degree(g, mode = "all")

  invisible(igraph::layout_as_bipartite(g, types = is_target))

  x <- ifelse(is_target, 1, 0)
  y <- numeric(length(vlayer))
  for (side in c(FALSE, TRUE)) {
    idx <- which(is_target == side)
    if (length(idx) > 0) y[idx] <- rank(deg[idx], ties.method = "first")
  }
  cbind(x, y)
}

#' Resolve `plot_network_layers()`/`plot_network_degeneracy()`'s per-node
#' colour grouping + palette for `colour_by`
#'
#' @return `list(group = character vector, same length/order as
#'   `nodes$name`; palette = named colour vector; legend_name =
#'   character(1))`.
#' @keywords internal
.network_layers_colour_info <- function(proj, conditions, nodes, colour_by) {
  switch(
    colour_by,
    layer = list(
      group = nodes$layer,
      palette = .network_layers_layer_colors(),
      legend_name = "Layer"
    ),
    node_type = list(
      group = ifelse(nodes$layer %in% c("compound", "target"), nodes$layer, "other"),
      palette = c(compound = "#e67e22", target = "#2980b9", other = "grey70"),
      legend_name = "Node type"
    ),
    module = .network_layers_module_colour_info(proj, conditions, nodes$name)
  )
}

#' Join [network_module_robustness()]'s per-node `network_module_membership`
#' table onto a node table for `colour_by = "module"`, with a
#' `"not clustered"` grey fallback level
#'
#' @description
#' Shared by [plot_network_layers()] and [plot_network_degeneracy()] (the
#' latter draws its base layout via the same `.network_layers_plot_data()`
#' pipeline, so a degeneracy pair's arc endpoints are coloured by the same
#' module assignments -- letting the reader see whether a degenerate pair
#' crosses a module boundary). Nodes present in the scope's graph but
#' absent from `network_module_membership` (e.g.
#' `network_module_robustness()` was never run, or the membership table
#' simply never covers a pathway/disease node -- it is only ever built
#' from the compound-target graph) fall into a grey `"not clustered"`
#' level rather than an `NA` colour. Pooled multi-condition scope: a node
#' built into more than one condition can carry a different `module_id`
#' per condition -- the first match (by row order in the table) wins; the
#' common case is a single condition in scope, where this is an exact
#' join.
#' @return `list(group, palette, legend_name = "Module")`, same shape as
#'   `.network_layers_colour_info()`.
#' @keywords internal
.network_layers_module_colour_info <- function(proj, conditions, node_names) {
  membership_all <- patliRResults(proj, "network_module_membership")
  if (is.null(membership_all) || nrow(membership_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_module_membership} entry in {.arg proj}.",
      "i" = "{.code colour_by = \"module\"} needs per-node module assignments -- run {.fn network_module_robustness} first."
    ))
  }
  mem <- membership_all[membership_all$condition %in% conditions, , drop = FALSE]
  mem <- mem[!duplicated(mem$node_id), , drop = FALSE]
  lookup <- stats::setNames(mem$module_id, mem$node_id)
  module_id <- unname(lookup[node_names])
  group <- ifelse(is.na(module_id), "not clustered", module_id)

  mod_levels <- sort(unique(group[group != "not clustered"]))
  palette <- stats::setNames(grDevices::hcl.colors(length(mod_levels), "Dark 3"), mod_levels)
  palette["not clustered"] <- "grey70"

  list(group = group, palette = palette, legend_name = "Module")
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
  nodes$tooltip <- sprintf("%s\n%s, degree %d\n%s: %s", nodes$label, nodes$layer, nodes$degree, pd$legend_name, nodes$colour_group)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                    linetype = .data$is_derived, alpha = .data$is_derived),
      colour = "grey75", linewidth = 0.15
    ) +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.2), guide = "none")

  ## Node grouping is carried on `fill` (with a fixed grey border), not
  ## `colour`, so an overlay that maps `colour` to a *continuous* value --
  ## plot_network_degeneracy()'s degeneracy arcs, coloured by z_score --
  ## gets its own independent scale instead of clashing with this discrete
  ## one (ggplot2 supports only one scale per aesthetic).
  if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    p <- p + ggiraph::geom_point_interactive(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$colour_group, size = .data$point_size,
                    tooltip = .data$tooltip, data_id = .data$name),
      shape = 21, colour = "grey30", stroke = 0.2
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$colour_group, size = .data$point_size),
      shape = 21, colour = "grey30", stroke = 0.2
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
    ggplot2::scale_fill_manual(values = pd$palette, name = pd$legend_name) +
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
