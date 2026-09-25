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
#' @section Splitting a large network into readable figures:
#' A condition with ~20 compounds and several hundred predicted targets,
#' most of them hit by a single compound, draws as a hairball whatever the
#' layout. Three optional filters (all off by default) cut it down:
#' `min_target_degree = 2` keeps only shared targets; `compound_ids` draws
#' one or a few compounds' sub-network; `module_id` draws one module of
#' [network_module_robustness()]. Every filtered figure gets its own PNG
#' and log row (the filter is encoded in the file name and in the log's
#' `subset` column), so one call per module or per compound builds a
#' series of small figures -- see the examples.
#'
#' @section Views -- readable figures for large networks (`view`):
#' Real conditions (~15-25 compounds, 400+ predicted targets, ~2000
#' edges, most targets shared by several compounds) cannot be drawn as a
#' readable force-directed graph. `view` picks one of five figures, each
#' answering one question:
#' * `"network"` -- the historical force-directed (or `layout =
#'   "bipartite"`) drawing of every node; `layout`, `colour_by`, `layers`
#'   and `max_pathways` apply. Best for small graphs or one module /
#'   compound (`module_id`, `compound_ids`).
#' * `"compound"` -- compound-centric summary. Compounds sit on a ring,
#'   ordered by the similarity (Jaccard) of their target sets. Only targets
#'   hit by at least `k` compounds are drawn, on concentric shells (the more
#'   compounds share a target, the closer to the centre, the larger and the
#'   darker it is); the other targets are collapsed into a per-compound
#'   badge (`"190 targets, +9 unique"`: total targets and targets no other
#'   compound hits). `k = min_shared`, or, when `NULL`, the smallest
#'   `k >= 2` that leaves at most `max_targets` targets. Edge width =
#'   prediction probability; the `top_hub_n` most shared targets are
#'   labelled. `colour_by = "module"` colours targets by module.
#' * `"bipartite"` -- compounds in a left column, targets in a right column,
#'   targets grouped by [network_module_robustness()] module when available
#'   (coloured band; compounds grouped and coloured by their own module)
#'   and sorted by how many compounds hit them; compounds ordered by the
#'   mean position of their targets (barycentric crossing reduction).
#'   Edges are nearly transparent so bundles carry the message. All
#'   targets are drawn unless `min_shared` is set.
#' * `"module"` -- one shaded region (padded convex hull) per module of
#'   [network_module_robustness()] (required), each laid out on its own
#'   with the module's compounds pinned on an inner ring; edges within a
#'   module in its colour, edges between modules faint grey; every
#'   compound labelled plus the targets with the highest
#'   `network_centrality` hub score in each module. Best per condition:
#'   with pooled conditions a node takes the module of the first condition
#'   listing it.
#' * `"pathway"` -- compound x enriched-term dot plot: the `max_pathways`
#'   most significant [network_enrich()] terms (restrict with
#'   `pathway_db`, e.g. `"kegg"`), dot size = number of the compound's
#'   targets annotated to the term, colour = share of the term's in-scope
#'   targets the compound hits; rows and columns clustered. Needs
#'   `clusterProfiler` + `org.Hs.eg.db`; always static.
#'
#' `view = "auto"` (default) keeps the `"network"` figure for graphs of at
#' most `auto_max_nodes` nodes (and whenever `layout` is given explicitly
#' or `layers` include pathways/diseases) and switches to `"compound"`
#' above that, saying so in the subtitle. The filters (`condition`,
#' `compound_ids`, `module_id`, `min_target_degree`) apply before every
#' view. In the simplified views `top_hub_n` is the number of target
#' labels (the most shared targets, spread across modules in the
#' `"bipartite"`/`"module"` views) and `colour_by` only matters for
#' `"compound"` (`"module"` colours the targets by module, anything else
#' by how many compounds share them). Colours are colour-blind safe (Okabe
#' & Ito, a module keeps its colour in every view; viridis for the
#' sequential scales); labels repel when `ggrepel` is installed.
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
#'   with the colour-blind-safe Okabe-Ito palette (`grDevices::hcl.colors(n,
#'   "Dark 3")` beyond eight modules) and a
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
#'   recommended -- see the section above for why). With `view =
#'   "pathway"`: the number of terms (rows) drawn.
#' @param pathway_db See `.network_target_pathway_edges()`, internal --
#'   restrict the pathway layer to specific `network_enrich()` `db`
#'   value(s) (e.g. `"kegg"` alone) before `max_pathways` even applies, or
#'   `NULL` (default) for every `db` combined. Also restricts the terms of
#'   `view = "pathway"`.
#' @param min_target_degree Integer >= 1, default `1` (every target, the
#'   historical behaviour). Drops targets hit by fewer than this many
#'   distinct compounds -- `2` keeps only the targets shared between
#'   compounds, which removes the degree-1 "dandelion" leaves that make a
#'   condition with a few hundred singleton targets unreadable. See the
#'   section on splitting large networks.
#' @param compound_ids `NULL` (default, every compound) or a character
#'   vector of `compound_id`s: draw only the sub-network of those compounds
#'   and their targets.
#' @param module_id `NULL` (default) or one or more `module_id`s from
#'   [network_module_robustness()]'s `network_module_membership`: draw only
#'   the compound/target nodes of those module(s).
#' @param view One of `"auto"` (default), `"network"`, `"compound"`,
#'   `"bipartite"`, `"module"` or `"pathway"` -- see the section on views.
#' @param min_shared `NULL` (default) or a single integer >= 1: the
#'   minimum number of compounds a target must be hit by to be drawn in the
#'   `"compound"`, `"bipartite"` and `"module"` views. `NULL` means
#'   automatic for `"compound"` (see `max_targets`) and every target for
#'   the other two.
#' @param max_targets Integer, default `40`: target budget of the
#'   `"compound"` view when `min_shared = NULL` -- the smallest sharing
#'   threshold `k >= 2` that leaves at most this many targets is used.
#' @param auto_max_nodes Number, default `300`: `view = "auto"` switches to
#'   the `"compound"` view above this many nodes (`Inf` never switches).
#' @param save Logical, default `TRUE`. If `TRUE`, also writes a PNG to
#'   `out_dir`.
#' @param out_dir Directory to write the PNG to (only used if
#'   `save = TRUE`). Defaults to `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()] for the saved PNG.
#'
#' @return A `ggplot` object (`engine = "static"`) or a `girafe` htmlwidget
#'   (`engine = "ggiraph"`). If `save = TRUE` (the default), the PNG path
#'   is also recorded in `patliRResults(proj, "network_layers_plot_log")`
#'   (one row per distinct `(condition, layout, colour_by, subset, view)`
#'   combination -- `"ALL"` for the pooled default `condition` -- upserted,
#'   see `.network_upsert()`; the filename encodes them too, e.g.
#'   `network_layers_ALL_fr_layer.png` for the `"network"` view and
#'   `network_layers_EFLO-S_compound_k12.png` for the others, so different
#'   views never overwrite each other's log row or PNG; `subset` is `""`
#'   unless `min_target_degree`/`compound_ids`/`module_id` filter the graph
#'   or the view adds its own threshold such as `"k12"`; `layout` is `NA`
#'   outside the `"network"` view), retrievable via `attr(result, "proj")`.
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
#' ## Readable views of a large condition (see the section on views):
#' plot_network_layers(proj, condition = "FLO-ET", view = "compound", save = FALSE)
#' plot_network_layers(proj, condition = "FLO-ET", view = "compound", min_shared = 3, save = FALSE)
#' plot_network_layers(proj, condition = "FLO-ET", view = "bipartite", save = FALSE)
#' # proj <- network_module_robustness(proj, condition = "FLO-ET")
#' # plot_network_layers(proj, condition = "FLO-ET", view = "module", save = FALSE)
#' # proj <- network_enrich(proj, condition = "FLO-ET", db = "kegg")
#' # plot_network_layers(proj, condition = "FLO-ET", view = "pathway", pathway_db = "kegg", save = FALSE)
#'
#' ## Colour by module (needs network_module_robustness() first):
#' # proj <- network_module_robustness(proj, condition = "FLO-ET")
#' # plot_network_layers(proj, condition = "FLO-ET", colour_by = "module", save = FALSE)
#'
#' ## Split a hairball: only targets shared by >= 2 compounds ...
#' plot_network_layers(proj, condition = "FLO-ET", min_target_degree = 2, save = FALSE)
#' ## ... or one figure per module (after network_module_robustness()) ...
#' # proj <- network_module_robustness(proj, condition = "FLO-ET")
#' # mem <- patliRResults(proj, "network_module_membership")
#' # for (m in sort(unique(stats::na.omit(mem$module_id[mem$condition == "FLO-ET"])))) {
#' #   proj <- attr(plot_network_layers(proj, condition = "FLO-ET", module_id = m,
#' #                                    colour_by = "module"), "proj")
#' # }
#' ## ... or one per compound:
#' # for (id in compounds(proj)$id) {
#' #   proj <- attr(plot_network_layers(proj, condition = "FLO-ET", compound_ids = id), "proj")
#' # }
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
                                 max_pathways = 30, pathway_db = NULL,
                                 min_target_degree = 1, compound_ids = NULL, module_id = NULL,
                                 view = c("auto", "network", "compound", "bipartite", "module", "pathway"),
                                 min_shared = NULL, max_targets = 40, auto_max_nodes = 300,
                                 save = TRUE, out_dir = NULL,
                                 width = 9, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  layout_given <- !missing(layout)
  layout <- match.arg(layout)
  colour_by <- match.arg(colour_by)
  view <- match.arg(view)
  unknown_layers <- setdiff(layers, c("compound", "target", "pathway", "disease"))
  if (length(unknown_layers) > 0) {
    cli::cli_abort("Unknown {.arg layers} value(s) {.val {unknown_layers}}; must be a subset of {.val {c('compound', 'target', 'pathway', 'disease')}}.")
  }
  if (!.network_layers_is_count(min_target_degree)) {
    cli::cli_abort("{.arg min_target_degree} must be a single integer >= 1.")
  }
  if (!is.null(min_shared) && !.network_layers_is_count(min_shared)) {
    cli::cli_abort("{.arg min_shared} must be {.code NULL} (automatic) or a single integer >= 1.")
  }
  if (!.network_layers_is_count(max_targets)) cli::cli_abort("{.arg max_targets} must be a single integer >= 1.")
  if (!is.numeric(top_hub_n) || length(top_hub_n) != 1L || !is.finite(top_hub_n) ||
      top_hub_n < 0 || top_hub_n != floor(top_hub_n)) {
    cli::cli_abort("{.arg top_hub_n} must be a non-negative integer.")
  }
  if (!is.numeric(auto_max_nodes) || length(auto_max_nodes) != 1L || is.na(auto_max_nodes) || auto_max_nodes < 1) {
    cli::cli_abort("{.arg auto_max_nodes} must be a single number >= 1 ({.code Inf} never simplifies).")
  }
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition) # condition NULL -> every built condition
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  g <- .network_layered_graph_multi(proj, conditions, pathway_db = pathway_db)
  g <- .network_layers_filter(proj, g, conditions, layers = layers, max_pathways = max_pathways)
  g <- .network_layers_subset(proj, g, conditions, compound_ids = compound_ids,
                              min_target_degree = min_target_degree, module_id = module_id)
  if (igraph::vcount(g) == 0) {
    cli::cli_abort(c(
      "Condition(s) {.val {conditions}} have zero nodes for {.arg layers} = {.val {layers}} -- nothing to plot.",
      "i" = if (!is.null(compound_ids) || !is.null(module_id) || min_target_degree > 1) "{.arg compound_ids}/{.arg module_id}/{.arg min_target_degree} removed every node; relax them."
    ))
  }
  subset_tag <- .network_layers_subset_tag(compound_ids, min_target_degree, module_id)
  title_suffix <- .network_layers_subset_title(scope_label, compound_ids, min_target_degree, module_id)

  ## `subset` and then `view` joined the log key after the first release --
  ## older logs get "" (the unfiltered figure) / "network" (the only view
  ## there was) so a default call still replaces its own row.
  proj <- .plot_log_backfill(proj, "network_layers_plot_log", "subset", "")
  proj <- .plot_log_backfill(proj, "network_layers_plot_log", "view", "network")
  key_cols <- c("condition", "layout", "colour_by", "subset", "view")

  resolved <- .network_layers_resolve_view(view, g, layout_given = layout_given, auto_max_nodes = auto_max_nodes)
  view <- resolved$view

  if (view == "network") {
    pd <- .network_layers_plot_data(proj, g, conditions, layout = layout, top_hub_n = top_hub_n, seed = seed, colour_by = colour_by)
    p <- .network_layers_ggplot(pd, engine = engine, title_suffix = title_suffix, fig_width = width)
    return(.plot_finish(
      proj, p,
      name = "network_layers_plot_log",
      filename = paste0("network_layers_", scope_label, "_", layout, "_", colour_by,
                        if (nzchar(subset_tag)) paste0("_", subset_tag), ".png"),
      log_row = data.frame(
        condition = scope_label, layout = layout, colour_by = colour_by, subset = subset_tag,
        view = "network", path = NA_character_, n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g),
        stringsAsFactors = FALSE
      ),
      key_cols = key_cols,
      engine = engine, save = save, out_dir = out_dir,
      width = width, height = height, dpi = dpi,
      static = if (engine == "static") p else .network_layers_ggplot(pd, engine = "static", title_suffix = title_suffix, fig_width = width)
    ))
  }

  ## Simplified views: built from the compound-target edges only.
  ct <- .network_layers_ct_table(proj, g, conditions)
  if (nrow(ct) == 0) {
    cli::cli_abort("No compound-target edge in scope -- {.code view = \"{view}\"} has nothing to draw.")
  }
  build <- function(eng) {
    switch(view,
      compound = .network_view_compound(proj, ct, conditions, min_shared = min_shared, max_targets = max_targets,
                                        top_hub_n = top_hub_n, colour_by = colour_by, seed = seed, engine = eng,
                                        title_suffix = title_suffix, note = resolved$note, fig_width = width),
      bipartite = .network_view_bipartite(proj, ct, conditions, min_shared = min_shared, top_hub_n = top_hub_n,
                                          engine = eng, title_suffix = title_suffix, fig_width = width),
      module = .network_view_module(proj, ct, conditions, min_shared = min_shared, top_hub_n = top_hub_n,
                                    seed = seed, engine = eng, title_suffix = title_suffix, fig_width = width),
      pathway = .network_view_pathway(proj, ct, conditions, pathway_db = pathway_db, max_pathways = max_pathways,
                                      title_suffix = title_suffix, fig_width = width)
    )
  }
  built <- build(engine)
  p <- built$plot
  full_tag <- paste(c(built$view_tag[nzchar(built$view_tag)], subset_tag[nzchar(subset_tag)]), collapse = "_")
  used_colour <- if (view == "compound") colour_by else NA_character_
  .plot_finish(
    proj, p,
    name = "network_layers_plot_log",
    filename = paste0("network_layers_", scope_label, "_", view,
                      if (!is.na(used_colour) && used_colour != "layer") paste0("_", used_colour),
                      if (nzchar(full_tag)) paste0("_", full_tag), ".png"),
    log_row = data.frame(
      condition = scope_label, layout = NA_character_, colour_by = used_colour, subset = full_tag,
      view = view, path = NA_character_, n_nodes = built$n_nodes, n_edges = built$n_edges,
      stringsAsFactors = FALSE
    ),
    key_cols = key_cols,
    engine = if (view == "pathway") "static" else engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi,
    static = if (engine == "static" || view == "pathway") p else build("static")$plot
  )
}

#' Is `x` a single whole number >= 1?
#' @keywords internal
.network_layers_is_count <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
}

#' Resolve `view = "auto"` to a concrete view
#'
#' @description
#' `"auto"` becomes `"compound"` (the compound-centric summary) when the
#' graph is a plain compound-target graph with more than `auto_max_nodes`
#' nodes and the caller did not ask for a specific force `layout`;
#' otherwise `"network"` (the historical force-directed figure). Any other
#' `view` is returned unchanged.
#' @return `list(view, note)` -- `note` is the subtitle sentence explaining
#'   an automatic switch, or `NULL`.
#' @keywords internal
.network_layers_resolve_view <- function(view, g, layout_given = FALSE, auto_max_nodes = 300) {
  if (view != "auto") return(list(view = view, note = NULL))
  n <- igraph::vcount(g)
  two_layer <- all(igraph::V(g)$layer %in% c("compound", "target"))
  if (two_layer && !layout_given && n > auto_max_nodes) {
    return(list(
      view = "compound",
      note = sprintf("Simplified view chosen automatically: the full graph has %d nodes (> %s); view = \"network\" draws all of them.",
                     n, format(auto_max_nodes))
    ))
  }
  list(view = "network", note = NULL)
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
  if (length(conditions) == 0) cli::cli_abort("{.arg condition} must select at least one condition.")
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

#' Split a hairball into a readable sub-network: restrict to some
#' compounds / one module, or to targets shared by several compounds
#'
#' @description
#' Applied after `.network_layers_filter()`. `compound_ids` drops every
#' other compound; `module_id` keeps only compound/target nodes that
#' [network_module_robustness()]'s `network_module_membership` assigns to
#' those module(s) in any condition in scope; `min_target_degree` drops
#' targets hit by fewer than that many distinct compounds (compound-target
#' edges only -- `2` keeps just the shared targets, which is what carries
#' the network's structure; singleton targets are the "dandelion" leaves).
#' Nodes left isolated by any of these are dropped too. Pathway/disease
#' nodes are kept while something still links to them (derived
#' compound-pathway/-disease edges included). With every argument at its
#' default the graph is returned unchanged.
#' @return An `igraph` induced subgraph.
#' @keywords internal
.network_layers_subset <- function(proj, g, conditions, compound_ids = NULL, min_target_degree = 1, module_id = NULL) {
  vname <- igraph::V(g)$name
  vlayer <- igraph::V(g)$layer
  drop <- rep(FALSE, length(vname))
  if (!is.null(compound_ids)) {
    drop <- drop | (vlayer == "compound" & !vname %in% compound_ids)
  }
  if (!is.null(module_id)) {
    membership_all <- patliRResults(proj, "network_module_membership")
    if (is.null(membership_all) || nrow(membership_all) == 0) {
      cli::cli_abort(c(
        "No {.val network_module_membership} entry in {.arg proj}.",
        "i" = "{.arg module_id} needs per-node module assignments -- run {.fn network_module_robustness} first."
      ))
    }
    mem <- membership_all[membership_all$condition %in% conditions & membership_all$module_id %in% module_id, , drop = FALSE]
    if (nrow(mem) == 0) {
      cli::cli_abort("No node of module(s) {.val {module_id}} in condition(s) {.val {conditions}}.")
    }
    drop <- drop | (vlayer %in% c("compound", "target") & !vname %in% mem$node_id)
  }
  changed <- any(drop)
  if (changed) g <- igraph::induced_subgraph(g, igraph::V(g)[!drop])

  if (min_target_degree > 1) {
    el <- igraph::as_data_frame(g, what = "edges")
    ct <- el[el$edge_kind == "compound_target", , drop = FALSE]
    n_compounds <- tapply(ct$from, ct$to, function(x) length(unique(x)))
    targets <- igraph::V(g)$name[igraph::V(g)$layer == "target"]
    n_hit <- n_compounds[targets]
    low <- targets[is.na(n_hit) | n_hit < min_target_degree]
    if (length(low) > 0) {
      g <- igraph::induced_subgraph(g, igraph::V(g)[!igraph::V(g)$name %in% low])
      changed <- TRUE
    }
  }
  ## Removing isolated vertices cannot isolate anything else -- one pass.
  if (changed) g <- igraph::induced_subgraph(g, igraph::V(g)[igraph::degree(g, mode = "all") > 0])
  g
}

#' File-name / log-key tag for `.network_layers_subset()`'s arguments
#' @return `""` when nothing is filtered, else e.g. `"mindeg2_mod-M1"`.
#' @keywords internal
.network_layers_subset_tag <- function(compound_ids = NULL, min_target_degree = 1, module_id = NULL) {
  parts <- character(0)
  if (!is.null(compound_ids)) parts <- c(parts, paste0("cmp-", substr(rlang::hash(sort(unique(compound_ids))), 1, 8)))
  if (!is.null(module_id)) parts <- c(parts, paste0("mod-", gsub("[^A-Za-z0-9_.-]+", "_", paste(sort(unique(module_id)), collapse = "+"))))
  if (min_target_degree > 1) parts <- c(parts, paste0("mindeg", min_target_degree))
  paste(parts, collapse = "_")
}

#' Plot-title suffix naming `.network_layers_subset()`'s active filters
#' @return `scope_label`, followed by e.g. `" (module M1; targets hit by >= 2 compounds)"`.
#' @keywords internal
.network_layers_subset_title <- function(scope_label, compound_ids = NULL, min_target_degree = 1, module_id = NULL) {
  parts <- character(0)
  if (!is.null(compound_ids)) parts <- c(parts, paste0(length(unique(compound_ids)), " selected compound(s)"))
  if (!is.null(module_id)) parts <- c(parts, paste0("module ", paste(sort(unique(module_id)), collapse = ", ")))
  if (min_target_degree > 1) parts <- c(parts, paste0("targets hit by >= ", min_target_degree, " compounds"))
  if (length(parts) == 0) scope_label else paste0(scope_label, " (", paste(parts, collapse = "; "), ")")
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
  if (!is.numeric(top_hub_n) || length(top_hub_n) != 1L || !is.finite(top_hub_n) ||
      top_hub_n < 0 || top_hub_n != floor(top_hub_n)) {
    cli::cli_abort("{.arg top_hub_n} must be a non-negative integer.")
  }
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
  hub_cut <- if (top_hub_n == 0) Inf else if (nrow(nodes) <= top_hub_n) -Inf else sort(nodes$degree, decreasing = TRUE)[top_hub_n]
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
#' module label includes its condition when more than one is pooled. The
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
  if (length(unique(conditions)) > 1) {
    mem$module_id <- ifelse(is.na(mem$module_id), NA_character_, paste0(mem$condition, ":", mem$module_id))
  }
  lookup <- stats::setNames(mem$module_id, mem$node_id)
  module_id <- unname(lookup[node_names])
  group <- ifelse(is.na(module_id), "not clustered", module_id)

  ## Colours are assigned over every module in scope (not only the ones
  ## drawn), so a module keeps its colour across views and subsets.
  palette <- .network_view_qual_palette(unique(c(stats::na.omit(mem$module_id), group[group != "not clustered"])))

  list(group = group, palette = palette, legend_name = "Module")
}

#' @keywords internal
.network_layers_layer_colors <- function() {
  c(compound = "#e67e22", target = "#2980b9", pathway = "#27ae60", disease = "#c0392b")
}

#' @keywords internal
.network_layers_ggplot <- function(pd, engine, title_suffix, fig_width = 9, edge_alpha = 1) {
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
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.2) * edge_alpha, guide = "none")

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

  ## Hub labels: long GC-MS names truncated for the drawn text only (the
  ## tooltip above keeps the full name), repelled off each other and off
  ## the hub points when ggrepel is installed -- hubs cluster in the dense
  ## core of a force-directed layout, where plain offset text piles up.
  hubs <- nodes[nodes$is_hub, , drop = FALSE]
  if (nrow(hubs) > 0) {
    hubs$short_label <- .plot_truncate(hubs$label)
    p <- p + .plot_text_layer(
      data = hubs,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$short_label),
      size = 3, fontface = "bold", colour = "grey15", show.legend = FALSE,
      repel_args = list(
        box.padding = 0.35, point.padding = 0.2, min.segment.length = 0.2,
        segment.colour = "grey40", max.overlaps = Inf,
        bg.colour = "white", bg.r = 0.12, seed = 1
      ),
      text_args = list(vjust = -1)
    )
  }

  layer_title <- paste(unique(nodes$layer[order(match(nodes$layer, c("compound", "target", "pathway", "disease")))]), collapse = "-")
  p +
    ggplot2::scale_fill_manual(values = pd$palette, name = pd$legend_name) +
    ## identity-scaled point sizes leave the legend keys at the tiny default
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::scale_size_identity(guide = "none") +
    ggplot2::coord_equal(xlim = pd$xlim, ylim = pd$ylim, clip = "off") +
    ggplot2::labs(
      title = .plot_wrap(paste0(layer_title, " network -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(
        "Node size = degree; labeled nodes are the top-degree hubs; dashed edges are derived (transitive) relationships",
        .plot_wrap_width(fig_width, 8)
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 12, 8, 12),
      legend.position = "right"
    )
}

## ---------------------------------------------------------------------------
## Simplified views (`view = "compound" / "bipartite" / "module" /
## "pathway"`). Real conditions have ~15-25 compounds and 400+ predicted
## targets with ~2000 edges; no force layout makes that readable. Each view
## below answers one question with a deterministic (or seeded) layout, and
## each is split into a pure data step (tested without pixels) and a
## drawing step.
## ---------------------------------------------------------------------------

#' Compound-target edge table of a (filtered) layered graph, with the
#' prediction probability attached
#'
#' @description
#' Only `edge_kind == "compound_target"` edges; `weight` is looked up in
#' `network_edges` (maximum across the conditions in scope when pooled;
#' `NA` when the table carries no weight).
#' @return `data.frame(compound_id, uniprot_id, weight)`, one row per pair.
#' @keywords internal
.network_layers_ct_table <- function(proj, g, conditions) {
  el <- igraph::as_data_frame(g, what = "edges")
  el <- unique(el[el$edge_kind == "compound_target", c("from", "to"), drop = FALSE])
  ct <- data.frame(compound_id = el$from, uniprot_id = el$to, stringsAsFactors = FALSE)
  ct$weight <- NA_real_
  edges_all <- patliRResults(proj, "network_edges")
  if (nrow(ct) > 0 && !is.null(edges_all) && "weight" %in% names(edges_all)) {
    e <- edges_all[edges_all$condition %in% conditions, , drop = FALSE]
    w <- tapply(e$weight, paste(e$compound_id, e$uniprot_id, sep = "\x01"),
                function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    ct$weight <- as.vector(w[paste(ct$compound_id, ct$uniprot_id, sep = "\x01")])
  }
  rownames(ct) <- NULL
  ct
}

#' Collapse a compound-target table to the targets shared by >= k compounds
#'
#' @description
#' The core of the simplified views. For every target, `n_compounds` is
#' the number of distinct compounds hitting it; targets with
#' `n_compounds >= k` are "drawn", the rest are only counted. With
#' `min_shared = NULL`, `k` is the smallest value >= 2 that leaves at most
#' `max_targets` drawn targets (the largest possible `k` if even that
#' leaves more; `1` when there is a single compound). Per compound, the
#' badge counts are `n_targets` (all its targets in scope), `n_unique`
#' (targets no other compound in scope hits), `n_drawn` and `n_hidden =
#' n_targets - n_drawn`.
#' @param ct `data.frame(compound_id, uniprot_id, weight)`, as
#'   `.network_layers_ct_table()` returns.
#' @param min_shared `NULL` (automatic `k`) or a single integer >= 1.
#' @param max_targets Target budget for the automatic `k`.
#' @return `list(k, n_compounds, targets = data.frame(uniprot_id,
#'   n_compounds, mean_weight, drawn), edges = the drawn rows of ct,
#'   badges = data.frame(compound_id, n_targets, n_unique, n_drawn,
#'   n_hidden))`. `targets` is sorted by `n_compounds`, then `mean_weight`
#'   (both decreasing), then id.
#' @keywords internal
.network_view_collapse <- function(ct, min_shared = NULL, max_targets = 40) {
  cmp <- sort(unique(ct$compound_id))
  n_cmp <- length(cmp)
  per_target <- tapply(ct$compound_id, ct$uniprot_id, function(x) length(unique(x)))
  targets <- data.frame(uniprot_id = names(per_target), n_compounds = as.integer(per_target),
                        stringsAsFactors = FALSE)
  mw <- tapply(ct$weight, ct$uniprot_id, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  targets$mean_weight <- as.vector(mw[targets$uniprot_id])

  if (is.null(min_shared)) {
    if (n_cmp <= 1 || nrow(targets) == 0) {
      k <- 1L
    } else {
      ks <- seq.int(2L, max(2L, max(targets$n_compounds)))
      n_ge <- vapply(ks, function(kk) sum(targets$n_compounds >= kk), integer(1))
      ok <- which(n_ge <= max_targets)
      k <- if (length(ok) > 0) ks[ok[1]] else max(ks)
    }
  } else {
    k <- as.integer(min_shared)
  }
  targets$drawn <- targets$n_compounds >= k
  mw_sort <- ifelse(is.na(targets$mean_weight), -Inf, targets$mean_weight)
  targets <- targets[order(-targets$n_compounds, -mw_sort, targets$uniprot_id), , drop = FALSE]
  rownames(targets) <- NULL

  edges <- ct[ct$uniprot_id %in% targets$uniprot_id[targets$drawn], , drop = FALSE]
  rownames(edges) <- NULL
  unique_targets <- targets$uniprot_id[targets$n_compounds == 1]
  count <- function(ids) as.integer(table(factor(ids, levels = cmp)))
  badges <- data.frame(
    compound_id = cmp,
    n_targets = count(ct$compound_id),
    n_unique = count(ct$compound_id[ct$uniprot_id %in% unique_targets]),
    n_drawn = count(edges$compound_id),
    stringsAsFactors = FALSE
  )
  badges$n_hidden <- badges$n_targets - badges$n_drawn
  list(k = k, n_compounds = n_cmp, targets = targets, edges = edges, badges = badges)
}

#' Order compounds so that compounds with similar target sets sit next to
#' each other
#'
#' @description
#' Average-linkage clustering of the Jaccard distance between compound
#' target sets (`stats::dist(method = "binary")`); the dendrogram leaf
#' order is used. Fewer than three compounds: decreasing target count.
#' Deterministic (no seed involved).
#' @return Character vector of compound ids.
#' @keywords internal
.network_view_compound_order <- function(ct) {
  cmp <- sort(unique(ct$compound_id))
  n_t <- as.integer(table(factor(ct$compound_id, levels = cmp)))
  if (length(cmp) < 3) return(cmp[order(-n_t, cmp)])
  m <- unclass(table(factor(ct$compound_id, levels = cmp), ct$uniprot_id)) > 0
  hc <- stats::hclust(stats::dist(m, method = "binary"), method = "average")
  cmp[hc$order]
}

#' Deterministic radial layout of the compound-centric view
#'
#' @description
#' Compounds sit on the unit circle in `compound_order`. Each drawn target
#' goes on a concentric shell whose radius encodes how many compounds hit
#' it (all of them: innermost shell; exactly `k`: outermost, radius
#' `r_max`). Within a shell, targets are ordered by the (probability
#' weighted) circular mean angle of the compounds hitting them and spaced
#' evenly, rotated to stay as close as possible to those angles -- a
#' target sits on the side of the compounds it links to, and no two
#' targets overlap.
#' @return `list(compounds = data.frame(name, x, y, angle), targets =
#'   col$targets[drawn] + x, y, shell_r)`.
#' @keywords internal
.network_view_compound_layout <- function(col, compound_order, r_min = 0.2, r_max = 0.78) {
  n <- length(compound_order)
  theta <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
  cmp <- data.frame(name = compound_order, x = cos(theta), y = sin(theta), angle = theta, stringsAsFactors = FALSE)
  tg <- col$targets[col$targets$drawn, , drop = FALSE]
  tg$x <- numeric(nrow(tg))
  tg$y <- numeric(nrow(tg))
  tg$shell_r <- numeric(nrow(tg))
  if (nrow(tg) == 0) return(list(compounds = cmp, targets = tg))

  ed <- col$edges
  w <- ifelse(is.na(ed$weight), 1, ed$weight)
  ci <- match(ed$compound_id, cmp$name)
  sx <- tapply(w * cos(cmp$angle[ci]), ed$uniprot_id, sum)
  sy <- tapply(w * sin(cmp$angle[ci]), ed$uniprot_id, sum)
  bary <- atan2(sy[tg$uniprot_id], sx[tg$uniprot_id])

  levels_n <- sort(unique(tg$n_compounds), decreasing = TRUE)
  n_top <- col$n_compounds
  n_low <- min(col$k, min(tg$n_compounds))
  shell_radius <- function(nc) {
    if (n_top == n_low) return((r_min + r_max) / 2)
    r_min + (r_max - r_min) * (n_top - nc) / (n_top - n_low)
  }
  for (lv in levels_n) {
    idx <- which(tg$n_compounds == lv)
    r <- shell_radius(lv)
    m <- length(idx)
    if (m == 1) {
      a <- bary[idx]
      if (r < 1e-9) a <- 0
    } else {
      o <- idx[order(bary[idx])]
      slots <- 2 * pi * (seq_len(m) - 1) / m
      offset <- atan2(sum(sin(bary[o] - slots)), sum(cos(bary[o] - slots)))
      a <- numeric(m)
      a[match(o, idx)] <- offset + slots
    }
    tg$x[idx] <- r * cos(a)
    tg$y[idx] <- r * sin(a)
    tg$shell_r[idx] <- r
  }
  list(compounds = cmp, targets = tg)
}

#' Two-column layout of the bipartite view
#'
#' @description
#' Targets on the right (`x = 1`), grouped by module (largest module on
#' top, `"not clustered"` last, a gap of `gap` slots between groups) and,
#' within a group, by decreasing `n_compounds`, then `mean_weight`.
#' Compounds on the left (`x = 0`), evenly spaced and ordered by the mean
#' height of their targets (barycentric ordering -- the classic
#' crossing-reduction heuristic for two-layer drawings; Sugiyama et al.
#' 1981), so each compound faces the block of targets it mostly hits.
#' When `compound_module` is given, compounds are first grouped in the
#' same module order as the targets (barycentric order within a group),
#' so a module's compounds face its target block. Compounds without any
#' drawn target go last within their group.
#' @param edges The drawn compound-target rows.
#' @param targets The drawn targets (`uniprot_id`, `n_compounds`,
#'   `mean_weight`).
#' @param compounds Every compound id to place.
#' @param module `NULL` or a named character vector (`uniprot_id` ->
#'   module); `NA` counts as `"not clustered"`.
#' @param compound_module `NULL` or a named character vector (compound id
#'   -> module), same convention.
#' @return `list(compounds = data.frame(name, bary, module, x, y), targets
#'   = targets + module, x, y)`, both sorted top to bottom.
#' @keywords internal
.network_view_bipartite_layout <- function(edges, targets, compounds, module = NULL, compound_module = NULL, gap = 3) {
  targets$module <- if (is.null(module)) "all" else unname(module[targets$uniprot_id])
  targets$module[is.na(targets$module)] <- "not clustered"
  mod_order <- .network_view_module_order(targets$module)
  mw <- ifelse(is.na(targets$mean_weight), -Inf, targets$mean_weight)
  o <- order(match(targets$module, mod_order), -targets$n_compounds, -mw, targets$uniprot_id)
  targets <- targets[o, , drop = FALSE]
  rownames(targets) <- NULL
  slot <- seq_len(nrow(targets)) + gap * (match(targets$module, mod_order) - 1)
  span <- max(1, max(c(slot, 1)) - 1)
  targets$x <- rep(1, nrow(targets))
  targets$y <- 1 - (slot - 1) / span

  ty <- stats::setNames(targets$y, targets$uniprot_id)
  bary <- tapply(ty[edges$uniprot_id], edges$compound_id, mean)
  cmp <- data.frame(name = compounds, bary = unname(bary[compounds]), stringsAsFactors = FALSE)
  cmp$module <- if (is.null(compound_module)) "all" else unname(compound_module[compounds])
  cmp$module[is.na(cmp$module)] <- "not clustered"
  cmod_order <- unique(c(mod_order, .network_view_module_order(cmp$module)))
  cmod_order <- c(setdiff(cmod_order, "not clustered"), intersect(cmod_order, "not clustered"))
  cmp <- cmp[order(match(cmp$module, cmod_order), -ifelse(is.na(cmp$bary), -Inf, cmp$bary), cmp$name), , drop = FALSE]
  rownames(cmp) <- NULL
  n <- nrow(cmp)
  cmp$x <- rep(0, n)
  cmp$y <- if (n == 1) 0.5 else 0.97 - 0.94 * (seq_len(n) - 1) / (n - 1)
  list(compounds = cmp, targets = targets)
}

#' Module display order: largest first, "not clustered" last
#' @keywords internal
.network_view_module_order <- function(module) {
  tab <- table(module)
  tab <- tab[order(-as.integer(tab), names(tab))]
  mods <- names(tab)
  c(setdiff(mods, "not clustered"), intersect(mods, "not clustered"))
}

#' Per-module layout for the module view: each module laid out on its own,
#' then the blocks packed row by row
#'
#' @description
#' Every module's nodes get a Fruchterman-Reingold layout of the edges
#' inside that module only, normalized to a disc whose area is
#' proportional to the module's node count. When `nodes` has a `layer`
#' column and a module holds two or more compounds, those compounds are
#' pinned on an inner ring (similar compounds adjacent) and only the
#' targets move -- otherwise compounds sharing most of their targets pile
#' up in the middle. Discs are packed in rows, largest module first.
#' Edges between modules therefore never pull layouts together -- the
#' grouping stays visible.
#' @param nodes `data.frame(name, module)`, optionally with `layer`.
#' @param edges `data.frame(from, to)`.
#' @return `list(nodes = nodes + x, y; centres = data.frame(module, x, y,
#'   radius))`.
#' @keywords internal
.network_view_module_layout <- function(nodes, edges, seed = 1) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  mods <- .network_view_module_order(nodes$module)
  nodes$x <- 0
  nodes$y <- 0
  blocks <- lapply(mods, function(m) {
    nm <- nodes$name[nodes$module == m]
    e <- edges[edges$from %in% nm & edges$to %in% nm, c("from", "to"), drop = FALSE]
    is_c <- if ("layer" %in% names(nodes)) nodes$layer[match(nm, nodes$name)] == "compound" else rep(FALSE, length(nm))
    nc <- sum(is_c)
    xy <- if (length(nm) == 1) {
      matrix(0, 1, 2)
    } else if (nc >= 2 && length(nm) > nc) {
      ## Several compounds share most of their targets, so a free force
      ## layout piles them on top of each other in the middle. Pin them on
      ## an inner ring (similar compounds adjacent) and let only the
      ## targets move.
      gm <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = nm, stringsAsFactors = FALSE))
      s <- sqrt(length(nm))
      c_names <- nm[is_c]
      c_order <- if (nrow(e) > 0) .network_view_compound_order(data.frame(compound_id = e$from, uniprot_id = e$to, stringsAsFactors = FALSE)) else character(0)
      c_order <- c(intersect(c_order, c_names), setdiff(c_names, c_order))
      ang <- pi / 2 - 2 * pi * (seq_len(nc) - 1) / nc
      ci <- match(c_order, nm)
      init <- matrix(stats::runif(2 * length(nm), -0.5 * s, 0.5 * s), ncol = 2)
      init[ci, ] <- cbind(0.45 * s * cos(ang), 0.45 * s * sin(ang))
      lo_x <- rep(-s, length(nm))
      hi_x <- rep(s, length(nm))
      lo_y <- lo_x
      hi_y <- hi_x
      lo_x[ci] <- hi_x[ci] <- init[ci, 1]
      lo_y[ci] <- hi_y[ci] <- init[ci, 2]
      igraph::layout_with_fr(gm, coords = init, niter = 1500, minx = lo_x, maxx = hi_x, miny = lo_y, maxy = hi_y)
    } else {
      gm <- igraph::graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = nm, stringsAsFactors = FALSE))
      igraph::layout_with_fr(gm, niter = 1500)
    }
    xy <- sweep(xy, 2, colMeans(xy))
    r <- max(sqrt(rowSums(xy^2)))
    if (r > 0) xy <- xy / r
    radius <- sqrt(length(nm))
    list(name = nm, xy = xy * radius, radius = radius)
  })
  radii <- vapply(blocks, `[[`, numeric(1), "radius")
  ## Small modules still get room for their two-line title.
  slot_r <- pmax(radii, 0.45 * max(radii))
  ncol <- ceiling(sqrt(length(mods)))
  pad <- 0.35 * max(radii) + 0.6
  centres <- data.frame(module = mods, x = 0, y = 0, radius = radii, stringsAsFactors = FALSE)
  y_cursor <- 0
  for (row_start in seq(1, length(mods), by = ncol)) {
    idx <- row_start:min(length(mods), row_start + ncol - 1)
    row_h <- max(slot_r[idx])
    x_cursor <- 0
    for (i in idx) {
      centres$x[i] <- x_cursor + slot_r[i]
      centres$y[i] <- y_cursor - row_h
      x_cursor <- x_cursor + 2 * slot_r[i] + pad
    }
    y_cursor <- y_cursor - 2 * row_h - pad
  }
  for (i in seq_along(blocks)) {
    at <- match(blocks[[i]]$name, nodes$name)
    nodes$x[at] <- blocks[[i]]$xy[, 1] + centres$x[i]
    nodes$y[at] <- blocks[[i]]$xy[, 2] + centres$y[i]
  }
  list(nodes = nodes, centres = centres)
}

#' Padded convex hull per group, for the module view's shaded regions
#'
#' @description
#' Each point is replaced by a small circle of radius `pad` (16 points)
#' before [grDevices::chull()], so single nodes and collinear groups still
#' get a visible rounded region.
#' @return `data.frame(group, x, y)` polygon vertices, in hull order.
#' @keywords internal
.network_view_hulls <- function(x, y, group, pad = 0.3) {
  ang <- seq(0, 2 * pi, length.out = 17)[-17]
  out <- lapply(unique(group), function(gr) {
    i <- which(group == gr)
    px <- as.vector(outer(x[i], pad * cos(ang), `+`))
    py <- as.vector(outer(y[i], pad * sin(ang), `+`))
    h <- grDevices::chull(px, py)
    data.frame(group = gr, x = px[h], y = py[h], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Soft module lookup: `NULL` when [network_module_robustness()] has not
#' been run, otherwise node id -> module (condition-prefixed when pooled)
#' @keywords internal
.network_view_module_lookup <- function(proj, conditions, node_names) {
  membership_all <- patliRResults(proj, "network_module_membership")
  if (is.null(membership_all) || nrow(membership_all) == 0 ||
      !any(membership_all$condition %in% conditions)) {
    return(NULL)
  }
  info <- .network_layers_module_colour_info(proj, conditions, node_names)
  out <- stats::setNames(info$group, node_names)
  out[out == "not clustered"] <- NA_character_
  attr(out, "palette") <- info$palette
  out
}

#' Colour-blind-safe qualitative palette (Okabe & Ito 2008), falling back
#' to `hcl.colors(n, "Dark 3")` beyond eight levels
#'
#' @description
#' Levels are sorted naturally (`"M2"` before `"M10"`) before colours are
#' assigned, so the same module gets the same colour whatever order the
#' caller lists them in; `"not clustered"` is always grey.
#' @keywords internal
.network_view_qual_palette <- function(levels) {
  levels <- unique(setdiff(levels, "not clustered"))
  levels <- levels[order(nchar(levels), levels)]
  oi <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#999933", "#F0E442")
  cols <- if (length(levels) <= length(oi)) oi[seq_along(levels)] else grDevices::hcl.colors(length(levels), "Dark 3")
  c(stats::setNames(cols, levels), `not clustered` = "grey70")
}

#' Colours shared by the simplified views
#' @keywords internal
.network_view_colours <- function() {
  c(compound = "#D55E00", target = "#0072B2", edge = "grey45")
}

#' Node layer that becomes a ggiraph interactive layer (tooltip = `tooltip`
#' column, `data_id` = `name`) when `engine = "ggiraph"`
#' @param aes_cols Named list: aesthetic -> column name (besides `x`/`y`).
#' @keywords internal
.network_view_points <- function(engine, data, aes_cols = list(), ...) {
  interactive <- identical(engine, "ggiraph") && requireNamespace("ggiraph", quietly = TRUE)
  cols <- c(list(x = "x", y = "y"), aes_cols)
  if (interactive) cols <- c(cols, list(tooltip = "tooltip", data_id = "name"))
  mapping <- ggplot2::aes(!!!lapply(cols, function(col) rlang::expr(.data[[!!col]])))
  if (interactive) {
    ggiraph::geom_point_interactive(data = data, mapping = mapping, ...)
  } else {
    ggplot2::geom_point(data = data, mapping = mapping, ...)
  }
}

#' Common theme of the simplified views (white background, no axes)
#' @keywords internal
.network_view_theme <- function(legend_position = "right") {
  ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30", lineheight = 1.1),
      plot.caption = ggplot2::element_text(size = 7, colour = "grey40", hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(8, 12, 8, 12),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.size = ggplot2::unit(0.45, "cm")
    )
}

#' Compound-centric summary view (`view = "compound"`)
#'
#' @description
#' Compounds on an outer ring (ordered by target-set similarity), the
#' targets shared by >= `k` compounds on concentric inner shells (more
#' shared = closer to the centre), everything else collapsed into a
#' per-compound badge. See `.network_view_collapse()` and
#' `.network_view_compound_layout()`.
#' @return `list(plot, view_tag, n_nodes, n_edges, data = list(collapse,
#'   compounds, targets))`.
#' @keywords internal
.network_view_compound <- function(proj, ct, conditions, min_shared = NULL, max_targets = 40, top_hub_n = 15,
                                   colour_by = "layer", seed = 1, engine = "static", title_suffix = "",
                                   note = NULL, fig_width = 9) {
  col <- .network_view_collapse(ct, min_shared = min_shared, max_targets = max_targets)
  lay <- .network_view_compound_layout(col, .network_view_compound_order(ct))
  cols <- .network_view_colours()

  cmp <- merge(lay$compounds, col$badges, by.x = "name", by.y = "compound_id", sort = FALSE)
  cmp <- cmp[match(lay$compounds$name, cmp$name), , drop = FALSE]
  cmp$label <- .plot_unique_labels(proj, conditions, cmp$name, "compound")
  cmp$short_label <- .plot_truncate(cmp$label, 26)
  cmp$badge <- sprintf("%d targets, +%d unique", cmp$n_targets, cmp$n_unique)
  cmp$tooltip <- sprintf("%s\n%d predicted targets\n%d hit by no other compound\n%d drawn (shared by >= %d)",
                         cmp$label, cmp$n_targets, cmp$n_unique, cmp$n_drawn, col$k)
  ## Labels (name + badge) sit just outside the ring, left-aligned on the
  ## right half and right-aligned on the left half (centred only for an
  ## exactly top/bottom compound). Near the top the two lines stack upwards
  ## from the anchor, near the bottom downwards, so they never cover a node.
  ca <- cos(cmp$angle)
  sa <- sin(cmp$angle)
  side <- ifelse(abs(ca) < 0.02, 0, sign(ca))
  cmp$lx <- 1.1 * ca
  cmp$ly <- 1.1 * sa
  cmp$hjust <- ifelse(side > 0, 0, ifelse(side < 0, 1, 0.5))
  band <- ifelse(sa > 0.45, 1, ifelse(sa < -0.45, -1, 0))
  cmp$name_vjust <- c(1.2, -0.15, -1.3)[band + 2]
  cmp$badge_vjust <- c(2.7, 1.25, -0.15)[band + 2]
  cmp$point_size <- 5 + 5 * sqrt(cmp$n_targets / max(cmp$n_targets))

  tg <- lay$targets
  tg$name <- tg$uniprot_id
  tg$label <- .plot_label_nodes(proj, conditions, tg$uniprot_id, rep("target", nrow(tg)))
  by_module <- colour_by == "module" && nrow(tg) > 0
  if (by_module) {
    info <- .network_layers_module_colour_info(proj, conditions, tg$uniprot_id)
    tg$fill_group <- info$group
    palette <- info$palette
  }
  tg$tooltip <- sprintf("%s (%s)\nhit by %d of %d compounds\nmean probability %.2f",
                        tg$label, tg$uniprot_id, tg$n_compounds, col$n_compounds, tg$mean_weight)

  ed <- col$edges
  ed$x <- cmp$x[match(ed$compound_id, cmp$name)]
  ed$y <- cmp$y[match(ed$compound_id, cmp$name)]
  ed$xend <- tg$x[match(ed$uniprot_id, tg$uniprot_id)]
  ed$yend <- tg$y[match(ed$uniprot_id, tg$uniprot_id)]
  ## Probabilities are clamped to [0.5, 1] (the usual import cut-off is 0.5)
  ## so the width legend is identical across figures.
  ed$w <- pmin(pmax(ifelse(is.na(ed$weight), 0.5, ed$weight), 0.5), 1)

  shells <- unique(tg[, c("n_compounds", "shell_r"), drop = FALSE])
  circle <- function(r, n = 180) {
    a <- seq(0, 2 * pi, length.out = n)
    data.frame(x = r * cos(a), y = r * sin(a))
  }
  guides <- do.call(rbind, c(
    list(cbind(circle(1), ring = "compounds")),
    lapply(seq_len(nrow(shells)), function(i) cbind(circle(shells$shell_r[i]), ring = paste0("s", i)))
  ))

  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = guides, ggplot2::aes(x = .data$x, y = .data$y, group = .data$ring),
                       colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_segment(
      data = ed,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linewidth = .data$w),
      colour = cols[["edge"]], alpha = min(0.35, max(0.08, 40 / max(1, nrow(ed))))
    ) +
    ggplot2::scale_linewidth(range = c(0.1, 0.9), limits = c(0.5, 1),
                             name = "Prediction\nprobability", breaks = c(0.5, 0.75, 1))
  if (nrow(tg) > 0) {
    p <- p + .network_view_points(engine, tg, list(size = "n_compounds", fill = if (by_module) "fill_group" else "n_compounds"),
                                  shape = 21, colour = "grey25", stroke = 0.3)
  }
  p <- p + .network_view_points(engine, cmp, list(), size = cmp$point_size, shape = 21,
                                fill = cols[["compound"]], colour = "white", stroke = 0.6)

  ## Target labels: the most shared drawn targets.
  lab_t <- utils::head(tg[order(-tg$n_compounds, -ifelse(is.na(tg$mean_weight), -Inf, tg$mean_weight)), , drop = FALSE], top_hub_n)
  if (nrow(lab_t) > 0) {
    p <- p + .plot_text_layer(
      data = lab_t, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.5, colour = "grey10", fontface = "italic", show.legend = FALSE,
      repel_args = list(box.padding = 0.25, point.padding = 0.15, min.segment.length = 0.1,
                        segment.colour = "grey50", segment.size = 0.2, max.overlaps = Inf,
                        bg.colour = "white", bg.r = 0.15, seed = 1),
      text_args = list(vjust = -0.9)
    )
  }
  p <- p +
    ggplot2::geom_text(data = cmp, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$short_label,
                                                hjust = .data$hjust, vjust = .data$name_vjust),
                       size = 2.9, fontface = "bold", colour = "grey10") +
    ggplot2::geom_text(data = cmp, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$badge,
                                                hjust = .data$hjust, vjust = .data$badge_vjust),
                       size = 2.4, colour = "grey35")

  ## One merged legend for "compounds sharing the target" (size + fill),
  ## unless the fill is taken by the module colours.
  n_lo <- if (nrow(tg) > 0) min(col$k, min(tg$n_compounds)) else col$k
  n_hi <- max(col$n_compounds, n_lo + 1)
  n_lo <- min(n_lo, n_hi - 1)
  size_breaks <- unique(round(seq(max(1, n_lo), n_hi, length.out = min(4, n_hi - n_lo + 1))))
  size_name <- "Compounds\nsharing target"
  p <- p + ggplot2::scale_size(range = c(2.6, 6), limits = c(n_lo, n_hi), breaks = size_breaks, name = size_name)
  if (by_module) {
    p <- p +
      ggplot2::scale_fill_manual(values = palette, name = "Target module") +
      ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "grey60"), order = 1),
                      fill = ggplot2::guide_legend(override.aes = list(size = 4), order = 2),
                      linewidth = ggplot2::guide_legend(override.aes = list(alpha = 0.7), order = 3))
  } else {
    p <- p +
      ggplot2::scale_fill_viridis_c(option = "D", direction = -1, end = 0.92, limits = c(n_lo, n_hi),
                                    breaks = size_breaks, name = size_name, guide = "legend") +
      ggplot2::guides(size = ggplot2::guide_legend(order = 1), fill = ggplot2::guide_legend(order = 1),
                      linewidth = ggplot2::guide_legend(override.aes = list(alpha = 0.7), order = 2))
  }
  p <- p + ggplot2::coord_equal(xlim = c(-1.78, 1.78), ylim = c(-1.26, 1.24), clip = "off")

  n_hidden_t <- sum(!col$targets$drawn)
  sub <- sprintf(paste0(
    "Ring: %d compounds, ordered by target-set similarity (Jaccard). Inside: the %d targets hit by >= %d compounds; ",
    "the more compounds share a target, the closer to the centre, the larger and the darker it is. The other %d targets ",
    "are not drawn: each compound's badge gives its total and its unique (hit by no other compound) targets. ",
    "Edge width = prediction probability."
  ), col$n_compounds, sum(col$targets$drawn), col$k, n_hidden_t)
  if (!is.null(note)) sub <- paste(note, sub)
  p <- p +
    ggplot2::labs(
      title = .plot_wrap(paste0("Compound-target summary -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(sub, .plot_wrap_width(fig_width, 8)),
      x = NULL, y = NULL
    ) +
    .network_view_theme(legend_position = "bottom")

  list(
    plot = p,
    view_tag = paste0("k", col$k),
    n_nodes = nrow(cmp) + nrow(tg), n_edges = nrow(ed),
    data = list(collapse = col, compounds = cmp, targets = tg)
  )
}

#' Bipartite two-column view (`view = "bipartite"`)
#'
#' @description
#' Compounds left, targets right, targets grouped by module when
#' [network_module_robustness()] has been run (coloured band + edge colour
#' = module), otherwise by how many compounds hit them. Edges are drawn
#' nearly transparent so that bundles, not single lines, carry the
#' message; only the `top_hub_n` most shared targets are labelled. With
#' `min_shared = k > 1` the targets hit by fewer compounds are dropped
#' (counted in the compound labels).
#' @return `list(plot, view_tag, n_nodes, n_edges, data = list(collapse,
#'   compounds, targets))`.
#' @keywords internal
.network_view_bipartite <- function(proj, ct, conditions, min_shared = NULL, top_hub_n = 15, engine = "static",
                                    title_suffix = "", fig_width = 9) {
  k <- if (is.null(min_shared)) 1L else as.integer(min_shared)
  col <- .network_view_collapse(ct, min_shared = k)
  tg0 <- col$targets[col$targets$drawn, , drop = FALSE]
  module <- .network_view_module_lookup(proj, conditions, tg0$uniprot_id)
  cmodule <- if (is.null(module)) NULL else .network_view_module_lookup(proj, conditions, col$badges$compound_id)
  lay <- .network_view_bipartite_layout(col$edges, tg0, col$badges$compound_id, module = module, compound_module = cmodule)
  cols <- .network_view_colours()

  cmp <- lay$compounds
  b <- col$badges[match(cmp$name, col$badges$compound_id), , drop = FALSE]
  cmp$n_targets <- b$n_targets
  cmp$n_unique <- b$n_unique
  cmp$label <- .plot_unique_labels(proj, conditions, cmp$name, "compound")
  cmp$short_label <- paste0(.plot_truncate(cmp$label, 30), "  (", cmp$n_targets, ")")
  cmp$tooltip <- sprintf("%s\n%d predicted targets, %d unique", cmp$label, cmp$n_targets, cmp$n_unique)

  tg <- lay$targets
  tg$name <- tg$uniprot_id
  tg$label <- .plot_label_nodes(proj, conditions, tg$uniprot_id, rep("target", nrow(tg)))
  tg$tooltip <- sprintf("%s (%s)\nhit by %d compounds\nmodule %s", tg$label, tg$uniprot_id, tg$n_compounds, tg$module)
  has_module <- !is.null(module)
  mod_levels <- .network_view_module_order(tg$module)
  if (has_module) {
    palette <- attr(module, "palette")
  } else {
    cmp$module <- "compound"
    palette <- c(all = cols[["target"]], compound = cols[["compound"]])
  }

  ed <- col$edges
  ed$y <- cmp$y[match(ed$compound_id, cmp$name)]
  ed$yend <- tg$y[match(ed$uniprot_id, tg$uniprot_id)]
  ed$module <- tg$module[match(ed$uniprot_id, tg$uniprot_id)]
  ed_alpha <- min(0.4, max(0.05, 110 / max(1, nrow(ed))))

  bands <- do.call(rbind, lapply(mod_levels, function(m) {
    yy <- tg$y[tg$module == m]
    data.frame(module = m, ymin = min(yy) - 0.004, ymax = max(yy) + 0.004, n = length(yy), stringsAsFactors = FALSE)
  }))
  bands$ymid <- (bands$ymin + bands$ymax) / 2
  bands$text <- paste0(bands$module, "\n", bands$n, " targets")

  ## Labels: the most shared targets of each module block (so every
  ## block gets named hubs), `top_hub_n` in total.
  per_block <- max(1L, ceiling(top_hub_n / length(mod_levels)))
  tg_rank <- tg[order(-tg$n_compounds, -ifelse(is.na(tg$mean_weight), -Inf, tg$mean_weight)), , drop = FALSE]
  lab_t <- do.call(rbind, lapply(mod_levels, function(m) utils::head(tg_rank[tg_rank$module == m, , drop = FALSE], per_block)))
  lab_t <- utils::head(lab_t[order(-lab_t$n_compounds), , drop = FALSE], top_hub_n)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = ed,
      ggplot2::aes(x = 0, y = .data$y, xend = 1, yend = .data$yend, colour = .data$module),
      linewidth = 0.25, alpha = ed_alpha
    ) +
    ggplot2::scale_colour_manual(values = palette, guide = "none")
  if (has_module) {
    p <- p +
      ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = 1.025, xmax = 1.05, ymin = .data$ymin, ymax = .data$ymax, fill = .data$module)) +
      ggplot2::geom_text(data = bands, ggplot2::aes(x = 1.065, y = .data$ymid, label = .data$text),
                         hjust = 0, size = 2.4, lineheight = 0.9, colour = "grey25")
  }
  p <- p +
    .network_view_points(engine, tg, list(size = "n_compounds", fill = "module"), shape = 21, colour = "grey20", stroke = 0.1) +
    .network_view_points(engine, cmp, list(fill = "module"), shape = 23, size = 2.5 + 3.5 * sqrt(cmp$n_targets / max(cmp$n_targets)),
                         colour = "grey10", stroke = 0.5) +
    ggplot2::geom_text(data = cmp, ggplot2::aes(x = -0.04, y = .data$y, label = .data$short_label),
                       hjust = 1, size = 2.7, colour = "grey10")
  if (nrow(lab_t) > 0) {
    lab_x <- if (has_module) 1.28 else 1.06
    p <- p + .plot_text_layer(
      data = lab_t, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 2.4, fontface = "italic", colour = "grey10", show.legend = FALSE, hjust = 0,
      repel_args = list(nudge_x = lab_x - 1, direction = "y", xlim = c(lab_x, NA), min.segment.length = 0,
                        segment.colour = "grey60", segment.size = 0.2, box.padding = 0.1,
                        max.overlaps = Inf, seed = 1),
      text_args = list(nudge_x = 0.02)
    )
  }
  p <- p +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_size_area(max_size = 2.4, name = "Compounds\nhitting target",
                             breaks = unique(round(seq(max(1, k), max(tg$n_compounds, k), length.out = 4)))) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "grey40", colour = "grey20"))) +
    ggplot2::scale_x_continuous(limits = c(-0.95, 1.55), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off")

  sub <- sprintf(paste0(
    "Left: %d compounds (diamonds; total predicted targets in brackets%s), ordered so each faces the targets it mostly hits. ",
    "Right: %d targets%s, %s. Edges are drawn transparent: dense bundles, not single lines, show where compounds converge. ",
    "Labelled: the %d targets hit by most compounds%s."
  ), nrow(cmp), if (has_module) "; colour and grouping = the compound's own module" else "", nrow(tg),
  if (k > 1) sprintf(" hit by >= %d compounds (%d others hidden)", k, sum(!col$targets$drawn)) else "",
  if (has_module) "grouped by network module (coloured band) and, within a module, by how many compounds hit them" else "sorted by how many compounds hit them",
  nrow(lab_t), if (has_module) " within each module" else "")
  p <- p +
    ggplot2::labs(
      title = .plot_wrap(paste0("Compound-target bipartite view -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(sub, .plot_wrap_width(fig_width, 8)),
      x = NULL, y = NULL
    ) +
    .network_view_theme()

  list(
    plot = p,
    view_tag = if (k > 1) paste0("k", k) else "",
    n_nodes = nrow(cmp) + nrow(tg), n_edges = nrow(ed),
    data = list(collapse = col, compounds = cmp, targets = tg)
  )
}

#' Module view (`view = "module"`)
#'
#' @description
#' One shaded region per [network_module_robustness()] module, each laid
#' out on its own (`.network_view_module_layout()`); edges inside a module
#' drawn in the module colour, edges between modules faint grey. Every
#' compound is labelled; per module, the targets with the highest
#' `network_centrality` hub score (degree when that table is missing) are
#' labelled, `top_hub_n` in total.
#' @return `list(plot, view_tag, n_nodes, n_edges, data = list(collapse,
#'   nodes, centres))`.
#' @keywords internal
.network_view_module <- function(proj, ct, conditions, min_shared = NULL, top_hub_n = 15, seed = 1,
                                 engine = "static", title_suffix = "", fig_width = 9) {
  k <- if (is.null(min_shared)) 1L else as.integer(min_shared)
  col <- .network_view_collapse(ct, min_shared = k)
  ed <- col$edges
  nodes <- data.frame(
    name = c(col$badges$compound_id, col$targets$uniprot_id[col$targets$drawn]),
    layer = c(rep("compound", nrow(col$badges)), rep("target", sum(col$targets$drawn))),
    stringsAsFactors = FALSE
  )
  info <- .network_layers_module_colour_info(proj, conditions, nodes$name) # aborts without membership
  nodes$module <- info$group
  lay <- .network_view_module_layout(nodes, data.frame(from = ed$compound_id, to = ed$uniprot_id, stringsAsFactors = FALSE), seed = seed)
  nodes <- lay$nodes
  mods <- .network_view_module_order(nodes$module)
  palette <- info$palette
  cols <- .network_view_colours()

  nodes$degree <- as.integer(table(factor(c(ed$compound_id, ed$uniprot_id), levels = nodes$name)))
  cen <- patliRResults(proj, "network_centrality")
  hub <- stats::setNames(rep(NA_real_, nrow(nodes)), nodes$name)
  if (!is.null(cen) && "hub_score" %in% names(cen)) {
    cc <- cen[cen$condition %in% conditions, , drop = FALSE]
    hs <- tapply(cc$hub_score, cc$node_id, max, na.rm = TRUE)
    hub[] <- unname(hs[nodes$name])
  }
  nodes$score <- ifelse(is.na(hub), nodes$degree / max(1, max(nodes$degree)), hub)
  nodes$label <- .plot_unique_labels(proj, conditions, nodes$name, nodes$layer)
  nodes$label[nodes$layer == "target"] <- .plot_label_nodes(proj, conditions, nodes$name[nodes$layer == "target"],
                                                             rep("target", sum(nodes$layer == "target")))
  nodes$tooltip <- sprintf("%s\n%s, module %s\n%d links", nodes$label, nodes$layer, nodes$module, nodes$degree)

  per_mod <- max(1L, ceiling(top_hub_n / max(1, length(mods))))
  tnodes <- nodes[nodes$layer == "target", , drop = FALSE]
  lab_t <- do.call(rbind, lapply(split(tnodes, tnodes$module), function(d) utils::head(d[order(-d$score, d$name), , drop = FALSE], per_mod)))
  if (top_hub_n == 0 || is.null(lab_t)) lab_t <- tnodes[0, , drop = FALSE]
  lab_c <- nodes[nodes$layer == "compound", , drop = FALSE]
  lab_c$short_label <- .plot_truncate(lab_c$label, 24)

  mi <- match(ed$compound_id, nodes$name)
  ti <- match(ed$uniprot_id, nodes$name)
  segs <- data.frame(x = nodes$x[mi], y = nodes$y[mi], xend = nodes$x[ti], yend = nodes$y[ti],
                     module = nodes$module[mi], within = nodes$module[mi] == nodes$module[ti], stringsAsFactors = FALSE)

  hulls <- .network_view_hulls(nodes$x, nodes$y, nodes$module, pad = 0.35)
  cen_df <- lay$centres
  cen_df$n_c <- vapply(cen_df$module, function(m) sum(nodes$module == m & nodes$layer == "compound"), integer(1))
  cen_df$n_t <- vapply(cen_df$module, function(m) sum(nodes$module == m & nodes$layer == "target"), integer(1))
  top_y <- tapply(hulls$y, hulls$group, max)
  cen_df$ly <- unname(top_y[cen_df$module]) + 0.1 + 0.08 * max(cen_df$radius)
  cen_df$text <- sprintf(if (nrow(cen_df) > 5) "%s\n%d compound%s, %d targets" else "%s: %d compound%s, %d targets",
                         cen_df$module, cen_df$n_c, ifelse(cen_df$n_c == 1, "", "s"), cen_df$n_t)

  nodes$size_value <- ifelse(nodes$layer == "compound", NA_real_, nodes$degree)
  tn <- nodes[nodes$layer == "target", , drop = FALSE]
  cn <- nodes[nodes$layer == "compound", , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$group),
                          alpha = 0.10, colour = NA) +
    ggplot2::geom_polygon(data = hulls, ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, colour = .data$group),
                          fill = NA, linewidth = 0.3, linetype = "22") +
    ggplot2::geom_segment(data = segs[!segs$within, , drop = FALSE],
                          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
                          colour = "grey70", linewidth = 0.12, alpha = 0.10) +
    ggplot2::geom_segment(data = segs[segs$within, , drop = FALSE],
                          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, colour = .data$module),
                          linewidth = 0.2, alpha = 0.35) +
    .network_view_points(engine, tn, list(fill = "module", size = "size_value"), shape = 21, colour = "white", stroke = 0.2) +
    .network_view_points(engine, cn, list(fill = "module"), shape = 23, size = 4.2, colour = "grey10", stroke = 0.5) +
    ggplot2::geom_text(data = cen_df, ggplot2::aes(x = .data$x, y = .data$ly, label = .data$text, colour = .data$module),
                       vjust = 0, size = 2.9, lineheight = 0.9, fontface = "bold", show.legend = FALSE)
  lab_all <- rbind(
    data.frame(x = lab_c$x, y = lab_c$y, label = lab_c$short_label, face = "bold", stringsAsFactors = FALSE),
    if (nrow(lab_t) > 0) data.frame(x = lab_t$x, y = lab_t$y, label = lab_t$label, face = "italic", stringsAsFactors = FALSE)
  )
  p <- p + .plot_text_layer(
    data = lab_all, mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label, fontface = .data$face),
    size = 2.5, colour = "grey10", show.legend = FALSE,
    repel_args = list(box.padding = 0.3, point.padding = 0.15, min.segment.length = 0.1, segment.colour = "grey40",
                      segment.size = 0.2, max.overlaps = Inf, bg.colour = "white", bg.r = 0.15, seed = 1),
    text_args = list(vjust = -0.9)
  )
  p <- p +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::scale_colour_manual(values = palette, guide = "none") +
    ggplot2::scale_size_area(max_size = 3.5, name = "Compounds\nhitting target") +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "grey50", colour = "white"))) +
    ggplot2::coord_equal(clip = "off")
  sub <- sprintf(paste0(
    "Nodes grouped by network module (Leiden, network_module_robustness()); diamonds = compounds (all labelled), circles = targets ",
    "(size = number of compounds hitting them). Edges within a module in its colour, between modules faint grey. ",
    "Labelled targets: the top %d per module by hub score%s."
  ), per_mod, if (k > 1) sprintf("; only targets hit by >= %d compounds are drawn", k) else "")
  p <- p +
    ggplot2::labs(
      title = .plot_wrap(paste0("Compound-target network by module -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(sub, .plot_wrap_width(fig_width, 8)),
      x = NULL, y = NULL
    ) +
    .network_view_theme()

  list(
    plot = p,
    view_tag = if (k > 1) paste0("k", k) else "",
    n_nodes = nrow(nodes), n_edges = nrow(ed),
    data = list(collapse = col, nodes = nodes, centres = cen_df)
  )
}

#' Compound x pathway cells for the pathway view
#'
#' @param ct Compound-target table (`compound_id`, `uniprot_id`), one row
#'   per pair.
#' @param term_sets Named list: term id -> UniProt ids annotated to it.
#' @return `data.frame(ID, compound_id, n_hit, n_term_targets, frac)` --
#'   one row per (term, compound) with at least one hit; `n_term_targets`
#'   counts the in-scope targets annotated to the term and `frac` is the
#'   share of them the compound hits.
#' @keywords internal
.network_view_pathway_cells <- function(ct, term_sets) {
  ct <- unique(ct[, c("compound_id", "uniprot_id")])
  rows <- lapply(names(term_sets), function(id) {
    hit <- ct[ct$uniprot_id %in% term_sets[[id]], , drop = FALSE]
    if (nrow(hit) == 0) return(NULL)
    n_term <- length(unique(hit$uniprot_id))
    tab <- table(hit$compound_id)
    data.frame(ID = id, compound_id = names(tab), n_hit = as.integer(tab), n_term_targets = n_term,
               frac = as.integer(tab) / n_term, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) {
    return(data.frame(ID = character(0), compound_id = character(0), n_hit = integer(0),
                      n_term_targets = integer(0), frac = numeric(0), stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Hierarchical-clustering order of the rows of a matrix (input order when
#' fewer than three rows)
#' @keywords internal
.network_view_hclust_order <- function(m) {
  if (nrow(m) < 3) return(seq_len(nrow(m)))
  stats::hclust(stats::dist(m), method = "average")$order
}

#' Compound x pathway dot view (`view = "pathway"`)
#'
#' @description
#' The `max_pathways` most significant `network_enrichment` terms of the
#' conditions in scope (optionally restricted to `pathway_db`) against the
#' compounds: dot size = number of the compound's targets annotated to the
#' term, dot colour = share of the term's in-scope targets the compound
#' hits. Rows and columns are ordered by hierarchical clustering so
#' compounds with a similar pathway profile sit together. Needs
#' `clusterProfiler` + `org.Hs.eg.db` (UniProt -> Entrez, the id type of
#' `network_enrichment$geneID`).
#' @return `list(plot, view_tag, n_nodes, n_edges, data = list(cells,
#'   terms))`.
#' @keywords internal
.network_view_pathway <- function(proj, ct, conditions, pathway_db = NULL, max_pathways = 30,
                                  title_suffix = "", fig_width = 9) {
  enr <- patliRResults(proj, "network_enrichment")
  if (is.null(enr) || nrow(enr) == 0 || !any(enr$condition %in% conditions)) {
    cli::cli_abort(c(
      "No {.val network_enrichment} rows for condition(s) {.val {conditions}}.",
      "i" = "{.code view = \"pathway\"} needs {.fn network_enrich} to be run first."
    ))
  }
  enr <- enr[enr$condition %in% conditions, , drop = FALSE]
  if (!is.null(pathway_db)) enr <- enr[enr$db %in% pathway_db, , drop = FALSE]
  if (nrow(enr) == 0) cli::cli_abort("No enrichment term of {.arg pathway_db} {.val {pathway_db}} in scope.")
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cli::cli_abort("{.code view = \"pathway\"} needs {.pkg clusterProfiler} and {.pkg org.Hs.eg.db} (Bioconductor) to map UniProt to Entrez ids.")
  }
  if (!is.numeric(max_pathways) || length(max_pathways) != 1L || is.na(max_pathways) || max_pathways < 1) {
    cli::cli_abort("{.arg max_pathways} must be a single number >= 1.")
  }
  enr <- enr[order(enr$p.adjust), , drop = FALSE]
  enr <- enr[!duplicated(enr$ID), , drop = FALSE]
  map <- tryCatch(
    clusterProfiler::bitr(unique(ct$uniprot_id), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE),
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) cli::cli_abort("No target UniProt id mapped to an Entrez id -- nothing to draw.")
  sets <- lapply(enr$geneID, function(s) unique(map$UNIPROT[map$ENTREZID %in% strsplit(s, "/", fixed = TRUE)[[1]]]))
  names(sets) <- enr$ID
  keep <- vapply(sets, length, integer(1)) > 0
  enr <- enr[keep, , drop = FALSE]
  sets <- sets[keep]
  top <- utils::head(seq_len(nrow(enr)), if (is.finite(max_pathways)) max_pathways else nrow(enr))
  enr <- enr[top, , drop = FALSE]
  sets <- sets[top]

  cells <- .network_view_pathway_cells(ct, sets)
  multi_db <- length(unique(enr$db)) > 1
  enr$label <- paste0(if (multi_db) paste0("[", toupper(enr$db), "] ") else "", .plot_truncate(enr$Description, 55))
  enr$label <- make.unique(enr$label, sep = " ")
  cmp_ids <- sort(unique(ct$compound_id))
  cmp_lab <- .plot_truncate(.plot_unique_labels(proj, conditions, cmp_ids, "compound"), 30)

  m <- matrix(0, nrow(enr), length(cmp_ids), dimnames = list(enr$ID, cmp_ids))
  m[cbind(match(cells$ID, enr$ID), match(cells$compound_id, cmp_ids))] <- cells$frac
  row_o <- .network_view_hclust_order(m)
  col_o <- .network_view_hclust_order(t(m))
  cells$term <- factor(enr$label[match(cells$ID, enr$ID)], levels = rev(enr$label[row_o]))
  cells$compound <- factor(cmp_lab[match(cells$compound_id, cmp_ids)], levels = cmp_lab[col_o])

  p <- ggplot2::ggplot(cells, ggplot2::aes(x = .data$compound, y = .data$term)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_hit, fill = .data$frac), shape = 21, colour = "grey30", stroke = 0.2) +
    ggplot2::scale_fill_viridis_c(option = "D", direction = -1, limits = c(0, 1), name = "Share of the\nterm's targets") +
    ggplot2::scale_size_area(max_size = 5.5, name = "Targets of the\ncompound in term") +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "grey60"))) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(
      title = .plot_wrap(paste0("Compounds x enriched pathways -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(sprintf(paste0(
        "The %d most significant %s terms (lowest adjusted p in scope; rows and columns clustered) against the compounds. ",
        "Dot size = how many of the compound's predicted targets are annotated to the term; colour = which share of ",
        "the term's targets in this network the compound hits."
      ), nrow(enr), paste(toupper(sort(unique(enr$db))), collapse = "/")),
      .plot_wrap_width(fig_width, 8)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1, size = 7.5),
      axis.text.y = ggplot2::element_text(size = 7.5),
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7)
    )

  db_tag <- if (!is.null(pathway_db)) paste0("_", gsub("[^A-Za-z0-9]+", "-", paste(sort(pathway_db), collapse = "+"))) else ""
  list(
    plot = p,
    view_tag = paste0("top", nrow(enr), db_tag),
    n_nodes = nrow(enr) + length(cmp_ids), n_edges = nrow(cells),
    data = list(cells = cells, terms = enr)
  )
}
