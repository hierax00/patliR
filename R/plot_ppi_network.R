#' @include AllGenerics.R internal.R network_build.R network_bowtie.R plot-helpers.R plot_network_layers.R
NULL

## The "compound -> protein -> protein (PPI) -> disease" picture: every
## protein the extract's compounds are predicted to hit, joined by the STRING
## protein-protein interactions among those proteins, with a disease gene set
## (disease_genes_fetch()/disease_genes_import()) overlaid. The other network
## plots draw either the bipartite compound-target graph
## (plot_network_layers()) or STRING distances as numbers
## (plot_proximity()/plot_synergy()); none shows which of the hit proteins
## actually talk to each other and where the disease genes sit among them.
##
## STRING is read straight from its flat files (links + aliases), cached under
## cacheDir(proj)/stringdb/ -- the same directory and files network_proximity()
## (through STRINGdb) and network_bowtie() already use -- so no STRINGdb object
## (and no STRING API version check over the network) is needed here.

#' Protein-protein interaction network of an extract's predicted targets,
#' with a disease overlay
#'
#' @description
#' Draws the proteins the compounds of one condition (or several, pooled) are
#' predicted to hit (`network_edges`), joined by the STRING protein-protein
#' interactions among them, and marks which of them belong to a disease gene
#' set (`disease_genes`). It is the step between the compound-target network
#' ([plot_network_layers()]) and the disease: which hit proteins interact
#' with each other, and whether the disease genes form a connected block
#' among them.
#'
#' With a few hundred targets the full graph is unreadable, so the figure is
#' simplified by default:
#'
#' * Targets are ranked by (1) being a disease gene, (2) the number of the
#'   extract's compounds predicted to hit them, (3) their STRING degree among
#'   the extract's targets, and the top `max_targets` are drawn.
#' * A target with no STRING partner (at `score_threshold`) among the drawn
#'   proteins is not drawn; its slot goes to the next target in the ranking.
#'   The subtitle states how many targets were left out and why.
#' * Only `top_n_labels` targets are labelled (gene symbols): the
#'   highest-ranked disease genes (at most 60% of the labels) and the
#'   best-connected other targets. Every node keeps its full name in the
#'   `ggiraph` tooltip.
#'
#' @section Views:
#' * `"compound_ppi"` (default): the targets and their PPI edges in the
#'   centre, the compounds as labelled diamonds on an outer ring, each
#'   joined to the drawn targets it hits by a thin line. A compound sits on
#'   the side of the ring nearest to its targets.
#' * `"ppi"`: only the targets and their PPI edges (with the disease
#'   overlay). Use it when the compound lines crowd the figure.
#'
#' @section Colours and sizes:
#' `colour_by = "disease"` fills disease genes in reddish purple and other
#' targets in sky blue. `colour_by = "module"` fills targets by their
#' [network_module_robustness()] module (`network_module_membership`,
#' colour-blind-safe Okabe-Ito palette, the same colours as
#' [plot_network_layers()]). In both cases a disease gene also gets a dark
#' purple ring, and an interaction between two disease genes is purple.
#' Node size is the number of compounds hitting the protein
#' (`size_by = "compounds"`) or its STRING degree in the drawn graph
#' (`size_by = "degree"`). PPI edge opacity follows the STRING combined
#' score.
#'
#' @section GPCRs (metabotropic receptors) in red:
#' Many plant metabolites are predicted to bind G protein-coupled
#' (metabotropic) receptors, and such receptors are expected to be hit by
#' several compounds. With `highlight = "gpcr"` (default) these proteins are
#' drawn as squares with a red outline, their labels in red, and every
#' interaction with at least one GPCR endpoint as a thicker red line; the
#' subtitle counts them. Red is used only for this overlay (the disease
#' overlay is purple, a fill plus a ring), so both can be shown together.
#'
#' A protein is a GPCR when `org.Hs.eg.db` annotates its UniProt accession
#' (or, failing that, its gene symbol) to GO:0004930 "G protein-coupled
#' receptor activity" or any descendant term, which includes the
#' metabotropic glutamate (mGluR, GO:0098988/GO:0001640) and GABA-B
#' (GO:0004965) receptors. GO:0008066 "glutamate receptor activity" is not
#' used because it also covers ionotropic receptors. Electronic (IEA)
#' annotations are kept: without them GO:0004930 covers only about 390 human
#' genes, fewer than half of the known human GPCRs. Genes supported only by
#' author statements (TAS/NAS/IC) are dropped, because those are mostly
#' pathway-level mis-annotations of non-receptors (e.g. PDGFRB, PPARG).
#' `highlight_set` adds any
#' UniProt accessions of your own (drawn the same way); with
#' `highlight = "none"` and `highlight_set` given, only your set is
#' highlighted. The classification of every predicted target is returned as
#' `attr(., "gpcr")`.
#'
#' @section Disease-gene enrichment:
#' The summary reports a one-sided hypergeometric test of the overlap
#' between the extract's targets and the disease genes, both mapped to STRING,
#' with the STRING proteome of `species` (every protein in STRING's alias
#' table) as the universe. Target prediction tools only cover a part of the
#' proteome (proteins with known ligands), which this universe does not
#' account for, so read the p-value as descriptive, not as a calibrated test.
#'
#' @section STRING data:
#' STRING's `protein.links` and `protein.aliases` flat files for `species`
#' and `version` are read from `cacheDir(proj)/stringdb/` (downloaded once
#' from stringdb-downloads.org when absent -- the same files
#' [network_proximity()] uses). The links filtered at `score_threshold` are
#' cached there as `ppi_links_<species>_<version>_<threshold>.rds`, so only
#' the first call pays for reading the full links file. UniProt accessions
#' that STRING does not know are counted as unmapped and never drawn.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease `NULL` (default, no overlay) or one `disease_id` present in
#'   `patliRResults(proj, "disease_genes")`, e.g. `"MONDO_0005044"`.
#' @param view `"compound_ppi"` (default) or `"ppi"`. See the section on views.
#' @param colour_by `"disease"` (default) or `"module"`. See the section on
#'   colours.
#' @param size_by `"compounds"` (default) or `"degree"`.
#' @param max_targets Maximum number of targets drawn (default `60`); `Inf`
#'   draws every target that has a STRING partner.
#' @param top_n_labels Number of targets that get a text label (default
#'   `30`); `0` labels none. Compounds are always labelled.
#' @param score_threshold Minimum STRING combined score (0-1000) of a PPI
#'   edge. Default `700` ("high confidence").
#' @param disease_partners Integer, default `0`. When `> 0` and `disease` is
#'   given, up to this many disease genes that are *not* predicted targets but
#'   interact with at least two drawn targets are added as hollow red circles
#'   -- the disease proteins the extract could reach through one PPI step.
#' @param disease_node Logical, default `FALSE`. When `TRUE` (and `disease` is
#'   given) a single labelled disease node is drawn beside the network and
#'   joined by dashed lines to every drawn disease gene.
#' @param highlight `"gpcr"` (default) marks G protein-coupled (metabotropic)
#'   receptors and their interactions in red; `"none"` turns this off. See
#'   the section on GPCRs.
#' @param highlight_set `NULL` (default) or a character vector of UniProt
#'   accessions to highlight in red as well (or instead, with
#'   `highlight = "none"`).
#' @param layout `"fr"` (Fruchterman-Reingold, default) or `"kk"`
#'   (Kamada-Kawai) for the PPI graph.
#' @param seed Integer seed for the layout (default `1`), so the figure is
#'   reproducible.
#' @param species NCBI taxonomy ID. Default `9606` (human).
#' @param version STRING version. Default `"12.0"`.
#' @param engine `"static"` (default) or `"ggiraph"` (tooltips on hover).
#'
#' @return A `ggplot` (or `girafe`) object. `attr(., "summary")` is a
#'   one-row `data.frame` (`n_targets`, `n_targets_mapped`, `n_ppi_edges`
#'   among all targets, `n_no_partner`, `n_targets_drawn`,
#'   `n_ppi_edges_drawn`, `n_compounds_drawn`, `n_disease_genes`,
#'   `n_disease_genes_mapped`, `n_disease_hit`, `n_disease_hit_drawn`,
#'   `n_disease_overlap_string`, `n_universe`, `fold_enrichment`,
#'   `p_hypergeom`, `highlight`, `n_gpcr_targets`, `n_gpcr_drawn`,
#'   `n_gpcr_edges_drawn`, `score_threshold`); `attr(., "gpcr")` a
#'   `data.frame` with one row per predicted target (`uniprot_id`,
#'   `gene_symbol`, `is_highlighted`, `source`, GO `evidence` codes,
#'   `n_compounds`, `drawn`); and `attr(., "data")` a list with the drawn
#'   `nodes`, `ppi_edges` (with a `gpcr_edge` flag), `compound_edges` and
#'   `compounds`. With `save = TRUE` a PNG is written, logged to
#'   `patliRResults(proj, "ppi_network_plot_log")` keyed by
#'   `(condition, disease_id, view, colour_by, variant)` -- `variant` is
#'   `"default"` or names the non-default overlays (`disease_partners`,
#'   `disease_node`, `highlight = "none"`, `highlight_set`), which also go
#'   into the file name -- and `attr(., "proj")` holds
#'   the updated project.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_load("my_project")
#' # needs network_build() and disease_genes_fetch() first
#' p <- plot_ppi_network(proj, condition = "EFLO-S", disease = "MONDO_0005044")
#' attr(p, "summary")
#' plot_ppi_network(proj, condition = "EFLO-S", disease = "MONDO_0005044",
#'                  view = "ppi", colour_by = "module", save = FALSE)
#' }
#'
#' @seealso [plot_network_layers()], [network_proximity()],
#'   [disease_genes_fetch()], [network_module_robustness()]
#' @export
plot_ppi_network <- function(proj, condition = NULL, disease = NULL,
                             view = c("compound_ppi", "ppi"),
                             colour_by = c("disease", "module"),
                             size_by = c("compounds", "degree"),
                             max_targets = 60, top_n_labels = 30, score_threshold = 700,
                             disease_partners = 0, disease_node = FALSE,
                             highlight = c("gpcr", "none"), highlight_set = NULL,
                             layout = c("fr", "kk"), seed = 1,
                             species = 9606, version = "12.0",
                             engine = c("static", "ggiraph"),
                             save = TRUE, out_dir = NULL, width = 10, height = 9, dpi = 300) {
  stopifnot(is(proj, "PatliRProject"))
  view <- match.arg(view)
  colour_by <- match.arg(colour_by)
  size_by <- match.arg(size_by)
  layout <- match.arg(layout)
  engine <- match.arg(engine)
  highlight <- match.arg(highlight)
  if (!is.null(highlight_set) && (!is.character(highlight_set) || anyNA(highlight_set))) {
    cli::cli_abort("{.arg highlight_set} must be {.code NULL} or a character vector of UniProt accessions.")
  }
  if (!is.numeric(max_targets) || length(max_targets) != 1 || is.na(max_targets) || max_targets < 2) {
    cli::cli_abort("{.arg max_targets} must be a single number >= 2 (or {.code Inf}).")
  }
  if (!is.numeric(top_n_labels) || length(top_n_labels) != 1 || is.na(top_n_labels) || top_n_labels < 0) {
    cli::cli_abort("{.arg top_n_labels} must be a single number >= 0.")
  }
  if (!is.numeric(score_threshold) || length(score_threshold) != 1 || is.na(score_threshold) ||
      score_threshold < 0 || score_threshold > 1000) {
    cli::cli_abort("{.arg score_threshold} must be a single number between 0 and 1000 (STRING combined score).")
  }
  if (!is.numeric(disease_partners) || length(disease_partners) != 1 || is.na(disease_partners) || disease_partners < 0) {
    cli::cli_abort("{.arg disease_partners} must be a single number >= 0.")
  }
  if (!is.logical(disease_node) || length(disease_node) != 1 || is.na(disease_node)) {
    cli::cli_abort("{.arg disease_node} must be {.code TRUE} or {.code FALSE}.")
  }
  if (!is.null(disease) && (!is.character(disease) || length(disease) != 1 || is.na(disease) || !nzchar(disease))) {
    cli::cli_abort("{.arg disease} must be {.code NULL} or a single {.field disease_id}.")
  }
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  ## ---- targets of the extract ---------------------------------------------
  edges_all <- patliRResults(proj, "network_edges")
  ct <- unique(edges_all[edges_all$condition %in% conditions & !is.na(edges_all$uniprot_id),
                         c("compound_id", "uniprot_id"), drop = FALSE])
  if (nrow(ct) == 0) cli::cli_abort("No compound-target edge in condition(s) {.val {conditions}}.")
  n_cmp <- table(ct$uniprot_id)
  targets <- data.frame(uniprot_id = names(n_cmp), n_compounds = as.integer(n_cmp), stringsAsFactors = FALSE)

  ## ---- disease genes --------------------------------------------------------
  dg <- NULL
  disease_name <- NA_character_
  if (!is.null(disease)) {
    dg_all <- patliRResults(proj, "disease_genes")
    if (is.null(dg_all) || nrow(dg_all) == 0) {
      cli::cli_abort(c(
        "No {.val disease_genes} entry in {.arg proj}.",
        "i" = "Run {.fn disease_genes_fetch} or {.fn disease_genes_import} first."
      ))
    }
    dg <- dg_all[dg_all$disease_id == disease & !is.na(dg_all$uniprot_id), , drop = FALSE]
    if (nrow(dg) == 0) {
      cli::cli_abort(c(
        "No gene associated with {.val {disease}} in {.val disease_genes}.",
        "i" = "Available disease IDs: {.val {unique(dg_all$disease_id)}}."
      ))
    }
    if ("disease_name" %in% names(dg)) {
      nm <- stats::na.omit(dg$disease_name)
      if (length(nm) > 0 && nzchar(nm[1])) disease_name <- as.character(nm[1])
    }
    dg <- dg[!duplicated(dg$uniprot_id), , drop = FALSE]
  }
  disease_ids <- if (is.null(dg)) character(0) else unique(dg$uniprot_id)
  targets$is_disease <- targets$uniprot_id %in% disease_ids

  ## ---- gene symbols -----------------------------------------------------------
  targets$label <- .ppi_symbols(proj, conditions, targets$uniprot_id, dg)

  ## ---- GPCR / highlighted proteins --------------------------------------------
  hl <- .ppi_highlight_table(targets$uniprot_id, targets$label, highlight, highlight_set)
  targets$is_gpcr <- hl$is_highlighted[match(targets$uniprot_id, hl$uniprot_id)]

  ## ---- STRING ---------------------------------------------------------------
  string <- .ppi_string_network(proj, unique(c(targets$uniprot_id, disease_ids)),
                                species = species, version = version, score_threshold = score_threshold)
  mapped <- string$map
  ppi_all <- string$edges
  targets$string_id <- unname(mapped[targets$uniprot_id])
  is_target <- ppi_all$uniprot_a %in% targets$uniprot_id & ppi_all$uniprot_b %in% targets$uniprot_id
  ppi_t <- ppi_all[is_target, , drop = FALSE]
  deg_all <- table(factor(c(ppi_t$uniprot_a, ppi_t$uniprot_b), levels = targets$uniprot_id))
  targets$ppi_degree_all <- as.integer(deg_all[targets$uniprot_id])

  ## ---- enrichment summary ---------------------------------------------------
  t_mapped <- unique(stats::na.omit(targets$string_id))
  d_mapped <- unique(stats::na.omit(unname(mapped[disease_ids])))
  enr <- .ppi_hypergeom(t_mapped, d_mapped, string$n_universe)

  ## ---- selection ------------------------------------------------------------
  sel <- .ppi_select_targets(targets, ppi_t, max_targets)
  drawn <- sel$drawn
  if (length(drawn) < 2) {
    cli::cli_abort(c(
      "No STRING interaction (combined score >= {score_threshold}) among the predicted targets of {.val {scope_label}}.",
      "i" = "{sum(!is.na(targets$string_id))} of {nrow(targets)} target(s) mapped to STRING; lower {.arg score_threshold} to see weaker interactions."
    ))
  }
  nodes <- targets[match(drawn, targets$uniprot_id), , drop = FALSE]
  nodes$kind <- "target"
  ppi <- ppi_t[ppi_t$uniprot_a %in% drawn & ppi_t$uniprot_b %in% drawn, , drop = FALSE]

  ## disease genes that are not targets but touch >= 2 drawn targets
  n_partners <- 0L
  if (!is.null(dg) && disease_partners > 0) {
    cand <- .ppi_disease_partners(ppi_all, drawn, setdiff(disease_ids, targets$uniprot_id), dg, disease_partners)
    if (nrow(cand) > 0) {
      n_partners <- nrow(cand)
      p_label <- .ppi_symbols(proj, conditions, cand$uniprot_id, dg)
      p_hl <- .ppi_highlight_table(cand$uniprot_id, p_label, highlight, highlight_set)
      nodes <- rbind(nodes, data.frame(
        uniprot_id = cand$uniprot_id, n_compounds = 0L, is_disease = TRUE, label = p_label,
        is_gpcr = p_hl$is_highlighted,
        string_id = unname(mapped[cand$uniprot_id]), ppi_degree_all = NA_integer_,
        kind = "disease_partner", stringsAsFactors = FALSE
      ))
      keep <- ppi_all$uniprot_a %in% nodes$uniprot_id & ppi_all$uniprot_b %in% nodes$uniprot_id
      ppi <- ppi_all[keep, , drop = FALSE]
      ## partner-partner edges would pull the partners into a block of their
      ## own; only their links to the drawn targets are the point
      both_partner <- ppi$uniprot_a %in% cand$uniprot_id & ppi$uniprot_b %in% cand$uniprot_id
      ppi <- ppi[!both_partner, , drop = FALSE]
    }
  }
  deg_drawn <- table(factor(c(ppi$uniprot_a, ppi$uniprot_b), levels = nodes$uniprot_id))
  nodes$degree <- as.integer(deg_drawn[nodes$uniprot_id])
  nodes$rank <- match(nodes$uniprot_id, drawn)

  ppi$gpcr_edge <- nodes$is_gpcr[match(ppi$uniprot_a, nodes$uniprot_id)] |
    nodes$is_gpcr[match(ppi$uniprot_b, nodes$uniprot_id)]

  ## ---- layout ---------------------------------------------------------------
  restore_rng <- .with_seed(as.integer(seed))
  xy <- .ppi_layout(nodes$uniprot_id, ppi, layout)
  restore_rng()
  nodes$x <- xy[, 1]
  nodes$y <- xy[, 2]

  ## ---- compounds ------------------------------------------------------------
  cmp <- NULL
  ce <- ct[ct$uniprot_id %in% nodes$uniprot_id[nodes$kind == "target"], , drop = FALSE]
  n_cmp_hidden <- 0L
  if (view == "compound_ppi") {
    all_cmp <- unique(ct$compound_id)
    cmp <- .ppi_compound_ring(ce, nodes, radius = 1.2)
    n_cmp_hidden <- length(setdiff(all_cmp, cmp$compound_id))
    cmp$n_targets <- as.integer(table(ct$compound_id)[cmp$compound_id])
    cmp$label <- .plot_unique_labels(proj, conditions, cmp$compound_id, "compound")
  }

  ## ---- summary --------------------------------------------------------------
  summary <- data.frame(
    condition = scope_label, disease_id = if (is.null(disease)) "none" else disease,
    view = view, n_targets = nrow(targets), n_targets_mapped = sum(!is.na(targets$string_id)),
    n_ppi_edges = nrow(ppi_t), n_no_partner = sel$n_no_partner,
    n_targets_drawn = sum(nodes$kind == "target"),
    n_ppi_edges_drawn = sum(ppi$uniprot_a %in% drawn & ppi$uniprot_b %in% drawn),
    n_disease_partners = n_partners,
    n_compounds_drawn = if (is.null(cmp)) 0L else nrow(cmp),
    n_disease_genes = length(disease_ids), n_disease_genes_mapped = length(d_mapped),
    n_disease_hit = sum(targets$is_disease),
    n_disease_hit_drawn = sum(nodes$is_disease & nodes$kind == "target"),
    n_disease_overlap_string = enr$k,
    n_universe = string$n_universe, fold_enrichment = enr$fold, p_hypergeom = enr$p,
    highlight = .ppi_highlight_name(highlight, highlight_set),
    n_gpcr_targets = sum(targets$is_gpcr), n_gpcr_drawn = sum(nodes$is_gpcr),
    n_gpcr_edges_drawn = sum(ppi$gpcr_edge),
    score_threshold = score_threshold, stringsAsFactors = FALSE
  )

  ## ---- plot -----------------------------------------------------------------
  p <- .ppi_ggplot(
    nodes = nodes, ppi = ppi, cmp = cmp, ce = ce, proj = proj, conditions = conditions,
    disease = disease, disease_name = disease_name, colour_by = colour_by, size_by = size_by,
    top_n_labels = top_n_labels, disease_node = disease_node && !is.null(disease),
    score_threshold = score_threshold, engine = engine, summary = summary,
    scope_label = scope_label, n_cmp_hidden = n_cmp_hidden, fig_width = width,
    hl_name = summary$highlight
  )

  disease_tag <- if (is.null(disease)) "none" else gsub("[^A-Za-z0-9_.-]+", "_", disease)
  ## non-default overlays get their own file and log row
  variant <- c(
    if (n_partners > 0) paste0("partners", n_partners) else NULL,
    if (disease_node && !is.null(disease)) "disease_node" else NULL,
    if (highlight == "none") "no_gpcr" else NULL,
    if (length(highlight_set) > 0) paste0("set_", substr(rlang::hash(sort(unique(highlight_set))), 1, 8)) else NULL
  )
  variant <- if (length(variant) == 0) "default" else paste(variant, collapse = "_")
  proj <- .plot_log_backfill(proj, "ppi_network_plot_log", "variant", "default")
  filename <- paste0("ppi_network_", scope_label, "_", disease_tag, "_", view,
                     if (colour_by == "module") "_module" else "",
                     if (variant != "default") paste0("_", variant) else "", ".png")
  log_row <- cbind(summary[, c("condition", "disease_id", "view")],
                   data.frame(colour_by = colour_by, variant = variant, path = NA_character_, stringsAsFactors = FALSE),
                   summary[, setdiff(names(summary), c("condition", "disease_id", "view"))])
  result <- .plot_finish(
    proj, p,
    name = "ppi_network_plot_log", filename = filename, log_row = log_row,
    key_cols = c("condition", "disease_id", "view", "colour_by", "variant"),
    engine = engine, save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
  attr(result, "summary") <- summary
  hl$gene_symbol <- targets$label[match(hl$uniprot_id, targets$uniprot_id)]
  hl$n_compounds <- targets$n_compounds[match(hl$uniprot_id, targets$uniprot_id)]
  hl$drawn <- hl$uniprot_id %in% nodes$uniprot_id[nodes$kind == "target"]
  attr(result, "gpcr") <- hl[, c("uniprot_id", "gene_symbol", "is_highlighted", "source", "evidence",
                                 "n_compounds", "drawn")]
  attr(result, "data") <- list(
    nodes = nodes, ppi_edges = ppi,
    compound_edges = if (is.null(cmp)) ce[0, , drop = FALSE] else ce,
    compounds = cmp
  )
  result
}

#' STRING PPI edges among a set of UniProt accessions
#'
#' @description
#' Maps `uniprot_ids` to STRING IDs through STRING's alias table
#' ([.network_stringdb_aliases_map()]) and returns the STRING links (at
#' `score_threshold`) among them, translated back to UniProt accessions.
#' When two accessions map to the same STRING protein, each gets the
#' protein's edges.
#' @return `list(edges = data.frame(uniprot_a, uniprot_b, score),
#'   map = named character (uniprot -> STRING id, `NA` when unmapped),
#'   n_universe = number of STRING proteins in the alias table)`.
#' @keywords internal
.ppi_string_network <- function(proj, uniprot_ids, species, version, score_threshold) {
  uniprot_ids <- unique(as.character(uniprot_ids))
  aliases <- .network_stringdb_aliases_map(proj, species, version)
  hit <- aliases[aliases$alias %in% uniprot_ids, , drop = FALSE]
  hit <- hit[!duplicated(hit$alias), , drop = FALSE]
  map <- stats::setNames(hit$string_protein_id[match(uniprot_ids, hit$alias)], uniprot_ids)
  n_universe <- length(unique(aliases$string_protein_id))

  links <- .ppi_string_links(proj, species, version, score_threshold)
  ids <- unique(stats::na.omit(map))
  links <- links[links$protein1 %in% ids & links$protein2 %in% ids & links$combined_score >= score_threshold, , drop = FALSE]
  empty <- data.frame(uniprot_a = character(0), uniprot_b = character(0), score = numeric(0), stringsAsFactors = FALSE)
  if (nrow(links) == 0) return(list(edges = empty, map = map, n_universe = n_universe))

  ## STRING id -> every accession mapped to it (a merge, so an accession
  ## sharing its STRING protein with another one gets the same edges)
  map_df <- data.frame(uniprot = names(map)[!is.na(map)], string = unname(map[!is.na(map)]), stringsAsFactors = FALSE)
  edges <- merge(links, stats::setNames(map_df, c("uniprot_a", "protein1")), by = "protein1")
  edges <- merge(edges, stats::setNames(map_df, c("uniprot_b", "protein2")), by = "protein2")
  edges <- data.frame(uniprot_a = edges$uniprot_a, uniprot_b = edges$uniprot_b,
                      score = as.numeric(edges$combined_score), stringsAsFactors = FALSE)
  edges <- edges[edges$uniprot_a != edges$uniprot_b, , drop = FALSE]
  ## undirected, one row per pair (highest score kept)
  swap <- edges$uniprot_a > edges$uniprot_b
  tmp <- edges$uniprot_a[swap]
  edges$uniprot_a[swap] <- edges$uniprot_b[swap]
  edges$uniprot_b[swap] <- tmp
  edges <- edges[order(-edges$score), , drop = FALSE]
  edges <- edges[!duplicated(paste(edges$uniprot_a, edges$uniprot_b)), , drop = FALSE]
  edges <- edges[order(edges$uniprot_a, edges$uniprot_b), , drop = FALSE]
  rownames(edges) <- NULL
  list(edges = edges, map = map, n_universe = n_universe)
}

#' STRING `protein.links` at a score threshold, cached
#'
#' @description
#' Reads `<species>.protein.links.v<version>.txt.gz` from
#' `cacheDir(proj)/stringdb/` (downloading it once when absent, like
#' [network_proximity()]'s `STRINGdb` object does), keeps the links with
#' `combined_score >= score_threshold`, one row per undirected pair, and
#' caches that table as `ppi_links_<species>_<version>_<threshold>.rds`. A
#' cached table at a lower threshold is filtered instead of re-reading the
#' flat file. Safe to delete: the flat file is the real cache.
#' @return `data.frame(protein1, protein2, combined_score)`.
#' @keywords internal
.ppi_string_links <- function(proj, species, version, score_threshold) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  prefix <- paste0("ppi_links_", species, "_", version, "_")
  cache_rds <- file.path(dir, paste0(prefix, score_threshold, ".rds"))
  if (file.exists(cache_rds)) return(readRDS(cache_rds))

  cached <- list.files(dir, pattern = paste0("^", gsub(".", "[.]", prefix, fixed = TRUE), "[0-9.]+[.]rds$"))
  cached_thr <- suppressWarnings(as.numeric(sub("[.]rds$", "", substring(cached, nchar(prefix) + 1))))
  usable <- which(!is.na(cached_thr) & cached_thr <= score_threshold)
  if (length(usable) > 0) {
    best <- usable[which.max(cached_thr[usable])]
    links <- readRDS(file.path(dir, cached[best]))
    return(links[links$combined_score >= score_threshold, , drop = FALSE])
  }

  raw_gz <- file.path(dir, paste0(species, ".protein.links.v", version, ".txt.gz"))
  if (!file.exists(raw_gz)) {
    url <- paste0("https://stringdb-downloads.org/download/protein.links.v", version, "/",
                  species, ".protein.links.v", version, ".txt.gz")
    utils::download.file(url, raw_gz, mode = "wb", quiet = FALSE)
  }
  links <- utils::read.table(gzfile(raw_gz), header = TRUE, sep = " ", stringsAsFactors = FALSE,
                             colClasses = c("character", "character", "integer"), comment.char = "", quote = "")
  required <- c("protein1", "protein2", "combined_score")
  if (!all(required %in% names(links))) {
    cli::cli_abort(c(
      "STRING's {.file protein.links} file for species {.val {species}}/v{version} lacks column(s) {.val {setdiff(required, names(links))}}.",
      "i" = "Columns found: {.val {names(links)}}."
    ))
  }
  links <- links[links$combined_score >= score_threshold & links$protein1 < links$protein2, required, drop = FALSE]
  rownames(links) <- NULL
  saveRDS(links, cache_rds)
  links
}

#' One-sided hypergeometric test of a target set / disease set overlap
#' @return `list(k, fold, p)`; `NA`s when either set is empty.
#' @keywords internal
.ppi_hypergeom <- function(targets, disease, n_universe) {
  n <- length(targets)
  K <- length(disease)
  k <- length(intersect(targets, disease))
  if (n == 0 || K == 0 || is.na(n_universe) || n_universe < max(n, K)) {
    return(list(k = k, fold = NA_real_, p = NA_real_))
  }
  list(k = k, fold = (k / n) / (K / n_universe),
       p = stats::phyper(k - 1, K, n_universe - K, n, lower.tail = FALSE))
}

#' Gene symbols for UniProt accessions
#'
#' @description
#' [.plot_label_nodes()] (org.Hs.eg.db, when installed), with the
#' `gene_symbol` column of the disease gene table as a fallback for
#' accessions it leaves unresolved.
#' @return Character vector, same length as `ids` (the accession itself when
#'   nothing resolves).
#' @keywords internal
.ppi_symbols <- function(proj, conditions, ids, dg = NULL) {
  if (length(ids) == 0) return(character(0))
  ## bitr()'s "x% of input gene IDs are fail to map" is expected here
  ## (unresolved IDs keep their accession) and only noise
  lab <- suppressWarnings(suppressMessages(.plot_label_nodes(proj, conditions, ids, "target")))
  if (!is.null(dg) && "gene_symbol" %in% names(dg)) {
    sym <- stats::setNames(as.character(dg$gene_symbol), dg$uniprot_id)[ids]
    unresolved <- lab == ids & !is.na(sym) & nzchar(sym)
    lab[unresolved] <- sym[unresolved]
  }
  unname(lab)
}

#' GO terms that define the GPCR (metabotropic receptor) highlight
#'
#' @description
#' `GO:0004930` "G protein-coupled receptor activity" and all its descendants
#' (queried through the `GOALL` key, which already includes every child
#' term). Metabotropic glutamate (`GO:0098988`, `GO:0001640`) and GABA-B
#' (`GO:0004965`) receptor activities are descendants of `GO:0004930`; they
#' are listed explicitly so the definition does not depend on the GO
#' release keeping them there. `GO:0008066` "glutamate receptor activity" is
#' deliberately *not* used: it also covers the ionotropic AMPA/NMDA/kainate
#' receptors.
#' @keywords internal
.ppi_gpcr_terms <- function() c("GO:0004930", "GO:0098988", "GO:0001640", "GO:0004965")

#' Human GPCRs (UniProt accessions and gene symbols) from org.Hs.eg.db
#'
#' @description
#' Every UniProt accession and gene symbol annotated in `org.Hs.eg.db` to one
#' of [.ppi_gpcr_terms()] (or a descendant), with the evidence codes
#' collapsed per accession. IEA (electronic) annotations are **kept**:
#' without them `GO:0004930` covers about 390 human genes, less than half of
#' the ~800 known human GPCRs -- most of the IEA-only genes are olfactory,
#' taste and orphan receptors whose 7-transmembrane assignment (InterPro) is
#' reliable. Genes whose only support is an author statement or curator
#' inference (TAS, NAS, IC, ND) are dropped: in practice these are
#' pathway-level mappings of non-receptors such as PDGFRB, PPARG or NPY
#' (17 genes in the Bioconductor 3.21 `org.Hs.eg.db`; LGR4/LGR5 are the only
#' true GPCRs lost, and can be added with `highlight_set`). Queried with [clusterProfiler::bitr()] (`fromType = "GOALL"`);
#' the result is cached for the session.
#' @return `data.frame(uniprot_id, symbol, evidence)`, or `NULL` when
#'   `clusterProfiler`/`org.Hs.eg.db` are not installed.
#' @keywords internal
.ppi_gpcr_go <- function() {
  if (!is.null(.ppi_cache$gpcr)) return(.ppi_cache$gpcr)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(NULL)
  }
  x <- tryCatch(
    suppressWarnings(suppressMessages(clusterProfiler::bitr(
      .ppi_gpcr_terms(), fromType = "GOALL", toType = c("UNIPROT", "SYMBOL", "EVIDENCEALL"),
      OrgDb = "org.Hs.eg.db", drop = TRUE
    ))),
    error = function(e) NULL
  )
  if (is.null(x) || nrow(x) == 0) return(NULL)
  x <- x[!is.na(x$UNIPROT) | !is.na(x$SYMBOL), , drop = FALSE]
  ## a gene supported only by author statements / curator inference is
  ## dropped: those are mostly pathway-level mappings of non-receptors
  ## (PDGFRB, PPARG, NPY, GNAT2, PAX8 ... in the 2024 org.Hs.eg.db)
  weak <- c("TAS", "NAS", "IC", "ND")
  strong <- unique(x$SYMBOL[!x$EVIDENCEALL %in% weak])
  x <- x[x$SYMBOL %in% strong, , drop = FALSE]
  ev <- tapply(x$EVIDENCEALL, x$UNIPROT, function(e) paste(sort(unique(stats::na.omit(e))), collapse = ","))
  sym <- tapply(x$SYMBOL, x$UNIPROT, function(s) s[1])
  out <- data.frame(uniprot_id = names(ev), symbol = as.character(sym[names(ev)]),
                    evidence = as.character(ev), stringsAsFactors = FALSE)
  ## symbols with no accession row still identify a GPCR by name
  extra <- setdiff(stats::na.omit(unique(x$SYMBOL)), out$symbol)
  if (length(extra) > 0) {
    out <- rbind(out, data.frame(uniprot_id = NA_character_, symbol = extra, evidence = NA_character_,
                                 stringsAsFactors = FALSE))
  }
  .ppi_cache$gpcr <- out
  out
}

#' @keywords internal
.ppi_cache <- new.env(parent = emptyenv())

#' Which proteins are highlighted (GPCRs and/or a user set)
#'
#' @description
#' `highlight = "gpcr"` marks a protein whose accession -- or, failing that,
#' gene symbol -- is in [.ppi_gpcr_go()]; `highlight_set` adds the listed
#' accessions. With `highlight = "none"` and no `highlight_set` nothing is
#' marked.
#' @return `data.frame(uniprot_id, is_highlighted, source, evidence)`, one
#'   row per `ids`. `source` is `"GO:0004930"` (GO annotation),
#'   `"highlight_set"` or `NA`.
#' @keywords internal
.ppi_highlight_table <- function(ids, symbols, highlight, highlight_set = NULL) {
  out <- data.frame(uniprot_id = ids, is_highlighted = FALSE, source = NA_character_,
                    evidence = NA_character_, stringsAsFactors = FALSE)
  if (highlight == "gpcr") {
    go <- .ppi_gpcr_go()
    if (is.null(go)) {
      cli::cli_warn(c(
        "{.code highlight = \"gpcr\"} needs {.pkg org.Hs.eg.db} (Bioconductor), not installed; no protein is highlighted.",
        "i" = "Pass the accessions yourself with {.arg highlight_set}."
      ))
    } else {
      by_acc <- match(ids, go$uniprot_id)
      by_sym <- match(symbols, go$symbol)
      hit <- !is.na(by_acc) | !is.na(by_sym)
      out$is_highlighted[hit] <- TRUE
      out$source[hit] <- "GO:0004930"
      out$evidence[hit] <- ifelse(!is.na(by_acc[hit]), go$evidence[by_acc[hit]], "symbol match")
    }
  }
  if (length(highlight_set) > 0) {
    user <- ids %in% highlight_set
    out$source[user] <- ifelse(out$is_highlighted[user], "GO:0004930+highlight_set", "highlight_set")
    out$is_highlighted[user] <- TRUE
  }
  out
}

#' Legend/log name of the highlighted set
#' @keywords internal
.ppi_highlight_name <- function(highlight, highlight_set) {
  user <- length(highlight_set) > 0
  if (highlight == "gpcr" && user) "GPCR / highlighted" else if (highlight == "gpcr") "GPCR" else if (user) "Highlighted" else "none"
}

#' Pick the targets to draw
#'
#' @description
#' Targets are ordered by disease-gene status, number of compounds, STRING
#' degree among all targets, then UniProt ID. Targets without any STRING
#' partner among the extract's targets are never drawn. The top
#' `max_targets` are taken; any of them with no partner *among the selected*
#' is replaced by the next in the order, repeated until every selected target
#' has a partner or the candidates run out.
#' @param targets `data.frame(uniprot_id, n_compounds, is_disease, ppi_degree_all)`.
#' @param ppi `data.frame(uniprot_a, uniprot_b)` among the targets.
#' @return `list(drawn = ordered character, n_no_partner = integer,
#'   n_replaced = integer)`.
#' @keywords internal
.ppi_select_targets <- function(targets, ppi, max_targets) {
  ord <- order(-targets$is_disease, -targets$n_compounds, -targets$ppi_degree_all, targets$uniprot_id)
  cand <- targets$uniprot_id[ord][targets$ppi_degree_all[ord] > 0]
  n_no_partner <- sum(targets$ppi_degree_all == 0)
  n_take <- min(length(cand), max_targets)
  selected <- cand[seq_len(n_take)]
  pool <- cand[-seq_len(n_take)]
  rejected <- character(0)
  repeat {
    sub <- ppi[ppi$uniprot_a %in% selected & ppi$uniprot_b %in% selected, , drop = FALSE]
    lonely <- setdiff(selected, c(sub$uniprot_a, sub$uniprot_b))
    if (length(lonely) == 0) break
    selected <- setdiff(selected, lonely)
    rejected <- c(rejected, lonely)
    if (length(pool) == 0) {
      ## no refill possible: only drop what is still isolated
      sub <- ppi[ppi$uniprot_a %in% selected & ppi$uniprot_b %in% selected, , drop = FALSE]
      selected <- intersect(selected, c(sub$uniprot_a, sub$uniprot_b))
      break
    }
    add <- pool[seq_len(min(length(pool), length(lonely)))]
    pool <- setdiff(pool, add)
    selected <- c(selected, add)
  }
  list(drawn = cand[cand %in% selected], n_no_partner = n_no_partner, n_replaced = length(unique(rejected)))
}

#' Disease genes that are not targets but interact with >= 2 drawn targets
#' @return `data.frame(uniprot_id, n_links, association_score)`, at most `k` rows.
#' @keywords internal
.ppi_disease_partners <- function(ppi_all, drawn, candidates, dg, k) {
  e <- ppi_all[(ppi_all$uniprot_a %in% drawn & ppi_all$uniprot_b %in% candidates) |
                 (ppi_all$uniprot_b %in% drawn & ppi_all$uniprot_a %in% candidates), , drop = FALSE]
  if (nrow(e) == 0) return(data.frame(uniprot_id = character(0), n_links = integer(0), association_score = numeric(0)))
  partner <- ifelse(e$uniprot_a %in% candidates, e$uniprot_a, e$uniprot_b)
  tab <- table(partner)
  out <- data.frame(uniprot_id = names(tab), n_links = as.integer(tab), stringsAsFactors = FALSE)
  score <- if ("association_score" %in% names(dg)) stats::setNames(dg$association_score, dg$uniprot_id) else NULL
  out$association_score <- if (is.null(score)) NA_real_ else as.numeric(score[out$uniprot_id])
  out <- out[out$n_links >= 2, , drop = FALSE]
  out <- out[order(-out$n_links, -ifelse(is.na(out$association_score), 0, out$association_score), out$uniprot_id), , drop = FALSE]
  utils::head(out, k)
}

#' Force-directed layout of the drawn PPI graph, scaled into the unit disc
#'
#' @description
#' Each connected component is laid out on its own and the components are
#' packed with [igraph::merge_coords()], so small components do not drift to
#' the edge of the figure. The result is centred, its short axis stretched
#' by the square root of the aspect ratio, scaled so that the farthest node
#' lies at radius 1, radially stretched (`r^0.6`) to open up the dense core
#' that force-directed layouts give hub-rich PPI graphs, and decluttered
#' ([.ppi_declutter()]). With
#' `layout = "kk"` edge lengths follow `1000 / STRING score`.
#' @return Numeric matrix, one row per `ids`, columns x/y.
#' @keywords internal
.ppi_layout <- function(ids, ppi, layout = "fr", min_dist = 0.15) {
  g <- igraph::graph_from_data_frame(
    data.frame(from = ppi$uniprot_a, to = ppi$uniprot_b, weight = ppi$score / 1000, stringsAsFactors = FALSE),
    vertices = data.frame(name = ids, stringsAsFactors = FALSE), directed = FALSE
  )
  lay_one <- function(h) {
    if (igraph::vcount(h) == 1) return(matrix(0, 1, 2))
    if (igraph::vcount(h) == 2) return(matrix(c(-0.5, 0.5, 0, 0), 2, 2))
    if (layout == "kk") {
      igraph::layout_with_kk(h, weights = 1 / igraph::E(h)$weight)
    } else {
      igraph::layout_with_fr(h, niter = 3000)
    }
  }
  comps <- igraph::decompose(g)
  if (length(comps) == 1) {
    xy <- lay_one(g)
  } else {
    lays <- lapply(comps, lay_one)
    merged <- igraph::merge_coords(comps, lays)
    comp_names <- unlist(lapply(comps, function(h) igraph::V(h)$name))
    xy <- merged[match(ids, comp_names), , drop = FALSE]
  }
  xy <- sweep(xy, 2, (apply(xy, 2, max) + apply(xy, 2, min)) / 2)
  ## an elongated layout leaves half of the disc empty: the short axis is
  ## stretched by the square root of the aspect ratio (a compromise between
  ## filling the disc and keeping the layout's shape)
  ext <- apply(abs(xy), 2, max)
  if (all(ext > 0)) {
    short <- which.min(ext)
    xy[, short] <- xy[, short] * sqrt(max(ext) / min(ext))
  }
  r <- sqrt(rowSums(xy^2))
  if (max(r) > 0) {
    xy <- xy / max(r)
    ## force-directed layouts pack hubs into a dense core; a radial power
    ## transform (r -> r^0.6) opens the core up without changing the angles
    r <- r / max(r)
    stretch <- ifelse(r > 0, r^0.6 / r, 1)
    xy <- xy * stretch
    xy <- .ppi_declutter(xy, min_dist = min_dist)
  }
  xy
}

#' Push apart nodes closer than `min_dist`, inside the unit disc
#'
#' @description
#' Rounds of pairwise collision resolution (each overlapping pair is moved
#' apart along the line joining it by half the missing distance); a node
#' pushed outside the unit disc is put back on its rim. Keeps the layout's
#' shape but stops hub clusters from drawing nodes on top of each other.
#' @return The adjusted coordinate matrix.
#' @keywords internal
.ppi_declutter <- function(xy, min_dist = 0.15, iter = 400) {
  n <- nrow(xy)
  if (n < 2 || min_dist <= 0) return(xy)
  for (k in seq_len(iter)) {
    dx <- outer(xy[, 1], xy[, 1], "-")
    dy <- outer(xy[, 2], xy[, 2], "-")
    d <- sqrt(dx^2 + dy^2)
    diag(d) <- Inf
    close <- d < min_dist
    if (!any(close)) break
    ## coincident points: separate along a deterministic direction
    same <- which(close & d < 1e-9, arr.ind = TRUE)
    if (nrow(same) > 0) {
      i <- same[same[, 1] < same[, 2], 1]
      xy[i, 1] <- xy[i, 1] + min_dist / 4
      next
    }
    push <- ifelse(close, (min_dist - d) / (2 * d), 0)
    xy[, 1] <- xy[, 1] + rowSums(push * dx)
    xy[, 2] <- xy[, 2] + rowSums(push * dy)
    r <- sqrt(rowSums(xy^2))
    out <- r > 1
    xy[out, ] <- xy[out, , drop = FALSE] / r[out]
  }
  xy
}

#' Place compounds evenly on a ring, each near the side of its targets
#'
#' @description
#' A compound's preferred angle is the direction of the centroid of its drawn
#' targets. Compounds are sorted by that angle and spaced evenly around the
#' ring, rotated to best match the preferred angles (circular mean of the
#' offsets), so labels never collide and lines stay short.
#' @return `data.frame(compound_id, n_drawn, angle, x, y)`.
#' @keywords internal
.ppi_compound_ring <- function(ce, nodes, radius = 1.2) {
  if (nrow(ce) == 0) {
    return(data.frame(compound_id = character(0), n_drawn = integer(0), angle = numeric(0),
                      x = numeric(0), y = numeric(0), stringsAsFactors = FALSE))
  }
  pos <- match(ce$uniprot_id, nodes$uniprot_id)
  cx <- tapply(nodes$x[pos], ce$compound_id, mean)
  cy <- tapply(nodes$y[pos], ce$compound_id, mean)
  n_drawn <- tapply(pos, ce$compound_id, length)
  pref <- atan2(cy, cx)
  ord <- order(pref, names(pref))
  m <- length(pref)
  base <- 2 * pi * (seq_len(m) - 1) / m
  d <- pref[ord] - base
  offset <- atan2(mean(sin(d)), mean(cos(d)))
  angle <- base + offset
  data.frame(compound_id = names(pref)[ord], n_drawn = as.integer(n_drawn[ord]), angle = angle,
             x = radius * cos(angle), y = radius * sin(angle), stringsAsFactors = FALSE)
}

#' Which drawn targets get a text label
#'
#' @description
#' Up to `n` targets. Disease genes (in selection order) take at most 60% of
#' the labels and the other targets -- ranked by PPI degree in the drawn
#' graph, then number of compounds -- the rest, so the hubs that are not
#' disease genes are named too; a group with fewer members than its share
#' passes the remainder to the other group. Disease partners are not counted.
#' Highlighted (GPCR) targets -- up to 15, most compounds first -- are
#' labelled in addition to the `n`.
#' @return Integer row indices into `nodes`.
#' @keywords internal
.ppi_label_pick <- function(nodes, n) {
  if (n <= 0) return(integer(0))
  tgt <- which(nodes$kind == "target")
  dis <- tgt[nodes$is_disease[tgt]]
  dis <- dis[order(nodes$rank[dis])]
  oth <- tgt[!nodes$is_disease[tgt]]
  oth <- oth[order(-nodes$degree[oth], -nodes$n_compounds[oth], nodes$rank[oth])]
  n_dis <- min(length(dis), if (length(oth) == 0) n else ceiling(0.6 * n))
  n_oth <- min(length(oth), n - n_dis)
  n_dis <- min(length(dis), n - n_oth)
  picked <- c(utils::head(dis, n_dis), utils::head(oth, n_oth))
  ## highlighted (GPCR) targets are the point of the red overlay: up to 15
  ## of them are labelled on top of the n
  if ("is_gpcr" %in% names(nodes)) {
    hl <- tgt[nodes$is_gpcr[tgt]]
    hl <- hl[order(-nodes$n_compounds[hl], -nodes$degree[hl], nodes$rank[hl])]
    picked <- unique(c(picked, utils::head(hl, 15)))
  }
  picked
}

#' Colours of [plot_ppi_network()]
#'
#' @description
#' Okabe-Ito based. Red is reserved for the GPCR / highlighted proteins and
#' their interactions; the disease overlay uses reddish purple (fill, ring,
#' disease-disease edges), so the two overlays stay distinguishable when
#' they are drawn together.
#' @keywords internal
.ppi_colours <- function() {
  c(compound = "#009E73", compound_text = "#00553D", target = "#56B4E9",
    disease = "#CC79A7", disease_ring = "#6A1B5A", disease_edge = "#8E4585", disease_text = "#5E1F52",
    gpcr = "#D7191C", gpcr_text = "#A50F15", edge = "grey30")
}

#' Build the ggplot for [plot_ppi_network()]
#' @keywords internal
.ppi_ggplot <- function(nodes, ppi, cmp, ce, proj, conditions, disease, disease_name, colour_by, size_by,
                        top_n_labels, disease_node, score_threshold, engine, summary, scope_label,
                        n_cmp_hidden, fig_width, hl_name = "GPCR") {
  pal <- .ppi_colours()
  has_disease <- !is.null(disease)
  has_hl <- hl_name != "none"
  dlabel <- if (is.na(disease_name)) disease else disease_name
  pos <- stats::setNames(seq_len(nrow(nodes)), nodes$uniprot_id)
  is_target <- nodes$kind == "target"
  has_cmp <- !is.null(cmp) && nrow(cmp) > 0

  ## sizes: computed here (not by a scale) so ggrepel gets the same point
  ## size to keep labels off the nodes
  size_val <- if (size_by == "compounds") nodes$n_compounds else nodes$degree
  size_val[!is_target] <- NA
  size_max <- max(c(size_val, 1), na.rm = TRUE)
  size_min <- 1
  s_lo <- 1.8
  s_hi <- 5.6
  scale_fun <- function(v) s_lo + (s_hi - s_lo) * sqrt(pmax(v - size_min, 0) / max(size_max - size_min, 1))
  nodes$pt_size <- ifelse(is_target, scale_fun(size_val), 2.6)
  size_title <- if (size_by == "compounds") "Compounds hitting" else "PPI partners drawn"
  size_breaks <- unique(round(pretty(c(size_min, size_max), n = 3)))
  size_breaks <- size_breaks[size_breaks >= size_min & size_breaks <= size_max]
  if (length(size_breaks) == 0) size_breaks <- size_max

  ## fill groups
  target_label <- "Other predicted target"
  disease_label <- if (has_disease) paste0("Disease gene (", .plot_truncate(dlabel, 34), ")") else NULL
  if (colour_by == "module") {
    info <- .network_layers_module_colour_info(proj, conditions, nodes$uniprot_id)
    nodes$fill_group <- info$group
    fill_values <- info$palette[names(info$palette) %in% unique(info$group)]
    fill_name <- "Module"
  } else {
    nodes$fill_group <- if (has_disease) ifelse(nodes$is_disease, disease_label, target_label) else "Predicted target"
    fill_values <- if (has_disease) {
      stats::setNames(c(pal[["disease"]], pal[["target"]]), c(disease_label, target_label))
    } else {
      c(`Predicted target` = pal[["target"]])
    }
    fill_name <- "Protein"
  }
  hl_level <- paste0(hl_name, " (red outline)")
  nodes$shape_group <- ifelse(nodes$is_gpcr, hl_level, "Other protein")
  nodes$shape_code <- ifelse(nodes$is_gpcr, 22, 21)

  nodes$tooltip <- sprintf(
    "%s (%s)\n%s%s\nhit by %d compound(s); %d PPI partner(s) drawn",
    nodes$label, nodes$uniprot_id,
    ifelse(nodes$kind == "disease_partner", "disease gene, not a predicted target",
           ifelse(nodes$is_disease, "predicted target, disease gene", "predicted target")),
    ifelse(nodes$is_gpcr, paste0("; ", hl_name), ""),
    nodes$n_compounds, nodes$degree
  )
  nodes$name <- nodes$uniprot_id

  ## PPI edges and their class
  ppi$x <- nodes$x[pos[ppi$uniprot_a]]
  ppi$y <- nodes$y[pos[ppi$uniprot_a]]
  ppi$xend <- nodes$x[pos[ppi$uniprot_b]]
  ppi$yend <- nodes$y[pos[ppi$uniprot_b]]
  edge_levels <- c(paste0("Involves a ", hl_name), "Between two disease genes", "Other interaction")
  disease_pair <- nodes$is_disease[pos[ppi$uniprot_a]] & nodes$is_disease[pos[ppi$uniprot_b]]
  ppi$edge_class <- ifelse(ppi$gpcr_edge, edge_levels[1], ifelse(disease_pair, edge_levels[2], edge_levels[3]))
  ppi <- ppi[order(ppi$gpcr_edge, ppi$score), , drop = FALSE]
  edge_values <- stats::setNames(c(pal[["gpcr"]], pal[["disease_edge"]], pal[["edge"]]), edge_levels)
  edge_present <- edge_levels[edge_levels %in% ppi$edge_class]

  p <- ggplot2::ggplot()

  ## compound -> target lines (under everything)
  if (has_cmp) {
    cpos <- match(ce$compound_id, cmp$compound_id)
    tpos <- pos[ce$uniprot_id]
    ce_seg <- data.frame(x = cmp$x[cpos], y = cmp$y[cpos], xend = nodes$x[tpos], yend = nodes$y[tpos])
    p <- p + ggplot2::geom_segment(
      data = ce_seg, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      colour = pal[["compound"]], alpha = 0.07, linewidth = 0.2
    )
  }

  ## disease node and its dashed links
  dnode <- NULL
  if (disease_node) {
    dg_nodes <- nodes[nodes$is_disease, , drop = FALSE]
    if (nrow(dg_nodes) > 0) {
      yr <- if (has_cmp) 1.55 else 1.3
      dnode <- data.frame(x = 0, y = -yr, label = .plot_wrap(dlabel, 30), stringsAsFactors = FALSE)
      dseg <- data.frame(x = dg_nodes$x, y = dg_nodes$y, xend = dnode$x, yend = dnode$y)
      p <- p + ggplot2::geom_curve(
        data = dseg, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        colour = pal[["disease"]], alpha = 0.45, linewidth = 0.3, linetype = "22", curvature = 0.12
      )
    }
  }

  ## PPI edges: opacity by STRING score, colour by class, GPCR edges thicker
  base <- ppi[!ppi$gpcr_edge, , drop = FALSE]
  red <- ppi[ppi$gpcr_edge, , drop = FALSE]
  edge_aes <- ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                           alpha = .data$score, colour = .data$edge_class)
  if (nrow(base) > 0) p <- p + ggplot2::geom_segment(data = base, mapping = edge_aes, linewidth = 0.55)
  if (nrow(red) > 0) p <- p + ggplot2::geom_segment(data = red, mapping = edge_aes, linewidth = 1.05)
  p <- p +
    ggplot2::scale_colour_manual(
      values = edge_values, breaks = edge_present, name = "Interaction",
      guide = ggplot2::guide_legend(override.aes = list(linewidth = 1.1, alpha = 1), order = 3, ncol = 1)
    ) +
    ggplot2::scale_alpha_continuous(
      range = c(0.3, 0.9), limits = c(score_threshold, 1000), name = "STRING score",
      guide = ggplot2::guide_legend(override.aes = list(linewidth = 1, colour = "grey30"), order = 4, ncol = 1)
    )

  ## disease ring (under the fill layer, larger stroke -> visible border)
  ring <- nodes[nodes$is_disease & is_target, , drop = FALSE]
  if (has_disease && nrow(ring) > 0) {
    p <- p + ggplot2::geom_point(
      data = ring, ggplot2::aes(x = .data$x, y = .data$y, size = .data$pt_size + 1.7),
      shape = ring$shape_code, fill = NA, colour = pal[["disease_ring"]], stroke = 1.1, show.legend = FALSE
    )
  }
  ## disease partners: hollow purple symbols
  partners <- nodes[!is_target, , drop = FALSE]
  if (nrow(partners) > 0) {
    p <- p + .network_view_points(engine, partners, list(size = "pt_size"),
                                  shape = partners$shape_code, fill = "white",
                                  colour = ifelse(partners$is_gpcr, pal[["gpcr"]], pal[["disease_ring"]]),
                                  stroke = 1.1, show.legend = FALSE)
  }
  tg <- nodes[is_target, , drop = FALSE]
  p <- p + .network_view_points(engine, tg, list(size = "pt_size", fill = "fill_group", shape = "shape_group"),
                                colour = ifelse(tg$is_gpcr, pal[["gpcr"]], "grey20"),
                                stroke = ifelse(tg$is_gpcr, 1.1, 0.3))

  ## compounds: diamonds on the ring, labels outside it
  if (has_cmp) {
    cmp$name <- cmp$compound_id
    cmp$tooltip <- sprintf("%s\n%d predicted targets, %d drawn", cmp$label, cmp$n_targets, cmp$n_drawn)
    p <- p + .network_view_points(engine, cmp, list(), shape = 23, size = 4.2,
                                  fill = pal[["compound"]], colour = "grey15", stroke = 0.4)
    ca <- cos(cmp$angle)
    sa <- sin(cmp$angle)
    cmp$lx <- cmp$x + 0.075 * ca
    cmp$ly <- cmp$y + 0.075 * sa
    ## labels run outward (left half right-aligned, right half left-aligned);
    ## only a compound exactly at the top/bottom is centred, and those near
    ## the top/bottom sit above/below their diamond, so two neighbours on
    ## either side of the bottom never write over each other
    cmp$hjust <- ifelse(abs(ca) < 0.05, 0.5, ifelse(ca > 0, 0, 1))
    cmp$vjust <- ifelse(abs(sa) > 0.9, ifelse(sa > 0, 0, 1), 0.5)
    cmp$short <- .plot_truncate(cmp$label, 30)
    p <- p + ggplot2::geom_text(
      data = cmp, ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$short,
                               hjust = .data$hjust, vjust = .data$vjust),
      size = 2.7, colour = pal[["compound_text"]], lineheight = 0.9
    )
  }

  ## target labels: every node is passed so the repulsion avoids unlabelled
  ## nodes too; only the picked ones get text
  lab <- nodes
  keep <- rep(FALSE, nrow(nodes))
  keep[.ppi_label_pick(nodes, top_n_labels)] <- TRUE
  keep[!is_target] <- TRUE
  lab$text <- ifelse(keep, lab$label, "")
  lab$face <- ifelse(lab$is_disease, "bold", "plain")
  lab$face[!is_target] <- "bold.italic"
  lab$col <- ifelse(lab$is_gpcr, pal[["gpcr_text"]], ifelse(lab$is_disease, pal[["disease_text"]], "grey10"))
  lab_map <- if (.plot_use_repel()) {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$text, fontface = .data$face, point.size = .data$pt_size)
  } else {
    ggplot2::aes(x = .data$x, y = .data$y, label = .data$text, fontface = .data$face)
  }
  p <- p + .plot_text_layer(
    data = lab, mapping = lab_map,
    colour = lab$col, size = 2.7, show.legend = FALSE,
    repel_args = list(box.padding = 0.3, point.padding = 0.15, min.segment.length = 0.25,
                      segment.colour = "grey45", segment.size = 0.25, max.overlaps = Inf,
                      bg.colour = "white", bg.r = 0.12, seed = 1, force = 2, force_pull = 0.8,
                      max.time = 5, max.iter = 50000),
    text_args = list(vjust = -1.1)
  )

  if (!is.null(dnode)) {
    p <- p +
      ggplot2::geom_point(data = dnode, ggplot2::aes(x = .data$x, y = .data$y),
                          shape = 21, size = 8, fill = pal[["disease"]], colour = pal[["disease_ring"]], stroke = 1.2) +
      ggplot2::geom_text(data = dnode, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                         vjust = 2.3, size = 3.2, fontface = "bold", colour = pal[["disease_text"]], lineheight = 0.9)
  }

  ## legends (one column each, side by side under the network)
  shape_values <- stats::setNames(c(22, 21), c(hl_level, "Other protein"))
  shape_present <- names(shape_values)[names(shape_values) %in% nodes$shape_group[is_target]]
  p <- p +
    ggplot2::scale_fill_manual(
      values = fill_values, breaks = names(fill_values), name = fill_name,
      guide = ggplot2::guide_legend(override.aes = list(size = 4, shape = 21, colour = "grey20", stroke = 0.3),
                                    order = 1, ncol = if (length(fill_values) > 6) 2 else 1)
    ) +
    ggplot2::scale_shape_manual(
      values = shape_values, breaks = shape_present, name = "Receptor class",
      guide = if (has_hl && any(nodes$is_gpcr[is_target])) {
        ggplot2::guide_legend(override.aes = list(size = 4, fill = "white", colour = c(pal[["gpcr"]], "grey20")[match(shape_present, names(shape_values))],
                                                  stroke = c(1.1, 0.3)[match(shape_present, names(shape_values))]),
                              order = 2, ncol = 1)
      } else "none"
    ) +
    ggplot2::scale_size_identity(guide = "legend", name = size_title,
                                 breaks = scale_fun(size_breaks), labels = size_breaks) +
    ggplot2::guides(size = ggplot2::guide_legend(
      title = size_title, order = 5, ncol = 1,
      override.aes = list(shape = 21, fill = "grey80", colour = "grey30", stroke = 0.3)
    ))

  ## frame: room for the compound labels on both sides
  lim_x <- if (has_cmp) 1.95 else 1.15
  lim_y_hi <- if (has_cmp) 1.35 else 1.1
  lim_y_lo <- if (!is.null(dnode)) -(abs(dnode$y) + 0.28) else -lim_y_hi
  p <- p +
    ggplot2::coord_equal(xlim = c(-lim_x, lim_x), ylim = c(lim_y_lo, lim_y_hi), clip = "off", expand = FALSE)

  ## titles
  s <- summary
  title <- paste0("Protein-protein interactions among predicted targets -- ", scope_label)
  if (has_disease) title <- paste0(title, " | ", dlabel)
  sub1 <- sprintf("%d of %d predicted targets drawn, joined by %d STRING interactions (combined score >= %s).",
                  s$n_targets_drawn, s$n_targets, s$n_ppi_edges_drawn, format(score_threshold))
  omitted <- s$n_targets - s$n_targets_drawn
  sub2 <- sprintf(
    "Not drawn: %d target(s) with no STRING partner among the extract's targets%s%s.",
    s$n_no_partner,
    if (s$n_targets_mapped < s$n_targets) sprintf(" (%d not in STRING)", s$n_targets - s$n_targets_mapped) else "",
    if (omitted - s$n_no_partner > 0) sprintf("; %d lower-ranked target(s) beyond max_targets", omitted - s$n_no_partner) else ""
  )
  sub3 <- if (has_disease) {
    sprintf(
      "Disease genes: %d of %d (STRING-mapped) are predicted targets, %d drawn; fold enrichment %s, hypergeometric p = %s (universe: %s STRING proteins).",
      s$n_disease_overlap_string, s$n_disease_genes_mapped, s$n_disease_hit_drawn,
      if (is.na(s$fold_enrichment)) "NA" else formatC(s$fold_enrichment, format = "f", digits = 2),
      if (is.na(s$p_hypergeom)) "NA" else format.pval(s$p_hypergeom, digits = 2, eps = 1e-300),
      format(s$n_universe, big.mark = ",")
    )
  } else NULL
  sub4 <- if (has_hl) {
    sprintf("%ss: %d of %d drawn proteins (%d of %d predicted targets); %d interaction(s) involving one drawn in red.",
            hl_name, s$n_gpcr_drawn, nrow(nodes), s$n_gpcr_targets, s$n_targets, s$n_gpcr_edges_drawn)
  } else NULL
  subtitle <- .plot_wrap(paste(c(sub1, sub2, sub3, sub4), collapse = "\n"), .plot_wrap_width(fig_width, 8.5))
  cap <- c(
    sprintf("Targets ranked by %snumber of compounds hitting them and STRING degree; each drawn target has >= 1 STRING partner in the figure.",
            if (has_disease) "disease-gene status, " else ""),
    if (has_cmp) "Diamonds: compounds, with thin lines to the drawn targets they are predicted to hit." else NULL,
    if (n_cmp_hidden > 0) sprintf("%d compound(s) hit none of the drawn targets and are not shown.", n_cmp_hidden) else NULL,
    if (has_disease) "Purple ring: disease gene." else NULL,
    if (nrow(partners) > 0) "Hollow symbols: disease genes that are not predicted targets but interact with >= 2 drawn targets." else NULL,
    if (has_hl && hl_name %in% c("GPCR", "GPCR / highlighted")) "GPCR: GO:0004930 (G protein-coupled receptor activity) and descendants, org.Hs.eg.db." else NULL
  )
  p <- p +
    ggplot2::labs(title = .plot_wrap(title, .plot_wrap_width(fig_width, 12)), subtitle = subtitle,
                  caption = .plot_wrap(paste(cap, collapse = " "), .plot_wrap_width(fig_width, 7.5))) +
    .network_view_theme(legend_position = "bottom") +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 14, 8, 14),
                   legend.box = "horizontal", legend.title.position = "top",
                   legend.box.just = "top", legend.spacing.x = ggplot2::unit(0.5, "cm"))
  p
}
