#' @include AllGenerics.R internal.R plot_ppi_network.R plot-helpers.R network_condition_compare.R
NULL

## Uses ggalluvial (Suggests), same rationale as plot_bowtie.R: curved flow
## ribbons between stacked strata are not reasonably hand-rollable with
## plain ggplot2 geoms.

#' Alluvial flow: compound -> GPCR -> G-protein family
#'
#' @description
#' Metabotropic (G protein-coupled) receptors are a large share of the
#' predicted targets for many plant extracts, and a single receptor's
#' downstream effect depends on *which* G protein it couples to -- this is
#' the piece [plot_ppi_network()] with `highlight = "gpcr"` does not show (it marks
#' GPCRs in the network, but not their coupling). This alluvial reads like
#' a small metro map: which compounds (left) reach which GPCRs (middle),
#' and which G-protein family (right, colour) each receptor signals
#' through.
#'
#' @section GPCR classification and coupling table (please read):
#' Receptors are the same GO:0004930-based set [plot_ppi_network()] uses
#' (see its "GPCRs" section for the electronic-annotation caveat -- a
#' handful of false positives are known, see below). G-protein coupling
#' comes from a small hand-curated table
#' (`inst/extdata/gpcr_gprotein_coupling.csv`), written from established
#' GPCR pharmacology (consistent with the IUPHAR/BPS Guide to
#' Pharmacology's *primary* coupling calls), **not** an automated lookup --
#' verify against GtoPdb before using this for a manuscript. It only
#' covers the receptors encountered while building it; a target classified
#' as a GPCR but absent from the table (a new receptor GO turns up, or one
#' of its 4 known false positives -- `HPGD`, `RORB`, `SPHK1`, `SPHK2`,
#' which are an enzyme, a nuclear receptor and two kinases mis-annotated
#' to GO:0004930 -- see the CSV's `note` column) is dropped from the flow
#' with a message, not silently mis-classified. Many GPCRs couple to more
#' than one G protein *in some tissue*; the table records one *primary*
#' family per receptor, not the full polypharmacology.
#'
#' @section Vascular effect (optional, `show_vascular_effect = TRUE`):
#' Adds a third flow into a coarse vasodilation/vasoconstriction/mixed/
#' positive inotropy-chronotropy call, but **only** for the receptors
#' where the literature actually documents a vascular or cardiac role
#' (mostly cardiovascular/renal receptors); every other receptor (most
#' chemokine, opioid, histamine, CNS ones) is shown as `"not established /
#' other tissue"` rather than guessing -- do not read that bucket as "no
#' effect", only as "not a mainstream vascular pharmacology fact".
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Optional disease id (as in [network_proximity()]/
#'   [targets_disease_filter()]). When given, a receptor that is *also* one
#'   of that disease's genes is marked with `*` after its label and
#'   reported in the subtitle -- it does not change which receptors are
#'   drawn.
#' @param top_n_receptors Integer, default `20`. Receptors are ranked by
#'   the number of compounds hitting them (ties broken by receptor name);
#'   with many predicted GPCRs the alluvial gets crowded.
#' @param show_vascular_effect Logical, default `FALSE`; see above.
#'
#' @return A `ggplot` object (`ggalluvial` has no `ggiraph` equivalent). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "gpcr_flow_plot_log")`. `attr(., "table")` holds
#'   every drawn `(compound_id, uniprot_id, gene_symbol, g_protein,
#'   vascular_effect, is_disease_gene)` combination.
#'
#' @seealso [plot_ppi_network()], [network_condition_compare()]
#' @export
plot_gpcr_flow <- function(proj, condition = NULL, disease = NULL, top_n_receptors = 20,
                           show_vascular_effect = FALSE,
                           save = TRUE, out_dir = NULL, width = 9, height = 7, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  .plot_require(extra = "ggalluvial")
  .pathway_check_count(top_n_receptors, "top_n_receptors", min = 1)
  .pathway_check_flag(show_vascular_effect, "show_vascular_effect")

  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c("No {.val network_edges} entry in {.arg proj}.", "i" = "Run {.fn network_build} first."))
  }
  ct <- unique(edges_all[edges_all$condition %in% conditions, c("compound_id", "uniprot_id"), drop = FALSE])

  gpcr <- .ppi_gpcr_go()
  if (is.null(gpcr)) {
    cli::cli_abort(c("{.fn plot_gpcr_flow} needs {.pkg org.Hs.eg.db} (Bioconductor) to classify GPCRs, not installed."))
  }
  ct <- ct[ct$uniprot_id %in% gpcr$uniprot_id, , drop = FALSE]
  if (nrow(ct) == 0) {
    cli::cli_abort("No predicted target for condition(s) {.val {conditions}} is classified as a GPCR (GO:0004930).")
  }

  coupling <- .gpcr_coupling_table()
  ct <- merge(ct, coupling, by = "uniprot_id", all.x = TRUE)
  n_unclassified <- sum(is.na(ct$g_protein))
  if (n_unclassified > 0) {
    cli::cli_inform(c("i" = paste0(
      n_unclassified, " GPCR-target row(s) dropped: not a real GPCR (a known GO:0004930 false positive) or ",
      "not yet in the curated coupling table (see {.fn plot_gpcr_flow}'s docs)."
    )))
  }
  ct <- ct[!is.na(ct$g_protein), , drop = FALSE]
  if (nrow(ct) == 0) {
    cli::cli_abort("No GPCR target in scope has a G-protein coupling call in the curated table.")
  }

  dz <- NULL
  if (!is.null(disease)) {
    dz <- tryCatch(.condition_compare_disease_genes(proj, disease), error = function(e) {
      cli::cli_warn(c(
        "!" = "Could not look up disease genes for {.val {disease}}: {conditionMessage(e)}",
        "i" = "No receptor will be marked as a disease gene -- this is \"not available\", not \"confirmed unrelated\"."
      ))
      NULL
    })
  }
  ct$is_disease_gene <- if (!is.null(dz)) ct$uniprot_id %in% dz$uniprot_id else FALSE

  ## ranked by number of compounds hitting the receptor, ties broken by gene symbol
  ## (not by accession, which has no biological meaning and would order arbitrarily)
  n_by_receptor <- stats::aggregate(compound_id ~ uniprot_id, ct, function(x) length(unique(x)))
  names(n_by_receptor)[2] <- "n_compounds"
  n_by_receptor$gene_symbol <- ct$gene_symbol[match(n_by_receptor$uniprot_id, ct$uniprot_id)]
  n_by_receptor <- n_by_receptor[order(-n_by_receptor$n_compounds, n_by_receptor$gene_symbol), , drop = FALSE]
  keep <- utils::head(n_by_receptor$uniprot_id, top_n_receptors)
  ct <- ct[ct$uniprot_id %in% keep, , drop = FALSE]

  ct$compound_label <- .plot_unique_labels(proj, conditions, ct$compound_id, "compound")
  ct$receptor_label <- paste0(ct$gene_symbol, ifelse(ct$is_disease_gene, "*", ""))
  if (show_vascular_effect) ct$vascular_effect[is.na(ct$vascular_effect)] <- "not established / other tissue"

  flow <- if (show_vascular_effect) {
    unique(ct[, c("compound_id", "compound_label", "receptor_label", "g_protein", "vascular_effect")])
  } else {
    unique(ct[, c("compound_id", "compound_label", "receptor_label", "g_protein")])
  }

  label_lookup <- stats::setNames(flow$compound_label, flow$compound_id)
  side_label <- function(stratum, x) {
    label <- unname(label_lookup[as.character(stratum)])
    ifelse(x == 1 & !is.na(label), .plot_truncate(label), ifelse(x == 1, as.character(stratum), ""))
  }
  left_labels <- .plot_truncate(unique(flow$compound_label))
  x_expand <- .bowtie_x_expansion(max(nchar(left_labels)), 6, width)
  compound_label_layer <- .plot_text_layer(
    stat = ggalluvial::StatStratum,
    mapping = ggplot2::aes(label = ggplot2::after_stat(side_label(.data$stratum, .data$x))),
    size = 2.8, hjust = 1, nudge_x = -(1 / 8 + 0.03),
    repel_args = list(direction = "y", min.segment.length = Inf, box.padding = 0.05,
                      point.padding = 0, force = 0.5, max.overlaps = Inf, seed = 1),
    text_args = list()
  )

  g_protein_colors <- c(
    "Gs" = "#e74c3c", "Gi/o" = "#2980b9", "Gq/11" = "#27ae60", "G12/13" = "#8e44ad"
  )

  axes <- if (show_vascular_effect) {
    ggplot2::aes(axis1 = .data$compound_id, axis2 = .data$receptor_label,
                axis3 = .data$g_protein, axis4 = .data$vascular_effect, y = 1)
  } else {
    ggplot2::aes(axis1 = .data$compound_id, axis2 = .data$receptor_label, axis3 = .data$g_protein, y = 1)
  }
  axis_names <- if (show_vascular_effect) c("compound", "receptor", "g_protein", "vascular_effect") else c("compound", "receptor", "g_protein")

  p <- ggplot2::ggplot(flow, axes) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = .data$g_protein), width = 1 / 6, alpha = 0.75) +
    ggalluvial::geom_stratum(width = 1 / 6, fill = "grey92", colour = "grey40") +
    compound_label_layer +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(ifelse(.data$x == 1, "", as.character(stratum)))),
                       size = 2.6) +
    ggplot2::scale_x_discrete(limits = axis_names, expand = ggplot2::expansion(add = x_expand)) +
    ggplot2::scale_fill_manual(values = g_protein_colors, name = "Primary\nG protein", na.translate = FALSE) +
    ggplot2::labs(
      title = paste0("Compound -> GPCR -> G protein -- ", scope_label),
      subtitle = .plot_wrap(paste0(
        "Flow = a predicted compound-receptor hit; colour = the receptor's primary G-protein family (curated, see docs). ",
        if (!is.null(disease)) sprintf("*: also a %s gene (%d of %d receptors shown). ", disease,
                                       sum(!duplicated(ct$uniprot_id) & ct$is_disease_gene), length(unique(ct$uniprot_id))) else "",
        "Top ", length(unique(flow$receptor_label)), " receptor(s) by number of compounds hitting them."
      ), .plot_wrap_width(width, 8)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_blank()
    )

  result <- .plot_finish(
    proj, p,
    name = "gpcr_flow_plot_log",
    filename = paste0("gpcr_flow_", scope_label, if (show_vascular_effect) "_vascular" else "", ".png"),
    log_row = data.frame(condition = scope_label, disease_id = disease %||% NA_character_,
                         n_receptors = length(unique(flow$receptor_label)), path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "disease_id"),
    engine = NULL, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- ct[, c("compound_id", "uniprot_id", "gene_symbol", "g_protein",
                                  intersect("vascular_effect", names(ct)), "is_disease_gene")]
  result
}

#' Read the curated GPCR -> primary G-protein (and, where established,
#' vascular effect) table bundled with the package
#' @return `data.frame(uniprot_id, gene_symbol, g_protein, vascular_effect)`;
#'   `g_protein`/`vascular_effect` are `NA` for the four known GO:0004930
#'   false positives (see the CSV's `note` column).
#' @keywords internal
.gpcr_coupling_table <- function() {
  path <- system.file("extdata", "gpcr_gprotein_coupling.csv", package = "patliR")
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d[, c("uniprot_id", "gene_symbol", "g_protein", "vascular_effect")]
}
