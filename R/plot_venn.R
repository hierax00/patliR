#' @include AllGenerics.R internal.R network_build.R targets_disease.R plot-helpers.R
NULL

## Uses ggVennDiagram (Suggests): proportional circle-intersection
## geometry is not worth hand-rolling when a small ggplot2-native package
## does it well. Same reasoning as ggalluvial in plot_bowtie.R.

#' Venn diagram: compound targets vs. disease-associated targets
#'
#' @description
#' Two-set overlap: `"compound_targets"` (every UniProt ID any compound
#' hits in `condition`'s `network_edges`) against `"disease_targets"`
#' (every UniProt ID associated with `disease` in `targets_disease`, from
#' [targets_disease_filter()]). The overlap is literally "this extract's
#' targets that are also known disease genes for this indication" -- a
#' quick sanity check before running the heavier, statistically rigorous
#' version of the same question, [network_proximity()].
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param disease Character EFO ID from `targets_disease$disease_id`. If
#'   `NULL` (default) and exactly one disease is present in
#'   `patliRResults(proj, "targets_disease")`, that one is used; otherwise
#'   required.
#' @param sets Character vector of exactly the two sets to compare,
#'   currently only `c("compound_targets", "disease_targets")` (default,
#'   and the only supported value) -- kept as an explicit argument, not
#'   hardcoded, so a third set (e.g. STRING bowtie core membership) can be
#'   added later without changing the function's contract.
#' @param label_top Integer, default `15`. Names the `label_top` shared
#'   targets (in both sets) hit by the most compounds in scope, as a
#'   ranked list beside the diagram with each target's compound count
#'   (gene symbols when `clusterProfiler` + `org.Hs.eg.db` are installed,
#'   UniProt IDs otherwise). `0` draws the bare diagram. The list needs the
#'   `patchwork` package; without it the names go into the caption. When
#'   `width` is not given, a labeled figure is saved 9 inches wide instead
#'   of 6.
#'
#' @return A `ggplot` object (`ggVennDiagram` builds on `ggplot2`, no
#'   separate `engine` argument; a `patchwork` -- also a `ggplot` -- when
#'   the shared-target list is drawn). `attr(., "shared_targets")` holds
#'   every shared target (`target_id`, `label`, `n_compounds`), ranked.
#'   The disease name (from [disease_genes_fetch()] /
#'   [targets_disease_profile()]) is shown next to its ID when available. If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "venn_plot_log")`.
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
#' plot_venn(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @export
plot_venn <- function(proj, condition = NULL, disease = NULL,
                       sets = c("compound_targets", "disease_targets"),
                       label_top = 15,
                       save = TRUE, out_dir = NULL, width = 6, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  if (!(is.numeric(label_top) && length(label_top) == 1 && !is.na(label_top) && label_top >= 0)) {
    cli::cli_abort("{.arg label_top} must be a single number >= 0.")
  }
  width_given <- !missing(width)
  if (!identical(sets, c("compound_targets", "disease_targets"))) {
    cli::cli_abort("{.arg sets} currently only supports {.val {c('compound_targets', 'disease_targets')}}.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_venn}.")
  }
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_venn} needs the {.pkg ggVennDiagram} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"ggVennDiagram\")} -- pure CRAN, ggplot2-based."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")

  edges_all <- patliRResults(proj, "network_edges")
  edges_scope <- edges_all[edges_all$condition %in% conditions, , drop = FALSE]
  compound_targets <- unique(edges_scope$uniprot_id)

  disease_all <- patliRResults(proj, "targets_disease")
  if (is.null(disease_all) || nrow(disease_all) == 0) {
    cli::cli_abort(c("No {.val targets_disease} entry in {.arg proj}.", "i" = "Run {.fn targets_disease_filter} first."))
  }
  if (is.null(disease)) {
    avail <- unique(disease_all$disease_id)
    if (length(avail) != 1) {
      cli::cli_abort("{.arg disease} is required when {.val targets_disease} has more than one disease_id; available: {.val {avail}}.")
    }
    disease <- avail
  }
  disease_targets <- unique(disease_all$target_id[disease_all$disease_id == disease])
  if (length(disease_targets) == 0) {
    cli::cli_abort("No target associated with {.val {disease}} in {.val targets_disease}.")
  }

  disease_label <- .plot_disease_label(proj, disease)
  disease_name <- .plot_disease_label(proj, disease, with_id = FALSE)

  venn_list <- list(compound_targets = compound_targets, disease_targets = disease_targets)
  ## readable set names (the disease name, wrapped), and room on both sides
  ## so the labels are not clipped
  names(venn_list) <- c("Compound targets", paste0("Disease targets\n(", .plot_wrap(disease_name, 26), ")"))
  title <- paste0("Compound targets vs. disease targets\n", .plot_wrap(disease_label, 60), " -- ", scope_label)
  p <- ggVennDiagram::ggVennDiagram(venn_list, label_alpha = 0) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 20, 10, 90)) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2980b9", name = "Count") +
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ## the overlap, named: shared targets ranked by how many compounds hit them
  shared <- .venn_shared_targets(edges_scope, intersect(compound_targets, disease_targets))
  shared$label <- if (nrow(shared) > 0) .plot_label_nodes(proj, conditions, shared$target_id, "target") else character(0)
  shared <- shared[, c("target_id", "label", "n_compounds"), drop = FALSE]
  if (label_top > 0 && nrow(shared) > 0) {
    top <- utils::head(shared, label_top)
    if (!width_given) width <- 9
    p <- .venn_add_shared_list(p, top, nrow(shared), title)
  }

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("venn_", scope_label, "_", disease, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
    n_overlap <- length(intersect(compound_targets, disease_targets))
    log_row <- data.frame(
      condition = scope_label, disease_id = disease, path = path,
      n_compound_targets = length(compound_targets), n_disease_targets = length(disease_targets),
      n_overlap = n_overlap, stringsAsFactors = FALSE
    )
    log_df <- .network_upsert(proj, "venn_plot_log", log_row, c("condition", "disease_id"))
    patliRResults(proj, "venn_plot_log") <- log_df
    .write_results_csv(proj, "venn_plot_log", log_df)
  }

  if (save) attr(p, "proj") <- proj
  attr(p, "shared_targets") <- shared
  p
}

#' Shared targets ranked by the number of compounds hitting them
#' @param edges `network_edges` rows in scope (`compound_id`, `uniprot_id`).
#' @param shared_ids UniProt IDs in both sets.
#' @return `data.frame(target_id, n_compounds)`, most-hit first (ties by ID).
#' @keywords internal
.venn_shared_targets <- function(edges, shared_ids) {
  if (length(shared_ids) == 0) {
    return(data.frame(target_id = character(0), n_compounds = integer(0), stringsAsFactors = FALSE))
  }
  e <- unique(edges[edges$uniprot_id %in% shared_ids, c("compound_id", "uniprot_id"), drop = FALSE])
  n <- table(factor(e$uniprot_id, levels = sort(unique(shared_ids))))
  out <- data.frame(target_id = names(n), n_compounds = as.integer(n), stringsAsFactors = FALSE)
  out <- out[order(-out$n_compounds, out$target_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Put the ranked shared-target list beside the Venn diagram
#' @description A `patchwork` of the diagram and a text column
#'   (`"1. AKT1 (12)"` ...); without `patchwork`, the list goes into the
#'   caption.
#' @keywords internal
.venn_add_shared_list <- function(p, top, n_shared, title) {
  lines <- sprintf("%d. %s (%d)", seq_len(nrow(top)), .plot_truncate(top$label, 22), top$n_compounds)
  heading <- paste0("Top ", nrow(top), " of ", n_shared, " shared targets\n(n compounds hitting each)")
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(p + ggplot2::labs(caption = .plot_wrap(paste0(sub("\n", " ", heading, fixed = TRUE), ": ",
                                                         paste(lines, collapse = ", ")), 90)))
  }
  n <- length(lines)
  side <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = n + 2.3, label = heading, hjust = 0, vjust = 1,
                      size = 3.3, fontface = "bold", colour = "grey15", lineheight = 0.95) +
    ggplot2::annotate("text", x = 0, y = rev(seq_len(n)), label = lines, hjust = 0,
                      size = 3.1, colour = "grey20") +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(min(0, n - 25), n + 2.5), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void()
  ## the title moves from the diagram to the whole figure, centred over both
  p <- p + ggplot2::labs(title = NULL)
  patchwork::wrap_plots(p, side, widths = c(2.3, 1)) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                             plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    )
}
