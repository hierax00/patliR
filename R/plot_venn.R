#' @include AllGenerics.R internal.R network_build.R targets_disease.R
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
#'
#' @return A `ggplot` object (`ggVennDiagram` builds on `ggplot2`, no
#'   separate `engine` argument). If `save = TRUE` (default), also writes
#'   a PNG and logs it to `patliRResults(proj, "venn_plot_log")`.
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
                       save = TRUE, out_dir = NULL, width = 6, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
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
  compound_targets <- unique(edges_all$uniprot_id[edges_all$condition %in% conditions])

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

  venn_list <- list(compound_targets = compound_targets, disease_targets = disease_targets)
  p <- ggVennDiagram::ggVennDiagram(venn_list, label_alpha = 0) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "#2980b9", name = "Count") +
    ggplot2::labs(title = paste0("Compound targets vs. ", disease, " targets -- ", scope_label))

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("venn_", scope_label, "_", disease, ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
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
  p
}
