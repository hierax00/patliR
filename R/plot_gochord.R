#' @include AllGenerics.R internal.R network_enrich.R
NULL

## Uses GOplot::GOChord() for the ribbon geometry. Called with `nlfc = 0`
## and no `logFC` column: target prediction has no expression fold-change,
## and a synthetic one would render a colored ring that looks like real
## per-gene magnitude but is not. GOChord()'s default `nlfc = 1` errors
## when no `logFC` column is present, so `0` must be passed explicitly.
## GOChord() returns a real ggplot object, so the `engine = c("static",
## "ggiraph")` contract is kept (girafe wraps the static picture; no
## per-ribbon tooltips).

#' Circular gene-pathway ribbon chord diagram (real `GOplot::GOChord()`)
#'
#' @description
#' The `top_n_terms` most significant [network_enrich()] pathway/GO terms
#' for `condition`, plus every gene (Entrez -> symbol, via the same
#' `org.Hs.eg.db` mapping [network_enrich()] itself uses) belonging to at
#' least one of them, rendered as a true ribbon chord diagram via
#' `GOplot::GOChord()` -- curved ribbons from each gene to every term it
#' belongs to, colored by term.
#'
#' @section No `logFC` ring:
#' `GOChord()` optionally colors an inner ring by a `logFC` column. This
#' package has no real per-gene fold-change to put there (target
#' prediction, not differential expression), so this function always calls
#' `GOChord(..., nlfc = 0)`: ribbons only, no inner ring, no `logFC` legend.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param db Passed through to filter `network_enrichment$db` if more than
#'   one enrichment database was run for `condition` (e.g. `"go"` vs.
#'   `"reactome"`); `NULL` (default) uses whichever is present, and errors
#'   if both are and `db` was not specified.
#' @param top_n_terms Integer, default `10`. Only the `top_n_terms` most
#'   significant terms (lowest `p.adjust`) are drawn -- with dozens of
#'   significant terms the ribbon diagram becomes unreadable.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "gochord_plot_log")`.
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
#' proj <- network_enrich(proj, condition = "FLO-ET", db = "go")
#' plot_gochord(proj, condition = "FLO-ET", top_n_terms = 8, save = FALSE)
#' }
#'
#' @export
plot_gochord <- function(proj, condition = NULL, db = NULL, top_n_terms = 10,
                          engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                          width = 14, height = 14, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n_terms), length(top_n_terms) == 1, top_n_terms >= 1)
  engine <- match.arg(engine)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_gochord}.")
  }
  if (!requireNamespace("GOplot", quietly = TRUE)) {
    cli::cli_abort("The {.pkg GOplot} package is required for {.fn plot_gochord}.")
  }
  if (engine == "ggiraph" && !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  conditions <- .network_resolve_conditions(proj, condition)
  if (length(conditions) != 1L) {
    cli::cli_abort("{.fn plot_gochord} needs a single {.arg condition} (enrichment terms are per-condition, pooling them is not meaningful).")
  }
  cond <- conditions[[1]]

  enrichment_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enrichment_all) || nrow(enrichment_all) == 0) {
    cli::cli_abort(c("No {.val network_enrichment} entry in {.arg proj}.", "i" = "Run {.fn network_enrich} first."))
  }
  enr <- enrichment_all[enrichment_all$condition == cond, , drop = FALSE]
  if (!is.null(db)) enr <- enr[enr$db == db, , drop = FALSE]
  avail_db <- unique(enr$db)
  if (length(avail_db) > 1) {
    cli::cli_abort("More than one {.val db} present for condition {.val {cond}} ({.val {avail_db}}); pass {.arg db} explicitly.")
  }
  if (nrow(enr) == 0) {
    cli::cli_abort("No {.val network_enrichment} rows for condition {.val {cond}}{if (!is.null(db)) paste0(' / db ', db) else ''}.")
  }
  enr <- enr[order(enr$p.adjust), , drop = FALSE]
  enr <- utils::head(enr, top_n_terms)

  membership <- do.call(rbind, lapply(seq_len(nrow(enr)), function(i) {
    entrez_ids <- strsplit(enr$geneID[i], "/")[[1]]
    data.frame(term_label = enr$Description[i], entrez_id = entrez_ids, stringsAsFactors = FALSE)
  }))
  gene_label <- .network_gochord_gene_labels(unique(membership$entrez_id))
  membership$gene_label <- gene_label[membership$entrez_id]
  membership$gene_label[is.na(membership$gene_label)] <- membership$entrez_id[is.na(membership$gene_label)]

  genes <- sort(unique(membership$gene_label))
  terms <- unique(membership$term_label) # already ordered by p.adjust from `enr`
  chord_matrix <- matrix(0L, nrow = length(genes), ncol = length(terms), dimnames = list(genes, terms))
  chord_matrix[cbind(membership$gene_label, membership$term_label)] <- 1L

  p <- GOplot::GOChord(
    chord_matrix,
    title = paste0("Gene-term membership -- ", cond, " (", enr$db[1], ")"),
    nlfc = 0
  )

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, paste0("gochord_", cond, "_", enr$db[1], ".png"))
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
    log_row <- data.frame(condition = cond, db = enr$db[1], path = path, n_terms = nrow(enr), stringsAsFactors = FALSE)
    log_df <- .network_upsert(proj, "gochord_plot_log", log_row, c("condition", "db"))
    patliRResults(proj, "gochord_plot_log") <- log_df
    .write_results_csv(proj, "gochord_plot_log", log_df)
  }

  result <- if (engine == "static") p else ggiraph::girafe(ggobj = p)
  if (save) attr(result, "proj") <- proj
  result
}

#' Entrez -> gene symbol lookup for `plot_gochord()`
#' @return Named character vector (`names` = Entrez IDs), `NA` for
#'   unmapped IDs or if `org.Hs.eg.db` is not installed.
#' @keywords internal
.network_gochord_gene_labels <- function(entrez_ids) {
  out <- stats::setNames(rep(NA_character_, length(entrez_ids)), entrez_ids)
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    return(out)
  }
  map <- tryCatch(
    clusterProfiler::bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = TRUE),
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) return(out)
  out[map$ENTREZID] <- map$SYMBOL
  out
}
