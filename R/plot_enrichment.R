#' @include AllGenerics.R internal.R network_enrich.R plot-helpers.R
NULL

#' Bubble plot of `network_enrich()` over-representation results
#'
#' @description
#' The standard `clusterProfiler::dotplot()`-style enrichment bubble: one
#' point per enriched term from `patliRResults(proj, "network_enrichment")`,
#' term name (`Description`) on `y`, `GeneRatio` (parsed from its `"k/n"`
#' string, e.g. `"3/20"` -> `0.15`) on `x`, point **size** = `Count` (genes
#' in the input set annotated to the term), point **colour** = `p.adjust`
#' (`scale_colour_gradient(low = "#c0392b", high = "#2980b9")` -- low
#' `p.adjust` = more significant = the "hot" red end, matching
#' [plot_proximity()]'s `TRUE`-is-red convention for "the thing worth
#' looking at"). Terms are ranked by `p.adjust` (most significant first)
#' to pick the `top_n` shown per `(condition, db)` group, then the `y`-axis
#' within that selection is ordered by `GeneRatio` -- the two orderings
#' answer different questions ("which terms are most confidently enriched"
#' vs. "which terms involve the largest fraction of the input set") and
#' `clusterProfiler::dotplot()` makes the same split. Faceted by `db` when
#' more than one database is in scope.
#'
#' @section `qvalue` may be absent:
#' `clusterProfiler` omits the `qvalue` column from its result whenever the
#' `qvalue` package's pi0 estimation fails, which is routine for the small
#' p-value vectors a 10-50-gene target set produces (see
#' [network_enrich()]'s own review notes). This plot does not use `qvalue`
#' at all, and does not require the column to be present.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param db Character vector of database names (`"go"`, `"reactome"`,
#'   `"kegg"`) to restrict to, or `NULL` (default) for every `db` present
#'   in scope.
#' @param top_n Integer, default `20`. Number of most-significant
#'   (lowest-`p.adjust`) terms to draw per `(condition, db)` group.
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object or `girafe` htmlwidget. If `save = TRUE`
#'   (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "enrichment_plot_log")`.
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
#' proj <- network_enrich(proj, db = "go") # needs org.Hs.eg.db + clusterProfiler
#' plot_enrichment(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @references
#' Yu, Wang, Han & He (2012), *OMICS* 16(5):284-287,
#' \doi{10.1089/omi.2011.0118} (the `clusterProfiler` dotplot form this
#' figure follows). Li S. (2021), *World J Tradit Chin Med* 7(1):146-154,
#' \doi{10.4103/wjtcm.wjtcm_11_21} (network-pharmacology reporting
#' guidance codifying the enrichment bubble as expected output).
#'
#' @export
plot_enrichment <- function(proj, condition = NULL, db = NULL, top_n = 20,
                             engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                             width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  enr_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enr_all) || nrow(enr_all) == 0) {
    cli::cli_abort(c("No {.val network_enrichment} entry in {.arg proj}.", "i" = "Run {.fn network_enrich} first."))
  }
  dat <- enr_all[enr_all$condition %in% conditions, , drop = FALSE]
  if (!is.null(db)) dat <- dat[dat$db %in% db, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_enrichment} rows for the requested condition(s)/db.")
  }
  db_label <- if (is.null(db)) "ALL" else paste(sort(unique(db)), collapse = "+")

  ## qvalue is not used by this plot and is not required to be present
  ## (see @section above); nothing to do here beyond documenting that we
  ## never touch it.

  dat$gene_ratio <- .enrichment_parse_ratio(dat$GeneRatio)
  unparsed <- is.na(dat$gene_ratio) & !is.na(dat$GeneRatio)
  if (any(unparsed)) {
    cli::cli_warn(c(
      "!" = "{sum(unparsed)} row(s) had a {.field GeneRatio} that could not be parsed as \"k/n\" and {?is/are} excluded: {.val {utils::head(unique(dat$GeneRatio[unparsed]), 5)}}."
    ))
  }
  dat <- dat[!is.na(dat$gene_ratio), , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No {.val network_enrichment} rows with a parseable {.field GeneRatio} for the requested condition(s)/db.")
  }

  multi_cond <- length(unique(dat$condition)) > 1
  multi_db <- length(unique(dat$db)) > 1
  dat$group <- paste(dat$condition, dat$db, sep = " / ")

  ## Select the top_n most-significant (lowest p.adjust, NA last) terms
  ## per (condition, db) group -- clusterProfiler::dotplot()'s own
  ## showCategory selection -- then order the y-axis within that
  ## selection by GeneRatio (a different question: "how large a fraction
  ## of the input set", not "how confidently enriched").
  dat <- do.call(rbind, lapply(split(dat, dat$group), function(d) {
    d <- d[order(d$p.adjust, na.last = TRUE), , drop = FALSE]
    utils::head(d, top_n)
  }))
  rownames(dat) <- NULL

  dat$term_label <- if (multi_cond) paste0(dat$condition, ": ", dat$Description) else dat$Description
  dat$term_label <- factor(dat$term_label, levels = unique(dat$term_label[order(dat$gene_ratio)]))
  dat$tooltip <- sprintf(
    "%s\nGeneRatio %s (%.2f), Count %d, p.adjust %.3g",
    dat$term_label, dat$GeneRatio, dat$gene_ratio, dat$Count, dat$p.adjust
  )

  p <- ggplot2::ggplot(dat, ggplot2::aes(
    x = .data$gene_ratio, y = .data$term_label, size = .data$Count, colour = .data$p.adjust
  ))
  p <- p + if (engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)) {
    ggiraph::geom_point_interactive(ggplot2::aes(tooltip = .data$tooltip, data_id = .data$term_label))
  } else {
    ggplot2::geom_point()
  }
  p <- p +
    ggplot2::scale_colour_gradient(low = "#c0392b", high = "#2980b9", name = "p.adjust") +
    ggplot2::scale_size(range = c(1.5, 6), name = "Count") +
    ggplot2::labs(
      title = paste0("Enrichment -- ", scope_label, " / ", db_label),
      subtitle = paste0(
        "Top ", top_n, " most-significant term(s) per (condition, db); x = GeneRatio, y ordered by GeneRatio;\n",
        "colour = p.adjust (red = more significant), size = Count"
      ),
      x = "GeneRatio", y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"), plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"))
  if (multi_db) p <- p + ggplot2::facet_wrap(~db, scales = "free_y")

  .plot_finish(
    proj, p,
    name = "enrichment_plot_log",
    filename = paste0("enrichment_", scope_label, "_", db_label, ".png"),
    log_row = data.frame(condition = scope_label, db = db_label, path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("condition", "db"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}

#' Parse a `clusterProfiler` `"k/n"` ratio string (`GeneRatio`/`BgRatio`)
#'
#' @description
#' `"3/20"` -> `0.15`. Returns `NA_real_` for anything that is not exactly
#' two `/`-separated non-negative integers, or whose denominator is `0`,
#' rather than erroring -- a malformed or missing ratio should drop that
#' row from the plot, not abort the whole call.
#'
#' @param x Character vector of `"k/n"` strings.
#' @return Numeric vector, same length as `x`.
#' @keywords internal
.enrichment_parse_ratio <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, "/", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) != 2L) return(NA_real_)
    num <- suppressWarnings(as.numeric(p[1]))
    den <- suppressWarnings(as.numeric(p[2]))
    if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
    num / den
  }, numeric(1))
}
