#' @include AllGenerics.R internal.R network_build.R network_enrich.R
NULL

## Renders KEGG pathway diagrams (pathview::pathview()) for a condition's
## top enrichKEGG() hits, coloring each gene node by a score from
## network_build()'s compound-target edges (`gene_score = "max_weight"` by
## default: the max import probability among compounds hitting that
## target).
##
## pathview() has NO argument for where it writes its output PNG --
## `kegg.native = TRUE` (the default) always writes
## `<pathway.id>.<out.suffix>.png` to the current working directory. So
## this function temporarily setwd()s into `out_dir` for the pathview()
## calls, restored via on.exit() even on error. `kegg.dir` is where its
## *input* KGML/PNG files are cached (only re-downloaded when absent),
## pointed at cacheDir(proj).

#' Render KEGG pathway diagrams (`pathview`) for a condition's top enriched
#' pathways, colored by target score
#'
#' @description
#' For one `condition` already enriched against KEGG by [network_enrich()]
#' (`db = "kegg"`), renders a native KEGG pathway PNG for each of its
#' `top_n_pathways` most significant pathways (or an explicit
#' `pathway_id` vector), using [pathview::pathview()]. Each pathway's
#' target/gene nodes are colored by `gene_score` -- a per-Entrez-ID summary
#' of that condition's own [network_build()] compound-target edges.
#'
#' @section Why this needs its own gene-to-Entrez mapping, not `geneID` off `network_enrichment`:
#' [network_enrich()]'s output stores the Entrez IDs it fed into
#' `enrichKEGG()` only as a single `/`-separated string per pathway
#' (`geneID`, `clusterProfiler`'s own format) -- convenient for display,
#' not for building a numeric-scored named vector. This function instead
#' re-maps `condition`'s own [network_build()] UniProt targets to Entrez
#' (`clusterProfiler::bitr()`, same call [network_enrich()] itself makes)
#' so it can attach a real per-gene score, not just presence/absence.
#'
#' @section Output files, not an in-memory plot object:
#' Unlike the `plot_*` family, `pathview()` does not return anything this
#' package would call a "plot object" -- it renders and writes PNG files as
#' its side effect. `network_pathview()` follows the `network_*` family's
#' own convention instead (input `proj` -> output `proj`, every result
#' written to disk and logged), not `plot_*`'s (input `proj` -> returned
#' `ggplot`/`girafe`).
#'
#' @section Heavy, optional dependencies:
#' Needs `pathview` (Bioconductor, **not** installed by default -- add via
#' `BiocManager::install("pathview")`) plus `clusterProfiler`/`org.Hs.eg.db`
#' (already needed by
#' [network_enrich()]). `pathview` downloads each pathway's KGML/PNG from
#' KEGG's REST API on first use (needs internet, like `db = "kegg"` in
#' [network_enrich()] itself) and caches them under `kegg_dir` afterwards
#' (default `cacheDir(proj)/kegg_pathview`) -- deleting that cache is
#' always safe, same principle as every other cache in `patliR`.
#'
#' @inheritParams network_build
#' @param condition Character scalar -- a single condition already enriched
#'   against KEGG by [network_enrich()]. Unlike most `network_*` functions,
#'   `NULL` (default) only works if exactly one condition has `db = "kegg"`
#'   enrichment results; otherwise pass it explicitly (same rule
#'   [plot_gochord()] uses, for the same reason -- pathway diagrams are
#'   inherently per-condition, pooling them is not meaningful).
#' @param gene_score `"max_weight"` (default -- the highest import
#'   probability among compounds hitting that target in `condition`),
#'   `"mean_weight"`, or
#'   `"n_compounds"` (count of distinct compounds hitting that target --
#'   not on a 0-1 scale, `limit` is adjusted automatically, see below).
#' @param pathway_id Character vector of KEGG pathway IDs (e.g.
#'   `"hsa04066"`, as returned in `network_enrichment$ID`), or `NULL`
#'   (default) to use the `top_n_pathways` most significant KEGG results
#'   for `condition` (lowest `p.adjust`). If given, every ID must already
#'   appear in `condition`'s `db = "kegg"` [network_enrich()] results.
#' @param top_n_pathways Integer, default `10`. Ignored if `pathway_id` is given.
#' @param low,mid,high Colors passed straight to `pathview::pathview()`'s
#'   `low`/`mid`/`high` (each internally wrapped as `list(gene = ...)`,
#'   `cpd.data` is never used here). Default `"white"`/`"yellow"`/`"red"`,
#'   not `pathview`'s own default (`"green"`/`"gray"`/`"red"`).
#' @param out_dir Directory the rendered PNGs are written to. Defaults to
#'   `file.path(projectDir(proj), "plots", "kegg_pathview")`.
#' @param kegg_dir Directory `pathview` reads/caches each pathway's
#'   downloaded `.xml`/`.png` from. Defaults to
#'   `file.path(cacheDir(proj), "kegg_pathview")`.
#'
#' @return The updated `proj`, with a `kegg_pathview_log` entry in
#'   [patliRResults()] (columns `condition`, `pathway_id`, `description`,
#'   `gene_score`, `n_genes_mapped`, `ok`, `path`, `message`), also written
#'   to `results/kegg_pathview_log.csv`. A pathway that fails to render
#'   (e.g. `pathview` cannot reach KEGG, or the pathway has no mappable
#'   nodes) gets `ok = FALSE` and a `message`, is logged
#'   (`"network_pathview_render_failed"`), and does **not** stop the other
#'   pathways in the same call from rendering.
#'
#' @references Luo, W. & Brouwer, C. (2013), "Pathview: an R/Bioconductor
#'   package for pathway-based data integration and visualization",
#'   *Bioinformatics* 29(14), 1830-1831. \doi{10.1093/bioinformatics/btt285}
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
#' proj <- network_enrich(proj, db = "kegg") # needs internet
#' proj <- network_pathview(proj, top_n_pathways = 3) # needs pathview + internet
#' patliRResults(proj, "kegg_pathview_log")
#' }
#'
#' @export
network_pathview <- function(proj, condition = NULL,
                              gene_score = c("max_weight", "mean_weight", "n_compounds"),
                              pathway_id = NULL, top_n_pathways = 10,
                              low = "white", mid = "yellow", high = "red",
                              out_dir = NULL, kegg_dir = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  gene_score <- match.arg(gene_score)
  stopifnot(is.numeric(top_n_pathways), length(top_n_pathways) == 1, top_n_pathways >= 1)
  if (!is.null(pathway_id)) stopifnot(is.character(pathway_id), length(pathway_id) >= 1)

  if (!requireNamespace("pathview", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_pathview} needs the {.pkg pathview} package (Bioconductor).",
      "i" = "Install it with {.code BiocManager::install(\"pathview\")}."
    ))
  }
  .network_enrich_check_deps("kegg")

  enrichment_all <- patliRResults(proj, "network_enrichment")
  if (is.null(enrichment_all) || nrow(enrichment_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_enrichment} entry in {.arg proj}.",
      "i" = "Run {.fn network_enrich} with {.code db = \"kegg\"} first."
    ))
  }
  kegg_all <- enrichment_all[enrichment_all$db == "kegg", , drop = FALSE]
  if (nrow(kegg_all) == 0) {
    cli::cli_abort(c(
      "No {.code db = \"kegg\"} rows in {.val network_enrichment}.",
      "i" = "Run {.fn network_enrich} with {.code db = \"kegg\"} first."
    ))
  }

  if (is.null(condition)) {
    kegg_conditions <- unique(kegg_all$condition)
    if (length(kegg_conditions) != 1) {
      cli::cli_abort(c(
        "{.fn network_pathview} needs a single {.arg condition} (pathway diagrams are per-condition).",
        "i" = "KEGG enrichment is available for {.val {kegg_conditions}}; pass one via {.arg condition}."
      ))
    }
    cond <- kegg_conditions
  } else {
    stopifnot(is.character(condition), length(condition) == 1)
    if (!condition %in% kegg_all$condition) {
      cli::cli_abort("No {.code db = \"kegg\"} {.val network_enrichment} rows for condition {.val {condition}}; available: {.val {unique(kegg_all$condition)}}.")
    }
    cond <- condition
  }

  enr <- kegg_all[kegg_all$condition == cond, , drop = FALSE]
  enr <- enr[order(enr$p.adjust), , drop = FALSE]

  if (!is.null(pathway_id)) {
    unknown <- setdiff(pathway_id, enr$ID)
    if (length(unknown) > 0) {
      cli::cli_abort("Pathway ID(s) {.val {unknown}} are not among condition {.val {cond}}'s {.code db = \"kegg\"} {.fn network_enrich} results.")
    }
    enr <- enr[match(pathway_id, enr$ID), , drop = FALSE]
  } else {
    enr <- utils::head(enr, top_n_pathways)
  }

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c("No {.val network_edges} entry in {.arg proj}.", "i" = "Run {.fn network_build} first."))
  }
  edges_cond <- edges_all[edges_all$condition == cond, , drop = FALSE]
  edges_cond <- edges_cond[!is.na(edges_cond$uniprot_id) & nzchar(edges_cond$uniprot_id), , drop = FALSE]
  if (nrow(edges_cond) == 0) {
    cli::cli_abort("No compound-target edges for condition {.val {cond}} in {.val network_edges}; nothing to score/color.")
  }

  map <- clusterProfiler::bitr(unique(edges_cond$uniprot_id), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
  unmapped <- setdiff(unique(edges_cond$uniprot_id), map$UNIPROT)
  if (length(unmapped) > 0) {
    proj <- .log_append(
      proj, step = "network_pathview", id = unmapped,
      message = "network_pathview_unmapped_target: no Entrez Gene ID found in org.Hs.eg.db for this UniProt ID; excluded from the gene score vector"
    )
  }
  edges_cond <- merge(edges_cond, map, by.x = "uniprot_id", by.y = "UNIPROT")
  if (nrow(edges_cond) == 0) {
    cli::cli_abort("None of condition {.val {cond}}'s targets mapped to an Entrez Gene ID; cannot build a gene score vector.")
  }

  gene_vector <- .network_pathview_gene_vector(edges_cond, gene_score)
  limit_gene <- if (gene_score == "n_compounds") c(0, max(gene_vector, 1)) else c(0, 1)

  if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots", "kegg_pathview")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(kegg_dir)) kegg_dir <- file.path(cacheDir(proj), "kegg_pathview")
  if (!dir.exists(kegg_dir)) dir.create(kegg_dir, recursive = TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir) # pathview() writes its native-KEGG PNG to the working
                 # directory unconditionally -- no output-path argument
                 # exists (verified against pathview's real source, see
                 # file header comment); this is the only way to control it.

  rows <- vector("list", nrow(enr))
  for (i in seq_len(nrow(enr))) {
    pid <- enr$ID[i]
    res <- tryCatch(
      pathview::pathview(
        gene.data = gene_vector, pathway.id = pid, species = "hsa",
        kegg.dir = kegg_dir, limit = list(gene = limit_gene),
        low = low, mid = mid, high = high
      ),
      error = function(e) e
    )
    ok <- !inherits(res, "error")
    if (!ok) {
      proj <- .log_append(
        proj, step = "network_pathview", id = pid,
        message = paste0("network_pathview_render_failed: ", conditionMessage(res))
      )
    }
    rows[[i]] <- data.frame(
      condition = cond, pathway_id = pid, description = enr$Description[i],
      gene_score = gene_score, n_genes_mapped = length(gene_vector), ok = ok,
      path = file.path(out_dir, paste0(pid, ".pathview.png")),
      message = if (ok) NA_character_ else conditionMessage(res),
      stringsAsFactors = FALSE
    )
  }
  setwd(old_wd)

  new_rows <- do.call(rbind, rows)
  rownames(new_rows) <- NULL
  new_rows <- .network_upsert(proj, "kegg_pathview_log", new_rows, c("condition", "pathway_id"))

  patliRResults(proj, "kegg_pathview_log") <- new_rows
  .write_results_csv(proj, "kegg_pathview_log", new_rows)
  proj <- .log_append(
    proj, step = "network_pathview", id = NA_character_,
    message = paste0("condition '", cond, "': ", sum(new_rows$condition == cond & new_rows$ok), "/",
                      sum(new_rows$condition == cond), " KEGG pathway diagram(s) rendered to ", out_dir)
  )
  .write_log_csv(proj)
  proj
}

#' Aggregate a condition's (Entrez-mapped) compound-target edges into one
#' named numeric vector, per `gene_score`
#' @return Named numeric vector (`names` = Entrez IDs).
#' @keywords internal
.network_pathview_gene_vector <- function(edges_cond, gene_score) {
  if (gene_score == "n_compounds") {
    agg <- stats::aggregate(compound_id ~ ENTREZID, edges_cond, function(x) length(unique(x)))
    return(stats::setNames(agg$compound_id, agg$ENTREZID))
  }
  fun <- if (gene_score == "max_weight") max else mean
  agg <- stats::aggregate(weight ~ ENTREZID, edges_cond, fun)
  stats::setNames(agg$weight, agg$ENTREZID)
}
