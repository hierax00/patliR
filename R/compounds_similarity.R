#' @include AllGenerics.R internal.R
NULL

## compounds_similarity() -- patliR_manual.md, "Planeado / no implementado"
## ("similitud química").
## fingerprint::fp.sim.matrix()'s real signature/return shape (a plain N x N
## symmetric matrix, diag = 1, no dimnames set) confirmed against its own
## source (fingerprint_3.5.10, R/matrix.R) before writing this, not
## assumed from documentation alone -- same discipline already applied to
## pheatmap()/GOplot::GOChord() (see R/plot_heatmap.R, R/plot_gochord.R).
## `fingerprint` is a hard dependency of `rcdk` (already an Imports of
## patliR), so it is always installed alongside it -- listed here in
## Suggests anyway (CRAN policy: a direct `pkg::fun()` call needs its
## package declared, even if guaranteed present transitively) and guarded
## with the same requireNamespace() pattern every other Suggests-gated
## function in this package uses.

#' Pairwise structural similarity between compounds
#'
#' @description
#' Computes a molecular fingerprint (via `rcdk::get.fingerprint()`) for
#' every compound and a pairwise similarity score (via
#' `fingerprint::fp.sim.matrix()`) between every pair -- the kind of
#' structural-similarity analysis behind Yıldırım et al. 2007's figures
#' explaining drug pairs that share a target (see `patliR_manual.md`,
#' family `network_*`). Purely structural (2D fingerprints); says nothing
#' about shared targets or biological activity on its own -- pair this
#' with [network_degeneracy()] (pathway overlap) or `network_edges`
#' (shared targets) for that.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id` to
#'   restrict to, or `NULL` (default) for every compound with a valid
#'   structure.
#' @param fingerprint_type One of `"standard"` (default, hashed path-based
#'   fingerprint), `"extended"`, `"circular"`, `"maccs"`, `"pubchem"` --
#'   passed straight to `rcdk::get.fingerprint(type = ...)`.
#' @param method One of `"tanimoto"` (default), `"dice"`, `"cosine"` --
#'   passed straight to `fingerprint::fp.sim.matrix(method = ...)`.
#'
#' @return The updated `proj`, with a `compounds_similarity` entry in
#'   [patliRResults()] (columns `compound_a`, `compound_b`, `similarity`,
#'   `fingerprint_type`, `method`; one row per unordered pair, self-pairs
#'   excluded), also written to `results/compounds_similarity.csv`.
#'   Compounds whose SMILES does not parse, or whose fingerprint could not
#'   be computed, are excluded and logged -- never silently dropped from
#'   view.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- compounds_similarity(proj)
#' patliRResults(proj, "compounds_similarity")
#' }
#'
#' @export
compounds_similarity <- function(proj, compound_ids = NULL,
                                  fingerprint_type = c("standard", "extended", "circular", "maccs", "pubchem"),
                                  method = c("tanimoto", "dice", "cosine")) {
  stopifnot(is(proj, "PatliRProject"))
  fingerprint_type <- match.arg(fingerprint_type)
  method <- match.arg(method)
  if (!requireNamespace("fingerprint", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn compounds_similarity} needs the {.pkg fingerprint} package.",
      "i" = "{.code install.packages(\"fingerprint\")} -- it is a dependency of {.pkg rcdk}, already required by this package, so it should already be installed in most setups."
    ))
  }

  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  if (nrow(cmp) < 2) {
    cli::cli_abort("Need at least 2 compounds to compute pairwise similarity; found {nrow(cmp)}.")
  }

  mols <- .parse_smiles_safe(cmp$smiles)
  parsed <- !vapply(mols, is.null, logical(1))
  unparsed <- cmp$id[!parsed]
  if (length(unparsed) > 0) {
    proj <- .log_append(
      proj, step = "compounds_similarity", id = unparsed,
      message = "compounds_similarity_excluded: SMILES could not be parsed by rcdk; excluded from the similarity matrix"
    )
  }
  cmp <- cmp[parsed, , drop = FALSE]
  mols <- mols[parsed]

  fps <- lapply(mols, function(m) tryCatch(rcdk::get.fingerprint(m, type = fingerprint_type), error = function(e) NULL))
  has_fp <- !vapply(fps, is.null, logical(1))
  unfingerprinted <- cmp$id[!has_fp]
  if (length(unfingerprinted) > 0) {
    proj <- .log_append(
      proj, step = "compounds_similarity", id = unfingerprinted,
      message = paste0("compounds_similarity_excluded: could not compute a '", fingerprint_type, "' fingerprint for this compound; excluded")
    )
  }
  cmp <- cmp[has_fp, , drop = FALSE]
  fps <- fps[has_fp]

  if (nrow(cmp) < 2) {
    cli::cli_abort("Fewer than 2 compounds with a computable fingerprint; cannot compute similarity.")
  }

  sim_matrix <- fingerprint::fp.sim.matrix(fps, method = method)
  rownames(sim_matrix) <- cmp$id
  colnames(sim_matrix) <- cmp$id

  pairs <- utils::combn(cmp$id, 2, simplify = FALSE)
  result <- do.call(rbind, lapply(pairs, function(p) {
    data.frame(
      compound_a = p[1], compound_b = p[2],
      similarity = sim_matrix[p[1], p[2]],
      fingerprint_type = fingerprint_type, method = method,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL

  patliRResults(proj, "compounds_similarity") <- result
  .write_results_csv(proj, "compounds_similarity", result)
  proj <- .log_append(
    proj, step = "compounds_similarity", id = NA_character_,
    message = paste0(
      "computed pairwise ", method, " similarity for ", nrow(cmp), " compounds (",
      nrow(result), " pairs), fingerprint_type='", fingerprint_type, "'"
    )
  )
  .write_log_csv(proj)
  proj
}
