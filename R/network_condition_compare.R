#' @include AllGenerics.R internal.R network_build.R network_proximity.R
NULL

#' Compare conditions (extracts) by how close their compounds sit to a disease
#'
#' @description
#' [network_proximity()] scores one compound at a time. This aggregates those
#' per-compound z-scores to a single number **per condition**, so several
#' extracts (e.g. different phenological stages, greenhouse vs. wild) can be
#' ranked by how strong a signal their whole compound set gives for a disease
#' -- answering "which extract has the best effect?" rather than "which
#' compound".
#'
#' Any `(condition, disease)` pair not yet in `patliRResults(proj,
#' "network_proximity")` is computed first with [network_proximity()], using
#' `...`. Already-computed pairs are reused as-is (no recomputation), so
#' calling this repeatedly while conditions are added one at a time is cheap.
#'
#' @inheritParams network_build
#' @param conditions Character vector of conditions to compare, or `NULL`
#'   (default) for every condition currently in `network_edges` (a pooled
#'   condition is included if present -- exclude it explicitly if that is
#'   not wanted).
#' @param disease Character EFO/MONDO/GO id, as used by [network_proximity()].
#' @param aggregate `"mean"` (default) or `"median"` z-score across the
#'   compounds present in a condition.
#' @param target_summary `TRUE` (default): also store a per-condition,
#'   per-target contribution table (see Value).
#' @param ... Passed to [network_proximity()] for any pair not yet computed
#'   (`score_threshold`, `n_random`, `seed`, `disease_genes`, ...).
#'
#' @return The updated `proj`. `patliRResults(proj,
#'   "network_condition_compare")` gains/updates rows keyed by
#'   `(condition, disease_id)`: `n_compounds_present` (in `network_edges`),
#'   `n_compounds_scored` (rows in `network_proximity` for this condition,
#'   whether or not each one's z-score came out finite), `n_compounds_finite`
#'   (of those, how many did -- this is what feeds `mean_z`/`median_z`/`sd_z`/`best_z`),
#'   `best_z` (most negative, i.e. closest single compound), `n_disease_genes_hit`,
#'   `n_disease_genes_total`, `disease_gene_coverage` (hit / total). Sorted
#'   with the most negative (closest-to-disease) `aggregate` z first.
#'   When `target_summary = TRUE`, also sets `patliRResults(proj,
#'   "network_condition_targets")`: one row per `(condition, disease_id,
#'   uniprot_id)` hit by that condition's compounds and also a disease gene,
#'   with `n_compounds` (how many of the condition's compounds hit it),
#'   `gene_symbol`, `association_score`, and, when available,
#'   `degree`/`hub_score` from [network_centrality()] -- the targets driving
#'   that condition's proximity signal, ranked by `n_compounds` then
#'   `association_score`.
#'
#' @seealso [network_proximity()], [plot_condition_compare()]
#' @export
network_condition_compare <- function(proj, conditions = NULL, disease,
                                       aggregate = c("mean", "median"),
                                       target_summary = TRUE, ...) {
  stopifnot(is(proj, "PatliRProject"))
  aggregate <- match.arg(aggregate)
  stopifnot(is.character(disease), length(disease) == 1, nzchar(disease))

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c("No {.val network_edges} entry in {.arg proj}.", "i" = "Run {.fn network_build} first."))
  }
  if (is.null(conditions)) conditions <- unique(edges_all$condition)
  conditions <- as.character(conditions)

  prox <- patliRResults(proj, "network_proximity")
  have <- if (is.null(prox)) character(0) else unique(prox$condition[prox$disease_id == disease])
  missing_cond <- setdiff(conditions, have)
  for (cond in missing_cond) {
    proj <- network_proximity(proj, condition = cond, disease = disease, ...)
  }
  prox <- patliRResults(proj, "network_proximity")
  prox <- prox[prox$condition %in% conditions & prox$disease_id == disease, , drop = FALSE]

  ## every condition being compared must share the same proximity provenance --
  ## mixing e.g. score_threshold = 400 for one condition and 700 (from a stale
  ## or differently-parameterised run) for another would rank them on
  ## different, incomparable null models without any visible warning
  provenance_cols <- c("disease_gene_source", "score_threshold", "species", "string_version")
  provenance_cols <- intersect(provenance_cols, names(prox))
  inconsistent <- Filter(function(col) length(unique(prox[[col]])) > 1, provenance_cols)
  if (length(inconsistent) > 0) {
    cli::cli_abort(c(
      "The conditions being compared have different {.val {inconsistent}} in their stored {.val network_proximity} rows.",
      "i" = paste0("Comparing conditions is only valid when every condition's proximity was computed the same way. ",
                   "Recompute the mismatched condition(s) with matching arguments, or drop them from {.arg conditions}.")
    ))
  }
  disease_gene_source <- if ("disease_gene_source" %in% names(prox) && nrow(prox) > 0) {
    prox$disease_gene_source[1]
  } else {
    "disease_genes"
  }

  dg <- .condition_compare_disease_genes(proj, disease, source = disease_gene_source)

  rows <- lapply(conditions, function(cond) {
    p <- prox[prox$condition == cond, , drop = FALSE]
    z <- p$z_score[is.finite(p$z_score)]
    ct <- edges_all[edges_all$condition == cond, , drop = FALSE]
    hit <- unique(ct$uniprot_id[ct$uniprot_id %in% dg$uniprot_id])
    data.frame(
      condition = cond, disease_id = disease,
      n_compounds_present = length(unique(ct$compound_id)),
      n_compounds_scored = nrow(p), n_compounds_finite = length(z),
      mean_z = if (length(z)) mean(z) else NA_real_,
      median_z = if (length(z)) stats::median(z) else NA_real_,
      sd_z = if (length(z) > 1) stats::sd(z) else NA_real_,
      best_z = if (length(z)) min(z) else NA_real_,
      n_disease_genes_hit = length(hit), n_disease_genes_total = nrow(dg),
      disease_gene_coverage = if (nrow(dg) > 0) length(hit) / nrow(dg) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  sort_col <- if (aggregate == "mean") result$mean_z else result$median_z
  result <- result[order(sort_col), , drop = FALSE]
  rownames(result) <- NULL

  result <- .network_upsert(proj, "network_condition_compare", result, c("condition", "disease_id"))
  patliRResults(proj, "network_condition_compare") <- result
  .write_results_csv(proj, "network_condition_compare", result)

  if (target_summary) {
    targets <- .condition_compare_targets(proj, conditions, disease, edges_all, dg)
    ## coarser key than the row itself (like network_proximity's own touched_keys, see
    ## its "Recompute the whole condition/disease" comment): a rerun must wipe every old
    ## uniprot_id row for a touched (condition, disease_id), including ones that no
    ## longer qualify (lost all edges, or their target is no longer a disease gene),
    ## not just overwrite the uniprot_ids the new run happens to still produce.
    cmp_touched <- data.frame(condition = conditions, disease_id = disease, stringsAsFactors = FALSE)
    targets <- .network_upsert(proj, "network_condition_targets", targets,
                               c("condition", "disease_id"), touched_keys = cmp_touched)
    patliRResults(proj, "network_condition_targets") <- targets
    .write_results_csv(proj, "network_condition_targets", targets)
  }

  .write_log_csv(proj)
  proj
}

#' One row per UniProt accession for a disease, from the source
#' [network_proximity()] actually used
#' @param source `"disease_genes"` or `"targets_disease"` (matches
#'   [network_proximity()]'s `disease_genes` argument and its stored
#'   `disease_gene_source` column) -- which table to read, not a free choice:
#'   reading the other one here would summarise a different disease module
#'   than the one the z-scores were actually computed against.
#' @return `data.frame(uniprot_id, gene_symbol, association_score)`, one row
#'   per `uniprot_id` (an accession appearing more than once keeps its
#'   highest `association_score`, so a duplicate never inflates a coverage
#'   denominator or duplicates a downstream join).
#' @keywords internal
.condition_compare_disease_genes <- function(proj, disease, source = "disease_genes") {
  out <- if (source == "targets_disease") {
    td <- patliRResults(proj, "targets_disease")
    if (is.null(td) || !any(td$disease_id == disease)) {
      cli::cli_abort(c(
        "No rows for {.val {disease}} in {.val targets_disease} (the source {.fn network_proximity} used).",
        "i" = "Run {.fn targets_disease_filter} first."
      ))
    }
    td <- td[td$disease_id == disease, , drop = FALSE]
    data.frame(uniprot_id = td$target_id, gene_symbol = NA_character_, association_score = NA_real_,
              stringsAsFactors = FALSE)
  } else {
    dg <- patliRResults(proj, "disease_genes")
    if (is.null(dg) || !any(dg$disease_id == disease)) {
      cli::cli_abort(c(
        "No rows for {.val {disease}} in {.val disease_genes} (the source {.fn network_proximity} used).",
        "i" = "Run {.fn disease_genes_fetch} first."
      ))
    }
    dg[dg$disease_id == disease, c("uniprot_id", "gene_symbol", "association_score"), drop = FALSE]
  }
  ## one row per accession: keep the highest association_score (NA loses to any number)
  out <- out[order(out$uniprot_id, -out$association_score, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out$uniprot_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @keywords internal
.condition_compare_targets <- function(proj, conditions, disease, edges_all, dg) {
  cent <- patliRResults(proj, "network_centrality")
  rows <- lapply(conditions, function(cond) {
    ct <- edges_all[edges_all$condition == cond, , drop = FALSE]
    ct <- ct[ct$uniprot_id %in% dg$uniprot_id, , drop = FALSE]
    if (nrow(ct) == 0) {
      return(data.frame(condition = character(0), disease_id = character(0), uniprot_id = character(0),
                        gene_symbol = character(0), association_score = numeric(0), n_compounds = integer(0),
                        degree = numeric(0), hub_score = numeric(0), stringsAsFactors = FALSE))
    }
    agg <- stats::aggregate(compound_id ~ uniprot_id, ct, function(x) length(unique(x)))
    names(agg) <- c("uniprot_id", "n_compounds")
    agg <- merge(agg, dg, by = "uniprot_id", all.x = TRUE)
    agg$degree <- NA_real_; agg$hub_score <- NA_real_
    if (!is.null(cent)) {
      cc <- cent[cent$condition == cond & cent$node_type == "target", , drop = FALSE]
      m <- match(agg$uniprot_id, cc$node_id)
      if ("degree" %in% names(cc)) agg$degree <- cc$degree[m]
      if ("hub_score" %in% names(cc)) agg$hub_score <- cc$hub_score[m]
    }
    agg <- agg[order(-agg$n_compounds, -agg$association_score), , drop = FALSE]
    data.frame(condition = cond, disease_id = disease, agg, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
