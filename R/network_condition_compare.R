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
#' @param permutation_test `FALSE` (default). A compound's proximity z-score
#'   does not depend on which condition it is grouped under (same targets,
#'   same disease module, same interactome) -- so a condition's `mean_z`/
#'   `median_z` is just the mean/median of the fixed per-compound scores of
#'   the compounds present in it, and asking "is this condition's aggregate
#'   unusual" is a question about that subset, not a fresh hypothesis test.
#'   When `TRUE`, each compound's z is first averaged over every row it
#'   contributes across `conditions` (collapsing the Monte Carlo noise
#'   [network_proximity()] adds on each rerun of the same compound), then
#'   each condition's observed mean/median is compared to `n_perm` random
#'   subsets of the same size drawn without replacement from that pooled,
#'   per-compound distribution -- an empirical two-sided p-value, adjusted
#'   across `conditions` with Benjamini-Hochberg. This is deliberately a
#'   descriptive add-on, not a filter: see Value.
#' @param n_perm Number of permutation draws when `permutation_test = TRUE`
#'   (default `10000`). Ignored otherwise.
#' @param seed `NULL` (default) or a single integer, for the permutation
#'   draws when `permutation_test = TRUE`. Logged either way (an
#'   auto-generated seed if `NULL`) and restored afterwards, like
#'   [network_proximity()]'s own `seed`. Ignored otherwise.
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
#'   `n_disease_genes_total`, `disease_gene_coverage` (hit / total),
#'   `disease_gene_source`, `score_threshold`, `species`, `string_version`,
#'   `network_type` (this call's resolved proximity provenance, stamped onto
#'   every row -- comparing conditions computed with different values here,
#'   whether in one call or across separate calls for the same disease, is
#'   what this function's own consistency check aborts on). Sorted
#'   with the most negative (closest-to-disease) `aggregate` z first. When
#'   `permutation_test = TRUE`, also gains `perm_null_mean`, `perm_null_sd`,
#'   `perm_p_value`, `perm_p_adjusted`, `perm_n_draws`, `perm_seed_used` --
#'   read these as "how surprising is this condition's subset of compounds",
#'   never as a per-compound significance claim, and note that a condition
#'   compared here with `conditions = NULL` or a small `conditions` vector
#'   only borrows a pool from those conditions, not from every condition
#'   that ever existed for this disease.
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
                                       target_summary = TRUE,
                                       permutation_test = FALSE, n_perm = 10000,
                                       seed = NULL, ...) {
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
  provenance_cols <- c("disease_gene_source", "score_threshold", "species", "string_version", "network_type")
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

  ## The check above only guards THIS call's own conditions. A separate,
  ## later call comparing a different subset of conditions for the SAME
  ## disease (e.g. proximity recomputed at a different STRING threshold or
  ## network_type for just one of them) would each pass their own internal
  ## check and still end up stored side by side -- silently combinable by
  ## plot_condition_compare()/plot_condition_trend(), which read every row
  ## for a disease at once. Guard the write itself against whatever is
  ## already stored for conditions this call is NOT about to replace.
  existing_cmp <- patliRResults(proj, "network_condition_compare")
  if (!is.null(existing_cmp) && nrow(existing_cmp) > 0) {
    other <- existing_cmp[existing_cmp$disease_id == disease & !existing_cmp$condition %in% conditions, , drop = FALSE]
    shared_cols <- intersect(provenance_cols, names(other))
    if (nrow(other) > 0 && length(shared_cols) > 0) {
      mism <- Filter(function(col) {
        ov <- unique(as.character(other[[col]])); pv <- unique(as.character(prox[[col]]))
        length(ov) > 1 || !identical(ov, pv)
      }, shared_cols)
      if (length(mism) > 0) {
        cli::cli_abort(c(
          "Condition(s) {.val {unique(other$condition)}}, already stored for {.val {disease}}, have different {.val {mism}} than this call.",
          "i" = "Comparing or plotting conditions together across separate network_condition_compare() calls is only valid when they share the same proximity provenance.",
          "i" = "Re-run the stored condition(s) with matching arguments, or restrict {.arg conditions} to just the ones computed the same way."
        ))
      }
    }
  }

  dg <- .condition_compare_disease_genes(proj, disease, source = disease_gene_source)

  ## resolved provenance (already validated as unique across `conditions`
  ## above), stamped onto every output row so a LATER call comparing a
  ## different condition subset for the same disease -- or a plotting
  ## consumer reading the whole stored table at once -- can detect a
  ## mismatch instead of silently combining incompatible interactomes.
  provenance_row <- stats::setNames(
    lapply(provenance_cols, function(col) if (nrow(prox) > 0) unique(as.character(prox[[col]]))[1] else NA_character_),
    provenance_cols
  )

  rows <- lapply(conditions, function(cond) {
    p <- prox[prox$condition == cond, , drop = FALSE]
    z <- p$z_score[is.finite(p$z_score)]
    ct <- edges_all[edges_all$condition == cond, , drop = FALSE]
    hit <- unique(ct$uniprot_id[ct$uniprot_id %in% dg$uniprot_id])
    out <- data.frame(
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
    for (col in names(provenance_row)) out[[col]] <- provenance_row[[col]]
    out
  })
  result <- do.call(rbind, rows)

  if (permutation_test) {
    stopifnot(is.numeric(n_perm), length(n_perm) == 1, is.finite(n_perm), n_perm >= 1)
    perm <- .condition_permutation_test(prox, edges_all, conditions,
                                        aggregate = aggregate, n_perm = as.integer(n_perm), seed = seed)
    result <- merge(result, perm, by = "condition", all.x = TRUE, sort = FALSE)
  }

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

#' Subset-permutation test: is a condition's aggregate z unusual relative to
#' a same-size random subset of the pooled per-compound z?
#' @param prox `network_proximity` rows already filtered to `conditions` and
#'   the disease being compared (as built by [network_condition_compare()]).
#' @param edges_all `network_edges` (all conditions -- filtered here per
#'   `cond`, matching how [network_condition_compare()] itself computes
#'   `n_compounds_present`).
#' @return `data.frame(condition, perm_null_mean, perm_null_sd,
#'   perm_p_value, perm_p_adjusted, perm_n_draws, perm_seed_used)`.
#' @keywords internal
.condition_permutation_test <- function(prox, edges_all, conditions, aggregate, n_perm, seed) {
  stat_fn <- if (aggregate == "mean") mean else stats::median
  finite <- prox[is.finite(prox$z_score), , drop = FALSE]
  ## one z per compound: averages away the Monte Carlo noise network_proximity()
  ## adds each time the same physical compound is scored in a different condition
  zc <- tapply(finite$z_score, finite$compound_id, mean)

  ## snapshot the caller's RNG state *before* touching it -- if `seed = NULL`,
  ## sample.int() below draws from (and advances) this same stream to pick
  ## used_seed, so capturing old_seed any later would restore a state one
  ## draw past the caller's real starting point
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  used_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else as.integer(seed)
  set.seed(used_seed)
  on.exit({
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  rows <- lapply(conditions, function(cond) {
    ids <- intersect(unique(edges_all$compound_id[edges_all$condition == cond]), names(zc))
    if (length(ids) == 0 || length(zc) == 0) {
      return(data.frame(condition = cond, perm_null_mean = NA_real_, perm_null_sd = NA_real_,
                        perm_p_value = NA_real_, stringsAsFactors = FALSE))
    }
    obs <- stat_fn(zc[ids])
    ## zc[sample.int(length(zc), length(ids))], NOT sample(zc, length(ids)):
    ## base R's sample() treats a numeric vector of length 1 as "sample from
    ## 1:that value" (?sample's documented surprise for a length-one x), so
    ## with a single-compound pool sample(zc, 1) would draw from 1:zc[[1]]
    ## instead of returning zc[[1]] itself, silently producing a null
    ## distribution unrelated to the real (degenerate, single-value) one.
    null <- replicate(n_perm, stat_fn(zc[sample.int(length(zc), length(ids))]))
    lo <- (1 + sum(null <= obs)) / (1 + n_perm)
    hi <- (1 + sum(null >= obs)) / (1 + n_perm)
    data.frame(condition = cond, perm_null_mean = mean(null), perm_null_sd = stats::sd(null),
              perm_p_value = min(1, 2 * min(lo, hi)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$perm_p_adjusted <- stats::p.adjust(out$perm_p_value, "BH")
  out$perm_n_draws <- n_perm
  out$perm_seed_used <- used_seed
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
