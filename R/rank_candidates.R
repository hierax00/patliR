#' @include AllGenerics.R internal.R network_build.R network_centrality.R network_hub_penalty.R network_proximity.R network_synergy.R network_module_robustness.R adme_filter.R
NULL

## Design: chilcuague-analysis/reviews/rank-candidates-design-spec.md (Milestone
## B). Combines up to 6 per-compound criteria via Robust Rank Aggregation
## (RRA, Kolde et al. 2012) -- 2 mandatory (adme, centrality), 4 optional
## and auto-detected from patliRResults(proj) (hub_penalty, proximity,
## synergy, module_robustness). RRA over a weighted linear sum of z-scored
## criteria because a criterion with a broken/skewed scale (e.g.
## betweenness_norm sitting at 0 for most nodes) cannot silently dominate --
## ranks are compared, not raw magnitudes. RobustRankAggreg::aggregateRanks()
## (CRAN, still current as of this writing) rather than hand-rolled: its
## null model is a closed-form beta-distribution order-statistic correction,
## not geometry -- the same "depend for algorithms" argument DESIGN.md
## already makes for `bipartite`'s Barber's Q_B optimisation (piece 13).
## Supply normalised ranks directly via `rmat`: rows are compounds and
## columns are criteria, with average ranks for ties and missing values
## imputed to relative rank 1. No ordered-list conversion is needed.
## Toxicity and bias_reweighted were both considered as criteria and
## deliberately dropped (user's call, 2026-09-13) -- not offered even as an
## opt-in override in this version.

#' Rank candidate compounds by Robust Rank Aggregation over network/ADME
#' criteria
#'
#' @description
#' Combines up to 6 per-compound criteria into one ranked table using
#' Robust Rank Aggregation (RRA; Kolde, Laur, Adler & Vilo 2012,
#' *Bioinformatics* 28(4):573-580). RRA compares each criterion's *rank*
#' ordering, not its raw values, so a criterion with a broken or skewed
#' scale cannot silently dominate the others -- the standard failure mode
#' of a weighted linear sum of z-scored criteria. Two criteria are
#' mandatory (`"adme"`, `"centrality"`); the rest (`"hub_penalty"`,
#' `"proximity"`, `"synergy"`, `"module_robustness"`) are auto-detected from
#' [patliRResults()] and silently degrade -- a project that never ran
#' [network_synergy()] still gets a ranking, just over fewer criteria, and
#' [projectLog()] records exactly which criteria were used vs skipped for
#' each condition.
#'
#' @section Criteria:
#' Direction is already normalised to "higher = better" internally for
#' ranking purposes; the *raw* value, in its own natural direction, is what
#' the output table actually shows.
#' \describe{
#'   \item{`"adme"` (mandatory)}{`crit_adme_pass_frac` -- fraction of
#'     [adme_filter()]'s selected rules a compound passes, from
#'     `adme_filtered`. Aborts if `adme_filtered` is absent (run
#'     [adme_local()] then [adme_filter()] first).}
#'   \item{`"centrality"` (mandatory)}{`crit_centrality` -- for each of a
#'     compound's own predicted targets, `rowMeans()` of whichever of
#'     `degree_norm` / `betweenness_norm` / `hub_score_component` are
#'     non-`NA` in `network_centrality`, then rolled up to the compound via
#'     `roll_up`. Aborts if `network_centrality` was run with
#'     `normalize = FALSE` (no normalised columns to build a composite
#'     from) -- re-run it with `normalize = TRUE`.}
#'   \item{`"hub_penalty"` (optional)}{`crit_hub_penalty` --
#'     `network_hub_penalty`'s `score_adjusted`, rolled up target->compound.}
#'   \item{`"proximity"` (optional)}{`crit_proximity_z` --
#'     `network_proximity`'s `z_score` (lower/more negative is better --
#'     the only criterion whose raw scale runs the opposite way; direction
#'     is flipped internally before ranking).}
#'   \item{`"synergy"` (optional)}{`crit_synergy_best` -- for each compound,
#'     the best (`max`) `synergy_score` among every partner it was scored
#'     against in `network_synergy` (its best achievable complementary
#'     pairing, not an average across all partners). `n_p2_partners` (count
#'     of `cheng_class == "P2"` partnerships) is kept as an *informational*
#'     column only -- not fed to the RRA, since it is a different view of
#'     the same underlying table as `crit_synergy_best` and including both
#'     would double-count one signal.}
#'   \item{`"module_robustness"` (optional)}{`crit_module_r_index` -- the
#'     `r_index` of the module a compound belongs to (direct join via
#'     `network_module_membership`, not an aggregation -- one compound has
#'     exactly one module per condition).}
#' }
#' `"proximity"` and `"synergy"` both key off a `disease_id`; if `disease`
#' is `NULL` and `network_proximity` holds results for exactly one disease
#' for a condition, that one is used automatically -- with more than one,
#' `disease` must be supplied explicitly (same policy [network_synergy()]
#' already uses for its own `disease` argument).
#'
#' @section Roll-up (target->compound, `roll_up=`):
#' `"weighted_mean"` (default) -- weighted by the [network_build()] edge
#' weight (import probability) for each `(compound_id, uniprot_id)` pair,
#' so a compound's score reflects its whole predicted target profile
#' weighted by confidence rather than being hijacked by one high-scoring
#' but low-probability predicted target. `"mean"` and `"max"` are also
#' available, always explicit -- never a hidden default.
#'
#' @section Missing data:
#' Two different situations get two different treatments. A whole
#' *optional* table absent from the project (e.g. [network_synergy()] never
#' run) drops that criterion entirely for the affected condition(s) --
#' logged, not imputed (imputing a constant discount for a criterion nobody
#' computed would be a silent no-op dressed up as data). A table that
#' exists but lacks a value for one particular compound (e.g. a mono-target
#' compound with no [network_synergy()] partner) gets **worst-rank
#' imputation** on that one criterion -- the standard convention for
#' incomplete rank lists under RRA. Every missing value receives rank `N`
#' (relative rank `1`), where `N` is the number of compounds in the condition,
#' so missing evidence cannot improve its rank. `n_criteria_used` records, per compound,
#' how many criteria it had a *real* (non-imputed) value for.
#'
#' @section Output columns not otherwise documented above:
#' `rank_<criterion>` -- one per `crit_*` column actually used for this
#' condition, 1 = best, direction already normalised. Exact ties share the
#' average of their occupied positions (possibly fractional), independent
#' of compound IDs and input order. These ranks divided by `N` are supplied
#' to RRA as a numeric matrix. `rra_score` --
#' `RobustRankAggreg::aggregateRanks()`'s output (lower = more robustly
#' top-ranked across every criterion used). `rra_rank` --
#' the rank of `rra_score` (ties -- common, because RRA scores saturate at 1
#' -- broken by the mean per-criterion rank; still `"min"` for a tie on
#' both), the table's primary sort key.
#' `pareto_front` -- integer tier (`1` = non-dominated) from a from-scratch
#' non-dominated sort (NSGA-II-style; O(n^2) per tier, fine at the
#' candidate-list scale this function targets -- tens to a few hundred
#' compounds) over every `crit_*` column actually used, with a missing
#' value treated as `-Inf` (worst), the same convention as the RRA
#' imputation. See [plot_rank()] for the companion visualisations.
#'
#' @inheritParams network_build
#' @param disease `NULL` (default, auto-resolve when unambiguous) or a
#'   `disease_id` string identifying which [network_proximity()] /
#'   [network_synergy()] run to use for the `"proximity"`/`"synergy"`
#'   criteria.
#' @param criteria `NULL` (default, auto-detect: use every criterion with
#'   data) or a character vector naming a subset of `"adme"`,
#'   `"centrality"`, `"hub_penalty"`, `"proximity"`, `"synergy"`,
#'   `"module_robustness"` to restrict to (the two mandatory ones must be
#'   included).
#' @param roll_up `"weighted_mean"` (default), `"mean"`, or `"max"` -- see
#'   the roll-up section above.
#' @param top_n How many top-`rra_rank` compounds `export` writes out.
#'   Default `15`.
#' @param export `"none"` (default), `"sdf"`, or `"smi"`. When not
#'   `"none"`, writes `results/rank_candidates_top<top_n>_<condition>.<ext>`
#'   per condition from the compounds' own SMILES (via the same
#'   `.parse_smiles_safe()` identity chain [compounds_similarity()] uses --
#'   no re-lookup, no new dependency). `"sdf"` uses
#'   `rcdk::write.molecules()`; `"smi"` writes plain `SMILES<tab>compound_id`
#'   lines (`write.molecules()` always emits an SDF mol block regardless of
#'   the file extension -- verified directly against the installed `rcdk`
#'   before relying on it here -- so `"smi"` is built by hand from
#'   `rcdk::get.smiles()`).
#'
#' @return The updated `proj`, with a `rank_candidates` entry in
#'   [patliRResults()] (columns `condition`, `compound_id`, one `crit_*` per
#'   criterion used anywhere in this call, one paired `rank_*`,
#'   `n_p2_partners`, `n_criteria_used`, `rra_score`, `rra_rank`,
#'   `pareto_front`), also written to `results/rank_candidates.csv`.
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
#' proj <- network_centrality(proj)
#' proj <- adme_local(proj)
#' proj <- adme_filter(proj)
#' proj <- rank_candidates(proj)
#' patliRResults(proj, "rank_candidates")
#' }
#'
#' @references
#' Kolde, Laur, Adler & Vilo (2012), *Bioinformatics* 28(4):573-580,
#' \doi{10.1093/bioinformatics/btr709}.
#'
#' @seealso [plot_rank()]
#' @export
rank_candidates <- function(proj, condition = NULL, disease = NULL,
                             criteria = NULL,
                             roll_up = c("weighted_mean", "mean", "max"),
                             top_n = 15,
                             export = c("none", "sdf", "smi")) {
  stopifnot(is(proj, "PatliRProject"))
  roll_up <- match.arg(roll_up)
  export <- match.arg(export)
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)

  if (!requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg RobustRankAggreg} package is required for {.fn rank_candidates}.",
      "i" = "{.code install.packages(\"RobustRankAggreg\")} -- implements Kolde et al. (2012)'s Robust Rank Aggregation, which {.fn rank_candidates} uses so no single criterion's scale can dominate the ranking."
    ))
  }

  all_criteria <- c("adme", "centrality", "hub_penalty", "proximity", "synergy", "module_robustness")
  mandatory <- c("adme", "centrality")
  if (is.null(criteria)) {
    criteria <- all_criteria
  } else {
    criteria <- match.arg(criteria, all_criteria, several.ok = TRUE)
    missing_mandatory <- setdiff(mandatory, criteria)
    if (length(missing_mandatory) > 0) {
      cli::cli_abort("{.arg criteria} must include the mandatory criterion/criteria {.val {missing_mandatory}}.")
    }
  }

  conditions <- .network_resolve_conditions(proj, condition)
  edges_all <- patliRResults(proj, "network_edges")

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    compound_ids <- sort(unique(edges_all$compound_id[edges_all$condition == cond]))
    if (length(compound_ids) < 2) {
      cli::cli_abort("Condition {.val {cond}} has fewer than 2 present compounds; {.fn rank_candidates} needs at least 2 to rank.")
    }
    edges_cond <- edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id", "weight"), drop = FALSE]

    crit_values <- list()
    used <- character(0)
    skipped <- character(0)
    n_p2_vec <- NULL

    if ("adme" %in% criteria) {
      v <- .rank_criterion_adme(proj, compound_ids)
      if (is.null(v)) {
        cli::cli_abort(c(
          "No {.val adme_filtered} results found; {.fn rank_candidates} requires {.fn adme_filter} to have run first.",
          "i" = "Run {.fn adme_local}, then {.fn adme_filter}, before {.fn rank_candidates}."
        ))
      }
      crit_values$crit_adme_pass_frac <- v
      used <- c(used, "adme")
    }

    if ("centrality" %in% criteria) {
      v <- .rank_criterion_centrality(proj, cond, compound_ids, edges_cond, roll_up)
      if (is.null(v)) {
        cli::cli_abort(c(
          "No {.val network_centrality} results for condition {.val {cond}}; {.fn rank_candidates} requires it.",
          "i" = "Run {.fn network_centrality} first."
        ))
      }
      crit_values$crit_centrality <- v
      used <- c(used, "centrality")
    }

    if ("hub_penalty" %in% criteria) {
      v <- .rank_criterion_hub_penalty(proj, cond, edges_cond, roll_up)
      if (is.null(v)) {
        skipped <- c(skipped, "hub_penalty")
      } else {
        crit_values$crit_hub_penalty <- v
        used <- c(used, "hub_penalty")
      }
    }

    resolved_disease <- NULL
    if (any(c("proximity", "synergy") %in% criteria)) {
      resolved_disease <- .rank_resolve_disease(proj, cond, disease)
    }

    if ("proximity" %in% criteria) {
      v <- if (is.null(resolved_disease)) NULL else .rank_criterion_proximity(proj, cond, resolved_disease, compound_ids)
      if (is.null(v)) {
        skipped <- c(skipped, "proximity")
      } else {
        crit_values$crit_proximity_z <- v
        used <- c(used, "proximity")
      }
    }

    if ("synergy" %in% criteria) {
      syn <- if (is.null(resolved_disease)) NULL else .rank_criterion_synergy(proj, cond, resolved_disease, compound_ids)
      if (is.null(syn)) {
        skipped <- c(skipped, "synergy")
      } else {
        crit_values$crit_synergy_best <- syn$best
        n_p2_vec <- syn$n_p2
        used <- c(used, "synergy")
      }
    }

    if ("module_robustness" %in% criteria) {
      v <- .rank_criterion_module_robustness(proj, cond, compound_ids)
      if (is.null(v)) {
        skipped <- c(skipped, "module_robustness")
      } else {
        crit_values$crit_module_r_index <- v
        used <- c(used, "module_robustness")
      }
    }

    if (length(crit_values) == 0) {
      cli::cli_abort("Condition {.val {cond}}: no usable criteria after auto-detection.")
    }

    mat_raw <- vapply(crit_values, function(v) unname(v[compound_ids]), numeric(length(compound_ids)))
    if (is.null(dim(mat_raw))) mat_raw <- matrix(mat_raw, nrow = length(compound_ids), dimnames = list(NULL, names(crit_values)))
    rownames(mat_raw) <- compound_ids

    dirs <- .rank_direction[colnames(mat_raw)]
    mat_dir <- sweep(mat_raw, 2, dirs, `*`)

    n_criteria_used <- rowSums(!is.na(mat_raw))

    rmat <- matrix(NA_real_, nrow(mat_dir), ncol(mat_dir), dimnames = dimnames(mat_dir))
    rank_cols <- vector("list", ncol(mat_dir))
    names(rank_cols) <- colnames(mat_dir)
    for (j in seq_len(ncol(mat_dir))) {
      ranks <- rank(-mat_dir[, j], ties.method = "average", na.last = "keep")
      ## Missing evidence always gets the worst rank, not an average of
      ## trailing positions: multiple missing values must not improve each
      ## other's ranks or acquire arbitrary ordering from compound IDs.
      ranks[is.na(ranks)] <- length(compound_ids)
      rank_cols[[j]] <- ranks
      rmat[, j] <- ranks / length(compound_ids)
    }

    rra <- RobustRankAggreg::aggregateRanks(
      rmat = rmat, method = "RRA"
    )
    rra_score <- unname(stats::setNames(rra$Score, rra$Name)[compound_ids])
    ## RRA scores saturate at 1 for every compound that is not significantly
    ## better than random, so many candidates can tie exactly. Ties are broken
    ## by the mean of their per-criterion ranks (still "min" for a tie on both).
    rra_rank <- .rank_rra_rank(rra_score, rowMeans(rmat))

    mat_pareto <- mat_dir
    mat_pareto[is.na(mat_pareto)] <- -Inf
    pareto_front <- .rank_pareto_front(mat_pareto)

    out_df <- data.frame(condition = cond, compound_id = compound_ids, stringsAsFactors = FALSE)
    for (col in colnames(mat_raw)) out_df[[col]] <- mat_raw[, col]
    for (col in names(rank_cols)) out_df[[paste0("rank_", sub("^crit_", "", col))]] <- rank_cols[[col]]
    out_df$n_p2_partners <- if (is.null(n_p2_vec)) NA_integer_ else unname(n_p2_vec[compound_ids])
    out_df$n_criteria_used <- n_criteria_used
    out_df$rra_score <- rra_score
    out_df$rra_rank <- rra_rank
    out_df$pareto_front <- pareto_front
    rownames(out_df) <- NULL

    rows[[cond]] <- out_df

    proj <- .log_append(
      proj, step = "rank_candidates", id = NA_character_,
      message = paste0(
        "condition '", cond, "': ranked ", length(compound_ids), " compound(s) using criteria [",
        paste(used, collapse = ", "), "]",
        if (length(skipped) > 0) paste0("; skipped (no data): [", paste(skipped, collapse = ", "), "]") else ""
      )
    )

    if (export != "none") {
      top_ids <- out_df$compound_id[order(out_df$rra_rank)][seq_len(min(top_n, nrow(out_df)))]
      path <- .rank_export_structures(proj, top_ids, cond, top_n, export)
      proj <- .log_append(
        proj, step = "rank_candidates", id = NA_character_,
        message = paste0("condition '", cond, "': wrote top ", length(top_ids), " structure(s) to '", path, "'")
      )
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(
    proj, "rank_candidates", result, "condition",
    touched_keys = data.frame(condition = conditions, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "rank_candidates") <- result
  .write_results_csv(proj, "rank_candidates", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
#' Rank by RRA score, breaking ties with the mean per-criterion rank
#'
#' @param rra_score Numeric vector, lower = better.
#' @param mean_rank Numeric vector of the same length, lower = better.
#' @return Integer ranks; identical on both keys share the smallest rank.
#' @keywords internal
.rank_rra_rank <- function(rra_score, mean_rank) {
  key <- paste(signif(rra_score, 12), signif(mean_rank, 12))
  ord <- order(rra_score, mean_rank)
  match(key, key[ord])
}

.rank_direction <- c(
  crit_adme_pass_frac = 1, crit_centrality = 1, crit_hub_penalty = 1,
  crit_proximity_z = -1, crit_synergy_best = 1, crit_module_r_index = 1
)

#' ADME pass fraction per compound from `adme_filtered` (long format)
#' @return Named numeric vector over `compound_ids`, or `NULL` if
#'   `adme_filtered` is absent.
#' @keywords internal
.rank_criterion_adme <- function(proj, compound_ids) {
  adme <- patliRResults(proj, "adme_filtered")
  if (is.null(adme) || nrow(adme) == 0) return(NULL)
  agg <- stats::aggregate(pass ~ compound_id, data = adme, FUN = function(x) mean(as.logical(x)))
  lookup <- stats::setNames(agg$pass, agg$compound_id)
  stats::setNames(unname(lookup[compound_ids]), compound_ids)
}

#' Target-level composite centrality (mean of available `*_norm` columns)
#' rolled up target->compound
#' @return Named numeric vector over `compound_ids`, or `NULL` if
#'   `network_centrality` has no rows for `cond`. Aborts if the table exists
#'   but has no normalised column (i.e. was run with `normalize = FALSE`).
#' @keywords internal
.rank_criterion_centrality <- function(proj, cond, compound_ids, edges_cond, roll_up) {
  cent <- patliRResults(proj, "network_centrality")
  if (is.null(cent) || nrow(cent) == 0) return(NULL)
  sub <- cent[cent$condition == cond & cent$node_type == "target", , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)

  norm_cols <- intersect(c("degree_norm", "betweenness_norm", "hub_score_component"), names(sub))
  if (length(norm_cols) == 0 || all(is.na(sub[norm_cols]))) {
    cli::cli_abort(c(
      "{.val network_centrality} for condition {.val {cond}} has no normalised column to build a composite from.",
      "i" = "Re-run {.fn network_centrality} with {.code normalize = TRUE}."
    ))
  }
  composite <- rowMeans(as.matrix(sub[, norm_cols, drop = FALSE]), na.rm = TRUE)
  composite[is.nan(composite)] <- NA_real_
  target_values <- stats::setNames(composite, sub$node_id)
  .rank_rollup_target_to_compound(target_values, edges_cond, compound_ids, roll_up)
}

#' `network_hub_penalty`'s `score_adjusted` rolled up target->compound
#' @return Named numeric vector over `compound_ids`, or `NULL` if
#'   `network_hub_penalty` has no rows for `cond`.
#' @keywords internal
.rank_criterion_hub_penalty <- function(proj, cond, edges_cond, roll_up) {
  hp <- patliRResults(proj, "network_hub_penalty")
  if (is.null(hp) || nrow(hp) == 0) return(NULL)
  sub <- hp[hp$condition == cond, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  target_values <- stats::setNames(sub$score_adjusted, sub$uniprot_id)
  compound_ids <- unique(edges_cond$compound_id)
  .rank_rollup_target_to_compound(target_values, edges_cond, compound_ids, roll_up)
}

#' Resolve which `disease_id` the `"proximity"`/`"synergy"` criteria should
#' use for one condition
#' @return A `disease_id` string, or `NULL` if `network_proximity` has no
#'   rows for `cond`. Aborts if `disease` is `NULL` and more than one
#'   `disease_id` is present.
#' @keywords internal
.rank_resolve_disease <- function(proj, cond, disease) {
  prox <- patliRResults(proj, "network_proximity")
  if (is.null(prox) || nrow(prox) == 0) return(NULL)
  sub <- prox[prox$condition == cond, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  diseases <- unique(sub$disease_id)

  if (!is.null(disease)) {
    if (!disease %in% diseases) return(NULL)
    return(disease)
  }
  if (length(diseases) > 1) {
    cli::cli_abort(c(
      "Condition {.val {cond}} has {.fn network_proximity} results for multiple diseases ({.val {diseases}}).",
      "i" = "Pass {.arg disease} to {.fn rank_candidates} to pick one."
    ))
  }
  diseases[1]
}

#' `network_proximity`'s `z_score`, direct (compound-level already)
#' @return Named numeric vector over `compound_ids`, or `NULL` if there are
#'   no rows for `(cond, disease)`.
#' @keywords internal
.rank_criterion_proximity <- function(proj, cond, disease, compound_ids) {
  prox <- patliRResults(proj, "network_proximity")
  sub <- prox[prox$condition == cond & prox$disease_id == disease, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)
  lookup <- stats::setNames(sub$z_score, sub$compound_id)
  stats::setNames(unname(lookup[compound_ids]), compound_ids)
}

#' Best (`max`) `synergy_score` among a compound's partners, plus its `P2`
#' partnership count
#' @return `list(best = <named numeric>, n_p2 = <named integer>)`, or `NULL`
#'   if there are no rows for `(cond, disease)` or every `synergy_score` is
#'   `NA`.
#' @keywords internal
.rank_criterion_synergy <- function(proj, cond, disease, compound_ids) {
  syn <- patliRResults(proj, "network_synergy")
  if (is.null(syn) || nrow(syn) == 0) return(NULL)
  sub <- syn[syn$condition == cond & syn$disease_id == disease, , drop = FALSE]
  if (nrow(sub) == 0) return(NULL)

  best <- stats::setNames(rep(NA_real_, length(compound_ids)), compound_ids)
  n_p2 <- stats::setNames(rep(0L, length(compound_ids)), compound_ids)
  for (cid in compound_ids) {
    as_a <- sub[sub$compound_a == cid, , drop = FALSE]
    as_b <- sub[sub$compound_b == cid, , drop = FALSE]
    scores <- c(as_a$synergy_score, as_b$synergy_score)
    scores <- scores[!is.na(scores)]
    if (length(scores) > 0) best[[cid]] <- max(scores)
    n_p2[[cid]] <- sum(as_a$cheng_class == "P2", na.rm = TRUE) + sum(as_b$cheng_class == "P2", na.rm = TRUE)
  }
  if (all(is.na(best))) return(NULL)
  list(best = best, n_p2 = n_p2)
}

#' A compound's module `r_index`, direct join via `network_module_membership`
#' @return Named numeric vector over `compound_ids`, or `NULL` if either
#'   table is absent/empty for `cond`, or every joined `r_index` is `NA`.
#' @keywords internal
.rank_criterion_module_robustness <- function(proj, cond, compound_ids) {
  membership <- patliRResults(proj, "network_module_membership")
  robustness <- patliRResults(proj, "network_module_robustness")
  if (is.null(membership) || is.null(robustness) || nrow(membership) == 0 || nrow(robustness) == 0) return(NULL)

  mem_sub <- membership[membership$condition == cond & membership$node_type == "compound", , drop = FALSE]
  rob_sub <- robustness[robustness$condition == cond, , drop = FALSE]
  if (nrow(mem_sub) == 0 || nrow(rob_sub) == 0) return(NULL)

  r_lookup <- stats::setNames(rob_sub$r_index, rob_sub$module_id)
  mod_of <- stats::setNames(mem_sub$module_id, mem_sub$node_id)

  out <- stats::setNames(rep(NA_real_, length(compound_ids)), compound_ids)
  for (cid in compound_ids) {
    mid <- unname(mod_of[cid])
    if (!is.na(mid) && mid %in% names(r_lookup)) out[[cid]] <- unname(r_lookup[[mid]])
  }
  if (all(is.na(out))) return(NULL)
  out
}

#' Roll up a target-level value to compound level via the compound's own
#' predicted targets
#'
#' @description
#' For each compound, looks up its targets and their [network_build()] edge
#' weights (import probability) in `edges_cond`, drops targets with no
#' value in `target_values`, and combines the rest per `roll_up`. A
#' compound with no targets carrying a value gets `NA` (worst-rank
#' imputation happens later, at the ranking step, not here).
#'
#' @param target_values Named numeric vector, names = `uniprot_id`.
#' @param edges_cond `data.frame(compound_id, uniprot_id, weight)` for one
#'   condition.
#' @return Named numeric vector over `compound_ids`.
#' @keywords internal
.rank_rollup_target_to_compound <- function(target_values, edges_cond, compound_ids, roll_up) {
  out <- stats::setNames(rep(NA_real_, length(compound_ids)), compound_ids)
  for (cid in compound_ids) {
    sub <- edges_cond[edges_cond$compound_id == cid, , drop = FALSE]
    if (nrow(sub) == 0) next
    v <- unname(target_values[sub$uniprot_id])
    w <- sub$weight
    keep <- !is.na(v)
    if (!any(keep)) next
    v <- v[keep]; w <- w[keep]
    out[[cid]] <- switch(roll_up,
      weighted_mean = if (sum(w) > 0) sum(v * w) / sum(w) else mean(v),
      mean = mean(v),
      max = max(v)
    )
  }
  out
}

#' From-scratch non-dominated sort (NSGA-II-style fronts)
#'
#' @description
#' `mat` must already be direction-normalised (every column higher = better)
#' with no `NA` (the caller substitutes `-Inf` for a missing value -- same
#' "missing = worst" convention as the RRA imputation). Front `1` = the
#' non-dominated set (no other row is at-least-as-good on every column and
#' strictly better on at least one); front `2` = non-dominated after
#' removing front `1`; etc. `O(n^2)` per front -- fine at the candidate-list
#' scale `rank_candidates()` targets (tens to a few hundred rows); not
#' meant for anything larger.
#'
#' @param mat Numeric matrix, rows = compounds, columns = criteria.
#' @return Integer vector, length `nrow(mat)`, the front/tier of each row
#'   (`1` = non-dominated).
#' @keywords internal
.rank_pareto_front <- function(mat) {
  n <- nrow(mat)
  front <- integer(n)
  remaining <- seq_len(n)
  tier <- 1L
  dominates <- function(i, j) all(mat[i, ] >= mat[j, ]) && any(mat[i, ] > mat[j, ])

  while (length(remaining) > 0) {
    non_dominated <- vapply(remaining, function(i) {
      !any(vapply(remaining, function(j) i != j && dominates(j, i), logical(1)))
    }, logical(1))
    front[remaining[non_dominated]] <- tier
    remaining <- remaining[!non_dominated]
    tier <- tier + 1L
  }
  front
}

#' Write the top-N ranked compounds' structures to an `.sdf`/`.smi` file
#'
#' @description
#' Reuses the same `.parse_smiles_safe()` identity chain
#' [compounds_similarity()] already relies on -- no new dependency, no
#' re-lookup of chemical identity. `rcdk::write.molecules()` always emits an
#' SDF mol block regardless of the output filename's extension (verified
#' directly against the installed `rcdk` before writing this) -- so
#' `export = "smi"` is built by hand from `rcdk::get.smiles()` instead of
#' delegating to `write.molecules()`.
#'
#' @return The path written to, invisibly.
#' @keywords internal
.rank_export_structures <- function(proj, compound_ids, condition, top_n, export) {
  cmp <- compounds(proj)
  cmp <- cmp[match(compound_ids, cmp$id), , drop = FALSE]
  mols <- .parse_smiles_safe(cmp$smiles)
  parsed <- !vapply(mols, is.null, logical(1))
  if (!any(parsed)) {
    cli::cli_abort("None of the top {length(compound_ids)} compound(s) for condition {.val {condition}} have a parseable SMILES.")
  }
  if (any(!parsed)) {
    cli::cli_warn("{sum(!parsed)} of {length(compound_ids)} top compound(s) for condition {.val {condition}} could not be parsed and are excluded from the export.")
  }
  mols <- mols[parsed]
  ids <- cmp$id[parsed]

  results_dir <- file.path(projectDir(proj), "results")
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  ext <- if (export == "sdf") "sdf" else "smi"
  path <- file.path(results_dir, paste0("rank_candidates_top", top_n, "_", condition, ".", ext))

  if (export == "sdf") {
    rcdk::write.molecules(mols, path, together = TRUE)
  } else {
    smiles_out <- vapply(mols, rcdk::get.smiles, character(1))
    writeLines(paste(smiles_out, ids), path)
  }
  invisible(path)
}
