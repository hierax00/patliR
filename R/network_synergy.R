#' @include AllGenerics.R internal.R network_build.R network_proximity.R
NULL

## Interim design until rank_candidates() ships (see ROADMAP.md).
## `pairs = "rank_top"` ranks candidates by network_proximity()'s z_score
## before pairing (avoids the combinatorial blow-up); `pairs = "all"` skips
## ranking and scores every pair, for small conditions. No permutation test
## on the synergy score itself -- network_proximity() already carries its
## own permutation z_score per compound.

#' Candidate synergistic compound pairs (target complementarity + joint
#' disease proximity)
#'
#' @description
#' For pairs of compounds present in a condition's network, combines two
#' signals into a `synergy_score`:
#'
#' - **Complementarity** (`1 - target_jaccard`): how different the two
#'   compounds' target sets are -- a pair that hits the same targets is
#'   redundant, not synergistic.
#' - **Joint disease proximity** (`-mean(z_score_a, z_score_b)`, from
#'   already-computed [network_proximity()] results for `disease`): how
#'   much closer than random *both* compounds' targets are to the disease
#'   gene module -- a pair only scores well here if both members are
#'   independently plausible against this disease, not just one of them.
#'
#' `synergy_score = complementarity * joint_closeness`: high only when a
#' pair is both structurally complementary *and* jointly close to the
#' disease module (the pattern the network-medicine literature associates
#' with real synergy, as opposed to two compounds that are simply
#' redundant, or complementary but biologically irrelevant to the disease
#' -- see `DESIGN.md`).
#'
#' @section Requires [network_proximity()] to have already run for `disease`:
#' Unlike [network_degeneracy()] (which only needs [network_build()]), pair
#' scoring here needs a per-compound proximity `z_score`, so this errors
#' clearly if [network_proximity()] has not been run yet for the
#' `(condition, disease)` combination, rather than silently treating missing
#' scores as zero.
#'
#' @section `pairs = "rank_top"` vs. `pairs = "all"`:
#' `"rank_top"` (default) ranks compounds by their [network_proximity()]
#' `z_score` (most negative = closest to the disease module = best) and
#' only forms pairs among the top `top_n`, avoiding the combinatorial cost
#' of scoring every possible pair -- this is an **interim** stand-in for the
#' not-yet-implemented `rank_candidates()` (see `ROADMAP.md`). `"all"`
#' scores every compound pair in the condition's network
#' regardless of ranking; pairs involving a compound with no
#' [network_proximity()] score (e.g. it had no targets mappable to STRING)
#' get `NA` in the proximity-derived columns, not a fabricated score.
#'
#' @inheritParams network_build
#' @param disease The same `disease_id` used in the [network_proximity()]
#'   call whose results this function reuses.
#' @param pairs `"rank_top"` (default) or `"all"`. See the section above.
#' @param top_n Only used when `pairs = "rank_top"`: how many top-ranked
#'   compounds (by `z_score`) to form pairs among. Default `10`.
#'
#' @return The updated `proj`, with a `network_synergy` entry in
#'   [patliRResults()] (columns `condition`, `disease_id`, `compound_a`,
#'   `compound_b`, `target_jaccard`, `complementarity`, `z_score_a`,
#'   `z_score_b`, `joint_closeness`, `synergy_score`, `pairs_mode`), also
#'   written to `results/network_synergy.csv`.
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
#' disease_id <- unique(patliRResults(proj, "targets_disease")$disease_id)[1]
#' proj <- network_proximity(proj, condition = "FLO-ET", disease = disease_id)
#' proj <- network_synergy(proj, condition = "FLO-ET", disease = disease_id)
#' patliRResults(proj, "network_synergy")
#' }
#'
#' @export
network_synergy <- function(proj, condition = NULL, disease,
                             pairs = c("rank_top", "all"), top_n = 10) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(disease), length(disease) == 1, nzchar(disease))
  pairs <- match.arg(pairs)
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 2)

  conditions <- .network_resolve_conditions(proj, condition)

  proximity_all <- patliRResults(proj, "network_proximity")
  if (is.null(proximity_all) || nrow(proximity_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_proximity} entry in {.arg proj}.",
      "i" = "Run {.fn network_proximity} for {.val {disease}} first."
    ))
  }

  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    prox <- proximity_all[proximity_all$condition == cond & proximity_all$disease_id == disease, , drop = FALSE]
    if (nrow(prox) == 0) {
      cli::cli_abort(c(
        "No {.fn network_proximity} results for condition {.val {cond}} and disease {.val {disease}}.",
        "i" = "Run {.code network_proximity(proj, condition = \"{cond}\", disease = \"{disease}\")} first."
      ))
    }
    z_of <- stats::setNames(prox$z_score, prox$compound_id)

    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    target_sets <- split(ct$uniprot_id, ct$compound_id)
    all_compounds <- unique(ct$compound_id)

    candidates <- if (pairs == "rank_top") {
      ranked <- prox[order(prox$z_score, na.last = TRUE), , drop = FALSE]
      ranked <- ranked[!is.na(ranked$z_score), , drop = FALSE]
      utils::head(unique(ranked$compound_id), top_n)
    } else {
      all_compounds
    }

    if (length(candidates) < 2) {
      rows[[cond]] <- .empty_network_synergy_row()
      proj <- .log_append(
        proj, step = "network_synergy", id = NA_character_,
        message = paste0("condition '", cond, "': fewer than 2 candidate compounds (pairs = '", pairs, "'), nothing to pair")
      )
      next
    }

    combos <- utils::combn(sort(candidates), 2, simplify = FALSE)
    pair_rows <- vector("list", length(combos))
    for (i in seq_along(combos)) {
      a <- combos[[i]][1]; b <- combos[[i]][2]
      ta <- target_sets[[a]]; tb <- target_sets[[b]]
      target_jaccard <- length(intersect(ta, tb)) / length(union(ta, tb))
      complementarity <- 1 - target_jaccard

      za <- z_of[[a]]; zb <- z_of[[b]]
      za <- if (is.null(za)) NA_real_ else za
      zb <- if (is.null(zb)) NA_real_ else zb
      joint_closeness <- if (!is.na(za) && !is.na(zb)) -mean(c(za, zb)) else NA_real_
      synergy_score <- if (!is.na(joint_closeness)) complementarity * joint_closeness else NA_real_

      pair_rows[[i]] <- data.frame(
        condition = cond, disease_id = disease, compound_a = a, compound_b = b,
        target_jaccard = target_jaccard, complementarity = complementarity,
        z_score_a = za, z_score_b = zb, joint_closeness = joint_closeness,
        synergy_score = synergy_score, pairs_mode = pairs, stringsAsFactors = FALSE
      )
    }

    rows[[cond]] <- do.call(rbind, pair_rows)
    n_scored <- sum(!is.na(rows[[cond]]$synergy_score))
    proj <- .log_append(
      proj, step = "network_synergy", id = NA_character_,
      message = paste0(
        "condition '", cond, "', disease '", disease, "' (pairs = '", pairs, "'): ",
        nrow(rows[[cond]]), " pair(s) evaluated, ", n_scored, " with a full synergy_score"
      )
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_synergy", result, c("condition", "disease_id", "compound_a", "compound_b"))

  patliRResults(proj, "network_synergy") <- result
  .write_results_csv(proj, "network_synergy", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.empty_network_synergy_row <- function() {
  data.frame(
    condition = character(0), disease_id = character(0), compound_a = character(0),
    compound_b = character(0), target_jaccard = double(0), complementarity = double(0),
    z_score_a = double(0), z_score_b = double(0), joint_closeness = double(0),
    synergy_score = double(0), pairs_mode = character(0), stringsAsFactors = FALSE
  )
}
