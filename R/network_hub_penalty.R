#' @include AllGenerics.R internal.R network_build.R
NULL

#' Promiscuity-adjusted target score per condition (hub penalty)
#'
#' @description
#' For every target node in a condition's graph, computes
#' \deqn{score_{adjusted} = degree_{raw} \times \log(N_{compounds} / degree_{raw})}
#' where `degree_raw` is how many of the condition's present compounds hit
#' that target, and `N_compounds` is the total number of compounds present
#' in that condition. This is a direct, self-contained computation from
#' [network_build()]'s graph -- it does **not** require [network_centrality()]
#' to have run first (deliberately: coupling this to another function's
#' specific `measures` choice would be a fragile dependency for a two-line
#' formula that only ever needs plain degree).
#'
#' A target hit by every compound in the condition (`degree_raw ==
#' N_compounds`) gets `score_adjusted = 0` -- maximally promiscuous, fully
#' discounted. A target hit by only one compound out of many keeps close to
#' its raw degree (`log(N_compounds / 1)` is large) -- same intuition as
#' TF-IDF, and the same log-ratio [bias_reweight()] will reuse later at the
#' reference-database level (`patliR_manual.md`, section 7) instead of at
#' this per-condition-network level.
#'
#' @inheritParams network_build
#'
#' @return The updated `proj`, with a `network_hub_penalty` entry in
#'   [patliRResults()] (columns `condition`, `uniprot_id`, `degree_raw`,
#'   `n_compounds_total`, `score_adjusted`), also written to
#'   `results/network_hub_penalty.csv`.
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
#' proj <- network_hub_penalty(proj)
#' patliRResults(proj, "network_hub_penalty")
#' }
#'
#' @export
network_hub_penalty <- function(proj, condition = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  conditions <- .network_resolve_conditions(proj, condition)

  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    g <- .network_graph(proj, cond)
    is_target <- igraph::V(g)$type
    target_ids <- igraph::V(g)$name[is_target]
    n_compounds_total <- sum(!is_target)

    if (length(target_ids) == 0) {
      rows[[cond]] <- .empty_network_hub_penalty_row()
      next
    }

    degree_raw <- igraph::degree(g, v = igraph::V(g)[is_target])
    score_adjusted <- ifelse(
      degree_raw > 0 & n_compounds_total > 0,
      degree_raw * log(n_compounds_total / degree_raw),
      NA_real_
    )

    rows[[cond]] <- data.frame(
      condition = cond, uniprot_id = target_ids, degree_raw = degree_raw,
      n_compounds_total = n_compounds_total, score_adjusted = score_adjusted,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_hub_penalty", result, "condition")

  patliRResults(proj, "network_hub_penalty") <- result
  .write_results_csv(proj, "network_hub_penalty", result)
  proj <- .log_append(
    proj, step = "network_hub_penalty", id = NA_character_,
    message = paste0("condition '", conditions, "': hub penalty computed for ",
                      vapply(conditions, function(cond) sum(result$condition == cond), integer(1)), " targets")
  )
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.empty_network_hub_penalty_row <- function() {
  data.frame(condition = character(0), uniprot_id = character(0), degree_raw = integer(0),
             n_compounds_total = integer(0), score_adjusted = double(0), stringsAsFactors = FALSE)
}
