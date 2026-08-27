#' @include AllGenerics.R internal.R refdb.R
NULL

## bias_audit() ships only the MAD-based "promiscuous compound/target"
## outlier flag, computed over the project's own reference_bioactivity.
## Categorical enrichment against a disease-category background is not
## implemented -- it needs a confirmed target -> MeSH/Disease-Ontology
## source, and this package never guesses an external API. Requesting it
## raises a clear "not implemented yet" error. See ROADMAP.md.

#' Audit database bias -- MAD-based "promiscuous" compound/target outlier
#' flag
#'
#' @description
#' Flags compounds and targets whose degree in the project's own reference
#' database (`reference_bioactivity`, from [refdb_build()]) is a statistical
#' outlier on the high side -- i.e. compounds tested against an unusually
#' large number of targets, or targets hit by an unusually large number of
#' compounds, relative to every other compound/target actually present in
#' `reference_bioactivity`. This is the mechanism behind the "homogeneity"
#' finding in Zhang, Q. (2026), *Frontiers in Pharmacology* 17:1748478 (the
#' same hub molecules/targets keep showing up across unrelated network
#' pharmacology studies) and quantified in Cell Genomics 2023 (PMC10363916):
#' in STRING, a node's degree alone predicts whether it is a known drug
#' target with AUC = 77.6% -- curation bias, not necessarily real biological
#' relevance.
#'
#' This function **never removes anything** -- it only categorizes. Use
#' [bias_reweight()] afterwards to get a log-ratio-adjusted score (never
#' replacing the raw one) for entries flagged `"promiscuo"`.
#'
#' @section Outlier rule:
#' Uses the median-absolute-deviation method of Leys, C. et al. (2013),
#' *Journal of Experimental Social Psychology* 49(4), 764-766 (`stats::mad()`
#' already applies the constant that makes it a consistent estimator of the
#' standard deviation under normality). Only the high side is flagged --
#' unlike a generic two-sided outlier test, an unusually *low* degree is not
#' a promiscuity/bias concern, just a rarely-touched compound/target.
#' Compounds and targets are each compared only against their own kind (a
#' compound's degree is "how many distinct targets it was tested against";
#' a target's degree is "how many distinct compounds tested against it" --
#' the two are not on the same scale).
#'
#' @section What this does not do yet:
#' `categories`/categorical enrichment against a disease-term background is
#' not implemented -- see `ROADMAP.md`. Passing a non-`NULL` `categories`
#' raises an informative error rather than pretending to succeed.
#'
#' @inheritParams compounds
#' @param check_homogeneity Logical, default `TRUE`. Must be `TRUE` in this
#'   version -- `categories` (categorical enrichment) is not implemented yet,
#'   so `FALSE` here would leave nothing to compute.
#' @param mad_threshold Single positive number, default `2.5` (the threshold
#'   used in most applications of Leys et al. 2013). A compound/target with
#'   `mad_score > mad_threshold` is flagged `"promiscuo"`.
#' @param categories `NULL` (default, the only implemented option). Any
#'   other value raises an informative "not implemented yet" error -- see
#'   the Description.
#'
#' @return The updated `proj`, with a `bias_homogeneity` entry in
#'   [patliRResults()] (columns `id`, `tipo` (`"compound"` or `"target"`),
#'   `frecuencia_global_refdb` (number of distinct targets/compounds this
#'   id is associated with across the whole `reference_bioactivity` table),
#'   `mad_score`, `categoria` (`"promiscuo"` or `"normal"`)), also written
#'   to `results/bias_homogeneity.csv`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- refdb_build(proj, sources = c("pubchem", "chembl")) # needs internet
#' proj <- bias_audit(proj)
#' patliRResults(proj, "bias_homogeneity")
#' }
#'
#' @seealso [bias_reweight()], [bias_report()].
#' @export
bias_audit <- function(proj, check_homogeneity = TRUE, mad_threshold = 2.5, categories = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(mad_threshold), length(mad_threshold) == 1, mad_threshold > 0)

  if (!is.null(categories)) {
    cli::cli_abort(c(
      "Categorical enrichment ({.arg categories}) is not implemented yet in this version of patliR.",
      "i" = "See {.file ROADMAP.md} -- it needs a confirmed MeSH/DO category source before coding it for real.",
      "i" = "{.fn bias_audit} with the default {.code categories = NULL} runs the homogeneity/outlier check only."
    ))
  }
  if (!isTRUE(check_homogeneity)) {
    cli::cli_abort("Nothing to audit: {.arg check_homogeneity} is {.val FALSE} and {.arg categories} is {.val NULL}.")
  }

  ref_bio <- patliRResults(proj, "reference_bioactivity")
  if (is.null(ref_bio) || nrow(ref_bio) == 0) {
    cli::cli_abort(c(
      "No {.val reference_bioactivity} entry in {.arg proj}.",
      "i" = "Run {.fn refdb_build} first -- {.fn bias_audit} needs it to compute how often each compound/target appears across the whole reference database."
    ))
  }

  compound_rows <- .bias_homogeneity_rows(ref_bio, "compound_id", "target_chembl_id", "compound", mad_threshold)
  target_rows <- .bias_homogeneity_rows(ref_bio, "target_chembl_id", "compound_id", "target", mad_threshold)
  homogeneity <- rbind(compound_rows, target_rows)
  rownames(homogeneity) <- NULL

  patliRResults(proj, "bias_homogeneity") <- homogeneity
  .write_results_csv(proj, "bias_homogeneity", homogeneity)

  n_promiscuous <- sum(homogeneity$categoria == "promiscuo")
  proj <- .log_append(
    proj, step = "bias_audit", id = NA_character_,
    message = paste0(
      "homogeneity check: ", n_promiscuous, " of ", nrow(homogeneity),
      " compound/target entries flagged 'promiscuo' (mad_threshold=", mad_threshold, ")"
    )
  )
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.bias_homogeneity_rows <- function(ref_bio, id_col, other_col, tipo, mad_threshold) {
  keep <- !is.na(ref_bio[[id_col]]) & nzchar(ref_bio[[id_col]])
  ref_bio <- ref_bio[keep, , drop = FALSE]
  if (nrow(ref_bio) == 0) return(.empty_bias_homogeneity_row())

  deg_list <- split(ref_bio[[other_col]], ref_bio[[id_col]])
  degree <- vapply(deg_list, function(x) length(unique(stats::na.omit(x))), integer(1))
  ids <- names(deg_list)

  med <- stats::median(degree)
  mad_val <- stats::mad(degree)
  mad_score <- if (mad_val > 0) (degree - med) / mad_val else rep(NA_real_, length(degree))
  categoria <- ifelse(!is.na(mad_score) & mad_score > mad_threshold, "promiscuo", "normal")

  data.frame(
    id = ids, tipo = tipo, frecuencia_global_refdb = as.integer(degree),
    mad_score = mad_score, categoria = categoria,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_bias_homogeneity_row <- function() {
  data.frame(
    id = character(0), tipo = character(0), frecuencia_global_refdb = integer(0),
    mad_score = double(0), categoria = character(0), stringsAsFactors = FALSE
  )
}

#' Log-ratio-adjusted score for "promiscuous" compounds/targets
#'
#' @description
#' For every entry flagged `"promiscuo"` by [bias_audit()], computes
#' \deqn{score_{adjusted} = degree \times \log(N_{universe} / degree)}
#' the same log-ratio formula [network_hub_penalty()] uses at the
#' per-condition-network level, applied here at the whole-reference-database
#' level instead: `degree` is `frecuencia_global_refdb` (from
#' [bias_audit()]), and `N_universe` is the total number of distinct
#' targets (for a compound entry) or distinct compounds (for a target
#' entry) across the whole `reference_bioactivity` table -- recovered
#' directly from [patliRResults()]`(proj, "bias_homogeneity")` (one row per
#' distinct id of each kind), with no need to re-read `reference_bioactivity`
#' itself. Entries not flagged `"promiscuo"` keep `score_ajustado ==
#' score_crudo` -- no penalty is applied where none was warranted.
#'
#' `score_ajustado` is always a **new column** -- `score_crudo` (the raw
#' `frecuencia_global_refdb`) is never overwritten or dropped, so any
#' downstream consumer can always fall back to the unadjusted count.
#'
#' @inheritParams compounds
#'
#' @return The updated `proj`, with a `bias_reweighted` entry in
#'   [patliRResults()] (columns `id`, `tipo`, `score_crudo`,
#'   `score_ajustado`, `categoria`), also written to
#'   `results/bias_reweighted.csv`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- refdb_build(proj, sources = c("pubchem", "chembl")) # needs internet
#' proj <- bias_audit(proj)
#' proj <- bias_reweight(proj)
#' patliRResults(proj, "bias_reweighted")
#' }
#'
#' @seealso [bias_audit()], [bias_report()], [network_hub_penalty()] for the
#'   same formula at per-condition-network scope.
#' @export
bias_reweight <- function(proj) {
  stopifnot(is(proj, "PatliRProject"))

  homog <- patliRResults(proj, "bias_homogeneity")
  if (is.null(homog) || nrow(homog) == 0) {
    cli::cli_abort(c(
      "No {.val bias_homogeneity} entry in {.arg proj}.",
      "i" = "Run {.fn bias_audit} first."
    ))
  }

  n_targets_total <- sum(homog$tipo == "target")
  n_compounds_total <- sum(homog$tipo == "compound")
  n_universe <- ifelse(homog$tipo == "compound", n_targets_total, n_compounds_total)
  degree <- homog$frecuencia_global_refdb

  score_adjusted <- ifelse(
    homog$categoria == "promiscuo" & degree > 0 & n_universe > 0,
    degree * log(n_universe / degree),
    degree
  )

  result <- data.frame(
    id = homog$id, tipo = homog$tipo,
    score_crudo = degree, score_ajustado = score_adjusted,
    categoria = homog$categoria,
    stringsAsFactors = FALSE
  )

  patliRResults(proj, "bias_reweighted") <- result
  .write_results_csv(proj, "bias_reweighted", result)

  n_adjusted <- sum(result$categoria == "promiscuo")
  proj <- .log_append(
    proj, step = "bias_reweight", id = NA_character_,
    message = paste0(n_adjusted, " of ", nrow(result), " entries adjusted (categoria == 'promiscuo')")
  )
  .write_log_csv(proj)
  proj
}

#' Summarize the "promiscuous" compounds/targets found by the bias audit --
#' never a removal
#'
#' @description
#' Read-only summary over [patliRResults()]`(proj, "bias_reweighted")`,
#' following the same contract as [tox_report()]: nothing is written to
#' `proj` or to disk, and a fixed disclaimer is always attached. The
#' `"promiscuo"` table is **always present** in the returned `summary` (zero
#' rows, with a warning, if [bias_audit()]/[bias_reweight()] have not run
#' yet) -- consistent with the rest of the package's "absence of a flag is
#' always visible, never just an omitted table" convention (see
#' [tox_safetyome()]).
#'
#' @inheritParams compounds
#'
#' @return A `list` with two elements:
#'   \itemize{
#'     \item `summary`: a `data.frame` (columns `id`, `tipo`, `score_crudo`,
#'       `score_ajustado`, `categoria`), restricted to `categoria ==
#'       "promiscuo"` and ordered by `score_ajustado` (ascending -- the most
#'       heavily discounted, i.e. most promiscuous, entries first).
#'     \item `note`: the fixed disclaimer string, always identical
#'       regardless of the data.
#'   }
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- refdb_build(proj, sources = c("pubchem", "chembl")) # needs internet
#' proj <- bias_audit(proj)
#' proj <- bias_reweight(proj)
#' report <- bias_report(proj)
#' report$summary
#' report$note
#' }
#'
#' @seealso [bias_audit()], [bias_reweight()].
#' @export
bias_report <- function(proj) {
  stopifnot(is(proj, "PatliRProject"))

  note <- paste(
    "A 'promiscuo' flag means this compound/target is a statistical outlier",
    "(MAD-based, Leys et al. 2013) in how often it appears across the whole",
    "reference database -- it does not mean the underlying association is",
    "wrong. score_ajustado is a log-ratio discount (same formula as",
    "network_hub_penalty()) applied only to flagged entries; it never",
    "replaces score_crudo, which is always kept alongside it."
  )

  reweighted <- patliRResults(proj, "bias_reweighted")
  if (is.null(reweighted) || nrow(reweighted) == 0) {
    cli::cli_warn(paste(
      "No {.val bias_reweighted} results found; run {.fn bias_audit} and",
      "{.fn bias_reweight} first."
    ))
    return(list(
      summary = data.frame(
        id = character(), tipo = character(), score_crudo = double(),
        score_ajustado = double(), categoria = character(), stringsAsFactors = FALSE
      ),
      note = note
    ))
  }

  promiscuous <- reweighted[reweighted$categoria == "promiscuo", , drop = FALSE]
  promiscuous <- promiscuous[order(promiscuous$score_ajustado), , drop = FALSE]
  rownames(promiscuous) <- NULL

  list(summary = promiscuous, note = note)
}
