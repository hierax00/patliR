#' @include AllGenerics.R internal.R network_build.R network_proximity.R
NULL

## Interim design until rank_candidates() ships (see ROADMAP.md).
## `pairs = "rank_top"` ranks candidates by network_proximity()'s z_score
## before pairing (avoids the combinatorial blow-up); `pairs = "all"` skips
## ranking and scores every pair, for small conditions.
##
## `separation = "network"` computes the real Menche (2015) topological
## separation s_AB on the STRING interactome and reports Cheng et al.
## (2019)'s P1..P6 classification. `separation = "jaccard"` keeps the cheap
## set-overlap proxy for users without STRINGdb.
##
## No permutation null is fitted for s_AB itself: Cheng et al. (2019)
## explicitly reject a z-score for the drug-drug relationship -- "each drug
## has only a small number of experimentally reported targets (on average
## 3) ... therefore, the randomization procedure is not producing a
## Gaussian distribution". network_proximity() already carries its own
## degree-preserving permutation z per compound, which is what the
## proximal_a / proximal_b predicates consume.

#' Cheng et al. (2019) drug-pair classification (Complementary Exposure)
#'
#' @description
#' For pairs of compounds present in a condition's network, computes the
#' topological network separation `s_AB` of Menche et al. (2015) between the
#' two compounds' STRING-mapped target sets and combines it with each
#' compound's individual disease proximity ([network_proximity()]) to place
#' the pair in one of Cheng, Kovacs & Barabasi (2019)'s six classes:
#'
#' | class | `separated` | `proximal_a` | `proximal_b` | meaning |
#' |---|---|---|---|---|
#' | `P1` Overlapping Exposure | FALSE | TRUE | TRUE | both hit D, same neighbourhood |
#' | `P2` **Complementary Exposure** | TRUE | TRUE | TRUE | both hit D, separate neighbourhoods |
#' | `P3` Indirect Exposure | FALSE | exactly one | | overlapping pair, one hits D |
#' | `P4` Single Exposure | TRUE | exactly one | | separated pair, one hits D |
#' | `P5` Non-Exposure | FALSE | FALSE | FALSE | overlapping pair, neither hits D |
#' | `P6` Independent Action | TRUE | FALSE | FALSE | everything separated |
#'
#' Cheng et al.'s finding is that **only `P2` correlates with therapeutic
#' efficacy**. `complementary_exposure` is the logical `cheng_class == "P2"`.
#'
#' Two classifications are returned, differing only in the proximity
#' predicate:
#'
#' - `cheng_class` -- **FDR-gated** (patliR's modification): `proximal_x =
#'   z_x < 0 & p_adjusted_x < alpha`. This is the stricter rule, and the one
#'   `complementary_exposure` and `synergy_score` are built on.
#' - `cheng_class_sign` -- **sign-only**, the rule of the paper: a compound
#'   is proximal when its proximity `z < 0`. `P1` = `s_AB < 0` and both `z`
#'   negative; `P2` = `s_AB >= 0` and both; `P3` / `P4` = `s_AB < 0` /
#'   `s_AB >= 0` and exactly one; `P5` / `P6` = `s_AB < 0` / `s_AB >= 0` and
#'   neither.
#'
#' The two differ whenever a compound has `z < 0` but its BH-adjusted p does
#' not clear `alpha`; report which rule a statement rests on.
#'
#' @section The separation statistic `s_AB` (Menche et al. 2015):
#' On the STRING largest connected component `G` with unweighted
#' shortest-path hop count `d(u,v)`:
#'
#' \deqn{s_{AB} = \langle d_{AB} \rangle - \tfrac{1}{2}(\langle d_{AA}
#'   \rangle + \langle d_{BB} \rangle)}
#'
#' where the within-set distance uses the "closest" convention
#' \eqn{\langle d_{AA} \rangle = |A|^{-1} \sum_{a \in A}
#'   \min_{a' \in A, a' \ne a} d(a, a')} and the between-set distance is the
#' pooled, symmetric form
#' \eqn{\langle d_{AB} \rangle = (|A| + |B|)^{-1} (\sum_{a} \min_{b} d(a,b)
#'   + \sum_{b} \min_{a} d(b,a))}. `s_AB < 0` means the two target modules
#' overlap topologically; `s_AB >= 0` (`separated = TRUE`) means they are
#' separated -- Menche et al.'s own sign convention, verbatim.
#'
#' No null model is fitted for `s_AB` (see the file header for Cheng et
#' al.'s small-target-set argument); it is used raw.
#'
#' @section Singletons (fewer than 2 mapped targets):
#' The within-set distance `<d_AA>` needs at least two members. When a
#' compound's STRING-mapped target set has fewer than 2, this function
#' still computes `d_aa = 0`, `s_ab` and `separated` (so the numbers stay
#' inspectable) and flags the row with `singleton_a` / `singleton_b`. That
#' `0` is **not neutral**: it makes `s_AB` systematically more positive, and
#' two identical single-target sets get `s_AB = 0`, i.e. `separated = TRUE`
#' and possibly `P2` for what is one and the same protein. Menche et al.'s
#' reference `separation.py` returns `nan` instead.
#'
#' Therefore, by default (`singleton = "na"`), the separation-based classes
#' `cheng_class`, `cheng_class_sign` and `complementary_exposure` are `NA`
#' whenever `singleton_a` or `singleton_b` is `TRUE`. `singleton = "zero"`
#' restores the old convention (classes computed from the `d_AA = 0`
#' separation); the flags are kept either way, the choice is recorded in
#' `singleton_policy`, and [plot_synergy()] never counts a singleton pair
#' as `P2`. `synergy_score` is `NA` for a singleton-flagged pair under both
#' policies. With predicted plant-metabolite targets a mapped set of size 1
#' is not rare.
#'
#' @section `alpha` is patliR's tightening, not Cheng's:
#' Cheng et al. (2019) use the drug-disease proximity z-score and never
#' state a p-value cutoff for the "overlaps the disease module" predicate.
#' `proximal_x = (z_x < 0 AND p_adjusted_x < alpha)` is a *stricter* patliR
#' choice; `cheng_class_sign` gives the paper's sign-only classification
#' alongside it. `p_adjusted` is [network_proximity()]'s Benjamini-Hochberg
#' value, computed across **all (condition, compound) tests of one
#' [network_proximity()] call**, not per pair and not per condition.
#'
#' The BH-adjusted p of the `i`-th smallest of `m` p-values is
#' `q_(i) = min_{j >= i} m * p_(j) / j`, so it depends on the rank as well
#' as on `m`: with 30 tests, `n_random = 100` and every empirical p at its
#' floor `1/101`, every `q` is `1/101 < 0.05`. The only situation in which
#' the gate is arithmetically unreachable is `1 / (n_random + 1) >= alpha`
#' (a BH-adjusted p is never below its raw p); that case warns. Separately,
#' when a condition has compounds with `z < 0` but none with
#' `p_adjusted < alpha` (an observed outcome, not a bound), this is logged
#' and reported as a message, since the FDR-gated classes can then only be
#' `P5`/`P6`. If a merged `network_proximity` table holds rows from calls
#' with different BH family sizes, a warning notes that their `p_adjusted`
#' values are not one coherent FDR family.
#'
#' @section Distance direction, interactome and comparability with Cheng et al.:
#' [network_proximity()]'s "closest" distance averages, over the
#' **compound's** mapped targets, the hop count to the nearest disease gene
#' (drug target -> nearest disease gene), as in the reference proximity
#' toolbox of Guney et al. (2016) / Cheng et al. Cheng et al.'s Eq. 1 as
#' printed averages the other way (disease gene -> nearest drug target).
#' The two are not symmetric, so `z` values are not interchangeable with a
#' disease-side implementation.
#'
#' The interactome is STRING at `score_threshold` (default `400`, "medium
#' confidence"), a **functional-association** network that includes
#' text-mining, co-expression and database-transferred edges -- not Cheng
#' et al.'s interactome of experimentally supported physical protein
#' interactions. Absolute `s_AB` values and class assignments are therefore
#' not directly comparable with the paper's; use them to compare pairs
#' within one patliR run. `score_threshold` is recorded in the output.
#'
#' @section Ranking quantities are ad hoc, not efficacy measures:
#' `complementarity` (`1 - target_jaccard`), `joint_closeness`
#' (`-max(z_a, z_b)`) and `synergy_score` (`complementarity *
#' joint_closeness`, gated on `complementary_exposure` in network mode and
#' on `both_proximal` in Jaccard mode) are patliR's **ranking** heuristics
#' for ordering candidate pairs. They are not validated measures of
#' efficacy, of pharmacological synergy (Bliss, Loewe, ...), or of
#' anything Cheng et al. tested.
#'
#' @section `separation = "jaccard"`:
#' Keeps the cheap set-overlap proxy for users without `STRINGdb`.
#' `d_aa` / `d_bb` / `d_ab` / `s_ab` / `separated` / `cheng_class` /
#' `cheng_class_sign` / `complementary_exposure` / `singleton_a` /
#' `singleton_b` are all `NA`, `separation_method = "jaccard"` records it,
#' and `synergy_score` is gated on `both_proximal` rather than
#' `complementary_exposure`. (Earlier versions defined `both_proximal` as
#' the raw-sign quantity now called `both_negative_z`; Jaccard-mode
#' `synergy_score` is now FDR-gated like the network mode.)
#' `target_jaccard` / `complementarity` are
#' computed on the **raw UniProt** target sets in *both* modes (never the
#' STRING-mapped subset), so they are comparable across modes;
#' `n_targets_a_mapped` / `n_targets_b_mapped` carry the mapped counts
#' separately.
#'
#' @section `proximal_*`, `both_proximal`, `both_negative_z` and their NA policies:
#' `proximal_x` is three-valued: `NA` when `z_x` is missing, otherwise
#' `z_x < 0 & p_adjusted_x < alpha` (falling back per-row to `z_x < 0`
#' alone, logged, when that row's `p_adjusted` is `NA`).
#'
#' `both_proximal` uses **the same predicates as `cheng_class`**:
#' `isTRUE(proximal_a) & isTRUE(proximal_b)`. It is two-valued -- a missing
#' `proximal_x` counts as not proximal, so a pair with an unscored compound
#' is `FALSE`, never `NA` -- and `both_proximal = TRUE` implies
#' `proximal_a = proximal_b = TRUE`. `both_negative_z` is the raw-sign
#' quantity (`z_a < 0 & z_b < 0`, missing `z` -> `FALSE`), i.e. the
#' predicate of the sign-only `cheng_class_sign`; it can be `TRUE` while
#' both `proximal_*` are `FALSE`.
#'
#' `cheng_class` is `NA` whenever any of `separated`, `proximal_a`,
#' `proximal_b` is `NA`; `cheng_class_sign` whenever `separated`, `z_a` or
#' `z_b` is `NA`; both are also `NA` for singleton pairs under the default
#' `singleton = "na"`.
#'
#' @section Requires [network_proximity()] to have already run for `disease`:
#' Pair scoring needs a per-compound proximity `z_score`, so this errors
#' clearly if [network_proximity()] has not been run for the
#' `(condition, disease)` combination. When `separation = "network"` the
#' `species` / `version` / `score_threshold` recorded by
#' [network_proximity()] must match this call's arguments, or `s_AB` and
#' `z` would be on different interactomes -- a mismatch aborts.
#'
#' @section `pairs = "rank_top"` vs. `pairs = "all"`:
#' `"rank_top"` (default) ranks compounds by their [network_proximity()]
#' `z_score` and only forms pairs among the top `top_n` -- an **interim**
#' stand-in for the not-yet-implemented `rank_candidates()`. `"all"` scores
#' every compound pair; pairs involving a compound with no
#' [network_proximity()] score get `NA` in the proximity-derived columns,
#' not a fabricated score.
#'
#' @inheritParams network_build
#' @inheritParams network_proximity
#' @param disease The same `disease_id` used in the [network_proximity()]
#'   call whose results this function reuses.
#' @param pairs `"rank_top"` (default) or `"all"`. See the section above.
#' @param top_n Only used when `pairs = "rank_top"`: how many top-ranked
#'   compounds (by `z_score`) to form pairs among. Default `10`.
#' @param separation `"network"` (default) for the real Menche (2015)
#'   `s_AB` on the STRING interactome, or `"jaccard"` for the set-overlap
#'   proxy. `"network"` aborts if `STRINGdb` is not installed.
#' @param alpha Proximity significance cut for the `proximal_a` /
#'   `proximal_b` predicates (patliR's tightening of Cheng et al., not
#'   theirs). Default `0.05`.
#' @param disease_gene_source Which [network_proximity()] rows to consume:
#'   `"disease_genes"` (default, the independent gene set) or the legacy,
#'   circular `"targets_disease"` (a warning is emitted -- the Cheng
#'   classification would inherit the Phase-1 circularity).
#' @param singleton `"na"` (default) or `"zero"`. What the separation-based
#'   classes (`cheng_class`, `cheng_class_sign`, `complementary_exposure`)
#'   do for a pair in which either compound has fewer than 2 STRING-mapped
#'   targets: `"na"` leaves them `NA`; `"zero"` computes them from the
#'   `d_AA = 0` convention (the pre-audit behaviour). See the Singletons
#'   section. Ignored when `separation = "jaccard"`.
#'
#' @return The updated `proj`, with a `network_synergy` entry in
#'   [patliRResults()] (columns `condition`, `disease_id`, `compound_a`,
#'   `compound_b`, `n_targets_a`, `n_targets_b`, `n_targets_a_mapped`,
#'   `n_targets_b_mapped`, `target_jaccard`, `complementarity`,
#'   `singleton_a`, `singleton_b`, `d_aa`, `d_bb`, `d_ab`, `s_ab`,
#'   `separated`, `z_score_a`, `z_score_b`, `p_adjusted_a`, `p_adjusted_b`,
#'   `proximal_a`, `proximal_b`, `both_proximal`, `both_negative_z`,
#'   `cheng_class` (FDR-gated), `cheng_class_sign` (sign-only, the paper's
#'   rule), `complementary_exposure`, `joint_closeness`, `synergy_score`,
#'   `separation_method`, `alpha`, `singleton_policy`, `pairs_mode`,
#'   `species`, `string_version`, `score_threshold`,
#'   `disease_gene_source`), also written to
#'   `results/network_synergy.csv`.
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
#' proj <- disease_genes_fetch(proj, disease = "type 2 diabetes mellitus")
#' disease_id <- unique(patliRResults(proj, "disease_genes")$disease_id)[1]
#' proj <- network_proximity(proj, condition = "FLO-ET", disease = disease_id)
#' proj <- network_synergy(proj, condition = "FLO-ET", disease = disease_id)
#' patliRResults(proj, "network_synergy")
#' }
#'
#' @references
#' Cheng, Kovacs & Barabasi (2019), *Nat Commun* 10:1197,
#' \doi{10.1038/s41467-019-09186-x}. Menche et al. (2015), *Science*
#' 347(6224):1257601, \doi{10.1126/science.1257601}.
#'
#' @export
network_synergy <- function(proj, condition = NULL, disease,
                             pairs = c("rank_top", "all"), top_n = 10,
                             separation = c("network", "jaccard"),
                             alpha = 0.05, species = 9606, version = "12.0",
                             score_threshold = 400,
                             disease_gene_source = c("disease_genes", "targets_disease"),
                             singleton = c("na", "zero")) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(disease), length(disease) == 1, nzchar(disease))
  pairs <- match.arg(pairs)
  separation <- match.arg(separation)
  disease_gene_source <- match.arg(disease_gene_source)
  singleton <- match.arg(singleton)
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 2)
  stopifnot(is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1)

  if (separation == "network" && !requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_synergy} with {.code separation = \"network\"} needs {.pkg STRINGdb}, not installed.",
      "i" = "Install it with {.code BiocManager::install(\"STRINGdb\")}, or pass {.code separation = \"jaccard\"} for the set-overlap proxy."
    ))
  }

  conditions <- .network_resolve_conditions(proj, condition)

  proximity_all <- patliRResults(proj, "network_proximity")
  if (is.null(proximity_all) || nrow(proximity_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_proximity} entry in {.arg proj}.",
      "i" = "Run {.fn network_proximity} for {.val {disease}} first."
    ))
  }

  ## Subset on disease_gene_source so a merged table holding both a
  ## "disease_genes" and a legacy "targets_disease" run for the same
  ## disease does not produce a duplicated z_of name (which setNames() +
  ## z_of[a] would silently resolve to the first match).
  ## `provenance_known` is FALSE for a pre-Phase-1 network_proximity table
  ## (no `disease_gene_source` column, or all-NA after a `.network_upsert()`
  ## back-fill). In that case we cannot subset on it, and the synergy output
  ## must record `disease_gene_source = NA` rather than fabricate a value.
  have_dgs <- "disease_gene_source" %in% names(proximity_all)
  provenance_known <- have_dgs &&
    any(!is.na(proximity_all$disease_gene_source))
  if (provenance_known) {
    present_sources <- unique(stats::na.omit(proximity_all$disease_gene_source))
    if (length(present_sources) > 1) {
      cli::cli_warn(c(
        "!" = "The {.val network_proximity} table holds rows from multiple {.field disease_gene_source} values ({.val {present_sources}}).",
        "i" = "Using {.val {disease_gene_source}}; pass {.arg disease_gene_source} to pick the other."
      ))
    }
    proximity_all <- proximity_all[
      !is.na(proximity_all$disease_gene_source) &
        proximity_all$disease_gene_source == disease_gene_source, , drop = FALSE]
  }
  effective_dgs <- if (provenance_known) disease_gene_source else NA_character_
  if (provenance_known && disease_gene_source == "targets_disease" &&
        nrow(proximity_all) > 0) {
    cli::cli_warn(c(
      "!" = "{.code disease_gene_source = \"targets_disease\"}: the Cheng classification inherits the Phase-1 circularity.",
      "i" = "{.fn targets_disease_filter} only annotates UniProt IDs already predicted as targets, so the {.field z_score} the {.field proximal_*} predicates consume is biased negative for reasons unrelated to topology."
    ))
  } else if (!provenance_known) {
    cli::cli_warn(c(
      "!" = "The {.val network_proximity} rows predate the {.field disease_gene_source} column.",
      "i" = "Cannot verify the disease module was independent of the compounds' predicted targets; the Phase-1 circularity may apply. {.field disease_gene_source} is recorded as {.val NA}."
    ))
  }

  string_db <- NULL
  lcc <- NULL
  g <- NULL
  g_names <- character(0)
  if (separation == "network") {
    string_db <- .network_stringdb(proj, species, version, score_threshold)
    lcc <- .network_string_lcc(proj, species, version, score_threshold, string_db = string_db)
    g <- lcc$graph
    g_names <- igraph::V(g)$name
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

    ## s_AB and z must be on the same interactome (Cheng et al. 2019). If
    ## network_proximity() recorded its interactome parameters, a mismatch
    ## with this call's arguments is a hard error.
    if (separation == "network") {
      if (all(c("species", "string_version", "score_threshold") %in% names(prox))) {
        ps <- unique(prox$species); pv <- unique(as.character(prox$string_version))
        pt <- unique(prox$score_threshold)
        ## `string_version` is character by contract on both sides
        ## (`.patliR_results_colclasses` pins it through the CSV round-trip),
        ## so a direct string compare is exact; `species` / `score_threshold`
        ## are genuine numbers.
        mism <- !isTRUE(all.equal(as.numeric(ps), as.numeric(species))) ||
          !identical(pv, as.character(version)) ||
          !isTRUE(all.equal(as.numeric(pt), as.numeric(score_threshold)))
        if (length(ps) > 1 || length(pv) > 1 || length(pt) > 1 || mism) {
          cli::cli_abort(c(
            "{.fn network_synergy}: the {.val network_proximity} rows for condition {.val {cond}} were computed on a different STRING interactome than this call.",
            "i" = "network_proximity: species {.val {ps}}, version {.val {pv}}, score_threshold {.val {pt}}.",
            "i" = "network_synergy: species {.val {species}}, version {.val {version}}, score_threshold {.val {score_threshold}}.",
            "i" = "Re-run {.fn network_proximity} with matching arguments, or match them here."
          ))
        }
      } else {
        cli::cli_warn(c(
          "!" = "The {.val network_proximity} rows predate the {.field species} / {.field string_version} / {.field score_threshold} provenance columns.",
          "i" = "Cannot verify {.field s_AB} and {.field z_score} are on the same STRING interactome."
        ))
      }
    }

    ## FDR-gate diagnostics. The BH-adjusted p of the i-th smallest of m
    ## p-values is q_(i) = min_{j >= i} m * p_(j) / j: it depends on the
    ## rank, so m / (n_random + 1) is NOT a floor (30 tests all at p = 1/101
    ## give q = 1/101 each). The only hard bound is q_(i) >= p_(i) >=
    ## 1 / (n_random + 1); beyond that, report what the stored p_adjusted
    ## values actually do.
    if ("n_tests_in_family" %in% names(prox)) {
      fam <- unique(stats::na.omit(prox$n_tests_in_family))
      if (length(fam) > 1) {
        cli::cli_warn(c(
          "!" = "The {.val network_proximity} rows for condition {.val {cond}} were assembled from calls with different BH family sizes ({.val {fam}}).",
          "i" = "Their {.field p_adjusted} values come from different FDR families, so the {.arg alpha} gate is not one coherent FDR threshold across the pair set."
        ))
      }
    }
    nr <- suppressWarnings(min(prox$n_random, na.rm = TRUE))
    if (is.finite(nr) && 1 / (nr + 1) >= alpha) {
      cli::cli_warn(c(
        "!" = "The proximity FDR gate cannot be cleared for condition {.val {cond}}: the empirical p has a floor of 1/({.arg n_random} + 1) = {signif(1 / (nr + 1), 3)} >= {.arg alpha} = {alpha}, and a BH-adjusted p is never below its raw p.",
        "i" = "Every compound with a {.field p_adjusted} is then not {.field proximal}, so {.field cheng_class} can only be {.val P5}/{.val P6}; {.field cheng_class_sign} (sign-only rule) is unaffected. Rerun {.fn network_proximity} with a larger {.arg n_random}."
      ))
    } else if ("p_adjusted" %in% names(prox)) {
      neg_z <- !is.na(prox$z_score) & prox$z_score < 0 & !is.na(prox$p_adjusted)
      if (any(neg_z) && !any(prox$p_adjusted[neg_z] < alpha)) {
        min_q <- min(prox$p_adjusted[neg_z])
        msg <- paste0(
          "condition '", cond, "': ", sum(neg_z), " compound(s) have z < 0 but none has BH p_adjusted < alpha = ",
          alpha, " (smallest ", signif(min_q, 3), "); the FDR-gated cheng_class can only be P5/P6 here, see cheng_class_sign for the sign-only rule"
        )
        cli::cli_inform(c("i" = msg))
        proj <- .log_append(proj, step = "network_synergy", id = NA_character_, message = msg)
      }
    }

    z_of <- stats::setNames(prox$z_score, prox$compound_id)
    padj_of <- if ("p_adjusted" %in% names(prox)) {
      stats::setNames(prox$p_adjusted, prox$compound_id)
    } else {
      NULL
    }

    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    target_sets <- lapply(split(ct$uniprot_id, ct$compound_id), unique)
    all_compounds <- unique(ct$compound_id)

    candidates <- if (pairs == "rank_top") {
      ranked <- prox[order(prox$z_score, na.last = TRUE), , drop = FALSE]
      ranked <- ranked[!is.na(ranked$z_score), , drop = FALSE]
      utils::head(unique(ranked$compound_id), top_n)
    } else {
      all_compounds
    }
    candidates <- intersect(candidates, all_compounds)

    if (length(candidates) < 2) {
      rows[[cond]] <- .empty_network_synergy_row()
      proj <- .log_append(
        proj, step = "network_synergy", id = NA_character_,
        message = paste0("condition '", cond, "': fewer than 2 candidate compounds (pairs = '", pairs, "'), nothing to pair")
      )
      next
    }

    ## network mode: one map() for every candidate target, then ONE
    ## distances() over U = the union of every candidate's mapped, in-graph
    ## targets. Every pair's <d_AA>/<d_BB>/<d_AB> is then a submatrix
    ## min/mean -- not a fresh BFS per pair (pre-review 12-U7).
    mapped_sets <- NULL
    D <- NULL
    d_aa_of <- NULL
    if (separation == "network") {
      cand_uni <- unique(unlist(target_sets[candidates], use.names = FALSE))
      map_df <- string_db$map(
        data.frame(uniprot_id = cand_uni, stringsAsFactors = FALSE),
        "uniprot_id", removeUnmappedRows = FALSE, quiet = TRUE
      )
      unmapped <- map_df$uniprot_id[is.na(map_df$STRING_id)]
      if (length(unmapped) > 0) {
        proj <- .log_append(
          proj, step = "network_synergy", id = NA_character_,
          message = paste0(
            "network_synergy_unmapped: ", length(unmapped), " of ", nrow(map_df),
            " UniProt ID(s) had no STRING_id for species ", species, " and were excluded (e.g. ",
            paste(utils::head(unmapped, 10), collapse = ", "), ")"
          )
        )
      }
      uni_to_string <- stats::setNames(map_df$STRING_id, map_df$uniprot_id)
      mapped_sets <- lapply(target_sets[candidates], function(u) {
        s <- unique(stats::na.omit(uni_to_string[unique(u)]))
        s[s %in% g_names]
      })
      names(mapped_sets) <- candidates
      U <- unique(unlist(mapped_sets, use.names = FALSE))
      D <- if (length(U) >= 1) {
        igraph::distances(g, v = U, to = U, weights = NA)
      } else {
        NULL
      }
      d_aa_of <- lapply(mapped_sets, function(s) .network_set_distance(g, s, D))
    }

    combos <- utils::combn(sort(candidates), 2, simplify = FALSE)
    pair_rows <- vector("list", length(combos))
    p_fallback_used <- FALSE
    disconnected_pairs <- character(0)
    n_singleton_na <- 0L

    for (i in seq_along(combos)) {
      a <- combos[[i]][1]; b <- combos[[i]][2]
      ta <- unique(target_sets[[a]]); tb <- unique(target_sets[[b]])
      n_targets_a <- length(ta); n_targets_b <- length(tb)
      inter <- length(intersect(ta, tb)); uni <- length(union(ta, tb))
      target_jaccard <- if (uni > 0) inter / uni else NA_real_
      complementarity <- 1 - target_jaccard

      za <- unname(z_of[a]); zb <- unname(z_of[b])
      pa <- if (is.null(padj_of)) NA_real_ else unname(padj_of[a])
      pb <- if (is.null(padj_of)) NA_real_ else unname(padj_of[b])
      if ((!is.na(za) && is.na(pa)) || (!is.na(zb) && is.na(pb))) p_fallback_used <- TRUE

      ## FDR-gated predicates (cheng_class) and raw-sign predicates
      ## (cheng_class_sign, the paper's rule). both_proximal is built from
      ## the former, both_negative_z from the latter; both map a missing
      ## predicate to FALSE.
      proximal_a <- .synergy_proximal(za, pa, alpha)
      proximal_b <- .synergy_proximal(zb, pb, alpha)
      negative_a <- if (is.na(za)) NA else za < 0
      negative_b <- if (is.na(zb)) NA else zb < 0
      both_proximal <- isTRUE(proximal_a) && isTRUE(proximal_b)
      both_negative_z <- isTRUE(negative_a) && isTRUE(negative_b)
      joint_closeness <- if (!is.na(za) && !is.na(zb)) -max(za, zb) else NA_real_

      if (separation == "network") {
        A <- mapped_sets[[a]]; B <- mapped_sets[[b]]
        n_a_mapped <- length(A); n_b_mapped <- length(B)
        singleton_a <- n_a_mapped < 2L
        singleton_b <- n_b_mapped < 2L
        singleton_pair <- singleton_a || singleton_b
        d_aa <- d_aa_of[[a]]; d_bb <- d_aa_of[[b]]
        if (n_a_mapped == 0L || n_b_mapped == 0L) {
          d_ab <- NA_real_; s_ab <- NA_real_
        } else {
          d_ab <- .network_between_set_distance(g, A, B, D)
          s_ab <- if (anyNA(c(d_ab, d_aa, d_bb))) NA_real_ else d_ab - (d_aa + d_bb) / 2
          if (is.na(s_ab)) disconnected_pairs <- c(disconnected_pairs, paste0(a, " ~ ", b))
        }
        separated <- if (is.na(s_ab)) NA else (s_ab >= 0)
        if (singleton_pair && singleton == "na") {
          ## d_AA = 0 is a convention, not a measurement: no class.
          cheng_class <- NA_character_
          cheng_class_sign <- NA_character_
          if (!is.na(s_ab)) n_singleton_na <- n_singleton_na + 1L
        } else {
          cheng_class <- .synergy_cheng_class(separated, proximal_a, proximal_b)
          cheng_class_sign <- .synergy_cheng_class(separated, negative_a, negative_b)
        }
        complementary_exposure <- if (is.na(cheng_class)) NA else (cheng_class == "P2")
        synergy_score <- if (isTRUE(complementary_exposure) && !singleton_pair) {
          complementarity * joint_closeness
        } else {
          NA_real_
        }
        n_a_mapped_out <- n_a_mapped; n_b_mapped_out <- n_b_mapped
      } else {
        n_a_mapped_out <- NA_integer_; n_b_mapped_out <- NA_integer_
        singleton_a <- NA; singleton_b <- NA
        d_aa <- NA_real_; d_bb <- NA_real_; d_ab <- NA_real_; s_ab <- NA_real_
        separated <- NA
        cheng_class <- NA_character_
        cheng_class_sign <- NA_character_
        complementary_exposure <- NA
        synergy_score <- if (both_proximal) complementarity * joint_closeness else NA_real_
      }

      pair_rows[[i]] <- data.frame(
        condition = cond, disease_id = disease, compound_a = a, compound_b = b,
        n_targets_a = n_targets_a, n_targets_b = n_targets_b,
        n_targets_a_mapped = n_a_mapped_out, n_targets_b_mapped = n_b_mapped_out,
        target_jaccard = target_jaccard, complementarity = complementarity,
        singleton_a = singleton_a, singleton_b = singleton_b,
        d_aa = d_aa, d_bb = d_bb, d_ab = d_ab, s_ab = s_ab, separated = separated,
        z_score_a = za, z_score_b = zb, p_adjusted_a = pa, p_adjusted_b = pb,
        proximal_a = proximal_a, proximal_b = proximal_b, both_proximal = both_proximal,
        both_negative_z = both_negative_z,
        cheng_class = cheng_class, cheng_class_sign = cheng_class_sign,
        complementary_exposure = complementary_exposure,
        joint_closeness = joint_closeness, synergy_score = synergy_score,
        separation_method = separation, alpha = alpha, singleton_policy = singleton,
        pairs_mode = pairs,
        species = as.numeric(species), string_version = as.character(version),
        score_threshold = as.numeric(score_threshold),
        disease_gene_source = effective_dgs,
        stringsAsFactors = FALSE
      )
    }

    rows[[cond]] <- do.call(rbind, pair_rows)
    n_p2 <- sum(rows[[cond]]$cheng_class == "P2", na.rm = TRUE)
    n_p2_sign <- sum(rows[[cond]]$cheng_class_sign == "P2", na.rm = TRUE)
    n_scored <- sum(!is.na(rows[[cond]]$synergy_score))
    proj <- .log_append(
      proj, step = "network_synergy", id = NA_character_,
      message = paste0(
        "condition '", cond, "', disease '", disease, "' (separation = '", separation,
        "', pairs = '", pairs, "'): ", nrow(rows[[cond]]), " pair(s), ", n_p2,
        " P2 (Complementary Exposure, FDR-gated), ", n_p2_sign, " P2 (sign-only rule), ",
        n_scored, " with a synergy_score"
      )
    )
    if (n_singleton_na > 0) {
      proj <- .log_append(
        proj, step = "network_synergy", id = NA_character_,
        message = paste0(
          "condition '", cond, "': ", n_singleton_na, " pair(s) involve a compound with < 2 mapped targets;",
          " cheng_class / cheng_class_sign / complementary_exposure set to NA (singleton = 'na')"
        )
      )
    }
    if (p_fallback_used) {
      proj <- .log_append(
        proj, step = "network_synergy", id = NA_character_,
        message = paste0("condition '", cond, "': some compounds had no p_adjusted; proximal_* fell back to the z-sign alone (alpha not applied) for those")
      )
    }
    if (length(disconnected_pairs) > 0) {
      proj <- .log_append(
        proj, step = "network_synergy", id = NA_character_,
        message = paste0(
          "condition '", cond, "': ", length(disconnected_pairs),
          " pair(s) had a non-finite topological distance (disconnected on the STRING LCC); s_ab set to NA: ",
          paste(utils::head(disconnected_pairs, 10), collapse = ", ")
        )
      )
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  ## Replace every pair for the recomputed condition/disease, including
  ## pairs involving compounds no longer present in the current network.
  touched_keys <- data.frame(condition = conditions, disease_id = disease, stringsAsFactors = FALSE)
  result <- .network_upsert(
    proj, "network_synergy", result,
    c("condition", "disease_id"), touched_keys = touched_keys
  )

  patliRResults(proj, "network_synergy") <- result
  .write_results_csv(proj, "network_synergy", result)
  .write_log_csv(proj)
  proj
}

#' Three-valued proximity predicate for the Cheng classification
#'
#' @description
#' `NA` when `z` is missing; otherwise `z < 0 & p_adj < alpha`, falling
#' back per-row to `z < 0` alone when that row's `p_adj` is `NA` (the
#' caller logs that `alpha` was not applied).
#' @keywords internal
.synergy_proximal <- function(z, p_adj, alpha) {
  if (is.na(z)) return(NA)
  if (is.na(p_adj)) return(z < 0)
  z < 0 && p_adj < alpha
}

#' Cheng, Kovacs & Barabasi (2019) P1..P6 drug-pair class
#'
#' @description
#' `NA_character_` whenever any of `separated`, `proximal_a`, `proximal_b`
#' is `NA`; otherwise the class from Cheng et al.'s Fig. 2 truth table.
#' `P2` (Complementary Exposure) is exactly `separated & proximal_a &
#' proximal_b`. [network_synergy()] calls it twice: with the FDR-gated
#' predicates (`cheng_class`) and with the raw `z < 0` signs
#' (`cheng_class_sign`).
#' @keywords internal
.synergy_cheng_class <- function(separated, proximal_a, proximal_b) {
  if (anyNA(c(separated, proximal_a, proximal_b))) return(NA_character_)
  n_prox <- sum(proximal_a, proximal_b)
  if (!separated && n_prox == 2L) {
    "P1"
  } else if (separated && n_prox == 2L) {
    "P2"
  } else if (!separated && n_prox == 1L) {
    "P3"
  } else if (separated && n_prox == 1L) {
    "P4"
  } else if (!separated && n_prox == 0L) {
    "P5"
  } else {
    "P6"
  }
}

#' Menche et al. (2015) within-set distance <d_AA> (closest convention)
#'
#' @description
#' `unique()`s `ids` first (igraph's C `distances()` rejects a `to=` with
#' duplicate vertices, and the deduped size is the denominator Menche uses).
#' Returns `0` when `length(unique(ids)) < 2` -- **patliR's** singleton
#' convention, not Menche's (`separation.py` returns `nan`); see
#' [network_synergy()]. Any non-finite distance between two set members
#' (only possible off the LCC) makes the whole quantity `NA_real_` rather
#' than silently shrinking the set.
#'
#' @param D Optional pre-computed distance matrix with dimnames covering
#'   `ids` (so a caller can compute one `distances()` per condition and
#'   every set/pair reads a submatrix).
#' @return Single numeric, or `NA_real_` / `0` per the conventions above.
#' @keywords internal
.network_set_distance <- function(g, ids, D = NULL) {
  ids <- unique(ids)
  if (length(ids) < 2L) return(0)
  d <- if (is.null(D)) {
    igraph::distances(g, v = ids, to = ids, weights = NA)
  } else {
    D[ids, ids, drop = FALSE]
  }
  ## catch Inf *and* NaN before the diagonal is blanked (the diagonal of a
  ## real distance matrix is 0, so it survives the finiteness test).
  if (any(!is.finite(d))) return(NA_real_)
  diag(d) <- NA_real_
  row_min <- apply(d, 1, function(r) if (all(is.na(r))) NA_real_ else min(r, na.rm = TRUE))
  if (anyNA(row_min)) return(NA_real_)
  mean(row_min)
}

#' Menche et al. (2015) between-set distance <d_AB> (pooled, symmetric)
#'
#' @description
#' One `igraph::distances()` matrix, `apply(d, 1, min)` walking A -> B and
#' `apply(d, 2, min)` walking B -> A, denominator `nrow(d) + ncol(d)` on
#' the **deduped** sets. Nodes in `A n B` contribute `d = 0` on both sides
#' (intended -- Menche's `separation.py` does the same). Any non-finite
#' entry makes the whole quantity `NA_real_` (do not copy
#' `.network_closest_distance()`'s silent set-shrinking).
#'
#' @inheritParams .network_set_distance
#' @return Single numeric, or `NA_real_` if either set is empty or any
#'   distance is non-finite.
#' @keywords internal
.network_between_set_distance <- function(g, a, b, D = NULL) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0L || length(b) == 0L) return(NA_real_)
  d <- if (is.null(D)) {
    igraph::distances(g, v = a, to = b, weights = NA)
  } else {
    D[a, b, drop = FALSE]
  }
  if (any(!is.finite(d))) return(NA_real_)
  sum(apply(d, 1, min), apply(d, 2, min)) / (nrow(d) + ncol(d))
}

#' Menche et al. (2015) network separation s_AB
#'
#' @description
#' `s_AB = <d_AB> - (<d_AA> + <d_BB>) / 2`. `NA_real_` if any of the three
#' component distances is `NA`. `s_AB >= 0` means the two modules are
#' topologically separated (Menche's sign convention).
#' @inheritParams .network_set_distance
#' @return Single numeric or `NA_real_`.
#' @keywords internal
.network_separation <- function(g, a, b, D = NULL) {
  d_ab <- .network_between_set_distance(g, a, b, D)
  d_aa <- .network_set_distance(g, a, D)
  d_bb <- .network_set_distance(g, b, D)
  if (anyNA(c(d_ab, d_aa, d_bb))) return(NA_real_)
  d_ab - (d_aa + d_bb) / 2
}

#' @keywords internal
.empty_network_synergy_row <- function() {
  data.frame(
    condition = character(0), disease_id = character(0),
    compound_a = character(0), compound_b = character(0),
    n_targets_a = integer(0), n_targets_b = integer(0),
    n_targets_a_mapped = integer(0), n_targets_b_mapped = integer(0),
    target_jaccard = double(0), complementarity = double(0),
    singleton_a = logical(0), singleton_b = logical(0),
    d_aa = double(0), d_bb = double(0), d_ab = double(0), s_ab = double(0),
    separated = logical(0),
    z_score_a = double(0), z_score_b = double(0),
    p_adjusted_a = double(0), p_adjusted_b = double(0),
    proximal_a = logical(0), proximal_b = logical(0), both_proximal = logical(0),
    both_negative_z = logical(0),
    cheng_class = character(0), cheng_class_sign = character(0),
    complementary_exposure = logical(0),
    joint_closeness = double(0), synergy_score = double(0),
    separation_method = character(0), alpha = double(0), singleton_policy = character(0),
    pairs_mode = character(0),
    species = double(0), string_version = character(0), score_threshold = double(0),
    disease_gene_source = character(0),
    stringsAsFactors = FALSE
  )
}
