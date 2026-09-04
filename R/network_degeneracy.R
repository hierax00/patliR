#' @include AllGenerics.R internal.R network_build.R network_layers.R
NULL

## "Degeneracy" in the systems-biology sense (Edelman & Gally, 2001, PNAS
## 98(24), 13763-8): two compounds are degenerate when they act on largely
## DIFFERENT targets but those targets converge on largely the SAME
## function -- redundant routes to the same effect. A pairwise compound
## property.
##
## annotation = "direct" (default) scores that convergence with GO semantic
## similarity (GOSemSim, Wang 2007) between the two compounds' targets'
## OWN GO annotations, against an annotation-count-matched permutation null.
## It needs only network_build() (+ the UniProt -> Entrez map and GOSemSim),
## NOT network_enrich() -- a gene's GO annotations do not move when the
## condition's compound list changes, which is what removes the
## enrichment-universe circularity that annotation = "enriched" (and the
## legacy annotation = "jaccard") carry (see the @section below).

#' Pairwise functional degeneracy between compounds (GO semantic similarity)
#'
#' @description
#' For each condition, and every pair of compounds present in that
#' condition's network ([network_build()]), computes:
#'
#' - `target_jaccard`: Jaccard similarity of the two compounds' target sets
#'   (structural overlap -- do they act on the same proteins?).
#' - `functional_similarity`: how strongly the two compounds' targets
#'   converge on the same biology. What it is depends on `annotation`:
#'   - `"direct"` (default) -- GO semantic similarity
#'     ([GOSemSim::clusterSim()]-equivalent) between the two compounds'
#'     Entrez-mapped target sets, scored on each protein's *own* GO
#'     annotations. No enrichment involved.
#'   - `"enriched"` -- [GOSemSim::mgoSim()] between the two compounds'
#'     *enriched* GO-term sets (from [network_enrich()], reconstructed via
#'     `.network_target_pathway_edges()`, internal). Carries the
#'     circularity described below.
#'   - `"jaccard"` -- the historical unweighted Jaccard of the two
#'     compounds' enriched-pathway sets (no GOSemSim). Regression path.
#' - `degeneracy_score = functional_similarity * (1 - target_jaccard)`:
#'   high only when two compounds are functionally similar
#'   (`functional_similarity` near 1) *despite* being structurally distinct
#'   (`target_jaccard` near 0) -- different elements, same function, in the
#'   spirit of degeneracy as Edelman & Gally (2001, *PNAS* 98(24), 13763-8)
#'   describe it qualitatively (their quantitative measure, Tononi et al.
#'   1999, is a different, information-theoretic construction). A pair that
#'   hits the same targets is redundant, not degenerate; a pair that shares
#'   neither targets nor function is unrelated. Only the different-targets /
#'   same-function quadrant scores high. The multiplicative form is a
#'   convenience scalarisation, not a derived result; read `target_jaccard`
#'   and `functional_similarity` together.
#'
#' @section The permutation null (`annotation = "direct"` only):
#' Proteins differ ~100-fold in how many GO BP terms they carry (median 6,
#' max > 270 in `org.Hs.eg.db`), and a random set of poorly-annotated
#' proteins scores low for reasons unrelated to biology. For each pair the
#' observed similarity is z-scored against `n_random` resamplings in which
#' each compound's gene set is replaced by an equal-size set of *distinct*
#' genes drawn from the same annotation-count bins (the same
#' consecutive-value >= `min_per_bin` binning [network_proximity()] uses for
#' degree, via the shared `.network_value_bins()` /
#' `.network_resample_matched()`). The resampling pool is `universe`:
#' `"project"` (default -- every imported target, so the null carries the
#' target predictor's family bias, which is what you want), `"condition"`
#' (this condition's targets only -- often too small), or `"genome"` (every
#' annotated gene -- the weakest null; capped for tractability and warned).
#'
#' `z_score` is the primary reported statistic (as in Guney/Menche/Cheng).
#' `p_empirical` is the **right-tail** permutation p
#' (`(1 + #{sim* >= sim_obs}) / (n_random + 1)`) -- high similarity is the
#' interesting direction, so the tail is flipped relative to
#' [network_proximity()], whose identically-named `p_adjusted` is a
#' *left*-tail p. `p_adjusted` is its Benjamini-Hochberg value across every
#' pair in the call. Because the GOSemSim score is quantised to 0.001, ties
#' at the observed value are common and inflate `p_empirical` (conservative);
#' `sd(sim*) == 0` is a realistic outcome for small pools and gives
#' `z_score = NA`. With `n_pairs / (n_random + 1) > 0.05` a BH-adjusted p
#' cannot resolve an isolated pair and a warning is emitted -- lead with
#' `z_score` there.
#'
#' @section Why the default does not use network_enrich():
#' In `annotation = "enriched"` mode a compound's functional profile is the
#' set of terms that were called significant for the **whole condition** --
#' a pool the compound's own targets helped define, shared by every compound
#' in the condition. Two compounds are then compared on membership in a
#' small, mutually-determined list, which inflates the similarity of any
#' pair sharing even one promiscuous target and makes the score depend on
#' how many *other* compounds were in the condition. `annotation = "direct"`
#' compares each compound's targets by their own GO annotations, which are
#' fixed properties of the proteins and do not move when the condition's
#' compound list changes. The permutation null corrects the scale but not
#' this dependence, so `"direct"` is the default. `measure = "Wang"` is
#' chosen for consistency with `network_enrich(simplify_go = TRUE)` (which
#' also uses Wang) and for cost, **not** because it is bias-free -- Wang
#' removes the term-level frequency dependence that Resnik/Lin/Rel/Jiang
#' carry, but not the gene-level one (well-studied proteins carry more
#' terms), and UniProt->Entrez and the GO annotations both come from
#' `org.Hs.eg.db`. It is the annotation-count-matched null that corrects the
#' curation bias.
#'
#' @section Requires network_enrich() only for `"enriched"` / `"jaccard"`:
#' `annotation = "direct"` (default) needs only [network_build()], the
#' `clusterProfiler` + `org.Hs.eg.db` UniProt->Entrez map, and `GOSemSim` --
#' it does not need [network_enrich()] and is usable one step after
#' [network_build()], like [network_centrality()]. `"enriched"` and
#' `"jaccard"` still require [network_enrich()] to have run for the
#' condition; `"jaccard"` does not need `GOSemSim`.
#'
#' @section Compound pairs with nothing to compare are omitted, not `NA`-filled:
#' In `"direct"` mode a pair is skipped (no row) when either compound has no
#' Entrez-mapped, GO-annotated target. In `"enriched"` / `"jaccard"` mode a
#' pair is skipped when neither compound has a target in a significant
#' pathway (`functional_similarity` would be `0/0`). The row *set* is
#' therefore mode-dependent; see `projectLog(proj)` for the skip counts.
#'
#' @inheritParams network_build
#' @param annotation How `functional_similarity` is computed: `"direct"`
#'   (default, GO semantic similarity on the targets' own annotations, no
#'   `network_enrich()`), `"enriched"` (GO semantic similarity on the
#'   enriched-term sets), or `"jaccard"` (historical enriched-pathway
#'   Jaccard).
#' @param ont GO ontology for `"direct"` / `"enriched"`: `"BP"` (default --
#'   process-level convergence, the concept degeneracy requires), `"MF"`
#'   (molecular action -- answers "do these compounds engage the same
#'   *kind* of protein", a different question), or `"CC"`.
#' @param measure GO semantic-similarity measure: `"Wang"` (default,
#'   `computeIC = FALSE`) or an IC measure (`"Resnik"`, `"Lin"`, `"Rel"`,
#'   `"Jiang"` -- each needs a separate `computeIC = TRUE` `godata` object,
#'   logged when built; Resnik is unbounded).
#' @param combine Set-to-set combination: `"BMA"` (default, best-match
#'   average), `"max"`, `"avg"`, or `"rcmax"`.
#' @param drop GO evidence codes excluded before scoring / before the
#'   annotation-count binning. Default `"IEA"` (electronic annotations --
#'   ~44% of `org.Hs.eg.db` BP annotations and the entire BP annotation of
#'   ~1,966 genes; Cheng et al. 2019 exclude them too). `NULL` keeps all.
#' @param universe Resampling pool for the permutation null: `"project"`
#'   (default), `"condition"`, or `"genome"`.
#' @param n_random Permutation-null resamplings for the z-score
#'   (`"direct"` mode). Default `200`; `>= 1` required. `n_random = 1`
#'   gives `sim_random_sd = NA`, `z_score = NA`, `p_empirical` in
#'   `{0.5, 1}`.
#' @param seed `NULL` (default, auto-generated and logged) or an integer,
#'   for the resampling. Does not leak into the caller's RNG state.
#' @param pathway_db See `.network_target_pathway_edges()`, internal --
#'   restrict `"enriched"` / `"jaccard"` pathway sets to specific
#'   `network_enrich()` `db` value(s), or `NULL` (default) for all.
#'
#' @return The updated `proj`, with a `network_degeneracy` entry in
#'   [patliRResults()] (columns `condition`, `compound_a`, `compound_b`,
#'   `n_targets_a`, `n_targets_b`, `n_genes_a_mapped`, `n_genes_b_mapped`,
#'   `n_pathways_a`, `n_pathways_b` (`NA` in `"direct"` mode),
#'   `target_jaccard`, `functional_similarity`, `sim_random_mean`,
#'   `sim_random_sd`, `z_score`, `p_empirical`, `p_adjusted` (right-tail;
#'   `NA` outside `"direct"` mode), `degeneracy_score`, and the per-call
#'   constants `annotation`, `ont`, `measure`, `combine`, `drop`,
#'   `universe`, `n_random`, `seed_used`), also written to
#'   `results/network_degeneracy.csv`. The pre-0.2.0 `pathway_jaccard`
#'   column is transparently renamed to `functional_similarity` on the
#'   first run.
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
#' # default "direct" mode -- no network_enrich() needed
#' proj <- network_degeneracy(proj, condition = "FLO-ET", seed = 1)
#' patliRResults(proj, "network_degeneracy")
#' }
#'
#' @export
network_degeneracy <- function(proj, condition = NULL,
                                annotation = c("direct", "enriched", "jaccard"),
                                ont = c("BP", "MF", "CC"),
                                measure = c("Wang", "Resnik", "Lin", "Rel", "Jiang"),
                                combine = c("BMA", "max", "avg", "rcmax"),
                                drop = "IEA",
                                universe = c("project", "condition", "genome"),
                                n_random = 200, seed = NULL, pathway_db = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  annotation <- match.arg(annotation)
  ont <- match.arg(ont)
  measure <- match.arg(measure)
  combine <- match.arg(combine)
  universe <- match.arg(universe)
  stopifnot(is.numeric(n_random), length(n_random) == 1, n_random >= 1)
  if (!is.null(drop)) stopifnot(is.character(drop))

  conditions <- .network_resolve_conditions(proj, condition)

  ## ---- mode-aware prerequisite guards (14-R3) --------------------------
  needs_enrich <- annotation %in% c("enriched", "jaccard")
  needs_cp     <- annotation %in% c("direct", "enriched", "jaccard")
  needs_gosem  <- annotation %in% c("direct", "enriched")

  if (needs_enrich) {
    enrichment_all <- patliRResults(proj, "network_enrichment")
    if (is.null(enrichment_all) || nrow(enrichment_all) == 0) {
      cli::cli_abort(c(
        "No {.val network_enrichment} entry in {.arg proj}.",
        "i" = "{.code annotation = {.val {annotation}}} needs the pathway layer -- run {.fn network_enrich} first, or use the default {.code annotation = \"direct\"}."
      ))
    }
    missing_enrich <- setdiff(conditions, unique(enrichment_all$condition))
    if (length(missing_enrich) > 0) {
      cli::cli_abort(c(
        "Condition(s) {.val {missing_enrich}} have no {.fn network_enrich} results.",
        "i" = "Run {.fn network_enrich} for {.val {missing_enrich}} first, or use {.code annotation = \"direct\"}."
      ))
    }
  }
  if (needs_cp &&
      (!requireNamespace("clusterProfiler", quietly = TRUE) ||
       !requireNamespace("org.Hs.eg.db", quietly = TRUE))) {
    cli::cli_abort(c(
      "{.fn network_degeneracy} needs {.pkg clusterProfiler} and {.pkg org.Hs.eg.db} (UniProt -> Entrez mapping), not installed.",
      "i" = "Install them with {.code BiocManager::install(c(\"clusterProfiler\", \"org.Hs.eg.db\"))}."
    ))
  }
  if (needs_gosem && !requireNamespace("GOSemSim", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.code annotation = {.val {annotation}}} needs {.pkg GOSemSim}, not installed.",
      "i" = "Install it with {.code BiocManager::install(\"GOSemSim\")}, or use {.code annotation = \"jaccard\"}."
    ))
  }

  used_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else as.integer(seed)
  restore_rng <- .with_seed(used_seed)
  on.exit(restore_rng(), add = TRUE)

  ## ---- GOSemSim godata + annotation index (direct/enriched) -----------
  sem_data <- NULL
  g2go <- NULL
  ann_counts <- NULL
  if (needs_gosem) {
    compute_ic <- measure != "Wang"
    sem_data <- .network_godata(proj, ont, computeIC = compute_ic)
    if (is.null(sem_data)) {
      cli::cli_abort(c(
        "{.fn GOSemSim::godata} could not be built for {.val {ont}} (computeIC = {compute_ic}).",
        "i" = "Check {.pkg org.Hs.eg.db} / {.pkg GO.db} are installed and loadable."
      ))
    }
    proj <- .log_append(
      proj, step = "network_degeneracy", id = NA_character_,
      message = paste0(
        "GOSemSim godata built/loaded: ont=", ont, ", measure=", measure,
        ", computeIC=", compute_ic, ", drop=", if (is.null(drop)) "none" else paste(drop, collapse = "+")
      )
    )
    ga <- sem_data@geneAnno
    if (!is.null(drop) && "EVIDENCE" %in% names(ga)) ga <- ga[!ga$EVIDENCE %in% drop, , drop = FALSE]
    g2go <- lapply(split(as.character(ga$GO), as.character(ga$ENTREZID)),
                   function(x) unique(x[!is.na(x)]))
    ann_counts <- lengths(g2go)
  }

  edges_all <- patliRResults(proj, "network_edges")

  ## project universe -> Entrez, once (direct mode)
  project_pool_entrez <- NULL
  if (annotation == "direct" && universe == "project") {
    ti <- patliRResults(proj, "targets_imported")
    proj_uni <- if (!is.null(ti) && nrow(ti) > 0) unique(ti$uniprot_id) else unique(edges_all$uniprot_id)
    project_pool_entrez <- .network_uniprot_to_entrez(proj_uni)
  }

  rows <- vector("list", length(conditions))
  names(rows) <- conditions
  n_skipped_total <- 0L

  for (cond in conditions) {
    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    all_compounds <- sort(unique(ct$compound_id))

    if (length(all_compounds) < 2) {
      rows[[cond]] <- .network_degeneracy_direct_empty()
      proj <- .log_append(
        proj, step = "network_degeneracy", id = NA_character_,
        message = paste0("condition '", cond, "': fewer than 2 compounds in the network, nothing to compare")
      )
      next
    }

    if (annotation == "direct") {
      rows[[cond]] <- .network_degeneracy_direct(
        proj, cond, ct, all_compounds, sem_data, g2go, ann_counts,
        measure, combine, universe, n_random, project_pool_entrez, edges_all
      )
      proj <- rows[[cond]]$proj
      rows[[cond]] <- rows[[cond]]$rows
    } else {
      rows[[cond]] <- .network_degeneracy_pathway(
        proj, cond, ct, all_compounds, annotation, sem_data, measure, combine, pathway_db
      )
      proj <- rows[[cond]]$proj
      n_skipped_total <- n_skipped_total + rows[[cond]]$n_skipped
      rows[[cond]] <- rows[[cond]]$rows
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  ## per-call constants
  n_call <- if (is.null(result)) 0L else nrow(result)
  if (n_call > 0) {
    result$annotation <- annotation
    result$ont <- if (annotation == "jaccard") NA_character_ else ont
    result$measure <- if (annotation == "jaccard") NA_character_ else measure
    result$combine <- if (annotation == "jaccard") NA_character_ else combine
    result$drop <- if (annotation == "jaccard" || is.null(drop)) NA_character_ else paste(drop, collapse = "+")
    result$universe <- if (annotation == "direct") universe else NA_character_
    result$n_random <- if (annotation == "direct") as.integer(n_random) else NA_integer_
    result$seed_used <- if (annotation == "direct") used_seed else NA_integer_

    ## right-tail BH across every pair in the call (direct mode only)
    if (annotation == "direct") {
      result$p_adjusted <- stats::p.adjust(result$p_empirical, method = "BH")
      if (n_call / (n_random + 1) > 0.05) {
        cli::cli_warn(c(
          "!" = "{.fn network_degeneracy}: {n_call} pair(s) vs n_random = {n_random}: a BH-adjusted p cannot resolve an isolated convergent pair (floor {signif(1/(n_random+1), 3)} x {n_call}).",
          "i" = "Lead with {.field z_score}; raise {.arg n_random} if you need {.field p_adjusted}."
        ))
      }
    }
  } else {
    result <- .empty_network_degeneracy_row()
  }

  ## ---- pathway_jaccard -> functional_similarity migration (14-U3) -----
  existing <- patliRResults(proj, "network_degeneracy")
  if (!is.null(existing) && "pathway_jaccard" %in% names(existing) &&
      !"functional_similarity" %in% names(existing)) {
    names(existing)[names(existing) == "pathway_jaccard"] <- "functional_similarity"
    patliRResults(proj, "network_degeneracy") <- existing
    proj <- .log_append(
      proj, step = "network_degeneracy", id = NA_character_,
      message = "migrated legacy column 'pathway_jaccard' -> 'functional_similarity'"
    )
  }

  result <- .network_upsert(
    proj, "network_degeneracy", result, "condition",
    touched_keys = data.frame(condition = conditions, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "network_degeneracy") <- result
  .write_results_csv(proj, "network_degeneracy", result)
  .write_log_csv(proj)
  proj
}

#' `annotation = "direct"` core: GO semantic similarity + permutation null
#' for one condition. Returns `list(rows, proj)`.
#' @keywords internal
.network_degeneracy_direct <- function(proj, cond, ct, all_compounds, sem_data,
                                        g2go, ann_counts, measure, combine, universe,
                                        n_random, project_pool_entrez, edges_all) {
  ## compound -> distinct Entrez (many-to-many UniProt -> Entrez)
  uni_all <- unique(ct$uniprot_id)
  u2e <- .network_uniprot_to_entrez_map(uni_all)
  compound_entrez <- lapply(split(ct$uniprot_id, ct$compound_id), function(u) {
    unique(stats::na.omit(unlist(u2e[unique(u)], use.names = FALSE)))
  })

  ## resampling pool
  pool <- switch(universe,
    project   = project_pool_entrez,
    condition = unique(stats::na.omit(unlist(u2e, use.names = FALSE))),
    genome    = names(g2go)
  )
  pool <- unique(as.character(pool))
  pool <- pool[pool %in% names(g2go)]        # annotated genes only (14-U2)

  ## cap the genome pool so the term x term matrix stays tractable
  pool_capped <- FALSE
  cap <- 1500L
  if (length(pool) > cap) {
    pool <- sample(pool, cap)
    pool_capped <- TRUE
  }

  if (length(pool) < 2) {
    proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
      message = paste0("condition '", cond, "': resampling pool (universe = '", universe,
                       "') has < 2 annotated genes; no rows"))
    return(list(rows = .empty_network_degeneracy_row(), proj = proj))
  }

  ## annotation-count bins over the pool (shared helper with proximity)
  pool_counts <- ann_counts[pool]
  bins <- .network_value_bins(pool_counts, min_per_bin = 100)
  bin_of_node <- integer(length(pool))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b
  if (length(bins) < 3) {
    cli::cli_warn(c(
      "!" = "{.fn network_degeneracy}: condition {.val {cond}} resampling pool (universe = {.val {universe}}, {length(pool)} gene{?s}) yields only {length(bins)} annotation-count bin{?s}.",
      "i" = "The permutation null is close to uniform sampling; a larger {.arg universe} gives a sharper null."
    ))
  }
  if (pool_capped) {
    cli::cli_warn(c(
      "!" = "{.fn network_degeneracy}: condition {.val {cond}} universe = {.val {universe}} pool capped at {cap} genes for the term-similarity matrix."
    ))
  }

  ## term x term similarity matrix over the pool's term union, ONCE
  terms_pool <- sort(unique(unlist(g2go[pool], use.names = FALSE)))
  sim_mat <- GOSemSim::termSim(terms_pool, terms_pool, semData = sem_data, method = measure)
  if (is.null(dim(sim_mat))) {
    sim_mat <- matrix(sim_mat, nrow = length(terms_pool), ncol = length(terms_pool),
                      dimnames = list(terms_pool, terms_pool))
  }

  ## genes that lost all annotation after the drop filter (log once)
  lost <- setdiff(unique(unlist(compound_entrez, use.names = FALSE)), names(g2go))
  if (length(lost) > 0) {
    proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
      message = paste0("condition '", cond, "': ", length(lost),
                       " Entrez gene(s) had no GO annotation after the drop filter and were excluded from scoring"))
  }

  pairs <- utils::combn(all_compounds, 2, simplify = FALSE)
  pair_rows <- vector("list", length(pairs))
  n_skipped <- 0L

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1]; b <- pairs[[i]][2]
    ea_all <- compound_entrez[[a]]; eb_all <- compound_entrez[[b]]
    if (is.null(ea_all)) ea_all <- character(0)
    if (is.null(eb_all)) eb_all <- character(0)
    ## scored sets: annotated AND in the resampling pool (14-U2)
    ea <- intersect(ea_all, pool)
    eb <- intersect(eb_all, pool)

    if (length(ea) == 0 || length(eb) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }

    tv_a <- unlist(g2go[ea], use.names = FALSE)
    tv_b <- unlist(g2go[eb], use.names = FALSE)
    sim_obs <- .network_go_bma(tv_a, tv_b, sim_mat, combine)

    sim_rand <- vapply(seq_len(n_random), function(j) {
      ra <- .network_resample_matched(ea, pool, bins, bin_of_node)
      rb <- .network_resample_matched(eb, pool, bins, bin_of_node)
      .network_go_bma(unlist(g2go[ra], use.names = FALSE),
                      unlist(g2go[rb], use.names = FALSE), sim_mat, combine, check = FALSE)
    }, numeric(1))
    sim_rand <- sim_rand[is.finite(sim_rand)]

    sim_mean <- if (length(sim_rand) > 0) mean(sim_rand) else NA_real_
    sim_sd   <- if (length(sim_rand) > 1) stats::sd(sim_rand) else NA_real_
    z_score  <- if (!is.na(sim_sd) && sim_sd > 0 && !is.na(sim_obs)) {
      (sim_obs - sim_mean) / sim_sd
    } else NA_real_
    p_emp <- if (length(sim_rand) > 0 && !is.na(sim_obs)) {
      (1 + sum(sim_rand >= sim_obs)) / (length(sim_rand) + 1)
    } else NA_real_

    tj <- length(intersect(ct$uniprot_id[ct$compound_id == a],
                           ct$uniprot_id[ct$compound_id == b])) /
          length(union(ct$uniprot_id[ct$compound_id == a],
                       ct$uniprot_id[ct$compound_id == b]))
    deg <- if (is.na(sim_obs)) NA_real_ else sim_obs * (1 - tj)

    pair_rows[[i]] <- data.frame(
      condition = cond, compound_a = a, compound_b = b,
      n_targets_a = length(ct$uniprot_id[ct$compound_id == a]),
      n_targets_b = length(ct$uniprot_id[ct$compound_id == b]),
      n_genes_a_mapped = length(ea_all), n_genes_b_mapped = length(eb_all),
      n_pathways_a = NA_integer_, n_pathways_b = NA_integer_,
      target_jaccard = tj, functional_similarity = sim_obs,
      sim_random_mean = sim_mean, sim_random_sd = sim_sd,
      z_score = z_score, p_empirical = p_emp, p_adjusted = NA_real_,
      degeneracy_score = deg,
      stringsAsFactors = FALSE
    )
  }

  pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
  out <- if (length(pair_rows) == 0) {
    .network_degeneracy_direct_empty()
  } else {
    do.call(rbind, pair_rows)
  }
  proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
    message = paste0("condition '", cond, "' (annotation = 'direct'): ", nrow(out),
                     " pair(s) scored, ", n_skipped,
                     " pair(s) skipped (a compound has no Entrez-mapped, annotated target), seed-driven null n_random = ", n_random))
  list(rows = out, proj = proj)
}

#' `annotation = "enriched"` / `"jaccard"` core for one condition.
#' Returns `list(rows, proj, n_skipped)`.
#' @keywords internal
.network_degeneracy_pathway <- function(proj, cond, ct, all_compounds, annotation,
                                         sem_data, measure, combine, pathway_db) {
  target_pathway <- .network_target_pathway_edges(proj, cond, unique(ct$uniprot_id), pathway_db = pathway_db)
  cp <- if (!is.null(target_pathway) && nrow(target_pathway) > 0) {
    unique(merge(ct, target_pathway, by = "uniprot_id")[, c("compound_id", "pathway_id")])
  } else {
    data.frame(compound_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE)
  }
  target_sets <- split(ct$uniprot_id, ct$compound_id)
  pathway_sets <- split(cp$pathway_id, cp$compound_id)

  pairs <- utils::combn(all_compounds, 2, simplify = FALSE)
  pair_rows <- vector("list", length(pairs))
  n_skipped <- 0L

  for (i in seq_along(pairs)) {
    a <- pairs[[i]][1]; b <- pairs[[i]][2]
    ta <- target_sets[[a]]; tb <- target_sets[[b]]
    pa <- pathway_sets[[a]]; pb <- pathway_sets[[b]]
    if (is.null(pa)) pa <- character(0)
    if (is.null(pb)) pb <- character(0)

    if (length(union(pa, pb)) == 0) {
      n_skipped <- n_skipped + 1L
      next
    }
    target_jaccard <- length(intersect(ta, tb)) / length(union(ta, tb))

    fs <- if (annotation == "jaccard") {
      length(intersect(pa, pb)) / length(union(pa, pb))
    } else {
      ## enriched: GO semantic similarity over the enriched-term sets
      go_a <- unique(pa[grepl("^GO:", pa)])
      go_b <- unique(pb[grepl("^GO:", pb)])
      if (length(go_a) == 0 || length(go_b) == 0) NA_real_
      else GOSemSim::mgoSim(go_a, go_b, semData = sem_data, measure = measure, combine = combine)
    }
    deg <- if (is.na(fs)) NA_real_ else fs * (1 - target_jaccard)

    pair_rows[[i]] <- data.frame(
      condition = cond, compound_a = a, compound_b = b,
      n_targets_a = length(ta), n_targets_b = length(tb),
      n_genes_a_mapped = NA_integer_, n_genes_b_mapped = NA_integer_,
      n_pathways_a = length(pa), n_pathways_b = length(pb),
      target_jaccard = target_jaccard, functional_similarity = fs,
      sim_random_mean = NA_real_, sim_random_sd = NA_real_,
      z_score = NA_real_, p_empirical = NA_real_, p_adjusted = NA_real_,
      degeneracy_score = deg,
      stringsAsFactors = FALSE
    )
  }

  pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
  out <- if (length(pair_rows) == 0) .network_degeneracy_direct_empty() else do.call(rbind, pair_rows)
  proj <- .log_append(proj, step = "network_degeneracy", id = NA_character_,
    message = paste0("condition '", cond, "' (annotation = '", annotation, "'): ", nrow(out),
                     " pair(s) scored, ", n_skipped,
                     " pair(s) skipped (neither compound has a target in a significant pathway)"))
  list(rows = out, proj = proj, n_skipped = n_skipped)
}

#' The per-pair columns, before the per-call constants are bound on.
#' @keywords internal
.network_degeneracy_direct_empty <- function() {
  data.frame(
    condition = character(0), compound_a = character(0), compound_b = character(0),
    n_targets_a = integer(0), n_targets_b = integer(0),
    n_genes_a_mapped = integer(0), n_genes_b_mapped = integer(0),
    n_pathways_a = integer(0), n_pathways_b = integer(0),
    target_jaccard = double(0), functional_similarity = double(0),
    sim_random_mean = double(0), sim_random_sd = double(0),
    z_score = double(0), p_empirical = double(0), p_adjusted = double(0),
    degeneracy_score = double(0),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_network_degeneracy_row <- function() {
  data.frame(
    condition = character(0), compound_a = character(0), compound_b = character(0),
    n_targets_a = integer(0), n_targets_b = integer(0),
    n_genes_a_mapped = integer(0), n_genes_b_mapped = integer(0),
    n_pathways_a = integer(0), n_pathways_b = integer(0),
    target_jaccard = double(0), functional_similarity = double(0),
    sim_random_mean = double(0), sim_random_sd = double(0),
    z_score = double(0), p_empirical = double(0), p_adjusted = double(0),
    degeneracy_score = double(0),
    annotation = character(0), ont = character(0), measure = character(0),
    combine = character(0), drop = character(0), universe = character(0),
    n_random = integer(0), seed_used = integer(0),
    stringsAsFactors = FALSE
  )
}

#' BMA (and max/avg/rcmax) over a pooled term-similarity submatrix,
#' reproducing `GOSemSim::combineScores()` exactly -- including the
#' all-NA row/column drop, the vector / one-row / one-column -> `max`
#' branch, and the `round(., 3)` quantisation. Duplicate terms in
#' `terms_a` / `terms_b` are kept (each contributes a row/column), so the
#' BMA is annotation-weighted the same way `clusterSim()` is.
#' @keywords internal
.network_go_bma <- function(terms_a, terms_b, sim_mat, combine = "BMA", check = TRUE) {
  if (check) {
    terms_a <- terms_a[!is.na(terms_a) & terms_a %in% rownames(sim_mat)]
    terms_b <- terms_b[!is.na(terms_b) & terms_b %in% colnames(sim_mat)]
  }
  if (length(terms_a) == 0 || length(terms_b) == 0) return(NA_real_)
  S <- sim_mat[terms_a, terms_b, drop = FALSE]
  if (!sum(!is.na(S))) return(NA_real_)

  .cs <- function(M) {
    if (is.vector(M) || nrow(M) == 1 || ncol(M) == 1) {
      if (combine == "avg") return(round(mean(M, na.rm = TRUE), 3))
      return(round(max(M, na.rm = TRUE), 3))
    }
    na <- is.na(M)
    rna <- rowSums(na) == ncol(M)
    if (any(rna)) { M <- M[!rna, , drop = FALSE]; na <- na[!rna, , drop = FALSE] }
    cna <- colSums(na) == nrow(M)
    if (any(cna)) M <- M[, !cna, drop = FALSE]
    if (is.vector(M) || nrow(M) == 1 || ncol(M) == 1) {
      if (combine == "avg") return(round(mean(M, na.rm = TRUE), 3))
      return(round(max(M, na.rm = TRUE), 3))
    }
    result <- switch(combine,
      avg   = mean(M, na.rm = TRUE),
      max   = max(M, na.rm = TRUE),
      rcmax = max(mean(apply(M, 1, max, na.rm = TRUE)),
                  mean(apply(M, 2, max, na.rm = TRUE))),
      ## BMA (rcmax.avg)
      sum(apply(M, 1, max, na.rm = TRUE), apply(M, 2, max, na.rm = TRUE)) / sum(dim(M))
    )
    round(result, 3)
  }
  .cs(S)
}

#' Build (and cache) a `GOSemSim::godata` object under `cacheDir(proj)`.
#'
#' @description
#' Factored so [network_enrich()] and [network_degeneracy()] share one
#' construction and one cache. Uses `annoDb = "org.Hs.eg.db"` (the
#' `OrgDb =` spelling is deprecated in GOSemSim >= 2.36). The cache key
#' carries the `org.Hs.eg.db` and `GO.db` package versions so a Bioconductor
#' annotation update is not silently served a stale corpus. Safe to delete:
#' the `.rds` is a ~1 MB rebuild of a ~30 s computation.
#'
#' @return A `GOSemSimDATA` object, or `NULL` if `GOSemSim` is unavailable
#'   or `godata()` errors (callers guard).
#' @keywords internal
.network_godata <- function(proj, ont, computeIC) {
  if (!requireNamespace("GOSemSim", quietly = TRUE)) return(NULL)
  ic_tag <- if (isTRUE(computeIC)) "ic" else "wang"
  org_ver <- tryCatch(as.character(utils::packageVersion("org.Hs.eg.db")), error = function(e) "NA")
  go_ver  <- tryCatch(as.character(utils::packageVersion("GO.db")), error = function(e) "NA")
  dir <- file.path(cacheDir(proj), "gosemsim")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  cache_rds <- file.path(dir, paste0("gosemsim_", ont, "_", ic_tag, "_", org_ver, "_", go_ver, ".rds"))
  if (file.exists(cache_rds)) {
    obj <- tryCatch(readRDS(cache_rds), error = function(e) NULL)
    if (!is.null(obj)) return(obj)
  }
  obj <- tryCatch(
    GOSemSim::godata(annoDb = "org.Hs.eg.db", keytype = "ENTREZID", ont = ont, computeIC = computeIC),
    error = function(e) NULL
  )
  if (!is.null(obj)) tryCatch(saveRDS(obj, cache_rds), error = function(e) NULL)
  obj
}

#' UniProt -> Entrez, returning the flat vector of distinct Entrez IDs.
#' @keywords internal
.network_uniprot_to_entrez <- function(uniprot_ids) {
  m <- .network_uniprot_to_entrez_map(uniprot_ids)
  unique(stats::na.omit(unlist(m, use.names = FALSE)))
}

#' UniProt -> Entrez as a named list (one character vector per input
#' accession; many-to-many). Uses the same `clusterProfiler::bitr` map as
#' [network_enrich()].
#' @keywords internal
.network_uniprot_to_entrez_map <- function(uniprot_ids) {
  uniprot_ids <- unique(uniprot_ids[!is.na(uniprot_ids) & nzchar(uniprot_ids)])
  if (length(uniprot_ids) == 0) return(list())
  map <- tryCatch(
    clusterProfiler::bitr(uniprot_ids, fromType = "UNIPROT", toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db", drop = TRUE),
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) {
    return(stats::setNames(vector("list", length(uniprot_ids)), uniprot_ids))
  }
  split(as.character(map$ENTREZID), map$UNIPROT)
}
