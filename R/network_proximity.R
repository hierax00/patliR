#' @include AllGenerics.R internal.R network_build.R network_module_robustness.R
NULL

## Network-based compound-disease proximity on the STRING PPI network. A
## per-condition compound-target bipartite graph has no real topological
## distance between two targets, so this needs a separate external PPI
## network.
##
## Method: Guney et al. (2016), Nat Commun 7, 10331 -- the "closest"
## proximity measure (their Eq. 2):
##   d(S,T) = (1/|S|) * sum_{s in S} min_{t in T} d(s,t)
## z-scored against a degree-preserving null model (S and T resampled from
## same-degree-bin nodes, n_random times). More negative z = closer than
## chance.

#' Network-based compound-disease proximity on the STRING interactome
#'
#' @description
#' For each compound present in a condition's network ([network_build()]),
#' computes the topological "closest" distance (Guney et al. 2016, *Nat
#' Commun* 7:10331) between that compound's target set and a disease's
#' associated gene set, on the STRING protein-protein interaction network
#' (`STRINGdb`), and z-scores it against a degree-preserving random null
#' model.
#'
#' @section Where the disease gene set comes from:
#' By default (`disease_genes = "disease_genes"`) the disease module `T` is
#' the independent, disease -> target gene set from [disease_genes_fetch()]
#' or [disease_genes_import()] -- built without reference to the compounds,
#' so the z-score reflects topology. The legacy
#' `disease_genes = "targets_disease"` path draws `T` from
#' [targets_disease_filter()], which only annotates UniProt IDs already
#' predicted as targets; `S` and `T` then overlap by construction and the
#' statistic degrades to a set-membership lookup. That path still runs but
#' emits a warning. The `n_overlap` column reports `|S ∩ T|` either way.
#'
#' @section STRINGdb is a heavy, optional dependency, not installed by default:
#' Like the `network_enrich()` Bioconductor annotation packages, `STRINGdb`
#' is kept in `Suggests`. The very first call downloads STRING's flat files
#' (protein list, aliases, interactions) for `species` under
#' `cacheDir(proj)` -- this needs internet access once; subsequent calls
#' reuse the downloaded files (STRINGdb's own on-disk caching, not
#' `patliR`-specific). This can be a genuinely large download (tens to a few
#' hundred MB for human at the full interactome). Install `STRINGdb` with
#' `BiocManager::install("STRINGdb")`.
#'
#' @section How target/gene IDs are mapped to the STRING network:
#' `STRINGdb$map()` joins the input UniProt IDs against STRING's own alias
#' table (`protein.aliases.v<version>/<species>...txt.gz`), which already
#' includes UniProt accessions as one of many alias types alongside gene
#' symbols and Ensembl IDs -- confirmed by reading `STRINGdb`'s own
#' `map()`/`get_aliases()` source rather than assumed: `map()` is a plain
#' `multi_map_df(my_data_frame, get_aliases(), id_col, "alias", "STRING_id")`
#' join, so UniProt accessions map directly with no special-casing needed.
#' UniProt IDs that fail to map (not in STRING's alias table for `species`)
#' are logged (`"network_proximity_unmapped"`) and excluded, never silently
#' dropped without a trace.
#'
#' @section Degree-preserving null model:
#' The STRING graph is first restricted to its largest connected component
#' (Guney et al. 2016; Menche et al. 2015). For each of the `n_random`
#' iterations, the compound's target set and the disease's gene set are
#' each replaced by an equal-size set of **distinct** STRING nodes drawn
#' from the same degree bins -- bins built from consecutive degree values,
#' each extended until it holds >= 100 nodes (Guney et al.'s Methods), not
#' fixed quantiles. The "closest" distance is recomputed on the resampled
#' pair. The disease side of the null is identical for every compound in a
#' `(condition, disease)`, so its `n_random` resampled sets are drawn once
#' and reused across compounds (Guney's reference implementation does the
#' same). Because each draw's `T'` is shared by every compound, the null
#' distances are computed with one multi-source breadth-first search per
#' draw (from `T'`) rather than one search per random source node -- an
#' exact reformulation (hop counts are integers), not an approximation, so
#' `d_random` is identical to the per-draw computation.
#' `z_score = (d_observed - mean(d_random)) / sd(d_random)`; a
#' strongly negative `z_score` means the compound's targets are
#' topologically closer to the disease genes than expected by chance given
#' their degree. `p_empirical` is the left-tail permutation p, and
#' `p_adjusted` its Benjamini-Hochberg value across every compound in the
#' call. `seed` (if supplied) is used only for this resampling and does not
#' leak into the caller's RNG state (same isolation pattern as
#' [network_module_robustness()]'s percolation step).
#'
#' @inheritParams network_build
#' @param disease A single `disease_id`. With `disease_genes = "disease_genes"`
#'   (the default) it must appear in
#'   `patliRResults(proj, "disease_genes")$disease_id` (i.e. it was fetched
#'   by [disease_genes_fetch()] or imported by [disease_genes_import()]).
#'   With `disease_genes = "targets_disease"` it must appear in
#'   `patliRResults(proj, "targets_disease")$disease_id` instead.
#' @param disease_genes Which slot the disease module `T` is built from.
#'   `"disease_genes"` (default): the independent, disease -> target gene set
#'   from [disease_genes_fetch()] / [disease_genes_import()]. `"targets_disease"`:
#'   the legacy path -- `T` is drawn from [targets_disease_filter()], which
#'   only ever annotates UniProt IDs already predicted as targets, so `S`
#'   and `T` overlap by construction and the z-score is biased (a warning
#'   naming this circularity is emitted; the run is not aborted). Rows
#'   record which was used in the `disease_gene_source` column.
#' @param species NCBI taxonomy ID passed to `STRINGdb$new()`. Default
#'   `9606` (human).
#' @param version STRING database version passed to `STRINGdb$new()`.
#'   Default `"12.0"`.
#' @param score_threshold Minimum STRING combined confidence score
#'   (0-1000) for an interaction to be loaded into the graph, passed to
#'   `STRINGdb$new()`. Default `400` (STRINGdb's own default -- "medium
#'   confidence").
#' @param n_random Number of degree-preserving random resamplings used to
#'   compute the null distribution for the z-score. Default `1000`.
#' @param seed `NULL` (default) or a single integer, for the random
#'   resampling. Logged either way (an auto-generated seed if `NULL`) so a
#'   run can be reproduced exactly.
#' @param store_null Logical, default `FALSE`. When `TRUE`, additionally
#'   writes every individual random draw's `d_random` value (the raw
#'   per-draw distance the loop already computes on the way to
#'   `d_random_mean`/`d_random_sd`/`z_score` -- nothing new is computed) to
#'   a separate `network_proximity_null` results slot, one row per
#'   `(compound, disease, draw)`. This is opt-in because it is large: at
#'   `n_random = 1000` and 30 compounds that is ~30,000 rows (~1 MB). It is
#'   the raw material [plot_proximity(view = "null")][plot_proximity()]
#'   needs to draw the permutation-histogram figure; when `FALSE` (default)
#'   the `network_proximity_null` slot is left untouched -- a prior
#'   `store_null = TRUE` run's stored draws are neither read nor wiped by a
#'   later `store_null = FALSE` call.
#'
#' @return The updated `proj`, with a `network_proximity` entry in
#'   [patliRResults()] (columns `condition`, `compound_id`, `disease_id`,
#'   `disease_gene_source`, `n_targets_mapped`, `n_disease_genes_mapped`,
#'   `n_overlap` (`|S ∩ T|` on the STRING-mapped sets, so the residual
#'   overlap driving `d_observed` toward 0 is always visible -- it equals
#'   `n_targets_mapped` exactly when `d_observed == 0`, i.e. `S ⊆ T`),
#'   `d_observed`, `d_random_mean`, `d_random_sd`, `z_score`, `p_empirical`,
#'   `n_random`, `seed_used`, `species`, `string_version`, `score_threshold`
#'   (the STRING interactome this run used -- recorded so [network_synergy()]
#'   can refuse to combine a `z` from one interactome with an `s_AB` from
#'   another), `n_tests_in_family` (the Benjamini-Hochberg family size for
#'   this call), `p_adjusted`), also written to
#'   `results/network_proximity.csv`. Compounds with zero targets mappable
#'   to the STRING network for `species` contribute no row (logged instead).
#'   When `store_null = TRUE`, also a `network_proximity_null` entry
#'   (columns `condition`, `compound_id`, `disease_id`,
#'   `disease_gene_source`, `draw` (integer, `1..n_random`), `d_random`
#'   (double, the resampled "closest" distance for that draw, `NA` when
#'   that draw's resampled pair had no finite path) -- `n_random` rows per
#'   compound that received a row in `network_proximity`), also written to
#'   `results/network_proximity_null.csv`.
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
#' proj <- network_proximity(
#'   proj, condition = "FLO-ET",
#'   disease = unique(patliRResults(proj, "disease_genes")$disease_id)[1]
#' ) # needs STRINGdb + internet on first call
#' patliRResults(proj, "network_proximity")
#' }
#'
#' @export
network_proximity <- function(proj, condition = NULL, disease,
                               disease_genes = c("disease_genes", "targets_disease"),
                               species = 9606, version = "12.0",
                               score_threshold = 400, n_random = 1000,
                               seed = NULL, store_null = FALSE) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(disease), length(disease) == 1, nzchar(disease))
  stopifnot(is.numeric(n_random), length(n_random) == 1, n_random >= 1)
  stopifnot(is.logical(store_null), length(store_null) == 1, !is.na(store_null))
  disease_genes <- match.arg(disease_genes)
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_proximity} needs {.pkg STRINGdb}, not installed.",
      "i" = "Install it with {.code BiocManager::install(\"STRINGdb\")}."
    ))
  }

  conditions <- .network_resolve_conditions(proj, condition)

  if (disease_genes == "disease_genes") {
    dg_all <- patliRResults(proj, "disease_genes")
    if (is.null(dg_all) || nrow(dg_all) == 0) {
      cli::cli_abort(c(
        "No {.val disease_genes} entry in {.arg proj}.",
        "i" = "Run {.fn disease_genes_fetch} or {.fn disease_genes_import} first.",
        "i" = "To use the legacy (circular) disease set from {.fn targets_disease_filter}, pass {.code disease_genes = \"targets_disease\"}."
      ))
    }
    disease_uniprot <- unique(dg_all$uniprot_id[dg_all$disease_id == disease & !is.na(dg_all$uniprot_id)])
    if (length(disease_uniprot) == 0) {
      cli::cli_abort(c(
        "No gene associated with {.val {disease}} in {.val disease_genes}.",
        "i" = "Available disease IDs: {.val {unique(dg_all$disease_id)}}.",
        "i" = "Run {.fn disease_genes_fetch} / {.fn disease_genes_import} for {.val {disease}} first."
      ))
    }
  } else {
    disease_all <- patliRResults(proj, "targets_disease")
    if (is.null(disease_all) || nrow(disease_all) == 0) {
      cli::cli_abort(c(
        "No {.val targets_disease} entry in {.arg proj}.",
        "i" = "Run {.fn targets_disease_filter} first."
      ))
    }
    disease_uniprot <- unique(disease_all$target_id[disease_all$disease_id == disease & !is.na(disease_all$target_id)])
    if (length(disease_uniprot) == 0) {
      cli::cli_abort(c(
        "No target associated with {.val {disease}} in {.val targets_disease}.",
        "i" = "Available disease IDs: {.val {unique(disease_all$disease_id)}}."
      ))
    }
    cli::cli_warn(c(
      "!" = "{.code disease_genes = \"targets_disease\"}: the disease module is drawn from the compounds' own predicted targets.",
      "i" = "{.fn targets_disease_filter} only annotates UniProt IDs already in {.val targets_imported}, so this disease set is {.strong circular} -- {.field S} and {.field T} overlap by construction, {.field d_observed} is dragged toward 0, and {.field z_score} is biased negative for reasons unrelated to topology.",
      "i" = "Run {.fn disease_genes_fetch} / {.fn disease_genes_import} and use the default {.code disease_genes = \"disease_genes\"} for an independent gene set."
    ))
  }

  used_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else as.integer(seed)
  restore_rng <- .with_seed(used_seed)
  on.exit(restore_rng(), add = TRUE)

  string_db <- .network_stringdb(proj, species, version, score_threshold)

  ## Guney et al. (2016) and Menche et al. (2015) both work on the
  ## interactome's largest connected component -- a target outside it has
  ## no finite distance to the disease module and would otherwise be
  ## silently dropped (biasing d_observed downward) or make min() warn.
  ## Extracted into .network_string_lcc() so network_synergy() gets a
  ## byte-identical graph -- s_AB and z must be on the same interactome.
  lcc <- .network_string_lcc(proj, species, version, score_threshold, string_db = string_db)
  g <- lcc$graph
  degree_all <- lcc$degree
  bins <- lcc$bins
  ## node_names / bin_of_node are hoisted out of
  ## .network_resample_matched() -- it is called ~2 * n_random per
  ## compound and both are invariant across every one of those calls.
  node_names <- names(degree_all)
  bin_of_node <- lcc$bin_of_node
  proj <- .log_append(
    proj, step = "network_proximity", id = NA_character_,
    message = paste0("STRING interactome restricted to its largest connected component (",
                      igraph::vcount(g), " nodes)")
  )

  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions
  ## Only populated when store_null = TRUE -- one element per condition,
  ## each the rbind of that condition's per-compound raw-draw rows.
  null_rows <- vector("list", length(conditions))
  names(null_rows) <- conditions

  for (cond in conditions) {
    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    compounds <- unique(ct$compound_id)
    target_sets <- split(ct$uniprot_id, ct$compound_id)

    all_uniprot <- unique(c(unlist(target_sets), disease_uniprot))
    map_df <- string_db$map(data.frame(uniprot_id = all_uniprot, stringsAsFactors = FALSE),
                             "uniprot_id", removeUnmappedRows = FALSE, quiet = TRUE)
    unmapped <- map_df$uniprot_id[is.na(map_df$STRING_id)]
    if (length(unmapped) > 0) {
      proj <- .log_append(
        proj, step = "network_proximity", id = NA_character_,
        message = paste0(
          "network_proximity_unmapped: ", length(unmapped), " of ", nrow(map_df),
          " UniProt ID(s) had no STRING_id for this species and were excluded (e.g. ",
          paste(utils::head(unmapped, 10), collapse = ", "), ")"
        )
      )
    }
    uni_to_string <- stats::setNames(map_df$STRING_id, map_df$uniprot_id)

    disease_string <- unique(stats::na.omit(uni_to_string[disease_uniprot]))
    if (length(disease_string) == 0) {
      slot_lab <- if (disease_genes == "disease_genes") "disease_genes" else "targets_disease"
      cli::cli_abort(c(
        "None of the {.val {disease}} gene(s) could be mapped to the STRING network.",
        "i" = "Check {.arg species} matches the organism the {.val {slot_lab}} set was built for."
      ))
    }

    ## The disease module T -- and therefore the degree-matched disease-side
    ## null -- is identical for every compound in this (condition, disease),
    ## so both the T -> STRING/LCC filter and the n_random null draws are
    ## hoisted out of the per-compound loop (spec 1.8 [SHOULD]).
    g_names <- igraph::V(g)$name
    target_string <- disease_string[disease_string %in% g_names]
    ## The null draws are held as integer vertex indices (node_names is
    ## V(g)$name in vertex order), not names: the batched distance step
    ## below indexes the graph directly, and a match() against ~17k names
    ## per draw would otherwise cost more than the draw itself.
    t_rand_list <- if (length(target_string) > 0) {
      target_idx <- match(target_string, node_names)
      lapply(seq_len(n_random), function(j)
        .network_resample_matched_idx(target_idx, bins, bin_of_node))
    } else {
      vector("list", n_random)
    }

    ## Pass 1 -- resolve every compound's S and draw its n_random
    ## degree-matched S' sets, in the same compound-then-draw order the
    ## draws were always taken in (nothing else in this loop touches the
    ## RNG), so a given seed yields exactly the same random sets as the
    ## historical one-draw-then-one-distances() loop.
    source_by_cp <- vector("list", length(compounds))
    s_rand_by_cp <- vector("list", length(compounds))
    for (i in seq_along(compounds)) {
      cp <- compounds[i]
      source_string <- unique(stats::na.omit(uni_to_string[target_sets[[cp]]]))
      source_string <- source_string[source_string %in% g_names]

      if (length(source_string) == 0 || length(target_string) == 0) {
        proj <- .log_append(
          proj, step = "network_proximity", id = cp,
          message = paste0("condition '", cond, "': no mappable/graph-present target(s) for this compound or the disease gene set; skipped")
        )
        next
      }
      source_by_cp[[i]] <- source_string
      source_idx <- match(source_string, node_names)
      s_rand_by_cp[[i]] <- lapply(seq_len(n_random), function(j)
        .network_resample_matched_idx(source_idx, bins, bin_of_node))
    }
    scored <- which(!vapply(source_by_cp, is.null, logical(1)))

    ## Pass 2 -- every compound's null distances at once: one multi-source
    ## BFS per draw (from that draw's T'), instead of |S'| BFS runs per
    ## (compound, draw). See .network_closest_distance_null().
    d_random_by_cp <- vector("list", length(compounds))
    if (length(scored) > 0) {
      d_random_by_cp[scored] <- .network_closest_distance_null(g, s_rand_by_cp[scored], t_rand_list)
    }

    cond_rows <- vector("list", length(compounds))
    cond_null_rows <- vector("list", length(compounds))
    for (i in scored) {
      cp <- compounds[i]
      source_string <- source_by_cp[[i]]

      n_overlap <- length(intersect(unique(source_string), unique(target_string)))
      d_observed <- .network_closest_distance(g, source_string, target_string)
      ## The raw, un-filtered per-draw distances -- kept before the
      ## is.finite() filter below, so store_null = TRUE can persist exactly
      ## what fed into d_random_mean/d_random_sd/z_score without changing
      ## how any of those are computed.
      d_random_raw <- d_random_by_cp[[i]]
      d_random <- d_random_raw[is.finite(d_random_raw)]

      if (store_null) {
        cond_null_rows[[i]] <- data.frame(
          condition = cond, compound_id = cp, disease_id = disease,
          disease_gene_source = disease_genes,
          draw = seq_len(n_random),
          d_random = ifelse(is.finite(d_random_raw), d_random_raw, NA_real_),
          stringsAsFactors = FALSE
        )
      }

      d_random_mean <- if (length(d_random) > 0) mean(d_random) else NA_real_
      d_random_sd <- if (length(d_random) > 1) stats::sd(d_random) else NA_real_
      z_score <- if (!is.na(d_random_sd) && d_random_sd > 0) {
        (d_observed - d_random_mean) / d_random_sd
      } else {
        NA_real_
      }
      ## left-tail permutation p: how often is a random module at least as
      ## close as the observed one. (1 + .) / (n + 1) is the standard
      ## never-zero estimator (Phipson & Smyth 2010).
      p_empirical <- if (length(d_random) > 0) {
        (1 + sum(d_random <= d_observed)) / (length(d_random) + 1)
      } else {
        NA_real_
      }

      cond_rows[[i]] <- data.frame(
        condition = cond, compound_id = cp, disease_id = disease,
        disease_gene_source = disease_genes,
        n_targets_mapped = length(source_string), n_disease_genes_mapped = length(target_string),
        n_overlap = n_overlap,
        d_observed = d_observed, d_random_mean = d_random_mean, d_random_sd = d_random_sd,
        z_score = z_score, p_empirical = p_empirical, n_random = length(d_random),
        seed_used = used_seed,
        species = as.numeric(species), string_version = as.character(version),
        score_threshold = as.numeric(score_threshold),
        n_tests_in_family = NA_integer_, p_adjusted = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    cond_rows <- cond_rows[!vapply(cond_rows, is.null, logical(1))]
    rows[[cond]] <- if (length(cond_rows) == 0) .empty_network_proximity_row() else do.call(rbind, cond_rows)
    if (store_null) {
      cond_null_rows <- cond_null_rows[!vapply(cond_null_rows, is.null, logical(1))]
      null_rows[[cond]] <- if (length(cond_null_rows) == 0) {
        .empty_network_proximity_null_row()
      } else {
        do.call(rbind, cond_null_rows)
      }
    }
    proj <- .log_append(
      proj, step = "network_proximity", id = NA_character_,
      message = paste0("condition '", cond, "', disease '", disease, "': proximity computed for ",
                        nrow(rows[[cond]]), " compound(s), seed ", used_seed)
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  ## p_adjusted / n_tests_in_family are part of the schema even for a
  ## zero-row result, so the CSV stays column-stable and .network_upsert()'s
  ## rbind never faces a column-count mismatch on a later zero-row rerun.
  ## BH across every (condition, compound) proximity test in this call; the
  ## family size is recorded so network_synergy() can detect a merged table
  ## whose rows came from calls of different family sizes and warn when the
  ## alpha gate is arithmetically unreachable (p_floor * m / k).
  if (nrow(result) > 0) {
    result$n_tests_in_family <- nrow(result)
    result$p_adjusted <- stats::p.adjust(result$p_empirical, method = "BH")
  }
  ## Recompute the whole condition/disease, including compounds that lost
  ## every edge and therefore no longer occur in the current edge table.
  touched_keys <- data.frame(condition = conditions, disease_id = disease, stringsAsFactors = FALSE)
  result <- .network_upsert(
    proj, "network_proximity", result,
    c("condition", "disease_id"), touched_keys = touched_keys
  )

  patliRResults(proj, "network_proximity") <- result

  ## store_null = FALSE (default): network_proximity_null is not read,
  ## written, or wiped -- any draws a previous store_null = TRUE call left
  ## there are untouched by this call.
  if (store_null) {
    null_result <- do.call(rbind, null_rows)
    rownames(null_result) <- NULL
    ## Re-running for this (condition, disease) replaces only that pair's
    ## stored draws -- keyed coarser than the main table (no compound_id)
    ## because every draw for every compound in scope is being recomputed
    ## together, unlike network_proximity's per-compound touched_keys.
    null_touched_keys <- data.frame(condition = conditions, disease_id = disease, stringsAsFactors = FALSE)
    null_result <- .network_upsert(
      proj, "network_proximity_null", null_result,
      c("condition", "disease_id"), touched_keys = null_touched_keys
    )
    patliRResults(proj, "network_proximity_null") <- null_result
    .write_results_csv(proj, "network_proximity_null", null_result)
  }

  .write_results_csv(proj, "network_proximity", result)
  .write_log_csv(proj)
  proj
}

#' Construct a `STRINGdb` object, caching its downloaded flat files under
#' `cacheDir(proj)` (STRINGdb's own on-disk caching, not `patliR`-specific
#' -- deleting `cacheDir(proj)/stringdb` just triggers a redownload, same
#' "safe to delete" cache principle as `.network_cache_path()`).
#' @return A `STRINGdb` reference-class instance.
#' @keywords internal
.network_stringdb <- function(proj, species, version, score_threshold) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  STRINGdb <- get("STRINGdb", envir = asNamespace("STRINGdb"))
  STRINGdb$new(version = version, species = species,
               score_threshold = score_threshold, input_directory = dir)
}

#' Contiguous-value bins over a numeric node property -- the pool each
#' random resampling draws from, per node, to preserve that property.
#'
#' @details
#' Generic over any per-node numeric value: [network_proximity()] passes
#' node degree (Guney et al. 2016, *Nat Commun* 7:10331 -- the
#' degree-preserving null), [network_degeneracy()] passes each gene's GO
#' annotation count. Bins are built from *consecutive values*, a bin
#' extended until it holds at least `min_per_bin` (default 100) nodes. This
#' keeps bins fine-grained where the values are dense and only coarsens in
#' a sparse tail -- unlike fixed quantile bins, where on a heavy-tailed
#' distribution a single decile spans an enormous range and a
#' high-value node can be swapped for a low-value one, destroying the
#' property preservation exactly where the bias it corrects for is
#' strongest.
#' @param values Named or unnamed numeric vector, one entry per node, in
#'   the node order the caller will align `bins` against.
#' @return A list of integer vectors (indices into `values`), one per bin.
#' @keywords internal
.network_value_bins <- function(values, min_per_bin = 100) {
  ## one sorted pass: rle() over the sorted values gives every unique
  ## value and its multiplicity at once, instead of an O(|V|) rescan per
  ## unique value.
  runs <- rle(sort(as.vector(values)))
  uniq_vals <- runs$values
  count_per_value <- runs$lengths

  ## walk up the unique values, closing a bin once it has >= min_per_bin
  bin_of_value <- integer(length(uniq_vals))
  b <- 1L
  count <- 0L
  for (i in seq_along(uniq_vals)) {
    bin_of_value[i] <- b
    count <- count + count_per_value[i]
    if (count >= min_per_bin) { b <- b + 1L; count <- 0L }
  }
  ## a trailing under-full bin is merged into the previous one
  if (count > 0 && b > 1L) bin_of_value[bin_of_value == b] <- b - 1L

  node_bin <- bin_of_value[match(values, uniq_vals)]
  split(seq_along(values), node_bin)
}

#' Replace `string_ids` with an equal-size set of *distinct* nodes drawn
#' from the same value bins ([.network_value_bins()]).
#'
#' @details
#' Nodes are grouped by their bin and, per bin, `k` distinct replacements
#' are drawn from that bin's pool -- so the returned set has exactly
#' `length(string_ids)` distinct nodes, matching the property-preserving
#' null (Guney et al.'s "a set of |S| proteins with matching degrees" for
#' [network_proximity()]; an annotation-count-matched gene set for
#' [network_degeneracy()]). Drawing each node independently with
#' replacement (and then `unique()`-ing downstream) would shrink the
#' random set, inflating its variance and shrinking `|z_score|` in a way
#' that varies with set size.
#'
#' `node_names` (the node identifiers, in the order `bins` indexes) and
#' `bin_of_node` (the bin index of every node, aligned to `node_names`)
#' are passed in rather than rebuilt here: this function is called
#' `~2 * n_random` times per pair and both are invariant across every call.
#' `string_ids` must all be in `node_names` (the resampling pool) -- an
#' absent id would silently corrupt the output.
#' @return Character vector of node IDs, length `length(string_ids)`.
#' @keywords internal
.network_resample_matched <- function(string_ids, node_names, bins, bin_of_node) {
  stopifnot(all(string_ids %in% node_names))
  node_names[.network_resample_matched_idx(match(string_ids, node_names), bins, bin_of_node)]
}

#' Index form of [.network_resample_matched()]
#'
#' @description
#' Same draw, same RNG consumption (one `sample.int()` per occupied bin, in
#' ascending bin order), but takes and returns integer positions into the
#' node vector `bins` indexes rather than node names -- so a caller making
#' thousands of draws against a ~17k-node graph (the [network_proximity()]
#' null) pays the name -> index `match()` once, not once per draw. For
#' index `idx` the two forms return `node_names[out]` and `out` for the
#' same random stream.
#' @param idx Integer positions (into the node vector `bins` /
#'   `bin_of_node` are aligned to) of the nodes to replace.
#' @return Integer vector of node positions, length `length(idx)`.
#' @keywords internal
.network_resample_matched_idx <- function(idx, bins, bin_of_node) {
  by_bin <- split(seq_along(idx), bin_of_node[idx])

  out <- integer(length(idx))
  for (bs in names(by_bin)) {
    slots <- by_bin[[bs]]
    pool <- bins[[as.integer(bs)]]
    k <- length(slots)
    ## distinct draw; if the bin is smaller than k (rare, only in the
    ## coarse tail) fall back to with-replacement for that bin.
    picks <- if (k <= length(pool)) {
      pool[sample.int(length(pool), k)]
    } else {
      pool[sample.int(length(pool), k, replace = TRUE)]
    }
    out[slots] <- picks
  }
  if (any(out == 0L)) {
    cli::cli_abort("{.fn .network_resample_matched}: {sum(out == 0L)} output slot(s) unfilled -- a `string_ids` entry has no bin.")
  }
  out
}

#' Guney et al. (2016) "closest" distance between two node sets on `g`
#'
#' @description
#' `unique()` on both `source_ids` and `target_ids` before calling
#' `igraph::distances()`: igraph's underlying C routine rejects a `to=`
#' argument with duplicate vertices. `.network_resample_matched()`
#' returns distinct nodes, but the observed sets can still contain a
#' repeat, and the "closest" measure only cares about the node *set*, so
#' deduping is harmless.
#'
#' @return A single numeric (unweighted shortest-path hops, `weights = NA`
#'   -- same rationale as [network_centrality()]: STRING's `combined_score`
#'   is a confidence, not a distance).
#' @keywords internal
.network_closest_distance <- function(g, source_ids, target_ids) {
  d <- igraph::distances(g, v = unique(source_ids), to = unique(target_ids), weights = NA)
  d[!is.finite(d)] <- NA_real_
  row_min <- apply(d, 1, min, na.rm = TRUE)
  mean(row_min[is.finite(row_min)])
}

#' Guney "closest" distance for every (compound, draw) of the proximity null
#'
#' @description
#' Batched, exact equivalent of calling [.network_closest_distance()] once
#' per random draw -- `.network_closest_distance(g, V(g)$name[s_rand[[k]][[j]]],
#' V(g)$name[t_rand[[j]]])` for every compound `k` and draw `j` -- which
#' costs `|S'|` breadth-first searches per `(k, j)` and dominated
#' [network_proximity()]'s run time on a real interactome (tens of
#' thousands of BFS runs per disease).
#'
#' @details
#' The "closest" measure only needs, for each source node `s`, the hop
#' distance from `s` to the *nearest* member of the draw's target set
#' `T'_j`: `min_{t in T'_j} d(s, t)`. That is a multi-source BFS from
#' `T'_j`, which gives the value for every node at once. It is run as a
#' single-source BFS from a virtual super-node: `g` is copied as a
#' directed graph (every undirected edge in both directions, matching
#' `distances()`'s default `mode = "all"`), and one extra vertex per draw
#' gets out-edges to that draw's `T'_j` and no in-edges -- so no path can
#' pass through another draw's super-node -- and `distances(mode = "out")`
#' from it minus 1 is `min_{t in T'_j} d(s, t)`, `Inf` where unreachable.
#' The random `T'_j` is shared by every compound (drawn once per
#' `(condition, disease)`), so the cost falls from
#' `sum_k |S_k| * n_random` BFS runs to `n_random`.
#'
#' Hop counts are integers, so the per-source minima are identical to the
#' old per-draw `distances()` matrix row minima; each draw's value is then
#' `mean()` over its unique, finite per-source minima -- the same
#' reduction, in the same function, as [.network_closest_distance()]
#' (a source with no finite path to `T'_j` is dropped; `NaN` if none has
#' one) -- so `d_random`, and every statistic built on it, is unchanged.
#'
#' Memory: the draws are processed in chunks so that the
#' `draws x needed-source-nodes` distance block stays at or below
#' `max_cells` doubles (default `1e7`, ~80 MB); only the columns of nodes
#' that actually occur in some `S'` are materialised.
#'
#' @param g The LCC `igraph` (undirected, unweighted hops).
#' @param s_rand List over compounds; each a list over the `n_random`
#'   draws of integer vertex indices (the random `S'`).
#' @param t_rand List over the `n_random` draws of integer vertex indices
#'   (the random `T'`, shared by every compound).
#' @param max_cells Upper bound on the cells of one chunk's distance block.
#' @return List parallel to `s_rand`, each a numeric vector of length
#'   `length(t_rand)` (the raw per-draw `d_random`, non-finite where the
#'   draw had no finite path).
#' @keywords internal
.network_closest_distance_null <- function(g, s_rand, t_rand, max_cells = 1e7) {
  n <- igraph::vcount(g)
  n_draws <- length(t_rand)
  out <- lapply(s_rand, function(x) rep(NA_real_, n_draws))
  if (n_draws == 0L || length(s_rand) == 0L) return(out)

  el <- igraph::as_edgelist(g, names = FALSE)
  base_edges <- as.vector(rbind(c(el[, 1], el[, 2]), c(el[, 2], el[, 1])))
  rm(el)

  ## only the columns of nodes that ever occur in some S' are needed
  needed <- sort(unique(unlist(s_rand, use.names = FALSE)))
  col_of <- integer(n)
  col_of[needed] <- seq_along(needed)
  ## unique() as in .network_closest_distance(): a with-replacement
  ## fallback draw can repeat a node, and the measure is over the node set.
  s_cols <- lapply(s_rand, function(draws) lapply(draws, function(s) col_of[unique(s)]))

  per_chunk <- as.integer(max(1, min(n_draws, floor(max_cells / length(needed)))))
  for (first in seq(1L, n_draws, by = per_chunk)) {
    js <- first:min(n_draws, first + per_chunk - 1L)
    super <- n + seq_along(js)
    t_sets <- lapply(t_rand[js], unique)
    super_edges <- as.vector(rbind(rep(super, lengths(t_sets)), unlist(t_sets, use.names = FALSE)))
    g2 <- igraph::make_graph(c(base_edges, super_edges), n = n + length(js), directed = TRUE)
    D <- igraph::distances(g2, v = super, to = needed, mode = "out", weights = NA)
    D <- D - 1
    rm(g2)
    for (k in seq_along(s_rand)) {
      out[[k]][js] <- vapply(seq_along(js), function(r) {
        d <- D[r, s_cols[[k]][[js[r]]]]
        mean(d[is.finite(d)])
      }, numeric(1))
    }
    rm(D)
  }
  out
}

#' @keywords internal
.empty_network_proximity_row <- function() {
  data.frame(
    condition = character(0), compound_id = character(0), disease_id = character(0),
    disease_gene_source = character(0),
    n_targets_mapped = integer(0), n_disease_genes_mapped = integer(0),
    n_overlap = integer(0),
    d_observed = double(0), d_random_mean = double(0), d_random_sd = double(0),
    z_score = double(0), p_empirical = double(0), n_random = integer(0),
    seed_used = integer(0),
    species = double(0), string_version = character(0), score_threshold = double(0),
    n_tests_in_family = integer(0), p_adjusted = double(0),
    stringsAsFactors = FALSE
  )
}

#' Zero-row constructor for the `network_proximity_null` slot (Phase-0
#' zero-row invariant -- see `.empty_network_proximity_row()` above)
#' @keywords internal
.empty_network_proximity_null_row <- function() {
  data.frame(
    condition = character(0), compound_id = character(0), disease_id = character(0),
    disease_gene_source = character(0),
    draw = integer(0), d_random = double(0),
    stringsAsFactors = FALSE
  )
}

#' STRING interactome largest connected component, with its degree bins
#'
#' @description
#' Factored out of [network_proximity()] so that it and [network_synergy()]
#' operate on a byte-identical graph -- the Menche (2015) separation `s_AB`
#' and the Guney (2016) proximity `z` are only comparable (Cheng et al.
#' 2019) when computed on the same interactome, at the same `species`,
#' `version` and `score_threshold`.
#'
#' Guney et al. (2016) and Menche et al. (2015) both restrict the STRING
#' graph to its largest connected component: a node outside it has no
#' finite distance to any module and would otherwise bias distances
#' downward or make `min()` warn.
#'
#' @details
#' Only the graph is cached, under
#' `cacheDir(proj)/stringdb/lcc_<species>_<version>_<threshold>.rds`. The
#' degree, the degree bins and the per-node bin index are recomputed from
#' it on every load -- a future `min_per_bin` parameter would otherwise be
#' served a stale cached binning. Like `.network_stringdb()`'s downloaded
#' flat files, this `.rds` is **safe to delete**: STRING's own flat files
#' under the same directory are the real cache, and deleting the `.rds`
#' just triggers one recompute of the component decomposition.
#'
#' @param string_db Optional pre-constructed `STRINGdb` instance (avoids a
#'   second `STRINGdb$new()` when the caller already built one); only used
#'   on a cache miss.
#' @return `list(graph, degree, bins, bin_of_node)` -- `graph` the LCC
#'   `igraph`, `degree` its named degree vector, `bins` the list of
#'   node-index vectors from [.network_value_bins()], `bin_of_node` the
#'   bin index of every node aligned to `names(degree)`.
#' @keywords internal
.network_string_lcc <- function(proj, species, version, score_threshold, string_db = NULL) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  cache_rds <- file.path(dir, paste0("lcc_", species, "_", version, "_", score_threshold, ".rds"))

  if (file.exists(cache_rds)) {
    g <- readRDS(cache_rds)
  } else {
    if (is.null(string_db)) string_db <- .network_stringdb(proj, species, version, score_threshold)
    g0 <- string_db$get_graph()
    comp <- igraph::components(g0)
    lcc_id <- which.max(comp$csize)
    g <- if (sum(comp$membership != lcc_id) > 0) {
      igraph::induced_subgraph(g0, which(comp$membership == lcc_id))
    } else {
      g0
    }
    saveRDS(g, cache_rds)
  }

  degree_all <- igraph::degree(g)
  bins <- .network_value_bins(degree_all)
  bin_of_node <- integer(length(degree_all))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b
  list(graph = g, degree = degree_all, bins = bins, bin_of_node = bin_of_node)
}
