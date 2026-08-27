#' @include AllGenerics.R internal.R network_build.R network_module_robustness.R
NULL

## network_proximity() -- network-based compound-disease proximity on the
## STRING protein-protein interaction network. STRINGdb chosen as the
## background interactome, confirmed with Uriel (chat, 2026-07-23): a
## per-condition compound-target *bipartite* graph (network_build()) has no
## real topological distance between two targets, so answering "how close
## is this compound's target set to this disease's gene set on the human
## interactome" needs a genuinely separate, external PPI network.
##
## Method: Guney, E., Menche, J., Vidal, M., & Barabasi, A.L. (2016).
## "Network-based in silico drug efficacy screening." Nature Communications,
## 7, 10331. https://doi.org/10.1038/ncomms10331 -- the "closest" proximity
## measure (their Eq. 2): d(S,T) = (1/|S|) * sum_{s in S} min_{t in T} d(s,t),
## z-scored against a degree-preserving null model (both S and T resampled
## from same-degree-bin nodes, n_random times), z = (d_obs - mean(d_rand)) /
## sd(d_rand). More negative z = closer than random chance = topologically
## plausible drug-disease relationship.
##
## STRINGdb API verified against the live Bioconductor docs (rdrr.io man
## pages + source of rstring.R, 2026-07-23) before writing this, per this
## project's "verify, don't assume" discipline -- see the roxygen sections
## below for exactly what was checked and why.

#' Network-based compound-disease proximity on the STRING interactome
#'
#' @description
#' For each compound present in a condition's network ([network_build()]),
#' computes the topological "closest" distance (Guney et al. 2016, *Nat
#' Commun* 7:10331) between that compound's target set and a disease's
#' associated gene set (from [targets_disease_filter()]), on the STRING
#' protein-protein interaction network (`STRINGdb`), and z-scores it against
#' a degree-preserving random null model.
#'
#' @section STRINGdb is a heavy, optional dependency, not installed by default:
#' Like the `network_enrich()` Bioconductor annotation packages, `STRINGdb`
#' is kept in `Suggests`. The very first call downloads STRING's flat files
#' (protein list, aliases, interactions) for `species` under
#' `cacheDir(proj)` -- this needs internet access once; subsequent calls
#' reuse the downloaded files (STRINGdb's own on-disk caching, not
#' `patliR`-specific). This can be a genuinely large download (tens to a few
#' hundred MB for human at the full interactome) -- see `TESTING_GUIDE.Rmd`,
#' section 0.
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
#' For each of the `n_random` iterations, every protein in the compound's
#' target set *and* every protein in the disease's gene set is independently
#' replaced by a randomly chosen STRING node from the same degree decile bin
#' (an approximation of Guney et al.'s degree-preserving resampling), and the
#' "closest" distance is recomputed on this resampled pair of sets. `z_score
#' = (d_observed - mean(d_random)) / sd(d_random)`; a strongly negative
#' `z_score` means the compound's targets are topologically closer to the
#' disease genes than expected by chance given their degree. `seed` (if
#' supplied) is used only for this resampling and does not leak into the
#' caller's RNG state (same isolation pattern as
#' [network_module_robustness()]'s percolation step).
#'
#' @inheritParams network_build
#' @param disease A single `disease_id` as it appears in
#'   `patliRResults(proj, "targets_disease")$disease_id` (i.e. already
#'   resolved by [targets_disease_filter()] for at least one target in this
#'   project).
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
#'
#' @return The updated `proj`, with a `network_proximity` entry in
#'   [patliRResults()] (columns `condition`, `compound_id`, `disease_id`,
#'   `n_targets_mapped`, `n_disease_genes_mapped`, `d_observed`,
#'   `d_random_mean`, `d_random_sd`, `z_score`, `n_random`, `seed_used`),
#'   also written to `results/network_proximity.csv`. Compounds with zero
#'   targets mappable to the STRING network for `species` contribute no row
#'   (logged instead).
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
#' proj <- targets_disease_filter(proj, disease = "type 2 diabetes mellitus")
#' proj <- network_proximity(
#'   proj, condition = "FLO-ET",
#'   disease = unique(patliRResults(proj, "targets_disease")$disease_id)[1]
#' ) # needs STRINGdb + internet on first call
#' patliRResults(proj, "network_proximity")
#' }
#'
#' @export
network_proximity <- function(proj, condition = NULL, disease,
                               species = 9606, version = "12.0",
                               score_threshold = 400, n_random = 1000,
                               seed = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(disease), length(disease) == 1, nzchar(disease))
  stopifnot(is.numeric(n_random), length(n_random) == 1, n_random >= 1)
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_proximity} needs {.pkg STRINGdb}, not installed.",
      "i" = "See {.file TESTING_GUIDE.Rmd}, section 0, for the {.fn BiocManager::install} chunk."
    ))
  }

  conditions <- .network_resolve_conditions(proj, condition)

  disease_all <- patliRResults(proj, "targets_disease")
  if (is.null(disease_all) || nrow(disease_all) == 0) {
    cli::cli_abort(c(
      "No {.val targets_disease} entry in {.arg proj}.",
      "i" = "Run {.fn targets_disease_filter} first."
    ))
  }
  disease_uniprot <- unique(disease_all$target_id[disease_all$disease_id == disease])
  if (length(disease_uniprot) == 0) {
    cli::cli_abort(c(
      "No target associated with {.val {disease}} in {.val targets_disease}.",
      "i" = "Available disease IDs: {.val {unique(disease_all$disease_id)}}."
    ))
  }

  used_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1) else as.integer(seed)
  restore_rng <- .with_seed(used_seed)
  on.exit(restore_rng(), add = TRUE)

  string_db <- .network_stringdb(proj, species, version, score_threshold)
  g <- string_db$get_graph()
  degree_all <- igraph::degree(g)
  bins <- .network_degree_bins(degree_all)

  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions

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
        proj, step = "network_proximity", id = unmapped,
        message = "network_proximity_unmapped: no STRING_id found for this UniProt ID in this species; excluded"
      )
    }
    uni_to_string <- stats::setNames(map_df$STRING_id, map_df$uniprot_id)

    disease_string <- unique(stats::na.omit(uni_to_string[disease_uniprot]))
    if (length(disease_string) == 0) {
      cli::cli_abort(c(
        "None of the {.val {disease}} target(s) could be mapped to the STRING network.",
        "i" = "Check {.arg species} matches the organism of {.val targets_disease}."
      ))
    }

    cond_rows <- vector("list", length(compounds))
    for (i in seq_along(compounds)) {
      cp <- compounds[i]
      source_string <- unique(stats::na.omit(uni_to_string[target_sets[[cp]]]))
      source_string <- source_string[source_string %in% igraph::V(g)$name]
      target_string <- disease_string[disease_string %in% igraph::V(g)$name]

      if (length(source_string) == 0 || length(target_string) == 0) {
        proj <- .log_append(
          proj, step = "network_proximity", id = cp,
          message = paste0("condition '", cond, "': no mappable/graph-present target(s) for this compound or the disease gene set; skipped")
        )
        cond_rows[[i]] <- NULL
        next
      }

      d_observed <- .network_closest_distance(g, source_string, target_string)
      d_random <- vapply(seq_len(n_random), function(j) {
        s_rand <- .network_resample_degree_matched(source_string, degree_all, bins)
        t_rand <- .network_resample_degree_matched(target_string, degree_all, bins)
        .network_closest_distance(g, s_rand, t_rand)
      }, numeric(1))
      d_random <- d_random[is.finite(d_random)]

      d_random_mean <- if (length(d_random) > 0) mean(d_random) else NA_real_
      d_random_sd <- if (length(d_random) > 1) stats::sd(d_random) else NA_real_
      z_score <- if (!is.na(d_random_sd) && d_random_sd > 0) {
        (d_observed - d_random_mean) / d_random_sd
      } else {
        NA_real_
      }

      cond_rows[[i]] <- data.frame(
        condition = cond, compound_id = cp, disease_id = disease,
        n_targets_mapped = length(source_string), n_disease_genes_mapped = length(target_string),
        d_observed = d_observed, d_random_mean = d_random_mean, d_random_sd = d_random_sd,
        z_score = z_score, n_random = length(d_random), seed_used = used_seed,
        stringsAsFactors = FALSE
      )
    }

    cond_rows <- cond_rows[!vapply(cond_rows, is.null, logical(1))]
    rows[[cond]] <- if (length(cond_rows) == 0) .empty_network_proximity_row() else do.call(rbind, cond_rows)
    proj <- .log_append(
      proj, step = "network_proximity", id = NA_character_,
      message = paste0("condition '", cond, "', disease '", disease, "': proximity computed for ",
                        nrow(rows[[cond]]), " compound(s), seed ", used_seed)
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_proximity", result, c("condition", "disease_id", "compound_id"))

  patliRResults(proj, "network_proximity") <- result
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

#' Degree decile bins over every node of the interactome -- the pool each
#' random resampling draws from, per node, to approximately preserve degree.
#' @return A list of integer vectors (node indices), one per bin.
#' @keywords internal
.network_degree_bins <- function(degree_all, n_bins = 10) {
  breaks <- stats::quantile(degree_all, probs = seq(0, 1, length.out = n_bins + 1), type = 1)
  breaks[1] <- breaks[1] - 1 ## make the lowest bin inclusive of the minimum
  bin_id <- cut(degree_all, breaks = unique(breaks), include.lowest = TRUE, labels = FALSE)
  split(seq_along(degree_all), bin_id)
}

#' Replace each node in `string_ids` with a random node from the same
#' degree bin (with replacement across nodes, i.e. two input nodes in the
#' same bin may draw the same replacement).
#' @return Character vector of STRING IDs, same length as `string_ids`.
#' @keywords internal
.network_resample_degree_matched <- function(string_ids, degree_all, bins) {
  node_names <- names(degree_all)
  bin_of_node <- integer(length(degree_all))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b

  vapply(string_ids, function(sid) {
    idx <- match(sid, node_names)
    b <- bin_of_node[idx]
    pool <- bins[[b]]
    ## sample(pool, 1) would misbehave if length(pool) == 1 -- sample()
    ## treats a length-1 numeric x as "sample from 1:x", not as "the pool
    ## is this one element" (a well-known base-R footgun). sample.int() on
    ## the pool's length, then indexing into pool, avoids that entirely.
    node_names[pool[sample.int(length(pool), 1)]]
  }, character(1), USE.NAMES = FALSE)
}

#' Guney et al. (2016) "closest" distance between two node sets on `g`
#'
#' @description
#' `unique()` on both `source_ids` and `target_ids` before calling
#' `igraph::distances()`: confirmed via a real error from Uriel's first
#' STRINGdb run (2026-07-23) that igraph's underlying C routine
#' (`vendor/cigraph/src/paths/unweighted.c`) rejects a `to=` argument with
#' duplicate vertices outright ("Target vertex list must not have any
#' duplicates"). Duplicates are expected here, not a caller bug:
#' `.network_resample_degree_matched()` draws each node independently *with
#' replacement* from its degree bin, so two different input nodes
#' frequently resample to the same replacement, especially for small
#' target sets (a handful of disease genes). Since the "closest" measure
#' only cares about the *set* of source/target nodes, deduping changes
#' nothing about the result -- it only removes a redundant row/column
#' `igraph::distances()` would otherwise refuse to compute.
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

#' @keywords internal
.empty_network_proximity_row <- function() {
  data.frame(
    condition = character(0), compound_id = character(0), disease_id = character(0),
    n_targets_mapped = integer(0), n_disease_genes_mapped = integer(0),
    d_observed = double(0), d_random_mean = double(0), d_random_sd = double(0),
    z_score = double(0), n_random = integer(0), seed_used = integer(0),
    stringsAsFactors = FALSE
  )
}
