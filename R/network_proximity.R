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
#' pair. `z_score = (d_observed - mean(d_random)) / sd(d_random)`; a
#' strongly negative `z_score` means the compound's targets are
#' topologically closer to the disease genes than expected by chance given
#' their degree. `p_empirical` is the left-tail permutation p, and
#' `p_adjusted` its Benjamini-Hochberg value across every compound in the
#' call. `seed` (if supplied) is used only for this resampling and does not
#' leak into the caller's RNG state (same isolation pattern as
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
#'   `d_random_mean`, `d_random_sd`, `z_score`, `p_empirical`,
#'   `p_adjusted`, `n_random`, `seed_used`), also written to
#'   `results/network_proximity.csv`. Compounds with zero targets mappable
#'   to the STRING network for `species` contribute no row (logged instead).
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
      "i" = "Install it with {.code BiocManager::install(\"STRINGdb\")}."
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

  ## Guney et al. (2016) and Menche et al. (2015) both work on the
  ## interactome's largest connected component -- a target outside it has
  ## no finite distance to the disease module and would otherwise be
  ## silently dropped (biasing d_observed downward) or make min() warn.
  comp <- igraph::components(g)
  lcc <- which.max(comp$csize)
  n_dropped <- sum(comp$membership != lcc)
  if (n_dropped > 0) {
    g <- igraph::induced_subgraph(g, which(comp$membership == lcc))
    proj <- .log_append(
      proj, step = "network_proximity", id = NA_character_,
      message = paste0("restricted STRING network to its largest connected component (",
                        igraph::vcount(g), " nodes; ", n_dropped, " dropped)")
    )
  }

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
        n_targets_mapped = length(source_string), n_disease_genes_mapped = length(target_string),
        d_observed = d_observed, d_random_mean = d_random_mean, d_random_sd = d_random_sd,
        z_score = z_score, p_empirical = p_empirical, n_random = length(d_random),
        seed_used = used_seed, stringsAsFactors = FALSE
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
  ## BH across every (condition, compound) proximity test in this call
  if (nrow(result) > 0 && "p_empirical" %in% names(result)) {
    result$p_adjusted <- stats::p.adjust(result$p_empirical, method = "BH")
  }
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

#' Degree bins over every node of the interactome -- the pool each random
#' resampling draws from, per node, to preserve degree.
#'
#' @details
#' Guney et al. (2016), *Nat Commun* 7:10331, Methods: bins are built from
#' *consecutive degree values*, a bin extended until it holds at least
#' `min_per_bin` (100) nodes. This keeps bins fine-grained at low degree
#' and only coarsens in the sparse high-degree tail -- unlike fixed
#' quantile bins, where a single decile on a scale-free graph spans an
#' enormous degree range and a hub can be swapped for a low-degree node,
#' destroying degree preservation exactly where the hub bias it corrects
#' for is strongest.
#' @return A list of integer vectors (node indices), one per bin.
#' @keywords internal
.network_degree_bins <- function(degree_all, min_per_bin = 100) {
  ord <- order(degree_all)
  degs_sorted <- degree_all[ord]
  uniq_degs <- unique(degs_sorted)

  ## walk up the unique degrees, closing a bin once it has >= min_per_bin
  bin_of_degree <- integer(length(uniq_degs))
  b <- 1L
  count <- 0L
  for (i in seq_along(uniq_degs)) {
    bin_of_degree[i] <- b
    count <- count + sum(degs_sorted == uniq_degs[i])
    if (count >= min_per_bin) { b <- b + 1L; count <- 0L }
  }
  ## a trailing under-full bin is merged into the previous one
  if (count > 0 && b > 1L) bin_of_degree[bin_of_degree == b] <- b - 1L

  node_bin <- bin_of_degree[match(degree_all, uniq_degs)]
  split(seq_along(degree_all), node_bin)
}

#' Replace `string_ids` with an equal-size set of *distinct* nodes drawn
#' from the same degree bins.
#'
#' @details
#' Nodes are grouped by their degree bin and, per bin, `k` distinct
#' replacements are drawn from that bin's pool -- so the returned set has
#' exactly `length(string_ids)` distinct nodes, matching Guney et al.'s
#' "a set of |S| proteins with matching degrees". Drawing each node
#' independently with replacement (and then `unique()`-ing downstream)
#' would shrink the random set, inflating its distance variance and
#' shrinking `|z_score|` in a way that varies with set size.
#' @return Character vector of STRING IDs, length `length(string_ids)`.
#' @keywords internal
.network_resample_degree_matched <- function(string_ids, degree_all, bins) {
  node_names <- names(degree_all)
  bin_of_node <- integer(length(degree_all))
  for (b in seq_along(bins)) bin_of_node[bins[[b]]] <- b

  idx <- match(string_ids, node_names)
  by_bin <- split(seq_along(idx), bin_of_node[idx])

  out <- character(length(string_ids))
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
    out[slots] <- node_names[picks]
  }
  out
}

#' Guney et al. (2016) "closest" distance between two node sets on `g`
#'
#' @description
#' `unique()` on both `source_ids` and `target_ids` before calling
#' `igraph::distances()`: igraph's underlying C routine rejects a `to=`
#' argument with duplicate vertices. `.network_resample_degree_matched()`
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

#' @keywords internal
.empty_network_proximity_row <- function() {
  data.frame(
    condition = character(0), compound_id = character(0), disease_id = character(0),
    n_targets_mapped = integer(0), n_disease_genes_mapped = integer(0),
    d_observed = double(0), d_random_mean = double(0), d_random_sd = double(0),
    z_score = double(0), p_empirical = double(0), n_random = integer(0),
    seed_used = integer(0),
    stringsAsFactors = FALSE
  )
}
