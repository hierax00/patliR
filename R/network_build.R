#' @include AllGenerics.R internal.R
NULL

## An igraph object does not fit in a CSV, so the real output is a plain
## long-format edge list (`network_edges`); one igraph per condition is
## cached under cacheDir(proj) for speed. `.network_graph()` (below) reads
## that cache or rebuilds it from `network_edges`, so deleting the cache is
## always safe. See DESIGN.md.

#' Build a compound-target network for each condition (network pharmacology
#' core)
#'
#' @description
#' Builds one bipartite compound-target graph per condition column of
#' [binarizedMatrix()] -- the same per-condition fork established by
#' [prep_binarize()] -- using the predicted/imported targets in
#' [patliRResults()]. Only compounds actually present (`1`) in a condition,
#' and their targets, go into that condition's graph.
#'
#' This is the foundation the rest of the `network_*` family builds on.
#'
#' @section Why `target_source = "imported"` is the default:
#' `targets_consensus()` and `targets_bipartite()` are not implemented yet
#' (see `ROADMAP.md`), so `"imported"` (from
#' [targets_import()]/[targets_import_batch()]) is the only source that
#' exists; the other two raise a clear "not implemented yet" error if
#' requested.
#'
#' @section Edges are enriched with disease association when available:
#' If [targets_disease_filter()] has already been run, each edge also gets a
#' `disease_association_score` attribute/column (`NA` if that target has no
#' recorded association, or if `targets_disease_filter()` was never run at
#' all) -- this is purely descriptive at this stage; `network_proximity()`
#' (not implemented yet) is what will actually use it topologically.
#'
#' @inheritParams compounds
#' @param condition Character vector of condition names (must match columns
#'   of [binarizedMatrix()]), or `NULL` (default) to (re)build every
#'   condition. Rebuilding a subset leaves previously-built conditions
#'   untouched in `patliRResults(proj, "network_edges")`.
#' @param target_source `"imported"` (default and only implemented source
#'   right now), `"consensus"`, or `"bipartite"`. See the section above.
#' @param min_score `NULL` (default, keep every imported target) or a single
#'   number in `[0, 1]`: drop target edges whose `probability` (from
#'   [targets_import()]) is below this threshold before building the graph.
#'
#' @return The updated `proj`, with a `network_edges` entry in
#'   [patliRResults()] (columns `condition`, `compound_id`, `uniprot_id`,
#'   `weight` (the import probability), `disease_association_score`), also
#'   written to `results/network_edges.csv`. One `igraph` object per
#'   (re)built condition is cached under [cacheDir()] for the rest of the
#'   `network_*` family to reuse -- see `.network_graph()`, internal.
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
#' patliRResults(proj, "network_edges")
#' }
#'
#' @export
network_build <- function(proj, condition = NULL,
                           target_source = c("imported", "consensus", "bipartite"),
                           min_score = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  target_source <- match.arg(target_source)

  if (target_source != "imported") {
    cli::cli_abort(c(
      "{.val {target_source}} is not implemented yet in this version of patliR.",
      "i" = "Only {.val imported} (from {.fn targets_import}/{.fn targets_import_batch}) is available until {.code targets_{target_source}()} ships -- see {.file ROADMAP.md}."
    ))
  }
  if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1 || min_score < 0 || min_score > 1)) {
    cli::cli_abort("{.arg min_score} must be a single number in [0, 1], or NULL.")
  }

  bin <- binarizedMatrix(proj)
  if (nrow(bin) == 0) {
    cli::cli_abort(c(
      "No binarized abundance matrix in {.arg proj}.",
      "i" = "Run {.fn prep_binarize} first -- {.fn network_build} builds one network per condition column of {.fn binarizedMatrix}."
    ))
  }
  all_conditions <- setdiff(names(bin), "compound_id")

  if (!is.null(condition)) {
    unknown <- setdiff(condition, all_conditions)
    if (length(unknown) > 0) {
      cli::cli_abort("Unknown condition(s) {.val {unknown}}; available: {.val {all_conditions}}.")
    }
    conditions <- condition
  } else {
    conditions <- all_conditions
  }

  imported <- patliRResults(proj, "targets_imported")
  if (is.null(imported) || nrow(imported) == 0) {
    cli::cli_abort(c(
      "No {.val targets_imported} entry in {.arg proj}.",
      "i" = "Run {.fn targets_import} or {.fn targets_import_batch} first."
    ))
  }
  if (!is.null(min_score)) {
    below <- is.na(imported$probability) | imported$probability < min_score
    imported <- imported[!below, , drop = FALSE]
  }

  disease <- patliRResults(proj, "targets_disease")

  edge_list <- vector("list", length(conditions))
  names(edge_list) <- conditions
  n_edges <- integer(length(conditions))
  names(n_edges) <- conditions

  for (cond in conditions) {
    present_ids <- unique(bin$compound_id[!is.na(bin[[cond]]) & bin[[cond]] == 1])
    edges <- imported[imported$compound_id %in% present_ids, c("compound_id", "uniprot_id", "probability"), drop = FALSE]
    names(edges)[names(edges) == "probability"] <- "weight"
    rownames(edges) <- NULL

    ## targets_imported has no (compound_id, uniprot_id) uniqueness
    ## guarantee; a duplicate would become a multi-edge and inflate a
    ## target's degree beyond the number of distinct compounds hitting it.
    ## Keep the highest-weight row per pair, log the rest.
    dup_key <- paste(edges$compound_id, edges$uniprot_id)
    if (nrow(edges) > 0 && anyDuplicated(dup_key) > 0) {
      ord <- order(dup_key, -edges$weight)
      edges <- edges[ord, , drop = FALSE]
      keep <- !duplicated(dup_key[ord])
      dropped_n <- sum(!keep)
      edges <- edges[keep, , drop = FALSE]
      rownames(edges) <- NULL
      proj <- .log_append(
        proj, step = "network_build", id = NA_character_,
        message = paste0("condition '", cond, "': ", dropped_n, " duplicate (compound_id, uniprot_id) edge(s) collapsed, highest weight kept")
      )
    }

    edges$disease_association_score <- rep(NA_real_, nrow(edges))
    if (!is.null(disease) && nrow(disease) > 0 && nrow(edges) > 0) {
      key <- paste(edges$compound_id, edges$uniprot_id)
      dkey <- paste(disease$compound_id, disease$target_id)
      hit <- match(key, dkey)
      edges$disease_association_score <- disease$association_score[hit]
    }

    g <- .network_build_igraph(present_ids, edges)
    saveRDS(g, .network_cache_path(proj, cond))

    edges$condition <- rep(cond, nrow(edges))
    edge_list[[cond]] <- edges[, c("condition", "compound_id", "uniprot_id", "weight", "disease_association_score")]
    n_edges[[cond]] <- nrow(edges)
  }

  new_edges <- do.call(rbind, edge_list)
  rownames(new_edges) <- NULL

  existing <- patliRResults(proj, "network_edges")
  if (!is.null(existing) && nrow(existing) > 0) {
    existing <- existing[!existing$condition %in% conditions, , drop = FALSE]
    new_edges <- rbind(existing, new_edges)
  }

  patliRResults(proj, "network_edges") <- new_edges
  .write_results_csv(proj, "network_edges", new_edges)

  proj <- .log_append(
    proj, step = "network_build", id = NA_character_,
    message = paste0("condition '", conditions, "': ", n_edges[conditions], " compound-target edges")
  )
  .write_log_csv(proj)
  proj
}

#' Build one condition's igraph object from a compound id set and edge table
#' @return An undirected `igraph` object; compound nodes have `type = FALSE`,
#'   target nodes have `type = TRUE` (the standard igraph bipartite
#'   convention). Isolated compounds (no surviving target edges) are still
#'   included as nodes.
#' @keywords internal
.network_build_igraph <- function(compound_ids, edges) {
  target_ids <- if (nrow(edges) > 0) setdiff(unique(edges$uniprot_id), compound_ids) else character(0)
  vertices <- data.frame(
    name = c(compound_ids, target_ids),
    type = c(rep(FALSE, length(compound_ids)), rep(TRUE, length(target_ids))),
    stringsAsFactors = FALSE
  )
  igraph::graph_from_data_frame(
    d = edges[, c("compound_id", "uniprot_id", "weight", "disease_association_score"), drop = FALSE],
    directed = FALSE, vertices = vertices
  )
}

#' Node "mode" labels for a bipartite compound-target graph
#'
#' @description
#' The shared accessor for the compound/target split every mode-aware
#' `network_*` function needs, so the `as.logical(V(g)$type)` idiom is not
#' re-derived (slightly differently each time) in
#' `.network_centrality_one()`, `network_hub_penalty()`,
#' `network_module_robustness()`, etc. `V(g)$type == FALSE` is a compound,
#' `TRUE` is a target -- the convention `.network_build_igraph()` sets.
#' @return `character` vector, one entry per vertex in `V(g)` order, each
#'   `"compound"` or `"target"` (`NA` for a vertex whose `type` is `NA`).
#' @keywords internal
.network_node_types <- function(g) {
  ty <- as.logical(igraph::V(g)$type)
  ifelse(is.na(ty), NA_character_, ifelse(ty, "target", "compound"))
}

#' Is `g` a genuine two-mode (bipartite) graph, per `network_build()`'s own
#' construction contract?
#'
#' @description
#' Checks that every vertex carries a non-`NA` `type` and that **no edge
#' connects two vertices of the same mode**. This is deliberately *not*
#' `igraph::bipartite_mapping()`: that function (a) throws on a named graph
#' with a non-bipartite structure instead of returning `FALSE`, (b) reports
#' `TRUE` for any 2-colourable graph -- including the tri/quadripartite
#' `.network_layered_graph()` (compounds and pathways share a colour) and
#' any tree -- and (c) ignores an existing `V(g)$type`, computing its own
#' 2-colouring, while `as_biadjacency_matrix()` trusts `V(g)$type` blindly.
#' The check here matches what the bipartite modularity path actually
#' requires: that `as_biadjacency_matrix(g)` is a faithful representation.
#' @return `logical(1)`.
#' @keywords internal
.network_is_bipartite <- function(g) {
  if (!"type" %in% igraph::vertex_attr_names(g)) return(FALSE)
  ty <- as.logical(igraph::V(g)$type)
  if (anyNA(ty)) return(FALSE)
  el <- igraph::as_edgelist(g, names = FALSE)
  nrow(el) == 0L || all(ty[el[, 1]] != ty[el[, 2]])
}

#' Path to the cached igraph object for one project condition
#' @keywords internal
.network_cache_path <- function(proj, condition) {
  file.path(cacheDir(proj), paste0("network_", .cache_key_slug(condition), ".rds"))
}

#' Get the igraph object for one condition -- the shared entry point every
#' other `network_*` function uses
#'
#' @description
#' Reads the cached `.rds` for `condition` if present; otherwise transparently
#' rebuilds it from `patliRResults(proj, "network_edges")` (falling back to
#' `results/network_edges.csv` on disk if `proj` was reloaded fresh via
#' [patliR_load()] and the in-memory results bag does not have it) and
#' repopulates the cache. Deleting the cache is always safe, same principle
#' as `refdb_rebuild_cache()`.
#' @return An `igraph` object.
#' @keywords internal
.network_graph <- function(proj, condition) {
  cache_path <- .network_cache_path(proj, condition)
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all)) {
    csv_path <- file.path(projectDir(proj), "results", "network_edges.csv")
    if (file.exists(csv_path)) {
      edges_all <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
    }
  }
  if (is.null(edges_all) || !condition %in% edges_all$condition) {
    cli::cli_abort(c(
      "No network built for condition {.val {condition}}.",
      "i" = "Run {.fn network_build} first."
    ))
  }

  edges <- edges_all[edges_all$condition == condition, , drop = FALSE]
  bin <- binarizedMatrix(proj)
  if (!condition %in% names(bin)) {
    cli::cli_abort("Condition {.val {condition}} not found in {.fn binarizedMatrix}; cannot rebuild its graph.")
  }
  present_ids <- bin$compound_id[!is.na(bin[[condition]]) & bin[[condition]] == 1]

  g <- .network_build_igraph(present_ids, edges)
  saveRDS(g, cache_path)
  g
}

#' Resolve/validate a `condition` argument against what `network_build()`
#' has actually built -- the shared validation every later `network_*`
#' function uses
#'
#' @return Character vector of condition names to operate on.
#' @keywords internal
.network_resolve_conditions <- function(proj, condition) {
  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_edges} entry in {.arg proj}.",
      "i" = "Run {.fn network_build} first."
    ))
  }
  built <- unique(edges_all$condition)
  if (is.null(condition)) {
    return(built)
  }
  unknown <- setdiff(condition, built)
  if (length(unknown) > 0) {
    cli::cli_abort("Condition(s) {.val {unknown}} were not built by {.fn network_build}; available: {.val {built}}.")
  }
  condition
}

#' Replace rows of an existing `patliRResults()` table that share a key with
#' `new_rows`, keeping every other existing row untouched
#'
#' @description
#' The shared "rebuild a subset without nuking the rest" pattern already
#' used ad hoc in [network_build()] -- factored out here so every later
#' `network_*` function (most of which are per-condition, some per
#' `(condition, method)` or `(condition, db)`) gets it for free and
#' consistently, instead of re-deriving slightly different upsert logic
#' per function.
#'
#' @param key_cols Character vector of column names that jointly identify a
#'   "slot" to replace (e.g. `"condition"`, or `c("condition", "method")`).
#' @param touched_keys `NULL` (default) or a data frame carrying the
#'   `key_cols` values the caller just recomputed. When supplied, exactly
#'   those slots are dropped from the existing table and replaced by
#'   whatever `new_rows` holds -- **including zero rows**, so a rerun that
#'   now produces nothing for a condition correctly removes that
#'   condition's stale rows instead of silently keeping them. When `NULL`,
#'   the touched slots are derived from `new_rows` itself (the historical
#'   behaviour, kept for callers not yet migrated): a zero-row `new_rows`
#'   then touches nothing and the previous table is returned unchanged.
#' @return `new_rows` with any previously-existing, non-overlapping rows of
#'   `patliRResults(proj, name)` prepended.
#' @keywords internal
.network_upsert <- function(proj, name, new_rows, key_cols, touched_keys = NULL) {
  existing <- patliRResults(proj, name)
  if (is.null(existing) || nrow(existing) == 0) {
    rownames(new_rows) <- NULL
    return(new_rows)
  }

  if (is.null(touched_keys)) {
    touched <- if (nrow(new_rows) > 0) {
      unique(do.call(paste, c(new_rows[key_cols], sep = "\r")))
    } else {
      character(0)
    }
  } else {
    touched_keys <- as.data.frame(touched_keys, stringsAsFactors = FALSE)
    missing_key <- setdiff(key_cols, names(touched_keys))
    if (length(missing_key) > 0) {
      cli::cli_abort(c(
        "{.fn .network_upsert}: {.arg touched_keys} is missing key column{?s} {.val {missing_key}}.",
        "i" = "It must carry every column named in {.arg key_cols}: {.val {key_cols}}."
      ))
    }
    touched <- unique(do.call(paste, c(touched_keys[key_cols], sep = "\r")))
  }

  key_existing <- do.call(paste, c(existing[key_cols], sep = "\r"))
  kept <- existing[!key_existing %in% touched, , drop = FALSE]

  if (nrow(new_rows) > 0) {
    only_new <- setdiff(names(new_rows), names(kept))
    only_old <- setdiff(names(kept), names(new_rows))
    shared   <- intersect(names(new_rows), names(kept))
    ## An all-NA column (common when a CSV round-trip reads an empty column
    ## back as logical) is compatible with anything; the numeric family
    ## (integer/double) rbinds without loss. Everything else must match.
    .col_compat <- function(a, b) {
      if (all(is.na(a)) || all(is.na(b))) return(TRUE)
      num <- c("integer", "numeric", "double")
      ca <- class(a)[1]; cb <- class(b)[1]
      if (ca %in% num && cb %in% num) return(TRUE)
      identical(ca, cb)
    }
    type_conflict <- shared[!vapply(shared, function(cn)
      .col_compat(kept[[cn]], new_rows[[cn]]), logical(1))]

    ## A pure column *addition* (new schema has columns the old table lacks,
    ## and nothing the other way, no type change on a shared column) is a
    ## forward migration: back-fill typed NA on the pre-existing rows and
    ## proceed, with one summary line. A column *removal* or a *type
    ## conflict* on a shared column still aborts -- those are not something
    ## this function can silently reconcile.
    if (length(only_old) > 0 || length(type_conflict) > 0) {
      cli::cli_abort(c(
        "{.fn .network_upsert}: the new {.val {name}} rows and the existing table are not reconcilable.",
        if (length(only_old) > 0) c("i" = "Only in the existing table (column removed): {.val {only_old}}") else NULL,
        if (length(type_conflict) > 0) c("i" = "Type changed on shared column{?s}: {.val {type_conflict}}") else NULL,
        "i" = "Delete {.file results/{name}.csv} and re-run the step to rebuild it under the current schema."
      ))
    }
    if (length(only_new) > 0) {
      for (cn in only_new) kept[[cn]] <- new_rows[[cn]][rep(NA_integer_, nrow(kept))]
      if (nrow(kept) > 0) {
        cli::cli_inform(c(
          "i" = "{.fn .network_upsert}: {.val {name}} gained column{?s} {.val {only_new}}; back-filled as NA on {nrow(kept)} pre-existing row{?s}."
        ))
      }
    }
    out <- rbind(kept[, names(new_rows), drop = FALSE], new_rows)
  } else {
    out <- kept
  }
  rownames(out) <- NULL
  out
}
