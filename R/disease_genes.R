#' @include AllGenerics.R internal.R network_build.R
NULL

## An independent disease gene set, fetched in the disease -> target
## direction (or imported from a curated list), stored in its own
## `disease_genes` results slot.
##
## Why this is separate from `targets_disease_filter()`:
## `targets_disease_filter()` only ever annotates UniProt IDs that are
## already in `targets_imported` -- so the set it produces is a subset of
## the compounds' own predicted targets. Using that as `network_proximity()`'s
## "disease module" makes the proximity z-score circular: S and T overlap by
## construction and the statistic measures set membership, not topology.
## `disease_genes_fetch()` / `disease_genes_import()` build a disease gene
## set that does not depend on the compounds at all.

#' Fetch a disease's associated gene set from Open Targets (disease -> target)
#'
#' @description
#' Queries the Open Targets Platform in the **disease -> target** direction
#' (`disease(efoId:){ associatedTargets }`) and stores the resulting gene
#' set in a `disease_genes` results slot, independent of any compound or of
#' `targets_imported`. This is what [network_proximity()] should build its
#' disease module `T` from: unlike [targets_disease_filter()] (which only
#' annotates UniProt IDs already predicted as targets), nothing here is a
#' function of the compounds, so the proximity z-score measures topology
#' rather than set membership.
#'
#' `disease` is resolved to a current Open Targets EFO/MONDO ID exactly the
#' way [targets_disease_filter()] does it -- free text goes through Open
#' Targets' `mapIds` full-text search (always resolves to the *current*
#' term), an ID-shaped string is looked up and validated directly (see
#' `.open_targets_resolve_disease()`, internal). Free text is the safer
#' default; ontology IDs get retired.
#'
#' @section Open Targets score vs. a curated disease-gene list:
#' Menche et al. (2015) and Guney et al. (2016) build the disease module
#' from curated lists (OMIM + GWAS Catalog). Open Targets' aggregated
#' `associationScore` (genetic association, known drug, literature, ...) is
#' a reasonable substitute but it is a continuous score, not a curated
#' membership call, so any `min_score` cut is a scientific choice and is
#' always logged. `disease_genes_import()` is the route for a genuinely
#' curated list.
#'
#' @inheritParams targets_disease_filter
#' @param disease A disease/phenotype name (e.g. `"Parkinson disease"`) or a
#'   currently-valid ontology identifier (EFO/MONDO/Orphanet/..., e.g.
#'   `"MONDO_0005180"`). Free text is the safer default. Required.
#' @param source Only `"open_targets"` is implemented. Kept explicit (rather
#'   than hardcoded) so a future alternative source raises an informative
#'   error instead of being silently ignored -- same pattern as
#'   [targets_disease_filter()] and [refdb_build()].
#' @param min_score `NULL` (default): keep every associated target that maps
#'   to at least one Swiss-Prot UniProt accession. A numeric threshold in
#'   `[0, 1]`: additionally drop rows whose `association_score` is below it
#'   (always logged, by UniProt ID and score, so nothing disappears
#'   silently). The default is `NULL` for consistency with
#'   [targets_disease_filter()].
#' @param fetch_mode `"warn_and_cache"` (default) or `"abort"`. See
#'   `.fetch_external()` (internal).
#'
#' @return The updated `proj`, with a `disease_genes` entry in
#'   [patliRResults()] (columns `disease_id`, `disease_name`, `uniprot_id`,
#'   `gene_symbol`, `association_score`, `source`, `fetched_at`), also
#'   written to `results/disease_genes.csv`. Re-fetching one disease
#'   replaces only that disease's rows (upsert keyed on `disease_id`).
#'   Associated targets with no `uniprot_swissprot` protein ID are logged
#'   (`"disease_genes_no_swissprot"`) and excluded; a target with more than
#'   one Swiss-Prot accession contributes one row per accession.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' proj <- disease_genes_fetch(proj, disease = "type 2 diabetes mellitus") # needs internet
#' patliRResults(proj, "disease_genes")
#' }
#'
#' @export
disease_genes_fetch <- function(proj, disease, source = c("open_targets"),
                                 min_score = NULL,
                                 fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"))
  source <- match.arg(source)
  fetch_mode <- match.arg(fetch_mode)

  if (missing(disease) || is.null(disease) || !is.character(disease) ||
      length(disease) != 1 || !nzchar(disease)) {
    cli::cli_abort(c(
      "{.arg disease} is required and must be a single non-empty string.",
      "i" = "This function fetches one disease at a time -- call it again for another."
    ))
  }
  if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1 ||
                              min_score < 0 || min_score > 1)) {
    cli::cli_abort("{.arg min_score} must be a single number in [0, 1], or NULL.")
  }

  efo_id <- .fetch_external(
    fetch_fun = function() .open_targets_resolve_disease(disease),
    cache_dir = cacheDir(proj),
    cache_key = paste0("opentargets_disease_", .cache_key_slug(disease)),
    mode = fetch_mode
  )
  if (is.null(efo_id)) {
    cli::cli_warn("Could not resolve {.val {disease}} to an Open Targets disease ID; nothing to do.")
    return(proj)
  }

  fetched <- .fetch_external(
    fetch_fun = function() .open_targets_disease_targets(efo_id),
    cache_dir = cacheDir(proj),
    cache_key = paste0("opentargets_disease_targets_", .cache_key_slug(efo_id)),
    mode = fetch_mode
  )
  if (is.null(fetched)) {
    cli::cli_warn("Could not fetch associated targets for {.val {disease}} ({.val {efo_id}}); nothing to do.")
    return(proj)
  }

  disease_name <- fetched$disease_name %||% NA_character_
  rows <- fetched$rows

  if (nrow(rows) > 0) {
    no_uniprot <- is.na(rows$uniprot_id)
    if (any(no_uniprot)) {
      proj <- .log_append(
        proj, step = "disease_genes_fetch", id = rows$ensembl_id[no_uniprot],
        message = paste0(
          "disease_genes_no_swissprot: Open Targets target '", rows$ensembl_id[no_uniprot],
          "' (", rows$gene_symbol[no_uniprot], ") has no uniprot_swissprot proteinId; excluded"
        )
      )
    }
    rows <- rows[!no_uniprot, , drop = FALSE]
  }

  fetched_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  new_rows <- data.frame(
    disease_id = rep(efo_id, nrow(rows)),
    disease_name = rep(disease_name, nrow(rows)),
    uniprot_id = as.character(rows$uniprot_id),
    gene_symbol = as.character(rows$gene_symbol),
    association_score = as.numeric(rows$association_score),
    source = rep("open_targets", nrow(rows)),
    fetched_at = rep(fetched_at, nrow(rows)),
    stringsAsFactors = FALSE
  )
  new_rows <- unique(new_rows)

  if (!is.null(min_score) && nrow(new_rows) > 0) {
    below <- !is.na(new_rows$association_score) & new_rows$association_score < min_score
    if (any(below)) {
      proj <- .log_append(
        proj, step = "disease_genes_fetch", id = new_rows$uniprot_id[below],
        message = paste0(
          "disease_genes_below_min_score: association_score ",
          format(round(new_rows$association_score[below], 4), nsmall = 4),
          " < min_score ", min_score, "; dropped"
        )
      )
    }
    new_rows <- new_rows[!below, , drop = FALSE]
  }

  proj <- .log_append(
    proj, step = "disease_genes_fetch", id = efo_id,
    message = paste0(
      "disease_genes_fetch: '", disease, "' resolved to '", efo_id, "' (", disease_name,
      "); kept ", nrow(new_rows), " Swiss-Prot gene(s); ",
      if (is.null(min_score)) "min_score = NULL (all kept)" else paste0("min_score = ", min_score)
    )
  )

  if (nrow(new_rows) == 0) new_rows <- .empty_disease_genes_row()

  result <- .network_upsert(
    proj, "disease_genes", new_rows, "disease_id",
    touched_keys = data.frame(disease_id = efo_id, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "disease_genes") <- result
  .write_results_csv(proj, "disease_genes", result)
  .write_log_csv(proj)
  proj
}

#' Import a curated disease gene list into the `disease_genes` slot
#'
#' @description
#' The manual counterpart to [disease_genes_fetch()], for curated
#' disease-gene lists exported from OMIM, the GWAS Catalog, DisGeNET, or a
#' hand-assembled table. Writes into the same `disease_genes` slot / schema
#' as [disease_genes_fetch()], with `source` recording provenance and
#' `association_score` allowed to be `NA` (a curated list carries no Open
#' Targets score). Follows the shape of the other `*_import()` functions in
#' the package ([adme_import()], [tox_import()], [targets_import()]).
#'
#' @section UniProt accessions, not gene symbols:
#' `table` must carry a column of UniProt accessions. A `gene_symbol` column
#' is accepted and stored if present, but `disease_genes_import()` does not
#' map symbols to accessions -- convert them first (UniProt's ID-mapping
#' tool, or `org.Hs.eg.db`). This keeps the import offline and deterministic
#' and avoids silently picking one accession per multi-mapping symbol.
#'
#' @inheritParams compounds
#' @param table A `data.frame`, or a path to a CSV file, with at least a
#'   UniProt-ID column. A `gene_symbol` column and an
#'   `association_score`/`score` column are used if present.
#' @param disease_id The disease identifier these genes belong to (e.g. an
#'   EFO/MONDO ID, or any stable label). Becomes the `disease_id` that
#'   [network_proximity()] is called with. Required.
#' @param disease_name Optional human-readable disease name, stored in the
#'   `disease_name` column (`NA` if not given).
#' @param source Free-text provenance string for the `source` column (e.g.
#'   `"OMIM"`, `"GWAS Catalog"`, `"DisGeNET_curated"`). Default `"manual"`.
#' @param uniprot_col,gene_symbol_col,score_col Optional explicit column
#'   names in `table`. When `NULL` (default) they are auto-detected
#'   (case/whitespace-insensitively) from the usual spellings
#'   (`uniprot_id`/`uniprot`/`accession`/...; `gene_symbol`/`symbol`/`gene`;
#'   `association_score`/`score`).
#'
#' @return The updated `proj`, with the curated rows upserted into the
#'   `disease_genes` slot (keyed on `disease_id`, so this replaces any
#'   previous rows for the same disease), also written to
#'   `results/disease_genes.csv`. Rows with a missing/empty UniProt ID are
#'   logged (`"disease_genes_import_no_uniprot"`) and excluded.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' curated <- data.frame(
#'   uniprot_id = c("P37840", "P05067", "P10636"),
#'   gene_symbol = c("SNCA", "APP", "MAPT")
#' )
#' proj <- disease_genes_import(
#'   proj, curated,
#'   disease_id = "MONDO_0005180", disease_name = "Parkinson disease",
#'   source = "curated_demo"
#' )
#' patliRResults(proj, "disease_genes")
#'
#' @export
disease_genes_import <- function(proj, table, disease_id, disease_name = NULL,
                                  source = "manual",
                                  uniprot_col = NULL, gene_symbol_col = NULL,
                                  score_col = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  if (missing(disease_id) || is.null(disease_id) || !is.character(disease_id) ||
      length(disease_id) != 1 || !nzchar(disease_id)) {
    cli::cli_abort("{.arg disease_id} is required and must be a single non-empty string.")
  }
  if (!is.character(source) || length(source) != 1 || !nzchar(source)) {
    cli::cli_abort("{.arg source} must be a single non-empty string.")
  }
  if (!is.null(disease_name) && (!is.character(disease_name) || length(disease_name) != 1)) {
    cli::cli_abort("{.arg disease_name} must be a single string, or NULL.")
  }

  raw <- if (is.data.frame(table)) {
    table
  } else if (is.character(table) && length(table) == 1 && file.exists(table)) {
    .read_csv_safe(table)
  } else {
    cli::cli_abort("{.arg table} must be a data frame or a path to an existing CSV file.")
  }
  if (nrow(raw) == 0) {
    cli::cli_abort("{.arg table} has no rows.")
  }

  find_col <- function(explicit, candidates) {
    if (!is.null(explicit)) return(.match_column_flexible(names(raw), explicit))
    for (cand in candidates) {
      hit <- .match_column_flexible(names(raw), cand)
      if (!is.na(hit)) return(hit)
    }
    NA_character_
  }
  uni_col <- find_col(uniprot_col, c("uniprot_id", "uniprot", "accession", "uniprot_accession", "swissprot"))
  sym_col <- find_col(gene_symbol_col, c("gene_symbol", "symbol", "gene", "approved_symbol"))
  sc_col  <- find_col(score_col, c("association_score", "score"))

  if (is.na(uni_col)) {
    cli::cli_abort(c(
      "No UniProt-ID column found in {.arg table}.",
      "i" = "Provide a column of UniProt accessions (e.g. named {.val uniprot_id}), or pass {.arg uniprot_col}.",
      if (!is.na(sym_col)) c(
        "x" = "A gene-symbol column ({.val {sym_col}}) is present, but {.fn disease_genes_import} does not map symbols to accessions."
      ) else NULL,
      if (!is.na(sym_col)) c(
        "i" = "Convert the symbols to UniProt accessions first (UniProt ID mapping, or {.pkg org.Hs.eg.db})."
      ) else NULL
    ))
  }

  uniprot_id <- trimws(as.character(raw[[uni_col]]))
  gene_symbol <- if (!is.na(sym_col)) as.character(raw[[sym_col]]) else rep(NA_character_, nrow(raw))
  association_score <- if (!is.na(sc_col)) suppressWarnings(as.numeric(raw[[sc_col]])) else rep(NA_real_, nrow(raw))

  missing_uni <- is.na(uniprot_id) | !nzchar(uniprot_id) | tolower(uniprot_id) == "na"
  if (any(missing_uni)) {
    proj <- .log_append(
      proj, step = "disease_genes_import", id = NA_character_,
      message = paste0("disease_genes_import_no_uniprot: row ", which(missing_uni),
                        " of the curated table has no UniProt ID; excluded")
    )
  }

  keep <- !missing_uni
  new_rows <- data.frame(
    disease_id = rep(disease_id, sum(keep)),
    disease_name = rep(disease_name %||% NA_character_, sum(keep)),
    uniprot_id = uniprot_id[keep],
    gene_symbol = gene_symbol[keep],
    association_score = association_score[keep],
    source = rep(source, sum(keep)),
    fetched_at = rep(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sum(keep)),
    stringsAsFactors = FALSE
  )
  new_rows <- unique(new_rows)

  proj <- .log_append(
    proj, step = "disease_genes_import", id = disease_id,
    message = paste0(
      "disease_genes_import: ", nrow(new_rows), " curated gene(s) for disease '",
      disease_id, "' (source '", source, "')",
      if (all(is.na(new_rows$association_score))) "; no association_score (curated list)" else ""
    )
  )

  if (nrow(new_rows) == 0) new_rows <- .empty_disease_genes_row()

  result <- .network_upsert(
    proj, "disease_genes", new_rows, "disease_id",
    touched_keys = data.frame(disease_id = disease_id, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "disease_genes") <- result
  .write_results_csv(proj, "disease_genes", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.empty_disease_genes_row <- function() {
  data.frame(
    disease_id = character(0), disease_name = character(0),
    uniprot_id = character(0), gene_symbol = character(0),
    association_score = double(0), source = character(0),
    fetched_at = character(0), stringsAsFactors = FALSE
  )
}

#' Fetch every target associated with a disease from Open Targets, in the
#' disease -> target direction, paginating on `count`.
#'
#' @return `list(disease_id, disease_name, rows)` where `rows` is a
#'   `data.frame` with columns `ensembl_id`, `gene_symbol`,
#'   `association_score`, `uniprot_id` -- one row per (target, Swiss-Prot
#'   accession), or one row with `uniprot_id = NA` for a target that has no
#'   `uniprot_swissprot` protein ID.
#' @keywords internal
.open_targets_disease_targets <- function(efo_id, page_size = 500L) {
  query_string <- "
    query DiseaseTargets($efoId: String!, $index: Int!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: $index, size: $size}) {
          count
          rows {
            score
            target { id approvedSymbol proteinIds { id source } }
          }
        }
      }
    }
  "

  index <- 0L
  disease_name <- NULL
  collected <- list()
  n_seen <- 0L
  repeat {
    data <- .open_targets_graphql(
      query_string = query_string,
      variables = list(efoId = efo_id, index = index, size = page_size)
    )
    dz <- data$disease
    if (is.null(dz)) {
      stop(paste0("Open Targets: disease(efoId = \"", efo_id, "\") returned no record."))
    }
    if (is.null(disease_name)) disease_name <- dz$name %||% NA_character_
    at <- dz$associatedTargets
    count <- at$count %||% 0L
    page_rows <- at$rows %||% list()
    if (length(page_rows) == 0) break
    collected <- c(collected, page_rows)
    n_seen <- n_seen + length(page_rows)
    index <- index + 1L
    if (n_seen >= count) break
    if (index > 100000L) break # unreachable safety stop
  }

  if (length(collected) == 0) {
    return(list(disease_id = efo_id, disease_name = disease_name,
                rows = data.frame(ensembl_id = character(0), gene_symbol = character(0),
                                   association_score = double(0), uniprot_id = character(0),
                                   stringsAsFactors = FALSE)))
  }

  row_frames <- lapply(collected, function(r) {
    tgt <- r$target
    ensembl_id <- tgt$id %||% NA_character_
    gene_symbol <- tgt$approvedSymbol %||% NA_character_
    score <- suppressWarnings(as.numeric(r$score %||% NA_real_))
    pids <- tgt$proteinIds %||% list()
    swissprot <- character(0)
    if (length(pids) > 0) {
      swissprot <- vapply(pids, function(p) {
        if (identical(p$source, "uniprot_swissprot")) p$id %||% NA_character_ else NA_character_
      }, character(1))
      swissprot <- unique(swissprot[!is.na(swissprot) & nzchar(swissprot)])
    }
    if (length(swissprot) == 0) swissprot <- NA_character_
    data.frame(
      ensembl_id = ensembl_id, gene_symbol = gene_symbol,
      association_score = score, uniprot_id = swissprot,
      stringsAsFactors = FALSE
    )
  })

  list(disease_id = efo_id, disease_name = disease_name,
       rows = do.call(rbind, row_frames))
}
