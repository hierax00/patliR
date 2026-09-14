#' @include AllGenerics.R internal.R
NULL

## target <-> disease association via the Open Targets GraphQL API v4 (no
## API key). Open Targets rather than GeneCards for license/scraping
## reasons. Optional -- never blocks the pipeline if `disease` is not
## supplied.

#' Filter/annotate predicted targets by their association to a disease
#'
#' @description
#' Looks up, for every distinct UniProt ID already present in
#' `patliRResults(proj, "targets_imported")`, its overall association score
#' to a single disease/phenotype in the Open Targets Platform (target-disease
#' association, aggregated across all of Open Targets' evidence types --
#' genetic association, known drug, literature, etc.). This is an optional,
#' enrichment-only step -- nothing else in `patliR` requires it to run.
#'
#' **Not a disease gene set for [network_proximity()].** This function only
#' annotates UniProt IDs that are already predicted targets of the compounds
#' under study, so its output is by construction a subset of those
#' compounds' own targets. Using it as the "disease module" would make
#' [network_proximity()]'s z-score circular (it would measure set membership,
#' not topology) -- use [disease_genes_fetch()] / [disease_genes_import()]
#' instead, which build the disease gene set independently.
#'
#' `targets_disease_filter()` never predicts targets itself and never touches
#' [compounds()] -- it only annotates the UniProt IDs already imported via
#' [targets_import()]/[targets_import_batch()]. Run one of those first.
#'
#' @inheritParams compounds
#' @param disease A disease or phenotype name (e.g. `"hypertension"`) or an
#'   existing, *currently valid* ontology identifier (EFO/MONDO/Orphanet/etc.,
#'   e.g. `"MONDO_0005044"`, hypertension's current MONDO ID). Free text is
#'   the safer default -- ontology IDs do get retired (e.g. hypertension's
#'   old EFO ID, `EFO_0000537`, was obsoleted in EFO 3.88.0 and replaced by
#'   `MONDO_0005044`; other catalogs, like GWAS Catalog, still link the old
#'   one, but Open Targets does not), and passing an obsolete ID fails with
#'   no error from the API itself -- just an empty result -- which this
#'   function surfaces as a clear "could not resolve" error rather than
#'   silently returning nothing (see `.open_targets_resolve_disease()`,
#'   internal). An ID-shaped string (`PREFIX_alphanumeric`) is looked up
#'   directly and validated against Open Targets; anything else (free text)
#'   is resolved via Open Targets' `mapIds` full-text search (same engine
#'   behind their web search box), which always resolves to the *current*
#'   term -- ambiguous free-text terms take the top-ranked hit
#'   (`mapIds` does not match on the ID string itself,
#'   confirmed against the live API -- this is why the two paths are
#'   separate rather than always going through `mapIds`). Required: this
#'   function only makes sense for one disease at a time, call it again for
#'   another.
#' @param source Only `"open_targets"` is implemented. Kept as an explicit
#'   argument (rather than hardcoding the source) so a future alternative
#'   source raises an informative error instead of silently being ignored --
#'   see [refdb_build()] for the same pattern with `"coconut"`.
#' @param min_score `NULL` (default): keep every UniProt ID that was
#'   successfully queried, whether or not Open Targets reports an
#'   association (rows with no association get `association_score = NA` and
#'   are dropped, since there is nothing left to rank -- this is always
#'   logged). A numeric threshold in `[0, 1]`: additionally drop rows whose
#'   `association_score` is below it (also always logged, by UniProt ID and
#'   score, so nothing disappears silently).
#' @param fetch_mode `"warn_and_cache"` (default) or `"abort"`. See
#'   `.fetch_external()` (internal): this is enrichment
#'   data, so failures do not stop the pipeline by default -- a UniProt ID
#'   that cannot be resolved/queried is logged and dropped, not fatal.
#'
#' @return The updated `proj`, with a `targets_disease` entry in
#'   [patliRResults()] (columns `compound_id`, `target_id` (the UniProt ID),
#'   `disease_id` (the resolved EFO ID), `association_score`, `evidence`
#'   (semicolon-separated `datatype=score` breakdown, or `NA` when Open
#'   Targets has no evidence at all for the pair)), also written to
#'   `results/targets_disease.csv`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- targets_import(
#'   proj,
#'   system.file("extdata", "import_targets", "Targets5280443.csv", package = "patliR"),
#'   platform = "superpred"
#' )
#' proj <- targets_disease_filter(proj, disease = "hypertension") # needs internet
#' patliRResults(proj, "targets_disease")
#' }
#'
#' @export
targets_disease_filter <- function(proj, disease, source = c("open_targets"),
                                    min_score = NULL,
                                    fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"))
  source <- match.arg(source)
  fetch_mode <- match.arg(fetch_mode)

  if (missing(disease) || is.null(disease) || !is.character(disease) || length(disease) != 1 || !nzchar(disease)) {
    cli::cli_abort(c(
      "{.arg disease} is required and must be a single non-empty string.",
      "i" = "This function only makes sense for one disease at a time -- call it again for another."
    ))
  }
  if (!is.null(min_score) && (!is.numeric(min_score) || length(min_score) != 1 || min_score < 0 || min_score > 1)) {
    cli::cli_abort("{.arg min_score} must be a single number in [0, 1], or NULL.")
  }

  imported <- patliRResults(proj, "targets_imported")
  if (is.null(imported) || nrow(imported) == 0) {
    cli::cli_abort(c(
      "No {.val targets_imported} entry in {.arg proj}.",
      "i" = "Run {.fn targets_import} or {.fn targets_import_batch} first."
    ))
  }

  efo_id <- .fetch_external(
    fetch_fun = function() .open_targets_resolve_disease(disease),
    cache_dir = cacheDir(proj), cache_key = paste0("opentargets_disease_", .cache_key_slug(disease)),
    mode = fetch_mode
  )
  if (is.null(efo_id)) {
    cli::cli_warn("Could not resolve {.val {disease}} to an Open Targets disease ID; nothing to do.")
    return(proj)
  }

  uniprot_ids <- sort(unique(imported$uniprot_id))
  rows <- vector("list", length(uniprot_ids))
  names(rows) <- uniprot_ids

  for (u in uniprot_ids) {
    ensembl_id <- .fetch_external(
      fetch_fun = function() .open_targets_map_id(u, entity = "target"),
      cache_dir = cacheDir(proj), cache_key = paste0("opentargets_target_", .cache_key_slug(u)),
      mode = fetch_mode
    )
    if (is.null(ensembl_id)) {
      proj <- .log_append(
        proj, step = "targets_disease_filter", id = u,
        message = paste0("targets_disease_unresolved_target: could not map UniProt ID '", u, "' to an Open Targets target ID")
      )
      next
    }

    assoc <- .fetch_external(
      fetch_fun = function() .open_targets_target_disease_score(ensembl_id, efo_id),
      cache_dir = cacheDir(proj),
      cache_key = paste0("opentargets_assoc_", .cache_key_slug(ensembl_id), "_", .cache_key_slug(efo_id)),
      mode = fetch_mode
    )
    if (is.null(assoc)) {
      proj <- .log_append(
        proj, step = "targets_disease_filter", id = u,
        message = paste0("targets_disease_unresolved_association: could not query Open Targets for '", u, "' <-> '", disease, "'")
      )
      next
    }

    rows[[u]] <- data.frame(
      target_id = u, disease_id = assoc$disease_id %||% efo_id,
      association_score = assoc$score, evidence = assoc$evidence %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  scores <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(scores)) scores <- .empty_targets_disease_scores()

  no_assoc <- is.na(scores$association_score)
  if (any(no_assoc)) {
    proj <- .log_append(
      proj, step = "targets_disease_filter", id = scores$target_id[no_assoc],
      message = "targets_disease_no_association: Open Targets reports no evidence for this target-disease pair; dropped"
    )
  }
  scores <- scores[!no_assoc, , drop = FALSE]

  if (!is.null(min_score) && nrow(scores) > 0) {
    below <- scores$association_score < min_score
    if (any(below)) {
      proj <- .log_append(
        proj, step = "targets_disease_filter", id = scores$target_id[below],
        message = paste0(
          "targets_disease_below_min_score: association_score ",
          format(round(scores$association_score[below], 4), nsmall = 4),
          " < min_score ", min_score
        )
      )
    }
    scores <- scores[!below, , drop = FALSE]
  }

  ## fan back out from one row per UniProt ID to one row per (compound, uniprot)
  merged <- merge(
    imported[, c("compound_id", "uniprot_id")], scores,
    by.x = "uniprot_id", by.y = "target_id", all = FALSE
  )
  result <- data.frame(
    compound_id = merged$compound_id, target_id = merged$uniprot_id,
    disease_id = merged$disease_id, association_score = merged$association_score,
    evidence = merged$evidence, stringsAsFactors = FALSE
  )
  result <- unique(result)

  patliRResults(proj, "targets_disease") <- result
  .write_results_csv(proj, "targets_disease", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.empty_targets_disease_scores <- function() {
  data.frame(target_id = character(), disease_id = character(),
             association_score = double(), evidence = character(), stringsAsFactors = FALSE)
}

#' Sanitize a free-text term into a safe, collision-resistant cache-file-name
#' fragment
#' @details
#' `gsub("[^A-Za-z0-9]+", "_", x)` alone collapses every run of punctuation
#' to a single `_`, so distinct raw strings that differ only in punctuation
#' -- `"FLO-ET"` vs `"FLO_ET"`, `"EFLO S"` vs `"EFLO-S"` -- collapse to the
#' *same* slug and silently share a cache file. A short hash of the raw
#' string (computed before the punctuation is collapsed, so it still
#' differs when the collapsed slug does not) is appended to disambiguate
#' them; no new dependency (`digest` is not otherwise used anywhere in this
#' package).
#' @return `character`, same length as `x`.
#' @keywords internal
.cache_key_slug <- function(x) {
  vapply(x, function(xi) {
    slug <- gsub("[^A-Za-z0-9]+", "_", xi)
    codes <- utf8ToInt(enc2utf8(xi))
    h <- if (length(codes) == 0) 0 else sum(as.numeric(codes) * seq_along(codes)) %% (2^31 - 1)
    paste0(slug, "_", sprintf("%08x", as.integer(h)))
  }, character(1), USE.NAMES = FALSE)
}

#' POST one query to the Open Targets Platform GraphQL API (v4, no API key)
#'
#' `disease_genes_fetch()` fires this ~`ceiling(count / 500)` times in a
#' paginated loop (20 requests for a large disease) against a rate-limiting
#' API, so it carries a User-Agent, a ~3 req/s throttle and exponential
#' backoff on 429/5xx -- same treatment as `.refdb_get_json()`.
#'
#' @return The `data` field of the parsed JSON response (a nested list).
#' @keywords internal
.open_targets_graphql <- function(query_string, variables = list()) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query Open Targets.")
  }
  req <- httr2::request("https://api.platform.opentargets.org/api/v4/graphql")
  req <- httr2::req_user_agent(req, paste0("patliR/", .patliR_version()))
  req <- httr2::req_throttle(req, capacity = 3, fill_time_s = 1)
  req <- httr2::req_retry(
    req, max_tries = 4,
    is_transient = function(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504),
    backoff = function(tries) 2^tries
  )
  req <- httr2::req_body_json(req, list(query = query_string, variables = variables))
  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_json(resp)
  if (!is.null(body$errors) && length(body$errors) > 0) {
    msgs <- vapply(body$errors, function(e) e$message %||% "unknown error", character(1))
    stop(paste("Open Targets GraphQL error:", paste(msgs, collapse = "; ")))
  }
  body$data
}

#' Resolve a free-text term or existing ID to a canonical Open Targets ID
#'
#' Uses the `mapIds` endpoint (the same full-text mapping that backs Open
#' Targets' own search box), restricted to one entity type at a time
#' (`"target"` or `"disease"`) so an ambiguous term cannot silently resolve
#' to the wrong entity kind. Takes the top-ranked hit -- `mapIds` returns
#' hits already ordered by relevance, same as their `search()` endpoint.
#'
#' @return A character scalar (the resolved ID), or throws if there is no hit.
#' @keywords internal
.open_targets_map_id <- function(term, entity = c("target", "disease")) {
  entity <- match.arg(entity)
  data <- .open_targets_graphql(
    query_string = "
      query MapIds($terms: [String!]!, $entities: [String!]) {
        mapIds(queryTerms: $terms, entityNames: $entities) {
          mappings { term hits { id name entity } }
        }
      }
    ",
    variables = list(terms = list(term), entities = list(entity))
  )
  hits <- data$mapIds$mappings[[1]]$hits
  if (length(hits) == 0) {
    stop(paste0("Open Targets: no ", entity, " matched the term '", term, "'"))
  }
  hits[[1]]$id
}

#' Resolve a `disease` argument to a canonical Open Targets EFO/MONDO ID
#'
#' @description
#' `mapIds()` (used by [.open_targets_map_id()]) is a *full-text* search over
#' disease names/synonyms, so it should not be relied on to match a raw
#' ontology ID string -- so an ID-shaped string (EFO/MONDO/Orphanet/HP/DOID,
#' the common `PREFIX_alnum` shape) is looked up directly via
#' `disease(efoId: ...)` first, which both resolves and validates it in one
#' call, rather than routed through `mapIds()`.
#'
#' @section A failure mode this guards against:
#' `disease(efoId = "EFO_0000537")` (hypertension's old EFO ID) returns
#' `NULL` with no GraphQL error at all -- not because direct lookup is
#' wrong, but because `EFO_0000537` is **obsolete** (retired in EFO 3.88.0,
#' replaced by `MONDO_0005044`; some catalogs still link the old ID, Open
#' Targets does not).
#' An earlier version of this function swallowed that `NULL`/failure into a
#' bare fallback and only ever surfaced `mapIds()`'s generic "no match"
#' error, hiding the real cause -- fixed to report both failure reasons (see
#' below). Practically: prefer free text (always resolves to the *current*
#' term via `mapIds()`) over a hand-typed ontology ID unless you know it is
#' still current. If the ID-shaped fast path finds nothing, this still
#' falls back to `mapIds()` rather than failing outright, in case the
#' string was coincidentally ID-shaped free text.
#'
#' @return A character scalar (the resolved EFO ID), or throws if nothing
#'   matches either path. The error message always includes *both* failure
#'   reasons when the ID-shaped fast path was attempted -- see the note
#'   below on why swallowing the first one to a bare `NULL` was itself a bug.
#'
#' @section Diagnostic history:
#' An earlier version of this function caught the direct-lookup failure with
#' `tryCatch(..., error = function(e) NULL)` and silently fell through to
#' `mapIds()` -- so when *both* paths failed, the only error the caller ever
#' saw was `mapIds()`'s generic "no disease matched", with the real reason
#' the direct lookup failed thrown away. That is exactly the kind of
#' silently-swallowed failure `patliR` otherwise takes care to avoid (see
#' `.fetch_external()`) -- fixed here to keep and report both messages.
#' @keywords internal
.open_targets_resolve_disease <- function(term) {
  direct_err <- NULL
  if (grepl("^[A-Za-z]+_[A-Za-z0-9]+$", term)) {
    direct <- tryCatch(.open_targets_disease_by_id(term), error = function(e) e)
    if (!inherits(direct, "error") && !is.null(direct)) return(direct)
    direct_err <- if (inherits(direct, "error")) conditionMessage(direct) else "disease(efoId = ...) returned no match"
  }

  mapped <- tryCatch(.open_targets_map_id(term, entity = "disease"), error = function(e) e)
  if (!inherits(mapped, "error")) return(mapped)

  if (is.null(direct_err)) {
    stop(conditionMessage(mapped))
  }
  stop(paste0(
    "Open Targets: could not resolve disease '", term, "' -- ",
    "direct ID lookup failed (", direct_err, "); ",
    "mapIds() fallback also failed (", conditionMessage(mapped), ")"
  ))
}

#' @keywords internal
.open_targets_disease_by_id <- function(efo_id) {
  data <- .open_targets_graphql(
    query_string = "query DiseaseById($id: String!) { disease(efoId: $id) { id } }",
    variables = list(id = efo_id)
  )
  data$disease$id
}

#' Query the overall target-disease association score for one target/disease
#' pair, using `Bs` to filter server-side to that one disease
#'
#' @return A `list(score, disease_id, evidence)`; `score` is `NA_real_` (and
#'   `evidence` is `NA_character_`) if Open Targets has no association at all
#'   for this pair.
#' @keywords internal
.open_targets_target_disease_score <- function(ensembl_id, efo_id) {
  data <- .open_targets_graphql(
    query_string = "
      query Assoc($ensemblId: String!, $efoIds: [String!]) {
        target(ensemblId: $ensemblId) {
          associatedDiseases(Bs: $efoIds) {
            rows {
              score
              disease { id }
              datatypeScores { id score }
            }
          }
        }
      }
    ",
    variables = list(ensemblId = ensembl_id, efoIds = list(efo_id))
  )
  rows <- data$target$associatedDiseases$rows
  if (length(rows) == 0) {
    return(list(score = NA_real_, disease_id = efo_id, evidence = NA_character_))
  }
  r <- rows[[1]]
  dt <- r$datatypeScores
  evidence <- if (length(dt) > 0) {
    paste(vapply(dt, function(x) paste0(x$id, "=", round(x$score, 3)), character(1)), collapse = "; ")
  } else {
    NA_character_
  }
  list(score = r$score, disease_id = r$disease$id %||% efo_id, evidence = evidence)
}
