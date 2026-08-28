#' @include AllGenerics.R internal.R
NULL

## reference_compounds.csv / reference_bioactivity.csv are stored
## long/relational (no list-columns) so they stay plain CSV. See DESIGN.md.

#' Build the local reference database of compound identity and bioactivity
#'
#' @description
#' `refdb_build()` populates `reference_compounds.csv` and
#' `reference_bioactivity.csv` inside [projectDir()] for the compounds
#' currently in [compounds()], pulling from public sources. This reference
#' database is what later families (`bias_*`, `network_degeneracy()`) treat
#' as "the background" -- it grows across an entire `patliR` install if you
#' point several projects at the same `project_dir` root, or stays
#' project-local otherwise; that choice is yours.
#'
#' Only sources with a real, confirmed REST API are fetched automatically
#' (`"pubchem"`, `"chembl"`). `"coconut"` is accepted as a source name for
#' forward-compatibility but is not implemented yet in this version --
#' request it once `coconut_fetch()` (natural-product family) ships; asking
#' for it now raises an informative error rather than pretending to
#' succeed.
#'
#' ChEMBL identity is resolved by a deterministic identity chain, never a
#' similarity search: PubChem CID -> standard InChIKey (batched) -> ChEMBL
#' `standard_inchi_key` exact match (batched) -> InChIKey-skeleton match
#' (first 14 characters, flagged as weaker in `match_type`) -> SMILES
#' `flexmatch` -> give up. When a query returns several molecules the one
#' kept is chosen by [.chembl_pick()] (own-parent, then has a preferred
#' name, then lowest ChEMBL id), never simply the first.
#'
#' @inheritParams compounds
#' @param sources Character vector, any of `"pubchem"`, `"chembl"`,
#'   `"coconut"`.
#' @param compound_ids Character vector of `compounds(proj)$id` to include,
#'   or `NULL` (default) for all compounds currently in `proj`.
#' @param fetch_mode `"warn_and_cache"` (default) or `"abort"`. See
#'   `.fetch_external()` (internal): this is enrichment
#'   data, so failures do not stop the pipeline by default.
#'
#' @return The updated `proj`, with `reference_compounds` and
#'   `reference_bioactivity` entries in [patliRResults()], also written to
#'   `reference_compounds.csv` / `reference_bioactivity.csv` inside
#'   [projectDir()]. `reference_compounds` carries `inchikey`, `match_type`,
#'   `parent_chembl_id` and `fetched_at`; `reference_bioactivity` carries
#'   `standard_relation`, `pchembl_value`, `target_organism`, `assay_type`
#'   and `data_validity_comment`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- refdb_build(proj, sources = c("pubchem", "chembl")) # needs internet
#' }
#'
#' @export
refdb_build <- function(proj, sources = c("pubchem", "chembl"),
                         compound_ids = NULL,
                         fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"))
  ## "coconut" is a valid name for forward-compat but not implemented; it is
  ## NOT in the default (a bare `refdb_build(proj)` used to abort because
  ## match.arg(several.ok = TRUE) returned all of the default vector).
  sources <- match.arg(sources, choices = c("pubchem", "chembl", "coconut"), several.ok = TRUE)
  fetch_mode <- match.arg(fetch_mode)

  if ("coconut" %in% sources) {
    cli::cli_abort(c(
      "{.val coconut} is not implemented yet in this version of patliR.",
      "i" = "It will be handled by {.fn coconut_fetch} once that family ships; drop it from {.arg sources} for now."
    ))
  }

  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  if (nrow(cmp) == 0) {
    cli::cli_warn("No matching compounds in {.arg proj}; nothing to build.")
    return(proj)
  }

  smiles_key <- ifelse(!is.na(cmp$canonical_smiles) & nzchar(cmp$canonical_smiles),
                       cmp$canonical_smiles, cmp$smiles)
  cid_key <- as.character(cmp$pubchem_id)
  fetched_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  ## Build ONLY this run's rows, then merge -- so a second refdb_build()
  ## (or one after patliR_load() has restored the CSVs into the results
  ## bag) replaces this run's compounds instead of duplicating every row.
  new_compounds <- .empty_reference_compounds()
  new_bioactivity <- .empty_reference_bioactivity()
  n_resolved <- 0L
  n_weak <- 0L

  ## ---- ChEMBL identity, resolved once for the whole batch ---------------
  identity <- NULL
  if ("chembl" %in% sources) {
    identity <- .fetch_external(
      fetch_fun = function() .chembl_resolve(smiles_key, cid_key),
      cache_dir = cacheDir(proj),
      ## keyed by chemical identity (CIDs + SMILES), NOT by the project-local
      ## C0001 ids -- re-prep in a different order can't return stale rows
      ## for the wrong compound.
      cache_key = paste0("refdb_identity_", .refdb_batch_key(cid_key, smiles_key)),
      mode = fetch_mode
    )
  }

  for (i in seq_len(nrow(cmp))) {
    id <- cmp$id[i]

    if ("pubchem" %in% sources) {
      pc <- .fetch_external(
        fetch_fun = function() .pubchem_lookup(cid = cid_key[i], smiles = smiles_key[i]),
        cache_dir = cacheDir(proj),
        cache_key = paste0("refdb_pubchem_", .refdb_identity_key(cid_key[i], smiles_key[i])),
        mode = fetch_mode
      )
      if (!is.null(pc)) {
        new_compounds <- rbind(new_compounds, data.frame(
          compound_id = id, source = "pubchem", external_id = pc$cid,
          name = pc$name %||% NA_character_, inchikey = pc$inchikey %||% NA_character_,
          match_type = NA_character_, parent_chembl_id = NA_character_,
          fetched_at = fetched_at, stringsAsFactors = FALSE
        ))
      }
    }

    if ("chembl" %in% sources && !is.null(identity)) {
      ch <- identity[identity$row == i, , drop = FALSE]
      if (nrow(ch) == 1 && !is.na(ch$chembl_id)) {
        n_resolved <- n_resolved + 1L
        if (!identical(ch$match_type, "inchikey_exact")) n_weak <- n_weak + 1L
        new_compounds <- rbind(new_compounds, data.frame(
          compound_id = id, source = "chembl", external_id = ch$chembl_id,
          name = ch$pref_name %||% NA_character_, inchikey = ch$inchikey %||% NA_character_,
          match_type = ch$match_type, parent_chembl_id = ch$parent_chembl_id %||% NA_character_,
          fetched_at = fetched_at, stringsAsFactors = FALSE
        ))
        bio <- .fetch_external(
          fetch_fun = function() .chembl_bioactivity(ch$chembl_id),
          cache_dir = cacheDir(proj),
          ## keyed by the resolved ChEMBL id, not the C0001 id.
          cache_key = paste0("refdb_chembl_bioactivity_", ch$chembl_id),
          mode = fetch_mode
        )
        if (!is.null(bio) && nrow(bio) > 0) {
          bio$compound_id <- id
          new_bioactivity <- rbind(new_bioactivity, bio[, names(new_bioactivity)])
        }
      }
    }
  }

  ref_compounds <- .refdb_merge(
    patliRResults(proj, "reference_compounds") %||% .empty_reference_compounds(),
    new_compounds, cmp$id, key = c("compound_id", "source"))
  ref_bioactivity <- .refdb_merge(
    patliRResults(proj, "reference_bioactivity") %||% .empty_reference_bioactivity(),
    new_bioactivity, cmp$id,
    key = c("compound_id", "target_chembl_id", "standard_type", "assay_chembl_id"))

  patliRResults(proj, "reference_compounds") <- ref_compounds
  patliRResults(proj, "reference_bioactivity") <- ref_bioactivity
  .write_results_csv(proj, "reference_compounds", ref_compounds)
  .write_results_csv(proj, "reference_bioactivity", ref_bioactivity)
  proj <- .log_append(
    proj, step = "refdb_build", id = NA_character_,
    message = paste0(
      "built reference DB for ", nrow(cmp), " compound(s): ",
      nrow(new_compounds), " identity row(s), ", nrow(new_bioactivity),
      " bioactivity row(s); ChEMBL resolved ", n_resolved, "/", nrow(cmp),
      " (", n_weak, " via a weaker match_type: skeleton/flexmatch)"
    )
  )
  .write_log_csv(proj)
  proj
}

#' Replace this run's rows in a reference table, keep every other
#' compound's, drop exact-key duplicates. Makes [refdb_build()] idempotent,
#' which the "every step is proj -> proj" contract needs. Tolerates an
#' `old` table written by an earlier schema (missing columns are filled
#' with `NA` so the two frames can be `rbind()`-ed).
#' @keywords internal
.refdb_merge <- function(old, new, touched_ids, key) {
  if (nrow(old) == 0) return(new)
  old <- old[!old$compound_id %in% touched_ids, , drop = FALSE]
  for (col in setdiff(names(new), names(old))) old[[col]] <- rep(NA, nrow(old))
  for (col in setdiff(names(old), names(new))) new[[col]] <- rep(NA, nrow(new))
  out <- rbind(old, new[, names(old), drop = FALSE])
  key <- intersect(key, names(out))
  if (length(key) > 0) out <- out[!duplicated(out[, key, drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Incrementally update the reference database for new compounds
#'
#' Same fetch logic as [refdb_build()], restricted to (and additive for)
#' the given compounds -- use this after adding new compounds with
#' [prep_compounds()] instead of rebuilding the whole reference database
#' from scratch.
#'
#' @inheritParams refdb_build
#' @param compound_ids Character vector of `compounds(proj)$id` to fetch.
#'   Required (unlike in [refdb_build()], there is no "all" default, to
#'   avoid accidentally re-fetching everything).
#'
#' @return The updated `proj`. See [refdb_build()].
#' @export
refdb_update <- function(proj, compound_ids,
                          sources = c("pubchem", "chembl"),
                          fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"), is.character(compound_ids), length(compound_ids) > 0)
  refdb_build(proj, sources = sources, compound_ids = compound_ids, fetch_mode = fetch_mode)
}

#' Rebuild the derived, performance-only reference database cache
#'
#' @description
#' Regenerates the `.rds` cache of the reference database inside
#' [cacheDir()] from `results/reference_compounds.csv` /
#' `results/reference_bioactivity.csv` in [projectDir()]. This cache is
#' purely a performance optimization for
#' downstream functions that repeatedly join against the reference
#' database (e.g. `bias_audit()`) -- it is always safe to delete, and this
#' function always regenerates it from the CSVs, never the other way
#' around.
#'
#' @inheritParams compounds
#' @return `invisible(NULL)`. Called for its side effect (writing the
#'   cache file).
#' @export
refdb_rebuild_cache <- function(proj) {
  stopifnot(is(proj, "PatliRProject"))
  results_dir <- file.path(projectDir(proj), "results")
  ref_compounds_path <- file.path(results_dir, "reference_compounds.csv")
  ref_bioactivity_path <- file.path(results_dir, "reference_bioactivity.csv")

  if (!file.exists(ref_compounds_path) && !file.exists(ref_bioactivity_path)) {
    cli::cli_warn("No reference database CSVs found in {.path {results_dir}}; run {.fn refdb_build} first.")
    return(invisible(NULL))
  }

  cache <- list(
    reference_compounds = if (file.exists(ref_compounds_path)) utils::read.csv(ref_compounds_path, stringsAsFactors = FALSE) else .empty_reference_compounds(),
    reference_bioactivity = if (file.exists(ref_bioactivity_path)) utils::read.csv(ref_bioactivity_path, stringsAsFactors = FALSE) else .empty_reference_bioactivity(),
    built_at = Sys.time()
  )
  saveRDS(cache, file.path(cacheDir(proj), "refdb_cache.rds"))
  invisible(NULL)
}

#' @keywords internal
.empty_reference_compounds <- function() {
  data.frame(compound_id = character(), source = character(), external_id = character(),
             name = character(), inchikey = character(), match_type = character(),
             parent_chembl_id = character(), fetched_at = character(),
             stringsAsFactors = FALSE)
}

#' @keywords internal
.empty_reference_bioactivity <- function() {
  data.frame(compound_id = character(), target_chembl_id = character(), target_name = character(),
             standard_type = character(), standard_value = double(), standard_units = character(),
             assay_chembl_id = character(), standard_relation = character(),
             pchembl_value = double(), target_organism = character(),
             assay_type = character(), data_validity_comment = character(),
             stringsAsFactors = FALSE)
}

## ---- HTTP plumbing (shared by all three refdb fetchers) -----------------

#' `User-Agent` string for patliR's outbound reference-database calls
#' @keywords internal
.refdb_user_agent <- function() {
  paste0("patliR/", .patliR_version(),
         " (network pharmacology; +https://github.com/hierax00/patliR)")
}

#' GET a URL and parse the JSON body, with a `User-Agent`, a ~3 req/s
#' throttle, and exponential-backoff retry on 429 / 5xx
#'
#' @description
#' The single choke point every PubChem/ChEMBL request in this file goes
#' through, so throttle/retry/User-Agent are applied once and tests can
#' mock the network by stubbing this one binding.
#' @param url Character scalar.
#' @param query Optional named list of query parameters; `httr2` handles
#'   the percent-encoding (this is how SMILES reaches `flexmatch` without
#'   ever touching the URL path).
#' @return The parsed JSON body (a list).
#' @keywords internal
.refdb_get_json <- function(url, query = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query external reference databases.")
  }
  req <- httr2::request(url)
  req <- httr2::req_user_agent(req, .refdb_user_agent())
  req <- httr2::req_throttle(req, capacity = 3, fill_time_s = 1)
  req <- httr2::req_retry(
    req, max_tries = 4,
    is_transient = function(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504),
    backoff = function(tries) 2^tries
  )
  if (length(query) > 0) req <- httr2::req_url_query(req, !!!query)
  httr2::resp_body_json(httr2::req_perform(req))
}

#' Short, order-independent digest of a set of chemical identifiers, used
#' to key the batch identity cache by *what* was asked for, not the run.
#' @keywords internal
.refdb_batch_key <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  substr(rlang::hash(paste(sort(unique(parts)), collapse = "|")), 1, 16)
}

#' Per-compound identity cache key: the CID if we have one, else a digest
#' of the SMILES; never the project-local id.
#' @keywords internal
.refdb_identity_key <- function(cid, smiles) {
  if (!is.na(cid) && nzchar(cid)) return(paste0("cid", cid))
  if (!is.na(smiles) && nzchar(smiles)) return(paste0("smi", substr(rlang::hash(smiles), 1, 16)))
  "none"
}

## ---- PubChem ------------------------------------------------------------

#' Look up a compound on PubChem, by CID if we have one, else by SMILES
#'
#' `patliR` has no locally computed InChIKey (see `.check_structures()`,
#' internal), so this prefers the exact `pubchem_id` (CID) already on the
#' compound when available (Scenario B input), and only falls back to a
#' SMILES-based lookup -- exact-structure but canonicalization-sensitive --
#' when there is no CID. Returns the standard InChIKey and connectivity
#' SMILES alongside the title so `reference_compounds` has a real
#' structure-identity column.
#' @keywords internal
.pubchem_lookup <- function(cid = NA_character_, smiles = NA_character_) {
  props <- c("Title", "InChIKey", "ConnectivitySMILES")
  if (!is.na(cid) && nzchar(cid)) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid,
                  "/property/", paste(props, collapse = ","), "/JSON")
    body <- .refdb_get_json(url)
  } else if (!is.na(smiles) && nzchar(smiles)) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/property/",
                  paste(props, collapse = ","), "/JSON")
    body <- .refdb_get_json(url, query = list(smiles = smiles))
  } else {
    stop("need either a PubChem CID or a SMILES to look up a compound on PubChem")
  }
  if (!is.null(body$Fault)) {
    stop("PubChem fault: ", body$Fault$Message %||% "unknown", " (", body$Fault$Code %||% "?", ")")
  }
  p <- body$PropertyTable$Properties[[1]]
  list(
    cid = as.character(p$CID),
    name = p$Title %||% NA_character_,
    inchikey = p$InChIKey %||% NA_character_,
    connectivity_smiles = p$ConnectivitySMILES %||% NA_character_
  )
}

#' Batched PubChem CID -> standard InChIKey (+ connectivity SMILES, title)
#' @param cids Character vector of CIDs (deduplicated internally).
#' @return Named list, `as.character(CID)` -> `list(cid, name, inchikey,
#'   connectivity_smiles)`.
#' @keywords internal
.pubchem_cids_to_identity <- function(cids, chunk = 100L) {
  cids <- unique(cids[!is.na(cids) & nzchar(cids)])
  out <- list()
  if (length(cids) == 0) return(out)
  groups <- split(cids, ceiling(seq_along(cids) / chunk))
  for (grp in groups) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
      paste(grp, collapse = ","),
      "/property/Title,InChIKey,ConnectivitySMILES/JSON"
    )
    body <- .refdb_get_json(url)
    if (!is.null(body$Fault)) {
      stop("PubChem fault: ", body$Fault$Message %||% "unknown")
    }
    for (p in body$PropertyTable$Properties) {
      out[[as.character(p$CID)]] <- list(
        cid = as.character(p$CID),
        name = p$Title %||% NA_character_,
        inchikey = p$InChIKey %||% NA_character_,
        connectivity_smiles = p$ConnectivitySMILES %||% NA_character_
      )
    }
  }
  out
}

## ---- ChEMBL identity chain --------------------------------------------

#' Resolve a batch of compounds to ChEMBL molecules by a deterministic
#' identity chain
#'
#' @description
#' Never a similarity search. The chain, applied to whatever is still
#' unresolved after each step:
#' \enumerate{
#'   \item PubChem CID -> standard InChIKey (batched), or a per-compound
#'     SMILES lookup when there is no CID;
#'   \item ChEMBL `molecule_structures__standard_inchi_key__in` exact
#'     match (batched) -> `match_type = "inchikey_exact"`;
#'   \item ChEMBL `..._startswith` on the 14-character InChIKey skeleton
#'     -> `match_type = "inchikey_skeleton"` (weaker: same 2D skeleton,
#'     possibly a different stereo/isotope/charge form);
#'   \item ChEMBL `molecule_structures__canonical_smiles__flexmatch` on the
#'     SMILES as a query parameter -> `match_type = "smiles_flexmatch"`.
#' }
#' When any query returns more than one molecule, [.chembl_pick()] chooses.
#'
#' @param smiles,cids Character vectors, same length (one entry per
#'   compound, in `compounds()` order). `NA`/`""` allowed.
#' @return A `data.frame` with one row per input compound: `row` (integer
#'   position), `inchikey`, `chembl_id`, `parent_chembl_id`, `pref_name`,
#'   `match_type`. `chembl_id` is `NA` for compounds that fell through the
#'   whole chain.
#' @keywords internal
.chembl_resolve <- function(smiles, cids) {
  n <- length(smiles)
  res <- data.frame(
    row = seq_len(n),
    inchikey = NA_character_, chembl_id = NA_character_,
    parent_chembl_id = NA_character_, pref_name = NA_character_,
    match_type = NA_character_, stringsAsFactors = FALSE
  )

  ## (1) InChIKeys
  cid_map <- .pubchem_cids_to_identity(cids)
  for (i in seq_len(n)) {
    if (!is.na(cids[i]) && nzchar(cids[i]) && !is.null(cid_map[[cids[i]]])) {
      res$inchikey[i] <- cid_map[[cids[i]]]$inchikey %||% NA_character_
    } else if (!is.na(smiles[i]) && nzchar(smiles[i])) {
      pc <- tryCatch(.pubchem_lookup(smiles = smiles[i]), error = function(e) NULL)
      if (!is.null(pc)) res$inchikey[i] <- pc$inchikey %||% NA_character_
    }
  }

  ## (2) exact InChIKey match, batched
  todo <- which(is.na(res$chembl_id) & !is.na(res$inchikey))
  if (length(todo) > 0) {
    exact <- .chembl_molecules_by_inchikey(unique(res$inchikey[todo]), filter = "exact")
    for (i in todo) {
      hit <- .chembl_pick(exact[!is.na(exact$query_inchikey) &
                                  exact$query_inchikey == res$inchikey[i], , drop = FALSE])
      if (!is.null(hit)) {
        res$chembl_id[i] <- hit$molecule_chembl_id
        res$parent_chembl_id[i] <- hit$parent_chembl_id
        res$pref_name[i] <- hit$pref_name
        res$match_type[i] <- "inchikey_exact"
      }
    }
  }

  ## (3) InChIKey-skeleton fallback (weaker), one query per skeleton
  todo <- which(is.na(res$chembl_id) & !is.na(res$inchikey))
  for (i in todo) {
    skel <- substr(res$inchikey[i], 1, 14)
    if (!nzchar(skel)) next
    hits <- .chembl_molecules_by_inchikey(skel, filter = "skeleton")
    hit <- .chembl_pick(hits)
    if (!is.null(hit)) {
      res$chembl_id[i] <- hit$molecule_chembl_id
      res$parent_chembl_id[i] <- hit$parent_chembl_id
      res$pref_name[i] <- hit$pref_name
      res$match_type[i] <- "inchikey_skeleton"
    }
  }

  ## (4) SMILES flexmatch, one query per compound
  todo <- which(is.na(res$chembl_id) & !is.na(smiles) & nzchar(smiles))
  for (i in todo) {
    hits <- .chembl_molecules_by_smiles_flexmatch(smiles[i])
    hit <- .chembl_pick(hits)
    if (!is.null(hit)) {
      res$chembl_id[i] <- hit$molecule_chembl_id
      res$parent_chembl_id[i] <- hit$parent_chembl_id
      res$pref_name[i] <- hit$pref_name
      res$match_type[i] <- "smiles_flexmatch"
    }
  }

  res
}

#' ChEMBL base URL for the data API.
#' @keywords internal
.chembl_base <- function() "https://www.ebi.ac.uk/chembl/api/data"

## `only=` field list for molecule identity queries (internal constant).
.chembl_molecule_only <- "molecule_chembl_id,pref_name,molecule_hierarchy,molecule_type,structure_type,molecule_structures"

#' One page of `molecule.json` -> list of flat molecule records
#' @keywords internal
.chembl_molecule_records <- function(query) {
  url <- paste0(.chembl_base(), "/molecule.json")
  recs <- list()
  repeat {
    body <- .refdb_get_json(url, query = query)
    query <- NULL # the `next` link already carries every parameter
    for (m in body$molecules %||% list()) {
      recs[[length(recs) + 1L]] <- list(
        molecule_chembl_id = m$molecule_chembl_id %||% NA_character_,
        pref_name = m$pref_name %||% NA_character_,
        parent_chembl_id = m$molecule_hierarchy$parent_chembl_id %||% m$molecule_chembl_id %||% NA_character_,
        molecule_type = m$molecule_type %||% NA_character_,
        structure_type = m$structure_type %||% NA_character_,
        standard_inchi_key = m$molecule_structures$standard_inchi_key %||% NA_character_
      )
    }
    nxt <- body$page_meta$`next`
    if (is.null(nxt) || is.na(nxt) || !nzchar(nxt)) break
    url <- paste0("https://www.ebi.ac.uk", nxt)
  }
  recs
}

#' ChEMBL molecule lookup by standard InChIKey (exact `__in`, batched) or
#' by 14-char skeleton (`__startswith`, one key at a time)
#'
#' @param inchikeys Character vector. For `filter = "exact"` these are full
#'   27-char standard InChIKeys and are queried in `__in` batches; for
#'   `filter = "skeleton"` each is a 14-char skeleton queried on its own.
#' @return A `data.frame` with `query_inchikey`, `molecule_chembl_id`,
#'   `pref_name`, `parent_chembl_id`, `molecule_type`, `structure_type`.
#' @keywords internal
.chembl_molecules_by_inchikey <- function(inchikeys, filter = c("exact", "skeleton"),
                                           chunk = 20L) {
  filter <- match.arg(filter)
  inchikeys <- unique(inchikeys[!is.na(inchikeys) & nzchar(inchikeys)])
  if (length(inchikeys) == 0) return(.empty_chembl_molecule_df())
  field <- if (filter == "exact") {
    "molecule_structures__standard_inchi_key__in"
  } else {
    "molecule_structures__standard_inchi_key__startswith"
  }
  groups <- if (filter == "exact") {
    split(inchikeys, ceiling(seq_along(inchikeys) / chunk))
  } else {
    as.list(inchikeys)
  }
  frames <- lapply(groups, function(grp) {
    q <- list(only = .chembl_molecule_only, limit = 1000L)
    q[[field]] <- if (filter == "exact") paste(grp, collapse = ",") else grp
    recs <- .chembl_molecule_records(q)
    if (length(recs) == 0) return(.empty_chembl_molecule_df())
    df <- do.call(rbind, lapply(recs, function(r) {
      data.frame(
        query_inchikey = r$standard_inchi_key,
        molecule_chembl_id = r$molecule_chembl_id, pref_name = r$pref_name,
        parent_chembl_id = r$parent_chembl_id, molecule_type = r$molecule_type,
        structure_type = r$structure_type, stringsAsFactors = FALSE
      )
    }))
    ## for a skeleton query every returned row shares the queried skeleton
    if (filter == "skeleton") df$query_inchikey <- grp
    df
  })
  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out
}

#' ChEMBL molecule lookup by SMILES flexmatch (SMILES as a query
#' parameter, never a path segment)
#' @keywords internal
.chembl_molecules_by_smiles_flexmatch <- function(smiles) {
  if (is.na(smiles) || !nzchar(smiles)) return(.empty_chembl_molecule_df())
  recs <- .chembl_molecule_records(list(
    molecule_structures__canonical_smiles__flexmatch = smiles,
    only = .chembl_molecule_only, limit = 1000L
  ))
  if (length(recs) == 0) return(.empty_chembl_molecule_df())
  out <- do.call(rbind, lapply(recs, function(r) {
    data.frame(
      query_inchikey = NA_character_,
      molecule_chembl_id = r$molecule_chembl_id, pref_name = r$pref_name,
      parent_chembl_id = r$parent_chembl_id, molecule_type = r$molecule_type,
      structure_type = r$structure_type, stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' @keywords internal
.empty_chembl_molecule_df <- function() {
  data.frame(query_inchikey = character(), molecule_chembl_id = character(),
             pref_name = character(), parent_chembl_id = character(),
             molecule_type = character(), structure_type = character(),
             stringsAsFactors = FALSE)
}

#' Deterministically choose one ChEMBL molecule when a query returns several
#'
#' @description
#' Order of preference: (1) the molecule is its own parent (not a salt,
#' mixture or other child form); (2) it has a `pref_name`; (3) lowest
#' numeric ChEMBL id. Never "the first row the API happened to return".
#' @param mols A `data.frame` as returned by the `.chembl_molecules_by_*`
#'   helpers, or `NULL`/0-row.
#' @return A one-row `list` (`molecule_chembl_id`, `pref_name`,
#'   `parent_chembl_id`, ...), or `NULL` if `mols` is empty.
#' @keywords internal
.chembl_pick <- function(mols) {
  if (is.null(mols) || nrow(mols) == 0) return(NULL)
  own_parent <- !is.na(mols$parent_chembl_id) &
    mols$parent_chembl_id == mols$molecule_chembl_id
  has_pref <- !is.na(mols$pref_name) & nzchar(mols$pref_name)
  num_id <- suppressWarnings(as.integer(sub("^CHEMBL", "", mols$molecule_chembl_id)))
  num_id[is.na(num_id)] <- .Machine$integer.max
  ord <- order(!own_parent, !has_pref, num_id)
  as.list(mols[ord[1], , drop = FALSE])
}

## ---- ChEMBL bioactivity ----------------------------------------------

## `target_pref_name` values that are not real biological targets and must
## not pollute `reference_bioactivity` (they would be counted as targets by
## `bias_audit()`). Internal constant.
.chembl_nontarget_names <- c(
  "No relevant target", "Unchecked", "Unspecified", "Non-molecular",
  "Log S", "LogD7.4", "LogP", "Solubility", "ADMET", "ADME"
)

#' Fetch every bioactivity row ChEMBL has for a molecule
#'
#' @description
#' Fully paginated via `page_meta$next` (no hard `limit`), field-filtered
#' with `only=`, non-target rows dropped, deterministically ordered.
#' `standard_relation` and `data_validity_comment` are kept as columns so a
#' censored / flagged measurement (`IC50 > 100000`, "Outside typical
#' range") is never stored as if it were an exact, clean value.
#'
#' @param chembl_id A single `molecule_chembl_id`.
#' @return A `data.frame` with the `reference_bioactivity` columns except
#'   `compound_id` (added by the caller).
#' @keywords internal
.chembl_bioactivity <- function(chembl_id) {
  only <- paste(c(
    "target_chembl_id", "target_pref_name", "target_organism",
    "standard_type", "standard_relation", "standard_value", "standard_units",
    "pchembl_value", "assay_chembl_id", "assay_type", "data_validity_comment"
  ), collapse = ",")
  query <- list(molecule_chembl_id = chembl_id, only = only, limit = 1000L)
  url <- paste0(.chembl_base(), "/activity.json")

  acts <- list()
  repeat {
    body <- .refdb_get_json(url, query = query)
    query <- NULL
    acts <- c(acts, body$activities %||% list())
    nxt <- body$page_meta$`next`
    if (is.null(nxt) || is.na(nxt) || !nzchar(nxt)) break
    url <- paste0("https://www.ebi.ac.uk", nxt)
  }

  empty <- .empty_reference_bioactivity()[, -1, drop = FALSE]
  if (length(acts) == 0) return(empty)

  rows <- lapply(acts, function(a) {
    data.frame(
      target_chembl_id = a$target_chembl_id %||% NA_character_,
      target_name = a$target_pref_name %||% NA_character_,
      standard_type = a$standard_type %||% NA_character_,
      standard_value = suppressWarnings(as.double(a$standard_value %||% NA)),
      standard_units = a$standard_units %||% NA_character_,
      assay_chembl_id = a$assay_chembl_id %||% NA_character_,
      standard_relation = a$standard_relation %||% NA_character_,
      pchembl_value = suppressWarnings(as.double(a$pchembl_value %||% NA)),
      target_organism = a$target_organism %||% NA_character_,
      assay_type = a$assay_type %||% NA_character_,
      data_validity_comment = a$data_validity_comment %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  ## drop non-target rows: no target id at all, or a placeholder name
  keep <- !is.na(out$target_chembl_id) & nzchar(out$target_chembl_id) &
    !(out$target_name %in% .chembl_nontarget_names)
  out <- out[keep, , drop = FALSE]
  if (nrow(out) == 0) return(empty)

  out <- out[order(out$target_chembl_id, out$standard_type, out$assay_chembl_id,
                   out$standard_value, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
