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
#'   [projectDir()].
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
refdb_build <- function(proj, sources = c("pubchem", "chembl", "coconut"),
                         compound_ids = NULL,
                         fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"))
  sources <- match.arg(sources, several.ok = TRUE)
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

  ref_compounds <- patliRResults(proj, "reference_compounds") %||% .empty_reference_compounds()
  ref_bioactivity <- patliRResults(proj, "reference_bioactivity") %||% .empty_reference_bioactivity()

  for (i in seq_len(nrow(cmp))) {
    id <- cmp$id[i]
    pubchem_cid <- cmp$pubchem_id[i]
    smiles_key <- cmp$canonical_smiles[i] %||% cmp$smiles[i]

    if ("pubchem" %in% sources) {
      pc <- .fetch_external(
        fetch_fun = function() .pubchem_lookup(cid = pubchem_cid, smiles = smiles_key),
        cache_dir = cacheDir(proj), cache_key = paste0("refdb_pubchem_", id), mode = fetch_mode
      )
      if (!is.null(pc)) {
        ref_compounds <- rbind(ref_compounds, data.frame(
          compound_id = id, source = "pubchem", external_id = pc$cid,
          name = pc$name %||% NA_character_, stringsAsFactors = FALSE
        ))
      }
    }

    if ("chembl" %in% sources) {
      ch <- .fetch_external(
        fetch_fun = function() .chembl_lookup(smiles_key),
        cache_dir = cacheDir(proj), cache_key = paste0("refdb_chembl_", id), mode = fetch_mode
      )
      if (!is.null(ch)) {
        ref_compounds <- rbind(ref_compounds, data.frame(
          compound_id = id, source = "chembl", external_id = ch$chembl_id,
          name = ch$pref_name %||% NA_character_, stringsAsFactors = FALSE
        ))
        bio <- .fetch_external(
          fetch_fun = function() .chembl_bioactivity(ch$chembl_id),
          cache_dir = cacheDir(proj), cache_key = paste0("refdb_chembl_bioactivity_", id), mode = fetch_mode
        )
        if (!is.null(bio) && nrow(bio) > 0) {
          bio$compound_id <- id
          ref_bioactivity <- rbind(ref_bioactivity, bio[, names(ref_bioactivity)])
        }
      }
    }
  }

  patliRResults(proj, "reference_compounds") <- ref_compounds
  patliRResults(proj, "reference_bioactivity") <- ref_bioactivity
  .write_results_csv(proj, "reference_compounds", ref_compounds)
  .write_results_csv(proj, "reference_bioactivity", ref_bioactivity)
  proj
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
             name = character(), stringsAsFactors = FALSE)
}

#' @keywords internal
.empty_reference_bioactivity <- function() {
  data.frame(compound_id = character(), target_chembl_id = character(), target_name = character(),
             standard_type = character(), standard_value = double(), standard_units = character(),
             assay_chembl_id = character(), stringsAsFactors = FALSE)
}

#' Look up a compound on PubChem, by CID if we have one, else by SMILES
#'
#' `patliR` has no true InChIKey (see `.check_structures()`, internal), so
#' this prefers the exact `pubchem_id` (CID) already on the compound when
#' available (Scenario B input), and only falls back to a SMILES-based
#' lookup -- which is exact-structure but canonicalization-sensitive -- when
#' there is no CID.
#' @keywords internal
.pubchem_lookup <- function(cid = NA_character_, smiles = NA_character_) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query PubChem.")
  }
  if (!is.na(cid) && nzchar(cid)) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/property/Title/JSON")
  } else if (!is.na(smiles) && nzchar(smiles)) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
      utils::URLencode(smiles, reserved = TRUE), "/property/Title/JSON"
    )
  } else {
    stop("need either a PubChem CID or a SMILES to look up a compound on PubChem")
  }
  resp <- httr2::req_perform(httr2::request(url))
  body <- httr2::resp_body_json(resp)
  props <- body$PropertyTable$Properties[[1]]
  list(cid = as.character(props$CID), name = props$Title)
}

#' Percent-encode a SMILES string for use as a ChEMBL URL path segment
#'
#' @description
#' ChEMBL's structure-search endpoints take the SMILES as a literal path
#' segment. Only `#`, `/`, and `%` need escaping (they change how the URL
#' itself is parsed); reserved-encoding `=`/`(`/`)` makes the endpoint
#' return `HTTP 500`. `%` must be escaped first, or it corrupts the `%23`/
#' `%2F` escapes this function then adds.
#' @return Character scalar, safe to paste into a ChEMBL path-segment URL.
#' @keywords internal
.chembl_url_encode_smiles <- function(smiles) {
  smiles <- gsub("%", "%25", smiles, fixed = TRUE) # must run first
  smiles <- gsub("#", "%23", smiles, fixed = TRUE)
  smiles <- gsub("/", "%2F", smiles, fixed = TRUE)
  utils::URLencode(smiles, reserved = FALSE) # handles spaces / non-ASCII, if any
}

#' Look up a compound on ChEMBL by SMILES, using ChEMBL's 100% similarity search
#'
#' @description
#' Uses `/similarity/<smiles>/100.json` (the
#' `molecule_structures__canonical_smiles__flexmatch` filter returns
#' `HTTP 500` on the current API). Similarity search matches the query
#' regardless of which toolkit canonicalized it, which is the cross-toolkit
#' problem patliR has without an InChIKey (see `.check_structures()`). A
#' 100% match can return more than one entry (e.g. a free base and its
#' salt); the first is used.
#' @keywords internal
.chembl_lookup <- function(smiles) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query ChEMBL.")
  }
  if (is.na(smiles) || !nzchar(smiles)) stop("need a SMILES to look up a compound on ChEMBL")
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/similarity/",
    .chembl_url_encode_smiles(smiles), "/100.json"
  )
  resp <- httr2::req_perform(httr2::request(url))
  body <- httr2::resp_body_json(resp)
  if (length(body$molecules) == 0) stop("no ChEMBL molecule found for this SMILES (100% similarity search)")
  m <- body$molecules[[1]]
  list(chembl_id = m$molecule_chembl_id, pref_name = m$pref_name)
}

#' @keywords internal
.chembl_bioactivity <- function(chembl_id) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query ChEMBL.")
  }
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=", chembl_id,
    "&limit=50"
  )
  resp <- httr2::req_perform(httr2::request(url))
  body <- httr2::resp_body_json(resp)
  acts <- body$activities
  if (length(acts) == 0) return(.empty_reference_bioactivity()[, -1])
  rows <- lapply(acts, function(a) {
    data.frame(
      target_chembl_id = a$target_chembl_id %||% NA_character_,
      target_name = a$target_pref_name %||% NA_character_,
      standard_type = a$standard_type %||% NA_character_,
      standard_value = suppressWarnings(as.double(a$standard_value %||% NA)),
      standard_units = a$standard_units %||% NA_character_,
      assay_chembl_id = a$assay_chembl_id %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
