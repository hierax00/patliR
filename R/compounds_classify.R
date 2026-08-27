#' @include AllGenerics.R internal.R
NULL

## compounds_classify() -- chemical-family classification via NPClassifier
## (Kim et al. 2021, J Nat Prod 84, 2795-2807), a free REST API purpose-
## built for natural products (GNPS/UCSD): given a SMILES it returns a
## pathway (Alkaloids, Amino acids and Peptides, Carbohydrates, Fatty
## acids, Polyketides, Shikimates and Phenylpropanoids, Terpenoids), a
## superclass, and a class. Chosen over the more general ClassyFire
## (broader chemical taxonomy, not natural-product-specific) because
## patliR's whole domain is natural product extracts -- NPClassifier's
## pathway categories are exactly "familias quimicas" the way a
## phytochemist would name them, decided with Uriel in chat (2026-08-11).
## Same .fetch_external()/httr2 pattern already established in refdb.R --
## no new Suggests needed.

#' Classify compounds into natural-product chemical families (NPClassifier)
#'
#' @description
#' Queries the free, no-auth NPClassifier REST API
#' (<https://npclassifier.gnps2.org>, Kim, H.W. et al. (2021), "NPClassifier:
#' A Deep Neural Network-Based Structural Classification Tool for Natural
#' Products", *J. Nat. Prod.* 84, 2795-2807, \doi{10.1021/acs.jnatprod.1c00399})
#' for every compound in [compounds()] (or a subset), by canonical SMILES.
#' This is what [plot_chemical_space()]'s `color_by = "family"` uses to
#' color points by chemical family (alkaloids, terpenoids, phenylpropanoids,
#' ...) -- the natural-product analogue of coloring by molecule category in
#' the classic "chemical space" figures (e.g. Reymond & Awale, 2012).
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound currently in [compounds()].
#' @param fetch_mode `"warn_and_cache"` (default) or `"abort"`. See
#'   `.fetch_external()` (internal): this is enrichment
#'   data, so failures do not stop the pipeline by default -- a compound
#'   that cannot be classified gets `pathway = NA` and a log entry, not a
#'   hard stop.
#'
#' @return The updated `proj`, with a `compounds_classified` entry in
#'   [patliRResults()] (columns `compound_id`, `pathway`, `superclass`,
#'   `class`, `isglycoside`, `source`, `fetch_date`), also written to
#'   `results/compounds_classified.csv`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- compounds_classify(proj) # needs internet
#' patliRResults(proj, "compounds_classified")
#' }
#'
#' @export
compounds_classify <- function(proj, compound_ids = NULL,
                                fetch_mode = c("warn_and_cache", "abort")) {
  stopifnot(is(proj, "PatliRProject"))
  fetch_mode <- match.arg(fetch_mode)

  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  if (nrow(cmp) == 0) {
    cli::cli_warn("No matching compounds in {.arg proj}; run {.fn prep_compounds} first. Nothing to do.")
    return(proj)
  }

  fetch_date <- as.character(Sys.Date())
  rows <- lapply(seq_len(nrow(cmp)), function(i) {
    id <- cmp$id[i]
    smiles_key <- cmp$canonical_smiles[i] %||% cmp$smiles[i]
    cls <- .fetch_external(
      fetch_fun = function() .npclassifier_lookup(smiles_key),
      cache_dir = cacheDir(proj), cache_key = paste0("npclassifier_", id), mode = fetch_mode
    )
    if (is.null(cls)) {
      return(data.frame(
        compound_id = id, pathway = NA_character_, superclass = NA_character_,
        class = NA_character_, isglycoside = NA, source = "npclassifier",
        fetch_date = fetch_date, stringsAsFactors = FALSE
      ))
    }
    data.frame(
      compound_id = id,
      pathway = cls$pathway %||% NA_character_,
      superclass = cls$superclass %||% NA_character_,
      class = cls$class %||% NA_character_,
      isglycoside = cls$isglycoside %||% NA,
      source = "npclassifier", fetch_date = fetch_date, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)

  failed <- is.na(result$pathway)
  if (any(failed)) {
    proj <- .log_append(
      proj, step = "compounds_classify", id = result$compound_id[failed],
      message = "NPClassifier lookup failed or returned no pathway; compound left unclassified"
    )
  }

  result <- .network_upsert(proj, "compounds_classified", result, "compound_id")
  patliRResults(proj, "compounds_classified") <- result
  .write_results_csv(proj, "compounds_classified", result)
  .write_log_csv(proj)
  proj
}

#' Look up a compound's NPClassifier pathway/superclass/class by SMILES
#'
#' @section Why this throttles and retries (2026-08-11):
#' On real data, `compounds_classify()` calls this once per compound in a
#' tight `lapply()` (tens to low hundreds of back-to-back requests). Tested
#' directly against the live endpoint: isolated requests -- including ones
#' with stereochemistry (`@`, `@@`) and percent-escaped characters -- always
#' return valid JSON, so this is not a SMILES-encoding bug. But real
#' natural-product runs (74/74 compounds, EFLO-S dataset) got a `200 OK`
#' whose body was HTML instead of the documented JSON, for essentially
#' every compound. That pattern -- works in isolation, fails under volume --
#' is the signature of a rate-limit or anti-bot response returned with a
#' success status rather than `429`/`503` (so httr2's default transient-
#' response detection, which only looks at status code, would never catch
#' it). NPClassifier/GNPS2 publish no documented rate limit, so there is no
#' official number to throttle to; one request per second is a conservative,
#' polite default. `req_throttle()`'s token bucket is shared across calls to
#' the same host within one R session, so this applies across the whole
#' `compounds_classify()` loop, not just within one call.
#' @return `list(pathway, superclass, class, isglycoside)`, the first
#'   (highest-confidence) entry of each field NPClassifier returns.
#' @keywords internal
.npclassifier_lookup <- function(smiles) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("the 'httr2' package is required to query NPClassifier.")
  }
  if (is.na(smiles) || !nzchar(smiles)) stop("need a SMILES to query NPClassifier")
  url <- paste0(
    "https://npclassifier.gnps2.org/classify?smiles=",
    utils::URLencode(smiles, reserved = TRUE)
  )

  req <- httr2::request(url)
  req <- httr2::req_user_agent(
    req, "patliR-Rpkg (network pharmacology research tool; https://github.com/)"
  )
  req <- httr2::req_throttle(req, capacity = 1, fill_time_s = 1, realm = "npclassifier.gnps2.org")
  req <- httr2::req_retry(
    req,
    max_tries = 4,
    is_transient = function(resp) {
      httr2::resp_status(resp) %in% c(429, 503) ||
        !identical(httr2::resp_content_type(resp), "application/json")
    },
    backoff = function(n) 2^n
  )

  resp <- httr2::req_perform(req)

  if (!identical(httr2::resp_content_type(resp), "application/json")) {
    body_preview <- tryCatch(
      substr(httr2::resp_body_string(resp), 1, 300),
      error = function(e) "<could not read response body>"
    )
    stop(
      "NPClassifier returned content-type '", httr2::resp_content_type(resp),
      "' instead of JSON, even after retrying with backoff. The endpoint ",
      "itself is known to return valid JSON for isolated requests (verified ",
      "2026-08-11), so this usually means it is rate-limiting or blocking ",
      "rapid sequential requests, not that this SMILES is malformed. ",
      "Response body preview: ", body_preview
    )
  }

  body <- httr2::resp_body_json(resp)
  list(
    pathway = if (length(body$pathway_results) > 0) body$pathway_results[[1]] else NA_character_,
    superclass = if (length(body$superclass_results) > 0) body$superclass_results[[1]] else NA_character_,
    class = if (length(body$class_results) > 0) body$class_results[[1]] else NA_character_,
    isglycoside = isTRUE(body$isglycoside)
  )
}
