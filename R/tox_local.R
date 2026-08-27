#' @include AllGenerics.R internal.R
NULL

#' Structural toxicity/interference alerts (PAINS, Brenk) -- local, no network
#'
#' @description
#' Flags each compound against published libraries of "structural alerts" --
#' substructures statistically associated with either assay interference
#' (PAINS) or generally undesirable screening-library behavior (Brenk) --
#' using SMARTS substructure matching via [rcdk::matches()]. This is the same
#' underlying mechanism (CDK's SMARTS matcher) that RDKit uses internally to
#' implement PAINS filtering in Python; no RDKit/`reticulate` dependency is
#' needed.
#'
#' @section This never produces a pass/fail verdict:
#' `tox_local()` has no `hard_cutoff` argument, unlike [adme_filter()], and
#' this is a deliberate design choice, not an oversight: a PAINS/Brenk flag
#' very often reflects real pharmacological activity (vasodilators,
#' vasoconstrictors, and other reactive-but-intentional chemotypes routinely
#' trip these filters) rather than a genuine liability. See [tox_report()]
#' for the one function in this family that summarizes results -- it always
#' carries a fixed disclaimer for exactly this reason. If you want to exclude
#' compounds, do so explicitly and manually from `patliRResults(proj,
#' "tox_local")`, with your own judgment about which alerts actually matter
#' for your extract/target combination.
#'
#' @section Coverage -- PAINS and Brenk, both the full published set:
#' `alert_sets = "pains"` uses the SMARTS list bundled at
#' `inst/extdata/pains_smarts.csv` -- the **complete 480-filter WEHI PAINS
#' set**, sourced verbatim from RDKit's `Data/Pains/wehi_pains.csv`
#' (BSD-3-Clause, which itself redistributes the substructure filters from
#' Baell, J.B. & Holloway, G.A. (2010), "New Substructure Filters for
#' Removal of Pan Assay Interference Compounds (PAINS) from Screening
#' Libraries...", *J. Med. Chem.* 53(7), 2719-2740,
#' \doi{10.1021/jm901137j}). An earlier version of this file only bundled a
#' 121-filter subset (frequency >= 3) because the automated fetch used to
#' source it kept truncating near the end of the file; the full set was
#' downloaded directly (no truncation) and confirmed complete.
#'
#' `alert_sets = "brenk"` uses `inst/extdata/brenk_smarts.csv` -- the
#' **complete 105-alert set** from Brenk, R. et al. (2008), "Lessons
#' Learnt from Assembling Screening Libraries for Drug Discovery for
#' Neglected Diseases", *ChemMedChem* 3, 435-444,
#' \doi{10.1002/cmdc.200700139}. Cross-validated two ways before bundling
#' (see `patliR_manual.md` for the full provenance/reproduction script):
#' count and alert names confirmed
#' against RDKit's own compiled `FilterCatalogs.BRENK` (105 entries), exact
#' SMARTS text taken from PatWalters/rd_filters' `alert_collection.csv`
#' (MIT; its `"Dundee"`-labelled rows -- Brenk's group was the University
#' of Dundee Drug Discovery Unit), every one of the 105 patterns
#' independently re-parsed with `rdkit.Chem.MolFromSmarts()` to confirm it
#' is well-formed. Not transcribed from memory.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound currently in [compounds()].
#' @param alert_sets Character vector, any of `"pains"`, `"brenk"`. Default
#'   `c("pains", "brenk")`.
#'
#' @return The updated `proj`, with a `tox_local` entry in [patliRResults()]
#'   (columns `compound_id`, `alert_set`, `alert_name`, `smarts`, `matched`),
#'   also written to `results/tox_local.csv`. Only alerts that actually
#'   matched at least one compound's structure trigger a log entry; the
#'   full result table (matched and non-matched) is always in
#'   `patliRResults()`.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- tox_local(proj, alert_sets = c("pains", "brenk"))
#' patliRResults(proj, "tox_local")
#' }
#'
#' @export
tox_local <- function(proj, compound_ids = NULL, alert_sets = c("pains", "brenk")) {
  stopifnot(is(proj, "PatliRProject"))
  alert_sets <- match.arg(alert_sets, c("pains", "brenk"), several.ok = TRUE)

  cmp <- .compounds_or_all(proj, compound_ids)
  if (nrow(cmp) == 0) {
    cli::cli_abort("No compounds found in {.arg proj}; run {.fn prep_compounds} first.")
  }

  alerts <- do.call(rbind, lapply(alert_sets, function(set) {
    tbl <- .load_alert_smarts(set)
    data.frame(alert_set = set, alert_name = tbl$name, smarts = tbl$smarts, stringsAsFactors = FALSE)
  }))

  mols <- .parse_smiles_safe(cmp$smiles)

  rows <- lapply(seq_len(nrow(cmp)), function(i) {
    id <- cmp$id[i]
    mol <- mols[[i]]
    if (is.null(mol)) {
      return(data.frame(
        compound_id = id, alert_set = alerts$alert_set, alert_name = alerts$alert_name,
        smarts = alerts$smarts, matched = NA, stringsAsFactors = FALSE
      ))
    }
    matched <- vapply(alerts$smarts, function(sm) {
      res <- tryCatch(rcdk::matches(sm, mol), error = function(e) NA)
      isTRUE(res[1])
    }, logical(1))
    data.frame(
      compound_id = id, alert_set = alerts$alert_set, alert_name = alerts$alert_name,
      smarts = alerts$smarts, matched = matched, stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)

  unparsed <- cmp$id[vapply(mols, is.null, logical(1))]
  if (length(unparsed) > 0) {
    proj <- .log_append(
      proj, step = "tox_local", id = unparsed,
      message = paste0("tox_local: SMILES could not be parsed by rcdk; all ", paste(alert_sets, collapse = "/"),
                        " alerts recorded as NA for this compound")
    )
  }

  result <- .network_upsert(proj, "tox_local", result, "compound_id")
  patliRResults(proj, "tox_local") <- result
  .write_results_csv(proj, "tox_local", result)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.compounds_or_all <- function(proj, compound_ids) {
  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  cmp
}

#' @keywords internal
.load_alert_smarts <- function(set) {
  switch(set,
    pains = .load_pains_smarts(),
    brenk = .load_brenk_smarts(),
    cli::cli_abort("Unknown {.arg alert_sets} value: {.val {set}}.")
  )
}

#' Load the bundled PAINS SMARTS subset
#'
#' @section `name` is not a unique key:
#' The upstream WEHI PAINS list encodes each filter's original identity as a
#' `regId` string like `"amino_acridine_A(1)"` -- base name plus a frequency
#' count in parens. Our loader splits that into separate `name`/`frequency`
#' columns (see the transcription script at the end of `TESTING_GUIDE.Rmd`),
#' which drops the guarantee of uniqueness: the real bundled data has the
#' base name `"amino_acridine_A"` twice, for two different SMARTS
#' (frequencies 1 and 46). `smarts` is always unique (it is the actual
#' structural identity); `(name, frequency)` together are also always
#' unique. Never assume `name` alone identifies a row -- see
#' `test-tox.R`'s test of this function for the confirmed invariants.
#' @return `data.frame(smarts, name, frequency)`.
#' @keywords internal
.load_pains_smarts <- function() {
  utils::read.csv(system.file("extdata", "pains_smarts.csv", package = "patliR"), stringsAsFactors = FALSE)
}

#' Load the bundled Brenk SMARTS set (105 alerts, Brenk et al. 2008)
#'
#' See the "Coverage" section of [tox_local()] for provenance. Errors
#' clearly (rather than fabricating SMARTS from memory) if
#' `inst/extdata/brenk_smarts.csv` is missing from this install -- should
#' not happen in a normal package install, only if `inst/extdata/` was
#' pruned somehow.
#' @return `data.frame(smarts, name)`.
#' @keywords internal
.load_brenk_smarts <- function() {
  path <- system.file("extdata", "brenk_smarts.csv", package = "patliR")
  if (!nzchar(path) || !file.exists(path)) {
    cli::cli_abort(c(
      "{.val inst/extdata/brenk_smarts.csv} not found in this {.pkg patliR} install.",
      "i" = "Re-install the package, or see {.file TESTING_GUIDE.Rmd} to regenerate it.",
      "i" = "Use {.code tox_local(proj, alert_sets = \"pains\")} in the meantime."
    ))
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}
