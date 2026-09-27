#' @include AllGenerics.R internal.R
NULL

#' Mark, and optionally remove, compounds by ADME rule
#'
#' @description
#' By default `adme_filter()` only *marks* pass/fail per rule -- it never
#' removes a compound silently. Set `hard_cutoff = TRUE` to actually drop
#' failing compounds from [compounds()] going forward, which always shows
#' you what would be removed (and why) first.
#'
#' @inheritParams compounds
#' @param rules Character vector, any of `"ro5"`, `"veber"`, `"ghose"`,
#'   `"egan"`, `"oprea"`, `"route"`. `"oprea"` is the set of property ranges
#'   (Oprea 2000: HBD 0-2, HBA 2-9, rotatable bonds 2-8, rings 1-4, inclusive)
#'   holding most of the drug-like compounds that paper surveyed; see [adme_local()]
#'   for every rule's exact criteria and citation. `"route"` expands to
#'   every `route_*` column [adme_local()] computed (whichever routes you
#'   asked for there).
#' @param source Currently only `"local"` (i.e. [adme_local()]'s output) is
#'   supported; `"imported"` will raise an informative error until a future
#'   version adds it.
#' @param hard_cutoff Logical (default `FALSE`). If `TRUE`, compounds that
#'   fail *any* selected rule are removed from [compounds()] (see `ask`);
#'   if `FALSE`, nothing is removed -- the return value only adds the
#'   `adme_filtered` marks.
#' @param ask Logical (default `interactive()`). Only used when
#'   `hard_cutoff = TRUE`. If `TRUE`, prompts before removing anything
#'   (`y`/`n`/`v` to view the list first). If `FALSE` (e.g. in a script),
#'   the cut is applied without prompting, but full detail is always
#'   written to [projectLog()] regardless.
#' @param ro5_max_violations `NULL` (default) uses [adme_local()]'s own
#'   `ro5_pass` as stored (zero violations -- all four Lipinski criteria).
#'   An integer 0-4 recomputes pass/fail here instead, from
#'   `adme_local()`'s `ro5_violations` column, without rerunning
#'   [adme_local()] -- e.g. `ro5_max_violations = 1` for the tolerance the
#'   original paper itself uses. Only applies when `"ro5"` is in `rules`;
#'   ignored otherwise. Needs `ro5_violations` in the `adme_local` results
#'   (added in this version -- an older run needs `adme_local()` rerun).
#'
#' @return The updated `proj`, with an `adme_filtered` entry in
#'   [patliRResults()] in long format (columns `compound_id`, `rule`,
#'   `pass`), also written to `results/adme_filtered.csv`. If
#'   `hard_cutoff = TRUE` and removal was confirmed (or `ask = FALSE`),
#'   [compounds()] is also updated and re-written to `01_compounds.csv`.
#'
#' @examples
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' proj <- adme_filter(proj, rules = c("ro5", "veber"))
#' patliRResults(proj, "adme_filtered")
#'
#' @export
adme_filter <- function(proj, rules = c("ro5", "veber", "ghose", "egan", "oprea", "route"),
                         source = c("local", "imported"),
                         hard_cutoff = FALSE, ask = interactive(),
                         ro5_max_violations = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  rules <- match.arg(rules, c("ro5", "veber", "ghose", "egan", "oprea", "route"), several.ok = TRUE)
  source <- match.arg(source)
  if (!is.null(ro5_max_violations)) {
    stopifnot(is.numeric(ro5_max_violations), length(ro5_max_violations) == 1,
             !is.na(ro5_max_violations), ro5_max_violations >= 0, ro5_max_violations <= 4)
  }

  if (source == "imported") {
    cli::cli_abort(c(
      "{.arg source = \"imported\"} is not implemented yet in this version of patliR.",
      "i" = "Use {.code source = \"local\"} (the default) with {.fn adme_local}'s output for now."
    ))
  }

  adme <- patliRResults(proj, "adme_local")
  if (is.null(adme) || nrow(adme) == 0) {
    cli::cli_abort("No {.val adme_local} results found; run {.fn adme_local} first.")
  }

  rule_cols <- unlist(lapply(rules, function(r) {
    if (r == "route") {
      grep("^route_", names(adme), value = TRUE)
    } else {
      paste0(r, "_pass")
    }
  }))
  rule_cols <- intersect(rule_cols, names(adme))
  if (length(rule_cols) == 0) {
    cli::cli_abort("None of the requested {.arg rules} have a matching column in {.val adme_local}; did you request a route that {.fn adme_local} was not run with?")
  }

  ## ro5_max_violations recomputes ro5_pass here (from the violation count adme_local()
  ## stored) rather than reusing the strict, zero-violation column as-is
  ro5_override <- NULL
  if (!is.null(ro5_max_violations) && "ro5_pass" %in% rule_cols) {
    if (!"ro5_violations" %in% names(adme)) {
      cli::cli_abort(c(
        "{.arg ro5_max_violations} needs an {.val ro5_violations} column in {.val adme_local}, not found.",
        "i" = "Rerun {.fn adme_local} (this version adds it) before using {.arg ro5_max_violations}."
      ))
    }
    ro5_override <- adme$ro5_violations <= ro5_max_violations
  }

  long <- do.call(rbind, lapply(rule_cols, function(col) {
    pass <- if (col == "ro5_pass" && !is.null(ro5_override)) ro5_override else adme[[col]]
    data.frame(compound_id = adme$compound_id, rule = col, pass = pass, stringsAsFactors = FALSE)
  }))

  patliRResults(proj, "adme_filtered") <- long
  .write_results_csv(proj, "adme_filtered", long)

  if (!hard_cutoff) {
    return(proj)
  }

  fails_any <- unique(long$compound_id[!is.na(long$pass) & !long$pass])
  unknown <- vapply(rule_cols, function(col) sum(is.na(adme[[col]])), integer(1))
  if (any(unknown > 0)) {
    ## a compound whose descriptors could not be computed is kept, not passed: say so
    cli::cli_inform(c("!" = "Compounds with no result (descriptors missing) are kept, not evaluated: {paste0(names(unknown)[unknown > 0], ' = ', unknown[unknown > 0], collapse = ', ')}."))
  }
  if (length(fails_any) == 0) {
    if (all(unknown == nrow(adme))) {
      cli::cli_inform("None of the selected rules could be evaluated (all results are NA); nothing to cut.")
    } else {
      cli::cli_inform("No compounds fail any of the selected rules; nothing to cut.")
    }
    return(proj)
  }

  cmp <- compounds(proj)
  n_total <- nrow(cmp)
  by_rule <- vapply(rule_cols, function(col) {
    pass <- if (col == "ro5_pass" && !is.null(ro5_override)) ro5_override else adme[[col]]
    sum(!pass, na.rm = TRUE)
  }, integer(1))
  summary_msg <- paste0(rule_cols, ": ", by_rule, collapse = ", ")

  proceed <- TRUE
  if (isTRUE(ask)) {
    proceed <- .ask_hard_cutoff(length(fails_any), n_total, summary_msg, fails_any, cmp)
  }

  proj <- .log_append(
    proj, step = "adme_filter", id = fails_any,
    message = paste0("hard_cutoff removal (rules: ", paste(rule_cols, collapse = ", "), "), ",
                      if (proceed) "applied" else "cancelled by user")
  )

  if (proceed) {
    compounds(proj) <- cmp[!cmp$id %in% fails_any, , drop = FALSE]
    .write_step_csv(proj, "01_compounds.csv", compounds(proj))
  }

  .write_log_csv(proj)
  proj
}

#' @keywords internal
.ask_hard_cutoff <- function(n_remove, n_total, summary_msg, fails_any, cmp) {
  cli::cli_inform("Failures by rule: {summary_msg}")
  repeat {
    answer <- tolower(trimws(readline(sprintf(
      "Remove %d of %d compounds that fail at least one selected rule? [y/n/v] ", n_remove, n_total
    ))))
    if (answer == "v") {
      print(cmp[cmp$id %in% fails_any, c("id", "name", "pubchem_id")])
      next
    }
    return(answer == "y")
  }
}
