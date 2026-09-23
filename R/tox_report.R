#' @include AllGenerics.R internal.R
NULL

#' Summarize toxicity results -- never a pass/fail verdict
#'
#' @description
#' Builds a per-compound summary combining [tox_local()]'s structural alerts
#' (PAINS/Brenk), [tox_safetyome()]'s target-level systemic-risk flags, and
#' [tox_import()]'s imported properties (hERG, DILI, Ames, ...), when any is
#' available. This function **never** collapses its output into a single
#' pass/fail flag, and always attaches a fixed disclaimer -- see the
#' Description of [tox_local()] for why: many compounds trip these filters
#' because of genuine pharmacological activity, not because they are
#' unsafe.
#'
#' @inheritParams compounds
#'
#' @return A `list` with two elements:
#'   \itemize{
#'     \item `summary`: a `data.frame` (one row per compound with any
#'       [tox_local()], [tox_safetyome()], and/or [tox_import()] results),
#'       columns `compound_id`, `n_pains_alerts`, `pains_alert_names`
#'       (comma-separated), `n_brenk_alerts`, `brenk_alert_names`
#'       (comma-separated), `n_safetyome_hits`,
#'       `safetyome_organ_systems` (comma-separated, deduplicated),
#'       `n_imported_properties`, `imported_properties` (comma-separated
#'       `property=value` pairs).
#'     \item `note`: the fixed disclaimer string, always identical
#'       regardless of the data (`"A toxicity/structural-alert flag can
#'       reflect real pharmacological activity, not just risk -- this is
#'       not a pass/fail verdict. Review each flagged compound
#'       individually."`).
#'   }
#'   `proj` itself is not modified (this function returns a `list`, not a
#'   `proj`), but the `summary` table is written to `results/tox_report.csv`
#'   inside [projectDir()] as a side effect -- the same plain-CSV
#'   source-of-truth every other family produces, so the per-compound
#'   toxicity/PAINS view survives the R session and is picked up by
#'   [patliR_load()] (as `patliRResults(proj, "tox_report")`) and by
#'   [patliR_export_llm()]. Nothing is written when there are no
#'   [tox_local()]/[tox_safetyome()]/[tox_import()] results to summarize.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- tox_local(proj, alert_sets = c("pains", "brenk"))
#' report <- tox_report(proj)
#' report$summary
#' report$note
#' # also on disk, as the durable per-compound view:
#' read.csv(file.path(projectDir(proj), "results", "tox_report.csv"))
#' }
#'
#' @export
tox_report <- function(proj) {
  stopifnot(is(proj, "PatliRProject"))

  note <- paste(
    "A toxicity/structural-alert flag can reflect real pharmacological",
    "activity, not just risk -- this is not a pass/fail verdict. Review",
    "each flagged compound individually."
  )

  local_res <- patliRResults(proj, "tox_local")
  safetyome_res <- patliRResults(proj, "tox_safetyome")
  imported_res <- patliRResults(proj, "tox_imported")

  if (all(vapply(list(local_res, safetyome_res, imported_res),
                 function(x) is.null(x) || nrow(x) == 0, logical(1)))) {
    cli::cli_warn(paste(
      "No {.val tox_local}, {.val tox_safetyome}, or {.val tox_imported}",
      "results found; run {.fn tox_local}, {.fn tox_safetyome}, and/or",
      "{.fn tox_import} first."
    ))
    return(list(
      summary = data.frame(
        compound_id = character(),
        n_pains_alerts = integer(), pains_alert_names = character(),
        n_brenk_alerts = integer(), brenk_alert_names = character(),
        n_safetyome_hits = integer(), safetyome_organ_systems = character(),
        n_imported_properties = integer(), imported_properties = character(),
        stringsAsFactors = FALSE
      ),
      note = note
    ))
  }

  ids <- unique(c(
    if (!is.null(local_res)) local_res$compound_id,
    if (!is.null(safetyome_res)) safetyome_res$compound_id,
    if (!is.null(imported_res)) imported_res$compound_id
  ))

  rows <- lapply(ids, function(id) {
    .alert_set_hits <- function(set) {
      if (is.null(local_res)) return(data.frame(alert_name = character()))
      local_res[local_res$compound_id == id & local_res$alert_set == set &
                   !is.na(local_res$matched) & local_res$matched, ]
    }
    pains_hits <- .alert_set_hits("pains")
    brenk_hits <- .alert_set_hits("brenk")

    safetyome_hits <- if (!is.null(safetyome_res)) {
      safetyome_res[safetyome_res$compound_id == id &
                       !is.na(safetyome_res$in_core_panel) & safetyome_res$in_core_panel, ]
    } else {
      data.frame(organ_system = character())
    }

    imp_rows <- if (!is.null(imported_res)) imported_res[imported_res$compound_id == id, ] else imported_res[0, ]
    imp_pairs <- if (!is.null(imported_res) && nrow(imp_rows) > 0) {
      paste(imp_rows$property, imp_rows$value, sep = "=", collapse = ", ")
    } else {
      NA_character_
    }

    data.frame(
      compound_id = id,
      n_pains_alerts = nrow(pains_hits),
      pains_alert_names = if (nrow(pains_hits) > 0) paste(pains_hits$alert_name, collapse = ", ") else NA_character_,
      n_brenk_alerts = nrow(brenk_hits),
      brenk_alert_names = if (nrow(brenk_hits) > 0) paste(brenk_hits$alert_name, collapse = ", ") else NA_character_,
      n_safetyome_hits = nrow(safetyome_hits),
      safetyome_organ_systems = if (nrow(safetyome_hits) > 0) {
        paste(unique(safetyome_hits$organ_system), collapse = ", ")
      } else {
        NA_character_
      },
      n_imported_properties = if (!is.null(imported_res)) nrow(imp_rows) else 0L,
      imported_properties = imp_pairs,
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, rows)
  .write_results_csv(proj, "tox_report", summary)

  list(summary = summary, note = note)
}
