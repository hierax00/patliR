#' patliR: Reproducible Network Pharmacology Analysis for Natural Product
#' Extracts
#'
#' @description
#' `patliR` takes a set of compounds (a pure isolate, a curated compound
#' list, or a raw GC-MS abundance matrix from one or more extracts) through a
#' full network pharmacology pipeline: compound identity and structure
#' validation, physicochemical/ADME profiling, structural toxicity alerting,
#' target prediction or import, compound-target-pathway network analysis,
#' database-bias auditing, candidate ranking, and a fully-logged report.
#'
#' @section Design principles:
#' Every exported function that transforms the pipeline state takes a
#' [`PatliRProject-class`] object and returns a new, updated one -- nothing
#' is modified in place. Every step also writes its result to a plain CSV
#' file inside the project directory, so the analysis can always be resumed
#' from disk with [patliR_load()], even in a brand new R session, without
#' depending on any R-specific binary format.
#'
#' @section Getting started:
#' See `vignette("patliR-intro", package = "patliR")` for a full walkthrough,
#' or start with [patliR_project()] and [prep_compounds()].
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom methods new validObject setClass setGeneric setMethod setValidity is
#' @importFrom stats quantile median mad
#' @importFrom utils packageVersion write.csv read.csv
#' @importFrom tools R_user_dir
#' @importFrom rlang .data
## usethis namespace: end
NULL
