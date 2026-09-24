#' @include AllGenerics.R internal.R
NULL

#' Generate and QC 2D structure depictions for every compound
#'
#' @description
#' Renders a 2D depiction (PNG) for every compound currently in
#' [compounds()] and records whether generation succeeded. This is a visual
#' QC step, additional to (and independent from) the structural validation
#' [prep_compounds()] already performs via `rcdk` -- it exists to catch
#' structures that parse correctly but depict in a degenerate or unexpected
#' way (disconnected fragments, wrong valence rendering, etc.), and to give
#' you a quick sanity-check image per compound before moving further down
#' the pipeline.
#'
#' @inheritParams compounds
#' @param engine `"rcdk"` (default): use `rcdk`'s own 2D depiction, no extra
#'   dependency. `"chemminer"` is not implemented and is rejected.
#' @param out_dir Directory to write the PNG files to. Defaults to
#'   `file.path(projectDir(proj), "structure2d")`.
#'
#' @return The updated `proj`, with a `structure2d_log` entry in
#'   [patliRResults()] (columns `id`, `generated_ok`, `engine_used`,
#'   `path`, `failure_reason`), also written to `results/structure2d_log.csv`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- prep_structure2d(proj)
#' patliRResults(proj, "structure2d_log")
#' }
#'
#' @export
prep_structure2d <- function(proj, engine = c("rcdk", "chemminer"), out_dir = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  engine <- match.arg(engine)
  if (engine == "chemminer") {
    cli::cli_abort("{.code engine = \"chemminer\"} is not implemented; use {.code engine = \"rcdk\"}.")
  }

  cmp <- compounds(proj)
  if (nrow(cmp) == 0) {
    cli::cli_warn("No compounds in {.arg proj} yet; run {.fn prep_compounds} first. Nothing to do.")
    return(proj)
  }

  if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "structure2d")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  results <- lapply(seq_len(nrow(cmp)), function(i) {
    id <- cmp$id[i]
    smi <- cmp$smiles[i]
    path <- file.path(out_dir, paste0(id, ".png"))
    failure <- tryCatch({
      .depict_2d(smi, path, engine = engine)
      NA_character_
    }, error = function(e) conditionMessage(e))
    ok <- is.na(failure)
    data.frame(
      id = id, generated_ok = ok, engine_used = engine,
      path = if (ok) path else NA_character_,
      failure_reason = failure,
      stringsAsFactors = FALSE
    )
  })
  log_df <- do.call(rbind, results)

  failed <- log_df[!log_df$generated_ok, "id"]
  if (length(failed) > 0) {
    proj <- .log_append(
      proj, step = "prep_structure2d", id = failed,
      message = paste0("2D depiction failed with engine '", engine, "': ",
                       log_df$failure_reason[!log_df$generated_ok])
    )
  }

  patliRResults(proj, "structure2d_log") <- log_df
  .write_results_csv(proj, "structure2d_log", log_df)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.depict_2d <- function(smiles, path, engine) {
  mol <- rcdk::parse.smiles(smiles)[[1]]
  if (is.null(mol)) stop("could not parse SMILES for depiction")
  if (engine == "rcdk") {
    img <- rcdk::view.image.2d(mol)
    grDevices::png(path, width = 300, height = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot(0:1, 0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
    graphics::rasterImage(img, 0, 0, 1, 1)
  } else {
    stop("ChemmineR depiction path not yet implemented")
  }
  invisible(path)
}
