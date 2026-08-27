#' @include AllGenerics.R internal.R prep_structure2d.R
NULL

## plot_structure2d() -- purely an assembly step. prep_structure2d()
## already renders one PNG per compound (rcdk/ChemmineR 2D depiction);
## this just arranges the already-rendered PNGs into one grid figure with
## compound_id/name captions, reusing patchwork (already added for
## plot_upset()) instead of hand-rolling base graphics::layout(). The one
## new thing needed is a way to read an existing PNG file back in as a
## raster to place inside a ggplot panel -- the `png` package, about as
## minimal as that gets (base R has no PNG reader of its own).

#' Grid of 2D structure depictions from `prep_structure2d()`
#'
#' @description
#' Arranges the PNG files [prep_structure2d()] already rendered per
#' compound into one grid figure, each panel captioned with the
#' compound's name (falling back to `compound_id`). Purely an assembly
#' step -- does not regenerate any depiction itself; run
#' [prep_structure2d()] first.
#'
#' @inheritParams compounds
#' @inheritParams plot_save_params
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound with a successful (`generated_ok`)
#'   `structure2d_log` entry.
#' @param ncol Integer, default `4`. Number of grid columns.
#'
#' @return A `patchwork` object (see [plot_upset()] for why `patchwork`,
#'   not a single `ggplot`). If `save = TRUE` (default), also writes a PNG
#'   and logs it to `patliRResults(proj, "structure2d_grid_log")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- prep_structure2d(proj)
#' plot_structure2d(proj, save = FALSE)
#' }
#'
#' @export
plot_structure2d <- function(proj, compound_ids = NULL, ncol = 4,
                              save = TRUE, out_dir = NULL, width = 10, height = 8, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(ncol), length(ncol) == 1, ncol >= 1)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_structure2d}.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_structure2d} needs the {.pkg patchwork} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"patchwork\")}."
    ))
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_structure2d} needs the {.pkg png} package (CRAN), not installed -- it is how the already-rendered depictions get read back in as images.",
      "i" = "{.code install.packages(\"png\")}."
    ))
  }

  log_df <- patliRResults(proj, "structure2d_log")
  if (is.null(log_df) || nrow(log_df) == 0) {
    cli::cli_abort(c("No {.val structure2d_log} entry in {.arg proj}.", "i" = "Run {.fn prep_structure2d} first."))
  }
  dat <- log_df[log_df$generated_ok & !is.na(log_df$path), , drop = FALSE]
  if (!is.null(compound_ids)) dat <- dat[dat$id %in% compound_ids, , drop = FALSE]
  if (nrow(dat) == 0) {
    cli::cli_abort("No successfully-generated 2D structures to plot for this selection.")
  }

  cmp <- compounds(proj)
  name_lookup <- stats::setNames(cmp$name, cmp$id)
  dat$label <- ifelse(is.na(name_lookup[dat$id]) | name_lookup[dat$id] == "", dat$id, name_lookup[dat$id])

  panels <- lapply(seq_len(nrow(dat)), function(i) {
    img <- tryCatch(png::readPNG(dat$path[i]), error = function(e) NULL)
    p <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = dat$label[i]) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 8, hjust = 0.5))
    if (!is.null(img)) {
      p <- p + ggplot2::annotation_raster(img, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
    } else {
      p <- p + ggplot2::annotate("text", x = 0, y = 0, label = "(image unreadable)", size = 2.5, colour = "grey50")
    }
    p
  })

  p <- patchwork::wrap_plots(panels, ncol = ncol)

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, "structure2d_grid.png")
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, limitsize = FALSE)
    grid_log <- data.frame(path = path, n_compounds = nrow(dat), stringsAsFactors = FALSE)
    patliRResults(proj, "structure2d_grid_log") <- grid_log
    .write_results_csv(proj, "structure2d_grid_log", grid_log)
  }

  if (save) attr(p, "proj") <- proj
  p
}
