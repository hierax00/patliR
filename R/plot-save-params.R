#' Shared parameters for the `plot_*` family's file-saving behavior
#'
#' Every `plot_*` function returns its plot object for further tweaking and,
#' by default, also writes it to the project's `plots/` directory. These
#' parameters control that save step and are documented once here, inherited
#' via `@inheritParams plot_save_params` elsewhere.
#'
#' @param save Logical, default `TRUE`. Whether to also write the plot to
#'   disk (in addition to returning it).
#' @param out_dir Directory to write the file(s) to when `save = TRUE`.
#'   `NULL` (default) means `file.path(projectDir(proj), "plots")`.
#' @param width,height Numeric, plot size in inches, passed to
#'   [ggplot2::ggsave()] (or the function's native device).
#' @param dpi Numeric resolution in dots per inch for raster output, passed
#'   to [ggplot2::ggsave()].
#'
#' @name plot_save_params
#' @keywords internal
NULL
