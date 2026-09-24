#' @include AllGenerics.R internal.R adme_filter.R
NULL

## Same hand-rolled ggplot2 + patchwork UpSet layout as plot_upset(), but
## over "which compounds pass which combination of adme_filter() rules"
## instead of "which targets are shared across which conditions". Kept
## separate from plot_upset() (different axes) rather than adding a
## `domain` switch. UpSet not Venn: 5 rules make a proportional Venn
## unreadable.

#' UpSet-style plot: drug-likeness/lead-likeness rule overlaps across compounds
#'
#' @description
#' Which compounds pass which *combination* of [adme_filter()]'s rules, and
#' how many -- e.g. "42 compounds pass Ro5+Veber+Ghose but fail Oprea
#' (not lead-like)". Needs [adme_filter()] to have been run first (any
#' `rules` selection with 2+ rules).
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param top_n Integer, default `15`. Only the `top_n` largest
#'   intersections are shown.
#'
#' @return A `patchwork` object (bar chart on top of the dot-matrix,
#'   `print()`-able and `ggsave()`-able like a single `ggplot`). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "adme_upset_log")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' proj <- adme_filter(proj, rules = c("ro5", "veber", "ghose", "egan", "oprea"))
#' plot_adme_upset(proj, save = FALSE)
#' }
#'
#' @export
plot_adme_upset <- function(proj, top_n = 15,
                             save = TRUE, out_dir = NULL, width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(top_n), length(top_n) == 1, top_n >= 1)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn plot_adme_upset}.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_adme_upset} needs the {.pkg patchwork} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"patchwork\")} -- pure CRAN, ggplot2 + grid only."
    ))
  }

  long <- patliRResults(proj, "adme_filtered")
  if (is.null(long) || nrow(long) == 0) {
    cli::cli_abort("No {.val adme_filtered} results found; run {.fn adme_filter} first.")
  }
  ## route_* columns are per-administration-route compatibility flags, not
  ## a drug-likeness/lead-likeness rule -- excluded from this plot on
  ## purpose (they answer a different question, "which route" not "how
  ## drug-like").
  long <- long[!grepl("^route_", long$rule), , drop = FALSE]
  rules <- sort(unique(long$rule))
  if (length(rules) < 2) {
    cli::cli_abort("{.fn plot_adme_upset} needs at least 2 rules in {.val adme_filtered} to compare; got {.val {rules}}.")
  }

  passing <- long[long$pass %in% TRUE, , drop = FALSE]
  combo <- if (nrow(passing) > 0) {
    stats::aggregate(rule ~ compound_id, passing,
                     function(x) paste(sort(unique(x)), collapse = "|"))
  } else {
    data.frame(compound_id = character(0), rule = character(0), stringsAsFactors = FALSE)
  }
  ## Compounds that pass none of the rules still form a real (empty-set)
  ## intersection -- keep them rather than silently dropping.
  all_ids <- unique(long$compound_id)
  none <- setdiff(all_ids, combo$compound_id)
  if (length(none) > 0) {
    combo <- rbind(combo, data.frame(compound_id = none, rule = "(none)", stringsAsFactors = FALSE))
  }

  sizes <- as.data.frame(table(combo$rule), stringsAsFactors = FALSE)
  names(sizes) <- c("combo", "n_compounds")
  sizes <- sizes[order(-sizes$n_compounds), , drop = FALSE]
  sizes <- utils::head(sizes, top_n)
  sizes$combo <- factor(sizes$combo, levels = sizes$combo)

  dot <- do.call(rbind, lapply(as.character(sizes$combo), function(cb) {
    members <- if (cb == "(none)") character(0) else strsplit(cb, "\\|")[[1]]
    data.frame(combo = cb, rule = rules, in_set = rules %in% members, stringsAsFactors = FALSE)
  }))
  dot$combo <- factor(dot$combo, levels = levels(sizes$combo))

  p_bar <- ggplot2::ggplot(sizes, ggplot2::aes(x = .data$combo, y = .data$n_compounds)) +
    ggplot2::geom_col(fill = "#27ae60") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_compounds), vjust = -0.4, size = 3) +
    ggplot2::labs(title = "Drug-likeness / lead-likeness rule intersections", x = NULL, y = "Compounds") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )

  p_dot <- ggplot2::ggplot(dot, ggplot2::aes(x = .data$combo, y = .data$rule)) +
    ggplot2::geom_line(data = dot[dot$in_set, , drop = FALSE], ggplot2::aes(group = .data$combo), colour = "grey30", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$in_set), size = 3) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey85"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())

  p <- patchwork::wrap_plots(p_bar, p_dot, ncol = 1, heights = c(2, 1))

  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, "adme_upset.png")
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
    log_row <- data.frame(path = path, n_rules = length(rules), n_intersections = nrow(sizes), stringsAsFactors = FALSE)
    patliRResults(proj, "adme_upset_log") <- log_row
    .write_results_csv(proj, "adme_upset_log", log_row)
  }

  if (save) attr(p, "proj") <- proj
  p
}
