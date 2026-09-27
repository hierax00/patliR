#' @include AllGenerics.R internal.R plot-helpers.R
NULL

## Uses ggalluvial (Suggests), same rationale as the other alluvial plots in this
## package. Unlike plot_bowtie()/plot_gpcr_flow()/plot_condition_disease_flow()
## (2-4 fixed axes), this one has an axis per *condition in a caller-given
## sequence*, and uses geom_flow() on long ("lode") data instead of
## geom_alluvium() on wide data, because a compound is not guaranteed to be
## present at every step -- geom_flow() draws that as the ribbon vanishing and
## (if it returns later) reappearing, which is exactly the reading wanted here.

#' River across an ordered sequence of conditions: how the extract's
#' composition changes from one stage to the next
#'
#' @description
#' One flow per compound, moving left to right through `conditions` **in the
#' order given** -- e.g. a phenological sequence (vegetative -> flowering ->
#' fruiting -> post-flowering) -- coloured by a chemical grouping
#' (`group_by`), with the ribbon narrowing, widening, or breaking as
#' compounds drop out of or appear in the presence/absence calls
#' ([prep_binarize()]) from one condition to the next. Reads as "how does
#' the extract's composition flow between conditions", rather than [plot_condition_compare()]'s
#' "which single condition is best" or [plot_condition_disease_flow()]'s
#' "which disease does each condition line up with".
#'
#' @inheritParams network_build
#' @param conditions Character vector, **the river's left-to-right order** --
#'   not resolved or sorted for you, because the whole point is a caller-
#'   chosen sequence (a phenological order, a treatment order, ...). At
#'   least two conditions, all present in [binarizedMatrix()].
#' @param group_by `"pathway"` (default), `"superclass"` or `"class"` --
#'   which [compounds_classify()] column colours the flows. `"pathway"` is
#'   usually the most legible (few, broad categories); `"class"` is the
#'   finest and typically needs a low `top_n_groups`.
#' @param top_n_groups Integer, default `8`. Groups beyond the `top_n_groups`
#'   most frequent (by total compound-condition presences across
#'   `conditions`) are pooled into `"Other"` -- with many small chemical
#'   classes the river gets illegible otherwise.
#' @inheritParams plot_save_params
#'
#' @return A `ggplot` (`ggalluvial` has no `ggiraph` equivalent). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "condition_flow_plot_log")`. `attr(., "table")`
#'   holds the exact `(condition, compound_id, group)` rows drawn.
#'
#' @seealso [plot_condition_compare()], [plot_condition_disease_flow()], [compounds_classify()]
#' @export
plot_condition_flow <- function(proj, conditions, group_by = c("pathway", "superclass", "class"),
                                top_n_groups = 8, save = TRUE, out_dir = NULL,
                                width = NULL, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  group_by <- match.arg(group_by)
  .plot_require(extra = "ggalluvial")
  if (!is.character(conditions) || length(conditions) < 2 || anyNA(conditions)) {
    cli::cli_abort("{.arg conditions} must be a character vector of at least 2 condition names, in the order to draw them.")
  }
  .pathway_check_count(top_n_groups, "top_n_groups", min = 1)

  bin <- binarizedMatrix(proj)
  if (is.null(bin)) cli::cli_abort(c("No binarized presence matrix in {.arg proj}.", "i" = "Run {.fn prep_binarize} first."))
  missing_cond <- setdiff(conditions, setdiff(names(bin), "compound_id"))
  if (length(missing_cond) > 0) {
    cli::cli_abort("Condition(s) {.val {missing_cond}} not found in {.fn binarizedMatrix}; available: {.val {setdiff(names(bin), 'compound_id')}}.")
  }

  cls <- patliRResults(proj, "compounds_classified")
  if (is.null(cls) || nrow(cls) == 0) {
    cli::cli_abort(c("No {.val compounds_classified} entry in {.arg proj}.", "i" = "Run {.fn compounds_classify} first."))
  }
  grp <- stats::setNames(cls[[group_by]], cls$compound_id)

  long <- do.call(rbind, lapply(conditions, function(cd) {
    present <- bin$compound_id[!is.na(bin[[cd]]) & bin[[cd]] == 1]
    if (length(present) == 0) return(NULL)
    g <- grp[present]
    data.frame(condition = cd, compound_id = present,
              group = ifelse(is.na(g), "Unclassified", g), stringsAsFactors = FALSE)
  }))
  if (is.null(long) || nrow(long) == 0) {
    cli::cli_abort("No compound is present (per {.fn binarizedMatrix}) in any of {.arg conditions}.")
  }

  totals <- sort(table(long$group), decreasing = TRUE)
  keep <- names(utils::head(totals, top_n_groups))
  long$group[!long$group %in% keep] <- "Other"
  long$condition <- factor(long$condition, levels = conditions)
  long$y <- 1

  groups_present <- if ("Other" %in% long$group) c(setdiff(sort(unique(long$group)), "Other"), "Other") else sort(unique(long$group))
  colors <- .plot_contrast_palette(length(groups_present))
  names(colors) <- groups_present
  if ("Other" %in% names(colors)) colors["Other"] <- "grey70"

  if (is.null(width)) width <- max(6, 1.3 * length(conditions) + 3)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$condition, stratum = .data$group,
                                          alluvium = .data$compound_id, y = .data$y)) +
    ggalluvial::geom_flow(ggplot2::aes(fill = .data$group), alpha = 0.7, width = 1 / 6) +
    ggalluvial::geom_stratum(width = 1 / 6, fill = "grey92", colour = "grey40") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(count)), size = 2.7) +
    ggplot2::scale_fill_manual(values = colors, name = tools::toTitleCase(group_by)) +
    ggplot2::labs(
      title = "How the extract's composition flows across conditions",
      subtitle = .plot_wrap(paste0(
        "One ribbon per compound, present (per prep_binarize()) at each condition in the order given; ",
        "a ribbon narrows/breaks where compounds drop out and widens/restarts where new ones appear. ",
        "Colour = ", group_by, "."
      ), .plot_wrap_width(width, 8)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey40"),
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  result <- .plot_finish(
    proj, p,
    name = "condition_flow_plot_log",
    filename = paste0("condition_flow_", paste(conditions, collapse = "-"), "_", group_by, ".png"),
    log_row = data.frame(conditions = paste(conditions, collapse = ">"), group_by = group_by,
                         n_compounds = length(unique(long$compound_id)), path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("conditions", "group_by"),
    engine = NULL, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- long[, c("condition", "compound_id", "group")]
  result
}
