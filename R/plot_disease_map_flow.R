#' @include AllGenerics.R internal.R plot-helpers.R plot_disease_map.R
NULL

## Same alluvial family as plot_condition_flow()/plot_condition_disease_flow() --
## see plot_condition_flow.R's header comment for the geom_flow() vs.
## geom_alluvium() rationale. This one generalises the single, pooled
## "Disease areas reached by the extract's targets" chart of plot_disease_map()
## to a caller-given ordered sequence of conditions: same 14-area taxonomy
## (disease_map_classes()), same target-disease association data
## (.disease_map_associations()), but read as a river across conditions
## instead of one static snapshot.

#' River across an ordered sequence of conditions: which disease areas the
#' extract's target set reaches, stage by stage
#'
#' @description
#' One flow per predicted target, moving left to right through `conditions`
#' **in the order given** (e.g. a phenological sequence), coloured by its
#' single most-associated disease area from [disease_map_classes()] -- the
#' same 14-area taxonomy and colour scheme as [plot_disease_map()], but
#' varying by condition instead of one pooled snapshot. A ribbon narrows,
#' widens or breaks as targets belonging to that area drop out of or appear
#' in [network_build()]'s edges from one condition to the next. Reads as
#' "which therapeutic areas does each stage's target set lean towards",
#' complementing [plot_condition_flow()]'s "how does the extract's chemistry
#' change" and [plot_condition_disease_flow()]'s "which (curated) disease
#' does each condition line up with best".
#'
#' @details
#' A target can be associated with diseases from more than one area (its
#' target-disease rows are pooled across `sources`, exactly as in
#' [plot_disease_map()]). To draw one ribbon per target, each is assigned a
#' single **primary** area: whichever area holds the most of its diseases
#' (ties broken by the higher association score) -- the same one-label-per-
#' unit simplification [plot_condition_flow()] makes for a compound's
#' chemical class. `attr(., "table")` also carries this assignment
#' (`target_id`, `disease_class`) so it can be checked or overridden upstream.
#'
#' @inheritParams network_build
#' @inheritParams plot_disease_map
#' @param conditions Character vector, **the river's left-to-right order** --
#'   not resolved or sorted for you, because the whole point is a caller-
#'   chosen sequence. At least two conditions, all already built by
#'   [network_build()], no duplicates.
#' @inheritParams plot_save_params
#'
#' @return A `ggplot` (`ggalluvial` has no `ggiraph` equivalent). If
#'   `save = TRUE` (default), also writes a PNG and logs it to
#'   `patliRResults(proj, "disease_map_flow_plot_log")`. `attr(., "table")`
#'   holds the exact `(condition, target_id, disease_class)` rows drawn.
#'
#' @seealso [plot_disease_map()], [disease_map_classes()],
#'   [plot_condition_flow()], [plot_condition_disease_flow()]
#' @export
plot_disease_map_flow <- function(proj, conditions, disease_classes = NULL,
                                  sources = c("profile", "disease_genes", "targets_disease"),
                                  min_score = 0.1, save = TRUE, out_dir = NULL,
                                  width = NULL, height = 7, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  .plot_require(extra = "ggalluvial")
  if (!is.character(conditions) || length(conditions) < 2 || anyNA(conditions)) {
    cli::cli_abort("{.arg conditions} must be a character vector of at least 2 condition names, in the order to draw them.")
  }
  if (anyDuplicated(conditions)) cli::cli_abort("{.arg conditions} must not contain duplicate condition names.")
  if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score) || min_score < 0 || min_score > 1) {
    cli::cli_abort("{.arg min_score} must be a single number in [0, 1].")
  }
  if (!is.character(sources) || length(sources) == 0 ||
        length(setdiff(sources, c("profile", "disease_genes", "targets_disease"))) > 0) {
    cli::cli_abort("{.arg sources} must be a non-empty subset of {.val {c('profile', 'disease_genes', 'targets_disease')}}.")
  }
  classes <- .disease_map_check_classes(disease_classes)

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c("No {.val network_edges} entry in {.arg proj}.", "i" = "Run {.fn network_build} first."))
  }
  missing_cond <- setdiff(conditions, unique(edges_all$condition))
  if (length(missing_cond) > 0) {
    cli::cli_abort("Condition(s) {.val {missing_cond}} not built by {.fn network_build}; available: {.val {unique(edges_all$condition)}}.")
  }

  ct_all <- edges_all[edges_all$condition %in% conditions, , drop = FALSE]
  assoc <- .disease_map_associations(proj, unique(ct_all$uniprot_id), sources = sources, min_score = min_score)
  if (nrow(assoc) == 0) {
    cli::cli_abort(c(
      "No target-disease association for the targets of condition(s) {.val {conditions}}.",
      "i" = "Run {.fn targets_disease_profile} (or {.fn disease_genes_fetch}) first, or lower {.arg min_score}."
    ))
  }
  assoc$disease_class <- .disease_map_classify(assoc$disease_name, classes)

  ## one primary area per target -- see @details
  n_dz <- stats::aggregate(disease_id ~ target_id + disease_class, assoc,
                           function(x) length(unique(x)))
  names(n_dz)[3] <- "n_diseases"
  ## na.action = na.pass: formula aggregate() defaults to na.omit, which
  ## drops every row with an NA response *before* grouping -- a
  ## (target_id, disease_class) whose association_score is NA for every row
  ## (e.g. a GO-sourced association, which carries no score at all) would
  ## then vanish entirely instead of reaching the `all(is.na(x))` branch
  ## below, and if EVERY row in `assoc` is unscored this aggregate() errors
  ## outright ("no rows to aggregate") rather than returning an empty result.
  best_sc <- stats::aggregate(association_score ~ target_id + disease_class, assoc,
                              function(x) if (all(is.na(x))) -Inf else max(x, na.rm = TRUE),
                              na.action = stats::na.pass)
  agg <- merge(n_dz, best_sc, by = c("target_id", "disease_class"), all.x = TRUE)
  agg$association_score[is.na(agg$association_score)] <- -Inf
  agg <- agg[order(agg$target_id, -agg$n_diseases, -agg$association_score), , drop = FALSE]
  primary <- agg[!duplicated(agg$target_id), c("target_id", "disease_class")]
  target_class <- stats::setNames(primary$disease_class, primary$target_id)

  long <- do.call(rbind, lapply(conditions, function(cd) {
    ## unique(): the same guard plot_condition_flow() applies to compound_id
    present <- unique(ct_all$uniprot_id[ct_all$condition == cd])
    present <- present[present %in% names(target_class)]
    if (length(present) == 0) return(NULL)
    data.frame(condition = cd, target_id = present,
              disease_class = unname(target_class[present]), stringsAsFactors = FALSE)
  }))
  if (is.null(long) || nrow(long) == 0) {
    cli::cli_abort("None of {.arg conditions}'s targets has a disease-area association in scope.")
  }
  long$condition <- factor(long$condition, levels = conditions)
  long$y <- 1

  classes_present <- intersect(c(setdiff(names(classes), "Other"), "Other"), unique(long$disease_class))
  colors <- .disease_map_palette(names(classes))[classes_present]

  if (is.null(width)) width <- max(8, 1.4 * length(conditions) + 3)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$condition, stratum = .data$disease_class,
                                          alluvium = .data$target_id, y = .data$y)) +
    ggalluvial::geom_flow(ggplot2::aes(fill = .data$disease_class), alpha = 0.7, width = 1 / 6) +
    ggalluvial::geom_stratum(width = 1 / 6, fill = "grey92", colour = "grey40") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(.plot_truncate(as.character(stratum), 20))),
                       size = 2.5) +
    ggplot2::scale_fill_manual(values = colors, name = "Disease area", breaks = classes_present) +
    ggplot2::labs(
      title = "Which disease areas does the extract's target set reach, condition by condition?",
      subtitle = .plot_wrap(paste0(
        "One ribbon per predicted target, present (per network_build()) at each condition in the order given; ",
        "coloured by its single most-associated disease area (disease_map_classes(); most of its diseases fall ",
        "there, ties broken by association score). A ribbon narrows/breaks where that area's targets drop out ",
        "of the condition and widens/restarts where new ones appear."
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
    name = "disease_map_flow_plot_log",
    filename = paste0("disease_map_flow_", paste(conditions, collapse = "-"), ".png"),
    log_row = data.frame(conditions = paste(conditions, collapse = ">"),
                         n_targets = length(unique(long$target_id)), n_classes = length(classes_present),
                         path = NA_character_, stringsAsFactors = FALSE),
    key_cols = c("conditions"),
    engine = NULL, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
  attr(result, "table") <- long[, c("condition", "target_id", "disease_class")]
  result
}
