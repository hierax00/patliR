#' @include AllGenerics.R internal.R network_build.R network_bowtie.R plot-helpers.R
NULL

## Arc diagram (Krzywinski et al. 2009 idiom) for target-target adjacency
## within one condition -- network-design-spec.md Sec 3.9. Reuses
## network_bowtie()'s STRING "actions" machinery
## (.network_stringdb(), .network_stringdb_actions_graph(),
## .network_actions_flag_true()) for the target-target edges rather than
## re-fetching STRING or re-deriving the directed/action-flag parsing --
## see network_bowtie.R's own docs for what the actions file is and why
## v11.0 is the default actions_version. network_bowtie.R is used
## unmodified here (it is also being read, unmodified, by the bow-tie
## family) -- nothing in this file changes its behaviour.
##
## Nodes sit on a single horizontal line, not a circle: the spec's own
## recommendation is "do the arc version; revisit circlize only if a user
## asks for the circular one." Order is by degree (within the
## target-target actions graph, restricted to this condition's targets),
## ascending left-to-right.
##
## Arcs are a hand-rolled quadratic Bezier via geom_path, not geom_curve
## (unlike plot_network_degeneracy()'s existing arc overlay): on a 2-D
## force-directed layout every degeneracy link spans roughly the same
## screen distance, so a constant geom_curve() curvature already looks
## consistent; here node spacing is 1-D and spans range from 1 (adjacent
## nodes) to n-1 (the two ends of the line), so a constant curvature would
## make a short-hop arc barely visible next to a towering long-hop one.
## The hand-rolled version makes "arc height scales with node span"
## explicit (.network_target_chord_bezier(), internal) instead of relying
## on geom_curve()'s implicit curvature-to-chord-length relationship.

#' Arc diagram of a condition's target-target STRING "actions" adjacency
#'
#' @description
#' network-design-spec.md Sec 3.9: nodes are a condition's targets, placed
#' along a single horizontal line (ordered by degree, ascending
#' left-to-right) rather than around a circle -- the spec's own
#' recommendation is "do the arc version, [...] revisit `circlize` only if
#' a user asks for the circular one." One arc is drawn above the line per
#' pair of targets connected in the STRING directed "actions" network
#' [network_bowtie()] already downloads and parses (same file, same
#' `actions_version` default `"11.0"`, same species default `9606`;
#' direction is dropped here -- this plot asks only "do these two targets
#' interact", not "which one acts on the other", which is what
#' [network_bowtie()]'s bow-tie decomposition is for). A subset of the
#' highest-degree targets (`top_n_labels`) is labeled via
#' [.plot_label_nodes()]; every other plotted target gets an unlabeled tick,
#' to keep a condition with dozens of targets legible.
#'
#' @section Arc colour -- when a per-edge STRING score is, and is not, available:
#' [network_bowtie()]'s cached actions graph
#' (`.network_stringdb_actions_graph()`, internal) keeps only
#' `item_id_a`/`item_id_b` when it parses the STRING actions flat file --
#' the file's own `score` column (0-999, per-interaction confidence) is
#' dropped, since [network_bowtie()] itself has no use for it, and
#' `network_bowtie.R` is used unmodified here. This function makes a
#' best-effort second, independent pass over the *same*, already
#' downloaded/cached raw actions flat file (same cache path, same
#' `is_directional`/`a_is_acting` filter, reusing
#' [network_bowtie()]'s own `.network_actions_flag_true()`) to recover
#' `score` when that raw file is still present on disk. When it is not
#' (e.g. only the parsed-graph `.rds` cache survived a cache cleanup, or in
#' every automated test, which mocks the actions graph directly and never
#' downloads anything), `actions_score_threshold` has no effect and every
#' arc is drawn in one constant colour -- exactly the fallback
#' network-design-spec.md Sec 3.9 itself anticipates ("otherwise a
#' constant colour").
#' A present but unreadable or invalid raw score file raises an error;
#' confidence filtering is never silently bypassed for a corrupt file.
#'
#' @inheritParams network_build
#' @inheritParams plot_save_params
#' @param actions_score_threshold Single number in `[0, 999]`, default
#'   `400` -- the same STRING confidence default [network_proximity()] /
#'   [network_synergy()] use for the *undirected* interactome (the actions
#'   file's own `score` column uses the same 0-999 scale). Only has an
#'   effect when a per-edge `score` could be recovered -- see the section
#'   above; a `cli_inform` names the fallback when it could not.
#' @param top_n_labels Integer, default `15`. The `top_n_labels`
#'   highest-degree targets (ties broken by `uniprot_id`, for a
#'   deterministic order) get a text label; every other plotted target
#'   gets an unlabeled tick mark instead.
#'   Labels are fanned out along the baseline (short leader lines back to
#'   their ticks) so neighbouring high-degree targets do not overprint, and
#'   names longer than 20 characters are shortened (full name in the
#'   `ggiraph` tooltip).
#' @param top_n_targets `NULL` (default, every target in scope) or an
#'   integer >= 2: draw only the `top_n_targets` targets with the most
#'   actions interactions among the condition's targets (ties by
#'   `uniprot_id`), and the arcs among them. A condition with a few hundred
#'   targets otherwise reduces the baseline to a smear of ticks (a
#'   `cli_inform` suggests this above 150 targets).
#' @param engine `"static"` (default) or `"ggiraph"`, save/out_dir/width/
#'   height/dpi -- same as [plot_network_layers()].
#'
#' @return A `ggplot` object (`engine = "static"`) or a `girafe` htmlwidget
#'   (`engine = "ggiraph"`). If `save = TRUE` (default), also writes a PNG
#'   and logs it to `patliRResults(proj, "target_chord_plot_log")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' abundance <- read.csv(
#'   system.file("extdata", "input_abundance_matrix.csv", package = "patliR"),
#'   check.names = FALSE
#' )
#' proj <- prep_binarize(proj, abundance)
#' proj <- targets_import_batch(
#'   proj,
#'   system.file("extdata", "import_targets", package = "patliR"),
#'   platform = "superpred"
#' )
#' proj <- network_build(proj)
#' # needs STRINGdb + internet on first call (downloads/caches the STRING
#' # "actions" flat file, same as network_bowtie())
#' plot_target_chord(proj, condition = "FLO-ET", save = FALSE)
#' }
#'
#' @references
#' Krzywinski et al. (2009), *Genome Res* 19(9):1639-1645,
#' \doi{10.1101/gr.092759.109}.
#'
#' @export
plot_target_chord <- function(proj, condition = NULL, actions_score_threshold = 400, top_n_labels = 15,
                               top_n_targets = NULL,
                               engine = c("static", "ggiraph"), save = TRUE, out_dir = NULL,
                               width = 8, height = 6, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(actions_score_threshold), length(actions_score_threshold) == 1,
            !is.na(actions_score_threshold), actions_score_threshold >= 0, actions_score_threshold <= 999)
  stopifnot(is.numeric(top_n_labels), length(top_n_labels) == 1, !is.na(top_n_labels), top_n_labels >= 0)
  stopifnot(is.null(top_n_targets) ||
              (is.numeric(top_n_targets) && length(top_n_targets) == 1 && !is.na(top_n_targets) && top_n_targets >= 2))
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn plot_target_chord} needs {.pkg STRINGdb} (for UniProt -> STRING_id mapping, same as {.fn network_bowtie}), not installed.",
      "i" = "Install it with {.code BiocManager::install(\"STRINGdb\")}."
    ))
  }
  engine <- match.arg(engine)
  engine <- .plot_require(engine)
  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  ## Same fixed defaults network_bowtie() itself uses -- not exposed here
  ## to keep this (explicitly lowest-priority) plot's signature to what
  ## network-design-spec.md Sec 3.9 asks for.
  species <- 9606
  actions_version <- "11.0"

  edges_all <- patliRResults(proj, "network_edges")
  targets <- sort(unique(edges_all$uniprot_id[edges_all$condition %in% conditions]))
  if (length(targets) < 2) {
    cli::cli_abort(c(
      "Condition(s) {.val {conditions}} have {length(targets)} target(s) in scope.",
      "i" = "{.fn plot_target_chord} needs at least 2 targets to draw any arc; widen {.arg condition} or check {.fn network_build}'s output."
    ))
  }

  string_db <- .network_stringdb(proj, species, actions_version, score_threshold = 0)
  map_df <- string_db$map(data.frame(uniprot_id = targets, stringsAsFactors = FALSE),
                           "uniprot_id", removeUnmappedRows = FALSE, quiet = TRUE)
  uni_to_string <- stats::setNames(map_df$STRING_id, map_df$uniprot_id)

  g <- .network_stringdb_actions_graph(proj, species, actions_version)$graph
  edge_df <- .target_chord_edges(uni_to_string, g)
  if (nrow(edge_df) == 0) {
    cli::cli_abort(c(
      "No STRING {.val actions} interaction found among the {length(targets)} target(s) of condition(s) {.val {conditions}}.",
      "i" = "This is a common, legitimate outcome -- the directed actions network ({.fn network_bowtie}'s own data source) is much sparser than STRING's full interactome; it is not an error in this plot."
    ))
  }

  scores <- .target_chord_actions_scores(proj, species, actions_version)
  has_score <- !is.null(scores)
  if (has_score) {
    n_before <- nrow(edge_df)
    edge_df <- merge(edge_df, scores, by = c("string_a", "string_b"), all.x = TRUE)
    edge_df <- edge_df[!is.na(edge_df$score) & edge_df$score >= actions_score_threshold, , drop = FALSE]
    if (nrow(edge_df) == 0) {
      cli::cli_abort(c(
        "Every one of the {n_before} STRING actions interaction(s) among these targets scored below {.arg actions_score_threshold} = {actions_score_threshold}.",
        "i" = "Lower {.arg actions_score_threshold}, or treat this as a legitimate 'no confident interaction at this threshold' result -- 400 is the package-wide STRING confidence default."
      ))
    }
  } else {
    edge_df$score <- NA_real_
    cli::cli_inform(c(
      "i" = "{.fn plot_target_chord}: STRING actions {.field score} is not available (the raw actions flat file is not on disk -- only the parsed graph is) -- drawing every arc in one colour and skipping {.arg actions_score_threshold} filtering."
    ))
  }

  n_targets_all <- length(targets)
  if (!is.null(top_n_targets)) {
    kept <- .target_chord_top_targets(targets, edge_df, top_n_targets)
    targets <- kept$targets
    edge_df <- kept$edge_df
  } else if (length(targets) > 150) {
    cli::cli_inform(c(
      "i" = "{.fn plot_target_chord}: {length(targets)} targets on one line -- most ticks will be indistinguishable; pass {.arg top_n_targets} (e.g. {.code top_n_targets = 60}) to keep only the most connected ones."
    ))
  }

  labels <- .plot_label_nodes(proj, conditions, targets, "target")
  nodes <- .target_chord_nodes(targets, edge_df, top_n_labels, labels)

  idx <- stats::setNames(nodes$x, nodes$uniprot_id)
  edge_df$x0 <- unname(idx[edge_df$uniprot_a])
  edge_df$x1 <- unname(idx[edge_df$uniprot_b])
  edge_df$edge_id <- seq_len(nrow(edge_df))
  arcs <- .network_target_chord_bezier(edge_df)

  p <- .network_target_chord_ggplot(nodes, arcs, has_score, engine, scope_label, fig_width = width, fig_height = height)

  ## top_n_targets joined the log key later -- older rows drew every target
  top_n_used <- if (is.null(top_n_targets)) NA_real_ else as.numeric(top_n_targets)
  proj <- .plot_log_backfill(proj, "target_chord_plot_log", "top_n_targets", NA_real_)
  .plot_finish(
    proj, p,
    name = "target_chord_plot_log",
    filename = paste0("target_chord_", scope_label, if (!is.na(top_n_used)) paste0("_top", top_n_used), ".png"),
    log_row = data.frame(
      condition = scope_label, top_n_targets = top_n_used, path = NA_character_,
      n_targets = length(targets), n_targets_in_scope = n_targets_all,
      n_edges = nrow(edge_df), actions_score_threshold = actions_score_threshold,
      score_available = has_score, stringsAsFactors = FALSE
    ),
    key_cols = c("condition", "top_n_targets"),
    engine = engine, save = save, out_dir = out_dir,
    width = width, height = height, dpi = dpi
  )
}

#' Undirected target-target edges from the STRING actions graph, restricted
#' to a set of UniProt targets
#'
#' @description
#' `uni_to_string` maps this condition's targets to STRING ids (possibly
#' `NA`, for a target absent from STRING altogether); `g` is
#' [network_bowtie()]'s cached directed actions graph for the whole
#' species. Direction is dropped and reciprocal/duplicate edges collapsed:
#' this plot asks only "do these two targets interact", not "which one acts
#' on the other". `string_a`/`string_b` are ordered independently of
#' `uniprot_a`/`uniprot_b` (by STRING id, not by whichever accession a
#' STRING id happens to have been mapped from) so
#' `.target_chord_actions_scores()`'s join key matches regardless of that
#' accession.
#'
#' @return `data.frame(uniprot_a, uniprot_b, string_a, string_b)`, zero rows
#'   (with the right columns) if fewer than 2 targets mapped onto `g`, or
#'   `g` has no edge between any two of them.
#' @keywords internal
.target_chord_edges <- function(uni_to_string, g) {
  mapped <- uni_to_string[!is.na(uni_to_string)]
  g_names <- igraph::V(g)$name
  in_graph <- mapped[mapped %in% g_names]
  if (length(unique(in_graph)) < 2) {
    return(.empty_target_chord_edges())
  }

  sub <- igraph::induced_subgraph(g, vids = unique(in_graph))
  ed <- igraph::as_data_frame(sub, what = "edges")
  if (nrow(ed) == 0) {
    return(.empty_target_chord_edges())
  }

  ## First uniprot_id wins if >1 accession maps to the same STRING id (a
  ## rare isoform collision) -- negligible for a plot whose whole point is
  ## "does an interaction exist", not exact accession bookkeeping.
  string_to_uni <- stats::setNames(names(in_graph), in_graph)
  ed$uniprot_a <- unname(string_to_uni[ed$from])
  ed$uniprot_b <- unname(string_to_uni[ed$to])
  ed <- ed[!is.na(ed$uniprot_a) & !is.na(ed$uniprot_b) & ed$uniprot_a != ed$uniprot_b, , drop = FALSE]
  if (nrow(ed) == 0) {
    return(.empty_target_chord_edges())
  }

  a_first <- ed$uniprot_a < ed$uniprot_b
  out <- data.frame(
    uniprot_a = ifelse(a_first, ed$uniprot_a, ed$uniprot_b),
    uniprot_b = ifelse(a_first, ed$uniprot_b, ed$uniprot_a),
    string_a = ifelse(ed$from < ed$to, ed$from, ed$to),
    string_b = ifelse(ed$from < ed$to, ed$to, ed$from),
    stringsAsFactors = FALSE
  )
  key <- paste(out$uniprot_a, out$uniprot_b, sep = "\r")
  out[!duplicated(key), , drop = FALSE]
}

#' @keywords internal
.empty_target_chord_edges <- function() {
  data.frame(uniprot_a = character(0), uniprot_b = character(0),
             string_a = character(0), string_b = character(0), stringsAsFactors = FALSE)
}

#' Best-effort per-edge STRING actions `score`, read directly from the
#' cached raw flat file
#'
#' @description
#' [network_bowtie()]'s `.network_stringdb_actions_graph()` (internal)
#' downloads and caches this exact file but drops its `score` column when
#' building the parsed graph. Modifying that shared function was out of
#' scope for this plot (see the "Arc colour" section of
#' `plot_target_chord()`'s own docs); this makes a second, independent pass
#' over the same cached file, when it is still on disk, reusing
#' [network_bowtie()]'s own `.network_actions_flag_true()` for the
#' direction filter rather than re-deriving it.
#'
#' @return `data.frame(string_a, string_b, score)` (undirected, deduplicated
#'   by keeping the strongest of any multiple action "modes" between the
#'   same pair), or `NULL` only if the raw file is not on disk. Unreadable
#'   or invalid files raise an error so confidence filtering is not bypassed.
#' @keywords internal
.target_chord_actions_scores <- function(proj, species, version) {
  raw_gz <- file.path(cacheDir(proj), "stringdb", paste0(species, ".protein.actions.v", version, ".txt.gz"))
  if (!file.exists(raw_gz)) return(NULL)

  actions <- tryCatch({
    con <- gzfile(raw_gz, "rt")
    on.exit(close(con), add = TRUE)
    utils::read.delim(con, stringsAsFactors = FALSE)
  }, error = function(e) {
    cli::cli_abort("Cannot read STRING actions scores from {.file {raw_gz}}; cannot apply {.arg actions_score_threshold}.", parent = e)
  }, warning = function(w) {
    cli::cli_abort("Cannot reliably read STRING actions scores from {.file {raw_gz}}; cannot apply {.arg actions_score_threshold}.", parent = w)
  })
  required <- c("item_id_a", "item_id_b", "is_directional", "a_is_acting", "score")
  if (!all(required %in% names(actions))) {
    cli::cli_abort("Invalid STRING actions score schema in {.file {raw_gz}}: missing {toString(setdiff(required, names(actions)))}; cannot apply {.arg actions_score_threshold}.")
  }

  keep <- .network_actions_flag_true(actions$is_directional) & .network_actions_flag_true(actions$a_is_acting)
  sc <- actions[keep, c("item_id_a", "item_id_b", "score")]
  sc$score <- suppressWarnings(as.numeric(sc$score))
  if (any(!is.finite(sc$score) | sc$score < 0 | sc$score > 999) ||
      anyNA(sc$item_id_a) || anyNA(sc$item_id_b) ||
      any(!nzchar(sc$item_id_a) | !nzchar(sc$item_id_b))) {
    cli::cli_abort("Invalid STRING actions scores or target IDs in {.file {raw_gz}}; cannot apply {.arg actions_score_threshold}.")
  }
  if (nrow(sc) == 0) {
    return(data.frame(string_a = character(0), string_b = character(0), score = numeric(0)))
  }

  a_first <- sc$item_id_a < sc$item_id_b
  sc$string_a <- ifelse(a_first, sc$item_id_a, sc$item_id_b)
  sc$string_b <- ifelse(a_first, sc$item_id_b, sc$item_id_a)
  ## A pair can have >1 action "mode" (activation, inhibition, ...), each
  ## with its own score -- keep the strongest (max) as this plot's one
  ## colour value per arc.
  stats::aggregate(score ~ string_a + string_b, data = sc, FUN = max)
}

#' Node placement + label selection for `plot_target_chord()`'s arc diagram
#'
#' @description
#' Every one of `targets` gets a node, placed by ascending degree (within
#' `edge_df`, i.e. only counting edges among the plotted targets), ties
#' broken alphabetically by `uniprot_id` -- deterministic given the same
#' input, no unseeded randomness anywhere in this plot. The `top_n_labels`
#' highest-degree targets (same tie-break) get their pre-computed `labels`
#' entry; the rest get `NA` (drawn as an unlabeled tick by the caller).
#'
#' @param targets Character vector of UniProt ids, one node each.
#' @param edge_df `data.frame` with (at least) `uniprot_a`, `uniprot_b`
#'   columns -- e.g. `.target_chord_edges()`'s output.
#' @param top_n_labels Integer >= 0.
#' @param labels Character vector, same length/order as `targets` -- e.g.
#'   [.plot_label_nodes()]'s output.
#' @return `data.frame(uniprot_id, x, degree, label)`, one row per target,
#'   ordered by `x` (ascending degree).
#' @keywords internal
.target_chord_nodes <- function(targets, edge_df, top_n_labels, labels) {
  n <- length(targets)
  deg <- stats::setNames(rep(0L, n), targets)
  if (nrow(edge_df) > 0) {
    tab <- table(c(edge_df$uniprot_a, edge_df$uniprot_b))
    deg[names(tab)] <- as.integer(tab)
  }

  ord <- order(deg[targets], targets)
  ordered <- targets[ord]
  nodes <- data.frame(
    uniprot_id = ordered, x = seq_len(n), degree = unname(deg[ordered]),
    stringsAsFactors = FALSE
  )

  label_of <- stats::setNames(labels, targets)
  rank_desc <- order(-nodes$degree, nodes$uniprot_id)
  keep_n <- max(0L, min(top_n_labels, n))
  shown <- if (keep_n > 0) rank_desc[seq_len(keep_n)] else integer(0)
  nodes$label <- NA_character_
  if (length(shown) > 0) {
    nodes$label[shown] <- unname(label_of[nodes$uniprot_id[shown]])
  }
  nodes
}

#' Keep only the `top_n` most connected targets for `plot_target_chord()`
#'
#' @description
#' Degree within `edge_df` (edges among the plotted targets), ties broken
#' by `uniprot_id`; edges with an end outside the kept set are dropped.
#' @return `list(targets, edge_df)`.
#' @keywords internal
.target_chord_top_targets <- function(targets, edge_df, top_n) {
  deg <- stats::setNames(rep(0L, length(targets)), targets)
  if (nrow(edge_df) > 0) {
    tab <- table(c(edge_df$uniprot_a, edge_df$uniprot_b))
    tab <- tab[names(tab) %in% targets]
    deg[names(tab)] <- as.integer(tab)
  }
  keep <- utils::head(targets[order(-deg, targets)], top_n)
  list(
    targets = sort(keep),
    edge_df = edge_df[edge_df$uniprot_a %in% keep & edge_df$uniprot_b %in% keep, , drop = FALSE]
  )
}

#' Spread label anchors along the baseline so no two are closer than
#' `min_gap`
#'
#' @description
#' `plot_target_chord()` orders nodes by ascending degree, so the labeled
#' (highest-degree) targets are always the rightmost, adjacent ticks --
#' drawn at their own x their labels pile into one smear. Each label is
#' moved the least needed to keep `min_gap` between neighbours (a forward
#' pass pushing right, then -- if that overshoots `upper` -- a backward pass
#' pulling left); order is preserved, so leader lines from node to label
#' never cross. If `min_gap` cannot fit at all, labels are spaced evenly
#' over `[lower, upper]`.
#' @param x Numeric node positions.
#' @return Numeric label positions, same order as `x`.
#' @keywords internal
.target_chord_label_x <- function(x, min_gap, lower, upper) {
  n <- length(x)
  if (n == 0) return(numeric(0))
  o <- order(x)
  y <- x[o]
  if ((n - 1) * min_gap >= upper - lower) {
    y <- seq(lower, upper, length.out = n)
  } else {
    y[1] <- max(y[1], lower)
    for (i in seq_len(n - 1) + 1L) y[i] <- max(y[i], y[i - 1] + min_gap)
    if (y[n] > upper) {
      y[n] <- upper
      for (i in rev(seq_len(n - 1))) y[i] <- min(y[i], y[i + 1] - min_gap)
    }
  }
  out <- numeric(n)
  out[o] <- y
  out
}

#' Quadratic-Bezier arc points for `plot_target_chord()`'s arc diagram
#'
#' @description
#' One arc per row of `edge_df` (needs numeric `x0`, `x1` and integer
#' `edge_id`, plus any extra columns to carry through unchanged to every
#' point of that edge's arc, e.g. `score`). Height is `height_scale *
#' abs(x1 - x0)` -- an explicit, deterministic function of node span,
#' rather than relying on `geom_curve()`'s implicit curvature-to-chord-
#' length relationship (the idiom [plot_network_degeneracy()]'s arc overlay
#' uses, on a 2-D force-directed layout where every edge spans roughly the
#' same screen distance; here spans range from 1, adjacent nodes, to n-1,
#' the two ends of the line, so a constant `geom_curve()` curvature would
#' draw a barely visible short-hop arc right next to a towering long-hop
#' one).
#'
#' @param edge_df Data frame with numeric `x0`, `x1` and integer `edge_id`.
#' @param n_points Points per arc (default `40`).
#' @param height_scale Arc height per unit of horizontal span (default
#'   `0.35`).
#' @return One row per (edge, point): every column of `edge_df` (constant
#'   within an edge, repeated `n_points` times) plus `x`, `y`.
#' @keywords internal
.network_target_chord_bezier <- function(edge_df, n_points = 40, height_scale = 0.35) {
  if (nrow(edge_df) == 0) {
    return(cbind(edge_df[FALSE, , drop = FALSE], x = numeric(0), y = numeric(0)))
  }
  t <- seq(0, 1, length.out = n_points)
  rows <- lapply(seq_len(nrow(edge_df)), function(i) {
    x0 <- edge_df$x0[i]
    x1 <- edge_df$x1[i]
    mid_x <- (x0 + x1) / 2
    h <- height_scale * abs(x1 - x0)
    x <- (1 - t)^2 * x0 + 2 * (1 - t) * t * mid_x + t^2 * x1
    y <- 2 * (1 - t) * t * h  # baseline y = 0 at both ends (t = 0, t = 1)
    cbind(edge_df[rep(i, n_points), , drop = FALSE], x = x, y = y)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
.network_target_chord_ggplot <- function(nodes, arcs, has_score, engine, title_suffix, fig_width = 8, fig_height = 6) {
  interactive <- engine == "ggiraph" && requireNamespace("ggiraph", quietly = TRUE)
  nodes$tooltip <- sprintf("%s\ndegree=%d", ifelse(is.na(nodes$label), nodes$uniprot_id, nodes$label), nodes$degree)

  p <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4)

  if (nrow(arcs) > 0) {
    p <- if (has_score) {
      p + ggplot2::geom_path(
        data = arcs, ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id, colour = .data$score),
        linewidth = 0.6, alpha = 0.75
      ) +
        ggplot2::scale_colour_gradient(low = "#aed6f1", high = "#1b4f72", name = "STRING actions\nscore", na.value = "grey70")
    } else {
      p + ggplot2::geom_path(
        data = arcs, ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id),
        colour = "#2980b9", linewidth = 0.6, alpha = 0.75
      )
    }
  }

  unlabeled <- nodes[is.na(nodes$label), , drop = FALSE]
  labeled <- nodes[!is.na(nodes$label), , drop = FALSE]
  if (nrow(unlabeled) > 0) {
    p <- if (interactive) {
      p + ggiraph::geom_point_interactive(
        data = unlabeled,
        ggplot2::aes(x = .data$x, y = 0, tooltip = .data$tooltip, data_id = .data$uniprot_id),
        shape = 3, size = 1.2, colour = "grey50"
      )
    } else {
      p + ggplot2::geom_point(data = unlabeled, ggplot2::aes(x = .data$x, y = 0), shape = 3, size = 1.2, colour = "grey50")
    }
  }

  ## Vertical room for the angled labels under the baseline, in data units:
  ## the arcs span y in [0, y_top]; a label of k characters at 60 degrees
  ## hangs ~k * 0.055 * sin(60) inches below it.
  n <- nrow(nodes)
  y_top <- if (nrow(arcs) > 0) max(arcs$y) else 1
  panel_w <- max(fig_width - 2, 2)  # minus legend + margins
  panel_h <- max(fig_height - 1.3, 1.5) # minus title + subtitle
  label_h_in <- if (nrow(labeled) > 0) max(nchar(.plot_truncate(labeled$label, 20))) * 0.055 * 0.87 + 0.25 else 0.1
  label_h_in <- min(label_h_in, 0.45 * panel_h)
  y_bottom <- -y_top * label_h_in / (panel_h - label_h_in)
  drop <- 0.08 * abs(y_bottom) / label_h_in # leader-line length, ~0.08 in

  if (nrow(labeled) > 0) {
    p <- if (interactive) {
      p + ggiraph::geom_point_interactive(
        data = labeled,
        ggplot2::aes(x = .data$x, y = 0, tooltip = .data$tooltip, data_id = .data$uniprot_id),
        shape = 16, size = 1.8, colour = "#c0392b"
      )
    } else {
      p + ggplot2::geom_point(data = labeled, ggplot2::aes(x = .data$x, y = 0), shape = 16, size = 1.8, colour = "#c0392b")
    }
    ## Labeled targets are the highest-degree ones, i.e. the rightmost
    ## adjacent ticks: spread their label anchors ~0.17 in apart
    ## (.target_chord_label_x()) and join each to its tick with a short
    ## leader line, instead of printing every label at its own tick.
    min_gap <- 0.17 * max(n - 1, 1) / panel_w
    labeled$label_x <- .target_chord_label_x(labeled$x, min_gap, lower = 1, upper = max(n, 1))
    labeled$short_label <- .plot_truncate(labeled$label, 20)
    p <- p +
      ggplot2::geom_segment(
        data = labeled, ggplot2::aes(x = .data$x, xend = .data$label_x, y = 0, yend = -drop),
        colour = "grey55", linewidth = 0.25
      ) +
      ggplot2::geom_text(
        data = labeled, ggplot2::aes(x = .data$label_x, y = -drop, label = .data$short_label),
        angle = 60, hjust = 1, vjust = 0.5, size = 2.6, colour = "grey15"
      )
  }

  p +
    ggplot2::coord_cartesian(ylim = c(y_bottom, 1.03 * y_top), expand = FALSE, clip = "off") +
    ggplot2::labs(
      title = paste0("Target-target STRING actions adjacency -- ", title_suffix),
      subtitle = .plot_wrap(paste0(
        "Nodes = condition target(s), ordered by degree (ascending, left to right); arcs = STRING directed 'actions' interactions ",
        "(network_bowtie()'s own data source, direction dropped here); top labeled by degree, rest ticked"
      ), .plot_wrap_width(fig_width, 7.5)),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40"),
      plot.title.position = "plot",
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}
