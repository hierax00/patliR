#' @include AllGenerics.R internal.R network_build.R
NULL

## network_module_robustness() -- module detection (dbscan::hdbscan()) +
## targeted-attack percolation (Schneider et al. 2011 R-index) per module,
## on network_build()'s per-condition graph. Verified dbscan::hdbscan()'s
## real signature against CRAN docs before writing this (2026-07-23):
## hdbscan(x, minPts, ...) where x is a numeric matrix or a dist object;
## there is no default for minPts (unlike some clustering functions), so it
## is exposed here as `min_module_size` rather than silently hardcoded.

#' Detect network modules and measure each module's robustness to targeted
#' node removal
#'
#' @description
#' For each condition's graph ([network_build()]), clusters nodes into
#' modules using HDBSCAN (`dbscan::hdbscan()`) over the graph's shortest-path
#' distance matrix (`igraph::distances()`) -- nodes that are topologically
#' close (few hops apart) cluster into the same module, regardless of
#' whether they are compounds or targets. For each module, runs a
#' **targeted-attack percolation**: repeatedly remove the current
#' highest-degree node from that module's induced subgraph (recalculating
#' degree after every removal -- the standard "malicious attack" strategy,
#' not random failure), tracking what fraction of the module stays in its
#' largest connected component at each step. Summarizes that fragmentation
#' curve as the **R-index** (Schneider, C.M. et al. 2011, *PNAS* 108(10),
#' 3838-41 -- see `patliR_manual.md`, section 9, where `plot_robustness()`
#' will visualize the curve this function computes): the mean largest-
#' component fraction across all removal steps, from `0` (module collapses
#' immediately) to `1` (module stays fully connected however much you
#' remove -- only possible for a single-node module).
#'
#' @section HDBSCAN "noise" points get their own module, tagged as such (2026-08-14):
#' HDBSCAN labels points that do not belong to any dense cluster as cluster
#' `0` ("noise"). Earlier versions **excluded** these nodes entirely from
#' module-level robustness output. On real, messier compound-target graphs
#' this silently dropped a meaningful fraction of nodes from every
#' downstream consumer of [patliRResults()]`(proj, "network_module_robustness")`
#' with no trace beyond a log line -- a node that never clustered cleanly is
#' not the same thing as a node that does not exist, and hiding it entirely
#' made "how many nodes does this condition actually cover" silently wrong.
#' Now every noise node in a component is collected into its own module
#' (`module_type = "noise"` in the summary table, vs. `"cluster"` for
#' HDBSCAN's real clusters and the small-graph/no-clusters fallback module),
#' and its robustness is still computed the same way as any other module --
#' a percolation curve over a scattered set of nodes is a legitimate (if
#' less structurally meaningful) thing to measure, and consistently
#' including it beats a second special case. Components that produce no
#' noise at all simply have no `"noise"` module row. How many nodes were
#' noise, and which, is always logged either way.
#'
#' @section Small graphs may not have real modules:
#' If a connected component has fewer than `2 * min_module_size` nodes,
#' HDBSCAN is not run on it at all (not enough points to distinguish
#' "cluster" from "noise" meaningfully) -- that **component** is treated as
#' a single module instead, and this is logged. This matters in practice:
#' `patliR`'s own bundled example data produces graphs this small.
#'
#' @section Disconnected graphs are clustered per connected component:
#' Real compound-target networks (unlike `patliR`'s small bundled example)
#' are frequently **disconnected** -- a target reachable only through
#' compounds absent from the current condition, or an isolated
#' compound-target pair with no path to the rest of the graph. HDBSCAN
#' requires a distance matrix, and `igraph::distances()` returns `Inf` for
#' any pair of nodes in different components; `dbscan::hdbscan()` does not
#' validate its input for non-finite values and crashes on it (confirmed
#' against its own source, 2026-08-11 -- not documented, just reproduced).
#' Rather than special-casing `Inf`, `.network_detect_modules()` first
#' splits the graph into its connected components (`igraph::components()`)
#' and runs the small-graph fallback and/or HDBSCAN **within each component
#' separately** -- cross-component node pairs were never meaningfully
#' clusterable together anyway, since no path connects them. Module IDs
#' (`M1`, `M2`, ...) are assigned sequentially across all components, and
#' the log message reports each component's outcome.
#'
#' @inheritParams network_build
#' @param clustering Only `"hdbscan"` is implemented.
#' @param min_module_size Minimum HDBSCAN cluster size (`dbscan::hdbscan()`'s
#'   `minPts`, which has no built-in default -- exposed here rather than
#'   silently hardcoded). Default `2`.
#' @param seed Integer seed for tie-breaking when multiple nodes share the
#'   current maximum degree during the attack sequence (`sample()` picks
#'   among ties) -- always set, exposed, and logged, even though the tie-
#'   breaking step is a small part of an otherwise deterministic algorithm,
#'   for exact reproducibility. `NULL` (default) draws and logs a fresh
#'   seed.
#'
#' @return The updated `proj`, with two entries in [patliRResults()]:
#'   `network_module_robustness` (columns `condition`, `module_id`,
#'   `module_type` (`"cluster"`/`"noise"` -- see "HDBSCAN 'noise' points get
#'   their own module" above), `n_nodes`, `r_index`, `seed_used`) and
#'   `network_robustness_curve`
#'   (columns `condition`, `module_id`, `n_removed`,
#'   `largest_component_fraction` -- the raw percolation curve behind each
#'   `r_index`, for `plot_robustness()` later). Both written to their
#'   matching `results/*.csv`.
#'
#' @examples
#' \donttest{
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
#' proj <- network_module_robustness(proj, condition = "FLO-ET")
#' patliRResults(proj, "network_module_robustness")
#' }
#'
#' @export
network_module_robustness <- function(proj, condition = NULL,
                                       clustering = c("hdbscan"),
                                       min_module_size = 2, seed = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  clustering <- match.arg(clustering)
  stopifnot(is.numeric(min_module_size), length(min_module_size) == 1, min_module_size >= 2)
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_module_robustness} needs the {.pkg dbscan} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"dbscan\")} -- pure CRAN, no Bioconductor/Java involved."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  stopifnot(is.numeric(seed), length(seed) == 1)

  summary_rows <- vector("list", length(conditions))
  curve_rows <- vector("list", length(conditions))
  names(summary_rows) <- conditions
  names(curve_rows) <- conditions

  for (cond in conditions) {
    g <- .network_graph(proj, cond)
    n_nodes <- igraph::vcount(g)

    if (n_nodes < 2) {
      proj <- .log_append(
        proj, step = "network_module_robustness", id = NA_character_,
        message = paste0("condition '", cond, "': fewer than 2 nodes, nothing to cluster/percolate")
      )
      summary_rows[[cond]] <- .empty_network_module_robustness_row()
      curve_rows[[cond]] <- .empty_network_robustness_curve_row()
      next
    }

    modules <- .network_detect_modules(g, min_module_size)
    proj <- .log_append(
      proj, step = "network_module_robustness", id = NA_character_,
      message = modules$log_message
    )

    cond_summary <- vector("list", length(modules$members))
    cond_curve <- vector("list", length(modules$members))
    for (i in seq_along(modules$members)) {
      module_nodes <- modules$members[[i]]
      module_id <- names(modules$members)[i]
      module_type <- modules$types[[i]]
      sub_g <- igraph::induced_subgraph(g, igraph::V(g)[module_nodes])

      perc <- .network_percolate(sub_g, seed = seed)
      cond_summary[[i]] <- data.frame(
        condition = cond, module_id = module_id, module_type = module_type,
        n_nodes = igraph::vcount(sub_g),
        r_index = perc$r_index, seed_used = seed, stringsAsFactors = FALSE
      )
      cond_curve[[i]] <- data.frame(
        condition = cond, module_id = module_id, n_removed = seq_along(perc$curve) - 1L,
        largest_component_fraction = perc$curve, stringsAsFactors = FALSE
      )
    }
    summary_rows[[cond]] <- do.call(rbind, cond_summary)
    curve_rows[[cond]] <- do.call(rbind, cond_curve)
  }

  summary_result <- do.call(rbind, summary_rows)
  rownames(summary_result) <- NULL
  summary_result <- .network_upsert(proj, "network_module_robustness", summary_result, "condition")

  curve_result <- do.call(rbind, curve_rows)
  rownames(curve_result) <- NULL
  curve_result <- .network_upsert(proj, "network_robustness_curve", curve_result, "condition")

  patliRResults(proj, "network_module_robustness") <- summary_result
  patliRResults(proj, "network_robustness_curve") <- curve_result
  .write_results_csv(proj, "network_module_robustness", summary_result)
  .write_results_csv(proj, "network_robustness_curve", curve_result)
  .write_log_csv(proj)
  proj
}

#' Cluster a graph's nodes into modules via HDBSCAN over shortest-path
#' distances, run separately per connected component (cross-component node
#' pairs have infinite graph distance and are never meaningfully
#' clusterable together -- see "Disconnected graphs are clustered per
#' connected component" above), falling back to "one module = whole
#' component" when a component is too small
#' @return `list(members = named list of character vectors of node names,
#'   types = character vector, same length/order as `members`, each
#'   `"cluster"` or `"noise"` -- see "HDBSCAN 'noise' points get their own
#'   module" above, log_message = character scalar)`.
#' @keywords internal
.network_detect_modules <- function(g, min_module_size) {
  comps <- igraph::components(g)
  n_components <- comps$no

  all_members <- list()
  all_types <- character(0)
  comp_messages <- character(n_components)

  for (comp_id in seq_len(n_components)) {
    comp_nodes <- names(comps$membership)[comps$membership == comp_id]
    sub_g <- igraph::induced_subgraph(g, igraph::V(g)[comp_nodes])
    comp_result <- .network_detect_modules_component(sub_g, min_module_size)
    all_members <- c(all_members, comp_result$members)
    all_types <- c(all_types, comp_result$types)
    comp_messages[comp_id] <- comp_result$log_message
  }

  names(all_members) <- paste0("M", seq_along(all_members))
  names(all_types) <- names(all_members)

  log_message <- if (n_components > 1) {
    paste0(
      n_components, " connected component(s) in graph (module detection run separately per ",
      "component, since cross-component nodes have no path between them): ",
      paste(comp_messages, collapse = " | ")
    )
  } else {
    comp_messages[1]
  }

  list(members = all_members, types = all_types, log_message = log_message)
}

#' Cluster a single connected component's nodes into modules via HDBSCAN
#' over shortest-path distances (always finite within one component),
#' falling back to "one module = whole component" when too small, and
#' collecting any HDBSCAN "noise" nodes into their own module rather than
#' dropping them
#' @return `list(members = named list of character vectors of node names,
#'   types = character vector, same length/order as `members`, each
#'   `"cluster"` or `"noise"`, log_message = character scalar)`.
#' @keywords internal
.network_detect_modules_component <- function(g, min_module_size) {
  node_names <- igraph::V(g)$name
  n_nodes <- length(node_names)

  if (n_nodes < 2 * min_module_size) {
    return(list(
      members = stats::setNames(list(node_names), "M1"),
      types = "cluster",
      log_message = paste0(
        "component with ", n_nodes, " node(s) (< 2 * min_module_size = ", 2 * min_module_size,
        "); treated as a single module instead of running HDBSCAN"
      )
    ))
  }

  dist_matrix <- igraph::distances(g)
  hdb <- dbscan::hdbscan(stats::as.dist(dist_matrix), minPts = min_module_size)

  cluster_ids <- sort(unique(hdb$cluster[hdb$cluster != 0]))
  noise_nodes <- node_names[hdb$cluster == 0]
  n_noise <- length(noise_nodes)

  if (length(cluster_ids) == 0) {
    ## Every point is HDBSCAN noise -- there is no real cluster here, so
    ## (unlike the mixed case below) the whole component becomes a single
    ## module, tagged "noise" rather than "cluster" since that is what it
    ## actually is.
    return(list(
      members = stats::setNames(list(node_names), "M1"),
      types = "noise",
      log_message = paste0(
        "HDBSCAN found no clusters in a ", n_nodes, "-node component (all noise); treated as a single noise module"
      )
    ))
  }

  members <- lapply(cluster_ids, function(cl) node_names[hdb$cluster == cl])
  names(members) <- paste0("M", cluster_ids)
  types <- rep("cluster", length(members))

  if (n_noise > 0) {
    members <- c(members, stats::setNames(list(noise_nodes), "Mnoise"))
    types <- c(types, "noise")
  }

  list(
    members = members,
    types = types,
    log_message = paste0(
      length(cluster_ids), " module(s) found via HDBSCAN in a ", n_nodes, "-node component (minPts = ",
      min_module_size, "); ", n_noise, " node(s) placed in a separate noise module (not dropped)"
    )
  )
}

#' Targeted-attack percolation on one (sub)graph: repeatedly remove the
#' current highest-degree node, track the largest-component fraction
#'
#' @description
#' `r_index` is Schneider et al. (2011), *PNAS* 108(10), 3838-3841's
#' robustness index: the mean of `s(Q)` (largest-component fraction) over
#' removal steps `Q = 1..N` -- **not** including the trivial, pre-removal
#' `s(0) = 1` state. `curve[1] == 1` (the untouched graph) is still kept in
#' the returned `curve` for plotting (so [plot_robustness()] can draw the
#' curve starting from "nothing removed yet"), but must be dropped before
#' averaging into `r_index` -- averaging it in was a real bug (found
#' 2026-08-25): it systematically inflates `r_index`, worst for exactly the
#' small modules this package's networks typically have after HDBSCAN
#' splits a component (e.g. a 1-node module used to report `r_index = 0.5`
#' via `mean(c(1, 0))`, when Schneider's own definition for `N = 1` is
#' `s(1) = 0`).
#' @return `list(curve = numeric vector, length vcount(g)+1, curve[1] == 1;
#'   r_index = mean(curve[-1]))`.
#' @keywords internal
.network_percolate <- function(g, seed) {
  n <- igraph::vcount(g)
  if (n == 1) {
    return(list(curve = c(1, 0), r_index = 0))
  }

  withr_seed <- .with_seed(seed)
  on.exit(withr_seed(), add = TRUE)

  remaining <- g
  curve <- numeric(n + 1)
  curve[1] <- 1 ## before removing anything, the whole module is "connected to itself"

  for (step in seq_len(n)) {
    deg <- igraph::degree(remaining)
    max_deg <- max(deg)
    candidates <- names(deg)[deg == max_deg]
    target <- if (length(candidates) > 1) sample(candidates, 1) else candidates
    remaining <- igraph::delete_vertices(remaining, target)

    if (igraph::vcount(remaining) == 0) {
      largest_cc <- 0
    } else {
      comp <- igraph::components(remaining)
      largest_cc <- max(comp$csize)
    }
    curve[step + 1] <- largest_cc / n
  }

  ## r_index averages steps 1..n only (Schneider's Q = 1..N) -- curve[1] is
  ## the pre-removal Q = 0 state and must not be included, see above.
  list(curve = curve, r_index = mean(curve[-1]))
}

#' Temporarily set the RNG seed, returning a closure that restores the
#' previous RNG state when called -- avoids `network_module_robustness()`
#' permanently perturbing the caller's random state
#' @keywords internal
.with_seed <- function(seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

#' @keywords internal
.empty_network_module_robustness_row <- function() {
  data.frame(condition = character(0), module_id = character(0), module_type = character(0),
             n_nodes = integer(0), r_index = double(0), seed_used = integer(0), stringsAsFactors = FALSE)
}

#' @keywords internal
.empty_network_robustness_curve_row <- function() {
  data.frame(condition = character(0), module_id = character(0), n_removed = integer(0),
             largest_component_fraction = double(0), stringsAsFactors = FALSE)
}
