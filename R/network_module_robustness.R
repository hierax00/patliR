#' @include AllGenerics.R internal.R network_build.R
NULL

## Module detection + per-module percolation robustness on network_build()'s
## per-condition graph. Three clustering backends:
##
##   "leiden"   -- igraph::cluster_leiden(objective_function = "modularity"),
##                 per connected component. The general default.
##   "bipartite"-- bipartite::computeModules() (Beckett 2016), maximising
##                 Barber (2007)'s bipartite modularity Q_B. Only valid on
##                 the genuinely two-mode compound-target graph.
##   "hdbscan"  -- dbscan::hdbscan() over the shortest-path distance matrix.
##                 A density heuristic on an integer metric; kept as a
##                 no-edge-density-assumption contrast, not the default.
##
## Then, per module, a percolation attack (targeted degree-descending and/or
## random) summarised as Schneider et al. (2011)'s R-index.

#' Detect network modules and measure each module's robustness to node
#' removal
#'
#' @description
#' For each condition's graph ([network_build()]), clusters nodes into
#' modules, then runs a **percolation attack** per module and reports each
#' module's **R-index** (Schneider, C.M. et al. 2011, *PNAS* 108(10),
#' 3838-41): the mean largest-connected-component fraction over removal
#' steps `Q = 1..N`. Higher = more robust; `0` is immediate collapse. The
#' theoretical maximum is `(N-1)/(2N) < 0.5` (Schneider's definition),
#' reached only by a module that never fragments until fully removed -- not
#' `1`. R-index depends strongly on module size, so it is **not comparable
#' across modules of different `n_nodes`** (a 1-node module scores `0`, a
#' 2-node one `0.25`, by arithmetic).
#'
#' @section Clustering backends (`clustering=`):
#' * `"leiden"` (default) -- [igraph::cluster_leiden()] with
#'   `objective_function = "modularity"`, run per connected component
#'   (reusing the same component loop as the other backends, so the module
#'   ID scheme and the `sum(n_nodes) == vcount(g)` invariant are shared).
#'   Makes no bipartiteness assumption and guarantees well-connected
#'   communities. `module_type` is `"cluster"` for every Leiden module.
#' * `"bipartite"` -- `bipartite::computeModules(method = "Beckett")`
#'   (Beckett 2016's DIRTLPAwb+), maximising Barber (2007)'s bipartite
#'   modularity `Q_B`. Standard (Newman) modularity's null model permits
#'   compound-compound and target-target edges, which cannot exist in a
#'   compound-target graph; Barber's null does not. **Only valid on the
#'   two-mode compound-target graph** -- aborts (via `.network_is_bipartite()`)
#'   otherwise. Needs the `bipartite` package (`Suggests`).
#' * `"hdbscan"` -- `dbscan::hdbscan()` over `igraph::distances(g, weights =
#'   NA)`. With `minPts = 2` every non-isolated node has core distance 1 and
#'   the mutual-reachability hierarchy degenerates to the single-linkage
#'   hierarchy of an integer metric with 3-4 distinct values, so the
#'   "density" it estimates is near-constant on a bipartite graph. Kept as
#'   an explicit contrast (it assumes nothing about edge density), not the
#'   default -- selecting it emits a one-time note. HDBSCAN "noise" nodes
#'   (label `0`) are collected into their own module (`module_type =
#'   "noise"`) rather than dropped.
#'
#' Clustering is run **unweighted**. For `"leiden"` the `weight` edge
#' attribute (the import probability) is stripped with
#' `igraph::delete_edge_attr()` before clustering -- *not* passed as
#' `weights = NA`, which would leave `strength()` silently reading the
#' attribute for the modularity null model and, on any `NA` probability,
#' return every node as its own module with no error. The probabilities are
#' also not calibrated across import platforms (a weighted `Q_B` would be a
#' modularity of prediction confidence, not of pharmacology), and every
#' other topological call in `network_*` is unweighted.
#'
#' @section Percolation attack (`attack=`):
#' * `"targeted"` (default) -- repeatedly remove the current highest-degree
#'   node from the module's induced subgraph, recomputing degree after every
#'   removal (the standard "malicious attack"). Deterministic apart from the
#'   `seed`-controlled tie-break among equal-degree nodes.
#' * `"random"` -- `n_random` independent uniformly-random removal orders
#'   (no degree recomputation -- "random failure"). The stored curve is the
#'   **mean** largest-component fraction across replicates at each step
#'   (`largest_component_sd` its SD); `r_index_random` is the mean over
#'   replicates of each replicate's R-index (equivalently
#'   `mean(mean_curve[-1])`, since the mean is linear). Each replicate draws
#'   under a derived sub-seed combining `seed`, the module index and the
#'   replicate index (computed overflow-safe, so a near-`.Machine$integer.max`
#'   seed does not wrap to `NA`).
#' * `"both"` -- both; `plot_robustness()` then draws the two-curve
#'   percolation figure (Albert, Jeong & Barabasi 2000).
#'
#' `r_index` always holds the targeted value (or `NA` if `attack =
#' "random"`); `r_index_random` holds the random baseline (`NA` unless
#' `attack` includes it). The interpretable quantity is `r_index_random -
#' r_index`: a module far below its random baseline is hub-dependent.
#'
#' @section Small and disconnected graphs:
#' Real compound-target networks are frequently disconnected;
#' `igraph::distances()` returns `Inf` across components. All three backends
#' run **per connected component** ([igraph::components()]); a component
#' with fewer than `min_component_size` nodes is treated as a single module
#' (no clustering run on it), logged. Module IDs (`M1`, `M2`, ...) are
#' assigned sequentially across components and are a **within-run label
#' only** -- `M3` in one run is unrelated to `M3` under a different `seed`,
#' `resolution`, `clustering`, or rebuilt graph.
#'
#' @inheritParams network_build
#' @param clustering `"leiden"` (default), `"bipartite"`, or `"hdbscan"` --
#'   see "Clustering backends" above. **The default changed from
#'   `"hdbscan"` in patliR 0.2.0**; the module partitions, and therefore
#'   every `r_index`, will differ from earlier versions.
#' @param attack `"targeted"` (default), `"random"`, or `"both"` -- see
#'   "Percolation attack" above.
#' @param resolution Leiden resolution (Reichardt-Bornholdt gamma,
#'   internally rescaled by `1 / sum(vertex_weights)`; *not* a CPM density
#'   threshold). Higher = more, smaller modules. Default `1`. Ignored unless
#'   `clustering = "leiden"`.
#' @param n_iterations Leiden iterations (Traag et al. 2019 recommend
#'   iterating to stability; igraph's own default is only `2`). Default `5`.
#'   Ignored unless `clustering = "leiden"`.
#' @param min_module_size HDBSCAN `minPts` (`dbscan::hdbscan()` has no
#'   built-in default). Default `2`. Ignored unless `clustering =
#'   "hdbscan"`.
#' @param min_component_size Connected components smaller than this are
#'   treated as a single module by every backend (the method-agnostic
#'   fallback). Default `3`.
#' @param n_random Number of random removal orders averaged for `attack =
#'   "random"`/`"both"`. Default `20`; **must be `>= 2`** for those attacks
#'   (a single draw has no SD and is not a baseline -- it inverts the
#'   targeted-vs-random conclusion for a non-trivial fraction of modules),
#'   so the baseline is always stored as an aggregate, never one replicate.
#' @param seed Integer seed covering (independently) the stochastic
#'   clustering step (Leiden / Beckett label propagation) and the
#'   percolation tie-break / random orders. Always set, exposed and logged.
#'   `NULL` (default) draws and logs a fresh seed.
#'
#' @return The updated `proj`, with three entries in [patliRResults()]:
#'   * `network_module_robustness` -- one row per module: `condition`,
#'     `module_id`, `module_type` (`"cluster"`/`"noise"`), `n_nodes`,
#'     `r_index` (targeted), `r_index_random`, `n_random`, `clustering`,
#'     `modularity` (Newman `Q` for Leiden via `igraph::modularity()`;
#'     node-weighted mean of per-component Barber `Q_B` for bipartite; `NA`
#'     for hdbscan -- repeated down the condition's rows), `resolution`,
#'     `method_detail`, `seed_used`.
#'   * `network_robustness_curve` -- the raw percolation curve behind each
#'     `r_index`: `condition`, `module_id`, `n_removed`,
#'     `largest_component_fraction`, `removal_strategy`
#'     (`"targeted"`/`"random"`), `largest_component_sd` (`NA` for
#'     targeted), `n_replicates` (`1` for targeted, `n_random` for random).
#'   * `network_module_membership` -- one row per node: `condition`,
#'     `node_id`, `node_type` (`"compound"`/`"target"`), `module_id`,
#'     `module_type`. To be consumed by a future `plot_network_layers()`
#'     module-colour mode (spec 3.2).
#'   All three written to their matching `results/*.csv`. Rerunning with a
#'   different `clustering` **replaces** the condition's rows (the
#'   `clustering` column records which method produced them); it does not
#'   accumulate.
#'
#' @references
#' Schneider, C.M., Moreira, A.A., Andrade, J.S., Havlin, S. & Herrmann,
#' H.J. (2011) Mitigation of malicious attacks on networks. *PNAS*
#' 108(10):3838-3841.
#'
#' Traag, V.A., Waltman, L. & van Eck, N.J. (2019) From Louvain to Leiden:
#' guaranteeing well-connected communities. *Sci Rep* 9:5233.
#'
#' Barber, M.J. (2007) Modularity and community detection in bipartite
#' networks. *Phys Rev E* 76:066102.
#'
#' Beckett, S.J. (2016) Improved community detection in weighted bipartite
#' networks. *R Soc Open Sci* 3:140536.
#'
#' Albert, R., Jeong, H. & Barabasi, A.-L. (2000) Error and attack
#' tolerance of complex networks. *Nature* 406:378-382.
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
#' proj <- network_module_robustness(proj, condition = "FLO-ET", attack = "both")
#' patliRResults(proj, "network_module_robustness")
#' }
#'
#' @export
network_module_robustness <- function(proj, condition = NULL,
                                       clustering = c("leiden", "bipartite", "hdbscan"),
                                       attack = c("targeted", "random", "both"),
                                       resolution = 1, n_iterations = 5L,
                                       min_module_size = 2, min_component_size = 3L,
                                       n_random = 20L, seed = NULL) {
  stopifnot(is(proj, "PatliRProject"))
  clustering <- match.arg(clustering)
  attack <- match.arg(attack)
  stopifnot(is.numeric(min_module_size), length(min_module_size) == 1, min_module_size >= 2)
  stopifnot(is.numeric(min_component_size), length(min_component_size) == 1, min_component_size >= 2)
  stopifnot(is.numeric(resolution), length(resolution) == 1, resolution > 0)
  stopifnot(is.numeric(n_iterations), length(n_iterations) == 1, n_iterations >= 1)
  stopifnot(is.numeric(n_random), length(n_random) == 1, n_random >= 1)
  n_iterations <- as.integer(n_iterations)
  n_random <- as.integer(n_random)
  if (attack %in% c("random", "both") && n_random < 2L) {
    cli::cli_abort(c(
      "{.arg n_random} must be >= 2 when {.arg attack} includes {.val random}.",
      "i" = "A single random removal order is not a baseline -- its per-step SD is undefined ({.field largest_component_sd} would be all {.val NA}). The default is 20."
    ))
  }

  if (clustering == "hdbscan") {
    cli::cli_inform(c(
      "i" = paste0(
        "{.arg clustering = \"hdbscan\"} is a density heuristic on an integer graph-distance ",
        "metric (it degenerates to single linkage at {.code minPts = 2}); {.arg clustering = ",
        "\"leiden\"} is the general default. See {.help network_module_robustness}."
      )
    ), .frequency = "once", .frequency_id = "patliR_hdbscan_clustering")
  }

  if (clustering == "hdbscan" && !requireNamespace("dbscan", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_module_robustness} with {.arg clustering = \"hdbscan\"} needs the {.pkg dbscan} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"dbscan\")} -- pure CRAN, no Bioconductor/Java involved."
    ))
  }
  if (clustering == "bipartite" && !requireNamespace("bipartite", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_module_robustness} with {.arg clustering = \"bipartite\"} needs the {.pkg bipartite} package (CRAN), not installed.",
      "i" = "{.code install.packages(\"bipartite\")} -- pure CRAN, no Bioconductor/Java involved. patliR only ever loads its namespace ({.fn requireNamespace}), which does not attach anything."
    ))
  }

  conditions <- .network_resolve_conditions(proj, condition)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  stopifnot(is.numeric(seed), length(seed) == 1)
  seed <- as.integer(seed)

  summary_rows <- vector("list", length(conditions))
  curve_rows <- vector("list", length(conditions))
  membership_rows <- vector("list", length(conditions))
  names(summary_rows) <- conditions
  names(curve_rows) <- conditions
  names(membership_rows) <- conditions

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
      membership_rows[[cond]] <- .empty_network_module_membership_row()
      next
    }

    if (clustering == "bipartite" && !.network_is_bipartite(g)) {
      cli::cli_abort(c(
        "{.arg clustering = \"bipartite\"} needs a genuine two-mode graph, but condition {.val {cond}}'s graph is not one.",
        "i" = "It has vertices with no {.code type} or an edge within one mode. Barber's bipartite modularity only applies to the compound-target graph {.fn network_build} builds, not the layered graph."
      ))
    }

    modules <- .network_detect_modules(
      g, clustering = clustering, min_module_size = min_module_size,
      resolution = resolution, n_iterations = n_iterations,
      min_component_size = min_component_size, seed = seed
    )
    proj <- .log_append(
      proj, step = "network_module_robustness", id = NA_character_,
      message = modules$log_message
    )
    if (!is.null(modules$warn_message)) {
      cli::cli_warn(.cli_escape(modules$warn_message))
    }

    node_type <- stats::setNames(.network_node_types(g), igraph::V(g)$name)

    cond_summary <- vector("list", length(modules$members))
    cond_curve <- vector("list", length(modules$members))
    cond_membership <- vector("list", length(modules$members))
    for (i in seq_along(modules$members)) {
      module_nodes <- modules$members[[i]]
      module_id <- names(modules$members)[i]
      module_type <- modules$types[[i]]
      sub_g <- igraph::induced_subgraph(g, igraph::V(g)[module_nodes])

      r_index_targeted <- NA_real_
      r_index_random <- NA_real_
      strat_curves <- list()

      if (attack %in% c("targeted", "both")) {
        perc <- .network_percolate(sub_g, seed = seed)
        r_index_targeted <- perc$r_index
        strat_curves$targeted <- data.frame(
          n_removed = seq_along(perc$curve) - 1L,
          largest_component_fraction = perc$curve,
          removal_strategy = "targeted",
          largest_component_sd = NA_real_,
          n_replicates = 1L, stringsAsFactors = FALSE
        )
      }
      if (attack %in% c("random", "both")) {
        ## Reserve enough seeds for every replicate: a fixed stride of
        ## 1000 reused the preceding module's draws when n_random > 1000.
        perc_r <- .network_percolate_random(sub_g, seed = .derive_seed(seed, max(1000, n_random) * i), n_random = n_random)
        r_index_random <- perc_r$r_index
        strat_curves$random <- data.frame(
          n_removed = seq_along(perc_r$curve) - 1L,
          largest_component_fraction = perc_r$curve,
          removal_strategy = "random",
          largest_component_sd = perc_r$curve_sd,
          n_replicates = perc_r$n_replicates, stringsAsFactors = FALSE
        )
      }

      cond_summary[[i]] <- data.frame(
        condition = cond, module_id = module_id, module_type = module_type,
        n_nodes = igraph::vcount(sub_g),
        r_index = r_index_targeted, r_index_random = r_index_random,
        n_random = if (attack %in% c("random", "both")) n_random else NA_integer_,
        clustering = clustering,
        modularity = modules$modularity,
        resolution = if (clustering == "leiden") resolution else NA_real_,
        method_detail = modules$method_detail,
        seed_used = seed, stringsAsFactors = FALSE
      )
      curve_df <- do.call(rbind, strat_curves)
      cond_curve[[i]] <- data.frame(
        condition = cond, module_id = module_id, curve_df,
        stringsAsFactors = FALSE, row.names = NULL
      )
      cond_membership[[i]] <- data.frame(
        condition = cond, node_id = module_nodes,
        node_type = unname(node_type[module_nodes]),
        module_id = module_id, module_type = module_type,
        stringsAsFactors = FALSE
      )
    }
    summary_rows[[cond]] <- do.call(rbind, cond_summary)
    curve_rows[[cond]] <- do.call(rbind, cond_curve)
    membership_rows[[cond]] <- do.call(rbind, cond_membership)
  }

  touched_keys <- data.frame(condition = conditions, stringsAsFactors = FALSE)

  summary_result <- do.call(rbind, summary_rows)
  rownames(summary_result) <- NULL
  summary_result <- .network_upsert(proj, "network_module_robustness", summary_result, "condition", touched_keys = touched_keys)

  curve_result <- do.call(rbind, curve_rows)
  rownames(curve_result) <- NULL
  curve_result <- .network_upsert(proj, "network_robustness_curve", curve_result, "condition", touched_keys = touched_keys)

  membership_result <- do.call(rbind, membership_rows)
  rownames(membership_result) <- NULL
  membership_result <- .network_upsert(proj, "network_module_membership", membership_result, "condition", touched_keys = touched_keys)

  patliRResults(proj, "network_module_robustness") <- summary_result
  patliRResults(proj, "network_robustness_curve") <- curve_result
  patliRResults(proj, "network_module_membership") <- membership_result
  .write_results_csv(proj, "network_module_robustness", summary_result)
  .write_results_csv(proj, "network_robustness_curve", curve_result)
  .write_results_csv(proj, "network_module_membership", membership_result)
  .write_log_csv(proj)
  proj
}

#' Cluster a graph's nodes into modules, run separately per connected
#' component (cross-component node pairs are never meaningfully clusterable
#' together), falling back to "one module = whole component" when a
#' component has fewer than `min_component_size` nodes
#'
#' @description
#' Dispatches to one of three per-component backends (`"leiden"`,
#' `"bipartite"`, `"hdbscan"`). The stochastic clustering step is wrapped in
#' `.with_seed(seed)` here -- independently of the per-module
#' `.with_seed(seed)` inside `.network_percolate()`, so `clustering =
#' "hdbscan"` (which uses no RNG) still reproduces its historical output
#' bit-for-bit.
#' @return `list(members = named list of character vectors of node names,
#'   types = character vector same length/order as `members`,
#'   log_message = character scalar, warn_message = character scalar or
#'   NULL, modularity = numeric scalar (Newman Q / mean Barber Q_B / NA),
#'   method_detail = character scalar)`.
#' @keywords internal
.network_detect_modules <- function(g, clustering, min_module_size = 2,
                                     resolution = 1, n_iterations = 5L,
                                     min_component_size = 3L, seed = NULL) {
  restore_seed <- .with_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  comps <- igraph::components(g)
  n_components <- comps$no

  all_members <- list()
  all_types <- character(0)
  comp_messages <- character(n_components)
  warn_messages <- character(0)
  qb_values <- numeric(0)
  qb_weights <- numeric(0)

  method_detail <- switch(
    clustering,
    leiden = "modularity",
    bipartite = "Beckett",
    hdbscan = paste0("hdbscan (minPts=", min_module_size, ")")
  )

  for (comp_id in seq_len(n_components)) {
    comp_nodes <- names(comps$membership)[comps$membership == comp_id]
    sub_g <- igraph::induced_subgraph(g, igraph::V(g)[comp_nodes])
    n_c <- length(comp_nodes)

    if (n_c < min_component_size) {
      cr <- list(
        members = stats::setNames(list(comp_nodes), "M1"), types = "cluster",
        log_message = paste0(
          "component with ", n_c, " node(s) (< min_component_size = ", min_component_size,
          "); treated as a single module instead of running ", clustering
        )
      )
    } else {
      cr <- switch(
        clustering,
        leiden = .network_detect_modules_leiden(sub_g, resolution, n_iterations),
        bipartite = .network_detect_modules_bipartite(sub_g, seed),
        hdbscan = .network_detect_modules_hdbscan(sub_g, min_module_size)
      )
    }

    all_members <- c(all_members, cr$members)
    all_types <- c(all_types, cr$types)
    comp_messages[comp_id] <- cr$log_message
    if (!is.null(cr$warn_message)) warn_messages <- c(warn_messages, cr$warn_message)
    if (!is.null(cr$qb) && !is.na(cr$qb)) {
      qb_values <- c(qb_values, cr$qb)
      qb_weights <- c(qb_weights, n_c)
    }
  }

  names(all_members) <- paste0("M", seq_along(all_members))
  names(all_types) <- names(all_members)

  modularity_value <- switch(
    clustering,
    hdbscan = NA_real_,
    bipartite = if (length(qb_values) == 0) NA_real_ else stats::weighted.mean(qb_values, qb_weights),
    leiden = {
      memb_vec <- integer(igraph::vcount(g))
      names(memb_vec) <- igraph::V(g)$name
      for (k in seq_along(all_members)) memb_vec[all_members[[k]]] <- k
      q <- igraph::modularity(g, memb_vec, weights = NA)
      if (is.nan(q)) NA_real_ else q
    }
  )

  log_message <- if (n_components > 1) {
    paste0(
      n_components, " connected component(s) in graph (module detection run separately per ",
      "component): ", paste(comp_messages, collapse = " | ")
    )
  } else {
    comp_messages[1]
  }

  list(
    members = all_members, types = all_types, log_message = log_message,
    warn_message = if (length(warn_messages) == 0) NULL else paste(warn_messages, collapse = " | "),
    modularity = modularity_value, method_detail = method_detail
  )
}

#' Leiden community detection on one connected component, unweighted
#'
#' @description
#' The `weight` edge attribute is **deleted** (not passed as `weights =
#' NA`): `igraph::cluster_leiden()` reads `E(g)$weight` via `strength()` for
#' the modularity null model even when `weights = NA` is passed, and any
#' `NA` weight then silently yields all singletons. `objective_function =
#' "modularity"` is set explicitly because igraph's default (`"CPM"`)
#' returns all singletons at `resolution = 1` on a sparse graph.
#' @return `list(members, types, log_message, warn_message, qb = NA_real_)`.
#' @keywords internal
.network_detect_modules_leiden <- function(g, resolution, n_iterations) {
  g2 <- if ("weight" %in% igraph::edge_attr_names(g)) igraph::delete_edge_attr(g, "weight") else g
  comm <- igraph::cluster_leiden(
    g2, objective_function = "modularity",
    resolution = resolution, n_iterations = n_iterations
  )
  memb <- igraph::membership(comm)
  node_names <- igraph::V(g)$name
  ids <- sort(unique(memb))
  members <- lapply(ids, function(k) node_names[memb == k])
  names(members) <- paste0("M", ids)

  warn_message <- NULL
  if (igraph::ecount(g) > 0 && length(members) == igraph::vcount(g)) {
    warn_message <- paste0(
      "Leiden (objective_function = 'modularity', resolution = ", resolution,
      ") returned every node as its own module on a ", igraph::vcount(g),
      "-node component -- the component may be near-tree; try a lower resolution"
    )
  }

  list(
    members = members, types = rep("cluster", length(members)),
    log_message = paste0(
      length(members), " Leiden module(s) in a ", igraph::vcount(g),
      "-node component (objective = modularity, resolution = ", resolution,
      ", n_iterations = ", n_iterations, ")"
    ),
    warn_message = warn_message, qb = NA_real_
  )
}

#' Bipartite (Barber Q_B) community detection on one connected component
#' @return `list(members, types, log_message, warn_message, qb)`.
#' @keywords internal
.network_detect_modules_bipartite <- function(g, seed) {
  node_names <- igraph::V(g)$name
  mat <- igraph::as_biadjacency_matrix(g, sparse = FALSE)

  ## bipartite::computeModules() errors on any biadjacency matrix with a
  ## dimension < 2 (a star / single-mode component), not only a zero
  ## dimension -- the intermediate web degrades to a vector and the
  ## moduleWeb S4 slot assignment refuses it. A star has no bipartite
  ## module structure anyway, so fall back to a single module and log.
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    return(list(
      members = stats::setNames(list(node_names), "M1"), types = "cluster",
      log_message = paste0(
        "component with ", nrow(mat), " compound(s) and ", ncol(mat),
        " target(s) (a star / single-mode component); treated as a single module ",
        "instead of running bipartite modularity"
      ),
      warn_message = NULL, qb = NA_real_
    ))
  }

  bm <- .network_bipartite_membership(mat, seed)
  memb <- bm$membership
  dropped <- names(memb)[is.na(memb)]
  ok <- memb[!is.na(memb)]
  ids <- sort(unique(ok))
  members <- lapply(ids, function(k) names(ok)[ok == k])
  names(members) <- paste0("M", ids)
  types <- rep("cluster", length(members))
  ## Within one connected component every compound/target has >= 1 edge, so
  ## bipartite::empty() drops nothing in practice -- but if it ever does,
  ## bucket the dropped nodes into a residual module so the
  ## sum(n_nodes) == vcount(g) invariant still holds.
  if (length(dropped) > 0) {
    members <- c(members, stats::setNames(list(dropped), "Miso"))
    types <- c(types, "isolated")
  }

  list(
    members = members, types = types,
    log_message = paste0(
      length(ids), " bipartite module(s) (Beckett) in a ", igraph::vcount(g),
      "-node component; Q_B = ", signif(bm$qb, 4),
      if (length(dropped) > 0) paste0("; ", length(dropped), " isolated node(s) in a residual module") else ""
    ),
    warn_message = NULL, qb = bm$qb
  )
}

#' Remap `bipartite::computeModules()` membership back to node names
#'
#' @description
#' `computeModules()` runs `bipartite::empty()` internally (drops zero-sum
#' rows/cols), so `module2constraints()` is indexed against the **emptied**
#' matrix. The remap therefore uses `mod@originalWeb`'s dimnames, not
#' `c(rownames(mat), colnames(mat))` of the matrix passed in -- getting this
#' wrong shifts every assignment past the first dropped row. Nodes dropped
#' by `empty()` get `module_id = NA`.
#' @return A list: `membership`, a named integer vector over
#'   `c(rownames(mat), colnames(mat))` with `NA` for `empty()`-dropped
#'   nodes; `qb`, Barber's Q_B (`mod@likelihood`); and `degenerate`, `TRUE`
#'   when `computeModules()` found no usable partition (all nodes then in
#'   one module).
#' @keywords internal
.network_bipartite_membership <- function(mat, seed) {
  restore_seed <- .with_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  mod <- suppressWarnings(bipartite::computeModules(mat, method = "Beckett"))
  memb <- bipartite::module2constraints(mod)
  ow <- mod@originalWeb
  all_names <- c(rownames(mat), colnames(mat))
  present <- c(rownames(ow), colnames(ow))
  out <- stats::setNames(rep(NA_integer_, length(all_names)), all_names)

  ## module2constraints() returns a *list* (not an integer vector) when a
  ## species matches zero or several modules -- i.e. there is no usable
  ## bipartite partition (Q_B ~ 0, e.g. a complete bipartite component).
  ## Treat that as a single module over the emptied web.
  if (!is.numeric(memb) || length(memb) != length(present)) {
    out[present] <- 1L
    return(list(membership = out, qb = mod@likelihood, degenerate = TRUE))
  }

  names(memb) <- present
  out[present] <- as.integer(memb)
  list(membership = out, qb = mod@likelihood, degenerate = FALSE)
}

#' HDBSCAN module detection on one connected component (over shortest-path
#' distances, always finite within a component), collecting HDBSCAN "noise"
#' nodes into their own module rather than dropping them
#' @return `list(members, types, log_message, warn_message = NULL, qb =
#'   NA_real_)`.
#' @keywords internal
.network_detect_modules_hdbscan <- function(g, min_module_size) {
  node_names <- igraph::V(g)$name
  n_nodes <- length(node_names)

  ## weights = NA: the edge `weight` attribute is the import probability,
  ## not a distance -- letting igraph use it as one would make confident
  ## predictions look "far" (and errors outright on any NA probability).
  dist_matrix <- igraph::distances(g, weights = NA)
  hdb <- dbscan::hdbscan(stats::as.dist(dist_matrix), minPts = min_module_size)

  cluster_ids <- sort(unique(hdb$cluster[hdb$cluster != 0]))
  noise_nodes <- node_names[hdb$cluster == 0]
  n_noise <- length(noise_nodes)

  if (length(cluster_ids) == 0) {
    return(list(
      members = stats::setNames(list(node_names), "M1"), types = "noise",
      log_message = paste0(
        "HDBSCAN found no clusters in a ", n_nodes, "-node component (all noise); treated as a single noise module"
      ),
      warn_message = NULL, qb = NA_real_
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
    members = members, types = types,
    log_message = paste0(
      length(cluster_ids), " module(s) found via HDBSCAN in a ", n_nodes, "-node component (minPts = ",
      min_module_size, "); ", n_noise, " node(s) placed in a separate noise module (not dropped)"
    ),
    warn_message = NULL, qb = NA_real_
  )
}

#' Targeted-attack percolation on one (sub)graph: repeatedly remove the
#' current highest-degree node, track the largest-component fraction
#'
#' @description
#' `r_index` is Schneider et al. (2011), *PNAS* 108(10), 3838-3841's
#' robustness index: the mean of `s(Q)` (largest-component fraction) over
#' removal steps `Q = 1..N` -- **not** including the trivial, pre-removal
#' `s(0) = 1` state. `curve[1] == 1` is kept in the returned `curve` for
#' plotting but must be dropped before averaging.
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

  list(curve = curve, r_index = mean(curve[-1]))
}

#' Random-failure percolation on one (sub)graph: remove nodes in a uniformly
#' random order (no degree recomputation), averaged over `n_random`
#' independent replicates
#'
#' @description
#' Schneider et al. (2011) always contrast targeted attack against random
#' failure. A single random draw has enough variance to invert the
#' comparison for a non-trivial fraction of modules, so the stored curve is
#' the mean over `n_random` replicates and `curve_sd` its per-step SD. Each
#' replicate `r` draws its order under `.with_seed(.derive_seed(seed, r))`
#' (the caller passes a per-module-distinct `seed`, so replicates are also
#' distinct across modules; `.derive_seed()` keeps the arithmetic clear of
#' integer overflow). `r_index` is `mean(mean_curve[-1])`, which by
#' linearity equals the mean over replicates of each replicate's Schneider
#' R-index.
#' @return `list(curve, curve_sd, r_index, n_replicates)`.
#' @keywords internal
.network_percolate_random <- function(g, seed, n_random) {
  n <- igraph::vcount(g)
  if (n == 1) {
    return(list(curve = c(1, 0), curve_sd = c(0, 0), r_index = 0, n_replicates = n_random))
  }

  node_names <- igraph::V(g)$name
  rep_curves <- matrix(NA_real_, nrow = n_random, ncol = n + 1)

  for (r in seq_len(n_random)) {
    draw_seed <- .with_seed(.derive_seed(seed, r))
    removal_order <- sample(node_names)
    draw_seed()

    remaining <- g
    curve <- numeric(n + 1)
    curve[1] <- 1
    for (step in seq_len(n)) {
      remaining <- igraph::delete_vertices(remaining, removal_order[step])
      if (igraph::vcount(remaining) == 0) {
        largest_cc <- 0
      } else {
        largest_cc <- max(igraph::components(remaining)$csize)
      }
      curve[step + 1] <- largest_cc / n
    }
    rep_curves[r, ] <- curve
  }

  mean_curve <- colMeans(rep_curves)
  sd_curve <- apply(rep_curves, 2, stats::sd)
  list(curve = mean_curve, curve_sd = sd_curve, r_index = mean(mean_curve[-1]), n_replicates = n_random)
}

#' Combine a base seed with an integer offset without integer overflow
#'
#' @description
#' `base + offset` is plain integer arithmetic in R and overflows to `NA`
#' (then `set.seed(NA)` errors) when `base` is near `.Machine$integer.max`
#' -- which `sample.int(.Machine$integer.max, 1)` can draw for the default
#' seed. The addition is done in double and wrapped back into valid
#' integer range, so a given top-level `seed` still maps deterministically
#' to the same per-replicate sub-seeds.
#' @keywords internal
.derive_seed <- function(base, offset) {
  as.integer((as.double(base) + as.double(offset)) %% .Machine$integer.max)
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
  data.frame(
    condition = character(0), module_id = character(0), module_type = character(0),
    n_nodes = integer(0), r_index = double(0), r_index_random = double(0),
    n_random = integer(0), clustering = character(0), modularity = double(0),
    resolution = double(0), method_detail = character(0), seed_used = integer(0),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_network_robustness_curve_row <- function() {
  data.frame(
    condition = character(0), module_id = character(0), n_removed = integer(0),
    largest_component_fraction = double(0), removal_strategy = character(0),
    largest_component_sd = double(0), n_replicates = integer(0),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_network_module_membership_row <- function() {
  data.frame(
    condition = character(0), node_id = character(0), node_type = character(0),
    module_id = character(0), module_type = character(0),
    stringsAsFactors = FALSE
  )
}
