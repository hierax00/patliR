#' @include AllGenerics.R internal.R network_build.R
NULL

## Shared boilerplate for the plot_*() family. Every network plot_*
## function repeats the same four blocks: the ggplot2/engine availability
## check, the "which conditions / pooled-vs-per-condition label" logic,
## the ggsave + log-upsert + girafe-wrap + attr(result, "proj") tail, and
## the node-label lookup. Those live here so a change to any one of them
## happens in a single place.

#' Check plotting dependencies and resolve the drawing engine
#'
#' @description
#' The `ggplot2` availability check plus the `engine = "ggiraph"` ->
#' `"static"` fallback that every `plot_*` function opens with. The name of
#' the calling function is picked up from the call stack so the error
#' message still reads `plot_<name>()`.
#'
#' @param engine `NULL` (the function has no `engine` argument, e.g.
#'   `plot_bowtie()`) or a single already-`match.arg()`-ed engine string.
#'   When `"ggiraph"` and `ggiraph` is not installed, a warning is emitted
#'   and `"static"` is returned instead.
#' @param extra Character vector of additional package names the plot needs
#'   (e.g. `"ggalluvial"`); each is checked with the same message shape as
#'   the `ggplot2` check.
#' @return The resolved engine (invisibly unchanged when `engine` is `NULL`
#'   or already `"static"`; possibly downgraded from `"ggiraph"`).
#' @keywords internal
.plot_require <- function(engine = NULL, extra = character()) {
  fn <- tryCatch(utils::tail(as.character(sys.call(-1L)[[1L]]), 1L), error = function(e) "this function")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("The {.pkg ggplot2} package is required for {.fn {fn}}.")
  }
  for (pkg in extra) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(c(
        "{.fn {fn}} needs the {.pkg {pkg}} package (CRAN), not installed.",
        "i" = "{.code install.packages(\"{pkg}\")} -- pure CRAN, ggplot2-based."
      ))
    }
  }
  if (!is.null(engine) && length(engine) == 1L && !is.na(engine) && engine == "ggiraph" &&
      !requireNamespace("ggiraph", quietly = TRUE)) {
    cli::cli_warn("The {.pkg ggiraph} package is not installed; falling back to {.val static}.")
    engine <- "static"
  }
  engine
}

#' Resolve a `plot_*` `condition` argument to a condition vector and a
#' scope label
#'
#' @description
#' `.network_resolve_conditions()` plus the `"ALL"` (pooled) vs.
#' `"cond1+cond2"` (explicit) label string that the `plot_*` family uses
#' for filenames, plot titles and the `condition` column of every plot log.
#'
#' @param proj A `PatliRProject`.
#' @param condition `NULL` (pool every built condition) or a character
#'   vector of condition names.
#' @return `list(conditions = <character>, scope_label = <character(1)>)`.
#' @keywords internal
.plot_scope <- function(proj, condition) {
  conditions <- .network_resolve_conditions(proj, condition)
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")
  list(conditions = conditions, scope_label = scope_label)
}

#' Save a `plot_*` figure, upsert its log row, and attach the updated
#' project
#'
#' @description
#' The tail every `plot_*` function shares: when `save = TRUE`, resolve
#' `out_dir`, `ggsave()` the figure, upsert a one-row log into the
#' `name` results slot (keyed by `key_cols`), write that slot's CSV; then
#' wrap the plot in `girafe()` for `engine = "ggiraph"` and stamp
#' `attr(result, "proj")` with the updated project.
#'
#' @param proj A `PatliRProject`.
#' @param save Logical; when `FALSE` nothing is written and `p` (or its
#'   `girafe` wrapper) is returned without the `proj` attribute.
#' @param out_dir `NULL` (default -> `file.path(projectDir(proj), "plots")`)
#'   or an explicit directory.
#' @param width,height,dpi Passed to [ggplot2::ggsave()].
#' @param p The plot object to return (a `ggplot`, or already a `girafe`
#'   is never passed -- the wrap happens here).
#' @param name Results-bag slot name for the log, e.g.
#'   `"centrality_plot_log"`.
#' @param filename Basename of the file to write inside `out_dir`.
#' @param log_row One-row `data.frame`. If it has a `path` column that
#'   value is overwritten with the resolved path (column position kept);
#'   otherwise a `path` column is appended.
#' @param key_cols Passed straight to [.network_upsert()].
#' @param engine `NULL` / `"static"` -> return `p`; `"ggiraph"` -> return
#'   `girafe(ggobj = p)`.
#' @param static Plot object actually handed to `ggsave()` -- defaults to
#'   `p`, but [plot_network_layers()] / [plot_network_degeneracy()] pass a
#'   separately rendered static plot because their interactive `p` does not
#'   rasterise well.
#' @return `p` (or its `girafe` wrapper), with `attr(., "proj")` set when
#'   `save = TRUE`.
#' @keywords internal
.plot_finish <- function(proj, p, name, filename, log_row, key_cols,
                         engine = NULL, save = TRUE, out_dir = NULL,
                         width = 8, height = 6, dpi = 150, static = p) {
  if (save) {
    if (is.null(out_dir)) out_dir <- file.path(projectDir(proj), "plots")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    path <- file.path(out_dir, filename)
    ggplot2::ggsave(path, static, width = width, height = height, dpi = dpi)
    log_row$path <- path
    log_df <- .network_upsert(proj, name, log_row, key_cols)
    patliRResults(proj, name) <- log_df
    .write_results_csv(proj, name, log_df)
  }
  result <- if (is.null(engine) || engine == "static") {
    p
  } else {
    ggiraph::girafe(ggobj = p, options = list(ggiraph::opts_tooltip(opacity = 0.9)))
  }
  if (save) attr(result, "proj") <- proj
  result
}

#' Look up display labels for a set of network nodes
#'
#' @description
#' Formerly `.network_layers_labels()` (it lived in `plot_network_layers.R`,
#' which four other plot files `@include`d purely to reach it). Resolves
#' compound IDs to `compounds()$name`, UniProt IDs to gene symbols (via
#' `clusterProfiler::bitr()` + `org.Hs.eg.db`, when both are installed),
#' and pathway IDs to `network_enrichment$Description`; anything that does
#' not resolve is left as its raw ID.
#'
#' @param proj A `PatliRProject`.
#' @param conditions Character vector of conditions in scope (used to
#'   restrict the pathway-description lookup).
#' @param ids Character vector of node IDs.
#' @param layer Character vector, recycled against `ids`: `"compound"`,
#'   `"target"`, `"pathway"` or `"disease"` per node.
#' @return Character vector, same length/order as `ids`.
#' @keywords internal
.plot_label_nodes <- function(proj, conditions, ids, layer) {
  nodes <- data.frame(name = ids, layer = layer, stringsAsFactors = FALSE)
  label <- nodes$name

  is_compound <- nodes$layer == "compound"
  if (any(is_compound)) {
    cmp <- compounds(proj)
    lookup <- stats::setNames(cmp$name, cmp$id)
    resolved <- lookup[nodes$name[is_compound]]
    label[is_compound] <- ifelse(is.na(resolved) | resolved == "", nodes$name[is_compound], resolved)
  }

  is_target <- nodes$layer == "target"
  if (any(is_target) && requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    map <- tryCatch(
      clusterProfiler::bitr(nodes$name[is_target], fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db", drop = TRUE),
      error = function(e) NULL
    )
    if (!is.null(map) && nrow(map) > 0) {
      lookup <- stats::setNames(map$SYMBOL, map$UNIPROT)
      resolved <- lookup[nodes$name[is_target]]
      label[is_target] <- ifelse(is.na(resolved), nodes$name[is_target], resolved)
    }
  }

  is_pathway <- nodes$layer == "pathway"
  if (any(is_pathway)) {
    enrichment_all <- patliRResults(proj, "network_enrichment")
    if (!is.null(enrichment_all)) {
      enr <- enrichment_all[enrichment_all$condition %in% conditions, , drop = FALSE]
      lookup <- stats::setNames(enr$Description, enr$ID)
      resolved <- lookup[nodes$name[is_pathway]]
      label[is_pathway] <- ifelse(is.na(resolved) | resolved == "", nodes$name[is_pathway], resolved)
    }
  }

  label
}
