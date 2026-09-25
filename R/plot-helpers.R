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
#' `"cond1+cond2_<hash>"` (explicit multiple) label string that the `plot_*` family uses
#' for filenames, plot titles and the `condition` column of every plot log.
#' Single-condition labels stay unchanged except literal `"ALL"`, which
#' also gets a hash to distinguish it from the pooled default.
#'
#' @param proj A `PatliRProject`.
#' @param condition `NULL` (pool every built condition) or a character
#'   vector of condition names.
#' @return `list(conditions = <character>, scope_label = <character(1)>)`.
#' @keywords internal
.plot_scope <- function(proj, condition) {
  conditions <- .network_resolve_conditions(proj, condition)
  if (length(conditions) == 0) cli::cli_abort("{.arg condition} must select at least one condition.")
  scope_label <- if (is.null(condition)) "ALL" else paste(conditions, collapse = "+")
  if (!is.null(condition) && (length(conditions) > 1 || identical(conditions, "ALL"))) {
    conditions <- sort(unique(conditions))
    scope_label <- paste0(paste(conditions, collapse = "+"), "_", rlang::hash(conditions))
  }
  list(conditions = conditions, scope_label = scope_label)
}

#' Disambiguate display labels without changing node identity
#' @keywords internal
.plot_unique_labels <- function(proj, conditions, ids, type) {
  unique_ids <- unique(ids)
  labels <- .plot_label_nodes(proj, conditions, unique_ids, type)
  duplicate <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate] <- paste0(labels[duplicate], " (", unique_ids[duplicate], ")")
  labels <- make.unique(labels)
  labels[match(ids, unique_ids)]
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
    ## bg = "white": theme_void() figures have a blank plot.background, and
    ## without an explicit device background the PNG comes out transparent
    ## (rendered black by most image viewers).
    ggplot2::ggsave(path, static, width = width, height = height, dpi = dpi, bg = "white")
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

#' Back-fill a new key column into an existing plot log
#'
#' @description
#' When a `plot_*` log gains a key column (e.g. `subset`, `view`), logs
#' written by an earlier version lack it and [.network_upsert()] could not
#' index them. The column is added with the value the old behaviour
#' corresponds to, so the default call still replaces its own old row.
#' @return `proj`, possibly with the `name` results slot updated.
#' @keywords internal
.plot_log_backfill <- function(proj, name, col, value) {
  log_df <- patliRResults(proj, name)
  if (!is.null(log_df) && nrow(log_df) > 0 && !col %in% names(log_df)) {
    log_df[[col]] <- value
    patliRResults(proj, name) <- log_df
  }
  proj
}

#' Wrap title / subtitle / caption text to a maximum line length
#'
#' @description
#' ggplot2 never wraps `labs()` text, so a long subtitle runs off the right
#' edge of the saved PNG. Each element of `x` is wrapped with
#' [strwrap()] at `width` characters; explicit `"\n"` line breaks already
#' in the text are kept (each line is wrapped on its own).
#'
#' @param x Character vector (`NULL`/`NA` pass through unchanged).
#' @param width Maximum characters per line -- typically
#'   `.plot_wrap_width()` of the figure's width.
#' @return Character vector, same length as `x`.
#' @keywords internal
.plot_wrap <- function(x, width = 100) {
  if (is.null(x)) return(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
    if (length(lines) == 0) return(s)
    paste(vapply(lines, function(l) {
      if (!nzchar(l)) return("")
      paste(unlist(lapply(strwrap(l, width = width), .plot_break_long, width = width)), collapse = "\n")
    }, character(1)), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

#' Break one over-long, space-free line (a systematic chemical name)
#' @description
#' [strwrap()] only breaks at spaces, so
#' `"2-(4-Ethenyl-4-methyl-3-(prop-1-en-2-yl)cyclohexyl)propan-2-ol"` stays
#' one line. Such a line is cut after the last `-`, `,` or `)` that keeps
#' the piece within `width` (hard cut at `width` when there is none).
#' @return Character vector of pieces, each at most `width` characters.
#' @keywords internal
.plot_break_long <- function(line, width) {
  out <- character(0)
  while (nchar(line) > width) {
    head <- substr(line, 1, width)
    cut <- max(c(0L, gregexpr("[-,)]", head)[[1]]))
    if (cut < width %/% 3) cut <- width
    out <- c(out, substr(line, 1, cut))
    line <- sub("^ +", "", substr(line, cut + 1, nchar(line)))
  }
  c(out, line)
}

#' Characters per line that fit across a saved figure
#'
#' @description
#' Rough capacity of one text line spanning the whole figure: an average
#' glyph of a sans-serif font is ~0.55 em wide, so `size` points of text
#' fit `72 / (0.55 * size)` characters per inch. `margin` (inches) covers
#' the plot margins. Deliberately a little conservative -- an early line
#' break costs nothing, a clipped subtitle loses information.
#'
#' @param fig_width Figure width in inches (the `width` passed to `ggsave()`).
#' @param size Font size in points of the text being wrapped.
#' @param margin Inches of the figure width not available to the text.
#' @return Integer, at least 20.
#' @keywords internal
.plot_wrap_width <- function(fig_width, size = 8, margin = 0.3) {
  max(20L, as.integer(floor((fig_width - margin) * 72 / (0.55 * size))))
}

#' Shorten long labels for on-plot text
#'
#' @description
#' Compound names from GC-MS libraries can run to 60+ characters
#' (`"Naphthalene, 1,2-dihydro-1,1,6-trimethyl-"`); drawn as point labels they
#' collide and run off the panel. Labels longer than `max_chars` are cut
#' and end in `"..."`. Only for drawn text -- tooltips and logs keep the
#' full name.
#'
#' @param x Character vector.
#' @param max_chars Maximum label length, including the `"..."`.
#' @return Character vector, same length as `x`.
#' @keywords internal
.plot_truncate <- function(x, max_chars = 28) {
  x <- as.character(x)
  long <- !is.na(x) & nchar(x) > max_chars
  x[long] <- paste0(substr(x[long], 1, max_chars - 3), "...")
  x
}

#' Is label repulsion available?
#'
#' @description
#' `ggrepel` is optional (Suggests). Label layers use
#' [ggrepel::geom_text_repel()] when it is installed and fall back to plain
#' [ggplot2::geom_text()] otherwise. `options(patliR.repel = FALSE)` forces
#' the fallback (reproducible figures across machines, or tests).
#' @return `TRUE`/`FALSE`.
#' @keywords internal
.plot_use_repel <- function() {
  isTRUE(getOption("patliR.repel", TRUE)) && requireNamespace("ggrepel", quietly = TRUE)
}

#' A text-label layer that repels when `ggrepel` is installed
#'
#' @description
#' Arguments in `...` go to either geom; `repel_args` only to
#' [ggrepel::geom_text_repel()] (e.g. `max.overlaps`, `box.padding`) and
#' `text_args` only to the [ggplot2::geom_text()] fallback (e.g. a `vjust`
#' offset that would fight the repulsion).
#' @return A ggplot2 layer.
#' @keywords internal
.plot_text_layer <- function(..., repel_args = list(), text_args = list()) {
  if (.plot_use_repel()) {
    do.call(ggrepel::geom_text_repel, c(list(...), repel_args))
  } else {
    do.call(ggplot2::geom_text, c(list(...), text_args))
  }
}

#' Resolve disease IDs to their human-readable names
#'
#' @description
#' Looks the IDs up in `disease_genes$disease_name` (from
#' [disease_genes_fetch()]), then `targets_disease_profile$disease_name`,
#' then `targets_disease$disease_name` when that column exists -- the first
#' non-empty name wins.
#' @param proj A `PatliRProject`.
#' @param ids Character vector of disease IDs (e.g. `"MONDO_0005044"`).
#' @return Named character vector (names = `ids`), `NA` where no name is
#'   recorded.
#' @keywords internal
.plot_disease_names <- function(proj, ids) {
  ids <- as.character(ids)
  out <- stats::setNames(rep(NA_character_, length(ids)), ids)
  for (slot in c("disease_genes", "targets_disease_profile", "targets_disease")) {
    todo <- is.na(out)
    if (!any(todo)) break
    tab <- tryCatch(patliRResults(proj, slot), error = function(e) NULL)
    if (is.null(tab) || nrow(tab) == 0 || !all(c("disease_id", "disease_name") %in% names(tab))) next
    tab <- tab[!is.na(tab$disease_name) & nzchar(tab$disease_name), c("disease_id", "disease_name"), drop = FALSE]
    lookup <- stats::setNames(as.character(tab$disease_name), tab$disease_id)
    out[todo] <- unname(lookup[ids[todo]])
  }
  out
}

#' Display label for a disease: `"name (ID)"`, the ID alone when no name is
#' recorded, the name alone when it already contains the ID
#' @param proj A `PatliRProject`.
#' @param ids Character vector of disease IDs.
#' @param with_id Logical; `FALSE` drops the `" (ID)"` suffix.
#' @return Character vector, same length as `ids`.
#' @keywords internal
.plot_disease_label <- function(proj, ids, with_id = TRUE) {
  ids <- as.character(ids)
  nm <- unname(.plot_disease_names(proj, ids))
  has_name <- !is.na(nm) & nzchar(nm)
  ## "inflammatory response (GO:0006954)" already carries its ID, in either
  ## the GO:0006954 or GO_0006954 spelling
  contains_id <- has_name & (mapply(grepl, ids, nm, fixed = TRUE) |
                               mapply(grepl, sub("_", ":", ids, fixed = TRUE), nm, fixed = TRUE))
  out <- ids
  out[has_name] <- if (with_id) {
    ifelse(contains_id[has_name], nm[has_name], paste0(nm[has_name], " (", ids[has_name], ")"))
  } else {
    nm[has_name]
  }
  out
}

#' Qualitative palette in which neighbouring colours never look alike
#'
#' @description
#' For figures that place categories side by side (e.g. the term arcs and
#' ribbons of [plot_gochord()]), where a hue ramp such as
#' [grDevices::rainbow()] makes adjacent categories nearly identical.
#' \describe{
#'   \item{`"contrast"`}{12 colours from the Okabe-Ito and Tol
#'     colour-blind-safe sets, ordered so consecutive colours differ as much
#'     as possible in lightness and hue (smallest CIELAB distance between
#'     neighbours > 80). Beyond 12, the cycle repeats darkened, then
#'     lightened, so every colour stays distinct.}
#'   \item{`"grey"`}{`n` greys (for black-and-white print) spread from near
#'     black to near white and interleaved dark, light, dark, light, so two
#'     neighbours are always about half the grey range apart.}
#'   \item{`"default"`}{[grDevices::rainbow()] -- the previous behaviour.}
#' }
#' @param n Number of colours.
#' @param palette `"contrast"`, `"grey"` or `"default"`.
#' @return Character vector of `n` distinct hex colours.
#' @keywords internal
.plot_contrast_palette <- function(n, palette = c("contrast", "grey", "default")) {
  palette <- match.arg(palette)
  n <- as.integer(n)
  if (is.na(n) || n <= 0) return(character(0))
  if (palette == "default") return(grDevices::rainbow(n))
  if (palette == "grey") {
    levels <- seq(0.12, 0.95, length.out = max(n, 2))[seq_len(n)]
    if (n == 1) levels <- 0.35
    half <- ceiling(n / 2)
    ord <- as.vector(rbind(seq_len(half), half + seq_len(half)))
    ord <- ord[ord <= n]
    return(grDevices::gray(levels[ord]))
  }
  ## order chosen to maximise the smallest CIELAB distance between
  ## consecutive colours (every neighbour pair differs by dE > 80)
  base <- c("#0072B2", "#D55E00", "#56B4E9", "#F0E442", "#882255", "#E69F00",
            "#CC79A7", "#117733", "#EE8866", "#009E73", "#332288", "#BBBBBB")
  mix <- function(cols, target, w) {
    m <- grDevices::col2rgb(cols) / 255
    t <- grDevices::col2rgb(target) / 255
    m <- m * (1 - w) + as.vector(t) * w
    grDevices::rgb(m[1, ], m[2, ], m[3, ])
  }
  cycles <- ceiling(n / length(base))
  out <- unlist(lapply(seq_len(cycles), function(k) {
    if (k == 1) return(base)
    ## alternate darker / lighter shades of the base cycle, further each round
    w <- min(0.2 + 0.15 * ((k - 2) %/% 2), 0.7)
    if (k %% 2 == 0) mix(base, "black", w) else mix(base, "white", w)
  }))
  out <- toupper(out[seq_len(n)])
  ## guarantee distinct values (GOplot's legend uses unique() of the colours)
  dup <- duplicated(out)
  while (any(dup)) {
    out[dup] <- toupper(mix(out[dup], "grey50", 0.1))
    dup <- duplicated(out)
  }
  out
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
