#' @include AllGenerics.R internal.R network_build.R network_proximity.R
NULL

## A bowtie decomposition needs a directed graph with cycles. The layered
## graph (network_layers.R) is a DAG (no cycles); STRING's regular graph
## has cycles but is undirected. So this uses STRING's directed "actions"
## channel (protein.actions.v<version>...txt.gz), which no STRINGdb R
## method exposes -- it is downloaded and parsed directly, cached under
## cacheDir(proj)/stringdb/. Only mechanistically-characterized
## interactions get a direction, so this is a much sparser graph than the
## regular "links" file.
##
## `actions_version` defaults to "11.0" (the actions flat file was
## discontinued after v11.0). The UniProt -> STRING_id mapping is done
## against that SAME release, not a newer one -- STRING re-derives its ENSP
## identifier space every release, so mixing versions would silently drop
## retired IDs into the "not in the actions network" bucket.
##
## The decomposition is a whole-network property, so it is computed once
## per (species, actions_version) and cached; each condition's targets are
## then annotated with which component they fall into.

#' Bowtie architecture of the STRING directed-action network, annotated per
#' condition's targets
#'
#' @description
#' Downloads STRING's directed "actions" file for `species` (activation/
#' inhibition/... interactions with known directionality -- a much sparser,
#' but genuinely *directed*, subset of the interactome `network_proximity()`
#' uses), decomposes it into the classic bowtie structure (Broder, A. et al.
#' 2000, "Graph structure in the Web," *Computer Networks* 33(1-6), 309-320),
#' and reports which component each of a condition's mapped targets falls
#' into:
#'
#' - **core**: the largest strongly connected component (every node can
#'   reach every other node by a directed path).
#' - **in_component**: nodes that can reach the core but cannot be reached
#'   from it (upstream regulators of the core).
#' - **out_component**: nodes reachable from the core but that cannot reach
#'   it back (downstream effectors of the core).
#' - **other**: everything else (disconnected, or connected only via
#'   "tendrils"/"tubes" that neither reach nor are reached by the core --
#'   Broder et al.'s full taxonomy distinguishes these further; collapsed
#'   here into one category, consistent with how this simplified 4-way
#'   split is normally reported in network-biology papers).
#'
#' @section This is a species-wide structural landmark, not a per-condition subgraph:
#' Unlike most of the `network_*` family, the bowtie decomposition itself
#' does not depend on `condition` at all -- it is a property of the whole
#' directed action network for `species`, computed once and cached (see
#' `.network_stringdb_actions_graph()`, internal). `condition` only
#' controls which compounds/targets get annotated with their (fixed)
#' component membership.
#'
#' @section `species`/`version` are documented once, not inherited:
#' `@inheritParams network_proximity` would also pull in that function's
#' `version` doc, which describes it as the live STRING release passed to
#' `STRINGdb$new()` -- true there, but `version` is deprecated and ignored
#' *here* (see its own `@param` below). Only `species` is actually shared;
#' it is documented explicitly instead, so the two functions' `version`
#' semantics never end up contradicting each other on the same page.
#'
#' @param proj,condition See [network_build()].
#' @param species NCBI taxonomy ID for the STRING actions file and the
#'   UniProt -> STRING_id mapping. Default `9606` (human).
#' @param actions_version STRING release for the directed "actions" file
#'   *and* the UniProt -> STRING_id mapping (see
#'   `.network_stringdb_actions_graph()`, internal). Defaults to `"11.0"`:
#'   STRING discontinued the "actions" flat file after v11.0.
#' @param actions_score_threshold Minimum STRING per-interaction confidence
#'   `score` (0-999, the actions file's own column -- distinct from the
#'   combined interactome score [network_proximity()] filters on) required
#'   for a directed action to become an edge. Default `400`, matching the
#'   rest of the package's STRING confidence convention. Without this
#'   filter every directed row is used regardless of confidence, so the
#'   core SCC size -- the headline bowtie number -- would be driven by
#'   STRING's lowest-confidence predictions.
#' @param version Deprecated and unused -- both the actions graph and the
#'   ID mapping now use `actions_version` (mixing releases dropped retired
#'   ENSP IDs into the "not in the actions network" bucket). Kept so old
#'   calls do not error.
#'
#' @return The updated `proj`, with two entries in [patliRResults()]:
#'   `network_bowtie` (columns `condition`, `compound_id`, `uniprot_id`,
#'   `string_id`, `bowtie_component` -- one row per (condition, compound,
#'   target). `bowtie_component` is `"core"` / `"in_component"` /
#'   `"out_component"` / `"other"` (tendrils, tubes, disconnected) for a
#'   target in the action graph, `"not_in_action_network"` for a target
#'   that maps to STRING but is absent from the sparse directed action
#'   graph, or `"unmapped"` for a target with no STRING_id at all -- none
#'   are dropped) and
#'   `network_bowtie_summary` (columns `species`, `version`,
#'   `actions_score_threshold`, `n_nodes`, `n_edges`, `n_core`,
#'   `n_in_component`, `n_out_component`, `n_other` -- one row per
#'   `(species, version)`; `version` here is `actions_version`).
#'   Both written to their matching `results/*.csv`. A core smaller than 50
#'   nodes, or under 1% of `n_nodes`, gets a `cli_warn` -- every downstream
#'   bowtie number is likely noise at that size.
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
#' proj <- network_bowtie(proj, condition = "FLO-ET") # needs STRINGdb + internet on first call
#' patliRResults(proj, "network_bowtie")
#' patliRResults(proj, "network_bowtie_summary")
#' }
#'
#' @export
network_bowtie <- function(proj, condition = NULL, species = 9606, version = "12.0",
                            actions_version = "11.0", actions_score_threshold = 400) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.numeric(actions_score_threshold), length(actions_score_threshold) == 1,
            !is.na(actions_score_threshold), actions_score_threshold >= 0, actions_score_threshold <= 999)
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_bowtie} needs {.pkg STRINGdb} (for UniProt -> STRING_id mapping), not installed.",
      "i" = "Install it with {.code BiocManager::install(\"STRINGdb\")}."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)

  actions <- .network_stringdb_actions_graph(proj, species, actions_version, score_threshold = actions_score_threshold)
  g <- actions$graph
  bowtie <- .network_bowtie_classify(g)

  summary_row <- data.frame(
    species = species, version = actions_version, actions_score_threshold = actions_score_threshold,
    n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g),
    n_core = sum(bowtie == "core"), n_in_component = sum(bowtie == "in_component"),
    n_out_component = sum(bowtie == "out_component"), n_other = sum(bowtie == "other"),
    stringsAsFactors = FALSE
  )

  ## A core this small (either in absolute terms, or relative to the whole
  ## actions network) makes every downstream bowtie number noise -- the
  ## in/out/other split is only meaningful relative to a real core.
  if (summary_row$n_core < 50 || summary_row$n_core < 0.01 * summary_row$n_nodes) {
    cli::cli_warn(c(
      "The bowtie core for species {.val {species}}/v{actions_version} (actions_score_threshold = {actions_score_threshold}) has only {summary_row$n_core} node(s) out of {summary_row$n_nodes}.",
      "i" = "Every downstream bowtie number (in/out/other component sizes, per-target classification) is likely noise at this size -- consider a lower {.arg actions_score_threshold}, or treat this run as inconclusive."
    ))
  }

  ## UniProt -> STRING_id mapping, read directly from STRING's own alias
  ## flat file (protein.aliases.v<actions_version>) rather than
  ## instantiating a second STRINGdb object at score_threshold = 0 (which
  ## would download the FULL interactome purely to reach an alias table
  ## that never depended on score in the first place) -- same pattern as
  ## the actions file itself: download once, cache, parse directly.
  ## Against the SAME STRING release as the actions graph (`actions_version`,
  ## not the pipeline-wide `version`): STRING re-derives its ENSP identifier
  ## space each release, so a v12.0 mapping checked against a v11.0 actions
  ## graph would silently miss retired IDs.
  aliases <- .network_stringdb_aliases_map(proj, species, actions_version)
  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    if (nrow(ct) == 0) {
      rows[[cond]] <- .empty_network_bowtie_row()
      next
    }

    uniprot_ids_needed <- unique(ct$uniprot_id)
    matched <- aliases[match(uniprot_ids_needed, aliases$alias), , drop = FALSE]
    uni_to_string <- stats::setNames(matched$string_protein_id, uniprot_ids_needed)
    unmapped <- uniprot_ids_needed[is.na(uni_to_string)]
    if (length(unmapped) > 0) {
      proj <- .log_append(
        proj, step = "network_bowtie", id = unmapped,
        message = "network_bowtie_unmapped: no STRING_id found for this UniProt ID in this species; component = 'unmapped'"
      )
    }

    string_id <- uni_to_string[ct$uniprot_id]
    component <- ifelse(
      is.na(string_id), "unmapped",
      ## a mapped target that is absent from the (much sparser) directed
      ## actions graph is NOT a bow-tie finding -- keep it distinct from a
      ## genuine tendril/tube node, which the classifier already folds into
      ## "other".
      ifelse(string_id %in% names(bowtie), bowtie[string_id], "not_in_action_network")
    )

    rows[[cond]] <- data.frame(
      condition = cond, compound_id = ct$compound_id, uniprot_id = ct$uniprot_id,
      string_id = unname(string_id), bowtie_component = component, stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  bowtie_touched <- unique(edges_all[edges_all$condition %in% conditions, c("condition", "compound_id", "uniprot_id"), drop = FALSE])
  result <- .network_upsert(
    proj, "network_bowtie", result, c("condition", "compound_id", "uniprot_id"),
    touched_keys = bowtie_touched
  )
  summary_result <- .network_upsert(
    proj, "network_bowtie_summary", summary_row, c("species", "version"),
    touched_keys = data.frame(species = species, version = actions_version, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "network_bowtie") <- result
  patliRResults(proj, "network_bowtie_summary") <- summary_result
  .write_results_csv(proj, "network_bowtie", result)
  .write_results_csv(proj, "network_bowtie_summary", summary_result)
  proj <- .log_append(
    proj, step = "network_bowtie", id = NA_character_,
    message = paste0(
      "species ", species, " v", actions_version, ": ", summary_row$n_nodes, " nodes / ", summary_row$n_edges,
      " directed action edges -- core=", summary_row$n_core, ", in=", summary_row$n_in_component,
      ", out=", summary_row$n_out_component, ", other=", summary_row$n_other
    )
  )
  .write_log_csv(proj)
  proj
}

#' Download (if not already cached) and parse STRING's directed "actions"
#' file, caching the parsed raw table under `cacheDir(proj)/stringdb/`
#'
#' @description
#' Just the download + parse + column-sanity-check, shared by
#' `.network_stringdb_actions_graph()` (below) and, independently,
#' `plot_target_chord()`'s own raw-score reader (same on-disk cache path,
#' `<species>.protein.actions.v<version>.txt.gz`, read separately there
#' because it needs the per-edge `score` after the graph has already been
#' built). Caching the parsed data frame (not just the downloaded `.gz`)
#' means a later call with a *different* `score_threshold` never re-parses
#' the flat file, only re-filters an in-memory data frame.
#' @return A `data.frame` with (at least) STRING's own `item_id_a`,
#'   `item_id_b`, `is_directional`, `a_is_acting` columns, and `score` when
#'   the release provides it.
#' @keywords internal
.network_stringdb_actions_raw <- function(proj, species, version) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  cache_rds <- file.path(dir, paste0("actions_raw_", species, "_", version, ".rds"))
  if (file.exists(cache_rds)) {
    return(readRDS(cache_rds))
  }

  raw_gz <- file.path(dir, paste0(species, ".protein.actions.v", version, ".txt.gz"))
  if (!file.exists(raw_gz)) {
    url <- paste0(
      "https://stringdb-downloads.org/download/protein.actions.v", version, "/",
      species, ".protein.actions.v", version, ".txt.gz"
    )
    utils::download.file(url, raw_gz, mode = "wb", quiet = FALSE)
  }

  actions <- utils::read.delim(gzfile(raw_gz), stringsAsFactors = FALSE)
  required <- c("item_id_a", "item_id_b", "is_directional", "a_is_acting")
  missing_cols <- setdiff(required, names(actions))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "STRING's {.file protein.actions} file for species {.val {species}}/v{version} is missing expected column(s) {.val {missing_cols}}.",
      "i" = "Columns found: {.val {names(actions)}}. STRING may have changed its file format -- check https://string-db.org/cgi/download before retrying."
    ))
  }
  saveRDS(actions, cache_rds)
  actions
}

#' Parse STRING's directed "actions" file into an `igraph`, optionally
#' filtered on the actions file's own per-interaction `score`
#'
#' @description
#' Only rows where both `is_directional` and `a_is_acting` are true (see
#' `.network_actions_flag_true()`, internal -- STRING encodes these as
#' `"t"`/`"f"` strings, not `1`/`0`) become edges (`item_id_a -> item_id_b`)
#' -- interactions STRING itself does not mark
#' as directional (e.g. plain physical `binding`) are excluded, since
#' including them as if they had a direction would fabricate information
#' the source data does not contain. `item_id_a`/`item_id_b` are already in
#' STRING's own `<species>.ENSPxxxxxxx` ID space -- the same space
#' `STRINGdb$get_graph()`/`$map()` use -- so no additional ID mapping is
#' needed to build this graph.
#'
#' @param score_threshold `NULL` (default) applies no confidence filter --
#'   every directed row becomes an edge, the historical behaviour every
#'   existing caller (`plot_target_chord()`, which does its own separate
#'   score-based arc filtering afterwards) still relies on. A number filters
#'   to `score >= score_threshold` first ([network_bowtie()] passes its
#'   `actions_score_threshold` here).
#' @return `list(graph = <directed igraph>)`.
#' @keywords internal
.network_stringdb_actions_graph <- function(proj, species, version, score_threshold = NULL) {
  actions <- .network_stringdb_actions_raw(proj, species, version)

  is_directional <- .network_actions_flag_true(actions$is_directional)
  a_is_acting <- .network_actions_flag_true(actions$a_is_acting)
  keep <- is_directional & a_is_acting

  score_note <- ""
  if (!is.null(score_threshold)) {
    if (!"score" %in% names(actions)) {
      cli::cli_abort(c(
        "{.arg score_threshold} = {score_threshold} was given but STRING's {.file protein.actions} file for species {.val {species}}/v{version} has no {.field score} column.",
        "i" = "Columns found: {.val {names(actions)}}. Pass {.code score_threshold = NULL} to skip confidence filtering, or check STRING's file format."
      ))
    }
    keep <- keep & !is.na(actions$score) & actions$score >= score_threshold
    score_note <- paste0(" and {.field score} >= ", score_threshold)
  }

  directed <- unique(actions[keep, c("item_id_a", "item_id_b")])
  names(directed) <- c("from", "to")
  if (nrow(directed) == 0) {
    cli::cli_abort(c(
      paste0(
        "STRING's {.file protein.actions} file for species {.val {species}}/v{version} yielded zero directed edges after filtering on {.field is_directional}/{.field a_is_acting}",
        score_note, "."
      ),
      "i" = "This almost always means STRING changed how those columns are encoded (seen so far: {.val t}/{.val f}, {.val 1}/{.val 0}), or the score threshold is too high for this release -- inspect the raw actions file under cacheDir(proj)/stringdb/ directly before retrying."
    ))
  }

  g <- igraph::graph_from_data_frame(directed, directed = TRUE)
  list(graph = g)
}

#' UniProt -> STRING_id mapping table, read directly from STRING's own
#' `protein.aliases.v<version>` flat file
#'
#' @description
#' `STRINGdb$map()` performs this exact join (`multi_map_df()` against
#' `get_aliases()`, confirmed by reading `STRINGdb`'s own source) but doing
#' it via a `STRINGdb` object means instantiating one -- and
#' `STRINGdb$new(score_threshold = 0, ...)` downloads/loads the FULL
#' interactome (hundreds of MB, minutes) just to reach an alias table that
#' never depended on score at all. This reads the alias flat file directly
#' instead, same download-once/cache/parse pattern as
#' `.network_stringdb_actions_raw()`.
#' @return A `data.frame(string_protein_id, alias)`, deduplicated, one row
#'   per (STRING ID, alias) pair -- `alias` includes UniProt accessions
#'   alongside gene symbols/Ensembl IDs/etc, so filtering to the accessions
#'   of interest is a plain `%in%`/`match()`, no special-casing needed.
#' @keywords internal
.network_stringdb_aliases_map <- function(proj, species, version) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  cache_rds <- file.path(dir, paste0("aliases_raw_", species, "_", version, ".rds"))
  if (file.exists(cache_rds)) {
    return(readRDS(cache_rds))
  }

  raw_gz <- file.path(dir, paste0(species, ".protein.aliases.v", version, ".txt.gz"))
  if (!file.exists(raw_gz)) {
    url <- paste0(
      "https://stringdb-downloads.org/download/protein.aliases.v", version, "/",
      species, ".protein.aliases.v", version, ".txt.gz"
    )
    utils::download.file(url, raw_gz, mode = "wb", quiet = FALSE)
  }

  aliases <- .network_read_aliases_file(raw_gz)
  required <- c("string_protein_id", "alias")
  missing_cols <- setdiff(required, names(aliases))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "STRING's {.file protein.aliases} file for species {.val {species}}/v{version} is missing expected column(s) {.val {missing_cols}}.",
      "i" = "Columns found: {.val {names(aliases)}}. STRING may have changed its file format -- check https://string-db.org/cgi/download before retrying."
    ))
  }
  aliases <- unique(aliases[, c("string_protein_id", "alias")])
  saveRDS(aliases, cache_rds)
  aliases
}

## Second real bug found on the first live run (2026-07-23, after fixing the
## v12.0/v11.0 URL issue above): STRING's `is_directional`/`a_is_acting`
## columns are encoded as the strings "t"/"f" (confirmed against
## biostars.org/p/9588100 and STRING's own column documentation), not
## numeric 0/1 -- `actions$is_directional == 1` therefore matched nothing
## (a character "t" is never `== 1`) and silently produced an empty graph
## (0 nodes/edges) instead of erroring, because
## `igraph::graph_from_data_frame()` on a zero-row edge list is not itself
## an error. Handled tolerantly here (covers "t"/"f", "true"/"false", "1"/
## "0", and native logical/numeric, case-insensitively) since STRING's own
## docs are not fully consistent about which encoding a given release uses.
#' @keywords internal
.network_actions_flag_true <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  tolower(trimws(as.character(x))) %in% c("t", "true", "1")
}

#' Classify every node of a directed graph into the simplified bowtie
#' structure (core / in_component / out_component / other)
#' @return A named character vector (names = node names, values = one of
#'   `"core"`, `"in_component"`, `"out_component"`, `"other"`).
#' @keywords internal
.network_bowtie_classify <- function(g) {
  node_names <- igraph::V(g)$name
  comp <- igraph::components(g, mode = "strong")
  core_idx <- which.max(comp$csize)
  core_nodes <- node_names[comp$membership == core_idx]

  result <- stats::setNames(rep("other", length(node_names)), node_names)
  result[core_nodes] <- "core"

  if (length(core_nodes) > 0) {
    rep_node <- core_nodes[1]
    can_reach_core <- is.finite(igraph::distances(g, v = igraph::V(g), to = rep_node, mode = "out")[, 1])
    reached_by_core <- is.finite(igraph::distances(g, v = rep_node, to = igraph::V(g), mode = "out")[1, ])

    result[can_reach_core & result != "core"] <- "in_component"
    result[reached_by_core & result != "core" & result != "in_component"] <- "out_component"
  }

  result
}

#' @keywords internal
.empty_network_bowtie_row <- function() {
  data.frame(condition = character(0), compound_id = character(0), uniprot_id = character(0),
             string_id = character(0), bowtie_component = character(0), stringsAsFactors = FALSE)
}

#' Read a STRING `protein.aliases` flat file
#'
#' The v11.0 file's header line is space-separated
#' (`## string_protein_id ## alias ## source ##`) while the data rows are
#' tab-separated, so `read.delim()` sees a one-column header over three-column
#' rows and fails with "more columns than column names". Newer releases use a
#' tab-separated `#string_protein_id<TAB>alias<TAB>source` header. This reads
#' both: the header is parsed on its own and the body is read headerless.
#'
#' @param path Path to the (gzipped) aliases file.
#' @return `data.frame(string_protein_id, alias, source)`.
#' @keywords internal
.network_read_aliases_file <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1)
  cols <- if (grepl("\t", header, fixed = TRUE)) {
    sub("^#+", "", strsplit(header, "\t", fixed = TRUE)[[1]])
  } else {
    trimws(gsub("#", " ", strsplit(header, "[[:space:]]+#+[[:space:]]*")[[1]]))
  }
  cols <- cols[nzchar(cols)]
  body <- utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(body) < length(cols)) cols <- cols[seq_len(ncol(body))]
  names(body)[seq_along(cols)] <- cols
  body
}
