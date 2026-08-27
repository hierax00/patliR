#' @include AllGenerics.R internal.R network_build.R network_proximity.R
NULL

## network_bowtie() -- confirmed with Uriel (chat, 2026-07-23) after two
## substrates turned out unusable for a classic bowtie decomposition:
##   1. The directed compound-target-pathway-disease layered graph
##      (network_layers.R) is a DAG by construction -- no cycles, so its
##      largest strongly connected component is always a single node
##      (see network_motifs.R's roxygen for that discovery).
##   2. STRINGdb's regular interaction graph (`get_graph()`, already used
##      by network_proximity()) has real cycles, but is UNDIRECTED --
##      bowtie's core/in-component/out-component split is a directed-graph
##      concept (core = largest *strongly connected* component), so an
##      undirected graph collapses "core" back to "the whole connected
##      component" -- the same trivial-collapse problem, from the opposite
##      direction (cycles without direction, instead of direction without
##      cycles).
##
## Real fix: STRING separately publishes a directed "actions" channel
## (protein.actions.v<version>/<species>...txt.gz -- confirmed against
## STRING's own file-format documentation and the 2025 STRING paper on
## directionality of regulation, PMC11701646, 2026-07-23) with `action`
## (activation/inhibition/...), `is_directional`, and `a_is_acting` columns
## marking which protein acts on which for interactions with a known
## mechanism. This file is NOT exposed by any `STRINGdb` R package method
## (confirmed against its full method list, rdrr.io, 2026-07-23) -- it is
## downloaded and parsed directly here, cached under the same
## `cacheDir(proj)/stringdb/` directory `.network_stringdb()` already uses
## (same "safe to delete, redownloads" cache principle). It covers a much
## sparser subset of the interactome than the regular "links" file (only
## mechanistically-characterized interactions get a direction), which is
## the real, unavoidable price of getting a genuine directed substrate.
##
## `actions_version` defaults to "11.0", NOT the pipeline-wide `version`
## (default "12.0"): confirmed by first-run failure (2026-07-23, 404 on
## protein.actions.v12.0) and then against STRING's own download page
## (string-db.org/cgi/download) and directory listings on
## stringdb-downloads.org -- the "actions" flat file was discontinued
## after v11.0 and was never republished for v11.5 or v12.0. `version`
## still defaults to "12.0" for the UniProt -> STRING_id mapping
## (`.network_stringdb()`, shared with the rest of the `network_*` family)
## since that mapping is unaffected and should stay current; STRING protein
## IDs (`<species>.ENSPxxxxxxx`) are stable Ensembl accessions across
## versions, so mixing a v12.0 mapping with a v11.0 actions graph does not
## introduce an ID-space mismatch.
##
## Since bowtie describes the architecture of a whole network (Broder et
## al. 2000's original web-graph bowtie was global, not computed per query)
## rather than a per-condition property, the decomposition itself is
## computed ONCE per (species, version) and cached -- then each condition's
## mapped targets are annotated with which bowtie component they fall into,
## the same "compute the shared substrate once, join per condition" pattern
## network_proximity() already uses for the interactome itself.

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
#' @inheritParams network_proximity
#' @param actions_version STRING release to use for the directed "actions"
#'   file specifically (see `.network_stringdb_actions_graph()`, internal).
#'   Defaults to `"11.0"`, independently of `version` (default `"12.0"`,
#'   used for the UniProt -> STRING_id mapping): STRING discontinued the
#'   "actions" flat file after v11.0 and never republished it for v11.5 or
#'   v12.0, so `version`'s default does not work for this file and the two
#'   are tracked separately.
#'
#' @return The updated `proj`, with two entries in [patliRResults()]:
#'   `network_bowtie` (columns `condition`, `compound_id`, `uniprot_id`,
#'   `string_id`, `bowtie_component` -- one row per (condition, compound,
#'   target); targets that could not be mapped to the STRING action network
#'   get `bowtie_component = "unmapped"`, not a dropped row) and
#'   `network_bowtie_summary` (columns `species`, `version`, `n_nodes`,
#'   `n_edges`, `n_core`, `n_in_component`, `n_out_component`, `n_other` --
#'   one row per `(species, version)`, describing the global decomposition
#'   regardless of which conditions were annotated; `version` here is
#'   `actions_version`, the release the actions graph was actually built
#'   from). Both written to their matching `results/*.csv`.
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
#' proj <- network_bowtie(proj, condition = "FLO-ET") # needs STRINGdb + internet on first call
#' patliRResults(proj, "network_bowtie")
#' patliRResults(proj, "network_bowtie_summary")
#' }
#'
#' @export
network_bowtie <- function(proj, condition = NULL, species = 9606, version = "12.0",
                            actions_version = "11.0") {
  stopifnot(is(proj, "PatliRProject"))
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_bowtie} needs {.pkg STRINGdb} (for UniProt -> STRING_id mapping), not installed.",
      "i" = "See {.file TESTING_GUIDE.Rmd}, section 0, for the {.fn BiocManager::install} chunk."
    ))
  }
  conditions <- .network_resolve_conditions(proj, condition)

  actions <- .network_stringdb_actions_graph(proj, species, actions_version)
  g <- actions$graph
  bowtie <- .network_bowtie_classify(g)

  summary_row <- data.frame(
    species = species, version = actions_version,
    n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g),
    n_core = sum(bowtie == "core"), n_in_component = sum(bowtie == "in_component"),
    n_out_component = sum(bowtie == "out_component"), n_other = sum(bowtie == "other"),
    stringsAsFactors = FALSE
  )

  string_db <- .network_stringdb(proj, species, version, score_threshold = 0)
  edges_all <- patliRResults(proj, "network_edges")
  rows <- vector("list", length(conditions))
  names(rows) <- conditions

  for (cond in conditions) {
    ct <- unique(edges_all[edges_all$condition == cond, c("compound_id", "uniprot_id")])
    if (nrow(ct) == 0) {
      rows[[cond]] <- .empty_network_bowtie_row()
      next
    }

    map_df <- string_db$map(data.frame(uniprot_id = unique(ct$uniprot_id), stringsAsFactors = FALSE),
                             "uniprot_id", removeUnmappedRows = FALSE, quiet = TRUE)
    unmapped <- map_df$uniprot_id[is.na(map_df$STRING_id)]
    if (length(unmapped) > 0) {
      proj <- .log_append(
        proj, step = "network_bowtie", id = unmapped,
        message = "network_bowtie_unmapped: no STRING_id found for this UniProt ID in this species; component = 'unmapped'"
      )
    }
    uni_to_string <- stats::setNames(map_df$STRING_id, map_df$uniprot_id)

    string_id <- uni_to_string[ct$uniprot_id]
    component <- ifelse(
      is.na(string_id), "unmapped",
      ifelse(string_id %in% names(bowtie), bowtie[string_id], "other")
    )

    rows[[cond]] <- data.frame(
      condition = cond, compound_id = ct$compound_id, uniprot_id = ct$uniprot_id,
      string_id = unname(string_id), bowtie_component = component, stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result <- .network_upsert(proj, "network_bowtie", result, c("condition", "compound_id", "uniprot_id"))
  summary_result <- .network_upsert(proj, "network_bowtie_summary", summary_row, c("species", "version"))

  patliRResults(proj, "network_bowtie") <- result
  patliRResults(proj, "network_bowtie_summary") <- summary_result
  .write_results_csv(proj, "network_bowtie", result)
  .write_results_csv(proj, "network_bowtie_summary", summary_result)
  proj <- .log_append(
    proj, step = "network_bowtie", id = NA_character_,
    message = paste0(
      "species ", species, " v", version, ": ", summary_row$n_nodes, " nodes / ", summary_row$n_edges,
      " directed action edges -- core=", summary_row$n_core, ", in=", summary_row$n_in_component,
      ", out=", summary_row$n_out_component, ", other=", summary_row$n_other
    )
  )
  .write_log_csv(proj)
  proj
}

#' Download (if not already cached) and parse STRING's directed "actions"
#' file into an `igraph`, caching the parsed graph under
#' `cacheDir(proj)/stringdb/`
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
#' @return `list(graph = <directed igraph>)`.
#' @keywords internal
.network_stringdb_actions_graph <- function(proj, species, version) {
  dir <- file.path(cacheDir(proj), "stringdb")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  cache_rds <- file.path(dir, paste0("actions_graph_", species, "_", version, ".rds"))
  if (file.exists(cache_rds)) {
    return(list(graph = readRDS(cache_rds)))
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

  is_directional <- .network_actions_flag_true(actions$is_directional)
  a_is_acting <- .network_actions_flag_true(actions$a_is_acting)
  directed <- actions[is_directional & a_is_acting, c("item_id_a", "item_id_b")]
  directed <- unique(directed)
  names(directed) <- c("from", "to")
  if (nrow(directed) == 0) {
    cli::cli_abort(c(
      "STRING's {.file protein.actions} file for species {.val {species}}/v{version} yielded zero directed edges after filtering on {.field is_directional}/{.field a_is_acting}.",
      "i" = "This almost always means STRING changed how those two columns are encoded (seen so far: {.val t}/{.val f}, {.val 1}/{.val 0}) -- inspect {.file {raw_gz}} directly before retrying."
    ))
  }

  g <- igraph::graph_from_data_frame(directed, directed = TRUE)
  saveRDS(g, cache_rds)
  list(graph = g)
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
