#' @include AllGenerics.R internal.R network_build.R network_enrich.R
NULL

## Adds KEGG's directed pathway topology (activation/inhibition/binding/...
## relations between genes WITHIN a pathway) to the network -- the piece
## ROADMAP.md flagged as deferred. network_pathview() only ever renders
## KEGG's own pre-made pathway images (pathview::pathview()); nothing in
## patliR before this file parsed KGML relations itself.
##
## Pattern follows network_bowtie()'s STRING precedent: fetch and parse the
## raw source file directly (KGML via KEGGREST::keggGet(id, "kgml")) rather
## than a heavier wrapper, cache the raw download, map IDs once in bulk
## (KEGGREST::keggConv("uniprot", species) -- one call for the whole
## species, not one per gene). New Suggests: KEGGREST (fetch) + xml2
## (parse) -- both lightweight; KEGGgraph/Rgraphviz were deliberately not
## used (heavier stack for the same handful of XML tags this file reads
## directly).
##
## KGML shape, confirmed live against a real pathway (hsa04151) before
## writing this parser:
##   <entry id="6" name="hsa:1977 hsa:253314 hsa:9470" type="gene" ...>
##     (an entry can bundle several KEGG gene IDs -- an ortholog/paralog box)
##   <entry id="12" type="group" ...>
##     <component id="10"/><component id="11"/>
##   </entry>
##     (a complex -- its own "genes" are the union of its components',
##     each of which is itself a type="gene" entry elsewhere in the file)
##   <relation entry1="95" entry2="43" type="PPrel">
##     <subtype name="activation" value="--&gt;"/>
##   </relation>
##     (a relation can carry more than one <subtype>; each becomes its own
##     output row rather than being collapsed to just the first)
## entry1/entry2 are KGML-internal entry IDs, not gene IDs -- must be
## resolved via the entry table first, then expanded (cartesian product)
## when either side bundles multiple genes (directly, or via a group).

#' Fetch and add KEGG's directed pathway topology to the network
#'
#' @description
#' For each KEGG pathway in scope, downloads its KGML (KEGG Markup
#' Language) file and parses the `<relation>` elements -- directed,
#' typed gene-gene relations (`PPrel` protein-protein, e.g. activation/
#' inhibition/binding; `GErel` gene expression; `ECrel` enzyme-catalysis)
#' -- into UniProt-keyed edges, including relations between `type="group"`
#' entries (protein complexes), resolved to the union of their member
#' genes. This is the piece [network_pathview()] does not provide: that
#' function only renders KEGG's own pre-made pathway diagrams (via the
#' `pathview` package), it never parses topology; nothing else in
#' `patliR` reads a KGML file before this function. `PCrel`
#' (protein-compound) relations are parsed but always dropped -- see
#' `relation_types` below for why.
#'
#' @section Where the pathway list comes from:
#' `pathway_ids = NULL` (default) pulls every KEGG pathway ID already found
#' significant by [network_enrich(db = "kegg")][network_enrich()] for the
#' condition(s) in scope -- so this function completes a network you have
#' already characterised, rather than dumping an unrelated bulk topology.
#' Pass `pathway_ids` explicitly (e.g. `"hsa04151"`) to bypass that and fetch
#' specific pathways regardless of enrichment.
#'
#' @section `restrict_to_network` (default `TRUE`):
#' A KGML relation is kept only when **both** its endpoints, after UniProt
#' mapping, are among the condition's own predicted targets
#' (`network_edges$uniprot_id`) -- directly usable to annotate the existing
#' compound-target graph with directionality/mechanism, not a disconnected
#' topology dump. `restrict_to_network = FALSE` keeps the full pathway
#' topology (every gene KEGG maps in that pathway, whether or not any
#' compound in this project happens to target it) -- useful for pathway-level
#' analysis independent of the current target predictions.
#'
#' @inheritParams network_build
#' @param pathway_ids `NULL` (default, see above) or a character vector of
#'   KEGG pathway IDs (e.g. `c("hsa04151", "hsa04010")`).
#' @param species KEGG organism code, default `"hsa"` (human) -- must match
#'   what [network_enrich()]/[network_pathview()] were run against.
#' @param relation_types Character vector of KGML relation `type`s to keep.
#'   Default `c("PPrel", "GErel", "ECrel")` -- every gene-gene relation type
#'   KGML defines. `"PCrel"` (protein-compound) is a valid value you can
#'   pass explicitly, but is not in the default: this function only ever
#'   resolves `<entry type="gene">` (and, as of the group-entry fix, the
#'   `type="group"` complexes built from them) to UniProt IDs, never
#'   `<entry type="compound">` -- so a `PCrel` relation's compound-side
#'   endpoint can never resolve and that relation is always dropped as an
#'   unmapped pair. Kept as an explicit opt-in (rather than removed
#'   outright) so a future patch that does add compound-entry resolution
#'   doesn't need to change this argument's contract, only what happens
#'   when it's used. `"maplink"` (link to another pathway diagram, not a
#'   biological relation) is never included.
#' @param restrict_to_network Logical, default `TRUE`. See the section above.
#'
#' @return The updated `proj`, with a `network_kegg_topology` entry in
#'   [patliRResults()] (columns `condition`, `pathway_id`, `pathway_title`,
#'   `from_kegg`, `to_kegg`, `from_uniprot`, `to_uniprot`, `relation_type`,
#'   `relation_subtype`, `relation_value` (KGML's own arrow-style code, e.g.
#'   `"-->"` activation, `"--|"` inhibition -- kept verbatim, not
#'   reinterpreted), `restrict_to_network`), also written to
#'   `results/network_kegg_topology.csv`. Entries KEGG could not map to a
#'   UniProt accession are logged and excluded, never silently dropped.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' # ... build the project, run network_build() and
#' # network_enrich(condition = "FLO-ET", db = "kegg") ...
#' proj <- network_kegg_topology(proj, condition = "FLO-ET") # needs internet
#' patliRResults(proj, "network_kegg_topology")
#' }
#'
#' @seealso [network_enrich()], [network_pathview()]
#' @export
network_kegg_topology <- function(proj, condition = NULL, pathway_ids = NULL,
                                   species = "hsa",
                                   relation_types = c("PPrel", "GErel", "ECrel"),
                                   restrict_to_network = TRUE) {
  stopifnot(is(proj, "PatliRProject"))
  stopifnot(is.character(species), length(species) == 1, nzchar(species))
  stopifnot(is.logical(restrict_to_network), length(restrict_to_network) == 1, !is.na(restrict_to_network))
  relation_types <- match.arg(relation_types, c("PPrel", "GErel", "ECrel", "PCrel"), several.ok = TRUE)
  if ("PCrel" %in% relation_types) {
    cli::cli_warn(c(
      "!" = "{.code relation_types} includes {.val PCrel}, but {.fn network_kegg_topology} never resolves {.code entry type=\"compound\"} to a UniProt ID.",
      "i" = "Every {.val PCrel} relation will be dropped as unmapped; see {.arg relation_types}' documentation."
    ))
  }

  if (!requireNamespace("KEGGREST", quietly = TRUE) || !requireNamespace("xml2", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.fn network_kegg_topology} needs {.pkg KEGGREST} and {.pkg xml2}, not both installed.",
      "i" = "Install with {.code BiocManager::install(\"KEGGREST\")} and {.code install.packages(\"xml2\")}."
    ))
  }

  conditions <- .network_resolve_conditions(proj, condition)
  edges_all <- patliRResults(proj, "network_edges")

  rows <- vector("list", length(conditions))
  names(rows) <- conditions
  uniprot_of_kegg <- NULL # built once, lazily, shared across every condition/pathway in this call

  for (cond in conditions) {
    ids <- pathway_ids
    if (is.null(ids)) {
      enrichment <- patliRResults(proj, "network_enrichment")
      if (!is.null(enrichment) && nrow(enrichment) > 0) {
        ids <- unique(enrichment$ID[enrichment$condition == cond & enrichment$db == "kegg"])
      }
      if (length(ids) == 0) {
        cli::cli_abort(c(
          "No {.arg pathway_ids} given and no KEGG-enriched pathways found for condition {.val {cond}}.",
          "i" = "Pass {.arg pathway_ids} explicitly, or run {.fn network_enrich}({.arg condition} = {.val {cond}}, {.arg db} = \"kegg\") first."
        ))
      }
    }

    target_ids <- unique(edges_all$uniprot_id[edges_all$condition == cond])
    pathway_rows <- vector("list", length(ids))
    names(pathway_rows) <- ids

    for (pid in ids) {
      kgml <- .fetch_external(
        fetch_fun = function() KEGGREST::keggGet(pid, "kgml"),
        cache_dir = file.path(cacheDir(proj), "kegg_topology"), cache_key = pid,
        mode = "warn_and_cache"
      )
      if (is.null(kgml) || !nzchar(kgml)) {
        proj <- .log_append(
          proj, step = "network_kegg_topology", id = pid,
          message = "network_kegg_topology_fetch_failed: could not download KGML for this pathway; skipped"
        )
        next
      }

      parsed <- tryCatch(.network_kegg_parse_kgml(kgml, relation_types), error = function(e) e)
      if (inherits(parsed, "error")) {
        proj <- .log_append(
          proj, step = "network_kegg_topology", id = pid,
          message = paste0("network_kegg_topology_parse_failed: ", conditionMessage(parsed))
        )
        next
      }
      if (nrow(parsed$relations) == 0) {
        pathway_rows[[pid]] <- NULL
        next
      }

      if (is.null(uniprot_of_kegg)) {
        conv <- KEGGREST::keggConv("uniprot", species)
        uniprot_of_kegg <- split(sub("^up:", "", unname(conv)), names(conv))
      }

      rel <- parsed$relations
      from_up <- uniprot_of_kegg[rel$from_kegg]
      to_up <- uniprot_of_kegg[rel$to_kegg]
      n_expand <- lengths(from_up) * lengths(to_up)
      unmapped <- lengths(from_up) == 0 | lengths(to_up) == 0
      if (any(unmapped)) {
        proj <- .log_append(
          proj, step = "network_kegg_topology", id = pid,
          message = paste0(
            "network_kegg_topology_unmapped: ", sum(unmapped), " of ", nrow(rel),
            " relation(s) in ", pid, " had a KEGG gene ID with no UniProt mapping; excluded"
          )
        )
      }
      keep <- !unmapped
      if (!any(keep)) { pathway_rows[[pid]] <- NULL; next }

      expanded <- do.call(rbind, lapply(which(keep), function(i) {
        expand.grid(
          from_kegg = rel$from_kegg[i], to_kegg = rel$to_kegg[i],
          from_uniprot = from_up[[i]], to_uniprot = to_up[[i]],
          relation_type = rel$relation_type[i], relation_subtype = rel$relation_subtype[i],
          relation_value = rel$relation_value[i],
          stringsAsFactors = FALSE
        )
      }))

      if (isTRUE(restrict_to_network)) {
        expanded <- expanded[expanded$from_uniprot %in% target_ids & expanded$to_uniprot %in% target_ids, , drop = FALSE]
        if (nrow(expanded) == 0) { pathway_rows[[pid]] <- NULL; next }
      }

      expanded <- unique(expanded)
      pathway_rows[[pid]] <- data.frame(
        condition = cond, pathway_id = pid, pathway_title = parsed$title, expanded,
        stringsAsFactors = FALSE
      )
    }

    pathway_rows <- pathway_rows[!vapply(pathway_rows, is.null, logical(1))]
    rows[[cond]] <- if (length(pathway_rows) == 0) {
      .empty_network_kegg_topology_row()
    } else {
      do.call(rbind, pathway_rows)
    }
    proj <- .log_append(
      proj, step = "network_kegg_topology", id = NA_character_,
      message = paste0(
        "condition '", cond, "': ", nrow(rows[[cond]]), " KEGG-topology edge(s) from ",
        length(pathway_rows), " of ", length(ids), " requested pathway(s), restrict_to_network = ", restrict_to_network
      )
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (is.null(result) || nrow(result) == 0) result <- .empty_network_kegg_topology_row()
  result$restrict_to_network <- restrict_to_network

  result <- .network_upsert(
    proj, "network_kegg_topology", result, "condition",
    touched_keys = data.frame(condition = conditions, stringsAsFactors = FALSE)
  )

  patliRResults(proj, "network_kegg_topology") <- result
  .write_results_csv(proj, "network_kegg_topology", result)
  .write_log_csv(proj)
  proj
}

#' Parse one pathway's KGML: the gene entries and the directed relations
#' between them (KGML entry IDs already resolved to KEGG gene IDs,
#' cartesian-expanded when either endpoint bundles more than one gene)
#'
#' @return `list(title, relations)`, `relations` a `data.frame(from_kegg,
#'   to_kegg, relation_type, relation_subtype, relation_value)` -- one row
#'   per (KEGG gene, KEGG gene) pair the relation expands to.
#' @keywords internal
.network_kegg_parse_kgml <- function(kgml_text, relation_types) {
  doc <- xml2::read_xml(kgml_text)
  title <- xml2::xml_attr(doc, "title")

  entries <- xml2::xml_find_all(doc, "//entry[@type='gene']")
  entry_genes <- stats::setNames(
    lapply(xml2::xml_attr(entries, "name"), function(x) strsplit(x, "\\s+")[[1]]),
    xml2::xml_attr(entries, "id")
  )

  ## `type="group"` entries (protein complexes) have no `name` of their
  ## own -- KGML lists their members as <component id="..."/> children,
  ## each referencing a type="gene" entry elsewhere in this same document.
  ## A relation naming a complex as entry1/entry2 (common for receptor
  ## complexes) would otherwise never resolve to any gene and be silently
  ## dropped by the `is_gene_pair` check below. Resolve each group to the
  ## union of its members' KEGG gene IDs and fold it into the same lookup
  ## table, keyed by the group's own entry id, so the rest of this
  ## function needs no group-specific branch.
  group_entries <- xml2::xml_find_all(doc, "//entry[@type='group']")
  if (length(group_entries) > 0) {
    group_genes <- stats::setNames(
      lapply(group_entries, function(g) {
        comp_ids <- xml2::xml_attr(xml2::xml_find_all(g, "component"), "id")
        ## A component id not found in entry_genes (e.g. a nested group,
        ## against spec but not worth erroring over) is dropped, not an
        ## error -- unlist() on a NULL element contributes nothing.
        unlist(entry_genes[comp_ids], use.names = FALSE)
      }),
      xml2::xml_attr(group_entries, "id")
    )
    group_genes <- group_genes[lengths(group_genes) > 0]
    entry_genes <- c(entry_genes, group_genes)
  }

  rel_nodes <- xml2::xml_find_all(doc, "//relation")
  if (length(rel_nodes) == 0) {
    return(list(title = title, relations = .network_kegg_empty_relations()))
  }

  rel_type <- xml2::xml_attr(rel_nodes, "type")
  keep <- rel_type %in% relation_types
  rel_nodes <- rel_nodes[keep]
  rel_type <- rel_type[keep]
  if (length(rel_nodes) == 0) {
    return(list(title = title, relations = .network_kegg_empty_relations()))
  }

  e1 <- xml2::xml_attr(rel_nodes, "entry1")
  e2 <- xml2::xml_attr(rel_nodes, "entry2")
  is_gene_pair <- e1 %in% names(entry_genes) & e2 %in% names(entry_genes)
  if (!any(is_gene_pair)) {
    return(list(title = title, relations = .network_kegg_empty_relations()))
  }

  ## A <relation> can carry more than one <subtype> (e.g. a PPrel with
  ## both "phosphorylation" and "activation" stacked) -- xml_find_first()
  ## used to keep only the first and silently drop the rest. One output
  ## row per subtype now (falling back to a single NA-subtype row when a
  ## relation has none, same as before), cartesian-expanded against the
  ## from/to gene sets exactly like the multi-gene-entry case already was.
  idx <- which(is_gene_pair)
  out <- do.call(rbind, lapply(idx, function(i) {
    subtype_nodes <- xml2::xml_find_all(rel_nodes[[i]], "subtype")
    if (length(subtype_nodes) == 0) {
      sub_name <- NA_character_
      sub_value <- NA_character_
    } else {
      sub_name <- xml2::xml_attr(subtype_nodes, "name")
      sub_value <- xml2::xml_attr(subtype_nodes, "value")
    }
    g <- expand.grid(
      from_kegg = entry_genes[[e1[i]]], to_kegg = entry_genes[[e2[i]]],
      subtype_idx = seq_along(sub_name), stringsAsFactors = FALSE
    )
    data.frame(
      from_kegg = g$from_kegg, to_kegg = g$to_kegg,
      relation_type = rel_type[i],
      relation_subtype = sub_name[g$subtype_idx], relation_value = sub_value[g$subtype_idx],
      stringsAsFactors = FALSE
    )
  }))
  list(title = title, relations = unique(out))
}

#' @keywords internal
.network_kegg_empty_relations <- function() {
  data.frame(
    from_kegg = character(0), to_kegg = character(0), relation_type = character(0),
    relation_subtype = character(0), relation_value = character(0), stringsAsFactors = FALSE
  )
}

#' @keywords internal
.empty_network_kegg_topology_row <- function() {
  data.frame(
    condition = character(0), pathway_id = character(0), pathway_title = character(0),
    from_kegg = character(0), to_kegg = character(0), from_uniprot = character(0), to_uniprot = character(0),
    relation_type = character(0), relation_subtype = character(0), relation_value = character(0),
    stringsAsFactors = FALSE
  )
}
