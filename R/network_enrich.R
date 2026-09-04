#' @include AllGenerics.R internal.R network_build.R
NULL

## ReactomePA::enrichPathway() only accepts Entrez gene IDs (no keyType
## argument), so this function maps UniProt -> Entrez once via
## clusterProfiler::bitr() and feeds Entrez to whichever backend was
## requested, rather than relying on per-function keyType support that is
## inconsistent across enrichGO()/enrichKEGG()/enrichPathway().

#' Pathway/ontology enrichment over each condition's target set
#'
#' @description
#' For each condition already built by [network_build()], maps its distinct
#' UniProt target IDs to Entrez Gene IDs (`org.Hs.eg.db`) and runs an
#' over-representation test against GO, Reactome, or KEGG gene sets. It does
#' not change the graph built by [network_build()]; it produces a separate
#' enrichment results table keyed by condition.
#'
#' @section Heavy, optional dependencies -- not installed by default:
#' Unlike `igraph` (a hard `Imports` dependency of the whole `network_*`
#' family since [network_build()]), the packages this function needs are
#' large Bioconductor annotation databases, kept in `Suggests` rather than
#' `Imports` so that installing `patliR` itself stays light. You need:
#' `org.Hs.eg.db` and `clusterProfiler` always (ID mapping, plus GO/KEGG
#' enrichment); `ReactomePA` additionally for `db = "reactome"`. Install
#' with `BiocManager::install(c("clusterProfiler", "org.Hs.eg.db",
#' "ReactomePA"))`. Calling this
#' without the right package installed raises a clear error naming exactly
#' which package is missing, rather than a cryptic "could not find
#' function".
#'
#' @section KEGG needs internet, GO/Reactome do not:
#' `org.Hs.eg.db` and `reactome.db`/`ReactomePA` ship their gene-set data
#' as local packages -- once installed, `db = "go"`/`db = "reactome"` never
#' touch the network. `db = "kegg"` is different:
#' `clusterProfiler::enrichKEGG()` queries KEGG's REST API live on every
#' call (there is no current local Bioconductor package for this -- the old
#' `KEGG.db` is deprecated). GO and Reactome are fully local once their
#' packages are installed.
#'
#' @inheritParams compounds
#' @param condition Character vector of condition names (must already have
#'   been built by [network_build()]), or `NULL` (default) for every built
#'   condition.
#' @param db `"reactome"` (default), `"go"`, or `"kegg"` -- one at a time;
#'   call again with a different `db` to accumulate results (existing
#'   `(condition, db)` pairs are replaced, not duplicated).
#' @param ont Only used when `db = "go"`: `"BP"` (default), `"MF"`, `"CC"`,
#'   or `"ALL"` -- passed straight to [clusterProfiler::enrichGO()].
#' @param pvalueCutoff,qvalueCutoff Passed straight to the underlying
#'   `clusterProfiler`/`ReactomePA` enrichment function (defaults `0.05`/
#'   `0.2`, same as their own defaults).
#' @param simplify_go Logical, default `TRUE`. Only used when `db = "go"`.
#'   GO terms are hierarchical and highly redundant -- a real target set
#'   routinely comes back with hundreds of "significant" GO terms, many of
#'   them near-duplicate parent/child terms covering almost the same genes.
#'   This redundancy is what makes [network_motifs()]'s derived
#'   `compound -> pathway` layer explode to compound out-degrees in the
#'   thousands, since every near-duplicate GO term adds its own derived
#'   edge. When `TRUE`, runs
#'   [clusterProfiler::simplify()] (semantic-similarity-based, via
#'   \pkg{GOSemSim}) on the `db = "go"` result before returning it -- this
#'   removes redundant terms by *similarity*, not by tightening
#'   `pvalueCutoff`/`qvalueCutoff` further, so it does not silently drop
#'   real signal the way a stricter cutoff would. If `GOSemSim` is not
#'   installed, or `simplify()` itself fails for this gene set (has
#'   happened for `ont = "ALL"` in some `clusterProfiler` versions -- it is
#'   designed for one ontology at a time), the unsimplified result is kept
#'   and the fallback is logged (`"network_enrich_simplify_failed"`), never
#'   a hard error -- this is enrichment quality, not a required step.
#' @param simplify_cutoff Single number in `(0, 1]`, default `0.7` (same as
#'   `clusterProfiler::simplify()`'s own default). Lower keeps fewer, more
#'   dissimilar terms (more aggressive collapsing); `1` keeps everything
#'   `simplify()` would otherwise consider a near-duplicate. Only used when
#'   `simplify_go = TRUE` and `db = "go"`.
#'
#' @return The updated `proj`, with a `network_enrichment` entry in
#'   [patliRResults()] (columns `condition`, `db`, `ID`, `Description`,
#'   `GeneRatio`, `BgRatio`, `pvalue`, `p.adjust`, `qvalue`, `geneID`
#'   (Entrez IDs, `/`-separated -- `clusterProfiler`'s own convention,
#'   preserved as-is rather than reformatted), `Count`), also written to
#'   `results/network_enrichment.csv`. UniProt IDs that could not be mapped
#'   to an Entrez ID are logged (`"network_enrich_unmapped_target"`) and
#'   excluded from that condition's gene set, never silently included as
#'   `NA`.
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
#' proj <- network_enrich(proj, db = "go") # needs org.Hs.eg.db + clusterProfiler
#' patliRResults(proj, "network_enrichment")
#' }
#'
#' @export
network_enrich <- function(proj, condition = NULL,
                            db = c("reactome", "go", "kegg"),
                            ont = c("BP", "MF", "CC", "ALL"),
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                            simplify_go = TRUE, simplify_cutoff = 0.7) {
  stopifnot(is(proj, "PatliRProject"))
  db <- match.arg(db)
  ont <- match.arg(ont)
  stopifnot(is.numeric(pvalueCutoff), length(pvalueCutoff) == 1, pvalueCutoff > 0, pvalueCutoff <= 1)
  stopifnot(is.numeric(qvalueCutoff), length(qvalueCutoff) == 1, qvalueCutoff > 0, qvalueCutoff <= 1)
  stopifnot(is.logical(simplify_go), length(simplify_go) == 1, !is.na(simplify_go))
  stopifnot(is.numeric(simplify_cutoff), length(simplify_cutoff) == 1, simplify_cutoff > 0, simplify_cutoff <= 1)

  .network_enrich_check_deps(db, simplify_go = simplify_go && db == "go")

  ## Built once, not once per condition: GOSemSim::godata() (called inside
  ## clusterProfiler::simplify() when semData isn't supplied) parses the
  ## whole GO DAG + org.Hs.eg.db's term-gene annotations, which is by far
  ## the slowest part of simplify() -- reusing it across every condition in
  ## this call turns an O(n_conditions) rebuild into a one-time cost. Shared
  ## construction + on-disk cache with network_degeneracy() via
  ## .network_godata() (annoDb=, not the deprecated OrgDb=); it already
  ## wraps godata() in tryCatch(..., NULL).
  go_sem_data <- if (db == "go" && simplify_go) {
    .network_godata(proj, ont = if (ont == "ALL") "BP" else ont, computeIC = FALSE)
  } else {
    NULL
  }

  edges_all <- patliRResults(proj, "network_edges")
  if (is.null(edges_all) || nrow(edges_all) == 0) {
    cli::cli_abort(c(
      "No {.val network_edges} entry in {.arg proj}.",
      "i" = "Run {.fn network_build} first."
    ))
  }
  built_conditions <- unique(edges_all$condition)

  if (!is.null(condition)) {
    unknown <- setdiff(condition, built_conditions)
    if (length(unknown) > 0) {
      cli::cli_abort("Condition(s) {.val {unknown}} were not built by {.fn network_build}; available: {.val {built_conditions}}.")
    }
    conditions <- condition
  } else {
    conditions <- built_conditions
  }

  result_list <- vector("list", length(conditions))
  names(result_list) <- conditions

  for (cond in conditions) {
    uniprot_ids <- unique(edges_all$uniprot_id[edges_all$condition == cond])
    uniprot_ids <- uniprot_ids[!is.na(uniprot_ids) & nzchar(uniprot_ids)]

    if (length(uniprot_ids) == 0) {
      proj <- .log_append(
        proj, step = "network_enrich", id = NA_character_,
        message = paste0("condition '", cond, "': no targets to enrich (empty gene set)")
      )
      result_list[[cond]] <- .empty_network_enrichment_row()
      next
    }

    map <- clusterProfiler::bitr(uniprot_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
    unmapped <- setdiff(uniprot_ids, map$UNIPROT)
    if (length(unmapped) > 0) {
      proj <- .log_append(
        proj, step = "network_enrich", id = unmapped,
        message = "network_enrich_unmapped_target: no Entrez Gene ID found in org.Hs.eg.db for this UniProt ID; excluded from the gene set"
      )
    }
    entrez_ids <- unique(map$ENTREZID)

    if (length(entrez_ids) == 0) {
      proj <- .log_append(
        proj, step = "network_enrich", id = NA_character_,
        message = paste0("condition '", cond, "': no UniProt ID mapped to Entrez; nothing to enrich")
      )
      result_list[[cond]] <- .empty_network_enrichment_row()
      next
    }

    enrich_result <- tryCatch(
      .network_enrich_run(db, entrez_ids, ont, pvalueCutoff, qvalueCutoff),
      error = function(e) e
    )
    if (inherits(enrich_result, "error")) {
      proj <- .log_append(
        proj, step = "network_enrich", id = NA_character_,
        message = paste0("condition '", cond, "', db '", db, "': enrichment call failed -- ", conditionMessage(enrich_result))
      )
      result_list[[cond]] <- .empty_network_enrichment_row()
      next
    }

    if (db == "go" && simplify_go && !is.null(enrich_result) && nrow(as.data.frame(enrich_result)) > 0) {
      if (ont == "ALL") {
        proj <- .log_append(
          proj, step = "network_enrich", id = NA_character_,
          message = paste0("condition '", cond, "': simplify_go skipped -- clusterProfiler::simplify() operates on one GO ontology at a time, not ont = \"ALL\"")
        )
      } else {
        simplified <- tryCatch(
          clusterProfiler::simplify(enrich_result, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min, semData = go_sem_data),
          error = function(e) e
        )
        if (inherits(simplified, "error")) {
          proj <- .log_append(
            proj, step = "network_enrich", id = NA_character_,
            message = paste0("network_enrich_simplify_failed: condition '", cond, "' -- keeping the unsimplified GO result -- ", .cli_escape(conditionMessage(simplified)))
          )
        } else {
          n_before <- nrow(as.data.frame(enrich_result))
          n_after <- nrow(as.data.frame(simplified))
          proj <- .log_append(
            proj, step = "network_enrich", id = NA_character_,
            message = paste0("condition '", cond, "': simplify_go reduced ", n_before, " GO term(s) to ", n_after, " (cutoff = ", simplify_cutoff, ")")
          )
          enrich_result <- simplified
        }
      }
    }

    df <- as.data.frame(enrich_result)
    if (nrow(df) == 0) {
      result_list[[cond]] <- .empty_network_enrichment_row()
    } else {
      df$condition <- cond
      df$db <- db
      result_list[[cond]] <- df[, c("condition", "db", "ID", "Description", "GeneRatio", "BgRatio",
                                     "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    }
  }

  new_rows <- do.call(rbind, result_list)
  rownames(new_rows) <- NULL

  existing <- patliRResults(proj, "network_enrichment")
  if (!is.null(existing) && nrow(existing) > 0) {
    touched <- paste(new_rows$condition, new_rows$db)
    existing_key <- paste(existing$condition, existing$db)
    existing <- existing[!existing_key %in% touched, , drop = FALSE]
    new_rows <- rbind(existing, new_rows)
  }

  patliRResults(proj, "network_enrichment") <- new_rows
  .write_results_csv(proj, "network_enrichment", new_rows)
  .write_log_csv(proj)
  proj
}

#' Check that the packages `network_enrich()` needs for `db` (and, if
#' requested, `simplify_go`) are installed
#' @keywords internal
.network_enrich_check_deps <- function(db, simplify_go = FALSE) {
  needed <- c("clusterProfiler", "org.Hs.eg.db")
  if (db == "reactome") needed <- c(needed, "ReactomePA")
  if (simplify_go) needed <- c(needed, "GOSemSim")

  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cli::cli_abort(c(
      "{.fn network_enrich} with {.code db = \"{db}\"}{if (simplify_go) ' and simplify_go = TRUE' else ''} needs the following package(s), not installed: {.val {missing}}.",
      "i" = "{.code BiocManager::install(c({paste(sprintf('\"%s\"', missing), collapse = ', ')}))}",
      "i" = "Or pass {.code simplify_go = FALSE} to skip GO-term redundancy reduction if only {.pkg GOSemSim} is missing."
    ))
  }
}

#' Run the actual enrichment call for one (already ID-mapped) gene set
#' @return An `enrichResult` object (from `clusterProfiler`/`ReactomePA`).
#' @keywords internal
.network_enrich_run <- function(db, entrez_ids, ont, pvalueCutoff, qvalueCutoff) {
  if (db == "reactome") {
    return(ReactomePA::enrichPathway(
      gene = entrez_ids, organism = "human",
      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff
    ))
  }
  if (db == "go") {
    return(clusterProfiler::enrichGO(
      gene = entrez_ids, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = ont,
      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff
    ))
  }
  ## db == "kegg" -- hits KEGG's live REST API, needs internet
  clusterProfiler::enrichKEGG(
    gene = entrez_ids, organism = "hsa",
    pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff
  )
}

#' A genuinely empty (0-row) result -- used when a condition has nothing to
#' enrich (empty/unmapped gene set) or the enrichment call itself failed.
#' Never carries a placeholder `condition`/`db` row: like the rest of
#' `patliR`, "nothing to report" means no row, not a row full of `NA`s --
#' the real story (why) lives in `projectLog(proj)`.
#' @keywords internal
.empty_network_enrichment_row <- function() {
  data.frame(
    condition = character(0), db = character(0), ID = character(0), Description = character(0),
    GeneRatio = character(0), BgRatio = character(0), pvalue = double(0),
    p.adjust = double(0), qvalue = double(0), geneID = character(0), Count = integer(0),
    stringsAsFactors = FALSE
  )
}
