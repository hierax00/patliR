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
#' @section Background universe -- why the default is not the whole genome:
#' Every over-representation test asks "is this term hit more than you would
#' expect against the background?" -- and "expect" is meaningless until the
#' background is specified (Boyle et al. 2004, *Bioinformatics* 20(18),
#' 3710-3715, the original GO::TermFinder framing). The foreground gene set
#' here is never a random sample of the genome: it is the reachable output
#' of a ligand-similarity target predictor (SuperPred/SwissTargetPrediction),
#' which can only return proteins with enough known ligands to build a
#' similarity model in the first place -- overwhelmingly GPCRs, kinases,
#' nuclear receptors, proteases, and transporters. Tested against all
#' ~20,000 genes annotated in `org.Hs.eg.db` (`clusterProfiler`'s own
#' default background), those families come out "significant" for *any*
#' input compound before a single biological difference between compounds
#' is considered -- Timmons, Szkop & Gallagher (2015), *Genome Biol* 16:186,
#' "Multiple sources of bias confound functional enrichment analysis of
#' global -omics data". `universe = "project"` (the default here) fixes
#' this by restricting the background to the project's own predictable
#' proteome -- every protein any compound in the project could plausibly
#' have been assigned as a target, not every protein in the genome -- so
#' the test asks the sharper question "is this term hit more than you'd
#' expect from *this predictor's* output," not "...from a random gene."
#' `universe = "genome"` is kept as an explicit, opt-in escape hatch back to
#' the old whole-genome/whole-pathway-database behaviour, for comparison or
#' for callers who have a specific reason to want it.
#'
#' **This is a breaking change for existing projects.** Before this
#' argument existed, every call implicitly ran with today's `"genome"`
#' behaviour. Re-running `network_enrich()` under the new `"project"`
#' default changes every p-value it reports, and therefore every downstream
#' result that consumes `network_enrichment` --
#' [network_degeneracy(annotation = "enriched")][network_degeneracy()] /
#' `annotation = "jaccard"`, [network_motifs()]'s pathway layer, and
#' [plot_network_layers()]/[plot_enrichment()]. See `NEWS.md`.
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
#' @param universe `"project"` (default) or `"genome"` -- the
#'   over-representation background; see the "Background universe" section
#'   above. `"project"`
#'   maps `unique(patliRResults(proj, "targets_imported")$uniprot_id)` (or,
#'   if that slot is absent, every target across every condition
#'   [network_build()] has built -- the same fallback and the same pool
#'   [network_degeneracy(universe = "project")][network_degeneracy()] uses)
#'   to Entrez once per call (not once per condition, since the pool does
#'   not depend on which condition is being enriched) and passes it as
#'   `universe =` to [clusterProfiler::enrichGO()]/
#'   [clusterProfiler::enrichKEGG()]/[ReactomePA::enrichPathway()].
#'   `"genome"` passes no `universe =` at all, so each function falls back
#'   to its own whole-genome/whole-pathway-database default. Recorded as a
#'   per-row `universe` column.
#' @param pvalueCutoff,qvalueCutoff Passed straight to the underlying
#'   `clusterProfiler`/`ReactomePA` enrichment function (defaults `0.05`/
#'   `0.2`, same as their own defaults).
#' @param pAdjustMethod One of `stats::p.adjust.methods` (default `"BH"`,
#'   same as `clusterProfiler`'s own default when this argument did not
#'   exist). Passed straight to [clusterProfiler::enrichGO()]/
#'   [clusterProfiler::enrichKEGG()]/[ReactomePA::enrichPathway()].
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
#'   `GeneRatio`, `BgRatio`, `pvalue`, `p.adjust`, `qvalue` (`NA` when
#'   `clusterProfiler` omitted it -- routine when the `qvalue` package's
#'   pi0 estimation fails on a small gene set, never a crash here),
#'   `geneID` (Entrez IDs, `/`-separated -- `clusterProfiler`'s own
#'   convention, preserved as-is rather than reformatted), `Count`,
#'   `ONTOLOGY` (only populated when `ont = "ALL"`; `NA` otherwise, kept
#'   for a stable schema across calls), `universe` (the resolved `universe`
#'   argument for that row, `"project"` or `"genome"`)), also written to
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
                            universe = c("project", "genome"),
                            pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                            pAdjustMethod = "BH",
                            simplify_go = TRUE, simplify_cutoff = 0.7) {
  stopifnot(is(proj, "PatliRProject"))
  db <- match.arg(db)
  ont <- match.arg(ont)
  universe <- match.arg(universe)
  pAdjustMethod <- match.arg(pAdjustMethod, choices = stats::p.adjust.methods)
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

  ## Background universe, built ONCE per call (not once per condition): it
  ## does not depend on which condition is being enriched. Same source and
  ## the same UniProt -> Entrez mapping pattern as
  ## network_degeneracy(universe = "project")'s `project_pool_entrez` --
  ## intentionally kept identical so the two `universe = "project"`
  ## concepts mean the same thing across the package. See the "Background
  ## universe" @section above for why this exists.
  universe_entrez <- NULL
  ## What actually ends up passed to the enrichment call -- may fall back
  ## to "genome" below even when the caller asked for "project"; the
  ## per-condition result rows must record this effective value, not the
  ## raw argument, or a fallback would be invisible in the durable table.
  universe_effective <- universe
  if (universe == "project") {
    ti <- patliRResults(proj, "targets_imported")
    proj_uni <- if (!is.null(ti) && nrow(ti) > 0) unique(ti$uniprot_id) else unique(edges_all$uniprot_id)
    universe_entrez <- .network_uniprot_to_entrez(proj_uni)
    if (length(universe_entrez) == 0) {
      cli::cli_warn(c(
        "!" = "{.fn network_enrich}: {.code universe = \"project\"} resolved to 0 Entrez-mapped gene(s) from {length(proj_uni)} UniProt accession(s).",
        "i" = "Falling back to no background restriction for this call (equivalent to {.code universe = \"genome\"})."
      ))
      universe_entrez <- NULL
      universe_effective <- "genome"
    } else {
      proj <- .log_append(
        proj, step = "network_enrich", id = NA_character_,
        message = paste0(
          "universe = 'project': background pool = ", length(universe_entrez),
          " Entrez gene(s) mapped from ", length(proj_uni), " distinct UniProt accession(s)"
        )
      )
    }
  }

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
      .network_enrich_run(db, entrez_ids, ont, pvalueCutoff, qvalueCutoff, pAdjustMethod, universe_entrez),
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
    result_list[[cond]] <- if (nrow(df) == 0) {
      .empty_network_enrichment_row()
    } else {
      .network_enrich_result_df(df, cond, db, universe_effective)
    }
  }

  new_rows <- do.call(rbind, result_list)
  rownames(new_rows) <- NULL

  new_rows <- .network_upsert(
    proj, "network_enrichment", new_rows, c("condition", "db"),
    touched_keys = data.frame(condition = conditions, db = db, stringsAsFactors = FALSE)
  )

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
#'
#' @param universe_entrez `NULL` (no background restriction -- each
#'   function's own whole-genome/whole-pathway-database default) or a
#'   character vector of Entrez Gene IDs to pass as `universe =`. Confirmed
#'   against the installed `clusterProfiler`/`ReactomePA` formals that all
#'   three underlying functions accept `universe =` in the same ID space as
#'   `gene =` (`enrichKEGG(keyType = "kegg")`'s KEGG-prefixed IDs are only
#'   used internally when querying KEGG's REST API -- both `gene =` and
#'   `universe =` are documented as "a vector of entrez gene id", and
#'   `DOSE::enricher_internal()`, which all three call into, intersects
#'   `universe` against its own background gene set without re-keying it).
#'   Passing `universe = NULL` explicitly is equivalent to omitting the
#'   argument (`enrichGO()` does `if (missing(universe)) universe <- NULL`;
#'   `enricher_internal()`'s own default is `universe = NULL` and it
#'   no-ops the background restriction when `is.null(universe)`), so this
#'   can always be passed rather than conditionally.
#' @return An `enrichResult` object (from `clusterProfiler`/`ReactomePA`).
#' @keywords internal
.network_enrich_run <- function(db, entrez_ids, ont, pvalueCutoff, qvalueCutoff,
                                 pAdjustMethod = "BH", universe_entrez = NULL) {
  if (db == "reactome") {
    return(ReactomePA::enrichPathway(
      gene = entrez_ids, organism = "human",
      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
      pAdjustMethod = pAdjustMethod, universe = universe_entrez
    ))
  }
  if (db == "go") {
    return(clusterProfiler::enrichGO(
      gene = entrez_ids, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", ont = ont,
      pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
      pAdjustMethod = pAdjustMethod, universe = universe_entrez
    ))
  }
  ## db == "kegg" -- hits KEGG's live REST API, needs internet
  clusterProfiler::enrichKEGG(
    gene = entrez_ids, organism = "hsa",
    pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
    pAdjustMethod = pAdjustMethod, universe = universe_entrez
  )
}

#' Build the per-row `network_enrichment` columns from one condition's
#' `as.data.frame(enrichResult)`, defensively.
#'
#' @description
#' Replaces a hard `df[, c(...)]` column subset, which threw `undefined
#' columns selected` whenever `clusterProfiler` omitted `qvalue` --
#' routine when the `qvalue` package's pi0 estimation fails on the small
#' p-value vectors a 10-50-gene target set produces (would otherwise abort
#' a run that had already spent minutes on GO/KEGG/Reactome). `ONTOLOGY` is
#' handled the same defensive way: `clusterProfiler` only adds it when
#' `ont = "ALL"`; keeping it present-but-`NA` otherwise (rather than
#' absent) matches the zero-row-invariant, stable-schema convention
#' `.empty_network_enrichment_row()` documents for this table.
#' @keywords internal
.network_enrich_result_df <- function(df, cond, db, universe) {
  qvalue_col   <- if ("qvalue" %in% names(df)) df$qvalue else rep(NA_real_, nrow(df))
  ontology_col <- if ("ONTOLOGY" %in% names(df)) df$ONTOLOGY else rep(NA_character_, nrow(df))
  data.frame(
    condition = cond, db = db,
    ID = df$ID, Description = df$Description,
    GeneRatio = df$GeneRatio, BgRatio = df$BgRatio,
    pvalue = df$pvalue, p.adjust = df$p.adjust, qvalue = qvalue_col,
    geneID = df$geneID, Count = df$Count,
    ONTOLOGY = ontology_col, universe = universe,
    stringsAsFactors = FALSE
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
    ONTOLOGY = character(0), universe = character(0),
    stringsAsFactors = FALSE
  )
}
