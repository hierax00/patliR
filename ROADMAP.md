# patliR roadmap

Designed but not yet built. Each item should become a GitHub issue once the
repo is public; this file is the interim list.

Everything the help pages document is implemented. The functions below are
design only — asking for one (e.g. `network_build(target_source =
"consensus")`) raises a clear "not implemented yet" error rather than
guessing.

## Target prediction

- **`targets_bipartite()` / `targets_consensus()`** — RWR/NBI link
  prediction over a compound-target bipartite graph from known bioactivity,
  then a consensus score across sources. Deferred: the cost is replicating
  ChEMBL's data curation well enough to trust an in-house predictor, not
  the algorithm. `targets_consensus()` needs at least two sources to
  combine, so it falls with the same decision.
- **`disease_genes_import()`** — manual alternative to
  `targets_disease_filter()`'s Open Targets call: import a disease-gene
  list into the same `targets_disease` slot, with a `source` column.

## Reference / natural-product databases

- **`coconut_fetch()`** — COCONUT 2.0 REST API, folded into
  `reference_*`. Endpoint not confirmed without guessing at the time of
  writing.
- **`tcm_import()`** — manual import for traditional-medicine databases
  (UniTCM, HERB) that lack a reliable API.

## Ranking and reporting

- **`rank_candidates()`** — combine adjusted centrality (hub penalty + bias
  reweight), ADME filtering, consensus targets, and modular robustness into
  one ranked table. Blocks `plot_rank()`.
- **`dock_prepare()` / `dock_parse()`** — compounds -> PDBQT + a search-box
  config (`box_center` required, no default); parse Vina/DiffDock output
  back to a data frame. No scoring of its own.
- **`report_generate()`** — one report per condition: ADME filtering,
  ranking, bias audit, modular robustness, full decision/parameter/seed
  log. Docking section only if `dock_results` exists.
- **`plot_rank()` / `plot_kegg_binding()`** — blocked by the two families
  above.

## Network

- **`network_kegg_complete()`** — add KEGG's directed topology to the
  network. Deferred on time budget.
- **`network_hub_penalty()` scoring** — `d * log(N/d)` is non-monotone
  (peaks at `d = N/e`), so a moderately promiscuous target can outrank a
  selective one. Decide whether "intermediate-specificity emphasis" is
  actually wanted or switch to a monotone selectivity weight.
- **`network_synergy()`** — replace `1 - Jaccard(targets)` with a real
  network separation `s_AB` (Menche et al. 2015) and report Cheng et al.
  (2019)'s Complementary Exposure classification, not just the scalar.
- **`network_module_robustness()` clustering** — HDBSCAN on integer graph
  distances degenerates to single-linkage; add modularity-based options
  (`cluster_leiden`, or bipartite modularity) via the `clustering=` arg.
- **`network_degeneracy()`** — swap `pathway_jaccard` for GO semantic
  similarity (GOSemSim, already a dependency); add a permutation null.
- **`network_centrality()`** — bipartite-aware normalisation so `degree` /
  `betweenness` are comparable across `node_type` (Borgatti & Everett 1997).
- **`bias_audit()` categorical enrichment** — the half of the original
  design with no code: MeSH / Disease Ontology enrichment vs. a
  STRING/DrugBank/reference background. Needs a confirmed category source.
- **AI interpretation module ("narrator")** — a curated narrative over a
  project's results. `patliR_export_llm()` is the raw-data-dump precursor.

## Extensions to functions that already exist

- **`plot_chemical_space()`** — overlay targets/proteins on the same space.
  A joint compound+target embedding is a genuinely different design
  question (what shared axis?), left for its own design pass.
- **`tox_safetyome()`** — swap the 500-gene core panel for the full
  ~11,300-gene catalogue if the core panel proves too narrow (only changes
  which CSV `.load_safetyome_core_panel()` reads).
- **`network_filter_proteome()`** — wire the filtered edge set back into
  `network_centrality()` / `network_module_robustness()` / etc., and
  support guided import from an external file. Needs a new scope axis
  parallel to `condition` in `.network_resolve_conditions()` /
  `.network_graph()`.
- **`compounds_similarity()`** — integrate with `plot_chemical_space()` /
  `network_synergy()` (currently standalone).
- **`adme_import()`** — give it the same `mapping_file` row-order bridge
  `tox_import()` now has, so its export/import loop closes without relying
  on cross-toolkit SMILES canonicalization either.

## Housekeeping

- Write `vignettes/patliR-intro.Rmd` and `paper.md` (JOSS) before any
  Bioconductor/JOSS submission.
- Reconcile "second run" semantics across families (`network_*` vs.
  `adme_local()` etc.), give `plot_*` consistent logging, and factor out
  the boilerplate repeated across the 16 `plot_*` functions.
- Audit `\donttest` examples: those needing network should be `\dontrun`,
  those needing Bioconductor Suggests should be `@examplesIf`.
