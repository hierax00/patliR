# patliR (development version)

## `network_module_robustness()` — clustering backends and a random-failure baseline

- **Breaking change to the numbers.** `network_module_robustness()`'s
  `clustering` argument now defaults to **`"leiden"`** (was `"hdbscan"`).
  Modules are detected by `igraph::cluster_leiden()` maximising Newman
  modularity, so every module partition — and therefore every `r_index` —
  differs from earlier versions. `clustering = "hdbscan"` still works
  (unchanged output for the same seed) but emits a one-time note that it is
  a density heuristic on an integer metric; it is not deprecated.
- **New `clustering = "bipartite"`** — `bipartite::computeModules(method =
  "Beckett")`, maximising Barber (2007)'s bipartite modularity `Q_B`. Only
  valid on the genuine two-mode compound–target graph (aborts otherwise).
  Adds `bipartite` to `Suggests`; note it attaches `sna` and `vegan`, which
  mask some `igraph` functions for the rest of the session.
- Clustering is run **unweighted**: for Leiden the `weight` edge attribute
  is *deleted* before clustering rather than passed as `weights = NA`
  (which would leave `strength()` reading it for the modularity null model
  and silently return all-singletons on any `NA` probability).
- **New `attack = c("targeted", "random", "both")`** (default
  `"targeted"`). `"random"` averages `n_random` (default 20) uniformly
  random removal orders; `network_robustness_curve` gains
  `removal_strategy`, `largest_component_sd`, `n_replicates`. The summary
  table gains `r_index_random` and `n_random` alongside the unchanged
  `r_index` (targeted).
- **New third output slot `network_module_membership`** — one row per node
  (`condition`, `node_id`, `node_type`, `module_id`, `module_type`),
  required by `plot_network_layers(colour_by = "module")`.
- The summary table gains `clustering`, `modularity` (Newman `Q` / mean
  Barber `Q_B` / `NA`), `resolution`, and `method_detail`. Rerunning with a
  different `clustering` **replaces** the condition's rows (the `clustering`
  column records which method produced them); a pre-0.2.0
  `network_module_robustness.csv` is back-filled with `NA`, not rejected.
- New arguments: `resolution` (1), `n_iterations` (5), `min_component_size`
  (3, the method-agnostic single-module fallback threshold), `n_random`
  (20).
- `plot_robustness()` now draws both removal strategies (targeted solid,
  random dashed with an SD ribbon), plots the x axis as *fraction removed*,
  and annotates `R_targ` / `R_rand` / their difference.

## Disease gene sets (`disease_genes_*`)

- **New `disease_genes_fetch()`** — fetches a disease's associated gene set
  from Open Targets in the **disease -> target** direction
  (`disease(efoId:){ associatedTargets{ rows{ target{ proteinIds } score } } }`,
  paginated on `count`), mapping `proteinIds` with `source ==
  "uniprot_swissprot"` to UniProt accessions (a target with several
  Swiss-Prot accessions contributes one row each; a target with none is
  logged and dropped). Stored in a new `disease_genes` results slot
  (`results/disease_genes.csv`; columns `disease_id`, `disease_name`,
  `uniprot_id`, `gene_symbol`, `association_score`, `source`, `fetched_at`),
  upserted per `disease_id`. `disease` is resolved to a current EFO/MONDO
  ID the same way `targets_disease_filter()` does it; `min_score` mirrors
  `targets_disease_filter()`'s semantics and defaults to `NULL` (keep all).
- **New `disease_genes_import()`** — the manual counterpart, for curated
  lists (OMIM / GWAS Catalog / DisGeNET exports). Takes a data frame or a
  CSV path with a UniProt-ID column (gene-symbol-only tables are rejected
  with a clear message -- convert to accessions first); `source` records
  provenance and `association_score` may be `NA`. Same slot / schema /
  per-`disease_id` upsert as `disease_genes_fetch()`. Implements the
  `disease_genes_import()` item from the roadmap.

## Network analysis (`network_*` / `plot_*`)

- **`network_proximity()` no longer builds its disease module from the
  compounds' own predicted targets.** It gains `disease_genes = c(
  "disease_genes", "targets_disease")` (default `"disease_genes"`): the
  default path builds the disease set `T` from the independent
  `disease_genes` slot (`disease_genes_fetch()` / `disease_genes_import()`),
  so the proximity z-score measures topology rather than set membership. It
  aborts if that slot is missing or empty for the requested disease. The
  legacy `disease_genes = "targets_disease"` path still runs but emits a
  `cli_warn` naming the **circular**ity (`targets_disease` only annotates
  UniProt IDs already predicted as targets, so `S` and `T` overlap by
  construction). Output gains two columns: `disease_gene_source`
  (`"disease_genes"` / `"targets_disease"`) and `n_overlap` = `|S ∩ T|` on
  the STRING-mapped sets (equals `n_targets_mapped` exactly when
  `d_observed == 0`, i.e. `S ⊆ T`). Both are declared in the empty-row
  constructor, so `results/network_proximity.csv` stays column-stable
  across zero-row reruns.

- **`network_centrality()` gains bipartite-aware normalisation
  (`normalize = TRUE`, the default).** For whichever `measures` were
  requested it adds `degree_norm` (degree over the opposite mode's size)
  and `betweenness_norm` (betweenness over the mode-specific maximum
  betweenness of a bipartite graph of these mode sizes), both in `[0, 1]`
  and comparable between compounds and targets of the same condition
  (Borgatti & Everett 1997, *Social Networks* 19:243; `B_max` cross-checked
  against NetworkX's `bipartite.betweenness_centrality` -- no factor-of-two
  applied, `igraph::betweenness()` already counts each unordered pair
  once). It also adds `hub_score_component` (HITS hub recomputed per
  connected component, rescaled to `max = 1` within each -- fixes the
  disconnected-graph case where every node outside the largest component
  got `hub_score` ~ 0) with a `component_id` label, and the per-condition
  constants `n_compounds` / `n_targets`. Edge cases: an empty mode gives
  `_norm = NA` (never `NaN`/`0`); `B_max == 0` (a mode of size 1) gives
  `NA` with a log line. `normalize = FALSE` reproduces the previous raw
  column set exactly. All new columns are declared in the empty-row
  constructor, so `results/network_centrality.csv` stays column-stable
  across zero-row reruns; `.network_upsert()` back-fills them on a project
  whose CSV predates this change. `_norm` values are **not** comparable
  across conditions of different size (documented like the R-index).
- **`plot_centrality()`**: `measure` accepts `"degree_norm"` /
  `"betweenness_norm"`; when `measure` is not supplied and
  `node_type = "both"` it now defaults to `"degree_norm"` (with a
  `cli_inform` naming the reason), falling back to `"degree"` with a
  warning on results produced before normalisation. The value axis states
  the normaliser (e.g. `/ |T| = 412`). Two nodes that resolve to the same
  gene symbol no longer collapse into one bar -- the accession is appended
  on collision.
- **`.network_upsert()` no longer keeps stale rows when a rerun produces
  zero rows for a condition.** It gains a `touched_keys` argument: callers
  now pass the key values they just recomputed, so a slot that now yields
  nothing (a proteome filter that matches no target, a condition with no
  motifs, a compound that lost its last STRING-mappable target, ...) has
  its previous rows dropped instead of silently surviving. Migrated every
  `network_*` caller (`network_build` excepted -- it never had the bug):
  `network_filter_proteome`, `network_motifs`, `network_degeneracy`,
  `network_synergy`, `network_proximity`, `network_centrality`,
  `network_hub_penalty`, `network_module_robustness`, `network_bowtie`,
  `network_pathview`. A pre-`rbind` column-set check now fails with a clear
  message instead of R's generic "numbers of columns of arguments do not
  match".
- **`network_proximity()`**: `p_adjusted` is now part of the results
  schema even for a zero-row result (declared in the empty-row
  constructor, assigned `NA_real_` before the BH block), so
  `results/network_proximity.csv` stays column-stable across reruns.
- **`network_synergy(pairs = "all")`** no longer errors with `subscript
  out of bounds` when a compound is absent from `network_proximity()`
  (e.g. it had no STRING-mappable target); such pairs get `NA` in the
  proximity-derived columns, as the docs already promised.
- **`plot_bowtie()`** now colours the `not_in_action_network` component
  (added to `network_bowtie()` in the previous round) -- previously the
  sixth, often largest, stratum rendered as an unnamed grey/NA flow.
- **`plot_robustness()`** caption corrected: `r_index` is bounded above by
  `(N-1)/(2N) < 0.5` (Schneider et al. 2011), not "close to 1 = robust".
- **`network_proximity()` performance**: the degree-bin lookup
  (`bin_of_node`, `node_names`) is hoisted out of
  `.network_resample_degree_matched()` and computed once per call rather
  than ~`2 * n_random` times per compound; `.network_degree_bins()`
  replaces its per-unique-degree rescan with a single `rle()` pass.
- **Internal**: the boilerplate the seven network `plot_*` functions
  repeated (ggplot2/engine check, scope label, `ggsave` + log upsert +
  `girafe` wrap + `attr(result, "proj")` tail, node-label lookup) is
  factored into `.plot_require()` / `.plot_scope()` / `.plot_finish()` /
  `.plot_label_nodes()` in `R/plot-helpers.R`. `.plot_label_nodes()` is the
  former `.network_layers_labels()`, moved out of `plot_network_layers.R`
  (which `plot_bowtie`/`plot_centrality`/`plot_proximity`/`plot_synergy`
  `@include`d only to reach it). No change to any figure, filename, log
  schema or return type.

## Phase 1 verification follow-ups (`disease_genes_*` / `network_proximity()`)

- **`disease_genes_fetch(min_score = )` now defaults to `0.4`, not `NULL`.**
  Open Targets' full associated-target list for a common disease runs to
  thousands of genes (type 2 diabetes: ~9,900) against a ~17k-node STRING
  LCC; at that size roughly half of any compound's targets fall inside the
  disease module by chance and the proximity z-score's dynamic range
  collapses. `0.4` brings common diseases to a Menche/Guney-scale module
  (~200-300 genes). `NULL` stays available as "keep everything". The
  kept/dropped counts and the cutoff are logged as a single summary line.
- **`disease_genes_import()` gains `map_symbols = FALSE`.** When `TRUE`, a
  symbol-only `table` is mapped through `org.Hs.eg.db`
  (`clusterProfiler::bitr()`, `SYMBOL -> UNIPROT`); every multi-mapping and
  every unmapped symbol is logged, all accessions of a multi-mapping symbol
  are kept, and `source` records that the mapping was applied. Symbol
  rejection stays the default (the mapping is many-to-many).
- **Per-gene logging replaced by summaries.** `disease_genes_no_swissprot`,
  `disease_genes_below_min_score` and `network_proximity_unmapped` wrote one
  `log.csv` row per gene (~10k rows per call against a real disease). Each
  is now a single count line (with the first ~10 IDs inline).
- **`.network_upsert()` migrates a legacy results CSV instead of aborting.**
  When `new_rows` has columns the existing table lacks and none the other
  way (a pure column addition, e.g. a pre-Phase-1 `network_proximity.csv`
  without `disease_gene_source` / `n_overlap`), the old rows are back-filled
  with typed `NA` and the upsert proceeds with a one-line `cli_inform`. A
  column *removal* or a *type conflict* on a shared column still aborts.
- **`network_proximity()` abort message** on the default path no longer
  tells the user to check `targets_disease` (it names the `disease_genes`
  set that was actually in use).
- **`network_proximity()` null efficiency.** The degree-matched disease-side
  null is drawn once per `(condition, disease)` and reused across compounds
  (it was redrawn per compound although `T` is identical); `target_string`
  and `V(g)$name` are hoisted out of the per-compound loop. Pure speedups.
- **`.open_targets_graphql()`** now carries a `User-Agent`, a ~3 req/s
  throttle and exponential-backoff retry on 429/5xx (mirroring
  `.refdb_get_json()`) -- `disease_genes_fetch()` drives a ~20-request
  paginated loop against a rate-limiting API.
- **`plot_proximity()`** warns and drops the circular `targets_disease`
  rows when a `network_proximity` table mixes both `disease_gene_source`
  values, so they cannot render as the most significant points.

## Reference database (`refdb_*`)

- **`refdb_build()` now resolves ChEMBL identity by a deterministic
  identity chain, not a similarity search.** The old
  `/similarity/<smiles>/100.json` call matched a 100% Tanimoto
  fingerprint, which is not structural identity and assigned wrong
  molecules (e.g. methyl palmitate resolved to methyl octanoate, limonene
  to the wrong enantiomer). The new `.chembl_resolve()` walks PubChem CID
  -> standard InChIKey (batched) -> ChEMBL
  `molecule_structures__standard_inchi_key__in` exact (batched) ->
  InChIKey-skeleton `__startswith` fallback -> SMILES `flexmatch` (SMILES
  passed as a query parameter). Multi-hit queries are resolved by
  `.chembl_pick()` (own-parent, then `pref_name`, then lowest ChEMBL id),
  never "the first molecule".
- **`reference_compounds`** gains `inchikey`, `match_type`
  (`inchikey_exact` / `inchikey_skeleton` / `smiles_flexmatch`),
  `parent_chembl_id` and `fetched_at`.
- **`reference_bioactivity`** gains `standard_relation`, `pchembl_value`,
  `target_organism`, `assay_type` and `data_validity_comment`; the hard
  `limit = 50` on the ChEMBL activity fetch is gone (full pagination via
  `page_meta$next`), non-target rows ("No relevant target", "Log S", rows
  with no `target_chembl_id`) are dropped, and rows are deterministically
  ordered. A censored measurement (`IC50 > 100000`) is no longer stored as
  if it were an exact value.
- All reference-database HTTP now goes through one choke point with a
  `User-Agent`, a ~3 req/s throttle and exponential-backoff retry on
  429/5xx. Identity cache keys are chemical (CID / InChIKey / SMILES
  digest), not the project-local `C0001` id, so re-preparing compounds in a
  different order can no longer return a stale, wrong compound's data.
- Removed `.chembl_url_encode_smiles()` (dead once SMILES moved off the URL
  path).

**Action required:** projects built with an earlier version must rerun
`refdb_build()` and rebuild any downstream `bias_*` results -- the previous
`reference_*.csv` files contain mis-assigned ChEMBL molecules.

# patliR 0.1.0

Initial public release. Core pipeline implemented and tested, pre-1.0.

- `prep_*` — compound import/validation (PubChem CID or SMILES),
  replicate-abundance binarization, 2D structure depiction.
- `refdb_*` — reference database sync (PubChem, ChEMBL).
- `compounds_classify()` — natural-product family classification via
  NPClassifier; `compounds_similarity()` — pairwise fingerprint
  similarity (Tanimoto/Dice/cosine).
- `adme_*` — local physicochemical drug-likeness/lead-likeness rules
  (Lipinski Ro5, Veber, Ghose, Egan, Oprea) and BOILED-Egg via `rcdk`,
  plus import from external ADME platforms, rule-based filtering, and
  SMILES export for platforms that only take pasted-in lists
  (`adme_export_smiles()`).
- `tox_*` — structural toxicity alerts (PAINS, 480 filters; Brenk, 105
  alerts), target-level systemic-risk screening against the Safetyome
  core panel, import from external toxicity platforms
  (`tox_export_smiles()` + `tox_import(mapping_file=)` for a row-order
  round trip), and reporting (`tox_report()` writes
  `results/tox_report.csv`).
- `targets_*` (partial) — target import from prediction platforms,
  disease-target association filtering.
- `network_*` — the systems-biology core: build, enrich (GO/Reactome/
  KEGG), pathview, centrality/hub penalty, module robustness, layers/
  motifs, degeneracy, proximity, synergy, bow-tie, and
  `network_filter_proteome()` (filter to a custom protein list).
- `plot_*` — 16 visualization functions covering ADME, 2D structures,
  and every `network_*` result.
- `launch_app()` — optional Shiny wizard orchestrating the pipeline
  step by step.
- `bias_*` (partial) — `bias_audit()`/`bias_reweight()`/`bias_report()`
  flag and discount statistical-outlier compounds/targets in the
  project's own reference database.
- `patliR_export_llm()` — flat-text project export for LLM analysis.

See [`ROADMAP.md`](ROADMAP.md) for what's designed but not yet built.
