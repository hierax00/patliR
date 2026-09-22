# patliR (development version)

## New `targets_disease_profile()` / `plot_disease_network()`, and `network_kegg_topology()`

- **New `targets_disease_profile()`.** Organizes every predicted target's
  Open Targets disease associations, resolving the disease **display name**
  (`targets_disease_filter()` stores only the EFO/MONDO ID). Two modes:
  `disease = NULL` (default, explore mode) takes each target's top
  `top_n_diseases` associated diseases with no disease fixed in advance, so
  you see the landscape of diseases a compound's target set is actually
  implicated in; `disease = "<name or ID>"` (single-disease mode, same
  resolution rules as `targets_disease_filter()`) queries every target
  against just that one disease. Writes a new `targets_disease_profile`
  slot, distinct from and independent of `targets_disease_filter()`'s
  `targets_disease` table.
- **New `plot_disease_network()`.** Draws the compound-target-disease
  network from `targets_disease_profile()`'s output: compounds and targets
  at the usual small point size (same palette as `plot_network_layers()`),
  disease nodes drawn substantially larger with a translucent convex-hull
  halo over every target associated with that disease -- reuses the exact
  `grDevices::chull()` + `geom_polygon()` pattern `plot_chemical_space()`
  already established for its per-family halos, rather than inventing a new
  one.
- **New `network_kegg_topology()`.** Closes the `network_kegg_complete()`
  ROADMAP item: fetches KGML for each condition's KEGG-enriched pathways
  (`network_enrich(db = "kegg")`, or an explicit `pathway_ids=`) and parses
  `<relation>` elements -- directed, typed gene-gene relations (PPrel/
  GErel/ECrel/PCrel: activation, inhibition, binding, phosphorylation, ...)
  -- into UniProt-keyed edges (`KEGGREST` fetch + bulk `keggConv()` ID
  mapping, `xml2` parse; both new `Suggests`, deliberately lightweight --
  `KEGGgraph`/Rgraphviz were not used for the same handful of XML tags this
  reads directly, same "depend for algorithms, not for a few tags" judgment
  as `bipartite`/`STRINGdb`). `restrict_to_network = TRUE` (default) keeps
  only relations between two of the condition's own predicted targets --
  directly usable to annotate the existing compound-target graph with
  directionality/mechanism, not a disconnected topology dump; `FALSE` keeps
  the pathway's full topology.

## New `report_generate()` (Milestone B)

- **New `report_generate()`.** One self-contained HTML file per condition
  (`reports/report_<condition>.html`): compounds present, ADME filtering,
  the database-bias audit, modular robustness, the `rank_candidates()`
  ranking, and the full `projectLog()`. Each section auto-includes when its
  source table exists and shows a one-line "not run" placeholder otherwise
  -- the same graceful-degradation convention `rank_candidates()` already
  established, not a hard prerequisite abort (only `network_build()` is
  required, to know which compounds are "present" per condition).
  Deliberately plain HTML built with base R string concatenation, **not**
  `rmarkdown`/`pandoc` -- this environment has no `pandoc`
  (`rmarkdown::pandoc_available()` is `FALSE` here), so an `Rmd`-based
  report could not even be rendered or tested; and per `DESIGN.md`'s own
  "minimal Suggests" philosophy, a handful of HTML tables does not need a
  heavier dependency anyway. Table cell values are HTML-escaped
  (`&`/`<`/`>`) before embedding.

## New `rank_candidates()` / `plot_rank()` (Milestone B)

- **New `rank_candidates()`.** Combines per-compound criteria into one
  ranked table via Robust Rank Aggregation (RRA; Kolde, Laur, Adler & Vilo
  2012, *Bioinformatics* 28(4):573-580) -- criteria are compared by *rank*,
  not raw value, so a criterion with a broken or skewed scale (e.g. a
  `betweenness_norm` that sits at 0 for most nodes) cannot silently
  dominate a weighted sum. Two criteria are mandatory: ADME pass fraction
  (from `adme_filter()`'s `adme_filtered`) and network centrality (a
  `rowMeans()` composite of `network_centrality()`'s available `*_norm`
  columns, rolled up target -> compound by the `network_build()` edge
  weight). Four more are auto-detected and silently skipped when absent,
  with the decision logged: hub penalty (`network_hub_penalty`), disease
  proximity (`network_proximity`'s `z_score`, direction-flipped since lower
  is better), synergy best-partner (`network_synergy`'s `max(synergy_score)`
  across a compound's partners), and module robustness (`r_index`, joined
  via `network_module_membership`). Toxicity alerts and
  `bias_reweight()`'s `score_ajustado` were both considered as criteria and
  deliberately dropped in this version. Needs the `RobustRankAggreg`
  package (new `Suggests`) -- namespace-guarded, same pattern as
  `bipartite`/`dbscan`. Output keeps every raw `crit_*` value, its
  direction-normalised `rank_*`, `rra_score`/`rra_rank`, and a
  from-scratch, all-criteria Pareto front tier -- the same "never one
  opaque number" convention `network_synergy()` already established.
  `export = "sdf"`/`"smi"` writes the top-N ranked compounds' structures
  (reuses `rcdk`/`.parse_smiles_safe()`, no new dependency).
- **New `plot_rank()`.** `view = "pareto"` (default): a new plain-`ggplot2`
  scatter (no dependency), coloured by the Pareto front `rank_candidates()`
  already computed, front-1 (non-dominated) compounds labelled by name.
  `view = "heatmap"`: compound x target matrix restricted to the top
  `rra_rank`ed compounds and their most relevant targets (breadth/weight/
  centrality) -- delegates to the same `pheatmap` rendering
  `plot_heatmap()` uses (factored out into `.plot_pheatmap_render()` in
  this release) rather than duplicating the geometry.
- Design rationale, the RRA-vs-weighted-sum-vs-desirability-functions
  comparison, and every roll-up/missing-data decision:
  `chilcuague-analysis/reviews/rank-candidates-design-spec.md`.

## `network_build()` / `network_bowtie()` / `network_filter_proteome()` / `network_pathview()` / `network_motifs()` — Phase 5 audit fixes (design spec §1.1, §1.6, §1.10, §1.11, §1.12)

- **BREAKING: `network_build()` no longer emits `disease_association_score`.**
  The old join keyed on `(compound_id, uniprot_id)` against `targets_disease`
  via `match()`, which only ever returns the *first* hit -- once
  `targets_disease_filter()` had accumulated more than one disease (the
  normal case for `network_proximity()`/`network_synergy()`), every edge
  silently carried an arbitrary one of them with no `disease_id` column
  recording which. There is no fix that keeps a single scalar column
  meaningful, so the column is gone from `network_edges`,
  `network_filtered_edges`, and their `.csv` files. Join
  `patliRResults(proj, "targets_disease")` yourself on
  `(compound_id, target_id, disease_id)` if you need this information.
- **`network_build()`/`.network_graph()` no longer serve a stale cache
  behind an edited CSV.** The cached `.rds` graph is stamped with a row
  count and a cheap content checksum for the condition's edges; a later
  `.network_graph()` call rebuilds instead of trusting the cache when a
  hand-edited `results/network_edges.csv` (or an in-memory
  `patliRResults()` write) disagrees with what was stamped.
- **`network_build()` now asserts the built graph is bipartite** (via the
  existing `.network_is_bipartite()` helper) and aborts with a clear
  message naming the colliding ID(s) if a `compound_id` ever equals a
  `uniprot_id` -- that ID would otherwise be silently typed as a compound
  and its target edge would become intra-mode.
- **`.cache_key_slug()` no longer collides on punctuation-only differences**
  (e.g. `"FLO-ET"` and `"FLO_ET"` used to map to the same cache file). A
  short hash of the raw string is now appended to the slug.
- **New `network_bowtie(actions_score_threshold = 400)`.** Previously every
  directed row of STRING's actions file was used regardless of STRING's own
  per-interaction confidence, so the core SCC size (the headline bow-tie
  number) was driven by the lowest-confidence predictions. The threshold
  used is recorded in `network_bowtie_summary`. A `cli_warn` now also fires
  when the resulting core is degenerate (< 50 nodes, or < 1% of the
  network). Internally, UniProt -> STRING_id mapping now reads STRING's
  `protein.aliases` flat file directly instead of instantiating a second
  `STRINGdb` object at `score_threshold = 0` (which downloaded the entire
  interactome just to reach an alias table that never depended on score).
- **New `network_filter_proteome(proteome_label = NULL)`.** Output rows now
  carry `proteome_label` (auto-derived from a hash of the sorted proteome
  when not supplied) and `n_proteome`; `proteome_label` is part of the
  upsert key alongside `condition`, so two different proteomes filtered
  against the same condition are both kept instead of one overwriting the
  other. A `cli_warn` now fires when `proteome` matches nothing in
  `network_edges` and fewer than half its entries look like UniProt
  accessions (the usual cause: gene symbols passed by mistake).
- **`network_pathview()`'s `ok` column now means the PNG actually exists**,
  not just "pathview() did not raise an R error" -- `pathview::pathview()`
  frequently warns and returns without writing anything (no mappable nodes,
  or a KGML/PNG download that silently returns an HTML error page); that
  case now gets `ok = FALSE`/`path = NA`, where it previously got
  `ok = TRUE` and a `path` pointing at a file that was never created.
  Targets with an `NA` `weight` are now explicitly dropped (and logged) from
  the `"max_weight"`/`"mean_weight"` colour vector instead of vanishing via
  `aggregate()`'s implicit `na.action = na.omit`; if that would leave zero
  targets scored, the function now aborts with a clear message instead of
  handing `pathview()` an empty vector. A `cli_warn` fires if any score
  exceeds 1 before the `[0, 1]` colour clamp (a tripwire for a platform
  reporting raw 0-100 percentages instead of `targets_import()`'s `[0, 1]`
  convention).
- **`network_motifs()`'s log line now reports the deduplicated feedback
  *cycle* count, not the row count.** A directed 3-cycle is still listed 3
  times in `network_motifs` output (once per starting node -- the
  row-emission shape is unchanged), but the `projectLog()` summary used to
  report that row count as if it were the number of cycles found (3x the
  true count for every real cycle).

## `network_proximity()` / `plot_proximity()` — opt-in raw null draws + a permutation-histogram view

- **New `network_proximity(store_null = FALSE)`.** When `TRUE`, additionally
  persists every individual random draw's `d_random` value (already computed
  on the way to `d_random_mean`/`d_random_sd`/`z_score` -- nothing new is
  computed) to a new `network_proximity_null` results slot, one row per
  `(compound, disease, draw)`. Opt-in because it is large (~1 MB at
  `n_random = 1000` / 30 compounds). `store_null = FALSE` (default) leaves any
  previously stored draws untouched, and the main `network_proximity` output
  is byte-identical either way.
- **New `plot_proximity(view = c("z", "null"))`.** `"z"` is the existing
  z-score lollipop, unchanged. `"null"` draws the literal permutation-
  histogram figure the literature uses (Guney et al. 2016 Fig. 1; Menche et
  al. 2015 SI) -- one facet per compound, its raw `d_random` draws
  histogrammed, `d_observed` overlaid, `z_score`/`p_adjusted` annotated.
  Needs `network_proximity(store_null = TRUE)`; aborts with a clear message
  naming it otherwise. New `top_n` argument (default `12`, `view = "null"`
  only) caps the number of facets, ranked by `|z_score|`.

## `plot_target_chord()` — new, arc diagram of target-target STRING actions

- **New `plot_target_chord()`** (spec 3.9, lowest priority of the Phase-4
  plot roster). An arc diagram (Krzywinski et al. 2009's idiom, plain
  `ggplot2` -- no `circlize` dependency) of a condition's targets placed
  along one line, with an arc above for every STRING *actions* interaction
  between two plotted targets, reusing [network_bowtie()]'s already-
  downloaded/cached actions graph rather than fetching anything new.
  `actions_score_threshold` (default `400`, matching the rest of the
  package's STRING confidence convention) recolours arcs by edge score when
  a per-edge score can be recovered from the cached raw actions file, and
  falls back to one constant arc colour otherwise. `top_n_labels` (default
  `15`) labels only the highest-degree targets; the rest get an unlabeled
  tick. Needs `STRINGdb`; aborts cleanly on fewer than 2 targets or zero
  actions-edges among the plotted set.

## `plot_network_layers()` / `plot_network_degeneracy()` — bipartite layout, module colouring, significance-based degeneracy filter

- **`plot_network_layers()` gains `layout = "bipartite"`** — a
  deterministic two-column layout (compounds on the left, targets on the
  right, each column ordered by degree so fewer edges cross; the standard
  affiliation-network rendering, Borgatti & Everett 1997). Only valid on a
  genuinely two-layer (compound-target) scope — `cli_abort`s naming the
  extra layer(s) if `layers` also includes `"pathway"`/`"disease"`. The
  existing `"fr"`/`"kk"`/`"drl"` force-directed layouts are unchanged
  (`"fr"` stays the default).
- **`plot_network_layers()` gains `colour_by = c("layer", "module",
  "node_type")`** (default `"layer"`, the historical behaviour, unchanged).
  `"module"` joins [network_module_robustness()]'s per-node
  `network_module_membership` table by `(condition, node_id)` and colours
  with a qualitative palette (`grDevices::hcl.colors(n, "Dark 3")`), with a
  grey `"not clustered"` level for any node absent from that table;
  `cli_abort`s if the table does not exist at all, naming
  `network_module_robustness()`. `"node_type"` colours by compound vs.
  target.
- **Node colour is now carried on `fill`, not `colour`** (a plain
  border replaces the old borderless point) — this is what lets
  `plot_network_degeneracy()`'s z-score-coloured arcs use their own
  independent colour scale without clashing with the node legend.
- **`network_layers_plot_log`'s upsert key is now `(condition, layout,
  colour_by)`** (was `condition` alone), and the saved filename encodes all
  three (e.g. `network_layers_ALL_fr_layer.png`) — different
  `layout`/`colour_by` views no longer overwrite each other's log row or
  PNG. **Breaking for anyone reading the PNG path by its old, condition-only
  filename.**
- **`plot_network_degeneracy()`'s default filter changes from an arbitrary
  score cutoff to significance** — **breaking**: the new `filter =
  c("p_adjusted", "score")` argument defaults to `"p_adjusted"`, drawing
  every compound pair with `p_adjusted < alpha` (`alpha` new, default
  `0.05`) instead of the old always-on `degeneracy_score >= min_degeneracy`
  (`min_degeneracy` default unchanged at `0.3`). `filter = "score"`
  reproduces the pre-revision behaviour exactly. `p_adjusted` is `NA` for
  every row when [network_degeneracy()] was run with `annotation %in%
  c("enriched", "jaccard")` (no permutation null exists for those modes,
  since piece 14) — in that case `filter = "p_adjusted"` (including the
  default) transparently falls back to `filter = "score"` and
  `cli_inform`s why.
- **Degeneracy arcs are now coloured by `z_score`** (bounded,
  significance-bearing) instead of sized by the unbounded
  `degeneracy_score`.
- **`plot_network_degeneracy()` inherits `layout`/`colour_by`** from
  `plot_network_layers()` (its scope is always two-layer, so `layout =
  "bipartite"` is always available; `colour_by = "module"` lets a reader
  see whether a degenerate pair's two compounds landed in the same module).
  Its own `network_degeneracy_plot_log` key/filename gains `layout`,
  `colour_by` and `filter` for the same reason as above.

## `network_synergy()` — real `s_AB` (Menche 2015) + Cheng P1–P6 classification

- **New `separation = c("network", "jaccard")`** (default `"network"`) —
  **breaking**: the default now requires `STRINGdb` and a one-time
  interactome download, and `cli_abort`s (naming `separation = "jaccard"`
  as the escape) when `STRINGdb` is absent, where the old behaviour was
  always the cheap set-overlap proxy. `"network"` computes the
  topological separation
  `s_AB = ⟨d_AB⟩ − (⟨d_AA⟩ + ⟨d_BB⟩)/2` (Menche et al. 2015) between the
  two compounds' STRING-mapped target sets on the interactome's largest
  connected component, and classifies every pair into Cheng, Kovács &
  Barabási (2019)'s `P1`–`P6` (`cheng_class`); `complementary_exposure`
  is `cheng_class == "P2"`. Absolute `s_AB` and the `P1`–`P6` boundary
  are not numerically comparable to Cheng's published figures (STRING at
  `score_threshold = 400` includes text-mining / co-expression edges);
  the classification is for relative ranking within a run. `"jaccard"`
  keeps the old set-overlap proxy for users without `STRINGdb` (all
  `s_AB`/Cheng columns `NA`, `synergy_score` stays on its historical
  `both_proximal` gate).
- **New arguments** `alpha` (proximity significance cut for the
  `proximal_a`/`proximal_b` predicates — patliR's tightening, *not*
  Cheng's, who use the z-score alone), `species` / `version` /
  `score_threshold` (must match the `network_proximity()` run — a
  mismatch now aborts), and `disease_gene_source` (default
  `"disease_genes"`; consuming the legacy circular `"targets_disease"`
  rows warns).
- **No permutation null for `s_AB`.** Cheng et al. explicitly reject a
  z-score for the drug–drug relationship (≈3 targets per drug → the
  randomisation is non-Gaussian); `s_AB` is used raw.
- **Singleton convention.** A STRING-mapped target set of size < 2 gets
  `⟨d_AA⟩ = 0` (patliR's choice; Menche's `separation.py` returns `nan`)
  and a `singleton_a` / `singleton_b` flag. This inflates `s_AB` toward
  `separated = TRUE`; flagged pairs keep a visible `cheng_class` but are
  excluded from the `synergy_score` scalar.
- **`cli_warn` when the proximity gate is arithmetically unreachable** —
  `n_tests_in_family / (n_random + 1) > alpha` means the BH-adjusted
  proximity p can never clear `alpha` and every pair silently falls into
  `P5`/`P6`.
- **Wider output** (pure column addition — `.network_upsert()` migrates a
  legacy `network_synergy.csv` in place): `n_targets_a`/`_b`,
  `n_targets_a_mapped`/`_b_mapped`, `singleton_a`/`_b`, `d_aa`, `d_bb`,
  `d_ab`, `s_ab`, `separated`, `p_adjusted_a`/`_b`, `proximal_a`/`_b`,
  `cheng_class`, `complementary_exposure`, `separation_method`, `alpha`,
  `species`, `string_version`, `disease_gene_source`. `target_jaccard` /
  `complementarity` are unchanged and still computed on the **raw
  UniProt** sets in both modes.
- **`network_proximity()`** gains `species`, `string_version`,
  `score_threshold` and `n_tests_in_family` columns (the STRING
  interactome and BH family size the run used), so `network_synergy()`
  can refuse to combine a `z` and an `s_AB` from different interactomes.
  Also pure column additions.
- The STRING LCC + degree-binning block is factored into
  `.network_string_lcc()`, shared by `network_proximity()` and
  `network_synergy()` so both operate on a byte-identical graph. Only the
  graph is cached (`cacheDir(proj)/stringdb/lcc_<species>_<version>_<thr>.rds`,
  safe to delete); `network_proximity()`'s `z_score` is unchanged.
- `network_synergy()` now `cli_warn`s when the `network_proximity()` rows
  it consumes were assembled from calls with different BH family sizes
  (`n_tests_in_family`), and records `disease_gene_source = NA` (with a
  warning) rather than fabricating `"disease_genes"` when the proximity
  table predates the provenance column.
- `plot_synergy()` draws a self-explanatory empty panel (instead of an
  all-points-dropped ggplot warning) when no pair is scored — the
  guaranteed outcome whenever the proximity BH gate is unreachable.
- CSV round-trip: results columns that are character *by contract* but
  hold numeric-looking values (`string_version` = `"12.0"`, and
  `condition`) are pinned via `colClasses` when `patliR_load()` /
  `.network_graph()` read them, so a save→load→rerun no longer trips
  `.network_upsert()`'s (correctly strict) column-type check. The
  short-lived character-vs-numeric tolerance in `.col_compat()` is
  reverted.

## `network_degeneracy()` — GO semantic similarity + a permutation null

- **New `annotation = c("direct", "enriched", "jaccard")`** (default
  `"direct"`). `"direct"` scores functional convergence as **GO semantic
  similarity** (`GOSemSim`, Wang/BMA — `GOSemSim::clusterSim()`-equivalent,
  reproduced exactly by a fast pooled-`termSim`-matrix path) between the two
  compounds' Entrez-mapped targets on each protein's *own* GO annotations.
  It needs only `network_build()` — **not `network_enrich()`** — which
  removes the enrichment-universe circularity the old `pathway_jaccard`
  carried. `"enriched"` is `GOSemSim::mgoSim()` over the enriched-term sets
  (keeps `network_enrich()` as a prerequisite and the circularity);
  `"jaccard"` is the historical enriched-pathway Jaccard (no `GOSemSim`).
- **Permutation null (`"direct"` mode).** The observed similarity is
  z-scored against `n_random` (default `200`) resamplings drawn from
  **GO-annotation-count bins** (the same `.network_value_bins()` consecutive
  binning `network_proximity()` uses for degree) over
  `universe = c("project", "condition", "genome")` (default `"project"`).
  `z_score` is the primary statistic; `p_empirical` / `p_adjusted` are the
  **right-tail** permutation p (opposite tail to `network_proximity()`'s
  identically-named left-tail `p_adjusted`), BH per call, and
  resolution-limited (a warning fires when `n_pairs / (n_random + 1) >
  0.05`). `n_random = 1` gives `sim_random_sd`/`z_score` `NA`.
- **`pathway_jaccard` renamed to `functional_similarity`.** An existing
  `results/network_degeneracy.csv` is migrated in place on the first run
  (the column is renamed before the upsert, and the new provenance columns
  back-filled `NA`) — no delete-and-rebuild. `degeneracy_score` is retained.
- **New arguments** `ont` (`"BP"` default), `measure` (`"Wang"` default; the
  IC measures need a separate `computeIC = TRUE` `godata`), `combine`
  (`"BMA"`), `drop` (`"IEA"` — excludes electronic annotations, ~44% of
  `org.Hs.eg.db` BP), `seed`. New output columns `n_genes_a_mapped` /
  `n_genes_b_mapped`, `sim_random_mean` / `sim_random_sd` / `z_score` /
  `p_empirical` / `p_adjusted`, and the per-call constants `annotation` /
  `ont` / `measure` / `combine` / `drop` / `universe` / `n_random` /
  `seed_used`. `n_pathways_a/b` are retained but `NA` in `"direct"` mode.
- **`godata` construction is now shared** between `network_enrich()` and
  `network_degeneracy()` via `.network_godata()` (uses the non-deprecated
  `annoDb =` spelling; caches under `cacheDir(proj)/gosemsim/` with the
  `org.Hs.eg.db` / `GO.db` versions in the key).
- The proximity resampling helpers `.network_degree_bins()` /
  `.network_resample_degree_matched()` are renamed
  `.network_value_bins()` / `.network_resample_matched()` (generalised
  from degree to any per-node numeric property; signatures unchanged).
  `network_proximity()` output is byte-identical.
- **Fixed: `universe = "genome"` could silently score zero pairs.** The
  random 1500-gene cap on the genome-wide resampling pool was applied
  *before* checking which genes a compound actually targets, so it could
  (and on the package's own fixture, did) evict a compound's own targets
  from the survivors -- every pair involving that compound then skipped,
  with a log message blaming a missing annotation rather than the cap.
  `.network_cap_pool_keep_observed()` now guarantees every observed gene
  survives the cap; only the remainder of the pool is randomly subsampled.
- `network_degeneracy(annotation = "enriched")` now warns (rather than
  silently scoring every pair `NA`) when the condition's enriched terms
  contain no GO IDs -- `functional_similarity` in that mode needs
  `GOSemSim::mgoSim()` over GO terms specifically.
- `drop` is recorded as `NA` for `annotation != "direct"` rows: the
  IEA-annotation filter only affects the `"direct"`-mode annotation index,
  not `"enriched"`'s `mgoSim()` call, so stamping `"IEA"` on an `"enriched"`
  row claimed a filter that had no effect on that row's number.
- `network_degeneracy()` now `cli_inform()`s, once per call, that
  `sim_random_*`/`z_score`/`p_*` are `NA` for every row when
  `annotation != "direct"` (no permutation null exists for `"enriched"`/
  `"jaccard"` -- see the `@section` on why).

## `network_module_robustness()` — clustering backends and a random-failure baseline

- **Breaking change to the numbers.** `network_module_robustness()`'s
  `clustering` argument now defaults to **`"leiden"`** (was `"hdbscan"`).
  Modules are detected by `igraph::cluster_leiden()` maximising Newman
  modularity, so every module partition — and therefore every `r_index` —
  differs from earlier versions. `clustering = "hdbscan"` still works and
  emits a one-time note that it is a density heuristic on an integer
  metric; it is not deprecated. Its module partitions and every `r_index`
  are unchanged for the same seed, but `module_type` on a 3-node component
  now reflects the backend's own verdict (`"noise"` rather than the old
  method-agnostic `"cluster"`): the single-module fallback threshold moved
  from `< 4` nodes to `< min_component_size` (default 3), so a 3-node
  component is now clustered rather than passed through.
- **New `clustering = "bipartite"`** — `bipartite::computeModules(method =
  "Beckett")`, maximising Barber (2007)'s bipartite modularity `Q_B`. Only
  valid on the genuine two-mode compound–target graph (aborts otherwise).
  A star / single-mode component within the graph falls back to a single
  module (`computeModules()` errors on a `1×k` biadjacency matrix). Adds
  `bipartite` to `Suggests` — loaded namespace-only, so it does **not**
  attach `sna`/`vegan` or mask any `igraph` function.
- Clustering is run **unweighted**: for Leiden the `weight` edge attribute
  is *deleted* before clustering rather than passed as `weights = NA`
  (which would leave `strength()` reading it for the modularity null model
  and silently return all-singletons on any `NA` probability).
- **New `attack = c("targeted", "random", "both")`** (default
  `"targeted"`). `"random"` averages `n_random` (default 20; **must be
  `>= 2`** when `attack` includes `"random"`) uniformly random removal
  orders under overflow-safe
  per-`(module, replicate)` sub-seeds; `network_robustness_curve` gains
  `removal_strategy`, `largest_component_sd`, `n_replicates`. The summary
  table gains `r_index_random` and `n_random` alongside the unchanged
  `r_index` (targeted).
- **New third output slot `network_module_membership`** — one row per node
  (`condition`, `node_id`, `node_type`, `module_id`, `module_type`), to be
  consumed by a future `plot_network_layers()` module-colour mode.
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
  once). It also adds `hub_score_component` (eigenvector centrality
  recomputed per connected component, each scaled to `max = 1` -- fixes the
  disconnected-graph case where every node outside the largest component
  got `hub_score` ~ 0) with a `component_id` label, and the per-condition
  constants `n_compounds` / `n_targets`. Edge cases: an empty mode gives
  `_norm = NA` (never `NaN`/`0`); `B_max == 0` (a mode of size 1) gives
  `NA` with a log line; an isolated node gives `hub_score_component = NA`.
  `normalize = FALSE` now emits the **same columns** as `normalize = TRUE`
  (the `_norm` / `component_id` / count columns just come out `NA`), so it
  no longer aborts `.network_upsert()` when run against a project that
  already holds a normalised table. All columns are declared in the
  empty-row constructor, so `results/network_centrality.csv` stays
  column-stable across zero-row reruns; `.network_upsert()` back-fills them
  on a project whose CSV predates this change. `_norm` values are **not**
  comparable across conditions of different size (documented like the
  R-index).
- **`hub_score` is now `igraph::eigen_centrality()`, not
  `igraph::hits_scores()$hub`.** On a two-mode graph the HITS hub is the
  principal eigenvector of a matrix with a degenerate top eigenspace, from
  which ARPACK returns an RNG-seeded arbitrary vector -- so `hub_score` (and
  the new `hub_score_component`) took different values on every run,
  breaking the pure-function contract. The Perron eigenvector is unique on
  each connected component, so `eigen_centrality()` is deterministic; on
  an undirected graph it is exactly the quantity the roxygen already
  claimed `hub_score` to be.
- **`plot_centrality()`**: `measure` accepts `"degree_norm"` /
  `"betweenness_norm"`; when `measure` is not supplied and
  `node_type = "both"` it now defaults to `"degree_norm"` (with a
  `cli_inform` naming the reason), falling back to `"degree"` with a
  warning when no `degree_norm` value is available for the conditions in
  scope (a pre-0.2.0 result, or one built with `normalize = FALSE`). The
  value axis states
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

## `plot_synergy()` — replaced with the Cheng Complementary-Exposure quadrant; new `plot_enrichment()`

- **`plot_synergy()` is replaced, not revised.** The old figure plotted
  `complementarity` x `joint_closeness` (a patliR invention) with
  `synergy_score` on colour/size. The new figure plots `z_score_a` x
  `z_score_b` (Cheng, Kovacs & Barabasi (2019) Fig. 2's own axes, re-ordered
  per-row so the more disease-proximal compound is always on `x` --
  display-only, `compound_a`/`compound_b` in
  `patliRResults(proj, "network_synergy")` are untouched), colour =
  `separated` (`s_AB >= 0`, a two-level manual scale, not a gradient),
  shape = `cheng_class` (`P1`-`P6` + `NA`), size = `abs(s_ab)` (`NA` ->
  smallest, never dropped). The `P2` (Complementary Exposure) region is
  shaded and annotated with its count per panel; the top-`synergy_score`
  `P2` pairs are labeled by compound name. **Breaking**: `top_n`'s default
  drops from `10` to `5` (the new label set is `P2`-only and much smaller).
  The empty-state guard now fires on "no row has a `cheng_class`" rather
  than "no row is scored" -- a pair can be classified without ever getting
  a `synergy_score`, and gating on the latter made the guard fire on almost
  every real run.
- **New companion figure**: `synergy_sab_<scope>.png`, a histogram of
  `s_AB` across every pair in scope with a `s_AB = 0` reference line --
  Cheng et al.'s Fig. 1 separation distribution. Logged separately to
  `patliRResults(proj, "synergy_sab_plot_log")`.
- **New `plot_enrichment()`** -- the `clusterProfiler::dotplot()`-style
  enrichment bubble that `network_enrich()` (a core function) previously
  had no matching plot for. Term (`Description`) on `y` ordered by
  `GeneRatio` (parsed from its `"k/n"` string, e.g. `"3/20"` -> `0.15`),
  `GeneRatio` on `x`, size = `Count`, colour = `p.adjust`
  (`scale_colour_gradient(low = "#c0392b", high = "#2980b9")`, matching the
  rest of the package's "red = the thing worth looking at" convention: low
  `p.adjust` = more significant = red). Selects the `top_n` (default `20`)
  most-significant terms per `(condition, db)`, facets by `db` when more
  than one is in scope. Does not require the `qvalue` column (absent
  whenever `clusterProfiler`'s pi0 estimation fails, routine for a small
  gene set). Logged to `patliRResults(proj, "enrichment_plot_log")`.

## `network_enrich()` — background universe (BREAKING), `pAdjustMethod`, and crash fixes

- **BREAKING: new default background universe.** `network_enrich()` gains a
  `universe = c("project", "genome")` argument, **defaulting to
  `"project"`**. Previously every call to `enrichGO()`/`enrichKEGG()`/
  `enrichPathway()` ran with no `universe =` at all, so the background was
  clusterProfiler's own whole-genome/whole-pathway-database default (all
  ~20,000 `org.Hs.eg.db` genes for GO, the whole KEGG/Reactome gene-set
  universe otherwise) -- appropriate for an unbiased gene list, but the
  foreground here is the reachable output of a ligand-similarity target
  predictor (SuperPred/SwissTargetPrediction), structurally skewed toward
  GPCRs, kinases, nuclear receptors, proteases, and transporters -- the
  only families with enough known ligands to build a similarity model from
  in the first place. Tested against a whole-genome background, those
  families come out "enriched" for *any* input compound before a single
  biological difference between compounds is considered (Timmons, Szkop &
  Gallagher 2015, *Genome Biol* 16:186, "Multiple sources of bias confound
  functional enrichment analysis of global -omics data"; Boyle et al. 2004,
  *Bioinformatics* 20(18):3710-3715, the original over-representation-test
  framing). `universe = "project"` restricts the background to the
  project's own predictable proteome --
  `unique(patliRResults(proj, "targets_imported")$uniprot_id)`, or every
  target across every built condition if that slot is absent -- mapped to
  Entrez once per call (not once per condition, since the pool does not
  depend on which condition is being enriched), the same pool and mapping
  `network_degeneracy(universe = "project")` already uses so the two
  concepts mean the same thing across the package. `universe = "genome"`
  restores the old, unrestricted behaviour as an explicit opt-out (for
  comparison, or when there is a specific reason to want it). The resolved
  universe is recorded per row in a new `universe` column.
  **This changes every p-value `network_enrich()` reports for every
  existing project that re-runs with the new default** -- and therefore
  every downstream `network_degeneracy(annotation = "enriched"/"jaccard")`,
  `network_motifs()`'s pathway layer, and `plot_network_layers()`/
  `plot_enrichment()` result that consumes `network_enrichment`. Re-run
  `network_enrich()` (and anything downstream of it) deliberately after
  upgrading -- do not assume the numbers carried over.
- **Fix: a missing `qvalue` column no longer crashes the run.**
  `clusterProfiler` omits `qvalue` when the `qvalue` package's pi0
  estimation fails -- routine on the small p-value vectors a 10-50-gene
  target set produces. The result table is now built column-by-column,
  with `NA` fill for `qvalue` (and for the new `ONTOLOGY` column, populated
  only when `ont = "ALL"`) instead of a hard `df[, c(...)]` subset that
  threw `undefined columns selected` and aborted a run that had already
  spent minutes on GO/KEGG/Reactome.
- **New `pAdjustMethod` argument** (default `"BH"`, matching
  `clusterProfiler`'s own internal default -- not previously exposed at
  all), passed through to `enrichGO()`/`enrichKEGG()`/`enrichPathway()`.
- **Internal: switched the `(condition, db)` replace-on-rerun logic to the
  shared `.network_upsert()`**, instead of a hand-rolled re-derivation --
  one code path, one set of tests, and it fixes a blind spot in the old
  logic (a condition whose enrichment came back empty never appeared in
  the touched-keys it derived, so a stale row from an earlier, non-empty
  run of that same `(condition, db)` could survive an empty rerun
  untouched). Legacy `network_enrichment.csv` files predating the
  `universe`/`ONTOLOGY` columns migrate through it cleanly (pure column
  addition, back-filled as `NA` on existing rows).
- **`network_degeneracy()`'s abort message**, shown when `annotation =
  "enriched"`/`"jaccard"` needs a `network_enrichment` entry that isn't
  there, no longer flatly says "run `network_enrich()` first" -- an empty
  or missing `network_enrichment` table is indistinguishable from "ran and
  found nothing" (by design: `patliR` never writes a placeholder row for
  "nothing to report"), so the message now also names that possibility and
  points at `projectLog()`.

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
