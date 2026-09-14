# What each function actually computes

A plain-language pass over every exported function: what goes in, what comes
out, **which operation the code performs** (does it average? filter? fit a
model? run a permutation test?), and the theory behind it. Written for
someone who knows the pharmacology/network science but does not read R.

Notation: `proj` is the project object every step takes and returns; the
real output is always a CSV in the project folder.

---

## Project object

| function | what it does |
|---|---|
| `patliR_project(dir)` | Creates an empty project rooted at a folder. No computation. |
| `patliR_load(dir)` | Rebuilds a project by **reading back the CSV files** a previous run wrote. No re-computation — the CSVs are the source of truth. |
| `compounds()`, `matrixRaw()`, `binarizedMatrix()`, `patliRResults()`, `projectLog()`, `projectDir()`, `cacheDir()` | Accessors. Return a stored table; compute nothing. |
| `patliR_export_llm(proj)` | Concatenates every table in the project into one plain-text file. String formatting only. |

---

## `prep_*` — getting compounds in

### `prep_compounds(proj, table, identifier)`
**In:** a table of name / CAS / PubChem CID / SMILES. **Out:** `compounds()` — one row per valid compound with an internal id and a canonical SMILES.
**Operation:** *validation and canonicalisation, no modelling.* Each SMILES is parsed by the CDK toolkit (`rcdk`); if it parses to a non-empty structure the compound is kept, otherwise it is logged and dropped. The stored `canonical_smiles` is CDK's canonical form (`Canonical` + aromatic flavour) — a deterministic re-writing of the same molecule, **not** an InChIKey. If a row has a CID but no SMILES, the SMILES is fetched from PubChem.
**Theory:** none — this is data hygiene. The canonical-SMILES choice is a compromise (`rcdk` has no InChIKey); see `DESIGN.md`.

### `prep_compound(proj, one_row, ...)`
Thin wrapper: same as `prep_compounds()` for a single row.

### `prep_binarize(proj, abundance_matrix)`
**In:** a GC-MS abundance matrix, columns named `R<n>-<CONDITION>` (replicate, condition). **Out:** `matrixRaw()` (per-condition averages) and `binarizedMatrix()` (0/1 per condition).
**Operation:** *average, then threshold.* For each condition the replicate columns are **averaged** (arithmetic mean, ignoring NAs). Then a **per-condition first quartile (Q1, the 25th percentile)** is computed over the non-zero averaged values, and a compound is called *present* (1) in that condition if its average is `> 0` **and** `>= Q1`, else *absent* (0). Columns that are already 0/1 are passed through with a majority rule instead.
**Theory:** the Q1 cut is a simple "drop the bottom quarter" detection threshold; no distributional model.

### `prep_as_condition(proj, condition, compound_ids)`
**In:** a plain compound list (no matrix). **Out:** `binarizedMatrix()` with one column, all `1` for the chosen compounds.
**Operation:** *labelling only.* Marks the whole list (or a subset) as one experimental condition so the network step has something to work on when there is no abundance matrix.

### `prep_structure2d(proj)`
**Operation:** renders one 2D structure image per compound (`rcdk` depiction) and records whether each rendered. Drawing, not analysis.

---

## `refdb_*` — the local reference database

### `refdb_build(proj, sources = c("pubchem", "chembl"))`
**In:** the compound list. **Out:** `reference_compounds.csv` (identity: which PubChem / ChEMBL entry each compound is) and `reference_bioactivity.csv` (measured activities pulled from ChEMBL).
**Operation:** *look-up and a deterministic identity chain, no modelling.* For each compound it queries PubChem (by CID). ChEMBL identity is resolved by a chain, never a similarity search: PubChem CID → standard InChIKey (batched) → ChEMBL `standard_inchi_key` exact match (batched) → InChIKey-skeleton (14-character) fallback → SMILES `flexmatch` → give up — each step only runs if the previous one found nothing, and every row records which `match_type` it landed on. A compound with no ChEMBL entry at all (common for natural-product metabolites) is logged, not an error. Re-running replaces that compound's rows rather than duplicating them.
**Theory:** none. This table is the "background" that `bias_audit()` and `network_degeneracy()` later treat as "the known universe" for these compounds.
**Fixed (commit `8cdfe36`):** ChEMBL matching used to run a 100%-Tanimoto *similarity* search, which is not the same as identity and could pick the wrong molecule for long chains or stereoisomers. The identity chain above replaced it; `match_type` on each `reference_compounds` row tells you how confident that particular match is (exact InChIKey is strongest, `smiles_flexmatch` weakest before giving up).

### `refdb_update(proj, compound_ids)`
Same as `refdb_build()` but only for the named compounds — use after adding compounds.

### `refdb_rebuild_cache(proj)`
Regenerates a fast binary copy of the reference CSVs. Pure caching; always safe to delete.

---

## `compounds_classify()` / `compounds_similarity()`

### `compounds_classify(proj)`
**In:** SMILES. **Out:** a natural-product *pathway / superclass / class* per compound (e.g. "Terpenoids", "Shikimates and Phenylpropanoids").
**Operation:** *external classifier call, no local computation.* Sends each SMILES to the NPClassifier web service (one request per second, with back-off) and stores the returned labels.
**Theory:** NPClassifier (Kim et al. 2021, *J. Nat. Prod.* 84:2795) is a neural classifier trained on a biosynthesis-oriented taxonomy of natural products. Chosen over the more general ClassyFire because this package's domain is plant extracts.

### `compounds_similarity(proj, method = "tanimoto")`
**In:** SMILES. **Out:** a long table of every compound pair with a similarity score.
**Operation:** *fingerprint + pairwise set similarity.* Computes a molecular fingerprint per compound (`rcdk::get.fingerprint()` — a bit-vector of substructure features), then for every pair computes **Tanimoto** (Jaccard on the bit sets), Dice, or cosine similarity.
**Theory:** Tanimoto similarity of 2D fingerprints — the standard cheminformatics measure of structural resemblance. It says nothing about shared targets or activity.

---

## `adme_*` — physicochemistry and drug-likeness

### `adme_local(proj)`
**In:** SMILES. **Out:** `adme_local.csv` — molecular descriptors plus five binary drug-/lead-likeness flags and a gut/brain-absorption call.
**Operation:** *descriptor calculation + fixed-threshold rules, no fitting.* Descriptors (MW, logP, H-bond donors/acceptors, TPSA, rotatable bonds, molar refractivity, atom count) are read straight from CDK (`rcdk`). Each rule is then a hard inequality on those descriptors:

| flag | rule (all must hold) | source |
|---|---|---|
| `ro5_pass` | MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 | Lipinski et al. 2001 |
| `veber_pass` | TPSA ≤ 140, rotatable bonds ≤ 10 | Veber et al. 2002 |
| `ghose_pass` | 160 ≤ MW ≤ 480, −0.4 ≤ logP ≤ 5.6, 20 ≤ atoms ≤ 70, 40 ≤ AMR ≤ 130 | Ghose et al. 1999 |
| `egan_pass` | logP ≤ 5.88, TPSA ≤ 131.6 | Egan et al. 2000 |
| `oprea_pass` | HBD < 2, 2 < HBA < 10, 2 < rot. bonds < 8, 1 < rings < 4 | Oprea 2000 (lead-likeness) |

`gi_absorption` / `bbb_permeant` come from a **point-in-ellipse test**: the compound's (TPSA, WLogP) point is checked against the two published BOILED-Egg ellipses (Daina & Zoete 2016) by ray-casting. Optional route flags (`route_oral` etc.) are further fixed inequalities.
**Theory:** all five rules are published "rule-of-thumb" boundaries of the physicochemical region where oral drugs tend to fall. They are heuristics, not classifiers — no training, no probability.

### `adme_filter(proj, rules)`
**In:** `adme_local()` results. **Out:** a long pass/fail table, and (only if you ask) a shorter compound list.
**Operation:** *filtering.* By default it just **records** which rules each compound passes. With `hard_cutoff = TRUE` it removes compounds that fail any selected rule — after showing you the count and (interactively) asking. Never silent.

### `adme_import(proj, csv, platform)`
**In:** a CSV exported from SwissADME / ADMETlab / pkCSM. **Out:** those properties in long form, matched to your compounds.
**Operation:** *column mapping and a join.* Renames the platform's columns to patliR's names and matches rows to compounds by canonical SMILES (or PubChem CID). No computation on the values.

### `adme_export_smiles(proj)`
Writes a plain SMILES list (one per line) plus a row-order "bridge" CSV, so a platform's export can be matched back by position. String output only.

---

## `tox_*` — toxicity / interference alerts

### `tox_local(proj, alert_sets = c("pains", "brenk"))`
**In:** SMILES. **Out:** `tox_local.csv` — one row per (compound, alert), with a `matched` flag.
**Operation:** *substructure pattern matching.* Each compound's structure is tested against every SMARTS pattern in the bundled PAINS (480) and Brenk (105) lists using CDK's SMARTS matcher (`rcdk::matches()`). A hit means the pattern is present somewhere in the molecule.
**Theory:** PAINS (Baell & Holloway 2010) are substructures statistically linked to *assay interference*; Brenk (2008) are substructures flagged as generally undesirable in screening libraries. **A hit is not a verdict** — many are also present in genuinely active chemotypes. The function deliberately has no cutoff.

### `tox_safetyome(proj)`
**In:** predicted targets (`targets_imported`). **Out:** for each target, whether it is one of ~500 "safety-critical" genes and which organ system it is tied to.
**Operation:** *ID mapping + table lookup.* Maps each target's UniProt ID to a gene symbol (`clusterProfiler::bitr` + `org.Hs.eg.db`), then checks membership in the bundled Safetyome core panel.
**Theory:** the panel is from Liu et al. (2026, *Tox. Sci.* 209(3)) — genes whose modulation is associated with systemic/organ toxicity. Again descriptive, not a pass/fail.

### `tox_import(proj, csv, platform, mapping_file)`
Same idea as `adme_import()` for hERG / DILI / Ames predictions. With `mapping_file` it matches by row position against the bridge CSV instead of by SMILES.

### `tox_export_smiles(proj)`
The toxicity-side counterpart of `adme_export_smiles()`.

### `tox_report(proj)`
**In:** whatever `tox_local` / `tox_safetyome` / `tox_import` produced. **Out:** one row per compound — PAINS count, Brenk count, safety-panel hits, imported properties — plus a fixed disclaimer, and `results/tox_report.csv`.
**Operation:** *aggregation only.* Counts and concatenates; runs no model and never collapses to a single score.

---

## `targets_*` — predicted protein targets

### `targets_import(proj, csv, platform)` / `targets_import_batch(proj, folder, platform)`
**In:** target-prediction exports (SuperPred / SwissTargetPrediction), one file per compound. **Out:** `targets_imported.csv` — compound, UniProt target, prediction probability.
**Operation:** *parsing and matching.* Reads the platform's columns, tolerant of small header differences, and matches each file to a compound by the PubChem CID in its file name. Probabilities given as `"96.5%"` strings are converted to `0.965`. No prediction is done here — patliR does not run a target model.
**Theory:** the upstream models (SuperPred etc.) use 2D/3D similarity to known ligands; patliR just ingests their output.

### `targets_disease_filter(proj, disease)`
**In:** the imported targets + a disease name. **Out:** `targets_disease.csv` — each target annotated with its association score to that one disease.
**Operation:** *external database query.* Resolves the disease name to an ontology ID and queries the Open Targets GraphQL API for each target's overall association score, optionally dropping targets below a score threshold.
**Theory:** Open Targets aggregates genetic, expression, pathway and literature evidence into a single 0–1 target–disease association score.
**Note:** this is a *target → disease* lookup — it only ever annotates UniProt IDs that are already predicted targets, so its output is a subset of the compounds' own targets. It is **not** a disease gene set for `network_proximity()`; use `disease_genes_*` for that.

---

## `disease_genes_*` — an independent disease gene set

`network_proximity()` needs a set of genes associated with the disease, assembled **without reference to the compounds under study**. Deriving it from `targets_disease_filter()` (above) makes the proximity z-score circular — the "disease module" would be a subset of the compounds' predicted targets, so the statistic would measure set membership, not topology. These two functions build the set independently and store it in a `disease_genes` slot (`results/disease_genes.csv`: `disease_id`, `disease_name`, `uniprot_id`, `gene_symbol`, `association_score`, `source`, `fetched_at`), upserted per `disease_id`.

### `disease_genes_fetch(proj, disease, min_score = 0.4)`
**In:** a disease name or ontology ID. **Out:** the `disease_genes` slot.
**Operation:** *external database query.* Resolves the disease to a current EFO/MONDO ID, then queries Open Targets in the **disease → target** direction (`disease(efoId:){ associatedTargets }`), paginating on `count`. Each associated target's `proteinIds` are filtered to `source == "uniprot_swissprot"` and mapped to UniProt accessions (a target with several Swiss-Prot accessions contributes one row each; one with none is excluded and counted in a summary log line). `min_score` drops targets whose aggregated Open Targets `association_score` is below the cut; kept/dropped counts are logged.
**Theory — and why `min_score` defaults to `0.4`.** Menche et al. (2015) and Guney et al. (2016) build the disease module from a curated list (OMIM + GWAS Catalog) of tens to low hundreds of genes. Open Targets' full associated-target list for a common disease is far larger — type 2 diabetes returns ~9,900 targets — and against the ~17k-node STRING largest connected component a module that big drags `d_observed` onto ≈ 0.5 by chance alone (about half of any compound's targets land inside `T` at random), which kills the z-score's dynamic range. `min_score = 0.4` brings common diseases to ~200–300 genes, a Menche-scale module. The Open Targets score aggregates several evidence types (genetic association, known drug, literature, …), so `0.4` is a pragmatic knob, not a principled biological threshold — tune it, or pass `min_score = NULL` to keep everything and filter `disease_genes` yourself. The aggregated score also folds in the *known-drug* channel, so a no-threshold module contains the targets of every drug approved for the disease — a second circularity if the compounds studied are themselves known drugs (not an issue for plant secondary metabolites).

### `disease_genes_import(proj, table, disease_id, map_symbols = FALSE)`
**In:** a curated list (data frame or CSV path) from OMIM / GWAS Catalog / DisGeNET / a hand-assembled table. **Out:** the same `disease_genes` slot.
**Operation:** *parsing and matching.* Requires a column of UniProt accessions; `association_score` may be absent (a curated list carries no Open Targets score). A symbol-only table is **rejected by default**, because `SYMBOL → UNIPROT` in `org.Hs.eg.db` is many-to-many and a silent one-per-symbol pick would put a nondeterministic mapping into the disease module. `map_symbols = TRUE` opts in to mapping a symbol column through `org.Hs.eg.db` (`clusterProfiler::bitr()`); every multi-mapping and every unmapped symbol is logged, all accessions of a multi-mapping symbol are kept, and `source` records that the mapping was applied.

---

## `network_*` — the network-pharmacology core

The substrate is one **bipartite compound–target graph per condition**:
compounds on one side, predicted targets on the other, an edge where a
compound is predicted to hit a target. Built by `network_build()`; every
other `network_*` function reads it.

### `network_build(proj)`
**In:** the binarized matrix + imported targets. **Out:** `network_edges.csv` — one row per (condition, compound, target) edge, with the prediction probability as `weight`.
**Operation:** *graph construction.* For each condition, keeps the compounds present in that condition and their target edges, de-duplicates repeated (compound, target) pairs (keeping the highest-probability one), and builds the graph. Every built graph is asserted bipartite (no compound_id may collide with a uniprot_id) before it is returned. A cached graph object is written for speed; the CSV is the truth — the cache carries a row-count + content checksum so a hand-edit to `network_edges.csv` is detected and the cache rebuilt, rather than silently served stale. `disease_association_score` (a per-edge join onto `targets_disease`) was removed: `targets_disease` accumulates rows across multiple diseases, and the old join's `match()` silently picked whichever disease happened to sort first, with no `disease_id` column recording which — a breaking schema change (see `NEWS.md`); join `targets_disease` yourself if you need it.

### `network_enrich(proj, condition, db = c("reactome", "go", "kegg"), universe = c("project", "genome"), pAdjustMethod = "BH")`
**In:** a condition's target set. **Out:** `network_enrichment.csv` — enriched GO terms / Reactome or KEGG pathways.
**Operation:** *over-representation test.* Maps the targets' UniProt IDs to Entrez IDs, then runs `clusterProfiler`'s **hypergeometric (Fisher) over-representation test** of that gene set against the chosen ontology, keeping terms below an adjusted-p cutoff (`p.adjust` via `pAdjustMethod`, default `"BH"`, newly exposed — previously hard-coded inside `clusterProfiler`). For GO, `simplify_go = TRUE` additionally collapses near-duplicate terms by **semantic similarity** (GOSemSim) at cutoff 0.7. The `qvalue` column is `NA` (not a crash) when `clusterProfiler` omits it, which happens routinely when the `qvalue` package's pi0 estimation fails on the small p-value vectors a 10–50-gene target set produces.
**Theory:** standard gene-set enrichment — "is this pathway hit more than you would expect if the target list were random?" (Boyle et al. 2004; Yu et al. 2012). "Than you would expect" is only meaningful once the background is fixed, and the foreground here is never a random gene sample: it is the reachable output of a ligand-similarity target predictor, structurally skewed toward GPCRs, kinases, nuclear receptors, proteases, and transporters — the only protein families with enough known ligands to build a similarity model from. Tested against clusterProfiler's whole-genome default (~20,000 `org.Hs.eg.db` genes), those families come out "enriched" for any input compound before a single biological difference between compounds is considered (Timmons, Szkop & Gallagher 2015, *Genome Biol* 16:186). **`universe = "project"` (default since this fix)** restricts the background to the project's own predictable proteome (`targets_imported$uniprot_id`, mapped to Entrez once per call) — the same background pool `network_degeneracy(universe = "project")` uses for its permutation null. `universe = "genome"` is kept as an explicit opt-out to the old, unrestricted behaviour. Changing the default changes every p-value the function reports; see `NEWS.md`.

### `network_centrality(proj, condition, measures = c("degree", "betweenness", "hub_score"), normalize = TRUE)`
**Out:** per-node centrality values.
**Operation:** *graph metrics, unweighted.* Degree = number of edges. Betweenness = fraction of shortest paths through the node (`igraph`). `hub_score` = `igraph::eigen_centrality(weights = NA)` — the principal (Perron) eigenvector of the adjacency matrix, scaled to `max = 1` (on this undirected graph it equals eigenvector centrality; the directional Kleinberg meaning does not apply). This replaces `igraph::hits_scores()$hub`, whose top eigenspace is degenerate on a two-mode graph so ARPACK returned an RNG-seeded arbitrary vector — the Perron eigenvector is unique per connected component and so deterministic. Edge weights are **ignored on purpose** — the weight is a prediction confidence, not a distance.
**Bipartite normalisation (`normalize = TRUE`, default):** the raw `degree` / `betweenness` columns are not comparable between compounds and targets (a compound's degree tops out at `|T|`, a target's at `|C|`, and `|T|` is usually far larger). `degree_norm` divides each node's degree by the *opposite* mode's size; `betweenness_norm` divides betweenness by the mode-specific maximum betweenness of a bipartite graph of these mode sizes (**Borgatti & Everett 1997**, *Social Networks* 19:243; the `B_max` formulas cross-checked against NetworkX's `bipartite.betweenness_centrality`). `igraph::betweenness()` on an undirected graph counts each unordered pair once, so no factor-of-two correction is applied. Both `_norm` columns land in `[0, 1]`. `hub_score_component` re-runs `eigen_centrality()` per connected component (each scaled to `max = 1`) so nodes outside the largest component are not all ≈ 0 — it is a *within-component relative* ranking; `component_id` labels the component, and an isolated node gets `NA`. `n_compounds` / `n_targets` are recorded per condition. `normalize = FALSE` still emits every column (the schema is `normalize`-independent) but leaves the `_norm` / `component_id` / count columns `NA`.
**Theory:** classic centrality (Freeman 1978); two-mode normalisation (Borgatti & Everett 1997). The `_norm` values are comparable between compounds and targets **of the same condition**, but — like the R-index of `network_module_robustness()` — **not across conditions of different size**, because `|C|` and `|T|` are per-condition constants. On a disconnected graph `betweenness_norm` is a mild under-estimate (cross-component pairs contribute no paths; `B_max` assumes connectivity).

### `network_hub_penalty(proj, condition)`
**Out:** per target, a promiscuity-adjusted score.
**Operation:** *a closed-form re-weighting, no fitting.* With `p` = (compounds hitting this target) / (compounds in the condition), the score is `N · p · log(1/p)`. It is `0` when a target is hit by every compound and **peaks at `p ≈ 0.37`** — it up-weights *intermediate*-specificity targets, and is **not monotone** in selectivity despite the name.
**Theory:** the shape is `N` times a Shannon surprisal term (TF-IDF-like in spirit, though not TF-IDF). Flagged for revision — see `ROADMAP.md`.

### `network_module_robustness(proj, condition, clustering, attack, seed)`
**Out:** the graph's modules, an **R-index** of how robust each module is to node removal, and a per-node module-membership table.
**Operation:** *cluster, then simulate node removal.*
*(1) Cluster*, per connected component, by one of three backends:
- `clustering = "leiden"` (**default since 0.2.0**) — `igraph::cluster_leiden()` maximising **Newman modularity** (`objective_function = "modularity"`, *not* igraph's `"CPM"` default, which returns all singletons on a sparse graph). Valid on any graph; makes no bipartiteness assumption and guarantees well-connected communities (Traag et al. 2019). The reported `modularity` is Newman's `Q`.
- `clustering = "bipartite"` — `bipartite::computeModules(method = "Beckett")` (Beckett 2016) maximising **Barber's bipartite modularity `Q_B`** (Barber 2007). *Only* valid on the genuine two-mode compound–target graph: Newman's null model allows compound–compound and target–target edges, which cannot exist here, so it over-reports modularity; Barber's null forbids them. Aborts on a non-bipartite graph; a star / single-mode connected component (`1×k` biadjacency matrix) falls back to a single module. The reported `modularity` is `Q_B` (node-weighted mean over components).
- `clustering = "hdbscan"` — `dbscan::hdbscan()` over the shortest-path-distance matrix. With `minPts = 2` every non-isolated node has core distance 1, so the mutual-reachability hierarchy **degenerates to the single-linkage hierarchy of an integer metric with 3–4 distinct values** — it estimates a near-constant "density" on a bipartite graph. Kept as a no-edge-density-assumption contrast, not the default; `modularity` is `NA`. Selecting it emits a one-time note.

Clustering is run **unweighted** — for Leiden the `weight` edge attribute is *deleted* before clustering (not passed as `weights = NA`, which leaves `strength()` reading it for the null model and silently returns all-singletons on any `NA` probability). The import probabilities are also not calibrated across prediction platforms, and every other topological `network_*` call is unweighted.

*(2) Attack.* `attack = "targeted"` (default) removes the current highest-degree node, recomputes degree, repeats — the "malicious attack". `attack = "random"` averages `n_random` uniformly-random removal orders (no recomputation) — "random failure"; the curve is stored as the replicate mean ± SD. `attack = "both"` runs both, giving the canonical two-curve percolation figure (Albert, Jeong & Barabási 2000).

*(3) R-index* = the **mean largest-component fraction over removal steps `Q = 1..N`** (Schneider et al. 2011). Higher = more robust; theoretical maximum `(N−1)/(2N) < 0.5`, **not** 1; **not comparable across modules of different size**. `r_index` holds the targeted value, `r_index_random` the random baseline; `r_index_random − r_index` is the interpretable quantity (a module far below its random baseline is hub-dependent). The seed independently covers the stochastic clustering step and the percolation tie-break / random orders.
**Theory:** targeted-attack percolation from network-robustness theory. `clustering` choice: Leiden always; bipartite only on the 2-mode compound–target graph (Barber's null vs Newman's); HDBSCAN degenerates to single linkage on the integer metric.
**References:** Schneider et al. 2011 *PNAS* 108(10):3838-3841; Traag et al. 2019 *Sci Rep* 9:5233; Barber 2007 *Phys Rev E* 76:066102; Beckett 2016 *R Soc Open Sci* 3:140536; Albert, Jeong & Barabási 2000 *Nature* 406:378-382.

### `network_motifs(proj, condition)`
**Out:** every feed-forward triple (compound → target → pathway/disease, with the compound also linked directly).
**Operation:** *enumeration by a table join, no statistics.* Lists all two-step paths in the directed layered graph and keeps those whose closing edge exists. A directed 3-cycle (feedback loop) is listed once per starting node — 3 rows per real cycle, kept deliberately (per-rotation visibility) — but the `projectLog()` summary line reports the deduplicated *cycle* count, not the row count.
**Theory:** the *shapes* are Milo et al. (2002)'s 3-node motifs, but **no null model is run**, and in this graph the closing edge is added by construction — so read the output as "which triples realise this pattern", not as statistically enriched motifs. Feedback loops are structurally impossible here (the graph has no cycles) and always come out zero.

### `network_degeneracy(proj, condition, annotation = "direct", ont = "BP", measure = "Wang", combine = "BMA", drop = "IEA", universe = "project", n_random = 200, seed)`
**Out:** per compound pair: `target_jaccard`, a `functional_similarity` (the renamed `pathway_jaccard`), its permutation `z_score` / right-tail `p_adjusted`, and `degeneracy_score = functional_similarity · (1 − target_jaccard)`.
**Operation.** `annotation = "direct"` (default) scores functional convergence as GO **semantic similarity** (GOSemSim; `GOSemSim::clusterSim()`-equivalent, best-match average, `round`ed to 0.001) between the two compounds' Entrez-mapped targets on each protein's *own* GO annotations — **no `network_enrich()` needed**, which removes the enrichment-universe circularity (a gene's GO annotations do not move when the condition's compound list changes). `measure = "Wang"` (Wang et al. 2007; `computeIC = FALSE`, shared `godata` cache with `network_enrich(simplify_go)`) for consistency and cost, **not** because it is bias-free; `drop = "IEA"` excludes electronic annotations (~44% of `org.Hs.eg.db` BP annotations, and the entire BP annotation of ~1,966 genes; Cheng et al. 2019 do the same). The observed similarity is z-scored against `n_random` resamplings in which each gene set is replaced by an equal-size distinct set drawn from the same **GO-annotation-count bins** (the `.network_value_bins()` consecutive-value binning `network_proximity()` uses for degree) over `universe` (`"project"` = every imported target, default). A fast permutation path — one pooled `GOSemSim::termSim` matrix, then matrix subsetting per permutation — makes ~10² permutations per pair cost seconds, not the ~60 h a per-permutation `clusterSim()` call would. `z_score` is the primary statistic; `p_empirical` / `p_adjusted` are the **right-tail** permutation p (high similarity is the interesting direction — the tail is flipped versus `network_proximity()`, whose identically-named `p_adjusted` is *left*-tail), BH-adjusted across the call, and resolution-limited (warned when `n_pairs / (n_random + 1) > 0.05`). `annotation = "enriched"` scores `GOSemSim::mgoSim()` over the enriched-term sets (keeps the circularity, no null); `annotation = "jaccard"` is the historical enriched-pathway Jaccard (no GOSemSim, no null). `n_pathways_a/b` are `NA` in `"direct"` mode; the row *set* is mode-dependent (a pair is skipped when it has nothing to compare).
**Theory:** "degeneracy" in the sense of Edelman & Gally (2001) — structurally different elements, same function. The multiplicative score is a convenience; in sparse data `target_jaccard ≈ 0` so it collapses toward `functional_similarity`. The annotation-count-matched permutation null is what corrects the curation bias. Read `target_jaccard` and `functional_similarity` together. Wang et al. 2007 *Bioinformatics* 23:1274; Schlicker et al. 2006 *BMC Bioinformatics* 7:302 (BMA).

### `network_proximity(proj, condition, disease, disease_genes = "disease_genes", n_random)`
**Out:** per compound: how topologically close its targets are to the disease's genes on the human interactome, as a z-score and an empirical p-value.
**Operation:** *a distance + a permutation test.* The disease module `T` is the disease's gene set from the `disease_genes` slot (`disease_genes_fetch()` / `disease_genes_import()`) — built independently of the compounds, so the statistic measures topology, not set membership. On the STRING protein–protein network (restricted to its largest connected component), the "closest" distance `d` = average, over the compound's STRING-mapped targets `S`, of the shortest-path hop count to the *nearest* gene in `T` (Guney et al. 2016, Eq. 2). Then `d` is recomputed `n_random` times with `S` and `T` each **replaced by random genes of matched degree** (bins built to hold ≥ 100 genes each); the disease side of the null is drawn once per `(condition, disease)` and reused across compounds. `z = (d_observed − mean(d_random)) / sd(d_random)`; `p_empirical` is the left-tail permutation p; `p_adjusted` is Benjamini–Hochberg across the compounds. A strongly negative z = closer than chance. Two columns record provenance: `disease_gene_source` (`"disease_genes"` or the legacy `"targets_disease"`) and `n_overlap` = |S ∩ T| on the STRING-mapped sets, which equals `n_targets_mapped` exactly when `d_observed == 0` (i.e. `S ⊆ T`).
**The circularity that was fixed.** Before this was rewired, `T` came from `targets_disease_filter()`, which only annotates UniProt IDs already in `targets_imported` — so `T` was a subset of the compounds' own predicted targets, `S ∩ T` was large by construction, `d_observed` was dragged to 0 by exact zeros, and `z` went strongly negative for reasons unrelated to topology. `disease_genes = "targets_disease"` still runs that path for backward comparison, but emits a warning naming the circularity and stamps `disease_gene_source = "targets_disease"` on the rows so the CSV records it.
**Theory:** network-medicine proximity (Guney/Menche/Barabási) — a drug is more likely relevant to a disease if its targets sit near the disease module in the interactome.

### `network_synergy(proj, condition, disease, separation = "network", alpha = 0.05)`
**Out:** per compound pair: the Menche (2015) topological separation `s_AB`, each compound's disease-proximity predicate, and the Cheng et al. (2019) class `P1`–`P6`.
**Operation:** *one interactome distance matrix per condition, then a truth table.* On the STRING largest connected component (the same graph `network_proximity()` uses — `.network_string_lcc()`, keyed by `species`/`version`/`score_threshold`, and a mismatch between the two calls aborts), a single `igraph::distances()` over the union of every candidate compound's STRING-mapped targets gives every pair's

  `s_AB = ⟨d_AB⟩ − (⟨d_AA⟩ + ⟨d_BB⟩)/2`

where `⟨d_AA⟩ = |A|⁻¹ Σ_{a∈A} min_{a'≠a} d(a,a')` (closest convention) and `⟨d_AB⟩ = (|A|+|B|)⁻¹ (Σ_a min_b d(a,b) + Σ_b min_a d(b,a))` is the pooled symmetric form (Menche et al. 2015, Eq. 2–3). `s_AB ≥ 0` ⇒ `separated`. No null model is fitted for `s_AB`: Cheng et al. explicitly reject a z-score for the drug–drug relationship because each drug has on average ~3 targets and "the randomization procedure is not producing a Gaussian distribution" — `s_AB` is used raw.

Each compound is `proximal` when `z_score < 0` **and** `p_adjusted < alpha` from `network_proximity()`. `alpha` is **patliR's tightening — Cheng et al. use the z-score alone and state no p-value cut.** The classification is the Cheng Fig. 2 truth table over `(separated, proximal_a, proximal_b)`: `P1` overlapping/both-proximal, **`P2` Complementary Exposure = separated & both-proximal (the only class Cheng associate with efficacy)**, `P3`/`P4` one proximal, `P5`/`P6` neither. `complementary_exposure = (cheng_class == "P2")`.

**Singleton convention.** When a compound's STRING-mapped target set has < 2 members, `⟨d_AA⟩ = 0` (patliR's choice — Menche's `separation.py` returns `nan`) and the row is flagged `singleton_a`/`singleton_b`. This inflates `s_AB`, biasing single-target compounds toward `separated = TRUE` (→ `P2`/`P4`/`P6`); flagged pairs keep a visible `cheng_class` but are dropped from the `synergy_score` scalar.

**BH-family caveat.** `network_proximity()`'s empirical p has a floor of `1/(n_random+1)` and its Benjamini–Hochberg family is per-call; if `n_tests_in_family/(n_random+1) > alpha` the proximity gate is unreachable, every pair falls to `P5`/`P6`, and `network_synergy()` warns.

**Interactome scope.** patliR's STRING graph at `score_threshold = 400` includes text-mining and co-expression edges, whereas Cheng et al. (2019) built their interactome from physical/experimental interactions only. So the absolute `s_AB` values and the `P1`–`P6` boundary here are **not numerically comparable** to Cheng's published figures; the classification is meaningful for *relative* ranking within one patliR run, not as a cross-study effect size. Raise `score_threshold` (e.g. to 700) for a higher-confidence, sparser graph closer to Cheng's.

`separation = "jaccard"` skips STRING entirely: `s_AB` and the Cheng columns are `NA`, `synergy_score` reverts to its historical `both_proximal` gate. `target_jaccard`/`complementarity` are always on the **raw UniProt** sets in both modes; `n_targets_*_mapped` carries the STRING-mapped counts. `synergy_score` (kept as an interim ranking scalar) `= complementarity · (−max(z_a, z_b))`, gated on `complementary_exposure` in network mode.
**Theory:** Menche et al. (2015) module separation; Cheng, Kovács & Barabási (2019) *Nat Commun* 10:1197 Complementary Exposure.

### `network_bowtie(proj, condition, actions_score_threshold = 400)`
**Out:** each target labelled with its position in a bow-tie decomposition of STRING's *directed regulatory* network: `core`, `in_component`, `out_component`, `other` (tendrils), `not_in_action_network`, or `unmapped`.
**Operation:** *strongly-connected-component analysis.* Downloads STRING's directed "actions" channel (which no STRINGdb function exposes), filters it to `actions$score >= actions_score_threshold` (previously used every directed row regardless of STRING's own confidence, so the core SCC size was driven by the lowest-confidence predictions), finds the largest strongly connected component = the **core**, then the nodes that can reach the core (**in**) and that the core can reach (**out**). A core under 50 nodes or under 1% of the network's nodes triggers a `cli_warn` (every downstream number is noise at that size). UniProt → STRING_id mapping reads STRING's `protein.aliases` flat file directly rather than instantiating a second `STRINGdb` object at `score_threshold = 0` (which downloaded the full interactome just for an alias table that never depended on score). Computed once per `(species, actions_version, actions_score_threshold)` and cached; conditions just get annotated.
**Theory:** Broder et al. (2000)'s web-graph bow-tie, applied to a regulatory network — the core is the mutually-reachable regulatory hub, in/out are upstream/downstream.

### `network_pathview(proj, condition)`
**Out:** rendered KEGG pathway diagrams, one PNG per top enriched pathway, with each gene node coloured by a target score.
**Operation:** *rendering.* Calls `pathview::pathview()`, which downloads each pathway's KEGG map and colours it. The gene score is, per target, the highest prediction probability among the compounds hitting it (`"max_weight"`). `ok` requires the expected PNG to actually exist on disk, not just the absence of an R error (`pathview()` can warn and return without writing anything — no mappable nodes, or a KGML/PNG download that silently returns an HTML error page). Targets with an `NA` weight are explicitly dropped from the colour vector and logged (`aggregate()`'s implicit `na.action = na.omit` used to do this silently); if every target's weight is `NA` the function aborts with a clear message instead of handing `pathview()` an empty vector. A `cli_warn` fires if any score exceeds 1 before the `[0, 1]` colour clamp — the clamp assumes `targets_import()`'s `"96.5%"` → `0.965` convention, which a platform reporting raw 0–100 percentages would silently violate.
**Theory:** none — visualisation of the enrichment result.

### `network_filter_proteome(proj, proteome, proteome_label = NULL)`
**Out:** `network_edges` filtered to a user-supplied list of UniProt IDs, plus `proteome_label` (auto-derived from a hash of the sorted proteome if not given) and `n_proteome`.
**Operation:** *a filter.* Row-subset, nothing else. Deliberately does not feed back into the other `network_*` functions. `proteome_label` is part of the upsert key alongside `condition`, so two different proteomes filtered against the same condition are both kept, not overwritten. A `cli_warn` fires when `proteome` matches nothing in `network_edges` and fewer than half its entries look like UniProt accessions — the usual cause is passing gene symbols by mistake.

---

## `bias_*` — database-bias audit

### `bias_audit(proj, mad_threshold = 2.5)`
**In:** `reference_bioactivity` (from `refdb_build()`). **Out:** per compound and per target, a "promiscuous" flag.
**Operation:** *a robust outlier test.* For each compound, "degree" = number of *distinct* targets it has measured activity against in the reference database (and symmetrically for targets). Then `mad_score = (degree − median) / MAD`, where **MAD is the median absolute deviation** (a median-based, outlier-resistant analogue of the standard deviation; `stats::mad`). A compound/target with `mad_score > 2.5` — i.e. unusually *heavily studied* — is flagged `"promiscuo"`. Only the high side is flagged; a low degree just means under-studied. Nothing is removed.
**Theory:** Leys et al. (2013) — use median ± k·MAD, not mean ± k·SD, for outlier detection when the distribution is skewed (bioactivity counts are very skewed). The point is to surface where a network result might be an artefact of *curation effort* rather than biology.

### `bias_reweight(proj)`
**Out:** for each flagged compound/target, an adjusted score.
**Operation:** *a closed-form discount.* `score_adjusted = degree · log(N_universe / degree)` — the same TF-IDF-style log-ratio as `network_hub_penalty()`, but over the whole reference database. Only `"promiscuo"` entries are discounted; the rest keep their raw score.

### `bias_report(proj)`
**Operation:** aggregation only. Lists the flagged entries, ordered by adjusted score, with a fixed caveat that a flag does not invalidate a result.

---

## `rank_candidates()` / `plot_rank()` — combining criteria into one ranking

### `rank_candidates(proj, condition, disease = NULL, criteria = NULL, roll_up = "weighted_mean", top_n = 15, export = "none")`
**In:** whichever of `adme_filtered` (mandatory), `network_centrality` (mandatory), `network_hub_penalty`, `network_proximity`, `network_synergy`, `network_module_membership`/`network_module_robustness` (all optional, auto-detected) exist for the condition. **Out:** one row per compound with every criterion's raw value, its direction-normalised rank, an aggregate score, and a Pareto tier.
**Operation:** *Robust Rank Aggregation (RRA), not a weighted sum.* Each criterion is direction-normalised so higher is always better internally (`network_proximity`'s `z_score` is the one exception — flipped, since more negative is better) and each compound is ordered best-to-worst on that one axis (a compound with no value for a criterion sorts last — worst-rank imputation, `order(..., na.last = TRUE)`). `RobustRankAggreg::aggregateRanks(glist, N, method = "RRA", full = TRUE)` (Kolde et al. 2012) then combines these *rank* orderings — not the raw values — into one `rra_score` (lower = more robustly top-ranked across every criterion used) via a closed-form order-statistic (beta distribution) correction. Comparing ranks rather than magnitudes is deliberate: a criterion whose raw scale is broken or skewed (e.g. `betweenness_norm` sitting at 0 for most nodes on a sparse bipartite graph) cannot silently dominate the way it would in a weighted linear sum of z-scored criteria. `pareto_front` is a separate, from-scratch non-dominated sort over the same direction-normalised criteria (front 1 = no other compound is at-least-as-good on every criterion and strictly better on one) — kept alongside `rra_rank` as a second, non-collapsed view, never the primary sort key.
**Roll-ups.** `network_centrality` and `network_hub_penalty` live at the target level; each is rolled up to its compound via `roll_up` (default `"weighted_mean"`, weighted by the `network_build()` edge weight/import probability — so one low-confidence but highly central predicted target cannot hijack a compound's score the way `"max"` would). The centrality criterion itself is `rowMeans()` of whichever of `degree_norm`/`betweenness_norm`/`hub_score_component` are non-`NA` per target (aborts if `network_centrality` was run with `normalize = FALSE`, since there is then nothing to build a composite from). `network_synergy` is pairwise; each compound takes `max(synergy_score)` across every partner it was scored against (its best achievable pairing, not an average) — `n_p2_partners` (count of `cheng_class == "P2"` partnerships) is kept as an informational column only, not fed to the RRA, since it is a different view of the same table and would double-count the signal. `network_module_robustness`'s `r_index` is a direct join via `network_module_membership` (one compound, one module, no aggregation).
**Missing data, two different policies.** A whole criterion's table absent from the project (e.g. `network_synergy()` never run) drops that criterion entirely for the affected condition — logged, not imputed. A table present but missing one compound's value gets the worst-rank treatment above, tracked per row via `n_criteria_used`.
**Deliberately excluded criteria (2026-09-13, user's call).** Structural toxicity alerts (`tox_local`) and `bias_reweight()`'s `score_ajustado` were both considered and dropped from this version — not even offered as an opt-in override.
**Export.** `export = "sdf"`/`"smi"` writes the top-`top_n` (by `rra_rank`) compounds' structures via `rcdk::write.molecules()`/`rcdk::get.smiles()` and the same `.parse_smiles_safe()` identity chain `compounds_similarity()` uses — no new dependency, no re-lookup.
**Theory:** Kolde, Laur, Adler & Vilo (2012), *Bioinformatics* 28(4):573-580 — Robust Rank Aggregation, originally for combining ranked gene lists across studies/platforms of differing, non-comparable scales.

### `plot_rank(proj, condition, view = c("pareto", "heatmap"), ...)`
**Operation, `view = "pareto"`:** a new plain-`ggplot2` scatter (no dependency) over `rank_candidates()`'s output, coloured by the already-computed `pareto_front`; `x`/`y` default to `rra_rank` and whichever optional criterion is present, in priority order `crit_synergy_best` > `crit_proximity_z` > `crit_module_r_index` > `crit_hub_penalty` (falling back to `crit_centrality`) — an illustrative 2D placement of a decision that used every criterion, not the decision itself, and freely overridable.
**Operation, `view = "heatmap"`:** restricts to the top `top_n_compounds` (by `rra_rank`) and, among their edges only, the `top_n_targets` most relevant targets (`target_relevance = "breadth"` — how many top compounds share it — default, `"weight"`, or `"centrality"`), then renders with the same `pheatmap` machinery `plot_heatmap(what = "compound_target")` uses (factored into a shared `.plot_pheatmap_render()` helper) — no new geometry, answers "which targets do the best compounds actually converge on" where the full every-compound-x-every-target matrix is too dense to read.
**Design rationale:** `chilcuague-analysis/reviews/rank-candidates-design-spec.md`.

### `report_generate(proj, condition, out_dir = NULL, top_n = 15)`
**Out:** one self-contained HTML file per condition (`reports/report_<condition>.html`) covering compounds present, ADME filtering, the database-bias audit, modular robustness, `rank_candidates()`'s ranking, and the full `projectLog()`.
**Operation:** *base-R string templating, not `rmarkdown`.* Each section checks its own source table in `patliRResults(proj)` and either renders it (as a plain HTML `<table>`, cell values HTML-escaped) or shows a one-line "has not been run" placeholder — nothing aborts except the one true prerequisite, `network_build()` (needed to know which compounds are "present" in a condition at all). Deliberately **not** an `Rmd`/`pandoc` report: this development environment has no `pandoc` installed (`rmarkdown::pandoc_available()` returns `FALSE`), so an `Rmd`-based report could not be rendered or tested here even though `rmarkdown`/`knitr` are already `Suggests` for the vignette; and a handful of HTML tables is squarely inside the "reimplement it, don't depend" half of the package's own `DESIGN.md` rule for anything that is layout rather than a non-trivial algorithm.
**Theory:** none — presentation of already-computed results.

---

## `plot_*` — figures

All 16 take `proj` and return a plot object; most also write a PNG. They
**compute nothing new** beyond layout — they read a `*_results` table and
draw it. A few do a small transform for display:

| function | display transform |
|---|---|
| `plot_admet_radar` | polar coordinates by hand (`x = r·sinθ`, `y = r·cosθ`); one panel per compound |
| `plot_boiled_egg` | scatter of (TPSA, WLogP) over the two published ellipses |
| `plot_chemical_space` | **PCA** (`prcomp`, scaled) or UMAP of the `adme_local` descriptors to 2 or 3 axes; convex-hull "halos" per chemical family |
| `plot_network_layers` | force-directed layout (`igraph::layout_with_fr`) of the layered graph; only hub nodes labelled |
| `plot_network_degeneracy` | the same layout + one arc per degenerate compound pair |
| `plot_heatmap` | `pheatmap` — hierarchical clustering on both axes |
| `plot_gochord` | `GOplot::GOChord` ribbon diagram, genes ↔ GO terms |
| `plot_bowtie` | alluvial flow, compound → bow-tie component |
| `plot_upset` / `plot_adme_upset` | hand-built UpSet (intersection bar + dot matrix) |
| `plot_venn` | 2-set proportional Venn (`ggVennDiagram`) |
| `plot_centrality` / `plot_robustness` / `plot_proximity` / `plot_synergy` | bar / curve / lollipop / scatter of the matching result table |
| `plot_structure2d` | grid of the pre-rendered 2D structure PNGs |

---

## `launch_app()`

Starts a Shiny wizard that walks the pipeline step by step. Orchestration
only — every button calls one of the functions above.
