# patliR 0.2.0

First public release. A pre-1.0 version: the core pipeline is implemented
and tested, and argument names may still change before 1.0.

## Compounds and chemistry

- `prep_compounds()` / `prep_compound()` — import and validate compounds by
  PubChem CID or SMILES, with stereo-aware deduplication and optional
  PubChem fetch of missing SMILES. `prep_binarize()` turns a replicate
  abundance matrix into presence/absence per condition (Q1 threshold, or
  `min_replicates`); `prep_as_condition()` makes a plain list one condition.
- `refdb_build()` / `refdb_update()` — local reference database of identity
  and bioactivity from PubChem and ChEMBL, resolving ChEMBL identity
  deterministically (InChIKey exact match, skeleton fallback, SMILES
  flexmatch).
- `compounds_classify()` (NPClassifier) and `compounds_similarity()`
  (fingerprint Tanimoto/Dice/cosine).
- `adme_local()` — Lipinski (with Lipinski's own HBA/HBD counts), Veber,
  Ghose, Egan, Oprea and BOILED-Egg; `adme_filter()` with a configurable
  Ro5 tolerance; import from and SMILES export to external ADME platforms.
- `tox_local()` — PAINS (480) and Brenk (105) structural alerts;
  `tox_safetyome()` against the Safetyome core panel; `tox_report()`.
  Alerts are flags for expert review, never a pass/fail verdict.

## Targets and diseases

- `targets_import()` / `targets_import_batch()` — predicted targets from
  SwissTargetPrediction, SuperPred and similar platforms.
- `targets_disease_filter()` / `targets_disease_profile()` — Open Targets
  disease associations per target.
- `disease_genes_fetch()` / `disease_genes_import()` — an independent
  disease gene module for network proximity.

## Networks

- `network_build()`, `network_enrich()` (GO/Reactome/KEGG against a
  project-level universe), `network_centrality()`, `network_hub_penalty()`,
  `network_module_robustness()`, `network_motifs()`, `network_degeneracy()`
  (GO semantic similarity with a permutation null), `network_proximity()`
  (Menche/Guney z-score), `network_synergy()` (`s_AB` and Cheng
  complementary-exposure classes), `network_bowtie()`, `network_pathview()`,
  `network_kegg_topology()`, `network_filter_proteome()` and
  `network_condition_compare()`.

## Bias, ranking and reporting

- `bias_audit()` / `bias_reweight()` / `bias_report()` — MAD-based
  detection and discounting of over-studied compounds and targets in the
  reference database.
- `rank_candidates()` — Robust Rank Aggregation over ADME and network
  criteria, with Pareto tiers and an SDF/SMILES export of the top
  candidates; `plot_rank()`.
- `report_generate()` — one self-contained HTML report per condition.
- `plot_*` — static and interactive figures for every result above.
- `patliR_load()` resumes a project from its CSVs; `patliR_export_llm()`
  writes a flat-text dump; `launch_app()` is an optional Shiny wizard.
