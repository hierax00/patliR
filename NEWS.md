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
  core panel, import from external toxicity platforms, and reporting.
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

See [`patliR_manual.md`](patliR_manual.md) (section "Planeado / no
implementado") for what's designed but not yet built.
