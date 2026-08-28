# patliR (development version)

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
