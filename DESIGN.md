# patliR design notes

Why the package is shaped the way it is. For what each function does, read
its help page (`?prep_compounds`, `?network_build`, ...); this document
only covers the cross-cutting decisions that are not obvious from any one
function.

## The project object is immutable; every step is a pure function

The pipeline is a chain of `proj <- step(proj, ...)` calls. `proj` is an S4
`PatliRProject`, never mutated in place (R6 was rejected for this reason).
Each step is an explicit, inspectable input -> output transformation, which
makes a run reconstructable after the fact and makes tests independent of
each other.

## CSV is the source of truth; any `.rds` is a disposable cache

Every step writes its output as plain CSV inside the project directory
(`01_compounds.csv`, `02_matrix_raw.csv`, `03_binarized.csv`,
`results/<name>.csv`, ...) in addition to returning the in-memory object.
`patliR_load(project_dir)` rebuilds `proj` from those CSVs, so an analysis
resumes in a fresh R session, or on another machine, with no dependency on
any R-specific binary format.

- The reference database (`results/reference_compounds.csv`,
  `results/reference_bioactivity.csv`) is stored long/relational, not with
  list-columns, so it stays 100% CSV.
- `network_build()` follows the same rule: the edge list
  (`results/network_edges.csv`) is the real output; the per-condition
  `igraph` object is cached separately under `cacheDir(proj)` purely for
  speed. Deleting the cache is always safe.

## External resources go through one shared failure policy

Every function that calls an external API (PubChem, ChEMBL, Open Targets,
KEGG, STRING, ...) goes through `.fetch_external()`, which has two modes:

- `abort` — fail loudly, no fallback, when continuing without the data
  makes no sense.
- `warn_and_cache` (default for enrichment steps) — use the local cache if
  present, otherwise warn and return `NULL`. Never breaks the pipeline.

Error text from a failed external call is passed through `.cli_escape()`
before being shown, because a raw JSON fragment echoed back by an API can
contain `{`/`}` that `cli` would otherwise try to evaluate as glue.

## No web scraping — the user exports, patliR imports

Platforms without an official API (SwissADME, ADMETlab, SuperPred, ...) are
never scraped. The user runs the platform, downloads its CSV export, and
`adme_import()` / `tox_import()` / `targets_import()` reconcile it against
`compounds(proj)` with guided column mapping. `adme_export_smiles()` closes
the loop for platforms that only take a pasted-in SMILES list: it writes
the list plus a bridge CSV (`row_order`, `compound_id`, ...) so the
platform's own export can be matched back by row order.

## Structure identity: CDK canonical SMILES, not InChIKey

Each compound carries `canonical_smiles` from
`rcdk::get.smiles(mol, smiles.flavors(c("Canonical", "UseAromaticSymbols")))`.
`rcdk` on CRAN has no InChIKey function. Consequence: a CDK canonical SMILES
is not guaranteed to match the canonical SMILES another toolkit (RDKit,
PubChem, ChEMBL) produces for the same molecule. Functions that reconcile
against an external source re-canonicalize with the same CDK flavor, or
prefer an exact ID (PubChem CID, ChEMBL structure search) when one is
available, rather than comparing raw SMILES.

## Two input scenarios

- **A** — a GC-MS abundance matrix with replicate columns named
  `R<n>-<CONDITION>` -> `prep_binarize()` (averages replicates per
  condition, then binarizes; `matrix_raw` is kept so the continuous value
  stays available).
- **B** — a curated table of name + PubChem CID + SMILES ->
  `prep_compounds()` directly.

`prep_binarize()` forks the rest of the pipeline per condition column; the
`network_*` family carries that same per-condition split through.

## `tox_*` and `bias_*` never produce a pass/fail verdict

A PAINS/Brenk structural alert, or a "promiscuous" bias flag, very often
reflects real pharmacological activity or an under-studied compound rather
than a genuine liability. `tox_local()` has no hard cutoff (unlike
`adme_filter()`); `tox_report()` and `bias_report()` always carry a fixed
disclaimer. Excluding a compound is left to the user's explicit judgment.

## `plot_*`: static first, minimal Suggests

All 16 `plot_*` functions default to a static, publication-resolution
image and add a `Suggests` package only when the geometry is unreasonable
to build by hand (`ggalluvial`, `ggVennDiagram`, `pheatmap`, `GOplot`).
`ggrepel`, `ggraph`, and `UpSetR`/`ComplexUpset` were deliberately
reimplemented in plain `ggplot2` instead of taken as dependencies. Each
network analysis gets its own view rather than one generic `plot_network()`
— compound-target-pathway graphs saturate visually too fast for that.

## Analysis `Suggests`: reimplement layout, depend for algorithms

The "reimplement it" rule above is about *layout* geometry. It does not
extend to non-trivial *algorithms*. `network_module_robustness(clustering
= "bipartite")` depends on `bipartite` (a `Suggests`) for Barber's `Q_B`
maximisation — an NP-hard optimisation with a decade of algorithmic
literature (Barber 2007; Beckett 2016), not a geometry anyone should
hand-roll. It is guarded by the standard `requireNamespace()` +
`cli_abort` pattern and is needed only for the non-default `clustering`
value. One real cost to note: `bipartite`'s `Depends:` attaches `sna` and
`vegan` to the search path, and `sna` masks `igraph`'s `degree()` /
`betweenness()` / `closeness()` for the rest of the session — an argument
for keeping `"bipartite"` a deliberate opt-in, never the default.
`dbscan` (for `clustering = "hdbscan"`) is a `Suggests` on the same
footing.

## Explicitly out of scope

- Parsing raw GC-MS instrument output (`prep_gcms()`) — that is
  metabolomics; patliR starts from a compound x condition matrix the user
  has already assembled.
- Molecular dynamics.
- A docking engine or 3D viewer of its own. `dock_prepare()` /
  `dock_parse()` (planned, see `ROADMAP.md`) stay engine-agnostic.
- Binding-pocket prediction (P2RANK/fpocket) — `box_center` is the user's
  to provide.
