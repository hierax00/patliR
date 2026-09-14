# patliR

![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![R >= 4.3.0](https://img.shields.io/badge/R-%3E%3D%204.3.0-blue.svg)
![Status: pre-1.0](https://img.shields.io/badge/status-pre--1.0-orange.svg)

**Reproducible network-pharmacology analysis for natural-product extracts.**
From a plain compound list (or a GC-MS abundance matrix) all the way to a
fully-logged compound–target–pathway network — every step a pure function,
every intermediate a plain CSV.

Developed at the Laboratorio de Investigación Química y Farmacológica de
Productos Naturales, UAQ.

---

## The pipeline

```mermaid
flowchart TD
    subgraph input [" "]
        L["compound list<br/><i>name · PubChem CID · SMILES</i>"]
        M["GC-MS abundance matrix<br/><i>replicates × conditions</i>"]
    end

    L --> PC["<b>prep_compounds</b><br/>validate structures"]
    M --> PB["<b>prep_binarize</b><br/>average replicates → presence/absence"]
    L -.list-only.-> PAC["<b>prep_as_condition</b><br/>list = one extract"]
    PC --> PAC

    subgraph chem ["identity & chemistry"]
        RDB["<b>refdb_build</b><br/>PubChem · ChEMBL"]
        CL["<b>compounds_classify</b><br/>NPClassifier family"]
        S2D["<b>prep_structure2d</b>"]
        SIM["<b>compounds_similarity</b>"]
    end
    PC --> RDB & CL & S2D & SIM

    subgraph pk ["physchem · ADME · toxicity"]
        AL["<b>adme_local</b><br/>Ro5 · Veber · Ghose · Egan · Oprea · BOILED-Egg"]
        AF["<b>adme_filter</b>"]
        AI["<b>adme_import</b> ← SwissADME / ADMETlab"]
        TL["<b>tox_local</b><br/>PAINS · Brenk"]
        TS["<b>tox_safetyome</b><br/>500-gene safety panel"]
        TR["<b>tox_report</b>"]
    end
    PC --> AL --> AF
    AI --> AF
    PC --> TL --> TR
    TS --> TR

    TR --> CSV[["triage CSV<br/><i>pick the shortlist here</i>"]]
    AF --> CSV

    CSV --> TI["<b>targets_import</b> ← SuperPred / SwissTarget"]
    PB --> NB
    PAC --> NB
    TI --> NB["<b>network_build</b><br/>compound–target graph / condition"]

    DG["<b>disease_genes_fetch</b><br/>independent disease module<br/>(Open Targets, min_score = 0.4)"]

    subgraph net ["network pharmacology"]
        NE["<b>network_enrich</b><br/>GO · Reactome · KEGG (universe = project)"]
        NC["<b>network_centrality</b> · <b>hub_penalty</b>"]
        NMR["<b>module_robustness</b><br/>R-index percolation"]
        NM["<b>network_motifs</b> · <b>degeneracy</b>"]
        NP["<b>network_proximity</b><br/>distance to disease module"]
        NS["<b>network_synergy</b><br/>Cheng P1-P6"]
        NBT["<b>network_bowtie</b>"]
    end
    NB --> NE --> NC --> NMR --> NM
    DG --> NP
    NB --> NP --> NS
    NB --> NBT

    NE & NC & NMR & NP --> PLOTS["<b>plot_*</b><br/>16 publication figures<br/>+ 3-axis chemical space"]
    TR --> BIAS["<b>bias_audit</b><br/>database-bias flag"]
    RDB --> BIAS

    NC & NMR & NP & NS & AF --> RC["<b>rank_candidates</b><br/>Robust Rank Aggregation"]
    RC --> PR["<b>plot_rank</b><br/>Pareto · heatmap"]
    RC --> RG["<b>report_generate</b><br/>one HTML per condition"]
```

Two entry points: a **curated compound list** (`prep_compounds`) or a **raw
abundance matrix** (`prep_binarize`). Everything downstream is shared. The
network layer keys off one graph per experimental *condition*; with a plain
list, `prep_as_condition()` makes the whole list a single condition.

## Install

```r
# install.packages("remotes")
remotes::install_github("hierax00/patliR")
```

Java (≥ 8) is required for the `rcdk`/`rJava` cheminformatics routines
(structure validation, PAINS/Brenk matching, descriptors). Everything else
is optional and pulled in only when you call a function that needs it —
Bioconductor packages for `network_enrich()` / `tox_safetyome()` /
`network_pathview()`, `STRINGdb` for `network_proximity()` / `network_bowtie()`,
`plotly` for the 3-D chemical space, and so on. Each help page lists its own
requirements.

## Quick start — list to triage table

```r
library(patliR)

proj <- patliR_project("chilcuague")
compounds_in <- read.csv("compounds.csv")          # name, CAS, PubChemCID, SMILES

proj <- prep_compounds(proj, compounds_in, identifier = "smiles")
proj <- refdb_build(proj, sources = c("pubchem", "chembl"))   # reference DB, on load
proj <- compounds_classify(proj)                    # natural-product family
proj <- adme_local(proj)                            # drug-likeness rules + BOILED-Egg
proj <- adme_filter(proj, rules = c("ro5", "veber", "ghose", "egan", "oprea"))
proj <- tox_local(proj, alert_sets = c("pains", "brenk"))

report <- tox_report(proj)                          # → results/tox_report.csv
report$summary                                      # one row per compound
```

## Quick start — shortlist to network

```r
keep <- c("C0001", "C0004", "C0009", "C0021")       # your picks from the triage

proj <- prep_as_condition(proj, condition = "Chilcuague", compound_ids = keep)
proj <- targets_import_batch(proj, "targets_superpred/", platform = "superpred")
proj <- network_build(proj)
proj <- network_enrich(proj, condition = "Chilcuague", db = "go")
proj <- network_centrality(proj, condition = "Chilcuague")
proj <- network_module_robustness(proj, condition = "Chilcuague", seed = 42)

plot_network_layers(proj, condition = "Chilcuague")
plot_chemical_space(proj, dims = 3, color_by = "family", engine = "plotly")

## close the loop: one ranked table + one report, combining everything above
proj <- rank_candidates(proj, condition = "Chilcuague", export = "sdf")
plot_rank(proj, condition = "Chilcuague", view = "pareto")
proj <- report_generate(proj, condition = "Chilcuague")
```

A worked end-to-end script against a real dataset lives in
[`chilcuague-analysis/`](https://github.com/hierax00/chilcuague-analysis).

## Function families

| family | what it does |
|---|---|
| `prep_*` | import & validate compounds, binarize an abundance matrix, 2-D depiction, list-as-condition |
| `refdb_*` | local reference DB of identity + bioactivity (PubChem, ChEMBL) |
| `compounds_classify()` / `compounds_similarity()` | NPClassifier family; pairwise fingerprint similarity |
| `adme_*` | local drug-/lead-likeness rules + BOILED-Egg; import from external platforms; rule filtering; SMILES export bridge |
| `tox_*` | PAINS/Brenk structural alerts; target-level safety panel; import; per-compound report — **never a pass/fail verdict** |
| `targets_*` | import predicted targets; disease-association filtering (Open Targets) |
| `disease_genes_*` | independent disease gene module for `network_proximity()` (Open Targets, or a curated import) |
| `network_*` | build, enrich, centrality/hub-penalty, module robustness, motifs, degeneracy, proximity, synergy, bow-tie, KEGG pathview, proteome filter |
| `plot_*` | 16 static/interactive figures for every result above |
| `bias_*` | MAD-based "promiscuous compound/target" flag against the reference DB |
| `rank_candidates()` / `plot_rank()` | Robust Rank Aggregation over ADME + network criteria into one ranked table, with Pareto/heatmap views |
| `report_generate()` | one self-contained HTML report per condition |
| `patliR_export_llm()` | flat-text dump of a whole project for an LLM to read |

Everything is documented on its own help page.
[`METHODS.md`](METHODS.md) explains, in plain language, what each function
actually computes and the theory behind it. For the cross-cutting design
decisions see [`DESIGN.md`](DESIGN.md); for what is designed but not yet
built see [`ROADMAP.md`](ROADMAP.md).

## Design in one paragraph

Every step is a pure transformation `proj <- step(proj, ...)`. `proj` is an
immutable S4 object, but the durable source of truth is always a plain CSV
in the project directory — `patliR_load()` rebuilds the whole project from
those CSVs in a fresh R session, on another machine, with no R-specific
binary format. Every step appends to a run log, so every filtering
decision, random seed, and data source is traceable after the fact.
External services (PubChem, ChEMBL, KEGG, STRING, …) all go through one
retry/cache wrapper and never break the pipeline when they fail.

## Bundled reference data

Small curated tables ship under `inst/extdata/` so the local steps work
offline:

- **PAINS** — 480 filters, Baell & Holloway (2010), *J. Med. Chem.* 53(7),
  2719–2740; verbatim from RDKit's `wehi_pains.csv` (BSD-3-Clause).
- **Brenk** — 105 alerts, Brenk et al. (2008), *ChemMedChem* 3, 435–444;
  from PatWalters/rd_filters (MIT), cross-checked against RDKit's
  `FilterCatalogs.BRENK`.
- **Safetyome core panel** — 500 genes, Liu et al. (2026), *Toxicological
  Sciences* 209(3), kfag021, Supplementary Table 4. Redistribution terms
  for that table were not independently confirmed at the time of writing.
- **BOILED-Egg** GIA/BBB ellipse boundaries — digitized from Daina & Zoete
  (2016) via PyBOILEDegg (GPL-3).

## License

MIT — see [`LICENSE`](LICENSE).
