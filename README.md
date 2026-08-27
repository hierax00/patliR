# patliR

![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)
![R >= 4.3.0](https://img.shields.io/badge/R-%3E%3D%204.3.0-blue.svg)
![Status: pre-1.0](https://img.shields.io/badge/status-pre--1.0%20%2F%20active%20development-orange.svg)

Reproducible network pharmacology analysis for natural product extracts and
compound mixtures, from a raw compound list or GC-MS abundance matrix all the
way to a fully-logged candidate ranking report.

Developed at the Laboratorio de Investigación Química y Farmacológica de
Productos Naturales, UAQ.

## Contents

- [Status](#status)
- [Design in one paragraph](#design-in-one-paragraph)
- [Install](#install-development-version)
- [Quick example](#quick-example)
- [Documentation](#documentation)
- [Bundled reference data](#bundled-reference-data)
- [License](#license)

## Status

Core pipeline implemented and tested, pre-1.0. Covers compound import and
validation (`prep_*`), reference-database sync (`refdb_*`), natural-product
family classification and molecular similarity, local + imported ADME
(`adme_*`) and toxicity (`tox_*`) screening, partial target
import/disease-association filtering (`targets_*`), the full
network-pharmacology core (`network_*`, 12 functions: build, enrich,
KEGG pathview, centrality/hub penalty, module robustness, layers/motifs,
degeneracy, proximity, synergy, bow-tie), 16 `plot_*` visualization
functions, a partial domain-bias audit (`bias_*`), an optional Shiny
wizard (`launch_app()`), and a couple of standalone utilities
(`network_filter_proteome()`, `patliR_export_llm()`).

Not yet implemented: `rank_*` (candidate prioritization), `dock_*`/
`report_*` (docking prep and final report generation), `coconut_*`/
`tcm_*` (COCONUT/TCM database import), `targets_bipartite()`/
`targets_consensus()`, KEGG-directed network completion, and an AI
narration module — see [`ROADMAP.md`](ROADMAP.md).

For what every function does — signature, output, design rationale — read
its help page. For the cross-cutting design decisions, see
[`DESIGN.md`](DESIGN.md).

## Design in one paragraph

Every pipeline step is a pure function: `proj <- some_step(proj, ...)`.
`proj` is an S4 `PatliRProject` object, but the durable source of truth is
always plain CSV files written to the project directory — the analysis can
resume in a brand new R session, on a different machine, without any
R-specific binary format, via `patliR_load()`. Every step also appends to a
run log (`projectLog(proj)`), so every filtering decision, random seed, and
data source is traceable after the fact.

## Install (development version)

```r
# install.packages("devtools")
devtools::install_local(".")
# or, while developing:
devtools::load_all(".")
```

Java (>= 8) is required for `rcdk`/`rJava` cheminformatics routines. Optional
`Suggests` packages unlock specific functions (Bioconductor packages for
`network_enrich()`/`tox_safetyome()`/`network_pathview()`, `ggalluvial`/
`ggVennDiagram`/`patchwork`/`umap` for parts of the `plot_*` family) — see
each function's own documentation (`?network_enrich`, `?network_pathview`,
`?plot_chemical_space`, ...) for what it needs.

## Quick example

```r
library(patliR)

proj <- patliR_project("my_analysis")
compound_list <- read.csv(system.file("extdata", "input_compound_list.csv", package = "patliR"))
proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
proj <- adme_local(proj)
proj <- tox_local(proj)

compounds(proj)
projectLog(proj)
```

## Documentation

- Each function has its own help page (`?prep_compounds`, `?network_build`, ...).
- [`DESIGN.md`](DESIGN.md) — the cross-cutting design decisions.
- [`ROADMAP.md`](ROADMAP.md) — what is designed but not yet built.

## Bundled reference data

A few functions ship small, curated reference tables under
`inst/extdata/` so they work offline, out of the box:

- **PAINS** (480 filters) — Baell & Holloway (2010), *J. Med. Chem.*
  53(7), 2719-2740. Transcribed verbatim from RDKit's
  `Data/Pains/wehi_pains.csv` (BSD-3-Clause).
- **Brenk** (105 alerts) — Brenk et al. (2008), *ChemMedChem* 3, 435-444.
  SMARTS text from PatWalters/rd_filters' `alert_collection.csv` (MIT),
  cross-validated against RDKit's own compiled `FilterCatalogs.BRENK`.
- **Safetyome core panel** (500 genes) — Liu et al. (2026),
  *Toxicological Sciences* 209(3), kfag021. Transcribed from the paper's
  Supplementary Table 4. Redistribution terms for this specific table
  were not independently confirmed at the time of writing.

## License

MIT. See [`LICENSE`](LICENSE).
