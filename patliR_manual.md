# patliR — Manual del Creador

Referencia única del proyecto: qué hace cada función implementada (por
familia, con firma input/output y el porqué de las decisiones de diseño
que no son obvias leyendo el código), y — en la sección final,
[Planeado / no implementado](#planeado--no-implementado) — todo lo
diseñado pero sin código todavía, incluidas las decisiones de alcance de
1.0 vs. 1.1. Antes eran dos documentos (`patliR_manual.md` + `ROADMAP.md`)
que se referenciaban mutuamente; fusionados en uno el 2026-08-23 para que
"qué hace patliR hoy" y "qué le falta" se lean como una sola cosa.

Este documento describe el estado actual del código y lo planeado para
adelante. Para **cómo** se llegó hasta acá — bugs reales encontrados
contra APIs/datos reales, decisiones tomadas sobre la marcha, fechas — hay
una bitácora de desarrollo cronológica que se mantiene fuera de este
repositorio público.

Última actualización: 2026-08-25.

---

## 0. Arquitectura general

### Objeto de proyecto: funcional, no R6

Se descartó R6 (mutable, `proj$fn()`) a favor de un objeto S4 inmutable
(`PatliRProject`), con funciones puras: `proj <- prep_compounds(proj, ...)`.

**Por qué**: coherencia con el principio rector de transparencia — cada
paso del pipeline es una transformación explícita e inspeccionable
input→output. Con mutación in-place es más difícil garantizar que un
reporte pueda reconstruir con exactitud qué pasó en cada paso. También es
más fácil de testear (sin estado compartido entre tests) y de pipear con
`|>`.

### Persistencia: CSV como fuente de verdad, RDS como caché derivada

- Cada paso del pipeline escribe su output a CSV dentro de la carpeta del
  proyecto (`01_compounds.csv`, `02_matrix_raw.csv`, `03_binarized.csv`,
  ...), además de devolver el objeto en memoria.
- `patliR_load(project_dir)` reconstruye `proj` leyendo esos CSVs de
  vuelta — el análisis se puede retomar sin depender de que la sesión de R
  siga viva.
- La base de referencia (`reference_compounds.csv`,
  `reference_bioactivity.csv`) vive en formato largo/relacional (no
  list-columns) para que sea 100% CSV y no dependa de tipos exclusivos de
  R.
- Cualquier `.rds` es una caché de performance, siempre regenerable desde
  los CSV. Nunca es requisito — el análisis debe poder continuar solo con
  los CSV. `network_*` sigue el mismo principio: guarda las aristas como
  CSV (`network_edges`) y cachea el `igraph` derivado por separado (ver
  familia `network_*` abajo).

**Por qué**: el análisis debe poder continuar en cualquier paso sin
depender de un formato binario específico de R — portabilidad y
reproducibilidad entre máquinas/sesiones.

### Manejo de red externa — patrón compartido

Helper interno `.fetch_external()`, reusado por toda función que pega a
una API/fuente externa, con dos modos de falla configurables:

- `abort`: falla explícita (`cli::cli_abort()`), sin fallback — para
  cuando no tiene sentido continuar sin el dato.
- `warn_and_cache`: usa caché local si existe; si no, `cli::cli_warn()` y
  sigue con `NULL` — nunca rompe el pipeline. Para enriquecimiento
  opcional.

El texto de `conditionMessage()` del error real (potencialmente texto
externo sin control, p. ej. un fragmento JSON crudo de una API que
falló) se pasa por `.cli_escape()` (fix 2026-08-16, ver `DEVLOG.md` --
bug real corriendo `compounds_classify()` contra Chilcuague: NPClassifier
devolvió un mensaje con `{`/`}` literales, y `cli::cli_warn()`/
`cli::cli_abort()` interpretan cualquier `{...}` de un bullet como una
expresión glue a evaluar, no texto literal -- eso tronaba con un error de
parseo que enmascaraba el error real) antes de insertarse en cualquier
bullet de `cli_warn()`/`cli_abort()`.

### Scraping de plataformas web: descartado

Se descartó automatizar plataformas sin API oficial (p. ej. SwissADME) vía
scraping — no es mantenible ni éticamente deseable. Patrón adoptado en su
lugar: el usuario corre la plataforma que prefiera en su navegador,
descarga el CSV que esa plataforma exporta, y `patliR` lo importa con
mapeo de columnas guiado (`adme_import()`, `tox_import()`,
`targets_import()`).

### Fuera de alcance (explícito)

- `prep_gcms()` — parsear salida cruda de instrumento GC-MS. Eso es
  metabolómica, un paquete/alcance aparte. `patliR` trabaja sobre matrices
  compuesto×condición ya armadas por el usuario.
- Dinámica molecular.
- Motor de docking propio, visualizador 3D — ver más abajo, 1.5, para lo
  que sí está en alcance (`dock_prepare()`/`dock_parse()`, 100% agnóstico
  al motor).
- Predicción de bolsillo de unión (P2RANK/fpocket) — `box_center` queda en
  manos del usuario.

---

## 1. Familia `prep_*` — entrada y normalización

```r
prep_compounds(
  data, identifier = c("pubchem", "smiles"),
  id_col = NULL, name_col = NULL, dedup = TRUE, cache_dir = proj$cache_dir
)
# Output: proj$compounds (id interno, pubchem_id, name, canonical_smiles, source)
#         proj$log (id, motivo: estructura_invalida | duplicado | sin_smiles_resuelto)
# Escribe 01_compounds.csv

prep_compound(data_row, identifier = c("pubchem","smiles"), ...)
# Wrapper delgado de un solo compuesto -> llama prep_compounds() con 1 fila.

prep_structure2d(proj, engine = c("chemminer","rcdk"), out_dir = NULL)
# QC visual/geométrico adicional (más allá de la validación de rcdk en prep_compounds).
# Output: proj$structure2d_log (id, generado_ok, engine_usado)

prep_binarize(
  matrix, id_col = "compound_id", replicate_sep = "_R", average_replicates = TRUE
)
# Input: compound_id + columnas "R<n>-<CONDICION>" (ej. R1-LEA-ET, R2-FLO-AQ
# -- réplica primero, la condición puede tener guiones adentro).
# Promedia réplicas por condición, luego Q1 por condición sobre el promedio.
# Output: proj$matrix_raw   (promedio por condición, SIN binarizar, se conserva)
#         proj$binarized    (0/1 por condición)
#         proj$log (condicion, ya_binaria, promedio_replicas, q1_usado)
# Escribe 02_matrix_raw.csv y 03_binarized.csv -- ninguna pisa a la otra.
```

**Identidad estructural: SMILES canónico de CDK, no InChIKey.** Cada
compuesto lleva `compounds$canonical_smiles`
(`rcdk::get.smiles(mol, smiles.flavors(c("Canonical", "UseAromaticSymbols")))`),
no un InChIKey — `rcdk::get.inchi.key()` no existe en el `rcdk` de CRAN
(ver `DEVLOG.md`). Costo aceptado: un SMILES canónico de CDK no
necesariamente coincide con el canónico de otro toolkit (RDKit, OpenBabel,
PubChem, ChEMBL) para la misma molécula — cada función que reconcilia
contra una fuente externa lo hace re-canonicalizando con el mismo flavor
de CDK, o prefiriendo un ID exacto (PubChem CID, endpoint `flexmatch` de
ChEMBL) cuando está disponible, en vez de comparar SMILES crudo.

**Por qué se conserva `matrix_raw`**: se necesita el valor continuo (no
solo presencia/ausencia) para poder colorear por concentración relativa
además de por confianza de predicción en visualizaciones futuras (ver
más abajo, 1.6) — nunca se descarta información de input que podría
reusarse.

Dos escenarios de input reconocidos: **Escenario A** (matriz GC-MS con
réplicas, formato `R<n>-<CONDICION>` de arriba) → `prep_binarize()`;
**Escenario B** (tabla curada nombre + PubChemCID + SMILES) →
`prep_compounds()` directo.

### Mini-familia `refdb_*` — infraestructura compartida

```r
refdb_build(proj, sources = c("pubchem","chembl"))
refdb_update(proj, compounds)
refdb_rebuild_cache(proj)   # regenera el .rds derivado desde los csv
```

Fuentes: PubChem + ChEMBL para identidad/bioactividad de compuesto
(`sources = "coconut"` da el mismo error explícito de "no implementado
todavía" que el resto de funciones con fuentes pendientes — ver
más abajo, 1.2). `refdb_build()`/PubChem prefiere `pubchem_id` (CID
exacto) sobre búsqueda por estructura; `refdb_build()`/ChEMBL usa el
endpoint `molecule_structures__canonical_smiles__flexmatch`, diseñado por
la propia API de ChEMBL para reconciliar SMILES entre toolkits distintos.

### `compounds_classify()` — familia química (NPClassifier)

```r
proj <- compounds_classify(proj, fetch_mode = "warn_and_cache")
# Output: proj$compounds_classified
#   (compound_id, pathway, superclass, class, isglycoside, source, fetch_date)
```

Clasifica cada compuesto por familia de producto natural vía
[NPClassifier](https://npclassifier.gnps2.org) (Kim, H.W. et al. 2021,
*J. Nat. Prod.* 84(11), 2795-2807, doi:10.1021/acs.jnatprod.1c00399), una
API REST gratuita y sin autenticación (`GET
https://npclassifier.gnps2.org/classify?smiles=...`) especializada en
productos naturales (se prefirió sobre ClassyFire por esa
especialización). Usa el mismo patrón `.fetch_external()`
(`abort`/`warn_and_cache`) y caché que `refdb_build()` — nunca vuelve a
pegarle a la API para un compuesto ya cacheado salvo que se le pida
explícitamente. `pathway` es la columna que usa
`plot_chemical_space(color_by = "family")` por default (ver sección 6).

**Nota real (2026-08-11)**: el endpoint de NPClassifier responde con
`text/html` (no JSON, status 200) bajo un patrón consistente con
rate-limit/anti-bot no documentado — `.npclassifier_lookup()` usa
`httr2::req_throttle()` (1 req/s) + `req_retry()` con backoff exponencial
que también reintenta cuando el content-type no es JSON, no solo en
429/503 (ver `DEVLOG.md`, 2026-08-11).

### `compounds_similarity()` — similitud molecular pareada

```r
compounds_similarity(
  proj, compound_ids = NULL,
  fingerprint_type = c("standard","extended","circular","maccs","pubchem"),
  method = c("tanimoto","dice","cosine")
)
# Output: proj$compounds_similarity -> tibble(compound_a, compound_b, similarity)
#   [formato largo pareado, mismo patrón que network_degeneracy()]
```

Implementada 2026-08-16. Huella molecular vía `rcdk::get.fingerprint()`,
similitud pareada vía `fingerprint::fp.sim.matrix()` — su firma real
(matriz N×N simétrica, diagonal en 1, sin dimnames) se confirmó contra el
código fuente real del paquete `fingerprint` (`R/matrix.R`), no solo su
documentación. `fingerprint` es dependencia dura de `rcdk` (ya en
`Imports`), agregado también a `Suggests` porque CRAN exige declarar
cualquier paquete al que se le hace `::` directo aunque llegue
transitivamente. Función nueva independiente — no integrada a
`plot_chemical_space()`/`network_synergy()` todavía; ambas siguen siendo
candidatas de integración futura sin bloquear esto (ver [Planeado / no
implementado](#planeado--no-implementado)).

---

## 2. Familia `adme_*` — ADME

Separada de `tox_*` (sección 3) porque toxicidad nunca debe funcionar como
filtro duro.

```r
adme_local(proj, compounds = NULL, routes = c("oral","topical","ophthalmic","injectable"))
# Ro5 + Veber + Ghose (completo, con AMR) + Egan + Oprea (lead-likeness) +
# BOILED-Egg (absorción GI / BBB) + reglas heurísticas por vía de
# administración. Todo vía rcdk, sin red.
# Output: proj$adme_local (compound_id, mw, logp, hbd, hba, tpsa, rotatable_bonds,
#                           amr, n_rings_approx, ro5_pass, veber_pass, ghose_pass,
#                           egan_pass, oprea_pass, gi_absorption, bbb_permeant,
#                           route_oral, route_topical, route_ophthalmic, route_injectable)

adme_import(proj, path, platform = c("swissadme","admetlab","other"), column_map = NULL)
# CSV exportado de SwissADME/ADMETlab/pkCSM/etc. column_map = NULL con
# platform="swissadme"/"admetlab" usa presets ya confirmados contra
# archivos de ejemplo reales:
#   swissadme: Molecule->pubchem_id, MW->mw, Consensus.Log.Po.w->logp, ...
#   admetlab:  smiles->smiles (join vía canonical_smiles, no pubchem_id --
#              no todas las plataformas exportan con el mismo identificador)
# Output: proj$adme_imported (compound_id, propiedad, valor, fuente, fecha_import)

adme_filter(
  proj, rules = c("ro5","veber","ghose","egan","oprea","route"), source = c("local","imported"),
  hard_cutoff = FALSE, ask = interactive()
)
# Default: solo marca pasa/no-pasa por regla (nunca elimina en silencio).
# "oprea" es la regla de lead-likeness (más angosta que drug-likeness,
# pensada para triage hit-to-lead, no para un candidato final).
# hard_cutoff = TRUE:
#   1. calcula cuántos se removerían y por qué regla, se lo muestra al usuario
#   2. ask=TRUE (interactivo) -> prompt estilo update.packages(): "¿Remover N de M? [y/n/v]"
#   3. ask=FALSE (script) -> aplica el corte igual, pero SIEMPRE deja detalle en el log
# Output: proj$adme_filtered (compound_id, regla, pasa, valor, umbral) [formato largo]

adme_export_smiles(proj, compound_ids = NULL, out_file = NULL)
# Exporta los compuestos del proyecto como lista plana de SMILES (uno por
# renglon, canonical_smiles si existe, si no smiles crudo) lista para
# pegar directo en SwissADME/plataformas similares que no aceptan CSV de
# entrada. out_file = NULL (default) solo devuelve el resultado en
# memoria; con out_file escribe <out_file>.txt (la lista) y
# <out_file>_map.csv (row_order, compound_id, name, smiles) -- row_order
# coincide con el orden de renglon en el .txt, que es el mismo orden en
# que la plataforma va a devolver su propio export (sin id propio).
# Cierra el ciclo: adme_export_smiles() -> pegar en la plataforma ->
# descargar su export -> adme_import().
# Output: invisible list(smiles_text, mapping)
```

Implementada 2026-08-25, pedido directo de Uriel en chat (no bug, no
diseño previo) al preparar Chilcuague completo para exportar a SwissADME.
Compuestos sin ningún SMILES se excluyen del export con aviso explícito
(nunca se manda una fila vacía a la plataforma en silencio). Archivo:
`R/adme_import.R`.

### GI absorption / BBB permeant: BOILED-Egg real

`gi_absorption`/`bbb_permeant` usan las elipses publicadas reales del
BOILED-Egg (Daina & Zoete 2016, *ChemMedChem* 11, 1117-1121), digitalizadas
desde [bfmilne/PyBOILEDegg](https://github.com/bfmilne/PyBOILEDegg) (GPL-3,
misma fuente), empaquetadas en `inst/extdata/boiled_egg_{gia,bbb}.csv`,
con test punto-en-polígono real. Única aproximación que queda: el eje
WLogP del modelo original (Wildman-Crippen, RDKit) no tiene equivalente
exacto en `rcdk`/CDK — se usa `rcdk::get.alogp()` (Ghose-Crippen, método
relacionado pero no idéntico) como proxy, documentado en `?adme_local`.

### Referencias — criterios físico-químicos por vía de administración

**Ro5 / oral** (base, Lipinski): MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10;
buena biodisponibilidad oral asociada a 0 < logP < 3. Lipinski, C.A. et
al. — regla de los 5.

**Inyectable**: Log D(pH 7.4) entre 1 y 3 para buen balance
solubilidad/permeabilidad.
[Lipophilicity – Cambridge MedChem Consulting](https://cambridgemedchemconsulting.com/lipophilicity/)

**Oftálmico** (Karami et al. 2022, *J Ocul Pharmacol Ther*): clogD(pH 7.4)
≤ 4.0, TPSA ≤ 250 Ų, ΔG(o/w) ≤ 20 kJ/mol (4.8 kcal/mol), solubilidad
(S_int y S_pH7.4) ≥ 1 μM.
[Eyes on Lipinski's Rule of Five (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8817695/)

**Tópico/dérmico**: MW óptimo ≲ 400 Da (rango aceptable 400–500 Da), logP
entre 1 y 3 (rango aceptable 1–4). La mayoría de moléculas activas no
cumple las tres condiciones simultáneamente — este flag es orientativo, no
una promesa de viabilidad de formulación.
[Predicting topical drug clearance from the skin (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7987642/)

**Nota de transparencia**: a diferencia de Ro5 (una sola regla ampliamente
citada), los criterios por vía de administración son heurísticos de
fuentes dispersas — cualquier reporte debe citar la referencia específica
de cada regla, nunca presentarlas con el mismo peso/autoridad que
Lipinski.

### Reglas drug-like/lead-like — el set completo y por qué se detuvo ahí

Cinco reglas pasa/no-pasa, cada una un criterio publicado real, nunca
inventado ni aproximado de memoria:

- **Lipinski Ro5**: Lipinski et al. (2001), *Adv Drug Deliv Rev* 46, 3-26.
  MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10.
- **Veber**: Veber et al. (2002), *J Med Chem* 45, 2615-2623. TPSA ≤ 140,
  rotables ≤ 10.
- **Ghose**: Ghose et al. (1999), *J Comb Chem* 1, 55-68. 160 ≤ MW ≤ 480,
  -0.4 ≤ logP ≤ 5.6, 20 ≤ átomos pesados ≤ 70, 40 ≤ AMR ≤ 130 (las cuatro
  condiciones, incluida la refractividad molar que una versión anterior
  de `adme_local()` no calculaba).
- **Egan**: Egan et al. (2000), *J Med Chem* 43, 3867-3877 (el modelo
  "egg" que el BOILED-Egg extiende). logP ≤ 5.88, TPSA ≤ 131.6.
- **Oprea (lead-likeness)**: Oprea (2000), *J Comput Aided Mol Des* 14,
  251-264. HBD < 2, 2 < HBA < 10, 2 < rotables < 8, 1 < anillos < 4 — la
  regla estándar de la literatura para "lead-like", más angosta que
  drug-likeness a propósito.

**Deliberadamente no implementadas**, decisión explícita con Uriel en
chat (2026-08-11): Hughes et al. (2008, *Bioorg Med Chem Lett* 18,
4872-4875), Ritchie & Macdonald (2009, *Drug Discov Today* 14, 1011-1020)
y Lovering et al. (2009, *J Med Chem* 52, 6752-6756) son hallazgos
**correlacionales** en sus papers originales (toxicidad/developability/
éxito clínico vs. una propiedad), no reglas pasa/no-pasa con un cutoff
que la fuente misma defina — implementarlas como corte duro habría
significado inventar un umbral. Muegge (2001, *J Med Chem* 44, 1841-1846)
es una regla pasa/no-pasa real pero sus umbrales exactos no se pudieron
reconfirmar contra la fuente primaria (paywall) — diferida, no
transcrita de memoria.

`n_rings_approx` (usado por `oprea_pass`) es una heurística regex sobre
pares de dígitos de cierre de anillo del SMILES, no percepción SSSR real
de CDK — `rcdk` no expone un wrapper de alto nivel para eso; ver el
código de `.smiles_ring_count_approx()` para el detalle exacto.

---

## 3. Familia `tox_*` — toxicidad

Familia aparte de ADME, sin opción de corte duro en la API. En experiencia
reportada, muchos compuestos salen "tóxicos" en estos modelos por cumplir
su efecto farmacológico real (vasodilatadores, vasoconstrictores...), no
por ser genuinamente indeseables. Que `tox_*` no tenga siquiera un
parámetro `hard_cutoff` (a diferencia de `adme_filter()`) es una decisión
de diseño, no solo de documentación.

```r
tox_local(proj, compounds = NULL, alert_sets = c("pains", "brenk"))
# Alertas estructurales vía coincidencia de subestructura (SMARTS) sobre
# las listas publicadas de PAINS y Brenk, usando rcdk::matches(). Ambos
# sets completos y empaquetados.
# Output: proj$tox_local (compound_id, alert_set, alert_name, smarts, matched)

tox_safetyome(proj, compounds = NULL)
# Anotación a nivel de BLANCO (no de compuesto): dado el uniprot_id de cada
# blanco predicho (proj$targets_imported), marca si es uno de los ~500
# genes del "core panel" de Safetyome (Liu et al. 2026) y para qué
# sistema de órganos MedDRA. Mapeo UniProt -> símbolo génico vía
# clusterProfiler::bitr() + org.Hs.eg.db (mismo mecanismo que ya usa
# network_enrich() para UniProt -> Entrez). Necesita proj$targets_imported
# ya poblado (targets_import()/targets_import_batch()).
# Output: proj$tox_safetyome (compound_id, uniprot_id, gene_symbol,
#   in_core_panel, organ_system, safetyome_scaled_score,
#   tau_score_percentile, conservation_score_percentile,
#   safetyome_scaled_score_percentile, median_tau_conservation_score,
#   median_all_scores, source)

tox_import(proj, path, platform = c("admetlab", "swissadme", "other"), column_map = NULL)
# hERG, hepatotoxicidad, etc. desde plataforma externa, mismo patrón de
# import que adme_import().
# Output: proj$tox_imported (compound_id, propiedad, valor, fuente, fecha_import)

tox_report(proj)
# Solo lectura, combina tox_local/tox_safetyome/tox_imported por
# compuesto. Nunca produce un pass/fail. Siempre incluye la nota fija:
# "un flag de toxicidad puede reflejar actividad farmacológica real, no solo riesgo"
```

**PAINS**: Baell, J.B. & Holloway, G.A. (2010), *J. Med. Chem.* 53(7),
2719-2740 — filtros SMARTS derivados de compuestos "frequent hitters" en
screenings bioquímicos (rodaninas, catecoles, quinonas entre los peores).
`inst/extdata/pains_smarts.csv` trae el set WEHI completo (480 filtros),
transcrito de `Data/Pains/wehi_pains.csv` de RDKit (BSD-3-Clause).
[New Substructure Filters for PAINS](https://research.monash.edu/en/publications/new-substructure-filters-for-removal-of-pan-assay-interference-co/)

**Brenk** (Brenk, R. et al. 2008, *ChemMedChem* 3, 435-444,
doi:10.1002/cmdc.200700139): lista de subestructuras indeseables para
librerías de screening (nitro = mutagénico, sulfatos/fosfatos = mala
farmacocinética, 2-halopiridinas/tioles = reactivos).
`inst/extdata/brenk_smarts.csv` trae el **set completo de 105 alertas**,
ya empaquetado. Provenencia (doble validación antes de empaquetar, ver
`DEVLOG.md` 2026-08-11 y `TESTING_GUIDE.Rmd` para el script de
reproducción): el conteo y los nombres de alerta se confirmaron contra el
`FilterCatalogs.BRENK` compilado de RDKit (105 entradas), el texto SMARTS
exacto se tomó de `alert_collection.csv` de PatWalters/rd_filters (MIT;
sus filas etiquetadas `"Dundee"` — el grupo de Brenk era la Dundee Drug
Discovery Unit), y las 105 patrones se re-parsearon de forma
independiente con `rdkit.Chem.MolFromSmarts()` para confirmar que son
válidos. Nada transcrito de memoria.
[TeachOpenCADD T003](https://projects.volkamerlab.org/teachopencadd/talktorials/T003_compound_unwanted_substructures.html)

**Safetyome** (Liu, X. et al. 2026, *Toxicological Sciences* 209(3),
kfag021, doi:10.1093/toxsci/kfag021): a diferencia de PAINS/Brenk (alertas
sobre la *estructura del compuesto*), esto anota el *blanco molecular* --
dado un `uniprot_id` predicho por `targets_*`, marca si corresponde a uno
de los ~500 genes del "core panel" que el propio pipeline del paper
identificó con la señal de riesgo sistémico multi-fenotipo más fuerte
(especificidad tisular + conservación evolutiva + densidad de evidencia de
genética humana, sobre 22 System Organ Class de MedDRA).
`inst/extdata/safetyome_core_panel.csv` viene de la Tabla Suplementaria 4
del paper (Uriel la descargó directo de Oxford Academic, ver carpeta
`Safetyome Data/` del proyecto) -- 500 genes únicos (507 filas; 7 genes
aparecen dos veces bajo dos sistemas de órganos distintos). El catálogo
completo del paper (~11,300 genes, Tabla Suplementaria 2) **no** se
empaquetó -- a esa escala deja de ser un panel curado de alerta y pasa a
ser "¿existe este gen en el genoma humano?" para casi cualquier lista de
blancos real; si hiciera falta más cobertura más adelante es un cambio de
config, no un rediseño (ver más abajo, sección Planeado / no implementado).

**Mecanismo técnico**: `rcdk` expone `matches(query, target,
return.matches = FALSE)`, búsqueda de subestructura vía SMARTS sobre CDK
(Java) — el mismo mecanismo que RDKit usa para PAINS en Python, ya
disponible nativo en una dependencia que el paquete ya usa para todo lo
demás. No hace falta RDKit ni `reticulate`.
[rcdk::matches() docs](https://www.rdocumentation.org/packages/rcdk/versions/3.4.3/topics/matches)

**Lo que no es calculable local sin modelo entrenado**: hERG,
hepatotoxicidad, mutagenicidad cuantitativa (Ames) — requieren modelos
entrenados sobre datasets curados, quedan exclusivamente en
`tox_import()`.

---

## 4. Familia `targets_*` — predicción de blancos (parcial)

Un predictor propio (estilo SwissTargetPrediction/SuperPred) queda fuera
de alcance por ahora: el costo no está en el algoritmo (regresión
logística sobre Morgan fingerprints, 1-2 días) sino en replicar con
exactitud la curación de ChEMBL (semanas).

```r
targets_import(
  proj, path, platform = c("swisstargetprediction","superpred","other"),
  target_col,                    # UniProt ID -- siempre obligatorio
  probability_col,               # siempre obligatorio
  confidence_col = NULL,         # opcional, no todas las plataformas lo dan
  id_from = c("filename","column"),
  compound_col = NULL            # solo si id_from = "column"
)
# id_from="filename" (default): estas plataformas solo aceptan un compuesto
# por corrida -- se pide bajar cada resultado como "Targets<pubchem_id>.csv"
# (SIN separador, ej. Targets5280443.csv -- convención heredada del
# cookbook de referencia) y el compound_id se extrae del nombre del archivo.
# platform = "superpred" trae preset de columnas confirmado contra archivos
# reales: target_col = "UniProt ID", probability_col = "Probability",
# confidence_col = "Model accuracy". "swisstargetprediction"/"other" no
# tienen preset -- target_col/probability_col siempre explícitos.
# Output: proj$targets_imported -> tibble(compound_id, uniprot_id, probability,
#                                          confidence, source, import_date)

targets_import_batch(proj, dir, platform, target_col, probability_col, confidence_col = NULL)
# Lee todos los "Targets*.csv" de una carpeta, valida el patrón de nombre,
# llama targets_import() por archivo, unifica en una sola tabla.
# Un archivo malo NUNCA tumba el resto del batch (fix 2026-08-16, ver
# DEVLOG.md -- bug real contra Chilcuague: un export de SuperPred sin
# columna "Probability" truena targets_import(), y antes eso hacía perder
# TODOS los archivos ya importados en la misma llamada) -- mismo patrón
# tryCatch-por-item que network_pathview() usa por vía de KEGG.
# targets_import() sigue abortando fuerte en una llamada directa (un solo
# archivo malo ahí SÍ debe frenar todo); solo el orquestador de batch se
# volvió resiliente.
# Output: proj$targets_import_batch_log -> tibble(path, pubchem_id,
#   compound_id, compound_name, ok, reason [imported/skipped_naming/
#   import_failed], message) -- UNA fila por archivo encontrado, éxito o
#   no, compound_name resuelto desde compounds(proj) para que un archivo
#   excluido sea identificable por nombre, no solo por CID crudo.

targets_disease_filter(
  proj, disease, source = "open_targets", min_score = NULL,
  fetch_mode = c("warn_and_cache", "abort")
)
# Anota cada UniProt ID distinto de targets_imported con su score de
# asociación a UNA enfermedad, vía Open Targets GraphQL v4
# (api.platform.opentargets.org/api/v4/graphql, sin API key). Un solo
# disease por llamada, a propósito.
# Output: proj$targets_disease -> tibble(compound_id, target_id, disease_id,
#                                         association_score, evidence)
```

**Resolución de IDs**: tanto `disease` (texto libre o ID de ontología ya
válido) como cada `uniprot_id` se resuelven a IDs canónicos de Open
Targets (EFO / Ensembl gene ID) vía `mapIds`, restringido a un tipo de
entidad por llamada. **Usar texto libre (`"hypertension"`), no un ID de
ontología a mano, salvo que sepas que está vigente** — IDs pueden quedar
obsoletos (p. ej. `EFO_0000537`, reemplazado por `MONDO_0005044` desde EFO
3.88.0); texto libre siempre resuelve al término vigente. Pares sin
ninguna asociación (`score = NA`) se descartan y quedan logueados; con
`min_score`, los pares por debajo también se descartan y logean por
separado. Cada corrida sobreescribe `proj$targets_disease` — guardar el
resultado en una variable propia si se quiere conservar una versión previa
antes de volver a llamar la función.

**Por qué Open Targets y no GeneCards**: GeneCards prohíbe scraping y
requiere licencia — mal candidato para dependencia default de un paquete
open-source. Open Targets tiene API GraphQL abierta y gratuita, diseñada
para asociación blanco↔enfermedad con score de evidencia.
[GeneCards Suite Terms of Use](https://www.lifemapsc.com/genecards-suite-terms-of-use/)

**Qué falta de esta familia**: `targets_bipartite()` (predicción propia
vía `netpredictor`, RWR/NBI) y `targets_consensus()` (fusión de fuentes) —
diseñadas, sin implementar, ver más abajo, 1.1. También
`disease_genes_import()` (alternativa manual a Open Targets vía exports de
GeneCards) — ver más abajo, 1.7.

---

## 5. Familia `network_*` — análisis de red (núcleo systems-biology)

12 funciones, la mayoría operan **por condición** (una red por columna de
`proj$binarized`, mismo fork que ya estableció `prep_binarize()`).

```r
network_build(proj, condition = NULL, target_source = c("imported","consensus","bipartite"), min_score = NULL)
# Grafo tripartito compuesto-blanco-vía vía igraph. condition=NULL corre
# todas. target_source default es "imported" (consensus/bipartite no
# existen todavía -- ver más abajo, 1.1, piden explícitos dan error claro).
# Deduplica aristas (compound_id, uniprot_id), se queda con la de mayor
# weight. Un igraph no es un CSV: guarda network_edges (tabla larga, CSV)
# como fuente de verdad, cachea un igraph por condición solo por
# performance -- borrar el cache siempre es seguro.
# Output: proj$network_edges -> tibble(condition, compound_id, uniprot_id,
#                                       weight, disease_association_score)

network_enrich(proj, condition = NULL, db = c("go","reactome","kegg"))
# Mapea UniProt -> Entrez (clusterProfiler::bitr() + org.Hs.eg.db) una sola
# vez, pasa Entrez a clusterProfiler::enrichGO()/enrichKEGG() o
# ReactomePA::enrichPathway() (solo acepta Entrez, sin keyType propio).
# db = "kegg" necesita internet (API REST de KEGG en cada llamada, no hay
# paquete de datos KEGG local vigente); "go"/"reactome" son 100% locales
# una vez instalados los paquetes de Bioconductor.
# Output: proj$network_enrichment -> tibble(condition, db, ID, Description,
#                                            geneID, p.adjust, ...)

network_pathview(proj, condition = NULL, gene_score = c("max_weight","mean_weight","n_compounds"), pathway_id = NULL, top_n_pathways = 10)
# Mapas KEGG renderizados (pathview::pathview()) para las top_n_pathways
# vías de network_enrich(db="kegg") de UNA condición, coloreadas por un
# score por blanco/Entrez derivado de network_edges (gene_score="max_weight"
# por defecto -- la probabilidad de import más alta entre los compuestos
# que tocan ese blanco en esa condición). Reproduce (no adivina) el paso
# "MAPAS DE KEGG AUTOMATIZADOS" que Juanjo ya hacía a mano por condición en
# su propio análisis pre-patliR (Farmacologia de Redes Juanjo/Redes/
# <condicion>/<condicion>.Rmd) -- misma llamada real a pathview() verificada
# ahí, mismos defaults de color (low="white", mid="yellow", high="red",
# distintos de los propios de pathview). A diferencia del resto de
# network_*, no hay salida "plot object" -- pathview() escribe PNG como
# efecto secundario (sin argumento de directorio de salida propio; esta
# función hace setwd() temporal a out_dir, restaurado con on.exit()) y
# cachea sus .xml/.png de entrada bajo cacheDir(proj)/kegg_pathview.
# Un pathway que falla en renderizar no frena a los demás (tryCatch por
# vía, igual que el loop original de Juanjo), se loguea y queda con ok=FALSE.
# Output: proj$kegg_pathview_log -> tibble(condition, pathway_id, description,
#                                           gene_score, n_genes_mapped, ok,
#                                           path, message)

network_centrality(proj, condition = NULL, measures = c("degree","betweenness","hub_score"))
# Para TODOS los nodos (compuestos y blancos por igual). Sin ponderar por
# weight (esa columna es probabilidad de predicción, no distancia/costo --
# ponderar de más invertiría el sentido).

network_hub_penalty(proj, condition = NULL)
# score_ajustado = centralidad_cruda * log(N_compuestos_total / N_compuestos_que_tocan_este_blanco)
# Autocontenida -- no depende de que network_centrality() haya corrido.

network_module_robustness(proj, condition = NULL, clustering = c("hdbscan"), min_module_size = 2, seed = NULL)
# Clustering: dbscan::hdbscan() sobre distancias de camino más corto.
# Grafos chicos (< 2*min_module_size nodos) tratan todo el grafo como un
# solo módulo, logueado. Percolación: ataque dirigido (remueve el nodo de
# mayor grado repetidamente), R-index = área bajo la curva de fragmentación
# (Schneider et al. 2011).

network_motifs(proj, condition = NULL)
# Sobre una capa dirigida nueva (compuesto -> blanco -> {vía, enfermedad},
# ver network_layers.R), NO sobre el grafo bipartito no dirigido de
# network_build(). Esa capa es un DAG por construcción -- nunca va a
# encontrar loops de feedback (esperado, no bug). Motivos enumerados por
# fuerza bruta (qué tripleta exacta de nodos), no igraph::motifs().

network_degeneracy(proj, condition = NULL)
# Agrupa compuestos por solapamiento de VÍAS (no de blancos directos) --
# "degeneracy" en sentido de biología de sistemas (Edelman & Gally 2001):
# blancos distintos, vías convergentes. Requiere network_enrich() ya
# corrido para la condición -- falla con mensaje claro si no.
# degeneracy_score = pathway_jaccard * (1 - target_jaccard).
# Output: proj$network_degeneracy -> tibble(condition, compound_a, compound_b,
#                                            target_jaccard, pathway_jaccard,
#                                            degeneracy_score)

network_proximity(proj, condition = NULL, disease = NULL, n_random = 1000, seed = NULL)
# Requiere proj$targets_disease. Distancia topológica "closest" (Guney et
# al. 2016) entre el módulo de blancos de CADA compuesto y el módulo de
# genes de la enfermedad, sobre el interactoma de STRINGdb (sin ponderar,
# weights = NA), contra modelo nulo de grado preservado.
# Output: proj$network_proximity -> tibble(condition, compound_id, disease_id,
#                                           d_observed, d_random_mean,
#                                           d_random_sd, z_score, seed_usado)

network_synergy(proj, condition = NULL, disease = NULL, pairs = c("rank_top","all"), top_n = 10)
# Requiere network_proximity() ya corrido para (condition, disease).
# synergy_score = complementarity * joint_closeness, donde complementarity
# = 1 - target_jaccard y joint_closeness = -promedio(z_score_a, z_score_b).
# pairs="rank_top" (ranking interino por z_score, ver más abajo, 1.4, para
# cuándo esto se reemplaza por rank_candidates()) o "all" (todos los pares
# sin rankear, para condiciones chicas).
# Output: proj$network_synergy -> tibble(condition, disease_id, compound_a,
#                                         compound_b, target_jaccard,
#                                         complementarity, z_score_a, z_score_b,
#                                         joint_closeness, synergy_score, pairs_mode)

network_bowtie(proj, condition = NULL, version = "12.0", actions_version = "11.0")
# Propiedad de TODA la red, no por condición -- calculado una vez por
# (species, version), cacheado. Sustrato: canal de acciones dirigidas de
# STRING (protein.actions, activation/inhibition con dirección conocida --
# STRINGdb no expone este archivo, se descarga y parsea directo). Solo
# is_directional & a_is_acting se vuelven aristas. Taxonomía de Broder et
# al. 2000 colapsada a 4 categorías.
# Output: proj$network_bowtie -> tibble(condition, uniprot_id, bowtie_component)
#         proj$network_bowtie_summary -> tibble(species, version, n_nodes, n_edges, ...)

network_filter_proteome(proj, proteome, condition = NULL)
# Filtro simple de network_edges a un vector de UniProt IDs propio (ej. el
# proteoma de una línea celular u organismo de interés) -- parte del
# alcance oficial de la propuesta FOPER 2026 (mes 5). Implementada
# 2026-08-16, alcance acotado a propósito ("sencillo, solo es un filter"):
# produce proj$network_filtered_edges para inspección/export directo, NO
# lo conecta de vuelta a network_centrality()/network_module_robustness()/
# etc. -- eso necesitaría un scope nuevo paralelo a "condition" en
# .network_resolve_conditions()/.network_graph(), cambio bastante más
# grande que lo pedido. Tampoco hace import guiado desde archivo externo
# (solo acepta un vector de UniProt IDs en R) -- ambas extensiones quedan
# como candidatas de 1.1 si hacen falta (ver Planeado / no implementado).
# Output: proj$network_filtered_edges -> mismo esquema que network_edges,
#                                          ya filtrado a `proteome`
```

### Innovaciones incorporadas — no están en ningún cookbook de network pharmacology típico

**`network_proximity()`** — de la escuela de "network medicine" (Barabási
et al.), no de network pharmacology clásico. Mide distancia topológica en
el interactoma entre el módulo de blancos del compuesto y el módulo de
genes de la enfermedad, contra distancia esperada por azar.
Guney, E., Menche, J., Vidal, M., Barabási, A.L. (2016). "Network-based in
silico drug efficacy screening." *Nat Commun* 7, 10331.
[Nature Communications](https://www.nature.com/articles/ncomms10331)

**`network_synergy()`** — extiende `network_degeneracy()` hacia pares con
blancos *complementarios* pero ambos cercanos al módulo de enfermedad, el
patrón que la literatura asocia con sinergia real (no redundancia).
Cheng, F. et al. (2019). "Network-based prediction of drug combinations."
*Nat Commun* 10, 1197.
[Nature Communications](https://www.nature.com/articles/s41467-019-09186-x)

**`network_pathview()`** — reproduce (verificado contra el `.Rmd` real de
Juanjo, no adivinado) el paso de mapas KEGG que su propio análisis
pre-`patliR` ya hacía a mano por condición.
Luo, W. & Brouwer, C. (2013). "Pathview: an R/Bioconductor package for
pathway-based data integration and visualization." *Bioinformatics*
29(14), 1830-1831. doi:10.1093/bioinformatics/btt285.

### Requisito para `plot_*` (resuelto en la sección 6)

Cada análisis de red (centralidad, bowtie, motifs, robustez modular)
necesita su propia vista, no un solo `plot_network()` genérico — redes de
compuesto-blanco-vía se saturan visualmente muy rápido.

---

## 6. Familia `plot_*` — visualización

16 funciones. Todas siguen el mismo patrón: `engine = c("static",
"ggiraph")` donde aplica, `save`/`out_dir`/`width`/`height`/`dpi`, log
propio en `patliRResults()`. `save = TRUE` por defecto, escribe PNG real
bajo `projectDir(proj)/plots/`; `save = FALSE` para explorar sin tocar
disco.

```r
plot_admet_radar(proj, compounds = NULL, engine = c("ggiraph","static"), save = TRUE, out_dir = NULL)
# Un ggplot/girafe por compuesto (lista nombrada por compound_id) -- no
# facet_wrap(), ilegible con 8+ compuestos. Trigonometría manual
# (x = r*sin(theta), y = r*cos(theta)) + coord_fixed(), NO coord_polar()
# (da lados curvos y desalineados). Ejes: LIPO (XLOGP3), SIZE (MW), POLAR
# (TPSA), INSOLU (Log S, ESOL de Delaney 2004), INSATU (Fsp3 aproximado
# por regex sobre SMILES -- sobreestima para sp2 no aromáticos, la pieza
# más floja del gráfico), FLEX (enlaces rotables). Zona rosa = rango
# óptimo publicado (Daina, Michielin & Zoete, 2017, Sci. Rep. 7:42717).
# Guarda dos variantes por compuesto: _labeled y _plain.
# proj actualizado viaja como attr(result, "proj") -- excepción a
# propósito, el valor de retorno principal acá es el gráfico, no proj.

plot_boiled_egg(proj, compounds = NULL, engine = c("ggiraph","static"), save = TRUE, out_dir = NULL)
# UN solo scatter con todos los compuestos (a diferencia del radar, esa
# comparación conjunta es el punto del gráfico). Elipses GIA/BBB reales
# (ver familia adme_*). Guarda boiled_egg_labeled.png (geom_text() simple,
# sin ggrepel -- etiquetas pueden pisarse con muchos compuestos juntos) y
# boiled_egg_plain.png.

plot_structure2d(proj, compounds = NULL, save = TRUE, out_dir = NULL)
# Grid de paneles con los PNG ya generados por prep_structure2d() (los
# relee vía png::readPNG(), R base no trae lector de PNG). Necesita
# patchwork + png.

plot_network_layers(proj, condition = NULL, layers = c("compound","target"), max_pathways = 30, top_hub_n = 15, layout = "fr")
# Red en capas (compound -> target -> {pathway, disease}) de
# .network_layered_graph()/.network_layered_graph_multi(). Estática
# (ggiraph + salida estática), no visNetwork/networkD3 -- las figuras de
# referencia son estáticas de alta resolución para publicación, no
# exploración en vivo. Sin ggraph -- layout con igraph::layout_with_fr()
# (niter = 2000), armado a mano con geom_segment()/geom_point()/geom_text().
# condition = NULL (default) agrupa TODAS las condiciones ya construidas
# en un solo grafo (.network_layered_graph_multi()); uno o más nombres de
# condición da el modo acotado. layers default es solo compound+target --
# pathway/disease son opt-in (con muchos términos GO/Reactome, si no se
# recortan por max_pathways, ahogan la estructura real -- ver DEVLOG.md).
# Etiquetas resueltas (nombre de compuesto, símbolo de gen, Description de
# vía) solo para los top_hub_n nodos de mayor grado. Aristas derivadas
# (transitivas) se dibujan punteadas.

plot_network_degeneracy(proj, condition = NULL, min_degeneracy = 0.3, ...)
# Reusa el layout completo de plot_network_layers() y suma un arco curvo
# por par de compuestos de network_degeneracy() con degeneracy_score >=
# min_degeneracy.

plot_centrality(proj, condition = NULL, measure = c("degree","betweenness","hub_score"))
# Barras de network_centrality(), facet por categoría específico/promiscuo
# cuando esté disponible bias_homogeneity (ver sección 7 / más abajo, 1.3).

plot_robustness(proj, condition = NULL)
# Curva de percolación + R-index de network_module_robustness().
# Schneider, C.M. et al. (2011). "Mitigation of malicious attacks on
# networks." PNAS 108(10), 3838-41.
# https://www.pnas.org/doi/10.1073/pnas.1009440108

plot_proximity(proj, condition = NULL, disease = NULL)
# Lollipop de z-scores por compuesto (network_proximity()), bandas de
# referencia en z = 0, +-1.96 -- NO un histograma de permutación: solo se
# persiste el resumen del modelo nulo (media/sd), no las n_random
# distancias re-muestreadas.

plot_synergy(proj, condition = NULL)
# Scatter complementarity x joint_closeness de network_synergy(), tamaño/
# color por synergy_score, top pares etiquetados.

plot_bowtie(proj, condition = NULL)
# Alluvial estático (ggalluvial) compound -> bowtie_component, no el
# Sankey/networkD3 interactivo del catálogo original -- mismo principio
# "estático primero, Suggests mínimo".

plot_venn(proj, condition = NULL, sets = c("compound_targets","disease_targets"))
# ggVennDiagram, 2 conjuntos -- compuesto ∩ enfermedad.

plot_upset(proj, condition = NULL)
# Barra + matriz de puntos con eje x compartido, armado a mano con
# patchwork (no UpSetR/ComplexUpset). Requiere >= 2 condiciones.

plot_heatmap(proj, condition = NULL, what = c("compound_target","compound_condition"),
  save = TRUE, out_dir = NULL, width = NULL, height = NULL)
# MIGRADO 2026-08-16 al pheatmap() real del cookbook (ver DEVLOG.md
# 2026-08-16) -- clustering jerárquico real en ambos ejes
# (desactivado solo, con aviso, en el eje con <2 filas/columnas -- hclust()
# no puede agrupar un solo elemento), valores impresos sobre la escala de
# color. Sin `engine`/dpi -- pheatmap es grid/base graphics, no ggplot2,
# ggiraph no tiene nada de qué agarrarse acá; save=TRUE ahora escribe PDF
# (no PNG) vía el propio argumento filename= de pheatmap() (confirmado
# contra su código fuente real, heatmap_motor(), no asumido).

plot_gochord(proj, condition = NULL, top_n_terms = 10, engine = c("static","ggiraph"))
# MIGRADO 2026-08-16 al GOplot::GOChord() real (cintas curvas, no líneas
# rectas) -- ver DEVLOG.md 2026-08-16. GOChord() resultó
# estar construido sobre ggplot2 por dentro (confirmado contra su código
# fuente real) -- el contrato engine=c("static","ggiraph") se conserva
# intacto, a diferencia de plot_heatmap(). SIEMPRE llama con nlfc = 0 (sin
# columna logFC): el cookbook usa un placeholder sintético
# (seq(-1,1,...), sin relación a datos reales) que a propósito NO se
# copió -- patliR no tiene fold-change real por gen (es predicción de
# blancos, no expresión diferencial), inventar uno hubiera sido engañoso.

plot_chemical_space(proj, condition = NULL, compound_ids = NULL,
  method = c("pca","umap"), color_by = "family", show_hulls = TRUE,
  seed = NULL, engine = c("ggiraph","static"), save = TRUE, out_dir = NULL)
# Reducción de dimensionalidad (PCA vía stats::prcomp() base R, o UMAP
# opcional vía el paquete umap) sobre los descriptores de adme_local()
# (mw, logp, hbd, hba, tpsa, rotatable_bonds, amr, fraction_csp3_approx,
# aromatic_proportion_approx, n_rings_approx). color_by = "family"
# (default) colorea por proj$compounds_classified$pathway -- aborta con
# mensaje claro si compounds_classify() no se ha corrido todavía;
# cualquier otra columna numérica/categórica de adme_local también sirve.
# show_hulls = TRUE dibuja un polígono convexo (grDevices::chull(), sin
# dependencia nueva) por grupo cuando color_by es categórico, más una
# etiqueta de texto en negritas en el centroide de cada familia (además
# de la leyenda). Estilo visual con theme_void() + ejes tipo flecha que
# cruzan el origen (geom_segment con arrow()), reconstruyendo a mano el
# estilo de las figuras de "chemical space" con ejes centrados en el
# origen que usa la literatura (Reymond, J.-L. y colegas; NO se copió
# ninguna figura, se reconstruyó el estilo). overlay_targets (proyectar
# blancos sobre el mismo espacio) quedó deliberadamente fuera -- ver más
# abajo, "Extensiones pendientes de funciones ya implementadas".

plot_adme_upset(proj, top_n = 15, save = TRUE, out_dir = NULL)
# UpSet armado a mano (mismo patrón que plot_upset(), sin UpSetR/
# ComplexUpset): barra de tamaño de intersección + matriz de puntos con
# eje x compartido vía patchwork, pero la "cosa que se intersecta" acá es
# el conjunto de reglas drug-like/lead-like que cada compuesto pasa
# (ro5/veber/ghose/egan/oprea de adme_filter()), no condiciones. Excluye
# las filas route_* (esas ya son booleanas post-decisión, no reglas
# individuales). top_n limita a las intersecciones más frecuentes.
```

### Dependencias

`Suggests`: `ggplot2`, `ggiraph` (motor base de toda la familia),
`ggalluvial` (solo `plot_bowtie()` — geometría de cintas curvas no es
razonable a mano), `ggVennDiagram` (solo `plot_venn()` — áreas
proporcionales), `patchwork` (`plot_upset()`, `plot_structure2d()`,
`plot_adme_upset()` — combinar paneles `ggplot` separados), `png` (solo
`plot_structure2d()`), `umap` (solo `plot_chemical_space(method =
"umap")` — PCA no necesita paquete nuevo, ya usa `stats::prcomp()`),
`pheatmap` (solo `plot_heatmap()`, real desde 2026-08-16), `GOplot` (solo
`plot_gochord()`, real desde 2026-08-16).

```r
install.packages(c("ggalluvial", "ggVennDiagram", "patchwork", "png", "umap", "pheatmap", "GOplot"))
```

**Deliberadamente sin agregar** (reimplementado a mano en `ggplot2` puro,
mismo principio de `Suggests` mínimo): `ggrepel` (`plot_boiled_egg()`),
`ggraph` (`plot_network_layers()`), `UpSetR`/`ComplexUpset`
(`plot_upset()`). `pheatmap`/`GOplot` salieron de esta lista el
2026-08-16 -- ver arriba, migración real ya hecha.

### Qué falta de esta familia

`plot_rank()` y `plot_kegg_binding()` — bloqueadas por `rank_candidates()`
y por la familia de docking, ninguna tiene sustrato real todavía. Ver
más abajo, 1.6.

---

## 7. Familia `bias_*` — auditoría de sesgo de dominio (parcial)

```r
bias_audit(proj, check_homogeneity = TRUE, mad_threshold = 2.5, categories = NULL)
# Homogeneity check autocontenido: para cada compuesto y cada blanco de
# proj$reference_bioactivity (refdb_build()), calcula su "grado" (para un
# compuesto: cuantos blancos distintos toca; para un blanco: cuantos
# compuestos distintos lo tocan) y lo compara contra la mediana +-
# mad_threshold*MAD de su propio tipo (Leys et al. 2013). Solo el lado
# alto se marca "promiscuo" -- un grado inusualmente BAJO no es un
# problema de sesgo, solo un compuesto/blanco poco estudiado. NUNCA
# elimina, solo categoriza. categories (enriquecimiento MeSH/DO vs.
# background STRING/DrugBank/reference_db) NO esta implementado --
# necesita una fuente externa confirmada antes de codearse (ver
# más abajo, 1.3); pedirlo explicito da error claro, no un resultado
# adivinado.
# Output: proj$bias_homogeneity -> tibble(id, tipo [compound/target],
#                                          frecuencia_global_refdb,
#                                          mad_score, categoria [promiscuo/normal])

bias_reweight(proj)
# Mismo log-ratio que network_hub_penalty() (degree * log(N_universo /
# degree)) pero a nivel de TODA la reference_db, no por condicion --
# N_universo es "total de blancos distintos" (para reponderar un
# compuesto) o "total de compuestos distintos" (para reponderar un
# blanco), leido directo de bias_homogeneity (una fila por id distinto de
# cada tipo), sin releer reference_bioactivity. Solo se aplica descuento a
# categoria=="promiscuo" -- el resto conserva score_ajustado == score_crudo,
# nunca se sobreescribe el crudo.
# Output: proj$bias_reweighted -> tibble(id, tipo, score_crudo, score_ajustado, categoria)

bias_report(proj)
# Solo lectura (mismo contrato que tox_report()): nada se escribe a proj
# ni a disco. Tabla de "promiscuos" (categoria=="promiscuo" de
# bias_reweighted, ordenada por score_ajustado ascendente), SIEMPRE
# presente (vacia + warning si bias_audit()/bias_reweight() no han
# corrido, nunca un error duro) + nota fija de que un flag "promiscuo" no
# invalida la asociacion, solo marca sesgo de curacion potencial.
```

Implementado 2026-08-16: ver `R/bias_audit.R` y la entrada de `DEVLOG.md`
2026-08-16. Deviaciones deliberadas de la firma tentativa original (más
abajo, 1.3): sin parámetro `condition` (el homogeneity check es una
propiedad de toda la reference_db, no por condición de red) ni
`background` (solo tiene sentido para el enriquecimiento categórico, que
no está implementado todavía) -- cargar parámetros sin uso real hasta que
esa parte exista hubiera sido código muerto desde el día uno.

**Sustento**: Zhang, Q. (2026), *Frontiers in Pharmacology* 17:1748478 —
análisis de 1,038 estudios publicados de network pharmacology, encuentra
"homogeneidad" (mismas moléculas/blancos hub top sin importar remedio o
enfermedad). Mecanismo cuantificado en Cell Genomics 2023 (PMC10363916):
en STRING, el grado de un nodo predice si es blanco de fármaco conocido
con AUC=77.6% — sesgo de curación, no relevancia biológica real.
Referencia técnica del método de outlier: Leys et al. 2013. Método
alternativo evaluado y no adoptado: MeSHOP (PMC3654871).

### Qué falta de esta familia

El enriquecimiento categórico (`categories = "mesh_disease"/"do"` vs.
`background = "string"/"drugbank"/"reference_db"`) -- la otra mitad del
diseño original de `bias_audit()`, ver más abajo, 1.3. Necesita que Uriel
confirme una fuente real de categorías de enfermedad por blanco (MeSH o
Disease Ontology, cubriendo TODAS las enfermedades, no solo la que
`targets_disease_filter()` trae por llamada) antes de codearse -- mismo
principio de "nunca adivinar una API externa" que ya se aplicó a
`coconut_fetch()`/`network_kegg_complete()`. Pedirlo explícitamente
(`categories` no `NULL`) da un error claro, no un resultado inventado.

También pendiente, ahora desbloqueado por esta familia:
`rank_candidates()` (más abajo, 1.4 -- ya puede usar `score_ajustado` de
`bias_reweight()` como una de sus cuatro entradas ponderadas) y el facet
específico/promiscuo de `plot_centrality()` (sección 6, nota sobre
`bias_homogeneity`).

---

## 8. `patliR_export_llm()` — exportación del proyecto para análisis por LLM

```r
patliR_export_llm(proj, out_file = NULL, max_rows = 200)
```

Implementada 2026-08-16, pedido directo de Uriel en chat sin diseño
previo. Serializa `compounds()`, cada entrada de `patliRResults()`, y
`projectLog()` completo a texto plano (tablas en formato CSV bajo
encabezados Markdown) — pensado para pegarse o subirse directo a un chat
de LLM, sin pasar por R. Deliberadamente genérico sobre `patliRResults()`
(itera lo que exista, no una lista fija por familia) en vez de un resumen
curado por sección — esa curaduría es exactamente el trabajo de
`report_generate()` (ver más abajo, 1.5, sin implementar) y del
"narrador" (2.7, sin diseñar): esta función es el volcado completo de
datos crudos para que el LLM razone, no un resumen pre-escrito. Nunca muta
`proj` (mismo contrato de solo lectura que `tox_report()`/`bias_report()`);
tablas largas se truncan (con aviso explícito, nunca en silencio) a
`max_rows` para no saturar la ventana de contexto del LLM.

---

## Planeado / no implementado

Todo lo de esta sección es **diseño, no implementación** — ninguna función
de acá tiene código todavía salvo que se indique lo contrario. Cada
familia/función tiene su propio espacio (firma tentativa, propósito,
output esperado, por qué), mismo formato que usan las secciones de arriba
para lo ya implementado, para que al implementarse sea una migración
directa, no una reescritura desde cero. (Esta sección vivía en un archivo
aparte, `ROADMAP.md`, hasta el 2026-08-23 — fusionada acá para que
"implementado" y "planeado" se lean como un solo documento.)

Fuentes: el diseño original del paquete (secciones que quedaron
"diseñadas pero no implementadas" tras el cierre de `network_*`/`plot_*`,
ver la bitácora de desarrollo interna), y la propuesta FOPER 2026
(documento interno del Laboratorio de Investigación Química y
Farmacológica de Productos Naturales, UAQ) más aclaraciones directas de
Uriel — la propuesta es "el spec final" en el sentido de que de ahí sigue
creciendo el alcance, no reemplaza el diseño ya cerrado.

### Estado general (alcance de 1.0 vs. 1.1)

Todo lo descrito en las secciones 1-8 de arriba está **implementado y
cerrado**. Decisiones de alcance de 1.0 confirmadas por Uriel (2026-08-16):
`targets_bipartite()`/`targets_consensus()` **descartados para 1.0** (1.1
de abajo); `tcm_import()` y los conectores regionales/validación contra
literatura de la propuesta FOPER **diferidos a 1.1** (1.2 y 2.8); `coconut_fetch()`
evaluado y también diferido a 1.1 (no se pudo confirmar el endpoint real
sin adivinar, 1.2); `rank_*`/`dock_*`/`report_*`/`plot_rank()`/
`plot_kegg_binding()` siguen sin implementar, bloqueados por dependencias
entre sí (1.4-1.6); `network_kegg_complete()` diferido a 1.1 por
presupuesto de tiempo (2.1). El enrichment categórico de `bias_audit()`
(la mitad del diseño original que no tiene código, ver sección 7) sigue
sin implementar — necesita una fuente MeSH/Disease Ontology confirmada
(1.3).

### 1.1 `targets_bipartite()` / `targets_consensus()` — descartado para 1.0 (2026-08-16)

Decisión explícita de Uriel: no se hará para 1.0. El costo real nunca
estuvo en el algoritmo (RWR/NBI son directos) sino en replicar la
curación de datos de ChEMBL con exactitud suficiente para que un
predictor propio sea confiable — inversión que no se justifica para esta
versión. `targets_consensus()` depende de tener al menos dos fuentes que
combinar (`"bipartite"` + `"imported"`); sin `targets_bipartite()` no hay
nada que consensuar, así que cae con la misma decisión.
`network_build(target_source = "consensus"|"bipartite")` sigue dando el
error de "no implementado todavía". Candidato a revisar para 1.1, no
antes.

```r
targets_bipartite(
  proj, compounds = NULL, network_source = c("chembl","reference_db"),
  method = c("rwr", "nbi"), restart_prob = 0.7,
  permutations = 999, seed = NULL
)
# Suggests: netpredictor (GPL-2, no se hereda al core). RWR/NBI sobre red
# bipartita compuesto-blanco desde bioactividad conocida (ChEMBL/refdb).
# netpredictor usa test de permutación para significancia -> seed real,
# se expone y se loguea (Seal & Wild, BMC Bioinformatics 2018).
# Output: proj$targets_bipartite -> tibble(compound_id, target_id, score,
#                                           p_value, method, seed_usado)

targets_consensus(proj, sources = c("bipartite","imported"), weights = NULL)
# weights = NULL -> promedio simple entre fuentes disponibles, expuesto.
# Output: proj$targets_consensus -> tibble(compound_id, target_id,
#                                           score_consensus, n_fuentes,
#                                           detalle_por_fuente)
```

### 1.2 `coconut_*` / `tcm_*` — diferido a 1.1 (2026-08-16)

`tcm_import()`: confirmado explícitamente por Uriel, no va en 1.0.

`coconut_fetch()`: evaluado si era "sencillo" de implementar ya. Veredicto
— no todavía, sin adivinar. La documentación Swagger de COCONUT
(`coconut.naturalproducts.net/api-documentation`) es una SPA que no se
pudo leer en el sandbox (JS renderizado, no HTML estático), y no se
encontró el repo del backend con las rutas reales (`Steinbeck-Lab/coconut`
en GitHub parece ser el frontend). La API "real y activa" sigue siendo un
hallazgo válido (viene del propio paper NAR 2025 de COCONUT 2.0
describiéndola, no de haber probado el endpoint en vivo), pero
implementar `coconut_fetch()` sin confirmar el path/parámetros/forma del
JSON de respuesta exactos sería exactamente el tipo de "adivinar una API
externa" que este proyecto evita a propósito (mismo principio ya aplicado
en `network_kegg_complete()` y en `refdb_build(sources = "coconut")`, que
hoy da error claro de "no implementado" en vez de fallar en silencio).
Candidato real para 1.1 — la siguiente sesión que lo retome debería
empezar confirmando el endpoint real antes de escribir código.

```r
coconut_fetch(proj, compounds = NULL)
# Suggests, httr2. API REST confirmada, activa (COCONUT 2.0, NAR 2025,
# compatible OpenAPI, releases mensuales). warn_and_cache si no hay conexión.
# Output: se integra a reference_bioactivity.csv / reference_compounds.csv

tcm_import(proj, path, platform = c("unitcm","herb2","other"), column_map = NULL)
# Sin API confiable confirmada para UniTCM ni HERB 2.0 -> import manual,
# mismo patrón que adme_import()/tox_import()/targets_import().
# Output: proj$tcm_imported -> tibble(compound_id, herb_name, formula_tcm,
#                                      fuente, fecha_import)
```

Fuentes evaluadas: COCONUT 2.0 (API real confirmada), UniTCM (sin API
pública confirmada), HERB 2.0 (más completo que UniTCM, tampoco API
confirmada, queda como plataforma de `tcm_import()`), NPASS (bioactividad
relevante pero sin REST limpio, descartada por ahora), LOTUS
(estructura-organismo, no bioactividad, no resuelve necesidad core),
THEOBROMA (preprint bioRxiv jun 2026, agregador con auditoría de licencia
por compuesto — vigilar, sin API estable todavía).

### 1.3 `bias_*` — enrichment categórico pendiente

Ver sección 7 arriba para lo ya implementado (homogeneity check completo)
y el "Sustento"/"Qué falta" ahí mismo para el detalle de qué falta y por
qué (necesita una fuente MeSH/Disease Ontology confirmada, cubriendo
TODAS las enfermedades, antes de codearse).

### 1.4 `rank_*` — priorización pre-docking

```r
rank_candidates(
  proj, condition = NULL,
  weights = list(centralidad = 0.25, adme = 0.25, targets = 0.25, robustez = 0.25),
  include_route_flags = TRUE
)
# Combina score_ajustado (network_hub_penalty + bias_reweight) + adme_filter
# + targets_consensus + network_module_robustness en un ranking compuesto.
# Output: proj$ranking -> tibble(compound_id, score_total, score_centralidad,
#                                 score_adme, score_targets, score_robustez,
#                                 targets_principales, route_oral,
#                                 route_topical, route_ophthalmic, route_injectable)
```

`route_*` siempre con el disclaimer fijo: "compatibilidad fisicoquímica
con la vía, NO es recomendación de formulación farmacéutica ni de
manufactura." Al implementarse, revisar `network_synergy(pairs =
"rank_top")`, que hoy usa un ranking interino basado en
`network_proximity()` a la espera de esta función (ver `DEVLOG.md`).

### 1.5 `dock_*` / `report_*`

```r
dock_prepare(
  proj, condition = NULL, compounds = NULL,   # default: top de rank_candidates()
  target_pdb, box_center, box_size = c(20, 20, 20),
  engine = c("vina", "diffdock"), out_dir = NULL
)
# Convierte compuestos a PDBQT (rcdk/openbabel) + config con la caja de
# búsqueda. box_center OBLIGATORIO, sin default -- no se adivina el sitio
# de unión (eso es predicción de bolsillo, P2RANK/fpocket, fuera de alcance).
# Output: proj$dock_config -> tibble(compound_id, pdbqt_path, config_path,
#                                     box_center, box_size, engine)

dock_parse(proj, results_dir, engine = c("vina", "diffdock"))
# Parsea salida externa (Vina/DiffDock) a dataframe R. Sin scoring propio.
# Output: proj$dock_results -> tibble(compound_id, pose_id, binding_affinity,
#                                      rmsd_lb, rmsd_ub, engine)
```

100% agnóstico al motor, confirmado también en la propuesta FOPER 2026:
patliR nunca corre docking, solo prepara inputs e importa resultados
externos. Fuera de alcance explícito: motor de docking en sí, visualizador
3D, dinámica molecular.

```r
report_generate(
  proj, condition = NULL, format = c("terminal", "html", "pdf", "latex", "csv"),
  sections = "all"
)
# Motor de datos único, renderizadores múltiples. Un reporte por condición.
# Secciones siempre presentes: filtrado ADME, ranking, auditoría de sesgo,
# robustez modular, log completo de decisiones/parámetros/semillas.
# Sección de docking: opcional, solo si proj$dock_results existe.
# Output: archivo(s) en cache_dir/reports/<condicion>/
```

Orden de desarrollo confirmado: `report_*` primero, `dock_*` después.

### 1.6 `plot_rank()` / `plot_kegg_binding()`

```r
plot_rank(proj, condition = NULL)                # ranking, target-anotado, facet promiscuo
plot_kegg_binding(proj, condition = NULL, color_by = c("confianza","concentracion_relativa"))
```

Bloqueadas por 1.4 (`rank_candidates()`) y por la familia de docking —
ninguna tiene sustrato real todavía para graficar. Nota ya anotada desde
`prep_binarize()`: se conserva `matrix_raw` (valor continuo, no solo
binarizado) específicamente para que `plot_kegg_binding()` pueda ofrecer
`color_by = "concentracion_relativa"` — falta definir exactamente contra
qué se normaliza (máximo de la condición, entre condiciones, z-score).

### 1.7 `disease_genes_import()` — alternativa manual a Open Targets

```r
disease_genes_import(proj, path, source = c("genecards","other"), column_map = NULL)
# Alternativa manual a targets_disease_filter(source="open_targets") --
# incorpora al mismo slot proj$targets_disease, columna `fuente` para que
# report_generate()/bias_report() muestren de dónde vino cada asociación.
# Output: proj$targets_disease -> tibble(gene_symbol, uniprot_id, disease,
#                                         relevance_score, fuente = "genecards_import")
```

Diseñada durante la revisión de higiene de 2026-07-23 (los archivos de
ejemplo `GeneCards{Hypertension,Inflammation}.csv` ya estaban
empaquetados como exports manuales del usuario, mismo principio ya
aplicado a SwissADME/SuperPred: importar un CSV que el usuario ya
descargó es válido aunque la fuente no tenga API abierta). Confirmado
contra `R/` que nunca se llegó a codear.

### 2.1 Completar la red con la topología dirigida de KEGG — diferido a 1.1 (2026-08-16)

Evaluado en la misma sesión que el filtro de proteoma/similitud (secciones
5 y 1) y descartado para 1.0 por presupuesto de tiempo, no por dificultad
conceptual: a diferencia de `network_filter_proteome()` (un filtro
autocontenido) y `compounds_similarity()` (fingerprints/Tanimoto, sin red
externa), esto necesita parsear KGML real de la API REST de KEGG (formato
XML dirigido, nunca ejercitado contra una respuesta real), decidir cacheo
local vs. pegar a la API en cada llamada, y tocar `plot_network_layers()`
para marcar visualmente nodo "propio" vs. "nativo de KEGG" — tres
preguntas de diseño reales, ninguna resoluble de forma responsable
apurado. Candidato prioritario para 1.1 — diseño ya detallado abajo, no
se perdió nada, solo no se implementó todavía.

```r
network_kegg_complete(proj, condition = NULL, pathway_id, mark_origin = TRUE)
# Firma tentativa, sujeta a cambio.
```

Hoy `network_enrich(db = "kegg")` solo trae el enrichment de vías (qué
vías están sobrerrepresentadas), no la topología interna de cada vía. Las
vías de KEGG ya vienen dirigidas (A activa/inhibe a B, no solo "A
interactúa con B" como STRING). Idea: para una vía seleccionada por el
usuario, traer el grafo dirigido completo de KEGG, marcar qué nodos son
"nuestros" (blancos predichos que ya están en la red de patliR) vs. nodos
nativos de KEGG que no estaban, y usar esa topología dirigida real para
trazar efectos downstream — ejemplo dado por Uriel: en `network_degeneracy()`
(o una variante), poder evaluar todo lo que lleva a un efecto puntual —
p. ej. vasoconstricción — siguiendo la vía dirigida real de KEGG, en vez
de solo las aristas ad-hoc que arma `.network_layered_graph()`/
`.network_layered_graph_multi()` hoy (que no distinguen dirección causal,
solo co-pertenencia a compuesto-blanco-vía).

Conecta directo con el flujo de "extensión de red" descrito en la
propuesta FOPER (módulo de red): el usuario selecciona una vía
significativa, patliR identifica qué proteínas de esa vía no están en la
red original, y vuelve a consultar STRING para completar las
interacciones faltantes — mismo patrón, pero con KEGG como fuente de
topología dirigida en vez de (o adicional a) STRING.

Preguntas de diseño abiertas: ¿el grafo dirigido de KEGG se cachea local o
se pega a la API REST en cada llamada (mismo caveat ya documentado para
`network_enrich(db="kegg")` — no hay paquete de datos KEGG vigente,
`KEGG.db` deprecado)? ¿Cómo se marca visualmente un nodo "propio" vs.
"nativo de KEGG" en `plot_network_layers()`?

**Distinto de `network_pathview()`** (sección 5, ya implementada): esa
función es solo *visualización* (pinta genes/blancos ya conocidos sobre el
mapa oficial de KEGG, no cambia ni extiende ningún grafo); esto es sobre
traer la *topología dirigida* de una vía de KEGG para extender el grafo de
patliR con ella. Dos piezas separadas, no se combinaron.

### 2.7 Módulo de interpretación asistido por IA ("narrador")

Estrictamente narra resultados ya calculados — nunca modifica el análisis
ni decide nada (no elige umbrales, no descarta compuestos, no cambia
parámetros). Todo prompt/respuesta se loguea y se declara en los reportes,
siguiendo reglas de la convocatoria FOPER. Sin diseño de implementación
todavía — qué modelo/API se usa, cómo se integra al logging existente de
`projectLog()`, qué pasa si el usuario no tiene acceso a un LLM (¿la GUI
sigue funcionando sin este módulo?).

### 2.8 Resto de alcance FOPER 2026 no cubierto arriba — diferido a 1.1 (2026-08-16)

Confirmado explícitamente por Uriel: queda pendiente para después de 1.0.

- **Conectores a bases regionales** — BIOFACQUIM, LANaPDB, contrastadas
  contra COCONUT (ver 1.2, `coconut_*`/`tcm_*`, ya en diseño).
- **Validación contra literatura** — casos de estudio contrastados contra
  la Biblioteca Digital de la Medicina Tradicional Mexicana (UNAM).
- **Cronograma de la propuesta**: 8 meses (septiembre 2026 - mayo 2027),
  hitos mensuales por módulo. Consultar
  [`docs/patliR_FOPER_2026.pdf`](docs/patliR_FOPER_2026.pdf) directamente
  para la fecha exacta de un hito puntual.

### Extensiones pendientes de funciones ya implementadas

Estas funciones tienen código real (ver sus secciones arriba); lo que
sigue son extensiones deliberadamente dejadas fuera de la primera versión,
no partes faltantes de un diseño incompleto:

- **`plot_chemical_space()` (sección 6) — overlay de blancos/proteínas.**
  La idea original (elipses o formas de blancos superpuestas sobre el
  mismo espacio) no se implementó: una proteína no tiene logP/TPSA en el
  mismo sentido que un compuesto, así que un embedding conjunto
  compuesto+blanco es una pregunta de diseño genuinamente distinta (¿qué
  eje compartido tendría sentido? ¿una matriz de perfiles de interacción
  compuesto-blanco reducida en vez de descriptores físico-químicos?) y no
  se quiso resolver a la ligera. El pedido de Uriel de colorear por
  familia (2026-08-11) ya está satisfecho; retomar el overlay de blancos
  como su propio punto de diseño si hace falta.
- **`tox_safetyome()` (sección 3) — catálogo completo (~11,300 genes) en
  vez del core panel (500).** Si el core panel se queda corto en la
  práctica, extender es cambiar qué CSV lee `.load_safetyome_core_panel()`,
  no un rediseño de `tox_safetyome()`.
- **`network_filter_proteome()` (sección 5) — conectar el filtrado de
  vuelta a `network_centrality()`/`network_module_robustness()`/etc., e
  import guiado desde archivo externo.** Necesitaría un scope nuevo
  paralelo a "condition" en `.network_resolve_conditions()`/
  `.network_graph()` — cambio bastante más grande que el filtro simple ya
  implementado.
- **`compounds_similarity()` (sección 1) — integración con
  `plot_chemical_space()`/`network_synergy()`.** Sigue siendo función
  independiente; ambas integraciones son candidatas futuras sin bloquear
  el uso standalone actual.

### Pendientes menores de higiene (no bloquean)

- Definir predicción de bolsillo (`box_center` en `dock_prepare()`, 1.5):
  queda en manos del usuario por ahora, sin wrapper propio a P2RANK/fpocket.
- La vignette (`vignettes/patliR-intro.Rmd`) y el `paper.md`/JOSS siguen
  sin escribir — pendientes antes de una eventual sumisión a Bioconductor.
- **Recordatorio operativo**: cada función nueva agregada necesita un
  `devtools::document()` corrido localmente antes de confiar en que
  `man/`/`NAMESPACE` reflejan el código real — no se puede correr `R`
  desde el entorno donde se escribe el código en estas sesiones de
  diseño, así que esto se le devuelve a Uriel como paso manual cada vez,
  no se asume hecho.
- Consolidación de arquitectura entre familias (semántica de "qué pasa en
  una segunda corrida" distinta entre `network_*`/`adme_local()` etc.,
  `plot_*` sin logging propio, boilerplate repetido en los 16 `plot_*`) —
  ver `DEVLOG.md`, 2026-08-16 (continuación 5), hallazgo de mayor
  prioridad para antes de escalar más `plot_*`/`network_*` en 1.1.

---

## Estructura de paquete

```
patliR/
├── DESCRIPTION
├── NAMESPACE            # generado por roxygen2
├── R/                   # por familia (prep_*.R, network_*.R, plot_*.R, etc.)
├── man/                 # docs generadas (load_all()/document())
├── tests/testthat/      # un test mínimo por función pública, sobre datos de ejemplo
├── vignettes/           # tutorial largo (pendiente, ver más abajo)
├── inst/extdata/        # datos de Sphaeralcea angustifolia como ejemplo empaquetado
├── docs/                # propuesta FOPER 2026 y otros documentos de referencia
├── README.md
├── LICENSE              # MIT
└── paper.md + paper.bib # artículo JOSS (pendiente)
```

Dataset de ejemplo (`inst/extdata/`): `input_single_compound.csv` (1 fila,
Apigenin — ejemplo de `prep_compound()`), `input_compound_list.csv` (8
compuestos — ejemplo de `prep_compounds()`, Escenario B),
`input_abundance_matrix.csv` (matriz GC-MS, columnas `R<n>-<CONDICION>`:
LEA-ET, LEA-AQ, FLO-ET, FLO-AQ, STE-ET — ejemplo de `prep_binarize()`,
Escenario A), `import_adme_swissadme.csv`/`import_tox_admetlab.csv`
(ejemplos de `adme_import()`/`tox_import()`),
`import_targets/Targets<CID>.csv` × 8 (ejemplo de
`targets_import_batch()`).
