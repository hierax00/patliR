# patliR -- Shiny wizard (ROADMAP.md, 2.6 `launch_app()`)
#
# This app is a thin orchestration layer: every action button below calls an
# existing, already-tested patliR function (prep_*/refdb_*/adme_*/tox_*/
# targets_*/network_*/plot_*) and stores the returned `proj` in a reactive
# value. No analysis logic lives here -- if a computation looks wrong, the
# bug is in the underlying R/*.R function, not in this file.
#
# Launched via patliR::launch_app(), which calls
# shiny::runApp(system.file("shiny", package = "patliR")). Can also be run
# directly with shiny::runApp() from this file's directory during
# development (needs `devtools::load_all()` run first in that session).

library(shiny)

has_dt <- requireNamespace("DT", quietly = TRUE)

render_table_ui <- function(id) {
  if (has_dt) DT::DTOutput(id) else tableOutput(id)
}
## `data_fn` MUST be a zero-arg function (not a bare expression) -- shiny's
## render*() functions capture their argument as an unevaluated expression
## and re-`eval()` it on every invalidation; passing a lazily-evaluated R
## promise through a wrapper function like this one would make that promise
## get forced (and memoized) only once, so the table would never update
## after the first render. A real function call has no such caching.
render_table_server <- function(output, id, data_fn) {
  if (has_dt) {
    output[[id]] <- DT::renderDT({ data_fn() }, options = list(pageLength = 10, scrollX = TRUE))
  } else {
    output[[id]] <- renderTable({ data_fn() })
  }
}

## ---- small shared helpers -------------------------------------------------

.notify_ok <- function(msg) showNotification(msg, type = "message", duration = 4)
.notify_err <- function(e) showNotification(conditionMessage(e), type = "error", duration = NULL)

## `run_step()` runs `fun` (a zero-arg closure calling the real patliR
## function), and on success replaces `rv$proj`; on failure it shows the
## real error message instead of crashing the app -- every patliR function
## already raises informative `cli::cli_abort()` errors, this just surfaces
## them in the UI instead of the R console.
run_step <- function(rv, fun, ok_msg) {
  tryCatch({
    rv$proj <- fun()
    .notify_ok(ok_msg)
    TRUE
  }, error = function(e) {
    .notify_err(e)
    FALSE
  })
}

condition_choices <- function(proj) {
  if (is.null(proj) || nrow(patliR::binarizedMatrix(proj)) == 0) return(character())
  setdiff(colnames(patliR::binarizedMatrix(proj)), "compound_id")
}

## ---- UI --------------------------------------------------------------

ui <- fluidPage(
  titlePanel("patliR -- asistente de análisis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Estado del proyecto"),
      verbatimTextOutput("project_status"),
      tags$hr(),
      h5("Registro (últimos eventos)"),
      render_table_ui("log_tail"),
      tags$hr(),
      downloadButton("dl_log", "Descargar log completo (CSV)")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "wizard",

        ## 1. Proyecto ---------------------------------------------------
        tabPanel("1. Proyecto",
          br(),
          fluidRow(
            column(6,
              wellPanel(
                h4("Crear proyecto nuevo"),
                textInput("new_project_dir", "Carpeta del proyecto", value = ""),
                textInput("new_cache_dir", "Carpeta de caché (opcional)", value = ""),
                actionButton("btn_new_project", "Crear", class = "btn-primary")
              )
            ),
            column(6,
              wellPanel(
                h4("Cargar proyecto existente"),
                textInput("load_project_dir", "Carpeta del proyecto", value = ""),
                textInput("load_cache_dir", "Carpeta de caché (opcional)", value = ""),
                actionButton("btn_load_project", "Cargar")
              )
            )
          )
        ),

        ## 2. Compuestos ---------------------------------------------------
        tabPanel("2. Compuestos",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("prep_compounds()"),
                fileInput("cmp_file", "Tabla de compuestos (CSV)", accept = ".csv"),
                radioButtons("cmp_identifier", "Identificador", c("PubChem CID" = "pubchem", "SMILES" = "smiles")),
                textInput("cmp_id_col", "Columna identificadora (vacío = default)", value = ""),
                textInput("cmp_name_col", "Columna de nombre (vacío = default)", value = ""),
                checkboxInput("cmp_dedup", "Eliminar duplicados (por SMILES canónico)", TRUE),
                selectInput("cmp_on_missing_smiles", "Si falta SMILES", c("abort", "fetch", "drop")),
                selectInput("cmp_fetch_mode", "Modo de búsqueda externa", c("warn_and_cache", "abort")),
                actionButton("btn_prep_compounds", "Importar compuestos", class = "btn-primary")
              )
            ),
            column(8, h4("compounds(proj)"), render_table_ui("cmp_table"))
          )
        ),

        ## 3. Matriz / binarizacion ----------------------------------------
        tabPanel("3. Matriz",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("prep_binarize()"),
                fileInput("mat_file", "Matriz de abundancia (CSV)", accept = ".csv"),
                textInput("mat_id_col", "Columna de identificador de compuesto", value = "Name"),
                checkboxInput("mat_average_replicates", "Promediar réplicas (columnas repetidas)", TRUE),
                numericInput("mat_q", "Cuantil de corte (binarización)", value = 0.25, min = 0.01, max = 0.99, step = 0.01),
                actionButton("btn_prep_binarize", "Procesar matriz", class = "btn-primary")
              )
            ),
            column(8,
              h4("binarizedMatrix(proj)"), render_table_ui("mat_table")
            )
          )
        ),

        ## 4. RefDB ----------------------------------------------------------
        tabPanel("4. Base de referencia",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("refdb_build()"),
                checkboxGroupInput("refdb_sources", "Fuentes", c("pubchem", "chembl", "coconut"), selected = "pubchem"),
                selectInput("refdb_fetch_mode", "Modo de búsqueda externa", c("warn_and_cache", "abort")),
                actionButton("btn_refdb_build", "Construir base de referencia", class = "btn-primary")
              )
            ),
            column(8, h4("Resultado"), verbatimTextOutput("refdb_status"))
          )
        ),

        ## 5. ADME -------------------------------------------------------
        tabPanel("5. ADME",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("adme_local()"),
                checkboxGroupInput("adme_routes", "Vias a evaluar",
                  c("oral", "topical", "ophthalmic", "injectable"),
                  selected = c("oral", "topical", "ophthalmic", "injectable")),
                actionButton("btn_adme_local", "Calcular ADME local (rcdk)", class = "btn-primary")
              ),
              wellPanel(
                h4("adme_import()"),
                fileInput("adme_import_file", "Tabla ADME (CSV)", accept = ".csv"),
                selectInput("adme_import_platform", "Plataforma", c("swissadme", "admetlab", "other")),
                actionButton("btn_adme_import", "Importar ADME externo")
              ),
              wellPanel(
                h4("adme_filter()"),
                checkboxGroupInput("adme_rules", "Reglas",
                  c("ro5", "veber", "ghose", "egan", "oprea", "route"),
                  selected = c("ro5", "veber")),
                selectInput("adme_filter_source", "Fuente de datos", c("local", "imported")),
                checkboxInput("adme_hard_cutoff", "Descartar compuestos que fallan (hard cutoff)", FALSE),
                actionButton("btn_adme_filter", "Filtrar por ADME")
              )
            ),
            column(8, h4("Resultado ADME"), render_table_ui("adme_table"))
          )
        ),

        ## 6. Toxicidad -----------------------------------------------------
        tabPanel("6. Toxicidad",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("tox_local()"),
                checkboxGroupInput("tox_alert_sets", "Alertas estructurales", c("pains", "brenk"), selected = c("pains", "brenk")),
                actionButton("btn_tox_local", "Evaluar alertas estructurales", class = "btn-primary")
              ),
              wellPanel(
                h4("tox_import()"),
                fileInput("tox_import_file", "Tabla de toxicidad (CSV)", accept = ".csv"),
                selectInput("tox_import_platform", "Plataforma", c("admetlab", "swissadme", "other")),
                actionButton("btn_tox_import", "Importar toxicidad externa")
              ),
              wellPanel(
                h4("tox_safetyome() / tox_report()"),
                actionButton("btn_tox_safetyome", "Anotar safetyome (blancos)"),
                br(), br(),
                actionButton("btn_tox_report", "Generar reporte de toxicidad", class = "btn-primary")
              )
            ),
            column(8, h4("Resultado"), render_table_ui("tox_table"))
          )
        ),

        ## 7. Blancos --------------------------------------------------------
        tabPanel("7. Blancos",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("targets_import()"),
                fileInput("tgt_file", "Predicción de blancos (CSV)", accept = ".csv"),
                selectInput("tgt_platform", "Plataforma", c("swisstargetprediction", "superpred", "other")),
                selectInput("tgt_id_from", "Compuesto identificado por", c("filename", "column")),
                actionButton("btn_targets_import", "Importar blancos (1 archivo)", class = "btn-primary")
              ),
              wellPanel(
                h4("targets_import_batch()"),
                textInput("tgt_batch_dir", "Carpeta con varios archivos (ruta en el servidor)", value = ""),
                selectInput("tgt_batch_platform", "Plataforma", c("swisstargetprediction", "superpred", "other")),
                actionButton("btn_targets_import_batch", "Importar blancos (carpeta)")
              ),
              wellPanel(
                h4("targets_disease_filter()"),
                textInput("tgt_disease", "Enfermedad / condición (Open Targets)", value = ""),
                numericInput("tgt_min_score", "Score mínimo (vacío/NA = sin filtro)", value = NA, min = 0, max = 1, step = 0.05),
                actionButton("btn_targets_disease", "Filtrar blancos por enfermedad")
              )
            ),
            column(8, h4("Resultado"), render_table_ui("tgt_table"))
          )
        ),

        ## 8. Red --------------------------------------------------------
        tabPanel("8. Red",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("network_build()"),
                helpText("Construye una red bipartita compuesto-blanco por condición, a partir de binarizedMatrix() + blancos importados."),
                selectizeInput("net_condition", "Condición(es) (vacío = todas)", choices = NULL, multiple = TRUE),
                numericInput("net_min_score", "Score mínimo de blanco (vacío/NA = sin filtro)", value = NA, min = 0, max = 1, step = 0.05),
                actionButton("btn_network_build", "Construir red", class = "btn-primary")
              )
            ),
            column(8, h4("network_edges"), render_table_ui("net_table"))
          )
        ),

        ## 9. Analisis de red -------------------------------------------------
        tabPanel("9. Análisis de red",
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h4("network_centrality()"),
                checkboxGroupInput("nc_measures", "Métricas", c("degree", "betweenness", "hub_score"), selected = c("degree", "betweenness", "hub_score")),
                actionButton("btn_network_centrality", "Calcular centralidad", class = "btn-primary")
              ),
              wellPanel(
                h4("network_module_robustness()"),
                numericInput("nmr_min_module_size", "Tamaño mínimo de módulo", value = 2, min = 2, step = 1),
                numericInput("nmr_seed", "Semilla (vacía = aleatoria)", value = NA, step = 1),
                actionButton("btn_network_module_robustness", "Evaluar robustez modular")
              ),
              wellPanel(
                h4("network_proximity() / network_synergy()"),
                textInput("nprox_disease", "Enfermedad (STRING + módulo de enfermedad)", value = ""),
                actionButton("btn_network_proximity", "Calcular proximidad"),
                br(), br(),
                actionButton("btn_network_synergy", "Calcular sinergia")
              ),
              wellPanel(
                h4("Otros"),
                actionButton("btn_network_hub_penalty", "network_hub_penalty()"),
                br(), br(),
                actionButton("btn_network_degeneracy", "network_degeneracy()"),
                br(), br(),
                actionButton("btn_network_bowtie", "network_bowtie()"),
                br(), br(),
                actionButton("btn_network_motifs", "network_motifs()")
              )
            ),
            column(8, h4("Resultado"), render_table_ui("analysis_table"))
          )
        ),

        ## 10. Graficos --------------------------------------------------
        tabPanel("10. Gráficos",
          br(),
          sidebarLayout(
            sidebarPanel(
              width = 3,
              selectInput("plot_choice", "Gráfico", choices = c(
                "Centralidad" = "centrality",
                "Robustez modular" = "robustness",
                "Proximidad" = "proximity",
                "Sinergia" = "synergy",
                "Degeneracy" = "degeneracy",
                "Bowtie" = "bowtie",
                "Capas de red" = "layers",
                "BOILED-Egg (ADME)" = "boiled_egg",
                "Radar ADMET" = "admet_radar",
                "Espacio químico" = "chemical_space",
                "UpSet ADME" = "adme_upset"
              )),
              actionButton("btn_plot_render", "Generar", class = "btn-primary")
            ),
            mainPanel(
              width = 9,
              plotOutput("plot_out", height = "600px")
            )
          )
        ),

        ## 11. Exportar ------------------------------------------------
        tabPanel("11. Log / Exportar",
          br(),
          h4("Registro completo del proyecto (projectLog)"),
          render_table_ui("full_log_table"),
          tags$hr(),
          h4("Entradas disponibles en patliRResults(proj)"),
          verbatimTextOutput("results_names")
        )
      )
    )
  )
)

## ---- server ------------------------------------------------------------

server <- function(input, output, session) {
  rv <- reactiveValues(proj = NULL)

  output$project_status <- renderText({
    if (is.null(rv$proj)) return("Ningún proyecto cargado todavía.")
    paste0(
      "Carpeta: ", patliR::projectDir(rv$proj), "\n",
      "Compuestos: ", nrow(patliR::compounds(rv$proj)), "\n",
      "Condiciones: ", ncol(patliR::binarizedMatrix(rv$proj))
    )
  })

  render_table_server(output, "log_tail", function() {
    req(rv$proj)
    log_df <- patliR::projectLog(rv$proj)
    utils::tail(log_df, 15)
  })

  output$dl_log <- downloadHandler(
    filename = function() "patliR_log.csv",
    content = function(file) {
      req(rv$proj)
      utils::write.csv(patliR::projectLog(rv$proj), file, row.names = FALSE)
    }
  )

  observe({
    updateSelectizeInput(session, "net_condition", choices = condition_choices(rv$proj), server = TRUE)
  })

  ## 1. Proyecto ------------------------------------------------------------
  observeEvent(input$btn_new_project, {
    req(nzchar(input$new_project_dir))
    cache_dir <- if (nzchar(input$new_cache_dir)) input$new_cache_dir else NULL
    run_step(rv, function() patliR::patliR_project(input$new_project_dir, cache_dir = cache_dir),
              "Proyecto creado.")
  })
  observeEvent(input$btn_load_project, {
    req(nzchar(input$load_project_dir))
    cache_dir <- if (nzchar(input$load_cache_dir)) input$load_cache_dir else NULL
    run_step(rv, function() patliR::patliR_load(input$load_project_dir, cache_dir = cache_dir),
              "Proyecto cargado.")
  })

  ## 2. Compuestos ------------------------------------------------------------
  observeEvent(input$btn_prep_compounds, {
    req(rv$proj, input$cmp_file)
    data <- utils::read.csv(input$cmp_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    id_col <- if (nzchar(input$cmp_id_col)) input$cmp_id_col else NULL
    name_col <- if (nzchar(input$cmp_name_col)) input$cmp_name_col else NULL
    run_step(rv, function() {
      patliR::prep_compounds(
        rv$proj, data,
        identifier = input$cmp_identifier, id_col = id_col, name_col = name_col,
        dedup = input$cmp_dedup, on_missing_smiles = input$cmp_on_missing_smiles,
        fetch_mode = input$cmp_fetch_mode
      )
    }, "Compuestos importados.")
  })
  render_table_server(output, "cmp_table", function() { req(rv$proj); patliR::compounds(rv$proj) })

  ## 3. Matriz ------------------------------------------------------------
  observeEvent(input$btn_prep_binarize, {
    req(rv$proj, input$mat_file)
    data <- utils::read.csv(input$mat_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    run_step(rv, function() {
      patliR::prep_binarize(
        rv$proj, data, id_col = input$mat_id_col,
        average_replicates = input$mat_average_replicates, q = input$mat_q
      )
    }, "Matriz procesada y binarizada.")
  })
  render_table_server(output, "mat_table", function() { req(rv$proj); patliR::binarizedMatrix(rv$proj) })

  ## 4. RefDB ------------------------------------------------------------
  observeEvent(input$btn_refdb_build, {
    req(rv$proj, length(input$refdb_sources) > 0)
    run_step(rv, function() {
      patliR::refdb_build(rv$proj, sources = input$refdb_sources, fetch_mode = input$refdb_fetch_mode)
    }, "Base de referencia construida.")
  })
  output$refdb_status <- renderPrint({
    req(rv$proj)
    utils::tail(patliR::projectLog(rv$proj)[patliR::projectLog(rv$proj)$step == "refdb_build", ], 20)
  })

  ## 5. ADME ------------------------------------------------------------
  observeEvent(input$btn_adme_local, {
    req(rv$proj)
    run_step(rv, function() patliR::adme_local(rv$proj, routes = input$adme_routes), "ADME local calculado.")
  })
  observeEvent(input$btn_adme_import, {
    req(rv$proj, input$adme_import_file)
    run_step(rv, function() {
      patliR::adme_import(rv$proj, input$adme_import_file$datapath, platform = input$adme_import_platform)
    }, "ADME importado.")
  })
  observeEvent(input$btn_adme_filter, {
    req(rv$proj, length(input$adme_rules) > 0)
    run_step(rv, function() {
      patliR::adme_filter(
        rv$proj, rules = input$adme_rules, source = input$adme_filter_source,
        hard_cutoff = input$adme_hard_cutoff, ask = FALSE
      )
    }, "Filtro ADME aplicado.")
  })
  render_table_server(output, "adme_table", function() {
    req(rv$proj)
    res <- patliR::patliRResults(rv$proj, "adme_local")
    if (is.null(res)) res <- patliR::patliRResults(rv$proj, "adme_imported")
    if (is.null(res)) res <- patliR::patliRResults(rv$proj, "adme_filtered")
    req(res)
    res
  })

  ## 6. Toxicidad ------------------------------------------------------------
  observeEvent(input$btn_tox_local, {
    req(rv$proj, length(input$tox_alert_sets) > 0)
    run_step(rv, function() patliR::tox_local(rv$proj, alert_sets = input$tox_alert_sets), "Alertas estructurales evaluadas.")
  })
  observeEvent(input$btn_tox_import, {
    req(rv$proj, input$tox_import_file)
    run_step(rv, function() {
      patliR::tox_import(rv$proj, input$tox_import_file$datapath, platform = input$tox_import_platform)
    }, "Toxicidad importada.")
  })
  observeEvent(input$btn_tox_safetyome, {
    req(rv$proj)
    run_step(rv, function() patliR::tox_safetyome(rv$proj), "Safetyome anotado.")
  })
  observeEvent(input$btn_tox_report, {
    req(rv$proj)
    run_step(rv, function() patliR::tox_report(rv$proj), "Reporte de toxicidad generado.")
  })
  render_table_server(output, "tox_table", function() {
    req(rv$proj)
    res <- patliR::patliRResults(rv$proj, "tox_report")
    if (is.null(res)) res <- patliR::patliRResults(rv$proj, "tox_local")
    req(res)
    res
  })

  ## 7. Blancos ------------------------------------------------------------
  observeEvent(input$btn_targets_import, {
    req(rv$proj, input$tgt_file)
    run_step(rv, function() {
      patliR::targets_import(
        rv$proj, input$tgt_file$datapath, platform = input$tgt_platform, id_from = input$tgt_id_from
      )
    }, "Blancos importados.")
  })
  observeEvent(input$btn_targets_import_batch, {
    req(rv$proj, nzchar(input$tgt_batch_dir))
    run_step(rv, function() {
      patliR::targets_import_batch(rv$proj, input$tgt_batch_dir, platform = input$tgt_batch_platform)
    }, "Blancos importados (carpeta).")
  })
  observeEvent(input$btn_targets_disease, {
    req(rv$proj, nzchar(input$tgt_disease))
    min_score <- if (is.na(input$tgt_min_score)) NULL else input$tgt_min_score
    run_step(rv, function() {
      patliR::targets_disease_filter(rv$proj, disease = input$tgt_disease, min_score = min_score)
    }, "Blancos filtrados por enfermedad.")
  })
  render_table_server(output, "tgt_table", function() {
    req(rv$proj)
    res <- patliR::patliRResults(rv$proj, "targets_imported")
    req(res)
    res
  })

  ## 8. Red ------------------------------------------------------------
  observeEvent(input$btn_network_build, {
    req(rv$proj)
    condition <- if (length(input$net_condition) == 0) NULL else input$net_condition
    min_score <- if (is.na(input$net_min_score)) NULL else input$net_min_score
    run_step(rv, function() {
      patliR::network_build(rv$proj, condition = condition, min_score = min_score)
    }, "Red construida.")
  })
  render_table_server(output, "net_table", function() {
    req(rv$proj)
    res <- patliR::patliRResults(rv$proj, "network_edges")
    req(res)
    res
  })

  ## 9. Analisis de red ------------------------------------------------------
  observeEvent(input$btn_network_centrality, {
    req(rv$proj, length(input$nc_measures) > 0)
    run_step(rv, function() patliR::network_centrality(rv$proj, measures = input$nc_measures), "Centralidad calculada.")
  })
  observeEvent(input$btn_network_module_robustness, {
    req(rv$proj)
    seed <- if (is.na(input$nmr_seed)) NULL else input$nmr_seed
    run_step(rv, function() {
      patliR::network_module_robustness(rv$proj, min_module_size = input$nmr_min_module_size, seed = seed)
    }, "Robustez modular calculada.")
  })
  observeEvent(input$btn_network_proximity, {
    req(rv$proj, nzchar(input$nprox_disease))
    run_step(rv, function() patliR::network_proximity(rv$proj, disease = input$nprox_disease), "Proximidad calculada.")
  })
  observeEvent(input$btn_network_synergy, {
    req(rv$proj, nzchar(input$nprox_disease))
    run_step(rv, function() patliR::network_synergy(rv$proj, disease = input$nprox_disease), "Sinergia calculada.")
  })
  observeEvent(input$btn_network_hub_penalty, {
    req(rv$proj)
    run_step(rv, function() patliR::network_hub_penalty(rv$proj), "Penalización de hubs calculada.")
  })
  observeEvent(input$btn_network_degeneracy, {
    req(rv$proj)
    run_step(rv, function() patliR::network_degeneracy(rv$proj), "Degeneracy calculada.")
  })
  observeEvent(input$btn_network_bowtie, {
    req(rv$proj)
    run_step(rv, function() patliR::network_bowtie(rv$proj), "Clasificación bowtie calculada.")
  })
  observeEvent(input$btn_network_motifs, {
    req(rv$proj)
    run_step(rv, function() patliR::network_motifs(rv$proj), "Motivos de red calculados.")
  })

  analysis_result_name <- reactiveVal("network_centrality")
  observeEvent(input$btn_network_centrality, analysis_result_name("network_centrality"))
  observeEvent(input$btn_network_module_robustness, analysis_result_name("network_module_robustness"))
  observeEvent(input$btn_network_proximity, analysis_result_name("network_proximity"))
  observeEvent(input$btn_network_synergy, analysis_result_name("network_synergy"))
  observeEvent(input$btn_network_hub_penalty, analysis_result_name("network_hub_penalty"))
  observeEvent(input$btn_network_degeneracy, analysis_result_name("network_degeneracy"))
  observeEvent(input$btn_network_bowtie, analysis_result_name("network_bowtie"))
  observeEvent(input$btn_network_motifs, analysis_result_name("network_motifs"))

  render_table_server(output, "analysis_table", function() {
    req(rv$proj)
    res <- patliR::patliRResults(rv$proj, analysis_result_name())
    req(res)
    res
  })

  ## 10. Graficos ------------------------------------------------------------
  plot_obj <- eventReactive(input$btn_plot_render, {
    req(rv$proj)
    switch(input$plot_choice,
      centrality      = patliR::plot_centrality(rv$proj, save = FALSE),
      robustness      = patliR::plot_robustness(rv$proj, save = FALSE),
      proximity       = patliR::plot_proximity(rv$proj, save = FALSE),
      synergy         = patliR::plot_synergy(rv$proj, save = FALSE),
      degeneracy      = patliR::plot_network_degeneracy(rv$proj, save = FALSE),
      bowtie          = patliR::plot_bowtie(rv$proj, save = FALSE),
      layers          = patliR::plot_network_layers(rv$proj, save = FALSE),
      boiled_egg      = patliR::plot_boiled_egg(rv$proj, save = FALSE),
      admet_radar     = patliR::plot_admet_radar(rv$proj, engine = "static", save = FALSE),
      chemical_space  = patliR::plot_chemical_space(rv$proj, save = FALSE),
      adme_upset      = patliR::plot_adme_upset(rv$proj, save = FALSE)
    )
  })
  output$plot_out <- renderPlot({
    tryCatch(print(plot_obj()), error = function(e) { .notify_err(e); NULL })
  })

  ## 11. Exportar ------------------------------------------------------------
  render_table_server(output, "full_log_table", function() { req(rv$proj); patliR::projectLog(rv$proj) })
  output$results_names <- renderPrint({
    req(rv$proj)
    names(patliR::patliRResults(rv$proj))
  })
}

shinyApp(ui, server)
