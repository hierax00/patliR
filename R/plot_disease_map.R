#' @include AllGenerics.R internal.R network_build.R plot-helpers.R plot_network_layers.R
NULL

## Compound -> target -> disease-area map. The question it answers is "which
## disease areas does this extract touch, and through which proteins": the
## hundreds of Open Targets disease terms of targets_disease_profile() are
## grouped into a dozen therapeutic areas (transparent keyword rules, no
## ontology needed offline), and the custom gene sets of disease_genes_fetch()
## (e.g. hypertension, inflammation) ride along as highlighted "focus"
## diseases. Same visual language as plot_network_layers(view = "compound"):
## white background, viridis for the shared-target scale, a colour-blind
## friendly qualitative palette for the disease areas, clean labels.
##
## Every view is split into a pure data step (tested without pixels) and a
## drawing step.

#' Compound-target-disease map: which disease areas an extract's predicted
#' targets are associated with
#'
#' @description
#' Classifies every disease associated with the extract's predicted target
#' proteins into a therapeutic area (cardiovascular, neoplasm, nervous
#' system, ...; see [disease_map_classes()]) and draws the compound ->
#' protein -> disease chain in one of three views:
#'
#' * `"class"` (default) -- three columns: compounds (diamonds, left), the
#'   most shared target proteins (circles, middle; size and viridis colour =
#'   number of compounds hitting the target, as in
#'   [plot_network_layers()]`(view = "compound")`) and the disease areas
#'   (right; size = number of the extract's targets associated with at least
#'   one disease of the area, over **all** targets in scope, not only the
#'   drawn ones). Edges target -> area are coloured by the area and their
#'   width is the association score (the best score of the target among the
#'   area's diseases). Focus diseases (see below) are drawn as squares right
#'   under their area, with dashed edges.
#' * `"disease"` -- the same three columns, but the right column holds the
#'   `top_n_diseases` individual diseases associated with the most targets
#'   (ties: higher mean association score), plus every focus disease,
#'   grouped and coloured by area.
#' * `"heatmap"` -- a tile map of compounds (or, with `heatmap_rows =
#'   "target"`, the drawn target proteins) x disease areas, plus one column
#'   per focus disease; the tile colour and number are the count of the
#'   compound's targets associated with the area (for target rows: the
#'   number of the area's diseases the target is associated with).
#'
#' @section Data used:
#' Target-disease associations are pooled from up to three results tables
#' (`sources`), restricted to the targets of `network_edges` in scope:
#' * `"profile"` -- [targets_disease_profile()] (Open Targets, each target's
#'   top diseases with their resolved names; the main source).
#' * `"disease_genes"` -- [disease_genes_fetch()]'s custom disease gene sets
#'   (e.g. hypertension `MONDO_0005044`, inflammation `GO_0006954`); a target
#'   in the set counts as associated (GO gene sets carry no score).
#' * `"targets_disease"` -- [targets_disease_filter()]'s per-target score for
#'   one chosen disease.
#' Rows with an association score below `min_score` are dropped (rows
#' without a score are kept); a target-disease pair found in several sources
#' is counted once, with its highest score. Disease names come from
#' `targets_disease_profile` or `disease_genes` (the raw ID when neither
#' names it).
#'
#' @section Focus diseases:
#' The diseases of `disease_genes` and `targets_disease` are the ones the
#' analysis was designed around, so by default (`focus = NULL`) they are the
#' "focus" diseases: always shown in the `"class"`/`"disease"` views
#' (whatever their rank), drawn as squares with a bold label, and given their
#' own columns in the heatmap. They also count towards their area (the
#' hypertension gene set is part of "Cardiovascular"). Pass a character
#' vector of `disease_id`s to choose them, or `character(0)` for none.
#'
#' @section Disease areas are keyword rules, on purpose:
#' No disease ontology is available offline, so a disease is assigned to the
#' first area (in list order) with a pattern matching its name
#' (case-insensitive Perl regular expressions); unmatched names go to
#' `"Other"`. The rules are plain data -- print [disease_map_classes()] to
#' read them, and pass a modified copy as `disease_classes`, e.g.
#' `modifyList(disease_map_classes(), list(Skin = c("skin", "derma")))`.
#' The classification table returned as `attr(result, "classification")`
#' lists every target-disease pair with its area, so any assignment can be
#' checked.
#'
#' @section Which targets are drawn:
#' Only proteins associated with at least one disease shown on the right are
#' candidates. With `min_shared = NULL` the `max_targets` candidates hit by
#' the most compounds are drawn (ties: more disease links, then higher mean
#' prediction probability); with `min_shared = k` every candidate hit by at
#' least `k` compounds is. Compounds none of whose targets is drawn are
#' listed at the bottom of the left column, without edges. The counts on the
#' right (targets, diseases) always refer to every target in scope.
#'
#' @param proj A `PatliRProject`, with [network_build()] and at least one of
#'   [targets_disease_profile()], [disease_genes_fetch()] or
#'   [targets_disease_filter()] already run.
#' @param condition `NULL` (default, every built condition pooled) or a
#'   character vector of condition names.
#' @param view One of `"class"` (default), `"disease"` or `"heatmap"` -- see
#'   the description.
#' @param top_n_diseases Integer >= 1, default `15`: number of individual
#'   diseases in the `"disease"` view (focus diseases come on top).
#' @param max_targets Integer >= 1, default `30`: number of target proteins
#'   drawn in the middle column (and heatmap rows for `heatmap_rows =
#'   "target"`) when `min_shared` is `NULL`.
#' @param min_shared `NULL` (default) or an integer >= 1: draw every
#'   candidate target hit by at least this many compounds instead.
#' @param top_hub_n Integer >= 0, default `30`: number of drawn targets that
#'   get a gene-symbol label (the most shared ones).
#' @param disease_classes `NULL` (default, [disease_map_classes()]) or a
#'   named list of character vectors of regular expressions -- see the
#'   section on disease areas.
#' @param focus `NULL` (default, automatic), a character vector of
#'   `disease_id`s, or `character(0)` -- see the section on focus diseases.
#' @param sources Subset of `c("profile", "disease_genes",
#'   "targets_disease")`: which association tables to use (all by default;
#'   missing tables are skipped).
#' @param min_score Number in `[0, 1]`, default `0.1`: minimum Open Targets
#'   association score of a target-disease pair (rows without a score are
#'   kept). The single-disease scores of `targets_disease` reach down to
#'   ~0.001, which would otherwise link almost every target to the focus
#'   disease.
#' @param heatmap_rows `"compound"` (default) or `"target"`: rows of the
#'   `"heatmap"` view.
#' @param save Logical, default `TRUE`: also write a PNG to `out_dir` and
#'   log it.
#' @param out_dir Directory for the PNG. Defaults to
#'   `file.path(projectDir(proj), "plots")`.
#' @param width,height,dpi Passed to [ggplot2::ggsave()].
#'
#' @return A `ggplot` object with attributes `"classification"`
#'   (`data.frame(target_id, target_label, disease_id, disease_name,
#'   disease_class, association_score, source, focus)`, one row per
#'   target-disease pair in scope), `"disease_summary"` (one row per
#'   disease: `disease_id, disease_name, disease_class, n_targets,
#'   mean_score, focus, rank`) and `"class_summary"` (one row per area:
#'   `disease_class, n_targets, n_diseases, mean_score, n_focus`). When
#'   `save = TRUE` the PNG (`disease_map_<condition>_<view>_<subset>.png`) is
#'   logged in `patliRResults(proj, "disease_map_plot_log")` (one row per
#'   `(condition, view, subset)`, upserted) and the updated project is
#'   `attr(result, "proj")`.
#'
#' @examples
#' \dontrun{
#' proj <- patliR_load("my_project") # network_build() + targets_disease_profile() run
#' p <- plot_disease_map(proj, condition = "EFLO-S")
#' plot_disease_map(proj, condition = "EFLO-S", view = "disease", top_n_diseases = 20)
#' plot_disease_map(proj, condition = "EFLO-S", view = "heatmap")
#' head(attr(p, "classification"))
#'
#' ## Adjust the keyword rules: move "stroke" from Cardiovascular to Nervous
#' cl <- disease_map_classes()
#' cl$Cardiovascular <- setdiff(cl$Cardiovascular, "stroke")
#' cl$`Nervous system & psychiatric` <- c(cl$`Nervous system & psychiatric`, "stroke")
#' plot_disease_map(proj, condition = "EFLO-S", disease_classes = cl)
#' }
#'
#' @seealso [disease_map_classes()], [targets_disease_profile()],
#'   [disease_genes_fetch()], [plot_network_layers()], [plot_disease_network()]
#' @export
plot_disease_map <- function(proj, condition = NULL, view = c("class", "disease", "heatmap"),
                             top_n_diseases = 15, max_targets = 30, min_shared = NULL, top_hub_n = 30,
                             disease_classes = NULL, focus = NULL,
                             sources = c("profile", "disease_genes", "targets_disease"),
                             min_score = 0.1, heatmap_rows = c("compound", "target"),
                             save = TRUE, out_dir = NULL, width = 11, height = 8.5, dpi = 150) {
  stopifnot(is(proj, "PatliRProject"))
  view <- match.arg(view)
  heatmap_rows <- match.arg(heatmap_rows)
  if (!.network_layers_is_count(top_n_diseases)) cli::cli_abort("{.arg top_n_diseases} must be a single integer >= 1.")
  if (!.network_layers_is_count(max_targets)) cli::cli_abort("{.arg max_targets} must be a single integer >= 1.")
  if (!is.null(min_shared) && !.network_layers_is_count(min_shared)) {
    cli::cli_abort("{.arg min_shared} must be {.code NULL} or a single integer >= 1.")
  }
  if (!is.numeric(top_hub_n) || length(top_hub_n) != 1L || !is.finite(top_hub_n) ||
      top_hub_n < 0 || top_hub_n != floor(top_hub_n)) {
    cli::cli_abort("{.arg top_hub_n} must be a non-negative integer.")
  }
  if (!is.numeric(min_score) || length(min_score) != 1L || is.na(min_score) || min_score < 0 || min_score > 1) {
    cli::cli_abort("{.arg min_score} must be a single number in [0, 1].")
  }
  if (!is.character(sources) || length(sources) == 0 ||
      length(setdiff(sources, c("profile", "disease_genes", "targets_disease"))) > 0) {
    cli::cli_abort("{.arg sources} must be a non-empty subset of {.val {c('profile', 'disease_genes', 'targets_disease')}}.")
  }
  if (!is.null(focus) && !is.character(focus)) {
    cli::cli_abort("{.arg focus} must be {.code NULL}, a character vector of disease IDs, or {.code character(0)}.")
  }
  classes <- .disease_map_check_classes(disease_classes)
  .plot_require()

  scope <- .plot_scope(proj, condition)
  conditions <- scope$conditions
  scope_label <- scope$scope_label

  ct <- .disease_map_ct(proj, conditions)
  if (nrow(ct) == 0) cli::cli_abort("No compound-target edge in condition(s) {.val {conditions}}.")
  assoc <- .disease_map_associations(proj, unique(ct$uniprot_id), sources = sources, min_score = min_score)
  if (nrow(assoc) == 0) {
    cli::cli_abort(c(
      "No target-disease association for the targets of condition(s) {.val {conditions}}.",
      "i" = "Run {.fn targets_disease_profile} (or {.fn disease_genes_fetch}) first, or lower {.arg min_score}."
    ))
  }
  assoc$disease_class <- .disease_map_classify(assoc$disease_name, classes)
  assoc$focus <- .disease_map_focus(proj, assoc, focus)
  sums <- .disease_map_summaries(assoc, names(classes))

  build <- switch(view,
    class = , disease = .disease_map_view_network(
      proj, ct, assoc, sums, conditions, view = view, top_n_diseases = top_n_diseases,
      max_targets = max_targets, min_shared = min_shared, top_hub_n = top_hub_n,
      palette = .disease_map_palette(names(classes)), title_suffix = scope_label, fig_width = width
    ),
    heatmap = .disease_map_view_heatmap(
      proj, ct, assoc, sums, conditions, rows = heatmap_rows, max_targets = max_targets,
      min_shared = min_shared, title_suffix = scope_label, fig_width = width
    )
  )
  p <- build$plot

  labels <- .plot_label_nodes(proj, conditions, unique(assoc$target_id), "target")
  classification <- data.frame(
    target_id = assoc$target_id,
    target_label = labels[match(assoc$target_id, unique(assoc$target_id))],
    disease_id = assoc$disease_id, disease_name = assoc$disease_name,
    disease_class = assoc$disease_class, association_score = assoc$association_score,
    source = assoc$source, focus = assoc$focus, stringsAsFactors = FALSE
  )
  classification <- classification[order(classification$target_id, classification$disease_class,
                                         classification$disease_name), , drop = FALSE]
  rownames(classification) <- NULL
  attr(p, "classification") <- classification
  attr(p, "disease_summary") <- sums$diseases
  attr(p, "class_summary") <- sums$classes

  .plot_finish(
    proj, p,
    name = "disease_map_plot_log",
    filename = paste0("disease_map_", scope_label, "_", view, if (nzchar(build$tag)) paste0("_", build$tag), ".png"),
    log_row = data.frame(
      condition = scope_label, view = view, subset = build$tag, path = NA_character_,
      n_compounds = build$n_compounds, n_targets = build$n_targets, n_right = build$n_right,
      n_edges = build$n_edges, n_diseases = nrow(sums$diseases),
      n_classes = sum(sums$classes$n_targets > 0), stringsAsFactors = FALSE
    ),
    key_cols = c("condition", "view", "subset"),
    save = save, out_dir = out_dir, width = width, height = height, dpi = dpi
  )
}

#' Default keyword rules grouping disease names into therapeutic areas
#'
#' @description
#' The rules [plot_disease_map()] uses to put each disease into one of
#' twelve therapeutic areas plus `"Other"`: a named list, one element per
#' area, each a character vector of case-insensitive Perl regular
#' expressions matched against the disease **name**. Areas are tried in list
#' order and the first match wins, so the order resolves overlaps: a lung
#' carcinoma is a neoplasm (not respiratory), a viral pneumonia is
#' infectious, "inflammatory bowel disease" is immune/inflammatory, the
#' hypertension gene set is cardiovascular and the inflammation GO term is
#' immune/inflammatory. Names matching nothing (e.g. traits such as
#' "smoking initiation", or syndromes named after people) go to `"Other"`.
#'
#' The areas loosely follow the top level of the EFO / Open Targets
#' "therapeutic area" hierarchy, merged into a readable number of groups:
#' Neoplasm; Infectious; Immune & inflammatory; Hematologic; Eye & ear;
#' Cardiovascular; Nervous system & psychiatric; Metabolic & endocrine;
#' Respiratory; Digestive & liver; Musculoskeletal; Renal & urogenital (with
#' reproduction and pregnancy); Skin.
#'
#' @return A named list of character vectors (regular expressions).
#' @examples
#' cl <- disease_map_classes()
#' names(cl)
#' cl$Cardiovascular
#' @seealso [plot_disease_map()]
#' @export
disease_map_classes <- function() {
  list(
    Neoplasm = c(
      "cancer", "carcinom", "neoplas", "tumou?r", "leuka?emia", "lymphoma", "myeloma", "melanoma",
      "sarcoma", "gliom", "blastoma", "adenoma", "mesothelioma", "malignan", "metasta", "mastocytosis",
      "myelofibrosis", "thrombocythemia", "polycythemia vera", "myeloproliferative", "leukoplakia",
      "hemangioma", "pheochromocytoma", "fibromatosis", "leiomyoma", "fibroma", "lipoma"
    ),
    Infectious = c(
      "infect", "\\bvir(al|us)", "virus", "bacteri", "fung", "sepsis", "septic", "covid", "\\bhiv\\b",
      "hepatitis [a-e]\\b", "pneumonia", "dengue", "influenza", "orthomyxo", "leishmania", "listeriosis",
      "diphtheria", "tubercul", "malaria", "abscess", "cellulitis", "candidiasis", "herpes", "zoster",
      "spotted fever", "mycosis", "protozoa", "severe acute respiratory syndrome", "epstein-barr",
      "parasit", "helminth", "lyme"
    ),
    `Immune & inflammatory` = c(
      "(?<!non-)autoimmun", "autoinflamm", "inflamm", "immunodeficien", "(?<![a-z-])immune", "lupus",
      "rheumatoid", "psoriatic arthritis", "juvenile idiopathic arthritis", "spondylitis", "vasculitis",
      "arteritis", "polyangiitis", "sj(o|ö)gren", "sarcoidosis", "graft versus host", "crohn", "colitis",
      "celiac", "allerg", "atopy", "hypersensitivity", "anaphyla", "agammaglobulin", "hyper-ige",
      "behcet", "hemophagocytic", "lymphoproliferative", "whim syndrome", "mhc class", "myasthenia gravis",
      "eosinophil", "lymphaden", "fever"
    ),
    Hematologic = c(
      "\\ban(a)?emi", "bleeding", "ha?emorrha", "ha?emophilia", "platelet", "thrombocyt", "thrombasthenia",
      "thrombophilia", "coagulation", "factor [ivx]+\\b", "prothrombin", "protein c deficiency",
      "ha?emochromatosis", "erythrocytosis", "polycyth", "neutropeni", "bone marrow", "splenomegaly",
      "hydrocytosis", "\\bblood", "iron deficiency", "warfarin", "methemoglobin", "thalass", "sickle",
      "ha?emoly", "pancytopeni"
    ),
    `Eye & ear` = c(
      "\\beye", "ocular", "retin", "macula", "glaucom", "cataract", "myopi", "uveit", "iridocycl",
      "keratitis", "\\boptic", "vision", "ophthalm", "cone dystrophy", "cone-rod", "strabism", "refractive",
      "anisometropia", "ectropion", "entropion", "s-cone", "goldmann-favre", "cornea", "conjunctiv",
      "blind", "keratoconus", "hearing", "deaf", "tinnitus", "otitis", "\\bear\\b", "vestibul"
    ),
    Cardiovascular = c(
      "cardi", "heart", "hypertens", "hypotens", "coronary", "myocard", "arrhythm", "atrial", "ventric",
      "aort", "aneurysm", "atheroscl", "arterioscl", "stroke", "vascul", "\\bvein", "venous", "varicose",
      "thrombo", "embol", "infarct", "qt (interval|syndrome)", "brugada", "tachycard", "bradycard",
      "atrioventricular", "wolff-parkinson", "telangiectasia", "arter(y|ial)", "angina", "raynaud",
      "loeys-dietz", "marfan"
    ),
    `Nervous system & psychiatric` = c(
      "neur", "nerv", "brain", "cerebr", "cerebell", "enceph", "alzheimer", "parkinson", "dementia",
      "epilep", "seizure", "ataxia", "paraplegia", "spastic", "dyston", "tremor", "migraine", "headache",
      "^pain$", "neuropathic pain", "pain (disorder|syndrome)", "insensitivity to pain", "erythermalgia",
      "intellectual", "developmental (delay|disability)", "learning disability", "autis", "attention deficit",
      "schizo", "psychot", "psychiat", "depress", "bipolar", "anxiety", "panic", "obsessive", "\\bmental",
      "stress disorder", "tourette", "insomnia", "narcolep", "cataplexy", "restless legs", "dependence",
      "substance", "opioid use", "addict", "melancholia", "behavio", "hemiplegia", "myoclon", "startle",
      "hyperekplexia", "lissenceph", "microcephal", "macroceph", "charcot-marie", "myasthen",
      "lambert-eaton", "isaacs", "kennedy disease", "friedreich", "lewy", "\\bcortical", "lipofuscinosis",
      "movement disorder", "motor", "dravet", "lennox", "sedation", "cognit", "huntington", "amyotrophic",
      "multiple sclerosis", "paralysis"
    ),
    `Metabolic & endocrine` = c(
      "diabet", "obes", "insulin", "glyc", "lipid", "lipodystrophy", "lipoprotein", "cholesterol", "cholan",
      "metabol", "thyroid", "goit(er|re)", "adrenal", "adreno", "cortisol", "cortisone", "aldosteron",
      "hypogonad", "growth hormone", "androgen", "aromatase", "gout", "hyperuric", "storage disease",
      "lysosomal", "gaucher", "niemann-pick", "mucopolysacchar", "mitochondri", "amyloidosis", "carnitine",
      "hypophosphat", "rickets", "vitamin", "overnutrition", "polyphagia", "aminoacid", "amino acid",
      "hypermethionin", "hyperoxaluria", "xanthinuria", "tryptophan", "trehalase", "fructose", "galactos",
      "glutamin", "cushing", "mineralocorticoid", "glucocorticoid", "short stature", "graves", "endocrin",
      "parathyroid", "pituitar", "acromegaly", "prader-willi", "alstrom"
    ),
    Respiratory = c(
      "lung", "pulmon", "respirat", "asthma", "bronch", "airway", "pneum", "rhinitis", "nasal", "sinusitis",
      "cough", "apnea", "apnoea", "cystic fibrosis", "emphysema", "kartagener", "pleur", "trache", "laryn"
    ),
    `Digestive & liver` = c(
      "liver", "hepat", "cirrhosis", "biliar", "cholang", "cholelith", "gallstone", "gallbladder", "pancrea",
      "gastr", "intestin", "bowel", "colon", "rect(al|um|o)", "esophag", "oesophag", "reflux", "ulcer",
      "digestive", "nausea", "vomiting", "diarrh", "constipation", "proctitis", "hirschsprung", "periodont",
      "tooth", "dental", "teeth", "\\boral\\b", "mucositis", "salivary", "malabsorption", "steato",
      "duoden", "trypsin", "enteropathy", "appendic"
    ),
    Musculoskeletal = c(
      "skelet", "\\bbone", "osteo", "arthr", "joint", "muscul", "muscle", "myopath", "myoton", "myosit",
      "dystrophy", "cartilage", "fracture", "sprain", "scoliosis", "spondyl", "\\bdis[ck]\\b", "back pain",
      "knee pain", "fibromyalgia", "tendin", "tendon", "fasciitis", "contracture", "carpal tunnel",
      "frozen shoulder", "genu (valgum|varum)", "\\blimbs?\\b", "\\bfoot\\b", "brachydactyly", "polydactyly",
      "dysostosis", "pterygium", "chondro", "ollier", "rheumatic", "cramp", "spasm", "temporomandibular",
      "tietze", "dupuytren", "dysplasia", "clubbing", "cleft"
    ),
    `Renal & urogenital` = c(
      "kidney", "renal", "nephr", "glomerul", "urin", "bladder", "ureth", "ureter", "urolith", "\\buro",
      "dysuria", "enuresis", "incontinence", "prostat", "testic", "testis", "ovar", "uter", "endometri",
      "cervicitis", "vas deferens", "fallopian", "infertility", "erectile", "sexual dysfunction",
      "dyspareunia", "peyronie", "menopaus", "menstru", "pregnan", "eclampsia", "placent", "premature birth",
      "preterm", "cesarean", "caesarean", "gestation", "fetal", "genital", "vagin", "vulv"
    ),
    Skin = c(
      "skin", "derma", "psoria", "eczema", "acne", "alopecia", "vitiligo", "kerato(derma|sis)", "keratolytic",
      "rosacea", "sebaceous", "seborrh", "hyperhidrosis", "erythema", "pigment", "piebald", "keloid",
      "cutaneous", "hair loss", "ichthyosis", "pemphig", "urticaria", "lichen", "wound", "prurit",
      "blister", "epiderm", "palmoplantar", "hidradenitis"
    )
  )
}

#' Validate / default the `disease_classes` argument
#' @return A named list of non-empty character vectors.
#' @keywords internal
.disease_map_check_classes <- function(disease_classes) {
  if (is.null(disease_classes)) return(disease_map_classes())
  ok <- is.list(disease_classes) && length(disease_classes) > 0 && !is.null(names(disease_classes)) &&
    all(nzchar(names(disease_classes))) && !anyDuplicated(names(disease_classes)) &&
    all(vapply(disease_classes, function(x) is.character(x) && length(x) > 0 && !anyNA(x), logical(1)))
  if (!ok) {
    cli::cli_abort(c(
      "{.arg disease_classes} must be a named list of non-empty character vectors (regular expressions), one per disease area.",
      "i" = "Start from {.code disease_map_classes()} and modify it."
    ))
  }
  bad <- character(0)
  for (nm in names(disease_classes)) {
    for (pat in disease_classes[[nm]]) {
      res <- tryCatch(suppressWarnings(grepl(pat, "x", perl = TRUE)), error = function(e) e)
      if (inherits(res, "error")) bad <- c(bad, pat)
    }
  }
  if (length(bad) > 0) cli::cli_abort("Invalid regular expression(s) in {.arg disease_classes}: {.val {bad}}.")
  disease_classes[names(disease_classes) != "Other"]
}

#' Assign each disease name to the first matching area
#'
#' @param disease_name Character vector.
#' @param classes Named list of regular-expression vectors, as
#'   [disease_map_classes()].
#' @return Character vector, same length as `disease_name`; `"Other"` when
#'   nothing matches (or the name is missing).
#' @keywords internal
.disease_map_classify <- function(disease_name, classes) {
  out <- rep("Other", length(disease_name))
  todo <- !is.na(disease_name) & nzchar(disease_name)
  for (nm in names(classes)) {
    if (!any(todo)) break
    pat <- paste0("(?:", paste(classes[[nm]], collapse = ")|(?:"), ")")
    hit <- todo
    hit[todo] <- grepl(pat, disease_name[todo], ignore.case = TRUE, perl = TRUE)
    out[hit] <- nm
    todo <- todo & !hit
  }
  out
}

#' Fixed colour per default area (so an area keeps its colour across views
#' and conditions); areas of a custom `disease_classes` beyond these get
#' `hcl.colors(, "Dark 3")`, `"Other"` is always grey
#' @keywords internal
.disease_map_palette <- function(class_names) {
  fixed <- c(
    Neoplasm = "#882255", Infectious = "#999933", `Immune & inflammatory` = "#EE7733",
    Hematologic = "#AA4499", `Eye & ear` = "#5FA8D3", Cardiovascular = "#CC3311",
    `Nervous system & psychiatric` = "#332288", `Metabolic & endocrine` = "#DDAA33",
    Respiratory = "#0077BB", `Digestive & liver` = "#117733", Musculoskeletal = "#8C564B",
    `Renal & urogenital` = "#44AA99", Skin = "#E7849B"
  )
  class_names <- setdiff(class_names, "Other")
  pal <- fixed[intersect(class_names, names(fixed))]
  extra <- setdiff(class_names, names(fixed))
  if (length(extra) > 0) pal <- c(pal, stats::setNames(grDevices::hcl.colors(length(extra), "Dark 3"), extra))
  c(pal[class_names], Other = "#A6A6A6")
}

#' Compound-target table of the conditions in scope (one row per pair,
#' maximum prediction probability across pooled conditions)
#' @return `data.frame(compound_id, uniprot_id, weight)`.
#' @keywords internal
.disease_map_ct <- function(proj, conditions) {
  e <- patliRResults(proj, "network_edges")
  e <- e[e$condition %in% conditions, , drop = FALSE]
  if (nrow(e) == 0) return(data.frame(compound_id = character(), uniprot_id = character(), weight = double()))
  if (!"weight" %in% names(e)) e$weight <- NA_real_
  key <- paste(e$compound_id, e$uniprot_id, sep = "\r")
  w <- tapply(e$weight, key, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  first <- !duplicated(key)
  ct <- data.frame(compound_id = e$compound_id[first], uniprot_id = e$uniprot_id[first],
                   weight = as.vector(w[key[first]]), stringsAsFactors = FALSE)
  ct <- ct[order(ct$compound_id, ct$uniprot_id), , drop = FALSE]
  rownames(ct) <- NULL
  ct
}

#' Pool the target-disease association tables for a set of targets
#'
#' @description
#' Reads `targets_disease_profile` (`"profile"`), `disease_genes` and
#' `targets_disease` (whichever are requested and exist), keeps the rows of
#' `targets`, drops rows scored below `min_score` (unscored rows are kept),
#' resolves disease names (profile first, then disease_genes, else the ID;
#' a trailing `" (GO:0006954)"` is dropped) and collapses a pair found in
#' several tables to one row with its highest score and the sources joined
#' by `"+"`.
#' @return `data.frame(target_id, disease_id, disease_name,
#'   association_score, source)`, one row per target-disease pair.
#' @keywords internal
.disease_map_associations <- function(proj, targets, sources = c("profile", "disease_genes", "targets_disease"),
                                      min_score = 0) {
  empty <- data.frame(target_id = character(), disease_id = character(), disease_name = character(),
                      association_score = double(), source = character(), stringsAsFactors = FALSE)
  parts <- list()
  name_tables <- list()
  prof <- if ("profile" %in% sources) patliRResults(proj, "targets_disease_profile") else NULL
  if (!is.null(prof) && nrow(prof) > 0) {
    name_tables$profile <- prof[, c("disease_id", "disease_name")]
    p <- prof[prof$target_id %in% targets, , drop = FALSE]
    parts$profile <- data.frame(target_id = p$target_id, disease_id = p$disease_id,
                                association_score = as.numeric(p$association_score), source = rep("profile", nrow(p)),
                                stringsAsFactors = FALSE)
  }
  dg <- if ("disease_genes" %in% sources) patliRResults(proj, "disease_genes") else NULL
  if (!is.null(dg) && nrow(dg) > 0) {
    if ("disease_name" %in% names(dg)) name_tables$dg <- dg[, c("disease_id", "disease_name")]
    d <- dg[dg$uniprot_id %in% targets, , drop = FALSE]
    sc <- if ("association_score" %in% names(d)) as.numeric(d$association_score) else rep(NA_real_, nrow(d))
    parts$dg <- data.frame(target_id = d$uniprot_id, disease_id = d$disease_id, association_score = sc,
                           source = rep("disease_genes", nrow(d)), stringsAsFactors = FALSE)
  }
  td <- if ("targets_disease" %in% sources) patliRResults(proj, "targets_disease") else NULL
  if (!is.null(td) && nrow(td) > 0) {
    d <- td[td$target_id %in% targets, , drop = FALSE]
    sc <- if ("association_score" %in% names(d)) as.numeric(d$association_score) else rep(NA_real_, nrow(d))
    parts$td <- data.frame(target_id = d$target_id, disease_id = d$disease_id, association_score = sc,
                           source = rep("targets_disease", nrow(d)), stringsAsFactors = FALSE)
  }
  if (length(parts) == 0) return(empty)
  a <- do.call(rbind, unname(parts))
  a <- a[!is.na(a$target_id) & !is.na(a$disease_id), , drop = FALSE]
  a <- a[is.na(a$association_score) | a$association_score >= min_score, , drop = FALSE]
  if (nrow(a) == 0) return(empty)

  key <- paste(a$target_id, a$disease_id, sep = "\r")
  score <- tapply(a$association_score, key, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  src <- tapply(a$source, key, function(x) paste(sort(unique(x)), collapse = "+"))
  first <- !duplicated(key)
  out <- data.frame(target_id = a$target_id[first], disease_id = a$disease_id[first],
                    association_score = as.vector(score[key[first]]), source = as.vector(src[key[first]]),
                    stringsAsFactors = FALSE)

  names_df <- do.call(rbind, c(unname(name_tables), list(data.frame(disease_id = character(), disease_name = character()))))
  names_df <- names_df[!is.na(names_df$disease_name) & nzchar(names_df$disease_name), , drop = FALSE]
  names_df <- names_df[!duplicated(names_df$disease_id), , drop = FALSE]
  nm <- names_df$disease_name[match(out$disease_id, names_df$disease_id)]
  nm <- ifelse(is.na(nm), out$disease_id, nm)
  out$disease_name <- sub("\\s*\\((GO|MONDO|EFO|HP|MP|Orphanet|DOID)[:_][0-9]+\\)$", "", nm)
  out <- out[order(out$target_id, out$disease_id), c("target_id", "disease_id", "disease_name", "association_score", "source")]
  rownames(out) <- NULL
  out
}

#' Which association rows belong to a focus disease
#' @param focus `NULL` (diseases of `disease_genes` / `targets_disease`), or
#'   a character vector of disease IDs.
#' @return Logical vector along `assoc`.
#' @keywords internal
.disease_map_focus <- function(proj, assoc, focus = NULL) {
  if (is.null(focus)) {
    ids <- character(0)
    for (nm in c("disease_genes", "targets_disease")) {
      tab <- patliRResults(proj, nm)
      if (!is.null(tab) && nrow(tab) > 0 && "disease_id" %in% names(tab)) ids <- c(ids, tab$disease_id)
    }
    focus <- unique(ids)
  }
  assoc$disease_id %in% focus
}

#' Per-disease and per-area summaries of an association table
#'
#' @description
#' Diseases are ranked by the number of distinct targets associated, then
#' the mean association score (unscored pairs ignored; a disease with no
#' score at all ranks after scored ones at equal count), then name. Areas:
#' distinct targets, distinct diseases, mean score, number of focus
#' diseases; every area of `class_names` plus `"Other"` gets a row (zero
#' when unused), sorted by `n_targets` with `"Other"` last.
#' @param assoc Association table with `disease_class` and `focus` columns.
#' @return `list(diseases, classes)`.
#' @keywords internal
.disease_map_summaries <- function(assoc, class_names) {
  mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  ids <- unique(assoc$disease_id)
  i <- match(ids, assoc$disease_id)
  dis <- data.frame(
    disease_id = ids, disease_name = assoc$disease_name[i], disease_class = assoc$disease_class[i],
    n_targets = as.integer(tapply(assoc$target_id, assoc$disease_id, function(x) length(unique(x)))[ids]),
    mean_score = as.vector(tapply(assoc$association_score, assoc$disease_id, mean_na)[ids]),
    focus = as.vector(tapply(assoc$focus, assoc$disease_id, any)[ids]),
    stringsAsFactors = FALSE
  )
  ms <- ifelse(is.na(dis$mean_score), -Inf, dis$mean_score)
  dis <- dis[order(-dis$n_targets, -ms, dis$disease_name, dis$disease_id), , drop = FALSE]
  dis$rank <- seq_len(nrow(dis))
  rownames(dis) <- NULL

  all_classes <- unique(c(setdiff(class_names, "Other"), "Other"))
  cls <- data.frame(
    disease_class = all_classes,
    n_targets = vapply(all_classes, function(k) length(unique(assoc$target_id[assoc$disease_class == k])), integer(1)),
    n_diseases = vapply(all_classes, function(k) length(unique(assoc$disease_id[assoc$disease_class == k])), integer(1)),
    mean_score = vapply(all_classes, function(k) mean_na(assoc$association_score[assoc$disease_class == k]), numeric(1)),
    n_focus = vapply(all_classes, function(k) length(unique(assoc$disease_id[assoc$disease_class == k & assoc$focus])), integer(1)),
    stringsAsFactors = FALSE
  )
  cls <- cls[order(cls$disease_class == "Other", -cls$n_targets, cls$disease_class), , drop = FALSE]
  rownames(cls) <- NULL
  list(diseases = dis, classes = cls)
}

#' Choose the target proteins drawn in the middle column
#'
#' @param ct Compound-target table (`compound_id, uniprot_id, weight`).
#' @param links `data.frame(target_id, node_id)` -- the target's links to
#'   the right-hand nodes on show; only linked targets are candidates.
#' @param max_targets,min_shared See [plot_disease_map()].
#' @return `data.frame(uniprot_id, n_compounds, n_links, mean_weight)` of
#'   the drawn targets, most shared first.
#' @keywords internal
.disease_map_select_targets <- function(ct, links, max_targets = 30, min_shared = NULL) {
  cand <- intersect(unique(ct$uniprot_id), unique(links$target_id))
  if (length(cand) == 0) {
    return(data.frame(uniprot_id = character(), n_compounds = integer(), n_links = integer(), mean_weight = double()))
  }
  sub <- ct[ct$uniprot_id %in% cand, , drop = FALSE]
  n_c <- tapply(sub$compound_id, sub$uniprot_id, function(x) length(unique(x)))
  mw <- tapply(sub$weight, sub$uniprot_id, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  lk <- links[links$target_id %in% cand, , drop = FALSE]
  n_l <- tapply(lk$node_id, lk$target_id, function(x) length(unique(x)))
  tg <- data.frame(uniprot_id = cand, n_compounds = as.integer(n_c[cand]), n_links = as.integer(n_l[cand]),
                   mean_weight = as.vector(mw[cand]), stringsAsFactors = FALSE)
  mws <- ifelse(is.na(tg$mean_weight), -Inf, tg$mean_weight)
  tg <- tg[order(-tg$n_compounds, -tg$n_links, -mws, tg$uniprot_id), , drop = FALSE]
  tg <- if (is.null(min_shared)) utils::head(tg, max_targets) else tg[tg$n_compounds >= min_shared, , drop = FALSE]
  rownames(tg) <- NULL
  tg
}

#' Right-hand nodes and target links of the class / disease views
#'
#' @description
#' `"class"`: one node per area with at least one associated target (sorted
#' by number of targets, `"Other"` last), each followed by that area's focus
#' diseases. `"disease"`: the `top_n` best-ranked diseases plus every focus
#' disease, grouped by area (areas in the order of the class summary),
#' focus diseases first within an area, then by rank. Link score = best
#' association score of the target among the node's diseases.
#' @return `list(nodes = data.frame(node_id, label, disease_class, type
#'   ("class"/"disease"/"focus"), n_targets, n_diseases, mean_score, group),
#'   links = data.frame(target_id, node_id, score, focus))`.
#' @keywords internal
.disease_map_right_nodes <- function(assoc, sums, view = c("class", "disease"), top_n = 15) {
  view <- match.arg(view)
  max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  dis <- sums$diseases
  foc <- dis[dis$focus, , drop = FALSE]
  class_order <- sums$classes$disease_class[sums$classes$n_targets > 0]
  focus_nodes <- if (nrow(foc) > 0) {
    data.frame(node_id = foc$disease_id, label = foc$disease_name, disease_class = foc$disease_class,
               type = "focus", n_targets = foc$n_targets, n_diseases = 1L, mean_score = foc$mean_score,
               stringsAsFactors = FALSE)
  } else NULL

  if (view == "class") {
    cls <- sums$classes[sums$classes$n_targets > 0, , drop = FALSE]
    nodes <- data.frame(node_id = cls$disease_class, label = cls$disease_class, disease_class = cls$disease_class,
                        type = "class", n_targets = cls$n_targets, n_diseases = cls$n_diseases,
                        mean_score = cls$mean_score, stringsAsFactors = FALSE)
    nodes <- rbind(nodes, focus_nodes)
    nodes <- nodes[order(match(nodes$disease_class, class_order), nodes$type != "class",
                         -nodes$n_targets, nodes$label), , drop = FALSE]
    key <- paste(assoc$target_id, assoc$disease_class, sep = "\r")
    first <- !duplicated(key)
    cl_links <- data.frame(target_id = assoc$target_id[first], node_id = assoc$disease_class[first],
                           score = as.vector(tapply(assoc$association_score, key, max_na)[key[first]]),
                           focus = FALSE, stringsAsFactors = FALSE)
    fa <- assoc[assoc$focus, , drop = FALSE]
    fo_links <- data.frame(target_id = fa$target_id, node_id = fa$disease_id, score = fa$association_score,
                           focus = rep(TRUE, nrow(fa)), stringsAsFactors = FALSE)
    links <- rbind(cl_links, fo_links)
  } else {
    top <- utils::head(dis[!dis$focus, , drop = FALSE], top_n)
    top <- data.frame(node_id = top$disease_id, label = top$disease_name, disease_class = top$disease_class,
                      type = rep("disease", nrow(top)), n_targets = top$n_targets, n_diseases = rep(1L, nrow(top)),
                      mean_score = top$mean_score, stringsAsFactors = FALSE)
    nodes <- rbind(top, focus_nodes)
    rk <- dis$rank[match(nodes$node_id, dis$disease_id)]
    nodes <- nodes[order(match(nodes$disease_class, class_order), nodes$type != "focus", rk), , drop = FALSE]
    sel <- assoc[assoc$disease_id %in% nodes$node_id, , drop = FALSE]
    links <- data.frame(target_id = sel$target_id, node_id = sel$disease_id, score = sel$association_score,
                        focus = sel$focus, stringsAsFactors = FALSE)
  }
  nodes$group <- nodes$disease_class
  rownames(nodes) <- NULL
  rownames(links) <- NULL
  list(nodes = nodes, links = links)
}

#' Three-column layout (compounds | targets | right nodes)
#'
#' @description
#' Right nodes keep their order, spaced evenly with a gap of `gap` slots
#' between groups (`group` column). Targets are ordered by the
#' score-weighted mean height of their right-hand nodes, compounds by the
#' mean height of their drawn targets (barycentric crossing reduction,
#' Sugiyama et al. 1981); compounds with no drawn target go last. Heights
#' run from 1 (top) to 0.
#' @return `list(compounds = data.frame(name, y, isolated), targets =
#'   targets + y, right = right + y)`.
#' @keywords internal
.disease_map_columns_layout <- function(ct_drawn, compounds, targets, right, links, gap = 0.6) {
  spread <- function(n) if (n <= 1) rep(0.5, n) else 1 - (seq_len(n) - 1) / (n - 1)
  if (nrow(right) > 0) {
    grp <- cumsum(c(TRUE, right$group[-1] != right$group[-nrow(right)]))
    slot <- seq_len(nrow(right)) - 1 + gap * (grp - 1)
    right$y <- if (max(slot) > 0) 1 - slot / max(slot) else rep(0.5, nrow(right))
  } else {
    right$y <- numeric(0)
  }
  ry <- stats::setNames(right$y, right$node_id)
  lk <- links[links$target_id %in% targets$uniprot_id & links$node_id %in% right$node_id, , drop = FALSE]
  w <- ifelse(is.na(lk$score), 0.5, lk$score)
  bary <- tapply(w * ry[lk$node_id], lk$target_id, sum) / tapply(w, lk$target_id, sum)
  tb <- as.vector(bary[targets$uniprot_id])
  tb[is.na(tb)] <- -1
  targets <- targets[order(-tb, -targets$n_compounds, targets$uniprot_id), , drop = FALSE]
  targets$y <- spread(nrow(targets))
  ty <- stats::setNames(targets$y, targets$uniprot_id)
  ce <- ct_drawn[ct_drawn$uniprot_id %in% targets$uniprot_id, , drop = FALSE]
  cb <- tapply(ty[ce$uniprot_id], ce$compound_id, mean)
  cmp <- data.frame(name = compounds, bary = as.vector(cb[compounds]), stringsAsFactors = FALSE)
  cmp$isolated <- is.na(cmp$bary)
  cmp <- cmp[order(cmp$isolated, -ifelse(is.na(cmp$bary), 0, cmp$bary), cmp$name), , drop = FALSE]
  cmp$y <- spread(nrow(cmp))
  rownames(cmp) <- NULL
  rownames(targets) <- NULL
  list(compounds = cmp, targets = targets, right = right)
}

#' Smooth (cubic Hermite, horizontal tangents) curves between column pairs
#' @return `data.frame(edge, x, y)` plus every column of `attrs` repeated
#'   per point.
#' @keywords internal
.disease_map_curves <- function(x0, y0, x1, y1, attrs = NULL, n = 32) {
  k <- length(x0)
  if (k == 0) return(NULL)
  t <- seq(0, 1, length.out = n)
  s <- t * t * (3 - 2 * t)
  out <- data.frame(edge = rep(seq_len(k), each = n),
                    x = rep(x0, each = n) + rep(x1 - x0, each = n) * rep(t, k),
                    y = rep(y0, each = n) + rep(y1 - y0, each = n) * rep(s, k))
  if (!is.null(attrs)) out <- cbind(out, attrs[rep(seq_len(k), each = n), , drop = FALSE])
  rownames(out) <- NULL
  out
}

#' `"class"` and `"disease"` views: compounds | targets | areas or diseases
#' @return `list(plot, tag, n_compounds, n_targets, n_right, n_edges, data)`.
#' @keywords internal
.disease_map_view_network <- function(proj, ct, assoc, sums, conditions, view = "class", top_n_diseases = 15,
                                      max_targets = 30, min_shared = NULL, top_hub_n = 30, palette = NULL,
                                      title_suffix = "", fig_width = 11) {
  rn <- .disease_map_right_nodes(assoc, sums, view = view, top_n = top_n_diseases)
  tg <- .disease_map_select_targets(ct, rn$links, max_targets = max_targets, min_shared = min_shared)
  compounds <- sort(unique(ct$compound_id))
  ct_drawn <- ct[ct$uniprot_id %in% tg$uniprot_id, , drop = FALSE]
  lay <- .disease_map_columns_layout(ct_drawn, compounds, tg, rn$nodes, rn$links)
  cmp <- lay$compounds
  tg <- lay$targets
  right <- lay$right
  links <- rn$links[rn$links$target_id %in% tg$uniprot_id, , drop = FALSE]

  ## Compound badges: all targets in scope, and how many carry any disease link.
  linked <- unique(assoc$target_id)
  cmp$n_targets <- as.integer(table(factor(ct$compound_id, levels = cmp$name)))
  cmp$n_linked <- as.integer(table(factor(ct$compound_id[ct$uniprot_id %in% linked], levels = cmp$name)))
  cmp$label <- .plot_truncate(.plot_unique_labels(proj, conditions, cmp$name, "compound"), 30)
  cmp$badge <- ifelse(cmp$isolated, sprintf("%d targets, none drawn", cmp$n_targets),
                      sprintf("%d targets, %d disease-linked", cmp$n_targets, cmp$n_linked))
  cmp$size <- 2.6 + 2.4 * sqrt(cmp$n_targets / max(1, cmp$n_targets))
  tg$label <- .plot_label_nodes(proj, conditions, tg$uniprot_id, rep("target", nrow(tg)))

  x_c <- 0
  x_t <- 1
  x_r <- 2
  ## Edges: compound -> target (thin grey), target -> area/disease (area colour).
  e1 <- ct_drawn
  e1 <- .disease_map_curves(rep(x_c, nrow(e1)), cmp$y[match(e1$compound_id, cmp$name)],
                            rep(x_t, nrow(e1)), tg$y[match(e1$uniprot_id, tg$uniprot_id)])
  links$disease_class <- right$disease_class[match(links$node_id, right$node_id)]
  links$score_w <- pmin(pmax(ifelse(is.na(links$score), 0.3, links$score), 0), 1)
  links$lty <- ifelse(links$focus, "22", "solid")
  links <- links[order(links$focus, links$score_w), , drop = FALSE]
  e2 <- .disease_map_curves(rep(x_t, nrow(links)), tg$y[match(links$target_id, tg$uniprot_id)],
                            rep(x_r, nrow(links)), right$y[match(links$node_id, right$node_id)],
                            attrs = links[, c("disease_class", "score_w", "lty", "focus")])

  max_n <- max(1, right$n_targets)
  right$size <- ifelse(right$type == "class", 3 + 7 * sqrt(right$n_targets / max_n),
                       2.6 + 5 * sqrt(right$n_targets / max_n))
  right$size[right$type == "focus"] <- pmax(right$size[right$type == "focus"], 3.5)
  right$name_lab <- ifelse(right$type == "focus", paste0(.plot_truncate(right$label, 34), "  [focus]"),
                           .plot_truncate(right$label, 38))
  score_txt <- ifelse(is.na(right$mean_score), "", sprintf(", mean score %.2f", right$mean_score))
  right$stat_lab <- ifelse(right$type == "class",
                           sprintf("%d targets, %d diseases", right$n_targets, right$n_diseases),
                           sprintf("%d target%s%s", right$n_targets, ifelse(right$n_targets == 1, "", "s"), score_txt))
  right$lx <- x_r + 0.035 + 0.0045 * right$size

  n_all_targets <- length(unique(ct$uniprot_id))
  n_linked_all <- sum(unique(ct$uniprot_id) %in% linked)
  n_lo <- if (nrow(tg) > 0) min(tg$n_compounds) else 1
  n_hi <- max(length(compounds), n_lo + 1)
  size_breaks <- unique(round(seq(n_lo, n_hi, length.out = min(4, n_hi - n_lo + 1))))
  size_name <- "Compounds\nsharing target"
  cols <- .network_view_colours()

  p <- ggplot2::ggplot()
  if (!is.null(e1)) {
    p <- p + ggplot2::geom_path(data = e1, ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge),
                                colour = "grey50", linewidth = 0.22,
                                alpha = min(0.35, max(0.1, 60 / max(1, nrow(ct_drawn)))))
  }
  if (!is.null(e2)) {
    p <- p + ggplot2::geom_path(
      data = e2,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge, colour = .data$disease_class,
                   linewidth = .data$score_w, linetype = .data$lty),
      alpha = 0.55, lineend = "round"
    )
  }
  p <- p +
    ggplot2::scale_linetype_identity() +
    ggplot2::scale_linewidth(range = c(0.15, 1.5), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                             name = "Association\nscore")
  ## Compounds: diamonds.
  p <- p + ggplot2::geom_point(data = cmp, ggplot2::aes(x = x_c, y = .data$y), shape = 23, size = cmp$size,
                               fill = ifelse(cmp$isolated, "grey75", cols[["compound"]]), colour = "white", stroke = 0.5)
  ## Targets: viridis circles, size = compounds sharing.
  if (nrow(tg) > 0) {
    p <- p + ggplot2::geom_point(data = tg, ggplot2::aes(x = x_t, y = .data$y, size = .data$n_compounds,
                                                         fill = .data$n_compounds),
                                 shape = 21, colour = "grey25", stroke = 0.3)
  }
  ## Right nodes: a white (or black, focus) backing point, then the coloured one.
  foc <- right[right$type == "focus", , drop = FALSE]
  nfoc <- right[right$type != "focus", , drop = FALSE]
  if (nrow(nfoc) > 0) {
    p <- p +
      ggplot2::geom_point(data = nfoc, ggplot2::aes(x = x_r, y = .data$y), shape = 16, size = nfoc$size + 0.9,
                          colour = "white") +
      ggplot2::geom_point(data = nfoc, ggplot2::aes(x = x_r, y = .data$y, colour = .data$disease_class),
                          shape = 16, size = nfoc$size)
  }
  if (nrow(foc) > 0) {
    p <- p +
      ggplot2::geom_point(data = foc, ggplot2::aes(x = x_r, y = .data$y), shape = 15, size = foc$size + 1.3,
                          colour = "grey10") +
      ggplot2::geom_point(data = foc, ggplot2::aes(x = x_r, y = .data$y, colour = .data$disease_class),
                          shape = 15, size = foc$size)
  }

  ## Labels.
  lab_t <- utils::head(tg[order(-tg$n_compounds, -ifelse(is.na(tg$mean_weight), -Inf, tg$mean_weight)), , drop = FALSE], top_hub_n)
  if (nrow(lab_t) > 0) {
    ## Targets are evenly spaced in their column, so a label pinned just
    ## right of its node never collides with another; the white box keeps
    ## it readable over the outgoing edges.
    lab_t$lx <- x_t + 0.04
    p <- p + ggplot2::geom_label(
      data = lab_t, ggplot2::aes(x = .data$lx, y = .data$y, label = .data$label),
      size = 2.35, colour = "grey10", fontface = "italic", hjust = 0, vjust = 0.5,
      fill = grDevices::adjustcolor("white", alpha.f = 0.85), label.size = 0,
      label.padding = ggplot2::unit(0.08, "lines"), label.r = ggplot2::unit(0.05, "lines"),
      show.legend = FALSE
    )
  }
  p <- p +
    ggplot2::geom_text(data = cmp, ggplot2::aes(x = x_c - 0.06, y = .data$y, label = .data$label),
                       hjust = 1, vjust = -0.15, size = 2.75, fontface = "bold",
                       colour = ifelse(cmp$isolated, "grey50", "grey10")) +
    ggplot2::geom_text(data = cmp, ggplot2::aes(x = x_c - 0.06, y = .data$y, label = .data$badge),
                       hjust = 1, vjust = 1.25, size = 2.2, colour = "grey40")
  if (nrow(right) > 0) {
    p <- p +
      ggplot2::geom_text(data = right, ggplot2::aes(x = .data$lx, y = .data$y, label = .data$name_lab),
                         hjust = 0, vjust = -0.15, size = ifelse(right$type == "class", 3.1, 2.75),
                         fontface = "bold", colour = "grey10") +
      ggplot2::geom_text(data = right, ggplot2::aes(x = .data$lx, y = .data$y, label = .data$stat_lab),
                         hjust = 0, vjust = 1.25, size = 2.2, colour = "grey35")
  }
  ## Column headers.
  heads <- data.frame(
    x = c(x_c, x_t, x_r), y = 1.075,
    label = c(sprintf("Compounds (%d)", length(compounds)),
              sprintf("Target proteins (%d of %d drawn)", nrow(tg), n_all_targets),
              if (view == "class") "Disease areas" else "Diseases"),
    stringsAsFactors = FALSE
  )
  p <- p + ggplot2::geom_text(data = heads, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                              size = 3.1, fontface = "bold", colour = "grey30",
                              hjust = c(1, 0.5, 0), nudge_x = c(0.04, 0, -0.04))

  p <- p +
    ggplot2::scale_size(range = c(1.8, 4.6), limits = c(n_lo, n_hi), breaks = size_breaks, name = size_name) +
    ggplot2::scale_fill_viridis_c(option = "D", direction = -1, end = 0.92, limits = c(n_lo, n_hi),
                                  breaks = size_breaks, name = size_name, guide = "legend") +
    ggplot2::scale_colour_manual(
      values = palette, name = "Disease area",
      guide = if (view == "disease") ggplot2::guide_legend(override.aes = list(shape = 16, size = 3.5, linewidth = 0, alpha = 1), order = 3, ncol = 2) else "none"
    ) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1), fill = ggplot2::guide_legend(order = 1),
                    linewidth = ggplot2::guide_legend(order = 2, override.aes = list(colour = "grey40", alpha = 0.8))) +
    ggplot2::scale_x_continuous(limits = c(-0.85, 2.95), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-0.04, 1.1), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off")

  n_iso <- sum(cmp$isolated)
  pick <- if (is.null(min_shared)) sprintf("the %d disease-linked targets hit by the most compounds", nrow(tg))
          else sprintf("the %d disease-linked targets hit by >= %d compounds", nrow(tg), min_shared)
  right_txt <- if (view == "class") {
    sprintf(paste0("Right: %d disease areas (%d diseases, keyword grouping of disease names); node size = extract targets ",
                   "associated with the area (all %d disease-linked targets, not only the drawn ones)."),
            sum(right$type == "class"), nrow(sums$diseases), n_linked_all)
  } else {
    sprintf("Right: the %d diseases associated with the most extract targets (of %d), coloured by disease area; node size = targets.",
            sum(right$type == "disease"), nrow(sums$diseases))
  }
  sub <- paste0(
    sprintf("Middle: %s (size and colour = compounds sharing the target). ", pick), right_txt,
    " Edge width = association score (Open Targets).",
    if (nrow(foc) > 0) " Squares and dashed edges: focus diseases (custom gene sets)." else "",
    if (n_iso > 0) sprintf(" %d compound%s with no drawn target listed at the bottom.", n_iso, if (n_iso == 1) "" else "s") else ""
  )
  title <- if (view == "class") "Disease areas reached by the extract's targets" else "Top diseases linked to the extract's targets"
  p <- p +
    ggplot2::labs(
      title = .plot_wrap(paste0(title, " -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
      subtitle = .plot_wrap(sub, .plot_wrap_width(fig_width, 8)),
      caption = .plot_wrap(sprintf(
        "%d of %d predicted targets have a disease association (min. score filter applied). Disease areas: disease_map_classes().",
        n_linked_all, n_all_targets), .plot_wrap_width(fig_width, 7)),
      x = NULL, y = NULL
    ) +
    .network_view_theme(legend_position = "bottom") +
    ggplot2::theme(legend.box = "horizontal", legend.box.just = "top",
                   plot.margin = ggplot2::margin(8, 14, 6, 14))

  tag <- paste0(if (view == "disease") paste0("top", top_n_diseases, "_") else "",
                if (is.null(min_shared)) paste0("t", max_targets) else paste0("k", min_shared))
  list(plot = p, tag = tag, n_compounds = length(compounds), n_targets = nrow(tg), n_right = nrow(right),
       n_edges = nrow(ct_drawn) + nrow(links),
       data = list(compounds = cmp, targets = tg, right = right, links = links))
}

#' Count table behind the heatmap view
#'
#' @description
#' `rows = "compound"`: for each compound and area (and focus disease), the
#' number of distinct targets of the compound associated with it.
#' `rows = "target"`: for each drawn target, the number of the area's
#' diseases it is associated with (1 for a focus disease it belongs to).
#' Columns: areas with at least one associated target in scope (sorted by
#' targets, `"Other"` last), then the focus diseases.
#' @return `data.frame(row_id, col_id, col_label, block, n)` (all
#'   combinations, zeros included).
#' @keywords internal
.disease_map_heatmap_counts <- function(ct, assoc, sums, rows = "compound", row_ids) {
  cls <- sums$classes$disease_class[sums$classes$n_targets > 0]
  foc <- sums$diseases[sums$diseases$focus, , drop = FALSE]
  cols <- data.frame(col_id = c(cls, foc$disease_id), col_label = c(cls, foc$disease_name),
                     block = c(rep("Disease area", length(cls)), rep("Focus disease", nrow(foc))),
                     stringsAsFactors = FALSE)
  if (rows == "compound") {
    pairs <- merge(ct[, c("compound_id", "uniprot_id")], assoc[, c("target_id", "disease_id", "disease_class", "focus")],
                   by.x = "uniprot_id", by.y = "target_id")
    a <- unique(data.frame(row_id = pairs$compound_id, col_id = pairs$disease_class, item = pairs$uniprot_id))
    f <- pairs[pairs$focus, , drop = FALSE]
    b <- unique(data.frame(row_id = f$compound_id, col_id = f$disease_id, item = f$uniprot_id))
  } else {
    s <- assoc[assoc$target_id %in% row_ids, , drop = FALSE]
    a <- unique(data.frame(row_id = s$target_id, col_id = s$disease_class, item = s$disease_id))
    f <- s[s$focus, , drop = FALSE]
    b <- unique(data.frame(row_id = f$target_id, col_id = f$disease_id, item = f$disease_id))
  }
  ab <- rbind(a, b)
  grid <- expand.grid(row_id = row_ids, col_id = cols$col_id, stringsAsFactors = FALSE)
  cnt <- table(factor(ab$row_id, levels = row_ids), factor(ab$col_id, levels = cols$col_id))
  grid$n <- as.integer(cnt[cbind(match(grid$row_id, row_ids), match(grid$col_id, cols$col_id))])
  grid$col_label <- cols$col_label[match(grid$col_id, cols$col_id)]
  grid$block <- cols$block[match(grid$col_id, cols$col_id)]
  grid
}

#' `"heatmap"` view: compounds (or targets) x disease areas
#' @return `list(plot, tag, n_compounds, n_targets, n_right, n_edges, data)`.
#' @keywords internal
.disease_map_view_heatmap <- function(proj, ct, assoc, sums, conditions, rows = "compound", max_targets = 30,
                                      min_shared = NULL, title_suffix = "", fig_width = 11) {
  compounds <- sort(unique(ct$compound_id))
  if (rows == "compound") {
    row_ids <- .network_view_compound_order(ct)
    row_lab <- .plot_truncate(.plot_unique_labels(proj, conditions, row_ids, "compound"), 34)
    n_tot <- as.integer(table(factor(ct$compound_id, levels = row_ids)))
    row_lab <- sprintf("%s (%d)", row_lab, n_tot)
    tg <- NULL
  } else {
    links <- data.frame(target_id = assoc$target_id, node_id = assoc$disease_class, stringsAsFactors = FALSE)
    tg <- .disease_map_select_targets(ct, links, max_targets = max_targets, min_shared = min_shared)
    row_ids <- tg$uniprot_id
    row_lab <- sprintf("%s (%d cpd)", .plot_label_nodes(proj, conditions, row_ids, rep("target", length(row_ids))),
                       tg$n_compounds)
  }
  if (length(row_ids) == 0) cli::cli_abort("No row to draw in the heatmap.")
  hm <- .disease_map_heatmap_counts(ct, assoc, sums, rows = rows, row_ids = row_ids)
  col_ids <- unique(hm$col_id)
  cls_n <- sums$classes$n_targets[match(col_ids, sums$classes$disease_class)]
  dis_n <- sums$diseases$n_targets[match(col_ids, sums$diseases$disease_id)]
  col_n <- ifelse(is.na(cls_n), dis_n, cls_n)
  col_lab <- .plot_wrap(.plot_truncate(hm$col_label[match(col_ids, hm$col_id)], 34), 16)
  col_lab <- sprintf("%s\n(%d)", col_lab, col_n)
  hm$col_f <- factor(hm$col_id, levels = col_ids, labels = make.unique(col_lab))
  hm$row_f <- factor(hm$row_id, levels = rev(row_ids), labels = rev(make.unique(row_lab)))
  hm$block <- factor(hm$block, levels = c("Disease area", "Focus disease"))
  hm$fill <- ifelse(hm$n > 0, hm$n, NA)
  hi <- max(1, hm$n)
  hm$txt_col <- ifelse(hm$n > 0.55 * hi, "white", "grey10")
  hm$txt <- ifelse(hm$n > 0, as.character(hm$n), "")
  fill_name <- if (rows == "compound") "Targets associated" else "Diseases associated"

  p <- ggplot2::ggplot(hm, ggplot2::aes(x = .data$col_f, y = .data$row_f)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$txt), colour = hm$txt_col, size = 2.3) +
    ggplot2::scale_fill_viridis_c(option = "D", direction = -1, end = 0.95, na.value = "grey95",
                                  limits = c(1, hi), name = fill_name) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$block), scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(position = "top", expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 50, hjust = 0, vjust = 0, size = 7.5, colour = "grey15",
                                              lineheight = 0.9),
      axis.text.y = ggplot2::element_text(size = 7.5, colour = "grey15"),
      strip.text = ggplot2::element_text(face = "bold", size = 8.5, colour = "grey25"),
      strip.placement = "outside",
      panel.spacing = ggplot2::unit(0.6, "lines"),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30", lineheight = 1.1),
      plot.caption = ggplot2::element_text(size = 7, colour = "grey40", hjust = 0),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(8, 40, 8, 10)
    )
  sub <- if (rows == "compound") {
    sprintf(paste0("Tile = number of the compound's predicted targets associated with at least one disease of the area ",
                   "(or with the focus disease). Rows: %d compounds (total targets in brackets), ordered by target-set ",
                   "similarity; columns: disease areas (extract targets associated, in brackets) and focus diseases."),
            length(row_ids))
  } else {
    sprintf(paste0("Tile = number of the area's diseases the target is associated with (1 = member of the focus gene set). ",
                   "Rows: the %d disease-linked targets hit by the most compounds (compounds in brackets)."),
            length(row_ids))
  }
  p <- p + ggplot2::labs(
    title = .plot_wrap(paste0("Disease-area map of the extract's targets -- ", title_suffix), .plot_wrap_width(fig_width, 13)),
    subtitle = .plot_wrap(sub, .plot_wrap_width(fig_width, 8)),
    caption = "Disease areas: keyword grouping of Open Targets disease names, see disease_map_classes()."
  )
  tag <- if (rows == "compound") "compound" else paste0("target_", if (is.null(min_shared)) paste0("t", max_targets) else paste0("k", min_shared))
  list(plot = p, tag = tag, n_compounds = length(compounds), n_targets = if (is.null(tg)) length(unique(ct$uniprot_id)) else nrow(tg),
       n_right = length(col_ids), n_edges = sum(hm$n > 0), data = list(counts = hm))
}
