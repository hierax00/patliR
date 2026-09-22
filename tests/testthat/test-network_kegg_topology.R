## network_kegg_topology()'s own KGML fetch (KEGGREST::keggGet) and the bulk
## UniProt mapping (KEGGREST::keggConv) are live-network calls, not mocked
## here -- same convention as targets_disease_filter()/targets_disease_profile()
## (see test-targets_disease.R). The function was exercised against the real
## KEGG REST API by hand before being committed (see NEWS.md); a real KGML
## response (hsa04151) was fetched live to ground the fixture below.
##
## .network_kegg_parse_kgml() itself is pure XML parsing with no network
## call, so it IS fully unit-tested here against a hand-built KGML string
## whose shape mirrors the real hsa04151 response exactly: an entry can
## bundle several KEGG gene IDs (an ortholog/paralog box), a relation's
## entry1/entry2 are KGML-internal entry IDs (not gene IDs) that must be
## resolved via the entry table, a non-gene entry (e.g. type="compound")
## must exclude its relation, and an unrequested relation type (here
## "maplink", a diagram cross-link, never a biological relation) must be
## dropped before the gene-pair check.

.kegg_fixture_kgml <- function() {
  '<?xml version="1.0"?>
<pathway name="path:hsa00000" org="hsa" number="00000" title="Test Pathway">
    <entry id="1" name="hsa:100 hsa:200" type="gene">
        <graphics name="GENE1, GENE2" type="rectangle" x="10" y="10" width="20" height="10"/>
    </entry>
    <entry id="2" name="hsa:300" type="gene">
        <graphics name="GENE3" type="rectangle" x="50" y="50" width="20" height="10"/>
    </entry>
    <entry id="3" name="cpd:C00001" type="compound">
        <graphics name="C00001" type="circle" x="80" y="80" width="8" height="8"/>
    </entry>
    <relation entry1="1" entry2="2" type="PPrel">
        <subtype name="activation" value="--&gt;"/>
    </relation>
    <relation entry1="2" entry2="3" type="PPrel">
        <subtype name="inhibition" value="--|"/>
    </relation>
    <relation entry1="1" entry2="2" type="maplink">
        <subtype name="compound" value="C99999"/>
    </relation>
</pathway>'
}

test_that(".network_kegg_parse_kgml() resolves multi-gene entries, drops non-gene endpoints and unrequested relation types", {
  testthat::skip_if_not_installed("xml2")
  parsed <- patliR:::.network_kegg_parse_kgml(.kegg_fixture_kgml(), c("PPrel", "GErel", "ECrel", "PCrel"))

  expect_equal(parsed$title, "Test Pathway")
  rel <- parsed$relations
  ## entry2="3" is type="compound", not a gene -- its PPrel relation is dropped;
  ## the maplink relation is dropped by relation_types before the gene check;
  ## only entry1=1 -> entry2=2 survives, cartesian-expanded over entry 1's two genes.
  expect_equal(nrow(rel), 2)
  expect_setequal(rel$from_kegg, c("hsa:100", "hsa:200"))
  expect_true(all(rel$to_kegg == "hsa:300"))
  expect_true(all(rel$relation_type == "PPrel"))
  expect_true(all(rel$relation_subtype == "activation"))
  expect_true(all(rel$relation_value == "-->"))
})

test_that(".network_kegg_parse_kgml() restricts to the requested relation_types", {
  testthat::skip_if_not_installed("xml2")
  parsed <- patliR:::.network_kegg_parse_kgml(.kegg_fixture_kgml(), "GErel")
  expect_equal(nrow(parsed$relations), 0)
})

test_that(".network_kegg_parse_kgml() returns an empty relations frame with the right columns when a pathway has no relations", {
  testthat::skip_if_not_installed("xml2")
  no_rel_kgml <- '<?xml version="1.0"?>
<pathway name="path:hsa00001" org="hsa" number="00001" title="No Relations">
    <entry id="1" name="hsa:1" type="gene"><graphics name="G1" type="rectangle"/></entry>
</pathway>'
  parsed <- patliR:::.network_kegg_parse_kgml(no_rel_kgml, c("PPrel", "GErel", "ECrel", "PCrel"))
  expect_equal(nrow(parsed$relations), 0)
  expect_setequal(names(parsed$relations), c("from_kegg", "to_kegg", "relation_type", "relation_subtype", "relation_value"))
})

test_that("network_kegg_topology() requires network_build() to have run first", {
  proj <- .test_project()
  proj <- prep_compounds(proj, .test_compound_list(), identifier = "pubchem")
  expect_error(network_kegg_topology(proj, pathway_ids = "hsa00000"), "network_edges")
})

test_that("network_kegg_topology() aborts clearly when no pathway_ids are given and there is no KEGG enrichment", {
  proj <- .network_stats_test_setup()
  expect_error(network_kegg_topology(proj, condition = "FLO-ET"), "network_enrich")
})

test_that("network_kegg_topology() validates species and restrict_to_network", {
  proj <- .network_stats_test_setup()
  expect_error(network_kegg_topology(proj, condition = "FLO-ET", pathway_ids = "hsa00000", species = character(0)))
  expect_error(network_kegg_topology(proj, condition = "FLO-ET", pathway_ids = "hsa00000", restrict_to_network = "yes"))
})
