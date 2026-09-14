#' @include AllGenerics.R internal.R
NULL

#' Compute physicochemical/ADME properties locally, no network required
#'
#' @description
#' Computes molecular descriptors with `rcdk` for every compound in
#' [compounds()] (or a subset) and derives: the Lipinski Rule of Five
#' (Ro5), the Veber rule, the Ghose filter, the Egan filter, the Oprea
#' lead-likeness rule, a BOILED-Egg-style estimate of passive GI
#' absorption / BBB permeation, and one physicochemical compatibility flag
#' per administration route.
#'
#' @section Drug-likeness vs. lead-likeness rules, and which ones patliR does *not* implement:
#' Five binary pass/fail rules are computed, each a real published
#' criterion:
#' \itemize{
#'   \item **Lipinski Ro5** (`ro5_pass`) -- Lipinski et al. (2001), *Adv.
#'     Drug Deliv. Rev.* 46, 3-26: MW <= 500, logP <= 5, HBD <= 5, HBA <= 10.
#'   \item **Veber** (`veber_pass`) -- Veber et al. (2002), *J. Med. Chem.*
#'     45, 2615-2623: TPSA <= 140, rotatable bonds <= 10.
#'   \item **Ghose** (`ghose_pass`) -- Ghose et al. (1999), *J. Comb.
#'     Chem.* 1, 55-68: 160 <= MW <= 480, -0.4 <= logP <= 5.6, 20 <= heavy
#'     atoms <= 70, 40 <= molar refractivity (AMR) <= 130. All four
#'     criteria now checked (an earlier version of this function omitted
#'     the AMR term).
#'   \item **Egan** (`egan_pass`) -- Egan et al. (2000), *J. Med. Chem.*
#'     43, 3867-3877 (the "egg" model that later BOILED-Egg extends):
#'     logP <= 5.88, TPSA <= 131.6.
#'   \item **Oprea lead-likeness** (`oprea_pass`) -- Oprea (2000), *J.
#'     Comput. Aided Mol. Des.* 14, 251-264: HBD < 2, 2 < HBA < 10,
#'     2 < rotatable bonds < 8, 1 < ring count < 4. This is the standard
#'     literature "lead-like properties" rule (narrower than drug-likeness,
#'     meant for hit-to-lead triage rather than a final candidate) -- see
#'     `adme_filter(rules = "oprea")`.
#' }
#' **Deliberately not implemented**: Hughes et al. (2008, *Bioorg. Med.
#' Chem. Lett.* 18, 4872-4875), Ritchie & Macdonald (2009, *Drug Discov.
#' Today* 14, 1011-1020) and Lovering et al. (2009, *J. Med. Chem.* 52,
#' 6752-6756) are *correlational findings*, not binary pass/fail rules the
#' original papers define with a cutoff -- a hard cutoff for them would
#' mean inventing a threshold the source literature does not give. Muegge
#' (2001, *J. Med. Chem.* 44, 1841-1846) is a real pass/fail rule but its
#' exact thresholds could not be re-confirmed against the (paywalled)
#' primary source -- deferred rather than transcribed from memory.
#'
#' @section BOILED-Egg (please read the WLogP caveat):
#' `gi_absorption`/`bbb_permeant` are now a real point-in-polygon test
#' against the published BOILED-Egg ellipses (Daina & Zoete, 2016,
#' *ChemMedChem* 11, 1117-1121), digitized in `inst/extdata/boiled_egg_*.csv`
#' -- see [plot_boiled_egg()] for the plot and the exact provenance note.
#' The one remaining approximation: the original model is defined on
#' **WLogP** (Wildman & Crippen's atom-contribution logP, as implemented in
#' RDKit); `rcdk`/CDK has no identical implementation, so `wlogp_proxy` uses
#' CDK's own ALogP (`rcdk::get.alogp()`, Ghose-Crippen-style, a different
#' but methodologically related atom-contribution method) instead. This can
#' shift a compound slightly relative to where SwissADME itself would place
#' it -- treat `gi_absorption`/`bbb_permeant` as "very likely right", not
#' "certified identical to SwissADME's own call".
#'
#' @section Route criteria and their references:
#' \itemize{
#'   \item **Oral** (Ro5, Lipinski): MW <= 500, logP <= 5, HBD <= 5, HBA <= 10.
#'   \item **Injectable**: approximate logD(pH 7.4) (here taken as the
#'     computed logP, since `patliR` does not model ionization) between 1
#'     and 3.
#'   \item **Ophthalmic** (Karami et al. 2022, *J Ocul Pharmacol Ther*):
#'     TPSA <= 250 sq. Angstrom and approximate clogD(pH 7.4) <= 4.0 (see
#'     the injectable caveat above -- same approximation applies).
#'   \item **Topical/dermal**: MW <= 500 (optimum <= 400) and logP between 1
#'     and 4 (optimum 1-3).
#' }
#' None of these route flags are a formulation or regulatory
#' recommendation -- treat them as a coarse physicochemical compatibility
#' screen, not a delivery-route decision.
#'
#' @inheritParams compounds
#' @param compound_ids Character vector of `compounds(proj)$id`, or `NULL`
#'   (default) for every compound currently in `proj`.
#' @param routes Character vector, any of `"oral"`, `"topical"`,
#'   `"ophthalmic"`, `"injectable"`.
#'
#' @section Radar-chart descriptors (please also read):
#' `fraction_csp3_approx` and `aromatic_proportion_approx` are computed with
#' a lightweight regex heuristic over the SMILES string itself (lowercase
#' aromatic atoms vs. uppercase `C`), **not** true CDK hybridization/
#' aromaticity perception -- `fraction_csp3_approx` therefore over-counts
#' non-aromatic sp2 carbons (e.g. in C=C or C=O) as if they were sp3. It is
#' good enough to place a compound roughly on [plot_admet_radar()]'s
#' INSATU axis, not for anything that needs a rigorous Fsp3. `logs_esol` is
#' the published Delaney (2004) ESOL equation applied to `logp`, `mw`,
#' `rotatable_bonds`, and `aromatic_proportion_approx` -- so it inherits
#' the same approximation for its aromatic-proportion term. `n_rings_approx`
#' counts SMILES ring-closure digit pairs (single digits and `%nn`
#' two-digit forms, outside `[...]` brackets so isotope/charge numbers
#' are not miscounted) -- deterministic per the SMILES specification for
#' well-formed strings, but not true CDK graph-theoretic SSSR ring
#' perception (`rcdk` has no high-level wrapper for that; see the code
#' comment where it is computed). Used by `oprea_pass` below.
#'
#' @return The updated `proj`, with an `adme_local` entry in
#'   [patliRResults()] (columns `compound_id`, `mw`, `logp`, `hbd`, `hba`,
#'   `tpsa`, `rotatable_bonds`, `n_atoms`, `amr`, `wlogp_proxy`,
#'   `fraction_csp3_approx`, `aromatic_proportion_approx`,
#'   `n_rings_approx`, `logs_esol`, `ro5_pass`, `veber_pass`, `ghose_pass`,
#'   `egan_pass`, `oprea_pass`, `gi_absorption`, `bbb_permeant`,
#'   `route_oral`, `route_topical`, `route_ophthalmic`,
#'   `route_injectable`), also written to `results/adme_local.csv`.
#'
#' @examples
#' \donttest{
#' proj <- patliR_project(tempfile("patliR_demo_"))
#' compound_list <- read.csv(
#'   system.file("extdata", "input_compound_list.csv", package = "patliR")
#' )
#' proj <- prep_compounds(proj, compound_list, identifier = "pubchem")
#' proj <- adme_local(proj)
#' patliRResults(proj, "adme_local")
#' }
#'
#' @export
adme_local <- function(proj, compound_ids = NULL,
                        routes = c("oral", "topical", "ophthalmic", "injectable")) {
  stopifnot(is(proj, "PatliRProject"))
  routes <- match.arg(routes, several.ok = TRUE)

  cmp <- compounds(proj)
  if (!is.null(compound_ids)) cmp <- cmp[cmp$id %in% compound_ids, , drop = FALSE]
  if (nrow(cmp) == 0) {
    cli::cli_warn("No matching compounds in {.arg proj}; run {.fn prep_compounds} first. Nothing to do.")
    return(proj)
  }

  desc <- .compute_adme_descriptors(cmp$smiles)
  failed <- is.na(desc$mw)
  if (any(failed)) {
    proj <- .log_append(
      proj, step = "adme_local", id = cmp$id[failed],
      message = "descriptor calculation failed (structure could not be re-parsed by rcdk)"
    )
  }

  out <- data.frame(
    compound_id = cmp$id,
    mw = desc$mw, logp = desc$logp, hbd = desc$hbd, hba = desc$hba,
    tpsa = desc$tpsa, rotatable_bonds = desc$rotatable_bonds,
    stringsAsFactors = FALSE
  )

  out$n_atoms <- desc$n_atoms

  approx <- .smiles_approx_descriptors(cmp$smiles, desc$n_atoms)
  out$fraction_csp3_approx <- approx$fraction_csp3_approx
  out$aromatic_proportion_approx <- approx$aromatic_proportion_approx
  out$n_rings_approx <- approx$n_rings_approx
  ## Delaney (2004) ESOL equation, J Chem Inf Comput Sci 44(3):1000-1005.
  out$logs_esol <- with(out,
    0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotatable_bonds - 0.74 * aromatic_proportion_approx
  )

  out$amr <- desc$amr

  out$ro5_pass   <- with(out, mw <= 500 & logp <= 5 & hbd <= 5 & hba <= 10)
  out$veber_pass <- with(out, tpsa <= 140 & rotatable_bonds <= 10)
  ## Ghose et al. 1999 -- now all four criteria, including molar
  ## refractivity (AMR); an earlier version of this function omitted AMR
  ## because it was not yet being computed.
  out$ghose_pass <- with(out, mw >= 160 & mw <= 480 & logp >= -0.4 & logp <= 5.6 &
    n_atoms >= 20 & n_atoms <= 70 & amr >= 40 & amr <= 130)
  ## Egan et al. 2000 ("egg" model).
  out$egan_pass <- with(out, logp <= 5.88 & tpsa <= 131.6)
  ## Oprea 2000 lead-likeness -- narrower than drug-likeness, meant for
  ## hit-to-lead triage. n_rings_approx is the SMILES ring-closure-digit
  ## heuristic (see .smiles_approx_descriptors()), not true CDK SSSR ring
  ## perception -- rcdk has no high-level ring-count wrapper, and this
  ## avoids guessing an unverified CDK descriptor column name (see the
  ## code comment there for exactly what it does and does not handle).
  out$oprea_pass <- with(out, hbd < 2 & hba > 2 & hba < 10 &
    rotatable_bonds > 2 & rotatable_bonds < 8 &
    n_rings_approx > 1 & n_rings_approx < 4)

  out$wlogp_proxy <- desc$wlogp_proxy
  egg <- .boiled_egg(out$tpsa, out$wlogp_proxy)
  out$gi_absorption <- egg$gi_absorption
  out$bbb_permeant  <- egg$bbb_permeant

  if ("oral" %in% routes) {
    out$route_oral <- out$ro5_pass
  }
  if ("injectable" %in% routes) {
    out$route_injectable <- with(out, logp >= 1 & logp <= 3)
  }
  if ("ophthalmic" %in% routes) {
    out$route_ophthalmic <- with(out, tpsa <= 250 & logp <= 4.0)
  }
  if ("topical" %in% routes) {
    out$route_topical <- with(out, mw <= 500 & logp >= 1 & logp <= 4)
  }

  out <- .network_upsert(proj, "adme_local", out, "compound_id")
  patliRResults(proj, "adme_local") <- out
  .write_results_csv(proj, "adme_local", out)
  .write_log_csv(proj)
  proj
}

#' @keywords internal
.compute_adme_descriptors <- function(smiles) {
  mols <- .parse_smiles_safe(smiles)
  empty_row <- data.frame(mw = NA_real_, logp = NA_real_, hbd = NA_integer_,
                           hba = NA_integer_, tpsa = NA_real_,
                           rotatable_bonds = NA_integer_, n_atoms = NA_integer_,
                           wlogp_proxy = NA_real_, amr = NA_real_)
  rows <- lapply(mols, function(m) {
    if (is.null(m)) return(empty_row)
    tryCatch({
      rcdk::convert.implicit.to.explicit(m)
      data.frame(
        mw  = .safe_mw(m),
        logp = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "XLogP"),
        hbd = as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor", "nHBDon")),
        hba = as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor", "nHBAcc")),
        tpsa = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor", "TopoPSA"),
        rotatable_bonds = as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor", "nRotB")),
        n_atoms = as.integer(tryCatch(rcdk::get.atom.count(m), error = function(e) NA_integer_)),
        wlogp_proxy = tryCatch(as.numeric(rcdk::get.alogp(m)[["ALogP"]]), error = function(e) NA_real_),
        ## Same ALOGPDescriptor rcdk::get.alogp() already calls, but read via
        ## eval.desc() for the AMR (molar refractivity) column it also
        ## returns alongside ALogP/ALogp2 -- needed for the Ghose filter's
        ## fourth criterion (Ghose et al. 1999).
        amr = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor", "AMR")
      )
    }, error = function(e) empty_row)
  })
  do.call(rbind, rows)
}

#' @keywords internal
.safe_mw <- function(mol) {
  tryCatch(rcdk::get.exact.mass(mol), error = function(e) NA_real_)
}

#' @keywords internal
.safe_desc <- function(mol, desc_name, col) {
  tryCatch({
    val <- rcdk::eval.desc(mol, desc_name, verbose = FALSE)
    as.numeric(val[[col]])
  }, error = function(e) NA_real_)
}

#' Rough SMILES-regex proxies for Fsp3, aromatic proportion, and ring count
#'
#' See the "Radar-chart descriptors" note in [adme_local()] -- these are
#' cheap heuristics over the raw SMILES string, not true CDK hybridization/
#' aromaticity/ring perception.
#' @keywords internal
.smiles_approx_descriptors <- function(smiles, n_atoms) {
  vals <- Map(function(s, n_heavy) {
    if (is.na(s) || !nzchar(s) || is.na(n_heavy) || n_heavy == 0) {
      return(c(fraction_csp3_approx = NA_real_, aromatic_proportion_approx = NA_real_,
               n_rings_approx = NA_real_))
    }
    aromatic_c     <- lengths(regmatches(s, gregexpr("c", s, fixed = TRUE)))
    aromatic_other <- lengths(regmatches(s, gregexpr("[nosp](?![a-z])", s, perl = TRUE)))
    aliphatic_c    <- lengths(regmatches(s, gregexpr("C(?![a-z])", s, perl = TRUE)))
    total_carbon <- aromatic_c + aliphatic_c
    c(
      fraction_csp3_approx = if (total_carbon > 0) aliphatic_c / total_carbon else NA_real_,
      aromatic_proportion_approx = (aromatic_c + aromatic_other) / n_heavy,
      n_rings_approx = .smiles_ring_count_approx(s)
    )
  }, smiles, n_atoms)
  as.data.frame(do.call(rbind, vals))
}

#' Count ring-closure digit pairs in a SMILES string (approximate ring count)
#'
#' @description
#' Per the SMILES specification, a ring bond is written as a digit (or a
#' `%nn` two-digit escape) appearing exactly twice in the string -- once at
#' each ring-closure atom. Counting closure tokens and halving gives the
#' ring count for well-formed SMILES. Bracket contents (`[...]`, e.g.
#' isotope labels like `[13C]` or explicit H counts) are stripped first so
#' those digits are never mistaken for ring closures.
#'
#' @section What this does not handle:
#' Only correct for well-formed SMILES that follow the digit-pair
#' convention; does not attempt real graph-theoretic SSSR ring perception
#' (`rcdk` has no high-level wrapper for that -- see the "Radar-chart
#' descriptors" note in [adme_local()] for why a CDK descriptor column was
#' not guessed instead). Good enough for `oprea_pass`'s ring-count
#' criterion, not for anything needing an exact SSSR count.
#'
#' @return Numeric scalar (ring count), `NA_real_` if `s` is empty.
#' @keywords internal
.smiles_ring_count_approx <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_real_)
  no_brackets <- gsub("\\[[^]]*\\]", "", s)
  two_digit <- regmatches(no_brackets, gregexpr("%\\d{2}", no_brackets))[[1]]
  remainder <- gsub("%\\d{2}", "", no_brackets)
  single_digit <- regmatches(remainder, gregexpr("\\d", remainder))[[1]]
  floor((length(two_digit) + length(single_digit)) / 2)
}

#' Real BOILED-Egg GI absorption / BBB permeation call (point-in-polygon)
#'
#' See the "BOILED-Egg" note in [adme_local()] for the WLogP-proxy caveat;
#' see `.point_in_polygon()` and `.load_boiled_egg_polygons()` (internal)
#' for the published ellipse boundaries themselves.
#' @keywords internal
.boiled_egg <- function(tpsa, wlogp) {
  poly <- .load_boiled_egg_polygons()
  has_val <- !is.na(tpsa) & !is.na(wlogp)
  gi <- rep(NA, length(tpsa))
  bbb <- rep(NA, length(tpsa))
  if (any(has_val)) {
    gi[has_val] <- .point_in_polygon(tpsa[has_val], wlogp[has_val], poly$gia$tpsa, poly$gia$wlogp)
    bbb[has_val] <- .point_in_polygon(tpsa[has_val], wlogp[has_val], poly$bbb$tpsa, poly$bbb$wlogp)
  }
  list(gi_absorption = ifelse(is.na(gi), NA_character_, ifelse(gi, "High", "Low")), bbb_permeant = bbb)
}
