#' @include AllGenerics.R internal.R
NULL

#' Compute physicochemical/ADME properties locally, no network required
#'
#' @description
#' Computes molecular descriptors with `rcdk` for every compound in
#' [compounds()] (or a subset) and derives: the Lipinski Rule of Five
#' (Ro5), the Veber rule, the Ghose filter, the Egan filter, the Oprea
#' property ranges, a BOILED-Egg-style estimate of passive GI
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
#'     45, 2615-2623: TPSA <= 140, rotatable bonds <= 10, where (as in the
#'     paper) a rotatable bond is a single non-ring bond to a non-terminal
#'     heavy atom **excluding amide C-N bonds**; `rotatable_bonds` is counted
#'     that way (CDK `RotatableBondsCountDescriptor`, `excludeAmides = TRUE`).
#'   \item **Ghose** (`ghose_pass`) -- Ghose et al. (1999), *J. Comb.
#'     Chem.* 1, 55-68: 160 <= MW <= 480, -0.4 <= logP <= 5.6, 20 <= atoms
#'     <= 70 (`n_atoms` counts **all atoms including hydrogens**, as the
#'     original atom-count criterion does), 40 <= molar refractivity (AMR)
#'     <= 130. All four criteria checked.
#'   \item **Egan** (`egan_pass`) -- Egan et al. (2000), *J. Med. Chem.*
#'     43, 3867-3877 (the "egg" model that later BOILED-Egg extends):
#'     logP <= 5.88, TPSA <= 131.6.
#'   \item **Oprea property ranges** (`oprea_pass`) -- Oprea (2000), *J.
#'     Comput. Aided Mol. Des.* 14, 251-264: 0 <= HBD <= 2, 2 <= HBA <= 9,
#'     2 <= rotatable bonds <= 8, 1 <= ring count <= 4 (inclusive ranges;
#'     the paper reports them as the ranges holding ~70% of the *drug-like*
#'     compounds it surveyed, not as a stricter lead-likeness rule -- an
#'     earlier version used narrower, exclusive limits). See
#'     `adme_filter(rules = "oprea")`.
#' }
#' **Descriptor conventions.** `mw` is the average molecular weight (natural
#' isotopic abundance, `rcdk::get.natural.mass()`), not the monoisotopic mass.
#' `logp` is CDK's XLogP and `hba`/`hbd` are CDK's donor/acceptor counts, so
#' the rules above are adapted screens on CDK descriptors, not bit-for-bit
#' reproductions of the papers' own logP/HBA definitions. `ro5_pass` requires
#' all four Lipinski criteria (zero violations); the original paper also
#' tolerates one violation. Structures are used as supplied (no salt
#' stripping or neutralisation).
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
#' but methodologically related atom-contribution method) instead. The two
#' logP values can differ enough to move a compound across an ellipse
#' boundary, and this substitution has not been validated here against
#' SwissADME's own calls -- treat `gi_absorption`/`bbb_permeant` as an
#' approximate screen, most reliable for compounds well inside or well
#' outside the ellipses.
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
#' @param wlogp_source `"cdk"` (default, no Python ever invoked): `wlogp_proxy`
#'   is CDK's ALogP, an atom-contribution method related to but not identical
#'   to Wildman & Crippen's WLogP (see the BOILED-Egg section). `"rdkit"`:
#'   compute the actual Wildman & Crippen (1999) WLogP with `RDKit::Crippen.MolLogP()`,
#'   via an optional `reticulate` + Python `rdkit` environment -- exactly the
#'   quantity BOILED-Egg was built on, no proxy. Needs `reticulate` (Suggests)
#'   and a `rdkit` Python package `reticulate` can resolve (installed once,
#'   cached by `reticulate`/`uv` after the first call -- the *first* call with
#'   `wlogp_source = "rdkit"` may download it, which needs internet). If it
#'   cannot be resolved, `adme_local()` warns and falls back to `"cdk"`
#'   automatically -- this argument never turns a working call into an error.
#' @param rdkit_qc Logical, default `FALSE`. Adds four **extra, informational**
#'   columns from the same optional RDKit -- never replacing or filtering the
#'   CDK-based ones, purely a second opinion to look at side by side:
#'   `hba_lipinski_rdkit`/`hbd_lipinski_rdkit` (RDKit's
#'   `CalcNumLipinskiHBA()`/`CalcNumLipinskiHBD()`, the literal N+O count
#'   Lipinski et al. 2001 define -- CDK's `hba`/`hbd` use CDK's own donor/
#'   acceptor rules instead, which can disagree at the margins, e.g. an N-O
#'   bond), `tpsa_rdkit` (RDKit's `TPSA()`, the same Ertl 2000 fragment
#'   method as CDK's `tpsa` -- implementation cross-check) and
#'   `fraction_csp3_rdkit`/`aromatic_proportion_rdkit` (RDKit's own
#'   hybridization/aromaticity perception -- CDK's `fraction_csp3_approx`/
#'   `aromatic_proportion_approx` already come from CDK's graph perception,
#'   not a proxy, but the two toolkits' *aromaticity models* (which ring
#'   systems count as aromatic) can differ on edge cases, e.g. some fused
#'   heterocycles -- a large disagreement is worth a manual look, not
#'   evidence either one is wrong). Same graceful degradation as
#'   `wlogp_source = "rdkit"`: unavailable RDKit warns once and every
#'   `_rdkit` column is `NA`, nothing else changes.
#'
#' @section Radar-chart descriptors (please also read):
#' `fraction_csp3_approx` is CDK's `FractionalCSP3Descriptor` (sp3 carbons /
#' all carbons, from perceived hybridization) and `aromatic_proportion_approx`
#' is CDK's aromatic-atom count (after aromaticity perception) divided by the
#' heavy-atom count, so neither depends on how the SMILES was written (the
#' `_approx` suffix is kept for column-name stability; an earlier version
#' used a SMILES-string regex heuristic that mistook C=C and C=O carbons for
#' sp3). `logs_esol` is
#' the published Delaney (2004) ESOL equation applied to `logp` (CDK XLogP,
#' not the predictor the original coefficients were fitted with), `mw`,
#' `rotatable_bonds`, and `aromatic_proportion_approx`. `n_rings_approx`
#' counts SMILES ring-closure digit pairs (single digits and `%nn`
#' two-digit forms, outside `[...]` brackets so isotope/charge numbers
#' are not miscounted) -- deterministic per the SMILES specification for
#' well-formed strings, but not true CDK graph-theoretic SSSR ring
#' perception (`rcdk` has no high-level wrapper for that; see the code
#' comment where it is computed). Used by `oprea_pass` below.
#'
#' @return The updated `proj`, with an `adme_local` entry in
#'   [patliRResults()] (columns `compound_id`, `mw`, `logp`, `hbd`, `hba`,
#'   `tpsa`, `rotatable_bonds`, `n_atoms`, `amr`, `wlogp_proxy`, `wlogp_source`
#'   (`"cdk_alogp_proxy"` or `"rdkit_wlogp"`, per compound),
#'   `fraction_csp3_approx`, `aromatic_proportion_approx` (plus, when
#'   `rdkit_qc = TRUE`, `hba_lipinski_rdkit`, `hbd_lipinski_rdkit`,
#'   `tpsa_rdkit`, `fraction_csp3_rdkit`, `aromatic_proportion_rdkit`),
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
                        routes = c("oral", "topical", "ophthalmic", "injectable"),
                        wlogp_source = c("cdk", "rdkit"), rdkit_qc = FALSE) {
  stopifnot(is(proj, "PatliRProject"))
  routes <- match.arg(routes, several.ok = TRUE)
  wlogp_source <- match.arg(wlogp_source)
  .pathway_check_flag(rdkit_qc, "rdkit_qc")

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
  ## a single descriptor can fail while the mass succeeds: log those too, per descriptor
  for (nm in c("logp", "hbd", "hba", "tpsa", "wlogp_proxy", "amr", "rotatable_bonds", "fraction_csp3", "n_arom_atoms")) {
    partial <- !failed & is.na(desc[[nm]])
    if (any(partial)) {
      proj <- .log_append(proj, step = "adme_local", id = cmp$id[partial],
                          message = paste0("descriptor '", nm, "' could not be computed (left NA)"))
    }
  }

  out <- data.frame(
    compound_id = cmp$id,
    mw = desc$mw, logp = desc$logp, hbd = desc$hbd, hba = desc$hba,
    tpsa = desc$tpsa, rotatable_bonds = desc$rotatable_bonds,
    stringsAsFactors = FALSE
  )

  out$n_atoms <- desc$n_atoms

  out$fraction_csp3_approx <- desc$fraction_csp3
  out$aromatic_proportion_approx <- ifelse(is.na(desc$n_heavy) | desc$n_heavy == 0, NA_real_,
                                           desc$n_arom_atoms / desc$n_heavy)
  out$n_rings_approx <- vapply(cmp$smiles, .smiles_ring_count_approx, numeric(1), USE.NAMES = FALSE)
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
  ## Oprea 2000 property ranges (inclusive; the paper's drug-like ranges).
  ## n_rings_approx is the SMILES ring-closure-digit count (see
  ## .smiles_ring_count_approx()), not true CDK SSSR ring perception.
  out$oprea_pass <- with(out, hbd >= 0 & hbd <= 2 & hba >= 2 & hba <= 9 &
    rotatable_bonds >= 2 & rotatable_bonds <= 8 &
    n_rings_approx >= 1 & n_rings_approx <= 4)

  out$wlogp_proxy <- desc$wlogp_proxy
  out$wlogp_source <- "cdk_alogp_proxy"
  rk <- NULL
  if (wlogp_source == "rdkit" || rdkit_qc) {
    rk <- .rdkit_descriptors(cmp$smiles)
    if (is.null(rk)) {
      cli::cli_warn(c(
        "!" = "{if (wlogp_source == 'rdkit') '{.arg wlogp_source} = {.val rdkit}' else '{.arg rdkit_qc} = TRUE'} needed a Python {.pkg rdkit} that could not be resolved.",
        "i" = "wlogp_source: falling back to the CDK ALogP proxy. rdkit_qc: every {.code _rdkit} column is {.val NA}."
      ))
    }
  }
  if (wlogp_source == "rdkit") {
    if (!is.null(rk)) {
      use <- !is.na(rk$wlogp)
      out$wlogp_proxy[use] <- rk$wlogp[use]
      out$wlogp_source[use] <- "rdkit_wlogp"
      if (any(!use)) {
        proj <- .log_append(proj, step = "adme_local", id = cmp$id[!use & !failed],
                            message = "RDKit could not parse this structure for WLogP; kept the CDK ALogP proxy")
      }
    }
  }
  if (rdkit_qc) {
    out$hba_lipinski_rdkit <- if (!is.null(rk)) rk$hba_lipinski else NA_integer_
    out$hbd_lipinski_rdkit <- if (!is.null(rk)) rk$hbd_lipinski else NA_integer_
    out$tpsa_rdkit <- if (!is.null(rk)) rk$tpsa else NA_real_
    out$fraction_csp3_rdkit <- if (!is.null(rk)) rk$fraction_csp3 else NA_real_
    out$aromatic_proportion_rdkit <- if (!is.null(rk)) rk$aromatic_proportion else NA_real_
  }
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
                           rotatable_bonds = NA_integer_, n_atoms = NA_integer_, n_heavy = NA_integer_,
                           wlogp_proxy = NA_real_, amr = NA_real_,
                           fraction_csp3 = NA_real_, n_arom_atoms = NA_integer_)
  rows <- lapply(mols, function(m) {
    if (is.null(m)) return(empty_row)
    tryCatch({
      ## aromaticity must be perceived before the descriptors are read (Kekule input);
      ## if perception throws, the aromatic count would come from unprepared flags -> NA
      arom_ok <- tryCatch({ rcdk::do.aromaticity(m); TRUE }, error = function(e) FALSE)
      rcdk::convert.implicit.to.explicit(m)
      data.frame(
        mw  = .safe_mw(m),
        logp = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor", "XLogP"),
        hbd = as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor", "nHBDon")),
        hba = as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor", "nHBAcc")),
        tpsa = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor", "TopoPSA"),
        rotatable_bonds = .rotatable_bonds_veber(m),
        fraction_csp3 = .safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor", "Fsp3"),
        n_arom_atoms = if (arom_ok) as.integer(.safe_desc(m, "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor", "naAromAtom")) else NA_integer_,
        n_atoms = as.integer(tryCatch(rcdk::get.atom.count(m), error = function(e) NA_integer_)),
        n_heavy = tryCatch(sum(vapply(rcdk::get.atoms(m), function(a) rcdk::get.symbol(a) != "H", logical(1))),
                          error = function(e) NA_integer_),
        wlogp_proxy = tryCatch(as.numeric(rcdk::get.alogp(m)), error = function(e) NA_real_),
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

#' Average molecular weight (natural isotopic abundance), not the monoisotopic mass
#' @keywords internal
.safe_mw <- function(mol) {
  tryCatch(rcdk::get.natural.mass(mol), error = function(e) NA_real_)
}

#' Rotatable bonds as Veber et al. (2002) define them
#' @description Single non-ring bonds to non-terminal heavy atoms, **excluding
#'   amide C-N bonds** (CDK `RotatableBondsCountDescriptor` with
#'   `includeTerminals = FALSE, excludeAmides = TRUE`, set through rJava --
#'   `rcdk::eval.desc()` only exposes the defaults, which keep amides).
#' @param mol An rcdk molecule.
#' @return Integer, `NA_integer_` on failure.
#' @keywords internal
.rotatable_bonds_veber <- function(mol) {
  tryCatch({
    d <- rJava::.jnew("org/openscience/cdk/qsar/descriptors/molecular/RotatableBondsCountDescriptor")
    params <- rJava::.jarray(list(rJava::.jnew("java/lang/Boolean", FALSE), rJava::.jnew("java/lang/Boolean", TRUE)),
                             contents.class = "java/lang/Object")
    rJava::.jcall(d, "V", "setParameters", params)
    ac <- rJava::.jcast(mol, "org/openscience/cdk/interfaces/IAtomContainer")
    dv <- rJava::.jcall(d, "Lorg/openscience/cdk/qsar/DescriptorValue;", "calculate", ac)
    ## CDK can hand back a value carrying an exception instead of throwing: treat that as a failure
    if (!rJava::is.jnull(rJava::.jcall(dv, "Ljava/lang/Exception;", "getException"))) return(NA_integer_)
    res <- rJava::.jcall(dv, "Lorg/openscience/cdk/qsar/result/IDescriptorResult;", "getValue")
    as.integer(rJava::.jcall(rJava::.jcast(res, "org/openscience/cdk/qsar/result/IntegerResult"), "I", "intValue"))
  }, error = function(e) NA_integer_)
}

#' @keywords internal
.safe_desc <- function(mol, desc_name, col) {
  tryCatch({
    val <- rcdk::eval.desc(mol, desc_name, verbose = FALSE)
    out <- as.numeric(val[[col]])
    ## Inf/NaN from a descriptor is a failure, not a value the rules should compare
    if (length(out) != 1L || !is.finite(out)) NA_real_ else out
  }, error = function(e) NA_real_)
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

#' RDKit descriptors, via an optional Python `rdkit` (batched, one Python call)
#'
#' @description
#' Everything [adme_local()] can optionally pull from RDKit, computed in a
#' single Python call (not one per molecule -- reticulate's R/Python boundary
#' has real per-call overhead): `wlogp` (`Crippen.MolLogP()`, a direct
#' implementation of Wildman & Crippen's (1999) atom-contribution logP --
#' the quantity BOILED-Egg (Daina & Zoete 2016) is actually defined on),
#' `hba_lipinski`/`hbd_lipinski` (`rdMolDescriptors.CalcNumLipinskiHBA()`/
#' `CalcNumLipinskiHBD()`, the literal N+O count Lipinski et al. 2001
#' define), `tpsa` (`Descriptors.TPSA()`, Ertl et al. 2000, same method as
#' CDK's `tpsa` -- a cross-check, not a different definition), and
#' `fraction_csp3`/`aromatic_proportion` (RDKit's own hybridization/
#' aromaticity perception, a cross-check against CDK's
#' `fraction_csp3_approx`/`aromatic_proportion_approx`). Every failure mode
#' (no `reticulate`, no Python, no `rdkit` package, offline on the first
#' call that would need to fetch it, one bad SMILES) is caught: a single
#' unusable molecule is `NA` in every column; RDKit unavailable at all
#' returns `NULL` for the whole call -- callers fall back / warn, never error.
#'
#' @param smiles Character vector of SMILES.
#' @return A `data.frame(wlogp, hba_lipinski, hbd_lipinski, tpsa,
#'   fraction_csp3, aromatic_proportion)`, one row per `smiles` (`NA` for
#'   any molecule RDKit could not parse), or `NULL` if RDKit could not be
#'   made available at all.
#' @keywords internal
.rdkit_descriptors <- function(smiles) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  empty <- data.frame(wlogp = NA_real_, hba_lipinski = NA_integer_, hbd_lipinski = NA_integer_,
                      tpsa = NA_real_, fraction_csp3 = NA_real_, aromatic_proportion = NA_real_)
  ok <- tryCatch({
    reticulate::py_require("rdkit")
    isTRUE(reticulate::py_module_available("rdkit"))
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(NULL)
  py_fun <- tryCatch({
    reticulate::py_run_string(paste(
      "def patliR_rdkit_batch(smiles_list):",
      "    from rdkit import Chem",
      "    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors",
      "    out = []",
      "    for smi in smiles_list:",
      "        m = None if smi is None else Chem.MolFromSmiles(smi)",
      "        if m is None:",
      "            out.append(dict(wlogp=None, hba_lipinski=None, hbd_lipinski=None,",
      "                            tpsa=None, fraction_csp3=None, aromatic_proportion=None))",
      "            continue",
      "        n_heavy = m.GetNumHeavyAtoms()",
      "        n_arom = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())",
      "        out.append(dict(",
      "            wlogp=Crippen.MolLogP(m),",
      "            hba_lipinski=rdMolDescriptors.CalcNumLipinskiHBA(m),",
      "            hbd_lipinski=rdMolDescriptors.CalcNumLipinskiHBD(m),",
      "            tpsa=Descriptors.TPSA(m),",
      "            fraction_csp3=Descriptors.FractionCSP3(m),",
      "            aromatic_proportion=(n_arom / n_heavy if n_heavy else None),",
      "        ))",
      "    return out",
      sep = "\n"
    ))$patliR_rdkit_batch
  }, error = function(e) NULL)
  if (is.null(py_fun)) return(NULL)
  ## one Python call for every molecule, not one call per molecule
  smi_in <- ifelse(is.na(smiles) | !nzchar(smiles), NA_character_, smiles)
  res <- tryCatch(py_fun(as.list(smi_in)), error = function(e) NULL)
  if (is.null(res)) return(NULL)
  rows <- lapply(res, function(r) as.data.frame(lapply(r, function(x) if (is.null(x)) NA else x)))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
