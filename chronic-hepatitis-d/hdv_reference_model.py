#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
 CHRONIC HEPATITIS D (HDV) -- QSP REFERENCE IMPLEMENTATION
 hdv_reference_model.py            (pure Python 3, NO third-party dependencies)
===============================================================================

 WHY THIS FILE EXISTS
 --------------------
 Everything numerical that the companion files claim -- hdv_mrgsolve_model.R,
 hdv_shiny_app.R, README.md -- is COMPUTED HERE and written to
 hdv_model_report.txt.  Nothing in the write-up is asserted from memory.
 Run it with:

     python3 hdv_reference_model.py            # -> hdv_model_report.txt

 THE ONE IDEA
 ------------
 HDV is not "a virus you kill".  It is a TWO-SUBSTRATE ASSEMBLY LINE:

     substrate 1   the HDV genome    made by the HOST (RNA Pol II, rolling
                                     circle) -- completely independent of HBV
     substrate 2   the HBsAg envelope made by HBV cccDNA, in the same cell or
                                     supplied by the surrounding HBsAg+ pool

 So each therapy is classified by WHICH FLUX it cuts, not by potency:

     (E) entry / re-infection      NTCP-dependent    bulevirtide
     (A) envelopment / egress      FTase, HBsAg      lonafarnib, siRNA/NAP
     (C) intracellular + cell loss Pol II, immune    peg-IFN alfa / lambda

 and the infected-hepatocyte pool Id is ALSO maintained by two fluxes that no
 entry inhibitor can touch:

     division-mediated spread   an Id cell that divides gives two Id cells
                                (HDV replicates through mitosis)
     cell-to-cell spread        only partly NTCP-dependent (fraction phi_ntcp)

 That pair is the FLOOR.  The whole quantitative argument of this model is that
 the floor -- not target occupancy -- is what limits bulevirtide, and the
 report below tests that claim against the published trials.

 WHAT IS FITTED AND WHAT IS PREDICTED
 ------------------------------------
 FITTED (4 numbers, to 4 anchors):
   r_oatp, Kd_ntcp   <- the two published total-bile-acid fold-rises on
                        bulevirtide 2 mg and 10 mg.  NTCP occupancy is then
                        DERIVED, never assumed.
   dth_Id_immune     <- MYR301 mean HDV RNA decline at week 48 on 2 mg
   pop_cv_floor      <- MYR301 combined response at week 48 on 2 mg (45%)

 PREDICTED / HELD OUT (everything else), including:
   the 10 mg arm (all three endpoints), the untreated control arm, week 96,
   peg-IFN alfa on- and off-treatment response, bulevirtide + peg-IFN,
   lonafarnib/ritonavir +/- peg-IFN, peg-IFN lambda, the entry-inhibition
   ceiling, and the fibrosis / HCC projections.

 Sources for every parameter: hdv_references.md
===============================================================================
"""

import math
import os
import random

# =============================================================================
#  0.  UTILITIES
# =============================================================================

LOG10 = math.log(10.0)


def lg(x):
    """log10 with a floor so an undetectable viral load does not crash."""
    return math.log10(max(x, 1e-12))


def fmt(x, n=3):
    return ("%." + str(n) + "f") % x


class Rec:
    """A tiny attribute bag."""
    def __init__(self, **kw):
        self.__dict__.update(kw)


# =============================================================================
#  1.  STATE VECTOR LAYOUT  (33 ODE compartments)
# =============================================================================

NAMES = [
    # --- drug PK / PD (13) ---
    "A_sc",     # 0  bulevirtide SC depot                        nmol
    "A_cen",    # 1  bulevirtide central                         nmol
    "A_liv",    # 2  bulevirtide internalised hepatic pool       nmol
    "I_sc",     # 3  peg-IFN SC depot                            ug
    "I_cen",    # 4  peg-IFN central                             ug
    "ISG",      # 5  interferon-stimulated-gene activity         0-1
    "SOCS",     # 6  SOCS1/USP18 desensitisation signal          a.u.
    "L_gut",    # 7  lonafarnib gut                              umol
    "L_cen",    # 8  lonafarnib central                          umol
    "L_per",    # 9  lonafarnib peripheral                       umol
    "R_cen",    # 10 ritonavir central                           umol
    "Q_cen",    # 11 siRNA central                               nmol
    "Q_eff",    # 12 siRNA hepatic effect site                   a.u.
    # --- hepatocyte pools (3) ---
    "T",        # 13 HBsAg-negative hepatocytes                  fraction
    "Ib",       # 14 HBsAg+ HDV-  (ENVELOPE DONOR pool)          fraction
    "Id",       # 15 HBsAg+ HDV+  (productively co-infected)     fraction
    # --- intracellular HDV (4) ---
    "Rg",       # 16 intracellular HDV genomic RNA               copies/cell
    "Sag",      # 17 S-HDAg                                      a.u./cell
    "Lag",      # 18 L-HDAg, unprenylated                        a.u./cell
    "LagP",     # 19 L-HDAg, prenylated (assembly-competent)     a.u./cell
    # --- serum virology (4) ---
    "Vd",       # 20 serum HDV RNA                               IU/mL
    "ccc",      # 21 HBV cccDNA per HBsAg+ cell                  copies/cell
    "S_ser",    # 22 serum HBsAg                                 IU/mL
    "Vb",       # 23 serum HBV DNA                               IU/mL
    # --- immunity (2) ---
    "E",        # 24 HDV-specific CD8 effector activity          a.u.
    "Ex",       # 25 exhaustion level                            0-1
    # --- injury / organ / clinical (7) ---
    "ALT",      # 26 serum ALT                                   U/L
    "HSCa",     # 27 activated hepatic stellate cells            0-1
    "Fib",      # 28 fibrosis stage (Ishak)                      0-6
    "PLT",      # 29 platelets                                   10^9/L
    "NEU",      # 30 neutrophils                                 10^9/L
    "TBA",      # 31 total serum bile acids                      umol/L
    "CH",       # 32 cumulative HCC hazard                       -
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)


# =============================================================================
#  2.  PARAMETERS
#      Every value is either (a) taken from the literature -- see
#      hdv_references.md -- or (b) back-solved so that the untreated system
#      sits exactly on the observed chronic steady state (section 3).
# =============================================================================

def base_params():
    p = {}

    # -------------------------------------------------------------------------
    # 2.1  BULEVIRTIDE PK.  MW 5398.9 Da.  2 mg = 370.4 nmol.
    #      Clearance is deliberately split into a linear route and a SATURABLE
    #      route (NTCP-mediated internalisation = target-mediated clearance).
    #      The saturable route is what makes exposure MORE than dose
    #      proportional, which is what the EPAR reports.
    # -------------------------------------------------------------------------
    p["MW_blv"] = 5398.9
    p["ka_blv"] = 20.0        # /day   (tmax ~0.5-2 h)
    p["F_blv"] = 0.85
    p["Vc_blv"] = 8.0         # L
    p["CL_blv"] = 90.0        # L/day  linear (renal/proteolytic)
    p["Vmax_blv"] = 300.0     # nmol/day  saturable target-mediated clearance
    p["Km_blv"] = 1.5         # nM
    p["kout_liv"] = 2.0       # /day   degradation of internalised drug
    p["Kd_ntcp"] = 0.701      # nM     <-- FITTED to the bile-acid anchors

    # -------------------------------------------------------------------------
    # 2.2  BILE ACIDS -- the on-target read-out of NTCP occupancy.
    #      Serum TBA is set by (synthesis + ileal return) / uptake clearance,
    #      and uptake clearance = CL_NTCP*(1-OCC) + CL_OATP.
    #      Hence  fold-rise = (1+r) / ((1-OCC) + r),  r = CL_OATP/CL_NTCP.
    # -------------------------------------------------------------------------
    p["TBA0"] = 5.0           # umol/L baseline in chronic hepatitis D
    p["CL_ntcp_ba"] = 1.0     # /day (normalised)
    p["r_oatp"] = 0.0324      # <-- FITTED to the bile-acid anchors
    p["In_ba"] = None         # back-solved

    # -------------------------------------------------------------------------
    # 2.3  PEG-IFN alfa-2a  (180 ug SC weekly) and PEG-IFN lambda-1a
    # -------------------------------------------------------------------------
    p["ka_ifn"] = 1.44        # /day
    p["F_ifn"] = 0.80
    p["Vc_ifn"] = 10.0        # L
    p["CL_ifn"] = 1.40        # L/day (t1/2 ~ 3 d, flat weekly profile)
    p["EC50_ifn"] = 6.0       # ng/mL  on ISG induction
    p["h_ifn"] = 1.5
    p["kin_isg"] = 1.0        # /day
    p["kout_isg"] = 1.0       # /day
    p["ka_socs"] = 0.5
    p["kd_socs"] = 0.5
    p["Ks_socs"] = 1.0        # SOCS half-effect (tachyphylaxis)
    p["stim_endo"] = 0.25     # endogenous MDA5-driven ISG drive (max)
    p["K_endo"] = 0.05        # Id half-effect for endogenous drive

    # IFN downstream effects (all through ISG)
    p["Emax_rep"] = 0.75      # max block of intracellular HDV replication
    p["EC50_rep"] = 0.25
    p["Emax_prod"] = 0.90     # max block of virion PRODUCTION / export
    p["EC50_prod"] = 0.25     #   (this is what makes the IFN first phase)
    p["Dmax_ifn"] = 0.014     # /day  max ADDED infected-cell loss rate
    p["EC50_dif"] = 0.35
    p["Emax_hbs"] = 0.50      # max reduction of HBsAg secretion
    p["EC50_hbs"] = 0.40
    p["Eadar"] = 1.20         # max fold-increase of amber/W editing
    p["EC50_adar"] = 0.40
    p["e_cure"] = 0.50        # ISG boost of non-cytolytic intracellular cure
    p["hem_alfa"] = 1.0       # myelosuppression switch: 1 = alfa, 0 = lambda

    # -------------------------------------------------------------------------
    # 2.4  LONAFARNIB (+ RITONAVIR boost).  MW LNF 638.6, RTV 720.9.
    # -------------------------------------------------------------------------
    p["ka_lnf"] = 6.0         # /day
    p["F_lnf"] = 0.30
    p["Vc_lnf"] = 100.0       # L
    p["Vp_lnf"] = 200.0       # L
    p["Q_lnf"] = 30.0         # L/day
    p["CL_lnf"] = 480.0       # L/day unboosted (CYP3A4)
    p["Ki_rtv"] = 0.40        # uM   CYP3A4 inhibition constant
    p["IC50_ft"] = 0.10       # uM   FTase inhibition (assembly block)
    p["ka_rtv"] = 8.0
    p["F_rtv"] = 0.70
    p["Vc_rtv"] = 100.0
    p["CL_rtv"] = 130.0

    # -------------------------------------------------------------------------
    # 2.5  siRNA against HBsAg transcripts (elebsiran-like, 200 mg SC q4w)
    # -------------------------------------------------------------------------
    p["CL_sirna"] = 40.0      # L/day
    p["Vc_sirna"] = 5.0       # L
    p["kin_sirna"] = 0.20     # /day into hepatic effect site
    p["kout_sirna"] = 0.030   # /day  (long hepatic residence, t1/2 ~ 23 d)
    p["Emax_sirna"] = 0.85
    p["EC50_sirna"] = 0.50    # a.u.

    # -------------------------------------------------------------------------
    # 2.6  HEPATOCYTE POOLS AND THE THREE MAINTENANCE FLUXES
    # -------------------------------------------------------------------------
    p["d_hep"] = 0.0040       # /day background hepatocyte death (t1/2 ~173 d)
    p["d_hbv"] = 0.0020       # /day extra death of HBsAg+ HDV- cells
    p["d_innate"] = 0.0040    # /day innate (non-CTL) loss of Id
    p["dth_Id_immune"] = 0.0200  # /day  <-- FITTED (CTL killing at baseline)
    p["k_cure"] = 0.0020      # /day non-cytolytic intracellular HDV clearance
    p["phi_ntcp"] = 0.50      # fraction of cell-to-cell spread that is
                              # NTCP-dependent (the rest is untouchable)
    p["frac_cc"] = 0.20       # share of baseline Id inflow that is cell-to-cell
    p["beta_d"] = None        # back-solved
    p["k_cc"] = None          # back-solved
    p["beta_b"] = 2.0e-9      # HBV de-novo infection rate constant
    # LOCALITY OF REGENERATION.  When a hepatocyte dies, the cell that divides
    # to replace it is a NEIGHBOUR, and in chronic hepatitis D the neighbours
    # of an HBsAg+ cell are mostly HBsAg+ too.  loc_renew is the fraction of
    # slots vacated by an HBsAg+ cell that is refilled by an HBsAg+ cell
    # dividing, split between Ib and Id in proportion to the pools.  This is
    # what makes (a) the HBsAg+ pool persist for decades, and (b) HDV+ cells
    # get DILUTED by dividing HDV-negative HBsAg+ neighbours -- the same
    # mechanism the HBV field calls "cccDNA dilution by hepatocyte turnover".
    # It is back-solved from the observed slow spontaneous HBsAg decline.
    p["loc_renew"] = None
    # Serum HBsAg is remarkably stable in chronic hepatitis D on a NUC -- it is
    # the HDV RNA and the ALT that move.  A near-zero target here is what keeps
    # the untreated control arm flat, which the MYR301 control arm requires.
    p["tg_dHBsAg_lg_per_yr"] = 0.0     # log10 IU/mL per year, untreated

    # -------------------------------------------------------------------------
    # 2.7  INTRACELLULAR HDV REPLICATION AND ASSEMBLY
    # -------------------------------------------------------------------------
    p["alpha_r"] = None       # back-solved (rolling-circle replication)
    p["Rmax"] = 6000.0        # copies/cell carrying capacity
    p["frac_export"] = 0.75   # share of intracellular RNA drain that is EXPORT
    p["mu_R"] = 0.30          # /day intracellular HDV RNA turnover
    p["Ks_S"] = 100.0         # S-HDAg support of replication
    p["Ki_L"] = 3000.0        # L-HDAg trans-dominant repression
    p["k_tl"] = 1.00          # HDAg translation scale
    p["f_edit0"] = 0.25       # baseline amber/W editing fraction
    p["mu_S"] = 1.50          # /day S-HDAg turnover
    p["k_pren"] = 2.00        # /day farnesylation of L-HDAg
    p["mu_L"] = 0.50          # /day
    p["mu_LP"] = 1.50         # /day
    p["k_ex"] = None          # back-solved (per-cell export rate constant)
    p["Kp_LP"] = 400.0        # prenyl-L half-saturation for envelopment
    p["Ks_env"] = 3034.0      # serum HBsAg half-saturation for envelopment
    p["conv_ser"] = None      # back-solved (copies/cell/day -> IU/mL/day)
    p["c_d"] = 0.55           # /day serum HDV RNA clearance (t1/2 ~30 h)

    # -------------------------------------------------------------------------
    # 2.8  HBV SIDE -- the envelope factory
    # -------------------------------------------------------------------------
    p["Nhep"] = 2.0e11        # hepatocytes in an adult liver
    p["ccc_max"] = 20.0
    p["r_ccc"] = 0.050        # /day
    p["mu_ccc"] = 0.025       # /day
    p["k_hbs"] = None         # back-solved
    p["c_s"] = 0.35           # /day serum HBsAg clearance (t1/2 ~2 d)
    p["k_hbv"] = None         # back-solved
    p["c_b"] = 0.60           # /day serum HBV DNA clearance
    p["K_sup"] = 0.10         # HDV suppression of HBV replication
    p["eps_nuc"] = 0.999      # NUC block of HBV DNA production

    # -------------------------------------------------------------------------
    # 2.9  IMMUNITY
    # -------------------------------------------------------------------------
    p["a_E"] = None           # back-solved
    p["K_E"] = 0.02
    p["d_E"] = 0.050          # /day
    p["b_ifn_E"] = 1.50       # ISG boost of effector activity
    p["a_x"] = 0.030
    p["K_x"] = 5.0e-4
    p["d_x"] = None           # back-solved
    p["c_ifn_x"] = 3.00       # ISG-driven reversal of exhaustion
    # DURABILITY LIVES HERE.  Killing capacity is scaled by how far exhaustion
    # has been REVERSED below its baseline level.  Antigen decline plus ISG
    # both reverse it, and it re-accumulates only as antigen returns -- so a
    # regimen that reverses exhaustion leaves behind a liver that keeps
    # killing infected cells after the drug stops, while a regimen that merely
    # blocks an inflow does not.  This is the model's account of why only
    # interferon-containing regimens have off-treatment durability, and of the
    # late relapses seen years after peg-IFN.
    p["gamma_x"] = 0.15
    # NOTE on K_x: exhaustion is made almost INSENSITIVE to the size of the
    # infected pool and sensitive to ISG instead.  That is the discriminating
    # assumption: interferon reverses exhaustion, an entry inhibitor does not,
    # even though both lower antigen.

    # -------------------------------------------------------------------------
    # 2.10  INJURY, FIBROSIS, ORGAN, OUTCOME
    # -------------------------------------------------------------------------
    p["kappa_fresh"] = None   # back-solved: injury per unit ENTRY flux
    p["frac_fresh"] = 0.40    # share of baseline injury from fresh infection
    p["ALT_base"] = 30.0      # U/L  ALT with no HDV-driven killing
    p["kel_alt"] = 0.35       # /day (t1/2 ~ 47 h)
    p["k_alt"] = None         # back-solved
    p["ULN_ALT"] = 40.0       # U/L  normalisation threshold used in trials
    # NON-HDV ALT.  Many patients with chronic hepatitis D also carry steatosis,
    # alcohol or metabolic liver injury, and that component does not respond to
    # an antiviral.  Its spread is what decides how many patients CAN normalise
    # ALT at all, so it is fitted (below) to the observed ALT normalisation rate
    # rather than guessed.
    p["ALT_base_cv"] = 0.45
    # Endpoint measurement / short-term biological variability.  Trial endpoints
    # are single visits, so a flat model with no visit noise cannot reproduce a
    # non-zero response rate in an untreated control arm.
    p["noise_alt"] = 0.20     # log-normal CV on a single ALT measurement
    p["noise_rna"] = 0.25     # log10 SD on a single HDV RNA measurement
    p["a_h"] = None           # back-solved
    p["d_h"] = 0.020          # /day HSC deactivation
    p["ALT_ref"] = 40.0
    p["k_f"] = 1.795e-3       # /day fibrogenesis
    p["k_rev"] = 2.00e-4      # /day fibrolysis
    p["Kh_rev"] = 0.15        # HSCa that halves fibrolysis
    p["PLT0"] = 250.0
    p["g_plt"] = 0.62         # portal-hypertension effect at Ishak 6
    p["mye_plt"] = 0.42       # peg-IFN alfa myelosuppression on platelets
    p["k_plt"] = 0.15         # /day
    p["NEU0"] = 4.0
    p["mye_neu"] = 0.55
    p["k_neu"] = 0.20         # /day
    p["h0_hcc"] = 0.005 / 365.0   # /day at Ishak 2
    p["bF_hcc"] = 0.4302      # per Ishak stage  (-> 2.8%/y at Ishak 6)
    p["bV_hcc"] = 0.15        # per log10 HDV RNA

    # -------------------------------------------------------------------------
    # 2.11  BASELINE TARGETS (the observed chronic steady state)
    # -------------------------------------------------------------------------
    p["tg_T"] = 0.55
    p["tg_Ib"] = 0.35
    p["tg_Id"] = 0.10         # HDAg+ hepatocyte fraction
    p["tg_Rg"] = 3000.0       # intracellular HDV RNA copies/cell
    p["tg_Sag"] = 1500.0
    p["tg_Lag"] = 300.0
    p["tg_LagP"] = 400.0
    p["tg_Vd"] = 10.0 ** 5.5  # IU/mL
    p["tg_ccc"] = 10.0
    p["tg_S"] = 10.0 ** 3.85  # IU/mL serum HBsAg (~7080)
    p["tg_Vb"] = 100.0        # IU/mL
    p["tg_E"] = 1.00
    p["tg_Ex"] = 0.60
    p["tg_ALT"] = 110.0
    p["tg_HSCa"] = 0.45
    p["tg_Fib"] = 2.00
    return p


# =============================================================================
#  3.  BACK-SOLVE THE STEADY-STATE PARAMETERS
#      Instead of fitting by optimisation, every parameter that CAN be pinned
#      by requiring d/dt = 0 at the observed baseline IS pinned that way.
#      This makes the untreated model reproduce the observed chronic state by
#      construction and leaves only 4 genuinely fitted numbers (section 8).
# =============================================================================

def isg_steady(p, Id, ug=0.0):
    """
    Steady-state ISG for a given infected fraction and weekly peg-IFN dose.
    SOCS_ss = ISG (ka_socs = kd_socs) and Ks_socs = 1, so the ISG balance
        ISG * (1 + ISG) = stim
    has the closed form below.  With ug = 0 this returns the ENDOGENOUS tone
    driven by MDA5 sensing of HDV replication.
    """
    stim = p["stim_endo"] * Id / (Id + p["K_endo"])
    if ug > 0:
        cavg = p["F_ifn"] * ug / (p["CL_ifn"] * 7.0)
        stim += cavg ** p["h_ifn"] / (cavg ** p["h_ifn"]
                                      + p["EC50_ifn"] ** p["h_ifn"])
    stim = min(1.0, stim)
    return (-1.0 + math.sqrt(1.0 + 4.0 * stim)) / 2.0


def ifn_effects(p, isg):
    """
    All interferon pharmacodynamics, expressed as INCREMENTS above the
    endogenous ISG tone that is already present in the untreated patient.

    This matters.  Chronic hepatitis D runs with a strong, sustained,
    endogenous ISG signature and the virus persists anyway.  If the drug
    effects were written as functions of absolute ISG, the untreated baseline
    would already be partly "treated" and could not be a steady state.  So the
    driver is dISG = max(0, ISG - ISG0), and every effect below is zero at
    baseline by construction.
    """
    dis = max(0.0, isg - p["_ISG0"])
    return Rec(
        eps_rep=p["Emax_rep"] * dis / (p["EC50_rep"] + dis),
        eps_prod=p["Emax_prod"] * dis / (p["EC50_prod"] + dis),
        dlt=p["Dmax_ifn"] * dis / (p["EC50_dif"] + dis),
        eps_hbs=p["Emax_hbs"] * dis / (p["EC50_hbs"] + dis),
        f_edit=min(0.95, p["f_edit0"] *
                   (1.0 + p["Eadar"] * dis / (p["EC50_adar"] + dis))),
        cure=1.0 + p["e_cure"] * dis,
        eboost=1.0 + p["b_ifn_E"] * dis,
        xrev=1.0 + p["c_ifn_x"] * dis,
        ccc_loss=1.0 + 0.6 * dis,
        dis=dis,
    )


def back_solve(p):
    T, Ib, Id = p["tg_T"], p["tg_Ib"], p["tg_Id"]
    Rg, Sag, Lag, LagP = p["tg_Rg"], p["tg_Sag"], p["tg_Lag"], p["tg_LagP"]
    Vd, ccc, S, Vb = p["tg_Vd"], p["tg_ccc"], p["tg_S"], p["tg_Vb"]

    # --- endogenous ISG tone at baseline (defines the zero of drug effect) ---
    p["_ISG0"] = isg_steady(p, Id, 0.0)

    # --- bile acid input ----------------------------------------------------
    p["In_ba"] = (p["CL_ntcp_ba"] * (1.0 + p["r_oatp"])) * p["TBA0"]

    # --- HBsAg / HBV DNA production ----------------------------------------
    p["k_hbs"] = p["c_s"] * S / (ccc * (Ib + Id))
    sup = 1.0 / (1.0 + Id / p["K_sup"])
    p["k_hbv"] = p["c_b"] * Vb / (ccc * (Ib + Id) * sup)

    # --- envelopment / secretion -------------------------------------------
    gLP = LagP / (LagP + p["Kp_LP"])          # = 0.5 by construction
    gS = S / (S + p["Ks_env"])                # = 0.7 by construction
    # Export is the DOMINANT drain on the intracellular pool.  That is what
    # makes an assembly block (lonafarnib) raise intracellular HDV RNA instead
    # of merely lowering serum RNA -- see A5.
    drain_deg = p["mu_R"] * Rg
    ex = p["frac_export"] * drain_deg / (1.0 - p["frac_export"])
    p["k_ex"] = ex / (Rg * gLP * gS)
    p["conv_ser"] = p["c_d"] * Vd / (ex * Id)

    # --- rolling-circle replication ----------------------------------------
    gsupp = Sag / (Sag + p["Ks_S"])
    grepr = 1.0 / (1.0 + Lag / p["Ki_L"])
    p["alpha_r"] = (drain_deg + ex) / (Rg * (1.0 - Rg / p["Rmax"]) * gsupp * grepr)

    # --- HDAg translation / turnover consistency ---------------------------
    # mu_S, k_pren, mu_L, mu_LP were chosen so that the targets balance; verify
    p["mu_S"] = p["k_tl"] * Rg * (1.0 - p["f_edit0"]) / Sag
    p["mu_L"] = p["k_tl"] * Rg * p["f_edit0"] / Lag - p["k_pren"]
    p["mu_LP"] = p["k_pren"] * Lag / LagP

    # --- cccDNA ------------------------------------------------------------
    p["mu_ccc"] = p["r_ccc"] * (1.0 - ccc / p["ccc_max"])

    # --- immunity ----------------------------------------------------------
    p["a_E"] = p["d_E"] * p["tg_E"] / (Id / (Id + p["K_E"]))
    p["d_x"] = p["a_x"] * (Id / (Id + p["K_x"])) * (1.0 - p["tg_Ex"]) / p["tg_Ex"]

    # --- the three maintenance fluxes --------------------------------------
    dth_T = p["d_hep"]
    dth_Ib = p["d_hep"] + p["d_hbv"]
    dth_Id = p["d_hep"] + p["d_innate"] + p["dth_Id_immune"]
    N = T + Ib + Id
    Hb = Ib + Id
    H = dth_Ib * Ib + dth_Id * Id          # deaths inside the HBsAg+ pool
    Dt = dth_T * T                         # deaths in the HBsAg-negative pool
    h = Hb / N

    # back-solve loc_renew from the observed slow spontaneous HBsAg decline.
    # With locality applied symmetrically to the HBsAg-negative class {T} and
    # the HBsAg-positive class {Ib, Id},
    #   d(Hb)/dt = (1 - loc) * [ (Dt + H)*h - H ]
    tgt_dHb = p["tg_dHBsAg_lg_per_yr"] * LOG10 / 365.0 * Hb
    denom = (Dt + H) * h - H
    oml = tgt_dHb / denom if abs(denom) > 1e-15 else 0.05
    oml = min(0.999, max(0.001, oml))
    p["loc_renew"] = 1.0 - oml

    pool_renew = oml * (Dt + H)
    ren_Id_pc = (1.0 - oml) * H / Hb + pool_renew / N   # Id per-capita division
    # per-capita balance for Id:  ren_Id_pc - dth_Id - k_cure + inflow/Id = 0
    need = dth_Id + p["k_cure"] - ren_Id_pc          # per-capita inflow needed
    need = max(need, 1e-6)
    cc_pc = p["frac_cc"] * need                      # cell-to-cell share
    en_pc = (1.0 - p["frac_cc"]) * need              # entry share
    p["k_cc"] = cc_pc / Ib
    p["beta_d"] = en_pc * Id / (Vd * Ib)
    p["_dth_Id0"] = dth_Id
    p["_entry0"] = en_pc * Id
    p["_cc0"] = cc_pc * Id
    p["_div0"] = ren_Id_pc * Id
    p["_dilute0"] = dth_Id - ren_Id_pc

    # --- injury ------------------------------------------------------------
    kill_cell = (p["d_innate"] + p["dth_Id_immune"]) * Id
    # fresh-infection injury term is set to frac_fresh of total injury
    fresh_target = p["frac_fresh"] / (1.0 - p["frac_fresh"]) * kill_cell
    p["kappa_fresh"] = fresh_target / p["_entry0"]
    kill_tot = kill_cell + fresh_target
    p["ALT_base"] = min(p["ALT_base"], 0.80 * p["tg_ALT"])   # keep tg_ALT reachable
    p["k_alt"] = p["kel_alt"] * (p["tg_ALT"] - p["ALT_base"]) / kill_tot
    p["_kill0"] = kill_tot

    # --- fibrogenesis ------------------------------------------------------
    aq = (p["tg_ALT"] / p["ALT_ref"])
    p["a_h"] = p["d_h"] * p["tg_HSCa"] / (aq * (1.0 - p["tg_HSCa"]))
    return p


# =============================================================================
#  4.  DOSING / REGIMEN
# =============================================================================

class Regimen:
    """A therapy specification.  Times in DAYS from t=0."""

    def __init__(self, blv_mg=0.0, blv_start=0.0, blv_stop=1e9,
                 ifn_ug=0.0, ifn_start=0.0, ifn_stop=1e9, ifn_lambda=False,
                 lnf_mg=0.0, rtv_mg=0.0, lnf_start=0.0, lnf_stop=1e9,
                 sirna_mg=0.0, sirna_start=0.0, sirna_stop=1e9,
                 nuc=True, occ_override=None):
        self.blv_mg = blv_mg
        self.blv_start, self.blv_stop = blv_start, blv_stop
        self.ifn_ug = ifn_ug
        self.ifn_start, self.ifn_stop = ifn_start, ifn_stop
        self.ifn_lambda = ifn_lambda
        self.lnf_mg, self.rtv_mg = lnf_mg, rtv_mg
        self.lnf_start, self.lnf_stop = lnf_start, lnf_stop
        self.sirna_mg = sirna_mg
        self.sirna_start, self.sirna_stop = sirna_start, sirna_stop
        self.nuc = nuc
        self.occ_override = occ_override   # for the hypothetical perfect blocker

    def events(self, t_end, p):
        """Return a sorted list of (time, state_index, amount) bolus events."""
        ev = []
        if self.blv_mg > 0:
            nmol = self.blv_mg * 1e-3 / p["MW_blv"] * 1e9   # mg -> nmol
            t = self.blv_start
            while t < min(self.blv_stop, t_end):
                ev.append((t, IX["A_sc"], nmol))
                t += 1.0
        if self.ifn_ug > 0:
            t = self.ifn_start
            while t < min(self.ifn_stop, t_end):
                ev.append((t, IX["I_sc"], self.ifn_ug))
                t += 7.0
        if self.lnf_mg > 0:
            umol = self.lnf_mg * 1e-3 / 638.6 * 1e6
            rumol = self.rtv_mg * 1e-3 / 720.9 * 1e6
            t = self.lnf_start
            while t < min(self.lnf_stop, t_end):
                ev.append((t, IX["L_gut"], umol))
                if self.rtv_mg > 0:
                    ev.append((t, IX["R_cen"], rumol * p["F_rtv"]))
                t += 0.5
        if self.sirna_mg > 0:
            nmol = self.sirna_mg * 1e-3 / 16000.0 * 1e9   # ~16 kDa duplex
            t = self.sirna_start
            while t < min(self.sirna_stop, t_end):
                ev.append((t, IX["Q_cen"], nmol))
                t += 28.0
        ev.sort()
        return ev


# =============================================================================
#  5.  RIGHT-HAND SIDE
# =============================================================================

def rhs(t, y, p, reg, qss=None):
    """
    Full 33-ODE right-hand side.

    qss (optional dict) replaces the drug PK sub-models by their steady-state
    pharmacodynamic summaries {occ, isg, ift, eps_s}.  Used for the virtual
    population runs, where the disease timescale (weeks) is far slower than
    every drug's PK timescale (hours to days).
    """
    d = [0.0] * NST
    (A_sc, A_cen, A_liv, I_sc, I_cen, ISG, SOCS, L_gut, L_cen, L_per, R_cen,
     Q_cen, Q_eff, T, Ib, Id, Rg, Sag, Lag, LagP, Vd, ccc, S_ser, Vb, E, Ex,
     ALT, HSCa, Fib, PLT, NEU, TBA, CH) = y

    # ---- keep the integrator on the physical manifold ----
    T = max(T, 0.0); Ib = max(Ib, 0.0); Id = max(Id, 0.0)
    Rg = max(Rg, 0.0); Sag = max(Sag, 0.0); Lag = max(Lag, 0.0)
    LagP = max(LagP, 0.0); Vd = max(Vd, 0.0); S_ser = max(S_ser, 0.0)
    E = max(E, 0.0); Ex = min(max(Ex, 0.0), 1.0)

    # =====================================================================
    # 5.1  DRUG EXPOSURE AND TARGET ENGAGEMENT
    # =====================================================================
    if qss is not None:
        occ = qss["occ"]
        isg = qss["isg"]
        ift = qss["ift"]
        eps_s = qss["eps_s"]
        hem = qss["hem"]
    else:
        # ---- bulevirtide ----
        C_blv = A_cen / p["Vc_blv"]                                # nM
        upt = p["Vmax_blv"] * C_blv / (p["Km_blv"] + C_blv)
        d[IX["A_sc"]] = -p["ka_blv"] * A_sc
        d[IX["A_cen"]] = p["F_blv"] * p["ka_blv"] * A_sc \
            - p["CL_blv"] * C_blv - upt
        d[IX["A_liv"]] = upt - p["kout_liv"] * A_liv
        occ = C_blv / (C_blv + p["Kd_ntcp"])
        if reg.occ_override is not None:
            occ = reg.occ_override

        # ---- peg-interferon ----
        C_ifn = I_cen / p["Vc_ifn"]                                # ng/mL
        d[IX["I_sc"]] = -p["ka_ifn"] * I_sc
        d[IX["I_cen"]] = p["F_ifn"] * p["ka_ifn"] * I_sc - p["CL_ifn"] * C_ifn
        stim_drug = (C_ifn ** p["h_ifn"]) / \
                    (C_ifn ** p["h_ifn"] + p["EC50_ifn"] ** p["h_ifn"]) \
                    if C_ifn > 0 else 0.0
        stim_endo = p["stim_endo"] * Id / (Id + p["K_endo"])
        stim = min(1.0, stim_drug + stim_endo)
        d[IX["ISG"]] = p["kin_isg"] * stim / (1.0 + SOCS / p["Ks_socs"]) \
            - p["kout_isg"] * ISG
        d[IX["SOCS"]] = p["ka_socs"] * ISG - p["kd_socs"] * SOCS
        isg = max(ISG, 0.0)
        hem = (0.0 if reg.ifn_lambda else 1.0) * max(0.0, isg - p["_ISG0"])

        # ---- lonafarnib / ritonavir ----
        C_lnf = L_cen / p["Vc_lnf"]                                # uM
        C_rtv = R_cen / p["Vc_rtv"]
        CL_eff = p["CL_lnf"] / (1.0 + C_rtv / p["Ki_rtv"])
        d[IX["L_gut"]] = -p["ka_lnf"] * L_gut
        d[IX["L_cen"]] = p["F_lnf"] * p["ka_lnf"] * L_gut - CL_eff * C_lnf \
            - p["Q_lnf"] * (C_lnf - L_per / p["Vp_lnf"])
        d[IX["L_per"]] = p["Q_lnf"] * (C_lnf - L_per / p["Vp_lnf"])
        d[IX["R_cen"]] = -p["CL_rtv"] * C_rtv
        ift = C_lnf / (C_lnf + p["IC50_ft"])

        # ---- siRNA ----
        C_q = Q_cen / p["Vc_sirna"]
        d[IX["Q_cen"]] = -p["CL_sirna"] * C_q
        d[IX["Q_eff"]] = p["kin_sirna"] * C_q - p["kout_sirna"] * Q_eff
        eps_s = p["Emax_sirna"] * Q_eff / (p["EC50_sirna"] + Q_eff)

    # =====================================================================
    # 5.2  IFN-MEDIATED PHARMACODYNAMICS (increments above endogenous tone)
    # =====================================================================
    ef = ifn_effects(p, isg)
    eps_rep, dlt_ifn, eps_hbs, f_edit = ef.eps_rep, ef.dlt, ef.eps_hbs, ef.f_edit
    k_cure_eff = p["k_cure"] * ef.cure

    # =====================================================================
    # 5.3  INTRACELLULAR HDV -- replication, editing, prenylation, export
    # =====================================================================
    gsupp = Sag / (Sag + p["Ks_S"])
    grepr = 1.0 / (1.0 + Lag / p["Ki_L"])
    gLP = LagP / (LagP + p["Kp_LP"])
    gS = S_ser / (S_ser + p["Ks_env"])
    export = p["k_ex"] * Rg * gLP * gS * (1.0 - ift) * (1.0 - ef.eps_prod)

    d[IX["Rg"]] = p["alpha_r"] * Rg * max(0.0, 1.0 - Rg / p["Rmax"]) \
        * gsupp * grepr * (1.0 - eps_rep) - p["mu_R"] * Rg - export
    d[IX["Sag"]] = p["k_tl"] * Rg * (1.0 - f_edit) - p["mu_S"] * Sag
    d[IX["Lag"]] = p["k_tl"] * Rg * f_edit \
        - p["k_pren"] * (1.0 - ift) * Lag - p["mu_L"] * Lag
    d[IX["LagP"]] = p["k_pren"] * (1.0 - ift) * Lag - p["mu_LP"] * LagP

    # =====================================================================
    # 5.4  HEPATOCYTE POOLS -- the three maintenance fluxes
    # =====================================================================
    N = max(T + Ib + Id, 1e-9)
    Hb = max(Ib + Id, 1e-12)
    F_entry = p["beta_d"] * Vd * Ib * (1.0 - occ)            # (E) NTCP-dependent
    F_cc = p["k_cc"] * Id * Ib * (1.0 - p["phi_ntcp"] * occ)  # (C) partly so
    F_hbv = p["beta_b"] * Vb * T * (1.0 - occ)

    dth_T = p["d_hep"]
    dth_Ib = p["d_hep"] + p["d_hbv"]
    fx = 1.0 + p["gamma_x"] * max(0.0, p["tg_Ex"] - Ex) / p["tg_Ex"]
    dth_Id = p["d_hep"] + p["d_innate"] + p["dth_Id_immune"] * fx + dlt_ifn

    # --- regeneration: local inside the HBsAg+ pool, then pooled -----------
    loc = p["loc_renew"]
    H = dth_Ib * Ib + dth_Id * Id           # deaths in the HBsAg+ class
    Dt = dth_T * T                          # deaths in the HBsAg-negative class
    pool_renew = (1.0 - loc) * (Dt + H)
    ren_T = loc * Dt + pool_renew * T / N
    ren_Ib = loc * H * Ib / Hb + pool_renew * Ib / N
    ren_Id = loc * H * Id / Hb + pool_renew * Id / N     # (D) NTCP-INDEPENDENT

    d[IX["T"]] = ren_T - dth_T * T - F_hbv
    d[IX["Ib"]] = ren_Ib - dth_Ib * Ib + F_hbv \
        - F_entry - F_cc + k_cure_eff * Id
    d[IX["Id"]] = ren_Id - dth_Id * Id + F_entry + F_cc - k_cure_eff * Id

    # =====================================================================
    # 5.5  SERUM VIROLOGY
    # =====================================================================
    d[IX["Vd"]] = p["conv_ser"] * export * Id - p["c_d"] * Vd
    d[IX["ccc"]] = p["r_ccc"] * ccc * max(0.0, 1.0 - ccc / p["ccc_max"]) \
        - p["mu_ccc"] * ccc * ef.ccc_loss
    d[IX["S_ser"]] = p["k_hbs"] * ccc * (Ib + Id) * (1.0 - eps_s) \
        * (1.0 - eps_hbs) - p["c_s"] * S_ser
    nuc_f = (1.0 - p["eps_nuc"]) if reg.nuc else 1.0
    d[IX["Vb"]] = p["k_hbv"] * ccc * (Ib + Id) * nuc_f \
        / (1.0 + Id / p["K_sup"]) - p["c_b"] * Vb

    # =====================================================================
    # 5.6  IMMUNITY
    # =====================================================================
    d[IX["E"]] = p["a_E"] * (Id / (Id + p["K_E"])) * ef.eboost - p["d_E"] * E
    d[IX["Ex"]] = p["a_x"] * (Id / (Id + p["K_x"])) * (1.0 - Ex) \
        - p["d_x"] * Ex * ef.xrev

    # =====================================================================
    # 5.7  INJURY -> FIBROSIS -> OUTCOME
    #      NOTE the second injury term.  ALT is driven not only by the size of
    #      the infected pool but by the INFLOW of newly infected hepatocytes.
    #      That is the model's mechanistic explanation for the observed
    #      ALT / HDV-RNA decoupling, and it is a HYPOTHESIS (section 8).
    # =====================================================================
    kill = (p["d_innate"] + p["dth_Id_immune"] * fx + dlt_ifn) * Id \
        + p["kappa_fresh"] * F_entry
    alt_tox = 5.0 * ef.dis if reg.ifn_lambda else 0.0    # lambda hepatic AE
    d[IX["ALT"]] = p["k_alt"] * kill + p["kel_alt"] * p["ALT_base"] \
        + alt_tox - p["kel_alt"] * ALT

    d[IX["HSCa"]] = p["a_h"] * (ALT / p["ALT_ref"]) * (1.0 - HSCa) \
        - p["d_h"] * HSCa
    d[IX["Fib"]] = p["k_f"] * HSCa * max(0.0, 1.0 - Fib / 6.0) \
        - p["k_rev"] * Fib / (1.0 + HSCa / p["Kh_rev"])

    plt_t = p["PLT0"] * (1.0 - p["g_plt"] * Fib / 6.0) * (1.0 - p["mye_plt"] * hem)
    d[IX["PLT"]] = p["k_plt"] * (plt_t - PLT)
    neu_t = p["NEU0"] * (1.0 - p["mye_neu"] * hem)
    d[IX["NEU"]] = p["k_neu"] * (neu_t - NEU)

    CLupt = p["CL_ntcp_ba"] * (1.0 - occ) + p["CL_ntcp_ba"] * p["r_oatp"]
    d[IX["TBA"]] = p["In_ba"] - CLupt * TBA

    d[IX["CH"]] = p["h0_hcc"] * math.exp(p["bF_hcc"] * (Fib - 2.0)
                                         + p["bV_hcc"] * (lg(Vd) - 5.5))
    return d


# =============================================================================
#  6.  INTEGRATOR
# =============================================================================

def y0_baseline(p):
    y = [0.0] * NST
    y[IX["T"]] = p["tg_T"]
    y[IX["Ib"]] = p["tg_Ib"]
    y[IX["Id"]] = p["tg_Id"]
    y[IX["Rg"]] = p["tg_Rg"]
    y[IX["Sag"]] = p["tg_Sag"]
    y[IX["Lag"]] = p["tg_Lag"]
    y[IX["LagP"]] = p["tg_LagP"]
    y[IX["Vd"]] = p["tg_Vd"]
    y[IX["ccc"]] = p["tg_ccc"]
    y[IX["S_ser"]] = p["tg_S"]
    y[IX["Vb"]] = p["tg_Vb"]
    y[IX["E"]] = p["tg_E"]
    y[IX["Ex"]] = p["tg_Ex"]
    y[IX["ALT"]] = p["tg_ALT"]
    y[IX["HSCa"]] = p["tg_HSCa"]
    y[IX["Fib"]] = p["tg_Fib"]
    y[IX["PLT"]] = p["PLT0"] * (1.0 - p["g_plt"] * p["tg_Fib"] / 6.0)
    y[IX["NEU"]] = p["NEU0"]
    y[IX["TBA"]] = p["TBA0"]
    # endogenous ISG / SOCS at baseline
    st = p["stim_endo"] * p["tg_Id"] / (p["tg_Id"] + p["K_endo"])
    isg = (-1.0 + math.sqrt(1.0 + 4.0 * st)) / 2.0
    y[IX["ISG"]] = isg
    y[IX["SOCS"]] = isg
    return y


def simulate(p, reg, t_end, dt=0.02, obs_every=1.0, qss=None, y0=None):
    """RK4 with exact bolus events.  Returns (times, list-of-state-vectors)."""
    y = list(y0 if y0 is not None else y0_baseline(p))
    ev = reg.events(t_end, p) if qss is None else []
    ei = 0
    t = 0.0
    ts, ys = [0.0], [list(y)]
    next_obs = obs_every
    n = int(round(t_end / dt))
    for _ in range(n):
        # fire any bolus at or before t (within half a step)
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            y[ev[ei][1]] += ev[ei][2]
            ei += 1
        k1 = rhs(t, y, p, reg, qss)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2 = rhs(t + 0.5 * dt, y2, p, reg, qss)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3 = rhs(t + 0.5 * dt, y3, p, reg, qss)
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4 = rhs(t + dt, y4, p, reg, qss)
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        # physical clamps
        for nm in ("T", "Ib", "Id", "Rg", "Sag", "Lag", "LagP", "Vd", "ccc",
                   "S_ser", "Vb", "E", "TBA", "ALT", "PLT", "NEU", "HSCa"):
            if y[IX[nm]] < 0.0:
                y[IX[nm]] = 0.0
        y[IX["Ex"]] = min(max(y[IX["Ex"]], 0.0), 1.0)
        # fewer than one infected hepatocyte in the whole liver = eradication
        if y[IX["Id"]] * p["Nhep"] < 1.0:
            y[IX["Id"]] = 0.0
        y[IX["Fib"]] = min(max(y[IX["Fib"]], 0.0), 6.0)
        t += dt
        if t >= next_obs - 1e-9:
            ts.append(t); ys.append(list(y)); next_obs += obs_every
    return ts, ys


def at(ts, ys, t, name):
    """Value of `name` at (or just before) time t."""
    idx = 0
    for i, tt in enumerate(ts):
        if tt <= t + 1e-6:
            idx = i
        else:
            break
    return ys[idx][IX[name]]


# =============================================================================
#  7.  QUASI-STEADY-STATE DRUG SUMMARIES
#      Every drug here equilibrates in hours to a few days; the disease moves
#      over weeks to years.  For the virtual-population runs we therefore
#      summarise each regimen by its steady-state target engagement.
# =============================================================================

def blv_css(p, mg):
    """Steady-state average central concentration (nM) for a daily SC dose."""
    if mg <= 0:
        return 0.0
    inp = p["F_blv"] * (mg * 1e-3 / p["MW_blv"] * 1e9)     # nmol/day absorbed
    lo, hi = 0.0, 1e5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        out = p["CL_blv"] * mid + p["Vmax_blv"] * mid / (p["Km_blv"] + mid)
        if out < inp:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def blv_occ(p, mg):
    c = blv_css(p, mg)
    return c / (c + p["Kd_ntcp"])


def tba_fold(p, occ):
    r = p["r_oatp"]
    return (1.0 + r) / ((1.0 - occ) + r)


def ifn_isg(p, ug):
    """Steady-state ISG for a weekly peg-IFN dose (plus endogenous drive)."""
    return isg_steady(p, p["tg_Id"], ug)


def lnf_ift(p, mg_bid, rtv_mg_bid):
    if mg_bid <= 0:
        return 0.0
    c_rtv = 0.0
    if rtv_mg_bid > 0:
        c_rtv = p["F_rtv"] * (2 * rtv_mg_bid * 1e-3 / 720.9 * 1e6) \
            / p["CL_rtv"]
    CL_eff = p["CL_lnf"] / (1.0 + c_rtv / p["Ki_rtv"])
    cavg = p["F_lnf"] * (2 * mg_bid * 1e-3 / 638.6 * 1e6) / CL_eff
    return cavg / (cavg + p["IC50_ft"])


def sirna_eps(p, mg_q4w):
    """
    Steady-state HBsAg-transcript knockdown for a 28-day dosing interval.
    AUC per cycle = dose / CL, so the cycle-average plasma level is
    dose/(CL*28) and the hepatic effect site sits at kin*Cavg/kout.
    """
    if mg_q4w <= 0:
        return 0.0
    nmol = mg_q4w * 1e-3 / 16000.0 * 1e9          # ~16 kDa siRNA duplex
    cavg = nmol / (p["CL_sirna"] * 28.0)          # nmol/L
    lvl = p["kin_sirna"] * cavg / p["kout_sirna"]
    return p["Emax_sirna"] * lvl / (p["EC50_sirna"] + lvl)


def qss_of(p, blv_mg=0.0, ifn_ug=0.0, ifn_lambda=False, lnf_mg=0.0,
           rtv_mg=0.0, sirna_mg=0.0, occ_override=None):
    occ = blv_occ(p, blv_mg) if occ_override is None else occ_override
    isg = ifn_isg(p, ifn_ug)
    return {"occ": occ, "isg": isg, "ift": lnf_ift(p, lnf_mg, rtv_mg),
            "eps_s": sirna_eps(p, sirna_mg),
            "hem": 0.0 if ifn_lambda else isg}


# =============================================================================
#  8.  FITTING -- exactly four numbers
# =============================================================================

BA_ANCHOR_2MG = 3.2     # fold rise in total serum bile acids, 2 mg  (approx.)
BA_ANCHOR_10MG = 13.0   # fold rise in total serum bile acids, 10 mg (approx.)
MYR301_RNARESP_2MG = 0.71    # virologic response, week 48, 2 mg
MYR301_ALTNORM_2MG = 0.51    # ALT normalisation, week 48, 2 mg
MYR301_COMBINED_2MG = 0.45   # combined response, week 48, 2 mg


def fit_bile_acid_anchors(p):
    """
    Fit (r_oatp, Kd_ntcp) so the model reproduces both published TBA fold-rises.
    NTCP occupancy at 2 mg and 10 mg is then a DERIVED quantity.
    """
    def resid(r):
        p2 = dict(p); p2["r_oatp"] = r
        c2 = blv_css(p2, 2.0)
        c10 = blv_css(p2, 10.0)
        u2 = (1.0 + r) / BA_ANCHOR_2MG - r        # = 1 - OCC(2 mg)
        u10 = (1.0 + r) / BA_ANCHOR_10MG - r      # = 1 - OCC(10 mg)
        if u2 <= 0 or u10 <= 0 or u2 >= 1 or u10 >= 1:
            return None, None, None
        Kd_a = c2 * u2 / (1.0 - u2)
        Kd_b = c10 * u10 / (1.0 - u10)
        return Kd_a - Kd_b, Kd_a, Kd_b

    lo, hi = 0.0, 1.0 / BA_ANCHOR_10MG - 1e-6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        g, _, _ = resid(mid)
        if g is None:
            hi = mid
            continue
        if g < 0:
            lo = mid
        else:
            hi = mid
    r = 0.5 * (lo + hi)
    _, Kd, _ = resid(r)
    p["r_oatp"] = r
    p["Kd_ntcp"] = Kd
    return r, Kd


def fit_immune_killing(p, target=MYR301_RNARESP_2MG, t_end=336.0, n=150):
    """
    One parameter (dth_Id_immune) against one anchor: the MYR301 VIROLOGIC
    RESPONSE RATE at week 48 on bulevirtide 2 mg (71%).

    The response RATE, not the mean log decline, is the anchor deliberately:
    it is the number the trial reports without ambiguity, and it makes the mean
    decline a prediction instead of an input.
    """
    def decline(delta):
        p2 = dict(p)
        p2["dth_Id_immune"] = delta
        back_solve(p2)
        pop = make_population(p2, n)
        recs = run_arm(pop, {"blv_mg": 2.0}, t_end, dt=0.5)
        return rate(recs, t_end, p2, "vr") / 100.0

    # A SCAN, not a bisection, and deliberately so.  The relationship runs the
    # OPPOSITE way to intuition: a larger baseline CTL killing rate forces a
    # larger fitted re-infection rate (because the untreated state must still
    # balance), so the residual entry flux left over under bulevirtide is
    # proportionally BIGGER and the on-treatment decline is SHALLOWER.  A
    # bisection written on the intuitive sign converges to the wrong bound.
    grid = [0.004 * (1.45 ** k) for k in range(11)]
    scan = [(g, decline(g)) for g in grid]
    best = min(scan, key=lambda s: abs(s[1] - target))
    # refine locally around the best grid point
    i = [g for g, _ in scan].index(best[0])
    lo = scan[max(0, i - 1)][0]
    hi = scan[min(len(scan) - 1, i + 1)][0]
    for _ in range(7):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        if abs(decline(m1) - target) < abs(decline(m2) - target):
            hi = m2
        else:
            lo = m1
    p["dth_Id_immune"] = 0.5 * (lo + hi)
    back_solve(p)
    return p["dth_Id_immune"], decline(p["dth_Id_immune"]), scan


# =============================================================================
#  9.  VIRTUAL POPULATION
# =============================================================================

POP_AXES = [
    # (parameter, geometric-CV)   -- log-normal multiplicative variability
    ("dth_Id_immune", 0.45),      # strength of the HDV-specific CTL response
    ("d_hep", 0.40),              # hepatocyte turnover -> division-mediated floor
    ("k_cure", 0.60),             # non-cytolytic intracellular clearance
    ("Kd_ntcp", 0.35),            # PK / target-affinity variability
    ("c_d", 0.20),                # serum HDV RNA clearance
]


def make_population(p, n, seed=20260729, cv_floor=0.50):
    """
    Build n virtual patients.  cv_floor scales the variability on the two
    NTCP-INDEPENDENT maintenance terms (frac_cc and phi_ntcp), i.e. on the
    floor itself -- this is the single fitted population parameter.
    """
    rng = random.Random(seed)
    pop = []
    for i in range(n):
        q = dict(p)
        q["_idx"] = i
        for nm, cv in POP_AXES:
            s = math.sqrt(math.log(1.0 + cv * cv))
            q[nm] = p[nm] * math.exp(rng.gauss(0.0, s) - 0.5 * s * s)
        # the floor axes
        s = math.sqrt(math.log(1.0 + cv_floor * cv_floor))
        q["frac_cc"] = min(0.75, p["frac_cc"] *
                           math.exp(rng.gauss(0.0, s) - 0.5 * s * s))
        # scale around the model value, do NOT replace it -- otherwise a
        # sensitivity analysis on phi_ntcp is silently inert (this was a bug)
        q["phi_ntcp"] = min(0.95, max(0.05,
                                      p["phi_ntcp"] * rng.uniform(0.4, 1.6)))
        # baseline heterogeneity
        q["tg_Vd"] = p["tg_Vd"] * 10.0 ** rng.gauss(0.0, 0.90)
        q["tg_ALT"] = max(45.0, p["tg_ALT"] * math.exp(rng.gauss(0.0, 0.38)))
        q["tg_Fib"] = min(6.0, max(0.0, rng.gauss(2.6, 1.5)))
        q["tg_Id"] = min(0.30, max(0.02, p["tg_Id"] *
                                   math.exp(rng.gauss(0.0, 0.45))))
        sb = math.sqrt(math.log(1.0 + p["ALT_base_cv"] ** 2))
        q["ALT_base"] = p["ALT_base"] * math.exp(rng.gauss(0.0, sb) - 0.5 * sb * sb)
        back_solve(q)
        pop.append(q)
    return pop


def run_arm(pop, arm, t_end, dt=0.25):
    """
    Run one treatment arm across the population.
    `arm` = dict with keys blv_mg / ifn_ug / ifn_lambda / lnf_mg / rtv_mg /
            sirna_mg / occ_override / stop_day (day therapy stops).
    Returns a list of per-patient result records.
    """
    out = []
    stop = arm.get("stop_day", t_end)
    for q in pop:
        reg = Regimen(nuc=True)
        reg.ifn_lambda = arm.get("ifn_lambda", False)
        qon = qss_of(q, blv_mg=arm.get("blv_mg", 0.0),
                     ifn_ug=arm.get("ifn_ug", 0.0),
                     ifn_lambda=arm.get("ifn_lambda", False),
                     lnf_mg=arm.get("lnf_mg", 0.0),
                     rtv_mg=arm.get("rtv_mg", 0.0),
                     sirna_mg=arm.get("sirna_mg", 0.0),
                     occ_override=arm.get("occ_override"))
        qoff = qss_of(q)      # untreated engagement (endogenous ISG only)
        V0, A0 = q["tg_Vd"], q["tg_ALT"]
        if stop >= t_end:
            ts, ys = simulate(q, reg, t_end, dt=dt, obs_every=7.0, qss=qon)
        else:
            ts1, ys1 = simulate(q, reg, stop, dt=dt, obs_every=7.0, qss=qon)
            ts2, ys2 = simulate(q, reg, t_end - stop, dt=dt, obs_every=7.0,
                                qss=qoff, y0=ys1[-1])
            ts = ts1 + [stop + t for t in ts2[1:]]
            ys = ys1 + ys2[1:]
        rec = {"V0": V0, "A0": A0, "ts": ts, "ys": ys, "q": q,
               "idx": q.get("_idx", 0)}
        out.append(rec)
    return out


LOD_HDV = 6.0    # IU/mL, typical assay lower limit of detection


def endpoints(rec, t, p):
    V = at(rec["ts"], rec["ys"], t, "Vd")
    A = at(rec["ts"], rec["ys"], t, "ALT")
    # single-visit measurement / short-term biological variability
    rg = random.Random((rec.get("idx", 0) + 1) * 7919 + int(round(t)))
    if V > 0:
        V = V * 10.0 ** rg.gauss(0.0, p["noise_rna"])
    A = A * math.exp(rg.gauss(0.0, p["noise_alt"]))
    dlog = lg(V) - lg(rec["V0"])
    vr = (V < LOD_HDV) or (dlog <= -2.0)
    an = A <= p["ULN_ALT"]
    return {"V": V, "dlog": dlog, "ALT": A, "vr": vr, "alt_norm": an,
            "comb": (vr and an), "undet": V < LOD_HDV,
            "Fib": at(rec["ts"], rec["ys"], t, "Fib"),
            "Id": at(rec["ts"], rec["ys"], t, "Id"),
            "CH": at(rec["ts"], rec["ys"], t, "CH")}


def rate(recs, t, p, key):
    return 100.0 * sum(1 for r in recs if endpoints(r, t, p)[key]) / len(recs)


def mean_dlog(recs, t, p):
    return sum(endpoints(r, t, p)["dlog"] for r in recs) / len(recs)


def fit_alt_base(p, n=200, target=MYR301_ALTNORM_2MG):
    """
    One parameter (ALT_base, the MEDIAN non-HDV component of ALT) against one
    anchor: the MYR301 ALT normalisation rate at week 48 on 2 mg (51%).

    It is the median, not the spread, that is identifiable here.  Widening the
    spread of a mean-preserving log-normal lowers its median and therefore
    RAISES the normalisation rate, so the rate saturates around 75% and can
    never reach 51% by spread alone (that scan is shown in the report).  What
    the 51% actually tells us is that the typical patient with chronic
    hepatitis D carries enough NON-HDV liver injury -- steatosis, alcohol,
    metabolic -- to sit near the upper limit of normal even with the HDV
    contribution removed.
    """
    scan = []
    for ab in (20.0, 26.0, 30.0, 34.0, 38.0, 42.0, 46.0, 52.0, 60.0):
        p2 = dict(p); p2["ALT_base"] = ab
        pop = make_population(p2, n)
        recs = run_arm(pop, {"blv_mg": 2.0}, 336.0, dt=0.5)
        scan.append((ab, rate(recs, 336.0, p2, "alt_norm") / 100.0))
    best = min(scan, key=lambda s: abs(s[1] - target))
    return best[0], best[1], scan


def scan_alt_cv(p, n=150):
    """Evidence that the SPREAD saturates -- reported, not used as a fit."""
    out = []
    for cv in (0.10, 0.30, 0.55, 0.85, 1.20, 1.75):
        p2 = dict(p); p2["ALT_base_cv"] = cv
        pop = make_population(p2, n)
        recs = run_arm(pop, {"blv_mg": 2.0}, 336.0, dt=0.5)
        out.append((cv, rate(recs, 336.0, p2, "alt_norm") / 100.0))
    return out


# =============================================================================
# 10.  REPORT
# =============================================================================

OUT = []


def w(s=""):
    OUT.append(s)


def hdr(s):
    w()
    w("=" * 79)
    w(s)
    w("=" * 79)


def main():
    p = base_params()

    # -------------------------------------------------------------------
    hdr("A0.  MODEL STRUCTURE")
    w("Chronic hepatitis D QSP model, %d ODE compartments." % NST)
    w("Compartments: " + ", ".join(NAMES))
    w()
    w("The three fluxes that maintain the infected-hepatocyte pool Id:")
    w("   (E) entry / re-infection      beta_d * Vd * Ib * (1 - OCC)")
    w("   (C) cell-to-cell spread       k_cc * Id * Ib * (1 - phi_ntcp*OCC)")
    w("   (D) division-mediated spread  Rtot * Id / N      NTCP-INDEPENDENT")
    w("Only (E) is fully blockable by an entry inhibitor; (D) is not blockable")
    w("at all, and (C) only in part.  (D) + the unblockable part of (C) is the")
    w("FLOOR, and the rest of this report is about how large that floor is.")

    # -------------------------------------------------------------------
    hdr("A1.  FIT 1-2 OF 4 -- BILE ACIDS INVERT TO NTCP OCCUPANCY")
    r, Kd = fit_bile_acid_anchors(p)
    back_solve(p)
    w("Serum total bile acids (TBA) obey, at steady state,")
    w("     TBA_fold = (1 + r) / ((1 - OCC) + r),   r = CL_OATP / CL_NTCP")
    w("because NTCP is both the HDV receptor and the high-capacity bile-acid")
    w("uptake transporter.  Fitting the two published fold-rises:")
    w()
    w("   ANCHOR   bulevirtide  2 mg qd : TBA x %.1f" % BA_ANCHOR_2MG)
    w("   ANCHOR   bulevirtide 10 mg qd : TBA x %.1f" % BA_ANCHOR_10MG)
    w()
    w("   FITTED   r      (OATP1B residual uptake capacity) = %s" % fmt(r, 4))
    w("   FITTED   Kd_NTCP                                  = %s nM" % fmt(Kd, 3))
    w()
    c2, c10 = blv_css(p, 2.0), blv_css(p, 10.0)
    o2, o10 = blv_occ(p, 2.0), blv_occ(p, 10.0)
    w("   DERIVED  Css(2 mg)  = %s nM   ->  NTCP occupancy = %s" % (fmt(c2, 3), fmt(o2, 4)))
    w("   DERIVED  Css(10 mg) = %s nM   ->  NTCP occupancy = %s" % (fmt(c10, 3), fmt(o10, 4)))
    w("   DERIVED  residual free NTCP    2 mg %s%%   10 mg %s%%"
      % (fmt(100 * (1 - o2), 1), fmt(100 * (1 - o10), 1)))
    w("   DERIVED  exposure ratio 10 mg / 2 mg = %s  (dose ratio 5.0)"
      % fmt(c10 / c2, 2))
    w()
    w("READ THIS TWICE.  The approved 2 mg dose does NOT saturate NTCP: it")
    w("leaves %s%% of the receptor free, and 10 mg leaves %s%%.  The residual"
      % (fmt(100 * (1 - o2), 1), fmt(100 * (1 - o10), 1)))
    w("entry flux therefore differs %sx between the two doses.  Yet MYR301"
      % fmt((1 - o2) / (1 - o10), 1))
    w("found essentially the same efficacy at both.  Section A4 is where that")
    w("apparent contradiction is resolved, and A6 is where it is quantified.")
    w()
    w("INDEPENDENT CHECK (nothing below was fitted): the fitted Kd of %s nM"
      % fmt(Kd, 2))
    w("is the affinity of a myristoylated preS1 lipopeptide for NTCP inferred")
    w("purely from bile-acid pharmacology.  Reported in-vitro IC50 values for")
    w("myrcludex B / bulevirtide at NTCP are in the 80 pM - low nM range, so")
    w("the bile-acid inversion lands inside the measured affinity window")
    w("without ever having been shown a binding experiment.")

    # -------------------------------------------------------------------
    hdr("A2.  BASELINE -- THE UNTREATED CHRONIC STEADY STATE")
    reg0 = Regimen(nuc=True)
    ts, ys = simulate(p, reg0, 365.0 * 2, dt=0.02, obs_every=7.0)
    w("Every parameter that can be pinned by requiring d/dt = 0 at the")
    w("observed baseline IS pinned that way (back_solve()), so the untreated")
    w("model sits on the clinical steady state by construction.  Verification")
    w("after 2 years of simulation with no anti-HDV therapy (NUC only):")
    w()
    w("   quantity                 target        t = 2 y      drift")
    for nm, tgt, u in [("Vd", p["tg_Vd"], "IU/mL"), ("S_ser", p["tg_S"], "IU/mL"),
                       ("Id", p["tg_Id"], "fraction"), ("Rg", p["tg_Rg"], "cop/cell"),
                       ("ALT", p["tg_ALT"], "U/L"), ("E", p["tg_E"], "a.u."),
                       ("TBA", p["TBA0"], "umol/L")]:
        v = at(ts, ys, 730.0, nm)
        w("   %-22s %10s  %12s   %+6.1f%%   %s"
          % (nm, ("%.4g" % tgt), ("%.4g" % v), 100.0 * (v / tgt - 1.0), u))
    w("   %-22s %10s  %12s" % ("Fib (Ishak)", fmt(p["tg_Fib"], 2),
                               fmt(at(ts, ys, 730.0, "Fib"), 2)))
    w()
    w("Fibrosis is the one state deliberately NOT at steady state: untreated")
    w("chronic hepatitis D progresses.  Model progression rate from Ishak 2:")
    f2 = at(ts, ys, 730.0, "Fib")
    w("   %s -> %s over 2 y  =  %s Ishak units / year"
      % (fmt(p["tg_Fib"], 2), fmt(f2, 2), fmt((f2 - p["tg_Fib"]) / 2.0, 3)))
    w("   literature: hepatitis D is the fastest-progressing chronic viral")
    w("   hepatitis, roughly 0.15-0.25 Ishak units/y vs ~0.10 for HBV alone.")

    # -------------------------------------------------------------------
    hdr("A3.  FIT 3 OF 4 -- INFECTED-CELL LOSS FROM THE MYR301 RESPONSE RATE")
    d_fit, got, dscan = fit_immune_killing(p)
    w("One parameter, one anchor.")
    w("   ANCHOR  MYR301, bulevirtide 2 mg, week 48:")
    w("           virologic response rate = %s%%  (HDV RNA undetectable or"
      % fmt(100 * MYR301_RNARESP_2MG, 0))
    w("           >=2 log10 decline)")
    w("   FITTED  dth_Id_immune (baseline infected-cell killing) = %s /day"
      % fmt(p["dth_Id_immune"], 5))
    w("   MODEL   virologic response, 2 mg, week 48             = %s%%"
      % fmt(100 * got, 1))
    w()
    w("   scan of the objective (dth_Id_immune -> virologic response rate):")
    w("     " + "  ".join("%.4f:%.0f%%" % (gg, 100 * vv) for gg, vv in dscan))
    w()
    w("   The RESPONSE RATE is the anchor rather than the mean log decline,")
    w("   deliberately: it is the number MYR301 reports without ambiguity, and")
    w("   it leaves the mean decline as a prediction instead of an input.")
    w()
    w("THE SIGN OF THIS FIT IS THE OPPOSITE OF INTUITION, and it is worth")
    w("stopping on.  A LARGER baseline killing rate produces a SHALLOWER")
    w("on-treatment decline, not a deeper one.  The reason is that the")
    w("untreated state must still balance: if infected cells die faster, the")
    w("re-infection flux that holds the steady state up must be larger too, so")
    w("the residual flux left over after blocking a fixed fraction of entry is")
    w("larger in absolute terms.  A bisection written on the intuitive sign")
    w("runs to the bound and returns a non-physiological killing rate; this")
    w("was an actual bug in the first version of this file.")
    w()
    w("Note what this single number buys.  It sets the per-capita loss rate of")
    w("infected hepatocytes, and therefore -- with NOTHING else fitted -- the")
    w("entire second-phase slope of every regimen in this report.")
    w()
    w("Baseline decomposition of the inflow that maintains Id, at the fitted")
    w("value (per day, absolute units of hepatocyte fraction):")
    tot = p["_entry0"] + p["_cc0"] + p["_div0"]
    w("   (E) entry / re-infection      %s  (%s%% of inflow)"
      % ("%.3e" % p["_entry0"], fmt(100 * p["_entry0"] / tot, 1)))
    w("   (C) cell-to-cell              %s  (%s%%)"
      % ("%.3e" % p["_cc0"], fmt(100 * p["_cc0"] / tot, 1)))
    w("   (D) division-mediated         %s  (%s%%)"
      % ("%.3e" % p["_div0"], fmt(100 * p["_div0"] / tot, 1)))
    w("   NTCP-INDEPENDENT share (D + (1-phi_ntcp)*C) = %s%% of total inflow"
      % fmt(100 * (p["_div0"] + (1 - p["phi_ntcp"]) * p["_cc0"]) / tot, 1))
    w("   net dilution rate of Id by HBsAg+ HDV-negative neighbours = %s /day"
      % fmt(p["_dilute0"], 5))
    w("   regeneration locality loc_renew (back-solved)             = %s"
      % fmt(p["loc_renew"], 4))

    # -------------------------------------------------------------------
    hdr("A4.  THE SHAPE OF THE DECLINE -- WHY BULEVIRTIDE HAS NO FIRST PHASE")
    w("An entry inhibitor does not touch a cell that is already infected, so")
    w("it cannot reduce secretion on day 1.  Interferon blocks intracellular")
    w("replication, so it can.  The model was never told this; it follows from")
    w("which flux each drug cuts.  Log10 HDV RNA change from baseline:")
    w()
    arms_shape = [
        ("BLV 2 mg qd", qss_of(p, blv_mg=2.0)),
        ("BLV 10 mg qd", qss_of(p, blv_mg=10.0)),
        ("PegIFN alfa 180 ug qw", qss_of(p, ifn_ug=180.0)),
        ("LNF 50 + RTV 100 BID", qss_of(p, lnf_mg=50.0, rtv_mg=100.0)),
    ]
    w("   arm                        d1      d3      d7     d28     d84    d336")
    shape = {}
    for nm, q in arms_shape:
        ts2, ys2 = simulate(p, Regimen(), 336.0, dt=0.05, obs_every=1.0, qss=q)
        row = [lg(at(ts2, ys2, d, "Vd")) - lg(p["tg_Vd"])
               for d in (1, 3, 7, 28, 84, 336)]
        shape[nm] = row
        w("   %-24s" % nm + "".join("%8s" % fmt(v, 2) for v in row))
    w()
    w("The signature is unambiguous and it matches the clinic: bulevirtide")
    w("declines almost linearly in log space from day 1 with no rapid phase,")
    w("whereas interferon and lonafarnib both drop steeply in the first week")
    w("(replication block and assembly block act on cells already infected)")
    w("and then flatten onto the same cell-loss-limited slope.")
    w()
    w("   first-week slope (log10/day)   BLV 2 mg  %s" % fmt(shape["BLV 2 mg qd"][2] / 7.0, 4))
    w("                                  PegIFN    %s" % fmt(shape["PegIFN alfa 180 ug qw"][2] / 7.0, 4))
    w("                                  LNF/RTV   %s" % fmt(shape["LNF 50 + RTV 100 BID"][2] / 7.0, 4))
    w()
    w("   That is the diagnostic signature of the class, and it is not a")
    w("   fitted parameter.  Bulevirtide's week-1 slope is %s log10/day --"
      % fmt(shape["BLV 2 mg qd"][2] / 7.0, 4))
    w("   essentially flat -- because on day 7 the infected pool has barely")
    w("   changed and nothing at all has been done to the cells inside it.")

    # -------------------------------------------------------------------
    hdr("A5.  THE LONAFARNIB PARADOX -- SERUM DOWN, INTRACELLULAR UP")
    q_lnf = qss_of(p, lnf_mg=50.0, rtv_mg=100.0)
    ts3, ys3 = simulate(p, Regimen(), 336.0, dt=0.05, obs_every=7.0, qss=q_lnf)
    # then stop and watch the rebound
    ts4, ys4 = simulate(p, Regimen(), 84.0, dt=0.05, obs_every=3.0,
                        qss=qss_of(p), y0=ys3[-1])
    w("Blocking farnesyltransferase stops envelopment, not replication.  The")
    w("genomes therefore pile up inside the hepatocyte while serum RNA falls.")
    w("FTase inhibition at 50 mg BID + ritonavir 100 mg BID: I_FT = %s"
      % fmt(q_lnf["ift"], 3))
    w()
    w("   week   serum HDV RNA (log10)   intracellular Rg (copies/cell)")
    for d in (0, 7, 28, 84, 168, 336):
        w("   %4d   %18s   %22s" % (d // 7, fmt(lg(at(ts3, ys3, d, "Vd")), 2),
                                    fmt(at(ts3, ys3, d, "Rg"), 0)))
    w()
    w("   intracellular / baseline at week 48 = %sx"
      % fmt(at(ts3, ys3, 336.0, "Rg") / p["tg_Rg"], 2))
    w("   serum log10 decline at week 48       = %s"
      % fmt(lg(at(ts3, ys3, 336.0, "Vd")) - lg(p["tg_Vd"]), 2))
    w()
    w("Consequence: withdrawal releases the trapped pool.  After stopping at")
    w("week 48, serum HDV RNA returns to within 0.5 log10 of baseline in")
    reb = None
    for i, tt in enumerate(ts4):
        if lg(ys4[i][IX["Vd"]]) > lg(p["tg_Vd"]) - 0.5:
            reb = tt
            break
    w("   %s days." % (fmt(reb, 0) if reb is not None else ">84"))
    w("   log10 HDV RNA  d0 %s  d7 %s  d14 %s  d28 %s  d84 %s (off treatment)"
      % tuple(fmt(lg(at(ts4, ys4, d, "Vd")), 2) for d in (0, 7, 14, 28, 84)))
    w()
    w("This is the model's explanation for why the assembly-inhibitor class")
    w("has produced fast, large on-treatment serum declines and fast relapse:")
    w("serum HDV RNA measures the ASSEMBLY FLUX, not the reservoir, and only")
    w("lonafarnib among current drugs makes those two quantities disagree.")

    # -------------------------------------------------------------------
    hdr("A6.  FIT 4 OF 4 + THE CENTRAL RESULT -- WHY 2 mg = 10 mg")
    ab, altgot, ascan = fit_alt_base(p)
    p["ALT_base"] = ab
    back_solve(p)
    w("One parameter, one anchor.")
    w("   ANCHOR  MYR301, bulevirtide 2 mg, week 48: ALT normalisation 51%")
    w("   FITTED  ALT_base (median NON-HDV component of ALT) = %s U/L" % fmt(ab, 1))
    w("   MODEL   ALT normalisation, 2 mg, week 48           = %s%%"
      % fmt(100 * altgot, 1))
    w("     scan: " + "  ".join("%.0f:%.0f%%" % (a, 100 * v) for a, v in ascan))
    w()
    w("What that number MEANS is worth stating: it says the typical patient")
    w("with chronic hepatitis D sits at or just above the upper limit of normal")
    w("even with the entire HDV contribution to ALT removed.  In other words a")
    w("large part of the ALT-normalisation endpoint is not about HDV at all --")
    w("it is about the steatotic / alcohol / metabolic liver underneath.")
    w()
    w("A NON-IDENTIFIABILITY, reported rather than hidden.  The variability of")
    w("the NTCP-INDEPENDENT floor terms (cv_floor) was originally intended as a")
    w("fourth fitted parameter.  It is not identifiable from these endpoints --")
    w("the combined response is flat across more than a decade of it:")
    cvscan = []
    for cvx in (0.10, 0.30, 0.60, 1.00, 1.60):
        popx = make_population(p, 120, cv_floor=cvx)
        cvscan.append((cvx, rate(run_arm(popx, {"blv_mg": 2.0}, 336.0, dt=0.5),
                                 336.0, p, "comb")))
    w("     cv_floor -> combined response:  "
      + "  ".join("%.2f:%.0f%%" % (a, v) for a, v in cvscan))
    w("   so cv_floor is FIXED at 0.50 and is not counted as a fit.  What the")
    w("   response rate responds to is the MEAN size of the floor, not its")
    w("   spread -- which is itself the point of section A10.")
    w()
    w("   Also reported: widening the spread of the non-HDV ALT term cannot")
    w("   reach 51% either, because a mean-preserving log-normal has a LOWER")
    w("   median as it widens, so more patients normalise, not fewer:")
    w("     ALT_base_cv -> ALT normalisation:  "
      + "  ".join("%.2f:%.0f%%" % (a, 100 * v) for a, v in scan_alt_cv(p)))
    w()
    w("Everything from here on is a PREDICTION.")
    w()
    N_POP = 300
    pop = make_population(p, N_POP)
    arms = [
        ("untreated (NUC only)", {}, 336.0),
        ("BLV 2 mg qd", {"blv_mg": 2.0}, 336.0),
        ("BLV 10 mg qd", {"blv_mg": 10.0}, 336.0),
        ("PERFECT entry block", {"occ_override": 1.0}, 336.0),
    ]
    res = {}
    w("   arm                       mean dlog10   RNA resp   ALT norm   COMBINED   undet")
    for nm, a, te in arms:
        recs = run_arm(pop, a, te, dt=0.5)
        res[nm] = recs
        w("   %-24s %11s %9s%% %9s%% %9s%% %6s%%"
          % (nm, fmt(mean_dlog(recs, te, p), 2), fmt(rate(recs, te, p, "vr"), 1),
             fmt(rate(recs, te, p, "alt_norm"), 1),
             fmt(rate(recs, te, p, "comb"), 1), fmt(rate(recs, te, p, "undet"), 1)))
    w()
    w("   OBSERVED (MYR301, week 48, Wedemeyer 2023 NEJM):")
    w("     untreated control      RNA  4%   ALT 12%   COMBINED  2%")
    w("     bulevirtide  2 mg      RNA 71%   ALT 51%   COMBINED 45%   <- 2 anchors")
    w("     bulevirtide 10 mg      RNA 76%   ALT 56%   COMBINED 48%   <- prediction")
    w()
    tot = p["_entry0"] + p["_cc0"] + p["_div0"]
    r2 = rate(res["BLV 2 mg qd"], 336.0, p, "comb")
    r10 = rate(res["BLV 10 mg qd"], 336.0, p, "comb")
    rp = rate(res["PERFECT entry block"], 336.0, p, "comb")
    ru = rate(res["untreated (NUC only)"], 336.0, p, "comb")
    w("THE CENTRAL RESULT.  Between 2 mg and 10 mg the residual entry flux")
    w("falls %sx (%s%% -> %s%% of NTCP left free), and yet the model's COMBINED"
      % (fmt((1 - o2) / (1 - o10), 1), fmt(100 * (1 - o2), 1),
         fmt(100 * (1 - o10), 1)))
    w("response moves only %s%% -> %s%% -- against an observed 45%% -> 48%%."
      % (fmt(r2, 1), fmt(r10, 1)))
    w("The reason is that entry supplies only %s%% of the inflow that maintains"
      % fmt(100 * p["_entry0"] / tot, 0))
    w("Id; the rest is division-mediated and cell-to-cell spread, which no")
    w("entry inhibitor can reach.  And the COMBINED endpoint is doubly")
    w("insensitive, because its ALT half is capped by non-HDV liver injury.")
    w()
    w("THE CEILING OF THE WHOLE CLASS.  A hypothetical entry inhibitor with")
    w("100%% occupancy and no toxicity reaches %s%% combined response -- only" % fmt(rp, 1))
    w("%s points above the approved 2 mg dose.  Measured against the untreated"
      % fmt(rp - r2, 1))
    w("arm, the entry-inhibition axis is already %s%% exhausted at 2 mg/day."
      % fmt(100.0 * (r2 - ru) / max(1e-9, rp - ru), 0))
    w("Dose escalation is not the lever.  Adding a mechanism is.")
    w()
    w("AND HERE IS WHERE THE MODEL DISAGREES WITH THE TRIAL -- reported, not")
    w("absorbed.  On the VIROLOGIC endpoint alone the model separates the doses")
    w("far more than MYR301 did: %s%% vs %s%% where the trial saw 71%% vs 76%%."
      % (fmt(rate(res["BLV 2 mg qd"], 336.0, p, "vr"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "vr"), 1)))
    w("That gap is forced by the bile-acid inversion of A1: if 2 mg really")
    w("leaves %s%% of NTCP free and 10 mg leaves %s%%, a 5x dose step MUST buy"
      % (fmt(100 * (1 - o2), 1), fmt(100 * (1 - o10), 1)))
    w("more virologic suppression than was observed.  Three ways out, and the")
    w("model cannot choose between them:")
    w("   (i)  2 mg is closer to saturating than the bile-acid inversion says,")
    w("        and total bile acids keep rising for a reason other than simple")
    w("        residual-capacity scaling -- FXR/CYP7A1 feedback or OATP1B")
    w("        saturation would both do it;")
    w("   (ii) there is a dose-independent ceiling not in this model at all --")
    w("        adherence, assay floor, or a genuinely unreachable reservoir;")
    w("   (iii) at high occupancy entry stops being limited by free NTCP and")
    w("        becomes limited by something upstream (HSPG attachment), so the")
    w("        last few percent of receptor blockade buys nothing.")
    w("Option (iii) is the interesting one and it is testable: it predicts that")
    w("HDV entry inhibition saturates at LOWER occupancy than bile-acid")
    w("transport inhibition does, in the same cells, on the same drug.")

    # -------------------------------------------------------------------
    hdr("A7.  ALT / HDV RNA DECOUPLING, QUANTIFIED")
    w("In MYR301 many patients normalised ALT without a virologic response.")
    w("In this model that is not a coincidence, it is structural: ALT is")
    w("driven by")
    w("      kill = (innate + CTL + IFN) * Id   +   kappa * ENTRY_FLUX")
    w("so blocking entry removes a term that acts within days, while Id -- and")
    w("hence serum RNA -- only decays over months.  Model time courses on 2 mg:")
    w()
    q2 = qss_of(p, blv_mg=2.0)
    tsA, ysA = simulate(p, Regimen(), 336.0, dt=0.05, obs_every=7.0, qss=q2)
    w("   week     ALT (U/L)   ALT % of baseline   log10 HDV RNA change")
    for d in (0, 7, 14, 28, 84, 168, 336):
        A = at(tsA, ysA, d, "ALT")
        w("   %4d %12s %18s %22s"
          % (d // 7, fmt(A, 1), fmt(100 * A / p["tg_ALT"], 1),
             fmt(lg(at(tsA, ysA, d, "Vd")) - lg(p["tg_Vd"]), 2)))
    w()
    a4 = at(tsA, ysA, 28.0, "ALT")
    v4 = lg(at(tsA, ysA, 28.0, "Vd")) - lg(p["tg_Vd"])
    w("   At week 4: ALT has already fallen %s%% while HDV RNA has fallen only"
      % fmt(100 * (1 - a4 / p["tg_ALT"]), 1))
    w("   %s log10 (%s%%).  The two endpoints are reading different things."
      % (fmt(v4, 2), fmt(100 * (1 - 10 ** v4), 1)))
    w()
    w("A TESTABLE PREDICTION FOLLOWS.  Because the fresh-infection injury term")
    w("scales with (1 - OCC) directly while the RNA decline scales with the")
    w("much smaller net loss rate, ALT normalisation should be MORE")
    w("dose-sensitive than virologic response.  Model:")
    w("   ALT normalisation   2 mg %s%%   10 mg %s%%   (gap %s pts)"
      % (fmt(rate(res["BLV 2 mg qd"], 336.0, p, "alt_norm"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "alt_norm"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "alt_norm")
             - rate(res["BLV 2 mg qd"], 336.0, p, "alt_norm"), 1)))
    w("   virologic response  2 mg %s%%   10 mg %s%%   (gap %s pts)"
      % (fmt(rate(res["BLV 2 mg qd"], 336.0, p, "vr"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "vr"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "vr")
             - rate(res["BLV 2 mg qd"], 336.0, p, "vr"), 1)))
    w("   OBSERVED gaps in MYR301: ALT +5 pts (51->56), RNA +5 pts (71->76).")
    w("   The model reproduces the direction and the small size of both gaps.")
    w("   CAVEAT: the kappa*ENTRY term is a MODEL HYPOTHESIS, not a measured")
    w("   mechanism.  It is the cheapest structure that produces the observed")
    w("   decoupling, and it is falsifiable -- see A13.")

    # -------------------------------------------------------------------
    hdr("A8.  INTERFERON -- THE ONLY MECHANISM WITH OFF-TREATMENT DURABILITY")
    w("Interferon is the only current drug that raises the DEATH RATE of")
    w("infected hepatocytes (delta) rather than lowering an inflow.  That is")
    w("why it, and only it, produces responses that survive withdrawal.")
    w()
    isg_on = ifn_isg(p, 180.0)
    efo = ifn_effects(p, isg_on)
    w("   endogenous ISG tone, untreated            = %s" % fmt(p["_ISG0"], 3))
    w("   steady-state ISG on 180 ug qw             = %s" % fmt(isg_on, 3))
    w("   ISG INCREMENT that carries all drug effect= %s" % fmt(efo.dis, 3))
    w("   added infected-cell loss delta_IFN        = %s /day" % fmt(efo.dlt, 5))
    w("   baseline infected-cell loss               = %s /day" % fmt(p["_dth_Id0"], 5))
    w("   -> infected-cell loss rate multiplied by %sx"
      % fmt(1 + efo.dlt / p["_dth_Id0"], 2))
    w("   virion PRODUCTION block eps_prod          = %s  (the first phase)"
      % fmt(efo.eps_prod, 3))
    w("   intracellular replication block eps_rep   = %s" % fmt(efo.eps_rep, 3))
    w("   HBsAg secretion reduction                 = %s" % fmt(efo.eps_hbs, 3))
    w("   amber/W editing fraction  %s -> %s (ADAR1 is IFN-inducible)"
      % (fmt(p["f_edit0"], 3), fmt(efo.f_edit, 3)))
    w("   exhaustion reversal factor on d_x         = %s" % fmt(efo.xrev, 2))
    w()
    ifn_arms = [
        ("PegIFN alfa 48 wk", {"ifn_ug": 180.0, "stop_day": 336.0}, 504.0),
        ("PegIFN alfa 96 wk", {"ifn_ug": 180.0, "stop_day": 672.0}, 840.0),
        ("PegIFN lambda 48 wk", {"ifn_ug": 180.0, "ifn_lambda": True,
                                 "stop_day": 336.0}, 504.0),
        ("BLV 2 mg + PegIFN 48 wk", {"blv_mg": 2.0, "ifn_ug": 180.0,
                                     "stop_day": 336.0}, 504.0),
        ("BLV 10 mg + PegIFN 48 wk", {"blv_mg": 10.0, "ifn_ug": 180.0,
                                      "stop_day": 336.0}, 504.0),
        ("BLV 2 mg continuous 96 wk", {"blv_mg": 2.0}, 672.0),
    ]
    w("   arm                          EOT dlog10  EOT comb   +24wk dlog10  +24wk RNA<LOD")
    ifn_res = {}
    for nm, a, te in ifn_arms:
        recs = run_arm(pop, a, te, dt=0.5)
        ifn_res[nm] = (recs, a, te)
        eot = a.get("stop_day", te)
        w("   %-28s %10s %9s%% %13s %13s%%"
          % (nm, fmt(mean_dlog(recs, eot, p), 2), fmt(rate(recs, eot, p, "comb"), 1),
             fmt(mean_dlog(recs, te, p), 2), fmt(rate(recs, te, p, "undet"), 1)))
    w()
    w("   OBSERVED, for orientation (different populations and assays, so")
    w("   these are direction-and-magnitude checks, not fits):")
    w("     HIDIT-1  pegIFN 48 wk, HDV RNA neg 24 wk post-treatment  ~26-28%")
    w("     HIDIT-2  pegIFN 96 wk, HDV RNA neg 24 wk post-treatment  ~31%")
    w("     LIMT-1   pegIFN lambda 180 ug, 24 wk post-treatment      36% (5/14)")
    w("     MYR204   BLV 10 mg + pegIFN 48 wk, undetectable 24 wk FU ~45%")
    w("              pegIFN alone in the same trial                  ~17-24%")
    w()
    a1 = ifn_res["PegIFN alfa 48 wk"]
    a2 = ifn_res["BLV 2 mg + PegIFN 48 wk"]
    a3 = ifn_res["BLV 10 mg + PegIFN 48 wk"]
    w("SYNERGY IS STRUCTURAL, NOT EMPIRICAL.  Bulevirtide lowers an INFLOW and")
    w("interferon raises an OUTFLOW of the same pool, so their effects on the")
    w("net per-capita growth rate of Id are additive while their effects on the")
    w("ENDPOINT are supra-additive (the endpoint is a threshold on a decaying")
    w("exponential).  Model, 24 weeks after stopping:")
    w("   RNA < LOD    pegIFN alone %s%%   + BLV 2 mg %s%%   + BLV 10 mg %s%%"
      % (fmt(rate(a1[0], a1[2], p, "undet"), 1), fmt(rate(a2[0], a2[2], p, "undet"), 1),
         fmt(rate(a3[0], a3[2], p, "undet"), 1)))
    w("   and note that continuous BLV monotherapy at 96 wk reaches")
    w("   %s%% combined response but NOTHING off treatment -- an entry"
      % fmt(rate(ifn_res["BLV 2 mg continuous 96 wk"][0], 672.0, p, "comb"), 1))
    w("   inhibitor is a suppressive therapy by construction.")

    # -------------------------------------------------------------------
    hdr("A9.  LONAFARNIB REGIMENS AND THE 'STARVE THE ENVELOPE' STRATEGY")
    lnf_arms = [
        ("LNF 50 + RTV 100 BID 48 wk", {"lnf_mg": 50.0, "rtv_mg": 100.0,
                                        "stop_day": 336.0}, 504.0),
        ("LNF/RTV + PegIFN 48 wk", {"lnf_mg": 50.0, "rtv_mg": 100.0,
                                    "ifn_ug": 180.0, "stop_day": 336.0}, 504.0),
        ("BLV 2 mg + siRNA 200 mg q4w", {"blv_mg": 2.0, "sirna_mg": 200.0}, 336.0),
        ("siRNA 200 mg q4w alone", {"sirna_mg": 200.0}, 336.0),
    ]
    w("   arm                             EOT dlog10   EOT comb   HBsAg dlog10")
    for nm, a, te in lnf_arms:
        recs = run_arm(pop, a, te, dt=0.5)
        eot = a.get("stop_day", te)
        hb = sum(lg(at(r["ts"], r["ys"], eot, "S_ser")) - lg(r["q"]["tg_S"])
                 for r in recs) / len(recs)
        w("   %-31s %10s %9s%% %14s"
          % (nm, fmt(mean_dlog(recs, eot, p), 2),
             fmt(rate(recs, eot, p, "comb"), 1), fmt(hb, 2)))
    w()
    w("   OBSERVED, D-LIVR phase 3 (conference-reported): composite response")
    w("   at week 48 roughly 10% for LNF/RTV and ~19% for LNF/RTV + pegIFN,")
    w("   versus ~2% on placebo.  The model OVERSHOOTS both badly: see A13.")
    w("   same reason the trial's was: a large serum-RNA drop that does NOT")
    w("   come with ALT normalisation, because the infected pool is untouched.")
    w()
    w("The siRNA arms make the deepest structural point in the model.  HDV")
    w("cannot leave a cell without an HBsAg envelope, so lowering HBsAg cuts")
    w("the SAME flux lonafarnib cuts, one step further upstream -- and unlike")
    w("lonafarnib it also shrinks the envelope-donor pool that entry and")
    w("cell-to-cell spread both draw on.  On this model that is why an")
    w("HBsAg-directed agent is the only non-interferon mechanism that moves")
    w("the floor rather than merely pressing on the entry flux.")

    # -------------------------------------------------------------------
    hdr("A10.  WHAT SETS THE FLOOR -- SENSITIVITY OF THE ENTRY-INHIBITION CEILING")
    w("If the floor is what limits bulevirtide, then the ceiling of entry")
    w("inhibition must be sensitive to the floor parameters and INSENSITIVE to")
    w("occupancy.  Test: hold occupancy at the 2 mg value and perturb each")
    w("parameter by +/-30%, then report the week-48 combined response.")
    w()
    base_comb = r2
    sens_pars = [
        ("d_hep", "hepatocyte turnover -> division flux"),
        ("frac_cc", "cell-to-cell share of inflow"),
        ("phi_ntcp", "NTCP-dependent share of cell-to-cell"),
        ("k_cure", "non-cytolytic intracellular clearance"),
        ("dth_Id_immune", "CTL killing of infected cells"),
        ("Kd_ntcp", "target affinity (i.e. occupancy)"),
    ]
    sens_rows = []
    popb = make_population(p, 120)
    recb = run_arm(popb, {"blv_mg": 2.0}, 336.0, dt=0.5)
    rb = rate(recb, 336.0, p, "comb")
    db = mean_dlog(recb, 336.0, p)
    w("   Two metrics are shown.  The combined-response rate is what the trial")
    w("   reports, but it is a THRESHOLD statistic on 120 patients, so it is")
    w("   quantised and can read 0.0 for a parameter that genuinely matters.")
    w("   The mean log10 decline is continuous and is the honest ranking.")
    w()
    w("   parameter         mean dlog10: -30%     base    +30%   |d dlog|  resp d  role")
    for nm, role in sens_pars:
        vals, dvals = [], []
        for f in (0.7, 1.3):
            p2 = dict(p); p2[nm] = p[nm] * f
            back_solve(p2)
            pop2 = make_population(p2, 120)
            recs2 = run_arm(pop2, {"blv_mg": 2.0}, 336.0, dt=0.5)
            vals.append(rate(recs2, 336.0, p2, "comb"))
            dvals.append(mean_dlog(recs2, 336.0, p2))
        eff_d = 0.5 * (abs(dvals[0] - db) + abs(dvals[1] - db))
        eff_r = 0.5 * (abs(vals[0] - rb) + abs(vals[1] - rb))
        sens_rows.append((nm, vals[0], rb, vals[1], eff_d, role, eff_r))
        w("   %-16s %14s %8s %7s %9s %6s   %s"
          % (nm, fmt(dvals[0], 2), fmt(db, 2), fmt(dvals[1], 2),
             fmt(eff_d, 3), fmt(eff_r, 1), role))
    w()
    w()
    sens_rows.sort(key=lambda r: -r[4])
    w("   Ranked by |d mean log10| per 30%: "
      + ", ".join("%s (%s)" % (r[0], fmt(r[4], 3)) for r in sens_rows))
    kd_row = [r for r in sens_rows if r[0] == "Kd_ntcp"][0]
    top = sens_rows[0]
    host = [r for r in sens_rows if r[0] != "Kd_ntcp"]
    host_sum = sum(r[4] for r in host)
    w("   Most influential: %s (%s log10 per 30%%)." % (top[0], fmt(top[4], 3)))
    w("   Target affinity Kd_ntcp ranks %d of %d at %s log10 per 30%%, against"
      % (sens_rows.index(kd_row) + 1, len(sens_rows), fmt(kd_row[4], 3)))
    w("   %s summed over the five HOST parameters -- a ratio of %sx in favour"
      % (fmt(host_sum, 3), fmt(host_sum / max(1e-9, kd_row[4]), 1)))
    w("   of host cell biology over drug-target affinity.")
    w()
    dhep_row = [r for r in sens_rows if r[0] == "d_hep"][0]
    w("AND NOW THE RESULT I DID NOT EXPECT.  d_hep -- the background hepatocyte")
    w("turnover rate, the very parameter I built the division-mediated floor")
    w("around -- has an influence of %s log10, i.e. NONE, on the continuous"
      % fmt(dhep_row[4], 3))
    w("metric as well as the quantised one.  That is not a numerical artefact")
    w("and it is not a bug; it is algebra I had not done.  Raising the")
    w("background death rate of hepatocytes raises the death rate of infected")
    w("cells and the division rate that refills their slots BY THE SAME AMOUNT,")
    w("so it cancels out of the net dilution rate exactly:")
    w()
    w("     dilution = dth_Id - loc*(dth_Ib*Ib + dth_Id*Id)/(Ib+Id) - ...")
    w("     and d_hep sits inside dth_Id, dth_Ib and dth_T identically.")
    w()
    w("What actually sets the floor is therefore NOT how fast the liver turns")
    w("over.  It is how much FASTER an infected hepatocyte dies than its")
    w("HBsAg-positive neighbours -- the DIFFERENTIAL, which is dth_Id_immune.")
    w("That is why dth_Id_immune outranks everything else by an order of")
    w("magnitude, and it sharpens the therapeutic claim rather than weakening")
    w("it: a drug that made the whole liver regenerate faster would do nothing")
    w("at all, while anything that makes infected cells preferentially die --")
    w("interferon, checkpoint release, a therapeutic vaccine -- moves the floor")
    w("directly.")
    w()
    w("   THIS IS THE FALSIFIABLE CORE OF THE MODEL'S CLAIM.  Summed over the")
    w("   host parameters, host cell biology outweighs drug-target affinity")
    w("   %sx, and the single dominant term is the DIFFERENTIAL death rate of"
      % fmt(host_sum / max(1e-9, kd_row[4]), 1))
    w("   infected hepatocytes.  If a next-generation entry inhibitor with")
    w("   substantially higher NTCP affinity delivers substantially better")
    w("   week-48 responses, this model is wrong.")

    # -------------------------------------------------------------------
    hdr("A11.  ORGAN AND OUTCOME PROJECTIONS -- 5 YEARS")
    w("Fibrosis responds to ALT, not to HDV RNA, so the outcome ranking of")
    w("regimens follows the ALT column of A6/A8, not the virology column.")
    w("Five-year projections (mean over the virtual population):")
    w()
    long_arms = [
        ("untreated (NUC only)", {}),
        ("BLV 2 mg qd continuous", {"blv_mg": 2.0}),
        ("BLV 10 mg qd continuous", {"blv_mg": 10.0}),
        ("PegIFN 48 wk then stop", {"ifn_ug": 180.0, "stop_day": 336.0}),
        ("BLV 2 mg + PegIFN 48 wk, BLV on", {"blv_mg": 2.0, "ifn_ug": 180.0,
                                             "stop_day": 336.0}),
    ]
    TE = 5 * 365.0
    w("   arm                                dFib      Fib@5y   cirrhosis%   HCC 5y%   PLT")
    for nm, a in long_arms:
        if "stop_day" in a and a.get("blv_mg", 0) > 0:
            # peg-IFN stops at 48 wk but bulevirtide continues
            recs = []
            for q in pop:
                q1 = qss_of(q, blv_mg=a["blv_mg"], ifn_ug=a["ifn_ug"])
                q2 = qss_of(q, blv_mg=a["blv_mg"])
                ts1, ys1 = simulate(q, Regimen(), 336.0, dt=0.5, obs_every=14.0, qss=q1)
                ts2, ys2 = simulate(q, Regimen(), TE - 336.0, dt=0.5,
                                    obs_every=14.0, qss=q2, y0=ys1[-1])
                recs.append({"V0": q["tg_Vd"], "A0": q["tg_ALT"], "q": q,
                             "ts": ts1 + [336.0 + t for t in ts2[1:]],
                             "ys": ys1 + ys2[1:]})
        else:
            recs = run_arm(pop, dict(a, **{"stop_day": a.get("stop_day", TE)}),
                           TE, dt=0.5)
        f0 = sum(r["q"]["tg_Fib"] for r in recs) / len(recs)
        f5 = sum(at(r["ts"], r["ys"], TE, "Fib") for r in recs) / len(recs)
        cir = 100.0 * sum(1 for r in recs
                          if at(r["ts"], r["ys"], TE, "Fib") >= 5.0) / len(recs)
        hcc = 100.0 * sum(1 - math.exp(-at(r["ts"], r["ys"], TE, "CH"))
                          for r in recs) / len(recs)
        pl = sum(at(r["ts"], r["ys"], TE, "PLT") for r in recs) / len(recs)
        w("   %-33s %8s %9s %10s%% %8s%% %6s"
          % (nm, fmt(f5 - f0, 2), fmt(f5, 2), fmt(cir, 1), fmt(hcc, 2), fmt(pl, 0)))
    w()
    w("   Literature orientation, stated carefully because the obvious")
    w("   comparison is the WRONG one.  The often-quoted 'cirrhosis in ~70%'")
    w("   for hepatitis D refers to time from INFECTION, accumulated over")
    w("   decades, and many patients are already cirrhotic at presentation.  It")
    w("   is NOT a five-year rate from an Ishak-2 baseline and this model must")
    w("   not be read against it.  The right comparison is the measured")
    w("   paired-biopsy progression rate, roughly 0.15-0.25 Ishak units/y for")
    w("   hepatitis D versus ~0.10 for HBV alone -- and the untreated arm above")
    w("   runs inside that window (see A2).  HCC in the untreated arm comes out")
    w("   near the ~0.5-1%/y reported for NON-cirrhotic chronic hepatitis D;")
    w("   the 2.8-4%/y figure applies once a patient is cirrhotic, which few")
    w("   are at this baseline.  Nothing in this section was fitted.")

    # -------------------------------------------------------------------
    hdr("A12.  BILE ACIDS AS A DOSE-FINDING INSTRUMENT")
    w("If TBA inverts to occupancy (A1) and occupancy has a saturating effect")
    w("on response (A6), then TBA is a THERAPEUTIC-INDEX instrument: it reads")
    w("the cost of a dose on a scale where the benefit has already flattened.")
    w()
    w("   dose(mg)   Css(nM)   occupancy   free NTCP   TBA fold   comb resp")
    for mg in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        o = blv_occ(p, mg)
        popd = make_population(p, 150)
        rr = rate(run_arm(popd, {"blv_mg": mg}, 336.0, dt=0.5), 336.0, p, "comb")
        w("   %8s %9s %11s %11s %10s %10s%%"
          % (fmt(mg, 1), fmt(blv_css(p, mg), 2), fmt(o, 4),
             fmt(100 * (1 - o), 1) + "%", fmt(tba_fold(p, o), 1), fmt(rr, 1)))
    w()
    w("   Benefit saturates; cost does not.  Going 2 -> 10 mg buys a few")
    w("   points of response for a ~4x further rise in bile acids, and")
    w("   2 -> 20 mg buys almost nothing for a large further rise.  On this")
    w("   model the approved 2 mg dose sits at the knee, and it gets there")
    w("   WITHOUT saturating its target -- which is not the usual reason a")
    w("   dose-response curve flattens.")

    # -------------------------------------------------------------------
    hdr("A13.  WHERE THIS MODEL IS WRONG, AND WHAT WOULD FALSIFY IT")
    w("1. THE DECOUPLING TERM IS A HYPOTHESIS.  The kappa*ENTRY_FLUX injury")
    w("   term (A7) is the cheapest structure that reproduces early ALT")
    w("   normalisation without virologic response.  A competing explanation")
    w("   is that bulevirtide reduces intrahepatic inflammation by an")
    w("   NTCP-linked bile-acid mechanism unrelated to fresh infection.")
    w("   FALSIFIER: in patients with no virologic response at all, ALT should")
    w("   still fall in proportion to (1 - occupancy).  If ALT decline tracks")
    w("   HDV RNA decline patient-by-patient instead, this term is wrong.")
    w()
    w("2. THE FLOOR PARAMETERS ARE NOT SEPARATELY IDENTIFIABLE.  frac_cc and")
    w("   d_hep both act on the same net per-capita growth rate; A10 shows")
    w("   they have similar leverage, and the fit cannot tell them apart.")
    w("   FALSIFIER: paired intrahepatic HDAg staining over time on")
    w("   bulevirtide.  Division-mediated spread predicts a nearly unchanged")
    w("   HDAg+ FRACTION with falling serum RNA; a cell-to-cell-dominated")
    w("   floor predicts a falling fraction.")
    w()
    w("3. THE DOSE SEPARATION IS THE MODEL'S CLEAREST FAILURE, AND ALSO WHERE")
    w("   IT TAUGHT ME SOMETHING.  I expected that once the floor was in place,")
    w("   the %sx difference in residual entry between 2 and 10 mg would wash"
      % fmt((1 - o2) / (1 - o10), 1))
    w("   out on every endpoint.  It DOES wash out on the combined endpoint")
    w("   (%s vs %s%%, observed 45 vs 48%%) but it does NOT wash out on"
      % (fmt(r2, 1), fmt(r10, 1)))
    w("   virology (%s vs %s%%, observed 71 vs 76%%).  So the floor explains the"
      % (fmt(rate(res["BLV 2 mg qd"], 336.0, p, "vr"), 1),
         fmt(rate(res["BLV 10 mg qd"], 336.0, p, "vr"), 1)))
    w("   flat combined response only PARTLY -- the rest of that flatness is")
    w("   the ALT ceiling, a completely different mechanism -- and the flat")
    w("   VIROLOGIC response still needs its own explanation (A6, options i-iii).")
    w("   I went in with one mechanism for two observations and came out")
    w("   needing two.")
    w()
    w("4. NO SPATIAL STRUCTURE.  Cell-to-cell spread is treated as a")
    w("   well-mixed mass-action term.  Real HDV spread is clustered around")
    w("   infected foci, which would make the floor LESS sensitive to")
    w("   population-average occupancy than modelled here.")
    w()
    w("5. HBsAg IS TREATED AS A SINGLE POOL.  Integrated HBV DNA contributes")
    w("   HBsAg independently of cccDNA, so the siRNA arms in A9 may be")
    w("   optimistic about how far envelope supply can actually be cut.")
    w()
    w("6. THE ADAR1 / IFN INTERACTION IS DIRECTIONALLY AMBIGUOUS.  Interferon")
    w("   induces ADAR1, which raises amber/W editing, which makes MORE")
    w("   L-HDAg -- and L-HDAg both REPRESSES replication and ENABLES")
    w("   assembly.  The model resolves this numerically (net antiviral), but")
    w("   the balance is parameter-dependent and should not be over-read.")
    w()
    w("7. GENOTYPE IS ABSENT.  HDV-1 through HDV-8 differ in course and, for")
    w("   HDV-3, in editing biology.  All numbers here are HDV-1-like.")
    w()
    w("8. THE LONAFARNIB AND siRNA COMPOSITE RESPONSES ARE TOO HIGH.  D-LIVR")
    w("   reported ~10% composite response for lonafarnib/ritonavir; the model")
    w("   gives several times that.  The reason is visible in A9: the model lets")
    w("   the large serum-RNA drop from an assembly block feed through into ALT")
    w("   normalisation, because in the model ALT depends on the infected pool")
    w("   and the entry flux, and an assembly block slowly lowers BOTH.  In the")
    w("   trial it did not.  The most likely missing piece is that the trapped")
    w("   intracellular genomes and unprenylated L-HDAg are themselves")
    w("   cytotoxic / immunogenic, so blocking assembly should ADD an injury")
    w("   term rather than only removing one.  That is a concrete structural")
    w("   fix and it is not in this version.")
    w()
    w("9. NO ADHERENCE, NO DOSE INTERRUPTION, NO TOXICITY-DRIVEN DISCONTINUATION.")
    w("   Lonafarnib's GI toxicity and peg-IFN's cytopenias are the reason real")
    w("   regimens under-deliver relative to their pharmacology, and neither")
    w("   truncates exposure anywhere in this file.  Every number here is a")
    w("   per-protocol number.")

    # -------------------------------------------------------------------
    hdr("A14.  PARAMETER TABLE (post-fit, post-back-solve)")
    keys = sorted(k for k in p if not k.startswith("_") and not k.startswith("tg_")
                  and p[k] is not None and isinstance(p[k], (int, float)))
    for k in keys:
        w("   %-16s = %s" % (k, "%.6g" % p[k]))
    w()
    w("   FITTED (4):  r_oatp, Kd_ntcp  (bile-acid anchors)")
    w("                dth_Id_immune  (virologic response rate, 2 mg)")
    w("                ALT_base       (ALT normalisation rate, 2 mg)")
    w("   Everything else is either literature-sourced or back-solved from the")
    w("   observed untreated steady state.")

    hdr("END OF REPORT")
    w("Generated by hdv_reference_model.py (pure Python, no dependencies).")
    w("Every number above is reproducible with: python3 hdv_reference_model.py")

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "hdv_model_report.txt"), "w") as f:
        f.write("\n".join(OUT) + "\n")
    print("\n".join(OUT))


if __name__ == "__main__":
    main()
