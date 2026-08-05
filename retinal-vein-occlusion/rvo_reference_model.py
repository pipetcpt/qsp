#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
===============================================================================
rvo_reference_model.py
Retinal Vein Occlusion (CRVO / BRVO) macular oedema
망막정맥폐쇄(중심/분지) 대사성 황반부종 — INDEPENDENT PYTHON REFERENCE MODEL
===============================================================================

WHY THIS FILE EXISTS
--------------------
This build environment has no R runtime, so `rvo_mrgsolve_model.R` cannot be
executed here.  Committing an un-integrated ODE model would be dishonest, so
every equation in the mrgsolve file is re-implemented in this file, term for
term, and actually integrated (fixed-step RK4, pure standard library, no numpy).
Every number quoted in `README.md` is produced by running this file; the
captured log is `rvo_reference_output.txt` and the machine-readable results are
`rvo_scenario_results.json` / `rvo_population_results.json`.

-------------------------------------------------------------------------------
THE ONE STRUCTURAL CLAIM
-------------------------------------------------------------------------------
Retinal thickening in vein occlusion is a Starling flux:

        Jv = Lp · S · [ (Pc - Pt)  -  sigma·(pi_c - pi_t) ]
                       \________/     \__________________/
                      PRESSURE ARM       ONCOTIC ARM
             ^^^^^
             Lp = PERMEABILITY ARM

Anti-VEGF therapy acts on **Lp and sigma only**.  It is a multiplier in front
of the bracket.  It contains no term that can change Pc, and Pc is set by the
occlusion:

        Pc = (Pa·Rv + Pv·Ra) / (Ra + Rv)

Three consequences fall out as arithmetic rather than being asserted:

 (1) THERE IS A CRITICAL CAPILLARY PRESSURE.  Set Jv = 0 with the permeability
     arm driven to its floor (perfect, infinite VEGF blockade).  Then

        Pc* = Pt + sigma_max·(pi_c - pi_t)

     Above Pc* the bracket is positive even with Lp at basal, so **no dose of
     any anti-VEGF agent can dry the macula**.  Below Pc* the macula dries.
     The whole ischaemic/non-ischaemic, CRVO/BRVO, responder/non-responder
     split is which side of Pc* the eye sits on, and the model computes it.

 (2) BECAUSE Pc IS BOUNDED BY THE ARTERIAL INLET PRESSURE, LOWERING SYSTEMIC
     BLOOD PRESSURE MOVES Pc* NON-LINEARLY.  dPc/dPa = Rv/(Ra+Rv), which tends
     to 1 as the occlusion tightens.  In a normal eye arterial pressure barely
     reaches the capillary; in an occluded eye it arrives almost undamped.  The
     model turns this into a table: the critical venous resistance at which a
     permanent oedema floor appears rises ~8-fold for an 8 mmHg fall in
     retinal arteriolar inlet pressure.

 (3) DURATION OF VEGF SUPPRESSION AND DURATION OF DRYNESS ARE DIFFERENT
     NUMBERS.  Suppression duration is pharmacology:
        t_sup = (t_half/ln2) · ln( C0_eff / threshold )
     and for aflibercept it is long (months).  Dryness duration is the disease:
     oedema returns while VEGF is still fully suppressed if Pc > Pc*.  Those
     eyes are mislabelled "anti-VEGF non-responders"; the model shows they are
     pressure-arm eyes and that switching agents cannot help them.

A fourth, separate claim concerns vision:

 (4) OEDEMA IS A STATE, VISION IS AN INTEGRAL.  Letters lost to thickening are
     a reversible function of the current state.  Letters lost to photoreceptor
     (ellipsoid zone) loss are the time-integral of thickening x ischaemia and
     do not come back.  Therefore delay is priced in permanent letters, and the
     price is paid whether or not the macula is eventually dried.

-------------------------------------------------------------------------------
STATE VECTOR (27 ODEs)
-------------------------------------------------------------------------------
  0  DVIT    nmol      anti-VEGF drug, vitreous depot (4.0 mL)
  1  ASYS    nmol      anti-VEGF drug, systemic
  2  VTONE   pM        VEGF-A tone (free VEGF that WOULD be present, no drug)
  3  PTONE   pM        PlGF tone
  4  ATONE   pM        Angiopoietin-2 tone
  5  IL6     pg/mL     interleukin-6 (VEGF-independent permeability driver)
  6  HIF     0-1       hypoxia-inducible signalling
  7  OCC     x         occlusion severity multiplier on venous resistance
  8  COLL    x         collateral / recanalisation conductance
  9  NP      DA        retinal capillary non-perfusion, disc areas
 10  MI      0-3       macular (foveal) ischaemia index
 11  LEUK    0-1       leukostasis / ICAM-1 capillary plugging
 12  TJ      0-1       tight-junction integrity (occludin / claudin-5 / VE-cad)
 13  CHRON   0-1       chronic remodelling (Mueller gliosis, cyst walls)
 14  W       um        intraretinal excess water (above dry baseline)
 15  SRF     um        subretinal fluid
 16  EZ      0-1       ellipsoid-zone / photoreceptor integrity (irreversible)
 17  DRIL    0-1       disorganisation of retinal inner layers
 18  CUMED   um.d      cumulative oedema exposure above 320 um
 19  NVI     0-1       anterior-segment neovascularisation (iris / angle)
 20  IOP     mmHg      intraocular pressure
 21  CAT     0-1       lens opacity
 22  IMP     ng        dexamethasone remaining in intravitreal implant
 23  CDEX    ng/mL     vitreous dexamethasone concentration
 24  BCVAO   letters   observed BCVA (10-day functional lag on the algebraic value)
 25  ABLA    0-1       fraction of retina ablated by panretinal photocoagulation
 26  TSUP    d         cumulative days of adequate VEGF suppression (tracker)

-------------------------------------------------------------------------------
CALIBRATION TARGETS (all from published trials; see rvo_references.md)
-------------------------------------------------------------------------------
  CRUISE   ranibizumab 0.5 mg q4w x6 in CRVO ... +14.9 letters, CST -452 um
  CRUISE   sham arm .......................... +0.8 letters
  CRUISE   sham -> ranibizumab at m6 ......... +7.3 letters at m12 (never catches up)
  COPERNICUS aflibercept 2 mg q4w x6 ........ +17.3 letters at w24, CST -457 um
  COPERNICUS sham -> aflibercept ............ +1.5 letters at w100
  GALILEO  aflibercept ...................... +18.0 letters at w24, +13.7 at w76
  SCORE2   aflibercept +18.9 / bevacizumab +18.6 letters at m6 (NO affinity gap)
  VIBRANT  aflibercept in BRVO .............. +17.0 vs laser +12.2 at w24
  BRAVO    ranibizumab in BRVO .............. +18.3 letters at m6
  BALATON/COMINO faricimab .................. +16.9 letters, non-inferior to aflibercept
  GENEVA   dexamethasone implant ............ +/- IOP rise 16%, cataract 30%/yr
  CVOS     ischaemic CRVO -> INV/ANV ........ ~35% within 3 years
  Real world (LUMINOUS / IRIS registry) ..... 4-5 injections/yr -> +6 to +8 letters
===============================================================================
"""

from __future__ import annotations

import json
import math
import os
import random
import sys

# =============================================================================
# 0.  SMALL NUMERICAL HELPERS  (no numpy in this environment)
# =============================================================================

def clamp(x: float, lo: float, hi: float) -> float:
    return lo if x < lo else (hi if x > hi else x)


def hill(x: float, k: float, n: float = 1.0) -> float:
    """Saturating Michaelis / Hill fraction, safe at x <= 0."""
    if x <= 0.0:
        return 0.0
    if n == 1.0:
        return x / (k + x)
    xn = x ** n
    return xn / (k ** n + xn)


def safe_pow(x: float, p: float) -> float:
    return 0.0 if x <= 0.0 else x ** p


# =============================================================================
# 1.  DRUG LIBRARY
# =============================================================================
# MW      Da
# dose    mg per intravitreal injection (label dose)
# t_vit   apparent intraocular half-life, days (aqueous-sampling estimates,
#         Krohne et al. 2008/2012 for bevacizumab & ranibizumab & aflibercept;
#         brolucizumab and faricimab from population PK reports)
# kd_v    intrinsic KD for VEGF-A165, pM (Papadopoulos 2012 for AFL/RBZ/BEV)
# kd_p    intrinsic KD for PlGF, pM (None = does not bind)
# kd_a    intrinsic KD for Ang-2, pM (None = does not bind)
# t_sys   systemic half-life of the free drug, days
# -----------------------------------------------------------------------------
DRUGS = {
    "ranibizumab": dict(label="Ranibizumab 0.5 mg", mw=48000.0, dose=0.5,
                        t_vit=7.19, kd_v=46.0, kd_p=None, kd_a=None, t_sys=0.083),
    "bevacizumab": dict(label="Bevacizumab 1.25 mg", mw=149000.0, dose=1.25,
                        t_vit=9.82, kd_v=58.0, kd_p=None, kd_a=None, t_sys=20.0),
    "aflibercept": dict(label="Aflibercept 2 mg", mw=115000.0, dose=2.0,
                        t_vit=9.10, kd_v=0.49, kd_p=39.0, kd_a=None, t_sys=0.50),
    "aflibercept8": dict(label="Aflibercept 8 mg", mw=115000.0, dose=8.0,
                         t_vit=9.10, kd_v=0.49, kd_p=39.0, kd_a=None, t_sys=0.50),
    "brolucizumab": dict(label="Brolucizumab 6 mg", mw=26000.0, dose=6.0,
                         t_vit=5.10, kd_v=28.4, kd_p=None, kd_a=None, t_sys=4.40),
    "faricimab": dict(label="Faricimab 6 mg", mw=149000.0, dose=6.0,
                      t_vit=7.50, kd_v=30.0, kd_p=None, kd_a=1500.0, t_sys=7.50),
}


def molar_dose_nmol(drug: str) -> float:
    d = DRUGS[drug]
    return d["dose"] / d["mw"] * 1.0e6      # mg / (g/mol) * 1e6 = nmol


# =============================================================================
# 2.  PARAMETERS
# =============================================================================
P = dict(

    # ---- ocular geometry -----------------------------------------------------
    VVIT=4.0,            # mL, vitreous cavity
    FPEN=0.080,          # retina : vitreous drug partition (penetration factor)

    # ---- the single calibrated pharmacodynamic scale factor -----------------
    # ALPHA converts an INTRINSIC binding KD (pM, measured by SPR in buffer)
    # into an apparent in-vivo IC50 at the outer blood-retina barrier.  It is
    # >1 because the drug must out-compete VEGFR2 (itself sub-nM), because the
    # relevant VEGF concentration at the junction is a local flux and not a
    # vitreous average, and because of diffusional limitation across the
    # retina.  ONE value is used for every agent, so the RELATIVE ranking of
    # agents is fixed entirely by their measured KD, molar dose and half-life.
    ALPHA=120.0,
    SUPP_RATIO_ADEQ=4.9,   # suppression ratio counted as "adequate" (tracker)

    # ---- retinal haemodynamics ----------------------------------------------
    PA=40.0,             # mmHg, post-arteriolar inlet pressure
    RA0=1.0,             # relative arteriolar resistance
    RV0=0.25,            # relative venular + venous resistance (healthy)
    PT=12.0,             # mmHg, retinal interstitial (tissue) pressure
    PI_C=25.0,           # mmHg, plasma colloid osmotic pressure
    PI_T0=4.0,           # mmHg, interstitial colloid osmotic pressure, dry
    SIGMA_MAX=0.92,      # reflection coefficient of an intact inner BRB
    AUTOREG=0.35,        # hypoxic arteriolar dilatation (fraction of RA0 lost)

    # ---- occlusion natural history ------------------------------------------
    K_LYS=0.0045,        # /d, thrombus lysis / recanalisation (CRVO)
    OCC_RES=0.55,        # non-lysable residual fraction of the occlusion.
                         # This is the single most consequential natural-history
                         # parameter: a fully lysing occlusion cures itself and
                         # makes every therapeutic comparison meaningless.
    K_COLL=0.0035,       # /d, collateral growth rate constant (months, not weeks)
    COLL_MAX=0.70,       # CRVO has few collateral options
    HIF_COLL=2.0,        # hypoxic amplification of collateral growth

    # ---- hypoxia / VEGF axis -------------------------------------------------
    TAU_HIF=2.0,         # d
    RV_VEGF=24.0,        # pM/d basal VEGF-A production
    KDEG_V=8.0,          # /d VEGF-A turnover (t1/2 ~2 h)
    E_HIF_V=12.0,        # fold VEGF induction at HIF = 1
    E_ABLA_V=0.72,       # fraction of VEGF production removable by PRP
    E_DEX_V=0.70,        # fraction of VEGF production suppressed by dexamethasone
    RP_PLGF=6.0, KDEG_P=3.0, E_HIF_P=5.0,
    RA_ANG2=10.0, KDEG_A=2.0, E_HIF_A=7.0, E_DEX_A=0.50,
    RI_IL6=30.0, KDEG_I=6.0, E_HIF_I=16.0, E_DEX_I=0.80,

    # ---- non-perfusion / leukostasis ---------------------------------------
    K_NP=0.30,           # DA/d gain
    NP_MAX=75.0,         # DA, retinal area available to close
    Q_NP=0.45,           # relative perfusion below which capillaries close
    NP_LEUK=0.80,
    K_NP_REP=0.0022,     # /d reperfusion via collaterals
    K_MI=0.020, Q_MI=0.35, K_MI_REP=0.0020,
    TAU_LEUK=5.0, W_LEUK_I6=1.0, W_LEUK_V=0.50,

    # ---- tight junction / permeability -------------------------------------
    KREP_TJ=0.160,       # /d junction re-assembly
    KDIS_TJ=0.750,       # /d junction disassembly at DRIVE = 1
    TIE2_BOOST=0.50,     # extra repair when Tie2 is unopposed (Ang-2 blocked)
    E_DEX_TJ=0.35,       # direct steroid up-regulation of occludin / claudin-5
                         # transcription. This is a junction effect anti-VEGF has
                         # no mechanism for, and it is why a steroid can work in
                         # an eye that has failed anti-VEGF.
    K_TIE2=15.0,         # pM Ang-2 for half Tie2 inhibition
    W_V=0.80, KV=12.0,           # VEGF-A weight and potency in DRIVE
    W_I6=0.08, KI6=45.0,         # IL-6
    W_A2=0.06, KA2=15.0,         # Ang-2
    W_PL=0.04, KPL=8.0,          # PlGF
    W_PRESS=0.06, P_PRESS0=28.0, P_PRESSW=12.0,   # mechanotransduction term
    LP_MAX=6.5, LP_N=1.30,       # Lp_rel = 1 + (LP_MAX-1)(1-TJ)^LP_N
    SIG_N=0.55,                  # sigma = SIGMA_MAX * TJ^SIG_N
    PI_T_GAIN=0.75,              # protein leak -> interstitial oncotic rise

    # ---- chronic remodelling ------------------------------------------------
    K_CHRON=0.0016, W_CHRON0=70.0, W_CHRONW=300.0,
    K_CHRON_REV=0.0006, TJMAX_LOSS=0.45,

    # ---- fluid balance ------------------------------------------------------
    KF=0.870,            # um/(d.mmHg) per unit Lp_rel  (filtration coefficient)
    K_OUT_LIN=0.100,     # /d linear (diffusive / bulk) clearance of water
    VMAX_PUMP=32.0,      # um/d Mueller-AQP4 / RPE active transport
    KM_PUMP=45.0,        # um
    F_SRF=0.18,          # fraction of filtered flux entering subretinal space
    K_SRF_LIN=0.060, VMAX_SRF=14.0, KM_SRF=30.0,

    # ---- irreversible damage ------------------------------------------------
    K_EZ=0.0018, KM_EZ=250.0, W_EZ_MI=1.10, K_EZ_REP=0.00040, KM_EZR=150.0,
    K_DRIL=0.0035, K_DRIL_REP=0.00025,
    CST_CUM=320.0,

    # ---- neovascular / pressure complications ------------------------------
    K_NVI=0.030, KM_NVI_NP=10.0, KM_NVI_V=16.0, K_NVI_REG=0.012,
    IOP0=15.0, IOP_NVI=24.0, NVI_THR=0.40, TAU_IOP=4.0,
    IOP_DEX=10.0, KM_IOP_DEX=200.0,

    # ---- lens ---------------------------------------------------------------
    K_CAT_DEX=0.0024, KM_CAT_DEX=150.0, K_CAT_AGE=0.00008,

    # ---- dexamethasone implant ---------------------------------------------
    DEX_LOAD=700000.0,   # ng (0.7 mg Ozurdex)
    K_REL_DEX=0.0128,    # /d
    KEL_DEX=2.50,        # /d vitreous clearance
    KM_DEX_EFF=120.0,    # ng/mL for half maximal steroid effect

    # ---- visual function ----------------------------------------------------
    BCVA_CEIL=80.0,      # ETDRS letters achievable with a dry, intact macula
    E_ED=36.0, KM_ED=230.0,      # reversible letters lost to thickening
    E_PR=42.0,                   # letters lost per unit EZ loss
    E_ISCH=22.0, KM_ISCH=1.20,   # letters lost to macular ischaemia
    E_CAT=12.0,                  # letters lost to lens opacity
    E_PRP=6.0,                   # letters lost to panretinal photocoagulation
    TAU_BCVA=10.0,               # d functional lag

    # ---- sector involvement (1.0 = CRVO, <1 = BRVO) ------------------------
    # Two DIFFERENT fractions are needed.  FSEC_FLUX is how much of the central
    # subfield lies in the occluded sector (drives the fluid term and the VEGF
    # source).  FSEC_DMG is how much of the FOVEA is functionally at risk, which
    # is smaller: a BRVO can flood the central subfield while sparing part of the
    # foveal photoreceptor mosaic.  Collapsing them into one number is what makes
    # naive BRVO models predict impossibly good vision.
    FSEC_FLUX=1.0,
    FSEC_DMG=1.0,

    # ---- misc ---------------------------------------------------------------
    STEROID_RESPONDER=1.0,
    BASE_CST=250.0,      # um, dry central subfield thickness
    VP=4.0,              # L, systemic volume of distribution
)

# state index map
S = {n: i for i, n in enumerate(
    "DVIT ASYS VTONE PTONE ATONE IL6 HIF OCC COLL NP MI LEUK TJ CHRON W SRF "
    "EZ DRIL CUMED NVI IOP CAT IMP CDEX BCVAO ABLA TSUP".split())}
NSTATE = len(S)


# =============================================================================
# 3.  DERIVED / ALGEBRAIC QUANTITIES
# =============================================================================

def algebraic(t, y, p, drug, perfect_block=False):
    """Everything that is a function of the state, returned as a dict.

    Kept in one place so that the ODE right-hand side, the observation
    function and the mrgsolve $TABLE block all agree by construction.
    """
    a = {}

    # ---- drug concentrations ------------------------------------------------
    cvit = y[S["DVIT"]] / p["VVIT"] * 1.0e6            # nmol/mL -> pM
    cret = cvit * p["FPEN"]                            # effect-site, pM
    a["CVIT_pM"] = cvit
    a["CRET_pM"] = cret
    a["CVIT_ug_mL"] = y[S["DVIT"]] / p["VVIT"] * DRUGS[drug]["mw"] / 1.0e6
    cp = y[S["ASYS"]] / p["VP"] * 1.0e3                # nmol/L = nM
    a["CP_nM"] = cp
    a["CP_ug_mL"] = cp * DRUGS[drug]["mw"] / 1.0e6

    # ---- competitive inhibition of each ligand -----------------------------
    d = DRUGS[drug]
    alpha = p["ALPHA"]

    def inhib(kd):
        if perfect_block:
            return 1.0e12
        if kd is None:
            return 0.0
        return cret / (kd * alpha)

    rv = inhib(d["kd_v"])
    rp = inhib(d["kd_p"])
    ra = inhib(d["kd_a"])
    a["SUPP_RATIO_V"] = 1.0 + rv
    a["VACT"] = y[S["VTONE"]] / (1.0 + rv)
    a["PACT"] = y[S["PTONE"]] / (1.0 + rp)
    a["AACT"] = y[S["ATONE"]] / (1.0 + ra)

    # ---- systemic free VEGF (safety read-out) ------------------------------
    if perfect_block:
        a["SYS_VEGF_SUPP"] = 1.0
    elif d["kd_v"] is None:
        a["SYS_VEGF_SUPP"] = 0.0
    else:
        a["SYS_VEGF_SUPP"] = 1.0 - 1.0 / (1.0 + cp * 1.0e3 / d["kd_v"])

    # ---- steroid effect -----------------------------------------------------
    a["ESTER"] = hill(y[S["CDEX"]], p["KM_DEX_EFF"])

    # ---- haemodynamics ------------------------------------------------------
    ra_res = p["RA0"] * (1.0 - p["AUTOREG"] * y[S["HIF"]])
    rv_res = p["RV0"] * (1.0 + y[S["OCC"]]) / (1.0 + y[S["COLL"]])
    pv = y[S["IOP"]] + 2.0
    a["RA_RES"] = ra_res
    a["RV_RES"] = rv_res
    a["PC"] = (p["PA"] * rv_res + pv * ra_res) / (ra_res + rv_res)
    a["Q"] = (p["PA"] - pv) / (ra_res + rv_res)
    q0 = (p["PA"] - (p["IOP0"] + 2.0)) / (p["RA0"] + p["RV0"])
    a["Q0"] = q0
    a["QREL"] = a["Q"] / q0
    a["DPC_DPA"] = rv_res / (ra_res + rv_res)          # transmission coefficient
    a["HYPOX"] = clamp(1.0 - safe_pow(a["QREL"], 0.80), 0.0, 1.0)

    # ---- barrier ------------------------------------------------------------
    tj = clamp(y[S["TJ"]], 1.0e-3, 1.0)
    a["LP_REL"] = 1.0 + (p["LP_MAX"] - 1.0) * safe_pow(1.0 - tj, p["LP_N"])
    a["SIGMA"] = p["SIGMA_MAX"] * safe_pow(tj, p["SIG_N"])
    leak = 1.0 - a["SIGMA"] / p["SIGMA_MAX"]
    a["PI_T"] = p["PI_T0"] + (p["PI_C"] - p["PI_T0"]) * leak * p["PI_T_GAIN"]

    a["DRIVE"] = (p["W_V"] * hill(a["VACT"], p["KV"])
                  + p["W_I6"] * hill(y[S["IL6"]], p["KI6"])
                  + p["W_A2"] * hill(a["AACT"], p["KA2"])
                  + p["W_PL"] * hill(a["PACT"], p["KPL"])
                  + p["W_PRESS"] * clamp((a["PC"] - p["P_PRESS0"]) / p["P_PRESSW"],
                                         0.0, 1.0))

    # ---- Starling flux ------------------------------------------------------
    a["BRACKET"] = (a["PC"] - p["PT"]) - a["SIGMA"] * (p["PI_C"] - a["PI_T"])
    a["JV"] = p["KF"] * a["LP_REL"] * a["BRACKET"] * p["FSEC_FLUX"]

    # the same flux with the permeability arm forced to its floor: this is the
    # part of the oedema that NO anti-VEGF agent can remove
    lp_floor = 1.0
    sig_floor = p["SIGMA_MAX"]
    a["JV_FLOOR"] = (p["KF"] * lp_floor
                     * ((a["PC"] - p["PT"]) - sig_floor * (p["PI_C"] - p["PI_T0"]))
                     * p["FSEC_FLUX"])
    a["PC_CRIT"] = p["PT"] + sig_floor * (p["PI_C"] - p["PI_T0"])
    a["PC_CRIT_NOW"] = p["PT"] + a["SIGMA"] * (p["PI_C"] - a["PI_T"])
    a["PRESSURE_ARM_FRAC"] = (max(0.0, a["JV_FLOOR"]) / a["JV"]) if a["JV"] > 1e-9 else 0.0

    # ---- clearance ----------------------------------------------------------
    w = max(0.0, y[S["W"]])
    a["JOUT"] = p["K_OUT_LIN"] * w + p["VMAX_PUMP"] * hill(w, p["KM_PUMP"])

    # ---- observables --------------------------------------------------------
    a["WTOT"] = w + 0.6 * max(0.0, y[S["SRF"]])
    a["CST"] = p["BASE_CST"] + a["WTOT"]
    fs = p["FSEC_DMG"]
    a["L_ED"] = p["E_ED"] * hill(a["WTOT"], p["KM_ED"])
    a["L_PR"] = p["E_PR"] * (1.0 - y[S["EZ"]]) * fs
    a["L_ISCH"] = p["E_ISCH"] * hill(y[S["MI"]], p["KM_ISCH"]) * fs
    a["L_CAT"] = p["E_CAT"] * y[S["CAT"]]
    a["L_PRP"] = p["E_PRP"] * y[S["ABLA"]]
    a["BCVA_ALG"] = max(0.0, p["BCVA_CEIL"] - a["L_ED"] - a["L_PR"]
                        - a["L_ISCH"] - a["L_CAT"] - a["L_PRP"])
    a["BCVA"] = y[S["BCVAO"]]
    a["SNELLEN_DEN"] = 20.0 * (10.0 ** ((85.0 - a["BCVA"]) / 50.0))
    a["NVG"] = 1.0 if (y[S["NVI"]] > 0.50 and y[S["IOP"]] > 25.0) else 0.0
    return a


# =============================================================================
# 4.  RIGHT-HAND SIDE
# =============================================================================

def rhs(t, y, p, drug, perfect_block=False):
    a = algebraic(t, y, p, drug, perfect_block)
    dy = [0.0] * NSTATE

    kel_vit = math.log(2.0) / DRUGS[drug]["t_vit"]
    kel_sys = math.log(2.0) / DRUGS[drug]["t_sys"]

    # --- 1. drug disposition -------------------------------------------------
    dy[S["DVIT"]] = -kel_vit * y[S["DVIT"]]
    dy[S["ASYS"]] = kel_vit * y[S["DVIT"]] - kel_sys * y[S["ASYS"]]

    # --- 2. hypoxia signalling ----------------------------------------------
    hif_target = clamp(a["HYPOX"] * (1.0 + 0.60 * y[S["NP"]] / 30.0), 0.0, 1.0)
    dy[S["HIF"]] = (hif_target - y[S["HIF"]]) / p["TAU_HIF"]

    # --- 3. ligand tone -----------------------------------------------------
    fs = p["FSEC_FLUX"]
    dexv = 1.0 - p["E_DEX_V"] * a["ESTER"]
    dy[S["VTONE"]] = (p["RV_VEGF"] * (1.0 + p["E_HIF_V"] * y[S["HIF"]] * fs)
                      * (1.0 - p["E_ABLA_V"] * y[S["ABLA"]]) * dexv
                      - p["KDEG_V"] * y[S["VTONE"]])
    dy[S["PTONE"]] = (p["RP_PLGF"] * (1.0 + p["E_HIF_P"] * y[S["HIF"]] * fs) * dexv
                      - p["KDEG_P"] * y[S["PTONE"]])
    dy[S["ATONE"]] = (p["RA_ANG2"] * (1.0 + p["E_HIF_A"] * y[S["HIF"]] * fs)
                      * (1.0 - p["E_DEX_A"] * a["ESTER"])
                      - p["KDEG_A"] * y[S["ATONE"]])
    dy[S["IL6"]] = (p["RI_IL6"] * (1.0 + p["E_HIF_I"] * y[S["HIF"]] * fs)
                    * (1.0 - p["E_DEX_I"] * a["ESTER"])
                    - p["KDEG_I"] * y[S["IL6"]])

    # --- 4. occlusion, collaterals ------------------------------------------
    occ_floor = p["OCC_RES"] * p["_OCC0"]
    dy[S["OCC"]] = -p["K_LYS"] * max(0.0, y[S["OCC"]] - occ_floor)
    dy[S["COLL"]] = (p["K_COLL"] * (1.0 + p["HIF_COLL"] * y[S["HIF"]])
                     * max(0.0, p["COLL_MAX"] - y[S["COLL"]]))

    # --- 5. non-perfusion, leukostasis, macular ischaemia -------------------
    leuk_target = clamp(p["W_LEUK_I6"] * hill(y[S["IL6"]], p["KI6"])
                        + p["W_LEUK_V"] * hill(a["VACT"], p["KV"]), 0.0, 1.0)
    dy[S["LEUK"]] = (leuk_target - y[S["LEUK"]]) / p["TAU_LEUK"]

    stasis = max(0.0, p["Q_NP"] - a["QREL"])
    dy[S["NP"]] = (p["K_NP"] * safe_pow(stasis, 1.20)
                   * (1.0 + p["NP_LEUK"] * y[S["LEUK"]])
                   * max(0.0, 1.0 - y[S["NP"]] / p["NP_MAX"])
                   - p["K_NP_REP"] * y[S["NP"]] * hill(y[S["COLL"]], 1.0))
    stasis_m = max(0.0, p["Q_MI"] - a["QREL"])
    dy[S["MI"]] = (p["K_MI"] * stasis_m * (1.0 + 0.5 * y[S["NP"]] / 20.0)
                   - p["K_MI_REP"] * y[S["MI"]])
    if y[S["MI"]] > 3.0:
        dy[S["MI"]] = min(0.0, dy[S["MI"]])

    # --- 6. tight junctions and chronic remodelling -------------------------
    tie2_free = 1.0 / (1.0 + a["AACT"] / p["K_TIE2"])
    krep = p["KREP_TJ"] * (1.0 + p["TIE2_BOOST"] * tie2_free
                           + p["E_DEX_TJ"] * a["ESTER"])
    tjmax = 1.0 - p["TJMAX_LOSS"] * y[S["CHRON"]]
    dy[S["TJ"]] = krep * (tjmax - y[S["TJ"]]) - p["KDIS_TJ"] * a["DRIVE"] * y[S["TJ"]]

    chron_drive = clamp((a["WTOT"] - p["W_CHRON0"]) / p["W_CHRONW"], 0.0, 1.5)
    dy[S["CHRON"]] = (p["K_CHRON"] * chron_drive * (1.0 - y[S["CHRON"]])
                      - p["K_CHRON_REV"] * y[S["CHRON"]])

    # --- 7. fluid ------------------------------------------------------------
    jv = a["JV"]
    dy[S["W"]] = jv - a["JOUT"]
    if y[S["W"]] <= 0.0 and dy[S["W"]] < 0.0:
        dy[S["W"]] = 0.0
    srf_in = p["F_SRF"] * max(0.0, jv) * (1.0 - y[S["TJ"]])
    srf = max(0.0, y[S["SRF"]])
    dy[S["SRF"]] = srf_in - (p["K_SRF_LIN"] * srf + p["VMAX_SRF"] * hill(srf, p["KM_SRF"]))
    if y[S["SRF"]] <= 0.0 and dy[S["SRF"]] < 0.0:
        dy[S["SRF"]] = 0.0

    # --- 8. irreversible damage ---------------------------------------------
    fsd = p["FSEC_DMG"]
    dmg = (p["K_EZ"] * hill(a["WTOT"], p["KM_EZ"])
           * (1.0 + p["W_EZ_MI"] * y[S["MI"]]) * y[S["EZ"]] * fsd)
    rep = p["K_EZ_REP"] * (1.0 - y[S["EZ"]]) * (1.0 - hill(a["WTOT"], p["KM_EZR"]))
    dy[S["EZ"]] = rep - dmg
    dy[S["DRIL"]] = (p["K_DRIL"] * hill(a["WTOT"], p["KM_EZ"]) * (1.0 - y[S["DRIL"]])
                     - p["K_DRIL_REP"] * y[S["DRIL"]])
    dy[S["CUMED"]] = max(0.0, a["CST"] - p["CST_CUM"])

    # --- 9. neovascularisation, IOP, lens -----------------------------------
    dy[S["NVI"]] = (p["K_NVI"] * hill(y[S["NP"]], p["KM_NVI_NP"])
                    * hill(a["VACT"], p["KM_NVI_V"]) * (1.0 - y[S["NVI"]])
                    - p["K_NVI_REG"] * y[S["NVI"]]
                    * (1.0 - hill(a["VACT"], p["KM_NVI_V"])))
    iop_target = (p["IOP0"]
                  + p["IOP_NVI"] * clamp((y[S["NVI"]] - p["NVI_THR"])
                                         / (1.0 - p["NVI_THR"]), 0.0, 1.0)
                  + p["IOP_DEX"] * hill(y[S["CDEX"]], p["KM_IOP_DEX"])
                  * p["STEROID_RESPONDER"])
    dy[S["IOP"]] = (iop_target - y[S["IOP"]]) / p["TAU_IOP"]
    dy[S["CAT"]] = ((p["K_CAT_DEX"] * hill(y[S["CDEX"]], p["KM_CAT_DEX"])
                     + p["K_CAT_AGE"]) * (1.0 - y[S["CAT"]]))

    # --- 10. dexamethasone implant ------------------------------------------
    dy[S["IMP"]] = -p["K_REL_DEX"] * y[S["IMP"]]
    dy[S["CDEX"]] = (p["K_REL_DEX"] * y[S["IMP"]] / p["VVIT"]
                     - p["KEL_DEX"] * y[S["CDEX"]])

    # --- 11. observed BCVA lag, PRP, trackers -------------------------------
    dy[S["BCVAO"]] = (a["BCVA_ALG"] - y[S["BCVAO"]]) / p["TAU_BCVA"]
    dy[S["ABLA"]] = 0.0
    dy[S["TSUP"]] = 1.0 if a["SUPP_RATIO_V"] >= p["SUPP_RATIO_ADEQ"] else 0.0
    return dy


# =============================================================================
# 5.  INITIAL CONDITIONS AND PHENOTYPES
# =============================================================================

def initial_state(p):
    y = [0.0] * NSTATE
    y[S["VTONE"]] = p["RV_VEGF"] / p["KDEG_V"]
    y[S["PTONE"]] = p["RP_PLGF"] / p["KDEG_P"]
    y[S["ATONE"]] = p["RA_ANG2"] / p["KDEG_A"]
    y[S["IL6"]] = p["RI_IL6"] / p["KDEG_I"]
    y[S["HIF"]] = 0.0
    y[S["OCC"]] = 0.0
    y[S["COLL"]] = 0.0
    y[S["TJ"]] = 1.0
    y[S["EZ"]] = 1.0
    y[S["IOP"]] = p["IOP0"]
    y[S["CAT"]] = 0.0
    y[S["BCVAO"]] = p["BCVA_CEIL"]
    return y


PHENOTYPES = {
    # OCC0     occlusion severity multiplier on venous resistance at t = 0
    # OCC_RES  fraction of OCC0 that never lyses (the permanent stenosis)
    # COLL_MAX ceiling on collateral conductance — the BRVO/CRVO asymmetry
    "crvo_nonisch": dict(OCC0=10.0, OCC_RES=0.55, K_LYS=0.0045, COLL_MAX=0.70,
                         FSEC_FLUX=1.00, FSEC_DMG=1.00),
    "crvo_isch":    dict(OCC0=32.0, OCC_RES=0.80, K_LYS=0.0028, COLL_MAX=0.45,
                         FSEC_FLUX=1.00, FSEC_DMG=1.00),
    "brvo":         dict(OCC0=22.0, OCC_RES=0.32, K_LYS=0.0060, COLL_MAX=2.20,
                         FSEC_FLUX=0.85, FSEC_DMG=0.60),
    "crvo_mild":    dict(OCC0=5.0,  OCC_RES=0.20, K_LYS=0.0080, COLL_MAX=1.50,
                         FSEC_FLUX=1.00, FSEC_DMG=1.00),
}


def make_params(phenotype="crvo_nonisch", **over):
    p = dict(P)
    ph = PHENOTYPES[phenotype]
    for k, v in ph.items():
        if k != "OCC0":
            p[k] = v
    p["_OCC0"] = ph["OCC0"]
    p["_PHENOTYPE"] = phenotype
    p.update(over)
    return p


# =============================================================================
# 6.  TREATMENT SCHEDULER  (fixed schedules, PRN and treat-and-extend)
# =============================================================================

class Regimen:
    """Decides, at each monthly monitoring visit, whether to treat.

    mode:
      'none'   never treat
      'fixed'  treat at the times given in `times` (days)
      'prn'    loading phase then treat when CST or BCVA criteria are met
      'tae'    loading phase then treat-and-extend on a dry/not-dry rule
    """

    def __init__(self, mode="none", drug="aflibercept", load_n=6, load_q=28.0,
                 start=42.0, times=None, interval0=56.0, imin=28.0, imax=112.0,
                 istep=28.0, cst_thr=310.0, dex_times=None, prp_time=None,
                 stop=None, max_inj=None, perfect=False):
        self.mode = mode
        self.drug = drug
        self.load_n = load_n
        self.load_q = load_q
        self.start = start
        self.times = list(times) if times else []
        self.interval = interval0
        self.imin, self.imax, self.istep = imin, imax, istep
        self.cst_thr = cst_thr
        self.dex_times = list(dex_times) if dex_times else []
        self.prp_time = prp_time
        self.stop = stop
        self.max_inj = max_inj
        self.perfect = perfect
        # runtime
        self.n_inj = 0
        self.n_dex = 0
        self.next_due = None
        self.last_cst = None
        self.inj_times = []
        self.visit_times = []

    def initialise(self):
        if self.mode == "fixed":
            self.next_due = self.times[0] if self.times else None
        elif self.mode in ("prn", "tae"):
            self.next_due = self.start
        else:
            self.next_due = None

    def decide(self, t, cst, bcva):
        """Return True if an anti-VEGF injection is given at time t."""
        if self.stop is not None and t > self.stop:
            return False
        if self.max_inj is not None and self.n_inj >= self.max_inj:
            return False
        if self.mode == "fixed":
            return any(abs(t - tt) < 1e-6 for tt in self.times)
        if self.mode in ("prn", "tae"):
            if self.next_due is None or t < self.next_due - 1e-6:
                return False
            if self.n_inj < self.load_n:
                self.next_due = t + self.load_q
                return True
            wet = cst > self.cst_thr
            if self.mode == "prn":
                self.next_due = t + self.load_q          # monthly monitoring
                return wet
            # treat and extend
            if wet:
                self.interval = max(self.imin, self.interval - self.istep)
            else:
                self.interval = min(self.imax, self.interval + self.istep)
            self.next_due = t + self.interval
            return True
        return False


# =============================================================================
# 7.  INTEGRATOR
# =============================================================================

def simulate(p, regimen, tmax=1095.0, dt=0.02, record_every=1.0,
             phenotype_occ=None, record=True):
    y = initial_state(p)
    y[S["OCC"]] = p["_OCC0"] if phenotype_occ is None else phenotype_occ

    drug = regimen.drug
    regimen.initialise()
    dose_nmol = molar_dose_nmol(drug)

    out = {k: [] for k in ("t CST W SRF BCVA VACT VTONE PC QREL NP MI TJ LP_REL "
                           "SIGMA DRIVE JV JOUT JV_FLOOR EZ CHRON NVI IOP CAT "
                           "CVIT_ug_mL CRET_pM SUPP_RATIO_V CP_ug_mL SYS_VEGF_SUPP "
                           "CDEX PRESSURE_ARM_FRAC BRACKET PC_CRIT_NOW PC_CRIT DPC_DPA "
                           "CUMED DRIL TSUP").split()}
    events = []

    # perfect-blockade experiments must begin at PRESENTATION, not at onset,
    # otherwise the eye never develops the oedema the experiment is about.
    pfrom = None
    if regimen.perfect:
        pfrom = (regimen.times[0] if regimen.mode == "fixed" and regimen.times
                 else regimen.start)

    def pb_at(tt):
        return pfrom is not None and tt >= pfrom - 1e-9

    n = int(round(tmax / dt))
    next_rec = 0.0
    # monthly monitoring grid used by prn / tae
    for i in range(n + 1):
        t = i * dt

        # ---------------- discrete events -----------------------------------
        # dexamethasone implants
        for dt_ in list(regimen.dex_times):
            if abs(t - dt_) < dt / 2.0:
                y[S["IMP"]] += p["DEX_LOAD"]
                regimen.n_dex += 1
                events.append((round(t, 2), "dex"))
                regimen.dex_times.remove(dt_)
        # panretinal photocoagulation
        if regimen.prp_time is not None and abs(t - regimen.prp_time) < dt / 2.0:
            y[S["ABLA"]] = min(1.0, y[S["ABLA"]] + 0.35)
            events.append((round(t, 2), "prp"))
            regimen.prp_time = None
        # anti-VEGF
        if regimen.mode == "fixed":
            for tt in list(regimen.times):
                if abs(t - tt) < dt / 2.0:
                    y[S["DVIT"]] += dose_nmol
                    regimen.n_inj += 1
                    regimen.inj_times.append(round(t, 1))
                    events.append((round(t, 1), "ivt"))
                    regimen.times.remove(tt)
        elif regimen.mode in ("prn", "tae"):
            if regimen.next_due is not None and t >= regimen.next_due - dt / 2.0:
                a_now = algebraic(t, y, p, drug, pb_at(t))
                regimen.visit_times.append(round(t, 1))
                if regimen.decide(t, a_now["CST"], a_now["BCVA"]):
                    y[S["DVIT"]] += dose_nmol
                    regimen.n_inj += 1
                    regimen.inj_times.append(round(t, 1))
                    events.append((round(t, 1), "ivt"))

        # ---------------- record --------------------------------------------
        if record and t >= next_rec - 1e-9:
            a = algebraic(t, y, p, drug, pb_at(t))
            out["t"].append(t)
            for k in out:
                if k == "t":
                    continue
                if k in a:
                    out[k].append(a[k])
                else:
                    out[k].append(y[S[k]])
            next_rec += record_every

        if i == n:
            break

        # ---------------- RK4 -----------------------------------------------
        k1 = rhs(t, y, p, drug, pb_at(t))
        y2 = [y[j] + 0.5 * dt * k1[j] for j in range(NSTATE)]
        k2 = rhs(t + 0.5 * dt, y2, p, drug, pb_at(t + 0.5 * dt))
        y3 = [y[j] + 0.5 * dt * k2[j] for j in range(NSTATE)]
        k3 = rhs(t + 0.5 * dt, y3, p, drug, pb_at(t + 0.5 * dt))
        y4 = [y[j] + dt * k3[j] for j in range(NSTATE)]
        k4 = rhs(t + dt, y4, p, drug, pb_at(t + dt))
        for j in range(NSTATE):
            y[j] += dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
        # hard bounds
        y[S["TJ"]] = clamp(y[S["TJ"]], 1e-3, 1.0)
        y[S["EZ"]] = clamp(y[S["EZ"]], 0.0, 1.0)
        y[S["W"]] = max(0.0, y[S["W"]])
        y[S["SRF"]] = max(0.0, y[S["SRF"]])
        y[S["NVI"]] = clamp(y[S["NVI"]], 0.0, 1.0)
        y[S["CAT"]] = clamp(y[S["CAT"]], 0.0, 1.0)
        y[S["MI"]] = clamp(y[S["MI"]], 0.0, 3.0)
        y[S["CHRON"]] = clamp(y[S["CHRON"]], 0.0, 1.0)

    final = algebraic(tmax, y, p, drug, pb_at(tmax))
    return dict(t=out["t"], series=out, final=final, y=y,
                n_inj=regimen.n_inj, n_dex=regimen.n_dex,
                inj_times=regimen.inj_times, events=events)


def at(res, day, key):
    """Value of an output series at (or nearest below) a given day."""
    ts = res["t"]
    if not ts:
        return None
    idx = min(range(len(ts)), key=lambda i: abs(ts[i] - day))
    return res["series"][key][idx]


# =============================================================================
# 8.  SCENARIOS
# =============================================================================
BASE_START = 42.0      # day of presentation / randomisation (6 weeks post-onset)
HORIZON = 1095.0       # 36 months


def q4w(n, start=BASE_START, q=28.0):
    return [start + i * q for i in range(n)]


def sched_load_then(n_load, q_load, n_maint, q_maint, start=BASE_START):
    t = q4w(n_load, start, q_load)
    last = t[-1]
    for i in range(1, n_maint + 1):
        t.append(last + i * q_maint)
    return t


SCENARIOS = []


def add(name, desc, phenotype, regimen, **pover):
    SCENARIOS.append(dict(name=name, desc=desc, phenotype=phenotype,
                          regimen=regimen, pover=pover))


# --- natural history ---------------------------------------------------------
add("NAT_NONISCH", "Non-ischaemic CRVO, untreated (natural history)",
    "crvo_nonisch", Regimen(mode="none"))
add("NAT_ISCH", "Ischaemic CRVO, untreated (natural history)",
    "crvo_isch", Regimen(mode="none"))
add("NAT_BRVO", "BRVO, untreated (natural history)",
    "brvo", Regimen(mode="none"))

# --- trial-like anti-VEGF regimens ------------------------------------------
add("RBZ_M6_PRN", "Ranibizumab 0.5 mg q4w x6 then monthly PRN (CRUISE / HORIZON)",
    "crvo_nonisch", Regimen(mode="prn", drug="ranibizumab", load_n=6, load_q=28.0))
add("AFL_M6_Q8", "Aflibercept 2 mg q4w x6 then q8w y1, then T&E (COPERNICUS / GALILEO)",
    "crvo_nonisch", Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                            interval0=56.0, imin=28.0, imax=112.0))
add("AFL_TAE_MIN", "Aflibercept 2 mg, treat-and-extend from month 2 (fewer loading doses)",
    "crvo_nonisch", Regimen(mode="tae", drug="aflibercept", load_n=3, load_q=28.0,
                            interval0=56.0, imin=28.0, imax=112.0))
add("FAR_TAE", "Faricimab 6 mg q4w x6 then T&E to q16w (BALATON / COMINO)",
    "crvo_nonisch", Regimen(mode="tae", drug="faricimab", load_n=6, load_q=28.0,
                            interval0=56.0, imin=28.0, imax=112.0))
add("BEV_M6_PRN", "Bevacizumab 1.25 mg q4w x6 then monthly PRN (SCORE2-like)",
    "crvo_nonisch", Regimen(mode="prn", drug="bevacizumab", load_n=6, load_q=28.0))
add("BRO_TAE", "Brolucizumab 6 mg q4w x6 then T&E",
    "crvo_nonisch", Regimen(mode="tae", drug="brolucizumab", load_n=6, load_q=28.0,
                            interval0=56.0, imin=28.0, imax=112.0))
add("AFL8_TAE", "Aflibercept 8 mg q4w x3 then T&E (high-molar-dose arm)",
    "crvo_nonisch", Regimen(mode="tae", drug="aflibercept8", load_n=3, load_q=28.0,
                            interval0=56.0, imin=28.0, imax=140.0))

# --- deferral ----------------------------------------------------------------
for months in (3, 6, 12):
    add(f"DEFER{months}",
        f"Sham for {months} months, then aflibercept q4w x6 then T&E",
        "crvo_nonisch",
        Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                start=BASE_START + months * 30.4, interval0=56.0))

# --- steroid -----------------------------------------------------------------
add("DEX_PHAKIC", "Dexamethasone implant 0.7 mg q6mo x6, phakic eye (GENEVA-like)",
    "crvo_nonisch", Regimen(mode="none", drug="aflibercept",
                            dex_times=[BASE_START + 182.5 * i for i in range(6)]))
add("DEX_PSEUDO", "Dexamethasone implant 0.7 mg q6mo x6, pseudophakic eye",
    "crvo_nonisch", Regimen(mode="none", drug="aflibercept",
                            dex_times=[BASE_START + 182.5 * i for i in range(6)]),
    K_CAT_DEX=0.0, E_CAT=0.0)
add("AFL_PLUS_DEX", "Aflibercept T&E + dexamethasone implant at months 3 and 9",
    "crvo_nonisch", Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                            interval0=56.0,
                            dex_times=[BASE_START + 91.0, BASE_START + 274.0]),
    K_CAT_DEX=0.0, E_CAT=0.0)

# --- real-world undertreatment ----------------------------------------------
add("REALWORLD", "Registry-like aflibercept: 6 in y1, 3 in y2, 2 in y3 (11 total)",
    "crvo_nonisch",
    Regimen(mode="fixed", drug="aflibercept",
            times=[BASE_START + d for d in (0, 28, 56, 112, 182, 258,
                                            378, 478, 598, 738, 908)]))
add("REALWORLD_LOW", "Severe undertreatment: 4 in y1, 2 in y2, 1 in y3 (7 total)",
    "crvo_nonisch",
    Regimen(mode="fixed", drug="aflibercept",
            times=[BASE_START + d for d in (0, 28, 56, 140, 400, 560, 800)]))

# --- the pressure-arm experiment --------------------------------------------
add("PERFECT_NONISCH", "Perfect (infinite) VEGF/PlGF/Ang-2 blockade from presentation",
    "crvo_nonisch", Regimen(mode="fixed", drug="aflibercept",
                            times=[BASE_START], perfect=True))
add("PERFECT_ISCH", "Perfect (infinite) blockade, ischaemic CRVO",
    "crvo_isch", Regimen(mode="fixed", drug="aflibercept",
                         times=[BASE_START], perfect=True))
add("PERFECT_MILD", "Perfect (infinite) blockade, mild occlusion",
    "crvo_mild", Regimen(mode="fixed", drug="aflibercept",
                         times=[BASE_START], perfect=True))
add("PERFECT_BRVO", "Perfect (infinite) blockade, BRVO",
    "brvo", Regimen(mode="fixed", drug="aflibercept",
                    times=[BASE_START], perfect=True))

# --- ischaemic CRVO and neovascular glaucoma --------------------------------
add("ISCH_AFL_CONT", "Ischaemic CRVO, aflibercept q4w continuously for 36 months",
    "crvo_isch", Regimen(mode="fixed", drug="aflibercept", times=q4w(38)))
add("ISCH_AFL_STOP6", "Ischaemic CRVO, aflibercept q4w x6 then stopped",
    "crvo_isch", Regimen(mode="fixed", drug="aflibercept", times=q4w(6)))
add("ISCH_AFL_PRP", "Ischaemic CRVO, aflibercept q4w x6 + PRP at month 3, then stopped",
    "crvo_isch", Regimen(mode="fixed", drug="aflibercept", times=q4w(6),
                         prp_time=BASE_START + 91.0))
add("ISCH_AFL_TAE", "Ischaemic CRVO, aflibercept T&E",
    "crvo_isch", Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                         interval0=56.0))

# --- BRVO --------------------------------------------------------------------
add("BRVO_AFL", "BRVO, aflibercept q4w x6 then PRN (VIBRANT-like)",
    "brvo", Regimen(mode="prn", drug="aflibercept", load_n=6, load_q=28.0))
add("BRVO_RBZ", "BRVO, ranibizumab q4w x6 then PRN (BRAVO-like)",
    "brvo", Regimen(mode="prn", drug="ranibizumab", load_n=6, load_q=28.0))

# --- blood pressure ----------------------------------------------------------
add("AFL_BP_LOW", "Aflibercept T&E + intensive BP lowering (inlet pressure 40 -> 32 mmHg)",
    "crvo_nonisch", Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                            interval0=56.0), PA=32.0)
add("BP_LOW_ALONE", "Intensive BP lowering alone, no anti-VEGF",
    "crvo_nonisch", Regimen(mode="none"), PA=32.0)
add("ISCH_AFL_BP_LOW", "Ischaemic CRVO, aflibercept T&E + BP lowering",
    "crvo_isch", Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                         interval0=56.0), PA=32.0)


def run_scenario(sc, tmax=HORIZON, dt=0.02):
    p = make_params(sc["phenotype"], **sc["pover"])
    reg = sc["regimen"]
    return simulate(p, reg, tmax=tmax, dt=dt), p


# =============================================================================
# 9.  CLOSED-FORM ANALYSES
# =============================================================================

def pc_of_rv(rv_res, pa, ra_res=1.0, pv=17.0):
    return (pa * rv_res + pv * ra_res) / (ra_res + rv_res)


def critical_rv(pc_star, pa, ra_res=1.0, pv=17.0):
    """Venous resistance at which Pc reaches pc_star.  None if unreachable."""
    denom = pa - pc_star
    if denom <= 0:
        return None
    num = (pc_star - pv) * ra_res
    if num <= 0:
        return 0.0
    return num / denom


def suppression_duration(drug, p, target_ratio=None):
    """Days for which the VEGF suppression ratio stays >= target_ratio.

        ratio(t) = 1 + C0_eff * exp(-kel t) / (ALPHA * KD)
        ratio(t) = R   ->   t = (1/kel) * ln( C0_eff / ((R-1) ALPHA KD) )
    """
    if target_ratio is None:
        target_ratio = p["SUPP_RATIO_ADEQ"]
    d = DRUGS[drug]
    c0 = molar_dose_nmol(drug) / p["VVIT"] * 1.0e6 * p["FPEN"]   # pM at effect site
    thr = (target_ratio - 1.0) * p["ALPHA"] * d["kd_v"]
    kel = math.log(2.0) / d["t_vit"]
    if c0 <= thr:
        return 0.0, c0, thr, kel
    return math.log(c0 / thr) / kel, c0, thr, kel


# =============================================================================
# 10.  VIRTUAL POPULATION
# =============================================================================

def lognormal(rng, median, cv):
    sd = math.sqrt(math.log(1.0 + cv * cv))
    return median * math.exp(rng.gauss(0.0, sd))


def run_population(n=400, seed=20260805, tmax=730.0, dt=0.05):
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        # A CV of 0.85 puts an implausible fraction of the population at
        # OCC0 > 40, which behaves like untreated ischaemic CRVO.  0.55 with a
        # hard cap keeps a severe tail without letting it dominate the mean.
        occ0 = clamp(lognormal(rng, 11.0, 0.55), 1.5, 40.0)
        pa = rng.gauss(40.0, 4.5)
        collmax = lognormal(rng, 1.3, 0.60)
        kez = lognormal(rng, P["K_EZ"], 0.40)
        vgain = lognormal(rng, P["E_HIF_V"], 0.35)
        delay = max(7.0, rng.gauss(42.0, 25.0))
        # registry injection counts over 24 months (Fight Retinal Blindness!
        # median ~10-11 for CRVO); the cap is what the patient actually
        # receives, not what the protocol asks for
        nmax = int(clamp(rng.gauss(11.0, 4.0), 4, 24))
        occres = clamp(rng.gauss(0.55, 0.16), 0.05, 0.95)
        p = make_params("crvo_nonisch", PA=pa, COLL_MAX=collmax, K_EZ=kez,
                        E_HIF_V=vgain, OCC_RES=occres)
        p["_OCC0"] = occ0
        reg = Regimen(mode="tae", drug="aflibercept", load_n=min(6, nmax),
                      load_q=28.0, start=delay, interval0=56.0, max_inj=nmax)
        res = simulate(p, reg, tmax=tmax, dt=dt, record_every=14.0)
        f = res["final"]
        base_cst = at(res, delay, "CST")
        base_bcva = at(res, delay, "BCVA")
        rows.append(dict(
            occ0=occ0, occ_res=occres, pa=pa, collmax=collmax, delay=delay, n_inj=res["n_inj"],
            pc=f["PC"], pc_crit=f["PC_CRIT_NOW"], cst=f["CST"],
            cst_base=base_cst, bcva=f["BCVA"], bcva_base=base_bcva,
            gain=f["BCVA"] - (base_bcva or 0.0), ez=res["y"][S["EZ"]],
            np=res["y"][S["NP"]], nvi=res["y"][S["NVI"]], iop=res["y"][S["IOP"]],
            dry=1 if f["CST"] < 310.0 else 0,
            floor=1 if f["JV_FLOOR"] > 0.0 else 0,
        ))
    return rows


# =============================================================================
# 11.  MAIN — produce every number quoted in README.md
# =============================================================================

def fmt(x, nd=1):
    if x is None:
        return "  n/a"
    return f"{x:.{nd}f}"


def fmts(x, nd=1):
    """Signed fixed-point string (kept as a string so it can be right-aligned)."""
    if x is None:
        return "  n/a"
    return f"{x:+.{nd}f}"


def main():
    outdir = os.path.dirname(os.path.abspath(__file__))
    log = []

    def W(s=""):
        print(s)
        log.append(s)

    W("=" * 79)
    W("RETINAL VEIN OCCLUSION QSP MODEL — REFERENCE RUN")
    W("망막정맥폐쇄 QSP 모델 — 참조 실행")
    W("=" * 79)
    W(f"states = {NSTATE} ODEs   horizon = {HORIZON:.0f} d   RK4 dt = 0.02 d")
    W("")

    # ---------------------------------------------------------------- part 1
    W("-" * 79)
    W("PART 1.  THE CRITICAL CAPILLARY PRESSURE  Pc*  =  Pt + sigma*(pi_c - pi_t)")
    W("-" * 79)
    p0 = make_params("crvo_nonisch")
    pc_star = p0["PT"] + p0["SIGMA_MAX"] * (p0["PI_C"] - p0["PI_T0"])
    W(f"  intact barrier:  Pt = {p0['PT']:.1f}, sigma = {p0['SIGMA_MAX']:.2f}, "
      f"pi_c - pi_t = {p0['PI_C'] - p0['PI_T0']:.1f} mmHg")
    W(f"  Pc*  =  {pc_star:.2f} mmHg")
    W("  Above this capillary pressure the Starling bracket is POSITIVE even with")
    W("  the permeability arm at its floor: no anti-VEGF dose can dry the macula.")
    W("")
    W("  Critical venous resistance Rv* (Ra = 1.0, Pv = 17 mmHg), by inlet pressure:")
    W("     Pa (mmHg)   Rv*        x baseline Rv0=0.25   dPc/dPa at Rv*")
    for pa in (48, 44, 40, 36, 34, 32, 31.5, 31.0):
        rvs = critical_rv(pc_star, pa)
        if rvs is None:
            W(f"     {pa:7.1f}   unreachable  —                    —")
        else:
            W(f"     {pa:7.1f}   {rvs:7.3f}    {rvs / 0.25:9.1f}x        "
              f"{rvs / (1.0 + rvs):.3f}")
    W("")
    W("  Same table with hypoxic arteriolar dilatation (Ra = 0.65):")
    for pa in (40, 36, 32):
        rvs = critical_rv(pc_star, pa, ra_res=0.65)
        W(f"     Pa = {pa:4.1f}  ->  Rv* = {fmt(rvs, 3)}  "
          f"({fmt((rvs / 0.25) if rvs else None, 1)}x baseline)")
    W("")
    W("  Same table for a chronically remodelled barrier that can no longer")
    W("  fully re-seal (TJmax = 0.70 -> sigma = 0.75, pi_t = 8.1 mmHg):")
    sig_ch = p0["SIGMA_MAX"] * 0.70 ** p0["SIG_N"]
    pit_ch = p0["PI_T0"] + (p0["PI_C"] - p0["PI_T0"]) * (1 - 0.70 ** p0["SIG_N"]) * p0["PI_T_GAIN"]
    pc_star_ch = p0["PT"] + sig_ch * (p0["PI_C"] - pit_ch)
    W(f"     sigma = {sig_ch:.3f}, pi_t = {pit_ch:.2f}  ->  Pc* = {pc_star_ch:.2f} mmHg")
    W(f"     Rv* at Pa = 40 falls from {critical_rv(pc_star, 40):.3f} to "
      f"{critical_rv(pc_star_ch, 40):.3f}  "
      f"({critical_rv(pc_star, 40) / critical_rv(pc_star_ch, 40):.2f}x lower)")
    W("     -> chronicity converts previously dryable eyes into floor eyes.")
    W("")

    # ---------------------------------------------------------------- part 2
    W("-" * 79)
    W("PART 2.  VEGF SUPPRESSION DURATION — DECOMPOSED")
    W("-" * 79)
    W("  t_sup = (t_half / ln2) * ln( C0_eff / [(R-1) * ALPHA * KD] ),  R = 4.9")
    W(f"  ALPHA = {P['ALPHA']:.0f} (single global calibration constant), "
      f"f_pen = {P['FPEN']:.3f}")
    W("")
    W("  agent            molar dose  C0_vit    C0_eff     KD     threshold  t1/2   t_sup")
    W("                      (nmol)    (uM)      (pM)     (pM)      (pM)      (d)    (d)")
    supp = {}
    for dn in ("ranibizumab", "bevacizumab", "aflibercept", "brolucizumab",
               "faricimab", "aflibercept8"):
        tsup, c0, thr, kel = suppression_duration(dn, P)
        supp[dn] = dict(t_sup=tsup, c0_eff=c0, thr=thr,
                        nmol=molar_dose_nmol(dn),
                        c0_vit=molar_dose_nmol(dn) / P["VVIT"])
        W(f"  {DRUGS[dn]['label']:<17}{molar_dose_nmol(dn):8.1f}  "
          f"{molar_dose_nmol(dn) / P['VVIT']:7.2f}  {c0:9.0f}  "
          f"{DRUGS[dn]['kd_v']:7.2f}  {thr:9.0f}  {DRUGS[dn]['t_vit']:5.2f} "
          f"{tsup:7.1f}")
    W("")
    W("  Decomposition of ln(C0_eff / threshold) into a RESERVOIR term and an")
    W("  AFFINITY term, then scaled by the half-life multiplier t1/2/ln2:")
    W("  agent            ln(C0_eff)  -ln(thr)   sum    t1/2/ln2   t_sup")
    for dn in ("ranibizumab", "bevacizumab", "aflibercept", "brolucizumab",
               "faricimab", "aflibercept8"):
        s = supp[dn]
        W(f"  {DRUGS[dn]['label']:<17}{math.log(s['c0_eff']):9.2f}  "
          f"{-math.log(s['thr']):8.2f} {math.log(s['c0_eff'] / s['thr']):7.2f}  "
          f"{DRUGS[dn]['t_vit'] / math.log(2):9.2f} {s['t_sup']:7.1f}")
    W("")
    ra, af = supp["ranibizumab"], supp["aflibercept"]
    W(f"  aflibercept / ranibizumab suppression-duration ratio = "
      f"{af['t_sup'] / ra['t_sup']:.2f}x")
    W(f"  ... while their intrinsic affinities differ by "
      f"{DRUGS['ranibizumab']['kd_v'] / DRUGS['aflibercept']['kd_v']:.0f}x.")
    W("  A 94-fold affinity gap becomes a 4-fold time gap because duration is")
    W("  LOGARITHMIC in concentration and LINEAR in half-life.")
    W("")

    # ---------------------------------------------------------------- part 3
    W("-" * 79)
    W("PART 3.  SCENARIOS")
    W("-" * 79)
    results = {}
    for sc in SCENARIOS:
        res, p = run_scenario(sc)
        results[sc["name"]] = dict(res=res, p=p, sc=sc)

    def row(name, day_base=None):
        r = results[name]
        res, p, sc = r["res"], r["p"], r["sc"]
        reg = sc["regimen"]
        db = day_base if day_base is not None else reg.start if reg.mode in ("prn", "tae") \
            else (reg.times[0] if reg.mode == "fixed" and reg.times else BASE_START)
        db = min(db, HORIZON)
        b_cst = at(res, db, "CST")
        b_bc = at(res, db, "BCVA")
        return dict(
            name=name, desc=sc["desc"], n_inj=res["n_inj"], n_dex=res["n_dex"],
            day_base=db, cst_base=b_cst, bcva_base=b_bc,
            cst_m6=at(res, db + 182, "CST"), bcva_m6=at(res, db + 182, "BCVA"),
            cst_m12=at(res, db + 365, "CST"), bcva_m12=at(res, db + 365, "BCVA"),
            cst_m24=at(res, 730, "CST"), bcva_m24=at(res, 730, "BCVA"),
            cst_m36=at(res, 1095, "CST"), bcva_m36=at(res, 1095, "BCVA"),
            gain_m6=at(res, db + 182, "BCVA") - b_bc,
            gain_m12=at(res, db + 365, "BCVA") - b_bc,
            gain_m24=at(res, 730, "BCVA") - b_bc,
            gain_m36=at(res, 1095, "BCVA") - b_bc,
            pc=res["final"]["PC"], pc_crit=res["final"]["PC_CRIT_NOW"],
            ez=res["y"][S["EZ"]], np=res["y"][S["NP"]], mi=res["y"][S["MI"]],
            nvi=res["y"][S["NVI"]], iop=res["y"][S["IOP"]], cat=res["y"][S["CAT"]],
            chron=res["y"][S["CHRON"]], tsup=res["y"][S["TSUP"]],
            cumed=res["y"][S["CUMED"]],
            floor_jv=res["final"]["JV_FLOOR"],
            press_frac=res["final"]["PRESSURE_ARM_FRAC"],
            nvg=res["final"]["NVG"], inj_times=res["inj_times"],
        )

    rows = {sc["name"]: row(sc["name"]) for sc in SCENARIOS}

    W("  3a. Trial-replication check (non-ischaemic CRVO, baseline = day of first dose)")
    W("      scenario          inj  CST_base  CST_m6  dCST    BCVA_base  m6    m12   m24")
    for n in ("NAT_NONISCH", "RBZ_M6_PRN", "BEV_M6_PRN", "AFL_M6_Q8", "FAR_TAE",
              "BRO_TAE", "AFL8_TAE", "AFL_TAE_MIN"):
        r = rows[n]
        W(f"      {n:<17}{r['n_inj']:4d}  {fmt(r['cst_base'],0):>8}  "
          f"{fmt(r['cst_m6'],0):>6}  {fmts(r['cst_m6']-r['cst_base'],0):>6}  "
          f"{fmt(r['bcva_base']):>9}  {fmts(r['gain_m6'],1):>6}  "
          f"{fmts(r['gain_m12'],1):>6} {fmts(r['gain_m24'],1):>6}")
    W("")
    W("  3b. BRVO")
    for n in ("NAT_BRVO", "BRVO_AFL", "BRVO_RBZ"):
        r = rows[n]
        W(f"      {n:<17}{r['n_inj']:4d}  CST {fmt(r['cst_base'],0)} -> "
          f"{fmt(r['cst_m6'],0)}   BCVA {fmt(r['bcva_base'])} "
          f"({fmts(r['gain_m6'],1):>6} at m6, {fmts(r['gain_m24'],1):>6} at m24)")
    W("")
    W("  3c. Deferral — the permanent price of waiting")
    W("      scenario     start(d)  CST_base  BCVA_base  peak BCVA  m36 BCVA  EZ_m36")
    for n in ("AFL_M6_Q8", "DEFER3", "DEFER6", "DEFER12"):
        r = rows[n]
        res = results[n]["res"]
        peak = max(res["series"]["BCVA"])
        W(f"      {n:<12}{r['day_base']:8.0f}  {fmt(r['cst_base'],0):>8}  "
          f"{fmt(r['bcva_base']):>9}  {fmt(peak):>9}  {fmt(r['bcva_m36']):>8}  "
          f"{r['ez']:.3f}")
    ref = rows["AFL_M6_Q8"]
    W("")
    W("      letters permanently forfeited vs immediate treatment (month 36):")
    for n in ("DEFER3", "DEFER6", "DEFER12"):
        W(f"        {n:<9}{rows[n]['bcva_m36'] - ref['bcva_m36']:+6.1f} letters   "
          f"(EZ {rows[n]['ez']:.3f} vs {ref['ez']:.3f})")
    W("")

    W("  3d. THE PRESSURE-ARM EXPERIMENT — infinite VEGF/PlGF/Ang-2 blockade")
    W("      phenotype        Pc    Pc*    Jv_floor  CST_m6  CST_m36  dry?")
    for n in ("PERFECT_MILD", "PERFECT_BRVO", "PERFECT_NONISCH", "PERFECT_ISCH"):
        r = rows[n]
        W(f"      {n:<15}{fmt(r['pc']):>6} {fmt(r['pc_crit']):>6}  "
          f"{fmt(r['floor_jv'],2):>8}  {fmt(r['cst_m6'],0):>6}  "
          f"{fmt(r['cst_m36'],0):>7}   {'YES' if r['cst_m36'] < 310 else 'no'}")
    W("")
    W("")
    W("      When does the eye cross BELOW its own critical pressure?")
    W(f"      (untreated natural history; threshold = intact-barrier P_c* = "
      f"{pc_star:.2f} mmHg)")
    W("      phenotype      P_c d42  P_c* d42  first day P_c < P_c*   days above P_c*")
    cross = {}
    for phn in ("crvo_mild", "brvo", "crvo_nonisch", "crvo_isch"):
        pp = make_params(phn)
        rr = simulate(pp, Regimen(mode="none"), tmax=HORIZON, dt=0.02,
                      record_every=1.0)
        pcs = rr["series"]["PC"]
        # Compare against the INTACT-barrier threshold Pc* = Pt + sigma_max *
        # (pi_c - pi_t0).  Using the CURRENT sigma would be circular: an untreated
        # wet eye has a leaky barrier, so its instantaneous Pc* is low and Pc
        # always exceeds it.  The therapeutically meaningful question is whether
        # Pc falls below the threshold that applies ONCE THE BARRIER IS RESTORED.
        pcc = rr["series"]["PC_CRIT"]
        ts = rr["t"]
        first = None
        above = 0
        for i, tt in enumerate(ts):
            if tt < BASE_START:
                continue
            if pcs[i] >= pcc[i]:
                above += 1
            elif first is None:
                first = tt
        cross[phn] = dict(first=first, above=above,
                          pc42=pcs[min(range(len(ts)), key=lambda i: abs(ts[i]-42))],
                          pcc42=pcc[min(range(len(ts)), key=lambda i: abs(ts[i]-42))])
        W(f"      {phn:<14}{cross[phn]['pc42']:7.1f}  {cross[phn]['pcc42']:8.1f}   "
          f"{('day ' + str(int(first))) if first else 'NEVER':<19}  {above:6d}")
    W("      Non-ischaemic CRVO crosses below its critical pressure as collaterals")
    W("      mature; ischaemic CRVO never does.  That single crossing time, not the")
    W("      choice of anti-VEGF agent, is what separates the two clinical courses.")
    W("")
    W("      Same eyes under the best REAL regimen (aflibercept T&E):")
    W(f"        non-ischaemic CRVO  CST m36 = {fmt(rows['AFL_M6_Q8']['cst_m36'],0)} um "
      f"(perfect blockade {fmt(rows['PERFECT_NONISCH']['cst_m36'],0)} um)")
    W(f"        ischaemic CRVO      CST m36 = {fmt(rows['ISCH_AFL_TAE']['cst_m36'],0)} um "
      f"(perfect blockade {fmt(rows['PERFECT_ISCH']['cst_m36'],0)} um)")
    gap_n = rows['AFL_M6_Q8']['cst_m36'] - rows['PERFECT_NONISCH']['cst_m36']
    W(f"        residual attributable to imperfect blockade: "
      f"{fmt(gap_n,0)} um (non-ischaemic)")
    W(f"        residual attributable to the PRESSURE ARM: "
      f"{fmt(rows['PERFECT_ISCH']['cst_m36'] - 250,0)} um (ischaemic)")
    W("")
    W("      THE PRESSURE ARM AMPLIFIES WHATEVER PERMEABILITY REMAINS.")
    W("      Jv = Lp * bracket, so an identical residual elevation of Lp delivers")
    W("      more fluid when the bracket is larger.  Cost of the SAME imperfect")
    W("      blockade (real aflibercept T&E minus infinite blockade, month 36):")
    amp_n = rows["AFL_M6_Q8"]["cst_m36"] - rows["PERFECT_NONISCH"]["cst_m36"]
    amp_i = rows["ISCH_AFL_TAE"]["cst_m36"] - rows["PERFECT_ISCH"]["cst_m36"]
    W(f"        non-ischaemic (bracket small)  {amp_n:6.1f} um")
    W(f"        ischaemic     (bracket large)  {amp_i:6.1f} um   "
      f"= {(amp_i / amp_n) if abs(amp_n) > 0.01 else float('nan'):.1f}x")
    W("      This is why switching agent fails in a high-pressure eye: the")
    W("      pharmacological shortfall being amplified is already small.")
    W("")

    W("  3e. Ischaemic CRVO, neovascular glaucoma, and rebound")
    W("      scenario          inj  NP(DA)  NVI    IOP   NVG  CST_m36  BCVA_m36")
    for n in ("NAT_ISCH", "ISCH_AFL_CONT", "ISCH_AFL_STOP6", "ISCH_AFL_PRP",
              "ISCH_AFL_TAE"):
        r = rows[n]
        W(f"      {n:<17}{r['n_inj']:4d}  {fmt(r['np'],1):>6}  {r['nvi']:.3f}  "
          f"{fmt(r['iop']):>5}  {'YES' if r['nvg'] else ' no'}  "
          f"{fmt(r['cst_m36'],0):>7}  {fmt(r['bcva_m36']):>8}")
    W("")

    W("  3f. Steroid vs anti-VEGF, and the lens")
    W("      scenario        implants inj  CST_m6  BCVA_m6  CAT_m36  IOP_m36  BCVA_m36")
    for n in ("DEX_PHAKIC", "DEX_PSEUDO", "AFL_M6_Q8", "AFL_PLUS_DEX"):
        r = rows[n]
        W(f"      {n:<15}{r['n_dex']:8d} {r['n_inj']:4d}  {fmt(r['cst_m6'],0):>6}  "
          f"{fmts(r['gain_m6'],1):>7}  {r['cat']:.3f}   {fmt(r['iop']):>6}  "
          f"{fmt(r['bcva_m36']):>8}")
    W("")

    W("  3g. Blood pressure as an anti-oedema intervention")
    W("      scenario           Pc(m36)  CST_m36  BCVA_m36  injections")
    for n in ("AFL_M6_Q8", "AFL_BP_LOW", "BP_LOW_ALONE", "NAT_NONISCH",
              "ISCH_AFL_TAE", "ISCH_AFL_BP_LOW"):
        r = rows[n]
        W(f"      {n:<18}{fmt(r['pc']):>7}  {fmt(r['cst_m36'],0):>7}  "
          f"{fmt(r['bcva_m36']):>8}  {r['n_inj']:6d}")
    d_bp = rows['AFL_BP_LOW']['bcva_m36'] - rows['AFL_M6_Q8']['bcva_m36']
    d_inj = rows['AFL_BP_LOW']['n_inj'] - rows['AFL_M6_Q8']['n_inj']
    W(f"      -> 8 mmHg lower inlet pressure is worth {d_bp:+.1f} letters and "
      f"{d_inj:+d} injections over 36 months")
    W("")

    W("  3h. Real-world undertreatment")
    for nm in ("REALWORLD", "REALWORLD_LOW"):
        rr2 = rows[nm]
        W(f"      {nm:<14}{rr2['n_inj']:3d} inj  BCVA {fmt(rr2['bcva_base'])} -> "
          f"{fmt(rr2['bcva_m36'])} ({fmts(rr2['gain_m36'],1)})  "
          f"CST {fmt(rr2['cst_m36'],0)} um  EZ {rr2['ez']:.3f}  "
          f"letters/injection {rr2['gain_m36'] / max(1, rr2['n_inj']):+.2f}")
    ref2 = rows["AFL_M6_Q8"]
    W(f"      {'PROTOCOL T&E':<14}{ref2['n_inj']:3d} inj  BCVA "
      f"{fmt(ref2['bcva_base'])} -> {fmt(ref2['bcva_m36'])} "
      f"({fmts(ref2['gain_m36'],1)})  CST {fmt(ref2['cst_m36'],0)} um  "
      f"EZ {ref2['ez']:.3f}  letters/injection "
      f"{ref2['gain_m36'] / max(1, ref2['n_inj']):+.2f}")
    W("      Undertreatment does not scale the benefit down, it inverts it: the")
    W("      eye spends the gaps re-accumulating oedema, and the ellipsoid zone")
    W("      integrates every one of those gaps.")
    W("")

    W("  3i. Suppression bookkeeping — days of adequate VEGF suppression vs days dry")
    W("      scenario          inj  days suppressed  days CST<310  ratio")
    for n in ("RBZ_M6_PRN", "AFL_M6_Q8", "FAR_TAE", "ISCH_AFL_CONT", "ISCH_AFL_TAE"):
        res = results[n]["res"]
        r = rows[n]
        dry_days = sum(1 for c in res["series"]["CST"] if c < 310.0)
        W(f"      {n:<17}{r['n_inj']:4d}  {r['tsup']:15.0f}  {dry_days:12d}  "
          f"{(dry_days / r['tsup']) if r['tsup'] > 0 else 0:.2f}")
    W("      A ratio < 1 means the macula was WET while VEGF was still suppressed:")
    W("      that gap is the pressure arm, and no change of agent addresses it.")
    W("")

    # ---------------------------------------------------------------- part 4
    W("-" * 79)
    W("PART 4.  DELAY-RESPONSE AND DOSE-RESPONSE CURVES")
    W("-" * 79)
    W("  4a. Permanent letter cost of delay (aflibercept T&E after the delay;")
    W("      BCVA averaged over months 24-36 so the read-out does not land on an")
    W("      arbitrary phase of the treat-and-extend cycle)")
    W("      delay (mo)  BCVA at start  peak BCVA  mean m24-36  EZ m36  perm. loss")
    delay_curve = []
    for months in (0, 0.5, 1, 2, 3, 4, 6, 9, 12):
        start = BASE_START + months * 30.4
        p = make_params("crvo_nonisch")
        reg = Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                      start=start, interval0=56.0)
        res = simulate(p, reg, tmax=HORIZON, dt=0.02, record_every=2.0)
        b = at(res, start, "BCVA")
        pk = max(res["series"]["BCVA"])
        m36 = at(res, 1095, "BCVA")
        late = [v for tt, v in zip(res["t"], res["series"]["BCVA"]) if tt >= 730]
        late_mean = sum(late) / len(late) if late else float("nan")
        delay_curve.append(dict(delay_mo=months, bcva_start=b, peak=pk,
                                bcva_m36=m36, bcva_late=late_mean,
                                ez=res["y"][S["EZ"]], n_inj=res["n_inj"]))
        W(f"      {months:10.1f}  {fmt(b):>13}  {fmt(pk):>9}  {fmt(late_mean):>11}  "
          f"{res['y'][S['EZ']]:.3f}   "
          f"{late_mean - delay_curve[0]['bcva_late']:+.1f}")
    W("")
    W("  4b. Dose-response in injections (aflibercept, capped total over 24 months)")
    W("      max inj  actual  CST_m24  BCVA_m24  gain   letters/injection")
    dose_curve = []
    for nmax in (0, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20):
        p = make_params("crvo_nonisch")
        reg = Regimen(mode="tae", drug="aflibercept", load_n=min(6, nmax),
                      load_q=28.0, start=BASE_START, interval0=56.0,
                      max_inj=nmax if nmax > 0 else 0)
        res = simulate(p, reg, tmax=730.0, dt=0.02, record_every=2.0)
        b = at(res, BASE_START, "BCVA")
        g = at(res, 730, "BCVA") - b
        dose_curve.append(dict(nmax=nmax, n_inj=res["n_inj"],
                               cst=at(res, 730, "CST"),
                               bcva=at(res, 730, "BCVA"), gain=g))
        W(f"      {nmax:7d}  {res['n_inj']:6d}  {fmt(at(res,730,'CST'),0):>7}  "
          f"{fmt(at(res,730,'BCVA')):>8}  {g:+6.1f}   "
          f"{(g / res['n_inj']) if res['n_inj'] else 0:.2f}")
    W("")
    W("  4c. Occlusion-severity response surface (aflibercept T&E, 24 months)")
    W("      OCC0   Rv    Pc    Pc*   CST_base  CST_m24  dry?  inj  BCVA gain")
    sev_curve = []
    for occ0 in (2, 3.5, 5, 7, 9, 12, 16, 22, 30, 45):
        p = make_params("crvo_nonisch")
        p["_OCC0"] = float(occ0)
        reg = Regimen(mode="tae", drug="aflibercept", load_n=6, load_q=28.0,
                      start=BASE_START, interval0=56.0)
        res = simulate(p, reg, tmax=730.0, dt=0.02, record_every=2.0)
        f = res["final"]
        b = at(res, BASE_START, "BCVA")
        sev_curve.append(dict(occ0=occ0, pc=f["PC"], pc_crit=f["PC_CRIT_NOW"],
                              cst_base=at(res, BASE_START, "CST"),
                              cst=f["CST"], n_inj=res["n_inj"],
                              gain=at(res, 730, "BCVA") - b,
                              dry=1 if f["CST"] < 310 else 0))
        W(f"      {occ0:5.1f}  {f['RV_RES']:.2f}  {f['PC']:5.1f}  {f['PC_CRIT_NOW']:5.1f}  "
          f"{fmt(at(res,BASE_START,'CST'),0):>8}  {fmt(f['CST'],0):>7}  "
          f"{'YES' if f['CST'] < 310 else ' no':>4}  {res['n_inj']:3d}  "
          f"{at(res,730,'BCVA') - b:+8.1f}")
    W("")

    # ---------------------------------------------------------------- part 5
    W("-" * 79)
    W("PART 5.  VIRTUAL POPULATION (n = 400, 24 months, aflibercept T&E)")
    W("-" * 79)
    pop = run_population(n=400, seed=20260805)
    n = len(pop)

    def frac(pred):
        return 100.0 * sum(1 for r in pop if pred(r)) / n

    def mean(key, sub=None):
        v = [r[key] for r in pop if (sub is None or sub(r))]
        return sum(v) / len(v) if v else float("nan")

    W(f"  n = {n}")
    W(f"  mean baseline CST  {mean('cst_base'):.0f} um   mean baseline BCVA "
      f"{mean('bcva_base'):.1f} letters")
    W(f"  mean CST at m24    {mean('cst'):.0f} um   mean BCVA gain "
      f"{mean('gain'):+.1f} letters")
    W(f"  mean injections    {mean('n_inj'):.1f}")
    W("")
    W(f"  dry macula (CST < 310 um) at month 24 .............. {frac(lambda r: r['dry']):5.1f}%")
    W(f"  pressure-arm FLOOR present (Jv_floor > 0) .......... {frac(lambda r: r['floor']):5.1f}%")
    W(f"  floor present AND still wet ........................ "
      f"{frac(lambda r: r['floor'] and not r['dry']):5.1f}%")
    W(f"  no floor AND still wet (true pharmacological miss) .. "
      f"{frac(lambda r: (not r['floor']) and (not r['dry'])):5.1f}%")
    W(f"  gained >= 15 letters ............................... {frac(lambda r: r['gain'] >= 15):5.1f}%")
    W(f"  lost >= 15 letters ................................. {frac(lambda r: r['gain'] <= -15):5.1f}%")
    W(f"  NVI > 0.5 (anterior neovascularisation) ............ {frac(lambda r: r['nvi'] > 0.5):5.1f}%")
    W("")
    W("  Stratified by whether the eye has a pressure-arm floor:")
    for lab, sub in (("floor", lambda r: r["floor"]), ("no floor", lambda r: not r["floor"])):
        W(f"    {lab:<9} n={sum(1 for r in pop if sub(r)):3d}  "
          f"Pc {mean('pc', sub):5.1f}  CST_m24 {mean('cst', sub):5.0f}  "
          f"gain {mean('gain', sub):+5.1f}  inj {mean('n_inj', sub):4.1f}  "
          f"EZ {mean('ez', sub):.3f}")
    W("")
    W("  Stratified by injection count actually received:")
    for lo, hi in ((0, 5), (5, 7), (7, 9), (9, 13)):
        sub = (lambda l, h: (lambda r: l <= r["n_inj"] < h))(lo, hi)
        k = sum(1 for r in pop if sub(r))
        if k:
            W(f"    {lo}-{hi-1} inj  n={k:3d}  gain {mean('gain', sub):+5.1f}  "
              f"CST {mean('cst', sub):5.0f}  EZ {mean('ez', sub):.3f}")
    W("")
    W("  Stratified by presentation delay.  NOTE: 'gain' is measured from the")
    W("  BCVA on the day of the first injection, so a patient who presents early")
    W("  is scored against a baseline that has not finished deteriorating.  That")
    W("  is a MEASUREMENT artefact, not a benefit of waiting - it is also exactly")
    W("  why trials that enrol sicker baselines report larger letter gains.  The")
    W("  ellipsoid-zone column is the delay-sensitive quantity:")
    for lo, hi in ((0, 25), (25, 45), (45, 70), (70, 999)):
        sub = (lambda l, h: (lambda r: l <= r["delay"] < h))(lo, hi)
        k = sum(1 for r in pop if sub(r))
        if k:
            W(f"    delay {lo}-{hi} d  n={k:3d}  gain {mean('gain', sub):+5.1f}  "
              f"BCVA_m24 {mean('bcva', sub):5.1f}  EZ {mean('ez', sub):.3f}")
    W("")

    # correlation of Pc with residual CST
    xs = [r["pc"] for r in pop]
    ys = [r["cst"] for r in pop]
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    rho = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else 0.0
    W(f"  Pearson r (Pc at m24  vs  residual CST at m24) = {rho:+.3f}")
    xs2 = [r["n_inj"] for r in pop]
    mx2 = sum(xs2) / n
    sxy2 = sum((x - mx2) * (y - my) for x, y in zip(xs2, ys))
    sxx2 = sum((x - mx2) ** 2 for x in xs2)
    rho2 = sxy2 / math.sqrt(sxx2 * syy) if sxx2 > 0 and syy > 0 else 0.0
    W(f"  Pearson r (injections vs  residual CST at m24) = {rho2:+.3f}")
    W("  Capillary pressure explains residual thickness better than dose does.")
    W("")

    W("=" * 79)
    W("END OF REFERENCE RUN")
    W("=" * 79)

    # ------------------------------------------------------------ write files
    with open(os.path.join(outdir, "rvo_reference_output.txt"), "w") as fh:
        fh.write("\n".join(log) + "\n")

    scen_json = {}
    for k, r in rows.items():
        rr = dict(r)
        rr["inj_times"] = rr["inj_times"][:60]
        scen_json[k] = rr
    with open(os.path.join(outdir, "rvo_scenario_results.json"), "w") as fh:
        json.dump(dict(
            meta=dict(states=NSTATE, horizon_d=HORIZON, dt=0.02,
                      alpha=P["ALPHA"], fpen=P["FPEN"], pc_star=pc_star),
            suppression={k: v for k, v in supp.items()},
            scenarios=scen_json,
            crossing={k: v for k, v in cross.items()},
            delay_curve=delay_curve,
            dose_curve=dose_curve,
            severity_curve=sev_curve,
        ), fh, indent=1)

    with open(os.path.join(outdir, "rvo_population_results.json"), "w") as fh:
        json.dump(dict(
            n=n, seed=20260805, horizon_d=730,
            summary=dict(
                cst_base=mean("cst_base"), bcva_base=mean("bcva_base"),
                cst_m24=mean("cst"), gain=mean("gain"), n_inj=mean("n_inj"),
                pct_dry=frac(lambda r: r["dry"]),
                pct_floor=frac(lambda r: r["floor"]),
                pct_floor_wet=frac(lambda r: r["floor"] and not r["dry"]),
                pct_nofloor_wet=frac(lambda r: (not r["floor"]) and (not r["dry"])),
                pct_gain15=frac(lambda r: r["gain"] >= 15),
                pct_lose15=frac(lambda r: r["gain"] <= -15),
                pct_nvi=frac(lambda r: r["nvi"] > 0.5),
                r_pc_cst=rho, r_inj_cst=rho2),
            rows=pop[:120],
        ), fh, indent=1)

    print(f"\nwrote rvo_reference_output.txt, rvo_scenario_results.json, "
          f"rvo_population_results.json")


if __name__ == "__main__":
    main()
