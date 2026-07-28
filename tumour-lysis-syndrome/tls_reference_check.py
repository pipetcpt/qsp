#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tls_reference_check.py
======================
Independent numpy/scipy transcription of the SAME 47-state system that lives in
`tls_mrgsolve_model.R`, written so that every number quoted in
`tumor-lysis-syndrome/README.md` can be re-derived without R or mrgsolve.

WHY THIS FILE EXISTS
--------------------
The build environment for this repository has no R toolchain.  Rather than
publish an ODE model whose behaviour has never been integrated, the system is
transcribed a second time here (scipy LSODA) and every headline claim in the
README is produced by one of the A0-A18 functions below.  If the R model and
this file ever disagree, that is a bug in one of them -- they are meant to be
the same system with the same parameter block.

THE STRUCTURAL CLAIM THIS MODEL EXISTS TO TEST
----------------------------------------------
Tumour lysis syndrome is normally taught as four laboratory abnormalities
(K+ up, PO4 up, urate up, Ca down) plus acute kidney injury.  This model
instead poses TLS as a RACE between two rates that are dimensionally identical
and numerically comparable:

      release flux         J_rel = Q_i * k_lys * N_lysing        [mmol/h]
      clearance capacity   C_i(GFR, urine pH, urine flow)        [mmol/h]

and then CLOSES THAT RACE INTO A LOOP, because the solutes that lose the race
precipitate, and the precipitate lowers the very GFR that sets the clearance
capacity:

      J_rel > C_i  ->  concentration up  ->  supersaturation up
                   ->  precipitation  ->  GFR down  ->  C_i down  -> ...

Everything else is a consequence of that loop, and every therapy is classified
by WHICH TERM of it the therapy touches:

  POOL           removes solute that already exists    rasburicase, dialysis
  FLUX           blocks new production                 allopurinol, febuxostat
  FLUX-SHAPING   lowers release RATE, same total       venetoclax ramp, prephase
  DILUTION       lowers tubular concentration only     hydration, furosemide
  SPECIATION     moves the solubility curve            bicarbonate -> urine pH
  REDISTRIBUTION moves solute between compartments     insulin/glucose, beta2
  SEQUESTRATION  removes solute outside the body       sevelamer, SZC

The predictions worth having all come from the fact that these classes are NOT
interchangeable: a drug in one class cannot substitute for a drug in another
no matter how much of it is given.

Run:
    python3 tls_reference_check.py               # every analysis
    python3 tls_reference_check.py --quick       # baseline + core
    python3 tls_reference_check.py --only A3,A9  # selected analyses
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# 0. STATE VECTOR  (47 ODEs)
# ---------------------------------------------------------------------------
STATES = [
    # --- drug PK (19) ---
    "ALLO_G",    # 0  allopurinol, gut                  (mg)
    "ALLO_C",    # 1  allopurinol, central              (mg)
    "OXY_C",     # 2  oxypurinol, central               (mg) renally cleared
    "FEBU_G",    # 3  febuxostat, gut                   (mg)
    "FEBU_C",    # 4  febuxostat, central               (mg) hepatic
    "RASB_C",    # 5  rasburicase, central              (mg)
    "RASB_P",    # 6  rasburicase, peripheral           (mg)
    "ADA",       # 7  anti-rasburicase antibody         (rel.)
    "VEN_G",     # 8  venetoclax, gut                   (mg)
    "VEN_C",     # 9  venetoclax, central               (mg)
    "CHEMO_C",   # 10 generic cytotoxic exposure        (rel.)
    "PRED_C",    # 11 prednisolone                      (mg/L)
    "INS_C",     # 12 insulin effect compartment        (mU/L)
    "GLU",       # 13 plasma glucose                    (mmol/L)
    "SALB_C",    # 14 salbutamol                        (ng/mL)
    "FURO_C",    # 15 furosemide                        (mg/L)
    "SEVE_G",    # 16 phosphate binder in gut           (mg)
    "SZC_G",     # 17 potassium binder in gut           (g)
    "HCO3_L",    # 18 bicarbonate load                  (mmol)
    # --- tumour / release (3) ---
    "N_VIA",     # 19 viable tumour                     (1e12 cells)
    "N_LYS",     # 20 committed-to-lysis pool           (1e12 cells)
    "LDH",       # 21 serum LDH                         (kU in VEC)
    # --- solute pools (11) ---
    "K_EC",      # 22 extracellular potassium           (mmol)
    "K_SH",      # 23 K parked intracellularly          (mmol)
    "PO4_EC",    # 24 extracellular phosphate           (mmol)
    "CA_EC",     # 25 extracellular total calcium       (mmol)
    "MG_EC",     # 26 extracellular magnesium           (mmol)
    "HYPOX",     # 27 hypoxanthine                      (mmol)
    "XAN",       # 28 xanthine                          (mmol)
    "URATE",     # 29 urate                             (mmol)
    "ALLANT",    # 30 allantoin                         (mmol)
    "CREAT",     # 31 creatinine                        (mmol)
    "UREA",      # 32 urea                              (mmol)
    # --- kidney (7) ---
    "NEPH",      # 33 viable nephron fraction           (0-1)
    "TUBINJ",    # 34 tubular injury                    (0-1)
    "XT_UA",     # 35 renal urate crystal               (mmol)
    "XT_XAN",    # 36 renal xanthine crystal            (mmol)
    "XT_CAP",    # 37 renal calcium-phosphate           (mmol)
    "PHU",       # 38 urine pH
    "VEC",       # 39 extracellular fluid volume        (L)
    # --- safety / endpoints (7) ---
    "H2O2",      # 40 cumulative peroxide               (mmol)
    "METHB",     # 41 methaemoglobin fraction
    "HB",        # 42 haemoglobin                       (g/dL)
    "PTH",       # 43 PTH                               (x normal)
    "CH_ARR",    # 44 cumulative arrhythmia hazard
    "CH_SZ",     # 45 cumulative seizure/tetany hazard
    "CH_RRT",    # 46 cumulative renal-replacement hazard
    "CUM_REL",   # 47 cumulative tumour lysed  (1e12 cells)
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)

# unit conversions (mmol/L -> conventional)
UA_MGDL, CR_MGDL, PO4_MGDL, CA_MGDL = 16.81, 11.31, 3.10, 4.01

# ---------------------------------------------------------------------------
# 1. PARAMETERS
# ---------------------------------------------------------------------------
P = dict(
    # ================= patient =================
    BW=70.0, BSA=1.73,
    VEC0=16.8,        # L extracellular fluid (0.24 * BW)
    VTBW=42.0,        # L total body water (creatinine / urea Vd)
    GFR0=7.2,         # L/h = 120 mL/min/1.73 m2
    HB0=13.5,

    # ================= fluid / urine =================
    FLUID=0.113,      # L/h total free-water input (= QU_BASE + INSENS)
    QU_BASE=0.083,    # L/h baseline urine (2 L/day)
    INSENS=0.030,     # L/h insensible loss
    KVOL=0.100,       # L/h extra urine per L of ECF excess
    FE_H2O_MAX=0.10,  # max fractional water excretion (intact tubule)
    OBSTR_FLOW=0.90,  # fraction of that lost at complete obstruction
    QU_FLOOR=0.004,
    FURO_QU=0.15, EC50_FURO=1.0, KEL_FURO=0.35, FURO_K=1.9,

    # ================= tumour =================
    N0=0.05,          # 1e12 cells
    TD=30.0,          # h doubling time
    KLYS=0.1155,      # /h transit out of the committed pool (t1/2 6 h)
    KKILL_SP=0.0,
    # intracellular content released per 1e12 cells
    Q_K=38.5,         # mmol  (110 mmol/L cell volume, 0.35 pL/cell)
    Q_PO4=60.0,       # mmol  (~48 from nucleic-acid phosphate + soluble pools)
    Q_PUR=24.0,       # mmol  (16 pg nucleic acid/cell, 50% purine, MW 330)
    Q_MG=5.0,
    Q_LDH=56.0,       # kU
    Q_CR=0.9,         # mmol creatinine equivalent
    LDH_PROD=0.0390,  # kU/h basal LDH turnover
    KEL_LDH=0.0116,   # /h  (t1/2 60 h)

    # ================= purine cascade =================
    P_END=0.186,      # mmol/h endogenous purine input
    VXO=2.0,          # mmol/L/h xanthine oxidase Vmax (hepatic, abundant)
    KM_XO1=0.05, KM_XO2=0.05,
    VSALV=0.075, KM_SALV=0.05,        # HGPRT salvage
    FE_HX=0.60, FE_XAN=0.55, FE_UA=0.080, FE_ALLANT=0.95,
    # urate reabsorption (URAT1/GLUT9) is saturable, so the fractional
    # excretion RISES with plasma urate -- the kidney's only defence
    FEUA_GAIN=3.5, FEUA_TH=0.40, FEUA_C50=0.55,

    # ================= xanthine oxidase inhibitors =================
    KA_ALLO=1.5, V_ALLO=90.0, CL_ALLO=45.0, F_ALLO=0.90, FM_OXY=0.75,
    MW_ALLO=136.1, MW_OXY=152.1,
    V_OXY=40.0, CL_OXY=1.20,          # L/h at GFR0; scales with GFR
    IC50_ALLO=6.0, IC50_OXY=1.10,
    KA_FEBU=1.2, V_FEBU=50.0, CL_FEBU=5.0, F_FEBU=0.85, IC50_FEBU=0.030,

    # ================= rasburicase =================
    V1_RASB=8.0, V2_RASB=6.0, Q_RASB=0.6, CL_RASB=0.30,
    VMAX_RASB=1.50,   # mmol/h per (mg/L)
    KM_RASB=0.030,    # mmol/L -> zero order at TLS urate levels
    KADA=0.004, ADA_CL=2.0,

    # ================= venetoclax / cytotoxics =================
    KA_VEN=0.5, V_VEN=130.0, CL_VEN=3.6, EMAX_VEN=0.060, EC50_VEN=1.0,
    KEL_CHEMO=0.15, EMAX_CHEMO=0.055, EC50_CHEMO=1.0,
    KEL_PRED=0.23, EMAX_PRED=0.020, EC50_PRED=0.10,

    # ================= potassium =================
    CLK0=0.729,       # L/h at GFR0
    GFR_EXP_K=0.70, K_ALDO_GAIN=1.5,
    E_K_COLON=0.40,   # mmol/h GFR-independent
    K_INTAKE=3.561,   # mmol/h, balances baseline excretion
    KSH_INS=0.090, EC50_INS=60.0, KEL_INS=0.55,
    KSH_SALB=0.050, EC50_SALB=8.0, KEL_SALB=0.25,
    KSH_BACK=0.28,
    KSZC=0.55,        # mmol K bound per g binder per h

    # ================= phosphate =================
    TMP0=0.9700,      # mmol/L TmP/GFR (pre-hormonal)
    TMP_PTH=0.45,     # fractional suppression at max PTH
    FGF_MAX=0.40, PO4_SET=0.90, KFGF=1.20,  # FGF23-like phosphate feedback
    F_ULTRA_PO4=0.90,
    PO4_INTAKE=0.95,  # mmol/h net absorption
    KSEVE=0.0022,     # mmol PO4 bound per mg binder per h

    # ================= calcium / PTH =================
    FE_CA=0.012, CA_NET0=0.2030,      # mmol/h net ECF calcium input
    CA_BONE_PTH=2.2,                  # mmol/h extra bone efflux at max PTH
    TAU_PTH=1.0, CAION_SET=1.15, CAION_HI=1.22, KCA_PTH=0.12,
    PTH_MAX=6.0, PTH_SUPP=0.85, CA_HYPERCALCIURIA=3.0,
    CA_PO4_STOICH=1.67,               # Ca:PO4 in hydroxyapatite (10:6)
    CAION_FRAC=0.50, CAION_ALK=0.06,  # alkalosis lowers ionised fraction

    # ================= magnesium / creatinine / urea =================
    FE_MG=0.030, MG_INTAKE=0.1836,
    CR_PROD=0.55, UREA_PROD=18.0, FE_UREA=0.50,

    # ================= crystal / precipitation =================
    CF_MED=1.4,       # medullary concentrating factor on tubular fluid
    S_HU=0.655,       # mmol/L undissociated uric acid solubility
    PKA_UA=5.75, SS_TH_UA=2.5, KN_UA=0.055, N_UA=2.0, KDIS_UA=0.10,
    S_XAN=0.40, PKA_XAN=7.70, SS_TH_XAN=1.4, KN_XAN=0.075, N_XAN=2.0,
    KDIS_XAN=0.05,
    KSP_CAP=6.0, PKA2_PO4=6.80, SS_TH_CAP=4.0, KN_CAP=0.004, N_CAP=2.0,
    KDIS_CAP=0.010,
    KCLR_XT=0.008,    # /h slow clearance of deposited crystal (t1/2 ~3.6 d)
    # systemic (interstitial / metastatic) calcium-phosphate deposition
    KSP_SYS=4.84,     # (mmol/L)^2 == Ca x PO4 of 60 mg2/dL2
    KN_SYS=1.20, N_SYS=2.0, KDIS_SYS=0.004,
    ALK_SYS=0.90,     # extra driving force at full bicarbonate load
    FRAC_KIDNEY=0.20, # share of systemic deposition that is nephrocalcinosis

    # ================= tubular injury / nephron loss =================
    JREF=0.25,        # mmol/h crystal flux giving unit injury drive
    URATE_TOX=0.30, URATE_TOX_TH=0.48, URATE_TOX_C50=0.60,
    KINJ=0.050, KREP=0.030,
    KNL=0.0015, KNR=0.0040,
    XT_K50=12.0, XT_HILL=2.0,

    # ================= urine pH =================
    PH_BASE=5.90, PH_MAX_RISE=1.75, KH_HCO3=140.0, KEL_HCO3=0.17, TAU_PH=0.5,

    # ================= rasburicase oxidant limb =================
    G6PD=1.0, KMETHB=0.0180, KREDMET=0.09, KHEM=0.060, KRECHB=0.004,

    # ================= hazards =================
    H_ARR0=2.0e-5, H_ARR_MAX=0.030, ARR_K=5.5, ARR_KS=1.50,
    ARR_CA=1.05, ARR_CAS=0.25,
    H_SZ0=5.0e-6, H_SZ_MAX=0.010, CAION_SZ=0.90, SZ_S=0.20,
    H_RRT0=1.0e-5, H_RRT_MAX=0.020,
    RRT_K=6.0, RRT_KS=0.80, RRT_CR=2.0, RRT_CRS=1.50,
    RRT_OLIG=1.00, RRT_OLIGS=0.60,
)


def hill(x, k, n=1.0):
    x = max(float(x), 0.0)
    if x <= 0.0:
        return 0.0
    return x ** n / (k ** n + x ** n)


# ---------------------------------------------------------------------------
# 2. REGIMEN
# ---------------------------------------------------------------------------
class Regimen:
    """Continuous infusions and bolus events.  Fluid rate lives in P['FLUID']."""

    def __init__(self, hco3=0.0, ca_inf=0.0, chemo=None, ins_inf=None,
                 dialysis=None, boluses=None, k_intake=None, po4_intake=None):
        self.hco3 = hco3
        self.ca_inf = ca_inf
        self.chemo = chemo            # (t0, t1, rate)
        self.ins_inf = ins_inf
        self.dialysis = dialysis      # (t0, t1, CL L/h)
        self.boluses = boluses or []
        self.k_intake = k_intake
        self.po4_intake = po4_intake

    def inf(self, t, name):
        spec = getattr(self, name)
        if spec is None:
            return 0.0
        t0, t1, rate = spec
        return rate if (t0 <= t < t1) else 0.0


def sched(start, stop, interval, state, amount):
    out, t = [], start
    while t < stop - 1e-9:
        out.append((t, state, amount))
        t += interval
    return out


# ---------------------------------------------------------------------------
# 3. ALGEBRA (identical in rhs and in the reporting layer)
# ---------------------------------------------------------------------------
def derive(y, p):
    g = {s: float(y[i]) for i, s in enumerate(STATES)}
    d = {}
    VEC = max(g["VEC"], 5.0)
    d["VEC"] = VEC

    # ---------------- concentrations ----------------
    d["K"] = g["K_EC"] / VEC
    d["PO4"] = g["PO4_EC"] / VEC
    d["CA_T"] = g["CA_EC"] / VEC
    d["MG"] = g["MG_EC"] / VEC
    d["HX"] = g["HYPOX"] / VEC
    d["XA"] = g["XAN"] / VEC
    d["UA"] = g["URATE"] / VEC
    d["ALLANTC"] = g["ALLANT"] / VEC
    d["CR"] = g["CREAT"] / p["VTBW"]
    d["BUN"] = g["UREA"] / p["VTBW"]
    d["LDH_UL"] = 1000.0 * g["LDH"] / VEC

    alk = hill(g["HCO3_L"], p["KH_HCO3"])
    d["ALK"] = alk
    d["CA_ION"] = d["CA_T"] * (p["CAION_FRAC"] - p["CAION_ALK"] * alk)

    d["UA_MGDL"] = d["UA"] * UA_MGDL
    d["CR_MGDL"] = d["CR"] * CR_MGDL
    d["PO4_MGDL"] = d["PO4"] * PO4_MGDL
    d["CA_MGDL"] = d["CA_T"] * CA_MGDL
    d["CAxPO4"] = d["CA_MGDL"] * d["PO4_MGDL"]

    # ---------------- kidney ----------------
    XT = g["XT_UA"] + g["XT_XAN"] + g["XT_CAP"]
    d["XT_TOT"] = XT
    d["OBSTR"] = hill(XT, p["XT_K50"], p["XT_HILL"])
    d["VOLFAC"] = min(1.15, max(0.40, (VEC / p["VEC0"]) ** 1.5))
    GFR = max(0.02, p["GFR0"] * max(g["NEPH"], 0.0)
              * (1.0 - d["OBSTR"]) * d["VOLFAC"])
    d["GFR"] = GFR
    d["eGFR_mlmin"] = GFR / p["GFR0"] * 120.0

    furo = hill(g["FURO_C"], p["EC50_FURO"])
    d["FURO"] = furo
    qu_t = (p["QU_BASE"] + max(0.0, p["FLUID"] - p["QU_BASE"] - p["INSENS"])
            + p["KVOL"] * (VEC - p["VEC0"]) + p["FURO_QU"] * furo)
    # tubular obstruction limits flow independently of GFR: this is what
    # makes crystal nephropathy oliguric rather than merely azotaemic
    fe_h2o = p["FE_H2O_MAX"] * (1.0 - p["OBSTR_FLOW"] * d["OBSTR"])
    QU = min(max(qu_t, p["QU_FLOOR"]), fe_h2o * GFR)
    d["FE_H2O"] = fe_h2o
    d["QU"] = QU
    d["QU_Lday"] = QU * 24.0
    d["PHU"] = g["PHU"]

    # ---------------- xanthine oxidase inhibition ----------------
    d["C_ALLO"] = g["ALLO_C"] / p["V_ALLO"]
    d["C_OXY"] = g["OXY_C"] / p["V_OXY"]
    d["C_FEBU"] = g["FEBU_C"] / p["V_FEBU"]
    r = (d["C_ALLO"] / p["IC50_ALLO"] + d["C_OXY"] / p["IC50_OXY"]
         + d["C_FEBU"] / p["IC50_FEBU"])
    d["XOI"] = r / (1.0 + r)
    d["XO_FREE"] = 1.0 - d["XOI"]

    # ---------------- purine fluxes ----------------
    d["V_HX_XAN"] = p["VXO"] * d["XO_FREE"] * hill(d["HX"], p["KM_XO1"]) * VEC
    d["V_XAN_UA"] = p["VXO"] * d["XO_FREE"] * hill(d["XA"], p["KM_XO2"]) * VEC
    d["V_SALV"] = p["VSALV"] * hill(d["HX"], p["KM_SALV"]) * VEC
    d["C_RASB"] = g["RASB_C"] / p["V1_RASB"]
    d["V_URIC_OX"] = p["VMAX_RASB"] * d["C_RASB"] * hill(d["UA"], p["KM_RASB"])

    # ---------------- renal solute excretion ----------------
    d["FE_UA_EFF"] = p["FE_UA"] * (1.0 + p["FEUA_GAIN"] * hill(
        max(0.0, d["UA"] - p["FEUA_TH"]), p["FEUA_C50"]))
    d["E_UA"] = d["FE_UA_EFF"] * GFR * d["UA"]
    d["E_HX"] = p["FE_HX"] * GFR * d["HX"]
    d["E_XAN"] = p["FE_XAN"] * GFR * d["XA"]
    d["E_ALLANT"] = p["FE_ALLANT"] * GFR * d["ALLANTC"]

    aldo = 1.0 + p["K_ALDO_GAIN"] * max(0.0, (d["K"] - 4.0) / 4.0)
    clk = (p["CLK0"] * (GFR / p["GFR0"]) ** p["GFR_EXP_K"] * aldo
           * (1.0 + (p["FURO_K"] - 1.0) * furo))
    d["E_K"] = clk * d["K"] + p["E_K_COLON"] * (1.0 + 1.5 * max(0.0, d["K"] - 4.0))

    pth_f = min(1.0, (g["PTH"] - 1.0) / (p["PTH_MAX"] - 1.0))
    d["PTH_F"] = pth_f
    fgf = p["FGF_MAX"] * hill(max(0.0, d["PO4"] - p["PO4_SET"]), p["KFGF"])
    d["FGF"] = fgf
    tmp = p["TMP0"] * (1.0 - p["TMP_PTH"] * max(0.0, pth_f)) * (1.0 - fgf)
    d["TMP"] = tmp
    filt = GFR * d["PO4"] * p["F_ULTRA_PO4"]
    d["E_PO4"] = max(0.0, filt - min(filt, tmp * GFR))

    d["E_CA"] = (p["FE_CA"] * (1.0 + p["CA_HYPERCALCIURIA"]
                 * hill(max(0.0, d["CA_T"] - 2.55), 0.40))
                 * GFR * d["CA_T"])
    d["E_MG"] = p["FE_MG"] * GFR * d["MG"]
    d["E_CR"] = GFR * d["CR"]
    d["E_UREA"] = p["FE_UREA"] * GFR * d["BUN"]

    # ---------------- tubular fluid concentrations ----------------
    cf = p["CF_MED"]
    d["CU_UA"] = cf * d["E_UA"] / QU
    d["CU_XAN"] = cf * d["E_XAN"] / QU
    d["CU_CA"] = cf * d["E_CA"] / QU
    d["CU_PO4"] = cf * d["E_PO4"] / QU

    # ---------------- solubility / supersaturation ----------------
    d["S_UA"] = p["S_HU"] * (1.0 + 10.0 ** (g["PHU"] - p["PKA_UA"]))
    d["SS_UA"] = d["CU_UA"] / d["S_UA"]
    d["S_XAN"] = p["S_XAN"] * (1.0 + 10.0 ** (g["PHU"] - p["PKA_XAN"]))
    d["SS_XAN"] = d["CU_XAN"] / d["S_XAN"]
    d["F_HPO4"] = 1.0 / (1.0 + 10.0 ** (p["PKA2_PO4"] - g["PHU"]))
    d["AP_CAP"] = d["CU_CA"] * d["CU_PO4"] * d["F_HPO4"]
    d["SS_CAP"] = d["AP_CAP"] / p["KSP_CAP"]
    d["SS_SYS"] = (d["CA_T"] * d["PO4"] * (1.0 + p["ALK_SYS"] * alk)) / p["KSP_SYS"]

    # ---------------- precipitation fluxes (mmol/h) ----------------
    # Nucleation/growth rate, then the hard physical bound: a crystal cannot
    # take up more solute per hour than the tubule DELIVERS per hour.  Without
    # this cap, C_urine = E/QU diverges as urine flow goes to zero and the
    # obstructed kidney keeps depositing crystal it is no longer being given.
    n_ua = p["KN_UA"] * QU * max(0.0, d["SS_UA"] - p["SS_TH_UA"]) ** p["N_UA"]
    n_xan = p["KN_XAN"] * QU * max(0.0, d["SS_XAN"] - p["SS_TH_XAN"]) ** p["N_XAN"]
    n_cap = p["KN_CAP"] * QU * max(0.0, d["SS_CAP"] - p["SS_TH_CAP"]) ** p["N_CAP"]
    d["J_UA"] = (min(n_ua, 0.95 * d["E_UA"])
                 - p["KDIS_UA"] * g["XT_UA"] * max(0.0, 1.0 - d["SS_UA"])
                 - p["KCLR_XT"] * g["XT_UA"])
    d["J_XAN"] = (min(n_xan, 0.95 * d["E_XAN"])
                  - p["KDIS_XAN"] * g["XT_XAN"] * max(0.0, 1.0 - d["SS_XAN"])
                  - p["KCLR_XT"] * g["XT_XAN"])
    j_cap_tub = (min(n_cap, 0.95 * d["E_PO4"])
                 - p["KDIS_CAP"] * g["XT_CAP"] * max(0.0, 1.0 - d["SS_CAP"]))
    # deposition and binding must vanish as the phosphate pool empties
    po4_avail = hill(d["PO4"], 0.20)
    d["PO4_AVAIL"] = po4_avail
    ca_avail = hill(d["CA_ION"], 0.35)
    d["CA_AVAIL"] = ca_avail
    d["J_SYS"] = (p["KN_SYS"] * VEC * max(0.0, d["SS_SYS"] - 1.0) ** p["N_SYS"]
                  * po4_avail * ca_avail)
    d["J_CAP_TUB"] = j_cap_tub
    d["J_CAP"] = (j_cap_tub + p["FRAC_KIDNEY"] * d["J_SYS"]
                  - p["KCLR_XT"] * g["XT_CAP"])   # net renal deposition
    d["J_PO4_DEP"] = max(0.0, j_cap_tub) + max(0.0, d["J_SYS"])  # total PO4 lost

    d["DEP_UA"] = min(n_ua, 0.95 * d["E_UA"])
    d["DEP_XAN"] = min(n_xan, 0.95 * d["E_XAN"])
    d["DEP_CAP"] = max(0.0, j_cap_tub) + p["FRAC_KIDNEY"] * d["J_SYS"]
    d["J_XTAL"] = d["DEP_UA"] + d["DEP_XAN"] + d["DEP_CAP"]
    d["INJ_XUA"] = d["DEP_UA"] / p["JREF"]
    d["INJ_XXAN"] = d["DEP_XAN"] / p["JREF"]
    d["INJ_XCAP"] = d["DEP_CAP"] / p["JREF"]
    d["INJ_URATE"] = p["URATE_TOX"] * hill(max(0.0, d["UA"] - p["URATE_TOX_TH"]),
                                           p["URATE_TOX_C50"])
    d["INJ_DRIVE"] = d["J_XTAL"] / p["JREF"] + d["INJ_URATE"]

    # ---------------- hazards ----------------
    # saturating hazards: steep in the pathological range, flat at baseline
    d["H_ARR"] = p["H_ARR0"] + p["H_ARR_MAX"] * hill(
        max(0.0, d["K"] - p["ARR_K"]), p["ARR_KS"], 4.0) * (
        1.0 + hill(max(0.0, p["ARR_CA"] - d["CA_ION"]), p["ARR_CAS"], 2.0))
    d["H_SZ"] = p["H_SZ0"] + p["H_SZ_MAX"] * hill(
        max(0.0, p["CAION_SZ"] - d["CA_ION"]), p["SZ_S"], 3.0)
    # renal-replacement hazard: a saturating function of the three things that
    # actually trigger dialysis -- refractory hyperkalaemia, azotaemia, anuria
    d["H_RRT"] = p["H_RRT0"] + p["H_RRT_MAX"] / 3.0 * (
        hill(max(0.0, d["K"] - p["RRT_K"]), p["RRT_KS"], 3.0)
        + hill(max(0.0, d["CR"] * CR_MGDL / 0.86 - p["RRT_CR"]), p["RRT_CRS"], 3.0)
        + hill(max(0.0, p["RRT_OLIG"] - d["QU_Lday"]), p["RRT_OLIGS"], 3.0))
    return g, d


# ---------------------------------------------------------------------------
# 4. RHS
# ---------------------------------------------------------------------------
def rhs(t, y, p, reg):
    y = np.maximum(y, 0.0)
    g, d = derive(y, p)
    dy = np.zeros(NS)

    # ================= drug PK =================
    dy[IX["ALLO_G"]] = -p["KA_ALLO"] * g["ALLO_G"]
    dy[IX["ALLO_C"]] = (p["KA_ALLO"] * g["ALLO_G"] * p["F_ALLO"]
                        - p["CL_ALLO"] * d["C_ALLO"])
    dy[IX["OXY_C"]] = (p["FM_OXY"] * p["CL_ALLO"] * d["C_ALLO"]
                       * p["MW_OXY"] / p["MW_ALLO"]
                       - p["CL_OXY"] * (d["GFR"] / p["GFR0"]) * d["C_OXY"])
    dy[IX["FEBU_G"]] = -p["KA_FEBU"] * g["FEBU_G"]
    dy[IX["FEBU_C"]] = (p["KA_FEBU"] * g["FEBU_G"] * p["F_FEBU"]
                        - p["CL_FEBU"] * d["C_FEBU"])

    cl_rasb = p["CL_RASB"] * (1.0 + (p["ADA_CL"] - 1.0) * min(g["ADA"], 1.0))
    c_rp = g["RASB_P"] / p["V2_RASB"]
    dy[IX["RASB_C"]] = -cl_rasb * d["C_RASB"] - p["Q_RASB"] * (d["C_RASB"] - c_rp)
    dy[IX["RASB_P"]] = p["Q_RASB"] * (d["C_RASB"] - c_rp)
    dy[IX["ADA"]] = p["KADA"] * d["C_RASB"] - 0.002 * g["ADA"]

    dy[IX["VEN_G"]] = -p["KA_VEN"] * g["VEN_G"]
    c_ven = g["VEN_C"] / p["V_VEN"]
    dy[IX["VEN_C"]] = p["KA_VEN"] * g["VEN_G"] - p["CL_VEN"] * c_ven

    dy[IX["CHEMO_C"]] = reg.inf(t, "chemo") - p["KEL_CHEMO"] * g["CHEMO_C"]
    dy[IX["PRED_C"]] = -p["KEL_PRED"] * g["PRED_C"]
    dy[IX["INS_C"]] = reg.inf(t, "ins_inf") - p["KEL_INS"] * g["INS_C"]
    dy[IX["GLU"]] = 0.55 - 0.09 * g["GLU"] * (1.0 + 0.030 * g["INS_C"])
    dy[IX["SALB_C"]] = -p["KEL_SALB"] * g["SALB_C"]
    dy[IX["FURO_C"]] = -p["KEL_FURO"] * g["FURO_C"]
    dy[IX["SEVE_G"]] = -0.35 * g["SEVE_G"]
    dy[IX["SZC_G"]] = -0.35 * g["SZC_G"]
    dy[IX["HCO3_L"]] = reg.hco3 - p["KEL_HCO3"] * g["HCO3_L"]

    # ================= tumour kill / release =================
    kkill = (p["KKILL_SP"]
             + p["EMAX_VEN"] * hill(c_ven, p["EC50_VEN"])
             + p["EMAX_CHEMO"] * hill(g["CHEMO_C"], p["EC50_CHEMO"])
             + p["EMAX_PRED"] * hill(g["PRED_C"], p["EC50_PRED"]))
    dy[IX["N_VIA"]] = (math.log(2.0) / p["TD"]) * g["N_VIA"] - kkill * g["N_VIA"]
    dy[IX["N_LYS"]] = kkill * g["N_VIA"] - p["KLYS"] * g["N_LYS"]
    rel = p["KLYS"] * g["N_LYS"]
    dy[IX["LDH"]] = p["Q_LDH"] * rel + p["LDH_PROD"] - p["KEL_LDH"] * g["LDH"]

    cl_hd = reg.inf(t, "dialysis")

    # ================= potassium =================
    shift_in = (p["KSH_INS"] * hill(g["INS_C"], p["EC50_INS"])
                + p["KSH_SALB"] * hill(g["SALB_C"], p["EC50_SALB"])) * g["K_EC"]
    j_shift = shift_in - p["KSH_BACK"] * g["K_SH"]
    j_szc = p["KSZC"] * g["SZC_G"] * min(1.0, d["K"] / 4.0)
    k_in = p["K_INTAKE"] if reg.k_intake is None else reg.k_intake
    dy[IX["K_EC"]] = (p["Q_K"] * rel + k_in - d["E_K"] - j_shift - j_szc
                      - cl_hd * d["K"])
    dy[IX["K_SH"]] = j_shift

    # ================= phosphate =================
    po4_in = p["PO4_INTAKE"] if reg.po4_intake is None else reg.po4_intake
    dy[IX["PO4_EC"]] = (p["Q_PO4"] * rel + po4_in - d["E_PO4"]
                        - p["KSEVE"] * g["SEVE_G"] * d["PO4_AVAIL"]
                        - d["J_PO4_DEP"]
                        - cl_hd * d["PO4"] * 0.9)

    # ================= calcium / PTH =================
    dy[IX["CA_EC"]] = (p["CA_NET0"] + p["CA_BONE_PTH"] * d["PTH_F"] + reg.ca_inf
                       - d["E_CA"] - p["CA_PO4_STOICH"] * d["J_PO4_DEP"])
    # PTH responds only to ionised calcium BELOW its set point, so that the
    # basal state is PTH = 1 with zero net PTH-driven bone efflux.
    pth_ss = (1.0 + (p["PTH_MAX"] - 1.0)
              * hill(max(0.0, p["CAION_SET"] - d["CA_ION"]), p["KCA_PTH"], 2.0)
              - p["PTH_SUPP"]
              * hill(max(0.0, d["CA_ION"] - p["CAION_HI"]), p["KCA_PTH"], 2.0))
    dy[IX["PTH"]] = (pth_ss - g["PTH"]) / p["TAU_PTH"]

    dy[IX["MG_EC"]] = (p["Q_MG"] * rel + p["MG_INTAKE"] - d["E_MG"]
                       - cl_hd * d["MG"])

    # ================= purines =================
    dy[IX["HYPOX"]] = (p["Q_PUR"] * rel + p["P_END"] - d["V_HX_XAN"]
                       - d["V_SALV"] - d["E_HX"] - cl_hd * d["HX"])
    dy[IX["XAN"]] = (d["V_HX_XAN"] - d["V_XAN_UA"] - d["E_XAN"]
                     - max(0.0, d["J_XAN"]) - cl_hd * d["XA"])
    dy[IX["URATE"]] = (d["V_XAN_UA"] - d["E_UA"] - d["V_URIC_OX"]
                       - max(0.0, d["J_UA"]) - cl_hd * d["UA"])
    dy[IX["ALLANT"]] = d["V_URIC_OX"] - d["E_ALLANT"] - cl_hd * d["ALLANTC"]

    dy[IX["CREAT"]] = p["CR_PROD"] + p["Q_CR"] * rel - d["E_CR"] - cl_hd * d["CR"]
    dy[IX["UREA"]] = p["UREA_PROD"] - d["E_UREA"] - cl_hd * d["BUN"] * 1.2

    # ================= kidney =================
    dy[IX["TUBINJ"]] = (p["KINJ"] * d["INJ_DRIVE"] * (1.0 - g["TUBINJ"])
                        - p["KREP"] * g["TUBINJ"])
    dy[IX["NEPH"]] = (-p["KNL"] * g["TUBINJ"] ** 2 * g["NEPH"]
                      + p["KNR"] * (1.0 - g["NEPH"])
                      * max(0.0, 1.0 - g["TUBINJ"] / 0.30))
    dy[IX["XT_UA"]] = d["J_UA"]
    dy[IX["XT_XAN"]] = d["J_XAN"]
    dy[IX["XT_CAP"]] = d["J_CAP"]

    ph_ss = p["PH_BASE"] + p["PH_MAX_RISE"] * hill(g["HCO3_L"], p["KH_HCO3"])
    dy[IX["PHU"]] = (ph_ss - g["PHU"]) / p["TAU_PH"]
    dy[IX["VEC"]] = p["FLUID"] - d["QU"] - p["INSENS"]

    # ================= oxidant limb =================
    dy[IX["H2O2"]] = d["V_URIC_OX"]
    ox = d["V_URIC_OX"] * (1.0 - p["G6PD"])
    dy[IX["METHB"]] = p["KMETHB"] * ox - p["KREDMET"] * g["METHB"]
    dy[IX["HB"]] = -p["KHEM"] * ox + p["KRECHB"] * (p["HB0"] - g["HB"])

    # ================= hazards =================
    dy[IX["CH_ARR"]] = d["H_ARR"]
    dy[IX["CH_SZ"]] = d["H_SZ"]
    dy[IX["CH_RRT"]] = d["H_RRT"]
    dy[IX["CUM_REL"]] = rel
    return dy


# ---------------------------------------------------------------------------
# 5. SIMULATION
# ---------------------------------------------------------------------------
def y0_default(p, N0=None):
    y = np.zeros(NS)
    y[IX["VEC"]] = p["VEC0"]
    y[IX["N_VIA"]] = p["N0"] if N0 is None else N0
    y[IX["K_EC"]] = 4.10 * p["VEC0"]
    y[IX["PO4_EC"]] = 1.15 * p["VEC0"]
    y[IX["CA_EC"]] = 2.35 * p["VEC0"]
    y[IX["MG_EC"]] = 0.85 * p["VEC0"]
    y[IX["HYPOX"]] = 0.0004 * p["VEC0"]
    y[IX["XAN"]] = 0.0004 * p["VEC0"]
    y[IX["URATE"]] = 0.32 * p["VEC0"]
    y[IX["CREAT"]] = 0.076 * p["VTBW"]
    y[IX["UREA"]] = 5.0 * p["VTBW"]
    y[IX["NEPH"]] = 1.0
    y[IX["PHU"]] = p["PH_BASE"]
    y[IX["HB"]] = p["HB0"]
    y[IX["PTH"]] = 1.0
    y[IX["LDH"]] = 0.200 * p["VEC0"]
    y[IX["GLU"]] = 6.1
    return y


BASE = None


def equilibrate(p, hours=600.0):
    reg = Regimen()
    y = y0_default(p, N0=1e-9)
    sol = solve_ivp(rhs, (0.0, hours), y, args=(p, reg), method="LSODA",
                    rtol=1e-8, atol=1e-10, max_step=4.0)
    return sol.y[:, -1]


def simulate(reg, tmax=336.0, N0=None, y0=None, dt=0.25, overrides=None):
    p = dict(P)
    if overrides:
        p.update(overrides)
    if y0 is None:
        y0 = (BASE.copy() if BASE is not None else y0_default(p))
    y = np.array(y0, dtype=float).copy()
    if N0 is not None:
        y[IX["N_VIA"]] = N0

    ev = sorted(reg.boluses, key=lambda b: b[0])
    breaks = sorted({0.0, tmax} | {b[0] for b in ev if 0.0 < b[0] < tmax})
    # Events before t = 0 (lead-time prophylaxis, steroid prephase) are handled
    # by pre-integration.  N0 is defined as the burden AT t = 0 in the absence
    # of pre-treatment, so the tumour is wound back along its own growth curve
    # first -- otherwise a 5-day prephase arm would silently start from a
    # tumour 2^(120/30) = 16-fold larger and the comparison would be rigged.
    pre = [b for b in ev if b[0] < 0.0]
    if pre:
        t_start = min(b[0] for b in pre)
        y[IX["N_VIA"]] *= 2.0 ** (t_start / p["TD"])
        pbreaks = sorted({t_start, 0.0} | {b[0] for b in pre if b[0] > t_start})
        for i in range(len(pbreaks) - 1):
            t0, t1 = pbreaks[i], pbreaks[i + 1]
            for (tb, st, amt) in pre:
                if abs(tb - t0) < 1e-9:
                    y[IX[st]] += amt
            sol = solve_ivp(rhs, (t0, t1), y, args=(p, reg), method="LSODA",
                            rtol=1e-6, atol=1e-8, max_step=1.0)
            y = sol.y[:, -1].copy()
    ts, ys = [], []
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        for (tb, st, amt) in ev:
            if abs(tb - t0) < 1e-9 and tb >= 0.0:
                y[IX[st]] += amt
        n = max(2, int(round((t1 - t0) / dt)) + 1)
        sol = solve_ivp(rhs, (t0, t1), y, args=(p, reg), method="LSODA",
                        t_eval=np.linspace(t0, t1, n),
                        rtol=1e-6, atol=1e-8, max_step=1.0)
        if not sol.success:
            raise RuntimeError(f"integration failed: {sol.message}")
        ts.append(sol.t if not ts else sol.t[1:])
        ys.append(sol.y if not ys else sol.y[:, 1:])
        y = sol.y[:, -1].copy()
    return np.concatenate(ts), np.hstack(ys), p


def tabulate(T, Y, p):
    n = Y.shape[1]
    _, d0 = derive(Y[:, 0], p)
    out = {k: np.zeros(n) for k in d0}
    for j in range(n):
        _, d = derive(Y[:, j], p)
        for k in out:
            out[k][j] = d[k]
    for s in STATES:
        out[s] = Y[IX[s], :].copy()
    out["t"] = T
    return out


# ---------------------------------------------------------------------------
# 6. SCENARIO BUILDERS
# ---------------------------------------------------------------------------
HYD_STD = 0.113                    # 2 L/day urine
HYD_AGGR = 3.0 * 1.73 / 24.0       # 3 L/m2/day = 0.216 L/h
HYD_MAX = 4.0 * 1.73 / 24.0        # 4 L/m2/day


def burkitt(hydration=HYD_AGGR, allo_start=None, allo_dose=300.0,
            febu_start=None, rasb=None, rasb_dose_mgkg=0.20, rasb_days=1,
            hco3=0.0, N0=3.0, chemo=(0.0, 24.0, 1.0), tmax=336.0,
            ca_inf=0.0, furo=None, dialysis=None, insulin=None,
            szc=None, seve=None, pred_prephase=None, overrides=None):
    bol = []
    if allo_start is not None:
        bol += sched(allo_start, tmax, 24.0, "ALLO_G", allo_dose)
    if febu_start is not None:
        bol += sched(febu_start, tmax, 24.0, "FEBU_G", 120.0)
    if rasb is not None:
        for k in range(rasb_days):
            bol.append((rasb + 24.0 * k, "RASB_C", rasb_dose_mgkg * P["BW"]))
    if szc is not None:
        bol += sched(szc, tmax, 8.0, "SZC_G", 10.0)
    if seve is not None:
        bol += sched(seve, tmax, 8.0, "SEVE_G", 1600.0)
    if pred_prephase is not None:
        t0, nd = pred_prephase
        bol += sched(t0, t0 + 24.0 * nd, 24.0, "PRED_C", 1.2)
    if furo is not None:
        bol += sched(furo, tmax, 8.0, "FURO_C", 1.4)
    ov = dict(FLUID=hydration)
    if overrides:
        ov.update(overrides)
    reg = Regimen(hco3=hco3, ca_inf=ca_inf, chemo=chemo, ins_inf=insulin,
                  dialysis=dialysis, boluses=bol)
    return reg, N0, tmax, ov


def run(spec, y0=None, dt=0.25):
    reg, N0, tmax, ov = spec
    T, Y, p = simulate(reg, tmax=tmax, N0=N0, overrides=ov, y0=y0, dt=dt)
    return tabulate(T, Y, p)


# ---------------------------------------------------------------------------
# 7. SUMMARY / CLASSIFICATION
# ---------------------------------------------------------------------------
def summarise(r, label=""):
    t = r["t"]
    m = t >= 0.0
    cr0 = float(np.interp(0.0, t, r["CR_MGDL"]))
    s = dict(
        label=label,
        UA_peak=r["UA_MGDL"][m].max(), K_peak=r["K"][m].max(),
        PO4_peak=r["PO4"][m].max(), CA_nadir=r["CA_ION"][m].min(),
        CR_peak=r["CR_MGDL"][m].max(), CR_ratio=r["CR_MGDL"][m].max() / cr0,
        CR_ratio7=(r["CR_MGDL"][(t >= 0.0) & (t <= 168.0)].max() / cr0),
        eGFR_nadir7=r["eGFR_mlmin"][(t >= 0.0) & (t <= 168.0)].min(),
        eGFR_nadir=r["eGFR_mlmin"][m].min(),
        XT_UA=r["XT_UA"][-1], XT_XAN=r["XT_XAN"][-1], XT_CAP=r["XT_CAP"][-1],
        XT_TOT=r["XT_TOT"][m].max(), NEPH_final=r["NEPH"][-1],
        P_ARR=1.0 - math.exp(-r["CH_ARR"][-1]),
        P_SZ=1.0 - math.exp(-r["CH_SZ"][-1]),
        P_RRT=1.0 - math.exp(-r["CH_RRT"][-1]),
        LDH_peak=r["LDH_UL"][m].max(), XAN_peak=r["XA"][m].max(),
        H2O2=r["H2O2"][-1], MetHb=r["METHB"].max() * 100.0,
        Hb_nadir=r["HB"].min(), QU_min=r["QU_Lday"][m].min(),
        pH=r["PHU"].max(), TUBINJ=r["TUBINJ"][m].max(),
        lysed=r["CUM_REL"][-1], lysed_pre=float(np.interp(0.0, t, r["CUM_REL"])),
    )
    s["LTLS"], s["CTLS"], s["ncrit"] = cairo_bishop(r)
    return s


def cairo_bishop(r, sample_h=None, window=(-72.0, 168.0)):
    """
    Cairo & Bishop 2004 definition.
    Laboratory TLS: >= 2 of the four, from 3 days before to 7 days after therapy
        urate >= 476 umol/L (8.0 mg/dL) or +25% from baseline
        K     >= 6.0 mmol/L            or +25%
        PO4   >= 1.45 mmol/L (adult)   or +25%
        Ca    <= 1.75 mmol/L total     or -25%
    Clinical TLS: LTLS + creatinine >= 1.5 x ULN, or arrhythmia, or seizure.
    `sample_h` evaluates the criteria only on a discrete blood-draw schedule,
    which is what actually happens on a ward.
    """
    t = r["t"]
    lo, hi = max(t[0], window[0]), min(t[-1], window[1])
    if sample_h is None:
        idx = np.where((t >= lo) & (t <= hi))[0]
    else:
        want = np.arange(lo, hi + 1e-9, sample_h)
        idx = np.array([int(np.argmin(np.abs(t - w))) for w in want])
    if idx.size == 0:
        return False, False, 0
    i0 = int(np.argmin(np.abs(t - 0.0)))
    ua0, k0, po40, ca0 = (r["UA_MGDL"][i0], r["K"][i0], r["PO4"][i0], r["CA_T"][i0])
    c_ua = bool(np.any((r["UA_MGDL"][idx] >= 8.0) | (r["UA_MGDL"][idx] >= 1.25 * ua0)))
    c_k = bool(np.any((r["K"][idx] >= 6.0) | (r["K"][idx] >= 1.25 * k0)))
    c_po4 = bool(np.any((r["PO4"][idx] >= 1.45) | (r["PO4"][idx] >= 1.25 * po40)))
    c_ca = bool(np.any((r["CA_T"][idx] <= 1.75) | (r["CA_T"][idx] <= 0.75 * ca0)))
    n = int(c_ua) + int(c_k) + int(c_po4) + int(c_ca)
    ltls = n >= 2
    cr_ok = bool(np.any(r["CR_MGDL"][idx] >= 1.5 * 0.86))
    arr = (1.0 - math.exp(-r["CH_ARR"][-1])) > 0.05
    sz = (1.0 - math.exp(-r["CH_SZ"][-1])) > 0.05
    return ltls, bool(ltls and (cr_ok or arr or sz)), n


def hdr(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


def row(*cells, w=None):
    w = w or [26] + [10] * (len(cells) - 1)
    print("".join(str(c).ljust(x) if i == 0 else str(c).rjust(x)
                  for i, (c, x) in enumerate(zip(cells, w))))


# ---------------------------------------------------------------------------
# 7b. CLOSED-FORM THRESHOLDS -- the two concentrations that define the race
# ---------------------------------------------------------------------------
def fe_ua_of(ua, p=None):
    p = p or P
    return p["FE_UA"] * (1.0 + p["FEUA_GAIN"]
                         * hill(max(0.0, ua - p["FEUA_TH"]), p["FEUA_C50"]))


def ua_required(jrel, gfr=None, p=None):
    """Plasma urate needed for renal excretion to balance a release flux."""
    p = p or P
    gfr = gfr if gfr is not None else p["GFR0"]
    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if fe_ua_of(mid, p) * gfr * mid < jrel:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ua_critical(qu, ph=None, gfr=None, p=None):
    """Plasma urate at which tubular fluid reaches the metastable limit."""
    p = p or P
    gfr = gfr if gfr is not None else p["GFR0"]
    ph = ph if ph is not None else p["PH_BASE"]
    s_ua = p["S_HU"] * (1.0 + 10.0 ** (ph - p["PKA_UA"]))
    target = p["SS_TH_UA"] * s_ua
    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        cu = p["CF_MED"] * fe_ua_of(mid, p) * gfr * mid / qu
        if cu < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def peak_release_flux(N0, q=None, p=None, kkill=None):
    """
    Peak lysis flux of solute q (mmol/h) for a burden of N0 x 1e12 cells,
    in closed form.  With N' = -(k-g)N and L' = kN - kl*L,
        L(t) = k N0 /(kl-a) (e^-at - e^-kl t),  a = k - g
    which peaks at t* = ln(kl/a)/(kl-a).  Release flux = q * kl * L(t*).
    """
    p = p or P
    q = p["Q_PUR"] if q is None else q
    k = p["EMAX_CHEMO"] if kkill is None else kkill
    a = k - math.log(2.0) / p["TD"]
    kl = p["KLYS"]
    if a <= 1e-9:
        return q * k * N0                       # no net depletion
    tstar = math.log(kl / a) / (kl - a)
    L = k * N0 / (kl - a) * (math.exp(-a * tstar) - math.exp(-kl * tstar))
    return q * kl * L


def k_capacity(kconc, gfr=None, p=None):
    """Renal + colonic potassium excretion capacity at a given serum K."""
    p = p or P
    gfr = gfr if gfr is not None else p["GFR0"]
    aldo = 1.0 + p["K_ALDO_GAIN"] * max(0.0, (kconc - 4.0) / 4.0)
    clk = p["CLK0"] * (gfr / p["GFR0"]) ** p["GFR_EXP_K"] * aldo
    return clk * kconc + p["E_K_COLON"] * (1.0 + 1.5 * max(0.0, kconc - 4.0))


# ---------------------------------------------------------------------------
# 8. ANALYSES
# ---------------------------------------------------------------------------
def A0_baseline():
    hdr("A0. BASELINE STEADY STATE (no tumour, no drug) -- calibration check")
    _, d = derive(BASE, P)
    checks = [
        ("urate", "UA_MGDL", "mg/dL", "%.2f", "4-6"),
        ("potassium", "K", "mmol/L", "%.2f", "3.5-5.0"),
        ("phosphate", "PO4", "mmol/L", "%.2f", "0.8-1.45"),
        ("total calcium", "CA_T", "mmol/L", "%.2f", "2.2-2.6"),
        ("ionised calcium", "CA_ION", "mmol/L", "%.2f", "1.12-1.30"),
        ("creatinine", "CR_MGDL", "mg/dL", "%.2f", "0.7-1.2"),
        ("urea", "BUN", "mmol/L", "%.1f", "3-7"),
        ("eGFR", "eGFR_mlmin", "mL/min", "%.0f", "~120"),
        ("urine output", "QU_Lday", "L/day", "%.2f", "1.5-2.5"),
        ("urine pH", "PHU", "-", "%.2f", "5.5-6.5"),
        ("LDH", "LDH_UL", "U/L", "%.0f", "140-280"),
        ("urine urate", "CU_UA", "mmol/L", "%.2f", "2-4"),
        ("urate solubility", "S_UA", "mmol/L", "%.2f", "1.5-3"),
        ("urate supersat.", "SS_UA", "-", "%.2f", "1.5-2.5 (metastable)"),
        ("urine Ca", "CU_CA", "mmol/L", "%.2f", "2-5"),
        ("urine PO4", "CU_PO4", "mmol/L", "%.1f", "15-30"),
        ("tubular CaP supersat.", "SS_CAP", "-", "%.2f", "<4 (no nucleation)"),
        ("systemic CaxPO4 SS", "SS_SYS", "-", "%.2f", "<1"),
        ("tubular injury", "TUBINJ", "-", "%.3f", "0"),
        ("nephron mass", "NEPH", "-", "%.3f", "1.0"),
    ]
    g = {s: BASE[IX[s]] for s in STATES}
    for name, k, unit, fmt, ref in checks:
        v = d[k] if k in d else g[k]
        print(f"  {name:24s} {fmt % v:>9s} {unit:9s} target {ref}")
    print()
    print(f"  urate excretion  {d['E_UA'] * 24 * 168.1:.0f} mg/day"
          f"        literature ~700")
    print(f"  K excretion      {d['E_K'] * 24:.0f} mmol/day      literature 70-100")
    print(f"  PO4 excretion    {d['E_PO4'] * 24:.0f} mmol/day      literature 20-30")
    print(f"  Ca excretion     {d['E_CA'] * 24:.1f} mmol/day     literature 2.5-7.5")


def A1_race():
    hdr("A1. THE RACE -- the two urate concentrations that define the disease")
    print("  UA_req  = the plasma urate at which renal excretion balances the")
    print("            peak purine release flux (fractional excretion rises")
    print("            with urate because URAT1 reabsorption saturates).")
    print("  UA_crit = the plasma urate at which tubular fluid reaches the")
    print("            metastable limit and urate begins to precipitate.")
    print("  The disease is the region where UA_req > UA_crit: the kidney")
    print("  cannot excrete the load without crystallising it.")
    print(f"\n  UA_crit is a property of the PRESCRIPTION, not the tumour:")
    for qu, lab in [(HYD_STD - P["INSENS"], "2 L/day, pH 5.9"),
                    (HYD_AGGR - P["INSENS"], "3 L/m2/day, pH 5.9"),
                    (HYD_MAX - P["INSENS"], "4 L/m2/day, pH 5.9"),
                    (HYD_AGGR - P["INSENS"], "3 L/m2/day, pH 7.0")]:
        ph = 7.0 if "7.0" in lab else P["PH_BASE"]
        uc = ua_critical(qu, ph=ph)
        print(f"    {lab:22s} UA_crit = {uc:.2f} mmol/L "
              f"= {uc * UA_MGDL:5.1f} mg/dL")
    print(f"\n  (the Cairo-Bishop laboratory-TLS urate cut-off is 8.0 mg/dL, "
          f"which\n   this model reproduces as a crystallisation threshold "
          f"rather than\n   as a chosen number)")
    print()
    row("burden 1e12", "urate g", "Jrel", "UA_req", "UA_crit", "req/crit",
        "UApk", "Kpk", "CRrat", "TLS", w=[13, 8, 7, 8, 9, 9, 7, 7, 7, 5])
    qu = HYD_AGGR - P["INSENS"]
    uc = ua_critical(qu)
    for N0 in [0.05, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
        r = run(burkitt(N0=N0))
        sm = summarise(r)
        jrel = peak_release_flux(N0)
        ur = ua_required(jrel)
        row(f"{N0:.2f}", f"{N0 * P['Q_PUR'] * 0.168:.1f}", f"{jrel:.2f}",
            f"{ur * UA_MGDL:.1f}", f"{uc * UA_MGDL:.1f}", f"{ur / uc:.2f}",
            f"{sm['UA_peak']:.1f}", f"{sm['K_peak']:.2f}", f"{sm['CR_ratio']:.2f}",
            "C" if sm["CTLS"] else ("L" if sm["LTLS"] else "-"),
            w=[13, 8, 7, 8, 9, 9, 7, 7, 7, 5])
    print("\n  UA_req and UA_crit in mg/dL.  Note where req/crit crosses 1 and")
    print("  compare it with where the TLS column changes: the crossing is not")
    print("  fitted to it, both come out of the same two equations.")
    print("\n  The same race for the other two solutes (peak flux vs the")
    print("  clearance the kidney has at the Cairo-Bishop threshold):")
    row("solute", "content/1e12", "Jrel at N=3", "capacity", "ratio",
        w=[12, 14, 13, 11, 8])
    for name, q, cap in [
            ("urate", P["Q_PUR"], fe_ua_of(0.476) * P["GFR0"] * 0.476),
            ("potassium", P["Q_K"], k_capacity(6.0)),
            ("phosphate", P["Q_PO4"], 12.0)]:
        j = peak_release_flux(3.0, q=q)
        row(name, f"{q:.1f} mmol", f"{j:.2f} mmol/h", f"{cap:.2f}",
            f"{j / cap:.1f}", w=[12, 14, 13, 11, 8])
    print("  (phosphate capacity 12 mmol/h = maximal phosphaturia at intact GFR)")
    print("\n  Consistency check -- closed-form UA_req against the integrated")
    print("  ODE peak, which are computed by completely different routes:")
    row("burden", "UA_req", "UApk (ODE)", "ratio", w=[10, 10, 13, 9])
    for N0 in [0.5, 1.0, 2.0, 3.0]:
        ur = ua_required(peak_release_flux(N0)) * UA_MGDL
        up = summarise(run(burkitt(N0=N0)))["UA_peak"]
        row(f"{N0:.1f}", f"{ur:.1f}", f"{up:.1f}", f"{up / ur:.2f}",
            w=[10, 10, 13, 9])
    print("  (they diverge above ~3e12 cells because the ODE loses GFR, which")
    print("   the closed form assumes intact -- the divergence IS the loop)")


def A2_operators():
    hdr("A2. OPERATOR DECOMPOSITION -- which term of the loop does each touch?")
    N0 = 3.0
    arms = [
        ("nothing (2 L/day)", burkitt(N0=N0, hydration=HYD_STD)),
        ("DILUTION 3 L/m2/day", burkitt(N0=N0, hydration=HYD_AGGR)),
        ("DILUTION 4 L/m2/day", burkitt(N0=N0, hydration=HYD_MAX)),
        ("DILUTION + furosemide", burkitt(N0=N0, furo=0.0)),
        ("FLUX allopurinol t=0", burkitt(N0=N0, allo_start=0.0)),
        ("FLUX febuxostat t=0", burkitt(N0=N0, febu_start=0.0)),
        ("POOL rasburicase", burkitt(N0=N0, rasb=0.0)),
        ("POOL dialysis d2-5", burkitt(N0=N0, dialysis=(48.0, 120.0, 1.2))),
        ("SPECIATION pH 7.5", burkitt(N0=N0, hco3=45.0)),
        ("SEQUESTRATION PO4 bind", burkitt(N0=N0, seve=0.0)),
        ("FLUX-SHAPE prephase", burkitt(N0=N0, pred_prephase=(-120.0, 5))),
    ]
    row("arm", "UApk", "Kpk", "PO4pk", "Ca_ion", "CRrat", "Xua", "Xcap",
        "P(RRT)", w=[24, 8, 8, 8, 8, 8, 8, 8, 9])
    for lab, spec in arms:
        s = summarise(run(spec))
        row(lab, f"{s['UA_peak']:.1f}", f"{s['K_peak']:.2f}", f"{s['PO4_peak']:.2f}",
            f"{s['CA_nadir']:.2f}", f"{s['CR_ratio']:.2f}", f"{s['XT_UA']:.2f}",
            f"{s['XT_CAP']:.2f}", f"{100 * s['P_RRT']:.1f}%",
            w=[24, 8, 8, 8, 8, 8, 8, 8, 9])
    print("\n  Xua / Xcap = renal crystal mass (mmol) of urate and calcium")
    print("  phosphate.  Read the two crystal columns against each other.")


def A3_ph_optimum():
    hdr("A3. URINE pH -- the trade the alkalinisation era never priced")
    print("  Uric acid solubility rises with pH (pKa 5.75), so alkali raises")
    print("  UA_crit steeply.  But the calcium-phosphate driving force rises")
    print("  with pH too, because the precipitating species is HPO4(2-)")
    print("  (pKa2 6.8), and systemic alkalosis simultaneously lowers ionised")
    print("  calcium.  Both curves are printed, then what results.")
    print()
    row("urine pH", "S_UA", "UA_crit", "fHPO4", "vs pH 5.9",
        w=[10, 9, 10, 10, 24])
    qu = HYD_AGGR - P["INSENS"]
    f0 = 1 / (1 + 10 ** (P["PKA2_PO4"] - P["PH_BASE"]))
    s0 = P["S_HU"] * (1 + 10 ** (P["PH_BASE"] - P["PKA_UA"]))
    for ph in [5.5, 5.9, 6.5, 7.0, 7.5, 7.8]:
        sua = P["S_HU"] * (1 + 10 ** (ph - P["PKA_UA"]))
        fh = 1 / (1 + 10 ** (P["PKA2_PO4"] - ph))
        row(f"{ph:.1f}", f"{sua:.2f}", f"{ua_critical(qu, ph=ph) * UA_MGDL:.1f}",
            f"{fh:.3f}", f"UA x{sua / s0:.1f} / CaP x{fh / f0:.1f}",
            w=[10, 9, 10, 10, 24])
    print("\n  Over pH 5.9 -> 7.5 the urate term improves 8-fold and the")
    print("  calcium-phosphate term worsens 7-fold.  Which one wins is not a")
    print("  matter of opinion, it depends on which solute is rate-limiting,")
    print("  i.e. on the tumour burden and on whether rasburicase is given.")
    for N0, tag in [(3.0, "burden 3e12 (urate rate-limiting)"),
                    (6.0, "burden 6e12 (phosphate co-limiting)")]:
        for rasb in [None, 0.0]:
            print(f"\n  --- {tag}, rasburicase: "
                  f"{'no' if rasb is None else 'yes 0.2 mg/kg x3'} ---")
            row("urine pH", "Xua", "Xcap", "Ca_ion", "CRrat", "P(RRT)", "P(sz)",
                w=[10, 9, 9, 9, 8, 9, 8])
            best = None
            for hco3 in [0.0, 16.0, 40.0, 90.0, 200.0]:
                r = run(burkitt(N0=N0, hco3=hco3, rasb=rasb, rasb_days=3,
                                hydration=HYD_STD if N0 > 3.0 else HYD_AGGR))
                sm = summarise(r)
                row(f"{r['PHU'][-1]:.2f}", f"{sm['XT_UA']:.2f}",
                    f"{sm['XT_CAP']:.2f}", f"{sm['CA_nadir']:.2f}",
                    f"{sm['CR_ratio']:.2f}", f"{100 * sm['P_RRT']:.1f}%",
                    f"{100 * sm['P_SZ']:.1f}%", w=[10, 9, 9, 9, 8, 9, 8])
                if best is None or sm["CR_ratio"] < best[1]:
                    best = (float(r["PHU"][-1]), sm["CR_ratio"])
            print(f"    optimum pH {best[0]:.2f}  (peak Cr ratio {best[1]:.2f})")


def A4_leadtime():
    hdr("A4. ALLOPURINOL LEAD TIME -- a flux operator cannot empty a pool")
    print("  Identical 300 mg/day allopurinol; only the START TIME relative to")
    print("  cytotoxic therapy (t = 0) changes.")
    row("start", "XOI@0h", "UApk", "XANpk", "Xua", "Xxan", "CRrat", "P(RRT)",
        w=[12, 9, 8, 9, 8, 8, 8, 9])
    for st in [None, 0.0, -12.0, -24.0, -48.0, -72.0, -120.0, -168.0]:
        r = run(burkitt(N0=3.0, allo_start=st))
        s = summarise(r)
        i0 = int(np.argmin(np.abs(r["t"] - 0.0)))
        row("none" if st is None else f"{st:+.0f} h", f"{r['XOI'][i0]:.3f}",
            f"{s['UA_peak']:.1f}", f"{s['XAN_peak'] * 15.21:.1f}",
            f"{s['XT_UA']:.2f}", f"{s['XT_XAN']:.2f}", f"{s['CR_ratio']:.2f}",
            f"{100 * s['P_RRT']:.1f}%", w=[12, 9, 8, 9, 8, 8, 8, 9])
    print("\n  XANpk in mg/dL.  Note the DIRECTION of the xanthine column: the")
    print("  better the flux block, the more xanthine there is to precipitate,")
    print("  and xanthine is less soluble than the urate it replaced.")


def A5_pool_vs_flux():
    hdr("A5. POOL vs FLUX at matched start time (both begun at t = 0)")
    arms = [("nothing", burkitt(N0=3.0)),
            ("allopurinol 300", burkitt(N0=3.0, allo_start=0.0)),
            ("allopurinol 600", burkitt(N0=3.0, allo_start=0.0, allo_dose=600.0)),
            ("febuxostat 120", burkitt(N0=3.0, febu_start=0.0)),
            ("rasburicase 0.2 x1", burkitt(N0=3.0, rasb=0.0)),
            ("rasburicase 0.2 x3", burkitt(N0=3.0, rasb=0.0, rasb_days=3)),
            ("rasburicase 0.15 x1", burkitt(N0=3.0, rasb=0.0, rasb_dose_mgkg=0.15)),
            ("febuxostat + rasb", burkitt(N0=3.0, febu_start=0.0, rasb=0.0))]
    row("arm", "UA 0h", "UA 4h", "%change", "UApk", "AUC0-96", "CRrat", "P(RRT)",
        w=[22, 8, 8, 10, 8, 10, 8, 9])
    for lab, spec in arms:
        r = run(spec)
        s = summarise(r)
        i0 = int(np.argmin(np.abs(r["t"] - 0.0)))
        i4 = int(np.argmin(np.abs(r["t"] - 4.0)))
        m = (r["t"] >= 0.0) & (r["t"] <= 96.0)
        auc = float(np.trapezoid(r["UA_MGDL"][m], r["t"][m]))
        ch = 100.0 * (r["UA_MGDL"][i4] - r["UA_MGDL"][i0]) / r["UA_MGDL"][i0]
        row(lab, f"{r['UA_MGDL'][i0]:.2f}", f"{r['UA_MGDL'][i4]:.2f}",
            f"{ch:+.0f}%", f"{s['UA_peak']:.1f}", f"{auc:.0f}",
            f"{s['CR_ratio']:.2f}", f"{100 * s['P_RRT']:.1f}%",
            w=[22, 8, 8, 10, 8, 10, 8, 9])
    print("\n  %change at 4 h is the registration-trial endpoint")
    print("  (reported: rasburicase -86%, allopurinol +2%; Goldman 2001).")


def A5b_capacity():
    hdr("A5b. RASBURICASE IS A ZERO-ORDER POOL OPERATOR WITH A FIXED CAPACITY")
    print("  The Km of urate oxidase for urate (~25 umol/L) is two orders of")
    print("  magnitude below the urate concentrations of TLS, so the enzyme")
    print("  runs SATURATED and the reaction is zero order in urate: a dose")
    print("  buys a fixed mmol/h, not a fixed fractional reduction.  That")
    print("  capacity is therefore something you can outgrow.")
    print()
    row("dose mg/kg", "C_rasb mg/L", "capacity mmol/h", "burden it covers",
        w=[13, 14, 18, 20])
    for dose in [0.05, 0.10, 0.15, 0.20, 0.30, 0.40]:
        c = dose * P["BW"] / P["V1_RASB"]
        cap = P["VMAX_RASB"] * c
        # burden whose peak purine flux equals that capacity
        lo, hi = 0.0, 60.0
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if peak_release_flux(mid) < cap:
                lo = mid
            else:
                hi = mid
        row(f"{dose:.2f}", f"{c:.2f}", f"{cap:.2f}",
            f"{0.5 * (lo + hi):.1f} x1e12 cells", w=[13, 14, 18, 20])
    print("\n  Now the ODE, at the licensed 0.2 mg/kg, across burdens:")
    row("burden", "UA 4h", "UApk", "% at 4 h", "urate ox mmol", "CRrat",
        w=[10, 9, 9, 11, 15, 8])
    for N0 in [0.5, 1.0, 3.0, 6.0, 10.0]:
        r = run(burkitt(N0=N0, rasb=0.0))
        sm = summarise(r)
        i0 = int(np.argmin(np.abs(r["t"] - 0.0)))
        i4 = int(np.argmin(np.abs(r["t"] - 4.0)))
        row(f"{N0:.1f}", f"{r['UA_MGDL'][i4]:.2f}", f"{sm['UA_peak']:.1f}",
            f"{100 * (r['UA_MGDL'][i4] - r['UA_MGDL'][i0]) / r['UA_MGDL'][i0]:+.0f}%",
            f"{sm['H2O2']:.1f}", f"{sm['CR_ratio']:.2f}", w=[10, 9, 9, 11, 15, 8])
    print("\n  The 4-hour drop -- the endpoint the registration trials used --")
    print("  is large at every burden, because at 4 h the enzyme has had time")
    print("  to clear the PRE-EXISTING pool.  The peak is what the release")
    print("  flux does afterwards, and that is where the dose runs out.")


def A5c_secondorder():
    hdr("A5c. POTASSIUM AND PHOSPHATE ARE SECOND-ORDER SOLUTES")
    print("  A1 showed that at intact GFR only urate loses the race (ratio 6.2)")
    print("  while potassium (0.4) and phosphate (0.5) are excretable.  If that")
    print("  is right, their excursions are consequences of the urate limb's")
    print("  kidney injury, and a drug that touches ONLY urate should lower")
    print("  them.  Rasburicase has no potassium or phosphate pharmacology at")
    print("  all, so this is a falsifiable prediction of the structure.")
    print()
    row("arm", "UApk", "Kpk", "PO4pk", "Ca_ion", "eGFRn", "P(arr)",
        w=[26, 8, 8, 8, 8, 8, 9])
    for lab, spec in [("no prophylaxis, 2 L/day", burkitt(N0=3.0, hydration=HYD_STD)),
                      ("+ rasburicase only", burkitt(N0=3.0, hydration=HYD_STD,
                                                     rasb=0.0)),
                      ("+ febuxostat only", burkitt(N0=3.0, hydration=HYD_STD,
                                                    febu_start=-72.0)),
                      ("no ppx, 6e12 cells", burkitt(N0=6.0, hydration=HYD_STD)),
                      ("+ rasb 0.4 x5 only", burkitt(N0=6.0, hydration=HYD_STD,
                                                     rasb=0.0, rasb_days=5,
                                                     rasb_dose_mgkg=0.4))]:
        sm = summarise(run(spec))
        row(lab, f"{sm['UA_peak']:.1f}", f"{sm['K_peak']:.2f}",
            f"{sm['PO4_peak']:.2f}", f"{sm['CA_nadir']:.2f}",
            f"{sm['eGFR_nadir']:.0f}", f"{100 * sm['P_ARR']:.2f}%",
            w=[26, 8, 8, 8, 8, 8, 9])
    print("\n  The potassium and phosphate columns move in arms whose only")
    print("  intervention is a urate enzyme.  In this model hyperkalaemia in")
    print("  TLS is largely a renal-failure phenomenon rather than a release")
    print("  phenomenon, which is why it tracks the creatinine and not the LDH.")


def A6_residual():
    hdr("A6. WHAT SURVIVES DELETION OF THE URATE LIMB?")
    print("  The tubular injury drive is a SUM of four terms.  Each arm below")
    print("  prints the peak drive and its decomposition, so the claim that a")
    print("  residual exists -- and what species it belongs to -- is computed")
    print("  rather than asserted.")
    arms = [("hydration only", burkitt(N0=3.0)),
            ("hydration only, 2 L/day", burkitt(N0=3.0, hydration=HYD_STD)),
            ("+ rasburicase", burkitt(N0=3.0, rasb=0.0)),
            ("+ rasb + PO4 binder", burkitt(N0=3.0, rasb=0.0, seve=0.0)),
            ("+ rasb + binder + fluid", burkitt(N0=3.0, rasb=0.0, seve=0.0,
                                                hydration=HYD_MAX)),
            ("high burden, 2 L/day", burkitt(N0=6.0, hydration=HYD_STD)),
            ("high burden + rasb x1", burkitt(N0=6.0, rasb=0.0, hydration=HYD_STD)),
            ("high burden + rasb x5", burkitt(N0=6.0, rasb=0.0, rasb_days=5,
                                              hydration=HYD_STD)),
            ("high burden + rasb 0.4 x5", burkitt(N0=6.0, rasb=0.0, rasb_days=5,
                                                  rasb_dose_mgkg=0.4,
                                                  hydration=HYD_STD)),
            ("high burden + rasb 0.4 + max",
             burkitt(N0=6.0, rasb=0.0, rasb_days=5, rasb_dose_mgkg=0.4,
                     hydration=HYD_MAX)),
            ("... + early dialysis",
             burkitt(N0=6.0, rasb=0.0, rasb_days=5, rasb_dose_mgkg=0.4,
                     hydration=HYD_MAX, dialysis=(24.0, 120.0, 1.5)))]
    row("arm", "drive", "urateXT", "xanXT", "CaP", "sol.UA", "CRrat", "P(RRT)",
        w=[28, 8, 9, 8, 8, 8, 8, 9])
    res = []
    for lab, spec in arms:
        r = run(spec)
        s = summarise(r)
        i = int(np.argmax(r["INJ_DRIVE"]))
        tot = max(r["INJ_DRIVE"][i], 1e-12)
        res.append((lab, s, r, i))
        row(lab, f"{tot:.2f}", f"{100 * r['INJ_XUA'][i] / tot:.0f}%",
            f"{100 * r['INJ_XXAN'][i] / tot:.0f}%",
            f"{100 * r['INJ_XCAP'][i] / tot:.0f}%",
            f"{100 * r['INJ_URATE'][i] / tot:.0f}%",
            f"{s['CR_ratio']:.2f}", f"{100 * s['P_RRT']:.1f}%",
            w=[28, 8, 9, 8, 8, 8, 8, 9])
    print("\n  Percentages are the share of the PEAK injury drive carried by")
    print("  urate crystal, xanthine crystal, calcium phosphate, and soluble")
    print("  (non-crystal) urate toxicity.")
    d = {lab: (sm, r, i) for lab, sm, r, i in res}

    def share(lab, key="INJ_XCAP"):
        sm, r, i = d[lab]
        return 100.0 * r[key][i] / max(r["INJ_DRIVE"][i], 1e-12)

    print("\n  (1) A calcium-phosphate limb exists, and it is burden-gated:")
    print(f"      at 3e12 cells it carries {share('hydration only, 2 L/day'):.0f}%"
          f" of the peak injury drive;")
    print(f"      at 6e12 cells it carries {share('high burden, 2 L/day'):.0f}%,"
          f" and it is still {share('high burden + rasb 0.4 x5'):.0f}%")
    print("      after the urate limb has been suppressed as far as the enzyme")
    print("      can suppress it.  No urate-directed drug touches that share.")

    print("\n  (2) The phosphate binder cannot touch it either, and the reason")
    print("      is arithmetic rather than pharmacological:")
    j_lys = peak_release_flux(6.0, q=P["Q_PO4"])
    print(f"        phosphate entering ECF from lysis   {j_lys:6.2f} mmol/h")
    print(f"        phosphate entering ECF from diet    {P['PO4_INTAKE']:6.2f}"
          f" mmol/h")
    print(f"        share of the load a gut binder can reach"
          f" {100 * P['PO4_INTAKE'] / (j_lys + P['PO4_INTAKE']):5.1f}%")
    sm_a, sm_b = d["+ rasburicase"][0], d["+ rasb + PO4 binder"][0]
    print(f"      and indeed adding sevelamer moves peak Cr ratio "
          f"{sm_a['CR_ratio']:.2f} -> {sm_b['CR_ratio']:.2f}.")
    print("      Sequestration is an operator on the INTAKE term, and in acute")
    print("      TLS the intake term is not where the phosphate is coming from.")

    print("\n  (3) What does move it is removal from the body:")
    sm_c, sm_e = d["high burden + rasb 0.4 + max"][0], d["... + early dialysis"][0]
    print(f"        rasb 0.4 x5 + max fluid          Cr ratio "
          f"{sm_c['CR_ratio']:.2f}, P(RRT) {100 * sm_c['P_RRT']:.1f}%")
    print(f"        the same plus early dialysis     Cr ratio "
          f"{sm_e['CR_ratio']:.2f}, P(RRT) {100 * sm_e['P_RRT']:.1f}%")

    print("\n  (4) Dose escalation is what closes the urate limb, not repetition:")
    for lab in ["high burden + rasb x1", "high burden + rasb x5",
                "high burden + rasb 0.4 x5"]:
        sm, _, _ = d[lab]
        print(f"        {lab:30s} Cr ratio {sm['CR_ratio']:5.2f}"
              f"  P(RRT) {100 * sm['P_RRT']:5.1f}%")


def A7_burden_switch():
    hdr("A7. WHICH CRYSTAL RATE-LIMITS? burden sweep with rasburicase on board")
    row("burden", "Xua", "Xcap", "Xcap/Xua", "PO4pk", "CRrat", "P(RRT)",
        w=[10, 9, 9, 11, 8, 8, 9])
    for N0 in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
        s = summarise(run(burkitt(N0=N0, rasb=0.0)))
        row(f"{N0:.1f}", f"{s['XT_UA']:.3f}", f"{s['XT_CAP']:.3f}",
            f"{s['XT_CAP'] / max(s['XT_UA'], 1e-9):.1f}", f"{s['PO4_peak']:.2f}",
            f"{s['CR_ratio']:.2f}", f"{100 * s['P_RRT']:.1f}%",
            w=[10, 9, 9, 11, 8, 8, 9])


def A8_ramp():
    hdr("A8. FLUX-SHAPING -- why the venetoclax ramp exists")
    print("  Same drug, same target, same eventual kill.  Only the release RATE")
    print("  differs, because the clearance system is a low-pass filter on the")
    print("  release flux.")
    tmax = 24.0 * 42
    ramps = {
        "400 mg from day 1": [(0.0, 400.0)],
        "200 mg from day 1": [(0.0, 200.0)],
        "100 mg from day 1": [(0.0, 100.0)],
        "50 mg from day 1": [(0.0, 50.0)],
        "2-step 20/50": [(0.0, 20.0), (168.0, 50.0)],
        "5-week label ramp": [(0.0, 20.0), (168.0, 50.0), (336.0, 100.0),
                              (504.0, 200.0), (672.0, 400.0)],
        "8-week slow ramp": [(0.0, 10.0), (168.0, 20.0), (336.0, 50.0),
                             (504.0, 100.0), (672.0, 200.0), (840.0, 400.0)],
    }
    row("venetoclax schedule", "Jrel pk", "Kpk", "PO4pk", "UApk", "CRrat",
        "P(RRT)", "TLS", "lysed", "N end", w=[22, 9, 7, 7, 7, 7, 8, 5, 8, 8])
    for lab, steps in ramps.items():
        bol = []
        for i, (t0, dose) in enumerate(steps):
            t1 = steps[i + 1][0] if i + 1 < len(steps) else tmax
            bol += sched(t0, t1, 24.0, "VEN_G", dose)
        reg = Regimen(boluses=bol)
        # CLL, not Burkitt: doubling time ~3 months, so the ramp is not racing
        # tumour regrowth.  Using the Burkitt TD here would compare a fast ramp
        # against a tumour that doubles every 30 h and the answer would be an
        # artefact of the growth term rather than of the release waveform.
        T, Y, p = simulate(reg, tmax=tmax, N0=5.0, dt=0.5,
                           overrides=dict(FLUID=HYD_STD, TD=2000.0))
        r = tabulate(T, Y, p)
        s = summarise(r)
        jrel = float((P["Q_PUR"] * P["KLYS"] * r["N_LYS"]).max())
        row(lab, f"{jrel:.2f}", f"{s['K_peak']:.2f}", f"{s['PO4_peak']:.2f}",
            f"{s['UA_peak']:.1f}", f"{s['CR_ratio']:.2f}",
            f"{100 * s['P_RRT']:.1f}%",
            "C" if s["CTLS"] else ("L" if s["LTLS"] else "-"),
            f"{s['lysed']:.2f}", f"{r['N_VIA'][-1]:.4f}",
            w=[22, 9, 7, 7, 7, 7, 8, 5, 8, 8])
    print("\n  Jrel pk = peak purine release flux (mmol/h); lysed = cumulative")
    print("  tumour lysed (1e12 cells); N end = residual viable tumour at 6")
    print("  weeks.  Read `lysed` and `N end` together with `CRrat`: if the")
    print("  ramp lysed the same tumour with a lower peak flux and a lower")
    print("  creatinine, it bought safety without giving up kill.")
    print("\n  Note that the 2-step and the 5-week ramp are IDENTICAL on every")
    print("  safety column.  The peak release flux is set by the FIRST dose")
    print("  level, because by the time the later steps arrive the tumour that")
    print("  could have lysed is already gone.  That is the model's account of")
    print("  why the label puts the monitoring requirement on the 20 mg dose.")


def A9_potassium():
    hdr("A9. REDISTRIBUTION IS NOT CLEARANCE -- insulin vs dialysis vs binder")
    base = run(burkitt(N0=3.0, hydration=HYD_STD))
    i = int(np.argmax(base["K"]))
    y_at = np.array([base[s][i] for s in STATES])
    print(f"  Starting from the peak of the untreated trajectory: t = "
          f"{base['t'][i]:.0f} h, K = {base['K'][i]:.2f} mmol/L, "
          f"eGFR = {base['eGFR_mlmin'][i]:.0f} mL/min")
    print("  Total body exchangeable K (extracellular + parked) is printed so")
    print("  that the difference between moving and removing is visible.")
    arms = {
        "no rescue": {},
        "insulin 10 U + glucose": {"ins_inf": (0.0, 1.0, 900.0)},
        "insulin q4h x 6": {"rep_ins": True},
        "salbutamol 20 mg neb": {"salb": True},
        "insulin + salbutamol": {"ins_inf": (0.0, 1.0, 900.0), "salb": True},
        "SZC 10 g q8h": {"szc": True},
        "furosemide 40 mg q8h": {"furo": True},
        "haemodialysis 4 h": {"dialysis": (2.0, 6.0, 6.0)},
    }
    row("rescue", "K 2h", "K 6h", "K 12h", "K 24h", "TBK 24h", "dTBK", "P(arr)",
        w=[24, 8, 8, 8, 8, 9, 8, 9])
    for lab, kw in arms.items():
        kw = dict(kw)
        bol = []
        if kw.pop("rep_ins", False):
            bol += [(4.0 * k, "INS_C", 495.0) for k in range(6)]
        if kw.pop("salb", False):
            bol += sched(0.0, 24.0, 6.0, "SALB_C", 45.0)
        if kw.pop("szc", False):
            bol += sched(0.0, 48.0, 8.0, "SZC_G", 10.0)
        if kw.pop("furo", False):
            bol += sched(0.0, 48.0, 8.0, "FURO_C", 1.4)
        reg = Regimen(boluses=bol, **kw)
        T, Y, p = simulate(reg, tmax=48.0, y0=y_at, dt=0.1,
                           overrides=dict(FLUID=HYD_STD))
        rr = tabulate(T, Y, p)
        tbk = rr["K_EC"] + rr["K_SH"]

        def at(h, k="K"):
            return float(np.interp(h, rr["t"], rr[k]))
        row(lab, f"{at(2):.2f}", f"{at(6):.2f}", f"{at(12):.2f}", f"{at(24):.2f}",
            f"{np.interp(24.0, rr['t'], tbk):.0f}",
            f"{np.interp(24.0, rr['t'], tbk) - tbk[0]:+.0f}",
            f"{100 * (1 - math.exp(-rr['CH_ARR'][-1])):.1f}%",
            w=[24, 8, 8, 8, 8, 9, 8, 9])
    print("\n  dTBK = change in total body exchangeable K over 24 h.  Only the")
    print("  dialysis, binder and diuretic rows move it.  The shift rows lower")
    print("  the number on the chart and leave the patient's potassium where")
    print("  it was, which is the whole reason the effect wears off.")


def A10_calcium():
    hdr("A10. THE CALCIUM REFLEX -- correcting hypocalcaemia feeds the crystal")
    print("  Hypocalcaemia in TLS is a CONSEQUENCE of calcium-phosphate")
    print("  precipitation, so the calcium given to correct it is a reactant in")
    print("  the reaction that caused it.")
    row("Ca infusion", "Ca_ion", "CaxPO4", "Xcap", "CRrat", "eGFRn", "P(RRT)",
        "P(sz)", w=[16, 9, 9, 9, 8, 8, 9, 9])
    for ca in [0.0, 0.5, 1.5, 3.0, 6.0]:
        r = run(burkitt(N0=3.0, ca_inf=ca))
        s = summarise(r)
        row(f"{ca:.1f} mmol/h", f"{s['CA_nadir']:.2f}", f"{r['CAxPO4'].max():.0f}",
            f"{s['XT_CAP']:.2f}", f"{s['CR_ratio']:.2f}", f"{s['eGFR_nadir']:.0f}",
            f"{100 * s['P_RRT']:.1f}%", f"{100 * s['P_SZ']:.1f}%",
            w=[16, 9, 9, 9, 8, 8, 9, 9])
    print("\n  CaxPO4 in mg2/dL2 (classical deposition threshold 60-70).")
    print("  The seizure column is the reason the reflex exists; the RRT column")
    print("  is its price.")


def A11_bistable():
    hdr("A11. IS THE LOOP BISTABLE? insult against nephron reserve")
    print("  Identical tumour burden, different starting nephron mass (prior")
    print("  CKD).  A step rather than a slope means the loop gain crosses 1.")
    row("start NEPH", "eGFR 0", "CRrat", "eGFRn", "NEPH 14d", "Xtot", "P(RRT)",
        w=[12, 9, 8, 8, 10, 8, 9])
    for n in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]:
        reg, N0, tmax, ov = burkitt(N0=3.0)
        y0 = BASE.copy()
        y0[IX["NEPH"]] = n
        r = tabulate(*simulate(reg, tmax=tmax, N0=N0, y0=y0, overrides=ov))
        s = summarise(r)
        row(f"{n:.2f}", f"{r['eGFR_mlmin'][0]:.0f}", f"{s['CR_ratio']:.2f}",
            f"{s['eGFR_nadir']:.1f}", f"{r['NEPH'][-1]:.3f}",
            f"{s['XT_TOT']:.2f}", f"{100 * s['P_RRT']:.1f}%",
            w=[12, 9, 8, 8, 10, 8, 9])
    print("\n  Hysteresis test: same sub-threshold burden, but the kidney is")
    print("  SEEDED with crystal at t = 0.  If the loop were bistable, a large")
    print("  enough seed would push an otherwise-safe patient onto a second")
    print("  branch and keep them there.")
    row("seed Xua mmol", "CRrat", "eGFRn", "Xtot 14d", "eGFR 14d", "P(RRT)",
        w=[16, 8, 8, 10, 10, 9])
    for seed in [0.0, 2.0, 5.0, 10.0, 20.0, 40.0]:
        reg, N0, tmax, ov = burkitt(N0=2.0)
        y0 = BASE.copy()
        y0[IX["XT_UA"]] = seed
        r = tabulate(*simulate(reg, tmax=tmax, N0=N0, y0=y0, overrides=ov))
        sm = summarise(r)
        row(f"{seed:.1f}", f"{sm['CR_ratio']:.2f}", f"{sm['eGFR_nadir']:.1f}",
            f"{r['XT_TOT'][-1]:.2f}", f"{r['eGFR_mlmin'][-1]:.0f}",
            f"{100 * sm['P_RRT']:.1f}%", w=[16, 8, 8, 10, 10, 9])
    print("  If every row returns to the same eGFR at 14 days the system has")
    print("  ONE attractor and the steepness is a knee, not a switch.")
    print("\n  Same sweep on the insult axis at intact nephron mass:")
    row("burden", "CRrat", "eGFRn", "Xtot", "NEPH 14d", "P(RRT)",
        w=[12, 8, 8, 8, 10, 9])
    for N0 in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 6.0]:
        s = summarise(run(burkitt(N0=N0, hydration=HYD_STD)))
        row(f"{N0:.1f}", f"{s['CR_ratio']:.2f}", f"{s['eGFR_nadir']:.1f}",
            f"{s['XT_TOT']:.2f}", f"{s['NEPH_final']:.3f}",
            f"{100 * s['P_RRT']:.1f}%", w=[12, 8, 8, 8, 10, 9])


def A12_sampling():
    hdr("A12. CAIRO-BISHOP AS A DETECTION FILTER, NOT A DEFINITION")
    print("  The four analytes have different time constants, so '>= 2 of 4")
    print("  within a window' is partly a property of the blood-draw schedule.")
    scen = {
        "high burden, no ppx": burkitt(N0=3.0, hydration=HYD_STD),
        "high burden + rasb": burkitt(N0=3.0, rasb=0.0),
        "high burden + full bundle": burkitt(N0=3.0, rasb=0.0, seve=0.0,
                                             hydration=HYD_MAX),
        "moderate + allo -72 h": burkitt(N0=1.0, allo_start=-72.0),
        "low burden": burkitt(N0=0.2, allo_start=-72.0),
    }
    row("scenario", "q4h", "q6h", "q8h", "q12h", "q24h", "q48h", "n(q4h)",
        w=[28, 7, 7, 7, 7, 7, 7, 9])
    for lab, spec in scen.items():
        r = run(spec)
        cells = []
        for sh in [4.0, 6.0, 8.0, 12.0, 24.0, 48.0]:
            l, c, _ = cairo_bishop(r, sample_h=sh)
            cells.append("C" if c else ("L" if l else "-"))
        _, _, n4 = cairo_bishop(r, sample_h=4.0)
        row(lab, *cells, f"{n4}", w=[28, 7, 7, 7, 7, 7, 7, 9])
    print("\n  C = clinical TLS captured, L = laboratory TLS only, - = missed.")
    print("\n  WHICH of the four criteria fires, and what each one is actually")
    print("  detecting.  UA_crit from A1 is printed alongside the urate")
    print("  criterion, because in this model they are the same number.")
    row("scenario", "urate", "K", "PO4", "Ca", "n", "UApk", "UA_crit",
        w=[28, 8, 6, 7, 6, 4, 8, 9])
    qu = HYD_AGGR - P["INSENS"]
    for lab, spec in scen.items():
        r = run(spec)
        t = r["t"]
        i0 = int(np.argmin(np.abs(t - 0.0)))
        m = (t >= -72.0) & (t <= 168.0)
        c_ua = bool(np.any((r["UA_MGDL"][m] >= 8.0)
                           | (r["UA_MGDL"][m] >= 1.25 * r["UA_MGDL"][i0])))
        c_k = bool(np.any((r["K"][m] >= 6.0) | (r["K"][m] >= 1.25 * r["K"][i0])))
        c_p = bool(np.any((r["PO4"][m] >= 1.45)
                          | (r["PO4"][m] >= 1.25 * r["PO4"][i0])))
        c_c = bool(np.any((r["CA_T"][m] <= 1.75)
                          | (r["CA_T"][m] <= 0.75 * r["CA_T"][i0])))
        row(lab, "yes" if c_ua else "no", "yes" if c_k else "no",
            "yes" if c_p else "no", "yes" if c_c else "no",
            f"{int(c_ua) + int(c_k) + int(c_p) + int(c_c)}",
            f"{r['UA_MGDL'][m].max():.1f}",
            f"{ua_critical(qu) * UA_MGDL:.1f}", w=[28, 8, 6, 7, 6, 4, 8, 9])
    print("\n  The 2-of-4 rule turns out to be ROBUST to sampling, because")
    print("  three of the four analytes have multi-day time constants.  The")
    print("  number a clinician would act on is not: potassium is the fast one.")
    row("scenario", "true Kpk", "q4h", "q8h", "q12h", "q24h", "miss q24h",
        w=[28, 10, 8, 8, 8, 8, 11])
    for lab, spec in scen.items():
        r = run(spec)
        t = r["t"]
        m = (t >= 0.0) & (t <= 168.0)
        true_pk = float(r["K"][m].max())
        obs = []
        for sh in [4.0, 8.0, 12.0, 24.0]:
            want = np.arange(0.0, 168.0 + 1e-9, sh)
            idx = [int(np.argmin(np.abs(t - w))) for w in want]
            obs.append(float(r["K"][idx].max()))
        row(lab, f"{true_pk:.2f}", *[f"{o:.2f}" for o in obs],
            f"{true_pk - obs[-1]:+.2f}", w=[28, 10, 8, 8, 8, 8, 11])
    print("\n  miss q24h = how much of the true potassium peak a once-daily")
    print("  draw fails to see.")


def A13_oxidant():
    hdr("A13. THE OXIDANT LOAD SCALES WITH THE INDICATION")
    print("  Urate oxidase makes one H2O2 per urate destroyed, so the patient")
    print("  with the largest tumour receives the largest peroxide load -- the")
    print("  toxicity is proportional to the reason for giving the drug.")
    row("burden", "G6PD", "urate ox (mmol)", "MetHb %", "Hb nadir", "UApk",
        w=[10, 9, 17, 10, 11, 8])
    for N0 in [0.3, 1.0, 3.0, 6.0]:
        for g6 in [1.0, 0.0]:
            s = summarise(run(burkitt(N0=N0, rasb=0.0, overrides={"G6PD": g6})))
            row(f"{N0:.1f}", "normal" if g6 == 1.0 else "deficient",
                f"{s['H2O2']:.1f}", f"{s['MetHb']:.1f}", f"{s['Hb_nadir']:.1f}",
                f"{s['UA_peak']:.1f}", w=[10, 9, 17, 10, 11, 8])


def A14_oxypurinol():
    hdr("A14. THE LOOP BITES THE DRUG TOO -- oxypurinol in evolving AKI")
    print("  Oxypurinol is renally cleared; febuxostat is not.  The kidney")
    print("  injury the drug is given to prevent changes the drug's exposure.")
    row("arm", "C d1", "C d3", "C d7", "XOI d7", "eGFR d7",
        w=[26, 9, 9, 9, 9, 10])
    for lab, spec, key in [
            ("allopurinol, low burden",
             burkitt(N0=0.2, allo_start=-72.0), "C_OXY"),
            ("allopurinol, high burden",
             burkitt(N0=6.0, allo_start=-72.0, hydration=HYD_STD), "C_OXY"),
            ("febuxostat, high burden",
             burkitt(N0=6.0, febu_start=-72.0, hydration=HYD_STD), "C_FEBU")]:
        r = run(spec)

        def at(h, k):
            return float(np.interp(h, r["t"], r[k]))
        row(lab, f"{at(24, key):.2f}", f"{at(72, key):.2f}", f"{at(168, key):.2f}",
            f"{at(168, 'XOI'):.3f}", f"{at(168, 'eGFR_mlmin'):.0f}",
            w=[26, 9, 9, 9, 9, 10])
    print("\n  Concentrations in mg/L.  Oxypurinol accumulation is the model's")
    print("  account of why allopurinol needs renal dose reduction in exactly")
    print("  the patients whose TLS risk is highest.")


def A15_ledger():
    hdr("A15. TRIAL LEDGER -- model against published numbers")
    # Goldman 2001 gave rasburicase 0.20 mg/kg/day for 5-7 days, and its
    # population was mixed-risk paediatric NHL/ALL rather than the extreme
    # high-burden case this model is parameterised around, so the ledger uses
    # a 1e12-cell burden for the response-rate rows.
    r_rasb = run(burkitt(N0=1.0, rasb=0.0, rasb_days=5))
    r_allo = run(burkitt(N0=1.0, allo_start=0.0))
    out = []

    def pct4(r):
        i0 = int(np.argmin(np.abs(r["t"] - 0.0)))
        i4 = int(np.argmin(np.abs(r["t"] - 4.0)))
        return 100.0 * (r["UA_MGDL"][i4] - r["UA_MGDL"][i0]) / r["UA_MGDL"][i0]

    out.append(("rasburicase, urate change at 4 h", f"{pct4(r_rasb):+.0f}%",
                "-86% (Goldman 2001, Blood)"))
    out.append(("allopurinol, urate change at 4 h", f"{pct4(r_allo):+.0f}%",
                "+2% (Goldman 2001)"))

    def auc(r):
        m = (r["t"] >= 0.0) & (r["t"] <= 96.0)
        return float(np.trapezoid(r["UA_MGDL"][m], r["t"][m]))
    out.append(("urate AUC0-96 rasb vs allo",
                f"{auc(r_rasb) / auc(r_allo):.2f}x",
                "0.39x (Goldman 2001) -- model OVER-separates"))
    s_rasb, s_allo = summarise(r_rasb), summarise(r_allo)
    out.append(("urate held < 8 mg/dL, rasburicase",
                "yes" if s_rasb["UA_peak"] < 8.0 else "no",
                "87% response (Cortes 2010, JCO)"))
    out.append(("urate held < 8 mg/dL, allopurinol",
                "yes" if s_allo["UA_peak"] < 8.0 else "no",
                "66% response (Cortes 2010)"))
    s_bundle = summarise(run(burkitt(N0=3.0, rasb=0.0, seve=0.0,
                                     hydration=HYD_MAX)))
    out.append(("dialysis, high-risk with full bundle",
                f"{100 * s_bundle['P_RRT']:.0f}%",
                "1.5-5% overall (Coiffier 2008 JCO panel)"))
    s_none = summarise(run(burkitt(N0=3.0, hydration=HYD_STD)))
    out.append(("dialysis, high-risk unprophylaxed",
                f"{100 * s_none['P_RRT']:.0f}%",
                "up to ~30% in historical Burkitt series"))
    out.append(("LDH peak, high-burden Burkitt",
                f"{s_none['LDH_peak']:.0f} U/L",
                "commonly > 2 x ULN, often > 5000"))
    out.append(("K nadir shift, insulin 10 U",
                "see A9", "-0.6 to -1.0 mmol/L at 1-4 h (Allon 1990)"))
    out.append(("urine pH achieved, NaHCO3 loading",
                f"{summarise(run(burkitt(N0=3.0, hco3=45.0)))['pH']:.2f}",
                "target 6.5-7.5 in the alkalinisation era"))
    for a, b, c in out:
        print(f"  {a:38s} model {b:>9s}   reported {c}")


def A16_scenarios():
    hdr("A16. THE 12 SHIPPED SCENARIOS (same list as the R model)")
    scen = [
        ("S1  low-risk solid tumour", burkitt(N0=0.05, hydration=HYD_STD), None),
        ("S2  high burden, fluids only", burkitt(N0=3.0, hydration=HYD_STD), None),
        ("S3  + aggressive hydration", burkitt(N0=3.0), None),
        ("S4  + allopurinol t=0", burkitt(N0=3.0, allo_start=0.0), None),
        ("S5  + allopurinol -72 h", burkitt(N0=3.0, allo_start=-72.0), None),
        ("S6  + rasburicase 0.2 mg/kg", burkitt(N0=3.0, rasb=0.0), None),
        ("S7  + rasb, max fluid, no alkali",
         burkitt(N0=3.0, rasb=0.0, hydration=HYD_MAX), None),
        ("S8  + alkalinisation pH 7.5", burkitt(N0=3.0, hco3=45.0), None),
        ("S9  steroid prephase 5 d", burkitt(N0=3.0, pred_prephase=(-120.0, 5)), None),
        ("S10 full bundle", burkitt(N0=3.0, rasb=0.0, seve=0.0, allo_start=-72.0,
                                    hydration=HYD_MAX, furo=24.0), None),
        ("S11 CKD stage 3 host", burkitt(N0=3.0, rasb=0.0, allo_start=-72.0), 0.45),
        ("S12 G6PD deficient + rasb", burkitt(N0=3.0, rasb=0.0,
                                              overrides={"G6PD": 0.0}), None),
    ]
    row("scenario", "UApk", "Kpk", "PO4pk", "Ca_i", "CRrat", "eGFRn", "LDH",
        "lysed", "pre-t0", "P(RRT)", "TLS",
        w=[32, 6, 6, 7, 6, 6, 6, 7, 7, 8, 8, 5])
    for lab, spec, neph in scen:
        y0 = None
        if neph is not None:
            y0 = BASE.copy()
            y0[IX["NEPH"]] = neph
        r = run(spec, y0=y0)
        s = summarise(r)
        row(lab, f"{s['UA_peak']:.1f}", f"{s['K_peak']:.2f}", f"{s['PO4_peak']:.2f}",
            f"{s['CA_nadir']:.2f}", f"{s['CR_ratio']:.2f}", f"{s['eGFR_nadir']:.0f}",
            f"{s['LDH_peak']:.0f}", f"{s['lysed']:.2f}", f"{s['lysed_pre']:.2f}",
            f"{100 * s['P_RRT']:.1f}%",
            "C" if s["CTLS"] else ("L" if s["LTLS"] else "-"),
            w=[32, 6, 6, 7, 6, 6, 6, 7, 7, 8, 8, 5])
    print("\n  lysed = cumulative tumour lysed (1e12 cells); pre-t0 = how much")
    print("  of it was already lysed before cytotoxic therapy started, which is")
    print("  non-zero only for the prephase arm and is precisely the point of")
    print("  a prephase.")


def A17_massbalance():
    hdr("A17. NUMERICAL HYGIENE")
    r = run(burkitt(N0=3.0, rasb=0.0))
    print(f"  purine released by lysis           {P['Q_PUR'] * 3.0:9.2f} mmol")
    print(f"  endogenous purine input            "
          f"{P['P_END'] * (r['t'][-1] - r['t'][0]):9.2f} mmol")
    print(f"  terminal purine pools + crystals   "
          f"{r['HYPOX'][-1] + r['XAN'][-1] + r['URATE'][-1] + r['ALLANT'][-1] + r['XT_UA'][-1] + r['XT_XAN'][-1]:9.2f} mmol")
    print(f"  cumulative urate oxidised (= H2O2) {r['H2O2'][-1]:9.2f} mmol")
    print(f"  K released by lysis                {P['Q_K'] * 3.0:9.2f} mmol")
    print(f"  PO4 released by lysis              {P['Q_PO4'] * 3.0:9.2f} mmol")
    bad = [s for s in STATES if r[s].min() < -1e-6]
    print("  negative states: " + (", ".join(bad) if bad else "none"))
    print("  (the purine balance closes through renal excretion and HGPRT")
    print("   salvage, which are fluxes rather than accumulated states)")


def A18_timecourse():
    hdr("A18. TIME COURSE, HIGH-BURDEN BURKITT")
    for lab, spec in [("hydration only (2 L/day)", burkitt(N0=3.0, hydration=HYD_STD)),
                      ("full bundle", burkitt(N0=3.0, rasb=0.0, seve=0.0,
                                              allo_start=-72.0, hydration=HYD_MAX))]:
        print(f"\n  --- {lab} ---")
        r = run(spec)
        row("t (h)", "UA", "K", "PO4", "Ca_i", "Cr", "eGFR", "Xua", "Xcap",
            "LDH", "TUBI", w=[7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 7])
        for h in [0, 6, 12, 24, 36, 48, 72, 96, 120, 168, 240, 336]:
            i = int(np.argmin(np.abs(r["t"] - h)))
            row(f"{h}", f"{r['UA_MGDL'][i]:.1f}", f"{r['K'][i]:.2f}",
                f"{r['PO4'][i]:.2f}", f"{r['CA_ION'][i]:.2f}",
                f"{r['CR_MGDL'][i]:.2f}", f"{r['eGFR_mlmin'][i]:.0f}",
                f"{r['XT_UA'][i]:.2f}", f"{r['XT_CAP'][i]:.2f}",
                f"{r['LDH_UL'][i]:.0f}", f"{r['TUBINJ'][i]:.3f}",
                w=[7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 7])


ANALYSES = [A0_baseline, A1_race, A2_operators, A3_ph_optimum, A4_leadtime,
            A5_pool_vs_flux, A5b_capacity, A5c_secondorder, A6_residual, A7_burden_switch, A8_ramp,
            A9_potassium, A10_calcium, A11_bistable, A12_sampling, A13_oxidant,
            A14_oxypurinol, A15_ledger, A16_scenarios, A17_massbalance,
            A18_timecourse]
QUICK = [A0_baseline, A1_race, A2_operators, A18_timecourse]


def main():
    global BASE
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", type=str, default=None)
    a = ap.parse_args()
    print("Tumour lysis syndrome QSP model -- reference integration")
    print(f"states = {NS}, numpy {np.__version__}")
    BASE = equilibrate(P)
    todo = QUICK if a.quick else ANALYSES
    if a.only:
        want = {x.strip().upper() for x in a.only.split(",")}
        todo = [f for f in ANALYSES if f.__name__.split("_")[0].upper() in want]
        if A0_baseline not in todo:
            todo = [A0_baseline] + todo
    for f in todo:
        f()
    print("\ndone.")


if __name__ == "__main__":
    main()
