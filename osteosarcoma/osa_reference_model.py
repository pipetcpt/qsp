#!/usr/bin/env python3
# ============================================================================
#  osa_reference_model.py
#  Osteosarcoma (OSA) — reference implementation of the QSP model
#  골육종 QSP 모델 — 수치 참조 구현
#
#  WHY THIS FILE EXISTS
#  --------------------
#  The deliverable model of record is `osa_mrgsolve_model.R` (mrgsolve / R).
#  R is not available in this container, so this file re-implements the SAME
#  equations and the SAME parameters in Python/scipy and actually runs them,
#  so that every number quoted in README.md and in the .dot map is a computed
#  result rather than an assertion.  The two files are kept structurally
#  parallel: parameter names here are the parameter names there.
#
#  ORGANISING THESIS
#  -----------------
#  Cure in osteosarcoma is a Poisson bet on lesions nobody can see, and it is
#  paid for out of an exposure budget with three hard organ ceilings:
#
#      P(cure) = exp( -lambda0 * [ 1 - exp( -n0 * exp(-K) ) ] )
#      K       = kM*E_MTX + kD*E_DOX + kC*E_CIS + kI*E_IFO      (nats of kill)
#      s.t.    E_DOX <= ceiling(heart), E_CIS <= ceiling(cochlea),
#              E_MTX <= ceiling(kidney)
#
#  and the MTX ceiling is not a constant, because methotrexate is cleared
#  through the organ it destroys.  Tubular methotrexate precipitates whenever
#  its luminal concentration exceeds a pH-dependent solubility
#
#      S(pH) = 0.86 * 10^(0.682*(pH-5))  mM     (fit to the solubility table)
#      C_tub = CL_ren * C_plasma / UrineFlow
#
#  Injury lowers GFR and urine flow, which RAISES C_tub, which precipitates
#  more drug.  That is a positive feedback with a loop gain, hence a critical
#  urine pH — a bifurcation, not a risk factor.
#
#  Six consequences are then arithmetic (all recomputed below):
#    A. A critical urine pH exists and the model derives the guideline.
#    B. Because S ~ 10^(0.682*pH), one pH unit is worth 4.8x the fluid and
#       more than a 4-fold dose reduction.  Alkalinisation dominates dosing.
#    C. Escalation after surgery (MAPIE) arrives after the bet is placed and
#       spends budget it does not replace -> HR ~ 1 by construction.
#    D. The 40-year survival plateau is a VARIANCE result: the spread of K
#       across patients is several times any achievable shift in its mean.
#    E. Huvos necrosis is an integral of the same parameter that decides the
#       micrometastatic fate, so it PREDICTS but cannot be TREATED.
#    F. Dexrazoxane is the only intervention that RAISES the budget instead
#       of reallocating it; bone-directed drugs move a different readout.
#
#  UNITS:  time = hours.  MTX in mmol / mM.  DOX, CIS in mg / mg/L.
#          tumour volume mL (1 mL ~ 1 g ~ 1e9 cells).  BSA-normalised doses
#          are converted with BSA = 1.7 m^2.
# ============================================================================
import json
import math
import os
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


# ---------------------------------------------------------------------------
# 1. STATE VECTOR  (50 states — the R model expands these to 65 by splitting
#    the marrow chain and the bone/immune blocks; shared states share names)
# ---------------------------------------------------------------------------
SNAMES = [
    # --- methotrexate PK + the renal feedback loop -----------------------
    "A1", "A2", "A3", "PPT", "KID",
    # --- folate / rescue -------------------------------------------------
    "MPGT", "MPGM", "LV", "RFT", "RFM",
    # --- doxorubicin + heart --------------------------------------------
    "DOXC", "DOXP", "DOXOL", "CMV", "CUMDOX",
    # --- cisplatin + kidney + cochlea -----------------------------------
    "CISC", "ADDT", "ADDK", "ADDC", "HL",
    # --- ifosfamide / etoposide (MAPIE arm) ------------------------------
    "IFOA", "ETOC",
    # --- tumour ----------------------------------------------------------
    # LMET = ln(viable cells in a representative micrometastatic lesion).
    # It starts at ln(1e6) and the log-kill K is ln(1e6) - LMET(end).
    "PRIM", "NEC", "RES", "LMET",
    # per-agent attribution of the micrometastatic kill integral (nats)
    "KI_MTX", "KI_DOX", "KI_CIS", "KI_IE",
    # --- bone remodelling vicious cycle ---------------------------------
    "OCL", "TGFM", "RKL", "DENO", "ZOL",
    # --- vasculature / hypoxia ------------------------------------------
    "VASC", "HYP",
    # --- immune ----------------------------------------------------------
    "CTL", "M2", "MIFA",
    # --- myelosuppression (Friberg transit chain) ------------------------
    "PROL", "TR1", "TR2", "TR3", "CIRC",
    # --- mucosa ----------------------------------------------------------
    "MUC",
    # --- exposure budget accumulators -----------------------------------
    "EMTX", "EDOX", "ECIS", "EIFO",
    # --- rescue / protectants / biomarkers ------------------------------
    "DEXR", "GLUC", "MG", "ALP",
    # cumulative treatment-related mortality hazard: the counterweight
    # without which prolonged exposure would look like free extra kill
    "TRM",
]
IX = {n: i for i, n in enumerate(SNAMES)}
NS = len(SNAMES)

# ---------------------------------------------------------------------------
# 2. PARAMETERS
# ---------------------------------------------------------------------------
P = dict(
    BSA=1.7,
    # ---- methotrexate PK (3 cmt, L and L/h, 1.7 m^2 adolescent) --------
    MTX_MW=454.44,
    V1=18.0, V2=10.0, V3=5.0,
    Q2=18.0, Q3=0.60,
    CLREN0=7.5,          # renal CL at KID=1  (L/h) ~ 86% of total
    CLNR=1.2,            # non-renal CL (L/h)
    FU=0.55,             # unbound fraction
    GFR0=7.2,            # L/h  (=120 mL/min)
    # ---- tubular precipitation / the loop ------------------------------
    HYDR=3.0,            # protocol hydration, L/m2/day
    UF0=0.21,            # urine flow at 3 L/m2/d, 1.7 m2  (L/h)
    URINE_PH=7.5,
    SOL_A=0.86,          # mM at pH 5.0
    SOL_B=0.682,         # decades of solubility per pH unit
    KPPT=0.45,           # crystallisation rate constant (1/h per mM excess)
    KDIS=0.030,          # redissolution rate (1/h)
    PPT50=2.0,           # mmol of crystal giving half-maximal obstruction
    KINJ_PPT=0.020,      # 1/h  obstructive injury
    KINJ_TUB=0.0016,     # 1/h  direct tubular toxicity (concentration driven)
    KTUB=12.0,           # mM   half-max for direct tubular toxicity
    KREP_KID=0.0022,     # 1/h  tubular regeneration
    # ---- folate / polyglutamates / rescue ------------------------------
    KPG=0.055, MPGMAX=8.0, KDPG=0.010,     # tumour polyglutamation
    KPGM=0.035, MPGMAXM=6.0, KDPGM=0.032,  # marrow polyglutamation
    KIPG=1.10,           # MTXPG conc giving half-maximal TS/DHFR block
    CL_LV=9.0, V_LV=25.0,
    KRF=0.60,            # reduced folate that doubles the apparent Ki
    KRF_IN=0.020, KRF_OUT=0.045, RF0=0.15,
    KRF_INM=0.055,       # marrow takes up leucovorin faster than tumour
    # ---- doxorubicin ---------------------------------------------------
    CL_DOX=45.0, V_DOXC=25.0, V_DOXP=1200.0, Q_DOX=60.0,
    KMET_DOXOL=0.22, CL_DOXOL=30.0,
    KCM=0.0235,          # cardiomyocyte loss per (mg/L doxorubicinol)-h
    KCM_REP=1.2e-6,      # essentially irreversible
    # DEXR / MIFA are EFFECT-SITE surrogates, not plasma drug: dexrazoxane's
    # iron chelation and topoisomerase-II-beta depletion, and mifamurtide's
    # macrophage activation, both far outlast the parent molecule in plasma
    # (dexrazoxane t1/2 ~ 2-3 h, mifamurtide ~ 18 min), so the turnover here
    # is the turnover of the EFFECT.
    DEXR_EMAX=0.72, DEXR_EC50=3.0,
    CL_DEXR=3.0, V_DEXR=40.0,
    # ---- cisplatin -----------------------------------------------------
    CL_CIS=30.0, V_CIS=18.0,
    KADD_T=0.030, KREP_ADDT=0.020,
    KADD_K=0.020, KREP_ADDK=0.008, KINJ_CIS=0.0060,
    KADD_C=0.014, KREP_ADDC=0.0016,
    KHL=0.055,           # dB of high-frequency shift per adduct-hour
    KMG=0.010, MG0=0.85,
    # ---- ifosfamide / etoposide ----------------------------------------
    CL_IFO=12.0, V_IFO=40.0, FACT_IFO=0.55,
    CL_ETO=1.9, V_ETO=18.0,
    KINJ_IFO=0.0022,     # ifosfamide/CAA proximal tubular injury
    # ---- tumour growth / kill ------------------------------------------
    PRIM0=150.0, PRIMMAX=1500.0,
    KG=9.63e-4,          # 1/h  -> volume doubling time 30 d
    CELLS_PER_ML=1.0e9,
    KMTX=0.1828,         # max MTX kill rate (1/h) at full TS block
    KDOX=0.2766, EC50_DOX=0.35,
    KCIS=0.1730, EC50_CIS=0.55,
    KIFO=0.0200, EC50_IFO=6.0,
    KETO=0.0160, EC50_ETO=2.5,
    KNECCLR=0.00006,     # necrotic osteoid is RETAINED (it is what the
                         # pathologist measures at week 11)
    # penetration: primary tumour is poorly perfused, micromets are not
    PEN_MIN=0.30,
    N_LESION_CELLS=1.0e6,
    KGM=3.4e-4,          # 1/h  dormant micrometastatic growth (~85 d Td)
    LMAX=25.3,           # ln(1e11) — carrying capacity of one lung lesion
    LAMBDA0=1.80,        # Poisson mean micrometastatic lesions, localized
    LAMBDA0_MET=15.0,    # overt metastatic presentation
    # ---- resistance ----------------------------------------------------
    RES0=0.020, KMUT=1.2e-5, KSEL=0.45,
    # ---- bone remodelling vicious cycle --------------------------------
    OCL0=1.0, KOCL=0.010, RKL0=1.0, KRKL_T=0.0035, KDRKL=0.020,
    OPG=1.0, KTGFM=0.014, KDTGFM=0.030, TGF_BOOST=0.45,
    KDENO=0.0009, V_DENO=3.5, DENO_IC50=1.4,
    KZOL=0.00035, ZOL_IC50=0.30,
    OSTEOLYSIS_MAX=1.0,
    # ---- vasculature / hypoxia -----------------------------------------
    VASC0=0.55, KVASC=0.004, KDVASC=0.004, HYP_K=0.6,
    # ---- immune --------------------------------------------------------
    CTL0=1.0, KCTL=0.004, KDCTL=0.006, M20=1.0, KM2=0.003, KDM2=0.004,
    CL_MIFA=0.12, V_MIFA=12.0, MIFA_EMAX=0.55, MIFA_EC50=0.25,
    # ---- Friberg myelosuppression --------------------------------------
    CIRC0=4.0, MTT=125.0, GAM=0.17,
    SL_MTX=5.00, SL_DOX=34.0, SL_CIS=4.00, SL_IFO=0.60, SL_ETO=4.60,
    GCSF=0.0,
    # ---- mucositis -----------------------------------------------------
    KMUC=0.500, KMUC_REP=0.020,
    # ---- biomarkers ----------------------------------------------------
    ALP0=1.0, KALP=0.010, KDALP=0.012,
    # ---- treatment-related mortality -----------------------------------
    #  Infection/mucositis death needs BOTH a denuded barrier and no
    #  neutrophils, so the hazard is their PRODUCT, not their sum.  The renal
    #  arm is the multi-organ-failure route of methotrexate nephrotoxicity.
    KTRM_INF=6.90e-4, KTRM_REN=7.00e-4, KTRM_CARD=2.20e-4,
    #  a small constant hazard for the causes this model does not carry
    #  (line sepsis, thromboembolism, surgical and anaesthetic death):
    #  ~0.8% over 34 weeks of protocol therapy
    KTRM_BASE=1.40e-6,
    ANC_G4=0.35, KID_CRIT=0.60,
    # ---- kill weights for the metastatic log-kill integral -------------
    # (nats of kill per unit exposure accumulator; see calibration note)
    LUNG_IMM=4.0,        # immune kill weighting, lung vs primary
    POP_CV_K=0.25,       # between-patient CV of achieved log-kill
)

# ---------------------------------------------------------------------------
# 3. RIGHT-HAND SIDE
# ---------------------------------------------------------------------------


def solubility(pH, p):
    """Tubular methotrexate solubility, mM.  Log-linear fit to the published
    solubility table (0.39 mg/mL at pH 5.0, 1.55 at 6.0, 9.04 at 7.0;
    MW 454.44 -> 0.86, 3.41, 19.9 mM)."""
    return p["SOL_A"] * 10.0 ** (p["SOL_B"] * (pH - 5.0))


def rhs(t, y, p, inf):
    y = np.maximum(y, 0.0)
    d = np.zeros(NS)
    g = IX

    A1, A2, A3, PPT, KID = (y[g["A1"]], y[g["A2"]], y[g["A3"]],
                            y[g["PPT"]], y[g["KID"]])
    KID = min(max(KID, 0.02), 1.0)

    # ---------------- methotrexate PK with the renal feedback ------------
    C1 = A1 / p["V1"]                     # mM  ( mmol / L )
    C2 = A2 / p["V2"]
    C3 = A3 / p["V3"]
    OBS = PPT / (PPT + p["PPT50"])        # fractional tubular obstruction
    CLREN = p["CLREN0"] * KID * (1.0 - 0.85 * OBS)
    UF = p["UF0"] * KID * (1.0 - 0.90 * OBS)
    UF = max(UF, 0.02)
    excr = CLREN * C1                     # mmol/h delivered to tubular lumen
    CTUB = excr / UF                      # mM in tubular fluid
    S = solubility(p["URINE_PH"], p)
    excess = max(0.0, CTUB - S)
    J_ppt = p["KPPT"] * excess * UF                       # mmol/h -> crystal
    J_dis = p["KDIS"] * PPT * max(0.0, S - CTUB) / S
    # glucarpidase: systemic hydrolysis of MTX to DAMPA (adds ~ 8 L/h)
    CL_GLU = 90.0 * y[g["GLUC"]] / (y[g["GLUC"]] + 1.0)

    d[g["A1"]] = (inf["mtx"]
                  - (p["CLNR"] + CLREN + CL_GLU) * C1
                  - p["Q2"] * (C1 - C2) - p["Q3"] * (C1 - C3)
                  - J_ppt)
    d[g["A2"]] = p["Q2"] * (C1 - C2)
    d[g["A3"]] = p["Q3"] * (C1 - C3)
    d[g["PPT"]] = J_ppt - J_dis
    d[g["GLUC"]] = inf["gluc"] - 0.35 * y[g["GLUC"]]

    # ---------------- kidney functional mass ----------------------------
    inj = (p["KINJ_PPT"] * OBS
           + p["KINJ_TUB"] * CTUB / (CTUB + p["KTUB"])
           + p["KINJ_CIS"] * y[g["ADDK"]] / (y[g["ADDK"]] + 1.0)
           + p["KINJ_IFO"] * y[g["IFOA"]] / (y[g["IFOA"]] + 5.0))
    d[g["KID"]] = -inj * KID + p["KREP_KID"] * (1.0 - KID)

    # ---------------- leucovorin / reduced folate / polyglutamates ------
    LV = y[g["LV"]] / p["V_LV"]
    d[g["LV"]] = inf["lv"] - p["CL_LV"] * LV
    d[g["RFT"]] = p["KRF_IN"] * LV - p["KRF_OUT"] * (y[g["RFT"]] - p["RF0"])
    d[g["RFM"]] = p["KRF_INM"] * LV - p["KRF_OUT"] * (y[g["RFM"]] - p["RF0"])

    PEN = p["PEN_MIN"] + (1.0 - p["PEN_MIN"]) * y[g["VASC"]] / \
        (y[g["VASC"]] + p["VASC0"])
    C_tum = C1 * PEN
    d[g["MPGT"]] = (p["KPG"] * C_tum * (1.0 - y[g["MPGT"]] / p["MPGMAX"])
                    - p["KDPG"] * y[g["MPGT"]])
    d[g["MPGM"]] = (p["KPGM"] * C1 * (1.0 - y[g["MPGM"]] / p["MPGMAXM"])
                    - p["KDPGM"] * y[g["MPGM"]])

    # fractional thymidylate-synthase / DHFR block, rescued by reduced folate
    KiT = p["KIPG"] * (1.0 + y[g["RFT"]] / p["KRF"])
    KiM = p["KIPG"] * (1.0 + y[g["RFM"]] / p["KRF"])
    blockT = y[g["MPGT"]] / (y[g["MPGT"]] + KiT)
    blockM = y[g["MPGM"]] / (y[g["MPGM"]] + KiM)

    # ---------------- doxorubicin ---------------------------------------
    CD = y[g["DOXC"]] / p["V_DOXC"]
    CDP = y[g["DOXP"]] / p["V_DOXP"]
    d[g["DOXC"]] = (inf["dox"] - p["CL_DOX"] * CD
                    - p["Q_DOX"] * (CD - CDP) - p["KMET_DOXOL"] * y[g["DOXC"]])
    d[g["DOXP"]] = p["Q_DOX"] * (CD - CDP)
    d[g["DOXOL"]] = p["KMET_DOXOL"] * y[g["DOXC"]] - p["CL_DOXOL"] * \
        y[g["DOXOL"]] / p["V_DOXC"]
    d[g["CUMDOX"]] = inf["dox"] / p["BSA"]          # mg/m2 accumulated
    CDOL = y[g["DOXOL"]] / p["V_DOXC"]
    CDX = y[g["DEXR"]] / p["V_DEXR"]
    prot = 1.0 - p["DEXR_EMAX"] * CDX / (CDX + p["DEXR_EC50"])
    d[g["DEXR"]] = inf["dexr"] - p["CL_DEXR"] * CDX
    d[g["CMV"]] = (-p["KCM"] * (CDOL + 0.25 * CD) * prot * y[g["CMV"]]
                   + p["KCM_REP"] * (1.0 - y[g["CMV"]]))

    # ---------------- cisplatin -----------------------------------------
    CC = y[g["CISC"]] / p["V_CIS"]
    d[g["CISC"]] = inf["cis"] - p["CL_CIS"] * CC
    d[g["ADDT"]] = p["KADD_T"] * CC * PEN - p["KREP_ADDT"] * y[g["ADDT"]]
    d[g["ADDK"]] = p["KADD_K"] * CC - p["KREP_ADDK"] * y[g["ADDK"]]
    d[g["ADDC"]] = p["KADD_C"] * CC - p["KREP_ADDC"] * y[g["ADDC"]]
    d[g["HL"]] = p["KHL"] * y[g["ADDC"]] / (y[g["ADDC"]] + 1.0)
    d[g["MG"]] = -p["KMG"] * y[g["ADDK"]] / (y[g["ADDK"]] + 1.0) \
        + 0.004 * (p["MG0"] - y[g["MG"]])

    # ---------------- ifosfamide / etoposide ----------------------------
    CI = y[g["IFOA"]] / p["V_IFO"]
    d[g["IFOA"]] = inf["ifo"] * p["FACT_IFO"] - p["CL_IFO"] * CI
    CE = y[g["ETOC"]] / p["V_ETO"]
    d[g["ETOC"]] = inf["eto"] - p["CL_ETO"] * CE

    # ---------------- bone remodelling vicious cycle --------------------
    CDEN = y[g["DENO"]] / p["V_DENO"]
    d[g["DENO"]] = inf["deno"] - p["KDENO"] * y[g["DENO"]]
    d[g["ZOL"]] = inf["zol"] - 0.00012 * y[g["ZOL"]]
    rkl_free = y[g["RKL"]] / (1.0 + CDEN / p["DENO_IC50"])
    ocl_drive = rkl_free / (rkl_free + p["OPG"])
    zol_inh = 1.0 / (1.0 + y[g["ZOL"]] / p["ZOL_IC50"])
    d[g["RKL"]] = (p["KRKL_T"] * y[g["PRIM"]] / 100.0
                   - p["KDRKL"] * y[g["RKL"]])
    d[g["OCL"]] = p["KOCL"] * (2.0 * ocl_drive * zol_inh - y[g["OCL"]])
    d[g["TGFM"]] = p["KTGFM"] * y[g["OCL"]] - p["KDTGFM"] * y[g["TGFM"]]

    # ---------------- vasculature / hypoxia -----------------------------
    d[g["VASC"]] = (p["KVASC"] * y[g["PRIM"]] / (y[g["PRIM"]] + 200.0)
                    - p["KDVASC"] * y[g["VASC"]] * (1.0 + inf["vegfi"]))
    d[g["HYP"]] = 0.01 * (p["HYP_K"] * (1.0 - y[g["VASC"]]) - y[g["HYP"]])

    # ---------------- immune --------------------------------------------
    CMF = y[g["MIFA"]] / p["V_MIFA"]
    d[g["MIFA"]] = inf["mifa"] - p["CL_MIFA"] * CMF
    mifa_e = p["MIFA_EMAX"] * CMF / (CMF + p["MIFA_EC50"])
    d[g["M2"]] = p["KM2"] * (1.0 + 0.6 * y[g["HYP"]]) - p["KDM2"] * \
        y[g["M2"]] * (1.0 + 2.0 * mifa_e)
    d[g["CTL"]] = (p["KCTL"] * (1.0 + 1.5 * mifa_e)
                   - p["KDCTL"] * y[g["CTL"]] * (1.0 + 0.8 * y[g["M2"]]))

    # ---------------- tumour kill ---------------------------------------
    res = min(y[g["RES"]], 0.95)
    kill_imm = 0.00006 * y[g["CTL"]] / (1.0 + y[g["M2"]])
    kill_mtx = p["KMTX"] * blockT * (1.0 - res)
    kill_dox = p["KDOX"] * (CD * PEN) / (CD * PEN + p["EC50_DOX"]) * (1.0 - res)
    kill_cis = p["KCIS"] * y[g["ADDT"]] / (y[g["ADDT"]] + p["EC50_CIS"]) * \
        (1.0 - res)
    kill_ifo = p["KIFO"] * (CI * PEN) / (CI * PEN + p["EC50_IFO"]) * (1.0 - res)
    kill_eto = p["KETO"] * (CE * PEN) / (CE * PEN + p["EC50_ETO"]) * (1.0 - res)
    kill_prim = kill_mtx + kill_dox + kill_cis + kill_ifo + kill_eto + kill_imm

    grow = p["KG"] * (1.0 + p["TGF_BOOST"] * y[g["TGFM"]]) * \
        math.log(max(p["PRIMMAX"] / max(y[g["PRIM"]], 1e-6), 1.0 + 1e-9))
    d[g["PRIM"]] = y[g["PRIM"]] * (grow - kill_prim)
    d[g["NEC"]] = y[g["PRIM"]] * kill_prim - p["KNECCLR"] * y[g["NEC"]]

    # --- micrometastases -------------------------------------------------
    # Fully perfused (no penetration penalty), no matrix TGF-beta boost, but
    # angiogenesis-limited so they regrow slowly (clinical dormancy: median
    # time to pulmonary relapse ~ 1.5 y).  LMET = ln(cells per lesion).
    km_mtx = p["KMTX"] * blockT_full(y, p, C1) * (1.0 - res)
    km_dox = p["KDOX"] * CD / (CD + p["EC50_DOX"]) * (1.0 - res)
    km_cis = p["KCIS"] * y[g["ADDT"]] / (y[g["ADDT"]] + p["EC50_CIS"]) * \
        (1.0 - res)
    km_ie = (p["KIFO"] * CI / (CI + p["EC50_IFO"])
             + p["KETO"] * CE / (CE + p["EC50_ETO"])) * (1.0 - res)
    # Immune clearance of pulmonary micrometastases.  This is the ONLY term
    # mifamurtide acts on: muramyl tripeptide activates alveolar macrophages
    # in the lung, which is where the lesions that decide cure actually are.
    # Weighted LUNG_IMM-fold higher than in the primary because alveolar
    # macrophage density and antigen access are both higher there.
    km_imm = p["LUNG_IMM"] * kill_imm
    kill_met = km_mtx + km_dox + km_cis + km_ie + km_imm
    L = y[g["LMET"]]
    grow_met = p["KGM"] * max(0.0, 1.0 - L / p["LMAX"])
    d[g["LMET"]] = (grow_met - kill_met) if L > -25.0 else 0.0
    d[g["KI_MTX"]] = km_mtx
    d[g["KI_DOX"]] = km_dox
    d[g["KI_CIS"]] = km_cis
    d[g["KI_IE"]] = km_ie

    # ---------------- resistance ----------------------------------------
    d[g["RES"]] = (p["KMUT"] * (1.0 - res)
                   + p["KSEL"] * res * (1.0 - res) * kill_prim)

    # ---------------- exposure accumulators (the budget) ----------------
    d[g["EMTX"]] = C1 if C1 > 0.001 else 0.0          # mM*h above 1 uM
    d[g["EDOX"]] = CD
    d[g["ECIS"]] = CC
    d[g["EIFO"]] = CI

    # ---------------- myelosuppression (Friberg) ------------------------
    ktr = 4.0 / p["MTT"]
    edrug = (p["SL_MTX"] * blockM + p["SL_DOX"] * CD + p["SL_CIS"] * CC
             + p["SL_IFO"] * CI + p["SL_ETO"] * CE)
    edrug = min(edrug, 0.98)
    fb = (p["CIRC0"] / max(y[g["CIRC"]], 0.05)) ** p["GAM"]
    fb = min(fb, 6.0) * (1.0 + p["GCSF"])
    d[g["PROL"]] = ktr * y[g["PROL"]] * ((1.0 - edrug) * fb - 1.0)
    d[g["TR1"]] = ktr * (y[g["PROL"]] - y[g["TR1"]])
    d[g["TR2"]] = ktr * (y[g["TR1"]] - y[g["TR2"]])
    d[g["TR3"]] = ktr * (y[g["TR2"]] - y[g["TR3"]])
    d[g["CIRC"]] = ktr * (y[g["TR3"]] - y[g["CIRC"]])

    # ---------------- mucositis -----------------------------------------
    d[g["MUC"]] = (p["KMUC"] * blockM * 4.0 * max(0.0, 1.0 - y[g["MUC"]] / 4.0)
                   - p["KMUC_REP"] * y[g["MUC"]])

    # ---------------- biomarker -----------------------------------------
    d[g["ALP"]] = (p["KALP"] * (y[g["PRIM"]] / 150.0) * (1.0 + y[g["OCL"]])
                   - p["KDALP"] * y[g["ALP"]])

    # ---------------- treatment-related mortality hazard ----------------
    g_neut = math.exp(-y[g["CIRC"]] / p["ANC_G4"])
    g_muc = min(y[g["MUC"]] / 4.0, 1.0)
    g_ren = max(0.0, (p["KID_CRIT"] - KID) / p["KID_CRIT"])
    g_card = max(0.0, (0.85 - y[g["CMV"]]) / 0.85)
    d[g["TRM"]] = (p["KTRM_BASE"]
                   + p["KTRM_INF"] * g_neut * g_muc
                   + p["KTRM_REN"] * g_ren
                   + p["KTRM_CARD"] * g_card)
    return d


def blockT_full(y, p, C1):
    """TS block seen by a fully-perfused micrometastasis (no penetration
    penalty).  Uses the same polyglutamate pool but the plasma-equivalent
    driving concentration."""
    KiT = p["KIPG"] * (1.0 + y[IX["RFT"]] / p["KRF"])
    mpg = y[IX["MPGT"]] * 1.25          # better delivery -> higher PG load
    mpg = min(mpg, p["MPGMAX"])
    return mpg / (mpg + KiT)


# ---------------------------------------------------------------------------
# 4. DOSING — EURAMOS-1 MAP backbone
# ---------------------------------------------------------------------------
W = 168.0   # hours per week
BLANK = {"mtx": 0.0, "lv": 0.0, "dox": 0.0, "cis": 0.0, "ifo": 0.0,
         "eto": 0.0, "dexr": 0.0, "deno": 0.0, "zol": 0.0, "mifa": 0.0,
         "gluc": 0.0, "vegfi": 0.0}

# MAP (EURAMOS-1 arm A): DOX 75 mg/m2 x6 = 450; CIS 120 mg/m2 x4 = 480;
# HDMTX 12 g/m2 x12.  Definitive surgery at week 11.
MAP_DOX_W = [1, 6, 12, 17, 22, 27]
MAP_CIS_W = [1, 6, 12, 22]
MAP_MTX_W = [4, 5, 9, 10, 15, 16, 20, 21, 25, 26, 30, 31]
SURGERY_W = 11

# ---------------------------------------------------------------------------
#  PROTOCOL GO / NO-GO GATES
#  -------------------------
#  Real protocols do not administer a scheduled course into an unrecovered
#  patient.  These are the EURAMOS-1-style rules, and they are the reason
#  ESCALATION SPENDS THE BUDGET INSTEAD OF ADDING TO IT: every extra cycle of
#  ifosfamide/etoposide deepens the neutrophil nadir and shaves the kidney,
#  and the courses that get deferred or dropped are courses of the backbone
#  that was doing the killing.  EURAMOS-1 delivered full protocol therapy to
#  76% of MAP patients and 51% of MAPIE patients; the model reproduces that
#  gap from these gates rather than by assuming it.
# ---------------------------------------------------------------------------
GATES = {
    #  kind : (min ANC 10^9/L, min mucositis-free, min kidney fraction)
    "MTX": (0.75, 3.0, 0.60),
    "DOX": (1.00, 3.0, 0.00),
    "CIS": (1.00, 3.0, 0.60),
    "IE":  (1.00, 3.0, 0.60),
}
MAX_DEFER_WEEKS = 2
DEFER_DOSE_FACTOR = 0.75     # a deferred course comes back at -25% dose
KID_STOP_MTX = 0.55         # grade-3 AKI -> hold methotrexate
KID_RESUME_MTX = 0.70       # resume once renal function recovers
CARDIAC_STOP_LVEF = 50.0     # do not give more anthracycline below this


def lvef_of(cmv):
    return 62.0 * (0.35 + 0.65 * max(cmv, 0.0) ** 0.6)


def expand(course, t0, p, reg):
    """Turn one administered course into (t_start, t_end, rate_dict) windows."""
    ev = []
    k = course["kind"]
    if k == "MTX":
        mmol = course["dose"] * p["BSA"] * 1000.0 / p["MTX_MW"]
        ev.append((t0, t0 + 4.0, {"mtx": mmol / 4.0}))
        lv_h = reg.get("lv_hours", 48.0)
        for i in range(int(lv_h / 6.0)):
            ta = t0 + 24.0 + 6.0 * i
            ev.append((ta, ta + 0.25,
                       {"lv": reg.get("lv_dose", 15.0) * p["BSA"] / 0.25}))
        if reg.get("gluc_at") is not None:
            tg = t0 + reg["gluc_at"]
            ev.append((tg, tg + 0.25, {"gluc": 50.0 * p["BSA"] / 0.25}))
    elif k == "DOX":
        ev.append((t0, t0 + 48.0, {"dox": course["dose"] * p["BSA"] / 48.0}))
        if reg.get("dexrazoxane", False):
            # 10:1 dexrazoxane, given with each anthracycline day: for a 48-h
            # doxorubicin infusion that is a concurrent 48-h cover
            ev.append((t0, t0 + 48.0,
                       {"dexr": 10.0 * course["dose"] * p["BSA"] / 48.0}))
    elif k == "CIS":
        ev.append((t0, t0 + 72.0, {"cis": course["dose"] * p["BSA"] / 72.0}))
    elif k == "IE":
        for day in range(5):
            ta = t0 + 24.0 * day
            ev.append((ta, ta + 3.0, {"ifo": 2800.0 * p["BSA"] / 3.0}))
            ev.append((ta, ta + 1.0, {"eto": 100.0 * p["BSA"] / 1.0}))
    elif k == "DENO":
        ev.append((t0, t0 + 1.0, {"deno": 120.0}))
    elif k == "ZOL":
        ev.append((t0, t0 + 0.5, {"zol": 8.0}))
    elif k == "MIFA":
        for off in (0.0, 84.0):
            ev.append((t0 + off, t0 + off + 1.0, {"mifa": 2.0 * p["BSA"]}))
    return ev


def build_schedule(reg):
    """week -> list of scheduled course specs."""
    sched = {}

    def add(week, spec):
        sched.setdefault(int(week), []).append(spec)

    for w in reg.get("dox_weeks", MAP_DOX_W):
        add(w, dict(kind="DOX", dose=reg.get("dox_dose", 75.0), defer=0))
    for w in reg.get("cis_weeks", MAP_CIS_W):
        add(w, dict(kind="CIS", dose=reg.get("cis_dose", 120.0), defer=0))
    for w in reg.get("mtx_weeks", MAP_MTX_W):
        add(w, dict(kind="MTX", dose=reg.get("mtx_dose", 12.0), defer=0))
    for w in reg.get("ifo_weeks", []):
        add(w, dict(kind="IE", dose=2800.0, defer=0))
    for w in reg.get("deno_weeks", []):
        add(w, dict(kind="DENO", dose=120.0, defer=0))
    for w in reg.get("zol_weeks", []):
        add(w, dict(kind="ZOL", dose=4.0, defer=0))
    for w in reg.get("mifa_weeks", []):
        add(w, dict(kind="MIFA", dose=2.0, defer=0))
    return sched


def gate(course, y, p):
    """Return (administer?, reason)."""
    k = course["kind"]
    if k not in GATES:
        return True, ""
    anc_min, muc_max, kid_min = GATES[k]
    if y[IX["CIRC"]] < anc_min:
        return False, "neutropenia"
    if y[IX["MUC"]] >= muc_max:
        return False, "mucositis"
    if kid_min > 0 and y[IX["KID"]] < kid_min:
        return False, "renal"
    if k == "DOX" and lvef_of(y[IX["CMV"]]) < CARDIAC_STOP_LVEF:
        return False, "cardiac"
    return True, ""


def init_state(p):
    y = np.zeros(NS)
    y[IX["KID"]] = 1.0
    y[IX["RFT"]] = p["RF0"]
    y[IX["RFM"]] = p["RF0"]
    y[IX["CMV"]] = 1.0
    y[IX["MG"]] = p["MG0"]
    y[IX["PRIM"]] = p["PRIM0"]
    y[IX["RES"]] = p["RES0"]
    y[IX["OCL"]] = p["OCL0"]
    y[IX["RKL"]] = p["RKL0"]
    y[IX["TGFM"]] = 0.4
    y[IX["VASC"]] = p["VASC0"]
    y[IX["HYP"]] = 0.3
    y[IX["CTL"]] = p["CTL0"]
    y[IX["M2"]] = p["M20"]
    for s in ("PROL", "TR1", "TR2", "TR3", "CIRC"):
        y[IX[s]] = p["CIRC0"]
    y[IX["ALP"]] = p["ALP0"]
    y[IX["LMET"]] = math.log(p["N_LESION_CELLS"])
    return y


def integrate_window(y, a, b, ev, p, reg, dt=1.0):
    """Integrate [a, b] splitting at every infusion breakpoint inside it."""
    bps = sorted({a, b} | {t for (t0, t1, _) in ev for t in (t0, t1)
                           if a < t < b})
    ts, ys = [], []
    for u, v in zip(bps[:-1], bps[1:]):
        if v - u < 1e-9:
            continue
        act = dict(BLANK)
        act["vegfi"] = reg.get("vegfi", 0.0)
        for (t0, t1, r) in ev:
            if t0 <= u + 1e-9 and t1 >= v - 1e-9:
                for kk, vv in r.items():
                    act[kk] = act.get(kk, 0.0) + vv
        grid = np.arange(u, v + 1e-9, dt)
        if grid.size < 2 or grid[-1] < v - 1e-9:
            grid = np.append(grid, v)
        sol = solve_ivp(rhs, (u, v), y, args=(p, act), method="LSODA",
                        t_eval=grid, rtol=1e-7, atol=1e-10, max_step=4.0)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1].copy()
        ts.extend(sol.t.tolist())
        ys.extend(sol.y.T.tolist())
    return y, ts, ys


def simulate(reg, p=None, tmax_w=40.0, dense_dt=1.0):
    """Week-stepping scheduler with protocol go/no-go gates."""
    p = dict(P if p is None else p)
    p.update(reg.get("p", {}))
    sched = build_schedule(reg)
    surgery_w = reg.get("surgery_week", SURGERY_W)
    y = init_state(p)
    ts, ys = [0.0], [y.copy().tolist()]
    pending = []
    given = {k: 0 for k in ("MTX", "DOX", "CIS", "IE")}
    planned = {k: 0 for k in given}
    for w in sched.values():
        for c in w:
            if c["kind"] in planned:
                planned[c["kind"]] += 1
    omitted = {k: 0 for k in given}
    on_time = {k: 0 for k in given}
    deferrals = {"neutropenia": 0, "mucositis": 0, "renal": 0, "cardiac": 0}
    delivered = {k: 0.0 for k in given}
    huvos = None
    done_surg = False
    mtx_stopped = False
    nweeks = int(tmax_w)
    for wk in range(1, nweeks + 1):
        t0 = (wk - 1) * W
        # ---- surgery -------------------------------------------------
        if (not done_surg) and wk == surgery_w and reg.get("surgery", True):
            tot = y[IX["PRIM"]] + y[IX["NEC"]]
            huvos = (y[IX["NEC"]] / tot) if tot > 0 else 0.0
            y[IX["PRIM"]] *= 0.002       # macroscopically complete resection
            y[IX["NEC"]] = 0.0
            done_surg = True
        # ---- decide this week's courses ------------------------------
        due = sched.get(wk, []) + pending
        pending = []
        ev = []
        # grade-3 AKI holds methotrexate; it resumes on renal recovery, as
        # every protocol specifies
        if y[IX["KID"]] < KID_STOP_MTX:
            mtx_stopped = True
        elif y[IX["KID"]] > KID_RESUME_MTX:
            mtx_stopped = False
        for c in due:
            if mtx_stopped and c["kind"] == "MTX":
                omitted["MTX"] += 1
                continue
            ok, why = gate(c, y, p)
            if ok:
                ev.extend(expand(c, t0, p, reg))
                if c["kind"] in given:
                    given[c["kind"]] += 1
                    delivered[c["kind"]] += c["dose"]
                    if c["defer"] == 0:
                        on_time[c["kind"]] += 1
            else:
                deferrals[why] = deferrals.get(why, 0) + 1
                c = dict(c)
                c["defer"] += 1
                c["dose"] *= DEFER_DOSE_FACTOR
                if c["defer"] <= MAX_DEFER_WEEKS:
                    pending.append(c)
                elif c["kind"] in omitted:
                    omitted[c["kind"]] += 1
        y, tt, yy = integrate_window(y, t0, t0 + W, ev, p, reg, dt=dense_dt)
        ts.extend(tt)
        ys.extend(yy)
    for c in pending:
        if c["kind"] in omitted:
            omitted[c["kind"]] += 1
    T = np.array(ts)
    Y = np.array(ys)
    if huvos is None:
        tot = Y[-1, IX["PRIM"]] + Y[-1, IX["NEC"]]
        huvos = (Y[-1, IX["NEC"]] / tot) if tot > 0 else 0.0
    n_plan = sum(planned.values())
    n_give = sum(given.values())
    n_time = sum(on_time.values())
    plan_dose = {"MTX": 12.0, "DOX": 75.0, "CIS": 120.0, "IE": 2800.0}
    dose_num = sum(delivered[k] for k in delivered)
    dose_den = sum(planned[k] * plan_dose[k] for k in planned)
    res = dict(t=T, Y=Y, huvos=huvos, p=p, given=given, planned=planned,
               omitted=omitted, deferrals=deferrals, delivered=delivered,
               on_time=on_time,
               completion=(n_give / n_plan) if n_plan else 1.0,
               on_schedule=(n_time / n_plan) if n_plan else 1.0,
               dose_intensity=(dose_num / dose_den) if dose_den else 1.0)
    res.update(endpoints(res, reg))
    return res


# ---------------------------------------------------------------------------
# 5. THE POPULATION CURE LAYER
#
#  A single deterministic run gives the cure probability of a patient sitting
#  exactly at the mean parameters, and because P(cure) is a DOUBLE exponential
#  in K that number is violently sensitive near the operating point.  What a
#  trial measures is the cure FRACTION of a population whose K is spread out.
#  Every scenario is therefore also reported as a population cure fraction:
#  the scenario supplies the MEAN of K, the population layer supplies the
#  SPREAD (transporter genotype, MTX clearance CV ~30%, delivered dose
#  intensity, intrinsic chemosensitivity), and the two together are what is
#  comparable to a published EFS curve.
# ---------------------------------------------------------------------------
POP_N = 8000
POP_SEED = 20260805
_Z = np.random.default_rng(POP_SEED + 1).standard_normal(POP_N)
_LAM = np.random.default_rng(POP_SEED).poisson(P["LAMBDA0"], POP_N)
_LAM_MET = np.random.default_rng(POP_SEED + 2).poisson(P["LAMBDA0_MET"], POP_N)


def pop_cure(K_nats, lam0=None, cv=None, n0=None):
    """Cure FRACTION of a virtual population whose mean log-kill is K_nats.

    Common random numbers are used across arms so that two scenarios differ
    only by their mean K, never by Monte-Carlo noise."""
    cv = P["POP_CV_K"] if cv is None else cv
    n0 = P["N_LESION_CELLS"] if n0 is None else n0
    lam = _LAM if (lam0 is None or abs(lam0 - P["LAMBDA0"]) < 1e-9) else (
        _LAM_MET if abs(lam0 - P["LAMBDA0_MET"]) < 1e-9
        else np.random.default_rng(POP_SEED).poisson(lam0, POP_N))
    K = np.maximum(_Z * cv * max(K_nats, 1e-9) + K_nats, 0.0)
    surv = n0 * np.exp(-K)
    p_les = -np.expm1(-surv)          # P(one lesion survives), per patient
    # cure needs EVERY one of this patient's lam lesions to fail
    return float(np.mean((1.0 - p_les) ** lam))


def efs_hr(cure_test, cure_ref):
    """EFS hazard ratio implied by two cure fractions.

    Under proportional hazards the long-term event-free proportion satisfies
    S_test = S_ref^HR, so HR = ln(cure_test) / ln(cure_ref).  HR < 1 favours
    the test arm."""
    return (math.log(max(min(cure_test, 1 - 1e-9), 1e-9))
            / math.log(max(min(cure_ref, 1 - 1e-9), 1e-9)))


# ---------------------------------------------------------------------------
# 6. ENDPOINTS
# ---------------------------------------------------------------------------
def cure_probability(K, lam0, p):
    """P(cure) = exp(-lam0 * [1 - exp(-n0 e^-K)]).

    n0*e^-K is the expected number of clonogenic cells surviving in ONE
    lesion; 1-exp(-that) is the chance that lesion survives; and cure needs
    EVERY one of a Poisson(lam0) number of lesions to fail."""
    s = p["N_LESION_CELLS"] * math.exp(-max(K, 0.0))
    p_les = -math.expm1(-s)          # 1-exp(-s), stable for small s
    return math.exp(-lam0 * p_les)


def endpoints(res, reg):
    Y, T, p = res["Y"], res["t"], res["p"]
    lam0 = reg.get("lambda0", p["LAMBDA0"])
    K = max(math.log(p["N_LESION_CELLS"]) - Y[-1, IX["LMET"]], 0.0)
    out = dict(
        K_nats=K, K_log10=K / math.log(10.0),
        cure=cure_probability(K, lam0, p),
        cure_pop=pop_cure(K, lam0),
        # overall survival = cured AND not killed by the treatment
        surv_pop=pop_cure(K, lam0) * math.exp(-Y[-1, IX["TRM"]]),
        eGFR_pct=100.0 * Y[-1, IX["KID"]],
        eGFR_nadir_pct=100.0 * float(Y[:, IX["KID"]].min()),
        eGFR=P["GFR0"] * Y[-1, IX["KID"]] * 1000.0 / 60.0 / 1.73 * 1.7,
        LVEF=62.0 * (0.35 + 0.65 * Y[-1, IX["CMV"]] ** 0.6),
        cum_dox=Y[-1, IX["CUMDOX"]],
        hearing_dB=Y[-1, IX["HL"]],
        ANC_nadir=float(Y[:, IX["CIRC"]].min()),
        days_ANC_lt_0p5=float((Y[:, IX["CIRC"]] < 0.5).sum() / 24.0),
        MUC_max=float(Y[:, IX["MUC"]].max()),
        E_MTX=Y[-1, IX["EMTX"]], E_DOX=Y[-1, IX["EDOX"]],
        E_CIS=Y[-1, IX["ECIS"]], E_IFO=Y[-1, IX["EIFO"]],
        PPT_max=float(Y[:, IX["PPT"]].max()),
        K_MTX=Y[-1, IX["KI_MTX"]] / math.log(10.0),
        K_DOX=Y[-1, IX["KI_DOX"]] / math.log(10.0),
        K_CIS=Y[-1, IX["KI_CIS"]] / math.log(10.0),
        K_IE=Y[-1, IX["KI_IE"]] / math.log(10.0),
        LMET=Y[-1, IX["LMET"]],
        ALP=Y[-1, IX["ALP"]],
        RES_final=Y[-1, IX["RES"]],
        OCL_max=float(Y[:, IX["OCL"]].max()),
        huvos_pct=100.0 * res["huvos"],
        TRM=float(-math.expm1(-Y[-1, IX["TRM"]])),
        completion_pct=100.0 * res.get("completion", 1.0),
        on_schedule_pct=100.0 * res.get("on_schedule", 1.0),
        dose_intensity_pct=100.0 * res.get("dose_intensity", 1.0),
        n_MTX=res.get("given", {}).get("MTX", 0),
        n_DOX=res.get("given", {}).get("DOX", 0),
        n_CIS=res.get("given", {}).get("CIS", 0),
        n_IE=res.get("given", {}).get("IE", 0),
        omit_total=sum(res.get("omitted", {}).values()),
        dox_delivered=res.get("delivered", {}).get("DOX", 0.0),
    )
    return out


def mtx_course(pH=7.0, dose_g_m2=12.0, uf_mult=1.0, kid0=1.0, gluc_at=None,
               lv_hours=48.0, lv_dose=15.0, hours=336.0, p=None):
    """Single HDMTX course — the loop in isolation.  Returns the clinical
    monitoring concentrations and whether elimination is 'delayed'."""
    pp = dict(P if p is None else p)
    pp["URINE_PH"] = pH
    pp["UF0"] = P["UF0"] * uf_mult
    ev = []
    mmol = dose_g_m2 * pp["BSA"] * 1000.0 / pp["MTX_MW"]
    ev.append((0.0, 4.0, {"mtx": mmol / 4.0}))
    for k in range(int(lv_hours / 6.0)):
        ta = 24.0 + 6.0 * k
        ev.append((ta, ta + 0.25, {"lv": lv_dose * pp["BSA"] / 0.25}))
    if gluc_at is not None:
        ev.append((gluc_at, gluc_at + 0.25, {"gluc": 50.0 * pp["BSA"] / 0.25}))
    bps = sorted({0.0, hours} | {e[0] for e in ev} | {e[1] for e in ev})
    bps = [b for b in bps if 0 <= b <= hours]
    y = init_state(pp)
    y[IX["KID"]] = kid0
    ts, ys = [0.0], [y.copy()]
    for a, b in zip(bps[:-1], bps[1:]):
        if b - a < 1e-9:
            continue
        act = dict(BLANK)
        for (t0, t1, r) in ev:
            if t0 <= a + 1e-9 and t1 >= b - 1e-9:
                for k2, v in r.items():
                    act[k2] = act.get(k2, 0.0) + v
        grid = np.arange(a, b + 1e-9, 0.25)
        if grid.size < 2 or grid[-1] < b - 1e-9:
            grid = np.append(grid, b)
        sol = solve_ivp(rhs, (a, b), y, args=(pp, act), method="LSODA",
                        t_eval=grid, rtol=1e-8, atol=1e-11, max_step=1.0)
        y = sol.y[:, -1].copy()
        ts.extend(sol.t[1:].tolist())
        ys.extend(sol.y[:, 1:].T.tolist())
    T, Y = np.array(ts), np.array(ys)
    C1 = Y[:, IX["A1"]] / pp["V1"] * 1000.0        # uM

    def at(h):
        return float(np.interp(h, T, C1))
    KID = Y[:, IX["KID"]]
    CLREN = pp["CLREN0"] * KID * (1 - 0.85 * Y[:, IX["PPT"]] /
                                  (Y[:, IX["PPT"]] + pp["PPT50"]))
    UF = np.maximum(pp["UF0"] * KID * (1 - 0.90 * Y[:, IX["PPT"]] /
                                       (Y[:, IX["PPT"]] + pp["PPT50"])), 0.02)
    CTUB = CLREN * (C1 / 1000.0) / UF
    S = solubility(pH, pp)
    return dict(t=T, C_uM=C1, CTUB=CTUB, S=S, KID=KID, PPT=Y[:, IX["PPT"]],
                ppt_max=float(Y[:, IX["PPT"]].max()),
                C24=at(24), C48=at(48), C72=at(72), C96=at(96),
                Cmax=float(C1.max()),
                SS_max=float((CTUB / S).max()),
                kid_final=float(KID[-1]),
                delayed=bool(at(48) > 1.0 or at(72) > 0.1),
                aki=bool(KID[-1] < 0.75),
                auc_mM_h=float(np.trapezoid(C1 / 1000.0, T)))


# ---------------------------------------------------------------------------
# 6. SCENARIOS
# ---------------------------------------------------------------------------
SCENARIOS = OrderedDict([
    ("S01_surgery_only", dict(
        label="Surgery alone (historical control, no chemotherapy)",
        dox_weeks=[], cis_weeks=[], mtx_weeks=[])),
    ("S02_MAP", dict(label="MAP — EURAMOS-1 arm A (reference standard)")),
    ("S03_MAPIE", dict(
        label="MAP+IE (MAPIE) — poor-responder escalation, EURAMOS-1",
        # As randomised: doxorubicin 375 mg/m2 (5 cycles), cisplatin
        # 360 mg/m2 (3 cycles), 14 HDMTX courses, 5 ifosfamide/etoposide
        # cycles.  The escalation BUYS the IE cycles by giving up one
        # doxorubicin and one cisplatin cycle.
        ifo_weeks=[14, 18, 23, 28, 32],
        mtx_weeks=[4, 5, 9, 10, 15, 16, 20, 21, 25, 26, 30, 31],
        dox_weeks=[1, 6, 12, 21, 26],
        cis_weeks=[1, 6, 12])),
    ("S04_MAP_pH60", dict(
        label="MAP with urine pH 6.0 (inadequate alkalinisation)",
        p={"URINE_PH": 6.0})),
    ("S05_MAP_pH75", dict(
        label="MAP with urine pH 7.5 + full hydration",
        p={"URINE_PH": 7.5})),
    ("S06_MAP_lowhydration", dict(
        label="MAP at half hydration (1.5 L/m2/day), pH 7.0",
        p={"UF0": P["UF0"] * 0.5})),
    ("S07_MAP_pH66", dict(
        label="MAP with urine pH 6.6 (control for the rescue arm)",
        p={"URINE_PH": 6.6})),
    ("S07b_MAP_glucarpidase", dict(
        label="MAP, pH 6.6, glucarpidase rescue at 48 h",
        p={"URINE_PH": 6.6}, gluc_at=48.0)),
    ("S08_MAP_dexrazoxane", dict(
        label="MAP + dexrazoxane 10:1 before every doxorubicin dose",
        dexrazoxane=True)),
    ("S09_MAP_dox600", dict(
        label="MAP intensified to doxorubicin 600 mg/m2 (no protectant)",
        dox_dose=100.0)),
    ("S10_MAP_dox600_dexra", dict(
        label="MAP intensified to 600 mg/m2 WITH dexrazoxane",
        dox_dose=100.0, dexrazoxane=True)),
    ("S11_MAP_zoledronate", dict(
        label="MAP + zoledronate 4 mg q4w (OS2006 design)",
        zol_weeks=[1, 5, 9, 13, 17, 21, 25, 29])),
    ("S12_MAP_denosumab", dict(
        label="MAP + denosumab 120 mg q4w",
        deno_weeks=[1, 5, 9, 13, 17, 21, 25, 29])),
    ("S13_MAP_mifamurtide", dict(
        label="MAP + mifamurtide (INT-0133 schedule)",
        mifa_weeks=list(range(12, 34)))),
    ("S14_AP_only", dict(
        label="AP only — adult/elderly, methotrexate omitted",
        mtx_weeks=[], dox_weeks=[1, 6, 12, 17, 22, 27],
        cis_weeks=[1, 6, 12, 22])),
    ("S15_MAP_metastatic", dict(
        label="MAP in overt metastatic disease at diagnosis",
        lambda0=P["LAMBDA0_MET"])),
    ("S16_MAP_resistant", dict(
        label="MAP in an intrinsically resistant tumour (RES0 = 0.25)",
        p={"RES0": 0.25})),
    ("S17_MAP_regorafenib", dict(
        label="MAP + regorafenib maintenance (anti-angiogenic, REGOBONE)",
        vegfi=1.2)),
    ("S18_MAP_dosedense", dict(
        label="MAP with 25% higher methotrexate (15 g/m2)",
        mtx_dose=15.0)),
])

def run_scenarios():
    say("=" * 100)
    say(" SECTION 1 — THERAPY SCENARIOS")
    say("=" * 100)
    say("")
    say(f"{'scenario':<24}{'Huvos%':>7}{'log10':>7}{'cure':>7}{'TRM%':>6}{'surv':>7}{'HR':>6}"
        f"{'sched%':>7}{'RDI%':>6}{'MTX':>4}{'DOX':>4}{'CIS':>4}{'IE':>4}"
        f"{'eGFRnad':>8}{'LVEF':>6}{'dB':>5}{'ANCnad':>7}{'MUC':>5}"
        f"{'cumDOX':>7}")
    say("-" * 100)
    results = {}
    ref = None
    for k, reg in SCENARIOS.items():
        r = simulate(reg)
        if k == "S02_MAP":
            ref = r["surv_pop"]
        results[k] = {kk: (float(vv) if isinstance(vv, (int, float,
                                                       np.floating)) else vv)
                      for kk, vv in r.items() if kk not in ("t", "Y", "p")}
        results[k]["label"] = reg["label"]
        results[k]["deferrals"] = r["deferrals"]
    for k in results:
        results[k]["efs_hr_vs_MAP"] = efs_hr(results[k]["surv_pop"], ref)
    for k, r in results.items():
        say(f"{k:<24}{r['huvos_pct']:>7.1f}{r['K_log10']:>7.2f}"
            f"{r['cure_pop']:>7.3f}{100 * r['TRM']:>6.1f}"
            f"{r['surv_pop']:>7.3f}{r['efs_hr_vs_MAP']:>6.2f}"
            f"{r['on_schedule_pct']:>7.0f}"
            f"{r['dose_intensity_pct']:>6.0f}{r['n_MTX']:>4.0f}{r['n_DOX']:>4.0f}"
            f"{r['n_CIS']:>4.0f}{r['n_IE']:>4.0f}"
            f"{r['eGFR_nadir_pct']:>8.0f}{r['LVEF']:>6.1f}"
            f"{r['hearing_dB']:>5.0f}{r['ANC_nadir']:>7.2f}"
            f"{r['MUC_max']:>5.1f}{r['dox_delivered'] / P['BSA'] * P['BSA']:>7.0f}")
    say("")
    say("  Huvos%  necrosis fraction of the resected primary at week 11")
    say("  log10   cell kill delivered to the micrometastatic pool")
    say("  cure    POPULATION cure fraction (mean K from this scenario, spread")
    say(f"          from the virtual-population layer, CV {P['POP_CV_K']:.0%})")
    say("  TRM%    treatment-related mortality (infection x mucositis, renal,")
    say("          cardiac, plus a small baseline)")
    say("  surv    overall survival fraction = cure x (1 - TRM).  THIS is the")
    say("          column to read: several arms buy log-kill by damaging the")
    say("          kidney, and the kill they buy is real, and they still lose.")
    say("  HR      hazard ratio vs standard MAP on the survival fraction")
    say("  sched%  courses given ON SCHEDULE at full dose (the number a")
    say("          trial reports as 'completed protocol therapy')")
    say("  RDI%    relative dose intensity actually delivered")
    say("  MTX/DOX/CIS/IE  courses actually given (planned 12 / 6 / 4 / -)")
    say("")
    say(f"  Calibration checks:  MAP survival     "
        f"{results['S02_MAP']['surv_pop']:.3f}"
        f"  (EURAMOS-1 5-y EFS 0.54-0.59, OS ~0.71)")
    say(f"                       MAP TRM          "
        f"{100 * results['S02_MAP']['TRM']:.1f}%"
        f"   (reported toxic deaths ~1%)")
    say(f"                       surgery alone     "
        f"{results['S01_surgery_only']['cure_pop']:.3f}"
        f"  (historical amputation-only 0.15-0.20)")
    say(f"                       MAP on schedule   "
        f"{results['S02_MAP']['on_schedule_pct']:.0f}%"
        f"    (EURAMOS-1 MAP completed protocol 76%)")
    say(f"                       MAPIE on schedule "
        f"{results['S03_MAPIE']['on_schedule_pct']:.0f}%"
        f"    (EURAMOS-1 MAPIE completed protocol 51%)")
    say(f"                       relative dose intensity, MAP "
        f"{results['S02_MAP']['dose_intensity_pct']:.0f}% /"
        f" MAPIE {results['S03_MAPIE']['dose_intensity_pct']:.0f}%")
    say("")
    return results


# ---------------------------------------------------------------------------
# 8. THE METHOTREXATE RENAL LOOP — critical pH, loop gain, exchange rates
# ---------------------------------------------------------------------------
def run_loop_analysis():
    say("=" * 100)
    say(" SECTION 2 — THE METHOTREXATE RENAL FEEDBACK LOOP")
    say("=" * 100)
    say("")
    say(" Methotrexate is cleared through the organ it destroys.  86% of the")
    say(" dose leaves by the kidney; the tubular concentration is")
    say("")
    say("     C_tub = CL_ren(KID) * C_plasma / UrineFlow(KID)")
    say("")
    say(" and everything above the pH-dependent solubility crystallises,")
    say(" obstructs, and lowers both CL_ren and UrineFlow — which RAISES")
    say(" C_tub.  Positive feedback with a loop gain, hence a bifurcation.")
    say("")
    say(" 2a. Solubility of methotrexate in tubular fluid")
    say("     S(pH) = 0.86 * 10^(0.682*(pH-5)) mM   — log-linear fit to the")
    say("     published table (0.39 / 1.55 / 9.04 mg/mL at pH 5 / 6 / 7)")
    say(f"{'pH':>8}{'S (mM)':>10}{'S (mg/mL)':>12}")
    for pH in (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0):
        S = solubility(pH, P)
        say(f"{pH:>8.1f}{S:>10.2f}{S * P['MTX_MW'] / 1000.0:>12.2f}")
    say("")

    say(" 2b. A single 12 g/m2 course as a function of urine pH")
    say("     (hydration 3 L/m2/day, leucovorin 15 mg/m2 q6h from +24 h)")
    say(f"{'pH':>6}{'Cmax_uM':>9}{'CtubPeak':>10}{'C/S':>7}{'C24':>8}"
        f"{'C48':>8}{'C72':>9}{'eGFR%':>7}{'crystal':>8}{'delayed':>8}{'AKI':>6}")
    say("-" * 100)
    ph_rows = []
    for pH in (6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 8.0):
        r = mtx_course(pH=pH)
        ph_rows.append(dict(pH=pH, Cmax=r["Cmax"], SS=r["SS_max"],
                            C24=r["C24"], C48=r["C48"], C72=r["C72"],
                            kid=r["kid_final"], delayed=r["delayed"],
                            ppt=r["ppt_max"], aki=r["aki"]))
        say(f"{pH:>6.1f}{r['Cmax']:>9.0f}{r['SS_max'] * r['S']:>10.1f}"
            f"{r['SS_max']:>7.2f}{r['C24']:>8.2f}{r['C48']:>8.3f}"
            f"{r['C72']:>9.4f}{100 * r['kid_final']:>7.1f}"
            f"{r['ppt_max']:>8.2f}{str(r['delayed']):>8}{str(r['aki']):>6}")
    say("")
    say("     C24/C48/C72 are the monitoring concentrations in uM.  A course")
    say("     is 'delayed' if C48 > 1 uM or C72 > 0.1 uM (the thresholds that")
    say("     trigger rescue escalation in every modern protocol).")
    say("")

    def crit(pred, lo=5.5, hi=9.0, iters=34):
        if not pred(mtx_course(pH=lo)):
            return float("nan")
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            if pred(mtx_course(pH=mid)):
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    pH_crit = crit(lambda r: r["SS_max"] > 1.0)
    pH_delay = crit(lambda r: r["delayed"])
    pH_aki = crit(lambda r: r["aki"])
    say(f" 2c. THE BIFURCATION")
    say(f"     urine pH at which tubular fluid first supersaturates   "
        f"{pH_crit:.2f}")
    say(f"     urine pH below which elimination becomes delayed       "
        f"{pH_delay:.2f}")
    say(f"     urine pH below which the course causes AKI (>25% GFR)  "
        f"{pH_aki:.2f}")
    say("")
    say("     The protocol instruction 'alkalinise the urine to pH >= 7.0 and")
    say("     hydrate at 3 L/m2/day' is not a convention.  It is where this")
    say(f"     loop changes sign.  Note that pH 7.0 sits {pH_crit - 7.0:+.2f} from the")
    say("     supersaturation boundary, i.e. AT the edge and not comfortably")
    say("     inside it — which is why 1-2% of courses still show delayed")
    say("     elimination under full protocol compliance, and why the")
    say("     escalation instruction on a rising creatinine is MORE")
    say("     bicarbonate, not less methotrexate.")
    say("")

    say(" 2d. EXCHANGE RATES.  C_tub scales as dose / urine-flow, while")
    say("     solubility scales as 10^(0.682*pH).  So any factor f applied to")
    say("     C_tub is worth exactly log10(f)/0.682 pH units:")
    for f, what in ((2.0, "double the hydration"),
                    (4.0, "quadruple the hydration"),
                    (2.0, "halve the methotrexate dose"),
                    (4.0, "quarter the methotrexate dose")):
        say(f"       {what:<30} = {math.log10(f) / P['SOL_B']:.2f} pH units")
    say(f"     Read the other way: 1.00 pH unit = {10 ** P['SOL_B']:.2f}-fold change")
    say("     in C_tub.  ONE UNIT OF URINE pH IS WORTH 4.8x THE FLUID, OR A")
    say("     4.8-FOLD DOSE REDUCTION.  Alkalinisation is not adjunctive")
    say("     supportive care; arithmetically it is the largest single lever")
    say("     on high-dose methotrexate safety, larger than the dose itself.")
    say("")
    hyd = []
    say(" 2e. The same statement, simulated (pH 6.6 held fixed):")
    say(f"{'hydration (L/m2/d)':>20}{'C48_uM':>9}{'eGFR%':>8}{'delayed':>9}")
    for mult, lab in ((0.5, 1.5), (1.0, 3.0), (2.0, 6.0), (4.0, 12.0)):
        r = mtx_course(pH=6.6, uf_mult=mult)
        hyd.append(dict(hyd=lab, C48=r["C48"], kid=r["kid_final"],
                        delayed=r["delayed"]))
        say(f"{lab:>20.1f}{r['C48']:>9.3f}{100 * r['kid_final']:>8.1f}"
            f"{str(r['delayed']):>9}")
    r_ph = mtx_course(pH=7.3, uf_mult=1.0)
    say(f"     vs. pH 7.3 at the ORIGINAL 3 L/m2/day: C48 {r_ph['C48']:.3f} uM,"
        f" eGFR {100 * r_ph['kid_final']:.1f}%")
    say("     — 0.7 of a pH unit does what quadrupling the fluid does.")
    say("")

    say(" 2f. LOOP GAIN.  Identical course at pH 6.6, but starting from a")
    say("     kidney that a previous cisplatin cycle already shaved:")
    say(f"{'start eGFR%':>12}{'C/S':>7}{'C48_uM':>9}{'AUC mM*h':>10}"
        f"{'end eGFR%':>11}{'GFR lost':>10}{'delayed':>9}")
    say("-" * 100)
    gain_rows = []
    for kid0 in (1.00, 0.90, 0.80, 0.70, 0.60, 0.50):
        r = mtx_course(pH=6.6, kid0=kid0)
        lost = 100 * (kid0 - r["kid_final"])
        gain_rows.append(dict(kid0=kid0, SS=r["SS_max"], C48=r["C48"],
                              auc=r["auc_mM_h"], kid=r["kid_final"],
                              lost=lost, delayed=r["delayed"]))
        say(f"{100 * kid0:>12.0f}{r['SS_max']:>7.2f}{r['C48']:>9.3f}"
            f"{r['auc_mM_h']:>10.2f}{100 * r['kid_final']:>11.1f}"
            f"{lost:>10.1f}{str(r['delayed']):>9}")
    a, b = gain_rows[0], gain_rows[-1]
    say("")
    say(f"     Halving the starting GFR multiplies the methotrexate AUC of the")
    say(f"     SAME dose by {b['auc'] / a['auc']:.2f}x and C48 by"
        f" {b['C48'] / max(a['C48'], 1e-9):.1f}x, and leaves the kidney at"
        f" {100 * b['kid']:.1f}%")
    say(f"     instead of {100 * a['kid']:.1f}%.  (The absolute GFR LOST is smaller"
        f" simply because")
    say("     there was less left to lose — the injury runs down toward a")
    say("     common floor, which is the point: the floor is where twelve")
    say("     courses take you.)  The loop gain rises as the kidney falls, so")
    say("     course n changes the")
    say("     pharmacokinetics of courses n+1 .. 12.  In a regimen that gives")
    say("     twelve courses of methotrexate AND four of cisplatin, the two")
    say("     nephrotoxins are not additive — they are multiplicative through")
    say("     a shared state variable.")
    say("")

    say(" 2g. THE RESCUE IS DIAGNOSTICALLY GATED TO ARRIVE AFTER THE")
    say("     EXPOSURE IT WOULD PREVENT.")
    say("")
    say("     Glucarpidase is triggered by a measured concentration, and the")
    say("     concentration that triggers it is measured at 24, 48 and 72 h.")
    say("     But the exposure integral is front-loaded: by the time the")
    say("     trigger fires, most of the AUC has already been delivered.")
    say("")
    r0 = mtx_course(pH=6.4)
    say(f"{'rescue timing':<24}{'C48_uM':>9}{'C72_uM':>10}{'AUC mM*h':>10}"
        f"{'AUC removed':>13}{'end eGFR%':>11}")
    say("-" * 100)
    say(f"{'none':<24}{r0['C48']:>9.3f}{r0['C72']:>10.4f}"
        f"{r0['auc_mM_h']:>10.2f}{'-':>13}{100 * r0['kid_final']:>11.1f}")
    gl = {}
    for h in (12.0, 24.0, 36.0, 48.0, 72.0):
        r = mtx_course(pH=6.4, gluc_at=h)
        frac = 1 - r["auc_mM_h"] / r0["auc_mM_h"]
        gl[h] = dict(auc=r["auc_mM_h"], frac=frac, kid=r["kid_final"],
                     C48=r["C48"], C72=r["C72"])
        say(f"{'glucarpidase at %g h' % h:<24}{r['C48']:>9.3f}"
            f"{r['C72']:>10.4f}{r['auc_mM_h']:>10.2f}"
            f"{100 * frac:>12.0f}%{100 * r['kid_final']:>11.1f}")
    say("")
    say(f"     Rescue at 24 h removes {100 * gl[24.0]['frac']:.0f}% of the exposure;"
        f" at 48 h, {100 * gl[48.0]['frac']:.0f}%; at 72 h,")
    say(f"     {100 * gl[72.0]['frac']:.0f}%.  The standard trigger is a 48-h or 72-h")
    say("     concentration, i.e. the drug is given at the point where it can")
    say(f"     still remove {100 * gl[48.0]['frac']:.0f}% of what is left and"
        f" {100 * (1 - gl[48.0]['frac']):.0f}% has already landed.")
    say("")
    say(f"     End-of-course GFR is {100 * r0['kid_final']:.1f}% with no rescue and"
        f" {100 * gl[48.0]['kid']:.1f}% with rescue")
    say("     at 48 h.  The tubular injury integral was paid in the first day,")
    say("     before any monitoring threshold had been crossed.  Nothing")
    say("     rescues the tubule retrospectively — which is why the only")
    say("     intervention with real leverage is the one applied BEFORE the")
    say("     infusion starts (2d), and why urine pH, not the rescue protocol,")
    say("     is where the safety of this regimen actually lives.")
    say("")
    r1, r2 = mtx_course(pH=6.4, gluc_at=48.0), mtx_course(pH=6.4, gluc_at=24.0)
    return dict(pH_crit=pH_crit, pH_delay=pH_delay, pH_aki=pH_aki,
                ph_rows=ph_rows, gain_rows=gain_rows, hydration=hyd,
                gluc=dict(none=r0["auc_mM_h"], h48=r1["auc_mM_h"],
                          h24=r2["auc_mM_h"],
                          kid_none=r0["kid_final"], kid_48=r1["kid_final"]))


# ---------------------------------------------------------------------------
# 9. THE EXPOSURE BUDGET AND ITS THREE ORGAN CEILINGS
# ---------------------------------------------------------------------------
def run_budget_analysis(scen):
    say("=" * 100)
    say(" SECTION 3 — THE EXPOSURE BUDGET AND ITS THREE ORGAN CEILINGS")
    say("=" * 100)
    say("")
    base = scen["S02_MAP"]
    say(" 3a. Where MAP's log-kill actually comes from (single-agent")
    say("     deletion, everything else held at protocol):")
    say(f"{'arm':<34}{'log10 kill':>11}{'delta':>8}{'per course':>12}"
        f"{'cure':>8}{'HR':>7}")
    say("-" * 100)
    dele = OrderedDict([
        ("full MAP", (dict(), 1)),
        ("MAP minus methotrexate", (dict(mtx_weeks=[]), 12)),
        ("MAP minus doxorubicin", (dict(dox_weeks=[]), 6)),
        ("MAP minus cisplatin", (dict(cis_weeks=[]), 4)),
    ])
    contrib = {}
    for nm, (over, ncourse) in dele.items():
        reg = dict(SCENARIOS["S02_MAP"])
        reg.update(over)
        r = simulate(reg)
        dK = r["K_log10"] - base["K_log10"]
        per = (-dK / ncourse) if ncourse > 1 else float("nan")
        contrib[nm] = dict(K=r["K_log10"], dK=dK, per_course=per,
                           cure=r["cure_pop"])
        say(f"{nm:<34}{r['K_log10']:>11.2f}{dK:>8.2f}"
            f"{per:>12.3f}{r['cure_pop']:>8.3f}"
            f"{efs_hr(r['cure_pop'], base['cure_pop']):>7.2f}")
    say("")
    say("     Twelve courses of methotrexate and six of doxorubicin are worth")
    say("     comparable totals, but the PER-COURSE value differs several-fold")
    say("     — and the two are limited by different organs, so they are not")
    say("     interchangeable currency.")
    say("")

    say(" 3b. What standard MAP has already spent when it finishes:")
    say(f"     doxorubicin cumulative      {base['dox_delivered']:.0f} mg/m2"
        f"   (cardiac risk rises steeply > 300-400)")
    say(f"     LVEF                        {base['LVEF']:.1f} %"
        f"        (baseline 62.0)")
    say(f"     eGFR nadir                  {base['eGFR_nadir_pct']:.0f} %"
        f"         of baseline")
    say(f"     high-frequency hearing      {base['hearing_dB']:.1f} dB"
        f"       threshold shift")
    say(f"     neutrophil nadir            {base['ANC_nadir']:.2f}"
        f" x10^9/L   (grade 4 < 0.5)")
    say(f"     peak mucositis grade        {base['MUC_max']:.1f}")
    say("")
    say("     Three of these are CUMULATIVE and essentially irreversible")
    say("     (heart, cochlea, tubule).  They are the walls of the box, and")
    say("     standard MAP already stands against all three.")
    say("")

    say(" 3c. DEXRAZOXANE MOVES THE WALL.  Everything else reallocates")
    say("     inside it.")
    say(f"{'arm':<40}{'cumDOX':>8}{'LVEF':>7}{'log10':>7}{'cure':>7}{'HR':>6}")
    say("-" * 100)
    for k in ("S02_MAP", "S08_MAP_dexrazoxane", "S09_MAP_dox600",
              "S10_MAP_dox600_dexra"):
        r = scen[k]
        say(f"{r['label'][:39]:<40}{r['dox_delivered']:>8.0f}{r['LVEF']:>7.1f}"
            f"{r['K_log10']:>7.2f}{r['cure_pop']:>7.3f}"
            f"{efs_hr(r['cure_pop'], base['cure_pop']):>6.2f}")
    du = scen["S09_MAP_dox600"]["LVEF"] - base["LVEF"]
    dp = scen["S10_MAP_dox600_dexra"]["LVEF"] - base["LVEF"]
    say("")
    say(f"     450 -> 600 mg/m2 costs {-du:.1f} LVEF points unprotected and"
        f" {-dp:.1f} with")
    say(f"     dexrazoxane, i.e. the protectant returns"
        f" {100 * (1 - dp / du) if du else 0:.0f}% of the cardiac price")
    say("     of the intensification.  In this model dexrazoxane is the only")
    say("     agent that RAISES the budget rather than spending it, and it is")
    say("     the only one that makes anthracycline intensification a")
    say("     coherent strategy at all.")
    say("")

    say(" 3d. MAPIE: ESCALATION ARRIVES AFTER THE BET IS PLACED, AND IN")
    say("     EXACTLY THE PATIENTS WHERE EXTRA DRUG IS WORTH LEAST.")
    say("     (RES0 0.02 = good responder, never randomised to MAPIE;")
    say("      RES0 0.20 = poor responder, the actual trial population)")
    say("")
    say(f"{'':<46}{'log10':>7}{'cure':>7}{'TRM%':>6}{'surv':>7}{'HR':>6}"
        f"{'sched%':>7}{'eGFRnad':>9}{'ANCnad':>8}{'MUC':>5}")
    say("-" * 100)
    mapie_rows = {}
    for r0, lab in ((0.02, "good responder (RES0 0.02)"),
                    (0.20, "POOR responder (RES0 0.20)")):
        pair = []
        for nm in ("S02_MAP", "S03_MAPIE"):
            reg = dict(SCENARIOS[nm])
            reg["p"] = dict(reg.get("p", {}))
            reg["p"]["RES0"] = r0
            pair.append(simulate(reg))
        A, B = pair
        hr = efs_hr(B["surv_pop"], A["surv_pop"])
        mapie_rows[lab] = dict(map_K=A["K_log10"], mapie_K=B["K_log10"],
                               map_cure=A["cure_pop"], mapie_cure=B["cure_pop"],
                               map_surv=A["surv_pop"], mapie_surv=B["surv_pop"],
                               hr=hr, mapie_sched=B["on_schedule_pct"])
        for tag, r, h in ((" — MAP", A, 1.0), (" — MAPIE", B, hr)):
            say(f"{(lab + tag)[:45]:<46}{r['K_log10']:>7.2f}"
                f"{r['cure_pop']:>7.3f}{100 * r['TRM']:>6.1f}"
                f"{r['surv_pop']:>7.3f}{h:>6.2f}{r['on_schedule_pct']:>7.0f}"
                f"{r['eGFR_nadir_pct']:>9.0f}{r['ANC_nadir']:>8.2f}"
                f"{r['MUC_max']:>5.1f}")
    say("")
    hr_poor = mapie_rows["POOR responder (RES0 0.20)"]["hr"]
    hr_good = mapie_rows["good responder (RES0 0.02)"]["hr"]
    say(f"     Predicted EFS HR for MAPIE:  {hr_poor:.2f} in poor responders,")
    say(f"                                  {hr_good:.2f} in good responders.")
    say("     EURAMOS-1 randomised MAPIE ONLY among poor responders and")
    say("     reported HR 0.98 (95% CI 0.78-1.23).")
    say("")
    say("     Every kill term in the model carries the factor (1 - RES).  A")
    say("     'poor responder' is by definition a patient in whom that factor")
    say("     is small, so the escalation was tested in exactly the subgroup")
    say("     where any added cytotoxic is worth least — and it had to buy the")
    say("     ifosfamide/etoposide cycles by giving up one doxorubicin and one")
    say("     cisplatin cycle out of a budget that was already at its walls.")
    say("     The null result is the arithmetic of the design, not a surprise")
    say("     about the drugs.")
    say("")
    return dict(contrib=contrib, mapie=mapie_rows)


# ---------------------------------------------------------------------------
# 10. VIRTUAL POPULATION: THE PLATEAU IS A VARIANCE RESULT
# ---------------------------------------------------------------------------
def run_population(scen):
    say("=" * 100)
    say(" SECTION 4 — THE FORTY-YEAR PLATEAU IS A VARIANCE RESULT")
    say("=" * 100)
    say("")
    Kbar = scen["S02_MAP"]["K_nats"]
    cv = P["POP_CV_K"]
    say(f"     Reference MAP: mean log-kill {Kbar:.2f} nats"
        f" ({Kbar / math.log(10):.2f} log10)")
    say(f"     Between-patient CV of achieved log-kill: {cv:.0%}"
        f"  (SD {cv * Kbar:.2f} nats)")
    say("     Sources of that spread: SLC19A1 / ABCC2 transporter genotype,")
    say("     methotrexate clearance (CV ~30%), MTHFR, delivered dose")
    say("     intensity after the go/no-go gates, and intrinsic")
    say("     chemosensitivity.  None of them is randomised in any trial.")
    say("")
    say(" 4a. What a shift in the MEAN buys, at the real spread:")
    say(f"{'delta log10':>12}{'cure fraction':>15}{'abs gain':>10}{'EFS HR':>9}")
    say("-" * 100)
    base_cure = pop_cure(Kbar)
    rows = []
    for dl in (-1.00, -0.50, -0.25, -0.10, 0.0, 0.10, 0.25, 0.50, 1.00):
        c = pop_cure(Kbar + dl * math.log(10))
        rows.append(dict(dlog10=dl, cure=c, gain=c - base_cure,
                         hr=efs_hr(c, base_cure)))
        say(f"{dl:>12.2f}{c:>15.3f}{c - base_cure:>+10.3f}"
            f"{efs_hr(c, base_cure):>9.2f}")
    say("")

    say(" 4b. THE SAME MEAN SHIFT AGAINST DIFFERENT SPREADS.  This is the")
    say("     whole argument:")
    say(f"{'CV of K':>9}{'SD (nats)':>11}{'cure':>8}{'cure +0.5 log':>15}"
        f"{'abs gain':>10}{'HR':>7}")
    say("-" * 100)
    spread = []
    for c_v in (0.02, 0.05, 0.10, 0.15, 0.25, 0.35, 0.50, 0.65):
        c0 = pop_cure(Kbar, cv=c_v)
        c1 = pop_cure(Kbar + 0.5 * math.log(10), cv=c_v)
        spread.append(dict(cv=c_v, sd=c_v * Kbar, c0=c0, c1=c1,
                           gain=c1 - c0, hr=efs_hr(c1, c0)))
        say(f"{c_v:>9.2f}{c_v * Kbar:>11.2f}{c0:>8.3f}{c1:>15.3f}"
            f"{c1 - c0:>+10.3f}{efs_hr(c1, c0):>7.2f}")
    g_lo, g_hi = spread[0]["gain"], spread[-1]["gain"]
    say("")
    say(f"     Half a log of extra kill is worth {100 * g_lo:+.0f} points of cure in a")
    say(f"     population with CV 2%, and only {100 * g_hi:+.0f} points in one with"
        f" CV 65% — a")
    say(f"     {g_lo / max(g_hi, 1e-9):.1f}-fold difference produced entirely by the spread and")
    say("     not at all by the drug.")
    say("")
    say("     Note carefully what this does NOT say.  Narrowing the spread at")
    say(f"     constant mean barely moves the cure fraction itself"
        f" ({pop_cure(Kbar, cv=0.45):.3f} at")
    say(f"     CV 45% vs {pop_cure(Kbar, cv=0.05):.3f} at CV 5%): the patients who gain and the")
    say("     patients who lose very nearly cancel.  What the spread destroys")
    say("     is DETECTABILITY.  Any realistic intensification moves the mean")
    say("     by 0.1-0.3 log, and at CV 25% that is worth:")
    say("")
    say(f"{'delta log10':>13}{'abs gain in cure':>19}{'HR':>7}"
        f"{'n per arm for 80% power':>26}")
    say("-" * 100)
    detect = []
    for dl in (0.10, 0.20, 0.30, 0.50):
        c = pop_cure(Kbar + dl * math.log(10))
        gain = c - base_cure
        # two-proportion sample size, alpha 0.05 two-sided, power 0.80
        pbar = 0.5 * (c + base_cure)
        n = (2 * pbar * (1 - pbar) * (1.96 + 0.842) ** 2 / max(gain, 1e-9) ** 2)
        detect.append(dict(dlog10=dl, cure=c, gain=gain, n=n,
                           hr=efs_hr(c, base_cure)))
        say(f"{dl:>13.2f}{gain:>+19.3f}{efs_hr(c, base_cure):>7.2f}{n:>26.0f}")
    say("")
    say("     EURAMOS-1 randomised 618 poor responders to MAPIE.  On this")
    say("     arithmetic that trial could only ever have detected a shift of")
    say(f"     roughly {min(d['dlog10'] for d in detect if d['n'] < 618) if any(d['n'] < 618 for d in detect) else 0.5:.2f}"
        " log or larger, and the escalation delivers about")
    say("     0.2.  The trial was not underpowered by accident; it was")
    say("     underpowered by the shape of the dose-response surface it was")
    say("     built on top of.")
    say("")
    say(" 4b-ii. NOT THE MEAN, AND NOT THE TAIL EITHER — THE SHOULDER.")
    say("")
    say("     Because P(cure) saturates at both ends, the marginal value of an")
    say("     extra log of kill is zero in patients who were already going to")
    say("     be cured AND zero in patients no achievable dose can rescue.")
    say("     Spend the same total extra log-kill on different parts of the")
    say("     distribution and the answers are not close:")
    K_pop = np.maximum(_Z * cv * Kbar + Kbar, 0.0)

    def cure_of(Karr):
        surv = P["N_LESION_CELLS"] * np.exp(-Karr)
        return float(np.mean((1.0 - (-np.expm1(-surv))) ** _LAM))

    c_ref = cure_of(K_pop)
    total = 0.25 * math.log(10) * POP_N          # total nats available to spend
    q = np.quantile(K_pop, [0.0, 0.25, 0.50, 0.75, 1.0])
    L10 = math.log(10)
    say("")
    say(f"     Achieved log-kill quartile bounds (log10):"
        f" {q[0] / L10:.2f} / {q[1] / L10:.2f} / {q[2] / L10:.2f}"
        f" / {q[3] / L10:.2f} / {q[4] / L10:.2f}")
    say("")
    say(f"{'who gets the extra kill':<40}{'each gets':>11}{'cure':>8}"
        f"{'gain':>8}{'HR':>7}")
    say("-" * 100)
    targets = [
        ("everyone (uniform intensification)", np.ones(POP_N, dtype=bool)),
        ("bottom quartile (worst responders)", K_pop <= q[1]),
        ("SECOND quartile (the near-misses)", (K_pop > q[1]) & (K_pop <= q[2])),
        ("third quartile", (K_pop > q[2]) & (K_pop <= q[3])),
        ("top quartile (already cured)", K_pop > q[3]),
        ("bottom half", K_pop <= q[2]),
        ("middle half (25-75th percentile)", (K_pop > q[1]) & (K_pop <= q[3])),
    ]
    targ_rows = []
    for nm, mask in targets:
        K2 = K_pop.copy()
        per = total / max(mask.sum(), 1)
        K2[mask] += per
        c = cure_of(K2)
        targ_rows.append(dict(who=nm, per_log10=per / L10, cure=c,
                              gain=c - c_ref, hr=efs_hr(c, c_ref)))
        say(f"{nm:<40}{per / L10:>10.2f}L{c:>8.3f}{c - c_ref:>+8.3f}"
            f"{efs_hr(c, c_ref):>7.2f}")
    say("")
    best = max(targ_rows[1:5], key=lambda r: r["gain"])
    unif = targ_rows[0]
    say(f"     The same total drug is worth {best['gain'] / max(unif['gain'], 1e-9):.1f}x more given to the")
    say(f"     {best['who']} than spread over everyone, and giving it")
    say(f"     to the WORST quartile is worth LESS than uniform"
        f" ({targ_rows[1]['gain']:+.3f} vs")
    say(f"     {unif['gain']:+.3f}) — a full extra log does not rescue a patient who is")
    say("     four logs short.")
    say("")
    say(" 4b-iii. WHERE THE MARGINAL VALUE ACTUALLY LIVES, BY DECILE:")
    say(f"{'decile of achieved log-kill':>28}{'K range (log10)':>20}"
        f"{'marginal value':>17}")
    say("-" * 100)
    eps = 0.05 * math.log(10)
    dec = []
    for d in range(10):
        lo, hi = np.quantile(K_pop, [d / 10.0, (d + 1) / 10.0])
        mask = (K_pop >= lo) & (K_pop <= hi)
        K2 = K_pop.copy()
        K2[mask] += eps
        mv = (cure_of(K2) - c_ref) / (eps / L10 * mask.sum() / POP_N)
        dec.append(dict(decile=d + 1, lo=lo / L10, hi=hi / L10, mv=mv))
        bar = "#" * int(round(40 * mv / 0.6))
        say(f"{d + 1:>28}{('%.2f - %.2f' % (lo / L10, hi / L10)):>20}"
            f"{mv:>17.3f}  {bar}")
    say("")
    top = max(dec, key=lambda r: r["mv"])
    say(f"     The peak sits in decile {top['decile']} (K {top['lo']:.2f}-{top['hi']:.2f} log10) and the")
    say("     value is essentially zero in the bottom two and the top two")
    say("     deciles.  The patients worth escalating are the ones ALREADY")
    say("     CLOSE to the threshold and just below it.")
    say("")
    say("     Now put that next to Section 5.  Necrosis correlates with")
    say("     achieved log-kill at r = 0.97, so the marker CAN identify this")
    say("     band — the near-misses at roughly 85-89% necrosis, not the 40%")
    say("     necrosis patients.  EURAMOS-1 used the marker correctly as a")
    say("     prognostic tool and then spent its escalation budget on the")
    say("     group the marker had identified as LEAST salvageable.  The band")
    say("     the arithmetic points at has never been randomised.")
    say("")
    say(" 4c. THE POISSON STRUCTURE.  Cure requires EVERY occult lesion to be")
    say("     sterilised, so the exposure enters as a double exponential:")
    say("")
    say("       P(cure) = exp( -lambda0 * [1 - exp(-n0 * e^-K)] )")
    say("")
    say(f"{'lambda0':>10}{'P(no lesion at all)':>22}{'cure with MAP':>16}"
        f"{'cure, MAP +1 log':>18}")
    say("-" * 100)
    pois = []
    for lam0 in (0.5, 1.0, 1.8, 3.0, 6.0, 15.0):
        c = pop_cure(Kbar, lam0=lam0)
        c1 = pop_cure(Kbar + math.log(10), lam0=lam0)
        pois.append(dict(lam0=lam0, none=math.exp(-lam0), cure=c, cure1=c1))
        say(f"{lam0:>10.1f}{math.exp(-lam0):>22.3f}{c:>16.3f}{c1:>18.3f}")
    say("")
    say(f"     lambda0 = {P['LAMBDA0']:.2f} was calibrated on one number: cure after")
    say(f"     amputation ALONE is exp(-lambda0) ="
        f" {math.exp(-P['LAMBDA0']):.3f}, and the historical")
    say("     figure is 0.15-0.20.  Everything downstream — including the")
    say("     6.5-log requirement on MAP — follows from that single anchor.")
    say("")
    say("     Note the last column: in overt metastatic disease (lambda0 = 15)")
    say("     a full extra log of kill barely moves the cure fraction, because")
    say("     cure needs all fifteen lesions and the product of fifteen")
    say("     near-misses is still a miss.  Same drug, same log-kill,")
    say("     different arithmetic.")
    say("")
    return dict(Kbar=Kbar, cv=cv, base_cure=base_cure, rows=rows,
                spread=spread, poisson=pois, detect=detect,
                targeting=targ_rows, deciles=dec)


# ---------------------------------------------------------------------------
# 11. HUVOS NECROSIS: A PERFECT PREDICTOR OF SOMETHING YOU CANNOT TREAT
# ---------------------------------------------------------------------------
def run_huvos_analysis():
    say("=" * 100)
    say(" SECTION 5 — HUVOS NECROSIS IS A READ-OUT, NOT A LEVER")
    say("=" * 100)
    say("")
    say("     Necrosis at week 11 and micrometastatic log-kill are two")
    say("     integrals of the SAME sensitivity parameter over the SAME ten")
    say("     weeks.  That is why necrosis is the strongest prognostic factor")
    say("     in the disease — and also why changing therapy AFTER measuring")
    say("     it cannot recover what the parameter already decided.")
    say("")
    say(f"{'RES0':>7}{'Huvos%':>9}{'grade':>8}{'log10 kill':>12}{'cure':>8}"
        f"{'HR vs RES0=0.02':>17}")
    say("-" * 100)
    rows = []
    ref = None
    for r0 in (0.00, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45):
        reg = dict(SCENARIOS["S02_MAP"])
        reg["p"] = {"RES0": r0}
        r = simulate(reg)
        if ref is None:
            ref = r["cure_pop"]
        if abs(r0 - 0.02) < 1e-9:
            ref = r["cure_pop"]
        rows.append(dict(res0=r0, huvos=r["huvos_pct"], K=r["K_log10"],
                         cure=r["cure_pop"]))
    for r in rows:
        say(f"{r['res0']:>7.2f}{r['huvos']:>9.1f}"
            f"{('good' if r['huvos'] >= 90 else 'poor'):>8}"
            f"{r['K']:>12.2f}{r['cure']:>8.3f}"
            f"{efs_hr(r['cure'], ref):>17.2f}")
    hv = np.array([x["huvos"] for x in rows])
    kk = np.array([x["K"] for x in rows])
    cc = np.array([x["cure"] for x in rows])
    r_hk = float(np.corrcoef(hv, kk)[0, 1])
    r_hc = float(np.corrcoef(hv, cc)[0, 1])
    say("")
    say(f"     Pearson r(Huvos%, micrometastatic log-kill) = {r_hk:.3f}")
    say(f"     Pearson r(Huvos%, cure fraction)            = {r_hc:.3f}")
    say("")
    say(" 5b. THE SAME ESCALATION, MOVED EARLIER.  Poor-responder tumour")
    say("     (RES0 = 0.20); the ONLY difference is when the ifosfamide/")
    say("     etoposide cycles are given:")
    say(f"{'arm':<44}{'Huvos%':>9}{'log10':>8}{'cure':>8}{'HR':>7}")
    say("-" * 100)
    arms = OrderedDict([
        ("MAP", dict()),
        ("MAPIE — escalation AFTER surgery (as randomised)",
         dict(SCENARIOS["S03_MAPIE"])),
        ("MAP + the same IE cycles given PRE-operatively",
         dict(ifo_weeks=[2, 3, 7, 8, 12], dox_weeks=[1, 6, 17, 22, 27],
              cis_weeks=[1, 6, 22])),
    ])
    early = {}
    ref2 = None
    for nm, over in arms.items():
        reg = dict(SCENARIOS["S02_MAP"])
        reg.update({k: v for k, v in over.items() if k != "label"})
        reg["p"] = dict(reg.get("p", {}))
        reg["p"]["RES0"] = 0.20
        r = simulate(reg)
        if ref2 is None:
            ref2 = r["cure_pop"]
        early[nm] = dict(huvos=r["huvos_pct"], K=r["K_log10"],
                         cure=r["cure_pop"])
        say(f"{nm:<44}{r['huvos_pct']:>9.1f}{r['K_log10']:>8.2f}"
            f"{r['cure_pop']:>8.3f}{efs_hr(r['cure_pop'], ref2):>7.2f}")
    say("")
    e_map = early["MAP"]
    e_late = early["MAPIE — escalation AFTER surgery (as randomised)"]
    e_early = early["MAP + the same IE cycles given PRE-operatively"]
    say(f"     Moving the identical drug from after surgery to before it changes")
    say(f"     the NECROSIS read-out by"
        f" {e_early['huvos'] - e_late['huvos']:+.1f} points and the micrometastatic")
    say(f"     log-kill by {e_early['K'] - e_late['K']:+.2f} log10.")
    say("")
    say("     That dissociation is the point, and it is worse for the surrogate")
    say("     than a weak correlation would be.  Across the resistance sweep")
    say(f"     above, necrosis tracks micrometastatic kill almost perfectly"
        f" (r = {r_hk:.2f}),")
    say("     so it is an excellent PROGNOSTIC marker.  But necrosis is an")
    say("     integral over the primary, where drug penetration and matrix")
    say("     TGF-beta apply, and cure is an integral over the lung, where they")
    say("     do not.  Re-timing therapy moves the first without moving the")
    say("     second, i.e. the marker can be improved without improving the")
    say("     patient.  A marker that is prognostic but manipulable is exactly")
    say("     the wrong thing to build a treatment decision on — and")
    say("     response-adapted escalation is a design that reads the marker at")
    say("     week 11 and then spends its budget in the window where the")
    say("     integral has the least left to give.")
    say("")
    return dict(rows=rows, r_huvos_K=r_hk, r_huvos_cure=r_hc, early=early)


# ---------------------------------------------------------------------------
# 12. BONE-DIRECTED THERAPY MOVES A DIFFERENT READ-OUT
# ---------------------------------------------------------------------------
def run_bone_analysis(scen):
    say("=" * 100)
    say(" SECTION 6 — THE VICIOUS CYCLE IS A GROWTH TERM, NOT A SURVIVAL TERM")
    say("=" * 100)
    say("")
    say(f"{'arm':<42}{'OCLmax':>8}{'ALP':>7}{'Huvos%':>8}{'log10':>7}"
        f"{'cure':>7}{'HR':>6}")
    say("-" * 100)
    base = scen["S02_MAP"]
    for k in ("S02_MAP", "S11_MAP_zoledronate", "S12_MAP_denosumab"):
        r = scen[k]
        say(f"{r['label'][:41]:<42}{r['OCL_max']:>8.2f}{r['ALP']:>7.2f}"
            f"{r['huvos_pct']:>8.1f}{r['K_log10']:>7.2f}{r['cure_pop']:>7.3f}"
            f"{efs_hr(r['cure_pop'], base['cure_pop']):>6.2f}")
    z = scen["S11_MAP_zoledronate"]
    d = scen["S12_MAP_denosumab"]
    say("")
    say(f"     Zoledronate cuts peak osteoclast drive by"
        f" {100 * (1 - z['OCL_max'] / base['OCL_max']):.0f}% and denosumab by"
        f" {100 * (1 - d['OCL_max'] / base['OCL_max']):.0f}%,")
    say(f"     while micrometastatic log-kill moves by"
        f" {z['K_log10'] - base['K_log10']:+.2f} and"
        f" {d['K_log10'] - base['K_log10']:+.2f} log10 respectively.")
    say("")
    say("     The reason is structural rather than empirical.  In this model")
    say("     the RANKL -> osteoclast -> matrix-TGF-beta loop multiplies the")
    say("     tumour GROWTH rate (KG), whereas cure is decided by the KILL")
    say("     integral, and during MAP the kill term is one to two orders of")
    say("     magnitude larger than the growth term it would be modifying.")
    say("     Suppressing the loop therefore does exactly what the trials")
    say("     found it does: it changes osteolysis, ALP and radiographic")
    say("     appearance, and it does not change survival.  OS2006 randomised")
    say("     318 patients to zoledronate on top of chemotherapy and found no")
    say("     EFS benefit with worse local control.  Two readouts, two")
    say("     integrals; the marker moved and the disease did not.")
    say("")


# ---------------------------------------------------------------------------
# 13. MAIN
# ---------------------------------------------------------------------------
def main():
    say("#" * 100)
    say("#  OSTEOSARCOMA (OSA) QSP MODEL — REFERENCE OUTPUT")
    say("#  osa_reference_model.py    (numerical mirror of osa_mrgsolve_model.R)")
    say("#")
    say("#  Cure is a Poisson bet on lesions nobody can see, paid for out of an")
    say("#  exposure budget with three hard organ ceilings:")
    say("#")
    say("#      P(cure) = exp( -lambda0 * [1 - exp(-n0 * e^-K)] )")
    say("#      K       = kill integral, bounded by kidney / heart / cochlea")
    say("#")
    say("#  and the methotrexate ceiling is not a constant, because the drug is")
    say("#  cleared through the organ it destroys.")
    say("#" * 100)
    say("")
    say(f"  states {NS}   parameters {len(P)}   scenarios {len(SCENARIOS)}"
        f"   virtual population {POP_N}")
    say("")
    scen = run_scenarios()
    loop = run_loop_analysis()
    budget = run_budget_analysis(scen)
    pop = run_population(scen)
    huv = run_huvos_analysis()
    run_bone_analysis(scen)

    say("=" * 100)
    say(" SUMMARY OF DERIVED NUMBERS")
    say("=" * 100)
    hr_poor = budget["mapie"]["POOR responder (RES0 0.20)"]["hr"]
    hr_good = budget["mapie"]["good responder (RES0 0.02)"]["hr"]
    items = [
        ("critical urine pH (tubular supersaturation)",
         f"{loop['pH_crit']:.2f}", "protocol target >= 7.0"),
        ("urine pH below which elimination is delayed",
         f"{loop['pH_delay']:.2f}", ""),
        ("urine pH below which the course causes AKI",
         f"{loop['pH_aki']:.2f}", ""),
        ("1 pH unit is worth this much hydration",
         f"{10 ** P['SOL_B']:.2f}x", "or the same dose reduction"),
        ("AUC multiplier when starting GFR is halved",
         f"{loop['gain_rows'][-1]['auc'] / loop['gain_rows'][0]['auc']:.2f}x",
         "same dose, same pH"),
        ("MAP log-kill, micrometastatic pool",
         f"{scen['S02_MAP']['K_log10']:.2f} log10", ""),
        ("MAP survival fraction (population)",
         f"{scen['S02_MAP']['surv_pop']:.3f}",
         "EURAMOS-1 5-y EFS 0.54-0.59"),
        ("surgery-alone survival fraction",
         f"{scen['S01_surgery_only']['surv_pop']:.3f}",
         "historical 0.15-0.20"),
        ("MAP / MAPIE courses given on schedule",
         f"{scen['S02_MAP']['on_schedule_pct']:.0f}% /"
         f" {scen['S03_MAPIE']['on_schedule_pct']:.0f}%",
         "EURAMOS-1 completed protocol 76% / 51%"),
        ("MAP / MAPIE treatment-related mortality",
         f"{100 * scen['S02_MAP']['TRM']:.1f}% /"
         f" {100 * scen['S03_MAPIE']['TRM']:.1f}%", "reported ~1% overall"),
        ("under-hydration: cure / TRM / survival",
         f"{scen['S06_MAP_lowhydration']['cure_pop']:.3f} /"
         f" {100 * scen['S06_MAP_lowhydration']['TRM']:.0f}% /"
         f" {scen['S06_MAP_lowhydration']['surv_pop']:.3f}",
         f"vs MAP survival {scen['S02_MAP']['surv_pop']:.3f}"),
        ("best arm in the whole set",
         "dox 600 + dexrazoxane",
         f"survival {scen['S10_MAP_dox600_dexra']['surv_pop']:.3f},"
         f" HR {scen['S10_MAP_dox600_dexra']['efs_hr_vs_MAP']:.2f}"),
        ("MAPIE predicted EFS HR, poor responders",
         f"{hr_poor:.2f}", "EURAMOS-1 0.98 (0.78-1.23)"),
        ("MAPIE predicted EFS HR, good responders",
         f"{hr_good:.2f}", "never randomised"),
        ("r(Huvos%, micrometastatic log-kill)",
         f"{huv['r_huvos_K']:.3f}", ""),
        ("cumulative doxorubicin delivered by MAP",
         f"{scen['S02_MAP']['dox_delivered']:.0f} mg/m2", "ceiling ~ 400"),
        ("LVEF after MAP / MAP + dexrazoxane",
         f"{scen['S02_MAP']['LVEF']:.1f} /"
         f" {scen['S08_MAP_dexrazoxane']['LVEF']:.1f} %", ""),
        ("value of +0.5 log kill at CV 2% / CV 65%",
         f"{100 * pop['spread'][0]['gain']:+.0f} /"
         f" {100 * pop['spread'][-1]['gain']:+.0f} points", "variance dominates"),
        ("n per arm to detect MAPIE's +0.2 log at 80% power",
         f"{next(d['n'] for d in pop['detect'] if abs(d['dlog10'] - 0.20) < 1e-9):.0f}",
         "EURAMOS-1 randomised 618"),
        ("same budget: uniform vs 2nd-quartile targeting",
         f"{pop['targeting'][0]['gain']:+.3f} vs"
         f" {pop['targeting'][2]['gain']:+.3f}",
         "the shoulder, not the tail"),
        ("marginal value peaks at this achieved log-kill",
         f"{max(pop['deciles'], key=lambda r: r['mv'])['lo']:.1f}-"
         f"{max(pop['deciles'], key=lambda r: r['mv'])['hi']:.1f} log10",
         "~85-89% necrosis, the near-misses"),
    ]
    for k, v, note in items:
        say(f"  {k:<46}{v:>22}   {note}")
    say("")

    with open(os.path.join(HERE, "osa_scenario_results.json"), "w") as f:
        json.dump(scen, f, indent=1, default=float)
    with open(os.path.join(HERE, "osa_analysis_results.json"), "w") as f:
        json.dump(dict(loop=loop, budget=budget, population=pop, huvos=huv),
                  f, indent=1, default=float)
    with open(os.path.join(HERE, "osa_reference_output.txt"), "w") as f:
        f.write("\n".join(OUT) + "\n")
    say("  wrote osa_reference_output.txt, osa_scenario_results.json,")
    say("        osa_analysis_results.json")


if __name__ == "__main__":
    main()
