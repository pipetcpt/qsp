#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pmdd_python_twin.py — dependency-free (stdlib-only) numerical twin of the
PMDD QSP model in `pmdd_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
Every quantitative statement in `README.md` is produced by running this file.
It implements the *same* 60-state ODE system, the *same* parameters and the
*same* dosing regimens as the mrgsolve model, using a fixed-step RK4
integrator and nothing but the Python standard library, so the numbers can be
re-checked on any machine without R, mrgsolve, numpy or scipy.

    python3 pmdd_python_twin.py             # all scenarios + the check table
    python3 pmdd_python_twin.py --csv out   # also dump per-scenario CSVs

THE STRUCTURAL CLAIM
--------------------
PMDD is modelled as a GAIN disorder, not a level disorder, because the level
story is already falsified in the clinic: luteal E2/P4/allopregnanolone (ALLO)
concentrations in PMDD are indistinguishable from those of asymptomatic
controls, yet the same physiological add-back reinstates symptoms in PMDD and
does nothing in controls.  Two structural commitments follow.

  (1) NON-MONOTONIC transduction.  The affective drive produced by
      neurosteroid load L at GABA-A is unimodal,

          NS_DRIVE = SENS * (L/KP) * exp(1 - L/KP)          [peak at L = KP]

      so the disease parameter is KP — WHERE THE PEAK SITS relative to the
      physiological luteal ALLO range — not the range itself.

  (2) A SLOW, SYMMETRIC CHANGE DETECTOR sets the gain SENS.  GABA-A subunit
      composition (delta from ALLO exposure history, alpha4 from the RATE OF
      CHANGE of ALLO and of total progestogen tone) moves on a timescale of
      days while ALLO itself turns over in minutes.  The detector is
      symmetric: the onset of exposure drives it as well as the withdrawal.

Five consequences are checked below.  Each disappears or inverts in the
`fals_*` runs, where SENS is frozen at 1, the transduction is made monotone
and every rate detector is replaced by a level detector:

  [P1] the symptom peak LAGS the ALLO peak by ~4-6 days (late luteal) instead
       of coinciding with it;
  [P2] dutasteride 2.5 mg helps and 0.5 mg does not — the top of the
       inverted U is locally FLAT, so a small ALLO reduction buys nothing;
  [P3] a SEDATIVE-dose exogenous neurosteroid helps (descending limb) while a
       SUB-sedative dose of the same agent makes things worse (ascending
       limb) — opposite signs from one drug, which no monotone model can give;
  [P4] a 4-day hormone-free interval beats a 7-day one, because the driver is
       the withdrawal RATE and not the trough LEVEL;
  [P5] E2/P4 add-back on top of GnRH-agonist suppression reinstates symptoms
       TRANSIENTLY and they then remit under continued, identical hormone
       exposure — the signature of adaptation, which a level detector cannot
       produce.
"""

import argparse
import math
import os
import sys

# ══════════════════════════════════════════════════════════════════════════
# 1. STATE VECTOR
# ══════════════════════════════════════════════════════════════════════════
STATES = [
    # --- HPO cycle engine (6) -------------------------------------------
    "FSH", "LH", "FOLL", "CL", "E2", "P4",
    # --- neurosteroidogenesis (5) ---------------------------------------
    "DHP", "ALLOP", "ALLOB", "ISOB", "ZUR",
    # --- GABA-A composition plasticity + change detectors (5) -----------
    "MALLO", "MPROG", "DELTA", "ALPHA4", "ME2",
    # --- serotonergic system (3) ----------------------------------------
    "SHT", "DES", "BDNF",
    # --- corticolimbic (3) ----------------------------------------------
    "AMY", "PFC", "REW",
    # --- HPA axis (3) ---------------------------------------------------
    "CRH", "ACTH", "CORT",
    # --- somatic / peripheral (4) ---------------------------------------
    "ALDO", "ECFV", "BREAST", "PGE",
    # --- DRSP symptom domains (11) --------------------------------------
    "S_IRR", "S_DEP", "S_ANX", "S_LAB", "S_ANH", "S_COG",
    "S_FAT", "S_APP", "S_SLP", "S_OVR", "S_PHY",
    # --- organ-system / safety (3) --------------------------------------
    "HF", "BMD", "ENDO",
    # --- drug PK (15) ---------------------------------------------------
    "SERT_D", "SERT_C", "DMS_C",          # sertraline depot/central/metabolite
    "FLX_C", "NFLX_C",                    # fluoxetine + norfluoxetine
    "DRSP_D", "DRSP_C",                   # drospirenone depot/central
    "LEUP_D", "LEUP_C",                   # leuprolide depot/central
    "SEPRA_D", "SEPRA_C",                 # sepranolone SC depot/central
    "DUTA_D", "DUTA_C",                   # dutasteride depot/central
    "ALPZ_D", "ALPZ_C",                   # alprazolam depot/central
    "UPA_C",                              # ulipristal acetate
    # --- pituitary GnRH-receptor downregulation (1) ---------------------
    "PITDS",
]
IX = {s: i for i, s in enumerate(STATES)}
NST = len(STATES)
CORE_DOMAINS = ("S_IRR", "S_DEP", "S_ANX", "S_LAB", "S_ANH")
ALL_DOMAINS = CORE_DOMAINS + ("S_COG", "S_FAT", "S_APP", "S_SLP", "S_OVR", "S_PHY")
DRSP_FLOOR = float(len(ALL_DOMAINS))          # all domains scored 1 -> 11
CORE_FLOOR = float(len(CORE_DOMAINS))         # -> 5

# ══════════════════════════════════════════════════════════════════════════
# 2. PARAMETERS
# ══════════════════════════════════════════════════════════════════════════
P = dict(
    # ---- cycle engine -------------------------------------------------
    TCYC=28.0,
    KFSH_IN=2.4, KFSH_OUT=0.55, KINH=0.30,
    KLH_IN=1.6, KLH_OUT=0.9, SURGE_AMP=26.0, SURGE_DAY=13.6, SURGE_SD=0.55,
    KFG=0.42, REC=0.030, KFSH_M=1.20, FSH_THR=0.90, FMAX=1.0, KFL=0.06,
    KOV=20.0, OVUL_DAY=14.0, OVUL_SD=0.45,
    KCL=0.50, KCL_G=1.00, CL_MAX=3.10,
    KLUT0=0.05, KLUT1=1.25, LUT_REG_DAY=24.0, LUT_REG_SD=0.75,
    KE2F=3520.0, E2_HILL=3.0, KE2L=60.0, KE2_OUT=2.2, E2_BASE=30.0,
    KP4=7.6, KP4_OUT=2.0, P4_BASE=0.25, KP4_FB=3.0,
    # ---- neurosteroidogenesis -----------------------------------------
    K5AR=3.0, KDHP=3.2, K3HSD=3.0, KALLOP=7.0, ADR_ALLO=8.26,
    KIN_B=7.2, KOUT_B=8.0, DENOVO=1.0, W_STRESS_NS=0.55,
    IC50_5AR1=110.0,           # ng/mL dutasteride at brain SRD5A1
    EMAX_HSD=0.55, EC50_HSD=30.0,    # SSRI shift of 3alpha-HSD (raises ALLO)
    KISO_OUT=5.0, KI_ISO=0.80,       # sepranolone antagonism — FITTED (see README)
    K_ZUR=0.83,                # exogenous neurosteroid analogue, t1/2 ~20 h
    # ---- GABA-A plasticity: the slow, SYMMETRIC change detector -------
    TAU_M=1.0, TAU_MP=1.0, TAU_D=5.0, TAU_A4=1.2,
    HD=1.25, KD=3.0,           # delta upregulation by ALLO exposure history
    HA=2.20, KA=0.55,          # alpha4 upregulation by the RATE of change
    KRISE=0.55,                # weight of a RISING signal vs a falling one
    W_PROG=0.85,               # weight of progestogen-tone change vs ALLO
    # GAMMA_D is NEGATIVE: delta-subunit upregulation under sustained exposure
    # is the TOLERANCE arm (chronic-progesterone pseudopregnancy data), while
    # withdrawal-induced alpha4 upregulation is the SENSITISING arm.  This is
    # what makes continuous, stable hormone exposure asymptomatic while a
    # change in exposure is symptomatic.
    GAMMA_D=-0.35, GAMMA_A4=1.50,    # PMDD gains (controls: -0.05 / 0.25)
    # ---- inverted-U transduction --------------------------------------
    KP=5.0,                    # PMDD: the peak sits inside the luteal range
    MONOTONE=0.0, RATE_OFF=0.0, MONO_GAIN=1.66,
    NS_HILL=4.0, NS50=1.50,
    BZEQ=0.12, W_BZ_AMY=1.30,  # alprazolam: BZ-site load + direct anxiolysis
    PREG_EQ=0.35,              # pregnanolone/THDOC additive load (nM)
    # ---- E2 withdrawal ------------------------------------------------
    TAU_ME2=6.0,
    # ---- serotonin ----------------------------------------------------
    KS_SYN=1.0, KS_OUT=1.0, OCC_MAX=0.85, EC50_SERT=13.0,
    POT_FLX=0.50, POT_DMS=0.05,
    WB=0.50, TAU_DES=12.0, KD2=1.6,
    KBDNF=1.10, TAU_BDNF=14.0,
    K5HT_NS=0.28,              # luteal 5-HT suppression by the NS drive
    # ---- corticolimbic ------------------------------------------------
    TAU_AMY=0.5, TAU_PFC=0.8, TAU_REW=1.2,
    W_NS_AMY=1.50, W_CORT_AMY=0.45, W_SHT_AMY=1.25, W_PFC_AMY=0.70,
    W_PROG_AMY=0.30,           # progestin-related mood side effect
    P_NS_PFC=0.34, P_E2_PFC=0.10, CBT=0.0,
    R_NS_REW=0.45, R_BDNF=0.30, R_WE2=0.15,
    # ---- HPA ----------------------------------------------------------
    TAU_CRH=0.2, TAU_ACTH=0.1, TAU_CORT=0.3,
    STRESS=0.30, K_ALLO_BRAKE=0.45, CORT0=1.0,
    # ---- somatic ------------------------------------------------------
    TAU_ALDO=1.0, A_P4=0.55,
    TAU_ECFV=1.5, ECF_GAIN=0.030,
    KI_MR_DRSP=9.0, KI_MR_SPIRO=1.0, KI_MR_P4=22.0,
    TAU_BR=1.5, B_E2=0.55, B_P4=0.55, E2REF=120.0, P4REF=11.0,
    TAU_PG=0.8, KPG=1.0, NSAID=0.0, UPA_PR_KI=8.0, KI_PR=9.0,
    # ---- symptom domains ----------------------------------------------
    TAU_S=0.6, BASE_SYM=0.08, CA_SUPP=0.0,
    # ---- organ / safety -----------------------------------------------
    TAU_HF=2.0, HF_MAX=9.0, E50_HF=38.0, HF_SSRI=0.45,
    KBMD_LOSS=0.040, E50_BMD=45.0, HYPO_REF=0.15, KBMD_REC=0.010,
    TAU_ENDO=4.0,
    # ---- PK -----------------------------------------------------------
    KA_SERT=14.4, V_SERT=1500.0, K_SERT=0.64, K_DMS=0.25, FM_DMS=0.35,
    V_FLX=2500.0, K_FLX=0.173, K_NFLX=0.077, FM_NFLX=0.60,
    KA_DRSP=9.0, V_DRSP=250.0, K_DRSP=0.533,
    EC50_SUPP=4.0, EE_E2EQ=72.0,
    KA_LEUP=0.075, V_LEUP=25.0, K_LEUP=2.8,
    K_DS_ON=1.60, K_DS_OFF=0.045, AGON_FLARE=2.6, EC50_LEUP=0.9,
    KA_SEPRA=2.2, V_SEPRA=90.0, K_SEPRA=4.5, SEPRA_BRAIN=0.55,
    KA_DUTA=6.0, V_DUTA=500.0, K_DUTA=0.0198,
    KA_ALPZ=12.0, V_ALPZ=70.0, K_ALPZ=1.386,
    K_UPA=0.35, EC50_UPA_SUPP=6.0,
    # ---- exogenous hormone add-back / surgery -------------------------
    ADDBACK_E2=0.0, ADDBACK_P4=0.0, ADDBACK_T0=0.0,
    OOPH=0.0, SPIRO=0.0,
)

# An asymptomatic control differs from PMDD in exactly THREE parameters:
# where the response peak sits (KP) and the two plasticity gains.
CONTROL_OVERRIDES = dict(GAMMA_D=-0.05, GAMMA_A4=0.25, KP=13.0)


# ══════════════════════════════════════════════════════════════════════════
# 3. DOSING
# ══════════════════════════════════════════════════════════════════════════
class Regimen(object):
    """A dosing regimen: a sorted list of (time, state index, amount) boluses."""

    def __init__(self):
        self.events = []

    def add(self, state, amt, start, stop, every, cyclic_window=None):
        """Repeated boluses; cyclic_window=(a,b) restricts to cycle days a<=cd<b."""
        t, n = start, 0
        while t < stop and n < 40000:
            if cyclic_window is None:
                self.events.append((t, IX[state], amt))
            else:
                cd = math.fmod(t, P["TCYC"]) + 1.0
                if cyclic_window[0] <= cd < cyclic_window[1]:
                    self.events.append((t, IX[state], amt))
            t += every
            n += 1
        return self

    def finalize(self):
        self.events.sort(key=lambda e: e[0])
        return self


# ══════════════════════════════════════════════════════════════════════════
# 4. HELPERS
# ══════════════════════════════════════════════════════════════════════════
def gauss(x, mu, sd):
    z = (x - mu) / sd
    return 0.0 if z * z > 60.0 else math.exp(-0.5 * z * z)


def sigmoid(x, x0, k):
    z = (x - x0) / k
    if z > 40.0:
        return 1.0
    if z < -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def clamp01(x):
    return 0.0 if x < 0.0 else (1.0 if x > 1.0 else x)


def pg_window(cd):
    """Perimenstrual prostaglandin window (cycle day ~27 through day ~3)."""
    if cd >= 26.5:
        return sigmoid(cd, 27.0, 0.5)
    if cd <= 3.5:
        return 1.0 - sigmoid(cd, 3.0, 0.6)
    return 0.0


# ══════════════════════════════════════════════════════════════════════════
# 5. THE ODE SYSTEM
# ══════════════════════════════════════════════════════════════════════════
def derivatives(t, y, p, out):
    cd = math.fmod(t, p["TCYC"]) + 1.0                     # cycle day, 1..29

    (FSH, LH, FOLL, CL, E2, P4,
     DHP, ALLOP, ALLOB, ISOB, ZUR,
     MALLO, MPROG, DELTA, ALPHA4, ME2,
     SHT, DES, BDNF,
     AMY, PFC, REW,
     CRH, ACTH, CORT,
     ALDO, ECFV, BREAST, PGE,
     S_IRR, S_DEP, S_ANX, S_LAB, S_ANH, S_COG,
     S_FAT, S_APP, S_SLP, S_OVR, S_PHY,
     HF, BMD, ENDO,
     SERT_D, SERT_C, DMS_C, FLX_C, NFLX_C,
     DRSP_D, DRSP_C, LEUP_D, LEUP_C, SEPRA_D, SEPRA_C,
     DUTA_D, DUTA_C, ALPZ_D, ALPZ_C, UPA_C, PITDS) = y

    # ─────────────────────────────────────────────────────────────── PK ──
    dSERT_D = -p["KA_SERT"] * SERT_D
    dSERT_C = p["KA_SERT"] * SERT_D / p["V_SERT"] * 1000.0 - p["K_SERT"] * SERT_C
    dDMS_C = p["FM_DMS"] * p["K_SERT"] * SERT_C - p["K_DMS"] * DMS_C
    dFLX_C = -p["K_FLX"] * FLX_C
    dNFLX_C = p["FM_NFLX"] * p["K_FLX"] * FLX_C - p["K_NFLX"] * NFLX_C
    dDRSP_D = -p["KA_DRSP"] * DRSP_D
    dDRSP_C = p["KA_DRSP"] * DRSP_D / p["V_DRSP"] * 1000.0 - p["K_DRSP"] * DRSP_C
    dLEUP_D = -p["KA_LEUP"] * LEUP_D
    dLEUP_C = p["KA_LEUP"] * LEUP_D / p["V_LEUP"] * 1000.0 - p["K_LEUP"] * LEUP_C
    dSEPRA_D = -p["KA_SEPRA"] * SEPRA_D
    dSEPRA_C = p["KA_SEPRA"] * SEPRA_D / p["V_SEPRA"] * 1000.0 - p["K_SEPRA"] * SEPRA_C
    dDUTA_D = -p["KA_DUTA"] * DUTA_D
    dDUTA_C = p["KA_DUTA"] * DUTA_D / p["V_DUTA"] * 1000.0 - p["K_DUTA"] * DUTA_C
    dALPZ_D = -p["KA_ALPZ"] * ALPZ_D
    dALPZ_C = p["KA_ALPZ"] * ALPZ_D / p["V_ALPZ"] * 1000.0 - p["K_ALPZ"] * ALPZ_C
    dUPA_C = -p["K_UPA"] * UPA_C
    dZUR = -p["K_ZUR"] * ZUR

    # ───────────────────────────────────── gonadotropin suppression ─────
    supp_oc = DRSP_C / (DRSP_C + p["EC50_SUPP"])
    supp_upa = UPA_C / (UPA_C + p["EC50_UPA_SUPP"])
    leup_norm = LEUP_C / (LEUP_C + p["EC50_LEUP"])
    dPITDS = p["K_DS_ON"] * leup_norm * (1.0 - PITDS) - p["K_DS_OFF"] * PITDS
    # GnRH agonist: flare first (PITDS still ~0), profound suppression later
    gnrh_eff = (1.0 - PITDS) * (1.0 + p["AGON_FLARE"] * leup_norm * (1.0 - PITDS))
    gnrh_eff *= (1.0 - supp_oc) * (1.0 - 0.55 * supp_upa)
    gnrh_eff /= (1.0 + max(P4, 0.0) / p["KP4_FB"])          # luteal P4 feedback
    gnrh_eff = max(gnrh_eff, 0.0)

    # ─────────────────────────────────────────── HPO cycle engine ───────
    dFSH = p["KFSH_IN"] * gnrh_eff / (1.0 + FOLL / p["KINH"]) - p["KFSH_OUT"] * FSH
    surge = p["SURGE_AMP"] * gauss(cd, p["SURGE_DAY"], p["SURGE_SD"]) \
        * clamp01((FOLL - 0.25) / 0.35) * gnrh_eff
    dLH = p["KLH_IN"] * gnrh_eff + surge - p["KLH_OUT"] * LH

    gonad = 0.0 if p["OOPH"] > 0.5 else 1.0
    ovul_flux = p["KOV"] * FOLL * gauss(cd, p["OVUL_DAY"], p["OVUL_SD"]) \
        * clamp01((LH - 2.5) / 6.0) * gonad
    # a fresh follicular cohort is recruited every cycle (REC) and then grows
    # autocatalytically; the FSH response saturates (KFSH_M) so that the very
    # high early-follicular FSH cannot drive explosive growth
    hfsh = max(0.0, FSH - p["FSH_THR"])
    hfsh = hfsh / (p["KFSH_M"] + hfsh)
    dFOLL = (p["REC"] + p["KFG"] * FOLL) * hfsh * (1.0 - FOLL / p["FMAX"]) * gonad \
        - p["KFL"] * FOLL - ovul_flux

    klut = p["KLUT0"] + p["KLUT1"] * sigmoid(cd, p["LUT_REG_DAY"], p["LUT_REG_SD"])
    # CL growth is gated to the post-ovulatory window so that a numerically
    # small CL remnant cannot regrow during the follicular phase
    w_lut = sigmoid(cd, 14.8, 0.5) * (1.0 - sigmoid(cd, 24.5, 0.7))
    grow = p["KCL_G"] * CL * (1.0 - CL / p["CL_MAX"]) * w_lut if CL > 0.02 else 0.0
    dCL = p["KCL"] * ovul_flux + grow - klut * CL

    ab = 1.0 if t >= p["ADDBACK_T0"] else 0.0
    dE2 = p["KE2F"] * FOLL ** p["E2_HILL"] + p["KE2L"] * CL \
        - p["KE2_OUT"] * (E2 - p["E2_BASE"]) + ab * p["ADDBACK_E2"]
    dP4 = p["KP4"] * CL - p["KP4_OUT"] * (P4 - p["P4_BASE"]) + ab * p["ADDBACK_P4"]

    # ethinylestradiol is not in the E2 state but it is estrogenic: every
    # estrogen-receptor-mediated readout below therefore sees E2_TOT
    E2_TOT = E2 + p["EE_E2EQ"] * (DRSP_C / (DRSP_C + 5.0))

    # ──────────────────────────────────── neurosteroidogenesis ──────────
    a5ar = 1.0 / (1.0 + DUTA_C / p["IC50_5AR1"])
    ssri_eq = SERT_C + p["POT_DMS"] * DMS_C + p["POT_FLX"] * (FLX_C + NFLX_C)
    f_hsd = 1.0 + p["EMAX_HSD"] * ssri_eq / (p["EC50_HSD"] + ssri_eq)

    dDHP = p["K5AR"] * a5ar * max(P4, 0.0) - p["KDHP"] * DHP
    dALLOP = p["K3HSD"] * f_hsd * DHP + p["ADR_ALLO"] * a5ar - p["KALLOP"] * ALLOP
    dALLOB = p["KIN_B"] * ALLOP \
        + p["DENOVO"] * (1.0 + p["W_STRESS_NS"] * p["STRESS"]) * a5ar \
        - p["KOUT_B"] * ALLOB
    dISOB = p["SEPRA_BRAIN"] * SEPRA_C - p["KISO_OUT"] * ISOB

    # ──────────────────── GABA-A composition: SYMMETRIC change detector ──
    # total progestogen tone = endogenous P4 + any exogenous progestin
    PROG = max(P4, 0.0) / p["P4REF"] + 1.4 * DRSP_C / (DRSP_C + p["KI_PR"])
    ALLO_EFF = ALLOB + ZUR          # exogenous analogues are sensed too
    dMALLO = (ALLO_EFF - MALLO) / p["TAU_M"]
    dMPROG = (PROG - MPROG) / p["TAU_MP"]
    if p["RATE_OFF"] > 0.5:
        chg, WE2 = 0.0, 0.0
    else:
        fall_a = max(0.0, (MALLO - ALLO_EFF) / (MALLO + 0.5))
        rise_a = max(0.0, (ALLO_EFF - MALLO) / (MALLO + 0.5))
        fall_p = max(0.0, (MPROG - PROG) / (MPROG + 0.25))
        rise_p = max(0.0, (PROG - MPROG) / (MPROG + 0.25))
        chg = (fall_a + p["KRISE"] * rise_a) \
            + p["W_PROG"] * (fall_p + p["KRISE"] * rise_p)
        WE2 = max(0.0, (ME2 - E2_TOT) / (ME2 + 20.0))
    delta_ss = 1.0 + p["HD"] * ALLOB / (p["KD"] + ALLOB)
    a4_ss = 1.0 + p["HA"] * chg / (p["KA"] + chg)
    dDELTA = (delta_ss - DELTA) / p["TAU_D"]
    dALPHA4 = (a4_ss - ALPHA4) / p["TAU_A4"]
    dME2 = (E2_TOT - ME2) / p["TAU_ME2"]

    # ───────────────────── inverted-U neurosteroid transduction ─────────
    sens = 1.0 + p["GAMMA_D"] * (DELTA - 1.0) + p["GAMMA_A4"] * (ALPHA4 - 1.0)
    bz = p["BZEQ"] * ALPZ_C / (1.0 + 0.60 * (ALPHA4 - 1.0))
    load = (ALLOB + p["PREG_EQ"] + bz + ZUR) / (1.0 + ISOB / p["KI_ISO"])
    if p["MONOTONE"] > 0.5:
        # Level detector.  MONO_GAIN is set so the falsified model reaches the
        # SAME peak drive as the base model on an untreated cycle, so the
        # comparison cannot be rigged by a difference of scale.
        ns_drive = p["MONO_GAIN"] * load / p["KP"]
    else:
        r = load / p["KP"]
        ns_drive = sens * r * math.exp(1.0 - r)
    ns_drive *= (1.0 - p["CA_SUPP"])
    h = p["NS_HILL"]
    eff = ns_drive ** h / (ns_drive ** h + p["NS50"] ** h)

    # ───────────────────────────────────────────────── serotonin ────────
    occ = p["OCC_MAX"] * ssri_eq / (ssri_eq + p["EC50_SERT"])
    brake = p["WB"] * max(0.0, SHT - 1.0) * (1.0 - DES)
    dSHT = p["KS_SYN"] / (1.0 + brake) - p["KS_OUT"] * (1.0 - occ) * SHT \
        - p["K5HT_NS"] * eff * SHT
    dDES = (min(1.0, p["KD2"] * max(0.0, SHT - 1.0)) - DES) / p["TAU_DES"]
    dBDNF = (p["KBDNF"] * max(0.0, SHT - 1.0) - BDNF) / p["TAU_BDNF"]

    # ─────────────────────────────────────────────────── HPA axis ───────
    allo_brake = p["K_ALLO_BRAKE"] * (ALLOB / (ALLOB + 2.0)) / sens
    dCRH = ((p["STRESS"] + 0.7) * (1.0 - allo_brake)
            / (1.0 + 0.8 * (CORT - 1.0)) - CRH) / p["TAU_CRH"]
    dACTH = (CRH - ACTH) / p["TAU_ACTH"]
    dCORT = (ACTH - CORT) / p["TAU_CORT"]

    # ────────────────────────────────────────────── corticolimbic ───────
    bz_occ = ALPZ_C / (ALPZ_C + 8.0) / (1.0 + 0.60 * (ALPHA4 - 1.0))
    prog_ex = DRSP_C / (DRSP_C + p["KI_PR"])
    amy_ss = 1.0 + p["W_NS_AMY"] * eff + p["W_CORT_AMY"] * (CORT / p["CORT0"] - 1.0) \
        - p["W_SHT_AMY"] * max(0.0, SHT - 1.0) - p["W_PFC_AMY"] * (PFC - 1.0) \
        - p["W_BZ_AMY"] * bz_occ + p["W_PROG_AMY"] * prog_ex
    pfc_ss = 1.0 - p["P_NS_PFC"] * eff \
        + p["P_E2_PFC"] * (E2_TOT / p["E2REF"] - 1.0) + p["CBT"]
    rew_ss = 1.0 - p["R_NS_REW"] * eff + p["R_BDNF"] * BDNF - p["R_WE2"] * WE2
    dAMY = (max(amy_ss, 0.2) - AMY) / p["TAU_AMY"]
    dPFC = (max(pfc_ss, 0.2) - PFC) / p["TAU_PFC"]
    dREW = (max(rew_ss, 0.2) - REW) / p["TAU_REW"]

    # ──────────────────────────────────────────────────── somatic ───────
    dALDO = (1.0 + p["A_P4"] * max(P4, 0.0) / p["P4REF"] - ALDO) / p["TAU_ALDO"]
    mr_block = 1.0 / (1.0 + DRSP_C / p["KI_MR_DRSP"] + p["SPIRO"] / p["KI_MR_SPIRO"]
                      + max(P4, 0.0) / p["KI_MR_P4"])
    dECFV = (1.0 + p["ECF_GAIN"] * ALDO * mr_block - ECFV) / p["TAU_ECFV"]
    pr_block = 1.0 / (1.0 + UPA_C / p["UPA_PR_KI"])
    br_ss = p["B_E2"] * E2_TOT / p["E2REF"] \
        + p["B_P4"] * (max(P4, 0.0) / p["P4REF"] + prog_ex) * pr_block
    dBREAST = (br_ss - BREAST) / p["TAU_BR"]
    dPGE = (p["KPG"] * pg_window(cd) * (1.0 - p["NSAID"])
            * (0.0 if p["OOPH"] > 0.5 else 1.0) - PGE) / p["TAU_PG"]

    # ───────────────────────────────────── DRSP symptom domains ─────────
    ca = 1.0 - (0.25 if p["CA_SUPP"] > 0.0 else 0.0)
    bloat = max(0.0, (ECFV - 1.0) / 0.030)
    cramp = max(0.0, PGE)
    sleep_q = 1.0 - clamp01(0.35 * eff + 0.030 * HF)
    d_amy = max(0.0, AMY - 1.0)
    d_pfc = max(0.0, 1.0 - PFC)
    d_rew = max(0.0, 1.0 - REW)
    sht_up = max(0.0, SHT - 1.0)
    sht_dn = max(0.0, 1.0 - SHT)

    f = {
        "S_IRR": 0.62 * d_amy + 0.22 * eff - 0.34 * sht_up - 0.18 * (PFC - 1.0),
        "S_ANX": 0.66 * d_amy + 0.18 * eff - 0.28 * sht_up,
        "S_LAB": 0.55 * d_amy + 0.30 * d_pfc + 0.15 * WE2 - 0.22 * sht_up,
        "S_DEP": 0.42 * eff + 0.38 * d_rew + 0.18 * WE2 - 0.30 * BDNF - 0.15 * sht_up,
        "S_ANH": 0.28 * eff + 0.62 * d_rew - 0.32 * BDNF,
        "S_COG": 0.52 * d_pfc + 0.26 * eff + 0.14 * (1.0 - sleep_q),
        "S_FAT": 0.34 * eff + 0.32 * (1.0 - sleep_q) + 0.22 * d_rew,
        "S_APP": (0.58 * eff + 0.30 * sht_dn - 0.20 * sht_up) * ca,
        "S_SLP": (1.0 - sleep_q) * 1.15,
        "S_OVR": 0.58 * eff + 0.30 * d_pfc - 0.20 * sht_up,
        "S_PHY": (0.34 * bloat + 0.34 * BREAST + 0.32 * cramp) * ca,
    }
    cur = dict(S_IRR=S_IRR, S_DEP=S_DEP, S_ANX=S_ANX, S_LAB=S_LAB, S_ANH=S_ANH,
               S_COG=S_COG, S_FAT=S_FAT, S_APP=S_APP, S_SLP=S_SLP, S_OVR=S_OVR,
               S_PHY=S_PHY)
    dS = {}
    for k, v in f.items():
        dS[k] = (1.0 + 5.0 * clamp01(v + p["BASE_SYM"]) - cur[k]) / p["TAU_S"]

    # ────────────────────────────────────────── organ / safety ──────────
    hf_ss = p["HF_MAX"] / (1.0 + (E2_TOT / p["E50_HF"]) ** 2) \
        * (1.0 - p["HF_SSRI"] * occ)
    dHF = (hf_ss - HF) / p["TAU_HF"]
    hypo = 1.0 / (1.0 + (E2_TOT / p["E50_BMD"]) ** 2)
    dBMD = -p["KBMD_LOSS"] * max(0.0, hypo - p["HYPO_REF"]) \
        + p["KBMD_REC"] * max(0.0, -BMD) \
        * max(0.0, p["HYPO_REF"] - hypo) / p["HYPO_REF"]
    dENDO = ((DRSP_C / (DRSP_C + 5.0)) - ENDO) / p["TAU_ENDO"]

    out[0] = dFSH; out[1] = dLH; out[2] = dFOLL; out[3] = dCL
    out[4] = dE2; out[5] = dP4
    out[6] = dDHP; out[7] = dALLOP; out[8] = dALLOB; out[9] = dISOB; out[10] = dZUR
    out[11] = dMALLO; out[12] = dMPROG; out[13] = dDELTA; out[14] = dALPHA4
    out[15] = dME2
    out[16] = dSHT; out[17] = dDES; out[18] = dBDNF
    out[19] = dAMY; out[20] = dPFC; out[21] = dREW
    out[22] = dCRH; out[23] = dACTH; out[24] = dCORT
    out[25] = dALDO; out[26] = dECFV; out[27] = dBREAST; out[28] = dPGE
    out[29] = dS["S_IRR"]; out[30] = dS["S_DEP"]; out[31] = dS["S_ANX"]
    out[32] = dS["S_LAB"]; out[33] = dS["S_ANH"]; out[34] = dS["S_COG"]
    out[35] = dS["S_FAT"]; out[36] = dS["S_APP"]; out[37] = dS["S_SLP"]
    out[38] = dS["S_OVR"]; out[39] = dS["S_PHY"]
    out[40] = dHF; out[41] = dBMD; out[42] = dENDO
    out[43] = dSERT_D; out[44] = dSERT_C; out[45] = dDMS_C
    out[46] = dFLX_C; out[47] = dNFLX_C
    out[48] = dDRSP_D; out[49] = dDRSP_C
    out[50] = dLEUP_D; out[51] = dLEUP_C
    out[52] = dSEPRA_D; out[53] = dSEPRA_C
    out[54] = dDUTA_D; out[55] = dDUTA_C
    out[56] = dALPZ_D; out[57] = dALPZ_C
    out[58] = dUPA_C; out[59] = dPITDS
    return out


# ══════════════════════════════════════════════════════════════════════════
# 6. INITIAL CONDITIONS + INTEGRATOR
# ══════════════════════════════════════════════════════════════════════════
def y0(p):
    y = [0.0] * NST
    y[IX["FSH"]] = 3.0
    y[IX["LH"]] = 1.8
    y[IX["FOLL"]] = 0.05
    y[IX["E2"]] = 40.0
    y[IX["P4"]] = 0.3
    y[IX["ALLOB"]] = 1.3
    y[IX["MALLO"]] = 1.3
    y[IX["MPROG"]] = 0.03
    y[IX["DELTA"]] = 1.4
    y[IX["ALPHA4"]] = 1.1
    y[IX["ME2"]] = 40.0
    for s in ("SHT", "AMY", "PFC", "REW", "CRH", "ACTH", "CORT", "ALDO", "ECFV"):
        y[IX[s]] = 1.0
    for s in ALL_DOMAINS:
        y[IX[s]] = 1.5
    return y


RECORD_KEYS = ["E2", "P4", "ALLOP", "ALLOB", "ISOB", "ZUR", "FOLL", "CL",
               "FSH", "LH", "DELTA", "ALPHA4", "SHT", "AMY", "PFC", "REW",
               "CORT", "ECFV", "BREAST", "HF", "BMD",
               "SERT_C", "DRSP_C", "DUTA_C", "LEUP_C", "SEPRA_C", "UPA_C",
               "ALPZ_C"] + list(ALL_DOMAINS)


def simulate(p, regimen, tend, dt=0.02, record_every=0.25):
    y = y0(p)
    ev = regimen.events if regimen else []
    ei = 0
    k1 = [0.0] * NST; k2 = [0.0] * NST; k3 = [0.0] * NST; k4 = [0.0] * NST
    ytmp = [0.0] * NST
    rec = {"t": [], "cd": [], "DRSP": [], "CORE": []}
    for kk in RECORD_KEYS:
        rec[kk] = []
    nrec = max(1, int(round(record_every / dt)))
    nstep = int(round(tend / dt))

    for i in range(nstep + 1):
        t = i * dt
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            y[ev[ei][1]] += ev[ei][2]
            ei += 1
        if i % nrec == 0:
            rec["t"].append(t)
            rec["cd"].append(math.fmod(t, p["TCYC"]) + 1.0)
            for kk in RECORD_KEYS:
                rec[kk].append(y[IX[kk]])
            rec["DRSP"].append(sum(y[IX[s]] for s in ALL_DOMAINS))
            rec["CORE"].append(sum(y[IX[s]] for s in CORE_DOMAINS))
        if i == nstep:
            break
        derivatives(t, y, p, k1)
        for j in range(NST):
            ytmp[j] = y[j] + 0.5 * dt * k1[j]
        derivatives(t + 0.5 * dt, ytmp, p, k2)
        for j in range(NST):
            ytmp[j] = y[j] + 0.5 * dt * k2[j]
        derivatives(t + 0.5 * dt, ytmp, p, k3)
        for j in range(NST):
            ytmp[j] = y[j] + dt * k3[j]
        derivatives(t + dt, ytmp, p, k4)
        for j in range(NST):
            y[j] += dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
        for j in range(NST):
            if j != IX["BMD"] and y[j] < 0.0:
                y[j] = 0.0
    return rec


# ══════════════════════════════════════════════════════════════════════════
# 7. SUMMARY METRICS
# ══════════════════════════════════════════════════════════════════════════
def cycle_window(rec, cyc, lo, hi, key, tcyc=28.0):
    out = []
    for i, t in enumerate(rec["t"]):
        if int(t // tcyc) == cyc and lo <= rec["cd"][i] < hi:
            out.append(rec[key][i])
    return out


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def metrics(rec, cyc=4, tcyc=28.0):
    """Summarise one cycle.

    Follicular window = cycle days 5-11 (symptom-free by definition of PMDD).
    Luteal window     = cycle days 21-28 (the trial-relevant window).
    The DRSP peak is located on a LUTEAL-ALIGNED axis (day 8 of `cyc` through
    day 6 of the next cycle, reported as days 29-34) so that a perimenstrual
    peak is not mis-dated to the beginning of the cycle.
    """
    m = {}

    def win(lo, hi, key):
        return cycle_window(rec, cyc, lo, hi, key, tcyc)

    m["fol"] = mean(win(5.0, 11.0, "DRSP"))
    m["lut_mean"] = mean(win(21.0, 29.0, "DRSP"))
    m["core_fol"] = mean(win(5.0, 11.0, "CORE"))
    m["core_lut"] = mean(win(21.0, 29.0, "CORE"))
    m["pct_inc"] = 100.0 * (m["core_lut"] - m["core_fol"]) / m["core_fol"]

    al_d, al_v, al_c = [], [], []
    for i, t in enumerate(rec["t"]):
        c, cdi = int(t // tcyc), rec["cd"][i]
        if c == cyc and cdi >= 8.0:
            al_d.append(cdi); al_v.append(rec["DRSP"][i]); al_c.append(rec["CORE"][i])
        elif c == cyc + 1 and cdi < 6.0:
            al_d.append(cdi + tcyc); al_v.append(rec["DRSP"][i])
            al_c.append(rec["CORE"][i])
    m["lut_peak"] = max(al_v) if al_v else float("nan")
    m["drsp_peak_day"] = al_d[al_v.index(max(al_v))] if al_v else float("nan")
    m["core_peak"] = max(al_c) if al_c else float("nan")
    m["amp"] = m["lut_peak"] - m["fol"]

    allo = win(1.0, 29.0, "ALLOB")
    cds = win(1.0, 29.0, "cd")
    m["allo_peak"] = max(allo) if allo else float("nan")
    m["allo_peak_day"] = cds[allo.index(max(allo))] if allo else float("nan")
    m["allo_fol"] = mean(win(5.0, 11.0, "ALLOB"))
    m["allop_peak"] = max(win(1.0, 29.0, "ALLOP") or [float("nan")])
    m["lag"] = m["drsp_peak_day"] - m["allo_peak_day"]

    e2, p4 = win(1.0, 29.0, "E2"), win(1.0, 29.0, "P4")
    m["e2_peak"] = max(e2) if e2 else float("nan")
    m["e2_mean"] = mean(e2)
    m["p4_peak"] = max(p4) if p4 else float("nan")
    m["p4_peak_day"] = cds[p4.index(max(p4))] if p4 else float("nan")
    m["hf"] = mean(win(1.0, 29.0, "HF"))
    # BMD change accrued during THIS cycle only (so that a run with a long
    # untreated lead-in is not charged for losses that predate the treatment)
    bmd_c = win(1.0, 29.0, "BMD")
    m["bmd_cyc"] = (bmd_c[-1] - bmd_c[0]) if bmd_c else float("nan")
    m["a4_peak"] = max(win(1.0, 29.0, "ALPHA4") or [float("nan")])
    m["delta_peak"] = max(win(1.0, 29.0, "DELTA") or [float("nan")])
    m["bmd"] = rec["BMD"][-1]
    return m


def core_reduction(ref, m):
    """% reduction of the luteal CORE-affective burden above the scale floor."""
    base = ref["core_lut"] - CORE_FLOOR
    if base <= 0:
        return float("nan")
    return 100.0 * (1.0 - (m["core_lut"] - CORE_FLOOR) / base)


def window_mean(rec, t0, t1, key="DRSP"):
    """Mean of `key` over absolute time [t0, t1) — used for the add-back flare."""
    v = [rec[key][i] for i, t in enumerate(rec["t"]) if t0 <= t < t1]
    return mean(v)


def luteal_reduction(ref, m):
    """% reduction of the luteal DRSP burden above the floor of the scale.

    This is the model analogue of the primary endpoint of PMDD trials (change
    in the luteal-phase DRSP score).  The floor is subtracted because a DRSP-11
    of 11 means "no symptoms at all", not "zero".
    """
    base = ref["lut_mean"] - DRSP_FLOOR
    if base <= 0:
        return float("nan")
    return 100.0 * (1.0 - (m["lut_mean"] - DRSP_FLOOR) / base)


# ══════════════════════════════════════════════════════════════════════════
# 8. SCENARIOS
# ══════════════════════════════════════════════════════════════════════════
NCYC = 6
TEND = NCYC * 28.0
LONG = 10 * 28.0          # dutasteride: t1/2 ~35 d needs a long lead-in
ADDBACK_END = 10 * 28.0   # 4 suppressed cycles, then 6 cycles of add-back


def base_params(control=False, **over):
    p = dict(P)
    if control:
        p.update(CONTROL_OVERRIDES)
    p.update(over)
    return p


def scenarios():
    """(name, label, params, regimen, t_end, cycle_to_evaluate)."""
    S = []
    add = S.append

    add(("pmdd_untreated", "Untreated PMDD (reference)",
         base_params(), Regimen().finalize(), TEND, 4))
    add(("control", "Asymptomatic control — identical hormones",
         base_params(control=True), Regimen().finalize(), TEND, 4))

    # ── SSRIs ────────────────────────────────────────────────────────────
    add(("sert50_luteal", "Sertraline 50 mg, luteal-only (cd 15-28)",
         base_params(),
         Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0, (15.0, 29.0)).finalize(),
         TEND, 4))
    add(("sert_cycle1", "Sertraline 50 mg luteal — FIRST treated cycle",
         base_params(),
         Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0, (15.0, 29.0)).finalize(),
         TEND, 0))
    add(("sert50_cont", "Sertraline 50 mg continuous",
         base_params(), Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0).finalize(),
         TEND, 4))
    add(("sert100_cont", "Sertraline 100 mg continuous",
         base_params(), Regimen().add("SERT_D", 100.0, 0.0, TEND, 1.0).finalize(),
         TEND, 4))
    add(("sert50_onset", "Sertraline 50 mg, symptom-onset dosing (cd 22-28)",
         base_params(),
         Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0, (22.0, 29.0)).finalize(),
         TEND, 4))
    add(("flx20_cont", "Fluoxetine 20 mg continuous",
         base_params(),
         Regimen().add("FLX_C", 20.0 * 1000.0 * 0.9 / 2500.0, 0.0, TEND, 1.0).finalize(),
         TEND, 4))

    # ── ovulation suppression ────────────────────────────────────────────
    add(("drsp_24_4", "Drospirenone 3 mg / EE 20 ug — 24/4",
         base_params(),
         Regimen().add("DRSP_D", 3.0, 0.0, TEND, 1.0, (1.0, 25.0)).finalize(),
         TEND, 4))
    add(("drsp_21_7", "Drospirenone 3 mg / EE 20 ug — 21/7",
         base_params(),
         Regimen().add("DRSP_D", 3.0, 0.0, TEND, 1.0, (1.0, 22.0)).finalize(),
         TEND, 4))
    add(("drsp_cont", "Drospirenone / EE — continuous, no HFI",
         base_params(), Regimen().add("DRSP_D", 3.0, 0.0, TEND, 1.0).finalize(),
         TEND, 4))

    # ── GnRH agonist ± add-back ──────────────────────────────────────────
    add(("leuprolide", "Leuprolide 3.75 mg IM q28d",
         base_params(), Regimen().add("LEUP_D", 3.75, 0.0, TEND, 28.0).finalize(),
         TEND, 4))
    ab_reg = Regimen().add("LEUP_D", 3.75, 0.0, ADDBACK_END, 28.0).finalize()
    ab_par = dict(ADDBACK_E2=210.0, ADDBACK_P4=20.0, ADDBACK_T0=112.0)
    add(("addback_early", "Leuprolide + E2/P4 add-back from d112 — first cycle",
         base_params(**ab_par), ab_reg, ADDBACK_END, 4))
    add(("addback_late", "the SAME run 4 cycles later (identical exposure)",
         base_params(**ab_par), ab_reg, ADDBACK_END, 8))
    add(("addback_ctrl", "Leuprolide + add-back in a CONTROL — first cycle",
         base_params(control=True, **ab_par), ab_reg, ADDBACK_END, 4))

    # ── neurosteroid-directed ────────────────────────────────────────────
    add(("sepranolone", "Sepranolone 10 mg SC q48h, luteal",
         base_params(),
         Regimen().add("SEPRA_D", 10.0, 0.0, TEND, 2.0, (15.0, 29.0)).finalize(),
         TEND, 4))
    add(("duta_2p5", "Dutasteride 2.5 mg/d",
         base_params(), Regimen().add("DUTA_D", 2.5, 0.0, LONG, 1.0).finalize(),
         LONG, 8))
    add(("duta_0p5", "Dutasteride 0.5 mg/d",
         base_params(), Regimen().add("DUTA_D", 0.5, 0.0, LONG, 1.0).finalize(),
         LONG, 8))
    add(("allo_high", "Exogenous neurosteroid, SEDATIVE dose (20 nM/d, cd 22-28)",
         base_params(),
         Regimen().add("ZUR", 20.0, 0.0, TEND, 1.0, (22.0, 29.0)).finalize(),
         TEND, 4))
    add(("allo_low", "the SAME agent, SUB-sedative dose (3 nM/d, cd 22-28)",
         base_params(),
         Regimen().add("ZUR", 3.0, 0.0, TEND, 1.0, (22.0, 29.0)).finalize(),
         TEND, 4))
    add(("alprazolam", "Alprazolam 0.25 mg TID, luteal",
         base_params(),
         Regimen().add("ALPZ_D", 0.25, 0.0, TEND, 1.0 / 3.0, (15.0, 29.0)).finalize(),
         TEND, 4))
    add(("ulipristal", "Ulipristal acetate 5 mg/d",
         base_params(),
         Regimen().add("UPA_C", 5.0 * 1000.0 * 0.8 / 500.0, 0.0, TEND, 1.0).finalize(),
         TEND, 4))

    # ── adjuncts, combination, surgery ───────────────────────────────────
    add(("calcium", "Calcium carbonate 1200 mg/d",
         base_params(CA_SUPP=0.18), Regimen().finalize(), TEND, 4))
    add(("spironolactone", "Spironolactone 100 mg/d luteal (somatic)",
         base_params(SPIRO=1.6), Regimen().finalize(), TEND, 4))
    add(("cbt", "CBT / mindfulness (prefrontal training)",
         base_params(CBT=0.18), Regimen().finalize(), TEND, 4))
    add(("sert_plus_drsp", "Sertraline 50 mg luteal + drospirenone 24/4",
         base_params(),
         Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0, (15.0, 29.0))
                  .add("DRSP_D", 3.0, 0.0, TEND, 1.0, (1.0, 25.0)).finalize(),
         TEND, 4))
    add(("oophorectomy", "Bilateral oophorectomy + transdermal E2",
         base_params(OOPH=1.0, ADDBACK_E2=170.0), Regimen().finalize(), TEND, 4))

    # ── falsification runs ───────────────────────────────────────────────
    fal = dict(MONOTONE=1.0, GAMMA_D=0.0, GAMMA_A4=0.0, RATE_OFF=1.0)
    add(("fals_untreated", "FALSIFY: level detector, no plasticity",
         base_params(**fal), Regimen().finalize(), TEND, 4))
    add(("fals_sert", "FALSIFY + sertraline 50 mg luteal",
         base_params(**fal),
         Regimen().add("SERT_D", 50.0, 0.0, TEND, 1.0, (15.0, 29.0)).finalize(),
         TEND, 4))
    add(("fals_duta_2p5", "FALSIFY + dutasteride 2.5 mg",
         base_params(**fal), Regimen().add("DUTA_D", 2.5, 0.0, LONG, 1.0).finalize(),
         LONG, 8))
    add(("fals_duta_0p5", "FALSIFY + dutasteride 0.5 mg",
         base_params(**fal), Regimen().add("DUTA_D", 0.5, 0.0, LONG, 1.0).finalize(),
         LONG, 8))
    add(("fals_24_4", "FALSIFY + drospirenone 24/4",
         base_params(**fal),
         Regimen().add("DRSP_D", 3.0, 0.0, TEND, 1.0, (1.0, 25.0)).finalize(),
         TEND, 4))
    add(("fals_21_7", "FALSIFY + drospirenone 21/7",
         base_params(**fal),
         Regimen().add("DRSP_D", 3.0, 0.0, TEND, 1.0, (1.0, 22.0)).finalize(),
         TEND, 4))
    add(("fals_allo_high", "FALSIFY + sedative-dose neurosteroid",
         base_params(**fal),
         Regimen().add("ZUR", 20.0, 0.0, TEND, 1.0, (22.0, 29.0)).finalize(),
         TEND, 4))
    add(("fals_allo_low", "FALSIFY + sub-sedative neurosteroid",
         base_params(**fal),
         Regimen().add("ZUR", 3.0, 0.0, TEND, 1.0, (22.0, 29.0)).finalize(),
         TEND, 4))
    add(("fals_addback_early", "FALSIFY + add-back, first cycle",
         base_params(**dict(ab_par, **fal)), ab_reg, ADDBACK_END, 4))
    add(("fals_addback_late", "FALSIFY + add-back, 4 cycles later",
         base_params(**dict(ab_par, **fal)), ab_reg, ADDBACK_END, 8))
    return S


# ══════════════════════════════════════════════════════════════════════════
# 9. MAIN
# ══════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=None, help="directory for per-scenario CSVs")
    ap.add_argument("--dt", type=float, default=0.02)
    args = ap.parse_args()

    res, lab, recs = {}, {}, {}
    print("=" * 79)
    print("PMDD QSP model — Python twin   (%d ODE states, RK4, dt = %.3f d)"
          % (NST, args.dt))
    print("=" * 79)
    print("  %-18s %8s %8s %8s %8s" %
          ("scenario", "DRSPfol", "DRSPlut", "peak", "cd(peak)"))
    for name, label, p, reg, tend, cyc in scenarios():
        rec = simulate(p, reg, tend, dt=args.dt)
        m = metrics(rec, cyc=cyc)
        res[name], lab[name], recs[name] = m, label, rec
        print("  %-18s %8.1f %8.1f %8.1f %8.1f   %s"
              % (name, m["fol"], m["lut_mean"], m["lut_peak"],
                 m["drsp_peak_day"], label[:36]))
        if args.csv:
            os.makedirs(args.csv, exist_ok=True)
            cols = ["t", "cd", "E2", "P4", "ALLOB", "DELTA", "ALPHA4", "SHT",
                    "AMY", "PFC", "REW", "HF", "BMD", "DRSP", "CORE"]
            with open(os.path.join(args.csv, name + ".csv"), "w") as fh:
                fh.write(",".join(cols) + "\n")
                for i in range(len(rec["t"])):
                    fh.write(",".join("%.5g" % rec[c][i] for c in cols) + "\n")

    ref, fal = res["pmdd_untreated"], res["fals_untreated"]

    def red(m, r=None):
        return luteal_reduction(r or ref, m)

    print("\n" + "-" * 79)
    print("PHYSIOLOGY OF THE UNTREATED CYCLE (identical in PMDD and controls)")
    print("-" * 79)
    print("  E2 peak                        %6.1f pg/mL  (target 180-320)" % ref["e2_peak"])
    print("  P4 peak                        %6.1f ng/mL  (target 8-16)" % ref["p4_peak"])
    print("  P4 peak cycle day              %6.1f        (target 20-23)" % ref["p4_peak_day"])
    print("  plasma ALLO luteal peak        %6.2f nM     (target 4-8)" % ref["allop_peak"])
    print("  brain ALLO follicular          %6.2f nM     (target 1-2)" % ref["allo_fol"])
    print("  brain ALLO luteal peak         %6.2f nM     (target 4-8)" % ref["allo_peak"])
    print("  brain ALLO peak cycle day      %6.1f        (target 20-24)" % ref["allo_peak_day"])
    print("  alpha4 peak / delta peak       %6.2f / %.2f x baseline"
          % (ref["a4_peak"], ref["delta_peak"]))
    print("\n  [P1] DRSP peak cycle day       %6.1f        (target 25-28)" % ref["drsp_peak_day"])
    print("  [P1] lag behind the ALLO peak  %6.1f d      (target 3-6 d)" % ref["lag"])
    print("  PMDD core luteal increase      %+6.0f %%      (DSM-5 needs >=30%%)"
          % ref["pct_inc"])
    print("  control core luteal increase   %+6.0f %%      (must stay <30%%)"
          % res["control"]["pct_inc"])
    print("  PMDD    DRSP-11 follicular -> luteal peak   %.1f -> %.1f"
          % (ref["fol"], ref["lut_peak"]))
    print("  control DRSP-11 follicular -> luteal peak   %.1f -> %.1f"
          % (res["control"]["fol"], res["control"]["lut_peak"]))

    print("\n" + "-" * 79)
    print("TREATMENT EFFECT — %% reduction of the luteal DRSP burden above floor")
    print("(the model analogue of the primary endpoint of PMDD trials)")
    print("-" * 79)
    for k in ["sert50_luteal", "sert_cycle1", "sert50_cont", "sert100_cont",
              "sert50_onset", "flx20_cont", "drsp_24_4", "drsp_21_7",
              "drsp_cont", "leuprolide", "addback_early", "addback_late",
              "sepranolone", "duta_2p5", "duta_0p5", "allo_high", "allo_low",
              "alprazolam", "ulipristal", "calcium", "spironolactone", "cbt",
              "sert_plus_drsp", "oophorectomy"]:
        print("  %-15s %+7.1f %%   %s" % (k, red(res[k]), lab[k][:45]))

    print("\n" + "-" * 79)
    print("THE FIVE PREDICTIONS vs THE FALSIFIED (LEVEL-DETECTOR) MODEL")
    print("-" * 79)
    print("  [P1] symptom-peak lag          %+6.1f d   ->  %+6.1f d"
          % (ref["lag"], fal["lag"]))
    print("  [P2] dutasteride 2.5 mg        %+6.1f %%   ->  %+6.1f %%"
          % (red(res["duta_2p5"]), red(res["fals_duta_2p5"], fal)))
    print("  [P2] dutasteride 0.5 mg        %+6.1f %%   ->  %+6.1f %%"
          % (red(res["duta_0p5"]), red(res["fals_duta_0p5"], fal)))
    print("  [P3] neurosteroid, sedative    %+6.1f %%   ->  %+6.1f %%"
          % (red(res["allo_high"]), red(res["fals_allo_high"], fal)))
    print("  [P3] neurosteroid, sub-dose    %+6.1f %%   ->  %+6.1f %%"
          % (red(res["allo_low"]), red(res["fals_allo_low"], fal)))
    print("  [P4] 24/4 minus 21/7           %+6.1f pt  ->  %+6.1f pt"
          % (red(res["drsp_24_4"]) - red(res["drsp_21_7"]),
             red(res["fals_24_4"], fal) - red(res["fals_21_7"], fal)))
    # [P5] is a statement about the TRANSITION, so it is measured on absolute
    # time: the first 14 days of add-back against days 112-126 of add-back,
    # under identical, unchanging hormone input.
    T0 = 112.0
    ab_flare = window_mean(recs["addback_early"], T0, T0 + 10.0)
    ab_adapt = window_mean(recs["addback_early"], T0 + 112.0, T0 + 126.0)
    ab_supp = window_mean(recs["addback_early"], T0 - 28.0, T0)
    fb_flare = window_mean(recs["fals_addback_early"], T0, T0 + 10.0)
    fb_adapt = window_mean(recs["fals_addback_early"], T0 + 112.0, T0 + 126.0)
    fb_supp = window_mean(recs["fals_addback_early"], T0 - 28.0, T0)
    print("  [P5] DRSP-11 on suppression      %6.1f    ->  %6.1f" % (ab_supp, fb_supp))
    print("  [P5] first 10 d of add-back      %6.1f    ->  %6.1f" % (ab_flare, fb_flare))
    print("  [P5] day 112-126 of add-back     %6.1f    ->  %6.1f" % (ab_adapt, fb_adapt))
    print("       PMDD: the flare FADES (%+.1f -> %+.1f above suppression);"
          % (ab_flare - ab_supp, ab_adapt - ab_supp))
    print("       level detector: it never fades (%+.1f -> %+.1f)"
          % (fb_flare - fb_supp, fb_adapt - fb_supp))
    print("  SSRI raises brain ALLO         %.2f -> %.2f nM (3alpha-HSD shift)"
          % (ref["allo_peak"], res["sert50_luteal"]["allo_peak"]))

    print("\n" + "-" * 79)
    print("SAFETY / ORGAN-SYSTEM READOUTS")
    print("-" * 79)
    print("  hot flushes/day  untreated %4.2f | leuprolide %4.2f | add-back %4.2f"
          % (ref["hf"], res["leuprolide"]["hf"], res["addback_late"]["hf"]))
    print("  mean E2 (pg/mL)  untreated %4.0f | leuprolide %4.0f | 24/4 %4.0f"
          % (ref["e2_mean"], res["leuprolide"]["e2_mean"], res["drsp_24_4"]["e2_mean"]))
    print("  lumbar BMD change: leuprolide %+5.2f %% over %d cycles | "
          "add-back %+5.2f %% | oophorectomy+E2 %+5.2f %%"
          % (res["leuprolide"]["bmd"], NCYC, res["addback_late"]["bmd"],
             res["oophorectomy"]["bmd"]))

    # ─────────────────────────────────────────────────────────── checks ──
    print("\n" + "=" * 79)
    print("SELF-CHECKS — every one of these is a claim made in README.md")
    print("=" * 79)
    checks = []

    def chk(name, cond, detail):
        checks.append((name, bool(cond), detail))

    chk("cycle: E2 peak 180-320 pg/mL", 180 <= ref["e2_peak"] <= 320,
        "%.0f" % ref["e2_peak"])
    chk("cycle: P4 peak 8-16 ng/mL on day 20-23",
        8 <= ref["p4_peak"] <= 16 and 20 <= ref["p4_peak_day"] <= 23,
        "%.1f ng/mL, day %.0f" % (ref["p4_peak"], ref["p4_peak_day"]))
    chk("cycle: plasma ALLO luteal peak 4-8 nM", 4 <= ref["allop_peak"] <= 8,
        "%.2f" % ref["allop_peak"])
    chk("cycle: brain ALLO 1-2 nM follicular and 4-8 nM luteal",
        1.0 <= ref["allo_fol"] <= 2.0 and 4 <= ref["allo_peak"] <= 8,
        "%.2f -> %.2f" % (ref["allo_fol"], ref["allo_peak"]))
    chk("[P1] symptom peak on cycle day 25-28", 25 <= ref["drsp_peak_day"] <= 28.5,
        "day %.1f" % ref["drsp_peak_day"])
    chk("[P1] symptom peak lags the ALLO peak by 3-6 d", 3.0 <= ref["lag"] <= 6.0,
        "%.1f d" % ref["lag"])
    chk("PMDD core luteal increase >=30% (DSM-5 criterion)", ref["pct_inc"] >= 30,
        "%.0f%%" % ref["pct_inc"])
    chk("control core luteal increase <30% on the SAME hormones",
        res["control"]["pct_inc"] < 30, "%.0f%%" % res["control"]["pct_inc"])
    chk("sertraline 50 mg luteal-only cuts the luteal burden 35-65%",
        35 <= red(res["sert50_luteal"]) <= 65, "%.1f%%" % red(res["sert50_luteal"]))
    chk("the SSRI benefit is already there in the FIRST treated cycle",
        red(res["sert_cycle1"]) >= 0.70 * red(res["sert50_luteal"]),
        "%.1f%% vs %.1f%% at steady state"
        % (red(res["sert_cycle1"]), red(res["sert50_luteal"])))
    chk("continuous SSRI beats luteal-only on the core affective domains",
        res["sert50_cont"]["core_lut"] < res["sert50_luteal"]["core_lut"],
        "core luteal %.1f vs %.1f"
        % (res["sert50_cont"]["core_lut"], res["sert50_luteal"]["core_lut"]))
    chk("drospirenone 24/4 cuts the luteal burden 30-70%",
        30 <= red(res["drsp_24_4"]) <= 70, "%.1f%%" % red(res["drsp_24_4"]))
    chk("[P4] 24/4 beats 21/7 by >=8 points",
        red(res["drsp_24_4"]) - red(res["drsp_21_7"]) >= 8.0,
        "%.1f pt" % (red(res["drsp_24_4"]) - red(res["drsp_21_7"])))
    chk("continuous dosing (no HFI) beats 24/4",
        red(res["drsp_cont"]) > red(res["drsp_24_4"]),
        "%.1f%% vs %.1f%%" % (red(res["drsp_cont"]), red(res["drsp_24_4"])))
    chk("leuprolide cuts the luteal burden >=60% (DRSP-11)",
        red(res["leuprolide"]) >= 60, "%.1f%%" % red(res["leuprolide"]))
    chk("leuprolide cuts the CORE affective burden >=78%",
        core_reduction(ref, res["leuprolide"]) >= 78,
        "%.1f%%" % core_reduction(ref, res["leuprolide"]))
    chk("[P5] add-back reinstates symptoms in PMDD (>=8 DRSP-11 points)",
        ab_flare - ab_supp >= 8.0, "%+.1f points above suppression"
        % (ab_flare - ab_supp))
    chk("[P5] the flare then FADES under identical continued exposure",
        (ab_adapt - ab_supp) < 0.4 * (ab_flare - ab_supp),
        "%+.1f -> %+.1f points" % (ab_flare - ab_supp, ab_adapt - ab_supp))
    chk("[P5] add-back does NOT make a control symptomatic",
        res["addback_ctrl"]["pct_inc"] < 30, "%.0f%%" % res["addback_ctrl"]["pct_inc"])
    chk("[P2] dutasteride 2.5 mg cuts the luteal burden >=25%",
        red(res["duta_2p5"]) >= 25, "%.1f%%" % red(res["duta_2p5"]))
    chk("[P2] dutasteride 0.5 mg does almost nothing (<12%)",
        red(res["duta_0p5"]) < 12, "%.1f%%" % red(res["duta_0p5"]))
    chk("[P3] a SEDATIVE-dose neurosteroid HELPS (>=25%)",
        red(res["allo_high"]) >= 25, "%.1f%%" % red(res["allo_high"]))
    chk("[P3] a SUB-sedative dose of the same agent HURTS (<0%)",
        red(res["allo_low"]) < 0.0, "%.1f%%" % red(res["allo_low"]))
    chk("sepranolone cuts the luteal burden 20-55%",
        20 <= red(res["sepranolone"]) <= 55, "%.1f%%" % red(res["sepranolone"]))
    chk("alprazolam gives a modest benefit (5-35%)",
        5 <= red(res["alprazolam"]) <= 35, "%.1f%%" % red(res["alprazolam"]))
    chk("calcium 1200 mg gives a small benefit (5-30%)",
        5 <= red(res["calcium"]) <= 30, "%.1f%%" % red(res["calcium"]))
    chk("spironolactone helps somatic items only (0-20% overall)",
        0 < red(res["spironolactone"]) < 20, "%.1f%%" % red(res["spironolactone"]))
    chk("SSRI + pill beats either treatment alone",
        red(res["sert_plus_drsp"]) > max(red(res["sert50_luteal"]),
                                         red(res["drsp_24_4"])),
        "%.1f%% vs %.1f / %.1f" % (red(res["sert_plus_drsp"]),
                                   red(res["sert50_luteal"]), red(res["drsp_24_4"])))
    chk("oophorectomy + E2 abolishes the CORE affective burden (>=80%)",
        core_reduction(ref, res["oophorectomy"]) >= 80,
        "%.1f%% core, %.1f%% DRSP-11"
        % (core_reduction(ref, res["oophorectomy"]), red(res["oophorectomy"])))
    chk("FALSIFY [P1]: the lag collapses below 2 d", fal["lag"] < 2.0,
        "%.1f d" % fal["lag"])
    chk("FALSIFY [P2]: dutasteride becomes dose-graded (0.5 mg >=35% of 2.5 mg)",
        red(res["fals_duta_0p5"], fal) >= 0.35 * max(1e-9, red(res["fals_duta_2p5"], fal)),
        "%.1f%% vs %.1f%%" % (red(res["fals_duta_0p5"], fal),
                              red(res["fals_duta_2p5"], fal)))
    chk("FALSIFY [P3]: both neurosteroid doses now hurt, monotonically",
        red(res["fals_allo_high"], fal) < red(res["fals_allo_low"], fal) < 0.0,
        "high %.0f%% < low %.0f%%" % (red(res["fals_allo_high"], fal),
                                      red(res["fals_allo_low"], fal)))
    chk("FALSIFY [P4]: the 24/4-vs-21/7 gap shrinks below 4 points",
        (red(res["fals_24_4"], fal) - red(res["fals_21_7"], fal)) < 4.0,
        "%.1f pt" % (red(res["fals_24_4"], fal) - red(res["fals_21_7"], fal)))
    chk("FALSIFY [P5]: the add-back response no longer fades",
        (fb_adapt - fb_supp) > 0.8 * (fb_flare - fb_supp),
        "%+.1f -> %+.1f points" % (fb_flare - fb_supp, fb_adapt - fb_supp))
    chk("FALSIFY: the SSRI benefit is much smaller than in the base model",
        red(res["fals_sert"], fal) < 0.6 * red(res["sert50_luteal"]),
        "%.1f%% vs %.1f%%" % (red(res["fals_sert"], fal), red(res["sert50_luteal"])))
    chk("leuprolide raises hot flushes by >=3/day",
        res["leuprolide"]["hf"] - ref["hf"] >= 3.0,
        "%.2f vs %.2f per day" % (res["leuprolide"]["hf"], ref["hf"]))
    chk("leuprolide loses 1.5-4% lumbar BMD over 6 cycles",
        -4.0 <= res["leuprolide"]["bmd"] <= -1.5, "%.2f%%" % res["leuprolide"]["bmd"])
    chk("add-back protects the skeleton (no further BMD loss per cycle)",
        res["addback_late"]["bmd_cyc"] > -0.10,
        "%+.3f %%/cycle vs %+.3f on leuprolide alone"
        % (res["addback_late"]["bmd_cyc"], res["leuprolide"]["bmd_cyc"]))

    npass = sum(1 for _, ok, _ in checks if ok)
    for name, ok, detail in checks:
        print("  [%s] %-60s %s" % ("PASS" if ok else "FAIL", name, detail))
    print("-" * 79)
    print("  %d/%d checks passed" % (npass, len(checks)))
    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
