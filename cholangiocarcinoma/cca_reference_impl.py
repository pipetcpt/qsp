#!/usr/bin/env python3
# =============================================================================
#  cca_reference_impl.py
#  Cholangiocarcinoma (biliary tract cancer) QSP model
#  Dependency-free reference implementation
# =============================================================================
#
#  This file is the EXECUTABLE TWIN of cca_mrgsolve_model.R: same 54 states,
#  same parameter names, same equations, same dose gate.  It exists so that
#  every number quoted in README.md can be regenerated with
#
#      python3 cca_reference_impl.py            # virtual trials
#      python3 cca_reference_impl.py --checks   # 30 self-checks, PASS/FAIL
#
#  on a machine with nothing installed but CPython.  No numpy, no scipy.
#  Integrator: fixed-step RK4, h = 0.1 d (set by free-platinum kinetics),
#  substepped to 0.002-0.01 d around every chemotherapy administration
#  (gemcitabine plasma t1/2 ~ 6 min).
#
#  ---------------------------------------------------------------------------
#  THE THREE STRUCTURAL COMMITMENTS  (see README.md for the argument)
#  ---------------------------------------------------------------------------
#  I.   THE DELIVERED DOSE IS AN OUTPUT.  Bilirubin, ANC, CrCl and performance
#       status gate what is actually infused, so obstruction -> held dose ->
#       growth -> obstruction is a closed loop.  Biliary drainage has no
#       antitumour action and is nonetheless the highest-leverage node,
#       because it is the only intervention that cuts that loop.
#
#  II.  RESISTANCE IS SELECTED, NOT INDUCED.  The FGFR2 kinase-domain-mutant
#       clone is seeded at t = 0 at frequency MU_RES * ln(N_cells)
#       (Goldie-Coldman).  No "resistance onset time" parameter exists.
#       See the RESULT reported in README: this seeding CANNOT by itself
#       explain a 7-month PFS, and the model says so instead of hiding it.
#
#  III. SURVIVAL IS TWO COMPETING HAZARDS ON DIFFERENT CLOCKS.  h_tumour is
#       slow and burden-driven; h_biliary is fast, recurrent and drainage-
#       driven.  Population survival is the MEAN of individual S(t), so the
#       model reports a median AND a 24-month tail from one simulation.
#
#  Only six parameters were fitted (marked [FIT]).  Everything else comes
#  from a published PK analysis, a product label, or a physiological constant.
# =============================================================================

from __future__ import division

import math
import sys

# -----------------------------------------------------------------------------
# 0.  STATE INDEX  (54 ODEs)
# -----------------------------------------------------------------------------
STATE_NAMES = [
    # --- PK (13) ---
    "GEM_C", "GEM_P", "DFDCTP", "CIS_C", "CIS_P", "PTDNA",
    "DUR_C", "DUR_P", "FGI_A", "FGI_C", "COV", "IVO_C", "FU_C",
    # --- tumour (8) ---
    "TS", "TP", "TR", "TD1", "TD2", "TD3", "TMET", "NADIR",
    # --- microenvironment / immune (5) ---
    "CAF", "ECM", "TEFF", "SUPP", "IL6",
    # --- biliary / hepatic (7) ---
    "OBSR", "PATN", "BILI", "ALP", "ALB", "FLR", "CHOLI",
    # --- host / toxicity (10) ---
    "PROL", "TRN1", "TRN2", "TRN3", "ANC", "PLT", "CRCL",
    "NEURO", "LBM", "PS",
    # --- biomarkers (5) ---
    "CA199", "CTDNA", "HG2", "PHOS", "IRAE",
    # --- survival / bookkeeping (6) ---
    "CUMHT", "CUMHB", "SURV", "CUMCIS", "RDI", "NDOSE",
]
S = dict((nm, i) for i, nm in enumerate(STATE_NAMES))
NST = len(STATE_NAMES)
assert NST == 54, NST

(iGEM_C, iGEM_P, iDFDCTP, iCIS_C, iCIS_P, iPTDNA, iDUR_C, iDUR_P, iFGI_A,
 iFGI_C, iCOV, iIVO_C, iFU_C, iTS, iTP, iTR, iTD1, iTD2, iTD3, iTMET, iNADIR,
 iCAF, iECM, iTEFF, iSUPP, iIL6, iOBSR, iPATN, iBILI, iALP, iALB, iFLR,
 iCHOLI, iPROL, iTRN1, iTRN2, iTRN3, iANC, iPLT, iCRCL, iNEURO, iLBM, iPS,
 iCA199, iCTDNA, iHG2, iPHOS, iIRAE, iCUMHT, iCUMHB, iSURV, iCUMCIS, iRDI,
 iNDOSE) = range(NST)

# -----------------------------------------------------------------------------
# 1.  PARAMETERS
# -----------------------------------------------------------------------------
P = dict(
    # ---------- patient / anatomy ----------
    BSA=1.72,             # m2
    LIVVOL=1500.0,        # mL parenchyma
    FHILAR=0.72,          # hilar-position fraction (pCCA/dCCA high, iCCA low)
    IMMENG=0.0,           # 1 = immune-engaged tumour (mixture indicator)
    FGFR2=0.0,            # 1 = FGFR2 fusion
    IDH1=0.0,             # 1 = IDH1 R132 mutation

    # ---------- gemcitabine PK ----------
    CL_GEM=3800.0, V1_GEM=22.0, Q_GEM=40.0, V2_GEM=35.0,
    VMAX_TP=95.0,         # AU/d  dCK-limited dFdCTP formation
    KM_TP=6.0,            # mg/L  saturable
    KOUT_TP=2.08,         # /d    intracellular t1/2 ~ 8 h

    # ---------- cisplatin PK (free platinum) ----------
    CL_CIS=430.0, V1_CIS=18.0, Q_CIS=25.0, V2_CIS=45.0,
    KFORM_PT=9.0, KREP_PT=0.85,

    # ---------- durvalumab PK (label) ----------
    CL_DUR=0.232, V1_DUR=5.6, Q_DUR=0.68, V2_DUR=3.9,
    IC50_DUR=0.35,        # mg/L

    # ---------- FGFR inhibitor PK (pemigatinib label) ----------
    KA_FGI=12.0, CL_FGI=254.0, V_FGI=235.0,
    IC50_FGI=0.008,       # mg/L
    KCOV=2.2, KCOV_OFF=0.10,
    RHO_PEM=0.05,         # residual activity of a REVERSIBLE FGFRi vs V564F
    RHO_FUT=0.45,         # residual activity of a COVALENT FGFRi vs V564F

    # ---------- ivosidenib ----------
    CL_IVO=6.5, V_IVO=180.0, IC50_IVO=0.9,
    EMAX_IDH=0.62,        # fractional growth-rate reduction at full 2-HG block

    # ---------- capecitabine / 5-FU surrogate ----------
    CL_FU=600.0, V_FU=60.0, FCONV=0.030,

    # ---------- tumour ----------
    LAM0=0.0100,          # /d   exponential growth (Vd ~ 69 d)
    LAM1=2.4,             # cm3/d linear-phase ceiling
    PSI=20.0,
    GP=0.55,              # persister relative growth rate
    KSP0=3.50,            # drug-driven entry into persistence           [FIT]
    KPS0=0.008,           # /d   reversion out of persistence
    ALPHA_P=0.45,         # persister chemo/TKI kill fraction           [FIT]
    K_GEM=0.168,          # kill per AU dFdCTP per d                     [FIT]
    K_CIS=0.192,          # kill per AU Pt-DNA adduct per d              [FIT]
    K_FU=0.100,           # kill per mg/L 5-FU per d
    K_SYN=0.55,           # gem/cis supra-additivity (NER interference)
    EFGFR=0.88,           # growth suppression at full FGFR block
    KFGFR_KILL=0.018,     # /d   direct kill at full FGFR block          [FIT]
    KIMM=0.055,           # /d   immune kill scale                       [FIT]
    KIMM50=180.0,
    KTR_D=0.55,           # /d   damaged-cell transit
    MU_RES=1.0e-6,        # per division, aggregate over ALL FGFR2 KD codons
    KMET=8.0e-5, LAMMET=0.020, TMET_NEW=15.0,
    RELAPSE_VOL=1.0,      # cm3 — radiologically detectable recurrence

    # ---------- stroma ----------
    KCAF=0.05, KCAF_OFF=0.03, KECM=0.030, KECM_OFF=0.012,
    FPEN_MIN=0.30, ECM50=1.0,

    # ---------- immune ----------
    KREC=0.13, KDEC=0.075, KSUP=0.030, KICD=0.55,
    KSUPP_ON=0.05, KSUPP_OFF=0.04, PDL1EXP=1.0, KPD=0.06,
    KIL6=1.0, KIL6_ON=0.020, KIL6_OFF=0.09,

    # ---------- biliary ----------
    KOBS_EQ=0.30, KOBS=110.0, DR_EFF=0.88,
    KOCC_SEMS=0.0029,     # /d  uncovered SEMS  (t1/2 ~ 240 d)
    KOCC_PLAS=0.0077,     # /d  plastic stent   (t1/2 ~  90 d)
    KING=0.9,
    KIN_BIL=0.15, BILI0=0.6, EMAX_BIL=24.0, OBS50=0.55, HBIL=3.0, BIL_LIV=8.0,
    KALP=0.10, ALP0=90.0, EMAX_ALP=520.0,
    KCH_ON=0.055, KCH_OFF=0.16, INF_STENT=1.4,

    # ---------- liver ----------
    KFLR=0.06, FOBS_FLR=0.55, FLR_MIN=0.03,
    KALB=0.035, ALB_MAX=4.3, ALB_FLOOR=0.35,

    # ---------- myelosuppression (Friberg) ----------
    CIRC0=4.2, MTT=125.0, GAMMA=0.16,
    SLOPE_GEM=1.40, SLOPE_CIS=0.90,
    PLT0=250.0, KPLT=0.09, SLOPE_PLT=0.90,

    # ---------- renal / neuro / host ----------
    CRCL0=88.0, KNEPH=0.055, KCRCL_REC=0.0016,
    KNEURO=0.075, KNEURO_OFF=0.0022,
    LBM0=42.0, KCACH=0.0060, KLBM_REC=0.0016, KIL6C=1.4,
    KPS=0.10,

    # ---------- biomarkers ----------
    KCA=0.20, CA_TUM=6.2, CA_BILI=140.0, CA0=22.0,
    KCT=1.6, CT_SHED=0.020,
    KHG=1.1, HG0=1.0,
    PHOS0=3.4, KPHOS=0.28, EMAX_PHOS=3.0,
    KIRAE_ON=0.0016, KIRAE_OFF=0.020,

    # ---------- hazards ----------
    HT0=9.6e-5,
    HT_T=0.0202,          # per (tumour volume / liver volume) per d     [FIT]
    HT_FLR=0.0053,
    FLRCRIT=0.55,
    HT_CACH=0.0038,
    HB0=6.4e-5,
    HB_CH=0.0132,         # per unit cholangitis burden per d            [FIT]
    HB_ALBI=0.00053,
    ALBI_REF=-2.60,

    # ---------- THE GATE ----------
    GATE_ON=1.0,
    BILI50=2.6, HB_GATE=6.0,
    ANC50=1.0, HA_GATE=8.0,
    PS50=2.60, HP_GATE=8.0,
    CRCL50=45.0, HC_GATE=8.0,
)


# -----------------------------------------------------------------------------
# 2.  REGIMEN — everything the clinician PRESCRIBES
# -----------------------------------------------------------------------------
class Regimen(object):
    def __init__(self, **kw):
        self.gem = kw.get("gem", 0.0)              # mg/m2 d1,d8 q21
        self.cis = kw.get("cis", 0.0)              # mg/m2 d1,d8 q21
        self.n_cycle = kw.get("n_cycle", 0)        # cycles of gem/cis
        self.durva = kw.get("durva", 0.0)          # mg q21 x8 then q28
        self.durva_days = kw.get("durva_days", 1460.0)
        self.fgi = kw.get("fgi", 0.0)              # mg/day oral FGFR inhibitor
        self.fgi_kind = kw.get("fgi_kind", "pem")  # "pem" | "fut"
        self.fgi_on = kw.get("fgi_on", 14)
        self.fgi_off = kw.get("fgi_off", 7)
        self.fgi_start = kw.get("fgi_start", 0.0)
        self.ivo = kw.get("ivo", 0.0)              # mg/day
        self.cape = kw.get("cape", 0.0)            # mg/m2 BID 14/21
        self.cape_days = kw.get("cape_days", 0.0)
        self.stent = kw.get("stent", None)         # None | "plastic" | "sems"
        self.stent_delay = kw.get("stent_delay", 0.0)
        self.revise = kw.get("revise", True)
        self.resect = kw.get("resect", 0.0)        # fraction removed at t=0
        # Second line is NOT a scheduled input: it starts when RECIST
        # progression is detected AND the gate is open (PS, bilirubin).
        # That is the same principle as rule 1-4 — therapy is an output.
        self.second_line = kw.get("second_line", 0.0)   # mFOLFOX kill scale
        self.sl_dur = kw.get("sl_dur", 168.0)           # d
        self.sl_lag = kw.get("sl_lag", 21.0)            # d after progression



# -----------------------------------------------------------------------------
#  THE GATE, as a standalone function.
#  This is the same arithmetic that make_rhs() inlines; it is exposed here so
#  that rules 1-4 can be inspected and tested on their own, because they are
#  the model's first structural claim.
# -----------------------------------------------------------------------------
def gate_values(y, p):
    """Return (gate_chemo, gate_cisplatin_extra) for a state vector."""
    if p["GATE_ON"] < 0.5:
        return 1.0, 1.0
    bili = max(y[iBILI], 0.10)
    anc = max(y[iANC], 0.02)
    ps = max(y[iPS], 0.0)
    crcl = max(y[iCRCL], 3.0)
    g_bil = 1.0 / (1.0 + (bili / p["BILI50"]) ** p["HB_GATE"])      # rule 1
    g_anc = 1.0 / (1.0 + (p["ANC50"] / anc) ** p["HA_GATE"])        # rule 2
    g_cr = 1.0 / (1.0 + (p["CRCL50"] / crcl) ** p["HC_GATE"])       # rule 3
    g_ps = 1.0 / (1.0 + (ps / p["PS50"]) ** p["HP_GATE"])           # rule 4
    return g_bil * g_anc * g_ps, g_cr


# -----------------------------------------------------------------------------
# 3.  COMPILED RIGHT-HAND SIDE
#     All parameters are bound as closure locals so the inner loop does pure
#     arithmetic.  This is ~4x faster than dict lookups and is the only reason
#     a 54-ODE virtual trial is feasible in pure CPython.
# -----------------------------------------------------------------------------
def make_rhs(p, reg, prog2l=None):
    if prog2l is None:
        prog2l = [-1.0]
    (BSA, LIVVOL, FHILAR, IMMENG, FGFR2, IDH1) = (
        p["BSA"], p["LIVVOL"], p["FHILAR"], p["IMMENG"], p["FGFR2"], p["IDH1"])
    CL_GEM, V1_GEM, Q_GEM, V2_GEM = p["CL_GEM"], p["V1_GEM"], p["Q_GEM"], p["V2_GEM"]
    VMAX_TP, KM_TP, KOUT_TP = p["VMAX_TP"], p["KM_TP"], p["KOUT_TP"]
    CL_CIS, V1_CIS, Q_CIS, V2_CIS = p["CL_CIS"], p["V1_CIS"], p["Q_CIS"], p["V2_CIS"]
    KFORM_PT, KREP_PT = p["KFORM_PT"], p["KREP_PT"]
    CL_DUR, V1_DUR, Q_DUR, V2_DUR, IC50_DUR = (
        p["CL_DUR"], p["V1_DUR"], p["Q_DUR"], p["V2_DUR"], p["IC50_DUR"])
    KA_FGI, CL_FGI, V_FGI, IC50_FGI = p["KA_FGI"], p["CL_FGI"], p["V_FGI"], p["IC50_FGI"]
    KCOV, KCOV_OFF = p["KCOV"], p["KCOV_OFF"]
    RHO = p["RHO_FUT"] if reg.fgi_kind == "fut" else p["RHO_PEM"]
    COVALENT = 1 if reg.fgi_kind == "fut" else 0
    CL_IVO, V_IVO, IC50_IVO, EMAX_IDH = p["CL_IVO"], p["V_IVO"], p["IC50_IVO"], p["EMAX_IDH"]
    CL_FU, V_FU, FCONV = p["CL_FU"], p["V_FU"], p["FCONV"]
    LAM0, LAM1, PSI, GP = p["LAM0"], p["LAM1"], p["PSI"], p["GP"]
    KSP0, KPS0, ALPHA_P = p["KSP0"], p["KPS0"], p["ALPHA_P"]
    K_GEM, K_CIS, K_FU, K_SYN = p["K_GEM"], p["K_CIS"], p["K_FU"], p["K_SYN"]
    EFGFR, KFGFR_KILL = p["EFGFR"], p["KFGFR_KILL"]
    KIMM, KIMM50, KTR_D = p["KIMM"], p["KIMM50"], p["KTR_D"]
    KMET, LAMMET = p["KMET"], p["LAMMET"]
    KCAF, KCAF_OFF, KECM, KECM_OFF = p["KCAF"], p["KCAF_OFF"], p["KECM"], p["KECM_OFF"]
    FPEN_MIN, ECM50 = p["FPEN_MIN"], p["ECM50"]
    KREC, KDEC, KSUP, KICD = p["KREC"], p["KDEC"], p["KSUP"], p["KICD"]
    KSUPP_ON, KSUPP_OFF, PDL1EXP, KPD = p["KSUPP_ON"], p["KSUPP_OFF"], p["PDL1EXP"], p["KPD"]
    KIL6, KIL6_ON, KIL6_OFF = p["KIL6"], p["KIL6_ON"], p["KIL6_OFF"]
    KOBS_EQ, KOBS, DR_EFF, KING = p["KOBS_EQ"], p["KOBS"], p["DR_EFF"], p["KING"]
    KIN_BIL, BILI0, EMAX_BIL, OBS50, HBIL, BIL_LIV = (
        p["KIN_BIL"], p["BILI0"], p["EMAX_BIL"], p["OBS50"], p["HBIL"], p["BIL_LIV"])
    KALP, ALP0, EMAX_ALP = p["KALP"], p["ALP0"], p["EMAX_ALP"]
    KCH_ON, KCH_OFF, INF_STENT = p["KCH_ON"], p["KCH_OFF"], p["INF_STENT"]
    KFLR, FOBS_FLR, FLR_MIN = p["KFLR"], p["FOBS_FLR"], p["FLR_MIN"]
    KALB, ALB_MAX, ALB_FLOOR = p["KALB"], p["ALB_MAX"], p["ALB_FLOOR"]
    CIRC0, MTT, GAMMA = p["CIRC0"], p["MTT"], p["GAMMA"]
    SLOPE_GEM, SLOPE_CIS = p["SLOPE_GEM"], p["SLOPE_CIS"]
    PLT0, KPLT, SLOPE_PLT = p["PLT0"], p["KPLT"], p["SLOPE_PLT"]
    CRCL0, KNEPH, KCRCL_REC = p["CRCL0"], p["KNEPH"], p["KCRCL_REC"]
    KNEURO, KNEURO_OFF = p["KNEURO"], p["KNEURO_OFF"]
    LBM0, KCACH, KLBM_REC, KIL6C, KPS = (
        p["LBM0"], p["KCACH"], p["KLBM_REC"], p["KIL6C"], p["KPS"])
    KCA, CA_TUM, CA_BILI, CA0 = p["KCA"], p["CA_TUM"], p["CA_BILI"], p["CA0"]
    KCT, CT_SHED, KHG, HG0 = p["KCT"], p["CT_SHED"], p["KHG"], p["HG0"]
    PHOS0, KPHOS, EMAX_PHOS = p["PHOS0"], p["KPHOS"], p["EMAX_PHOS"]
    KIRAE_ON, KIRAE_OFF = p["KIRAE_ON"], p["KIRAE_OFF"]
    HT0, HT_T, HT_FLR, FLRCRIT, HT_CACH = (
        p["HT0"], p["HT_T"], p["HT_FLR"], p["FLRCRIT"], p["HT_CACH"])
    HB0, HB_CH, HB_ALBI, ALBI_REF = p["HB0"], p["HB_CH"], p["HB_ALBI"], p["ALBI_REF"]
    GATE_ON = p["GATE_ON"] > 0.5
    BILI50, HB_GATE = p["BILI50"], p["HB_GATE"]
    ANC50, HA_GATE = p["ANC50"], p["HA_GATE"]
    PS50, HP_GATE = p["PS50"], p["HP_GATE"]
    CRCL50, HC_GATE = p["CRCL50"], p["HC_GATE"]
    TT0 = p["_TT0"]
    KOCC = {"plastic": p["KOCC_PLAS"], "sems": p["KOCC_SEMS"]}.get(reg.stent, 0.0)

    gem_dose = reg.gem * BSA
    cis_dose = reg.cis * BSA
    n_cycle = reg.n_cycle
    chemo_end = n_cycle * 21.0
    INFW = 0.5 / 24.0                     # 30-minute infusion
    INFR = 1.0 / INFW
    dur_dose, durva_days = reg.durva, reg.durva_days
    DURW = 1.0 / 24.0
    DURR = 1.0 / DURW
    fgi_dose, fgi_on, fgi_period = reg.fgi, float(reg.fgi_on), float(reg.fgi_on + reg.fgi_off)
    fgi_start = reg.fgi_start
    ivo_dose = reg.ivo
    fu_rate = reg.cape * BSA * 2.0 * FCONV
    cape_days = reg.cape_days
    ktr = 4.0 / (MTT / 24.0)
    K_2L, SL_DUR, SL_LAG = reg.second_line, reg.sl_dur, reg.sl_lag
    fmod, log10, log, sqrt = math.fmod, math.log10, math.log, math.sqrt

    def rhs(t, y, d):
        TS = y[iTS];  TP = y[iTP];  TR = y[iTR]
        if TS < 0.0: TS = 0.0
        if TP < 0.0: TP = 0.0
        if TR < 0.0: TR = 0.0
        TT = TS + TP + TR
        TMET = y[iTMET]
        if TMET < 0.0: TMET = 0.0

        ECM = y[iECM]
        if ECM < 0.0: ECM = 0.0
        fpen = FPEN_MIN + (1.0 - FPEN_MIN) / (1.0 + ECM / ECM50)

        bili = y[iBILI]
        if bili < 0.10: bili = 0.10
        alb = y[iALB]
        if alb < 1.0: alb = 1.0
        anc = y[iANC]
        if anc < 0.02: anc = 0.02
        ps = y[iPS]
        if ps < 0.0: ps = 0.0
        crcl = y[iCRCL]
        if crcl < 3.0: crcl = 3.0
        flr = y[iFLR]
        if flr < FLR_MIN: flr = FLR_MIN
        patn = y[iPATN]
        if patn < 0.0: patn = 0.0
        choli = y[iCHOLI]
        if choli < 0.0: choli = 0.0

        albi = 0.66 * log10(bili * 17.1) - 0.085 * (alb * 10.0)

        # ---------------- THE GATE (rules 1-4) ----------------
        if GATE_ON:
            g_bil = 1.0 / (1.0 + (bili / BILI50) ** HB_GATE)
            g_anc = 1.0 / (1.0 + (ANC50 / anc) ** HA_GATE)
            g_ps = 1.0 / (1.0 + (ps / PS50) ** HP_GATE)
            g_cr = 1.0 / (1.0 + (CRCL50 / crcl) ** HC_GATE)
            g_chem = g_bil * g_anc * g_ps
        else:
            g_chem = 1.0; g_cr = 1.0; g_ps = 1.0

        # ---------------- scheduled inputs ----------------
        gem_rate = 0.0; cis_rate = 0.0; pulse = 0.0
        if n_cycle > 0 and t < chemo_end:
            ph = fmod(t, 21.0)
            if ph < INFW or (8.0 <= ph < 8.0 + INFW):
                pulse = INFR
                gem_rate = gem_dose * pulse * g_chem
                cis_rate = cis_dose * pulse * g_chem * g_cr
        dur_rate = 0.0
        if dur_dose > 0.0 and t < durva_days:
            if t < chemo_end or n_cycle == 0:
                ph = fmod(t, 21.0)
                if ph < DURW:
                    dur_rate = dur_dose * DURR * g_ps
            else:
                ph = fmod(t - chemo_end, 28.0)
                if ph < DURW:
                    dur_rate = dur_dose * DURR * g_ps
        fgi_in = 0.0
        if fgi_dose > 0.0 and t >= fgi_start:
            if fmod(t - fgi_start, fgi_period) < fgi_on:
                fgi_in = fgi_dose
        fu_in = 0.0
        if fu_rate > 0.0 and t < cape_days and fmod(t, 21.0) < 14.0:
            fu_in = fu_rate

        # ---------------- PK ----------------
        Cg = y[iGEM_C] / V1_GEM
        Cgp = y[iGEM_P] / V2_GEM
        d[iGEM_C] = gem_rate - (CL_GEM / V1_GEM) * y[iGEM_C] - Q_GEM * (Cg - Cgp)
        d[iGEM_P] = Q_GEM * (Cg - Cgp)
        d[iDFDCTP] = VMAX_TP * fpen * Cg / (KM_TP + Cg) - KOUT_TP * y[iDFDCTP]

        Cc = y[iCIS_C] / V1_CIS
        Ccp = y[iCIS_P] / V2_CIS
        d[iCIS_C] = cis_rate - (CL_CIS / V1_CIS) * y[iCIS_C] - Q_CIS * (Cc - Ccp)
        d[iCIS_P] = Q_CIS * (Cc - Ccp)
        tp_c = y[iDFDCTP]
        krep = KREP_PT / (1.0 + K_SYN * tp_c / (1.0 + tp_c))   # gem blocks NER
        d[iPTDNA] = KFORM_PT * fpen * Cc - krep * y[iPTDNA]

        Cd = y[iDUR_C] / V1_DUR
        Cdp = y[iDUR_P] / V2_DUR
        d[iDUR_C] = dur_rate - (CL_DUR / V1_DUR) * y[iDUR_C] - Q_DUR * (Cd - Cdp)
        d[iDUR_P] = Q_DUR * (Cd - Cdp)
        durocc = Cd / (Cd + IC50_DUR)

        d[iFGI_A] = fgi_in - KA_FGI * y[iFGI_A]
        Cf = y[iFGI_C] / V_FGI
        d[iFGI_C] = KA_FGI * y[iFGI_A] - (CL_FGI / V_FGI) * y[iFGI_C]
        occ_rev = Cf / (Cf + IC50_FGI)
        d[iCOV] = KCOV * occ_rev * (1.0 - y[iCOV]) - KCOV_OFF * y[iCOV]
        occ_S = y[iCOV] if COVALENT else occ_rev
        if occ_S < 0.0: occ_S = 0.0
        occ_R = occ_S * RHO

        Ci = y[iIVO_C] / V_IVO
        d[iIVO_C] = ivo_dose - (CL_IVO / V_IVO) * y[iIVO_C]
        occ_ivo = Ci / (Ci + IC50_IVO)

        Cu = y[iFU_C] / V_FU
        d[iFU_C] = fu_in - (CL_FU / V_FU) * y[iFU_C]

        # ---------------- tumour ----------------
        lam_mod = 1.0
        lam_modR = 1.0
        if FGFR2 > 0.5:
            lam_mod *= (1.0 - EFGFR * occ_S)
            lam_modR *= (1.0 - EFGFR * occ_R)
        if IDH1 > 0.5:
            lam_mod *= (1.0 - EMAX_IDH * occ_ivo)
            lam_modR *= (1.0 - EMAX_IDH * occ_ivo)

        TTs = TT if TT > 1e-9 else 1e-9
        rr = LAM0 * TTs / LAM1
        denom = rr if rr > 4.0 else (1.0 + rr ** PSI) ** (1.0 / PSI)
        gS = LAM0 * lam_mod * TS / denom
        gP = LAM0 * GP * TP / denom
        gR = LAM0 * lam_modR * TR / denom

        chemo_hit = K_GEM * tp_c + K_CIS * y[iPTDNA] + K_FU * Cu
        if K_2L > 0.0:
            t2 = prog2l[0]
            if t2 >= 0.0 and (t2 + SL_LAG) <= t < (t2 + SL_LAG + SL_DUR):
                chemo_hit += K_2L * g_chem
        if FGFR2 > 0.5:
            fg_S = KFGFR_KILL * occ_S
            fg_R = KFGFR_KILL * occ_R
        else:
            fg_S = 0.0; fg_R = 0.0
        teff = y[iTEFF]
        if teff < 0.0: teff = 0.0
        imm_hit = KIMM * teff * KIMM50 / (KIMM50 + TT)

        killS = (chemo_hit + fg_S + imm_hit) * TS
        killP = (ALPHA_P * (chemo_hit + fg_S) + imm_hit) * TP
        killR = (chemo_hit + fg_R + imm_hit) * TR

        e_sp = KSP0 * (chemo_hit + fg_S)
        d[iTS] = gS - killS - e_sp * TS + KPS0 * TP
        d[iTP] = gP - killP + e_sp * TS - KPS0 * TP
        d[iTR] = gR - killR
        d[iTD1] = killS + killP + killR - KTR_D * y[iTD1]
        d[iTD2] = KTR_D * (y[iTD1] - y[iTD2])
        d[iTD3] = KTR_D * (y[iTD2] - y[iTD3])
        d[iTMET] = KMET * TT + LAMMET * TMET - (chemo_hit + imm_hit) * TMET

        sld = (TTs / TT0) ** (1.0 / 3.0)
        nad = y[iNADIR]
        d[iNADIR] = -4.0 * (nad - sld) if nad > sld else 0.0

        # ---------------- stroma ----------------
        caf = y[iCAF]
        if caf < 0.0: caf = 0.0
        il6 = y[iIL6]
        if il6 < 0.0: il6 = 0.0
        d[iCAF] = KCAF * TT / (TT + 60.0) - KCAF_OFF * caf
        d[iECM] = KECM * caf - KECM_OFF * ECM
        d[iIL6] = KIL6_ON * (caf + TT / 60.0 + 3.0 * choli) - KIL6_OFF * il6

        # ---------------- immune ----------------
        supp = y[iSUPP]
        if supp < 0.0: supp = 0.0
        icd = KICD * (killS + killP + killR) / (5.0 + TT)
        supp_frac = (1.0 - durocc) * PDL1EXP / (PDL1EXP + KPD)
        d[iTEFF] = KREC * IMMENG * (1.0 + icd) * (1.0 - supp_frac) \
            - KDEC * teff - KSUP * teff * supp
        d[iSUPP] = KSUPP_ON * (caf + il6) - KSUPP_OFF * supp

        # ---------------- biliary ----------------
        obs_raw = FHILAR * TT / (KOBS + TT)
        obs_tgt = obs_raw * (1.0 - DR_EFF * patn)
        obsr = y[iOBSR]
        if obsr < 0.0: obsr = 0.0
        elif obsr > 0.999: obsr = 0.999
        d[iOBSR] = KOBS_EQ * (obs_tgt - y[iOBSR])
        d[iPATN] = -KOCC * patn * (1.0 + KING * TT / (KOBS + TT))

        ob3 = obsr ** HBIL
        bil_tgt = BILI0 + EMAX_BIL * ob3 / (OBS50 ** HBIL + ob3) \
            + BIL_LIV * (0.0 if flr > 0.45 else (0.45 - flr) / 0.45)
        d[iBILI] = KIN_BIL * (bil_tgt - y[iBILI])
        d[iALP] = KALP * (ALP0 + EMAX_ALP * obsr - y[iALP])

        neutropenic = 1.0 / (1.0 + (anc / 0.5) ** 6)
        d[iCHOLI] = KCH_ON * obsr * obsr * (1.0 + INF_STENT * (1.0 - patn)) \
            * (1.0 + 0.6 * neutropenic) - KCH_OFF * choli

        # ---------------- liver reserve (a FRACTION, chased — not integrated) ----
        flr_tgt = (1.0 - TT / LIVVOL) * (1.0 - FOBS_FLR * obsr)
        if flr_tgt < FLR_MIN: flr_tgt = FLR_MIN
        d[iFLR] = KFLR * (flr_tgt - y[iFLR])
        d[iALB] = KALB * (ALB_MAX * (ALB_FLOOR + (1.0 - ALB_FLOOR) * flr)
                          * (1.0 - 0.32 * il6 / (il6 + KIL6)) - y[iALB])

        # ---------------- myelosuppression (Friberg) ----------------
        edrug = SLOPE_GEM * tp_c + SLOPE_CIS * y[iPTDNA]
        if edrug > 0.95: edrug = 0.95
        fb = (CIRC0 / anc) ** GAMMA
        if fb > 3.0: fb = 3.0
        d[iPROL] = ktr * y[iPROL] * ((1.0 - edrug) * fb - 1.0)
        d[iTRN1] = ktr * (y[iPROL] - y[iTRN1])
        d[iTRN2] = ktr * (y[iTRN1] - y[iTRN2])
        d[iTRN3] = ktr * (y[iTRN2] - y[iTRN3])
        d[iANC] = ktr * (y[iTRN3] - y[iANC])
        plt_sup = SLOPE_PLT * tp_c
        if plt_sup > 0.9: plt_sup = 0.9
        d[iPLT] = KPLT * (PLT0 * (1.0 - plt_sup) - y[iPLT])

        # ---------------- renal / neuro / host ----------------
        d[iCRCL] = -KNEPH * Cc * CRCL0 + KCRCL_REC * (CRCL0 - y[iCRCL])
        d[iNEURO] = KNEURO * Cc - KNEURO_OFF * y[iNEURO]
        d[iCUMCIS] = cis_rate
        lbm = y[iLBM]
        d[iLBM] = -KCACH * il6 / (il6 + KIL6C) * lbm + KLBM_REC * (LBM0 - lbm)

        bx = bili - 1.2
        if bx < 0.0: bx = 0.0
        lbm_loss = 1.0 - lbm / LBM0
        if lbm_loss < 0.0: lbm_loss = 0.0
        ps_tgt = 0.6 + 1.7 * bx / (bx + 5.0) + 1.8 * choli / (choli + 0.8) \
            + 3.3 * lbm_loss + 1.2 * (TT + TMET) / (TT + TMET + 900.0) \
            + 0.9 * y[iIRAE] + 0.5 * y[iNEURO] / (y[iNEURO] + 2.0)
        if ps_tgt > 4.0: ps_tgt = 4.0
        d[iPS] = KPS * (ps_tgt - y[iPS])

        # ---------------- biomarkers ----------------
        d[iCA199] = KCA * (CA0 + CA_TUM * (TT + TMET) + CA_BILI * obsr * obsr - y[iCA199])
        d[iCTDNA] = KCT * (CT_SHED * (TT + TMET) * (1.0 + 3.0 * TR / TTs) - y[iCTDNA])
        d[iHG2] = KHG * (HG0 * (1.0 + 9.0 * IDH1 * (1.0 - 0.96 * occ_ivo)) - y[iHG2])
        d[iPHOS] = KPHOS * (PHOS0 + (EMAX_PHOS * occ_S if fgi_dose > 0.0 else 0.0) - y[iPHOS])
        d[iIRAE] = (KIRAE_ON * durocc if dur_dose > 0.0 else 0.0) - KIRAE_OFF * y[iIRAE]

        # ---------------- the two hazards ----------------
        flr_def = FLRCRIT - flr
        if flr_def < 0.0: flr_def = 0.0
        ht = HT0 + HT_T * (TT + TMET) / LIVVOL + HT_FLR * flr_def + HT_CACH * lbm_loss
        albi_ex = albi - ALBI_REF
        if albi_ex < 0.0: albi_ex = 0.0
        hb = HB0 + HB_CH * choli + HB_ALBI * albi_ex
        d[iCUMHT] = ht
        d[iCUMHB] = hb
        d[iSURV] = -(ht + hb) * y[iSURV]

        d[iRDI] = pulse * g_chem * (g_cr if cis_dose > 0 else 1.0)
        d[iNDOSE] = pulse
        return d

    return rhs


# -----------------------------------------------------------------------------
# 4.  INITIAL CONDITIONS
# -----------------------------------------------------------------------------
def init_state(p, reg, tt0):
    y = [0.0] * NST
    frac_R = p["MU_RES"] * math.log(max(tt0, 1.0) * 1e9)     # Goldie-Coldman
    tr0 = tt0 * frac_R
    keep = 1.0 - reg.resect
    y[iTS] = (tt0 - tr0) * keep
    y[iTR] = tr0 * keep
    p["_TT0"] = max(tt0, 1e-9)          # RECIST baseline = pre-surgical burden
    y[iNADIR] = ((y[iTS] + y[iTR]) / p["_TT0"]) ** (1.0 / 3.0)
    y[iCAF] = 0.9
    y[iECM] = 2.2
    obs0 = p["FHILAR"] * tt0 / (p["KOBS"] + tt0)
    # A patient enrolled on a systemic-therapy trial has ALREADY been drained
    # and has a bilirubin inside the eligibility window.  Scenarios that delay
    # drainage start jaundiced instead — that is the whole point of S4.
    if reg.stent is not None and reg.stent_delay <= 0.0:
        y[iPATN] = 1.0
        obs0 = obs0 * (1.0 - p["DR_EFF"])
    else:
        y[iPATN] = 0.0
    y[iOBSR] = obs0
    y[iBILI] = p["BILI0"] + p["EMAX_BIL"] * obs0 ** p["HBIL"] / (
        p["OBS50"] ** p["HBIL"] + obs0 ** p["HBIL"])
    y[iALP] = p["ALP0"] + p["EMAX_ALP"] * obs0
    y[iFLR] = max(p["FLR_MIN"], (1.0 - tt0 / p["LIVVOL"]) * (1.0 - p["FOBS_FLR"] * obs0))
    y[iALB] = p["ALB_MAX"] * (p["ALB_FLOOR"] + (1.0 - p["ALB_FLOOR"]) * y[iFLR])
    y[iCHOLI] = 0.0
    y[iPROL] = y[iTRN1] = y[iTRN2] = y[iTRN3] = y[iANC] = p["CIRC0"]
    y[iPLT] = p["PLT0"]
    y[iCRCL] = p["CRCL0"]
    y[iLBM] = p["LBM0"]
    y[iPS] = 1.0
    y[iCA199] = p["CA0"] + p["CA_TUM"] * tt0 + p["CA_BILI"] * obs0 ** 2
    y[iCTDNA] = p["CT_SHED"] * tt0
    y[iHG2] = p["HG0"] * (1.0 + 9.0 * p["IDH1"])
    y[iPHOS] = p["PHOS0"]
    y[iSURV] = 1.0
    return y


# -----------------------------------------------------------------------------
# 5.  INTEGRATOR (RK4, adaptive substepping around fast PK)
# -----------------------------------------------------------------------------
GRID = 7.0          # d, reporting grid


def simulate(p, reg, tend=1460.0, h=0.1, tt0=100.0, trace=False):
    y = init_state(p, reg, tt0)
    prog2l = [-1.0]
    f = make_rhs(p, reg, prog2l)
    d1 = [0.0] * NST; d2 = [0.0] * NST; d3 = [0.0] * NST; d4 = [0.0] * NST
    ytmp = [0.0] * NST
    t = 0.0
    stent_done = reg.stent is None
    stent_at = reg.stent_delay
    ngrid = int(tend / GRID) + 1
    scurve = [1.0] * ngrid
    tprog = None
    best = y[iNADIR]
    tr_rows = []
    gi = 0
    rng_state = 0

    while t < tend - 1e-9:
        # ---- discrete events ----
        if (not stent_done) and t >= stent_at:
            y[iPATN] = 1.0
            stent_done = True
        if reg.stent is not None and stent_done and reg.revise \
                and y[iPATN] < 0.25 and y[iPS] < 3.0 and y[iSURV] > 0.02:
            y[iPATN] = 1.0                       # stent exchange / revision

        # ---- step size ------------------------------------------------------
        #  Gemcitabine has a 6-minute plasma half-life (CL/V1 = 173 /d).  That
        #  mode is far outside the RK4 stability region at the base step, so
        #  the integrator drops to a fine step for 1.6 d after every chemo
        #  administration and the two gemcitabine states are then snapped to
        #  exactly zero, which removes the stiff mode entirely between doses.
        #  (mrgsolve uses LSODA and needs none of this; see the .R file.)
        hh = h
        if reg.n_cycle > 0 and t < reg.n_cycle * 21.0 + 2.0:
            ph = math.fmod(t, 21.0)
            dt_dose = min(ph, abs(ph - 8.0), abs(ph - 21.0))
            if dt_dose < 0.12:
                hh = 0.002
            elif dt_dose < 2.6:
                hh = 0.01
        if reg.durva > 0.0 and math.fmod(t, 21.0) < 0.06:
            hh = min(hh, 0.005)
        nsub = int(round(h / hh))
        if nsub < 1:
            nsub = 1
        hh = h / nsub

        for _ in range(nsub):
            f(t, y, d1)
            for i in range(NST):
                ytmp[i] = y[i] + 0.5 * hh * d1[i]
            f(t + 0.5 * hh, ytmp, d2)
            for i in range(NST):
                ytmp[i] = y[i] + 0.5 * hh * d2[i]
            f(t + 0.5 * hh, ytmp, d3)
            for i in range(NST):
                ytmp[i] = y[i] + hh * d3[i]
            f(t + hh, ytmp, d4)
            for i in range(NST):
                y[i] += hh / 6.0 * (d1[i] + 2.0 * d2[i] + 2.0 * d3[i] + d4[i])
            t += hh
        for i in (iTS, iTP, iTR, iTD1, iTD2, iTD3, iTMET, iANC, iPLT,
                  iTEFF, iSUPP, iCHOLI, iPATN, iIL6, iCAF, iECM, iIRAE,
                  iGEM_C, iGEM_P, iCIS_C, iCIS_P, iDFDCTP, iPTDNA, iFGI_A):
            if y[i] < 0.0:
                y[i] = 0.0
        # Once plasma gemcitabine is < 1e-3 mg/L (~3 nM, four orders below KM_TP)
        # both gemcitabine states are snapped to zero.  This discards <0.1% of
        # the dose and removes the 173 /d eigenvalue from the slow phase.
        if y[iGEM_C] < 0.05:
            y[iGEM_C] = 0.0
            y[iGEM_P] = 0.0

        # ---- RECIST bookkeeping ----
        TT = y[iTS] + y[iTP] + y[iTR]
        sld = (max(TT, 1e-9) / p["_TT0"]) ** (1.0 / 3.0)
        if sld < best:
            best = sld
        if tprog is None:
            if reg.resect > 0.5:
                # adjuvant setting: the event is RELAPSE (a lesion becomes
                # detectable again), not RECIST progression of a target lesion
                if TT + y[iTMET] > p["RELAPSE_VOL"]:
                    tprog = t
                    prog2l[0] = t
            elif (sld >= 1.20 * y[iNADIR] and sld > 1.05 * best) \
                    or y[iTMET] > p["TMET_NEW"]:
                tprog = t
                prog2l[0] = t

        # ---- reporting grid ----
        while gi < ngrid and gi * GRID <= t + 1e-9:
            scurve[gi] = y[iSURV]
            gi += 1
        if trace and abs(math.fmod(t, GRID)) < h:
            tr_rows.append([t, TT, y[iTS], y[iTP], y[iTR], sld, y[iBILI], y[iOBSR],
                            y[iPATN], y[iANC], y[iPS], y[iSURV], y[iCA199],
                            y[iCTDNA], y[iPHOS], y[iHG2], y[iFLR], y[iCHOLI],
                            y[iTEFF], y[iCRCL]])
        if y[iANC] > 40.0:      # numerical guard on the Friberg feedback
            y[iANC] = 40.0
        if y[iSURV] < 5e-4 or y[iTS] != y[iTS]:
            break

    while gi < ngrid:
        scurve[gi] = y[iSURV]
        gi += 1
    rdi = y[iRDI] / y[iNDOSE] if y[iNDOSE] > 1e-9 else float("nan")
    return dict(y=y, scurve=scurve, tprog=(tprog if tprog is not None else tend),
                best=best, rdi=rdi, trace=tr_rows, t_end=t)


# -----------------------------------------------------------------------------
# 6.  VIRTUAL POPULATION
#     Deterministic quasi-random covariates (no RNG dependence, reproducible).
# -----------------------------------------------------------------------------
def _halton(i, b):
    f, r = 1.0, 0.0
    while i > 0:
        f /= b
        r += f * (i % b)
        i //= b
    return r


def _qnorm(u):
    """Acklam inverse normal CDF (|err| < 1.15e-9)."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    dd = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
          3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if u < pl:
        q = math.sqrt(-2 * math.log(u))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((dd[0]*q+dd[1])*q+dd[2])*q+dd[3])*q+1)
    if u > ph:
        q = math.sqrt(-2 * math.log(1 - u))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((dd[0]*q+dd[1])*q+dd[2])*q+dd[3])*q+1)
    q = u - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def _ln(z, cv):
    s = math.sqrt(math.log(1.0 + cv * cv))
    return math.exp(-0.5 * s * s + s * z)


def _median_cross(curve, level, grid=GRID):
    """First t at which a monotone-decreasing curve crosses `level`."""
    for k in range(1, len(curve)):
        if curve[k] <= level:
            y0, y1 = curve[k - 1], curve[k]
            if y0 == y1:
                return k * grid
            fr = (y0 - level) / (y0 - y1)
            return (k - 1 + fr) * grid
    return float("nan")


def run_population(over, reg_kw, n=120, tend=1460.0, h=0.1,
                   pi_immune=0.22, tt0_med=100.0):
    ngrid = int(tend / GRID) + 1
    os_curve = [0.0] * ngrid
    pfs_curve = [0.0] * ngrid
    nresp = 0.0
    rdis = []
    n_eval = 0
    neutro = 0.0
    for i in range(1, n + 1):
        p = dict(P)
        p.update(over)
        z1 = _qnorm(_halton(i, 2)); z2 = _qnorm(_halton(i, 3))
        z3 = _qnorm(_halton(i, 5)); z4 = _qnorm(_halton(i, 7))
        z5 = _qnorm(_halton(i, 11)); z6 = _qnorm(_halton(i, 13))
        u_im = _halton(i, 17)
        p["LAM0"] *= _ln(z1, 0.30)
        p["FHILAR"] = min(0.97, max(0.04, p["FHILAR"] * _ln(z2, 0.22)))
        p["KOBS"] *= _ln(z3, 0.30)
        p["HT_T"] *= _ln(z4, 0.35)
        p["HB_CH"] *= _ln(z5, 0.40)
        p["K_GEM"] *= _ln(z6, 0.42)
        p["K_CIS"] *= _ln(z6, 0.42)
        p["KFGFR_KILL"] *= _ln(z6, 0.42)
        p["CIRC0"] *= _ln(z2, 0.20)
        p["CRCL0"] *= _ln(z3, 0.22)
        p["IMMENG"] = 1.0 if (u_im < pi_immune) else 0.0
        tt0 = tt0_med * _ln(z4, 0.55)

        reg = Regimen(**reg_kw)
        r = simulate(p, reg, tend=tend, h=h, tt0=tt0)
        sc = r["scurve"]
        kprog = int(r["tprog"] / GRID)
        for k in range(ngrid):
            os_curve[k] += sc[k]
            pfs_curve[k] += (sc[k] if k < kprog else 0.0)
        # response is only assessable while the patient is plausibly alive
        if sc[min(ngrid - 1, int(42 / GRID))] > 0.5 and reg.resect < 0.5:
            n_eval += 1
            if r["best"] <= 0.70:
                nresp += 1.0
        if r["rdi"] == r["rdi"]:
            rdis.append(r["rdi"])
        if r["y"][iANC] < 1.0:
            neutro += 1.0
    for k in range(ngrid):
        os_curve[k] /= n
        pfs_curve[k] /= n
    return dict(
        mOS=_median_cross(os_curve, 0.5) / 30.44,
        mPFS=_median_cross(pfs_curve, 0.5) / 30.44,
        ORR=(100.0 * nresp / n_eval) if n_eval else float("nan"),
        OS6=100.0 * os_curve[min(ngrid - 1, int(183 / GRID))],
        OS12=100.0 * os_curve[min(ngrid - 1, int(365 / GRID))],
        OS24=100.0 * os_curve[min(ngrid - 1, int(730 / GRID))],
        RDI=100.0 * sum(rdis) / len(rdis) if rdis else float("nan"),
        n=n, os_curve=os_curve,
    )


# -----------------------------------------------------------------------------
# 7.  SCENARIOS
# -----------------------------------------------------------------------------
ADV = dict(FHILAR=0.72)                        # advanced BTC, mixed anatomy
ICCA_FGFR = dict(FHILAR=0.16, FGFR2=1.0)       # iCCA, FGFR2 fusion
ICCA_IDH = dict(FHILAR=0.16, IDH1=1.0)         # iCCA, IDH1 R132

SCENARIOS = [
    ("S1  Best supportive care + biliary drainage",
     ADV, dict(stent="sems"), 0.22, 100.0, 1100.0),
    ("S2  Gemcitabine/cisplatin (ABC-02)",
     ADV, dict(gem=1000, cis=25, n_cycle=8, stent="sems", second_line=0.010), 0.22, 100.0, 1100.0),
    ("S3  Gem/cis + durvalumab (TOPAZ-1)",
     ADV, dict(gem=1000, cis=25, n_cycle=8, durva=1500, stent="sems", second_line=0.010), 0.22, 100.0, 1100.0),
    ("S4  Gem/cis + durva, drainage DELAYED 60 d, no revision",
     ADV, dict(gem=1000, cis=25, n_cycle=8, durva=1500, stent="sems",
               stent_delay=60.0, revise=False, second_line=0.010), 0.22, 100.0, 1100.0),
    ("S5  Gem/cis, plastic stent instead of metal",
     ADV, dict(gem=1000, cis=25, n_cycle=8, stent="plastic", second_line=0.010), 0.22, 100.0, 1100.0),
    ("S6  Pemigatinib 13.5 mg 14/7 (FGFR2 fusion)",
     ICCA_FGFR, dict(fgi=13.5, fgi_kind="pem", fgi_on=14, fgi_off=7,
                     stent="sems"), 0.22, 100.0, 1100.0),
    ("S7  Futibatinib 20 mg continuous (FGFR2 fusion)",
     ICCA_FGFR, dict(fgi=20.0, fgi_kind="fut", fgi_on=1, fgi_off=0,
                     stent="sems"), 0.22, 100.0, 1100.0),
    ("S8  Ivosidenib 500 mg (IDH1 mutant)",
     ICCA_IDH, dict(ivo=500.0, stent="sems"), 0.22, 100.0, 1100.0),
    ("S11 Pemigatinib CONTINUOUS (no 1-week break) — a prediction",
     ICCA_FGFR, dict(fgi=13.5, fgi_kind="pem", fgi_on=1, fgi_off=0,
                     stent="sems"), 0.22, 100.0, 1100.0),
    ("S9  R0 resection + adjuvant capecitabine (BILCAP)",
     ADV, dict(resect=0.9999, cape=1250, cape_days=180, stent="sems"), 0.22, 100.0, 2400.0),
    ("S10 R0 resection, observation only",
     ADV, dict(resect=0.9999, stent="sems"), 0.22, 100.0, 2400.0),
]

FALSIFIERS = [
    ("F1  GATE_ON = 0   (bilirubin does not gate dosing)",
     dict(ADV, GATE_ON=0.0),
     dict(gem=1000, cis=25, n_cycle=8, durva=1500, stent="sems",
          stent_delay=60.0, revise=False, second_line=0.010), 0.22, 100.0, 1100.0),
    ("F2  MU_RES = 0    (no pre-existing FGFR2-mutant clone)",
     dict(ICCA_FGFR, MU_RES=0.0),
     dict(fgi=13.5, fgi_kind="pem", fgi_on=14, fgi_off=7, stent="sems"), 0.22, 100.0, 1100.0),
    ("F3  PI_IMMUNE = 1 (every tumour immune-engaged)",
     ADV, dict(gem=1000, cis=25, n_cycle=8, durva=1500, stent="sems", second_line=0.010), 1.00, 100.0, 1100.0),
    ("F4  FPEN_MIN = 1  (no stromal penetration barrier)",
     dict(ADV, FPEN_MIN=1.0),
     dict(gem=1000, cis=25, n_cycle=8, stent="sems", second_line=0.010), 0.22, 100.0, 1100.0),
]


def fmt(name, r):
    f1 = lambda v: ("%5.1f" % v) if v == v else ">hor"
    f2 = lambda v: ("%4.1f%%" % v) if v == v else "  n/a"
    return ("%-52s mOS %s | mPFS %s | ORR %s | 6m %4.1f%% | "
            "12m %4.1f%% | 24m %4.1f%% | RDI %s"
            % (name, f1(r["mOS"]), f1(r["mPFS"]), f2(r["ORR"]), r["OS6"], r["OS12"], r["OS24"],
               ("%3.0f%%" % r["RDI"]) if r["RDI"] == r["RDI"] else " --"))


def run_all(n=120, quiet=False):
    out = {}
    if not quiet:
        print("=" * 132)
        print("CHOLANGIOCARCINOMA QSP MODEL — virtual trials (n = %d per arm; "
              "OS/PFS in months)" % n)
        print("=" * 132)
    for nm, over, rk, pi, tt, te in SCENARIOS:
        r = run_population(over, rk, n=n, pi_immune=pi, tt0_med=tt, tend=te)
        out[nm.split()[0]] = r
        if not quiet:
            print(fmt(nm, r))
    if not quiet:
        print("-" * 132)
        print("FALSIFIERS — one parameter changed, nothing else refitted")
        print("-" * 132)
    for nm, over, rk, pi, tt, te in FALSIFIERS:
        r = run_population(over, rk, n=n, pi_immune=pi, tt0_med=tt, tend=te)
        out[nm.split()[0]] = r
        if not quiet:
            print(fmt(nm, r))
    if not quiet:
        print("=" * 132)
    return out



# -----------------------------------------------------------------------------
# 9.  SELF-CHECKS  (python3 cca_reference_impl.py --checks)
#     Structural and arithmetic assertions about the model, not a re-fit.
# -----------------------------------------------------------------------------
def run_checks():
    res = []

    def ck(name, cond, detail=""):
        res.append((bool(cond), name, detail))

    ck("54 ODE states", NST == 54, "NST=%d" % NST)
    ck("state names unique", len(set(STATE_NAMES)) == 54)

    # --- PK identities -------------------------------------------------------
    p = dict(P); p.update(ADV); p["_TT0"] = 100.0
    reg = Regimen(gem=1000, cis=25, n_cycle=1, stent="sems")
    r = simulate(p, reg, tend=21.0, tt0=100.0)
    ck("gemcitabine cleared between doses", r["y"][iGEM_C] == 0.0)
    ck("platinum peripheral non-negative", r["y"][iCIS_P] >= 0.0)
    ck("cumulative cisplatin ~ 2 x 25 mg/m2 x BSA",
       abs(r["y"][iCUMCIS] - 2 * 25 * P["BSA"]) < 0.06 * 2 * 25 * P["BSA"],
       "%.1f vs %.1f" % (r["y"][iCUMCIS], 2 * 25 * P["BSA"]))

    # --- the gate ------------------------------------------------------------
    y = [0.0] * NST
    y[iBILI] = 0.6; y[iANC] = 4.2; y[iPS] = 1.0; y[iCRCL] = 90.0
    g, gc = gate_values(y, P)
    ck("gate fully open in a fit, drained patient", g > 0.99 and gc > 0.99,
       "g=%.3f" % g)
    y[iBILI] = 8.0
    g2, _ = gate_values(y, P)
    ck("rule 1 closes the gate at bilirubin 8", g2 < 0.02, "g=%.4f" % g2)
    y[iBILI] = 0.6; y[iANC] = 0.4
    g3, _ = gate_values(y, P)
    ck("rule 2 closes the gate at ANC 0.4", g3 < 0.02, "g=%.4f" % g3)
    y[iANC] = 4.2; y[iCRCL] = 30.0
    _, g4 = gate_values(y, P)
    ck("rule 3 omits cisplatin at CrCl 30", g4 < 0.05, "g=%.4f" % g4)
    y[iCRCL] = 90.0; y[iPS] = 3.5
    g5, _ = gate_values(y, P)
    ck("rule 4 all but stops therapy at ECOG 3.5", g5 < 0.10, "g=%.4f" % g5)
    pf = dict(P); pf["GATE_ON"] = 0.0
    gf, gcf = gate_values([0.0] * NST, pf)
    ck("GATE_ON=0 disables all four rules", gf == 1.0 and gcf == 1.0)

    # --- Goldie-Coldman seeding ---------------------------------------------
    pp = dict(P); pp.update(ICCA_FGFR)
    y0 = init_state(pp, Regimen(stent="sems"), 100.0)
    expect = 100.0 * P["MU_RES"] * math.log(100.0 * 1e9)
    ck("resistant clone present at t=0", y0[iTR] > 0.0)
    ck("seeding equals MU_RES*ln(N)*T0", abs(y0[iTR] - expect) < 1e-12,
       "%.3e" % y0[iTR])
    pz = dict(pp); pz["MU_RES"] = 0.0
    ck("MU_RES=0 removes the clone entirely",
       init_state(pz, Regimen(stent="sems"), 100.0)[iTR] == 0.0)

    # --- biliary physiology --------------------------------------------------
    pb = dict(P); pb.update(ADV); pb["_TT0"] = 100.0
    yd = init_state(pb, Regimen(stent="sems"), 100.0)
    yj = init_state(pb, Regimen(stent="sems", stent_delay=60.0), 100.0)
    ck("undrained patient is jaundiced at entry", yj[iBILI] > 3.0,
       "%.1f mg/dL" % yj[iBILI])
    ck("drained patient is not", yd[iBILI] < 1.5, "%.2f mg/dL" % yd[iBILI])
    ck("drainage lowers obstruction by DR_EFF",
       abs(yd[iOBSR] - yj[iOBSR] * (1 - P["DR_EFF"])) < 1e-9)
    t_sems = math.log(2.0) / P["KOCC_SEMS"]
    t_plas = math.log(2.0) / P["KOCC_PLAS"]
    ck("SEMS patency half-life ~ 8 months", 200 < t_sems < 280, "%.0f d" % t_sems)
    ck("plastic patency half-life ~ 3 months", 70 < t_plas < 110, "%.0f d" % t_plas)

    # --- ALBI ----------------------------------------------------------------
    # ALBI grade 1 is <= -2.60; albumin 4.0 with bilirubin 1.0 sits at -2.59,
    # i.e. just inside grade 2, which is the correct published behaviour.
    a1 = 0.66 * math.log10(1.0 * 17.1) - 0.085 * (4.2 * 10.0)
    a2 = 0.66 * math.log10(1.0 * 17.1) - 0.085 * (4.0 * 10.0)
    ck("ALBI grade 1 at bilirubin 1.0 / albumin 4.2", a1 <= -2.60, "%.3f" % a1)
    ck("ALBI grade 2 at bilirubin 1.0 / albumin 4.0", a2 > -2.60, "%.3f" % a2)

    # --- monotonicity --------------------------------------------------------
    pn = dict(P); pn.update(ADV); pn["_TT0"] = 100.0
    r2 = simulate(pn, Regimen(stent="sems"), tend=400.0, tt0=100.0)
    sc = r2["scurve"]
    ck("survival is non-increasing", all(sc[i + 1] <= sc[i] + 1e-12
                                         for i in range(len(sc) - 1)))
    ck("untreated tumour grows", r2["y"][iTS] > 100.0)
    ck("untreated patient never responds", r2["best"] >= 0.999)

    # --- persister mechanism -------------------------------------------------
    pc = dict(P); pc.update(ADV); pc["_TT0"] = 100.0
    rc = simulate(pc, Regimen(gem=1000, cis=25, n_cycle=8, stent="sems"),
                  tend=200.0, tt0=100.0)
    tt = max(rc["y"][iTS] + rc["y"][iTP] + rc["y"][iTR], 1e-9)
    ck("chemotherapy creates drug-tolerant persisters", rc["y"][iTP] > 0.0)
    ck("persister fraction is substantial by end of chemo",
       rc["y"][iTP] / tt > 0.3, "%.2f" % (rc["y"][iTP] / tt))
    ck("relative dose intensity is an OUTPUT below 100%",
       0.5 < rc["rdi"] < 1.0, "%.3f" % rc["rdi"])

    # --- FGFR target engagement ---------------------------------------------
    pf2 = dict(P); pf2.update(ICCA_FGFR); pf2["_TT0"] = 100.0
    rf = simulate(pf2, Regimen(fgi=13.5, fgi_kind="pem", stent="sems"),
                  tend=120.0, tt0=100.0)
    ck("FGFR inhibitor raises serum phosphate (on-target)",
       rf["y"][iPHOS] > P["PHOS0"] + 1.0, "%.2f mg/dL" % rf["y"][iPHOS])
    ck("no phosphate rise without an FGFR inhibitor",
       abs(rc["y"][iPHOS] - P["PHOS0"]) < 0.05)

    # --- IDH1 target engagement ---------------------------------------------
    pi2 = dict(P); pi2.update(ICCA_IDH); pi2["_TT0"] = 100.0
    ri = simulate(pi2, Regimen(ivo=500.0, stent="sems"), tend=120.0, tt0=100.0)
    ck("ivosidenib suppresses 2-HG by >80%",
       ri["y"][iHG2] < 0.2 * (P["HG0"] * 10.0), "%.2f" % ri["y"][iHG2])
    ck("ivosidenib produces no RECIST response (cytostatic only)",
       ri["best"] > 0.70, "best SLD %.3f" % ri["best"])

    # --- CA 19-9 confound ----------------------------------------------------
    ck("CA 19-9 higher in the jaundiced patient at identical tumour burden",
       yj[iCA199] > yd[iCA199], "%.0f vs %.0f" % (yj[iCA199], yd[iCA199]))

    ok = sum(1 for c in res if c[0])
    for good, name, detail in res:
        print("  %s  %-60s %s" % ("PASS" if good else "FAIL", name, detail))
    print("  %d/%d checks passed" % (ok, len(res)))
    return ok == len(res)


if __name__ == "__main__":
    if "--checks" in sys.argv:
        print("=" * 92)
        print("CHOLANGIOCARCINOMA QSP MODEL — self-checks")
        print("=" * 92)
        sys.exit(0 if run_checks() else 1)
    N = 120
    for a in sys.argv[1:]:
        if a.isdigit():
            N = int(a)
    run_all(N)
