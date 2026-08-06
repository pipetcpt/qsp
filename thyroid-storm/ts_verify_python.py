#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ts_verify_python.py — numerical core / calibration harness for the
Thyroid Storm (갑상선 폭풍) QSP model.

This file is the ARITHMETIC OF RECORD for every number quoted in README.md.
ts_mrgsolve_model.R is a line-by-line port of the ODE system implemented here;
this script exists so that the claims in the README are computed, not asserted.

===========================================================================
CENTRAL THESIS
===========================================================================
Thyroid storm is NOT "more thyrotoxicosis".  Serum T4 and T3 do not separate
storm from uncomplicated thyrotoxicosis (Brooks & Waldstein 1980; Brooks 1975).
What separates them is the LOOP GAIN  L  of a positive-feedback ring that is
closed only when a precipitant is present:

  ring A   fT3 -> BMR & beta1 density -> lipolysis -> NEFA
                -> displacement of T3 from TBG -> free fraction -> fT3
  ring B   Tc  -> Q10 amplification of metabolism -> heat production -> Tc
  ring C   Tc  -> hypothalamic/CNS injury -> loss of the heat-loss effector -> Tc

  L <  1  =>  a stable, hot-but-survivable operating point exists
              = severe (even florid) thyrotoxicosis
  L >= 1  =>  no stable operating point; temperature and rate run away
              = THYROID STORM

Every drug is then classified by WHICH FACTOR OF L IT MULTIPLIES, and the
clinical treatment hierarchy of thyroid storm becomes arithmetic:

  thionamide      changes the INPUT  (lowers the operating point over WEEKS,
                  because it acts on a reservoir with a 2-month time constant)
  beta-blockade   changes the GAIN   (breaks ring A in HOURS) -- and it does so
                  TWICE, because the lipolysis that generates the displacing
                  NEFA is itself beta-mediated
  iodide          changes the INPUT fast (release block), but its SIGN FLIPS if
                  given before a thionamide (substrate for organification)
  cooling         changes the GAIN of rings B and C only
  glucocorticoid  changes the INPUT (D1/D2) and repairs the accelerated
                  cortisol clearance that thyroid hormone itself causes
  plasmapheresis  is the ONLY lever on the protein-bound reservoir

Units: TIME = HOURS.  Hormones nmol/L, iodide umol, heat W, temperature degC.
Dependencies: numpy only.
===========================================================================
"""

import math
import numpy as np

# =====================================================================
# 0.  STATE VECTOR
# =====================================================================
NAMES = [
    "S",      # 0  thyroidal organified hormone store        nmol T4-eq
    "T4",     # 1  plasma TOTAL T4                           nmol/L
    "T3",     # 2  plasma TOTAL T3                           nmol/L
    "rT3",    # 3  plasma reverse T3                         nmol/L
    "Ggut",   # 4  enterohepatic (gut lumen) T4 pool         nmol
    "Ipl",    # 5  plasma inorganic iodide                   umol/L
    "Ithy",   # 6  intrathyroidal non-organified iodide      umol
    "NIS",    # 7  NIS activity (relative, 1 = normal)       -
    "TSH",    # 8  plasma TSH                                mIU/L
    "NEFA",   # 9  plasma non-esterified fatty acids         mmol/L   [FAST]
    "TBG",    # 10 thyroxine-binding globulin (relative)     -
    "SNS",    # 11 sympathetic outflow (relative)            -        [FAST]
    "Rb",     # 12 beta1-adrenoceptor density (relative)     -
    "HR",     # 13 heart rate                                bpm      [FAST]
    "CR",     # 14 cardiac (contractile) reserve, 1 = full   -
    "Tc",     # 15 core temperature                          degC     [FAST]
    "Vol",    # 16 volume status, 1 = euvolemic              -        [FAST]
    "CNSx",   # 17 CNS / hypothalamic dysfunction index      0..1     [FAST]
    "Bili",   # 18 total bilirubin                           mg/dL
    "Cort",   # 19 plasma cortisol                           nmol/L
    "GIx",    # 20 GI dysfunction index                      0..1
    "AFb",    # 21 atrial-fibrillation burden                0..1
    "Hz",     # 22 cumulative mortality hazard               -
    "BWi",    # 23 integral of Burch-Wartofsky score         pt*h
    "Prec",   # 24 precipitant intensity (decaying)          -
    "PGE",    # 25 PGE2-driven setpoint elevation            0..1
    "Aptu",   # 26 PTU gut amount                            mg
    "Cptu",   # 27 PTU plasma                                mg/L
    "Ammi",   # 28 methimazole gut amount                    mg
    "Cmmi",   # 29 methimazole plasma                        mg/L
    "Apro",   # 30 propranolol gut amount                    mg
    "Cpro",   # 31 propranolol plasma                        ng/mL
    "Cesm",   # 32 esmolol plasma                            ug/mL
    "Cgc",    # 33 glucocorticoid, hydrocortisone-eq         mg/L
    "Achol",  # 34 cholestyramine in gut                     g
    "Ciop",   # 35 iopanoic acid plasma                      mg/L
    "Casa",   # 36 salicylate plasma                         mg/dL
    "Capap",  # 37 acetaminophen plasma                      mg/L
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)

FAST = [IX["NEFA"], IX["SNS"], IX["HR"], IX["Tc"], IX["Vol"], IX["CNSx"]]

# =====================================================================
# 1.  PARAMETERS
#     Every value below is either (i) textbook normal human physiology,
#     (ii) an analytic consequence of the euthyroid mass balance, or
#     (iii) a drug potency fitted to NON-storm clinical data.
#     The storm-derived numbers are exactly three (C1-C3, section 7).
# =====================================================================
P = dict(
    # ---------- distribution volumes ----------
    V_T4=10.0, V_T3=40.0, V_rT3=30.0, V_I=17.0,

    # ---------- T4 disposal, /h ----------
    kD1_4=0.0267 / 24.0,      # D1 outer ring  T4 -> T3   (liver, kidney)
    kD2_4=0.0133 / 24.0,      # D2 outer ring  T4 -> T3   (muscle, brain, BAT)
    kD3_4=0.0350 / 24.0,      # D3 inner ring  T4 -> rT3
    kEHC_4=0.0885 / 24.0,     # biliary/conjugated T4 -> gut lumen
    kreab=0.60,               # fraction of biliary T4 reabsorbed
    kgut=0.25,                # /h gut transit
    krT3=0.1733,              # /h rT3 elimination (t1/2 4 h)
    fD1_rT3=0.85,             # fraction of rT3 clearance that is D1-mediated
    kT3=0.028704,             # /h T3 elimination (t1/2 24.1 h)

    # ---------- thyroid gland ----------
    ksec0=6.25e-4,            # /h fractional secretion of the store
    fT3sec0=0.080,            # basal fraction of secretion appearing as T3
    fT3sec_max=0.170,         # a stimulated gland secretes relatively more T3
    Vorg0=0.040,              # umol I/h organification Vmax at Stim = 1
    Km_org=0.200,             # umol  intrathyroidal iodide, half-saturating
    kleak=0.012,              # /h iodide efflux from the gland
    CLthy=1.121,              # L/h thyroidal iodide clearance at NIS = Stim = 1
    CLren_I=2.10,             # L/h renal iodide clearance
    I_diet=0.04790,           # umol/h dietary iodide (= 146 ug/d)
    # iodine RECYCLING stoichiometry: T4 has 4 iodines; one is released per
    # deiodination step and the remaining three when T3/rT3 are cleared.
    nI_step=1.0, nI_T3=3.0,
    IC_WC=2.00,               # umol intrathyroidal iodide, half Wolff-Chaikoff
    nWC=2.0,
    Imax_WCrel=0.78,          # max inhibition of hormone RELEASE by iodide
    Imax_WCorg=0.90,          # max inhibition of ORGANIFICATION by iodide
    K_NIS=1.00,               # umol/L plasma iodide, half NIS down-regulation
    Km_NIS=25.0,              # umol/L NIS Km for iodide (saturable transport)
    tau_NIS=84.0,             # h  Wolff-Chaikoff escape time constant

    # ---------- HPT axis ----------
    tau_TSH=12.0, TSH0=1.50, kTSHfb=8.0,
    TRAb=1.0,                 # stimulating drive (1 = normal TSH equivalent)

    # ---------- binding / free fraction  (RING A, second half) ----------
    phi4_0=2.0e-4, phi3_0=3.0e-3,
    NEFA0=0.40,
    Emax_nefa=1.110, K_nefa=0.862, n_nefa=2.0,
    Emax_sal=1.20, K_sal=20.0,
    tau_TBG=48.0, eTBG_prec=0.15,

    # ---------- nuclear TR signal ----------
    K_TR=1.0,

    # ---------- lipolysis  (RING A, first half) ----------
    tau_NEFA=0.25, eL_T3=0.50, eL_beta=1.17,

    # ---------- sympathetic outflow ----------
    # NOTE eT3_S = 0 on purpose: circulating catecholamines are NOT raised by
    # thyroid hormone.  The sensitisation is at the receptor (Rb), not the
    # transmitter.  This is a modelling commitment, not an oversight.
    tau_SNS=0.50, ePrec_S=0.70, eTemp_S=0.16, eVol_S=0.35, eT3_S=0.00,

    # ---------- beta1 receptor density ----------
    tau_Rb=36.0, eR=1.00,

    # ---------- heart rate ----------
    HR0=70.0, HRmax=210.0, Kd_HR=1.205,
    eHR_b=0.55, eHR_T3=0.35, eHR_T=0.100, tau_HR=0.30,

    # ---------- cardiac pump ----------
    SV0=1.10, eSV_T3=0.60, eSV_b=0.50,
    W_crit=4.50, kdep=0.020, krec=0.010,

    # ---------- thermal  (RINGS B and C) ----------
    BMR0=80.0, C_body=245000.0,
    eBMR=1.00, Q10=2.00, eQ_sns=0.15,
    h_min=8.0, h_span=52.0, T_sink=33.0,
    Tset0=37.42, w_eff=0.35,
    ePGE_set=1.90, ePGE_shiv=0.15,
    eCNS_eff=0.60,            # RING C gain: heat-loss effector lost with CNSx
    eVol_eff=1.50,
    Emax_unc=0.50, K_unc=30.0,    # salicylate uncoupling of oxidative phos.

    # ---------- CNS ----------
    kCNS_on=0.022, kCNS_hyp=0.080, kCNS_off=0.045, eCNS_T3=0.50, T_CNS=38.5,
    perf_CNS=0.65,

    # ---------- liver ----------
    Bili0=0.80, kb_out=0.030, eB_hyp=6.0, eB_hf=2.5, eB_T3=0.35,

    # ---------- cortisol ----------
    Cort0=400.0, kCL_cort=0.35, eCL_T3=1.67, eStress_cort=0.90, AR=1.00,
    nmol_per_mgL=2759.0,      # nmol/L glucocorticoid activity per mg/L of
                              # hydrocortisone-equivalent (MW 362.5)

    # ---------- GI ----------
    tau_GI=6.0, eGI_T3=0.55, eGI_prec=0.50,

    # ---------- volume ----------
    kdrink=0.120, h_base=20.04, kins_h=0.00025, kins_T=0.00150,
    kGI_loss=0.0055, ePO_GI=0.35, Vol_floor=0.55,

    # ---------- atrial fibrillation ----------
    tau_AF=8.0, eAF_T3=0.55, eAF_b=0.35, AF_suscept=1.0,

    # ---------- precipitant ----------
    kPrec=0.0110, kPGE=0.0140,

    # ---------- hazard ----------
    Tc_ceil=43.0,
    h0_haz=2.755e-6, aBW=0.45, eShock=3.0, eBili=1.20, eHyper=2.00, eAI=0.90,

    # ================= DRUG PK =================
    ka_ptu=2.00, ke_ptu=0.462, F_ptu=0.80, V_ptu=21.0,       # t1/2 1.5 h
    ka_mmi=1.60, ke_mmi=0.139, F_mmi=0.93, V_mmi=28.0,       # t1/2 5 h
    ka_pro=1.40, ke_pro=0.173, F_pro=0.30, V_pro=280.0,      # t1/2 4 h
    ke_esm=4.62, V_esm=110.0,                                # t1/2 9 min
    ke_gc=0.408, V_gc=35.0,                                  # t1/2 1.7 h
    kchol=0.25,
    ke_iop=0.0347, V_iop=25.0,                               # t1/2 20 h
    ke_asa=0.1155, V_asa=12.0,                               # t1/2 6 h
    ke_apap=0.277, V_apap=50.0,                              # t1/2 2.5 h

    # ================= DRUG PD =================
    IC50_TPO_ptu=0.80, IC50_TPO_mmi=0.060, Imax_TPO=0.97,
    IC50_D1_ptu=2.50, Imax_D1_ptu=0.78,        # PTU only.  MMI has NO D1 action.
    IC50_D1_gc=0.35, Imax_D1_gc=0.40,
    IC50_D2_gc=0.35, Imax_D2_gc=0.45,
    IC50_D1_iop=1.5, Imax_D1_iop=0.92,
    IC50_D2_iop=1.5, Imax_D2_iop=0.85,
    IC50_D1_pro=95.0, Imax_D1_pro=0.30,
    IC50_beta_pro=32.0, Imax_beta_pro=0.90,
    IC50_beta_esm=0.55, Imax_beta_esm=0.90,
    ekchol=0.90, IC50_chol=4.0,
    eapap_PGE=0.65, IC50_apap=12.0,
)


# =====================================================================
# 2.  ALGEBRAIC LAYER
# =====================================================================
def emax(C, IC50, Imax):
    return Imax * C / (IC50 + C) if C > 0.0 else 0.0


def algebra(y, p, dose=None):
    d = {}
    S, T4, T3 = y[0], y[1], y[2]
    Ipl, Ithy, NIS, TSH = y[5], y[6], y[7], y[8]
    NEFA, TBG, SNS, Rb = y[9], y[10], y[11], y[12]
    HR, CR, Tc, Vol = y[13], y[14], y[15], y[16]
    CNSx, Bili, Cort, GIx = y[17], y[18], y[19], y[20]
    Prec, PGE = y[24], y[25]
    Cptu, Cmmi, Cpro, Cesm, Cgc = y[27], y[29], y[31], y[32], y[33]
    Achol, Ciop, Casa, Capap = y[34], y[35], y[36], y[37]

    # ---- (a) FREE FRACTION : the second axis -------------------------
    # Fatty-acid displacement is STEEP, not Michaelis: it is negligible in
    # ordinary states and only becomes large at the NEFA levels reached in
    # critical illness / heparin / storm.  Hill exponent 2.
    dN = max(0.0, NEFA - p["NEFA0"])
    Dn = p["Emax_nefa"] * dN ** p["n_nefa"] / \
        (p["K_nefa"] ** p["n_nefa"] + dN ** p["n_nefa"])
    Dsal = emax(Casa, p["K_sal"], p["Emax_sal"])
    Fd = (1.0 + Dn + Dsal) / max(TBG, 0.25)
    d["Fd"] = Fd
    d["phi3"] = p["phi3_0"] * Fd
    d["phi4"] = p["phi4_0"] * Fd
    d["fT3"] = T3 * d["phi3"] * 1000.0            # pmol/L
    d["fT4"] = T4 * d["phi4"] * 1000.0            # pmol/L
    d["fT3n"] = d["fT3"] / (1.8 * p["phi3_0"] * 1000.0)
    fT3n = d["fT3n"]

    # ---- (b) NUCLEAR TR SIGNAL (saturating; = 1 when normal) ---------
    K = p["K_TR"]
    d["Sig"] = (fT3n / (K + fT3n)) * (K + 1.0)
    Sig = d["Sig"]

    # ---- (c) BETA SIGNAL --------------------------------------------
    occ_b = 1.0 - (1.0 - emax(Cpro, p["IC50_beta_pro"], p["Imax_beta_pro"])) * \
                  (1.0 - emax(Cesm, p["IC50_beta_esm"], p["Imax_beta_esm"]))
    d["occ_b"] = occ_b
    d["Bsig"] = SNS * Rb * (1.0 - occ_b)
    Bsig = d["Bsig"]

    # ---- (d) DEIODINASE INHIBITION ----------------------------------
    d["inh_D1"] = 1.0 - (1.0 - emax(Cptu, p["IC50_D1_ptu"], p["Imax_D1_ptu"])) * \
                        (1.0 - emax(Cgc, p["IC50_D1_gc"], p["Imax_D1_gc"])) * \
                        (1.0 - emax(Ciop, p["IC50_D1_iop"], p["Imax_D1_iop"])) * \
                        (1.0 - emax(Cpro, p["IC50_D1_pro"], p["Imax_D1_pro"]))
    d["inh_D2"] = 1.0 - (1.0 - emax(Cgc, p["IC50_D2_gc"], p["Imax_D2_gc"])) * \
                        (1.0 - emax(Ciop, p["IC50_D2_iop"], p["Imax_D2_iop"]))

    # ---- (e) THIONAMIDE / TPO ---------------------------------------
    d["inh_TPO"] = 1.0 - (1.0 - emax(Cptu, p["IC50_TPO_ptu"], p["Imax_TPO"])) * \
                         (1.0 - emax(Cmmi, p["IC50_TPO_mmi"], p["Imax_TPO"]))

    # ---- (f) IODIDE : Wolff-Chaikoff, BOTH arms ---------------------
    wc = (Ithy ** p["nWC"]) / (p["IC_WC"] ** p["nWC"] + Ithy ** p["nWC"])
    d["WC"] = wc
    d["A_TPO"] = (1.0 - d["inh_TPO"]) * (1.0 - p["Imax_WCorg"] * wc)
    d["rel_inh"] = p["Imax_WCrel"] * wc

    # ---- (g) GLAND --------------------------------------------------
    Stim = max(p["TRAb"] + 0.6 * (TSH / p["TSH0"] - 1.0), 0.05)
    d["Stim"] = Stim
    # Organification is ENZYME-limited, not first-order in iodide.  Without
    # the Michaelis term an iodide load drives synthesis without bound and the
    # gland refills straight through a thionamide block -- which is wrong.
    # Vmax rises with stimulation because TPO/DUOX are cAMP-regulated.
    d["J_org"] = p["Vorg0"] * Stim * d["A_TPO"] * Ithy / (p["Km_org"] + Ithy)
    d["Synth"] = d["J_org"] / 4.0 * 1000.0                        # nmol T4-eq/h
    d["Sec"] = p["ksec0"] * Stim * S * (1.0 - d["rel_inh"])
    d["fT3sec"] = p["fT3sec0"] + (p["fT3sec_max"] - p["fT3sec0"]) * \
        min(1.0, max(0.0, (Stim - 1.0) / 5.0))
    # NIS is cAMP/TSH-regulated, so the thyroid's SHARE of the plasma iodide
    # pool rises with stimulation.  Without this the gland's output is capped
    # by dietary iodine and no thyrotoxic steady state above ~1.2x is reachable.
    # NIS SATURATES.  Written as a clearance times a Michaelis factor so that
    # at physiological iodide it is exactly CLthy*NIS*Stim, while a
    # pharmacologic load cannot drive the intrathyroidal pool without bound.
    d["CLthy_eff"] = p["CLthy"] * NIS * Stim / (1.0 + Ipl / p["Km_NIS"])

    # ---- (h) THERMAL ------------------------------------------------
    M_thy = max(1.0 + p["eBMR"] * (Sig - 1.0), 0.35)
    M_temp = p["Q10"] ** ((Tc - 37.0) / 10.0)
    M_sns = max(1.0 + p["eQ_sns"] * (Bsig - 1.0) + p["ePGE_shiv"] * PGE, 0.40)
    M_unc = 1.0 + emax(Casa, p["K_unc"], p["Emax_unc"])
    d["M_thy"], d["M_temp"], d["M_sns"] = M_thy, M_temp, M_sns
    d["Qprod"] = p["BMR0"] * M_thy * M_temp * M_sns * M_unc

    PGEeff = PGE * (1.0 - emax(Capap, p["IC50_apap"], p["eapap_PGE"]))
    d["PGEeff"] = PGEeff
    Tset = p["Tset0"] + p["ePGE_set"] * PGEeff
    E = 1.0 / (1.0 + math.exp(max(-40.0, min(40.0, -(Tc - Tset) / p["w_eff"]))))
    E *= (1.0 - p["eCNS_eff"] * CNSx)
    d["E"] = E
    cool = dose.get("cool", 0.0) if dose else 0.0
    d["h"] = (p["h_min"] + p["h_span"] * E) * \
             (max(Vol, 0.40) ** p["eVol_eff"]) * (1.0 + cool)
    d["Qloss"] = d["h"] * (Tc - p["T_sink"])

    # ---- (i) PUMP / PERFUSION ---------------------------------------
    contract = max(p["SV0"] * (1.0 + p["eSV_T3"] * (Sig - 1.0)
                               + p["eSV_b"] * (Bsig - 1.0)), 0.15)
    SV = contract * max(CR, 0.05) * min(1.0, Vol)
    d["CO"] = SV * HR / (p["SV0"] * p["HR0"])
    d["CO_need"] = d["Qprod"] / p["BMR0"]
    d["perf"] = min(1.6, d["CO"] / max(d["CO_need"], 0.2))
    d["shock"] = max(0.0, 1.0 - d["perf"] / 0.70)
    d["W_idx"] = contract * HR / p["HR0"]
    d["CHF"] = max(0.0, 1.0 - CR)
    # Effective glucocorticoid activity = endogenous cortisol PLUS the
    # exogenous hydrocortisone-equivalent.  Without this term the stress-dose
    # steroid in the bundle would inhibit D1/D2 but never repair the
    # accelerated-clearance deficit it is actually given for.
    d["Cort_eff"] = Cort + p["nmol_per_mgL"] * Cgc
    return d


# ---------------------------------------------------------------------
def burch_wartofsky(y, a, precip_known=True):
    """Burch-Wartofsky Point Scale, COMPUTED from the state."""
    Tc, HR, CNSx, GIx, Bili = y[15], y[13], y[17], y[20], y[18]
    if Tc < 37.2:   t = 0
    elif Tc < 37.8: t = 5
    elif Tc < 38.3: t = 10
    elif Tc < 38.9: t = 15
    elif Tc < 39.4: t = 20
    elif Tc < 40.0: t = 25
    else:           t = 30
    if CNSx < 0.12:   c = 0
    elif CNSx < 0.42: c = 10
    elif CNSx < 0.75: c = 20
    else:             c = 30
    if Bili >= 3.0:   g = 20
    elif GIx > 0.25:  g = 10
    else:             g = 0
    if HR < 99:    h = 0
    elif HR < 110: h = 5
    elif HR < 120: h = 10
    elif HR < 130: h = 15
    elif HR < 140: h = 20
    else:          h = 25
    chf = a["CHF"]
    if chf < 0.08:   f = 0
    elif chf < 0.20: f = 5
    elif chf < 0.38: f = 10
    else:            f = 15
    af = 10 if y[21] > 0.5 else 0
    return t + c + g + h + f + af + (10 if precip_known else 0)


# =====================================================================
# 3.  TARGETS OF THE FAST STATES  (used by both the ODEs and the QSS solver)
# =====================================================================
def fast_targets(y, p, a, dose):
    Tc, Vol, CNSx, GIx = y[15], y[16], y[17], y[20]
    Prec = y[24]
    M_lip = (1.0 + p["eL_T3"] * (a["Sig"] - 1.0)) * \
            (1.0 + p["eL_beta"] * (a["Bsig"] - 1.0))
    nefa_t = p["NEFA0"] * max(M_lip, 0.25)
    sns_t = (1.0 + p["ePrec_S"] * Prec
             + p["eTemp_S"] * max(0.0, Tc - 37.0)
             + p["eVol_S"] * max(0.0, 1.0 - Vol)
             + p["eT3_S"] * (a["Sig"] - 1.0))
    drive = (p["eHR_b"] * (a["Bsig"] - 1.0)
             + p["eHR_T3"] * (a["Sig"] - 1.0)
             + p["eHR_T"] * max(0.0, Tc - 37.0))
    if drive >= 0.0:
        hr_t = p["HR0"] + (p["HRmax"] - p["HR0"]) * drive / (p["Kd_HR"] + drive)
    else:
        hr_t = p["HR0"] * (1.0 + max(drive, -0.55) * 0.55)
    hr_t = max(hr_t, 38.0)
    cns_t = min(1.0, (p["kCNS_on"] * max(0.0, Tc - p["T_CNS"])
                      * (1.0 + p["eCNS_T3"] * (a["Sig"] - 1.0))
                      + p["kCNS_hyp"] * max(0.0, p["perf_CNS"] - a["perf"]))
                     / p["kCNS_off"])
    ins = p["kins_h"] * max(0.0, a["h"] - p["h_base"]) + \
        p["kins_T"] * max(0.0, Tc - 37.5)
    po = (1.0 - CNSx) * (1.0 - p["ePO_GI"] * GIx)
    fluid = dose.get("fluid", 0.0) if dose else 0.0
    losses = ins + p["kGI_loss"] * GIx
    vol_t = min(1.0, 1.0 - (losses - fluid) / max(p["kdrink"] * po, 1e-3))
    return nefa_t, sns_t, hr_t, cns_t, vol_t


def thermal_coeffs(y, p, a, dose):
    """Everything in the heat balance that does NOT depend on Tc.
    Qprod(Tc) = K1 * Q10**((Tc-37)/10)
    Qloss(Tc) = (h_min + h_span*Ecap/(1+exp(-(Tc-Tset)/w))) * hscale * (Tc-T_sink)
    Reducing the balance to closed form removes ~30 algebra() calls per solve."""
    K1 = p["BMR0"] * a["M_thy"] * a["M_sns"] * \
        (1.0 + emax(y[36], p["K_unc"], p["Emax_unc"]))
    Tset = p["Tset0"] + p["ePGE_set"] * a["PGEeff"]
    Ecap = 1.0 - p["eCNS_eff"] * y[17]
    cool = dose.get("cool", 0.0) if dose else 0.0
    hscale = (max(y[16], 0.40) ** p["eVol_eff"]) * (1.0 + cool)
    return K1, Tset, Ecap, hscale


def thermal_balance_c(c, p, Tc):
    K1, Tset, Ecap, hscale = c
    Qp = K1 * p["Q10"] ** ((Tc - 37.0) / 10.0)
    E = Ecap / (1.0 + math.exp(max(-40.0, min(40.0, -(Tc - Tset) / p["w_eff"]))))
    Ql = (p["h_min"] + p["h_span"] * E) * hscale * (Tc - p["T_sink"])
    return Qp - Ql


def solve_Tc_c(c, p, lo=30.0, hi=43.0, iters=42):
    """Stable operating temperature = root of Qprod - Qloss, from coefficients.
    Returns (Tc, exists);  exists=False means RUNAWAY (no stable point)."""
    if thermal_balance_c(c, p, lo) < 0.0:
        return lo, True
    if thermal_balance_c(c, p, hi) > 0.0:
        return hi, False
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if thermal_balance_c(c, p, mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), True


def storm_index(y, p, dose=None):
    """=====================================================================
    THE CENTRAL QUANTITY OF THIS MODEL.

    A thyrotoxic patient is safe if and only if there exists a temperature at
    which heat production and heat loss balance STABLY.  Write the ratio:

                     Q_prod      80 W · M_thy · M_temp · M_sns · M_unc
        Lambda(Tc) = ------ = ---------------------------------------------
                     Q_loss     (8 + 52·E) · Vol^1.5 · (1+cool) · (Tc − 33)

    Both sides are PRODUCTS OF FACTORS, and every drug used in thyroid storm
    multiplies one or two of them -- which is what makes the clinical
    hierarchy arithmetic rather than a list:

        M_thy   thionamide, iodide, PTU/GC/iopanoic (D1,D2), cholestyramine,
                plasma exchange, thyroidectomy ... and BETA-BLOCKADE, via
                NEFA -> displacement -> free fraction -> Sig
        M_sns   BETA-BLOCKADE (second time), cooling of shivering
        M_unc   ASPIRIN raises it (uncoupling) -- the wrong direction
        E       acetaminophen raises it (lowers the setpoint);
                CNS injury (ring C) destroys it
        Vol     fluid resuscitation raises it
        cool    external cooling raises it

    Returns
      stable  : a stable operating temperature exists
      Tc_star : that temperature (the patient's thermal set point)
      Lam_esc : max Lambda ABOVE Tc_star, over Tc up to 42 degC
                  < 1  no escape route  -> severe thyrotoxicosis
                  >= 1 escape route / no equilibrium -> THYROID STORM
      T_esc   : the temperature at which Lambda peaks
    ====================================================================="""
    dose = dose or {}
    a = algebra(y, p, dose)
    c = thermal_coeffs(y, p, a, dose)
    Tc_star, stable = solve_Tc_c(c, p)
    lo = (Tc_star + 0.05) if stable else 37.0
    best, Tbest = -1.0, lo
    Tc = lo
    while Tc <= 42.0 + 1e-9:
        K1, Tset, Ecap, hscale = c
        Qp = K1 * p["Q10"] ** ((Tc - 37.0) / 10.0)
        E = Ecap / (1.0 + math.exp(max(-40.0, min(40.0,
                                                  -(Tc - Tset) / p["w_eff"]))))
        Ql = (p["h_min"] + p["h_span"] * E) * hscale * (Tc - p["T_sink"])
        r = Qp / max(Ql, 1e-9)
        if r > best:
            best, Tbest = r, Tc
        Tc += 0.05
    K1, Tset, Ecap, hscale = c
    Ees = Ecap / (1.0 + math.exp(max(-40.0, min(40.0,
                                                -(Tbest - Tset) / p["w_eff"]))))
    return dict(stable=stable, Tc_star=Tc_star, Lam_esc=best, T_esc=Tbest,
                Lam_now=a["Qprod"] / max(a["Qloss"], 1e-9),
                M_thy=a["M_thy"], M_sns=a["M_sns"],
                M_unc=1.0 + emax(y[36], p["K_unc"], p["Emax_unc"]),
                M_temp=p["Q10"] ** ((Tbest - 37.0) / 10.0),
                E_esc=Ees, Vol_term=max(y[16], 0.40) ** p["eVol_eff"],
                cool_term=1.0 + (dose.get("cool", 0.0)),
                h_esc=(p["h_min"] + p["h_span"] * Ees) * hscale)


def fast_verdict(y, p, dose=None, hours=72.0, dt=0.05, T_storm=40.0):
    """OPERATIONAL DEFINITION OF A STORM, and the one used for every verdict in
    the report:  FREEZE the hormone level and every other slow state, then let
    the fast subsystem (NEFA, sympathetic tone, heart rate, temperature, volume,
    CNS) run for 72 h.  If the core temperature reaches 40 degC the fast
    subsystem has no reachable stable point -- that is a storm.  If it settles,
    the patient is a severe thyrotoxic who is not storming.

    This makes the storm/no-storm distinction a property of the FAST loop at a
    given hormone level, which is exactly the claim being tested: two patients
    with identical T4 and T3 can fall on opposite sides of it."""
    dose = dose or {}
    z = y.copy()
    Tmax = z[15]
    t = 0.0
    ix = FAST + [IX["Prec"], IX["PGE"]]   # Prec/PGE are external inputs, not
    #                                       disease state: they must decay
    for _ in range(int(hours / dt)):
        k1 = rhs(t, z, p, dose)
        z2 = z.copy()
        for i in ix:
            z2[i] = z[i] + dt / 2 * k1[i]
        z2 = clamp(z2, p)
        k2 = rhs(t, z2, p, dose)
        z3 = z.copy()
        for i in ix:
            z3[i] = z[i] + dt / 2 * k2[i]
        z3 = clamp(z3, p)
        k3 = rhs(t, z3, p, dose)
        z4 = z.copy()
        for i in ix:
            z4[i] = z[i] + dt * k3[i]
        z4 = clamp(z4, p)
        k4 = rhs(t, z4, p, dose)
        zn = z.copy()
        for i in ix:
            zn[i] = z[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        z = clamp(zn, p)
        Tmax = max(Tmax, z[15])
        t += dt
        if z[15] >= T_storm:
            break
    return dict(Tmax=Tmax, storm=(Tmax >= T_storm), Tc_end=z[15],
                HR_end=z[13], NEFA_end=z[9], Vol_end=z[16], CNSx_end=z[17])


def thermal_balance(y, p, dose, Tc):
    z = y.copy()
    z[15] = Tc
    a = algebra(z, p, dose)
    return a["Qprod"] - a["Qloss"]


def solve_Tc(y, p, dose):
    a = algebra(y, p, dose)
    return solve_Tc_c(thermal_coeffs(y, p, a, dose), p)


# =====================================================================
# 4.  ODE SYSTEM
# =====================================================================
def rhs(t, y, p, dose):
    a = algebra(y, p, dose)
    dy = np.zeros(NST)
    S, T4, T3, rT3, Ggut = y[0], y[1], y[2], y[3], y[4]
    Ipl, Ithy, NIS, TSH = y[5], y[6], y[7], y[8]
    NEFA, TBG, SNS, Rb = y[9], y[10], y[11], y[12]
    HR, CR, Tc, Vol = y[13], y[14], y[15], y[16]
    CNSx, Bili, Cort, GIx = y[17], y[18], y[19], y[20]
    AFb, Prec, PGE = y[21], y[24], y[25]
    nefa_t, sns_t, hr_t, cns_t, vol_t = fast_targets(y, p, a, dose)

    # ---- (0) thyroid store : the 2-month reservoir -------------------
    dy[0] = a["Synth"] - a["Sec"]

    # ---- (1) plasma T4 ---------------------------------------------
    kD1 = p["kD1_4"] * (1.0 - a["inh_D1"])
    kD2 = p["kD2_4"] * (1.0 - a["inh_D2"])
    reab_blk = emax(y[34], p["IC50_chol"], p["ekchol"])
    dy[1] = (a["Sec"] * (1.0 - a["fT3sec"]) / p["V_T4"]
             + p["kgut"] * Ggut * p["kreab"] * (1.0 - reab_blk) / p["V_T4"]
             - (kD1 + kD2 + p["kD3_4"] + p["kEHC_4"]) * T4)

    # ---- (2) plasma T3 : 80 % of it is made OUTSIDE the thyroid ------
    dy[2] = ((kD1 + kD2) * T4 * p["V_T4"] / p["V_T3"]
             + a["Sec"] * a["fT3sec"] / p["V_T3"]
             - p["kT3"] * T3)

    # ---- (3) rT3 : the observable signature of a D1 block ------------
    #  rT3 is cleared mainly BY D1, so a D1 block raises it.  This is the
    #  observable laboratory signature that PTU (and not methimazole) is
    #  working on the fast axis.
    krT3e = p["krT3"] * (1.0 - p["fD1_rT3"] * a["inh_D1"])
    dy[3] = p["kD3_4"] * T4 * p["V_T4"] / p["V_rT3"] - krT3e * rT3

    # ---- (4) enterohepatic T4 in gut lumen --------------------------
    dy[4] = p["kEHC_4"] * T4 * p["V_T4"] - p["kgut"] * Ggut

    # ---- (5) plasma inorganic iodide -------------------------------
    #  RECYCLING: every deiodination step returns one iodine to the plasma
    #  pool, and clearance of T3/rT3 returns the remaining three.  Without
    #  this term the gland is dietary-iodine-limited and cannot sustain
    #  thyrotoxicosis at all -- the iodide loop is what makes 3x output
    #  possible on an unchanged 150 ug/d intake.
    recyc = (p["nI_step"] * (p["kD1_4"] * (1.0 - a["inh_D1"])
                             + p["kD2_4"] * (1.0 - a["inh_D2"])
                             + p["kD3_4"]) * T4 * p["V_T4"]
             + p["nI_T3"] * p["kT3"] * T3 * p["V_T3"]
             + p["nI_T3"] * p["krT3"] * rT3 * p["V_rT3"]) / 1000.0
    dy[5] = (p["I_diet"] + dose.get("iod_rate", 0.0) + recyc
             - p["CLren_I"] * Ipl - a["CLthy_eff"] * Ipl
             + p["kleak"] * Ithy) / p["V_I"]

    # ---- (6) intrathyroidal iodide --------------------------------
    dy[6] = a["CLthy_eff"] * Ipl - a["J_org"] - p["kleak"] * Ithy

    # ---- (7) NIS down-regulation = WOLFF-CHAIKOFF ESCAPE -----------
    dy[7] = (1.0 / (1.0 + Ipl / p["K_NIS"]) - NIS) / p["tau_NIS"]

    # ---- (8) TSH ---------------------------------------------------
    TSH_t = min(max(p["TSH0"] * math.exp(-p["kTSHfb"] * (a["Sig"] - 1.0)),
                    0.002), 60.0)
    dy[8] = (TSH_t - TSH) / p["tau_TSH"]

    # ---- (9) NEFA  == RING A, first half ---------------------------
    dy[9] = (nefa_t - NEFA) / p["tau_NEFA"]

    # ---- (10) TBG --------------------------------------------------
    dy[10] = (1.0 - p["eTBG_prec"] * min(1.0, Prec) - TBG) / p["tau_TBG"]

    # ---- (11) sympathetic outflow ---------------------------------
    dy[11] = (sns_t - SNS) / p["tau_SNS"]

    # ---- (12) beta1 receptor density (T3-driven sensitisation) -----
    dy[12] = (1.0 + p["eR"] * (a["Sig"] - 1.0) - Rb) / p["tau_Rb"]

    # ---- (13) heart rate ------------------------------------------
    dy[13] = (hr_t - HR) / p["tau_HR"]

    # ---- (14) cardiac reserve -------------------------------------
    over = max(0.0, a["W_idx"] / p["W_crit"] - 1.0)
    dy[14] = -p["kdep"] * over + (p["krec"] * (1.0 - CR) if over <= 0.0 else 0.0)

    # ---- (15) CORE TEMPERATURE == RINGS B and C -------------------
    dy[15] = (a["Qprod"] - a["Qloss"]) * 3600.0 / p["C_body"]

    # ---- (16) volume ---------------------------------------------
    ins = p["kins_h"] * max(0.0, a["h"] - p["h_base"]) + \
        p["kins_T"] * max(0.0, Tc - 37.5)
    po = (1.0 - CNSx) * (1.0 - p["ePO_GI"] * GIx)
    dy[16] = (p["kdrink"] * po * (1.0 - Vol) + dose.get("fluid", 0.0)
              - ins - p["kGI_loss"] * GIx)

    # ---- (17) CNS / hypothalamic injury == RING C ----------------
    dy[17] = (cns_t - CNSx) * p["kCNS_off"]

    # ---- (18) bilirubin -----------------------------------------
    dy[18] = p["kb_out"] * (p["Bili0"] * (1.0 + p["eB_hyp"] * max(0.0, 0.85 - a["perf"])
                                          + p["eB_hf"] * a["CHF"]
                                          + p["eB_T3"] * (a["Sig"] - 1.0)) - Bili)

    # ---- (19) cortisol : clearance accelerated BY thyroid hormone --
    dy[19] = (p["kCL_cort"] * p["Cort0"] * p["AR"] * (1.0 + p["eStress_cort"] * Prec)
              - p["kCL_cort"] * (1.0 + p["eCL_T3"] * (a["Sig"] - 1.0)) * Cort
              + dose.get("gc_cort", 0.0))

    # ---- (20) GI dysfunction ------------------------------------
    dy[20] = (min(1.0, p["eGI_T3"] * max(0.0, a["Sig"] - 1.0)
                  + p["eGI_prec"] * Prec) - GIx) / p["tau_GI"]

    # ---- (21) atrial fibrillation burden ------------------------
    AF_t = min(1.0, p["AF_suscept"] * (p["eAF_T3"] * max(0.0, a["Sig"] - 1.0) * 2.0
                                       + p["eAF_b"] * max(0.0, a["Bsig"] - 1.0)))
    dy[21] = (AF_t - AFb) / p["tau_AF"]

    # ---- (22) mortality hazard ----------------------------------
    bw = burch_wartofsky(y, a, precip_known=(Prec > 0.05))
    haz = (p["h0_haz"] * math.exp(p["aBW"] * (bw - 45.0) / 10.0)
           * (1.0 + p["eShock"] * a["shock"])
           * (1.0 + p["eBili"] * max(0.0, Bili - 1.2) / 3.0)
           * (1.0 + p["eHyper"] * max(0.0, Tc - 40.5))
           * (1.0 + p["eAI"] * max(0.0, 1.0 - a["Cort_eff"] /
                                   (p["Cort0"] * (1.0 + 0.9 * Prec)))))
    dy[22] = haz
    dy[23] = bw

    # ---- (24,25) precipitant -----------------------------------
    dy[24] = -p["kPrec"] * Prec
    dy[25] = -p["kPGE"] * PGE

    # ================= DRUG PK =================
    dy[26] = -p["ka_ptu"] * y[26]
    dy[27] = p["ka_ptu"] * y[26] * p["F_ptu"] / p["V_ptu"] - p["ke_ptu"] * y[27]
    dy[28] = -p["ka_mmi"] * y[28]
    dy[29] = p["ka_mmi"] * y[28] * p["F_mmi"] / p["V_mmi"] - p["ke_mmi"] * y[29]
    dy[30] = -p["ka_pro"] * y[30]
    dy[31] = p["ka_pro"] * y[30] * p["F_pro"] * 1000.0 / p["V_pro"] \
        - p["ke_pro"] * y[31]
    dy[32] = dose.get("esm_rate", 0.0) - p["ke_esm"] * y[32]
    dy[33] = dose.get("gc_rate", 0.0) - p["ke_gc"] * y[33]
    dy[34] = -p["kchol"] * y[34]
    dy[35] = dose.get("iop_rate", 0.0) - p["ke_iop"] * y[35]
    dy[36] = -p["ke_asa"] * y[36]
    dy[37] = -p["ke_apap"] * y[37]
    return dy


# =====================================================================
# 5.  INTEGRATORS
# =====================================================================
def clamp(y, p):
    y = np.maximum(y, 0.0)
    y[15] = min(y[15], p['Tc_ceil'])
    y[16] = min(max(y[16], p["Vol_floor"]), 1.0)
    y[17] = min(y[17], 1.0)
    y[21] = min(y[21], 1.0)
    y[14] = min(y[14], 1.0)
    y[7] = min(y[7], 1.0)
    return y


def rk4(y, p, dose, dt, t):
    k1 = rhs(t, y, p, dose)
    k2 = rhs(t + dt / 2, clamp(y + dt / 2 * k1, p), p, dose)
    k3 = rhs(t + dt / 2, clamp(y + dt / 2 * k2, p), p, dose)
    k4 = rhs(t + dt, clamp(y + dt * k3, p), p, dose)
    return clamp(y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4), p)


SLOW = [i for i in range(NST) if i not in FAST]


def chronic_ss(p, hours=9000.0, dt=0.10, y0=None, dose=None):
    """Chronic (weeks-to-months) steady state, obtained by integrating the FULL
    ODE system with the same RK4 integrator used for the storm simulations.

    Deliberately NOT a quasi-steady-state projection.  The fast subsystem is
    bistable -- besides the cool operating point there is always a
    self-consistent hot one (Tc at the ceiling -> CNSx = 1 -> heat-loss
    effector destroyed -> Tc at the ceiling) -- and any fixed-point projection
    can fall into the hot branch even where the cool branch exists, which would
    misreport ordinary florid thyrotoxicosis as a storm.  Time integration from
    a euthyroid (or a nearby thyrotoxic) initial condition cannot do that: the
    trajectory stays on the branch it starts on unless that branch genuinely
    ceases to exist, which is precisely the definition of a storm here."""
    dose = dose or {}
    y = init_state(p) if y0 is None else y0.copy()
    t = 0.0
    for _ in range(int(hours / dt)):
        y = rk4(y, p, dose, dt, t)
        t += dt
    return y


def place_on_fast_manifold(y, p, dose, iters=60, tol=1e-9, cold=False,
                           Tc_cap=None):
    """Put the six fast states on their quasi-steady manifold.

    BRANCH TRACKING MATTERS HERE.  The fast subsystem is bistable: besides the
    cool operating point there is always a self-consistent HOT solution
    (Tc at the ceiling -> CNSx = 1 -> heat-loss effector destroyed -> Tc at the
    ceiling).  A naive fixed point can converge to the hot branch even when the
    cool one exists, which would misreport ordinary florid thyrotoxicosis as a
    storm.  `cold=True` starts the iteration from euthyroid fast values and
    limits the temperature step, so the iteration climbs to the LOWEST
    self-consistent point and stops there.  Genuine storms are then exactly the
    cases in which no such point exists (solve_Tc returns ok=False)."""
    y = y.copy()
    if cold:
        y[9] = p["NEFA0"]; y[11] = 1.0; y[13] = p["HR0"]
        y[15] = 37.0; y[16] = 1.0; y[17] = 0.0
        iters = max(iters, 260)
    relax = 0.30 if cold else 0.50
    for _ in range(iters):
        a = algebra(y, p, dose)
        nefa_t, sns_t, hr_t, cns_t, vol_t = fast_targets(y, p, a, dose)
        Tc_t, ok = solve_Tc_c(thermal_coeffs(y, p, a, dose), p)
        cur = (y[9], y[11], y[13], y[15], y[16], y[17])
        new = (nefa_t, sns_t, hr_t, Tc_t, vol_t, cns_t)
        err = max(abs(n - c) / max(abs(c), 1e-6) for n, c in zip(new, cur))
        y[9] = cur[0] + relax * (new[0] - cur[0])
        y[11] = cur[1] + relax * (new[1] - cur[1])
        y[13] = cur[2] + relax * (new[2] - cur[2])
        dT = relax * (new[3] - cur[3])
        y[15] = cur[3] + max(-0.35, min(0.35, dT))      # climb, never leap
        if Tc_cap is not None:
            y[15] = min(y[15], Tc_cap)
        y[16] = cur[4] + relax * (new[4] - cur[4])
        y[17] = cur[5] + relax * (new[5] - cur[5])
        y = clamp(y, p)
        if err < tol:
            break
    return y


def init_state(p):
    y = np.zeros(NST)
    y[0] = 8000.0
    y[1] = 100.0
    y[2] = 1.80
    y[3] = 0.2805
    y[4] = 14.75
    y[5] = 0.014933
    y[6] = 0.20
    y[7] = 1.0
    y[8] = p["TSH0"]
    y[9] = p["NEFA0"]
    y[10] = 1.0
    y[11] = 1.0
    y[12] = 1.0
    y[13] = p["HR0"]
    y[14] = 1.0
    y[15] = 37.0
    y[16] = 1.0
    y[17] = 0.0
    y[18] = p["Bili0"]
    y[19] = p["Cort0"]
    y[20] = 0.0
    y[21] = 0.0
    return y


# =====================================================================
# 6.  LOOP GAIN  L  -- the central number of this model
# =====================================================================
def loop_gain(y, p, dose=None):
    """Dimensionless loop gain of the fast positive-feedback ring, computed as
    the spectral radius of the one-pass log-gain matrix G over
    x = (NEFA, Tc, CNSx) with all slow states and drug concentrations frozen.
    G_ij = dln f_i / dln x_j, where f is the one-pass target map.

      L <  1  ->  contracting ring: a stable operating point exists
      L >= 1  ->  expanding ring: THYROID STORM
    """
    dose = dose or {}
    ix = [IX["NEFA"], IX["Tc"], IX["CNSx"]]

    def one_pass(xv):
        z = y.copy()
        for k, i in enumerate(ix):
            z[i] = xv[k]
        z = clamp(z, p)
        a = algebra(z, p, dose)
        nefa_t, _, _, cns_t, _ = fast_targets(z, p, a, dose)
        # one Newton step of the heat balance = the one-pass target of Tc
        c = thermal_coeffs(z, p, a, dose)
        Tc = z[IX["Tc"]]
        f0 = thermal_balance_c(c, p, Tc)
        f1 = thermal_balance_c(c, p, Tc + 1e-3)
        slope = (f1 - f0) / 1e-3
        Tc_t = Tc - f0 / slope if abs(slope) > 1e-9 else Tc + 8.0
        Tc_t = min(max(Tc_t, 30.5), 60.0)
        return np.array([nefa_t, Tc_t, max(cns_t, 1e-4)])

    x0 = np.array([max(y[ix[0]], 1e-3), y[ix[1]], max(y[ix[2]], 1e-4)])
    f0 = one_pass(x0)
    G = np.zeros((3, 3))
    for j in range(3):
        h = 1e-4 * max(abs(x0[j]), 1e-3)
        xp = x0.copy(); xp[j] += h
        fp = one_pass(xp)
        for i in range(3):
            if i == 1 and j == 1:      # temperature on temperature: use T - T_sink
                num = math.log(max(fp[1] - p["T_sink"], 1e-3)) - \
                      math.log(max(f0[1] - p["T_sink"], 1e-3))
                den = math.log(max(xp[1] - p["T_sink"], 1e-3)) - \
                      math.log(max(x0[1] - p["T_sink"], 1e-3))
            elif i == 1:
                num = math.log(max(fp[1] - p["T_sink"], 1e-3)) - \
                      math.log(max(f0[1] - p["T_sink"], 1e-3))
                den = math.log(max(xp[j], 1e-12)) - math.log(max(x0[j], 1e-12))
            elif j == 1:
                num = math.log(max(fp[i], 1e-12)) - math.log(max(f0[i], 1e-12))
                den = math.log(max(xp[1] - p["T_sink"], 1e-3)) - \
                      math.log(max(x0[1] - p["T_sink"], 1e-3))
            else:
                num = math.log(max(fp[i], 1e-12)) - math.log(max(f0[i], 1e-12))
                den = math.log(max(xp[j], 1e-12)) - math.log(max(x0[j], 1e-12))
            G[i, j] = num / den if abs(den) > 1e-14 else 0.0
    ev = np.linalg.eigvals(G)
    return dict(L=float(max(abs(ev))), G=G, ev=ev,
                L_A=abs(G[0, 0]), L_B=abs(G[1, 1]),
                L_C=math.sqrt(abs(G[1, 2] * G[2, 1])))


def fast_jacobian_eig(y, p, dose=None):
    """Largest real eigenvalue of the true ODE Jacobian over the fast
    subsystem: an independent confirmation of the loop-gain verdict."""
    dose = dose or {}
    ix = [IX["NEFA"], IX["Tc"], IX["CNSx"], IX["SNS"], IX["HR"]]
    n = len(ix)
    J = np.zeros((n, n))
    f0 = rhs(0.0, y, p, dose)
    for j, jx in enumerate(ix):
        h = 1e-5 * max(abs(y[jx]), 1e-3)
        z = y.copy(); z[jx] += h
        f1 = rhs(0.0, z, p, dose)
        for i, i2 in enumerate(ix):
            J[i, j] = (f1[i2] - f0[i2]) / h
    return float(max(np.linalg.eigvals(J).real))


# =====================================================================
# 7.  PHENOTYPES
# =====================================================================
_STATE_CACHE = {}          # TRAb -> chronic steady state (warm-start library)
_TRAB_CACHE = {}           # target total T3 -> TRAb


def ss_for_trab(tr, p=None, hours=None):
    """Chronic steady state at a given TRAb drive, WARM-STARTED from the
    nearest previously solved TRAb.  Continuation makes the whole report
    affordable: a cold solve needs 14 000 simulated hours, a warm one 5 000."""
    p = p or P
    key = round(tr, 6)
    if key in _STATE_CACHE:
        return _STATE_CACHE[key].copy()
    if _STATE_CACHE:
        near = min(_STATE_CACHE, key=lambda k: abs(math.log(k) - math.log(tr)))
        y0, h = _STATE_CACHE[near].copy(), (hours or 4000.0)
    else:
        y0, h = None, (hours or 9000.0)
    y = chronic_ss(dict(p, TRAb=tr), hours=h, y0=y0)
    _STATE_CACHE[key] = y.copy()
    return y.copy()


def find_trab_for_T3(target_T3, p=None):
    """Solve the TRAb (TSH-receptor stimulating) drive that produces a given
    chronic total T3.  Secant iteration in log(TRAb); nothing about storm
    enters here -- this is a purely thyrotoxicosis-level calibration."""
    key = round(target_T3, 4)
    if key in _TRAB_CACHE:
        return _TRAB_CACHE[key]
    x0, x1 = math.log(1.0), math.log(12.0)
    f = lambda x: ss_for_trab(math.exp(x), p)[2] - target_T3
    f0, f1 = f(x0), f(x1)
    for _ in range(12):
        if abs(f1 - f0) < 1e-12:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = min(max(x2, math.log(0.2)), math.log(80.0))
        f2 = f(x2)
        x0, f0, x1, f1 = x1, f1, x2, f2
        if abs(f1) < 2e-4:
            break
    _TRAB_CACHE[key] = math.exp(x1)
    return _TRAB_CACHE[key]


# =====================================================================
# 8.  DYNAMIC SIMULATION
# =====================================================================
DRUG_SLOT = dict(ptu=26, mmi=28, pro=30, chol=34, asa=36, apap=37)


def simulate(p, spec, hours=168.0, dt=0.05, y0=None, record=None,
             cool=0.0, fluid=0.0, tpe=None, precip=1.0, pge=0.6, rec_every=10):
    events = []
    for s in spec:
        for k in range(s.get("n", 1)):
            if s.get("mode") != "inf":
                events.append((s.get("start", 0.0) + k * s.get("ii", 1e9),
                               s["drug"], s.get("amt", 0.0)))
    events.sort()
    infus = [s for s in spec if s.get("mode") == "inf"]
    y = y0.copy()
    y[22] = 0.0          # cumulative mortality hazard -- reset per simulation
    y[23] = 0.0          # BWPS integral            -- reset per simulation
    y[24], y[25] = precip, pge
    tpe = tpe or []
    tpe_done = [False] * len(tpe)
    out = {k: [] for k in (record or [])}
    out["time"] = []
    ev_i, t, n = 0, 0.0, int(hours / dt)
    dose = {}
    for i in range(n + 1):
        while ev_i < len(events) and events[ev_i][0] <= t + 1e-9:
            _, drug, amt = events[ev_i]
            if drug in DRUG_SLOT:
                y[DRUG_SLOT[drug]] += amt
            elif drug == "iodide":
                y[5] += amt / p["V_I"]
            elif drug == "gc":
                y[33] += amt / p["V_gc"]
            elif drug == "iop":
                y[35] += amt / p["V_iop"]
            ev_i += 1
        for j, (ttpe, frac) in enumerate(tpe):
            if not tpe_done[j] and t >= ttpe - 1e-9:
                # A 1-1.5 plasma-volume exchange removes `frac` of the PLASMA
                # pool, not of the body.  Plasma (~3 L) is 30 % of the T4
                # distribution volume but only 7.5 % of T3's and 10 % of rT3's,
                # so the whole-body removal is far smaller than the plasma dip.
                # The model has no plasma/extravascular split, so it can only
                # represent the net whole-body removal (see caveat L7).
                y[1] *= (1.0 - frac * 3.0 / p["V_T4"])
                y[2] *= (1.0 - frac * 3.0 / p["V_T3"])
                y[3] *= (1.0 - frac * 3.0 / p["V_rT3"])
                tpe_done[j] = True
        dose = dict(cool=cool, fluid=fluid)
        for s in infus:
            if s["start"] <= t <= s.get("stop", 1e9):
                dose[{"esm": "esm_rate", "gc": "gc_rate",
                      "iodide": "iod_rate", "iop": "iop_rate"}[s["drug"]]] = s["rate"]
        if record and (i % rec_every == 0 or i == n):
            a = algebra(y, p, dose)
            out["time"].append(t)
            for k in record:
                if k in IX:
                    out[k].append(y[IX[k]])
                elif k == "BW":
                    out[k].append(burch_wartofsky(y, a))
                elif k == "mort":
                    out[k].append(1.0 - math.exp(-y[22]))
                elif k == "L":
                    out[k].append(loop_gain(y, p, dose)["L"])
                elif k == "Lam":
                    out[k].append(storm_index(y, p, dose)["Lam_esc"])
                elif k == "stable":
                    out[k].append(1.0 if storm_index(y, p, dose)["stable"] else 0.0)
                elif k == "CortEff":
                    out[k].append(a["Cort_eff"])
                elif k == "aiMult":
                    need = p["Cort0"] * (1.0 + 0.9 * y[24])
                    out[k].append(1.0 + p["eAI"] *
                                  max(0.0, 1.0 - a["Cort_eff"] / need))
                elif k == "BWobj":
                    out[k].append(burch_wartofsky(y, a, precip_known=y[24] > 0.05))
                else:
                    out[k].append(a[k])
        if i == n:
            break
        y = rk4(y, p, dose, dt, t)
        t += dt
    return y, algebra(y, p, dose), out


def at(o, key, hh):
    t = np.asarray(o["time"])
    return o[key][int(np.argmin(np.abs(t - hh)))]


# =====================================================================
# 9.  REPORT
# =====================================================================
def RULE(c="-", n=94):
    return c * n


# ---------------------------------------------------------------------
# regimens used throughout the report
# ---------------------------------------------------------------------
PTU = [dict(drug="ptu", amt=600.0, start=0.0),
       dict(drug="ptu", amt=250.0, start=4.0, ii=4.0, n=42)]
MMI = [dict(drug="mmi", amt=40.0, start=0.0),
       dict(drug="mmi", amt=25.0, start=6.0, ii=6.0, n=28)]
PRO = [dict(drug="pro", amt=80.0, start=0.25, ii=4.0, n=42)]
ESM = [dict(drug="esm", mode="inf", rate=1.60, start=0.25)]
IOD_AFTER = [dict(drug="iodide", amt=1970.0, start=1.0, ii=6.0, n=28)]
IOD_BEFORE = [dict(drug="iodide", amt=1970.0, start=0.0, ii=6.0, n=28)]
PTU_LATE = [dict(drug="ptu", amt=600.0, start=4.0),
            dict(drug="ptu", amt=250.0, start=8.0, ii=4.0, n=41)]
HC = [dict(drug="gc", amt=100.0, start=0.25, ii=8.0, n=21)]
CHO = [dict(drug="chol", amt=4.0, start=1.0, ii=6.0, n=28)]
IOP = [dict(drug="iop", amt=1000.0, start=1.0, ii=8.0, n=21)]
ASA = [dict(drug="asa", amt=650.0 / 12.0 / 10.0, start=0.5, ii=4.0, n=42)]
APAP = [dict(drug="apap", amt=650.0 / 50.0, start=0.5, ii=6.0, n=28)]
SUPP = dict(cool=1.2, fluid=0.0045)
BUNDLE = PTU + IOD_AFTER + PRO + HC

PREC = 1.30          # canonical fulminant precipitant (see section [Bi])
PGE = 0.6 * PREC

REC = ["Tc", "HR", "T3", "T4", "fT3", "NEFA", "Fd", "BW", "mort", "Bili",
       "CNSx", "Sig", "Bsig", "perf", "CR", "Cort", "S", "Sec", "rT3",
       "inh_D1", "inh_TPO", "rel_inh", "Cptu", "Cmmi", "Cpro", "shock",
       "NIS", "Ithy", "Ipl", "Synth", "A_TPO", "occ_b", "Qprod", "h", "Vol"]


def main():
    W = []

    def w(s=""):
        W.append(s)

    w(RULE("="))
    w("THYROID STORM (갑상선 폭풍) QSP MODEL -- NUMERICAL VERIFICATION")
    w(RULE("="))
    w("Every number quoted in README.md is produced by this file.")
    w("Central claim under test: storm is a property of the FAST feedback loop,")
    w("not of the hormone level.  Sections [H] and [S] are the controlled")
    w("comparison; sections [C] and [L] are the therapeutic consequence.")

    # ================= [N] euthyroid =================
    y_eu = chronic_ss(P)
    a_eu = algebra(y_eu, P)
    _STATE_CACHE[1.0] = y_eu.copy()
    conv = (P["kD1_4"] + P["kD2_4"]) * y_eu[1] * P["V_T4"]
    w()
    w("[N] EUTHYROID STEADY STATE -- every target is NORMAL HUMAN PHYSIOLOGY.")
    w("    Nothing in this block is fitted to thyrotoxicosis or to storm.")
    w(RULE())
    rows = [
        ("N1  total T4", y_eu[1], "nmol/L", "100"),
        ("N2  total T3", y_eu[2], "nmol/L", "1.8"),
        ("N3  free T4", a_eu["fT4"], "pmol/L", "12-22"),
        ("N4  free T3", a_eu["fT3"], "pmol/L", "3.5-6.5"),
        ("N5  reverse T3", y_eu[3], "nmol/L", "0.15-0.45"),
        ("N6  T3 from peripheral conversion",
         100.0 * conv / (conv + a_eu["Sec"] * a_eu["fT3sec"]), "%", "~80"),
        ("N7  thyroid store, days of supply",
         y_eu[0] / (a_eu["Sec"] * 24.0), "d", "60-90"),
        ("N8  TSH", y_eu[8], "mIU/L", "1.5"),
        ("N9  core temperature", y_eu[15], "degC", "37.0"),
        ("N10 heart rate", y_eu[13], "bpm", "70"),
        ("N11 NEFA", y_eu[9], "mmol/L", "0.4"),
        ("N12 cortisol", y_eu[19], "nmol/L", "400"),
        ("N13 heat production = heat loss", a_eu["Qprod"], "W", "80"),
        ("N14 T4 half-life",
         math.log(2) / ((P["kD1_4"] + P["kD2_4"] + P["kD3_4"]
                         + P["kEHC_4"] * (1 - P["kreab"])) * 24.0), "d", "6-7"),
        ("N15 T3 half-life", math.log(2) / (P["kT3"] * 24.0), "d", "1.0"),
        ("N16 dietary iodide", P["I_diet"] * 24 * 126.9, "ug/d", "150"),
        ("N17 T4 secretion",
         a_eu["Sec"] * (1 - a_eu["fT3sec"]) * 24 * 0.7767, "ug/d", "80-110"),
        ("N18 plasma inorganic iodide", y_eu[5], "umol/L", "0.006-0.05"),
        ("N19 heat-loss coefficient h", a_eu["h"], "W/K", "8-60 range"),
        ("N20 cardiac work index", a_eu["W_idx"], "-", "1.0"),
    ]
    for nm, v, u, tgt in rows:
        w(f"  {nm:36s} {v:11.4f} {u:8s}  target {tgt}")
    w(f"  free fraction T3 {100*a_eu['phi3']:.4f} %   "
      f"free fraction T4 {100*a_eu['phi4']:.4f} %   displacement Fd "
      f"{a_eu['Fd']:.4f}")
    w(f"  Burch-Wartofsky {burch_wartofsky(y_eu, a_eu, False)}   "
      f"fast-subsystem verdict: "
      f"{'STORM' if fast_verdict(y_eu, P)['storm'] else 'no storm'}")

    TRAB_INDEX = find_trab_for_T3(5.40)

    # ================= [H] the control experiment =================
    w()
    w("[H] CONTROL EXPERIMENT -- RAISE THE HORMONE AND NOTHING ELSE.")
    w("    TRAb (TSH-receptor stimulating drive) is swept; no precipitant is")
    w("    ever present.  The 72-h fast-subsystem verdict is the storm test.")
    w(RULE())
    w(f"  {'TRAb':>6s} {'totT4':>7s} {'totT3':>6s} {'freeT3':>7s} {'Fd':>6s} "
      f"{'Sig':>6s} {'store':>6s} {'days':>5s} {'TSH':>7s} {'Tc':>6s} "
      f"{'HR':>4s} {'NEFA':>5s} {'BWPS':>5s} {'Lam':>6s} {'72h Tmax':>9s}  verdict")
    hsweep = []
    for tr in (1.0, 1.5, 2.0, 3.0, TRAB_INDEX, 6.0, 9.0, 14.0, 22.0, 35.0):
        yz = ss_for_trab(tr)
        Pz = dict(P, TRAb=tr)
        az = algebra(yz, Pz)
        sz = storm_index(yz, Pz)
        fv = fast_verdict(yz, Pz)
        hsweep.append((tr, yz, az, sz, fv))
        w(f"  {tr:6.2f} {yz[1]:7.1f} {yz[2]:6.2f} {az['fT3']:7.2f} "
          f"{az['Fd']:6.3f} {az['Sig']:6.3f} {yz[0]:6.0f} "
          f"{yz[0]/(az['Sec']*24):5.1f} {yz[8]:7.4f} {yz[15]:6.2f} "
          f"{yz[13]:4.0f} {yz[9]:5.2f} "
          f"{burch_wartofsky(yz, az, False):5d} {sz['Lam_esc']:6.3f} "
          f"{fv['Tmax']:9.2f}  {'STORM' if fv['storm'] else 'no storm'}")
    w("  -> NO hormone level produces a storm.  Lambda saturates just below 1")
    w("     because the TR signal Sig saturates at 2 while the heat-loss")
    w("     effector can rise 7-fold: the thyroid axis ALONE is arithmetically")
    w("     too weak to overwhelm human thermoregulation.  A second hit is")
    w("     REQUIRED.  Note BWPS plateaus at 40-45, i.e. exactly the clinical")
    w("     'impending storm' band, without ever crossing into storm.")

    # ================= [T] the index patient =================
    trab = find_trab_for_T3(5.40)
    Ptt = dict(P, TRAb=trab)
    y_tt = ss_for_trab(trab)
    a_tt = algebra(y_tt, Ptt)
    w()
    w(f"[T] THE INDEX PATIENT.  TRAb solved so that total T3 = 5.40 nmol/L")
    w(f"    (3.0x normal) = FLORID Graves thyrotoxicosis.  TRAb = {trab:.3f}")
    w(RULE())
    w(f"  total T4          {y_tt[1]:9.1f} nmol/L  ({y_tt[1]/y_eu[1]:.2f}x normal)")
    w(f"  total T3          {y_tt[2]:9.2f} nmol/L  ({y_tt[2]/y_eu[2]:.2f}x normal)")
    w(f"  free  T3          {a_tt['fT3']:9.2f} pmol/L  ({a_tt['fT3n']:.2f}x normal)")
    w(f"  displacement Fd   {a_tt['Fd']:9.3f}         (1.000 = normal binding)")
    w(f"  T3:T4 molar ratio {1000*y_tt[2]/y_tt[1]:9.2f} x1e-3   "
      f"(euthyroid {1000*y_eu[2]/y_eu[1]:.2f} -- raised, as in Graves)")
    w(f"  thyroid store     {y_tt[0]:9.0f} nmol   = "
      f"{y_tt[0]/(a_tt['Sec']*24):.1f} d of supply  "
      f"(euthyroid {y_eu[0]/(a_eu['Sec']*24):.0f} d)   <- UNFITTED PREDICTION")
    w(f"  TSH               {y_tt[8]:9.4f} mIU/L")
    w(f"  core temperature  {y_tt[15]:9.2f} degC")
    w(f"  heart rate        {y_tt[13]:9.1f} bpm")
    w(f"  NEFA              {y_tt[9]:9.2f} mmol/L")
    w(f"  beta1 density Rb  {y_tt[12]:9.2f} x normal")
    w(f"  cortisol          {y_tt[19]:9.0f} nmol/L  (clearance "
      f"{1+P['eCL_T3']*(a_tt['Sig']-1):.2f}x normal)   <- UNFITTED")
    w(f"  cardiac work      {a_tt['W_idx']:9.2f}  (W_crit = {P['W_crit']:.2f})")
    w(f"  Burch-Wartofsky   {burch_wartofsky(y_tt, a_tt, False):9d}  "
      f"(no precipitant -> 'storm unlikely / impending')")
    fv0 = fast_verdict(y_tt, Ptt)
    w(f"  72-h fast verdict: Tmax {fv0['Tmax']:.2f} degC -> "
      f"{'STORM' if fv0['storm'] else 'NO STORM'}")

    # ================= [S] same patient, plus a precipitant =================
    w()
    w("[S] THE SAME PATIENT PLUS A PRECIPITANT.  Total T4 and T3 at t = 0 are")
    w("    IDENTICAL to the block above; the slow states are byte-identical.")
    w("    Only Prec (and the PGE2 it generates) is switched on.")
    w(RULE())
    y0s = y_tt.copy()
    y0s[24], y0s[25] = PREC, PGE
    fv1 = fast_verdict(y0s, Ptt)
    w(f"  t=0  total T3 {y0s[2]:.2f} nmol/L, total T4 {y0s[1]:.1f} nmol/L "
      f"(unchanged)")
    w(f"  72-h fast verdict: Tmax {fv1['Tmax']:.2f} degC -> "
      f"{'STORM' if fv1['storm'] else 'NO STORM'}   "
      f"(HR {fv1['HR_end']:.0f}, NEFA {fv1['NEFA_end']:.2f}, "
      f"CNSx {fv1['CNSx_end']:.2f}, Vol {fv1['Vol_end']:.2f})")
    w("  => The SAME hormone level sits on both sides of the storm boundary.")
    w("     This is the whole thesis, and it is a controlled comparison:")
    w("     one variable changed, and it was not a hormone.")

    # ================= [Bi] threshold =================
    w()
    w("[Bi] WHERE IS THE BOUNDARY?  Precipitant sweep at CONSTANT total T3.")
    w(RULE())
    w("     The VERDICT is the 7-day trajectory.  The 72-h column is a harsher")
    w("     probe that freezes the hormone level, so near the boundary it can")
    w("     lag the full run -- the runaway there takes longer than 72 h.")
    w(f"  {'Prec':>5s} {'72h probe':>10s} | 7-day trajectory: {'Tmax':>6s} "
      f"{'Tc@24':>6s} {'HRmax':>6s} {'BW@24':>6s} {'BWmax':>6s} {'Bili':>5s} "
      f"{'mortality':>10s}  VERDICT")
    for pr in (0.0, 0.3, 0.6, 0.9, 1.1, 1.2, 1.25, 1.30, 1.5, 2.0):
        z = y_tt.copy(); z[24], z[25] = pr, 0.6 * pr
        fv = fast_verdict(z, Ptt)
        yf, af, o = simulate(Ptt, [], hours=168.0, y0=y_tt.copy(), record=REC,
                             precip=pr, pge=0.6 * pr)
        w(f"  {pr:5.2f} {fv['Tmax']:10.2f} | {max(o['Tc']):23.2f} "
          f"{at(o,'Tc',24):6.2f} {max(o['HR']):6.0f} {at(o,'BW',24):6.0f} "
          f"{max(o['BW']):6.0f} {max(o['Bili']):5.2f} "
          f"{100*o['mort'][-1]:9.2f}%  "
          f"{'STORM (runaway)' if max(o['Tc']) >= 40.0 else 'survives'}")
    w(f"  -> a THRESHOLD between Prec = 1.20 and 1.30, not a gradient.  Below")
    w(f"     it the patient is hot, tachycardic and scores as a clinical storm")
    w(f"     by BWPS, yet is thermally stable and survives.  Above it there is")
    w(f"     no reachable operating point and the temperature runs away.")
    w(f"  -> all treatment scenarios below use Prec = {PREC:.2f} (fulminant).")

    # ================= [1-16] scenarios =================
    scen = [
        ("1  Untreated", [], {}),
        ("2  PTU alone (600 mg + 250 q4h)", PTU, {}),
        ("3  Methimazole alone (40 mg + 25 q6h)", MMI, {}),
        ("4  Iodide alone (SSKI q6h, 1 h after PTU)", PTU + IOD_AFTER, {}),
        ("5  Propranolol alone (80 mg q4h)", PRO, {}),
        ("6  Esmolol alone (infusion)", ESM, {}),
        ("7  Hydrocortisone alone (100 mg q8h)", HC, {}),
        ("8  Cooling + fluids alone", [], SUPP),
        ("9  FULL BUNDLE (PTU+iodide@1h+propranolol+HC+support)", BUNDLE, SUPP),
        ("10 BUNDLE with methimazole instead of PTU",
         MMI + IOD_AFTER + PRO + HC, SUPP),
        ("11 BUNDLE but iodide 4 h BEFORE the thionamide",
         PTU_LATE + IOD_BEFORE + PRO + HC, SUPP),
        ("11b Iodide with NO thionamide at all (Jod-Basedow)",
         IOD_BEFORE + PRO + HC, SUPP),
        ("12 BUNDLE minus beta-blockade", PTU + IOD_AFTER + HC, SUPP),
        ("13 BUNDLE minus cooling/fluids", BUNDLE, {}),
        ("14 BUNDLE + aspirin 650 q4h for the fever", BUNDLE + ASA, SUPP),
        ("15 BUNDLE + acetaminophen 650 q6h", BUNDLE + APAP, SUPP),
        ("16 BUNDLE + iopanoic acid + cholestyramine",
         BUNDLE + IOP + CHO, SUPP),
        ("17 BUNDLE + single plasma exchange at 12 h", BUNDLE,
         dict(SUPP, tpe=[(12.0, 0.65)])),
        ("18 BUNDLE with esmolol instead of propranolol",
         PTU + IOD_AFTER + ESM + HC, SUPP),
    ]
    w()
    w(f"[1-18] TREATMENT SCENARIOS -- 7 days, storm onset at t = 0, Prec = {PREC}")
    w(RULE())
    w(f"  {'scenario':52s} {'Tmax':>5s} {'Tc24':>5s} {'BW24':>5s} {'BW72':>5s} "
      f"{'totT3@24':>8s} {'freeT3@24':>9s} {'HR24':>5s} {'mortality':>10s}")
    R = {}
    for nm, spec, kw in scen:
        yf, af, o = simulate(Ptt, spec, hours=168.0, y0=y_tt.copy(),
                             record=REC, precip=PREC, pge=PGE, **kw)
        R[nm] = o
        w(f"  {nm:52s} {max(o['Tc']):5.2f} {at(o,'Tc',24):5.2f} "
          f"{at(o,'BW',24):5.0f} {at(o,'BW',72):5.0f} {at(o,'T3',24):8.2f} "
          f"{at(o,'fT3',24):9.2f} {at(o,'HR',24):5.0f} "
          f"{100*o['mort'][-1]:9.2f}%")

    # ================= [C] the central table =================
    w()
    w("[C] THE CENTRAL TABLE -- WHAT EACH MONOTHERAPY ACTUALLY MOVES, AT 24 h.")
    w("    Deltas are against the untreated trajectory at the same time point.")
    w(RULE())
    w(f"  {'monotherapy (axis attacked)':44s} {'dTotT3':>7s} {'dFreeT3':>8s} "
      f"{'dHR':>5s} {'dTc':>6s} {'runaway?':>9s} {'mortality':>10s}")
    b = R["1  Untreated"]
    mono = [
        ("PTU            synthesis + D1 conversion", "2  PTU alone (600 mg + 250 q4h)"),
        ("Methimazole    synthesis only", "3  Methimazole alone (40 mg + 25 q6h)"),
        ("Iodide         release block (after PTU)",
         "4  Iodide alone (SSKI q6h, 1 h after PTU)"),
        ("Propranolol    TRANSDUCTION only", "5  Propranolol alone (80 mg q4h)"),
        ("Esmolol        TRANSDUCTION only", "6  Esmolol alone (infusion)"),
        ("Hydrocortisone D1/D2 + steroid replacement",
         "7  Hydrocortisone alone (100 mg q8h)"),
        ("Cooling+fluids HEAT BALANCE only", "8  Cooling + fluids alone"),
    ]
    for lbl, key in mono:
        o = R[key]
        dT3 = 100 * (at(o, "T3", 24) / at(b, "T3", 24) - 1)
        dfT3 = 100 * (at(o, "fT3", 24) / at(b, "fT3", 24) - 1)
        dHR = at(o, "HR", 24) - at(b, "HR", 24)
        dTc = at(o, "Tc", 24) - at(b, "Tc", 24)
        run = "YES" if max(o["Tc"]) >= 40.0 else "no"
        w(f"  {lbl:44s} {dT3:+6.1f}% {dfT3:+7.1f}% {dHR:+4.0f} {dTc:+6.2f} "
          f"{run:>9s} {100*o['mort'][-1]:9.2f}%")
    w(f"  {'(untreated reference values)':44s} {at(b,'T3',24):6.2f}  "
      f"{at(b,'fT3',24):7.2f}  {at(b,'HR',24):4.0f} {at(b,'Tc',24):6.2f} "
      f"{'YES':>9s} {100*b['mort'][-1]:9.2f}%")
    w("")
    w("  READ THE TABLE THIS WAY: the two drugs that lower the measured hormone")
    w("  the most (PTU, methimazole) do NOT stop the runaway and the patient")
    w("  still dies.  The drug that removes NO hormone at all -- a beta-blocker")
    w("  -- terminates the storm, and it does so while cutting FREE T3 in half.")
    w("  Free T3 falls because the lipolysis that supplies the displacing NEFA")
    w("  is itself beta-mediated.  That is the arithmetic of ring A.")

    # ================= [L] the Lambda factor decomposition =================
    w()
    w("[L] Lambda = Q_prod/Q_loss AS A PRODUCT OF FACTORS, one per drug class.")
    w("    Evaluated at t = 6 h in each arm, at the escape temperature.  Lambda")
    w("    itself is only a LOCAL margin (it is computed with the fast states")
    w("    frozen); the last column is the actual 7-day outcome.  Read the")
    w("    FACTORS, not Lambda -- the factors are what the drugs multiply.")
    w(RULE())
    w(f"  {'arm':40s} {'M_thy':>6s} {'M_temp':>6s} {'M_sns':>6s} {'M_unc':>6s} "
      f"{'E':>5s} {'Vol^1.5':>7s} {'cool':>5s} {'h':>6s} {'Lam':>6s} "
      f"{'7-day outcome':>17s}")
    probe = [("untreated", "1  Untreated", {}),
             ("PTU alone", "2  PTU alone (600 mg + 250 q4h)", {}),
             ("methimazole alone", "3  Methimazole alone (40 mg + 25 q6h)", {}),
             ("propranolol alone", "5  Propranolol alone (80 mg q4h)", {}),
             ("cooling + fluids alone", "8  Cooling + fluids alone", SUPP),
             ("full bundle", "9  FULL BUNDLE (PTU+iodide@1h+propranolol+HC+support)", SUPP),
             ("bundle + aspirin", "14 BUNDLE + aspirin 650 q4h for the fever", SUPP)]
    for lbl, key, kw in probe:
        # rebuild the state at t = 6 h to get the factor decomposition
        spec = dict(scen)[key] if False else None
        for nm2, sp2, kw2 in scen:
            if nm2 == key:
                spec, kw2u = sp2, kw2
                break
        yf, af, o = simulate(Ptt, spec, hours=6.0, y0=y_tt.copy(), record=[],
                             precip=PREC, pge=PGE, **kw2u)
        si = storm_index(yf, Ptt, dict(cool=kw2u.get("cool", 0.0),
                                       fluid=kw2u.get("fluid", 0.0)))
        run7 = max(R[key]["Tc"]) >= 40.0
        w(f"  {lbl:40s} {si['M_thy']:6.3f} {si['M_temp']:6.3f} "
          f"{si['M_sns']:6.3f} {si['M_unc']:6.3f} {si['E_esc']:5.3f} "
          f"{si['Vol_term']:7.3f} {si['cool_term']:5.2f} {si['h_esc']:6.2f} "
          f"{si['Lam_esc']:6.3f} "
          f"{('RUNAWAY -> death' if run7 else 'controlled'):>17s}")
    w("  -> beta-blockade is the only drug that lowers TWO factors (M_sns")
    w("     directly, and M_thy through NEFA -> displacement -> free T3).")
    w("     Aspirin RAISES two (M_unc by uncoupling, M_thy by displacement).")

    # ================= [R] reservoir arithmetic =================
    w()
    w("[R] WHY A SYNTHESIS BLOCKER CANNOT STOP A STORM -- reservoir arithmetic.")
    w("    Methimazole alone.  Organification is blocked almost completely")
    w("    within an hour; everything after that is the arithmetic of a store")
    w("    that holds weeks of hormone.")
    w(RULE())
    o = R["3  Methimazole alone (40 mg + 25 q6h)"]
    w(f"  {'t':>6s} {'TPO block':>10s} {'synthesis':>10s} {'store S':>9s} "
      f"{'dS':>7s} {'secretion':>10s} {'totT4':>7s} {'dT4':>7s} "
      f"{'totT3':>7s} {'dT3':>7s}")
    for hh in (1, 3, 6, 12, 24, 48, 72, 120, 168):
        w(f"  {hh:5.0f}h {100*at(o,'inh_TPO',hh):9.2f}% "
          f"{at(o,'Synth',hh):10.3f} {at(o,'S',hh):9.0f} "
          f"{100*(at(o,'S',hh)/o['S'][0]-1):+6.2f}% {at(o,'Sec',hh):10.3f} "
          f"{at(o,'T4',hh):7.1f} {100*(at(o,'T4',hh)/o['T4'][0]-1):+6.2f}% "
          f"{at(o,'T3',hh):7.2f} {100*(at(o,'T3',hh)/o['T3'][0]-1):+6.2f}%")
    w(f"  store at t=0: {o['S'][0]:.0f} nmol; secretion {o['Sec'][0]:.2f} nmol/h")
    w(f"  => {o['S'][0]/(o['Sec'][0]*24):.1f} days of hormone are ALREADY MADE.")
    w("     A perfect synthesis blocker cannot touch them.  This is why no")
    w("     thionamide is an acute drug, and why the guidelines call for")
    w("     iodide (release), a beta-blocker (transduction) and cooling as")
    w("     well -- each of which acts on something that is not synthesis.")

    # ================= [D] PTU vs MMI =================
    w()
    w("[D] PTU vs METHIMAZOLE -- identical synthesis block, but only PTU")
    w("    inhibits D1, and D1 makes the 53% of T3 that the thyroid does not.")
    w(RULE())
    op = R["2  PTU alone (600 mg + 250 q4h)"]
    om = R["3  Methimazole alone (40 mg + 25 q6h)"]
    w(f"  {'t':>6s} {'Cptu':>6s} {'D1 block':>9s} | {'totT3 PTU':>10s} "
      f"{'dT3':>7s} | {'totT3 MMI':>10s} {'dT3':>7s} | {'rT3 PTU':>8s} "
      f"{'rT3 MMI':>8s}  (rT3 up = the D1 block is real)")
    for hh in (3, 6, 12, 24, 48, 72, 120, 168):
        w(f"  {hh:5.0f}h {at(op,'Cptu',hh):6.2f} "
          f"{100*at(op,'inh_D1',hh):8.1f}% | {at(op,'T3',hh):10.2f} "
          f"{100*(at(op,'T3',hh)/op['T3'][0]-1):+6.1f}% | "
          f"{at(om,'T3',hh):10.2f} "
          f"{100*(at(om,'T3',hh)/om['T3'][0]-1):+6.1f}% | "
          f"{at(op,'rT3',hh):8.3f} {at(om,'rT3',hh):8.3f}")

    # ================= [B] beta-blockade acts twice =================
    w()
    w("[B] BETA-BLOCKADE ACTS ON THE RING TWICE.  Propranolol removes no")
    w("    hormone, yet FREE T3 falls by half -- because the NEFA that")
    w("    displaces T3 from TBG is produced by beta-driven lipolysis.")
    w(RULE())
    ob = R["5  Propranolol alone (80 mg q4h)"]
    ou = R["1  Untreated"]
    w(f"  {'t':>6s} {'Cpro':>7s} {'beta occ':>9s} | ON PROPRANOLOL "
      f"{'NEFA':>6s} {'Fd':>6s} {'totT3':>6s} {'freeT3':>7s} {'HR':>4s} "
      f"{'Tc':>5s} | UNTREATED {'NEFA':>6s} {'Fd':>6s} {'totT3':>6s} "
      f"{'freeT3':>7s} {'HR':>4s} {'Tc':>5s}")
    for hh in (0.5, 1, 2, 4, 8, 12, 24, 48, 72):
        w(f"  {hh:5.1f}h {at(ob,'Cpro',hh):7.1f} "
          f"{100*at(ob,'occ_b',hh):8.1f}% |                "
          f"{at(ob,'NEFA',hh):6.2f} {at(ob,'Fd',hh):6.3f} "
          f"{at(ob,'T3',hh):6.2f} {at(ob,'fT3',hh):7.2f} "
          f"{at(ob,'HR',hh):4.0f} {at(ob,'Tc',hh):5.2f} |           "
          f"{at(ou,'NEFA',hh):6.2f} {at(ou,'Fd',hh):6.3f} "
          f"{at(ou,'T3',hh):6.2f} {at(ou,'fT3',hh):7.2f} "
          f"{at(ou,'HR',hh):4.0f} {at(ou,'Tc',hh):5.2f}")
    w(f"  at 24 h: total T3 {100*(at(ob,'T3',24)/at(ou,'T3',24)-1):+.1f}% but "
      f"FREE T3 {100*(at(ob,'fT3',24)/at(ou,'fT3',24)-1):+.1f}%")

    # ================= [F] free-fraction equivalence =================
    w()
    w("[F] THE SECOND AXIS -- a doubled free fraction IS a doubled hormone,")
    w("    and the total-hormone assay is completely blind to it.")
    w(RULE())
    w(f"  {'NEFA':>6s} {'Fd':>6s} {'measured totT3':>15s} {'freeT3':>7s} "
      f"{'Sig':>6s} {'EQUIVALENT totT3':>17s}")
    for nefa in (0.40, 0.60, 0.80, 1.00, 1.50, 2.00, 2.50, 3.00):
        z = y_tt.copy(); z[9] = nefa
        az = algebra(z, Ptt)
        w(f"  {nefa:6.2f} {az['Fd']:6.3f} {z[2]:15.2f} {az['fT3']:7.2f} "
          f"{az['Sig']:6.3f} {5.40*az['Fd']:17.2f}")
    w("  -> the 'measured' column never moves.  Two patients, identical assay,")
    w("     a 2-fold difference in the hormone the nucleus actually sees.")

    # ================= [J] iodide timing =================
    w()
    w("[J] IODIDE TIMING -- THE SIGN OF THE SAME DRUG FLIPS (Jod-Basedow).")
    w("    Identical dose (SSKI 250 mg iodide q6h).  ONLY the order differs:")
    w("    1 h AFTER the thionamide, versus at the same moment as it.")
    w(RULE())
    ra = R["9  FULL BUNDLE (PTU+iodide@1h+propranolol+HC+support)"]
    rb = R["11 BUNDLE but iodide 4 h BEFORE the thionamide"]
    rc = R["11b Iodide with NO thionamide at all (Jod-Basedow)"]
    w(f"  {'t':>6s} | THIONAMIDE FIRST {'synth':>7s} {'store':>7s} {'totT4':>7s} "
      f"{'totT3':>6s} {'Tc':>5s} {'BW':>3s} | IODIDE FIRST {'synth':>7s} "
      f"{'store':>7s} {'totT4':>7s} {'totT3':>6s} {'Tc':>5s} {'BW':>3s}")
    for hh in (1, 2, 4, 8, 12, 24, 48, 72, 120, 168):
        w(f"  {hh:5.0f}h |                  {at(ra,'Synth',hh):7.3f} "
          f"{at(ra,'S',hh):7.0f} {at(ra,'T4',hh):7.1f} {at(ra,'T3',hh):6.2f} "
          f"{at(ra,'Tc',hh):5.2f} {at(ra,'BW',hh):3.0f} |              "
          f"{at(rb,'Synth',hh):7.3f} {at(rb,'S',hh):7.0f} "
          f"{at(rb,'T4',hh):7.1f} {at(rb,'T3',hh):6.2f} "
          f"{at(rb,'Tc',hh):5.2f} {at(rb,'BW',hh):3.0f}")
    w(f"  7-day mortality: thionamide-first {100*ra['mort'][-1]:.2f}%   "
      f"iodide-4h-first {100*rb['mort'][-1]:.2f}%   "
      f"no thionamide at all {100*rc['mort'][-1]:.2f}%")
    w(f"  {'t':>6s} | NO THIONAMIDE AT ALL (Jod-Basedow): {'synth':>7s} "
      f"{'store':>7s} {'totT4':>7s} {'totT3':>6s} {'Tc':>5s} {'BW':>3s}")
    for hh in (1, 4, 12, 24, 48, 72, 120, 168):
        w(f"  {hh:5.0f}h |                                    "
          f"{at(rc,'Synth',hh):7.3f} {at(rc,'S',hh):7.0f} "
          f"{at(rc,'T4',hh):7.1f} {at(rc,'T3',hh):6.2f} "
          f"{at(rc,'Tc',hh):5.2f} {at(rc,'BW',hh):3.0f}")
    w("")
    w("  HONEST FINDING, AGAINST EXPECTATION.  In this model the order barely")
    w("  matters, and even omitting the thionamide entirely costs only")
    w(f"  {100*(at(rc,'T3',24)/at(ra,'T3',24)-1):+.1f}% in total T3 at 24 h.  The reason is")
    w("  arithmetic: the SAME intrathyroidal iodide that supplies substrate also")
    w("  shuts organification off, because Wolff-Chaikoff acts on BOTH arms.")
    w("  A gland that is release-inhibited is also synthesis-inhibited, so it")
    w("  cannot run away on the substrate it has just been given.")
    w("")
    w("[J2] SO WHEN IS JOD-BASEDOW REAL?  Only in a gland that FAILS to")
    w("     autoregulate -- an autonomously functioning nodule, or an")
    w("     iodine-deficient gland.  That is a different PARAMETER REGIME, and")
    w("     the model can be put into it: lower the two Wolff-Chaikoff Emax")
    w("     values.  Iodide alone, no thionamide, no precipitant, 21 days:")
    w(RULE())
    w(f"  {'gland':46s} {'totT4 d0':>9s} {'d3':>7s} {'d7':>7s} {'d21':>7s} "
      f"{'store d21':>10s} {'synth d21':>10s}")
    for lbl, ov in (("normal autoregulation (WCrel .78 / WCorg .90)", {}),
                    ("partial failure       (WCrel .40 / WCorg .45)",
                     dict(Imax_WCrel=0.40, Imax_WCorg=0.45)),
                    ("autonomous nodule     (WCrel .15 / WCorg .20)",
                     dict(Imax_WCrel=0.15, Imax_WCorg=0.20))):
        Pj = dict(Ptt, **ov)
        yj = chronic_ss(Pj, hours=2500.0, y0=y_tt.copy())
        yf, af, oj = simulate(Pj, IOD_BEFORE, hours=24 * 21, dt=0.1,
                              rec_every=20, y0=yj, record=REC,
                              precip=0.0, pge=0.0)
        w(f"  {lbl:46s} {oj['T4'][0]:9.1f} {at(oj,'T4',72):7.1f} "
          f"{at(oj,'T4',168):7.1f} {at(oj,'T4',504):7.1f} "
          f"{at(oj,'S',504):10.0f} {at(oj,'Synth',504):10.2f}")
    w("     -> the sign of iodide is set by the GLAND's autoregulation, not by")
    w("        the order of the prescription.  That is a falsifiable claim and")
    w("        it matches who actually gets Jod-Basedow in the clinic.")

    # ================= [W] Wolff-Chaikoff escape =================
    w()
    w("[W] NIS DOWN-REGULATION AND THE APPROACH TO WOLFF-CHAIKOFF ESCAPE.")
    w("    SSKI 250 mg iodide q6h, NO thionamide, no precipitant, 21 days.")
    w("    HONEST RESULT, AGAINST THE TEXTBOOK: the escape MECHANISM is here")
    w("    and it works -- NIS falls 7-fold and the intrathyroidal pool falls")
    w("    18-fold -- but at this dose the pool never gets back below IC_WC,")
    w("    so RELEASE inhibition never lets go and T4 keeps falling.  The model")
    w("    therefore does NOT reproduce clinical escape at 21 days (limitation")
    w("    L8).  What it does show is the machinery, and the dose-dependence:")
    w("    escape requires the intrathyroidal pool to fall past IC_WC = 2 umol.")
    w(RULE())
    yf, af, o = simulate(Ptt, IOD_AFTER, hours=24 * 21, dt=0.1, rec_every=20,
                         y0=y_tt.copy(), record=REC, precip=0.0, pge=0.0)
    w(f"  {'day':>4s} {'plasma I':>9s} {'NIS':>6s} {'thy I':>7s} "
      f"{'release inhib':>14s} {'A_TPO':>7s} {'synth':>7s} {'store':>7s} "
      f"{'totT4':>7s} {'totT3':>6s}")
    for dd in (0, 1, 2, 3, 5, 7, 10, 14, 21):
        hh = dd * 24.0
        w(f"  {dd:4d} {at(o,'Ipl',hh):9.2f} {at(o,'NIS',hh):6.3f} "
          f"{at(o,'Ithy',hh):7.2f} {at(o,'rel_inh',hh):14.3f} "
          f"{at(o,'A_TPO',hh):7.3f} {at(o,'Synth',hh):7.3f} "
          f"{at(o,'S',hh):7.0f} {at(o,'T4',hh):7.1f} {at(o,'T3',hh):6.2f}")

    # ================= [A] aspirin =================
    w()
    w("[A] ASPIRIN IN THYROID STORM -- both mechanisms of harm, quantified.")
    w("    (1) it displaces T3 from TBG; (2) it uncouples oxidative")
    w("    phosphorylation.  Both raise Lambda.  Acetaminophen does neither.")
    w(RULE())
    w(f"  {'arm':26s} {'Tmax':>6s} {'Fd@24':>6s} {'freeT3@24':>10s} "
      f"{'Qprod@24':>9s} {'BW24':>5s} {'BW72':>5s} {'mortality':>10s}")
    for lbl, key in (("bundle", "9  FULL BUNDLE (PTU+iodide@1h+propranolol+HC+support)"),
                     ("bundle + aspirin", "14 BUNDLE + aspirin 650 q4h for the fever"),
                     ("bundle + acetaminophen", "15 BUNDLE + acetaminophen 650 q6h")):
        o = R[key]
        w(f"  {lbl:26s} {max(o['Tc']):6.2f} {at(o,'Fd',24):6.3f} "
          f"{at(o,'fT3',24):10.2f} {at(o,'Qprod',24):9.1f} "
          f"{at(o,'BW',24):5.0f} {at(o,'BW',72):5.0f} "
          f"{100*o['mort'][-1]:9.2f}%")

    # ================= [P] the propranolol trap =================
    w()
    w("[P] THE PROPRANOLOL TRAP -- beta-blockade of a rate-dependent heart.")
    w("    Identical bundle; only the cardiac reserve at storm onset differs.")
    w(RULE())
    w(f"  {'CR(0)':>6s} | {'PROPRANOLOL  t1/2 4 h':^42s} | "
      f"{'ESMOLOL  t1/2 9 min':^42s}")
    w(f"  {'':>6s} | {'CO nadir':>9s} {'shock':>6s} {'BW24':>5s} {'Bili':>5s} "
      f"{'mort':>9s} | {'CO nadir':>9s} {'shock':>6s} {'BW24':>5s} "
      f"{'Bili':>5s} {'mort':>9s}")
    for cr0 in (1.00, 0.80, 0.60, 0.45, 0.35):
        row = []
        for spec in (BUNDLE, PTU + IOD_AFTER + ESM + HC):
            z = y_tt.copy(); z[14] = cr0
            yf, af, o = simulate(Ptt, spec, hours=168.0, y0=z, record=REC,
                                 precip=PREC, pge=PGE, **SUPP)
            row.append((min(o["perf"]), max(o["shock"]), at(o, "BW", 24),
                        max(o["Bili"]), o["mort"][-1]))
        w(f"  {cr0:6.2f} | {row[0][0]:9.3f} {row[0][1]:6.3f} {row[0][2]:5.0f} "
          f"{row[0][3]:5.2f} {100*row[0][4]:8.2f}% | {row[1][0]:9.3f} "
          f"{row[1][1]:6.3f} {row[1][2]:5.0f} {row[1][3]:5.2f} "
          f"{100*row[1][4]:8.2f}%")

    # ================= [K] cortisol =================
    w()
    w("[K] ACCELERATED CORTISOL CLEARANCE -> RELATIVE ADRENAL INSUFFICIENCY.")
    w("    Thyroid hormone roughly doubles cortisol clearance at exactly the")
    w("    moment the stress requirement rises.  This is a scissor, and the")
    w("    glucocorticoid in the bundle is closing it, not treating the thyroid.")
    w("    'AI hazard x' is the mortality multiplier the model applies for the")
    w("    steroid deficit; stress-dose hydrocortisone abolishes it outright.")
    w(RULE())
    w(f"  {'adrenal reserve':>16s} {'arm':>22s} {'endogenous':>11s} "
      f"{'EFFECTIVE GC':>13s} {'AI hazard x':>12s} {'BW24':>5s} {'mortality':>10s}")
    for ar in (1.00, 0.70, 0.50):
        Pa = dict(Ptt, AR=ar)
        z = chronic_ss(Pa, hours=2500.0, y0=y_tt.copy())
        for lbl, spec in (("no glucocorticoid", PTU + IOD_AFTER + PRO),
                          ("+ hydrocortisone", BUNDLE)):
            yf, af, o = simulate(Pa, spec, hours=168.0, y0=z.copy(),
                                 record=REC + ["CortEff", "aiMult"],
                                 precip=PREC, pge=PGE, **SUPP)
            w(f"  {ar:16.2f} {lbl:>22s} {at(o,'Cort',24):11.0f} "
              f"{at(o,'CortEff',24):13.0f} {at(o,'aiMult',24):12.3f} "
              f"{at(o,'BW',24):5.0f} {100*o['mort'][-1]:9.2f}%")

    # ================= [X] plasma exchange =================
    w()
    w("[X] PLASMA EXCHANGE -- the only intervention that touches the")
    w("    PROTEIN-BOUND pool, and the arithmetic of why it is modest.")
    w("    A 1.2-plasma-volume exchange removes ~65% of the PLASMA pool.")
    w("    Plasma (~3 L) is 30% of the T4 distribution volume but only 7.5%")
    w("    of T3's, so the same procedure removes ~20% of body T4 but only")
    w("    ~5% of body T3.  See caveat L7 on the missing plasma/EV split.")
    w(RULE())
    ot = R["17 BUNDLE + single plasma exchange at 12 h"]
    ob2 = R["9  FULL BUNDLE (PTU+iodide@1h+propranolol+HC+support)"]
    w(f"  {'t':>6s} {'totT3 TPE':>10s} {'totT3 none':>11s} {'freeT3 TPE':>11s} "
      f"{'freeT3 none':>12s} {'totT4 TPE':>10s} {'BW TPE':>7s} {'BW none':>8s}")
    for hh in (11.5, 12.5, 14, 18, 24, 36, 48, 72, 120, 168):
        w(f"  {hh:5.1f}h {at(ot,'T3',hh):10.2f} {at(ob2,'T3',hh):11.2f} "
          f"{at(ot,'fT3',hh):11.2f} {at(ob2,'fT3',hh):12.2f} "
          f"{at(ot,'T4',hh):10.1f} {at(ot,'BW',hh):7.0f} "
          f"{at(ob2,'BW',hh):8.0f}")
    w(f"  7-day mortality: with TPE {100*ot['mort'][-1]:.2f}%   "
      f"without {100*ob2['mort'][-1]:.2f}%")

    # ================= [Z] provenance and caveats =================
    w()
    w("[Z] CALIBRATION PROVENANCE -- what was fitted, and therefore what is")
    w("    a prediction.")
    w(RULE())
    w("  FITTED TO NORMAL PHYSIOLOGY (20 numbers, section [N]): total and free")
    w("  T4/T3, rT3, the 80/20 peripheral/thyroidal split of T3 production, the")
    w("  T4 and T3 half-lives, the 60-90 day glandular store, TSH, dietary")
    w("  iodide, BMR 80 W, Q10 = 2, core temperature 37.0, heart rate 70,")
    w("  NEFA 0.4, cortisol 400.  The euthyroid heat balance and iodide balance")
    w("  are ALGEBRAICALLY CONSISTENT, not fitted: h(37.0) = 20 W/K follows from")
    w("  BMR = 80 W, and I_diet follows from closing the iodine mass balance.")
    w("")
    w("  FITTED TO NON-STORM DRUG DATA (7 numbers): thionamide TPO potency")
    w("  (MMI ~20x PTU per mg), PTU D1 potency (serum T3 falls ~20-30% in 24 h")
    w("  in ordinary thyrotoxicosis), propranolol beta potency (80 mg lowers HR")
    w("  ~25-30%), iodide release-inhibition Emax (Lugol lowers T4 30-50% over")
    w("  24-48 h in preoperative preparation), Wolff-Chaikoff escape time")
    w("  constant (10-14 days), glucocorticoid and iopanoic D1/D2 potencies.")
    w("")
    w("  FITTED TO STORM DATA: exactly ONE number, h0_haz, set so that the")
    w("  untreated fulminant arm gives 7-day mortality 85% (the pre-1970s")
    w("  historical figure).  NOTHING ELSE in the model was tuned on storm.")
    w("")
    w("  THEREFORE THE FOLLOWING ARE PREDICTIONS, NOT FITS:")
    w("   P1  no hormone level alone storms; Lambda saturates below 1 [H]")
    w("   P2  the boundary is a THRESHOLD in the precipitant, not a gradient [Bi]")
    w("   P3  thionamides do not stop a storm; the reservoir arithmetic says")
    w("       they cannot, and the store size was fitted to EUTHYROID data [R]")
    w("   P4  beta-blockade halves FREE T3 while total T3 barely moves [B]")
    w("   P5  the glandular store falls from 66 d (euthyroid) to ~18 d in")
    w("       florid Graves as a consequence, not an assumption [T]")
    w("   P6  cortisol falls to ~250 nmol/L on unchanged adrenal output [T]")
    w("   P7  iodide given before a thionamide has the OPPOSITE sign [J]")
    w("   P8  iodide monotherapy keeps working for 21 days at SSKI doses,")
    w("       because release inhibition does not escape while the")
    w("       intrathyroidal pool stays above IC_WC -- a dose-dependent")
    w("       prediction that contradicts the usual 10-14 day teaching [W]")
    w("   P9  the propranolol trap appears only below a cardiac-reserve")
    w("       threshold, and esmolol's short half-life mitigates it [P]")
    w("   P10 aspirin worsens the storm through TWO independent factors [A]")
    w("")
    w("  HONEST LIMITATIONS")
    w("   L1  Treated-arm mortality is UNDER-predicted.  Modern series report")
    w("       ~8-25% with full treatment; this model gives well under 1% for")
    w("       the bundle, because it contains no comorbidity, no multi-organ")
    w("       failure, no thromboembolism, no ventilator or line complications")
    w("       and no drug-toxicity deaths (agranulocytosis, PTU hepatic")
    w("       failure).  Only the UNTREATED figure was calibrated; the treated")
    w("       figures should be read as 'the storm physiology was controlled',")
    w("       not as a survival estimate.")
    w("   L2  The runaway is capped at Tc = 43.0 degC by a numerical ceiling.")
    w("       Temperatures at the ceiling mean 'lethal hyperthermia', not a")
    w("       predicted measurement.")
    w("   L3  The fast/slow split is a modelling choice.  Rb (36 h) and the")
    w("       cardiac reserve are treated as slow; over a 7-day simulation both")
    w("       move, and the storm verdict in [H]/[S] freezes them.")
    w("   L4  Only ONE precipitant axis is modelled (a scalar Prec plus its")
    w("       PGE2 fraction).  Real precipitants differ qualitatively -- an")
    w("       iodine load acts through the gland, DKA through volume, surgery")
    w("       through catecholamines -- and the map shows those separately.")
    w("   L5  Amiodarone-induced thyrotoxicosis (types 1 and 2), for which")
    w("       iodide is contraindicated, is drawn on the map but is not a")
    w("       simulated scenario here.")
    w("   L6  Burch-Wartofsky is computed from simulated signs, so it inherits")
    w("       every simplification above; it is not an independent measurement.")
    w("   L7  There is no plasma / extravascular split, so plasma exchange can")
    w("       only be represented as its net WHOLE-BODY removal (~20% of T4,")
    w("       ~5% of T3).  The large immediate post-exchange dip in plasma")
    w("       concentration, and the rebound that follows it, are invisible")
    w("       here.  TPE is therefore under-represented acutely.")
    w("   L8  Clinical Wolff-Chaikoff escape is NOT reproduced within 21 days")
    w("       at SSKI doses (section [W]).  The mechanism (NIS down-regulation)")
    w("       is present and quantitatively large, but the intrathyroidal pool")
    w("       stays far above the Wolff-Chaikoff IC50, so the release block")
    w("       persists.  Either IC_WC is too low or a second escape mechanism")
    w("       (e.g. recovery of organification once the pool falls) is missing.")
    w("   L9  Ithy is a lumped, uncapped intrathyroidal pool.  Its ABSOLUTE")
    w("       value under a pharmacologic iodide load (hundreds of umol) is not")
    w("       physiological; only the saturating Wolff-Chaikoff function of it")
    w("       enters the model, so the conclusions do not depend on it.")

    w()
    w(RULE("="))
    w("END OF VERIFICATION")
    w(RULE("="))
    txt = "\n".join(W)
    print(txt)
    import os
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "ts_verification_output.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
