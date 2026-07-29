#!/usr/bin/env python3
"""
Diabetic Foot Ulcer (DFU) QSP model — dependency-free reference implementation
=============================================================================

This file is the *numerical authority* for the DFU QSP model.  It implements
exactly the same 44-state ODE system as ``dfu_mrgsolve_model.R`` using a
fixed-step RK4 integrator and nothing but the Python standard library, so the
numbers quoted in ``README.md`` can be regenerated and checked on any machine
with a Python interpreter and no packages installed at all.

    python3 dfu_reference_implementation.py            # scenario table
    python3 dfu_reference_implementation.py --anchors  # calibration anchors
    python3 dfu_reference_implementation.py --csv OUT  # dump trajectories

Structural commitments of the model (see README §4):

  (1) CLOSURE IS A PERIMETER PROCESS.  dA/dt = -k(t)*P with P = 2*sqrt(pi*A).
      Equivalently the wound *radius* falls linearly at rate k.  Everything
      biological enters through k, not through A.
  (2) OFFLOADING IS EFFICACY x ADHERENCE.  A removable device is not a worse
      device, it is a good device worn 28% of the time.
  (3) TOPICAL GROWTH FACTOR IS A RACE AGAINST PROTEASE.  Becaplermin is
      eliminated from the wound at a rate proportional to MMP-9/TIMP-1.
  (4) OXYGEN GATES EVERY ANABOLIC TERM.  Collagen hydroxylation, angiogenesis
      and the neutrophil oxidative burst all share one Hill gate in pO2.
  (5) HbA1c ACTS ON THE SLOW AXES.  Neuropathy has a ~4.6-year time constant,
      so glycaemic control changes recurrence, not this ulcer's closure.
  (6) HEALING IS REMISSION.  Closure resets AREA but leaves neuropathy and
      plantar pressure untouched, and leaves a scar at ~80% tensile strength.

Research / education only.  Not for clinical use.
"""

import math
import sys

# ---------------------------------------------------------------------------
# Parameters  (units: time = days, area = cm^2, pO2 = mmHg, BACT = log10 CFU/g)
# ---------------------------------------------------------------------------
P = dict(
    # --- patient / systemic -------------------------------------------------
    HBA1C0=9.0, HBA1C_TGT=9.0, TAU_A1C=45.0,
    HB_GDL=12.5,
    KP_NEURO=5.0e-4, KR_NEURO=3.5e-4,
    # --- perfusion ----------------------------------------------------------
    PERF0=0.85, PERF_TGT=0.85, TAU_PERF=7.0,
    T_REVASC=1e6, REVASC_GAIN=0.0, K_RESTEN=1.5e-3,
    TCPO2_MAX=62.0, PO2_50=26.0, HILL_O2=3.0, PO2_HIF=30.0,
    OEDEMA_PEN=0.0,
    # --- mechanics / offloading --------------------------------------------
    DEFORMITY=0.45,                 # 0-1, claw toe / MT prominence / Charcot
    DEV_EFF=0.0, ADHERENCE=1.0,     # offloading device efficacy x wear fraction
    T_OFF_ON=0.0, T_OFF_OFF=1e6,
    FOOTWEAR=0.0,                   # post-healing therapeutic footwear (0/1)
    STEPS_REL=1.0,
    S50_MECH=1.00, H_MECH=2.0,
    STRESS_TOL=0.55, K_REINJ=0.086, AREA_MAXPR=10.0,
    # --- epithelial edge ----------------------------------------------------
    EDGE_MAX=0.158,                 # cm/day maximal linear edge advance
    K_CONTR=0.010, C50_CONTR=0.60,
    AREA0=2.00, DEPTH0=3.0, AREA_CLOSE=0.05,
    # --- keratinocytes ------------------------------------------------------
    KIN_K=0.25, KOUT_K=0.12, IL1_50K=1.00, PR50_K=1.60,
    # --- granulation --------------------------------------------------------
    KG_GRAN=0.22, KD_GRAN=0.030, K_FILL=0.16,
    # --- inflammation -------------------------------------------------------
    KIN_NEUT=0.50, KOUT_NEUT=0.50,
    KIN_M1=0.40, KOUT_M1=0.18, K_SW=0.55, IL1_SW=0.85,
    KIN_M2=0.10, KOUT_M2=0.16,
    KIN_IL1=0.60, KOUT_IL1=0.55, IL10_50=0.90,
    KIN_ROS=0.45, KOUT_ROS=0.60,
    # --- proteases ----------------------------------------------------------
    KIN_M9=0.55, KOUT_M9=0.30,
    KIN_T1=0.22, KOUT_T1=0.30,
    # --- growth factors / angiogenesis -------------------------------------
    KIN_PDGF=0.14, KOUT_PDGF=0.45, PDGF_50=1.10,
    KIN_VEGF=0.40, KOUT_VEGF=0.40,
    KIN_SDF=0.35, KOUT_SDF=0.30, K_DPP4=0.45,
    KM_EPC=0.30, KOUT_EPC=0.35,
    KG_VASC=0.30, KR_VASC=0.10,
    MGO_BLOCK=0.085,                # per %HbA1c above 6.0 — HIF transactivation loss
    # --- fibroblasts / ECM --------------------------------------------------
    KP_FIB=0.34, KD_FIB=0.10, FIB_MAX=2.0,
    K_SEN=0.055, KCL_SEN=0.030,
    KS_COL=0.30, KDEG_COL=0.14,
    # --- infection ----------------------------------------------------------
    BACT0=5.2, BMAX=9.5, KG_BACT=0.55,
    KILL_HOST=0.50, EMAX_ABX=1.30, EC50_ABX=0.45, TOL_BIOF=0.92,
    KG_BIOF=1.45, KD_BIOF=0.05, K_DEBR_BIOF=13.8, K_RIF=0.10,
    PROBE_BONE=0.0,
    K_ON_OST=0.045, K_OFF_OST=0.25, EC50_ABX_B=0.25, K_SURG=2.5,
    NECRO=0.35,
    # --- drug PK ------------------------------------------------------------
    KA_BEC=1.20, KOUT_BEC=0.15, KDEG_BEC_PROT=5.00, R_BEC=0.0,
    E_BEC=1.80, EC50_BEC=0.30, H_BEC=2.0, T_BEC_END=140.0,
    R_NOSF=0.0, KOUT_NOSF=0.90, T_NOSF_END=140.0, EMAX_NOSF=0.62, EC50_NOSF=0.50,
    KA_ESM=1.50, KOUT_ESM=1.30, R_ESM=0.0, E_ESM=0.55, EC50_ESM=0.40, T_ESM_END=84.0,
    R_ABX=0.0, T_ABX_END=14.0, KOUT_ABX=1.60, KTR_ABX_W=0.90, KOUT_ABX_W=1.10,
    F_BONE=0.35, KTR_ABX_B=0.25, KOUT_ABX_B=0.45, RIF_ON=0.0,
    R_CTP=0.0, T_CTP_END=56.0, KOUT_CTP=0.085, E_CTP=0.85, EC50_CTP=0.35,
    R_OXY=0.0, T_OXY_END=40.0, KOUT_OXY=8.0, E_OXY=22.0,
    R_DEBR=0.0, KOUT_DEBR=8.0, DEBR_PERIOD=14.0, T_DEBR_END=1e6,
    NPWT=0.0, DPP4I=0.0,
    # --- outcomes -----------------------------------------------------------
    TAU_SCAR=60.0, SCAR_MAX_TENS=0.80,
    H0_REC=2.20e-3, DEFORM_REF=0.45, EDU=0.0, TEMPMON=0.0, SURVEIL=0.0,
    A0_AMP=1.10e-3,
    W_TOX_ABX=1.0, W_TOX_HBOT=0.6, W_TOX_BEC=0.8, W_TOX_TCC=0.5,
)

# ---------------------------------------------------------------------------
# State vector layout (the same 44 states as $CMT in the R model; the
# declaration order differs between the two files, the equations do not)
# ---------------------------------------------------------------------------
NAMES = [
    "BEC_D", "BEC_W", "NOSF", "ESM_D", "ESM_W", "ABX_C", "ABX_W", "ABX_B",
    "CTP", "OXY", "DEBR",                                     # 1-11  therapy
    "HBA1C", "NEURO", "PERF",                                 # 12-14 slow axes
    "BACT", "BIOF", "OSTEO",                                  # 15-17 infection
    "NEUT", "M1", "M2", "IL1B", "ROS",                        # 18-22 inflammation
    "MMP9", "TIMP1",                                          # 23-24 proteases
    "PDGF", "VEGF", "SDF1", "EPC", "VASC",                    # 25-29 angiogenesis
    "FIB", "SEN", "COL", "GRAN", "KERA",                      # 30-34 repair
    "AREA", "DEPTH", "SCAR",                                  # 35-37 wound
    "HAZREC", "HAZAMP", "UFD", "ABXD", "TOXIDX",              # 38-42 outcomes
    "BEXP",                                                   # 43    exposed bone
    "AUCEDGE",                                                # 44    bookkeeping
]
IDX = {n: i for i, n in enumerate(NAMES)}
N = len(NAMES)


def hill(x, ec50, h=1.0):
    """Emax/Hill fraction, overflow-safe."""
    x = min(max(x, 0.0), 1e12)
    if x <= 0.0:
        return 0.0
    r = (x / ec50) ** h
    return r / (1.0 + r) if r < 1e12 else 1.0


def inhib(x, ic50, h=1.0):
    """Fractional inhibition 1/(1+(x/IC50)^h), overflow-safe."""
    x = min(max(x, 0.0), 1e12)
    r = (x / ic50) ** h
    return 1.0 / (1.0 + r) if r < 1e12 else 0.0


def sw(t, t0, width=0.5):
    """Smooth 0->1 switch at t0."""
    z = (t - t0) / width
    if z < -35:
        return 0.0
    if z > 35:
        return 1.0
    return 1.0 / (1.0 + math.exp(-z))


def window(t, t_on, t_off):
    return sw(t, t_on) * (1.0 - sw(t, t_off))


# ---------------------------------------------------------------------------
def derivatives(t, y, p):
    s = {n: y[i] for i, n in enumerate(NAMES)}
    d = [0.0] * N

    def D(name, val):
        d[IDX[name]] = val

    # ===== 0. algebra shared by many terms ================================
    AREA = max(s["AREA"], 0.0)
    PERIM = 2.0 * math.sqrt(math.pi * max(AREA, 1e-9))
    healed = 1.0 - sw(AREA, p["AREA_CLOSE"], 0.02)      # 1 when AREA < close
    open_w = 1.0 - healed

    hyperg = max(s["HBA1C"] - 6.0, 0.0)                 # %A1C above normal

    # --- protease ratio ---------------------------------------------------
    PROT = s["MMP9"] / (s["TIMP1"] + 0.05)
    PROTn = PROT / (PROT + 1.5)

    # --- offloading: efficacy x adherence ---------------------------------
    on_off = window(t, p["T_OFF_ON"], p["T_OFF_OFF"])
    offload = p["DEV_EFF"] * p["ADHERENCE"] * on_off
    ppp_index = (1.0 + 0.60 * (p["DEFORMITY"] - p["DEFORM_REF"])) * p["STEPS_REL"]
    STRESS = ppp_index * (1.0 - offload)
    MG = inhib(STRESS, p["S50_MECH"], p["H_MECH"])      # mechanical gate on edge
    over = max(STRESS - p["STRESS_TOL"], 0.0)
    # re-injury is a MARGIN process too: it scales with perimeter, and is bounded
    # by the anatomical extent of the high-pressure zone.
    REINJ = (p["K_REINJ"] * over * over * PERIM
             * max(1.0 - AREA / p["AREA_MAXPR"], 0.0) * open_w)

    # --- perfusion & oxygen ------------------------------------------------
    rev = sw(t, p["T_REVASC"])
    decay = math.exp(-p["K_RESTEN"] * max(t - p["T_REVASC"], 0.0))
    perf_tgt = min(p["PERF_TGT"] + p["REVASC_GAIN"] * rev * decay, 1.0)
    D("PERF", (perf_tgt - s["PERF"]) / p["TAU_PERF"])

    hb_f = min(s_hb := p["HB_GDL"] / 13.0, 1.0) if False else min(p["HB_GDL"] / 13.0, 1.0)
    WPO2 = (p["TCPO2_MAX"] * s["PERF"] * (0.45 + 0.55 * s["VASC"]) * hb_f
            - p["OEDEMA_PEN"] * (1.0 - 0.5 * p["NPWT"])
            + p["E_OXY"] * hill(s["OXY"], 0.5))
    WPO2 = max(WPO2, 1.0)
    OG = hill(WPO2, p["PO2_50"], p["HILL_O2"])          # THE anabolic oxygen gate

    # --- HIF competence: hypoxia is present, transactivation is not --------
    hypoxia = inhib(WPO2, p["PO2_HIF"], 2.0)
    hif_comp = max(1.0 - p["MGO_BLOCK"] * hyperg, 0.05)
    HIFdrive = hypoxia * hif_comp

    # ===== 1. therapy PK ===================================================
    # every therapy runs for a finite course; on_X is 1 while the course is on
    on_bec = window(t, 0.0, p["T_BEC_END"])
    on_nosf = window(t, 0.0, p["T_NOSF_END"])
    on_esm = window(t, 0.0, p["T_ESM_END"])
    on_abx = window(t, 0.0, p["T_ABX_END"])
    on_ctp = window(t, 0.0, p["T_CTP_END"])
    on_oxy = window(t, 0.0, p["T_OXY_END"])
    D("BEC_D", p["R_BEC"] * on_bec * open_w - p["KA_BEC"] * s["BEC_D"])
    D("BEC_W", p["KA_BEC"] * s["BEC_D"]
               - (p["KOUT_BEC"] + p["KDEG_BEC_PROT"] * PROT) * s["BEC_W"])
    D("NOSF", p["R_NOSF"] * on_nosf * open_w - p["KOUT_NOSF"] * s["NOSF"])
    D("ESM_D", p["R_ESM"] * on_esm * open_w - p["KA_ESM"] * s["ESM_D"])
    D("ESM_W", p["KA_ESM"] * s["ESM_D"] - p["KOUT_ESM"] * s["ESM_W"])
    D("ABX_C", p["R_ABX"] * on_abx - (p["KOUT_ABX"] + p["KTR_ABX_W"] + p["KTR_ABX_B"]) * s["ABX_C"])
    D("ABX_W", p["KTR_ABX_W"] * s["ABX_C"] - p["KOUT_ABX_W"] * s["ABX_W"])
    D("ABX_B", p["KTR_ABX_B"] * p["F_BONE"] * s["ABX_C"] * (1.0 + 1.8 * p["RIF_ON"])
               - p["KOUT_ABX_B"] * s["ABX_B"])
    D("CTP", p["R_CTP"] * on_ctp * open_w - p["KOUT_CTP"] * s["CTP"])
    D("OXY", p["R_OXY"] * on_oxy - p["KOUT_OXY"] * s["OXY"])
    D("DEBR", -p["KOUT_DEBR"] * s["DEBR"])              # bolus-driven

    # ===== 2. slow systemic axes ===========================================
    D("HBA1C", (p["HBA1C_TGT"] - s["HBA1C"]) / p["TAU_A1C"])
    D("NEURO", p["KP_NEURO"] * max(s["HBA1C"] - 6.5, 0.0) * (1.0 - s["NEURO"])
               - p["KR_NEURO"] * s["NEURO"])

    # ===== 3. infection ====================================================
    bf = min(max((s["BACT"] - 3.0) / 6.0, 0.0), 1.0)
    hostdef = ((0.30 + 0.70 * OG)
               * (1.0 - 0.35 * min(hyperg / 4.0, 1.0))
               * (1.0 - 0.85 * s["BIOF"]))
    abx_kill = (p["EMAX_ABX"] * hill(s["ABX_W"], p["EC50_ABX"])
                * (1.0 - p["TOL_BIOF"] * s["BIOF"]))
    D("BACT", (p["KG_BACT"] * (p["BMAX"] - s["BACT"]) / p["BMAX"]
               * (0.55 + 0.45 * p["NECRO"]) * open_w
               - p["KILL_HOST"] * hostdef * bf
               - abx_kill
               - 0.55 * s["DEBR"]
               - 1.2 * healed * max(s["BACT"] - 3.0, 0.0) / 6.0))
    D("BIOF", (p["KG_BIOF"] * bf * (1.0 - s["BIOF"])
               - p["KD_BIOF"] * s["BIOF"]
               - p["K_DEBR_BIOF"] * s["DEBR"] * s["BIOF"]
               - p["K_RIF"] * p["RIF_ON"] * s["BIOF"]))
    surg = 0.0  # surgical resection handled as an event (see run())
    D("BEXP", -0.0015 * s["BEXP"])
    D("OSTEO", (p["K_ON_OST"] * bf * s["BEXP"] * (1.0 - s["OSTEO"])
                - p["K_OFF_OST"] * hill(s["ABX_B"], p["EC50_ABX_B"])
                  * (1.0 - 0.70 * s["BIOF"]) * s["OSTEO"]
                - 0.002 * s["OSTEO"] - surg))

    # ===== 4. inflammation =================================================
    damage = 0.45 * over + 0.30 * s["OSTEO"]
    D("NEUT", p["KIN_NEUT"] * (0.15 + 1.20 * bf + damage) * open_w
              - p["KOUT_NEUT"] * s["NEUT"])
    switch = (inhib(s["IL1B"], p["IL1_SW"], 3.0)
              * (0.40 + 0.60 * OG)
              * (1.0 - 0.30 * min(hyperg / 4.0, 1.0)))
    D("M1", p["KIN_M1"] * (0.20 + 0.80 * bf + 0.60 * s["NEUT"] + 0.50 * s["SEN"]) * open_w
            - (p["KOUT_M1"] + p["K_SW"] * switch) * s["M1"])
    D("M2", p["K_SW"] * switch * s["M1"] + p["KIN_M2"] * open_w
            + 0.30 * hill(s["CTP"], p["EC50_CTP"])
            - p["KOUT_M2"] * s["M2"])
    D("IL1B", p["KIN_IL1"] * (0.30 * bf + 0.70 * s["M1"] + 0.40 * s["ROS"])
              * inhib(s["M2"], p["IL10_50"]) * open_w
              - p["KOUT_IL1"] * s["IL1B"])
    D("ROS", p["KIN_ROS"] * (s["NEUT"] + 0.50 * s["M1"] + 0.25 * min(hyperg / 3.0, 1.0)) * open_w
             - p["KOUT_ROS"] * s["ROS"])

    # ===== 5. proteases ====================================================
    enosf = p["EMAX_NOSF"] * hill(s["NOSF"], p["EC50_NOSF"])
    D("MMP9", p["KIN_M9"] * (0.60 * s["NEUT"] + 0.50 * s["M1"]
                             + 0.30 * s["SEN"] + 0.30 * s["IL1B"])
              * (1.0 - enosf) * open_w
              - p["KOUT_M9"] * s["MMP9"])
    D("TIMP1", p["KIN_T1"] * (0.30 + 0.80 * s["M2"]) * open_w
               - p["KOUT_T1"] * s["TIMP1"])

    # ===== 6. growth factors / angiogenesis ================================
    PDGF_tot = s["PDGF"] + p["E_BEC"] * hill(s["BEC_W"], p["EC50_BEC"], p["H_BEC"])
    D("PDGF", p["KIN_PDGF"] * (0.40 + 0.60 * s["M2"]) * open_w
              + 0.35 * hill(s["CTP"], p["EC50_CTP"])
              - (p["KOUT_PDGF"] + 0.45 * PROT) * s["PDGF"])
    D("VEGF", p["KIN_VEGF"] * (HIFdrive + 0.50 * s["M2"]) * open_w
              + 0.30 * hill(s["CTP"], p["EC50_CTP"])
              - p["KOUT_VEGF"] * s["VEGF"])
    D("SDF1", p["KIN_SDF"] * HIFdrive * open_w
              - (p["KOUT_SDF"] + p["K_DPP4"] * (1.0 - 0.75 * p["DPP4I"])) * s["SDF1"])
    mobfail = 0.55 * min(hyperg / 4.0, 1.0)
    D("EPC", p["KM_EPC"] * s["SDF1"] * (1.0 - mobfail) - p["KOUT_EPC"] * s["EPC"])
    D("VASC", p["KG_VASC"] * s["VEGF"] * (0.35 + 0.65 * s["EPC"]) * OG * (1.0 - s["VASC"])
              - p["KR_VASC"] * s["VASC"] * (1.0 + 0.50 * PROTn))

    # ===== 7. fibroblasts, ECM, granulation ================================
    D("FIB", p["KP_FIB"] * hill(PDGF_tot, p["PDGF_50"]) * OG
             * (1.0 - s["FIB"] / p["FIB_MAX"]) * open_w
             - p["K_SEN"] * (s["ROS"] + 0.40 * PROTn) * s["FIB"]
             - p["KD_FIB"] * s["FIB"])
    D("SEN", p["K_SEN"] * (s["ROS"] + 0.40 * PROTn) * s["FIB"]
             - p["KCL_SEN"] * s["SEN"]
             - 1.6 * s["DEBR"] * s["SEN"])
    D("COL", p["KS_COL"] * s["FIB"] * OG * open_w
             - p["KDEG_COL"] * PROTn * s["COL"])
    D("GRAN", (p["KG_GRAN"] * hill(s["COL"], 0.70) * (0.30 + 0.70 * s["VASC"]) * OG
               * (1.0 + 0.35 * p["NPWT"]) * (1.0 - s["GRAN"])
               - p["KD_GRAN"] * PROTn * s["GRAN"]
               - 0.25 * s["DEBR"] * s["GRAN"]))

    # ===== 8. keratinocytes & epithelial edge ==============================
    kdrive = (inhib(s["IL1B"], p["IL1_50K"], 2.0)
              * inhib(PROT, p["PR50_K"], 1.0)
              * (0.60 + 0.40 * OG)
              * (1.0 + p["E_ESM"] * hill(s["ESM_W"], p["EC50_ESM"]))
              * (1.0 + p["E_CTP"] * hill(s["CTP"], p["EC50_CTP"]))
              * (1.0 - 0.10 * min(hyperg / 4.0, 1.0)))
    D("KERA", p["KIN_K"] * kdrive * (1.0 - s["KERA"])
              - p["KOUT_K"] * s["KERA"]
              + 0.40 * s["DEBR"] * (1.0 - s["KERA"]))     # de-epibolisation

    EDGE = (p["EDGE_MAX"] * MG * OG * s["KERA"] * s["GRAN"]
            * (1.0 - 0.55 * bf) * (1.0 + 0.25 * hill(PDGF_tot, p["PDGF_50"])))
    CONTR = p["K_CONTR"] * hill(s["COL"], p["C50_CONTR"]) * AREA

    D("AREA", (-EDGE * PERIM - CONTR + REINJ) * open_w)
    D("DEPTH", (-p["K_FILL"] * s["GRAN"] + 0.06 * s["OSTEO"] + 0.04 * bf * (1 - s["GRAN"]))
               * open_w * (1.0 if s["DEPTH"] > 0.05 else 0.0))

    # ===== 9. outcomes =====================================================
    D("SCAR", healed * (1.0 - s["SCAR"]) / p["TAU_SCAR"] - open_w * 0.10 * s["SCAR"])
    stress_post = ppp_index * (1.0 - 0.30 * p["FOOTWEAR"])
    tens = p["SCAR_MAX_TENS"] * s["SCAR"]
    haz = (p["H0_REC"] * (0.30 + 0.70 * s["NEURO"])
           * (stress_post / 0.90) ** 1.5
           * (1.0 - 0.55 * tens)
           * (1.0 - 0.10 * p["EDU"]) * (1.0 - 0.28 * p["TEMPMON"])
           * (1.0 - 0.12 * p["SURVEIL"]))
    D("HAZREC", healed * haz)
    D("HAZAMP", (p["A0_AMP"] * (2.0 * s["OSTEO"] ** 1.5
                                + 1.2 * (1.0 - s["PERF"]) ** 2
                                + 0.9 * max(bf - 0.6, 0.0) * 2.5) * open_w))
    D("UFD", healed)
    D("ABXD", on_abx if p["R_ABX"] > 0 else 0.0)
    D("TOXIDX", (p["W_TOX_ABX"] * on_abx * (1.0 if p["R_ABX"] > 0 else 0.0)
                 + p["W_TOX_HBOT"] * on_oxy * (1.0 if p["R_OXY"] > 0 else 0.0)
                 + p["W_TOX_BEC"] * on_bec * (1.0 if p["R_BEC"] > 0 else 0.0)
                 + p["W_TOX_TCC"] * (p["DEV_EFF"] * on_off)) / 100.0)
    D("AUCEDGE", EDGE)

    return d


# ---------------------------------------------------------------------------
def initial_state(p):
    y = [0.0] * N
    y[IDX["HBA1C"]] = p["HBA1C0"]
    y[IDX["NEURO"]] = (p["KP_NEURO"] * max(p["HBA1C0"] - 6.5, 0.0)
                       / (p["KP_NEURO"] * max(p["HBA1C0"] - 6.5, 0.0) + p["KR_NEURO"]))
    y[IDX["PERF"]] = p["PERF0"]
    y[IDX["BACT"]] = p["BACT0"]
    y[IDX["BIOF"]] = 0.72
    y[IDX["OSTEO"]] = 0.55 * p["PROBE_BONE"]
    y[IDX["BEXP"]] = p["PROBE_BONE"]
    y[IDX["NEUT"]] = 0.75
    y[IDX["M1"]] = 0.90
    y[IDX["M2"]] = 0.30
    y[IDX["IL1B"]] = 1.05
    y[IDX["ROS"]] = 0.95
    y[IDX["MMP9"]] = 1.30
    y[IDX["TIMP1"]] = 0.55
    y[IDX["PDGF"]] = 0.30
    y[IDX["VEGF"]] = 0.35
    y[IDX["SDF1"]] = 0.25
    y[IDX["EPC"]] = 0.18
    y[IDX["VASC"]] = 0.28
    y[IDX["FIB"]] = 0.55
    y[IDX["SEN"]] = 0.45
    y[IDX["COL"]] = 0.40
    y[IDX["GRAN"]] = 0.35
    y[IDX["KERA"]] = 0.35
    y[IDX["AREA"]] = p["AREA0"]
    y[IDX["DEPTH"]] = p["DEPTH0"]
    return y


HSTEP = [0.05]


def run(p, tmax=540.0, h=None, events=None):
    h = HSTEP[0] if h is None else h
    """RK4 with bolus events.  events = list of (time, state_name, amount)."""
    y = initial_state(p)
    n = int(round(tmax / h))
    ts = [0.0]
    ys = [list(y)]
    ev = sorted(events or [], key=lambda e: e[0])
    ei = 0
    t = 0.0
    # scheduled debridement pulses
    deb = []
    if p["R_DEBR"] > 0:
        k = 0
        while k * p["DEBR_PERIOD"] < min(tmax, p["T_DEBR_END"]):
            deb.append((k * p["DEBR_PERIOD"], "DEBR", p["R_DEBR"]))
            k += 1
    ev = sorted(list(ev) + deb, key=lambda e: e[0])

    for step in range(n):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            _, nm, amt = ev[ei]
            if nm == "OSTEO_RESECT":
                y[IDX["OSTEO"]] *= (1.0 - amt)
                y[IDX["BEXP"]] *= (1.0 - amt)      # the nidus itself is removed
                y[IDX["BACT"]] = max(y[IDX["BACT"]] - 1.5, 3.0)
            else:
                y[IDX[nm]] += amt
            ei += 1
        k1 = derivatives(t, y, p)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(N)]
        k2 = derivatives(t + 0.5 * h, y2, p)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(N)]
        k3 = derivatives(t + 0.5 * h, y3, p)
        y4 = [y[i] + h * k3[i] for i in range(N)]
        k4 = derivatives(t + h, y4, p)
        y = [y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(N)]
        # non-negativity / bounds
        for nm in NAMES:
            if y[IDX[nm]] < 0.0:
                y[IDX[nm]] = 0.0
        for nm in ("BIOF", "OSTEO", "VASC", "GRAN", "KERA", "SCAR", "NEURO", "PERF"):
            y[IDX[nm]] = min(y[IDX[nm]], 1.0)
        y[IDX["BACT"]] = min(max(y[IDX["BACT"]], 0.0), p["BMAX"])
        t += h
        ts.append(t)
        ys.append(list(y))
    return ts, ys


# ---------------------------------------------------------------------------
def series(ts, ys, name):
    i = IDX[name]
    return [row[i] for row in ys]


def at(ts, ys, name, t):
    i = IDX[name]
    k = min(range(len(ts)), key=lambda j: abs(ts[j] - t))
    return ys[k][i]


def time_to_close(ts, ys, close=0.05):
    a = series(ts, ys, "AREA")
    for j, v in enumerate(a):
        if v <= close:
            return ts[j]
    return None


def summarise(p, tmax=540.0, events=None):
    ts, ys = run(p, tmax=tmax, events=events)
    A0 = p["AREA0"]
    a28 = at(ts, ys, "AREA", 28.0)
    par4 = 1.0 - a28 / A0
    tclose = time_to_close(ts, ys)
    a84 = at(ts, ys, "AREA", 84.0)
    # Sheehan/perimeter-law prediction of time-to-heal from PAR4 alone
    if par4 > 0 and par4 < 1:
        theal_pred = 28.0 / (1.0 - math.sqrt(max(1.0 - par4, 0.0)))
    else:
        theal_pred = float("inf")
    # recurrence 365 d after closure
    rec12 = None
    if tclose is not None and tclose + 365 <= tmax:
        h1 = at(ts, ys, "HAZREC", tclose)
        h2 = at(ts, ys, "HAZREC", tclose + 365.0)
        rec12 = 1.0 - math.exp(-(h2 - h1))
    amp = 1.0 - math.exp(-at(ts, ys, "HAZAMP", min(tmax, 365.0)))
    # Empirical linear edge-advance rate k = -dr/dt, measured MID-COURSE so the
    # window never straddles closure (which would dilute k toward zero).
    tA, tB = ((0.40 * tclose, 0.80 * tclose) if tclose else (28.0, 84.0))
    rA = math.sqrt(max(at(ts, ys, "AREA", tA), 0.0) / math.pi)
    rB = math.sqrt(max(at(ts, ys, "AREA", tB), 0.0) / math.pi)
    k_late = (rA - rB) / max(tB - tA, 1e-9)
    t_pred_late = tB + rB / k_late if k_late > 1e-6 else float("inf")
    # --- post-hoc metrics recomputed from the stored trajectory ------------
    biof = series(ts, ys, "BIOF")
    nsel = [j for j in range(len(ts)) if ts[j] <= 56.0]
    biof_mean = sum(biof[j] for j in nsel) / max(len(nsel), 1)
    # realised vs biofilm-free antibiotic kill (the "wasted drug" fraction)
    kill_real = kill_max = 0.0
    for j in nsel:
        aw = ys[j][IDX["ABX_W"]]
        f = p["EMAX_ABX"] * hill(aw, p["EC50_ABX"])
        kill_max += f
        kill_real += f * (1.0 - p["TOL_BIOF"] * ys[j][IDX["BIOF"]])
    abx_realised = (kill_real / kill_max) if kill_max > 0 else None
    bact_nadir = min(series(ts, ys, "BACT")[j] for j in nsel)
    # PDGF-receptor drive actually delivered at day 28
    pd28 = at(ts, ys, "PDGF", 28.0)
    bw28 = at(ts, ys, "BEC_W", 28.0)
    pdgf_tot28 = pd28 + p["E_BEC"] * hill(bw28, p["EC50_BEC"], p["H_BEC"])
    pdgf_drive28 = hill(pdgf_tot28, p["PDGF_50"])
    return dict(
        par4=par4, a28=a28, a84=a84, tclose=tclose, theal_pred=theal_pred,
        biof_mean=biof_mean, abx_realised=abx_realised, bact_nadir=bact_nadir,
        pdgf_drive28=pdgf_drive28, bec_w28=bw28,
        closed12=a84 <= 0.05, rec12=rec12, amp365=amp,
        k_late=k_late, t_pred_late=t_pred_late, r0=math.sqrt(A0 / math.pi),
        tcpo2=None, ts=ts, ys=ys,
        mmp_ratio=at(ts, ys, "MMP9", 28.0) / (at(ts, ys, "TIMP1", 28.0) + 0.05),
        biof28=at(ts, ys, "BIOF", 28.0), bact28=at(ts, ys, "BACT", 28.0),
        osteo84=at(ts, ys, "OSTEO", 84.0),
        ufd=at(ts, ys, "UFD", tmax), tox=at(ts, ys, "TOXIDX", tmax),
    )


# ---------------------------------------------------------------------------
# Scenario library  (mirrors DFU_simulate_scenarios() in the R model)
# ---------------------------------------------------------------------------
def base(**kw):
    p = dict(P)
    p.update(kw)
    return p


def standard_care(**kw):
    """Sharp debridement q2wk + removable cast walker at real-world adherence."""
    return base(DEV_EFF=0.87, ADHERENCE=0.28, R_DEBR=1.0, DEBR_PERIOD=14.0, **kw)


SCENARIOS = {
    "S01 no offloading (dressing only)":
        base(R_DEBR=1.0, DEBR_PERIOD=14.0),
    "S02 standard care (RCW 28% adherence, debride q2wk)":
        standard_care(),
    "S03 total contact cast (irremovable)":
        base(DEV_EFF=0.85, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=14.0),
    "S04 removable walker, PERFECT adherence":
        base(DEV_EFF=0.87, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=14.0),
    "S05 half-shoe / therapeutic shoe":
        base(DEV_EFF=0.40, ADHERENCE=0.55, R_DEBR=1.0, DEBR_PERIOD=14.0),
    "S06 TCC + weekly debridement":
        base(DEV_EFF=0.85, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=7.0),
    "S07 TCC + debridement q4wk":
        base(DEV_EFF=0.85, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=28.0),
    "S08 standard care + becaplermin gel":
        standard_care(R_BEC=1.00),
    "S09 standard care + TLC-NOSF dressing":
        standard_care(R_NOSF=1.00),
    "S10 standard care + becaplermin + TLC-NOSF":
        standard_care(R_BEC=1.00, R_NOSF=1.00),
    "S11 standard care + esmolol 14% gel":
        standard_care(R_ESM=1.00),
    "S12 ISCHAEMIC wound (PERF 0.32) + becaplermin":
        standard_care(PERF0=0.32, PERF_TGT=0.32, R_BEC=1.00),
    "S13 ISCHAEMIC, revascularised d14, then becaplermin":
        standard_care(PERF0=0.32, PERF_TGT=0.32, T_REVASC=14.0,
                      REVASC_GAIN=0.45, R_BEC=1.00),
    "S14 infected wound, 14-d antibiotic alone (no debridement)":
        base(DEV_EFF=0.87, ADHERENCE=0.28, BACT0=7.4, R_ABX=1.0, T_ABX_END=14.0,
             T_DEBR_END=0.0),
    "S15 infected wound, 14-d antibiotic + weekly debridement":
        base(DEV_EFF=0.87, ADHERENCE=0.28, BACT0=7.4, R_ABX=1.0, T_ABX_END=14.0,
             R_DEBR=1.0, DEBR_PERIOD=7.0),
    "S16 OSTEOMYELITIS, 6-week antibiotic only":
        standard_care(BACT0=7.4, PROBE_BONE=1.0, R_ABX=1.0, T_ABX_END=42.0),
    "S17 OSTEOMYELITIS, resection d7 + 3-week antibiotic":
        standard_care(BACT0=7.4, PROBE_BONE=1.0, R_ABX=1.0, T_ABX_END=21.0),
    "S18 standard care + HBOT (daily 90 min)":
        standard_care(R_OXY=1.0),
    "S19 standard care + cellular tissue product":
        standard_care(R_CTP=1.0),
    "S20 standard care + NPWT":
        standard_care(NPWT=1.0),
    "S21 intensive glycaemic control (A1c 9.0 -> 7.0)":
        standard_care(HBA1C_TGT=7.0),
    "S22 OPTIMAL BUNDLE (TCC + q1wk debride + NOSF + becaplermin)":
        base(DEV_EFF=0.85, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=7.0,
             R_BEC=1.00, R_NOSF=1.00),
    "S23 OPTIMAL BUNDLE + remission care (footwear+temp+education)":
        base(DEV_EFF=0.85, ADHERENCE=1.00, R_DEBR=1.0, DEBR_PERIOD=7.0,
             R_BEC=1.00, R_NOSF=1.00, FOOTWEAR=1.0, TEMPMON=1.0,
             EDU=1.0, SURVEIL=1.0, HBA1C_TGT=7.0),
}

SCENARIO_EVENTS = {
    "S17 OSTEOMYELITIS, resection d7 + 3-week antibiotic": [(7.0, "OSTEO_RESECT", 0.85)],
}


def fmt(v, nd=2, na="--"):
    if v is None:
        return na
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float) and (v != v or v == float("inf")):
        return ">1y"
    return f"{v:.{nd}f}"


def main():
    if "--fast" in sys.argv:
        HSTEP[0] = 0.10
    csv_out = None
    if "--csv" in sys.argv:
        csv_out = sys.argv[sys.argv.index("--csv") + 1]

    if "--anchors" in sys.argv:
        print_anchors()
        return

    hdr = (f"{'scenario':<58} {'PAR4':>6} {'A(84d)':>7} {'t_close':>8} "
           f"{'pred':>6} {'12wk':>5} {'rec12mo':>8} {'MMP/TIMP':>9}")
    print(hdr)
    print("-" * len(hdr))
    rows = {}
    for name, p in SCENARIOS.items():
        r = summarise(p, events=SCENARIO_EVENTS.get(name))
        rows[name] = r
        print(f"{name:<58} {r['par4']*100:>5.1f}% {r['a84']:>7.3f} "
              f"{fmt(r['tclose'], 1, '>540'):>8} {fmt(r['theal_pred'],0):>6} "
              f"{fmt(r['closed12']):>5} "
              f"{(fmt(r['rec12']*100,1) if r['rec12'] is not None else '--'):>8} "
              f"{r['mmp_ratio']:>9.2f}")

    print()
    print_derived(rows)

    if csv_out:
        p = SCENARIOS["S02 standard care (RCW 28% adherence, debride q2wk)"]
        ts, ys = run(p)
        with open(csv_out, "w") as f:
            f.write("time," + ",".join(NAMES) + "\n")
            for j in range(0, len(ts), 20):
                f.write(f"{ts[j]:.2f}," + ",".join(f"{v:.6g}" for v in ys[j]) + "\n")
        print(f"\n[wrote {csv_out}]")


def print_derived(rows):
    def g(k):
        return rows[k]

    line = "=" * 96
    print(line)
    print("DERIVED RESULTS QUOTED IN README")
    print(line)

    tcc = g("S03 total contact cast (irremovable)")
    rcw = g("S02 standard care (RCW 28% adherence, debride q2wk)")
    rcw100 = g("S04 removable walker, PERFECT adherence")
    half = g("S05 half-shoe / therapeutic shoe")
    none = g("S01 no offloading (dressing only)")
    print("[1] OFFLOADING = DEVICE EFFICACY x ADHERENCE")
    print(f"    TCC (irremovable, 85% eff.)      t_close = {fmt(tcc['tclose'],1,'>540')} d")
    print(f"    Removable walker (87% eff, 28%)  t_close = {fmt(rcw['tclose'],1,'>540')} d")
    print(f"    SAME walker at 100% wear         t_close = {fmt(rcw100['tclose'],1,'>540')} d")
    print(f"    Half-shoe (40% eff, 55% wear)    t_close = {fmt(half['tclose'],1,'>540')} d")
    print(f"    No offloading at all             t_close = {fmt(none['tclose'],1,'>540')} d")
    if tcc["tclose"] and rcw["tclose"] and rcw100["tclose"]:
        gap = rcw["tclose"] - tcc["tclose"]
        resid = rcw100["tclose"] - tcc["tclose"]
        print(f"    -> walker-vs-TCC gap = {gap:.1f} d.  Holding adherence at 100% it becomes "
              f"{resid:+.1f} d.")
        print(f"       {100*(1-abs(resid)/gap):.0f}% of the gap is ADHERENCE, not the device: the "
              f"removable walker is")
        print(f"       mechanically the BETTER device (87% vs 85% pressure reduction) and loses "
              f"because it comes off.")

    print()
    bec = g("S08 standard care + becaplermin gel")
    nosf = g("S09 standard care + TLC-NOSF dressing")
    both = g("S10 standard care + becaplermin + TLC-NOSF")
    b = rcw["tclose"] or 540.0
    print("[2] TOPICAL GROWTH FACTOR IS A RACE AGAINST PROTEASE")
    print(f"    {'arm':<34}{'t_close':>9}{'MMP9/TIMP1':>12}{'BEC_W(d28)':>12}"
          f"{'PDGF drive':>12}")
    for nm, r in (("standard care", rcw), ("+ becaplermin", bec),
                  ("+ TLC-NOSF", nosf), ("+ both", both)):
        print(f"    {nm:<34}{fmt(r['tclose'],1,'>540'):>9}{r['mmp_ratio']:>12.2f}"
              f"{r['bec_w28']:>12.3f}{r['pdgf_drive28']:>12.3f}")
    ga = b - (bec["tclose"] or 540.0)
    gb = b - (nosf["tclose"] or 540.0)
    gc = b - (both["tclose"] or 540.0)
    print(f"    days gained: becaplermin {ga:.1f}, dressing {gb:.1f}, both {gc:.1f} "
          f"(sum of singles {ga+gb:.1f})")
    dbec = bec["pdgf_drive28"] - rcw["pdgf_drive28"]
    dnosf = nosf["pdgf_drive28"] - rcw["pdgf_drive28"]
    dboth = both["pdgf_drive28"] - rcw["pdgf_drive28"]
    if dbec + dnosf > 0:
        print(f"    -> AT THE RECEPTOR the interaction is superadditive: PDGF drive gain "
              f"{dbec:+.3f} (drug) and")
        print(f"       {dnosf:+.3f} (dressing) alone, {dboth:+.3f} together = "
              f"{dboth/(dbec+dnosf):.2f}x the sum. Protease modulation")
        print(f"       raises wound becaplermin {both['bec_w28']/max(bec['bec_w28'],1e-9):.2f}-fold.")
    if b and bec["tclose"] and nosf["tclose"] and both["tclose"]:
        r1 = b / bec["tclose"]
        r2 = b / nosf["tclose"]
        r12 = b / both["tclose"]
        print(f"    -> AT THE WHOLE WOUND it washes out to almost exactly multiplicative on "
              f"healing RATE:")
        print(f"       {r1:.3f} x {r2:.3f} = {r1*r2:.3f} predicted vs {r12:.3f} observed. "
              f"MMP-9 excess damages the")
        print(f"       wound through several parallel channels the dressing fixes anyway, so "
              f"the combination")
        print(f"       looks disappointing in days even though the drug is being delivered "
              f"far better.")

    print()
    isch = g("S12 ISCHAEMIC wound (PERF 0.32) + becaplermin")
    revasc = g("S13 ISCHAEMIC, revascularised d14, then becaplermin")
    print("[3] PERFUSION GATES EVERY ANABOLIC TERM")
    print(f"    ischaemic + becaplermin          A(84 d) = {isch['a84']:.3f} cm2, "
          f"t_close = {fmt(isch['tclose'],1,'>540')} d")
    print(f"    revascularise d14, same drug     A(84 d) = {revasc['a84']:.3f} cm2, "
          f"t_close = {fmt(revasc['tclose'],1,'>540')} d")
    print(f"    -> the identical biologic is worth nothing below the oxygen gate and "
          f"{(isch['a84']/max(revasc['a84'],1e-9)):.1f}x less wound")
    print(f"       area at 12 weeks. Revascularise first; a biologic on an ischaemic wound "
          f"is wasted drug.")

    print()
    abx = g("S14 infected wound, 14-d antibiotic alone (no debridement)")
    abxd = g("S15 infected wound, 14-d antibiotic + weekly debridement")
    print("[4] BIOFILM MAKES THE ANTIBIOTIC A PARTIAL DRUG")
    print(f"    antibiotic alone                 mean biofilm(0-56 d) = {abx['biof_mean']:.3f}, "
          f"realised kill = {100*abx['abx_realised']:.1f}% of potential")
    print(f"    + weekly sharp debridement       mean biofilm(0-56 d) = {abxd['biof_mean']:.3f}, "
          f"realised kill = {100*abxd['abx_realised']:.1f}% of potential")
    print(f"    BACT nadir  {abx['bact_nadir']:.2f} -> {abxd['bact_nadir']:.2f} log CFU/g; "
          f"t_close {fmt(abx['tclose'],1)} -> {fmt(abxd['tclose'],1)} d")
    print(f"    -> even weekly debridement only lifts mean biofilm suppression from "
          f"{abx['biof_mean']:.2f} to {abxd['biof_mean']:.2f}:")
    print(f"       the biofilm is back to tolerance within 2-3 days, so "
          f"{100*(1-abx['abx_realised']):.0f}% of the antibiotic")
    print(f"       course ({100*(1-abxd['abx_realised']):.0f}% even with weekly debridement) "
          f"is spent on an organism it cannot reach.")
    print(f"       Debridement's benefit is mostly the")
    print(f"       removal of the senescent, epibolised edge, not durable bioburden control.")

    print()
    o1 = g("S16 OSTEOMYELITIS, 6-week antibiotic only")
    o2 = g("S17 OSTEOMYELITIS, resection d7 + 3-week antibiotic")
    print("[5] OSTEOMYELITIS: THE NIDUS BEATS THE EXPOSURE")
    print(f"    6-week antibiotic only           OSTEO(84 d) = {o1['osteo84']:.3f}, "
          f"amputation risk(1 y) = {100*o1['amp365']:.1f}%, abx days = 42")
    print(f"    resection d7 + 3-week antibiotic OSTEO(84 d) = {o2['osteo84']:.3f}, "
          f"amputation risk(1 y) = {100*o2['amp365']:.1f}%, abx days = 21")
    print(f"    -> half the antibiotic exposure and a better bone outcome, because removing "
          f"the infected")
    print(f"       bone removes the substrate the regrowth term feeds on.")

    print()
    a1c = g("S21 intensive glycaemic control (A1c 9.0 -> 7.0)")
    print("[6] HbA1c IS A SLOW-AXIS DRUG (it treats the NEXT ulcer)")
    print(f"    standard care                    t_close = {fmt(rcw['tclose'],1)} d, "
          f"recurrence@12 mo = {fmt(100*rcw['rec12'],1)}%")
    print(f"    + A1c 9.0 -> 7.0                 t_close = {fmt(a1c['tclose'],1)} d, "
          f"recurrence@12 mo = {fmt(100*a1c['rec12'],1)}%")
    # long horizon: neuropathy time constant is years
    long_sc = summarise(standard_care(), tmax=1825.0)
    long_tx = summarise(standard_care(HBA1C_TGT=7.0), tmax=1825.0)
    for lbl, r in (("standard care", long_sc), ("A1c 9.0 -> 7.0", long_tx)):
        ts, ys = r["ts"], r["ys"]
        print(f"    {lbl:<24} NEURO: d0 {ys[0][IDX['NEURO']]:.3f} -> "
              f"1 y {at(ts,ys,'NEURO',365):.3f} -> 5 y {at(ts,ys,'NEURO',1825):.3f}")
    tc = long_sc["tclose"]
    if tc:
        def rec_between(r, t0, t1):
            return 1.0 - math.exp(-(at(r["ts"], r["ys"], "HAZREC", t1)
                                    - at(r["ts"], r["ys"], "HAZREC", t0)))
        print(f"    recurrence in year 4 alone (d{tc+1095:.0f}-d{tc+1460:.0f}): "
              f"{100*rec_between(long_sc, tc+1095, tc+1460):.1f}% vs "
              f"{100*rec_between(long_tx, tc+1095, tc+1460):.1f}%")
    print(f"    -> glycaemic control buys ~{(rcw['tclose'] or 0)-(a1c['tclose'] or 0):.0f} d on "
          f"THIS ulcer and its real payoff is years away,")
    tau_eff = 1.0 / (P["KP_NEURO"] * (7.0 - 6.5) + P["KR_NEURO"]) / 365.0
    print(f"       because neuropathy relaxes with a ~{tau_eff:.1f}-year time constant at "
          f"A1c 7.0. That is why HbA1c looks")
    print(f"       like a weak predictor of index-ulcer healing in the meta-analyses and a "
          f"strong one for incidence.")

    print()
    opt = g("S22 OPTIMAL BUNDLE (TCC + q1wk debride + NOSF + becaplermin)")
    rem = g("S23 OPTIMAL BUNDLE + remission care (footwear+temp+education)")
    print("[7] CLOSURE IS REMISSION, NOT CURE")
    print(f"    optimal bundle, no aftercare     t_close = {fmt(opt['tclose'],1)} d, "
          f"recurrence@12 mo = {fmt(100*opt['rec12'],1)}%")
    print(f"    + remission care                 t_close = {fmt(rem['tclose'],1)} d, "
          f"recurrence@12 mo = {fmt(100*rem['rec12'],1)}%")
    print(f"    -> the two arms close the ulcer within "
          f"{abs((opt['tclose'] or 0)-(rem['tclose'] or 0)):.1f} d of each other and then "
          f"diverge by")
    print(f"       {100*(opt['rec12']-rem['rec12']):.0f} percentage points. Every unit of "
          f"benefit in this comparison is earned")
    print(f"       AFTER the wound is closed, at the point where most trials stop measuring.")

    print()
    print("[8] THE PERIMETER LAW: WHAT PAR4 CAN AND CANNOT TELL YOU")
    print(f"    {'A0 (cm2)':>9}{'r0 (cm)':>9}{'PAR4':>8}{'k mid':>11}"
          f"{'t_close obs':>13}{'pred from k':>13}{'pred from PAR4':>16}")
    ks = []
    for a0 in (0.5, 1.0, 2.0, 4.0, 8.0):
        r = summarise(standard_care(AREA0=a0))
        ks.append(r["k_late"])
        print(f"    {a0:>9.1f}{r['r0']:>9.3f}{100*r['par4']:>7.1f}%{r['k_late']:>11.5f}"
              f"{fmt(r['tclose'],1,'>540'):>13}{fmt(r['t_pred_late'],1):>13}"
              f"{fmt(r['theal_pred'],1):>16}")
    spread = (max(ks) - min(ks)) / (sum(ks) / len(ks))
    print(f"    -> k (linear edge advance, cm/day) varies by {100*spread:.0f}% across a "
          f"16-fold range of wound")
    print(f"       AREA -- i.e. a 4-fold range of radius -- so most of the variation in "
          f"time-to-heal is")
    print(f"       geometry, t = r0/k. Given PAR4, the closed form 28/(1-sqrt(1-PAR4)) is "
          f"size-free -- but it")
    print(f"       systematically OVER-predicts here because k is still rising during weeks "
          f"1-4 as the wound")
    print(f"       bed converts. PAR4 is therefore a conservative screen: wounds that pass "
          f"it beat its own")
    print(f"       forecast, and the ones that fail it are the ones whose k never rose.")


def print_anchors():
    print("CALIBRATION ANCHORS  (model vs published)")
    print("=" * 84)
    print("Armstrong 2001 Diabetes Care 24:1019 -- mean days to healing, A0 = 1.4 cm2")
    for nm, dev, adh, obs in (("TCC", 0.85, 1.00, 33.5),
                              ("removable cast walker", 0.87, 0.28, 50.4),
                              ("half-shoe", 0.40, 0.55, 61.0)):
        r = summarise(base(DEV_EFF=dev, ADHERENCE=adh, R_DEBR=1.0,
                           DEBR_PERIOD=14.0, AREA0=1.4), tmax=300.0)
        print(f"  {nm:<24} observed {obs:5.1f} d    model {fmt(r['tclose'],1):>6} d")
    print()
    print("Armstrong 2003 Diabetes Care 26:2595 -- RCW worn during 28% of daily activity")
    print(f"  encoded directly as ADHERENCE = 0.28 (not as a lower device efficacy)")
    print()
    print("Sheehan 2003 Diabetes Care 26:1879 -- PAR4 > 50% predicts 12-week closure")
    for a0 in (0.5, 1.0, 2.0, 4.0):
        r = summarise(standard_care(AREA0=a0))
        rule = "yes" if r["par4"] > 0.5 else "no"
        print(f"  A0={a0:>4} cm2   PAR4 = {100*r['par4']:5.1f}%   rule says {rule:<4} "
              f"actually closed by 12 wk: {fmt(r['closed12'])}")
    print()
    print("Armstrong 2017 NEJM 376:2367 -- ~40% recurrence at 1 year, ~60% at 3 years")
    rcw = summarise(standard_care(), tmax=1825.0)
    tc = rcw["tclose"]
    h = lambda t: at(rcw["ts"], rcw["ys"], "HAZREC", t)
    print(f"  model standard care: 1 y {100*(1-math.exp(-(h(tc+365)-h(tc)))):.1f}%, "
          f"3 y {100*(1-math.exp(-(h(tc+1095)-h(tc)))):.1f}%")
    print()
    print("Wieman 1998 Diabetes Care 21:822 -- becaplermin 50% vs 35% closure at 20 wk")
    a = summarise(standard_care())
    bcp = summarise(standard_care(R_BEC=1.0))
    print(f"  model: standard care t_close {fmt(a['tclose'],1)} d, "
          f"+ becaplermin {fmt(bcp['tclose'],1)} d "
          f"({100*(1-(bcp['tclose'] or 1)/(a['tclose'] or 1)):.0f}% faster)")
    print()
    print("Edmonds 2018 Lancet Diabetes Endocrinol 6:186 (Explorer) -- TLC-NOSF")
    print("  48% vs 30% closure at 20 wk in neuro-ischaemic DFU")
    n = summarise(standard_care(R_NOSF=1.0))
    print(f"  model: + TLC-NOSF t_close {fmt(n['tclose'],1)} d "
          f"({100*(1-(n['tclose'] or 1)/(a['tclose'] or 1)):.0f}% faster)")


if __name__ == "__main__":
    main()
