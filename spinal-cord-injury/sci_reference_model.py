#!/usr/bin/env python3
"""
Traumatic Spinal Cord Injury (SCI) — QSP reference implementation
================================================================================
Dependency-free (standard library only) re-implementation of the 40-ODE system
in ``sci_mrgsolve_model.R``.  Its purpose is NOT to replace mrgsolve but to make
every quantitative claim in ``README.md`` reproducible without an R
installation, and to serve as a cross-check that the R model's equations do what
the text says they do.

    python3 sci_reference_model.py            # full report
    python3 sci_reference_model.py --quick    # scenarios only (skips sweeps)

The integrator is a fixed-step RK4 with a step schedule (fine during the acute
cascade, coarse in the chronic phase).  All equations, parameter names and
default values are kept character-for-character identical to the [PARAM]/[ODE]
blocks of the R model so the two can be diffed by eye.

Time unit: DAYS throughout (3 h = 0.125 d).  Author: Claude Code Routine (QSP
Disease Model Library).  Research/education only — not for clinical use.
"""

import math
import sys

# =============================================================================
# 1. PARAMETERS  (mirrors [PARAM] @annotated in sci_mrgsolve_model.R)
# =============================================================================

P = dict(
    # ---- patient / injury descriptors -------------------------------------
    LEVEL_IDX=5.0,     # C1-C8 = 1-8, T1-T12 = 9-20, L1-L5 = 21-25 (C5 default)
    AIS_INIT=1.0,      # 1=A, 2=B, 3=C, 4=D  (initial severity grade)
    DISRUPT=-1.0,      # override: mechanically disrupted tract fraction (<0 = from AIS)
    COMP0=0.60,        # residual mechanical cord compression (0-1) before surgery
    SUBLESIONAL_LMN=1.0,  # 1 = suprasacral (reflex arc intact), <1 = conus/cauda
    # ---- interventions (timing switches) ---------------------------------
    TDECOMP=1.0,       # d, time of surgical decompression
    DECOMP_EFF=0.95,   # fraction of compression relieved by surgery
    VASO=0.0,          # 1 = MAP >= 85-90 mmHg protocol active
    VASO_GAIN=20.0,    # mmHg, MAP support achieved by vasopressors
    TVASO_END=7.0,     # d, duration of the MAP protocol
    REHAB=0.0,         # 0-1 activity-based rehabilitation intensity
    TREHAB_START=14.0, # d
    TREHAB_END=180.0,  # d
    FES=0.0,           # 0-1 functional electrical stimulation / FES cycling
    ESTIM=0.0,         # 1 = epidural electrical stimulation active
    TESTIM_START=120.0,
    CHASE=0.0,         # 0-1 chondroitinase ABC effect (preclinical)
    ROCKI=0.0,         # 0-1 Rho-ROCK inhibitor effect
    ANTIMUSC=0.0,      # 0-1 antimuscarinic / onabotulinumtoxinA bladder effect
    # ---- methylprednisolone PK/PD (NASCIS II/III) ------------------------
    KOUT_MP=6.0,       # 1/d  (t1/2 ~ 2.8 h)
    KTR_MP=4.0,
    EC50_MP=8.0,       # au (mg/kg-equivalent central amount)
    EMAX_MP_LP=0.60,   # max fractional inhibition of lipid peroxidation
    EMAX_MP_CYT=0.50,  # max fractional cytokine-production suppression
    EMAX_MP_NEU=0.45,  # max fractional neutrophil-recruitment suppression
    MP_AUC50=40.0,     # au*d, cumulative exposure giving half-maximal complication index
    # ---- riluzole PK/PD (RISCIS) -----------------------------------------
    KA_RIL=8.0,
    KOUT_RIL=1.4,      # 1/d  (t1/2 ~ 12 h)
    EC50_RIL=45.0,
    EMAX_RIL_REL=0.55, # max fractional glutamate-release inhibition
    K_RIL_UPTAKE=0.60, # fractional gain in glutamate clearance
    # ---- minocycline PK/PD ----------------------------------------------
    KOUT_MINO=1.05,    # 1/d  (t1/2 ~ 16 h)
    EC50_MINO=400.0,
    EMAX_MINO_M1=0.40, # max fractional M1 activation blockade
    EMAX_MINO_CASP=0.25,
    K_MINO_SWITCH=0.50,
    # ---- glibenclamide PK/PD (SUR1-TRPM4) --------------------------------
    KOUT_GLY=2.5,
    EC50_GLY=0.5,
    EMAX_GLY=0.55,     # max fractional SUR1-TRPM4-dependent edema blockade
    # ---- anti-Nogo-A antibody (intrathecal, NISCI) -----------------------
    KOUT_NOGO=0.15,    # 1/d (CSF/tissue residence, t1/2 ~ 4.6 d)
    EC50_NOGO=0.60,
    EMAX_NOGO=0.80,    # max fractional Nogo-A neutralization
    # ---- baclofen / pregabalin (symptomatic) -----------------------------
    KOUT_BAC=4.0,
    EC50_BAC=25.0,
    EMAX_BAC=0.65,     # max fractional spasticity reduction
    K_BAC_WEAK=0.12,   # functional-weakness penalty scalar
    KOUT_PGB=2.7,
    EC50_PGB=150.0,
    EMAX_PGB=0.45,     # max fractional neuropathic-pain reduction
    # ---- hemodynamics ----------------------------------------------------
    MAP_BASE=90.0,
    SHOCK_MAG=25.0,    # mmHg MAP drop from neurogenic shock (lesion >= T6)
    K_SHOCK=0.20,      # 1/d resolution of neurogenic shock
    KOUT_MAP=12.0,
    ISP_BASE=7.0,      # mmHg
    K_ISP_EDEMA=18.0,  # mmHg per unit edema
    K_ISP_COMP=20.0,   # mmHg per unit residual compression
    SCPP50=65.0,       # mmHg, perfusion pressure at half-maximal ischemia
    SCPP_SLOPE=8.0,
    K_ISCH_COMP=0.40,  # direct (non-perfusion) ischemic effect of compression
    KOUT_ISCH=6.0,
    # ---- acute injury driver --------------------------------------------
    KOUT_INJ=2.0,      # 1/d decay of primary necrosis/DAMP release (t1/2 ~ 8.3 h)
    # ---- edema -----------------------------------------------------------
    K_ED_ISCH=0.60,
    K_ED_INJ=0.85,
    K_ED_INFL=0.10,
    KOUT_EDEMA=0.35,
    # ---- excitotoxicity --------------------------------------------------
    K_GLU_INJ=110.0,
    K_GLU_ISCH=45.0,
    KOUT_GLU=24.0,     # 1/d (t1/2 ~ 42 min)
    K_CA_GLU=1.00,
    KOUT_CA=8.0,
    # ---- oxidative stress ------------------------------------------------
    K_ROS_CA=0.90,
    K_ROS_ISCH=1.30,
    K_ROS_NEUT=0.80,
    K_ROS_IRON=1.20,
    KOUT_ROS=6.0,
    # ---- inflammation ----------------------------------------------------
    CYTO_REF=50.0,     # au, reference cytokine scale
    K_CYTO_DAMP=350.0,
    K_CYTO_M1=8.0,
    K_CYTO_NEUT=6.0,
    KOUT_CYTO=2.5,
    K_CYTO_M2=0.35,
    K_NEUT_IN=1.30,
    KOUT_NEUT=1.20,
    K_M1_DAMP=0.85,
    K_M1_CYTO=0.20,
    K_M1_DEBRIS=0.40,
    KOUT_M1=0.15,
    K_SWITCH=0.060,
    K_M2_SUPPRESS=0.35,  # M2-mediated resolution of the M1 phenotype
    K_M2_BASE=0.45,
    KOUT_M2=0.10,
    # ---- cell death ------------------------------------------------------
    K_AP_ROS=0.72,
    K_AP_CA=0.95,
    K_AP_CYTO=0.24,
    K_AP_ISCH=0.40,
    KOUT_APOP=0.50,
    APOP_THR=0.050,    # sub-threshold caspase activity is physiological turnover
    ISCH_THR=0.040,    # sub-threshold perfusion deficit is tolerated
    K_DEATH_N=0.030,
    K_NEC_ISCH=0.120,
    K_DEATH_O=0.050,
    K_OX_O=0.030,
    K_REMYEL=0.022,
    K_DEB_GEN=14.0,
    K_DEB_NEC=0.20,
    K_DEB_CLEAR=0.60,
    KOUT_DEB=0.020,
    DEB50=0.45,
    OLIG_MAX=0.92,
    K_AX_APOP=0.017,
    K_AX_ISCH=0.070,
    MYEL_FLOOR=0.35,
    K_ABLOCK=0.08,     # 1/d resolution of spinal shock / acute conduction block
    K_CB_ABLOCK=0.45,  # depth of the initial reversible conduction block
    K_CB_EDEMA=0.35,   # reversible conduction block from cord edema
    K_CB_ISCH=0.25,    # reversible conduction block from ischemia
    CB_FLOOR=0.25,
    # ---- scar ------------------------------------------------------------
    K_GFAP_CYTO=2.00,
    K_GFAP_INJ=1.50,
    KOUT_GFAP=0.012,
    K_CSPG_GFAP=0.012,
    KOUT_CSPG=0.015,
    CSPG50=0.55,       # CSPG level halving plasticity drive
    # ---- lesion cavity (imaging surrogate) -------------------------------
    K_CAV_APOP=0.035,
    K_CAV_NEUT=0.012,
    K_CAV_ISCH=0.030,
    KOUT_CAV=0.0004,
    CAV_MAX=1.0,       # lesion index saturates at full cross-sectional cavitation
    CAV_ML=1.8,        # mL of cavity per unit lesion index (reporting scale only)
    K_CAP_LOSS=0.45,   # permanent loss of remyelination capacity per unit oligo death
    # ---- plasticity ------------------------------------------------------
    K_CRIT=0.015,      # 1/d intrinsic closure of the critical period
    K_CRIT_SCAR=0.90,  # scar-accelerated closure
    K_PL_REHAB=0.055,
    K_PL_SPONT=0.010,
    K_PL_NOGO=1.10,
    K_PL_ROCK=0.60,
    KOUT_PLAST=0.0008,
    PLAST_MAX=1.0,
    G_PLAST=1.00,      # gain of plasticity on effective descending drive
    G_ESTIM=0.35,      # gain of epidural stimulation on residual drive
    # ---- muscle / motor mapping -----------------------------------------
    K_ATRO=0.020,
    K_ATRO_REC=0.030,
    ATRO_MAX=0.60,
    ATRO_FLOOR=0.10,   # atrophy that loading cannot reverse without voluntary drive
    K_ATRO_MOT=0.25,
    TH50=0.22,         # effective descending drive at half-maximal motor score
    HILL=2.2,
    KOUT_MOTOR=0.040,
    # ---- reflex / spasticity / pain / bladder ---------------------------
    K_REFLEX=0.050,
    ASH_MAX=2.8,
    KOUT_SPAST=0.15,
    K_NP_SENS=0.28,
    K_NP_ECTOPIC=1.30,
    K_NP_REFLEX=0.45,
    NP_MAX=9.0,
    KOUT_NP=0.08,
    BL_MAX=1.0,
    KOUT_BL=0.06,
    KOUT_ADT=24.0,
    K_AD_SBP=60.0,     # mmHg SBP surge at full trigger, mature reflex, >= T6
    # ---- respiratory / bone ---------------------------------------------
    KOUT_RESP=0.50,
    K_RESP_EDEMA=30.0,
    KOUT_BMD=0.0022,
    BMD_FLOOR=0.60,
)

# state ordering (identical to [CMT] in the R model)
NAMES = [
    "MP_CENT", "MP_PERIPH", "RIL_DEPOT", "RIL_CENT", "MINO_CENT", "GLY_CENT",
    "NOGO_ITH", "BAC_CENT", "PGB_CENT", "MP_AUC", "SHOCK", "MAP", "INJ",
    "ISCHEMIA", "EDEMA", "GLU", "CAI", "ROS", "NEUT", "CYTO", "M1", "M2",
    "APOP", "NEURON", "OLIG", "AXON", "GFAP", "CSPG", "CAVITY", "CRIT",
    "PLAST", "ATRO", "MOTOR", "REFLEX", "SPAST", "NPAIN", "BLADDER",
    "AD_TRIG", "RESP", "BMD", "DEBRIS", "OLIG_CAP", "ABLOCK",
]
IDX = {n: i for i, n in enumerate(NAMES)}
N = len(NAMES)

AIS_DISRUPT = {1: 0.94, 2: 0.86, 3: 0.70, 4: 0.45}


# =============================================================================
# 2. DERIVED CONSTANTS AND INITIAL CONDITIONS  (mirrors [MAIN])
# =============================================================================

def derived(p):
    """Level- and severity-dependent constants computed once per subject."""
    lvl = p["LEVEL_IDX"]
    d = {}
    d["above_t6"] = 1.0 if lvl <= 14.0 else 0.0        # T6 == index 14
    d["is_cervical"] = 1.0 if lvl <= 8.0 else 0.0
    # ISNCSCI points preserved rostral to the lesion (of 100 total)
    if lvl <= 4.0:
        d["above_score"] = 0.0
    elif lvl <= 8.0:
        d["above_score"] = (lvl - 4.0) * 10.0          # C5->10 ... C8->40
    else:
        d["above_score"] = 50.0                        # all upper-limb myotomes
    # respiratory penalty (percentage points of vital capacity)
    if lvl <= 4.0:
        d["resp_pen"] = 70.0
    elif lvl <= 8.0:
        d["resp_pen"] = 45.0
    elif lvl <= 14.0:
        d["resp_pen"] = 20.0
    else:
        d["resp_pen"] = 5.0
    d["isch_base"] = 1.0 / (1.0 + math.exp(
        (p["MAP_BASE"] - p["ISP_BASE"] - p["SCPP50"]) / p["SCPP_SLOPE"]))
    disrupt = p["DISRUPT"]
    if disrupt < 0:
        disrupt = AIS_DISRUPT[int(round(p["AIS_INIT"]))]
    d["disrupt"] = disrupt
    return d


def initial(p, d):
    y = [0.0] * N
    y[IDX["SHOCK"]] = d["above_t6"]
    y[IDX["MAP"]] = p["MAP_BASE"] - p["SHOCK_MAG"] * d["above_t6"]
    y[IDX["INJ"]] = d["disrupt"]
    y[IDX["NEURON"]] = 1.0 - 0.5 * d["disrupt"]
    y[IDX["OLIG"]] = p["OLIG_MAX"]
    y[IDX["OLIG_CAP"]] = p["OLIG_MAX"]
    y[IDX["AXON"]] = 1.0 - d["disrupt"]
    y[IDX["CRIT"]] = 1.0
    y[IDX["ABLOCK"]] = 1.0
    y[IDX["RESP"]] = 100.0 - d["resp_pen"] * (1.0 - 0.6 * hill_of(conn_static(y, p), p))
    y[IDX["BMD"]] = 1.0
    y[IDX["MOTOR"]] = motor_target_static(y, p, d)
    return y


def hill_of(conn, p):
    if conn <= 0:
        return 0.0
    c = conn ** p["HILL"]
    return c / (p["TH50"] ** p["HILL"] + c)


def conn_static(y, p):
    """Effective descending drive at t=0 (acute conduction block included)."""
    myel = p["MYEL_FLOOR"] + (1.0 - p["MYEL_FLOOR"]) * y[IDX["OLIG"]]
    cb = max(p["CB_FLOOR"], 1.0 - p["K_CB_ABLOCK"] * y[IDX["ABLOCK"]])
    return min(0.999, y[IDX["AXON"]] * myel * cb)


def motor_target_static(y, p, d):
    return (d["above_score"]
            + (100.0 - d["above_score"]) * hill_of(conn_static(y, p), p))


# =============================================================================
# 3. RIGHT-HAND SIDE  (mirrors [ODE])
# =============================================================================

def sw(t, t0, width=0.02):
    """Smooth 0->1 switch at t0 (keeps the system differentiable)."""
    z = (t - t0) / width
    if z > 30:
        return 1.0
    if z < -30:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def rhs(t, y, p, d, infusions):
    (MP_CENT, MP_PERIPH, RIL_DEPOT, RIL_CENT, MINO_CENT, GLY_CENT, NOGO_ITH,
     BAC_CENT, PGB_CENT, MP_AUC, SHOCK, MAP, INJ, ISCHEMIA, EDEMA, GLU, CAI,
     ROS, NEUT, CYTO, M1, M2, APOP, NEURON, OLIG, AXON, GFAP, CSPG, CAVITY,
     CRIT, PLAST, ATRO, MOTOR, REFLEX, SPAST, NPAIN, BLADDER, AD_TRIG, RESP,
     BMD, DEBRIS, OLIG_CAP, ABLOCK) = y

    dx = [0.0] * N

    # ---------- intervention switches ----------
    comp = p["COMP0"] * (1.0 - p["DECOMP_EFF"] * sw(t, p["TDECOMP"]))
    vaso_on = p["VASO"] * (1.0 - sw(t, p["TVASO_END"]))
    rehab_on = p["REHAB"] * sw(t, p["TREHAB_START"]) * (1.0 - sw(t, p["TREHAB_END"]))
    estim_on = p["ESTIM"] * sw(t, p["TESTIM_START"])
    fes_on = p["FES"] * sw(t, p["TREHAB_START"])

    # ---------- drug effect (Emax) terms ----------
    mp_lp = p["EMAX_MP_LP"] * MP_CENT / (p["EC50_MP"] + MP_CENT)
    mp_cyt = p["EMAX_MP_CYT"] * MP_CENT / (p["EC50_MP"] + MP_CENT)
    mp_neu = p["EMAX_MP_NEU"] * MP_CENT / (p["EC50_MP"] + MP_CENT)
    ril_rel = p["EMAX_RIL_REL"] * RIL_CENT / (p["EC50_RIL"] + RIL_CENT)
    ril_up = p["K_RIL_UPTAKE"] * RIL_CENT / (p["EC50_RIL"] + RIL_CENT)
    mino_m1 = p["EMAX_MINO_M1"] * MINO_CENT / (p["EC50_MINO"] + MINO_CENT)
    mino_casp = p["EMAX_MINO_CASP"] * MINO_CENT / (p["EC50_MINO"] + MINO_CENT)
    gly_blk = p["EMAX_GLY"] * GLY_CENT / (p["EC50_GLY"] + GLY_CENT)
    nogo_blk = p["EMAX_NOGO"] * NOGO_ITH / (p["EC50_NOGO"] + NOGO_ITH)
    bac_eff = p["EMAX_BAC"] * BAC_CENT / (p["EC50_BAC"] + BAC_CENT)
    pgb_eff = p["EMAX_PGB"] * PGB_CENT / (p["EC50_PGB"] + PGB_CENT)

    # ---------- drug PK ----------
    dx[IDX["MP_CENT"]] = -p["KOUT_MP"] * MP_CENT - p["KTR_MP"] * MP_CENT + p["KTR_MP"] * MP_PERIPH
    dx[IDX["MP_PERIPH"]] = p["KTR_MP"] * MP_CENT - p["KTR_MP"] * MP_PERIPH
    dx[IDX["RIL_DEPOT"]] = -p["KA_RIL"] * RIL_DEPOT
    dx[IDX["RIL_CENT"]] = p["KA_RIL"] * RIL_DEPOT - p["KOUT_RIL"] * RIL_CENT
    dx[IDX["MINO_CENT"]] = -p["KOUT_MINO"] * MINO_CENT
    dx[IDX["GLY_CENT"]] = -p["KOUT_GLY"] * GLY_CENT
    dx[IDX["NOGO_ITH"]] = -p["KOUT_NOGO"] * NOGO_ITH
    dx[IDX["BAC_CENT"]] = -p["KOUT_BAC"] * BAC_CENT
    dx[IDX["PGB_CENT"]] = -p["KOUT_PGB"] * PGB_CENT
    dx[IDX["MP_AUC"]] = MP_CENT           # cumulative exposure (au*d)

    # ---------- hemodynamics: MAP -> SCPP -> ischemia ----------
    dx[IDX["SHOCK"]] = -p["K_SHOCK"] * SHOCK
    map_target = p["MAP_BASE"] - p["SHOCK_MAG"] * SHOCK + p["VASO_GAIN"] * vaso_on
    dx[IDX["MAP"]] = p["KOUT_MAP"] * (map_target - MAP)

    isp = p["ISP_BASE"] + p["K_ISP_EDEMA"] * EDEMA + p["K_ISP_COMP"] * comp
    scpp = MAP - isp
    z = max(-30.0, min(30.0, (scpp - p["SCPP50"]) / p["SCPP_SLOPE"]))
    # ischemic burden is a DEFICIT relative to intact perfusion, so that an
    # uninjured cord (SCPP = MAP_BASE - ISP_BASE) sits at exactly zero
    isch_target = 1.0 / (1.0 + math.exp(z)) - d["isch_base"]
    isch_target = min(1.0, max(0.0, isch_target)) + p["K_ISCH_COMP"] * comp
    isch_target = min(1.0, isch_target)
    dx[IDX["ISCHEMIA"]] = p["KOUT_ISCH"] * (isch_target - ISCHEMIA)

    # ---------- acute injury driver ----------
    dx[IDX["INJ"]] = -p["KOUT_INJ"] * INJ
    dx[IDX["ABLOCK"]] = -p["K_ABLOCK"] * ABLOCK

    # ---------- edema (the compression-edema-ischemia vicious cycle) ----------
    ed_drive = (p["K_ED_ISCH"] * ISCHEMIA + p["K_ED_INJ"] * INJ
                + p["K_ED_INFL"] * (CYTO / p["CYTO_REF"])) * (1.0 - gly_blk)
    dx[IDX["EDEMA"]] = ed_drive * (1.0 - EDEMA) - p["KOUT_EDEMA"] * EDEMA

    # ---------- excitotoxicity ----------
    glu_prod = (p["K_GLU_INJ"] * INJ + p["K_GLU_ISCH"] * ISCHEMIA) * (1.0 - ril_rel)
    dx[IDX["GLU"]] = glu_prod - p["KOUT_GLU"] * GLU * (1.0 + ril_up)
    dx[IDX["CAI"]] = p["K_CA_GLU"] * GLU - p["KOUT_CA"] * CAI

    # ---------- oxidative stress ----------
    ros_prod = (p["K_ROS_CA"] * CAI + p["K_ROS_ISCH"] * ISCHEMIA
                + p["K_ROS_NEUT"] * NEUT + p["K_ROS_IRON"] * INJ)
    dx[IDX["ROS"]] = ros_prod * (1.0 - mp_lp) - p["KOUT_ROS"] * ROS

    # ---------- inflammation ----------
    cyto_prod = (p["K_CYTO_DAMP"] * INJ + p["K_CYTO_M1"] * M1
                 + p["K_CYTO_NEUT"] * NEUT) * (1.0 - mp_cyt)
    dx[IDX["CYTO"]] = cyto_prod - p["KOUT_CYTO"] * CYTO - p["K_CYTO_M2"] * M2 * CYTO
    dx[IDX["NEUT"]] = (p["K_NEUT_IN"] * (CYTO / p["CYTO_REF"]) * (1.0 - mp_neu)
                       - p["KOUT_NEUT"] * NEUT)
    m1_drive = (p["K_M1_DAMP"] * INJ + p["K_M1_CYTO"] * (CYTO / p["CYTO_REF"])
                + p["K_M1_DEBRIS"] * DEBRIS) * (1.0 - mino_m1)
    switch = p["K_SWITCH"] * M1 * (1.0 + p["K_MINO_SWITCH"] * mino_m1)
    dx[IDX["M1"]] = (m1_drive - p["KOUT_M1"] * M1 - switch
                     - p["K_M2_SUPPRESS"] * M2 * M1)
    dx[IDX["M2"]] = switch + p["K_M2_BASE"] * DEBRIS - p["KOUT_M2"] * M2

    # ---------- cell death ----------
    ap_drive = (p["K_AP_ROS"] * ROS + p["K_AP_CA"] * CAI
                + p["K_AP_CYTO"] * (CYTO / p["CYTO_REF"])
                + p["K_AP_ISCH"] * ISCHEMIA) * (1.0 - mino_casp)
    dx[IDX["APOP"]] = ap_drive - p["KOUT_APOP"] * APOP

    # only supra-threshold death signalling destroys tissue: without this the
    # model would drain axons for ever at the (tiny) chronic steady state
    ap_eff = max(0.0, APOP - p["APOP_THR"])
    isch_eff = max(0.0, ISCHEMIA - p["ISCH_THR"])

    # myelin/cell debris: generated by death, cleared by M2 phagocytosis
    death_flux = ((p["K_DEATH_O"] * ap_eff + p["K_OX_O"] * ROS) * OLIG
                  + (p["K_DEATH_N"] * ap_eff
                     + p["K_NEC_ISCH"] * isch_eff * isch_eff) * NEURON)
    dx[IDX["DEBRIS"]] = (p["K_DEB_GEN"] * death_flux
                         + p["K_DEB_NEC"] * isch_eff * isch_eff
                         - p["K_DEB_CLEAR"] * M2 * DEBRIS - p["KOUT_DEB"] * DEBRIS)
    dx[IDX["NEURON"]] = -(p["K_DEATH_N"] * ap_eff
                          + p["K_NEC_ISCH"] * isch_eff * isch_eff) * NEURON
    dx[IDX["OLIG"]] = (-(p["K_DEATH_O"] * ap_eff + p["K_OX_O"] * ROS) * OLIG
                       + p["K_REMYEL"] * (1.0 / (1.0 + DEBRIS / p["DEB50"]))
                       * max(0.0, OLIG_CAP - OLIG))
    # remyelination can never restore more than the surviving oligodendrocyte
    # lineage capacity: a permanent ceiling that only ever falls
    dx[IDX["OLIG_CAP"]] = -p["K_CAP_LOSS"] * (p["K_DEATH_O"] * ap_eff
                                              + p["K_OX_O"] * ROS) * OLIG
    dx[IDX["AXON"]] = -(p["K_AX_APOP"] * ap_eff
                        + p["K_AX_ISCH"] * isch_eff * isch_eff) * AXON

    # ---------- glial scar ----------
    dx[IDX["GFAP"]] = (p["K_GFAP_CYTO"] * (CYTO / p["CYTO_REF"])
                       + p["K_GFAP_INJ"] * INJ - p["KOUT_GFAP"] * GFAP)
    dx[IDX["CSPG"]] = (p["K_CSPG_GFAP"] * GFAP * (1.0 - p["CHASE"])
                       - p["KOUT_CSPG"] * CSPG)

    # ---------- lesion cavity (MRI surrogate) ----------
    cav_drive = (p["K_CAV_APOP"] * ap_eff + p["K_CAV_NEUT"] * NEUT
                 + p["K_CAV_ISCH"] * isch_eff * isch_eff)
    dx[IDX["CAVITY"]] = (cav_drive * (1.0 - CAVITY / p["CAV_MAX"])
                         - p["KOUT_CAV"] * CAVITY)

    # ---------- plasticity ----------
    dx[IDX["CRIT"]] = -p["K_CRIT"] * CRIT * (1.0 + p["K_CRIT_SCAR"] * CSPG)
    pl_drive = ((p["K_PL_REHAB"] * rehab_on + p["K_PL_SPONT"]) * CRIT
                / (1.0 + CSPG / p["CSPG50"])
                * (1.0 + p["K_PL_NOGO"] * nogo_blk + p["K_PL_ROCK"] * p["ROCKI"]))
    dx[IDX["PLAST"]] = (pl_drive * (1.0 - PLAST / p["PLAST_MAX"])
                        - p["KOUT_PLAST"] * PLAST)

    # ---------- effective descending drive & motor score ----------
    myel = p["MYEL_FLOOR"] + (1.0 - p["MYEL_FLOOR"]) * OLIG
    # reversible conduction block: the acute exam is worse than the anatomy
    cb = max(p["CB_FLOOR"], 1.0 - p["K_CB_ABLOCK"] * ABLOCK
             - p["K_CB_EDEMA"] * EDEMA - p["K_CB_ISCH"] * ISCHEMIA)
    conn = AXON * myel * cb * (1.0 + p["G_PLAST"] * PLAST + p["G_ESTIM"] * estim_on)
    conn = min(0.999, conn)
    hill = hill_of(conn, p)

    dx[IDX["ATRO"]] = (p["K_ATRO"] * max(0.0, 1.0 - rehab_on - fes_on - hill)
                       * (p["ATRO_MAX"] - ATRO)
                       - p["K_ATRO_REC"] * (rehab_on + fes_on)
                       * max(0.0, ATRO - p["ATRO_FLOOR"] * (1.0 - hill)))

    motor_target = (d["above_score"]
                    + (100.0 - d["above_score"]) * hill * (1.0 - p["K_ATRO_MOT"] * ATRO)
                    - p["K_BAC_WEAK"] * bac_eff * (100.0 - d["above_score"]))
    motor_target = max(d["above_score"] * 0.5, motor_target)
    dx[IDX["MOTOR"]] = p["KOUT_MOTOR"] * (motor_target - MOTOR)

    # ---------- reflex reorganization -> spasticity / bladder / dysreflexia ----
    reflex_target = p["SUBLESIONAL_LMN"] * max(0.0, 1.0 - 0.85 * hill)
    dx[IDX["REFLEX"]] = p["K_REFLEX"] * (reflex_target - REFLEX)
    spast_target = p["ASH_MAX"] * REFLEX * (1.0 - bac_eff)
    dx[IDX["SPAST"]] = p["KOUT_SPAST"] * (spast_target - SPAST)

    np_drive = (p["K_NP_SENS"] * (CYTO / p["CYTO_REF"] + 0.6 * M1)
                + p["K_NP_ECTOPIC"] * 4.0 * conn * (1.0 - conn)
                + p["K_NP_REFLEX"] * REFLEX)
    np_target = p["NP_MAX"] * (np_drive / (1.0 + np_drive)) * (1.0 - pgb_eff)
    dx[IDX["NPAIN"]] = p["KOUT_NP"] * (np_target - NPAIN)

    bl_target = p["BL_MAX"] * REFLEX * p["SUBLESIONAL_LMN"] * (1.0 - p["ANTIMUSC"])
    dx[IDX["BLADDER"]] = p["KOUT_BL"] * (bl_target - BLADDER)
    dx[IDX["AD_TRIG"]] = -p["KOUT_ADT"] * AD_TRIG

    # ---------- respiratory / bone ----------
    resp_target = (100.0 - d["resp_pen"] * (1.0 - 0.6 * hill)
                   - p["K_RESP_EDEMA"] * EDEMA * d["is_cervical"])
    dx[IDX["RESP"]] = p["KOUT_RESP"] * (resp_target - RESP)

    weight_bearing = min(1.0, hill + 0.5 * rehab_on + 0.4 * fes_on)
    dx[IDX["BMD"]] = -p["KOUT_BMD"] * (1.0 - weight_bearing) * (BMD - p["BMD_FLOOR"])

    # ---------- infusions (zero-order input) ----------
    for cmt, rate, t0, t1 in infusions:
        if t0 <= t < t1:
            dx[IDX[cmt]] += rate

    return dx


# =============================================================================
# 4. INTEGRATOR
# =============================================================================

def step_schedule(t):
    if t < 3.0:
        return 0.0005
    if t < 30.0:
        return 0.005
    return 0.02


def simulate(p_over=None, boluses=(), infusions=(), end=365.0, record=None):
    """RK4 with a step schedule. boluses = [(cmt, amt, time), ...]."""
    p = dict(P)
    if p_over:
        p.update(p_over)
    d = derived(p)
    y = initial(p, d)

    ev = sorted(boluses, key=lambda e: e[2])
    ei = 0
    rec_times = sorted(record) if record else []
    ri = 0
    out = {}
    t = 0.0

    # boluses at t = 0
    while ei < len(ev) and ev[ei][2] <= 0.0:
        y[IDX[ev[ei][0]]] += ev[ei][1]
        ei += 1
    if rec_times and rec_times[0] <= 0.0:
        out[rec_times[0]] = list(y)
        ri = 1

    while t < end - 1e-12:
        h = step_schedule(t)
        # do not step over the next bolus / recording time / infusion edge
        limits = [end]
        if ei < len(ev):
            limits.append(ev[ei][2])
        if ri < len(rec_times):
            limits.append(rec_times[ri])
        for _c, _r, t0, t1 in infusions:
            if t0 > t:
                limits.append(t0)
            if t1 > t:
                limits.append(t1)
        h = min(h, max(1e-9, min(limits) - t))

        k1 = rhs(t, y, p, d, infusions)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(N)]
        k2 = rhs(t + 0.5 * h, y2, p, d, infusions)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(N)]
        k3 = rhs(t + 0.5 * h, y3, p, d, infusions)
        y4 = [y[i] + h * k3[i] for i in range(N)]
        k4 = rhs(t + h, y4, p, d, infusions)
        y = [y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(N)]
        # non-negativity guard (mirrors mrgsolve's practical behaviour)
        for i in range(N):
            if y[i] < 0.0 and NAMES[i] not in ("MAP", "MOTOR", "RESP"):
                y[i] = 0.0
        t += h

        while ei < len(ev) and ev[ei][2] <= t + 1e-12:
            y[IDX[ev[ei][0]]] += ev[ei][1]
            ei += 1
        while ri < len(rec_times) and rec_times[ri] <= t + 1e-12:
            out[rec_times[ri]] = list(y)
            ri += 1

    res = {"p": p, "d": d, "y_end": y, "rec": out}
    res["conn_end"] = conn_of(y, p, 0.0)
    res["hill_end"] = hill_of(res["conn_end"], p)
    return res


def conn_of(y, p, estim_on=0.0):
    myel = p["MYEL_FLOOR"] + (1.0 - p["MYEL_FLOOR"]) * y[IDX["OLIG"]]
    c = y[IDX["AXON"]] * myel * (1.0 + p["G_PLAST"] * y[IDX["PLAST"]]
                                 + p["G_ESTIM"] * estim_on)
    return min(0.999, c)


def get(res, name, t=None):
    y = res["y_end"] if t is None else res["rec"][t]
    return y[IDX[name]]


# =============================================================================
# 5. REGIMEN BUILDERS  (identical dosing to the R scenarios)
# =============================================================================

def mp_regimen(start_h, hours=23.0):
    """NASCIS: 30 mg/kg bolus then 5.4 mg/kg/h for 23 h (or 47 h)."""
    t0 = start_h / 24.0
    inf_amt = 5.4 * hours
    return ([("MP_CENT", 30.0, t0)],
            [("MP_CENT", inf_amt / (hours / 24.0), t0, t0 + hours / 24.0)])


def riluzole_regimen(start_h=12.0, days=14.0):
    """100 mg q12h for 24 h then 50 mg q12h to day 14 (RISCIS)."""
    b = []
    t = start_h / 24.0
    b.append(("RIL_DEPOT", 100.0, t))
    b.append(("RIL_DEPOT", 100.0, t + 0.5))
    k = 2
    while t + 0.5 * k <= days:
        b.append(("RIL_DEPOT", 50.0, t + 0.5 * k))
        k += 1
    return b


def minocycline_regimen(start_h=12.0, days=7.0):
    b = []
    t = start_h / 24.0
    k = 0
    while t + 0.5 * k <= days:
        b.append(("MINO_CENT", 800.0 if k < 4 else 400.0, t + 0.5 * k))
        k += 1
    return b


def glibenclamide_regimen(start_h=6.0, days=3.0):
    b = []
    t = start_h / 24.0
    k = 0
    while t + 0.25 * k <= days:
        b.append(("GLY_CENT", 1.0, t + 0.25 * k))
        k += 1
    return b


def antinogo_regimen(start_d=14.0, n=6, interval=7.0):
    return [("NOGO_ITH", 1.0, start_d + i * interval) for i in range(n)]


def baclofen_regimen(start_d=30.0, end_d=365.0):
    b, t = [], start_d
    while t <= end_d:
        b.append(("BAC_CENT", 20.0, t))
        t += 1.0 / 3.0
    return b


def pregabalin_regimen(start_d=30.0, end_d=365.0):
    b, t = [], start_d
    while t <= end_d:
        b.append(("PGB_CENT", 150.0, t))
        t += 0.5
    return b


# =============================================================================
# 6. SCENARIOS
# =============================================================================

BASE = dict(LEVEL_IDX=5.0, AIS_INIT=1.0, COMP0=0.60, TDECOMP=1.0)


def scenario(name, over=None, boluses=(), infusions=(), end=365.0):
    p = dict(BASE)
    if over:
        p.update(over)
    res = simulate(p, boluses=list(boluses), infusions=list(infusions), end=end,
                   record=[0.0, 1.0, 3.0, 7.0, 14.0, 30.0, 90.0, 180.0, 365.0])
    res["name"] = name
    return res


def build_scenarios():
    s = []
    mp3_b, mp3_i = mp_regimen(3.0)
    mp9_b, mp9_i = mp_regimen(9.0)
    mp48_b, mp48_i = mp_regimen(3.0, hours=47.0)

    s.append(scenario("1. Natural history (AIS A C5, no decompression, no rehab)",
                      dict(TDECOMP=999.0, REHAB=0.0)))
    s.append(scenario("2. Standard care: decompression 24 h + rehab 0.5",
                      dict(TDECOMP=1.0, REHAB=0.5)))
    s.append(scenario("3. Early decompression < 12 h + rehab",
                      dict(TDECOMP=0.5, REHAB=0.5)))
    s.append(scenario("4. Late decompression 72 h + rehab",
                      dict(TDECOMP=3.0, REHAB=0.5)))
    s.append(scenario("5. MP within 3 h (NASCIS II) + std care",
                      dict(REHAB=0.5), mp3_b, mp3_i))
    s.append(scenario("6. MP started at 9 h (outside window) + std care",
                      dict(REHAB=0.5), mp9_b, mp9_i))
    s.append(scenario("7. MP 48-h infusion (NASCIS III) + std care",
                      dict(REHAB=0.5), mp48_b, mp48_i))
    s.append(scenario("8. Riluzole 14 d (RISCIS) + std care",
                      dict(REHAB=0.5), riluzole_regimen()))
    s.append(scenario("9. Minocycline 7 d + std care",
                      dict(REHAB=0.5), minocycline_regimen()))
    s.append(scenario("10. Glibenclamide 3 d (SUR1-TRPM4) + std care",
                      dict(REHAB=0.5), glibenclamide_regimen()))
    s.append(scenario("11. MAP >= 85 mmHg x 7 d + std care",
                      dict(REHAB=0.5, VASO=1.0)))
    s.append(scenario("12. Bundle: decompression 8 h + MAP protocol + riluzole",
                      dict(TDECOMP=8.0 / 24.0, VASO=1.0, REHAB=0.5),
                      riluzole_regimen(start_h=8.0)))
    s.append(scenario("13. Anti-Nogo-A + intensive rehab from d14",
                      dict(REHAB=1.0, TDECOMP=1.0), antinogo_regimen()))
    s.append(scenario("14. Intensive rehab EARLY (d14) vs late - early arm",
                      dict(REHAB=1.0, TREHAB_START=14.0)))
    s.append(scenario("15. Intensive rehab LATE (d120)",
                      dict(REHAB=1.0, TREHAB_START=120.0, TREHAB_END=286.0)))
    s.append(scenario("16. Symptomatic: baclofen + pregabalin from d30",
                      dict(REHAB=0.5),
                      baclofen_regimen() + pregabalin_regimen()))
    s.append(scenario("17. Incomplete AIS C, standard care + rehab",
                      dict(AIS_INIT=3.0, REHAB=0.5)))
    s.append(scenario("18. AIS C + epidural stimulation from d120",
                      dict(AIS_INIT=3.0, REHAB=0.5, ESTIM=1.0)))
    return s


# =============================================================================
# 7. REPORT
# =============================================================================

def hdr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def report_scenarios(scen):
    hdr("A. SCENARIO SUMMARY (day 365 unless noted)")
    print(f"{'scenario':<58}{'MOTOR':>7}{'CAVITY':>8}{'AXON':>7}{'CONN':>7}{'ASH':>6}{'NRS':>6}")
    print("-" * 99)
    for r in scen:
        print(f"{r['name'][:57]:<58}"
              f"{get(r,'MOTOR'):>7.1f}{get(r,'CAVITY'):>8.3f}"
              f"{get(r,'AXON'):>7.3f}{r['conn_end']:>7.3f}"
              f"{get(r,'SPAST'):>6.2f}{get(r,'NPAIN'):>6.2f}")

    hdr("B. ACUTE-PHASE TRAJECTORY (scenario 2, standard care)")
    r = scen[1]
    print(f"{'day':>6}{'MAP':>7}{'ISCH':>7}{'EDEMA':>7}{'GLU':>7}{'ROS':>7}"
          f"{'CYTO':>8}{'NEUT':>7}{'M1':>7}{'APOP':>7}{'AXON':>7}{'MOTOR':>7}")
    for t in [0.0, 1.0, 3.0, 7.0, 14.0, 30.0, 90.0, 365.0]:
        y = r["rec"][t]
        print(f"{t:>6.0f}{y[IDX['MAP']]:>7.1f}{y[IDX['ISCHEMIA']]:>7.3f}"
              f"{y[IDX['EDEMA']]:>7.3f}{y[IDX['GLU']]:>7.2f}{y[IDX['ROS']]:>7.2f}"
              f"{y[IDX['CYTO']]:>8.1f}{y[IDX['NEUT']]:>7.2f}{y[IDX['M1']]:>7.2f}"
              f"{y[IDX['APOP']]:>7.2f}{y[IDX['AXON']]:>7.3f}{y[IDX['MOTOR']]:>7.1f}")


def report_window(quick=False, ais=1.0, tag="AIS A"):
    hdr(f"C. THERAPEUTIC WINDOW: METHYLPREDNISOLONE START TIME ({tag} C5)")
    ref = scenario("ref", dict(REHAB=0.5, AIS_INIT=ais))
    m_ref = get(ref, "MOTOR")
    cav_ref = get(ref, "CAVITY")
    print(f"reference (no MP): MOTOR {m_ref:.2f}, CAVITY {cav_ref:.3f}")
    print(f"\n{'start (h)':>10}{'MOTOR':>9}{'dMOTOR':>9}{'CAVITY':>9}"
          f"{'dCAV %':>9}{'MP_AUC':>9}{'compl.idx':>11}")
    rows = []
    for h in ([1, 3, 6, 9, 12, 24] if not quick else [3, 9]):
        b, i = mp_regimen(float(h))
        r = scenario(f"MP@{h}h", dict(REHAB=0.5, AIS_INIT=ais), b, i)
        auc = get(r, "MP_AUC")
        ci = auc / (P["MP_AUC50"] + auc)
        dcav = 100.0 * (get(r, "CAVITY") - cav_ref) / cav_ref
        rows.append((h, get(r, "MOTOR"), get(r, "MOTOR") - m_ref,
                     get(r, "CAVITY"), dcav, auc, ci))
        print(f"{h:>10}{rows[-1][1]:>9.2f}{rows[-1][2]:>+9.2f}"
              f"{rows[-1][3]:>9.3f}{dcav:>+9.1f}{auc:>9.1f}{ci:>11.3f}")
    if len(rows) > 2:
        first, last = rows[0], rows[-1]
        print(f"\n  -> benefit at {first[0]} h = {first[2]:+.2f} motor points; "
              f"at {last[0]} h = {last[2]:+.2f}")
        print(f"  -> the complication index is IDENTICAL "
              f"({first[6]:.3f} vs {last[6]:.3f}) because it depends on dose, "
              f"not timing")
        # 48-h regimen
        b, i = mp_regimen(3.0, hours=47.0)
        r48 = scenario("MP48", dict(REHAB=0.5, AIS_INIT=ais), b, i)
        auc48 = get(r48, "MP_AUC")
        ci48 = auc48 / (P["MP_AUC50"] + auc48)
        print(f"  -> 48-h regimen started at 3 h: MOTOR {get(r48,'MOTOR'):.2f} "
              f"({get(r48,'MOTOR')-m_ref:+.2f}), complication index {ci48:.3f} "
              f"({ci48/rows[1][6]:.2f}x the 24-h regimen)")
    return rows


def report_decompression(quick=False, ais=1.0, tag="AIS A"):
    hdr(f"D. THERAPEUTIC WINDOW: TIME TO SURGICAL DECOMPRESSION ({tag} C5)")
    print(f"{'t_decomp (h)':>13}{'MOTOR':>9}{'AXON':>8}{'CAVITY':>9}"
          f"{'peak ISCH':>11}{'peak EDEMA':>12}")
    out = []
    for h in ([4, 8, 12, 24, 48, 72, 999] if not quick else [8, 24, 72]):
        p = dict(BASE); p.update(dict(REHAB=0.5, TDECOMP=h / 24.0, AIS_INIT=ais))
        rec = [i * 0.25 for i in range(0, 25)] + [365.0]
        r = simulate(p, end=365.0, record=rec)
        pk_i = max(r["rec"][t][IDX["ISCHEMIA"]] for t in rec)
        pk_e = max(r["rec"][t][IDX["EDEMA"]] for t in rec)
        y = r["rec"][365.0]
        out.append((h, y[IDX["MOTOR"]], y[IDX["AXON"]], y[IDX["CAVITY"]], pk_i, pk_e))
        lbl = "never" if h == 999 else f"{h}"
        print(f"{lbl:>13}{y[IDX['MOTOR']]:>9.2f}{y[IDX['AXON']]:>8.3f}"
              f"{y[IDX['CAVITY']]:>9.3f}{pk_i:>11.3f}{pk_e:>12.3f}")
    if len(out) > 3:
        print(f"\n  -> 4 h vs 24 h: {out[0][1]-out[3][1]:+.2f} motor points; "
              f"24 h vs 72 h: {out[3][1]-out[5][1]:+.2f}")
    return out


def report_threshold():
    hdr("E. SURROGATE-ENDPOINT DISSOCIATION: MOTOR = f(CONN) IS A STEEP HILL")
    p = dict(P)
    print(f"TH50 = {p['TH50']}, HILL = {p['HILL']}\n")
    print(f"{'CONN':>7}{'hill':>8}{'MOTOR (C5, above=10)':>24}"
          f"{'d MOTOR per +0.01 CONN':>25}")
    for c in [0.02, 0.03, 0.05, 0.08, 0.11, 0.15, 0.22, 0.30, 0.45, 0.60, 0.80]:
        h1 = hill_of(c, p)
        h2 = hill_of(c + 0.01, p)
        print(f"{c:>7.2f}{h1:>8.3f}{10 + 90*h1:>24.1f}{90*(h2-h1):>25.2f}")
    print("\n  A 2-fold improvement in spared drive buys very different amounts of")
    print("  function depending on where the patient starts:")
    for c0 in [0.03, 0.06, 0.11, 0.20]:
        m0 = 10 + 90 * hill_of(c0, p)
        m1 = 10 + 90 * hill_of(2 * c0, p)
        print(f"    CONN {c0:.2f} -> {2*c0:.2f} (+100%):  MOTOR "
              f"{m0:5.1f} -> {m1:5.1f}  ({m1-m0:+.1f} points)")
    print("\n  ISNCSCI total-motor measurement noise is ~5 points, so a doubling of")
    print("  spared drive is UNDETECTABLE in the most severe patients and large in")
    print("  the mid-range - the arithmetic behind enrolling AIS B/C, not AIS A.")


def report_severity_interaction():
    hdr("F. SAME DRUG, DIFFERENT GRADE: RILUZOLE ACROSS BASELINE SEVERITY")
    print(f"{'AIS':>5}{'no drug MOTOR':>15}{'riluzole MOTOR':>16}"
          f"{'dMOTOR':>9}{'dCAVITY %':>11}")
    for ais, lbl in [(1, "A"), (2, "B"), (3, "C"), (4, "D")]:
        r0 = scenario("x", dict(AIS_INIT=float(ais), REHAB=0.5))
        r1 = scenario("y", dict(AIS_INIT=float(ais), REHAB=0.5),
                      riluzole_regimen())
        dcav = 100.0 * (get(r1, "CAVITY") - get(r0, "CAVITY")) / get(r0, "CAVITY")
        print(f"{lbl:>5}{get(r0,'MOTOR'):>15.2f}{get(r1,'MOTOR'):>16.2f}"
              f"{get(r1,'MOTOR')-get(r0,'MOTOR'):>+9.2f}{dcav:>+11.1f}")
    print("\n  The tissue-level effect (cavity) is nearly grade-independent; the")
    print("  motor-score effect is not. A trial powered on AIS A patients can be")
    print("  mechanistically successful and clinically negative at the same time.")


def report_two_clocks(ais=3.0, tag="AIS C"):
    hdr(f"G. TWO CLOCKS: NEUROPROTECTION (h) vs PLASTICITY (weeks-months) [{tag}]")
    early = scenario("early", dict(REHAB=1.0, AIS_INIT=ais,
                                   TREHAB_START=14.0, TREHAB_END=180.0))
    late = scenario("late", dict(REHAB=1.0, AIS_INIT=ais,
                                 TREHAB_START=120.0, TREHAB_END=286.0))
    none = scenario("none", dict(REHAB=0.0, AIS_INIT=ais))
    print(f"{'arm':<34}{'PLAST':>8}{'CRIT(d365)':>12}{'CONN':>8}{'MOTOR':>8}")
    for lbl, r in [("no rehab", none), ("rehab d14-180 (same dose)", early),
                   ("rehab d120-286 (same dose)", late)]:
        print(f"{lbl:<34}{get(r,'PLAST'):>8.3f}{get(r,'CRIT'):>12.4f}"
              f"{r['conn_end']:>8.3f}{get(r,'MOTOR'):>8.2f}")
    print(f"\n  Identical rehabilitation dose, 106 days apart: "
          f"{get(early,'MOTOR')-get(late,'MOTOR'):+.2f} motor points.")
    print("  The window that closes here is the CSPG/scar-gated critical period,")
    print("  not the acute cascade - a second clock with its own units.")


def report_complications():
    hdr("H. RECOVERY CREATES NEW PATHOLOGY: REFLEX RETURN")
    r = scenario("std", dict(REHAB=0.5, LEVEL_IDX=12.0),  # T4 lesion
                 end=365.0)
    print("T4 AIS A lesion, standard care:")
    print(f"{'day':>6}{'REFLEX':>9}{'ASH':>7}{'NRS':>7}{'BLADDER':>9}{'MOTOR':>8}")
    for t in [1.0, 7.0, 14.0, 30.0, 90.0, 180.0, 365.0]:
        y = r["rec"][t]
        print(f"{t:>6.0f}{y[IDX['REFLEX']]:>9.3f}{y[IDX['SPAST']]:>7.2f}"
              f"{y[IDX['NPAIN']]:>7.2f}{y[IDX['BLADDER']]:>9.3f}"
              f"{y[IDX['MOTOR']]:>8.1f}")
    # autonomic dysreflexia challenge at day 7 vs day 180
    hdr("I. AUTONOMIC DYSREFLEXIA CHALLENGE (bladder distension, T4 lesion)")
    print(f"{'trigger day':>13}{'REFLEX':>9}{'MAP':>8}{'peak SBP':>10}{'surge':>8}")
    for tday in [7.0, 30.0, 90.0, 180.0]:
        p = dict(BASE); p.update(dict(LEVEL_IDX=12.0, REHAB=0.5))
        rec = [tday + i * 0.005 for i in range(0, 40)]
        rr = simulate(p, boluses=[("AD_TRIG", 1.0, tday)], end=tday + 1.0,
                      record=rec)
        best = 0.0
        refl = 0.0
        mapv = 0.0
        for t in rec:
            y = rr["rec"][t]
            sbp = 1.15 * y[IDX["MAP"]] + P["K_AD_SBP"] * y[IDX["AD_TRIG"]] * \
                rr["d"]["above_t6"] * y[IDX["REFLEX"]] * (1 + 0.5 * y[IDX["BLADDER"]])
            if sbp > best:
                best, refl, mapv = sbp, y[IDX["REFLEX"]], y[IDX["MAP"]]
        print(f"{tday:>13.0f}{refl:>9.3f}{mapv:>8.1f}{best:>10.1f}"
              f"{best-1.15*mapv:>8.1f}")
    print("\n  The same noxious stimulus is harmless in the flaccid first week and")
    print("  a hypertensive emergency once the reflex arc has reorganized.")


def report_symptomatic():
    hdr("J. SYMPTOMATIC PHARMACOLOGY IS NOT FREE (baclofen trade-off)")
    r0 = scenario("no sympt", dict(REHAB=0.5))
    r1 = scenario("baclofen", dict(REHAB=0.5), baclofen_regimen())
    r2 = scenario("bac+pgb", dict(REHAB=0.5),
                  baclofen_regimen() + pregabalin_regimen())
    print(f"{'arm':<24}{'ASH':>7}{'NRS':>7}{'MOTOR':>8}")
    for lbl, r in [("none", r0), ("baclofen 20 mg TID", r1),
                   ("+ pregabalin 300 mg/d", r2)]:
        print(f"{lbl:<24}{get(r,'SPAST'):>7.2f}{get(r,'NPAIN'):>7.2f}"
              f"{get(r,'MOTOR'):>8.2f}")
    print("\n  Ashworth falls and pain falls, but the modelled motor score falls too:")
    print("  the antispastic effect and the weakness share one mechanism.")


def report_organ_systems():
    hdr("K. ORGAN-SYSTEM CONSEQUENCES AT 1 YEAR")
    rows = [
        ("C4 AIS A, no rehab", dict(LEVEL_IDX=4.0, REHAB=0.0)),
        ("C4 AIS A, rehab+FES", dict(LEVEL_IDX=4.0, REHAB=1.0, FES=1.0)),
        ("C5 AIS A, no rehab", dict(LEVEL_IDX=5.0, REHAB=0.0)),
        ("C5 AIS A, rehab+FES", dict(LEVEL_IDX=5.0, REHAB=1.0, FES=1.0)),
        ("T10 AIS A, rehab+FES", dict(LEVEL_IDX=18.0, REHAB=1.0, FES=1.0)),
    ]
    print(f"{'arm':<26}{'RESP %VC':>10}{'BMD':>8}{'ATRO':>8}{'MOTOR':>8}")
    for lbl, over in rows:
        o = dict(REHAB=0.5); o.update(over)
        r = scenario(lbl, o)
        print(f"{lbl:<26}{get(r,'RESP'):>10.1f}{get(r,'BMD'):>8.3f}"
              f"{get(r,'ATRO'):>8.3f}{get(r,'MOTOR'):>8.1f}")
    print("\n  Sublesional BMD loss is driven by absent weight-bearing, so it is the")
    print("  one endpoint that FES/standing improves without any change in CONN.")


def main():
    quick = "--quick" in sys.argv
    print(__doc__.split("\n\n")[0])
    print("Traumatic Spinal Cord Injury QSP — reference implementation report")
    print(f"{N} ODE states, RK4, step schedule "
          f"0.0005/0.005/0.02 d (acute/subacute/chronic)")

    scen = build_scenarios()
    report_scenarios(scen)
    report_threshold()
    report_window(quick)
    report_window(quick, ais=3.0, tag="AIS C")
    report_decompression(quick)
    report_decompression(quick, ais=3.0, tag="AIS C")
    report_severity_interaction()
    report_two_clocks()
    report_complications()
    report_symptomatic()
    report_organ_systems()
    print("\n" + "=" * 78)
    print("END OF REPORT — research/education only, not for clinical use.")
    print("=" * 78)


if __name__ == "__main__":
    main()
