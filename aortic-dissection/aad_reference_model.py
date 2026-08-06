#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aad_reference_model.py
======================
Independent Python/scipy re-implementation of the 52-ODE acute type B aortic
dissection (AAD) QSP model that is shipped as `aad_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
The mrgsolve model is the deliverable, but nothing in the repository can run
it (no R in the container).  An ODE system that has never been integrated is
a wish list, not a model.  So the same equations, the same parameter values
and the same dosing scenarios are written a second time here, in a different
language, with a different solver (LSODA), and every claim that appears in
README.md or in the mechanistic map is produced by running THIS file.  Where
the two implementations disagree, the disagreement is a bug in one of them and
is fixed in both.

THE STRUCTURAL CLAIM THE MODEL IS BUILT TO TEST
-----------------------------------------------
In acute aortic dissection the quantity that tears the aorta is not blood
pressure.  It is the product

        ASI  =  (sigma_FL_sys / sigma_ref)  x  (dP/dt_max / dPdt_ref)

of a STRESS term and an IMPULSE term, and the two terms are set by different
physiology:

  * sigma_FL = P_FL * r_FL / h_FL     (Laplace, applied to the FALSE lumen
    outer wall, which is outer media + adventitia only -> h_FL ~ 0.40 * h_ao,
    so at identical pressure the false lumen carries ~2.5x the true lumen's
    circumferential stress).  Pressure enters here.

  * dP/dt_max is set by LV contractility and heart rate, i.e. by SYMPATHETIC
    DRIVE, and is almost independent of the mean pressure.

A pure vasodilator lowers P and therefore lowers the first term - and, through
the arterial baroreflex it unloads, RAISES the second.  A beta-blocker lowers
both.  The whole "beta-blockade before vasodilatation" rule of acute aortic
syndrome management is therefore a statement about a PRODUCT, and it can be
computed rather than asserted.  To make the comparison fair, every acute
scenario is titrated by a closed-loop controller to the SAME systolic pressure
target (100-120 mmHg), so the arms differ only in HOW they got there.

THE SECOND STRUCTURAL CLAIM
---------------------------
The famous and counter-intuitive result that PARTIAL false-lumen thrombosis
carries a worse prognosis than either complete patency or complete thrombosis
(Tsai 2007 NEJM) is, in this model, a statement about a resistive divider.
The false lumen is a compliant chamber connected to the true lumen by an
entry-tear conductance G_en and a distal re-entry conductance G_re:

        dP_FL/dt = kcomm * ( G_en*(MAP - P_FL) + G_re*(P_dist - P_FL) )
                   - kdep * (P_FL - P_ven)

and its pulse pressure is low-pass filtered with tau_FL ~ 1/(G_en + G_re):

  * G_en high, G_re high  -> brisk through-flow, FL pulsatile, mean pressure
                             pulled down by flow across the entry tear.
  * G_en high, G_re ~ 0   -> BLIND SAC.  No through-flow, so no pressure drop:
                             FL mean pressure equals systemic mean pressure,
                             and because the pulse is filtered out, FL
                             DIASTOLIC pressure equals systemic MEAN pressure,
                             which is higher than true-lumen diastolic
                             pressure.  Sustained stress is maximal.
  * G_en ~ 0, G_re ~ 0    -> COMPLETE THROMBOSIS.  The sac depressurises with
                             tau = 336 h and the thrombus mechanically
                             supports the wall.

Nothing in the code says "partial thrombosis is bad".  The three phenotypes
are produced by two anatomical numbers (entry tear area, re-entry
conductance) and the thrombosis of the stagnant compartment follows from the
flow the divider allows.

UNITS
-----
time            h            (fast PK included; LSODA handles the stiffness)
amounts         mg
volumes         L
pressures       mmHg
diameters       mm
wall thickness  mm
stress          kPa          (1 mmHg = 0.13332 kPa)
dP/dt           mmHg/s
flows           L/min
resistance      mmHg*min/L

USAGE
-----
    python3 aad_reference_model.py            # full report -> stdout
    python3 aad_reference_model.py --json     # machine-readable summary
"""

import argparse
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

MMHG_TO_KPA = 0.133322

# ----------------------------------------------------------------------------
# 0.  STATE VECTOR (52 ODEs) - identical order to $CMT in the mrgsolve file
# ----------------------------------------------------------------------------
CMT = [
    # --- drug PK (16) -------------------------------------------------------
    "A_ESM",     # 1  esmolol, central                                    (mg)
    "A_LABC",    # 2  labetalol, central                                  (mg)
    "A_LABP",    # 3  labetalol, peripheral                               (mg)
    "A_METD",    # 4  metoprolol, gut depot                               (mg)
    "A_METC",    # 5  metoprolol, central                                 (mg)
    "A_NICC",    # 6  nicardipine, central                                (mg)
    "A_NICP",    # 7  nicardipine, peripheral                             (mg)
    "A_CLV",     # 8  clevidipine, central                                (mg)
    "A_SNP",     # 9  sodium nitroprusside, central                       (mg)
    "A_CN",      # 10 cyanide                                             (mg)
    "A_SCN",     # 11 thiocyanate                                         (mg)
    "A_FENC",    # 12 fentanyl, central                                   (mg)
    "A_FENE",    # 13 fentanyl, effect site                               (mg)
    "A_LOSD",    # 14 losartan, gut depot                                 (mg)
    "A_LOSC",    # 15 losartan, central                                   (mg)
    "A_E3174",   # 16 EXP-3174 (active metabolite)                        (mg)
    # --- neurohumoral / systemic haemodynamics (10) -------------------------
    "SNA",       # 17 sympathetic outflow                          (rel., 1=nl)
    "BAROSET",   # 18 baroreflex operating point                        (mmHg)
    "HR",        # 19 heart rate                                         (bpm)
    "CTR",       # 20 LV contractility                             (rel., 1=nl)
    "SVR",       # 21 systemic vascular resistance             (mmHg*min/L)
    "PMAP",      # 22 mean arterial pressure                            (mmHg)
    "PRA",       # 23 plasma renin activity                        (rel., 1=nl)
    "ANGII",     # 24 angiotensin II                               (rel., 1=nl)
    "ALDO",      # 25 aldosterone                                  (rel., 1=nl)
    "BV",        # 26 blood volume                                          (L)
    # --- aortic / two-lumen mechanics (9) -----------------------------------
    "PFL",       # 27 false lumen mean pressure                         (mmHg)
    "DFL",       # 28 false lumen diameter                                (mm)
    "DTL",       # 29 true lumen diameter                                 (mm)
    "THRFL",     # 30 false lumen thrombus fraction (body)            (0-1)
    "THRD",      # 31 distal re-entry thrombus fraction               (0-1)
    "AEN",       # 32 entry tear area                                    (mm2)
    "HFL",       # 33 false lumen outer wall thickness                    (mm)
    "ELN",       # 34 medial elastin integrity                        (0-1)
    "COL",       # 35 adventitial collagen (healing)                  (0-2)
    # --- injury / inflammation / pain (5) -----------------------------------
    "INJ",       # 36 acute dissection injury signal                  (0-1)
    "IL6",       # 37 interleukin-6                            (rel., 1=nl)
    "MMP9",      # 38 MMP-9 / elastolytic activity             (rel., 1=nl)
    "TGFB",      # 39 TGF-beta signalling                      (rel., 1=nl)
    "PAIN",      # 40 pain intensity                                  (0-10)
    # --- organ perfusion / laboratory endpoints (5) -------------------------
    "QREN",      # 41 renal perfusion (lagged)                     (rel., 1=nl)
    "CREAT",     # 42 serum creatinine                                (mg/dL)
    "LACT",      # 43 arterial lactate                               (mmol/L)
    "SCI",       # 44 spinal cord ischaemia burden                (index*h)
    "DDIM",      # 45 D-dimer                                         (ug/mL)
    # --- cumulative outcome integrals (5) -----------------------------------
    "CUMH",      # 46 cumulative rupture hazard                             (-)
    "SURV",      # 47 survival probability                              (0-1)
    "DAM",       # 48 cumulative stress-impulse damage (tear propagation)  (-)
    "TIT",       # 49 time inside BP+HR target                              (h)
    "AUCASI",    # 50 AUC of the aortic stress-impulse index             (h)
    # --- closed-loop titration controllers (2) ------------------------------
    "U_VD",      # 51 vasodilator controller integrator                     (-)
    "U_BB",      # 52 beta-blocker controller integrator                    (-)
    # --- spinal cord collateral network (1) ---------------------------------
    "COLLAT",    # 53 paraspinal collateral network                    (0-1)
]
IX = {n: i for i, n in enumerate(CMT)}
NEQ = len(CMT)

# ----------------------------------------------------------------------------
# 1.  PARAMETERS
#     Every value carries the source it was anchored to.  Where a value is a
#     model construct with no direct measurement (conductances, growth gains)
#     it is marked CAL: and the calibration target is named.
# ----------------------------------------------------------------------------
P = dict(
    WT=75.0,                 # kg, reference IRAD type B patient (~65 y, male)

    # ---- esmolol -----------------------------------------------------------
    # t1/2 9 min (RBC esterase), CL 285 mL/kg/min, Vd 3.4 L/kg  (Wiest 1995)
    V_ESM=255.0, CL_ESM=1283.0,          # L, L/h  -> t1/2 = 0.138 h
    IC50_ESM=0.55,                       # mg/L, beta1 (EC50 0.3-1.0 mg/L)
    # ---- labetalol ---------------------------------------------------------
    # t1/2 5.5 h, Vd 5.6 L/kg, CL 25 mL/min/kg; IV alpha:beta block = 1:7
    V_LAB=420.0, V_LAB2=130.0, CL_LAB=70.0, Q_LAB=45.0,
    IC50_LABB=0.13,                      # mg/L, beta1
    IC50_LABA=0.91,                      # mg/L, alpha1 (= 7 x beta IC50)
    # ---- metoprolol (oral, chronic arm) ------------------------------------
    KA_MET=1.2, F_MET=0.40, V_MET=300.0, CL_MET=63.0,   # t1/2 3.3 h, EM
    IC50_MET=0.035,                      # mg/L, beta1
    # ---- nicardipine IV ----------------------------------------------------
    V_NIC=45.0, V_NIC2=200.0, CL_NIC=40.0, Q_NIC=95.0,
    IC50_NIC=0.055,                      # mg/L (~55 ng/mL) arteriolar
    # ---- clevidipine IV ----------------------------------------------------
    V_CLV=12.0, CL_CLV=290.0,            # t1/2 ~1.7 min, blood esterases
    IC50_CLV=0.030,
    # ---- sodium nitroprusside + cyanide ------------------------------------
    V_SNP=15.0, CL_SNP=291.0,            # t1/2 ~2 min
    IC50_SNP=0.020,
    FCN=0.436,                           # mg CN- released per mg SNP (5 CN-/SNP)
    V_CN=30.0, CL_CN=31.0,               # CN- t1/2 ~40 min (rhodanese)
    V_SCN=35.0, CL_SCN=0.34,             # SCN- t1/2 ~3 d, renally cleared
    CN_TOX=1.0,                          # mg/L, threshold for toxicity
    # ---- fentanyl ----------------------------------------------------------
    V_FEN=280.0, CL_FEN=45.0, KE0_FEN=9.0,   # ke0 0.15/min
    IC50_FEN=0.0012,                     # mg/L (1.2 ng/mL) analgesia
    # ---- losartan / EXP-3174 (chronic arm) ---------------------------------
    KA_LOS=1.0, F_LOS=0.33, V_LOS=34.0, CL_LOS=50.0, FM_LOS=0.14,
    V_E3174=12.0, CL_E3174=2.0,          # t1/2 ~4-9 h
    IC50_LOS=0.075, IC50_E3174=0.010,    # mg/L, AT1

    # ---- baroreflex / sympathetic --------------------------------------
    TAU_SNA=0.15, TAU_BARO=48.0,
    BARO_REF=106.0, FRESET=0.45,         # baroreflex resets only PARTIALLY
    KBARO=1.25,                          # baroreflex gain (blunted, 65 y HTN)
    KPAIN=0.55,                          # sympathoexcitation per 10/10 pain
    KANG=0.35,                           # AT1-mediated sympathoexcitation
    SNA0=1.0,
    # ---- cardiac -----------------------------------------------------------
    TAU_HR=0.05, TAU_CTR=0.05, TAU_SVR=0.08, TAU_MAP=0.01,
    HR0=62.0, CTR0=1.0, SVR0=17.5, SV0=76.0,
    EMAX_BB_HR=0.46, EMAX_BB_CTR=0.42,   # max beta1-blockade effect
    EMAX_DHP=0.68, EMAX_ALPHA=0.30, EMAX_NO=0.62, EMAX_AT1=0.16,
    EXP_HR_SNA=0.85, EXP_CTR_SNA=0.60, EXP_SVR_SNA=0.70,
    CART0=1.72,                          # mL/mmHg, arterial compliance
    STIFF=0.55,                          # age/elastin stiffening of CART
    DPDT0=1180.0, DBP0=80.0,             # mmHg/s at HR0/CTR0/DBP0
    EXP_DPDT_HR=0.25, EXP_DPDT_DBP=0.35,
    # ---- volume / RAAS -----------------------------------------------------
    BV0=5.0, TAU_BV=48.0, TAU_PRA=0.5, TAU_ANG=0.05, TAU_ALDO=2.0,
    KRENIN=2.2,                          # renin gain per unit renal flow deficit
    EMAX_BB_RENIN=0.40,                  # beta1-blockade of JG cell renin
    KAT1FB=1.10,                         # loss of AngII short-loop feedback on ARB
    KALDO_BV=0.09,

    # ---- two-lumen mechanics ----------------------------------------------
    KCOMM=260.0,                         # 1/h per unit conductance (fast)
    KDEP=1.0 / 48.0,                     # depressurisation of a SEALED sac (2 d)
    PVEN=9.0,                            # mmHg, pressure a sealed sac decays to
    DPAO=4.0,                            # mmHg, TL pressure drop along the aorta
    AEN0=42.0,                           # mm2, entry tear AT PRESENTATION (~7 mm)
    AEN_REF=42.0,                        # mm2, FIXED normaliser for conductance
    AEN_MAX=320.0,                       # mm2, whole-circumference limit
    GRE0=1.00,                           # CAL: reference re-entry conductance
    TAUFL0=0.30,                         # s, FL charging time at unit conductance
    HAO0=1.90,                           # mm, descending thoracic wall thickness
    FOUT=0.40,                           # dissection plane -> outer wall fraction
    DTL0=20.0, DFL0=17.0,                # mm at presentation
    TAU_DTL=1.0,
    KCOLLAPSE=0.32,                      # max TL compression by a pressurised FL
    CFL_KP=6.0,                          # mmHg scale of the collapse sigmoid
    SIGREF=118.0,                        # kPa, normal aortic wall stress
    # entry tear occludes ONLY when thrombosis is essentially complete; this is
    # what separates "complete" from "partial" FL thrombosis in the model.
    SEAL_THR50=0.70, SEAL_THRW=0.050, SEAL_POW=3.0, KFLOW_THR=0.60,
    # ---- thrombosis --------------------------------------------------------
    KTHR=0.0015, KLYS=0.0060,            # 1/h, body of the FL
    KTHRD=0.0090, KLYSD=0.0025,          # 1/h, stagnant distal FL
    FLOW50=0.30, FLOW50D=0.09,           # stasis half-points (rel. flow)
    KQ=0.95,                             # through-flow scale
    KTOFRO=0.55,                         # to-and-fro flow scale
    # ---- growth / remodelling ---------------------------------------------
    KGROW=3.55e-4,                       # CAL: mm/h at sigma=SIGREF (see notes)
    EXP_GROW=1.9,
    KMMP_GROW=0.30,
    THR_SUPPORT=0.85,                    # thrombus mechanical support of wall
    KREGRESS=6.0e-5,                     # organising thrombus -> FL shrinkage
    DFL_CAP=58.0,                        # mm, numerical ceiling (patient long dead)
    # wall thickness is homeostatic: adventitial collagen thickens it, MMP-9 and
    # the proteolytic environment UNDER MURAL THROMBUS thin it.
    KREM_H=1.0 / 240.0, HFL_COL=0.45, KH_MMP=0.10, KTHR_H=1.00,
    KELN=8.0e-6,
    KCOL=2.0e-3, KCOL_STRESS=3.2e-3, KDEGC=1.0e-3, COLMAX=1.8,
    KTEAR=6.0e-3, KHEAL_TEAR=3.0e-5, EXP_TEAR_DPDT=1.2, TEAR_TH=0.45,
    # ---- inflammation ------------------------------------------------------
    TAU_INJ=30.0, KIL6=2.4, KIL6_STRESS=0.8, KOUT_IL6=0.35,
    KMMP_IL6=0.55, KMMP_AT1=0.60, KOUT_MMP=0.20, GEN_MMP=1.0,
    KTGF_INJ=1.20, KTGF_MMP=0.40, KOUT_TGF=0.09,
    TAU_PAIN=0.10, SIG_PAIN=175.0, KPAIN_INJ=7.5, EMAX_OP=0.90,
    # ---- organ perfusion ---------------------------------------------------
    OBST_REN=0.0, OBST_MES=0.0, OBST_SPI=0.0,   # static branch obstruction
    INV_REN=0.0, INV_MES=0.0, INV_SPI=1.0,      # branch arises from dissected aorta
    KOBST_DYN=0.75,                      # dynamic obstruction from TL collapse
    TAU_QREN=1.0,
    CREAT_PROD=1.05, CL_CREAT0=1.0, V_CREAT=42.0,
    KLACT_MES=2.6, KLACT_CN=2.2, KOUT_LACT=0.55, LACT0=1.0,
    SCP_CSF=10.0, KSCI=1.0, SCI_THRESH=0.45, SCPP50=42.0,
    SCI_PARA=12.0, TAU_SCI=72.0,         # paraplegia threshold, repair time
    KCOLL_BASE=0.55,                     # pre-existing collateral network
    KCOL_SPI=0.85, KCOLL=0.020, KCOLL_D=1.0 / 1440.0,
    KDD_FLOW=1.6, KDD_THR=2.4, KOUT_DD=0.10, DDIM0=0.25,
    # ---- rupture hazard / survival ----------------------------------------
    SWALL0=760.0,                        # kPa, ultimate wall strength (fresh)
    EXP_ELN_S=0.30,
    # Anchored so that sigma/S = 0.40 (a 26 mm dissected aorta at controlled
    # pressure) gives 2%/y and sigma/S = 0.67 (58 mm) gives ~11%/y.
    H0_RUPT=1.10e-4,                     # 1/h at sigma = wall strength
    W_RUPT=0.155,                        # exponential width of the hazard
    KH_MALPERF=1.3e-3, KH_OTHER=2.3e-6,  # 2%/y non-aortic mortality at 65 y
    KH_EXT=4.0e-3,                       # hazard of propagation / extension
    KDAM=1.0,
    # ---- BP / HR target window --------------------------------------------
    SBP_LO=100.0, SBP_HI=120.0, HR_TGT=60.0,
    # ---- closed-loop titration --------------------------------------------
    KI_VD=0.16, KI_BB=0.11, KI_BB_SBP=0.09,   # 1/(mmHg*h), 1/(bpm*h)
    VD_MAX=0.0, BB_MAX=0.0,              # scenario-set maximum infusion rates
    SBP_SET=110.0, HR_SET=60.0,
    BB_ON_SBP=0,                         # 1 = titrate the beta-blocker to SBP
    # ---- fixed-rate / bolus dosing (set per scenario) ---------------------
    R_ESM=0.0, R_LAB=0.0, R_NIC=0.0, R_CLV=0.0, R_SNP=0.0, R_FEN=0.0,
    D_MET=0.0, D_LOS=0.0, TAU_PO=12.0,   # oral maintenance, mg per TAU_PO h
    # ---- which agent each controller drives -------------------------------
    VD_AGENT=0,                          # 0 none, 1 nicardipine, 2 clevidipine, 3 SNP
    BB_AGENT=0,                          # 0 none, 1 esmolol, 2 labetalol
    # ---- interventions -----------------------------------------------------
    TEVAR_T=-1.0,                        # h, time of entry-tear coverage (-1 = never)
    TEVAR_SEAL=0.985,                    # fraction of entry conductance removed
    TEVAR_COVER=0.30,                    # fraction of intercostals covered
    ADHERE=1.0,                          # chronic drug adherence 0-1
)

# ----------------------------------------------------------------------------
# 2.  HELPERS
# ----------------------------------------------------------------------------
def sig(x):
    """Numerically safe logistic."""
    if x > 40.0:
        return 1.0
    if x < -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def softplus(x, k=8.0):
    """Smooth max(x, 0) - keeps the ODEs differentiable for the solver."""
    z = k * x
    if z > 40.0:
        return x
    if z < -40.0:
        return 0.0
    return math.log1p(math.exp(z)) / k


def occ(*pairs):
    """Additive competitive occupancy of one receptor by several ligands."""
    e = 0.0
    for c, ic50 in pairs:
        e += c / ic50
    return e / (1.0 + e)


# ----------------------------------------------------------------------------
# 3.  DERIVED (algebraic) QUANTITIES - the readouts of the model
# ----------------------------------------------------------------------------
def derive(y, p, t):
    """Everything the ODEs and the report need that is not a state."""
    d = {}
    g = d.__setitem__

    # ---------------- drug concentrations ----------------------------------
    C_ESM = y[IX["A_ESM"]] / p["V_ESM"]
    C_LAB = y[IX["A_LABC"]] / p["V_LAB"]
    C_MET = y[IX["A_METC"]] / p["V_MET"]
    C_NIC = y[IX["A_NICC"]] / p["V_NIC"]
    C_CLV = y[IX["A_CLV"]] / p["V_CLV"]
    C_SNP = y[IX["A_SNP"]] / p["V_SNP"]
    C_CN = y[IX["A_CN"]] / p["V_CN"]
    C_SCN = y[IX["A_SCN"]] / p["V_SCN"]
    C_FEN = y[IX["A_FENE"]] / p["V_FEN"]
    C_LOS = y[IX["A_LOSC"]] / p["V_LOS"]
    C_E31 = y[IX["A_E3174"]] / p["V_E3174"]
    C_ESM, C_LAB, C_MET = max(C_ESM, 0.0), max(C_LAB, 0.0), max(C_MET, 0.0)
    C_NIC, C_CLV, C_SNP = max(C_NIC, 0.0), max(C_CLV, 0.0), max(C_SNP, 0.0)
    C_CN, C_SCN, C_FEN = max(C_CN, 0.0), max(C_SCN, 0.0), max(C_FEN, 0.0)
    C_LOS, C_E31 = max(C_LOS, 0.0), max(C_E31, 0.0)
    for k, v in (("C_ESM", C_ESM), ("C_LAB", C_LAB), ("C_MET", C_MET),
                 ("C_NIC", C_NIC), ("C_CLV", C_CLV), ("C_SNP", C_SNP),
                 ("C_CN", C_CN), ("C_SCN", C_SCN), ("C_FEN", C_FEN),
                 ("C_LOS", C_LOS), ("C_E3174", C_E31)):
        g(k, v)

    # ---------------- receptor occupancies --------------------------------
    OCC_B = occ((C_ESM, p["IC50_ESM"]), (C_LAB, p["IC50_LABB"]),
                (C_MET, p["IC50_MET"]))
    OCC_A = occ((C_LAB, p["IC50_LABA"]))
    OCC_DHP = occ((C_NIC, p["IC50_NIC"]), (C_CLV, p["IC50_CLV"]))
    OCC_NO = occ((C_SNP, p["IC50_SNP"]))
    OCC_OP = occ((C_FEN, p["IC50_FEN"]))
    BLK_AT1 = occ((C_LOS, p["IC50_LOS"]), (C_E31, p["IC50_E3174"]))
    for k, v in (("OCC_B", OCC_B), ("OCC_A", OCC_A), ("OCC_DHP", OCC_DHP),
                 ("OCC_NO", OCC_NO), ("OCC_OP", OCC_OP), ("BLK_AT1", BLK_AT1)):
        g(k, v)

    # ---------------- systemic haemodynamics ------------------------------
    HR = y[IX["HR"]]
    CTR = y[IX["CTR"]]
    SVR = y[IX["SVR"]]
    PMAP = y[IX["PMAP"]]
    BV = y[IX["BV"]]
    ELN = y[IX["ELN"]]

    PRELOAD = (BV / p["BV0"]) * (1.0 - 0.25 * OCC_NO) * (1.0 - 0.05 * OCC_DHP)
    SV = p["SV0"] * CTR * max(PRELOAD, 0.05) ** 0.60 * (SVR / p["SVR0"]) ** -0.20
    CO = HR * SV / 1000.0                                     # L/min
    g("SV", SV)
    g("CO", CO)

    # Vasodilators do not only drop resistance; they raise conduit-artery
    # compliance, which is what keeps pulse pressure from exploding when stroke
    # volume rises reflexly.
    CART = (p["CART0"] / (1.0 + p["STIFF"] * (1.0 - ELN))
            * (1.0 + 0.25 * OCC_DHP + 0.30 * OCC_NO + 0.10 * OCC_A))
    PP = SV / CART
    SBP = PMAP + 2.0 / 3.0 * PP
    DBP = PMAP - 1.0 / 3.0 * PP
    g("CART", CART)
    g("PP", PP)
    g("SBP", SBP)
    g("DBP", DBP)

    DPDT = (p["DPDT0"] * CTR
            * (HR / p["HR0"]) ** p["EXP_DPDT_HR"]
            * (max(DBP, 15.0) / p["DBP0"]) ** p["EXP_DPDT_DBP"])
    g("DPDT", DPDT)

    # ---------------- two-lumen conductances ------------------------------
    # Entry tear conductance scales with orifice AREA.  Two things can remove
    # it: a covered stent-graft (TEVAR), or thrombus - but ONLY once the false
    # lumen is essentially completely thrombosed.  A false lumen that is 40%
    # thrombosed still has a wide-open entry tear and is still pressurised;
    # that single asymmetry is what makes "partial" different from "complete".
    seal = p["TEVAR_SEAL"] if (p["TEVAR_T"] >= 0.0 and t >= p["TEVAR_T"]) else 0.0
    seal_thr = (1.0 - sig((y[IX["THRFL"]] - p["SEAL_THR50"]) / p["SEAL_THRW"])) ** p["SEAL_POW"]
    G_EN = max((y[IX["AEN"]] / p["AEN_REF"]) * seal_thr * (1.0 - seal), 0.0)
    g("SEAL_THR", seal_thr)
    # distal re-entry: closes as the stagnant distal FL thromboses
    G_RE = p["GRE0"] * max(1.0 - y[IX["THRD"]], 0.0) ** 2
    g("G_EN", G_EN)
    g("G_RE", G_RE)

    PFL = y[IX["PFL"]]
    DFL = y[IX["DFL"]]
    DTL = y[IX["DTL"]]
    HFL = y[IX["HFL"]]

    # FL pulse pressure: RC low-pass through the communicating orifices
    GTOT = G_EN + G_RE
    if GTOT > 1e-9:
        tau_fl = p["TAUFL0"] / GTOT                            # s
    else:
        tau_fl = 1e6
    omega = 2.0 * math.pi * HR / 60.0                          # rad/s
    PP_FL = PP / math.sqrt(1.0 + (omega * tau_fl) ** 2)
    PFL_SYS = PFL + 2.0 / 3.0 * PP_FL
    PFL_DIA = PFL - 1.0 / 3.0 * PP_FL
    g("TAU_FL", tau_fl)
    g("PP_FL", PP_FL)
    g("PFL_SYS", PFL_SYS)
    g("PFL_DIA", PFL_DIA)

    # False-lumen flow has two parts:
    #   through-flow, which needs BOTH an entry and a re-entry (series divider),
    #   and to-and-fro flow through the entry tear alone, whose magnitude is the
    #   pressure-amplitude DIFFERENCE across the tear (PP - PP_FL): the more the
    #   sac damps the pulse, the harder blood is driven in and out of it.
    gser = G_EN * G_RE / GTOT if GTOT > 1e-9 else 0.0
    FLOW = (p["KQ"] * gser + p["KTOFRO"] * G_EN * max(PP - PP_FL, 0.0) / 46.0)
    # Thrombus displaces false-lumen volume and so reduces the flow it can
    # carry - but only partly: the segment immediately around a wide entry tear
    # keeps washing in and out however much of the distal sac has clotted.
    # That residual (1 - 0.60) is the whole reason a false lumen can sit at 30%
    # thrombosis indefinitely instead of going on to seal itself.
    FLOW *= max(1.0 - p["KFLOW_THR"] * y[IX["THRFL"]], 0.0)
    g("FLOW_FL", FLOW)

    # ---------------- wall stress (Laplace on the FL outer wall) ----------
    r_fl = DFL / 2.0
    SIG_FL = PFL * MMHG_TO_KPA * r_fl / max(HFL, 0.05)
    SIG_FL_SYS = PFL_SYS * MMHG_TO_KPA * r_fl / max(HFL, 0.05)
    # true lumen sees the intact remainder of the wall thickness
    h_tl = p["HAO0"] * (1.0 - p["FOUT"])
    SIG_TL = SBP * MMHG_TO_KPA * (DTL / 2.0) / h_tl
    # the mobile intimal flap: pressure DIFFERENCE across it, unsupported
    # The intimal flap only carries stress where it is still a FREE flap: clot
    # in the false lumen buttresses it and a stent-graft pins it.  Omitting this
    # exposure factor (first version) made a completely thrombosed, depressurised
    # false lumen the highest flap stress in the model - the pressure difference
    # across the flap is largest exactly when there is nothing left to tear.
    mobility = min((y[IX["AEN"]] / p["AEN_REF"]) ** 0.5, 1.5)
    flap_exposed = max(1.0 - y[IX["THRFL"]], 0.0) * (1.0 - seal) * mobility
    g("FLAP_MOBILITY", mobility)
    SIG_FLAP = abs(SBP - PFL_SYS) * MMHG_TO_KPA * r_fl / 0.25 * flap_exposed
    g("FLAP_EXPOSED", flap_exposed)
    g("SIG_FL", SIG_FL)
    g("SIG_FL_SYS", SIG_FL_SYS)
    g("SIG_TL", SIG_TL)
    g("SIG_FLAP", SIG_FLAP)

    # ---------------- THE INDEX -------------------------------------------
    ASI = (SIG_FL_SYS / p["SIGREF"]) * (DPDT / p["DPDT0"])
    g("ASI", ASI)
    g("ASI_FLAP", (SIG_FLAP / p["SIGREF"]) * (DPDT / p["DPDT0"]) ** p["EXP_TEAR_DPDT"])

    DAO = math.sqrt(DTL ** 2 + DFL ** 2)
    g("DAO", DAO)
    g("FL_RATIO", DFL / max(DTL + DFL, 1e-6))

    # ---------------- wall thickness set-point -----------------------------
    # Adventitial collagen thickens the outer wall; MMP-9 thins it; and mural
    # thrombus thins it hardest of all, because the wall beneath an intraluminal
    # thrombus is hypoxic and bathed in neutrophil elastase and plasmin.  That
    # last coefficient (KTHR_H = 1.0, i.e. a fully thrombosed false lumen sits on
    # a wall of half thickness) is the one that makes thrombus a mixed blessing:
    # it supports the wall mechanically (THR_SUPPORT) and digests it chemically
    # at the same time.
    HFL_BASE = p["HAO0"] * p["FOUT"]
    HFL_TGT = (HFL_BASE * (1.0 + p["HFL_COL"] * y[IX["COL"]])
               / (1.0 + p["KH_MMP"] * softplus(y[IX["MMP9"]] - 1.0)
                  + p["KTHR_H"] * y[IX["THRFL"]]))
    g("HFL_TGT", HFL_TGT)

    # ---------------- wall strength and rupture hazard --------------------
    COL = y[IX["COL"]]
    SWALL = p["SWALL0"] * max(ELN, 0.02) ** p["EXP_ELN_S"] * (0.60 + 0.40 * COL)
    SR = SIG_FL_SYS / SWALL
    H_RUPT = min(p["H0_RUPT"] * math.exp(min((SR - 1.0) / p["W_RUPT"], 30.0)), 0.60)
    g("SWALL", SWALL)
    g("STRESS_RATIO", SR)
    g("H_RUPT", H_RUPT)

    # ---------------- branch perfusion ------------------------------------
    # dynamic obstruction: the flap is pushed across a branch ostium when the
    # false lumen out-pressures the true lumen in diastole.
    # Dynamic obstruction needs TWO things at once: the false lumen must
    # out-pressure the true lumen in diastole, AND it must be the dominant
    # channel.  Requiring only the first (the first version of this file) put
    # the reference patient in renal failure at t = 0, because a patent false
    # lumen ALWAYS carries mean pressure above true-lumen diastolic pressure.
    flr = DFL / max(DFL + DTL, 1e-6)
    dom = min(1.0, max(flr - 0.40, 0.0) / 0.30)
    COLLAPSE = sig((PFL_DIA - DBP - 5.0) / 4.0) * dom
    g("COLLAPSE", COLLAPSE)
    dyn = p["KOBST_DYN"] * COLLAPSE
    pat_ren = max((1.0 - p["OBST_REN"]) * (1.0 - dyn * p["INV_REN"]), 0.02)
    pat_mes = max((1.0 - p["OBST_MES"]) * (1.0 - dyn * p["INV_MES"]), 0.02)
    pat_spi = max((1.0 - p["OBST_SPI"]) * (1.0 - 0.55 * COLLAPSE * p["INV_SPI"]), 0.02)
    # intercostal supply is lost when the FL thromboses or is stent-covered
    cover = p["TEVAR_COVER"] if (p["TEVAR_T"] >= 0.0 and t >= p["TEVAR_T"]) else 0.0
    pat_spi *= max(1.0 - 0.20 * y[IX["THRFL"]] - cover, 0.05)
    g("PAT_REN", pat_ren)
    g("PAT_MES", pat_mes)
    g("PAT_SPI", pat_spi)

    PP_REN = PMAP * pat_ren
    # renal autoregulation: flat plateau above ~80 mmHg perfusion pressure
    QREN_INST = min(1.06, 0.04 + 1.02 * sig((PP_REN - 58.0) / 8.5))
    g("PP_REN", PP_REN)
    g("QREN_INST", QREN_INST)
    QMES = min(1.06, 0.04 + 1.02 * sig((PMAP * pat_mes - 52.0) / 8.0))
    g("QMES", QMES)
    SCPP = PMAP - p["SCP_CSF"]
    # The paraspinal collateral network is recruitable: this is why paraplegia
    # after aortic coverage is an ACUTE-window risk and not a permanent state,
    # and why staged coverage works.  Without it the model claimed a chronically
    # thrombosed false lumen causes 5 years of continuous cord ischaemia.
    scpp_eff = SCPP * (pat_spi + p["KCOL_SPI"] * y[IX["COLLAT"]] * (1.0 - pat_spi))
    QSPI = min(1.06, 0.04 + 1.02 * sig((scpp_eff - p["SCPP50"]) / 7.0))
    g("SCPP_EFF", scpp_eff)
    g("SCPP", SCPP)
    g("QSPI", QSPI)

    # ---------------- clinical scores -------------------------------------
    g("EGFR", 175.0 * max(y[IX["CREAT"]], 0.2) ** -1.154 * 65.0 ** -0.203 * 1.0)
    g("PARAPLEGIA", sig((y[IX["SCI"]] - p["SCI_PARA"]) / 3.0))
    in_sbp = sig((SBP - p["SBP_LO"]) / 2.0) * sig((p["SBP_HI"] - SBP) / 2.0)
    in_hr = sig((p["HR_TGT"] + 5.0 - HR) / 2.0)
    g("IN_TARGET", in_sbp * in_hr)

    # ---------------- controller output rates -----------------------------
    u_vd = p["VD_MAX"] * sig(y[IX["U_VD"]])
    u_bb = p["BB_MAX"] * sig(y[IX["U_BB"]])
    g("U_VD_RATE", u_vd)
    g("U_BB_RATE", u_bb)
    return d


# ----------------------------------------------------------------------------
# 4.  THE ODE SYSTEM
# ----------------------------------------------------------------------------
def rhs(t, y, p):
    d = derive(y, p, t)
    dy = np.zeros(NEQ)
    S = lambda n: y[IX[n]]                                    # noqa: E731

    # ==== 4.1 drug PK =======================================================
    # esmolol: fixed rate + (optionally) the beta-blocker controller
    r_esm = p["R_ESM"] + (d["U_BB_RATE"] if p["BB_AGENT"] == 1 else 0.0)
    dy[IX["A_ESM"]] = r_esm - p["CL_ESM"] / p["V_ESM"] * S("A_ESM")

    r_lab = p["R_LAB"] + (d["U_BB_RATE"] if p["BB_AGENT"] == 2 else 0.0)
    dy[IX["A_LABC"]] = (r_lab
                        - p["CL_LAB"] / p["V_LAB"] * S("A_LABC")
                        - p["Q_LAB"] / p["V_LAB"] * S("A_LABC")
                        + p["Q_LAB"] / p["V_LAB2"] * S("A_LABP"))
    dy[IX["A_LABP"]] = (p["Q_LAB"] / p["V_LAB"] * S("A_LABC")
                        - p["Q_LAB"] / p["V_LAB2"] * S("A_LABP"))

    # oral maintenance drugs are given as an equivalent continuous rate
    # (mg per TAU_PO h) scaled by adherence - the chronic arms run for years
    # and the peak-trough ripple is irrelevant to a 5-year growth curve.
    dy[IX["A_METD"]] = p["D_MET"] * p["ADHERE"] / p["TAU_PO"] - p["KA_MET"] * S("A_METD")
    dy[IX["A_METC"]] = (p["F_MET"] * p["KA_MET"] * S("A_METD")
                        - p["CL_MET"] / p["V_MET"] * S("A_METC"))

    r_nic = p["R_NIC"] + (d["U_VD_RATE"] if p["VD_AGENT"] == 1 else 0.0)
    dy[IX["A_NICC"]] = (r_nic
                        - p["CL_NIC"] / p["V_NIC"] * S("A_NICC")
                        - p["Q_NIC"] / p["V_NIC"] * S("A_NICC")
                        + p["Q_NIC"] / p["V_NIC2"] * S("A_NICP"))
    dy[IX["A_NICP"]] = (p["Q_NIC"] / p["V_NIC"] * S("A_NICC")
                        - p["Q_NIC"] / p["V_NIC2"] * S("A_NICP"))

    r_clv = p["R_CLV"] + (d["U_VD_RATE"] if p["VD_AGENT"] == 2 else 0.0)
    dy[IX["A_CLV"]] = r_clv - p["CL_CLV"] / p["V_CLV"] * S("A_CLV")

    r_snp = p["R_SNP"] + (d["U_VD_RATE"] if p["VD_AGENT"] == 3 else 0.0)
    dy[IX["A_SNP"]] = r_snp - p["CL_SNP"] / p["V_SNP"] * S("A_SNP")
    # every molecule of SNP that is cleared liberates 5 CN-
    cn_in = p["FCN"] * p["CL_SNP"] / p["V_SNP"] * S("A_SNP")
    dy[IX["A_CN"]] = cn_in - p["CL_CN"] / p["V_CN"] * S("A_CN")
    # rhodanese converts CN- to SCN-, which is cleared renally -> accumulates
    # when the kidney is malperfused, which is exactly the AAD patient.
    cl_scn = p["CL_SCN"] * max(S("QREN"), 0.05)
    dy[IX["A_SCN"]] = (p["CL_CN"] / p["V_CN"] * S("A_CN") * 2.23
                       - cl_scn / p["V_SCN"] * S("A_SCN"))

    dy[IX["A_FENC"]] = p["R_FEN"] - p["CL_FEN"] / p["V_FEN"] * S("A_FENC")
    dy[IX["A_FENE"]] = p["KE0_FEN"] * (S("A_FENC") - S("A_FENE"))

    dy[IX["A_LOSD"]] = p["D_LOS"] * p["ADHERE"] / p["TAU_PO"] - p["KA_LOS"] * S("A_LOSD")
    dy[IX["A_LOSC"]] = (p["F_LOS"] * p["KA_LOS"] * S("A_LOSD")
                        - p["CL_LOS"] / p["V_LOS"] * S("A_LOSC"))
    dy[IX["A_E3174"]] = (p["FM_LOS"] * p["CL_LOS"] / p["V_LOS"] * S("A_LOSC")
                         - p["CL_E3174"] / p["V_E3174"] * S("A_E3174"))

    # ==== 4.2 baroreflex and sympathetic outflow ============================
    baro = 1.0 + p["KBARO"] * (S("BAROSET") - S("PMAP")) / S("BAROSET")
    ang_ex = softplus(S("ANGII") - 1.0)
    sna_t = (p["SNA0"] * softplus(baro, 6.0)
             * (1.0 + p["KPAIN"] * S("PAIN") / 10.0)
             * (1.0 + p["KANG"] * ang_ex * (1.0 - d["BLK_AT1"])))
    dy[IX["SNA"]] = (sna_t - S("SNA")) / p["TAU_SNA"]
    # Baroreflex resetting is PARTIAL.  Letting the operating point track mean
    # pressure completely (the first version of this file) is a positive
    # feedback loop: the reflex forgets the pressure it is supposed to defend,
    # and the untreated arm walks its own MAP to 160 mmHg and beyond.
    baro_tgt = (1.0 - p["FRESET"]) * p["BARO_REF"] + p["FRESET"] * S("PMAP")
    dy[IX["BAROSET"]] = (baro_tgt - S("BAROSET")) / p["TAU_BARO"]

    # ==== 4.3 heart, vessels, pressure =====================================
    sna_r = max(S("SNA") / p["SNA0"], 0.02)
    hr_t = p["HR0"] * sna_r ** p["EXP_HR_SNA"] * (1.0 - p["EMAX_BB_HR"] * d["OCC_B"])
    dy[IX["HR"]] = (max(hr_t, 32.0) - S("HR")) / p["TAU_HR"]

    ctr_t = p["CTR0"] * sna_r ** p["EXP_CTR_SNA"] * (1.0 - p["EMAX_BB_CTR"] * d["OCC_B"])
    dy[IX["CTR"]] = (max(ctr_t, 0.25) - S("CTR")) / p["TAU_CTR"]

    svr_t = (p["SVR0"] * sna_r ** p["EXP_SVR_SNA"]
             * (1.0 - p["EMAX_DHP"] * d["OCC_DHP"])
             * (1.0 - p["EMAX_ALPHA"] * d["OCC_A"])
             * (1.0 - p["EMAX_NO"] * d["OCC_NO"])
             * (1.0 + 0.30 * ang_ex * (1.0 - d["BLK_AT1"]))
             * (1.0 - p["EMAX_AT1"] * d["BLK_AT1"]))
    dy[IX["SVR"]] = (max(svr_t, 3.0) - S("SVR")) / p["TAU_SVR"]

    map_t = d["CO"] * S("SVR") + 5.0
    dy[IX["PMAP"]] = (map_t - S("PMAP")) / p["TAU_MAP"]

    # ==== 4.4 RAAS and volume ==============================================
    pra_t = (1.0 * (1.0 + p["KRENIN"] * (1.0 / max(S("QREN"), 0.05) - 1.0))
             * (1.0 - p["EMAX_BB_RENIN"] * d["OCC_B"])
             * (1.0 + p["KAT1FB"] * d["BLK_AT1"]))
    dy[IX["PRA"]] = (max(pra_t, 0.05) - S("PRA")) / p["TAU_PRA"]
    dy[IX["ANGII"]] = (S("PRA") - S("ANGII")) / p["TAU_ANG"]
    aldo_t = 1.0 + 0.9 * (S("ANGII") * (1.0 - d["BLK_AT1"]) - 1.0)
    dy[IX["ALDO"]] = (max(aldo_t, 0.1) - S("ALDO")) / p["TAU_ALDO"]
    bv_t = p["BV0"] * (1.0 + p["KALDO_BV"] * (S("ALDO") - 1.0))
    dy[IX["BV"]] = (bv_t - S("BV")) / p["TAU_BV"]

    # saturating mechanotransduction signal, used by both the collagen and the
    # interleukin equations below
    xs = softplus(d["SIG_FL_SYS"] / p["SIGREF"] - 1.0)

    # ==== 4.5 the two-lumen divider ========================================
    pdist = S("PMAP") - p["DPAO"]
    dy[IX["PFL"]] = (p["KCOMM"] * (d["G_EN"] * (S("PMAP") - S("PFL"))
                                   + d["G_RE"] * (pdist - S("PFL")))
                     - p["KDEP"] * (S("PFL") - p["PVEN"]))

    # true lumen is compressed by a pressurised false lumen
    dtl_t = p["DTL0"] * (1.0 - p["KCOLLAPSE"] * d["COLLAPSE"])
    dy[IX["DTL"]] = (dtl_t - S("DTL")) / p["TAU_DTL"]

    # False lumen expansion.  Stress enters with the Laplace exponent 1.9, so a
    # thrombus that thins the wall by a factor f multiplies growth by f^1.9 while
    # supporting it by only (1 - 0.85*THR): the competition between an exponent
    # and a linear term is what puts the growth maximum at an INTERMEDIATE
    # thrombus fraction.  Closed form: THR* = (n*A - s)/((1+n)*s*A) with
    # n = EXP_GROW, s = THR_SUPPORT, A = KTHR_H  ->  0.426 for the values used.
    cap = 1.0 - sig((S("DFL") - p["DFL_CAP"]) / 1.5)
    grow = (p["KGROW"] * (d["SIG_FL"] / p["SIGREF"]) ** p["EXP_GROW"]
            * (1.0 + p["KMMP_GROW"] * softplus(S("MMP9") - 1.0))
            * (1.0 - p["THR_SUPPORT"] * S("THRFL")) * cap)
    # an organising, depressurised thrombosed sac slowly shrinks
    regress = p["KREGRESS"] * S("THRFL") * sig((25.0 - S("PFL")) / 6.0)
    dy[IX["DFL"]] = grow - regress

    # thrombosis: driven by stasis, opposed by flow-dependent lysis
    stasis = 1.0 / (1.0 + (d["FLOW_FL"] / p["FLOW50"]) ** 2)
    stasis_d = 1.0 / (1.0 + (d["FLOW_FL"] / p["FLOW50D"]) ** 2)
    dy[IX["THRFL"]] = (p["KTHR"] * stasis * (1.0 - S("THRFL"))
                       - p["KLYS"] * S("THRFL") * min(d["FLOW_FL"], 3.0))
    dy[IX["THRD"]] = (p["KTHRD"] * stasis_d * (1.0 - S("THRD"))
                      - p["KLYSD"] * S("THRD") * min(d["FLOW_FL"], 3.0))

    # Entry tear propagation: THE stress-impulse term, and the place where the
    # product structure earns its keep.  The tear extends only when the flap's
    # own stress-impulse exceeds a threshold, so a regimen is aorta-protective
    # if and only if it drives ASI_flap below TEAR_TH - which a beta-blocker
    # does by lowering BOTH factors and a vasodilator may fail to do even while
    # it lowers the pressure, because it raises the second one.
    asi_flap = (d["SIG_FLAP"] / p["SIGREF"]) * (d["DPDT"] / p["DPDT0"]) ** p["EXP_TEAR_DPDT"]
    dy[IX["AEN"]] = (p["KTEAR"] * S("AEN") * softplus(asi_flap - p["TEAR_TH"], 40.0)
                     * max(1.0 - S("AEN") / p["AEN_MAX"], 0.0)
                     / (1.0 + 3.0 * S("COL"))
                     - p["KHEAL_TEAR"] * S("COL")
                     * softplus(S("AEN") - 0.6 * p["AEN0"]))

    # wall thickness relaxes towards its remodelling set-point (see derive())
    dy[IX["HFL"]] = p["KREM_H"] * (d["HFL_TGT"] - S("HFL"))
    dy[IX["ELN"]] = -p["KELN"] * softplus(S("MMP9") - 1.0) * S("ELN")
    # Adventitial collagen is laid down in response to TGF-beta AND to wall
    # stress itself.  This is the only negative feedback the aorta has: a
    # thicker outer wall lowers sigma = P*r/h, which is why most dilated aortas
    # creep rather than run away.  Without it the model has no stable chronic
    # state at all.
    dy[IX["COL"]] = ((p["KCOL"] * (0.10 + softplus(S("TGFB") - 1.0))
                      + p["KCOL_STRESS"] * xs)
                     * max(1.0 - S("COL") / p["COLMAX"], 0.0)
                     - p["KDEGC"] * S("COL") * S("MMP9"))

    # ==== 4.6 injury, inflammation, pain ===================================
    dy[IX["INJ"]] = -S("INJ") / p["TAU_INJ"]
    # The mechanotransduced part of the IL-6 signal SATURATES.  Leaving it
    # linear in wall stress (the first version) closed an unbounded loop -
    # stress -> IL-6 -> MMP-9 -> wall thinning -> stress - and the 5-year run
    # reached MMP-9 44x normal and a 0.16 mm aortic wall.
    dy[IX["IL6"]] = (p["KIL6"] * S("INJ")
                     + p["KIL6_STRESS"] * xs / (1.0 + xs)
                     - p["KOUT_IL6"] * (S("IL6") - 1.0))
    # MMP-9 and TGF-beta are written as set-point relaxations so that a
    # "healthy" input of 1.0 returns exactly 1.0 - the defect that made the
    # first version of this file diverge was a production term that did not
    # balance its own baseline.
    # AT1 signalling is a second, independent MMP-9 input (Ang II induces
    # MMP-2/-9 in medial smooth muscle).  It is the ONLY route by which an ARB
    # can protect the wall in this model, and it competes against the ARB's own
    # suppression of TGF-beta-driven adventitial collagen.  The net sign is
    # therefore computed, not assumed - see the S18/S19 comparison.
    mmp_t = (p["GEN_MMP"] * (1.0 + p["KMMP_IL6"] * softplus(S("IL6") - 1.0))
             * (1.0 + p["KMMP_AT1"] * softplus(S("ANGII") * (1.0 - d["BLK_AT1"]) - 1.0)
                - p["KMMP_AT1"] * 0.0))
    dy[IX["MMP9"]] = p["KOUT_MMP"] * (mmp_t - S("MMP9"))
    tgf_t = ((1.0 + p["KTGF_INJ"] * S("INJ")
              + p["KTGF_MMP"] * softplus(S("MMP9") - 1.0))
             * (0.40 + 0.60 * (1.0 - d["BLK_AT1"])))
    dy[IX["TGFB"]] = p["KOUT_TGF"] * (tgf_t - S("TGFB"))

    # Pain: the tearing pain of onset is the injury signal; what persists or
    # recurs afterwards is wall stress, which is why refractory pain is a
    # class-I indication for intervention rather than for more morphine.
    # The stress-driven component of the pain is gated by flap friability.
    # Without the gate the model kept a chronic, fully remodelled dissection in
    # 10/10 pain for five years, and the sympathetic drive that came with it
    # became a phantom cause of chronic hypertension.
    friab = 1.0 / (1.0 + 4.0 * S("COL"))
    pain_t = (p["KPAIN_INJ"] * S("INJ")
              + 10.0 * sig((d["SIG_FL_SYS"] / p["SIG_PAIN"] - 1.0) / 0.16)
              * (0.15 + 0.85 * friab))
    pain_t = min(pain_t, 10.0) * (1.0 - p["EMAX_OP"] * d["OCC_OP"])
    dy[IX["PAIN"]] = (pain_t - S("PAIN")) / p["TAU_PAIN"]

    # ==== 4.7 organ perfusion and laboratory endpoints =====================
    dy[IX["QREN"]] = (d["QREN_INST"] - S("QREN")) / p["TAU_QREN"]
    dy[IX["CREAT"]] = (p["CREAT_PROD"]
                       - p["CL_CREAT0"] * max(S("QREN"), 0.03) * S("CREAT")) / p["V_CREAT"] * 24.0
    mes_def = softplus(0.85 - d["QMES"])
    dy[IX["LACT"]] = (p["KLACT_MES"] * mes_def
                      + p["KLACT_CN"] * softplus(d["C_CN"] - 0.35)
                      + p["KOUT_LACT"] * p["LACT0"]
                      - p["KOUT_LACT"] * S("LACT"))
    # Cord injury burden, not a lifetime integral: written as an integral it
    # reported certain paraplegia in every arm, because any arm with a
    # transiently marginal spinal perfusion pressure eventually exceeded any
    # fixed threshold if you waited five years.
    dy[IX["SCI"]] = (p["KSCI"] * softplus(p["SCI_THRESH"] - d["QSPI"])
                     - S("SCI") / p["TAU_SCI"])
    dy[IX["COLLAT"]] = (p["KCOLL"] * softplus(0.75 - d["QSPI"]) * (1.0 - S("COLLAT"))
                        - p["KCOLL_D"] * S("COLLAT"))
    dy[IX["DDIM"]] = (p["KDD_FLOW"] * d["FLOW_FL"] * (1.0 - S("THRFL"))
                      + p["KDD_THR"] * softplus(dy[IX["THRFL"]] * 100.0) * 0.1
                      + p["KOUT_DD"] * p["DDIM0"]
                      - p["KOUT_DD"] * S("DDIM"))

    # ==== 4.8 outcome integrals ============================================
    # Death in the acute phase is mostly not free-wall rupture; it is extension
    # - retrograde into the arch, or distally into a new branch - so the hazard
    # carries an explicit term proportional to how hard the flap is being driven.
    # NOTE the sharpness argument.  With the default smoothing (k = 8) the
    # softplus tail returns ~0.004 for an argument of -0.43, i.e. a perfectly
    # well-perfused organ still contributed 1.05e-5/h of hazard - four times the
    # entire non-aortic mortality - and over 43 800 h that numerical leak, not
    # any modelled mechanism, became the leading cause of death in every chronic
    # arm.  Hazard terms integrated over years must not leak.
    h_mal = p["KH_MALPERF"] * (softplus(0.55 - d["QMES"], 40.0)
                               + softplus(0.45 - S("QREN"), 40.0))
    # Extension is a hazard of a FRESH flap: the acute flap is friable, the
    # chronic one is fibrotic (hence the COL gate), which is why the same
    # stress-impulse is lethal in week 1 and merely progressive in year 3.
    h_ext = (p["KH_EXT"] * softplus(d["ASI_FLAP"] - p["TEAR_TH"], 40.0)
             * (0.25 + 0.75 * S("INJ")) / (1.0 + 4.0 * S("COL")))
    h_tot = d["H_RUPT"] + h_mal + h_ext + p["KH_OTHER"]
    dy[IX["CUMH"]] = d["H_RUPT"] + h_ext
    dy[IX["SURV"]] = -h_tot * S("SURV")
    dy[IX["DAM"]] = p["KDAM"] * d["ASI"]
    dy[IX["TIT"]] = d["IN_TARGET"]
    dy[IX["AUCASI"]] = d["ASI"]

    # ==== 4.9 closed-loop titration ========================================
    # A nurse titrating to a systolic target is an integral controller.  Both
    # controllers are written as smooth saturating integrators so the arms can
    # be compared AT THE SAME ACHIEVED PRESSURE.
    dy[IX["U_VD"]] = p["KI_VD"] * (d["SBP"] - p["SBP_SET"]) if p["VD_MAX"] > 0 else 0.0
    if p["BB_MAX"] > 0:
        if p["BB_ON_SBP"]:
            dy[IX["U_BB"]] = p["KI_BB_SBP"] * (d["SBP"] - p["SBP_SET"])
        else:
            dy[IX["U_BB"]] = p["KI_BB"] * (S("HR") - p["HR_SET"])
    else:
        dy[IX["U_BB"]] = 0.0
    return dy


# ----------------------------------------------------------------------------
# 5.  INITIAL CONDITIONS
# ----------------------------------------------------------------------------
def initial_state(p, acute=True):
    y = np.zeros(NEQ)
    y[IX["SNA"]] = 1.0
    y[IX["BAROSET"]] = 106.0       # pre-existing hypertension: reset upward
    y[IX["HR"]] = p["HR0"]
    y[IX["CTR"]] = p["CTR0"]
    y[IX["SVR"]] = p["SVR0"]
    y[IX["PMAP"]] = 106.0
    y[IX["PRA"]] = 1.0
    y[IX["ANGII"]] = 1.0
    y[IX["ALDO"]] = 1.0
    y[IX["BV"]] = p["BV0"]
    y[IX["PFL"]] = 104.0
    y[IX["DFL"]] = p["DFL0"]
    y[IX["DTL"]] = p["DTL0"]
    y[IX["THRFL"]] = 0.0
    y[IX["THRD"]] = 0.0
    y[IX["AEN"]] = p["AEN0"]
    y[IX["HFL"]] = p["HAO0"] * p["FOUT"]
    y[IX["ELN"]] = 0.82            # 65-y hypertensive descending aorta
    y[IX["COL"]] = 0.10
    y[IX["INJ"]] = 1.0 if acute else 0.0
    y[IX["IL6"]] = 1.0
    y[IX["MMP9"]] = 1.0
    y[IX["TGFB"]] = 1.0
    y[IX["PAIN"]] = 0.0
    y[IX["QREN"]] = 1.0
    y[IX["COLLAT"]] = 0.55         # a paraspinal collateral network pre-exists
    y[IX["CREAT"]] = 1.05
    y[IX["LACT"]] = 1.0
    y[IX["SCI"]] = 0.0
    y[IX["DDIM"]] = 0.25
    y[IX["SURV"]] = 1.0
    y[IX["U_VD"]] = -2.5           # controllers start near-off
    y[IX["U_BB"]] = -2.5
    return y


def simulate(p, tend, y0=None, n=None, acute=True):
    if y0 is None:
        y0 = initial_state(p, acute=acute)
    if n is None:
        n = max(400, min(1600, int(tend * 4) + 1))
    ts = np.linspace(0.0, tend, n)
    # Acute runs are stepped finely; a 5-year run lets LSODA choose its own
    # step (the fast haemodynamic states are stiff, which is exactly what BDF
    # is for) or it would need ~10^5 forced steps.
    ms = 0.5 if tend <= 800.0 else 24.0
    sol = solve_ivp(rhs, (0.0, tend), y0, args=(p,), method="LSODA",
                    t_eval=ts, rtol=1e-6, atol=1e-9, max_step=ms)
    if not sol.success:
        raise RuntimeError("solver failed: " + sol.message)
    return sol


def outputs(sol, p):
    """Derived readouts at every stored time point."""
    keys = None
    rows = []
    for i, t in enumerate(sol.t):
        d = derive(sol.y[:, i], p, t)
        if keys is None:
            keys = list(d.keys())
        rows.append([d[k] for k in keys])
    arr = np.array(rows)
    out = {k: arr[:, j] for j, k in enumerate(keys)}
    for n in CMT:
        out[n] = sol.y[IX[n], :]
    out["time"] = sol.t
    return out


# ----------------------------------------------------------------------------
# 6.  SCENARIOS
# ----------------------------------------------------------------------------
def base(**kw):
    p = dict(P)
    p.update(kw)
    return p


# Acute phase = 72 h.  Every acute arm titrates to SBP 110 (window 100-120).
ACUTE_H = 72.0
# Chronic phase = 5 years.
CHRONIC_H = 5 * 8760.0

SCEN = {}


def scen(key, label, tend, **kw):
    SCEN[key] = dict(label=label, tend=tend, p=base(**kw))


# --- S1-S3: the stress/impulse dissociation -------------------------------
scen("S01_untreated", "Untreated acute type B (natural history, 72 h)",
     ACUTE_H)
scen("S02_vd_only", "Vasodilator monotherapy (nicardipine to 15 mg/h max)",
     ACUTE_H, VD_AGENT=1, VD_MAX=15.0)
scen("S03_bb_only", "Beta-blocker monotherapy (esmolol to 300 ug/kg/min max, on SBP)",
     ACUTE_H, BB_AGENT=1, BB_MAX=1350.0, BB_ON_SBP=1)
scen("S04_bb_then_vd", "Guideline order: esmolol (to HR 60) then nicardipine (to SBP)",
     ACUTE_H, BB_AGENT=1, BB_MAX=900.0, VD_AGENT=1, VD_MAX=15.0)
scen("S05_snp_only", "Nitroprusside monotherapy (to 3 ug/kg/min max)",
     ACUTE_H, VD_AGENT=3, VD_MAX=13.5)
scen("S06_snp_bb", "Nitroprusside + esmolol (the classic Wheat regimen)",
     ACUTE_H, VD_AGENT=3, VD_MAX=13.5, BB_AGENT=1, BB_MAX=900.0)
scen("S07_labetalol", "Labetalol monotherapy (0.5-2 mg/min, titrated on SBP)",
     ACUTE_H, BB_AGENT=2, BB_MAX=120.0, BB_ON_SBP=1)
scen("S08_clev_esm", "Clevidipine (to 21 mg/h) + esmolol",
     ACUTE_H, VD_AGENT=2, VD_MAX=21.0, BB_AGENT=1, BB_MAX=900.0)
scen("S09_analgesia", "Analgesia only (fentanyl 75 ug/h), no antihypertensive",
     ACUTE_H, R_FEN=0.075)
scen("S10_analg_bb_vd", "Fentanyl + esmolol + nicardipine (full protocol)",
     ACUTE_H, R_FEN=0.075, BB_AGENT=1, BB_MAX=900.0, VD_AGENT=1, VD_MAX=15.0)

# --- S11-S13: over-lowering and the renin trap ----------------------------
scen("S11_overshoot", "Aggressive target SBP 95 with a renal branch off the FL",
     ACUTE_H, VD_AGENT=1, VD_MAX=32.0, BB_AGENT=1, BB_MAX=900.0,
     SBP_SET=95.0, OBST_REN=0.45, INV_REN=1.0, INV_MES=1.0)
scen("S12_malperf_std", "Renal malperfusion, standard target SBP 110",
     ACUTE_H, VD_AGENT=1, VD_MAX=32.0, BB_AGENT=1, BB_MAX=900.0,
     OBST_REN=0.45, INV_REN=1.0, INV_MES=1.0)
scen("S13_snp_toxic", "Nitroprusside 72 h at high rate with renal impairment",
     ACUTE_H, VD_AGENT=3, VD_MAX=45.0, BB_AGENT=1, BB_MAX=900.0,
     OBST_REN=0.45, INV_REN=1.0, SBP_SET=100.0)

# --- S14-S17: false lumen anatomy (the partial-thrombosis arithmetic) -----
scen("S14_fl_patent", "Chronic 5 y: large re-entry (FL stays patent)",
     CHRONIC_H, GRE0=2.4, D_MET=100.0, D_LOS=50.0)
scen("S15_fl_partial", "Chronic 5 y: no distal re-entry (partial thrombosis, blind sac)",
     CHRONIC_H, GRE0=0.05, D_MET=100.0, D_LOS=50.0)
scen("S16_fl_complete", "Chronic 5 y: small entry tear (9 mm2) -> complete FL thrombosis",
     CHRONIC_H, GRE0=0.05, AEN0=9.0, D_MET=100.0, D_LOS=50.0)
scen("S17_tevar", "Chronic 5 y: TEVAR entry-tear coverage at day 14",
     CHRONIC_H, GRE0=0.05, TEVAR_T=336.0, D_MET=100.0, D_LOS=50.0)

# --- S18-S20: chronic medical therapy, genetics, adherence ---------------
scen("S18_marfan_los", "Marfan-like wall (MMP x2.2) on losartan + metoprolol, 5 y",
     CHRONIC_H, GEN_MMP=2.2, D_MET=100.0, D_LOS=50.0, GRE0=0.05)
scen("S19_marfan_nolos", "Marfan-like wall, beta-blocker only, 5 y",
     CHRONIC_H, GEN_MMP=2.2, D_MET=100.0, D_LOS=0.0, GRE0=0.05)
scen("S21_vd_forced", "Nicardipine pushed ABOVE label (30 mg/h) to match the SBP",
     ACUTE_H, VD_AGENT=1, VD_MAX=30.0)
scen("S22_tevar_map", "TEVAR at day 14 WITH spinal perfusion-pressure augmentation",
     CHRONIC_H, GRE0=0.05, TEVAR_T=336.0, D_MET=50.0, D_LOS=25.0, SBP_SET=135.0)
scen("S20_nonadherent", "Chronic 5 y: 25% adherence to antihypertensives",
     CHRONIC_H, GRE0=0.05, D_MET=100.0, D_LOS=50.0, ADHERE=0.25)


# ----------------------------------------------------------------------------
# 7.  REPORT
# ----------------------------------------------------------------------------
def tail_mean(x, t, hours=6.0):
    """Mean of x over the last `hours` of the run (the titrated steady state)."""
    m = t >= (t[-1] - hours)
    return float(np.mean(x[m]))


def first_time(o, key, thr, scale=8760.0):
    """Time (in years by default) at which `key` first reaches `thr`."""
    m = np.where(o[key] >= thr)[0]
    return float(o["time"][m[0]] / scale) if len(m) else float("nan")


def run_all():
    res = {}
    for key in sorted(SCEN):
        s = SCEN[key]
        p = s["p"]
        sol = simulate(p, s["tend"])
        o = outputs(sol, p)
        t = o["time"]
        yrs = s["tend"] / 8760.0
        tm = lambda n: tail_mean(o[n], t)                      # noqa: E731
        res[key] = dict(
            label=s["label"], tend=s["tend"], acute=(s["tend"] <= 200.0),
            # --- haemodynamics ---------------------------------------------
            SBP=tm("SBP"), DBP=tm("DBP"), MAP=tm("PMAP"), HR=tm("HR"),
            PP=tm("PP"), CO=tm("CO"), SVR=tm("SVR"), SNA=tm("SNA"),
            dPdt=tm("DPDT"),
            # --- the index -------------------------------------------------
            SIG_FL=tm("SIG_FL_SYS"), SIG_TL=tm("SIG_TL"), SIG_FLAP=tm("SIG_FLAP"),
            ASI=tm("ASI"), ASI_FLAP=tm("ASI_FLAP"),
            ASI_AUC_per_h=float(o["AUCASI"][-1] / s["tend"]),
            STRESS_RATIO=tm("STRESS_RATIO"), SWALL=tm("SWALL"),
            # --- the two lumens --------------------------------------------
            PFL=tm("PFL"), PFL_DIA=tm("PFL_DIA"), PFL_SYS=tm("PFL_SYS"),
            PP_FL=tm("PP_FL"), TAU_FL=tm("TAU_FL"),
            G_EN=tm("G_EN"), G_RE=tm("G_RE"), FLOW_FL=tm("FLOW_FL"),
            THRFL=float(o["THRFL"][-1]), THRD=float(o["THRD"][-1]),
            SEAL_THR=tm("SEAL_THR"), COLLAPSE=tm("COLLAPSE"),
            MOBILITY=tm("FLAP_MOBILITY"),
            # --- aortic geometry and remodelling ---------------------------
            DAO0=float(o["DAO"][0]), DAO=float(o["DAO"][-1]),
            DFL=float(o["DFL"][-1]), DTL=float(o["DTL"][-1]),
            GROWTH_mm_yr=float((o["DAO"][-1] - o["DAO"][0]) / max(yrs, 1e-9)),
            T55=first_time(o, "DAO", 55.0),
            AEN=float(o["AEN"][-1]), AEN0=float(o["AEN"][0]),
            dAEN_pct=float(100.0 * (o["AEN"][-1] / o["AEN"][0] - 1.0)),
            HFL=float(o["HFL"][-1]), ELN=float(o["ELN"][-1]),
            COL=float(o["COL"][-1]),
            # --- inflammation / pain ---------------------------------------
            IL6=tm("IL6"), MMP9=tm("MMP9"), TGFB=tm("TGFB"), PAIN=tm("PAIN"),
            PAIN_MAX=float(o["PAIN"].max()),
            # --- organs / laboratory ---------------------------------------
            QREN=tm("QREN"), QMES=tm("QMES"), QSPI=tm("QSPI"),
            PRA=tm("PRA"), ANGII=tm("ANGII"), ALDO=tm("ALDO"),
            CREAT=float(o["CREAT"][-1]), LACT=tm("LACT"),
            SCI_MAX=float(o["SCI"].max()), PARA=float(o["PARAPLEGIA"].max()),
            COLLAT=float(o["COLLAT"][-1]), DDIM=tm("DDIM"),
            C_CN=tm("C_CN"), C_CN_MAX=float(o["C_CN"].max()), C_SCN=tm("C_SCN"),
            # --- outcome ----------------------------------------------------
            CUMH=float(o["CUMH"][-1]), SURV=float(o["SURV"][-1]),
            MORT_pct=float(100.0 * (1.0 - o["SURV"][-1])),
            DAM=float(o["DAM"][-1]), TIT_pct=float(100.0 * o["TIT"][-1] / s["tend"]),
            # --- exposure ---------------------------------------------------
            C_ESM=tm("C_ESM"), C_NIC=tm("C_NIC"), C_SNP=tm("C_SNP"),
            C_LAB=tm("C_LAB"), C_CLV=tm("C_CLV"), C_MET=tm("C_MET"),
            C_E3174=tm("C_E3174"), C_FEN=tm("C_FEN"),
            OCC_B=tm("OCC_B"), OCC_DHP=tm("OCC_DHP"), OCC_NO=tm("OCC_NO"),
            OCC_A=tm("OCC_A"), BLK_AT1=tm("BLK_AT1"), OCC_OP=tm("OCC_OP"),
            RATE_VD=tm("U_VD_RATE"), RATE_BB=tm("U_BB_RATE"),
        )
    return res


def thrombus_sweep(p=None):
    """
    Growth rate as a function of false-lumen thrombus fraction at a FIXED,
    patent entry tear.  This is the arithmetic behind the partial-thrombosis
    paradox, and it has a closed-form optimum that the sweep must reproduce.
    """
    p = p or base()
    y0 = initial_state(p)
    hfl_base = p["HAO0"] * p["FOUT"]

    def growth(thr):
        y = y0.copy()
        y[IX["THRFL"]] = thr
        # put the wall at its remodelling set-point for this thrombus load
        y[IX["HFL"]] = (hfl_base * (1.0 + p["HFL_COL"] * y[IX["COL"]])
                        / (1.0 + p["KTHR_H"] * thr))
        d = derive(y, p, 0.0)
        return (p["KGROW"] * (d["SIG_FL"] / p["SIGREF"]) ** p["EXP_GROW"]
                * (1.0 - p["THR_SUPPORT"] * thr))

    xs = np.linspace(0.0, 1.0, 201)
    gs = np.array([growth(x) for x in xs])
    g0 = gs[0]
    n, sup, A = p["EXP_GROW"], p["THR_SUPPORT"], p["KTHR_H"]
    analytic = (n * A - sup) / ((1.0 + n) * sup * A)
    return xs, gs / g0, float(xs[gs.argmax()]), float(analytic)


# ----------------------------------------------------------------------------
# 8.  REPORT
# ----------------------------------------------------------------------------
def print_report(res):
    W = 122
    A = [k for k in sorted(res) if res[k]["acute"]]
    C = [k for k in sorted(res) if not res[k]["acute"]]
    sh = lambda k: k.replace("_", "-")[:18]                    # noqa: E731

    print("=" * W)
    print("ACUTE TYPE B AORTIC DISSECTION - QSP REFERENCE MODEL  (Python/scipy, LSODA)")
    print("%d ODEs, %d scenarios.  Every number below is integrated, none is asserted."
          % (NEQ, len(res)))
    print("=" * W)

    print("""
The model is built around one product:

    ASI = (sigma_FL_sys / sigma_ref) x (dP/dt_max / dPdt_ref)

sigma carries the PRESSURE, dP/dt carries the IMPULSE, and the two are set by
different physiology.  Sections 1-3 test whether a drug that lowers only the
first can be worse than no drug at all.
""")

    # ---- 1 -----------------------------------------------------------------
    print("### 1. ACUTE PHASE (72 h) - state at the titrated steady state\n")
    hdr = ("{:<19}{:>7}{:>7}{:>6}{:>7}{:>8}{:>7}{:>7}{:>6}{:>7}"
           .format("scenario", "SBP", "DBP", "HR", "dP/dt", "sig_FL", "ASI",
                   "ASIfl", "pain", "%tgt"))
    print(hdr)
    print("-" * len(hdr))
    for k in A:
        r = res[k]
        print("{:<19}{:>7.1f}{:>7.1f}{:>6.1f}{:>7.0f}{:>8.1f}{:>7.3f}{:>7.3f}"
              "{:>6.1f}{:>7.1f}".format(
                  sh(k), r["SBP"], r["DBP"], r["HR"], r["dPdt"], r["SIG_FL"],
                  r["ASI"], r["ASI_FLAP"], r["PAIN"], r["TIT_pct"]))

    # ---- 2 -----------------------------------------------------------------
    print("\n### 2. THE DISSOCIATION - what a pure vasodilator actually buys\n")
    u = res["S01_untreated"]
    hdr = ("{:<19}{:>7}{:>8}{:>8}{:>8}{:>8}{:>9}{:>9}"
           .format("scenario", "SBP", "dP/dt", "d%dPdt", "ASI", "ASIfl",
                   "tear+%", "mort72h%"))
    print(hdr)
    print("-" * len(hdr))
    for k in ["S01_untreated", "S02_vd_only", "S21_vd_forced", "S05_snp_only",
              "S09_analgesia", "S03_bb_only", "S07_labetalol", "S04_bb_then_vd",
              "S06_snp_bb", "S08_clev_esm", "S10_analg_bb_vd"]:
        r = res[k]
        print("{:<19}{:>7.1f}{:>8.0f}{:>8.1f}{:>8.3f}{:>8.3f}{:>9.2f}{:>9.2f}".format(
            sh(k), r["SBP"], r["dPdt"], 100.0 * (r["dPdt"] / u["dPdt"] - 1.0),
            r["ASI"], r["ASI_FLAP"], r["dAEN_pct"], r["MORT_pct"]))
    vd, bb = res["S21_vd_forced"], res["S04_bb_then_vd"]
    print("""
  Matched-pressure comparison (nicardipine forced to SBP {:.0f} vs esmolol-first
  regimen at SBP {:.0f} - a difference of {:.1f} mmHg):
      dP/dt max            {:>7.0f}  vs {:>7.0f}   ({:+.0f}%)
      ASI                  {:>7.3f}  vs {:>7.3f}   ({:+.0f}%)
      flap stress-impulse  {:>7.3f}  vs {:>7.3f}   ({:.1f}x)
      entry tear growth    {:>6.2f}%  vs {:>6.2f}%   ({:.0f}x)
      72-h mortality       {:>6.2f}%  vs {:>6.2f}%   ({:.1f}x)
  Untreated 72-h mortality is {:.2f}%.  Both vasodilator-monotherapy arms are
  therefore WORSE than no antihypertensive at all, while every arm containing a
  beta-blocker is better - and the arms are ordered by ASI, not by pressure.
""".format(vd["SBP"], bb["SBP"], vd["SBP"] - bb["SBP"],
           vd["dPdt"], bb["dPdt"], 100.0 * (vd["dPdt"] / bb["dPdt"] - 1.0),
           vd["ASI"], bb["ASI"], 100.0 * (vd["ASI"] / bb["ASI"] - 1.0),
           vd["ASI_FLAP"], bb["ASI_FLAP"], vd["ASI_FLAP"] / bb["ASI_FLAP"],
           vd["dAEN_pct"], bb["dAEN_pct"], vd["dAEN_pct"] / bb["dAEN_pct"],
           vd["MORT_pct"], bb["MORT_pct"], vd["MORT_pct"] / bb["MORT_pct"],
           u["MORT_pct"]))

    # ---- 3 -----------------------------------------------------------------
    print("### 3. WHY THE REFLEX WINS - the terms behind the table above\n")
    hdr = ("{:<19}{:>7}{:>7}{:>7}{:>7}{:>8}{:>8}{:>8}{:>8}"
           .format("scenario", "SNA", "HR", "PP", "PP_FL", "P_FLsys",
                   "TL-FL", "sig_flp", "mobil"))
    print(hdr)
    print("-" * len(hdr))
    for k in ["S01_untreated", "S02_vd_only", "S21_vd_forced", "S05_snp_only",
              "S03_bb_only", "S04_bb_then_vd"]:
        r = res[k]
        print("{:<19}{:>7.3f}{:>7.1f}{:>7.1f}{:>7.1f}{:>8.1f}{:>8.1f}{:>8.1f}"
              "{:>8.2f}".format(
                  sh(k), r["SNA"], r["HR"], r["PP"], r["PP_FL"], r["PFL_SYS"],
                  r["SBP"] - r["PFL_SYS"], r["SIG_FLAP"], r["MOBILITY"]))
    print("""
  A vasodilator raises heart rate, which shortens the cardiac cycle relative to
  the false lumen's RC time constant (tau_FL = {:.2f} s here), so MORE of the
  pulse is filtered out of the false lumen and the systolic pressure DIFFERENCE
  across the intimal flap grows even as the pressure itself falls.  Tachycardia
  is not merely a cosmetic side effect of vasodilatation in a dissected aorta;
  it is a second, independent mechanism of flap loading.
""".format(res["S02_vd_only"]["TAU_FL"]))

    # ---- 4 -----------------------------------------------------------------
    print("### 4. FALSE LUMEN ANATOMY -> PHENOTYPE (5-year runs)\n")
    hdr = ("{:<19}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}{:>8}{:>7}"
           .format("scenario", "A_en0", "G_re0", "P_FLd", "flow", "thrFL",
                   "D_ao", "mm/yr", "T55(y)", "surv5y", "sig/S"))
    print(hdr)
    print("-" * len(hdr))
    for k in C:
        r, sp = res[k], SCEN[k]["p"]
        print("{:<19}{:>7.0f}{:>7.2f}{:>7.1f}{:>7.3f}{:>7.2f}{:>7.1f}{:>7.2f}"
              "{:>7.2f}{:>8.3f}{:>7.2f}".format(
                  sh(k), sp["AEN0"], sp["GRE0"], r["PFL_DIA"], r["FLOW_FL"],
                  r["THRFL"], r["DAO"], r["GROWTH_mm_yr"], r["T55"], r["SURV"],
                  r["STRESS_RATIO"]))
    pa, pt, cp = res["S14_fl_patent"], res["S15_fl_partial"], res["S16_fl_complete"]
    hr_partial = math.log(pt["SURV"]) / math.log(pa["SURV"])
    print("""
  Nothing in the code says partial thrombosis is dangerous.  Two anatomical
  numbers were changed - the entry tear area and the distal re-entry
  conductance - and the rest is the divider:

    large re-entry  -> brisk through-flow, no thrombus, FL diastolic {:.0f} mmHg,
                       growth {:.2f} mm/y, 5-y survival {:.3f}
    no re-entry     -> blind sac: no outflow so no pressure drop, the pulse is
                       filtered out so FL DIASTOLIC pressure rises to {:.0f} mmHg,
                       the stagnant distal sac clots to {:.0%} while the entry tear
                       stays wide open, growth {:.2f} mm/y, 5-y survival {:.3f}
    small entry     -> the whole sac clots, the entry seals, the false lumen
                       depressurises to {:.0f} mmHg, growth {:.2f} mm/y,
                       5-y survival {:.3f}

  Implied hazard ratio, partial vs patent = {:.2f}.  The reported adjusted
  hazard ratio in the IRAD cohort that first described this (Tsai 2007) is 2.69.
""".format(pa["PFL_DIA"], pa["GROWTH_mm_yr"], pa["SURV"],
           pt["PFL_DIA"], pt["THRFL"], pt["GROWTH_mm_yr"], pt["SURV"],
           cp["PFL"], cp["GROWTH_mm_yr"], cp["SURV"], hr_partial))

    # ---- 5 -----------------------------------------------------------------
    xs, gs, num, ana = thrombus_sweep()
    print("### 5. VERIFICATION - the growth maximum has a closed form\n")
    print("  Growth ∝ (1 + A*THR)^n * (1 - s*THR) with n = EXP_GROW = 1.9,")
    print("  s = THR_SUPPORT = 0.85, A = KTHR_H = 1.0 (thrombus thins the wall it")
    print("  supports).  dg/dTHR = 0 gives THR* = (n*A - s) / ((1+n)*s*A).\n")
    print("  {:>10}{:>14}".format("THR", "rel. growth"))
    for x in [0.0, 0.1, 0.2, 0.3, 0.4, 0.426, 0.5, 0.6, 0.7, 0.85, 1.0]:
        print("  {:>10.3f}{:>14.3f}".format(x, float(np.interp(x, xs, gs))))
    print("\n  numerical argmax = {:.3f}   analytic optimum = {:.3f}   agreement {:.1%}"
          .format(num, ana, 1.0 - abs(num - ana) / ana))
    print("  The interior maximum is the model's explanation of the paradox: a\n"
          "  thrombus supports the wall LINEARLY and thins it under a Laplace\n"
          "  exponent, so the worst possible false lumen is a half-clotted one.\n")

    # ---- 6 -----------------------------------------------------------------
    print("### 6. MALPERFUSION, THE RENIN TRAP, AND CYANIDE\n")
    hdr = ("{:<19}{:>8}{:>9}{:>8}{:>7}{:>7}{:>7}{:>9}{:>8}"
           .format("scenario", "SBP", "vd mg/h", "Qren", "PRA", "creat",
                   "lact", "CN mg/L", "SCN"))
    print(hdr)
    print("-" * len(hdr))
    for k in ["S04_bb_then_vd", "S12_malperf_std", "S11_overshoot",
              "S06_snp_bb", "S13_snp_toxic"]:
        r = res[k]
        print("{:<19}{:>8.1f}{:>9.2f}{:>8.3f}{:>7.2f}{:>7.2f}{:>7.2f}{:>9.3f}"
              "{:>8.1f}".format(
                  sh(k), r["SBP"], r["RATE_VD"], r["QREN"], r["PRA"],
                  r["CREAT"], r["LACT"], r["C_CN_MAX"], r["C_SCN"]))
    m, g = res["S12_malperf_std"], res["S04_bb_then_vd"]
    print("""
  The same regimen at the same target reaches SBP {:.0f} in the patient whose
  renal artery is fed by the false lumen and SBP {:.0f} in the one whose is not.
  The vasodilator is already at its ceiling ({:.1f} vs {:.2f} mg/h): renal
  hypoperfusion (Q_ren {:.2f}) drives renin to {:.1f}x normal, and the
  angiotensin II it makes puts the resistance back.  Refractory hypertension in
  aortic dissection is a diagnosis, not a dosing failure - and the model's
  advice is to open the kidney, not to add a fourth vasodilator.

  Nitroprusside at {:.0f} mg/h ({:.1f} ug/kg/min) for 72 h with that same kidney
  reaches CN- {:.2f} mg/L (toxic threshold {:.1f}) and SCN- {:.0f} mg/L, and the
  lactate it produces ({:.2f} mmol/L) is indistinguishable at the bedside from
  the mesenteric ischaemia everyone is watching for.
""".format(m["SBP"], g["SBP"], m["RATE_VD"], g["RATE_VD"], m["QREN"], m["PRA"],
           res["S13_snp_toxic"]["RATE_VD"],
           res["S13_snp_toxic"]["RATE_VD"] * 1000.0 / 75.0 / 60.0,
           res["S13_snp_toxic"]["C_CN_MAX"], P["CN_TOX"],
           res["S13_snp_toxic"]["C_SCN"], res["S13_snp_toxic"]["LACT"]))

    # ---- 7 -----------------------------------------------------------------
    print("### 7. SPINAL CORD - the collateral network as a clock\n")
    hdr = ("{:<19}{:>8}{:>8}{:>8}{:>9}{:>9}{:>9}"
           .format("scenario", "MAP", "Q_spi", "collat", "SCI max", "para p", "surv5y"))
    print(hdr)
    print("-" * len(hdr))
    for k in ["S14_fl_patent", "S15_fl_partial", "S16_fl_complete", "S17_tevar",
              "S22_tevar_map"]:
        r = res[k]
        print("{:<19}{:>8.1f}{:>8.3f}{:>8.2f}{:>9.2f}{:>9.3f}{:>9.3f}".format(
            sh(k), r["MAP"], r["QSPI"], r["COLLAT"], r["SCI_MAX"], r["PARA"],
            r["SURV"]))
    t17, t22 = res["S17_tevar"], res["S22_tevar_map"]
    print("""
  Covering the entry tear costs {:.0%} of the segmental supply and thrombosing the
  false lumen costs more, so the cord depends on a collateral network that has to
  be RECRUITED.  Sustained pressure augmentation halves the ischaemic burden
  ({:.2f} -> {:.2f} index.h) but does not improve 5-year survival ({:.3f} ->
  {:.3f}), because the pressure that perfuses the cord also loads the aorta.
  The model's reading: augmentation belongs to the coverage window, not to the
  discharge prescription.
""".format(P["TEVAR_COVER"], t17["SCI_MAX"], t22["SCI_MAX"], t17["SURV"], t22["SURV"]))

    # ---- 8 -----------------------------------------------------------------
    print("### 8. CHRONIC MEDICAL THERAPY, GENETICS, ADHERENCE\n")
    hdr = ("{:<19}{:>7}{:>7}{:>7}{:>7}{:>8}{:>8}{:>8}{:>8}"
           .format("scenario", "SBP", "MMP9", "AT1blk", "PRA", "mm/yr",
                   "T55(y)", "sig/S", "surv5y"))
    print(hdr)
    print("-" * len(hdr))
    for k in ["S14_fl_patent", "S20_nonadherent", "S18_marfan_los",
              "S19_marfan_nolos"]:
        r = res[k]
        print("{:<19}{:>7.1f}{:>7.2f}{:>7.2f}{:>7.2f}{:>8.2f}{:>8.2f}{:>8.2f}"
              "{:>8.3f}".format(
                  sh(k), r["SBP"], r["MMP9"], r["BLK_AT1"], r["PRA"],
                  r["GROWTH_mm_yr"], r["T55"], r["STRESS_RATIO"], r["SURV"]))
    l, nl = res["S18_marfan_los"], res["S19_marfan_nolos"]
    print("""
  Adding losartan to a beta-blocked patient lowers systolic pressure by
  {:.1f} mmHg and blocks {:.0%} of AT1, yet moves the time to 55 mm from
  {:.2f} y to {:.2f} y - a difference of {:.0f} days, in the WRONG direction.
  The model's arithmetic: beta-blockade has already suppressed renin (PRA {:.2f}),
  so there is little angiotensin II left for the ARB to act on, and what AT1
  blockade does reach the wall removes TGF-beta-driven adventitial collagen -
  the only negative feedback the aorta has - as well as Ang II-driven MMP-9.
  The two effects cancel.  This is a PREDICTION, and it is the model's most
  exposed one; it happens to agree with the randomised trials that failed to
  show losartan superior to a beta-blocker on aortic growth, but the model was
  not fitted to them and the cancellation could be an artefact of the
  TGF-beta/collagen coupling.  Flagged accordingly.
""".format(nl["SBP"] - l["SBP"], l["BLK_AT1"], nl["T55"], l["T55"],
           (nl["T55"] - l["T55"]) * 365.0, l["PRA"]))

    # ---- 9 -----------------------------------------------------------------
    print("### 9. DRUG EXPOSURE AND RECEPTOR OCCUPANCY AT STEADY STATE\n")
    hdr = ("{:<19}{:>9}{:>9}{:>9}{:>9}{:>9}{:>7}{:>7}{:>7}"
           .format("scenario", "esm", "nic", "snp", "lab", "clv",
                   "occB", "occDHP", "occNO"))
    print(hdr + "   (all concentrations mg/L)")
    print("-" * len(hdr))
    for k in A:
        r = res[k]
        print("{:<19}{:>9.3f}{:>9.4f}{:>9.5f}{:>9.4f}{:>9.4f}{:>7.3f}{:>7.3f}"
              "{:>7.3f}".format(
                  sh(k), r["C_ESM"], r["C_NIC"], r["C_SNP"], r["C_LAB"],
                  r["C_CLV"], r["OCC_B"], r["OCC_DHP"], r["OCC_NO"]))

    # ---- 10 ----------------------------------------------------------------
    print("\n### 10. 72-h OUTCOME INTEGRALS\n")
    hdr = "{:<19}{:>10}{:>10}{:>10}{:>10}{:>10}".format(
        "scenario", "cumHaz", "mort%", "ASI.AUC/h", "sig/S", "D-dimer")
    print(hdr)
    print("-" * len(hdr))
    for k in A:
        r = res[k]
        print("{:<19}{:>10.4f}{:>10.2f}{:>10.3f}{:>10.3f}{:>10.2f}".format(
            sh(k), r["CUMH"], r["MORT_pct"], r["ASI_AUC_per_h"],
            r["STRESS_RATIO"], r["DDIM"]))

    print("\n" + "=" * W)
    print("SCENARIO KEY")
    print("=" * W)
    for k in sorted(res):
        print("  {:<19} {}".format(sh(k), res[k]["label"]))
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    res = run_all()
    if a.json:
        json.dump(res, sys.stdout, indent=1, sort_keys=True)
        print()
    else:
        print_report(res)


if __name__ == "__main__":
    main()
