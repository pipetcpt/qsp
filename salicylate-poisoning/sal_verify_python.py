#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sal_verify_python.py
====================
Independent Python/scipy re-implementation of the salicylate (aspirin)
poisoning QSP model distributed as `sal_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
The mrgsolve model is the deliverable; this file is the *check* on it.  Every
ODE is re-typed here from the mechanism rather than copied across, the two
implementations are required to agree, and the trajectories are then compared
against published human data.  Defects the check exposes are fixed in both
implementations and recorded in `sal_verification_output.txt` rather than
quietly tuned away.

THE THESIS THE MODEL IS BUILT TO TEST
-------------------------------------
Salicylate poisoning is graded clinically by a number - the total plasma
salicylate concentration - which is separated from the quantity that actually
does the killing (salicylate inside brain cells) by two multiplications that
the assay cannot see:

      C_brain  =  C_total  x  fu(C_total, pH)  x  f_n(pH_plasma)/f_n(pH_brain)
                             \_____________/     \___________________________/
                              free fraction:      pH partition of a weak acid
                              albumin saturates   (pKa 3.0).  Brain intra-
                              -> fu rises ~8x     cellular pH is BUFFERED, so
                              from 150 to 1000    it moves far less than plasma
                              mg/L                pH does.

Both multipliers rise as the patient deteriorates and neither appears on the
laboratory report.  Everything the model is asked to reproduce - the failure of
the Done nomogram, the lethality of chronic salicylism at "reassuring" levels,
the catastrophe of intubating a hyperventilating salicylate patient, the
primacy of bicarbonate over saline diuresis, and the potassium-dependence of
urinary alkalinisation - is meant to FOLLOW from those two factors rather than
to be asserted separately.

RUN
---
    python3 sal_verify_python.py
"""

import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

MW_SAL = 138.12          # salicylic acid,      g/mol
MW_ASA = 180.16          # acetylsalicylic acid, g/mol

# ============================================================================
# PARAMETERS
# ============================================================================
P0 = dict(
    BW      = 70.0,      # kg

    # --- absorption ------------------------------------------------------
    KA_ST   = 0.55,      # 1/h    stomach -> intestine (gastric emptying)
    KA_GUT  = 1.60,      # 1/h    intestine -> portal
    KA_CONC = 0.030,     # 1/h    concretion / bezoar dissolution
    F_CONC  = 0.0,       # fraction of the dose forming a concretion
    F_ORAL  = 0.90,      # oral bioavailability
    K_HYD   = 2.10,      # 1/h    systemic hydrolysis ASA -> salicylate
    K_AC    = 6.0,       # 1/h/g  adsorption of luminal drug onto charcoal
    KAC_OUT = 0.35,      # 1/h    charcoal transit out of the gut

    # --- distribution -----------------------------------------------------
    VC      = 5.0,       # L      central (measured) space
    VP      = 28.0,      # L      peripheral tissue water
    VCNS    = 1.20,      # L      brain water
    Q_P     = 45.0,      # L/h    central <-> peripheral exchange (free drug)
    PS_CNS  = 0.95,      # L/h    blood-brain PS product for the neutral acid
    KP_T    = 1.50,      # tissue trapping factor at the reference pH pair

    # --- plasma protein binding: two classes of albumin sites -------------
    B1      = 232.5,     # mg/L   high-affinity capacity  (fitted, see notes)
    KD1     = 11.67,     # mg/L
    B2      = 491.4,     # mg/L   low-affinity capacity
    KD2     = 380.8,     # mg/L
    ALB     = 40.0,      # g/L    albumin (capacities scale with this)
    KD_PH   = 2.00,      # acidaemia weakens binding, per pH unit below 7.40

    # --- ionisation --------------------------------------------------------
    PKA     = 3.00,      # salicylic acid

    # --- hepatic metabolism (Levy capacity-limited) -----------------------
    VMAX_SU = 100.0,     # mg/h   salicylurate (glycine conjugation)
    KM_SU   = 34.0,      # mg/L
    VMAX_PG = 27.0,      # mg/h   salicyl phenolic glucuronide
    KM_PG   = 550.0,     # mg/L
    CL_AG   = 0.045,     # L/h    salicyl acyl glucuronide (first order)
    CL_GA   = 0.025,     # L/h    gentisic acid (first order)
    GLY_T   = 9.0,       # h      glycine pool turnover time
    GLY_USE = 0.0030,    # fractional glycine drain per 100 mg/h salicylurate

    # --- renal --------------------------------------------------------------
    GFR0    = 7.2,       # L/h    (120 mL/min)
    A_REAB  = 260.0,     # L/h    tubular non-ionic back-diffusion constant
    HCO3_TH = 25.0,      # mmol/L renal bicarbonate threshold
    TH_K    = 0.22,      # threshold rise per mmol/L of hypokalaemia
    TH_VOL  = 0.55,      # threshold rise per unit fractional ECF contraction
    TH_CO2  = 0.20,      # mmol/L threshold shift per mmHg of PaCO2 (proximal
                         # HCO3 reabsorption varies directly with PCO2, so the
                         # patient's own hyperventilation LOWERS the threshold)
    UO_B    = 0.050,     # L/h    basal urine output
    UO_G    = 1.20,      # L/h    per L of ECF expansion
    UK_B    = 20.0,      # mmol/L basal urinary potassium
    UK_HCO3 = 0.35,      # mmol/L urinary K per mmol/L urinary bicarbonate
    UPH_OFF = 1.35,      # mmol/L non-bicarbonate buffer floor of the urine

    # --- acid-base -----------------------------------------------------------
    VBIC    = 35.0,      # L      bicarbonate space (0.5 L/kg)
    TAU_CO2 = 0.12,      # h      lung time constant
    EMAX_R  = 1.65,      # maximal fractional rise in alveolar ventilation
    EC50_R  = 22.0,      # mg/L   brain salicylate for half-maximal drive
    K_CHEM  = 7.0,       # ventilatory drive per unit fall in brain pH
    KF_ON   = 0.030,     # 1/h    respiratory muscle fatigue accrual
    KF_OFF  = 0.10,      # 1/h    fatigue recovery
    NAG0    = 2.9,       # mmol/h dietary fixed acid load
    NAG_C   = 20.0,      # mmol/h extra renal acid excretion per pH unit below 7.4
    PH_CAP  = 7.60,      # arterial pH at which the NaHCO3 infusion is held
    PH_CAPW = 0.015,     # width of that titration switch

    # --- intracellular pH: buffer bases and buffer powers --------------------
    # Brain and muscle intracellular pH are NOT computed from a fixed
    # intracellular bicarbonate.  They are computed from a buffer base that
    # follows plasma bicarbonate slowly, together with a non-bicarbonate buffer
    # power - which is why an acute rise in PaCO2 moves plasma pH far more than
    # it moves intracellular pH, and therefore drives the partition.
    BB_RATIO = 0.536,    # brain buffer base / plasma bicarbonate at equilibrium
    TAU_BB   = 5.8,      # h  (t1/2 4 h) brain equilibration with plasma HCO3
    BETA_B   = 35.0,     # mmol/L per pH unit, brain non-bicarbonate buffers
    PHREF_B  = 7.05,     # reference brain intracellular pH
    DPCO2_B  = 8.0,      # mmHg brain tissue PCO2 above arterial
    TB_RATIO = 0.458,    # muscle buffer base / plasma bicarbonate
    TAU_TB   = 4.3,      # h
    BETA_T   = 30.0,     # mmol/L per pH unit, muscle
    PHREF_T  = 7.00,
    DPCO2_T  = 6.0,      # mmHg

    # --- mitochondrial uncoupling and organic acids ---------------------------
    EC50_U  = 200.0,     # mg/L tissue salicylate, half-maximal uncoupling
    K_LAC   = 26.0,      # mmol/h lactate production at full uncoupling
    CL_LAC  = 0.35,      # 1/h
    K_KET   = 12.0,      # mmol/h ketoacid production at full uncoupling
    CL_KET  = 0.15,      # 1/h

    # --- potassium and volume --------------------------------------------------
    VK      = 250.0,     # L      apparent K space (100 mmol ~ 0.4 mmol/L)
    K_DIET  = 3.6,       # mmol/h
    KPH_SH  = 3.0,       # mmol/L plasma K per pH unit (transcellular shift)
    VECF0   = 14.0,      # L
    R_INTAKE = 0.110,    # L/h    oral/maintenance water intake
    INSENS  = 0.045,     # L/h    insensible loss
    SWEAT   = 0.16,      # L/h    at full uncoupling
    GFR_VOL = 3.0,       # GFR sensitivity to fractional ECF contraction

    # --- brain glucose -----------------------------------------------------------
    GLC_P   = 5.5,       # mmol/L plasma glucose (raised by dextrose)
    TMAX_G  = 60.0,      # mmol/L/h symmetric GLUT1 transport
    KT_G    = 6.0,       # mmol/L
    CMR0    = 18.0,      # mmol/L/h cerebral glucose consumption
    CMR_U   = 0.50,      # extra consumption at full uncoupling

    # --- thermogenesis, lung ------------------------------------------------------
    K_TEMP  = 0.90, CL_TEMP = 0.50,
    K_LUNG  = 0.055, CL_LUNG = 0.030, EC50_LUNG = 150.0,

    # --- toxicodynamics -------------------------------------------------------------
    EC50_CNS = 95.0,     # mg/L brain salicylate, half-maximal CNS depression
    HILL_CNS = 2.2,
    CNS_THR  = 45.0,     # mg/L brain: threshold for cumulative injury
    GLC_THR  = 0.80,     # mmol/L brain glucose below which injury accelerates

    # --- interventions (driven by the scenario scheduler) -----------------------------
    R_BIC   = 0.0,       # mmol/h sodium bicarbonate infusion
    R_FLUID = 0.0,       # L/h    crystalloid (including bicarbonate vehicle)
    R_KCL   = 0.0,       # mmol/h potassium replacement
    CL_HD   = 0.0,       # L/h    extracorporeal clearance (0 = off)
    HD_BIC  = 0.0,       # mmol/h bicarbonate gained from the dialysate
    VENT    = 0.0,       # 0 spontaneous, 1 controlled mechanical ventilation
    VA_SET  = 1.0,       # relative alveolar ventilation when VENT = 1
    SED     = 0.0,       # fractional respiratory depression (sedatives)
    AGE_F   = 0.0,       # 0 young, 1 elderly (lung-injury susceptibility)
)

IDX = dict(
    AST=0, AGUT=1, ACONC=2, AASA=3, ACENT=4, APER=5, ACNS=6,
    ASU=7, APG=8, AAG=9, AGA=10, AUR=11, AHD=12,
    HCO3=13, PACO2=14, BBB=15, BBT=16, LAC=17, KET=18, KB=19,
    VECF=20, GLY=21, TEMP=22, GLUB=23, LUNG=24, CNSI=25,
    AAC=26, FATIG=27,
)
NST = 28


# ============================================================================
# ALGEBRA  (identical formulations are used in the mrgsolve $GLOBAL block)
# ============================================================================
def free_conc(ctot, ph, p):
    """Unbound salicylate (mg/L) from total plasma concentration.

    bound(Cf) = B1*Cf/(K1+Cf) + B2*Cf/(K2+Cf), capacities scaled to albumin and
    dissociation constants widened by acidaemia.  Solved by Newton iteration -
    the same routine is used in the C++ side of the mrgsolve model.
    """
    if ctot <= 1e-12:
        return 0.0
    sc = p['ALB'] / 40.0
    b1, b2 = p['B1'] * sc, p['B2'] * sc
    m = 1.0 + p['KD_PH'] * max(0.0, 7.40 - ph)
    k1, k2 = p['KD1'] * m, p['KD2'] * m
    cf = ctot * 0.3 + 1e-6
    for _ in range(60):
        f = cf + b1 * cf / (k1 + cf) + b2 * cf / (k2 + cf) - ctot
        d = 1.0 + b1 * k1 / (k1 + cf) ** 2 + b2 * k2 / (k2 + cf) ** 2
        step = f / d
        cf_new = cf - step
        if cf_new <= 0.0:
            cf_new = cf * 0.5
        if abs(cf_new - cf) < 1e-10 * max(1.0, cf):
            cf = cf_new
            break
        cf = cf_new
    return min(max(cf, 0.0), ctot)


def f_neutral(ph, p):
    """Fraction present as the un-ionised, membrane-permeant acid."""
    return 1.0 / (1.0 + 10.0 ** (ph - p['PKA']))


def ph_from(hco3, pco2):
    return 6.10 + math.log10(max(hco3, 0.30) / (0.0301 * max(pco2, 5.0)))


def ph_intracell(buf, pco2, beta, phref):
    """Intracellular pH of a compartment with a non-bicarbonate buffer power.

    [HCO3-]i(pH) = BUF - beta*(pH - phref)          (closed-system buffering)
    pH           = 6.1 + log10([HCO3-]i / (0.0301*PCO2))

    Solved by Newton.  The consequence that matters clinically: raising PCO2
    from 25 to 60 mmHg moves plasma pH by ~0.35 but intracellular pH by <0.1,
    so the plasma-to-cell partition of a weak acid swings sharply.
    """
    x = phref
    for _ in range(40):
        h = buf - beta * (x - phref)
        if h < 0.20:
            h = 0.20
        g = 6.10 + math.log10(h / (0.0301 * max(pco2, 5.0)))
        fx = g - x
        dfx = -0.434294 * beta / h - 1.0
        step = fx / dfx
        xn = x - step
        if xn < 5.0:
            xn = 5.0
        if xn > 8.5:
            xn = 8.5
        if abs(xn - x) < 1e-10:
            x = xn
            break
        x = xn
    return x


def urine_ph(u_hco3, p):
    return min(8.20, max(4.60,
               6.10 + math.log10((u_hco3 + p['UPH_OFF']) / 1.354)))


# ============================================================================
# DERIVED QUANTITIES
# ============================================================================
def derive(y, p):
    d = {}
    ctot = max(y[IDX['ACENT']], 0.0) / p['VC']
    hco3 = max(y[IDX['HCO3']], 1.0)
    paco2 = max(y[IDX['PACO2']], 5.0)
    ph = ph_from(hco3, paco2)
    cf = free_conc(ctot, ph, p)
    fu = cf / ctot if ctot > 1e-9 else free_conc(1.0, ph, p)

    phb = ph_intracell(max(y[IDX['BBB']], 1.0), paco2 + p['DPCO2_B'],
                       p['BETA_B'], p['PHREF_B'])
    pht = ph_intracell(max(y[IDX['BBT']], 1.0), paco2 + p['DPCO2_T'],
                       p['BETA_T'], p['PHREF_T'])
    ccns = max(y[IDX['ACNS']], 0.0) / p['VCNS']
    ctis = max(y[IDX['APER']], 0.0) / p['VP']

    kp = max(0.8, y[IDX['KB']] - p['KPH_SH'] * (ph - 7.40))

    vecf = max(y[IDX['VECF']], 6.0)
    contract = max(0.0, (p['VECF0'] - vecf) / p['VECF0'])
    gfrf = min(1.10, max(0.12, 1.0 - p['GFR_VOL'] * contract))
    gfr = p['GFR0'] * gfrf
    uo = min(0.80, max(0.004,
                       (p['UO_B'] + p['UO_G'] * (vecf - p['VECF0'])) * gfrf))

    thr = (p['HCO3_TH'] * (1.0 + p['TH_K'] * max(0.0, 4.0 - kp)
                           + p['TH_VOL'] * contract)
           + p['TH_CO2'] * (paco2 - 40.0))
    j_uhco3 = max(0.0, gfr * (hco3 - thr))
    u_hco3 = j_uhco3 / uo
    uph = urine_ph(u_hco3, p)

    fn_u = f_neutral(uph, p)
    frac_reab = 1.0 - math.exp(-p['A_REAB'] * fn_u / uo)
    cl_renal = gfr * fu * (1.0 - frac_reab)

    u = ctis / (p['EC50_U'] + ctis)

    # --- partition coefficients (the heart of the model) -----------------
    fn_p, fn_b, fn_t = f_neutral(ph, p), f_neutral(phb, p), f_neutral(pht, p)
    d['kp_brain'] = fn_p / fn_b                # brain water : free plasma
    d['kp_tissue'] = fn_p / fn_t
    d['ratio_meas'] = ccns / ctot if ctot > 1e-9 else 0.0   # brain : TOTAL plasma

    # --- ventilation ------------------------------------------------------
    drive_sal = p['EMAX_R'] * ccns / (p['EC50_R'] + ccns)
    drive_ac = min(1.5, max(0.0, p['K_CHEM'] * (p['PHREF_B'] - phb)))
    va_spont = ((1.0 + drive_sal + drive_ac)
                * (1.0 - 0.85 * y[IDX['FATIG']])
                * (1.0 - 0.45 * y[IDX['LUNG']])
                * (1.0 - p['SED']))
    va = max(p['VA_SET'] if p['VENT'] > 0.5 else va_spont, 0.15)
    vco2 = 1.0 + 0.60 * u
    paco2_ss = min(120.0, max(8.0, 40.0 * vco2 / va))

    d['bic_gate'] = 1.0 / (1.0 + math.exp((ph - p['PH_CAP']) / p['PH_CAPW']))
    d.update(ctot=ctot, cfree=cf, fu=fu, ph=ph, phb=phb, pht=pht, ccns=ccns,
             ctis=ctis, hco3=hco3, paco2=paco2, kp=kp, gfr=gfr, gfrf=gfrf,
             uo=uo, uph=uph, u_hco3=u_hco3, frac_reab=frac_reab,
             cl_renal=cl_renal, uncouple=u, va=va, va_spont=va_spont,
             paco2_ss=paco2_ss, thr=thr, j_uhco3=j_uhco3, contract=contract)

    # --- elimination fluxes -------------------------------------------------
    gly = max(0.05, y[IDX['GLY']])
    d['v_su'] = p['VMAX_SU'] * gly * ctot / (p['KM_SU'] + ctot)
    d['v_pg'] = p['VMAX_PG'] * ctot / (p['KM_PG'] + ctot)
    d['v_ag'] = p['CL_AG'] * ctot
    d['v_ga'] = p['CL_GA'] * ctot
    d['v_ren'] = cl_renal * ctot
    d['v_hd'] = p['CL_HD'] * ctot
    tot = d['v_su'] + d['v_pg'] + d['v_ag'] + d['v_ga'] + d['v_ren'] + d['v_hd']
    d['cl_tot'] = tot / ctot if ctot > 1e-9 else 0.0
    d['f_renal'] = d['v_ren'] / tot if tot > 1e-9 else 0.0

    # --- clinical readouts ---------------------------------------------------
    d['anion_gap'] = 12.0 + (y[IDX['LAC']] - 1.0) + (y[IDX['KET']] - 0.10) \
        + ctot / MW_SAL
    d['cns_score'] = 100.0 * ccns ** p['HILL_CNS'] / (
        p['EC50_CNS'] ** p['HILL_CNS'] + ccns ** p['HILL_CNS'])
    d['tinnitus'] = 100.0 * ccns ** 2 / (18.0 ** 2 + ccns ** 2)
    d['vd_app'] = ((y[IDX['ACENT']] + y[IDX['APER']] + y[IDX['ACNS']]) / ctot
                   if ctot > 1e-6 else p['VC'])
    d['burden'] = (y[IDX['ACENT']] + y[IDX['APER']] + y[IDX['ACNS']]
                   + (y[IDX['AST']] + y[IDX['AGUT']] + y[IDX['ACONC']]
                      + y[IDX['AASA']]) * MW_SAL / MW_ASA)
    sev = (0.45 * d['cns_score'] / 100.0
           + 0.25 * min(1.0, max(0.0, (7.35 - ph) / 0.30))
           + 0.15 * min(1.0, max(0.0, (p['GLC_THR'] - y[IDX['GLUB']]) / p['GLC_THR']))
           + 0.15 * min(1.0, max(0.0, (y[IDX['TEMP']] - 37.5) / 2.5)))
    d['severity'] = 100.0 * min(1.0, sev)
    return d


# ============================================================================
# RIGHT-HAND SIDE
# ============================================================================
def rhs(t, y, p):
    d = derive(y, p)
    dy = np.zeros(NST)

    # ---- gut -------------------------------------------------------------
    ac = max(y[IDX['AAC']], 0.0)
    kbind = p['K_AC'] * ac / (1.0 + p['K_AC'] * ac / 3.0)
    dy[IDX['AST']] = -p['KA_ST'] * y[IDX['AST']] - kbind * y[IDX['AST']]
    dy[IDX['ACONC']] = -p['KA_CONC'] * y[IDX['ACONC']]
    dy[IDX['AGUT']] = (p['KA_ST'] * y[IDX['AST']]
                       + p['KA_CONC'] * y[IDX['ACONC']]
                       - p['KA_GUT'] * y[IDX['AGUT']]
                       - kbind * y[IDX['AGUT']])
    dy[IDX['AAC']] = -p['KAC_OUT'] * ac

    absorbed = p['F_ORAL'] * p['KA_GUT'] * y[IDX['AGUT']]     # mg/h as ASA
    dy[IDX['AASA']] = absorbed - p['K_HYD'] * y[IDX['AASA']]
    to_sal = p['K_HYD'] * y[IDX['AASA']] * MW_SAL / MW_ASA    # mg/h salicylate

    # ---- disposition -------------------------------------------------------
    kp_t = p['KP_T'] * d['kp_tissue'] / (f_neutral(7.40, p) / f_neutral(7.00, p))
    j_per = p['Q_P'] * (d['cfree'] - y[IDX['APER']] / (p['VP'] * max(kp_t, 1e-3)))
    j_cns = p['PS_CNS'] * (d['cfree'] - d['ccns'] / max(d['kp_brain'], 1e-3))

    dy[IDX['APER']] = j_per
    dy[IDX['ACNS']] = j_cns
    dy[IDX['ACENT']] = (to_sal - j_per - j_cns
                        - d['v_su'] - d['v_pg'] - d['v_ag'] - d['v_ga']
                        - d['v_ren'] - d['v_hd'])
    dy[IDX['ASU']] = d['v_su']
    dy[IDX['APG']] = d['v_pg']
    dy[IDX['AAG']] = d['v_ag']
    dy[IDX['AGA']] = d['v_ga']
    dy[IDX['AUR']] = d['v_ren']
    dy[IDX['AHD']] = d['v_hd']
    dy[IDX['GLY']] = ((1.0 - y[IDX['GLY']]) / p['GLY_T']
                      - p['GLY_USE'] * d['v_su'] / 100.0)

    # ---- organic acids and bicarbonate --------------------------------------
    dlac = p['K_LAC'] * d['uncouple'] / p['VBIC'] - p['CL_LAC'] * (y[IDX['LAC']] - 1.0)
    dket = p['K_KET'] * d['uncouple'] / p['VBIC'] - p['CL_KET'] * (y[IDX['KET']] - 0.10)
    dy[IDX['LAC']] = dlac
    dy[IDX['KET']] = dket

    sal_acid = to_sal / MW_SAL                                # mmol/h of acid
    renal_comp = p['NAG_C'] * max(0.0, 7.40 - d['ph']) * d['gfrf']
    # Clinical titration rule (Proudfoot 2004): the bicarbonate infusion is held
    # once the arterial pH reaches the ceiling.  This is what converts a
    # potassium problem into a therapeutic failure rather than into alkalaemia.
    dy[IDX['HCO3']] = (-dlac - dket
                       + (p['R_BIC'] * d['bic_gate'] + p['HD_BIC'] + renal_comp
                          - sal_acid - d['j_uhco3']
                          + p['NAG0'] * (d['gfrf'] - 1.0)) / p['VBIC'])

    # ---- ventilation ----------------------------------------------------------
    dy[IDX['PACO2']] = (d['paco2_ss'] - y[IDX['PACO2']]) / p['TAU_CO2']
    dy[IDX['FATIG']] = (p['KF_ON'] * max(0.0, d['va_spont'] - 2.4)
                        * (1.0 - y[IDX['FATIG']])
                        - p['KF_OFF'] * y[IDX['FATIG']])

    # ---- intracellular buffer bases -------------------------------------------
    dy[IDX['BBB']] = (p['BB_RATIO'] * d['hco3'] - y[IDX['BBB']]) / p['TAU_BB']
    dy[IDX['BBT']] = (p['TB_RATIO'] * d['hco3'] - y[IDX['BBT']]) / p['TAU_TB']

    # ---- potassium and volume ----------------------------------------------------
    u_k = min(120.0, p['UK_B'] + p['UK_HCO3'] * d['u_hco3'])
    dy[IDX['KB']] = (p['K_DIET'] + p['R_KCL'] - d['uo'] * u_k) / p['VK']
    dy[IDX['VECF']] = (p['R_INTAKE'] + p['R_FLUID'] - d['uo'] - p['INSENS']
                       - p['SWEAT'] * d['uncouple'])

    # ---- brain glucose ---------------------------------------------------------------
    gp, gb = p['GLC_P'], max(y[IDX['GLUB']], 0.0)
    dy[IDX['GLUB']] = (p['TMAX_G'] * (gp / (p['KT_G'] + gp) - gb / (p['KT_G'] + gb))
                       - p['CMR0'] * (1.0 + p['CMR_U'] * d['uncouple']))

    # ---- temperature, lung, cumulative CNS injury -----------------------------------
    dy[IDX['TEMP']] = p['K_TEMP'] * d['uncouple'] - p['CL_TEMP'] * (y[IDX['TEMP']] - 37.0)
    dy[IDX['LUNG']] = (p['K_LUNG'] * (1.0 + p['AGE_F'])
                       * d['ccns'] / (p['EC50_LUNG'] + d['ccns'])
                       * (1.0 - y[IDX['LUNG']]) - p['CL_LUNG'] * y[IDX['LUNG']])
    dy[IDX['CNSI']] = (max(0.0, d['ccns'] - p['CNS_THR']) / 100.0
                       * (1.0 + 2.0 * max(0.0, (p['GLC_THR'] - gb) / p['GLC_THR'])))
    return dy


# ============================================================================
# INITIAL CONDITIONS AND SCENARIO ENGINE
# ============================================================================
def y0(p, hco3=24.0, paco2=40.0):
    y = np.zeros(NST)
    y[IDX['HCO3']] = hco3
    y[IDX['PACO2']] = paco2
    y[IDX['BBB']] = p['BB_RATIO'] * hco3
    y[IDX['BBT']] = p['TB_RATIO'] * hco3
    y[IDX['LAC']] = 1.0
    y[IDX['KET']] = 0.10
    y[IDX['KB']] = 4.0
    y[IDX['VECF']] = p['VECF0']
    y[IDX['GLY']] = 1.0
    y[IDX['TEMP']] = 37.0
    y[IDX['GLUB']] = 1.30
    return y


class Sched:
    """Piecewise-constant intervention schedule.

    Items are (time, {key: value}).  A key that names a parameter changes that
    parameter from that time on.  Two special keys act instantaneously on the
    state: 'AAC' adds grams of activated charcoal, '+HCO3' adds a bicarbonate
    bolus expressed in mmol.
    """

    def __init__(self, items=None):
        self.items = sorted(items or [], key=lambda x: x[0])

    def times(self):
        return [t for t, _ in self.items]


def _params_at(p, sched, t):
    q = dict(p)
    for td, upd in sched.items:
        if td <= t + 1e-9:
            for k, v in upd.items():
                if k not in ('AAC', '+HCO3'):
                    q[k] = v
    return q


def simulate(dose_g=0.0, tmax=48.0, p=None, sched=None, doses=None,
             form='plain', n=1400, hco3_0=24.0, paco2_0=40.0):
    p = dict(P0 if p is None else p)
    sched = sched or Sched()
    if doses is None:
        doses = [(0.0, dose_g)] if dose_g > 0 else []

    y = y0(p, hco3_0, paco2_0)
    if form == 'enteric':
        p['KA_ST'], p['KA_GUT'] = 0.16, 0.40
    elif form == 'massive':
        p['KA_ST'], p['KA_GUT'] = 0.40, 1.20

    evts = sorted(set([0.0] + [t for t, _ in doses] + sched.times() + [tmax]))
    ts, ys = [], []
    for i in range(len(evts) - 1):
        t0, t1 = evts[i], evts[i + 1]
        for td, g in doses:
            if abs(td - t0) < 1e-9 and g > 0:
                mg = g * 1000.0
                y[IDX['ACONC']] += mg * p['F_CONC']
                y[IDX['AST']] += mg * (1.0 - p['F_CONC'])
        for td, upd in sched.items:
            if abs(td - t0) < 1e-9:
                for k, v in upd.items():
                    if k == 'AAC':
                        y[IDX['AAC']] += v
                    elif k == '+HCO3':
                        y[IDX['HCO3']] += v / p['VBIC']
                    else:
                        p[k] = v
        if t1 <= t0:
            continue
        nn = max(5, int(n * (t1 - t0) / tmax))
        sol = solve_ivp(rhs, (t0, t1), y, args=(p,),
                        t_eval=np.linspace(t0, t1, nn),
                        method='LSODA', rtol=1e-7, atol=1e-9, max_step=0.25)
        if not sol.success:
            raise RuntimeError("integration failed: " + sol.message)
        ts.append(sol.t if i == 0 else sol.t[1:])
        ys.append(sol.y if i == 0 else sol.y[:, 1:])
        y = sol.y[:, -1].copy()

    T = np.concatenate(ts)
    Y = np.concatenate(ys, axis=1)
    D = {}
    for j in range(Y.shape[1]):
        dd = derive(Y[:, j], _params_at(p, sched, T[j]))
        for k, v in dd.items():
            D.setdefault(k, []).append(v)
    D = {k: np.array(v) for k, v in D.items()}
    for k, i in IDX.items():
        D[k] = Y[i, :]
    D['time'] = T
    return D


# ============================================================================
# REPORT HELPERS
# ============================================================================
LINE = "=" * 78
OK, BAD = "  [OK]  ", "  [FAIL]"
_results = []


def check(name, value, lo, hi, unit="", note=""):
    good = (lo is None or value >= lo) and (hi is None or value <= hi)
    rng = (f"target {lo:g}-{hi:g}{unit}" if lo is not None and hi is not None
           else (f"target >{lo:g}{unit}" if lo is not None else ""))
    print(f"{OK if good else BAD} {name:<48s} = {value:>9.4g} {unit:<7s} {rng} {note}")
    _results.append((name, good))
    return good


def info(name, value, unit=""):
    print(f"        {name:<48s} = {value:>9.4g} {unit}")


def hdr(t):
    print("\n" + LINE + "\n" + t + "\n" + LINE)


def peak(D, k):
    return float(np.max(D[k]))


def at(D, k, t):
    return float(np.interp(t, D['time'], D[k]))


def halflife(D, k, t0, t1):
    c0, c1 = at(D, k, t0), at(D, k, t1)
    if c0 <= 0 or c1 <= 0 or c1 >= c0:
        return float('nan')
    return (t1 - t0) * math.log(2.0) / math.log(c0 / c1)


# Standard adult regimens -----------------------------------------------------
def rx_supportive(t=4.0):
    """Crystalloid at 125 mL/h - what every poisoned patient actually gets."""
    return [(t, dict(R_FLUID=0.125, R_KCL=4.0))]


def rx_alkalinise(t=4.0, kcl=10.0):
    """Standard regimen (Proudfoot 2004 / AAPCC): 1.5 mmol/kg NaHCO3 bolus, then
    1 L of 150 mmol/L over 4 h, then the same solution at 250 mL/h."""
    return [(t, {'+HCO3': 100.0, 'R_BIC': 75.0, 'R_FLUID': 0.25, 'R_KCL': kcl}),
            (t + 4.0, {'R_BIC': 37.5, 'R_FLUID': 0.25, 'R_KCL': kcl})]


# ============================================================================
# MAIN REPORT
# ============================================================================
def main():
    print(LINE)
    print("SALICYLATE POISONING QSP MODEL - INDEPENDENT VERIFICATION (Python/scipy)")
    print(LINE)
    print("28 ODEs re-implemented from mechanism, integrated with LSODA.")
    print("Every target below is published human data cited in sal_references.md.")

    # ----------------------------------------------------------------
    hdr("1. STATIC ALGEBRA - the two multipliers the assay cannot see")
    # ----------------------------------------------------------------
    print("(a) albumin saturation")
    for ct in (50, 150, 300, 500, 800, 1000):
        info(f"free fraction at {ct} mg/L, pH 7.40",
             100 * free_conc(ct, 7.40, P0) / ct, "%")
    fu150 = free_conc(150, 7.40, P0) / 150
    fu500 = free_conc(500, 7.40, P0) / 500
    fu800 = free_conc(800, 7.40, P0) / 800
    check("free fraction, therapeutic 150 mg/L", 100 * fu150, 5.0, 12.0, "%",
          "lit. 8-10%")
    check("free fraction, 500 mg/L", 100 * fu500, 20.0, 40.0, "%", "lit. 25-40%")
    check("free fraction, 800 mg/L", 100 * fu800, 30.0, 55.0, "%", "lit. >30%")
    check("free-fraction fold-rise, 150 -> 800 mg/L", fu800 / fu150, 3.5, 8.0, "x",
          "multiplier 1")
    check("acidaemia (pH 7.10) further raises free fraction",
          (free_conc(500, 7.10, P0) / 500) / fu500, 1.10, 1.6, "x")

    print("\n(b) pH partition of a weak acid, pKa 3.0")
    info("un-ionised fraction at pH 7.40", f_neutral(7.40, P0) * 1e6, "ppm")
    info("un-ionised fraction at pH 7.10", f_neutral(7.10, P0) * 1e6, "ppm")
    kp_norm = f_neutral(7.40, P0) / f_neutral(7.02, P0)
    check("brain:free-plasma partition at pH 7.40 / brain 7.02", kp_norm,
          0.30, 0.55, "", "= 10^(pHb - pHp)")

    print("\n(c) intracellular buffering - why CO2 and HCO3 are not symmetric")
    bb0 = P0['BB_RATIO'] * 24.0
    for pco2 in (20, 25, 40, 60, 80):
        info(f"brain pHi at PaCO2 {pco2} mmHg (buffer base fixed)",
             ph_intracell(bb0, pco2 + P0['DPCO2_B'], P0['BETA_B'], P0['PHREF_B']))
    php_lo = ph_from(24.0, 25.0)
    php_hi = ph_from(24.0, 60.0)
    phb_lo = ph_intracell(bb0, 25 + 8, P0['BETA_B'], P0['PHREF_B'])
    phb_hi = ph_intracell(bb0, 60 + 8, P0['BETA_B'], P0['PHREF_B'])
    info("plasma pH swing, PaCO2 25 -> 60", php_lo - php_hi)
    info("brain pHi swing over the same change", phb_lo - phb_hi)
    check("plasma pH moves more than brain pHi does",
          (php_lo - php_hi) / (phb_lo - phb_hi), 2.5, 12.0, "x",
          "31P-MRS: brain pHi is defended")
    part_lo = f_neutral(php_lo, P0) / f_neutral(phb_lo, P0)
    part_hi = f_neutral(php_hi, P0) / f_neutral(phb_hi, P0)
    check("brain partition rises when PaCO2 is 'normalised'",
          part_hi / part_lo, 1.5, 4.0, "x", "multiplier 2 - the lethal one")

    # ----------------------------------------------------------------
    hdr("2. THERAPEUTIC PHARMACOKINETICS")
    # ----------------------------------------------------------------
    D = simulate(dose_g=0.65, tmax=24)
    check("single 650 mg aspirin: Cmax salicylate", peak(D, 'ctot'), 25.0, 70.0,
          "mg/L", "lit. 30-60")
    check("single 650 mg: apparent t1/2 (4-10 h)", halflife(D, 'ctot', 4, 10),
          2.0, 5.0, "h", "lit. 2-4.5")
    check("apparent Vd, therapeutic", at(D, 'vd_app', 8) / P0['BW'], 0.10, 0.25,
          "L/kg", "lit. 0.15-0.2")
    info("urine pH, untreated normal subject", at(D, 'uph', 8))

    doses = [(t, 0.65) for t in np.arange(0, 120, 4)]
    Dc = simulate(doses=doses, tmax=132)
    css = at(Dc, 'ctot', 120)
    check("aspirin 3.9 g/day x 5 d: steady-state salicylate", css, 180.0, 340.0,
          "mg/L", "lit. 150-300 anti-inflammatory")
    check("fraction excreted unchanged, acid urine", 100 * at(Dc, 'f_renal', 120),
          1.0, 25.0, "%", "lit. 2-30%, pH-dependent")
    info("urine output at steady state", 1000 * at(Dc, 'uo', 120), "mL/h")
    check("capacity-limited: 6x the dose gives >6x the level",
          css / peak(D, 'ctot'), 3.0, 14.0, "x")
    info("free salicylate at steady state", at(Dc, 'cfree', 120), "mg/L")
    info("brain salicylate at steady state", at(Dc, 'ccns', 120), "mg/L")
    info("tinnitus score at steady state", at(Dc, 'tinnitus', 120), "/100")

    # ----------------------------------------------------------------
    hdr("3. ACUTE MASSIVE INGESTION (30 g) ON SUPPORTIVE CARE")
    # ----------------------------------------------------------------
    A = simulate(dose_g=30.0, tmax=72, form='massive',
                 sched=Sched(rx_supportive()))
    check("peak plasma salicylate", peak(A, 'ctot'), 500.0, 1000.0, "mg/L",
          "lit. severe >700")
    check("PaCO2 at 6 h (respiratory alkalosis comes first)", at(A, 'paco2', 6),
          18.0, 34.0, "mmHg", "lit. 20-30")
    check("arterial pH at 6 h is alkalaemic", at(A, 'ph', 6), 7.42, 7.60, "",
          "lit. mixed disturbance, pH high early")
    check("nadir bicarbonate", float(np.min(A['hco3'])), 8.0, 20.0, "mmol/L")
    check("peak anion gap", peak(A, 'anion_gap'), 18.0, 40.0, "mmol/L")
    check("apparent t1/2 over 24-48 h", halflife(A, 'ctot', 24, 48), 15.0, 45.0,
          "h", "lit. 20-36 in overdose")
    check("peak apparent Vd", peak(A, 'vd_app') / P0['BW'], 0.20, 0.55, "L/kg",
          "lit. rises with acidaemia")
    check("neuroglycopenia despite normal plasma glucose",
          float(np.min(A['GLUB'])), 0.0, 0.95, "mmol/L", "plasma glucose fixed 5.5")
    info("peak brain salicylate", peak(A, 'ccns'), "mg/L")
    info("peak temperature", peak(A, 'TEMP'), "degC")
    info("peak free salicylate", peak(A, 'cfree'), "mg/L")
    info("fraction of the dose excreted unchanged by 72 h",
         100 * A['AUR'][-1] / (30000 * MW_SAL / MW_ASA * P0['F_ORAL']), "%")

    NOFL = simulate(dose_g=30.0, tmax=72, form='massive')
    check("withholding fluid collapses renal clearance",
          min(1e4, at(A, 'cl_renal', 24) / max(at(NOFL, 'cl_renal', 24), 1e-9)),
          3.0, None, "x", "volume depletion -> aciduria + low GFR")
    info("urine output at 24 h, fluids", 1000 * at(A, 'uo', 24), "mL/h")
    info("urine output at 24 h, no fluids", 1000 * at(NOFL, 'uo', 24), "mL/h")
    info("plasma salicylate at 48 h, fluids", at(A, 'ctot', 48), "mg/L")
    info("plasma salicylate at 48 h, no fluids", at(NOFL, 'ctot', 48), "mg/L")

    # ----------------------------------------------------------------
    hdr("4. THE CENTRAL EXPERIMENT - intubating to a 'normal' PaCO2")
    # ----------------------------------------------------------------
    print("One switch at t = 8 h, nothing else changed.")
    print("  Arm A  ventilator set to a conventional PaCO2   (VA_SET = 1.0)")
    print("  Arm B  minute ventilation matched to the patient's own drive")
    base = simulate(dose_g=30.0, tmax=36, form='massive',
                    sched=Sched(rx_supportive()))
    va_own = at(base, 'va_spont', 8.0)
    armA = simulate(dose_g=30.0, tmax=36, form='massive',
                    sched=Sched(rx_supportive() + [(8.0, dict(VENT=1.0, VA_SET=1.0))]))
    armB = simulate(dose_g=30.0, tmax=36, form='massive',
                    sched=Sched(rx_supportive() + [(8.0, dict(VENT=1.0, VA_SET=va_own))]))
    info("the patient's own alveolar ventilation at 8 h", va_own, "x normal")
    print(f"   {'arm':<32s}{'pH':>6s}{'PaCO2':>7s}{'plasma':>9s}{'brain':>8s}"
          f"{'br/pl':>8s}{'CNS':>7s}")
    for nm, S in (("not intubated", base), ("intubated, PaCO2 normalised", armA),
                  ("intubated, ventilation matched", armB)):
        print(f"   {nm:<32s}{at(S,'ph',9):6.2f}{at(S,'paco2',9):7.1f}"
              f"{at(S,'ctot',9):9.1f}{at(S,'ccns',9):8.1f}"
              f"{at(S,'ratio_meas',9):8.3f}{at(S,'cns_score',9):7.1f}")
    print("   time course after the switch (arm A / arm B):")
    print(f"   {'t (h)':>6s}{'plasma A':>10s}{'plasma B':>10s}{'brain A':>9s}"
          f"{'brain B':>9s}{'ratio':>8s}")
    for t in (8.0, 9.0, 12.0, 18.0, 24.0, 36.0):
        print(f"   {t:6.0f}{at(armA,'ctot',t):10.1f}{at(armB,'ctot',t):10.1f}"
              f"{at(armA,'ccns',t):9.1f}{at(armB,'ccns',t):9.1f}"
              f"{at(armA,'ccns',t)/at(armB,'ccns',t):8.2f}")
    check("brain:plasma ratio jumps within 1 h of the switch",
          at(armA, 'ratio_meas', 9) / at(armB, 'ratio_meas', 9), 1.35, 3.0, "x")
    check("...while the MEASURED plasma level falls",
          at(armA, 'ctot', 9) / at(armB, 'ctot', 9), None, 0.95, "x",
          "the lab number improves as the patient is poisoned")
    check("brain salicylate nonetheless rises at 1 h",
          at(armA, 'ccns', 9) / at(armB, 'ccns', 9), 1.05, None, "x")
    check("brain salicylate excess at 24 h",
          at(armA, 'ccns', 24) / at(armB, 'ccns', 24), 1.15, None, "x",
          "renal clearance also collapses in acid urine")
    check("cumulative CNS injury at 36 h is higher in arm A",
          at(armA, 'CNSI', 36) / max(at(armB, 'CNSI', 36), 1e-6), 1.15, None, "x")
    sedo = simulate(dose_g=30.0, tmax=36, form='massive',
                    sched=Sched(rx_supportive() + [(8.0, dict(SED=0.45))]))
    check("sedation alone (45% respiratory depression) raises brain level",
          at(sedo, 'ccns', 12) / at(base, 'ccns', 12), 1.05, None, "x",
          "an opioid co-ingestion is the same lesion")

    # ----------------------------------------------------------------
    hdr("5. URINARY ALKALINISATION AND ITS POTASSIUM FAILURE MODE")
    # ----------------------------------------------------------------
    N = simulate(dose_g=30.0, tmax=48, form='massive', sched=Sched(rx_supportive()))
    S = simulate(dose_g=30.0, tmax=48, form='massive',
                 sched=Sched([(4.0, dict(R_FLUID=0.25, R_KCL=10.0))]))
    B = simulate(dose_g=30.0, tmax=48, form='massive', sched=Sched(rx_alkalinise()))
    Bk = simulate(dose_g=30.0, tmax=48, form='massive',
                  sched=Sched(rx_alkalinise(kcl=0.0)))
    print(f"   {'arm':<24s}{'uPH12':>7s}{'CLren':>9s}{'plasma24':>10s}"
          f"{'brain24':>9s}{'K24':>6s}{'CNSI48':>8s}")
    for nm, X in (("supportive only", N), ("saline diuresis", S),
                  ("NaHCO3 + K", B), ("NaHCO3, no K", Bk)):
        print(f"   {nm:<24s}{at(X,'uph',12):7.2f}"
              f"{1000*at(X,'cl_renal',12)/60:9.1f}{at(X,'ctot',24):10.1f}"
              f"{at(X,'ccns',24):9.1f}{at(X,'kp',24):6.2f}{at(X,'CNSI',48):8.2f}")
    check("alkalinisation raises urine pH above 7.5", at(B, 'uph', 12), 7.5, 8.2)
    check("renal clearance gain vs supportive care",
          at(B, 'cl_renal', 12) / max(at(N, 'cl_renal', 12), 1e-9), 4.0, 60.0, "x",
          "lit. 10-20x")
    check("NaHCO3 beats equal-volume saline diuresis",
          at(B, 'cl_renal', 12) / max(at(S, 'cl_renal', 12), 1e-9), 2.0, 40.0, "x",
          "Prescott 1982")
    print("   the potassium failure mode, hour by hour (no-K arm):")
    print(f"   {'t (h)':>6s}{'K':>7s}{'HCO3':>7s}{'pH':>7s}{'gate':>7s}"
          f"{'uPH':>7s}{'CLren':>8s}")
    for t in (6.0, 12.0, 18.0, 24.0, 36.0, 48.0):
        print(f"   {t:6.0f}{at(Bk,'kp',t):7.2f}{at(Bk,'hco3',t):7.1f}"
              f"{at(Bk,'ph',t):7.2f}{at(Bk,'bic_gate',t):7.2f}"
              f"{at(Bk,'uph',t):7.2f}{1000*at(Bk,'cl_renal',t)/60:8.1f}")
    def t_cross(X, key, thr, tmax=48.0):
        m = X['time'] <= tmax
        w = np.where(X[key][m] >= thr)[0]
        return float(X['time'][m][w[0]]) if len(w) else float('nan')
    tk = t_cross(B, 'uph', 7.5)
    tnk = t_cross(Bk, 'uph', 7.5)
    info("hours to reach urine pH 7.5, with K", tk, "h")
    info("hours to reach urine pH 7.5, without K", tnk, "h")
    check("hypokalaemia delays urinary alkalinisation", tnk - tk, 4.0, None, "h",
          "it is a delay, not an abolition - and it costs bicarbonate")
    check("plasma HCO3 needed to alkalinise the urine without K",
          at(Bk, 'hco3', tnk) - at(B, 'hco3', tk), 2.0, None, "mmol/L",
          "the price of not replacing potassium")
    check("plasma K without replacement", at(Bk, 'kp', 24), 1.5, 3.4, "mmol/L")
    check("renal clearance lost at 12 h by not replacing K",
          at(B, 'cl_renal', 12) / max(at(Bk, 'cl_renal', 12), 1e-9), 1.5, None, "x")
    check("salicylate excreted by 24 h, with-K / no-K",
          at(B, 'AUR', 24) / max(at(Bk, 'AUR', 24), 1e-9), 1.05, None, "x",
          "modest in mass terms - but see the brain concentration below")
    check("brain salicylate at 48 h: no-K / with-K",
          at(Bk, 'ccns', 48) / at(B, 'ccns', 48), 1.05, 12.0, "x",
          "the therapy is only as good as the potassium")
    check("NaHCO3 lowers brain level before the urine responds",
          at(N, 'ccns', 4.6) / at(B, 'ccns', 4.6), 1.02, 2.0, "x",
          "plasma-pH arm of the same drug")

    # ----------------------------------------------------------------
    hdr("6. HAEMODIALYSIS - and the cost of waiting")
    # ----------------------------------------------------------------
    HE = simulate(dose_g=30.0, tmax=48, form='massive',
                  sched=Sched(rx_alkalinise() + [(6.0, dict(CL_HD=6.0, HD_BIC=45.0))]))
    HL = simulate(dose_g=30.0, tmax=48, form='massive',
                  sched=Sched(rx_alkalinise() + [(18.0, dict(CL_HD=6.0, HD_BIC=45.0))]))
    check("HD halves the plasma level within ~4 h",
          at(HE, 'ctot', 6) / at(HE, 'ctot', 10), 1.6, 6.0, "x")
    check("early HD gives less cumulative CNS injury than late HD",
          at(HL, 'CNSI', 48) / max(at(HE, 'CNSI', 48), 1e-6), 1.2, None, "x")
    check("HD clearance exceeds maximal alkaline-urine clearance",
          6.0 / max(at(B, 'cl_renal', 12), 1e-9), 1.5, None, "x",
          "why EXTRIP recommends HD in severe cases")
    info("cumulative CNS injury, HD at 6 h", at(HE, 'CNSI', 48))
    info("cumulative CNS injury, HD at 18 h", at(HL, 'CNSI', 48))
    info("salicylate removed by the dialyser (early)", at(HE, 'AHD', 48) / 1000, "g")

    # ----------------------------------------------------------------
    hdr("7. CHRONIC SALICYLISM - toxic at a 'reassuring' number")
    # ----------------------------------------------------------------
    pc = dict(P0)
    pc['GFR0'] = 3.6          # 60 mL/min, age-related
    pc['ALB'] = 28.0          # hypoalbuminaemia
    pc['AGE_F'] = 1.0
    chronic = [(t, 1.3) for t in np.arange(0, 336, 6)]       # 5.2 g/day, 14 d
    thirst = [(168.0, dict(R_INTAKE=0.055))]                 # stops drinking d7
    C = simulate(doses=chronic, tmax=340, p=pc, sched=Sched(thirst))
    print(f"   day 7  plasma {at(C,'ctot',168):5.0f} mg/L  brain {at(C,'ccns',168):5.1f}"
          f"  pH {at(C,'ph',168):.2f}  CNS {at(C,'cns_score',168):4.0f}")
    print(f"   day 14 plasma {at(C,'ctot',336):5.0f} mg/L  brain {at(C,'ccns',336):5.1f}"
          f"  pH {at(C,'ph',336):.2f}  CNS {at(C,'cns_score',336):4.0f}")
    c_chr = at(C, 'ctot', 336)
    # find the acute trajectory point with the same measured plasma level
    m = A['time'] < 30
    j = int(np.argmin(np.abs(A['ctot'][m] - c_chr)))
    print(f"   ACUTE patient at the same measured level "
          f"({A['ctot'][j]:.0f} mg/L, t = {A['time'][j]:.1f} h): "
          f"brain {A['ccns'][j]:.1f} mg/L, CNS {A['cns_score'][j]:.0f}")
    check("chronic brain:plasma ratio exceeds acute at the same level",
          at(C, 'ratio_meas', 336) / max(A['ratio_meas'][j], 1e-9), 1.3, 60.0, "x",
          "equilibration + low albumin + acidaemia")
    check("chronic plasma level sits below the classic 'severe' cut-off",
          c_chr, 150.0, 700.0, "mg/L", "lit. chronic toxicity at 150-450")
    check("but the chronic CNS score is high", at(C, 'cns_score', 336), 25.0, None,
          "/100")
    info("chronic patient free salicylate", at(C, 'cfree', 336), "mg/L")
    info("chronic patient albumin", pc['ALB'], "g/L")
    info("chronic patient urine pH", at(C, 'uph', 336))

    # ----------------------------------------------------------------
    hdr("8. THE DONE NOMOGRAM, EVALUATED RATHER THAN ASSUMED")
    # ----------------------------------------------------------------
    print("Done (1960) grades severity from the plasma level and the hours since")
    print("ingestion.  Both numbers are read off the simulation and compared with")
    print("the brain concentration the model computes.")

    def done_zone(c_mgdl, t_h):
        t_h = min(max(t_h, 6.0), 60.0)
        sev = 130.0 * 0.5 ** ((t_h - 6.0) / 20.0)
        mod = 65.0 * 0.5 ** ((t_h - 6.0) / 20.0)
        mild = 32.5 * 0.5 ** ((t_h - 6.0) / 20.0)
        return ("severe" if c_mgdl >= sev else "moderate" if c_mgdl >= mod
                else "mild" if c_mgdl >= mild else "asymptomatic")

    def brain_zone(cb):
        return ("severe" if cb > 90 else "moderate" if cb > 55 else
                "mild" if cb > 25 else "asymptomatic")

    rows = [("acute 30 g, 6 h", A, 6.0, 6.0),
            ("acute 30 g, 24 h", A, 24.0, 24.0),
            ("intubated, PaCO2 40, 9 h", armA, 9.0, 9.0),
            ("ventilation matched, 9 h", armB, 9.0, 9.0),
            ("chronic, day 14", C, 336.0, 60.0),
            ("enteric-coated, 3 h", None, 3.0, 3.0)]
    pe = dict(P0); pe['F_CONC'] = 0.45
    E = simulate(dose_g=30.0, tmax=96, form='enteric', p=pe,
                 sched=Sched(rx_supportive()))
    mism = 0
    for nm, X, t, tnom in rows:
        X = E if X is None else X
        c = at(X, 'ctot', t) / 10.0
        cb = at(X, 'ccns', t)
        z, real = done_zone(c, tnom), brain_zone(cb)
        bad = z != real
        mism += bad
        print(f"   {nm:<27s}{c:7.1f} mg/dL  Done: {z:<12s} brain {cb:6.1f} "
              f"-> {real:<12s}{'   <-- MISMATCH' if bad else ''}")
    check("the nomogram misclassifies several real situations", mism, 2, None,
          "cases", "it reads the wrong variable")

    # ----------------------------------------------------------------
    hdr("9. ENTERIC-COATED ASPIRIN AND CONCRETION")
    # ----------------------------------------------------------------
    tpk = float(E['time'][int(np.argmax(E['ctot']))])
    check("time of peak plasma level", tpk, 8.0, 60.0, "h", "lit. peaks >12 h late")
    check("the 3 h level understates the eventual peak",
          peak(E, 'ctot') / at(E, 'ctot', 3.0), 2.0, 60.0, "x",
          "a single early level is uninterpretable")
    EA = simulate(dose_g=30.0, tmax=96, form='enteric', p=pe,
                  sched=Sched(rx_supportive()
                              + [(2.0, dict(AAC=50.0)), (6.0, dict(AAC=50.0)),
                                 (10.0, dict(AAC=50.0))]))
    check("multi-dose charcoal lowers the peak", peak(E, 'ctot') / peak(EA, 'ctot'),
          1.05, 4.0, "x")

    # ----------------------------------------------------------------
    hdr("10. MASS BALANCE AND NUMERICAL INTEGRITY")
    # ----------------------------------------------------------------
    Z = simulate(dose_g=30.0, tmax=96, form='massive',
                 sched=Sched(rx_alkalinise()), n=5000)
    given = 30000.0 * (MW_SAL / MW_ASA) * P0['F_ORAL']
    accounted = (Z['ACENT'][-1] + Z['APER'][-1] + Z['ACNS'][-1]
                 + Z['ASU'][-1] + Z['APG'][-1] + Z['AAG'][-1] + Z['AGA'][-1]
                 + Z['AUR'][-1] + Z['AHD'][-1]
                 + (Z['AST'][-1] + Z['AGUT'][-1] + Z['ACONC'][-1]) * P0['F_ORAL']
                 * MW_SAL / MW_ASA
                 + Z['AASA'][-1] * MW_SAL / MW_ASA)
    err = 100.0 * abs(accounted - given) / given
    check("salicylate mass-balance error", err, 0.0, 1.0, "%")
    check("no negative concentrations", float(np.min(Z['ctot'])), 0.0, None, "mg/L")
    allsc = (A, B, Bk, S, N, HE, HL, C, E, armA, armB, NOFL)
    check("pH stays physiological in every scenario",
          float(min(np.min(X['ph']) for X in allsc)), 6.60, None)
    check("brain pHi stays physiological in every scenario",
          float(min(np.min(X['phb']) for X in allsc)), 6.60, None)

    # ----------------------------------------------------------------
    hdr("11. THE WHOLE ARGUMENT IN ONE TABLE")
    # ----------------------------------------------------------------
    print("The measured number, and what it is actually worth:")
    print(f"   {'situation':<32s}{'plasma':>8s}{'fu':>7s}{'pH':>7s}{'pHi(br)':>9s}"
          f"{'brain':>8s}{'brain/plasma':>14s}")
    tab = [("acute 30 g, 6 h", A, 6.0),
           ("acute 30 g, 30 h", A, 30.0),
           ("intubated to PaCO2 40", armA, 9.0),
           ("chronic elderly, day 14", C, 336.0)]
    for nm, X, t in tab:
        print(f"   {nm:<32s}{at(X,'ctot',t):8.0f}{at(X,'fu',t):7.2f}"
              f"{at(X,'ph',t):7.2f}{at(X,'phb',t):9.2f}{at(X,'ccns',t):8.1f}"
              f"{at(X,'ratio_meas',t):14.3f}")
    ratios = [at(X, 'ratio_meas', t) for _, X, t in tab]
    check("spread of brain:plasma ratio across situations",
          max(ratios) / min(ratios), 1.5, 20.0, "x",
          "one plasma number, many different doses to the brain")

    # ----------------------------------------------------------------
    hdr("SUMMARY")
    # ----------------------------------------------------------------
    npass = sum(1 for _, g in _results if g)
    print(f"{npass} / {len(_results)} checks passed")
    for nm, g in _results:
        if not g:
            print(f"   FAILED: {nm}")
    return 0 if npass == len(_results) else 1


if __name__ == "__main__":
    sys.exit(main())
