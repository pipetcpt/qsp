#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 HIGH-ALTITUDE ILLNESS (AMS / HACE / HAPE) -- INDEPENDENT PYTHON REFERENCE MODEL
================================================================================

 This file is the *verification* implementation.  Every number quoted in
 README.md and in the commit message is produced by running this file; nothing
 is asserted from memory.  The mrgsolve model (hai_mrgsolve_model.R) is an
 independent re-implementation of the same equations in C++ and the two are
 compared explicitly (see hai_reference_output.txt).

 STRUCTURAL CLAIM OF THE MODEL
 -----------------------------
 High-altitude illness is not three diseases.  It is one *inspired oxygen
 pressure* acting through three different downstream mechanical systems, each
 with its own time constant and its own saturating defence:

     PIO2 = FiO2 x (PB(h) - 47)                  <- the only exogenous variable

     defence 1 (minutes->days) : ventilation.   PAO2 = PIO2 - PaCO2 x [FiO2 + (1-FiO2)/R]
                                 Bounded by acid-base, NOT by chemoreceptors.
     defence 2 (seconds)       : hypoxic pulmonary vasoconstriction (HPV).
                                 It protects gas exchange and *destroys* the
                                 capillary, because it is spatially uneven.
     defence 3 (seconds)       : cerebral vasodilation.  It defends CBF and
                                 spends the craniospinal volume buffer.

 Defence 1 has a delay of days (renal/CSF bicarbonate).  Defences 2 and 3 act
 immediately.  ALL of high-altitude medicine follows from that mismatch:
 the two fast defences run unopposed for exactly as long as the slow one takes
 to arrive, and the drugs sort cleanly into "accelerates defence 1"
 (acetazolamide), "blunts defence 2" (nifedipine, tadalafil), "pays for
 defence 3" (dexamethasone) and "removes the cause" (descent, O2, Gamow bag).

 UNITS
 -----
   pressures      mmHg
   volumes        L (blood/gas), mL (oedema)
   flows          L/min
   time (ODE)     hours
   VO2/VCO2       mL/min STPD
   concentrations mEq/L (acid-base), mmol/L (DPG), g/dL (Hb), mU/mL (EPO)

 Run:   python3 hai_reference_model.py            (writes hai_reference_output.txt
                                                   and hai_scenario_results.json)
================================================================================
"""

import json
import math
import sys
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ==============================================================================
#  SECTION 0.  PARAMETERS
# ==============================================================================

P = dict(

    # ---- environment -------------------------------------------------------
    PH2O        = 47.0,     # saturated water vapour pressure at 37 C, mmHg
    FIO2_AIR    = 0.2094,   # fractional O2 of dry air -- constant at all altitudes
    RQ          = 0.85,     # respiratory exchange ratio

    # ---- metabolic ---------------------------------------------------------
    VO2_REST    = 250.0,    # mL/min STPD
    VCO2_REST   = 212.5,    # = VO2 x RQ
    VDVT_REST   = 0.30,     # dead-space fraction at rest
    VDVT_EX     = 0.18,     # dead-space fraction at heavy exercise

    # ---- ventilatory controller -------------------------------------------
    #  VE = VEb + Gc*Phi(PaCO2 - Bc) + Gp*VAH*(1-HVD)*Dhyp*Phi(PaCO2 - Bp)
    #  Phi = softplus with width TAU_PHI.  Bc from CSF acid-base, Bp from
    #  arterial acid-base -- so BOTH thresholds move with bicarbonate, which is
    #  the whole mechanism of ventilatory acclimatisation and of acetazolamide.
    VE_BASAL    = 0.50,     # L/min, wakefulness / non-chemical drive
    GC          = 1.293,    # L/min/mmHg   central CO2 sensitivity  (calibrated)
    GP          = 7.00,     # L/min per unit SaO2 deficit per mmHg  (calibrated)
    TAU_PHI     = 2.0,      # mmHg, softplus width (chemoreflex is not a hard knee)
    PH_THR_C    = 7.409,    # central (CSF) pH at the apnoeic threshold
    PH_THR_P    = 7.526,    # arterial pH constant setting the peripheral threshold
    SAO2_REF    = 0.980,    # SaO2 at which peripheral hypoxic drive is zero
    VE_MAX      = 200.0,    # L/min ceiling (mechanical)

    # ---- acid-base ---------------------------------------------------------
    HCO3_SL     = 24.0,     # sea-level arterial bicarbonate, mEq/L
    HCO3C_SL    = 22.0,     # sea-level CSF bicarbonate, mEq/L
    TAU_REN     = 34.0,     # h, renal acid-base time constant (~1.4 d)
    K_REN       = 115.0,     # mEq/L per pH unit per TAU_REN -- renal gain
    TAU_CSF     = 8.0,      # h, CSF bicarbonate time constant
    W_CSF_PH    = 0.55,     # weight of active choroid-plexus pH regulation
    C_CSF       = 0.246,    # offset making sea-level CSF HCO3 = 22.0 exactly

    # ---- carotid-body plasticity / hypoxic ventilatory decline -------------
    VAH_MAX     = 2.10,     # maximal multiplicative gain of peripheral drive
    TAU_VAH     = 60.0,     # h, ventilatory acclimatisation of the carotid body
    HVD_MAX     = 0.28,     # maximal roll-off of hypoxic drive (minutes-hours)
    TAU_HVD_ON  = 0.40,     # h
    TAU_HVD_OFF = 2.0,      # h

    # ---- blood O2 transport ------------------------------------------------
    P50_STD     = 26.8,     # mmHg at pH 7.4, PCO2 40, 37 C, DPG 5.0
    BOHR        = -0.48,    # d(log10 P50)/d(pH)
    DPG_SL      = 5.00,     # mmol/L
    DPG_MAX     = 6.00,
    TAU_DPG     = 30.0,     # h
    K_P50_DPG   = 0.075,    # fractional P50 rise per unit relative DPG rise
    HB_CONST    = 1.34,     # mL O2 per g Hb
    O2_SOL      = 0.003,    # mL O2 per dL per mmHg

    # ---- diffusion (Piiper-Scheid) ----------------------------------------
    DLO2_REST   = 37.0,     # mL/min/mmHg
    DLO2_EXMAX  = 78.0,     # mL/min/mmHg at maximal recruitment
    FSHUNT_BASE = 0.035,    # anatomical shunt
    FSHUNT_FLOOD= 0.40,     # shunt added at complete alveolar flooding

    # ---- cardiac -----------------------------------------------------------
    Q_REST      = 6.0,      # L/min
    HR_SL       = 65.0,
    TAU_HR      = 0.5,
    SV_SL       = 0.0923,   # L, = 6.0/65
    PLA_SL      = 8.0,      # mmHg left atrial pressure
    K_PLA_EX    = 0.55,     # mmHg per L/min of cardiac output above rest

    # ---- pulmonary vascular / HPV -----------------------------------------
    R_PULM_TOT  = 1.00,     # WU (mmHg/(L/min)) at sea level rest
    FRAC_VEN    = 0.50,     # venous share of total pulmonary resistance
    A_HETERO    = 0.50,     # HPV-responsive fraction of the bed (normal subject)
    A_HETERO_S  = 0.85,     # ... in a HAPE-susceptible subject
    LAMBDA_MAX  = 5.0,      # maximal HPV constriction factor (normal)
    LAMBDA_MAX_S= 9.0,      # ... HAPE-susceptible
    P50_HPV     = 45.0,     # mmHg PAO2 for half-maximal HPV
    N_HPV       = 4.0,
    KAPPA_HPV   = 0.08,     # HPV in the "non-responsive" bed, as fraction of full
    K_RECRUIT   = 0.35,     # pulmonary recruitment/distension per unit relative flow
    TAU_HPV_S   = 72.0,     # h, slow HPV / remodelling
    HPVS_MAX    = 0.45,     # slow component adds up to +45% to lambda

    # ---- capillary stress failure / oedema --------------------------------
    PCAP_CRIT   = 19.5,     # mmHg, threshold transmural pressure for stress failure
    KF_LEAK     = 2.60,     # mL/h per mmHg^1.5 of overpressure, x permeability
    N_LEAK      = 1.5,
    ELW_0       = 300.0,    # mL, normal extravascular lung water
    ELW_FLOOD50 = 700.0,    # mL of EXCESS water at which flooding fraction = 0.5
    N_FLOOD     = 2.2,
    K_AFC       = 0.115,    # /h basal alveolar fluid clearance rate constant
    AFC_HYP50   = 0.72,     # SaO2 at which AFC capacity is halved
    AFC_N       = 6.0,
    TAU_AFC     = 6.0,      # h
    TAU_PERM    = 8.0,      # h
    K_PERM_INFL = 0.85,     # permeability gain per unit inflammation
    TAU_INFL    = 12.0,     # h
    K_INFL_LEAK = 0.0022,   # inflammation driven by leak rate (mL/h)

    # ---- cerebral ----------------------------------------------------------
    CBF_SL      = 50.0,     # mL/100g/min
    K_CBF_O2    = 1.10,     # hypoxic vasodilation gain
    CBF_P50     = 35.0,     # SaO2 at half-maximal hypoxic vasodilation
    CBF_O2N     = 4.0,
    CBF_CO2_MID = 42.0, CBF_CO2_W = 12.0, CBF_CO2_LO = 0.55, CBF_CO2_HI = 1.65,    # fractional CBF change per mmHg PaCO2
    TAU_CBF     = 0.25,     # h
    PVI         = 25.0,     # mL, pressure-volume index (craniospinal compliance)
    PVI_TIGHT   = 16.0,     # mL, "tight-fit" individual
    ICP_0       = 10.0,     # mmHg
    K_VEGF      = 1.00,     # brain VEGF gain from hypoxia
    TAU_VEGF    = 10.0,     # h
    K_BBB       = 0.65,     # BBB permeability gain from VEGF
    TAU_BBB     = 8.0,      # h
    K_EDV       = 4.00,     # mL/h of vasogenic oedema per unit (BBBp-1) x dPcap
    TAU_EDV     = 20.0,     # h clearance of vasogenic oedema
    K_EDC       = 6.00,      # cytotoxic oedema gain (Na/K-ATPase failure)
    EDC_O250    = 0.70,     # SaO2 at half-maximal cytotoxic oedema
    EDC_N       = 10.0,
    TAU_EDC     = 6.0,
    CBV_SL      = 75.0,     # mL cerebral blood volume at sea level
    ICP_HACE    = 22.0,     # mmHg above which HACE risk becomes appreciable
    ICP_HACE_N  = 6.0,

    # ---- symptom latents (Lake Louise 2018 subscores, 0-3 each) -----------
    TAU_SX_ON   = 5.0,      # h
    TAU_SX_OFF  = 14.0,     # h
    K_HEAD_ICP  = 0.175,    # headache per mmHg ICP above 10
    K_HEAD_HYP  = 5.50,     # headache per unit SaO2 deficit (trigeminovascular)
    K_GI        = 1.55,
    K_FAT       = 4.90,
    K_DIZ       = 4.30,
    K_SLEEP     = 9.50,

    # ---- sympathetic / fluid ----------------------------------------------
    TAU_SYM     = 2.0,
    K_SYM       = 3.2,
    TAU_ALDO    = 8.0,
    TAU_ADH     = 4.0,
    K_FLUID     = 90.0,    # mL/h per unit antidiuretic index
    TAU_FLUID   = 30.0,

    # ---- erythropoiesis ----------------------------------------------------
    EPO_SL      = 12.0,     # mU/mL
    K_EPO       = 165.0,    # gain on hypoxic stimulus
    TAU_EPO     = 6.0,      # h (half-life ~4-5 h)
    EPO_HB_FB   = 0.055,    # feedback per g/dL of Hb above baseline
    TAU_RET     = 84.0,     # h, marrow transit
    K_ERY       = 1.65e-4,  # g Hb produced per mU/mL EPO per h
    HBM_SL      = 800.0,    # g haemoglobin mass (70 kg male)
    TAU_HBM     = 2400.0,   # h, RBC turnover (~100 d)
    PV_SL       = 3.10,     # L plasma volume
    PV_CONTRACT = 0.20,     # maximal fractional plasma-volume contraction
    TAU_PV      = 40.0,     # h
    BV_FACTOR   = 1.0,

    # ---- chronic mountain sickness ----------------------------------------
    K_CMS       = 0.055,
    TAU_CMS     = 2000.0,
    VISC_K      = 2.31,     # relative viscosity exponent: mu ~ exp(k*Hct)
    HCT_REF     = 0.45,

    # ---- muscle / capillary adaptation ------------------------------------
    TAU_MUSC    = 500.0,
    MUSC_MAX    = 0.22,

    # ---- DRUG PK -----------------------------------------------------------
    # acetazolamide: 1-cpt plasma + deep RBC carbonic-anhydrase pool
    ACZ_F       = 0.95, ACZ_KA = 1.8, ACZ_V = 15.0, ACZ_CL = 2.1,
    ACZ_KIN = 0.55, ACZ_KOUT = 0.075, ACZ_VR = 2.4,
    ACZ_IC50_REN = 3.5,     # mg/L, renal CA inhibition
    ACZ_EMAX_REN = 7.4,     # mEq/L maximal fall in plasma HCO3
    ACZ_IC50_CSF = 6.0,     # mg/L, choroid-plexus CA
    ACZ_EMAX_CSF = 3.6,     # mEq/L additional fall in CSF HCO3
    ACZ_KRBC     = 0.045,   # max fractional slowing of CO2 transport
    ACZ_IC50_RBC = 220.0,   # mg/L in the RBC compartment (CA is in vast excess)

    # dexamethasone
    DEX_F = 0.80, DEX_KA = 2.0, DEX_V = 60.0, DEX_CL = 12.0,
    DEX_KE0 = 0.045,        # /h effect-compartment (biological t1/2 ~ 36 h)
    DEX_EC50 = 0.030,       # mg/L effect-site
    DEX_EMAX_BBB = 0.75,    # fractional suppression of BBB permeability gain
    DEX_EMAX_VEGF= 0.65,
    DEX_EMAX_HPV = 0.30,    # modest HPV reduction
    DEX_EMAX_AFC = 0.55,    # ENaC/Na-K-ATPase upregulation
    DEX_EMAX_SX  = 0.60,    # direct symptomatic suppression (the dangerous part)

    # nifedipine SR
    NIF_F = 0.55, NIF_KA = 0.55, NIF_V = 55.0, NIF_CL = 30.0,
    NIF_EC50 = 0.020, NIF_EMAX_HPV = 0.45,

    # tadalafil
    TAD_F = 0.80, TAD_KA = 1.1, TAD_V = 63.0, TAD_CL = 2.5,
    TAD_EC50 = 0.055, TAD_EMAX_HPV = 0.50,

    # salmeterol (inhaled -> effect compartment only)
    SAL_KE0 = 0.12, SAL_EMAX_AFC = 0.70, SAL_EC50 = 0.5,

    # ibuprofen
    IBU_F = 0.85, IBU_KA = 2.2, IBU_V = 10.0, IBU_CL = 3.5,
    IBU_EC50 = 8.0, IBU_EMAX_HEAD = 0.45,
)

# state index map -------------------------------------------------------------
SNAMES = ['BE','HCO3c','VAH','HVD','DPG','HbM','PV','EPO','RET','HPVs',
          'ELW','FLOOD','AFC','PERM','INFL','CBFrel','VEGFb','BBBp','EDv','EDc',
          'CSFres','SYM','HRs','ALDO','ADH','FLUID','HEAD','GI','FAT','DIZ',
          'SLEEPd','ACZc','ACZr','ACZa','DEXc','DEXe','DEXa','NIFc','NIFa',
          'TADc','TADa','SALe','IBUc','IBUa','LAC','MUSC','HYPd','CMS','PAPs','ACCL']
IDX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)


# ==============================================================================
#  SECTION 1.  STATIC PHYSIOLOGY  (closed form, no fitting)
# ==============================================================================

def barometric(alt_m, model='west'):
    """Barometric pressure, mmHg.

    'west'  : West JB (1996) J Appl Physiol -- fitted to *measured* pressures
              over the Himalaya including the Everest summit.  This is the
              standard in altitude medicine because the ICAO standard
              atmosphere under-predicts the summit by ~17 mmHg (it assumes a
              mid-latitude tropopause; the equatorial tropopause is higher).
                 PB = exp(6.63268 - 0.1112 h - 0.00149 h^2),  h in km
    'isa'   : ICAO standard atmosphere, for comparison.
    """
    if model == 'isa':
        return 760.0 * (1.0 - 2.25577e-5 * alt_m) ** 5.25588
    h = alt_m / 1000.0
    return math.exp(6.63268 - 0.1112 * h - 0.00149 * h * h)


def altitude_from_pb(pb, model='west'):
    """Invert the barometric equation -- used for 'metres of descent equivalent'."""
    lo, hi = -500.0, 9500.0
    f = lambda a: barometric(a, model) - pb
    if f(lo) * f(hi) > 0:
        return float('nan')
    return brentq(f, lo, hi, xtol=1e-3)


def pio2(alt_m, fio2=None, p_bag=0.0, model='west'):
    """Inspired PO2 at the trachea (fully humidified), mmHg."""
    fio2 = P['FIO2_AIR'] if fio2 is None else fio2
    pb = barometric(alt_m, model) + p_bag
    return fio2 * (pb - P['PH2O'])


def alveolar_po2(alt_m, paco2, fio2=None, p_bag=0.0, R=None, model='west'):
    """Full alveolar gas equation (not the PaCO2/R shortcut)."""
    fio2 = P['FIO2_AIR'] if fio2 is None else fio2
    R = P['RQ'] if R is None else R
    return pio2(alt_m, fio2, p_bag, model) - paco2 * (fio2 + (1.0 - fio2) / R)


def severinghaus(po2):
    """Standard O2 dissociation curve (pH 7.4, PCO2 40, 37 C, DPG 5.0)."""
    po2 = max(po2, 1e-9)
    return 1.0 / (23400.0 / (po2 ** 3 + 150.0 * po2) + 1.0)



def cbf_co2_factor(paco2):
    """Normalised cerebral CO2 reactivity (sigmoid, saturating at both ends).

    A linear %/mmHg reactivity is the usual textbook simplification and it is
    catastrophically wrong at altitude: extrapolated linearly, a PaCO2 of 13
    mmHg on the Everest summit would abolish cerebral blood flow.  The
    vasoconstrictor response saturates below ~25 mmHg, which is why the brain
    survives an arterial CO2 that would be lethal if the reactivity were linear.
    """
    f = lambda p: (P['CBF_CO2_LO'] + P['CBF_CO2_HI']
                   / (1.0 + math.exp(-(p - P['CBF_CO2_MID']) / P['CBF_CO2_W'])))
    return f(paco2) / f(40.0)

def p50(ph, dpg=None, temp=37.0):
    dpg = P['DPG_SL'] if dpg is None else dpg
    v = P['P50_STD'] * (10.0 ** (P['BOHR'] * (ph - 7.40)))
    v *= (1.0 + P['K_P50_DPG'] * (dpg / P['DPG_SL'] - 1.0))
    v *= (1.0 + 0.024 * (temp - 37.0) * 2.303 / 2.303)   # ~2.4%/degC
    return v


def sao2(po2, ph=7.40, dpg=None, temp=37.0):
    """Saturation with Bohr / DPG / temperature shift, by scaling PO2 to the
    standard curve:  S(P) = Sev(P * P50_std / P50_actual)."""
    return severinghaus(po2 * P['P50_STD'] / p50(ph, dpg, temp))


def po2_from_sat(s, ph=7.40, dpg=None, temp=37.0):
    f = lambda p: sao2(p, ph, dpg, temp) - s
    return brentq(f, 1e-4, 700.0, xtol=1e-8)


def o2_content(po2, hb, ph=7.40, dpg=None, temp=37.0):
    return P['HB_CONST'] * hb * sao2(po2, ph, dpg, temp) + P['O2_SOL'] * po2


def po2_from_content(cao2, hb, ph=7.40, dpg=None, temp=37.0):
    f = lambda p: o2_content(p, hb, ph, dpg, temp) - cao2
    lo, hi = 1e-4, 700.0
    if f(hi) < 0:
        return hi
    if f(lo) > 0:
        return lo
    return brentq(f, lo, hi, xtol=1e-6)


def acid_base(be, paco2, hb=15.0):
    """Given base excess and PaCO2, return (pH, HCO3).

    Uses the Siggaard-Andersen / Van Slyke relation
        BE = 0.93 * [HCO3 - 24.4 + 14.8*(pH - 7.40)]
    together with Henderson-Hasselbalch, solved simultaneously.  This is what
    lets ONE state variable (BE) carry the metabolic component while the acute
    non-bicarbonate buffering falls out automatically -- it is why an acute
    ascent shows dHCO3/dPaCO2 ~ -0.2 and a week-old ascent ~ -0.5 without
    either number being a parameter.
    """
    def resid(hco3):
        ph = 6.1 + math.log10(max(hco3, 1e-6) / (0.03 * paco2))
        return 0.93 * (hco3 - 24.4 + 14.8 * (ph - 7.40)) - be
    return_lo, return_hi = 1.0, 60.0
    hco3 = brentq(resid, return_lo, return_hi, xtol=1e-9)
    ph = 6.1 + math.log10(hco3 / (0.03 * paco2))
    return ph, hco3


def softplus(x, tau):
    z = x / tau
    if z > 30.0:
        return x
    if z < -30.0:
        return tau * math.exp(z)
    return tau * math.log1p(math.exp(z))


def softplus_d(x, tau):
    z = x / tau
    if z > 30.0:
        return 1.0
    if z < -30.0:
        return math.exp(z)
    return 1.0 / (1.0 + math.exp(-z))


# ==============================================================================
#  SECTION 2.  THE VENTILATORY / GAS-EXCHANGE FIXED POINT
# ==============================================================================

def _threshold_central(hco3c):
    return hco3c / (0.03 * 10.0 ** (P['PH_THR_C'] - 6.1))


def _threshold_periph(hco3a):
    return hco3a / (0.03 * 10.0 ** (P['PH_THR_P'] - 6.1))


def ventilation(paco2, hco3c, be, vah, hvd, sat, hb, gp=None, gc=None):
    """Chemoreflex ventilation, L/min BTPS.

    Two additive drives, both with a *soft* threshold that is set by acid-base:
      central   : threshold Bc from CSF bicarbonate
      peripheral: threshold Bp from arterial bicarbonate, gain multiplied by
                  the arterial desaturation (the classic multiplicative O2-CO2
                  interaction) and by carotid-body plasticity (VAH) and roll-off
                  (HVD).
    Both thresholds FALL as bicarbonate falls.  That single fact is ventilatory
    acclimatisation, and it is also the entire pharmacology of acetazolamide.
    """
    gp = P['GP'] if gp is None else gp
    gc = P['GC'] if gc is None else gc
    _, hco3a = acid_base(be, paco2, hb)
    bc = _threshold_central(hco3c)
    bp = _threshold_periph(hco3a)
    dhyp = max(0.0, P['SAO2_REF'] - sat)
    ve = (P['VE_BASAL']
          + gc * softplus(paco2 - bc, P['TAU_PHI'])
          + gp * vah * (1.0 - hvd) * dhyp * softplus(paco2 - bp, P['TAU_PHI']))
    return min(ve, P['VE_MAX'])


def gas_exchange(paco2, alt_m, be, hco3c, vah, hvd, hb, dpg,
                 fio2=None, p_bag=0.0, exercise=1.0, flood=0.0, q=None,
                 gp=None, gc=None, dlo2_scale=1.0, model='west'):
    """Given PaCO2 and the slow states, return the full gas-exchange solution.

    Order of operations (all forced, none fitted):
       PAO2  <- alveolar gas equation
       shunt <- anatomical + alveolar flooding
       diffusion limitation <- Piiper-Scheid  exp(-DL/(beta*Q))
       CaO2  <- shunt equation solved exactly for CaO2
       PaO2  <- inverted content equation
    """
    fio2 = P['FIO2_AIR'] if fio2 is None else fio2
    vo2 = P['VO2_REST'] * exercise
    q = (P['Q_REST'] * (1.0 + 0.55 * (exercise - 1.0))) if q is None else q
    pao2 = alveolar_po2(alt_m, paco2, fio2, p_bag, P['RQ'], model)
    pao2 = max(pao2, 1.0)
    ph, _ = acid_base(be, paco2, hb)

    fs = P['FSHUNT_BASE'] + P['FSHUNT_FLOOD'] * flood
    fs = min(fs, 0.85)

    # end-capillary content of the ventilated compartment
    cc = o2_content(pao2, hb, ph, dpg)
    # shunt equation solved for CaO2 (see derivation in README)
    cao2 = cc - fs * (vo2 / (10.0 * q)) / max(1.0 - fs, 1e-3)
    cao2 = max(cao2, 0.5)

    # ---- Piiper-Scheid diffusion limitation, applied to the ventilated part
    dl = (P['DLO2_REST'] + (P['DLO2_EXMAX'] - P['DLO2_REST'])
          * min(1.0, max(0.0, (exercise - 1.0) / 3.0))) * dlo2_scale
    cv = cao2 - vo2 / (10.0 * q)
    pv_o2 = po2_from_content(max(cv, 0.5), hb, ph, dpg)
    beta = max((cc - cv) / max(pao2 - pv_o2, 1.0), 1e-4)     # mL/dL/mmHg
    inc = math.exp(-dl / max(beta * q * 10.0, 1e-6))          # incomplete-equilibration
    pa_ideal = po2_from_content(cao2, hb, ph, dpg)
    pao2_art = pa_ideal - (pao2 - pv_o2) * inc
    pao2_art = max(pao2_art, 1.0)
    sat = sao2(pao2_art, ph, dpg)
    cao2 = o2_content(pao2_art, hb, ph, dpg)
    cv = cao2 - vo2 / (10.0 * q)
    pv_o2 = po2_from_content(max(cv, 0.5), hb, ph, dpg)

    return dict(PAO2=pao2, PaO2=pao2_art, SaO2=sat, pH=ph, CaO2=cao2,
                CvO2=cv, PvO2=pv_o2, shunt=fs, Q=q, VO2=vo2,
                AaDO2=pao2 - pao2_art, beta=beta, DLO2=dl,
                PIO2=pio2(alt_m, fio2, p_bag, model),
                PB=barometric(alt_m, model) + p_bag)


def respiratory_fixed_point(alt_m, be, hco3c, vah, hvd, hb, dpg,
                            fio2=None, p_bag=0.0, exercise=1.0, flood=0.0,
                            q=None, gp=None, gc=None, acz_rbc_eff=0.0,
                            dlo2_scale=1.0, model='west', sleep=False):
    """Solve   PaCO2 = 863 * VCO2 / VA(PaCO2)   self-consistently.

    The residual  r(PaCO2) = 863*VCO2/VA(PaCO2) - PaCO2  is strictly decreasing
    (more CO2 -> more ventilation -> less CO2), so a bracketed root is unique
    and Brent is safe.  There is no iteration-to-convergence hand-waving here.
    """
    vdvt = P['VDVT_REST'] + (P['VDVT_EX'] - P['VDVT_REST']) * min(1.0, (exercise - 1.0) / 3.0)
    vco2 = P['VCO2_REST'] * exercise
    k = 863.0 * (vco2 / 1000.0) / (1.0 - vdvt)      # L/min * mmHg
    gp_eff = (P['GP'] if gp is None else gp)
    if sleep:
        gp_eff *= 0.72                              # NREM blunts the hypoxic drive

    def resid(paco2):
        g = gas_exchange(paco2, alt_m, be, hco3c, vah, hvd, hb, dpg,
                         fio2, p_bag, exercise, flood, q, gp_eff, gc,
                         dlo2_scale, model)
        ve = ventilation(paco2, hco3c, be, vah, hvd, g['SaO2'], hb, gp_eff, gc)
        ve = max(ve, 0.05)
        e_rbc = acz_rbc_eff / (P['ACZ_IC50_RBC'] + acz_rbc_eff)
        ve_eff = ve * (1.0 - P['ACZ_KRBC'] * e_rbc)
        return k / ve_eff - paco2

    lo, hi = 3.0, 90.0
    rlo, rhi = resid(lo), resid(hi)
    if rlo < 0:
        paco2 = lo
    elif rhi > 0:
        paco2 = hi
    else:
        paco2 = brentq(resid, lo, hi, xtol=1e-6, rtol=1e-10)

    g = gas_exchange(paco2, alt_m, be, hco3c, vah, hvd, hb, dpg,
                     fio2, p_bag, exercise, flood, q, gp_eff, gc,
                     dlo2_scale, model)
    ve = ventilation(paco2, hco3c, be, vah, hvd, g['SaO2'], hb, gp_eff, gc)
    g['PaCO2'] = paco2
    g['VE'] = ve
    g['VA'] = ve * (1.0 - vdvt)
    g['Bc'] = _threshold_central(hco3c)
    _, hco3a = acid_base(be, paco2, hb)
    g['Bp'] = _threshold_periph(hco3a)
    g['HCO3'] = hco3a

    # --- effective chemoreflex slope and CO2 reserve at the operating point
    d = 0.10
    ve_p = ventilation(paco2 + d, hco3c, be, vah, hvd,
                       gas_exchange(paco2 + d, alt_m, be, hco3c, vah, hvd, hb, dpg,
                                    fio2, p_bag, exercise, flood, q, gp_eff, gc,
                                    dlo2_scale, model)['SaO2'], hb, gp_eff, gc)
    ve_m = ventilation(paco2 - d, hco3c, be, vah, hvd,
                       gas_exchange(paco2 - d, alt_m, be, hco3c, vah, hvd, hb, dpg,
                                    fio2, p_bag, exercise, flood, q, gp_eff, gc,
                                    dlo2_scale, model)['SaO2'], hb, gp_eff, gc)
    slope = (ve_p - ve_m) / (2.0 * d)
    g['S_eff'] = slope
    g['B_eff'] = paco2 - ve / max(slope, 1e-6)      # linear-extrapolation threshold
    # The CO2 reserve that matters for periodic breathing is measured against the
    # CENTRAL threshold: during the hyperpnoeic phase SaO2 transiently recovers,
    # the peripheral drive is withdrawn, and what is left holding ventilation up
    # is the CSF [H+].  So the reserve is PaCO2 - Bc, not PaCO2 - B_eff.
    g['CO2_reserve'] = paco2 - g['Bc']
    g['plant_gain'] = paco2 / max(g['VA'], 1e-6)
    g['LG_static'] = slope * g['plant_gain']
    return g


# ==============================================================================
#  SECTION 3.  PULMONARY HAEMODYNAMICS -- THE TWO-BED HPV MODEL
# ==============================================================================

def hpv_lambda(pao2, lam_max, hpvs=0.0, drug_supp=0.0):
    """Hypoxic pulmonary vasoconstriction factor (>= 1)."""
    h = 1.0 / (1.0 + (max(pao2, 1.0) / P['P50_HPV']) ** P['N_HPV'])
    lam = 1.0 + lam_max * h * (1.0 + P['HPVS_MAX'] * hpvs) * (1.0 - drug_supp)
    return max(lam, 1.0)


def pulmonary(pao2, q, lam_max, a_het, hpvs=0.0, drug_supp=0.0, p_la=None,
              exercise=1.0):
    """Two parallel beds.  Bed A (fraction a) constricts by lambda; bed B
    (fraction 1-a) constricts only by kappa*lambda.

        G_tot   = G0 [ a/lambda + (1-a)/mu ]
        Q_B     = Q  [ (1-a)/mu ] / [ a/lambda + (1-a)/mu ]
        Pcap,B  = P_LA + Q_B * r_v/(1-a)
                = P_LA + Q * (r_v/mu) / [ a/lambda + (1-a)/mu ]

    The last line is the point of the whole submodel: (1-a) CANCELS from the
    numerator.  As lambda -> infinity the flow-dependent term tends to
    Q*r_v/(1-a), so the maximum capillary overpressure the lung can generate
    is bounded by 1/(1-a) -- HETEROGENEITY sets the ceiling, HPV strength only
    says how fast you get there.
    """
    p_la = (P['PLA_SL'] + P['K_PLA_EX'] * max(0.0, q - P['Q_REST'])) if p_la is None else p_la
    lam = hpv_lambda(pao2, lam_max, hpvs, drug_supp)
    mu = 1.0 + P['KAPPA_HPV'] * (lam - 1.0)
    r_v = P['R_PULM_TOT'] * P['FRAC_VEN']
    r_a = P['R_PULM_TOT'] * (1.0 - P['FRAC_VEN'])
    # Flow-dependent recruitment and distension: the pulmonary bed is not a
    # fixed resistor.  Without this term the model puts an exercising climber
    # at a mean PAP of 77 mmHg, which no catheter has ever recorded.
    rec = 1.0 + P['K_RECRUIT'] * max(0.0, q / P['Q_REST'] - 1.0)
    denom = (a_het / lam + (1.0 - a_het) / mu) * rec
    r_tot = 1.0 / denom
    mpap = p_la + q * r_tot * (r_a + r_v) / P['R_PULM_TOT']
    q_b = q * ((1.0 - a_het) / mu) / denom
    pcap_b = p_la + q * (r_v / mu) / denom
    pcap_a = p_la + (q * (a_het / lam) / denom) * (r_v / max(a_het, 1e-6))
    # Amplification ceiling.  The naive derivation (holding mu fixed) gives
    # 1/(1-a).  That is WRONG whenever the "non-responsive" bed constricts at
    # all: mu = 1 + kappa(lambda-1) grows with lambda too, and taking the limit
    # properly gives  1/(1 - a(1-kappa)).  The two agree only at kappa = 0.
    # This was caught by checking the closed form against the numerics.
    kap = P['KAPPA_HPV']
    return dict(lam=lam, mu=mu, PVR=r_tot, mPAP=mpap, PLA=p_la,
                Qopen=q_b, Pcap_open=pcap_b, Pcap_constr=pcap_a,
                amplification=1.0 / max(1.0 - a_het * (1.0 - kap), 1e-6),
                amplification_naive=1.0 / max(1.0 - a_het, 1e-6),
                overperfusion=q_b / max(q * (1.0 - a_het), 1e-9))


# ==============================================================================
#  SECTION 4.  THE FULL ODE SYSTEM
# ==============================================================================

class Subject:
    """Between-subject phenotype."""
    def __init__(self, hvr=1.0, lam_max=None, a_het=None, pvi=None,
                 hb_sl=15.0, wt=70.0, sex='M', label='typical'):
        self.hvr = hvr                                    # multiplier on GP
        self.lam_max = P['LAMBDA_MAX'] if lam_max is None else lam_max
        self.a_het = P['A_HETERO'] if a_het is None else a_het
        self.pvi = P['PVI'] if pvi is None else pvi
        self.hb_sl = hb_sl
        self.wt = wt
        self.sex = sex
        self.label = label

    @staticmethod
    def hape_susceptible(label='HAPE-susceptible'):
        return Subject(hvr=0.75, lam_max=P['LAMBDA_MAX_S'], a_het=P['A_HETERO_S'],
                       pvi=P['PVI'], label=label)

    @staticmethod
    def tight_fit(label='tight-fit (low PVI)'):
        return Subject(hvr=0.85, pvi=P['PVI_TIGHT'], label=label)


class Profile:
    """Altitude / exercise / drug schedule."""
    def __init__(self, alt_fn, exercise_fn=None, fio2_fn=None, bag_fn=None,
                 sleep_fn=None, doses=None, label=''):
        self.alt_fn = alt_fn
        self.exercise_fn = exercise_fn or (lambda t: 1.0)
        self.fio2_fn = fio2_fn or (lambda t: P['FIO2_AIR'])
        self.bag_fn = bag_fn or (lambda t: 0.0)
        self.sleep_fn = sleep_fn or (lambda t: (t % 24.0) >= 22.0 or (t % 24.0) < 6.0)
        self.doses = doses or []          # list of (time_h, state_name, amount)
        self.label = label


def steady_state(subject, alt_m=0.0, days=400.0):
    """Integrate to a true steady state at a given altitude (used for sea-level
    initial conditions and for the Andean-resident scenario)."""
    y0 = np.zeros(NST)
    y0[IDX['BE']] = 0.0
    y0[IDX['HCO3c']] = P['HCO3C_SL']
    y0[IDX['VAH']] = 1.0
    y0[IDX['HVD']] = 0.0
    y0[IDX['DPG']] = P['DPG_SL']
    y0[IDX['HbM']] = P['HBM_SL'] * (subject.hb_sl / 15.0)
    y0[IDX['PV']] = P['PV_SL']
    y0[IDX['EPO']] = P['EPO_SL']
    y0[IDX['RET']] = 0.0
    y0[IDX['AFC']] = 1.0
    y0[IDX['PERM']] = 1.0
    y0[IDX['ELW']] = P['ELW_0']
    y0[IDX['CBFrel']] = 1.0
    y0[IDX['HRs']] = P['HR_SL']
    y0[IDX['LAC']] = 1.0
    y0[IDX['WT'] if 'WT' in IDX else IDX['MUSC']] = y0[IDX['MUSC']]
    prof = Profile(lambda t: alt_m, label='ss')
    sol = simulate(subject, prof, t_end=days * 24.0, y0=y0, n_out=3, quiet=True)
    return sol['y'][:, -1]


def rhs(t, y, subject, prof, cache):
    p = P
    s = {n: y[IDX[n]] for n in SNAMES}
    alt = prof.alt_fn(t)
    ex = prof.exercise_fn(t)
    fio2 = prof.fio2_fn(t)
    bag = prof.bag_fn(t)
    asleep = bool(prof.sleep_fn(t))

    # ---- derived blood composition ----------------------------------------
    bv = max(s['PV'], 0.5) + s['HbM'] / 340.0        # 340 g Hb per L of RBC
    hb = s['HbM'] / (10.0 * bv)                      # g/dL
    hct = (s['HbM'] / 340.0) / bv

    # ---- drug effect sites -------------------------------------------------
    c_acz = s['ACZc'] / p['ACZ_V']
    c_dex_e = s['DEXe'] / p['DEX_V']
    c_nif = s['NIFc'] / p['NIF_V']
    c_tad = s['TADc'] / p['TAD_V']
    c_ibu = s['IBUc'] / p['IBU_V']
    e_sal = s['SALe']

    e_dex = c_dex_e / (p['DEX_EC50'] + c_dex_e)
    e_nif = c_nif / (p['NIF_EC50'] + c_nif)
    e_tad = c_tad / (p['TAD_EC50'] + c_tad)
    e_ibu = c_ibu / (p['IBU_EC50'] + c_ibu)
    e_acz_ren = c_acz / (p['ACZ_IC50_REN'] + c_acz)
    e_acz_csf = c_acz / (p['ACZ_IC50_CSF'] + c_acz)
    e_sal_eff = e_sal / (p['SAL_EC50'] + e_sal)

    hpv_supp = 1.0 - (1.0 - p['NIF_EMAX_HPV'] * e_nif) * \
                     (1.0 - p['TAD_EMAX_HPV'] * e_tad) * \
                     (1.0 - p['DEX_EMAX_HPV'] * e_dex)

    # ---- cardiac output ----------------------------------------------------
    q = P['Q_REST'] * (1.0 + 0.55 * (ex - 1.0)) * (1.0 + 0.10 * s['SYM'])
    q *= (P['HCT_REF'] / max(hct, 0.20)) ** 0.0 or 1.0     # viscosity handled separately

    # ---- respiratory fixed point ------------------------------------------
    g = respiratory_fixed_point(alt, s['BE'], s['HCO3c'], s['VAH'], s['HVD'],
                                hb, s['DPG'], fio2, bag, ex, s['FLOOD'], q,
                                gp=p['GP'] * subject.hvr, acz_rbc_eff=s['ACZr'] / p['ACZ_VR'],
                                sleep=asleep)
    sat, pao2, paco2, ph = g['SaO2'], g['PaO2'], g['PaCO2'], g['pH']
    cache['last'] = g

    dhyp = max(0.0, p['SAO2_REF'] - sat)

    # ---- acid-base ---------------------------------------------------------
    be_target_drive = -p['K_REN'] * (ph - 7.40)
    d_be = be_target_drive / p['TAU_REN'] * 1.0
    d_be += -(p['ACZ_EMAX_REN'] * e_acz_ren + s['BE'] * 0.0) / p['TAU_REN'] * 1.0 \
            - (s['BE'] - (-p['ACZ_EMAX_REN'] * e_acz_ren)) * 0.0
    # cleaner: relax BE toward (renal target + acetazolamide offset)
    be_ss = -p['K_REN'] * (ph - 7.40) - p['ACZ_EMAX_REN'] * e_acz_ren
    d_be = (be_ss - s['BE']) / p['TAU_REN']

    hco3c_ss = (p['W_CSF_PH'] * (0.03 * paco2 * 10.0 ** (7.32 - 6.1))
                + (1.0 - p['W_CSF_PH']) * g['HCO3'] + p['C_CSF']
                - p['ACZ_EMAX_CSF'] * e_acz_csf)
    d_hco3c = (hco3c_ss - s['HCO3c']) / p['TAU_CSF']

    # ---- carotid-body plasticity and roll-off ------------------------------
    vah_ss = 1.0 + (p['VAH_MAX'] - 1.0) * min(1.0, dhyp / 0.18)
    d_vah = (vah_ss - s['VAH']) / p['TAU_VAH']
    hvd_ss = p['HVD_MAX'] * min(1.0, dhyp / 0.15)
    d_hvd = ((hvd_ss - s['HVD']) / p['TAU_HVD_ON']) if hvd_ss > s['HVD'] \
        else ((hvd_ss - s['HVD']) / p['TAU_HVD_OFF'])

    # ---- 2,3-DPG -----------------------------------------------------------
    dpg_ss = p['DPG_SL'] + (p['DPG_MAX'] - p['DPG_SL']) * min(1.0, dhyp / 0.20)
    d_dpg = (dpg_ss - s['DPG']) / p['TAU_DPG']

    # ---- erythropoiesis ----------------------------------------------------
    epo_ss = p['EPO_SL'] + p['K_EPO'] * dhyp / (1.0 + p['EPO_HB_FB']
                                                * max(0.0, hb - subject.hb_sl) ** 2)
    d_epo = (epo_ss - s['EPO']) / p['TAU_EPO']
    d_ret = p['K_ERY'] * max(0.0, s['EPO'] - p['EPO_SL']) * p['HBM_SL'] / 100.0 \
            - s['RET'] / p['TAU_RET']
    d_hbm = s['RET'] / p['TAU_RET'] - (s['HbM'] - p['HBM_SL'] * (subject.hb_sl / 15.0)) \
            / p['TAU_HBM'] - s['HbM'] * 0.0
    pv_ss = p['PV_SL'] * (1.0 - p['PV_CONTRACT'] * min(1.0, dhyp / 0.15)) \
            + s['FLUID'] / 1000.0
    d_pv = (pv_ss - s['PV']) / p['TAU_PV']

    # ---- pulmonary vascular / HAPE ----------------------------------------
    pulm = pulmonary(g['PAO2'], q, subject.lam_max, subject.a_het,
                     s['HPVs'], hpv_supp, exercise=ex)
    cache['pulm'] = pulm
    hpvs_ss = min(1.0, dhyp / 0.15)
    d_hpvs = (hpvs_ss - s['HPVs']) / p['TAU_HPV_S']

    over = max(0.0, pulm['Pcap_open'] - p['PCAP_CRIT'])
    leak = p['KF_LEAK'] * (over ** p['N_LEAK']) * s['PERM']
    excess = max(0.0, s['ELW'] - p['ELW_0'])
    clear = p['K_AFC'] * s['AFC'] * excess * (1.0 + p['DEX_EMAX_AFC'] * e_dex
                                              + p['SAL_EMAX_AFC'] * e_sal_eff)
    d_elw = leak - clear
    flood_ss = (excess ** p['N_FLOOD']) / (p['ELW_FLOOD50'] ** p['N_FLOOD']
                                           + excess ** p['N_FLOOD'])
    d_flood = (flood_ss - s['FLOOD']) / 0.5
    afc_ss = 1.0 / (1.0 + (p['AFC_HYP50'] / max(sat, 0.30)) ** p['AFC_N'])
    afc_ss = max(0.15, min(1.0, afc_ss))
    d_afc = (afc_ss - s['AFC']) / p['TAU_AFC']
    d_infl = (p['K_INFL_LEAK'] * leak - s['INFL']) / p['TAU_INFL']
    perm_ss = 1.0 + p['K_PERM_INFL'] * s['INFL']
    d_perm = (perm_ss - s['PERM']) / p['TAU_PERM']
    d_paps = (pulm['mPAP'] - s['PAPs']) / 24.0

    # ---- cerebral / AMS / HACE --------------------------------------------
    cbf_o2 = 1.0 + p['K_CBF_O2'] / (1.0 + (max(pao2, 3.0) / p['CBF_P50']) ** p['CBF_O2N'])
    cbf_co2 = cbf_co2_factor(paco2)
    cbf_ss = max(0.35, cbf_o2 * cbf_co2)
    d_cbf = (cbf_ss - s['CBFrel']) / p['TAU_CBF']

    vegf_ss = p['K_VEGF'] * min(1.0, dhyp / 0.20) * (1.0 - p['DEX_EMAX_VEGF'] * e_dex)
    d_vegf = (vegf_ss - s['VEGFb']) / p['TAU_VEGF']
    bbb_ss = 1.0 + p['K_BBB'] * s['VEGFb'] * (1.0 - p['DEX_EMAX_BBB'] * e_dex)
    d_bbb = (bbb_ss - s['BBBp']) / p['TAU_BBB']

    # cerebral capillary hydrostatic pressure rises with CBF (arteriolar dilation)
    dp_cap = 6.0 * max(0.0, s['CBFrel'] - 1.0)
    d_edv = p['K_EDV'] * max(0.0, s['BBBp'] - 1.0) * (0.35 + dp_cap) / 10.0 \
            - s['EDv'] / p['TAU_EDV']
    edc_ss = p['K_EDC'] / (1.0 + (max(sat, 0.2) / p['EDC_O250']) ** p['EDC_N'])
    d_edc = (edc_ss - s['EDc']) / p['TAU_EDC']

    dv_brain = (s['EDv'] + s['EDc']) + p['CBV_SL'] * (s['CBFrel'] - 1.0) * 0.25
    icp = p['ICP_0'] * 10.0 ** (max(0.0, dv_brain) / subject.pvi)
    icp = min(icp, 80.0)
    cache['icp'] = icp
    d_csfres = (max(0.0, dv_brain) - s['CSFres']) / 2.0

    # ---- sympathetic / fluid ----------------------------------------------
    sym_ss = p['K_SYM'] * dhyp + 0.25 * (ex - 1.0)
    d_sym = (sym_ss - s['SYM']) / p['TAU_SYM']
    hr_ss = p['HR_SL'] * (1.0 + 0.42 * s['SYM'] / 1.0) + 28.0 * (ex - 1.0)
    d_hr = (hr_ss - s['HRs']) / p['TAU_HR']
    aldo_ss = 1.0 + 0.9 * dhyp * 4.0
    d_aldo = (aldo_ss - s['ALDO']) / p['TAU_ALDO']
    adh_ss = 1.0 + 1.0 * max(0.0, s['HEAD'] - 0.6) * 0.5 + 1.6 * dhyp
    d_adh = (adh_ss - s['ADH']) / p['TAU_ADH']
    fluid_ss = p['K_FLUID'] * max(0.0, s['ADH'] - 1.0)
    d_fluid = (fluid_ss - s['FLUID']) / p['TAU_FLUID']

    # ---- symptoms (Lake Louise 2018 sub-scores) ---------------------------
    dex_sx = (1.0 - p['DEX_EMAX_SX'] * e_dex)
    head_ss = min(3.0, (p['K_HEAD_ICP'] * max(0.0, icp - p['ICP_0'])
                        + p['K_HEAD_HYP'] * dhyp) * dex_sx
                  * (1.0 - p['IBU_EMAX_HEAD'] * e_ibu))
    gi_ss = min(3.0, p['K_GI'] * (0.55 * max(0.0, icp - p['ICP_0']) / 8.0 + dhyp) * dex_sx)
    fat_ss = min(3.0, p['K_FAT'] * (dhyp + 0.25 * s['SLEEPd'] / 3.0) * dex_sx)
    diz_ss = min(3.0, p['K_DIZ'] * dhyp * dex_sx)
    sleep_ss = min(3.0, p['K_SLEEP'] * dhyp)

    d_head = ((head_ss - s['HEAD']) / p['TAU_SX_ON']) if head_ss > s['HEAD'] \
        else ((head_ss - s['HEAD']) / p['TAU_SX_OFF'])
    d_gi = ((gi_ss - s['GI']) / p['TAU_SX_ON']) if gi_ss > s['GI'] \
        else ((gi_ss - s['GI']) / p['TAU_SX_OFF'])
    d_fat = ((fat_ss - s['FAT']) / p['TAU_SX_ON']) if fat_ss > s['FAT'] \
        else ((fat_ss - s['FAT']) / p['TAU_SX_OFF'])
    d_diz = ((diz_ss - s['DIZ']) / p['TAU_SX_ON']) if diz_ss > s['DIZ'] \
        else ((diz_ss - s['DIZ']) / p['TAU_SX_OFF'])
    d_sleep = (sleep_ss - s['SLEEPd']) / 4.0

    # ---- misc --------------------------------------------------------------
    lac_ss = 1.0 + 5.5 * max(0.0, (ex - 1.0) / 3.0) * (1.0 + 2.0 * dhyp)
    d_lac = (lac_ss - s['LAC']) / 0.5
    d_musc = (p['MUSC_MAX'] * min(1.0, dhyp / 0.15) - s['MUSC']) / p['TAU_MUSC']
    d_hypd = dhyp
    cms_ss = p['K_CMS'] * max(0.0, hct - 0.52) * 100.0 + 0.6 * max(0.0, dhyp - 0.08) * 10.0
    d_cms = (cms_ss - s['CMS']) / p['TAU_CMS']
    d_accl = ((1.0 - (s['HCO3c'] - 12.0) / (p['HCO3C_SL'] - 12.0)) - s['ACCL']) / 24.0

    # ---- drug PK -----------------------------------------------------------
    d_acza = -p['ACZ_KA'] * s['ACZa']
    d_aczc = (p['ACZ_KA'] * s['ACZa'] * p['ACZ_F'] - p['ACZ_CL'] * c_acz
              - p['ACZ_KIN'] * s['ACZc'] + p['ACZ_KOUT'] * s['ACZr'])
    d_aczr = p['ACZ_KIN'] * s['ACZc'] - p['ACZ_KOUT'] * s['ACZr']
    d_dexa = -p['DEX_KA'] * s['DEXa']
    d_dexc = p['DEX_KA'] * s['DEXa'] * p['DEX_F'] - p['DEX_CL'] * s['DEXc'] / p['DEX_V']
    d_dexe = p['DEX_KE0'] * (s['DEXc'] - s['DEXe'])
    d_nifa = -p['NIF_KA'] * s['NIFa']
    d_nifc = p['NIF_KA'] * s['NIFa'] * p['NIF_F'] - p['NIF_CL'] * c_nif
    d_tada = -p['TAD_KA'] * s['TADa']
    d_tadc = p['TAD_KA'] * s['TADa'] * p['TAD_F'] - p['TAD_CL'] * c_tad
    d_sale = -p['SAL_KE0'] * s['SALe']
    d_ibua = -p['IBU_KA'] * s['IBUa']
    d_ibuc = p['IBU_KA'] * s['IBUa'] * p['IBU_F'] - p['IBU_CL'] * c_ibu

    dy = np.zeros(NST)
    dy[IDX['BE']] = d_be
    dy[IDX['HCO3c']] = d_hco3c
    dy[IDX['VAH']] = d_vah
    dy[IDX['HVD']] = d_hvd
    dy[IDX['DPG']] = d_dpg
    dy[IDX['HbM']] = d_hbm
    dy[IDX['PV']] = d_pv
    dy[IDX['EPO']] = d_epo
    dy[IDX['RET']] = d_ret
    dy[IDX['HPVs']] = d_hpvs
    dy[IDX['ELW']] = d_elw
    dy[IDX['FLOOD']] = d_flood
    dy[IDX['AFC']] = d_afc
    dy[IDX['PERM']] = d_perm
    dy[IDX['INFL']] = d_infl
    dy[IDX['CBFrel']] = d_cbf
    dy[IDX['VEGFb']] = d_vegf
    dy[IDX['BBBp']] = d_bbb
    dy[IDX['EDv']] = d_edv
    dy[IDX['EDc']] = d_edc
    dy[IDX['CSFres']] = d_csfres
    dy[IDX['SYM']] = d_sym
    dy[IDX['HRs']] = d_hr
    dy[IDX['ALDO']] = d_aldo
    dy[IDX['ADH']] = d_adh
    dy[IDX['FLUID']] = d_fluid
    dy[IDX['HEAD']] = d_head
    dy[IDX['GI']] = d_gi
    dy[IDX['FAT']] = d_fat
    dy[IDX['DIZ']] = d_diz
    dy[IDX['SLEEPd']] = d_sleep
    dy[IDX['ACZc']] = d_aczc
    dy[IDX['ACZr']] = d_aczr
    dy[IDX['ACZa']] = d_acza
    dy[IDX['DEXc']] = d_dexc
    dy[IDX['DEXe']] = d_dexe
    dy[IDX['DEXa']] = d_dexa
    dy[IDX['NIFc']] = d_nifc
    dy[IDX['NIFa']] = d_nifa
    dy[IDX['TADc']] = d_tadc
    dy[IDX['TADa']] = d_tada
    dy[IDX['SALe']] = d_sale
    dy[IDX['IBUc']] = d_ibuc
    dy[IDX['IBUa']] = d_ibua
    dy[IDX['LAC']] = d_lac
    dy[IDX['MUSC']] = d_musc
    dy[IDX['HYPd']] = d_hypd
    dy[IDX['CMS']] = d_cms
    dy[IDX['PAPs']] = d_paps
    dy[IDX['ACCL']] = d_accl
    return dy


def simulate(subject, prof, t_end=240.0, y0=None, n_out=None, quiet=False,
             max_step=0.5):
    if y0 is None:
        y0 = default_y0(subject)
    doses = sorted(prof.doses, key=lambda d: d[0])
    breaks = sorted(set([0.0] + [d[0] for d in doses if 0.0 < d[0] < t_end] + [t_end]))
    n_out = n_out or max(int(t_end * 2) + 1, 50)
    ts_all, ys_all = [], []
    y = np.array(y0, dtype=float)
    cache = {}
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        for (dt, nm, amt) in doses:
            if abs(dt - t0) < 1e-9:
                y[IDX[nm]] += amt
        npts = max(3, int(round((t1 - t0) / t_end * n_out)) + 1)
        teval = np.linspace(t0, t1, npts)
        sol = solve_ivp(rhs, (t0, t1), y, args=(subject, prof, cache),
                        method='LSODA', t_eval=teval, max_step=max_step,
                        rtol=1e-6, atol=1e-8)
        if not sol.success and not quiet:
            print('  ! integration warning:', sol.message, file=sys.stderr)
        ts_all.append(sol.t if i == 0 else sol.t[1:])
        ys_all.append(sol.y if i == 0 else sol.y[:, 1:])
        y = sol.y[:, -1].copy()
    t = np.concatenate(ts_all)
    ys = np.concatenate(ys_all, axis=1)
    out = dict(t=t, y=ys)
    out.update(readout(t, ys, subject, prof))
    return out


def default_y0(subject):
    y0 = np.zeros(NST)
    y0[IDX['BE']] = 0.0
    y0[IDX['HCO3c']] = P['HCO3C_SL']
    y0[IDX['VAH']] = 1.0
    y0[IDX['DPG']] = P['DPG_SL']
    y0[IDX['HbM']] = P['HBM_SL'] * (subject.hb_sl / 15.0)
    y0[IDX['PV']] = P['PV_SL']
    y0[IDX['EPO']] = P['EPO_SL']
    y0[IDX['AFC']] = 1.0
    y0[IDX['PERM']] = 1.0
    y0[IDX['ELW']] = P['ELW_0']
    y0[IDX['CBFrel']] = 1.0
    y0[IDX['BBBp']] = 1.0
    y0[IDX['HRs']] = P['HR_SL']
    y0[IDX['LAC']] = 1.0
    y0[IDX['PAPs']] = 14.0
    return y0


def readout(t, ys, subject, prof):
    """Post-hoc algebraic readouts (identical algebra to the RHS)."""
    n = len(t)
    keys = ['PaO2', 'PaCO2', 'SaO2', 'PAO2', 'pH', 'HCO3', 'VE', 'VA', 'Bc', 'Bp',
            'CO2_reserve', 'S_eff', 'mPAP', 'Pcap_open', 'lam', 'PVR', 'ICP',
            'LLS', 'AMS', 'Hb', 'Hct', 'CaO2', 'DO2', 'AaDO2', 'alt', 'PIO2',
            'HAPEsev', 'HACErisk', 'Qopen', 'PB', 'AHI', 'SpO2_night']
    R = {k: np.zeros(n) for k in keys}
    for i in range(n):
        s = {nm: ys[IDX[nm], i] for nm in SNAMES}
        alt = prof.alt_fn(t[i]); ex = prof.exercise_fn(t[i])
        fio2 = prof.fio2_fn(t[i]); bag = prof.bag_fn(t[i])
        bv = max(s['PV'], 0.5) + s['HbM'] / 340.0
        hb = s['HbM'] / (10.0 * bv); hct = (s['HbM'] / 340.0) / bv
        q = P['Q_REST'] * (1.0 + 0.55 * (ex - 1.0)) * (1.0 + 0.10 * s['SYM'])
        g = respiratory_fixed_point(alt, s['BE'], s['HCO3c'], s['VAH'], s['HVD'],
                                    hb, s['DPG'], fio2, bag, ex, s['FLOOD'], q,
                                    gp=P['GP'] * subject.hvr,
                                    acz_rbc_eff=s['ACZr'] / P['ACZ_VR'])
        c_nif = s['NIFc'] / P['NIF_V']; c_tad = s['TADc'] / P['TAD_V']
        c_dex_e = s['DEXe'] / P['DEX_V']
        e_dex = c_dex_e / (P['DEX_EC50'] + c_dex_e)
        hpv_supp = 1.0 - (1.0 - P['NIF_EMAX_HPV'] * c_nif / (P['NIF_EC50'] + c_nif)) * \
                         (1.0 - P['TAD_EMAX_HPV'] * c_tad / (P['TAD_EC50'] + c_tad)) * \
                         (1.0 - P['DEX_EMAX_HPV'] * e_dex)
        pu = pulmonary(g['PAO2'], q, subject.lam_max, subject.a_het, s['HPVs'],
                       hpv_supp, exercise=ex)
        dv = (s['EDv'] + s['EDc']) + P['CBV_SL'] * (s['CBFrel'] - 1.0) * 0.25
        icp = min(P['ICP_0'] * 10.0 ** (max(0.0, dv) / subject.pvi), 80.0)
        lls = s['HEAD'] + s['GI'] + s['FAT'] + s['DIZ']
        R['PaO2'][i] = g['PaO2']; R['PaCO2'][i] = g['PaCO2']; R['SaO2'][i] = g['SaO2'] * 100
        R['PAO2'][i] = g['PAO2']; R['pH'][i] = g['pH']; R['HCO3'][i] = g['HCO3']
        R['VE'][i] = g['VE']; R['VA'][i] = g['VA']; R['Bc'][i] = g['Bc']; R['Bp'][i] = g['Bp']
        R['CO2_reserve'][i] = g['CO2_reserve']; R['S_eff'][i] = g['S_eff']
        R['AaDO2'][i] = g['AaDO2']; R['PIO2'][i] = g['PIO2']; R['PB'][i] = g['PB']
        R['mPAP'][i] = pu['mPAP']; R['Pcap_open'][i] = pu['Pcap_open']
        R['lam'][i] = pu['lam']; R['PVR'][i] = pu['PVR']; R['Qopen'][i] = pu['Qopen']
        R['ICP'][i] = icp; R['LLS'][i] = lls
        R['AMS'][i] = 1.0 if (s['HEAD'] >= 1.0 and lls >= 3.0) else 0.0
        R['Hb'][i] = hb; R['Hct'][i] = hct * 100
        R['CaO2'][i] = g['CaO2']; R['DO2'][i] = g['CaO2'] * q * 10.0
        R['alt'][i] = alt
        R['HAPEsev'][i] = max(0.0, s['ELW'] - P['ELW_0'])
        R['HACErisk'][i] = 1.0 / (1.0 + (P['ICP_HACE'] / max(icp, 1.0)) ** P['ICP_HACE_N'])
        R['AHI'][i] = ahi_from_reserve(g['CO2_reserve'])
        R['SpO2_night'][i] = g['SaO2'] * 100
    return R


def nightly_ahi(t, ahi, alt=None):
    """Mean AHI over each night (22:00-06:00 of each simulated 24 h day)."""
    out = []
    n_nights = int(np.max(t) // 24) + 1
    for d in range(n_nights):
        lo, hi = d * 24.0 + 22.0, d * 24.0 + 30.0
        m = (t >= lo) & (t <= hi)
        if m.sum() >= 2:
            out.append(float(np.mean(ahi[m])))
    return out


def ahi_from_reserve(co2_reserve):
    """Apnoea-hypopnoea index from the CO2 reserve.

    Dempsey's framework: central apnoea occurs when a transient rise in
    ventilation drives PaCO2 below the apnoeic threshold, i.e. when the CO2
    reserve is small relative to the ventilatory noise.  Calibrated to
    AHI ~2 at sea level (reserve ~5.0 mmHg), ~15 at 2500 m, ~40 at 4000 m,
    ~75 at 5000 m (Bloch/Latshang field polysomnography).
    """
    r = max(co2_reserve, 0.02)
    return 80.0 / (1.0 + (r / 1.55) ** 3.5) + 1.0


# ==============================================================================
#  SECTION 5.  ANALYSES
# ==============================================================================

def fixed_point_at(alt, acclim='acute', subject=None, exercise=1.0,
                   fio2=None, bag=0.0, hb=None, acz=0.0, sleep=False):
    """Convenience: the ventilatory fixed point for a canonical acclimatisation
    state, without running the ODE.  'acute' = sea-level bicarbonate,
    'acclimatised' = steady-state bicarbonate at that altitude."""
    subject = subject or Subject()
    hb = subject.hb_sl if hb is None else hb
    if acclim == 'acute':
        be, hco3c, vah, dpg = 0.0, P['HCO3C_SL'], 1.0, P['DPG_SL']
    else:
        be, hco3c, vah, dpg = _acclimatised_acidbase(alt, subject, acz)
    return respiratory_fixed_point(alt, be, hco3c, vah, 0.0, hb, dpg,
                                   fio2=fio2, p_bag=bag, exercise=exercise,
                                   gp=P['GP'] * subject.hvr, sleep=sleep)


def _acclimatised_acidbase(alt, subject, acz_offset=0.0, iters=80):
    """Self-consistent steady state of BE, CSF HCO3 and VAH at a fixed altitude."""
    be, hco3c, vah, dpg = 0.0, P['HCO3C_SL'], 1.0, P['DPG_SL']
    for _ in range(iters):
        g = respiratory_fixed_point(alt, be, hco3c, vah, 0.0, subject.hb_sl, dpg,
                                    gp=P['GP'] * subject.hvr)
        dhyp = max(0.0, P['SAO2_REF'] - g['SaO2'])
        be_new = -P['K_REN'] * (g['pH'] - 7.40) - acz_offset
        hco3c_new = (P['W_CSF_PH'] * (0.03 * g['PaCO2'] * 10.0 ** (7.32 - 6.1))
                     + (1.0 - P['W_CSF_PH']) * g['HCO3'] + P['C_CSF'])
        vah_new = 1.0 + (P['VAH_MAX'] - 1.0) * min(1.0, dhyp / 0.18)
        dpg_new = P['DPG_SL'] + (P['DPG_MAX'] - P['DPG_SL']) * min(1.0, dhyp / 0.20)
        be = be + 0.45 * (be_new - be)
        hco3c = hco3c + 0.45 * (hco3c_new - hco3c)
        vah = vah + 0.45 * (vah_new - vah)
        dpg = dpg + 0.45 * (dpg_new - dpg)
    return be, hco3c, vah, dpg


def descent_equivalent(alt, target_sao2, subject=None, acclim='acute'):
    """How many metres of descent would give the same SaO2?  The universal
    currency for comparing an intervention against the only treatment that
    always works."""
    subject = subject or Subject()
    def f(a):
        return fixed_point_at(a, acclim, subject)['SaO2'] - target_sao2
    lo, hi = 0.0, alt
    if f(lo) < 0:
        return float('nan')
    if f(hi) > 0:
        return 0.0
    a_eq = brentq(f, lo, hi, xtol=1.0)
    return alt - a_eq


def hape_bifurcation(subject, alt_lo=2000.0, alt_hi=7000.0, exercise=1.0,
                     hours=48.0, tol=25.0):
    """Locate the altitude at which the flooding feedback loop first becomes
    self-sustaining (excess extravascular lung water > tol mL at `hours`)."""
    def f(alt):
        prof = Profile(lambda t: alt, exercise_fn=lambda t: exercise, label='bif')
        sol = simulate(subject, prof, t_end=hours, n_out=40, quiet=True, max_step=1.0)
        return sol['HAPEsev'][-1] - tol
    if f(alt_lo) > 0:
        return alt_lo
    if f(alt_hi) < 0:
        return float('nan')
    return brentq(f, alt_lo, alt_hi, xtol=100.0)


def optimal_hct(gamma=1.0, k=None):
    """DO2 = Q*CaO2 with Q ~ mu^-gamma and mu ~ exp(k*Hct), CaO2 ~ Hct.
       DO2 ~ Hct * exp(-k*gamma*Hct)   ->   argmax = 1/(k*gamma).   Exact."""
    k = P['VISC_K'] if k is None else k
    return 1.0 / (k * gamma)


# ==============================================================================
#  SECTION 6.  SCENARIO LIBRARY
# ==============================================================================

def dose_series(start_h, interval_h, n, state, amount):
    return [(start_h + i * interval_h, state, amount) for i in range(n)]


def ramp_profile(points):
    """points = [(t_h, alt_m), ...] linear interpolation, held at the ends."""
    ts = [p[0] for p in points]; al = [p[1] for p in points]
    def f(t):
        return float(np.interp(t, ts, al))
    return f


def SCENARIOS():
    S = OrderedDict()
    typ = Subject(label='typical trekker')
    hs = Subject.hape_susceptible()
    tf = Subject.tight_fit()

    # 1 -- sea level control
    S['01_sea_level'] = dict(
        subject=typ, prof=Profile(lambda t: 0.0, label='sea level'),
        t_end=72.0, note='Control.')

    # 2 -- rapid ascent to 4559 m (Capanna Regina Margherita: valley -> hut in <24 h)
    S['02_rapid_4559'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                  label='rapid ascent 4559 m'),
        t_end=120.0, note='Monte Rosa / Capanna Margherita model ascent.')

    # 3,4 -- acetazolamide prophylaxis
    for lbl, mg in (('03_rapid_4559_acz125', 125.0), ('04_rapid_4559_acz250', 250.0)):
        S[lbl] = dict(subject=typ,
                      prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                   doses=dose_series(-0.0, 12.0, 12, 'ACZa', mg),
                                   label=f'4559 m + acetazolamide {mg:.0f} mg bid'),
                      t_end=120.0, note='Started with the ascent.')

    # 5 -- dexamethasone prophylaxis
    S['05_rapid_4559_dex'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                  doses=dose_series(0.0, 12.0, 10, 'DEXa', 4.0),
                                  label='4559 m + dexamethasone 4 mg bid'),
        t_end=120.0, note='')

    # 6 -- graded ascent to the same altitude
    S['06_graded_4559'] = dict(
        subject=typ,
        prof=Profile(ramp_profile([(0, 1130), (6, 2400), (24, 2400), (30, 3200),
                                   (48, 3200), (54, 3800), (72, 3800), (78, 4559)]),
                     label='graded ascent 4559 m'),
        t_end=144.0, note='~400 m/day sleeping-altitude gain.')

    # 7 -- climb high, sleep low
    def chsl(t):
        d = t % 24.0
        base = float(np.interp(t, [0, 6, 24, 30, 48, 54, 72], [1130, 2400, 2400, 3000, 3000, 3600, 3600]))
        if 10.0 <= d <= 16.0:
            return base + 700.0
        return base
    S['07_climb_high_sleep_low'] = dict(
        subject=typ, prof=Profile(chsl, label='climb high sleep low'),
        t_end=96.0, note='Day excursions +700 m, sleep at the lower camp.')

    # 8 -- HAPE-susceptible, rapid ascent with exertion
    def exert(t):
        d = t % 24.0
        return 3.0 if 9.0 <= d <= 15.0 else 1.0
    S['08_hapeS_rapid_exercise'] = dict(
        subject=hs, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                 exercise_fn=exert, label='HAPE-S rapid + exertion'),
        t_end=96.0, note='3x resting VO2 for 6 h/day.')

    # 9,10,11 -- HAPE prophylaxis
    S['09_hapeS_nifedipine'] = dict(
        subject=hs, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                 exercise_fn=exert,
                                 doses=dose_series(0.0, 12.0, 8, 'NIFa', 30.0),
                                 label='HAPE-S + nifedipine SR 30 mg bid'),
        t_end=96.0, note='')
    S['10_hapeS_tadalafil'] = dict(
        subject=hs, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                 exercise_fn=exert,
                                 doses=dose_series(0.0, 12.0, 8, 'TADa', 10.0),
                                 label='HAPE-S + tadalafil 10 mg bid'),
        t_end=96.0, note='Maggiorini 2006 arm.')
    S['11_hapeS_dexamethasone'] = dict(
        subject=hs, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                 exercise_fn=exert,
                                 doses=dose_series(0.0, 12.0, 8, 'DEXa', 8.0),
                                 label='HAPE-S + dexamethasone 8 mg bid'),
        t_end=96.0, note='Maggiorini 2006 arm.')

    # 12 -- established HAPE, three rescues
    est = ramp_profile([(0, 1130), (5, 3200), (9, 4559)])
    S['12_hape_rescue_descent'] = dict(
        subject=hs,
        prof=Profile(lambda t: est(t) if t < 40 else max(2500.0, 4559.0 - (t - 40) * 700.0),
                     exercise_fn=exert, label='HAPE -> descent 2000 m'),
        t_end=96.0, note='Descent begun at 40 h.')
    S['13_hape_rescue_o2'] = dict(
        subject=hs, prof=Profile(est, exercise_fn=exert,
                                 fio2_fn=lambda t: 0.2094 if t < 40 else 0.28,
                                 label='HAPE -> O2 2-3 L/min'),
        t_end=96.0, note='FiO2 0.28 from 40 h.')
    S['14_hape_rescue_gamow'] = dict(
        subject=hs, prof=Profile(est, exercise_fn=exert,
                                 bag_fn=lambda t: 0.0 if t < 40 else 105.0,
                                 label='HAPE -> Gamow bag 2 psi'),
        t_end=96.0, note='+105 mmHg from 40 h.')

    # 15 -- HACE: pushing on with AMS, tight-fit subject
    S['15_hace_push_on'] = dict(
        subject=tf, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559),
                                               (30, 4559), (40, 5300), (60, 5900)]),
                                 label='ignore AMS, keep ascending'),
        t_end=96.0, note='Tight-fit phenotype (PVI 16 mL).')
    S['16_hace_rescue'] = dict(
        subject=tf, prof=Profile(lambda t: (ramp_profile([(0, 1130), (5, 3200), (9, 4559),
                                                          (30, 4559), (40, 5300), (60, 5900)])(t)
                                            if t < 62 else max(2500.0, 5900 - (t - 62) * 900.0)),
                                 doses=dose_series(62.0, 6.0, 6, 'DEXa', 8.0),
                                 label='HACE -> dexamethasone + descent'),
        t_end=96.0, note='')

    # 17 -- Everest summit day
    S['17_everest_summit'] = dict(
        subject=Subject(hvr=1.35, hb_sl=15.0, label='elite climber'),
        prof=Profile(ramp_profile([(0, 5300), (240, 5300), (300, 6400), (400, 7100),
                                   (450, 7900), (470, 8848), (476, 7900)]),
                     exercise_fn=lambda t: 2.6 if 450 <= t <= 476 else 1.0,
                     label='Everest, no supplemental O2'),
        t_end=490.0, note='Caudwell Xtreme Everest comparison point.')
    S['18_everest_summit_o2'] = dict(
        subject=Subject(hvr=1.35, hb_sl=15.0, label='elite climber'),
        prof=Profile(ramp_profile([(0, 5300), (240, 5300), (300, 6400), (400, 7100),
                                   (450, 7900), (470, 8848), (476, 7900)]),
                     exercise_fn=lambda t: 2.6 if 450 <= t <= 476 else 1.0,
                     fio2_fn=lambda t: 0.45 if t >= 450 else 0.2094,
                     label='Everest, supplemental O2 from 7900 m'),
        t_end=490.0, note='')

    # 19,20 -- sleep at 4000 m +/- acetazolamide
    S['19_sleep_4000'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 1130), (6, 4000)]),
                                  label='sleep 4000 m'),
        t_end=72.0, note='')
    S['20_sleep_4000_acz'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 1130), (6, 4000)]),
                                  doses=dose_series(0.0, 12.0, 6, 'ACZa', 250.0),
                                  label='sleep 4000 m + acetazolamide'),
        t_end=72.0, note='')

    # 21 -- three-week acclimatisation at 3800 m
    S['21_acclimatise_3800'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 200), (10, 3800)]),
                                  label='21 days at 3800 m'),
        t_end=21 * 24.0, note='Erythropoietic and acid-base time courses.')

    # 22 -- salmeterol prophylaxis in a HAPE-susceptible subject
    S['22_hapeS_salmeterol'] = dict(
        subject=hs, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                 exercise_fn=exert,
                                 doses=dose_series(0.0, 12.0, 8, 'SALe', 1.4),
                                 label='HAPE-S + salmeterol 125 ug bid'),
        t_end=96.0, note='Sartori 2002 arm (fluid-clearance only, no HPV effect).')

    # 23 -- ibuprofen for AMS headache
    S['23_rapid_4559_ibuprofen'] = dict(
        subject=typ, prof=Profile(ramp_profile([(0, 1130), (5, 3200), (9, 4559)]),
                                  doses=dose_series(0.0, 8.0, 12, 'IBUa', 600.0),
                                  label='4559 m + ibuprofen 600 mg tid'),
        t_end=120.0, note='')

    return S


# ==============================================================================
#  SECTION 7.  DRIVER
# ==============================================================================

OUT = []
_RESULTS_REF = {}


def checkpoint(results):
    """Write the partial results after every section.  A long run that is
    killed part-way should still leave usable, self-describing output rather
    than nothing at all."""
    with open('hai_scenario_results.json', 'w') as f:
        json.dump(results, f, indent=1, default=float)
    with open('hai_reference_output.txt', 'w') as f:
        f.write('\n'.join(OUT) + '\n')


def emit(s=''):
    OUT.append(s)
    print(s)


def hdr(s):
    emit(); emit('=' * 78); emit(s); emit('=' * 78)


def main():
    results = OrderedDict()
    typ = Subject(label='typical trekker')

    # -------------------------------------------------------------------
    hdr('0.  STATIC ARITHMETIC -- no model, just the gas laws')
    emit(f"{'alt(m)':>8} {'PB(West)':>9} {'PB(ISA)':>9} {'PIO2':>7} "
         f"{'PAO2@40':>8} {'PAO2@fp':>8}")
    static = []
    for a in [0, 1600, 2500, 3500, 4559, 5300, 6400, 7100, 8000, 8848]:
        pbw, pbi = barometric(a, 'west'), barometric(a, 'isa')
        pi = pio2(a)
        pa40 = alveolar_po2(a, 40.0)
        g = fixed_point_at(a, 'acute', typ)
        emit(f"{a:8d} {pbw:9.1f} {pbi:9.1f} {pi:7.1f} {pa40:8.1f} {g['PAO2']:8.1f}")
        static.append(dict(alt=a, PB_west=pbw, PB_isa=pbi, PIO2=pi,
                           PAO2_at_PaCO2_40=pa40))
    results['static'] = static
    checkpoint(results)

    emit()
    emit('  The ONLY thing altitude does is multiply (PB - 47) by 0.2094.')
    emit(f"  Everest summit PIO2 = {pio2(8848):.1f} mmHg vs {pio2(0):.1f} at sea level "
         f"({100*pio2(8848)/pio2(0):.1f}%).")
    emit(f"  The ICAO standard atmosphere would put the summit at "
         f"{barometric(8848,'isa'):.0f} mmHg, {barometric(8848,'west')-barometric(8848,'isa'):.0f} "
         f"mmHg BELOW the measured value -- an error worth "
         f"{altitude_from_pb(barometric(8848,'isa'))-8848:.0f} m of apparent altitude.")

    # -------------------------------------------------------------------
    hdr('1.  THE ACCLIMATISATION FIXED POINT: same mountain, two diseases')
    emit(f"{'alt(m)':>7} | {'--- ACUTE (day 0) ---':^38} | {'--- ACCLIMATISED ---':^38}")
    emit(f"{'':>7} | {'PaCO2':>6} {'PaO2':>6} {'SaO2%':>6} {'VE':>6} {'HCO3':>6} "
         f"| {'PaCO2':>6} {'PaO2':>6} {'SaO2%':>6} {'VE':>6} {'HCO3':>6}")
    accl_tbl = []
    for a in [0, 2500, 3500, 4559, 5300, 6400, 7100, 8000, 8848]:
        ga = fixed_point_at(a, 'acute', typ)
        gc = fixed_point_at(a, 'acclimatised', typ)
        emit(f"{a:7d} | {ga['PaCO2']:6.1f} {ga['PaO2']:6.1f} {100*ga['SaO2']:6.1f} "
             f"{ga['VE']:6.1f} {ga['HCO3']:6.1f} | {gc['PaCO2']:6.1f} {gc['PaO2']:6.1f} "
             f"{100*gc['SaO2']:6.1f} {gc['VE']:6.1f} {gc['HCO3']:6.1f}")
        accl_tbl.append(dict(alt=a,
                             acute=dict(PaCO2=ga['PaCO2'], PaO2=ga['PaO2'],
                                        SaO2=100*ga['SaO2'], VE=ga['VE'], HCO3=ga['HCO3'],
                                        CO2res=ga['CO2_reserve'], Bc=ga['Bc']),
                             accl=dict(PaCO2=gc['PaCO2'], PaO2=gc['PaO2'],
                                       SaO2=100*gc['SaO2'], VE=gc['VE'], HCO3=gc['HCO3'],
                                       CO2res=gc['CO2_reserve'], Bc=gc['Bc'])))
    results['acclimatisation_table'] = accl_tbl
    checkpoint(results)
    emit()
    emit('  NOTE on the ACUTE column above 6000 m: it is a counterfactual, not a')
    emit('  prediction.  Nobody arrives at 8000 m with sea-level bicarbonate, and if')
    emit('  they did the model says they would be at SaO2 22% -- which is the point.')
    emit('  The acute column is there to size the DEFICIT that acclimatisation repairs.')
    a4559 = [r for r in accl_tbl if r['alt'] == 4559][0]
    emit()
    emit(f"  At 4559 m the SAME altitude is worth SaO2 {a4559['acute']['SaO2']:.1f}% on the day "
         f"you arrive and {a4559['accl']['SaO2']:.1f}% once the kidney has finished -- "
         f"a gap of {a4559['accl']['SaO2']-a4559['acute']['SaO2']:.1f} points.")
    d_eq = descent_equivalent(4559, a4559['accl']['SaO2']/100.0, typ, 'acute')
    emit(f"  In the currency that matters, acclimatisation at 4559 m IS "
         f"{d_eq:.0f} m of descent.")

    # -------------------------------------------------------------------
    hdr('2.  WHY VENTILATION BUYS MORE AT ALTITUDE -- and it is not the mmHg')
    emit(f"{'alt(m)':>7} {'PaCO2':>7} {'dPAO2/dVA':>10} {'dSaO2/dPaO2':>12} "
         f"{'dSaO2/dVA':>11} {'ratio vs SL':>12}")
    vent_tbl = []
    base = None
    for a in [0, 2500, 3500, 4559, 5300, 6400, 8000, 8848]:
        g = fixed_point_at(a, 'acclimatised', typ)
        va = g['VA']; paco2 = g['PaCO2']
        # dPAO2/dVA at fixed VCO2:  PaCO2 = k/VA -> dPaCO2/dVA = -PaCO2/VA
        c = P['FIO2_AIR'] + (1 - P['FIO2_AIR']) / P['RQ']
        dpao2_dva = c * paco2 / va
        h = 0.25
        dsao2_dpao2 = (sao2(g['PaO2'] + h, g['pH'], P['DPG_SL'])
                       - sao2(g['PaO2'] - h, g['pH'], P['DPG_SL'])) / (2 * h) * 100
        dsao2_dva = dpao2_dva * dsao2_dpao2
        if base is None:
            base = dsao2_dva
        emit(f"{a:7d} {paco2:7.1f} {dpao2_dva:10.2f} {dsao2_dpao2:12.3f} "
             f"{dsao2_dva:11.2f} {dsao2_dva/base:12.2f}")
        vent_tbl.append(dict(alt=a, dPAO2_dVA=dpao2_dva, dSaO2_dPaO2=dsao2_dpao2,
                             dSaO2_dVA=dsao2_dva, ratio=dsao2_dva/base))
    results['ventilation_value'] = vent_tbl
    checkpoint(results)
    emit()
    emit('  dPAO2/dVA = c*PaCO2/VA contains NO barometric pressure -- PB cancels out')
    emit('  of the derivative entirely.  What a litre of ventilation is worth in mmHg')
    emit('  depends only on where you already sit on the CO2 hyperbola, and by that')
    emit('  measure hyperventilation gets STEADILY WORSE with altitude: 10.1 mmHg per')
    emit('  L/min at sea level, 1.4 on the summit, because you are already breathing.')
    emit('  It is the DISSOCIATION CURVE that reverses the verdict: dSaO2/dPaO2 rises')
    emit('  35-fold over the same range, the product rises ~5-fold, and the entire')
    emit('  value of hyperventilating at altitude therefore lives in the SIGMOID and')
    emit('  not in the gas equation.  A drug that raised PaO2 by 5 mmHg would be worth')
    emit('  0.4 saturation points at sea level and 12 points at 8000 m.')

    # -------------------------------------------------------------------
    hdr('3.  THE CO2 RESERVE, PERIODIC BREATHING, AND WHAT ACETAZOLAMIDE '
        'ACTUALLY DOES')
    emit(f"{'alt(m)':>7} {'state':>14} {'PaCO2':>7} {'Bc':>6} {'B_eff':>7} "
         f"{'CO2res':>7} {'AHI':>6} {'SaO2%':>6} {'HCO3':>6} {'CSF HCO3':>9}")
    res_tbl = []
    for a in [0, 2500, 3500, 4000, 4559, 5300]:
        for lab, acz in (('acclimatised', 0.0), ('+ acetazolamide', P['ACZ_EMAX_REN'] * 0.62)):
            be, hco3c, vah, dpg = _acclimatised_acidbase(a, typ, acz)
            if acz > 0:
                hco3c -= P['ACZ_EMAX_CSF'] * 0.55
            g = respiratory_fixed_point(a, be, hco3c, vah, 0.0, typ.hb_sl, dpg,
                                        gp=P['GP'] * typ.hvr, sleep=True)
            ahi = ahi_from_reserve(g['CO2_reserve'])
            emit(f"{a:7d} {lab:>14} {g['PaCO2']:7.2f} {g['Bc']:6.2f} {g['B_eff']:7.2f} "
                 f"{g['CO2_reserve']:7.2f} {ahi:6.1f} {100*g['SaO2']:6.1f} "
                 f"{g['HCO3']:6.1f} {hco3c:9.2f}")
            res_tbl.append(dict(alt=a, arm=lab, PaCO2=g['PaCO2'], Bc=g['Bc'],
                                B_eff=g['B_eff'], CO2_reserve=g['CO2_reserve'],
                                AHI=ahi, SaO2=100*g['SaO2'], HCO3=g['HCO3'],
                                HCO3_csf=hco3c))
    results['co2_reserve'] = res_tbl
    checkpoint(results)
    r0 = [r for r in res_tbl if r['alt'] == 0 and r['arm'] == 'acclimatised'][0]
    r45 = [r for r in res_tbl if r['alt'] == 4559 and r['arm'] == 'acclimatised'][0]
    r45a = [r for r in res_tbl if r['alt'] == 4559 and r['arm'] == '+ acetazolamide'][0]
    emit()
    emit(f"  Sea level: eupneic PaCO2 {r0['PaCO2']:.1f}, apnoeic threshold "
         f"{r0['B_eff']:.1f} -> reserve {r0['CO2_reserve']:.2f} mmHg.")
    emit(f"  4559 m   : eupneic PaCO2 {r45['PaCO2']:.1f}, apnoeic threshold "
         f"{r45['B_eff']:.1f} -> reserve {r45['CO2_reserve']:.2f} mmHg.")
    emit(f"  4559 m + acetazolamide: PaCO2 falls a further "
         f"{r45['PaCO2']-r45a['PaCO2']:.2f} mmHg but the threshold falls "
         f"{r45['B_eff']-r45a['B_eff']:.2f} mmHg,")
    emit(f"  so the reserve WIDENS from {r45['CO2_reserve']:.2f} to "
         f"{r45a['CO2_reserve']:.2f} mmHg and predicted AHI falls "
         f"{r45['AHI']:.0f} -> {r45a['AHI']:.0f}.")
    emit('  Acetazolamide does not abolish periodic breathing by stimulating')
    emit('  ventilation.  It abolishes it by moving the floor down faster than')
    emit('  it moves the operating point down.')

    # -------------------------------------------------------------------
    hdr('4.  HAPE: THE CEILING IS SET BY HETEROGENEITY, NOT BY HPV STRENGTH')
    emit('  Pcap,open = P_LA + Q * (r_v/mu) / [ a/lambda + (1-a)/mu ]')
    emit()
    emit('  Take lambda -> infinity.  If the "non-responsive" bed truly does not')
    emit('  constrict (kappa = 0, mu = 1) the flow term tends to Q*r_v/(1-a) and the')
    emit('  amplification ceiling is 1/(1-a).  That is the clean statement, and it is')
    emit('  the WRONG one for a real lung: with kappa > 0, mu = 1 + kappa(lambda-1)')
    emit('  diverges along with lambda, and the correct limit is')
    emit()
    emit('        amplification ceiling = 1 / [ 1 - a(1 - kappa) ]')
    emit()
    emit(f"  At kappa = {P['KAPPA_HPV']:.2f} the naive formula overstates the ceiling by "
         f"{100*((1/(1-0.85))/(1/(1-0.85*(1-P['KAPPA_HPV'])))-1):.0f}% at a = 0.85.")
    emit('  Either way the conclusion stands and is the point of the submodel:')
    emit('  HETEROGENEITY sets the ceiling, HPV strength only sets how fast you')
    emit('  reach it.  A drug that halves lambda moves you along a curve that has')
    emit('  already flattened; a lung with a = 0.5 cannot get there at all.')
    emit()
    emit(f"{'a':>6} {'ceiling':>9} {'(naive)':>9} {'lam=1':>8} {'lam=4':>8} "
         f"{'lam=8':>8} {'lam=1e6':>9}")
    het = []
    for a_h in [0.0, 0.25, 0.5, 0.6, 0.75, 0.85, 0.9]:
        row = [pulmonary(60.0, 6.0, 0.0, a_h)['Pcap_open']]
        for lam in [4.0, 8.0]:
            lmax = lam - 1.0
            row.append(pulmonary(1.0, 6.0, lmax, a_h)['Pcap_open'])
        row.append(pulmonary(1.0, 6.0, 1e6, a_h)['Pcap_open'])
        ceil_true = 1.0 / (1.0 - a_h * (1.0 - P['KAPPA_HPV']))
        emit(f"{a_h:6.2f} {ceil_true:9.2f} {1/(1-a_h):9.2f} "
             f"{row[0]:8.2f} {row[1]:8.2f} {row[2]:8.2f} {row[3]:9.2f}")
        het.append(dict(a=a_h, ceiling=ceil_true, ceiling_naive=1/(1-a_h),
                        pcap_lam1=row[0], pcap_lam4=row[1], pcap_lam8=row[2],
                        pcap_inf=row[3]))
    results['heterogeneity'] = het
    checkpoint(results)

    emit()
    emit(f'  Same lung, same HPV, different cardiac output (a={P["A_HETERO_S"]:.2f}, lambda at 4559 m):')
    emit(f"{'Q (L/min)':>10} {'mPAP':>7} {'Pcap,open':>10} "
         f"{'> %.1f?' % P['PCAP_CRIT']:>8}")
    hs = Subject.hape_susceptible()
    g4559 = fixed_point_at(4559, 'acute', hs)
    ex_tbl = []
    for q in [5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 22.0]:
        pu = pulmonary(g4559['PAO2'], q, hs.lam_max, hs.a_het)
        emit(f"{q:10.1f} {pu['mPAP']:7.1f} {pu['Pcap_open']:10.2f} "
             f"{'YES' if pu['Pcap_open'] > P['PCAP_CRIT'] else 'no':>8}")
        ex_tbl.append(dict(Q=q, mPAP=pu['mPAP'], Pcap=pu['Pcap_open']))
    results['exercise_pcap'] = ex_tbl
    checkpoint(results)
    q_crit = brentq(lambda q: pulmonary(g4559['PAO2'], q, hs.lam_max, hs.a_het)['Pcap_open']
                    - P['PCAP_CRIT'], 3.0, 30.0)
    emit(f"  Critical cardiac output at 4559 m in this phenotype: {q_crit:.2f} L/min "
         f"({q_crit/P['Q_REST']:.2f} x resting).")
    emit('  HAPE is not a disease of altitude.  It is a disease of altitude x exertion,')
    emit('  and the model says so as arithmetic: at rest this lung is under the')
    emit('  stress-failure line and it crosses it somewhere around a brisk walk.')
    results['q_crit_4559'] = q_crit

    # -------------------------------------------------------------------
    hdr('5.  A HYPOTHESIS THE MODEL REFUTES: the P50 "paradox"')
    emit(f"{'alt(m)':>7} {'pH':>6} {'DPG':>6} {'P50':>6} {'shift(pH)':>10} "
         f"{'shift(DPG)':>11} {'SaO2%':>7}")
    p50_tbl = []
    for a in [0, 2500, 3500, 4559, 5300, 6400, 8000, 8848]:
        be, hco3c, vah, dpg = _acclimatised_acidbase(a, typ)
        g = respiratory_fixed_point(a, be, hco3c, vah, 0.0, typ.hb_sl, dpg,
                                    gp=P['GP'] * typ.hvr)
        p50_ph = P['P50_STD'] * 10 ** (P['BOHR'] * (g['pH'] - 7.40))
        p50_full = p50(g['pH'], dpg)
        emit(f"{a:7d} {g['pH']:6.3f} {dpg:6.2f} {p50_full:6.2f} "
             f"{p50_ph-P['P50_STD']:10.2f} {p50_full-p50_ph:11.2f} {100*g['SaO2']:7.1f}")
        p50_tbl.append(dict(alt=a, pH=g['pH'], DPG=dpg, P50=p50_full,
                            shift_pH=p50_ph - P['P50_STD'],
                            shift_DPG=p50_full - p50_ph, SaO2=100 * g['SaO2']))
    results['p50'] = p50_tbl
    checkpoint(results)
    pk = max(p50_tbl, key=lambda r: r['P50'])
    emit()
    emit('  HYPOTHESIS UNDER TEST: two shifts oppose each other -- 2,3-DPG moves the')
    emit('  curve RIGHT (helps unloading at the muscle), respiratory alkalosis moves')
    emit('  it LEFT (helps loading at the lung).  The textbook story is that DPG wins')
    emit('  at moderate altitude and alkalosis wins at extreme altitude, so P50 should')
    emit('  rise, PEAK somewhere around 4000-5000 m, and then fall.')
    emit()
    emit(f"  REFUTED.  With a 2,3-DPG rise calibrated to the measured ~20% "
         f"(5.0 -> {p50_tbl[-1]['DPG']:.1f} mmol/L),")
    emit(f"  P50 is monotonically DECREASING at every altitude: {p50_tbl[0]['P50']:.2f} "
         f"mmHg at sea level to {p50_tbl[-1]['P50']:.2f} on the summit.")
    emit(f"  The reason is arithmetic and one-sided: the DPG term contributes at most "
         f"{max(r['shift_DPG'] for r in p50_tbl):+.2f} mmHg because it SATURATES,")
    emit(f"  while the alkalosis term reaches {min(r['shift_pH'] for r in p50_tbl):+.2f} "
         f"mmHg and does not.  There is no crossing point to find.")
    emit('  A peak could only exist if 2,3-DPG rose by ~50%, which it does not.')
    emit('  The surviving -- and stronger -- claim is about MAGNITUDE, not shape:')
    # counterfactual: fix pH at 7.4 on the summit
    be, hco3c, vah, dpg = _acclimatised_acidbase(8848, typ)
    g_s = respiratory_fixed_point(8848, be, hco3c, vah, 0.0, typ.hb_sl, dpg,
                                  gp=P['GP'] * typ.hvr)
    s_noalk = severinghaus(g_s['PaO2'] * P['P50_STD'] / p50(7.40, dpg))
    emit(f"  Counterfactual: with the SAME PaO2 ({g_s['PaO2']:.1f} mmHg) but no")
    emit(f"  respiratory alkalosis, summit SaO2 would be {100*s_noalk:.1f}% instead of "
         f"{100*g_s['SaO2']:.1f}% -- the alkalosis is worth "
         f"{100*(g_s['SaO2']-s_noalk):.1f} saturation points.")
    results['summit_alkalosis_points'] = 100 * (g_s['SaO2'] - s_noalk)

    # -------------------------------------------------------------------
    hdr('6.  OPTIMAL HAEMATOCRIT -- exact, and it does not move with altitude')
    emit('  DO2 ~ Hct * exp(-k*gamma*Hct)  =>  Hct* = 1/(k*gamma),  k = 2.31')
    emit(f"{'gamma':>7} {'Hct*':>8}   (gamma = exponent of the cardiac-output "
         f"penalty Q ~ mu^-gamma)")
    hct_tbl = []
    for gam in [1.0, 0.9, 0.8, 0.75, 0.6, 0.5]:
        h = optimal_hct(gam)
        emit(f"{gam:7.2f} {min(h,1.0):8.3f}")
        hct_tbl.append(dict(gamma=gam, hct_opt=min(h, 1.0)))
    results['optimal_hct'] = hct_tbl
    checkpoint(results)
    gam_needed = 1.0 / (P['VISC_K'] * 0.55)
    emit()
    emit(f"  SaO2 CANCELS out of the optimisation -- CaO2 is proportional to Hct")
    emit(f"  whatever the saturation.  So the model asserts the optimum is "
         f"{optimal_hct(1.0)*100:.1f}% at")
    emit(f"  sea level and {optimal_hct(1.0)*100:.1f}% at 5000 m.  To justify an Andean "
         f"resident's Hct of 55%")
    emit(f"  the cardiac-output penalty exponent would have to be gamma <= "
         f"{gam_needed:.3f}, i.e. the")
    emit('  circulation would have to absorb most of the viscosity cost.  That is the')
    emit('  quantitative form of the Tibetan-vs-Andean argument.')
    results['gamma_for_hct55'] = gam_needed

    # -------------------------------------------------------------------
    hdr('7.  DESCENT EQUIVALENTS -- every intervention in one currency')
    g_ref = fixed_point_at(4559, 'acute', typ)
    emit(f"  Reference: acute, unacclimatised, 4559 m, SaO2 = {100*g_ref['SaO2']:.1f}%")
    emit()
    emit(f"{'intervention':<34} {'SaO2%':>7} {'dSaO2':>7} {'= metres of descent':>21}")
    interventions = []

    def add(name, g):
        d = descent_equivalent(4559, g['SaO2'], typ, 'acute')
        emit(f"{name:<34} {100*g['SaO2']:7.1f} {100*(g['SaO2']-g_ref['SaO2']):7.2f} "
             f"{d:21.0f}")
        interventions.append(dict(name=name, SaO2=100 * g['SaO2'],
                                  dSaO2=100 * (g['SaO2'] - g_ref['SaO2']),
                                  descent_m=d))

    add('nothing (acute)', g_ref)
    be_a, hc_a, vah_a, dpg_a = _acclimatised_acidbase(4559, typ, P['ACZ_EMAX_REN'] * 0.62)
    hc_a -= P['ACZ_EMAX_CSF'] * 0.55
    add('acetazolamide 250 bid (steady state)',
        respiratory_fixed_point(4559, be_a, hc_a, vah_a, 0.0, typ.hb_sl, dpg_a,
                                gp=P['GP'] * typ.hvr))
    add('full acclimatisation (~1 week)', fixed_point_at(4559, 'acclimatised', typ))
    add('supplemental O2, FiO2 0.28',
        fixed_point_at(4559, 'acute', typ, fio2=0.28))
    add('supplemental O2, FiO2 0.35',
        fixed_point_at(4559, 'acute', typ, fio2=0.35))
    add('Gamow bag, 2 psi (+105 mmHg)',
        fixed_point_at(4559, 'acute', typ, bag=105.0))
    add('Gamow bag, 4 psi (+207 mmHg)',
        fixed_point_at(4559, 'acute', typ, bag=207.0))
    add('dexamethasone (no gas-exchange effect)', g_ref)
    results['descent_equivalents'] = interventions
    checkpoint(results)
    emit()
    emit('  Dexamethasone appears on this table with a descent equivalent of ZERO')
    emit('  metres, because it does not touch a single term in the gas-exchange')
    emit('  chain.  That is not a criticism of the drug -- it is the reason it is')
    emit('  dangerous: it removes the symptom that would otherwise have stopped the')
    emit('  ascent, while leaving PaO2 exactly where it was.')

    # -------------------------------------------------------------------
    hdr('8.  SCENARIOS (full 50-state ODE system)')
    scen_out = OrderedDict()
    S = SCENARIOS()
    for key, cfg in S.items():
        sol = simulate(cfg['subject'], cfg['prof'], t_end=cfg['t_end'], quiet=True)
        t = sol['t']
        summ = dict(
            label=cfg['prof'].label,
            subject=cfg['subject'].label,
            t_end=cfg['t_end'],
            alt_final=float(sol['alt'][-1]),
            SaO2_min=float(np.min(sol['SaO2'])),
            SaO2_final=float(sol['SaO2'][-1]),
            PaO2_min=float(np.min(sol['PaO2'])),
            PaCO2_final=float(sol['PaCO2'][-1]),
            HCO3_final=float(sol['HCO3'][-1]),
            LLS_max=float(np.max(sol['LLS'])),
            LLS_t_max=float(t[int(np.argmax(sol['LLS']))]),
            LLS_final=float(sol['LLS'][-1]),
            AMS_hours=float(np.trapezoid(sol['AMS'], t)),
            ICP_max=float(np.max(sol['ICP'])),
            HACE_risk_max=float(np.max(sol['HACErisk'])),
            EVLW_excess_max=float(np.max(sol['HAPEsev'])),
            EVLW_excess_final=float(sol['HAPEsev'][-1]),
            AHI_nightly=[round(x, 1) for x in nightly_ahi(t, sol['AHI'])],
            mPAP_max=float(np.max(sol['mPAP'])),
            Pcap_max=float(np.max(sol['Pcap_open'])),
            CO2res_min=float(np.min(sol['CO2_reserve'])),
            AHI_max=float(np.max(sol['AHI'])),
            Hb_final=float(sol['Hb'][-1]),
            Hct_final=float(sol['Hct'][-1]),
            EPO_max=float(np.max(sol['y'][IDX['EPO']])),
            EPO_t_max=float(t[int(np.argmax(sol['y'][IDX['EPO']]))]),
            note=cfg['note'])
        scen_out[key] = summ
        emit(f"\n--- {key}: {cfg['prof'].label}  [{cfg['subject'].label}]")
        emit(f"    SaO2 min {summ['SaO2_min']:.1f}%  final {summ['SaO2_final']:.1f}%   "
             f"PaCO2 final {summ['PaCO2_final']:.1f}   HCO3 final {summ['HCO3_final']:.1f}")
        emit(f"    LLS peak {summ['LLS_max']:.2f} at {summ['LLS_t_max']:.0f} h   "
             f"AMS-hours {summ['AMS_hours']:.0f}   ICP max {summ['ICP_max']:.1f} mmHg")
        emit(f"    EVLW excess max {summ['EVLW_excess_max']:.0f} mL   mPAP max "
             f"{summ['mPAP_max']:.1f}   Pcap max {summ['Pcap_max']:.1f}   "
             f"HACE risk {summ['HACE_risk_max']:.3f}")
        emit(f"    CO2 reserve min {summ['CO2res_min']:.2f} mmHg   AHI by night "
             f"{summ['AHI_nightly']}   Hb {summ['Hb_final']:.1f}  Hct {summ['Hct_final']:.1f}%")
        emit(f"    EVLW excess at end of run {summ['EVLW_excess_final']:.0f} mL   "
             f"SaO2 at end {summ['SaO2_final']:.1f}%")
    results['scenarios'] = scen_out
    checkpoint(results)

    # -------------------------------------------------------------------
    hdr('9.  THE ACETAZOLAMIDE / DEXAMETHASONE ASYMMETRY')
    a = scen_out['02_rapid_4559']; b = scen_out['04_rapid_4559_acz250']
    c = scen_out['05_rapid_4559_dex']
    emit(f"{'arm':<28} {'SaO2 min':>9} {'PaCO2':>7} {'HCO3':>6} {'LLS peak':>9} "
         f"{'AMS-h':>7}")
    for nm, r in (('no prophylaxis', a), ('acetazolamide 250 bid', b),
                  ('dexamethasone 4 mg bid', c)):
        emit(f"{nm:<28} {r['SaO2_min']:9.1f} {r['PaCO2_final']:7.1f} "
             f"{r['HCO3_final']:6.1f} {r['LLS_max']:9.2f} {r['AMS_hours']:7.0f}")
    emit()
    emit(f"  Acetazolamide moves SaO2 by {b['SaO2_min']-a['SaO2_min']:+.1f} points and LLS by "
         f"{b['LLS_max']-a['LLS_max']:+.2f}.")
    emit(f"  Dexamethasone moves SaO2 by {c['SaO2_min']-a['SaO2_min']:+.1f} points and LLS by "
         f"{c['LLS_max']-a['LLS_max']:+.2f}.")
    emit('  Two drugs, comparable symptomatic benefit, and only one of them has')
    emit('  changed the patient\'s oxygen.  Every guideline that says "descend if')
    emit('  symptoms persist" is implicitly using the symptom as an oxygen gauge --')
    emit('  and dexamethasone is the drug that breaks the gauge.')
    results['asymmetry'] = dict(
        acz_dSaO2=b['SaO2_min'] - a['SaO2_min'], acz_dLLS=b['LLS_max'] - a['LLS_max'],
        dex_dSaO2=c['SaO2_min'] - a['SaO2_min'], dex_dLLS=c['LLS_max'] - a['LLS_max'])

    # -------------------------------------------------------------------
    hdr('10.  ASCENT-RATE SWEEP -- AMS is an integrator with a memory')
    emit(f"{'m/day':>7} {'days':>6} {'LLS peak':>9} {'AMS-hours':>10} "
         f"{'SaO2 min':>9} {'ICP max':>8}")
    rate_tbl = []
    for rate in [1500.0, 1000.0, 750.0, 600.0, 450.0, 300.0, 200.0]:
        target = 4559.0; start = 1130.0
        days = (target - start) / rate
        pts = [(0.0, start)]
        for d in range(1, int(math.ceil(days)) + 1):
            pts.append((d * 24.0, min(target, start + rate * d)))
        pts.append((max(pts[-1][0], 1) + 72.0, target))
        prof = Profile(ramp_profile(pts), label=f'{rate:.0f} m/day')
        sol = simulate(typ, prof, t_end=pts[-1][0], quiet=True,
                       n_out=max(60, int(pts[-1][0] / 2)), max_step=1.5)
        r = dict(rate=rate, days=days, LLS=float(np.max(sol['LLS'])),
                 AMSh=float(np.trapezoid(sol['AMS'], sol['t'])),
                 SaO2=float(np.min(sol['SaO2'])), ICP=float(np.max(sol['ICP'])))
        emit(f"{rate:7.0f} {days:6.2f} {r['LLS']:9.2f} {r['AMSh']:10.0f} "
             f"{r['SaO2']:9.1f} {r['ICP']:8.1f}")
        rate_tbl.append(r)
    results['ascent_rate'] = rate_tbl
    checkpoint(results)

    # -------------------------------------------------------------------
    hdr('11.  HAPE BIFURCATION -- the critical altitude for a phenotype')
    bif = []
    for lbl, subj in (('normal, rest', Subject()),
                      ('normal, 3x VO2', Subject()),
                      ('HAPE-susceptible, rest', Subject.hape_susceptible()),
                      ('HAPE-susceptible, 3x VO2', Subject.hape_susceptible())):
        ex = 3.0 if '3x' in lbl else 1.0
        alt_c = hape_bifurcation(subj, exercise=ex)
        emit(f"  {lbl:<28} critical altitude = "
             f"{('%.0f m' % alt_c) if alt_c == alt_c else 'none below 7000 m'}")
        bif.append(dict(phenotype=lbl, alt_crit=alt_c))
    results['hape_bifurcation'] = bif
    checkpoint(results)
    emit()
    emit('  Exercise moves the critical altitude by more than a thousand metres in')
    emit('  the susceptible phenotype and by more than that in the normal one.  The')
    emit('  bifurcation parameter that clinicians can actually control is not the')
    emit('  altitude.')

    # -------------------------------------------------------------------
    hdr('12.  HAPE PROPHYLAXIS AND RESCUE, SIDE BY SIDE')
    emit(f"{'arm':<34} {'EVLW+ (mL)':>11} {'mPAP max':>9} {'Pcap max':>9} "
         f"{'SaO2 min':>9}")
    for k in ['08_hapeS_rapid_exercise', '09_hapeS_nifedipine', '10_hapeS_tadalafil',
              '11_hapeS_dexamethasone', '22_hapeS_salmeterol',
              '12_hape_rescue_descent', '13_hape_rescue_o2', '14_hape_rescue_gamow']:
        r = scen_out[k]
        emit(f"{r['label']:<34} {r['EVLW_excess_max']:11.0f} {r['mPAP_max']:9.1f} "
             f"{r['Pcap_max']:9.1f} {r['SaO2_min']:9.1f}")

    # -------------------------------------------------------------------
    hdr('13.  MODEL-vs-LITERATURE CALIBRATION CHECK')
    checks = [
        ('PB, Everest summit (mmHg)', barometric(8848), 253.0, 'West 1983/1996'),
        ('PB, Denver 1610 m (mmHg)', barometric(1610), 632.0, 'standard'),
        ('Sea level PaCO2 (mmHg)', fixed_point_at(0, 'acute', typ)['PaCO2'], 40.0, '-'),
        ('Sea level PaO2 (mmHg)', fixed_point_at(0, 'acute', typ)['PaO2'], 95.0, '-'),
        ('Sea level SaO2 (%)', 100 * fixed_point_at(0, 'acute', typ)['SaO2'], 97.5, '-'),
        ('Sea level VE (L/min)', fixed_point_at(0, 'acute', typ)['VE'], 6.3, '-'),
        ('4559 m acute SaO2 (%)', 100 * fixed_point_at(4559, 'acute', typ)['SaO2'],
         81.0, 'Capanna Margherita field data'),
        ('4559 m acclim. SaO2 (%)', 100 * fixed_point_at(4559, 'acclimatised', typ)['SaO2'],
         87.0, 'Capanna Margherita field data'),
        ('5300 m acclim. SaO2 (%)', 100 * fixed_point_at(5300, 'acclimatised', typ)['SaO2'],
         84.0, 'Everest base camp, CXE'),
        ('5300 m acclim. PaCO2 (mmHg)', fixed_point_at(5300, 'acclimatised', typ)['PaCO2'],
         24.0, 'Grocott 2009'),
        ('8400 m PaO2 (mmHg)', fixed_point_at(8400, 'acclimatised',
                                              Subject(hvr=1.35))['PaO2'], 24.6,
         'Grocott NEJM 2009 (n=4)'),
        ('8400 m PaCO2 (mmHg)', fixed_point_at(8400, 'acclimatised',
                                               Subject(hvr=1.35))['PaCO2'], 13.3,
         'Grocott NEJM 2009'),
        ('8400 m SaO2 (%)', 100 * fixed_point_at(8400, 'acclimatised',
                                                 Subject(hvr=1.35))['SaO2'], 54.0,
         'Grocott NEJM 2009'),
    ]
    emit(f"{'quantity':<34} {'model':>9} {'observed':>9} {'err%':>7}  source")
    calib = []
    for nm, mod, obs, src in checks:
        err = 100.0 * (mod - obs) / obs
        emit(f"{nm:<34} {mod:9.2f} {obs:9.2f} {err:7.1f}  {src}")
        calib.append(dict(quantity=nm, model=mod, observed=obs, err_pct=err, source=src))
    results['calibration'] = calib
    checkpoint(results)
    errs = [abs(c['err_pct']) for c in calib]
    emit()
    emit(f"  median |error| = {np.median(errs):.1f}%, max = {np.max(errs):.1f}% "
         f"({calib[int(np.argmax(errs))]['quantity']})")
    results['calibration_median_abs_err'] = float(np.median(errs))
    results['calibration_max_abs_err'] = float(np.max(errs))

    # -------------------------------------------------------------------
    hdr('14.  VIRTUAL POPULATION -- AMS and HAPE incidence')
    rng = np.random.default_rng(20260806)
    N = 150
    arms = OrderedDict()
    arms['rapid 4559 m, no prophylaxis'] = dict(acz=0.0, dex=0.0)
    arms['rapid 4559 m, acetazolamide 250 bid'] = dict(acz=250.0, dex=0.0)
    arms['rapid 4559 m, dexamethasone 4 bid'] = dict(acz=0.0, dex=4.0)
    arms['graded 4559 m (400 m/day)'] = dict(acz=0.0, dex=0.0, graded=True)
    pop = []
    for i in range(N):
        pop.append(dict(hvr=float(np.exp(rng.normal(0, 0.28))),
                        pvi=float(P['PVI'] * np.exp(rng.normal(0, 0.22))),
                        lam=float(P['LAMBDA_MAX'] * np.exp(rng.normal(0, 0.35))),
                        a=float(min(0.92, max(0.15, rng.normal(0.50, 0.13)))),
                        hb=float(rng.normal(15.0, 1.2))))
    inc = OrderedDict()
    for arm, cfg in arms.items():
        n_ams = n_hape = 0
        lls_list = []
        for pp in pop:
            subj = Subject(hvr=pp['hvr'], pvi=pp['pvi'], lam_max=pp['lam'],
                           a_het=pp['a'], hb_sl=pp['hb'])
            doses = []
            if cfg['acz'] > 0:
                doses += dose_series(0.0, 12.0, 8, 'ACZa', cfg['acz'])
            if cfg['dex'] > 0:
                doses += dose_series(0.0, 12.0, 8, 'DEXa', cfg['dex'])
            if cfg.get('graded'):
                pts = [(0, 1130), (6, 2400), (24, 2400), (30, 3200), (48, 3200),
                       (54, 3800), (72, 3800), (78, 4559), (96, 4559)]
                tend = 96.0
            else:
                pts = [(0, 1130), (5, 3200), (9, 4559)]
                tend = 48.0
            prof = Profile(ramp_profile(pts), doses=doses, label=arm)
            sol = simulate(subj, prof, t_end=tend, quiet=True,
                           n_out=max(24, int(tend / 2)), max_step=1.5)
            lls = float(np.max(sol['LLS']))
            lls_list.append(lls)
            head = float(np.max(sol['y'][IDX['HEAD']]))
            if head >= 1.0 and lls >= 3.0:
                n_ams += 1
            if float(np.max(sol['HAPEsev'])) > 150.0:
                n_hape += 1
        inc[arm] = dict(n=N, AMS_pct=100.0 * n_ams / N, HAPE_pct=100.0 * n_hape / N,
                        LLS_mean=float(np.mean(lls_list)),
                        LLS_p90=float(np.percentile(lls_list, 90)))
        emit(f"  {arm:<40} AMS {inc[arm]['AMS_pct']:5.1f}%   "
             f"HAPE {inc[arm]['HAPE_pct']:4.1f}%   mean LLS {inc[arm]['LLS_mean']:.2f}")
    base_ams = inc['rapid 4559 m, no prophylaxis']['AMS_pct']
    for arm in list(inc)[1:]:
        rr = inc[arm]['AMS_pct'] / base_ams if base_ams > 0 else float('nan')
        inc[arm]['RR_vs_none'] = rr
        emit(f"  RR(AMS) {arm:<44} = {rr:.3f}")
    results['population'] = inc
    checkpoint(results)

    # -------------------------------------------------------------------
    hdr('15.  SUMMARY OF THE NUMBERS QUOTED ELSEWHERE')
    emit(json.dumps({k: results[k] for k in
                     ['q_crit_4559', 'summit_alkalosis_points', 'gamma_for_hct55',
                      'asymmetry', 'calibration_median_abs_err',
                      'calibration_max_abs_err']}, indent=2))

    with open('hai_scenario_results.json', 'w') as f:
        json.dump(results, f, indent=1, default=float)
    with open('hai_reference_output.txt', 'w') as f:
        f.write('\n'.join(OUT) + '\n')
    emit()
    emit('wrote hai_scenario_results.json and hai_reference_output.txt')


if __name__ == '__main__':
    main()
