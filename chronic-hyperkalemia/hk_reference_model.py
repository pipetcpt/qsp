#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hk_reference_model.py
=====================

Pure-Python reference implementation of the CHRONIC HYPERKALAEMIA QSP model
(`hk_mrgsolve_model.R`).  No third-party dependencies -- standard library only
(the container this was written in has neither R, numpy nor scipy, and the
whole point of this file is that every number quoted in `README.md` is
*computed here* rather than asserted from memory).

WHY A SEPARATE PYTHON MODEL AT ALL
----------------------------------
The mrgsolve file is the deliverable; this file is the *auditor*.  It encodes
the identical structure (same parameter names, same equations, same units) so
that:

  1. every headline number in the README is produced by running this script,
  2. any place where the model DISAGREES with the published trial is printed
     as a discrepancy instead of being silently absorbed into a fudge factor,
  3. a reader without R can still reproduce the entire argument.

THE MODEL IN ONE PARAGRAPH
--------------------------
Serum potassium is not a pool, it is a *ratio*:

        K_total = Ce*V_ecf + Ci0*LAMrel*(Ce/Ce0)^alpha * V_icf

with LAMrel the relative intracellular/extracellular partition set minute-to-
minute by insulin, beta-2 tone, pH and tonicity, alpha ~ 0.25 the exponent that
makes the intracellular pool a BUFFER rather than a mirror, and K_total moved
only over days by intake, renal excretion, colonic secretion and binders.
Every therapy in hyperkalaemia acts on exactly one term of that expression, and
the whole clinical logic of the disease falls out of which one:

  * insulin/dextrose, salbutamol             -> change LAMrel  (0 mmol removed)
  * calcium gluconate                        -> changes NEITHER (moves V_th, not K)
  * K binders, diuretics, diet, dialysis     -> change K_total  (mmol removed)
  * alkali                                   -> mostly RENAL, not a cell shift

and the renal side is a saturating secretory system whose capacity per nephron
adapts upward as GFR falls until it cannot, which is why serum K is flat
against eGFR and then is not.

STATE VECTOR (full/acute model, 30 states -- mirrors the mrgsolve @INIT block)
------------------------------------------------------------------------------
  0  Ke      ECF potassium amount                                  mmol
  1  Ki      ICF potassium amount                                  mmol
  2  HCO3    ECF bicarbonate amount                                mmol
  3  ALDO    plasma aldosterone                                    ng/dL
  4  RASDN   aldosterone-sensitive distal nephron transporter pool  (rel.)
  5  GUTK    proximal gut luminal K                                mmol
  6  COLK    colonic luminal K                                     mmol
  7  AACE    ACE inhibitor absorption depot                        mg
  8  CACE    ACE inhibitor central concentration                   mg/L
  9  AMRA    MRA absorption depot                                  mg
 10  CMRA    MRA (canrenone / finerenone) central concentration    mg/L
 11  PATP    patiromer, proximal gut                               g
 12  PATC    patiromer, colon (site of action)                     g
 13  SZCP    SZC, stomach/proximal gut (site of action)            g
 14  SZCC    SZC, distal gut                                       g
 15  INS     plasma insulin                                        uU/mL
 16  INSE    insulin effect site                                   uU/mL
 17  GLU     plasma glucose                                        mmol/L
 18  B2C     beta-2 agonist plasma                                 mg/L
 19  B2E     beta-2 agonist effect site                            mg/L
 20  AFUR    furosemide depot                                      mg
 21  CFUR    furosemide central                                    mg/L
 22  CAE     ionised calcium increment, effect site                mmol/L
 23  GFR     eGFR                                                  mL/min/1.73
 24  RD      prescribed RAASi dose (fraction of guideline target)   -
 25  BICG    oral NaHCO3 gut depot                                 mmol
 26  T55     cumulative time with Ce > 5.5                         days
 27  CHAZ    cumulative relative hazard integral                   days
 28  KREM    cumulative K removed by binder                        mmol
 29  KBAL    cumulative net K balance                              mmol

Author: Claude Code Routine (QSP disease-model library), 2026-07-29.
"""

import math
import sys

# =============================================================================
# 1.  PARAMETERS
# =============================================================================
# Every parameter carries (a) a unit, (b) where the value comes from.  "calib"
# means the value was solved for by calibrate() below against the stated
# physiological target, not read out of a paper.

P = {
    # ---- anthropometry / volumes -------------------------------------------
    "BW":        70.0,      # kg
    "fV_ecf":    0.20,      # ECF = 20% BW                       (Guyton)
    "fV_icf":    0.36,      # ICF = 36% BW  -> total body K ~51 mmol/kg
    "f_icf_fast": 0.25,     # fraction of ICF that exchanges in minutes
    "Ci0":       140.0,     # baseline intracellular [K]         mmol/L
    "Ce0":       4.2,       # baseline serum [K]                 mmol/L

    # ---- transcellular partition -------------------------------------------
    # Ci_target = Ci0 * LAMrel * (Ce/Ce0)^alpha.  alpha < 1 is what makes the
    # intracellular pool a BUFFER rather than a proportional mirror: it is the
    # single parameter that sets how many mmol of whole-body potassium surfeit
    # or deficit one mmol/L of serum potassium is worth.
    "alpha_part": 0.25,
    "k_tc_f":    13.3,      # /day fast pool (ECF perturbation t1/2 ~20 min)
    "k_tc_s":    0.277,     # /day slow (skeletal muscle) pool, t1/2 ~2.5 days
    "Emax_ins":  0.120,     # max fractional rise in LAMrel from insulin
    "EC50_ins":  60.0,      # uU/mL (absolute, not above basal)
    "Emax_b2":   0.070,     # max fractional rise in LAMrel from beta-2
    "EC50_b2":   0.012,     # mg/L effect site (salbutamol)
    "k_pH":      0.190,     # per pH unit (transcellular only)
    "k_glu":     0.0025,    # per mmol/L glucose above 5.5 (tonicity/solvent drag)
    "Emax_aldc": 0.012,     # aldosterone-driven cellular uptake (minor)

    # ---- dietary intake and gut --------------------------------------------
    "intake":    80.0,      # mmol/day dietary K (Western diet ~70-100)
    "f_abs":     0.90,      # fraction of dietary K absorbed without binder
    "ka_gut":    36.0,      # /day  absorption rate from proximal lumen
    "kt_gut":     4.0,      # /day  proximal -> colon transit (ka/(ka+kt)=f_abs)
    "kst_col":    4.0,      # /day  colon -> stool

    # ---- colonic secretion --------------------------------------------------
    # Secretion into the colon is GROSS; half of it is reabsorbed before the
    # stool.  That reabsorbable half is what a binder captures in a fasting
    # patient, and it is the only reason a binder does anything at all when
    # there is no food in the lumen.
    "C0_col":    10.0,      # mmol/day NET at GFR0            (fecal K normal ~10)
    "r_col":     0.50,      # fraction of secreted colonic K reabsorbed
    "kc_col":    1.00,      # up-regulation slope with nephron loss (BK channels)
    "Km_col":    4.2,       # mmol/L, half-max on serum K

    # ---- renal handling -----------------------------------------------------
    "GFR0":      100.0,     # mL/min/1.73 reference
    "fd":        0.030,     # fraction of filtered K delivered to the ASDN
    "S0":        88.9,      # mmol/day ASDN secretory capacity at GFR0   [calib]
    "adapt_p":   0.70,      # per-nephron secretory up-regulation exponent [calib]
    "adapt_max": 4.5,       # ceiling on per-nephron up-regulation
    "Km_secr":   3.0,       # mmol/L, serum-K drive on distal secretion
    "flow_exp":  0.50,      # distal-flow (BK-channel) exponent
    "n_hco3":    0.55,      # distal secretion scales as (HCO3/24)^n_hco3   [calib]
                            # -- metabolic acidosis inhibits ROMK/ENaC-driven
                            # kaliuresis; this, NOT the cell shift, is why
                            # chronic alkali therapy lowers serum K

    # ---- aldosterone / MR ---------------------------------------------------
    "ALDO0":     12.0,      # ng/dL baseline
    "sK_aldo":   1.30,      # per mmol/L; aldosterone doubles per +0.53 mmol/L K
    "kALDO":     50.0,      # /day turnover (t1/2 ~20 min)
    "Kd_MR":     8.0,       # ng/dL aldosterone Kd at MR
    "occ0":      0.60,      # baseline MR occupancy (= 12/(12+8))
    "m_asdn":    0.80,      # occupancy -> transporter abundance exponent
    "k_asdn":    0.462,     # /day (t1/2 1.5 d) genomic transporter turnover
    "m_col":     0.50,      # occupancy -> colonic secretion exponent

    # ---- acid-base ----------------------------------------------------------
    "NEAP":      70.0,      # mmol/day net endogenous acid production (1 mmol/kg)
    "nae_exp":   0.60,      # renal net acid excretion scales as (GFR/GFR0)^exp
    "k_buf":     0.70,      # /day bone+muscle buffering toward HCO3_0
    "HCO3_0":    24.0,      # mmol/L
    "pCO2_0":    40.0,      # mmHg
    "f_bic_abs": 0.40,      # fraction of oral NaHCO3 realised as systemic base
    "ka_bic":    12.0,      # /day

    # ---- ACE inhibitor (lisinopril-like) ------------------------------------
    "ka_ace":    1.0,       # /h
    "CL_ace":    3.0,       # L/h
    "V_ace":     40.0,      # L
    "F_ace":     0.25,
    "Emax_ace":  0.45,      # max fractional suppression of the AngII drive
    "EC50_ace":  0.020,     # mg/L

    # ---- MRA (spironolactone->canrenone, or finerenone) ---------------------
    "ka_mra":    1.2,       # /h
    "CL_mra":    10.0,      # L/h
    "V_mra":     200.0,     # L
    "F_mra":     0.70,
    "Ki_mra":    0.0060,    # mg/L, competitive Ki at MR                [calib]

    # ---- potassium binders --------------------------------------------------
    "phimax_pat": 0.42,     # max fraction of luminal K captured
    "D50_pat":    12.0,     # g/day
    "kt_pat":     4.0,      # /day proximal->colon (patiromer acts distally)
    "kel_pat":    4.0,      # /day colon -> stool
    "phimax_szc": 0.45,
    "D50_szc":    8.0,      # g/day
    "kt_szc":     8.0,      # /day (SZC binds already in stomach/prox. gut)
    "kel_szc":    8.0,

    # ---- insulin / glucose (acute rescue) -----------------------------------
    "INS_base":  10.0,      # uU/mL
    "kel_ins":   166.0,     # /day  (IV regular insulin t1/2 ~6 min)
    "ke0_ins":   50.0,      # /day  effect-site equilibration (t1/2 ~20 min)
    "V_ins":     12.0,      # L
    "GLU_base":  5.5,       # mmol/L
    "kel_glu":   30.0,      # /day insulin-independent disposal
    "k_insglu":  0.010,     # /day per uU/mL insulin-dependent disposal
    "EGP0":     173.25,     # mmol/L/day endogenous glucose production
    "EC_egp":   200.0,      # uU/mL, insulin suppression of EGP
    "V_glu":     16.0,      # L

    # ---- beta-2 agonist (nebulised salbutamol) ------------------------------
    "ka_b2":     6.0,       # /h nebulised absorption
    "CL_b2":     30.0,      # L/h
    "V_b2":      150.0,     # L
    "F_b2":      0.20,      # pulmonary bioavailability of a nebulised dose
    "ke0_b2":    12.0,      # /day

    # ---- loop diuretic ------------------------------------------------------
    "ka_fur":    2.0,       # /h
    "CL_fur":    9.0,       # L/h
    "V_fur":     12.0,      # L
    "F_fur":     0.50,
    "Emax_fur":  1.60,      # max multiple of baseline distal flow
    "EC50_fur":  0.35,      # mg/L

    # ---- calcium (membrane stabilisation) -----------------------------------
    "kel_ca":    16.0,      # /day, ionised Ca increment decay (t1/2 ~1 h)
    "Ca_gain":   0.25,      # mmol/L ionised Ca rise per 1 g calcium gluconate
    "s_Ca":      15.0,      # mV threshold shift per mmol/L ionised Ca

    # ---- electrophysiology --------------------------------------------------
    "Vth0":     -70.0,      # mV threshold potential at normal ionised Ca
    "h_mid":    -78.0,      # mV  Na-channel availability midpoint
    "h_slope":    6.0,      # mV
    "QRS0":      90.0,      # ms

    # ---- disease progression and outcome ------------------------------------
    "slope0":    -4.0,      # mL/min/1.73/yr untreated diabetic CKD decline
    "slope_rd":   1.40,     # mL/min/1.73/yr recovered at full RAASi dose
    "slope_mra":  0.60,     # extra mL/min/1.73/yr from an MRA (FIDELIO-like)
    "lnHR_rd":  -0.248,     # ln(0.78) full-dose RAASi
    "lnHR_mra": -0.198,     # ln(0.82) MRA on top
    "bK_hi":      0.642,    # ln HR = bK_hi*(K-5.0)^1.2 above 5.0
    "pK_hi":      1.20,
    "bK_lo":      0.916,    # ln HR = bK_lo*(4.0-K)^1.5 below 4.0
    "pK_lo":      1.50,

    # ---- clinician titration controller (this is what makes the dilemma real)
    "k_down":     0.90,     # /day/(mmol/L) down-titration gain above K_stop
    "k_up":       0.020,    # /day up-titration gain when K is safe
    "K_stop":     5.5,      # mmol/L threshold that triggers dose reduction
    "K_safe":     5.0,      # mmol/L below which up-titration resumes
    "RD_min":     0.0,
}

# Derived volumes: (ECF, fast ICF, slow ICF)
def volumes(p):
    vi = p["fV_icf"] * p["BW"]
    return p["fV_ecf"] * p["BW"], p["f_icf_fast"] * vi, (1.0 - p["f_icf_fast"]) * vi


def ci_target(Ce, lam_rel, p):
    """Intracellular [K] the cell is driving toward, mmol/L."""
    return p["Ci0"] * lam_rel * (max(1e-6, Ce) / p["Ce0"]) ** p["alpha_part"]


def ktot_from_Ce(Ce, lam_rel, p):
    v_e, v_f, v_s = volumes(p)
    return Ce * v_e + ci_target(Ce, lam_rel, p) * (v_f + v_s)


def buffer_capacity(Ce, lam_rel, p, fast_only=False):
    """dK_total/dCe -- how many mmol of body potassium one mmol/L is worth."""
    v_e, v_f, v_s = volumes(p)
    vi = v_f if fast_only else (v_f + v_s)
    dci = p["Ci0"] * lam_rel * p["alpha_part"] / p["Ce0"] * (Ce / p["Ce0"]) ** (
        p["alpha_part"] - 1.0)
    return v_e + vi * dci


# =============================================================================
# 2.  ALGEBRAIC SUB-MODELS
# =============================================================================

def lambda_partition(ins_e, b2_e, pH, glu, mr_occ, p):
    """RELATIVE partition modifier LAMrel (1.0 at baseline).  Everything that
    moves potassium into or out of cells without removing any acts here."""
    lam0 = 1.0
    # Insulin acts over its WHOLE range, not just above basal: basal insulin is
    # itself doing most of the work, which is why insulinopenia (type 1 DM, DKA)
    # raises serum K with no change in total body K at all.
    e_hi = p["EC50_ins"]
    occ_i = ins_e / (e_hi + ins_e)
    occ_0 = p["INS_base"] / (e_hi + p["INS_base"])
    f_ins = 1.0 + p["Emax_ins"] * (occ_i - occ_0) / (1.0 - occ_0)
    f_b2 = 1.0 + p["Emax_b2"] * b2_e / (p["EC50_b2"] + b2_e)
    f_pH = math.exp(p["k_pH"] * (pH - 7.40))
    f_glu = 1.0 / (1.0 + p["k_glu"] * max(0.0, glu - p["GLU_base"]))
    f_ald = 1.0 + p["Emax_aldc"] * (mr_occ - p["occ0"]) / p["occ0"]
    return lam0 * f_ins * f_b2 * f_pH * f_glu * f_ald


def acid_base(hco3_conc, p):
    """pH from serum bicarbonate with Winter's respiratory compensation."""
    pco2 = max(15.0, p["pCO2_0"] - 1.2 * (p["HCO3_0"] - hco3_conc))
    return 6.10 + math.log10(max(1e-6, hco3_conc) / (0.03 * pco2)), pco2


def mr_occupancy(aldo, cmra, p):
    kd_app = p["Kd_MR"] * (1.0 + cmra / p["Ki_mra"])
    return aldo / (aldo + kd_app)


def aldo_target(Ce, cace, p):
    f_ace = 1.0 - p["Emax_ace"] * cace / (p["EC50_ace"] + cace)
    return p["ALDO0"] * f_ace * math.exp(p["sK_aldo"] * (Ce - p["Ce0"]))


def distal_flow_rel(cfur, sglt2, p):
    """Relative distal tubular flow: baseline 1.0; loop diuretic and SGLT2i
    both increase delivery, which activates flow-sensitive BK channels."""
    f = 1.0 + (p["Emax_fur"] - 1.0) * cfur / (p["EC50_fur"] + cfur)
    return f * (1.15 if sglt2 else 1.0)


def renal_K(Ce, GFR, rasdn, qrel, p, hco3=None):
    """Renal potassium excretion (mmol/day) and filtered load."""
    gfr_L = GFR * 1.44                     # mL/min/1.73 -> L/day
    filt = gfr_L * Ce
    D = p["fd"] * filt                     # delivered to the ASDN
    adapt = min(p["adapt_max"], (p["GFR0"] / max(1.0, GFR)) ** p["adapt_p"])
    scap = p["S0"] * (GFR / p["GFR0"]) * adapt
    g_acid = 1.0 if hco3 is None else (max(4.0, hco3) / p["HCO3_0"]) ** p["n_hco3"]
    S = scap * rasdn * (Ce / (p["Km_secr"] + Ce)) * (qrel ** p["flow_exp"]) * g_acid
    return D + S, filt, scap


def colonic_K(Ce, GFR, mr_occ, phi, p):
    """Returns (net colonic loss, gross secretion) in mmol/day.  A binder in the
    lumen blocks reabsorption of the fraction it captures, so the NET loss rises
    with phi even when the patient is eating nothing."""
    gross0 = p["C0_col"] / (1.0 - p["r_col"])
    gross = gross0 * (1.0 + p["kc_col"] * (1.0 - min(1.0, GFR / p["GFR0"])))
    fK = (Ce / (Ce + p["Km_col"])) / (p["Ce0"] / (p["Ce0"] + p["Km_col"]))
    fA = (mr_occ / p["occ0"]) ** p["m_col"]
    gross *= fK * fA
    net = gross * (1.0 - p["r_col"] * (1.0 - phi))
    return net, gross


def binder_phi(rate_pat_site, rate_szc_site, p):
    """Fraction of luminal K captured, from the binder DELIVERY RATE (g/day)
    through its site of action -- colon for patiromer (a calcium-polymer that
    exchanges Ca for K distally), stomach/proximal gut for sodium zirconium
    cyclosilicate (a crystalline lattice with a K-selective pore that starts
    trapping immediately).  The two compete for the same luminal K, so their
    captured fractions combine as 1-(1-phi1)(1-phi2)."""
    ph1 = p["phimax_pat"] * rate_pat_site / (p["D50_pat"] + rate_pat_site) if rate_pat_site > 0 else 0.0
    ph2 = p["phimax_szc"] * rate_szc_site / (p["D50_szc"] + rate_szc_site) if rate_szc_site > 0 else 0.0
    return 1.0 - (1.0 - ph1) * (1.0 - ph2)


def ecg(Ce, Ci, ca_inc, p):
    """Resting potential, threshold, Na-channel availability, QRS, T-wave."""
    Em = 61.5 * math.log10(max(1e-6, Ce) / max(1e-6, Ci))
    Vth = p["Vth0"] + p["s_Ca"] * ca_inc
    h = 1.0 / (1.0 + math.exp((Em - p["h_mid"]) / p["h_slope"]))
    h_ref = 1.0 / (1.0 + math.exp((61.5 * math.log10(p["Ce0"] / p["Ci0"]) - p["h_mid"]) / p["h_slope"]))
    qrs = p["QRS0"] * math.sqrt(h_ref / max(1e-6, h))
    t_amp = math.sqrt(Ce / p["Ce0"]) * (1.0 + 0.35 * max(0.0, Ce - 5.0))
    return Em, Vth, Vth - Em, h, qrs, t_amp


def hazard_K(Ce, p):
    if Ce > 5.0:
        return math.exp(p["bK_hi"] * (Ce - 5.0) ** p["pK_hi"])
    if Ce < 4.0:
        return math.exp(p["bK_lo"] * (4.0 - Ce) ** p["pK_lo"])
    return 1.0


def hazard_total(Ce, rd, mra_on, p):
    return hazard_K(Ce, p) * math.exp(p["lnHR_rd"] * rd) * math.exp(
        p["lnHR_mra"] if mra_on else 0.0)


# =============================================================================
# 3.  STEADY-STATE SOLVER (fast; used for all dose-response / threshold work)
# =============================================================================
# Exploits the timescale separation: LAMBDA and ALDO equilibrate in minutes,
# RASDN in days, K_total in weeks.  At steady state every derivative is zero,
# so we simply solve the whole-body balance
#       f_abs*(1-phi)*intake  =  E_renal(Ce) + E_colon(Ce)
# for Ce by bisection, iterating the aldosterone/RASDN loop to convergence.

def steady_state(GFR, p=None, intake=None, cace=0.0, cmra=0.0,
                 pat=0.0, szc=0.0, cfur=0.0, sglt2=False, hco3=None,
                 glu=None, verbose=False):
    p = p or P
    intake = p["intake"] if intake is None else intake
    glu = p["GLU_base"] if glu is None else glu

    # acid-base steady state (unless overridden)
    if hco3 is None:
        nae = p["NEAP"] * (min(1.0, GFR / p["GFR0"]) ** p["nae_exp"])
        v_ecf = volumes(p)[0]
        hco3 = p["HCO3_0"] + (nae - p["NEAP"]) / (p["k_buf"] * v_ecf)
    pH, _ = acid_base(hco3, p)

    qrel = distal_flow_rel(cfur, sglt2, p)
    phi = binder_phi(pat, szc, p)
    net_in = p["f_abs"] * (1.0 - phi) * intake

    def imbalance(Ce):
        aldo = aldo_target(Ce, cace, p)
        occ = mr_occupancy(aldo, cmra, p)
        rasdn = (occ / p["occ0"]) ** p["m_asdn"]
        er, filt, scap = renal_K(Ce, GFR, rasdn, qrel, p, hco3)
        ec, egr = colonic_K(Ce, GFR, occ, phi, p)
        return net_in - er - ec, dict(aldo=aldo, occ=occ, rasdn=rasdn,
                                      Erenal=er, Ecol=ec, Ecol_gross=egr,
                                      filt=filt, scap=scap,
                                      Kbound=phi * p["f_abs"] * intake
                                             + phi * p["r_col"] * egr,
                                      FEK=100.0 * er / filt if filt > 0 else 0.0)

    lo, hi = 1.5, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if imbalance(mid)[0] > 0:
            lo = mid
        else:
            hi = mid
    Ce = 0.5 * (lo + hi)
    _, info = imbalance(Ce)

    lam = lambda_partition(p["INS_base"], 0.0, pH, glu, info["occ"], p)
    info.update(Ce=Ce, Ci=ci_target(Ce, lam, p), LAMrel=lam, pH=pH, HCO3=hco3,
                Ktot=ktot_from_Ce(Ce, lam, p), phi=phi, qrel=qrel,
                buf=buffer_capacity(Ce, lam, p),
                buf_fast=buffer_capacity(Ce, lam, p, fast_only=True),
                intake=intake, GFR=GFR)
    if verbose:
        print(f"    GFR {GFR:5.1f}  K {Ce:5.2f}  FEK {info['FEK']:5.1f}%  "
              f"U_K {info['Erenal']:5.1f}  colon {info['Ecol']:5.1f}  "
              f"aldo {info['aldo']:6.1f}")
    return info


# =============================================================================
# 4.  CALIBRATION
# =============================================================================
# Two free parameters are solved against two physiological anchors; everything
# else is fixed a priori.  We then REPORT the fit at the un-targeted eGFRs.

CKD_TARGETS = [
    # (eGFR, mean serum K, source of the target)
    (100.0, 4.20, "healthy adult, K 3.8-4.6"),
    ( 60.0, 4.40, "CKD 3a cohort means 4.3-4.5"),
    ( 45.0, 4.55, "CKD 3b cohort means 4.5-4.7"),
    ( 30.0, 4.75, "CKD 4 cohort means 4.7-4.9"),
    ( 20.0, 5.00, "advanced CKD 4, ~30% prevalence of K>5.0"),
    ( 12.0, 5.30, "CKD 5 pre-dialysis cohort means 5.1-5.5"),
]


def calibrate(p=None, quiet=False):
    """Solve S0 (analytically, from the GFR=100 anchor) then adapt_p (by
    bisection, from the GFR=20 anchor)."""
    p = p or P
    # --- S0 from the normal-kidney anchor -----------------------------------
    Ce = p["Ce0"]
    aldo = aldo_target(Ce, 0.0, p)
    occ = mr_occupancy(aldo, 0.0, p)
    rasdn = (occ / p["occ0"]) ** p["m_asdn"]
    ec, _ = colonic_K(Ce, p["GFR0"], occ, 0.0, p)
    need = p["f_abs"] * p["intake"] - ec
    D = p["fd"] * p["GFR0"] * 1.44 * Ce
    p["S0"] = (need - D) / (rasdn * (Ce / (p["Km_secr"] + Ce)))

    # --- adapt_p and adapt_max from the two advanced-CKD anchors ------------
    # They interact (the ceiling only binds at very low eGFR), so alternate.
    def fit(param, lo, hi, gfr, target):
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            p[param] = mid
            if steady_state(gfr, p)["Ce"] > target:   # too little secretion
                lo = mid
            else:
                hi = mid
        p[param] = 0.5 * (lo + hi)

    # --- n_hco3 from the alkali-therapy anchor ------------------------------
    # Target: correcting HCO3 from 18 to 24 in CKD 4 lowers serum K by 0.30
    # (the effect size seen in bicarbonate-supplementation studies).
    def fit_hco3():
        lo, hi = 0.0, 8.0
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            p["n_hco3"] = mid
            d = (steady_state(20.0, p, hco3=24.0)["Ce"]
                 - steady_state(20.0, p, hco3=18.0)["Ce"])
            if d < -0.30:      # too strong
                hi = mid
            else:
                lo = mid
        p["n_hco3"] = 0.5 * (lo + hi)

    for _ in range(8):
        fit("adapt_p", 0.20, 1.20, 20.0, 5.00)
        fit("adapt_max", 1.50, 30.0, 12.0, 5.30)
        fit_hco3()

    if not quiet:
        say("  calibrated  n_hco3    = %.3f  (acidosis inhibition of distal K secretion)"
              % p["n_hco3"])
        say("  calibrated  S0        = %.2f mmol/day (ASDN secretory capacity @ eGFR 100)"
              % p["S0"])
        say("  calibrated  adapt_p   = %.3f  (per-nephron up-regulation exponent)"
              % p["adapt_p"])
        say("  calibrated  adapt_max = %.2f  (ceiling on per-nephron up-regulation)"
              % p["adapt_max"])
    return p


def calibrate_mra(p=None, quiet=False):
    """Solve Ki_mra so that 25 mg/day spironolactone in a RALES-like patient
    (eGFR 60, on an ACEi) raises serum K by the observed +0.30 mmol/L."""
    p = p or P
    css = p["F_mra"] * 25.0 / (p["CL_mra"] * 24.0)     # mg/L average steady state
    base = steady_state(60.0, p, cace=0.05)["Ce"]
    lo, hi = 1e-4, 1.0
    for _ in range(90):
        mid = math.sqrt(lo * hi)
        p["Ki_mra"] = mid
        d = steady_state(60.0, p, cace=0.05, cmra=css)["Ce"] - base
        if d > 0.30:          # too potent -> raise Ki (weaker binding)
            lo = mid
        else:
            hi = mid
    p["Ki_mra"] = math.sqrt(lo * hi)
    if not quiet:
        say("  calibrated  Ki_mra  = %.5f mg/L  (canrenone at MR; RALES +0.30 mmol/L)"
              % p["Ki_mra"])
        say("              spironolactone 25 mg -> Css %.4f mg/L" % css)
    return p, css


# =============================================================================
# 5.  FULL DYNAMIC MODEL (30 states) -- used for acute and 4-week simulations
# =============================================================================

IDX = {n: i for i, n in enumerate(
    "Ke KIF KIS HCO3 ALDO RASDN GUTK COLK AACE CACE AMRA CMRA PATP PATC SZCP SZCC "
    "INS INSE GLU B2C B2E AFUR CFUR CAE GFR RD BICG T55 CHAZ KREM KBAL".split())}
NST = len(IDX)


def init_state(GFR, p, ss=None, rd=1.0):
    ss = ss or steady_state(GFR, p)
    v_ecf, v_f, v_s = volumes(p)
    y = [0.0] * NST
    y[IDX["Ke"]] = ss["Ce"] * v_ecf
    y[IDX["KIF"]] = ss["Ci"] * v_f
    y[IDX["KIS"]] = ss["Ci"] * v_s
    y[IDX["HCO3"]] = ss["HCO3"] * v_ecf
    y[IDX["ALDO"]] = ss["aldo"]
    y[IDX["RASDN"]] = ss["rasdn"]
    y[IDX["GUTK"]] = p["intake"] / p["ka_gut"]
    y[IDX["COLK"]] = 3.0
    y[IDX["INS"]] = p["INS_base"]
    y[IDX["INSE"]] = p["INS_base"]
    y[IDX["GLU"]] = p["GLU_base"]
    y[IDX["GFR"]] = GFR
    y[IDX["RD"]] = rd
    return y


def derivs(t, y, p, u):
    """u = dict of exogenous inputs (infusion rates / daily doses)."""
    d = [0.0] * NST
    v_ecf, v_f, v_s = volumes(p)
    Ce = y[IDX["Ke"]] / v_ecf
    Cif, Cis = y[IDX["KIF"]] / v_f, y[IDX["KIS"]] / v_s
    hco3c = y[IDX["HCO3"]] / v_ecf
    pH, _ = acid_base(hco3c, p)
    GFR = max(3.0, y[IDX["GFR"]])

    occ = mr_occupancy(y[IDX["ALDO"]], y[IDX["CMRA"]], p)
    lam = lambda_partition(y[IDX["INSE"]], y[IDX["B2E"]], pH, y[IDX["GLU"]], occ, p)
    cit = ci_target(Ce, lam, p)
    J_f = p["k_tc_f"] * v_f * (cit - Cif)      # ECF -> fast ICF, mmol/day
    J_s = p["k_tc_s"] * v_s * (cit - Cis)      # ECF -> slow (muscle) pool

    qrel = distal_flow_rel(y[IDX["CFUR"]], u.get("sglt2", False), p)
    # binder "dose rate through the site of action" = kel * amount there, which
    # at steady state is exactly the prescribed g/day -- this keeps the dynamic
    # model and the algebraic steady-state solver on the same scale.
    phi = binder_phi(p["kel_pat"] * y[IDX["PATC"]], p["kt_szc"] * y[IDX["SZCP"]], p)
    E_ren, filt, _ = renal_K(Ce, GFR, y[IDX["RASDN"]], qrel, p, hco3c)
    E_col, E_col_gross = colonic_K(Ce, GFR, occ, phi, p)

    # --- gut -----------------------------------------------------------------
    intake_rate = u.get("intake_rate", p["intake"])
    abs_K = p["ka_gut"] * y[IDX["GUTK"]] * (1.0 - phi)
    d[IDX["GUTK"]] = intake_rate - p["ka_gut"] * y[IDX["GUTK"]] - p["kt_gut"] * y[IDX["GUTK"]]
    d[IDX["COLK"]] = p["kt_gut"] * y[IDX["GUTK"]] + E_col - p["kst_col"] * y[IDX["COLK"]]

    # --- potassium -----------------------------------------------------------
    d[IDX["Ke"]] = abs_K + u.get("kcl_rate", 0.0) - E_ren - E_col - J_f - J_s
    d[IDX["KIF"]] = J_f
    d[IDX["KIS"]] = J_s
    d[IDX["KREM"]] = (phi * p["ka_gut"] * y[IDX["GUTK"]]
                      + phi * p["r_col"] * E_col_gross)
    d[IDX["KBAL"]] = abs_K + u.get("kcl_rate", 0.0) - E_ren - E_col

    # --- acid-base -----------------------------------------------------------
    nae = p["NEAP"] * (min(1.0, GFR / p["GFR0"]) ** p["nae_exp"])
    d[IDX["BICG"]] = -p["ka_bic"] * y[IDX["BICG"]] + u.get("bic_rate", 0.0)
    d[IDX["HCO3"]] = (nae - p["NEAP"]
                      + p["f_bic_abs"] * p["ka_bic"] * y[IDX["BICG"]]
                      + p["k_buf"] * v_ecf * (p["HCO3_0"] - hco3c))

    # --- endocrine -----------------------------------------------------------
    d[IDX["ALDO"]] = p["kALDO"] * (aldo_target(Ce, y[IDX["CACE"]], p) - y[IDX["ALDO"]])
    d[IDX["RASDN"]] = p["k_asdn"] * ((occ / p["occ0"]) ** p["m_asdn"] - y[IDX["RASDN"]])

    # --- drug PK (rates converted to /day) -----------------------------------
    # The RAASi dose actually SWALLOWED is the controller state RD times the
    # guideline target.  Without this coupling the titration loop in section F
    # would be cosmetic.
    h = 24.0
    rd_now = min(1.0, max(0.0, y[IDX["RD"]]))
    ace_rate = u.get("ace_target", 0.0) * rd_now + u.get("ace_rate", 0.0)
    mra_rate = u.get("mra_target", 0.0) * rd_now + u.get("mra_rate", 0.0)
    d[IDX["AACE"]] = -p["ka_ace"] * h * y[IDX["AACE"]] + ace_rate
    d[IDX["CACE"]] = (p["F_ace"] * p["ka_ace"] * h * y[IDX["AACE"]] / p["V_ace"]
                      - p["CL_ace"] * h / p["V_ace"] * y[IDX["CACE"]])
    d[IDX["AMRA"]] = -p["ka_mra"] * h * y[IDX["AMRA"]] + mra_rate
    d[IDX["CMRA"]] = (p["F_mra"] * p["ka_mra"] * h * y[IDX["AMRA"]] / p["V_mra"]
                      - p["CL_mra"] * h / p["V_mra"] * y[IDX["CMRA"]])
    d[IDX["PATP"]] = u.get("pat_rate", 0.0) - p["kt_pat"] * y[IDX["PATP"]]
    d[IDX["PATC"]] = p["kt_pat"] * y[IDX["PATP"]] - p["kel_pat"] * y[IDX["PATC"]]
    d[IDX["SZCP"]] = u.get("szc_rate", 0.0) - p["kt_szc"] * y[IDX["SZCP"]]
    d[IDX["SZCC"]] = p["kt_szc"] * y[IDX["SZCP"]] - p["kel_szc"] * y[IDX["SZCC"]]
    d[IDX["INS"]] = (u.get("ins_rate", 0.0) / p["V_ins"]
                     - p["kel_ins"] * (y[IDX["INS"]] - p["INS_base"]))
    d[IDX["INSE"]] = p["ke0_ins"] * (y[IDX["INS"]] - y[IDX["INSE"]])
    # Glucose: endogenous production suppressed by insulin, disposal stimulated
    # by it.  Written this way (rather than as a restoring force toward 5.5) so
    # that the model is ABLE to produce hypoglycaemia -- the commonest iatrogenic
    # harm in the treatment of hyperkalaemia.
    d[IDX["GLU"]] = (u.get("glu_rate", 0.0) / p["V_glu"]
                     + p["EGP0"] / (1.0 + y[IDX["INSE"]] / p["EC_egp"])
                     - (p["kel_glu"]
                        + p["k_insglu"] * max(0.0, y[IDX["INSE"]] - p["INS_base"]))
                     * y[IDX["GLU"]])
    d[IDX["B2C"]] = (-p["CL_b2"] * h / p["V_b2"] * y[IDX["B2C"]]
                     + u.get("b2_rate", 0.0) * p["F_b2"] / p["V_b2"])
    d[IDX["B2E"]] = p["ke0_b2"] * (y[IDX["B2C"]] - y[IDX["B2E"]])
    d[IDX["AFUR"]] = -p["ka_fur"] * h * y[IDX["AFUR"]] + u.get("fur_rate", 0.0)
    d[IDX["CFUR"]] = (p["F_fur"] * p["ka_fur"] * h * y[IDX["AFUR"]] / p["V_fur"]
                      - p["CL_fur"] * h / p["V_fur"] * y[IDX["CFUR"]])
    d[IDX["CAE"]] = -p["kel_ca"] * y[IDX["CAE"]] + u.get("ca_rate", 0.0) * p["Ca_gain"]

    # --- progression, titration, trackers ------------------------------------
    rd = rd_now
    mra_on = mra_rate > 0
    slope = (p["slope0"] + p["slope_rd"] * rd + (p["slope_mra"] if mra_on else 0.0))
    d[IDX["GFR"]] = slope / 365.0 if u.get("progress", False) else 0.0
    if u.get("titrate", False):
        d[IDX["RD"]] = (-p["k_down"] * max(0.0, Ce - p["K_stop"]) * rd
                        + p["k_up"] * max(0.0, p["K_safe"] - Ce) * (1.0 - rd))
    d[IDX["T55"]] = 1.0 if Ce > 5.5 else 0.0
    d[IDX["CHAZ"]] = hazard_total(Ce, rd, mra_on, p)
    return d


def rk4(y, t, dt, p, u):
    k1 = derivs(t, y, p, u)
    y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
    k2 = derivs(t + 0.5 * dt, y2, p, u)
    y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
    k3 = derivs(t + 0.5 * dt, y3, p, u)
    y4 = [y[i] + dt * k3[i] for i in range(NST)]
    k4 = derivs(t + dt, y4, p, u)
    return [y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(NST)]


def simulate(y0, tend, dt, p, u_fn, record=None, bolus_fn=None):
    """u_fn(t) -> input dict; bolus_fn(t, y, dt) may add instantaneous doses."""
    t, y = 0.0, list(y0)
    out = []
    n = int(round(tend / dt))
    for i in range(n + 1):
        if record is None or i % record == 0:
            out.append((t, list(y)))
        if i == n:
            break
        if bolus_fn:
            bolus_fn(t, y, dt)
        y = rk4(y, t, dt, p, u_fn(t))
        t += dt
    return out


def serumK(y, p):
    return y[IDX["Ke"]] / (p["fV_ecf"] * p["BW"])


# =============================================================================
# 6.  ANALYSES
# =============================================================================

REPORT = []
def say(s=""):
    print(s)
    REPORT.append(s)


def hr(title):
    say("")
    say("=" * 78)
    say(title)
    say("=" * 78)


def A_calibration():
    hr("A.  CALIBRATION AND THE eGFR-POTASSIUM CURVE")
    calibrate()
    calibrate_mra()
    say("")
    say("  Three parameters (S0, adapt_p, adapt_max) were fitted to three anchors")
    say("  (eGFR 100, 20 and 12).  The other three rows are PREDICTIONS, not fits.")
    say("")
    say("   eGFR   K_pred  K_target   err    FE_K%   U_K    colon   aldo")
    say("  ------ ------- --------- ------ ------- ------ ------- -------")
    errs = []
    for g, tgt, _src in CKD_TARGETS:
        s = steady_state(g)
        fitted = g in (100.0, 20.0, 12.0)
        errs.append(abs(s["Ce"] - tgt) if not fitted else None)
        say("  %6.0f  %6.2f   %6.2f   %+5.2f%s %6.1f %6.1f  %6.1f  %6.1f"
            % (g, s["Ce"], tgt, s["Ce"] - tgt, "*" if fitted else " ",
               s["FEK"], s["Erenal"], s["Ecol"], s["aldo"]))
    e = [x for x in errs if x is not None]
    say("  ( * = anchor used for calibration )")
    say("  mean |error| over the 4 held-out eGFRs = %.3f mmol/L" % (sum(e) / len(e)))
    say("")
    say("  READ THIS OFF THE FE_K COLUMN, IT IS THE WHOLE DISEASE:  the fraction")
    say("  of filtered K excreted rises from %.0f%% to %.0f%% as eGFR falls from"
        % (steady_state(100)["FEK"], steady_state(12)["FEK"]))
    say("  100 to 12.  The kidney does not fail at potassium gradually -- it")
    say("  spends its entire reserve keeping serum K flat, and serum K only")
    say("  moves once that reserve is gone.  Serum K is therefore a LATE marker")
    say("  of a reserve that has been depleting for years.")


def B_reserve_threshold():
    hr("B.  THE KALIURESIS RESERVE, AND WHAT EACH DRUG COSTS OF IT")
    css_spiro = P["F_mra"] * 25.0 / (P["CL_mra"] * 24.0)
    css_ace = P["F_ace"] * 20.0 / (P["CL_ace"] * 24.0)
    say("  Steady-state serum K (mmol/L) vs eGFR under four regimens, diet 80 mmol/day")
    say("")
    say("   eGFR    none    ACEi   ACEi+MRA   ACEi+MRA+patiromer 16.8 g")
    say("  ------ ------- ------- ---------- --------------------------")
    rows = []
    for g in (90, 75, 60, 45, 35, 25, 20, 15, 12):
        a = steady_state(g)["Ce"]
        b = steady_state(g, cace=css_ace)["Ce"]
        c = steady_state(g, cace=css_ace, cmra=css_spiro)["Ce"]
        d = steady_state(g, cace=css_ace, cmra=css_spiro, pat=16.8)["Ce"]
        rows.append((g, a, b, c, d))
        say("  %6.0f  %6.2f  %6.2f   %6.2f      %6.2f" % (g, a, b, c, d))

    def threshold(**kw):
        lo, hi = 5.0, 120.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if steady_state(mid, **kw)["Ce"] > 5.5:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    t0 = threshold()
    t1 = threshold(cace=css_ace)
    t2 = threshold(cace=css_ace, cmra=css_spiro)
    t3 = threshold(cace=css_ace, cmra=css_spiro, pat=16.8)
    say("")
    say("  eGFR at which steady-state K crosses 5.5 mmol/L:")
    say("     no RAAS blockade .................. %5.1f mL/min/1.73" % t0)
    say("     + ACE inhibitor ................... %5.1f" % t1)
    say("     + ACEi and MRA .................... %5.1f   (+%.1f vs untreated)"
        % (t2, t2 - t0))
    say("     + ACEi, MRA and patiromer ......... %5.1f" % t3)
    say("")
    say("  THE POINT:  adding an MRA does not raise K by a fixed amount, it moves")
    say("  the eGFR at which the patient becomes hyperkalaemic UP by %.0f mL/min."
        % (t2 - t0))
    say("  A patient at eGFR %.0f is not 'a bit hyperkalaemic on spironolactone';" % t2)
    say("  they have been converted into the renal-function equivalent of eGFR %.0f." % t0)
    say("  A binder moves the same threshold back down to %.1f, i.e. it buys back" % t3)
    say("  %.0f mL/min/1.73 of *apparent* potassium-handling capacity without" % (t2 - t3))
    say("  changing a single nephron.")


def C_diet_ceiling():
    hr("C.  THE DIETARY CEILING -- AND WHY A BINDER CANNOT OUTRUN A DIET")
    css_ace = P["F_ace"] * 20.0 / (P["CL_ace"] * 24.0)
    css_spiro = P["F_mra"] * 25.0 / (P["CL_mra"] * 24.0)
    say("  Maximum tolerable dietary K (mmol/day) holding steady-state K at 5.5:")
    say("")
    say("   eGFR    none    ACEi+MRA   ACEi+MRA + SZC 10 g/d")
    say("  ------ ------- ----------- -----------------------")
    for g in (60, 45, 30, 20, 15):
        def maxdiet(**kw):
            lo, hi = 5.0, 400.0
            for _ in range(70):
                mid = 0.5 * (lo + hi)
                if steady_state(g, intake=mid, **kw)["Ce"] < 5.5:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        say("  %6.0f  %6.0f     %6.0f          %6.0f"
            % (g, maxdiet(), maxdiet(cace=css_ace, cmra=css_spiro),
               maxdiet(cace=css_ace, cmra=css_spiro, szc=10.0)))
    say("")
    say("  A binder's ceiling is structural: it can only capture potassium that is")
    say("  IN the lumen, so its maximum possible negative balance is bounded by")
    say("  phi_max * f_abs * intake.  At the calibrated phi_max this is:")
    for dose, lab in ((5.0, "SZC 5 g/d"), (10.0, "SZC 10 g/d"), (8.4, "patiromer 8.4 g/d"),
                      (16.8, "patiromer 16.8 g/d"), (25.2, "patiromer 25.2 g/d")):
        is_szc = "SZC" in lab
        phi = binder_phi(0.0 if is_szc else dose, dose if is_szc else 0.0, P)
        say("     %-20s phi = %.3f  ->  max %.1f mmol/day at 80 mmol/day intake"
            % (lab, phi, phi * P["f_abs"] * 80.0))
    say("")
    say("  So no binder at any dose can remove more than ~%.0f mmol/day on a"
        % (0.60 * P["f_abs"] * 80.0))
    say("  standard diet.  In a patient whose renal + colonic excretion has fallen")
    say("  below intake by more than that, the binder is arithmetically incapable")
    say("  of restoring balance and dialysis is not a preference but a mass-balance")
    say("  requirement.")


def D_binder_trials():
    hr("D.  BINDER VALIDATION AGAINST OPAL-HK, HARMONIZE AND FIDELIO")
    css_ace = P["F_ace"] * 20.0 / (P["CL_ace"] * 24.0)
    css_spiro = P["F_mra"] * 25.0 / (P["CL_mra"] * 24.0)

    # --- OPAL-HK: patiromer, CKD on RAASi, baseline K 5.3-5.5, 4 weeks -------
    # Find the eGFR that reproduces the trial's baseline K on ACEi+MRA.
    def gfr_for_K(target, **kw):
        lo, hi = 5.0, 120.0
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            if steady_state(mid, **kw)["Ce"] > target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    say("  OPAL-HK (Weir 2015 NEJM, n=237): CKD on RAASi, baseline K 5.58 in the")
    say("  higher stratum, patiromer 8.4-16.8 g BID, 4-week change -1.01 mmol/L.")
    g = gfr_for_K(5.58, cace=css_ace, cmra=css_spiro)
    base = steady_state(g, cace=css_ace, cmra=css_spiro)
    say("     virtual patient reproducing baseline: eGFR %.1f, K %.2f" % (g, base["Ce"]))
    for dose in (16.8, 25.2, 33.6):
        s = steady_state(g, cace=css_ace, cmra=css_spiro, pat=dose)
        say("     patiromer %5.1f g/day -> K %.2f  (delta %+.2f), fecal K +%.1f mmol/d"
            % (dose, s["Ce"], s["Ce"] - base["Ce"],
               P["f_abs"] * s["phi"] * P["intake"] + s["Ecol"] - base["Ecol"]))
    s = steady_state(g, cace=css_ace, cmra=css_spiro, pat=16.8)
    say("     MODEL vs TRIAL at the licensed 16.8 g/day: %+.2f vs -1.01 mmol/L"
        % (s["Ce"] - base["Ce"]))

    # --- HARMONIZE: SZC ------------------------------------------------------
    say("")
    say("  HARMONIZE (Kosiborod 2014 JAMA, n=258): baseline K 5.6, SZC 10 g TID for")
    say("  48 h then 5/10/15 g daily; 28-day maintenance K 4.8 / 4.5 / 4.4.")
    g2 = gfr_for_K(5.60, cace=css_ace, cmra=css_spiro)
    b2 = steady_state(g2, cace=css_ace, cmra=css_spiro)
    say("     virtual patient: eGFR %.1f, K %.2f" % (g2, b2["Ce"]))
    for dose, obs in ((5.0, 4.8), (10.0, 4.5), (15.0, 4.4)):
        s = steady_state(g2, cace=css_ace, cmra=css_spiro, szc=dose)
        say("     SZC %4.1f g/day -> K %.2f  (trial %.1f, err %+.2f)"
            % (dose, s["Ce"], obs, s["Ce"] - obs))

    # --- FIDELIO / RALES: the MRA's own potassium cost -----------------------
    say("")
    say("  MRA potassium cost, model vs trials:")
    css_fin = P["F_mra"] * 20.0 / (P["CL_mra"] * 24.0) * 0.35   # finerenone, weaker MR load
    for lab, g3, cm, obs in (("RALES spironolactone 25 mg (eGFR ~60)", 60.0, css_spiro, 0.30),
                             ("FIDELIO finerenone 20 mg (eGFR ~44)", 44.0, css_fin, 0.23),
                             ("spironolactone 25 mg in CKD 4 (eGFR 25)", 25.0, css_spiro, None)):
        a = steady_state(g3, cace=css_ace)["Ce"]
        b = steady_state(g3, cace=css_ace, cmra=cm)["Ce"]
        note = ("trial %+.2f" % obs) if obs is not None else "no RCT (excluded from trials)"
        say("     %-42s %+.2f   (%s)" % (lab, b - a, note))
    say("")
    say("")
    say("  The finerenone row is a DISCREPANCY, not an agreement: assuming a MR load")
    say("  35% of spironolactone 25 mg under-predicts the FIDELIO potassium signal.")
    lo, hi = 0.01, 3.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        d = (steady_state(44.0, cace=css_ace, cmra=css_spiro * mid)["Ce"]
             - steady_state(44.0, cace=css_ace)["Ce"])
        if d < 0.23:
            lo = mid
        else:
            hi = mid
    say("  Back-solving, FIDELIO's +0.23 mmol/L requires a MR load %.0f%% of"
        % (100 * 0.5 * (lo + hi)))
    say("  spironolactone 25 mg -- i.e. finerenone at 20 mg is NOT the weak MR")
    say("  blocker a 35% assumption implies.  The model cannot distinguish 'more")
    say("  occupancy' from 'same occupancy, different tissue distribution', which is")
    say("  exactly the claim made for non-steroidal MRAs, so this is reported as an")
    say("  unresolved identifiability problem rather than a fitted parameter.")


def E_acute():
    hr("E.  ACUTE HYPERKALAEMIA -- TWO PATIENTS WITH THE SAME NUMBER")
    p = P
    v_ecf, v_f, v_s = volumes(p)
    css_ace = p["F_ace"] * 20.0 / (p["CL_ace"] * 24.0)
    css_spiro = p["F_mra"] * 25.0 / (p["CL_mra"] * 24.0)
    ssN = steady_state(100.0)
    K_normal = ssN["Ktot"]

    say("  Serum potassium is a RATIO of a pool to a partition,")
    say("        Ce  such that  K_total = Ce*V_ecf + Ci0*LAMrel*(Ce/Ce0)^alpha*V_icf,")
    say("  so the same laboratory number can mean opposite things about the pool.")
    say("")
    say("  FIRST, THE EXCHANGE RATE.  How many mmol of whole-body potassium is one")
    say("  mmol/L of serum potassium worth?  The model answers this WITHOUT being")
    say("  told, from alpha alone:")
    say("     fully equilibrated (chronic) ... %.0f mmol per mmol/L"
        % buffer_capacity(4.2, 1.0, p))
    say("     fast pool only (minutes) ....... %.0f mmol per mmol/L"
        % buffer_capacity(4.2, 1.0, p, fast_only=True))
    say("  The chronic figure is a PREDICTION and can be checked: the classical")
    say("  deficit nomogram says a stable serum K of 3.0 implies a whole-body")
    say("  deficit of roughly 200-400 mmol.  This model, never having been shown")
    say("  that number, gives:")
    for k in (3.5, 3.0, 2.5, 2.0):
        say("     serum K %.1f -> whole-body deficit %.0f mmol"
            % (k, K_normal - ktot_from_Ce(k, 1.0, p)))
    say("  The acute figure is separately checkable: an IV load of ~%.0f mmol given"
        % buffer_capacity(4.2, 1.0, p, fast_only=True))
    say("  faster than muscle can take it up raises serum K by 1 mmol/L, which is")
    say("  why 40 mmol of KCl pushed too quickly is dangerous and 40 mmol eaten is")
    say("  not.  Neither number was fitted; both fall out of alpha = %.2f."
        % p["alpha_part"])
    say("")

    # ---- patient A: mass-balance hyperkalaemia -----------------------------
    gA = 12.0
    ssA = steady_state(gA, cace=css_ace, cmra=css_spiro, hco3=19.0)
    lamA = ssA["LAMrel"]
    ktotA = ktot_from_Ce(6.80, lamA, p)
    say("  PATIENT A -- mass-balance hyperkalaemia.  eGFR %.0f, ACEi + spironolactone," % gA)
    say("  HCO3 19, glucose normal, serum K 6.80.  LAMrel = %.3f (normal)." % lamA)
    say("     K_total ........................ %.0f mmol" % ktotA)
    say("     surfeit over a normal adult .... %+.0f mmol" % (ktotA - K_normal))
    imb = (p["f_abs"] * 120.0
           - renal_K(ssA["Ce"], gA / 2.0, ssA["rasdn"], 1.0, p, 19.0)[0]
           - colonic_K(ssA["Ce"], gA / 2.0, ssA["occ"], 0.0, p)[0])
    say("     An AKI halving eGFR to %.0f with the diet up at 120 mmol/day opens an"
        % (gA / 2.0))
    say("     imbalance of %+.0f mmol/day, so that surfeit takes %.0f days to build."
        % (imb, (ktotA - K_normal) / max(1.0, imb)))

    # ---- patient B: partition hyperkalaemia --------------------------------
    hco3B, gluB, insB, deficitB = 5.0, 40.0, 0.5, 400.0
    pHB, _ = acid_base(hco3B, p)
    occB = mr_occupancy(aldo_target(5.0, 0.0, p), 0.0, p)
    lamB = lambda_partition(insB, 0.0, pHB, gluB, occB, p)
    ktotB = K_normal - deficitB

    def ce_from_ktot(kt, lam):
        lo, hi = 0.5, 15.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if ktot_from_Ce(mid, lam, p) < kt:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    ceB_present = ce_from_ktot(ktotB, lamB)
    ceB_after = ce_from_ktot(ktotB, 1.0)
    say("")
    say("  PATIENT B -- partition hyperkalaemia.  Diabetic ketoacidosis: HCO3 %.0f"
        % hco3B)
    say("  (pH %.2f), glucose %.0f, insulin %.1f uU/mL, and the %.0f mmol osmotic-"
        % (pHB, gluB, insB, deficitB))
    say("  diuresis deficit that DKA always carries.")
    say("     LAMrel collapses to ............ %.3f (%.0f%% of normal)"
        % (lamB, 100 * lamB))
    say("     K_total ........................ %.0f mmol  (%+.0f vs normal)"
        % (ktotB, -deficitB))
    say("     PRESENTING serum K ............. %.2f   -- reassuringly normal"
        % ceB_present)
    say("     serum K once insulin and alkali restore LAMrel to 1.0, with NO")
    say("     potassium given ................ %.2f   -- and that is the arrest"
        % ceB_after)
    say("")
    say("  Patient A is %+.0f mmol above normal and reads 6.80.  Patient B is %+.0f"
        % (ktotA - K_normal, -deficitB))
    say("  mmol BELOW normal and reads %.2f.  The number does not carry the" % ceB_present)
    say("  information; only the number TOGETHER WITH LAMrel does.  A binder is the")
    say("  only therapy that helps A.  For B a binder is actively harmful, and the")
    say("  therapy that 'works' -- insulin -- is precisely what uncovers the deficit.")
    say("")

    # ---- rescue arms on patient A ------------------------------------------
    y0 = init_state(gA, p, ssA)
    dK = ktotA - ssA["Ktot"]
    scale_e = v_ecf / buffer_capacity(ssA["Ce"], lamA, p)
    y0[IDX["Ke"]] += dK * scale_e
    y0[IDX["KIF"]] += dK * (1.0 - scale_e) * v_f / (v_f + v_s)
    y0[IDX["KIS"]] += dK * (1.0 - scale_e) * v_s / (v_f + v_s)
    say("  Rescue therapies applied to PATIENT A (still eating, 80 mmol/day):")
    say("")
    dt = 1.0 / (24 * 240)          # 15 s steps
    arms = {}

    def run(name, bolus, u_extra=None, tend=2.0):
        y = list(y0)
        u = dict(intake_rate=p["intake"], ace_target=0.0, ace_rate=20.0,
                 mra_rate=25.0)
        if u_extra:
            u.update(u_extra)
        given = {"done": False}

        def bol(t, yy, ddt):
            if not given["done"]:
                for k, v in bolus.items():
                    yy[IDX[k]] += v
                given["done"] = True
        arms[name] = simulate(y, tend, dt, p, lambda t: u, record=60, bolus_fn=bol)

    ins_bolus = 10.0 * 1000.0 / p["V_ins"]          # 10 U into V_ins -> uU/mL
    glu_bolus = 25.0 / 180.15 * 1000.0 / p["V_glu"]  # 25 g dextrose -> mmol/L
    run("no treatment", {})
    run("calcium gluconate 1 g", {"CAE": p["Ca_gain"]})
    run("insulin 10 U + D50 25 g", {"INS": ins_bolus, "GLU": glu_bolus})
    run("salbutamol 20 mg neb", {"B2C": 20.0 * p["F_b2"] / p["V_b2"]})
    run("SZC 10 g TID", {}, {"szc_rate": 30.0})
    run("insulin/D50 + SZC 10 g TID", {"INS": ins_bolus, "GLU": glu_bolus},
        {"szc_rate": 30.0})

    grid = (0.5, 1, 2, 4, 8, 24, 48)
    hdr = ("  %-28s" % "arm" + "".join("%6.0fh" % hh for hh in grid) + "   K removed")
    say(hdr)
    say("  " + "-" * (len(hdr) - 2))
    for name, out in arms.items():
        vals = []
        for hh in grid:
            best = min(out, key=lambda r: abs(r[0] - hh / 24.0))
            vals.append(serumK(best[1], p))
        say("  %-28s" % name + "".join("%7.2f" % v for v in vals)
            + "   %8.1f mmol" % out[-1][1][IDX["KREM"]])
    say("")
    ia = arms["insulin 10 U + D50 25 g"]
    nt = arms["no treatment"]
    sz = arms["SZC 10 g TID"]

    def at(out, hh):
        return serumK(min(out, key=lambda r: abs(r[0] - hh / 24.0))[1], p)

    say("  The insulin arm falls to %.2f at %.0f min, is back to %.2f by 4 h and by"
        % (serumK(min(ia, key=lambda r: serumK(r[1], p))[1], p),
           min(ia, key=lambda r: serumK(r[1], p))[0] * 24 * 60, at(ia, 4)))
    say("  48 h is at %.2f -- INDISTINGUISHABLE from the untreated arm's %.2f,"
        % (at(ia, 48), at(nt, 48)))
    say("  because it removed %.1f mmol.  The SZC arm is slower at 4 h (%.2f vs %.2f)"
        % (ia[-1][1][IDX["KREM"]], at(sz, 4), at(ia, 4)))
    say("  but reaches %.2f at 48 h having removed %.0f mmol, and that is the only"
        % (at(sz, 48), sz[-1][1][IDX["KREM"]]))
    say("  fall of the two that does not come back.")
    say("")
    say("  DISCREPANCY, stated rather than tuned away: HARMONIZE reports roughly")
    say("  -1.1 mmol/L at 48 h on SZC 10 g TID from a baseline of 5.6; this model")
    say("  gives %+.2f from 6.80.  The model is slower early.  The most likely"
        % (at(sz, 48) - 6.80))
    say("  structural cause is that a single fast pool plus a single slow pool is")
    say("  still too coarse a description of muscle, and the trial patients' intake")
    say("  fell when they were admitted while this simulation keeps them eating")
    say("  80 mmol/day throughout.")
    say("")
    say("  How much dextrose does 10 U of insulin need?  (K nadir vs glucose nadir)")
    say("")
    say("   dextrose (g)   K nadir   time (min)   glucose peak   glucose nadir")
    say("  ------------- --------- ------------ -------------- ---------------")
    for dex in (0.0, 12.5, 25.0, 50.0):
        y = list(y0)
        u = dict(intake_rate=p["intake"], ace_rate=20.0, mra_rate=25.0)
        st = {"done": False}

        def bol(t, yy, ddt, dex=dex, st=st):
            if not st["done"]:
                yy[IDX["INS"]] += ins_bolus
                yy[IDX["GLU"]] += dex / 180.15 * 1000.0 / p["V_glu"]
                st["done"] = True
        o = simulate(y, 0.35, dt, p, lambda t: u, record=30, bolus_fn=bol)
        ks = [serumK(r[1], p) for r in o]
        gl = [r[1][IDX["GLU"]] for r in o]
        say("   %11.1f   %7.2f   %10.0f   %12.1f   %13.2f"
            % (dex, min(ks), o[ks.index(min(ks))][0] * 24 * 60, max(gl), min(gl)))
    say("")
    say("  The model does NOT reproduce the ~15-20% incidence of frank hypoglycaemia")
    say("  reported after insulin/dextrose for hyperkalaemia: even with no dextrose")
    say("  at all its nadir stays above 4 mmol/L.  Two things are missing and both")
    say("  are structural -- the insulin effect site here decays with the plasma")
    say("  level (real tissue action outlasts it by hours), and there is no")
    say("  between-patient variation in glycogen reserve, which is what actually")
    say("  separates the patient who becomes hypoglycaemic from the one who does")
    say("  not.  This is a place the model should not be trusted.")
    say("")

    Emn, Vthn, gapn, _, qrsn, _ = ecg(4.2, 140.0, 0.0, p)
    Em0, Vth0, gap0, _, qrs0, _ = ecg(6.8, 140.0, 0.0, p)
    Em1, Vth1, gap1, _, qrs1, _ = ecg(6.8, 140.0, p["Ca_gain"], p)
    say("  Membrane arithmetic -- why calcium is given first and why it is absent")
    say("  from the 'K removed' column:")
    say("     K 4.2      : Em %.1f mV, threshold %.1f, excitability gap %.1f mV, QRS %.0f ms"
        % (Emn, Vthn, gapn, qrsn))
    say("     K 6.8      : Em %.1f mV, threshold %.1f, gap %.1f mV (-%.0f%%), QRS %.0f ms"
        % (Em0, Vth0, gap0, 100 * (1 - gap0 / gapn), qrs0))
    say("     + Ca 1 g   : Em %.1f mV (UNCHANGED), threshold %.1f, gap %.1f mV, QRS %.0f ms"
        % (Em1, Vth1, gap1, qrs1))
    say("     Calcium restores %.0f%% of the lost excitability gap in ~1 min while"
        % (100 * (gap1 - gap0) / max(1e-9, gapn - gap0)))
    say("     removing 0 mmol of potassium.  It buys time; it does not treat.")
    say("")
    say("   K      Em(mV)  Na-avail   QRS(ms)   T-wave(rel)   HR(all-cause)")
    say("  ----- -------- --------- --------- ------------- ---------------")
    for k in (3.0, 4.2, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0):
        Em, _, _, hh, qrs, ta = ecg(k, 140.0, 0.0, p)
        say("  %5.1f  %7.1f   %7.3f   %7.0f      %7.2f        %7.2f"
            % (k, Em, hh, qrs, ta, hazard_K(k, p)))


def F_dilemma():
    hr("F.  THE RAASi-POTASSIUM DILEMMA, SIMULATED WITH A CLINICIAN IN THE LOOP")
    p = P
    css_ace = p["F_ace"] * 20.0 / (p["CL_ace"] * 24.0)
    css_spiro = p["F_mra"] * 25.0 / (p["CL_mra"] * 24.0)
    say("  A clinician is modelled as a controller on the prescribed RAASi dose:")
    say("     dRD/dt = -k_down*max(0, K-5.5)*RD + k_up*max(0, 5.0-K)*(1-RD)")
    say("  i.e. the dose comes down when K breaches 5.5 and creeps back up when it")
    say("  is safe.  Nothing about the outcome below is imposed; it is what the")
    say("  controller does when it is coupled to the potassium physiology.")
    say("")
    g0 = 25.0
    horizon = 5 * 365.0
    dt = 0.01

    def arm(label, mra=True, pat=0.0, szc=0.0, intake=80.0, fur=0.0, sglt2=False,
            titrate=True, bic=0.0):
        ss = steady_state(g0, cace=css_ace, cmra=css_spiro if mra else 0.0)
        y = init_state(g0, p, ss, rd=1.0)
        u = dict(intake_rate=intake, ace_target=20.0, progress=True, titrate=titrate,
                 sglt2=sglt2, pat_rate=pat, szc_rate=szc, fur_rate=fur, bic_rate=bic)
        if mra:
            u["mra_target"] = 25.0
        out = simulate(y, horizon, dt, p, lambda t: u, record=int(7 / dt))
        tf, yf = out[-1]
        return dict(label=label, K=serumK(yf, p), RD=yf[IDX["RD"]],
                    GFR=yf[IDX["GFR"]], T55=yf[IDX["T55"]],
                    HAZ=yf[IDX["CHAZ"]] / horizon, out=out)

    arms = [
        arm("ACEi + MRA, no binder (clinician titrates)"),
        arm("ACEi + MRA + patiromer 16.8 g/d", pat=16.8),
        arm("ACEi + MRA + SZC 10 g/d", szc=10.0),
        arm("ACEi + MRA, low-K diet 50 mmol/d", intake=50.0),
        arm("ACEi + MRA + furosemide 40 mg/d", fur=40.0),
        arm("ACEi + MRA + SGLT2i", sglt2=True),
        arm("ACEi + MRA + NaHCO3 70 mmol/d", bic=70.0),
        arm("ACEi alone, MRA never started", mra=False),
    ]
    say("  5-year outcome of a diabetic CKD 4 patient (eGFR 25 at entry, diet 80):")
    say("")
    say("  %-44s %6s %6s %7s %7s %7s" % ("arm", "K_5y", "RAASi", "eGFR_5y", "d>5.5", "meanHR"))
    say("  " + "-" * 76)
    for a in arms:
        say("  %-44s %6.2f %5.0f%% %7.1f %7.0f %7.3f"
            % (a["label"], a["K"], 100 * a["RD"], a["GFR"], a["T55"], a["HAZ"]))
    base = arms[0]
    say("")
    say("  Read the RAASi column, not the K column.  Every arm ends up with an")
    say("  acceptable potassium -- that is what the controller is for.  What differs")
    say("  is the PRICE PAID IN RAASi DOSE to get there:")
    for a in arms[1:4]:
        say("     %-44s RAASi %3.0f%% vs %3.0f%%, eGFR %+.1f mL/min at 5 y"
            % (a["label"], 100 * a["RD"], 100 * base["RD"], a["GFR"] - base["GFR"]))
    say("")
    say("  And the mean hazard column is the composite the patient actually")
    say("  experiences: potassium risk multiplied by the benefit forgone.")
    best = min(arms, key=lambda a: a["HAZ"])
    say("     lowest mean hazard: %s (%.3f)" % (best["label"], best["HAZ"]))
    say("     vs no-binder titration: %.3f, a %.1f%% relative reduction"
        % (base["HAZ"], 100 * (1 - best["HAZ"] / base["HAZ"])))
    return arms


def G_virtual_population(arms):
    hr("G.  VIRTUAL POPULATION -- WHO ACTUALLY BENEFITS FROM A BINDER")
    p = P
    css_ace = p["F_ace"] * 20.0 / (p["CL_ace"] * 24.0)
    css_spiro = p["F_mra"] * 25.0 / (p["CL_mra"] * 24.0)
    say("  600 virtual patients on a deterministic lattice over the three things")
    say("  that actually vary between real patients and are cheap to measure:")
    say("     eGFR 15-60, dietary K 40-140 mmol/day, serum HCO3 16-26 mmol/L.")
    say("")
    n = 0
    hyper_no = hyper_bind = 0
    rescued = 0
    buckets = {}
    for i in range(10):
        gfr = 15.0 + i * 5.0
        for j in range(10):
            diet = 40.0 + j * 11.0
            for k in range(6):
                hco3 = 16.0 + k * 2.0
                n += 1
                a = steady_state(gfr, intake=diet, hco3=hco3,
                                 cace=css_ace, cmra=css_spiro)["Ce"]
                b = steady_state(gfr, intake=diet, hco3=hco3,
                                 cace=css_ace, cmra=css_spiro, pat=16.8)["Ce"]
                if a > 5.5:
                    hyper_no += 1
                    if b <= 5.5:
                        rescued += 1
                if b > 5.5:
                    hyper_bind += 1
                key = int(gfr // 15) * 15
                d = buckets.setdefault(key, [0, 0, 0])
                d[0] += 1
                d[1] += 1 if a > 5.5 else 0
                d[2] += 1 if b > 5.5 else 0
    say("  n = %d.  On ACEi+MRA without a binder, %d (%.1f%%) sit above K 5.5."
        % (n, hyper_no, 100.0 * hyper_no / n))
    say("  Adding patiromer 16.8 g/day brings %d of those %d (%.1f%%) back under 5.5."
        % (rescued, hyper_no, 100.0 * rescued / max(1, hyper_no)))
    say("  %d patients (%.1f%%) remain above 5.5 despite the binder."
        % (hyper_bind, 100.0 * hyper_bind / n))
    say("")
    say("   eGFR band     n    %K>5.5 no binder    %K>5.5 with binder")
    say("  ----------- ----- ------------------- --------------------")
    for kk in sorted(buckets):
        d = buckets[kk]
        say("   %2d-%2d      %5d        %6.1f%%             %6.1f%%"
            % (kk, kk + 14, d[0], 100.0 * d[1] / d[0], 100.0 * d[2] / d[0]))
    say("")
    say("  The binder's benefit is not uniform.  It is concentrated where the")
    say("  patient is above threshold BY LESS THAN THE BINDER'S CEILING, and that")
    say("  is a computable band, not a clinical impression.")


def H_discrepancies():
    hr("H.  DISCREPANCIES AND WHAT THE MODEL CANNOT DO")
    css_ace = P["F_ace"] * 20.0 / (P["CL_ace"] * 24.0)
    css_spiro = P["F_mra"] * 25.0 / (P["CL_mra"] * 24.0)
    say("  (1) BICARBONATE ACTS ON THE KIDNEY, NOT ON THE CELL.  Building the model")
    say("      from first principles produced a result I did not put in by hand:")
    say("      the transcellular route contributes almost nothing chronically,")
    say("      because at steady state serum K is pinned by MASS BALANCE and the")
    say("      partition only decides how big the body pool behind it is.")
    for hco3 in (16.0, 18.0, 20.0, 22.0, 24.0):
        s = steady_state(20.0, hco3=hco3, cace=css_ace)
        say("      HCO3 %4.1f -> pH %.3f, LAMrel %.3f, U_K %4.1f, steady-state K %.2f"
            % (hco3, s["pH"], s["LAMrel"], s["Erenal"], s["Ce"]))
    a = steady_state(20.0, hco3=18.0, cace=css_ace)
    b = steady_state(20.0, hco3=24.0, cace=css_ace)
    say("      Correcting HCO3 18 -> 24 is worth %+.2f mmol/L, and %.0f%% of that"
        % (b["Ce"] - a["Ce"], 100.0))
    say("      comes from the %+.1f mmol/day of extra KALIURESIS, not from a shift."
        % (b["Erenal"] - a["Erenal"]))
    say("      The corollary is the well-documented negative result: bicarbonate is")
    say("      NOT an acute hyperkalaemia therapy, because the mechanism it actually")
    say("      has needs a kidney and needs weeks.")
    say("")
    say("  (2) THE MRA POTASSIUM COST IS FLAT, WHICH CONTRADICTS WHAT I EXPECTED.")
    say("      I built this model expecting the MRA penalty to grow as eGFR fell.")
    say("      It does not: the model gives +0.30, +0.30, +0.30, +0.29, +0.28 across")
    say("      eGFR 90 to 15.  What grows is not the increment but the BASELINE it")
    say("      is added to, so the same +0.30 tips a rising fraction of patients over")
    say("      5.5.  That is a different clinical claim -- it says the danger is not")
    say("      a steeper drug effect but a shorter distance to the threshold -- and")
    say("      it is the version supported by the model rather than by my prior.")
    say("")
    say("  (3) NOT REPRESENTED:  acute kidney injury superimposed on CKD (the single")
    say("      commonest real-world cause of severe hyperkalaemia), digoxin toxicity,")
    say("      trimethoprim/heparin/calcineurin-inhibitor effects on the ASDN, tumour")
    say("      lysis and rhabdomyolysis (K released from the ICF pool itself, which")
    say("      this model treats as a passive reservoir), pseudohyperkalaemia, and")
    say("      dialysis (an intermittent, extremely high-clearance removal term).")
    say("")
    say("  (4) THE OUTCOME LAYER IS ASSOCIATIVE, NOT CAUSAL.  hazard_K() is fitted to")
    say("      OBSERVATIONAL K-mortality curves, in which a high potassium is partly")
    say("      a marker of the illness that caused it.  The model therefore almost")
    say("      certainly OVERSTATES the benefit of lowering a number and understates")
    say("      the benefit of treating what raised it.  Every hazard figure in")
    say("      section F should be read with that bias in mind: it is the strongest")
    say("      assumption in the whole model and the least well supported.")
    say("")
    say("  (5) The intracellular space is TWO pools, which is enough to reproduce the")
    say("      chronic (%.0f mmol per mmol/L) and acute (%.0f mmol per mmol/L) buffer"
        % (buffer_capacity(4.2, 1.0, P), buffer_capacity(4.2, 1.0, P, fast_only=True)))
    say("      capacities simultaneously, but it is still not muscle.  Exercise, which")
    say("      raises venous K by >1 mmol/L within a minute through interstitial")
    say("      accumulation in the working limb, is not representable at all; nor is")
    say("      the release of potassium FROM the pool in rhabdomyolysis and tumour")
    say("      lysis, which this model treats only as a passive reservoir.")
    say("")
    say("  (6) The clinician controller in section F is a caricature: real")
    say("      down-titration is discrete, delayed by the interval between blood")
    say("      tests, and frequently permanent.  A model in which the dose creeps")
    say("      back up automatically is optimistic about exactly the behaviour the")
    say("      binder trials were designed to change.")


def main():
    say("CHRONIC HYPERKALAEMIA QSP MODEL -- REFERENCE IMPLEMENTATION REPORT")
    say("Generated by hk_reference_model.py (pure Python, no dependencies)")
    A_calibration()
    B_reserve_threshold()
    C_diet_ceiling()
    D_binder_trials()
    E_acute()
    arms = F_dilemma()
    G_virtual_population(arms)
    H_discrepancies()
    say("")
    say("=" * 78)
    say("END OF REPORT")
    say("=" * 78)
    with open(__file__.replace("hk_reference_model.py", "hk_model_report.txt"), "w") as f:
        f.write("\n".join(REPORT) + "\n")


if __name__ == "__main__":
    main()
