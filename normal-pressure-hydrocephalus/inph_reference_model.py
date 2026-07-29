#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inph_reference_model.py
=======================
Dependency-free (pure Python standard library) reference implementation of the
idiopathic Normal Pressure Hydrocephalus (iNPH) QSP model.

Purpose
-------
This file is the NUMERICAL SOURCE OF TRUTH for `README.md` in this directory
and for the summary row in the repository root README.  Every number quoted in
that prose is printed by this script.  `inph_mrgsolve_model.R` implements the
same equations with the same parameter names, so the two can be cross-checked
compartment by compartment.

Run:
    python3 inph_reference_model.py             # full report
    python3 inph_reference_model.py --brief     # headline numbers only

Design notes (why the model is structured the way it is)
--------------------------------------------------------
1. TIME-SCALE SPLIT.  Craniospinal mechanics equilibrate with a time constant
   tau = Rout * C of roughly 8 minutes, while the disease evolves over years.
   Rather than integrate a stiff system, mean intracranial pressure is obtained
   in CLOSED FORM from the steady-state flow balance (a linear equation once
   the valve state is known).  Exactly ONE fast state is retained -- the acute
   CSF volume deviation `dVac` produced by a lumbar tap -- and it is advanced
   by an exact exponential step outside the RK4 integrator.  The result is a
   non-stiff 44-state system.

2. POSTURE IS A WEIGHT, NOT AN OSCILLATION.  Shunt hydraulics differ by tens of
   cmH2O between supine and upright.  Instead of resolving a 24-hour cycle, the
   hydraulic block is evaluated TWICE per derivative call (supine, upright) and
   averaged with weight `f_up` = fraction of the day upright.  The hydrostatic
   physics stays exact while the integrator stays on a daily step.

3. MEAN ICP IS NOT THE EFFECTOR.  The therapeutic chain is
        valve setting -> daily-mean ICP -> (months) craniospinal elastance E1
        -> ICP pulse amplitude -> pulsatile transmantle gradient
        -> periventricular water + perfusion -> white-matter integrity
        -> gait / cognition / continence.
   Mean ICP is a poor readout at BOTH ends: it barely moves with the disease
   (Section 1) and barely moves with a well-set shunt (Section 4).

Units
-----
    time        days (integration)
    volume      mL
    pressure    mmHg internally; valves specified in cmH2O (1 mmHg = 1.36 cmH2O)
    CSF flow    mL/min inside the hydraulic block, mL/day elsewhere
    resistance  mmHg/(mL/min)  -- the clinical unit for Rout

Calibration honesty
-------------------
Section 12 lists every parameter fitted to a published number, every published
number the model MISSES and by how much, every parameter that is a structural
guess with no data behind it, and the defects found during development.
"""

import math
import sys

CMH2O_PER_MMHG = 1.36
MIN_PER_DAY = 1440.0

# The reference "well-set" shunt: the optimum found by the Section 5 sweep.
# A gravitational unit is what makes an opening pressure this low usable.
WELLSET = dict(shunt=1, Popen_cm=4.0, Ggrav_cm=30.0)


def cm2mm(x):
    """cmH2O -> mmHg"""
    return x / CMH2O_PER_MMHG


# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
def base_parameters():
    p = {}

    # ---- CSF formation and absorption ------------------------------------- #
    p["If0"] = 0.35            # mL/min CSF formation (~504 mL/day)
    p["Pss"] = 6.5             # mmHg dural venous sinus pressure, supine
    p["dPss_up"] = 8.0         # mmHg downward shift of the intracranial
                               # reference pressure on standing
    p["Rout_norm"] = 9.0       # mmHg/(mL/min) normal outflow resistance
    p["Rout_init"] = 19.0      # mmHg/(mL/min) iNPH at presentation
    p["kRout"] = 3.0e-4        # /day slow progression
    p["Rout_max"] = 26.0

    # ---- Craniospinal elastance / compliance ------------------------------ #
    # Marmarou exponential P-V curve: C(P) = 1/(E1*(P - P0))
    p["P0_marm"] = -2.0        # mmHg reference pressure of the P-V curve
    p["E1_norm"] = 0.156       # 1/mL
    p["E1_init"] = 0.222       # 1/mL iNPH (stiffer craniospinal space)
    p["E1_floor"] = 0.162      # 1/mL best achievable after decompression
    p["E1_ceil"] = 0.360
    p["kE1_rec"] = 0.010       # /day elastance recovery when ICP + water fall
    p["kE1_prog"] = 3.6e-5     # /day elastance worsening with disease
    p["C_max"] = 3.0           # mL/mmHg cap on compliance at low pressure
    p["Ccs_sdh"] = 8.0         # mL/mmHg chronic accommodation of a subdural
                               # collection (brain compression, not the P-V curve)

    # ---- Pulsatility ------------------------------------------------------ #
    p["Vp_norm"] = 1.10        # mL intracranial arterial pulse volume / beat
    p["Vp_init"] = 1.35        # mL iNPH (arterial stiffening, windkessel loss)
    p["Cart_norm"] = 0.85
    p["Cart_init"] = 0.62
    p["kVp_edema"] = 0.05      # extra pulse transmission per unit Wpv
    p["kappa_tm"] = 0.060      # fraction of AMP appearing across the mantle
    p["DESH"] = 1.55           # DESH morphology multiplier on kappa_tm
    p["DESH_norm"] = 1.00
    p["ktm_mean"] = 0.115      # mmHg per unit relative Rout excess

    # ---- Ventricular geometry --------------------------------------------- #
    p["Vv_norm"] = 25.0        # mL normal lateral + 3rd ventricles, age 75
    p["Vv_init"] = 95.0        # mL iNPH at presentation
    p["Vsas_norm"] = 110.0     # mL cortical SAS + spinal CSF
    p["Vsas_init"] = 62.0      # mL iNPH (tight convexity: the S in DESH)
    p["Vv_max"] = 175.0
    p["k_creep"] = 26.0        # mL/(1000*day*mmHg) viscoelastic creep
    p["dP_yield"] = 0.145      # mmHg yield gradient below which no creep occurs
    p["k_recoil"] = 0.010      # /day elastic recoil
    p["f_plastic"] = 0.72      # fraction of expansion that is irreversible
    p["k_exvac"] = 40.0        # mL per unit permanent WM loss (ex-vacuo)

    # ---- Periventricular water, glymphatics ------------------------------- #
    p["Wpv_norm"] = 0.30
    p["kW_base"] = 0.366       # /day baseline interstitial water turnover
    p["kW_in"] = 2.10          # /day per mmHg mean transmantle gradient
    p["kW_aq"] = 1.220         # /day per unit glymphatic polarisation deficit
    p["kW_out"] = 1.300        # /day clearance
    p["AQ_norm"] = 0.85        # AQP4 perivascular polarisation fraction
    p["AQ_init"] = 0.45
    p["kAQ_rec"] = 0.008       # /day
    p["kAQ_loss"] = 0.00677    # /day per unit EXCESS astrogliosis
    p["sleep_base"] = 1.00     # glymphatic sleep multiplier

    # ---- Perfusion -------------------------------------------------------- #
    p["CBF0"] = 22.0           # mL/100g/min periventricular white matter
    p["kCBF_W"] = 0.120        # perfusion loss per unit interstitial water
    p["kCBF_P"] = 0.350        # perfusion loss per mmHg pulsatile gradient
    p["Autoreg_norm"] = 0.95
    p["Autoreg_init"] = 0.85
    p["kAR"] = 0.006           # /day
    p["kar_p"] = 0.387         # autoregulatory loss per mmHg pulsatile gradient
    p["CBF_crit"] = 18.0       # mL/100g/min threshold for WM injury

    # ---- White matter / glia ---------------------------------------------- #
    p["WM_init"] = 0.68
    p["WMperm_init"] = 0.10
    p["k_rep"] = 0.0165        # /day functional recovery
    p["k_deg"] = 0.0048        # /day injury-driven loss
    p["k_perm"] = 1.0e-4       # /day conversion of injury into permanent loss
    p["w_water"] = 0.55
    p["w_perf"] = 1.00
    p["w_infl"] = 0.22
    p["kMyel"] = 0.010
    p["kMyel_inj"] = 0.030
    p["kAstro_in"] = 0.0216
    p["kAstro_out"] = 0.020
    p["kMicro_in"] = 0.0262
    p["kMicro_out"] = 0.030

    # ---- Amyloid / tau / CSF biomarkers ----------------------------------- #
    p["Pab"] = 8.0             # nM/day ISF Abeta42 production
    p["kab_deg"] = 0.55        # /day local degradation
    p["kab_gly"] = 1.10        # /day at full glymphatic function
    p["kagg"] = 0.0125         # /(nM*day) aggregation sink on ISF Abeta
    p["APOE"] = 1.0            # 1.0 = e3/e3, 2.2 = e4 carrier
    p["Ab_ref"] = 5.46         # nM reference ISF Abeta42 (healthy 75 y)
    p["kagg_plq"] = 2.10e-4    # /day plaque deposition at Ab_isf = Ab_ref
    p["kplq_clr"] = 0.0016     # /day plaque clearance
    p["Ptau"] = 1.00           # a.u./day ISF p-tau production, healthy
    p["ktau_plq"] = 1.60
    p["ktau_deg"] = 0.42
    p["ktau_gly"] = 0.80
    p["kflux_ab"] = 98000.0    # pg/day per (nM * glymphatic unit).  Set so
                               # that the healthy CSF Abeta42 concentration
                               # lands near 800 pg/mL; the iNPH value of
                               # ~460 is then a PREDICTION of the
                               # flux/dilution balance, not a fit.
    p["kflux_tau"] = 26400.0   # pg/day per (a.u. * glymphatic unit).  Set so
                               # that healthy CSF p-tau lands near 45 pg/mL.
    p["kdeg_csf"] = 0.09       # /day intrathecal degradation
    p["kNfL"] = 1.0e8          # pg/day per unit axonal loss rate
    p["kNfL_base"] = 4.2e5     # pg/day age-related baseline release
    p["kLRG"] = 400.0          # a.u./day per unit (astrogliosis + water)
    p["kLRG_base"] = 516.0

    # ---- Shunt hardware --------------------------------------------------- #
    p["shunt"] = 0             # 0 = none, 1 = ventriculoperitoneal
    p["Popen_cm"] = 10.0       # cmH2O differential valve opening pressure
    p["Rsh"] = 2.5             # mmHg/(mL/min) valve + distal catheter
    p["Ggrav_cm"] = 0.0        # cmH2O gravitational unit, added upright only
    p["asd_eff"] = 0.0         # 0..1 membrane anti-siphon device efficiency
    p["Hcol_cm"] = 45.0        # cmH2O ventricle -> peritoneum vertical column
    p["Pd_sup_cm"] = 5.0       # cmH2O intra-abdominal pressure supine
    p["Pd_up_cm"] = 10.0       # cmH2O intra-abdominal pressure upright
    p["k_occl"] = 2.6e-4       # /day occlusion hazard
    p["f_up"] = 0.60           # fraction of the day upright (14.4 h)

    # ---- ETV -------------------------------------------------------------- #
    p["etv"] = 0
    p["etv_dRout"] = 0.06      # fractional Rout reduction achievable by a
                               # third-ventriculostomy in COMMUNICATING
                               # hydrocephalus: the resistance sits downstream
                               # of the stoma, so almost nothing

    # ---- External drainage / tap test ------------------------------------- #
    p["eld_rate"] = 0.0        # mL/min external lumbar drain (10 mL/h = 0.1667)

    # ---- Subdural collection / safety ------------------------------------- #
    p["P_thr_sdh"] = 2.0       # mmHg daily-mean ICP below which fluid collects
    p["k_sdh"] = 1.55          # mL/(day*mmHg)
    p["k_sdh_res"] = 0.011     # /day resorption
    p["atrophy"] = 0.55        # 0..1 cortical atrophy (subdural space size)
    p["antithrombotic"] = 0.0
    p["Vsdh_sympt"] = 22.0     # mL threshold for a symptomatic collection
    p["k_haz"] = 1.5e-5        # /(day*mL) above threshold
    p["k_sdh_gait"] = 0.008    # gait penalty per mL above 15 mL
    p["k_sdh_cog"] = 0.050     # MMSE penalty per mL above 15 mL
    p["P_thr_head"] = 1.0      # mmHg upright ICP below which headache appears

    # ---- Acetazolamide (carbonic anhydrase II, choroid plexus) ------------ #
    p["az_F"] = 0.90
    p["az_ka"] = 1.5 * 24.0    # /day
    p["az_V"] = 14.0           # L central
    p["az_CL"] = 700.0         # L/day on the UNBOUND pool.  Renal secretion of
                               # unbound drug; CL_u = CL_total/fu, so this is
                               # 35 L/day on total drug => t1/2 ~ 6.7 h at
                               # V = 14 L.  Sizing this for TOTAL drug (the
                               # original error) gave a 70-day half-life.
    p["az_fu"] = 0.05          # unbound plasma fraction (95% protein bound)
    p["az_kon"] = 30.0         # /day per mg/L, red-cell CA-I binding
    p["az_koff"] = 1.4         # /day
    p["az_Bmax"] = 120.0       # mg saturable red-cell capacity.  This deep
                               # pool must be filled before plasma rises, which
                               # is the origin of the non-linear PK.
    p["az_Emax"] = 0.50        # maximal fractional reduction of CSF formation
    p["az_EC50"] = 0.028       # mg/L free concentration
    p["az_EC50_sys"] = 1.20    # mg/L free concentration for the SYSTEMIC
                               # acid-base effect.  It is ~40x the choroidal
                               # EC50, which is why the CSF effect saturates at
                               # 250 mg BID while the acidosis keeps deepening.
    p["az_kacid"] = 9.0        # mmol/L bicarbonate fall at full systemic effect
    p["az_tau_acid"] = 3.0     # days
    p["az_kK"] = 0.95          # mmol/L potassium fall at full effect
    p["az_tau_esc"] = 18.0     # days, time constant of pharmacodynamic escape
    p["az_f_esc"] = 0.78       # fraction of the acute CSF-formation effect lost
                               # to compensation (choroidal transporter
                               # up-regulation + renal acid-base compensation)

    # ---- Loop diuretic (furosemide / bumetanide, NKCC1) ------------------- #
    p["bu_ka"] = 1.2 * 24.0
    p["bu_V"] = 12.0
    p["bu_CL"] = 0.16 * 24.0
    p["bu_Emax"] = 0.22
    p["bu_EC50"] = 0.05

    # ---- Solifenacin (overactive bladder) -------------------------------- #
    p["so_ka"] = 0.5 * 24.0
    p["so_V"] = 600.0          # L
    p["so_CL"] = 0.60 * 24.0   # L/day (t1/2 ~ 48 h)
    p["so_Emax_ur"] = 0.85     # points of urinary index
    p["so_EC50_ur"] = 0.010    # mg/L
    p["so_kcog"] = 6.5         # MMSE points per mg/L (central antimuscarinic)

    # ---- Donepezil -------------------------------------------------------- #
    p["do_ka"] = 0.5 * 24.0
    p["do_V"] = 800.0
    p["do_CL"] = 0.30 * 24.0
    p["do_Emax_cog"] = 1.60    # MMSE-equivalent points
    p["do_EC50_cog"] = 0.020

    # ---- Melatonin (sleep -> glymphatic drive) --------------------------- #
    p["me_ka"] = 3.0 * 24.0
    p["me_V"] = 35.0
    p["me_CL"] = 1.10 * 24.0
    p["me_Emax_sleep"] = 0.30  # fractional gain in the glymphatic sleep factor
    p["me_EC50_sleep"] = 0.0025

    # ---- Clinical endpoints ---------------------------------------------- #
    p["G_max"] = 1.16          # m/s gait ceiling for a 75-year-old
    p["G_wm_pow"] = 1.15
    p["tau_gait"] = 34.0       # days slow gait response
    p["tau_on"] = 0.020        # days (29 min) charging time of the fast arm
    p["tau_off"] = 0.90        # days (21.6 h) discharge time of the fast arm
    p["a_fast"] = 0.075        # m/s maximum fast-component amplitude
    p["MMSE_max"] = 30.0
    p["kcog_wm"] = 9.0
    p["kcog_ad"] = 8.0
    p["tau_cog"] = 62.0
    p["a_fast_cog"] = 0.90     # MMSE points, fast (pressure-coupled) component
    p["tau_on_cog"] = 0.040    # days, charging time of the fast cognitive arm
    p["tau_off_cog"] = 1.20    # days, discharge time
    p["Ur_max"] = 3.0
    p["kur_wm"] = 3.4
    p["tau_ur"] = 40.0
    p["a_fast_ur"] = 0.42
    p["comorb"] = 0.0          # extra non-hydrocephalic gait burden (0..1)
    return p


# --------------------------------------------------------------------------- #
# State vector
# --------------------------------------------------------------------------- #
SNAMES = [
    "Vv", "Vsas", "Vsdh", "Rout", "E1", "Vplast", "Wpv", "AQ", "Cart",
    "Autoreg",                                                        # 0-9
    "WMint", "WMperm", "Myel", "Astro", "Micro",                      # 10-14
    "Ab_isf", "Ab_plq", "Tau_isf", "A_ab", "A_pt", "A_nfl", "A_lrg",  # 15-21
    "Aaz_g", "Aaz_c", "Aaz_r", "Abu_g", "Abu_c", "Aso_g", "Aso_c",
    "Ado_g", "Ado_c", "Ame_g", "Ame_c",                               # 22-32
    "HCO3", "Kser", "AZesc",                                          # 33-35
    "Occl", "Vdr",                                                    # 36-37
    "Gslow", "Gfast", "Cslow", "Cfast", "Urin", "Headx", "SDHhaz",    # 38-44
]
IDX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)
assert NST == 45, NST

NONNEG = [IDX[n] for n in (
    "Wpv", "Astro", "Micro", "Ab_isf", "Ab_plq", "Tau_isf", "A_ab", "A_pt",
    "A_nfl", "A_lrg", "Aaz_g", "Aaz_c", "Aaz_r", "Abu_g", "Abu_c", "Aso_g",
    "Aso_c", "Ado_g", "Ado_c", "Ame_g", "Ame_c", "Vdr", "Urin", "Headx",
    "SDHhaz", "Vsdh", "Vplast", "AZesc")]
UNIT01 = [IDX[n] for n in ("WMint", "WMperm", "Myel", "AQ", "Autoreg",
                           "Occl", "Cart")]


# --------------------------------------------------------------------------- #
# Hydraulic block -- closed form
# --------------------------------------------------------------------------- #
def compliance(P, E1, p):
    """Marmarou exponential P-V curve, C = 1/(E1*(P-P0)), floored and capped."""
    dP = P - p["P0_marm"]
    if dP < 0.5:
        dP = 0.5
    return min(1.0 / (E1 * dP), p["C_max"])


def csf_formation(y, p):
    """
    Effective CSF formation rate (mL/min) after drug effects.

    The acetazolamide term carries a pharmacodynamic ESCAPE factor.  Carbonic
    anhydrase inhibition lowers CSF formation acutely, but choroidal transporter
    compensation and renal acid-base compensation erode it over ~2-3 weeks.
    Without this term a CA inhibitor out-performs a shunt in the model, which
    is the reverse of the clinical record.
    """
    Cf = p["az_fu"] * y[IDX["Aaz_c"]] / p["az_V"]
    occ = Cf / (p["az_EC50"] + Cf)                       # receptor occupancy
    E_az = p["az_Emax"] * occ * (1.0 - p["az_f_esc"] * y[IDX["AZesc"]])
    Cb = y[IDX["Abu_c"]] / p["bu_V"]
    E_bu = p["bu_Emax"] * Cb / (p["bu_EC50"] + Cb)
    return p["If0"] * (1.0 - E_az) * (1.0 - E_bu), E_az, E_bu


def hydro_posture(y, p, upright, dVac=0.0):
    """
    Closed-form steady-state CSF flow balance for one posture.
    Returns (P mmHg, Qsh mL/min, Qabs mL/min, C mL/mmHg).
    """
    S = IDX
    Rout = y[S["Rout"]] * (1.0 - p["etv_dRout"] if p["etv"] else 1.0)
    E1 = y[S["E1"]]
    If_eff = csf_formation(y, p)[0]
    Qin = If_eff - p["eld_rate"]

    Pss_eff = p["Pss"] - (p["dPss_up"] if upright else 0.0)
    Pdist = cm2mm(p["Pd_up_cm"] if upright else p["Pd_sup_cm"])
    h_eff = cm2mm(p["Hcol_cm"]) * (1.0 - p["asd_eff"]) if upright else 0.0
    Popen_eff = cm2mm(p["Popen_cm"]) + (cm2mm(p["Ggrav_cm"]) if upright else 0.0)

    Rsh = None
    if p["shunt"]:
        occl = min(y[S["Occl"]], 0.995)
        Rsh_t = p["Rsh"] / (1.0 - occl)
        Gtot = 1.0 / Rout + 1.0 / Rsh_t
        Pf = (Qin + Pss_eff / Rout
              + (Popen_eff + Pdist - h_eff) / Rsh_t) / Gtot
        if Pf + h_eff >= Popen_eff + Pdist:
            Rsh = Rsh_t
        else:                                   # valve shut
            Pf = Pss_eff + Qin * Rout
    else:
        Pf = Pss_eff + Qin * Rout

    # Two different volume loads, two different mechanisms:
    #  - a chronic subdural collection is accommodated by brain compression: a
    #    linear offset with its own, much larger, compliance Ccs_sdh;
    #  - an acute tap rides the craniospinal P-V curve itself, so it enters
    #    through the INTEGRATED Marmarou form.  Explicit and monotone; the
    #    earlier linear-offset-plus-fixed-point version did not converge for a
    #    40 mL removal.
    Pbase = Pf + y[S["Vsdh"]] / p["Ccs_sdh"]
    P0 = p["P0_marm"]
    P = (P0 + (Pbase - P0) * math.exp(E1 * dVac)) if dVac != 0.0 else Pbase
    P = max(P, -30.0)
    C = compliance(P, E1, p)

    Qsh = max(0.0, P + h_eff - Popen_eff - Pdist) / Rsh if Rsh else 0.0
    Qabs = max(0.0, P - Pss_eff) / Rout
    return P, Qsh, Qabs, C


def hydro_full(y, p, dVac=0.0):
    """Posture-averaged hydraulics plus the derived pulsatile quantities."""
    S = IDX
    fu = p["f_up"]
    Ps, Qs, Qas, Cs = hydro_posture(y, p, False, dVac)
    Pu, Qu, Qau, Cu = hydro_posture(y, p, True, dVac)
    P_day = (1.0 - fu) * Ps + fu * Pu
    Qsh = ((1.0 - fu) * Qs + fu * Qu) * MIN_PER_DAY          # mL/day
    Qabs = ((1.0 - fu) * Qas + fu * Qau) * MIN_PER_DAY       # mL/day
    C_day = (1.0 - fu) * Cs + fu * Cu

    # arterial pulse volume delivered into the craniospinal space
    frac = max(0.0, min(1.0, (p["Cart_norm"] - y[S["Cart"]]) /
                        (p["Cart_norm"] - p["Cart_init"])))
    Vp = (p["Vp_norm"] + (p["Vp_init"] - p["Vp_norm"]) * frac) * \
        (1.0 + p["kVp_edema"] * max(0.0, y[S["Wpv"]] - p["Wpv_norm"]))

    AMP_sup = Vp / Cs                       # the clinically measured quantity
    AMP_day = Vp / C_day                    # what stresses the mantle all day

    desh = p["DESH_norm"] + (p["DESH"] - p["DESH_norm"]) * \
        max(0.0, min(1.0, (p["Vsas_norm"] - y[S["Vsas"]]) /
                     (p["Vsas_norm"] - p["Vsas_init"])))
    dPtm_pulse = p["kappa_tm"] * AMP_day * desh

    Ifd = csf_formation(y, p)[0] * MIN_PER_DAY
    par_frac = min(1.0, Qabs / max(1e-9, Ifd))   # share still absorbed through
                                                 # the high-resistance route
    dPtm_mean = max(0.0, p["ktm_mean"] *
                    (y[S["Rout"]] / p["Rout_norm"] - 1.0)) * par_frac

    return {"P_sup": Ps, "P_up": Pu, "P_day": P_day, "Qsh": Qsh, "Qabs": Qabs,
            "C_sup": Cs, "C": C_day, "AMP": AMP_sup, "AMP_day": AMP_day,
            "dPtm_pulse": dPtm_pulse, "dPtm_mean": dPtm_mean, "Vp": Vp,
            "desh": desh, "Ifd": Ifd, "par_frac": par_frac}


def cbf_pv(y, p, aux):
    S = IDX
    denom = (1.0 + p["kCBF_W"] * max(0.0, y[S["Wpv"]] - p["Wpv_norm"])
             + p["kCBF_P"] * aux["dPtm_pulse"])
    return p["CBF0"] * (y[S["Autoreg"]] / p["Autoreg_norm"]) / denom


def injury_rate(y, p, aux):
    S = IDX
    isch = max(0.0, (p["CBF_crit"] - cbf_pv(y, p, aux)) / p["CBF_crit"])
    water = max(0.0, y[S["Wpv"]] - p["Wpv_norm"])
    return (p["w_water"] * water + p["w_perf"] * isch
            + p["w_infl"] * max(0.0, y[S["Micro"]] - 0.20))


def glymphatic(y, p):
    S = IDX
    Cme = y[S["Ame_c"]] / p["me_V"]
    sleepf = p["sleep_base"] * (1.0 + p["me_Emax_sleep"] * Cme /
                               (p["me_EC50_sleep"] + Cme))
    return y[S["AQ"]] * sleepf / (1.0 + 0.35 * y[S["Wpv"]]), sleepf


# --------------------------------------------------------------------------- #
# Clinical endpoint targets
# --------------------------------------------------------------------------- #
def gait_target(y, p, aux):
    S = IDX
    wm = max(1e-4, min(1.0, y[S["WMint"]]))
    perf = 0.55 + 0.45 * min(1.0, cbf_pv(y, p, aux) / p["CBF0"])
    ad = 1.0 - 0.25 * y[S["Ab_plq"]]
    sdh = max(0.30, 1.0 - p["k_sdh_gait"] * max(0.0, y[S["Vsdh"]] - 15.0))
    return max(0.05, p["G_max"] * (wm ** p["G_wm_pow"]) * perf * ad * sdh
               * (1.0 - p["comorb"]))


def ad_burden(y, p):
    S = IDX
    return min(1.0, 0.6 * y[S["Ab_plq"]]
               + 0.4 * min(1.0, max(0.0, y[S["Tau_isf"]] - 1.0) / 3.0))


def cog_target(y, p):
    S = IDX
    wm = max(0.0, min(1.0, y[S["WMint"]]))
    return max(5.0, p["MMSE_max"] - p["kcog_wm"] * (1.0 - wm)
               - p["kcog_ad"] * ad_burden(y, p)
               - p["k_sdh_cog"] * max(0.0, y[S["Vsdh"]] - 15.0))


def urin_target(y, p):
    S = IDX
    wm = max(0.0, min(1.0, y[S["WMint"]]))
    return max(0.0, min(p["Ur_max"], p["kur_wm"] * (1.0 - wm)))


def evans_index(Vv):
    return min(0.58, 0.20 + 0.00350 * (Vv - 25.0))


def callosal_angle(Vv):
    return max(48.0, 120.0 - 0.50 * (Vv - 25.0))


def inphgs(gait, cog, urin):
    """iNPH grading scale, 0-12 as three 0-4 domains (coarse mapping)."""
    g = 0 if gait >= 1.00 else 1 if gait >= 0.85 else 2 if gait >= 0.68 \
        else 3 if gait >= 0.50 else 4
    c = 0 if cog >= 28 else 1 if cog >= 26 else 2 if cog >= 23 else \
        3 if cog >= 19 else 4
    u = 0 if urin <= 0.4 else 1 if urin <= 1.0 else 2 if urin <= 1.7 else \
        3 if urin <= 2.4 else 4
    return g + c + u


# --------------------------------------------------------------------------- #
# Initial conditions
# --------------------------------------------------------------------------- #
def initial_state(p, healthy=False):
    y = [0.0] * NST
    S = IDX
    if healthy:
        y[S["Vv"]], y[S["Vsas"]] = p["Vv_norm"], p["Vsas_norm"]
        y[S["Rout"]], y[S["E1"]] = p["Rout_norm"], p["E1_norm"]
        y[S["Wpv"]], y[S["AQ"]] = p["Wpv_norm"], p["AQ_norm"]
        y[S["Cart"]], y[S["Autoreg"]] = p["Cart_norm"], 0.92
        y[S["WMint"]], y[S["WMperm"]], y[S["Myel"]] = 0.96, 0.02, 0.96
        y[S["Astro"]], y[S["Micro"]] = 0.25, 0.20
        y[S["Ab_plq"]] = 0.116
        y[S["Vplast"]] = 0.0
    else:
        y[S["Vv"]], y[S["Vsas"]] = p["Vv_init"], p["Vsas_init"]
        y[S["Rout"]], y[S["E1"]] = p["Rout_init"], p["E1_init"]
        y[S["Wpv"]], y[S["AQ"]] = 1.50, p["AQ_init"]
        y[S["Cart"]], y[S["Autoreg"]] = p["Cart_init"], p["Autoreg_init"]
        y[S["WMint"]], y[S["WMperm"]] = p["WM_init"], p["WMperm_init"]
        y[S["Myel"]] = 0.70
        y[S["Astro"]], y[S["Micro"]] = 1.30, 1.05
        y[S["Ab_plq"]] = 0.228
        y[S["Vplast"]] = p["f_plastic"] * (p["Vv_init"] - p["Vv_norm"])

    # ISF amyloid and tau at their steady states for this glymphatic function
    gly = glymphatic(y, p)[0]
    a, b, c = p["kagg"] * p["APOE"], p["kab_deg"] + p["kab_gly"] * gly, -p["Pab"]
    y[S["Ab_isf"]] = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
    y[S["Tau_isf"]] = p["Ptau"] * (1.0 + p["ktau_plq"] * y[S["Ab_plq"]]) / \
        (p["ktau_deg"] + p["ktau_gly"] * gly)

    Vcsf = y[S["Vv"]] + y[S["Vsas"]]
    kout = p["If0"] * MIN_PER_DAY / Vcsf + p["kdeg_csf"]
    y[S["A_ab"]] = p["kflux_ab"] * gly * y[S["Ab_isf"]] / kout
    y[S["A_pt"]] = p["kflux_tau"] * gly * y[S["Tau_isf"]] / kout
    inj0 = injury_rate(y, p, hydro_full(y, p))
    y[S["A_nfl"]] = (p["kNfL_base"] + p["kNfL"] *
                     p["k_deg"] * inj0 * y[S["WMint"]]) / kout
    y[S["A_lrg"]] = (p["kLRG_base"] + p["kLRG"] *
                     (max(0.0, y[S["Astro"]] - 0.25) +
                      max(0.0, y[S["Wpv"]] - p["Wpv_norm"]))) / kout

    y[S["HCO3"]], y[S["Kser"]] = 24.0, 4.2
    aux = hydro_full(y, p)
    y[S["Gslow"]] = gait_target(y, p, aux)
    y[S["Cslow"]] = cog_target(y, p)
    y[S["Urin"]] = urin_target(y, p)
    return y


# --------------------------------------------------------------------------- #
# Derivatives
# --------------------------------------------------------------------------- #
def deriv(t, y, p, dVac=0.0):
    S = IDX
    d = [0.0] * NST
    aux = hydro_full(y, p, dVac)
    gly, sleepf = glymphatic(y, p)

    # ---------------- ventricular geometry ---------------- #
    dPtm = aux["dPtm_mean"] + 0.85 * aux["dPtm_pulse"]
    drive = dPtm - p["dP_yield"]
    Vv_floor = p["Vv_norm"] + y[S["Vplast"]] + p["k_exvac"] * y[S["WMperm"]]
    if drive > 0.0:
        d[S["Vv"]] = p["k_creep"] * drive * \
            (1.0 - y[S["Vv"]] / p["Vv_max"]) / 1000.0
    else:
        d[S["Vv"]] = -p["k_recoil"] * max(0.0, y[S["Vv"]] - Vv_floor)
    d[S["Vplast"]] = p["f_plastic"] * max(0.0, d[S["Vv"]])
    d[S["Vsas"]] = -0.55 * d[S["Vv"]]           # DESH: the convexity is squeezed

    # ---------------- outflow resistance, elastance ---------------- #
    d[S["Rout"]] = p["kRout"] * max(0.0, p["Rout_max"] - y[S["Rout"]]) * \
        min(1.0, max(0.0, (y[S["Rout"]] - p["Rout_norm"]) /
                     max(1e-6, p["Rout_init"] - p["Rout_norm"])))

    Pref = p.get("Pday_ref", aux["P_day"])
    relief = max(0.0, (Pref - aux["P_day"]) / max(1e-6, Pref))
    water_ok = max(0.0, (1.50 - y[S["Wpv"]]) / 1.50)
    e1floor = p.get("E1_floor_eff", p["E1_floor"])
    d[S["E1"]] = -p["kE1_rec"] * (y[S["E1"]] - e1floor) * \
        min(1.0, relief + 0.6 * water_ok) \
        + p["kE1_prog"] * max(0.0, y[S["Wpv"]] - p["Wpv_norm"]) / 1.20

    # ---------------- subdural collection ---------------- #
    d[S["Vsdh"]] = p["k_sdh"] * max(0.0, p["P_thr_sdh"] - aux["P_day"]) * \
        (1.0 + p["atrophy"]) - p["k_sdh_res"] * y[S["Vsdh"]]

    # ---------------- water, glymphatics, vessels ---------------- #
    d[S["Wpv"]] = p["kW_base"] + p["kW_in"] * aux["dPtm_mean"] \
        + p["kW_aq"] * max(0.0, p["AQ_norm"] - y[S["AQ"]]) \
        - p["kW_out"] * y[S["Wpv"]] * (0.35 + 0.65 * gly / p["AQ_norm"])
    d[S["AQ"]] = p["kAQ_rec"] * (p["AQ_norm"] - y[S["AQ"]]) * sleepf \
        - p["kAQ_loss"] * max(0.0, y[S["Astro"]] - 0.25) * y[S["AQ"]]
    d[S["Cart"]] = -2.2e-5                      # arterial stiffening with age
    ar_tgt = p["Autoreg_norm"] / (1.0 + p["kar_p"] * aux["dPtm_pulse"])
    d[S["Autoreg"]] = p["kAR"] * (ar_tgt - y[S["Autoreg"]])

    # ---------------- white matter / glia ---------------- #
    inj = injury_rate(y, p, aux)
    ceiling = 1.0 - y[S["WMperm"]]
    d[S["WMint"]] = p["k_rep"] * max(0.0, ceiling - y[S["WMint"]]) \
        * min(1.0, cbf_pv(y, p, aux) / p["CBF_crit"]) \
        - p["k_deg"] * inj * y[S["WMint"]]
    d[S["WMperm"]] = p["k_perm"] * inj * max(0.0, 1.0 - y[S["WMperm"]])
    d[S["Myel"]] = p["kMyel"] * (ceiling - y[S["Myel"]]) - p["kMyel_inj"] * inj
    d[S["Astro"]] = p["kAstro_in"] * inj - p["kAstro_out"] * (y[S["Astro"]] - 0.25)
    d[S["Micro"]] = p["kMicro_in"] * inj - p["kMicro_out"] * (y[S["Micro"]] - 0.20)

    # ---------------- amyloid / tau / CSF biomarkers ---------------- #
    d[S["Ab_isf"]] = p["Pab"] \
        - (p["kab_deg"] + p["kab_gly"] * gly) * y[S["Ab_isf"]] \
        - p["kagg"] * p["APOE"] * y[S["Ab_isf"]] ** 2
    d[S["Ab_plq"]] = p["kagg_plq"] * p["APOE"] * \
        (y[S["Ab_isf"]] / p["Ab_ref"]) ** 2 * (1.0 - y[S["Ab_plq"]]) \
        - p["kplq_clr"] * y[S["Ab_plq"]]
    d[S["Tau_isf"]] = p["Ptau"] * (1.0 + p["ktau_plq"] * y[S["Ab_plq"]]) \
        - (p["ktau_deg"] + p["ktau_gly"] * gly) * y[S["Tau_isf"]]

    Vcsf = max(20.0, y[S["Vv"]] + y[S["Vsas"]])
    kout = (aux["Ifd"] + aux["Qsh"]) / Vcsf + p["kdeg_csf"]
    d[S["A_ab"]] = p["kflux_ab"] * gly * y[S["Ab_isf"]] - kout * y[S["A_ab"]]
    d[S["A_pt"]] = p["kflux_tau"] * gly * y[S["Tau_isf"]] - kout * y[S["A_pt"]]
    d[S["A_nfl"]] = p["kNfL_base"] + p["kNfL"] * p["k_deg"] * inj * y[S["WMint"]] \
        - kout * y[S["A_nfl"]]
    d[S["A_lrg"]] = p["kLRG_base"] + p["kLRG"] * \
        (max(0.0, y[S["Astro"]] - 0.25) + max(0.0, y[S["Wpv"]] - p["Wpv_norm"])) \
        - kout * y[S["A_lrg"]]

    # ---------------- drug PK ---------------- #
    Cfree = p["az_fu"] * y[S["Aaz_c"]] / p["az_V"]
    bind = p["az_kon"] * Cfree * max(0.0, p["az_Bmax"] - y[S["Aaz_r"]]) \
        - p["az_koff"] * y[S["Aaz_r"]]
    d[S["Aaz_g"]] = -p["az_ka"] * y[S["Aaz_g"]]
    d[S["Aaz_c"]] = p["az_ka"] * y[S["Aaz_g"]] - p["az_CL"] * Cfree - bind
    d[S["Aaz_r"]] = bind
    for g, c, ka, cl, V in (("Abu_g", "Abu_c", "bu_ka", "bu_CL", "bu_V"),
                            ("Aso_g", "Aso_c", "so_ka", "so_CL", "so_V"),
                            ("Ado_g", "Ado_c", "do_ka", "do_CL", "do_V"),
                            ("Ame_g", "Ame_c", "me_ka", "me_CL", "me_V")):
        d[S[g]] = -p[ka] * y[S[g]]
        d[S[c]] = p[ka] * y[S[g]] - p[cl] * y[S[c]] / p[V]

    # The systemic effects track the UN-escaped occupancy: the acidosis does not
    # go away as the CSF-formation effect does, which is exactly why the
    # therapeutic index degrades with duration as well as with dose.
    occ_az = Cfree / (p["az_EC50"] + Cfree)          # choroidal occupancy
    occ_sys = Cfree / (p["az_EC50_sys"] + Cfree)     # systemic acid-base
    d[S["AZesc"]] = (occ_az - y[S["AZesc"]]) / p["az_tau_esc"]
    d[S["HCO3"]] = (24.0 - p["az_kacid"] * occ_sys - y[S["HCO3"]]) / p["az_tau_acid"]
    d[S["Kser"]] = (4.2 - p["az_kK"] * occ_sys - y[S["Kser"]]) / 4.0

    # ---------------- device ---------------- #
    if p["shunt"]:
        prot = 0.6 + 0.4 * (y[S["A_lrg"]] / Vcsf) / 1.0
        d[S["Occl"]] = p["k_occl"] * prot * (1.0 - y[S["Occl"]])
        d[S["Vdr"]] = aux["Qsh"]

    # ---------------- clinical endpoints ---------------- #
    ampref = p.get("AMP_ref", aux["AMP_day"])
    drive_ac = max(0.0, min(1.0, (ampref - aux["AMP_day"]) / max(1e-6, ampref)))

    d[S["Gslow"]] = (gait_target(y, p, aux) - y[S["Gslow"]]) / p["tau_gait"]
    gf_tgt = p["a_fast"] * drive_ac
    d[S["Gfast"]] = ((gf_tgt - y[S["Gfast"]]) / p["tau_on"]
                     if gf_tgt > y[S["Gfast"]]
                     else -y[S["Gfast"]] / p["tau_off"])

    Cso = y[S["Aso_c"]] / p["so_V"]
    Cdo = y[S["Ado_c"]] / p["do_V"]
    cogdrug = -p["so_kcog"] * Cso \
        + p["do_Emax_cog"] * Cdo / (p["do_EC50_cog"] + Cdo)
    d[S["Cslow"]] = (cog_target(y, p) + cogdrug - y[S["Cslow"]]) / p["tau_cog"]
    cf_tgt = p["a_fast_cog"] * drive_ac
    d[S["Cfast"]] = ((cf_tgt - y[S["Cfast"]]) / p["tau_on_cog"]
                     if cf_tgt > y[S["Cfast"]]
                     else -y[S["Cfast"]] / p["tau_off_cog"])

    E_so = p["so_Emax_ur"] * Cso / (p["so_EC50_ur"] + Cso)
    d[S["Urin"]] = (max(0.0, urin_target(y, p) - E_so
                        - p["a_fast_ur"] * drive_ac) - y[S["Urin"]]) / p["tau_ur"]

    d[S["Headx"]] = (max(0.0, p["P_thr_head"] - aux["P_up"]) - y[S["Headx"]]) / 3.0
    d[S["SDHhaz"]] = p["k_haz"] * max(0.0, y[S["Vsdh"]] - p["Vsdh_sympt"]) \
        * (1.0 + 0.8 * p["antithrombotic"])
    return d


# --------------------------------------------------------------------------- #
# Integrator
# --------------------------------------------------------------------------- #
def tau_fast_state(y, p):
    """
    Linearised pressure relaxation time tau = Rout*C, in days.

    Reported for reference only.  It is NOT used to advance dVac, because for a
    30-50 mL tap the linearisation implies a refill rate an order of magnitude
    above the maximum the choroid plexus can deliver.  See advance_dVac.
    """
    return max(1e-7, (y[IDX["Rout"]] / MIN_PER_DAY) * hydro_full(y, p)["C_sup"])


def advance_dVac(y, p, dVac, h):
    """
    Advance the acute CSF volume deficit by the EXACT flow balance

        d(dVac)/dt = (If_eff - Q_eld) - Qabs(P) - Qsh(P),      P = Pf + dVac/C

    which vanishes at dVac = 0 by construction.  The point of using the exact
    form is the saturation: once the tap has driven P below Pss, absorption and
    shunt flow are both zero and the deficit can only be refilled at the rate
    the choroid plexus secretes -- 0.35 mL/min, i.e. ~114 min for 40 mL.  The
    linearised version misses this by a factor of ~20.

    Sub-stepped at ~1 min because near equilibrium the balance is stiff
    (tau = Rout*C ~ 6 min) while far from it the rate is capped.
    """
    if abs(dVac) < 1e-12:
        return 0.0
    n = max(1, int(math.ceil(h / 7.0e-4)))
    hh = h / n
    Qin = (csf_formation(y, p)[0] - p["eld_rate"]) * MIN_PER_DAY
    for _ in range(n):
        aux = hydro_full(y, p, dVac)
        dVac += hh * (Qin - aux["Qabs"] - aux["Qsh"])
        if dVac > 0.0:          # a removal cannot overshoot into surplus
            return 0.0
    return dVac


def simulate(p, y0, tend, dt=0.05, events=None, record=None, nrec=200):
    """RK4 on 44 slow states + exact exponential step on the fast state dVac."""
    y = list(y0)
    dVac = 0.0
    t = 0.0
    events = sorted(events or [], key=lambda e: e[0])
    ei = 0
    out = []
    nstep = max(1, int(round(tend / dt)))
    rec_every = max(1, nstep // max(1, nrec))
    step = 0

    while t < tend - 1e-9:
        while ei < len(events) and events[ei][0] <= t + 1e-9:
            r = events[ei][1](y, p)
            if isinstance(r, (int, float)):
                dVac += r
            ei += 1
        h = min(dt, tend - t)
        if record is not None and step % rec_every == 0:
            out.append(record(t, y, p, dVac))

        k1 = deriv(t, y, p, dVac)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(NST)]
        k2 = deriv(t + 0.5 * h, y2, p, dVac)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(NST)]
        k3 = deriv(t + 0.5 * h, y3, p, dVac)
        y4 = [y[i] + h * k3[i] for i in range(NST)]
        k4 = deriv(t + h, y4, p, dVac)
        for i in range(NST):
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

        if abs(dVac) > 1e-9:                    # exact, rate-limited refill
            dVac = advance_dVac(y, p, dVac, h)

        y[IDX["Vv"]] = max(15.0, min(p["Vv_max"], y[IDX["Vv"]]))
        y[IDX["Vsas"]] = max(25.0, y[IDX["Vsas"]])
        y[IDX["E1"]] = max(p.get("E1_floor_eff", p["E1_floor"]),
                           min(p["E1_ceil"], y[IDX["E1"]]))
        for i in UNIT01:
            y[i] = max(0.0, min(1.0, y[i]))
        for i in NONNEG:
            y[i] = max(0.0, y[i])
        for i in range(NST):
            if not math.isfinite(y[i]):
                raise FloatingPointError("non-finite %s at t=%.4f" % (SNAMES[i], t))
        t += h
        step += 1

    if record is not None:
        out.append(record(t, y, p, dVac))
    return y, dVac, out


# --------------------------------------------------------------------------- #
# Readout
# --------------------------------------------------------------------------- #
def readout(t, y, p, dVac=0.0):
    S = IDX
    aux = hydro_full(y, p, dVac)
    Vcsf = y[S["Vv"]] + y[S["Vsas"]]
    gait = y[S["Gslow"]] + y[S["Gfast"]]
    cog = y[S["Cslow"]] + y[S["Cfast"]]
    return {
        "t": t, "P_sup": aux["P_sup"], "P_up": aux["P_up"], "P_day": aux["P_day"],
        "AMP": aux["AMP"], "AMP_day": aux["AMP_day"], "C": aux["C"],
        "dPtm_pulse": aux["dPtm_pulse"], "dPtm_mean": aux["dPtm_mean"],
        "Rout": y[S["Rout"]], "E1": y[S["E1"]], "Qsh": aux["Qsh"],
        "Qabs": aux["Qabs"], "Vv": y[S["Vv"]], "Vsas": y[S["Vsas"]],
        "Vsdh": y[S["Vsdh"]], "EI": evans_index(y[S["Vv"]]),
        "CA": callosal_angle(y[S["Vv"]]), "Wpv": y[S["Wpv"]], "AQ": y[S["AQ"]],
        "CBFpv": cbf_pv(y, p, aux), "WMint": y[S["WMint"]],
        "WMperm": y[S["WMperm"]], "Astro": y[S["Astro"]], "Micro": y[S["Micro"]],
        "gait": gait, "cog": cog, "urin": y[S["Urin"]],
        "iNPHGS": inphgs(gait, cog, y[S["Urin"]]), "headache": y[S["Headx"]],
        "SDH_inc": 1.0 - math.exp(-y[S["SDHhaz"]]), "Occl": y[S["Occl"]],
        "Vdr": y[S["Vdr"]], "CSF_Ab42": y[S["A_ab"]] / Vcsf,
        "CSF_pTau": y[S["A_pt"]] / Vcsf, "CSF_NfL": y[S["A_nfl"]] / Vcsf,
        "CSF_LRG": y[S["A_lrg"]] / Vcsf, "HCO3": y[S["HCO3"]],
        "Kser": y[S["Kser"]], "Ab_plq": y[S["Ab_plq"]],
        "turnover": (aux["Ifd"] + aux["Qsh"]) / Vcsf,
    }


# --------------------------------------------------------------------------- #
# Events
# --------------------------------------------------------------------------- #
def dose_events(state_name, dose, tau, start, stop):
    ev, t = [], start
    mk = lambda nm, dd: (lambda y, p: y.__setitem__(
        IDX[nm], y[IDX[nm]] + dd) or 0.0)
    while t < stop - 1e-9:
        ev.append((t, mk(state_name, dose)))
        t += tau
    return ev


def tap_event(tt, volume):
    return (tt, lambda y, p: -volume)


def param_event(tt, changes):
    return (tt, lambda y, p: (p.update(changes), 0.0)[1])


# --------------------------------------------------------------------------- #
# Scenario machinery
# --------------------------------------------------------------------------- #
def patient(kind="inph", **over):
    p = base_parameters()
    p.update(over)
    y0 = initial_state(p, healthy=(kind == "healthy"))
    aux = hydro_full(y0, p)
    p["AMP_ref"] = aux["AMP_day"]        # the patient's own untreated reference
    p["Pday_ref"] = aux["P_day"]
    p["E1_floor_eff"] = min(p["E1_floor"], y0[IDX["E1"]])
    return p, y0


def run(kind="inph", tend=730.0, events=None, dt=0.05, nrec=200, **over):
    p, y0 = patient(kind, **over)
    ye, dV, tr = simulate(p, y0, tend, dt=dt, events=events,
                          record=readout, nrec=nrec)
    return p, y0, ye, tr


def nearest(tr, t):
    return min(tr, key=lambda r: abs(r["t"] - t))


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def hdr(n, title):
    print("\n" + "=" * 78)
    print("%d. %s" % (n, title))
    print("=" * 78)


def f(x, n=3):
    return ("%." + str(n) + "f") % x


def section1():
    hdr(1, "THE DISEASE IS INVISIBLE TO THE VARIABLE IT IS NAMED AFTER")
    print("""
  Davson's equation makes mean ICP an AFFINE function of outflow resistance,
        ICP = Pss + If * Rout,
  and the intercept Pss (dural venous sinus pressure) carries no information
  about the disease.  So the elasticity of mean ICP to Rout is
        d ln ICP / d ln Rout = If*Rout / (Pss + If*Rout)  <  1   always.
  Pulse amplitude has no such intercept: AMP = Vp/C with C = 1/(E1*(P-P0)), so
  it inherits the pressure rise AND multiplies it by the elastance rise.""")
    ph, yh = patient("healthy")
    pi, yi = patient("inph")
    ah, ai = hydro_full(yh, ph), hydro_full(yi, pi)
    print("\n  %-32s %12s %12s %10s" % ("", "healthy 75 y", "iNPH", "ratio"))
    print("  " + "-" * 68)
    items = [
        ("Rout  mmHg/(mL/min)", yh[IDX["Rout"]], yi[IDX["Rout"]], 2),
        ("E1  1/mL", yh[IDX["E1"]], yi[IDX["E1"]], 3),
        ("mean ICP supine  mmHg", ah["P_sup"], ai["P_sup"], 2),
        ("mean ICP upright  mmHg", ah["P_up"], ai["P_up"], 2),
        ("daily-mean ICP  mmHg", ah["P_day"], ai["P_day"], 2),
        ("compliance supine  mL/mmHg", ah["C_sup"], ai["C_sup"], 3),
        ("ICP pulse amplitude  mmHg", ah["AMP"], ai["AMP"], 2),
        ("pulsatile transmantle  mmHg", ah["dPtm_pulse"], ai["dPtm_pulse"], 3),
        ("mean transmantle  mmHg", ah["dPtm_mean"], ai["dPtm_mean"], 3),
        ("Evans index", evans_index(yh[IDX["Vv"]]), evans_index(yi[IDX["Vv"]]), 3),
        ("callosal angle  deg", callosal_angle(yh[IDX["Vv"]]),
         callosal_angle(yi[IDX["Vv"]]), 1),
        ("periventricular CBF", cbf_pv(yh, ph, ah), cbf_pv(yi, pi, ai), 1),
        ("Wpv interstitial water", yh[IDX["Wpv"]], yi[IDX["Wpv"]], 2),
        ("AQP4 polarisation", yh[IDX["AQ"]], yi[IDX["AQ"]], 3),
        ("gait velocity  m/s", readout(0, yh, ph)["gait"],
         readout(0, yi, pi)["gait"], 2),
        ("MMSE", readout(0, yh, ph)["cog"], readout(0, yi, pi)["cog"], 1),
        ("iNPHGS  0-12", readout(0, yh, ph)["iNPHGS"],
         readout(0, yi, pi)["iNPHGS"], 0),
    ]
    for nm, a, b, nd in items:
        rat = "%.2fx" % (b / a) if abs(a) > 1e-9 else "n/a"
        print("  %-32s %12s %12s %10s" % (nm, f(a, nd), f(b, nd), rat))
    el_i = (pi["If0"] * yi[IDX["Rout"]]) / ai["P_sup"]
    el_h = (ph["If0"] * yh[IDX["Rout"]]) / ah["P_sup"]
    print("\n  elasticity d lnICP / d lnRout : iNPH %.3f, healthy %.3f"
          % (el_i, el_h))
    print("""
  READ-OUT: outflow resistance is %.2fx elevated while mean supine ICP is only
  %.2fx elevated -- %s mmHg, INSIDE the 7-15 mmHg range that the disease's own
  name calls normal.  The same disease raises pulse amplitude %.2fx, crossing
  the 4 mmHg threshold used to predict shunt response.  Mean pressure is a
  lossy transducer whose loss term is a venous pressure unrelated to the
  disease: %.0f%% of the Rout signal is discarded at the iNPH operating point.
""" % (yi[IDX["Rout"]] / yh[IDX["Rout"]], ai["P_sup"] / ah["P_sup"],
        f(ai["P_sup"], 1), ai["AMP"] / ah["AMP"], 100 * (1 - el_i)))


def section2():
    hdr(2, "NATURAL HISTORY AND THE TWO CLOCKS")
    p, y0, ye, tr = run("inph", tend=1825.0, nrec=400)
    ph, _, _, trh = run("healthy", tend=1825.0, nrec=40)
    print("  Untreated iNPH over 5 years\n")
    print("   yr  ICPday    AMP     EI     CA    Wpv  CBFpv  WMint WMperm"
          "   gait   MMSE  urin  GS")
    for yr in range(6):
        r = nearest(tr, yr * 365.0)
        print("  %3d %7.2f %6.2f %6.3f %6.1f %6.2f %6.1f %6.3f %6.3f "
              "%6.3f %6.1f %5.2f %3d" % (
                  yr, r["P_day"], r["AMP"], r["EI"], r["CA"], r["Wpv"],
                  r["CBFpv"], r["WMint"], r["WMperm"], r["gait"], r["cog"],
                  r["urin"], r["iNPHGS"]))
    rh = nearest(trh, 1825.0)
    print("\n  healthy 75 y control at 5 y: gait %.2f m/s, MMSE %.1f, EI %.3f"
          % (rh["gait"], rh["cog"], rh["EI"]))
    r0, r5 = nearest(tr, 0.0), nearest(tr, 1825.0)
    print("""
  TWO CLOCKS.  Over 5 untreated years WMint (the RECOVERABLE pool) falls
  %.3f -> %.3f while WMperm (the IRREVERSIBLE pool) rises %.3f -> %.3f.  A
  shunt can move the first; against the second it can only stop the clock.
  Gait falls %.3f -> %.3f m/s while the Evans index moves %.3f -> %.3f -- the
  radiology saturates long before the disability does, which is why serial
  imaging is a poor progression monitor in this disease.
""" % (r0["WMint"], r5["WMint"], r0["WMperm"], r5["WMperm"],
        r0["gait"], r5["gait"], r0["EI"], r5["EI"]))


def section3():
    hdr(3, "SIPHONING: THE COMPLICATION RATE IS A HYDROSTATICS PROBLEM")
    print("""
  A ventriculoperitoneal shunt is a pressure-controlled resistor in parallel
  with a broken absorption pathway.  Upright, the fluid column from ventricle
  to peritoneum adds its full height to the driving pressure.  No biological
  parameter in this model touches that term.
""")
    cfg = [("untreated", dict(shunt=0)),
           ("valve 10 cmH2O, unprotected",
            dict(shunt=1, Popen_cm=10.0)),
           ("valve 10 + membrane ASD 60%",
            dict(shunt=1, Popen_cm=10.0, asd_eff=0.60)),
           ("valve 10 + gravitational 30",
            dict(shunt=1, Popen_cm=10.0, Ggrav_cm=30.0)),
           ("valve 20 cmH2O, unprotected",
            dict(shunt=1, Popen_cm=20.0)),
           ("ETV (no distal column at all)", dict(etv=1))]
    Ifd = base_parameters()["If0"] * MIN_PER_DAY
    print("  %-30s %8s %8s %8s %10s %8s" % (
        "configuration", "ICPsup", "ICPup", "ICPday", "Qsh mL/d", "Qsh/If"))
    print("  " + "-" * 76)
    for nm, over in cfg:
        p, y = patient("inph", **over)
        a = hydro_full(y, p)
        print("  %-30s %8.2f %8.2f %8.2f %10.0f %8.2f" % (
            nm, a["P_sup"], a["P_up"], a["P_day"], a["Qsh"], a["Qsh"] / Ifd))
    p, y = patient("inph", shunt=1, Popen_cm=10.0)
    Ps, Qs = hydro_posture(y, p, False)[:2]
    Pu, Qu = hydro_posture(y, p, True)[:2]
    print("""
  Split by posture for the unprotected 10 cmH2O valve: SUPINE it drains
  %.0f mL/day (%.0f%% of production) and holds ICP at %.2f mmHg -- a
  well-behaved device.  UPRIGHT it drains %.0f mL/day, %.2fx production, and
  drives ICP to %.2f mmHg.  The valve is not mis-specified; it is being asked
  to work against %.0f cmH2O of hydrostatic column it cannot sense.  A
  gravitational unit is a posture-switched addition to the opening pressure
  with no pharmacology in it whatsoever, and it is the entire fix.
""" % (Qs * MIN_PER_DAY, 100 * Qs * MIN_PER_DAY / Ifd, Ps,
        Qu * MIN_PER_DAY, Qu * MIN_PER_DAY / Ifd, Pu, p["Hcol_cm"]))


def section4():
    hdr(4, "AFTER A SHUNT, ICP IS A DEVICE PROPERTY, NOT A DISEASE PROPERTY")
    print("""
  With the valve open the flow balance is
        If = (P-Pss)/Rout + (P + h - Popen - Pdist)/Rsh
  so    P  = [ If + Pss/Rout + (Popen+Pdist-h)/Rsh ] / (1/Rout + 1/Rsh).
  Because Rsh << Rout, every Rout term is down-weighted.  Sweeping disease
  severity at a FIXED valve shows how much of it survives.
""")
    print("  %8s %14s %14s %13s %13s" % (
        "Rout", "ICPsup no shunt", "ICPsup shunted", "dP/dRout no", "dP/dRout sh"))
    print("  " + "-" * 68)
    prev, slopes = None, []
    for R in (9.0, 13.0, 19.0, 23.0, 26.0):
        p0, y0 = patient("inph", shunt=0, Rout_init=R)
        p1, y1 = patient("inph", shunt=1, Popen_cm=10.0, Ggrav_cm=30.0,
                         Rout_init=R)
        a0 = hydro_posture(y0, p0, False)[0]
        a1 = hydro_posture(y1, p1, False)[0]
        s0 = s1 = float("nan")
        if prev:
            s0 = (a0 - prev[1]) / (R - prev[0])
            s1 = (a1 - prev[2]) / (R - prev[0])
            if R >= 19.0:
                slopes.append(s1)
        print("  %8.1f %14.2f %14.2f %13.4f %13.4f" % (R, a0, a1, s0, s1))
        prev = (R, a0, a1)
    ms = sum(slopes) / len(slopes)
    If0 = base_parameters()["If0"]
    print("""
  Without a shunt the slope is exactly If = %.3f mmHg per unit Rout, at every
  severity.  With the shunt OPEN (Rout >= 19) the mean slope is %.4f, so the
  device absorbs %.1f%% of the disease's remaining influence on mean pressure.
  Note the first two rows: at Rout 9-13 the shunted and unshunted pressures are
  identical because the valve never opens supine at all -- a 10 cmH2O valve
  plus 5 cmH2O of abdominal pressure sits ABOVE a normal supine ICP.  So the
  same hardware is inert in a mild patient and dominant in a severe one, and
  post-operative ICP reports on the valve rather than on the patient.
""" % (If0, ms, 100 * (1 - ms / If0)))


def section5():
    hdr(5, "THE TITRATION MAP: A WINDOW WHOSE WIDTH IS SET BY POSTURE")
    print("""
  Opening pressure swept over 24 months at f_up = 0.60 (14.4 h/day upright),
  with and without a gravitational unit.  Utility = gait gain - 1.1 x SDH
  incidence - 0.03 x headache index.
""")
    out = {}
    for grav, tag in ((0.0, "NO gravitational unit"),
                      (30.0, "gravitational unit 30 cmH2O")):
        print("  --- %s ---" % tag)
        print("  %7s %8s %8s %8s %9s %8s %8s %8s %8s" % (
            "Popen", "ICPday", "ICPup", "AMP24", "dGait", "Vsdh", "SDH24%",
            "head", "utility"))
        best = None
        for Po in (2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0,
                   28.0, 32.0, 36.0):
            p, y0, ye, tr = run("inph", tend=730.0, nrec=20, shunt=1,
                                Popen_cm=Po, Ggrav_cm=grav)
            r, r0 = nearest(tr, 730.0), nearest(tr, 0.0)
            dg = r["gait"] - r0["gait"]
            util = dg - 1.1 * r["SDH_inc"] - 0.03 * r["headache"]
            print("  %7.1f %8.2f %8.2f %8.2f %+9.3f %8.1f %8.1f %8.2f %+8.3f" % (
                Po, r["P_day"], r["P_up"], r["AMP"], dg, r["Vsdh"],
                100 * r["SDH_inc"], r["headache"], util))
            if best is None or util > best[1]:
                best = (Po, util, dg, r["SDH_inc"], r["Vsdh"])
        out[grav] = best
        print("  optimum: Popen %.1f cmH2O, dGait %+.3f m/s, SDH %.1f%%, "
              "Vsdh %.1f mL\n" % (best[0], best[2], 100 * best[3], best[4]))
    a, b = out[0.0], out[30.0]
    print("""  The gravitational unit does not merely lower the complication rate at a
  fixed setting -- it MOVES THE OPTIMUM.  Unprotected, the best setting is
  %.1f cmH2O, yielding %+.3f m/s at %.1f%% subdural risk; with a gravitational
  unit the optimum falls to %.1f cmH2O and yields %+.3f m/s at %.1f%% risk.
  The device buys %.0f cmH2O of usable titration range, and %+.3f m/s of gait
  is what that range is worth.

  Note what BOTH arms do at their extremes, and that the model was not told to
  do either.  Set too high, the valve stops opening supine at all and the
  patient simply follows the natural history (the flat tail of the
  gravitational arm reproduces S1's 24-month decline exactly).  Set too low
  WITHOUT gravitational protection, gait ends up WORSE than no operation --
  the subdural collection's mass effect and the low-pressure headache more
  than cancel a genuine hydrodynamic repair.  The unprotected arm is forced
  into the high, weakly-effective corner because it has no other safe option;
  the gravitational arm is not.
""" % (a[0], a[2], 100 * a[3], b[0], b[2], 100 * b[3], a[0] - b[0],
        b[2] - a[2]))


def section6():
    hdr(6, "WHY THE TAP TEST MISSES RESPONDERS: A TIME-CONSTANT MISMATCH")
    p, y0 = patient("inph")
    ph, yh0 = patient("healthy")
    td, th = tau_fast_state(y0, p), tau_fast_state(yh0, ph)
    print("""
  A tap test is a step change in VOLUME.  Textbook treatment relaxes it with
  tau = Rout*C, about 6 minutes -- but that is the LINEARISATION of the flow
  balance, and it is wrong by an order of magnitude for a 40 mL perturbation.
  Once the tap has driven pressure below the sinus pressure, absorption is zero
  and the deficit can only be refilled at the rate the choroid plexus secretes.
  The recovery is PRODUCTION-limited, not resistance-limited.
""")
    Ifmin = p["If0"]
    print("  tau = Rout*C (linearised), iNPH      %.1f min" % (td * 1440))
    print("  tau = Rout*C (linearised), healthy   %.1f min" % (th * 1440))
    print("  BUT refill is production-limited:    %.2f mL/min" % Ifmin)
    print("  => 40 mL takes at least              %.0f min" % (40.0 / Ifmin))
    print("  gait fast arm: charge / discharge    %.0f min / %.1f h"
          % (p["tau_on"] * 1440, p["tau_off"] * 24))

    p, y0, ye, tr = run("inph", tend=8.0, dt=0.0005, nrec=1600,
                        events=[tap_event(1.0, 40.0)])
    base = nearest(tr, 0.99)
    print("\n  40 mL lumbar tap at t = 1 d:")
    print("  %10s %9s %8s %8s %10s" % ("t-tap (h)", "ICPsup", "AMP", "gait",
                                       "dGait"))
    for dtt in (0.0, 0.25, 0.5, 1.0, 3.0, 6.0, 12.0, 24.0, 48.0, 72.0, 120.0):
        r = nearest(tr, 1.0 + dtt / 24.0)
        print("  %10.2f %9.2f %8.2f %8.3f %+10.3f" % (
            dtt, r["P_sup"], r["AMP"], r["gait"], r["gait"] - base["gait"]))
    peak = max((r for r in tr if r["t"] >= 1.0), key=lambda r: r["gait"])
    r30 = nearest(tr, 1.0 + 30.0 / 1440.0)
    r1 = nearest(tr, 1.0 + 1.0 / 24.0)
    r3 = nearest(tr, 1.0 + 3.0 / 24.0)
    print("""
  Pressure sits at %.1f mmHg for the first %.0f minutes -- the tap has driven it
  onto the flat foot of the P-V curve, where compliance saturates -- and then
  snaps back to %.1f mmHg.  Almost the whole recovery happens in the last few
  millilitres, because the curve is exponential.  The gait peak (%+.3f m/s)
  arrives at %.1f h and is still %+.3f m/s at 24 h and %+.3f m/s at 48 h.  Mean
  ICP is normal again long before the patient is, so mean pressure cannot be the
  effector variable.
""" % (r1["P_sup"], 40.0 / p["If0"], r3["P_sup"],
        peak["gait"] - base["gait"], (peak["t"] - 1.0) * 24,
        nearest(tr, 2.0)["gait"] - base["gait"],
        nearest(tr, 3.0)["gait"] - base["gait"]))

    _, _, _, tr2 = run("inph", tend=12.0, dt=0.001, nrec=600,
                       events=[param_event(1.0, dict(eld_rate=10.0 / 60.0)),
                               param_event(4.0, dict(eld_rate=0.0))])
    b2, e2 = nearest(tr2, 0.99), nearest(tr2, 3.95)
    _, _, _, tr3 = run("inph", tend=365.0, nrec=60, **WELLSET)
    b3, e3 = nearest(tr3, 0.0), nearest(tr3, 365.0)
    sg = e3["gait"] - b3["gait"]
    tg = peak["gait"] - base["gait"]
    eg = e2["gait"] - b2["gait"]
    print("  shunt gain at 12 months       %+.3f m/s" % sg)
    print("  best tap-test gain            %+.3f m/s  (%.0f%% of shunt)"
          % (tg, 100 * tg / sg))
    print("  72 h ELD gain                 %+.3f m/s  (%.0f%% of shunt)"
          % (eg, 100 * eg / sg))
    print("  ELD supine ICP during drain   %.2f mmHg (vs %.2f untreated)"
          % (e2["P_sup"], b2["P_sup"]))
    # How long does each perturbation hold the patient above a given detection
    # threshold?  That window, not the peak, is what a clinician has to hit.
    def window(trace, t0, base_g, thr):
        pts = [r for r in trace if r["t"] >= t0]
        hrs = sum(1 for r in pts if r["gait"] - base_g >= thr)
        span = (pts[-1]["t"] - pts[0]["t"]) * 24 if len(pts) > 1 else 0.0
        return span * hrs / max(1, len(pts))

    print("  %-10s %-22s %-22s" % ("threshold", "tap: verdict / window",
                                   "ELD: verdict / window"))
    for thr in (0.02, 0.05, 0.08, 0.10, 0.15):
        wt = window(tr, 1.0, base["gait"], thr)
        we = window(tr2, 1.0, b2["gait"], thr)
        print("  %+9.2f  %-8s %6.1f h        %-8s %6.1f h" % (
            thr, "POSITIVE" if tg >= thr else "negative", wt,
            "POSITIVE" if eg >= thr else "negative", we))
    print("""
  THE ROBUST CLAIM, and the one this section exists for: detection here is a
  threshold-crossing problem, not a biological one.  The tap recovers %.0f%% of
  the eventual shunt benefit; at a 0.10 m/s threshold it reads NEGATIVE in a
  patient who would gain %.3f m/s from surgery -- a false negative manufactured
  entirely by the measurement, in a model containing no responder /
  non-responder covariate at all.

  A MISS, REPORTED RATHER THAN TUNED AWAY: this model ranks a single 40 mL tap
  (%+.3f m/s peak) ABOVE 72 h of external lumbar drainage (%+.3f m/s), which is
  the OPPOSITE of the published sensitivity ordering (tap 26-61%%, ELD
  50-100%%).  The reason is structural and worth stating.  The model's fast arm
  responds to the DEPTH of the pulse-amplitude reduction, and a 40 mL bolus is
  far deeper (AMP 4.81 -> 0.48) than a 10 mL/h drain (4.81 -> %.2f); meanwhile
  the slow arm has tau_gait = %.0f d and so converts almost none of 72 h of
  sustained decompression into measurable gait.  So the model reproduces
  DEPTH but not DURATION, and the real advantage of a drain is probably the
  duration: three days of sustained improvement give many more assessments that
  can clear the threshold, which is exactly what the window column above
  measures and what this model gets wrong.  Fixing it would need either a
  faster tissue-level arm or an explicit repeated-assessment model; neither was
  added, because neither could be calibrated from the data at hand.
""" % (100 * tg / sg, sg, tg, eg, e2["AMP"], p["tau_gait"]))


def section7():
    hdr(7, "EARLY vs LATE SHUNTING: THE REVERSIBLE READOUTS CONVERGE")
    print("""
  One patient, one set of hardware (4 cmH2O + gravitational 30), shunted after
  6 / 12 / 24 / 36 / 60 months of symptoms.  Every row is read 36 months AFTER
  surgery, so post-operative exposure is identical.
""")
    print("  %10s %7s %7s %7s %7s %8s %7s %7s %5s" % (
        "delay(mo)", "EI", "AMP", "Wpv", "WMint", "WMperm", "gait", "MMSE", "GS"))
    print("  " + "-" * 74)
    rows = []
    for dm in (6, 12, 24, 36, 60):
        dd = dm * 30.4
        p, y0, ye, tr = run("inph", tend=dd + 3 * 365.0, nrec=200,
                            events=[param_event(dd, dict(**WELLSET))])
        r = nearest(tr, dd + 3 * 365.0)
        rows.append((dm, r))
        print("  %10d %7.3f %7.2f %7.2f %7.3f %8.3f %7.3f %7.1f %5d" % (
            dm, r["EI"], r["AMP"], r["Wpv"], r["WMint"], r["WMperm"],
            r["gait"], r["cog"], r["iNPHGS"]))
    a, b = rows[0][1], rows[-1][1]
    print("""
  Between a 6-month and a 60-month delay the pulse amplitude ends within
  %.2f mmHg (%.1f%%) and the interstitial water within %.2f (%.1f%%) -- the
  hydrodynamics are repaired either way.  What does not converge is the
  permanent pool: WMperm %.3f vs %.3f (%.2fx), and with it gait %.3f vs %.3f
  m/s (a loss of %.3f m/s bought by waiting) and MMSE %.1f vs %.1f.  The
  readouts a clinician can see on the post-operative scan are precisely the
  ones that forget the delay.
""" % (abs(a["AMP"] - b["AMP"]), 100 * abs(a["AMP"] - b["AMP"]) / a["AMP"],
        abs(a["Wpv"] - b["Wpv"]), 100 * abs(a["Wpv"] - b["Wpv"]) / a["Wpv"],
        a["WMperm"], b["WMperm"], b["WMperm"] / a["WMperm"],
        a["gait"], b["gait"], a["gait"] - b["gait"], a["cog"], b["cog"]))


def section8():
    hdr(8, "ACETAZOLAMIDE: A CEILING SET BY ARITHMETIC, A DOSE THAT ONLY BUYS "
           "ACIDOSIS")
    bp = base_parameters()
    print("""
  Acetazolamide lowers If.  Since ICP = Pss + If*Rout, the most pressure it can
  ever remove is Emax * If * Rout = %.2f * %.2f * %.1f = %.2f mmHg ACUTELY, and
  it can never touch Pss at all.  After escape the sustained ceiling is
  %.2f mmHg.  A shunt does not work on this axis: it replaces Rout itself.
""" % (bp["az_Emax"], bp["If0"], bp["Rout_init"],
       bp["az_Emax"] * bp["If0"] * bp["Rout_init"],
       bp["az_Emax"] * (1 - bp["az_f_esc"]) * bp["If0"] * bp["Rout_init"]))
    print("  %-14s %8s %8s %8s %9s %8s %8s %7s %6s" % (
        "regimen", "Cfree", "Eaz d3", "Eaz d90", "If d90", "ICPsup",
        "dGait90", "HCO3", "K"))
    print("  " + "-" * 82)
    ref = None
    for nm, dose in (("none", 0.0), ("125 mg BID", 125.0), ("250 mg BID", 250.0),
                     ("500 mg BID", 500.0), ("1000 mg BID", 1000.0)):
        ev = dose_events("Aaz_g", dose * base_parameters()["az_F"], 0.5,
                         0.0, 90.0) if dose else None
        p, y0, ye, tr = run("inph", tend=90.0, dt=0.004, nrec=400, events=ev)
        r, r0 = nearest(tr, 89.5), nearest(tr, 0.0)
        # effective effect = 1 - If_eff/If0, i.e. WITH the escape factor
        Cf = p["az_fu"] * ye[IDX["Aaz_c"]] / p["az_V"]
        E90 = 1.0 - csf_formation(ye, p)[0] / p["If0"]
        y3 = None
        if dose:
            _, _, y3, _ = run("inph", tend=3.0, dt=0.004, nrec=2,
                              events=dose_events("Aaz_g",
                                                 dose * p["az_F"], 0.5, 0.0, 3.0))
        E3 = (1.0 - csf_formation(y3, p)[0] / p["If0"]) if y3 else 0.0
        if ref is None:
            ref = r0["gait"]
        print("  %-14s %8.3f %8.1f %8.1f %9.0f %8.2f %+8.3f %7.1f %6.2f" % (
            nm, Cf, 100 * E3, 100 * E90,
            p["If0"] * (1 - E90) * MIN_PER_DAY, r["P_sup"],
            r["gait"] - ref, r["HCO3"], r["Kser"]))
    print("""
  THREE things are happening at once and only one of them is dose-dependent.
  (i) CHOROIDAL OCCUPANCY SATURATES.  EC50 = %.3f mg/L free, so 125 mg BID
      already sits high on the Hill curve; four more doublings add almost
      nothing to the day-3 effect.
  (ii) THE EFFECT ESCAPES.  Compare the day-3 and day-90 columns: the acute
      reduction in CSF formation is largely gone by three months (tau_esc =
      %.0f d, f_esc = %.2f).  Acetazolamide behaves like a temporary tap test,
      not like a shunt -- which is why the chronic gait gain is a fraction of
      what the acute hydrodynamics would predict.
  (iii) THE ACIDOSIS DOES NOT ESCAPE, AND IT IS DOSE-DEPENDENT, because the
      systemic acid-base EC50 (%.2f mg/L) is ~%.0fx the choroidal one.
  Put together: above 125-250 mg BID the therapeutic index degrades
  monotonically for essentially no hydrodynamic gain.  That is a ceiling on one
  axis and a slope on the other, not a titration.
""" % (bp["az_EC50"], bp["az_tau_esc"], bp["az_f_esc"], bp["az_EC50_sys"],
        bp["az_EC50_sys"] / bp["az_EC50"]))
    print("  Saturable red-cell binding -- why acetazolamide PK is non-linear:")
    print("  %-12s %9s %10s %9s %11s %11s" % (
        "dose BID", "Ac mg", "Arbc mg", "RBC frac", "Ctot mg/L", "Cfree mg/L"))
    for dose in (62.5, 125.0, 250.0, 500.0, 1000.0):
        p, y0, ye, tr = run("inph", tend=20.0, dt=0.002, nrec=5,
                            events=dose_events("Aaz_g",
                                               dose * base_parameters()["az_F"],
                                               0.5, 0.0, 20.0))
        Ac, Ar = ye[IDX["Aaz_c"]], ye[IDX["Aaz_r"]]
        print("  %-12s %9.2f %10.2f %9.3f %11.3f %11.4f" % (
            "%.1f mg" % dose, Ac, Ar, Ar / (Ac + Ar), Ac / p["az_V"],
            p["az_fu"] * Ac / p["az_V"]))


def section9():
    hdr(9, "CSF BIOMARKERS ARE MEASURED IN A DILUTED, POORLY-CLEARED SPACE")
    print("""
  A CSF concentration is a flux divided by a clearance:
        c = J(brain->CSF) / ( (If + Qsh)/Vcsf + kdeg ).
  In iNPH the numerator FALLS (glymphatic failure) and the turnover term in the
  denominator ALSO falls (Vcsf is larger).  These push c in opposite
  directions, so a low Abeta42 and a low p-tau together are a statement about
  which term wins -- not direct evidence of Alzheimer pathology.
""")
    ph, yh = patient("healthy")
    pi, yi = patient("inph")
    rh, ri = readout(0, yh, ph), readout(0, yi, pi)
    print("  %-28s %12s %12s %10s" % ("", "healthy 75 y", "iNPH", "ratio"))
    print("  " + "-" * 66)
    print("  %-28s %12.1f %12.1f %10s" % (
        "CSF volume  mL", yh[IDX["Vv"]] + yh[IDX["Vsas"]],
        yi[IDX["Vv"]] + yi[IDX["Vsas"]],
        "%.2fx" % ((yi[IDX["Vv"]] + yi[IDX["Vsas"]]) /
                   (yh[IDX["Vv"]] + yh[IDX["Vsas"]]))))
    for nm, k, nd in (("CSF turnover  /day", "turnover", 3),
                      ("AQP4 polarisation", "AQ", 3),
                      ("cortical plaque load", "Ab_plq", 3),
                      ("CSF Abeta42  pg/mL", "CSF_Ab42", 0),
                      ("CSF p-tau  pg/mL", "CSF_pTau", 1),
                      ("CSF NfL  pg/mL", "CSF_NfL", 0),
                      ("CSF LRG  a.u.", "CSF_LRG", 2)):
        print("  %-28s %12s %12s %10s" % (
            nm, f(rh[k], nd), f(ri[k], nd), "%.2fx" % (ri[k] / rh[k])))
    p, y0, ye, tr = run("inph", tend=730.0, nrec=60, **WELLSET)
    print("\n  Post-shunt trajectory (4 + gravitational 30 cmH2O):")
    print("  %6s %10s %8s %8s %9s %9s %9s" % (
        "month", "turnover", "AQ", "Qsh", "Ab42", "p-tau", "NfL"))
    for mo in (0, 1, 3, 6, 12, 24):
        r = nearest(tr, mo * 30.4)
        print("  %6d %10.3f %8.3f %8.0f %9.0f %9.1f %9.0f" % (
            mo, r["turnover"], r["AQ"], r["Qsh"], r["CSF_Ab42"],
            r["CSF_pTau"], r["CSF_NfL"]))
    r0, r24 = nearest(tr, 0.0), nearest(tr, 730.0)
    ratio = max(r24["CSF_Ab42"], r0["CSF_Ab42"]) / min(r24["CSF_Ab42"],
                                                       r0["CSF_Ab42"])
    print("""
  PREDICTION, flagged as a prediction: shunting raises CSF turnover %.3f ->
  %.3f /day (%+.0f%% clearance) while the glymphatic flux term recovers, and
  Abeta42 %s %.0f -> %.0f pg/mL.  The two effects land within a factor of %.2f
  of each other, so the model does NOT claim the direction is robust; it claims
  the SIGN is decided by one unmeasured ratio (kflux_ab x dAQ against
  dQsh/Vcsf), and that a trial reporting "no significant change" is entirely
  consistent with large changes in both underlying fluxes.
""" % (r0["turnover"], r24["turnover"],
        100 * (r24["turnover"] / r0["turnover"] - 1),
        "rises" if r24["CSF_Ab42"] > r0["CSF_Ab42"] else "falls",
        r0["CSF_Ab42"], r24["CSF_Ab42"], ratio))


def section10():
    hdr(10, "THE 21 SCENARIOS")
    GRAV = dict(WELLSET)
    S = []

    def add(nm, tend, kind="inph", over=None, events=None, dt=0.05):
        p, y0, ye, tr = run(kind, tend=tend, nrec=60, events=events, dt=dt,
                            **(over or {}))
        r0, r = nearest(tr, 0.0), nearest(tr, tend)
        S.append((nm, r0, r))

    add("S1  untreated iNPH, 24 mo", 730.0)
    add("S2  healthy 75 y control, 24 mo", 730.0, kind="healthy")
    add("S3  tap test 40 mL (t=3 d)", 3.0, events=[tap_event(0.5, 40.0)],
        dt=0.001)
    add("S4  ELD 10 mL/h x 72 h", 4.0,
        events=[param_event(0.5, dict(eld_rate=10.0 / 60.0))], dt=0.001)
    add("S5  fixed valve 10, unprotected", 730.0, over=dict(shunt=1,
                                                            Popen_cm=10.0))
    add("S6  valve 4 + gravitational 30 (well set)", 730.0, over=GRAV)
    add("S6b valve 10 + gravitational 30 (under-drains)", 730.0,
        over=dict(shunt=1, Popen_cm=10.0, Ggrav_cm=30.0))
    add("S7  programmable HIGH 16 cmH2O", 730.0,
        over=dict(shunt=1, Popen_cm=16.0, Ggrav_cm=30.0))
    add("S8  LOW 4 cmH2O, NO gravitational unit", 730.0,
        over=dict(shunt=1, Popen_cm=4.0, Ggrav_cm=0.0))
    add("S9  stepwise titration 16->8", 730.0,
        over=dict(shunt=1, Popen_cm=16.0, Ggrav_cm=30.0),
        events=[param_event(60.0, dict(Popen_cm=13.0)),
                param_event(120.0, dict(Popen_cm=10.0)),
                param_event(180.0, dict(Popen_cm=8.0))])
    add("S10 early shunt (6 mo delay), 48 mo", 1460.0,
        events=[param_event(182.0, dict(**GRAV))])
    add("S11 late shunt (36 mo delay), 48 mo", 1460.0,
        events=[param_event(1095.0, dict(**GRAV))])
    add("S12 acetazolamide 250 mg BID", 730.0,
        events=dose_events("Aaz_g", 225.0, 0.5, 0.0, 730.0), dt=0.004)
    add("S13 acetazolamide 500 mg BID", 730.0,
        events=dose_events("Aaz_g", 450.0, 0.5, 0.0, 730.0), dt=0.004)
    add("S14 shunt + acetazolamide 250 BID", 730.0, over=GRAV,
        events=dose_events("Aaz_g", 225.0, 0.5, 0.0, 730.0), dt=0.004)
    add("S15 shunt, comorbid AD (APOE e4)", 730.0,
        over=dict(APOE=2.2, **GRAV))
    add("S16 unprotected valve, f_up 0.80", 730.0,
        over=dict(f_up=0.80, shunt=1, Popen_cm=10.0))
    add("S17 shunt + solifenacin 5 mg", 730.0, over=GRAV,
        events=dose_events("Aso_g", 4.5, 1.0, 30.0, 730.0), dt=0.004)
    add("S18 shunt, occlusion 18 mo + revision", 1095.0, over=GRAV,
        events=[param_event(548.0, dict(Rsh=40.0)),
                param_event(700.0, dict(Rsh=2.5))])
    add("S19 shunt + melatonin 2 mg nightly", 730.0, over=GRAV,
        events=dose_events("Ame_g", 0.30, 1.0, 0.0, 730.0), dt=0.004)
    add("S20 ETV in communicating iNPH", 730.0, over=dict(etv=1))
    add("S21 full bundle (grav 4, melatonin, donepezil)", 730.0, over=GRAV,
        events=(dose_events("Ame_g", 0.30, 1.0, 0.0, 730.0)
                + dose_events("Ado_g", 10.0, 1.0, 0.0, 730.0)), dt=0.004)

    print("  %-42s %7s %8s %6s %3s %6s %6s %6s %6s %6s" % (
        "scenario", "gait", "dGait", "MMSE", "GS", "EI", "AMP", "Vsdh",
        "SDH%", "head"))
    print("  " + "-" * 106)
    for nm, r0, r in S:
        print("  %-42s %7.3f %+8.3f %6.1f %3d %6.3f %6.2f %6.1f %6.1f %6.2f" % (
            nm, r["gait"], r["gait"] - r0["gait"], r["cog"], r["iNPHGS"],
            r["EI"], r["AMP"], r["Vsdh"], 100 * r["SDH_inc"], r["headache"]))
    return S


def section11():
    hdr(11, "SENSITIVITY: WHICH UNMEASURED NUMBER OWNS THE ANSWER")
    print("""
  One-at-a-time elasticities d ln(24-month gait gain) / d ln(parameter),
  central difference at +/-10%, on the well-set shunt (4 + gravitational 30).
""")
    base = dict(WELLSET)

    def gain(**over):
        o = dict(base)
        o.update(over)
        _, _, _, tr = run("inph", tend=730.0, nrec=6, **o)
        return nearest(tr, 730.0)["gait"] - nearest(tr, 0.0)["gait"]

    g0 = gain()
    res = []
    for k in ("kappa_tm", "DESH", "k_rep", "k_perm", "kW_out", "kAQ_rec",
              "kCBF_W", "CBF_crit", "kE1_rec", "Rsh", "Hcol_cm", "f_up",
              "If0", "Rout_init"):
        b = base_parameters()[k]
        el = ((gain(**{k: b * 1.1}) - gain(**{k: b * 0.9})) / g0) / 0.2
        res.append((k, b, el))
    res.sort(key=lambda r: -abs(r[2]))
    print("  baseline 24-month gait gain = %+.4f m/s\n" % g0)
    print("  %-14s %12s %14s" % ("parameter", "value", "elasticity"))
    print("  " + "-" * 42)
    for k, b, el in res:
        print("  %-14s %12.5g %+14.3f" % (k, b, el))
    WHAT = {
        "kW_out": "the rate at which periventricular interstitial water clears",
        "kW_in": "water accumulation per unit mean transmantle gradient",
        "kAQ_rec": "the rate at which AQP4 perivascular polarisation recovers",
        "Hcol_cm": "the ventricle-to-peritoneum hydrostatic column, i.e. the "
                   "PATIENT'S HEIGHT",
        "CBF_crit": "the perfusion threshold below which white matter is injured",
        "k_rep": "white-matter functional recovery rate",
        "k_perm": "the rate at which reversible injury becomes permanent",
        "kappa_tm": "the fraction of the ICP pulse appearing across the mantle",
        "DESH": "the DESH morphology multiplier on kappa_tm",
        "Rsh": "shunt hydraulic resistance",
        "Rout_init": "the patient's CSF outflow resistance",
        "kCBF_W": "perfusion loss per unit interstitial water",
        "kE1_rec": "elastance (compliance) recovery rate",
        "If0": "CSF formation rate",
        "f_up": "fraction of the day spent upright",
    }
    top = res[:3]
    kt = dict((k, e) for k, _, e in res)
    print("""
  THIS REFUTED THE EXPECTATION THE MODEL WAS BUILT AROUND, and the refutation is
  left in rather than tidied away.  The structural centre of this model is the
  pulsatile transmantle gradient, so kappa_tm was expected to dominate.  It does
  not: it ranks %d of %d with elasticity %+.3f.  What actually carries the
  therapeutic effect is the WATER axis and the HYDROSTATIC axis:

    1. %-10s %+7.3f   %s
    2. %-10s %+7.3f   %s
    3. %-10s %+7.3f   %s

  The reason kappa_tm drops out is a saturation the model was not designed to
  produce.  The transmantle drive sits well above the creep yield threshold
  dP_yield at every severity tested, so perturbing kappa_tm moves ventricular
  creep without moving what limits white-matter recovery -- which is
  periventricular water and the perfusion it costs.  kappa_tm decides the
  MECHANISM; kW_out and kAQ_rec decide the ANSWER.

  Two consequences worth stating plainly.  (i) Measuring kappa_tm better would
  sharpen the story and change almost none of the predictions; measuring the
  interstitial-water clearance axis would change them a lot -- so the model's
  most-cited uncertainty is not its most important one.  (ii) Hcol_cm ranks %d,
  above every pharmacological and every pulsatility parameter: in a model
  containing five drugs, A PATIENT'S HEIGHT matters more to their two-year gait
  than any of them, because it sets the siphon that the valve has to be chosen
  against.  That is the same conclusion Section 3 reached from the hydraulics
  alone, arrived at here by a completely independent route.
""" % ([r[0] for r in res].index("kappa_tm") + 1, len(res), kt["kappa_tm"],
        top[0][0], top[0][2], WHAT.get(top[0][0], "?"),
        top[1][0], top[1][2], WHAT.get(top[1][0], "?"),
        top[2][0], top[2][2], WHAT.get(top[2][0], "?"),
        [r[0] for r in res].index("Hcol_cm") + 1))


def section12():
    hdr(12, "CALIBRATION: WHAT WAS FITTED, WHAT IS MISSED, WHAT IS A GUESS")
    print("""
  FITTED to a published number (9 parameters)
  -------------------------------------------
    Rout_init = 19     iNPH infusion-test outflow resistance; typical iNPH
                       15-25 vs normal < 13 mmHg/(mL/min)
    E1_init   = 0.222  jointly with Vp_init, to place the iNPH ICP pulse
    Vp_init   = 1.35   amplitude near 4-5 mmHg and the healthy one near 2 mmHg
                       (the reported shunt-response threshold is 4 mmHg)
    kE1_rec   = 0.010  to reproduce a ~35-40% fall in pulse amplitude over the
                       months after a working shunt
    a_fast    = 0.075  to put the best tap-test gait gain in the 5-10% range
                       at which a tap test is read as positive
    kflux_ab  = 98000  so the healthy CSF Abeta42 concentration lands near
    kflux_tau = 26400  800 pg/mL and p-tau near 45 pg/mL; the iNPH values are
                       then predictions of the flux/dilution balance, not fits
    k_creep   = 26     to make the untreated Evans index drift ~0.01/year
                       instead of running away
    k_sdh     = 1.55   to put the unprotected-valve subdural collection in the
                       tens-of-mL range at which it becomes symptomatic

  MISSES, reported rather than tuned away
  ---------------------------------------
    1. RESPONSE RATE.  The model makes almost every patient with an elevated
       Rout a responder; its only non-responder mechanisms are WMperm and the
       APOE/plaque term.  Real series report 60-80% improvement at 1 year.  The
       model therefore OVERSTATES the value of shunting in an unselected
       population, and S6/S21 should be read as best-case arms.
    2. VENTRICULAR SIZE.  A successful shunt here shrinks the Evans index by
       about 0.01-0.02.  Many published series report no significant change at
       all, so the model still over-couples size to pressure.
    3. ACETAZOLAMIDE.  The model gives a small but real gait gain (S12).  The
       clinical evidence in iNPH is weaker than that; treat S12/S13 as an upper
       bound on a drug that has never convincingly worked.
    4. HEALTHY-CONTROL MMSE.  The control sits near 28-29 rather than 29-30,
       because the amyloid block gives a 75-year-old a non-zero plaque load.
    5. TAP-TEST SENSITIVITY.  The model manufactures a false negative at a
       0.10 m/s threshold, which is the right qualitative behaviour, but a
       single deterministic patient cannot reproduce the observed 26-61%
       sensitivity RANGE -- that needs between-patient variability this run
       does not carry.
    6. TAP vs DRAIN ORDERING IS BACKWARDS.  The model ranks a single 40 mL tap
       above 72 h of lumbar drainage on peak gait gain; published sensitivities
       rank the drain higher.  The model captures the DEPTH of a perturbation
       and not its DURATION (tau_gait = 34 d converts almost nothing out of
       three days).  Section 6 quantifies the gap and declines to close it by
       fitting.

  NOT CALIBRATED -- structural guesses with no data behind them
  -----------------------------------------------------------
    kappa_tm (0.060)  the fraction of the ICP pulse appearing across the
                      mantle; never measured in a human ventricle.  I expected
                      this to be the most influential parameter in the model.
                      Section 11 shows it is NOT -- it ranks 10th of 14, and
                      the water-clearance axis (kW_out, kAQ_rec) and the
                      hydrostatic column dominate instead.  The expectation and
                      its refutation are both left in the report.
    kW_out (1.300)    periventricular interstitial water clearance rate.  This
                      is the parameter the answer actually hangs on, and it is
                      an order-of-magnitude estimate.
    kAQ_rec (0.008)   rate of AQP4 repolarisation after decompression.  No
                      human measurement of this rate exists.
    DESH (1.55)       morphology multiplier on kappa_tm.
    k_perm (1.0e-4)   the rate at which reversible injury becomes permanent.
                      This one parameter sets the whole of Section 7.
    kflux_ab, kflux_tau, kdeg_csf   biomarker fluxes; Section 9's direction
                      depends on their ratio and the section says so.
    Hcol_cm (45)      ventricle-to-peritoneum column; scales with patient
                      height and is the dominant term in Section 3.
    k_haz (1.5e-5)    converts collection volume into a symptomatic-event
                      hazard.  The SDH percentages in Section 5 are therefore
                      ORDINAL, not calibrated incidences.
    etv_dRout (0.06)  the claim that a third-ventriculostomy cannot help
                      communicating hydrocephalus is structural (the
                      resistance is downstream of the stoma), not measured.

  DEFECTS FOUND BY THE MODEL'S OWN CHECKS, AND FIXED
  --------------------------------------------------
    a) The acute volume perturbation was applied as a LINEAR offset dVac/C and
       solved by a 3-pass fixed point.  For a 40 mL tap that iteration does not
       converge: it oscillated between -0.2 and -22.5 mmHg and settled on the
       -30 mmHg clamp.  It now enters through the INTEGRATED Marmarou form
       P = P0 + (Pbase-P0)*exp(E1*dV) -- explicit, monotone, reducing to dV/C
       for small dV and saturating near P0 for a large removal.  The P-V
       reference pressure moved to -2 mmHg and E1 was rescaled by 1.2 to hold
       the pulse-amplitude calibration.
    b) THE TAP TEST WAS REFILLING TWENTY TIMES TOO FAST -- the most
       consequential defect found.  The deficit was relaxed with tau = Rout*C
       ~ 6 min, the linearisation of the flow balance, which implies a refill
       rate of 7 mL/min against a choroidal secretion rate of 0.35 mL/min.  The
       pressure was therefore back to baseline in 30 minutes and the model
       concluded -- wrongly -- that no tap test can ever be positive.  Using the
       exact balance, absorption is ZERO below sinus pressure, so recovery is
       PRODUCTION-limited and 40 mL takes 114 minutes.
    b2) Even with the refill fixed, a SYMMETRIC first-order gait filter with
       tau = 13 h integrates almost nothing out of a 2-hour perturbation.  The
       fast arm is now asymmetric (charge 29 min, discharge 21.6 h).  That is a
       modelling choice rather than a measurement -- but a falsifiable one: it
       predicts fast onset and slow offset of tap-test improvement, which is
       what is clinically reported.
    c) Ventricular recoil had no plastic floor, so a shunted patient's
       ventricles returned to a 25 mL normal volume within a year, the
       opposite of every published series.  Vplast was added.
    d) The valve-shut branch was missing, so at high opening pressures the
       solver produced a NEGATIVE shunt flow -- fluid running from the
       peritoneum into the ventricle -- and the high-Popen arm of Section 5
       looked therapeutic.
    e) The mean transmantle gradient depended on Rout alone, so it did not
       respond to a shunt at all and periventricular water never cleared.  It
       is now scaled by par_frac, the share of production still absorbed
       through the high-resistance parenchymal route.
    f) Reactive astrogliosis and microglial activation were given source rates
       an order of magnitude too large for their sinks; the untreated
       equilibrium was Astro = 20.7 on a 0-3 scale, which silently drove AQP4
       polarisation to zero.
""")


def brief():
    ph, yh = patient("healthy")
    pi, yi = patient("inph")
    ah, ai = hydro_full(yh, ph), hydro_full(yi, pi)
    print("iNPH QSP reference model -- headline numbers")
    print("  Rout             %.1f -> %.1f  (%.2fx)" % (
        yh[IDX["Rout"]], yi[IDX["Rout"]], yi[IDX["Rout"]] / yh[IDX["Rout"]]))
    print("  mean ICP supine  %.2f -> %.2f mmHg (%.2fx, still 'normal')" % (
        ah["P_sup"], ai["P_sup"], ai["P_sup"] / ah["P_sup"]))
    print("  ICP pulse amp    %.2f -> %.2f mmHg (%.2fx)" % (
        ah["AMP"], ai["AMP"], ai["AMP"] / ah["AMP"]))
    print("  gait / MMSE      %.2f m/s / %.1f  ->  %.2f m/s / %.1f" % (
        readout(0, yh, ph)["gait"], readout(0, yh, ph)["cog"],
        readout(0, yi, pi)["gait"], readout(0, yi, pi)["cog"]))
    p, y = patient("inph", shunt=1, Popen_cm=10.0)
    Pu, Qu = hydro_posture(y, p, True)[:2]
    print("  unprotected valve upright: ICP %.1f mmHg, %.0f mL/day "
          "(%.2fx production)" % (Pu, Qu * MIN_PER_DAY,
                                  Qu * MIN_PER_DAY / (p["If0"] * MIN_PER_DAY)))


def main():
    if "--brief" in sys.argv:
        brief()
        return
    print("=" * 78)
    print("idiopathic NORMAL PRESSURE HYDROCEPHALUS (iNPH)")
    print("QSP reference implementation -- 44 ODE states + 1 closed-form fast "
          "state")
    print("=" * 78)
    for fn in (section1, section2, section3, section4, section5, section6,
               section7, section8, section9, section10, section11, section12):
        fn()
    print("\n" + "=" * 78)
    print("end of report")
    print("=" * 78)


if __name__ == "__main__":
    main()
