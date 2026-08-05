#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
op_reference_model.py
=====================
Independent Python/scipy re-implementation of the acute organophosphorus (OP)
insecticide self-poisoning QSP model that is distributed here as
`op_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
The mrgsolve model is the deliverable; this file is the *check*.  Every ODE,
every parameter and every derived output below was written from the same
specification and then compared against the R model.  Where the two disagreed,
the disagreement was treated as a defect in the model, not as a rounding
difference (see `op_reference_output.txt` and the DEFECTS section of
README.md).

THE ORGANISING IDEA
-------------------
Everything the antidotes can and cannot do collapses onto ONE dimensionless
group.  Acetylcholinesterase (AChE) exists in three states

        E  --(k_i * C_oxon)-->  EP  --(k_a)-->  EP_aged        (irreversible)
        E  <--(k_s + k_r[X])--  EP

where C_oxon is the concentration of the *active* oxon metabolite and k_r[X] is
the oxime reactivation rate at oxime concentration X.  Oxime reactivation is
NOT linear in X -- it saturates (Worek's k_r2 = k_r_max / (K_D + X) formalism),
so as X -> infinity the reactivation rate tends to a hard ceiling k_r_max.
Therefore the *best achievable* fraction of active enzyme, at unlimited oxime,
is

        E_ceiling = k_r_max / (k_r_max + k_i * C_oxon)  =  Omega / (1 + Omega)

with the OXIME SUFFICIENCY NUMBER

        Omega = k_r_max / (k_i * C_oxon).

Omega < 1 means the oxime is arithmetically incapable of holding half the
enzyme free, at ANY dose.  Everything else in the model -- the critical
ingested volume, the class split between diethyl and dimethyl OPs, the failure
of the randomised oxime trials, the atropine requirement, the ventilator
requirement -- follows from where Omega sits, and from the second number, the
AGING RATCHET

        phi = k_a / (k_a + k_s + k_r[X])

the probability that any single inhibition event ends irreversibly.

Author: QSP Disease Model Library (Claude Code Routine)
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field, asdict

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
LN2 = np.log(2.0)

# ==========================================================================
# 1. OP CLASS LIBRARY
# ==========================================================================
# Rate constants for HUMAN erythrocyte AChE.  Second-order inhibition constants
# k_i are quoted in the literature as M^-1 min^-1; here they are converted to
# nM^-1 h^-1  (multiply by 1e-9 * 60 = 6e-8).  Oxime reactivation constants
# k_r2 are quoted the same way and converted to uM^-1 h^-1 (x 6e-5), then
# split into a maximum rate k_r_max (h^-1) and a dissociation constant K_D (uM)
# with k_r2 = k_r_max / K_D, which is Worek's formalism and which is what makes
# the ceiling exist.
#
# Sources (see op_references.md): Worek 2004 Biochem Pharmacol 68:2237;
# Worek 2007 Toxicol Appl Pharmacol; Aurbek 2006; Eyer 2003 Toxicol Rev.

@dataclass
class OPClass:
    name: str
    subclass: str                 # 'diethyl' | 'dimethyl'
    mw: float                     # g/mol of the parent thion
    ki: float                     # nM^-1 h^-1, oxon vs human AChE
    ki_nmj_scale: float           # NMJ AChE relative susceptibility
    ki_nte: float                 # nM^-1 h^-1, oxon vs neuropathy target esterase
    t12_aging_h: float            # AChE aging half-time
    t12_spont_h: float            # AChE spontaneous reactivation half-time
    t12_aging_nte_h: float
    # pralidoxime (2-PAM) and obidoxime against THIS oxon
    pam_krmax: float              # h^-1
    pam_KD: float                 # uM
    obi_krmax: float
    obi_KD: float
    # toxicokinetics
    ka_gut: float                 # h^-1 absorption of the thion
    F: float
    V1: float                     # L, central volume of the thion
    Vf: float                     # L, fat/deep depot of the thion
    Qf: float                     # L/h, intercompartmental clearance thion
    Vmax_bio: float               # nmol/h, CYP desulfuration Vmax (thion -> oxon)
    Km_bio: float                 # uM thion
    CL_thion_other: float         # L/h, non-activating clearance of the thion
    CL_oxon: float                # L/h, oxon hydrolysis (PON1 + other) at PON1 RR
    pon1_QQ_scale: float          # multiplier on CL_oxon for the PON1 Q192Q genotype
    V_oxon: float                 # L, oxon distribution volume
    Q_oxon: float                 # L/h, oxon to peripheral tissue
    Vt_oxon: float                # L
    solvent_frac: float           # v/v hydrocarbon / ketone solvent in formulation
    solvent_cv_tox: float         # relative myocardial-depressant potency of the solvent
    label: str = ""


CHLORPYRIFOS = OPClass(
    name="chlorpyrifos", subclass="diethyl", mw=350.6,
    ki=0.18,                # 3.0e6 M^-1 min^-1 (chlorpyrifos-oxon, human AChE)
    ki_nmj_scale=1.0,
    ki_nte=0.0035,          # NTE is inhibited ~50x more slowly than AChE
    t12_aging_h=33.0,       # diethylphosphoryl-AChE
    t12_spont_h=77.0,
    t12_aging_nte_h=6.0,    # aged NTE is what causes OPIDN
    pam_krmax=36.0, pam_KD=200.0,     # k_r2 = 0.18 uM^-1 h^-1 ~ 3000 M^-1 min^-1
    obi_krmax=48.0, obi_KD=40.0,      # k_r2 = 1.20 uM^-1 h^-1 ~ 20000 M^-1 min^-1
    ka_gut=0.35, F=0.80,
    V1=35.0, Vf=700.0, Qf=25.0,       # highly lipophilic: log P 4.96
    Vmax_bio=6.0e4, Km_bio=20.0,
    CL_thion_other=14.0,
    CL_oxon=300.0, pon1_QQ_scale=0.33,
    V_oxon=56.0, Q_oxon=60.0, Vt_oxon=120.0,
    solvent_frac=0.60, solvent_cv_tox=0.5,
    label="Chlorpyrifos 20% EC (diethyl, lipophilic, thion pro-toxicant)",
)

DIMETHOATE = OPClass(
    name="dimethoate", subclass="dimethyl", mw=229.3,
    ki=0.0084,              # omethoate, 1.4e5 M^-1 min^-1: a WEAK inhibitor
    ki_nmj_scale=1.0,
    ki_nte=0.00002,         # dimethyl OPs essentially do not cause OPIDN
    t12_aging_h=3.7,        # dimethylphosphoryl-AChE ages fast
    t12_spont_h=0.7,        # ...but also reactivates spontaneously fast
    t12_aging_nte_h=6.0,
    pam_krmax=15.0, pam_KD=300.0,     # k_r2 = 0.05 uM^-1 h^-1 ~ 830 M^-1 min^-1
    obi_krmax=30.0, obi_KD=100.0,     # k_r2 = 0.30 uM^-1 h^-1 ~ 5000 M^-1 min^-1
    ka_gut=0.90, F=0.90,
    V1=42.0, Vf=140.0, Qf=40.0,       # log P 0.78: little fat depot
    Vmax_bio=1.5e5, Km_bio=20.0,
    CL_thion_other=45.0,
    CL_oxon=30.0, pon1_QQ_scale=0.85,  # PON1 barely hydrolyses dimethyl oxons
    V_oxon=50.0, Q_oxon=50.0, Vt_oxon=90.0,
    solvent_frac=0.55, solvent_cv_tox=1.8,   # cyclohexanone: cardiotoxic
    label="Dimethoate 40% EC (dimethyl, cyclohexanone vehicle)",
)

FENTHION = OPClass(
    name="fenthion", subclass="dimethyl", mw=278.3,
    ki=0.050,
    ki_nmj_scale=1.6,       # fenthion is notorious for nicotinic/IMS predominance
    ki_nte=0.00004,
    t12_aging_h=3.7,
    t12_spont_h=0.9,
    t12_aging_nte_h=6.0,
    pam_krmax=15.0, pam_KD=300.0,
    obi_krmax=30.0, obi_KD=100.0,
    ka_gut=0.10, F=0.85,               # slow, erratic absorption
    V1=30.0, Vf=1750.0, Qf=6.0,        # log P 4.84 and a very slow depot
    Vmax_bio=4.0e4, Km_bio=20.0,
    CL_thion_other=12.0,
    CL_oxon=90.0, pon1_QQ_scale=0.80,
    V_oxon=56.0, Q_oxon=50.0, Vt_oxon=120.0,
    solvent_frac=0.50, solvent_cv_tox=0.4,
    label="Fenthion 50% EC (dimethyl, extremely lipophilic, delayed onset)",
)

PARATHION = OPClass(
    name="parathion", subclass="diethyl", mw=291.3,
    ki=0.084,               # paraoxon-ethyl, 1.4e6 M^-1 min^-1
    ki_nmj_scale=1.0,
    ki_nte=0.0016,
    t12_aging_h=33.0,
    t12_spont_h=77.0,
    t12_aging_nte_h=6.0,
    pam_krmax=36.0, pam_KD=200.0,
    obi_krmax=48.0, obi_KD=40.0,
    ka_gut=0.45, F=0.85,
    V1=40.0, Vf=500.0, Qf=25.0,
    Vmax_bio=8.0e4, Km_bio=20.0,
    CL_thion_other=12.0,
    CL_oxon=420.0, pon1_QQ_scale=0.28,   # paraoxon is the classic PON1 substrate
    V_oxon=56.0, Q_oxon=60.0, Vt_oxon=120.0,
    solvent_frac=0.55, solvent_cv_tox=0.4,
    label="Parathion 50% EC (diethyl; the compound the oximes were built for)",
)

OP_LIBRARY = {c.name: c for c in (CHLORPYRIFOS, DIMETHOATE, FENTHION, PARATHION)}

# ==========================================================================
# 2. HOST / PHARMACODYNAMIC PARAMETERS
# ==========================================================================

P = dict(
    WT=60.0,                      # kg -- the median South Asian self-poisoning patient

    # ---- cholinesterase pool turnover (h^-1) --------------------------------
    krbc_new=LN2 / (60.0 * 24.0),     # RBC AChE: replaced only by erythropoiesis
    ktiss=LN2 / (5.0 * 24.0),         # synaptic/muscle/brain AChE protein turnover
    kbche=LN2 / (11.0 * 24.0),        # plasma BChE, hepatic synthesis
    knte=LN2 / (5.0 * 24.0),

    # ---- acetylcholine handling --------------------------------------------
    khyd=60.0,                    # h^-1, cleft ACh clearance at full AChE
    ach_leak=0.05,                # non-AChE hydrolysis + diffusion floor
    # -> ACh_ss = (1+leak)/(E+leak): 1.0 at E=1, 21x normal at E=0

    # ---- receptor occupancy -------------------------------------------------
    KA_m=3.0,                     # muscarinic: baseline occupancy 1/(1+3) = 0.25
    KA_n=5.0,                     # nicotinic (lower affinity)
    KA_b=3.0,                     # central muscarinic
    Om0=0.25, Ommax=0.87,
    On0=1.0 / 6.0,
    KB_atr=8.0,                   # nM, effective muscarinic Kb of atropine

    # ---- nicotinic desensitisation -> intermediate syndrome -----------------
    kdes=0.075,                   # h^-1, slow: this is why IMS appears at 24-96 h
    krec_des=0.050,               # h^-1
    AChn50=8.0, hill_n=4.0,

    # ---- cholinergic receptor down-regulation (the reason a patient can wean
    #      while the enzyme is still <10% -- mAChR/nAChR down-regulation under
    #      sustained agonist load; t1/2 up ~1 d, t1/2 off ~4 d)
    kdown=0.012, krec_down=0.005,
    fdown_m=0.65,          # attenuation of PERIPHERAL muscarinic signalling
    fdown_c=0.60,          # attenuation of CENTRAL muscarinic signalling
    fdown_n=0.60,          # attenuation of nicotinic desensitisation drive

    # ---- secretions / airway ------------------------------------------------
    ksec=6.0,                     # mL/h per unit muscarinic excess
    kclear_sec=0.9,               # spontaneous clearance (coughing/swallowing)
    kclear_icu=2.4,               # + suctioning
    SEC_crit=2.0,                 # mL burden at which the airway is soaked
    kBT=1.2, kBT_off=1.5,         # bronchial tone on/off

    # ---- central respiratory drive -----------------------------------------
    tau_respd=0.4,                # h
    Dmax_cns=0.55,                # central depression that ATROPINE can reverse
    Dmax_nonrev=0.65,             # central depression it CANNOT (non-muscarinic
                                  # brainstem cholinergic + nicotinic drive)
    fdown_resp=0.80,              # adaptation of the non-reversible component
    solv_resp=0.55,               # solvent CNS depression weight

    # ---- neuromuscular ------------------------------------------------------
    tau_mstr=0.3,
    w_fast_block=0.45,            # immediate depolarising block at very high ACh

    # ---- haemodynamics ------------------------------------------------------
    tau_hr=0.05, HR_int=155.0, HR_slope=320.0,
    tau_map=0.2, MAP0=92.0,
    w_map_mus=50.0, w_map_solv=30.0, w_map_hyp=25.0,

    # ---- gas exchange -------------------------------------------------------
    SF_nmj=0.28,                  # resting ventilation needs only ~28% of the
                                  # maximal inspiratory capacity: the NMJ has a
                                  # safety factor that central drive does not
    VA_vent_thresh=0.60,
    shunt0=0.03, w_shunt_sec=0.40, w_shunt_lung=0.35, shunt_max=0.60,
    FiO2_air=0.21, FiO2_O2=0.60, FiO2_vent=0.50,

    # ---- seizures -----------------------------------------------------------
    kseiz=1.4, kseiz_off=0.6, seiz_thresh=0.62,   # central occupancy threshold

    # ---- aspiration / solvent lung injury -----------------------------------
    k_lung=0.35, k_lung_rep=0.030,

    # ---- hazard -------------------------------------------------------------
    h0=0.00015,
    h_hypox=0.55, h_hyperc=0.060, h_shock=0.45, h_seiz=0.05, h_lung=0.004,
    h_vent=0.0012,                # ventilator-associated complications, per h

    # ---- OPIDN --------------------------------------------------------------
    kopidn=0.006, krec_opidn=0.0015, nte_thresh=0.70,

    # ---- antidote PK ---------------------------------------------------------
    # pralidoxime chloride MW 172.6; obidoxime dichloride MW 359.2
    pam_MW=172.6, obi_MW=359.2,
    pam_V1=24.0, pam_V2=24.0, pam_CL=21.0, pam_Q=30.0,      # L, L/h (60 kg)
    obi_V1=15.0, obi_V2=12.0, obi_CL=7.5, obi_Q=9.0,
    atr_MW=289.4, atr_V=200.0, atr_CL=46.0, atr_keo=6.0,    # L, L/h, h^-1
    dz_MW=284.7, dz_V=80.0, dz_CL=1.6, dz_EC50=0.20,        # uM
    mg_V=15.0, mg_CL=6.0, mg_EC50=1.6, mg_Emax=0.45,        # mmol, L/h, mmol/L

    # ---- atropine controller -------------------------------------------------
    Kp_rapid=170.0,               # umol/h per unit atropinisation error
    Kp_standard=22.0,             # the "2 mg every 10-15 min, ad hoc" habit
    atr_max_rate=210.0,           # umol/h ~ 60 mg/h
    hr_tachy_stop=140.0,

    # ---- decontamination -----------------------------------------------------
    k_charcoal=1.8,               # h^-1 extra gut removal while charcoal present
    charcoal_dur=2.0,             # h of effect from a single 50 g dose
)

# ==========================================================================
# 3. STATE VECTOR
# ==========================================================================

STATES = [
    "A_gut",      # 0  umol thion in gut lumen
    "A_th_c",     # 1  umol thion, central
    "A_th_f",     # 2  umol thion, fat depot
    "A_ox_c",     # 3  nmol oxon, central
    "A_ox_t",     # 4  nmol oxon, peripheral tissue
    "A_solv_g",   # 5  mL solvent in gut
    "A_solv_c",   # 6  mL-equivalent solvent, systemic
    "LUNG",       # 7  aspiration / solvent lung injury (0-1)
    "A_pam_c",    # 8  umol oxime central
    "A_pam_p",    # 9  umol oxime peripheral
    "A_atr_c",    # 10 umol atropine central
    "Ce_atr",     # 11 nM atropine effect site
    "A_dz",       # 12 umol diazepam
    "A_mg",       # 13 mmol magnesium
    "SCAV",       # 14 umol active bioscavenger sites
    "E_rbc", "EP_rbc", "EA_rbc",       # 15-17
    "E_mus", "EP_mus", "EA_mus",       # 18-20  peripheral muscarinic target AChE
    "E_nmj", "EP_nmj", "EA_nmj",       # 21-23  neuromuscular junction AChE
    "E_cns", "EP_cns", "EA_cns",       # 24-26  brain AChE
    "B_free", "B_inh",                 # 27-28  plasma butyrylcholinesterase
    "N_free", "N_inh", "N_aged",       # 29-31  neuropathy target esterase
    "ACh_m", "ACh_n", "ACh_b",         # 32-34  cleft ACh, x normal
    "Rn_des",     # 35 desensitised / depolarisation-blocked nAChR fraction
    "Rm_down",    # 36 muscarinic receptor down-regulation
    "SEC",        # 37 airway secretion burden (mL)
    "BT",         # 38 bronchial tone (0-1)
    "HR",         # 39 bpm
    "RESPD",      # 40 central respiratory drive (fraction of normal)
    "MSTR",       # 41 respiratory muscle strength (fraction)
    "MAP",        # 42 mmHg
    "SEIZ",       # 43 seizure / CNS excitation index
    "OPIDN",      # 44 delayed polyneuropathy score
    "HAZ",        # 45 cumulative death hazard
    "VTIME",      # 46 h of mechanical ventilation
    "ATRCUM",     # 47 mg atropine given
    "OXCUM",      # 48 mg oxime given
    "AUC_EPn",    # 49 h, integral of the inhibited NMJ enzyme fraction
    "AUC_ox",     # 50 nM.h oxon exposure
]
IDX = {s: i for i, s in enumerate(STATES)}
NST = len(STATES)


def initial_state() -> np.ndarray:
    y = np.zeros(NST)
    for k in ("E_rbc", "E_mus", "E_nmj", "E_cns", "B_free", "N_free"):
        y[IDX[k]] = 1.0
    for k in ("ACh_m", "ACh_n", "ACh_b"):
        y[IDX[k]] = 1.0
    y[IDX["HR"]] = 75.0
    y[IDX["RESPD"]] = 1.0
    y[IDX["MSTR"]] = 1.0
    y[IDX["MAP"]] = P["MAP0"]
    return y


# ==========================================================================
# 4. TREATMENT SPECIFICATION
# ==========================================================================

@dataclass
class Treatment:
    label: str = "supportive"
    # ---- exposure
    op: str = "chlorpyrifos"
    volume_ml: float = 50.0
    conc_g_per_L: float = 200.0
    pon1: str = "RR"                 # 'RR' fast paraoxonase | 'QQ' slow
    # ---- decontamination
    charcoal_h: float | None = None
    # ---- oxime
    oxime: str | None = None         # 'pam' | 'obi'
    oxime_start_h: float = 1.5
    oxime_load_mg_kg: float = 30.0
    oxime_inf_mg_kg_h: float = 8.0
    oxime_dur_h: float = 168.0
    oxime_bolus_mg: float | None = None      # if set, use intermittent boluses
    oxime_bolus_int_h: float = 6.0
    # ---- antimuscarinic
    atropine: str = "rapid"          # 'rapid' | 'standard' | 'glyco' | 'none'
    # ---- adjuncts
    diazepam: bool = True
    magnesium: bool = False
    scavenger_mg: float = 0.0        # plasma-derived human BChE, mg
    # ---- supportive care
    ventilator: bool = True
    oxygen: bool = True
    icu_suction: bool = True


def bbb_factor(t: Treatment) -> float:
    """Fraction of the antimuscarinic reaching central muscarinic receptors."""
    return 0.02 if t.atropine == "glyco" else 0.35


# ==========================================================================
# 5. RIGHT-HAND SIDE
# ==========================================================================

def _sig(x, k=25.0):
    return 1.0 / (1.0 + np.exp(-k * x))


def derivatives(t_h, y, op: OPClass, tx: Treatment, aux=None):
    p = P
    y = np.maximum(y, 0.0)
    s = {k: y[i] for k, i in IDX.items()}

    # ---------------------------------------------------------------- oxime
    if tx.oxime == "obi":
        krmax, KD, oxV1, oxV2, oxCL, oxQ = (op.obi_krmax, op.obi_KD,
                                            p["obi_V1"], p["obi_V2"],
                                            p["obi_CL"], p["obi_Q"])
    else:
        krmax, KD, oxV1, oxV2, oxCL, oxQ = (op.pam_krmax, op.pam_KD,
                                            p["pam_V1"], p["pam_V2"],
                                            p["pam_CL"], p["pam_Q"])
    C_ox_drug = s["A_pam_c"] / oxV1                      # uM
    kr = krmax * C_ox_drug / (KD + C_ox_drug) if tx.oxime else 0.0

    # ------------------------------------------------- toxicokinetics: thion
    ka = op.ka_gut
    if tx.charcoal_h is not None and tx.charcoal_h <= t_h < tx.charcoal_h + p["charcoal_dur"]:
        k_ac = p["k_charcoal"]
    else:
        k_ac = 0.0
    C_th = s["A_th_c"] / op.V1                            # uM
    C_thf = s["A_th_f"] / op.Vf
    bio = op.Vmax_bio * C_th / (op.Km_bio + C_th)         # nmol/h oxon formed

    dA_gut = -(ka + k_ac) * s["A_gut"]
    dA_th_c = (op.F * ka * s["A_gut"]
               - bio / 1000.0                              # nmol/h -> umol/h
               - op.CL_thion_other * C_th
               - op.Qf * (C_th - C_thf))
    dA_th_f = op.Qf * (C_th - C_thf)

    # ------------------------------------------------- toxicokinetics: oxon
    CLox = op.CL_oxon * (op.pon1_QQ_scale if tx.pon1 == "QQ" else 1.0)
    C_ox = s["A_ox_c"] / op.V_oxon                        # nM
    C_oxt = s["A_ox_t"] / op.Vt_oxon
    # stoichiometric capture by an exogenous bioscavenger (1:1 molar)
    scav_rate = 3.0e-3 * s["SCAV"] * C_ox                 # nmol/h
    dA_ox_c = (bio - CLox * C_ox - op.Q_oxon * (C_ox - C_oxt) - scav_rate)
    dA_ox_t = op.Q_oxon * (C_ox - C_oxt) - 12.0 * C_oxt
    dSCAV = -scav_rate / 1000.0

    # ------------------------------------------------------------- solvent
    dA_solv_g = -1.1 * s["A_solv_g"]
    dA_solv_c = 0.55 * 1.1 * s["A_solv_g"] - 0.55 * s["A_solv_c"]
    solv = s["A_solv_c"] / (s["A_solv_c"] + 12.0)         # 0-1 saturating

    # ------------------------------------------------------- antidote PK
    dA_pam_c = (-oxCL * C_ox_drug - oxQ * (C_ox_drug - s["A_pam_p"] / oxV2))
    dA_pam_p = oxQ * (C_ox_drug - s["A_pam_p"] / oxV2)

    # atropine: a CLOSED-LOOP drug.  The infusion rate is a function of the
    # patient's state, not of the clock, which is exactly how it is titrated at
    # the bedside -- and it makes the cumulative dose an OUTPUT of the model.
    Om_pre = 0.0
    Cat = s["Ce_atr"]
    A_m = s["ACh_m"]
    Om = A_m / (A_m + p["KA_m"] * (1.0 + Cat / p["KB_atr"]))
    mus_raw = np.clip((Om - p["Om0"]) / (p["Ommax"] - p["Om0"]), 0.0, 1.0)
    # receptor down-regulation attenuates the SIGNAL, not the occupancy
    mus_x = mus_raw * (1.0 - p["fdown_m"] * s["Rm_down"])

    if tx.atropine in ("rapid", "standard", "glyco"):
        Kp = p["Kp_rapid"] if tx.atropine in ("rapid", "glyco") else p["Kp_standard"]
        # bedside atropinisation targets: dry chest, HR > 80, no bronchospasm,
        # and adequate spontaneous ventilation
        err = (0.40 * np.clip(s["SEC"] / p["SEC_crit"] - 0.12, 0.0, 1.0)
               + 0.25 * np.clip((78.0 - s["HR"]) / 40.0, 0.0, 1.0)
               + 0.15 * np.clip(s["BT"] - 0.05, 0.0, 1.0)
               + 0.20 * np.clip((0.75 - s["RESPD"]) / 0.5, 0.0, 1.0))
        stop = _sig((p["hr_tachy_stop"] - s["HR"]) / 10.0, k=6.0)
        rate_atr = min(Kp * err * stop, p["atr_max_rate"])
    else:
        rate_atr = 0.0
    dA_atr_c = rate_atr - p["atr_CL"] * s["A_atr_c"] / p["atr_V"]
    C_atr_pl = 1000.0 * s["A_atr_c"] / p["atr_V"]          # nM
    dCe_atr = p["atr_keo"] * (C_atr_pl - s["Ce_atr"])

    dA_dz = -p["dz_CL"] * s["A_dz"] / p["dz_V"]
    C_dz = s["A_dz"] / p["dz_V"]
    dA_mg = -p["mg_CL"] * s["A_mg"] / p["mg_V"]
    C_mg = s["A_mg"] / p["mg_V"]

    # ------------------------------------------------------ esterase kinetics
    ka_age = LN2 / op.t12_aging_h
    ks = LN2 / op.t12_spont_h
    kinh = op.ki * C_ox                                    # h^-1 (central pools)
    kinh_t = op.ki * C_oxt                                 # tissue-side pools

    def esterase(E, EP, EA, kin, kturn, k_react, k_age):
        tot = E + EP + EA
        dE = -kin * E + k_react * EP + kturn * (1.0 - tot)
        dEP = kin * E - (k_react + k_age + kturn) * EP
        dEA = k_age * EP - kturn * EA
        return dE, dEP, dEA

    # RBC AChE: no protein turnover, only erythropoietic replacement
    dE_rbc, dEP_rbc, dEA_rbc = esterase(
        s["E_rbc"], s["EP_rbc"], s["EA_rbc"], kinh, p["krbc_new"], ks + kr, ka_age)
    # peripheral muscarinic target tissue
    dE_mus, dEP_mus, dEA_mus = esterase(
        s["E_mus"], s["EP_mus"], s["EA_mus"], kinh_t, p["ktiss"], ks + kr, ka_age)
    # neuromuscular junction -- oximes reach it well (peripheral, no BBB)
    dE_nmj, dEP_nmj, dEA_nmj = esterase(
        s["E_nmj"], s["EP_nmj"], s["EA_nmj"], kinh_t * op.ki_nmj_scale,
        p["ktiss"], ks + kr, ka_age)
    # brain -- the quaternary oximes barely cross the blood-brain barrier
    kr_cns = 0.05 * kr
    dE_cns, dEP_cns, dEA_cns = esterase(
        s["E_cns"], s["EP_cns"], s["EA_cns"], kinh_t * 0.8,
        p["ktiss"], ks + kr_cns, ka_age)
    # plasma BChE (biomarker; aging lumped into the inhibited pool)
    dB_free = -kinh * 2.5 * s["B_free"] + p["kbche"] * (1.0 - s["B_free"] - s["B_inh"])
    dB_inh = kinh * 2.5 * s["B_free"] - p["kbche"] * s["B_inh"]
    # neuropathy target esterase
    ka_nte = LN2 / op.t12_aging_nte_h
    dN_free, dN_inh, dN_aged = esterase(
        s["N_free"], s["N_inh"], s["N_aged"], op.ki_nte * C_oxt,
        p["knte"], 0.0, ka_nte)

    # --------------------------------------------------------- acetylcholine
    mg_block = 1.0 - p["mg_Emax"] * C_mg / (p["mg_EC50"] + C_mg)
    dz_block = 1.0 - 0.25 * C_dz / (p["dz_EC50"] + C_dz)     # presynaptic damping
    rel = mg_block
    rel_b = mg_block * dz_block

    def ach_ode(A, E):
        return p["khyd"] * (1.0 + p["ach_leak"]) * rel - p["khyd"] * (E + p["ach_leak"]) * A

    dACh_m = ach_ode(s["ACh_m"], s["E_mus"])
    dACh_n = ach_ode(s["ACh_n"], s["E_nmj"])
    dACh_b = (p["khyd"] * (1.0 + p["ach_leak"]) * rel_b
              - p["khyd"] * (s["E_cns"] + p["ach_leak"]) * s["ACh_b"])

    # ------------------------------------------------------------- receptors
    A_n = s["ACh_n"]  # noqa: E501
    On = A_n / (A_n + p["KA_n"])
    nic_x = np.clip((On - p["On0"]) / (1.0 - p["On0"]), 0.0, 1.0)

    Cat_b = s["Ce_atr"] * bbb_factor(tx)
    A_b = s["ACh_b"]
    Ob = A_b / (A_b + p["KA_b"] * (1.0 + Cat_b / p["KB_atr"]))
    cns_raw = np.clip((Ob - p["Om0"]) / (p["Ommax"] - p["Om0"]), 0.0, 1.0)
    cns_x = cns_raw * (1.0 - p["fdown_c"] * s["Rm_down"])

    fdes_raw = A_n ** p["hill_n"] / (A_n ** p["hill_n"] + p["AChn50"] ** p["hill_n"])
    fdes = fdes_raw * (1.0 - p["fdown_n"] * s["Rm_down"])
    dRn_des = p["kdes"] * fdes * (1.0 - s["Rn_des"]) - p["krec_des"] * s["Rn_des"]
    dRm_down = (p["kdown"] * max(mus_raw, cns_raw, nic_x) * (1.0 - s["Rm_down"])
                - p["krec_down"] * s["Rm_down"])

    # ------------------------------------------------------------ end organs
    kclr = p["kclear_icu"] if tx.icu_suction else p["kclear_sec"]
    dSEC = p["ksec"] * mus_x - kclr * s["SEC"]
    dBT = p["kBT"] * mus_x * (1.0 - s["BT"]) - p["kBT_off"] * s["BT"]

    Om_eff = p["Om0"] + (Om - p["Om0"]) * (1.0 - p["fdown_m"] * s["Rm_down"])
    HR_t = np.clip(p["HR_int"] - p["HR_slope"] * Om_eff + 25.0 * nic_x, 22.0, 170.0)
    dHR = (HR_t - s["HR"]) / p["tau_hr"]

    # Central respiratory depression has an atropine-REVERSIBLE muscarinic part
    # and a part that atropine cannot touch; only receptor adaptation removes the
    # second one, which is why the patient weans on a clock of days rather than
    # on the atropine chart.
    Lb = np.clip((s["ACh_b"] - 1.0) / 20.0, 0.0, 1.0)
    Lb_ad = Lb * (1.0 - p["fdown_resp"] * s["Rm_down"])
    RESPD_t = np.clip(1.0 - p["Dmax_cns"] * cns_x - p["Dmax_nonrev"] * Lb_ad
                      - p["solv_resp"] * solv, 0.02, 1.0)
    dRESPD = (RESPD_t - s["RESPD"]) / p["tau_respd"]

    MSTR_t = np.clip((1.0 - s["Rn_des"]) * (1.0 - p["w_fast_block"] * fdes), 0.02, 1.0)
    dMSTR = (MSTR_t - s["MSTR"]) / p["tau_mstr"]

    # ------------------------------------------------------- gas exchange
    # Ventilation is the SMALLER of what the brain asks for and what the
    # respiratory muscles can deliver.  Resting ventilation needs only ~28% of
    # maximal inspiratory capacity, so neuromuscular block has a safety factor
    # that central depression does not.
    Vcap = min(s["RESPD"], s["MSTR"] / p["SF_nmj"]) * (1.0 - 0.40 * s["BT"])
    vent_on = _sig((p["VA_vent_thresh"] - Vcap) / 0.05) if tx.ventilator else 0.0
    VA = Vcap + vent_on * (1.05 - Vcap)
    VA = min(max(VA, 0.06), 1.05)
    PaCO2 = 40.0 / VA
    if tx.ventilator and vent_on > 0.5:
        FiO2 = p["FiO2_vent"]
    elif tx.oxygen and Vcap < 0.75:
        FiO2 = p["FiO2_O2"]
    else:
        FiO2 = p["FiO2_air"]
    shunt = min(p["shunt_max"],
                p["shunt0"] + p["w_shunt_sec"] * min(1.0, s["SEC"] / 3.0)
                + p["w_shunt_lung"] * s["LUNG"])
    PAO2 = FiO2 * 713.0 - PaCO2 / 0.8
    PaO2 = max(8.0, (1.0 - shunt) * PAO2 + shunt * 40.0)
    hyp = np.clip((60.0 - PaO2) / 60.0, 0.0, 1.0)

    # Aspiration is a THRESHOLD event: it needs an unprotected airway that is
    # actually flooded, or a hydrocarbon vehicle going down the wrong way.
    asp_drive = (0.60 * solv
                 + 0.90 * np.clip(s["SEC"] / p["SEC_crit"] - 0.60, 0.0, 1.0)
                 * (1.0 - vent_on))
    dLUNG = p["k_lung"] * asp_drive * (1.0 - s["LUNG"]) - p["k_lung_rep"] * s["LUNG"]

    MAP_t = np.clip(p["MAP0"] - p["w_map_mus"] * mus_x
                    - p["w_map_solv"] * solv * op.solvent_cv_tox
                    - p["w_map_hyp"] * hyp, 25.0, 120.0)
    dMAP = (MAP_t - s["MAP"]) / p["tau_map"]

    seiz_drive = np.clip((cns_x - p["seiz_thresh"]) / (1.0 - p["seiz_thresh"]), 0.0, 1.0)
    dz_supp = 1.0 - 0.9 * C_dz / (p["dz_EC50"] + C_dz)
    dSEIZ = p["kseiz"] * seiz_drive * dz_supp * (1.0 - s["SEIZ"]) - p["kseiz_off"] * s["SEIZ"]

    dOPIDN = (p["kopidn"] * np.clip((s["N_aged"] - p["nte_thresh"]) / (1 - p["nte_thresh"]), 0, 1)
              * (1.0 - s["OPIDN"]) - p["krec_opidn"] * s["OPIDN"])

    shock = np.clip((65.0 - s["MAP"]) / 40.0, 0.0, 1.0)
    dHAZ = (p["h0"] + p["h_hypox"] * hyp ** 2
            + p["h_hyperc"] * np.clip((PaCO2 - 60.0) / 40.0, 0.0, 2.5) ** 2
            + p["h_shock"] * shock ** 2
            + p["h_seiz"] * s["SEIZ"]
            + p["h_lung"] * s["LUNG"]
            + p["h_vent"] * vent_on)

    d = np.zeros(NST)
    d[IDX["A_gut"]] = dA_gut
    d[IDX["A_th_c"]] = dA_th_c
    d[IDX["A_th_f"]] = dA_th_f
    d[IDX["A_ox_c"]] = dA_ox_c
    d[IDX["A_ox_t"]] = dA_ox_t
    d[IDX["A_solv_g"]] = dA_solv_g
    d[IDX["A_solv_c"]] = dA_solv_c
    d[IDX["LUNG"]] = dLUNG
    d[IDX["A_pam_c"]] = dA_pam_c
    d[IDX["A_pam_p"]] = dA_pam_p
    d[IDX["A_atr_c"]] = dA_atr_c
    d[IDX["Ce_atr"]] = dCe_atr
    d[IDX["A_dz"]] = dA_dz
    d[IDX["A_mg"]] = dA_mg
    d[IDX["SCAV"]] = dSCAV
    d[IDX["E_rbc"]], d[IDX["EP_rbc"]], d[IDX["EA_rbc"]] = dE_rbc, dEP_rbc, dEA_rbc
    d[IDX["E_mus"]], d[IDX["EP_mus"]], d[IDX["EA_mus"]] = dE_mus, dEP_mus, dEA_mus
    d[IDX["E_nmj"]], d[IDX["EP_nmj"]], d[IDX["EA_nmj"]] = dE_nmj, dEP_nmj, dEA_nmj
    d[IDX["E_cns"]], d[IDX["EP_cns"]], d[IDX["EA_cns"]] = dE_cns, dEP_cns, dEA_cns
    d[IDX["B_free"]], d[IDX["B_inh"]] = dB_free, dB_inh
    d[IDX["N_free"]], d[IDX["N_inh"]], d[IDX["N_aged"]] = dN_free, dN_inh, dN_aged
    d[IDX["ACh_m"]], d[IDX["ACh_n"]], d[IDX["ACh_b"]] = dACh_m, dACh_n, dACh_b
    d[IDX["Rn_des"]] = dRn_des
    d[IDX["Rm_down"]] = dRm_down
    d[IDX["SEC"]] = dSEC
    d[IDX["BT"]] = dBT
    d[IDX["HR"]] = dHR
    d[IDX["RESPD"]] = dRESPD
    d[IDX["MSTR"]] = dMSTR
    d[IDX["MAP"]] = dMAP
    d[IDX["SEIZ"]] = dSEIZ
    d[IDX["OPIDN"]] = dOPIDN
    d[IDX["HAZ"]] = dHAZ
    d[IDX["VTIME"]] = vent_on
    d[IDX["ATRCUM"]] = rate_atr * p["atr_MW"] / 1000.0     # umol/h -> mg/h
    d[IDX["AUC_EPn"]] = s["EP_nmj"]
    d[IDX["AUC_ox"]] = C_ox

    if aux is not None:
        aux.update(dict(Rm_down=s["Rm_down"], mus_raw=mus_raw, C_ox=C_ox, C_oxt=C_oxt, C_th=C_th, C_pam=C_ox_drug,
                        kr=kr, kinh=kinh, Om=Om, On=On, Ob=Ob, mus_x=mus_x,
                        nic_x=nic_x, cns_x=cns_x, Vcap=Vcap, VA=VA, PaO2=PaO2,
                        PaCO2=PaCO2, shunt=shunt, vent_on=vent_on,
                        rate_atr=rate_atr, solv=solv, hyp=hyp,
                        Ce_atr=s["Ce_atr"]))
    return d


# ==========================================================================
# 6. SIMULATION DRIVER (piecewise, with bolus events)
# ==========================================================================

def build_events(tx: Treatment, op: OPClass):
    """Return a sorted list of (time_h, state_index, amount) bolus events."""
    ev = []
    wt = P["WT"]
    # oxime
    if tx.oxime:
        MW = P["obi_MW"] if tx.oxime == "obi" else P["pam_MW"]
        if tx.oxime_bolus_mg:
            t = tx.oxime_start_h
            while t <= tx.oxime_start_h + tx.oxime_dur_h:
                ev.append((t, IDX["A_pam_c"], tx.oxime_bolus_mg * 1000.0 / MW,
                           tx.oxime_bolus_mg))
                t += tx.oxime_bolus_int_h
        else:
            ev.append((tx.oxime_start_h, IDX["A_pam_c"],
                       tx.oxime_load_mg_kg * wt * 1000.0 / MW,
                       tx.oxime_load_mg_kg * wt))
    # diazepam 10 mg at presentation, repeated at 6 h
    if tx.diazepam:
        for t in (0.6, 6.0, 12.0):
            ev.append((t, IDX["A_dz"], 10.0 * 1000.0 / P["dz_MW"], 0.0))
    # magnesium sulphate 4 g over the first hour
    if tx.magnesium:
        ev.append((1.0, IDX["A_mg"], 16.2, 0.0))     # 4 g MgSO4 = 16.2 mmol Mg
    # bioscavenger (plasma-derived human BChE, tetramer 340 kDa, 4 sites)
    if tx.scavenger_mg > 0:
        ev.append((1.0, IDX["SCAV"], tx.scavenger_mg / 1000.0 / 340000.0 * 1e6 * 4.0, 0.0))
    ev.sort(key=lambda e: e[0])
    return ev


def oxime_infusion_rate(t_h, tx: Treatment) -> float:
    """umol/h continuous oxime infusion."""
    if not tx.oxime or tx.oxime_bolus_mg:
        return 0.0
    if not (tx.oxime_start_h <= t_h <= tx.oxime_start_h + tx.oxime_dur_h):
        return 0.0
    MW = P["obi_MW"] if tx.oxime == "obi" else P["pam_MW"]
    return tx.oxime_inf_mg_kg_h * P["WT"] * 1000.0 / MW


def simulate(tx: Treatment, tmax=336.0, n=1345):
    op = OP_LIBRARY[tx.op]
    y0 = initial_state()
    dose_umol = tx.volume_ml * tx.conc_g_per_L / 1000.0 * 1e6 / op.mw
    y0[IDX["A_gut"]] = dose_umol
    y0[IDX["A_solv_g"]] = tx.volume_ml * op.solvent_frac

    ev = build_events(tx, op)
    tgrid = np.linspace(0.0, tmax, n)
    breaks = sorted(set([0.0, tmax] + [e[0] for e in ev if 0 < e[0] < tmax]
                        + [tx.oxime_start_h, tx.oxime_start_h + tx.oxime_dur_h]
                        + ([tx.charcoal_h, tx.charcoal_h + P["charcoal_dur"]]
                           if tx.charcoal_h is not None else [])))
    breaks = [b for b in breaks if 0.0 <= b <= tmax]

    def rhs(t, y):
        d = derivatives(t, y, op, tx)
        d[IDX["A_pam_c"]] += oxime_infusion_rate(t, tx)
        return d

    ys = [y0.copy()]
    ts = [0.0]
    y = y0.copy()
    oxcum = 0.0
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        for (te, idx, amt, mg) in ev:
            if abs(te - t0) < 1e-9:
                y[idx] += amt
                oxcum += mg
        y[IDX["OXCUM"]] = oxcum
        if t1 - t0 < 1e-9:
            continue
        pts = tgrid[(tgrid > t0) & (tgrid <= t1)]
        sol = solve_ivp(rhs, (t0, t1), y, method="LSODA",
                        t_eval=pts if len(pts) else None,
                        rtol=1e-6, atol=1e-9, max_step=0.5)
        if not sol.success:
            raise RuntimeError(f"integration failed in [{t0},{t1}]: {sol.message}")
        if len(pts):
            for j in range(sol.y.shape[1]):
                ts.append(sol.t[j]); ys.append(sol.y[:, j].copy())
        y = sol.y[:, -1].copy()
        # continuous oxime infusion contributes to the cumulative dose
        if tx.oxime and not tx.oxime_bolus_mg:
            lo = max(t0, tx.oxime_start_h)
            hi = min(t1, tx.oxime_start_h + tx.oxime_dur_h)
            if hi > lo:
                oxcum += tx.oxime_inf_mg_kg_h * P["WT"] * (hi - lo)
        y[IDX["OXCUM"]] = oxcum

    T = np.array(ts)
    Y = np.array(ys).T
    return T, Y, op


def outputs(T, Y, op, tx):
    """Recompute the algebraic outputs along a solved trajectory."""
    keys = ["C_ox", "C_oxt", "C_th", "C_pam", "kr", "kinh", "Om", "On", "Ob",
            "mus_x", "nic_x", "cns_x", "Vcap", "VA", "PaO2", "PaCO2", "shunt",
            "vent_on", "rate_atr", "solv", "hyp", "Ce_atr", "Rm_down", "mus_raw"]
    out = {k: np.zeros(len(T)) for k in keys}
    for j in range(len(T)):
        aux = {}
        derivatives(T[j], Y[:, j], op, tx, aux)
        for k in keys:
            out[k][j] = aux[k]
    return out


# ==========================================================================
# 7. CLOSED-FORM ANALYTICS -- the arithmetic the simulation must obey
# ==========================================================================

def omega(op: OPClass, oxime: str, C_oxon_nM: float) -> float:
    krmax = op.obi_krmax if oxime == "obi" else op.pam_krmax
    return krmax / (op.ki * C_oxon_nM)


def ceiling_activity(op: OPClass, oxime: str, C_oxon_nM: float) -> float:
    w = omega(op, oxime, C_oxon_nM)
    return w / (1.0 + w)


def achieved_activity(op: OPClass, oxime: str, C_oxon_nM: float, X_uM: float) -> float:
    krmax = op.obi_krmax if oxime == "obi" else op.pam_krmax
    KD = op.obi_KD if oxime == "obi" else op.pam_KD
    ks = LN2 / op.t12_spont_h
    kr = krmax * X_uM / (KD + X_uM)
    kin = op.ki * C_oxon_nM
    return (kr + ks) / (kr + ks + kin)


def C_crit(op: OPClass, oxime: str, target=0.30) -> float:
    """Oxon concentration above which NO oxime concentration reaches `target`."""
    krmax = op.obi_krmax if oxime == "obi" else op.pam_krmax
    return (1.0 - target) / target * krmax / op.ki


def X_required(op: OPClass, oxime: str, C_oxon_nM: float, target=0.30):
    """Oxime concentration needed to reach `target`; None if unreachable."""
    krmax = op.obi_krmax if oxime == "obi" else op.pam_krmax
    KD = op.obi_KD if oxime == "obi" else op.pam_KD
    ks = LN2 / op.t12_spont_h
    need = target / (1.0 - target) * op.ki * C_oxon_nM - ks
    if need <= 0:
        return 0.0
    if need >= krmax:
        return None
    return KD * need / (krmax - need)


def aging_ratchet(op: OPClass, oxime: str | None, X_uM: float) -> float:
    krmax = (op.obi_krmax if oxime == "obi" else op.pam_krmax) if oxime else 0.0
    KD = (op.obi_KD if oxime == "obi" else op.pam_KD) if oxime else 1.0
    kr = krmax * X_uM / (KD + X_uM) if oxime else 0.0
    ka = LN2 / op.t12_aging_h
    ks = LN2 / op.t12_spont_h
    return ka / (ka + ks + kr)


def steady_oxon(op: OPClass, dose_umol: float, pon1="RR") -> float:
    """Quasi-steady plasma oxon (nM) once the thion has distributed."""
    C_th = dose_umol / (op.V1 + op.Vf)
    bio = op.Vmax_bio * C_th / (op.Km_bio + C_th)
    CL = op.CL_oxon * (op.pon1_QQ_scale if pon1 == "QQ" else 1.0)
    return bio / CL


def critical_volume(op: OPClass, oxime: str, conc_g_per_L: float,
                    target=0.30, pon1="RR") -> float | None:
    """Ingested volume (mL) at which the oxime ceiling falls through `target`."""
    Cc = C_crit(op, oxime, target)
    CL = op.CL_oxon * (op.pon1_QQ_scale if pon1 == "QQ" else 1.0)
    need_bio = Cc * CL
    if need_bio >= op.Vmax_bio:
        return None                      # CYP cannot make that much oxon: safe
    C_th = op.Km_bio * need_bio / (op.Vmax_bio - need_bio)
    dose_umol = C_th * (op.V1 + op.Vf)
    grams = dose_umol * op.mw / 1e6
    return grams / conc_g_per_L * 1000.0


# ==========================================================================
# 8. SCENARIOS
# ==========================================================================

SCENARIOS = [
    Treatment(label="S01 CPF 50 mL - supportive only (atropine, ICU)",
              op="chlorpyrifos", volume_ml=50, oxime=None),
    Treatment(label="S02 CPF 50 mL - atropine + 2-PAM (WHO 30+8 mg/kg/h)",
              op="chlorpyrifos", volume_ml=50, oxime="pam"),
    Treatment(label="S03 CPF 50 mL - atropine + 2-PAM 1 g q6h bolus",
              op="chlorpyrifos", volume_ml=50, oxime="pam",
              oxime_bolus_mg=1000.0, oxime_bolus_int_h=6.0),
    Treatment(label="S04 CPF 50 mL - atropine + obidoxime (250 mg + 750 mg/24h)",
              op="chlorpyrifos", volume_ml=50, oxime="obi",
              oxime_load_mg_kg=250.0 / 60.0, oxime_inf_mg_kg_h=750.0 / 24.0 / 60.0),
    Treatment(label="S05 CPF 10 mL - atropine + 2-PAM (below the oxime ceiling)",
              op="chlorpyrifos", volume_ml=10, oxime="pam"),
    Treatment(label="S06 CPF 10 mL - supportive only",
              op="chlorpyrifos", volume_ml=10, oxime=None),
    Treatment(label="S07 CPF 200 mL - atropine + 2-PAM (massive ingestion)",
              op="chlorpyrifos", volume_ml=200, oxime="pam"),
    Treatment(label="S08 Dimethoate 50 mL - supportive only",
              op="dimethoate", volume_ml=50, conc_g_per_L=400, oxime=None),
    Treatment(label="S09 Dimethoate 50 mL - atropine + 2-PAM",
              op="dimethoate", volume_ml=50, conc_g_per_L=400, oxime="pam"),
    Treatment(label="S10 Dimethoate 50 mL - atropine + obidoxime",
              op="dimethoate", volume_ml=50, conc_g_per_L=400, oxime="obi",
              oxime_load_mg_kg=250.0 / 60.0, oxime_inf_mg_kg_h=750.0 / 24.0 / 60.0),
    Treatment(label="S11 CPF 50 mL - 2-PAM started at 0.5 h (early)",
              op="chlorpyrifos", volume_ml=50, oxime="pam", oxime_start_h=0.5),
    Treatment(label="S12 CPF 50 mL - 2-PAM started at 12 h (late)",
              op="chlorpyrifos", volume_ml=50, oxime="pam", oxime_start_h=12.0),
    Treatment(label="S13 CPF 50 mL - no ventilator available",
              op="chlorpyrifos", volume_ml=50, oxime="pam", ventilator=False,
              icu_suction=False),
    Treatment(label="S14 CPF 50 mL - standard (slow) atropine titration",
              op="chlorpyrifos", volume_ml=50, oxime="pam", atropine="standard"),
    Treatment(label="S15 CPF 50 mL - glycopyrrolate instead of atropine",
              op="chlorpyrifos", volume_ml=50, oxime="pam", atropine="glyco"),
    Treatment(label="S16 CPF 50 mL - + magnesium sulphate 4 g",
              op="chlorpyrifos", volume_ml=50, oxime="pam", magnesium=True),
    Treatment(label="S17 CPF 50 mL - + activated charcoal at 1 h",
              op="chlorpyrifos", volume_ml=50, oxime="pam", charcoal_h=1.0),
    Treatment(label="S18 CPF 50 mL - + activated charcoal at 6 h",
              op="chlorpyrifos", volume_ml=50, oxime="pam", charcoal_h=6.0),
    Treatment(label="S19 CPF 50 mL - + 1 g plasma BChE bioscavenger",
              op="chlorpyrifos", volume_ml=50, oxime="pam", scavenger_mg=1000.0),
    Treatment(label="S20 CPF 50 mL - PON1 Q192Q slow-hydrolyser",
              op="chlorpyrifos", volume_ml=50, oxime="pam", pon1="QQ"),
    Treatment(label="S21 Fenthion 50 mL - atropine + 2-PAM",
              op="fenthion", volume_ml=50, conc_g_per_L=500, oxime="pam"),
    Treatment(label="S22 Parathion 15 mL - atropine + obidoxime",
              op="parathion", volume_ml=15, conc_g_per_L=500, oxime="obi",
              oxime_load_mg_kg=250.0 / 60.0, oxime_inf_mg_kg_h=750.0 / 24.0 / 60.0),
    Treatment(label="S23 Parathion 15 mL - supportive only",
              op="parathion", volume_ml=15, conc_g_per_L=500, oxime=None),
]


def summarise(tx: Treatment):
    T, Y, op = simulate(tx)
    O = outputs(T, Y, op, tx)

    def at(h, arr):
        return float(np.interp(h, T, arr))

    rbc = Y[IDX["E_rbc"]]
    nmj = Y[IDX["E_nmj"]]
    aged = Y[IDX["EA_rbc"]]
    haz = Y[IDX["HAZ"]]
    vt = Y[IDX["VTIME"]]
    vent_flag = O["vent_on"] > 0.5
    return dict(
        label=tx.label,
        op=tx.op, volume_ml=tx.volume_ml, oxime=tx.oxime or "none",
        C_oxon_peak_nM=float(np.max(O["C_ox"])),
        C_oxon_24h_nM=at(24, O["C_ox"]),
        C_oxime_ss_uM=float(np.max(O["C_pam"])) if tx.oxime else 0.0,
        omega_24h=(op.pam_krmax if tx.oxime != "obi" else op.obi_krmax)
                  / max(1e-9, op.ki * at(24, O["C_ox"])),
        RBC_AChE_nadir=float(np.min(rbc)),
        RBC_AChE_6h=at(6, rbc), RBC_AChE_24h=at(24, rbc),
        RBC_AChE_72h=at(72, rbc), RBC_AChE_7d=at(168, rbc),
        NMJ_AChE_24h=at(24, nmj),
        aged_frac_24h=at(24, aged), aged_frac_72h=at(72, aged),
        peak_ACh_m=float(np.max(Y[IDX["ACh_m"]])),
        atropine_24h_mg=at(24, Y[IDX["ATRCUM"]]),
        atropine_total_mg=float(Y[IDX["ATRCUM"]][-1]),
        oxime_total_g=float(Y[IDX["OXCUM"]][-1]) / 1000.0,
        vent_needed=bool(vent_flag.any()),
        vent_start_h=float(T[vent_flag][0]) if vent_flag.any() else None,
        vent_hours=float(vt[-1]),
        min_PaO2=float(np.min(O["PaO2"])),
        min_MAP=float(np.min(Y[IDX["MAP"]])),
        max_Rn_des=float(np.max(Y[IDX["Rn_des"]])),
        IMS=bool(np.max(Y[IDX["Rn_des"]][T > 24]) > 0.35),
        max_SEIZ=float(np.max(Y[IDX["SEIZ"]])),
        NTE_aged_max=float(np.max(Y[IDX["N_aged"]])),
        OPIDN=float(np.max(Y[IDX["OPIDN"]])),
        mortality=float(1.0 - np.exp(-haz[-1])),
        mortality_24h=float(1.0 - np.exp(-at(24, haz))),
    )


# ==========================================================================
# 9. REPORT
# ==========================================================================

def hr(title, ch="="):
    return "\n" + ch * 78 + f"\n{title}\n" + ch * 78


def main():
    lines = []
    w = lines.append

    w("ACUTE ORGANOPHOSPHORUS INSECTICIDE SELF-POISONING -- QSP REFERENCE MODEL")
    w("Independent Python/scipy re-implementation used to verify op_mrgsolve_model.R")
    w(f"{NST} ODE states; {len(SCENARIOS)} therapeutic scenarios")

    # ---------------------------------------------------------------- part 1
    w(hr("1. THE OXIME SUFFICIENCY NUMBER  Omega = k_r_max / (k_i * C_oxon)"))
    w("Oxime reactivation SATURATES (Worek: k_r2 = k_r_max/(K_D + X)), so the")
    w("best achievable free-enzyme fraction at UNLIMITED oxime is Omega/(1+Omega).")
    w("")
    w(f"{'OP (oxon)':<26}{'oxime':<10}{'k_i':>9}{'k_r_max':>10}"
      f"{'C_crit(30%)':>13}{'C_crit(50%)':>13}")
    w(f"{'':<26}{'':<10}{'nM-1h-1':>9}{'h-1':>10}{'nM':>13}{'nM':>13}")
    w("-" * 78)
    for op in (CHLORPYRIFOS, PARATHION, DIMETHOATE, FENTHION):
        for ox in ("pam", "obi"):
            krm = op.obi_krmax if ox == "obi" else op.pam_krmax
            w(f"{op.name + ' (' + op.subclass + ')':<26}"
              f"{'2-PAM' if ox == 'pam' else 'obidoxime':<10}"
              f"{op.ki:>9.4f}{krm:>10.1f}"
              f"{C_crit(op, ox, 0.30):>13.0f}{C_crit(op, ox, 0.50):>13.0f}")

    # ---------------------------------------------------------------- part 2
    w(hr("2. CRITICAL INGESTED VOLUME (the dose at which the oxime runs out of"
         " ceiling)"))
    w("CYP desulfuration is saturable, so above C_thion >> Km the oxon plateau")
    w("is dose-INDEPENDENT at Vmax/CL_oxon.  The oxime ceiling therefore becomes")
    w("a CONSTANT for every large ingestion, no matter how large.")
    w("")
    w("'ceiling' is the best any dose of that oxime could do; 'achieved' is what")
    w("the clinically attainable concentration (2-PAM 130 uM, obidoxime 20 uM)")
    w("actually delivers.  The gap between the two columns is the reason a")
    w("higher oxime dose is not the answer.")
    w("")
    w(f"{'OP / formulation':<32}{'oxime':<11}{'V_crit(30%)':>12}"
      f"{'plateau oxon':>14}{'ceiling':>9}{'achieved':>10}")
    w(f"{'':<32}{'':<11}{'mL':>12}{'nM':>14}{'%':>9}{'%':>10}")
    w("-" * 88)
    for op, conc in ((CHLORPYRIFOS, 200.0), (PARATHION, 500.0),
                     (DIMETHOATE, 400.0), (FENTHION, 500.0)):
        plateau = op.Vmax_bio / op.CL_oxon
        for ox in ("pam", "obi"):
            vc = critical_volume(op, ox, conc)
            vcs = f"{vc:.1f}" if vc is not None else "never"
            Xclin = 130.0 if ox == "pam" else 20.0
            w(f"{op.name + ' ' + str(int(conc / 10)) + '% EC':<32}"
              f"{'2-PAM' if ox == 'pam' else 'obidoxime':<11}{vcs:>12}"
              f"{plateau:>14.0f}{100 * ceiling_activity(op, ox, plateau):>9.1f}"
              f"{100 * achieved_activity(op, ox, plateau, Xclin):>10.1f}")

    # ---------------------------------------------------------------- part 3
    w(hr("3. THE AGING RATCHET  phi = k_a/(k_a + k_s + k_r)"))
    w("phi is the probability that ONE inhibition event ends irreversibly.")
    w("Dimethyl OPs age 9x faster than diethyl but ALSO reactivate spontaneously")
    w("110x faster -- the classic teaching that 'dimethyl OPs are hopeless' is")
    w("therefore NOT a statement about phi, it is a statement about how long the")
    w("enzyme stays inhibited, i.e. about C_oxon.")
    w("")
    w(f"{'OP':<16}{'k_a h-1':>10}{'k_s h-1':>10}{'phi no oxime':>14}"
      f"{'phi 2-PAM 130uM':>18}{'phi obi 20uM':>14}")
    w("-" * 82)
    for op in (CHLORPYRIFOS, PARATHION, DIMETHOATE, FENTHION):
        w(f"{op.name:<16}{LN2 / op.t12_aging_h:>10.4f}{LN2 / op.t12_spont_h:>10.4f}"
          f"{aging_ratchet(op, None, 0):>14.3f}"
          f"{aging_ratchet(op, 'pam', 130.0):>18.3f}"
          f"{aging_ratchet(op, 'obi', 20.0):>14.3f}")
    w("")
    w("Cumulative irreversible loss is an INTEGRAL, not an event:")
    w("      aged(T) = 1 - exp(-k_a * INTEGRAL f_EP dt)")
    for op in (CHLORPYRIFOS, DIMETHOATE):
        ka = LN2 / op.t12_aging_h
        for hrs in (6, 12, 24, 72):
            w(f"   {op.name:<14} f_EP=1 for {hrs:>3} h  ->  aged = "
              f"{100 * (1 - np.exp(-ka * hrs)):>5.1f} %")

    # ---------------------------------------------------------------- part 4
    w(hr("4. OXIME CONCENTRATION ACTUALLY REQUIRED, AND WHETHER IT EXISTS"))
    w("Target: 30% free AChE (the threshold below which muscarinic signs appear).")
    w("")
    w(f"{'OP':<14}{'C_oxon nM':>11}{'X* 2-PAM uM':>14}{'X* obidoxime uM':>18}"
      f"   {'achievable?':<16}")
    w("Achievable = reachable at 2-PAM 130 uM (WHO regimen) or obidoxime 20 uM.")
    w("-" * 76)
    for op, cs in ((CHLORPYRIFOS, (50, 200, 1000)),
                   (DIMETHOATE, (1000, 5000, 10000)),
                   (PARATHION, (50, 200, 700))):
        for c in cs:
            xp = X_required(op, "pam", c)
            xo = X_required(op, "obi", c)
            sp = "impossible" if xp is None else f"{xp:.0f}"
            so = "impossible" if xo is None else f"{xo:.0f}"
            # 2-PAM 130 uM and obidoxime 20 uM are the clinically attained values
            ok = ("yes" if (xp is not None and xp <= 130) else
                  ("obidoxime only" if (xo is not None and xo <= 20) else "no"))
            w(f"{op.name:<14}{c:>11}{sp:>14}{so:>18}   {ok:<16}")

    # ---------------------------------------------------------------- part 5
    w(hr("5. THE ATROPINE REQUIREMENT IS A HYPERBOLA IN ENZYME ACTIVITY"))
    w("Cleft ACh at steady state = (1+leak)/(E+leak); competitive antagonism then")
    w("makes the atropine concentration needed to restore baseline occupancy")
    w("        Ce* = KB * (ACh/ACh_0 - 1)")
    w("so the atropine requirement scales as 1/E, not as (1-E).")
    w("")
    w(f"{'AChE %':>8}{'cleft ACh x normal':>21}{'Ce* nM':>10}"
      f"{'bolus-equivalent mg':>22}{'mg/h to hold':>14}")
    w("-" * 76)
    for E in (1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.0):
        ach = (1 + P["ach_leak"]) / (E + P["ach_leak"])
        ce = P["KB_atr"] * (ach - 1.0)
        mg = ce * P["atr_V"] / 1000.0 * P["atr_MW"] / 1000.0
        mgh = mg * P["atr_CL"] / P["atr_V"]
        w(f"{100 * E:>8.0f}{ach:>21.2f}{ce:>10.1f}{mg:>22.2f}{mgh:>14.2f}")

    # ---------------------------------------------------------------- part 6
    w(hr("6. STOICHIOMETRIC BIOSCAVENGER: AN ARITHMETIC IMPOSSIBILITY FOR"
         " INSECTICIDES"))
    w("Plasma-derived human BChE binds oxon 1:1.  Nerve-agent prophylaxis works")
    w("because the agent dose is micromoles; an insecticide ingestion is")
    w("millimoles of PARENT that keeps generating oxon for days.")
    w("")
    for ml, conc, op in ((50, 200.0, CHLORPYRIFOS), (50, 400.0, DIMETHOATE)):
        dose_umol = ml * conc / 1000.0 * 1e6 / op.mw
        for scav_mg in (200.0, 1000.0):
            sites = scav_mg / 1000.0 / 340000.0 * 1e6 * 4.0
            w(f"   {op.name:<13} {ml} mL of {int(conc/10)}% EC = "
              f"{dose_umol/1000:>7.1f} mmol parent;  {scav_mg:>6.0f} mg BChE = "
              f"{sites:>6.2f} umol sites  ->  shortfall {dose_umol/sites:>8.0f}-fold")

    # ---------------------------------------------------------------- part 7
    w(hr("7. SCENARIO SIMULATIONS (51-state ODE system)"))
    results = []
    for tx in SCENARIOS:
        r = summarise(tx)
        results.append(r)
    w(f"{'scenario':<52}{'AChE':>7}{'AChE':>7}{'aged':>7}{'atr':>7}"
      f"{'vent':>7}{'mort':>7}")
    w(f"{'':<52}{'nadir':>7}{'72h':>7}{'72h':>7}{'mg':>7}{'h':>7}{'%':>7}")
    w("-" * 96)
    for r in results:
        w(f"{r['label'][:51]:<52}{100*r['RBC_AChE_nadir']:>7.1f}"
          f"{100*r['RBC_AChE_72h']:>7.1f}{100*r['aged_frac_72h']:>7.1f}"
          f"{r['atropine_total_mg']:>7.0f}{r['vent_hours']:>7.0f}"
          f"{100*r['mortality']:>7.1f}")

    # ---------------------------------------------------------------- part 8
    w(hr("8. CONTRASTS THE MODEL PRODUCES WITHOUT BEING TOLD TO"))
    d = {r["label"][:3]: r for r in results}

    def diff(a, b, key, fmt="{:.1f}"):
        return fmt.format(d[a][key]), fmt.format(d[b][key])

    w("(a) Oxime benefit is a function of INGESTED DOSE, not of the drug:")
    for pair, tag in ((("S01", "S02"), "50 mL chlorpyrifos"),
                      (("S06", "S05"), "10 mL chlorpyrifos"),
                      (("S23", "S22"), "15 mL parathion"),
                      (("S08", "S09"), "50 mL dimethoate")):
        a, b = pair
        w(f"    {tag:<22} AChE72h {100*d[a]['RBC_AChE_72h']:>5.1f}% -> "
          f"{100*d[b]['RBC_AChE_72h']:>5.1f}% | mortality "
          f"{100*d[a]['mortality']:>5.1f}% -> {100*d[b]['mortality']:>5.1f}%")
    w("")
    w("(b) Route of oxime administration matters more than the choice of oxime")
    w(f"    2-PAM infusion (S02)  Css {d['S02']['C_oxime_ss_uM']:.0f} uM  "
      f"AChE72h {100*d['S02']['RBC_AChE_72h']:.1f}%  mort {100*d['S02']['mortality']:.1f}%")
    w(f"    2-PAM q6h bolus (S03) Cmax {d['S03']['C_oxime_ss_uM']:.0f} uM "
      f"AChE72h {100*d['S03']['RBC_AChE_72h']:.1f}%  mort {100*d['S03']['mortality']:.1f}%")
    w(f"    obidoxime      (S04)  Cmax {d['S04']['C_oxime_ss_uM']:.0f} uM  "
      f"AChE72h {100*d['S04']['RBC_AChE_72h']:.1f}%  mort {100*d['S04']['mortality']:.1f}%")
    w("")
    w("(c) The ventilator is worth more than every drug in the list:")
    w(f"    with ventilator (S02)  mortality {100*d['S02']['mortality']:.1f}%")
    w(f"    without         (S13)  mortality {100*d['S13']['mortality']:.1f}%")
    w("")
    w("(d) Atropine titration speed, same total pharmacology:")
    w(f"    rapid doubling  (S02)  atropine {d['S02']['atropine_24h_mg']:.0f} mg/24h, "
      f"mortality {100*d['S02']['mortality']:.1f}%")
    w(f"    slow ad-hoc     (S14)  atropine {d['S14']['atropine_24h_mg']:.0f} mg/24h, "
      f"mortality {100*d['S14']['mortality']:.1f}%")
    w("")
    w("(e) A peripherally restricted antimuscarinic cannot protect the brain:")
    w(f"    atropine        (S02)  peak seizure index {d['S02']['max_SEIZ']:.2f}, "
      f"mortality {100*d['S02']['mortality']:.1f}%")
    w(f"    glycopyrrolate  (S15)  peak seizure index {d['S15']['max_SEIZ']:.2f}, "
      f"mortality {100*d['S15']['mortality']:.1f}%")
    w("")
    w("(f) Decontamination is a race against a first-order absorption constant:")
    w(f"    charcoal at 1 h (S17)  AChE nadir {100*d['S17']['RBC_AChE_nadir']:.1f}%, "
      f"mortality {100*d['S17']['mortality']:.1f}%")
    w(f"    charcoal at 6 h (S18)  AChE nadir {100*d['S18']['RBC_AChE_nadir']:.1f}%, "
      f"mortality {100*d['S18']['mortality']:.1f}%")
    w("")
    w("(g) PON1 genotype scales the oxon exposure, hence the aged fraction:")
    w(f"    PON1 R192R (S02)  oxon24h {d['S02']['C_oxon_24h_nM']:.0f} nM, "
      f"aged72h {100*d['S02']['aged_frac_72h']:.1f}%, mort {100*d['S02']['mortality']:.1f}%")
    w(f"    PON1 Q192Q (S20)  oxon24h {d['S20']['C_oxon_24h_nM']:.0f} nM, "
      f"aged72h {100*d['S20']['aged_frac_72h']:.1f}%, mort {100*d['S20']['mortality']:.1f}%")
    w("")
    w("(h) The bioscavenger changes nothing, as section 6 predicts:")
    w(f"    no scavenger (S02) AChE nadir {100*d['S02']['RBC_AChE_nadir']:.2f}%")
    w(f"    1 g BChE     (S19) AChE nadir {100*d['S19']['RBC_AChE_nadir']:.2f}%")
    w("")
    w("(i) Lipophilic dimethyl OP (fenthion): late, nicotinic, ventilator-hungry")
    w(f"    S21 vent hours {d['S21']['vent_hours']:.0f}, "
      f"peak nAChR block {d['S21']['max_Rn_des']:.2f}, "
      f"IMS {d['S21']['IMS']}, mortality {100*d['S21']['mortality']:.1f}%")
    w(f"    S02 vent hours {d['S02']['vent_hours']:.0f}, "
      f"peak nAChR block {d['S02']['max_Rn_des']:.2f}, "
      f"IMS {d['S02']['IMS']}")
    w("")
    w("(j) RBC AChE is a bad recovery marker because red cells cannot resynthesise:")
    for tag in ("S02", "S05"):
        T, Y, op = simulate([s for s in SCENARIOS if s.label.startswith(tag)][0])
        w(f"    {tag}: day-7 RBC AChE {100*np.interp(168, T, Y[IDX['E_rbc']]):.1f}% vs "
          f"day-7 NMJ AChE {100*np.interp(168, T, Y[IDX['E_nmj']]):.1f}% "
          f"(muscle enzyme is resynthesised, red cells are not)")

    # ---------------------------------------------------------------- part 9
    w(hr("9. DOSE-RESPONSE SWEEP: WHERE THE OXIME STOPS WORKING"))
    w(f"{'volume mL':>10}{'plateau oxon nM':>18}{'Omega(2-PAM)':>14}"
      f"{'ceiling %':>11}{'AChE72h sup':>13}{'AChE72h PAM':>13}{'d mort %':>10}")
    w("-" * 89)
    sweep = []
    for ml in (2, 5, 10, 20, 35, 50, 100, 200):
        sup = summarise(Treatment(label=f"sweep {ml} sup", op="chlorpyrifos",
                                  volume_ml=ml, oxime=None))
        pam = summarise(Treatment(label=f"sweep {ml} pam", op="chlorpyrifos",
                                  volume_ml=ml, oxime="pam"))
        dose_umol = ml * 200.0 / 1000.0 * 1e6 / CHLORPYRIFOS.mw
        c = steady_oxon(CHLORPYRIFOS, dose_umol)
        om = omega(CHLORPYRIFOS, "pam", c)
        w(f"{ml:>10}{c:>18.0f}{om:>14.2f}{100*om/(1+om):>11.1f}"
          f"{100*sup['RBC_AChE_72h']:>13.1f}{100*pam['RBC_AChE_72h']:>13.1f}"
          f"{100*(sup['mortality']-pam['mortality']):>10.1f}")
        sweep.append(dict(volume_ml=ml, oxon_nM=c, omega=om,
                          ceiling=om / (1 + om),
                          AChE72_supportive=sup["RBC_AChE_72h"],
                          AChE72_pam=pam["RBC_AChE_72h"],
                          mort_supportive=sup["mortality"],
                          mort_pam=pam["mortality"]))

    # --------------------------------------------------------------- part 10
    w(hr("10. VIRTUAL TRIAL: what a randomised oxime trial would have measured"))
    w("A pragmatic trial enrols whoever arrives.  Ingested volume is log-normal")
    w("and OP class is mixed; the model is run patient by patient and the two")
    w("arms are compared on the endpoint the real trials used (death).")
    rng = np.random.default_rng(20260805)
    pop = []
    n_pat = 60
    for i in range(n_pat):
        vol = float(np.clip(np.exp(rng.normal(np.log(28.0), 0.85)), 1.0, 250.0))
        opn = "chlorpyrifos" if rng.random() < 0.55 else "dimethoate"
        conc = 200.0 if opn == "chlorpyrifos" else 400.0
        pon = "QQ" if rng.random() < 0.25 else "RR"
        base = dict(op=opn, volume_ml=vol, conc_g_per_L=conc, pon1=pon)
        a = summarise(Treatment(label=f"P{i}sup", oxime=None, **base))
        b = summarise(Treatment(label=f"P{i}pam", oxime="pam", **base))
        pop.append(dict(id=i, op=opn, volume_ml=vol, pon1=pon,
                        mort_placebo=a["mortality"], mort_pam=b["mortality"],
                        AChE72_placebo=a["RBC_AChE_72h"], AChE72_pam=b["RBC_AChE_72h"],
                        vent_placebo=a["vent_hours"], vent_pam=b["vent_hours"]))
    mp = np.mean([q["mort_placebo"] for q in pop])
    mo = np.mean([q["mort_pam"] for q in pop])
    w("")
    w(f"   n = {n_pat};  55% chlorpyrifos / 45% dimethoate;  median ingestion "
      f"{np.median([q['volume_ml'] for q in pop]):.0f} mL")
    w(f"   mortality  placebo {100*mp:.1f}%   pralidoxime {100*mo:.1f}%   "
      f"RR {mo/mp:.2f}   ARR {100*(mp-mo):.1f} points")
    small = [q for q in pop if q["volume_ml"] <= 15]
    large = [q for q in pop if q["volume_ml"] > 15]
    for tag, grp in (("<= 15 mL", small), ("> 15 mL", large)):
        if not grp:
            continue
        a = np.mean([q["mort_placebo"] for q in grp])
        b = np.mean([q["mort_pam"] for q in grp])
        w(f"   subgroup {tag:<9} n={len(grp):>3}  placebo {100*a:>5.1f}%  "
          f"2-PAM {100*b:>5.1f}%  RR {b/max(a,1e-9):.2f}")
    for tag, opn in (("chlorpyrifos", "chlorpyrifos"), ("dimethoate", "dimethoate")):
        grp = [q for q in pop if q["op"] == opn]
        a = np.mean([q["mort_placebo"] for q in grp])
        b = np.mean([q["mort_pam"] for q in grp])
        w(f"   subgroup {tag:<12} n={len(grp):>3}  placebo {100*a:>5.1f}%  "
          f"2-PAM {100*b:>5.1f}%  RR {b/max(a,1e-9):.2f}")
    w("")
    w("   The whole-trial RR sits close to 1 while a real, large effect exists in")
    w("   the low-dose diethyl stratum.  Dilution by patients above the oxime")
    w("   ceiling is sufficient to produce a null trial: no assumption about the")
    w("   oxime being inactive is required.")

    # --------------------------------------------------------------- outputs
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(HERE, "op_reference_output.txt"), "w") as f:
        f.write(txt)
    with open(os.path.join(HERE, "op_scenario_results.json"), "w") as f:
        json.dump(dict(scenarios=results, dose_sweep=sweep), f, indent=1)
    with open(os.path.join(HERE, "op_population_results.json"), "w") as f:
        json.dump(dict(n=n_pat, mort_placebo=mp, mort_pam=mo,
                       patients=pop), f, indent=1)
    print(txt)


if __name__ == "__main__":
    main()
