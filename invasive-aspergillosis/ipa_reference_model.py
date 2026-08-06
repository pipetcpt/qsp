#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ipa_reference_model.py
======================
Independent Python/scipy re-implementation of the invasive pulmonary
aspergillosis (IPA) QSP model that is written in mrgsolve C++ in
`ipa_mrgsolve_model.R`.

WHY A SECOND IMPLEMENTATION EXISTS
----------------------------------
Every quantitative claim in `README.md` is produced by THIS file.  The
mrgsolve model and this file were written from the same equation sheet but
in two different languages, and any number that the two disagree on is a
defect in one of them.  Writing the system twice is the only cheap way to
catch the class of error that a single implementation cannot see: an
amount used where a concentration was meant, a rate constant with the wrong
sign, a Hill term that saturates at the wrong end.

THE MODEL IN ONE PARAGRAPH
--------------------------
Hyphal biomass grows exponentially and is removed by exactly two things:
recruited neutrophils and a fungicidal drug.  Triazoles are NOT in that
list -- they multiply the growth term instead (ergosterol depletion slows
elongation; it does not lyse hyphae).  Therefore in a host with no
neutrophils, an azole with maximal growth inhibition Imax < 1 leaves
dB/dt = kgrow*(1-Imax)*B > 0 and the burden still rises, only slower.  That
single structural asymmetry, plus a perfusion term that falls as the lesion
angioinvades (so the drug that must reach the hyphae is delivered through
the vessels the hyphae are destroying), generates most of what the clinical
literature reports: the dominance of neutrophil recovery over drug choice,
a critical time-to-treatment beyond which no regimen sterilises, and a
galactomannan trajectory that moves for reasons unrelated to burden.

STATE VECTOR: 53 ODEs.  Units: time h, concentration mg/L, fungal burden in
CFU-equivalents (CFUe), ANC in 10^9/L, cytokines in pg/mL, GM in index units.
"""

import json
import math
import sys
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# 1. STATE VECTOR
# ----------------------------------------------------------------------------
STATES = [
    # --- voriconazole (VRC): oral depot, 2-cpt plasma, ELF, brain -----------
    "VRC_gut", "VRC_c", "VRC_p", "VRC_elf", "VRC_brain",
    # --- isavuconazole (ISA) ------------------------------------------------
    "ISA_gut", "ISA_c", "ISA_p", "ISA_elf",
    # --- posaconazole (POS) -------------------------------------------------
    "POS_gut", "POS_c", "POS_p", "POS_elf",
    # --- liposomal amphotericin B (AMB) ------------------------------------
    "AMB_c", "AMB_p", "AMB_elf", "AMB_brain",
    # --- echinocandin (ECH: anidulafungin / caspofungin) --------------------
    "ECH_c", "ECH_p", "ECH_elf",
    # --- co-medication ------------------------------------------------------
    "TAC_c",      # tacrolimus, ng/mL  (CYP3A4 victim of the azole)
    "GCS_c",      # corticosteroid, prednisone-equivalent mg/L
    "GCSF_c",     # filgrastim, ng/mL
    # --- fungus -------------------------------------------------------------
    "CON",        # resting/swollen conidia, CFUe
    "GERM",       # germlings (committed, not yet hyphal), CFUe
    "HYW",        # hyphal biomass, azole-susceptible (wild type), CFUe
    "HYR",        # hyphal biomass, azole-resistant (cyp51A TR34/L98H), CFUe
    "BRB",        # cerebral fungal burden, CFUe
    "BIO",        # galactosaminogalactan / hyphal matrix, relative units
    # --- lesion / vasculature ----------------------------------------------
    "VLES",       # lesion volume, mL
    "NEC",        # infarcted (non-perfused) lesion volume, mL
    "PERF",       # lesion perfusion index, 0-1
    "ANG",        # angioinvasion burden / vessel wall damage, relative
    # --- host defence -------------------------------------------------------
    "MAR",        # marrow granulocyte reserve, relative
    "ANC",        # blood absolute neutrophil count, 10^9/L
    "NLES",       # neutrophils recruited into lesion, relative
    "MAC",        # functional alveolar macrophages, relative
    "IFNG",       # IFN-gamma / Th1 axis, pg/mL
    "IL6",        # IL-6, pg/mL
    "CRP",        # C-reactive protein, mg/L
    "IL10",       # IL-10 (immunoparalysis arm), pg/mL
    # --- biomarkers ---------------------------------------------------------
    "GMles",      # galactomannan in lesion, index units x mL
    "GMser",      # serum galactomannan optical density index
    "BDG",        # serum (1,3)-beta-D-glucan, pg/mL
    "PCRs",       # serum Aspergillus DNA, copies/mL
    # --- organ / toxicity ---------------------------------------------------
    "ALT",        # U/L
    "SCR",        # serum creatinine, mg/dL
    "KMG",        # renal K+/Mg2+ wasting index, relative
    "QTC",        # QTc prolongation above baseline, ms
    "NEUROE",     # visual / neuropsychiatric effect compartment, relative
    "PFR",        # PaO2/FiO2 ratio, mmHg
    "TEMP",       # core temperature, degC
    "HAZ",        # cumulative mortality hazard (dimensionless)
]
IDX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)

# ----------------------------------------------------------------------------
# 2. PARAMETERS
# ----------------------------------------------------------------------------
P = OrderedDict(
    # ===================== voriconazole =====================================
    # 2-cpt, Michaelis-Menten hepatic elimination (CYP2C19/CYP3A4).
    # Calibrated to 4 mg/kg IV q12h -> Cmin ~2-3 mg/L, AUC24 ~40-50 mg*h/L
    # in a CYP2C19 normal metaboliser (Pascual 2008; Dolton 2012).
    VRC_ka=1.10, VRC_F=0.96, VRC_Vc=50.0, VRC_Vp=150.0, VRC_Q=25.0,
    VRC_Vmax=47.8,          # mg/h, CYP2C19 path (saturable, genotype-scaled)
    VRC_Km=3.00,            # mg/L  -- nonlinearity: dose 2x -> AUC >2x
    VRC_CLlin=2.50,         # L/h, CYP3A4/FMO3 path (linear, genotype-neutral)
    VRC_fu=0.42,            # 58% protein bound
    VRC_Kp_elf=7.1,         # ELF:plasma total ~7 (Crandon 2009)
    VRC_keq_elf=0.60,       # 1/h
    VRC_Kp_brain=0.50,      # CSF:plasma ~0.5 (Lutsar 2003)
    VRC_keq_brain=0.10,
    VRC_infl_Emax=0.62,     # CRP-driven suppression of CYP2C19 (Veringa 2017)
    VRC_infl_EC50=105.0,    # mg/L CRP
    # ===================== isavuconazole ====================================
    # ISAVUCONAZONIUM -> isavuconazole. 200 mg q8h x6 then 200 mg q24h.
    # Target: Cav ~4 mg/L, AUC24 ~100 mg*h/L (Desai 2016; SECURE).
    ISA_ka=0.90, ISA_F=0.98, ISA_Vc=90.0, ISA_Vp=360.0, ISA_Q=8.0,
    ISA_CL=1.75, ISA_fu=0.010,
    ISA_Kp_elf=1.20, ISA_keq_elf=0.50, ISA_Kp_brain=0.11,
    # ===================== posaconazole =====================================
    # Delayed-release tablet 300 mg (loading bid d1) -> Cav ~1.6 mg/L.
    POS_ka=0.30, POS_F=0.54, POS_Vc=190.0, POS_Vp=290.0, POS_Q=6.0,
    POS_CL=4.50, POS_fu=0.015,
    POS_Kp_elf=1.60, POS_keq_elf=0.40, POS_Kp_brain=0.03,
    # ===================== liposomal amphotericin B =========================
    # 3 mg/kg -> Cmax ~83 mg/L, AUC24 ~555 mg*h/L (AmBisome label; Walsh 2001).
    AMB_Vc=2.00, AMB_Vp=10.0, AMB_Q=0.90, AMB_CL=0.38,
    AMB_Kp_elf=0.42, AMB_keq_elf=0.06,     # slow lung accumulation
    AMB_Kp_brain=0.045, AMB_keq_brain=0.02,
    # ===================== echinocandin =====================================
    # Anidulafungin 200 mg load / 100 mg qd: Cmax ~7, AUC24 ~110 mg*h/L.
    ECH_Vc=30.0, ECH_Vp=20.0, ECH_Q=1.0, ECH_CL=0.90,
    ECH_Kp_elf=0.50, ECH_keq_elf=0.10, ECH_Kp_brain=0.02,
    # ===================== co-medication ====================================
    TAC_Vc=1000.0, TAC_CL=30.0,      # CL/F; 3 mg bid -> trough ~8 ng/mL
    TAC_KI_VRC=1.30,        # mg/L voriconazole giving 50% CYP3A4 inhibition
    TAC_KI_ISA=10.0,        # isavuconazole is a much weaker 3A4 inhibitor
    TAC_KI_POS=0.45,        # posaconazole is the strongest of the three
    # GCS_c is an EFFECT compartment, not plasma prednisolone: its k_out is set
    # by the ~22 h genomic effect half-life of glucocorticoids on the phagocyte
    # oxidative burst, not by the 2-3 h plasma half-life of prednisolone.
    GCS_CL=1.40, GCS_Vc=45.0,
    GCSF_CL=0.55, GCSF_Vc=4.0,
    # ===================== fungal growth ====================================
    KGROW=0.095,            # 1/h -> 0.99 log10/day net in a naive neutropenic host
    BMAX=3.0e10,            # CFUe carrying capacity of a lobe
    KGERM=0.085,            # conidium -> germling, 1/h
    KMAT=0.220,             # germling -> hypha, 1/h
    KMUT=2.0e-8,            # cyp51A resistance per replication
    KGROW_BRAIN=0.100,      # slower growth in brain parenchyma
    BMAX_BRAIN=1.0e9,       # carrying capacity of cerebral parenchyma (CFUe)
    # ===================== host killing =====================================
    KHOST_MAX=0.320,        # 1/h maximal RECRUITED neutrophil killing
    KHOST_BASE=0.220,       # 1/h constitutive alveolar surveillance at normal ANC
    KN50=0.55,              # NLES giving half-maximal killing
    KMAC_CON=0.150,         # macrophage conidial killing, 1/h per unit MAC
    GCS_IC50=0.900,         # mg/L prednisone-eq: 50% loss of phagocyte burst
    BIO_SHIELD=0.55,        # max fraction of neutrophil killing blocked by GAG
    KBIO=0.60,
    KBIO_ON=2.2e-11, KBIO_OFF=0.020,
    # ===================== drug PD ==========================================
    # Azoles: GROWTH inhibition (fungistatic). Imax < 1 by construction.
    IMAX_AZOLE=0.955,
    VRC_EC50_MIC=3.78,      # EC50_site = coef * MIC (mg/L)
    ISA_EC50_MIC=0.525,
    POS_EC50_MIC=2.33,
    HILL_AZ=2.0,
    # Amphotericin B: KILLING (fungicidal).
    AMB_EMAX=0.185,         # 1/h  (calibrated to AmBiLoad 12-week survival)
    AMB_EC50_MIC=3.40,      # mg/L per (MIC/1.0)
    HILL_AMB=2.4,
    # Echinocandin: apical growth inhibition only + paradoxical attenuation.
    ECH_IMAX=0.520,
    ECH_EC50=0.30,          # mg/L
    ECH_CPAR=12.0,          # mg/L above which the paradoxical effect appears
    HILL_ECH=1.6,
    # MICs (mg/L)
    MIC_VRC=0.5, MIC_ISA=1.0, MIC_POS=0.12, MIC_AMB=1.0,
    MIC_VRC_R=8.0, MIC_ISA_R=8.0, MIC_POS_R=0.75, MIC_AMB_R=1.0,
    # ===================== lesion / perfusion ===============================
    KVLES=2.0e-10,          # mL per CFUe per h -> ~1e8 CFUe per mL of lesion
    KRESOLVE_BASE=0.020,    # 1/h passive resolution
    KRESOLVE=0.010,         # 1/h extra, neutrophil-dependent
    KANG=1.2e-12,           # angioinvasion index per CFUe per h
    KANG_OFF=0.0090,
    KPERF_LOSS=0.020,       # max perfusion loss rate (1/h) -- now saturating
    KANG50=0.35,            # angioinvasion index at half-maximal perfusion loss
    KPERF_REC=0.012,
    PERF_GAMMA=1.0,         # C_site = C_elf * (DELIV_FLOOR + (1-DELIV_FLOOR)*PERF^gamma)
    DELIV_FLOOR=0.25,       # diffusive delivery to the viable rim, perfusion-independent
    KNEC=0.0025, KNEC_OFF=0.0018,
    # ===================== haematopoiesis ===================================
    KPROL=0.052, KTR=0.052, KOUT_ANC=0.056,
    ANC_SS=4.0,
    KMARG=0.011,            # margination into the lesion, per h
    KMARG_EC50=1.0e5,       # CFUe giving half-maximal neutrophil recruitment
    GCSF_EMAX=3.6, GCSF_EC50=6.0,
    KNLES_IN=0.020, KNLES_OUT=0.030,
    # ===================== innate / cytokines ===============================
    MAC_SS=1.0, KMAC_TURN=0.010,
    KIFNG_ON=0.9, KIFNG_OFF=0.10,
    KIL6_MAX=40.0, KIL6_EC50=2.0e9, KIL6_OFF=0.16,
    KCRP_ON=0.014, KCRP_OFF=0.0145,   # CRP t1/2 ~19 h
    KIL10_ON=0.25, KIL10_OFF=0.13,
    # ===================== biomarkers =======================================
    # GM is released by ACTIVELY ELONGATING hyphae and by lysis; it is not
    # proportional to standing biomass.  This is the modelling decision that
    # makes GM and burden dissociate under a fungistatic drug.
    KGM_GROWTH=3.5e-10,      # index*mL per CFUe produced per h of growth flux
    KGM_LYSIS=3.90e-10,
    KGM_ABS=0.045,          # lesion -> serum, scaled by perfusion
    KGM_DEG=0.030,
    KGM_EL=0.030,           # serum GM elimination, 1/h (t1/2 ~23 h)
    GM_VD=1.0,
    KBDG_ON=2.0e-8, KBDG_OFF=0.028,
    KPCR_ON=1.0e-6, KPCR_OFF=0.075,
    # ===================== organ / toxicity =================================
    ALT_BASE=25.0, KALT_OFF=0.011,
    ALT_VRC_SLOPE=1.55,     # U/L per h per mg/L above threshold
    ALT_VRC_THR=4.0,        # mg/L; proxy for the trough >4-5.5 hepatotoxicity signal
    SCR_BASE=0.90, KSCR_OFF=0.020,
    SCR_AMB_SLOPE=3.0e-4,   # mg/dL per h per mg/L above threshold
    SCR_AMB_THR=18.0,       # L-AmB plasma conc above which tubular injury runs
    KMG_ON=1.4e-4, KMG_OFF=0.020,
    QTC_VRC=3.4, QTC_POS=5.0, KQTC_OFF=0.30,
    NEURO_ON=0.062, NEURO_OFF=0.10, NEURO_THR=4.0,
    PFR_BASE=400.0, KPFR_REC=0.010, KPFR_LOSS=0.030,
    TEMP_BASE=36.8, KTEMP=0.10, TEMP_GAIN=0.030,
    # ===================== mortality hazard =================================
    HAZ0=2.4e-5,
    HAZ_BURDEN=2.20e-4,     # per log10 CFUe above 4 per h
    HAZ_PERF=1.00e-3,       # per (1-PERF) per h
    HAZ_HYPOX=1.10e-3,      # per unit (1 - PFR/400) per h
    HAZ_BRAIN=1.20e-4,      # per log10 cerebral CFUe per h
    HAZ_NEUTRO=2.50e-4,     # per h while ANC < 0.5
    HAZ_ALT=1.5e-6,
    HAZ_SCR=2.6e-4,
    # ===================== host phenotype switches ==========================
    CYP2C19=1.0,            # multiplier on VRC_Vmax
    STEROID_ON=0.0,
    RESIST_FRAC=0.0,        # fraction of the inoculum already TR34/L98H
    T_ANC_REC=240.0,        # h at which marrow recovery starts
    INOCULUM=1.0e3,         # conidia deposited in the alveoli
    ANC0=0.05,              # ANC at inoculation
    ADJ_REDUCE_IS=0.0,      # 1 = immunosuppression reduced at diagnosis
)

# CYP2C19 diplotype -> Vmax multiplier (Weiss 2009; Owusu Obeng 2014)
CYP2C19 = {
    "UM (*17/*17)": 1.85,
    "RM (*1/*17)": 1.38,
    "NM (*1/*1)": 1.00,
    "IM (*1/*2)": 0.62,
    "PM (*2/*2)": 0.30,
}


# ----------------------------------------------------------------------------
# 3. DOSING
# ----------------------------------------------------------------------------
class Regimen:
    """A list of bolus doses and of zero-order infusions."""

    def __init__(self):
        self.boluses = []      # (time_h, state_name, amount)
        self.infusions = []    # (start_h, end_h, state_name, rate_amount_per_h)

    def bolus(self, t, cmt, amt):
        self.boluses.append((float(t), cmt, float(amt)))
        return self

    def infuse(self, t, cmt, amt, dur):
        self.infusions.append((float(t), float(t + dur), cmt, amt / float(dur)))
        return self

    def multi_bolus(self, t0, cmt, amt, ii, n):
        for k in range(int(n)):
            self.bolus(t0 + k * ii, cmt, amt)
        return self

    def multi_infuse(self, t0, cmt, amt, dur, ii, n):
        for k in range(int(n)):
            self.infuse(t0 + k * ii, cmt, amt, dur)
        return self

    def event_times(self):
        ts = {b[0] for b in self.boluses}
        for s, e, _, _ in self.infusions:
            ts.add(s)
            ts.add(e)
        return sorted(ts)

    def infusion_rate(self, t, cmt):
        r = 0.0
        for s, e, c, rate in self.infusions:
            if c == cmt and s <= t < e:
                r += rate
        return r


def hill(c, ec50, n):
    if c <= 0.0:
        return 0.0
    if c > 1e6 * ec50:          # avoid c**n overflowing for burdens ~1e10
        return 1.0
    cn = c ** n
    return cn / (ec50 ** n + cn)


# ----------------------------------------------------------------------------
# 4. RIGHT-HAND SIDE
# ----------------------------------------------------------------------------
def derivs(t, y, p, reg):
    y = np.maximum(y, 0.0)
    d = np.zeros(NS)
    g = lambda n: y[IDX[n]]

    # ---------------- concentrations ---------------------------------------
    Cvrc = g("VRC_c") / p["VRC_Vc"]
    Cisa = g("ISA_c") / p["ISA_Vc"]
    Cpos = g("POS_c") / p["POS_Vc"]
    Camb = g("AMB_c") / p["AMB_Vc"]
    Cech = g("ECH_c") / p["ECH_Vc"]
    Ctac = g("TAC_c") / p["TAC_Vc"] * 1000.0          # ng/mL
    Cgcs = g("GCS_c") / p["GCS_Vc"]
    Cgcsf = g("GCSF_c") / p["GCSF_Vc"] * 1000.0       # ng/mL

    Evrc = g("VRC_elf")
    Eisa = g("ISA_elf")
    Epos = g("POS_elf")
    Eamb = g("AMB_elf")
    Eech = g("ECH_elf")

    # ---------------- lesion perfusion gates drug delivery -----------------
    PERF = min(max(g("PERF"), 1e-4), 1.0)
    fdel = p["DELIV_FLOOR"] + (1.0 - p["DELIV_FLOOR"]) * PERF ** p["PERF_GAMMA"]
    Svrc, Sisa, Spos, Samb, Sech = (Evrc * fdel, Eisa * fdel, Epos * fdel,
                                    Eamb * fdel, Eech * fdel)

    # ---------------- PK: voriconazole (Michaelis-Menten) ------------------
    infl = 1.0 - p["VRC_infl_Emax"] * hill(g("CRP"), p["VRC_infl_EC50"], 1.0)
    vmax = p["VRC_Vmax"] * p["CYP2C19"] * infl
    cl_vrc = (vmax * Cvrc / (p["VRC_Km"] + Cvrc)               # CYP2C19, mg/h
              + p["VRC_CLlin"] * Cvrc)                          # CYP3A4/FMO3
    q_vrc = p["VRC_Q"] * (Cvrc - g("VRC_p") / p["VRC_Vp"])
    d[IDX["VRC_gut"]] = -p["VRC_ka"] * g("VRC_gut")
    d[IDX["VRC_c"]] = (p["VRC_F"] * p["VRC_ka"] * g("VRC_gut")
                       + reg.infusion_rate(t, "VRC_c") - cl_vrc - q_vrc)
    d[IDX["VRC_p"]] = q_vrc
    d[IDX["VRC_elf"]] = p["VRC_keq_elf"] * (p["VRC_Kp_elf"] * Cvrc - Evrc)
    d[IDX["VRC_brain"]] = p["VRC_keq_brain"] * (p["VRC_Kp_brain"] * Cvrc
                                                - g("VRC_brain"))

    # ---------------- PK: isavuconazole -------------------------------------
    q_isa = p["ISA_Q"] * (Cisa - g("ISA_p") / p["ISA_Vp"])
    d[IDX["ISA_gut"]] = -p["ISA_ka"] * g("ISA_gut")
    d[IDX["ISA_c"]] = (p["ISA_F"] * p["ISA_ka"] * g("ISA_gut")
                       + reg.infusion_rate(t, "ISA_c")
                       - p["ISA_CL"] * Cisa - q_isa)
    d[IDX["ISA_p"]] = q_isa
    d[IDX["ISA_elf"]] = p["ISA_keq_elf"] * (p["ISA_Kp_elf"] * Cisa - Eisa)

    # ---------------- PK: posaconazole --------------------------------------
    q_pos = p["POS_Q"] * (Cpos - g("POS_p") / p["POS_Vp"])
    d[IDX["POS_gut"]] = -p["POS_ka"] * g("POS_gut")
    d[IDX["POS_c"]] = (p["POS_F"] * p["POS_ka"] * g("POS_gut")
                       + reg.infusion_rate(t, "POS_c")
                       - p["POS_CL"] * Cpos - q_pos)
    d[IDX["POS_p"]] = q_pos
    d[IDX["POS_elf"]] = p["POS_keq_elf"] * (p["POS_Kp_elf"] * Cpos - Epos)

    # ---------------- PK: liposomal amphotericin B --------------------------
    q_amb = p["AMB_Q"] * (Camb - g("AMB_p") / p["AMB_Vp"])
    d[IDX["AMB_c"]] = (reg.infusion_rate(t, "AMB_c")
                       - p["AMB_CL"] * Camb - q_amb)
    d[IDX["AMB_p"]] = q_amb
    d[IDX["AMB_elf"]] = p["AMB_keq_elf"] * (p["AMB_Kp_elf"] * Camb - Eamb)
    d[IDX["AMB_brain"]] = p["AMB_keq_brain"] * (p["AMB_Kp_brain"] * Camb
                                                - g("AMB_brain"))

    # ---------------- PK: echinocandin --------------------------------------
    q_ech = p["ECH_Q"] * (Cech - g("ECH_p") / p["ECH_Vp"])
    d[IDX["ECH_c"]] = (reg.infusion_rate(t, "ECH_c")
                       - p["ECH_CL"] * Cech - q_ech)
    d[IDX["ECH_p"]] = q_ech
    d[IDX["ECH_elf"]] = p["ECH_keq_elf"] * (p["ECH_Kp_elf"] * Cech - Eech)

    # ---------------- PK: tacrolimus, with azole CYP3A4 inhibition ----------
    # Competitive CYP3A4 inhibition: each azole contributes C/Ki independently.
    inh = 1.0 / (1.0 + Cvrc / p["TAC_KI_VRC"] + Cisa / p["TAC_KI_ISA"]
                 + Cpos / p["TAC_KI_POS"])
    d[IDX["TAC_c"]] = (reg.infusion_rate(t, "TAC_c")
                       - p["TAC_CL"] * inh * (g("TAC_c") / p["TAC_Vc"]))

    d[IDX["GCS_c"]] = reg.infusion_rate(t, "GCS_c") - p["GCS_CL"] * Cgcs
    d[IDX["GCSF_c"]] = reg.infusion_rate(t, "GCSF_c") - p["GCSF_CL"] * Cgcsf / 1000.0

    # ---------------- host defence competence -------------------------------
    fster = 1.0 / (1.0 + Cgcs / p["GCS_IC50"])
    shield = 1.0 - p["BIO_SHIELD"] * hill(g("BIO"), p["KBIO"], 1.0)
    Neff = hill(g("NLES"), p["KN50"], 1.0) * fster * shield
    surveil = min(g("ANC") / p["ANC_SS"], 1.5) * fster * shield
    kill_host = p["KHOST_MAX"] * Neff + p["KHOST_BASE"] * surveil

    # ---------------- drug effects at the lesion ----------------------------
    # Azoles act on GROWTH; each has its own EC50 scaled by the isolate MIC.
    def az_inhib(mic_w, mic_r):
        iw = (p["IMAX_AZOLE"] * hill(Svrc, p["VRC_EC50_MIC"] * mic_w[0], p["HILL_AZ"])
              + p["IMAX_AZOLE"] * hill(Sisa, p["ISA_EC50_MIC"] * mic_w[1], p["HILL_AZ"])
              + p["IMAX_AZOLE"] * hill(Spos, p["POS_EC50_MIC"] * mic_w[2], p["HILL_AZ"]))
        ir = (p["IMAX_AZOLE"] * hill(Svrc, p["VRC_EC50_MIC"] * mic_r[0], p["HILL_AZ"])
              + p["IMAX_AZOLE"] * hill(Sisa, p["ISA_EC50_MIC"] * mic_r[1], p["HILL_AZ"])
              + p["IMAX_AZOLE"] * hill(Spos, p["POS_EC50_MIC"] * mic_r[2], p["HILL_AZ"]))
        return min(iw, p["IMAX_AZOLE"]), min(ir, p["IMAX_AZOLE"])

    Iw, Ir = az_inhib((p["MIC_VRC"], p["MIC_ISA"], p["MIC_POS"]),
                      (p["MIC_VRC_R"], p["MIC_ISA_R"], p["MIC_POS_R"]))

    # Echinocandin: apical growth inhibition, attenuated above CPAR
    par = 1.0 / (1.0 + (Sech / p["ECH_CPAR"]) ** 3)
    Iech = p["ECH_IMAX"] * hill(Sech, p["ECH_EC50"], p["HILL_ECH"]) * par
    Iw = min(1.0 - (1.0 - Iw) * (1.0 - Iech), 0.995)
    Ir = min(1.0 - (1.0 - Ir) * (1.0 - Iech), 0.995)

    # Amphotericin B: killing
    kill_amb_w = p["AMB_EMAX"] * hill(Samb, p["AMB_EC50_MIC"] * p["MIC_AMB"],
                                      p["HILL_AMB"])
    kill_amb_r = p["AMB_EMAX"] * hill(Samb, p["AMB_EC50_MIC"] * p["MIC_AMB_R"],
                                      p["HILL_AMB"])

    # ---------------- fungal dynamics ---------------------------------------
    HY = g("HYW") + g("HYR")
    logistic = max(0.0, 1.0 - HY / p["BMAX"])

    kgerm = p["KGERM"] * (1.0 - 0.55 * hill(g("MAC") * fster, 0.7, 1.0))
    conid_kill = p["KMAC_CON"] * g("MAC") * fster
    d[IDX["CON"]] = -kgerm * g("CON") - conid_kill * g("CON")
    d[IDX["GERM"]] = (kgerm * g("CON") - p["KMAT"] * g("GERM")
                      - 0.6 * kill_host * g("GERM") - conid_kill * 0.35 * g("GERM"))

    # Extinction gate: a continuous state can carry 1e-30 "CFUe" and regrow
    # from it. Nothing below one organism can divide.
    gate_w = hill(g("HYW"), 1.0, 4.0)
    gate_r = hill(g("HYR"), 1.0, 4.0)
    grow_w = p["KGROW"] * (1.0 - Iw) * logistic * gate_w
    grow_r = p["KGROW"] * (1.0 - Ir) * logistic * gate_r
    mut_flux = p["KMUT"] * grow_w * g("HYW")
    lyse_w = (kill_host + kill_amb_w) * g("HYW")
    lyse_r = (kill_host + kill_amb_r) * g("HYR")

    d[IDX["HYW"]] = (p["KMAT"] * g("GERM") * (1.0 - p["RESIST_FRAC"])
                     + grow_w * g("HYW") - lyse_w - mut_flux)
    d[IDX["HYR"]] = (p["KMAT"] * g("GERM") * p["RESIST_FRAC"]
                     + grow_r * g("HYR") - lyse_r + mut_flux)

    # cerebral dissemination: seeding is proportional to angioinvasion
    seed = 2.6e-8 * g("ANG") * HY
    kill_brain = (kill_host * 0.25
                  + p["AMB_EMAX"] * hill(g("AMB_brain"),
                                         p["AMB_EC50_MIC"] * p["MIC_AMB"],
                                         p["HILL_AMB"]))
    grow_brain = p["KGROW_BRAIN"] * (
        1.0 - min(p["IMAX_AZOLE"] * hill(g("VRC_brain"),
                                         p["VRC_EC50_MIC"] * p["MIC_VRC"],
                                         p["HILL_AZ"])
                  + p["IMAX_AZOLE"] * hill(Cisa * p["ISA_Kp_brain"],
                                           p["ISA_EC50_MIC"] * p["MIC_ISA"],
                                           p["HILL_AZ"]),
                  p["IMAX_AZOLE"]))
    lg_brain = max(0.0, 1.0 - g("BRB") / p["BMAX_BRAIN"])
    d[IDX["BRB"]] = (seed * lg_brain
                     + grow_brain * lg_brain * hill(g("BRB"), 1.0, 4.0) * g("BRB")
                     - kill_brain * g("BRB"))

    d[IDX["BIO"]] = p["KBIO_ON"] * HY - p["KBIO_OFF"] * g("BIO")

    # ---------------- lesion, angioinvasion, perfusion ----------------------
    d[IDX["VLES"]] = (p["KVLES"] * HY
                      - (p["KRESOLVE_BASE"]
                         + p["KRESOLVE"] * hill(g("NLES"), 0.5, 1.0)) * g("VLES"))
    d[IDX["ANG"]] = p["KANG"] * HY - p["KANG_OFF"] * g("ANG")
    heal = hill(g("NLES"), 0.5, 1.0) * (1.0 - hill(HY, 1.0e7, 1.0))
    d[IDX["PERF"]] = (-p["KPERF_LOSS"] * hill(g("ANG"), p["KANG50"], 1.0) * PERF
                      + p["KPERF_REC"] * (1.0 - PERF) * (0.25 + heal))
    # NEC is a SUBSET of VLES -- production works on the still-perfused rim,
    # which is what keeps the infarcted volume below the lesion volume.
    d[IDX["NEC"]] = (p["KNEC"] * (1.0 - PERF) * max(g("VLES") - g("NEC"), 0.0)
                     - p["KNEC_OFF"] * g("NEC"))

    # ---------------- haematopoiesis ----------------------------------------
    marrow_on = 1.0 if t >= p["T_ANC_REC"] else 0.0
    egcsf = p["GCSF_EMAX"] * hill(Cgcsf, p["GCSF_EC50"], 1.0)
    d[IDX["MAR"]] = (p["KPROL"] * marrow_on * (1.0 + egcsf) * p["ANC_SS"]
                     - p["KTR"] * g("MAR"))
    marg = p["KMARG"] * hill(HY, p["KMARG_EC50"], 1.0) * g("ANC")
    d[IDX["ANC"]] = p["KTR"] * g("MAR") - p["KOUT_ANC"] * g("ANC") - marg
    d[IDX["NLES"]] = (p["KNLES_IN"] * marg * 100.0 * fster
                      - p["KNLES_OUT"] * g("NLES"))
    d[IDX["MAC"]] = p["KMAC_TURN"] * (p["MAC_SS"] - g("MAC"))

    # ---------------- cytokines ---------------------------------------------
    d[IDX["IFNG"]] = p["KIFNG_ON"] * hill(HY, 1.0e6, 1.0) * fster - p["KIFNG_OFF"] * g("IFNG")
    d[IDX["IL6"]] = (p["KIL6_MAX"] * hill(HY, p["KIL6_EC50"], 1.0)
                     * (0.3 + 0.7 * fster) - p["KIL6_OFF"] * g("IL6"))
    d[IDX["CRP"]] = p["KCRP_ON"] * g("IL6") - p["KCRP_OFF"] * g("CRP")
    d[IDX["IL10"]] = p["KIL10_ON"] * hill(g("IL6"), 250.0, 1.0) - p["KIL10_OFF"] * g("IL10")

    # ---------------- biomarkers --------------------------------------------
    growth_flux = grow_w * g("HYW") + grow_r * g("HYR")
    lysis_flux = lyse_w + lyse_r
    d[IDX["GMles"]] = (p["KGM_GROWTH"] * growth_flux
                       + p["KGM_LYSIS"] * lysis_flux
                       - p["KGM_ABS"] * PERF * g("GMles")
                       - p["KGM_DEG"] * g("GMles"))
    d[IDX["GMser"]] = (p["KGM_ABS"] * PERF * g("GMles") / p["GM_VD"]
                       - p["KGM_EL"] * g("GMser"))
    d[IDX["BDG"]] = p["KBDG_ON"] * (growth_flux + lysis_flux) - p["KBDG_OFF"] * g("BDG")
    d[IDX["PCRs"]] = p["KPCR_ON"] * (0.25 * growth_flux + lysis_flux) - p["KPCR_OFF"] * g("PCRs")

    # ---------------- organ toxicity ----------------------------------------
    vrc_ex = max(0.0, Cvrc - p["ALT_VRC_THR"])
    d[IDX["ALT"]] = (p["ALT_VRC_SLOPE"] * vrc_ex
                     + 0.9 * max(0.0, Cisa - 8.0)
                     - p["KALT_OFF"] * (g("ALT") - p["ALT_BASE"]))
    amb_ex = max(0.0, Camb - p["SCR_AMB_THR"])
    d[IDX["SCR"]] = (p["SCR_AMB_SLOPE"] * amb_ex
                     - p["KSCR_OFF"] * (g("SCR") - p["SCR_BASE"]))
    d[IDX["KMG"]] = p["KMG_ON"] * Camb - p["KMG_OFF"] * g("KMG")
    d[IDX["QTC"]] = (p["QTC_VRC"] * Cvrc + p["QTC_POS"] * Cpos
                     - p["KQTC_OFF"] * g("QTC"))
    d[IDX["NEUROE"]] = (p["NEURO_ON"] * max(0.0, Cvrc - p["NEURO_THR"])
                        - p["NEURO_OFF"] * g("NEUROE"))
    # shunt fraction scales with consolidated + infarcted volume; written on
    # PFR/PFR_BASE so oxygenation degrades toward a floor instead of to zero
    d[IDX["PFR"]] = (p["KPFR_REC"] * (p["PFR_BASE"] - g("PFR"))
                     - p["KPFR_LOSS"] * (g("VLES") + 0.8 * g("NEC"))
                     * g("PFR") / p["PFR_BASE"])

    d[IDX["TEMP"]] = (p["KTEMP"] * (p["TEMP_BASE"] - g("TEMP"))
                      + p["TEMP_GAIN"] * hill(g("IL6"), 90.0, 1.0) * fster * 10.0)

    # ---------------- mortality hazard --------------------------------------
    lb = math.log10(max(HY, 1.0))
    lbb = math.log10(max(g("BRB"), 1.0))
    haz = (p["HAZ0"]
           + p["HAZ_BURDEN"] * max(0.0, lb - 4.0)
           + p["HAZ_PERF"] * (1.0 - PERF)
           + p["HAZ_HYPOX"] * max(0.0, 1.0 - g("PFR") / p["PFR_BASE"])
           + p["HAZ_BRAIN"] * max(0.0, lbb - 2.0)
           + (p["HAZ_NEUTRO"] if g("ANC") < 0.5 else 0.0)
           + p["HAZ_ALT"] * max(0.0, g("ALT") - 120.0)
           + p["HAZ_SCR"] * max(0.0, g("SCR") - 1.5))
    d[IDX["HAZ"]] = haz
    return d


# ----------------------------------------------------------------------------
# 5. INITIAL CONDITIONS & SIMULATION DRIVER
# ----------------------------------------------------------------------------
def y0_from(p):
    y = np.zeros(NS)
    y[IDX["CON"]] = p["INOCULUM"]
    y[IDX["PERF"]] = 1.0
    y[IDX["MAC"]] = p["MAC_SS"]
    y[IDX["ANC"]] = p["ANC0"]
    y[IDX["MAR"]] = p["ANC0"] * p["KOUT_ANC"] / p["KTR"]
    y[IDX["ALT"]] = p["ALT_BASE"]
    y[IDX["SCR"]] = p["SCR_BASE"]
    y[IDX["PFR"]] = p["PFR_BASE"]
    y[IDX["TEMP"]] = p["TEMP_BASE"]
    return y


def simulate(par_over=None, reg=None, tmax=2016.0, dt=2.0):
    """Integrate the 53-state system. tmax default = 12 weeks."""
    p = OrderedDict(P)
    if par_over:
        p.update(par_over)
    reg = reg or Regimen()

    y = y0_from(p)
    tgrid = np.arange(0.0, tmax + dt / 2, dt)
    out = np.zeros((len(tgrid), NS))
    out[0] = y

    # break integration at every dose / infusion boundary
    breaks = sorted(set([0.0, tmax] + [e for e in reg.event_times() if 0 < e < tmax]))
    bol = {}
    for tt, cmt, amt in reg.boluses:
        bol.setdefault(round(tt, 6), []).append((cmt, amt))

    ti = 0
    for a, b in zip(breaks[:-1], breaks[1:]):
        for cmt, amt in bol.get(round(a, 6), []):
            y[IDX[cmt]] += amt
        seg = tgrid[(tgrid > a) & (tgrid <= b)]
        # ALWAYS evaluate at the right-hand break itself, appending it if it is
        # not already an output point, so the state handed to the next segment
        # is the state at `b` and no sub-grid interval is dropped.
        if seg.size and abs(seg[-1] - b) < 1e-9:
            teval = seg
        else:
            teval = np.append(seg, b)
        sol = solve_ivp(derivs, (a, b), y, args=(p, reg), method="LSODA",
                        t_eval=teval, rtol=1e-8, atol=1e-10, max_step=0.5)
        if not sol.success:
            raise RuntimeError(f"integration failed on [{a},{b}]: {sol.message}")
        ys = np.asarray(sol.y)
        if seg.size:
            out[ti + 1: ti + 1 + seg.size] = ys[:, :seg.size].T
            ti += seg.size
        y = ys[:, -1].copy()
    for cmt, amt in bol.get(round(tmax, 6), []):
        y[IDX[cmt]] += amt
    return tgrid, np.maximum(out, 0.0)


# ----------------------------------------------------------------------------
# 6. DERIVED OUTPUTS
# ----------------------------------------------------------------------------
def burden_log(out):
    hy = out[:, IDX["HYW"]] + out[:, IDX["HYR"]]
    return np.log10(np.maximum(hy, 1.0))


def survival(out):
    return np.exp(-out[:, IDX["HAZ"]])


def auc(t, c):
    return float(np.trapezoid(c, t)) if hasattr(np, "trapezoid") else float(np.trapz(c, t))


def summarise(tag, t, out, p=None):
    hy = out[:, IDX["HYW"]] + out[:, IDX["HYR"]]
    lb = np.log10(np.maximum(hy, 1.0))
    s = survival(out)
    i42 = int(np.argmin(np.abs(t - 42 * 24)))
    i84 = int(np.argmin(np.abs(t - 84 * 24)))
    return OrderedDict(
        scenario=tag,
        log10_burden_peak=round(float(lb.max()), 2),
        t_peak_burden_d=round(float(t[int(np.argmax(lb))] / 24), 1),
        log10_burden_d14=round(float(lb[int(np.argmin(np.abs(t - 336)))]), 2),
        log10_burden_d42=round(float(lb[i42]), 2),
        log10_burden_d84=round(float(lb[i84]), 2),
        GM_peak=round(float(out[:, IDX["GMser"]].max()), 2),
        t_GM_peak_d=round(float(t[int(np.argmax(out[:, IDX["GMser"]]))] / 24), 1),
        GM_d14=round(float(out[int(np.argmin(np.abs(t - 336))), IDX["GMser"]]), 2),
        PERF_min=round(float(out[:, IDX["PERF"]].min()), 3),
        lesion_peak_mL=round(float(out[:, IDX["VLES"]].max()), 1),
        ALT_max=round(float(out[:, IDX["ALT"]].max()), 0),
        SCr_max=round(float(out[:, IDX["SCR"]].max()), 2),
        mortality_d42_pct=round(100 * (1 - float(s[i42])), 1),
        mortality_d84_pct=round(100 * (1 - float(s[i84])), 1),
    )


# ----------------------------------------------------------------------------
# 7. REGIMEN BUILDERS
# ----------------------------------------------------------------------------
WT = 70.0


def reg_vrc_iv(start_d=3.0, days=42, load_mg_kg=6.0, maint_mg_kg=4.0):
    r = Regimen()
    t0 = start_d * 24
    r.infuse(t0, "VRC_c", load_mg_kg * WT, 2.0)
    r.infuse(t0 + 12, "VRC_c", load_mg_kg * WT, 2.0)
    r.multi_infuse(t0 + 24, "VRC_c", maint_mg_kg * WT, 2.0, 12.0, days * 2)
    return r


def reg_isa(start_d=3.0, days=42):
    r = Regimen()
    t0 = start_d * 24
    r.multi_infuse(t0, "ISA_c", 200.0, 1.0, 8.0, 6)
    r.multi_infuse(t0 + 48, "ISA_c", 200.0, 1.0, 24.0, days)
    return r


def reg_amb(start_d=3.0, days=42, mgkg=3.0):
    r = Regimen()
    r.multi_infuse(start_d * 24, "AMB_c", mgkg * WT, 2.0, 24.0, days)
    return r


def reg_ech(start_d=3.0, days=42):
    r = Regimen()
    t0 = start_d * 24
    r.infuse(t0, "ECH_c", 200.0, 1.5)
    r.multi_infuse(t0 + 24, "ECH_c", 100.0, 1.5, 24.0, days)
    return r


def reg_pos(start_d=0.0, days=90):
    r = Regimen()
    t0 = start_d * 24
    r.multi_bolus(t0, "POS_gut", 300.0, 12.0, 2)
    r.multi_bolus(t0 + 24, "POS_gut", 300.0, 24.0, days)
    return r


def add_steroid(r, mg_per_day=60.0, days=60, start_d=0.0):
    for k in range(days):
        r.infuse(start_d * 24 + k * 24, "GCS_c", mg_per_day, 1.0)
    return r


def add_gcsf(r, start_d, days=10, mcg=300.0):
    for k in range(days):
        r.infuse((start_d + k) * 24, "GCSF_c", mcg / 1000.0, 1.0)
    return r


def add_tac(r, days=90, mg_bid=3.0):
    for k in range(days * 2):
        r.bolus(k * 12.0, "TAC_c", mg_bid)
    return r


def merge(*regs):
    out = Regimen()
    for r in regs:
        out.boluses += r.boluses
        out.infusions += r.infusions
    return out


# ----------------------------------------------------------------------------
# 8. ANALYSES
# ----------------------------------------------------------------------------
RESULTS = OrderedDict()
LINES = []


def say(s=""):
    LINES.append(s)
    print(s)


def hdr(s):
    say("")
    say("=" * 78)
    say(s)
    say("=" * 78)


def main():
    # ---------------------------------------------------------------- PK QC
    hdr("A. PHARMACOKINETIC QUALIFICATION (drug alone, no infection)")
    say(f"{'drug / genotype':38s} {'Cmin,ss':>9s} {'Cmax,ss':>9s} {'AUC24,ss':>10s}")
    say("-" * 78)
    pk_qc = OrderedDict()

    for gt, mult in CYP2C19.items():
        t, o = simulate({"CYP2C19": mult, "INOCULUM": 0.0},
                        reg_vrc_iv(start_d=0.0, days=14), tmax=14 * 24, dt=0.25)
        c = o[:, IDX["VRC_c"]] / P["VRC_Vc"]
        m = (t >= 12 * 24) & (t <= 13 * 24)
        cmin, cmax = float(c[m].min()), float(c[m].max())
        a24 = auc(t[m], c[m])
        pk_qc[f"VRC 4 mg/kg q12h IV, {gt}"] = dict(Cmin=round(cmin, 2),
                                                   Cmax=round(cmax, 2),
                                                   AUC24=round(a24, 1))
        say(f"{'VRC 4 mg/kg q12h IV, ' + gt:38s} {cmin:9.2f} {cmax:9.2f} {a24:10.1f}")

    for tag, reg, cmt, vol in (
        ("ISA 200 mg q24h IV", reg_isa(0.0, 14), "ISA_c", P["ISA_Vc"]),
        ("POS DR-tab 300 mg q24h PO", reg_pos(0.0, 14), "POS_c", P["POS_Vc"]),
        ("L-AmB 3 mg/kg q24h", reg_amb(0.0, 14, 3.0), "AMB_c", P["AMB_Vc"]),
        ("L-AmB 10 mg/kg q24h", reg_amb(0.0, 14, 10.0), "AMB_c", P["AMB_Vc"]),
        ("Anidulafungin 100 mg q24h", reg_ech(0.0, 14), "ECH_c", P["ECH_Vc"]),
    ):
        t, o = simulate({"INOCULUM": 0.0}, reg, tmax=14 * 24, dt=0.25)
        c = o[:, IDX[cmt]] / vol
        m = (t >= 12 * 24) & (t <= 13 * 24)
        cmin, cmax, a24 = float(c[m].min()), float(c[m].max()), auc(t[m], c[m])
        pk_qc[tag] = dict(Cmin=round(cmin, 2), Cmax=round(cmax, 2), AUC24=round(a24, 1))
        say(f"{tag:38s} {cmin:9.2f} {cmax:9.2f} {a24:10.1f}")

    RESULTS["A_pk_qualification"] = pk_qc

    say("")
    say("NOTE: the VRC row set is the model's nonlinearity in one column --")
    r_um = pk_qc["VRC 4 mg/kg q12h IV, UM (*17/*17)"]["AUC24"]
    r_pm = pk_qc["VRC 4 mg/kg q12h IV, PM (*2/*2)"]["AUC24"]
    say(f"  AUC24 spans {r_um:.1f} (UM) to {r_pm:.1f} (PM) mg*h/L on an IDENTICAL")
    say(f"  mg/kg dose: a {r_pm / r_um:.1f}-fold genotype-driven exposure range.")
    RESULTS["A_cyp2c19_auc_fold"] = round(r_pm / r_um, 2)

    # -------------------------------------------------- dose-proportionality
    hdr("B. VORICONAZOLE IS SUPRAPROPORTIONAL: 1.5x DOSE IS NOT 1.5x EXPOSURE")
    say(f"{'dose (mg/kg q12h)':>20s} {'AUC24':>10s} {'ratio vs 3':>12s} {'dose ratio':>12s}")
    say("-" * 78)
    prop = []
    base = None
    for mgkg in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        t, o = simulate({"INOCULUM": 0.0},
                        reg_vrc_iv(0.0, 14, mgkg, mgkg), tmax=14 * 24, dt=0.25)
        c = o[:, IDX["VRC_c"]] / P["VRC_Vc"]
        m = (t >= 12 * 24) & (t <= 13 * 24)
        a = auc(t[m], c[m])
        if mgkg == 3.0:
            base = a
        prop.append((mgkg, round(a, 1)))
    for mgkg, a in prop:
        say(f"{mgkg:20.1f} {a:10.1f} {a / base:12.2f} {mgkg / 3.0:12.2f}")
    RESULTS["B_dose_proportionality"] = [dict(mg_kg=m_, AUC24=a_) for m_, a_ in prop]
    r8 = [a for m_, a in prop if m_ == 8.0][0]
    say("")
    say(f"  Doubling 4 -> 8 mg/kg multiplies AUC24 by "
        f"{r8 / [a for m_, a in prop if m_ == 4.0][0]:.2f}, not 2.00.")
    RESULTS["B_fold_4to8"] = round(r8 / [a for m_, a in prop if m_ == 4.0][0], 2)

    # ------------------------------------------------- the structural result
    hdr("C. THE CENTRAL STRUCTURAL RESULT: A FUNGISTATIC DRUG CANNOT CLEAR "
        "A HOST\n   WITH NO NEUTROPHILS")
    say("dB/dt = kgrow*(1 - I_azole)*B - (k_host*N_eff + k_cidal)*B")
    say(f"kgrow = {P['KGROW']:.3f}/h, Imax_azole = {P['IMAX_AZOLE']:.3f}")
    floor = P["KGROW"] * (1 - P["IMAX_AZOLE"])
    say(f"=> with N_eff = 0 and a perfectly dosed azole, the floor of the growth")
    say(f"   term is {floor:.5f}/h = {floor * 24 / math.log(10):.3f} log10/day.")
    say(f"   Doubling time at maximal azole effect: {math.log(2) / floor:.0f} h "
        f"({math.log(2) / floor / 24:.1f} days).")
    say(f"   The minimum k_host needed to make dB/dt < 0 is {floor:.5f}/h,")
    say(f"   i.e. N_eff > {floor / P['KHOST_MAX']:.4f} "
        f"({100 * floor / P['KHOST_MAX']:.1f}% of maximal neutrophil function).")
    RESULTS["C_azole_growth_floor_per_h"] = round(floor, 6)
    RESULTS["C_azole_growth_floor_log10_per_day"] = round(floor * 24 / math.log(10), 4)
    RESULTS["C_min_Neff_for_clearance"] = round(floor / P["KHOST_MAX"], 4)

    # ------------------------------------------------------------ scenarios
    hdr("D. TREATMENT SCENARIOS (inoculum on day 0, therapy from day 3, "
        "12-week horizon)")
    scen = OrderedDict()

    def run(tag, over=None, reg=None):
        o_ = dict(over or {})
        t, y = simulate(o_, reg, tmax=84 * 24, dt=4.0)
        scen[tag] = summarise(tag, t, y)
        return t, y

    run("S01 untreated, persistent neutropenia", {"T_ANC_REC": 1e9})
    run("S02 untreated, ANC recovery d10", {"T_ANC_REC": 240.0})
    run("S03 VRC NM, persistent neutropenia", {"T_ANC_REC": 1e9}, reg_vrc_iv())
    run("S04 VRC NM, ANC recovery d10", {"T_ANC_REC": 240.0}, reg_vrc_iv())
    run("S05 VRC UM (*17/*17), ANC rec d10",
        {"CYP2C19": CYP2C19["UM (*17/*17)"], "T_ANC_REC": 240.0}, reg_vrc_iv())
    run("S06 VRC PM (*2/*2), ANC rec d10",
        {"CYP2C19": CYP2C19["PM (*2/*2)"], "T_ANC_REC": 240.0}, reg_vrc_iv())
    run("S07 VRC UM + TDM dose escalation",
        {"CYP2C19": CYP2C19["UM (*17/*17)"], "T_ANC_REC": 240.0},
        reg_vrc_iv(maint_mg_kg=4.0 * 1.9))
    run("S08 isavuconazole, ANC rec d10", {"T_ANC_REC": 240.0}, reg_isa())
    run("S09 L-AmB 3 mg/kg, ANC rec d10", {"T_ANC_REC": 240.0}, reg_amb(mgkg=3.0))
    run("S10 L-AmB 10 mg/kg, ANC rec d10", {"T_ANC_REC": 240.0}, reg_amb(mgkg=10.0))
    run("S11 anidulafungin mono, ANC rec d10", {"T_ANC_REC": 240.0}, reg_ech())
    run("S12 VRC + anidulafungin, ANC rec d10", {"T_ANC_REC": 240.0},
        merge(reg_vrc_iv(), reg_ech()))
    run("S13 VRC vs TR34/L98H (MIC 8)",
        {"T_ANC_REC": 240.0, "RESIST_FRAC": 1.0}, reg_vrc_iv())
    run("S14 TR34/L98H salvaged with L-AmB d10",
        {"T_ANC_REC": 240.0, "RESIST_FRAC": 1.0},
        merge(reg_vrc_iv(days=7), reg_amb(start_d=10.0, days=35)))
    run("S15 steroid host (no neutropenia), VRC",
        {"T_ANC_REC": 0.0, "ANC0": 4.0},
        add_steroid(reg_vrc_iv(), 60.0, 60))
    run("S16 steroid host, VRC + steroid taper d14",
        {"T_ANC_REC": 0.0, "ANC0": 4.0},
        add_steroid(reg_vrc_iv(), 60.0, 14))
    run("S17 VRC + G-CSF from d3", {"T_ANC_REC": 240.0},
        add_gcsf(reg_vrc_iv(), 3.0, 12))
    run("S18 late start d10, VRC", {"T_ANC_REC": 240.0}, reg_vrc_iv(start_d=10.0))
    run("S19 CNS disease, VRC", {"T_ANC_REC": 240.0}, reg_vrc_iv())
    run("S20 CNS disease, L-AmB", {"T_ANC_REC": 240.0}, reg_amb(mgkg=5.0))
    run("S21 transplant, tacrolimus + VRC", {"T_ANC_REC": 0.0, "ANC0": 3.0},
        merge(reg_vrc_iv(), add_tac(Regimen()), add_steroid(Regimen(), 20.0, 84)))

    cols = ("log10_burden_d14", "log10_burden_d42", "GM_peak", "GM_d14",
            "PERF_min", "mortality_d42_pct", "mortality_d84_pct")
    say(f"{'scenario':40s}" + "".join(f"{c.replace('log10_burden_', 'logB'):>12s}"
                                      for c in cols))
    say("-" * 124)
    for k, v in scen.items():
        say(f"{k:40s}" + "".join(f"{v[c]:12.2f}" for c in cols))
    RESULTS["D_scenarios"] = scen

    # ------------------------------------------------------ time-to-therapy
    hdr("E. TIME-TO-TREATMENT BIFURCATION (VRC, NM, ANC recovery d14)")
    say(f"{'start day':>10s} {'logB d28':>10s} {'logB d84':>10s} "
        f"{'PERF min':>10s} {'mort d84 %':>12s}")
    say("-" * 78)
    delay = []
    for sd in (1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14):
        t, o = simulate({"T_ANC_REC": 336.0}, reg_vrc_iv(start_d=float(sd)),
                        tmax=84 * 24, dt=2.0)
        lb = burden_log(o)
        i28 = int(np.argmin(np.abs(t - 28 * 24)))
        rec = dict(start_day=sd,
                   logB_d28=round(float(lb[i28]), 2),
                   logB_d84=round(float(lb[-1]), 2),
                   PERF_min=round(float(o[:, IDX["PERF"]].min()), 3),
                   mort_d84=round(100 * (1 - float(survival(o)[-1])), 1))
        delay.append(rec)
        say(f"{sd:10d} {rec['logB_d28']:10.2f} {rec['logB_d84']:10.2f} "
            f"{rec['PERF_min']:10.3f} {rec['mort_d84']:12.1f}")
    RESULTS["E_start_delay"] = delay
    lo = [r for r in delay if r["mort_d84"] < 50.0]
    hi = [r for r in delay if r["mort_d84"] >= 50.0]
    if lo and hi:
        a = max(r["start_day"] for r in lo)
        b = min(r["start_day"] for r in hi)
        ma = [r["mort_d84"] for r in delay if r["start_day"] == a][0]
        mb = [r["mort_d84"] for r in delay if r["start_day"] == b][0]
        say("")
        say(f"  This is a BIFURCATION, not a gradient. Days 1-{a} form a graded arm")
        say(f"  ({delay[0]['mort_d84']:.1f}% -> {ma:.1f}% mortality); day {b} jumps to {mb:.1f}%.")
        say(f"  One day of delay ({a} -> {b}) is worth {mb - ma:.1f} mortality points,")
        say(f"  more than the entire preceding {a - 1}-day window ({ma - delay[0]['mort_d84']:.1f} points).")
        RESULTS["E_bifurcation"] = dict(last_recoverable_day=a, first_failing_day=b,
                                        mort_at_last=ma, mort_at_first=mb,
                                        one_day_cost=round(mb - ma, 1))

    # --------------------------------------------------- ANC recovery timing
    hdr("F. DRUG CHOICE vs NEUTROPHIL RECOVERY: WHICH LEVER IS BIGGER?")
    say(f"{'ANC recovery day':>18s} {'no drug':>10s} {'VRC':>10s} "
        f"{'ISA':>10s} {'L-AmB 3':>10s} {'VRC+ECH':>10s}")
    say("-" * 78)
    lever = []
    for rd in (5, 10, 14, 21, 28, 1e9):
        row = {"anc_rec_day": (None if rd > 1e8 else rd)}
        vals = []
        for tag, reg in (("none", None), ("VRC", reg_vrc_iv()),
                         ("ISA", reg_isa()), ("AMB", reg_amb(mgkg=3.0)),
                         ("VRC+ECH", merge(reg_vrc_iv(), reg_ech()))):
            t, o = simulate({"T_ANC_REC": (1e9 if rd > 1e8 else rd * 24.0)},
                            reg, tmax=84 * 24, dt=4.0)
            m = round(100 * (1 - float(survival(o)[-1])), 1)
            row[tag] = m
            vals.append(m)
        lever.append(row)
        lab = "never" if rd > 1e8 else f"{int(rd)}"
        say(f"{lab:>18s} " + "".join(f"{v:10.1f}" for v in vals))
    RESULTS["F_lever"] = lever
    # The honest comparison is the spread AMONG ACTIVE DRUGS at a fixed host,
    # against the spread across hosts at a fixed drug. Including the untreated
    # arm would measure "drug vs no drug", which is not the clinical question.
    drugs = ("VRC", "ISA", "AMB", "VRC+ECH")
    per_host = {}
    for row in lever:
        vals = [row[k] for k in drugs]
        lab = "never" if row["anc_rec_day"] is None else f"day {int(row['anc_rec_day'])}"
        per_host[lab] = round(max(vals) - min(vals), 1)
    d_drug = max(per_host.values())
    host_spans = {k: round(max(r[k] for r in lever) - min(r[k] for r in lever), 1)
                  for k in drugs}
    d_host = max(host_spans.values())
    say("")
    say("  Spread AMONG THE FOUR ACTIVE REGIMENS, at each fixed host:")
    for k, v in per_host.items():
        say(f"    ANC recovery {k:10s} : {v:5.1f} mortality points")
    say("  Spread ACROSS HOSTS, at each fixed regimen:")
    for k, v in host_spans.items():
        say(f"    {k:10s} : {v:5.1f} mortality points")
    say("")
    say(f"  Largest drug lever {d_drug:.1f} points; largest host lever {d_host:.1f} points")
    say(f"  -- the marrow is worth {d_host / max(d_drug, 1e-9):.1f}x the prescription.")
    RESULTS["F_drug_spread_by_host"] = per_host
    RESULTS["F_host_spread_by_drug"] = host_spans
    RESULTS["F_drug_effect_points"] = d_drug
    RESULTS["F_host_effect_points"] = d_host
    RESULTS["F_host_over_drug"] = round(d_host / max(d_drug, 1e-9), 1)

    # ------------------------------------------------------------ GM paradox
    hdr("G. GALACTOMANNAN MOVES FOR TWO DIFFERENT REASONS")
    say("GM production = k_growth * (growth flux) + k_lysis * (killing flux).")
    say("A fungistatic azole removes the FIRST term; a cidal polyene inflates")
    say("the SECOND before it removes it.  Week-1 GM therefore behaves")
    say("differently under the two drugs at comparable burden.")
    say("")
    say(f"{'arm':26s} {'GM d5':>8s} {'GM d7':>8s} {'GM d10':>8s} {'GM d14':>8s} "
        f"{'GM d21':>8s} {'logB d10':>9s} {'logB d14':>9s} {'logB d21':>9s}")
    say("-" * 110)
    gmrows = []
    for tag, reg in (("no therapy", None), ("voriconazole", reg_vrc_iv()),
                     ("L-AmB 3 mg/kg", reg_amb(mgkg=3.0)),
                     ("anidulafungin", reg_ech())):
        t, o = simulate({"T_ANC_REC": 1e9}, reg, tmax=24 * 24, dt=0.5)
        gm = o[:, IDX["GMser"]]
        lb = burden_log(o)
        pick_ = lambda dd: float(gm[int(np.argmin(np.abs(t - dd * 24)))])
        pb = lambda dd: float(lb[int(np.argmin(np.abs(t - dd * 24)))])
        row = dict(arm=tag, GM_d3=round(pick_(3), 2), GM_d5=round(pick_(5), 2),
                   GM_d7=round(pick_(7), 2), GM_d10=round(pick_(10), 2),
                   GM_d14=round(pick_(14), 2), GM_d21=round(pick_(21), 2),
                   logB_d5=round(pb(5), 2), logB_d7=round(pb(7), 2),
                   logB_d10=round(pb(10), 2), logB_d14=round(pb(14), 2),
                   logB_d21=round(pb(21), 2))
        gmrows.append(row)
        say(f"{tag:26s} {row['GM_d5']:8.2f} {row['GM_d7']:8.2f} {row['GM_d10']:8.2f} "
            f"{row['GM_d14']:8.2f} {row['GM_d21']:8.2f} {row['logB_d10']:9.2f} "
            f"{row['logB_d14']:9.2f} {row['logB_d21']:9.2f}")
    RESULTS["G_gm_kinetics"] = gmrows
    nt = [r for r in gmrows if r["arm"] == "no therapy"][0]
    ec = [r for r in gmrows if r["arm"] == "anidulafungin"][0]
    vr = [r for r in gmrows if r["arm"] == "voriconazole"][0]
    say("")
    say("  Two consequences of writing GM release on the FLUX and not the STOCK:")
    say("")
    say(f"  (1) UNTREATED GM FALLS WHILE THE PATIENT GETS WORSE. No therapy: GM")
    say(f"      {nt['GM_d10']:.2f} on day 10 at burden {nt['logB_d10']:.2f}, then")
    say(f"      {nt['GM_d14']:.2f} on day 14 at burden {nt['logB_d14']:.2f}. The burden")
    say(f"      did not fall -- it hit carrying capacity, so the GROWTH flux went to")
    say(f"      zero and GM production stopped. A falling GM is not a treated infection.")
    back = (ec["GM_d14"] > nt["GM_d14"]) and (ec["logB_d14"] < nt["logB_d14"])
    say("")
    say(f"  (2) DAY-14 RANKING IS INVERTED: anidulafungin GM {ec['GM_d14']:.2f} at burden")
    say(f"      {ec['logB_d14']:.2f}, untreated GM {nt['GM_d14']:.2f} at burden {nt['logB_d14']:.2f}.")
    say(f"      {'The higher GM belongs to the LOWER burden.' if back else 'Ranking preserved here.'}")
    say("")
    say(f"  (3) On voriconazole GM never crosses 0.5 (peak {max(vr['GM_d5'], vr['GM_d7'], vr['GM_d10'], vr['GM_d14'], vr['GM_d21']):.2f})")
    say(f"      while the burden sits at 10^{vr['logB_d14']:.2f}: a GM-negative live infection,")
    say(f"      which is the model's version of the reduced GM sensitivity under")
    say(f"      mould-active therapy reported by Marr 2005.")
    RESULTS["G_untreated_gm_falls_at_capacity"] = dict(
        gm_d10=nt["GM_d10"], logB_d10=nt["logB_d10"],
        gm_d14=nt["GM_d14"], logB_d14=nt["logB_d14"])
    RESULTS["G_day14_ranking_inverted"] = bool(back)
    RESULTS["G_vrc_gm_negative_live_infection"] = dict(
        gm_max=max(vr["GM_d5"], vr["GM_d7"], vr["GM_d10"], vr["GM_d14"], vr["GM_d21"]),
        logB_d14=vr["logB_d14"], logB_d21=vr["logB_d21"])

    # ----------------------------------------------------- perfusion feedback
    hdr("H. A HYPOTHESIS OF MY OWN, TESTED AND REFUTED:\n"
        "   'LATE TREATMENT FAILS BECAUSE ANGIOINVASION BLOCKS DRUG DELIVERY'")
    say("The map draws a positive feedback loop -- hyphae destroy the vessels")
    say("that deliver the drug that would kill the hyphae -- and I expected it")
    say("to BE the time-to-treatment cliff. The way to find out is to delete it:")
    say("set DELIV_FLOOR = 1 so the drug reaches the hyphae at full ELF")
    say("concentration no matter how thrombosed the lesion is, and re-run.")
    say("")
    say(f"{'start day':>10s} {'peak logB':>11s} {'peak logB':>11s} "
        f"{'mort d84':>10s} {'mort d84':>10s}")
    say(f"{'':>10s} {'(feedback)':>11s} {'(deleted)':>11s} "
        f"{'(feedback)':>10s} {'(deleted)':>10s}")
    say("-" * 78)
    fb = []
    for sd in (3, 4, 5, 6, 7, 8, 10, 14):
        _, o1 = simulate({"T_ANC_REC": 336.0}, reg_vrc_iv(start_d=float(sd)),
                         tmax=84 * 24, dt=4.0)
        _, o0 = simulate({"T_ANC_REC": 336.0, "DELIV_FLOOR": 1.0},
                         reg_vrc_iv(start_d=float(sd)), tmax=84 * 24, dt=4.0)
        r = dict(start_day=sd,
                 peak_on=round(float(burden_log(o1).max()), 2),
                 peak_off=round(float(burden_log(o0).max()), 2),
                 mort_on=round(100 * (1 - float(survival(o1)[-1])), 1),
                 mort_off=round(100 * (1 - float(survival(o0)[-1])), 1))
        fb.append(r)
        say(f"{sd:10d} {r['peak_on']:11.2f} {r['peak_off']:11.2f} "
            f"{r['mort_on']:10.1f} {r['mort_off']:10.1f}")
    RESULTS["H_delivery_ablation"] = fb
    on = [r for r in fb if r["mort_on"] >= 50.0]
    off = [r for r in fb if r["mort_off"] >= 50.0]
    b_on = min((r["start_day"] for r in on), default=None)
    b_off = min((r["start_day"] for r in off), default=None)
    say("")
    say(f"  Cliff with the feedback : first failing start day = {b_on}")
    say(f"  Cliff without it        : first failing start day = {b_off}")
    if b_on == b_off:
        say("  => REFUTED. Deleting perfusion-limited delivery does not move the")
        say("     cliff at all. The time-to-treatment threshold is NOT a delivery")
        say("     problem; it is set by how large the burden already is when a")
        say("     GROWTH-RATE inhibitor takes over, relative to what that class of")
        say("     drug can hold until the marrow returns. The angioinvasion loop is")
        say("     real in the model and it drives perfusion, oxygenation and the")
        say("     hazard -- but it is not the mechanism of the cliff, and the map")
        say("     should not be read as claiming that it is.")
    else:
        say(f"  => The feedback moves the cliff by {b_off - b_on} day(s); it is")
        say("     contributory but not the whole mechanism.")
    RESULTS["H_cliff_with_feedback"] = b_on
    RESULTS["H_cliff_without_feedback"] = b_off
    RESULTS["H_hypothesis_refuted"] = bool(b_on == b_off)

    # ---------------------------------------------------------- MIC / target
    hdr("I. PK/PD TARGET ATTAINMENT: AUC/MIC AND THE MIC BREAKPOINT")
    say(f"{'MIC (mg/L)':>11s} {'fAUC24/MIC':>12s} {'logB d42':>10s} "
        f"{'mort d84 %':>12s}")
    say("-" * 78)
    microws = []
    t, o = simulate({"INOCULUM": 0.0}, reg_vrc_iv(0.0, 14), tmax=14 * 24, dt=0.25)
    c = o[:, IDX["VRC_c"]] / P["VRC_Vc"]
    m = (t >= 12 * 24) & (t <= 13 * 24)
    auc24_vrc = auc(t[m], c[m])
    fauc = auc24_vrc * P["VRC_fu"]
    for mic in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
        t, o = simulate({"T_ANC_REC": 336.0, "MIC_VRC": mic},
                        reg_vrc_iv(), tmax=84 * 24, dt=4.0)
        lb = burden_log(o)
        i42 = int(np.argmin(np.abs(t - 42 * 24)))
        r = dict(MIC=mic, fAUC_MIC=round(fauc / mic, 1),
                 logB_d42=round(float(lb[i42]), 2),
                 mort_d84=round(100 * (1 - float(survival(o)[-1])), 1))
        microws.append(r)
        say(f"{mic:11.3f} {r['fAUC_MIC']:12.1f} {r['logB_d42']:10.2f} "
            f"{r['mort_d84']:12.1f}")
    RESULTS["I_mic_ladder"] = microws
    RESULTS["I_fAUC24_vrc"] = round(fauc, 1)

    # ------------------------------------------------------ resistance window
    hdr("J. RESISTANCE SELECTION UNDER AZOLE MONOTHERAPY")
    say(f"{'arm':30s} {'logB_R d42':>12s} {'R fraction d42':>16s} "
        f"{'logB_total d42':>16s}")
    say("-" * 84)
    resrows = []
    for tag, over, reg in (
        ("no therapy", {"T_ANC_REC": 336.0}, None),
        ("VRC monotherapy", {"T_ANC_REC": 336.0}, reg_vrc_iv()),
        ("VRC, persistent neutropenia", {"T_ANC_REC": 1e9}, reg_vrc_iv()),
        ("L-AmB 3 mg/kg", {"T_ANC_REC": 336.0}, reg_amb(mgkg=3.0)),
        ("VRC + anidulafungin", {"T_ANC_REC": 336.0}, merge(reg_vrc_iv(), reg_ech())),
    ):
        t, o = simulate(over, reg, tmax=84 * 24, dt=4.0)
        i42 = int(np.argmin(np.abs(t - 42 * 24)))
        hw, hr = o[i42, IDX["HYW"]], o[i42, IDX["HYR"]]
        tot = max(hw + hr, 1.0)
        r = dict(arm=tag, logB_R=round(float(np.log10(max(hr, 1e-30))), 2),
                 R_frac=float(f"{hr / tot:.3e}"),
                 logB_tot=round(float(np.log10(tot)), 2))
        resrows.append(r)
        say(f"{tag:30s} {r['logB_R']:12.2f} {r['R_frac']:16.3e} {r['logB_tot']:16.2f}")
    RESULTS["J_resistance"] = resrows

    # ------------------------------------------------------ virtual population
    hdr("K. VIRTUAL POPULATION (n=320): WHAT PREDICTS DEATH?")
    rng = np.random.default_rng(20260806)
    n = 320
    geno_keys = list(CYP2C19.keys())
    geno_p = [0.04, 0.26, 0.40, 0.25, 0.05]
    rows = []
    for i in range(n):
        gk = geno_keys[int(rng.choice(len(geno_keys), p=geno_p))]
        anc_rec = float(rng.choice([5, 8, 10, 14, 21, 28, 1e9],
                                   p=[0.10, 0.16, 0.20, 0.20, 0.16, 0.12, 0.06]))
        start = float(rng.choice([2, 3, 4, 6, 9], p=[0.20, 0.30, 0.25, 0.15, 0.10]))
        mic = float(rng.choice([0.25, 0.5, 1.0, 2.0, 8.0],
                               p=[0.30, 0.42, 0.18, 0.06, 0.04]))
        steroid = bool(rng.random() < 0.35)
        over = {"CYP2C19": CYP2C19[gk],
                "T_ANC_REC": (1e9 if anc_rec > 1e8 else anc_rec * 24.0),
                "MIC_VRC": mic,
                "INOCULUM": float(10 ** rng.normal(3.0, 0.35))}
        reg = reg_vrc_iv(start_d=start)
        if steroid:
            reg = add_steroid(reg, 40.0, 45)
        t, o = simulate(over, reg, tmax=84 * 24, dt=12.0)
        rows.append(dict(geno=gk, anc_rec=anc_rec, start=start, mic=mic,
                         steroid=steroid,
                         mort=100 * (1 - float(survival(o)[-1])),
                         logB84=float(burden_log(o)[-1])))
        if (i + 1) % 100 == 0:
            print(f"    ... {i + 1}/{n}", file=sys.stderr)

    overall = float(np.mean([r["mort"] for r in rows]))
    say(f"Overall predicted 12-week mortality: {overall:.1f}%")
    say("")

    def strat(name, keyfn, order=None):
        say(f"  by {name}:")
        groups = OrderedDict()
        for r in rows:
            groups.setdefault(keyfn(r), []).append(r["mort"])
        keys = order or sorted(groups.keys(), key=lambda z: str(z))
        out_ = OrderedDict()
        for k in keys:
            if k not in groups:
                continue
            v = groups[k]
            out_[str(k)] = dict(n=len(v), mortality=round(float(np.mean(v)), 1))
            say(f"    {str(k):24s} n={len(v):4d}  mortality {np.mean(v):5.1f}%")
        return out_

    pop = OrderedDict()
    pop["by_anc_recovery"] = strat(
        "neutrophil recovery day",
        lambda r: ("never" if r["anc_rec"] > 1e8 else f"day {int(r['anc_rec'])}"),
        ["day 5", "day 8", "day 10", "day 14", "day 21", "day 28", "never"])
    pop["by_genotype"] = strat("CYP2C19 genotype", lambda r: r["geno"], geno_keys)
    pop["by_start_day"] = strat("therapy start day",
                                lambda r: f"day {int(r['start'])}",
                                ["day 2", "day 3", "day 4", "day 6", "day 9"])
    pop["by_mic"] = strat("voriconazole MIC", lambda r: f"MIC {r['mic']}",
                          ["MIC 0.25", "MIC 0.5", "MIC 1.0", "MIC 2.0", "MIC 8.0"])
    pop["by_steroid"] = strat("corticosteroid exposure",
                              lambda r: "steroid" if r["steroid"] else "no steroid",
                              ["no steroid", "steroid"])
    pop["overall_mortality_pct"] = round(overall, 1)
    RESULTS["K_virtual_population"] = pop

    # rank the covariates by spread
    spreads = OrderedDict()
    for lbl, blk in (("neutrophil recovery", pop["by_anc_recovery"]),
                     ("therapy start day", pop["by_start_day"]),
                     ("voriconazole MIC", pop["by_mic"]),
                     ("CYP2C19 genotype", pop["by_genotype"]),
                     ("corticosteroid", pop["by_steroid"])):
        vals = [v["mortality"] for v in blk.values()]
        spreads[lbl] = round(max(vals) - min(vals), 1)
    say("")
    say("  Covariate spread in predicted 12-week mortality (percentage points):")
    for k, v in sorted(spreads.items(), key=lambda kv: -kv[1]):
        say(f"    {k:24s} {v:6.1f}")
    RESULTS["K_covariate_spread"] = spreads

    # --------------------------------------------------------------- TDM value
    hdr("L. WHAT TDM IS WORTH (VRC, ANC recovery d14, MIC 0.5)")
    say(f"{'strategy':34s} {'AUC24':>10s} {'logB d42':>10s} {'ALT max':>9s} "
        f"{'neuro AUC':>11s} {'mort d84':>10s}")
    say("-" * 90)
    tdm = []
    for tag, over, mgkg in (
        ("UM, fixed 4 mg/kg", {"CYP2C19": CYP2C19["UM (*17/*17)"]}, 4.0),
        ("UM, TDM-escalated 7.6 mg/kg", {"CYP2C19": CYP2C19["UM (*17/*17)"]}, 7.6),
        ("NM, fixed 4 mg/kg", {"CYP2C19": CYP2C19["NM (*1/*1)"]}, 4.0),
        ("PM, fixed 4 mg/kg", {"CYP2C19": CYP2C19["PM (*2/*2)"]}, 4.0),
        ("PM, TDM-reduced 1.8 mg/kg", {"CYP2C19": CYP2C19["PM (*2/*2)"]}, 1.8),
    ):
        o_ = dict(over)
        o_["T_ANC_REC"] = 336.0
        t, o = simulate(o_, reg_vrc_iv(maint_mg_kg=mgkg), tmax=84 * 24, dt=2.0)
        c = o[:, IDX["VRC_c"]] / P["VRC_Vc"]
        m = (t >= 20 * 24) & (t <= 21 * 24)
        a = auc(t[m], c[m])
        i42 = int(np.argmin(np.abs(t - 42 * 24)))
        r = dict(strategy=tag, AUC24=round(a, 1),
                 logB_d42=round(float(burden_log(o)[i42]), 2),
                 ALT_max=round(float(o[:, IDX["ALT"]].max()), 0),
                 neuro_AUC=round(auc(t, o[:, IDX["NEUROE"]]), 1),
                 mort_d84=round(100 * (1 - float(survival(o)[-1])), 1))
        tdm.append(r)
        say(f"{tag:34s} {r['AUC24']:10.1f} {r['logB_d42']:10.2f} {r['ALT_max']:9.0f} "
            f"{r['neuro_AUC']:11.1f} {r['mort_d84']:10.1f}")
    RESULTS["L_tdm"] = tdm

    # ------------------------------------------------------------------- DDI
    hdr("M. AZOLE-TACROLIMUS INTERACTION (transplant scenario)")
    say(f"{'arm':40s} {'TAC trough d14 (ng/mL)':>24s} {'fold vs alone':>16s}")
    say("-" * 84)
    ddi = []
    base_tac = None
    for tag, over, reg in (
        ("tacrolimus alone", {"T_ANC_REC": 0.0, "ANC0": 3.0}, add_tac(Regimen())),
        ("tacrolimus + voriconazole", {"T_ANC_REC": 0.0, "ANC0": 3.0},
         merge(add_tac(Regimen()), reg_vrc_iv(start_d=0.0))),
        ("tacrolimus + isavuconazole", {"T_ANC_REC": 0.0, "ANC0": 3.0},
         merge(add_tac(Regimen()), reg_isa(start_d=0.0))),
        ("tacrolimus 1/3 dose + voriconazole", {"T_ANC_REC": 0.0, "ANC0": 3.0},
         merge(add_tac(Regimen(), mg_bid=1.0), reg_vrc_iv(start_d=0.0))),
    ):
        t, o = simulate(over, reg, tmax=20 * 24, dt=0.5)
        i = int(np.argmin(np.abs(t - 14 * 24)))
        tr = float(o[i, IDX["TAC_c"]] / P["TAC_Vc"] * 1000.0)
        if base_tac is None:
            base_tac = tr
        ddi.append(dict(arm=tag, trough=round(tr, 2), fold=round(tr / base_tac, 2)))
        say(f"{tag:40s} {tr:24.2f} {tr / base_tac:16.2f}")
    RESULTS["M_ddi"] = ddi

    # ------------------------------------------------------------------- CNS
    hdr("N. CEREBRAL ASPERGILLOSIS: PENETRATION DECIDES, NOT POTENCY")
    say(f"{'arm':34s} {'brain Cavg':>12s} {'Kp':>7s} {'log brain burden d42':>22s}")
    say("-" * 84)
    cns = []
    for tag, reg, kp, key in (
        ("voriconazole 4 mg/kg q12h", reg_vrc_iv(), P["VRC_Kp_brain"], "VRC_brain"),
        ("L-AmB 5 mg/kg", reg_amb(mgkg=5.0), P["AMB_Kp_brain"], "AMB_brain"),
        ("isavuconazole 200 mg qd", reg_isa(), P["ISA_Kp_brain"], None),
        ("anidulafungin 100 mg qd", reg_ech(), P["ECH_Kp_brain"], None),
        ("no therapy", None, 0.0, None),
    ):
        t, o = simulate({"T_ANC_REC": 336.0}, reg, tmax=84 * 24, dt=4.0)
        i42 = int(np.argmin(np.abs(t - 42 * 24)))
        if key:
            cb = float(np.mean(o[(t > 20 * 24), IDX[key]]))
        elif tag.startswith("isavu"):
            cb = float(np.mean(o[(t > 20 * 24), IDX["ISA_c"]] / P["ISA_Vc"]) * kp)
        elif tag.startswith("anidula"):
            cb = float(np.mean(o[(t > 20 * 24), IDX["ECH_c"]] / P["ECH_Vc"]) * kp)
        else:
            cb = 0.0
        brb = float(o[i42, IDX["BRB"]])
        lbb = float(np.log10(brb)) if brb >= 1.0 else None
        cns.append(dict(arm=tag, brain_Cavg=round(cb, 3), Kp=kp,
                        brain_burden_d42=brb,
                        log_brain_burden_d42=(round(lbb, 2) if lbb is not None else "cleared")))
        shown = f"{lbb:22.2f}" if lbb is not None else f"{'cleared (<1 CFUe)':>22s}"
        say(f"{tag:34s} {cb:12.3f} {kp:7.3f} {shown}")
    RESULTS["N_cns"] = cns

    # ------------------------------------------------------------- mass check
    hdr("O. NUMERICAL QUALITY CONTROL")
    t, o = simulate({"T_ANC_REC": 240.0}, reg_vrc_iv(), tmax=84 * 24, dt=1.0)
    neg = int(np.sum(o < -1e-6))
    say(f"negative states encountered : {neg}")
    say(f"max |dPERF| out of [0,1]    : "
        f"{max(0.0, float(o[:, IDX['PERF']].max() - 1.0)):.2e}")
    say(f"survival monotone decreasing: {bool(np.all(np.diff(survival(o)) <= 1e-12))}")
    # tight vs loose tolerance
    t2, o2 = simulate({"T_ANC_REC": 240.0}, reg_vrc_iv(), tmax=84 * 24, dt=4.0)
    lb1 = burden_log(o)[int(np.argmin(np.abs(t - 42 * 24)))]
    lb2 = burden_log(o2)[int(np.argmin(np.abs(t2 - 42 * 24)))]
    say(f"grid-independence logB(d42) : {lb1:.4f} (dt=1h) vs {lb2:.4f} (dt=4h), "
        f"delta {abs(lb1 - lb2):.2e}")
    RESULTS["O_qc"] = dict(negative_states=neg,
                           grid_delta_logB_d42=round(float(abs(lb1 - lb2)), 6),
                           survival_monotone=bool(np.all(np.diff(survival(o)) <= 1e-12)))

    # ------------------------------------------------------------------ write
    with open("ipa_scenario_results.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    with open("ipa_reference_output.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")
    say("")
    say("wrote ipa_scenario_results.json and ipa_reference_output.txt")


if __name__ == "__main__":
    main()
