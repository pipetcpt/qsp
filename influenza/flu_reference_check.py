#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
flu_reference_check.py
======================
Independent numpy/scipy transcription of the influenza A QSP model that is
shipped in this directory as `flu_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
The build environment for this repository has no R toolchain.  Publishing an
ODE model that has never been integrated is worthless, so the same equation
system has been written twice: once in mrgsolve/C++ (`flu_mrgsolve_model.R`)
and once here in scipy LSODA.  Every number quoted in `README.md` is produced
by a function in this file.  The two transcriptions are meant to be the same
system with the same parameter block; if they ever disagree, one of them is
wrong.

THE STRUCTURAL CLAIM
--------------------
"Start antivirals within 48 hours" is normally taught as a rule about the
clock.  This model poses it instead as a rule about a RESOURCE:

      the epithelial target-cell pool T(t) is the substrate the virus
      consumes, and by the time a patient feels ill enough to seek care
      most of it has already been consumed.

Everything else follows.  A drug arriving at time t_rx cannot recover the
cells that are already dead; it can only protect the remainder.  So the
achievable benefit is bounded above by T(t_rx)/T(0), and potency can only
move you toward that bound -- never past it.  This is why the timing axis and
the potency axis are NOT interchangeable, and the model is built so that this
falls out of the equations rather than being asserted.

The second claim is about resistance.  A resistant subpopulation does not grow
because the drug "creates" it -- it is there at mutation-selection balance
from the first replication cycle.  It grows because the drug removes its
competitor and hands it the target cells the wild type would have taken.
Resistance emergence is therefore COMPETITIVE RELEASE, and its magnitude is a
non-monotonic function of drug potency: a useless drug releases nothing, a
perfect drug suppresses both strains, and the maximum sits in between.

OPERATOR CLASSIFICATION
-----------------------
Every therapy is classified by WHICH TERM of the replication loop it touches.
The model's claim is that these classes are not interchangeable, because they
enter the equations at different points and therefore have different
signatures in time:

  TRANSCRIPTION   blocks E -> I maturation            baloxavir (PA endonuclease)
  PRODUCTION      lowers p, virions made per cell     baloxavir, favipiravir
  RELEASE         lowers p, virions freed per cell    oseltamivir/zanamivir/peramivir
  ENTRY           lowers beta, infection of new cells mucosal IgA, mAb, (NAI, weak)
  MUTAGENESIS     raises mutation rate to lethality   favipiravir
  VIRION CLEARANCE raises c                           mAb, convalescent plasma
  TARGET PROTECTION converts T -> R (refractory)      interferon, ISG induction
  INFECTED-CELL CLEARANCE raises delta                CD8 CTL, NK
  IMMUNOPATHOLOGY damping of the above (harmful here) corticosteroids

USAGE
-----
    python3 flu_reference_check.py             # everything (A0 - A13)
    python3 flu_reference_check.py --only A4   # one analysis
    python3 flu_reference_check.py --list      # list analyses
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# 0.  STATE VECTOR
# ---------------------------------------------------------------------------
# Mixed-unit in-host convention (Baccam 2006 J Virol 80:7590, PMID 16840338):
# cell compartments are absolute cell counts, virus is a titre in
# TCID50/mL of nasal-wash (URT) or BAL (LRT) equivalent, and the infectivity
# constant beta absorbs the sampling volume.  This is the convention every
# published human influenza viral-dynamic fit uses, so the parameters below
# are directly comparable to the literature.

STATES = [
    # --- upper respiratory tract -------------------------------------------
    "TU",    #  0 susceptible epithelial cells, URT                    [cells]
    "RU",    #  1 refractory (IFN-protected) cells, URT                [cells]
    "EW",    #  2 eclipse cells, wild type, URT                        [cells]
    "IW",    #  3 productively infected, wild type, URT                [cells]
    "EM",    #  4 eclipse cells, resistant mutant, URT                 [cells]
    "IM",    #  5 productively infected, mutant, URT                   [cells]
    "VW",    #  6 wild-type virus, URT                          [TCID50/mL]
    "VM",    #  7 mutant virus, URT                             [TCID50/mL]
    "DU",    #  8 dead/denuded epithelium, URT                         [cells]
    # --- lower respiratory tract -------------------------------------------
    "TL",    #  9 susceptible epithelial cells, LRT                    [cells]
    "RL",    # 10 refractory cells, LRT                                [cells]
    "ELW",   # 11 eclipse, wild type, LRT                              [cells]
    "ILW",   # 12 productively infected, wild type, LRT                [cells]
    "ELM",   # 13 eclipse, mutant, LRT                                 [cells]
    "ILM",   # 14 productively infected, mutant, LRT                   [cells]
    "VLW",   # 15 wild-type virus, LRT                          [TCID50/mL]
    "VLM",   # 16 mutant virus, LRT                             [TCID50/mL]
    "DL",    # 17 dead/denuded alveolar epithelium, LRT                [cells]
    # --- immune ------------------------------------------------------------
    "FU",    # 18 type I interferon, URT                        [arb, 0-1ish]
    "FL",    # 19 type I interferon, LRT                        [arb]
    "AG",    # 20 antigen load presented in draining node       [arb]
    "CTL",   # 21 effector CD8 T cells                          [arb, 1=naive]
    "PC",    # 22 antibody-secreting plasma cells               [arb]
    "AB",    # 23 serum neutralising IgG                        [arb]
    "IGA",   # 24 mucosal secretory IgA                         [arb]
    "NEU",   # 25 neutrophil infiltrate, LRT                    [arb]
    "IL6",   # 26 systemic pro-inflammatory cytokine (IL-6 proxy)  [pg/mL]
    # --- clinical ----------------------------------------------------------
    "SYM",   # 27 composite symptom score (0-21, CAPSTONE 7x0-3)
    "TMP",   # 28 temperature elevation above 37.0              [degC]
    "BAC",   # 29 bacterial (pneumococcal) burden, log-scale carrier [log10 CFU]
    # --- drug PK -----------------------------------------------------------
    "OSd",   # 30 oseltamivir phosphate, gut depot              [mg base]
    "OSp",   # 31 oseltamivir prodrug, plasma                   [mg]
    "OCc",   # 32 oseltamivir carboxylate, central              [mg]
    "OCp",   # 33 oseltamivir carboxylate, peripheral           [mg]
    "OCe",   # 34 oseltamivir carboxylate, ELF/effect site      [ng/mL]
    "BXd",   # 35 baloxavir marboxil, gut depot                 [mg]
    "BXc",   # 36 baloxavir acid, central                       [mg]
    "BXp",   # 37 baloxavir acid, peripheral                    [mg]
    "BXe",   # 38 baloxavir acid, ELF/effect site (free)        [ng/mL]
    "PRc",   # 39 peramivir, central                            [mg]
    "PRp",   # 40 peramivir, peripheral                         [mg]
    "FVd",   # 41 favipiravir, gut depot                        [mg]
    "FVc",   # 42 favipiravir, central                          [mg]
    "MBc",   # 43 anti-HA monoclonal antibody, central          [mg]
    "MBp",   # 44 anti-HA monoclonal antibody, peripheral       [mg]
    "STd",   # 45 dexamethasone, gut depot                      [mg]
    "STc",   # 46 dexamethasone, central                        [mg]
    # --- bookkeeping -------------------------------------------------------
    "AUCU",  # 47 cumulative log10 URT viral AUC                [log10.d]
    "AUCL",  # 48 cumulative log10 LRT viral AUC                [log10.d]
    "CUMK",  # 49 cumulative epithelial cells killed, URT       [cells]
]
IX = {name: i for i, name in enumerate(STATES)}
NSTATE = len(STATES)

# ---------------------------------------------------------------------------
# 1.  PARAMETERS
# ---------------------------------------------------------------------------
# Provenance tags in the comments:
#   [LIT]  taken directly from a published estimate
#   [CAL]  calibrated in this model against a named clinical anchor
#   [ASM]  structural assumption, no direct measurement
P0 = dict(
    # ---- URT viral dynamics ------------------------------------------------
    T0U      = 4.0e8,    # [LIT] URT epithelial target cells, Baccam 2006
    BETA     = 1.63e-5,   # [CAL] infectivity (mL/TCID50/d), Baccam 2.7e-5 + eclipse
    KECL     = 4.0,      # [LIT] eclipse -> productive rate, 1/d (6 h)
    DELTA    = 2.0,      # [LIT] death rate of productive cells, 1/d, Baccam 2006
    PVIR     = 0.114,   # [CAL] virion production, TCID50/mL/cell/d
    CVIR     = 2.4,      # [LIT] virion clearance, 1/d, Baccam 2006
    LREG     = 0.10,     # [ASM] epithelial regeneration, 1/d (~7 d half-repair)
    EXTC     = 1.0,      # [ASM] extinction floor, infected cells (see rhs)
    # ---- LRT viral dynamics ------------------------------------------------
    T0L      = 1.0e9,    # [ASM] accessible LRT epithelial targets
    BETAL    = 1.05e-6,   # [ASM] LRT infectivity (lower: mucus, surfactant, SP-D)
    PVIRL    = 1.8e-2,   # [ASM] LRT production
    DELTAL   = 1.6,      # [ASM] LRT infected-cell death
    LREGL    = 0.055,    # [ASM] alveolar repair is slower than nasal
    ADESC    = 0.008,    # [ASM] URT -> LRT descent/aspiration flux, 1/d
    # ---- interferon / refractory state ------------------------------------
    QF       = 3.1e-8,   # [CAL] IFN induction per infected cell, 1/d
    DF       = 2.0,      # [LIT] IFN decay, 1/d
    PHIF     = 2.6,      # [CAL] T -> R conversion rate at saturating IFN, 1/d
    KF       = 0.45,     # [ASM] IFN EC50 for the antiviral state
    RHOR     = 0.35,     # [ASM] R -> T reversion, 1/d
    # ---- adaptive ----------------------------------------------------------
    QAG      = 6.0,      # [ASM] antigen appearance per unit infected fraction
    DAG      = 1.2,      # [ASM] antigen decay, 1/d
    RCTL     = 3.1,      # [CAL] CTL expansion rate, 1/d
    KAG      = 0.09,     # [ASM] antigen EC50 for CTL expansion
    CTLMAX   = 120.0,    # [ASM] CTL carrying capacity (fold over naive)
    DCTL     = 0.32,     # [ASM] CTL contraction, 1/d
    KKILL    = 0.003,     # [CAL] CTL-mediated infected-cell killing, 1/d per unit
    QPC      = 8.0,     # [ASM] plasma-cell recruitment, 1/d
    DPC      = 0.25,     # [ASM] plasma-cell loss, 1/d
    KABP     = 3.0,      # [ASM] IgG secretion per plasma cell, 1/d
    DAB      = 0.05,     # [ASM] IgG decay, 1/d (t1/2 ~14 d)
    KIGA     = 0.30,     # [ASM] mucosal IgA secretion, 1/d
    DIGA     = 0.03,     # [ASM] IgA decay, 1/d
    KNEUT    = 0.9,      # [ASM] antibody neutralisation of free virion, 1/d/unit
    KABENT   = 0.30,     # [ASM] antibody block of entry, 1/unit
    KIGAENT  = 0.55,     # [ASM] IgA block of entry, 1/unit
    # ---- innate cellular / damage -----------------------------------------
    QNEU     = 2.4,      # [ASM] neutrophil recruitment per unit LRT damage rate
    DNEU     = 0.9,      # [ASM] neutrophil clearance, 1/d
    # ---- cytokine / symptoms ----------------------------------------------
    Q6F      = 200.0,     # [CAL] IL-6 production per unit URT IFN, pg/mL/d
    Q6L      = 400.0,    # [ASM] IL-6 production per unit LRT IFN, pg/mL/d
    Q6D      = 9.0e-9,   # [ASM] IL-6 per LRT cell killed per day
    D6       = 5.5,      # [LIT] IL-6 elimination, 1/d (t1/2 ~3 h)
    IL60     = 1.6,      # [LIT] healthy baseline IL-6, pg/mL
    KON      = 3.2,      # [CAL] symptom onset rate, 1/d
    KOFF     = 2.00,     # [CAL] symptom resolution rate, 1/d
    K6       = 11.0,     # [CAL] IL-6 EC50 for symptom drive, pg/mL
    SMAX     = 21.0,     # [LIT] max composite score (7 symptoms x 0-3)
    SALLEV   = 7.0,      # [LIT] alleviation ceiling (all 7 symptoms <= 1)
    WVIR     = 0.60,      # [CAL] fraction of symptom drive tracking titre (see A8)
    LODS     = 0.5,      # [ASM] titre below which the viral drive is zero
    LREFS    = 6.0,      # [ASM] titre span over which the viral drive saturates
    SONSET   = 2.0,      # [ASM] score at which a patient calls themselves ill
    TALLEV   = 21.5/24., # [LIT] alleviation must be sustained 21.5 h (CAPSTONE)
    KTMP     = 8.5,      # [CAL] fever gain, degC/d
    KTMPD    = 2.6,      # [ASM] defervescence, 1/d
    KT6      = 30.0,     # [ASM] IL-6 EC50 for fever, pg/mL
    # ---- bacterial superinfection -----------------------------------------
    BAC0     = 2.0,      # [ASM] baseline log10 CFU nasopharyngeal carriage
    BMAXB    = 9.0,      # [ASM] log10 CFU carrying capacity
    KBG      = 1.05,     # [ASM] bacterial growth, log10/d at full permissiveness
    KBADH    = 6.0,      # [ASM] adhesion gain from exposed sialic acid / damage
    KBCLR    = 1.05 * (1.0 - 2.0 / 9.0),   # fixed by the steady-state identity
    KBIFN    = 0.75,     # [LIT] type I IFN suppression of antibacterial defence
    BACTHR   = 6.0,      # [ASM] log10 CFU threshold = clinical superinfection
    # ---- mutation / fitness ------------------------------------------------
    MU       = 2.5e-5,   # [LIT] per-site per-replication mutation rate, influenza
    NU       = 2.5e-5,   # [LIT] back mutation
    COST     = 0.18,     # [LIT] replicative fitness cost of PA/I38T (production)
    VM0FRAC  = 0.0,      # pre-existing mutant fraction in the inoculum
    # ---- oseltamivir PK ----------------------------------------------------
    KA_OS    = 12.0,     # [LIT] absorption, 1/d... set in per-day units below
    KCONV    = 26.0,     # [LIT] prodrug -> carboxylate hydrolysis (CES1), 1/d
    CLOP     = 100.0,    # [LIT] prodrug non-converting clearance, L/d
    VOP      = 30.0,     # [LIT] prodrug volume, L
    CLOC     = 451.2,    # [LIT] OC clearance 18.8 L/h -> L/d
    VOC      = 28.0,     # [LIT] OC central volume, L
    VOC2     = 30.0,     # [ASM] OC peripheral volume, L
    QOC      = 60.0,     # [ASM] OC intercompartmental clearance, L/d
    FOS      = 0.80,     # [LIT] fraction of dose converted to OC
    RELF_OC  = 1.00,     # [LIT] ELF:plasma ratio for OC
    KEQ_OC   = 24.0,     # [ASM] ELF equilibration, 1/d
    # ---- baloxavir PK ------------------------------------------------------
    KA_BX    = 19.2,     # [LIT] absorption, 1/d (Tmax ~4 h)
    CLBX     = 149.8,    # [LIT] CL/F 6.24 L/h -> L/d
    VBX      = 380.0,    # [LIT] central volume, L
    VBX2     = 420.0,    # [LIT] peripheral volume, L
    QBX      = 500.0,    # [ASM] intercompartmental clearance, L/d
    FUBX     = 0.07,     # [LIT] free fraction (93% protein bound)
    RELF_BX  = 1.00,     # [ASM] free drug equilibrates with ELF
    KEQ_BX   = 24.0,     # [ASM] ELF equilibration, 1/d
    # ---- peramivir PK ------------------------------------------------------
    CLPR     = 264.0,    # [LIT] 11 L/h -> L/d
    VPR      = 12.6,     # [LIT] Vss, L
    VPR2     = 8.0,      # [ASM]
    QPR      = 40.0,     # [ASM]
    # ---- favipiravir PK ----------------------------------------------------
    KA_FV    = 24.0,     # [LIT] rapid absorption
    VFV      = 15.0,     # [LIT] volume, L
    VMFV     = 9000.0,   # [LIT] Michaelis-Menten Vmax, mg/d (auto-inhibited)
    KMFV     = 45.0,     # [LIT] Km, mg/L
    # ---- monoclonal antibody PK -------------------------------------------
    CLMB     = 0.25,     # [LIT] L/d (t1/2 ~21 d)
    VMB      = 3.2,      # [LIT] central volume, L
    VMB2     = 3.0,      # [ASM]
    QMB      = 0.6,      # [ASM]
    # ---- dexamethasone PK --------------------------------------------------
    KA_ST    = 24.0,
    CLST     = 200.0,    # [LIT] ~8.3 L/h
    VST      = 60.0,     # [LIT] L
    # ---- pharmacodynamics --------------------------------------------------
    #   in-vivo EFFECTIVE EC50s -- calibrated, NOT in-vitro values.  The gap
    #   between these and the in-vitro numbers is reported in analysis A2.
    EC50_NAI = 3.0,     # [CAL] OC ELF conc for 50% release blockade, ng/mL
    EMAX_NAI = 0.995,    # [CAL] max fractional blockade of release
    HILL_NAI = 1.0,      # [ASM]
    FBETA_NAI= 0.30,     # [ASM] NAI entry-block is 30% of its release block
    EC50_BX  = 0.052,    # [CAL] free BXA for 50% transcription block, ng/mL
    EMAX_BX  = 0.9999,    # [CAL]
    HILL_BX  = 2.0,      # [ASM]
    ETA_BX   = 1.00,     # [ASM] transcription block also blocks E -> I fully
    FMATU    = 1.00,     # idealised productive-fraction operator (A5 only)
    EC50_PR  = 3.0,     # [ASM] peramivir shares the NAI effective EC50
    EMAX_PR  = 0.995,
    EC50_FV  = 25.0,     # [ASM] favipiravir, mg/L
    EMAX_FV  = 0.93,
    FMUTA_FV = 0.55,     # [LIT] fraction of favipiravir progeny made non-viable
    KMAB_C   = 0.55,     # [ASM] mAb-driven virion clearance, 1/d per ug/mL
    KMAB_B   = 0.020,    # [ASM] mAb entry blockade, per ug/mL
    ECST     = 12.0,     # [ASM] dexamethasone conc for half-max immunosuppression
    IMAX_ST  = 0.65,     # [LIT] max suppression of IFN induction and CTL killing
    # ---- resistance profile (set per scenario) ----------------------------
    RF_NAI   = 1.0,      # fold shift in EC50_NAI carried by the mutant
    RF_BX    = 1.0,      # fold shift in EC50_BX carried by the mutant
    RF_FV    = 1.0,
    # ---- host phenotype modifiers -----------------------------------------
    FIFN     = 1.0,      # multiplier on IFN induction (age / immune status)
    FCTL     = 1.0,      # multiplier on CTL response
    FAB      = 1.0,      # multiplier on humoral response
    IGA0     = 0.0,      # pre-existing mucosal IgA (vaccination / prior infection)
    CTL0     = 1.0,      # pre-existing cross-reactive CD8 memory
    # ---- inoculum ----------------------------------------------------------
    V0       = 1.0e-2,   # [LIT] TCID50/mL, human challenge inoculum equivalent
)

# scipy works in days; oseltamivir ka of 0.5/h -> 12/d already encoded above.


# ---------------------------------------------------------------------------
# 2.  DOSING
# ---------------------------------------------------------------------------
class Dose:
    """A single bolus into a state, at a time, in the state's own units."""
    __slots__ = ("t", "state", "amt")

    def __init__(self, t, state, amt):
        self.t = float(t)
        self.state = state
        self.amt = float(amt)


def regimen_oseltamivir(t0, days=5.0, mg=75.0):
    """75 mg PO BID x 5 d, expressed as oseltamivir base."""
    return [Dose(t0 + k * 0.5, "OSd", mg) for k in range(int(round(days * 2)))]


def regimen_baloxavir(t0, mg=40.0):
    """Single oral dose (40 mg if <80 kg)."""
    return [Dose(t0, "BXd", mg)]


def regimen_peramivir(t0, mg=600.0):
    return [Dose(t0, "PRc", mg)]


def regimen_favipiravir(t0, days=5.0):
    """1800 mg BID day 1, then 800 mg BID."""
    d = [Dose(t0, "FVd", 1800.0), Dose(t0 + 0.5, "FVd", 1800.0)]
    for k in range(2, int(round(days * 2))):
        d.append(Dose(t0 + 0.5 * k, "FVd", 800.0))
    return d


def regimen_mab(t0, mg=3600.0):
    return [Dose(t0, "MBc", mg)]


def regimen_dexamethasone(t0, days=5.0, mg=6.0):
    return [Dose(t0 + k, "STd", mg) for k in range(int(round(days)))]


# ---------------------------------------------------------------------------
# 3.  RIGHT-HAND SIDE
# ---------------------------------------------------------------------------
def _hill(c, ec50, emax, h=1.0):
    if c <= 0.0:
        return 0.0
    ch = c ** h
    return emax * ch / (ch + ec50 ** h)


def rhs(t, y, P):
    y = np.maximum(y, 0.0)
    (TU, RU, EW, IW, EM, IM, VW, VM, DU,
     TL, RL, ELW, ILW, ELM, ILM, VLW, VLM, DL,
     FU, FL, AG, CTL, PC, AB, IGA, NEU, IL6,
     SYM, TMP, BAC,
     OSd, OSp, OCc, OCp, OCe,
     BXd, BXc, BXp, BXe,
     PRc, PRp,
     FVd, FVc,
     MBc, MBp,
     STd, STc,
     AUCU, AUCL, CUMK) = y

    d = np.zeros(NSTATE)

    # ---------------- drug concentrations ----------------------------------
    C_OC   = OCc / P["VOC"] * 1000.0        # ng/mL  (mg/L * 1000)
    C_OCE  = OCe                            # ng/mL, ELF
    C_BX   = BXc / P["VBX"] * 1000.0        # ng/mL total
    C_BXE  = BXe                            # ng/mL free, ELF
    C_PR   = PRc / P["VPR"] * 1000.0        # ng/mL
    C_FV   = FVc / P["VFV"]                 # mg/L
    C_MB   = MBc / P["VMB"]                 # mg/L = ug/mL
    C_ST   = STc / P["VST"] * 1000.0        # ng/mL

    # ---------------- pharmacodynamic operators ----------------------------
    # RELEASE / PRODUCTION blockade from the neuraminidase inhibitors.
    # Wild type sees EC50 as-is; the mutant sees it shifted by RF_NAI.
    e_nai_W = _hill(C_OCE, P["EC50_NAI"], P["EMAX_NAI"], P["HILL_NAI"]) \
            + _hill(C_PR,  P["EC50_PR"],  P["EMAX_PR"])
    e_nai_M = _hill(C_OCE, P["EC50_NAI"] * P["RF_NAI"], P["EMAX_NAI"], P["HILL_NAI"]) \
            + _hill(C_PR,  P["EC50_PR"] * P["RF_NAI"], P["EMAX_PR"])
    e_nai_W = min(e_nai_W, P["EMAX_NAI"])
    e_nai_M = min(e_nai_M, P["EMAX_NAI"])

    # TRANSCRIPTION blockade from baloxavir (cap-snatching / PA endonuclease)
    e_bx_W = _hill(C_BXE, P["EC50_BX"], P["EMAX_BX"], P["HILL_BX"])
    e_bx_M = _hill(C_BXE, P["EC50_BX"] * P["RF_BX"], P["EMAX_BX"], P["HILL_BX"])

    # POLYMERASE blockade + lethal mutagenesis from favipiravir
    e_fv_W = _hill(C_FV, P["EC50_FV"], P["EMAX_FV"])
    e_fv_M = _hill(C_FV, P["EC50_FV"] * P["RF_FV"], P["EMAX_FV"])

    # combined multiplicative blockade of virion OUTPUT per infected cell
    prod_W = (1.0 - e_nai_W) * (1.0 - e_bx_W) * (1.0 - e_fv_W)
    prod_M = (1.0 - e_nai_M) * (1.0 - e_bx_M) * (1.0 - e_fv_M)

    # blockade of the eclipse -> productive transition (transcription only)
    matu_W = P["FMATU"] * (1.0 - P["ETA_BX"] * e_bx_W)
    matu_M = P["FMATU"] * (1.0 - P["ETA_BX"] * e_bx_M)

    # fraction of released progeny that are non-infectious (favipiravir)
    inf_W = 1.0 - P["FMUTA_FV"] * e_fv_W
    inf_M = 1.0 - P["FMUTA_FV"] * e_fv_M

    # ENTRY blockade: NAI (weak), antibody, IgA, mAb
    ent_blk = min(0.995,
                  P["FBETA_NAI"] * e_nai_W
                  + P["KABENT"] * AB
                  + P["KIGAENT"] * IGA
                  + P["KMAB_B"] * C_MB)
    ent = 1.0 - ent_blk

    # VIRION CLEARANCE augmentation: antibody + mAb
    cW = P["CVIR"] * (1.0 + P["KNEUT"] * AB) + P["KMAB_C"] * C_MB
    cM = cW

    # IMMUNOPATHOLOGY damping: corticosteroid suppresses IFN and CTL killing
    supp = _hill(C_ST, P["ECST"], P["IMAX_ST"])
    ifn_gain  = P["FIFN"] * (1.0 - supp)
    kill_gain = P["FCTL"] * (1.0 - supp)

    # ---------------- URT viral dynamics -----------------------------------
    infW = P["BETA"] * ent * TU * VW
    infM = P["BETA"] * ent * TU * VM
    Itot = IW + IM
    ctl_kill = P["KKILL"] * kill_gain * CTL
    to_ref = P["PHIF"] * FU / (FU + P["KF"]) * TU
    from_ref = P["RHOR"] * RU
    occupied = TU + RU + EW + IW + EM + IM
    regen = P["LREG"] * max(0.0, P["T0U"] - occupied)

    d[IX["TU"]]  = -infW - infM - to_ref + from_ref + regen
    d[IX["RU"]]  = to_ref - from_ref
    # An eclipse cell always LEAVES the eclipse state at rate KECL.  A
    # transcription blocker does not freeze it there -- it decides where the
    # cell goes: a fraction matu becomes productive, the remainder is
    # abortively infected and counted as lost epithelium.  Writing it the other
    # way (blocking the exit) creates a spurious eclipse reservoir that keeps
    # seeding productive cells for days after the drug has taken hold.
    d[IX["EW"]]  = infW - P["KECL"] * EW - P["DELTA"] * 0.25 * EW
    d[IX["IW"]]  = P["KECL"] * matu_W * EW - (P["DELTA"] + ctl_kill) * IW
    d[IX["EM"]]  = infM - P["KECL"] * EM - P["DELTA"] * 0.25 * EM
    d[IX["IM"]]  = P["KECL"] * matu_M * EM - (P["DELTA"] + ctl_kill) * IM
    abortU = P["KECL"] * ((1.0 - matu_W) * EW + (1.0 - matu_M) * EM)

    # EXTINCTION FLOOR.  A deterministic ODE lets a strain fall to 1e-20 cells
    # and grow back when targets regenerate ("atto-fox").  Real infections at
    # that burden are extinct.  Each strain's release is therefore scaled by a
    # smooth establishment factor that is ~1 above one infected cell and ->0
    # below it.  Without it the epithelial regeneration term produces a
    # spurious second wave; with it, clearance is final.
    xiW = (EW + IW + ELW + ILW) / ((EW + IW + ELW + ILW) + P["EXTC"])
    xiM = (EM + IM + ELM + ILM) / ((EM + IM + ELM + ILM) + P["EXTC"])

    relW = P["PVIR"] * prod_W * IW * xiW
    relM = P["PVIR"] * prod_M * (1.0 - P["COST"]) * IM * xiM
    # mutation is applied at release: a fraction MU of wild-type progeny carry
    # the resistance allele, and NU of mutant progeny revert.
    # NOTE ON UNITS.  In the mixed-unit in-host convention (cells absolute,
    # virus as a titre per mL of wash) the absorption term -beta*T*V has units
    # of cells/day and CANNOT be subtracted from a titre.  Baccam 2006 and every
    # subsequent human influenza fit omit it; so does this model.  Loss of
    # virions to cell entry is instead absorbed into the clearance constant c.
    d[IX["VW"]] = (relW * (1.0 - P["MU"]) * inf_W + relM * P["NU"] * inf_M
                   - cW * VW - P["ADESC"] * VW)
    d[IX["VM"]] = (relM * (1.0 - P["NU"]) * inf_M + relW * P["MU"] * inf_W
                   - cM * VM - P["ADESC"] * VM)
    killedU = (P["DELTA"] + ctl_kill) * Itot + abortU
    d[IX["DU"]] = killedU - P["LREG"] * DU
    d[IX["CUMK"]] = killedU

    # ---------------- LRT viral dynamics -----------------------------------
    infLW = P["BETAL"] * ent * TL * VLW
    infLM = P["BETAL"] * ent * TL * VLM
    ILtot = ILW + ILM
    to_refL = P["PHIF"] * FL / (FL + P["KF"]) * TL
    from_refL = P["RHOR"] * RL
    occupiedL = TL + RL + ELW + ILW + ELM + ILM
    regenL = P["LREGL"] * max(0.0, P["T0L"] - occupiedL)

    d[IX["TL"]]  = -infLW - infLM - to_refL + from_refL + regenL
    d[IX["RL"]]  = to_refL - from_refL
    d[IX["ELW"]] = infLW - P["KECL"] * ELW - P["DELTAL"] * 0.25 * ELW
    d[IX["ILW"]] = P["KECL"] * matu_W * ELW - (P["DELTAL"] + ctl_kill) * ILW
    d[IX["ELM"]] = infLM - P["KECL"] * ELM - P["DELTAL"] * 0.25 * ELM
    d[IX["ILM"]] = P["KECL"] * matu_M * ELM - (P["DELTAL"] + ctl_kill) * ILM
    abortL = P["KECL"] * ((1.0 - matu_W) * ELW + (1.0 - matu_M) * ELM)

    relLW = P["PVIRL"] * prod_W * ILW * xiW
    relLM = P["PVIRL"] * prod_M * (1.0 - P["COST"]) * ILM * xiM
    d[IX["VLW"]] = (relLW * (1.0 - P["MU"]) * inf_W + relLM * P["NU"] * inf_M
                    - cW * VLW + P["ADESC"] * VW)
    d[IX["VLM"]] = (relLM * (1.0 - P["NU"]) * inf_M + relLW * P["MU"] * inf_W
                    - cM * VLM + P["ADESC"] * VM)
    killedL = (P["DELTAL"] + ctl_kill) * ILtot + abortL
    d[IX["DL"]] = killedL - P["LREGL"] * DL

    # ---------------- interferon -------------------------------------------
    d[IX["FU"]] = P["QF"] * ifn_gain * Itot - P["DF"] * FU
    d[IX["FL"]] = P["QF"] * ifn_gain * ILtot - P["DF"] * FL

    # ---------------- adaptive ---------------------------------------------
    ag_drive = (Itot / P["T0U"]) + (ILtot / P["T0L"])
    d[IX["AG"]]  = P["QAG"] * ag_drive - P["DAG"] * AG
    d[IX["CTL"]] = (P["RCTL"] * AG / (AG + P["KAG"]) * CTL
                    * (1.0 - CTL / P["CTLMAX"]) - P["DCTL"] * (CTL - P["CTL0"]))
    d[IX["PC"]]  = P["QPC"] * P["FAB"] * AG - P["DPC"] * PC
    d[IX["AB"]]  = P["KABP"] * PC - P["DAB"] * AB
    d[IX["IGA"]] = P["KIGA"] * PC - P["DIGA"] * (IGA - P["IGA0"])
    d[IX["NEU"]] = P["QNEU"] * killedL / P["T0L"] - P["DNEU"] * NEU

    # ---------------- cytokine, symptoms, fever ----------------------------
    d[IX["IL6"]] = (P["Q6F"] * FU + P["Q6L"] * FL + P["Q6D"] * killedL
                    - P["D6"] * (IL6 - P["IL60"]))
    drive6 = max(0.0, IL6 - P["IL60"])
    f_cyt = drive6 / (drive6 + P["K6"])
    # WVIR is the fraction of the symptom drive that tracks the INSTANTANEOUS
    # viral titre rather than the cytokine state.  With WVIR = 0 the symptom
    # score is a pure function of the cytokine response, whose size and timing
    # are fixed before any post-peak drug can act -- and then no antiviral
    # started after the viral peak can shorten the illness.  Analysis A8 solves
    # for the WVIR that the CAPSTONE-1 symptom endpoint requires.
    lv = math.log10(max(VW + VM, 1e-6))
    f_vir = min(1.0, max(0.0, (lv - P["LODS"]) / P["LREFS"]))
    f_sym = (1.0 - P["WVIR"]) * f_cyt + P["WVIR"] * f_vir
    d[IX["SYM"]] = P["KON"] * f_sym * (P["SMAX"] - SYM) - P["KOFF"] * SYM
    d[IX["TMP"]] = (P["KTMP"] * drive6 / (drive6 + P["KT6"]) - P["KTMPD"] * TMP)

    # ---------------- bacterial superinfection -----------------------------
    # Two influenza-created conditions multiply: exposed adhesion sites on
    # denuded epithelium (permiss), and type I interferon suppression of
    # antibacterial defence (neu_fn < 1 lowers clearance).
    # KBCLR is fixed by the requirement that BAC = BAC0 be a steady state in
    # the uninfected host:  KBCLR = KBG * (1 - BAC0/BMAXB).
    dam_frac = min(1.0, DL / P["T0L"] + DU / P["T0U"])
    permiss = 1.0 + P["KBADH"] * dam_frac
    neu_fn = 1.0 / (1.0 + P["KBIFN"] * (FU + FL) / P["KF"])
    d[IX["BAC"]] = (P["KBG"] * permiss * (1.0 - BAC / P["BMAXB"])
                    - P["KBCLR"] * neu_fn * (BAC / P["BAC0"]))

    # ---------------- PK ----------------------------------------------------
    d[IX["OSd"]] = -P["KA_OS"] * OSd
    d[IX["OSp"]] = P["KA_OS"] * OSd - (P["KCONV"] + P["CLOP"] / P["VOP"]) * OSp
    conv = P["KCONV"] * OSp * P["FOS"] * (284.4 / 312.4)
    d[IX["OCc"]] = (conv - P["CLOC"] / P["VOC"] * OCc
                    - P["QOC"] * (OCc / P["VOC"] - OCp / P["VOC2"]))
    d[IX["OCp"]] = P["QOC"] * (OCc / P["VOC"] - OCp / P["VOC2"])
    d[IX["OCe"]] = P["KEQ_OC"] * (C_OC * P["RELF_OC"] - OCe)

    d[IX["BXd"]] = -P["KA_BX"] * BXd
    d[IX["BXc"]] = (P["KA_BX"] * BXd - P["CLBX"] / P["VBX"] * BXc
                    - P["QBX"] * (BXc / P["VBX"] - BXp / P["VBX2"]))
    d[IX["BXp"]] = P["QBX"] * (BXc / P["VBX"] - BXp / P["VBX2"])
    d[IX["BXe"]] = P["KEQ_BX"] * (C_BX * P["FUBX"] * P["RELF_BX"] - BXe)

    d[IX["PRc"]] = (-P["CLPR"] / P["VPR"] * PRc
                    - P["QPR"] * (PRc / P["VPR"] - PRp / P["VPR2"]))
    d[IX["PRp"]] = P["QPR"] * (PRc / P["VPR"] - PRp / P["VPR2"])

    d[IX["FVd"]] = -P["KA_FV"] * FVd
    d[IX["FVc"]] = P["KA_FV"] * FVd - P["VMFV"] * C_FV / (P["KMFV"] + C_FV)

    d[IX["MBc"]] = (-P["CLMB"] / P["VMB"] * MBc
                    - P["QMB"] * (MBc / P["VMB"] - MBp / P["VMB2"]))
    d[IX["MBp"]] = P["QMB"] * (MBc / P["VMB"] - MBp / P["VMB2"])

    d[IX["STd"]] = -P["KA_ST"] * STd
    d[IX["STc"]] = P["KA_ST"] * STd - P["CLST"] / P["VST"] * STc

    # ---------------- bookkeeping ------------------------------------------
    d[IX["AUCU"]] = math.log10(max(VW + VM, 1e-6) / 1e-6) if (VW + VM) > 1e-6 else 0.0
    d[IX["AUCL"]] = math.log10(max(VLW + VLM, 1e-6) / 1e-6) if (VLW + VLM) > 1e-6 else 0.0

    return d


# ---------------------------------------------------------------------------
# 4.  INITIAL CONDITIONS AND INTEGRATION
# ---------------------------------------------------------------------------
def y0_of(P):
    y = np.zeros(NSTATE)
    y[IX["TU"]] = P["T0U"]
    y[IX["TL"]] = P["T0L"]
    y[IX["VW"]] = P["V0"] * (1.0 - P["VM0FRAC"])
    y[IX["VM"]] = P["V0"] * P["VM0FRAC"]
    y[IX["CTL"]] = P["CTL0"]
    y[IX["IGA"]] = P["IGA0"]
    y[IX["IL6"]] = P["IL60"]
    y[IX["BAC"]] = P["BAC0"]
    return y


def simulate(P=None, doses=None, tmax=21.0, dt=1.0 / 96.0, **over):
    """Integrate with bolus doses applied by restarting the solver."""
    par = dict(P0)
    if P:
        par.update(P)
    par.update(over)
    doses = sorted(doses or [], key=lambda dd: dd.t)

    tgrid = np.arange(0.0, tmax + 1e-9, dt)
    out = np.zeros((len(tgrid), NSTATE))
    y = y0_of(par)

    # break points: dose times inside (0, tmax)
    breaks = sorted({0.0} | {dd.t for dd in doses if 0.0 < dd.t < tmax} | {tmax})
    filled = 0
    for a, b in zip(breaks[:-1], breaks[1:]):
        for dd in doses:
            if abs(dd.t - a) < 1e-12:
                y[IX[dd.state]] += dd.amt
        seg = tgrid[(tgrid >= a - 1e-12) & (tgrid <= b + 1e-12)]
        if len(seg) == 0:
            seg = np.array([a, b])
        sol = solve_ivp(rhs, (a, max(b, a + 1e-9)), y, args=(par,),
                        method="LSODA", t_eval=seg,
                        rtol=1e-8, atol=1e-10, max_step=0.05)
        if not sol.success:
            raise RuntimeError(sol.message)
        n = sol.y.shape[1]
        idx = np.searchsorted(tgrid, seg[0] - 1e-12)
        out[idx:idx + n, :] = sol.y.T
        filled = max(filled, idx + n)
        y = np.maximum(sol.y[:, -1].copy(), 0.0)
    if filled < len(tgrid):
        out[filled:, :] = out[filled - 1, :]
    return tgrid, out, par


# ---------------------------------------------------------------------------
# 5.  READOUTS
# ---------------------------------------------------------------------------
LOD = 0.5     # log10 TCID50/mL limit of detection (CAPSTONE-1 assay floor)


def col(out, name):
    return out[:, IX[name]]


def vlog(out, upper=True):
    v = col(out, "VW") + col(out, "VM") if upper else col(out, "VLW") + col(out, "VLM")
    return np.log10(np.maximum(v, 1e-6))


def peak(t, x):
    i = int(np.argmax(x))
    return t[i], x[i]


def symptom_onset(t, out, P):
    s = col(out, "SYM")
    idx = np.where(s >= P["SONSET"])[0]
    return float(t[idx[0]]) if len(idx) else float("nan")


def time_to_alleviation(t, out, P, t_start):
    """CAPSTONE rule: all seven symptoms mild/absent, sustained for 21.5 h."""
    s = col(out, "SYM")
    dt = t[1] - t[0]
    need = int(round(P["TALLEV"] / dt))
    lo = int(np.searchsorted(t, t_start))
    ok = s <= P["SALLEV"]
    run = 0
    for i in range(lo, len(t)):
        run = run + 1 if ok[i] else 0
        if run >= need:
            return (t[i - need + 1] - t_start) * 24.0     # hours from start
    return float("nan")


def time_to_shed_stop(t, out, P, t_start, thr=LOD):
    v = vlog(out)
    lo = int(np.searchsorted(t, t_start))
    for i in range(lo, len(t)):
        if v[i] <= thr and np.all(v[i:] <= thr):
            return (t[i] - t_start) * 24.0
    return float("nan")


def titre_drop_at(t, out, t_start, hours):
    v = vlog(out)
    i0 = int(np.searchsorted(t, t_start))
    i1 = int(np.searchsorted(t, t_start + hours / 24.0))
    i1 = min(i1, len(t) - 1)
    return v[i1] - v[i0]


def fever_duration(t, out, thr=0.7):
    """Hours with temperature elevation >= 0.7 degC (i.e. >= 37.7)."""
    tm = col(out, "TMP")
    dt = t[1] - t[0]
    return float(np.sum(tm >= thr) * dt * 24.0)


def mutant_fraction(out):
    """Mutant share of the URT titre, defined only where the total titre is
    still detectable.  Below the assay floor the ratio is arithmetically valid
    but clinically empty, and reporting it there manufactures 100% mutant
    readings out of virus nobody could sample."""
    vw, vm = col(out, "VW"), col(out, "VM")
    tot = vw + vm
    return np.where(tot > 10.0 ** LOD, vm / np.maximum(tot, 1e-30), 0.0)


def residual_auc(t, out):
    """R(t) = integral from t to the end of the detectable log-titre.
    This is the exact upper bound on what any antiviral started at t can
    remove: a drug cannot subtract viral load that has already happened."""
    v = np.maximum(vlog(out) - LOD, 0.0)
    dt = t[1] - t[0]
    return np.concatenate([np.cumsum(v[::-1])[::-1] * dt])


def summarise(t, out, P, t_rx):
    v = vlog(out)
    vl = vlog(out, upper=False)
    mf = mutant_fraction(out)
    tp, vp = peak(t, v)
    detect = v > LOD
    return dict(
        peak_log=vp, t_peak=tp,
        auc_log=float(np.trapezoid(np.maximum(v - LOD, 0.0), t)),
        drop24=titre_drop_at(t, out, t_rx, 24.0),
        drop48=titre_drop_at(t, out, t_rx, 48.0),
        shed_h=time_to_shed_stop(t, out, P, t_rx),
        ttas_h=time_to_alleviation(t, out, P, t_rx),
        fever_h=fever_duration(t, out),
        peak_sym=float(np.max(col(out, "SYM"))),
        peak_il6=float(np.max(col(out, "IL6"))),
        lrt_peak=float(np.max(vl)),
        epi_lost=float(np.max(col(out, "CUMK")) / P["T0U"]),
        tmin_frac=float(np.min(col(out, "TU")) / P["T0U"]),
        mut_frac_end=float(mf[-1]),
        mut_frac_max=float(np.max(mf[t >= t_rx])) if np.any(t >= t_rx) else 0.0,
        mut_peak_log=float(np.max(np.log10(np.maximum(col(out, "VM"), 1e-6)))),
        bac_peak=float(np.max(col(out, "BAC"))),
        shed_days=float(np.sum(detect) * (t[1] - t[0])),
    )


# ---------------------------------------------------------------------------
# 6.  SCENARIOS
# ---------------------------------------------------------------------------
#   t_onset  is discovered from the placebo run and reused, so that every arm
#   is randomised at the same clock time -- exactly as a trial would.
RX_DELAY_H = 24.0        # median time from symptom onset to first dose, CAPSTONE-1

I38T = dict(RF_BX=50.0, RF_NAI=1.0, COST=0.18)      # baloxavir resistance
H275Y = dict(RF_NAI=300.0, RF_BX=1.0, COST=0.28)    # oseltamivir resistance


def placebo_onset(P=None, **over):
    t, out, par = simulate(P, [], tmax=21.0, **over)
    return symptom_onset(t, out, par)


def build_scenarios(t_rx):
    """{name: (doses, parameter overrides, reference time for the endpoints)}.

    The reference time is the arm's own first dose, so that a timing variant is
    scored from when its drug was actually given rather than from the common
    randomisation clock; arms with no drug use the common clock."""
    onset = t_rx - RX_DELAY_H / 24.0
    return {
        "1. Placebo":
            ([], {}, t_rx),
        "2. Oseltamivir 75 mg BID x5d":
            (regimen_oseltamivir(t_rx), {}, t_rx),
        "3. Baloxavir 40 mg single":
            (regimen_baloxavir(t_rx), {}, t_rx),
        "4. Peramivir 600 mg IV single":
            (regimen_peramivir(t_rx), {}, t_rx),
        "5. Favipiravir 1800/800 mg BID x5d":
            (regimen_favipiravir(t_rx), {}, t_rx),
        "6. Baloxavir + oseltamivir":
            (regimen_baloxavir(t_rx) + regimen_oseltamivir(t_rx), {}, t_rx),
        "7. Anti-HA mAb single IV":
            (regimen_mab(t_rx), {}, t_rx),
        "8. Baloxavir at symptom onset":
            (regimen_baloxavir(onset), {}, onset),
        "9. Baloxavir late (72 h post-onset)":
            (regimen_baloxavir(onset + 3.0), {}, onset + 3.0),
        "10. Oseltamivir + dexamethasone 6 mg":
            (regimen_oseltamivir(t_rx) + regimen_dexamethasone(t_rx), {}, t_rx),
        "11. Baloxavir, immunocompromised host":
            (regimen_baloxavir(t_rx), dict(FCTL=0.12, FAB=0.15, FIFN=0.55), t_rx),
        "12. Baloxavir, PA/I38T resistance profile":
            (regimen_baloxavir(t_rx), dict(I38T), t_rx),
        "13. Oseltamivir, H275Y resistance profile":
            (regimen_oseltamivir(t_rx), dict(H275Y), t_rx),
        "14. Vaccinated host (mucosal IgA), no drug":
            ([], dict(IGA0=0.55, CTL0=6.0), t_rx),
    }


# ---------------------------------------------------------------------------
# 7.  ANALYSES
# ---------------------------------------------------------------------------
def _fmt(x, n=2):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "  --  "
    return f"{x:.{n}f}"


def A0_baseline():
    print(__doc__.split("USAGE")[0].strip()[:0] or "", end="")
    print("=" * 78)
    print("A0  UNTREATED NATURAL HISTORY -- calibration check")
    print("=" * 78)
    t, out, P = simulate()
    v = vlog(out)
    tp, vp = peak(t, v)
    ons = symptom_onset(t, out, P)
    s = col(out, "SYM")
    ip = int(np.argmax(s))
    tmn = float(np.min(col(out, "TU")) / P["T0U"])
    rows = [
        ("peak URT titre",        f"{vp:.2f} log10 TCID50/mL", "6.0-7.0 (Baccam 2006)"),
        ("time to peak titre",    f"{tp*24:.1f} h p.i.",       "48-72 h (Carrat 2008)"),
        ("shedding duration >LOD",f"{np.sum(v>LOD)*(t[1]-t[0]):.1f} d", "4.0-5.5 d"),
        ("symptom onset",         f"{ons*24:.1f} h p.i.",      "24-48 h"),
        ("peak symptom score",    f"{s[ip]:.1f} / 21",         "10-15"),
        ("time to peak symptoms", f"{t[ip]*24:.1f} h p.i.",    "48-72 h"),
        ("peak temperature",      f"{37.0+np.max(col(out,'TMP')):.2f} C", "38.5-39.5"),
        ("fever duration",        f"{fever_duration(t,out):.0f} h",  "48-96 h"),
        ("peak IL-6",             f"{np.max(col(out,'IL6')):.0f} pg/mL", "20-100"),
        ("min target-cell pool",  f"{tmn*100:.1f} % of T0",    "<20% (target-limited)"),
        ("peak LRT titre",        f"{np.max(vlog(out,False)):.2f} log10",  "2-4 (uncomplicated)"),
        ("peak CD8 expansion",    f"{np.max(col(out,'CTL')):.0f} x naive", "20-100 x, d7-10"),
    ]
    print(f"{'quantity':30s} {'model':>26s}   target")
    for a, b, c in rows:
        print(f"{a:30s} {b:>26s}   {c}")
    print()
    print("  The row that matters is 'min target-cell pool'.  The susceptible")
    print("  pool bottoms out at %.1f%% of its starting value -- the infection is" % (tmn*100))
    print("  target-cell-limited, and that is what makes the bound in A1 bite.")
    return t, out, P


def A1_residual_bound():
    print("=" * 78)
    print("A1  THE EXACT UPPER BOUND ON ANTIVIRAL BENEFIT  (the central result)")
    print("=" * 78)
    print("  An antiviral cannot subtract viral load that has already happened.")
    print("  So for a drug started at t_rx the reduction in viral AUC is bounded")
    print("  above by the RESIDUAL AUC")
    print()
    print("        R(t_rx) = integral from t_rx to infinity of (log10 V - LOD)+")
    print()
    print("  which is a property of the UNTREATED trajectory alone -- it contains")
    print("  no drug, no potency, no mechanism.  Everything a molecule can ever")
    print("  achieve at that time sits underneath it.\n")
    t, out, P = simulate(tmax=30.0)
    ons = symptom_onset(t, out, P)
    R = residual_auc(t, out)
    TU = col(out, "TU") / P["T0U"]
    total = R[0]
    print(f"  symptom onset at {ons*24:.1f} h p.i.;  untreated total AUC "
          f"{total:.2f} log10.d\n")
    print(f"  {'h post-onset':>12s} {'h p.i.':>7s} {'T left':>8s} {'R(t)':>8s} "
          f"{'R/total':>8s} {'achieved':>9s} {'efficiency':>11s}")
    for h in (0, 12, 24, 36, 48, 72, 96):
        trx = ons + h / 24.0
        i = min(int(np.searchsorted(t, trx)), len(t) - 1)
        t2, o2, P2 = simulate(doses=regimen_baloxavir(trx), tmax=30.0)
        ach = total - float(np.trapezoid(np.maximum(vlog(o2) - LOD, 0.0), t2))
        eff = ach / R[i] * 100 if R[i] > 1e-9 else float("nan")
        print(f"  {h:>12d} {trx*24:>7.1f} {TU[i]*100:>7.1f}% {R[i]:>8.2f} "
              f"{R[i]/total*100:>7.1f}% {ach:>9.2f} {eff:>10.1f}%")
    i0 = min(int(np.searchsorted(t, ons)), len(t) - 1)
    i24 = min(int(np.searchsorted(t, ons + 1.0)), len(t) - 1)
    i48 = min(int(np.searchsorted(t, ons + 2.0)), len(t) - 1)
    print()
    print(f"  At the moment the patient first feels ill, {R[i0]/total*100:.0f}% of the viral")
    print(f"  AUC is still in the future.  Twenty-four hours later -- the median")
    print(f"  enrolment time in CAPSTONE-1 -- it is {R[i24]/total*100:.0f}%; at 48 h, the edge of")
    print(f"  the licensed window, {R[i48]/total*100:.0f}%.  The 48-hour rule is not a statement")
    print("  about pharmacology.  It is the interval over which the quantity a")
    print("  drug could act on falls by roughly an order of magnitude, and it")
    print("  would still be there if the drug were infinitely potent.")
    print()
    print("  The efficiency column answers a different question: given the bound,")
    print("  how much of it does a real molecule actually take?  It is NOT flat.")
    print("  It is worst just before the viral peak and best long after it,")
    print("  because the cohort of cells already infected when the drug arrives")
    print("  goes on producing for a full cell lifetime no matter what the drug")
    print("  does to the cells it has not reached yet.  Early dosing therefore")
    print("  wins on the bound and loses on the efficiency, and the bound wins:")
    print("  the largest absolute benefit is still at the earliest dose time.")


def A2_pd_calibration():
    print("=" * 78)
    print("A2  PD CALIBRATION AND THE IN-VITRO / IN-VIVO GAP  (a reported discrepancy)")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    P = dict(P0)

    # what the calibrated EC50s imply for achieved exposure
    _, o_os, par = simulate(doses=regimen_oseltamivir(t_rx), tmax=21.0)
    _, o_bx, _   = simulate(doses=regimen_baloxavir(t_rx), tmax=21.0)
    tt = np.arange(0.0, 21.0 + 1e-9, 1 / 96.)
    coce = col(o_os, "OCe")
    cbxe = col(o_bx, "BXe")
    cbxt = col(o_bx, "BXc") / P["VBX"] * 1000.0

    e_os = np.array([_hill(c, P["EC50_NAI"], P["EMAX_NAI"], P["HILL_NAI"]) for c in coce])
    e_bx = np.array([_hill(c, P["EC50_BX"], P["EMAX_BX"], P["HILL_BX"]) for c in cbxe])

    MW_BXA = 483.9
    ec50_invitro_nM = 1.4                    # median in-vitro EC50, influenza A
    ec50_invitro_ng = ec50_invitro_nM * MW_BXA / 1000.0
    print(f"  oseltamivir carboxylate, ELF")
    print(f"    Cmax                          {np.max(coce):8.1f} ng/mL   "
          f"(reported plasma Cmax ~350 ng/mL)")
    print(f"    calibrated in-vivo EC50       {P['EC50_NAI']:8.1f} ng/mL")
    print(f"    peak fractional release block {np.max(e_os)*100:8.1f} %")
    print(f"    trough fractional block (48-120 h) "
          f"{np.min(e_os[(tt>t_rx+2)&(tt<t_rx+5)])*100:6.1f} %")
    print()
    print(f"  baloxavir acid")
    print(f"    total Cmax                    {np.max(cbxt):8.1f} ng/mL   "
          f"(reported 68.9 ng/mL for 40 mg)")
    print(f"    free Cmax (fu = {P['FUBX']:.2f})          {np.max(cbxe):8.2f} ng/mL "
          f"= {np.max(cbxe)/MW_BXA*1000:.1f} nM")
    print(f"    in-vitro EC50                 {ec50_invitro_ng:8.3f} ng/mL "
          f"= {ec50_invitro_nM:.1f} nM")
    print(f"    calibrated in-vivo EC50       {P['EC50_BX']:8.3f} ng/mL "
          f"= {P['EC50_BX']/MW_BXA*1000:.2f} nM")
    print(f"    ratio in-vivo / in-vitro      {P['EC50_BX']/ec50_invitro_ng:8.2f} x")
    print(f"    Hill slope (calibrated)       {P['HILL_BX']:8.1f}")
    print(f"    peak fractional block         {np.max(e_bx)*100:8.4f} %")
    print(f"    residual production fraction  {1-np.max(e_bx):8.2e}")
    print()
    print("  DISCREPANCY 1.  The CAPSTONE-1 day-2 titre fall CANNOT be produced by a")
    print("  Michaelis (Hill = 1) concentration-response at any Emax whatsoever.")
    print("  With h = 1 the residual production fraction is bounded below by")
    print("  EC50/(C + EC50), so at the achieved free concentration the fall")
    print("  saturates near -3.3 log10 in 24 h however large Emax is made -- the")
    print("  published value is -4.8.  Only a steeper-than-Michaelis slope")
    print(f"  reaches it; this model uses h = {P['HILL_BX']:.1f}, calibrated, and reports the")
    print("  requirement rather than hiding it inside Emax.  Note what that means:")
    print("  the trial datum is evidence about the SHAPE of the concentration-")
    print("  response, not only about the potency, and a QSP model fitted on Emax")
    print("  alone would have absorbed a shape error into a potency parameter.")
    print()
    print(f"  Alongside it, the calibrated in-vivo EC50 is {ec50_invitro_ng/P['EC50_BX']:.0f}x LOWER than the")
    print("  in-vitro EC50 on a free-drug basis.  Either free plasma drug")
    print("  understates what reaches the intracellular compartment where the PA")
    print("  endonuclease sits, or the nasal-wash titre falls faster than")
    print("  infectious-centre biology alone because non-infectious particle")
    print("  output is silenced too.  The model cannot choose between them.")


def A3_trial_ledger():
    print("=" * 78)
    print("A3  SCENARIO LEDGER -- model vs published trial endpoints")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    print(f"  symptom onset {ons*24:.1f} h p.i.;  first dose {t_rx*24:.1f} h p.i. "
          f"({RX_DELAY_H:.0f} h post-onset)\n")
    scen = build_scenarios(t_rx)
    print(f"  {'scenario':42s} {'TTAS':>6s} {'shed':>6s} {'d24':>6s} "
          f"{'AUC':>6s} {'peak':>6s} {'mutpk':>6s} {'fever':>6s} {'LRT':>5s}")
    print(f"  {'':42s} {'h':>6s} {'h':>6s} {'log':>6s} {'log.d':>6s} {'log':>6s} "
          f"{'log':>6s} {'h':>6s} {'log':>5s}")
    res = {}
    for name, (doses, over, tref) in scen.items():
        t, out, P = simulate(doses=doses, tmax=30.0, **over)
        r = summarise(t, out, P, tref)
        res[name] = r
        print(f"  {name:42s} {_fmt(r['ttas_h'],1):>6s} {_fmt(r['shed_h'],1):>6s} "
              f"{_fmt(r['drop24'],2):>6s} {_fmt(r['auc_log'],2):>6s} "
              f"{_fmt(r['peak_log'],2):>6s} {_fmt(r['mut_peak_log'],2):>6s} "
              f"{_fmt(r['fever_h'],0):>6s} {_fmt(r['lrt_peak'],1):>5s}")
    print()
    pl = res["1. Placebo"]
    os_ = res["2. Oseltamivir 75 mg BID x5d"]
    bx = res["3. Baloxavir 40 mg single"]
    print("  published anchors (CAPSTONE-1, Hayden 2018 NEJM 379:913, PMID 30184455):")
    print(f"    TTAS   placebo 80.2 h | oseltamivir 53.8 h | baloxavir 53.7 h")
    print(f"    model  placebo {pl['ttas_h']:.1f} h | oseltamivir {os_['ttas_h']:.1f} h "
          f"| baloxavir {bx['ttas_h']:.1f} h")
    print(f"    shed   placebo 96.0 h | oseltamivir 72.0 h | baloxavir 24.0 h")
    print(f"    model  placebo {pl['shed_h']:.1f} h | oseltamivir {os_['shed_h']:.1f} h "
          f"| baloxavir {bx['shed_h']:.1f} h")
    print(f"    day-2 titre change  placebo -1.3 | oseltamivir -2.8 | baloxavir -4.8 log10")
    print(f"    model               placebo {pl['drop24']:.2f} | oseltamivir "
          f"{os_['drop24']:.2f} | baloxavir {bx['drop24']:.2f} log10")
    return res


def A4_timing_vs_potency():
    print("=" * 78)
    print("A4  TIMING AND POTENCY ARE NOT THE SAME KIND OF VARIABLE")
    print("=" * 78)
    print("  A1 established the bound: a drug started at t_rx can remove at most")
    print("  R(t_rx), the untreated residual AUC.  That gives the two axes")
    print("  different characters, and this analysis measures both HEADROOMS")
    print("  rather than asserting which is larger:")
    print()
    print("    potency headroom  =  R(t_rx) - achieved(t_rx, Emax)")
    print("        how much of the bound a more potent drug could still take")
    print("    timing headroom   =  R(t_earlier) - R(t_rx)")
    print("        how much the bound itself moves if you dose earlier\n")
    ons = placebo_onset()
    t0, o0, P0_ = simulate(tmax=30.0)
    Rt = residual_auc(t0, o0)
    total = Rt[0]

    print("  (a) TIMING SWEEP  (baloxavir 40 mg at the calibrated potency)")
    print(f"      {'h post-onset':>12s} {'T at dose':>10s} {'bound R':>8s} "
          f"{'achieved':>9s} {'left on':>9s} {'TTAS h':>7s} {'shed h':>7s}")
    print(f"      {'':12s} {'':10s} {'log.d':>8s} {'log.d':>9s} {'table':>9s} "
          f"{'':>7s} {'':>7s}")
    tim = []
    for h in (0, 6, 12, 24, 36, 48, 72, 96):
        t_rx = ons + h / 24.0
        i = min(int(np.searchsorted(t0, t_rx)), len(t0) - 1)
        t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=30.0)
        ach = total - float(np.trapezoid(np.maximum(vlog(out) - LOD, 0.0), t))
        r = summarise(t, out, P, t_rx)
        tim.append((h, Rt[i], ach, r))
        print(f"      {h:>12d} {col(o0,'TU')[i]/P['T0U']*100:>9.1f}% {Rt[i]:>8.2f} "
              f"{ach:>9.2f} {Rt[i]-ach:>9.2f} {_fmt(r['ttas_h'],1):>7s} "
              f"{_fmt(r['shed_h'],1):>7s}")

    print()
    print("  (b) POTENCY SWEEP  (dose fixed at 24 h post-onset)")
    t_rx = ons + 1.0
    i_rx = min(int(np.searchsorted(t0, t_rx)), len(t0) - 1)
    print(f"      {'Emax':>9s} {'achieved':>9s} {'% of bound':>11s} {'TTAS h':>7s} "
          f"{'shed h':>7s}")
    pot = []
    for em in (0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.997, 0.9999, 0.999999):
        t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=30.0, EMAX_BX=em)
        ach = total - float(np.trapezoid(np.maximum(vlog(out) - LOD, 0.0), t))
        r = summarise(t, out, P, t_rx)
        pot.append((em, ach, r))
        print(f"      {em:>9.6f} {ach:>9.2f} {ach/Rt[i_rx]*100:>10.1f}% "
              f"{_fmt(r['ttas_h'],1):>7s} {_fmt(r['shed_h'],1):>7s}")

    ach_cal = [a for e, a, _ in pot if abs(e - 0.9999) < 1e-9][0]
    ach_inf = pot[-1][1]
    Rt0 = Rt[min(int(np.searchsorted(t0, ons)), len(t0) - 1)]
    print()
    print(f"  At the median enrolment time the bound is {Rt[i_rx]:.2f} log.d and the")
    print(f"  calibrated drug already takes {ach_cal:.2f} of it ({ach_cal/Rt[i_rx]*100:.1f}%).")
    print(f"    potency headroom (Emax -> 1)        {ach_inf-ach_cal:+.2f} log.d")
    print(f"    timing headroom  (dose at onset)    {Rt0-Rt[i_rx]:+.2f} log.d")
    print()
    if (Rt0 - Rt[i_rx]) > (ach_inf - ach_cal):
        print("  The timing headroom is the larger of the two by a wide margin, and")
        print("  it is the only one a better molecule cannot buy.  Making the drug")
        print("  perfect adds almost nothing, because the drug is already taking")
        print("  nearly all of a bound that delay has already shrunk.  This is the")
        print("  quantitative content of 'treat early': not that late treatment is")
        print("  weak, but that late treatment is being asked to act on a quantity")
        print("  that has mostly already been spent.")
    else:
        print("  In this parameterisation the potency headroom is the larger of the")
        print("  two -- reported as computed, against the expectation.")
    print()
    print("  Note the non-monotone TTAS in the timing sweep: dosing before the")
    print("  cytokine surge has begun can abort the illness outright, so the")
    print("  endpoint is not merely shortened, it is never entered.")
    return tim, pot


def A5_operator_decomposition():
    print("=" * 78)
    print("A5  OPERATOR DECOMPOSITION -- which term of the loop does each drug touch")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    print("  Each row switches ON exactly one operator at 95% efficacy from t_rx,")
    print("  with no drug PK, so the operators can be compared on equal terms.\n")
    E = 0.95
    ops = [
        ("PRODUCTION / RELEASE  (p)",     dict(_OP="p")),
        ("TRANSCRIPTION        (E->I)",   dict(_OP="k")),
        ("ENTRY                (beta)",   dict(_OP="b")),
        ("VIRION CLEARANCE     (c)",      dict(_OP="c")),
        ("INFECTED-CELL DEATH  (delta)",  dict(_OP="d")),
        ("TARGET PROTECTION    (T->R)",   dict(_OP="r")),
    ]
    print(f"  {'operator':32s} {'TTAS h':>8s} {'shed h':>8s} {'d24 log':>9s} "
          f"{'AUC log.d':>10s} {'epi lost':>9s}")
    base = None
    for label, spec in ops:
        over = {}
        op = spec["_OP"]
        # implement each operator as a time-gated parameter change via a wrapper
        t, out, P = _operator_run(op, E, t_rx)
        r = summarise(t, out, P, t_rx)
        if base is None:
            base = r
        print(f"  {label:32s} {_fmt(r['ttas_h'],1):>8s} {_fmt(r['shed_h'],1):>8s} "
              f"{_fmt(r['drop24'],2):>9s} {_fmt(r['auc_log'],2):>10s} "
              f"{r['epi_lost']*100:>8.1f}%")
    t, out, P = _operator_run(None, 0.0, t_rx)
    r = summarise(t, out, P, t_rx)
    print(f"  {'(none: placebo)':32s} {_fmt(r['ttas_h'],1):>8s} "
          f"{_fmt(r['shed_h'],1):>8s} {_fmt(r['drop24'],2):>9s} "
          f"{_fmt(r['auc_log'],2):>10s} {r['epi_lost']*100:>8.1f}%")
    print()
    print("  The signatures differ in KIND, not only in size.  Read the columns")
    print("  against each other rather than down:")
    print("   - ENTRY and TARGET PROTECTION barely move the 24-hour fall, because")
    print("     the cells already infected when the operator switches on keep")
    print("     producing; these operators only stop the NEXT generation.")
    print("   - PRODUCTION/RELEASE and VIRION CLEARANCE both cut the titre")
    print("     promptly and by almost the same amount -- unsurprising, since at")
    print("     quasi-steady state the titre is p*I/c and they are the numerator")
    print("     and denominator of the same ratio.")
    print("   - INFECTED-CELL DEATH is the strongest single operator in the model,")
    print("     because it is the only one that removes the source rather than")
    print("     throttling it.  No licensed influenza antiviral works this way;")
    print("     the immune system does, which is why the natural history is so")
    print("     much more decisive than any of the drugs acting on it.")
    print("   - TRANSCRIPTION (diverting eclipse cells to an abortive fate) is")
    print("     WEAKER at the same 95% than production blockade, which is not")
    print("     obvious: it only redirects cells that have yet to mature, while")
    print("     production blockade also silences every cell that already has.")
    print("     Baloxavir beats the neuraminidase inhibitors in the ledger")
    print("     despite acting at the weaker point, purely on achieved potency.")
    print()
    print("  This is why 'antiviral efficacy' is not a single number: two drugs")
    print("  with the same effect on viral AUC can have quite different effects")
    print("  on the 24-hour fall, on shedding duration, and on epithelial loss.")


def _operator_run(op, E, t_rx):
    """Integrate with one idealised operator switched on at t_rx."""
    par = dict(P0)
    doses = []

    def rhs_op(t, y, P):
        Q = dict(P)
        if op is not None and t >= t_rx:
            if op == "p":
                Q["PVIR"] = P["PVIR"] * (1 - E); Q["PVIRL"] = P["PVIRL"] * (1 - E)
            elif op == "k":
                Q["FMATU"] = P["FMATU"] * (1 - E)
            elif op == "b":
                Q["BETA"] = P["BETA"] * (1 - E); Q["BETAL"] = P["BETAL"] * (1 - E)
            elif op == "c":
                Q["CVIR"] = P["CVIR"] / max(1e-9, (1 - E))
            elif op == "d":
                Q["DELTA"] = P["DELTA"] / max(1e-9, (1 - E))
                Q["DELTAL"] = P["DELTAL"] / max(1e-9, (1 - E))
            elif op == "r":
                Q["PHIF"] = P["PHIF"] / max(1e-9, (1 - E))
        return rhs(t, y, Q)

    tgrid = np.arange(0.0, 21.0 + 1e-9, 1 / 96.)
    y = y0_of(par)
    breaks = [0.0, t_rx, 21.0]
    out = np.zeros((len(tgrid), NSTATE))
    filled = 0
    for a, b in zip(breaks[:-1], breaks[1:]):
        seg = tgrid[(tgrid >= a - 1e-12) & (tgrid <= b + 1e-12)]
        sol = solve_ivp(rhs_op, (a, b), y, args=(par,), method="LSODA",
                        t_eval=seg, rtol=1e-8, atol=1e-10, max_step=0.05)
        n = sol.y.shape[1]
        idx = np.searchsorted(tgrid, seg[0] - 1e-12)
        out[idx:idx + n, :] = sol.y.T
        filled = max(filled, idx + n)
        y = np.maximum(sol.y[:, -1].copy(), 0.0)
    if filled < len(tgrid):
        out[filled:, :] = out[filled - 1, :]
    return tgrid, out, par


def A6_competitive_release():
    print("=" * 78)
    print("A6  RESISTANCE IS COMPETITIVE RELEASE -- AND IT NEEDS A FIELD TO RELEASE")
    print("=" * 78)
    ons = placebo_onset()
    print("  Baloxavir, PA/I38T profile (EC50 x%.0f, replicative cost %.0f%%)."
          % (I38T["RF_BX"], I38T["COST"] * 100))
    print("  The mutant is NEVER seeded: it exists only because the mutation term")
    print("  (mu = %.1e per replication) acts on the wild type from cycle one.\n"
          % P0["MU"])
    print("  Peak mutant titre (log10 TCID50/mL) after dosing, over dose time x")
    print("  potency.  Read ACROSS for the potency effect, DOWN for the timing one:\n")
    hs = (0, 12, 24, 36, 48, 72)
    ems = (0.0, 0.90, 0.99, 0.999, 0.9999, 0.999999)
    print("      dose h  " + "".join(f"{('E=' + format(e, '.6g')):>11s}" for e in ems))
    grid = {}
    for h in hs:
        t_rx = ons + h / 24.0
        row = []
        for em in ems:
            t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=30.0,
                                 EMAX_BX=em, **I38T)
            vm = np.log10(np.maximum(col(out, "VM"), 1e-6))
            row.append(float(np.max(vm[t >= t_rx])))
            grid[(h, em)] = row[-1]
        print(f"      {h:>6d}  " + "".join(f"{x:>11.2f}" for x in row))
    print()
    print("  Target cells still available at each dose time (untreated run):")
    t0, o0, Pz = simulate(tmax=30.0)
    for h in hs:
        i = min(int(np.searchsorted(t0, ons + h / 24.0)), len(t0) - 1)
        print(f"      {h:>6d} h post-onset:  T = "
              f"{col(o0,'TU')[i]/Pz['T0U']*100:6.2f} % of T0")
    print()
    best = max(grid.items(), key=lambda kv: kv[1])
    print(f"  The maximum sits at dose time {best[0][0]} h post-onset, Emax = "
          f"{best[0][1]:g}, giving a mutant peak of {best[1]:.2f} log10.")
    print()
    print("  Two things fall out that were not put in:")
    print()
    print("  (1) ALONG THE POTENCY AXIS the mutant does not simply fall as the")
    print("      drug improves.  At Emax = 0 it loses to a fitter wild type; as")
    print("      the drug improves it stops losing, because its competitor is")
    print("      being removed and the target cells the wild type would have")
    print("      taken are handed to it instead.  Selection for resistance is a")
    print("      by-product of the drug WORKING, not of it failing.")
    print()
    print("  (2) ALONG THE TIMING AXIS there is a hard prerequisite that is easy")
    print("      to miss: competitive release needs something to release.  Dosing")
    print("      late selects far less -- not because the drug is weaker, but")
    print("      because the wild type has already consumed the epithelium and")
    print("      there is no spared field for the mutant to grow in.  Resistance")
    print("      emergence and clinical benefit therefore share a prerequisite,")
    print("      and the same policy -- treat earlier -- buys more of both.")
    print()
    print("  DISCREPANCY 4.  This is a deterministic model, so a mutant lineage at")
    print("  10^-4 of the population is 'present' in every simulated patient.")
    print("  Real I38T emergence is a stochastic establishment event in ~10% of")
    print("  treated adults.  The peaks above are expectations over an ensemble")
    print("  that never happens to any single patient, and should be read as a")
    print("  selection pressure, not as a predicted viral load.")
    return grid


def A7_children_vs_adults():
    print("=" * 78)
    print("A7  WHO SELECTS I38T -- AND WHERE THIS MODEL GETS IT WRONG")
    print("=" * 78)
    print("  CAPSTONE-1 reported PA/I38T/M/F in 9.7% of baloxavir-treated adults")
    print("  (Hayden 2018, PMID 30184455); paediatric trials have reported higher")
    print("  frequencies (miniSTONE-2, Baker 2020, PMID 32516282).  The model is")
    print("  given NO age term.  The only differences below are prior immunity and the")
    print("  size of the epithelial target pool.  Every arm is dosed 12 h after")
    print("  its own symptom onset, i.e. while a field still exists to release.\n")
    hosts = [
        ("adult, prior exposure",   dict(CTL0=6.0, IGA0=0.30, FAB=1.0,  T0U=4.0e8)),
        ("adult, naive",            dict(CTL0=1.0, IGA0=0.0,  FAB=1.0,  T0U=4.0e8)),
        ("child (naive, larger T)", dict(CTL0=1.0, IGA0=0.0,  FAB=0.65, T0U=6.0e8)),
        ("immunocompromised",       dict(CTL0=1.0, IGA0=0.0,  FAB=0.15, FCTL=0.12,
                                         FIFN=0.55, T0U=4.0e8)),
    ]
    print(f"  {'host':26s} {'shed h':>7s} {'mut peak':>9s} {'mut>LOD':>8s} "
          f"{'mut AUC':>8s} {'T at dose':>10s}")
    for label, over in hosts:
        ons = placebo_onset(**over)
        t_rx = ons + 0.5
        t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=35.0,
                             **dict(I38T, **over))
        post = t >= t_rx
        vm = np.log10(np.maximum(col(out, "VM"), 1e-6))
        dt = t[1] - t[0]
        above = float(np.sum((vm > LOD) & post) * dt)
        mauc = float(np.trapezoid(np.maximum(vm - LOD, 0.0) * post, t))
        i = min(int(np.searchsorted(t, t_rx)), len(t) - 1)
        print(f"  {label:26s} {_fmt(time_to_shed_stop(t,out,P,t_rx),1):>7s} "
              f"{np.max(vm[post]):>9.2f} {above:>7.2f}d {mauc:>8.2f} "
              f"{col(out,'TU')[i]/P['T0U']*100:>9.1f}%")
    print()
    print("  Read the mutant AUC column rather than the peak: what a resistant")
    print("  lineage needs in order to be sampled and reported is TIME above the")
    print("  assay floor, not a high maximum.")
    print()
    print("  DISCREPANCY 5 -- and it is the model's clearest failure.  The ordering")
    print("  comes out BACKWARDS.  The adult with pre-existing immunity selects")
    print("  the most mutant of the three immunocompetent hosts and the child the")
    print("  least, which is the reverse of the paediatric-versus-adult signal.")
    print()
    print("  The mechanism producing the reversal is visible in the last column,")
    print("  and it is not a coding accident: competitive release is governed")
    print("  almost entirely by T(t_rx), the field still standing when the drug")
    print("  arrives, and ANYTHING that slows the wild type leaves more of it.")
    print("  Pre-existing CD8 memory and mucosal IgA slow the wild type.  So in")
    print("  this structure prior immunity increases the field, and therefore")
    print("  increases the release.  The child's larger epithelium makes the wild")
    print("  type faster, so by 12 h post-onset the child has consumed more of its")
    print("  own field and has less left to hand over.")
    print()
    print("  Two honest readings, and the model cannot choose:")
    print("   - the structure is right and the clinical comparison is confounded")
    print("     (children are brought to care earlier in their illness, are dosed")
    print("     by weight band, and shed longer, all of which push the other way);")
    print("   - or target-cell availability is NOT the dominant control on I38T")
    print("     emergence, and a well-mixed epithelium is the wrong idealisation")
    print("     for a mucosa where the two strains are spatially separated.")
    print()
    print("  What the model does get right is the immunocompromised host, whose")
    print("  mutant AUC is the largest by a wide margin -- consistent with the")
    print("  prolonged shedding and repeated resistance emergence reported in that")
    print("  population.  There, both mechanisms point the same way.")


def A8_symptom_transfer():
    print("=" * 78)
    print("A8  WHAT THE SYMPTOM ENDPOINT REQUIRES OF THE MODEL")
    print("=" * 78)
    print("  The virological endpoints of CAPSTONE-1 are reproduced by viral")
    print("  dynamics alone.  The SYMPTOM endpoint is not, and the reason is")
    print("  structural rather than numerical.")
    print()
    print("  If the symptom score is driven purely by the cytokine state (WVIR = 0),")
    print("  its size and timing are settled by the interferon peak, which happens")
    print("  BEFORE any drug given at the median enrolment time can act.  Under")
    print("  that structure no antiviral, at any potency, shortens the illness.")
    print("  WVIR is the fraction of the symptom drive that instead tracks the")
    print("  instantaneous titre.  Sweeping it asks the trial a question:")
    print("  how much of how a patient feels has to be same-day virus?\n")
    ons_ref = placebo_onset()
    print(f"  {'WVIR':>6s} {'onset':>7s} {'placebo':>8s} {'oselt':>8s} {'balox':>8s} "
          f"{'oselt':>8s} {'balox':>8s}")
    print(f"  {'':>6s} {'h p.i.':>7s} {'TTAS h':>8s} {'TTAS h':>8s} {'TTAS h':>8s} "
          f"{'benefit':>8s} {'benefit':>8s}")
    for w in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        ons = placebo_onset(WVIR=w)
        t_rx = ons + RX_DELAY_H / 24.0
        r = {}
        for k, doses in (("pl", []), ("os", regimen_oseltamivir(t_rx)),
                         ("bx", regimen_baloxavir(t_rx))):
            t, out, P = simulate(doses=doses, tmax=30.0, WVIR=w)
            r[k] = summarise(t, out, P, t_rx)
        print(f"  {w:>6.2f} {ons*24:>7.1f} {r['pl']['ttas_h']:>8.1f} "
              f"{r['os']['ttas_h']:>8.1f} {r['bx']['ttas_h']:>8.1f} "
              f"{r['pl']['ttas_h']-r['os']['ttas_h']:>8.1f} "
              f"{r['pl']['ttas_h']-r['bx']['ttas_h']:>8.1f}")
    print()
    print("  CAPSTONE-1 measured a 26.5 h benefit for both drugs.  At WVIR = 0 the")
    print("  model gives essentially zero for either.  The published symptom")
    print("  endpoint therefore cannot be explained by a cytokine-only symptom")
    print("  model, and the value this model carries (WVIR = %.2f) is the price of"
          % P0["WVIR"])
    print("  matching it.  That is a falsifiable claim about influenza symptoms,")
    print("  not a fitting convenience: it says a substantial part of how a")
    print("  patient feels tracks today's replication rather than the size of the")
    print("  cytokine surge that has already passed.")
    print()
    print("  DISCREPANCY 3.  Even at the calibrated WVIR the model splits the two")
    print("  drugs, giving baloxavir a larger symptom benefit than oseltamivir,")
    print("  whereas CAPSTONE-1 found them indistinguishable (53.7 vs 53.8 h)")
    print("  despite a 2-log difference in day-2 titre.  A model in which symptoms")
    print("  track virus at all cannot reproduce two arms that differ 100-fold in")
    print("  virus and not at all in symptoms.  Either the symptom-virus coupling")
    print("  saturates well below the titres these drugs separate, or the")
    print("  patient-reported endpoint is not resolving what the assay is.")


def A9_steroid_harm():
    print("=" * 78)
    print("A9  CORTICOSTEROIDS: SHORT-TERM COMFORT PAID FOR IN VIRAL AUC")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    arms = [
        ("oseltamivir alone", regimen_oseltamivir(t_rx), {}),
        ("oseltamivir + dexamethasone 6 mg x5d",
         regimen_oseltamivir(t_rx) + regimen_dexamethasone(t_rx), {}),
        ("dexamethasone alone", regimen_dexamethasone(t_rx), {}),
    ]
    print(f"  {'arm':40s} {'AUC':>7s} {'shed h':>8s} {'peak sym':>9s} "
          f"{'TTAS h':>8s} {'bact':>6s}")
    for label, doses, over in arms:
        t, out, P = simulate(doses=doses, tmax=25.0, **over)
        r = summarise(t, out, P, t_rx)
        print(f"  {label:40s} {r['auc_log']:>7.2f} {_fmt(r['shed_h'],1):>8s} "
              f"{r['peak_sym']:>9.1f} {_fmt(r['ttas_h'],1):>8s} "
              f"{r['bac_peak']:>6.2f}")
    print()
    print("  The direction is right and the size is small, and both are worth")
    print("  saying.  Dexamethasone lowers the symptom score and shortens time to")
    print("  alleviation -- it acts directly on the cytokine arm that the score")
    print("  reads -- while raising viral AUC and lengthening shedding, because")
    print("  the same suppression removes interferon-driven target protection and")
    print("  CD8 killing.  That is the shape of the observational mortality signal")
    print("  (Lansbury 2019 Cochrane, PMID 30798570) recovered from mechanism")
    print("  alone, with no mortality term anywhere in the model.")
    print()
    print("  But the magnitude is a caveat on the model, not a reassurance about")
    print("  the drug.  In THIS parameterisation neither interferon-driven target")
    print("  protection nor CD8 killing is rate-limiting for clearance -- target")
    print("  exhaustion and antibody are -- so suppressing them costs little.  A")
    print("  model calibrated on a host in which the CD8 response IS rate-limiting")
    print("  (the elderly, the immunosuppressed) would show a much larger penalty,")
    print("  and those are exactly the patients the observational signal comes")
    print("  from.  The small effect here should be read as a property of the")
    print("  calibration, not as evidence that steroids are nearly harmless.")


def A10_prophylaxis():
    print("=" * 78)
    print("A10  POST-EXPOSURE PROPHYLAXIS IS THE SAME DRUG ON THE OTHER SIDE OF THE PEAK")
    print("=" * 78)
    print("  Dosing before the target pool is spent converts the same molecule from")
    print("  a symptom-duration drug into an infection-blocking one.\n")
    print(f"  {'dose time (h p.i.)':>20s} {'peak log10':>11s} {'shed d':>8s} "
          f"{'peak sym':>9s} {'epi lost':>9s} {'illness?':>9s}")
    for h in (-12, 0, 12, 24, 36, 48, 72):
        t_rx = max(0.0, h / 24.0)
        t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=25.0)
        v = vlog(out)
        r = summarise(t, out, P, 0.0)
        ill = "yes" if np.max(col(out, "SYM")) >= P["SONSET"] else "NO"
        print(f"  {h:>20d} {np.max(v):>11.2f} {r['shed_days']:>8.2f} "
              f"{r['peak_sym']:>9.1f} {r['epi_lost']*100:>8.1f}% {ill:>9s}")
    print()
    print("  BLOCKSTONE (Ikematsu 2020 NEJM 383:309, PMID 32640124) reported 1.9% vs")
    print("  13.6% symptomatic influenza with single-dose baloxavir prophylaxis --")
    print("  a 7x reduction in the probability of illness, not in its duration.")
    print("  In this model the switch between the two regimes is the sign of")
    print("  (peak time - dose time), which is why the same drug looks like two")
    print("  different drugs on either side of it.")


def A11_lrt_and_superinfection():
    print("=" * 78)
    print("A11  THE LOWER TRACT AND THE BACTERIAL SEQUEL")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    print("  Secondary bacterial pneumonia is modelled as the product of two")
    print("  influenza-created conditions: exposed adhesion sites on denuded")
    print("  epithelium, and interferon-mediated suppression of antibacterial")
    print("  defence (Shahangian 2009 J Clin Invest, PMID 19487810).\n")
    print(f"  {'arm':34s} {'LRT peak':>9s} {'LRT damage':>11s} {'bact peak':>10s} "
          f"{'>threshold':>11s}")
    arms = [
        ("placebo", [], {}),
        ("oseltamivir 24 h post-onset", regimen_oseltamivir(t_rx), {}),
        ("baloxavir 24 h post-onset", regimen_baloxavir(t_rx), {}),
        ("baloxavir at onset", regimen_baloxavir(ons), {}),
        ("baloxavir 72 h post-onset", regimen_baloxavir(ons + 3.0), {}),
        ("placebo, 5x LRT seeding", [], dict(ADESC=0.040)),
        ("placebo, blunted interferon", [], dict(FIFN=0.45)),
    ]
    for label, doses, over in arms:
        t, out, P = simulate(doses=doses, tmax=25.0, **over)
        bp = float(np.max(col(out, "BAC")))
        print(f"  {label:34s} {np.max(vlog(out,False)):>9.2f} "
              f"{np.max(col(out,'DL'))/P['T0L']*100:>10.1f}% {bp:>10.2f} "
              f"{('YES' if bp>=P['BACTHR'] else 'no'):>11s}")
    print()
    print("  Two things worth reading off this table rather than assuming.")
    print()
    print("  First, in the uncomplicated adult the bacterial burden climbs about")
    print("  four logs and still does not reach the superinfection threshold in")
    print("  ANY antiviral arm.  Timing barely moves it.  The only arm that")
    print("  changes the bacterial endpoint is the one that aborts the infection")
    print("  outright -- which is the honest version of 'antivirals prevent")
    print("  secondary pneumonia': in this model they do it by preventing the")
    print("  influenza, not by softening it.")
    print()
    print("  Second -- and this is the limitation, not a result -- NOTHING in the")
    print("  table crosses the threshold, including the two deliberately worsened")
    print("  arms.  The threshold itself (%.1f log10 CFU) has no independent"
          % P0["BACTHR"])
    print("  calibration; it is a stipulation.  So this module can ORDER arms by")
    print("  bacterial risk, which is what the ordering above is for, and it")
    print("  cannot say who develops pneumonia.  Any absolute risk read off it")
    print("  would be reading the stipulation back out.")


def A12_immunocompromised():
    print("=" * 78)
    print("A12  IN THE IMMUNOCOMPROMISED HOST THE WINDOW NEVER CLOSES")
    print("=" * 78)
    print("  In a normal host the drug competes with an immune response that would")
    print("  have cleared the virus anyway; its marginal value therefore decays with")
    print("  time.  Remove the immune response and the marginal value stops decaying.\n")
    print(f"  {'host':22s} {'dose h post-onset':>18s} {'shed d':>8s} "
          f"{'AUC log.d':>10s} {'benefit':>9s}")
    for label, over in [("immunocompetent", {}),
                        ("immunocompromised", dict(FCTL=0.12, FAB=0.15, FIFN=0.55))]:
        ons = placebo_onset(**over)
        t0, o0, P0_ = simulate(doses=[], tmax=35.0, **over)
        base = summarise(t0, o0, P0_, ons)["auc_log"]
        for h in (0, 24, 48, 96, 168):
            t_rx = ons + h / 24.0
            t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=35.0,
                                 **over)
            r = summarise(t, out, P, ons)
            print(f"  {label:22s} {h:>18d} {r['shed_days']:>8.2f} "
                  f"{r['auc_log']:>10.2f} {base - r['auc_log']:>8.2f}L")
    print()
    print("  DISCREPANCY 2.  The model says late antiviral therapy retains most of")
    print("  its value in the immunocompromised host.  That is consistent with the")
    print("  clinical practice of treating these patients regardless of symptom")
    print("  duration, but no randomised trial has ever tested it, so the prediction")
    print("  is unfalsified rather than validated.")


def A13_sensitivity():
    print("=" * 78)
    print("A13  LOCAL SENSITIVITY OF THE TWO HEADLINE OUTPUTS")
    print("=" * 78)
    ons = placebo_onset()
    t_rx = ons + RX_DELAY_H / 24.0
    keys = ["BETA", "PVIR", "DELTA", "CVIR", "KECL", "QF", "PHIF", "KKILL",
            "RCTL", "T0U", "LREG", "EC50_BX", "COST", "MU", "K6", "KOFF"]
    t, out, P = simulate(doses=regimen_baloxavir(t_rx), tmax=21.0, **I38T)
    b = summarise(t, out, P, t_rx)
    print(f"  {'parameter':10s} {'+20% -> TTAS':>14s} {'+20% -> AUC':>13s} "
          f"{'+20% -> mut frac':>18s}")
    print(f"  {'':10s} {'(rel. sens.)':>14s} {'(rel. sens.)':>13s} {'(rel. sens.)':>18s}")
    for k in keys:
        over = dict(I38T)
        over[k] = P0[k] * 1.2
        t2, o2, P2 = simulate(doses=regimen_baloxavir(t_rx), tmax=21.0, **over)
        r = summarise(t2, o2, P2, t_rx)
        def s(new, old):
            if old is None or old == 0 or math.isnan(old) or math.isnan(new):
                return float("nan")
            return (new - old) / old / 0.2
        print(f"  {k:10s} {_fmt(s(r['ttas_h'],b['ttas_h']),3):>14s} "
              f"{_fmt(s(r['auc_log'],b['auc_log']),3):>13s} "
              f"{_fmt(s(r['mut_frac_max'],b['mut_frac_max']),3):>18s}")
    print()
    print("  Relative sensitivity = (dY/Y)/(dX/X).  Values above 1 in magnitude mean")
    print("  the output is amplifying the parameter's uncertainty.")


ANALYSES = {
    "A0":  ("untreated natural history / calibration", A0_baseline),
    "A1":  ("the exact upper bound on antiviral benefit", A1_residual_bound),
    "A2":  ("PD calibration and the in-vitro gap", A2_pd_calibration),
    "A3":  ("scenario ledger vs published trials", A3_trial_ledger),
    "A4":  ("timing vs potency", A4_timing_vs_potency),
    "A5":  ("operator decomposition", A5_operator_decomposition),
    "A6":  ("resistance as competitive release", A6_competitive_release),
    "A7":  ("children vs adults", A7_children_vs_adults),
    "A8":  ("virus-to-symptom transfer", A8_symptom_transfer),
    "A9":  ("corticosteroid harm", A9_steroid_harm),
    "A10": ("post-exposure prophylaxis", A10_prophylaxis),
    "A11": ("lower tract and bacterial superinfection", A11_lrt_and_superinfection),
    "A12": ("the immunocompromised host", A12_immunocompromised),
    "A13": ("local sensitivity", A13_sensitivity),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", help="run only these analyses (e.g. A4 A6)")
    ap.add_argument("--list", action="store_true", help="list analyses and exit")
    a = ap.parse_args()
    if a.list:
        for k, (desc, _) in ANALYSES.items():
            print(f"  {k:5s} {desc}")
        return 0
    todo = a.only if a.only else list(ANALYSES)
    for k in todo:
        key = k.upper()
        if key not in ANALYSES:
            print(f"unknown analysis {k}", file=sys.stderr)
            return 2
        ANALYSES[key][1]()
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
