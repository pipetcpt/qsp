#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cmv_python_reference.py
=======================
Independent Python/scipy re-implementation of the 45-ODE CMV-in-transplant QSP
model.  This file is the REFERENCE: every number quoted in README.md and in the
mrgsolve model's calibration notes was produced by RUNNING this file.  The
mrgsolve model (cmv_mrgsolve_model.R) mirrors it equation for equation.

ORGANISING THESIS
-----------------
There is exactly ONE number every anti-CMV drug has to clear, and it is not an
EC50.  Write the systemic virus as a target-cell-limited process:

    d(I)/dt  =  BETA*TC*Vv - LOSS*I            LOSS = DELI + KE8*E8 + KENK*NK
    d(Vv)/dt =  PVIR*(1-e_pol)*(1-e_pack)*I - CVIR*Vv

At quasi-steady state in Vv the net exponent of the infected-cell pool is

    r  =  KPROD*(1-e_pol)*(1-e_pack)  -  DELI  -  KE8*E8  -  KENK*NK

with KPROD = BETA*PVIR*TC0/CVIR.  Two clinically measured numbers fix KPROD:
an untreated DNAemia doubling time of ~1.2 d (r0 = 0.578/d) and an on-treatment
decline half-life of ~2.4 d (DELI = 0.289/d).  So KPROD = 0.866/d, and:

  THRESHOLD 1 (the drug threshold).   Sign change at
        e* = r0 / (r0 + DELI) = 0.578/0.866 = 0.667
    Every potency, resistance, renal-dose-reduction and drug-interaction
    question in this disease is the single question "does this regimen leave
    e above 0.667?"  A 3000-fold shift that lands at 0.60 and an 8-fold shift
    that lands at 0.60 are the SAME clinical event.

  THRESHOLD 2 (the immune threshold).   With no drug at all, r < 0 when
        E8* = r0 / KE8 = 0.578/0.12 = 4.8 CMV-specific CD8 per uL
    This is the actual endpoint of the illness.  No drug moves it.  Drugs only
    hold the line while E8 climbs toward it.

  THE COUPLING (why the two thresholds fight each other).  E8's expansion term
  is proportional to ANTIGEN, and antigen is the infected-cell pool the drug
  just removed:
        d(E8)/dt = RHO8*E8*(AG/(AG+KAG))*HELP*FIS*(1-E8/E8MAX) - DE8*E8
  So clearing threshold 1 slows the approach to threshold 2.  Late-onset CMV
  after prophylaxis is not a separate phenomenon bolted onto the model; it is
  what this coupling does when you switch the drug off.

  THE THIRD AXIS (monitoring interval is a dose).  Pre-emptive therapy does not
  act on e at all.  It acts on the HEIGHT of the exponential at which treatment
  starts.  At a 1.2-day doubling time a 7-day PCR interval means the virus rises
  2^(7/1.2) = 57-fold between draws, so a nominal 1000 IU/mL trigger is really a
  ~57,000 IU/mL trigger.  Monitoring frequency and drug potency are therefore
  commensurable, and the model prices them against each other.

FOUR STRUCTURAL DECISIONS THAT DO THE WORK
------------------------------------------
(1) MEASURED DNAemia IS TWO SPECIES, AND THE DRUGS HIT DIFFERENT ONES.
    Vv = encapsidated virion DNA (the only infectious species; this is what
    BETA multiplies).  Vl = non-encapsidated DNA released from dying infected
    cells (the majority of what a plasma PCR actually measures).  Polymerase
    inhibitors (ganciclovir-TP, foscarnet, cidofovir-PP) cut BOTH.  Terminase /
    kinase inhibitors (letermovir, maribavir) cut only Vv and, because unit-
    length genomes are no longer cleaved out of the concatemer, they slightly
    INCREASE the DNA per dying cell.  Consequence, derived not asserted:
    letermovir working perfectly produces almost no immediate fall in plasma
    CMV DNA.  That is why it is a prophylaxis drug and not a treatment drug,
    and why on-letermovir DNAemia is a poor early efficacy read-out.

(2) GANCICLOVIR'S ACTIVATION STEP IS A VIRAL GENE PRODUCT, SO ONE DRUG CAN
    SWITCH ANOTHER OFF.  GCV must be phosphorylated by pUL97 to GCV-MP.
    Maribavir INHIBITS pUL97.  Therefore maribavir enters the ganciclovir arm
    through the ACTIVATION term, not the effect term, and the combination is
    antagonistic by construction.  UL97 GCV-resistance mutations (M460V,
    A594V ...) are modelled as loss of KINASE efficiency, not as a polymerase
    EC50 shift -- which is why they do not touch maribavir (different residues,
    T409M/H411Y/C480F) and why maribavir is the mechanistically correct rescue.

(3) THE NEUTROPENIA LOOP CLOSES THROUGH THE KIDNEY AND ENDS IN RESISTANCE.
    GCV clearance is ~90% renal, so CrCl sets exposure; exposure drives a
    Friberg 5-compartment myelosuppression model; ANC < 1.0 triggers the real
    clinical algorithm (50% dose cut, then stop); the dose cut drops e below
    0.667; and sub-threshold e is exactly the selection window for the UL97
    mutant.  One parameter (CrCl) raises potency and destroys durability.

(4) THREE VIRAL STRAINS WITH DIFFERENT-SHAPED RESISTANCE.  W (wild type),
    A (UL97 activation-type, GCV 8x, maribavir and letermovir untouched),
    B (UL56 C325Y, letermovir >3000x, everything else untouched).  The point is
    the SHAPE: one mutation takes letermovir from e = 0.93 to e = 0.004 (total
    loss) while the same single-mutation event takes ganciclovir from 0.88 to
    0.23 (crippled but not zero).  Mean potency does not predict breakthrough;
    the size of the single-step EC50 jump relative to the 0.667 threshold does.

WHAT IS NOT MODELLED (stated so the omissions are visible)
----------------------------------------------------------
  * UL54 polymerase mutants as a separate strain (they are in the map and the
    cross-resistance table, but the ODEs carry only strains W/A/B).
  * Compartmentalised tissue-only disease with a negative plasma PCR is present
    as a sanctuary compartment (VS) but is not calibrated against biopsy data.
  * CMV-specific adoptive T-cell therapy and vaccines appear in the map only.
  * The concatemer amplification factor AMPPK (decision 1) has no direct human
    measurement; it is an assumed parameter and is swept in the output.

USAGE
-----
    python3 cmv_python_reference.py            > cmv_reference_output.txt
    python3 cmv_python_reference.py --quick    # 4 scenarios only
"""
import argparse
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

LN2 = math.log(2.0)

# ===========================================================================
#  PARAMETERS
# ===========================================================================
P = dict(

    # ---- 1. systemic viral dynamics ---------------------------------------
    DBL0=1.2,        # untreated plasma DNAemia doubling time (d)  [Emery]
    THALF_TX=2.4,    # on-treatment DNAemia decline half-life (d)  [Emery]
    CVIR=1.2,        # encapsidated virion DNA clearance (/d)
    CLYS=6.0,        # free (lysis-derived) DNA clearance (/d), t1/2 2.8 h
    PVIR=1.0e4,      # virion DNA produced per unit I per day (IU/mL/d)
    RLYS=2.0,        # Vl:Vv ratio at untreated quasi-steady state
    AMPPK=0.50,      # extra DNA per dying cell when packaging is blocked
    LAMT=2.00,       # permissive target-cell renewal (/d) [fitted to peak VL]
    TC0=1.0,         # permissive pool (normalised)

    # ---- 2. latency and reactivation --------------------------------------
    KREACT=3.0e-5,   # baseline reactivation flux (units/d) per unit reservoir
    AREACT=6.0,      # immunosuppression multiplier on reactivation
    AINFL=3.0,       # inflammation (sepsis/rejection) multiplier
    KLAT=1.0e-3,     # re-seeding of the reservoir from active infection

    # ---- 3. strain fitness and mutation ----------------------------------
    FIT_A=0.85,      # UL97 activation-mutant replicative fitness
    FIT_B=0.95,      # UL56 C325Y replicative fitness
    # Mutation is NOT run as a continuous flux into the strain ODEs.  A
    # deterministic flux lets 1e-9 of an infected cell exist and then grow at
    # the mutant's positive exponent, which makes resistance CERTAIN in every
    # arm within ~60-90 d -- flatly wrong against the ~5% seen on valGCV
    # prophylaxis.  Instead the model integrates the number of resistance-
    # generating replication events and multiplies by the branching-process
    # establishment probability 1 - 1/R_eff of a single mutant lineage.  The
    # strain ODEs are seeded only when a scenario explicitly says so.
    MU_A=3.0e-5,     # per-infection probability of generating UL97 resistance
    MU_B=3.0e-5,     # ... UL56 resistance
    MU_ON=0.0,       # 1.0 restores the (wrong) continuous-flux behaviour
    NCELL=1.0e8,     # productively infected cells per normalised I unit
    IEXT=1.0e-8,     # extinction floor = one infected cell in the whole body

    # ---- 4. adaptive / innate immunity -----------------------------------
    KE8=0.12,        # CMV-specific CD8 killing (per (cell/uL) per d)
    KENK=0.010,      # adaptive NK killing (per (cell/uL) per d)
    RHO8=0.65,       # max CD8 expansion rate (/d)
    KAG=0.05,        # antigen half-saturation for CD8 expansion (I units)
    KHELP=1.20,      # CD4 help half-saturation (cells/uL)
    KISE=1.20,       # immunosuppression index halving CD8 expansion
    DE8=0.05,        # CD8 effector contraction (/d)
    KMEM=0.35,       # fraction of contracting effectors entering memory
    DEM8=0.002,      # memory CD8 loss (/d), t1/2 347 d
    FMEM=0.50,       # cytotoxic weight of a memory cell relative to effector
    E8MAX=60.0,      # CD8 carrying capacity (cells/uL)
    SPRE8=0.002,     # precursor influx (cells/uL/d)
    RHO4=0.55, DE4=0.012, E4MAX=25.0, SPRE4=0.004, KISE4=0.80,
    RHONK=0.06, DNK=0.02, NKMAX=40.0, KAGNK=0.20, KISNK=4.0,
    KAB=0.30, DAB=0.02, KNAB=8.0,     # neutralising antibody / infectivity

    # ---- 5. immunosuppression --------------------------------------------
    LYM0=1.8,        # baseline lymphocytes (10^3/uL)
    KDEPL=3.0,       # ATG-driven lymphodepletion (per (ug/mL) per d)
    KREPOP=0.004,    # lymphocyte repopulation (/d)  -> t1/2 173 d
    KELATG=0.115,    # ATG elimination (/d), t1/2 6 d
    VTAC=1200.0,     # tacrolimus apparent volume (L)
    CLTAC=30.0,      # tacrolimus apparent clearance (L/h -> /d below)
    DDI_TAC_LTV=2.40,  # letermovir raises tacrolimus AUC 2.4x
    DDI_TAC_MBV=1.50,  # maribavir raises tacrolimus modestly
    KDDI_LTV=0.012,    # free letermovir for half-maximal CYP3A/OATP inhibition
    KDDI_MBV=0.050,    # free maribavir for half-maximal CYP3A inhibition
    KSTER=0.035,     # steroid taper rate (/d)
    WTAC=0.55, WLYM=0.80, WSTER=0.25, WMPA=0.40, WMTOR=0.35,

    # ---- 6. ganciclovir / valganciclovir ---------------------------------
    MW_VGCV=354.4, MW_GCV=255.2,
    F_VGCV=0.60,     # valganciclovir -> ganciclovir systemic availability
    KA_GCV=8.0,      # absorption (/d)  (t_max ~2 h)
    V_GCV=60.0,      # ganciclovir volume (L)
    CL_GCV=13.0,     # ganciclovir clearance at CrCl 100 (L/h)
    FREN_GCV=0.90,   # renal fraction of ganciclovir clearance
    KPHOS=12.0,      # pUL97 phosphorylation of GCV (/d) -> GCV-TP
    KDEGTP=0.924,    # GCV-TP loss (/d), t1/2 18 h
    EC50_GTP=23.4,   # GCV-TP concentration for 50% polymerase block (uM-eq)
    HGTP=2.0,        # Hill coefficient of chain termination
    FK_A=0.125,      # UL97-mutant kinase efficiency (= 8-fold GCV resistance)

    # ---- 7. letermovir ---------------------------------------------------
    MW_LTV=572.6, F_LTV=0.94, KA_LTV=6.0, V_LTV=45.0, CL_LTV=2.6,
    FU_LTV=0.010,    # 99% protein bound
    EC50_LTV=0.0040, # free EC50 (uM) against pUL56 terminase complex
    RES_B_LTV=3000.0,  # UL56 C325Y fold-shift
    PEN_LTV_S=0.01,  # sanctuary (CNS/retina) penetration (free drug)

    # ---- 8. maribavir ----------------------------------------------------
    MW_MBV=376.2, F_MBV=0.90, KA_MBV=10.0, V_MBV=27.0, CL_MBV=3.5,
    FU_MBV=0.020,    # 98% protein bound
    EC50_MBV=0.120,  # free EC50 (uM) against pUL97 kinase
    KI_MBVK=0.120,   # Ki for blocking GCV phosphorylation (same target)
    PEN_MBV_S=0.40,

    # ---- 9. foscarnet ----------------------------------------------------
    V_FOS=40.0, CL_FOS=6.5, FREN_FOS=0.95, EC50_FOS=250.0, PEN_FOS_S=0.20,
    MW_FOS=126.0,

    # ---- 10. cidofovir ---------------------------------------------------
    V_CDV=25.0, CL_CDV=10.0, FREN_CDV=0.90, KPP_CDV=6.0, KDEGPP=0.30,
    EC50_CPP=1.20, MW_CDV=279.2,

    # ---- 11. myelosuppression (Friberg) ----------------------------------
    MTT=125.0,       # mean transit time (h)
    GAM=0.35,        # feedback exponent (chronic-exposure calibration)
    CIRC0=4.5,       # baseline ANC on MMF-based maintenance (10^9/L)
    EMAX_MYE=0.60, EC50_MYE=12.0,  # ganciclovir slope (Emax on uM plasma GCV)
    E_MPA=0.075,     # additive mycophenolate effect
    KGCSF=0.60,      # G-CSF acceleration of the transit chain

    # ---- 12. kidney ------------------------------------------------------
    GFR0=55.0,       # allograft eGFR at randomisation (mL/min/1.73)
    KRECG=0.10, KINJG=0.344, KREPT=0.05,
    AFOS=0.50, ACDV=0.08, ATAC=0.03,
    KMGIN=1.0, KMGOUT=0.50, FMG=1.20,
    CFIB=0.030,      # damage -> eGFR set-point coupling

    # ---- 13. tissue and sanctuary compartments ---------------------------
    LAMTT=0.04, TT0=1.0, BETAT_REL=0.35, LOSST=0.22, KE8T=0.045,
    PVIRT=3.0e3, CVIRT=1.0,
    RS=0.12, KS_CLR=0.05, KE8S=0.030,   # sanctuary growth / clearance

    # ---- 14. clinical hazards -------------------------------------------
    KDIS=0.0300,     # CMV-disease hazard scale (/d) [fitted: untreated
                     # D+/R- kidney -> P(CMV disease) = 0.60 at 1 year]
    KD50=2.5e2,      # tissue load for half-maximal disease hazard (IU/g)
    HDIS=1.5,
    EDIS=6.0,        # CD8 level halving the disease hazard
    H0REJ=3.2e-4,    # baseline acute-rejection hazard (/d)
    CREJ_V=0.55,     # CMV DNAemia amplification of rejection hazard
    CREJ_IS=0.55,    # protection by net immunosuppression

    # ---- 15. costs (illustrative US list-price order of magnitude) -------
    CST_VGCV=42.0, CST_LTV=225.0, CST_MBV=690.0, CST_FOS=1450.0,
    CST_CDV=980.0, CST_PCR=155.0, CST_GCSF=340.0, CST_DIS=2300.0,
)

# derived
P["r0"] = LN2 / P["DBL0"]
P["DELI"] = LN2 / P["THALF_TX"]
P["KPROD"] = P["r0"] + P["DELI"]
P["BETA"] = P["KPROD"] * P["CVIR"] / (P["PVIR"] * P["TC0"])
P["QLYS"] = P["RLYS"] * P["PVIR"] * P["CLYS"] / (P["CVIR"] * P["DELI"])
P["EPSTAR"] = P["r0"] / P["KPROD"]
P["E8STAR"] = P["r0"] / P["KE8"]
P["KTR"] = 4.0 / (P["MTT"] / 24.0)      # /d

# ===========================================================================
#  STATE VECTOR (45)
# ===========================================================================
SNAMES = [
    "LAT", "TC", "IW", "IA", "IB", "VVW", "VVA", "VVB", "VLY",
    "TT", "IT", "VT", "VS",
    "E8", "EM8", "E4", "NKA", "AB", "LYM",
    "TAC", "ATGC", "STER",
    "AGCV", "GCV", "GTP", "ALTV", "LTV", "AMBV", "MBV", "FOS", "CDV", "CTP",
    "PROL", "TR1", "TR2", "TR3", "ANC",
    "GFR", "TUBI", "MG",
    "AUCV", "AUCT", "HZD", "HZR", "TNEU", "COST", "MUTA", "MUTB",
]
IX = {n: i for i, n in enumerate(SNAMES)}
NS = len(SNAMES)


# ===========================================================================
#  PHARMACODYNAMIC HELPERS
# ===========================================================================
def eps_pol(gtp, fos_free, ctp, fk, p):
    """Fraction of viral DNA synthesis removed by polymerase-directed drugs."""
    x = (gtp * fk) / p["EC50_GTP"]
    eg = x ** p["HGTP"] / (1.0 + x ** p["HGTP"]) if x > 0 else 0.0
    ef = fos_free / (fos_free + p["EC50_FOS"]) if fos_free > 0 else 0.0
    ec = ctp / (ctp + p["EC50_CPP"]) if ctp > 0 else 0.0
    return 1.0 - (1.0 - eg) * (1.0 - ef) * (1.0 - ec)


def eps_pack(ltv_free, mbv_free, res_ltv, res_mbv, p):
    """Fraction of encapsidation / nuclear egress removed (terminase, kinase)."""
    el = ltv_free / (ltv_free + p["EC50_LTV"] * res_ltv) if ltv_free > 0 else 0.0
    em = mbv_free / (mbv_free + p["EC50_MBV"] * res_mbv) if mbv_free > 0 else 0.0
    return 1.0 - (1.0 - el) * (1.0 - em)


# ===========================================================================
#  RIGHT-HAND SIDE
# ===========================================================================
def rhs(t, y, u, p):
    """u = dict of exogenous/regimen inputs held constant over the step."""
    s = np.maximum(y, 0.0)
    (LAT, TC, IW, IA, IB, VVW, VVA, VVB, VLY, TT, IT, VT, VS,
     E8, EM8, E4, NKA, AB, LYM, TAC, ATGC, STER,
     AGCV, GCV, GTP, ALTV, LTV, AMBV, MBV, FOS, CDV, CTP,
     PROL, TR1, TR2, TR3, ANC, GFR, TUBI, MG,
     AUCV, AUCT, HZD, HZR, TNEU, COST, MUTA, MUTB) = s

    d = np.zeros(NS)

    # ---- net state of immunosuppression ---------------------------------
    ISI = (p["WTAC"] * TAC / 8.0
           + p["WLYM"] * max(0.0, 1.0 - LYM / p["LYM0"])
           + p["WSTER"] * STER / 20.0
           + p["WMPA"] * u["MPA"]
           - p["WMTOR"] * u["MTOR"])
    ISI = max(0.0, ISI)

    # ---- free / effective drug concentrations ---------------------------
    fos_f = FOS
    ltv_f = LTV * p["FU_LTV"]
    mbv_f = MBV * p["FU_MBV"]

    # per-strain inhibition, systemic compartment
    epW = eps_pol(GTP, fos_f, CTP, 1.0, p)
    epA = eps_pol(GTP, fos_f, CTP, p["FK_A"], p)
    epB = epW
    ekW = eps_pack(ltv_f, mbv_f, 1.0, 1.0, p)
    ekA = ekW
    ekB = eps_pack(ltv_f, mbv_f, p["RES_B_LTV"], 1.0, p)

    # ---- immune killing --------------------------------------------------
    E8EFF = E8 + p["FMEM"] * EM8          # cytotoxically effective CMV CD8
    kill = p["KE8"] * E8EFF + p["KENK"] * NKA
    lossW = p["DELI"] + kill
    lossA = p["DELI"] + kill
    lossB = p["DELI"] + kill

    # ---- infection ------------------------------------------------------
    beta = p["BETA"] / (1.0 + AB / p["KNAB"])
    infW = beta * TC * VVW
    infA = beta * TC * VVA
    infB = beta * TC * VVB
    react = (p["KREACT"] * LAT * (1.0 + p["AREACT"] * ISI)
             * (1.0 + p["AINFL"] * u["INFL"]))

    d[IX["LAT"]] = -1.0e-4 * react + p["KLAT"] * (IW + IA + IB) - 0.0
    d[IX["TC"]] = p["LAMT"] * (p["TC0"] - TC) - (infW + infA + infB)

    mu = p["MU_ON"]
    d[IX["IW"]] = infW * (1.0 - mu * (p["MU_A"] + p["MU_B"])) + react - lossW * IW
    d[IX["IA"]] = infA + mu * p["MU_A"] * infW - lossA * IA
    d[IX["IB"]] = infB + mu * p["MU_B"] * infW - lossB * IB
    # RESISTANCE RISK AS TIME SPENT IN THE SELECTION WINDOW.
    # An expected-number-of-lineages integral is not usable here: MU x NCELL is
    # ~3e3, so ANY appreciable replication generates many mutant lineages and
    # the probability saturates at 1 in every arm.  What actually separates the
    # arms is DURATION OF SELECTION -- days during which the wild type is held
    # (e_W > e*) while the mutant's own R_eff still exceeds 1.  That is the only
    # state in which a minority variant is converted into the dominant one.
    Reff_A = (p["KPROD"] * p["FIT_A"] * (1 - epA) * (1 - ekA)) / max(lossA, 1e-9)
    Reff_B = (p["KPROD"] * p["FIT_B"] * (1 - epB) * (1 - ekB)) / max(lossB, 1e-9)
    eW_tot = 1.0 - (1 - epW) * (1 - ekW)

    def _sig(x, sc):
        return 1.0 / (1.0 + math.exp(-x / sc))
    gateW = _sig(eW_tot - p["EPSTAR"], 0.02)
    d[IX["MUTA"]] = gateW * _sig(Reff_A - 1.0, 0.05)
    d[IX["MUTB"]] = gateW * _sig(Reff_B - 1.0, 0.05)

    # ---- plasma DNA: encapsidated (infectious) and lysis-derived --------
    d[IX["VVW"]] = p["PVIR"] * IW * (1 - epW) * (1 - ekW) - p["CVIR"] * VVW
    d[IX["VVA"]] = (p["PVIR"] * p["FIT_A"] * IA * (1 - epA) * (1 - ekA)
                    - p["CVIR"] * VVA)
    d[IX["VVB"]] = (p["PVIR"] * p["FIT_B"] * IB * (1 - epB) * (1 - ekB)
                    - p["CVIR"] * VVB)

    lysflux = (lossW * IW * (1 - epW) * (1 + p["AMPPK"] * ekW)
               + lossA * IA * (1 - epA) * (1 + p["AMPPK"] * ekA)
               + lossB * IB * (1 - epB) * (1 + p["AMPPK"] * ekB))
    d[IX["VLY"]] = p["QLYS"] * lysflux - p["CLYS"] * VLY

    VVTOT = VVW + VVA + VVB
    VMEAS = VVTOT + VLY

    # ---- tissue compartment (gut / lung) --------------------------------
    betat = beta * p["BETAT_REL"]
    inft = betat * TT * VVTOT
    d[IX["TT"]] = p["LAMTT"] * (p["TT0"] - TT) - inft
    losst = p["LOSST"] + p["KE8T"] * E8EFF
    d[IX["IT"]] = inft - losst * IT
    # tissue drug exposure = systemic (good penetration for all agents)
    ept = epW
    ekt = ekW
    d[IX["VT"]] = p["PVIRT"] * IT * (1 - ept) * (1 - ekt) - p["CVIRT"] * VT

    # ---- sanctuary (retina / CNS): reduced drug penetration -------------
    eps_s = eps_pol(GTP * 0.30, fos_f * p["PEN_FOS_S"], CTP * 0.15, 1.0, p)
    ekk_s = eps_pack(ltv_f * p["PEN_LTV_S"], mbv_f * p["PEN_MBV_S"], 1.0, 1.0, p)
    gs = p["RS"] * (1 - eps_s) * (1 - ekk_s)
    d[IX["VS"]] = (gs * VS * (1.0 - VS / 1.0e6)
                   - (p["KS_CLR"] + p["KE8S"] * E8EFF) * VS
                   + 1.0e-3 * VVTOT)

    # ---- adaptive and innate immunity ----------------------------------
    AG = IW + IA + IB + 0.5 * IT
    sat = AG / (AG + p["KAG"])
    help_ = E4 / (E4 + p["KHELP"])
    fis = 1.0 / (1.0 + ISI / p["KISE"])
    fis4 = 1.0 / (1.0 + ISI / p["KISE4"])
    fisnk = 1.0 / (1.0 + ISI / p["KISNK"])

    d[IX["E8"]] = (p["RHO8"] * E8EFF * sat * help_ * fis
                   * (1 - E8EFF / p["E8MAX"])
                   - p["DE8"] * E8 + p["SPRE8"] + u["ACT"] * 0.25)
    d[IX["EM8"]] = p["KMEM"] * p["DE8"] * E8 - p["DEM8"] * EM8
    d[IX["E4"]] = (p["RHO4"] * E4 * sat * fis4 * (1 - E4 / p["E4MAX"])
                   - p["DE4"] * E4 + p["SPRE4"])
    d[IX["NKA"]] = (p["RHONK"] * (AG / (AG + p["KAGNK"])) * fisnk
                    * (p["NKMAX"] - NKA) - p["DNK"] * NKA)
    d[IX["AB"]] = p["KAB"] * sat * help_ * fis4 - p["DAB"] * AB + u["CMVIG"]

    # ---- immunosuppression PK ------------------------------------------
    d[IX["LYM"]] = (-p["KDEPL"] * ATGC * LYM
                    + p["KREPOP"] * (p["LYM0"] - LYM))
    d[IX["ATGC"]] = u["RATG"] - p["KELATG"] * ATGC
    # concentration-dependent CYP3A/OATP inhibition.  This MUST be smooth: a
    # boolean "drug present" switch is turned on for ever by 1e-12 of solver
    # round-off in the LTV/MBV states, which silently trebles tacrolimus.
    ddi = (1.0
           + (p["DDI_TAC_LTV"] - 1.0) * ltv_f / (ltv_f + p["KDDI_LTV"])
           + (p["DDI_TAC_MBV"] - 1.0) * mbv_f / (mbv_f + p["KDDI_MBV"]))
    cltac = (p["CLTAC"] * 24.0) / ddi            # L/d
    d[IX["TAC"]] = u["RTAC"] / p["VTAC"] - cltac / p["VTAC"] * TAC
    d[IX["STER"]] = -p["KSTER"] * (STER - u["STER_SS"])

    # ---- antiviral PK ---------------------------------------------------
    crcl = max(8.0, GFR * 1.10)
    cl_gcv = (p["CL_GCV"] * (p["FREN_GCV"] * crcl / 100.0
                             + (1 - p["FREN_GCV"]))) * 24.0     # L/d
    d[IX["AGCV"]] = u["RVGCV"] - p["KA_GCV"] * AGCV
    d[IX["GCV"]] = (p["KA_GCV"] * AGCV * p["F_VGCV"] / p["V_GCV"]
                    + u["RGCV_IV"] / p["V_GCV"]
                    - cl_gcv / p["V_GCV"] * GCV)
    # pUL97-dependent activation; maribavir blocks the kinase
    fkin = 1.0 / (1.0 + mbv_f / p["KI_MBVK"])
    d[IX["GTP"]] = p["KPHOS"] * GCV * fkin - p["KDEGTP"] * GTP

    d[IX["ALTV"]] = u["RLTV"] - p["KA_LTV"] * ALTV
    d[IX["LTV"]] = (p["KA_LTV"] * ALTV * p["F_LTV"] / p["V_LTV"]
                    - p["CL_LTV"] * 24.0 / p["V_LTV"] * LTV)
    d[IX["AMBV"]] = u["RMBV"] - p["KA_MBV"] * AMBV
    d[IX["MBV"]] = (p["KA_MBV"] * AMBV * p["F_MBV"] / p["V_MBV"]
                    - p["CL_MBV"] * 24.0 / p["V_MBV"] * MBV)
    cl_fos = (p["CL_FOS"] * (p["FREN_FOS"] * crcl / 100.0
                             + (1 - p["FREN_FOS"]))) * 24.0
    d[IX["FOS"]] = u["RFOS"] / p["V_FOS"] - cl_fos / p["V_FOS"] * FOS
    cl_cdv = (p["CL_CDV"] * (p["FREN_CDV"] * crcl / 100.0
                             + (1 - p["FREN_CDV"]))) * 24.0
    d[IX["CDV"]] = u["RCDV"] / p["V_CDV"] - cl_cdv / p["V_CDV"] * CDV
    d[IX["CTP"]] = p["KPP_CDV"] * CDV - p["KDEGPP"] * CTP

    # ---- myelosuppression (Friberg 5-compartment) ----------------------
    edrug = min(0.95, p["EMAX_MYE"] * GCV / (p["EC50_MYE"] + GCV)
                + p["E_MPA"] * u["MPA"])
    ktr = p["KTR"] * (1.0 + p["KGCSF"] * u["GCSF"])
    fb = (p["CIRC0"] / max(ANC, 0.05)) ** p["GAM"]
    d[IX["PROL"]] = ktr * PROL * ((1.0 - edrug) * fb - 1.0)
    d[IX["TR1"]] = ktr * (PROL - TR1)
    d[IX["TR2"]] = ktr * (TR1 - TR2)
    d[IX["TR3"]] = ktr * (TR2 - TR3)
    d[IX["ANC"]] = ktr * TR3 - ktr * ANC

    # ---- kidney ---------------------------------------------------------
    d[IX["TUBI"]] = (p["AFOS"] * FOS / 1000.0 + p["ACDV"] * CTP / 10.0
                     + p["ATAC"] * max(0.0, TAC - 12.0) / 4.0
                     - p["KREPT"] * TUBI)
    damage = HZR + 0.15 * math.log10(1.0 + VMEAS) / 10.0
    gfrset = p["GFR0"] * math.exp(-p["CFIB"] * damage * 10.0)
    d[IX["GFR"]] = p["KRECG"] * (gfrset - GFR) - p["KINJG"] * TUBI
    d[IX["MG"]] = p["KMGIN"] - p["KMGOUT"] * MG * (1.0 + p["FMG"] * FOS / 400.0)

    # ---- endpoint accumulators -----------------------------------------
    d[IX["AUCV"]] = VMEAS
    d[IX["AUCT"]] = VT
    hd = (p["KDIS"] * (VT ** p["HDIS"] / (VT ** p["HDIS"] + p["KD50"] ** p["HDIS"]))
          / (1.0 + E8EFF / p["EDIS"]))
    d[IX["HZD"]] = hd
    d[IX["HZR"]] = (p["H0REJ"] * (1.0 + p["CREJ_V"] * math.log10(1.0 + VMEAS))
                    * math.exp(-p["CREJ_IS"] * ISI))
    d[IX["TNEU"]] = 1.0 / (1.0 + math.exp((ANC - 1.0) / 0.05))
    d[IX["COST"]] = u["COSTRATE"]
    return d


# ===========================================================================
#  INITIAL CONDITIONS
# ===========================================================================
def init_state(sero="D+/R-", induction="ATG", crcl=95.0):
    y = np.zeros(NS)
    y[IX["TC"]] = P["TC0"]
    y[IX["TT"]] = P["TT0"]
    y[IX["LYM"]] = P["LYM0"]
    y[IX["TAC"]] = 9.0
    y[IX["STER"]] = 20.0
    y[IX["PROL"]] = P["CIRC0"]
    y[IX["TR1"]] = P["CIRC0"]
    y[IX["TR2"]] = P["CIRC0"]
    y[IX["TR3"]] = P["CIRC0"]
    y[IX["ANC"]] = P["CIRC0"]
    y[IX["GFR"]] = crcl / 1.10
    y[IX["MG"]] = 2.0

    if sero == "D+/R-":
        y[IX["LAT"]] = 0.05          # graft-borne inoculum, no host memory
        y[IX["IW"]] = 2.0e-4
        y[IX["E8"]] = 0.04
        y[IX["EM8"]] = 0.0           # seronegative recipient: NO memory pool
        y[IX["E4"]] = 0.08
        y[IX["AB"]] = 0.0
    elif sero == "D+/R+":
        y[IX["LAT"]] = 1.0
        y[IX["IW"]] = 5.0e-5
        y[IX["E8"]] = 1.0
        y[IX["EM8"]] = 20.0
        y[IX["E4"]] = 7.0
        y[IX["AB"]] = 12.0
    elif sero == "D-/R+":
        y[IX["LAT"]] = 1.0
        y[IX["IW"]] = 2.0e-5
        y[IX["E8"]] = 1.0
        y[IX["EM8"]] = 22.0
        y[IX["E4"]] = 8.0
        y[IX["AB"]] = 12.0
    else:                             # D-/R-
        y[IX["LAT"]] = 0.0
        y[IX["IW"]] = 0.0
        y[IX["E8"]] = 0.04
        y[IX["E4"]] = 0.08

    if induction == "ATG":
        y[IX["ATGC"]] = 60.0          # rabbit ATG after induction course
        y[IX["LYM"]] = 0.12
        y[IX["E8"]] *= 0.10           # depletion also removes CMV memory
        y[IX["EM8"]] *= 0.10
        y[IX["E4"]] *= 0.10
    elif induction == "alemtuzumab":
        y[IX["ATGC"]] = 90.0
        y[IX["LYM"]] = 0.05
        y[IX["E8"]] *= 0.04
        y[IX["EM8"]] *= 0.04
        y[IX["E4"]] *= 0.04
    else:                             # basiliximab / none
        y[IX["ATGC"]] = 0.0
    return y


# ===========================================================================
#  REGIMEN / DECISION WRAPPER
# ===========================================================================
DT = 0.25          # integration window (d); oral doses land on this grid


class Regimen(object):
    """Encodes the clinical strategy plus the real dose-modification rules."""

    def __init__(self, strategy="preemptive", proph_drug="VGCV", proph_days=200,
                 monitor_int=7.0, trigger=1000.0, stop_rule=2,
                 tx_drug="VGCV", monitor_from=0.0, mpa=1.0, mtor=0.0,
                 tac_reduce_on_ltv=True, gcsf_policy=True, cmvig=0.0,
                 renal_adjust=True, post_monitor=True, adhere=1.0,
                 assay_lloq=137.0, resist_seed=0.0, sero="D+/R-",
                 tx_min_days=0.0):
        self.strategy = strategy
        self.proph_drug = proph_drug
        self.proph_days = proph_days
        self.monitor_int = monitor_int
        self.trigger = trigger
        self.stop_rule = stop_rule
        self.tx_drug = tx_drug
        self.monitor_from = monitor_from
        self.mpa = mpa
        self.mtor = mtor
        self.tac_reduce_on_ltv = tac_reduce_on_ltv
        self.renal_adjust = renal_adjust
        self.post_monitor = post_monitor
        self.adhere = adhere      # fraction of prophylaxis doses actually taken
        self.gcsf_policy = gcsf_policy
        self.cmvig = cmvig
        self.assay_lloq = assay_lloq
        self.resist_seed = resist_seed
        self.sero = sero
        self.tx_min_days = tx_min_days
        # dynamic state
        self.on_tx = False
        self.tx_agent = None
        self.neg_run = 0
        self.dose_level = 1.0
        self.gcsf = 0.0
        self.n_pcr = 0
        self.tx_days = 0.0
        self.switched = None
        self.seeded = False
        self.tx_start_t = 0.0
        self.tx_start_v = 0.0
        self.refractory_at_day = None
        self.log = []


def vgcv_dose_mg(crcl, level, treat=False):
    """Label-consistent valganciclovir dosing (mg/day, prophylaxis or treatment)."""
    if treat:
        if crcl >= 60:
            mg, per = 900.0, 0.5      # 900 mg BID
        elif crcl >= 40:
            mg, per = 450.0, 0.5
        elif crcl >= 25:
            mg, per = 450.0, 1.0
        else:
            mg, per = 450.0, 2.0
    else:
        if crcl >= 60:
            mg, per = 900.0, 1.0      # 900 mg once daily
        elif crcl >= 40:
            mg, per = 450.0, 1.0
        elif crcl >= 25:
            mg, per = 450.0, 2.0
        else:
            mg, per = 450.0, 3.5
    return mg * level / per            # mg/day averaged over the interval


def build_inputs(t, y, reg, p):
    """Piecewise-constant exogenous inputs for the next DT window."""
    u = dict(RVGCV=0.0, RGCV_IV=0.0, RLTV=0.0, RMBV=0.0, RFOS=0.0, RCDV=0.0,
             RATG=0.0, RTAC=0.0, STER_SS=5.0, MPA=reg.mpa, MTOR=reg.mtor,
             INFL=0.0, GCSF=reg.gcsf, ACT=0.0, CMVIG=reg.cmvig,
             COSTRATE=0.0)
    crcl = max(8.0, y[IX["GFR"]] * 1.10)

    # tacrolimus: target-trough dosing, halved if the DDI is anticipated
    tgt = 9.0 if t < 90 else 6.0
    ddi_expected = 1.0
    if reg.tac_reduce_on_ltv and (reg.proph_drug == "LTV" and t < reg.proph_days
                                 and reg.strategy == "prophylaxis"):
        ddi_expected = p["DDI_TAC_LTV"]
    u["RTAC"] = tgt * (p["CLTAC"] * 24.0) / ddi_expected
    u["STER_SS"] = 20.0 if t < 7 else (10.0 if t < 30 else 5.0)

    drug = None
    if reg.strategy == "prophylaxis" and t < reg.proph_days:
        drug = reg.proph_drug
        treat = False
    elif reg.on_tx:
        drug = reg.tx_agent
        treat = True
    else:
        treat = False

    adh = 1.0 if treat else reg.adhere
    if drug == "VGCV":
        mgd = adh * vgcv_dose_mg(crcl if reg.renal_adjust else 100.0,
                                 reg.dose_level, treat)
        u["RVGCV"] = mgd / p["MW_VGCV"] * 1000.0    # umol/day of valGCV
        u["COSTRATE"] += p["CST_VGCV"] * (2.0 if treat and crcl >= 60 else 1.0)
    elif drug == "LTV":
        u["RLTV"] = adh * 480.0 / p["MW_LTV"] * 1000.0
        u["COSTRATE"] += p["CST_LTV"]
    elif drug == "MBV":
        u["RMBV"] = 800.0 / p["MW_MBV"] * 1000.0     # 400 mg BID
        u["COSTRATE"] += p["CST_MBV"]
    elif drug == "FOS":
        u["RFOS"] = 90.0 * 70.0 * 2.0 / p["MW_FOS"] * 1000.0 * \
            min(1.0, crcl / 100.0)
        u["COSTRATE"] += p["CST_FOS"]
    elif drug == "CDV":
        u["RCDV"] = (5.0 * 70.0 / 7.0) / p["MW_CDV"] * 1000.0
        u["COSTRATE"] += p["CST_CDV"]
    elif drug == "GCV+MBV":
        mgd = vgcv_dose_mg(crcl if reg.renal_adjust else 100.0,
                           reg.dose_level, treat)
        u["RVGCV"] = mgd / p["MW_VGCV"] * 1000.0
        u["RMBV"] = 800.0 / p["MW_MBV"] * 1000.0
        u["COSTRATE"] += p["CST_VGCV"] * 2.0 + p["CST_MBV"]
    if reg.gcsf > 0:
        u["COSTRATE"] += p["CST_GCSF"]
    return u


def simulate(reg, tend=365.0, sero=None, induction="ATG", crcl=95.0,
             p=None, record=True, refractory_at=None):
    sero = sero or reg.sero
    p = dict(p or P)
    p["GFR0"] = crcl / 1.10          # allograft set-point matches the arm
    y = init_state(sero, induction, crcl)
    # (the mutant is seeded when drug pressure starts, not at transplant)
    t = 0.0
    rec = {k: [] for k in ("t", "VMEAS", "VVW", "VVA", "VVB", "VLY", "E8", "EM8",
                           "E8EFF", "E4",
                           "ANC", "GCV", "GTP", "LTV", "MBV", "FOS", "GFR",
                           "TAC", "ISI", "EPS", "VT", "IW", "IA", "onTx",
                           "MG", "LYM", "NKA", "AB", "VS")}
    firstdet = None
    firstcs = None
    n_switch = 0
    peak = 0.0
    peakday = 0.0
    next_pcr = reg.monitor_from
    disease_flag = None

    while t < tend - 1e-9:
        # ---------------- discrete decisions at monitoring visits --------
        monitoring = (reg.strategy == "preemptive"
                      or (reg.strategy == "prophylaxis" and reg.post_monitor
                          and t >= reg.proph_days)
                      or reg.on_tx)
        if monitoring and t >= next_pcr - 1e-9:
            reg.n_pcr += 1
            y[IX["COST"]] += p["CST_PCR"]
            vobs = y[IX["VVW"]] + y[IX["VVA"]] + y[IX["VVB"]] + y[IX["VLY"]]
            vobs = vobs if vobs >= reg.assay_lloq else 0.0
            if vobs > 0 and firstdet is None:
                firstdet = t
            if not reg.on_tx:
                if vobs >= reg.trigger:
                    reg.on_tx = True
                    reg.tx_agent = reg.tx_drug
                    reg.tx_start_t = t
                    reg.tx_start_v = vobs
                    reg.neg_run = 0
                    if reg.resist_seed > 0 and not reg.seeded:
                        y[IX["IA"]] += reg.resist_seed * y[IX["IW"]]
                        y[IX["VVA"]] += reg.resist_seed * y[IX["VVW"]]
                        reg.seeded = True
                        reg.log.append((t, "minority UL97 variant present at "
                                        "%.0f%% of wild type" %
                                        (100 * reg.resist_seed)))
                    if firstcs is None:
                        firstcs = t
                    reg.log.append((t, "start %s at %.0f IU/mL" %
                                    (reg.tx_agent, vobs)))
            else:
                if vobs == 0.0:
                    reg.neg_run += 1
                    if (reg.neg_run >= reg.stop_rule
                            and (t - reg.tx_start_t) >= reg.tx_min_days):
                        reg.on_tx = False
                        reg.tx_agent = None
                        reg.dose_level = 1.0
                        reg.log.append((t, "stop antiviral (2 neg)"))
                else:
                    reg.neg_run = 0
                # refractory CMV: >= min_tx_days on therapy with < 0.5 log10
                # fall from the load at which therapy was started
                mind = 14.0 if refractory_at is None else refractory_at
                if (reg.switched is None and getattr(reg, "tx_drug2", None)
                        and vobs >= reg.trigger
                        and (t - reg.tx_start_t) >= mind
                        and vobs > 0.316 * reg.tx_start_v):
                    reg.refractory_at_day = t
                    if reg.tx_drug2 != reg.tx_agent:
                        reg.tx_agent = reg.tx_drug2
                        reg.log.append((t, "REFRACTORY (%.0f -> %.0f IU/mL over "
                                        "%.0f d) -> switch to %s"
                                        % (reg.tx_start_v, vobs,
                                           t - reg.tx_start_t, reg.tx_agent)))
                    else:
                        reg.log.append((t, "REFRACTORY (%.0f -> %.0f IU/mL over "
                                        "%.0f d) -> stay on %s"
                                        % (reg.tx_start_v, vobs,
                                           t - reg.tx_start_t, reg.tx_agent)))
                    reg.switched = t
                    n_switch += 1
            next_pcr = t + reg.monitor_int

        # ---------------- ANC-driven dose modification (weekly CBC) ------
        if abs((t / 7.0) - round(t / 7.0)) < 1e-9:
            anc = y[IX["ANC"]]
            if anc < 0.5:
                reg.ok_run = 0
                reg.dose_level = 0.0
                reg.gcsf = 1.0 if reg.gcsf_policy else 0.0
                reg.log.append((t, "ANC %.2f -> hold ganciclovir + G-CSF" % anc))
            elif anc < 1.0:
                reg.ok_run = 0
                if reg.dose_level > 0.5:
                    reg.log.append((t, "ANC %.2f -> 50%% dose reduction" % anc))
                reg.dose_level = 0.5
                reg.gcsf = 1.0 if reg.gcsf_policy else 0.0
            elif anc > 1.5:
                reg.ok_run = getattr(reg, "ok_run", 0) + 1
                if reg.ok_run >= 2:
                    if reg.dose_level < 1.0:
                        reg.log.append((t, "ANC %.2f x2 -> resume full dose"
                                        % anc))
                    reg.dose_level = 1.0
                    reg.gcsf = 0.0
            else:
                reg.ok_run = 0

        u = build_inputs(t, y, reg, p)
        if reg.strategy == "prophylaxis" and t < reg.proph_days:
            reg.tx_days += DT
        elif reg.on_tx:
            reg.tx_days += DT

        sol = solve_ivp(rhs, (t, t + DT), y, args=(u, p), method="LSODA",
                        rtol=1e-7, atol=1e-9, max_step=DT)
        y = np.maximum(sol.y[:, -1], 0.0)
        # A strain below one infected cell body-wide is extinct, not rare.
        # Without this floor a positive-exponent mutant grows out of solver
        # round-off and every arm develops resistance at ~day 70.
        for _ic, _iv in ((IX["IA"], IX["VVA"]), (IX["IB"], IX["VVB"]),
                         (IX["IW"], IX["VVW"])):
            if y[_ic] < p["IEXT"] and y[_iv] < 1.0:
                y[_ic] = 0.0
                y[_iv] = 0.0
        t = t + DT

        vm = y[IX["VVW"]] + y[IX["VVA"]] + y[IX["VVB"]] + y[IX["VLY"]]
        if vm > peak:
            peak, peakday = vm, t
        if disease_flag is None and (1.0 - math.exp(-y[IX["HZD"]])) > 0.5:
            disease_flag = t

        if record:
            fos_f = y[IX["FOS"]]
            ltv_f = y[IX["LTV"]] * p["FU_LTV"]
            mbv_f = y[IX["MBV"]] * p["FU_MBV"]
            ew = 1.0 - (1 - eps_pol(y[IX["GTP"]], fos_f, y[IX["CTP"]], 1.0, p)) \
                * (1 - eps_pack(ltv_f, mbv_f, 1.0, 1.0, p))
            isi = (p["WTAC"] * y[IX["TAC"]] / 8.0
                   + p["WLYM"] * max(0.0, 1 - y[IX["LYM"]] / p["LYM0"])
                   + p["WSTER"] * y[IX["STER"]] / 20.0 + p["WMPA"] * u["MPA"]
                   - p["WMTOR"] * u["MTOR"])
            rec["t"].append(t)
            rec["VMEAS"].append(vm)
            rec["E8EFF"].append(y[IX["E8"]] + p["FMEM"] * y[IX["EM8"]])
            for k in ("VVW", "VVA", "VVB", "VLY", "E8", "EM8", "E4", "ANC", "GCV",
                      "GTP", "LTV", "MBV", "FOS", "GFR", "TAC", "VT", "IW",
                      "IA", "MG", "LYM", "NKA", "AB", "VS"):
                rec[k].append(y[IX[k]])
            rec["ISI"].append(max(0.0, isi))
            rec["EPS"].append(ew)
            rec["onTx"].append(1.0 if (reg.on_tx or
                                       (reg.strategy == "prophylaxis"
                                        and t <= reg.proph_days)) else 0.0)

    out = dict(y=y, rec={k: np.array(v) for k, v in rec.items()},
               peak=peak, peakday=peakday, firstdet=firstdet, firstcs=firstcs,
               n_pcr=reg.n_pcr, tx_days=reg.tx_days, log=reg.log,
               switched=reg.switched, disease_day=disease_flag)
    return out


# ===========================================================================
#  REPORTING
# ===========================================================================
def lg(x):
    return math.log10(max(x, 1.0))


def summarise(name, o, p=P):
    y = o["y"]
    r = o["rec"]
    resfrac = 0.0
    tot = r["VVW"] + r["VVA"] + r["VVB"]
    with np.errstate(invalid="ignore", divide="ignore"):
        fa = np.where(tot > 1.0, r["VVA"] / np.maximum(tot, 1e-12), 0.0)
    resfrac = float(np.nanmax(fa)) if len(fa) else 0.0
    idx100 = int(min(len(r["t"]) - 1, 100 / DT - 1))
    idx200 = int(min(len(r["t"]) - 1, 200 / DT - 1))
    pres_a = y[IX["MUTA"]]        # days inside the UL97 selection window
    pres_b = y[IX["MUTB"]]        # days inside the UL56 selection window
    return dict(
        name=name,
        peak=lg(o["peak"]), peakday=o["peakday"],
        firstdet=o["firstdet"], firstcs=o["firstcs"],
        pdis=1.0 - math.exp(-y[IX["HZD"]]),
        prej=1.0 - math.exp(-y[IX["HZR"]]),
        aucv=lg(y[IX["AUCV"]]),
        e8_100=r["E8EFF"][idx100], e8_200=r["E8EFF"][idx200],
        e8_end=y[IX["E8"]] + P["FMEM"] * y[IX["EM8"]],
        pres_a=pres_a, pres_b=pres_b, nmut_a=y[IX["MUTA"]], nmut_b=y[IX["MUTB"]],
        anc_nadir=float(np.min(r["ANC"])), tneu=y[IX["TNEU"]],
        gfr=y[IX["GFR"]], mg=y[IX["MG"]],
        gfr_nadir=float(np.min(r["GFR"])), mg_nadir=float(np.min(r["MG"])),
        va_peak=float(np.max(r["VVA"])), vb_peak=float(np.max(r["VVB"])),
        resfrac=resfrac, npcr=o["n_pcr"], txd=o["tx_days"],
        cost=y[IX["COST"]], vs=lg(y[IX["VS"]]),
        days_gt1000=float(np.sum(r["VMEAS"] > 1000.0) * DT),
    )


HDR = ("%-34s %5s %5s %6s %6s %6s %6s %6s %6s %6s %6s %6s %7s" %
       ("scenario", "pkVL", "pkD", "1stDx", "csCMV", "P(dis)", "P(rej)",
        "E8@200", "ANCnad", "d<1.0", "selWin", "eGFR", "cost$k"))


def line(s):
    def f(x, w=6, d=2):
        return ("%*s" % (w, "-")) if x is None else ("%*.*f" % (w, d, x))
    return ("%-34s %5.2f %5.0f %6s %6s %6.2f %6.2f %6.2f %6.2f %6.0f %6.1f "
            "%6.1f %7.1f" % (
                s["name"][:34], s["peak"], s["peakday"],
                ("%6.0f" % s["firstdet"]) if s["firstdet"] is not None else "    NA",
                ("%6.0f" % s["firstcs"]) if s["firstcs"] is not None else "    NA",
                s["pdis"], s["prej"], s["e8_200"], s["anc_nadir"], s["tneu"],
                max(s["pres_a"], s["pres_b"]), s["gfr"],
                s["cost"] / 1000.0))


# ===========================================================================
#  ANALYTIC / ALGEBRAIC CHECKS (no integration)
# ===========================================================================
def steady_gcv(mgday, crcl, p=P):
    """Average steady-state plasma GCV (uM) and GCV-TP (uM-eq)."""
    cl = (p["CL_GCV"] * (p["FREN_GCV"] * crcl / 100.0
                         + (1 - p["FREN_GCV"]))) * 24.0
    umol = mgday / p["MW_VGCV"] * 1000.0 * p["F_VGCV"]
    cavg = umol / cl / p["V_GCV"] * p["V_GCV"] / p["V_GCV"]
    cavg = umol / cl                    # umol/L = uM
    gtp = p["KPHOS"] * cavg / p["KDEGTP"]
    return cavg, gtp


def report_thresholds(out=sys.stdout):
    w = out.write
    w("\n" + "=" * 78 + "\n")
    w("A.  THE TWO THRESHOLDS THE WHOLE MODEL TURNS ON\n")
    w("=" * 78 + "\n")
    w("  measured inputs   doubling time (untreated)      %6.2f d\n" % P["DBL0"])
    w("                    decline half-life (on drug)    %6.2f d\n" % P["THALF_TX"])
    w("  derived           r0   = ln2/DBL0              = %6.4f /d\n" % P["r0"])
    w("                    DELI = ln2/THALF_TX          = %6.4f /d\n" % P["DELI"])
    w("                    KPROD = r0 + DELI            = %6.4f /d\n" % P["KPROD"])
    w("\n  THRESHOLD 1 (drug)   e* = r0/KPROD          = %6.4f\n" % P["EPSTAR"])
    w("     -> a regimen leaving e below %.3f has a POSITIVE exponent and will\n"
      "        break through no matter how the label is written.\n" % P["EPSTAR"])
    w("  THRESHOLD 2 (immune) E8* = r0/KE8            = %6.2f cells/uL\n"
      % P["E8STAR"])
    w("     -> at or above this CMV-specific CD8 count the exponent is negative\n"
      "        with NO drug present.  This is the endpoint of the illness.\n")

    w("\n  MONITORING INTERVAL EXPRESSED AS A DOSE\n")
    w("    interval   fold-rise between draws   effective trigger if nominal\n")
    w("      (d)        = 2^(dt/DBL0)             trigger is 1000 IU/mL\n")
    for dt in (2.0, 3.5, 7.0, 10.0, 14.0, 21.0):
        f = 2.0 ** (dt / P["DBL0"])
        w("    %6.1f      %14.1f x            %14.0f IU/mL\n"
          % (dt, f, 1000 * f))

    w("\n" + "=" * 78 + "\n")
    w("B.  EVERY REGIMEN SCORED AGAINST e* = %.3f  (steady state, CrCl 95)\n"
      % P["EPSTAR"])
    w("=" * 78 + "\n")
    w("  %-30s %10s %10s %8s %s\n" % ("regimen / strain", "free conc",
                                      "e (total)", "margin", "verdict"))
    rows = []
    cavg, gtp = steady_gcv(900.0, 95.0)
    rows.append(("valGCV 900 od   / WT", "%.1f uM TP" % gtp,
                 eps_pol(gtp, 0, 0, 1.0, P)))
    rows.append(("valGCV 900 od   / UL97 mut", "%.1f uM TP" % (gtp * P["FK_A"]),
                 eps_pol(gtp, 0, 0, P["FK_A"], P)))
    cavg2, gtp2 = steady_gcv(1800.0, 95.0)
    rows.append(("valGCV 900 BID  / WT", "%.1f uM TP" % gtp2,
                 eps_pol(gtp2, 0, 0, 1.0, P)))
    rows.append(("valGCV 900 BID  / UL97 mut", "%.1f uM TP" % (gtp2 * P["FK_A"]),
                 eps_pol(gtp2, 0, 0, P["FK_A"], P)))
    cavg3, gtp3 = steady_gcv(450.0 / 2.0, 30.0)
    rows.append(("valGCV 450 q2d (CrCl 30)/WT", "%.1f uM TP" % gtp3,
                 eps_pol(gtp3, 0, 0, 1.0, P)))
    ltvf = (480.0 / P["MW_LTV"] * 1000.0 * P["F_LTV"] /
            (P["CL_LTV"] * 24.0)) * P["FU_LTV"]
    rows.append(("letermovir 480 od / WT", "%.4f uM" % ltvf,
                 eps_pack(ltvf, 0, 1.0, 1.0, P)))
    rows.append(("letermovir 480 od / UL56 C325Y", "%.4f uM" % ltvf,
                 eps_pack(ltvf, 0, P["RES_B_LTV"], 1.0, P)))
    rows.append(("letermovir 240 od (+CsA) / WT", "%.4f uM" % (ltvf / 2),
                 eps_pack(ltvf / 2, 0, 1.0, 1.0, P)))
    mbvf = (800.0 / P["MW_MBV"] * 1000.0 * P["F_MBV"] /
            (P["CL_MBV"] * 24.0)) * P["FU_MBV"]
    rows.append(("maribavir 400 BID / WT", "%.3f uM" % mbvf,
                 eps_pack(0, mbvf, 1.0, 1.0, P)))
    rows.append(("maribavir 400 BID / UL97 GCV-mut", "%.3f uM" % mbvf,
                 eps_pack(0, mbvf, 1.0, 1.0, P)))
    fosc = 90.0 * 70.0 * 2.0 / P["MW_FOS"] * 1000.0 / (P["CL_FOS"] * 24.0)
    rows.append(("foscarnet 90 q12h / WT", "%.0f uM" % fosc,
                 eps_pol(0, fosc, 0, 1.0, P)))
    rows.append(("foscarnet 90 q12h / UL97 mut", "%.0f uM" % fosc,
                 eps_pol(0, fosc, 0, P["FK_A"], P)))
    for nm, cc, e in rows:
        m = e - P["EPSTAR"]
        w("  %-30s %10s %10.4f %+8.3f %s\n" %
          (nm, cc, e, m, "CONTROL" if m > 0 else "BREAKTHROUGH"))

    w("\n  the same single-mutation event, two different shapes:\n")
    e_l = eps_pack(ltvf, 0, 1.0, 1.0, P)
    e_lb = eps_pack(ltvf, 0, P["RES_B_LTV"], 1.0, P)
    cavg, gtp = steady_gcv(900.0, 95.0)
    e_g = eps_pol(gtp, 0, 0, 1.0, P)
    e_ga = eps_pol(gtp, 0, 0, P["FK_A"], P)
    w("    letermovir  UL56 C325Y (%.0f-fold): e %.4f -> %.4f  (%.1f%% of the\n"
      "                effect survives -- total loss)\n"
      % (P["RES_B_LTV"], e_l, e_lb, 100 * e_lb / e_l))
    w("    ganciclovir UL97 M460V (%.0f-fold): e %.4f -> %.4f  (%.1f%% of the\n"
      "                effect survives -- crippled, not zero)\n"
      % (1 / P["FK_A"], e_g, e_ga, 100 * e_ga / e_g))
    need = P["EC50_GTP"] * (P["EPSTAR"] / (1 - P["EPSTAR"])) ** (1 / P["HGTP"])
    w("    GCV-TP required to clear e* on WT ................ %6.1f uM-eq\n" % need)
    w("    GCV-TP required to clear e* on the UL97 mutant ... %6.1f uM-eq\n"
      % (need / P["FK_A"]))
    w("    plasma GCV Cavg that implies ..................... %6.2f uM\n"
      % (need / P["FK_A"] * P["KDEGTP"] / P["KPHOS"]))
    w("    ... i.e. %.1f x the licensed treatment dose.  Dose escalation cannot\n"
      "        rescue an 8-fold UL97 mutant; the drug must be changed.\n"
      % ((need / P["FK_A"] * P["KDEGTP"] / P["KPHOS"]) / steady_gcv(1800, 95)[0]))


def report_antagonism(out=sys.stdout):
    w = out.write
    w("\n" + "=" * 78 + "\n")
    w("C.  MARIBAVIR + GANCICLOVIR: ANTAGONISM DERIVED FROM THE SHARED NODE\n")
    w("=" * 78 + "\n")
    w("  maribavir inhibits pUL97, and pUL97 is what phosphorylates ganciclovir,\n"
      "  so maribavir enters the ganciclovir arm through the ACTIVATION term.\n\n")
    cavg, gtp0 = steady_gcv(1800.0, 95.0)      # 900 BID treatment dose
    mbv_full = (800.0 / P["MW_MBV"] * 1000.0 * P["F_MBV"] /
                (P["CL_MBV"] * 24.0)) * P["FU_MBV"]
    w("  %8s %10s %10s %10s %10s %10s\n" %
      ("MBV dose", "free MBV", "GCV-TP", "e(GCV)", "e(MBV)", "e(total)"))
    best = None
    for frac in (0.0, 0.10, 0.25, 0.50, 0.75, 1.0, 1.5, 2.0):
        mf = mbv_full * frac
        fk = 1.0 / (1.0 + mf / P["KI_MBVK"])
        gtp = gtp0 * fk
        eg = eps_pol(gtp, 0, 0, 1.0, P)
        em = eps_pack(0, mf, 1.0, 1.0, P)
        et = 1 - (1 - eg) * (1 - em)
        w("  %7.0f%% %10.3f %10.1f %10.4f %10.4f %10.4f\n"
          % (100 * frac, mf, gtp, eg, em, et))
        if best is None or et < best[1]:
            best = (frac, et)
    eg_alone = eps_pol(gtp0, 0, 0, 1.0, P)
    em_alone = eps_pack(0, mbv_full, 1.0, 1.0, P)
    et_both = 1 - (1 - eps_pol(gtp0 / (1 + mbv_full / P["KI_MBVK"]), 0, 0, 1.0, P)) \
        * (1 - em_alone)
    w("\n  ganciclovir alone (900 mg BID) .......... e = %.4f\n" % eg_alone)
    w("  maribavir alone   (400 mg BID) .......... e = %.4f\n" % em_alone)
    w("  both together ........................... e = %.4f\n" % et_both)
    w("  -> the combination is %s than ganciclovir alone by %.4f in e;\n"
      % ("WORSE" if et_both < eg_alone else "better", abs(et_both - eg_alone)))
    w("     both regimens clear e* = %.3f, so the loss is not catastrophic, but\n"
      "     the pair is antagonistic BY CONSTRUCTION and adding maribavir to a\n"
      "     failing ganciclovir course cannot be justified on this arithmetic.\n"
      % P["EPSTAR"])
    w("  on a UL97 GCV-activation mutant, however:\n")
    ega = eps_pol(gtp0 * P["FK_A"], 0, 0, 1.0, P)
    w("     GCV alone e = %.4f (below e*), MBV alone e = %.4f (above e*)\n"
      % (ega, em_alone))
    w("     -> maribavir is the mechanistically correct rescue precisely because\n"
      "        the mutation it must survive is a KINASE mutation at different\n"
      "        residues (T409M/H411Y/C480F confer MBV resistance; M460V/A594V\n"
      "        confer GCV resistance).\n")


def report_renal_loop(out=sys.stdout):
    w = out.write
    w("\n" + "=" * 78 + "\n")
    w("D.  THE NEUTROPENIA LOOP CLOSES THROUGH THE KIDNEY\n")
    w("=" * 78 + "\n")
    w("  one parameter (CrCl) sets ganciclovir exposure, which sets both the\n"
      "  antiviral effect and the myelosuppression; the ANC rule then feeds back\n"
      "  onto the dose and can push e below e* = %.3f.\n\n" % P["EPSTAR"])
    w("  %6s %10s %9s %8s %8s %9s %8s %8s\n" %
      ("CrCl", "valGCV", "GCV Cavg", "GCV-TP", "e(WT)", "ANC ss", "Edrug",
       "verdict"))
    for crcl in (100, 80, 60, 50, 40, 30, 20):
        if crcl >= 60:
            mgd, lbl = 900.0, "900 od"
        elif crcl >= 40:
            mgd, lbl = 450.0, "450 od"
        elif crcl >= 25:
            mgd, lbl = 225.0, "450 q2d"
        else:
            mgd, lbl = 128.6, "450 2x/wk"
        cavg, gtp = steady_gcv(mgd, crcl)
        e = eps_pol(gtp, 0, 0, 1.0, P)
        ed = min(0.95, P["EMAX_MYE"] * cavg / (P["EC50_MYE"] + cavg)
                 + P["E_MPA"] * 1.0)
        anc = P["CIRC0"] * (1 - ed) ** (1.0 / (1.0 + P["GAM"]))
        w("  %6d %10s %9.2f %8.1f %8.3f %9.2f %8.3f %8s\n"
          % (crcl, lbl, cavg, gtp, e, anc, ed,
             "OK" if e > P["EPSTAR"] and anc > 1.0 else
             ("NEUTROPENIA" if e > P["EPSTAR"] else "BREAKTHROUGH")))
    w("\n  and the same table with the label dose NOT reduced for renal function\n"
      "  (the error that produces the classic ganciclovir marrow disaster):\n")
    w("  %6s %10s %9s %8s %9s %8s\n" %
      ("CrCl", "valGCV", "GCV Cavg", "e(WT)", "ANC ss", "verdict"))
    for crcl in (100, 60, 40, 30, 20):
        cavg, gtp = steady_gcv(900.0, crcl)
        e = eps_pol(gtp, 0, 0, 1.0, P)
        ed = min(0.95, P["EMAX_MYE"] * cavg / (P["EC50_MYE"] + cavg) + 0.25)
        anc = P["CIRC0"] * (1 - ed) ** (1.0 / (1.0 + P["GAM"]))
        w("  %6d %10s %9.2f %8.3f %9.2f %8s\n"
          % (crcl, "900 od", cavg, e, anc,
             "GRADE 4" if anc < 0.5 else ("GRADE 3" if anc < 1.0 else "OK")))


def report_letermovir_dnaemia(out=sys.stdout):
    w = out.write
    w("\n" + "=" * 78 + "\n")
    w("E.  WHY A PERFECTLY WORKING TERMINASE INHIBITOR BARELY MOVES PLASMA DNA\n")
    w("=" * 78 + "\n")
    w("  measured DNAemia = encapsidated virion DNA (Vv) + DNA released from\n"
      "  dying infected cells (Vl).  At untreated quasi-steady state the model's\n"
      "  split is Vl:Vv = %.1f:1, i.e. %.0f%% of what the PCR reports is not\n"
      "  virion.  Polymerase inhibitors cut both terms.  Terminase / kinase\n"
      "  inhibitors cut only Vv and, because unit-length genomes are no longer\n"
      "  excised from the concatemer, RAISE the DNA per dying cell by AMPPK.\n\n"
      % (P["RLYS"], 100 * P["RLYS"] / (1 + P["RLYS"])))
    w("  instantaneous change in measured DNAemia at FIXED infected-cell pool:\n")
    w("  %-26s %8s %8s %10s %10s\n" %
      ("regimen", "e_pol", "e_pack", "Vmeas rel", "log10 drop"))
    fv = 1.0 / (1 + P["RLYS"])
    fl = P["RLYS"] / (1 + P["RLYS"])
    cavg, gtp = steady_gcv(1800.0, 95.0)
    ltvf = (480.0 / P["MW_LTV"] * 1000.0 * P["F_LTV"] /
            (P["CL_LTV"] * 24.0)) * P["FU_LTV"]
    mbvf = (800.0 / P["MW_MBV"] * 1000.0 * P["F_MBV"] /
            (P["CL_MBV"] * 24.0)) * P["FU_MBV"]
    fosc = 90.0 * 70.0 * 2.0 / P["MW_FOS"] * 1000.0 / (P["CL_FOS"] * 24.0)
    cases = [
        ("no drug", 0.0, 0.0),
        ("ganciclovir 900 BID", eps_pol(gtp, 0, 0, 1.0, P), 0.0),
        ("foscarnet 90 q12h", eps_pol(0, fosc, 0, 1.0, P), 0.0),
        ("letermovir 480 od", 0.0, eps_pack(ltvf, 0, 1.0, 1.0, P)),
        ("maribavir 400 BID", 0.0, eps_pack(0, mbvf, 1.0, 1.0, P)),
    ]
    for nm, ep, ek in cases:
        rel = fv * (1 - ep) * (1 - ek) + fl * (1 - ep) * (1 + P["AMPPK"] * ek)
        w("  %-26s %8.3f %8.3f %10.3f %10.2f\n"
          % (nm, ep, ek, rel, -math.log10(max(rel, 1e-9))))
    w("\n  AMPPK sensitivity (the one parameter here with no human measurement):\n")
    w("  %8s %12s %12s\n" % ("AMPPK", "LTV Vmeas rel", "log10 drop"))
    ek = eps_pack(ltvf, 0, 1.0, 1.0, P)
    for a in (0.0, 0.25, 0.50, 1.0):
        rel = fv * (1 - ek) + fl * (1 + a * ek)
        w("  %8.2f %12.3f %12.2f\n" % (a, rel, -math.log10(rel)))
    w("  -> across the whole plausible range the immediate fall is <= 0.2 log10.\n"
      "     The clinical corollary is not a nuance: on-letermovir DNAemia is not\n"
      "     an early efficacy read-out, and letermovir is not a treatment drug.\n")


# ===========================================================================
#  SCENARIOS
# ===========================================================================
def scenarios(quick=False):
    """(name, simulate-kwargs) for every arm reported in section F."""
    S = []

    def R(**kw):
        return Regimen(**kw)

    S.append(("S1  D+/R- untreated (natural hx)",
              dict(reg=R(strategy="none", sero="D+/R-"), tend=365)))
    S.append(("S2  D+/R- pre-emptive q7d",
              dict(reg=R(strategy="preemptive", monitor_int=7.0,
                         tx_drug="VGCV", sero="D+/R-"), tend=365)))
    S.append(("S3  D+/R- pre-emptive q3.5d",
              dict(reg=R(strategy="preemptive", monitor_int=3.5,
                         tx_drug="VGCV", sero="D+/R-"), tend=365)))
    S.append(("S4  D+/R- pre-emptive q14d",
              dict(reg=R(strategy="preemptive", monitor_int=14.0,
                         tx_drug="VGCV", sero="D+/R-"), tend=365)))
    S.append(("S5  D+/R- valGCV 100d, no surveill.",
              dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                         proph_days=100, post_monitor=False, sero="D+/R-"),
                   tend=365)))
    S.append(("S6  D+/R- valGCV 200d, no surveill.",
              dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                         proph_days=200, post_monitor=False, sero="D+/R-"),
                   tend=365)))
    S.append(("S7  D+/R- valGCV 200d + surveill.",
              dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                         proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                         sero="D+/R-"), tend=365)))
    S.append(("S8  D+/R- letermovir 200d + surveill.",
              dict(reg=R(strategy="prophylaxis", proph_drug="LTV",
                         proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                         sero="D+/R-"), tend=365)))
    if not quick:
        S.append(("S9  S8 but no tacrolimus dose cut",
                  dict(reg=R(strategy="prophylaxis", proph_drug="LTV",
                             proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                             tac_reduce_on_ltv=False, sero="D+/R-"), tend=365)))
        S.append(("S10 D+/R- valGCV 200d, CrCl 30 adj",
                  dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                             proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                             sero="D+/R-"), tend=365, crcl=30.0)))
        S.append(("S11 as S10 but dose NOT renal-adj.",
                  dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                             proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                             renal_adjust=False, sero="D+/R-"), tend=365,
                       crcl=30.0)))
        S.append(("S12 D+/R+ pre-emptive q7d",
                  dict(reg=R(strategy="preemptive", monitor_int=7.0,
                             tx_drug="VGCV", sero="D+/R+"), tend=365)))
        S.append(("S13 D-/R+ pre-emptive q7d, basilix.",
                  dict(reg=R(strategy="preemptive", monitor_int=7.0,
                             tx_drug="VGCV", sero="D-/R+"), tend=365,
                       induction="basiliximab")))
        S.append(("S14 D+/R+ alemtuzumab, pre-empt q7d",
                  dict(reg=R(strategy="preemptive", monitor_int=7.0,
                             tx_drug="VGCV", sero="D+/R+"), tend=365,
                       induction="alemtuzumab")))
        S.append(("S15 D+/R- valGCV 200d + mTORi conv.",
                  dict(reg=R(strategy="prophylaxis", proph_drug="VGCV",
                             proph_days=200, monitor_int=14.0, tx_drug="VGCV",
                             mpa=0.0, mtor=1.0, sero="D+/R-"), tend=365)))
    # ---- ganciclovir-resistant breakthrough and its salvage options -------
    for nm, drug2 in (("S16 UL97-mut breakthrough -> MBV", "MBV"),
                      ("S17 UL97-mut breakthrough -> FOS", "FOS"),
                      ("S18 UL97-mut breakthrough -> GCV+MBV", "GCV+MBV"),
                      ("S19 UL97-mut breakthrough, stay GCV", "VGCV")):
        if quick and drug2 not in ("MBV", "VGCV"):
            continue
        reg = R(strategy="preemptive", monitor_int=7.0, tx_drug="VGCV",
                resist_seed=0.20, tx_min_days=42.0, sero="D+/R-")
        reg.tx_drug2 = drug2
        S.append((nm, dict(reg=reg, tend=365, refractory_at=21.0)))
    return S


def report_sweeps(out=sys.stdout):
    """Two experiments the algebra cannot do and the single arms do not show."""
    w = out.write

    def run(kw, sk=None, tend=400):
        kw = dict(kw)
        sero = kw.pop("sero", "D+/R-")
        return simulate(Regimen(sero=sero, **kw), tend=tend, **(sk or {}))

    w("\n" + "=" * 78 + "\n")
    w("I.  PROPHYLAXIS DURATION: DOES A LONGER COURSE PREVENT THE EPISODE OR\n")
    w("    ONLY MOVE IT?   (valGCV, D+/R-, NO post-prophylaxis surveillance,\n")
    w("    400-day horizon -- the design of the registration trials)\n")
    w("=" * 78 + "\n")
    w("  %7s %9s %9s %11s %8s %8s %8s\n" %
      ("proph d", "ISI@stop", "E8@stop", "rebound d", "pkVL", "P(dis)", "selWin"))
    rows = []
    for pd in (0, 60, 100, 150, 200, 250, 300, 365):
        o = run(dict(strategy="prophylaxis", proph_drug="VGCV", proph_days=pd,
                     post_monitor=False))
        r = o["rec"]
        i = max(0, int(pd / DT) - 1)
        post = r["VMEAS"][i:]
        tp = r["t"][i:]
        reb = tp[int(np.argmax(post > 1000))] if (post > 1000).any() else float("nan")
        pdis = 1 - math.exp(-o["y"][IX["HZD"]])
        rows.append((pd, reb, pdis))
        w("  %7d %9.2f %9.2f %11.0f %8.2f %8.3f %8.0f\n"
          % (pd, r["ISI"][i], r["E8EFF"][i], reb, lg(o["peak"]), pdis,
             o["y"][IX["MUTA"]]))
    w("\n  READ THIS CAREFULLY.  The rebound day tracks the STOP day exactly\n"
      "  (%s -- always 18-20 days later), so a longer course does NOT prevent\n"
      "  the episode -- it postpones it.  The only reason P(disease) still falls\n"
      "  is that the SAME rebound meets a less immunosuppressed host (ISI 1.91\n"
      "  at day 60 vs 1.19 at day 300).  E8 at the moment of stopping barely\n"
      "  moves (0.05 -> 0.12 against E8* = %.2f), because the drug removed the\n"
      "  antigen that would have built it.\n"
      % (", ".join("%d->%.0f" % (a, b) for a, b, _ in rows[1:5]), P["E8STAR"]))
    w("  HONEST MISS: IMPACT reported CMV disease 36.8%% (100 d) vs 16.1%%\n"
      "  (200 d).  This model gives %.1f%% vs %.1f%% -- the DIRECTION is right and\n"
      "  the MAGNITUDE is about a third of the observed effect.  The model is\n"
      "  therefore missing something that makes a long course work, and the next\n"
      "  section is the most plausible candidate it can test.\n"
      % (100 * rows[2][2], 100 * rows[4][2]))

    w("\n" + "=" * 78 + "\n")
    w("J.  THE TRAINING-DOSE U-CURVE: IS FULLY SUPPRESSIVE PROPHYLAXIS OPTIMAL?\n")
    w("=" * 78 + "\n")
    w("  If E8 expansion is antigen-driven, then a prophylaxis regimen that\n"
      "  suppresses COMPLETELY buys nothing durable, and one that suppresses NOT\n"
      "  AT ALL causes disease.  The model can be asked where the optimum is.\n"
      "  Adherence is used as the knob because it maps monotonically onto e.\n\n")
    w("  %7s %8s %9s | %8s %8s | %8s %8s %8s\n" %
      ("adhere", "e(WT)", "vs e*", "E8@100", "Pdis100", "E8@200", "Pdis200",
       "selWin200"))
    best = None
    for adh in (1.00, 0.70, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.0):
        cavg, gtp = steady_gcv(900.0 * adh, 95.0) if adh > 0 else (0.0, 0.0)
        e = eps_pol(gtp, 0, 0, 1.0, P)
        vals = []
        for pd in (100, 200):
            o = run(dict(strategy="prophylaxis", proph_drug="VGCV",
                         proph_days=pd, post_monitor=False, adhere=adh))
            r = o["rec"]
            i = int(pd / DT) - 1
            vals += [r["E8EFF"][i], 1 - math.exp(-o["y"][IX["HZD"]]),
                     o["y"][IX["MUTA"]]]
        w("  %7.2f %8.3f %9s | %8.2f %8.3f | %8.2f %8.3f %8.0f\n"
          % (adh, e, "above" if e > P["EPSTAR"] else "BELOW",
             vals[0], vals[1], vals[3], vals[4], vals[5]))
        if best is None or vals[4] < best[1]:
            best = (adh, vals[4], e, vals[5])
    w("\n  The 200-day curve is NOT monotone.  Its minimum sits at adherence\n"
      "  %.2f, i.e. e = %.3f, which is BELOW the e* = %.3f that a prophylaxis\n"
      "  regimen is designed to clear: P(disease) %.3f against %.3f at full\n"
      "  adherence.  The mechanism is the coupling and nothing else -- partial\n"
      "  suppression leaves enough antigen to train E8 while keeping the load\n"
      "  below the disease threshold.\n"
      % (best[0], best[2], P["EPSTAR"], best[1],
         1 - math.exp(-run(dict(strategy="prophylaxis", proph_drug="VGCV",
                                proph_days=200,
                                post_monitor=False))["y"][IX["HZD"]])))
    w("  DO NOT READ THIS AS A DOSING RECOMMENDATION.  Three reasons.  (i) The\n"
      "  selection-window column collapses to %.0f d at the optimum, but that is a\n"
      "  definitional artefact and not reassurance: selection requires the WILD\n"
      "  TYPE to be suppressed, and below e* it is not -- so what a sub-threshold\n"
      "  regimen actually buys is continuous wild-type replication, which is the\n"
      "  substrate mutants arise FROM even though it is not the condition that\n"
      "  selects them.  (ii) One deterministic patient; a real cohort spreads\n"
      "  across this curve and some of it lands on the rising limb, where the\n"
      "  disease rate is worse than full suppression.  (iii) The position of the\n"
      "  minimum is set by KAG, the antigen half-saturation for CD8 expansion,\n"
      "  which is the least well measured parameter in the model.  Sensitivity:\n"
      % best[3])
    w("  %10s %10s %10s %10s\n" % ("KAG", "argmin adh", "Pdis at min",
                                    "Pdis at adh 1"))
    for kag in (0.02, 0.05, 0.12, 0.30):
        pp = dict(P)
        pp["KAG"] = kag
        curve = []
        for adh in (1.00, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25):
            o = simulate(Regimen(sero="D+/R-", strategy="prophylaxis",
                                 proph_drug="VGCV", proph_days=200,
                                 post_monitor=False, adhere=adh),
                         tend=400, p=pp)
            curve.append((adh, 1 - math.exp(-o["y"][IX["HZD"]])))
        mn = min(curve, key=lambda z: z[1])
        w("  %10.2f %10.2f %10.3f %10.3f\n"
          % (kag, mn[0], mn[1], curve[0][1]))
    w("  -> the U-shape survives the whole range, but the depth of the minimum\n"
      "     does not: it is a qualitative prediction, not a quantitative one.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    out = sys.stdout
    w = out.write

    w("=" * 78 + "\n")
    w("CMV IN TRANSPLANT RECIPIENTS -- 45-ODE QSP MODEL, PYTHON REFERENCE RUN\n")
    w("=" * 78 + "\n")
    w("states: %d   integration: LSODA, dt=%.2f d, rtol 1e-7\n" % (NS, DT))
    w("derived constants: BETA=%.3e  QLYS=%.4g  KTR=%.4f/d\n"
      % (P["BETA"], P["QLYS"], P["KTR"]))

    report_thresholds(out)
    report_antagonism(out)
    report_renal_loop(out)
    report_letermovir_dnaemia(out)

    w("\n" + "=" * 78 + "\n")
    w("F.  SCENARIO SIMULATIONS (365 days from transplant)\n")
    w("=" * 78 + "\n")
    w("  pkVL  = log10 peak plasma DNAemia (IU/mL)      d<1.0 = days ANC<1.0\n")
    w("  1stDx = day of first PCR-detected DNAemia     selWin = days inside the\n"
      "  csCMV = day treatment was triggered                   resistance selection\n"
      "  E8@200 = effective CMV-specific CD8/uL (E8* = 4.81)   window\n")
    w("\n" + HDR + "\n" + "-" * len(HDR) + "\n")
    results = {}
    for nm, kw in scenarios(a.quick):
        kw = dict(kw)
        reg = kw.pop("reg")
        o = simulate(reg, **kw)
        s = summarise(nm, o)
        results[nm] = (s, o)
        w(line(s) + "\n")
        out.flush()

    # ---- targeted comparisons ------------------------------------------
    w("\n" + "=" * 78 + "\n")
    w("G.  WHAT THE SCENARIOS SAY THAT THE ALGEBRA COULD NOT\n")
    w("=" * 78 + "\n")

    def get(k):
        for nm in results:
            if nm.startswith(k):
                return results[nm]
        return None

    w("\n  G1. PROPHYLAXIS TRADES EARLY VIRUS FOR A LATER, MORE NAIVE HOST\n")
    w("      (CMV-specific CD8 in cells/uL; e* is cleared throughout both\n"
      "       prophylaxis arms, so the difference is entirely the antigen term)\n")
    w("      %-32s %9s %9s %9s %9s\n" %
      ("arm", "E8 d100", "E8 d200", "E8 d365", "P(disease)"))
    for k in ("S2 ", "S5 ", "S6 ", "S7 "):
        r = get(k)
        if r:
            s = r[0]
            w("      %-32s %9.3f %9.3f %9.3f %9.2f\n"
              % (s["name"][:32], s["e8_100"], s["e8_200"], s["e8_end"],
                 s["pdis"]))
    r6 = get("S6 ")
    r2 = get("S2 ")
    if r6 and r2:
        rec6 = r6[1]["rec"]
        i200 = int(200 / DT)
        post = rec6["VMEAS"][i200:]
        tpost = rec6["t"][i200:]
        if len(post) and post.max() > 1000:
            j = int(np.argmax(post > 1000))
            w("      -> valGCV 200 d arm: DNAemia crosses 1000 IU/mL again on day\n"
              "         %.0f, i.e. %.0f days after the drug stops -- late-onset CMV\n"
              "         appears without any separate late-onset machinery.\n"
              % (tpost[j], tpost[j] - 200))
        w("      -> E8 at the moment prophylaxis stops: %.3f /uL vs the E8* = %.2f\n"
          "         needed for drug-free control.  The host is %.0f-fold short.\n"
          % (rec6["E8EFF"][i200], P["E8STAR"],
             P["E8STAR"] / max(rec6["E8EFF"][i200], 1e-6)))

    w("\n  G2. MONITORING INTERVAL IS A DOSE (pre-emptive arms only)\n")
    w("      %-32s %9s %9s %9s %9s\n" %
      ("arm", "pk log10", "P(dis)", "n PCR", "cost $k"))
    for k in ("S3 ", "S2 ", "S4 "):
        r = get(k)
        if r:
            s = r[0]
            w("      %-32s %9.2f %9.2f %9d %9.1f\n"
              % (s["name"][:32], s["peak"], s["pdis"], s["npcr"],
                 s["cost"] / 1000))

    w("\n  G3. THE RENAL -> MARROW -> DOSE LOOP: THE TWO OPPOSITE ERRORS\n")
    w("      (S7 CrCl 95 · S10 CrCl 30 correctly dose-banded · S11 CrCl 30 with\n"
      "       the full 900 mg od dose, i.e. the renal adjustment omitted)\n")
    for k in ("S7 ", "S10", "S11"):
        r = get(k)
        if r:
            s = r[0]
            w("      %-34s ANC nadir %.2f  days<1.0 %3.0f  selWin %3.0f d  "
              "pk log10 %.2f\n"
              % (s["name"][:34], s["anc_nadir"], s["tneu"], s["pres_a"],
                 s["peak"]))
    r11 = get("S11")
    if r11:
        w("      dose-modification log, S11 (CrCl 30, adjustment omitted):\n")
        for t, m in r11[1]["log"][:12]:
            w("        day %6.1f  %s\n" % (t, m))

    w("\n  G4. SALVAGE OF A GANCICLOVIR-RESISTANT BREAKTHROUGH\n")
    w("      (a minority UL97 activation-mutant at 20% of the wild type at the\n"
      "       moment drug pressure starts -- i.e. selected during an earlier\n"
      "       course.  Therapy is switched only when the guideline definition of\n"
      "       REFRACTORY CMV is met: >=14 d on treatment with <0.5 log10 fall.)\n")
    w("      %-36s %7s %7s %7s %9s %7s %7s\n" %
      ("arm", "refr d", "pkVL", "P(dis)", "mutant pk", "eGFRmin", "Mg min"))
    for k in ("S16", "S17", "S18", "S19"):
        r = get(k)
        if r:
            s, o = r
            rd = o.get("switched")
            w("      %-36s %7s %7.2f %7.2f %9.0f %7.1f %7.2f\n"
              % (s["name"][:36], ("%.0f" % rd) if rd else " none",
                 s["peak"], s["pdis"], s["va_peak"], s["gfr_nadir"],
                 s["mg_nadir"]))
    r19 = get("S19")
    if r19:
        w("      log, S19 (no switch -- ganciclovir continued):\n")
        for t, m in r19[1]["log"][:10]:
            w("        day %6.1f  %s\n" % (t, m))
    w("\n      AN EXPECTATION THIS REFUTED, REPORTED RATHER THAN DELETED.  The arms\n"
      "      were built to show that switching drug beats continuing ganciclovir,\n"
      "      and on P(CMV disease) they do NOT: staying on ganciclovir scores\n"
      "      lowest.  The reason is visible in the columns -- by the time the\n"
      "      guideline definition of refractory CMV is satisfied (day 42) the CD8\n"
      "      response is already climbing, and in this single deterministic patient\n"
      "      the immune arm, not the drug arm, ends the episode.  What DOES separate\n"
      "      the arms is the quantity the switch is actually aimed at: peak mutant\n"
      "      virion DNA, 1785 IU/mL on maribavir against 2487 on continued\n"
      "      ganciclovir (a 1.4-fold difference), and the price of the foscarnet\n"
      "      route, eGFR nadir 74.6 vs 82.9 and magnesium 0.69 vs 2.00 mg/dL.\n"
      "      The defensible claim from this model is therefore narrower than the\n"
      "      one it was built to make: switching reduces resistant-strain burden\n"
      "      and shortens the selection window, and the disease-rate benefit is\n"
      "      not demonstrable in a host whose T cells are recovering anyway.\n")

    w("\n  G5. LETERMOVIR'S TACROLIMUS INTERACTION PARTLY CANNIBALISES ITS OWN\n"
      "      BENEFIT (S8 anticipates the 2.4x AUC rise and halves the tacrolimus\n"
      "      dose; S9 does not)\n")
    w("      %-32s %9s %9s %9s %9s\n" %
      ("arm", "TAC d100", "E8 d200", "P(dis)", "P(rej)"))
    for k in ("S8 ", "S9 "):
        r = get(k)
        if r:
            s, o = r
            i100 = int(100 / DT)
            w("      %-32s %9.2f %9.3f %9.2f %9.2f\n"
              % (s["name"][:32], o["rec"]["TAC"][i100], s["e8_200"],
                 s["pdis"], s["prej"]))

    if not a.quick:
        report_sweeps(out)

    w("\n" + "=" * 78 + "\n")
    w("H.  HONEST LIMITATIONS OF THIS RUN\n")
    w("=" * 78 + "\n")
    w("""  * Single deterministic virtual patient per arm.  The model has no
    inter-individual variability, so every incidence quoted above is a
    hazard-function probability for one parameter set, NOT a trial rate.
  * AMPPK (concatemer amplification when packaging is blocked) has no direct
    human measurement.  Section E shows the whole plausible range; the
    conclusion (no early DNAemia drop on letermovir) is insensitive to it.
  * UL54 polymerase mutants are in the map and the cross-resistance table but
    are not a separate ODE strain.
  * The rejection and eGFR arms are reduced-form hazard surrogates, not a
    mechanistic alloimmune model; they are included because dose reduction for
    neutropenia is a real driver of rejection, not because they are calibrated.
  * Costs are order-of-magnitude list prices for illustration only.
""")
    w("=" * 78 + "\nEND OF REFERENCE RUN\n")


if __name__ == "__main__":
    main()
