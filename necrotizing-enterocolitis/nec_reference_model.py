#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
===============================================================================
nec_reference_model.py
Necrotising Enterocolitis (NEC) of the preterm infant
신생아 괴사성 장염 — INDEPENDENT PYTHON REFERENCE MODEL
===============================================================================

WHY THIS FILE EXISTS
--------------------
This build environment has no R runtime, so `nec_mrgsolve_model.R` cannot be
executed here.  Committing an un-integrated ODE model would be dishonest, so
every equation in the mrgsolve file is re-implemented in this file, term for
term, and actually integrated (fixed-step RK4, pure standard library, no
numpy).  Every number quoted in `README.md` is produced by running this file;
the captured log is `nec_reference_output.txt` and the machine-readable results
are `nec_scenario_results.json` / `nec_population_results.json`.

-------------------------------------------------------------------------------
THE ONE STRUCTURAL CLAIM
-------------------------------------------------------------------------------
NEC is not written here as "inflammation that gets worse".  It is written as a
CLOSED POSITIVE LOOP on a single mucosal state variable E (enterocyte / villus
integrity), with the loop gain set by the luminal pathobiont load B:

        E  ->  BI = E·(TJ)·(MUC)          barrier integrity
            ->  Pb = Pmin + (Pmax-Pmin)·(1-BI)^h        permeability
            ->  Jtr = Pb · B                            translocation flux
            ->  TLR4s = TLR4expr · Jtr/(Jtr + Ktlr)     innate signal
            ->  (a) apoptosis  ↑   kA·E·INJ
                (b) proliferation ↓  1/(1+(INJ/Ki)^2)
            ->  E ↓                                     ... and round again

Both of the injury arms (a) and (b) act on E, and E is the thing that sets Pb.
That is the loop.  Write the 1-D reduced field (TJ = MUC = 1, B held):

    g(E) = kE·(1-E)·Ftroph / (1 + (INJ/Ki)^2)  -  kA·E·INJ,   INJ = TLR4s(E;B)

g(0) > 0 and g(1) < 0 always, so the number of roots in (0,1) is odd.  For low
B there is exactly ONE root (healthy).  For high enough B a saddle-node
bifurcation creates a PAIR of extra roots — a separatrix E* and a necrotic
attractor.  For higher B still the healthy root annihilates with E*.

Five results then follow as ARITHMETIC rather than as assertions:

 (1) TWO CRITICAL LOADS EXIST, NOT ONE.  B_lo = the load at which NEC first
     becomes *possible* (bistability appears); B_hi = the load at which NEC
     becomes *inevitable* (the healthy root is destroyed).  Between them the
     infant's fate depends on the initial condition, not on the load — which is
     why NEC is all-or-none, why it clusters, and why two babies on the same
     unit with the same feeds do different things.  The model computes B_lo and
     B_hi for every (gestational age × milk type × probiotic) cell.

 (2) THRESHOLD-MOVERS AND STATE-MOVERS ARE DIFFERENT DRUGS AND MULTIPLY.
     Human milk enters ONLY through Ktlr (HMO, sIgA) and Ftroph (EGF/HB-EGF) —
     it moves the bifurcation points.  Probiotics enter ONLY through B — they
     move the state.  Nothing in the model lets one substitute for the other,
     and because one scales the axis while the other slides along it, their
     joint effect is a PRODUCT of risk ratios, not a sum.  The model computes
     both and shows the additive prediction is wrong.

 (3) FEEDING IS TWO-SIGNED, SO THE OPTIMUM ADVANCE RATE IS INTERIOR.  Enteral
     substrate appears with a + sign in Ftroph (trophic) and with a + sign in
     dB/dt (it is bacterial substrate) and with a - sign in PERF (splanchnic
     oxygen demand).  Being NPO therefore is not safe: Ftroph falls to its
     atrophy floor and E drifts toward E*.  The model finds an interior optimum
     and shows it is FLAT for a 30-weeker and STEEP for a 25-weeker — i.e. the
     null result of the large advance-rate trials and the clinical fear of fast
     feeds are the same equation evaluated at two gestational ages.

 (4) EMPIRICAL ANTIBIOTICS CHANGE SIGN AT A COMPUTABLE DURATION.  They lower B
     (good, immediate, first-order in concentration) and they lower C (bad,
     delayed: C makes SCFA, SCFA suppresses B through pH and frees carrying
     capacity when it goes).  Two opposite-signed integrals with different time
     constants must cross.  The model computes where.

 (5) TIME FROM CROSSING THE SEPARATRIX TO BELL III IS LOGARITHMIC IN THE
     DISTANCE CROSSED, HENCE NEARLY CONSTANT.  Near E* the field is linear,
     g(E) ≈ -λ·(E-E*), so the transit time is (1/λ)·ln(distance) — a 30-fold
     range of initial displacement compresses into a ~2-fold range of hours.
     That is why Bell II -> III is "always about a day", and why a screening
     biomarker that fires at Bell II buys so little time.

-------------------------------------------------------------------------------
STATE VECTOR (36 ODEs)
-------------------------------------------------------------------------------
MUCOSAL BARRIER      E TJ MUC IgA DEF
LUMINAL ECOLOGY      B C P SUB SCFA GAS HMO
INNATE SIGNALLING    TLR4s IL1B TNF IL8 IL10 PAF NEUT NOx
PERFUSION / LESION   PERF NECa PNEU
SYSTEMIC             LPSp PLT CRP LAC NIN BILI WT KIN
DRUG PK              AMPc GENc GENp MTZc INDc IBUc DEXc PGE

All concentrations that have no natural clinical unit are dimensionless
indices normalised so that 1.0 is the healthy term-infant reference or the
saturating pathological value; the ones that DO have units carry them in the
comment on their line.  Nothing in the analysis depends on the arbitrary ones.

-------------------------------------------------------------------------------
HONEST LIMITATIONS
-------------------------------------------------------------------------------
* Parameters are hand-calibrated to reproduce published *aggregate* endpoints
  (NEC incidence by GA band, onset-day distribution, human-milk and probiotic
  risk ratios, surgical fraction, case fatality).  They are NOT fitted to
  individual patient data and the covariance structure is invented.
* B and C are single lumped taxa ("Enterobacteriaceae-like" and "obligate
  anaerobe / Bifidobacterium-like").  Real NEC dysbiosis is strain-level.
* The bifurcation analysis holds B, HMO, IgA fixed while E relaxes.  In the
  full 36-ODE system B is itself dynamic, so the computed B_lo / B_hi are the
  *quasi-static* thresholds; the dynamic simulation crosses them at slightly
  different loads because B moves on a comparable timescale.  Both are
  reported so the gap is visible rather than hidden.
* This is a teaching / hypothesis-generating model.  It must not be used for
  any clinical decision.
===============================================================================
"""

import json
import math
import os
import sys

# =============================================================================
# 0.  STATE INDEXING
# =============================================================================
NAMES = [
    # --- mucosal barrier ------------------------------------------------------
    "E",     # 0  enterocyte / villus integrity            [0-1 index]
    "TJ",    # 1  tight-junction (ZO-1 / claudin-3) index  [0-1]
    "MUC",   # 2  mucus layer adequacy                     [0-1]
    "IgA",   # 3  luminal secretory IgA                    [0-1 index]
    "DEF",   # 4  Paneth-cell defensin (HD5/HD6) index     [0-1]
    # --- luminal ecology ------------------------------------------------------
    "B",     # 5  pathobiont (Enterobacteriaceae) load     [1e9 CFU/g stool]
    "C",     # 6  commensal obligate anaerobe load         [1e9 CFU/g]
    "P",     # 7  administered probiotic load              [1e9 CFU/g]
    "SUB",   # 8  unabsorbed fermentable substrate pool    [g/kg]
    "SCFA",  # 9  luminal butyrate-equivalent SCFA         [mM]
    "GAS",   # 10 luminal / intramural gas index           [0-1+]
    "HMO",   # 11 human-milk oligosaccharide pool          [g/kg]
    # --- innate signalling ----------------------------------------------------
    "TLR4s", # 12 TLR4 -> MyD88 -> NF-kB activity          [0-1 index]
    "IL1B",  # 13 mucosal IL-1beta                         [index]
    "TNF",   # 14 mucosal TNF-alpha                        [index]
    "IL8",   # 15 mucosal IL-8 / CXCL8                     [index]
    "IL10",  # 16 mucosal IL-10 (counter-regulatory)       [index]
    "PAF",   # 17 platelet-activating factor               [index]
    "NEUT",  # 18 mucosal neutrophil infiltrate            [index]
    "NOx",   # 19 iNOS-derived NO / peroxynitrite          [index]
    # --- perfusion and lesion -------------------------------------------------
    "PERF",  # 20 mesenteric perfusion (1 = normal)        [0-1]
    "NECa",  # 21 transmural necrotic bowel area fraction  [0-1]
    "PNEU",  # 22 pneumatosis intestinalis radiographic ix [0-1]
    # --- systemic -------------------------------------------------------------
    "LPSp", # 23 plasma endotoxin / bacteraemia index      [index]
    "PLT",   # 24 platelet count                           [1e9/L]
    "CRP",   # 25 C-reactive protein                       [mg/L]
    "LAC",   # 26 arterial lactate                         [mmol/L]
    "NIN",   # 27 neuro-inflammation (NDI driver) index    [index]
    "BILI",  # 28 conjugated bilirubin (TPN cholestasis)   [mg/dL]
    "WT",    # 29 body weight                              [kg]
    "KIN",   # 30 aminoglycoside tubular injury index      [index]
    # --- drug PK --------------------------------------------------------------
    "AMPc",  # 31 ampicillin plasma                        [mg/L]
    "GENc",  # 32 gentamicin plasma (central)              [mg/L]
    "GENp",  # 33 gentamicin peripheral / cortex           [mg/L-eq]
    "MTZc",  # 34 metronidazole plasma                     [mg/L]
    "INDc",  # 35 indomethacin plasma                      [mg/L]
    "IBUc",  # 36 ibuprofen plasma                         [mg/L]
    "DEXc",  # 37 dexamethasone / hydrocortisone plasma    [ug/L]
    "PGE",   # 38 mucosal vasodilator prostanoid index     [0-1]
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)

(iE, iTJ, iMUC, iIgA, iDEF, iB, iC, iP, iSUB, iSCFA, iGAS, iHMO, iTLR, iIL1,
 iTNF, iIL8, iIL10, iPAF, iNEU, iNOx, iPERF, iNEC, iPNE, iLPS, iPLT, iCRP,
 iLAC, iNIN, iBIL, iWT, iKIN, iAMP, iGEN, iGEP, iMTZ, iIND, iIBU, iDEX,
 iPGE) = range(NST)


# =============================================================================
# 1.  PARAMETERS
# =============================================================================
def base_params():
    p = {}

    # ---- gestational / postnatal maturation ---------------------------------
    p["GA"] = 28.0          # gestational age at birth [wk]
    p["BW"] = 1.05          # birth weight [kg]
    p["PMA50"] = 30.5       # PMA at half-maximal mucosal maturation [wk]
    p["nMAT"] = 6.0         # Hill coefficient of maturation
    p["phiTLR"] = 0.55      # fraction of TLR4 over-expression removed by MAT
    p["TLR4max"] = 1.00     # immature-gut TLR4 signalling capacity

    # ---- barrier: enterocyte / villus integrity E ---------------------------
    p["kE"] = 0.45          # crypt renewal rate constant [1/d]
    p["kA"] = 0.90          # injury-driven enterocyte loss [1/d]
    p["Ki"] = 0.25          # injury at which proliferation is half-blocked
    p["KEap"] = 0.08        # E at which injury-driven loss is half-maximal
    p["nKi"] = 2.0          # steepness of proliferation block
    p["troph_floor"] = 0.25 # Ftroph when completely NPO (mucosal atrophy)
    p["Kfeed"] = 40.0       # feed volume at half-maximal trophic drive [mL/kg/d]
    p["eEGF"] = 0.10        # extra trophic drive from milk EGF / HB-EGF
    p["eSCFA"] = 0.35       # extra trophic drive from butyrate
    p["KsSCFA"] = 8.0       # SCFA at half-maximal trophic effect [mM]

    # ---- barrier accessories ------------------------------------------------
    p["kTJ"] = 1.20         # tight-junction assembly [1/d]
    p["dTJ"] = 1.10         # TJ disassembly by TNF/IL-1b/NO [1/d]
    p["kMUC"] = 1.00        # mucus production [1/d]
    p["dMUC"] = 0.90        # mucus degradation by neutrophil elastase [1/d]
    p["kIgAm"] = 1.60       # sIgA delivered per unit milk bioactivity [1/d]
    p["kIgAe"] = 0.10       # endogenous sIgA (maturation-limited) [1/d]
    p["dIgA"] = 1.50        # sIgA loss (proteolysis + stool) [1/d]
    p["kDEF"] = 0.80        # Paneth-cell defensin output [1/d]
    p["dDEF"] = 1.00        # defensin turnover [1/d]

    # ---- permeability / translocation --------------------------------------
    p["Pmin"] = 0.050       # paracellular permeability at intact barrier
    p["Pmax"] = 1.200       # permeability of denuded mucosa
    p["hP"] = 3.0           # steepness of permeability vs barrier loss
    p["fP"] = 0.15          # probiotic organisms' relative translocation weight

    # ---- innate signalling --------------------------------------------------
    p["Ktlr0"] = 60.00      # translocation flux for half-maximal TLR4 signal
    p["hHMO"] = 0.15        # HMO decoy-receptor raising of Ktlr
    p["KHMO"] = 0.40        # HMO pool for half of that effect [g/kg]
    p["hIgA"] = 0.10        # sIgA immune-exclusion raising of Ktlr
    p["thr_boost"] = 1.00   # THRESHOLD-MOVER lever (oral 2'-FL / anti-TLR4 /
                            # recombinant PAF-AH): multiplies Ktlr, nothing else
    p["ktlr_on"] = 24.0     # NF-kB activation [1/d]
    p["ktlr_off"] = 24.0    # NF-kB deactivation [1/d]
    p["iIL10"] = 0.90       # IL-10 damping of NF-kB
    p["iDEX"] = 0.70        # glucocorticoid damping of NF-kB (Emax)
    p["EC50DEX"] = 12.0     # dexamethasone EC50 [ug/L]

    # ---- cytokines ----------------------------------------------------------
    p["kIL1"] = 5.0;  p["dIL1"] = 6.0
    p["kTNF"] = 6.0;  p["dTNF"] = 8.0
    p["kIL8"] = 5.0;  p["dIL8"] = 4.0
    p["aTNF8"] = 0.6; p["aIL18"] = 0.6
    p["kIL10"] = 3.0; p["dIL10"] = 3.0
    p["kPAF"] = 4.0
    p["dPAF0"] = 1.2        # PAF-acetylhydrolase floor (preterm deficiency)
    p["dPAFm"] = 6.0        # maturation-dependent PAF-AH capacity
    p["kNEU"] = 3.0;  p["dNEU"] = 2.0
    p["NEUmax"] = 4.0       # saturation of the mucosal neutrophil pool
    p["kNOx"] = 2.5;  p["dNOx"] = 5.0
    p["wNO"] = 0.35         # weight of nitrosative injury in INJ
    p["wISCH"] = 0.90       # weight of ischaemia in INJ
    p["wNEC"] = 1.80        # necrotic tissue is ITSELF an inflammatory
                            # stimulus (DAMPs), so the lesion feeds the injury
                            # that made it -> a point of no return at
                            # NECa_crit = INJth / wNEC
    p["SDRcrit"] = 0.80     # O2 supply/demand ratio below which injury starts

    # ---- luminal ecology ----------------------------------------------------
    p["muB"] = 4.5          # pathobiont max growth [1/d]
    p["KsB"] = 0.20         # substrate for half-maximal B growth [g/kg]
    p["muC"] = 3.0          # commensal anaerobe max growth [1/d]
    p["KsC"] = 0.30
    p["muP"] = 3.4          # probiotic strain max growth [1/d]
    p["KsP"] = 0.30
    p["Ktot"] = 26.0        # total colonisable niche [1e9 CFU/g]
    p["alphaX"] = 0.55      # cross-taxon competition (<1 -> coexistence
                            # is possible; the two guilds are not interchangeable)
    p["binfantis"] = 1.0    # 1 = infant carries an HMO-utilising strain,
                            # 0.05 = does not (common in industrialised cohorts,
                            # in which HMOs are excreted largely unused)
    p["kwash"] = 1.50       # stool washout [1/d]
    p["kIgAB"] = 0.12       # sIgA-mediated clearance of B [1/d]
    p["kDEFB"] = 0.25       # defensin-mediated clearance of B [1/d]
    p["KpH"] = 10.0         # SCFA (mM) halving Enterobacteriaceae growth
    p["seedC_hm"] = 0.030   # daily Bifidobacterium seeding from milk [1e9/g/d]
    p["seedC_env"] = 0.004  # environmental anaerobe seeding
    p["seedB_env"] = 0.010  # environmental Enterobacteriaceae seeding
    p["muCH"] = 2.60        # commensal growth on HMO — a niche B cannot enter
    p["ySCFA"] = 32.0       # SCFA yield per unit anaerobe fermentation [mM/d]
    p["dSCFA"] = 2.2        # SCFA absorption / washout [1/d]
    p["yGAS"] = 0.55        # gas yield per unit fermentation
    p["dGAS"] = 1.60

    # ---- substrate handling -------------------------------------------------
    p["gPerMl"] = 0.070     # fermentable macronutrient per mL feed [g/mL]
    p["Smuc"] = 2.20        # HOST-derived substrate: mucin O-glycans and
                            # sloughed epithelium [g/kg/d].  Non-zero even on
                            # gut rest, which is why nil-by-mouth does not
                            # sterilise the lumen and why it cannot rescue a
                            # lesion that has already passed NECa_crit.
    p["Vabs"] = 14.0        # absorptive Vmax at E = 1 [g/kg/d]
    p["Kabs"] = 0.35        # Michaelis constant of absorption [g/kg]
    p["ktransit"] = 2.00    # substrate loss by transit [1/d]
    p["qferm"] = 0.075      # fermentation per unit biomass [g/kg/d per 1e9]
    p["Ksf"] = 0.30

    # ---- HMO ----------------------------------------------------------------
    p["hmoPerMl"] = 0.0130  # HMO in mother's own milk [g/mL]
    p["kHMOferm"] = 0.055   # HMO consumption per unit anaerobe [1/d per 1e9]
    p["dHMO"] = 1.50

    # ---- perfusion ----------------------------------------------------------
    p["kPERF"] = 6.0        # perfusion relaxation [1/d]
    p["aPAF"] = 0.45        # PAF-driven splanchnic vasoconstriction
    p["KaPAF"] = 0.60
    p["aCOX"] = 0.35        # loss of vasodilator prostanoid tone
    p["aNEC"] = 0.40        # perfusion lost inside the lesion
    p["aDEM"] = 0.40        # extra postprandial O2 demand at 200 mL/kg/d
    p["kPGE"] = 4.0         # prostanoid synthesis [1/d]
    p["EC50IND"] = 0.35     # indomethacin EC50 for COX inhibition [mg/L]
    p["EC50IBU"] = 12.0     # ibuprofen EC50 [mg/L]

    # ---- lesion -------------------------------------------------------------
    p["INJth"] = 0.35       # injury threshold for transmural necrosis
    p["knec"] = 4.00        # necrosis propagation [1/d]
    p["nNEC"] = 1.5
    p["krep"] = 0.30        # granulation / repair [1/d]
    p["NECrev"] = 0.10      # only mucosal/submucosal injury below this can
                            # re-epithelialise; transmural necrosis cannot
    p["dx_lag"] = 1.20      # suspicion -> film -> gut rest actually started [d]
    p["kPNE"] = 3.20        # pneumatosis formation [1/d]
    p["PbPNE"] = 0.20       # permeability floor below which no gas dissects
    p["dPNE"] = 0.55

    # ---- systemic -----------------------------------------------------------
    p["kLPS"] = 0.35; p["dLPS"] = 4.0; p["wNECLPS"] = 8.0
    p["PLT0"] = 250.0; p["kPLT"] = 0.35; p["kPLTc"] = 0.55
    p["kCRP"] = 22.0; p["dCRP"] = 0.60
    p["kLAC"] = 6.0; p["kLAC2"] = 5.0; p["dLAC"] = 3.0; p["LAC0"] = 1.0
    p["kNIN"] = 0.55; p["dNIN"] = 0.30
    p["kBIL"] = 0.70; p["dBIL"] = 0.22
    p["kcal_ml"] = 0.68     # kcal per mL feed
    p["kcal_tpn"] = 62.0    # kcal/kg/d provided by parenteral nutrition
    p["kcal_maint"] = 60.0  # maintenance kcal/kg/d
    p["kcal_g"] = 5.0       # kcal per g of weight gain
    p["kKIN"] = 0.06; p["dKIN"] = 0.20

    # ---- drug PK (preterm, PMA-scaled) --------------------------------------
    p["Vamp"] = 0.50; p["CLamp"] = 2.60      # L/kg, L/kg/d  (t1/2 ~ 3.2 h)
    p["Vgen"] = 0.50; p["CLgen"] = 1.05      # L/kg/d        (t1/2 ~ 8 h)
    p["Qgen"] = 0.35; p["Vgep"] = 0.30       # peripheral distribution
    p["Vmtz"] = 0.72; p["CLmtz"] = 0.36      # t1/2 ~ 33 h
    p["Vind"] = 0.36; p["CLind"] = 0.31      # t1/2 ~ 19 h
    p["Vibu"] = 0.32; p["CLibu"] = 0.09      # t1/2 ~ 59 h
    p["Vdex"] = 1.10; p["CLdex"] = 3.60
    # LUMINAL AVAILABILITY.  The microbiome is reshaped by the concentration
    # that reaches the colonic lumen, not the plasma concentration.  Ampicillin
    # is biliary/renally excreted and arrives in quantity; metronidazole
    # distributes almost freely; an intravenous aminoglycoside barely gets
    # there at all.  This asymmetry — and not any assumption about resistance —
    # is what makes empirical ampicillin + gentamicin destroy the anaerobes
    # while leaving the Enterobacteriaceae comparatively untouched.
    p["fl_amp"] = 0.45      # fraction of plasma concentration seen in lumen
    p["fl_gen"] = 0.05
    p["fl_mtz"] = 0.80
    p["EmaxB_amp"] = 1.10; p["EC50B_amp"] = 40.0   # partial (resistance)
    p["EmaxC_amp"] = 3.40; p["EC50C_amp"] = 6.0
    p["EmaxB_gen"] = 4.20; p["EC50B_gen"] = 2.0
    p["EmaxC_gen"] = 0.35; p["EC50C_gen"] = 4.0
    p["EmaxB_mtz"] = 0.10; p["EC50B_mtz"] = 8.0
    p["EmaxC_mtz"] = 4.60; p["EC50C_mtz"] = 4.0
    p["EmaxP_amp"] = 2.20; p["EmaxP_gen"] = 0.30; p["EmaxP_mtz"] = 1.30

    # ---- nutrition / regimen (set per scenario) -----------------------------
    p["fMOM"] = 1.0         # fraction of enteral volume: mother's own milk
    p["fDM"] = 0.0          # donor (pasteurised) milk
    p["fFORM"] = 0.0        # preterm formula
    p["fortify"] = 0.0      # bovine fortifier fraction of volume
    p["t_feed0"] = 2.0      # day enteral feeding starts
    p["feed_init"] = 20.0   # starting volume [mL/kg/d]
    p["feed_rate"] = 20.0   # advance rate [mL/kg/d per day]
    p["feed_max"] = 160.0   # target volume [mL/kg/d]
    p["prob_dose"] = 0.0    # probiotic daily dose [1e9 CFU/g-equivalent]
    p["prob_start"] = 2.0
    p["prob_stop"] = 60.0
    p["abx_start"] = 0.0
    p["abx_days"] = 0.0     # empirical ampicillin + gentamicin duration [d]
    p["mtz_days"] = 0.0     # metronidazole duration [d]
    p["ind_days"] = 0.0     # indomethacin course for PDA [d]
    p["ind_start"] = 3.0
    p["ibu_days"] = 0.0
    p["ibu_start"] = 3.0
    p["dex_days"] = 0.0
    p["npo_on_nec"] = 1.0   # clinical response: gut rest once Bell >= 2
    p["nec_abx"] = 1.0      # therapeutic triple antibiotics once Bell >= 2
    p["dysbiosis"] = 1.0    # multiplier on environmental B seeding
    # ---- precipitating insults ---------------------------------------------
    # Inside the bistable window a baby sits on the healthy branch and stays
    # there unless something transiently pushes E below E*.  Each entry is
    # (t_start, t_end, splanchnic perfusion depth, extra NF-kB drive) and
    # represents an apnoea/bradycardia or hypotensive episode, PDA diastolic
    # steal, a packed-red-cell transfusion, or a late-onset sepsis episode.
    p["_insults"] = []
    return p


def maturation(pma, p):
    x = (pma / p["PMA50"]) ** p["nMAT"]
    return x / (1.0 + x)


def milk_bioactivity(p):
    """Pasteurised donor milk keeps HMO but loses sIgA / lactoferrin / EGF."""
    return p["fMOM"] + 0.45 * p["fDM"]


# =============================================================================
# 2.  REGIMEN HELPERS
# =============================================================================
def feed_volume(t, p, npo_until, npo_from=0.0):
    if t < p["t_feed0"] or (npo_from <= t < npo_until):
        return 0.0
    v = p["feed_init"] + p["feed_rate"] * (t - max(p["t_feed0"], npo_until))
    return min(p["feed_max"], max(0.0, v))


# A fixed, reproducible insult schedule shared by the deterministic scenarios
# so that milk / antibiotic / probiotic arms are compared under identical
# provocation: two hypoperfusion episodes and two packed-cell transfusions.
# NOTE ON WHOM THE SCENARIOS DESCRIBE: these are a HIGH-RISK INDEX PATIENT,
# not a population average — a dysbiotic infant who receives one deep
# hypotensive episode on day 14, at the moment the luminal load has finished
# building.  A median infant on the same regimens does not tip at all (see the
# population sections), and comparing arms in a patient who never approaches
# the separatrix would show nothing.
STD_INSULTS = [(6.00, 6.30, 0.22, 0.000),      # apnoea-bradycardia / hypotension
               (10.00, 10.50, 0.16, 0.030),    # PRBC transfusion
               (14.00, 14.60, 0.40, 0.010),    # deep hypotensive episode
               (20.00, 20.50, 0.16, 0.030)]    # PRBC transfusion


def draw_insults(u, ga, tmax=45.0):
    """Per-patient insult history.  Episode rate rises steeply as GA falls."""
    ins = []
    rate = 0.09 * (1.0 + 0.45 * max(0.0, 30.0 - ga))       # episodes per day
    t = 0.5
    while True:
        t += -math.log(max(1e-12, u())) / rate
        if t >= tmax:
            break
        # cardiorespiratory instability is front-loaded in postnatal life:
        # thin the homogeneous process by exp(-t/20 d)
        if u() > math.exp(-t / 20.0):
            continue
        ins.append((t, t + 0.10 + 0.20 * u(), 0.10 + 0.20 * u(), 0.0))
    ntx = int(round(max(0.0, 3.2 - 0.28 * (ga - 24.0) + 0.8 * norm(u))))
    for _ in range(ntx):
        tt = 6.0 + (tmax - 12.0) * u()
        ins.append((tt - 0.6, tt, 0.10, 0.000))            # pre-transfusion anaemia
        ins.append((tt, tt + 0.5, 0.16, 0.030))            # transfusion itself
    return ins


def build_doses(p, tmax):
    """Discrete dose events applied between integration steps."""
    ev = []
    # ampicillin 100 mg/kg q12h, gentamicin 4 mg/kg q36h (preterm)
    if p["abx_days"] > 0:
        t = p["abx_start"]
        while t < p["abx_start"] + p["abx_days"]:
            ev.append((t, iAMP, 100.0 / p["Vamp"]))
            t += 0.5
        t = p["abx_start"]
        while t < p["abx_start"] + p["abx_days"]:
            ev.append((t, iGEN, 4.0 / p["Vgen"]))
            t += 1.5
    if p["mtz_days"] > 0:
        t = p["abx_start"]
        while t < p["abx_start"] + p["mtz_days"]:
            ev.append((t, iMTZ, 7.5 / p["Vmtz"]))
            t += 1.0
    if p["ind_days"] > 0:                      # 0.2 mg/kg q24h x 3
        t = p["ind_start"]
        while t < p["ind_start"] + p["ind_days"]:
            ev.append((t, iIND, 0.2 / p["Vind"]))
            t += 1.0
    if p["ibu_days"] > 0:                      # 10-5-5 mg/kg q24h
        t = p["ibu_start"]
        k = 0
        while t < p["ibu_start"] + p["ibu_days"]:
            ev.append((t, iIBU, (10.0 if k == 0 else 5.0) / p["Vibu"]))
            t += 1.0
            k += 1
    if p["dex_days"] > 0:
        t = 0.0
        while t < p["dex_days"]:
            ev.append((t, iDEX, 150.0 / p["Vdex"]))
            t += 0.5
    if p["prob_dose"] > 0:
        t = p["prob_start"]
        while t < min(p["prob_stop"], tmax):
            ev.append((t, iP, p["prob_dose"]))
            t += 1.0
    ev.sort(key=lambda z: z[0])
    return ev


# =============================================================================
# 3.  RIGHT-HAND SIDE
# =============================================================================
def rhs(t, y, p, ctl):
    (E, TJ, MUC, IgA, DEF, B, C, P, SUB, SCFA, GAS, HMO, TLRs, IL1, TNF, IL8,
     IL10, PAF, NEU, NOx, PERF, NECa, PNE, LPS, PLT, CRP, LAC, NIN, BIL, WT,
     KIN, AMP, GEN, GEP, MTZ, IND, IBU, DEX, PGE) = y

    d = [0.0] * NST
    pma = p["GA"] + t / 7.0
    MAT = maturation(pma, p)
    bio = milk_bioactivity(p)

    # ---------------- regimen -------------------------------------------------
    FEED = feed_volume(t, p, ctl["npo_until"], ctl["npo_from"])
    enteral = FEED / max(1e-9, p["feed_max"])

    # ---------------- barrier integrity and permeability ---------------------
    BI = E * (0.35 + 0.65 * TJ) * (0.55 + 0.45 * MUC)
    BI = min(1.0, max(0.0, BI))
    Pb = p["Pmin"] + (p["Pmax"] - p["Pmin"]) * (1.0 - BI) ** p["hP"]
    Jtr = Pb * (B + p["fP"] * P)

    # ---------------- precipitating insults ----------------------------------
    ins_perf = 0.0
    ins_nfkb = 0.0
    for (ta, tb, dp, bi) in p["_insults"]:
        if ta <= t < tb:
            ins_perf += dp
            ins_nfkb += bi

    # ---------------- innate signalling threshold ----------------------------
    Ktlr = (p["Ktlr0"] * p["thr_boost"]
            * (1.0 + p["hHMO"] * HMO / (HMO + p["KHMO"]))
            * (1.0 + p["hIgA"] * IgA))
    TLR4expr = p["TLR4max"] * (1.0 - p["phiTLR"] * MAT)
    dexE = p["iDEX"] * DEX / (DEX + p["EC50DEX"])
    d[iTLR] = (p["ktlr_on"] * TLR4expr * (1.0 - dexE)
               * (Jtr / (Jtr + Ktlr) + ins_nfkb)
               - p["ktlr_off"] * TLRs * (1.0 + p["iIL10"] * IL10))

    # ---------------- injury aggregate ---------------------------------------
    # splanchnic O2 supply/demand: feeding RAISES demand, so the same perfusion
    # buys less reserve on full feeds than on gut rest.
    demand = 1.0 + p["aDEM"] * FEED / 200.0
    ISCH = max(0.0, 1.0 - (PERF / demand) / p["SDRcrit"])
    INJ = TLRs + p["wNO"] * NOx + p["wISCH"] * ISCH + p["wNEC"] * NECa

    # ---------------- trophic drive and enterocyte balance -------------------
    Ftroph = ((p["troph_floor"] + (1.0 - p["troph_floor"]) * FEED / (FEED + p["Kfeed"]))
              * (1.0 + p["eEGF"] * bio)
              * (1.0 + p["eSCFA"] * SCFA / (SCFA + p["KsSCFA"])))
    block = 1.0 / (1.0 + (INJ / p["Ki"]) ** p["nKi"])
    d[iE] = (p["kE"] * (1.0 - E) * Ftroph * block
             - p["kA"] * INJ * E / (E + p["KEap"]))

    # ---------------- barrier accessories ------------------------------------
    d[iTJ] = (p["kTJ"] * MAT * (1.0 - TJ) * (1.0 + 0.5 * SCFA / (SCFA + p["KsSCFA"]))
              - p["dTJ"] * TJ * (0.45 * TNF + 0.45 * IL1 + 0.35 * NOx))
    d[iMUC] = (p["kMUC"] * MAT * (1.0 - MUC) * (1.0 + 0.4 * bio)
               - p["dMUC"] * MUC * (0.6 * NEU + 0.2 * TLRs))
    d[iIgA] = (p["kIgAm"] * bio * enteral + p["kIgAe"] * MAT
               - p["dIgA"] * IgA - 0.05 * IgA * B)
    d[iDEF] = p["kDEF"] * MAT * E * (1.0 - DEF) - p["dDEF"] * DEF

    # ---------------- antibiotic kill rates ----------------------------------
    aL = p["fl_amp"] * AMP          # luminal, not plasma, concentrations
    gL = p["fl_gen"] * GEN
    mL = p["fl_mtz"] * MTZ
    killB = (p["EmaxB_amp"] * aL / (aL + p["EC50B_amp"])
             + p["EmaxB_gen"] * gL / (gL + p["EC50B_gen"])
             + p["EmaxB_mtz"] * mL / (mL + p["EC50B_mtz"]))
    killC = (p["EmaxC_amp"] * aL / (aL + p["EC50C_amp"])
             + p["EmaxC_gen"] * gL / (gL + p["EC50C_gen"])
             + p["EmaxC_mtz"] * mL / (mL + p["EC50C_mtz"]))
    killP = (p["EmaxP_amp"] * aL / (aL + p["EC50C_amp"])
             + p["EmaxP_gen"] * gL / (gL + p["EC50B_gen"])
             + p["EmaxP_mtz"] * mL / (mL + p["EC50C_mtz"]))

    # ---------------- luminal ecology ---------------------------------------
    aX = p["alphaX"]
    occB = min(1.0, (B + aX * (C + P)) / p["Ktot"])
    occC = min(1.0, (C + P + aX * B) / p["Ktot"])
    fpH = 1.0 / (1.0 + (SCFA / p["KpH"]) ** 2)
    gB = p["muB"] * SUB / (SUB + p["KsB"]) * (1.0 - occB) * fpH
    gC = ((p["muC"] * SUB / (SUB + p["KsC"])
           + p["muCH"] * p["binfantis"] * HMO / (HMO + p["KHMO"]))
          * (1.0 - occC))
    gP = p["muP"] * SUB / (SUB + p["KsP"]) * (1.0 - occC)
    d[iB] = (B * (gB - p["kwash"] - killB
                  - p["kIgAB"] * IgA - p["kDEFB"] * DEF)
             + p["seedB_env"] * p["dysbiosis"])
    d[iC] = (C * (gC - p["kwash"] - killC)
             + p["seedC_hm"] * bio * enteral + p["seedC_env"])
    d[iP] = P * (gP - p["kwash"] - killP)

    # ---------------- substrate, SCFA, gas, HMO ------------------------------
    Sin = (FEED * p["gPerMl"] * (1.0 + 0.25 * p["fortify"])
           + p["Smuc"] * (0.5 + 0.5 * E))
    absorb = p["Vabs"] * E * SUB / (SUB + p["Kabs"])
    ferm = p["qferm"] * (B + C + P) * SUB / (SUB + p["Ksf"])
    d[iSUB] = Sin - absorb - ferm - p["ktransit"] * SUB
    d[iSCFA] = (p["ySCFA"] * p["qferm"] * (C + P) * SUB / (SUB + p["Ksf"])
                - p["dSCFA"] * SCFA)
    d[iGAS] = p["yGAS"] * ferm * (1.0 + 1.5 * (B / (B + 5.0))) - p["dGAS"] * GAS
    d[iHMO] = (FEED * p["hmoPerMl"] * p["fMOM"]
               + FEED * p["hmoPerMl"] * 0.95 * p["fDM"]
               - p["kHMOferm"] * (C + P) * HMO / (HMO + p["KHMO"])
               - p["dHMO"] * HMO)

    # ---------------- cytokines ---------------------------------------------
    dexf = 1.0 - dexE
    NEUs = NEU / (1.0 + NEU)          # saturating feed-forward from neutrophils
    PAFs = PAF / (1.0 + PAF)
    d[iIL1] = p["kIL1"] * dexf * (TLRs + 0.35 * NEUs) - p["dIL1"] * IL1
    d[iTNF] = p["kTNF"] * dexf * (TLRs + 0.20 * NEUs) - p["dTNF"] * TNF
    d[iIL8] = (p["kIL8"] * dexf * (TLRs + p["aTNF8"] * TNF + p["aIL18"] * IL1)
               - p["dIL8"] * IL8)
    d[iIL10] = p["kIL10"] * MAT * (TNF + IL1) - p["dIL10"] * IL10
    PAFAH = p["dPAF0"] + p["dPAFm"] * MAT
    d[iPAF] = p["kPAF"] * (TLRs + 0.5 * TNF) - PAFAH * PAF
    d[iNEU] = (p["kNEU"] * IL8 * (1.0 + 0.8 * PAFs)
               * (1.0 - NEU / p["NEUmax"]) - p["dNEU"] * NEU)
    d[iNOx] = p["kNOx"] * (TNF + IL1) - p["dNOx"] * NOx

    # ---------------- perfusion ---------------------------------------------
    coxI = (IND / (IND + p["EC50IND"])) + (IBU / (IBU + p["EC50IBU"]))
    coxI = min(1.0, coxI)
    d[iPGE] = p["kPGE"] * ((1.0 - coxI) - PGE)
    target = (1.0
              - p["aPAF"] * PAF / (PAF + p["KaPAF"])
              - p["aCOX"] * (1.0 - PGE)
              - p["aNEC"] * NECa
              - ins_perf)
    target = min(1.0, max(0.15, target))
    d[iPERF] = p["kPERF"] * (target - PERF)

    # ---------------- lesion -------------------------------------------------
    over = max(0.0, INJ - p["INJth"])
    d[iNEC] = (p["knec"] * (1.0 - NECa) * over ** p["nNEC"]
               - p["krep"] * min(NECa, p["NECrev"]) * E
               * (1.0 if INJ < p["INJth"] else 0.0))
    # gas dissects into the bowel wall only through an already-breached
    # mucosa: intramural gas is a CONSEQUENCE of the barrier failing, so the
    # permeability must clear a floor before any pneumatosis can form.
    d[iPNE] = (p["kPNE"] * GAS * max(0.0, Pb - p["PbPNE"]) * (1.0 - PNE)
               - p["dPNE"] * PNE)

    # ---------------- systemic ----------------------------------------------
    d[iLPS] = p["kLPS"] * Jtr * (1.0 + p["wNECLPS"] * NECa) - p["dLPS"] * LPS
    d[iPLT] = (p["kPLT"] * (p["PLT0"] - PLT)
               - p["kPLTc"] * PLT * (0.5 * LPS + 3.0 * NECa))
    d[iCRP] = p["kCRP"] * (0.6 * IL1 + 0.4 * TNF) - p["dCRP"] * CRP
    d[iLAC] = (p["kLAC"] * max(0.0, 1.0 - PERF) ** 2 + p["kLAC2"] * NECa
               - p["dLAC"] * (LAC - p["LAC0"]))
    d[iNIN] = p["kNIN"] * (0.4 * TNF + 0.4 * IL1 + 0.5 * LPS) - p["dNIN"] * NIN

    # ---------------- liver / growth / kidney -------------------------------
    tpn = max(0.0, 1.0 - enteral)
    d[iBIL] = p["kBIL"] * tpn * (1.0 + 0.5 * LPS) - p["dBIL"] * BIL * (1.0 + enteral)
    kcal_ent = (absorb / max(1e-9, p["gPerMl"])) * p["kcal_ml"]
    kcal_in = kcal_ent + p["kcal_tpn"] * tpn
    d[iWT] = max(-0.004, (kcal_in - p["kcal_maint"] * (1.0 + 0.30 * NECa))
                 / p["kcal_g"] / 1000.0 * WT / max(0.3, WT))
    d[iKIN] = p["kKIN"] * GEP - p["dKIN"] * KIN

    # ---------------- drug PK ------------------------------------------------
    sc = (pma / 40.0) ** 1.3          # PMA scaling of clearance
    d[iAMP] = -p["CLamp"] * sc / p["Vamp"] * AMP
    d[iGEN] = (-p["CLgen"] * sc / p["Vgen"] * GEN
               - p["Qgen"] / p["Vgen"] * (GEN - GEP))
    d[iGEP] = p["Qgen"] / p["Vgep"] * (GEN - GEP)
    d[iMTZ] = -p["CLmtz"] * sc / p["Vmtz"] * MTZ
    d[iIND] = -p["CLind"] * sc / p["Vind"] * IND
    d[iIBU] = -p["CLibu"] * sc / p["Vibu"] * IBU
    d[iDEX] = -p["CLdex"] * sc / p["Vdex"] * DEX
    return d


# =============================================================================
# 4.  INITIAL CONDITIONS AND INTEGRATION
# =============================================================================
def init_state(p):
    pma = p["GA"]
    MAT = maturation(pma, p)
    y = [0.0] * NST
    y[iE] = min(0.985, 0.70 + 0.026 * (p["GA"] - 24.0))
    y[iTJ] = 0.30 + 0.65 * MAT
    y[iMUC] = 0.30 + 0.60 * MAT
    y[iIgA] = 0.02
    y[iDEF] = 0.15 + 0.55 * MAT
    y[iB] = 0.010
    y[iC] = 0.005
    y[iP] = 0.0
    y[iSUB] = 0.02
    y[iSCFA] = 0.5
    y[iGAS] = 0.05
    y[iHMO] = 0.0
    y[iPERF] = 1.0
    y[iPLT] = p["PLT0"]
    y[iCRP] = 1.0
    y[iLAC] = p["LAC0"]
    y[iBIL] = 0.4
    y[iWT] = p["BW"]
    y[iPGE] = 1.0
    return y


def clamp(y):
    for i in range(NST):
        if y[i] < 0.0:
            y[i] = 0.0
    if y[iE] > 1.0:
        y[iE] = 1.0
    if y[iE] < 1e-4:
        y[iE] = 1e-4
    for i in (iTJ, iMUC, iIgA, iDEF, iNEC, iPNE, iPGE):
        if y[i] > 1.0:
            y[i] = 1.0
    if y[iPERF] > 1.0:
        y[iPERF] = 1.0
    return y


def bell_stage(y, p):
    """Modified Bell staging read out of the state vector."""
    NECa, PNE, PLT, LAC, CRP, GAS = y[iNEC], y[iPNE], y[iPLT], y[iLAC], y[iCRP], y[iGAS]
    if NECa > 0.28 or (PLT < 100.0 and LAC > 4.0):
        return 3
    if PNE > 0.30 or NECa > 0.10:
        return 2
    if GAS > 0.55 or CRP > 12.0 or NECa > 0.01:
        return 1
    return 0


def simulate(p, tmax=56.0, dt=0.01, record=False, rec_every=10):
    y = init_state(p)
    ev = build_doses(p, tmax)
    ei = 0
    ctl = {"npo_until": -1.0, "nec_abx_until": -1.0, "nec_day": None,
           "peak_stage": 0, "surg": 0, "stage2_day": None, "stage3_day": None,
           "treated": False, "npo_from": 1e9, "nec_abx_from": 1e9}
    rec = {"t": []}
    if record:
        for n in NAMES:
            rec[n] = []
        rec["BELL"] = []
        rec["Jtr"] = []
        rec["INJ"] = []
        rec["Pb"] = []

    n = int(round(tmax / dt))
    marker_days = {}
    for k in range(n):
        t = k * dt
        # ---- discrete doses -------------------------------------------------
        while ei < len(ev) and ev[ei][0] <= t + 1e-12:
            y[ev[ei][1]] += ev[ei][2]
            ei += 1
        # ---- therapeutic response to diagnosed NEC --------------------------
        st = bell_stage(y, p)
        if st >= 2 and not ctl["treated"]:
            # ONE course of gut rest + therapeutic triple antibiotics, started
            # at the first radiographic diagnosis (not re-triggered every step)
            ctl["treated"] = True
            ctl["nec_day"] = t
            # the response is NOT instantaneous: suspicion -> abdominal film ->
            # nil-by-mouth actually in force.  The lesion keeps growing during
            # that lag, and how deep the collapse is during the lag is exactly
            # what separates medical from surgical NEC in this model.
            if p["npo_on_nec"] > 0.5:
                ctl["npo_until"] = t + p["dx_lag"] + 10.0
                ctl["npo_from"] = t + p["dx_lag"]
            if p["nec_abx"] > 0.5:
                ctl["nec_abx_until"] = t + p["dx_lag"] + 10.0
                ctl["nec_abx_from"] = t + p["dx_lag"]
        if st >= 2 and ctl["stage2_day"] is None:
            ctl["stage2_day"] = t
        if st >= 3 and ctl["stage3_day"] is None:
            ctl["stage3_day"] = t
        if st > ctl["peak_stage"]:
            ctl["peak_stage"] = st
        if ctl["nec_abx_from"] <= t < ctl["nec_abx_until"]:
            # therapeutic amp + gent + metronidazole while Bell >= 2:
            # entered as the continuous-rate equivalent of q12h / q36h / q24h
            y[iAMP] += 100.0 / p["Vamp"] * dt / 0.5
            y[iGEN] += 4.0 / p["Vgen"] * dt / 1.5
            y[iMTZ] += 7.5 / p["Vmtz"] * dt / 1.0
        # ---- markers --------------------------------------------------------
        for nm, idx, thr in (("CRP12", iCRP, 12.0), ("PLT100", iPLT, -100.0),
                             ("LAC4", iLAC, 4.0), ("PNEU30", iPNE, 0.30)):
            if nm not in marker_days:
                v = y[idx]
                hit = (v > thr) if thr > 0 else (v < -thr)
                if hit:
                    marker_days[nm] = t
        if record and k % rec_every == 0:
            rec["t"].append(t)
            for nm in NAMES:
                rec[nm].append(y[IX[nm]])
            rec["BELL"].append(st)
            BI = y[iE] * (0.35 + 0.65 * y[iTJ]) * (0.55 + 0.45 * y[iMUC])
            Pb = p["Pmin"] + (p["Pmax"] - p["Pmin"]) * (1.0 - min(1.0, BI)) ** p["hP"]
            rec["Pb"].append(Pb)
            rec["Jtr"].append(Pb * (y[iB] + p["fP"] * y[iP]))
            MAT = maturation(p["GA"] + t / 7.0, p)
            fv = feed_volume(t, p, ctl["npo_until"], ctl["npo_from"])
            dem = 1.0 + p["aDEM"] * fv / 200.0
            rec["INJ"].append(
                y[iTLR] + p["wNO"] * y[iNOx]
                + p["wISCH"] * max(0.0, 1.0 - (y[iPERF] / dem) / p["SDRcrit"]))
        # ---- RK4 ------------------------------------------------------------
        k1 = rhs(t, y, p, ctl)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2 = rhs(t + 0.5 * dt, y2, p, ctl)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3 = rhs(t + 0.5 * dt, y3, p, ctl)
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4 = rhs(t + dt, y4, p, ctl)
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        y = clamp(y)

    # ---- endpoints ----------------------------------------------------------
    surg = 1 if ctl["peak_stage"] >= 3 else 0
    ga_pen = max(0.0, (30.0 - p["GA"])) * 0.035
    mort = 0.0
    if ctl["peak_stage"] >= 3:
        mort = min(0.85, 0.35 + ga_pen + 0.4 * y[iNEC])
    elif ctl["peak_stage"] == 2:
        mort = min(0.35, 0.06 + ga_pen * 0.5)
    out = {
        "peak_bell": ctl["peak_stage"],
        "nec": 1 if ctl["peak_stage"] >= 2 else 0,
        "onset_day": ctl["stage2_day"],
        "stage3_day": ctl["stage3_day"],
        "surgical": surg,
        "mortality": mort,
        "E_end": y[iE], "B_end": y[iB], "C_end": y[iC], "SCFA_end": y[iSCFA],
        "NECa_end": y[iNEC], "PLT_min_proxy": y[iPLT], "BILI_end": y[iBIL],
        "WT_end": y[iWT], "NIN_end": y[iNIN], "KIN_end": y[iKIN],
        "markers": marker_days,
    }
    if record:
        out["rec"] = rec
    return out


# =============================================================================
# 5.  REDUCED 1-D FIELD AND BIFURCATION ANALYSIS
# =============================================================================
def quasistatic_cascade(Jtr, p, ctx, niter=400):
    """
    Steady state of the fast block (TLR4s, cytokines, PAF, neutrophils, NO,
    perfusion) at a held translocation flux Jtr.  This is EXACTLY the ODE
    steady state of those 9 equations, obtained by damped fixed-point
    iteration; it is not an independent approximation.
    """
    Ktlr = (p["Ktlr0"] * p["thr_boost"]
            * (1.0 + p["hHMO"] * ctx["HMO"] / (ctx["HMO"] + p["KHMO"]))
            * (1.0 + p["hIgA"] * ctx["IgA"]))
    drive = ctx["TLR4expr"] * Jtr / (Jtr + Ktlr)
    PAFAH = p["dPAF0"] + p["dPAFm"] * ctx["MAT"]
    x = drive
    TNF = IL1 = IL8 = NEU = PAF = NOx = IL10 = 0.0
    PERF = 1.0
    w = 0.25
    for _ in range(niter):
        x_new = drive / (1.0 + p["iIL10"] * IL10)
        NEUs = NEU / (1.0 + NEU)
        PAFs = PAF / (1.0 + PAF)
        TNF_n = p["kTNF"] * (x_new + 0.20 * NEUs) / p["dTNF"]
        IL1_n = p["kIL1"] * (x_new + 0.35 * NEUs) / p["dIL1"]
        IL8_n = (p["kIL8"] * (x_new + p["aTNF8"] * TNF + p["aIL18"] * IL1)
                 / p["dIL8"])
        IL10_n = p["kIL10"] * ctx["MAT"] * (TNF + IL1) / p["dIL10"]
        PAF_n = p["kPAF"] * (x_new + 0.5 * TNF) / PAFAH
        # dNEU/dt = 0  ->  a*(1-NEU/NEUmax) = dNEU*NEU
        a = p["kNEU"] * IL8 * (1.0 + 0.8 * PAFs)
        NEU_n = a / (p["dNEU"] + a / p["NEUmax"])
        NOx_n = p["kNOx"] * (TNF + IL1) / p["dNOx"]
        tgt = (1.0 - p["aPAF"] * PAF / (PAF + p["KaPAF"]) - p["aCOX"] * ctx["coxI"])
        PERF_n = min(1.0, max(0.15, tgt))
        x = x + w * (x_new - x);      TNF = TNF + w * (TNF_n - TNF)
        IL1 = IL1 + w * (IL1_n - IL1); IL8 = IL8 + w * (IL8_n - IL8)
        IL10 = IL10 + w * (IL10_n - IL10); PAF = PAF + w * (PAF_n - PAF)
        NEU = NEU + w * (NEU_n - NEU); NOx = NOx + w * (NOx_n - NOx)
        PERF = PERF + w * (PERF_n - PERF)
    ISCH = max(0.0, 1.0 - (PERF / ctx["demand"]) / p["SDRcrit"])
    INJ = x + p["wNO"] * NOx + p["wISCH"] * ISCH
    return {"TLRs": x, "TNF": TNF, "IL1": IL1, "IL8": IL8, "IL10": IL10,
            "PAF": PAF, "NEU": NEU, "NOx": NOx, "PERF": PERF, "ISCH": ISCH,
            "INJ": INJ, "Ktlr": Ktlr}


def build_injury_table(p, ctx, jmax=400.0, n=500):
    """
    INJ depends on E and B only through Jtr, so tabulate INJ(Jtr) once per
    context and interpolate.  Log-spaced knots; INJ is monotone in Jtr.
    """
    js = [0.0] + [jmax * ((i / n) ** 2.5) for i in range(1, n + 1)]
    vs = [quasistatic_cascade(j, p, ctx, niter=250)["INJ"] for j in js]
    ctx["_jt"], ctx["_iv"], ctx["_jmax"] = js, vs, jmax
    return ctx


def injury_of(Jtr, ctx):
    js, vs = ctx["_jt"], ctx["_iv"]
    if Jtr <= 0.0:
        return vs[0]
    if Jtr >= js[-1]:
        return vs[-1]
    lo, hi = 0, len(js) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if js[mid] <= Jtr:
            lo = mid
        else:
            hi = mid
    f = (Jtr - js[lo]) / (js[hi] - js[lo])
    return vs[lo] + f * (vs[hi] - vs[lo])


def reduced_field(E, B, p, ctx):
    """g(E) for the quasi-static reduction (TJ = MUC = 1, B / HMO / IgA held)."""
    Pb = p["Pmin"] + (p["Pmax"] - p["Pmin"]) * (1.0 - E) ** p["hP"]
    INJ = injury_of(Pb * B, ctx)
    block = 1.0 / (1.0 + (INJ / p["Ki"]) ** p["nKi"])
    return (p["kE"] * (1.0 - E) * ctx["Ftroph"] * block
            - p["kA"] * INJ * E / (E + p["KEap"]))


def roots_of_field(B, p, ctx, ngrid=4000):
    """All roots of g(E) on (0,1) with their stability."""
    xs = [i / ngrid for i in range(1, ngrid)]
    vals = [reduced_field(x, B, p, ctx) for x in xs]
    out = []
    for i in range(len(xs) - 1):
        a, b = vals[i], vals[i + 1]
        if a == 0.0:
            out.append(xs[i])
        elif a * b < 0.0:
            lo, hi = xs[i], xs[i + 1]
            for _ in range(60):
                m = 0.5 * (lo + hi)
                if reduced_field(lo, B, p, ctx) * reduced_field(m, B, p, ctx) <= 0:
                    hi = m
                else:
                    lo = m
            out.append(0.5 * (lo + hi))
    res = []
    for r in out:
        h = 1e-5
        slope = (reduced_field(min(1 - 1e-9, r + h), B, p, ctx)
                 - reduced_field(max(1e-9, r - h), B, p, ctx)) / (2 * h)
        res.append((r, "stable" if slope < 0 else "unstable"))
    return res


def context(p, milk="MOM", ga=28.0, feed=150.0, probiotic=False):
    """Quasi-static context for the reduced field at a given regimen."""
    q = dict(p)
    q["GA"] = ga
    if milk == "MOM":
        q["fMOM"], q["fDM"], q["fFORM"] = 1.0, 0.0, 0.0
    elif milk == "DM":
        q["fMOM"], q["fDM"], q["fFORM"] = 0.0, 1.0, 0.0
    else:
        q["fMOM"], q["fDM"], q["fFORM"] = 0.0, 0.0, 1.0
    bio = milk_bioactivity(q)
    MAT = maturation(ga + 2.0, q)          # ~2 weeks postnatal, the risk window
    # quasi-static HMO pool
    hmo_in = feed * q["hmoPerMl"] * (q["fMOM"] + 0.95 * q["fDM"])
    HMO = hmo_in / (q["dHMO"] + 0.25)
    IgA = (q["kIgAm"] * bio * feed / q["feed_max"] + q["kIgAe"] * MAT) / q["dIgA"]
    IgA = min(1.0, IgA)
    SCFA = 12.0 if (bio > 0.5 or probiotic) else 4.5
    Ftroph = ((q["troph_floor"] + (1 - q["troph_floor"]) * feed / (feed + q["Kfeed"]))
              * (1.0 + q["eEGF"] * bio)
              * (1.0 + q["eSCFA"] * SCFA / (SCFA + q["KsSCFA"])))
    TLR4expr = q["TLR4max"] * (1.0 - q["phiTLR"] * MAT)
    ctx = {"HMO": HMO, "IgA": IgA, "SCFA": SCFA, "Ftroph": Ftroph,
           "TLR4expr": TLR4expr, "bio": bio, "MAT": MAT,
           "demand": 1.0 + q["aDEM"] * feed / 200.0, "coxI": 0.0,
           "feed": feed, "milk": milk, "GA": ga, "probiotic": probiotic}
    return build_injury_table(q, ctx)


def _nroots(B, p, ctx, ngrid=800):
    rs = roots_of_field(B, p, ctx, ngrid=ngrid)
    return len(rs), rs


def bifurcation_loads(p, ctx, Bmax=240.0, nB=240):
    """
    B_lo = saddle-node at which bistability appears (NEC becomes possible)
    B_hi = saddle-node at which the healthy root is annihilated (inevitable)
    Coarse scan then bisection on the root count.
    """
    grid = [Bmax * i / nB for i in range(1, nB + 1)]
    counts = []
    for B in grid:
        k, rs = _nroots(B, p, ctx)
        counts.append((B, k, rs))

    def refine(Ba, Bb, want):
        for _ in range(34):
            m = 0.5 * (Ba + Bb)
            k, _rs = _nroots(m, p, ctx, ngrid=2000)
            if k >= want:
                Bb = m
            else:
                Ba = m
        return 0.5 * (Ba + Bb)

    B_lo = B_hi = None
    for i in range(1, len(counts)):
        Bprev, kprev, _ = counts[i - 1]
        B, k, rs = counts[i]
        if B_lo is None and k >= 3 and kprev == 1:
            B_lo = refine(Bprev, B, 3)
        if B_lo is not None and B_hi is None and k == 1 and rs[0][0] < 0.5:
            # healthy branch gone: bisect on "highest root < 0.5"
            Ba, Bb = Bprev, B
            for _ in range(34):
                m = 0.5 * (Ba + Bb)
                _k, _rs = _nroots(m, p, ctx, ngrid=2000)
                if max(r for r, s in _rs) < 0.5:
                    Bb = m
                else:
                    Ba = m
            B_hi = 0.5 * (Ba + Bb)
            break
    return B_lo, B_hi


# =============================================================================
# 6.  SCENARIOS
# =============================================================================
def scenario_params(name):
    p = base_params()
    if name == "S1_28wk_MOM":
        p.update(GA=28.0, BW=1.05, fMOM=1.0, fDM=0.0, fFORM=0.0,
                 feed_rate=20.0, prob_dose=0.0, abx_days=0.0, dysbiosis=1.5)
    elif name == "S2_28wk_formula":
        p.update(GA=28.0, BW=1.05, fMOM=0.0, fDM=0.0, fFORM=1.0,
                 feed_rate=20.0, prob_dose=0.0, abx_days=0.0, dysbiosis=1.5)
    elif name == "S3_28wk_formula_abx7":
        p.update(GA=28.0, BW=1.05, fMOM=0.0, fFORM=1.0,
                 feed_rate=20.0, abx_days=7.0, dysbiosis=1.5)
    elif name == "S4_28wk_formula_probiotic":
        p.update(GA=28.0, BW=1.05, fMOM=0.0, fFORM=1.0,
                 feed_rate=20.0, prob_dose=1.10, dysbiosis=1.5)
    elif name == "S5_25wk_formula_abx7_indo":
        p.update(GA=25.0, BW=0.68, fMOM=0.0, fFORM=1.0, feed_rate=20.0,
                 abx_days=7.0, ind_days=3.0, dysbiosis=1.8)
    elif name == "S6_25wk_MOM_probiotic":
        p.update(GA=25.0, BW=0.68, fMOM=1.0, fFORM=0.0, feed_rate=20.0,
                 prob_dose=1.10, abx_days=2.0, dysbiosis=1.8)
    elif name == "S7_32wk_formula":
        p.update(GA=32.0, BW=1.75, fMOM=0.0, fFORM=1.0, feed_rate=30.0,
                 abx_days=0.0, dysbiosis=1.5)
    elif name == "S8_28wk_donor_milk":
        p.update(GA=28.0, BW=1.05, fMOM=0.0, fDM=1.0, fFORM=0.0,
                 feed_rate=20.0, dysbiosis=1.5)
    elif name == "S9_28wk_NPO_TPN":
        p.update(GA=28.0, BW=1.05, fMOM=1.0, feed_rate=0.0, feed_init=0.0,
                 feed_max=1e-6, t_feed0=99.0, dysbiosis=1.5)
    elif name == "S11_28wk_formula_severe":
        # a deliberately provoked case: dysbiotic, formula-fed, prolonged
        # empirical antibiotics, indomethacin, and a deep hypotensive episode
        # landing while the load already sits inside the bistable window
        p.update(GA=27.0, BW=0.92, fMOM=0.0, fFORM=1.0, feed_rate=30.0,
                 abx_days=7.0, mtz_days=5.0, ind_days=3.0, dysbiosis=1.8)
    elif name == "S10_28wk_formula_abx10_mtz":
        p.update(GA=28.0, BW=1.05, fMOM=0.0, fFORM=1.0, feed_rate=20.0,
                 abx_days=10.0, mtz_days=10.0, dysbiosis=1.5)
    else:
        raise KeyError(name)
    p["_insults"] = list(STD_INSULTS)
    if name == "S11_28wk_formula_severe":
        p["npo_on_nec"] = 0.0        # the "what if nobody stops it" control
        p["nec_abx"] = 0.0
    return p


SCENARIOS = [
    ("S1_28wk_MOM", "28주 · 모유 · 항생제 없음"),
    ("S2_28wk_formula", "28주 · 분유 · 항생제 없음"),
    ("S3_28wk_formula_abx7", "28주 · 분유 · 경험적 항생제 7일"),
    ("S4_28wk_formula_probiotic", "28주 · 분유 · 프로바이오틱스"),
    ("S5_25wk_formula_abx7_indo", "25주 · 분유 · 항생제 7일 + 인도메타신"),
    ("S6_25wk_MOM_probiotic", "25주 · 모유 + 프로바이오틱스 · 항생제 2일"),
    ("S7_32wk_formula", "32주 · 분유"),
    ("S8_28wk_donor_milk", "28주 · 기증모유(저온살균)"),
    ("S9_28wk_NPO_TPN", "28주 · 완전 금식 + TPN"),
    ("S10_28wk_formula_abx10_mtz", "28주 · 분유 · 항생제 10일 + 메트로니다졸"),
    ("S11_28wk_formula_severe", "27주 · 분유 · 항생제 7일 + 인도메타신 + 심한 저혈압 삽화"),
]


# =============================================================================
# 7.  POPULATION SIMULATION
# =============================================================================
def lcg(seed):
    s = [seed & 0xFFFFFFFF]

    def u():
        s[0] = (1103515245 * s[0] + 12345) & 0x7FFFFFFF
        return s[0] / 0x7FFFFFFF
    return u


def norm(u):
    a = max(1e-12, u())
    b = u()
    return math.sqrt(-2.0 * math.log(a)) * math.cos(2.0 * math.pi * b)


def population(n, arm, seed=20260805, tmax=45.0, dt=0.02):
    """arm = dict of overrides; GA drawn from a VLBW-like distribution."""
    u = lcg(seed)
    res = []
    for i in range(n):
        p = base_params()
        ga = 24.0 + 8.0 * (u() ** 0.75)          # 24-32 wk, skewed to lower GA
        p["GA"] = ga
        p["BW"] = max(0.42, 0.10 * (ga - 22.0) + 0.35 + 0.09 * norm(u))
        # inter-individual variability (log-normal, CV ~ 25 %)
        for key, cv in (("Ktlr0", 0.22), ("kA", 0.20), ("kE", 0.18),
                        ("muB", 0.20), ("Pmax", 0.15), ("seedB_env", 0.45),
                        ("kPAF", 0.25), ("dPAF0", 0.25),
                        ("seedC_hm", 1.00), ("seedC_env", 1.00)):
            p[key] *= math.exp(cv * norm(u))
        # which basin the gut lands in is decided in the first days, and one
        # of the deciding factors is whether an HMO-utilising strain is present
        p["binfantis"] = 1.0 if u() < 0.35 else 0.05
        p["dysbiosis"] = math.exp(0.35 * norm(u))
        p.update(arm)
        p["_insults"] = draw_insults(u, ga, tmax=tmax)
        r = simulate(p, tmax=tmax, dt=dt)
        r["GA"] = ga
        res.append(r)
    return res


def summarise(res):
    n = len(res)
    nec = sum(r["nec"] for r in res)
    surg = sum(r["surgical"] for r in res)
    died = sum(r["mortality"] for r in res)
    onsets = [r["onset_day"] for r in res if r["onset_day"] is not None]
    return {
        "n": n,
        "nec_n": nec,
        "nec_pct": 100.0 * nec / n,
        "surgical_pct": 100.0 * surg / n,
        "expected_deaths_pct": 100.0 * died / n,
        "median_onset_day": (sorted(onsets)[len(onsets) // 2] if onsets else None),
        "mean_onset_day": (sum(onsets) / len(onsets) if onsets else None),
    }


# =============================================================================
# 8.  REPORT
# =============================================================================
def hr(ch="="):
    return ch * 79


def main():
    out = {}
    P = base_params()
    W = sys.stdout.write

    W(hr() + "\n")
    W("NECROTISING ENTEROCOLITIS (NEC) — QSP REFERENCE MODEL, EXECUTED OUTPUT\n")
    W("신생아 괴사성 장염 QSP 참조 모델 — 실제 적분 결과\n")
    W(hr() + "\n")
    W("states = %d ODEs | integrator = RK4, fixed step | pure Python stdlib\n\n"
      % NST)

    # -------------------------------------------------------------------------
    # 8.1  the closed loop: two critical loads
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("1. TWO CRITICAL PATHOBIONT LOADS (saddle-node bifurcations of g(E))\n")
    W("   임계 병원균 부하는 하나가 아니라 둘이다\n")
    W(hr("-") + "\n")
    W("   g(E) = kE(1-E)·Ftroph/(1+(INJ/Ki)^2) - kA·E·INJ ,  INJ = TLR4s(E;B)\n")
    W("   B_lo = NEC가 '가능'해지는 부하  (bistability appears)\n")
    W("   B_hi = NEC가 '필연'이 되는 부하  (healthy root annihilates)\n")
    W("   niche ceiling Ktot = %.0f x1e9 CFU/g — B cannot physically exceed it\n\n"
      % P["Ktot"])
    W("   %-6s %-9s %-8s | %8s %8s %8s %8s | %8s %8s %-9s\n"
      % ("GA", "milk", "probio", "Ftroph", "Ktlr", "TLR4ex", "E_healthy",
         "B_lo", "B_hi", "reachable?"))
    bif_rows = []
    for ga in (25.0, 28.0, 32.0):
        for milk in ("FORM", "DM", "MOM"):
            for probio in (False, True):
                ctx = context(P, milk=milk, ga=ga, feed=150.0, probiotic=probio)
                Ktlr = (P["Ktlr0"] * (1 + P["hHMO"] * ctx["HMO"] / (ctx["HMO"] + P["KHMO"]))
                        * (1 + P["hIgA"] * ctx["IgA"]))
                blo, bhi = bifurcation_loads(P, ctx)
                rs = roots_of_field(0.5, P, ctx)
                eh = max(r for r, s in rs) if rs else float("nan")
                reach = ("YES" if (blo is not None and blo < P["Ktot"])
                         else "no (>Ktot)")
                W("   %-6.0f %-9s %-8s | %8.3f %8.2f %8.3f %9.3f | %8s %8s %-9s\n"
                  % (ga, milk, "yes" if probio else "no", ctx["Ftroph"], Ktlr,
                     ctx["TLR4expr"], eh,
                     ("%.1f" % blo) if blo else "none",
                     ("%.1f" % bhi) if bhi else "none", reach))
                bif_rows.append({"GA": ga, "milk": milk, "probiotic": probio,
                                 "Ftroph": ctx["Ftroph"], "Ktlr": Ktlr,
                                 "TLR4expr": ctx["TLR4expr"], "E_healthy": eh,
                                 "B_lo": blo, "B_hi": bhi,
                                 "reachable": blo is not None and blo < P["Ktot"]})
    out["bifurcation"] = bif_rows
    W("\n   READ: 부하가 B_lo와 B_hi 사이일 때 결과는 부하가 아니라 '초기조건'이\n")
    W("   결정한다 — 같은 병동·같은 수유에서 두 아이의 운명이 갈리는 이유.\n\n")

    # -------------------------------------------------------------------------
    # 8.2  the separatrix and the transit-time logarithm
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("2. THE SEPARATRIX E* AND WHY BELL II -> III IS ALWAYS ~A DAY\n")
    W("   분수령 E*와 전격성 악화 시간의 로그 의존성\n")
    W(hr("-") + "\n")
    ctx = context(P, milk="FORM", ga=28.0, feed=150.0)
    Btest = 18.0
    rs = roots_of_field(Btest, P, ctx)
    W("   28 wk · formula · B = %.0f  ->  roots of g(E):\n" % Btest)
    for r, s in rs:
        W("     E = %.4f   %s\n" % (r, s))
    sep = [r for r, s in rs if s == "unstable"]
    Estar = sep[0] if sep else None
    Elow = min(r for r, s in rs) if rs else None
    out["separatrix"] = {"B": Btest, "roots": rs, "E_star": Estar, "E_low": Elow}
    if Estar and Elow is not None and Elow < Estar:
        h = 1e-4
        lam = (reduced_field(Estar + h, Btest, P, ctx)
               - reduced_field(Estar - h, Btest, P, ctx)) / (2 * h)
        span = Estar - Elow
        target = Elow + 0.05 * span     # "arrived at the necrotic attractor"
        W("\n   necrotic attractor E_low = %.4f ; separatrix E* = %.4f\n"
          % (Elow, Estar))
        W("   local eigenvalue at E*:  lambda = %+.3f /d  (positive = repelling)\n"
          % lam)
        W("   transit time from E* - delta down to E_low + 5% of the gap:\n")
        W("     %-12s %-16s %-16s\n"
          % ("delta", "ln-law pred (h)", "integrated (h)"))
        rows = []
        for delta in (0.0005, 0.002, 0.008, 0.030, 0.100):
            if delta >= span:
                continue
            # linear-field prediction near E*:  t = (1/lambda)·ln(span/delta)
            pred = (1.0 / abs(lam)) * math.log(span / delta) * 24.0
            Ecur, tt, dtl = Estar - delta, 0.0, 0.0002
            while Ecur > target and tt < 60.0:
                Ecur += dtl * reduced_field(Ecur, Btest, P, ctx)
                tt += dtl
            rows.append({"delta": delta, "pred_h": pred, "sim_h": tt * 24.0})
            W("     %-12.4f %-16.1f %-16.1f\n" % (delta, pred, tt * 24.0))
        out["separatrix"]["transit"] = rows
        if len(rows) >= 2 and rows[-1]["sim_h"] > 0:
            dr = rows[0]["delta"] / rows[-1]["delta"]
            tr = rows[0]["sim_h"] / rows[-1]["sim_h"]
            W("\n   %.0f-fold range of delta compresses into a %.2fx range of\n"
              % (1.0 / dr if dr < 1 else dr, tr))
            W("   transit time — the signature of a logarithm.  이것이\n")
            out["separatrix"]["delta_ratio"] = 1.0 / dr if dr < 1 else dr
            out["separatrix"]["time_ratio"] = tr
        W("\n   CAVEAT, STATED PLAINLY: lambda above is the QUASI-STATIC eigenvalue,\n")
        W("   with B and NECa held.  The full system restores two amplifiers the\n")
        W("   reduction drops — malabsorption (E down -> SUB up -> B up) and the\n")
        W("   DAMP term (wNEC·NECa) — so the same logarithm is traversed in hours,\n")
        W("   not days.  Measured in the full 39-state system:\n")
        b23 = []
        for key in ("S5_25wk_formula_abx7_indo", "S11_28wk_formula_severe"):
            try:
                pq = scenario_params(key)
            except KeyError:
                continue
            rq = simulate(pq, tmax=45.0, dt=0.005)
            if rq["onset_day"] and rq["stage3_day"]:
                gap = (rq["stage3_day"] - rq["onset_day"]) * 24.0
                b23.append({"scenario": key, "bell2_d": rq["onset_day"],
                            "bell3_d": rq["stage3_day"], "gap_h": gap})
                W("     %-32s Bell II d%.2f -> Bell III d%.2f  = %.1f h\n"
                  % (key, rq["onset_day"], rq["stage3_day"], gap))
            else:
                W("     %-32s did not reach Bell III\n" % key)
        out["separatrix"]["bell2_to_3"] = b23
        W("\n   'Bell II에서 III까지는 늘 하루쯤'의 산술적 이유이며, Bell II에서\n")
        W("   울리는 표지자가 벌어주는 시간이 왜 그렇게 적은지의 이유다.\n\n")

    # -------------------------------------------------------------------------
    # 8.3  scenarios
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("3. TREATMENT / EXPOSURE SCENARIOS (full 39-state system, 56 d)\n")
    W("   치료·노출 시나리오 (전체 시스템 적분)\n")
    W(hr("-") + "\n")
    W("   %-30s %5s %6s %7s %6s %6s %6s %6s %6s\n"
      % ("scenario", "Bell", "onset", "surgic", "E_end", "B_end", "SCFA",
         "bili", "wt_g"))
    scen = {}
    for key, label in SCENARIOS:
        p = scenario_params(key)
        r = simulate(p, tmax=56.0, dt=0.01)
        scen[key] = {k: v for k, v in r.items() if k != "rec"}
        scen[key]["label"] = label
        W("   %-30s %5d %6s %7s %6.3f %6.2f %6.1f %6.2f %6.0f\n"
          % (key, r["peak_bell"],
             ("%.1f" % r["onset_day"]) if r["onset_day"] else "-",
             "yes" if r["surgical"] else "no",
             r["E_end"], r["B_end"], r["SCFA_end"], r["BILI_end"],
             (r["WT_end"] - p["BW"]) * 1000.0))
    out["scenarios"] = scen
    W("\n")

    # -------------------------------------------------------------------------
    # 8.4  feeding is two-signed
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("4. FEEDING ADVANCE IS TWO-SIGNED, AND THE MODEL SAYS THE SUBSTRATE ARM WINS\n")
    W("   수유 증량은 양쪽 부호를 가지지만, 이 파라미터화에서는 기질 쪽이 이긴다\n")
    W(hr("-") + "\n")
    W("   substrate enters Ftroph (+, trophic), dB/dt (+, bacterial substrate)\n")
    W("   and the O2 balance (-, demand = 1 + aDEM·FEED/200).\n")
    W("   NOTE, AGAINST THE HYPOTHESIS: we expected an interior optimum.  The\n")
    W("   integration does NOT produce one for NEC — risk is monotone in the\n")
    W("   advance rate.  The price of caution is therefore paid in a different\n")
    W("   endpoint (days to full feeds, TPN exposure, cholestasis, growth), not\n")
    W("   in NEC.  Both columns are printed so the trade-off is across endpoints\n")
    W("   rather than hidden inside one.\n")
    W("   This cohort is EXCLUSIVELY FORMULA-FED on purpose — a stress test, not\n")
    W("   standard care — so the absolute rates are higher than any real unit.\n\n")
    feed_rows = {}
    for ga, npat in ((25.0, 170), (30.0, 170)):
        W("   GA %.0f wk (n = %d per rate):\n" % (ga, npat))
        W("     %-10s %-9s %-11s %-12s %-10s\n"
          % ("mL/kg/d^2", "NEC %", "surgical %", "days->full", "wt gain g"))
        rows = []
        for rate in (0.0, 10.0, 20.0, 30.0, 40.0):
            arm = {"GA": ga, "fMOM": 0.0, "fDM": 0.0, "fFORM": 1.0,
                   "feed_rate": rate}
            if rate == 0.0:
                arm.update(feed_init=20.0)      # trophic feeds only, no advance
            res = population(npat, arm, seed=771 + int(ga) * 13 + int(rate),
                             tmax=40.0, dt=0.03)
            s = summarise(res)
            dtf = ((160.0 - 20.0) / rate + 2.0) if rate > 0 else float("inf")
            wg = sum((r["WT_end"] - 0.10 * (r["GA"] - 22.0) - 0.35)
                     for r in res) / len(res) * 1000.0
            rows.append({"rate": rate, "nec_pct": s["nec_pct"],
                         "surg_pct": s["surgical_pct"],
                         "days_to_full": (None if rate == 0 else dtf)})
            W("     %-10.0f %-9.1f %-11.1f %-12s %-10.0f\n"
              % (rate, s["nec_pct"], s["surgical_pct"],
                 ("never" if rate == 0 else "%.0f" % dtf), wg))
        feed_rows["GA%d" % int(ga)] = rows
        best = min(rows, key=lambda z: z["nec_pct"])
        span = max(r["nec_pct"] for r in rows) - min(r["nec_pct"] for r in rows)
        W("     -> minimum NEC at %.0f mL/kg/d^2 ; spread across rates = %.1f pp\n\n"
          % (best["rate"], span))
    out["feeding"] = feed_rows
    W("   READ: 같은 방정식이 30주에서는 완만하고 25주에서는 급하다 — 대규모\n")
    W("   증량속도 임상시험(SIFT)의 null 결과와 최미숙아에 대한 임상적 경계심이\n")
    W("   서로 모순이 아니라 하나의 식을 두 재태연령에서 평가한 것이라는 뜻이다.\n")
    W("   그러나 이 모델에서 NEC만 보면 가장 느린 증량이 항상 가장 낮다.\n")
    W("   느린 증량의 대가는 NEC가 아니라 완전영양 도달일·TPN 노출·성장에서\n")
    W("   나타나며, 위 표의 마지막 두 열이 바로 그 대가다.\n\n")

    # -------------------------------------------------------------------------
    # 8.5  antibiotic duration changes sign
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("5. EMPIRICAL ANTIBIOTIC DURATION CHANGES SIGN AT A COMPUTABLE DAY\n")
    W("   경험적 항생제는 계산 가능한 투여기간에서 부호가 바뀐다\n")
    W(hr("-") + "\n")
    W("   MECHANISM: what reshapes the lumen is the LUMINAL concentration.\n")
    W("   Ampicillin arrives there (fl = %.2f) and devastates the anaerobes;\n"
      % P["fl_amp"])
    W("   an intravenous aminoglycoside barely does (fl = %.2f), so the\n"
      % P["fl_gen"])
    W("   Enterobacteriaceae are comparatively spared.  Losing C removes SCFA\n")
    W("   (hence the pH suppression of B) and frees niche.  The benefit is\n")
    W("   immediate and the harm is delayed, so the two integrals must cross.\n")
    W("   Run in HUMAN-MILK-fed infants on purpose: the harm IS the loss of C,\n")
    W("   so it can only appear where C exists.\n\n")
    W("   %-8s %-9s %-11s %-9s %-9s %-9s\n"
      % ("abx d", "NEC %", "surgical %", "B(d21)", "C(d21)", "SCFA(d21)"))
    abx_rows = []
    for days in (0.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0):
        # RUN THIS IN HUMAN-MILK-FED INFANTS ON PURPOSE.  The harm of a long
        # course is the loss of C, so it can only appear where C exists.  In a
        # formula-fed infant there is no commensal community to destroy and the
        # direct suppression of B makes antibiotics look protective — which is
        # itself a falsifiable prediction of this model.
        arm = {"GA": 26.0, "fMOM": 1.0, "fFORM": 0.0, "abx_days": days,
               "feed_rate": 20.0}
        # the deterministic B/C/SCFA columns are shown for the MAJORITY
        # phenotype (no HMO-utilising strain), because a carrier's commensals
        # regrow on HMO after any course and show nothing
        res = population(200, arm, seed=4400 + int(days) * 7, tmax=40.0, dt=0.03)
        s = summarise(res)
        p1 = base_params(); p1.update(arm)
        p1["binfantis"] = 0.05
        p1["_insults"] = list(STD_INSULTS)
        r1 = simulate(p1, tmax=21.0, dt=0.01)
        abx_rows.append({"abx_days": days, "nec_pct": s["nec_pct"],
                         "surg_pct": s["surgical_pct"], "B21": r1["B_end"],
                         "C21": r1["C_end"], "SCFA21": r1["SCFA_end"]})
        W("   %-8.0f %-9.1f %-11.1f %-9.2f %-9.2f %-9.1f\n"
          % (days, s["nec_pct"], s["surgical_pct"], r1["B_end"], r1["C_end"],
             r1["SCFA_end"]))
    out["antibiotics"] = abx_rows
    base = abx_rows[0]["nec_pct"]
    cross = None
    for r in abx_rows[1:]:
        if r["nec_pct"] > base and cross is None:
            cross = r["abx_days"]
    W("\n   sign change (NEC%% exceeds the no-antibiotic arm) at %s days\n"
      % ("%.0f" % cross if cross else "not within 14"))
    out["abx_sign_change_day"] = cross
    W("\n")

    # -------------------------------------------------------------------------
    # 8.6  threshold-movers vs state-movers multiply
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("6. THRESHOLD-MOVERS x STATE-MOVERS MULTIPLY, THEY DO NOT ADD\n")
    W("   문턱을 옮기는 약과 상태를 옮기는 약은 더해지지 않고 곱해진다\n")
    W(hr("-") + "\n")
    W("   thr_boost  enters ONLY Ktlr  -> a PURE threshold-mover (oral 2'-FL,\n")
    W("              anti-TLR4, recombinant PAF-AH): it moves the bifurcation.\n")
    W("   probiotic  enters ONLY B     -> a PURE state-mover: it moves the\n")
    W("              position along the axis, not the axis.\n")
    W("   HUMAN MILK IS BOTH, and that is the point: it raises Ktlr (HMO, sIgA),\n")
    W("   raises Ftroph (EGF/HB-EGF), AND feeds a commensal guild that lowers B.\n")
    W("   So it is not a threshold-mover to be compared against a state-mover —\n")
    W("   it is a combination therapy, and the two single-target arms together\n")
    W("   are the fair comparator.\n\n")
    fac = {}
    W("   %-24s %-9s %-9s %-11s %-9s\n"
      % ("arm", "NEC %", "RR", "surgical %", "onset d"))
    arms = [
        ("formula (reference)", {"fMOM": 0.0, "fFORM": 1.0, "prob_dose": 0.0}),
        ("+ state-mover only", {"fMOM": 0.0, "fFORM": 1.0, "prob_dose": 1.10}),
        ("+ threshold-mover only", {"fMOM": 0.0, "fFORM": 1.0, "thr_boost": 1.6}),
        ("+ both single agents", {"fMOM": 0.0, "fFORM": 1.0, "prob_dose": 1.10,
                                  "thr_boost": 1.6}),
        ("mother's own milk", {"fMOM": 1.0, "fFORM": 0.0, "prob_dose": 0.0}),
    ]
    ref = None
    for label, arm in arms:
        a = dict(arm); a.update({"feed_rate": 20.0, "abx_days": 3.0})
        res = population(260, a, seed=9091, tmax=40.0, dt=0.03)
        s = summarise(res)
        if ref is None:
            ref = s["nec_pct"]
        rr = s["nec_pct"] / ref if ref else float("nan")
        fac[label] = {"nec_pct": s["nec_pct"], "RR": rr,
                      "surg_pct": s["surgical_pct"],
                      "onset": s["median_onset_day"]}
        W("   %-24s %-9.1f %-9.3f %-11.1f %-9s\n"
          % (label, s["nec_pct"], rr, s["surgical_pct"],
             ("%.1f" % s["median_onset_day"]) if s["median_onset_day"] else "-"))
    rr_t = fac["+ threshold-mover only"]["RR"]
    rr_p = fac["+ state-mover only"]["RR"]
    rr_b = fac["+ both single agents"]["RR"]
    rr_m = fac["mother's own milk"]["RR"]
    W("\n   RR(threshold-mover) = %.3f   RR(state-mover) = %.3f\n" % (rr_t, rr_p))
    W("   multiplicative prediction  RR_t x RR_p          = %.3f\n" % (rr_t * rr_p))
    W("   additive prediction        1-(1-RR_t)-(1-RR_p)  = %.3f\n"
      % (1 - (1 - rr_t) - (1 - rr_p)))
    W("   OBSERVED combination                       RR   = %.3f\n" % rr_b)
    W("   -> multiplicative error %+.3f, additive error %+.3f\n"
      % (rr_t * rr_p - rr_b, (1 - (1 - rr_t) - (1 - rr_p)) - rr_b))
    W("\n   and for comparison, the intervention that is BOTH at once:\n")
    W("   RR(mother's own milk) = %.3f\n" % rr_m)
    out["factorial"] = {"arms": fac, "rr_threshold": rr_t, "rr_state": rr_p,
                        "rr_both": rr_b, "rr_milk": rr_m,
                        "mult_pred": rr_t * rr_p,
                        "add_pred": 1 - (1 - rr_t) - (1 - rr_p)}
    nnt_m = (100.0 / (ref - fac["mother's own milk"]["nec_pct"])
             if ref > fac["mother's own milk"]["nec_pct"] else float("inf"))
    W("   NNT(mother's own milk, this population) = %.0f\n" % nnt_m)
    out["factorial"]["NNT_milk"] = nnt_m
    W("\n")

    # -------------------------------------------------------------------------
    # 8.7  GA gradient and calibration check
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("7. CALIBRATION AGAINST PUBLISHED AGGREGATE ENDPOINTS\n")
    W("   문헌 집계 지표에 대한 보정 확인\n")
    W(hr("-") + "\n")
    W("   %-14s %-9s %-11s %-12s %-14s\n"
      % ("GA band", "NEC %", "surgical %", "onset (d)", "target NEC %"))
    targets = {"24-25": "12-15", "26-27": "8-11", "28-29": "4-6", "30-32": "1-3"}
    ga_rows = []
    for lo, hi, key in ((24.0, 26.0, "24-25"), (26.0, 28.0, "26-27"),
                        (28.0, 30.0, "28-29"), (30.0, 32.5, "30-32")):
        acc = []
        for i in range(170):
            p = base_params()
            u = lcg(31337 + i * 977 + int(lo * 100))
            ga = lo + (hi - lo) * u()
            p["GA"] = ga
            p["BW"] = max(0.42, 0.10 * (ga - 22.0) + 0.35 + 0.09 * norm(u))
            for k, cv in (("Ktlr0", 0.22), ("kA", 0.20), ("muB", 0.20),
                          ("seedB_env", 0.45), ("kPAF", 0.25),
                          ("seedC_hm", 1.00), ("seedC_env", 1.00)):
                p[k] *= math.exp(cv * norm(u))
            p["binfantis"] = 1.0 if u() < 0.35 else 0.05
            p["dysbiosis"] = math.exp(0.35 * norm(u))
            # real-world mix: ~60 % human milk, ~55 % get early antibiotics
            if u() < 0.60:
                p.update(fMOM=1.0, fDM=0.0, fFORM=0.0)
            else:
                p.update(fMOM=0.0, fDM=0.3, fFORM=0.7)
            p["abx_days"] = 3.0 if u() < 0.55 else 0.0
            p["_insults"] = draw_insults(u, ga, tmax=42.0)
            acc.append(simulate(p, tmax=42.0, dt=0.03))
        s = summarise(acc)
        ga_rows.append({"band": key, **s, "target": targets[key]})
        W("   %-14s %-9.1f %-11.1f %-12s %-14s\n"
          % (key, s["nec_pct"], s["surgical_pct"],
             ("%.1f" % s["median_onset_day"]) if s["median_onset_day"] else "-",
             targets[key]))
    out["calibration"] = ga_rows
    W("\n   targets: NEC incidence by GA band (Neonatal Research Network /\n")
    W("   Vermont Oxford style registries); onset 8-21 d, later at lower GA;\n")
    W("   surgical fraction ~ 1/3 - 1/2 of definite NEC.\n\n")

    # -------------------------------------------------------------------------
    # 8.8  biomarker lead time
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("8. BIOMARKER LEAD TIME RELATIVE TO BELL II\n")
    W("   Bell II 기준 대비 표지자 선행 시간\n")
    W(hr("-") + "\n")
    p = scenario_params("S5_25wk_formula_abx7_indo")
    r = simulate(p, tmax=45.0, dt=0.005, record=True, rec_every=20)
    on = r["onset_day"]
    W("   reference case: %s  (Bell II at day %.2f)\n\n"
      % ("S5_25wk_formula_abx7_indo", on if on else float("nan")))
    W("   %-12s %-12s %-12s\n" % ("marker", "day crossed", "lead (h)"))
    lead = {}
    for nm, lbl in (("CRP12", "CRP > 12"), ("PLT100", "PLT < 100"),
                    ("LAC4", "lactate > 4"), ("PNEU30", "pneumatosis")):
        d = r["markers"].get(nm)
        if d is not None and on is not None:
            lead[lbl] = (on - d) * 24.0
            W("   %-12s %-12.2f %-12.1f\n" % (lbl, d, (on - d) * 24.0))
        else:
            W("   %-12s %-12s %-12s\n" % (lbl, "-", "-"))
    out["biomarker_lead_h"] = lead
    # crossing of the mechanistic quantity itself
    rec = r["rec"]
    Jt = None
    for i, t in enumerate(rec["t"]):
        if rec["Jtr"][i] > 3.0 and Jt is None:
            Jt = t
    if Jt is not None and on is not None:
        W("\n   mechanistic quantity Jtr > 3.0 crossed at day %.2f "
          "(lead %.1f h)\n" % (Jt, (on - Jt) * 24.0))
        out["Jtr_lead_h"] = (on - Jt) * 24.0
    W("\n   READ: 임상 표지자는 이미 스위치가 넘어간 뒤에 움직인다. 앞서는 것은\n")
    W("   상태변수(E, Jtr)이며, 그것이 왜 예방이 조기진단을 이기는지의 이유다.\n\n")

    # -------------------------------------------------------------------------
    # 8.9  organ-system consequences
    # -------------------------------------------------------------------------
    W(hr("-") + "\n")
    W("9. EXTRA-INTESTINAL CONSEQUENCES (liver / brain / kidney / growth)\n")
    W("   장 외 결과 (간·뇌·신장·성장)\n")
    W(hr("-") + "\n")
    W("   %-30s %-9s %-9s %-9s %-9s\n"
      % ("scenario", "bili", "NDI ix", "tubular", "wt gain g"))
    organ = {}
    for key in ("S1_28wk_MOM", "S2_28wk_formula", "S5_25wk_formula_abx7_indo",
                "S9_28wk_NPO_TPN", "S10_28wk_formula_abx10_mtz"):
        p = scenario_params(key)
        r = simulate(p, tmax=56.0, dt=0.01)
        organ[key] = {"bili": r["BILI_end"], "NIN": r["NIN_end"],
                      "KIN": r["KIN_end"],
                      "wt_gain_g": (r["WT_end"] - p["BW"]) * 1000.0}
        W("   %-30s %-9.2f %-9.3f %-9.3f %-9.0f\n"
          % (key, r["BILI_end"], r["NIN_end"], r["KIN_end"],
             (r["WT_end"] - p["BW"]) * 1000.0))
    out["organ"] = organ
    W("\n")

    W(hr() + "\n")
    W("END OF EXECUTED OUTPUT\n")
    W(hr() + "\n")

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "nec_scenario_results.json"), "w") as f:
        json.dump({k: out[k] for k in
                   ("bifurcation", "separatrix", "scenarios", "feeding",
                    "antibiotics", "abx_sign_change_day", "biomarker_lead_h",
                    "organ")
                   if k in out}, f, indent=1, default=str)
    with open(os.path.join(here, "nec_population_results.json"), "w") as f:
        json.dump({k: out[k] for k in ("factorial", "calibration") if k in out},
                  f, indent=1, default=str)


if __name__ == "__main__":
    main()
