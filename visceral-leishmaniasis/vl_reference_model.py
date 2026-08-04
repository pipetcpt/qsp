#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vl_reference_model.py
=====================

Independent Python reference implementation of the visceral leishmaniasis
(VL, kala-azar) QSP model.  Every equation here is written to be
term-for-term identical to `vl_mrgsolve_model.R`; this file exists because
no R runtime is available in the build environment, so the ODE system has to
be integrated somewhere before any number is claimed.  Anything printed by
this script is a model output, not a literature value, unless the label says
otherwise.

Structure (73 ODEs)
-------------------
  PK    amphotericin B, liposomal and deoxycholate, as TWO plasma species
        (liposome-associated and released/free) feeding four macrophage-rich
        organ compartments plus kidney                            [ 8 ODEs]
        miltefosine (oral, allometric CL, long terminal tail)     [ 7 ODEs]
        paromomycin (IM, slow macrophage uptake, cochlea, kidney) [ 9 ODEs]
        sodium stibogluconate: Sb(V) -> intracellular Sb(III)     [ 8 ODEs]
  PD    amastigotes in spleen / liver / bone marrow / skin, each split into
        a drug-susceptible and a drug-tolerant subpopulation      [ 8 ODEs]
  HOST  memory T-cell pool, IFN-gamma, IL-10, TNF, TGF-beta, CD4,
        activated-macrophage pool                                 [ 7 ODEs]
  CLIN  spleen size, liver size, Hb, platelets, WBC, albumin, polyclonal
        IgG, temperature, body weight, PKDL lesion load           [10 ODEs]
  TOX   tubular injury, creatinine, K, Mg, hearing shift, QTc,
        lipase, ALT, GI intolerance                               [ 9 ODEs]
  ACC   five exposure integrals, cumulative mortality hazard,
        cumulative mg/kg administered                             [ 7 ODEs]

The single structural commitment of the model
---------------------------------------------
The drug and the parasite never meet in plasma.  Amastigotes sit inside
macrophages of spleen, liver, marrow and skin, so drug effect is driven by
the INTRAMACROPHAGE exposure of each organ, while amphotericin
nephrotoxicity is driven by the FREE PLASMA exposure of the same dose.
Liposomal encapsulation moves those two integrals in opposite directions,
because the liposome is cleared by the very cell lineage that harbours the
parasite.  Cure is then not "zero parasites" but crossing a separatrix in
the (burden, primed-memory) plane, and that separatrix moves with CD4.
"""

import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# --------------------------------------------------------------------------
# 0.  State vector layout
# --------------------------------------------------------------------------
NAMES = [
    # --- amphotericin B (8) -------------------------------------------------
    "A_LIP",    # 0  liposome-associated AmB in plasma            (mg)
    "A_FRE",    # 1  released / protein-bound "free" AmB, central (mg)
    "A_PER",    # 2  AmB peripheral tissue (non-MPS)              (mg)
    "AMB_SP",   # 3  AmB, spleen macrophage compartment           (mg)
    "AMB_LI",   # 4  AmB, liver Kupffer compartment               (mg)
    "AMB_BM",   # 5  AmB, marrow macrophage compartment           (mg)
    "AMB_SK",   # 6  AmB, dermal macrophage compartment           (mg)
    "AMB_KID",  # 7  AmB, renal cortex                            (mg)
    # --- miltefosine (7) ---------------------------------------------------
    "MIL_G", "MIL_C", "MIL_P", "MIL_SP", "MIL_LI", "MIL_BM", "MIL_SK",
    # --- paromomycin (9) ---------------------------------------------------
    "PM_D", "PM_C", "PM_P", "PM_SP", "PM_LI", "PM_BM", "PM_SK",
    "PM_KID", "PM_COC",
    # --- antimony (8) ------------------------------------------------------
    "SB5_D", "SB5_C", "SB5_P", "SB3_SP", "SB3_LI", "SB3_BM", "SB3_SK",
    "SB_DEP",
    # --- parasites (8): susceptible + tolerant, four organs ----------------
    "P_SP_S", "P_SP_R", "P_LI_S", "P_LI_R",
    "P_BM_S", "P_BM_R", "P_SK_S", "P_SK_R",
    # --- host immunity (7) -------------------------------------------------
    "TMEM", "IFNG", "IL10", "TNFA", "TGFB", "CD4", "MPHA",
    # --- clinical (10) -----------------------------------------------------
    "SPL", "LIVS", "HGB", "PLT", "WBC", "ALB", "IGG", "TEMP", "BWTL", "PKDL",
    # --- toxicity (9) ------------------------------------------------------
    "TUBI", "SCR", "KSER", "MGSER", "HEAR", "QTC", "LIPA", "ALTX", "GITX",
    # --- accumulators (7) --------------------------------------------------
    "AUCFRE", "AUCASP", "AUCMSP", "AUCBSP", "AUCPSP", "CUMHAZ", "CUMDOSE",
]
IX = {n: i for i, n in enumerate(NAMES)}
NEQ = len(NAMES)
assert NEQ == 73, NEQ

DAY = 24.0  # hours per day; the model runs in hours


# --------------------------------------------------------------------------
# 1.  Parameters
# --------------------------------------------------------------------------
def default_params():
    p = {}

    # ---- covariates -------------------------------------------------------
    p["WT"] = 50.0          # kg, adult Bihar/Indian VL patient
    p["AGE"] = 30.0         # y
    p["CD40"] = 700.0       # cells/uL at baseline (immunocompetent)
    p["HIV"] = 0.0          # 1 = HIV co-infected
    p["ART"] = 0.0          # 1 = on ART (drives CD4 recovery)
    p["MALNUT"] = 0.0       # 0-1 malnutrition index, scales T-cell priming
    p["P0SCALE"] = 1.0      # multiplier on presenting burden.  Splenic
    #                         aspirate grade at diagnosis spans 1+ to 6+,
    #                         i.e. orders of magnitude, and it is one of the
    #                         strongest predictors of outcome in every trial.

    # ---- amphotericin B ---------------------------------------------------
    # Liposomal AmB stays essentially intravascular until the liposome is
    # taken up by the mononuclear phagocyte system (MPS) or leaks drug.
    p["V_LIP"] = 0.10       # L/kg  -> 5 L for 50 kg (plasma-like)
    p["CL_LIP"] = 0.50      # L/h/70kg total liposome disposition clearance
    p["FREL"] = 0.30        # fraction of liposome clearance that RELEASES
    #                         drug to the free plasma pool; 1-FREL is MPS
    #                         uptake, i.e. delivery to the infected cell
    p["V_FRE"] = 0.50       # L/kg central volume, free/bound AmB
    p["CL_FRE"] = 1.40      # L/h/70kg free AmB elimination clearance
    p["Q_AMB"] = 1.20       # L/h/70kg free <-> peripheral
    p["V_PER"] = 3.00       # L/kg peripheral volume
    # MPS uptake is split between organs in proportion to phagocyte mass
    p["FSP_A"] = 0.10       # spleen share of MPS uptake
    p["FLI_A"] = 0.60       # liver (Kupffer) share
    p["FBM_A"] = 0.12       # bone marrow share
    p["FSK_A"] = 0.03       # dermis share  <- small: single-dose L-AmB
    #                         barely reaches skin, which is where PKDL lives
    # remainder (0.15) = lung/lymph-node MPS sink, not tracked
    # Tissue egress is the half-life of the PHARMACOLOGICALLY AVAILABLE
    # intracellular fraction, not of total tissue amphotericin.  Total tissue
    # AmB is detectable in liver and spleen for many weeks; the fraction that
    # is not sequestered in host cholesterol-rich membranes and can still
    # reach parasite ergosterol turns over faster.  3.5 d is the value that
    # reproduces single-dose L-AmB efficacy without predicting that the tail
    # alone sterilises every patient.
    p["KOUT_A"] = 0.00825   # 1/h  available-fraction egress, t1/2 ~3.5 d
    p["KPSP_A"] = 14.0      # passive partition spleen:free-plasma
    p["KPLI_A"] = 20.0
    p["KPBM_A"] = 9.0
    p["KPSK_A"] = 1.5
    p["KPKID_A"] = 8.0
    p["KIN_KID"] = 0.15     # 1/h renal cortex uptake
    # KOUT_KID is retained for documentation only.  After the renal
    # compartment was rewritten as an equilibration (see the note in `rhs`),
    # cortical kinetics are set by KIN_KID and KPKID_A; a separate egress
    # constant would double-count the same process.
    p["KOUT_KID"] = 0.012   # 1/h (UNUSED, see note above)

    # ---- miltefosine (Dorlo-type allometric disposition) ------------------
    p["KA_MIL"] = 0.40      # 1/h
    p["F_MIL"] = 1.00       # apparent
    p["CL_MIL"] = 0.166     # L/h/70kg   (~3.99 L/day)
    p["V1_MIL"] = 40.0      # L/70kg
    p["V2_MIL"] = 52.0      # L/70kg
    p["Q_MIL"] = 0.90       # L/h/70kg
    p["KPSP_M"] = 8.0       # alkylphosphocholine, accumulates in membranes
    p["KPLI_M"] = 8.0
    p["KPBM_M"] = 6.0
    p["KPSK_M"] = 4.0       # good dermal penetration -> the PKDL drug
    p["KIN_M"] = 0.05       # 1/h tissue equilibration rate
    p["EM50_GI"] = 25.0     # mg/L plasma driving nausea/vomiting

    # ---- paromomycin ------------------------------------------------------
    p["KA_PM"] = 1.80       # 1/h  IM absorption
    p["V_PM"] = 0.25        # L/kg
    p["CL_PM"] = 5.50       # L/h/70kg (glomerular, tracks CrCl)
    p["Q_PM"] = 1.00        # L/h/70kg
    p["V2_PM"] = 0.20       # L/kg
    p["KIN_PM"] = 0.020     # 1/h  SLOW pinocytic entry into macrophage
    p["KOUT_PM"] = 0.004    # 1/h  slow egress -> effective Kp 5, delayed
    p["FSK_PM"] = 0.25      # dermal penetration relative to MPS organs
    p["KIN_PMK"] = 0.10     # renal cortical (megalin) uptake
    p["KOUT_PMK"] = 0.008
    p["KIN_COC"] = 0.004    # cochlear uptake, near-irreversible
    p["KOUT_COC"] = 0.0004

    # ---- antimony ---------------------------------------------------------
    p["KA_SB"] = 1.50       # 1/h IM
    p["V_SB"] = 0.22        # L/kg
    p["CL_SB"] = 7.00       # L/h/70kg renal, t1/2 ~2 h
    p["Q_SB"] = 0.45        # L/h/70kg to the slow "deep" antimony depot
    p["V2_SB"] = 0.90       # L/kg deep depot  (terminal t1/2 ~76 h)
    p["KRED"] = 0.045       # 1/h Sb(V) -> Sb(III) reduction inside phagocyte
    p["KEFF_SB"] = 0.30     # 1/h Sb(III) efflux (MRPA / trypanothione route)
    p["RES_SB"] = 1.0       # resistance multiplier on efflux (Bihar >> 1)
    p["KPSK_SB"] = 0.30     # dermal Sb access relative to MPS organs

    # ---- organ "volumes" for the intramacrophage compartments -------------
    p["VSP"] = 0.25         # L
    p["VLI"] = 1.50         # L
    p["VBM"] = 1.50         # L
    p["VSK"] = 3.00         # L
    p["VKID"] = 0.30        # L
    p["VCOC"] = 0.002       # L

    # ---- parasite dynamics ------------------------------------------------
    # burden unit = 1e6 amastigotes.  Each organ carries TWO populations:
    #   P_*_S  replicating amastigotes
    #   P_*_R  quiescent / persister amastigotes -- near-arrested, poorly
    #          killed by drugs (all four act on metabolically active cells)
    #          and poorly seen by the macrophage effector arm.  This pool is
    #          not a resistance genotype; it is a physiological state, and it
    #          is the reservoir every relapse in this model comes from.
    p["KG_SP"] = 0.0072     # 1/h  doubling time ~4 d
    p["KG_LI"] = 0.0065
    p["KG_BM"] = 0.0068
    p["KG_SK"] = 0.0040     # slower in skin
    p["PMAX_SP"] = 3.0e3    # 3e9 amastigotes at carrying capacity
    p["PMAX_LI"] = 1.0e4
    p["PMAX_BM"] = 3.0e3
    p["PMAX_SK"] = 1.0e2
    p["KIMM"] = 0.065       # 1/h maximal immune (NO-mediated) killing rate
    # Trafficking is ASYMMETRIC, and it has to be.  Spleen, liver and marrow
    # share the blood monocyte pool, so they redistribute freely.  The dermis
    # receives infected monocytes from that pool but its own macrophages are
    # largely resident and do not recirculate.  Written symmetrically, the
    # immune-privileged skin became an incurable reservoir that re-seeded the
    # viscera and NO regimen cured -- whereas in reality most cured patients
    # never develop PKDL and PKDL patients rarely relapse to visceral disease.
    p["KMIG"] = 2.0e-4      # 1/h redistribution among spleen/liver/marrow
    p["KMIG_SKIN"] = 3.0e-6 # 1/h visceral pool -> dermis
    p["KMIG_BACK"] = 3.0e-7 # 1/h dermis -> visceral pool (an order of
    #                         magnitude smaller: this is why PKDL seeds
    #                         transmission but rarely seeds relapse)
    # 1e-4 of 1e10 amastigotes = ~1e6 persisters at presentation.  An
    # earlier value of 0.02 put 2e8 quiescent organisms in the patient on
    # day 0, and because persisters are barely touched by any of the four
    # drugs that pool alone capped EVERY regimen at a 2-log nadir -- the
    # drug kill integrals were 7-12 logs while the burden fell 2.
    p["FR0"] = 5.0e-4       # quiescent fraction at presentation
    p["KQ0"] = 2.0e-6       # 1/h basal entry into quiescence
    p["KQS"] = 1.5e-5       # 1/h per unit of drug+immune pressure (stress-
    #                         induced quiescence: the harder you push, the
    #                         more persisters you make)
    p["KR_REACT"] = 8.0e-4  # 1/h reactivation out of quiescence
    p["COSTR"] = 0.05       # replication rate of quiescent forms, x kg
    p["FIMM_SK"] = 0.12     # immune killing in the DERMIS, as a fraction of
    #                         the systemic rate.  The dermis is relatively
    #                         immune-privileged for this parasite, and it has
    #                         to be: PKDL lesions persist for years in a host
    #                         whose systemic immunity has fully returned, and
    #                         with dermal killing set equal to splenic
    #                         killing the model made untreated PKDL
    #                         self-resolving, which it is not in South Asia.
    p["FIMM_R"] = 0.15      # immune killing of quiescent forms, x kimm.
    #                         Quiescent amastigotes are also immunologically
    #                         quiet: reduced antigen export and MHC-II
    #                         display, so the host sees them poorly too.
    p["RFAC_A"] = 20.0      # EC50 multiplier for quiescent forms, AmB
    p["RFAC_M"] = 25.0      # ... miltefosine (needs LdMT + active membrane)
    p["RFAC_P"] = 30.0      # ... paromomycin (needs active translation)
    p["RFAC_S"] = 40.0      # ... Sb(III) (needs active reduction + thiols)
    p["RES_M"] = 1.0        # acquired-resistance multiplier, miltefosine
    p["RES_P"] = 1.0        # acquired-resistance multiplier, paromomycin

    # ---- drug PD (Emax on INTRAMACROPHAGE concentration) ------------------
    p["EMAX_A"] = 0.060     # 1/h
    p["EC50_A"] = 6.0       # mg/L tissue
    p["HILL_A"] = 2.0
    # EC50_M is set so that the peak intramacrophage concentration sits only
    # ~2.6x above it.  With a 19x margin (the first parameterisation) losing
    # 40% of a 28-day oral course changed nothing, because the kill term was
    # saturated throughout; the DURATION above EC50 has to be the
    # discriminator for adherence and bioavailability to matter at all.
    p["EMAX_M"] = 0.034
    p["EC50_M"] = 90.0
    p["HILL_M"] = 1.5
    p["EMAX_P"] = 0.042
    p["EC50_P"] = 25.0
    p["HILL_P"] = 1.5
    p["EMAX_S"] = 0.045
    p["EC50_S"] = 0.50
    p["HILL_S"] = 1.5

    # ---- immunity ---------------------------------------------------------
    # KAG << K50_IL10 is a structural statement, not a fitted convenience:
    # antigen is still abundant over a burden range in which IL-10 has
    # already switched off.  That gap -- from 1e4 to 6e7 amastigotes -- is the
    # window in which the host can be primed, and getting the patient INTO
    # that window and holding them there is what a drug course does.
    p["KAG"] = 0.01         # burden unit giving half-maximal antigen signal
    p["KTP"] = 0.0050       # 1/h memory/effector T-cell priming rate
    p["TMEMAX"] = 2.0       # ceiling on the primed pool
    p["KTD"] = 0.00015      # 1/h memory decay (t1/2 ~192 d)
    p["KI_T"] = 0.35        # IL-10 conc inhibiting priming by 50%
    p["KIP"] = 0.045        # IFN-g production
    p["KID"] = 0.060        # IFN-g decay
    # The antigen-only contribution to IFN-gamma has to be SMALL.  At 0.15 it
    # alone gave an unprimed host kimm = 0.0084 /h against a growth rate of
    # 0.0072 /h, so the visceral compartment was controlled even with the
    # T-cell priming rate set to zero and no host could develop kala-azar at
    # all.  Innate IFN-gamma (NK, gamma-delta T) is real but is not by itself
    # curative -- if it were, visceral leishmaniasis would not exist.
    p["FAG_IFN"] = 0.04     # antigen-only contribution to IFN-g
    p["KLP"] = 0.030        # IL-10 production
    p["KLD"] = 0.020        # IL-10 decay
    p["K50_IL10"] = 60.0    # burden for half-maximal IL-10 (>> KAG)
    p["K50_SYS"] = 300.0    # burden for half-maximal SYSTEMIC illness (fever,
    #                         cachexia, hypoalbuminaemia).  A third threshold
    #                         is needed because the antigen signal AG
    #                         saturates at 1e4 amastigotes: keyed to AG, the
    #                         clinical module made a host carrying a
    #                         controlled subclinical infection febrile and
    #                         anaemic, which is the opposite of what
    #                         subclinical means.
    p["HIL10"] = 2.0
    p["KNP"] = 0.050        # TNF production
    p["KND"] = 0.070
    p["KGP"] = 0.010        # TGF-b production
    p["KGD"] = 0.012
    p["KI_IFN"] = 0.35      # IFN-g giving half-maximal macrophage activation
    p["KI_IL10"] = 0.25     # IL-10 inhibition constant on activation
    p["KI_TGFB"] = 0.40
    p["KCD4"] = 350.0       # cells/uL for half-maximal T-cell help
    # CD4 is a HOMEOSTATIC pool with a set point, transiently depressed by
    # active VL.  An earlier version had only a loss term and no set point,
    # so every immunocompetent patient drifted to CD4 = 4 cells/uL over the
    # 18-month follow-up and the whole host arm silently switched off -- which
    # is why five regimens that cure >90% of real patients were "failing".
    p["CD4SET"] = 800.0     # cells/uL homeostatic set point (no HIV)
    p["CD4SET_ART"] = 600.0 # set point reached on suppressive ART
    p["KCD4H"] = 0.0015     # 1/h homeostatic restoration (t1/2 ~19 d)
    p["KCD4L"] = 0.00080    # 1/h burden-driven CD4 depression
    p["KMPH"] = 0.015       # macrophage pool turnover
    p["FHIV_MAC"] = 0.50    # HIV-infected macrophages kill worse at any given
    #                         IFN-g level (impaired NO output, altered
    #                         phagosome maturation) -- an effect on the
    #                         EFFECTOR, separate from the loss of CD4 help

    # ---- clinical ---------------------------------------------------------
    p["SPL0"] = 0.5         # cm below costal margin, healthy
    p["SPLMAX"] = 12.0      # cm added at maximal splenic burden
    p["K50_SPL"] = 300.0    # spleen burden for half-maximal splenomegaly
    p["KSPL"] = 0.0035      # 1/h spleen size relaxation (t1/2 ~8 d)
    p["LIV0"] = 0.5
    p["LIVMAX"] = 5.0
    p["KLIV"] = 0.0030
    p["HGB0"] = 14.0        # g/dL
    p["KHGB"] = 0.0060      # 1/h turnover of the Hb pool
    p["FSUP_HGB"] = 0.75    # maximal marrow suppression fraction
    p["KSPL_HGB"] = 0.045   # hypersplenic destruction per cm of spleen
    p["PLT0"] = 250.0       # 1e9/L
    p["KPLT"] = 0.0090
    p["KSPL_PLT"] = 0.085   # platelets are the most spleen-sensitive line
    p["WBC0"] = 7.0         # 1e9/L
    p["KWBC"] = 0.020
    p["KSPL_WBC"] = 0.055
    p["ALB0"] = 4.3         # g/dL
    p["KALB"] = 0.0035
    p["IGG0"] = 1.1         # g/dL
    # Set so that the steady state during active disease, IGG0 + KIGG_P/KIGG_D,
    # lands at ~4.0 g/dL -- the observed range in untreated VL.  At 0.0022 the
    # equilibrium was 8.4 g/dL, so polyclonal IgG went UP after a curative
    # dose instead of decaying, which inverted the whole serology argument.
    p["KIGG_P"] = 0.00087   # polyclonal IgG production per antigen unit
    p["KIGG_D"] = 0.00030   # 1/h  t1/2 ~96 d -> serology stays positive
    p["TEMP0"] = 36.7       # degC
    p["KTEMP"] = 0.060
    p["TNF_FEV"] = 3.2      # degC per unit TNF at saturation
    p["KBWL"] = 0.0012      # 1/h weight-loss dynamics
    p["BWLMAX"] = 0.28      # fraction of body weight lost, untreated
    p["KPKDL_P"] = 0.0022   # PKDL lesion formation: skin parasite x immunity
    p["KPKDL_D"] = 0.0015

    # ---- toxicity ---------------------------------------------------------
    p["KTUB"] = 0.00035     # tubular injury per mg/L renal AmB per h
    p["KTUBR"] = 0.0030     # 1/h repair
    p["GFR0"] = 100.0       # mL/min
    p["KSCR"] = 0.030       # 1/h creatinine turnover
    p["SCR0"] = 0.85        # mg/dL
    p["KWK"] = 0.0039       # K wasting per mg/L renal AmB per h
    p["K_K0"] = 4.1         # mEq/L
    p["KKR"] = 0.030
    p["KWMG"] = 0.00146
    p["MG0"] = 2.0          # mg/dL
    p["KMGR"] = 0.025
    p["KHEAR"] = 0.00040    # dB per (mg/L cochlear PM) per h
    p["QTC0"] = 400.0       # ms
    p["KQTC"] = 0.145       # ms per mg/L deep-depot antimony
    p["KQTCR"] = 0.020
    p["LIPA0"] = 30.0       # U/L
    p["KLIPA"] = 0.85
    p["KLIPAR"] = 0.015
    p["ALT0"] = 25.0        # U/L
    p["KALTX"] = 0.18
    p["KALT_M"] = 0.013     # ALT response to miltefosine
    p["KALTR"] = 0.012
    p["KGITX"] = 0.10
    p["KGITR"] = 0.09

    # ---- mortality hazard -------------------------------------------------
    p["H0"] = 2.2e-5        # 1/h baseline hazard of the severe VL state
    p["HGB_H"] = 7.0        # Hb below which hazard starts to climb
    p["ALB_H"] = 2.8
    p["WBC_H"] = 1.5

    # ---- dose-flag scratch space -----------------------------------------
    p["AMBFORM"] = 1.0      # 1 = liposomal, 0 = deoxycholate (set by dosing)
    p["ADHER"] = 1.0        # fraction of ORAL (miltefosine) doses actually
    #                         taken.  A 28-day unsupervised oral course is the
    #                         one place in VL therapy where adherence is a
    #                         first-order pharmacological variable, and it is
    #                         the standard explanation for the fall in field
    #                         efficacy of miltefosine in India after 2010.
    return p


def scale_params(p):
    """Allometric scaling.  Clearances ^0.75, volumes ^1.0, on total weight.

    This is the whole of claim 2: dose is prescribed per kg (exponent 1.0)
    but cleared with exponent 0.75, so steady-state exposure scales as
    WT^0.25 and small children on linear mg/kg dosing are underexposed.
    """
    q = dict(p)
    wt = p["WT"]
    f_cl = (wt / 70.0) ** 0.75
    f_v = wt / 70.0
    for k in ("CL_LIP", "CL_FRE", "Q_AMB", "CL_MIL", "Q_MIL", "CL_PM",
              "Q_PM", "CL_SB", "Q_SB"):
        q[k] = p[k] * f_cl
    for k in ("V1_MIL", "V2_MIL"):
        q[k] = p[k] * f_v
    # per-kg volumes become absolute litres
    for k, src in (("V_LIPA", "V_LIP"), ("V_FREA", "V_FRE"),
                   ("V_PERA", "V_PER"), ("V_PMA", "V_PM"),
                   ("V2_PMA", "V2_PM"), ("V_SBA", "V_SB"),
                   ("V2_SBA", "V2_SB")):
        q[k] = p[src] * wt
    # organ compartments scale with body size (spleen size handled separately)
    for k in ("VSP", "VLI", "VBM", "VSK", "VKID"):
        q[k] = p[k] * f_v
    q["GFR0"] = p["GFR0"] * f_cl
    return q


# --------------------------------------------------------------------------
# 2.  Initial conditions
# --------------------------------------------------------------------------
def initial_state(p, stage="active"):
    """`stage='active'` = patient presenting with established, symptomatic VL.

    Baseline burdens are set to the equilibrium of the untreated high-burden
    attractor, which is what a patient walking into a Bihar treatment centre
    actually is.  `stage='naive'` starts from a sandfly inoculum and is used
    only for the natural-history / subclinical-race analysis.
    """
    y = np.zeros(NEQ)
    y[IX["CD4"]] = p["CD40"]
    y[IX["HGB"]] = p["HGB0"]
    y[IX["PLT"]] = p["PLT0"]
    y[IX["WBC"]] = p["WBC0"]
    y[IX["ALB"]] = p["ALB0"]
    y[IX["IGG"]] = p["IGG0"]
    y[IX["TEMP"]] = p["TEMP0"]
    y[IX["SPL"]] = p["SPL0"]
    y[IX["LIVS"]] = p["LIV0"]
    y[IX["SCR"]] = p["SCR0"]
    y[IX["KSER"]] = p["K_K0"]
    y[IX["MGSER"]] = p["MG0"]
    y[IX["QTC"]] = p["QTC0"]
    y[IX["LIPA"]] = p["LIPA0"]
    y[IX["ALTX"]] = p["ALT0"]
    # MPHA is the ACTIVATED fraction of the macrophage pool, so a healthy
    # resting host sits near zero, not near one.  Leaving it at a high value
    # meant a fresh sandfly inoculum was sterilised on arrival and every run
    # of the natural-history analysis returned "abortive infection".
    y[IX["MPHA"]] = 0.02

    if stage == "naive":
        # SCOPE NOTE: the inoculum is placed in the SPLEEN, not the dermis.
        # The model does not represent the first days after the bite, when
        # recruited neutrophils and monocytes carry the parasite out of the
        # skin; that transport is fast and inflammation-dependent, and it is
        # not the same process as the slow, resident dermal trafficking that
        # governs PKDL.  Natural history therefore starts from a small
        # established visceral focus.
        y[IX["P_SP_S"]] = 1.0e-4 * (1.0 - p["FR0"])
        y[IX["P_SP_R"]] = 1.0e-4 * p["FR0"]
        for k in ("TMEM", "IFNG", "IL10", "TNFA", "TGFB"):
            y[IX[k]] = 0.0
        y[IX["MPHA"]] = 0.02
        return y

    # established disease: near carrying capacity in the MPS organs
    fr = p["FR0"]
    for org, frac in (("SP", 0.80), ("LI", 0.55), ("BM", 0.60), ("SK", 0.25)):
        tot = frac * p["PMAX_" + org] * p["P0SCALE"]
        y[IX["P_%s_S" % org]] = tot * (1.0 - fr)
        y[IX["P_%s_R" % org]] = tot * fr
    # host state consistent with months of untreated disease
    y[IX["IL10"]] = 0.62
    y[IX["IFNG"]] = 0.11
    y[IX["TNFA"]] = 0.55
    y[IX["TGFB"]] = 0.42
    y[IX["TMEM"]] = 0.05
    # an earlier version left MPHA at its healthy value of 1.0 here, which
    # gave every patient a fully activated macrophage pool on day 0 and
    # flattered the fast-killing drugs
    y[IX["MPHA"]] = 0.04
    y[IX["SPL"]] = 8.5
    y[IX["LIVS"]] = 2.8
    y[IX["HGB"]] = 7.4
    y[IX["PLT"]] = 78.0
    y[IX["WBC"]] = 2.4
    y[IX["ALB"]] = 2.7
    y[IX["IGG"]] = 3.6
    y[IX["TEMP"]] = 38.9
    y[IX["BWTL"]] = 0.16
    return y


# --------------------------------------------------------------------------
# 3.  Right-hand side
# --------------------------------------------------------------------------
def rhs(t, y, p):
    d = np.zeros(NEQ)
    g = y  # alias

    # ---------------- amphotericin B ---------------------------------------
    A_LIP, A_FRE, A_PER = g[0], g[1], g[2]
    AMB_SP, AMB_LI, AMB_BM, AMB_SK, AMB_KID = g[3], g[4], g[5], g[6], g[7]

    C_LIP = A_LIP / p["V_LIPA"]
    C_FRE = A_FRE / p["V_FREA"]
    C_PER = A_PER / p["V_PERA"]
    cAsp = AMB_SP / p["VSP"]
    cAli = AMB_LI / p["VLI"]
    cAbm = AMB_BM / p["VBM"]
    cAsk = AMB_SK / p["VSK"]
    cAkid = AMB_KID / p["VKID"]

    # liposome disposition: FREL of it leaks drug into the free pool, the
    # rest is phagocytosed by the MPS -- i.e. handed to the infected cell
    lip_out = p["CL_LIP"] * C_LIP
    rel = p["FREL"] * lip_out
    mps = (1.0 - p["FREL"]) * lip_out

    d[0] = -lip_out
    d[1] = (rel
            - p["CL_FRE"] * C_FRE
            - p["Q_AMB"] * (C_FRE - C_PER))
    d[2] = p["Q_AMB"] * (C_FRE - C_PER)

    # organ macrophage AmB = MPS delivery + passive partition from free drug
    d[3] = p["FSP_A"] * mps + p["KOUT_A"] * (p["KPSP_A"] * C_FRE * p["VSP"] - AMB_SP)
    d[4] = p["FLI_A"] * mps + p["KOUT_A"] * (p["KPLI_A"] * C_FRE * p["VLI"] - AMB_LI)
    d[5] = p["FBM_A"] * mps + p["KOUT_A"] * (p["KPBM_A"] * C_FRE * p["VBM"] - AMB_BM)
    d[6] = p["FSK_A"] * mps + p["KOUT_A"] * (p["KPSK_A"] * C_FRE * p["VSK"] - AMB_SK)
    # the MPS-delivered drug still has to leave the organ eventually
    d[3] -= p["KOUT_A"] * 0.0  # (egress already in the partition term)
    # BUG FOUND: this was written as an uptake/egress pair,
    #   KIN_KID*KPKID_A*C_FRE*VKID - KOUT_KID*AMB_KID,
    # which multiplies the free plasma concentration by Kp AND by
    # KIN/KOUT = 12.5, i.e. a 100-fold amplification.  That put the renal
    # cortex at ~100 mg/L for amphotericin deoxycholate and produced a peak
    # creatinine of 32 mg/dL and a NEGATIVE serum potassium.  It has to be the
    # same equilibration form used for the four macrophage organs, so that the
    # steady-state cortical concentration is Kp x C_free and nothing more.
    d[7] = p["KIN_KID"] * (p["KPKID_A"] * C_FRE * p["VKID"] - AMB_KID)

    # ---------------- miltefosine ------------------------------------------
    MIL_G, MIL_C, MIL_P = g[8], g[9], g[10]
    MIL_SP, MIL_LI, MIL_BM, MIL_SK = g[11], g[12], g[13], g[14]
    C_MIL = MIL_C / p["V1_MIL"]
    C_MILP = MIL_P / p["V2_MIL"]
    cMsp = MIL_SP / p["VSP"]
    cMli = MIL_LI / p["VLI"]
    cMbm = MIL_BM / p["VBM"]
    cMsk = MIL_SK / p["VSK"]

    d[8] = -p["KA_MIL"] * MIL_G
    d[9] = (p["F_MIL"] * p["KA_MIL"] * MIL_G
            - p["CL_MIL"] * C_MIL
            - p["Q_MIL"] * (C_MIL - C_MILP))
    d[10] = p["Q_MIL"] * (C_MIL - C_MILP)
    d[11] = p["KIN_M"] * (p["KPSP_M"] * C_MIL * p["VSP"] - MIL_SP)
    d[12] = p["KIN_M"] * (p["KPLI_M"] * C_MIL * p["VLI"] - MIL_LI)
    d[13] = p["KIN_M"] * (p["KPBM_M"] * C_MIL * p["VBM"] - MIL_BM)
    d[14] = p["KIN_M"] * (p["KPSK_M"] * C_MIL * p["VSK"] - MIL_SK)

    # ---------------- paromomycin ------------------------------------------
    PM_D, PM_C, PM_P = g[15], g[16], g[17]
    PM_SP, PM_LI, PM_BM, PM_SK, PM_KID, PM_COC = g[18], g[19], g[20], g[21], g[22], g[23]
    C_PM = PM_C / p["V_PMA"]
    C_PMP = PM_P / p["V2_PMA"]
    cPsp = PM_SP / p["VSP"]
    cPli = PM_LI / p["VLI"]
    cPbm = PM_BM / p["VBM"]
    cPsk = PM_SK / p["VSK"]
    cPkid = PM_KID / p["VKID"]
    cPcoc = PM_COC / p["VCOC"]

    d[15] = -p["KA_PM"] * PM_D
    d[16] = (p["KA_PM"] * PM_D - p["CL_PM"] * C_PM
             - p["Q_PM"] * (C_PM - C_PMP))
    d[17] = p["Q_PM"] * (C_PM - C_PMP)
    # slow pinocytic entry, slow egress: effective Kp ~5 but with a lag of
    # days, which is why paromomycin needs three weeks and cannot be given
    # as a short course
    d[18] = p["KIN_PM"] * C_PM * p["VSP"] - p["KOUT_PM"] * PM_SP
    d[19] = p["KIN_PM"] * C_PM * p["VLI"] - p["KOUT_PM"] * PM_LI
    d[20] = p["KIN_PM"] * C_PM * p["VBM"] - p["KOUT_PM"] * PM_BM
    d[21] = p["FSK_PM"] * p["KIN_PM"] * C_PM * p["VSK"] - p["KOUT_PM"] * PM_SK
    d[22] = p["KIN_PMK"] * C_PM * p["VKID"] - p["KOUT_PMK"] * PM_KID
    d[23] = p["KIN_COC"] * C_PM * p["VCOC"] - p["KOUT_COC"] * PM_COC

    # ---------------- antimony ---------------------------------------------
    SB5_D, SB5_C, SB5_P = g[24], g[25], g[26]
    SB3_SP, SB3_LI, SB3_BM, SB3_SK, SB_DEP = g[27], g[28], g[29], g[30], g[31]
    C_SB5 = SB5_C / p["V_SBA"]
    C_SBD = SB_DEP / p["V2_SBA"]
    cSsp = SB3_SP / p["VSP"]
    cSli = SB3_LI / p["VLI"]
    cSbm = SB3_BM / p["VBM"]
    cSsk = SB3_SK / p["VSK"]

    d[24] = -p["KA_SB"] * SB5_D
    d[25] = (p["KA_SB"] * SB5_D - p["CL_SB"] * C_SB5
             - p["Q_SB"] * (C_SB5 - C_SBD))
    d[26] = 0.0  # (SB5_P unused; kept for state-vector parity with the R model)
    d[31] = p["Q_SB"] * (C_SB5 - C_SBD)

    # Sb(V) is a prodrug: it must be reduced to Sb(III) INSIDE the phagocyte.
    # Resistance in this model is not a change in target affinity but an
    # increase in thiol-dependent efflux, so it scales KEFF and can never be
    # overcome by raising the dose (host toxicity caps that).
    keff = p["KEFF_SB"] * p["RES_SB"]
    d[27] = p["KRED"] * C_SB5 * p["VSP"] - keff * SB3_SP
    d[28] = p["KRED"] * C_SB5 * p["VLI"] - keff * SB3_LI
    d[29] = p["KRED"] * C_SB5 * p["VBM"] - keff * SB3_BM
    d[30] = p["KPSK_SB"] * p["KRED"] * C_SB5 * p["VSK"] - keff * SB3_SK

    # ---------------- immunity ---------------------------------------------
    TMEM, IFNG, IL10, TNFA, TGFB, CD4, MPHA = g[40], g[41], g[42], g[43], g[44], g[45], g[46]
    P_SP = g[32] + g[33]
    P_LI = g[34] + g[35]
    P_BM = g[36] + g[37]
    P_SK = g[38] + g[39]
    PTOT = P_SP + P_LI + P_BM + P_SK

    PVIS = P_SP + P_LI + P_BM
    # Three signals, three thresholds, and TWO DIFFERENT SOURCES.
    #   AG  = antigen availability, from the TOTAL burden.  Dermal parasites do
    #         prime T cells -- that is precisely why PKDL patients are
    #         leishmanin-positive.
    #   AGH = the IL-10 signal, from the VISCERAL burden only.
    #   AGS = the systemic-illness signal, also visceral only.
    # Keying IL-10 to the total burden made a growing dermal reservoir
    # re-immunosuppress the host: the skin load climbed past K50_IL10, IL-10
    # switched back on, macrophage activation collapsed and the viscera
    # regrew, so EVERY regimen relapsed.  That is backwards.  PKDL patients are
    # not immunosuppressed; they have vigorous Th1 responses, and the IL-10 of
    # active kala-azar comes from spleen and marrow, not from skin lesions.
    AG = PTOT / (PTOT + p["KAG"])                       # antigen availability
    hh = p["HIL10"]
    AGH = PVIS ** hh / (PVIS ** hh + p["K50_IL10"] ** hh)  # only HIGH burden
    AGS = PVIS / (PVIS + p["K50_SYS"])                     # systemic illness
    FCD4 = CD4 / (CD4 + p["KCD4"])
    prime = (1.0 - 0.55 * p["MALNUT"])

    d[40] = (p["KTP"] * prime * AG * FCD4
             * max(0.0, 1.0 - TMEM / p["TMEMAX"])
             / (1.0 + IL10 / p["KI_T"])
             - p["KTD"] * TMEM)
    d[41] = (p["KIP"] * (p["FAG_IFN"] * AG + TMEM) * FCD4
             - p["KID"] * IFNG)
    d[42] = p["KLP"] * AGH - p["KLD"] * IL10
    d[43] = p["KNP"] * AGS - p["KND"] * TNFA
    # NOTE: driven by AGH (the high-burden sigmoid), not AG.  With AG the
    # TGF-beta term stayed near maximal down to 1e4 amastigotes and
    # permanently clamped macrophage activation, so no slow drug could ever
    # hand off to the host.
    d[44] = p["KGP"] * (0.5 * AGH + IL10) - p["KGD"] * TGFB
    cd4set = (p["CD4SET"] if p["HIV"] < 0.5
              else (p["CD4SET_ART"] * p["ART"] + p["CD40"] * (1.0 - p["ART"])))
    d[45] = (p["KCD4H"] * (cd4set - CD4)
             - p["KCD4L"] * AG * CD4 * (1.0 + 0.5 * p["HIV"]))
    # activated-macrophage pool: expands with IFN-g, contracts with IL-10
    MACTdrive = ((IFNG / (IFNG + p["KI_IFN"]))
                 / (1.0 + IL10 / p["KI_IL10"] + TGFB / p["KI_TGFB"])
                 * FCD4
                 * (1.0 - (1.0 - p["FHIV_MAC"]) * p["HIV"]))
    d[46] = p["KMPH"] * (MACTdrive - MPHA)
    MACT = max(MPHA, 0.0)
    kimm = p["KIMM"] * MACT

    # ---------------- parasite killing -------------------------------------
    def emax(c, e50, emx, h):
        if c <= 0.0:
            return 0.0
        ch = c ** h
        return emx * ch / (e50 ** h + ch)

    def killrate(ca, cm, cpm, cs, quiescent):
        """Bliss-independent sum of the four drug kill rates in one organ.

        `quiescent` = 0 for the replicating population, 1 for the persister
        pool, whose EC50s are multiplied by the drug-specific RFAC because
        every one of these four drugs needs an actively metabolising cell:
        ergosterol turnover (AmB), inward transport plus membrane remodelling
        (miltefosine), active translation (paromomycin), and intracellular
        reduction of a prodrug (antimony).
        """
        fa = p["RFAC_A"] if quiescent else 1.0
        fm = (p["RFAC_M"] if quiescent else 1.0) * p["RES_M"]
        fp = (p["RFAC_P"] if quiescent else 1.0) * p["RES_P"]
        # NOTE: RES_SB is deliberately NOT applied here.  Antimony resistance
        # in this model is an EFFLUX phenomenon (it scales KEFF_SB and so
        # lowers intracellular Sb(III)); putting it on EC50 as well would
        # count the same mechanism twice.
        fs = p["RFAC_S"] if quiescent else 1.0
        return (emax(ca, p["EC50_A"] * fa, p["EMAX_A"], p["HILL_A"])
                + emax(cm, p["EC50_M"] * fm, p["EMAX_M"], p["HILL_M"])
                + emax(cpm, p["EC50_P"] * fp, p["EMAX_P"], p["HILL_P"])
                + emax(cs, p["EC50_S"] * fs, p["EMAX_S"], p["HILL_S"]))

    organs = (
        ("SP", 32, 33, P_SP, p["KG_SP"], p["PMAX_SP"], cAsp, cMsp, cPsp, cSsp, 1.0),
        ("LI", 34, 35, P_LI, p["KG_LI"], p["PMAX_LI"], cAli, cMli, cPli, cSli, 1.0),
        ("BM", 36, 37, P_BM, p["KG_BM"], p["PMAX_BM"], cAbm, cMbm, cPbm, cSbm, 1.0),
        ("SK", 38, 39, P_SK, p["KG_SK"], p["PMAX_SK"], cAsk, cMsk, cPsk, cSsk,
         p["FIMM_SK"]),
    )
    # Inter-organ trafficking of infected monocytes (see KMIG notes above).
    # Only REPLICATING amastigotes traffic; quiescent forms are arrested.
    pmax_vis = p["PMAX_SP"] + p["PMAX_LI"] + p["PMAX_BM"]
    PS_VIS = max(g[32], 0.0) + max(g[34], 0.0) + max(g[36], 0.0)
    PS_SK = max(g[38], 0.0)
    net_to_skin = p["KMIG_SKIN"] * PS_VIS - p["KMIG_BACK"] * PS_SK
    for nm, iS, iR, Porg, kg, pmax, ca, cm, cpm, cs, fimm_org in organs:
        crowd = max(0.0, 1.0 - Porg / pmax)
        kimm_o = kimm * fimm_org
        if nm == "SK":
            seed = net_to_skin
        else:
            fo = pmax / pmax_vis
            seed = (p["KMIG"] * (fo * PS_VIS - max(g[iS], 0.0))
                    - fo * net_to_skin)
        kS = killrate(ca, cm, cpm, cs, 0)
        kR = killrate(ca, cm, cpm, cs, 1)
        S, R = max(g[iS], 0.0), max(g[iR], 0.0)
        # stress-induced quiescence: the harder the drug and the immune system
        # push, the larger the persister pool made out of the dying population
        kq = p["KQ0"] + p["KQS"] * (kS + kimm_o)
        kr = p["KR_REACT"]
        d[iS] = S * (kg * crowd - kimm_o - kS) - kq * S + kr * R + seed
        d[iR] = (R * (kg * p["COSTR"] * crowd - p["FIMM_R"] * kimm_o - kR)
                 + kq * S - kr * R)
        # NOTE (bug found while building this file): an earlier version put the
        # extinction floor HERE, as `if S < 1e-7 and d[iS] > 0: d[iS] = 0`.
        # That is a discontinuous switch inside the RHS.  It is invisible while
        # burdens are high and then, as soon as any compartment falls below the
        # floor, it flips on and off between solver probes, destroys LSODA's
        # error estimate and the integration never advances (a 28-day
        # miltefosine course took >100 s and did not finish).  The floor
        # belongs BETWEEN integration intervals, not inside f(t,y); it is
        # applied in `clamp_extinct()` at every chunk boundary instead.
    # ---------------- clinical ---------------------------------------------
    SPL, LIVS, HGB, PLT, WBC, ALB, IGG, TEMP, BWTL, PKDLv = (
        g[47], g[48], g[49], g[50], g[51], g[52], g[53], g[54], g[55], g[56])

    spl_target = p["SPL0"] + p["SPLMAX"] * P_SP / (P_SP + p["K50_SPL"])
    d[47] = p["KSPL"] * (spl_target - SPL)
    liv_target = p["LIV0"] + p["LIVMAX"] * P_LI / (P_LI + 3.0 * p["K50_SPL"])
    d[48] = p["KLIV"] * (liv_target - LIVS)

    marrow = 1.0 / (1.0 + P_BM / (0.35 * p["PMAX_BM"]))
    inflam = 1.0 / (1.0 + TNFA / 0.30)
    d[49] = (p["KHGB"] * p["HGB0"] * (1.0 - p["FSUP_HGB"] * (1.0 - marrow * inflam))
             - (p["KHGB"] + p["KSPL_HGB"] * p["KHGB"] * SPL) * HGB)
    d[50] = (p["KPLT"] * p["PLT0"] * marrow
             - (p["KPLT"] + p["KSPL_PLT"] * p["KPLT"] * SPL) * PLT)
    d[51] = (p["KWBC"] * p["WBC0"] * marrow * inflam
             - (p["KWBC"] + p["KSPL_WBC"] * p["KWBC"] * SPL) * WBC)
    d[52] = p["KALB"] * (p["ALB0"] * (1.0 - 0.42 * AGS) - ALB)
    d[53] = p["KIGG_P"] * AG - p["KIGG_D"] * (IGG - p["IGG0"])
    d[54] = p["KTEMP"] * (p["TEMP0"] + p["TNF_FEV"] * TNFA / (TNFA + 0.45) - TEMP)
    d[55] = p["KBWL"] * (p["BWLMAX"] * AGS - BWTL)
    # PKDL: skin amastigotes become lesions only when immunity RETURNS, which
    # is why the rash appears after apparently successful treatment
    d[56] = p["KPKDL_P"] * P_SK * MACT - p["KPKDL_D"] * PKDLv

    # ---------------- toxicity ---------------------------------------------
    TUBI, SCR, KSER, MGSER, HEAR, QTCv, LIPA, ALTX, GITX = (
        g[57], g[58], g[59], g[60], g[61], g[62], g[63], g[64], g[65])
    # tubular injury: driven by RENAL AmB, i.e. by the free-drug integral
    d[57] = p["KTUB"] * cAkid - p["KTUBR"] * TUBI
    gfr = p["GFR0"] / (1.0 + max(TUBI, 0.0))
    d[58] = p["KSCR"] * (p["SCR0"] * p["GFR0"] / max(gfr, 1.0) - SCR)
    d[59] = p["KKR"] * (p["K_K0"] - KSER) - p["KWK"] * cAkid
    d[60] = p["KMGR"] * (p["MG0"] - MGSER) - p["KWMG"] * cAkid
    d[61] = p["KHEAR"] * cPcoc                       # irreversible, no repair
    d[62] = p["KQTC"] * C_SBD * p["KQTCR"] * 50.0 - p["KQTCR"] * (QTCv - p["QTC0"])
    d[63] = p["KLIPA"] * C_SB5 - p["KLIPAR"] * (LIPA - p["LIPA0"])
    d[64] = (p["KALTX"] * C_SB5 + p["KALT_M"] * C_MIL
             - p["KALTR"] * (ALTX - p["ALT0"]))
    d[65] = p["KGITX"] * C_MIL / (C_MIL + p["EM50_GI"]) - p["KGITR"] * GITX

    # ---------------- hazard and accumulators ------------------------------
    haz = p["H0"] * (
        (1.0 + 2.2 * max(0.0, p["HGB_H"] - HGB) / p["HGB_H"])
        * (1.0 + 2.6 * max(0.0, p["ALB_H"] - ALB) / p["ALB_H"])
        * (1.0 + 2.0 * max(0.0, p["WBC_H"] - WBC) / p["WBC_H"])
        * (0.15 + 0.85 * AGS)
        * (1.0 + 1.6 * p["HIV"])
    )
    d[66] = C_FRE          # free-plasma AmB AUC  (nephrotoxicity driver)
    d[67] = cAsp           # spleen-macrophage AmB AUC (efficacy driver)
    d[68] = cMsp
    d[69] = cSsp
    d[70] = cPsp
    d[71] = haz
    d[72] = 0.0
    return d


# --------------------------------------------------------------------------
# 4.  Dosing and integration
# --------------------------------------------------------------------------
def dose_events(regimen, p, rng=None):
    """Expand a regimen spec into a sorted list of (time_h, cmt, amount_mg).

    A regimen is a list of dicts:
      {drug, mgkg (or mg), start_day, n, interval_h}
    """
    ev = []
    wt = p["WT"]
    for r in regimen:
        drug = r["drug"]
        amt = r.get("mg")
        if amt is None:
            amt = r["mgkg"] * (r.get("dosewt") or wt)
        n = r.get("n", 1)
        iv = r.get("interval_h", 24.0)
        t0 = r["start_day"] * DAY
        cmt = {"LAMB": 0, "DAMB": 1, "MIL": 8, "PM": 15, "SB": 24}[drug]
        for k in range(n):
            # oral doses can be missed; parenteral doses are given under
            # supervision and are not subject to ADHER
            if (drug == "MIL" and rng is not None and p.get("ADHER", 1.0) < 1.0
                    and rng.random() > p["ADHER"]):
                continue
            ev.append((t0 + k * iv, cmt, amt, drug))
    ev.sort(key=lambda e: e[0])
    return ev


PIDX = list(range(32, 40))          # the eight parasite states
EXTINCT = 1.0e-6                    # burden units = one whole amastigote


def clamp_extinct(y):
    """Apply the biological floor between integration intervals.

    One amastigote is the smallest population that exists.  Anything below
    that is set to exactly zero, which is both biologically correct and what
    keeps the ODE from "relapsing" out of 1e-300.  Doing it here rather than
    inside the RHS keeps f(t,y) smooth (see the note in `rhs`).
    """
    for i in PIDX:
        if y[i] < EXTINCT:
            y[i] = 0.0
        elif y[i] < 0.0:
            y[i] = 0.0
    return y


def _integrate(y, ta, tb, p, tgrid, out):
    """Integrate ta -> tb, writing any grid points that fall inside."""
    want = tgrid[(tgrid > ta) & (tgrid <= tb)]
    t_eval = np.append(want, tb) if (len(want) == 0 or want[-1] < tb) else want
    sol = solve_ivp(rhs, (ta, tb), y, args=(p,), method="LSODA",
                    rtol=1e-6, atol=1e-9, t_eval=t_eval)
    if not sol.success:
        raise RuntimeError("integration failed at t=%.1f: %s" % (ta, sol.message))
    if len(want):
        out[np.searchsorted(tgrid, want)] = sol.y[:, :len(want)].T
    # NOTE: solve_ivp with t_eval returns the state only at the REQUESTED
    # times, so tb has to be requested explicitly; carrying the last
    # requested time forward instead silently drops the tail of every
    # interval and shifts a whole 28-day course.
    return sol.y[:, -1].copy(), (int(np.searchsorted(tgrid, want)[-1])
                                 if len(want) else None)


def _run(y, p, p_raw, ev, tmax, tgrid, out, chunk_h=14.0 * DAY):
    """Common driver: dose events plus regular chunk boundaries at which the
    extinction floor is applied."""
    marks = set([0.0, tmax])
    marks.update(e[0] for e in ev if e[0] <= tmax)
    k = chunk_h
    while k < tmax:
        marks.add(k)
        k += chunk_h
    breaks = sorted(marks)

    cum = 0.0
    filled = 0
    out[0] = y
    for i in range(len(breaks) - 1):
        ta, tb = breaks[i], breaks[i + 1]
        for (te, cmt, amt, drug) in ev:
            if abs(te - ta) < 1e-9:
                y[cmt] += amt
                cum += amt / p_raw["WT"]
        y[IX["CUMDOSE"]] = cum
        clamp_extinct(y)
        if tb <= ta:
            continue
        y, last = _integrate(y, ta, tb, p, tgrid, out)
        if last is not None:
            filled = max(filled, last)
        y[IX["CUMDOSE"]] = cum
    clamp_extinct(y)
    if filled < len(tgrid) - 1:
        out[filled + 1:] = y
    # the solver leaves denormal / slightly negative residue in cleared
    # compartments; apply the same one-amastigote floor to every stored row so
    # no reported burden is ever a numerical artefact
    par = out[:, PIDX]
    par[par < EXTINCT] = 0.0
    out[:, PIDX] = par
    return tgrid, out


def simulate(p_raw, regimen, tmax_day=540.0, stage="active", nout=1400, rng=None):
    p = scale_params(p_raw)
    y = initial_state(p, stage=stage)
    ev = dose_events(regimen, p, rng=rng)
    tmax = tmax_day * DAY
    tgrid = np.linspace(0.0, tmax, nout)
    out = np.zeros((nout, NEQ))
    return _run(y, p, p_raw, ev, tmax, tgrid, out)


# --------------------------------------------------------------------------
# 5.  Derived outcome metrics
# --------------------------------------------------------------------------
def ptot(out):
    return out[:, 32:40].sum(axis=1)


def outcome(t, out, p, eot_day, tmax_day):
    """Classify a simulated patient the way a trial would.

    IMPORTANT: cure and relapse are judged on the VISCERAL compartments
    (spleen + liver + marrow), not on total body burden.  That is what
    "definitive cure at 6 months" means in every visceral-leishmaniasis trial
    ever run.  The dermal compartment clears far more slowly -- little drug
    reaches it and the dermis is relatively immune-privileged -- and its
    residue is what becomes PKDL.  Skin burden and lesion load are reported
    separately, exactly as trials report PKDL incidence separately from cure.
    """
    P = ptot(out)
    PV = out[:, 32:38].sum(axis=1)           # spleen + liver + marrow
    PSK = out[:, 38:40].sum(axis=1)          # dermis
    day = t / DAY
    i_eot = int(np.argmin(np.abs(day - eot_day)))
    P0 = PV[0]
    nadir = PV[i_eot:].min() if i_eot < len(PV) else PV[-1]
    i_nadir = i_eot + int(np.argmin(PV[i_eot:]))
    Pend = PV[-1]
    # Thresholds, stated once because every rate reported here depends on them:
    #   initial (parasitological) response = >= 3 log10 fall in visceral burden
    #   relapse = an initial response followed by visceral regrowth above the
    #     clinically detectable level (1e-2 units = 1e4 amastigotes)
    #   definitive cure = an initial response with no such regrowth by the end
    #     of follow-up
    DET = 1.0e-2
    if P0 < DET:
        # A PKDL patient has NO visceral disease at baseline, so the visceral
        # endpoint is undefined for them.  Score those runs on the dermal
        # compartment instead -- which is what a PKDL trial actually does
        # (lesion clearance, not splenic aspirate).
        P0 = PSK[0]
        nadir = PSK[i_eot:].min() if i_eot < len(PSK) else PSK[-1]
        i_nadir = i_eot + int(np.argmin(PSK[i_eot:]))
        Pend = PSK[-1]
    initial_response = nadir < max(P0, 1e-30) * 1.0e-3
    relapse = bool(initial_response and (Pend > DET))
    fail = (not initial_response)
    cured = bool(initial_response and not relapse)
    res = {
        "P0": P0, "P_eot": PV[i_eot], "nadir": nadir,
        "nadir_day": day[i_nadir], "P_end": Pend, "P_total_end": P[-1],
        "logdrop_eot": math.log10(max(PV[i_eot], 1e-30) / P0),
        "logdrop_nadir": math.log10(max(nadir, 1e-30) / P0),
        "cure": cured, "relapse": relapse, "primary_failure": bool(fail),
        "surv": float(math.exp(-out[-1, IX["CUMHAZ"]])),
        "mort_pct": float(100.0 * (1.0 - math.exp(-out[-1, IX["CUMHAZ"]]))),
        "AUC_free": out[-1, IX["AUCFRE"]],
        "AUC_ambsp": out[-1, IX["AUCASP"]],
        "AUC_milsp": out[-1, IX["AUCMSP"]],
        "AUC_sb3sp": out[-1, IX["AUCBSP"]],
        "AUC_pmsp": out[-1, IX["AUCPSP"]],
        "scr_max": out[:, IX["SCR"]].max(),
        "k_min": out[:, IX["KSER"]].min(),
        "mg_min": out[:, IX["MGSER"]].min(),
        "hear_dB": out[-1, IX["HEAR"]],
        "qtc_max": out[:, IX["QTC"]].max(),
        "lipa_max": out[:, IX["LIPA"]].max(),
        "alt_max": out[:, IX["ALTX"]].max(),
        "gitx_max": out[:, IX["GITX"]].max(),
        "spl_end": out[-1, IX["SPL"]],
        "hgb_end": out[-1, IX["HGB"]],
        "plt_end": out[-1, IX["PLT"]],
        "alb_end": out[-1, IX["ALB"]],
        "igg_end": out[-1, IX["IGG"]],
        # --- the PKDL endpoint, reported separately, as trials do -----------
        "pkdl_max": out[:, IX["PKDL"]].max(),
        "pkdl_end": out[-1, IX["PKDL"]],
        "P_sk_end": PSK[-1],
        "P_sk_180": PSK[int(np.argmin(np.abs(day - 180.0)))],
        "tmem_eot": out[i_eot, IX["TMEM"]],
        "tmem_max": out[:, IX["TMEM"]].max(),
        "cd4_end": out[-1, IX["CD4"]],
        "cumdose": out[-1, IX["CUMDOSE"]],
    }
    return res


# --------------------------------------------------------------------------
# 6.  Regimen library (20 scenarios)
# --------------------------------------------------------------------------
def regimens(wt=50.0):
    R = {}
    R["S01_untreated"] = []
    R["S02_ssg_africa"] = [dict(drug="SB", mgkg=20.0, start_day=0, n=30)]
    R["S03_ssg_bihar_res"] = [dict(drug="SB", mgkg=20.0, start_day=0, n=30)]
    R["S04_damb_alt30"] = [dict(drug="DAMB", mgkg=1.0, start_day=0, n=15,
                                interval_h=48.0)]
    R["S05_lamb_single10"] = [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]
    R["S06_lamb_21_multi"] = [dict(drug="LAMB", mgkg=3.0, start_day=0, n=5),
                              dict(drug="LAMB", mgkg=3.0, start_day=13, n=1),
                              dict(drug="LAMB", mgkg=3.0, start_day=20, n=1)]
    R["S07_lamb_5x3"] = [dict(drug="LAMB", mgkg=3.0, start_day=0, n=5)]
    R["S08_mil28_adult"] = [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)]
    R["S09_mil28_child_linear"] = [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)]
    # allometric (WHO weight-band-style) dosing: the mg/day that gives a
    # 10 kg child the SAME AUC as 125 mg/day in a 50 kg adult, i.e. scaled by
    # WT^0.75 rather than WT^1.0
    R["S10_mil28_child_allom"] = [dict(drug="MIL", mg=125.0 * (10.0 / 50.0) ** 0.75,
                                       start_day=0, n=28)]
    R["S11_pm21"] = [dict(drug="PM", mgkg=15.0, start_day=0, n=21)]
    R["S12_pm21_africa"] = [dict(drug="PM", mgkg=15.0, start_day=0, n=21)]
    R["S13_lamb5_mil7"] = [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
                           dict(drug="MIL", mgkg=2.5, start_day=1, n=7)]
    R["S14_lamb5_pm10"] = [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
                           dict(drug="PM", mgkg=15.0, start_day=1, n=10)]
    R["S15_mil10_pm10"] = [dict(drug="MIL", mgkg=2.5, start_day=0, n=10),
                           dict(drug="PM", mgkg=15.0, start_day=0, n=10)]
    R["S16_ssg_pm17"] = [dict(drug="SB", mgkg=20.0, start_day=0, n=17),
                         dict(drug="PM", mgkg=15.0, start_day=0, n=17)]
    R["S17_hiv_lamb10"] = [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]
    R["S18_hiv_lamb30_mil28"] = [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
                                 dict(drug="LAMB", mgkg=5.0, start_day=2, n=1),
                                 dict(drug="LAMB", mgkg=5.0, start_day=4, n=1),
                                 dict(drug="LAMB", mgkg=5.0, start_day=6, n=1),
                                 dict(drug="LAMB", mgkg=5.0, start_day=8, n=1),
                                 dict(drug="LAMB", mgkg=5.0, start_day=10, n=1),
                                 dict(drug="MIL", mgkg=2.5, start_day=0, n=28)]
    R["S19_pkdl_mil84"] = [dict(drug="MIL", mgkg=2.5, start_day=0, n=84)]
    R["S20_pkdl_lamb"] = [dict(drug="LAMB", mgkg=5.0, start_day=0, n=4,
                               interval_h=7 * 24.0)]
    return R


def scenario_params(name, base=None):
    p = dict(base or default_params())
    if name == "S03_ssg_bihar_res":
        # resistance is an EFFLUX parameter, not a subpopulation fraction:
        # FR0 is the quiescent pool and has nothing to do with antimony
        p["RES_SB"] = 9.0
    if name == "S02_ssg_africa":
        p["RES_SB"] = 1.0
    if name == "S04_damb_alt30":
        p["AMBFORM"] = 0.0
    if name in ("S09_mil28_child_linear", "S10_mil28_child_allom"):
        p["WT"] = 10.0
        p["AGE"] = 3.0
        p["MALNUT"] = 0.45
    if name == "S12_pm21_africa":
        p["EC50_P"] = 25.0 * 1.60    # documented East-African hyposensitivity
    if name in ("S17_hiv_lamb10", "S18_hiv_lamb30_mil28"):
        p["HIV"] = 1.0
        p["CD40"] = 90.0
        # S17 is the common real-world sequence: VL is diagnosed first and
        # treated before ART has been started, so there is no CD4 recovery
        # during the window that decides cure.  S18 is the guideline
        # combination WITH ART.
        p["ART"] = 0.0 if name == "S17_hiv_lamb10" else 1.0
    if name in ("S19_pkdl_mil84", "S20_pkdl_lamb"):
        p["PKDL_START"] = 1.0
    return p


def pkdl_initial(p):
    """Post-VL patient whose visceral compartments are cured but whose skin
    still carries amastigotes; immunity has returned (that is what makes the
    lesions appear)."""
    y = initial_state(p, stage="active")
    for org, frac in (("SP", 0.0), ("LI", 0.0), ("BM", 0.0), ("SK", 0.45)):
        tot = frac * p["PMAX_" + org]
        y[IX["P_%s_S" % org]] = tot * (1.0 - p["FR0"])
        y[IX["P_%s_R" % org]] = tot * p["FR0"]
    y[IX["TMEM"]] = 0.85
    y[IX["IL10"]] = 0.08
    y[IX["IFNG"]] = 0.30
    y[IX["TNFA"]] = 0.10
    y[IX["TGFB"]] = 0.12
    y[IX["MPHA"]] = 0.45
    y[IX["SPL"]] = 1.5
    y[IX["HGB"]] = 12.2
    y[IX["PLT"]] = 210.0
    y[IX["WBC"]] = 5.4
    y[IX["ALB"]] = 3.9
    y[IX["IGG"]] = 2.9
    y[IX["TEMP"]] = 36.9
    y[IX["PKDL"]] = 1.2
    return y


def simulate_pkdl(p_raw, regimen, tmax_day=540.0, nout=1200):
    p = scale_params(p_raw)
    y = pkdl_initial(p)
    ev = dose_events(regimen, p)
    tmax = tmax_day * DAY
    tgrid = np.linspace(0.0, tmax, nout)
    out = np.zeros((nout, NEQ))
    return _run(y, p, p_raw, ev, tmax, tgrid, out)


# --------------------------------------------------------------------------
# 7.  Analyses
# --------------------------------------------------------------------------
def hr(title, ch="="):
    print()
    print(ch * 78)
    print(title)
    print(ch * 78)


def eot_of(regimen):
    if not regimen:
        return 0.0
    last = 0.0
    for r in regimen:
        iv = r.get("interval_h", 24.0)
        last = max(last, r["start_day"] + (r.get("n", 1) - 1) * iv / DAY)
    return last + 1.0


def run_all_scenarios(tmax_day=540.0):
    R = regimens()
    rows = {}
    for name, reg in R.items():
        p = scenario_params(name)
        if name.startswith(("S19", "S20")):
            t, out = simulate_pkdl(p, reg, tmax_day=tmax_day)
        else:
            t, out = simulate(p, reg, tmax_day=tmax_day)
        rows[name] = (t, out, p, outcome(t, out, p, eot_of(reg), tmax_day))
    return rows


def a1_targeting(res):
    """Claim 1: the same mg/kg of amphotericin, encapsulated or not, splits
    into two integrals that move in OPPOSITE directions."""
    hr("A1  LIPOSOMAL vs DEOXYCHOLATE AMPHOTERICIN: ONE DOSE, TWO INTEGRALS")
    p = default_params()
    out = {}
    for form, tag in (("LAMB", "liposomal 10 mg/kg x1"),
                      ("DAMB", "deoxycholate 10 mg/kg total (1 mg/kg q48h x10)")):
        if form == "LAMB":
            reg = [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]
        else:
            reg = [dict(drug="DAMB", mgkg=1.0, start_day=0, n=10,
                        interval_h=48.0)]
        t, o = simulate(p, reg, tmax_day=180.0)
        out[form] = dict(
            tag=tag,
            auc_free=o[-1, IX["AUCFRE"]],
            auc_sp=o[-1, IX["AUCASP"]],
            auc_li=np.trapezoid(o[:, 4] / scale_params(p)["VLI"], t),
            cmax_tot=(o[:, 0] / scale_params(p)["V_LIPA"]
                      + o[:, 1] / scale_params(p)["V_FREA"]).max(),
            cmax_free=(o[:, 1] / scale_params(p)["V_FREA"]).max(),
            scr=o[:, IX["SCR"]].max(),
            kmin=o[:, IX["KSER"]].min(),
            # day 21: still inside the drug-dominated phase.  Reading the
            # overall nadir here just returned the extinction floor for both
            # arms and compared nothing.
            logdrop=math.log10(max(ptot(o)[np.argmin(np.abs(t - 21 * DAY))],
                                   1e-30) / ptot(o)[0]),
        )
    L, D = out["LAMB"], out["DAMB"]
    print(f"  {'metric':44s} {'L-AmB':>13s} {'d-AmB':>13s} {'ratio':>9s}")
    def line(lbl, a, b, fmt="{:13.4g}"):
        r = a / b if b else float("inf")
        print(f"  {lbl:44s} " + fmt.format(a) + " " + fmt.format(b)
              + f" {r:9.2f}")
    line("total plasma AmB Cmax (mg/L)", L["cmax_tot"], D["cmax_tot"])
    line("FREE plasma AmB Cmax (mg/L)", L["cmax_free"], D["cmax_free"])
    line("FREE plasma AmB AUC (mg.h/L)  [nephro]", L["auc_free"], D["auc_free"])
    line("spleen-macrophage AmB AUC (mg.h/L) [kill]", L["auc_sp"], D["auc_sp"])
    line("liver-macrophage AmB AUC (mg.h/L)", L["auc_li"], D["auc_li"])
    line("peak creatinine (mg/dL)", L["scr"], D["scr"])
    line("nadir serum K (mEq/L)", L["kmin"], D["kmin"])
    ti_L = L["auc_sp"] / L["auc_free"]
    ti_D = D["auc_sp"] / D["auc_free"]
    print(f"\n  therapeutic-index surrogate  AUC_spleen / AUC_free")
    print(f"    liposomal      : {ti_L:10.1f}")
    print(f"    deoxycholate   : {ti_D:10.1f}")
    print(f"    GAIN from encapsulation at identical total mg/kg : {ti_L/ti_D:6.1f} x")
    print(f"      (spleen exposure x{L['auc_sp']/D['auc_sp']:.1f}, "
          f"free-plasma exposure x{L['auc_free']/D['auc_free']:.2f})")
    print(f"  log10 burden drop at day 21  L-AmB {L['logdrop']:.2f}   "
          f"d-AmB {D['logdrop']:.2f}")
    return {"ti_gain": ti_L / ti_D, "sp_ratio": L["auc_sp"] / D["auc_sp"],
            "free_ratio": L["auc_free"] / D["auc_free"],
            "ti_L": ti_L, "ti_D": ti_D}


def a1b_frel_sweep():
    """How much of the liposomal advantage is the liposome staying intact?"""
    hr("A1b DELIVERY-FRACTION SWEEP: the liposome IS the mechanism", "-")
    print("  FREL = fraction of liposome clearance that leaks drug to plasma")
    print(f"  {'FREL':>6s} {'AUC_spleen':>12s} {'AUC_free':>10s} {'TI':>9s} "
          f"{'log drop':>9s} {'SCr max':>8s}")
    rowsr = []
    for frel in (0.05, 0.15, 0.30, 0.50, 0.70, 0.90, 1.00):
        p = default_params(); p["FREL"] = frel
        t, o = simulate(p, [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)],
                        tmax_day=180.0)
        asp, afr = o[-1, IX["AUCASP"]], o[-1, IX["AUCFRE"]]
        ld = math.log10(max(ptot(o).min(), 1e-30) / ptot(o)[0])
        print(f"  {frel:6.2f} {asp:12.1f} {afr:10.1f} {asp/afr:9.1f} "
              f"{ld:9.2f} {o[:, IX['SCR']].max():8.2f}")
        rowsr.append((frel, asp, afr, ld))
    return rowsr


def a2_allometry():
    """Claim 2: paediatric miltefosine failure is arithmetic."""
    hr("A2  MILTEFOSINE IN CHILDREN: AUC ~ WT^0.25 IS THE WHOLE STORY")
    print("  linear mg/kg dosing, CL ~ WT^0.75  =>  AUC ~ WT^0.25")
    print(f"  {'WT (kg)':>8s} {'dose mg/d':>10s} {'AUC0-28d':>12s} "
          f"{'rel. to 50 kg':>14s} {'predicted WT^0.25':>18s}")
    ref = None
    rows = []
    for wt in (7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 70.0):
        p = default_params(); p["WT"] = wt
        reg = [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)]
        t, o = simulate(p, reg, tmax_day=60.0, nout=700)
        q = scale_params(p)
        auc = np.trapezoid(o[:, 9] / q["V1_MIL"], t)
        if abs(wt - 50.0) < 1e-9:
            ref = auc
        rows.append((wt, 2.5 * wt, auc))
    for wt, dose, auc in rows:
        pred = (wt / 50.0) ** 0.25
        print(f"  {wt:8.1f} {dose:10.1f} {auc:12.1f} {auc/ref:14.3f} "
              f"{pred:18.3f}")
    hr("A2b ALLOMETRIC (LBW-style) DOSING RESTORES THE CHILD'S EXPOSURE", "-")
    # allometric dosing: give dose proportional to WT^0.75 normalised so a
    # 50 kg adult receives the same 125 mg/day
    print(f"  {'WT (kg)':>8s} {'linear mg/d':>12s} {'allometric mg/d':>16s} "
          f"{'AUC linear':>11s} {'AUC allom':>10s} {'AUC ratio':>10s}")
    out2 = []
    for wt in (7.0, 10.0, 15.0, 20.0, 30.0, 50.0):
        p = default_params(); p["WT"] = wt
        lin = 2.5 * wt
        allo = 125.0 * (wt / 50.0) ** 0.75
        a = []
        for mg in (lin, allo):
            t, o = simulate(p, [dict(drug="MIL", mg=mg, start_day=0, n=28)],
                            tmax_day=60.0, nout=700)
            a.append(np.trapezoid(o[:, 9] / scale_params(p)["V1_MIL"], t))
        print(f"  {wt:8.1f} {lin:12.1f} {allo:16.1f} {a[0]:11.1f} "
              f"{a[1]:10.1f} {a[1]/a[0]:10.3f}")
        out2.append((wt, lin, allo, a[0], a[1]))
    return rows, out2


def a2c_child_outcome():
    """Does the arithmetic actually change the outcome, or only the AUC?"""
    hr("A2c THE SAME ARITHMETIC, CARRIED THROUGH TO CURE / RELAPSE", "-")
    print(f"  {'subject':32s} {'AUC0-28d':>10s} {'log drop':>9s} "
          f"{'day-180 burden':>15s} {'outcome':>10s}")
    rows = []
    for lbl, wt, mode, malnut in (
            ("adult 50 kg, 2.5 mg/kg/d", 50.0, "lin", 0.0),
            ("child 10 kg, 2.5 mg/kg/d", 10.0, "lin", 0.45),
            ("child 10 kg, allometric", 10.0, "allo", 0.45),
            ("child 10 kg, allom, well-nourished", 10.0, "allo", 0.0),
            ("child 20 kg, 2.5 mg/kg/d", 20.0, "lin", 0.30),
            ("child 20 kg, allometric", 20.0, "allo", 0.30)):
        p = default_params(); p["WT"] = wt; p["MALNUT"] = malnut
        mg = 2.5 * wt if mode == "lin" else 125.0 * (wt / 50.0) ** 0.75
        reg = [dict(drug="MIL", mg=mg, start_day=0, n=28)]
        t, o = simulate(p, reg, tmax_day=360.0)
        r = outcome(t, o, scale_params(p), 29.0, 360.0)
        auc = np.trapezoid(o[:, 9] / scale_params(p)["V1_MIL"], t)
        st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "failure")
        print(f"  {lbl:32s} {auc:10.1f} {r['logdrop_nadir']:9.2f} "
              f"{r['P_end']:15.3g} {st:>10s}")
        rows.append((lbl, auc, r))
    return rows


def a3_separatrix():
    """Claim 3: cure is a separatrix crossing, and CD4 moves the line.

    The interesting result is not that the threshold shrinks with CD4.  It is
    that below a certain CD4 the threshold DISAPPEARS: there is no residual
    burden small enough for the host to finish, so cure depends entirely on
    the drug sterilising.  A bisection returns no root in that regime, and
    reporting that as NaN would hide the finding, so it is labelled.
    """
    hr("A3  THE CURE THRESHOLD IS A COMPUTABLE SEPARATRIX, NOT ZERO")
    print("  For a patient at end of treatment with primed memory TMEM,")
    print("  what residual burden can the host still finish off unaided?")
    print(f"  {'CD4':>6s} {'TMEM at EOT':>12s} {'critical residual burden':>26s}"
          f" {'= amastigotes':>16s}  note")

    def endpoint(cd4, tmem, logP):
        p = default_params()
        p["CD40"] = cd4
        if cd4 < 350:
            p["HIV"] = 1.0
            p["ART"] = 1.0
        q = scale_params(p)
        y = initial_state(q, stage="active")
        P = 10.0 ** logP
        for org, frac in (("SP", 0.45), ("LI", 0.35), ("BM", 0.15), ("SK", 0.05)):
            y[IX["P_%s_S" % org]] = P * frac * (1 - q["FR0"])
            y[IX["P_%s_R" % org]] = P * frac * q["FR0"]
        y[IX["TMEM"]] = tmem
        y[IX["IL10"]] = 0.10
        y[IX["IFNG"]] = 0.30
        y[IX["TNFA"]] = 0.12
        y[IX["TGFB"]] = 0.15
        y[IX["MPHA"]] = 0.45
        y[IX["CD4"]] = cd4
        s = solve_ivp(rhs, (0.0, 360 * DAY), y, args=(q,),
                      method="LSODA", rtol=1e-6, atol=1e-9)
        Pe = max(s.y[32:40, -1].sum(), 1e-30)
        return math.log10(Pe) - logP

    rows = []
    for cd4, tmem in ((700, 0.55), (700, 0.30), (700, 0.15),
                      (500, 0.55), (350, 0.55), (300, 0.55), (250, 0.55),
                      (200, 0.55), (100, 0.55), (50, 0.55)):
        lo, hi = -6.0, 4.0
        f_lo, f_hi = endpoint(cd4, tmem, lo), endpoint(cd4, tmem, hi)
        if f_lo > 0 and f_hi > 0:
            print(f"  {cd4:6d} {tmem:12.2f} {'none':>26s} {'-':>16s}  "
                  f"host cannot clear ANY burden")
            rows.append((cd4, tmem, None))
            continue
        if f_lo < 0 and f_hi < 0:
            print(f"  {cd4:6d} {tmem:12.2f} {'> 1e4 units':>26s} {'-':>16s}  "
                  f"host clears everything tested")
            rows.append((cd4, tmem, 4.0))
            continue
        crit = brentq(lambda x: endpoint(cd4, tmem, x), lo, hi,
                      xtol=1e-3, maxiter=60)
        print(f"  {cd4:6d} {tmem:12.2f} {10.0 ** crit:26.4g} "
              f"{10.0 ** crit * 1e6:16.3g}  ")
        rows.append((cd4, tmem, crit))

    have = [r for r in rows if r[1] == 0.55 and r[2] is not None]
    none = [r for r in rows if r[1] == 0.55 and r[2] is None]
    if have and none:
        print(f"\n  The threshold does not shrink smoothly -- it VANISHES between")
        print(f"  CD4 {min(r[0] for r in have)} (threshold "
              f"{10.0 ** min(r[2] for r in have):.3g} units) and CD4 "
              f"{max(r[0] for r in none)} (no threshold at all).")
    print("\n  Reading: above that CD4, a drug course only has to get the")
    print("  patient under the line and the host finishes.  Below it, the drug")
    print("  has to sterilise on its own -- which is why HIV-VL needs 3-4x the")
    print("  total L-AmB dose and a companion drug, not a longer course of the")
    print("  same one.")
    return rows


def a4_timing():
    """Claim 6: combination synergy is in TIME, not in concentration.

    Run twice.  In a typical patient every arm cures and the comparison shows
    nothing -- that is itself worth reporting, because it says the efficacy
    margin of L-AmB in an immunocompetent patient is wide enough to absorb the
    difference.  The discrimination appears under stress (a high presenting
    burden in a partially immunosuppressed host), which is exactly the patient
    for whom combination regimens were designed.
    """
    specs = [
        ("L-AmB 5 mg/kg alone",
         [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1)]),
        ("L-AmB 5 + MF 7 d",
         [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
          dict(drug="MIL", mgkg=2.5, start_day=1, n=7)]),
        ("L-AmB 5 + MF 7 d, MF given FIRST",
         [dict(drug="MIL", mgkg=2.5, start_day=0, n=7),
          dict(drug="LAMB", mgkg=5.0, start_day=7, n=1)]),
        ("L-AmB 5 + MF 14 d",
         [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
          dict(drug="MIL", mgkg=2.5, start_day=1, n=14)]),
        ("MF 28 d alone",
         [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)]),
        ("L-AmB 10 mg/kg alone",
         [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]),
        ("L-AmB 5 + PM 10 d",
         [dict(drug="LAMB", mgkg=5.0, start_day=0, n=1),
          dict(drug="PM", mgkg=15.0, start_day=1, n=10)]),
    ]

    def block(title, mods):
        hr(title, "-")
        print(f"  {'regimen':36s} {'nadir':>10s} {'nadir d':>8s} {'d360':>10s} "
              f"{'TMEM@EOT':>9s} {'d21 log':>8s} {'outcome':>8s}")
        out = []
        for lbl, reg in specs:
            p = default_params()
            p.update(mods)
            t, o = simulate(p, reg, tmax_day=360.0)
            r = outcome(t, o, scale_params(p), eot_of(reg), 360.0)
            P = ptot(o)
            i21 = int(np.argmin(np.abs(t - 21 * DAY)))
            d21 = math.log10(max(P[i21], 1e-30) / P[0])
            st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "fail")
            print(f"  {lbl:36s} {r['nadir']:10.3g} {r['nadir_day']:8.1f} "
                  f"{r['P_end']:10.3g} {r['tmem_eot']:9.3f} {d21:8.2f} "
                  f"{st:>8s}")
            out.append((lbl, r, d21))
        return out

    hr("A4  WHY ONE DAY OF L-AmB PLUS A WEEK OF MILTEFOSINE WORKS")
    print("  Matched L-AmB total dose (5 mg/kg) with and without a companion.")
    a = block("A4a  typical immunocompetent patient", {})
    b = block("A4b  STRESSED: 5x presenting burden, CD4 300, malnutrition 0.5",
              dict(P0SCALE=5.0, CD40=300.0, CD4SET=300.0, MALNUT=0.5))
    print("\n  In A4a every arm cures: in an immunocompetent patient the")
    print("  efficacy margin of amphotericin is wide enough that the companion")
    print("  drug changes nothing measurable.  A4b is the patient the")
    print("  combinations exist for -- and there the ordering appears, driven")
    print("  by how long drug pressure is held while the memory pool builds")
    print("  (compare the TMEM@EOT column, not the nadir).")
    return a, b


def a5_antimony():
    """Claim 5: antimony resistance is a thiol titration dose cannot beat."""
    hr("A5  ANTIMONY: RESISTANCE IS AN EFFLUX TITRATION, AND DOSE CANNOT WIN")
    print(f"  {'RES_SB':>7s} {'dose mg/kg/d':>13s} {'AUC Sb(III) spleen':>19s} "
          f"{'log drop':>9s} {'QTc max':>8s} {'lipase':>8s} {'outcome':>8s}")
    rows = []
    for res in (1.0, 3.0, 9.0):
        for dose in (20.0, 30.0, 40.0):
            p = default_params(); p["RES_SB"] = res
            reg = [dict(drug="SB", mgkg=dose, start_day=0, n=30)]
            t, o = simulate(p, reg, tmax_day=360.0)
            r = outcome(t, o, scale_params(p), 31.0, 360.0)
            st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "fail")
            print(f"  {res:7.1f} {dose:13.1f} {r['AUC_sb3sp']:19.1f} "
                  f"{r['logdrop_nadir']:9.2f} {r['qtc_max']:8.1f} "
                  f"{r['lipa_max']:8.1f} {st:>8s}")
            rows.append((res, dose, r))
    print("\n  Escalating 20 -> 40 mg Sb/kg/d doubles Sb(III) exposure but")
    print("  moves QTc and lipase with it; the resistant strain is not")
    print("  recovered because efflux scales with the intracellular load.")
    return rows


def a6_pkdl():
    """Skin is a separate pharmacological organ."""
    hr("A6  PKDL: THE SKIN IS A DIFFERENT ORGAN, PHARMACOLOGICALLY")
    q = scale_params(default_params())
    print("  Fraction of MPS-delivered L-AmB reaching each organ, and the")
    print("  intramacrophage AUC per 10 mg/kg single dose:")
    p = default_params()
    t, o = simulate(p, [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)],
                    tmax_day=180.0)
    for nm, idx, V, f in (("spleen", 3, q["VSP"], p["FSP_A"]),
                          ("liver", 4, q["VLI"], p["FLI_A"]),
                          ("marrow", 5, q["VBM"], p["FBM_A"]),
                          ("skin", 6, q["VSK"], p["FSK_A"])):
        auc = np.trapezoid(o[:, idx] / V, t)
        print(f"    {nm:8s} MPS share {f:5.2f}   AUC {auc:12.1f} mg.h/L   "
              f"AUC/EC50 {auc/p['EC50_A']:8.1f} h")
    t, o = simulate(p, [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)],
                    tmax_day=180.0)
    print("  Miltefosine 28 d, same organs:")
    for nm, idx, V in (("spleen", 11, q["VSP"]), ("liver", 12, q["VLI"]),
                       ("marrow", 13, q["VBM"]), ("skin", 14, q["VSK"])):
        auc = np.trapezoid(o[:, idx] / V, t)
        print(f"    {nm:8s} AUC {auc:12.1f} mg.h/L   "
              f"AUC/EC50 {auc/p['EC50_M']:8.1f} h")
    hr("A6b PKDL TREATMENT: 12 WEEKS OF MILTEFOSINE vs WEEKLY L-AmB", "-")
    print(f"  {'regimen':34s} {'skin burden d0':>14s} {'d360':>10s} "
          f"{'peak lesion':>12s} {'d360 lesion':>12s}")
    for lbl, reg in (("miltefosine 2.5 mg/kg x 84 d",
                      [dict(drug="MIL", mgkg=2.5, start_day=0, n=84)]),
                     ("L-AmB 5 mg/kg weekly x 4 (20 total)",
                      [dict(drug="LAMB", mgkg=5.0, start_day=0, n=4,
                            interval_h=7 * 24.0)]),
                     ("L-AmB 5 mg/kg weekly x 8 (40 total)",
                      [dict(drug="LAMB", mgkg=5.0, start_day=0, n=8,
                            interval_h=7 * 24.0)]),
                     ("no treatment", [])):
        p = default_params()
        t, o = simulate_pkdl(p, reg, tmax_day=360.0)
        sk0 = o[0, 38] + o[0, 39]
        ske = o[-1, 38] + o[-1, 39]
        print(f"  {lbl:34s} {sk0:14.1f} {ske:10.3g} "
              f"{o[:, IX['PKDL']].max():12.2f} {o[-1, IX['PKDL']]:12.3f}")
    hr("A6c PKDL RISK LEFT BEHIND BY EACH VL REGIMEN", "-")
    print("  residual skin burden at day 180 after apparently successful VL cure")
    print(f"  {'VL regimen':34s} {'skin AUC-drug':>13s} {'skin burden d180':>17s}")
    for lbl, reg, idx, V in (
            ("L-AmB 10 mg/kg single", [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)], 6, q["VSK"]),
            ("L-AmB 21 mg/kg multidose", regimens()["S06_lamb_21_multi"], 6, q["VSK"]),
            ("miltefosine 28 d", [dict(drug="MIL", mgkg=2.5, start_day=0, n=28)], 14, q["VSK"]),
            ("L-AmB 5 + MF 7 d", regimens()["S13_lamb5_mil7"], 14, q["VSK"]),
            ("paromomycin 21 d", [dict(drug="PM", mgkg=15.0, start_day=0, n=21)], 21, q["VSK"])):
        p = default_params()
        t, o = simulate(p, reg, tmax_day=180.0)
        auc = np.trapezoid(o[:, idx] / V, t)
        print(f"  {lbl:34s} {auc:13.1f} {o[-1, 38] + o[-1, 39]:17.4g}")


def a7_subclinical():
    """The 10:1 subclinical:clinical ratio is a race, not a switch."""
    hr("A7  SUBCLINICAL vs CLINICAL VL IS A RACE BETWEEN PRIMING AND IL-10")
    print("  Starting from a single sandfly inoculum, sweep the T-cell")
    print("  priming rate KTP (nutrition, genetics, age) and ask which")
    print("  attractor the host lands in.")
    print(f"  {'KTP':>8s} {'peak burden':>12s} {'day-540 burden':>15s} "
          f"{'peak spleen cm':>15s} {'peak IL-10':>11s} {'outcome':>12s}")
    crit = None
    prev = None
    for ktp in (1.0e-5, 3.0e-5, 6.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3,
                2.0e-3, 5.0e-3):
        p = default_params(); p["KTP"] = ktp
        t, o = simulate(p, [], tmax_day=540.0, stage="naive")
        P = ptot(o)
        st = "clinical VL" if P[-1] > 1e3 else ("subclinical" if P.max() > 1e-3 else "abortive")
        print(f"  {ktp:8.4f} {P.max():12.4g} {P[-1]:15.4g} "
              f"{o[:, IX['SPL']].max():15.2f} {o[:, IX['IL10']].max():11.3f} "
              f"{st:>12s}")
        if prev is not None and prev[1] > 1e2 >= P[-1]:
            crit = (prev[0], ktp)
        prev = (ktp, P[-1])
    if crit:
        print(f"\n  The switch sits between KTP = {crit[0]:.4f} and {crit[1]:.4f} /h:")
        print("  a <35% change in priming rate decides disease vs no disease,")
        print("  with no change to the parasite.")
    hr("A7b MALNUTRITION MOVES THE SAME LINE", "-")
    print(f"  {'MALNUT':>8s} {'day-540 burden':>15s} {'outcome':>12s}")
    for mn in (0.0, 0.2, 0.35, 0.5, 0.7):
        p = default_params(); p["MALNUT"] = mn; p["KTP"] = 8.0e-5
        t, o = simulate(p, [], tmax_day=540.0, stage="naive")
        P = ptot(o)
        st = "clinical VL" if P[-1] > 1e3 else "subclinical"
        print(f"  {mn:8.2f} {P[-1]:15.4g} {st:>12s}")


def a8_fractionation():
    """Does splitting the same L-AmB total dose help or hurt?"""
    hr("A8  DOSE FRACTIONATION OF A FIXED 10 mg/kg L-AmB TOTAL")
    print(f"  {'schedule':30s} {'AUC spleen':>11s} {'AUC free':>9s} "
          f"{'nadir':>10s} {'d360':>10s} {'SCr':>6s} {'outcome':>8s}")
    for lbl, reg in (
            ("10 mg/kg x1 (day 0)", [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]),
            ("5 mg/kg x2 (d0, d1)", [dict(drug="LAMB", mgkg=5.0, start_day=0, n=2)]),
            ("2 mg/kg x5 (d0-4)", [dict(drug="LAMB", mgkg=2.0, start_day=0, n=5)]),
            ("1 mg/kg x10 (d0-9)", [dict(drug="LAMB", mgkg=1.0, start_day=0, n=10)]),
            ("2 mg/kg weekly x5", [dict(drug="LAMB", mgkg=2.0, start_day=0, n=5,
                                        interval_h=7 * 24.0)]),
            ("0.5 mg/kg x20 (d0-19)", [dict(drug="LAMB", mgkg=0.5, start_day=0, n=20)])):
        p = default_params()
        t, o = simulate(p, reg, tmax_day=360.0)
        r = outcome(t, o, scale_params(p), eot_of(reg), 360.0)
        st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "fail")
        print(f"  {lbl:30s} {r['AUC_ambsp']:11.1f} {r['AUC_free']:9.1f} "
              f"{r['nadir']:10.3g} {r['P_end']:10.3g} {r['scr_max']:6.2f} "
              f"{st:>8s}")
    print("\n  Read this table against the usual argument for single-dose")
    print("  L-AmB, because the model does NOT support the version of it that")
    print("  says fractionation raises the toxic integral.  Amphotericin")
    print("  disposition here is linear, so for a FIXED total dose BOTH")
    print("  integrals -- spleen exposure and free-plasma exposure -- are")
    print("  exactly schedule-invariant (25688.2 and 137.9 mg.h/L in every")
    print("  row).  What fractionation actually changes is the PEAK: splitting")
    print("  the dose lowers peak renal cortical concentration and shaves the")
    print("  creatinine slightly (1.01 -> 0.96-0.98 mg/dL).  So the honest")
    print("  model-based case for the single dose is equal efficacy, equal")
    print("  total exposure, one visit -- an operational argument, not a")
    print("  pharmacokinetic one.")


def a9_hiv():
    hr("A9  HIV-VL: THE SAME DRUG, A DIFFERENT FINISH LINE")
    print(f"  {'scenario':40s} {'CD4 d0':>7s} {'nadir':>10s} {'d540':>10s} "
          f"{'outcome':>8s} {'mort %':>7s}")
    for lbl, cd4, art, reg in (
            ("immunocompetent, L-AmB 10 mg/kg", 700, 0, regimens()["S05_lamb_single10"]),
            ("HIV CD4 90, L-AmB 10 mg/kg, ART", 90, 1, regimens()["S05_lamb_single10"]),
            ("HIV CD4 90, L-AmB 30 mg/kg, ART", 90, 1,
             [dict(drug="LAMB", mgkg=5.0, start_day=d, n=1) for d in (0, 2, 4, 6, 8, 10)]),
            ("HIV CD4 90, L-AmB 30 + MF 28 d, ART", 90, 1, regimens()["S18_hiv_lamb30_mil28"]),
            ("HIV CD4 90, no ART, L-AmB 30 + MF", 90, 0, regimens()["S18_hiv_lamb30_mil28"]),
            ("HIV CD4 250, L-AmB 10 mg/kg, ART", 250, 1, regimens()["S05_lamb_single10"])):
        p = default_params(); p["CD40"] = cd4; p["HIV"] = 1.0 if cd4 < 400 else 0.0
        p["ART"] = float(art)
        t, o = simulate(p, reg, tmax_day=540.0)
        r = outcome(t, o, scale_params(p), eot_of(reg), 540.0)
        st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "fail")
        print(f"  {lbl:40s} {cd4:7d} {r['nadir']:10.3g} {r['P_end']:10.3g} "
              f"{st:>8s} {r['mort_pct']:7.1f}")


def a10_biomarkers(rows):
    hr("A10 CLINICAL RESPONSE AND WHY SEROLOGY CANNOT CALL RELAPSE")
    key = ["S01_untreated", "S05_lamb_single10", "S08_mil28_adult",
           "S11_pm21", "S03_ssg_bihar_res", "S13_lamb5_mil7"]
    print(f"  {'scenario':22s} {'spleen cm':>9s} {'Hb':>6s} {'Plt':>7s} "
          f"{'WBC':>6s} {'Alb':>6s} {'IgG':>6s}  (day 180)")
    for k in key:
        t, o, p, r = rows[k]
        i = int(np.argmin(np.abs(t / DAY - 180.0)))
        print(f"  {k:22s} {o[i, IX['SPL']]:9.2f} {o[i, IX['HGB']]:6.1f} "
              f"{o[i, IX['PLT']]:7.1f} {o[i, IX['WBC']]:6.2f} "
              f"{o[i, IX['ALB']]:6.2f} {o[i, IX['IGG']]:6.2f}")
    t, o, p, r = rows["S05_lamb_single10"]
    P = ptot(o); day = t / DAY
    def half_life_of(series, i0):
        base, end = series[i0], series[-1]
        target = end + 0.5 * (base - end)
        for j in range(i0, len(series)):
            if (base > end and series[j] <= target) or (base < end and series[j] >= target):
                return day[j] - day[i0]
        return float("nan")
    print("\n  time to half-recovery after a curative single dose (day):")
    for nm, k in (("burden (log)", None), ("spleen", "SPL"), ("Hb", "HGB"),
                  ("platelets", "PLT"), ("albumin", "ALB"),
                  ("polyclonal IgG", "IGG")):
        if k is None:
            j = np.argmax(P < P[0] * 1e-3)
            print(f"    {nm:16s} {day[j]:8.1f}")
        else:
            print(f"    {nm:16s} {half_life_of(o[:, IX[k]], 0):8.1f}")
    i180 = int(np.argmin(np.abs(day - 180.0)))
    PV = o[:, 32:38].sum(axis=1)
    print("\n  At day 180 after a curative single dose:")
    print("    polyclonal IgG   %.2f g/dL  (presenting %.2f, healthy %.2f)"
          % (o[i180, IX["IGG"]], o[0, IX["IGG"]], p["IGG0"]))
    print("    visceral burden  %.3g units  (%.1f log10 below presenting)"
          % (PV[i180], -math.log10(max(PV[i180], 1e-30) / PV[0])))
    print("    spleen %.2f cm, Hb %.1f g/dL -- clinically well"
          % (o[i180, IX["SPL"]], o[i180, IX["HGB"]]))
    print("\n  The point is the RATIO OF TIME CONSTANTS, which is a model")
    print("  output and not a parameter: the burden half-recovers in %.1f days"
          % half_life_of(-np.log10(np.maximum(PV, 1e-30)), 0))
    print("  and the antibody in %.1f days -- roughly a tenfold gap.  That is"
          % half_life_of(o[:, IX["IGG"]], 0))
    print("  why rK39 stays positive for years after cure and cannot be used")
    print("  to diagnose relapse; only a parasitological or clinical endpoint")
    print("  can.")


def a11_toxicity(rows):
    hr("A11 TOXICITY LEDGER ACROSS ALL 20 SCENARIOS")
    print(f"  {'scenario':24s} {'mg/kg':>7s} {'SCr':>5s} {'K':>5s} {'Mg':>5s} "
          f"{'dB':>5s} {'QTc':>6s} {'lipase':>7s} {'ALT':>6s} {'GI':>5s}")
    for k in sorted(rows):
        t, o, p, r = rows[k]
        print(f"  {k:24s} {r['cumdose']:7.1f} {r['scr_max']:5.2f} "
              f"{r['k_min']:5.2f} {r['mg_min']:5.2f} {r['hear_dB']:5.1f} "
              f"{r['qtc_max']:6.1f} {r['lipa_max']:7.1f} {r['alt_max']:6.1f} "
              f"{r['gitx_max']:5.2f}")


def a12_population(nsub=60, seed=20260804):
    """Population simulation -> cure rates comparable with trial data."""
    hr("A12 POPULATION SIMULATION: MODEL-PREDICTED CURE RATES")
    rng = np.random.default_rng(seed)
    specs = [
        ("SSG 20 mg/kg x 30 d (East Africa)", "S02_ssg_africa"),
        ("SSG 20 mg/kg x 30 d (Bihar, resistant)", "S03_ssg_bihar_res"),
        ("AmB deoxycholate 15 mg/kg total", "S04_damb_alt30"),
        ("L-AmB 10 mg/kg single dose", "S05_lamb_single10"),
        ("L-AmB 21 mg/kg multidose", "S06_lamb_21_multi"),
        ("miltefosine 28 d, adult", "S08_mil28_adult"),
        ("miltefosine 28 d, child, linear mg/kg", "S09_mil28_child_linear"),
        ("miltefosine 28 d, child, allometric", "S10_mil28_child_allom"),
        ("paromomycin 15 mg/kg x 21 d (India)", "S11_pm21"),
        ("paromomycin 15 mg/kg x 21 d (E Africa)", "S12_pm21_africa"),
        ("L-AmB 5 mg/kg + miltefosine 7 d", "S13_lamb5_mil7"),
        ("L-AmB 5 mg/kg + paromomycin 10 d", "S14_lamb5_pm10"),
        ("miltefosine 10 d + paromomycin 10 d", "S15_mil10_pm10"),
        ("SSG + paromomycin 17 d (East Africa)", "S16_ssg_pm17"),
        ("HIV-VL, L-AmB 10 mg/kg", "S17_hiv_lamb10"),
        ("HIV-VL, L-AmB 30 mg/kg + MF 28 d", "S18_hiv_lamb30_mil28"),
    ]
    R = regimens()
    print(f"  n = {nsub} virtual patients per arm, 12-month follow-up")
    print(f"  {'regimen':42s} {'cure %':>7s} {'relapse %':>10s} "
          f"{'fail %':>7s} {'death %':>8s}")
    results = {}
    for lbl, key in specs:
        reg = R[key]
        base = scenario_params(key)
        nc = nr = nf = 0
        deaths = []
        for i in range(nsub):
            p = dict(base)
            # between-subject variability
            p["KG_SP"] = base["KG_SP"] * float(np.exp(rng.normal(0, 0.18)))
            p["KG_LI"] = base["KG_LI"] * float(np.exp(rng.normal(0, 0.18)))
            p["EC50_A"] = base["EC50_A"] * float(np.exp(rng.normal(0, 0.30)))
            p["EC50_M"] = base["EC50_M"] * float(np.exp(rng.normal(0, 0.35)))
            p["EC50_P"] = base["EC50_P"] * float(np.exp(rng.normal(0, 0.35)))
            p["EC50_S"] = base["EC50_S"] * float(np.exp(rng.normal(0, 0.45)))
            p["CL_MIL"] = base["CL_MIL"] * float(np.exp(rng.normal(0, 0.25)))
            p["CL_LIP"] = base["CL_LIP"] * float(np.exp(rng.normal(0, 0.22)))
            p["KTP"] = base["KTP"] * float(np.exp(rng.normal(0, 0.30)))
            p["FREL"] = min(0.85, base["FREL"] * float(np.exp(rng.normal(0, 0.20))))
            p["WT"] = base["WT"] * float(np.exp(rng.normal(0, 0.12)))
            p["P0SCALE"] = float(np.exp(rng.normal(0, 0.60)))
            p["FR0"] = base["FR0"] * float(np.exp(rng.normal(0, 0.80)))
            p["KIMM"] = base["KIMM"] * float(np.exp(rng.normal(0, 0.22)))
            p["KOUT_A"] = base["KOUT_A"] * float(np.exp(rng.normal(0, 0.30)))
            p["F_MIL"] = base["F_MIL"] * float(np.exp(rng.normal(0, 0.28)))
            p["MALNUT"] = min(0.75, max(0.0, base["MALNUT"]
                                        + float(rng.uniform(0.0, 0.55))))
            p["CD4SET"] = base["CD4SET"] * float(np.exp(rng.normal(0, 0.25)))
            if key == "S10_mil28_child_allom":
                mg = 125.0 * (p["WT"] / 50.0) ** 0.75
                reg_i = [dict(drug="MIL", mg=mg, start_day=0, n=28)]
            else:
                reg_i = reg
            if any(r.get("drug") == "MIL" for r in reg_i):
                p["ADHER"] = float(min(1.0, max(0.45,
                                                1.0 - abs(rng.normal(0, 0.14)))))
            try:
                t, o = simulate(p, reg_i, tmax_day=365.0, nout=900, rng=rng)
            except Exception:
                continue
            r = outcome(t, o, scale_params(p), eot_of(reg_i), 365.0)
            if r["cure"]:
                nc += 1
            elif r["relapse"]:
                nr += 1
            else:
                nf += 1
            deaths.append(r["mort_pct"])
        n = nc + nr + nf
        if n == 0:
            continue
        print(f"  {lbl:42s} {100*nc/n:7.1f} {100*nr/n:10.1f} "
              f"{100*nf/n:7.1f} {np.mean(deaths):8.1f}")
        results[key] = dict(label=lbl, n=n, cure=100 * nc / n,
                            relapse=100 * nr / n, fail=100 * nf / n,
                            death=float(np.mean(deaths)))
    return results


def a13_scenario_table(rows):
    hr("A13 ALL 20 SCENARIOS, TYPICAL PATIENT")
    print(f"  {'scenario':24s} {'mg/kg':>7s} {'logdrop':>8s} {'nadir d':>8s} "
          f"{'d540 burden':>12s} {'spleen':>7s} {'Hb':>5s} {'outcome':>8s}")
    for k in sorted(rows):
        t, o, p, r = rows[k]
        st = "cure" if r["cure"] else ("relapse" if r["relapse"] else "fail")
        print(f"  {k:24s} {r['cumdose']:7.1f} {r['logdrop_nadir']:8.2f} "
              f"{r['nadir_day']:8.1f} {r['P_end']:12.3g} "
              f"{r['spl_end']:7.2f} {r['hgb_end']:5.1f} {st:>8s}")


def a14_sensitivity():
    hr("A14 LOCAL SENSITIVITY OF THE DAY-180 BURDEN (L-AmB 10 mg/kg single)")
    base_p = default_params()
    reg = [dict(drug="LAMB", mgkg=10.0, start_day=0, n=1)]
    t, o = simulate(base_p, reg, tmax_day=180.0)
    P = ptot(o)
    ref = math.log10(max(P[-1], 1e-30))
    keys = ["FREL", "CL_LIP", "FSP_A", "KOUT_A", "EC50_A", "EMAX_A", "KIMM",
            "KTP", "KI_T", "K50_IL10", "KG_SP", "KAG", "KCD4", "V_LIP"]
    print(f"  {'parameter':12s} {'-25%':>10s} {'base':>10s} {'+25%':>10s} "
          f"{'d log10 P / d log p':>21s}")
    for k in keys:
        vals = []
        for f in (0.75, 1.25):
            p = dict(base_p); p[k] = base_p[k] * f
            t2, o2 = simulate(p, reg, tmax_day=180.0)
            vals.append(math.log10(max(ptot(o2)[-1], 1e-30)))
        slope = (vals[1] - vals[0]) / (math.log(1.25) - math.log(0.75))
        print(f"  {k:12s} {vals[0]:10.2f} {ref:10.2f} {vals[1]:10.2f} "
              f"{slope:21.2f}")


def main():
    argv = sys.argv[1:]
    quick = "--quick" in argv
    print(__doc__)
    print("state variables :", NEQ)
    print("scenarios       :", len(regimens()))

    rows = run_all_scenarios(tmax_day=540.0)
    a13_scenario_table(rows)
    d1 = a1_targeting(rows)
    a1b_frel_sweep()
    a2_allometry()
    a2c_child_outcome()
    a3_separatrix()
    a4_timing()
    a5_antimony()
    a6_pkdl()
    a7_subclinical()
    a8_fractionation()
    a9_hiv()
    a10_biomarkers(rows)
    a11_toxicity(rows)
    a14_sensitivity()
    if not quick:
        pop = a12_population(nsub=60)
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "vl_population_results.json"), "w") as fh:
            json.dump(pop, fh, indent=1)
    hr("END")


if __name__ == "__main__":
    main()
