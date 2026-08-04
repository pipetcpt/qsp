#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pws_reference_model.py — Prader-Willi Syndrome (PWS) QSP model
================================================================
Dependency-free (pure standard-library) RK4 reference implementation of the
same 64-state ODE system written in `pws_mrgsolve_model.R`.  Its purpose is to
make every number quoted in README.md reproducible with `python3` alone, and to
serve as the arithmetic check on the mrgsolve translation.

RUN
---
    python3 pws_reference_model.py            > pws_reference_output.txt
    python3 pws_reference_model.py --fast     # coarser dt, for quick edits

THE FIVE STRUCTURAL CHOICES THIS FILE MAKES
-------------------------------------------
(1) ONE ENZYME, MANY SUBSTRATES.  The paternal 15q11-q13 lesion is not written
    as a list of hormone deficits.  It is written as one scalar — PC1/3
    (PCSK1) activity, `PC13` — feeding five prohormone branches that each have
    their own *escape ratio*

        eps_i = d_i * Km_i / kcat_i        (precursor escape / processing)

    In the sub-saturating regime the steady-state product flux of branch i is
    flux_i/S_i = 1/(1 + eps_i/PC13), so

        L_i  = 1 - (1 + eps_i)/(1 + eps_i/PC13)          [Eq. A]  product loss
        R_i  = (1 + eps_i)/(PC13 + eps_i)                [Eq. B]  precursor gain

    L_i rises with eps_i while R_i FALLS with eps_i: the branches that lose
    product are exactly the branches that do NOT accumulate precursor.  And
    because flux + escape = S at steady state, any assay that cross-reacts
    with precursor and product measures S and is EXACTLY blind to PC13.

(2) SATIETY IS A HARMONIC MEAN AND THE MELANOCORTIN ARM IS IN SERIES.
    SAT = 1/sum_i(w_i/x_i) — the weakest arm rate-limits.  alpha-MSH is not
    its own arm: it is a saturating INPUT GAIN on the PVN-oxytocin arm,
    because MC4R satiety signalling is relayed through PVN oxytocin neurons.
    PWS blocks pro-oxytocin processing, i.e. DOWNSTREAM of MC4R, so MC4R
    agonism is bounded above by the very quantity PWS lacks.  Oxytocin
    receptor agonism, by contrast, acts on BOTH sides of the synapse — it
    raises the satiety arm AND presynaptically inhibits AgRP.  That asymmetry
    is the whole difference between the carbetocin and setmelanotide results,
    and it needs no drug-specific fitting.

(3) FOOD-SEEKING IS BISTABLE.  `SEEK` carries a Hill-4 self-reinforcement
    term, so at food-secure inputs it has two stable states separated by a
    separatrix, and at free-access inputs the LOW state is annihilated
    (saddle-node).  Food security does not lower the drive: it restores the
    existence of a low branch.  Drugs then decide which branch you sit on.

(4) TWO CLOCKS ON THE AIRWAY.  Growth hormone reaches the upper airway twice:
    fast (IGF-1 -> lymphoid/adenotonsillar hypertrophy, tau ~ 20 d, with a
    100-d adaptation that involutes it) and slow (lean mass -> respiratory
    muscle strength, and falling fat mass; tau 75-200 d).  The transient
    worsening of AHI in the first weeks of GH is an emergent consequence of
    tau_fast < tau_slow, not a fitted curve.

(5) INTAKE IS ANCHORED TO WHAT IS PROVIDED, NOT TO WHAT IS SPENT.  Caregivers
    serve age-normative portions, so EI is referenced to EIREQ(age) of a
    reference child.  A control is therefore in energy balance by
    construction, while PWS — whose expenditure is low because its lean mass
    and activity are low — is in surplus on a NORMAL intake.  This is what
    makes weight gain precede hyperphagia (Miller phase 2a) without any
    appetite parameter changing.

TIME SCALE
----------
The model runs on a DAY scale.  Every species whose true plasma half-life is
sub-hour (ghrelin, insulin, GH, LH) is a day-scale pool: its fast kinetics are
pre-averaged, which is why no rate constant exceeds ~14 /d.  This is legitimate
for endpoints (fat mass, HQ-CT, IGF-1 SDS, AHI, HbA1c) that integrate over
weeks, and it is what keeps RK4 stable at dt = 0.125 d over a 30-year run.
Drug concentrations are consequently APPARENT day-scale averages calibrated to
reproduce the correct 24-h AUC, not the true Cmax.

UNITS
-----
time d | mass kg | height cm | energy kcal/d | AG,UAG pg/mL | PYY,GLP-1 pmol/L
leptin ng/mL | GH ug/L | IGF-1 ng/mL | glucose mg/dL | insulin uU/mL | HbA1c %
adiponectin ug/mL | testosterone ng/dL | AHI events/h | Cobb deg.  Neuropeptide
tones and prohormone pools are dimensionless, normalized so the non-PWS
reference sits at 1.0.
"""

import math
import sys

# ----------------------------------------------------------------------------
# 0.  STATE VECTOR
# ----------------------------------------------------------------------------
STATES = [
    # --- prohormone pools (PC1/3 substrates) -------------------------------
    "PPOMC", "PPOXT", "PPGHRH", "PPINS", "PPGHR",            # 0-4
    # --- neuropeptide products / gut-brain signals -------------------------
    "AMSH", "OXT", "AGRP", "AG", "UAG", "PYY", "GLP1E",      # 5-11
    # --- adiposity, energy, growth -----------------------------------------
    "LEP", "FM", "LFM", "GUT", "HT", "ACT", "BA",            # 12-18
    # --- food-seeking / behaviour ------------------------------------------
    "SEEK", "HQ", "BEH",                                     # 19-21
    # --- somatotropic axis -------------------------------------------------
    "GHRH", "SST", "GH", "IGF1", "IGFBP3",                   # 22-26
    # --- glucose / insulin -------------------------------------------------
    "GLU", "INS", "HBA1C", "ADPN",                           # 27-30
    # --- gonadal axis / bone -----------------------------------------------
    "GNRH", "LH", "SEX", "BMD",                              # 31-34
    # --- upper airway / respiration ---------------------------------------
    "LYMPH", "ADAPT", "RMS", "AHI",                          # 35-38
    # --- musculoskeletal / thermoregulation --------------------------------
    "TONE", "COBB", "TEMP",                                  # 39-41
    # --- drug PK -----------------------------------------------------------
    "AGHD", "AGHC",      # somatropin daily SC        42-43
    "ALGD", "ALGC",      # long-acting weekly GH      44-45
    "ACBD", "ACBC",      # carbetocin intranasal      46-47
    "ADZD", "ADZC",      # diazoxide choline ER       48-49
    "ASMD", "ASMC",      # setmelanotide SC           50-51
    "ASGD", "ASGC",      # semaglutide SC             52-53
    "AOCD", "AOCC",      # octreotide LAR             54-55
    "ALVD", "ALVC",      # livoletide SC              56-57
    "ATSD", "ATSC",      # testosterone IM            58-59
    "AMFD", "AMFC",      # metformin PO               60-61
    # --- diagnostics -------------------------------------------------------
    "CUMEI", "CUMTEE",                                       # 62-63
]
IX = {n: i for i, n in enumerate(STATES)}
NST = len(STATES)

# ----------------------------------------------------------------------------
# 1.  PARAMETERS
# ----------------------------------------------------------------------------
P = dict(

    # ==== 1.1  genetic lesion layer ======================================
    LES=1.0,        # 1 = PWS (paternal 15q11-q13 expression absent), 0 = control
    DPC13=0.60,     # fractional loss of hypothalamic/islet PC1/3 activity
    FOXTN=0.75,     # surviving PVN oxytocin neuron fraction (Swaab 1995 -42%)
    FPYYS=0.60,     # blunted post-prandial PYY response
    FVAG=0.70,      # reduced vagal afferent / gastric distension gain
    FGHRC=2.20,     # increased gastric X/A (ghrelin) cell secretory drive
    FADPN=1.45,     # relative hyperadiponectinaemia for a given fat mass
    FHYPO=0.22,     # hypogonadotropic GnRH amplitude factor
    FSOM=0.75,      # somatotroph secretory reserve factor
    DMKRN3=1.20,    # y — MKRN3 loss advances the pubertal gate
    FTONE0=0.60,    # neonatal hypotonia factor
    FSUB=1.00,      # subtype modifier on behaviour (DEL 1.15 / mUPD 0.90)

    # ==== 1.2  prohormone processing (Eq. A / Eq. B) ======================
    # KC_i normalized to 1/d so that eps_i == KD_i (the escape ratio)
    KCPOMC=1.0, EPSPOMC=1.50,      # POMC       -> alpha-MSH
    KCPOXT=1.0, EPSPOXT=6.00,      # pro-OXT    -> oxytocin      (high escape)
    KCPGHRH=1.0, EPSPGHRH=1.00,    # pro-GHRH   -> GHRH
    KCPINS=1.0, EPSPINS=0.15,      # proinsulin -> insulin       (low escape)
    KCPGHR=1.0, EPSPGHR=0.10,      # proghrelin -> ghrelin       (low escape)
    SPOMC=1.0, SPOXT=1.0, SPGHRH=1.0, SPINS=1.0, SPGHR=1.0,
    KLEPP=0.45,     # leptin drive on POMC transcription
    KLEPR0=12.0,    # ng/mL — CNS leptin signal half-saturation

    # ==== 1.3  neuropeptide tones (normalized: control = 1.0) =============
    KSAMSH=6.25, KDAMSH=2.50,      # KSAMSH = KDAMSH*(1+EPSPOMC)
    KSOXT=49.0, KDOXT=7.00,        # KSOXT  = KDOXT *(1+EPSPOXT)
    SAGRP=3.0, KDAGRP=1.0,         # AgRP/NPY orexigenic tone
    KINHM=0.00,                    # alpha-MSH does NOT inhibit AgRP neurons:
                                   # melanocortin action on satiety is POST-
                                   # synaptic, at the PVN, not at the arcuate
    KINHO=2.00,                    # oxytocin inhibition of AgRP (presynaptic)
    KLEPA=25.0,                    # ng/mL — leptin suppression of AgRP
    KFASTA=0.40,                   # fasting amplification of AgRP
    X1REF=0.285,                   # reference PWS oxytocin-arm value (BEH term)
    X1CTL=1.000,                   # control oxytocin-arm value (relay integrity)

    # ==== 1.4  ghrelin ====================================================
    KSAG=12400.0,   # pg/mL/d per unit acyl processing flux
    FACYL=0.22,     # GOAT acylation fraction
    FACYLP=1.35,    # raised acylated:unacylated ghrelin ratio in PWS
    KCLAG=4.0, KDEACYL=0.55,
    KSUAG=10500.0, KCLUAG=4.0,
    KIINS=18.0,     # uU/mL — insulin suppression of ghrelin secretion
    KFASTG=0.45,
    KAG50=300.0,    # pg/mL — GHS-R1a occupancy half-saturation (SATURATION!)
    KGHRD=0.10,     # weight of the ghrelin arm in orexigenic drive

    # ==== 1.5  gut peptides / gastric handling ============================
    SPYY=18.0, KDPYY=2.0, PYY0=9.0,
    SGLP1=26.0, KDGLP1=2.6, GLP10=10.0,
    KEMPT=5.0,      # /d — day-scale gastric emptying
    FGE=0.70,       # PWS delayed gastric emptying

    # ==== 1.6  leptin =====================================================
    KSLEP=1.37, KDLEP=1.6,

    # ==== 1.7  energy balance / body composition ==========================
    # REE = REEref(age) * (WREEL*LFM/LFMnorm + WREEF*FM/FMnorm)
    WREEL=0.88, WREEF=0.12,
    PALX=0.62,      # activity thermogenesis as fraction of REE at ACT = 1
    ACTREF=1.000,   # reference activity level (defines EIREQ)
    FDIT=0.10,      # diet-induced thermogenesis
    KGHRMR=0.16,    # GH signal -> resting energy expenditure
    ACT0=1.00, TAUACT=30.0, KFATACT=0.45, KFA=22.0,
    RHOF=9440.0, RHOL=1800.0, FMMIN=1.2,
    KSEXPART=0.35,   # gonadal-steroid contribution to the lean-mass target
    EIBF=0.933,     # appetite set-point as fraction of the requirement
    KTITR=1.20,     # caregiver titration gain toward the weight target
    KEIS=0.55,      # realized food-seeking -> intake
    KEIH=0.35,      # hyperphagia score -> intake
    KLEPEI=0.25,    # adiposity negative feedback on intake
    EICAPF=1.40,    # absolute physiological ceiling on EI / EIREQ
    TAULFM=180.0,   # d — lean mass tracks its developmental target
    WLIGF0=0.78,    # lean target: floor of the IGF-1 term (control = 1.0)
    LIGFCAP=1.25,   # lean target: IGF-1 drive is capped here
    WLTON0=0.80,    # lean target: floor of the muscle-tone term
    FMODCAP=1.05,   # lean target: overall ceiling
    KLHT=0.35,      # kg lean mass per cm of height gain (diagnostic only)
    GVFL=0.72,      # GH-independent floor of growth velocity
    GVCAP=1.45,     # ceiling on the IGF-1 drive to growth velocity
    KSEXH=380.0,    # ng/dL — sex-steroid half-constant for the growth spurt
    BACLOSE=17.2, KBASEX=0.30,

    # ==== 1.8  food-seeking bistability ===================================
    KSON=0.120, KSELF=3.10, KSHALF=0.450, NSELF=4.0, KSOFF=2.20,
    RFL=0.55,       # floor of the reinforcement factor on self-sensitization:
                    # cue-driven food-seeking is maintained partly by the food
                    # it actually obtains, so removing ACCESS (not just cues)
                    # is what can destroy the latched upper state

    W1=0.44, W2=0.16, W3=0.13, W4=0.14, W5=0.13,
    KREL=0.20,      # MC4R -> PVN relay half-saturation (the series gain)
    NAGRP=1.00,
    HQMAX=36.0, WHQS=0.35, WHQD=0.65, KHQ=2.75, NHQ=1.60, TAUHQ=14.0,
    AGEHP=6.50, WAGEHP=1.40,
    AGEHP4=32.0, WAGEHP4=6.0, FHP4=0.22,

    # ==== 1.9  behaviour ==================================================
    BEH0=34.0, TAUBEH=45.0, KBFR=0.78, KBOXT=0.16,

    # ==== 1.10 somatotropic axis ==========================================
    KSGHRH=4.0, KDGHRH=2.0,
    TAUSST=5.0, KSSTI=0.55, KSSTG=0.60,
    KCLGH=5.5,
    KIGFFB=260.0,   # ng/mL — IGF-1 feedback on GH secretion
    KFMGH=26.0,     # kg — fat mass suppression of GH secretion
    KGIGF=3.00,     # saturating GH -> IGF-1 half-constant (GHREL units)
    KDIGF=1.15,
    KBP3=2.9, KDBP3=1.05,
    FPOTGH=0.63,    # potency of a SMOOTH day-scale somatropin exposure
                    # relative to pulsatile pituitary GH.  Pulsatility is more
                    # IGF-1-efficient per unit AUC, so a flattened day-scale
                    # surrogate must carry a factor < 1 to land on the observed
                    # +1.5 to +2 SDS of IGF-1 at 0.035 mg/kg/d.

    # ==== 1.11 glucose / insulin ==========================================
    KEGP=2100.0, KIEGP=14.0, KGABS=0.055, KGU0=6.5, KGUI=0.93,
    VG=2.6, KGGLU=105.0, KSINS=88.0, VINS=1.0, KCLINS=6.0,
    SI0=1.00, KADPNSI=0.42, KFMSI=0.85, KFS=24.0, KGHSI=0.30,
    TAUA1C=35.0, KSADPN=17.8, KDADPN=1.1, KFADP=30.0,

    # ==== 1.12 gonadal axis / bone ========================================
    BAPUB=11.5, WPUB=1.10, SGNRH=1.0, TAUGN=8.0,
    KLH=4.6, KDLH=2.4, KSEX=360.0, KDSEX=1.15,
    KBMDU=0.00135, KBMDR=0.00120, KBMDL=0.35, KBMDI=0.22, KBMDH=0.55,

    # ==== 1.13 upper airway ===============================================
    KLYU=0.075, KLYD=0.052, LMAX=2.30, KADL=3.4, TAUAD=100.0,
    RMS0=1.00, TAURMS=75.0, KRMSIGF=0.28, NRMS=0.60,
    AHI0=1.6, KAHIL=7.0, KAHIF=11.0, FLR0=0.35, KAHIR=6.5,
    KAHIC=6.0, TAUAHI=7.0,

    # ==== 1.14 muscle tone / scoliosis / thermoregulation =================
    TAUTONE=60.0, KTIGF=0.15, TONEA=0.78, FTONEP=0.78,
    KCOB=0.058, NCOB=1.30, TAUTEMP=420.0,

    # ==== 1.15 drug PK (apparent, day-scale; AUC-calibrated) ==============
    KAGH=5.5, KEGH=1.4, VGH=2.60,          # somatropin  (V in L/kg)
    KALG=0.42, KELG=0.55, VLG=0.95,        # long-acting weekly GH (L/kg)
    KACB=14.0, KECB=8.0, VCB=32.0, FCB=0.055,   # carbetocin intranasal
    KADZ=1.5, KEDZ=0.594, VDZ=0.62,        # diazoxide choline ER (L/kg)
    KASM=3.0, KESM=1.45, VSM=52.0,         # setmelanotide
    KASG=0.30, KESG=0.099, VSG=12.5,       # semaglutide
    KAOC=0.075, KEOC=0.34, VOC=22.0,       # octreotide LAR
    KALV=6.0, KELV=2.8, VLV=0.20,          # livoletide (L/kg)
    KATS=0.35, KETS=0.90, VTS=1400.0,      # testosterone enanthate
    KAMF=6.0, KEMF=2.9, VMF=63.0,          # metformin

    # ==== 1.16 drug PD ====================================================
    ECBMAX=1.35, ECB50=0.55,      # carbetocin OXTR agonism        (ng/mL)
    EV1AMX=1.30, EV1A50=2.40,     # carbetocin V1a counter-effect  (ng/mL)
    EDZMAX=0.40, EDZ50=15.0,      # diazoxide -> AgRP K-ATP        (ug/mL)
    EDZIMX=0.70, EDZI50=22.0,     # diazoxide -> beta-cell K-ATP   (the AE)
    ESMMAX=3.20, ESM50=9.0,       # setmelanotide MC4R agonism     (ng/mL)
    ESGGE=1.20, ESG50=18.0,       # semaglutide gastric emptying   (ng/mL)
    ESGSAT=1.60, ESGEI=0.22, ESGSI=0.30,
    ESGAM=0.45,                   # GLP-1R agonist inhibition of AgRP neurons
    EOCG=3.10, EOC50=1.40,        # octreotide ghrelin suppression (ng/mL)
    EOCGH=2.60, EOCINS=0.80,
    KLVAG=0.008,                  # livoletide AG antagonism       (per ng/mL)
    EMFSI=0.22, EMF50=1.10, EMFEI=0.035,
)

# ----------------------------------------------------------------------------
# 2.  REFERENCE (non-PWS) DEVELOPMENTAL SCAFFOLD
# ----------------------------------------------------------------------------
_AGES = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
         13, 14, 15, 16, 17, 18, 20, 25, 30, 40, 60]
_LFMN = [2.7, 5.4, 7.2, 9.6, 11.4, 13.0, 14.5, 16.2, 18.0, 20.0, 22.2, 24.8,
         27.5, 30.5, 34.5, 39.0, 43.5, 47.0, 49.5, 51.0, 53.0, 55.0, 56.0,
         56.0, 54.0]
_FMN = [0.5, 2.3, 2.9, 3.2, 3.3, 3.4, 3.6, 4.0, 4.5, 5.2, 6.0, 6.8, 7.6,
        8.2, 8.5, 8.8, 9.2, 9.8, 10.5, 11.2, 12.5, 14.0, 15.5, 18.0, 20.0]
_HTN = [50, 68, 76, 87, 96, 103, 110, 116, 122, 128, 133, 138, 144, 150,
        157, 164, 170, 174, 176, 177, 177, 177, 177, 177, 176]
_IGFN = [50, 55, 60, 68, 76, 85, 95, 110, 124, 140, 160, 180, 215, 260,
         310, 350, 348, 340, 310, 280, 250, 210, 190, 160, 120]
_GHN = [3.6, 3.5, 3.4, 3.2, 3.1, 3.0, 3.0, 2.9, 2.9, 2.9, 3.0, 3.1, 3.4,
        3.6, 4.0, 4.2, 4.0, 3.8, 3.4, 3.0, 2.4, 1.8, 1.5, 1.1, 0.8]
_LYN = [0.55, 0.70, 0.85, 1.05, 1.20, 1.28, 1.30, 1.26, 1.18, 1.10, 1.02,
        0.96, 0.90, 0.86, 0.82, 0.79, 0.77, 0.75, 0.74, 0.73, 0.72, 0.70,
        0.70, 0.70, 0.70]
_CAPN = [1.00, 0.85, 0.62, 0.42, 0.30, 0.23, 0.18, 0.15, 0.12, 0.10, 0.09,
         0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05,
         0.05, 0.06, 0.08]
_SUCKP = [0.42, 0.50, 0.66, 0.86, 0.96, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
          1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
          1.00, 1.00, 1.00]
_GVB = [36, 20, 13, 10, 8.0, 7.2, 6.5, 6.2, 6.0, 5.6, 5.0, 4.6, 4.0, 3.4,
        2.6, 1.8, 1.0, 0.45, 0.15, 0.05, 0, 0, 0, 0, 0]
_GVS = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2, 0.8, 1.8, 2.6, 3.8, 4.6, 4.4, 3.2,
        1.6, 0.7, 0.2, 0, 0, 0, 0, 0]
# reference sex-steroid exposure (SEXrel of a non-PWS individual) by age
_SXR = [0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.05,
        0.10, 0.24, 0.46, 0.70, 0.86, 0.94, 0.98, 1.00, 1.00, 1.00, 1.00,
        1.00, 0.95, 0.80]


def _interp(tab, a):
    if a <= _AGES[0]:
        return tab[0]
    if a >= _AGES[-1]:
        return tab[-1]
    lo, hi = 0, len(_AGES) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if _AGES[mid] <= a:
            lo = mid
        else:
            hi = mid
    f = (a - _AGES[lo]) / (_AGES[hi] - _AGES[lo])
    return tab[lo] + f * (tab[hi] - tab[lo])


def LFMNORM(a):
    return _interp(_LFMN, a)


def FMNORM(a):
    return _interp(_FMN, a)


def HTNORM(a):
    return _interp(_HTN, a)


def IGFNORM(a):
    return _interp(_IGFN, a)


def GHNORM(a):
    return _interp(_GHN, a)


def LYNORM(a):
    return _interp(_LYN, a)


def CAPNORM(a):
    return _interp(_CAPN, a)


def SUCKPWS(a):
    return _interp(_SUCKP, a)


def GVBASE(a):
    return _interp(_GVB, a)


def GVSPURT(a):
    return _interp(_GVS, a)


def SEXRELREF(a):
    return _interp(_SXR, a)


def SCHOFIELD(a, bw):
    """Schofield resting energy expenditure, kcal/d, by age band and weight."""
    if a < 3.0:
        return max(90.0, 60.9 * bw - 54.0)
    if a < 10.0:
        return 22.7 * bw + 495.0
    if a < 18.0:
        return 17.5 * bw + 651.0
    if a < 30.0:
        return 15.3 * bw + 679.0
    return 11.6 * bw + 879.0


def REEREF(a):
    """Reference REE for age — Schofield at reference body weight."""
    return SCHOFIELD(a, LFMNORM(a) + FMNORM(a))


def EIREQREF(a):
    """Age-normative energy requirement of a reference individual, kcal/d."""
    return REEREF(a) * (1.0 + P["PALX"] * P["ACTREF"]) / (1.0 - P["FDIT"])


def LEPREF(a):
    return P["KSLEP"] * FMNORM(a) / P["KDLEP"]


def GUTREF(a):
    return EIREQREF(a) / P["KEMPT"]


def HTSD(a):
    return 1.9 + 0.031 * HTNORM(a)


def IGFSD(a):
    return 0.32 * IGFNORM(a)


# ----------------------------------------------------------------------------
# 3.  ALGEBRAIC LESION LAYER (Eq. A / Eq. B)
# ----------------------------------------------------------------------------
def pc13_activity(p):
    return 1.0 - p["DPC13"] * p["LES"]


def branch_loss(eps, pc13):
    return 1.0 - (1.0 + eps) / (1.0 + eps / pc13)


def branch_precursor(eps, pc13):
    return (1.0 + eps) / (pc13 + eps)


# ----------------------------------------------------------------------------
# 4.  CONTROLS
# ----------------------------------------------------------------------------
class Ctl(object):
    """The food environment, as four separable knobs.

    avail  — portions offered, as a multiple of the age-normative requirement
    cue    — external food-cue / access exposure driving food-seeking
    titr   — 0 = portions are NOT adjusted to the growth chart
             1 = caregivers titrate portions toward a weight target
    bwtgtr — the weight-for-age the family/clinic is actually holding
    fosa   — upper-airway vulnerability multiplier
    """
    __slots__ = ("avail", "cue", "fosa", "titr", "bwtgtr")

    def __init__(self, avail=1.0, cue=1.0, fosa=1.0, titr=0.0, bwtgtr=1.15):
        self.avail = avail
        self.cue = cue
        self.fosa = fosa
        self.titr = titr
        self.bwtgtr = bwtgtr


# ----------------------------------------------------------------------------
# 5.  CORE ALGEBRA — shared by derivs() and outputs()
# ----------------------------------------------------------------------------
def height_age(ht):
    """Age at which a reference individual has this height (weight-for-height)."""
    if ht <= _HTN[0]:
        return 0.0
    for i in range(len(_AGES) - 1):
        if _HTN[i] <= ht <= _HTN[i + 1]:
            if _HTN[i + 1] == _HTN[i]:
                return _AGES[i]
            f = (ht - _HTN[i]) / (_HTN[i + 1] - _HTN[i])
            return _AGES[i] + f * (_AGES[i + 1] - _AGES[i])
    return 18.0


def algebra(t, y, p, ctl, age0):
    """Every algebraic quantity the ODEs need, in dependency order.

    ORDER MATTERS and is not arbitrary: the melanocortin/oxytocin relay arm
    (x1) has to be known before the energy block, because the adiposity
    negative feedback on intake is TRANSDUCED by that relay, and the realized
    intake has to be known before the remaining satiety arms, because meal
    load is what drives distension, PYY and GLP-1.
    """
    a = age0 + t / 365.25
    A = {"a": a}
    FM = max(y[IX["FM"]], p["FMMIN"])
    LFM = max(y[IX["LFM"]], 1.0)
    BW = FM + LFM
    A["FM"], A["LFM"], A["BW"] = FM, LFM, BW
    pc13 = pc13_activity(p)
    A["pc13"] = pc13
    pws = p["LES"] > 0.5

    # ==== (1) drug concentrations — apparent, day-scale, AUC-calibrated ====
    A["CGH"] = y[IX["AGHC"]] / (p["VGH"] * BW) * 1000.0        # ug/L
    A["CLG"] = y[IX["ALGC"]] / (p["VLG"] * BW) * 1000.0        # ug/L
    A["CCB"] = y[IX["ACBC"]] / p["VCB"] * 1000.0               # ng/mL
    A["CDZ"] = y[IX["ADZC"]] / (p["VDZ"] * BW)                 # ug/mL
    A["CSM"] = y[IX["ASMC"]] / p["VSM"] * 1000.0               # ng/mL
    A["CSG"] = y[IX["ASGC"]] / p["VSG"] * 1000.0               # ng/mL
    A["COC"] = y[IX["AOCC"]] / p["VOC"] * 1000.0               # ng/mL
    A["CLV"] = y[IX["ALVC"]] / (p["VLV"] * BW) * 1000.0        # ng/mL
    A["CTS"] = y[IX["ATSC"]] / p["VTS"] * 100000.0             # ng/dL
    A["CMF"] = y[IX["AMFC"]] / p["VMF"]                        # ug/mL

    # ==== (2) receptor-level drug effects =================================
    CCB = A["CCB"]
    A["ECB"] = (p["ECBMAX"] * CCB / (p["ECB50"] + CCB)
                - p["EV1AMX"] * CCB / (p["EV1A50"] + CCB))   # BIPHASIC
    A["EDZ"] = p["EDZMAX"] * A["CDZ"] / (p["EDZ50"] + A["CDZ"])
    A["EDZI"] = p["EDZIMX"] * A["CDZ"] / (p["EDZI50"] + A["CDZ"])
    A["ESM"] = p["ESMMAX"] * A["CSM"] / (p["ESM50"] + A["CSM"])
    A["ESGG"] = p["ESGGE"] * A["CSG"] / (p["ESG50"] + A["CSG"])
    A["ESGS"] = p["ESGSAT"] * A["CSG"] / (p["ESG50"] + A["CSG"])
    A["ESGA"] = p["ESGAM"] * A["CSG"] / (p["ESG50"] + A["CSG"])
    A["EOCg"] = p["EOCG"] * A["COC"] / (p["EOC50"] + A["COC"])
    A["EOCgh"] = p["EOCGH"] * A["COC"] / (p["EOC50"] + A["COC"])
    A["EOCin"] = p["EOCINS"] * A["COC"] / (p["EOC50"] + A["COC"])
    A["EMF"] = p["EMFSI"] * A["CMF"] / (p["EMF50"] + A["CMF"])

    # ==== (3) somatotropic axis (reference-normalized) ====================
    ghn = GHNORM(a)
    fbI = (1.0 + IGFNORM(a) / p["KIGFFB"]) / (1.0 + y[IX["IGF1"]] / p["KIGFFB"])
    fbF = (1.0 + FMNORM(a) / p["KFMGH"]) / (1.0 + FM / p["KFMGH"])
    fbS = (1.0 + p["KSSTG"]) / (1.0 + p["KSSTG"] * max(0.05, y[IX["SST"]]))
    A["GHsec"] = (p["KCLGH"] * ghn * y[IX["GHRH"]]
                  * (p["FSOM"] if pws else 1.0)
                  * fbI * fbF * fbS / (1.0 + A["EOCgh"]))
    A["GHREL"] = (y[IX["GH"]] + p["FPOTGH"] * (A["CGH"] + A["CLG"])) / ghn
    A["GHX"] = max(0.0, A["GHREL"] - 1.0) / (1.0 + max(0.0, A["GHREL"] - 1.0))
    g = max(1e-4, A["GHREL"])
    # GH -> IGF-1 SATURATES (GH-binding protein and receptor limited), which is
    # why a 4-fold rise in GH exposure buys ~2 SDS of IGF-1 and not 6.
    A["IGFT"] = IGFNORM(a) * (g / (p["KGIGF"] + g)) * (p["KGIGF"] + 1.0)
    A["IGFdrv"] = y[IX["IGF1"]] / max(1.0, IGFNORM(a))
    A["exc"] = min(1.5, max(0.0, A["IGFdrv"] - 1.0))

    # ==== (4) the melanocortin -> PVN-oxytocin relay arm ==================
    # alpha-MSH is NOT a parallel satiety arm.  MC4R satiety signalling is
    # relayed through PVN oxytocin neurons, so it enters as a SATURATING INPUT
    # GAIN on the oxytocin arm and is bounded above by oxytocin availability.
    MCin = y[IX["AMSH"]] + A["ESM"]
    A["MCin"] = MCin
    A["relay"] = (MCin / (p["KREL"] + MCin)) * (p["KREL"] + 1.0)
    A["x1"] = max(1e-4, y[IX["OXT"]] * A["relay"] * (1.0 + A["ECB"]))
    # RELAY INTEGRITY.  Leptin's anorexigenic signal is transduced by this same
    # arcuate -> PVN relay, so the adiposity negative-feedback loop is
    # attenuated in PWS BY CONSTRUCTION: functional leptin resistance with
    # normal leptin levels and normal receptors, and no extra parameter.
    A["RI"] = min(1.20, A["x1"] / p["X1CTL"])

    # ==== (5) energy balance and what is actually eaten ===================
    REE = (REEREF(a) * (p["WREEL"] * LFM / LFMNORM(a)
                        + p["WREEF"] * FM / FMNORM(a))
           * (1.0 + p["KGHRMR"] * A["GHX"]))
    A["REE"] = REE
    fx = max(0.0, FM - FMNORM(a))          # only EXCESS fat costs mobility
    A["ACTt"] = (p["ACT0"] * y[IX["TONE"]] / p["TONEA"]
                 * (1.0 - p["KFATACT"] * fx / (fx + p["KFA"])))
    A["AT"] = REE * p["PALX"] * y[IX["ACT"]]
    EIREQ = EIREQREF(a)
    A["EIREQ"] = EIREQ
    LEP = y[IX["LEP"]]
    offs = (p["EIBF"] + p["KEIS"] * y[IX["SEEK"]]
            + p["KEIH"] * y[IX["HQ"]] / p["HQMAX"]
            - p["KLEPEI"] * A["RI"] * (LEP / LEPREF(a) - 1.0))
    offs *= (1.0 - p["ESGEI"] * A["ESGS"] / p["ESGSAT"]
             - p["EMFEI"] * A["EMF"] / p["EMFSI"])
    A["EIdrv"] = max(0.0, EIREQ * offs)
    # What is OFFERED.  titr = 0: portions are age-normative and take no notice
    # of the child's weight — the phase-2a situation, and the reason a NORMAL
    # intake is a surplus for a low-expenditure child.  titr = 1: caregivers
    # titrate portions toward a weight-FOR-HEIGHT target, the way a PWS clinic
    # does.  Then body WEIGHT is pinned and COMPOSITION is decided by the
    # lean-mass target — which is how growth hormone lowers %fat while barely
    # moving BMI.
    ha = height_age(y[IX["HT"]])
    BWTGT = ctl.bwtgtr * (LFMNORM(ha) + FMNORM(ha))
    A["BWTGT"] = BWTGT
    titr = 1.0 + ctl.titr * p["KTITR"] * (BWTGT / BW - 1.0)
    A["titr"] = min(1.50, max(0.45, titr))
    A["EIcap"] = ctl.avail * EIREQ * A["titr"]
    suck = SUCKPWS(a) if pws else 1.0
    A["EI"] = max(0.0, suck * min(A["EIdrv"], A["EIcap"], p["EICAPF"] * EIREQ))
    A["ASK"] = EIREQ * (p["EIBF"] + p["KEIS"] * y[IX["SEEK"]]
                        + p["KEIH"] * y[IX["HQ"]] / p["HQMAX"])
    A["FRUST"] = min(1.0, max(0.0, (A["ASK"] - A["EI"]) / max(1.0, A["ASK"])))
    A["REINF"] = 1.0 - A["FRUST"]
    A["DIT"] = p["FDIT"] * A["EI"]
    A["TEE"] = REE + A["AT"] + A["DIT"]
    A["ESUFF"] = min(1.0, max(0.45, A["EI"] / max(1.0, A["TEE"])))

    # ==== (6) meal load — the day-scale stand-in for meal-related signals ==
    # Note this is the INTAKE ratio, not the mean gastric content: on a day
    # scale the mean gastric pool cannot distinguish meal-related distension
    # from chronic overfill, and using it would turn PWS's delayed gastric
    # emptying into a permanent satiety bonus, which is not what PWS is.
    ml = min(1.30, max(0.50, A["EI"] / EIREQ))
    ml *= (1.0 + 0.25 * A["ESGG"] / p["ESGGE"])
    A["mealload"] = ml
    A["FASTs"] = max(0.35, 1.0 + p["KFASTA"] * (1.0 - ml))
    A["FASTg"] = max(0.35, 1.0 + p["KFASTG"] * (1.0 - ml))

    # ==== (7) prohormone synthesis and processing fluxes ==================
    lep_sat = ((LEP / (p["KLEPR0"] + LEP))
               / (LEPREF(a) / (p["KLEPR0"] + LEPREF(a))))
    A["SPOMC"] = p["SPOMC"] * (1.0 + p["KLEPP"] * (lep_sat - 1.0))
    A["SPOXT"] = p["SPOXT"] * (p["FOXTN"] if pws else 1.0)
    A["SPGHRH"] = p["SPGHRH"]
    phiG = y[IX["GLU"]] ** 2 / (p["KGGLU"] ** 2 + y[IX["GLU"]] ** 2)
    A["phiG"] = phiG
    A["SPINS"] = p["SPINS"] * phiG / 0.3959
    A["SPGHR"] = (p["SPGHR"] * (p["FGHRC"] if pws else 1.0) * A["FASTg"]
                  / (1.0 + y[IX["INS"]] / p["KIINS"]) / (1.0 + A["EOCg"]))
    A["fPOMC"] = p["KCPOMC"] * pc13 * y[IX["PPOMC"]]
    A["fPOXT"] = p["KCPOXT"] * pc13 * y[IX["PPOXT"]]
    A["fPGHRH"] = p["KCPGHRH"] * pc13 * y[IX["PPGHRH"]]
    A["fPINS"] = p["KCPINS"] * pc13 * y[IX["PPINS"]]
    A["fPGHR"] = p["KCPGHR"] * pc13 * y[IX["PPGHR"]]

    # ==== (8) the satiety integrator — a HARMONIC MEAN ====================
    # 1/SAT = sum(w_i/x_i): the weakest arm rate-limits, so an intact arm
    # cannot be over-driven to rescue a broken one.
    A["x2"] = max(1e-4, (p["FVAG"] if pws else 1.0) * ml)
    A["x3"] = max(1e-4, (p["FPYYS"] if pws else 1.0) * y[IX["PYY"]] / p["PYY0"])
    A["x4"] = max(1e-4, y[IX["GLP1E"]] / p["GLP10"] + A["ESGS"])
    A["x5"] = max(1e-4, (max(0.15, LEP) / LEPREF(a)) ** 0.28
                  * (max(1.0, y[IX["INS"]]) / 11.0) ** 0.22)
    A["SAT"] = 1.0 / (p["W1"] / A["x1"] + p["W2"] / A["x2"] + p["W3"] / A["x3"]
                      + p["W4"] / A["x4"] + p["W5"] / A["x5"])

    # ==== (9) orexigenic drive — the ghrelin arm SATURATES =================
    AGeff = y[IX["AG"]] / (1.0 + p["KLVAG"] * A["CLV"])
    A["AGeff"] = AGeff
    ghr = AGeff / (p["KAG50"] + AGeff)
    ghr0 = 350.0 / (p["KAG50"] + 350.0)
    A["GHRARM"] = ghr / ghr0
    A["gateHP"] = (1.0 / (1.0 + math.exp(-(a - p["AGEHP"]) / p["WAGEHP"]))
                   * (1.0 - p["FHP4"]
                      / (1.0 + math.exp(-(a - p["AGEHP4"]) / p["WAGEHP4"]))))
    A["DRVi"] = max(0.0, (max(0.02, y[IX["AGRP"]]) ** p["NAGRP"])
                    * (1.0 + p["KGHRD"] * (A["GHRARM"] - 1.0)) * A["gateHP"])
    A["DRVe"] = A["DRVi"] * ctl.cue
    A["SELFG"] = p["KSELF"] * (p["RFL"] + (1.0 - p["RFL"]) * A["REINF"])

    # ==== (10) growth and the lean-mass developmental target ==============
    SEXtot = y[IX["SEX"]] + A["CTS"]
    A["SEXtot"] = SEXtot
    A["SEXrel"] = ((SEXtot / (SEXtot + p["KSEXH"]))
                   / (620.0 / (620.0 + p["KSEXH"])))
    A["sexfac"] = min(1.30, A["SEXrel"] / max(0.03, SEXRELREF(a)))
    A["epi"] = 1.0 / (1.0 + math.exp((y[IX["BA"]] - p["BACLOSE"]) / 0.35))
    # growth velocity has a large GH-INDEPENDENT floor (GVFL); this is why
    # untreated PWS loses ~20 cm of adult height and not 60.
    gvf = p["GVFL"] + (1.0 - p["GVFL"]) * min(p["GVCAP"], max(0.15, A["IGFdrv"]))
    A["gvf"] = gvf
    A["GV"] = ((GVBASE(a) * gvf + GVSPURT(a) * A["sexfac"] * (gvf ** 0.5))
               * A["ESUFF"] * A["epi"])
    # Each lean-target factor equals 1.0 for the non-PWS reference at every
    # age, so the control tracks LFMNORM(a) by construction and every
    # deviation is a named mechanism.  Growth hormone and testosterone reach
    # lean mass THROUGH these factors — the model gives them no private path.
    fmod = ((p["WLIGF0"] + (1.0 - p["WLIGF0"]) * min(p["LIGFCAP"], A["IGFdrv"]))
            * (p["WLTON0"] + (1.0 - p["WLTON0"]) * y[IX["TONE"]] / p["TONEA"])
            * (A["ESUFF"] ** 0.5)
            * (1.0 - 0.5 * p["KSEXPART"] + 0.5 * p["KSEXPART"] * A["sexfac"]))
    A["fmod"] = min(p["FMODCAP"], fmod)
    A["LFMT"] = LFMNORM(a) * A["fmod"]
    return A


# ----------------------------------------------------------------------------
# 6.  DERIVATIVES
# ----------------------------------------------------------------------------
def derivs(t, y, p, ctl, age0):
    d = [0.0] * NST
    A = algebra(t, y, p, ctl, age0)
    a = A["a"]
    pws = p["LES"] > 0.5
    pc13 = A["pc13"]

    # --- prohormone pools ------------------------------------------------
    d[0] = A["SPOMC"] - (p["KCPOMC"] * pc13 + p["EPSPOMC"]) * y[0]
    d[1] = A["SPOXT"] - (p["KCPOXT"] * pc13 + p["EPSPOXT"]) * y[1]
    d[2] = A["SPGHRH"] - (p["KCPGHRH"] * pc13 + p["EPSPGHRH"]) * y[2]
    d[3] = A["SPINS"] - (p["KCPINS"] * pc13 + p["EPSPINS"]) * y[3]
    d[4] = A["SPGHR"] - (p["KCPGHR"] * pc13 + p["EPSPGHR"]) * y[4]

    # --- neuropeptide tones ---------------------------------------------
    d[5] = p["KSAMSH"] * A["fPOMC"] - p["KDAMSH"] * y[5]
    d[6] = p["KSOXT"] * A["fPOXT"] - p["KDOXT"] * y[6]
    # oxytocin receptor agonism acts on BOTH sides of the synapse:
    #   postsynaptic -> satiety arm x1 (in algebra())
    #   presynaptic  -> AgRP inhibition (here)
    oxt_eff = y[6] * (1.0 + A["ECB"])
    lepf = ((1.0 + LEPREF(a) / p["KLEPA"])
            / (1.0 + A["RI"] * y[IX["LEP"]] / p["KLEPA"]))
    d[7] = (p["SAGRP"] * A["FASTs"] * lepf / (1.0 + A["EDZ"] + A["ESGA"])
            - p["KDAGRP"] * y[7]
            * (1.0 + p["KINHM"] * y[5] + p["KINHO"] * oxt_eff))

    # --- ghrelin ----------------------------------------------------------
    fac = min(0.60, p["FACYL"] * (p["FACYLP"] if pws else 1.0))
    d[8] = p["KSAG"] * A["fPGHR"] * fac - (p["KCLAG"] + p["KDEACYL"]) * y[8]
    d[9] = (p["KSUAG"] * A["fPGHR"] * (1.0 - fac) + p["KDEACYL"] * y[8]
            - p["KCLUAG"] * y[9])

    # --- gut peptides -----------------------------------------------------
    d[10] = p["SPYY"] * A["mealload"] - p["KDPYY"] * y[10]
    d[11] = p["SGLP1"] * A["mealload"] - p["KDGLP1"] * y[11]

    # --- leptin -----------------------------------------------------------
    d[12] = p["KSLEP"] * A["FM"] - p["KDLEP"] * y[12]

    # --- body composition -------------------------------------------------
    # Lean mass follows a DEVELOPMENTAL TARGET modulated by the IGF-1 drive,
    # muscle tone, energy sufficiency and the anabolic (GH / sex-steroid)
    # partition factor.  Fat mass is then the ENERGY BUFFER: it takes
    # whatever the energy balance leaves after lean accretion is paid for.
    # Written this way, energy is exactly conserved and PWS's lean deficit is
    # DERIVED (low IGF-1 x low tone), not imposed.
    d[14] = (A["LFMT"] - y[14]) / p["TAULFM"]
    d[13] = (A["EI"] - A["TEE"] - p["RHOL"] * d[14]) / p["RHOF"]

    # --- gut pool, height, activity, bone age ------------------------------
    d[15] = (A["EI"] - p["KEMPT"] * (p["FGE"] if pws else 1.0)
             / (1.0 + A["ESGG"]) * y[15])
    d[16] = A["GV"] / 365.25
    d[17] = (A["ACTt"] - y[17]) / p["TAUACT"]
    d[18] = (1.0 / 365.25) * (1.0 + p["KBASEX"] * min(1.2, A["SEXrel"]))

    # --- food-seeking bistability ----------------------------------------
    sh = y[19] ** p["NSELF"]
    d[19] = (p["KSON"] * A["DRVe"] * (1.0 - y[19])
             + A["SELFG"] * sh / (p["KSHALF"] ** p["NSELF"] + sh) * (1.0 - y[19])
             - p["KSOFF"] * y[19] * A["SAT"])

    # --- hyperphagia score (HQ-CT) ---------------------------------------
    dh = A["DRVi"] ** p["NHQ"]
    HQt = p["HQMAX"] * (p["WHQS"] * y[19]
                        + p["WHQD"] * dh / (p["KHQ"] ** p["NHQ"] + dh))
    d[20] = (HQt - y[20]) / p["TAUHQ"]

    # --- behaviour --------------------------------------------------------
    FRUST = A["FRUST"]
    BEHt = (p["BEH0"] * p["FSUB"] * (1.0 + p["KBFR"] * FRUST)
            * (1.0 - p["KBOXT"] * (A["x1"] / p["X1REF"] - 1.0))
            * (0.55 + 0.45 * A["gateHP"]))
    d[21] = (max(0.0, BEHt) - y[21]) / p["TAUBEH"]

    # --- somatotropic axis ------------------------------------------------
    d[22] = p["KSGHRH"] * A["fPGHRH"] - p["KDGHRH"] * y[22]
    d[23] = (1.0 + p["KSSTI"] * (A["IGFdrv"] - 1.0) - y[23]) / p["TAUSST"]
    d[24] = A["GHsec"] - p["KCLGH"] * y[24]
    liver = 0.55 + 0.45 * A["ESUFF"]
    d[25] = p["KDIGF"] * (A["IGFT"] * liver - y[25])
    d[26] = p["KBP3"] * (max(1e-3, A["GHREL"]) ** 0.65) - p["KDBP3"] * y[26]

    # --- glucose / insulin ------------------------------------------------
    d[30] = (p["KSADPN"] * (p["FADPN"] if pws else 1.0)
             / (1.0 + A["FM"] / p["KFADP"]) - p["KDADPN"] * y[30])
    SI = (p["SI0"] * (1.0 + p["KADPNSI"] * (y[30] / 11.0 - 1.0))
          / (1.0 + p["KFMSI"] * A["FM"] / (A["FM"] + p["KFS"]))
          / (1.0 + p["KGHSI"] * A["GHX"])
          * (1.0 + A["EMF"] + p["ESGSI"] * A["ESGS"] / p["ESGSAT"]))
    SI = max(0.15, SI)
    A["SI"] = SI
    d[28] = (p["KSINS"] * A["fPINS"] / (1.0 + A["EDZI"]) / (1.0 + A["EOCin"])
             / p["VINS"] - p["KCLINS"] * y[28])
    EGP = p["KEGP"] / (1.0 + y[28] / p["KIEGP"])
    Ra = p["KGABS"] * p["KEMPT"] * (p["FGE"] if pws else 1.0) * y[15]
    Rd = (p["KGU0"] + p["KGUI"] * y[28] * SI) * y[27]
    d[27] = (EGP + Ra - Rd) / p["VG"]
    d[29] = ((y[27] * 1.14 + 46.7) / 28.7 - y[29]) / p["TAUA1C"]

    # --- gonadal axis / bone ----------------------------------------------
    bapub = p["BAPUB"] - (p["DMKRN3"] if pws else 0.0)
    gateP = 1.0 / (1.0 + math.exp(-(y[18] - bapub) / p["WPUB"]))
    d[31] = ((p["SGNRH"] * gateP * (p["FHYPO"] if pws else 1.0) - y[31])
             / p["TAUGN"])
    d[32] = p["KLH"] * y[31] - p["KDLH"] * y[32]
    d[33] = p["KSEX"] * y[32] - p["KDSEX"] * y[33]
    d[34] = (p["KBMDU"] * (max(0.02, A["SEXtot"]) / 620.0) ** 0.5
             * (1.0 + p["KBMDL"] * A["LFM"] / LFMNORM(a))
             * (1.0 + p["KBMDI"] * A["exc"])
             - p["KBMDR"] * y[34]
             * (1.0 + p["KBMDH"] * max(0.0, 1.0 - A["SEXtot"] / 620.0)))

    # --- upper airway (two clocks) ----------------------------------------
    d[36] = (A["exc"] - y[36]) / p["TAUAD"]
    d[35] = (p["KLYU"] * A["exc"] * (p["LMAX"] - y[35])
             - p["KLYD"] * (y[35] - LYNORM(a)) * (1.0 + p["KADL"] * y[36]))
    RMSt = ((A["LFM"] / LFMNORM(a)) ** p["NRMS"]) * (1.0 + p["KRMSIGF"] * A["exc"])
    d[37] = (p["RMS0"] * RMSt - y[37]) / p["TAURMS"]
    AHIobs = (p["KAHIL"] * max(0.0, y[35] - LYNORM(a))
              + p["KAHIF"] * max(0.0, A["FM"] / A["LFM"] - p["FLR0"]))
    AHIt = (p["AHI0"] + ctl.fosa * AHIobs - p["KAHIR"] * (y[37] - 1.0)
            + p["KAHIC"] * CAPNORM(a) * p["LES"] * (1.0 - 0.45 * A["GHX"]))
    d[38] = (max(0.3, AHIt) - y[38]) / p["TAUAHI"]

    # --- tone / scoliosis / thermoregulation ------------------------------
    # hypotonia improves through infancy but never resolves in PWS
    tone0 = (p["TONEA"] * (p["FTONE0"]
                           + (p["FTONEP"] - p["FTONE0"]) * min(1.0, a / 5.0))
             if pws else p["TONEA"])
    TONEt = (tone0 * (1.0 + p["KTIGF"] * A["exc"])
             * (A["LFM"] / LFMNORM(a)) ** 0.40)
    d[39] = (TONEt - y[39]) / p["TAUTONE"]
    d[40] = (p["KCOB"] * (A["GV"] / 365.25)
             * (p["TONEA"] / max(0.25, y[39])) ** p["NCOB"])
    TEMPT = p["LES"] * (0.35 + 0.65 * math.exp(-a / 6.0))
    d[41] = (TEMPT - y[41]) / p["TAUTEMP"]

    # --- drug PK -----------------------------------------------------------
    d[42] = -p["KAGH"] * y[42]
    d[43] = p["KAGH"] * y[42] - p["KEGH"] * y[43]
    d[44] = -p["KALG"] * y[44]
    d[45] = p["KALG"] * y[44] - p["KELG"] * y[45]
    d[46] = -p["KACB"] * y[46]
    d[47] = p["KACB"] * y[46] * p["FCB"] - p["KECB"] * y[47]
    d[48] = -p["KADZ"] * y[48]
    d[49] = p["KADZ"] * y[48] - p["KEDZ"] * y[49]
    d[50] = -p["KASM"] * y[50]
    d[51] = p["KASM"] * y[50] - p["KESM"] * y[51]
    d[52] = -p["KASG"] * y[52]
    d[53] = p["KASG"] * y[52] - p["KESG"] * y[53]
    d[54] = -p["KAOC"] * y[54]
    d[55] = p["KAOC"] * y[54] - p["KEOC"] * y[55]
    d[56] = -p["KALV"] * y[56]
    d[57] = p["KALV"] * y[56] - p["KELV"] * y[57]
    d[58] = -p["KATS"] * y[58]
    d[59] = p["KATS"] * y[58] - p["KETS"] * y[59]
    d[60] = -p["KAMF"] * y[60]
    d[61] = p["KAMF"] * y[60] - p["KEMF"] * y[61]

    d[62] = A["EI"]
    d[63] = A["TEE"]
    return d


# ----------------------------------------------------------------------------
# 7.  OUTPUTS
# ----------------------------------------------------------------------------
def outputs(t, y, p, ctl, age0):
    A = algebra(t, y, p, ctl, age0)
    a = A["a"]
    o = dict(A)
    o["age"] = a
    o["PBF"] = 100.0 * A["FM"] / A["BW"]
    o["HT"] = y[IX["HT"]]
    o["BMI"] = A["BW"] / (y[IX["HT"]] / 100.0) ** 2
    o["HTSDS"] = (y[IX["HT"]] - HTNORM(a)) / HTSD(a)
    o["IGF1"] = y[IX["IGF1"]]
    o["IGFSDS"] = (y[IX["IGF1"]] - IGFNORM(a)) / IGFSD(a)
    o["AG"] = y[IX["AG"]]
    o["UAG"] = y[IX["UAG"]]
    o["AGUAG"] = y[IX["AG"]] / max(1e-6, y[IX["UAG"]])
    o["AGUAGeff"] = A["AGeff"] / max(1e-6, y[IX["UAG"]])
    o["HQ"] = y[IX["HQ"]]
    o["SEEK"] = y[IX["SEEK"]]
    o["BEH"] = y[IX["BEH"]]
    o["AHI"] = y[IX["AHI"]]
    o["LYMPH"] = y[IX["LYMPH"]]
    o["RMS"] = y[IX["RMS"]]
    o["GLU"] = y[IX["GLU"]]
    o["INS"] = y[IX["INS"]]
    o["HBA1C"] = y[IX["HBA1C"]]
    o["ADPN"] = y[IX["ADPN"]]
    o["HOMA"] = y[IX["GLU"]] * y[IX["INS"]] / 405.0
    o["BMD"] = y[IX["BMD"]]
    o["COBB"] = y[IX["COBB"]]
    o["TONE"] = y[IX["TONE"]]
    o["OXT"] = y[IX["OXT"]]
    o["AMSH"] = y[IX["AMSH"]]
    o["AGRP"] = y[IX["AGRP"]]
    o["GH"] = y[IX["GH"]]
    o["LH"] = y[IX["LH"]]
    o["SEX"] = A["SEXtot"]
    o["BA"] = y[IX["BA"]]
    o["TEEfr"] = A["TEE"] / A["EIREQ"]
    o["EIfr"] = A["EI"] / A["EIREQ"]

    o["WTSDS"] = ((A["BW"] - (LFMNORM(a) + FMNORM(a)))
                  / (0.13 * (LFMNORM(a) + FMNORM(a))))
    o["PHASE"] = miller_phase(a)
    return o


def miller_phase(a):
    if a < 0.75:
        return "1a"
    if a < 2.1:
        return "1b"
    if a < 4.5:
        return "2a"
    if a < 8.0:
        return "2b"
    if a < 30.0:
        return "3"
    return "4"


# ----------------------------------------------------------------------------
# 8.  INTEGRATOR + DOSING
# ----------------------------------------------------------------------------
class Regimen(object):
    def __init__(self, cmt, amt, interval, start=0.0, stop=1e12, per_kg=False,
                 n_per_day=1):
        self.cmt = IX[cmt]
        self.amt = amt
        self.interval = interval
        self.start = start
        self.stop = stop
        self.per_kg = per_kg
        self.n_per_day = n_per_day

    def times(self, tmax):
        ts = []
        t = self.start
        step = self.interval / self.n_per_day
        while t < min(tmax, self.stop) + 1e-9:
            ts.append(t)
            t += step
        return ts


class Sim(object):
    def __init__(self, p, y0, age0):
        self.p = dict(p)
        self.y = list(y0)
        self.age0 = age0
        self.t = 0.0
        self.trace = []
        self.last_ctl = Ctl()

    def snapshot(self):
        s = Sim(self.p, self.y, self.age0)
        s.t = self.t
        s.last_ctl = self.last_ctl
        return s

    def advance(self, duration, dt, ctl=None, regimens=(), record=0.0, tag=""):
        p = self.p
        ctl = ctl or Ctl()
        self.last_ctl = ctl
        if duration <= 0.0:
            return self
        events = []
        for r in regimens:
            for tt in r.times(duration):
                events.append((round(tt / dt) * dt, r))
        events.sort(key=lambda e: e[0])
        ei = 0
        nsteps = int(round(duration / dt))
        rec_every = int(round(record / dt)) if record > 0 else 0
        t0 = self.t
        for k in range(nsteps):
            tk = k * dt
            while ei < len(events) and events[ei][0] <= tk + 1e-9:
                _, r = events[ei]
                bw = max(2.0, self.y[IX["FM"]] + self.y[IX["LFM"]])
                self.y[r.cmt] += r.amt * (bw if r.per_kg else 1.0)
                ei += 1
            self.y = rk4(t0 + tk, self.y, dt, p, ctl, self.age0)
            if rec_every and (k + 1) % rec_every == 0:
                o = outputs(t0 + tk + dt, self.y, p, ctl, self.age0)
                o["t"] = t0 + tk + dt
                o["tag"] = tag
                self.trace.append(o)
        self.t += duration
        return self

    def out(self, ctl=None):
        return outputs(self.t, self.y, self.p, ctl or self.last_ctl, self.age0)


def rk4(t, y, h, p, ctl, age0):
    k1 = derivs(t, y, p, ctl, age0)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(NST)]
    k2 = derivs(t + 0.5 * h, y2, p, ctl, age0)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(NST)]
    k3 = derivs(t + 0.5 * h, y3, p, ctl, age0)
    y4 = [y[i] + h * k3[i] for i in range(NST)]
    k4 = derivs(t + h, y4, p, ctl, age0)
    out = [0.0] * NST
    for i in range(NST):
        v = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        out[i] = v if v > 0.0 else 0.0
    out[IX["SEEK"]] = min(1.0, out[IX["SEEK"]])
    out[IX["TONE"]] = min(1.20, max(0.05, out[IX["TONE"]]))
    return out


# ----------------------------------------------------------------------------
# 9.  INITIAL CONDITIONS (birth)
# ----------------------------------------------------------------------------
def init_birth(p):
    pc = pc13_activity(p)
    fox = p["FOXTN"] if p["LES"] > 0.5 else 1.0
    y = [0.0] * NST
    y[IX["PPOMC"]] = 1.0 / (pc + p["EPSPOMC"])
    y[IX["PPOXT"]] = fox / (pc + p["EPSPOXT"])
    y[IX["PPGHRH"]] = 1.0 / (pc + p["EPSPGHRH"])
    y[IX["PPINS"]] = 1.0 / (pc + p["EPSPINS"])
    y[IX["PPGHR"]] = 1.0 / (pc + p["EPSPGHR"])
    y[IX["AMSH"]] = p["KSAMSH"] / p["KDAMSH"] * pc * y[IX["PPOMC"]]
    y[IX["OXT"]] = p["KSOXT"] / p["KDOXT"] * pc * y[IX["PPOXT"]]
    y[IX["AGRP"]] = 1.0
    y[IX["AG"]] = 350.0
    y[IX["UAG"]] = 1250.0
    y[IX["PYY"]] = p["PYY0"]
    y[IX["GLP1E"]] = p["GLP10"]
    y[IX["LEP"]] = LEPREF(0.0)
    y[IX["FM"]] = 0.5
    y[IX["LFM"]] = 2.7
    y[IX["GUT"]] = GUTREF(0.0)
    y[IX["HT"]] = 50.0
    y[IX["ACT"]] = 0.60
    y[IX["SEEK"]] = 0.05
    y[IX["HQ"]] = 1.0
    y[IX["BEH"]] = 8.0
    y[IX["GHRH"]] = p["KSGHRH"] / p["KDGHRH"] * pc * y[IX["PPGHRH"]]
    y[IX["SST"]] = 1.0
    y[IX["GH"]] = GHNORM(0.0)
    y[IX["IGF1"]] = 50.0
    y[IX["IGFBP3"]] = 2.0
    y[IX["GLU"]] = 82.0
    y[IX["INS"]] = 8.0
    y[IX["HBA1C"]] = 4.6
    y[IX["ADPN"]] = 12.0
    y[IX["SEX"]] = 8.0
    y[IX["BMD"]] = 0.42
    y[IX["LYMPH"]] = LYNORM(0.0)
    y[IX["RMS"]] = 1.0
    y[IX["AHI"]] = 2.0
    y[IX["TONE"]] = (p["TONEA"] * p["FTONE0"]) if p["LES"] > 0.5 else p["TONEA"]
    y[IX["COBB"]] = 2.0
    y[IX["TEMP"]] = 1.0 if p["LES"] > 0.5 else 0.0
    return y


def pars(**kw):
    q = dict(P)
    q.update(kw)
    return q


DT_LONG = 0.125
DT_DRUG = 1.0 / 24.0
GHRX = lambda: [Regimen("AGHD", 0.035, 1.0, per_kg=True)]


def CT(avail, cue, titr=1.0, bwtgtr=1.15, fosa=1.0):
    return Ctl(avail=avail, cue=cue, titr=titr, bwtgtr=bwtgtr, fosa=fosa)


# the four food environments used throughout the report
ENV_CTRL = dict(avail=1.40, cue=1.00, titr=0.0)     # non-PWS, self-regulated
ENV_FREE = dict(avail=1.40, cue=1.00, titr=0.0)     # PWS, unmanaged
ENV_MGMT = dict(avail=1.00, cue=0.28, titr=1.0)     # standard PWS management
ENV_TIGHT = dict(avail=0.90, cue=0.10, titr=1.0)    # strict environment


def run_to(les=1.0, avail=1.0, cue=0.28, age=12.0, dt=None, gh_from=None,
           avail_early=1.40, titr=1.0, bwtgtr=1.15, record=0.0, **kw):
    """Birth -> `age`; free access until 2 y, then (avail, cue); optional GH."""
    dt = dt or DT_LONG
    p = pars(LES=les, **kw)
    s = Sim(p, init_birth(p), 0.0)
    s.advance(2.0 * 365.25, dt, Ctl(avail=avail_early, cue=1.0, titr=0.0),
              record=record, tag="infancy")
    if gh_from is None:
        s.advance((age - 2.0) * 365.25, dt, CT(avail, cue, titr, bwtgtr),
                  record=record, tag="childhood")
    else:
        if gh_from > 2.0:
            s.advance((gh_from - 2.0) * 365.25, dt, CT(avail, cue, titr, bwtgtr),
                      record=record, tag="pre-GH")
        s.advance((age - max(2.0, gh_from)) * 365.25, dt,
                  CT(avail, cue, titr, bwtgtr), regimens=GHRX(),
                  record=record, tag="on-GH")
    return s


# ----------------------------------------------------------------------------
# 10.  REPORT
# ----------------------------------------------------------------------------
def hr(c="-", n=78):
    return c * n


def seek_rhs_g(S, DRVe, SAT, selfg):
    sh = S ** P["NSELF"]
    return (P["KSON"] * DRVe * (1.0 - S)
            + selfg * sh / (P["KSHALF"] ** P["NSELF"] + sh) * (1.0 - S)
            - P["KSOFF"] * S * SAT)


def fixed_points_g(DRVe, SAT, selfg, n=8001):
    """All roots of dSEEK/dt on [0,1], each tagged stable / SADDLE."""
    fps = []
    prev = seek_rhs_g(0.0, DRVe, SAT, selfg)
    for i in range(1, n):
        S = i / (n - 1.0)
        cur = seek_rhs_g(S, DRVe, SAT, selfg)
        if (prev < 0.0) != (cur < 0.0):
            lo, hi = (i - 1) / (n - 1.0), S
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if ((seek_rhs_g(lo, DRVe, SAT, selfg) < 0.0)
                        != (seek_rhs_g(mid, DRVe, SAT, selfg) < 0.0)):
                    hi = mid
                else:
                    lo = mid
            Sf = 0.5 * (lo + hi)
            stab = ("stable" if seek_rhs_g(Sf + 1e-4, DRVe, SAT, selfg) < 0.0
                    else "SADDLE")
            fps.append((Sf, stab))
        prev = cur
    return fps


def main():
    fast = "--fast" in sys.argv
    dt = 0.25 if fast else DT_LONG
    dtd = 1.0 / 12.0 if fast else DT_DRUG
    out = []
    W = out.append
    W("=" * 78)
    W("PRADER-WILLI SYNDROME (PWS) QSP MODEL — REFERENCE OUTPUT")
    W("pws_reference_model.py · pure-Python RK4 · %d states · dt = %.4f d"
      % (NST, dt))
    W("=" * 78)

    # ================= [1] the lesion, algebraically ====================
    pc = pc13_activity(P)
    W("")
    W("[1] PC1/3 BRANCH ALGEBRA (Eq. A / Eq. B) at PC13 = %.3f" % pc)
    W(hr())
    W("%-15s %7s %9s %10s %10s %10s" % ("branch", "eps_i", "loss L_i",
                                        "product", "precursor", "total IR"))
    br = {}
    for nm, eps in (("POMC>aMSH", P["EPSPOMC"]), ("proOXT>OXT", P["EPSPOXT"]),
                    ("proGHRH>GHRH", P["EPSPGHRH"]),
                    ("proins>insulin", P["EPSPINS"]),
                    ("proghrelin>AG", P["EPSPGHR"])):
        L = branch_loss(eps, pc)
        R = branch_precursor(eps, pc)
        prod = 1.0 / (1.0 + eps / pc)
        prod_n = 1.0 / (1.0 + eps)
        W("%-15s %7.2f %9.4f %10.4f %10.4f %10.6f"
          % (nm, eps, L, 1.0 - L, R,
             (prod + (1.0 - prod)) / (prod_n + (1.0 - prod_n))))
        br[nm] = (L, R)
    W("")
    W("  product loss rises with eps_i while precursor accumulation FALLS")
    W("  with eps_i -> the branches that lose product are exactly the")
    W("  branches that do NOT accumulate precursor.  Last column: an assay")
    W("  cross-reacting with precursor + product measures the synthesis rate")
    W("  and is EXACTLY blind to PC13 (1.000000 for every branch), which is")
    W("  why plasma/CSF 'oxytocin' studies in PWS disagree with each other.")
    W("")
    W("  proinsulin : insulin ratio vs control  = %.3f x"
      % (br["proins>insulin"][1] / (1.0 - br["proins>insulin"][0])))
    W("  branch failure ordering, from ONE number each (not five fits):")
    W("    %s" % " > ".join("%s (%.0f%%)" % (k, 100 * v[0])
                            for k, v in sorted(br.items(),
                                               key=lambda kv: -kv[1][0])))

    # ================= [2] natural history ==============================
    W("")
    W("[2] NATURAL HISTORY at 25 y")
    W(hr())
    ctrl = run_to(les=0.0, age=25.0, dt=dt, **ENV_CTRL)
    pfree = run_to(les=1.0, age=25.0, dt=dt, **ENV_FREE)
    psec = run_to(les=1.0, age=25.0, dt=dt, **ENV_MGMT)
    pghs = run_to(les=1.0, age=25.0, gh_from=1.0, dt=dt, **ENV_MGMT)
    W("%-26s %6s %6s %6s %6s %6s %6s %6s %6s"
      % ("", "BW", "BMI", "%fat", "HT", "HTsds", "HQ-CT", "AHI", "TEE/rq"))
    for nm, s, c in (("non-PWS control", ctrl, CT(1.40, 1.0, 0.0)),
                     ("PWS free access", pfree, CT(1.40, 1.0, 0.0)),
                     ("PWS food-secure", psec, CT(1.00, 0.28, 1.0)),
                     ("PWS secure + GH from 1y", pghs, CT(1.00, 0.28, 1.0))):
        o = s.out(c)
        W("%-26s %6.1f %6.1f %6.1f %6.1f %6.2f %6.1f %6.1f %6.3f"
          % (nm, o["BW"], o["BMI"], o["PBF"], o["HT"], o["HTSDS"], o["HQ"],
             o["AHI"], o["TEEfr"]))
    W("")
    W("  energy arithmetic at 8 y")
    W("  %-18s %7s %7s %6s %8s %8s %8s"
      % ("", "LFM", "REE", "ACT", "TEE", "req", "TEE/req"))
    for nm, les, env, ghf in (("control", 0.0, ENV_CTRL, None),
                              ("PWS, no GH", 1.0, ENV_MGMT, None),
                              ("PWS, on GH", 1.0, ENV_MGMT, 1.0)):
        s8 = run_to(les=les, age=8.0, gh_from=ghf, dt=dt, **env)
        o = s8.out(CT(env["avail"], env["cue"], env["titr"]))
        W("  %-18s %7.1f %7.0f %6.2f %8.0f %8.0f %8.3f"
          % (nm, o["LFM"], o["REE"], s8.y[IX["ACT"]], o["TEE"], o["EIREQ"],
             o["TEEfr"]))
    W("")
    W("  the 60-70% of predicted total energy expenditure reported in PWS is")
    W("  not a parameter here: it is (low lean mass) x (low activity).")

    # ================= [3] Miller phases ================================
    W("")
    W("[3] MILLER NUTRITIONAL PHASES — weight gain PRECEDES hyperphagia")
    W(hr())
    s = Sim(pars(LES=1.0), init_birth(pars(LES=1.0)), 0.0)
    W("%6s %6s %7s %7s %7s %7s %7s %7s %7s"
      % ("age", "phase", "BW", "%fat", "EI/req", "wtSDS", "HQ-CT", "SEEK", "AG"))
    prev = 0.0
    cinf = CT(1.00, 1.00, 0.0)   # normal portions, NOT titrated
    for tgt in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 14.0,
                20.0):
        s.advance((tgt - prev) * 365.25, dt, cinf)
        prev = tgt
        o = s.out(cinf)
        W("%6.1f %6s %7.1f %7.1f %7.3f %7.2f %7.1f %7.3f %7.0f"
          % (o["age"], o["PHASE"], o["BW"], o["PBF"], o["EIfr"], o["WTSDS"],
             o["HQ"], o["SEEK"], o["AG"]))
    W("")
    W("  EI/req stays BELOW 1 through phase 1a (oro-motor limitation), then")
    W("  crosses it; weight-for-age SDS turns positive well before HQ-CT")
    W("  rises, because the expenditure lesion makes a normal intake a")
    W("  surplus.  The hyperphagia gate itself is a developmental clock, so")
    W("  the phase sequence needs no change in any appetite parameter.")

    # ================= [4] satiety arms ================================
    W("")
    W("[4] SATIETY INTEGRATOR — harmonic mean; MC arm in series with OXT arm")
    W(hr())
    ref = run_to(les=1.0, age=12.0, gh_from=1.0, dt=dt, **ENV_MGMT)
    c12 = CT(1.00, 0.28, 1.0)
    o = ref.out(c12)
    cref = run_to(les=0.0, age=12.0, dt=dt, **ENV_CTRL).out(CT(1.40, 1.0, 0.0))
    W("  %-32s %9s %9s" % ("arm", "control", "PWS @12y"))
    for k, lab in (("x1", "PVN-oxytocin relay   (w=0.44)"),
                   ("x2", "vagal / distension   (w=0.16)"),
                   ("x3", "PYY3-36 -> Y2        (w=0.13)"),
                   ("x4", "GLP-1 -> GLP1R       (w=0.14)"),
                   ("x5", "leptin / insulin CNS (w=0.13)")):
        W("  %-32s %9.4f %9.4f" % (lab, cref[k], o[k]))
    for k, lab in (("SAT", "SAT (harmonic mean)"), ("OXT", "oxytocin tone"),
                   ("AMSH", "alpha-MSH tone"), ("AGRP", "AgRP tone"),
                   ("DRVi", "internal orexigenic drive")):
        W("  %-32s %9.4f %9.4f" % (lab, cref[k], o[k]))
    W("")

    def satx(x1=None, x2=None, x3=None, x4=None, x5=None):
        return 1.0 / (P["W1"] / (x1 or o["x1"]) + P["W2"] / (x2 or o["x2"])
                      + P["W3"] / (x3 or o["x3"]) + P["W4"] / (x4 or o["x4"])
                      + P["W5"] / (x5 or o["x5"]))

    base = o["SAT"]
    mcsat = o["OXT"] * ((3.72 / (P["KREL"] + 3.72)) * (P["KREL"] + 1.0))
    W("  %-44s %8s %8s" % ("perturbation", "SAT", "dSAT"))
    for lab, v in (("baseline", base),
                   ("MC4R arm driven to saturation", satx(x1=mcsat)),
                   ("OXTR arm x2.1 (carbetocin-like)", satx(x1=o["x1"] * 2.1)),
                   ("PYY arm normalized", satx(x3=1.0)),
                   ("GLP-1 arm doubled", satx(x4=2.0 * o["x4"])),
                   ("ALL non-oxytocin arms normalized",
                    satx(x2=1.0, x3=1.0, x4=1.0, x5=1.0))):
        W("  %-44s %8.4f %8.4f" % (lab, v, v - base))
    W("")
    W("  normalizing every arm EXCEPT the oxytocin arm leaves SAT at %.4f of"
      % satx(x2=1.0, x3=1.0, x4=1.0, x5=1.0))
    W("  a possible 1.0.  A harmonic mean cannot be rescued by its healthy")
    W("  terms, and MC4R agonism enters THROUGH the broken term — bounded")
    W("  above by the oxytocin availability PWS does not have.")

    # ================= [5] ghrelin saturation ===========================
    W("")
    W("[5] THE GHRELIN ARM IS ALREADY SATURATED")
    W(hr())
    W("  %10s %11s %11s" % ("AG pg/mL", "occupancy", "arm/control"))
    ghr0 = 350.0 / (P["KAG50"] + 350.0)
    for ag in (200, 350, 500, 700, o["AG"], 1000, 1400):
        occ = ag / (P["KAG50"] + ag)
        W("  %10.0f %11.4f %11.4f" % (ag, occ, occ / ghr0))
    W("")
    W("  model PWS AG = %.0f pg/mL = %.2f x control; AG:UAG = %.3f vs %.3f"
      % (o["AG"], o["AG"] / cref["AG"], o["AGUAG"], cref["AGUAG"]))
    W("  that %.2f-fold elevation buys only %+.1f%% of receptor arm, so"
      % (o["AG"] / cref["AG"], 100.0 * (o["GHRARM"] - 1.0)))
    W("  abolishing it entirely cannot return more than %.1f%% of the drive."
      % (100.0 * P["KGHRD"] * (o["GHRARM"] - 1.0)))
    W("")
    W("  NOTE the sign structure: PC1/3 deficiency REDUCES proghrelin")
    W("  processing by %.1f%% and yet PWS is hyperghrelinaemic, because the"
      % (100.0 * branch_loss(P["EPSPGHR"], pc)))
    W("  same lesion throttles proinsulin processing and insulin is")
    W("  ghrelin's brake.  The lesion reaches ghrelin twice, with opposite")
    W("  signs, and the net sign is decided by which loop has more gain.")

    # ================= [6] bistability ==================================
    W("")
    W("[6] FOOD-SEEKING BISTABILITY — food security restores the low branch")
    W(hr())

    DRVi, SATp = o["DRVi"], o["SAT"]
    SELFG = o["SELFG"]

    def fixed_points(DRVe, SAT):
        return fixed_points_g(DRVe, SAT, SELFG)
    W("  PWS @12y: DRVi = %.4f, SAT = %.4f, reinforcement gain = %.4f"
      % (DRVi, SATp, SELFG))
    W("  control : internal drive DRVi = %.4f, SAT = %.4f"
      % (cref["DRVi"], cref["SAT"]))
    W("")
    W("  %-14s %8s   %s" % ("CUE", "DRVe", "fixed points of dSEEK/dt"))
    for cue in (1.00, 0.60, 0.45, 0.35, 0.28, 0.15):
        fps = fixed_points(DRVi * cue, SATp)
        W("  CUE = %.2f    %8.3f   %s"
          % (cue, DRVi * cue,
             "  ".join("%.4f(%s)" % (f, st) for f, st in fps)))
    lo_c, hi_c = 0.02, 3.0
    for _ in range(60):
        mid = 0.5 * (lo_c + hi_c)
        if len(fixed_points(DRVi * mid, SATp)) >= 3:
            lo_c = mid
        else:
            hi_c = mid
    W("")
    W("  SADDLE-NODE at CUE* = %.4f.  Above it the low branch does not exist"
      % (0.5 * (lo_c + hi_c)))
    W("  and the patient is FORCED onto the food-seeking state — no dose of")
    W("  anything can hold a state the vector field has deleted.  This is")
    W("  the model's account of why environmental food security outperforms")
    W("  every pharmacology in the panel below, and is not commensurable")
    W("  with it: it changes which states exist, not how deep they are.")
    W("  Control fixed points (DRVe = %.3f, SAT = %.3f): %s"
      % (cref["DRVi"], cref["SAT"],
         "  ".join("%.4f(%s)" % (f, st)
                   for f, st in fixed_points_g(cref["DRVi"], cref["SAT"],
                                               cref["SELFG"]))))

    # ================= [7] growth hormone ===============================
    W("")
    W("[7] GROWTH HORMONE — one drug, four clocks")
    W(hr())
    W("%-22s %6s %6s %6s %6s %6s %6s %6s %6s"
      % ("PWS @14y", "HT", "HTsds", "%fat", "LFM", "IGF1", "IGFsds", "HQ-CT",
         "Cobb"))
    for nm, ghf in (("no GH", None), ("GH from 8.0 y", 8.0),
                    ("GH from 4.0 y", 4.0), ("GH from 1.0 y", 1.0)):
        s14 = run_to(les=1.0, age=14.0, gh_from=ghf, dt=dt, **ENV_MGMT)
        oo = s14.out(c12)
        W("%-22s %6.1f %6.2f %6.1f %6.1f %6.0f %6.2f %6.1f %6.1f"
          % (nm, oo["HT"], oo["HTSDS"], oo["PBF"], oo["LFM"], oo["IGF1"],
             oo["IGFSDS"], oo["HQ"], oo["COBB"]))
    W("")
    W("  GH moves height / lean / fat and does NOT move HQ-CT: the")
    W("  somatotropic and satiety branches share the lesion, not a path.")
    W("  Scoliosis: GH raises growth velocity (worse) and muscle tone")
    W("  (better) and the two nearly cancel — which is what the trials found.")

    W("")
    W("  THE AIRWAY WINDOW — AHI after starting GH at 6 y")
    base6 = run_to(les=1.0, age=6.0, dt=dt, **ENV_MGMT)
    for fosa, label in ((1.0, "typical airway"), (1.9, "vulnerable airway")):
        W("  --- %s (FOSA = %.1f) ---" % (label, fosa))
        W("  %6s %8s %8s %8s %8s %8s"
          % ("week", "IGFsds", "LYMPH", "RMS", "%fat", "AHI"))
        s = base6.snapshot()
        cc = CT(1.00, 0.28, 1.0, fosa=fosa)
        prevw = 0
        for wk in (0, 2, 4, 6, 8, 12, 26, 52, 104):
            if wk > prevw:
                s.advance((wk - prevw) * 7.0, dtd, cc, regimens=GHRX())
                prevw = wk
            oo = s.out(cc)
            W("  %6d %8.2f %8.3f %8.3f %8.1f %8.2f"
              % (wk, oo["IGFSDS"], oo["LYMPH"], oo["RMS"], oo["PBF"],
                 oo["AHI"]))
    W("")
    W("  the peak sits at weeks 4-8 and is not fitted: tau_LYMPH ~ 20 d <")
    W("  tau_RMS = %.0f d < tau_fat (months).  This is the mechanistic"
      % P["TAURMS"])
    W("  content of 'polysomnography before, and 6-8 weeks after, starting")
    W("  growth hormone', and of not starting GH during a respiratory")
    W("  infection: airway vulnerability multiplies exactly the window and")
    W("  not the plateau, so the same drug is dangerous for 8 weeks and")
    W("  protective thereafter.")

    # ================= [8] the drug panel ===============================
    W("")
    W("[8] DRUG PANEL — biomarker vs endpoint (PWS aged 12 y, food-secure,")
    W("    every arm carrying background growth hormone)")
    W(hr())
    base12 = ref
    b = base12.out(c12)
    ghmt = GHRX()[0]

    def trial(name, weeks, regimens, cue=0.28, bg_gh=True):
        s = base12.snapshot()
        cc = CT(1.00, cue, 1.0)
        rg = list(regimens) + ([ghmt] if bg_gh else [])
        s.advance(weeks * 7.0, dtd, cc, regimens=rg)
        return name, s.out(cc), weeks

    rows = [trial("placebo", 13, [])]
    rows.append(trial("carbetocin 3.2 mg TID", 8,
                      [Regimen("ACBD", 3.2, 1.0, n_per_day=3)]))
    rows.append(trial("carbetocin 9.6 mg TID", 8,
                      [Regimen("ACBD", 9.6, 1.0, n_per_day=3)]))
    rows.append(trial("DCCR 5.1 mg/kg QD", 13,
                      [Regimen("ADZD", 5.1, 1.0, per_kg=True)]))
    rows.append(trial("setmelanotide 3 mg QD", 12, [Regimen("ASMD", 3.0, 1.0)]))
    rows.append(trial("livoletide 60 ug/kg QD", 12,
                      [Regimen("ALVD", 0.060, 1.0, per_kg=True)]))
    rows.append(trial("octreotide LAR 30 mg q28d", 12,
                      [Regimen("AOCD", 30.0, 28.0)]))
    rows.append(trial("octreotide LAR, GH stopped", 12,
                      [Regimen("AOCD", 30.0, 28.0)], bg_gh=False))
    rows.append(trial("GH stopped (comparator)", 12, [], bg_gh=False))
    rows.append(trial("semaglutide 2.4 mg QW", 13, [Regimen("ASGD", 2.4, 7.0)]))
    rows.append(trial("metformin 1500 mg QD", 13,
                      [Regimen("AMFD", 1500.0, 1.0)]))
    W("%-27s %7s %7s %8s %7s %7s %7s"
      % ("arm", "dHQ-CT", "dAG %", "dAG:UAG", "dIGF1%", "d%fat", "dHbA1c"))
    for nm, oo, wk in rows:
        W("%-27s %7.2f %7.1f %8.1f %7.1f %7.2f %7.2f"
          % (nm, oo["HQ"] - b["HQ"], 100.0 * (oo["AG"] / b["AG"] - 1.0),
             100.0 * (oo["AGUAGeff"] / b["AGUAGeff"] - 1.0),
             100.0 * (oo["IGF1"] / b["IGF1"] - 1.0), oo["PBF"] - b["PBF"],
             oo["HBA1C"] - b["HBA1C"]))
    pl = rows[0][1]
    W("")
    W("  placebo-corrected dHQ-CT next to the biomarker each arm 'fixed':")
    for nm, oo, wk in rows[1:]:
        W("    %-28s %+7.2f   (AG %+7.1f%%, AG:UAG %+7.1f%%)"
          % (nm, (oo["HQ"] - b["HQ"]) - (pl["HQ"] - b["HQ"]),
             100.0 * (oo["AG"] / b["AG"] - 1.0),
             100.0 * (oo["AGUAGeff"] / b["AGUAGeff"] - 1.0)))
    W("")
    W("  the ghrelin columns and the HQ-CT column are close to ORTHOGONAL.")
    W("  octreotide is the extreme case: it fixes the biomarker best, moves")
    W("  the endpoint least, and — once exogenous growth hormone is not there")
    W("  to mask it — collapses the very axis PWS already lacks:")
    W("    octreotide with GH stopped   dIGF1 %+7.1f %%"
      % (100.0 * (rows[7][1]["IGF1"] / b["IGF1"] - 1.0)))
    W("    GH stopped alone             dIGF1 %+7.1f %%"
      % (100.0 * (rows[8][1]["IGF1"] / b["IGF1"] - 1.0)))
    W("  a strictly negative trade that a ghrelin-centred reading of the")
    W("  disease cannot see.")

    # ================= [9] carbetocin optimum ===========================
    W("")
    W("[9] CARBETOCIN — the inverted dose-response is ANALYTIC")
    W(hr())
    E1, K1, E2, K2 = P["ECBMAX"], P["ECB50"], P["EV1AMX"], P["EV1A50"]
    W("  E(C) = E1*C/(K1+C) - E2*C/(K2+C)   OXTR agonism minus V1a effect")
    W("  E1 = %.2f  K1 = %.2f ng/mL      E2 = %.2f  K2 = %.2f ng/mL"
      % (E1, K1, E2, K2))
    num = K1 * math.sqrt(E2 * K2) - K2 * math.sqrt(E1 * K1)
    den = math.sqrt(E1 * K1) - math.sqrt(E2 * K2)
    Cstar = num / den
    Estar = E1 * Cstar / (K1 + Cstar) - E2 * Cstar / (K2 + Cstar)
    W("  dE/dC = 0  =>  C* = [K1*sqrt(E2K2) - K2*sqrt(E1K1)]"
      " / [sqrt(E1K1) - sqrt(E2K2)]")
    W("  C* = %.4f ng/mL,  E(C*) = %.4f,  E(inf) = %.4f"
      % (Cstar, Estar, E1 - E2))
    W("")
    W("  %11s %11s %9s %9s %9s"
      % ("mg TID", "Cavg ng/mL", "E(C)", "E/E(C*)", "dHQ-CT"))
    for dose in (0.8, 1.6, 3.2, 4.8, 6.4, 9.6, 14.4):
        cavg = 3.0 * dose * P["FCB"] / (P["VCB"] * P["KECB"]) * 1000.0
        ee = E1 * cavg / (K1 + cavg) - E2 * cavg / (K2 + cavg)
        s = base12.snapshot()
        s.advance(8 * 7.0, dtd, c12,
                  regimens=[Regimen("ACBD", dose, 1.0, n_per_day=3), ghmt])
        W("  %11.1f %11.3f %9.4f %9.4f %9.2f"
          % (dose, cavg, ee, ee / Estar, s.out(c12)["HQ"] - b["HQ"]))
    c32 = 3.0 * 3.2 * P["FCB"] / (P["VCB"] * P["KECB"]) * 1000.0
    W("")
    W("  3.2 mg TID lands at Cavg %.2f ng/mL, i.e. %.1f x the analytic optimum"
      % (c32, c32 / Cstar))
    W("  C* but still %.0f%% of its EFFECT; 9.6 mg TID is %.1f x C* and only"
      % (100.0 * (E1*c32/(K1+c32) - E2*c32/(K2+c32)) / Estar, 3.0 * c32 / Cstar))
    W("  %.0f%% of the peak effect.  The model therefore reproduces"
      % (100.0 * (E1*(3*c32)/(K1+3*c32) - E2*(3*c32)/(K2+3*c32)) / Estar))
    W("  the CARE-PWS dose ordering (3.2 mg > 9.6 mg) from V1a cross-")
    W("  activation alone, with NO dose-specific parameter, and predicts")
    W("  that the useful dose range is bounded on BOTH sides.")

    # ================= [10] DCCR ========================================
    W("")
    W("[10] DCCR — the K-ATP channel that helps and the one that hurts")
    W(hr())
    W("  %8s %10s %9s %9s %8s %8s"
      % ("mg/kg", "C ug/mL", "E(AgRP)", "E(bcell)", "dHQ-CT", "dHbA1c"))
    for dose in (1.5, 2.5, 3.5, 5.1, 7.0):
        s = base12.snapshot()
        s.advance(13 * 7.0, dtd, c12,
                  regimens=[Regimen("ADZD", dose, 1.0, per_kg=True), ghmt])
        oo = s.out(c12)
        cc = oo["CDZ"]
        W("  %8.1f %10.2f %9.4f %9.4f %8.2f %8.3f"
          % (dose, cc, P["EDZMAX"] * cc / (P["EDZ50"] + cc),
             P["EDZIMX"] * cc / (P["EDZI50"] + cc), oo["HQ"] - b["HQ"],
             oo["HBA1C"] - b["HBA1C"]))
    W("")
    W("  efficacy and hyperglycaemia are the SAME channel in two tissues, so")
    W("  their ratio is set once by EDZ50/EDZI50 = %.3f; no dose separates"
      % (P["EDZ50"] / P["EDZI50"]))
    W("  them and only a tissue-selective opener could.")
    W("")
    W("  baseline-severity interaction (the DESTINY-PWS subgroup result):")
    W("  %-30s %8s %8s %8s" % ("subgroup", "HQ base", "SEEK", "dHQ vs pbo"))
    for lab, cue in (("moderate (CUE 0.28)", 0.28), ("severe (CUE 0.60)", 0.60)):
        cc2 = CT(1.00, cue, 1.0)
        s0 = base12.snapshot(); s0.advance(13 * 7.0, dtd, cc2, regimens=[ghmt])
        s1 = base12.snapshot()
        s1.advance(13 * 7.0, dtd, cc2,
                   regimens=[Regimen("ADZD", 5.1, 1.0, per_kg=True), ghmt])
        W("  %-30s %8.2f %8.3f %8.2f"
          % (lab, s0.out(cc2)["HQ"], s0.out(cc2)["SEEK"],
             s1.out(cc2)["HQ"] - s0.out(cc2)["HQ"]))
    W("")
    W("  the same dose is larger in the severe stratum because there the")
    W("  patient sits nearer the upper branch, so lowering the drive can")
    W("  cross the separatrix.  A subgroup effect can be a bifurcation")
    W("  rather than noise — and if it is, it will replicate.")

    # ================= [11] the latch ===================================
    W("")
    W("[11] THE LATCH — what it takes to leave the food-seeking state")
    W(hr())
    base10 = run_to(les=1.0, age=10.0, gh_from=1.0, dt=dt, **ENV_MGMT)
    latched = base10.snapshot()
    latched.advance(2.0 * 365.25, dt, CT(1.40, 1.00, 0.0), regimens=GHRX())
    lo_ = latched.out(CT(1.40, 1.00, 0.0))
    W("  a patient who spent 10-12 y with free access is at SEEK = %.4f."
      % lo_["SEEK"])
    W("  Standard management is then restored and held for 2 years:")
    W("  %-26s %8s %8s %8s %8s"
      % ("", "SEEK", "HQ-CT", "%fat", "REINF"))
    rec = latched.snapshot()
    prevm = 0
    for mo in (0, 3, 6, 12, 24):
        if mo > prevm:
            rec.advance((mo - prevm) * 30.44, dt, CT(1.00, 0.28, 1.0),
                        regimens=GHRX())
            prevm = mo
        oo = rec.out(CT(1.00, 0.28, 1.0))
        W("  %-26s %8.4f %8.2f %8.1f %8.3f"
          % ("month %d" % mo, oo["SEEK"], oo["HQ"], oo["PBF"], oo["REINF"]))
    W("")
    W("  restoring standard management removes the fat but NOT the state.")
    W("  That is the model's argument for pre-emptive food security: the")
    W("  cheap moment to keep a patient off the upper branch is before they")
    W("  ever reach it, because afterwards standard care no longer suffices.")
    W("")
    W("  what does suffice — 6 months of each intervention, from month 24:")
    W("  %-34s %8s %8s %8s" % ("", "SEEK", "HQ-CT", "branch"))
    for lab, env, extra in (
            ("standard management (continued)", ENV_MGMT, []),
            ("+ carbetocin 1.6 mg TID", ENV_MGMT,
             [Regimen("ACBD", 1.6, 1.0, n_per_day=3)]),
            ("+ carbetocin 3.2 mg TID", ENV_MGMT,
             [Regimen("ACBD", 3.2, 1.0, n_per_day=3)]),
            ("+ DCCR 5.1 mg/kg", ENV_MGMT,
             [Regimen("ADZD", 5.1, 1.0, per_kg=True)]),
            ("+ semaglutide 2.4 mg QW", ENV_MGMT,
             [Regimen("ASGD", 2.4, 7.0)]),
            ("strict environment (CUE 0.10)", ENV_TIGHT, []),
            ("strict environment + carbetocin", ENV_TIGHT,
             [Regimen("ACBD", 3.2, 1.0, n_per_day=3)])):
        cc = CT(env["avail"], env["cue"], env["titr"])
        sx = rec.snapshot()
        sx.advance(182.0, dtd, cc, regimens=GHRX() + extra)
        oo = sx.out(cc)
        W("  %-34s %8.4f %8.2f %8s"
          % (lab, oo["SEEK"], oo["HQ"],
             "LOW" if oo["SEEK"] < 0.25 else "upper"))
    W("")
    W("")
    W("  NOTHING in the panel leaves the upper branch.  Rather than assert")
    W("  that some dose would, solve for what it would take.  Hold the")
    W("  latched patient's drive fixed and ask which satiety level")
    W("  annihilates the upper fixed point:")
    ol = rec.out(CT(1.00, 0.28, 1.0))
    DRVe_l, SELFG_l = ol["DRVi"] * 0.28, ol["SELFG"]

    def n_fp(SAT):
        return len(fixed_points_g(DRVe_l, SAT, SELFG_l))

    lo_s, hi_s = ol["SAT"], 4.0
    for _ in range(60):
        mid = 0.5 * (lo_s + hi_s)
        if n_fp(mid) >= 2:
            lo_s = mid
        else:
            hi_s = mid
    SATcrit = 0.5 * (lo_s + hi_s)
    W("    current satiety            SAT   = %.4f" % ol["SAT"])
    W("    upper state disappears at  SAT_c = %.4f  (%.2f x current)"
      % (SATcrit, SATcrit / ol["SAT"]))
    # translate into the oxytocin-arm gain that would be required
    others = (P["W2"] / ol["x2"] + P["W3"] / ol["x3"] + P["W4"] / ol["x4"]
              + P["W5"] / ol["x5"])
    need = 1.0 / SATcrit - others
    if need <= 0:
        W("    the oxytocin arm ALONE cannot reach SAT_c: even x1 -> infinity")
        W("    leaves 1/SAT >= %.4f, i.e. SAT <= %.4f.  The other four arms")
        W("    have to be repaired as well." % (others, 1.0 / others))
    else:
        gain = (P["W1"] / need) / ol["x1"]
        W("    required oxytocin-arm gain          = %.2f x" % gain)
        W("    best achievable with carbetocin     = %.2f x  (1 + E(C*))"
          % (1.0 + 0.4924))
        W("    shortfall                           = %.2f x" % (gain / 1.4924))
    W("")
    W("  This is the model's concrete drug-development target: it is not a")
    W("  potency problem, it is a STATE problem, and the HQ-CT reductions the")
    W("  panel delivers happen entirely WITHIN the latched state.  A trial")
    W("  powered on mean HQ-CT can succeed while every patient stays latched.")

    # ================= [12] puberty =====================================
    W("")
    W("[12] PUBERTY — MKRN3 loss advances the gate, the hypothalamic lesion")
    W("     lowers its amplitude: early onset, low amplitude, non-progressive")
    W(hr())
    W("  %6s %9s %9s %9s %9s %9s"
      % ("age", "ctrl LH", "PWS LH", "ctrl T", "PWS T", "PWS+Trx"))
    sc = Sim(pars(LES=0.0), init_birth(pars(LES=0.0)), 0.0)
    sp = Sim(pars(LES=1.0), init_birth(pars(LES=1.0)), 0.0)
    st = Sim(pars(LES=1.0), init_birth(pars(LES=1.0)), 0.0)
    prev = 0.0
    cc1 = CT(1.40, 1.0, 0.0)
    for tgt in (8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 16.0, 18.0):
        dur = (tgt - prev) * 365.25
        sc.advance(dur, dt, cc1)
        sp.advance(dur, dt, c12, regimens=GHRX())
        rx = GHRX() + ([Regimen("ATSD", 100.0, 14.0)] if tgt > 14.0 else [])
        st.advance(dur, dt, c12, regimens=rx)
        prev = tgt
        W("  %6.1f %9.2f %9.2f %9.0f %9.0f %9.0f"
          % (tgt, sc.y[IX["LH"]], sp.y[IX["LH"]], sc.out(cc1)["SEX"],
             sp.out(c12)["SEX"], st.out(c12)["SEX"]))
    W("")
    W("  bone age at 18 y: control %.2f, PWS %.2f — the epiphyses stay open"
      % (sc.y[IX["BA"]], sp.y[IX["BA"]]))
    W("  longer in PWS, which is why late growth hormone still buys height,")
    W("  and why sex-steroid replacement shortens that window as it works.")

    # ================= [13] metabolic paradox ===========================
    W("")
    W("[13] THE METABOLIC PARADOX — matched adiposity, unmatched insulin")
    W(hr())
    ob = run_to(les=0.0, age=25.0, dt=dt, EIBF=1.68, KLEPEI=0.05,
                **ENV_CTRL)
    W("%-26s %6s %6s %6s %6s %6s %6s %6s"
      % ("", "BMI", "%fat", "INS", "HOMA", "ADPN", "HbA1c", "GLU"))
    for nm, s2, cc in (("PWS, free access", pfree, CT(1.40, 1.0, 0.0)),
                       ("non-PWS, hyperphagic", ob, CT(1.40, 1.0, 0.0))):
        oo = s2.out(cc)
        W("%-26s %6.1f %6.1f %6.1f %6.2f %6.1f %6.2f %6.0f"
          % (nm, oo["BMI"], oo["PBF"], oo["INS"], oo["HOMA"], oo["ADPN"],
             oo["HBA1C"], oo["GLU"]))
    W("")
    W("  the lesion that starves the satiety arms also throttles proinsulin")
    W("  processing, so PWS reaches a given adiposity with LESS insulin and")
    W("  MORE adiponectin.  Same fat mass, different metabolic phenotype,")
    W("  from one parameter rather than two.")

    # ================= [14] behaviour trade-off =========================
    W("")
    W("[14] THE MANAGEMENT TRADE-OFF, MADE EXPLICIT")
    W(hr())
    W("  %-30s %8s %8s %8s %8s"
      % ("environment", "%fat", "HQ-CT", "BEH", "frustr."))
    for lab, env in (("free access", ENV_FREE),
                     ("normal portions, no titration",
                      dict(avail=1.00, cue=0.55, titr=0.0)),
                     ("standard PWS management", ENV_MGMT),
                     ("strict environment", ENV_TIGHT)):
        s2 = run_to(les=1.0, age=14.0, gh_from=1.0, dt=dt, **env)
        cc = CT(env["avail"], env["cue"], env["titr"])
        oo = s2.out(cc)
        W("  %-30s %8.1f %8.2f %8.1f %8.3f"
          % (lab, oo["PBF"], oo["HQ"], oo["BEH"], oo["FRUST"]))
    W("")
    W("  restriction buys adiposity and costs behaviour, because the")
    W("  unrealized drive IS the frustration term.  Carbetocin is the only")
    W("  arm in the panel that lowers both at once, since oxytocin enters")
    W("  the drive and the behaviour term through the same receptor.")

    # ================= [15] numerical integrity =========================
    W("")
    W("[15] NUMERICAL INTEGRITY")
    W(hr())
    s = run_to(les=1.0, age=12.0, gh_from=1.0, dt=dt, **ENV_MGMT)
    s2 = run_to(les=1.0, age=12.0, gh_from=1.0,
                dt=dt / 2.0, **ENV_MGMT)
    o1, o2 = s.out(c12), s2.out(c12)
    W("  step-halving check, PWS on GH, birth -> 12 y")
    W("  %-8s %13s %13s %11s" % ("output", "dt=%.4f" % dt,
                                 "dt=%.4f" % (dt / 2.0), "rel. diff"))
    worst = 0.0
    for k in ("PBF", "HQ", "IGF1", "AG", "AHI", "HT", "SAT", "SEEK", "LFM",
              "HBA1C"):
        r = abs(o1[k] - o2[k]) / max(1e-9, abs(o2[k]))
        worst = max(worst, r)
        W("  %-8s %13.5f %13.5f %11.2e" % (k, o1[k], o2[k], r))
    W("  worst relative drift on halving dt: %.2e" % worst)
    W("")
    W("=" * 78)
    W("end of reference output")
    W("=" * 78)
    print("\n".join(out))


if __name__ == "__main__":
    main()
