#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cio_reference_model.py
=======================
Independent, self-contained re-implementation of the cisplatin-induced
ototoxicity (CIO) QSP model shipped here as `cio_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
No R/mrgsolve runtime is available in the build environment, so every
structural claim made in README.md is derived here instead, from the SAME
equations, and the numbers printed by this script are the numbers quoted in
README.md.  The R file is a transliteration of what follows: same parameter
names, same state names, same ordering, same closed-form read-outs.

MODEL AT A GLANCE
-----------------
73 ODEs:
  * 25 systemic / cochlear-global states — cisplatin PK, kidney and GFR,
    tumour platinum and Pt-DNA adducts, stria vascularis, perilymph,
    endocochlear potential, cochlear inflammation, sodium thiosulfate PK
    (systemic and intratympanic), N-acetylcysteine
  * 8 tonotopic bands x 6 states each — labile (redox-active) platinum,
    bound (inert but persistent) platinum, glutathione, outer hair cells,
    inner hair cells, spiral ganglion neurons

The eight bands sit at their Greenwood positions for 0.25, 0.5, 1, 2, 4, 6,
8 and 12.5 kHz, so "the audiogram" is read off a spatial gradient rather
than fitted eight times.

THE FOUR STRUCTURAL CHOICES THAT DO ALL THE WORK
------------------------------------------------
1. Cochlear platinum is a STOCK, not a concentration.  PT_j has an efflux
   half-life of years, so the damaging variable is an integral that does not
   reset between cycles and does not stop when treatment does.

2. The cochlea is reached through a CASCADE, not directly: free plasma
   platinum -> stria vascularis -> perilymph -> hair cell.  Free plasma
   platinum has a ~23 min half-life; perilymph platinum peaks ~8 h later.
   Everything about the 6-hour thiosulfate rule is that one lag, and the
   rescue is possible only because thiosulfate (a small diffusible anion)
   does NOT need the transporter-mediated strial route that platinum does.

3. Glutathione is consumed at a SATURATING rate by the oxidant flux, so each
   band has a critical platinum load PT* = KSYN*GMAX_j/(KCON*KROS) above
   which its reserve collapses.  Ototoxicity is therefore a threshold that
   sweeps apex-ward through the cochlea as the cumulative dose rises — not a
   dose-response curve fitted per frequency.

4. Tonotopic vulnerability is that threshold divided by the local dose rate,
   PT*_j / uptake_j, i.e. the PRODUCT of a reserve gradient and an uptake
   gradient.  Which of the two carries the audiogram slope is measured in
   Part 6, not assumed.

Run:  python3 cio_reference_model.py        (also writes cio_reference_output.txt)
"""

import math
import os
import sys
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# 0. Tonotopy — Greenwood positions of the eight audiometric test frequencies
# ----------------------------------------------------------------------------
# Human Greenwood map: f(Hz) = 165.4 * (10^(2.1*x) - 0.88), x = 0 at the apex
# and x = 1 at the base.  Inverting gives each band's normalised position.
FREQS_HZ = np.array([250.0, 500.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 12500.0])
FREQ_LABEL = ["0.25k", "0.5k", "1k", "2k", "4k", "6k", "8k", "12.5k"]
XPOS = np.log10(FREQS_HZ / 165.4 + 0.88) / 2.1
NB = len(FREQS_HZ)
_RTOL_OVERRIDE = None
I_SPEECH = [1, 2, 3, 4]          # 0.5/1/2/4 kHz — the speech-frequency average

# ----------------------------------------------------------------------------
# 1. Parameters
# ----------------------------------------------------------------------------
P = OrderedDict(
    # ---- patient -----------------------------------------------------------
    BSA=1.73,          # m2
    GFR0=100.0,        # mL/min
    AGE=8.0,           # years — scales cochlear antioxidant reserve
    AGEK=4.0,          # years to half-maximal maturation of that reserve
    AGEMIN=0.62,       # reserve floor in the neonate (fraction of adult)

    # ---- cisplatin PK, free (ultrafilterable) platinum ----------------------
    # Free platinum leaving plasma is excreted, protein-bound or already
    # covalently bound in tissue; none of it returns as reactive drug, so the
    # peripheral compartment is a terminal sink (K21 = 0).
    VC=20.0,           # L
    KBIND=1.30,        # 1/h, irreversible binding to plasma protein
    CLRSL=0.054,       # L/h per mL/min GFR  ->  CLR = CLRSL * GFR
    K12=0.30, K21=0.0,  # 1/h, free Pt -> tissue (irreversible)
    KELB=0.0060,       # 1/h, elimination of protein-bound Pt (t1/2 ~ 4.8 d)
    MWCIS=300.05,      # g/mol

    # ---- kidney ------------------------------------------------------------
    KKID=0.55,         # 1/h, OCT2-mediated cortical uptake (tracer, no plasma sink)
    KKOUT=0.030,       # 1/h, cortical platinum efflux
    KTUB=0.0200,       # 1/h, tubular injury formation
    KT50=9.0,          # cortical Pt at half-maximal injury
    KREPT=0.0060,      # 1/h, tubular repair
    KLOSS=0.0120,      # 1/h, reversible GFR loss per unit injury
    KREC=0.0060,       # 1/h, recovery of reversible GFR loss
    KPERM=0.00016,     # 1/h, permanent nephron loss per unit injury

    # ---- tumour ------------------------------------------------------------
    QT=0.90,           # 1/h, perfusion-limited exchange at TUMPERF = 1
    TUMPERF=1.0,       # 1.0 well-perfused localised; 0.15 disseminated
    KPTUM=1.0,         # tissue:plasma partition
    KADF=0.55,         # 1/h, Pt-DNA adduct formation
    KREPAIR=0.018,     # 1/h, nucleotide-excision repair of adducts
    KKILL=0.045,       # log-kill per adduct-hour
    FADCARB=0.25,      # adduct-forming potency of carboplatin relative to cisplatin
    FCARB=0.060,       # cochlear uptake of carboplatin relative to cisplatin, per mole

    # ---- cochlear route: plasma -> stria -> perilymph -----------------------
    KINBL=0.085,       # 1/h, free plasma Pt -> stria vascularis (CTR1/OCT2)
    KSTRO=0.050,       # 1/h, strial platinum efflux
    KSTP=0.150,        # 1/h, stria -> perilymph
    KOUTPE=0.120,      # 1/h, perilymph platinum clearance
    FURO=0.0,          # 0/1 loop diuretic co-administration
    FUROF=1.45,         # fold increase in KINBL on furosemide

    # ---- cochlear route: perilymph -> hair cell ----------------------------
    KUPT=0.025,        # 1/h, perilymph -> organ of Corti
    BUPT=0.95,         # base-ward uptake gradient (log-slope over x)
    # Two intracellular platinum pools.  Only the LABILE pool is redox-active;
    # the BOUND pool is what a mass-spectrometry study measures decades later.
    TLAB_D=60.0,      # d, half-life of the labile, redox-active pool
    FBND=0.75,         # fraction of labile platinum that ends up covalently bound
    TRET_D=730.0,      # d, retention half-life of the bound (inert) pool

    # ---- redox -------------------------------------------------------------
    GSH0=1.00,         # apical glutathione capacity (normalised)
    BGSH=1.35,         # base-ward decline in reserve (log-slope over x)
    KSYN=0.045,        # 1/h, glutathione resynthesis
    KCON=0.02720,      # 1/h per unit oxidant flux, saturating GSH consumption
    KMG=0.050,         # Michaelis constant of that consumption
    KROS=1.00,         # oxidant flux per unit retained platinum
    RBASE=0.010,       # constitutive oxidant flux
    KSCV=1.00,         # scavenging strength of glutathione
    KBOX=0.030,        # non-glutathione floor of antioxidant defence
    NOISE=0.0,         # 0/1 concurrent high-level noise exposure
    FNOISE=0.16,       # oxidant-flux multiplier from noise
    AMGLY=0.0,         # 0/1 concurrent aminoglycoside
    FAG=0.26,          # oxidant-flux multiplier from aminoglycoside

    # ---- cell death --------------------------------------------------------
    KAP=0.00220,       # 1/h, maximal hair-cell death rate
    OXC=9.0,           # oxidative stress at half-maximal death rate
    HILL=1.70,         # steepness of the death function
    FIHC=0.30,         # inner hair cells are more resistant than outer
    FSGN=0.10,         # direct spiral-ganglion sensitivity
    KDEAFF=0.0016,     # 1/h, secondary ganglion loss after inner-hair-cell loss
    GINFL=0.85,        # inflammation multiplier on the death rate

    # ---- stria vascularis and endocochlear potential -----------------------
    KSV=0.0060,        # 1/h, strial injury rate
    KSVREP=0.00008,    # 1/h, strial repair
    OXCSV=0.80,         # strial stress at half-maximal injury
    EP0=85.0,          # mV, normal endocochlear potential
    KEP=0.25,          # 1/h, EP equilibration with strial function
    GEP=0.55,          # exponent linking EP to cochlear amplifier gain

    # ---- cochlear inflammation --------------------------------------------
    KINF=0.010,        # 1/h, induction by the damage signal
    KINFOUT=0.030,     # 1/h, resolution

    # ---- audiometry --------------------------------------------------------
    TSO=48.0,          # dB, shift at complete loss of the cochlear amplifier
    GTS=2.20,          # curvature of amplifier loss -> dB
    TSI=42.0,          # dB, further shift from inner-hair-cell (dead-region) loss
    TSMAX=110.0,       # dB, audiometer ceiling

    # ---- sodium thiosulfate ------------------------------------------------
    # Thiosulfate is a small hydrophilic anion and reaches perilymph directly
    # rather than through the transporter-mediated strial route platinum uses.
    VCS=20.0,          # L
    KELS=1.00,         # 1/h (t1/2 ~ 42 min)
    K12S=0.55, K21S=0.40,
    KINBLS=0.010,      # 1/h, plasma -> perilymph
    KOUTPES=0.500,     # 1/h, perilymph clearance of thiosulfate
    K2STS=0.006,        # 1/(uM*h), 2nd-order platinum + thiosulfate inactivation
    K2LOC=0.002,       # 1/uM, local block of hair-cell platinum uptake
    BSTS=0.45,         # maximal fractional boost of glutathione capacity
    KSTS=60.0,         # uM, half-maximal concentration for that boost
    KITA=0.55,         # 1/h, round-window absorption of an intratympanic dose
    FRW=0.005,         # fraction of a middle-ear dose that crosses the round window
    VPERI=1.6e-4,      # L, human perilymph volume
    KITOUT=1.20,       # 1/h, clearance of intratympanically-derived thiosulfate
    BIT=2.20,          # base-ward gradient of intratympanic delivery
    MWSTS=248.2,       # g/mol, sodium thiosulfate pentahydrate

    # ---- N-acetylcysteine --------------------------------------------------
    VCN=30.0,
    KELN=0.35,
    KINBLN=0.030,
    KOUTPEN=0.30,
    BNAC=0.40,
    KNAC=40.0,
    MWNAC=163.2,
)

GLOBAL_STATES = [
    "CISC", "CISB", "CISP", "KIDPT", "TUBI", "GFRL", "GPERM",
    "TUMPT", "TUMAD", "TUMLK", "STRPT", "PERI", "SVF", "EP", "INFL",
    "STSC", "STSP", "STSPE", "STSIT", "STSTU", "ITD", "NACC", "NACPE",
    "CUMD", "PTAUC",
]
NG = len(GLOBAL_STATES)
IDX = {n: i for i, n in enumerate(GLOBAL_STATES)}
I_PT, I_PTB, I_GSH, I_OHC, I_IHC, I_SGN = (NG, NG + NB, NG + 2 * NB, NG + 3 * NB,
                                           NG + 4 * NB, NG + 5 * NB)
NSTATE = NG + 6 * NB            # 25 + 48 = 73


# ----------------------------------------------------------------------------
# 2. Right-hand side
# ----------------------------------------------------------------------------
def gradients(p):
    upt = np.exp(p["BUPT"] * (XPOS - 0.5))
    agef = p["AGEMIN"] + (1.0 - p["AGEMIN"]) * p["AGE"] / (p["AGEK"] + p["AGE"])
    agef = agef / (p["AGEMIN"] + (1.0 - p["AGEMIN"]) * 25.0 / (p["AGEK"] + 25.0))
    gmax = p["GSH0"] * np.exp(-p["BGSH"] * XPOS) * agef
    itg = np.exp(p["BIT"] * (XPOS - 0.5))
    return upt, gmax, itg


def y0(p):
    y = np.zeros(NSTATE)
    _, gmax, _ = gradients(p)
    y[IDX["SVF"]] = 1.0
    y[IDX["EP"]] = p["EP0"]
    y[I_GSH:I_GSH + NB] = gmax
    y[I_OHC:I_OHC + NB] = 1.0
    y[I_IHC:I_IHC + NB] = 1.0
    y[I_SGN:I_SGN + NB] = 1.0
    return y


def rhs(t, y, p, upt, gmax, itg):
    d = np.zeros_like(y)
    g = IDX
    pos = lambda v: np.maximum(v, 0.0)

    CISC, CISP = pos(y[g["CISC"]]), pos(y[g["CISP"]])
    CISB = pos(y[g["CISB"]])
    KIDPT, TUBI = pos(y[g["KIDPT"]]), min(max(y[g["TUBI"]], 0.0), 1.0)
    GFRL, GPERM = pos(y[g["GFRL"]]), pos(y[g["GPERM"]])
    TUMPT, TUMAD = pos(y[g["TUMPT"]]), pos(y[g["TUMAD"]])
    STRPT, PERI = pos(y[g["STRPT"]]), pos(y[g["PERI"]])
    SVF = min(max(y[g["SVF"]], 0.0), 1.0)
    EP, INFL = pos(y[g["EP"]]), pos(y[g["INFL"]])
    STSC, STSP = pos(y[g["STSC"]]), pos(y[g["STSP"]])
    STSPE, STSIT = pos(y[g["STSPE"]]), pos(y[g["STSIT"]])
    STSTU, ITD = pos(y[g["STSTU"]]), pos(y[g["ITD"]])
    NACC, NACPE = pos(y[g["NACC"]]), pos(y[g["NACPE"]])

    PT = pos(y[I_PT:I_PT + NB])
    PTB = pos(y[I_PTB:I_PTB + NB])
    GSH = pos(y[I_GSH:I_GSH + NB])
    OHC = np.clip(y[I_OHC:I_OHC + NB], 0.0, 1.0)
    IHC = np.clip(y[I_IHC:I_IHC + NB], 0.0, 1.0)
    SGN = np.clip(y[I_SGN:I_SGN + NB], 0.0, 1.0)

    # -- renal function ------------------------------------------------------
    floss = min(0.85, GFRL + GPERM)
    GFR = p["GFR0"] * (1.0 - floss)

    # -- cisplatin PK --------------------------------------------------------
    d[g["CISC"]] = (-p["KBIND"] * CISC - (p["CLRSL"] * GFR / p["VC"]) * CISC
                    - p["K12"] * CISC + p["K21"] * CISP
                    - p["K2STS"] * CISC * STSC)
    d[g["CISB"]] = p["KBIND"] * CISC - p["KELB"] * CISB
    d[g["CISP"]] = p["K12"] * CISC - p["K21"] * CISP

    # -- kidney and nephrotoxicity ------------------------------------------
    d[g["KIDPT"]] = p["KKID"] * CISC - p["KKOUT"] * KIDPT
    d[g["TUBI"]] = (p["KTUB"] * KIDPT / (p["KT50"] + KIDPT) * (1.0 - TUBI)
                    - p["KREPT"] * TUBI)
    d[g["GFRL"]] = p["KLOSS"] * TUBI - p["KREC"] * GFRL
    d[g["GPERM"]] = p["KPERM"] * TUBI

    # -- tumour --------------------------------------------------------------
    qt = p["QT"] * p["TUMPERF"]
    d[g["TUMPT"]] = (qt * (CISC - TUMPT / p["KPTUM"]) - p["KADF"] * TUMPT
                     - p["K2STS"] * TUMPT * STSTU)
    d[g["TUMAD"]] = p["KADF"] * p["FADC"] * TUMPT - p["KREPAIR"] * TUMAD
    d[g["TUMLK"]] = p["KKILL"] * TUMAD
    d[g["STSTU"]] = qt * (STSC - STSTU) - p["K2STS"] * TUMPT * STSTU

    # -- cochlear platinum: plasma -> stria -> perilymph -> hair cell --------
    kin_bl = p["KINBL"] * (1.0 + p["FURO"] * (p["FUROF"] - 1.0))
    # Systemic thiosulfate arrives through the blood and is well mixed, so it
    # acts on the whole perilymph pool.  An intratympanic bolus does not: it
    # enters at the round window and is cleared from scala tympani far faster
    # than it can equilibrate along the duct, so it is a BAND-LOCAL term,
    # weighted by the base-ward delivery gradient itg.
    STSLOC = STSPE + STSIT * itg
    uptake = p["KUPT"] * upt * PERI / (1.0 + p["K2LOC"] * STSLOC)
    # NB: thiosulfate reaches plasma and perilymph, NOT the intracellular
    # strial depot — so the stria keeps refilling perilymph after the rescue
    # has washed out, which is what caps the achievable protection.
    d[g["STRPT"]] = kin_bl * CISC - (p["KSTRO"] + p["KSTP"]) * STRPT
    d[g["PERI"]] = (p["KSTP"] * STRPT - p["KOUTPE"] * PERI
                    - p["K2STS"] * PERI * STSPE - uptake.mean())
    ktot = math.log(2.0) / (p["TLAB_D"] * 24.0)
    kbnd, klab = p["FBND"] * ktot, (1.0 - p["FBND"]) * ktot
    kout_c = math.log(2.0) / (p["TRET_D"] * 24.0)
    d[I_PT:I_PT + NB] = uptake - (kbnd + klab) * PT
    d[I_PTB:I_PTB + NB] = kbnd * PT - kout_c * PTB

    # -- redox: saturating glutathione consumption gives a per-band threshold -
    env = 1.0 + p["FNOISE"] * p["NOISE"] + p["FAG"] * p["AMGLY"]
    flux = (p["KROS"] * PT + p["RBASE"]) * env
    boost = (1.0 + p["BSTS"] * STSLOC / (p["KSTS"] + STSLOC)
             + p["BNAC"] * NACPE / (p["KNAC"] + NACPE))
    d[I_GSH:I_GSH + NB] = (p["KSYN"] * (gmax * boost - GSH)
                           - p["KCON"] * flux * GSH / (p["KMG"] + GSH))
    OX = flux / (p["KSCV"] * GSH + p["KBOX"])

    # -- cell death ----------------------------------------------------------
    hz = p["KAP"] * OX ** p["HILL"] / (p["OXC"] ** p["HILL"] + OX ** p["HILL"])
    hz = hz * (1.0 + p["GINFL"] * INFL)
    d[I_OHC:I_OHC + NB] = -hz * OHC
    d[I_IHC:I_IHC + NB] = -p["FIHC"] * hz * IHC
    d[I_SGN:I_SGN + NB] = -p["FSGN"] * hz * SGN - p["KDEAFF"] * (1.0 - IHC) * SGN

    # -- stria vascularis and endocochlear potential -------------------------
    ox_sv = p["KROS"] * STRPT * env / (p["KSCV"] * p["GSH0"] + p["KBOX"])
    d[g["SVF"]] = (-p["KSV"] * ox_sv ** p["HILL"]
                   / (p["OXCSV"] ** p["HILL"] + ox_sv ** p["HILL"]) * SVF
                   + p["KSVREP"] * (1.0 - SVF))
    d[g["EP"]] = p["KEP"] * (p["EP0"] * SVF - EP)

    # -- cochlear inflammation ----------------------------------------------
    d[g["INFL"]] = p["KINF"] * (hz.mean() / p["KAP"]) * (1.0 - INFL) - p["KINFOUT"] * INFL

    # -- sodium thiosulfate --------------------------------------------------
    d[g["STSC"]] = (-p["KELS"] * STSC - p["K12S"] * STSC + p["K21S"] * STSP
                    - p["K2STS"] * CISC * STSC)
    d[g["STSP"]] = p["K12S"] * STSC - p["K21S"] * STSP
    d[g["STSPE"]] = (p["KINBLS"] * STSC - p["KOUTPES"] * STSPE
                     - p["K2STS"] * PERI * STSPE)
    d[g["ITD"]] = -p["KITA"] * ITD
    d[g["STSIT"]] = (p["KITA"] * ITD * p["FRW"] / p["VPERI"] - p["KITOUT"] * STSIT
                     - p["K2STS"] * PERI * STSIT)

    # -- N-acetylcysteine ----------------------------------------------------
    d[g["NACC"]] = -p["KELN"] * NACC - p["KINBLN"] * NACC
    d[g["NACPE"]] = p["KINBLN"] * NACC - p["KOUTPEN"] * NACPE

    # -- trackers ------------------------------------------------------------
    d[g["PTAUC"]] = PT.mean()
    return d


# ----------------------------------------------------------------------------
# 3. Read-outs
# ----------------------------------------------------------------------------
def threshold_shift(y, p):
    OHC = np.clip(y[I_OHC:I_OHC + NB], 0.0, 1.0)
    IHC = np.clip(y[I_IHC:I_IHC + NB], 0.0, 1.0)
    EP = max(y[IDX["EP"]], 0.0)
    amp = OHC * (EP / p["EP0"]) ** p["GEP"]
    ts = p["TSO"] * (1.0 - amp) ** p["GTS"] + p["TSI"] * (1.0 - IHC)
    return np.minimum(ts, p["TSMAX"])


def brock_grade(ts):
    """Brock: 40 dB criterion at 8 / 4 / 2 / 1 kHz, marching apex-ward.
    Pre-treatment hearing is assumed normal, so shift == threshold."""
    if ts[2] >= 40: return 4
    if ts[3] >= 40: return 3
    if ts[4] >= 40: return 2
    if ts[6] >= 40: return 1
    return 0


def siop_grade(ts):
    """SIOP Boston Ototoxicity Scale — 20 dB criterion."""
    if ts[2] > 20: return 4
    if ts[3] > 20: return 3
    if ts[4] > 20: return 2
    if max(ts[5], ts[6], ts[7]) > 20: return 1
    return 0


def asha_positive(ts):
    """ASHA: >=20 dB shift at one frequency, or >=10 dB at two adjacent."""
    return bool(np.any(ts >= 20.0) or np.any((ts[:-1] >= 10.0) & (ts[1:] >= 10.0)))


def ctcae_grade(ts):
    """CTCAE v5.0 adult hearing-loss grade, from the shift at 2-8 kHz."""
    hi = max(ts[3], ts[4], ts[5], ts[6])
    spe = float(np.mean(ts[I_SPEECH]))
    if spe >= 40: return 4
    if spe >= 25: return 3
    if hi >= 25: return 2
    if hi >= 15: return 1
    return 0


def pta_speech(ts):
    return float(np.mean(ts[I_SPEECH]))


def snr_loss(y):
    """Speech-in-noise penalty (dB SNR) from inner-hair-cell / ganglion loss —
    the component a pure-tone audiogram cannot see."""
    IHC = np.clip(y[I_IHC:I_IHC + NB], 0, 1)
    SGN = np.clip(y[I_SGN:I_SGN + NB], 0, 1)
    return float(12.0 * (1.0 - np.mean(IHC[I_SPEECH] * SGN[I_SPEECH])))


def tinnitus_prob(y, p):
    ts = threshold_shift(y, p)
    z = -2.2 + 0.055 * float(np.mean(ts[4:]))
    return 1.0 / (1.0 + math.exp(-z))


# ----------------------------------------------------------------------------
# 4. Dosing engine
# ----------------------------------------------------------------------------
def simulate(p, cycles=6, cycle_h=21 * 24.0, dose_mgm2=100.0, ndays=1, tinf=1.0,
             sts=None, sts_gm2=20.0, sts_delay=6.0,
             it_sts=None, it_umol=161.0, nac=None, nac_mg=0.0,
             follow_h=0.0, dense=None, carbo=False):
    """`sts`, `it_sts`, `nac` are lists of 0-based cycle indices at which the
    protectant is given; None means never."""
    p = dict(p)
    p.setdefault("FADC", 1.0)
    if carbo:
        # Carboplatin: one leaving group is a stable dicarboxylate, so per mole
        # it aquates far more slowly — less cochlear uptake, fewer adducts.
        p.update(KINBL=p["KINBL"] * p["FCARB"], KKID=p["KKID"] * 0.25,
                 MWCIS=371.25, FADC=p["FADCARB"])

    y = y0(p)
    upt, gmax, itg = gradients(p)
    dose_umol = dose_mgm2 * p["BSA"] / p["MWCIS"] * 1e3 / ndays
    sts_umol = sts_gm2 * p["BSA"] / p["MWSTS"] * 1e3     # mmol -> umol below
    sts_umol = sts_gm2 * p["BSA"] / p["MWSTS"] * 1e6

    total_h = cycles * cycle_h + follow_h
    events = []
    for c in range(cycles):
        t0c = c * cycle_h
        for dday in range(ndays):
            ts_ = t0c + dday * 24.0
            events.append((ts_, ts_ + tinf, dose_umol / p["VC"] / tinf, "CISC"))
        anchor = t0c + (ndays - 1) * 24.0 + tinf
        if sts is not None and c in sts:
            events.append((anchor + sts_delay, anchor + sts_delay + 0.25,
                           sts_umol / p["VCS"] / 0.25, "STSC"))
        if it_sts is not None and c in it_sts:
            events.append((anchor + sts_delay, anchor + sts_delay + 0.05,
                           it_umol / 0.05, "ITD"))
        if nac is not None and c in nac:
            events.append((anchor + sts_delay, anchor + sts_delay + 0.5,
                           nac_mg / p["MWNAC"] * 1e3 / p["VCN"] / 0.5, "NACC"))

    bps = sorted(set([0.0, total_h] + [e[0] for e in events] + [e[1] for e in events]
                     + [c * cycle_h for c in range(cycles + 1)]))
    bps = [b for b in bps if 0.0 <= b <= total_h + 1e-9]

    rec_t, rec_y = [0.0], [y.copy()]
    for a, b in zip(bps[:-1], bps[1:]):
        if b - a < 1e-12:
            continue
        rates = np.zeros(NSTATE)
        for (e0, e1, r, st) in events:
            if e0 <= a + 1e-9 and b <= e1 + 1e-9:
                rates[IDX[st]] += r

        def f(t, yy, _r=rates):
            return rhs(t, yy, p, upt, gmax, itg) + _r

        # NB: t_eval must always contain the interval end, otherwise solve_ivp
        # returns its last REQUESTED point as sol.y[:, -1] and the state carried
        # into the next interval is silently truncated.
        teval = None
        if dense is not None:
            teval = np.array(sorted(set([tt for tt in dense if a < tt < b] + [b])))
        sol = solve_ivp(f, (a, b), y, method="LSODA",
                        rtol=(_RTOL_OVERRIDE or 1e-7), atol=1e-10,
                        t_eval=teval, max_step=(b - a))
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1].copy()
        if teval is not None:
            for k in range(sol.t.size):
                rec_t.append(sol.t[k]); rec_y.append(sol.y[:, k].copy())
        else:
            rec_t.append(b); rec_y.append(y.copy())
        assert abs(sol.t[-1] - b) < 1e-8, "solver did not reach interval end"
    return dict(p=p, y=y, t=np.array(rec_t), Y=np.array(rec_y).T,
                cum_dose=cycles * dose_mgm2 * ndays if ndays > 1 else cycles * dose_mgm2)


def report(res, name="x"):
    p, y = res["p"], res["y"]
    ts = threshold_shift(y, p)
    return dict(
        name=name, ts=ts, pta=pta_speech(ts), brock=brock_grade(ts),
        siop=siop_grade(ts), asha=asha_positive(ts), ctcae=ctcae_grade(ts),
        snr=snr_loss(y), tinn=tinnitus_prob(y, p),
        ohc=np.clip(y[I_OHC:I_OHC + NB], 0, 1),
        ihc=np.clip(y[I_IHC:I_IHC + NB], 0, 1),
        gsh=y[I_GSH:I_GSH + NB], pt=y[I_PT:I_PT + NB],
        ptb=y[I_PTB:I_PTB + NB], ep=y[IDX["EP"]],
        gfr=p["GFR0"] * (1 - min(0.85, y[IDX["GFRL"]] + y[IDX["GPERM"]])),
        logkill=y[IDX["TUMLK"]], adduct=y[IDX["TUMAD"]], y=y)


# ----------------------------------------------------------------------------
# 5. Output helpers
# ----------------------------------------------------------------------------
OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def bl(tag, v, fmt="{:7.1f}"):
    return tag.ljust(32) + " ".join(fmt.format(x) for x in v)


def hdr():
    return " " * 32 + " ".join("{:>7s}".format(f) for f in FREQ_LABEL)


# ============================================================================
def main():
    say("=" * 104)
    say("CISPLATIN-INDUCED OTOTOXICITY (CIO) QSP MODEL — reference implementation")
    say("73 ODEs = 25 systemic/cochlear-global + 8 tonotopic bands x 6 states")
    say("=" * 104)
    say()
    upt, gmax, itg = gradients(P)
    say(hdr() + "     unit")
    say(bl("frequency", FREQS_HZ / 1000.0, "{:7.2f}") + "     kHz")
    say(bl("Greenwood x (0 apex, 1 base)", XPOS, "{:7.3f}"))
    say(bl("platinum uptake gradient", upt, "{:7.3f}"))
    say(bl("glutathione reserve GMAX", gmax, "{:7.3f}"))
    say(bl("intratympanic delivery", itg, "{:7.3f}"))
    ptstar = P["KSYN"] * gmax / (P["KCON"] * P["KROS"])
    say(bl("critical platinum load PT*", ptstar, "{:7.3f}"))
    say(bl("PT*/uptake (relative)", ptstar / upt / (ptstar / upt)[6], "{:7.2f}"))
    say()
    say("  base/apex uptake ratio  = {:.2f}".format(upt[-1] / upt[0]))
    say("  apex/base reserve ratio = {:.2f}".format(gmax[0] / gmax[-1]))
    say("  apex/base vulnerability = {:.2f}  (the product of the two)".format(
        (gmax[0] / upt[0]) / (gmax[-1] / upt[-1])))
    say()

    # ---------------------------------------------------------------- PART 1
    say("-" * 104)
    say("PART 1 — PHARMACOKINETICS: THE LAG THE WHOLE MODEL RESTS ON")
    say("-" * 104)
    dense = np.unique(np.concatenate([np.linspace(0, 72, 1441),
                                      np.linspace(72, 21 * 24, 300)]))
    r1 = simulate(P, cycles=1, dose_mgm2=100.0, dense=dense)
    t, Y = r1["t"], r1["Y"]
    cisc, stro, peri = Y[IDX["CISC"]], Y[IDX["STRPT"]], Y[IDX["PERI"]]
    say("  free plasma platinum   Cmax {:7.2f} uM  at {:5.2f} h".format(
        cisc.max(), t[cisc.argmax()]))
    m = (t > 1.2) & (t < 3.0)
    say("  free plasma platinum   t1/2 {:7.1f} min".format(
        math.log(2) / (-np.polyfit(t[m], np.log(cisc[m]), 1)[0]) * 60))
    say("  stria vascularis Pt    Cmax {:7.3f}     at {:5.2f} h".format(
        stro.max(), t[stro.argmax()]))
    say("  perilymph platinum     Cmax {:7.3f} uM  at {:5.2f} h  "
        "(peri:plasma Cmax {:.4f})".format(peri.max(), t[peri.argmax()],
                                           peri.max() / cisc.max()))
    m = (t > 30) & (t < 90)
    say("  perilymph platinum     t1/2 {:7.1f} h".format(
        math.log(2) / (-np.polyfit(t[m], np.log(np.maximum(peri[m], 1e-14)), 1)[0])))
    say()

    def ahead(tv, cv, thr):
        tot = np.trapezoid(cv, tv)
        m = tv >= thr
        return float(np.trapezoid(cv[m], tv[m]) / tot)

    # tumour adduct formation still ahead at time thr, both perfusion states
    ad_ahead = {}
    for perf in (1.0, 0.15):
        rr = simulate(dict(P, TUMPERF=perf), cycles=1, dose_mgm2=100.0, dense=dense)
        form = P["KADF"] * rr["Y"][IDX["TUMPT"]]
        ad_ahead[perf] = [ahead(rr["t"], form, th) for th in
                          [0, 1, 2, 4, 6, 8, 12, 24]]

    say("  WHAT IS STILL AHEAD OF YOU AT TIME t (fraction of each integral, %)")
    say("  {:>6s} {:>16s} {:>16s} {:>20s} {:>20s}".format(
        "t (h)", "free plasma Pt", "cochlear uptake", "tumour adducts",
        "tumour adducts"))
    say("  {:>6s} {:>16s} {:>16s} {:>20s} {:>20s}".format(
        "", "", "", "well-perfused", "poorly perfused"))
    for i, thr in enumerate([0, 1, 2, 4, 6, 8, 12, 24]):
        say("  {:6.1f} {:16.5f} {:16.1f} {:20.2f} {:20.2f}".format(
            thr, ahead(t, cisc, thr) * 100, ahead(t, peri, thr) * 100,
            ad_ahead[1.0][i] * 100, ad_ahead[0.15][i] * 100))
    say()
    say("  Read the 6 h row.  A systemic rescue given then is downstream of all but")
    say("  {:.4f}% of the free plasma platinum and {:.2f}% of the adduct formation in a".format(
        ahead(t, cisc, 6.0) * 100, ad_ahead[1.0][4] * 100))
    say("  well-perfused tumour, but still upstream of {:.0f}% of the platinum the".format(
        ahead(t, peri, 6.0) * 100))
    say("  cochlea has yet to take up.  In a poorly perfused tumour the adduct")
    say("  figure is {:.1f}%, i.e. {:.0f}x larger — the same protocol is no longer free.".format(
        ad_ahead[0.15][4] * 100, ad_ahead[0.15][4] / ad_ahead[1.0][4]))
    say()

    # ---------------------------------------------------------------- PART 2
    say("-" * 104)
    say("PART 2 — SCENARIO LIBRARY (20 scenarios)")
    say("-" * 104)
    S = OrderedDict()
    S["S01 adult 100 mg/m2 q3w x6"] = (dict(P), dict(cycles=6, dose_mgm2=100.0))
    S["S02 paediatric 80 mg/m2 q3w x6"] = (dict(P, AGE=3.0), dict(cycles=6, dose_mgm2=80.0))
    S["S03 S02 + STS 20 g/m2 at 6 h"] = (dict(P, AGE=3.0),
                                         dict(cycles=6, dose_mgm2=80.0, sts=list(range(6))))
    S["S04 S02 + STS at 0 h"] = (dict(P, AGE=3.0),
                                 dict(cycles=6, dose_mgm2=80.0, sts=list(range(6)),
                                      sts_delay=0.0))
    S["S05 S02 + intratympanic STS"] = (dict(P, AGE=3.0),
                                        dict(cycles=6, dose_mgm2=80.0,
                                             it_sts=list(range(6))))
    S["S06 S02 + systemic NAC 150 mg/kg"] = (dict(P, AGE=3.0),
                                             dict(cycles=6, dose_mgm2=80.0,
                                                  nac=list(range(6)), nac_mg=2250.0))
    S["S07 600 mg/m2 as 20 x5d q3w"] = (dict(P), dict(cycles=6, dose_mgm2=100.0, ndays=5))
    S["S08 carboplatin 560 mg/m2 x6"] = (dict(P), dict(cycles=6, dose_mgm2=560.0,
                                                       carbo=True))
    S["S08b myeloablative carbo x3"] = (dict(P), dict(cycles=3, dose_mgm2=1500.0,
                                                      carbo=True))
    S["S09 S01 + aminoglycoside"] = (dict(P, AMGLY=1.0), dict(cycles=6, dose_mgm2=100.0))
    S["S10 S01 + noise 85 dBA"] = (dict(P, NOISE=1.0), dict(cycles=6, dose_mgm2=100.0))
    S["S11 S01 + furosemide"] = (dict(P, FURO=1.0), dict(cycles=6, dose_mgm2=100.0))
    S["S12 S02 in a 1-year-old"] = (dict(P, AGE=1.0), dict(cycles=6, dose_mgm2=80.0))
    S["S13 S02 in a 16-year-old"] = (dict(P, AGE=16.0), dict(cycles=6, dose_mgm2=80.0))
    S["S14 S01 then 24 mo off-treatment"] = (dict(P), dict(cycles=6, dose_mgm2=100.0,
                                                           follow_h=24 * 730.0))
    S["S15 S01, GFR clamped"] = (dict(P, KLOSS=0.0, KPERM=0.0),
                                 dict(cycles=6, dose_mgm2=100.0))
    S["S16 S01, uniform reserve"] = (dict(P, BGSH=0.0), dict(cycles=6, dose_mgm2=100.0))
    S["S17 S01, uniform uptake"] = (dict(P, BUPT=0.0), dict(cycles=6, dose_mgm2=100.0))
    S["S18 S01 + STS at 6 h"] = (dict(P), dict(cycles=6, dose_mgm2=100.0,
                                               sts=list(range(6))))
    S["S19 S01 in CKD (GFR 45)"] = (dict(P, GFR0=45.0), dict(cycles=6, dose_mgm2=100.0))
    S["S20 low-dose 40 mg/m2 weekly x9"] = (dict(P), dict(cycles=9, dose_mgm2=40.0,
                                                          cycle_h=7 * 24.0))

    res = OrderedDict((k, report(simulate(pp, **kw), k)) for k, (pp, kw) in S.items())

    say()
    say("Threshold shift (dB HL) at end of treatment")
    say(hdr() + "  |    PTA Brock SIOP CTCAE ASHA   SNR  tinn")
    for k, r in res.items():
        say(k.ljust(32) + " ".join("{:7.1f}".format(x) for x in r["ts"]) +
            "  | {:6.1f} {:5d} {:4d} {:5d} {:>4s} {:5.1f} {:5.2f}".format(
                r["pta"], r["brock"], r["siop"], r["ctcae"],
                "yes" if r["asha"] else "no", r["snr"], r["tinn"]))
    say()
    say("Supporting state at end of treatment")
    say("{:32s}{:>8s}{:>8s}{:>10s}{:>10s}{:>10s}{:>12s}".format(
        "", "EP mV", "GFR", "OHC 8k", "OHC 4k", "OHC 1k", "tumour LK"))
    for k, r in res.items():
        say("{:32s}{:8.1f}{:8.1f}{:10.3f}{:10.3f}{:10.3f}{:12.2f}".format(
            k, r["ep"], r["gfr"], r["ohc"][6], r["ohc"][4], r["ohc"][2], r["logkill"]))
    say()

    # ---------------------------------------------------------------- PART 3
    say("-" * 104)
    say("PART 3 — THE THIOSULFATE DELAY SWEEP: SIOPEL-6 AND ACCL0431 IN ONE PICTURE")
    say("-" * 104)
    say("Same drug, same dose, same schedule.  Only the delay changes — and the")
    say("cochlea and the tumour read that delay through different lags.")
    say()
    for perf, tag in [(1.0, "well-perfused localised tumour (SIOPEL-6-like)"),
                      (0.15, "poorly-perfused disseminated tumour (ACCL0431-like)")]:
        say("  " + tag)
        base = report(simulate(dict(P, AGE=3.0, TUMPERF=perf), cycles=6, dose_mgm2=80.0))
        say("  {:>8s}{:>10s}{:>10s}{:>8s}{:>8s}{:>12s}{:>14s}{:>12s}".format(
            "delay h", "PTA dB", "8 kHz dB", "Brock", "SIOP", "tumour LK",
            "LK lost %", "otoprot %"))
        for dly in [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]:
            r = report(simulate(dict(P, AGE=3.0, TUMPERF=perf), cycles=6,
                                dose_mgm2=80.0, sts=list(range(6)), sts_delay=dly))
            say("  {:8.1f}{:10.1f}{:10.1f}{:8d}{:8d}{:12.2f}{:14.1f}{:12.1f}".format(
                dly, r["pta"], r["ts"][6], r["brock"], r["siop"], r["logkill"],
                (base["logkill"] - r["logkill"]) / base["logkill"] * 100,
                (base["pta"] - r["pta"]) / base["pta"] * 100))
        say("    reference, no thiosulfate: PTA {:.1f} dB, 8 kHz {:.1f} dB, "
            "log-kill {:.2f}".format(base["pta"], base["ts"][6], base["logkill"]))
        say()

    # ---------------------------------------------------------------- PART 4
    say("-" * 104)
    say("PART 4 — CUMULATIVE DOSE IS THE EXPOSURE VARIABLE; THE GRADES ARE A STAIRCASE")
    say("-" * 104)
    say("{:>8s}{:>12s}{:>9s}{:>9s}{:>9s}{:>9s}{:>8s}{:>7s}{:>7s}".format(
        "cycles", "cum mg/m2", "0.5k", "1k", "2k", "4k", "8k", "Brock", "SIOP"))
    for nc in range(1, 11):
        r = report(simulate(P, cycles=nc, dose_mgm2=100.0))
        say("{:8d}{:12.0f}{:9.1f}{:9.1f}{:9.1f}{:9.1f}{:8.1f}{:7d}{:7d}".format(
            nc, nc * 100.0, r["ts"][1], r["ts"][2], r["ts"][3], r["ts"][4],
            r["ts"][6], r["brock"], r["siop"]))
    say()
    say("Same 600 mg/m2, three schedules — fractionation can only help to the extent")
    say("that the death function is supralinear in instantaneous oxidative stress.")
    say("{:>34s}{:>10s}{:>10s}{:>8s}".format("schedule", "PTA dB", "8 kHz", "Brock"))
    for tag, kw in [("100 mg/m2 q3w x6", dict(cycles=6, dose_mgm2=100.0)),
                    ("20 mg/m2 d1-5 q3w x6", dict(cycles=6, dose_mgm2=100.0, ndays=5)),
                    ("50 mg/m2 q10.5d x12", dict(cycles=12, dose_mgm2=50.0,
                                                 cycle_h=10.5 * 24))]:
        r = report(simulate(P, **kw))
        say("{:>34s}{:10.2f}{:10.2f}{:8d}".format(tag, r["pta"], r["ts"][6], r["brock"]))
    say()
    say("{:>10s}{:>16s}{:>16s}{:>16s}".format(
        "HILL", "bolus PTA dB", "5-day PTA dB", "benefit dB"))
    for h in [1.0, 1.4, 1.7, 2.2, 3.0]:
        pp = dict(P, HILL=h)
        a = report(simulate(pp, cycles=6, dose_mgm2=100.0))["pta"]
        b = report(simulate(pp, cycles=6, dose_mgm2=100.0, ndays=5))["pta"]
        say("{:10.1f}{:16.2f}{:16.2f}{:16.3f}".format(h, a, b, a - b))
    say()

    # ---------------------------------------------------------------- PART 5
    say("-" * 104)
    say("PART 5 — LATE PROGRESSION BELONGS TO THE LABILE POOL, NOT TO THE PLATINUM")
    say("          THAT A MASS SPECTROMETER STILL FINDS DECADES LATER")
    say("-" * 104)
    say("Cochlear platinum is split into a LABILE, redox-active pool (half-life")
    say("TLAB_D) and a BOUND, inert pool (half-life TRET_D) that carries almost all")
    say("of the retained mass.  Only the first drives damage.  Sweeping each in turn")
    say("makes that a falsifiable claim rather than a modelling convenience.")
    say()
    say("  (a) labile half-life TLAB_D  —  bound pool fixed at 2 years")
    say("  {:>16s}{:>11s}{:>11s}{:>11s}{:>11s}{:>14s}{:>13s}".format(
        "TLAB_D", "EOT PTA", "+6 mo", "+12 mo", "+24 mo", "PTA drift dB", "8 kHz drift"))
    for tl, lab in [(7.0, "1 week"), (30.0, "1 month"), (60.0, "2 months"),
                    (120.0, "4 months"), (365.0, "1 year")]:
        pp = dict(P, TLAB_D=tl)
        v, hi = [], []
        for fu in [0.0, 24 * 182.0, 24 * 365.0, 24 * 730.0]:
            rr = report(simulate(pp, cycles=6, dose_mgm2=100.0, follow_h=fu))
            v.append(rr["pta"]); hi.append(rr["ts"][6])
        say("  {:>16s}{:11.2f}{:11.2f}{:11.2f}{:11.2f}{:14.2f}{:13.2f}".format(
            lab, v[0], v[1], v[2], v[3], v[3] - v[0], hi[3] - hi[0]))
    say()
    say("  (b) bound-pool half-life TRET_D  —  labile pool fixed at 60 days")
    say("  {:>16s}{:>11s}{:>11s}{:>14s}{:>22s}".format(
        "TRET_D", "EOT PTA", "+24 mo", "PTA drift dB", "bound Pt at +24 mo"))
    for tr, lab in [(30.0, "1 month"), (365.0, "1 year"), (730.0, "2 years"),
                    (7300.0, "20 years")]:
        pp = dict(P, TRET_D=tr)
        a0 = report(simulate(pp, cycles=6, dose_mgm2=100.0))
        a2 = report(simulate(pp, cycles=6, dose_mgm2=100.0, follow_h=24 * 730.0))
        say("  {:>16s}{:11.2f}{:11.2f}{:14.2f}{:22.4f}".format(
            lab, a0["pta"], a2["pta"], a2["pta"] - a0["pta"], a2["ptb"].mean()))
    say()
    say("  The audiogram does not move when TRET_D is changed by a factor of 240,")
    say("  while the measured cochlear platinum moves by the same factor.  The model")
    say("  therefore predicts that retained platinum is a MARKER of past exposure and")
    say("  not the agent of continuing damage — and that an intervention aimed at")
    say("  chelating it years later would change nothing.")
    say()
    say("Brock / SIOP grade over off-treatment follow-up (adult 600 mg/m2)")
    say("{:>14s}{:>10s}{:>10s}{:>10s}{:>8s}{:>7s}".format(
        "months after", "PTA dB", "4 kHz dB", "8 kHz dB", "Brock", "SIOP"))
    for mo in [0, 3, 6, 12, 24, 48]:
        r = report(simulate(P, cycles=6, dose_mgm2=100.0, follow_h=24 * 30.4 * mo))
        say("{:14d}{:10.2f}{:10.2f}{:10.2f}{:8d}{:7d}".format(
            mo, r["pta"], r["ts"][4], r["ts"][6], r["brock"], r["siop"]))
    say()

    # ---------------------------------------------------------------- PART 6
    say("-" * 104)
    say("PART 6 — WHAT CARRIES THE AUDIOGRAM SLOPE: RESERVE OR DOSE?")
    say("-" * 104)
    full = report(simulate(P, cycles=6, dose_mgm2=100.0))
    nores = report(simulate(dict(P, BGSH=0.0), cycles=6, dose_mgm2=100.0))
    nodose = report(simulate(dict(P, BUPT=0.0), cycles=6, dose_mgm2=100.0))
    neither = report(simulate(dict(P, BGSH=0.0, BUPT=0.0), cycles=6, dose_mgm2=100.0))
    say(hdr())
    say(bl("full model", full["ts"]))
    say(bl("uniform reserve (BGSH=0)", nores["ts"]))
    say(bl("uniform uptake  (BUPT=0)", nodose["ts"]))
    say(bl("both uniform", neither["ts"]))
    say()
    sl = lambda r: r["ts"][6] - r["ts"][2]
    s0 = sl(full)
    say("  audiogram slope, 8 kHz minus 1 kHz (dB)")
    say("    full model                {:8.1f}".format(s0))
    say("    reserve gradient removed  {:8.1f}   ({:5.1f}% of the slope gone)".format(
        sl(nores), (s0 - sl(nores)) / s0 * 100))
    say("    uptake gradient removed   {:8.1f}   ({:5.1f}% of the slope gone)".format(
        sl(nodose), (s0 - sl(nodose)) / s0 * 100))
    say("    both removed              {:8.1f}".format(sl(neither)))
    say()

    # ---------------------------------------------------------------- PART 7
    say("-" * 104)
    say("PART 7 — NEPHROTOXICITY IS AN AMPLIFIER OF OTOTOXICITY, NOT A PARALLEL HARM")
    say("-" * 104)
    say("{:>28s}{:>10s}{:>10s}{:>10s}{:>14s}".format(
        "", "GFR end", "PTA dB", "8 kHz dB", "mean cochl Pt"))
    for tag, pp in [("with renal feedback", dict(P)),
                    ("GFR clamped at baseline", dict(P, KLOSS=0.0, KPERM=0.0)),
                    ("pre-existing CKD, GFR 45", dict(P, GFR0=45.0)),
                    ("CKD 45 + GFR clamped", dict(P, GFR0=45.0, KLOSS=0.0, KPERM=0.0))]:
        r = report(simulate(pp, cycles=6, dose_mgm2=100.0))
        say("{:>28s}{:10.1f}{:10.2f}{:10.2f}{:14.4f}".format(
            tag, r["gfr"], r["pta"], r["ts"][6], r["pt"].mean()))
    say()

    # ---------------------------------------------------------------- PART 8
    say("-" * 104)
    say("PART 8 — ROUTE: A DELIVERY GRADIENT THAT MATCHES THE VULNERABILITY GRADIENT")
    say("-" * 104)
    pk = dict(P, AGE=3.0)
    r_no = report(simulate(pk, cycles=6, dose_mgm2=80.0))
    r_sy = report(simulate(pk, cycles=6, dose_mgm2=80.0, sts=list(range(6))))
    r_it = report(simulate(pk, cycles=6, dose_mgm2=80.0, it_sts=list(range(6))))
    say(hdr())
    say(bl("no protectant       dB", r_no["ts"]))
    say(bl("systemic STS at 6 h dB", r_sy["ts"]))
    say(bl("intratympanic STS   dB", r_it["ts"]))
    say(bl("systemic protection %",
           (r_no["ts"] - r_sy["ts"]) / np.maximum(r_no["ts"], 1e-9) * 100))
    say(bl("intratymp protection %",
           (r_no["ts"] - r_it["ts"]) / np.maximum(r_no["ts"], 1e-9) * 100))
    say()
    dn = np.linspace(0, 72, 2881)
    a_sys = simulate(pk, cycles=1, dose_mgm2=80.0, sts=[0], dense=dn)
    a_it = simulate(pk, cycles=1, dose_mgm2=80.0, it_sts=[0], dense=dn)
    auc = lambda r, st: float(np.trapezoid(r["Y"][IDX[st]], r["t"]))
    say("  one cycle, thiosulfate exposure (uM.h)")
    say("    {:<28s}{:>16s}{:>18s}{:>18s}".format(
        "route", "plasma AUC", "tumour AUC", "peak perilymph"))
    say("    {:<28s}{:16.0f}{:18.1f}{:18.1f}".format(
        "systemic 20 g/m2 at 6 h", auc(a_sys, "STSC"), auc(a_sys, "STSTU"),
        a_sys["Y"][IDX["STSPE"]].max()))
    say("    {:<28s}{:16.4f}{:18.4f}{:18.1f}".format(
        "intratympanic 161 umol", auc(a_it, "STSC"), auc(a_it, "STSTU"),
        a_it["Y"][IDX["STSIT"]].max()))
    say("  The intratympanic route reaches a perilymph peak {:.0f}x higher than the".format(
        a_it["Y"][IDX["STSIT"]].max() / max(a_sys["Y"][IDX["STSPE"]].max(), 1e-9)))
    say("  systemic route while putting nothing measurable into the tumour, which is")
    say("  the only way to escape the efficacy trade-off of Part 3 altogether.")
    say()

    # ---------------------------------------------------------------- PART 9
    say("-" * 104)
    say("PART 9 — RISK-FACTOR DECOMPOSITION (adult 600 mg/m2 backbone)")
    say("-" * 104)
    base = report(simulate(P, cycles=6, dose_mgm2=100.0))
    say("{:>36s}{:>10s}{:>10s}{:>8s}{:>7s}{:>7s}".format(
        "", "PTA dB", "d PTA", "8 kHz", "Brock", "SIOP"))
    say("{:>36s}{:10.2f}{:>10s}{:8.1f}{:7d}{:7d}".format(
        "reference", base["pta"], "-", base["ts"][6], base["brock"], base["siop"]))
    rows = [
        ("age 1 y", dict(P, AGE=1.0), {}),
        ("age 3 y", dict(P, AGE=3.0), {}),
        ("age 16 y", dict(P, AGE=16.0), {}),
        ("concurrent aminoglycoside", dict(P, AMGLY=1.0), {}),
        ("concurrent noise 85 dBA", dict(P, NOISE=1.0), {}),
        ("furosemide", dict(P, FURO=1.0), {}),
        ("pre-existing CKD (GFR 45)", dict(P, GFR0=45.0), {}),
        ("all five risk factors", dict(P, AGE=1.0, AMGLY=1.0, NOISE=1.0, FURO=1.0,
                                       GFR0=45.0), {}),
        ("STS at 6 h", dict(P), dict(sts=list(range(6)))),
        ("intratympanic STS", dict(P), dict(it_sts=list(range(6)))),
        ("STS at 6 h + all five risks", dict(P, AGE=1.0, AMGLY=1.0, NOISE=1.0,
                                             FURO=1.0, GFR0=45.0),
         dict(sts=list(range(6)))),
    ]
    for tag, pp, kw in rows:
        r = report(simulate(pp, cycles=6, dose_mgm2=100.0, **kw))
        say("{:>36s}{:10.2f}{:10.2f}{:8.1f}{:7d}{:7d}".format(
            tag, r["pta"], r["pta"] - base["pta"], r["ts"][6], r["brock"], r["siop"]))
    say()
    say("  additivity check on the speech-frequency average (dB):")
    singles = 0.0
    for tag, pp, kw in rows[:1] + rows[3:7]:
        singles += report(simulate(pp, cycles=6, dose_mgm2=100.0, **kw))["pta"] - base["pta"]
    combo = report(simulate(dict(P, AGE=1.0, AMGLY=1.0, NOISE=1.0, FURO=1.0, GFR0=45.0),
                            cycles=6, dose_mgm2=100.0))["pta"] - base["pta"]
    say("    sum of the five single-factor effects  {:7.2f}".format(singles))
    say("    the five applied together             {:7.2f}".format(combo))
    say("    excess over additivity                {:7.2f}  ({:.2f}x)".format(
        combo - singles, combo / singles))
    say()

    # ---------------------------------------------------------------- PART 10
    say("-" * 104)
    say("PART 10 — NUMERICAL SANITY AND MASS BALANCE")
    say("-" * 104)
    rr = simulate(P, cycles=2, dose_mgm2=100.0, dense=np.linspace(0, 42 * 24, 1200))
    Y = rr["Y"]
    bad = [GLOBAL_STATES[i] for i in range(NG) if Y[i].min() < -1e-6]
    bad += ["band state %d" % i for i in range(NG, NSTATE) if Y[i].min() < -1e-6]
    say("  states going materially negative : {}".format(bad if bad else "none"))
    for nm, b in [("OHC", I_OHC), ("IHC", I_IHC), ("SGN", I_SGN)]:
        say("  {:4s} stays within [0,1]            : {}  (min {:.6f}, max {:.6f})".format(
            nm, bool(Y[b:b + NB].min() >= -1e-9 and Y[b:b + NB].max() <= 1 + 1e-9),
            Y[b:b + NB].min(), Y[b:b + NB].max()))
    say("  EP stays within [0, EP0]         : {}".format(
        bool(Y[IDX["EP"]].min() >= -1e-9 and Y[IDX["EP"]].max() <= P["EP0"] + 1e-6)))
    say("  GSH never negative               : {}  (min {:.6f})".format(
        bool(Y[I_GSH:I_GSH + NB].min() >= -1e-9), Y[I_GSH:I_GSH + NB].min()))
    say("  platinum given per 100 mg/m2 dose : {:.1f} umol; free Cmax {:.2f} uM".format(
        100.0 * P["BSA"] / P["MWCIS"] * 1e3, Y[IDX["CISC"]].max()))
    say()
    say("  CLOSED-FORM CHECK OF THE GLUTATHIONE THRESHOLD")
    say("  At the critical labile load PT* = KSYN*GMAX/(KCON*KROS) the quasi-steady")
    say("  glutathione balance KSYN*(GMAX-G) = KCON*KROS*PT**G/(KMG+G) reduces to")
    say("  G*^2 + KMG*G* - GMAX*KMG = 0, i.e. G* = (-KMG + sqrt(KMG^2+4*GMAX*KMG))/2,")
    say("  with NO free parameter.  The two crossings below are computed from")
    say("  different state variables and must agree if the approximation holds.")
    say("  {:>7s}{:>10s}{:>10s}{:>10s}{:>13s}{:>13s}{:>9s}".format(
        "band", "GMAX", "PT*", "G*", "cyc PT>PT*", "cyc G<G*", "G*/GMAX"))
    perc = simulate(P, cycles=12, dose_mgm2=100.0,
                    dense=np.arange(0, 12 * 21 * 24 + 1, 6.0))
    upt2, gmax2, _ = gradients(P)
    for j in range(NB):
        ptst = P["KSYN"] * gmax2[j] / (P["KCON"] * P["KROS"])
        gst = (-P["KMG"] + math.sqrt(P["KMG"] ** 2 + 4 * gmax2[j] * P["KMG"])) / 2
        ptser, gser = perc["Y"][I_PT + j], perc["Y"][I_GSH + j]
        cp = np.where(ptser >= ptst)[0]
        cg = np.where(gser <= gst)[0]
        f = lambda idx: perc["t"][idx[0]] / (21 * 24.0) if idx.size else float("nan")
        say("  {:>7s}{:10.3f}{:10.3f}{:10.4f}{:13.2f}{:13.2f}{:9.3f}".format(
            FREQ_LABEL[j], gmax2[j], ptst, gst, f(cp), f(cg), gst / gmax2[j]))
    say()
    say("  Bands whose reserve never collapses within 12 cycles show nan in both")
    say("  columns — the same bands, by both criteria.  That agreement is the test.")
    say()
    say("  solver-tolerance check, S01 PTA (dB) at rtol 1e-6 / 1e-7 / 1e-8:")
    vals = []
    for rt in [1e-6, 1e-7, 1e-8]:
        global _RTOL_OVERRIDE
        _RTOL_OVERRIDE = rt
        vals.append(report(simulate(P, cycles=6, dose_mgm2=100.0))["pta"])
    _RTOL_OVERRIDE = None
    say("    {:.5f} / {:.5f} / {:.5f}".format(*vals))
    say()
    say("=" * 104)
    say("END OF REFERENCE RUN")
    say("=" * 104)

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "cio_reference_output.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
