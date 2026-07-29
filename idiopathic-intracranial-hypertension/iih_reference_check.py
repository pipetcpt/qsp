#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
iih_reference_check.py
======================

Dependency-free (Python standard library only) reference transcription of the
44-state IIH QSP model held in ``iih_mrgsolve_model.R``.

Why this file exists
--------------------
No R / mrgsolve interpreter is available in the build environment, so every
number quoted in ``README.md`` and in the commit message is produced *here*,
by integrating the same equations with an explicit RK4 stepper, rather than
being recalled or asserted.  The R file is a transcription of this file (same
parameter names, same right-hand sides, same units); ``ANALYSIS 0`` below
prints a derivative fingerprint at an awkward state so the two
implementations can be compared term by term.

The structural claim the whole model is built on
------------------------------------------------
Intracranial pressure in IIH is *not* a quantity that rises because something
is pumped in faster than it drains.  It is the FIXED POINT of the Davson
relation

        ICP = P_sss  +  R_out * F_form                                    (1)

closed through a mechanically collapsible transverse sinus:

        P_sss = P_cv + Q_v * ( R_ts0 + R_tsf + R_tsc(ICP) )               (2)

Because the sinus is a thin-walled tube sitting *inside* the cranium, its
transmural pressure is (sinus pressure - ICP); when ICP exceeds sinus
pressure by more than the wall can support (PCRIT) the lumen narrows, R_tsc
rises, P_sss rises, and ICP rises further.  (1)+(2) is therefore a loop, and
the loop has a gain

        g = d P_sss / d ICP  =  Q_v * d R_tsc / d ICP                     (3)

Everything interesting follows from (3) rather than from any drug:

  * every intervention's effect on ICP is multiplied by 1/(1-g);
  * when g -> 1 the fixed point is destroyed by a saddle-node (fulminant IIH);
  * between the two folds the same patient has TWO stable pressures, so a
    transient can move the state without moving any parameter, and a
    treatment can be curative by RESETTING the state rather than by being
    continued;
  * cerebral blood flow multiplies R_ts in (2), so a drug that vasodilates
    (acetazolamide, via its own metabolic acidosis) pays a penalty
    proportional to the stenosis it is being amplified by.

Units
-----
time      day        (the mrgsolve model uses the same unit)
pressure  mmHg       (1 mmHg = 1.35951 cmH2O; the clinical literature reports
                      lumbar opening pressure in cmH2O, so both are printed)
flow      mL/min     (converted inside the ODE by *1440)
volume    mL
RNFL      micrometre
weight    kg

Author: Claude Code Routine (QSP disease-model library)
"""

import math
import random

# --------------------------------------------------------------------------
# 0.  small helpers
# --------------------------------------------------------------------------

MMHG_PER_CMH2O = 1.0 / 1.35951
CMH2O = 1.35951          # multiply mmHg by this to get cmH2O


def cm(p_mmhg):
    return p_mmhg * CMH2O


def sig(x):
    """numerically safe logistic"""
    if x > 40.0:
        return 1.0
    if x < -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def pos(x):
    return x if x > 0.0 else 0.0


# --------------------------------------------------------------------------
# 1.  parameters
# --------------------------------------------------------------------------

def params(**over):
    p = dict(

        # ---- cerebral venous outflow / transverse sinus -------------------
        QV0     = 600.0,     # mL/min through the transverse sinuses
        RTS0    = 0.00333,   # mmHg.min/mL  -> 2.0 mmHg normal gradient
        RTSCMAX = 0.0200,    # mmHg.min/mL  -> 12.0 mmHg max compressive add
        PCRIT   = 10.95,     # mmHg  wall support of the sinus (SUSCEPTIBILITY)
        PSTIFF  = 4.14,      # mmHg  steepness of the collapse (SUSCEPTIBILITY)
        FSEG    = 0.50,      # fraction of the sinus gradient upstream of the
                             # collapsible segment (0 = collapse driven by the
                             # jugular pressure, 1 = by the sagittal sinus)
        TAURTS  = 0.0097,    # day (14 min) collapse time constant
        KREM    = 0.0016,    # /day  compression -> fixed remodelled stenosis
        KREV    = 0.0025,    # /day  reversal of remodelling
        RTSFMAX = 0.0100,    # mmHg.min/mL  ceiling of the fixed component
        KRESTEN = 0.00035,   # /day  loss of stent patency (restenosis)

        # ---- CSF hydrodynamics -------------------------------------------
        FFORM0  = 0.350,     # mL/min  choroid-plexus secretion (504 mL/day)
        # CSF outflow resistance is NOT an independent lesion: the arachnoid
        # granulations and the meningeal/cervical lymphatics discharge INTO the
        # venous compartment, so a raised outflow pressure narrows the egress
        # route.  R_out therefore carries a venous-pressure-dependent term,
        # which is the second (amplified) route by which body habitus acts.
        ROUTB   = 0.00,      # mmHg.min/mL  pressure-independent component
        KROUTP  = 1.69,      # mmHg.min/mL per mmHg of central venous pressure
        TAUROUT = 1.00,      # day
        E1      = 0.200,     # /mL  elastance coefficient (C = 1/(E1*ICP))
        RSHUNT  = 2.5,       # mmHg.min/mL  shunt resistance
        PSHUNT  = 7.4,       # mmHg  shunt opening pressure (10 cmH2O)
        KSHFAIL = 0.0011,    # /day  shunt obstruction hazard

        # ---- body habitus -> central venous pressure ----------------------
        HT      = 1.65,      # m
        TAUWT   = 60.0,      # day   weight trajectory time constant (diet)
        TAUWTB  = 25.0,      # day   weight trajectory after bariatric surgery
        KIAP    = 0.350,     # mmHg per BMI unit above BMIREF
        BMIREF  = 22.0,
        KPCV    = 1.00,      # transmission of IAP to central venous pressure
        PCV0    = 4.50,      # mmHg  central venous pressure at BMIREF

        # ---- choroid plexus secretion pharmacology ------------------------
        TAUSEC   = 0.150,    # day  secretory-machinery adaptation (3.6 h)
        EMAXCA   = 0.420,    # max fractional fall in F_form (CA inhibition)
        EC50ACZ  = 12.00,    # mg/L acetazolamide plasma
        EC50TPM  = 9.00,     # mg/L topiramate (weak CA inhibitor)
        POTTPM   = 0.45,     # relative CA potency of topiramate
        FWTTPM   = 0.100,    # max fractional weight loss on topiramate
        ED50WTPM = 90.0,     # mg/day for half-maximal weight loss
        EMAXGLP  = 0.300,    # max fractional fall (GLP-1R agonism)
        EC50EX   = 0.055,    # ng/mL exenatide
        EC50SEM  = 12.0,     # ng/mL semaglutide
        EMAXFUR  = 0.200,    # max fractional fall (NKCC1)
        EC50FUR  = 0.25,     # mg/L furosemide
        EMAXGC   = 0.100,    # glucocorticoid drive on secretion (11b-HSD1)
        EC50GC   = 1.00,
        KANDR    = 0.0120,   # adiposity/androgen drive on secretion, per BMI unit
                             # above 25 (NOT blockable by 11b-HSD1 inhibition)

        # ---- systemic consequences of carbonic-anhydrase inhibition -------
        HCO30    = 25.0,     # mEq/L
        FHCO3    = 0.280,    # max fractional fall in bicarbonate
        TAUHCO3  = 3.00,     # day
        KACID    = 0.550,    # CBF gain per fractional fall in HCO3
        TAUQV    = 0.050,    # day

        # ---- 11beta-HSD1 / cortisol tone ---------------------------------
        TAUCORT  = 0.500,
        KCORTBMI = 0.035,    # per BMI unit above 25
        EMAXAZD  = 0.850,    # AZD4017 target engagement
        EC50AZD  = 55.0,     # ng/mL

        # ---- optic nerve head --------------------------------------------
        TAUONS   = 0.020,    # day  sheath equilibration with ICP
        KFEN     = 22.0,     # /day fenestration conductance (ONSF)
        PORB     = 6.0,      # mmHg orbital tissue pressure
        IOP0     = 15.5,     # mmHg
        KIOPBMI  = 0.100,    # mmHg per BMI unit above 25
        EMAXIOP  = 2.50,     # mmHg max IOP fall from CA inhibition
        TAUIOP   = 0.100,    # day
        TLPDTHR  = 3.00,     # mmHg translaminar gradient for axoplasmic stasis
        TAUSTAS  = 2.00,     # day
        SSTAS    = 2.20,     # mmHg logistic width of the stasis switch
        KSW      = 1.20,     # um/day per mmHg supra-threshold gradient
        RNFLSMAX = 320.0,    # um  ceiling of prelaminar swelling
        KRES     = 0.100,    # /day resorption of swelling
        AXON0    = 98.0,     # um  RNFL axonal thickness
        AXONFL   = 38.0,     # um  floor (glial/residual)
        KAXL     = 0.00120,  # /day axon loss (stasis-driven)
        KISCH    = 0.00600,  # /day axon loss (ischaemic, TLPD > TLPDISCH)
        TLPDISCH = 22.0,     # mmHg
        SMD      = 0.300,    # dB per um of lost axon
        SREV     = 2.00,     # dB max reversible field loss
        TAUMD    = 10.0,     # day

        # ---- headache ----------------------------------------------------
        ICPHTHR  = 13.0,     # mmHg  pressure above which ICP drives headache
        W1HA     = 0.550,
        W2HA     = 0.300,
        W3HA     = 0.250,
        KSENS    = 0.120,    # /day  central sensitisation
        KDESENS  = 0.045,    # /day  desensitisation
        FMOHDES  = 0.600,    # fractional block of desensitisation by MOH
        KMOH     = 0.030,    # /day
        KMOHOFF  = 0.020,    # /day
        TAUMHD   = 3.00,     # day
        TAUCGRP  = 1.00,
        KCG      = 0.700,
        KCGICP   = 0.030,

        # ---- other symptoms ----------------------------------------------
        TINNREF  = 12.0,     # mmHg gradient giving mid-scale tinnitus
        TAUTINN  = 0.200,
        KDIPL    = 0.020,
        KDIPLR   = 0.060,
        ICPDIPL  = 26.0,     # mmHg  threshold for abducens stretch

        # ---- tolerability / adherence ------------------------------------
        EC50PAR  = 6.00,     # mg/L acetazolamide -> paraesthesia
        EC50COG  = 5.50,     # mg/L topiramate    -> cognitive slowing
        TAUAE    = 0.500,
        TAUADH   = 7.00,
        ADHMAX   = 0.950,
        KADHPAR  = 0.450,
        KADHCOG  = 0.500,
        TAUSTONE = 10.0,

        # ---- PK ----------------------------------------------------------
        # acetazolamide: F~1, Vc 14 L, CL 46 L/day (t1/2 ~5 h), deep saturable
        # red-cell carbonic-anhydrase binding
        KAACZ  = 28.8,   VCACZ = 14.0,  CLACZ = 46.0,
        KONACZ = 2.60,   KOFFACZ = 0.90, RBCMAX = 130.0,
        # topiramate: F 0.8, Vc 55 L, t1/2 21 h
        KATPM  = 24.0,   VCTPM = 55.0,  CLTPM = 43.0,  FTPM = 0.80,
        # furosemide: Vc 12 L, t1/2 1.5 h
        KAFUR  = 36.0,   VCFUR = 12.0,  CLFUR = 130.0, FFUR = 0.55,
        # exenatide: SC, t1/2 2.4 h, Vc 28 L
        KAEX   = 8.30,   VCEX  = 28.0,  CLEX  = 190.0,
        # semaglutide: SC weekly, t1/2 165 h, Vc 12.5 L
        KASEM  = 0.55,   VCSEM = 12.5,  CLSEM = 0.87,
        # prednisolone: Vc 35 L, t1/2 3 h
        KAPRED = 24.0,   VCPRED = 35.0, CLPRED = 190.0,
        # AZD4017: Vc 60 L, t1/2 6 h
        KAAZD  = 20.0,   VCAZD = 60.0,  CLAZD = 165.0, FAZD = 0.65,
    )
    p.update(over)
    return p


# --------------------------------------------------------------------------
# 2.  state vector
# --------------------------------------------------------------------------

SNAMES = [
    # PK (13)
    "AACZG", "AACZP", "AACZR", "ATPMG", "ATPMC", "AFURC",
    "AEXSC", "AEXC", "ASEMS", "ASEMC", "APRED", "AAZDG", "AAZDC",
    # secretion / CSF (5)
    "FSEC", "HCO3", "QVREL", "ROUT", "SHUNTP",
    # venous / mechanical (4)
    "ICP", "RTSC", "RTSF", "STENT",
    # body (4)
    "WT", "IAP", "PCV", "CORTL",
    # optic nerve (7)
    "PONS", "IOP", "STASIS", "RNFLS", "AXON", "MD", "CUMEXP",
    # headache / symptoms (6)
    "STG", "CGRP", "MOH", "MHD", "TINN", "DIPL",
    # tolerability (5)
    "PARES", "COGSL", "ADH", "KSTONE", "ONSFS",
]
IX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)
assert NST == 44, NST


# --------------------------------------------------------------------------
# 3.  the algebra of the venous loop  (shared by ODE and steady-state code)
# --------------------------------------------------------------------------

def sinus(p, icp, rtsf, stent, pcv, qv, rtsc=None):
    """Self-consistent transverse-sinus solution at a given ICP.

    Returns (rtsc, rtot, gradient, psss, pseg).
    If ``rtsc`` is given it is treated as the current dynamic state and the
    *target* is returned instead (used by the ODE).  If it is None the
    quasi-steady value is found by damped fixed-point iteration (used by the
    steady-state / continuation analyses).
    """
    rtscmax = p["RTSCMAX"] * (1.0 - stent)
    rfix = p["RTS0"] + rtsf * (1.0 - stent)

    def target(r):
        rtot = rfix + r
        grad = qv * rtot
        pseg = pcv + (1.0 - p["FSEG"]) * grad   # pressure at the collapsible segment
        drive = (icp - pseg - p["PCRIT"]) / p["PSTIFF"]
        return rtscmax * sig(drive)

    if rtsc is None:
        # target(r) - r is strictly decreasing in r (a larger compressive
        # resistance raises the mid-segment pressure, which relieves the
        # collapse), so the inner fixed point is unique and bisectable.
        a, b = 0.0, rtscmax
        if rtscmax <= 0.0:
            r = 0.0
        else:
            for _ in range(44):
                m = 0.5 * (a + b)
                if target(m) - m > 0.0:
                    a = m
                else:
                    b = m
            r = 0.5 * (a + b)
        rtsc_out = r
    else:
        rtsc_out = target(rtsc)          # target for the dynamic state
        r = rtsc

    rtot = rfix + r
    grad = qv * rtot
    psss = pcv + grad
    pseg = pcv + (1.0 - p["FSEG"]) * grad
    return rtsc_out, rtot, grad, psss, pseg


def steady_icp(p, pcv=None, rout=None, fform=None, rtsf=0.0, stent=0.0,
               qvrel=1.0, shuntp=0.0, bracket=(0.5, 140.0), n=4000,
               branch=None):
    """All roots of the Davson fixed-point equation.

    Returns a list of (icp, stability) with stability -1 = stable,
    +1 = unstable (saddle).  ``branch`` picks 'low'/'high'/'nearest'.
    """
    if pcv is None:
        pcv = p["PCV0"]
    if rout is None:
        rout = p["ROUT0"]
    if fform is None:
        fform = p["FFORM0"]
    qv = p["QV0"] * qvrel

    def G(icp):
        _, _, _, psss, _ = sinus(p, icp, rtsf, stent, pcv, qv, rtsc=None)
        fabs = pos(icp - psss) / rout
        fsh = shuntp * pos(icp - p["PSHUNT"]) / p["RSHUNT"]
        return fform - fabs - fsh

    lo, hi = bracket
    xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
    gs = [G(x) for x in xs]
    roots = []
    for i in range(n):
        if gs[i] == 0.0:
            roots.append(xs[i])
        elif gs[i] * gs[i + 1] < 0.0:
            a, b = xs[i], xs[i + 1]
            fa = gs[i]
            for _ in range(80):
                m = 0.5 * (a + b)
                fm = G(m)
                if fa * fm <= 0.0:
                    b = m
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))
    out = []
    for r in roots:
        # dV/dt = G ; stable if dG/dICP < 0
        d = (G(r + 1e-4) - G(r - 1e-4)) / 2e-4
        out.append((r, -1 if d < 0 else +1))
    if branch is None:
        return out
    st = [r for r, s in out if s < 0]
    if not st:
        return None
    if branch == "low":
        return min(st)
    if branch == "high":
        return max(st)
    return st


def loop_gain(p, icp, rtsf=0.0, stent=0.0, pcv=None, qvrel=1.0, h=1e-3):
    """g = dP_sss/dICP at the operating point (equation 3)."""
    if pcv is None:
        pcv = p["PCV0"]
    qv = p["QV0"] * qvrel
    _, _, _, p1, _ = sinus(p, icp + h, rtsf, stent, pcv, qv, rtsc=None)
    _, _, _, p0, _ = sinus(p, icp - h, rtsf, stent, pcv, qv, rtsc=None)
    return (p1 - p0) / (2 * h)


# --------------------------------------------------------------------------
# 4.  regimen description
# --------------------------------------------------------------------------

def regimen(**over):
    r = dict(
        acz=0.0,        # mg/day prescribed
        tpm=0.0,        # mg/day
        fur=0.0,        # mg/day
        ex=0.0,         # ug/day (exenatide 10 ug BID = 20)
        sem=0.0,        # mg/week
        pred=0.0,       # mg/day prednisolone-equivalent
        azd=0.0,        # mg/day AZD4017
        start=0.0,      # day therapy starts
        stop=1e9,       # day therapy stops
        diet=0.0,       # fractional weight loss target (diet/lifestyle)
        diet_start=0.0,
        bariatric=0.0,  # fractional weight loss target (surgery)
        bar_day=1e9,
        stent_day=1e9,
        shunt_day=1e9,
        onsf_day=1e9,
        lp=(),          # tuple of days on which a 25 mL LP is done
        analg=0.25,     # analgesic exposure 0-1 (drives medication overuse)
        citrate=0.0,    # potassium citrate co-prescription 0/1
    )
    r.update(over)
    return r


def _on(ev, t):
    return 1.0 if (ev["start"] <= t < ev["stop"]) else 0.0


# --------------------------------------------------------------------------
# 5.  the right-hand side (transcribed 1:1 into the mrgsolve $ODE block)
# --------------------------------------------------------------------------

def derivs(t, y, p, ev):
    d = [0.0] * NST
    g = IX

    AACZG = y[g["AACZG"]]; AACZP = y[g["AACZP"]]; AACZR = y[g["AACZR"]]
    ATPMG = y[g["ATPMG"]]; ATPMC = y[g["ATPMC"]]; AFURC = y[g["AFURC"]]
    AEXSC = y[g["AEXSC"]]; AEXC = y[g["AEXC"]]
    ASEMS = y[g["ASEMS"]]; ASEMC = y[g["ASEMC"]]
    APRED = y[g["APRED"]]; AAZDG = y[g["AAZDG"]]; AAZDC = y[g["AAZDC"]]
    FSEC = y[g["FSEC"]]; HCO3 = y[g["HCO3"]]; QVREL = y[g["QVREL"]]
    ROUT = y[g["ROUT"]]; SHUNTP = y[g["SHUNTP"]]
    ICP = y[g["ICP"]]; RTSC = y[g["RTSC"]]; RTSF = y[g["RTSF"]]
    STENT = y[g["STENT"]]
    WT = y[g["WT"]]; IAP = y[g["IAP"]]; PCV = y[g["PCV"]]; CORTL = y[g["CORTL"]]
    PONS = y[g["PONS"]]; IOP = y[g["IOP"]]; STASIS = y[g["STASIS"]]
    RNFLS = y[g["RNFLS"]]; AXON = y[g["AXON"]]; MD = y[g["MD"]]
    STG = y[g["STG"]]; CGRP = y[g["CGRP"]]; MOH = y[g["MOH"]]; MHD = y[g["MHD"]]
    TINN = y[g["TINN"]]; DIPL = y[g["DIPL"]]
    PARES = y[g["PARES"]]; COGSL = y[g["COGSL"]]; ADH = y[g["ADH"]]
    ONSFS = y[g["ONSFS"]]

    on = _on(ev, t)

    # ---- 5.1 PK ---------------------------------------------------------
    # continuous prescribed-dose input scaled by adherence: this is the QSP
    # abstraction that lets the tolerability -> adherence -> exposure loop
    # close inside the ODE system (see README, "continuous dosing note").
    d[g["AACZG"]] = on * ev["acz"] * ADH - p["KAACZ"] * AACZG
    CACZ = AACZP / p["VCACZ"]                                   # mg/L
    rbc_in = p["KONACZ"] * CACZ * (p["RBCMAX"] - AACZR)
    rbc_out = p["KOFFACZ"] * AACZR
    d[g["AACZP"]] = (p["KAACZ"] * AACZG
                     - p["CLACZ"] / p["VCACZ"] * AACZP
                     - rbc_in + rbc_out)
    d[g["AACZR"]] = rbc_in - rbc_out

    d[g["ATPMG"]] = on * ev["tpm"] * p["FTPM"] * ADH - p["KATPM"] * ATPMG
    CTPM = ATPMC / p["VCTPM"]
    d[g["ATPMC"]] = p["KATPM"] * ATPMG - p["CLTPM"] / p["VCTPM"] * ATPMC

    CFUR = AFURC / p["VCFUR"]
    d[g["AFURC"]] = (on * ev["fur"] * p["FFUR"] * ADH
                     - p["CLFUR"] / p["VCFUR"] * AFURC)

    d[g["AEXSC"]] = on * ev["ex"] - p["KAEX"] * AEXSC
    CEX = AEXC / p["VCEX"]                                      # ng/mL (ug/L)
    d[g["AEXC"]] = p["KAEX"] * AEXSC - p["CLEX"] / p["VCEX"] * AEXC

    d[g["ASEMS"]] = on * ev["sem"] / 7.0 * 1000.0 - p["KASEM"] * ASEMS  # ug/day
    CSEM = ASEMC / p["VCSEM"]
    d[g["ASEMC"]] = p["KASEM"] * ASEMS - p["CLSEM"] / p["VCSEM"] * ASEMC

    CPRED = APRED / p["VCPRED"]
    d[g["APRED"]] = on * ev["pred"] - p["CLPRED"] / p["VCPRED"] * APRED

    d[g["AAZDG"]] = on * ev["azd"] * p["FAZD"] - p["KAAZD"] * AAZDG
    CAZD = AAZDC / p["VCAZD"] * 1000.0                          # ng/mL
    d[g["AAZDC"]] = p["KAAZD"] * AAZDG - p["CLAZD"] / p["VCAZD"] * AAZDC

    BMI = WT / (p["HT"] ** 2)

    # ---- 5.2 carbonic-anhydrase engagement ------------------------------
    XCA = CACZ / p["EC50ACZ"] + p["POTTPM"] * CTPM / p["EC50TPM"]
    INHCA = p["EMAXCA"] * XCA / (1.0 + XCA)          # fraction of F_form
    OCCCA = XCA / (1.0 + XCA)                        # 0-1 target occupancy

    EGLP = p["EMAXGLP"] * (CEX / (p["EC50EX"] + CEX)
                           + CSEM / (p["EC50SEM"] + CSEM))
    EGLP = min(EGLP, p["EMAXGLP"])
    EFUR = p["EMAXFUR"] * CFUR / (p["EC50FUR"] + CFUR)
    EGC = p["EMAXGC"] * CORTL / (p["EC50GC"] + CORTL)

    EANDR = p["KANDR"] * pos(BMI - 25.0)
    FSECT = ((1.0 - INHCA) * (1.0 - EGLP) * (1.0 - EFUR)
             * (1.0 + EGC) * (1.0 + EANDR))
    d[g["FSEC"]] = (FSECT - FSEC) / p["TAUSEC"]
    FFORM = p["FFORM0"] * FSEC

    # ---- 5.3 systemic acidosis -> cerebral blood flow -------------------
    d[g["HCO3"]] = (p["HCO30"] * (1.0 - p["FHCO3"] * OCCCA) - HCO3) / p["TAUHCO3"]
    QVRELT = 1.0 + p["KACID"] * (p["HCO30"] - HCO3) / p["HCO30"]
    d[g["QVREL"]] = (QVRELT - QVREL) / p["TAUQV"]
    QV = p["QV0"] * QVREL

    # ---- 5.4 cortisol tone (11beta-HSD1) --------------------------------
    EAZD = p["EMAXAZD"] * CAZD / (p["EC50AZD"] + CAZD)
    CORTT = ((1.0 + p["KCORTBMI"] * pos(BMI - 25.0)) * (1.0 - EAZD)
             + 0.18 * CPRED)
    d[g["CORTL"]] = (CORTT - CORTL) / p["TAUCORT"]

    # ---- 5.5 the venous loop --------------------------------------------
    RTSCT, RTOT, GRAD, PSSS, PSEG = sinus(p, ICP, RTSF, STENT, PCV, QV,
                                          rtsc=RTSC)
    d[g["RTSC"]] = (RTSCT - RTSC) / p["TAURTS"]
    d[g["RTSF"]] = (p["KREM"] * (RTSC / max(p["RTSCMAX"], 1e-9))
                    * (p["RTSFMAX"] - RTSF) * (1.0 - STENT)
                    - p["KREV"] * RTSF)
    d[g["STENT"]] = -p["KRESTEN"] * STENT

    # ---- 5.6 ICP ---------------------------------------------------------
    FABS = pos(ICP - PSSS) / ROUT
    FSH = SHUNTP * pos(ICP - p["PSHUNT"]) / p["RSHUNT"]
    FLP = 2.5 if any(dd <= t < dd + 25.0 / 2.5 / 1440.0 for dd in ev["lp"]) else 0.0
    d[g["ICP"]] = (FFORM - FABS - FSH - FLP) * 1440.0 * p["E1"] * ICP
    d[g["ROUT"]] = (p["ROUTB"] + p["KROUTP"] * PCV - ROUT) / p["TAUROUT"]
    d[g["SHUNTP"]] = -p["KSHFAIL"] * SHUNTP

    # ---- 5.7 body habitus ------------------------------------------------
    WT0 = ev.get("_wt0", WT)
    fdiet = ev["diet"] if t >= ev["diet_start"] else 0.0
    fbar = ev["bariatric"] if t >= ev["bar_day"] else 0.0
    fster = 0.0016 * ev["pred"] * on          # steroid weight gain
    ftpm = (p["FWTTPM"] * ev["tpm"] / (p["ED50WTPM"] + ev["tpm"])) * on
    fglpw = (0.055 * ev["ex"] / (14.0 + ev["ex"])
             + 0.150 * ev["sem"] / (1.4 + ev["sem"])) * on
    if fbar > 0.0:
        WTT = WT0 * (1.0 - fbar - ftpm - fglpw + fster)
        tauw = p["TAUWTB"]
    else:
        WTT = WT0 * (1.0 - fdiet - ftpm - fglpw + fster)
        tauw = p["TAUWT"]
    d[g["WT"]] = (WTT - WT) / tauw
    IAPT = p["KIAP"] * pos(BMI - p["BMIREF"])
    d[g["IAP"]] = (IAPT - IAP) / 1.0
    d[g["PCV"]] = (p["PCV0"] + p["KPCV"] * IAP - PCV) / 1.0

    # ---- 5.8 optic nerve head -------------------------------------------
    d[g["ONSFS"]] = ((1.0 if t >= ev["onsf_day"] else 0.0) - ONSFS) / 0.5
    d[g["PONS"]] = ((ICP - PONS) / p["TAUONS"]
                    - ONSFS * p["KFEN"] * (PONS - p["PORB"]))
    IOPT = (p["IOP0"] + p["KIOPBMI"] * pos(BMI - 25.0)
            - p["EMAXIOP"] * OCCCA)
    d[g["IOP"]] = (IOPT - IOP) / p["TAUIOP"]
    TLPD = PONS - IOP
    STAST = sig((TLPD - p["TLPDTHR"]) / p["SSTAS"])
    d[g["STASIS"]] = (STAST - STASIS) / p["TAUSTAS"]
    d[g["RNFLS"]] = (p["KSW"] * pos(TLPD - p["TLPDTHR"])
                     * (1.0 - RNFLS / p["RNFLSMAX"])
                     * (AXON / p["AXON0"])
                     - p["KRES"] * RNFLS)
    d[g["AXON"]] = -(p["KAXL"] * STASIS ** 3
                     + p["KISCH"] * pos(TLPD - p["TLPDISCH"])) \
        * (AXON - p["AXONFL"])
    MDT = -(p["SMD"] * (p["AXON0"] - AXON)
            + p["SREV"] * min(1.0, pos(TLPD - p["TLPDTHR"]) / 15.0))
    d[g["MD"]] = (MDT - MD) / p["TAUMD"]
    d[g["CUMEXP"]] = pos(TLPD - p["TLPDTHR"])

    # ---- 5.9 headache ----------------------------------------------------
    DHA = (p["W1HA"] * pos(ICP - p["ICPHTHR"]) / 10.0
           + p["W2HA"] * pos(GRAD - 4.0) / 10.0
           + p["W3HA"] * pos(CGRP - 1.0))
    d[g["STG"]] = (p["KSENS"] * DHA * (1.0 - STG)
                   - p["KDESENS"] * STG * (1.0 - p["FMOHDES"] * MOH))
    d[g["CGRP"]] = (1.0 + p["KCG"] * STG
                    + p["KCGICP"] * pos(ICP - 15.0) - CGRP) / p["TAUCGRP"]
    d[g["MOH"]] = (p["KMOH"] * ev["analg"] * (1.0 - MOH)
                   - p["KMOHOFF"] * MOH * (1.0 - ev["analg"]))
    MHDT = 30.0 * min(1.0, 0.04 + 0.40 * STG + 0.20 * MOH + 0.26 * DHA)
    d[g["MHD"]] = (MHDT - MHD) / p["TAUMHD"]

    # ---- 5.10 tinnitus / diplopia ---------------------------------------
    TINNT = 10.0 * min(1.0, (GRAD / p["TINNREF"]) ** 1.5)
    d[g["TINN"]] = (TINNT - TINN) / p["TAUTINN"]
    d[g["DIPL"]] = (p["KDIPL"] * pos(ICP - p["ICPDIPL"]) * (1.0 - DIPL)
                    - p["KDIPLR"] * DIPL)

    # ---- 5.11 tolerability / adherence ----------------------------------
    d[g["PARES"]] = ((CACZ / (p["EC50PAR"] + CACZ)
                      + 0.45 * CTPM / (p["EC50COG"] + CTPM)) - PARES) / p["TAUAE"]
    d[g["COGSL"]] = ((CTPM / (p["EC50COG"] + CTPM)) - COGSL) / p["TAUAE"]
    ADHT = p["ADHMAX"] * (1.0 - p["KADHPAR"] * min(1.0, PARES)
                          - p["KADHCOG"] * min(1.0, COGSL))
    ADHT = max(0.05, ADHT)
    d[g["ADH"]] = (ADHT - ADH) / p["TAUADH"]
    d[g["KSTONE"]] = ((OCCCA * (1.0 - 0.60 * ev["citrate"])) - y[g["KSTONE"]]) \
        / p["TAUSTONE"]

    return d


# --------------------------------------------------------------------------
# 6.  initial condition = the patient's own untreated fixed point
# --------------------------------------------------------------------------

def fsec_base(p, bmi):
    """untreated secretory drive: adiposity-linked glucocorticoid (11b-HSD1)
    and androgen tone on the choroid plexus"""
    cortl = 1.0 + p["KCORTBMI"] * pos(bmi - 25.0)
    egc = p["EMAXGC"] * cortl / (p["EC50GC"] + cortl)
    eandr = p["KANDR"] * pos(bmi - 25.0)
    return (1.0 + egc) * (1.0 + eandr)


def init_state(p, bmi=38.0, rout=None, branch="high"):
    y = [0.0] * NST
    g = IX
    wt = bmi * p["HT"] ** 2
    iap = p["KIAP"] * pos(bmi - p["BMIREF"])
    pcv = p["PCV0"] + p["KPCV"] * iap
    rout = (p["ROUTB"] + p["KROUTP"] * pcv) if rout is None else rout
    fs0 = fsec_base(p, bmi)

    # steady venous / CSF fixed point
    roots = steady_icp(p, pcv=pcv, rout=rout, fform=p["FFORM0"] * fs0)
    stable = [r for r, s in roots if s < 0]
    if not stable:
        icp = 45.0
    elif branch == "high":
        icp = max(stable)
    else:
        icp = min(stable)
    rtsc, rtot, grad, psss, pseg = sinus(p, icp, 0.0, 0.0, pcv,
                                        p["QV0"], rtsc=None)

    y[g["FSEC"]] = fs0
    y[g["HCO3"]] = p["HCO30"]
    y[g["QVREL"]] = 1.0
    y[g["ROUT"]] = rout
    y[g["SHUNTP"]] = 0.0
    y[g["ICP"]] = icp
    y[g["RTSC"]] = rtsc
    y[g["RTSF"]] = 0.0
    y[g["STENT"]] = 0.0
    y[g["WT"]] = wt
    y[g["IAP"]] = iap
    y[g["PCV"]] = pcv
    y[g["CORTL"]] = 1.0 + p["KCORTBMI"] * pos(bmi - 25.0)
    y[g["PONS"]] = icp
    y[g["IOP"]] = p["IOP0"] + p["KIOPBMI"] * pos(bmi - 25.0)
    tlpd = y[g["PONS"]] - y[g["IOP"]]
    y[g["STASIS"]] = sig((tlpd - p["TLPDTHR"]) / p["SSTAS"])
    y[g["AXON"]] = p["AXON0"]
    y[g["RNFLS"]] = 0.0
    y[g["MD"]] = 0.0
    y[g["STG"]] = 0.0
    y[g["CGRP"]] = 1.0
    y[g["MOH"]] = 0.0
    y[g["MHD"]] = 1.5
    y[g["TINN"]] = 0.0
    y[g["DIPL"]] = 0.0
    y[g["ADH"]] = p["ADHMAX"]
    return y


def burn_in(p, y, ev, days, dt=0.004):
    """Let the slow disease states (papilledema, sensitisation, remodelling)
    reach the values they would have at presentation."""
    return simulate(p, y, ev, days, dt=dt, out_every=None)[1]


# --------------------------------------------------------------------------
# 7.  integrator
# --------------------------------------------------------------------------

def simulate(p, y0, ev, tend, dt=0.002, out_every=1.0, t0=0.0):
    ev = dict(ev)
    ev.setdefault("_wt0", y0[IX["WT"]])
    y = list(y0)
    t = t0
    rec = []
    nstep = int(round((tend - t0) / dt))
    nout = max(1, int(round(out_every / dt))) if out_every else None
    if nout:
        rec.append((t, list(y)))
    for k in range(nstep):
        # explicit events
        if ev["stent_day"] <= t < ev["stent_day"] + dt:
            y[IX["STENT"]] = 1.0
            y[IX["RTSC"]] = 0.0
        if ev["shunt_day"] <= t < ev["shunt_day"] + dt:
            y[IX["SHUNTP"]] = 1.0
        k1 = derivs(t, y, p, ev)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2 = derivs(t + 0.5 * dt, y2, p, ev)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3 = derivs(t + 0.5 * dt, y3, p, ev)
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4 = derivs(t + dt, y4, p, ev)
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        # hard physical floors
        y[IX["ICP"]] = max(y[IX["ICP"]], 0.5)
        y[IX["AXON"]] = min(max(y[IX["AXON"]], p["AXONFL"]), p["AXON0"])
        y[IX["RNFLS"]] = max(y[IX["RNFLS"]], 0.0)
        for nm in ("STASIS", "STG", "MOH", "DIPL", "STENT", "SHUNTP", "ONSFS"):
            y[IX[nm]] = min(max(y[IX[nm]], 0.0), 1.0)
        t += dt
        if nout and (k + 1) % nout == 0:
            rec.append((t, list(y)))
    return rec, y


def outputs(p, y):
    """derived (non-state) quantities, identical to the mrgsolve $TABLE block"""
    g = IX
    qv = p["QV0"] * y[g["QVREL"]]
    _, rtot, grad, psss, pseg = sinus(p, y[g["ICP"]], y[g["RTSF"]],
                                      y[g["STENT"]], y[g["PCV"]], qv,
                                      rtsc=y[g["RTSC"]])
    tlpd = y[g["PONS"]] - y[g["IOP"]]
    rnfl = y[g["AXON"]] + y[g["RNFLS"]]
    return dict(
        ICP=y[g["ICP"]], ICPCM=cm(y[g["ICP"]]),
        GRAD=grad, PSSS=psss, PCV=y[g["PCV"]],
        FFORM=p["FFORM0"] * y[g["FSEC"]],
        TLPD=tlpd, RNFL=rnfl, AXON=y[g["AXON"]], RNFLS=y[g["RNFLS"]],
        FRISEN=min(5.0, 5.0 * y[g["RNFLS"]] / (y[g["RNFLS"]] + 70.0)),
        MD=y[g["MD"]], MHD=y[g["MHD"]], TINN=y[g["TINN"]],
        DIPL=y[g["DIPL"]], WT=y[g["WT"]], BMI=y[g["WT"]] / p["HT"] ** 2,
        IOP=y[g["IOP"]], HCO3=y[g["HCO3"]], ADH=y[g["ADH"]],
        PARES=y[g["PARES"]], CACZ=y[g["AACZP"]] / p["VCACZ"],
        G=loop_gain(p, y[g["ICP"]], y[g["RTSF"]], y[g["STENT"]],
                    y[g["PCV"]], y[g["QVREL"]]),
    )


# --------------------------------------------------------------------------
# 8.  reporting helpers
# --------------------------------------------------------------------------

LINES = []
OUTFILE = "iih_reference_output.txt"
_FH = [None]


def say(s=""):
    print(s)
    LINES.append(s)
    if _FH[0] is None:
        _FH[0] = open(OUTFILE, "w")
    _FH[0].write(s + "\n")
    _FH[0].flush()


def head(n, title):
    say()
    say("=" * 78)
    say("ANALYSIS %s -- %s" % (n, title))
    say("=" * 78)


# --------------------------------------------------------------------------
# 9.  steady-state convenience wrappers used by the analyses
# --------------------------------------------------------------------------

def pcv_of(p, bmi):
    return p["PCV0"] + p["KPCV"] * p["KIAP"] * pos(bmi - p["BMIREF"])


def rout_of(p, bmi):
    return p["ROUTB"] + p["KROUTP"] * pcv_of(p, bmi)


def ss(p, bmi=38.0, fmul=1.0, rmul=1.0, pcvadd=0.0, rtsf=0.0, stent=0.0,
       qvrel=1.0, branch="high", nscan=1800):
    """steady ICP and the decomposition of equation (1)+(2)"""
    pcv = pcv_of(p, bmi) + pcvadd
    rout = rout_of(p, bmi) * rmul
    fform = p["FFORM0"] * fsec_base(p, bmi) * fmul
    roots = steady_icp(p, pcv=pcv, rout=rout, fform=fform, rtsf=rtsf,
                       stent=stent, qvrel=qvrel, n=nscan)
    stable = [r for r, s in roots if s < 0]
    if not stable:
        return None
    icp = max(stable) if branch == "high" else min(stable)
    _, rtot, grad, psss, pseg = sinus(p, icp, rtsf, stent, pcv,
                                      p["QV0"] * qvrel, rtsc=None)
    g = loop_gain(p, icp, rtsf, stent, pcv, qvrel)
    return dict(ICP=icp, GRAD=grad, PSSS=psss, PCV=pcv, ROUT=rout,
                FFORM=fform, G=g, AMP=1.0 / (1.0 - g) if g < 1 else float("inf"),
                NROOT=len(stable))


def presentation(p, bmi=38.0, days=120.0, ev=None):
    """untreated disease allowed to establish papilledema / sensitisation"""
    y0 = init_state(p, bmi=bmi)
    _, y = simulate(p, y0, ev or regimen(), days, dt=0.002, out_every=None)
    return y


# --------------------------------------------------------------------------
# ANALYSIS 0 -- derivative fingerprint (for cross-checking the R model)
# --------------------------------------------------------------------------

def analysis0():
    head(0, "derivative fingerprint at a deliberately awkward state")
    p = params()
    y = presentation(p, 38.0, 60.0)
    # load every drug at once, mid-stent, mid-shunt, post-ONSF, LP running
    y[IX["AACZG"]] = 180.0;  y[IX["AACZP"]] = 320.0; y[IX["AACZR"]] = 90.0
    y[IX["ATPMG"]] = 40.0;   y[IX["ATPMC"]] = 260.0
    y[IX["AFURC"]] = 6.0
    y[IX["AEXSC"]] = 3.1;    y[IX["AEXC"]] = 1.4
    y[IX["ASEMS"]] = 640.0;  y[IX["ASEMC"]] = 210.0
    y[IX["APRED"]] = 22.0
    y[IX["AAZDG"]] = 130.0;  y[IX["AAZDC"]] = 4.2
    y[IX["STENT"]] = 0.55;   y[IX["SHUNTP"]] = 0.70; y[IX["ONSFS"]] = 0.80
    y[IX["RTSF"]] = 0.0045;  y[IX["MOH"]] = 0.42
    ev = regimen(acz=2000, tpm=100, fur=40, ex=20, sem=2.4, pred=20, azd=800,
                 diet=0.05, analg=0.7, citrate=1.0, lp=(0.0,))
    ev["_wt0"] = y[IX["WT"]]
    d = derivs(0.001, y, p, ev)
    say("t = 0.001 d, all seven drugs on board, stent 0.55, shunt 0.70,")
    say("ONSF 0.80, lumbar puncture draining.  State / derivative pairs:")
    say("")
    say("%-8s %18s %18s" % ("state", "value", "d/dt"))
    for nm in SNAMES:
        say("%-8s %18.10g %18.10g" % (nm, y[IX[nm]], d[IX[nm]]))
    say("")
    say("sum(|d/dt|) = %.10g" % sum(abs(x) for x in d))


# --------------------------------------------------------------------------
# ANALYSIS 1 -- the loop, decomposed
# --------------------------------------------------------------------------

def analysis1():
    head(1, "ICP is a fixed point, not an accumulation: term-by-term decomposition")
    p = params()
    say("Davson:  ICP = P_cv + Q_v*(R_ts0 + R_tsc(ICP)) + R_out*F_form")
    say("")
    say("%-6s %8s %8s %8s %8s %8s %7s %7s" %
        ("BMI", "ICP(cm)", "P_cv", "grad", "Rout*F", "Rout", "g", "1/(1-g)"))
    for bmi in (22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 48):
        r = ss(p, bmi)
        say("%-6.0f %8.1f %8.2f %8.2f %8.2f %8.2f %7.3f %7.2f" %
            (bmi, cm(r["ICP"]), r["PCV"], r["GRAD"],
             r["ROUT"] * r["FFORM"], r["ROUT"], r["G"], r["AMP"]))
    say("")
    say("The three terms are not interchangeable.  P_cv appears BOTH in the")
    say("pressure balance and (through P_seg) in the collapse drive, so its")
    say("net coefficient is close to 1 -- it is NOT amplified.  Anything")
    say("entering through R_out*F_form is amplified by 1/(1-g).")
    say("")
    r38 = ss(p, 38.0)
    icp0 = r38["ICP"]
    dp = 1.0
    a_pcv = (ss(p, 38.0, pcvadd=dp)["ICP"] - icp0) / dp
    # 1 mmHg delivered through the CSF term
    dfrac = dp / (r38["ROUT"] * r38["FFORM"])
    a_csf = (ss(p, 38.0, fmul=1.0 - dfrac)["ICP"] - icp0) / (-dp)
    say("open-loop delivery of exactly 1.00 mmHg, measured at BMI 38:")
    say("  through P_cv (weight, posture, jugular)   dICP = %+.3f mmHg" % a_pcv)
    say("  through R_out*F_form (drug, outflow)      dICP = %+.3f mmHg" % a_csf)
    say("  ratio = %.2f   (1/(1-g) = %.2f)" % (a_csf / a_pcv, r38["AMP"]))
    say("")
    say("sensitivity to FSEG, the position along the stenosis at which the")
    say("collapse is driven (0.5 = mid-segment, the default):")
    say("%-7s %9s %9s %9s %9s" % ("FSEG", "ICP(cm)", "g", "P_cv coef", "CSF coef"))
    for fs in (0.0, 0.25, 0.50, 0.75, 1.0):
        q = params(FSEG=fs)
        r = ss(q, 38.0)
        i0 = r["ICP"]
        ap = (ss(q, 38.0, pcvadd=1.0)["ICP"] - i0)
        df = 1.0 / (r["ROUT"] * r["FFORM"])
        ac = (ss(q, 38.0, fmul=1.0 - df)["ICP"] - i0) / (-1.0)
        say("%-7.2f %9.1f %9.3f %9.3f %9.3f" % (fs, cm(i0), r["G"], ap, ac))
    say("")
    say("The unamplified status of body weight is a consequence of FSEG=0.5;")
    say("if the collapse were driven at the jugular end (FSEG=0) weight would")
    say("be amplified too.  This is the model's single largest structural")
    say("uncertainty and it is recorded here rather than hidden.")


# --------------------------------------------------------------------------
# ANALYSIS 2 -- one mechanical parameter generates the responder spectrum
# --------------------------------------------------------------------------

def analysis2():
    head(2, "identical target engagement, different patients: 1/(1-g)")
    say("Every row below receives EXACTLY the same molecular effect: a 30 %")
    say("reduction in choroid-plexus secretion (the GLP-1R effect size).")
    say("Nothing pharmacological differs between rows.  Only PCRIT -- how")
    say("much negative transmural pressure the transverse sinus wall can")
    say("carry before it narrows -- differs.")
    say("")
    say("%-8s %9s %8s %8s %10s %10s %9s" %
        ("PCRIT", "ICP(cm)", "grad", "g", "1/(1-g)", "dICP(cm)", "n roots"))
    for pc in (16.0, 14.0, 13.0, 12.0, 11.5, 10.95, 10.5, 10.0, 9.5, 9.0, 8.5):
        q = params(PCRIT=pc)
        b = ss(q, 38.0)
        if b is None:
            say("%-8.2f %9s -- no stable fixed point (fulminant)" % (pc, ""))
            continue
        t = ss(q, 38.0, fmul=0.70)
        say("%-8.2f %9.1f %8.2f %8.3f %10.2f %10.2f %9d" %
            (pc, cm(b["ICP"]), b["GRAD"], b["G"], b["AMP"],
             cm(t["ICP"] - b["ICP"]), b["NROOT"]))
    say("")
    say("The clinical dichotomy 'acetazolamide responder / non-responder'")
    say("needs no pharmacogenetics: the same 30 % secretion cut buys")
    say("2.6 cmH2O in a stiff-sinus patient and >8 cmH2O in a compliant one.")


# --------------------------------------------------------------------------
# ANALYSIS 3 -- folds, bistability, and treat-then-withdraw
# --------------------------------------------------------------------------

def analysis3():
    head(3, "where the fixed point folds -- and the finding that the"
            " trial-calibrated patient is NOT there")
    say("The loop can in principle destroy its own low-pressure solution.")
    say("Whether it does is decided by one number: the maximum slope of the")
    say("collapse map, which for a logistic collapse is")
    say("")
    say("    h'max = FSEG * Q_v * RTSCMAX / (4 * PSTIFF)")
    say("")
    say("and multiplicity requires h'max > 1.  At the calibrated values")
    p = params()
    hmax = p["FSEG"] * p["QV0"] * p["RTSCMAX"] / (4.0 * p["PSTIFF"])
    say("    FSEG %.2f, Q_v %.0f, RTSCMAX %.4f, PSTIFF %.2f  ->  h'max = %.3f"
        % (p["FSEG"], p["QV0"], p["RTSCMAX"], p["PSTIFF"], hmax))
    say("")
    say("so h'max < 1 and the calibrated patient has EXACTLY ONE stable")
    say("pressure at every BMI.  This is a result, not an omission: the")
    say("parameters were fixed by three trial anchors before this analysis was")
    say("run, and they place the typical patient safely below the fold.")
    say("")
    say("(a) how compliant does the sinus have to be before the fold appears?")
    say("    PSTIFF is swept; PCRIT is re-solved each time so that the")
    say("    presenting pressure stays inside the clinical range.")
    say("")
    say("%-9s %9s %30s %6s %9s" %
        ("PSTIFF", "h'max", "stable ICP at BMI 38 (cmH2O)", "n", "g(upper)"))
    for ps in (4.14, 3.00, 2.20, 1.80, 1.50, 1.40, 1.20, 1.00, 0.80, 0.60):
        hm = p["FSEG"] * p["QV0"] * p["RTSCMAX"] / (4.0 * ps)
        # keep the upper branch near 35 cmH2O by re-solving PCRIT
        best, bpc = None, None
        for k in range(0, 45):
            pc = 4.0 + 0.25 * k
            q = params(PSTIFF=ps, PCRIT=pc)
            r = ss(q, 38.0, nscan=400)
            if r is None:
                continue
            e = abs(cm(r["ICP"]) - 35.1)
            if best is None or e < best:
                best, bpc = e, pc
        q = params(PSTIFF=ps, PCRIT=bpc)
        pcv = pcv_of(q, 38.0); rout = rout_of(q, 38.0)
        ff = q["FFORM0"] * fsec_base(q, 38.0)
        roots = steady_icp(q, pcv=pcv, rout=rout, fform=ff, n=3000)
        st = [r for r, sg in roots if sg < 0]
        gu = loop_gain(q, max(st), 0.0, 0.0, pcv, 1.0) if st else float("nan")
        say("%-9.2f %9.3f %30s %6d %9.3f" %
            (ps, hm, ", ".join("%.1f" % cm(x) for x in st), len(st), gu))
    say("")
    say("(b) continuation in the CSF-side parameter K = R_out*F_form for a")
    say("    COMPLIANT-SINUS patient (PSTIFF 1.0), showing the two folds and")
    say("    therefore the hysteresis loop")
    say("")
    best, bpc = None, None
    for k in range(0, 45):
        pc = 4.0 + 0.25 * k
        q = params(PSTIFF=1.0, PCRIT=pc)
        r = ss(q, 38.0, nscan=400)
        if r is None:
            continue
        e = abs(cm(r["ICP"]) - 35.1)
        if best is None or e < best:
            best, bpc = e, pc
    q = params(PSTIFF=1.0, PCRIT=bpc)
    say("    (PSTIFF 1.0, PCRIT solved to %.2f)" % bpc)
    say("%-8s %10s %34s %6s" % ("fmul", "K(mmHg)", "stable ICP (cmH2O)", "n"))
    fm = 1.30
    lowfold = highfold = None
    prev = None
    while fm >= 0.30:
        pcv = pcv_of(q, 38.0); rout = rout_of(q, 38.0)
        ff = q["FFORM0"] * fsec_base(q, 38.0) * fm
        roots = steady_icp(q, pcv=pcv, rout=rout, fform=ff, n=3200)
        st = [r for r, sg in roots if sg < 0]
        say("%-8.2f %10.2f %34s %6d" %
            (fm, rout * ff, ", ".join("%.1f" % cm(x) for x in st), len(st)))
        if prev is not None:
            if prev == 1 and len(st) > 1:
                highfold = fm
            if prev > 1 and len(st) == 1:
                lowfold = fm
        prev = len(st)
        fm -= 0.05
    say("")
    say("    upper fold near fmul %s, lower fold near fmul %s"
        % (highfold, lowfold))
    say("")
    say("(c) treat for 6 months, then STOP.  Same drug, same BMI, two")
    say("    sinus phenotypes.")
    say("")
    say("%-34s %9s %9s %10s %11s" %
        ("patient", "presn", "on ACZ", "+90 d off", "+365 d off"))
    for lab, kw in (("calibrated (PSTIFF 4.14) monostable", {}),
                    ("compliant sinus (PSTIFF 1.5)",
                     dict(PSTIFF=1.5, PCRIT=9.05)),
                    ("compliant sinus (PSTIFF 1.0)",
                     dict(PSTIFF=1.0, PCRIT=bpc))):
        q2 = params(**kw)
        yp = presentation(q2, 38.0, 120.0)
        a = outputs(q2, yp)["ICPCM"]
        _, y1 = simulate(q2, yp, regimen(acz=2000), 180.0, dt=0.002,
                         out_every=None)
        b = outputs(q2, y1)["ICPCM"]
        _, y2 = simulate(q2, y1, regimen(acz=0), 90.0, dt=0.002, out_every=None)
        c = outputs(q2, y2)["ICPCM"]
        _, y3 = simulate(q2, y2, regimen(acz=0), 275.0, dt=0.002,
                         out_every=None)
        dd = outputs(q2, y3)["ICPCM"]
        say("%-34s %9.1f %9.1f %10.1f %11.1f" % (lab, a, b, c, dd))
    say("")
    say("Read the last two columns.  In the monostable phenotype the pressure")
    say("returns to where it started, so the drug is a parameter shift that")
    say("must be maintained -- which is what is observed in ordinary IIH and")
    say("why relapse on withdrawal is the rule.  A fold, when it exists, is")
    say("what would make a finite course curative, and the model says the")
    say("sinus has to be roughly three times more compliant than the")
    say("trial-calibrated value before that becomes possible.  That is a")
    say("falsifiable statement about WHICH patients could ever stop.")


# ANALYSIS 4 -- four therapeutic targets are four terms of one equation
# --------------------------------------------------------------------------

def analysis4():
    head(4, "four targets, one equation")
    p = params()
    b = ss(p, 38.0)
    say("baseline BMI 38: ICP %.1f cmH2O, gradient %.2f mmHg, g %.3f" %
        (cm(b["ICP"]), b["GRAD"], b["G"]))
    say("")
    say("%-34s %10s %10s %9s %9s" %
        ("intervention (matched -20 % where", "dICP(cm)", "grad", "g after",
         "tinnitus"))
    say("%-34s %10s %10s %9s %9s" %
        ("  applicable)", "", "(mmHg)", "", "(0-10)"))
    rows = [
        ("F_form  -20 % (CA inhibition)", ss(p, 38.0, fmul=0.80)),
        ("R_out   -20 % (shunt / lymphatic)", ss(p, 38.0, rmul=0.80)),
        ("P_cv    -20 % (weight, posture)",
         ss(p, 38.0, pcvadd=-0.20 * b["PCV"])),
        ("R_ts -> 0 (venous sinus stent)", ss(p, 38.0, stent=1.0)),
        ("Q_v +15 % (acidosis vasodilation)", ss(p, 38.0, qvrel=1.15)),
    ]
    for lab, r in rows:
        ti = 10.0 * min(1.0, (r["GRAD"] / p["TINNREF"]) ** 1.5)
        say("%-34s %10.2f %10.2f %9.3f %9.2f" %
            (lab, cm(r["ICP"] - b["ICP"]), r["GRAD"], r["G"], ti))
    say("")
    ti0 = 10.0 * min(1.0, (b["GRAD"] / p["TINNREF"]) ** 1.5)
    say("baseline pulsatile tinnitus %.2f/10" % ti0)
    say("")
    say("Stenting is the only row that changes g: it does not push on a")
    say("drive, it removes the amplifier.  That is why it also abolishes")
    say("pulsatile tinnitus immediately (tinnitus reads the GRADIENT, not")
    say("the pressure) while a drug that lowers ICP by a similar amount")
    say("leaves a residual gradient and residual tinnitus.")


# --------------------------------------------------------------------------
# ANALYSIS 5 -- calibration against the published trials
# --------------------------------------------------------------------------

def analysis5():
    head(5, "calibration: simulated vs published")
    p = params()
    yp = presentation(p, 38.0, 120.0)
    o = outputs(p, yp)
    say("presenting patient: BMI %.1f, ICP %.1f cmH2O, gradient %.1f mmHg,"
        % (o['BMI'], o['ICPCM'], o['GRAD']))
    say("  R_out %.1f mmHg.min/mL, RNFL %.0f um, Frisen %.2f, MD %.2f dB,"
        % (yp[IX['ROUT']], o['RNFL'], o['FRISEN'], o['MD']))
    say("  headache %.0f d/month, pulsatile tinnitus %.1f/10, g = %.2f"
        % (o['MHD'], o['TINN'], o['G']))
    say("")
    arms = [
        ("untreated control", regimen(), 180.0, None),
        ("diet -15.7 % (Sinclair 2010)", regimen(diet=0.157), 90.0, -8.4),
        ("diet -3.3 % (IIHTT placebo arm)", regimen(diet=0.033), 180.0, None),
        ("ACZ 2 g/d + diet (IIHTT active)", regimen(acz=2000, diet=0.070),
         180.0, None),
        ("ACZ 4 g/d + diet (IIHTT target)", regimen(acz=4000, diet=0.070),
         180.0, None),
        ("exenatide 10 ug BID (Mitchell 2023)", regimen(ex=20.0), 84.0, -5.6),
        ("AZD4017 400 mg BID (Markey 2020)", regimen(azd=800.0), 84.0, +0.3),
        ("bariatric surgery (IIH:WT 12 mo)",
         regimen(bariatric=0.21, bar_day=0.0), 365.0, None),
        ("topiramate 150 mg/d (Celebisoy 2007)", regimen(tpm=150), 180.0, None),
        ("furosemide 80 mg/d (adjunct)", regimen(fur=80), 180.0, None),
        ("venous sinus stent", regimen(stent_day=0.0), 180.0, None),
        ("CSF shunt (10 cmH2O valve)", regimen(shunt_day=0.0), 180.0, None),
        ("ONSF (right eye)", regimen(onsf_day=0.0), 180.0, None),
    ]
    say("%-37s %5s %8s %8s %7s %7s %7s %6s %6s" %
        ("arm", "day", "ICP(cm)", "dICP", "grad", "RNFL", "Frisen", "MD",
         "HA d/m"))
    base = o["ICPCM"]
    res = {}
    for lab, ev, dur, target in arms:
        _, y = simulate(p, yp, ev, dur, dt=0.002, out_every=None)
        q = outputs(p, y)
        res[lab] = q
        say("%-37s %5.0f %8.1f %+8.2f %7.2f %7.0f %7.2f %6.2f %6.1f" %
            (lab, dur, q["ICPCM"], q["ICPCM"] - base, q["GRAD"], q["RNFL"],
             q["FRISEN"], q["MD"], q["MHD"]))
    say("")
    say("published anchors and the simulated equivalents")
    say("  Sinclair 2010 (15.7 kg diet, 3 mo): ICP -8.4 cmH2O")
    say("      simulated %+.2f cmH2O"
        % (res["diet -15.7 % (Sinclair 2010)"]["ICPCM"] - base))
    say("  Mitchell 2023 exenatide, 12 wk: ICP -5.6 cmH2O")
    say("      simulated %+.2f cmH2O vs concurrent control %+.2f"
        % (res["exenatide 10 ug BID (Mitchell 2023)"]["ICPCM"] - base,
           res["exenatide 10 ug BID (Mitchell 2023)"]["ICPCM"]
           - res["untreated control"]["ICPCM"]))
    say("  Markey 2020 AZD4017, 12 wk: no significant ICP change")
    say("      simulated %+.2f cmH2O"
        % (res["AZD4017 400 mg BID (Markey 2020)"]["ICPCM"] - base))
    d1 = res["ACZ 2 g/d + diet (IIHTT active)"]["ICPCM"]
    d0 = res["diet -3.3 % (IIHTT placebo arm)"]["ICPCM"]
    say("  IIHTT 2014, 6 mo: ICP difference acetazolamide - placebo"
        " -4.4 cmH2O,")
    say("      papilledema grade difference -0.70, MD difference +0.71 dB")
    say("      simulated ICP difference %+.2f cmH2O, grade %+.2f, MD %+.2f dB"
        % (d1 - d0,
           res["ACZ 2 g/d + diet (IIHTT active)"]["FRISEN"]
           - res["diet -3.3 % (IIHTT placebo arm)"]["FRISEN"],
           res["ACZ 2 g/d + diet (IIHTT active)"]["MD"]
           - res["diet -3.3 % (IIHTT placebo arm)"]["MD"]))
    say("  IIH:WT 2021, bariatric vs diet at 12 mo: ICP difference -6.0 cmH2O")
    say("      simulated bariatric %+.2f cmH2O from presentation"
        % (res["bariatric surgery (IIH:WT 12 mo)"]["ICPCM"] - base))


# --------------------------------------------------------------------------
# ANALYSIS 6 -- acetazolamide: one molecule, three signs
# --------------------------------------------------------------------------

def analysis6():
    head(6, "acetazolamide has three actions with three different signs")
    say("(i)  choroid-plexus carbonic anhydrase -> F_form down    (helps ICP)")
    say("(ii) metabolic acidosis -> cerebral vasodilation -> Q_v up, and Q_v")
    say("     MULTIPLIES the stenosis resistance                  (raises ICP)")
    say("(iii) ciliary-epithelium carbonic anhydrase -> IOP down, and the")
    say("     translaminar gradient is ICP - IOP                  (harms the disc)")
    say("")
    p = params()
    yp = presentation(p, 38.0, 120.0)
    base = outputs(p, yp)
    say("%-9s %8s %8s %7s %7s %7s %7s %7s %7s %7s" %
        ("mg/day", "Cacz", "adher", "HCO3", "Qv rel", "ICP(cm)", "IOP",
         "TLPD", "Frisen", "MD"))
    doses = (0, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000)
    store = {}
    for dz in doses:
        _, y = simulate(p, yp, regimen(acz=dz), 180.0, dt=0.002, out_every=None)
        q = outputs(p, y)
        store[dz] = q
        say("%-9d %8.1f %8.2f %7.1f %7.3f %7.1f %7.2f %7.2f %7.2f %7.2f" %
            (dz, q["CACZ"], q["ADH"], q["HCO3"],
             y[IX["QVREL"]], q["ICPCM"], q["IOP"], q["TLPD"], q["FRISEN"],
             q["MD"]))
    say("")
    say("counterfactual dissection at 2000 mg/day (each arm switched off in")
    say("turn, everything else identical):")
    say("")
    variants = [
        ("full model", params()),
        ("no acidosis-vasodilation arm (KACID=0)", params(KACID=0.0)),
        ("no ocular arm (EMAXIOP=0)", params(EMAXIOP=0.0)),
        ("no adherence loss (KADHPAR=0)", params(KADHPAR=0.0)),
        ("neither penalty arm", params(KACID=0.0, EMAXIOP=0.0)),
    ]
    say("%-40s %9s %8s %8s %8s" %
        ("variant", "ICP(cm)", "TLPD", "Frisen", "MD"))
    for lab, q in variants:
        yq = presentation(q, 38.0, 120.0)
        b = outputs(q, yq)
        _, y = simulate(q, yq, regimen(acz=2000), 180.0, dt=0.002,
                        out_every=None)
        r = outputs(q, y)
        say("%-40s %9.2f %8.2f %8.2f %8.2f" %
            (lab, r["ICPCM"] - b["ICPCM"], r["TLPD"] - b["TLPD"],
             r["FRISEN"] - b["FRISEN"], r["MD"] - b["MD"]))
    say("")
    say("the acidosis penalty is proportional to the stenosis it is amplified")
    say("by, so acetazolamide's NET ICP benefit is not monotone in disease")
    say("severity:")
    say("")
    say("%-9s %9s %9s %10s %10s" %
        ("PCRIT", "grad", "g", "dICP(cm)", "penalty%"))
    for pc in (14.0, 12.0, 11.5, 10.95, 10.5, 10.0, 9.6):
        q = params(PCRIT=pc)
        yq = presentation(q, 38.0, 120.0)
        b = outputs(q, yq)
        _, y = simulate(q, yq, regimen(acz=2000), 180.0, dt=0.002,
                        out_every=None)
        r = outputs(q, y)
        q2 = params(PCRIT=pc, KACID=0.0)
        yq2 = presentation(q2, 38.0, 120.0)
        b2 = outputs(q2, yq2)
        _, y2 = simulate(q2, yq2, regimen(acz=2000), 180.0, dt=0.002,
                         out_every=None)
        r2 = outputs(q2, y2)
        d, d2 = r["ICPCM"] - b["ICPCM"], r2["ICPCM"] - b2["ICPCM"]
        say("%-9.2f %9.2f %9.3f %10.2f %10.1f" %
            (pc, b["GRAD"], b["G"], d, 100.0 * (1.0 - d / d2) if d2 else 0.0))


# --------------------------------------------------------------------------
# ANALYSIS 7 -- the measured RNFL is not the surviving nerve
# --------------------------------------------------------------------------

def analysis7():
    head(7, "OCT retinal nerve fibre layer is a SUM of two states with"
            " opposite meanings")
    say("RNFL(measured) = surviving axons + prelaminar swelling.")
    say("Swelling requires living axons to swell, so a dying disc and a")
    say("recovering disc pass through the SAME measured thickness.")
    say("")
    p = params()
    say("(a) two patients followed to an identical measured RNFL")
    say("%-6s %11s %11s %11s %8s %8s | %11s %11s %11s %8s %8s" %
        ("day", "RNFL", "axons", "swelling", "MD", "ICP",
         "RNFL", "axons", "swelling", "MD", "ICP"))
    say("%-6s %57s | %s" % ("", "  A: treated at presentation (ACZ 2 g + diet)",
                            "  B: untreated for 12 months"))
    pa = params(); pb = params()
    ya = presentation(pa, 38.0, 120.0)
    yb = presentation(pb, 38.0, 120.0)
    eva = regimen(acz=2000, diet=0.10)
    evb = regimen()
    for k in range(9):
        oa, ob = outputs(pa, ya), outputs(pb, yb)
        say("%-6.0f %11.1f %11.1f %11.1f %8.2f %8.1f | %11.1f %11.1f %11.1f"
            " %8.2f %8.1f" %
            (45 * k, oa["RNFL"], oa["AXON"], oa["RNFLS"], oa["MD"],
             oa["ICPCM"], ob["RNFL"], ob["AXON"], ob["RNFLS"], ob["MD"],
             ob["ICPCM"]))
        if k < 8:
            _, ya = simulate(pa, ya, eva, 45.0, dt=0.002, out_every=None)
            _, yb = simulate(pb, yb, evb, 45.0, dt=0.002, out_every=None)
    say("")
    say("(b) the injury endpoint is an INTEGRAL of supra-threshold")
    say("    translaminar pressure, so delay beats dose")
    say("")
    say("%-34s %9s %9s %9s %9s" %
        ("policy (2 y horizon)", "MD(dB)", "axons", "RNFL", "exposure"))
    horizon = 730.0
    for lab, delay, ev in (
            ("full therapy at week 0", 0.0, regimen(acz=2000, diet=0.10)),
            ("full therapy at week 2", 14.0, regimen(acz=2000, diet=0.10)),
            ("full therapy at week 4", 28.0, regimen(acz=2000, diet=0.10)),
            ("full therapy at week 8", 56.0, regimen(acz=2000, diet=0.10)),
            ("full therapy at week 12", 84.0, regimen(acz=2000, diet=0.10)),
            ("full therapy at week 26", 182.0, regimen(acz=2000, diet=0.10)),
            ("HALF therapy at week 0", 0.0, regimen(acz=500, diet=0.05)),
            ("stent at week 0", 0.0, regimen(stent_day=0.0)),
            ("stent at week 12", 84.0, regimen(stent_day=84.0)),
            ("no therapy", 0.0, regimen())):
        q = params()
        y = presentation(q, 38.0, 120.0)
        e = dict(ev); e["start"] = delay
        if e["stent_day"] < 1e8:
            e["stent_day"] = delay
        e["diet_start"] = delay
        _, y = simulate(q, y, e, horizon, dt=0.002, out_every=None)
        r = outputs(q, y)
        say("%-34s %9.2f %9.1f %9.1f %9.0f" %
            (lab, r["MD"], r["AXON"], r["RNFL"], y[IX["CUMEXP"]]))
    say("")
    say("(c) severe phenotype (PCRIT 9.3, BMI 44).  The ischaemic axon-loss")
    say("    term only switches on above a translaminar gradient of %.0f mmHg,"
        % params()["TLPDISCH"])
    say("    which this model reaches only above an opening pressure of about")
    say("    53 cmH2O -- so read the row below as showing that even at 46")
    say("    cmH2O the loss is still the slow stasis-driven kind, and that")
    say("    what separates the options over three weeks is how fast each one")
    say("    removes the gradient, not which injury mechanism it blocks.")
    say("%-30s %9s %9s %9s %9s" % ("", "ICP(cm)", "TLPD", "axons", "MD"))
    q = params(PCRIT=9.3)
    y = presentation(q, 44.0, 30.0)
    for lab, ev, dur in (("presentation (day 30)", None, 0.0),
                         ("+21 d medical therapy", regimen(acz=4000), 21.0),
                         ("+21 d stent instead", regimen(stent_day=0.0), 21.0),
                         ("+21 d shunt instead", regimen(shunt_day=0.0), 21.0),
                         ("+21 d untreated", regimen(), 21.0)):
        if ev is None:
            r = outputs(q, y)
        else:
            _, yy = simulate(q, y, ev, dur, dt=0.002, out_every=None)
            r = outputs(q, yy)
        say("%-30s %9.1f %9.2f %9.1f %9.2f" %
            (lab, r["ICPCM"], r["TLPD"], r["AXON"], r["MD"]))


# --------------------------------------------------------------------------
# ANALYSIS 8 -- optic nerve sheath fenestration decouples, it does not treat
# --------------------------------------------------------------------------

def analysis8():
    head(8, "ONSF changes the optic nerve's boundary condition, not the"
            " patient's pressure")
    p = params()
    yp = presentation(p, 38.0, 120.0)
    b = outputs(p, yp)
    say("%-32s %8s %8s %8s %8s %8s %8s %8s" %
        ("", "ICP(cm)", "PONS", "TLPD", "RNFL", "MD", "HA d/m", "tinnitus"))
    say("%-32s %8.1f %8.2f %8.2f %8.1f %8.2f %8.1f %8.2f" %
        ("presentation", b["ICPCM"], yp[IX["PONS"]], b["TLPD"], b["RNFL"],
         b["MD"], b["MHD"], b["TINN"]))
    for lab, ev in (("ONSF alone", regimen(onsf_day=0.0)),
                    ("ACZ 2 g/d alone", regimen(acz=2000)),
                    ("stent alone", regimen(stent_day=0.0)),
                    ("shunt alone", regimen(shunt_day=0.0)),
                    ("ONSF + ACZ 2 g/d",
                     regimen(onsf_day=0.0, acz=2000))):
        _, y = simulate(p, yp, ev, 180.0, dt=0.002, out_every=None)
        q = outputs(p, y)
        say("%-32s %8.1f %8.2f %8.2f %8.1f %8.2f %8.1f %8.2f" %
            (lab, q["ICPCM"], y[IX["PONS"]], q["TLPD"], q["RNFL"], q["MD"],
             q["MHD"], q["TINN"]))
    say("")
    say("Rank the rows by mean deviation and by headache and you get two")
    say("different orders.  ONSF collapses the translaminar gradient from")
    say("9.8 to 4.1 mmHg and the measured RNFL from 152 to 100 um while")
    say("leaving ICP, headache days and tinnitus exactly where they were,")
    say("because it acts on PONS and PONS is the only place the optic nerve")
    say("reads pressure.  Acetazolamide does the opposite: it takes 3.7 cmH2O")
    say("off the pressure and yet leaves the disc worse than ONSF does,")
    say("because the same carbonic anhydrase inhibition lowers IOP and the")
    say("disc reads the DIFFERENCE.  ONSF also does nothing for the fellow")
    say("eye, which in this model is not a limitation of the operation but an")
    say("identity: the fellow nerve still sees the undiminished ICP.")


# --------------------------------------------------------------------------
# ANALYSIS 9 -- the lumbar puncture as measurement and as therapy
# --------------------------------------------------------------------------

def analysis9():
    head(9, "the lumbar puncture: a 30-minute pressure effect and a"
            " multi-day symptom effect")
    p = params()
    yp = presentation(p, 38.0, 120.0)
    say("(a) fast integration (dt = 0.86 s) of a 25 mL removal")
    say("    relaxation time tau = C * R_out = 1/(E1*ICP) * R_out")
    icp0 = yp[IX["ICP"]]
    tau = (1.0 / (p["E1"] * icp0)) * yp[IX["ROUT"]]
    say("    C = %.3f mL/mmHg at ICP %.1f mmHg  ->  tau = %.1f min"
        % (1.0 / (p["E1"] * icp0), icp0, tau))
    say("")
    ev = regimen(lp=(0.0,))
    rec, _ = simulate(p, yp, ev, 0.30, dt=1.0e-5, out_every=0.005)
    say("%-10s %10s %10s %10s" % ("min", "ICP(cm)", "grad", "TLPD"))
    for t, y in rec[::4]:
        q = outputs(p, y)
        say("%-10.1f %10.1f %10.2f %10.2f" % (t * 1440.0, q["ICPCM"],
                                              q["GRAD"], q["TLPD"]))
    say("")
    say("(b) the same LP followed for 30 days, with the headache states")
    rec, _ = simulate(p, yp, ev, 30.0, dt=0.002, out_every=1.0)
    say("%-8s %10s %10s %10s %10s" % ("day", "ICP(cm)", "HA d/m", "STG", "CGRP"))
    for t, y in rec[::3]:
        q = outputs(p, y)
        say("%-8.0f %10.1f %10.1f %10.3f %10.3f" %
            (t, q["ICPCM"], q["MHD"], y[IX["STG"]], y[IX["CGRP"]]))
    say("")
    say("Note the shape.  tau printed above is the relaxation time AT the")
    say("presenting pressure; because compliance is 1/(E1*ICP) it lengthens")
    say("as the pressure falls, so recovery is not a single exponential: the")
    say("pressure undershoots to a few cmH2O within half an hour and is back")
    say("within 1 cmH2O of baseline at about two hours.  Over 30 days the")
    say("headache states do not move at all.  Whatever relief a diagnostic")
    say("LP provides for days is therefore not a pressure effect in this")
    say("model, and 'the LP helped' is not evidence that pressure is the")
    say("driver -- which is what the therapeutic-LP literature reports.")
    say("")
    say("(c) a single opening pressure is a sample of a moving quantity.")
    say("    Posture, CO2 and the cardiac cycle each move it; below, the")
    say("    diagnostic threshold of 25 cmH2O is applied to patients whose")
    say("    TRUE mean ICP is known.")
    say("")
    say("%-12s %12s %12s %14s %14s" %
        ("true(cm)", "pulse amp", "meas 5th", "meas 95th", "P(<25 cmH2O)"))
    for bmi in (28, 30, 32, 34, 36, 38, 42):
        r = ss(p, bmi)
        icp = r["ICP"]
        # pulse amplitude on the exponential P-V curve: dP = E1*P*dV_pulse
        amp = p["E1"] * icp * 0.55
        # measurement noise: posture/CO2/anxiety, sd 2.2 mmHg, plus pulse
        sd = math.sqrt(2.2 ** 2 + (amp / 2.0) ** 2)
        lo, hi = icp - 1.645 * sd, icp + 1.645 * sd
        z = (25.0 * MMHG_PER_CMH2O - icp) / sd
        pr = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        say("%-12.1f %12.2f %12.1f %14.1f %14.2f" %
            (cm(icp), amp, cm(lo), cm(hi), pr))


# --------------------------------------------------------------------------
# ANALYSIS 10 -- headache does not read ICP
# --------------------------------------------------------------------------

def analysis10():
    head(10, "why two thirds of patients still have headache after the"
             " pressure is fixed")
    p = params()
    say("%-42s %9s %9s %8s %8s %8s" %
        ("", "ICP(cm)", "HA d/m", "STG", "MOH", "tinnitus"))
    for lab, analg, ev in (
            ("presentation, low analgesic use", 0.15, None),
            ("presentation, high analgesic use", 0.85, None),
            ("shunt, low analgesic use", 0.15, regimen(shunt_day=0.0)),
            ("shunt, high analgesic use", 0.85, regimen(shunt_day=0.0)),
            ("shunt + analgesic withdrawal at day 90", 0.85, "withdraw"),
            ("stent, high analgesic use", 0.85, regimen(stent_day=0.0)),
            ("ACZ 2 g/d, high analgesic use", 0.85, regimen(acz=2000))):
        q = params()
        y = presentation(q, 38.0, 180.0, ev=regimen(analg=analg))
        if ev is None:
            r = outputs(q, y)
            say("%-42s %9.1f %9.1f %8.3f %8.3f %8.2f" %
                (lab, r["ICPCM"], r["MHD"], y[IX["STG"]], y[IX["MOH"]],
                 r["TINN"]))
            continue
        if ev == "withdraw":
            e1 = regimen(shunt_day=0.0, analg=0.85)
            _, y = simulate(q, y, e1, 90.0, dt=0.002, out_every=None)
            e2 = regimen(analg=0.05)
            _, y = simulate(q, y, e2, 180.0, dt=0.002, out_every=None)
        else:
            e = dict(ev); e["analg"] = analg
            _, y = simulate(q, y, e, 270.0, dt=0.002, out_every=None)
        r = outputs(q, y)
        say("%-42s %9.1f %9.1f %8.3f %8.3f %8.2f" %
            (lab, r["ICPCM"], r["MHD"], y[IX["STG"]], y[IX["MOH"]], r["TINN"]))
    say("")
    say("The sensitisation state STG has its own hysteresis: medication")
    say("overuse blocks its decay (FMOHDES), so the SAME normalised pressure")
    say("leaves a different headache burden depending on a variable no")
    say("pressure measurement contains.")


# --------------------------------------------------------------------------
# ANALYSIS 11 -- a virtual cohort: which patients does each policy help?
# --------------------------------------------------------------------------

def analysis11(n=48, seed=20260728):
    head(11, "virtual cohort: allocation, not efficacy")
    rnd = random.Random(seed)
    say("%d virtual patients.  Randomised: sinus wall support PCRIT" % n)
    say("(N(11.0, 0.9)), BMI (N(38, 5)), analgesic use U(0.1, 0.9) and")
    say("diagnostic delay (exponential, mean 70 d, capped at 300 d).")
    say("Each patient is simulated under five policies to month 14.")
    say("")
    pol = [("no therapy", lambda d: regimen(start=d)),
           ("ACZ 2 g/d", lambda d: regimen(acz=2000, start=d)),
           ("weight programme -12 %",
            lambda d: regimen(diet=0.12, diet_start=d, start=d)),
           ("ACZ + weight programme",
            lambda d: regimen(acz=2000, diet=0.12, diet_start=d, start=d)),
           ("stent if gradient > 8 mmHg",
            lambda d: regimen(stent_day=d, start=d))]
    tot = {k: [] for k, _ in pol}
    gsave = []
    for i in range(n):
        pc = max(9.15, min(13.5, rnd.gauss(11.0, 0.9)))
        bmi = max(27.0, min(52.0, rnd.gauss(38.0, 5.0)))
        analg = rnd.uniform(0.1, 0.9)
        delay = min(300.0, rnd.expovariate(1.0 / 70.0))
        q = params(PCRIT=pc)
        y0 = init_state(q, bmi=bmi)
        ev0 = regimen(analg=analg)
        _, ypres = simulate(q, y0, ev0, max(1.0, delay), dt=0.006,
                            out_every=None)
        b = outputs(q, ypres)
        gsave.append((b["G"], b["GRAD"], b["ICPCM"], delay, bmi, pc))
        for lab, mk in pol:
            e = mk(0.0); e["analg"] = analg
            if lab.startswith("stent") and b["GRAD"] <= 8.0:
                e = regimen(analg=analg)      # not eligible
            _, y = simulate(q, ypres, e, 425.0 - max(1.0, delay), dt=0.006,
                            out_every=None)
            r = outputs(q, y)
            tot[lab].append((r["MD"] - b["MD"], r["ICPCM"] - b["ICPCM"],
                             r["MHD"] - b["MHD"], r["AXON"], r["MD"],
                             y[IX["CUMEXP"]], r["ICPCM"]))
    def mean(v):
        return sum(v) / len(v)
    say("%-30s %9s %9s %9s %9s %9s" %
        ("policy", "final MD", "dMD(dB)", "final ICP", "dHA d/m", "axons"))
    for lab, _ in pol:
        v = tot[lab]
        say("%-30s %9.2f %9.2f %9.1f %9.2f %9.1f" %
            (lab, mean([x[4] for x in v]), mean([x[0] for x in v]),
             mean([x[6] for x in v]), mean([x[2] for x in v]),
             mean([x[3] for x in v])))
    say("")
    say("stratified by baseline loop gain g (ACZ + weight programme):")
    say("%-18s %6s %10s %10s %10s" % ("g band", "n", "dMD", "dICP", "mean g"))
    bands = [(0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 0.70), (0.70, 1.01)]
    v = tot["ACZ + weight programme"]
    for lo, hi in bands:
        idx = [i for i in range(n) if lo <= gsave[i][0] < hi]
        if not idx:
            continue
        say("%-18s %6d %10.2f %10.2f %10.3f" %
            ("%.2f-%.2f" % (lo, hi), len(idx),
             mean([v[i][0] for i in idx]), mean([v[i][1] for i in idx]),
             mean([gsave[i][0] for i in idx])))
    say("")
    say("stratified by diagnostic delay (ACZ + weight programme):")
    say("    (dMD is measured from each patient's OWN presentation, which for")
    say("     a late-diagnosed patient is already a damaged baseline -- so read")
    say("     the FINAL columns, not the deltas.)")
    say("%-14s %5s %10s %10s %11s %11s" %
        ("delay (d)", "n", "final MD", "dMD", "final axons", "exposure"))
    for lo, hi in ((0, 21), (21, 60), (60, 120), (120, 300.1)):
        idx = [i for i in range(n) if lo <= gsave[i][3] < hi]
        if not idx:
            continue
        say("%-14s %5d %10.2f %10.2f %11.1f %11.0f" %
            ("%d-%d" % (lo, hi), len(idx), mean([v[i][4] for i in idx]),
             mean([v[i][0] for i in idx]), mean([v[i][3] for i in idx]),
             mean([v[i][5] for i in idx])))
    say("")
    say("is the MEASURABLE venous gradient a usable proxy for the")
    say("UNMEASURABLE loop gain that actually sets drug response?")
    say("")

    def corr(a, b):
        ma, mb = mean(a), mean(b)
        num = sum((a[i] - ma) * (b[i] - mb) for i in range(len(a)))
        den = math.sqrt(sum((x - ma) ** 2 for x in a)
                        * sum((x - mb) ** 2 for x in b))
        return num / den if den else float("nan")

    gs = [x[0] for x in gsave]
    gr = [x[1] for x in gsave]
    ic = [x[2] for x in gsave]
    say("  r(g, transverse sinus gradient) = %+.3f" % corr(gs, gr))
    say("  r(g, opening pressure)          = %+.3f" % corr(gs, ic))
    say("  r(gradient, opening pressure)   = %+.3f" % corr(gr, ic))
    say("")
    say("  The answer is NO, and the sign is the interesting part.  ANALYSIS 1")
    say("  showed that g is NOT monotone in disease severity: it peaks around")
    say("  BMI 36-38 and falls again as the collapse saturates.  In a cohort")
    say("  spanning that peak, the gradient and the gain therefore run in")

    say("  OPPOSITE directions over part of the range, and the linear")
    say("  correlation comes out negative.  The clinically attractive shortcut")
    say("  -- 'measure the gradient, infer how much the drug will do' -- is")
    say("  refuted by the same structure that makes the gradient worth")
    say("  measuring for stent selection.  Two different questions, two")
    say("  different variables.")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    say("=" * 78)
    say("IDIOPATHIC INTRACRANIAL HYPERTENSION -- QSP reference computation")
    say("44 states, RK4, dt = 0.002 day unless stated; pure Python stdlib")
    say("=" * 78)
    p = params()
    say("")
    say("default parameter set (the 'typical' susceptible patient):")
    for k in sorted(p):
        say("  %-10s = %s" % (k, p[k]))
    analysis0()
    analysis1()
    analysis2()
    analysis3()
    analysis4()
    analysis5()
    analysis6()
    analysis7()
    analysis8()
    analysis9()
    analysis10()
    analysis11()
    say()
    say("=" * 78)
    say("END")
    say("=" * 78)
    if _FH[0] is not None:
        _FH[0].close()


if __name__ == "__main__":
    main()
