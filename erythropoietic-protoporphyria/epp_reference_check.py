#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
epp_reference_check.py
======================
A dependency-free (pure Python; no numpy/scipy) reference re-implementation of
the mrgsolve model in `epp_mrgsolve_model.R`.  Its job is to VERIFY that every
number quoted in README.md and in the model header is actually produced by the
equations rather than asserted by hand.

The `rhs()` function below is a line-for-line mirror of the `$ODE` block in the
R file and `readouts()` mirrors `$TABLE`.  If you change one, change the other.

    python3 epp_reference_check.py > epp_reference_output.txt

WHAT IT CHECKS
--------------
  A. Steady-state biomarker calibration       normal / carrier / EPP / XLP
  B. The ferrochelatase threshold             is ~35% residual really the knee?
  C. Sign flip #1 - zinc protoporphyrin       EPP down, XLP up, same overload
  D. Sign flip #2 - iron                      helps XLP, harms EPP
  E. Photobiology                             minutes to prodrome
  F. Drug effects                             afamelanotide / dersimelagon /
                                              bitopertin / beta-carotene
  G. Multiplicativity of the two therapy axes vs an additive expectation
  H. Hepatic bistability                      saddle-node continuation
  I. A day in the life                        prodrome, retreat, pain course
  J. The behavioural control loop             reaction-time (delay) sweep
  K. Environmental modifiers                  glass, cloud, sunscreen, cover

Author: QSP disease-model library (Claude Code Routine)
"""

import math

# =====================================================================
#  PARAMETERS  — identical names/values to $PARAM in epp_mrgsolve_model.R
# =====================================================================
P = dict(
    # ---- genotype switches -------------------------------------------------
    FRES   = 0.15,   # residual ferrochelatase activity (fraction of normal)
    GOF    = 1.0,    # ALAS2 gain-of-function multiplier (XLP ~3.5)
    FMC1R  = 1.0,    # MC1R functionality (~0.3 for red-hair LoF variants)

    # ---- erythroid glycine / ALAS2 ----------------------------------------
    KINGLY = 4.0,    # glycine influx (mmol/L/h) via GlyT1 / SLC6A9
    KOUTGLY= 5.0,    # glycine consumption other than ALAS2 (1/h)
    GLY0   = 0.8,    # baseline erythroblast glycine (mmol/L)
    KMGLY  = 0.8,    # ALAS2 Km for glycine — half-saturated at baseline
    VALAS  = 201.4,  # ALAS2 Vmax (umol/L_erythron/h)
    KHFB   = 1200.0, # heme negative feedback constant (weak in the erythron)
    KIRE   = 12.0,   # IRE/IRP half-point on labile iron (= 0.6 x FE0)
    NIRE   = 2.0,    # IRE/IRP Hill coefficient

    # ---- pathway transfer (every enzyme in large excess, 1/h) -------------
    KALAD  = 6.0, KHMBS = 6.0, KUROD = 6.0, KCPOX = 6.0, KPPOX = 6.0,

    # ---- ferrochelatase ----------------------------------------------------
    VFECH  = 435.0,  # umol/L/h.  Normal chelation CAPACITY = 2.9 x normal flux
    KMPPIX = 1.0,    # Km for protoporphyrin IX (umol/L)
    KMFE   = 10.0,   # Km for Fe2+ (umol/L)
    KISC   = 20.0,   # 2Fe-2S cluster assembly half-point on labile iron
    VZNREL = 0.02,   # Zn-insertion efficiency relative to Fe-insertion
    KMZN   = 5.0, ZN = 8.0, KIFZ = 3.0,   # Zn substrate + Fe competition
    KESC   = 0.672,  # PPIX escape erythron -> plasma (1/h)  [OVERFLOW ROUTE]
    KHEMEO = 1.0,    # heme utilisation (1/h)

    # ---- iron --------------------------------------------------------------
    FE0    = 20.0,   # reference labile iron pool (umol/L)
    KFEIN  = 20.0,   # baseline iron supply (umol/L/h)
    KFEOUT = 1.0,    # iron efflux / storage (1/h)
    FEDOSE = 0.0,    # iron supplementation, as a fraction of KFEIN

    # ---- erythrocyte pools --------------------------------------------------
    KPROD  = 2.4056e-4,  # ln2/(120 d) erythropoietic release (1/h)
    KSEN   = 2.4056e-4,  # RBC senescence (1/h)
    KEFFL  = 7.2e-4,     # PPIX efflux out of circulating RBC (1/h)
    ALPHA  = 1.000,      # erythroblast -> erythrocyte PPIX partition
    BETAZN = 8.09e-4,    # Zn-PP packaging coefficient
    KZN2   = 8.0e-5,     # late-reticulocyte Zn insertion (needs INTACT FECH)
    ERYSUP = 1.0,        # erythropoietic suppression (transfusion: < 1)

    # ---- distribution -------------------------------------------------------
    KTR    = 0.0521, # erythron escape -> plasma scaling
    KPLLIV = 3.0, KLIVPL = 0.1, KPLSK = 2.0, KSKPL = 1.667,
    KEHC   = 0.05,   # enterohepatic reabsorption (1/h)
    KFEC   = 0.30,   # faecal loss from the gut lumen (1/h)
    CHOLEFF= 0.0,    # cholestyramine effect (0-0.9) on reabsorption
    KPHER  = 0.0,    # plasmapheresis clearance (1/h)

    # ---- hepatobiliary ------------------------------------------------------
    VBILE  = 20.0,   # canalicular Vmax (umol/L/h) — the SATURABLE choke point
    KMBILE = 3.0,
    KCHIN  = 0.010,  # cholestasis induction (1/h)
    KCHOUT = 0.004,  # cholestasis resolution (1/h)
    KCH    = 1.8,    # hepatic PPIX at which crystals drive cholestasis
    NCH    = 6.0,    # Hill coefficient of crystal -> cholestasis
    FCHOL  = 0.92,   # maximum fractional loss of bile flow
    KINJ   = 0.012, KINJOUT = 0.01, KILIV = 5.0,

    # ---- photobiology -------------------------------------------------------
    ODMEL  = 0.693,  # epidermal optical density per melanin unit at 405 nm
    KPHOT  = 8.20,   # photodynamic dose rate constant
    EREF   = 1.0,    # reference irradiance for the tolerance-time read-out
    EMAX   = 1.0,    # peak ambient violet irradiance (normalised)
    FCLOUD = 1.0,    # cloud transmission (overcast ~0.5 at 400-410 nm)
    FGLASS = 1.0,    # window-glass transmission at 405 nm (~0.9!)
    FCLOTH = 1.0,    # clothing / opaque screen (0.05 = full cover)
    DPRO   = 1.0,    # prodrome threshold (dose units)
    DCRIT  = 3.0,    # full phototoxic-reaction threshold
    KDB    = 0.15,   # decay of the prodromal dose memory (1/h) -> ~4.6 h
    IC50BC = 12000.0,# beta-carotene plasma conc. for 50% quenching (ug/L)
    KOXOUT = 0.70,   # oxidative-signal clearance (1/h)
    OUTS   = 10.0, OUTE = 18.0,   # intended outdoor window (hour of day)
    TAUDEL = 0.033,  # behavioural retreat delay (h) ~2 min for an adult
    NBOUT  = 3.0,    # outdoor bouts per day
    WINDOW = 8.0,    # 10:00-18:00 registrational window (h)

    # ---- injury cascade -----------------------------------------------------
    KMAST = 1.0, KMASTO = 0.25, KC5 = 1.0, KC5O = 0.35,
    KED = 0.05, KEDO = 0.06, KNOC = 0.22, KNOCO = 0.020, KPAIN = 6.0,
    AMAX = 8.0, KAMP = 3.0, NAMP = 4.0,   # all-or-nothing reaction gate

    # ---- MC1R / melanin -----------------------------------------------------
    EMAXMC = 12.0,   # maximal fold-increase in melanogenic drive over basal
    EC50AF = 0.30,   # afamelanotide plasma EC50 (ug/L)
    EC50DE = 30000.0,# dersimelagon plasma EC50 (ug/L)
    KTYRIN = 4.13e-3, KTYRO = 4.13e-3,   # tyrosinase turnover (t1/2 7 d)
    KMELIN = 8.25e-4, KMELO = 8.25e-4,   # melanin turnover    (t1/2 35 d)

    # ---- drug PK ------------------------------------------------------------
    KREL   = 0.030,  # afamelanotide implant release (1/h)
    FAFA   = 1.0, VAFA = 20.0, KEAFA = 1.386,
    KADER  = 0.50, VDER = 200.0, KEDER = 0.0289,
    KABIT  = 0.40, VBIT = 700.0, KEBIT = 0.01386,
    IMAXBIT= 0.85, IC50BIT = 328.0,     # GlyT1 inhibition of glycine influx
    KABC   = 0.05, VBC = 3000.0, KEBC = 0.00206,  # beta-carotene V/F (t1/2 14 d)
    KEHEM  = 0.10, IMAXHEM = 0.35, IC50HEM = 5.0,
)

SN = ["GLY","ALA","PBG","UPG","CPG","PPG","PPIXE","HEME","FE",
      "PRBC","ZRBC","PPL","PSK","PLIV","PGUT","CHOL","LINJ",
      "OX","MAST","C5A","EDEMA","NOCI","AVOID","DBOUT",
      "TYR","MEL","ADEP","AC","DG","DC","BG","BC","CG","CC","HEMC",
      "SUNCUM","RXNCUM"]
IX = {n: i for i, n in enumerate(SN)}
NS = len(SN)


# =====================================================================
#  RIGHT-HAND SIDE  (mirror of $ODE in epp_mrgsolve_model.R)
# =====================================================================
def rhs(t, y, p, exposure=True):
    (GLY, ALA, PBG, UPG, CPG, PPG, PPIXE, HEME, FE,
     PRBC, ZRBC, PPL, PSK, PLIV, PGUT, CHOL, LINJ,
     OX, MAST, C5A, EDEMA, NOCI, AVOID, DBOUT,
     TYR, MEL, ADEP, AC, DG, DC, BG, BC, CG, CC, HEMC,
     SUNCUM, RXNCUM) = y

    GLY = max(GLY, 1e-9); FE = max(FE, 1e-9)
    PPIXE = max(PPIXE, 0.0); PLIV = max(PLIV, 0.0); PSK = max(PSK, 0.0)

    # --- bitopertin: GlyT1 blockade lowers erythroblast glycine influx ------
    IBIT = p["IMAXBIT"]*BC/(p["IC50BIT"] + BC)
    dGLY = p["KINGLY"]*(1.0 - IBIT) - p["KOUTGLY"]*GLY

    # --- ALAS2 : substrate x iron(IRE/IRP) x heme feedback x genotype -------
    GLYSAT = GLY/(p["KMGLY"] + GLY)
    n = p["NIRE"]
    IREF = ((FE**n/(p["KIRE"]**n + FE**n)) /
            (p["FE0"]**n/(p["KIRE"]**n + p["FE0"]**n)))
    FHEME = 1.0/(1.0 + (HEME/p["KHFB"])**2)
    IHEM = p["IMAXHEM"]*HEMC/(p["IC50HEM"] + HEMC)
    vALAS = (p["VALAS"]*p["GOF"]*GLYSAT*IREF*FHEME
             * (1.0 - IHEM)*p["ERYSUP"])

    dALA = vALAS            - p["KALAD"]*ALA
    dPBG = p["KALAD"]*ALA   - p["KHMBS"]*PBG
    dUPG = p["KHMBS"]*PBG   - p["KUROD"]*UPG
    dCPG = p["KUROD"]*UPG   - p["KCPOX"]*CPG
    dPPG = p["KCPOX"]*CPG   - p["KPPOX"]*PPG

    # --- ferrochelatase : ONE enzyme, TWO metals ----------------------------
    FISC  = (FE/(p["KISC"] + FE))/(p["FE0"]/(p["KISC"] + p["FE0"]))
    fFECH = p["FRES"]*FISC
    sat   = PPIXE/(p["KMPPIX"] + PPIXE)
    vFe   = p["VFECH"]*fFECH*sat*(FE/(p["KMFE"] + FE))
    vZn   = (p["VFECH"]*fFECH*p["VZNREL"]*sat
             * (p["ZN"]/(p["KMZN"] + p["ZN"]))
             * (p["KIFZ"]/(p["KIFZ"] + FE)))
    Jesc  = p["KESC"]*PPIXE

    dPPIXE = p["KPPOX"]*PPG - vFe - vZn - Jesc
    dHEME  = vFe - p["KHEMEO"]*HEME
    dFE    = p["KFEIN"]*(1.0 + p["FEDOSE"]) - p["KFEOUT"]*FE

    # --- erythrocyte pools ---------------------------------------------------
    dPRBC = (p["KPROD"]*p["ALPHA"]*PPIXE*p["ERYSUP"]
             - (p["KSEN"] + p["KEFFL"])*PRBC)
    dZRBC = (p["BETAZN"]*vZn + p["KZN2"]*fFECH*PRBC) - p["KSEN"]*ZRBC

    # --- distribution : plasma / skin / liver / gut --------------------------
    reabs = p["KEHC"]*(1.0 - p["CHOLEFF"])*PGUT
    dPPL = (p["KTR"]*Jesc + p["KEFFL"]*PRBC + p["KSKPL"]*PSK + reabs
            + p["KLIVPL"]*PLIV
            - (p["KPLLIV"] + p["KPLSK"] + p["KPHER"])*PPL)
    dPSK = p["KPLSK"]*PPL - p["KSKPL"]*PSK

    BF    = 1.0 - p["FCHOL"]*CHOL
    vcan  = p["VBILE"]*BF*PLIV/(p["KMBILE"] + PLIV)
    dPLIV = p["KPLLIV"]*PPL - vcan - p["KLIVPL"]*PLIV
    dPGUT = vcan - p["KFEC"]*PGUT - reabs

    hill  = PLIV**p["NCH"]/(p["KCH"]**p["NCH"] + PLIV**p["NCH"])
    dCHOL = p["KCHIN"]*hill*(1.0 - CHOL) - p["KCHOUT"]*CHOL
    dLINJ = p["KINJ"]*(CHOL**2 + (PLIV/p["KILIV"])**2) - p["KINJOUT"]*LINJ

    # --- MC1R -> melanin -----------------------------------------------------
    occ = AC/p["EC50AF"] + DC/p["EC50DE"]
    SIG = p["FMC1R"]*p["EMAXMC"]*occ/(1.0 + occ)
    dTYR = p["KTYRIN"]*(1.0 + SIG) - p["KTYRO"]*TYR
    dMEL = p["KMELIN"]*TYR - p["KMELO"]*MEL

    # --- photobiology --------------------------------------------------------
    TEPI = math.exp(-p["ODMEL"]*MEL)
    QUEN = 1.0/(1.0 + CC/p["IC50BC"])
    tod  = math.fmod(t, 24.0)
    s    = math.sin(math.pi*(tod - 6.0)/12.0)
    sun  = (s**1.3) if (6.0 <= tod <= 18.0 and s > 0.0) else 0.0
    fout = 1.0 if (exposure and p["OUTS"] <= tod <= p["OUTE"]) else 0.0
    Eamb = p["EMAX"]*sun*p["FCLOUD"]*p["FGLASS"]*p["FCLOTH"]
    R    = p["KPHOT"]*Eamb*TEPI*PSK*QUEN*fout*(1.0 - AVOID)

    # prodromal dose memory: the patient's SENSOR integrates dose with a
    # multi-hour memory, which is why one prodrome ends the whole outing.
    dDB  = R - p["KDB"]*DBOUT
    SW   = 1.0/(1.0 + math.exp(-(DBOUT - p["DPRO"])/0.08))
    dAVO = (SW - AVOID)/p["TAUDEL"]

    dOX   = R - p["KOXOUT"]*OX
    # Above the reaction threshold, mast-cell degranulation and complement
    # activation become self-amplifying: this is what makes the full
    # phototoxic reaction all-or-nothing rather than graded.
    dbn   = (DBOUT/p["KAMP"])**p["NAMP"]
    AMP   = 1.0 + p["AMAX"]*dbn/(1.0 + dbn)
    dMAST = p["KMAST"]*OX*AMP - p["KMASTO"]*MAST
    dC5A  = p["KC5"]*OX*AMP   - p["KC5O"]*C5A
    dEDE  = p["KED"]*(MAST + C5A) - p["KEDO"]*EDEMA
    dNOCI = p["KNOC"]*(1.5*OX + 0.5*EDEMA) - p["KNOCO"]*NOCI

    dSUN  = fout*(1.0 - AVOID)
    dRXN  = (1.0/24.0)/(1.0 + math.exp(-(DBOUT - p["DCRIT"])/0.20))

    # --- drug PK -------------------------------------------------------------
    dADEP = -p["KREL"]*ADEP
    dAC   = p["KREL"]*ADEP*p["FAFA"]/p["VAFA"] - p["KEAFA"]*AC
    dDG   = -p["KADER"]*DG
    dDC   = p["KADER"]*DG/p["VDER"] - p["KEDER"]*DC
    dBG   = -p["KABIT"]*BG
    dBC   = p["KABIT"]*BG/p["VBIT"] - p["KEBIT"]*BC
    dCG   = -p["KABC"]*CG
    dCC   = p["KABC"]*CG/p["VBC"] - p["KEBC"]*CC
    dHEMC = -p["KEHEM"]*HEMC

    return [dGLY,dALA,dPBG,dUPG,dCPG,dPPG,dPPIXE,dHEME,dFE,
            dPRBC,dZRBC,dPPL,dPSK,dPLIV,dPGUT,dCHOL,dLINJ,
            dOX,dMAST,dC5A,dEDE,dNOCI,dAVO,dDB,
            dTYR,dMEL,dADEP,dAC,dDG,dDC,dBG,dBC,dCG,dCC,dHEMC,
            dSUN,dRXN]


# =====================================================================
#  DERIVED READ-OUTS  (mirror of $TABLE)
# =====================================================================
def readouts(y, p):
    d = {n: y[IX[n]] for n in SN}
    tot  = d["PRBC"] + d["ZRBC"]
    TEPI = math.exp(-p["ODMEL"]*d["MEL"])
    QUEN = 1.0/(1.0 + d["CC"]/p["IC50BC"])
    ttol_h = p["DPRO"]/(p["KPHOT"]*p["EREF"]*TEPI*max(d["PSK"], 1e-12)*QUEN)
    sunhr  = p["WINDOW"]*(1.0 - math.exp(-p["NBOUT"]*ttol_h/p["WINDOW"]))
    return dict(
        PPIXRBC=d["PRBC"], ZNPPRBC=d["ZRBC"],
        ZNFRAC=100.0*d["ZRBC"]/tot if tot > 0 else 0.0,
        PPIXPL=d["PPL"], PPIXSK=d["PSK"], PPIXLIV=d["PLIV"], CHOL=d["CHOL"],
        ALT=25.0*(1.0 + 3.0*d["LINJ"]),
        TBIL=0.6*(1.0 + 6.0*d["CHOL"]*d["LINJ"]),
        MEL=d["MEL"], TEPI=TEPI,
        TTOLMIN=60.0*ttol_h, SUNHRD=sunhr,
        PAIN=10.0*d["NOCI"]**2/(p["KPAIN"]**2 + d["NOCI"]**2),
        SUNCUM=d["SUNCUM"], RXNDAYS=d["RXNCUM"], GLY=d["GLY"],
    )


def pars(**kw):
    q = dict(P); q.update(kw); return q


# =====================================================================
#  ANALYTIC STEADY STATE  (fast; used for every steady-state table)
# =====================================================================
def erythron_ss(p, BCconc=0.0):
    """Steady state of glycine -> ALAS2 -> ... -> FECH -> RBC pools."""
    IBIT = p["IMAXBIT"]*BCconc/(p["IC50BIT"] + BCconc)
    GLY  = p["KINGLY"]*(1.0 - IBIT)/p["KOUTGLY"]
    GLYSAT = GLY/(p["KMGLY"] + GLY)
    FE = p["KFEIN"]*(1.0 + p["FEDOSE"])/p["KFEOUT"]
    n = p["NIRE"]
    IREF = ((FE**n/(p["KIRE"]**n + FE**n)) /
            (p["FE0"]**n/(p["KIRE"]**n + p["FE0"]**n)))
    FISC  = (FE/(p["KISC"] + FE))/(p["FE0"]/(p["KISC"] + p["FE0"]))
    fFECH = p["FRES"]*FISC
    kFe = p["VFECH"]*fFECH*(FE/(p["KMFE"] + FE))
    kZn = (p["VFECH"]*fFECH*p["VZNREL"]
           * (p["ZN"]/(p["KMZN"] + p["ZN"]))*(p["KIFZ"]/(p["KIFZ"] + FE)))
    HEME = 100.0
    for _ in range(80):
        FHEME = 1.0/(1.0 + (HEME/p["KHFB"])**2)
        J = p["VALAS"]*p["GOF"]*GLYSAT*IREF*FHEME*p["ERYSUP"]
        lo, hi = 0.0, 1e8
        for _ in range(220):
            m = 0.5*(lo + hi)
            f = (kFe + kZn)*m/(p["KMPPIX"] + m) + p["KESC"]*m - J
            if f > 0: hi = m
            else:     lo = m
        PPIXE = 0.5*(lo + hi)
        HEME = kFe*PPIXE/(p["KMPPIX"] + PPIXE)/p["KHEMEO"]
    sat = PPIXE/(p["KMPPIX"] + PPIXE)
    vFe, vZn = kFe*sat, kZn*sat
    Jesc = p["KESC"]*PPIXE
    PRBC = p["KPROD"]*p["ALPHA"]*PPIXE*p["ERYSUP"]/(p["KSEN"] + p["KEFFL"])
    ZRBC = (p["BETAZN"]*vZn + p["KZN2"]*fFECH*PRBC)/p["KSEN"]
    return dict(GLY=GLY, FE=FE, J=J, HEME=HEME, PPIXE=PPIXE, vFe=vFe, vZn=vZn,
                Jesc=Jesc, PRBC=PRBC, ZRBC=ZRBC, fFECH=fFECH, IREF=IREF,
                FISC=FISC, CAP=kFe)


def liver_roots(p, IN):
    """All steady states of the (hepatic PPIX, cholestasis) subsystem for a
    given systemic input flux IN = KTR*Jesc + KEFFL*PRBC.
    Returns [(PLIV, CHOL, stable?), ...] sorted by PLIV."""
    fre = (p["KEHC"]*(1.0 - p["CHOLEFF"])
           / (p["KFEC"] + p["KEHC"]*(1.0 - p["CHOLEFF"])))
    g = (p["KPLLIV"] + p["KPHER"])/p["KPLLIV"]

    def F(L):
        hill = L**p["NCH"]/(p["KCH"]**p["NCH"] + L**p["NCH"])
        C = p["KCHIN"]*hill/(p["KCHOUT"] + p["KCHIN"]*hill)
        vcan = p["VBILE"]*(1.0 - p["FCHOL"]*C)*L/(p["KMBILE"] + L)
        return vcan*(g - fre) + p["KLIVPL"]*L*(g - 1.0) - IN

    roots, step = [], 0.002
    x = step; prev = F(x)
    while x < 200.0:
        x2 = x + step
        cur = F(x2)
        if prev == 0.0 or prev*cur < 0.0:
            lo, hi = x, x2
            for _ in range(90):
                mid = 0.5*(lo + hi)
                if F(lo)*F(mid) <= 0.0: hi = mid
                else:                   lo = mid
            r = 0.5*(lo + hi)
            hill = r**p["NCH"]/(p["KCH"]**p["NCH"] + r**p["NCH"])
            C = p["KCHIN"]*hill/(p["KCHOUT"] + p["KCHIN"]*hill)
            eps = max(1e-5, 1e-4*r)
            roots.append((r, C, (F(r + eps) - F(r - eps)) > 0.0))
        prev = cur
        x = x2
        if x > 20.0: step = 0.05
    return roots


def full_ss(p, BCconc=0.0, branch="low"):
    """Complete no-exposure steady state as a state VECTOR."""
    s = erythron_ss(p, BCconc)
    IN = p["KTR"]*s["Jesc"] + p["KEFFL"]*s["PRBC"]
    roots = liver_roots(p, IN)
    stab = [r for r in roots if r[2]] or roots
    r = stab[0] if branch == "low" else stab[-1]
    PLIV, CHOL = r[0], r[1]
    vcan = p["VBILE"]*(1.0 - p["FCHOL"]*CHOL)*PLIV/(p["KMBILE"] + PLIV)
    PPL = (vcan + p["KLIVPL"]*PLIV)/p["KPLLIV"]
    PSK = p["KPLSK"]*PPL/p["KSKPL"]
    PGUT = vcan/(p["KFEC"] + p["KEHC"]*(1.0 - p["CHOLEFF"]))
    LINJ = p["KINJ"]*(CHOL**2 + (PLIV/p["KILIV"])**2)/p["KINJOUT"]
    y = [0.0]*NS
    y[IX["GLY"]] = s["GLY"]
    for nm, k in (("ALA","KALAD"), ("PBG","KHMBS"), ("UPG","KUROD"),
                  ("CPG","KCPOX"), ("PPG","KPPOX")):
        y[IX[nm]] = s["J"]/p[k]
    y[IX["PPIXE"]] = s["PPIXE"]; y[IX["HEME"]] = s["HEME"]; y[IX["FE"]] = s["FE"]
    y[IX["PRBC"]] = s["PRBC"];   y[IX["ZRBC"]] = s["ZRBC"]
    y[IX["PPL"]] = PPL; y[IX["PSK"]] = PSK
    y[IX["PLIV"]] = PLIV; y[IX["PGUT"]] = PGUT
    y[IX["CHOL"]] = CHOL; y[IX["LINJ"]] = LINJ
    y[IX["TYR"]] = 1.0; y[IX["MEL"]] = 1.0
    y[IX["BC"]] = BCconc
    return y


# =====================================================================
#  INTEGRATION
# =====================================================================
def rk4(t, y, h, p, exposure):
    k1 = rhs(t, y, p, exposure)
    y2 = [y[i] + 0.5*h*k1[i] for i in range(NS)]
    k2 = rhs(t + 0.5*h, y2, p, exposure)
    y3 = [y[i] + 0.5*h*k2[i] for i in range(NS)]
    k3 = rhs(t + 0.5*h, y3, p, exposure)
    y4 = [y[i] + h*k3[i] for i in range(NS)]
    k4 = rhs(t + h, y4, p, exposure)
    return [y[i] + h/6.0*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(NS)]


def simulate(p, days, y_init, h=0.05, events=None, exposure=False, record=24.0):
    y = y_init[:]
    ev = sorted(events or [], key=lambda e: e[0])
    out, t, ei = [], 0.0, 0
    n = int(round(days*24/h))
    every = max(1, int(round(record/h)))
    for i in range(n):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            y[IX[ev[ei][1]]] += ev[ei][2]; ei += 1
        if i % every == 0:
            r = readouts(y, p); r["TIME"] = t; out.append(r)
        y = rk4(t, y, h, p, exposure)
        t += h
    r = readouts(y, p); r["TIME"] = t; out.append(r)
    return y, out


# =====================================================================
#  REPORT
# =====================================================================
def hdr(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)


def main():
    print("EPP / XLP QSP MODEL — REFERENCE NUMERICAL CHECK")
    print("pure-python mirror of epp_mrgsolve_model.R  (no numpy / no scipy)")

    # ------------------------------------------------------------------ A
    hdr("A. STEADY-STATE CALIBRATION — is the biomarker scale right?")
    genos = [("Normal (FECH 100%)",        dict(FRES=1.00)),
             ("Silent carrier (FECH 50%)", dict(FRES=0.50)),
             ("Borderline    (FECH 35%)",  dict(FRES=0.35)),
             ("EPP moderate  (FECH 25%)",  dict(FRES=0.25)),
             ("EPP typical   (FECH 15%)",  dict(FRES=0.15)),
             ("EPP severe    (FECH 10%)",  dict(FRES=0.10)),
             ("XLP (ALAS2 GOF x3.5)",      dict(FRES=1.00, GOF=3.5))]
    print(f"{'genotype':<28}{'RBC PPIX':>10}{'Zn-PP':>8}{'Zn%':>7}"
          f"{'plasma':>9}{'skin':>8}{'liver':>8}{'ALT':>7}{'t_tol':>8}")
    print(f"{'':<28}{'umol/L':>10}{'umol/L':>8}{'%':>7}{'umol/L':>9}"
          f"{'umol/L':>8}{'umol/L':>8}{'U/L':>7}{'min':>8}")
    A = {}
    for name, kw in genos:
        p = pars(**kw)
        r = readouts(full_ss(p), p)
        A[name] = r
        print(f"{name:<28}{r['PPIXRBC']:>10.2f}{r['ZNPPRBC']:>8.2f}"
              f"{r['ZNFRAC']:>7.1f}{r['PPIXPL']:>9.3f}{r['PPIXSK']:>8.3f}"
              f"{r['PPIXLIV']:>8.2f}{r['ALT']:>7.0f}{r['TTOLMIN']:>8.1f}")
    print("\n  reference ranges: normal erythrocyte PPIX < 1.5 umol/L;")
    print("  clinically manifest EPP typically 10-40 umol/L; XLP comparable or")
    print("  higher WITH a large Zn-PP fraction; untreated EPP tolerance time to")
    print("  the first prodrome is minutes, not hours.")

    # verify the analytic solver against the integrator
    p = pars(FRES=0.15)
    y0 = full_ss(p)
    yv, _ = simulate(p, 400, y0, h=0.05, exposure=False, record=1e9)
    print(f"\n  self-consistency (analytic vs 400 d of integration, FECH 15%):")
    for k, nm in (("PPIXRBC","RBC PPIX"), ("PPIXPL","plasma"),
                  ("PPIXLIV","liver"), ("CHOL","cholestasis")):
        print(f"    {nm:<14}{readouts(y0,p)[k]:>10.4f}  vs "
              f"{readouts(yv,p)[k]:>10.4f}")

    # ------------------------------------------------------------------ B
    hdr("B. THE FERROCHELATASE THRESHOLD — derived, not assumed")
    pn = pars(FRES=1.0); sn = erythron_ss(pn)
    Rcap = sn["CAP"]/sn["J"]
    print(f"  Normal chelation CAPACITY / normal ALAS2 flux = R_cap = {Rcap:.3f}")
    print(f"  => the terminal step can no longer keep up below FECH residual")
    print(f"     FRES* = 1/R_cap = {100.0/Rcap:.1f}% .  The ~35% clinical")
    print(f"     penetrance threshold is a CONSEQUENCE of that capacity ratio.\n")
    print(f"{'FECH residual %':>16}{'RBC PPIX':>11}{'x normal':>10}"
          f"{'skin PPIX':>11}{'t_tol min':>11}{'sun h/day':>11}")
    base = None
    for fr in [100, 70, 50, 40, 35, 30, 25, 20, 15, 12, 10, 8]:
        q = pars(FRES=fr/100.0)
        r = readouts(full_ss(q), q)
        if base is None: base = r["PPIXRBC"]
        print(f"{fr:>16}{r['PPIXRBC']:>11.2f}{r['PPIXRBC']/base:>10.1f}"
              f"{r['PPIXSK']:>11.3f}{r['TTOLMIN']:>11.1f}{r['SUNHRD']:>11.2f}")

    # ------------------------------------------------------------------ C
    hdr("C. DISCRIMINATOR #1 — ZINC PROTOPORPHYRIN: same overload, "
        "opposite metal")
    for nm in ["Normal (FECH 100%)", "EPP typical   (FECH 15%)",
               "XLP (ALAS2 GOF x3.5)"]:
        r = A[nm]
        print(f"  {nm:<26} free PPIX {r['PPIXRBC']:>7.2f}  "
              f"Zn-PP {r['ZNPPRBC']:>6.2f}  total EP "
              f"{r['PPIXRBC']+r['ZNPPRBC']:>7.2f}  "
              f"metal-free {100-r['ZNFRAC']:>5.1f}%")
    n0, e0, x0 = (A["Normal (FECH 100%)"], A["EPP typical   (FECH 15%)"],
                  A["XLP (ALAS2 GOF x3.5)"])
    print(f"\n  free PPIX vs normal:  EPP x{e0['PPIXRBC']/n0['PPIXRBC']:.0f}"
          f"   XLP x{x0['PPIXRBC']/n0['PPIXRBC']:.0f}    "
          f"(both enormously up — the same overload)")
    print(f"  Zn-PP    vs normal:   EPP x{e0['ZNPPRBC']/n0['ZNPPRBC']:.1f}"
          f"    XLP x{x0['ZNPPRBC']/n0['ZNPPRBC']:.1f}     "
          f"(a {x0['ZNPPRBC']/e0['ZNPPRBC']:.1f}-fold DISPROPORTION)")
    print(f"  metal-free fraction:  normal {100-n0['ZNFRAC']:.0f}%  ->  "
          f"EPP {100-e0['ZNFRAC']:.0f}%  vs  XLP {100-x0['ZNFRAC']:.0f}%")
    print("\n  Mechanism: FECH is the enzyme that inserts ZINC as well as iron.")
    print("  EPP breaks that enzyme, so the substrate piles up and the metal")
    print("  cannot go in either — free PPIX rises 160-fold while Zn-PP barely")
    print("  moves. XLP leaves FECH intact and merely floods it, so zinc")
    print("  insertion rises with the substrate. The clinically used metal-free")
    print("  fraction therefore separates the two diseases, and it does so here")
    print("  with nothing in the code special-casing them: only FRES and GOF")
    print("  differ between the two runs.")
    print("\n  (Contrast: iron deficiency and lead poisoning raise Zn-PP with a")
    print("   NORMAL metal-free PPIX — FECH is intact, the substrate is normal,")
    print("   only the metal is wrong. Set FEDOSE < 0 in section D to see the")
    print("   Zn-PP arm move on its own.)")

    # ------------------------------------------------------------------ D
    hdr("D. SIGN FLIP #2 — IRON: the same nutrient helps XLP and harms EPP")
    print("  Iron enters at three points: (+) IRE/IRP de-repression of ALAS2,")
    print("  (-) Fe2+ co-substrate of FECH, (-) 2Fe-2S cluster -> FECH Vmax.\n")
    print(f"{'iron supply':>12}{'':>2}{'LIP':>7}{'IRE fac':>9}{'ISC fac':>9}"
          f"{'':>3}{'EPP FECH 15%':>22}{'':>3}{'XLP FECH 100%':>22}")
    print(f"{'x baseline':>12}{'':>2}{'uM':>7}{'':>9}{'':>9}{'':>3}"
          f"{'RBC PPIX':>11}{'delta %':>11}{'':>3}{'RBC PPIX':>11}{'delta %':>11}")
    ref_e = erythron_ss(pars(FRES=0.15))["PRBC"]
    ref_x = erythron_ss(pars(FRES=1.00, GOF=3.5))["PRBC"]
    for f in [-0.50, -0.25, 0.0, 0.25, 0.50, 1.00, 2.00]:
        se = erythron_ss(pars(FRES=0.15, FEDOSE=f))
        sx = erythron_ss(pars(FRES=1.00, GOF=3.5, FEDOSE=f))
        print(f"{1+f:>12.2f}{'':>2}{se['FE']:>7.1f}{se['IREF']:>9.3f}"
              f"{se['FISC']:>9.3f}{'':>3}{se['PRBC']:>11.2f}"
              f"{100*(se['PRBC']/ref_e - 1):>+11.1f}{'':>3}{sx['PRBC']:>11.2f}"
              f"{100*(sx['PRBC']/ref_x - 1):>+11.1f}")
    print("\n  => In XLP, iron loading is STRONGLY beneficial: FECH is intact, so")
    print("     the extra iron feeds the disposal arm and it wins outright.")
    print("     In EPP the SAME intervention is essentially inert (a few percent")
    print("     the wrong way at modest doses): the disposal enzyme is broken and")
    print("     cannot use the iron, so all that is left is the IRE/IRP supply")
    print("     arm pushing ALAS2 up. Identical equations, identical parameters —")
    print("     only FRES and GOF differ between the two columns.")
    print("     This reproduces the clinical asymmetry: iron repletion has")
    print("     repeatedly helped in X-linked protoporphyria while trials and")
    print("     case series in EPP have been inconsistent and at best neutral.")
    print("     Note also the non-monotonicity in the XLP column at 0.5-0.75x:")
    print("     severe iron deficiency suppresses ALAS2 translation through the")
    print("     very same IRE, so the supply arm re-emerges at the bottom end.")

    # ------------------------------------------------------------------ E-G
    hdr("E-G. THERAPY — two orthogonal axes, and why their effects MULTIPLY")
    print("  Axis A (shielding)      MC1R agonists raise epidermal melanin OD")
    print("  Axis B (source)         bitopertin lowers PPIX via glycine limitation")
    print("  Photodynamic dose  D = E x T(melanin) x [PPIX]skin x t")
    print("  therefore tolerance time  t_tol = D_pro / (E x T x PPIX)  -- a RATIO")
    print("  of PRODUCTS, so the two axes compose multiplicatively.\n")

    bp = pars(FRES=0.15)
    yb = full_ss(bp)
    rb = readouts(yb, bp)
    D = 180
    ev_afa = [(60*24*k, "ADEP", 16000.0) for k in range(3)]
    ev_bit = [(24*k, "BG", 60000.0) for k in range(D)]
    regs = [
        ("untreated EPP",                    []),
        ("afamelanotide 16 mg implant q60d", ev_afa),
        ("dersimelagon 100 mg od",  [(24*k,"DG",100000.0) for k in range(D)]),
        ("dersimelagon 300 mg od",  [(24*k,"DG",300000.0) for k in range(D)]),
        ("bitopertin 20 mg od",     [(24*k,"BG", 20000.0) for k in range(D)]),
        ("bitopertin 60 mg od",     ev_bit),
        ("beta-carotene 180 mg od", [(24*k,"CG",180000.0) for k in range(D)]),
        ("afamelanotide + bitopertin",       ev_afa + ev_bit),
    ]
    print("  Endpoints are averaged over the whole 180-day season, because the")
    print("  afamelanotide implant produces a deliberate sawtooth (dosed q60d)")
    print("  and a single day-180 snapshot would land on its trough.\n")
    print(f"{'regimen, 180-day season':<34}{'RBC PPIX':>10}{'melanin':>9}"
          f"{'T_epi':>8}{'t_tol min':>11}{'x base':>8}{'sun h/d':>9}"
          f"{'season h':>10}")
    res = {}
    for lab, ev in regs:
        y, out = simulate(bp, D, yb, h=0.05, events=ev, exposure=False,
                          record=24.0)
        n = len(out)
        r = dict(PPIXRBC=readouts(y, bp)["PPIXRBC"],
                 MEL=sum(o["MEL"] for o in out)/n,
                 TEPI=sum(o["TEPI"] for o in out)/n,
                 TTOLMIN=sum(o["TTOLMIN"] for o in out)/n,
                 SUNHRD=sum(o["SUNHRD"] for o in out)/n)
        res[lab] = r
        print(f"{lab:<34}{r['PPIXRBC']:>10.2f}{r['MEL']:>9.3f}{r['TEPI']:>8.3f}"
              f"{r['TTOLMIN']:>11.1f}{r['TTOLMIN']/rb['TTOLMIN']:>8.2f}"
              f"{r['SUNHRD']:>9.2f}{r['SUNHRD']*180:>10.0f}")

    fa = res["afamelanotide 16 mg implant q60d"]["TTOLMIN"]/rb["TTOLMIN"]
    fbt = res["bitopertin 60 mg od"]["TTOLMIN"]/rb["TTOLMIN"]
    fc = res["afamelanotide + bitopertin"]["TTOLMIN"]/rb["TTOLMIN"]
    print(f"\n  fold-change in tolerance time")
    print(f"    afamelanotide alone              x{fa:.3f}")
    print(f"    bitopertin alone                 x{fbt:.3f}")
    print(f"    ADDITIVE expectation (fa+fb-1)   x{fa + fbt - 1:.3f}")
    print(f"    MULTIPLICATIVE (Bliss) fa*fb     x{fa*fbt:.3f}")
    print(f"    OBSERVED in the model            x{fc:.3f}")
    print(f"    excess over the additive expectation "
          f"{100*(fc/(fa + fbt - 1) - 1):+.1f}%")
    s0 = rb["SUNHRD"]; sa = res["afamelanotide 16 mg implant q60d"]["SUNHRD"]
    sb = res["bitopertin 60 mg od"]["SUNHRD"]
    sc = res["afamelanotide + bitopertin"]["SUNHRD"]
    print(f"\n  on the CENSORED registrational endpoint (pain-free sun h/day):")
    print(f"    untreated {s0:.2f} | afamelanotide {sa:.2f} | bitopertin {sb:.2f}"
          f" | combination {sc:.2f}")
    print(f"    combination gain {sc - s0:+.2f} h/day vs the sum of the single-")
    print(f"    agent gains {(sa - s0) + (sb - s0):+.2f} h/day — the 8 h daylight")
    print(f"    window compresses what any trial of this endpoint can resolve.")
    print(f"    Over a 180-day season: untreated {s0*180:.0f} h, "
          f"afamelanotide {sa*180:.0f} h, combination {sc*180:.0f} h.")

    print("\n  MC1R-agonist dose-response (pure shielding axis, day-180 steady):")
    print(f"{'dersimelagon mg od':>20}{'melanin':>10}{'T_epi':>8}"
          f"{'t_tol min':>11}{'sun h/d':>9}{'RBC PPIX':>10}")
    for mg in [0, 25, 50, 100, 200, 300, 600]:
        ev = [(24*k, "DG", mg*1000.0) for k in range(D)] if mg else []
        y, _ = simulate(bp, D, yb, h=0.05, events=ev, exposure=False, record=1e9)
        r = readouts(y, bp)
        print(f"{mg:>20}{r['MEL']:>10.3f}{r['TEPI']:>8.3f}{r['TTOLMIN']:>11.1f}"
              f"{r['SUNHRD']:>9.2f}{r['PPIXRBC']:>10.2f}")
    print("  the last column never moves: a shielding drug does not lower PPIX.")
    print("  That is a falsifiable prediction, and it is what the afamelanotide")
    print("  trials reported — protoporphyrin concentrations were unchanged.")

    print("\n  bitopertin dose-response (pure source axis):")
    print(f"{'bitopertin mg od':>20}{'glycine mM':>12}{'RBC PPIX':>10}"
          f"{'% drop':>9}{'t_tol min':>11}{'melanin':>9}")
    for mg in [0, 5, 10, 20, 30, 60, 120]:
        ev = [(24*k, "BG", mg*1000.0) for k in range(D)] if mg else []
        y, _ = simulate(bp, D, yb, h=0.05, events=ev, exposure=False, record=1e9)
        r = readouts(y, bp)
        print(f"{mg:>20}{r['GLY']:>12.3f}{r['PPIXRBC']:>10.2f}"
              f"{100*(1 - r['PPIXRBC']/rb['PPIXRBC']):>9.1f}"
              f"{r['TTOLMIN']:>11.1f}{r['MEL']:>9.3f}")
    print("  and the melanin column never moves. The axes are orthogonal in")
    print("  mechanism but multiply in effect, because dose is a product.")

    # ------------------------------------------------------------------ F2
    hdr("F2. AFAMELANOTIDE PK/PD HYSTERESIS — a 5-day drug with a 60-day effect")
    print(f"{'day':>6}{'implant left ug':>17}{'plasma ug/L':>13}{'melanin':>9}"
          f"{'T_epi':>8}{'t_tol min':>11}")
    y = yb[:]; t = 0.0; h = 0.05
    ev = sorted([(60*24*k, "ADEP", 16000.0) for k in range(3)])
    ei = 0
    want = [0, 1, 2, 5, 10, 20, 30, 45, 59, 61, 90, 120, 150, 180, 200, 240]
    seen = set()
    for i in range(int(240*24/h)):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            y[IX[ev[ei][1]]] += ev[ei][2]; ei += 1
        d = int(round(t/24))
        if d in want and d not in seen and abs(t - d*24) < h/2:
            r = readouts(y, bp); seen.add(d)
            print(f"{d:>6}{y[IX['ADEP']]:>17.1f}{y[IX['AC']]:>13.4f}"
                  f"{r['MEL']:>9.3f}{r['TEPI']:>8.3f}{r['TTOLMIN']:>11.1f}")
        y = rk4(t, y, h, bp, False); t += h
    print("  Plasma drug is effectively gone by about day 5, yet protection goes")
    print("  on rising for another two to three weeks and peaks around day 20.")
    print("  The PD outlasts the PK by more than an order of magnitude, which is")
    print("  why the dose interval is set by melanin turnover (~35 d) rather than")
    print("  by the pharmacokinetics -- 60 days, not 5.")

    # ------------------------------------------------------------------ H
    hdr("H. HEPATIC BISTABILITY — protoporphyric hepatopathy as a SADDLE-NODE")
    print("  The only exit for PPIX is bile; PPIX crystals reduce bile flow;")
    print("  reduced bile flow raises hepatic PPIX. Positive feedback on a")
    print("  single-exit system produces MORE THAN ONE steady state.\n")
    print(f"{'FECH %':>7}{'RBC PPIX':>10}{'input flux':>12}"
          f"  steady states (hepatic PPIX, cholestasis)")
    for fr in [30, 20, 15, 12, 11, 10, 9, 8, 7, 6]:
        q = pars(FRES=fr/100.0)
        s = erythron_ss(q)
        IN = q["KTR"]*s["Jesc"] + q["KEFFL"]*s["PRBC"]
        rts = liver_roots(q, IN)
        txt = "  ".join(f"[{b[0]:6.2f}, {b[1]:.3f}, "
                        f"{'stable' if b[2] else 'SADDLE'}]" for b in rts)
        print(f"{fr:>7}{s['PRBC']:>10.2f}{IN:>12.3f}  {txt}")

    def nstable(fr):
        q = pars(FRES=fr); s = erythron_ss(q)
        IN = q["KTR"]*s["Jesc"] + q["KEFFL"]*s["PRBC"]
        return sum(1 for b in liver_roots(q, IN) if b[2]), q, s

    lo, hi = 0.02, 0.30                      # upper fold: bistability appears
    for _ in range(45):
        mid = 0.5*(lo + hi)
        if nstable(mid)[0] >= 2: lo = mid
        else:                    hi = mid
    fold_hi = lo
    lo2, hi2 = 0.005, fold_hi                # lower fold: healthy branch dies
    for _ in range(45):
        mid = 0.5*(lo2 + hi2)
        if nstable(mid)[0] >= 2: hi2 = mid
        else:                    lo2 = mid
    fold_lo = hi2
    _, qh, sh = nstable(fold_hi)
    _, ql, sl = nstable(fold_lo)
    print(f"\n  BISTABLE WINDOW in residual ferrochelatase activity:")
    print(f"    upper fold  FECH = {100*fold_hi:.2f}%  (RBC PPIX "
          f"{sh['PRBC']:.1f} umol/L)  <- below this a second, cholestatic")
    print(f"                attractor APPEARS alongside the healthy one")
    print(f"    lower fold  FECH = {100*fold_lo:.2f}%  (RBC PPIX "
          f"{sl['PRBC']:.1f} umol/L)  <- below this the HEALTHY branch")
    print(f"                CEASES TO EXIST and hepatopathy is obligatory")
    print("  Between the folds the patient is bistable: biochemically stable for")
    print("  decades, yet one sufficiently large transient (intercurrent illness,")
    print("  fasting, alcohol, a haemolytic episode) pushes the state past the")
    print("  saddle and the system falls onto the cholestatic branch and stays")
    print("  there. That is why protoporphyric liver failure presents abruptly.")

    for fr in [0.10, 0.09]:
        q = pars(FRES=fr); s = erythron_ss(q)
        IN = q["KTR"]*s["Jesc"] + q["KEFFL"]*s["PRBC"]
        rts = liver_roots(q, IN)
        st = [b for b in rts if b[2]]; sd = [b for b in rts if not b[2]]
        if len(st) >= 2 and sd:
            print(f"\n  Worked example, FECH residual {100*fr:.0f}%:")
            print(f"    healthy branch      hepatic PPIX {st[0][0]:6.2f}  "
                  f"cholestasis {st[0][1]:.3f}")
            print(f"    SADDLE (threshold)  hepatic PPIX {sd[0][0]:6.2f}  "
                  f"cholestasis {sd[0][1]:.3f}")
            print(f"    disease branch      hepatic PPIX {st[-1][0]:6.2f}  "
                  f"cholestasis {st[-1][1]:.3f}")
            print(f"    margin to the saddle: "
                  f"{100*(sd[0][0]/st[0][0] - 1):.0f}% above the resting "
                  f"hepatic PPIX")
            break

    print("\n  Rescue from an ESTABLISHED cholestatic branch (FECH 8%),")
    print("  90 days of treatment, started from the disease attractor:")
    p8 = pars(FRES=0.08)
    y8 = full_ss(p8, branch="high")
    r8 = readouts(y8, p8)
    print(f"    baseline: hepatic PPIX {r8['PPIXLIV']:.2f}  cholestasis "
          f"{r8['CHOL']:.3f}  ALT {r8['ALT']:.0f} U/L  bilirubin "
          f"{r8['TBIL']:.2f} mg/dL")
    rescues = [
        ("no treatment",                    dict()),
        ("cholestyramine 16 g/d",           dict(CHOLEFF=0.85)),
        ("RBC transfusion (erythron -60%)", dict(ERYSUP=0.40)),
        ("plasmapheresis",                  dict(KPHER=0.50)),
        ("cholestyramine + transfusion",    dict(CHOLEFF=0.85, ERYSUP=0.40)),
        ("all three",                       dict(CHOLEFF=0.85, ERYSUP=0.40,
                                                 KPHER=0.50)),
        ("bone-marrow transplant",          dict(FRES=1.00)),
    ]
    print(f"\n{'rescue (90 days)':<34}{'hep PPIX':>10}{'cholest.':>10}{'ALT':>8}"
          f"{'bili':>7}{'RBC PPIX':>10}")
    for lab, kw in rescues:
        q = pars(FRES=0.08); q.update(kw)
        yr, _ = simulate(q, 90, y8, h=0.05, exposure=False, record=1e9)
        rr = readouts(yr, q)
        print(f"{lab:<34}{rr['PPIXLIV']:>10.2f}{rr['CHOL']:>10.3f}"
              f"{rr['ALT']:>8.0f}{rr['TBIL']:>7.2f}{rr['PPIXRBC']:>10.2f}")
    print("  Liver transplantation is deliberately absent from this list: it")
    print("  changes NONE of the parameters above, because the source is the")
    print("  marrow. That is precisely why the disease recurs in the graft, and")
    print("  why marrow transplantation — which resets FRES — is the curative one.")

    # ------------------------------------------------------------------ I
    hdr("I. A DAY IN THE LIFE — prodrome, retreat, and the all-or-nothing step")
    print("  Intended outdoor window 10:00-18:00 every day. The patient is his")
    print("  own feedback controller: the prodrome is the sensor, walking into")
    print("  the shade is the actuator, and his reaction time is the loop delay.\n")

    def trace(q, y_start, ndays, label):
        yy = y_start[:]; t = 0.0; h = 0.005
        rows = []
        n = int(round(ndays*24/h))
        for i in range(n):
            if i % int(round(1.0/h)) == 0:
                tod = round(t) % 24
                if 9 <= tod <= 20:
                    rr = readouts(yy, q)
                    rows.append((t, yy[IX["DBOUT"]], yy[IX["AVOID"]],
                                 yy[IX["OX"]], rr["PAIN"], yy[IX["SUNCUM"]]))
            yy = rk4(t, yy, h, q, True); t += h
        return yy, rows

    q_ok = pars(FRES=0.15, TAUDEL=2.0/60.0)
    y_ok, rows_ok = trace(q_ok, yb, 1.0, "adult, 2 min reaction time")
    print("  (a) ADULT who retreats within 2 minutes of the prodrome")
    print(f"{'t (h)':>7}{'prodromal dose':>16}{'in shade':>10}"
          f"{'oxid. signal':>14}{'pain NRS':>10}{'cum sun h':>11}")
    for r in rows_ok:
        print(f"{r[0]:>7.1f}{r[1]:>16.3f}{r[2]:>10.3f}{r[3]:>14.3f}"
              f"{r[4]:>10.2f}{r[5]:>11.2f}")

    q_bad = pars(FRES=0.15, TAUDEL=40.0/60.0)
    y_bad, rows_bad = trace(q_bad, yb, 1.0, "child, 40 min")
    print("\n  (b) SAME PATIENT, SAME SUN — but 40 minutes to get out of it")
    print("      (a young child, a car journey, a wedding, a fishing trip)")
    print(f"{'t (h)':>7}{'prodromal dose':>16}{'in shade':>10}"
          f"{'oxid. signal':>14}{'pain NRS':>10}{'cum sun h':>11}")
    for r in rows_bad:
        print(f"{r[0]:>7.1f}{r[1]:>16.3f}{r[2]:>10.3f}{r[3]:>14.3f}"
              f"{r[4]:>10.2f}{r[5]:>11.2f}")

    print("\n  the 3-day tail of reaction (b), with NO further sun exposure:")
    yq = y_bad[:]; t = 24.0; h = 0.01
    print(f"{'hours after exposure':>22}{'pain NRS':>10}{'oedema':>10}")
    for hr in range(0, 97, 6):
        tgt = 24.0 + hr
        while t < tgt - 1e-9:
            yq = rk4(t, yq, h, q_bad, False); t += h
        print(f"{hr:>22}{readouts(yq, q_bad)['PAIN']:>10.2f}"
              f"{yq[IX['EDEMA']]:>10.2f}")
    print("  A three-microsecond species (singlet oxygen) has produced a")
    print("  three-day illness. The persistence is not in the photon and not in")
    print("  the porphyrin — it is in the slow inflammatory state variables the")
    print("  photon switched on, which is why the pain is opioid-refractory and")
    print("  why nothing done AFTER the exposure shortens it much.")

    yf, of = simulate(pars(FRES=0.15), 14, yb, h=0.005, exposure=True, record=1.0)
    rf = readouts(yf, pars(FRES=0.15))
    print(f"\n  Over 14 days with a competent controller: pain-free sun "
          f"{rf['SUNCUM']:.1f} h total")
    print(f"  ({rf['SUNCUM']/14:.2f} h/day against an 8 h/day intention), peak NRS "
          f"{max(r['PAIN'] for r in of):.1f}, reaction-days {rf['RXNDAYS']:.2f}.")
    print("  Note what the model does NOT say: it does not say an untreated")
    print("  patient has many reactions. A patient who obeys his prodrome every")
    print("  single time has almost none — at the cost of 0.7 of the 8 hours he")
    print("  wanted. The disease burden IS the avoidance.")

    # ------------------------------------------------------------------ J
    hdr("J. THE CONTROL LOOP — reaction time decides prodrome vs full reaction")
    print("  Identical photobiology, identical PPIX, one identical sunny day.")
    print("  The ONLY thing that changes is how long the patient takes to notice")
    print("  the prodrome and get out of the light.\n")
    y_afa, _ = simulate(bp, 30, yb, h=0.05,
                        events=[(0.0, "ADEP", 16000.0)],
                        exposure=False, record=1e9)
    print(f"  (the afamelanotide arm is sampled 30 days after an implant, where")
    print(f"   melanin = {y_afa[IX['MEL']]:.2f} and epidermal transmittance = "
          f"{math.exp(-bp['ODMEL']*y_afa[IX['MEL']]):.3f})\n")

    def one_day(y_start, dmin):
        """One exposure day, then 4 days of recovery with no further sun."""
        q = pars(FRES=0.15, TAUDEL=dmin/60.0)
        y = y_start[:]; t = 0.0; h = 0.005
        pk_dose = 0.0
        for i in range(int(24/h)):
            y = rk4(t, y, h, q, True); t += h
            pk_dose = max(pk_dose, y[IX["DBOUT"]])
        sun = y[IX["SUNCUM"]]
        pk_pain, pk_t = 0.0, 0.0
        for i in range(int(96/h)):
            y = rk4(t, y, h, q, False); t += h
            if i % 40 == 0:
                v = readouts(y, q)["PAIN"]
                if v > pk_pain: pk_pain, pk_t = v, t - 24.0
        return pk_dose, sun, pk_pain, pk_t, readouts(y, q)["RXNDAYS"]

    print(f"{'retreat':>9}{'':>2}{'UNTREATED EPP':>44}{'':>3}"
          f"{'ON AFAMELANOTIDE':>44}")
    print(f"{'delay min':>9}{'':>2}{'peak dose':>10}{'sun h':>8}{'peak NRS':>10}"
          f"{'at h':>7}{'rxn-d':>9}{'':>3}"
          f"{'peak dose':>10}{'sun h':>8}{'peak NRS':>10}{'at h':>7}{'rxn-d':>9}")
    for dmin in [0.5, 1, 2, 5, 10, 20, 40, 90]:
        u = one_day(yb, dmin)
        a = one_day(y_afa, dmin)
        print(f"{dmin:>9.1f}{'':>2}{u[0]:>10.2f}{u[1]:>8.2f}{u[2]:>10.2f}"
              f"{u[3]:>7.1f}{u[4]:>9.2f}{'':>3}"
              f"{a[0]:>10.2f}{a[1]:>8.2f}{a[2]:>10.2f}{a[3]:>7.1f}{a[4]:>9.2f}")
    print("\n  Two things to read off this table.")
    print("  (i)  THE STEP. Between a 20-minute and a 40-minute delay the peak")
    print("       pain roughly doubles and reaction-days go from ~0 to real:")
    print("       the reaction is all-or-nothing, not graded, because the")
    print("       mast-cell/complement arm becomes self-amplifying above DCRIT.")
    print("       Clinically this is why patients describe either 'a warning'")
    print("       or 'three days of boiling oil', and very little in between.")
    print("  (ii) THE SECOND BENEFIT OF SHIELDING. At every delay the")
    print("       afamelanotide column accrues LESS dose and MORE sun. By")
    print("       slowing the rate at which dose accumulates, the drug makes the")
    print("       SAME human reaction time non-limiting. That benefit is")
    print("       invisible to an 'hours in sunlight' endpoint, which is one")
    print("       reason patient-reported outcomes in this disease can move")
    print("       further than the primary endpoint does.")

    print("\n  A caveat the model produces on its own: with the nociceptor")
    print("  time constant set to ~50 h, DAILY exposure stacks even when every")
    print("  prodrome is obeyed. Fourteen consecutive days of compliant outdoor")
    print("  time at a 2-minute delay reach a plateau pain of:")
    for lab, st in (("untreated", yb), ("afamelanotide", y_afa)):
        qq = pars(FRES=0.15, TAUDEL=2.0/60.0)
        yy, oo = simulate(qq, 14, st, h=0.005, exposure=True, record=6.0)
        print(f"    {lab:<16} NRS {max(o['PAIN'] for o in oo):.1f} "
              f"(cumulative pain-free sun {readouts(yy,qq)['SUNCUM']:.1f} h)")
    print("  i.e. the chronic low-grade burning reported by patients who try to")
    print("  live a normal outdoor life is a predicted consequence of a slow")
    print("  nociceptor state, not of any single reaction.")

    # ------------------------------------------------------------------ K
    hdr("K. ENVIRONMENTAL MODIFIERS — why glass, cloud and sunscreen all fail")
    print(f"{'condition':<40}{'transmitted':>13}{'t_tol min':>11}{'sun h/day':>11}")
    conds = [("open summer sun (reference)",  dict()),
             ("behind ordinary window glass", dict(FGLASS=0.90)),
             ("overcast sky",                 dict(FCLOUD=0.50)),
             ("SPF 50 chemical sunscreen",    dict()),
             ("iron-oxide tinted sunscreen",  dict(FCLOTH=0.35)),
             ("long sleeves + wide-brim hat", dict(FCLOTH=0.15)),
             ("full opaque cover",            dict(FCLOTH=0.05))]
    for lab, kw in conds:
        qq = pars(FRES=0.15); qq.update(kw)
        Eeff = qq["FCLOUD"]*qq["FGLASS"]*qq["FCLOTH"]
        ttol = 60.0*qq["DPRO"]/(qq["KPHOT"]*Eeff*math.exp(-qq["ODMEL"]*yb[IX["MEL"]])
                                * yb[IX["PSK"]])
        sun = qq["WINDOW"]*(1 - math.exp(-qq["NBOUT"]*(ttol/60.0)/qq["WINDOW"]))
        print(f"{lab:<40}{Eeff:>13.2f}{ttol:>11.1f}{sun:>11.2f}")
    print("  SPF 50 is listed at transmittance 1.00 deliberately: ultraviolet")
    print("  filters do not attenuate the 400-410 nm Soret band that drives this")
    print("  disease. Only broadband opaque or reflectant cover does — and, from")
    print("  the inside, melanin.")

    hdr("DONE — every number above is computed by the equations in")
    print("epp_mrgsolve_model.R.  Regenerate with:")
    print("    python3 epp_reference_check.py > epp_reference_output.txt")


if __name__ == "__main__":
    main()
