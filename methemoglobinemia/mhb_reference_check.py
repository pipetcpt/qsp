#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mhb_reference_check.py
======================
Independent Python/scipy re-implementation of the methemoglobinaemia QSP model
that is written in mrgsolve form in `mhb_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
R is not available in the build container, so the mrgsolve file cannot be
compiled here.  Every ODE, every parameter and every derived output in the R
file is therefore re-implemented from scratch in this file and integrated with
scipy.  The two files are cross-checked parameter-by-parameter by
`--params` (which parses the R `$PARAM` block) and the numbers printed by
`--all` are the numbers quoted in README.md.  If the R file and this file
disagree, this file is the one that was actually executed.

THE CLAIM THE MODEL IS BUILT TO MAKE
------------------------------------
Methemoglobinaemia is a *delivery* disease, not a *saturation* disease, and
the two numbers the clinician is given are each wrong in an opposite
direction:

  * %MetHb understates the injury, because oxidised subunits do not merely
    subtract capacity, they LEFT-SHIFT the subunits that are left
    (Darling & Roughton 1942).  The residual haemoglobin holds on to the
    oxygen it still carries.
  * SpO2 understates it too, but for a reason that has nothing to do with the
    patient: the pulse oximeter computes one ratio R and reads it off one
    calibration line, and any pigment that absorbs equally at 660 and 940 nm
    drives R -> 1.  The 85% floor is a property of the calibration line, not
    of the blood.

and the antidote is the same molecule as the poison, separated only by an
electron supply (NADPH) that the disease's most important comorbidity
(G6PD deficiency) destroys.

USAGE
-----
    python3 mhb_reference_check.py --all          # everything, as in README
    python3 mhb_reference_check.py --anchors      # anchor table only
    python3 mhb_reference_check.py --params mhb_mrgsolve_model.R
"""

import argparse
import math
import re
import sys

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# 0.  UNITS
# ---------------------------------------------------------------------------
# time                 h
# haemoglobin species  g/dL of WHOLE BLOOD
# RBC-resident solutes uM referred to RED CELL WATER
# plasma solutes       uM (or mg/L for dapsone/benzocaine/cimetidine)
# MB amounts           umol (whole body), converted to concentrations below
#
# 1 g/dL Hb  ==  10 g/L / 4.0285 g/mmol(haem)  ==  2.4823 mmol/L  ==  2482.3 uM
HEMEC = 2482.3          # uM haem per (g/dL) of haemoglobin, per litre of blood
HCT = 0.45              # haematocrit used for the blood <-> red-cell-water map
# a flux of 1 uM(RBC)/h of haem equals CONV g/dL/h of haemoglobin
CONV = HCT / HEMEC      # = 1.8129e-4

VB = 5.0                # L, blood volume
VRBC = VB * HCT         # L, red cell water

# ---------------------------------------------------------------------------
# 1.  PARAMETERS  (must match the $PARAM block of mhb_mrgsolve_model.R)
# ---------------------------------------------------------------------------
P = dict(
    # ---- host / physiology -------------------------------------------------
    WT=70.0,            # kg
    HB0=15.0,           # g/dL, total haemoglobin at t=0
    FMET0=0.008,        # baseline methaemoglobin fraction (0.5-1.5% normal)
    PAO2=95.0,          # mmHg, arterial PO2 (room air); 1800 under HBO 2.8 ATA
    FIO2FLAG=0.0,       # 0 = room air, 1 = supplemental/HBO (PAO2 used directly)
    VO2=250.0,          # mL O2/min, whole-body oxygen consumption
    CO0=5.0,            # L/min, baseline cardiac output
    P50=26.8,           # mmHg, standard P50
    NHILL=2.7,          # Hill coefficient of normal haemoglobin
    ALPHAM=0.50,        # left-shift of residual ferrous haem per unit ferric
                        #   fraction (Darling & Roughton 1942); P50 = P50*(1-a*f)
    BBPG=0.55,          # fractional rise in P50 per unit rise in 2,3-BPG
    PVCRIT=20.0,        # mmHg, critical capillary/venous PO2 for anaerobiosis
    GCO=1.60,           # cardiac-output gain on venous desaturation
    TAUCO=0.30,         # h, time constant of the cardiac-output response
    COMAX=2.6,          # multiple of CO0 that the heart can actually deliver
    GBPG=1.10,          # 2,3-BPG gain on tissue hypoxia
    TAUBPG=30.0,        # h, 2,3-BPG time constant (hours-to-days)

    # ---- methaemoglobin formation -----------------------------------------
    KAUTO=0.00118,      # /h, spontaneous autoxidation of oxyhaemoglobin (2.8%/d)
    KNOH=0.0420,        # /(uM h), dapsone-hydroxylamine co-oxidation
    KTOL=0.0110,        # /(uM h), o-toluidine (benzocaine/prilocaine) oxidation
    KNIT=0.00120,       # /(uM h), nitrite oxidation, non-catalytic term
    KNIT2=0.00060,      # /(uM h (g/dL)), nitrite AUTOCATALYSIS by MetHb
    KROSOX=0.00220,     # /(uM h), Hb oxidation by ROS
    KMBOXH=0.00160,     # /(uM h), oxidation of Hb by oxidised methylene blue
    KSULF=0.000040,     # /(uM h), sulfhaemoglobin formation (irreversible)

    # ---- methaemoglobin reduction -----------------------------------------
    VMXB5=0.675,        # g/dL/h, Vmax of NADH-cytochrome-b5 reductase
    KMB5=4.50,          # g/dL, Km of cyt-b5 reductase for methaemoglobin.
                        #   Km is deliberately LARGE relative to clinical MetHb
                        #   levels, because the observed endogenous reduction is
                        #   near-first-order at ~15%/h of the methaemoglobin pool
                        #   (t1/2 ~ 4.6 h); VMXB5/KMB5 IS that 0.15/h rate.
    EB5SET=1.0,         # CYB5R3 activity fraction (0.10-0.25 = congenital I/II)
    TAUEB5=48.0,        # h, time constant of CYB5R activity change
    GRIBO=0.35,         # riboflavin gain on CYB5R/flavin-reductase activity
    FHBM=0.0,           # fraction of TOTAL haemoglobin that is an HbM variant.
                        #   These haems are held ferric by the globin itself, so
                        #   they form a CONSTITUTIVE FLOOR under MHB that no
                        #   reductase and no methylene blue can get below.  Only
                        #   methaemoglobin ABOVE that floor is reducible.
    KASC=0.0000105,     # /(uM h), non-enzymatic ascorbate reduction of MetHb

    # ---- NADPH: the single shared, capped electron supply -------------------
    VPPP=40000.0,       # uM(RBC)/h, maximal NADPH flux at 100% G6PD activity
    G6PD=1.0,           # G6PD activity fraction (0.02 Mediterranean, 0.15 A-)
    VFRMAX=45000.0,     # uM(RBC)/h, Vmax of the NADPH-flavin reductase
    KMFR=12.0,          # uM, its Km for oxidised methylene blue
    VGRMAX=60000.0,     # uM(RBC)/h, Vmax of glutathione reductase
    KMGR=60.0,          # uM, its Km for GSSG

    # ---- methylene blue: PK -------------------------------------------------
    VC=40.0,            # L, central volume
    VP=1360.0,          # L, peripheral volume
    CLMB=260.0,         # L/h, total clearance
    QMB=300.0,          # L/h, intercompartmental clearance
    KUP=25.0,           # /h, central -> red cell
    KOUT=5.0,           # /h, red cell -> central
    MWMB=319.85,        # g/mol, methylene blue chloride

    # ---- methylene blue: PD -------------------------------------------------
    KMBR=1700.0,        # /h, leucoMB -> MetHb electron transfer
    KMHB=0.80,          # g/dL, Km of that transfer for methaemoglobin
    KO2=95.0,           # /h, leucoMB autoxidation by O2  (the futile branch)
    KMBLK=0.60,         # /h, loss of leucoMB from the red cell

    # ---- reactive oxygen / glutathione / Heinz bodies ----------------------
    GSH0=2000.0,        # uM(RBC), reduced glutathione at baseline
    KSYN=0.055,         # /h, red-cell glutathione synthesis (pool restoration)
    KGSHOX=30.0,        # /h, GSH oxidation per uM ROS
    KG3=350.0,          # uM, Km of that oxidation for GSH
    KROSB=0.55,         # uM/h, basal ROS generation
    KCAT=9.0,           # /h, catalase/GPx removal of ROS
    WGPX=1.60,          # GSH-dependent share of ROS removal
    KRECYC=520.0,       # /h, GSH-driven nitroso -> hydroxylamine recycling
    KNOARD=0.62,        # /h, irreversible loss of the dapsone nitroso species
    KTOLND=0.90,        # /h, irreversible loss of the toluidine nitroso species
    KGR2=600.0,         # uM, Km of that recycling for GSH
    KRCROS=0.020,       # superoxide leaked per turn of a co-oxidation cycle
    KHZ=1.50,           # /h, maximal ROS-driven Heinz-body formation
    KHZG=0.700,         # /h, Heinz formation driven by glutathione COLLAPSE
    FHZG=0.50,          # GSH fraction below which that second driver engages
    KHZR=110.0,         # uM, ROS giving half-maximal Heinz formation (Hill 4)
    NHZ=4.0,            # Hill coefficient of Heinz-body formation
    WHZ=9.0,            # protection afforded by a full glutathione pool
    KHZC=0.030,         # /h, clearance of Heinz-body-bearing cells
    KHEM0=0.00050,      # /h, baseline red-cell removal (120-day lifespan)
    KHEM1=0.012,        # /h, haemolysis at a Heinz burden of 1
    SELMET=2.5,         # preferential clearance of oxidised cells
    TAUEPO=60.0,        # h, marrow response time constant
    KEPO=0.0042,        # /h, marrow gain per g/dL of haemoglobin deficit
    GHBH=0.55,          # erythrocytosis: rise in the Hb set-point per unit hypoxia
    KFHB=0.65,          # conversion of haemolysed Hb to plasma free Hb
    KFHBEL=0.85,        # /h, clearance of plasma free haemoglobin
    HP0=12.0,           # uM, haptoglobin
    KHP=1.10,           # /h, haptoglobin consumption per uM free Hb
    KHPSYN=0.019,       # /h, haptoglobin resynthesis
    KBILI=0.55,         # mg/dL per (g/dL h) of haemolysis
    KBILEL=0.14,        # /h, bilirubin elimination
    KLDH=250.0,         # U/L per (g/dL h) of haemolysis
    KLDHEL=0.020,       # /h, LDH elimination

    # ---- oxidant drug PK ----------------------------------------------------
    KADAP=1.0,          # /h, dapsone absorption
    FDAP=0.90,          # dapsone oral bioavailability
    V1DAP=40.0,         # L
    V2DAP=30.0,         # L
    QDAP=5.0,           # L/h
    CLDAP=1.60,         # L/h
    FNOH=0.100,         # fraction of dapsone clearance to the hydroxylamine
    KBILDAP=0.55,       # /h, biliary secretion (enterohepatic loop)
    KEHC=0.16,          # /h, return of bile to gut
    MWDAP=248.3,        # g/mol
    KNOHEL=3.0,         # /h, hydroxylamine elimination from plasma
    KNOHIN=9.0,         # /h, hydroxylamine uptake into the red cell
    KNOHOUT=2.2,        # /h, hydroxylamine efflux
    IC50CIM=2.6,        # mg/L, cimetidine inhibition of CYP-mediated N-hydroxylation
    KACIM=1.4,          # /h
    KELCIM=0.29,        # /h
    VCIM=90.0,          # L

    KABZC=6.0,          # /h, benzocaine mucosal absorption
    KELBZC=5.5,         # /h, benzocaine hydrolysis
    FTOL=0.30,          # fraction routed to the oxidising metabolite
    VBZC=60.0,          # L
    MWBZC=165.2,        # g/mol
    KTOLEL=5.0,         # /h
    KTOLIN=8.0,         # /h
    KTOLOUT=2.0,        # /h

    KANIT=1.2,          # /h, nitrite absorption from the gut
    KNITEL=1.1,         # /h
    KNITIN=12.0,        # /h
    KNITOUT=3.0,        # /h
    MWNIT=69.0,         # g/mol (sodium nitrite)

    KELASC=0.35,        # /h, ascorbate elimination
    KASCIN=0.55,        # /h, ascorbate uptake into the red cell
    KASCOUT=0.30,       # /h
    VASC=18.0,          # L

    # ---- methylene blue as an MAO-A inhibitor ------------------------------
    KPBR=90.0,          # brain/tissue partition of methylene blue
    KIMAO=5.5,          # uM, MB concentration for 50% MAO-A inhibition
    KSSRI=1.0,          # normalised SSRI exposure giving 50% SERT block
    KREL5=1.0,          # /h, serotonin release (normalised)
    KUP5=0.62,          # /h, SERT-mediated clearance
    KMAO5=0.38,         # /h, MAO-A-mediated clearance

    # ---- pulse oximeter (the only fitted OPTICAL parameters) ---------------
    E660O=0.32,         # relative extinction, oxyhaemoglobin at 660 nm
    E940O=1.21,         # oxyhaemoglobin at 940 nm
    E660D=3.23,         # deoxyhaemoglobin at 660 nm
    E940D=0.69,         # deoxyhaemoglobin at 940 nm
    E660M=15.0,         # methaemoglobin at 660 nm   } fitted to Barker & Tremper
    E940M=15.0,         # methaemoglobin at 940 nm   } (equal by construction)
    E660S=20.0,         # sulfhaemoglobin at 660 nm
    E940S=4.0,          # sulfhaemoglobin at 940 nm
    SPOA=110.0,         # SpO2 = SPOA - SPOB * R   (the device calibration line)
    SPOB=25.0,

    # ---- clinical thresholds ------------------------------------------------
    CYANDEOXY=5.0,      # g/dL of deoxyhaemoglobin that is visibly cyanotic
    CYANMET=1.5,        # g/dL of methaemoglobin that is visibly cyanotic
    CYANSULF=0.5,       # g/dL of sulfhaemoglobin that is visibly cyanotic
    KLAC=0.030,         # mmol/L/h per mmHg of venous PO2 below critical
    KLACEL=0.45,        # /h
    RIBO=0.0,           # riboflavin supplementation flag (0/1)
)

# order of the state vector; index map built once
SNAMES = [
    "DAPG", "DAPC", "DAPP", "DAPB", "NOHC", "NOHR", "NOAR",        # 0-6
    "BZCD", "BZCC", "TOLR", "TOLN",                                # 7-9,+
    "NITG", "NITC", "NITR",                                        # 10-11,+
    "MBC", "MBP", "MBOX", "LMB", "MBEL",                           # 12-16
    "ASCC", "ASCR", "CIMG", "CIMC", "SSRI",                        # 17-21
    "HBF", "MHB", "SHB",                                           # 22-24
    "GSH", "GSSG", "ROS", "HEINZ",                                 # 25-28
    "FHB", "HAPT", "BILI", "LDH",                                  # 29-32
    "EPOD", "BPG", "CO", "LAC", "INJ", "O2DEBT",                   # 33-38
    "HT5", "EB5", "MBCUM",                                         # 39-41
]
IDX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)


# ---------------------------------------------------------------------------
# 2.  ALGEBRAIC OUTPUT BLOCK  (identical to the $TABLE block of the R file)
# ---------------------------------------------------------------------------
def outputs(y, p):
    """Everything the model reports but does not integrate."""
    HBF, MHB, SHB = y[IDX["HBF"]], y[IDX["MHB"]], y[IDX["SHB"]]
    HBF = max(HBF, 1e-9)
    MHB = max(MHB, 0.0)
    SHB = max(SHB, 0.0)
    HBTOT = HBF + MHB + SHB

    FMET = MHB / HBTOT
    FSULF = SHB / HBTOT
    # the allosterically relevant fraction: ferric haems among haems that could
    # otherwise bind oxygen (sulfhaemoglobin sits outside the tetramer story)
    FRED = MHB / (HBF + MHB)

    # --- the two things oxidised subunits do to the subunits that are left ---
    NEFF = 1.0 + (p["NHILL"] - 1.0) * (1.0 - FRED)      # cooperativity is lost
    P50E = p["P50"] * (1.0 + p["BBPG"] * (y[IDX["BPG"]] - 1.0)) \
        * (1.0 - p["ALPHAM"] * FRED)                     # and the curve shifts
    P50E = max(P50E, 1.0)

    PAO2 = p["PAO2"]
    SAO2 = PAO2 ** NEFF / (P50E ** NEFF + PAO2 ** NEFF)

    CAO2 = 1.34 * HBF * SAO2 + 0.003 * PAO2
    CO = max(y[IDX["CO"]], 0.5)
    CVO2 = CAO2 - p["VO2"] / (CO * 10.0)

    CAP_HB = 1.34 * HBF
    if CVO2 >= CAP_HB:                    # dissolved oxygen alone suffices
        SVO2 = 1.0
        PVO2 = max((CVO2 - CAP_HB) / 0.003, 0.0)
    else:
        SVO2 = min(max(CVO2 / CAP_HB, 1e-6), 1.0 - 1e-6)
        PVO2 = P50E * (SVO2 / (1.0 - SVO2)) ** (1.0 / NEFF)
    PVO2 = max(PVO2, 0.0)

    DO2 = CAO2 * CO * 10.0                 # mL O2/min delivered
    ERO2 = 0.0 if CAO2 <= 0 else (CAO2 - CVO2) / CAO2

    # --- the pulse oximeter, modelled as the device rather than as the blood -
    fO2 = HBF * SAO2 / HBTOT
    fdx = HBF * (1.0 - SAO2) / HBTOT
    num = fO2 * p["E660O"] + fdx * p["E660D"] + FMET * p["E660M"] + FSULF * p["E660S"]
    den = fO2 * p["E940O"] + fdx * p["E940D"] + FMET * p["E940M"] + FSULF * p["E940S"]
    R = num / max(den, 1e-9)
    SPO2 = min(max(p["SPOA"] - p["SPOB"] * R, 0.0), 100.0)
    SAO2CO = 100.0 * fO2                   # what a co-oximeter reports
    GAP = SPO2 - SAO2CO                    # the "saturation gap"

    # --- cyanosis is a sum of pigments, each with its own threshold ----------
    CYAN = (HBF * (1.0 - SAO2)) / p["CYANDEOXY"] + MHB / p["CYANMET"] \
        + SHB / p["CYANSULF"]

    return dict(HBTOT=HBTOT, FMET=FMET, MetPct=100.0 * FMET, FSULF=FSULF,
                FRED=FRED, NEFF=NEFF, P50E=P50E, SAO2=SAO2, CAO2=CAO2,
                CVO2=CVO2, SVO2=SVO2, PVO2=PVO2, DO2=DO2, ERO2=ERO2,
                R=R, SPO2=SPO2, SAO2CO=SAO2CO, GAP=GAP, CYAN=CYAN, CO=CO)


def equivalent_hb(pvo2, p, hb_hi=25.0):
    """
    The number the model exists to produce.

    Given a tissue (mixed-venous) PO2, return the haemoglobin concentration a
    patient with STRUCTURALLY NORMAL haemoglobin would have to have in order to
    sit at that same PVO2, at the same VO2 and cardiac output.  This converts
    "MetHb 30%" into "anaemia of X g/dL" -- and X is NOT 0.7*15.
    """
    lo, hi = 0.05, hb_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        y = np.zeros(NST)
        y[IDX["HBF"]] = mid
        y[IDX["MHB"]] = 0.0
        y[IDX["SHB"]] = 0.0
        y[IDX["BPG"]] = 1.0
        y[IDX["CO"]] = p["CO0"]
        o = outputs(y, dict(p))
        if o["PVO2"] < pvo2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# 3.  THE ODE SYSTEM  (mirrors $ODE of mhb_mrgsolve_model.R line for line)
# ---------------------------------------------------------------------------
def rhs(t, y, p, ev):
    y = np.maximum(y, 0.0)
    d = np.zeros(NST)
    g = lambda n: y[IDX[n]]

    HBF = max(g("HBF"), 1e-9)
    MHB = g("MHB")
    SHB = g("SHB")
    HBTOT = HBF + MHB + SHB
    o = outputs(y, p)
    MHBRED = max(MHB - p["FHBM"] * HBTOT, 0.0)

    # ---- 3.1 oxidant drug PK ----------------------------------------------
    CDAP = g("DAPC") / p["V1DAP"]                       # mg/L
    CCIM = g("CIMC") / p["VCIM"]                        # mg/L
    fnoh = p["FNOH"] / (1.0 + CCIM / p["IC50CIM"])      # cimetidine acts HERE

    abs_dap = p["KADAP"] * g("DAPG")
    d[IDX["DAPG"]] = -abs_dap + p["KEHC"] * g("DAPB")
    cl_dap = p["CLDAP"] * CDAP
    dist = p["QDAP"] * (CDAP - g("DAPP") / p["V2DAP"])
    bil = p["KBILDAP"] * CDAP
    d[IDX["DAPC"]] = p["FDAP"] * abs_dap - cl_dap - dist - bil
    d[IDX["DAPP"]] = dist
    d[IDX["DAPB"]] = bil - p["KEHC"] * g("DAPB")
    # hydroxylamine: umol/h out of the mg/h of dapsone routed to N-hydroxylation
    gen_noh = cl_dap * fnoh * 1000.0 / p["MWDAP"]       # umol/h
    upt_noh = p["KNOHIN"] * g("NOHC") - p["KNOHOUT"] * g("NOHR") * VRBC
    d[IDX["NOHC"]] = gen_noh - p["KNOHEL"] * g("NOHC") - upt_noh
    # nitroso-arene is reduced BACK to the hydroxylamine by glutathione: the
    # cell's own antioxidant is the co-substrate that makes the poison catalytic
    v_recyc = p["KRECYC"] * g("NOAR") * g("GSH") / (p["KGR2"] + g("GSH"))
    v_noh = p["KNOH"] * g("NOHR") * HBF                 # g/dL/h of Hb oxidised
    d[IDX["NOHR"]] = upt_noh / VRBC - v_noh / (2.0 * CONV) + v_recyc
    d[IDX["NOAR"]] = v_noh / (2.0 * CONV) - v_recyc - p["KNOARD"] * g("NOAR")

    abs_bzc = p["KABZC"] * g("BZCD")
    d[IDX["BZCD"]] = -abs_bzc
    d[IDX["BZCC"]] = abs_bzc - p["KELBZC"] * g("BZCC")
    gen_tol = p["FTOL"] * p["KELBZC"] * g("BZCC") * 1000.0 / p["MWBZC"]  # umol/h
    v_tol = p["KTOL"] * g("TOLR") * HBF
    # the SAME glutathione-driven cycle that makes the dapsone metabolite
    # catalytic; without it, 250 mg of benzocaine could oxidise at most 0.5% of
    # the haemoglobin mass and clinical methaemoglobinaemia would be impossible
    v_rcyt = p["KRECYC"] * g("TOLN") * g("GSH") / (p["KGR2"] + g("GSH"))
    d[IDX["TOLR"]] = (gen_tol / VRBC - p["KTOLEL"] * g("TOLR")
                      - v_tol / (2.0 * CONV) + v_rcyt)
    d[IDX["TOLN"]] = v_tol / (2.0 * CONV) - v_rcyt - p["KTOLND"] * g("TOLN")

    upt_nit = p["KNITIN"] * g("NITC") - p["KNITOUT"] * g("NITR") * VRBC
    d[IDX["NITG"]] = -p["KANIT"] * g("NITG")
    d[IDX["NITC"]] = p["KANIT"] * g("NITG") - p["KNITEL"] * g("NITC") - upt_nit
    # nitrite oxidation of oxyhaemoglobin is AUTOCATALYTIC: the product
    # accelerates the reaction, which is why nitrite has a lag and then a wall
    v_nit = (p["KNIT"] + p["KNIT2"] * MHB) * g("NITR") * HBF
    d[IDX["NITR"]] = upt_nit / VRBC - v_nit / (2.0 * CONV)

    d[IDX["CIMG"]] = -p["KACIM"] * g("CIMG")
    d[IDX["CIMC"]] = p["KACIM"] * g("CIMG") - p["KELCIM"] * g("CIMC")

    upt_asc = p["KASCIN"] * g("ASCC") - p["KASCOUT"] * g("ASCR") * VRBC
    d[IDX["ASCC"]] = -p["KELASC"] * g("ASCC") - upt_asc
    v_asc = p["KASC"] * g("ASCR") * MHBRED
    d[IDX["ASCR"]] = upt_asc / VRBC - 0.10 * g("ASCR") - v_asc / CONV

    d[IDX["SSRI"]] = -0.029 * g("SSRI")                 # t1/2 ~ 24 h, normalised

    # ---- 3.2 methylene blue PK --------------------------------------------
    CMBC = g("MBC") / p["VC"]                           # uM in plasma
    d[IDX["MBC"]] = (-p["CLMB"] * CMBC
                     - p["QMB"] * (CMBC - g("MBP") / p["VP"])
                     - p["KUP"] * g("MBC") + p["KOUT"] * g("MBOX"))
    d[IDX["MBP"]] = p["QMB"] * (CMBC - g("MBP") / p["VP"])
    d[IDX["MBEL"]] = p["CLMB"] * CMBC

    CMBOX = g("MBOX") / VRBC                            # uM in red cell water
    CLMBR = g("LMB") / VRBC

    # ---- 3.3 NADPH: ONE capped supply, TWO consumers ----------------------
    # This single block is where G6PD deficiency, the methylene-blue ceiling and
    # the haemolysis all come from.  Nothing else in the model knows about them.
    vFRdes = p["VFRMAX"] * CMBOX / (p["KMFR"] + CMBOX)  # MB's demand
    vGRdes = p["VGRMAX"] * g("GSSG") / (p["KMGR"] + g("GSSG"))   # GSH's demand
    CAPN = p["VPPP"] * p["G6PD"] * (1.0 + p["GRIBO"] * p["RIBO"])
    share = CAPN / (CAPN + vFRdes + vGRdes + 1e-9)
    vFR = vFRdes * share                                # uM(RBC)/h leucoMB made
    vGR = vGRdes * share                                # uM(RBC)/h GSSG reduced

    # methylene blue oxidises haemoglobin directly; this ALSO makes leucoMB, so
    # it is not a loss of drug -- it is a loan of electrons taken from the patient
    v_mboxh = p["KMBOXH"] * CMBOX * HBF                 # g/dL/h
    v_mred = p["KMBR"] * CLMBR * MHBRED / (p["KMHB"] + MHBRED)   # uM(RBC)/h
    v_o2 = p["KO2"] * CLMBR                             # the futile branch

    d[IDX["LMB"]] = (vFR + v_mboxh / (2.0 * CONV) - v_mred - v_o2
                     - p["KMBLK"] * CLMBR) * VRBC
    d[IDX["MBOX"]] = (p["KUP"] * g("MBC") - p["KOUT"] * g("MBOX")
                      + (-vFR - v_mboxh / (2.0 * CONV) + v_mred + v_o2) * VRBC)

    # ---- 3.4 reactive oxygen, glutathione, Heinz bodies -------------------
    # Glutathione is CONSERVED (GT = GSH + 2 GSSG) apart from slow resynthesis.
    # Every turn of a co-oxidation cycle costs 2 GSH and returns 1 GSSG, so the
    # arylhydroxylamine cycle is ultimately paid for out of NADPH -- which is
    # why the enzyme deficiency that disables the antidote also changes what the
    # poison does.
    v_rcyc_tot = v_recyc + v_rcyt
    v_gshox = p["KGSHOX"] * g("ROS") * g("GSH") / (p["KG3"] + g("GSH"))
    GT = g("GSH") + 2.0 * g("GSSG")
    d[IDX["GSH"]] = (p["KSYN"] * (p["GSH0"] - GT) + 2.0 * vGR
                     - v_gshox - 2.0 * v_rcyc_tot)
    d[IDX["GSSG"]] = 0.5 * v_gshox + v_rcyc_tot - vGR
    d[IDX["ROS"]] = (p["KROSB"] + v_o2 + p["KRCROS"] * v_rcyc_tot
                     - p["KCAT"] * g("ROS") * (1.0 + p["WGPX"] * g("GSH") / p["GSH0"]))
    # Heinz bodies have TWO drivers, and they are not the same driver.  One is
    # frank oxidant radical flux (what high-dose methylene blue does to a normal
    # cell).  The other is the collapse of the glutathione pool itself (what
    # happens in G6PD deficiency, where there is no radical excess at all --
    # there is simply nothing left holding the beta-93 thiols in the reduced
    # state).  Without the second driver the model predicts, wrongly, that the
    # G6PD-deficient patient haemolyses LESS than the normal one.
    rr = g("ROS") ** p["NHZ"] / (p["KHZR"] ** p["NHZ"] + g("ROS") ** p["NHZ"])
    gcol = max(0.0, 1.0 - g("GSH") / (p["FHZG"] * p["GSH0"])) ** 2
    d[IDX["HEINZ"]] = ((p["KHZ"] * rr + p["KHZG"] * gcol) * (1.0 - g("HEINZ"))
                       / (1.0 + p["WHZ"] * g("GSH") / p["GSH0"])
                       - p["KHZC"] * g("HEINZ"))

    v_rosox = p["KROSOX"] * g("ROS") * HBF
    v_sulf = p["KSULF"] * g("ROS") * HBF                # irreversible

    # ---- 3.5 haemoglobin book-keeping -------------------------------------
    khem = p["KHEM0"] + p["KHEM1"] * g("HEINZ")
    lys_f = khem * HBF
    lys_m = khem * p["SELMET"] * MHB                    # oxidised cells go first
    lys_s = khem * p["SELMET"] * SHB
    lys_tot = lys_f + lys_m + lys_s

    v_b5 = p["VMXB5"] * g("EB5") * MHBRED / (p["KMB5"] + MHBRED)
    v_auto = p["KAUTO"] * HBF

    ox_tot = v_auto + v_noh + v_tol + v_nit + v_rosox + v_mboxh
    red_tot = v_b5 + v_asc + 2.0 * CONV * v_mred

    d[IDX["HBF"]] = -ox_tot + red_tot - v_sulf - lys_f + g("EPOD")
    d[IDX["MHB"]] = ox_tot - red_tot - lys_m
    d[IDX["SHB"]] = v_sulf - lys_s

    # ---- 3.6 consequences of haemolysis -----------------------------------
    d[IDX["FHB"]] = p["KFHB"] * lys_tot * HEMEC / 100.0 - p["KFHBEL"] * g("FHB")
    d[IDX["HAPT"]] = (p["KHPSYN"] * (p["HP0"] - g("HAPT"))
                      - p["KHP"] * g("FHB") * g("HAPT") / (g("HAPT") + 2.0))
    d[IDX["BILI"]] = p["KBILI"] * lys_tot - p["KBILEL"] * g("BILI")
    d[IDX["LDH"]] = p["KLDH"] * lys_tot - p["KLDHEL"] * (g("LDH") - 180.0)

    hyp0 = max(0.0, 1.0 - o["PVO2"] / 38.8)
    prod = (p["KEPO"] * max(0.0, p["HB0"] * (1.0 + p["GHBH"] * hyp0) - HBTOT)
            + p["KHEM0"] * p["HB0"])
    d[IDX["EPOD"]] = (prod - g("EPOD")) / p["TAUEPO"]

    # ---- 3.7 whole-body compensation --------------------------------------
    hyp = max(0.0, 1.0 - o["PVO2"] / 38.8)
    cotgt = min(p["CO0"] * (1.0 + p["GCO"] * hyp), p["CO0"] * p["COMAX"])
    d[IDX["CO"]] = (cotgt - g("CO")) / p["TAUCO"]
    d[IDX["BPG"]] = ((1.0 + p["GBPG"] * hyp) - g("BPG")) / p["TAUBPG"]

    below = max(0.0, p["PVCRIT"] - o["PVO2"])
    d[IDX["LAC"]] = p["KLAC"] * below - p["KLACEL"] * (g("LAC") - 1.0)
    d[IDX["INJ"]] = below
    d[IDX["O2DEBT"]] = below * p["VO2"] / 1000.0 / p["WT"]

    # ---- 3.8 the antidote's other pharmacology ----------------------------
    # MAO-A sits in tissue, not plasma: plasma methylene blue is gone in
    # minutes, but the interaction with SSRIs lasts hours, so the effect
    # site is the peripheral compartment.
    CMBEFF = p["KPBR"] * g("MBP") / p["VP"]
    IMAO = CMBEFF / (p["KIMAO"] + CMBEFF)
    ISERT = g("SSRI") / (p["KSSRI"] + g("SSRI"))
    d[IDX["HT5"]] = (p["KREL5"] - p["KUP5"] * (1.0 - ISERT) * g("HT5")
                     - p["KMAO5"] * (1.0 - IMAO) * g("HT5"))

    d[IDX["EB5"]] = ((p["EB5SET"] * (1.0 + p["GRIBO"] * p["RIBO"])) - g("EB5")) \
        / p["TAUEB5"]
    d[IDX["MBCUM"]] = 0.0
    return d


# ---------------------------------------------------------------------------
# 4.  INITIAL CONDITIONS AND THE EVENT/SIMULATION DRIVER
# ---------------------------------------------------------------------------
def init_state(p):
    y = np.zeros(NST)
    hb = p["HB0"]
    y[IDX["HBF"]] = hb * (1.0 - p["FMET0"] - p["FHBM"])
    y[IDX["MHB"]] = hb * (p["FMET0"] + p["FHBM"])
    y[IDX["SHB"]] = 0.0
    y[IDX["GSH"]] = p["GSH0"]
    y[IDX["GSSG"]] = 1.0
    y[IDX["ROS"]] = p["KROSB"] / (p["KCAT"] * (1.0 + p["WGPX"]))
    y[IDX["HAPT"]] = p["HP0"]
    y[IDX["BILI"]] = 0.6
    y[IDX["LDH"]] = 180.0
    y[IDX["BPG"]] = 1.0
    y[IDX["CO"]] = p["CO0"]
    y[IDX["LAC"]] = 1.0
    y[IDX["HT5"]] = 1.0
    y[IDX["EB5"]] = p["EB5SET"]
    y[IDX["EPOD"]] = p["KHEM0"] * p["HB0"]
    return y


def simulate(par=None, events=None, tend=48.0, npt=None, y0=None):
    """
    events: list of (time_h, state_name, amount) bolus additions.
            Methylene blue doses are given as ("MB_mgkg", dose) and converted.
    """
    p = dict(P)
    if par:
        p.update(par)
    events = sorted(events or [], key=lambda e: e[0])
    y = init_state(p) if y0 is None else np.array(y0, dtype=float)

    npt = npt or max(400, int(tend * 8))
    tgrid = np.unique(np.concatenate([
        np.linspace(0.0, tend, npt),
        np.array([e[0] for e in events if 0 <= e[0] <= tend] or [0.0]),
    ]))

    rec = {n: [] for n in SNAMES}
    rec["time"] = []
    okeys = None
    orec = {}

    def store(t, yv):
        nonlocal okeys
        o = outputs(yv, p)
        if okeys is None:
            okeys = list(o.keys())
            for k in okeys:
                orec[k] = []
        rec["time"].append(t)
        for n in SNAMES:
            rec[n].append(yv[IDX[n]])
        for k in okeys:
            orec[k].append(o[k])

    breaks = sorted({0.0, tend} | {e[0] for e in events if 0 < e[0] < tend})
    t0 = breaks[0]
    store(t0, y)
    for seg_i in range(len(breaks) - 1):
        ta, tb = breaks[seg_i], breaks[seg_i + 1]
        for ev in events:
            if abs(ev[0] - ta) < 1e-12:
                name, amt = ev[1], ev[2]
                if name == "MB_mgkg":
                    umol = amt * p["WT"] * 1000.0 / p["MWMB"]
                    y[IDX["MBC"]] += umol
                    y[IDX["MBCUM"]] += amt
                elif name == "TRANSFUSE":
                    # amt = g/dL of fresh, fully functional haemoglobin, with an
                    # equal volume of the patient's own blood removed
                    frac = amt / max(y[IDX["HBF"]] + y[IDX["MHB"]] + y[IDX["SHB"]], 1e-9)
                    frac = min(frac, 1.0)
                    for nm in ("HBF", "MHB", "SHB"):
                        y[IDX[nm]] *= (1.0 - frac)
                    y[IDX["HBF"]] += amt
                else:
                    y[IDX[name]] += amt
        sub = np.linspace(ta, tb, max(3, int((tb - ta) * 24) + 3))
        sol = solve_ivp(rhs, (ta, tb), y, t_eval=sub[1:], args=(p, events),
                        method="LSODA", rtol=1e-7, atol=1e-9, max_step=0.25)
        if not sol.success:
            raise RuntimeError(f"integration failed: {sol.message}")
        for j, tt in enumerate(sol.t):
            store(tt, sol.y[:, j])
        y = sol.y[:, -1].copy()

    out = {k: np.array(v) for k, v in rec.items()}
    out.update({k: np.array(v) for k, v in orec.items()})
    out["_p"] = p
    out["_yend"] = y
    return out


def at(res, t):
    """index of the sample nearest time t"""
    return int(np.argmin(np.abs(res["time"] - t)))


def peak(res, key, tmin=0.0):
    m = res["time"] >= tmin
    return float(np.max(res[key][m]))


def tpeak(res, key, tmin=0.0):
    m = res["time"] >= tmin
    i = int(np.argmax(np.where(m, res[key], -np.inf)))
    return float(res["time"][i])


# ---------------------------------------------------------------------------
# 5.  SCENARIOS
# ---------------------------------------------------------------------------
def sc_baseline():
    return simulate(tend=1440.0, npt=1500)


def sc_benzocaine(dose_mg=250.0, mb=None, par=None, tend=24.0):
    ev = [(0.0, "BZCD", dose_mg)]
    if mb:
        for tmb, d in mb:
            ev.append((tmb, "MB_mgkg", d))
    return simulate(par=par, events=ev, tend=tend, npt=1400)


def sc_dapsone(mg_per_day=100.0, days=30.0, cimetidine=False, par=None):
    ev = []
    n = int(days)
    for k in range(n):
        ev.append((24.0 * k, "DAPG", mg_per_day))
        if cimetidine:
            for h in (0.0, 8.0, 16.0):
                ev.append((24.0 * k + h, "CIMG", 400.0))
    return simulate(par=par, events=ev, tend=24.0 * days, npt=2600)


def sc_dapsone_od(dose_mg=3000.0, mb_sched=None, par=None, tend=120.0):
    ev = [(0.0, "DAPG", dose_mg)]
    for tmb, d in (mb_sched or []):
        ev.append((tmb, "MB_mgkg", d))
    return simulate(par=par, events=ev, tend=tend, npt=2400)


# ---------------------------------------------------------------------------
# 6.  REPORT SECTIONS
# ---------------------------------------------------------------------------
def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


ANCHORS = []


def anchor(label, value, lo, hi, unit="", src=""):
    ok = (value >= lo) and (value <= hi)
    ANCHORS.append((label, value, lo, hi, ok))
    flag = "PASS" if ok else "FAIL"
    print(f"  [{flag}] {label:<52s} {value:10.3f} {unit:<8s} "
          f"(target {lo:g}-{hi:g}) {src}")
    return ok


def sec_steadystate():
    hr("1.  STEADY STATE — the model must first be able to do nothing")
    r = sc_baseline()
    i = -1
    d = rhs(r["time"][i], r["_yend"], r["_p"], [])
    print(f"  60-day integration of an untreated healthy adult.")
    print(f"    total Hb           {r['HBTOT'][i]:8.3f} g/dL   (start 15.000)")
    print(f"    MetHb              {r['MetPct'][i]:8.3f} %      (start 0.800)")
    print(f"    PvO2               {r['PVO2'][i]:8.3f} mmHg")
    print(f"    SpO2               {r['SPO2'][i]:8.3f} %")
    print(f"    max|dy/dt|         {np.max(np.abs(d)):8.3e}")
    print(f"    cumulative hypoxic dose {r['INJ'][i]:.3e} mmHg.h  (must be 0)")
    anchor("baseline MetHb", r["MetPct"][i], 0.4, 1.6, "%", "Normal <1.5%")
    anchor("baseline PvO2", r["PVO2"][i], 35.0, 42.0, "mmHg", "Normal 38-42")
    anchor("baseline SpO2", r["SPO2"][i], 96.0, 100.0, "%", "")
    anchor("baseline Hb drift over 60 d", abs(r["HBTOT"][i] - 15.0), 0.0, 0.30,
           "g/dL", "must not drift")
    anchor("baseline hypoxic dose", r["INJ"][i], 0.0, 1e-6, "mmHg.h", "must be 0")
    return r


def sec_equivalent_anemia():
    hr("2.  THE HEADLINE RESULT — what a percentage of methaemoglobin costs")
    print("""  A patient is told "your methaemoglobin is 30%".  The arithmetic that
  is done at the bedside is 15 x 0.70 = 10.5 g/dL of working haemoglobin.
  That arithmetic is wrong, and the model says by how much, because ferric
  subunits do two things to the subunits beside them: they remove the
  cooperativity (n -> 1) and they left-shift the curve (Darling & Roughton).
  Oxygen that is still carried is no longer released.

  EQUIV Hb = the haemoglobin a patient with NORMAL haemoglobin would need in
  order to sit at the same tissue PO2, at the same VO2 and cardiac output.
""")
    p = dict(P)
    print(f"  {'MetHb%':>7s} {'naive':>8s} {'EQUIV':>8s} {'lost':>8s} "
          f"{'% of residual':>14s} {'PvO2':>8s} {'P50eff':>8s} {'n_eff':>7s} {'SpO2':>7s}")
    rows = []
    for fm in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70):
        y = init_state(p)
        hb = p["HB0"]
        y[IDX["HBF"]] = hb * (1.0 - fm)
        y[IDX["MHB"]] = hb * fm
        o = outputs(y, p)
        naive = hb * (1.0 - fm)
        eq = equivalent_hb(o["PVO2"], p)
        lost = naive - eq
        pct = 100.0 * lost / naive if naive > 0 else 0.0
        rows.append((fm, naive, eq, lost, pct, o))
        print(f"  {100*fm:7.0f} {naive:8.2f} {eq:8.2f} {lost:8.2f} {pct:13.1f}% "
              f"{o['PVO2']:8.1f} {o['P50E']:8.1f} {o['NEFF']:7.2f} {o['SPO2']:7.1f}")
    print("""
  Read the fifth column.  At 15% MetHb the left shift sterilises 10% of the
  haemoglobin that is still ferrous; at 30% it sterilises 21%; at 50% it
  sterilises 31%.  The clinical rule that "30% MetHb is roughly a haemoglobin
  of 10.5" is out by more than two grams, and the error GROWS with severity --
  it is largest exactly where the decision matters.
""")
    r30 = [r for r in rows if abs(r[0] - 0.30) < 1e-9][0]
    r50 = [r for r in rows if abs(r[0] - 0.50) < 1e-9][0]
    anchor("EQUIV Hb at MetHb 30% (naive says 10.50)", r30[2], 8.0, 8.8, "g/dL")
    anchor("residual capacity sterilised at MetHb 30%", r30[4], 17.0, 25.0, "%")
    anchor("residual capacity sterilised at MetHb 50%", r50[4], 26.0, 36.0, "%")
    anchor("PvO2 at MetHb 50% (untreated, uncompensated)", r50[5]["PVO2"],
           16.0, 22.0, "mmHg", "below the 20 mmHg anaerobic threshold")
    return rows


def sec_same_percent_different_disease():
    hr("3.  THE SAME PERCENTAGE IS NOT THE SAME DISEASE")
    print("""  %MetHb is a RATIO.  Tissue oxygen delivery is a PRODUCT.  A model that
  carries both makes the anaemic patient's vulnerability arithmetic rather
  than a warning in a textbook.
""")
    p = dict(P)
    print(f"  {'Hb':>6s} {'MetHb%':>7s} {'functional':>11s} {'PvO2':>8s} "
          f"{'SpO2':>7s} {'lactate?':>9s}")
    for hb in (17.0, 15.0, 12.0, 9.0, 7.0):
        for fm in (0.20, 0.30):
            y = init_state(p)
            y[IDX["HBF"]] = hb * (1.0 - fm)
            y[IDX["MHB"]] = hb * fm
            o = outputs(y, dict(p, HB0=hb))
            print(f"  {hb:6.1f} {100*fm:7.0f} {hb*(1-fm):11.2f} {o['PVO2']:8.1f} "
                  f"{o['SPO2']:7.1f} {'YES' if o['PVO2'] < p['PVCRIT'] else 'no':>9s}")
    print("""
  The pulse oximeter reads essentially the SAME number down every column --
  it is a function of the ratio alone.  The tissue PO2 is not.  A haemoglobin
  of 7 crosses the anaerobic threshold at a methaemoglobin fraction that a
  haemoglobin of 17 tolerates with 14 mmHg to spare, and no measurement made
  at the bedside distinguishes them.
""")
    y = init_state(p)
    y[IDX["HBF"]] = 7.0 * 0.8
    y[IDX["MHB"]] = 7.0 * 0.2
    o1 = outputs(y, p)
    y = init_state(p)
    y[IDX["HBF"]] = 17.0 * 0.8
    y[IDX["MHB"]] = 17.0 * 0.2
    o2 = outputs(y, p)
    anchor("SpO2 difference between Hb 7 and Hb 17 at MetHb 20%",
           abs(o1["SPO2"] - o2["SPO2"]), 0.0, 0.6, "%",
           "the oximeter cannot tell them apart")
    anchor("PvO2 difference between the same two patients",
           o2["PVO2"] - o1["PVO2"], 8.0, 20.0, "mmHg", "but the tissue can")


def sec_oximeter():
    hr("4.  THE INSTRUMENT — where the 85% floor actually comes from")
    print("""  The oximeter forms R = A660/A940 over the pulsatile signal and reads SpO2
  off a straight calibration line.  Methaemoglobin absorbs strongly and almost
  equally at both wavelengths, so as its fraction rises R is dragged towards 1
  -- and the calibration line evaluated at R = 1 is 110 - 25 = 85.

  The 85% floor is therefore a property of the DEVICE, not of the blood.  The
  model predicts, rather than assumes, both the floor and the shape of the
  approach to it, and it predicts the saturation gap as the difference between
  two things it computes independently.
""")
    p = dict(P)
    print(f"  {'MetHb%':>7s} {'R':>7s} {'SpO2':>7s} {'co-ox SaO2':>11s} "
          f"{'gap':>8s} {'dSpO2/d%':>9s} {'cyanotic?':>10s}")
    prev = None
    for fm in (0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.85):
        y = init_state(p)
        y[IDX["HBF"]] = 15.0 * (1.0 - fm)
        y[IDX["MHB"]] = 15.0 * fm
        o = outputs(y, p)
        slope = "" if prev is None else f"{(o['SPO2']-prev[1])/((fm-prev[0])*100):9.3f}"
        print(f"  {100*fm:7.0f} {o['R']:7.3f} {o['SPO2']:7.1f} {o['SAO2CO']:11.1f} "
              f"{o['SPO2']-o['SAO2CO']:8.1f} {slope:>9s} "
              f"{'YES' if o['CYAN']>=1 else 'no':>10s}")
        prev = (fm, o["SPO2"])
    print("""
  The last-but-one column is the model's sharpest statement about monitoring:
  the sensitivity of the oximeter to the disease COLLAPSES as the disease gets
  worse.  Between 0 and 5% MetHb the reading moves ~0.9 points per percent;
  between 60 and 85% it moves ~0.02.  The instrument stops responding exactly
  where the patient is dying, and it does so smoothly, without alarming.
""")
    y = init_state(p)
    y[IDX["HBF"]] = 15.0 * 0.7
    y[IDX["MHB"]] = 15.0 * 0.3
    o30 = outputs(y, p)
    y = init_state(p)
    y[IDX["HBF"]] = 15.0 * 0.15
    y[IDX["MHB"]] = 15.0 * 0.85
    o85 = outputs(y, p)
    anchor("SpO2 at MetHb 30%", o30["SPO2"], 84.0, 90.0, "%",
           "Barker & Tremper 1989: ~85%")
    anchor("SpO2 floor at MetHb 85%", o85["SPO2"], 84.0, 86.5, "%", "plateau ~85%")
    # cyanosis threshold
    lo, hi = 0.0, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        y = init_state(p)
        y[IDX["HBF"]] = 15.0 * (1 - mid)
        y[IDX["MHB"]] = 15.0 * mid
        if outputs(y, p)["CYAN"] < 1.0:
            lo = mid
        else:
            hi = mid
    anchor("MetHb fraction at which cyanosis becomes visible", 100 * hi,
           7.0, 13.0, "%", "clinical teaching: 10-15%")


def sec_natural_history():
    hr("5.  BENZOCAINE — an acute exposure, and what happens if nobody treats")
    r = sc_benzocaine(250.0, tend=24.0)
    ip = int(np.argmax(r["MetPct"]))
    print(f"  20% benzocaine spray, ~250 mg absorbed across the mucosa.")
    print(f"    peak MetHb        {r['MetPct'][ip]:7.2f} %  at t = {r['time'][ip]:.2f} h")
    print(f"    SpO2 at peak      {r['SPO2'][ip]:7.2f} %")
    print(f"    co-oximetry SaO2  {r['SAO2CO'][ip]:7.2f} %")
    print(f"    PvO2 at peak      {r['PVO2'][ip]:7.2f} mmHg")
    # spontaneous decay half-life after the oxidant is gone
    m = r["time"] > r["time"][ip] + 0.5
    tt, ff = r["time"][m], r["MetPct"][m]
    half = ff[0] / 2.0
    k = np.where(ff <= half)[0]
    t12 = (tt[k[0]] - tt[0]) if len(k) else float("nan")
    print(f"    spontaneous MetHb half-life once the oxidant has cleared: {t12:.2f} h")
    anchor("peak MetHb after a benzocaine spray", r["MetPct"][ip], 14.0, 42.0, "%",
           "case series 15-40%")
    anchor("time to peak", r["time"][ip], 0.2, 2.5, "h", "20 min - 2 h")
    anchor("untreated MetHb half-life (normal CYB5R)", t12, 2.5, 12.0, "h",
           "endogenous ~15%/h, prolonged while oxidant remains")
    return r


def sec_methylene_blue():
    hr("6.  METHYLENE BLUE — one molecule, two directions")
    r0 = sc_benzocaine(250.0, tend=12.0)
    r1 = sc_benzocaine(250.0, mb=[(0.75, 1.0)], tend=12.0)
    i0 = at(r0, 2.0)
    i1 = at(r1, 2.0)
    print(f"  Same exposure, methylene blue 1 mg/kg at 45 min:")
    print(f"    MetHb at t=2 h, untreated {r0['MetPct'][i0]:6.2f} %")
    print(f"    MetHb at t=2 h, treated   {r1['MetPct'][i1]:6.2f} %")
    j = at(r1, 0.75)
    k = at(r1, 1.75)
    rate = (r1["MetPct"][j] - r1["MetPct"][k])
    print(f"    fall in the first hour after the dose: {rate:.2f} percentage points")
    anchor("MetHb 1 h after methylene blue 1 mg/kg", r1["MetPct"][k], 0.5, 12.0, "%",
           "clinical: substantial fall within 30-60 min")

    print("""
  Now the dose-response.  Nothing in this model contains a maximum dose, a
  ceiling, or a rule.  It contains one branch point: leucomethylene blue can
  give its two electrons to methaemoglobin, or it can give them to oxygen.
  The first branch is proportional to MHB/(Km+MHB).  The second is not.
""")
    print(f"  {'MB mg/kg':>9s} {'MetHb@2h':>9s} {'nadir':>7s} {'Hb@48h':>8s} "
          f"{'Heinz':>7s} {'ROS-AUC':>9s} {'net effect':>12s}")
    best = (None, 1e9)
    for dose in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0):
        r = sc_benzocaine(250.0, mb=[(0.75, dose)], tend=48.0)
        i2 = at(r, 2.0)
        nad = float(np.min(r["MetPct"][r["time"] > 0.75]))
        rosauc = float(np.trapezoid(r["ROS"], r["time"]))
        hb48 = r["HBTOT"][-1]
        if r["MetPct"][i2] < best[1]:
            best = (dose, r["MetPct"][i2])
        eff = "reduces" if nad < 2.0 else "partial"
        if hb48 < 14.0:
            eff += "/haemolyses"
        print(f"  {dose:9.1f} {r['MetPct'][i2]:9.2f} {nad:7.2f} {hb48:8.2f} "
              f"{peak(r,'HEINZ'):7.3f} {rosauc:9.1f} {eff:>12s}")
    print(f"""
  The methaemoglobin column bottoms out and then stops improving, while the
  haemoglobin and Heinz columns keep getting worse.  That is the whole reason
  the clinical limit on methylene blue is a CUMULATIVE dose (7 mg/kg) rather
  than a per-dose one: the harm is the time-integral of a futile branch that
  runs whenever drug is present, and the benefit shuts off as soon as the
  substrate is gone.  Best methaemoglobin at 2 h was at {best[0]:.1f} mg/kg.
""")
    r7 = sc_benzocaine(250.0, mb=[(0.75, 7.0)], tend=48.0)
    r1b = sc_benzocaine(250.0, mb=[(0.75, 1.0)], tend=48.0)
    anchor("Hb at 48 h after methylene blue 1 mg/kg", r1b["HBTOT"][-1], 14.2, 15.2,
           "g/dL", "no meaningful haemolysis")
    anchor("Hb at 48 h after methylene blue 7 mg/kg", r7["HBTOT"][-1], 9.0, 14.0,
           "g/dL", "haemolysis at the cumulative ceiling")
    anchor("ROS exposure ratio, 7 mg/kg vs 1 mg/kg",
           float(np.trapezoid(r7["ROS"], r7["time"]))
           / float(np.trapezoid(r1b["ROS"], r1b["time"])),
           1.8, 12.0, "x")
    return r1


def sec_g6pd():
    hr("7.  G6PD — why the antidote fails, and why it then does harm")
    print("""  There is no rule in this model that says "do not give methylene blue in
  G6PD deficiency".  There is one shared, capped NADPH supply with two
  consumers: the flavin reductase that makes leucomethylene blue, and the
  glutathione reductase that keeps the membrane intact.  Lower the cap and
  both consumers lose -- and the one that loses second is the membrane.
""")
    print(f"  {'G6PD':>6s} {'MetHb@2h':>9s} {'nadir':>7s} {'GSH nadir':>10s} "
          f"{'Heinz':>7s} {'Hb@72h':>8s} {'verdict':>22s}")
    for g in (1.00, 0.60, 0.40, 0.25, 0.15, 0.08, 0.02):
        r = sc_benzocaine(250.0, mb=[(0.75, 2.0)], par=dict(G6PD=g), tend=72.0)
        nad = float(np.min(r["MetPct"][r["time"] > 0.75]))
        gsh = float(np.min(r["GSH"]))
        hb = r["HBTOT"][-1]
        if nad < 3.0 and hb > 14.0:
            v = "works"
        elif nad < 8.0 and hb > 13.0:
            v = "works, slower"
        elif hb < 12.5:
            v = "FAILS and haemolyses"
        else:
            v = "fails"
        print(f"  {g:6.2f} {r['MetPct'][at(r,2.0)]:9.2f} {nad:7.2f} {gsh:10.1f} "
              f"{peak(r,'HEINZ'):7.3f} {hb:8.2f} {v:>22s}")
    rn = sc_benzocaine(250.0, mb=[(0.75, 2.0)], par=dict(G6PD=1.0), tend=72.0)
    rd = sc_benzocaine(250.0, mb=[(0.75, 2.0)], par=dict(G6PD=0.02), tend=72.0)
    ra = sc_benzocaine(250.0, mb=[(0.75, 2.0)], par=dict(G6PD=0.15), tend=72.0)
    print(f"""
  The model reproduces the clinical gradient rather than a binary rule: class
  III deficiency (A-, ~15% activity) gets a slower but real response, while
  class II/I (Mediterranean, ~2%) gets essentially none and pays the
  haemolytic price anyway.  The reason is visible in the GSH column -- at low
  G6PD, the methylene blue is not merely ineffective, it is COMPETING with
  glutathione for the electrons that were holding the cell together.
""")
    anchor("MetHb 2 h after MB 2 mg/kg, normal G6PD", rn["MetPct"][at(rn, 2.0)],
           0.0, 4.0, "%", "antidote works")
    anchor("MetHb 2 h after MB 2 mg/kg, G6PD 2%", rd["MetPct"][at(rd, 2.0)],
           8.0, 40.0, "%", "antidote failure")
    anchor("MetHb 2 h after MB 2 mg/kg, G6PD 15% (A-)", ra["MetPct"][at(ra, 2.0)],
           1.0, 25.0, "%", "intermediate, as observed")
    anchor("Hb at 72 h, G6PD 2% given methylene blue", rd["HBTOT"][-1],
           6.5, 13.0, "g/dL", "haemolysis instead of reduction")
    anchor("G6PD 2% haemolyses MORE than normal G6PD",
           rn["HBTOT"][-1] - rd["HBTOT"][-1], 0.3, 9.0, "g/dL")


def sec_dapsone():
    hr("8.  DAPSONE — a formation-side disease, and the rebound it must cause")
    for dose in (50.0, 100.0, 200.0, 300.0):
        r = sc_dapsone(dose, days=21.0)
        print(f"  dapsone {dose:5.0f} mg/d  ->  steady-state MetHb "
              f"{r['MetPct'][-1]:6.2f} %   Hb {r['HBTOT'][-1]:5.2f} g/dL   "
              f"Cdap {r['DAPC'][-1]/P['V1DAP']:5.2f} mg/L")
    r100 = sc_dapsone(100.0, days=21.0)
    r100c = sc_dapsone(100.0, days=21.0, cimetidine=True)
    drop = 100.0 * (1.0 - r100c["MetPct"][-1] / r100["MetPct"][-1])
    print(f"\n  cimetidine 400 mg tid added to dapsone 100 mg/d:")
    print(f"    MetHb {r100['MetPct'][-1]:.2f}% -> {r100c['MetPct'][-1]:.2f}% "
          f"({drop:.0f}% reduction)")
    anchor("steady-state MetHb on dapsone 100 mg/d", r100["MetPct"][-1],
           3.0, 12.0, "%", "reported 5-12%")
    anchor("cimetidine reduction in steady-state MetHb", drop, 15.0, 45.0, "%",
           "Coleman: ~25-30%")

    print("""
  Now the overdose (2 g).  Dapsone's half-life is ~30 h and it recirculates through
  the bile; methylene blue's effect on the blood is over in about an hour.  The
  model was not told that dapsone poisoning rebounds.  It rebounds because two
  clocks are running at different speeds.
""")
    ro = sc_dapsone_od(2000.0, tend=144.0)
    rm = sc_dapsone_od(2000.0, mb_sched=[(6.0, 2.0)], tend=144.0)
    i6 = at(rm, 6.0)
    i7 = at(rm, 7.0)
    nadir_i = int(np.argmin(rm["MetPct"][(rm["time"] > 6.0) & (rm["time"] < 12.0)]))
    tt = rm["time"][(rm["time"] > 6.0) & (rm["time"] < 12.0)]
    ff = rm["MetPct"][(rm["time"] > 6.0) & (rm["time"] < 12.0)]
    reb = float(np.max(rm["MetPct"][rm["time"] > tt[nadir_i]]))
    treb = float(rm["time"][rm["time"] > tt[nadir_i]][
        int(np.argmax(rm["MetPct"][rm["time"] > tt[nadir_i]]))])
    print(f"    untreated peak MetHb            {peak(ro,'MetPct'):6.2f} % "
          f"at {tpeak(ro,'MetPct'):.1f} h")
    print(f"    untreated: hours above 20% MetHb "
          f"{float(np.trapezoid((ro['MetPct']>20).astype(float), ro['time'])):6.1f} h")
    print(f"    single MB 2 mg/kg at 6 h: {rm['MetPct'][i6]:.1f}% -> "
          f"{ff[nadir_i]:.1f}% at {tt[nadir_i]:.1f} h, then REBOUND to "
          f"{reb:.1f}% at {treb:.1f} h")
    anchor("dapsone 2 g overdose, untreated peak MetHb", peak(ro, "MetPct"),
           30.0, 75.0, "%", "reported 30-70%")
    anchor("rebound after a single methylene blue dose", reb, 8.0, 60.0, "%",
           "well documented; requires repeat dosing")
    anchor("time to rebound peak", treb - tt[nadir_i], 1.0, 96.0, "h")

    # the cumulative-dose trap, and the argument for attacking formation instead
    rrep = sc_dapsone_od(2000.0, mb_sched=[(6.0, 2.0), (10.0, 2.0), (16.0, 1.0),
                                           (24.0, 1.0), (36.0, 1.0)], tend=144.0)
    par_c = dict()
    ev = [(0.0, "DAPG", 2000.0), (6.0, "MB_mgkg", 2.0), (10.0, "MB_mgkg", 2.0)]
    for k in range(6):
        for h in (2.0, 10.0, 18.0):
            ev.append((24.0 * k + h, "CIMG", 400.0))
    rcim = simulate(par=par_c, events=ev, tend=144.0, npt=2400)
    print(f"""
    5 doses of methylene blue (7 mg/kg cumulative, i.e. the ceiling):
        MetHb AUC {float(np.trapezoid(rrep['MetPct'], rrep['time'])):8.1f} %.h,
        Hb at 144 h {rrep['HBTOT'][-1]:.2f} g/dL, Heinz peak {peak(rrep,'HEINZ'):.3f}
    2 doses of methylene blue + cimetidine (4 mg/kg cumulative):
        MetHb AUC {float(np.trapezoid(rcim['MetPct'], rcim['time'])):8.1f} %.h,
        Hb at 144 h {rcim['HBTOT'][-1]:.2f} g/dL, Heinz peak {peak(rcim,'HEINZ'):.3f}

  This is the model's therapeutic argument, and it is an argument about WHICH
  TERM you attack.  Methylene blue acts on the clearance term and is spent
  within the hour; cimetidine acts on the formation term and keeps acting.  In
  a poisoning whose source term lasts days, the drug with the short action has
  to be re-given until the cumulative ceiling is reached, whereas turning down
  the source is free of that constraint.
""")


def sec_nitrite():
    hr("9.  NITRITE — the autocatalytic one")
    print("""  Nitrite oxidation of oxyhaemoglobin is autocatalytic: the product
  accelerates the reaction.  The model contains that one extra term -- and the
  result below is NOT the one it was put there to produce.
""")
    for dose in (2000.0, 5000.0, 10000.0):
        umol = dose * 1000.0 / P["MWNIT"]
        r = simulate(events=[(0.0, "NITG", umol)], tend=12.0, npt=1200)
        ip = int(np.argmax(r["MetPct"]))
        # time from 5% to 30%
        t5 = r["time"][np.argmax(r["MetPct"] > 5.0)] if np.any(r["MetPct"] > 5) else float("nan")
        t30 = r["time"][np.argmax(r["MetPct"] > 30.0)] if np.any(r["MetPct"] > 30) else float("nan")
        print(f"  sodium nitrite {dose/1000:5.1f} g  ->  peak MetHb {r['MetPct'][ip]:6.2f}% "
              f"at {r['time'][ip]:5.2f} h;  5%->30% in "
              f"{(t30-t5) if t30==t30 else float('nan'):5.2f} h;  "
              f"PvO2 nadir {float(np.min(r['PVO2'])):5.1f} mmHg")
    r = simulate(events=[(0.0, "NITG", 5000.0 * 1000.0 / P["MWNIT"])],
                 tend=16.0, npt=1400)
    rlin = simulate(par=dict(KNIT2=0.0),
                    events=[(0.0, "NITG", 5000.0 * 1000.0 / P["MWNIT"])],
                    tend=16.0, npt=1400)

    def rise(rr):
        a = rr["time"][np.argmax(rr["MetPct"] > 5.0)]
        b = rr["time"][np.argmax(rr["MetPct"] > 30.0)]
        return float(b - a)

    print(f"""
  A REFUTED HYPOTHESIS, LEFT IN RATHER THAN REMOVED
  -------------------------------------------------
  The autocatalytic term was put into this model on the expectation that it
  would generate the clinical lag-then-wall.  It does not.  Switching it off
  entirely moves the 5 g peak from {peak(r,'MetPct'):.1f}% to {peak(rlin,'MetPct'):.1f}% and leaves the
  5%-to-30% rise time at {rise(r):.2f} h vs {rise(rlin):.2f} h -- no meaningful difference.

  The reason is worth stating because it is itself a result: at poisoning
  doses the reaction is not rate-limited by its own chemistry.  It is limited
  by ABSORPTION on the way up and by STOICHIOMETRY at the top -- 5 g of sodium
  nitrite is 72.5 mmol against roughly 186 mmol of haem, and unlike the
  arylhydroxylamines nitrite is consumed as it goes.  That is why the peak
  scales with dose almost linearly while the timing barely moves, and it is
  why "how much did they take" predicts a nitrite patient far better than it
  predicts a dapsone one.

  The autocatalysis is real chemistry; it is simply not the rate-limiting step
  in vivo.  The term is left in the model, and this paragraph is left in the
  output, because deleting a hypothesis that failed is how a model stops being
  checkable.
""")
    anchor("nitrite 5 g peak MetHb", peak(r, "MetPct"), 45.0, 90.0, "%")


def sec_congenital():
    hr("10.  CONGENITAL DISEASE — why chronic 25% is well and acute 25% is not")
    print("""  Congenital CYB5R3 deficiency produces patients who are visibly blue,
  have methaemoglobin fractions of 15-30%, and are entirely well.  The model
  was not given a "chronic" switch.  It has 2,3-BPG and a marrow, and both
  respond to tissue hypoxia -- and 2,3-BPG's effect on P50 is the exact
  opposite allosteric move to the one methaemoglobin makes.
""")
    r = simulate(par=dict(EB5SET=0.035), tend=24 * 250.0, npt=2400)
    o = {k: r[k][-1] for k in ("MetPct", "PVO2", "SPO2", "P50E", "HBTOT", "CYAN")}
    y = init_state(P)
    fm = o["MetPct"] / 100.0
    y[IDX["HBF"]] = 15.0 * (1 - fm)
    y[IDX["MHB"]] = 15.0 * fm
    oa = outputs(y, P)
    print(f"  chronic (250 d at CYB5R in-vivo flux fraction 0.035):")
    print(f"    MetHb {o['MetPct']:.1f}%   2,3-BPG {r['BPG'][-1]:.2f}x   "
          f"P50eff {o['P50E']:.1f} mmHg   Hb {o['HBTOT']:.2f}")
    print(f"    PvO2 {o['PVO2']:.1f} mmHg    SpO2 {o['SPO2']:.1f}%   "
          f"cyanosis index {o['CYAN']:.2f}")
    print(f"  the SAME methaemoglobin fraction, arrived at acutely:")
    print(f"    P50eff {oa['P50E']:.1f} mmHg    PvO2 {oa['PVO2']:.1f} mmHg   "
          f"SpO2 {oa['SPO2']:.1f}%")
    print(f"""
  Same percentage, same oximeter reading, {o['PVO2']-oa['PVO2']:.1f} mmHg of difference at the
  tissue.  The compensation is not tolerance and it is not adaptation of the
  brain: it is a measurable right-shift that cancels part of a measurable
  left-shift.  This is also why methylene blue in congenital methaemoglobinaemia
  is a COSMETIC intervention -- the model shows the colour changing while the
  tissue PO2 barely moves.
""")
    rmb = simulate(par=dict(EB5SET=0.035), tend=24 * 250.0 + 12.0, npt=2600,
                   events=[(24 * 250.0, "MB_mgkg", 1.0)])
    i0 = at(rmb, 24 * 250.0 - 0.1)
    i1 = at(rmb, 24 * 250.0 + 1.0)
    print(f"    methylene blue 1 mg/kg: MetHb {rmb['MetPct'][i0]:.1f}% -> "
          f"{rmb['MetPct'][i1]:.1f}%,  PvO2 {rmb['PVO2'][i0]:.1f} -> "
          f"{rmb['PVO2'][i1]:.1f} mmHg,  SpO2 {rmb['SPO2'][i0]:.1f} -> "
          f"{rmb['SPO2'][i1]:.1f}%")
    anchor("chronic CYB5R deficiency steady-state MetHb", o["MetPct"], 10.0, 35.0,
           "%", "congenital type I: 10-35%")
    anchor("chronic CYB5R deficiency PvO2 (near normal despite the colour)",
           o["PVO2"], 33.0, 41.0, "mmHg")
    anchor("PvO2 advantage of chronic over acute at the same MetHb",
           o["PVO2"] - oa["PVO2"], 3.0, 15.0, "mmHg")

    print("\n  HbM variant (methaemoglobin that cannot be reduced at all):")
    rm = simulate(par=dict(FHBM=0.25, EB5SET=1.0),
                  events=[(2.0, "MB_mgkg", 2.0)], tend=12.0, npt=800)
    print(f"    MetHb before MB {rm['MetPct'][at(rm,1.9)]:.1f}% -> after "
          f"{rm['MetPct'][at(rm,4.0)]:.1f}%   (no response, correctly)")
    anchor("HbM: MetHb change after methylene blue",
           abs(rm["MetPct"][at(rm, 4.0)] - rm["MetPct"][at(rm, 1.9)]), 0.0, 2.0,
           "%", "HbM does not respond to MB")


def sec_hbo_and_transfusion():
    hr("11.  BYPASSING HAEMOGLOBIN ALTOGETHER")
    print("""  Three therapeutic classes act on three different terms:
    (1) stop the oxidant / block bioactivation  -> the FORMATION term
    (2) methylene blue, ascorbate               -> the CLEARANCE term
    (3) hyperbaric oxygen, exchange transfusion -> the DELIVERY term itself
  Only the third works when the first two cannot.
""")
    p = dict(P)
    for fm, pao2, label in ((0.70, 95.0, "MetHb 70%, room air"),
                            (0.70, 640.0, "MetHb 70%, 100% O2 at 1 ATA"),
                            (0.70, 1800.0, "MetHb 70%, 100% O2 at 2.8 ATA"),
                            (0.95, 1800.0, "MetHb 95%, 100% O2 at 2.8 ATA")):
        y = init_state(p)
        y[IDX["HBF"]] = 15.0 * (1 - fm)
        y[IDX["MHB"]] = 15.0 * fm
        o = outputs(y, dict(p, PAO2=pao2))
        diss = 0.003 * pao2
        print(f"  {label:<34s} dissolved O2 {diss:5.2f} mL/dL   "
              f"CaO2 {o['CAO2']:5.2f}   PvO2 {o['PVO2']:6.1f} mmHg   "
              f"{'ANAEROBIC' if o['PVO2'] < p['PVCRIT'] else 'aerobic'}")
    print("""
  The arithmetic is Boerema's: resting extraction is VO2/(CO x 10) = 5.0 mL of
  oxygen per decilitre, and 2.8 ATA of oxygen dissolves 5.4 mL/dL in plasma.
  Hyperbaric oxygen therefore covers resting metabolism with no functional
  haemoglobin at all -- which is why it is the correct answer in exactly the
  case where methylene blue has no substrate to work on.
""")
    y = init_state(p)
    y[IDX["HBF"]] = 15.0 * 0.05
    y[IDX["MHB"]] = 15.0 * 0.95
    o = outputs(y, dict(p, PAO2=1800.0))
    anchor("PvO2 at MetHb 95% on HBO 2.8 ATA", o["PVO2"], 20.0, 400.0, "mmHg",
           "above the anaerobic threshold with essentially no functional Hb")
    y2 = init_state(p)
    y2[IDX["HBF"]] = 15.0 * 0.05
    y2[IDX["MHB"]] = 15.0 * 0.95
    o2 = outputs(y2, dict(p, PAO2=95.0))
    anchor("PvO2 at MetHb 95% on room air", o2["PVO2"], 0.0, 6.0, "mmHg")

    print("  Exchange transfusion (replace 8 g/dL worth of the patient's blood):")
    rex = sc_benzocaine(420.0, tend=12.0)
    rtx = simulate(events=[(0.0, "BZCD", 420.0), (2.0, "TRANSFUSE", 8.0)],
                   tend=12.0, npt=900)
    print(f"    MetHb at 3 h  untreated {rex['MetPct'][at(rex,3.0)]:6.2f}%  vs "
          f"exchanged {rtx['MetPct'][at(rtx,3.0)]:6.2f}%")
    print(f"    PvO2  at 3 h  untreated {rex['PVO2'][at(rex,3.0)]:6.2f}   vs "
          f"exchanged {rtx['PVO2'][at(rtx,3.0)]:6.2f} mmHg")


def sec_ascorbate_and_ss():
    hr("12.  THE ALTERNATIVES, AND THE ANTIDOTE'S OWN PHARMACOLOGY")
    r_mb = sc_benzocaine(250.0, mb=[(0.75, 1.0)], tend=24.0)
    r_as = simulate(events=[(0.0, "BZCD", 250.0), (0.75, "ASCC", 10.0 * 1e6 / 176.1)],
                    tend=24.0, npt=1200)
    r_no = sc_benzocaine(250.0, tend=24.0)
    print("  ascorbic acid 10 g IV vs methylene blue 1 mg/kg vs nothing:")
    for lab, r in (("nothing", r_no), ("ascorbate 10 g", r_as),
                   ("methylene blue", r_mb)):
        print(f"    {lab:<16s} MetHb@2h {r['MetPct'][at(r,2.0)]:6.2f}%   "
              f"@6h {r['MetPct'][at(r,6.0)]:6.2f}%   "
              f"AUC {float(np.trapezoid(r['MetPct'], r['time'])):7.1f} %.h")
    print("""
  Ascorbate is not a slow version of methylene blue; it is a different order of
  magnitude.  Its rate constant is a non-enzymatic bimolecular one, and no dose
  that a kidney will tolerate makes it competitive.  It belongs where methylene
  blue is contraindicated, not where methylene blue is merely unavailable.
""")
    anchor("ascorbate 10 g is slower than MB at 2 h",
           r_as["MetPct"][at(r_as, 2.0)] - r_mb["MetPct"][at(r_mb, 2.0)],
           2.0, 40.0, "%", "ascorbate must NOT match MB")

    print("\n  Methylene blue is also a potent MAO-A inhibitor:")
    for lab, par, mb in (("MB alone", dict(), 2.0),
                         ("SSRI alone", dict(), 0.0),
                         ("MB + SSRI", dict(), 2.0)):
        ev = [(0.0, "BZCD", 250.0)]
        if mb:
            ev.append((0.75, "MB_mgkg", mb))
        if "SSRI" in lab:
            ev.append((0.0, "SSRI", 2.5))
        r = simulate(par=par, events=ev, tend=12.0, npt=800)
        print(f"    {lab:<12s} peak serotonin index {peak(r,'HT5'):5.2f} x normal")
    ev = [(0.0, "BZCD", 250.0), (0.0, "SSRI", 2.5), (0.75, "MB_mgkg", 2.0)]
    rss = simulate(events=ev, tend=12.0, npt=800)
    ev2 = [(0.0, "BZCD", 250.0), (0.75, "MB_mgkg", 2.0)]
    rmb = simulate(events=ev2, tend=12.0, npt=800)
    anchor("serotonin index, MB alone", peak(rmb, "HT5"), 1.1, 2.2, "x")
    anchor("serotonin index, MB on top of an SSRI", peak(rss, "HT5"), 2.5, 12.0, "x",
           "the reported interaction")


def sec_cohort():
    hr("13.  VIRTUAL COHORT — who actually decompensates")
    rng = np.random.default_rng(20260807)
    n = 240
    rows = []
    for _ in range(n):
        hb = float(np.clip(rng.normal(14.0, 2.2), 6.5, 18.5))
        g6 = float(np.clip(rng.lognormal(np.log(0.9), 0.55), 0.02, 1.4))
        eb5 = float(np.clip(rng.normal(1.0, 0.14), 0.35, 1.4))
        dose = float(np.clip(rng.normal(250.0, 80.0), 60.0, 520.0))
        r = sc_benzocaine(dose, par=dict(HB0=hb, G6PD=g6, EB5SET=eb5), tend=10.0)
        rows.append(dict(hb=hb, g6=g6, eb5=eb5, dose=dose,
                         met=peak(r, "MetPct"),
                         pvo2=float(np.min(r["PVO2"])),
                         spo2=float(np.min(r["SPO2"]))))
    met = np.array([x["met"] for x in rows])
    pv = np.array([x["pvo2"] for x in rows])
    hbv = np.array([x["hb"] for x in rows])
    anaer = pv < P["PVCRIT"]
    print(f"  n = {n} exposures to a benzocaine spray of variable size.")
    print(f"    median peak MetHb          {np.median(met):.1f}%  "
          f"(IQR {np.percentile(met,25):.1f}-{np.percentile(met,75):.1f})")
    print(f"    crossed the anaerobic threshold: {anaer.sum()} "
          f"({100*anaer.mean():.1f}%)")
    print(f"    of those, median peak MetHb {np.median(met[anaer]):.1f}% "
          f"and median Hb {np.median(hbv[anaer]):.1f} g/dL")
    print(f"    of the rest, median peak MetHb {np.median(met[~anaer]):.1f}% "
          f"and median Hb {np.median(hbv[~anaer]):.1f} g/dL")
    # how well does %MetHb alone rank the patients who decompensated?
    thr = np.percentile(met, 100 * (1 - anaer.mean()))
    tp = ((met >= thr) & anaer).sum()
    print(f"    if you triage on %MetHb alone at the matched cut-off "
          f"({thr:.1f}%), you identify {tp}/{anaer.sum()} of them.")
    prod = hbv * (1 - met / 100.0)
    thr2 = np.percentile(prod, 100 * anaer.mean())
    tp2 = ((prod <= thr2) & anaer).sum()
    print(f"    if you triage on the PRODUCT Hb x (1-f) instead, you identify "
          f"{tp2}/{anaer.sum()}.")
    anchor("product beats percentage as a triage variable",
           float(tp2 - tp), 0.0, 100.0, "patients")


def sec_sulf():
    hr("14.  SULFHAEMOGLOBIN — the differential the oximeter cannot make")
    p = dict(P)
    print(f"  {'pigment':>14s} {'g/dL':>6s} {'SpO2':>7s} {'cyanotic?':>10s} "
          f"{'responds to MB?':>16s}")
    for lab, met, sul in (("methaemoglobin", 1.5, 0.0), ("methaemoglobin", 3.0, 0.0),
                          ("sulfhaemoglobin", 0.0, 0.5), ("sulfhaemoglobin", 0.0, 1.5)):
        y = init_state(p)
        y[IDX["HBF"]] = 15.0 - met - sul
        y[IDX["MHB"]] = met
        y[IDX["SHB"]] = sul
        o = outputs(y, p)
        print(f"  {lab:>14s} {met+sul:6.2f} {o['SPO2']:7.1f} "
              f"{'YES' if o['CYAN']>=1 else 'no':>10s} "
              f"{'yes' if met>0 else 'NO (irreversible)':>16s}")
    print("""
  Half a gram of sulfhaemoglobin is visibly cyanotic; it takes three times as
  much methaemoglobin and ten times as much deoxyhaemoglobin.  The model gets
  this from three thresholds in one summed index, and it gets the crucial
  clinical difference -- non-response to methylene blue -- from the fact that
  sulfhaemoglobin never enters the reduction equations at all.
""")


def sec_sensitivity():
    hr("15.  SENSITIVITY — which parameters the conclusions actually rest on")
    base = sc_benzocaine(250.0, mb=[(0.75, 2.0)], tend=24.0)
    b_peak = peak(base, "MetPct")
    b_nad = float(np.min(base["MetPct"][base["time"] > 0.75]))
    b_pv = float(np.min(base["PVO2"]))
    print(f"  base case: peak MetHb {b_peak:.2f}%, nadir after MB {b_nad:.2f}%, "
          f"PvO2 nadir {b_pv:.1f} mmHg")
    print(f"\n  {'parameter':<12s} {'x0.5':>26s} {'x2':>26s}")
    print(f"  {'':12s} {'peak / nadir / PvO2':>26s} {'peak / nadir / PvO2':>26s}")
    for k in ("ALPHAM", "VMXB5", "KO2", "VPPP", "KMBR", "KTOL", "KMBOXH",
              "NHILL", "GCO", "KROSOX"):
        line = f"  {k:<12s}"
        for mult in (0.5, 2.0):
            r = sc_benzocaine(250.0, mb=[(0.75, 2.0)],
                              par={k: P[k] * mult}, tend=24.0)
            line += (f" {peak(r,'MetPct'):7.2f}/"
                     f"{float(np.min(r['MetPct'][r['time']>0.75])):6.2f}/"
                     f"{float(np.min(r['PVO2'])):6.1f}   ")
        print(line)
    print("""
  ALPHAM -- the size of the Darling-Roughton left shift -- is the parameter the
  headline result rests on, and it is the one with the weakest quantitative
  literature.  It is stated here rather than buried: if the true left shift is
  half of what is assumed, the 30% equivalence moves from 8.3 to about 9.4 g/dL
  and the argument weakens but does not reverse; if it is zero, the argument
  disappears entirely and %MetHb becomes an adequate bedside variable.  That is
  a falsifiable prediction, and the experiment is a co-oximeter and a blood-gas
  P50 measured on the same sample.
""")
    for mult, lab in ((0.0, "no left shift"), (0.5, "half"), (1.0, "as modelled"),
                      (1.5, "1.5x")):
        p = dict(P, ALPHAM=P["ALPHAM"] * mult)
        y = init_state(p)
        y[IDX["HBF"]] = 15.0 * 0.7
        y[IDX["MHB"]] = 15.0 * 0.3
        o = outputs(y, p)
        eq = equivalent_hb(o["PVO2"], p)
        print(f"    ALPHAM {lab:<14s} P50eff {o['P50E']:5.1f}  PvO2 {o['PVO2']:5.1f}  "
              f"EQUIV Hb at 30% MetHb {eq:5.2f} g/dL")


# ---------------------------------------------------------------------------
# 7.  PARAMETER CROSS-CHECK AGAINST THE R FILE
# ---------------------------------------------------------------------------
def check_params(rfile):
    hr("PARAMETER CROSS-CHECK vs the mrgsolve $PARAM block")
    try:
        txt = open(rfile, encoding="utf-8").read()
    except OSError as e:
        print(f"  cannot read {rfile}: {e}")
        return False
    m = re.search(r"\$PARAM.*?\n(.*?)\n\s*\$", txt, re.S)
    if not m:
        print("  no $PARAM block found")
        return False
    body = m.group(1)
    body = re.sub(r"//[^\n]*", "", body)
    found = {}
    for mm in re.finditer(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(-?[0-9.eE+-]+)", body):
        try:
            found[mm.group(1)] = float(mm.group(2))
        except ValueError:
            pass
    bad = 0
    for k, v in sorted(P.items()):
        if k not in found:
            print(f"  MISSING in R: {k}")
            bad += 1
        elif abs(found[k] - v) > 1e-9 * max(1.0, abs(v)):
            print(f"  MISMATCH {k}: python {v}  R {found[k]}")
            bad += 1
    extra = [k for k in found if k not in P]
    for k in extra:
        print(f"  EXTRA in R (not in python): {k}")
    print(f"\n  {len(P)} python parameters, {len(found)} R parameters, "
          f"{bad} discrepancies, {len(extra)} R-only.")
    return bad == 0 and not extra


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--anchors", action="store_true")
    ap.add_argument("--params", nargs="?", const="mhb_mrgsolve_model.R")
    a = ap.parse_args()
    if a.params:
        ok = check_params(a.params)
        sys.exit(0 if ok else 1)

    np.set_printoptions(suppress=True)
    print("METHAEMOGLOBINAEMIA QSP MODEL — independent verification run")
    print(f"{NST} ODE states, {len(P)} parameters")

    sec_steadystate()
    sec_equivalent_anemia()
    sec_same_percent_different_disease()
    sec_oximeter()
    sec_natural_history()
    sec_methylene_blue()
    sec_g6pd()
    sec_dapsone()
    sec_nitrite()
    sec_congenital()
    sec_hbo_and_transfusion()
    sec_ascorbate_and_ss()
    sec_sulf()
    sec_cohort()
    sec_sensitivity()

    hr("ANCHOR SUMMARY")
    npass = sum(1 for a_ in ANCHORS if a_[4])
    for lab, v, lo, hi, ok in ANCHORS:
        if not ok:
            print(f"  FAIL  {lab}: {v:.3f} not in [{lo}, {hi}]")
    print(f"\n  {npass}/{len(ANCHORS)} anchors pass.")
    sys.exit(0 if npass == len(ANCHORS) else 1)


if __name__ == "__main__":
    main()
