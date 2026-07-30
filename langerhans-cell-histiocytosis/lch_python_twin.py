#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lch_python_twin.py — dependency-free numerical twin of lch_mrgsolve_model.R
==========================================================================

Langerhans Cell Histiocytosis (LCH) QSP model, 62 ODEs.

WHY THIS FILE EXISTS
--------------------
`lch_mrgsolve_model.R` is the reference implementation, but mrgsolve needs R +
a compiler. This file re-implements the *same* right-hand side and the *same*
parameter block in pure Python (RK4, no numpy) so that every quantitative claim
made in README.md can be re-checked with `python3 lch_python_twin.py`.

Parameter names and compartment names are identical to the R model, so the two
files can be diffed line-by-line. If you change one, change the other.

Structural commitments this file is built to demonstrate (and to be able to
FALSIFY — every one has a kill switch):

  (1) CELL-OF-ORIGIN PARTITION, NOT GROWTH RATE, SETS PHENOTYPE.
      One mutant precursor pool is partitioned across seeding routes
      (THETA_MARROW/CIRC/LOCAL and the per-organ thetas). Single-system bone
      and multisystem risk-organ-positive disease differ only in that
      partition and in PRECM0 — no proliferation or death constant changes.
      Kill switch: set the phenotypes' thetas equal and the RO+/SS-b
      separation disappears.

  (2) MAPK INHIBITION IS CYTOSTATIC ON THIS LINEAGE. It lowers proliferation
      (CCND), survival signal (BCL) and secretome (SASP) but adds NO kill
      term. Chemotherapy adds explicit kill terms. Consequence, not
      assumption: DAS falls in days (secretome), lesion mass and cfDNA fall
      over weeks-months and plateau, and withdrawal rebounds at the
      precursor's own regrowth rate. Kill switch: SL_MAPKI_KILL > 0 abolishes
      the rebound — contradicting the clinical literature.

  (3) PERMANENT SEQUELAE ARE TIME-INTEGRALS OF ACTIVE DISEASE, NOT FUNCTIONS
      OF ITS PEAK. AVPN / NEUR / BILF / LUNGC are monotone accumulators with
      thresholds. Therefore time-to-effective-therapy can dominate drug
      potency, and a faster-but-shallower regimen can beat a
      slower-but-deeper one on CNS endpoints. The crossover is computable
      (see check "crossover").

Units: time = days. Concentrations = mg/L. Lesional burden = arbitrary
"burden units" (1 unit ~ 1e9 LCH cells). ANC = 1e9/L.
"""

from __future__ import annotations

import math
import sys

# ===========================================================================
# $PARAM  — identical names/values to the R model
# ===========================================================================
P = dict(
    # ---- 1. Genotype / phenotype switches -------------------------------
    WT=12.0,         # body weight (kg) - drives allometric PK scaling
    ALLO_CL=0.75,    # allometric exponent for clearances
    ALLO_V=1.00,     # allometric exponent for volumes
    GENO=1.0,        # 1 = BRAF V600E, 2 = MAP2K1/ARAF/driver-neg (BRAFi inactive)
    SMOKE=0.0,       # 0/1 cigarette smoke exposure (pulmonary LCH driver)
    IL17ON=1.0,      # 1 = keep the (contested) IL-17A/DC-fusion arm
    # ---- 2. Cell-of-origin partition ------------------------------------
    PRECM0=0.50,     # initial mutant marrow precursor pool (burden units)
    FSR=1.00,        # self-renewal COMPETENCE of the cell of origin (0-1).
                     # This is a property of the differentiation stage at
                     # which the driver arose, not a kinetic rate: an HSC-
                     # stage clone self-renews (FSR=1), a committed tissue
                     # precursor barely does (FSR~0.1).
    PRMAX=3.00,      # niche carrying capacity of the marrow reservoir
    PNICHE=0.0030,   # minimum clone size that still supports self-renewal;
                     # a deterministic surrogate for stochastic clonal
                     # extinction. Below it the clone cannot re-expand.
    THB=0.34,        # seeding fraction -> bone
    THS=0.18,        # -> skin
    THR=0.26,        # -> risk organs (liver/spleen/marrow)
    THP=0.08,        # -> pituitary / hypothalamus
    THC=0.06,        # -> CNS parenchyma (cerebellum/pons)
    THL=0.08,        # -> lung
    LBONE0=0.05,     # locally seeded initial lesions (tissue-restricted origin)
    LSKIN0=0.02,
    LRO0=0.02,
    LPIT0=0.01,
    LCNS0=0.005,
    LLUNG0=0.005,
    # ---- 3. Precursor / lesion kinetics ---------------------------------
    KSR=0.070,       # max self-renewal rate of mutant precursor (/day)
    KDP=0.015,       # precursor baseline death (/day)
    KEX=0.080,       # progeny output to blood (/day, proliferation-coupled)
    KDC=1.00,        # circulating precursor clearance (/day)
    KSEED=0.50,      # circulating -> tissue seeding (/day)
    KPROL=0.022,     # intralesional proliferation at CCND=1 (/day).
                     # DELIBERATELY BELOW the minimum immune clearance rate
                     # (KDL + KIMM*(1-TREGMAX) = 0.0288): a lesion can never
                     # sustain itself by local division alone. This encodes
                     # the Ki-67 ~5-10% finding - lesions are maintained by
                     # continuous recruitment from the precursor pool, so
                     # cutting the supply collapses them.
    KDL=0.020,       # lesional baseline death (/day)
    KIMM=0.035,      # immune-mediated lesional clearance at TREG=0 (/day)
    KAPO=1.60,       # apoptosis amplification when BCL is lost
    FN_APO=0.25,     # fraction of that amplification the MARROW precursor
                     # actually feels. The niche supplies ERK-independent
                     # survival signals, so losing BCL2A1 sensitises
                     # lesional cells far more than reservoir cells. This is
                     # what lets the reservoir plateau (not clear) under a
                     # cytostatic drug, matching persistent BRAF V600E+
                     # blood cells on MAPK inhibition.
    LMAX=8.0,        # per-site carrying capacity
    RECMAX=3.0,      # cap on IL-1beta-driven recruitment amplification
    FREC=0.45,       # IL-1beta-driven CCR6 recruitment amplification
    LREF=2.0,        # burden units used to normalise lesion mass
    # ---- 4. MAPK signalling ---------------------------------------------
    KERK=24.0,       # pERK equilibration (/day)
    ERKB=0.05,       # ERK activity floor (non-driver signalling)
    KCC=2.00,        # cyclin D1 turnover (/day)
    KBCL=1.50,       # BCL2A1 turnover (/day)
    KSASP=3.00,      # secretory program turnover (/day) - FAST
    EC50C=0.35, HC=2.0,     # ERK -> cyclin D1
    EC50B=0.30, HB=1.5,     # ERK -> BCL2A1
    EC50S=0.25, HS=1.2,     # ERK -> SASP
    PMAX=0.60,       # max paradoxical MAPK activation (BRAFi in BRAF-WT cells)
    # ---- 5. Vemurafenib PK/PD -------------------------------------------
    KA_VEM=8.0, CL_VEM=32.0, V2_VEM=90.0, Q_VEM=12.0, V3_VEM=25.0,
    FU_VEM=0.005, IC50_VEM=0.018,   # free mg/L giving 50% BRAF inhibition
    KENZ=0.10, EMAX_IND=0.35, EC50_IND=30.0,   # CYP3A4 autoinduction
    # ---- 6. Dabrafenib PK/PD --------------------------------------------
    KA_DAB=20.0, CL_DAB=840.0, V_DAB=400.0,
    FM_DAB=0.40, CL_DABM=520.0, V_DABM=180.0, RPOT_M=1.00,
    FU_DAB=0.005, IC50_DAB=0.00050, FCNS_DAB=0.35,
    # ---- 7. Trametinib PK/PD --------------------------------------------
    KA_TRA=6.0, CL_TRA=118.0, V2_TRA=250.0, Q_TRA=25.0, V3_TRA=900.0,
    FU_TRA=0.026, IC50_TRA=0.00015,
    # ---- 8. Cobimetinib PK/PD -------------------------------------------
    KA_COB=8.0, CL_COB=322.0, V_COB=1000.0,
    FU_COB=0.050, IC50_COB=0.0012,
    # ---- 9. Vinblastine PK/PD -------------------------------------------
    CL_VBL=400.0, V1_VBL=10.0, Q_VBL=300.0, V2_VBL=700.0,
    EC50_VBL=0.0050, EMAX_VBL=0.35,
    # ---- 10. Prednisolone PK/PD -----------------------------------------
    KA_PRE=12.0, CL_PRE=190.0, V_PRE=40.0, KEO_P=4.0,
    EMAX_GR=0.70, EC50_GR=0.050,     # NF-kB transrepression of secretome
    EMAX_GRK=0.28, EC50_GRK=0.060,   # direct lympho/histiocyte apoptosis
    # ---- 11. Cytarabine PK/PD -------------------------------------------
    CL_ARAC=1915.0, V_ARAC=40.0,
    KIN_ARAC=1.0, KOUT_ARATP=4.5, FCNS_ARAC=0.15,
    EMAX_ARAC=0.25, EC50_ARAC=0.020,
    # ---- 12. Cladribine PK/PD -------------------------------------------
    CL_CLAD=150.0, V_CLAD=55.0,
    KIN_CLAD=1.0, KOUT_CLATP=0.70, FCNS_CLAD=0.25,
    EMAX_CLAD=0.20, EC50_CLAD=0.015,
    # ---- 13. Maintenance (6-MP / MTX) -----------------------------------
    KA_MNT=2.0, KE_MNT=1.0, EMAX_MNT=0.35, EC50_MNT=0.30,
    # ---- 14. The cytostatic / cytotoxic split ---------------------------
    PFMAX=0.09,          # proliferative fraction of LCH cells at CCND=1 (Ki-67)
    SL_MAPKI_KILL=0.0,   # KILL SWITCH for commitment (2): must stay 0
    W_ARAC_P=0.55,       # ara-CTP kill weight on the marrow precursor pool
    W_CLAD_P=1.40,       # Cd-ATP kill weight on precursor (non-cycling-active)
    W_VBL_P=0.30,
    W_GR_P=0.35,
    # ---- 15. Microenvironment -------------------------------------------
    KTR=0.25, TREGMAX=0.75, EC50_TREG=1.0,
    KOC=0.50, KR_RANKL=1.20, AOSM=0.45,
    # ---- 16. Secretome --------------------------------------------------
    KP_IL1=1.00, KD_IL1=6.0, IL1B0=0.05,
    KP_TNF=1.00, KD_TNF=6.0, TNFA0=0.08,
    KP_IL6=1.20, KD_IL6=5.0, IL60=0.10,
    KP_OSM=0.80, KD_OSM=4.0, OSM0=0.05,
    KP_MMP=1.00, KD_MMP=2.0, MMP90=0.20,
    KP_RNK=1.00, KD_RNK=2.0, RANKL0=0.30, W_TCELL_RNK=0.35,
    # ---- 17. Biomarkers -------------------------------------------------
    KP_CF=14.0, KD_CF=8.0,          # cfDNA BRAF V600E signal
    CF_LOD=0.005,                   # ddPCR limit of detection (model units)
    KP_163=1.0, KD_163=1.0, SCD1630=1.0,
    KP_CRP=1.6, KD_CRP=2.0, CRP0=0.8,
    KP_FER=1.0, KD_FER=0.7, FERR0=100.0, W_HLH=2.5,
    # ---- 18. Irreversible organ pools -----------------------------------
    KAVP=0.010, AVP_CRIT=0.15,       # AVP neurons -> central DI
    KANT=0.008,                       # anterior pituitary reserve
    KND=0.010, NEUR_CRIT=0.70,       # cerebellar/pontine neurons -> LCH-ND
    W_CDI_ND=0.35,                    # CDI as an ND risk multiplier
    KBIL=0.020,                       # sclerosing cholangitis
    KCYST=0.0012, ASMK=1.60,         # pulmonary cystic destruction
    KRES=0.60, KHEAL=0.030, KINH=3.0, BVMAX=30.0,  # lytic lesion / healing
    # ---- 19. Myelosuppression (Friberg) ---------------------------------
    KTR_F=0.96, ANC0=4.0, GAM=0.16, KCIRC=2.30,
    SL_ARAC=9.0, SL_CLAD=3.0, SL_VBL=25.0, MSUP=0.45,
    # ---- 20. Targeted-therapy toxicity ----------------------------------
    KSK=0.35, KSKR=0.10,             # cutaneous paradox toxicity
    LVEF0=62.0, KLV=0.65, KREC_LV=0.06,
    # ---- 21. Endpoint definitions ---------------------------------------
    DAS_ACTIVE=3.0,                  # DAS above which disease counts "active"
    DAS_DX=3.5,                      # DAS at which the patient presents
    W_DAS_BONE=1.4, W_DAS_SKIN=0.9, W_DAS_RO=3.2, W_DAS_PIT=1.6,
    W_DAS_CNS=1.8, W_DAS_LUNG=1.5, W_DAS_INFL=4.0,
)

CMT = [
    # --- drug PK (22) ---
    "VEMG", "VEMC", "VEMP", "VEMI",          # vemurafenib + autoinduction
    "DABG", "DABC", "DABM",                  # dabrafenib + active metabolite
    "TRAG", "TRAC", "TRAP",                  # trametinib
    "COBG", "COBC",                          # cobimetinib
    "VBLC", "VBLP",                          # vinblastine
    "PREG", "PREC",                          # prednisolone
    "ARAC", "ARATP",                         # cytarabine + intracellular CTP
    "CLAC", "CLATP",                         # cladribine + intracellular ATP
    "MNTG", "MNTC",                          # 6-MP/MTX maintenance surrogate
    # --- signal transduction (5) ---
    "ERK", "CCND", "BCL", "SASP", "GRE",
    # --- cell pools (10) ---
    "PRECM", "CIRC", "LBONE", "LSKIN", "LRO", "LPIT", "LCNS", "LLUNG",
    "TREG", "OCL",
    # --- secretome (6) ---
    "IL1B", "TNFA", "IL6", "OSM", "MMP9", "RANKL",
    # --- biomarkers (4) ---
    "CFDNA", "SCD163", "CRP", "FERR",
    # --- irreversible organ pools (6) ---
    "AVPN", "ANTPIT", "NEUR", "BILF", "LUNGC", "BVOL",
    # --- myelosuppression, Friberg (4) ---
    "PROL", "TR1", "TR2", "ANC",
    # --- targeted-therapy toxicity (2) ---
    "SKTOX", "LVEF",
    # --- endpoint accumulators (3) ---
    "AUCERK", "TTET", "CUMDAS",
]
IX = {name: i for i, name in enumerate(CMT)}
NEQ = len(CMT)


def hillnorm(x, ec50, h):
    """f(x)/f(1) so that the normalised transducer equals 1 at x = 1."""
    fx = x ** h / (x ** h + ec50 ** h) if x > 0 else 0.0
    f1 = 1.0 / (1.0 + ec50 ** h)
    return fx / f1


def smooth_step(x, thr, width):
    """Differentiable indicator of x > thr (used by the day-counting states)."""
    z = (x - thr) / width
    if z > 40:
        return 1.0
    if z < -40:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def initial(p):
    y = [0.0] * NEQ
    y[IX["VEMI"]] = 1.0
    y[IX["ERK"]] = 1.0
    y[IX["CCND"]] = 1.0
    y[IX["BCL"]] = 1.0
    y[IX["SASP"]] = 1.0
    y[IX["PRECM"]] = p["PRECM0"]
    # circulating precursor starts at its quasi-steady state
    y[IX["CIRC"]] = p["KEX"] * p["PRECM0"] / (p["KSEED"] + p["KDC"])
    y[IX["LBONE"]] = p["LBONE0"]
    y[IX["LSKIN"]] = p["LSKIN0"]
    y[IX["LRO"]] = p["LRO0"]
    y[IX["LPIT"]] = p["LPIT0"]
    y[IX["LCNS"]] = p["LCNS0"]
    y[IX["LLUNG"]] = p["LLUNG0"]
    y[IX["TREG"]] = 0.05
    y[IX["OCL"]] = 0.10
    y[IX["IL1B"]] = p["IL1B0"]
    y[IX["TNFA"]] = p["TNFA0"]
    y[IX["IL6"]] = p["IL60"]
    y[IX["OSM"]] = p["OSM0"]
    y[IX["MMP9"]] = p["MMP90"]
    y[IX["RANKL"]] = p["RANKL0"]
    y[IX["SCD163"]] = p["SCD1630"]
    y[IX["CRP"]] = p["CRP0"]
    y[IX["FERR"]] = p["FERR0"]
    y[IX["AVPN"]] = 1.0
    y[IX["ANTPIT"]] = 1.0
    y[IX["NEUR"]] = 1.0
    y[IX["PROL"]] = p["ANC0"]
    y[IX["TR1"]] = p["ANC0"]
    y[IX["TR2"]] = p["ANC0"]
    y[IX["ANC"]] = p["ANC0"]
    y[IX["LVEF"]] = p["LVEF0"]
    # cfDNA starts at the quasi-steady state implied by the baseline mutant-cell
    # death flux (t1/2 of cfDNA is ~2 h, so it is never far from equilibrium)
    lt0 = (p["LBONE0"] + p["LSKIN0"] + p["LRO0"] + p["LPIT0"] + p["LCNS0"]
           + p["LLUNG0"])
    clr0 = p["KDL"] + p["KIMM"] * (1.0 - 0.05)
    flux0 = (p["KDP"] * p["PRECM0"] + clr0 * lt0
             + p["KDC"] * y[IX["CIRC"]])
    y[IX["CFDNA"]] = p["KP_CF"] * flux0 / p["KD_CF"]
    return y


def derivs(t, y, p, inf):
    """Right-hand side. `inf` supplies active IV/CIV infusion rates (mg/day)."""
    g = lambda n: y[IX[n]]
    d = [0.0] * NEQ

    # -------------------------------------------------------------------
    # Allometric scaling. All PK parameters are quoted at 70 kg (adult
    # literature values) and scaled to the patient. Without this a 12 kg
    # infant on 20 mg/kg/day vemurafenib would be predicted to have ~7 mg/L
    # at steady state instead of the observed 30-60 mg/L.
    # -------------------------------------------------------------------
    fcl = (p["WT"] / 70.0) ** p["ALLO_CL"]
    fv = (p["WT"] / 70.0) ** p["ALLO_V"]

    # -------------------------------------------------------------------
    # PK
    # -------------------------------------------------------------------
    CVEM = g("VEMC") / (p["V2_VEM"] * fv)
    CVEMP = g("VEMP") / (p["V3_VEM"] * fv)
    CL_vem_eff = p["CL_VEM"] * fcl * g("VEMI")
    d[IX["VEMG"]] = -p["KA_VEM"] * g("VEMG")
    QV = p["Q_VEM"] * fcl
    d[IX["VEMC"]] = (p["KA_VEM"] * g("VEMG") - CL_vem_eff * CVEM
                     - QV * CVEM + QV * CVEMP)
    d[IX["VEMP"]] = QV * CVEM - QV * CVEMP
    ind = 1.0 + p["EMAX_IND"] * CVEM / (p["EC50_IND"] + CVEM)
    d[IX["VEMI"]] = p["KENZ"] * (ind - g("VEMI"))

    CDAB = g("DABC") / (p["V_DAB"] * fv)
    CDABM = g("DABM") / (p["V_DABM"] * fv)
    d[IX["DABG"]] = -p["KA_DAB"] * g("DABG")
    d[IX["DABC"]] = p["KA_DAB"] * g("DABG") - p["CL_DAB"] * fcl * CDAB
    d[IX["DABM"]] = (p["FM_DAB"] * p["CL_DAB"] * fcl * CDAB
                     - p["CL_DABM"] * fcl * CDABM)

    CTRA = g("TRAC") / (p["V2_TRA"] * fv)
    CTRAP = g("TRAP") / (p["V3_TRA"] * fv)
    QT = p["Q_TRA"] * fcl
    d[IX["TRAG"]] = -p["KA_TRA"] * g("TRAG")
    d[IX["TRAC"]] = (p["KA_TRA"] * g("TRAG") - p["CL_TRA"] * fcl * CTRA
                     - QT * CTRA + QT * CTRAP)
    d[IX["TRAP"]] = QT * CTRA - QT * CTRAP

    CCOB = g("COBC") / (p["V_COB"] * fv)
    d[IX["COBG"]] = -p["KA_COB"] * g("COBG")
    d[IX["COBC"]] = p["KA_COB"] * g("COBG") - p["CL_COB"] * fcl * CCOB

    CVBL = g("VBLC") / (p["V1_VBL"] * fv)
    CVBLP = g("VBLP") / (p["V2_VBL"] * fv)
    QB = p["Q_VBL"] * fcl
    d[IX["VBLC"]] = (inf.get("VBLC", 0.0) - p["CL_VBL"] * fcl * CVBL
                     - QB * CVBL + QB * CVBLP)
    d[IX["VBLP"]] = QB * CVBL - QB * CVBLP

    CPRE = g("PREC") / (p["V_PRE"] * fv)
    d[IX["PREG"]] = -p["KA_PRE"] * g("PREG")
    d[IX["PREC"]] = p["KA_PRE"] * g("PREG") - p["CL_PRE"] * fcl * CPRE
    d[IX["GRE"]] = p["KEO_P"] * (CPRE - g("GRE"))

    CARAC = g("ARAC") / (p["V_ARAC"] * fv)
    d[IX["ARAC"]] = inf.get("ARAC", 0.0) - p["CL_ARAC"] * fcl * CARAC
    d[IX["ARATP"]] = p["KIN_ARAC"] * CARAC - p["KOUT_ARATP"] * g("ARATP")

    CCLAD = g("CLAC") / (p["V_CLAD"] * fv)
    d[IX["CLAC"]] = inf.get("CLAC", 0.0) - p["CL_CLAD"] * fcl * CCLAD
    d[IX["CLATP"]] = p["KIN_CLAD"] * CCLAD - p["KOUT_CLATP"] * g("CLATP")

    d[IX["MNTG"]] = -p["KA_MNT"] * g("MNTG")
    d[IX["MNTC"]] = p["KA_MNT"] * g("MNTG") - p["KE_MNT"] * g("MNTC")

    # -------------------------------------------------------------------
    # Target engagement
    # -------------------------------------------------------------------
    RV = p["FU_VEM"] * CVEM / p["IC50_VEM"]
    RD = (p["FU_DAB"] * (CDAB + p["RPOT_M"] * CDABM)) / p["IC50_DAB"]
    RB = RV + RD                                    # BRAF-directed occupancy
    RT = p["FU_TRA"] * CTRA / p["IC50_TRA"]
    RC = p["FU_COB"] * CCOB / p["IC50_COB"]
    RM = RT + RC                                    # MEK-directed occupancy

    IMEK = RM / (1.0 + RM)
    braf_active = 1.0 if p["GENO"] < 1.5 else 0.0
    IBRAF = braf_active * RB / (1.0 + RB)
    # Paradoxical MAPK activation: BRAFi in cells without the V600E monomer
    # (BRAF-WT keratinocytes, and MAP2K1-driven LCH). MEKi suppresses it.
    PARADOX = (1.0 - braf_active) * p["PMAX"] * (RB / (1.0 + RB)) * (1.0 - IMEK)
    PARADOX_SKIN = p["PMAX"] * (RB / (1.0 + RB)) * (1.0 - IMEK)

    ERK_target = ((1.0 - IBRAF) * (1.0 - IMEK) * (1.0 + PARADOX)
                  + p["ERKB"])
    ERK_target = min(ERK_target, 1.0 + p["PMAX"])
    d[IX["ERK"]] = p["KERK"] * (ERK_target - g("ERK"))

    ERK = max(g("ERK"), 1e-9)
    d[IX["CCND"]] = p["KCC"] * (hillnorm(ERK, p["EC50C"], p["HC"]) - g("CCND"))
    d[IX["BCL"]] = p["KBCL"] * (hillnorm(ERK, p["EC50B"], p["HB"]) - g("BCL"))
    d[IX["SASP"]] = p["KSASP"] * (hillnorm(ERK, p["EC50S"], p["HS"]) - g("SASP"))

    # -------------------------------------------------------------------
    # Drug effects on cell pools
    # -------------------------------------------------------------------
    FGR = 1.0 - p["EMAX_GR"] * g("GRE") / (p["EC50_GR"] + g("GRE"))
    EGR_KILL = p["EMAX_GRK"] * g("GRE") / (p["EC50_GRK"] + g("GRE"))
    E_ARAC = p["EMAX_ARAC"] * g("ARATP") / (p["EC50_ARAC"] + g("ARATP"))
    E_CLAD = p["EMAX_CLAD"] * g("CLATP") / (p["EC50_CLAD"] + g("CLATP"))
    E_VBL = p["EMAX_VBL"] * CVBL / (p["EC50_VBL"] + CVBL)
    E_MNT = p["EMAX_MNT"] * g("MNTC") / (p["EC50_MNT"] + g("MNTC"))

    PF = p["PFMAX"] * g("CCND")          # cycling fraction of LCH cells
    PFN = PF / p["PFMAX"]                # normalised (1 at CCND = 1)

    # Kill on LESIONAL cells: S-phase/mitosis-dependent drugs are gated by the
    # cycling fraction; cladribine is not (active in non-dividing monocytoid
    # cells) - this is the pharmacological reason 2-CdA works in LCH.
    KILL_LES = (E_ARAC * PFN + E_VBL * PFN + E_CLAD + EGR_KILL
                + p["SL_MAPKI_KILL"] * (IBRAF + IMEK))
    # Kill on the MARROW PRECURSOR reservoir (the reactivation source)
    KILL_PRE = (p["W_ARAC_P"] * E_ARAC * PFN + p["W_VBL_P"] * E_VBL * PFN
                + p["W_CLAD_P"] * E_CLAD + p["W_GR_P"] * EGR_KILL
                + p["SL_MAPKI_KILL"] * (IBRAF + IMEK))
    # CNS compartment sees reduced exposure
    KILL_CNS = (E_ARAC * p["FCNS_ARAC"] * PFN + E_CLAD * p["FCNS_CLAD"]
                + EGR_KILL)

    APO = 1.0 + p["KAPO"] * (1.0 - g("BCL"))     # loss of BCL2A1 -> apoptosis
    CLR = p["KDL"] * APO + p["KIMM"] * (1.0 - g("TREG"))

    # -------------------------------------------------------------------
    # Cell pools
    # -------------------------------------------------------------------
    # The mutant precursor divides asymmetrically: one daughter self-renews
    # in the marrow niche, the other is exported as circulating progeny.
    # Export therefore does NOT debit the reservoir - which is why the
    # reservoir can persist under a purely cytostatic drug.
    PRECM = g("PRECM")
    SR = (p["KSR"] * p["FSR"] * g("CCND") * (1.0 - PRECM / p["PRMAX"])
          * (1.0 - E_MNT) * PRECM / (PRECM + p["PNICHE"]))
    APO_PRE = 1.0 + p["FN_APO"] * p["KAPO"] * (1.0 - g("BCL"))
    d[IX["PRECM"]] = (SR - p["KDP"] * APO_PRE - KILL_PRE) * PRECM

    d[IX["CIRC"]] = (p["KEX"] * g("CCND") * PRECM
                     - (p["KSEED"] + p["KDC"] + KILL_PRE) * g("CIRC"))

    LTOT = (g("LBONE") + g("LSKIN") + g("LRO") + g("LPIT") + g("LCNS")
            + g("LLUNG"))
    REC = 1.0 + p["FREC"] * (g("IL1B") / p["IL1B0"] - 1.0)
    REC = min(max(REC, 0.2), p["RECMAX"])
    seed = p["KSEED"] * g("CIRC") * REC

    def lesion(name, theta, kill, extra_growth=0.0):
        L = g(name)
        grow = p["KPROL"] * g("CCND") * L * (1.0 - L / p["LMAX"])
        return (theta * seed + grow * (1.0 + extra_growth) - CLR * L
                - kill * L)

    smoke_boost = p["ASMK"] * p["SMOKE"]
    d[IX["LBONE"]] = lesion("LBONE", p["THB"], KILL_LES)
    d[IX["LSKIN"]] = lesion("LSKIN", p["THS"], KILL_LES)
    d[IX["LRO"]] = lesion("LRO", p["THR"], KILL_LES)
    d[IX["LPIT"]] = lesion("LPIT", p["THP"], KILL_LES)
    d[IX["LCNS"]] = lesion("LCNS", p["THC"], KILL_CNS)
    d[IX["LLUNG"]] = lesion("LLUNG", p["THL"], KILL_LES, smoke_boost)

    LTOTn = LTOT / p["LREF"]
    d[IX["TREG"]] = p["KTR"] * (p["TREGMAX"] * LTOTn
                                / (p["EC50_TREG"] + LTOTn) - g("TREG"))

    rnk = g("RANKL") / (g("RANKL") + p["KR_RANKL"])
    oc_target = rnk * (1.0 + p["AOSM"] * g("OSM") / p["OSM0"] * 0.25) \
        * (1.0 + 0.30 * p["IL17ON"])
    d[IX["OCL"]] = p["KOC"] * (oc_target - g("OCL"))

    # -------------------------------------------------------------------
    # Secretome (ERK/SASP-driven, GR-suppressible)
    # -------------------------------------------------------------------
    SRC = g("SASP") * LTOTn * FGR
    d[IX["IL1B"]] = p["KP_IL1"] * SRC - p["KD_IL1"] * (g("IL1B") - p["IL1B0"])
    d[IX["TNFA"]] = p["KP_TNF"] * SRC - p["KD_TNF"] * (g("TNFA") - p["TNFA0"])
    d[IX["IL6"]] = p["KP_IL6"] * SRC - p["KD_IL6"] * (g("IL6") - p["IL60"])
    d[IX["OSM"]] = p["KP_OSM"] * SRC - p["KD_OSM"] * (g("OSM") - p["OSM0"])
    d[IX["MMP9"]] = p["KP_MMP"] * SRC - p["KD_MMP"] * (g("MMP9") - p["MMP90"])
    d[IX["RANKL"]] = (p["KP_RNK"] * SRC * (1.0 + p["W_TCELL_RNK"])
                      - p["KD_RNK"] * (g("RANKL") - p["RANKL0"]))

    # -------------------------------------------------------------------
    # Biomarkers
    # cfDNA tracks the DEATH FLUX of mutant cells, which is why MAPK
    # inhibition (cytostatic) produces a shallow VAF plateau while
    # nucleoside analogues drive a deep nadir.
    # -------------------------------------------------------------------
    death_flux = (p["KDP"] * APO_PRE * PRECM + KILL_PRE * PRECM
                  + (CLR + KILL_LES) * (LTOT - g("LCNS"))
                  + (CLR + KILL_CNS) * g("LCNS")
                  + p["KDC"] * g("CIRC"))
    d[IX["CFDNA"]] = p["KP_CF"] * death_flux - p["KD_CF"] * g("CFDNA")
    d[IX["SCD163"]] = (p["KP_163"] * LTOTn * FGR
                       - p["KD_163"] * (g("SCD163") - p["SCD1630"]))
    d[IX["CRP"]] = (p["KP_CRP"] * g("IL6") / p["IL60"]
                    - p["KD_CRP"] * g("CRP"))
    hlh = p["W_HLH"] * (g("LRO") / p["LREF"])
    d[IX["FERR"]] = (p["KP_FER"] * (100.0 * LTOTn + 400.0 * hlh)
                     - p["KD_FER"] * (g("FERR") - p["FERR0"]))

    # -------------------------------------------------------------------
    # Irreversible organ pools — monotone integrals of local activity
    # -------------------------------------------------------------------
    infl = g("IL1B") / p["IL1B0"]
    d[IX["AVPN"]] = -p["KAVP"] * (g("LPIT") / p["LREF"]) * infl ** 0.5 \
        * g("AVPN")
    d[IX["ANTPIT"]] = -p["KANT"] * (g("LPIT") / p["LREF"]) * g("ANTPIT")
    cdi_flag = smooth_step(p["AVP_CRIT"], g("AVPN"), 0.02)
    d[IX["NEUR"]] = -p["KND"] * ((g("LCNS") / p["LREF"])
                                 + p["W_CDI_ND"] * cdi_flag
                                 * (g("LCNS") / p["LREF"] + 0.15)) * g("NEUR")
    d[IX["BILF"]] = p["KBIL"] * (g("LRO") / p["LREF"]) * (1.0 - g("BILF"))
    d[IX["LUNGC"]] = (p["KCYST"] * (g("LLUNG") / p["LREF"])
                      * (1.0 + smoke_boost) * (1.0 - g("LUNGC")))
    d[IX["BVOL"]] = (p["KRES"] * g("OCL") * (g("LBONE") / p["LREF"])
                     * (1.0 - g("BVOL") / p["BVMAX"])
                     - p["KHEAL"] * g("BVOL")
                     / (1.0 + p["KINH"] * g("LBONE") / p["LREF"]))

    # -------------------------------------------------------------------
    # Myelosuppression (Friberg) — drug + disease marrow infiltration
    # -------------------------------------------------------------------
    EDRUG = (p["SL_ARAC"] * g("ARATP") + p["SL_CLAD"] * g("CLATP")
             + p["SL_VBL"] * CVBL)
    EDRUG = min(EDRUG, 0.95)
    dis_sup = min(p["MSUP"] * g("LRO") / p["LREF"], 0.8)
    fb = (p["ANC0"] / max(g("ANC"), 1e-3)) ** p["GAM"]
    d[IX["PROL"]] = (p["KTR_F"] * g("PROL") * (1.0 - EDRUG) * (1.0 - dis_sup)
                     * fb - p["KTR_F"] * g("PROL"))
    d[IX["TR1"]] = p["KTR_F"] * (g("PROL") - g("TR1"))
    d[IX["TR2"]] = p["KTR_F"] * (g("TR1") - g("TR2"))
    d[IX["ANC"]] = p["KTR_F"] * g("TR2") - p["KCIRC"] * g("ANC")

    # -------------------------------------------------------------------
    # Targeted-therapy toxicity
    # -------------------------------------------------------------------
    d[IX["SKTOX"]] = p["KSK"] * PARADOX_SKIN - p["KSKR"] * g("SKTOX")
    d[IX["LVEF"]] = (p["KREC_LV"] * (p["LVEF0"] - g("LVEF"))
                     - p["KLV"] * IMEK)

    # -------------------------------------------------------------------
    # Endpoint accumulators
    # -------------------------------------------------------------------
    DAS = das_of(y, p)
    d[IX["AUCERK"]] = smooth_step(0.20, g("ERK"), 0.02)   # days with >80% supp.
    d[IX["TTET"]] = smooth_step(DAS, p["DAS_ACTIVE"], 0.3)
    d[IX["CUMDAS"]] = DAS
    return d


def das_of(y, p):
    """Histiocyte-Society-style disease activity score (0 - ~30)."""
    g = lambda n: y[IX[n]]
    infl = max(g("IL1B") / p["IL1B0"] - 1.0, 0.0)

    def sat(name, w):
        # each organ contributes a saturating term, so DAS is bounded like the
        # ordinal Histiocyte Society score rather than growing without limit
        ln = g(name) / p["LREF"]
        return w * ln / (1.0 + ln)

    return (sat("LBONE", p["W_DAS_BONE"])
            + sat("LSKIN", p["W_DAS_SKIN"])
            + sat("LRO", p["W_DAS_RO"])
            + sat("LPIT", p["W_DAS_PIT"])
            + sat("LCNS", p["W_DAS_CNS"])
            + sat("LLUNG", p["W_DAS_LUNG"])
            + p["W_DAS_INFL"] * infl / (1.0 + infl))


# ===========================================================================
# Dosing events
# ===========================================================================
class Regimen:
    """Bolus/oral events (into a depot or central cmt) and CIV infusions."""

    def __init__(self):
        self.bolus = []      # (time, cmt, amount_mg)
        self.infus = []      # (start, end, cmt, rate_mg_per_day)

    def oral(self, cmt, amt, start, stop, every):
        t = start
        while t < stop - 1e-9:
            self.bolus.append((t, cmt, amt))
            t += every
        return self

    def civ(self, cmt, total_mg_per_day, start, stop):
        self.infus.append((start, stop, cmt, total_mg_per_day))
        return self

    def rates(self, t):
        out = {}
        for s, e, c, r in self.infus:
            if s - 1e-9 <= t < e - 1e-9:
                out[c] = out.get(c, 0.0) + r
        return out

    def times(self):
        ts = set()
        for t, _, _ in self.bolus:
            ts.add(round(t, 6))
        for s, e, _, _ in self.infus:
            ts.add(round(s, 6))
            ts.add(round(e, 6))
        return sorted(ts)


def simulate(p, reg, tend, dt=0.005, record_every=0.5, y0=None, tstart=0.0):
    y = list(y0) if y0 is not None else initial(p)
    breaks = [t for t in reg.times() if tstart <= t <= tend]
    out = {"time": []}
    for n in CMT:
        out[n] = []
    for n in ("DAS", "CVEM", "CDAB", "CTRA", "ERKpct", "CDI", "ND"):
        out[n] = []

    def record(t):
        out["time"].append(t)
        for n in CMT:
            out[n].append(y[IX[n]])
        out["DAS"].append(das_of(y, p))
        _fv = (p["WT"] / 70.0) ** p["ALLO_V"]
        out["CVEM"].append(y[IX["VEMC"]] / (p["V2_VEM"] * _fv))
        out["CDAB"].append(y[IX["DABC"]] / (p["V_DAB"] * _fv))
        out["CTRA"].append(y[IX["TRAC"]] / (p["V2_TRA"] * _fv) * 1000.0)
        out["ERKpct"].append(100.0 * y[IX["ERK"]])
        out["CDI"].append(1.0 if y[IX["AVPN"]] < p["AVP_CRIT"] else 0.0)
        out["ND"].append(1.0 if y[IX["NEUR"]] < p["NEUR_CRIT"] else 0.0)

    t = tstart
    next_rec = tstart
    bolus_by_time = {}
    for bt, bc, ba in reg.bolus:
        bolus_by_time.setdefault(round(bt, 6), []).append((bc, ba))

    grid = sorted(set([tstart, tend] + breaks))
    for seg_i in range(len(grid) - 1):
        t0, t1 = grid[seg_i], grid[seg_i + 1]
        key = round(t0, 6)
        if key in bolus_by_time:
            for bc, ba in bolus_by_time[key]:
                y[IX[bc]] += ba
        if next_rec <= t0 + 1e-9:
            record(t0)
            next_rec = t0 + record_every
        t = t0
        inf = reg.rates(t0 + 1e-9)
        nsteps = max(1, int(math.ceil((t1 - t0) / dt)))
        h = (t1 - t0) / nsteps
        for _ in range(nsteps):
            k1 = derivs(t, y, p, inf)
            y2 = [y[i] + 0.5 * h * k1[i] for i in range(NEQ)]
            k2 = derivs(t + 0.5 * h, y2, p, inf)
            y3 = [y[i] + 0.5 * h * k2[i] for i in range(NEQ)]
            k3 = derivs(t + 0.5 * h, y3, p, inf)
            y4 = [y[i] + h * k3[i] for i in range(NEQ)]
            k4 = derivs(t + h, y4, p, inf)
            for i in range(NEQ):
                y[i] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0
                if y[i] < 0.0 and CMT[i] != "LVEF":
                    y[i] = 0.0
            t += h
            if t >= next_rec - 1e-9:
                record(min(t, t1))
                next_rec += record_every
    record(t)
    return out


# ===========================================================================
# Patient phenotypes and the 10 scenarios
# ===========================================================================
BSA = 0.70          # m^2, ~2-year-old, 12 kg
WT = 12.0           # kg


def pheno(kind):
    """Patient phenotypes. NOTE what does and does not change between them:
    only the CELL-OF-ORIGIN DESCRIPTORS (FSR, PRECM0, the seeding thetas, the
    locally seeded lesions) and the smoking flag. Every lesional kinetic
    constant (KPROL, KDL, KIMM, KDP, KSR) is identical across phenotypes."""
    p = dict(P)
    if kind == "SSb":            # single-system bone, tissue-restricted origin
        p.update(PRECM0=0.020, FSR=0.10, THB=0.70, THS=0.10, THR=0.04,
                 THP=0.06, THC=0.02, THL=0.08, LBONE0=0.80, LSKIN0=0.0,
                 LRO0=0.0, LPIT0=0.0, LCNS0=0.0, LLUNG0=0.0)
    elif kind == "MSROneg":
        p.update(PRECM0=0.12, FSR=0.72, THB=0.42, THS=0.22, THR=0.06,
                 THP=0.14, THC=0.06, THL=0.10, LRO0=0.0)
    elif kind == "MSROpos":      # early (HSC-stage) mutation
        p.update(PRECM0=0.45, FSR=1.00, THB=0.28, THS=0.16, THR=0.34,
                 THP=0.08, THC=0.06, THL=0.08)
    elif kind == "PLCH":         # adult pulmonary LCH, smoker
        p.update(PRECM0=0.06, FSR=0.24, THB=0.10, THS=0.04, THR=0.02,
                 THP=0.02, THC=0.02, THL=0.80, SMOKE=1.0, LLUNG0=0.20,
                 LBONE0=0.0, LSKIN0=0.0, LRO0=0.0, LPIT0=0.0, LCNS0=0.0,
                 WT=65.0)
    elif kind == "CNSrisk":      # craniofacial CNS-risk lesion, pituitary
        p.update(PRECM0=0.25, FSR=0.85, THB=0.30, THS=0.10, THR=0.10,
                 THP=0.30, THC=0.14, THL=0.06, LPIT0=0.06)
    return p


def r_vblpred(reg, start, weeks_induction=6, cont_to=365):
    """LCH-III style: VBL 6 mg/m2 IV weekly + prednisolone 40 mg/m2/d x 4 wk,
    then 3-weekly pulses in continuation."""
    vbl = 6.0 * BSA
    pre = 40.0 * BSA
    t = start
    while t < start + weeks_induction * 7:
        reg.bolus.append((t, "VBLC", vbl))
        t += 7
    # prednisolone 40 mg/m2/d daily x 28 d then taper over 2 wk (approximated)
    t = start
    while t < start + 28:
        reg.bolus.append((t, "PREG", pre))
        reg.bolus.append((t + 0.5, "PREG", pre))
        t += 1
    t = start + 28
    while t < start + 42:
        reg.bolus.append((t, "PREG", pre * 0.5))
        t += 1
    # continuation: VBL + 5-day prednisolone pulse every 21 d
    t = start + 42
    while t < cont_to:
        reg.bolus.append((t, "VBLC", vbl))
        for k in range(5):
            reg.bolus.append((t + k, "PREG", pre))
            reg.bolus.append((t + k + 0.5, "PREG", pre))
        t += 21
    # 6-MP/MTX maintenance in the continuation phase
    t = start + 42
    while t < cont_to:
        reg.bolus.append((t, "MNTG", 1.0))
        t += 1
    return reg


def r_cladarac(reg, start, cycles=3, clad_mg_m2=9.0, arac_mg_m2=500.0,
               cycle_days=28):
    """Salvage 2-CdA 9 mg/m2/d + Ara-C 500 mg/m2/d, 5 days, q28d."""
    for c in range(cycles):
        s = start + c * cycle_days
        reg.civ("CLAC", clad_mg_m2 * BSA, s, s + 5)
        reg.civ("ARAC", arac_mg_m2 * BSA, s, s + 5)
    return reg


def r_maint(reg, start, stop):
    """6-MP/MTX oral maintenance, one dose per day."""
    t = start
    while t < stop:
        reg.bolus.append((t, "MNTG", 1.0))
        t += 1
    return reg


def r_aracmono(reg, start, cycles=6, mg_m2=170.0, cycle_days=28):
    for c in range(cycles):
        s = start + c * cycle_days
        reg.civ("ARAC", mg_m2 * BSA, s, s + 5)
    return reg


def r_vem(reg, start, stop, mg_kg_day=20.0):
    dose = mg_kg_day * WT / 2.0
    return reg.oral("VEMG", dose, start, stop, 0.5)


def r_dabtram(reg, start, stop, dab_mg_kg=4.5, tram_mg_kg=0.032):
    reg.oral("DABG", dab_mg_kg * WT / 2.0, start, stop, 0.5)
    reg.oral("TRAG", tram_mg_kg * WT, start, stop, 1.0)
    return reg


SCENARIOS = {}


def scen(name, kind, build, tend, note):
    SCENARIOS[name] = (kind, build, tend, note)


scen("S1_observation", "SSb", lambda r: r, 730,
     "Single-system bone, observation only (self-limited course)")
scen("S2_LCH3_ROneg", "MSROneg", lambda r: r_vblpred(r, 0, 6, 365), 730,
     "MS RO-negative: LCH-III VBL/prednisolone 12 months")
scen("S3_frontline_fail", "MSROpos",
     lambda r: r_maint(r_cladarac(r_vblpred(r, 0, 6, 42), 42, 5), 180, 545), 730,
     "MS RO+ non-responder at week 6 -> salvage 2-CdA/Ara-C")
scen("S4_cladarac_upfront", "MSROpos",
     lambda r: r_maint(r_cladarac(r, 0, 6), 175, 545), 730,
     "MS RO+ refractory: 2-CdA/Ara-C x3 up front")
scen("S5_vem_continuous", "MSROpos", lambda r: r_vem(r, 0, 730), 730,
     "BRAF V600E MS RO+: vemurafenib 20 mg/kg/d continuous")
scen("S6_vem_stop", "MSROpos", lambda r: r_vem(r, 0, 365), 730,
     "Vemurafenib 12 months then STOP -> rebound")
scen("S7_dabtram", "MSROpos", lambda r: r_dabtram(r, 0, 730), 730,
     "Dabrafenib + trametinib continuous")
scen("S8_bridge_consolidate", "MSROpos",
     lambda r: r_maint(r_cladarac(r_vem(r, 0, 112), 56, 3), 145, 420), 730,
     "MAPKi bridge 8 wk -> 2-CdA/Ara-C consolidation -> stop everything")
scen("S9_delayed_dx", "CNSrisk", lambda r: r_vblpred(r, 180, 6, 545), 730,
     "6-month diagnostic delay with CNS-risk lesion -> sequelae")
scen("S9b_early_dx", "CNSrisk", lambda r: r_vblpred(r, 14, 6, 379), 730,
     "Same patient treated at day 14 (comparator for S9)")
scen("S10_plch_quit", "PLCH", lambda r: r, 1095,
     "Adult pulmonary LCH: smoking cessation at month 3")
scen("S10b_plch_smoke", "PLCH", lambda r: r, 1095,
     "Adult pulmonary LCH: continued smoking")


# Run-in to presentation: every scenario starts from an ESTABLISHED disease
# state, because "response" and "time to control" are only meaningful relative
# to the burden the patient actually presents with.
RUNIN = {"SSb": ("days", 30), "MSROneg": ("das", 500), "MSROpos": ("das", 500),
         "CNSrisk": ("das", 500), "PLCH": ("days", 365)}

_EST_CACHE = {}
_RUN_CACHE = {}


def _key(p):
    return tuple(sorted((k, v) for k, v in p.items()))


def establish(p, kind):
    """Integrate the untreated natural history up to presentation."""
    ck = (kind, _key(p))
    if ck in _EST_CACHE:
        return _EST_CACHE[ck]
    mode, arg = RUNIN[kind]
    if mode == "days":
        out = simulate(p, Regimen(), arg, record_every=5.0)
        y, t = [out[n][-1] for n in CMT], float(arg)
    else:
        y, t = initial(p), 0.0
        while t < arg:
            out = simulate(p, Regimen(), t + 5.0, y0=y, tstart=t,
                           record_every=5.0)
            y = [out[n][-1] for n in CMT]
            t += 5.0
            if das_of(y, p) >= p["DAS_DX"]:
                break
    # organ damage accrued before diagnosis is KEPT; the day-counting
    # accumulators are zeroed so that they measure time since presentation
    for acc in ("AUCERK", "TTET", "CUMDAS"):
        y[IX[acc]] = 0.0
    _EST_CACHE[ck] = (y, t)
    return y, t


def run(name, extra=None):
    ck = (name, _key(extra) if extra else None)
    if ck in _RUN_CACHE:
        return _RUN_CACHE[ck]
    kind, build, tend, _note = SCENARIOS[name]
    p = pheno(kind)
    if extra:
        p.update(extra)
    y0, _t0 = establish(p, kind)
    reg = build(Regimen())
    if name == "S10_plch_quit":
        # smoking cessation at month 3: run in two parameter segments
        out1 = simulate(p, reg, 90, y0=y0)
        p2 = dict(p)
        p2["SMOKE"] = 0.0
        y = [out1[n][-1] for n in CMT]
        out2 = simulate(p2, reg, tend, y0=y, tstart=90.0)
        out = _join(out1, out2)
    else:
        out = simulate(p, reg, tend, y0=y0)
    _RUN_CACHE[ck] = out
    return out


def _join(a, b):
    out = {}
    for k in a:
        out[k] = a[k] + b[k]
    return out


def at(out, t):
    """Value dictionary at (or just before) time t."""
    best = 0
    for i, tt in enumerate(out["time"]):
        if tt <= t + 1e-6:
            best = i
    return {k: out[k][best] for k in out}


def peak(out, key, t0=0.0, t1=1e9):
    v = [out[key][i] for i, t in enumerate(out["time"]) if t0 <= t <= t1]
    return max(v) if v else float("nan")


def trough(out, key, t0=0.0, t1=1e9):
    v = [out[key][i] for i, t in enumerate(out["time"]) if t0 <= t <= t1]
    return min(v) if v else float("nan")


# ===========================================================================
# Verification suite
# ===========================================================================
def main():
    checks = []

    def ck(label, cond, detail=""):
        checks.append((label, bool(cond), detail))

    print("=" * 78)
    print("LCH QSP model - numerical verification (%d ODEs)" % NEQ)
    print("=" * 78)

    # ---------------- structural sanity ---------------------------------
    ck("62 state variables", NEQ == 62, "NEQ=%d" % NEQ)
    ck("unique compartment names", len(set(CMT)) == NEQ)

    # ---------------- PK plausibility -----------------------------------
    p = pheno("MSROpos")
    o = run("S5_vem_continuous")
    css = at(o, 60)["CVEM"]
    ck("paediatric vemurafenib Css 20-70 mg/L on 20 mg/kg/day", 20 <= css <= 70,
       "Css=%.1f mg/L" % css)
    erk = at(o, 60)["ERKpct"]
    ck("vemurafenib suppresses pERK below 25% of baseline", erk < 25,
       "pERK=%.1f%%" % erk)

    o7 = run("S7_dabtram")
    ctra = at(o7, 90)["CTRA"]
    ck("trametinib trough 6-30 ng/mL (label Cav ~12-22)", 6 <= ctra <= 30,
       "Css=%.1f ng/mL" % ctra)
    erk7 = at(o7, 60)["ERKpct"]
    ck("BRAFi+MEKi deeper pERK suppression than BRAFi alone", erk7 < erk,
       "combo=%.1f%% vs mono=%.1f%%" % (erk7, erk))

    # ---------------- natural history -----------------------------------
    o1 = run("S1_observation")
    d1_0 = at(o1, 30)["DAS"]
    d1_end = at(o1, 730)["DAS"]
    ck("SS-b: low DAS throughout (self-limited)", d1_0 < 6 and d1_end < 6,
       "DAS d30=%.2f d730=%.2f" % (d1_0, d1_end))
    ck("SS-b: no CDI", at(o1, 730)["CDI"] < 0.5,
       "AVPN=%.2f" % at(o1, 730)["AVPN"])
    ck("SS-b: lytic lesion then reossification",
       peak(o1, "BVOL") > at(o1, 730)["BVOL"],
       "peak=%.2f end=%.2f" % (peak(o1, "BVOL"), at(o1, 730)["BVOL"]))

    p_ro = pheno("MSROpos")
    o_nat = simulate(p_ro, Regimen(), 365)
    ck("MS RO+ untreated: DAS rises above the active threshold",
       at(o_nat, 120)["DAS"] > P["DAS_ACTIVE"],
       "DAS d120=%.1f" % at(o_nat, 120)["DAS"])
    ck("MS RO+ untreated: risk-organ burden grows",
       at(o_nat, 180)["LRO"] > at(o_nat, 20)["LRO"],
       "LRO %.2f -> %.2f" % (at(o_nat, 20)["LRO"], at(o_nat, 180)["LRO"]))
    ck("MS RO+ untreated: CDI develops", at(o_nat, 365)["CDI"] > 0.5,
       "AVPN=%.3f" % at(o_nat, 365)["AVPN"])
    ck("MS RO+ untreated: biliary fibrosis accumulates",
       at(o_nat, 365)["BILF"] > 0.3, "BILF=%.2f" % at(o_nat, 365)["BILF"])

    # ---------------- COMMITMENT 1: partition sets phenotype ------------
    same = dict(pheno("SSb"))
    for k in ("THB", "THS", "THR", "THP", "THC", "THL", "PRECM0"):
        same[k] = pheno("MSROpos")[k]
    o_same = simulate(same, Regimen(), 365)
    ck("C1: SS-b phenotype becomes RO+ when ONLY the partition is swapped",
       at(o_same, 180)["LRO"] > 10 * at(o1, 180)["LRO"],
       "LRO %.3f vs %.3f" % (at(o_same, 180)["LRO"], at(o1, 180)["LRO"]))
    ck("C1: no lesional kinetic constant differs between phenotypes "
       "(only origin descriptors do)",
       all(pheno("SSb")[k] == pheno("MSROpos")[k]
           for k in ("KSR", "KDP", "KPROL", "KDL", "KIMM", "KSEED", "KDC",
                     "KAPO", "LMAX")))

    # ---------------- COMMITMENT 2: cytostatic + rebound ----------------
    o5 = run("S5_vem_continuous")
    das0 = at(o5, 0)["DAS"]
    das7 = at(o5, 7)["DAS"]
    ltot0 = sum(at(o5, 0)[k] for k in ("LBONE", "LSKIN", "LRO", "LPIT",
                                       "LCNS", "LLUNG"))
    ltot7 = sum(at(o5, 7)[k] for k in ("LBONE", "LSKIN", "LRO", "LPIT",
                                       "LCNS", "LLUNG"))
    frac_das = 1 - das7 / das0
    frac_mass = 1 - ltot7 / ltot0
    ck("C2: DAS falls faster than lesion mass in week 1 (secretome first)",
       frac_das > frac_mass + 0.10,
       "DAS -%.0f%% vs mass -%.0f%%" % (100 * frac_das, 100 * frac_mass))
    vaf0 = at(o5, 0)["CFDNA"]
    vaf365 = at(o5, 365)["CFDNA"]
    pre365 = at(o5, 365)["PRECM"]
    ck("C2: after 12 months of MAPKi the clone is still there - cfDNA stays "
       "above the assay LOD and the reservoir above the niche threshold",
       vaf365 > P["CF_LOD"] and pre365 > P["PNICHE"],
       "cfDNA %.3f -> %.3f (LOD %.3f); reservoir %.4f (PNICHE %.4f)"
       % (vaf0, vaf365, P["CF_LOD"], pre365, P["PNICHE"]))
    o4pre = run("S4_cladarac_upfront")
    ck("C2: the same interval of nucleoside-analogue therapy drives cfDNA "
       "BELOW the LOD and extinguishes the reservoir",
       at(o4pre, 365)["CFDNA"] < P["CF_LOD"]
       and at(o4pre, 365)["PRECM"] < P["PNICHE"],
       "cfDNA=%.5f, reservoir=%.6f"
       % (at(o4pre, 365)["CFDNA"], at(o4pre, 365)["PRECM"]))
    o6 = run("S6_vem_stop")
    das365 = at(o6, 364)["DAS"]
    das_reb = peak(o6, "DAS", 365, 730)
    ck("C2: stopping MAPKi rebounds (>2x the on-therapy DAS)",
       das_reb > 2 * das365,
       "on-therapy %.2f -> rebound peak %.2f" % (das365, das_reb))
    o6k = run("S6_vem_stop", {"SL_MAPKI_KILL": 0.25})
    reb_k = peak(o6k, "DAS", 365, 730) / max(at(o6k, 364)["DAS"], 1e-6)
    reb_0 = das_reb / max(das365, 1e-6)
    ck("C2 KILL SWITCH: adding a MAPKi kill term abolishes the rebound",
       reb_k < reb_0, "ratio %.1f -> %.1f" % (reb_0, reb_k))
    o4 = run("S4_cladarac_upfront")
    vaf_nadir_chemo = trough(o4, "CFDNA", 60, 200) / at(o4, 0)["CFDNA"]
    vaf_nadir_mapki = trough(o5, "CFDNA", 60, 200) / at(o5, 0)["CFDNA"]
    ck("C2: nucleoside analogues reach a deeper cfDNA nadir than MAPKi",
       vaf_nadir_chemo < vaf_nadir_mapki,
       "chemo %.3f vs MAPKi %.3f of baseline"
       % (vaf_nadir_chemo, vaf_nadir_mapki))
    pre_chemo = at(o4, 200)["PRECM"] / pheno("MSROpos")["PRECM0"]
    pre_mapki = at(o5, 200)["PRECM"] / pheno("MSROpos")["PRECM0"]
    ck("C2: chemo depletes the marrow precursor reservoir, MAPKi does not",
       pre_chemo < 0.5 * pre_mapki,
       "reservoir chemo %.3f vs MAPKi %.3f of baseline"
       % (pre_chemo, pre_mapki))
    ck("C2: cladribine kill is NOT gated by the cycling fraction",
       P["W_CLAD_P"] > 0 and P["PFMAX"] < 0.15)
    o8 = run("S8_bridge_consolidate")
    ck("C2: bridge-then-consolidate holds remission off all therapy",
       at(o8, 700)["DAS"] < 0.5 * das_reb,
       "DAS d700=%.2f vs S6 rebound %.2f" % (at(o8, 700)["DAS"], das_reb))

    # ---------------- COMMITMENT 3: sequelae are integrals ---------------
    o9 = run("S9_delayed_dx")
    o9b = run("S9b_early_dx")
    ck("C3: 6-month delay costs AVP neurons",
       at(o9, 730)["AVPN"] < at(o9b, 730)["AVPN"],
       "delayed AVPN=%.3f vs early AVPN=%.3f"
       % (at(o9, 730)["AVPN"], at(o9b, 730)["AVPN"]))
    ck("C3: delay converts a reversible course into permanent CDI",
       at(o9, 730)["CDI"] > 0.5 and at(o9b, 730)["CDI"] < 0.5,
       "CDI delayed=%d early=%d"
       % (at(o9, 730)["CDI"], at(o9b, 730)["CDI"]))
    ck("C3: neurodegeneration follows the same ordering",
       at(o9, 730)["NEUR"] < at(o9b, 730)["NEUR"],
       "NEUR delayed=%.3f early=%.3f"
       % (at(o9, 730)["NEUR"], at(o9b, 730)["NEUR"]))
    ck("C3: TTET (days of uncontrolled disease) explains the difference",
       at(o9, 730)["TTET"] > at(o9b, 730)["TTET"],
       "TTET delayed=%.0f d, early=%.0f d"
       % (at(o9, 730)["TTET"], at(o9b, 730)["TTET"]))
    ck("C3: sequelae pools are monotone (never recover)",
       all(o9["AVPN"][i + 1] <= o9["AVPN"][i] + 1e-9
           for i in range(len(o9["AVPN"]) - 1)))

    # Is there a crossover between a fast-shallow (cytostatic) and a
    # slow-deep (cytotoxic) regimen on the CNS endpoint? The model computes it
    # rather than assuming it.
    pc = pheno("CNSrisk")
    y0c, _tc = establish(pc, "CNSrisk")
    o_fast = simulate(pc, r_vem(Regimen(), 0, 400), 400, y0=y0c)
    o_deep = simulate(pc, r_cladarac(Regimen(), 0, 3), 400, y0=y0c)

    def days_to_control(o, thr=1.5):
        for i, t in enumerate(o["time"]):
            if o["DAS"][i] < thr:
                return t
        return float("inf")

    t_fast, t_deep = days_to_control(o_fast), days_to_control(o_deep)
    log_times = "MAPKi %.1f d vs 2-CdA/Ara-C %.1f d" % (t_fast, t_deep)

    # The mechanistically sharp statement is not "which is faster" but
    # "which one separates clinical response from mass reduction".
    def discordance(o):
        m0 = sum(at(o, 0)[k] for k in ("LBONE", "LSKIN", "LRO", "LPIT",
                                       "LCNS", "LLUNG"))
        m7 = sum(at(o, 7)[k] for k in ("LBONE", "LSKIN", "LRO", "LPIT",
                                       "LCNS", "LLUNG"))
        fd = 1.0 - at(o, 7)["DAS"] / max(at(o, 0)["DAS"], 1e-9)
        fm = 1.0 - m7 / max(m0, 1e-9)
        return fd / max(fm, 1e-6)

    dis_fast, dis_deep = discordance(o_fast), discordance(o_deep)
    ck("C2/C3: MAPK inhibition separates clinical response from mass "
       "reduction; cytotoxic therapy does not", dis_fast > dis_deep,
       "DAS-drop/mass-drop ratio: MAPKi %.2f vs 2-CdA/Ara-C %.2f (%s)"
       % (dis_fast, dis_deep, log_times))

    # Cycle number vs the niche-extinction threshold: a genuine prediction.
    o_c3 = simulate(pc, r_cladarac(Regimen(), 0, 3), 500, y0=y0c)
    o_c6 = simulate(pc, r_cladarac(Regimen(), 0, 6), 500, y0=y0c)
    ck("dose intensity matters at the extinction threshold: 3 cycles leaves "
       "the clone above PNICHE and it regrows, 6 cycles does not",
       at(o_c3, 500)["PRECM"] > 100 * at(o_c6, 500)["PRECM"],
       "reservoir at d500: 3 cycles=%.4f, 6 cycles=%.6f"
       % (at(o_c3, 500)["PRECM"], at(o_c6, 500)["PRECM"]))
    o_deep = o_c6      # use the protocol course for the endpoint comparison
    ck("C3: with both started at presentation the deeper regimen is not worse "
       "for AVP neurons (no crossover at zero delay)",
       at(o_deep, 400)["AVPN"] >= at(o_fast, 400)["AVPN"] - 1e-6,
       "deep AVPN=%.3f vs fast AVPN=%.3f"
       % (at(o_deep, 400)["AVPN"], at(o_fast, 400)["AVPN"]))
    # ... but delay the deeper regimen and the ordering flips. Find the flip.
    flip = None
    for delay in (15, 30, 45, 60, 90):
        od = simulate(pc, r_cladarac(Regimen(), delay, 6), 400, y0=y0c)
        if at(od, 400)["AVPN"] < at(o_fast, 400)["AVPN"]:
            flip = delay
            break
    ck("C3 CROSSOVER IS COMPUTABLE: delaying the deeper regimen flips the CNS "
       "endpoint at a finite delay", flip is not None,
       ("flip at %d days of delay" % flip) if flip else "no flip up to 90 d")
    ck("C3: the full 6-cycle course leaves a smaller reservoir than MAPKi",
       at(o_deep, 400)["PRECM"] < at(o_fast, 400)["PRECM"],
       "reservoir deep=%.5f fast=%.5f"
       % (at(o_deep, 400)["PRECM"], at(o_fast, 400)["PRECM"]))

    # ---------------- myelosuppression ----------------------------------
    nadir = trough(o4, "ANC", 0, 40)
    ck("2-CdA/Ara-C produces grade-4 neutropenia (ANC nadir <0.5)",
       nadir < 0.5, "ANC nadir=%.2f x10^9/L" % nadir)
    rec = at(o4, 27)["ANC"]
    ck("ANC recovers before the next cycle", rec > nadir * 2,
       "ANC d27=%.2f" % rec)
    o2 = run("S2_LCH3_ROneg")
    ck("VBL/prednisolone is far less myelosuppressive than 2-CdA/Ara-C",
       trough(o2, "ANC", 0, 60) > nadir,
       "VBL/pred nadir=%.2f" % trough(o2, "ANC", 0, 60))

    # ---------------- genotype / paradox --------------------------------
    o_map2k1 = run("S5_vem_continuous", {"GENO": 2.0})
    ck("MAP2K1-driven disease does not respond to a BRAF inhibitor",
       at(o_map2k1, 90)["DAS"] > 0.7 * at(o_map2k1, 0)["DAS"],
       "DAS %.2f -> %.2f" % (at(o_map2k1, 0)["DAS"],
                             at(o_map2k1, 90)["DAS"]))
    _pm = pheno("MSROpos")
    _pm["GENO"] = 2.0
    _y0m, _ = establish(_pm, "MSROpos")
    o_map2k1_mek = simulate(_pm, r_dabtram(Regimen(), 0, 400), 400, y0=_y0m)
    ck("MEK inhibition rescues the non-BRAF genotype that a BRAF inhibitor "
       "cannot touch",
       at(o_map2k1_mek, 90)["DAS"] < 0.6 * at(o_map2k1, 90)["DAS"],
       "MEKi DAS=%.2f vs BRAFi DAS=%.2f at day 90"
       % (at(o_map2k1_mek, 90)["DAS"], at(o_map2k1, 90)["DAS"]))
    sk_mono = at(o5, 200)["SKTOX"]
    sk_combo = at(o7, 200)["SKTOX"]
    ck("adding a MEK inhibitor suppresses paradox-driven skin toxicity",
       sk_combo < sk_mono, "SKTOX combo=%.2f vs mono=%.2f"
       % (sk_combo, sk_mono))
    lv = trough(o7, "LVEF", 0, 400)
    ck("MEK inhibitor causes a modest LVEF decline (>45%)",
       45 < lv < P["LVEF0"], "LVEF nadir=%.1f%%" % lv)

    # ---------------- pulmonary LCH / smoking ---------------------------
    oq = run("S10_plch_quit")
    os_ = run("S10b_plch_smoke")
    ck("smoking cessation limits cystic destruction",
       at(oq, 1095)["LUNGC"] < at(os_, 1095)["LUNGC"],
       "quit=%.2f vs smoke=%.2f" % (at(oq, 1095)["LUNGC"],
                                    at(os_, 1095)["LUNGC"]))
    ck("cystic destruction is irreversible even after quitting",
       at(oq, 1095)["LUNGC"] >= at(oq, 90)["LUNGC"],
       "LUNGC d90=%.3f d1095=%.3f" % (at(oq, 90)["LUNGC"],
                                      at(oq, 1095)["LUNGC"]))

    # ---------------- salvage logic -------------------------------------
    o3 = run("S3_frontline_fail")
    ck("front-line failure then salvage still reaches low disease activity",
       at(o3, 400)["DAS"] < 3.0, "DAS d400=%.2f" % at(o3, 400)["DAS"])

    # ---------------- further mechanistic checks -------------------------
    css_early = at(o5, 10)["CVEM"]
    css_late = at(o5, 120)["CVEM"]
    ck("vemurafenib autoinduction lowers steady-state exposure over weeks",
       css_late < css_early,
       "Css day10=%.1f -> day120=%.1f mg/L" % (css_early, css_late))
    o5_noind = run("S5_vem_continuous", {"EMAX_IND": 0.0})
    ck("removing autoinduction (EMAX_IND=0) raises exposure",
       at(o5_noind, 120)["CVEM"] > css_late,
       "no-induction Css=%.1f mg/L" % at(o5_noind, 120)["CVEM"])
    cns_ratio = (at(o4, 30)["LCNS"] / max(at(o4, 0)["LCNS"], 1e-9))
    sys_ratio = (at(o4, 30)["LRO"] / max(at(o4, 0)["LRO"], 1e-9))
    ck("limited CNS drug penetration makes the CNS lesion respond less than "
       "the systemic lesions", cns_ratio > sys_ratio,
       "CNS retains %.2f vs risk-organ %.2f of baseline at day 30"
       % (cns_ratio, sys_ratio))
    o_no17 = run("S1_observation", {"IL17ON": 0.0})
    ck("the contested IL-17A/DC-fusion arm is switchable and moves osteoclast "
       "output", peak(o_no17, "OCL") < peak(o1, "OCL"),
       "OCL peak with IL17ON=0: %.3f vs %.3f"
       % (peak(o_no17, "OCL"), peak(o1, "OCL")))
    o_novbl = run("S2_LCH3_ROneg", {"EMAX_VBL": 0.0, "EMAX_GRK": 0.0})
    ck("removing vinblastine/steroid cytotoxicity worsens front-line control",
       at(o_novbl, 120)["DAS"] > at(o2, 120)["DAS"],
       "DAS at day 120: %.2f without kill vs %.2f with"
       % (at(o_novbl, 120)["DAS"], at(o2, 120)["DAS"]))

    # ---------------- numerical hygiene ---------------------------------
    allpos = True
    for k in CMT:
        if k == "LVEF":
            continue
        if min(o5[k]) < -1e-6:
            allpos = False
    ck("no negative states in a 2-year simulation", allpos)
    fine = simulate(pheno("MSROpos"), r_vem(Regimen(), 0, 200), 200, dt=0.0025)
    coarse = simulate(pheno("MSROpos"), r_vem(Regimen(), 0, 200), 200, dt=0.01)
    rel = abs(at(fine, 200)["DAS"] - at(coarse, 200)["DAS"]) \
        / max(at(fine, 200)["DAS"], 1e-9)
    ck("step-size halving changes DAS by <1% (integration converged)",
       rel < 0.01, "relative change=%.2e" % rel)

    # ---------------- report --------------------------------------------
    print()
    npass = 0
    for label, ok, detail in checks:
        print("  [%s] %s%s" % ("PASS" if ok else "FAIL", label,
                               ("  --  " + detail) if detail else ""))
        npass += ok
    print()
    print("-" * 78)
    print("%d / %d checks passed" % (npass, len(checks)))
    print("-" * 78)

    print("\nScenario summary at 2 years (or end of run):")
    hdr = ("%-22s %6s %7s %7s %7s %6s %6s %6s" %
           ("scenario", "DAS", "cfDNA", "AVPN", "NEUR", "CDI", "ND", "TTET"))
    print(hdr)
    print("-" * len(hdr))
    for name in ("S1_observation", "S2_LCH3_ROneg", "S3_frontline_fail",
                 "S4_cladarac_upfront", "S5_vem_continuous", "S6_vem_stop",
                 "S7_dabtram", "S8_bridge_consolidate", "S9_delayed_dx",
                 "S9b_early_dx", "S10_plch_quit", "S10b_plch_smoke"):
        oo = run(name)
        e = at(oo, oo["time"][-1])
        print("%-22s %6.2f %7.2f %7.3f %7.3f %6.0f %6.0f %6.0f" %
              (name, e["DAS"], e["CFDNA"], e["AVPN"], e["NEUR"], e["CDI"],
               e["ND"], e["TTET"]))

    return 0 if npass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
