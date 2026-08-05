#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
icp_reference_model.py
======================
Dependency-free reference implementation of the Intrahepatic Cholestasis of
Pregnancy (ICP) QSP model.  This file exists because there is no R runtime in
the build container: every equation in `icp_mrgsolve_model.R` is first written
and integrated HERE (fixed-step RK4, pure Python, no numpy) so that the
published mrgsolve model is a translation of a system that has actually been
run, not a system that has only been typed.

WHY A SPECIES-RESOLVED MODEL
----------------------------
Almost every quantitative claim in the ICP literature is about *total* serum
bile acid (TBA), because that is what the routine 3alpha-hydroxysteroid
dehydrogenase assay reports.  But the two things clinicians actually care
about -- fetal death and itch -- are not driven by the total.  They are driven
by (i) the hydrophobic, cardiotoxic subset of the fetal bile acid pool and
(ii) an autotaxin/LPA axis that is only loosely coupled to bile acids at all.
A model whose state variable is "TBA" therefore cannot, even in principle,
reproduce the two central puzzles of the field:

  * the stillbirth risk threshold sits at TBA >= 100 umol/L and is *sharp*
    (Ovadia 2019: 0.13% / 0.28% / 3.44% for <40 / 40-99 / >=100 umol/L), and
  * ursodeoxycholic acid lowers TBA by ~half and yet did not improve perinatal
    outcome in the PITCHES randomised trial (Chappell 2019).

So the model carries five bile acid species (CA, CDCA, DCA, LCA, UDCA) through
six maternal compartments and three fetal ones, plus a hyocholic-acid sink for
CYP3A4 6alpha-hydroxylation.  Both puzzles then stop being puzzles and become
arithmetic; see `icp_calibration.py` and README.md.

STATE VECTOR (78 ODEs)
----------------------
  0-4    HL[5]   hepatocyte bile acid          (umol, effective Vd VHL)
  5-9    BI[5]   biliary tree + gallbladder    (umol)
  10-14  IL[5]   small-intestinal lumen        (umol)
  15-19  CO[5]   colonic lumen                 (umol)
  20-24  SP[5]   maternal systemic plasma      (umol)
  25-29  FP[5]   fetal plasma + ECF            (umol)
  30-34  FL[5]   fetal liver                   (umol)
  35     AF      amniotic fluid total BA       (umol)
  36     MEC     meconium bile acid (cumulative, umol)
  37     CYP7A1  cholesterol 7alpha-hydroxylase (relative, 1 = normal)
  38     SHP     small heterodimer partner      (relative)
  39     FGF19   ileal FGF19 -> portal          (relative)
  40     E2      estradiol                      (nmol/L)
  41     P4      progesterone                   (nmol/L)
  42     P4S     progesterone sulfates          (nmol/L)
  43     E2G     estradiol-17beta-glucuronide   (nmol/L)
  44     ATX     autotaxin                      (relative, 1 = normal)
  45     LPA     lysophosphatidic acid          (relative)
  46     ITCHC   central itch / scratch state   (relative)
  47     SLEEP   sleep disturbance index
  48     ROS     hepatocyte oxidative stress
  49     ALT     alanine aminotransferase       (U/L)
  50     GJ      fetal cardiac Cx43 coupling    (relative, 1 = normal)
  51     FCA     fetal cardiomyocyte Ca index   (relative, 1 = normal)
  52     VASO    chorionic vasoconstriction
  53     TROPH   trophoblast damage
  54     HYP     fetal hypoxia index
  55     OTR     myometrial oxytocin receptor   (relative)
  56     UTA     uterine activity index
  57     VITK    maternal vitamin K status      (relative, 1 = normal)
  58     PIVKA   under-carboxylated factor II
  59     HSB     cumulative stillbirth hazard
  60     HPT     cumulative spontaneous preterm birth hazard
  61     HMEC    cumulative meconium-passage hazard
  62     UDDEP   UDCA gut depot                 (umol)
  63     RIFDEP  rifampicin gut depot           (mg)
  64     RIFC    rifampicin central             (mg)
  65     RIFP    rifampicin peripheral          (mg)
  66     CYP3A   CYP3A4 induction               (relative, 1 = normal)
  67     HCA     hyocholic/hyodeoxycholic pool  (umol, maternal plasma)
  68     CHOLL   cholestyramine luminal         (g)
  69     SAMC    S-adenosylmethionine central   (mg)
  70     ODEV    IBAT inhibitor luminal         (umol)
  71     NTXC    naltrexone central             (mg)
  72     AHC     antihistamine central          (mg)
  73     BETC    betamethasone central          (mg)
  74     SURF    fetal surfactant / lung maturity (relative)
  75     SULT    SULT2A1 induction              (relative)
  76     MRP4    basolateral MRP3/4/OSTab        (relative)
  77     TDOSE   cumulative UDCA dose bookkeeping (umol, diagnostic only)

Run:  python3 icp_reference_model.py            -> full scenario table
      python3 icp_reference_model.py --quick     -> 6 core scenarios only
"""

import math
import sys

# ---------------------------------------------------------------------------
# species bookkeeping
# ---------------------------------------------------------------------------
SPECIES = ["CA", "CDCA", "DCA", "LCA", "UDCA"]
NS = 5
iCA, iCDCA, iDCA, iLCA, iUDCA = 0, 1, 2, 3, 4

# Cytotoxicity / detergency weights.  Ordering and spacing follow the
# hydrophobicity (critical micellar concentration) series and the fetal
# cardiomyocyte data: LCA > DCA > CDCA >> CA >>> UDCA.  UDCA is not zero
# because taurine-conjugated UDCA is measurably (weakly) detergent, but it is
# two orders below LCA, and it is the reason UDCA can lower the assay without
# lowering the toxic load.
WTOX = [0.20, 0.60, 0.72, 1.00, 0.02]
# FXR agonist potency (CDCA is the physiological ligand; UDCA is ~inert)
WFXR = [0.30, 1.00, 0.60, 0.50, 0.02]

# offsets
oHL, oBI, oIL, oCO, oSP, oFP, oFL = 0, 5, 10, 15, 20, 25, 30
(iAF, iMEC, iCYP7A1, iSHP, iFGF19, iE2, iP4, iP4S, iE2G, iATX, iLPA, iITCHC,
 iSLEEP, iROS, iALT, iGJ, iFCA, iVASO, iTROPH, iHYP, iOTR, iUTA, iVITK,
 iPIVKA, iHSB, iHPT, iHMEC, iUDDEP, iRIFDEP, iRIFC, iRIFP, iCYP3A, iHCA,
 iCHOLL, iSAMC, iODEV, iNTXC, iAHC, iBETC, iSURF, iSULT, iMRP4,
 iTDOSE) = range(35, 78)
NST = 78

MW_UDCA = 392.57   # g/mol
MW_RIF = 822.9

# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------
P = dict(
    # ---- volumes -----------------------------------------------------------
    # VHL is an EFFECTIVE hepatocyte distribution volume, not cytosolic water:
    # intracellular bile acids are largely bound to 3alpha-HSD/FABP-class
    # proteins, so the transport-relevant free concentration corresponds to a
    # much larger apparent volume.  Using cytosolic water (1.5 L) makes the
    # hepatocyte turnover ~500/day and the system numerically stiff for no
    # physiological gain.
    VHL=15.0, VBI=0.06, VIL=1.0, VCO=0.5,
    # VSP is the APPARENT volume of distribution of total bile acid, not the
    # plasma volume: >95% of circulating bile acid is albumin-bound and the
    # bound pool equilibrates with interstitial fluid within minutes, so the
    # measured serum concentration behaves as amount/Vd with Vd ~ 4x plasma
    # volume.  Using plasma volume instead puts a 7-minute half-life on a
    # 4 L compartment and makes the system needlessly stiff.
    VSP=18.0,

    # ---- synthesis and its feedback ---------------------------------------
    VSYN=520.0,          # umol/d de-novo synthesis at CYP7A1 = 1
    FCA=0.55,            # CYP8B1 branch: fraction to cholic acid
    KS_C7=1.0, KD_C7=1.0, KSHP=1.0, KF19=1.2,
    KS_SHP=2.2, KD_SHP=2.2, EC50_FXR=120.0,
    KS_F19=3.0, KD_F19=3.0, EC50_FXRI=900.0,

    # ---- MRP2/ABCC2: the canalicular route for UDCA -----------------------
    # UDCA and its conjugates are exported canalicularly largely by MRP2, not
    # BSEP.  Modelling UDCA as a BSEP substrate inside the shared denominator
    # forces an unphysiological choice: either UDCA cannot reach bile (and
    # accumulates to plasma concentrations of 60+ umol/L that are never seen),
    # or it displaces the endogenous species and RAISES their serum level,
    # which is the opposite of the observed effect of the drug.  Giving UDCA
    # its own non-competing route resolves both.
    VMRP2=24000.0, KM_MRP2=55.0, EMRP2P=0.60,

    # ---- BSEP (ABCB11) ----------------------------------------------------
    # BSEP operates well below saturation in health (denB ~ 1.9): if it is
    # placed near Vmax the model amplifies a 1.3-fold fall in transport
    # capacity into a 6-fold rise in hepatocyte load, which overshoots the
    # clinical spread of ICP by an order of magnitude.
    VBSEP=74000.0,
    KMB=[52.0, 42.0, 42.0, 28.0, 1e9],
    ETRAF=0.07, KTRAF=6.0,        # UDCA-driven canalicular carrier insertion
    ESAM=0.18, KSAM=9.0,          # SAMe: membrane fluidity / transmethylation

    # ---- sex-steroid inhibition of BSEP -----------------------------------
    KI_E2G=11.0,      # nmol/L, estradiol-17beta-glucuronide (cis-inhibition)
    KI_P4S=140.0,     # nmol/L, progesterone sulfates (trans-inhibition)

    # ---- basolateral escape valve (MRP3/MRP4/OSTalpha-beta) ---------------
    KBL=39.0,                       # L/d at MRP4 = 1
    SBL=[1.00, 0.90, 0.85, 0.60, 0.18],
    KS_M4=1.0, KD_M4=1.0, EM4=7.0, KM4=0.62, HM4=2.5, EM4P=0.25,

    # ---- hepatic uptake / first pass --------------------------------------
    EMAXUP=[0.990, 0.985, 0.980, 0.970, 0.930],  # intrinsic first-pass
    JHALF=800000.0,     # umol/d, load at which first-pass efficiency halves
    KHLSAT=300.0,       # umol/L effective hepatocyte load that trans-inhibits
    HHLSAT=2.0,         # Hill: NTCP downregulation is a threshold, not linear
    KI_RIF_OATP=12.0,   # mg/L rifampicin (total), OATP1B1/NTCP inhibition
    CLHSYS=600.0,       # L/d hepatic clearance of systemic (recirculating) BA

    # ---- bile flow and MDR3 ----------------------------------------------
    KBILE=7.0,
    VPC=6000.0,         # umol/d biliary phosphatidylcholine at MDR3 = 1
    BSPC0=2.8,          # normal biliary bile salt : phosphatidylcholine ratio

    # ---- intestine --------------------------------------------------------
    VASBT=20500.0,
    KMA=[300.0, 260.0, 260.0, 200.0, 5000.0],
    KPASS=[0.004, 0.010, 0.026, 0.030, 0.007],   # /d jejunal passive
    KTR=0.42,           # /d ileum -> colon
    KSEQ=0.45, KCHOL=3.0,
    FBIND=[1.00, 1.05, 1.05, 1.10, 0.95],
    KI_ODEV=0.9,        # umol luminal odevixibat equivalent, ASBT inhibition
    KDH_CA=1.10, KDH_CD=1.25, KDH_UD=0.20,   # /d microbial 7-dehydroxylation
    # Hepatic LCA detoxification: SULT2A1 sulfation (plus CYP3A4 6alpha-OH,
    # handled separately).  Without this route the model has no way to dispose
    # of the lithocholate that gut bacteria make from UDCA, and it then
    # overstates the hydrophobicity shift that UDCA therapy causes.  The same
    # enzyme is PXR-inducible, which is part of why rifampicin is protective.
    VLCAS=900.0, KM_LCAS=6.0,
    KCOL=[0.06, 0.10, 0.55, 0.60, 0.08],     # /d colonic passive absorption
    KFEC=1.05,

    # ---- renal ------------------------------------------------------------
    CLREN=[0.170, 0.160, 0.145, 0.095, 3.000],   # L/d
    ESULT=4.0,
    KS_S=1.0, KD_S=1.0, EPXR_S=0.85,

    # ---- sex steroids -----------------------------------------------------
    E2BASE=20.0, E2K=0.0627, KOUT_E2=2.77,
    P4BASE=200.0, P4K=0.0550, KOUT_P4=1.39,
    KF_E2G=0.185, KOUT_E2G=3.10, EUGT=0.15,
    KF_P4S=0.088, KOUT_P4S=0.55,

    # ---- placental transfer ----------------------------------------------
    # VP is the maximal fetal->maternal transfer capacity of the
    # syncytiotrophoblast.  Its SATURATION is the mechanism this model
    # proposes for the 100 umol/L stillbirth threshold, so KMP matters more
    # than any other single parameter here.
    VP=88.0,                                     # umol/d
    KMP=[26.0, 20.0, 20.0, 15.0, 55.0],          # umol/L
    TRANSIN=0.30,        # maternal-side occupancy of the same carrier
    PS=[0.550, 0.750, 0.850, 0.950, 0.150],      # L/d, BIDIRECTIONAL passive
    SF=13.0,             # umol/d fetal de-novo synthesis
    FSF=[0.60, 0.40, 0.0, 0.0, 0.0],
    CLFL=0.55, KFLO=2.2, KFB=0.55,               # fetal liver in/out/bile
    KUR=0.200, KSW=1.30, KIMA=0.55,              # urine / swallow / intramembr.

    # ---- fetal myocardium -------------------------------------------------
    IMAXGJ=0.85, IC50GJ=30.0, HGJ=1.60,
    KS_GJ=1.0, KD_GJ=1.0,
    ECA=2.50, KCA=45.0, KS_CA=1.0, KD_CA=1.0,

    # ---- placental vasculature / fetal hypoxia ---------------------------
    EV=1.0, KV=40.0, KS_V=1.0, KD_V=1.0,
    KT=55.0, KS_T=1.0, KD_T=1.0,
    KS_H=1.0, KD_H=1.0, WT=0.65,

    # ---- hazards (the ONLY constants fitted to outcome epidemiology) -----
    HSB0=1.857e-5,   # /d background stillbirth hazard  (fitted, Ovadia <40)
    HSBSC=4.36e-4,   # /d scale on arrhythmia index     (fitted, Ovadia 40-99)
    HN=2.0,          # exponent on arrhythmia index     (fitted, Ovadia >=100)
    KHY=0.55,
    HPT0=5.4e-4, KOT=2.20, KUT=1.00,
    HM0=3.7e-3, HMSC=2.52e-2, KMG=25.0,

    # ---- myometrium -------------------------------------------------------
    KS_O=1.0, KD_O=1.0, EO=1.60, KO=22.0,
    KS_U=1.0, KD_U=1.0,

    # ---- pruritus (autotaxin / LPA axis) ---------------------------------
    # Autotaxin is driven mainly by the sex-steroid load and only weakly, and
    # SATURABLY, by circulating bile acid.  This is not a convenience: serum
    # autotaxin does not fall with UDCA but does fall with rifampicin and with
    # nasobiliary drainage, itch severity correlates poorly with total bile
    # acid, and in a large minority of women itch PRECEDES the biochemical
    # abnormality by weeks.  A model that drives itch from TBA reproduces none
    # of those four observations.
    KS_A=1.0, KD_A=1.0, EA_E2=1.05, KA_E2=48.0,
    EA_P4S=3.00, KA_P4S=200.0, HP4S=2.0, EA_CH=0.35, KA_CH=3.0,
    ERIFA=0.62,
    KS_L=1.0, KD_L=1.0,
    # High apparent cooperativity downstream of LPA is the LPAR -> MRGPRX4 ->
    # TRPV1 -> spinal GRPR cascade plus the scratch-itch loop (KSENS), not a
    # single cooperative binding event.
    IMAXI=1.0, ILP50=3.20, HLP=4.0, WBA=0.30, KBA=25.0,
    KSENS=0.90, KSC=0.45,
    KS_I=1.0, KD_I=1.0, ENTX=0.30, EAH=0.22, VSCALE=4.81,
    KS_SL=1.0, KD_SL=1.0,

    # ---- hepatocellular injury -------------------------------------------
    KS_R=1.0, KD_R=1.0, KR=65.0, HR=3.0, WBSPC=0.62, EUD_ROS=0.42, KUD_ROS=7.0,
    ESAM_R=0.55, WRIF_ROS=0.28, KRIF_ROS=3.0,
    KS_ALT=20.0, KD_ALT=1.0, EALT=9.0,

    # ---- vitamin K / coagulation -----------------------------------------
    MIC50=1136.0, KIN_K=1.0, KD_K=1.0, DIETK=1.0, PIVKA0=0.278,
    KS_PI=1.0, KD_PI=1.0, KVK=0.62, WINR=0.42,

    # ---- fetal lung maturity ---------------------------------------------
    KS_SURF=1.0, KD_SURF=1.0, EBET=1.30, KBET=1.8,

    # ---- drug PK ----------------------------------------------------------
    KAUD=8.0, KTRUD=3.2, FUD=0.42,
    KA_RIF=15.0, F_RIF=0.90, CL_RIF=170.0, VC_RIF=45.0, Q_RIF=30.0, VP_RIF=20.0,
    EMAX3A=4.2, EC503A=1.6, KS_3A=0.75, KD_3A=0.75,
    V6A=260.0, F6A=[0.15, 1.00, 0.30, 0.90, 0.05],
    CLHCA=60.0,
    KCLR_CHOL=4.5,
    CL_SAM=48.0, V_SAM=22.0,
    KODEV=6.0,
    CL_NTX=190.0, V_NTX=1350.0, EC50_NTX=1.3,
    CL_AH=95.0, V_AH=340.0, EC50_AH=9.0,
    CL_BET=13.0, V_BET=42.0,
)

# genetic / constitutional susceptibility defaults
GEN0 = dict(gBSEP=1.0, gMDR3=1.0, gFIC1=1.0, gSULT=1.0, gFXR=1.0, gP=1.0,
            twin=1.0)


# ---------------------------------------------------------------------------
# gestational covariates
# ---------------------------------------------------------------------------
def fetal_weight(ga):
    """Hadlock-style log-quadratic fetal weight, kg."""
    lg = 0.578 + 0.332 * ga - 0.00354 * ga * ga
    return math.exp(lg) / 1000.0


def af_volume(ga):
    """Amniotic fluid volume, L (peak ~0.95 L near 33 wk)."""
    return 0.80 * math.exp(-((ga - 33.0) / 9.0) ** 2) + 0.15


def rds_risk(ga, surf):
    """Respiratory distress / NICU-for-respiratory-cause probability.

    Anchors (singleton, no antenatal steroid): 34 wk 0.22, 35 wk 0.12,
    36 wk 0.052, 37 wk 0.021, 38 wk 0.009, 39 wk 0.004, 40 wk 0.003.
    SURF (betamethasone) shifts the curve, it does not remove it.
    """
    base = 0.0028 + 1.0 / (1.0 + math.exp(1.15 * (ga - 33.6)))
    return base / (1.0 + 1.55 * max(0.0, surf - 1.0))


def nicu_risk(ga, surf, mec, hyp):
    """Any NICU admission >=4 h -- the PITCHES composite component that
    dominates the trial's primary endpoint."""
    r = 0.020 + 1.0 / (1.0 + math.exp(0.82 * (ga - 35.3)))
    r = r / (1.0 + 0.85 * max(0.0, surf - 1.0))
    r = r + 0.16 * mec + 0.22 * min(1.0, hyp / 1.2)
    return min(0.98, r)


# ---------------------------------------------------------------------------
# right-hand side
# ---------------------------------------------------------------------------
def derivs(t, y, p, ctl):
    d = [0.0] * NST
    ga = ctl["ga0"] + t / 7.0
    delivered = 1.0 if t >= ctl["tdel"] else 0.0
    PL = 0.0 if delivered else 1.0            # placental function switch
    g = ctl["gen"]

    fw = fetal_weight(ga)
    VFP = max(0.05, 0.25 * fw)
    VFL = max(0.01, 0.035 * fw)
    VAF = af_volume(ga)

    HL = y[oHL:oHL + NS]
    BI = y[oBI:oBI + NS]
    IL = y[oIL:oIL + NS]
    CO = y[oCO:oCO + NS]
    SP = y[oSP:oSP + NS]
    FP = y[oFP:oFP + NS]
    FL = y[oFL:oFL + NS]

    CHL = [max(0.0, HL[i]) / p["VHL"] for i in range(NS)]
    CIL = [max(0.0, IL[i]) / p["VIL"] for i in range(NS)]
    CSP = [max(0.0, SP[i]) / p["VSP"] for i in range(NS)]
    CFP = [max(0.0, FP[i]) / VFP for i in range(NS)]

    # ---------------- sex steroids -----------------------------------------
    E2tgt = p["E2BASE"] * math.exp(p["E2K"] * (ga - 20.0)) * g["twin"]
    P4tgt = p["P4BASE"] * math.exp(p["P4K"] * (ga - 20.0)) * g["twin"]
    E2, P4, P4S, E2G = y[iE2], y[iP4], y[iP4S], y[iE2G]
    d[iE2] = p["KOUT_E2"] * (E2tgt * PL - E2)
    d[iP4] = p["KOUT_P4"] * (P4tgt * PL - P4)
    CYP3A = max(1e-6, y[iCYP3A])
    # rifampicin co-induces UGT1A1 -> more E2-17G.  This is a real, and
    # partly self-defeating, wrinkle of using rifampicin in ICP and it is
    # kept in the model rather than hidden.
    d[iE2G] = p["KF_E2G"] * E2 * (1.0 + p["EUGT"] * (CYP3A - 1.0)) \
        - p["KOUT_E2G"] * E2G
    SULT = max(1e-6, y[iSULT])
    d[iP4S] = p["KF_P4S"] * P4 * g["gSULT"] * SULT - p["KOUT_P4S"] * P4S

    # steroid inhibitory pressure on BSEP
    SIP = E2G / p["KI_E2G"] + P4S / p["KI_P4S"]
    IBSEP = 1.0 / (1.0 + SIP)

    # ---------------- hepatic synthesis and its feedback -------------------
    FXRW = sum(WFXR[i] * CHL[i] for i in range(NS))
    FXRH = g["gFXR"] * FXRW / (p["EC50_FXR"] + FXRW)
    SHP, FGF19 = max(0.0, y[iSHP]), max(0.0, y[iFGF19])
    d[iSHP] = p["KS_SHP"] * FXRH - p["KD_SHP"] * SHP
    ILBA = sum(WFXR[i] * CIL[i] for i in range(NS))
    FXRI = g["gFXR"] * ILBA / (p["EC50_FXRI"] + ILBA)
    d[iFGF19] = p["KS_F19"] * FXRI - p["KD_F19"] * FGF19
    C7 = max(0.0, y[iCYP7A1])
    d[iCYP7A1] = p["KS_C7"] / (1.0 + (SHP / p["KSHP"]) ** 2
                               + FGF19 / p["KF19"]) - p["KD_C7"] * C7
    SYN = p["VSYN"] * C7
    S = [SYN * p["FCA"], SYN * (1.0 - p["FCA"]), 0.0, 0.0, 0.0]

    # ---------------- canalicular export (BSEP) ---------------------------
    SAMC = max(0.0, y[iSAMC])
    TRAF = 1.0 + p["ETRAF"] * CHL[iUDCA] / (p["KTRAF"] + CHL[iUDCA])
    SAMEF = 1.0 + p["ESAM"] * SAMC / (p["KSAM"] + SAMC)
    # BSEP carries the four endogenous species and competes only among them
    denB = 1.0 + sum(CHL[j] / p["KMB"][j] for j in range(4))
    gB = g["gBSEP"] * g["gFIC1"]
    JB = [p["VBSEP"] * gB * IBSEP * TRAF * SAMEF
          * (CHL[i] / p["KMB"][i]) / denB for i in range(4)] + [0.0]
    # UDCA: MRP2, independent of BSEP and up-regulated by PXR (rifampicin)
    JB[iUDCA] = p["VMRP2"] * (1.0 + p["EMRP2P"] * (CYP3A - 1.0)) \
        * CHL[iUDCA] / (p["KM_MRP2"] + CHL[iUDCA])

    # ---------------- LCA detoxification (SULT2A1 sulfation) --------------
    JLCAS = p["VLCAS"] * SULT * CHL[iLCA] / (p["KM_LCAS"] + CHL[iLCA])

    # ---------------- 6alpha-hydroxylation (CYP3A4) -----------------------
    J6A = [p["V6A"] * CYP3A * p["F6A"][i] * CHL[i] / 100.0 for i in range(NS)]

    # ---------------- basolateral escape valve ----------------------------
    MRP4 = max(0.0, y[iMRP4])
    # MRP3/MRP4/OSTalpha-beta are induced 5-20x in human cholestasis, but
    # they are near-silent in health; a proportional FXR term cannot express
    # both, so the induction is written as a Hill in FXR occupancy.
    rm4 = (FXRH / p["KM4"]) ** p["HM4"]
    d[iMRP4] = p["KS_M4"] * (1.0 + p["EM4"] * rm4 / (1.0 + rm4)
                             + p["EM4P"] * (CYP3A - 1.0)) \
        - p["KD_M4"] * MRP4
    JBL = [p["KBL"] * MRP4 * p["SBL"][i] * CHL[i] for i in range(NS)]

    # ---------------- intestinal handling ---------------------------------
    CHOLL = max(0.0, y[iCHOLL])
    ODEV = max(0.0, y[iODEV])
    denA = 1.0 + sum(CIL[j] / p["KMA"][j] for j in range(NS)) \
        + ODEV / p["KI_ODEV"]
    JASBT = [p["VASBT"] * (CIL[i] / p["KMA"][i]) / denA
             + p["KPASS"][i] * IL[i] for i in range(NS)]
    fseq = p["KSEQ"] * CHOLL / (p["KCHOL"] + CHOLL)
    JSEQ = [fseq * p["FBIND"][i] * IL[i] for i in range(NS)]

    # ---------------- colonic microbial transformation --------------------
    rCA = p["KDH_CA"] * max(0.0, CO[iCA])
    rCD = p["KDH_CD"] * max(0.0, CO[iCDCA])
    rUD = p["KDH_UD"] * max(0.0, CO[iUDCA])
    JCOL = [p["KCOL"][i] * max(0.0, CO[i]) for i in range(NS)]

    # ---------------- hepatic first pass ----------------------------------
    ABS = [JASBT[i] + JCOL[i] for i in range(NS)]
    ABStot = sum(ABS)
    CRIF = max(0.0, y[iRIFC]) / p["VC_RIF"]
    CHLtot = sum(CHL)
    # Sinusoidal uptake efficiency.  NTCP/OATP1B1 are transcriptionally
    # DOWN-regulated once the hepatocyte is loaded (a protective response), so
    # the term is a threshold in CHLtot, not a proportionality: first pass is
    # essentially complete in health and degrades only once cholestasis is
    # established.  This is what turns a linear fall in BSEP capacity into a
    # supralinear rise in maternal plasma bile acid.
    updep = 1.0 + ABStot / p["JHALF"] + CRIF / p["KI_RIF_OATP"] \
        + (CHLtot / p["KHLSAT"]) ** p["HHLSAT"]
    E = [p["EMAXUP"][i] / updep for i in range(NS)]
    JU = [E[i] * ABS[i] for i in range(NS)]        # -> hepatocyte
    JSPILL = [(1.0 - E[i]) * ABS[i] for i in range(NS)]   # -> systemic
    JUS = [p["CLHSYS"] / updep * CSP[i] for i in range(NS)]   # systemic uptake

    # ---------------- renal -----------------------------------------------
    JREN = [p["CLREN"][i] * (1.0 + p["ESULT"] * (SULT - 1.0)) * CSP[i]
            for i in range(NS)]
    d[iSULT] = p["KS_S"] * (1.0 + p["EPXR_S"] * (CYP3A - 1.0)) \
        - p["KD_S"] * SULT

    # ---------------- placental transfer ----------------------------------
    denP = 1.0 + sum(CFP[j] / p["KMP"][j] for j in range(NS)) \
        + p["TRANSIN"] * sum(CSP[j] / p["KMP"][j] for j in range(NS))
    JF2M = [p["VP"] * g["gP"] * PL * (CFP[i] / p["KMP"][i]) / denP
            for i in range(NS)]
    # The non-carrier component of placental transfer is DIFFUSIVE and
    # therefore bidirectional.  Writing it as a one-way maternal->fetal influx
    # (as several published ICP models do) leaves the fetal compartment with no
    # concentration-independent escape once the carrier saturates, and the
    # fetal pool then grows without bound.  With the gradient term present the
    # fetal concentration saturates just above the maternal one, which is what
    # cord-blood series actually show in very severe ICP.
    JPASS = [p["PS"][i] * PL * (CSP[i] - CFP[i]) for i in range(NS)]

    # ---------------- fetal liver / gut / amniotic fluid ------------------
    JFLIN = [p["CLFL"] * CFP[i] for i in range(NS)]
    JFLOUT = [p["KFLO"] * max(0.0, FL[i]) for i in range(NS)]
    JFBILE = [p["KFB"] * max(0.0, FL[i]) for i in range(NS)]
    FPtot = sum(CFP)
    JUR = p["KUR"] * FPtot * VFP
    AF = max(0.0, y[iAF])
    d[iAF] = JUR - (p["KSW"] + p["KIMA"]) * AF
    d[iMEC] = sum(JFBILE) * PL

    # ---------------- species mass balances -------------------------------
    for i in range(NS):
        d[oHL + i] = S[i] + JU[i] + JUS[i] - JB[i] - JBL[i] - J6A[i]
        d[oBI + i] = JB[i] - p["KBILE"] * max(0.0, BI[i])
        d[oIL + i] = p["KBILE"] * max(0.0, BI[i]) - JASBT[i] - JSEQ[i] \
            - p["KTR"] * max(0.0, IL[i])
        d[oCO + i] = p["KTR"] * max(0.0, IL[i]) - JCOL[i] \
            - p["KFEC"] * max(0.0, CO[i])
        d[oSP + i] = JSPILL[i] + JBL[i] - JUS[i] - JREN[i] \
            + JF2M[i] - JPASS[i]
        d[oFP + i] = JPASS[i] - JF2M[i] + JFLOUT[i] - JFLIN[i]
        d[oFL + i] = p["SF"] * p["FSF"][i] * PL + JFLIN[i] - JFLOUT[i] \
            - JFBILE[i]
    # fetal urine leaves fetal plasma proportionally
    for i in range(NS):
        d[oFP + i] -= p["KUR"] * CFP[i] * VFP
    d[oHL + iLCA] -= JLCAS
    # microbial 7-dehydroxylation in the colon
    d[oCO + iCA] -= rCA
    d[oCO + iDCA] += rCA
    d[oCO + iCDCA] -= rCD
    d[oCO + iLCA] += rCD
    d[oCO + iUDCA] -= rUD
    d[oCO + iLCA] += rUD
    # 6alpha-hydroxylation product leaves the five-species system
    d[iHCA] = sum(J6A) + JLCAS - p["CLHCA"] * max(0.0, y[iHCA]) / p["VSP"]

    # ---------------- derived exposure indices ----------------------------
    FCL = sum(WTOX[i] * CFP[i] for i in range(NS))     # fetal cytotoxic load
    MCL = sum(WTOX[i] * CSP[i] for i in range(NS))     # maternal cytotoxic load
    HSTRESS = sum(WTOX[i] * CHL[i] for i in range(NS))

    # ---------------- fetal myocardium ------------------------------------
    GJ, FCA = max(0.0, y[iGJ]), max(0.0, y[iFCA])
    r = (FCL / p["IC50GJ"]) ** p["HGJ"] if FCL > 0 else 0.0
    GJinh = p["IMAXGJ"] * r / (1.0 + r)
    d[iGJ] = p["KS_GJ"] * (1.0 - GJinh) - p["KD_GJ"] * GJ
    d[iFCA] = p["KS_CA"] * (1.0 + p["ECA"] * FCL / (p["KCA"] + FCL)) \
        - p["KD_CA"] * FCA
    ARRI = max(0.0, 1.0 - GJ) * FCA

    # ---------------- placental vasculature, fetal hypoxia ---------------
    VASO, TROPH, HYP = max(0.0, y[iVASO]), max(0.0, y[iTROPH]), max(0.0, y[iHYP])
    d[iVASO] = p["KS_V"] * p["EV"] * FCL / (p["KV"] + FCL) - p["KD_V"] * VASO
    d[iTROPH] = p["KS_T"] * FCL / (p["KT"] + FCL) - p["KD_T"] * TROPH
    d[iHYP] = p["KS_H"] * (VASO + p["WT"] * TROPH) - p["KD_H"] * HYP

    # ---------------- hazards ---------------------------------------------
    # Gating: "stillbirth" is only defined from 24 weeks, so the hazard is not
    # allowed to accrue before then -- otherwise the fitted background hazard
    # absorbs 4 weeks of pre-viable gestation and no longer means what the
    # epidemiology means by it.
    gate_sb = 1.0 if ga >= 24.0 else 0.0
    hSB = p["HSB0"] + p["HSBSC"] * (ARRI ** p["HN"]) * (1.0 + p["KHY"] * HYP)
    d[iHSB] = hSB * PL * gate_sb
    OTR, UTA = max(0.0, y[iOTR]), max(0.0, y[iUTA])
    d[iOTR] = p["KS_O"] * (1.0 + p["EO"] * MCL / (p["KO"] + MCL)) \
        - p["KD_O"] * OTR
    d[iUTA] = p["KS_U"] * max(0.0, OTR - 1.0) - p["KD_U"] * UTA
    gate_pt = 1.0 if 24.0 <= ga < 37.0 else 0.0
    hPT = p["HPT0"] * (1.0 + p["KOT"] * max(0.0, OTR - 1.0)) \
        * (1.0 + p["KUT"] * UTA)
    d[iHPT] = hPT * PL * gate_pt
    # Meconium passage in utero is a near-term phenomenon (gut motility and
    # TGR5 expression both mature late), so the hazard is gated on GA rather
    # than allowed to accrue from the second trimester.
    gate_mec = 1.0 / (1.0 + math.exp(-1.2 * (ga - 35.0)))
    hMEC = (p["HM0"] + p["HMSC"] * FCL / (p["KMG"] + FCL) * (1.0 + HYP)) \
        * gate_mec
    d[iHMEC] = hMEC * PL

    # ---------------- pruritus --------------------------------------------
    ATX, LPA, ITCHC = max(0.0, y[iATX]), max(0.0, y[iLPA]), max(0.0, y[iITCHC])
    rifatx = 1.0 / (1.0 + p["ERIFA"] * (CYP3A - 1.0))
    rp4s = (P4S / p["KA_P4S"]) ** p["HP4S"] if P4S > 0 else 0.0
    d[iATX] = p["KS_A"] * (1.0 + p["EA_E2"] * E2 / (p["KA_E2"] + E2)
                           + p["EA_P4S"] * rp4s / (1.0 + rp4s)
                           + p["EA_CH"] * MCL / (p["KA_CH"] + MCL)) * rifatx \
        - p["KD_A"] * ATX
    d[iLPA] = p["KS_L"] * ATX - p["KD_L"] * LPA
    NTXocc = max(0.0, y[iNTXC]) / p["V_NTX"] * 1000.0
    NTXocc = NTXocc / (p["EC50_NTX"] + NTXocc)
    AHocc = max(0.0, y[iAHC]) / p["V_AH"] * 1000.0
    AHocc = AHocc / (p["EC50_AH"] + AHocc)
    lr = (LPA / p["ILP50"]) ** p["HLP"] if LPA > 0 else 0.0
    itch_raw = p["IMAXI"] * lr / (1.0 + lr) + p["WBA"] * MCL / (p["KBA"] + MCL)
    sens = 1.0 + p["KSENS"] * ITCHC / (p["KSC"] + ITCHC)
    d[iITCHC] = p["KS_I"] * itch_raw * sens * (1.0 - p["ENTX"] * NTXocc) \
        * (1.0 - p["EAH"] * AHocc) - p["KD_I"] * ITCHC
    d[iSLEEP] = p["KS_SL"] * ITCHC - p["KD_SL"] * max(0.0, y[iSLEEP])

    # ---------------- hepatocellular injury -------------------------------
    PCFLUX = p["VPC"] * g["gMDR3"]
    # UDCA is excluded from the biliary bile-salt : phosphatidylcholine ratio:
    # it is not detergent at biliary concentrations, and this is precisely why
    # it protects the ABCB4-deficient canaliculus.  Counting it would make the
    # drug look as if it worsened the very lesion it is given to treat.
    BSPC = sum(JB[i] for i in range(4)) / max(PCFLUX, 1.0)
    ROS = max(0.0, y[iROS])
    prot_ud = 1.0 - p["EUD_ROS"] * CHL[iUDCA] / (p["KUD_ROS"] + CHL[iUDCA])
    rr = (HSTRESS / p["KR"]) ** p["HR"] if HSTRESS > 0 else 0.0
    # Rifampicin carries its own recognised hepatotoxic signal, so the drug is
    # not allowed to drive ALT to normal purely by depleting CDCA.
    d[iROS] = p["KS_R"] * (rr / (1.0 + rr)
                           + p["WBSPC"] * max(0.0, BSPC / p["BSPC0"] - 1.0)
                           + p["WRIF_ROS"] * CRIF / (p["KRIF_ROS"] + CRIF)) \
        * prot_ud - p["KD_R"] * ROS \
        * (1.0 + p["ESAM_R"] * SAMC / (p["KSAM"] + SAMC))
    d[iALT] = p["KS_ALT"] * (1.0 + p["EALT"] * ROS) \
        - p["KD_ALT"] * max(0.0, y[iALT])

    # ---------------- vitamin K -------------------------------------------
    ILtot = sum(CIL)
    MICELLE = ILtot / (p["MIC50"] + ILtot)
    VITK = max(1e-6, y[iVITK])
    d[iVITK] = p["KIN_K"] * (p["DIETK"] + ctl["vkdose"]) * MICELLE / 0.5 \
        - p["KD_K"] * VITK
    d[iPIVKA] = p["KS_PI"] / (1.0 + (VITK / p["KVK"]) ** 2) \
        - p["KD_PI"] * max(0.0, y[iPIVKA])

    # ---------------- fetal lung maturity ---------------------------------
    BETC = max(0.0, y[iBETC])
    d[iSURF] = p["KS_SURF"] * (1.0 + p["EBET"] * BETC / p["V_BET"]
                               / (p["KBET"] + BETC / p["V_BET"])) \
        - p["KD_SURF"] * max(0.0, y[iSURF])

    # ---------------- drug PK ---------------------------------------------
    UDDEP = max(0.0, y[iUDDEP])
    d[iUDDEP] = ctl["ud_rate"] - (p["KAUD"] + p["KTRUD"]) * UDDEP
    # absorbed UDCA is presented to the liver; unabsorbed reaches the colon
    d[oHL + iUDCA] += p["KAUD"] * UDDEP * p["FUD"]
    d[oSP + iUDCA] += p["KAUD"] * UDDEP * p["FUD"] * (1.0 - p["FUD"]) * 0.25
    d[oCO + iUDCA] += p["KTRUD"] * UDDEP
    d[iTDOSE] = ctl["ud_rate"]

    RIFDEP = max(0.0, y[iRIFDEP])
    d[iRIFDEP] = ctl["rif_rate"] - p["KA_RIF"] * RIFDEP
    d[iRIFC] = p["KA_RIF"] * RIFDEP * p["F_RIF"] \
        - p["CL_RIF"] / p["VC_RIF"] * max(0.0, y[iRIFC]) \
        - p["Q_RIF"] / p["VC_RIF"] * max(0.0, y[iRIFC]) \
        + p["Q_RIF"] / p["VP_RIF"] * max(0.0, y[iRIFP])
    d[iRIFP] = p["Q_RIF"] / p["VC_RIF"] * max(0.0, y[iRIFC]) \
        - p["Q_RIF"] / p["VP_RIF"] * max(0.0, y[iRIFP])
    d[iCYP3A] = p["KS_3A"] * (1.0 + p["EMAX3A"] * CRIF
                              / (p["EC503A"] + CRIF)) - p["KD_3A"] * CYP3A

    d[iCHOLL] = ctl["chol_rate"] - p["KCLR_CHOL"] * CHOLL
    d[iSAMC] = ctl["sam_rate"] - p["CL_SAM"] / p["V_SAM"] * SAMC
    d[iODEV] = ctl["odev_rate"] - p["KODEV"] * ODEV
    d[iNTXC] = ctl["ntx_rate"] - p["CL_NTX"] / p["V_NTX"] * max(0.0, y[iNTXC])
    d[iAHC] = ctl["ah_rate"] - p["CL_AH"] / p["V_AH"] * max(0.0, y[iAHC])
    d[iBETC] = ctl["bet_rate"](t) - p["CL_BET"] / p["V_BET"] * BETC

    return d


# ---------------------------------------------------------------------------
# observations
# ---------------------------------------------------------------------------
def observe(t, y, p, ctl):
    ga = ctl["ga0"] + t / 7.0
    fw = fetal_weight(ga)
    VFP = max(0.05, 0.25 * fw)
    VFL = max(0.01, 0.035 * fw)
    CSP = [max(0.0, y[oSP + i]) / p["VSP"] for i in range(NS)]
    CFP = [max(0.0, y[oFP + i]) / VFP for i in range(NS)]
    CHL = [max(0.0, y[oHL + i]) / p["VHL"] for i in range(NS)]
    HCA = max(0.0, y[iHCA]) / p["VSP"]
    # The routine clinical assay is 3alpha-hydroxysteroid dehydrogenase based:
    # it counts UDCA and hyocholic acid along with the endogenous species.
    tba_assay = sum(CSP) + HCA
    tba_endo = sum(CSP) - CSP[iUDCA]
    FCL = sum(WTOX[i] * CFP[i] for i in range(NS))
    MCL = sum(WTOX[i] * CSP[i] for i in range(NS))
    ARRI = max(0.0, 1.0 - y[iGJ]) * max(0.0, y[iFCA])
    return dict(
        ga=ga, TBA=tba_assay, TBA_endo=tba_endo, UDCAp=CSP[iUDCA], HCA=HCA,
        CORD=sum(CFP), FCL=FCL, MCL=MCL, ARRI=ARRI,
        VAS=min(10.0, p["VSCALE"] * max(0.0, y[iITCHC])),
        ALT=y[iALT], ATX=y[iATX], LPA=y[iLPA],
        GJ=y[iGJ], FCA=y[iFCA], HYP=y[iHYP], VASO=y[iVASO],
        SB=1.0 - math.exp(-max(0.0, y[iHSB])),
        PTB=1.0 - math.exp(-max(0.0, y[iHPT])),
        MEC=1.0 - math.exp(-max(0.0, y[iHMEC])),
        INR=1.0 + p["WINR"] * max(0.0, y[iPIVKA] - p["PIVKA0"]),
        VITK=y[iVITK], SURF=y[iSURF], CYP3A=y[iCYP3A],
        E2=y[iE2], P4S=y[iP4S], E2G=y[iE2G],
        CSP=CSP, CFP=CFP, CHL=CHL,
        FGUT=sum(WTOX[i] * max(0.0, y[oFL + i]) for i in range(NS)) / VFL,
    )


# ---------------------------------------------------------------------------
# initial state
# ---------------------------------------------------------------------------
def initial_state(p):
    y = [0.0] * NST
    # a plausible pre-pregnancy enterohepatic pool, refined by the burn-in
    for i, amt in zip(range(NS), [520.0, 380.0, 210.0, 30.0, 8.0]):
        y[oHL + i] = amt
    for i, amt in zip(range(NS), [1500.0, 1050.0, 560.0, 80.0, 20.0]):
        y[oBI + i] = amt
    for i, amt in zip(range(NS), [1050.0, 730.0, 400.0, 60.0, 15.0]):
        y[oIL + i] = amt
    for i, amt in zip(range(NS), [180.0, 120.0, 320.0, 130.0, 8.0]):
        y[oCO + i] = amt
    for i, amt in zip(range(NS), [6.0, 4.0, 3.6, 0.6, 0.4]):
        y[oSP + i] = amt
    for i, amt in zip(range(NS), [1.4, 0.9, 0.35, 0.06, 0.03]):
        y[oFP + i] = amt
    for i, amt in zip(range(NS), [1.0, 0.7, 0.25, 0.04, 0.02]):
        y[oFL + i] = amt
    y[iAF] = 1.0
    y[iCYP7A1] = 0.6
    y[iSHP] = 0.6
    y[iFGF19] = 0.7
    y[iE2] = 20.0
    y[iP4] = 200.0
    y[iP4S] = 30.0
    y[iE2G] = 1.2
    y[iATX] = 1.6
    y[iLPA] = 1.6
    y[iITCHC] = 0.05
    y[iROS] = 0.2
    y[iALT] = 20.0
    y[iGJ] = 1.0
    y[iFCA] = 1.0
    y[iOTR] = 1.0
    y[iVITK] = 1.0
    y[iPIVKA] = 0.3
    y[iCYP3A] = 1.0
    y[iSURF] = 1.0
    y[iSULT] = 1.0
    y[iMRP4] = 1.0
    return y


# ---------------------------------------------------------------------------
# scenario control
# ---------------------------------------------------------------------------
def make_ctl(gen=None, ga0=20.0, tdel_ga=39.0, udca_mg=0.0, udca_start=30.0,
             rif_mg=0.0, rif_start=32.0, chol_g=0.0, chol_start=32.0,
             sam_mg=0.0, sam_start=32.0, odev_umol=0.0, odev_start=32.0,
             ntx_mg=0.0, ntx_start=32.0, ah_mg=0.0, ah_start=32.0,
             bet_ga=None, vkdose=0.0, vk_start=34.0):
    """Build a control dict.  Oral regimens are entered as daily input rates
    into the relevant depot; the mrgsolve version uses discrete dosing events.
    Because every reported endpoint integrates over weeks, the two agree to
    within the width of the plotted lines (checked in icp_calibration.py)."""
    g = dict(GEN0)
    if gen:
        g.update(gen)
    ctl = dict(gen=g, ga0=ga0, tdel=(tdel_ga - ga0) * 7.0,
               ud_mg=udca_mg, ud_t0=(udca_start - ga0) * 7.0,
               rif_mg=rif_mg, rif_t0=(rif_start - ga0) * 7.0,
               chol_g=chol_g, chol_t0=(chol_start - ga0) * 7.0,
               sam_mg=sam_mg, sam_t0=(sam_start - ga0) * 7.0,
               odev=odev_umol, odev_t0=(odev_start - ga0) * 7.0,
               ntx_mg=ntx_mg, ntx_t0=(ntx_start - ga0) * 7.0,
               ah_mg=ah_mg, ah_t0=(ah_start - ga0) * 7.0,
               vk_mg=vkdose, vk_t0=(vk_start - ga0) * 7.0,
               bet_t0=None if bet_ga is None else (bet_ga - ga0) * 7.0)
    ctl["ud_rate"] = 0.0
    ctl["rif_rate"] = 0.0
    ctl["chol_rate"] = 0.0
    ctl["sam_rate"] = 0.0
    ctl["odev_rate"] = 0.0
    ctl["ntx_rate"] = 0.0
    ctl["ah_rate"] = 0.0
    ctl["vkdose"] = 0.0
    ctl["bet_rate"] = lambda t: 0.0
    return ctl


def set_rates(ctl, t):
    """Turn the scheduled regimens on/off at the current time."""
    on = lambda t0: (t0 is not None and t >= t0 and t < ctl["tdel"])
    ctl["ud_rate"] = (ctl["ud_mg"] / MW_UDCA * 1000.0) if on(ctl["ud_t0"]) else 0.0
    ctl["rif_rate"] = ctl["rif_mg"] if on(ctl["rif_t0"]) else 0.0
    ctl["chol_rate"] = ctl["chol_g"] if on(ctl["chol_t0"]) else 0.0
    ctl["sam_rate"] = ctl["sam_mg"] if on(ctl["sam_t0"]) else 0.0
    ctl["odev_rate"] = ctl["odev"] if on(ctl["odev_t0"]) else 0.0
    ctl["ntx_rate"] = ctl["ntx_mg"] if on(ctl["ntx_t0"]) else 0.0
    ctl["ah_rate"] = ctl["ah_mg"] if on(ctl["ah_t0"]) else 0.0
    ctl["vkdose"] = ctl["vk_mg"] if on(ctl["vk_t0"]) else 0.0
    if ctl["bet_t0"] is None:
        ctl["bet_rate"] = lambda tt: 0.0
    else:
        t0 = ctl["bet_t0"]
        # 12 mg IM x 2, 24 h apart, delivered as two 6-hour inputs
        def br(tt, t0=t0):
            r = 0.0
            for k in (0.0, 1.0):
                if t0 + k <= tt < t0 + k + 0.25:
                    r += 12.0 / 0.25
            return r
        ctl["bet_rate"] = br


# ---------------------------------------------------------------------------
# integrator
# ---------------------------------------------------------------------------
def simulate(ctl, p=P, dt=0.01, tend=None, burnin=140.0, record_every=1.0):
    """Burn in on a fixed GA-20 background, then run the pregnancy forward."""
    y = initial_state(p)
    # ---- burn-in: hold GA at ga0 so the enterohepatic pool equilibrates ----
    bctl = dict(ctl)
    bctl["tdel"] = 1e9
    bctl["ud_rate"] = bctl["rif_rate"] = bctl["chol_rate"] = 0.0
    bctl["sam_rate"] = bctl["odev_rate"] = bctl["ntx_rate"] = 0.0
    bctl["ah_rate"] = bctl["vkdose"] = 0.0
    bctl["bet_rate"] = lambda t: 0.0
    bfix = dict(bctl)
    nb = int(burnin / 0.02)
    for k in range(nb):
        y = rk4(0.0, y, 0.02, p, bfix, freeze_ga=True)

    if tend is None:
        tend = ctl["tdel"] + 21.0     # follow 3 weeks postpartum
    n = int(round(tend / dt))
    out = []
    next_rec = 0.0
    for k in range(n + 1):
        t = k * dt
        if t >= next_rec - 1e-9:
            out.append(observe(t, y, p, ctl))
            next_rec += record_every
        if k == n:
            break
        set_rates(ctl, t)
        y = rk4(t, y, dt, p, ctl)
    return out, y


def rk4(t, y, h, p, ctl, freeze_ga=False):
    if freeze_ga:
        tt = 0.0
    else:
        tt = t
    k1 = derivs(tt, y, p, ctl)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(NST)]
    k2 = derivs(tt + 0.5 * h, y2, p, ctl)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(NST)]
    k3 = derivs(tt + 0.5 * h, y3, p, ctl)
    y4 = [y[i] + h * k3[i] for i in range(NST)]
    k4 = derivs(tt + h, y4, p, ctl)
    return [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
            for i in range(NST)]


# ---------------------------------------------------------------------------
# endpoint assembly
# ---------------------------------------------------------------------------
def endpoints(rows, ctl, p=P):
    tdel_ga = ctl["ga0"] + ctl["tdel"] / 7.0
    pre = [r for r in rows if r["ga"] <= tdel_ga + 1e-6]
    at = pre[-1]
    post = [r for r in rows if r["ga"] > tdel_ga]
    peak = max(pre, key=lambda r: r["TBA"])
    onset = None
    for r in pre:
        if r["TBA"] >= 10.0:
            onset = r["ga"]
            break
    # postpartum normalisation: first day with assay TBA < 10
    resolve = None
    for r in post:
        if r["TBA"] < 10.0:
            resolve = (r["ga"] - tdel_ga) * 7.0
            break
    sb = at["SB"]
    ptb_spont = at["PTB"] if tdel_ga >= 37.0 else 1.0
    iatro_pt = 1.0 if tdel_ga < 37.0 else 0.0
    mec = at["MEC"]
    rds = rds_risk(tdel_ga, at["SURF"])
    nicu = nicu_risk(tdel_ga, at["SURF"], mec, at["HYP"])
    # PITCHES composite: perinatal death OR preterm birth <37 wk OR NICU >=4 h
    ptb_any = min(1.0, iatro_pt + (0.0 if iatro_pt else ptb_spont))
    comp = 1.0 - (1.0 - sb) * (1.0 - ptb_any) * (1.0 - nicu)
    pph = 0.045 * (1.0 + 1.35 * max(0.0, at["INR"] - 1.05) / 0.15)
    return dict(
        ga_del=tdel_ga, TBA_del=at["TBA"], TBA_endo_del=at["TBA_endo"],
        TBA_peak=peak["TBA"], TBA_peak_ga=peak["ga"],
        TBA_endo_peak=peak["TBA_endo"],
        onset_ga=onset, resolve_d=resolve,
        UDCAp=at["UDCAp"], HCA=at["HCA"],
        CORD=at["CORD"], FCL=at["FCL"], FCL_peak=peak["FCL"],
        ARRI=at["ARRI"], GJ=at["GJ"], HYP=at["HYP"],
        VAS=at["VAS"], VAS_peak=max(r["VAS"] for r in pre),
        ATX=at["ATX"], ALT=at["ALT"], ALT_peak=max(r["ALT"] for r in pre),
        INR=at["INR"], VITK=at["VITK"],
        SB=sb, PTB_spont=ptb_spont, MEC=mec, RDS=rds, NICU=nicu,
        COMPOSITE=comp, PPH=pph, SURF=at["SURF"], CYP3A=at["CYP3A"],
    )


# ---------------------------------------------------------------------------
# scenario library
# ---------------------------------------------------------------------------
# Severity strata are DEFINED by the total bile acid band they land in (the
# same way the literature defines them: >=10 diagnostic, >=40 severe,
# >=100 very severe), not by a genotype chosen in advance.  These three
# susceptibility vectors were selected by scanning gBSEP/gMDR3/gSULT until
# the simulated term TBA fell in the middle of each band.
MILD = dict(gBSEP=0.735, gMDR3=0.62, gSULT=1.95)     # term TBA ~ 25
SEVERE = dict(gBSEP=0.650, gMDR3=0.50, gSULT=2.40)   # term TBA ~ 55
VSEVERE = dict(gBSEP=0.520, gMDR3=0.35, gSULT=3.10)  # term TBA ~ 120


def scenarios():
    S = []
    A = S.append
    A(("01 정상 임신 (reference)", make_ctl()))
    A(("02 ICP 경증 무치료", make_ctl(gen=MILD)))
    A(("03 ICP 중증 무치료", make_ctl(gen=SEVERE)))
    A(("04 ICP 최중증 (TBA>=100) 무치료", make_ctl(gen=VSEVERE)))
    A(("05 중증 + UDCA 500mg bid @30wk", make_ctl(gen=SEVERE, udca_mg=1000)))
    A(("06 중증 + UDCA 750mg bid @30wk", make_ctl(gen=SEVERE, udca_mg=1500)))
    A(("07 중증 + UDCA 늦은 시작 @34wk",
       make_ctl(gen=SEVERE, udca_mg=1000, udca_start=34.0)))
    A(("08 최중증 + UDCA 500mg bid", make_ctl(gen=VSEVERE, udca_mg=1000)))
    A(("09 중증 + 리팜피신 300mg bid", make_ctl(gen=SEVERE, rif_mg=600)))
    A(("10 중증 + UDCA + 리팜피신",
       make_ctl(gen=SEVERE, udca_mg=1000, rif_mg=600)))
    A(("11 최중증 + UDCA + 리팜피신",
       make_ctl(gen=VSEVERE, udca_mg=1000, rif_mg=600)))
    A(("12 중증 + 콜레스티라민 16g/d", make_ctl(gen=SEVERE, chol_g=16.0)))
    A(("13 중증 + 콜레스티라민 + 비타민K",
       make_ctl(gen=SEVERE, chol_g=16.0, vkdose=1.6)))
    A(("14 중증 + SAMe 1000mg IV", make_ctl(gen=SEVERE, sam_mg=1000)))
    A(("15 중증 + UDCA + SAMe",
       make_ctl(gen=SEVERE, udca_mg=1000, sam_mg=1000)))
    A(("16 최중증 + IBAT 억제제 (가설)", make_ctl(gen=VSEVERE, odev_umol=4.0)))
    A(("17 중증 + IBAT 억제제 (가설)", make_ctl(gen=SEVERE, odev_umol=4.0)))
    A(("18 중증 + 날트렉손 50mg", make_ctl(gen=SEVERE, ntx_mg=50.0)))
    A(("19 중증 + 항히스타민 12mg/d", make_ctl(gen=SEVERE, ah_mg=12.0)))
    A(("20 쌍태 임신 무치료", make_ctl(gen=dict(MILD, twin=1.55))))
    A(("21 쌍태 임신 + UDCA",
       make_ctl(gen=dict(MILD, twin=1.55), udca_mg=1000)))
    A(("22 ABCB11 V444A 동형접합",
       make_ctl(gen=dict(gBSEP=0.55, gMDR3=0.95, gSULT=2.2))))
    A(("23 ABCB4 중증 변이 (MDR3 30%)",
       make_ctl(gen=dict(gBSEP=0.92, gMDR3=0.30, gSULT=2.2))))
    A(("24 ATP8B1 변이", make_ctl(gen=dict(gFIC1=0.68, gMDR3=0.80, gSULT=2.4))))
    A(("25 NR1H4(FXR) 기능저하",
       make_ctl(gen=dict(gFXR=0.55, gBSEP=0.85, gMDR3=0.70, gSULT=2.6))))
    A(("26 태반 수송능 저하 (gP 0.6)",
       make_ctl(gen=dict(SEVERE, gP=0.60))))
    A(("27 최중증 + 36주 분만", make_ctl(gen=VSEVERE, tdel_ga=36.0)))
    A(("28 최중증 + 36주 분만 + 베타메타손",
       make_ctl(gen=VSEVERE, tdel_ga=36.0, bet_ga=35.0)))
    A(("29 최중증 + 37주 분만", make_ctl(gen=VSEVERE, tdel_ga=37.0)))
    A(("30 최중증 + 38주 분만", make_ctl(gen=VSEVERE, tdel_ga=38.0)))
    A(("31 최중증 + 40주 분만", make_ctl(gen=VSEVERE, tdel_ga=40.0)))
    A(("32 중증 + 37주 분만", make_ctl(gen=SEVERE, tdel_ga=37.0)))
    A(("33 중증 + 38주 분만", make_ctl(gen=SEVERE, tdel_ga=38.0)))
    A(("34 경증 + 39주 분만", make_ctl(gen=MILD, tdel_ga=39.0)))
    A(("35 최대 병용 (UDCA+RIF+IBAT+36주)",
       make_ctl(gen=VSEVERE, udca_mg=1500, rif_mg=600, odev_umol=4.0,
                tdel_ga=36.0, bet_ga=35.0)))
    return S


# ---------------------------------------------------------------------------
def main():
    quick = "--quick" in sys.argv
    S = scenarios()
    if quick:
        S = [S[i] for i in (0, 1, 2, 3, 4, 8)]
    hdr = ("{:<34s} {:>7s} {:>7s} {:>7s} {:>6s} {:>6s} {:>6s} {:>6s} {:>6s} "
           "{:>7s} {:>7s} {:>7s} {:>7s} {:>7s}")
    print(hdr.format("scenario", "TBApk", "TBAdel", "endo", "cord", "FCL",
                     "ARRI", "VAS", "ALT", "SB%", "PTB%", "MEC%", "NICU%",
                     "COMP%"))
    print("-" * 150)
    results = {}
    for name, ctl in S:
        rows, _ = simulate(ctl)
        e = endpoints(rows, ctl)
        results[name] = e
        print(hdr.format(
            name[:34], f"{e['TBA_peak']:.1f}", f"{e['TBA_del']:.1f}",
            f"{e['TBA_endo_del']:.1f}", f"{e['CORD']:.1f}", f"{e['FCL']:.1f}",
            f"{e['ARRI']:.3f}", f"{e['VAS']:.1f}", f"{e['ALT']:.0f}",
            f"{100*e['SB']:.3f}", f"{100*e['PTB_spont']:.1f}",
            f"{100*e['MEC']:.1f}", f"{100*e['NICU']:.1f}",
            f"{100*e['COMPOSITE']:.1f}"))
    print()
    print("onset / resolution / vitamin K:")
    for name, ctl in S:
        e = results[name]
        print("  {:<34s} onset GA {:>6s}  postpartum<10 {:>6s} d  INR {:.2f}"
              "  VITK {:.2f}  UDCAp {:.1f}  HCA {:.1f}".format(
                  name[:34],
                  "n/a" if e["onset_ga"] is None else f"{e['onset_ga']:.1f}",
                  "n/a" if e["resolve_d"] is None else f"{e['resolve_d']:.0f}",
                  e["INR"], e["VITK"], e["UDCAp"], e["HCA"]))
    return results


if __name__ == "__main__":
    main()
