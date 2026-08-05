#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gsd1a_reference_model.py
========================
Independent Python/scipy re-implementation of the GSD Ia (von Gierke disease)
QSP model that is shipped as `gsd1a_mrgsolve_model.R`.

WHY A SECOND IMPLEMENTATION EXISTS
-----------------------------------
The mrgsolve model is the deliverable; this file is the *check* on it.  Every
number quoted in README.md and in the .dot map is produced here, by
integration, not by assertion.  Writing the same 47 ODEs twice — once in C++
inside mrgsolve, once in Python — is the cheapest way to catch the class of
error that a QSP model is most prone to: an equation that is dimensionally
wrong, or a parameter that silently makes a clinically catastrophic drug look
beneficial.  Defects this file caught during construction are logged at the
bottom of README.md.

THE ORGANISING IDEA
-------------------
Glucose-6-phosphatase (G6PC1, 17q21.31) is the single terminal step shared by
BOTH routes of endogenous glucose production: glycogenolysis and
gluconeogenesis both converge on glucose-6-phosphate, and only G6Pase can
release free glucose from it.  Losing it therefore does not merely reduce
hepatic glucose output — it converts the liver from the body's glucose SOURCE
into an obligate glucose SINK, because glucokinase keeps running in the
forward direction with nothing to reverse it.

Everything else in GSD Ia is the overflow.  G6P that cannot leave as glucose
leaves as:
    glycogen   (G6P allosterically activates glycogen synthase)  -> hepatomegaly
    lactate    (glycolysis)                                      -> acidosis
    lipid      (DNL via acetyl-CoA / ChREBP)                     -> steatosis, hyperTG
    ribose-5-P (pentose phosphate pathway) -> PRPP -> purines    -> hyperuricaemia

So the model is written as ONE branch-point mass balance on hepatic G6P, and
the five clinical syndromes fall out of the branch fractions rather than being
five separately-parameterised submodels.

UNITS
-----
time            h
glucose         mmol/L (plasma), mmol/h (flux)
G6P, glycogen   umol/g liver
lactate         mmol/L
urate           mg/dL
triglyceride    mg/dL (plasma), mg/g (liver)
conversion      1 mg/kg/min glucose = 0.33333 mmol/kg/h
"""

import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# 0. STATE VECTOR
# ---------------------------------------------------------------------------
STATES = [
    # --- carbohydrate delivery (the "drug" in dietary therapy) ---------------
    "AST",       # 0  slow starch remaining, gut lumen        (mmol glucose-eq)
    "AGL",       # 1  free glucose available for absorption   (mmol)
    "ACOL",      # 2  starch escaping to colon (resistant)    (mmol glucose-eq)
    # --- systemic glucose ----------------------------------------------------
    "Gp",        # 3  plasma glucose                          (mmol/L)
    # --- hepatic carbon ------------------------------------------------------
    "G6P",       # 4  hepatic glucose-6-phosphate             (umol/g liver)
    "Glyc",      # 5  hepatic glycogen (glucosyl units)       (umol/g liver)
    "LV",        # 6  liver volume                            (mL)
    # --- lactate / redox / acid-base ----------------------------------------
    "Lac",       # 7  blood lactate                           (mmol/L)
    "Pyr",       # 8  blood pyruvate                          (mmol/L)
    "HCO3",      # 9  plasma bicarbonate                      (mmol/L)
    # --- lipid ---------------------------------------------------------------
    "MalCoA",    # 10 hepatic malonyl-CoA (CPT1 brake)        (rel. to normal)
    "TGliv",     # 11 hepatic triglyceride                    (mg/g liver)
    "TGpl",      # 12 plasma triglyceride                     (mg/dL)
    "FFA",       # 13 plasma free fatty acid                  (mmol/L)
    "KB",        # 14 3-hydroxybutyrate                       (mmol/L)
    # --- hormones ------------------------------------------------------------
    "Ins",       # 15 plasma insulin                          (pmol/L)
    "Gcg",       # 16 plasma glucagon                         (pmol/L)
    "Epi",       # 17 plasma epinephrine                      (nmol/L)
    "Cort",      # 18 plasma cortisol                         (nmol/L)
    "GH",        # 19 growth hormone tone                     (ug/L)
    "IGF1",      # 20 IGF-1                                   (nmol/L)
    # --- purine / urate ------------------------------------------------------
    "PRPP",      # 21 hepatic PRPP pool                       (rel. to normal)
    "UA",        # 22 serum urate                             (mg/dL)
    "ALLOg",     # 23 allopurinol, gut                        (mg)
    "ALLOc",     # 24 allopurinol, plasma                     (mg/L)
    "OXY",       # 25 oxypurinol, plasma (active)             (mg/L)
    # --- kidney --------------------------------------------------------------
    "GFRrel",    # 26 GFR relative to age-normal              (-)
    "Rdam",      # 27 accumulated nephron damage              (-)
    "UACR",      # 28 urine albumin:creatinine                (mg/g)
    "CitU",      # 29 urinary citrate                         (rel. to normal)
    "NCa",       # 30 nephrocalcinosis burden                 (-)
    # --- long-term liver -----------------------------------------------------
    "HCA",       # 31 hepatocellular adenoma burden           (rel. volume)
    "HCChaz",    # 32 cumulative HCC hazard                   (-)
    # --- growth / bone -------------------------------------------------------
    "HtSDS",     # 33 height standard-deviation score         (SDS)
    "BMDz",      # 34 lumbar BMD Z-score                      (Z)
    # --- symptom / exposure integrators -------------------------------------
    "TimeHypo",  # 35 cumulative time with Gp < 3.9           (h)
    "TimeNGC",   # 36 cumulative time with cerebral FAI < thr (h)
    "AUCdef",    # 37 cumulative cerebral fuel deficit        (mmol)
    # --- enzyme-replacement / gene therapy ----------------------------------
    "AAVvg",     # 38 transduced-hepatocyte vector load       (rel.)
    "MRNA",      # 39 LNP-mRNA derived G6Pase protein         (rel.)
    "G6Pact",    # 40 EXPRESSED fractional G6Pase activity    (frac of normal)
    "AntiAAV",   # 41 anti-capsid neutralising antibody       (rel.)
    # --- GSD Ib branch (SLC37A4 / G6PT) -------------------------------------
    "AG15",      # 42 plasma 1,5-anhydroglucitol              (ug/mL)
    "AG6Pn",     # 43 neutrophil 1,5-AG6P                     (rel.)
    "ANC",       # 44 absolute neutrophil count               (10^9/L)
    "EMPAc",     # 45 empagliflozin plasma                    (nmol/L)
    # --- adjunctive drugs ----------------------------------------------------
    "ACEIc",     # 46 ACE inhibitor plasma (lisinopril-like)  (ug/L)
]
IX = {s: i for i, s in enumerate(STATES)}
NST = len(STATES)


# ---------------------------------------------------------------------------
# 1. PARAMETERS
# ---------------------------------------------------------------------------
def base_params(age_y=1.0, bw=10.0, genotype="Ia", resid=0.0):
    """
    Baseline parameter set.

    age_y  chronological age (years) -- sets glucose demand per kg, brain
           fraction, liver fraction and the remaining liver growth that will
           dilute an AAV episome.
    bw     body weight (kg)
    resid  residual G6Pase activity from the patient's own genotype, as a
           fraction of normal.  True null (e.g. p.Arg83Cys/p.Gln347Ter) ~0;
           p.Gly188Arg and some Japanese/Chinese alleles retain a few %.
    """
    p = {}
    p["age"] = age_y
    p["BW"] = bw
    p["genotype"] = genotype        # "Ia" | "Ib" | "control"
    p["resid"] = resid

    # -- whole-body glucose demand -------------------------------------------
    # Glucose utilisation rate falls steeply with age because it is dominated
    # by brain, and brain:body mass falls.  Bier 1977 / Haymond: newborn 8,
    # infant 6-7, child 4-5, adolescent 3, adult 2.0-2.2 mg/kg/min.
    p["GUR_mgkgmin"] = 2.10 + 4.60 * math.exp(-age_y / 7.0)
    # brain share of that demand
    p["fCNS"] = 0.42 + 0.36 * math.exp(-age_y / 6.0)
    # glucose distribution volume
    p["Vg"] = 0.20 * bw                                    # L
    p["Vlac"] = 0.50 * bw                                  # L
    p["Vpl"] = 0.045 * bw                                  # L (plasma)
    p["Vua"] = 0.30 * bw                                   # L (urate space)
    # liver mass (g): 3.6 % of BW in infancy -> 2.4 % in adults
    p["fLW"] = 0.024 + 0.014 * math.exp(-age_y / 4.0)
    p["LW"] = p["fLW"] * bw * 1000.0

    GUR = p["GUR_mgkgmin"] * bw / 3.0                      # mmol/h, total
    p["GUR_tot"] = GUR
    p["Dcns"] = p["fCNS"] * GUR                            # mmol glucose-eq/h
    p["Uper0"] = 0.30 * (1 - p["fCNS"]) * GUR              # insulin-independent
    p["Si"] = 0.70 * (1 - p["fCNS"]) * GUR / 60.0          # per pmol/L insulin

    # -- cerebral fuel --------------------------------------------------------
    # Lactate is a real cerebral fuel.  van Hall 2009 / Boumezbeur 2010: at
    # plasma lactate ~7 mM lactate supplies ~55-60 % of brain oxidative
    # metabolism.  This is not a detail -- it is why a GSD Ia child can be at
    # 1.8 mmol/L glucose and asymptomatic, and it is modelled explicitly.
    p["VlacCNS"] = 0.90                                    # x Dcns, max share
    p["KlacCNS"] = 4.0                                     # mmol/L
    p["Kb"] = 0.70                                         # mmol/L, GLUT1 app.
    p["FAIthr"] = 0.830                                    # neuroglycopenia

    # -- gut / starch ---------------------------------------------------------
    p["kdis"] = 0.45          # 1/h  uncooked cornstarch amylolysis
    p["Fabs"] = 0.75          # fraction reaching portal vein as glucose
    p["ka"] = 3.0             # 1/h  absorption of luminal free glucose
    p["kcol"] = 0.20          # 1/h  colonic disposal (SCFA), not glucose
    p["Rdrip"] = 0.0          # mmol/h mean enteral glucose delivery
    # Chronic runs cannot afford to integrate 30 years of q4h boluses, so the
    # day-cycle is represented by its mean and its AMPLITUDE.  That distinction
    # is not cosmetic: two regimens with the SAME daily carbohydrate but
    # different pulsatility are different diseases, because the overnight
    # trough is where the glycogen unloads into lactate.
    p["RdripAmp"] = 0.0       # 0 = continuous drip, ~0.9 = 3 large meals
    p["RdripPer"] = 24.0      # h

    # -- hepatic glucose uptake (glucokinase) --------------------------------
    # Never reverses in GSD Ia; this is the futile-cycle term that makes the
    # liver a net consumer of glucose at all times.  Calibrated to ~1.5-2.0
    # mg/kg/min hepatic glucose uptake at euglycaemia with insulin present.
    p["Vgk"] = 0.28 * GUR
    p["Kgk"] = 8.0
    p["ngk"] = 1.7
    p["gkIns"] = 0.60         # fractional GK activation by insulin (GKRP)

    # -- G6Pase --------------------------------------------------------------
    # Vmax is deliberately well ABOVE basal EGP: a normal liver runs G6Pase at
    # a small fraction of capacity.  That headroom is what makes the restored-
    # activity dose-response saturate so steeply in gene therapy.
    p["VmaxG6P"] = 6.2 * GUR                               # mmol/h at a = 1
    p["KmG6P"] = 1.5                                       # umol/g
    # extrahepatic (renal + intestinal) G6Pase, also absent in Ia
    p["fEGPextra"] = 0.10

    # -- glycogen -------------------------------------------------------------
    # kgp is calibrated so that a fully activated liver runs glycogenolysis at
    # ~45-50 % of total glucose demand, i.e. a normal infant liver empties its
    # glycogen in ~10 h -- which is the observed fasting tolerance.
    p["Glyc0"] = 310.0        # umol glucosyl/g  (~50 mg/g, normal)
    p["kgp"] = 0.048          # 1/h phosphorylase, x activation
    p["Vgs"] = 1.30 * GUR     # mmol/h
    p["Kgs"] = 1.10           # umol/g, Hill 2 (G6P activates GS)
    # Residual EGP calibration target: isotopic studies in GSD I children put
    # endogenous glucose production at ~1.0-1.5 mg/kg/min (vs 5-6 normal).  The
    # ONLY two routes that bypass G6Pase are the alpha-1,6 branch-point glucose
    # released by the debranching enzyme and lysosomal acid-maltase.
    p["fdeb"] = 0.100         # free glucose from alpha-1,6 debranching
    p["klys"] = 5.0e-3        # 1/h lysosomal (GAA) glycogen -> free glucose
    p["GlycMax"] = 950.0      # umol/g, physical ceiling

    # -- glycolysis / PPP -----------------------------------------------------
    p["Vpfk"] = 0.95 * GUR
    p["Kpfk"] = 3.0           # umol/g
    p["Vppp"] = 0.10 * GUR
    p["Kppp"] = 3.0

    # -- gluconeogenesis ------------------------------------------------------
    p["Vgng"] = 0.55 * GUR
    p["Kgng"] = 2.0           # mmol/L lactate
    p["Sub0"] = 1.2           # mmol/L alanine + glycerol equivalent supply
    p["fGngLac"] = 0.55       # fraction of gluconeogenic carbon from lactate

    # -- lactate --------------------------------------------------------------
    p["Vlacox"] = 0.75 * GUR  # mmol/h peripheral (non-CNS) oxidation
    p["Klacox"] = 3.5
    p["klacren"] = 0.55 * bw / 10.0
    p["LacThrRen"] = 5.0
    p["LacProdPer"] = 0.50 * GUR
    p["fLacOut"] = 0.80       # glycolytic pyruvate leaving as lactate vs PDH
    p["KacidPFK"] = 9.0       # mmol/L HCO3, acidotic inhibition of PFK-1
    p["LPratio"] = 14.0       # near-equilibrium LDH; NORMAL in GSD I

    # -- acid-base ------------------------------------------------------------
    p["HCO3n"] = 24.0
    p["kbuf"] = 1.2
    p["cLac_HCO3"] = 1.25

    # -- lipid ----------------------------------------------------------------
    p["kdnl"] = 0.55
    p["kmc"] = 1.6            # 1/h malonyl-CoA turnover
    p["Kmc"] = 1.0            # CPT1 inhibition constant (beta-oxidation)
    p["KmcKet"] = 0.60        # CPT1 inhibition constant (ketogenesis)
    p["Vvldl"] = 0.625        # mg/g/h  ApoB-limited VLDL-TG export Vmax
    p["KvldlTG"] = 12.0       # mg/g,   half-maximal hepatic TG for export
    p["Vlpl"] = 110.0         # mg/dL/h saturable (LPL) clearance
    p["Klpl"] = 500.0
    p["klinTG"] = 0.035       # 1/h non-saturable remnant/spillover clearance
    p["fLPLdef"] = 0.52       # LPL activity in GSD Ia (Forget 1974, Levy 1988)
    p["kbox"] = 0.030
    p["kffaTG"] = 0.22        # hepatic FFA re-esterification -> TG
    p["klip"] = 0.55          # lipolysis rate scaler
    p["kffau"] = 1.9
    p["kket"] = 3.0
    p["kkbu"] = 1.3

    # -- hormones -------------------------------------------------------------
    p["Imax"] = 420.0
    p["Kins"] = 7.0
    p["hins"] = 2.4
    p["kIns"] = 4.2
    p["Gcgmax"] = 42.0
    p["Kgcg"] = 4.3
    p["hgcg"] = 3.0
    p["kGcg"] = 5.0
    p["Epimax"] = 3.2
    p["Kepi"] = 3.2
    p["hepi"] = 6.0
    p["kEpi"] = 8.0
    p["Cortmax"] = 620.0
    p["Kcort"] = 3.1
    p["kCort"] = 0.45
    p["GHmax"] = 9.0
    p["Kgh"] = 3.4
    p["kGH"] = 0.9
    p["IGF1n"] = 25.0
    p["kIGF"] = 0.03

    # -- purine / urate -------------------------------------------------------
    p["UAprod0"] = 2.6 * bw / 10.0         # mg/h baseline
    p["kprpp"] = 2.2
    p["kprppUse"] = 2.2
    p["nPRPP"] = 0.45                      # compressive: pool << flux ratio
    p["JpppRef"] = None                    # set below, = control-liver PPP flux
    p["kATPdeg"] = 0.15
    p["G6Pth_ATP"] = 1.2                   # umol/g, Pi-trapping threshold
    p["UAcl0"] = 0.115 * bw / 10.0         # L/h effective urate clearance
    p["KiLacUA"] = 7.0                     # mmol/L, URAT1 exchange competition
    p["KiKBUA"] = 3.0
    p["EmaxAllo"] = 0.62
    p["IC50oxy"] = 5.5                     # mg/L oxypurinol
    p["kaAllo"] = 2.2
    p["kAllo"] = 0.55                      # allopurinol elimination 1/h
    p["fOxy"] = 0.76                       # fraction converted to oxypurinol
    p["kOxy"] = 0.031                      # 1/h  (t1/2 ~22 h)
    p["VdAllo"] = 1.6 * bw / 10.0

    # -- kidney ---------------------------------------------------------------
    p["khf"] = 0.0016
    p["kdam"] = 8.5e-7
    p["krep"] = 4.0e-6
    p["wLac"] = 0.20
    p["wHF"] = 0.85
    p["wUA"] = 0.045
    p["wTG"] = 0.00035
    p["kUACR"] = 0.004
    p["UACRmax"] = 1800.0
    p["EmaxACEI"] = 0.55
    p["IC50ACEI"] = 12.0                   # ug/L
    p["kaACEI"] = 1.1
    p["kACEI"] = 0.058                     # t1/2 ~12 h
    p["kCit"] = 0.05
    p["kNCa"] = 2.0e-5

    # -- long-term liver ------------------------------------------------------
    p["kHCAini"] = 3.1e-7
    p["kHCAgro"] = 1.4e-5
    p["HCAmax"] = 0.35                     # adenoma fraction of liver volume
    p["nMCI"] = 2.6                        # supralinearity of control index
    p["kHCC"] = 2.4e-6

    # -- growth / bone --------------------------------------------------------
    # Written as relaxation towards a control-dependent TARGET, not as a free
    # integral.  A height SDS that integrates a constant negative rate for 30
    # years reaches -8, which is not a number a human being can have; and
    # linear growth stops at skeletal maturity whatever the metabolic state.
    p["kgro"] = 6.0e-5                     # 1/h relaxation of height SDS
    p["HtSDStar0"] = -0.40                 # attainable SDS on perfect control
    p["HtSDSmci"] = 2.60                   # SDS lost at maximal MCI
    p["ageMature"] = 18.0                  # years
    p["kbmd"] = 4.0e-6                     # 1/h relaxation of BMD Z
    p["BMDztar0"] = 0.20
    p["BMDzAcid"] = 2.60
    p["BMDzMci"] = 1.60

    # -- gene / mRNA therapy --------------------------------------------------
    p["epsAAV"] = 1.0         # fractional activity per unit transduced load
    p["kexpr"] = 0.020        # 1/h expression onset (t1/2 ~35 h)
    p["ksil"] = 6.0e-6        # 1/h slow episome loss beyond growth dilution
    p["kmrna"] = 0.0072       # 1/h  LNP-mRNA protein decay (t1/2 ~4 d)
    p["epsMRNA"] = 1.0
    p["kAbUp"] = 0.02
    p["kAb"] = 2.0e-5
    # remaining liver growth from current age to adult, as a multiple
    p["LWadult"] = 0.024 * 70.0 * 1000.0
    p["tauGrow"] = 5.5 * 8766.0            # h, liver growth time constant

    # -- GSD Ib ---------------------------------------------------------------
    # 1,5-AG is a CONCENTRATION, so its input term must not scale with body
    # weight the way an amount would.  Normal plasma 1,5-AG is 10-25 ug/mL.
    p["AG15in"] = 0.85                     # ug/mL/h  dietary 1,5-AG input
    p["kAG15"] = 0.055                     # 1/h renal handling (SGLT4/SLC5A9)
    p["kAGup"] = 0.030                     # neutrophil uptake -> 1,5-AG6P
    p["kAG6Pout"] = 0.200                  # G6PT-dependent efflux (absent in Ib)
    p["kAG6Pdeg"] = 0.020                  # slow G6PT-independent loss
    p["fG6PT"] = 1.0                       # 1 = Ia/normal, 0 = Ib null
    p["ANCprod"] = 3.9                     # 10^9/L/h scaled
    p["KAG6P"] = 1.0
    p["hAG6P"] = 3.2
    p["kANC"] = 1.0
    p["EMPAemax"] = 3.4                    # fold rise in 1,5-AG clearance
    p["EMPAec50"] = 60.0                   # nmol/L
    p["kaEMPA"] = 1.4
    p["kEMPA"] = 0.055
    p["GCSF"] = 0.0

    # -- therapy switches -----------------------------------------------------
    # reference PPP flux of a NORMAL liver (G6P ~ 0.25 umol/g) -- the PRPP
    # pool is expressed relative to this
    p["JpppRef"] = p["Vppp"] * (0.25 / (p["Kppp"] + 0.25))

    p["fibrate"] = 0.0
    p["statin"] = 0.0
    p["citrate"] = 0.0
    p["illness"] = 0.0        # 0-1 catabolic stress / vomiting
    p["chronic_scale"] = 1.0  # accelerate slow states in long-horizon runs
    return p


def set_genotype(p, g):
    """Ia = G6PC1 null; Ib = SLC37A4 (transporter) null; control = healthy."""
    p["genotype"] = g
    if g == "control":
        p["resid"] = 1.0
        p["fG6PT"] = 1.0
        p["fLPLdef"] = 1.0
    elif g == "Ia":
        p["fG6PT"] = 1.0
    elif g == "Ib":
        # The transporter defect blocks G6P ENTRY to the ER, so the catalytic
        # subunit is intact but starved of substrate: same metabolic phenotype,
        # usually a shade milder, PLUS the neutrophil lesion.
        p["resid"] = max(p["resid"], 0.03)
        p["fG6PT"] = 0.0
    return p


# ---------------------------------------------------------------------------
# 2. ALGEBRAIC HELPERS
# ---------------------------------------------------------------------------
def hill(x, K, n=1.0):
    x = max(x, 0.0)
    return x ** n / (K ** n + x ** n) if (x > 0 or K > 0) else 0.0


def cerebral_fuel(Gp, Lac, p):
    """
    Returns (lactate-derived supply, glucose need, glucose actually taken,
             fuel adequacy index).  All in glucose-equivalents mmol/h.
    """
    D = p["Dcns"]
    sup_lac = p["VlacCNS"] * D * hill(Lac, p["KlacCNS"])
    need_glc = max(D - sup_lac, 0.0)
    sat = hill(Gp, p["Kb"])
    got_glc = need_glc * sat
    FAI = (sup_lac + got_glc) / D if D > 0 else 1.0
    return sup_lac, need_glc, got_glc, FAI


def isofuel_glucose(Lac, p, FAI_target=None):
    """
    Invert the cerebral fuel balance: at blood lactate `Lac`, what plasma
    glucose delivers the same cerebral ATP flux as `FAI_target`?
    This is the quantitative statement of "GSD Ia patients tolerate glucose
    concentrations that would seize a normal child".
    """
    if FAI_target is None:
        FAI_target = p["FAIthr"]
    D = 1.0
    fl = p["VlacCNS"] * hill(Lac, p["KlacCNS"])
    need = D - fl
    if need <= 1e-9:
        return 0.0                      # lactate alone covers the brain
    sat = (FAI_target - fl) / need
    if sat <= 0:
        return 0.0
    if sat >= 1:
        return float("inf")
    return sat * p["Kb"] / (1.0 - sat)


def liver_growth_rate(t, p):
    """Fractional liver growth per hour at the patient's current age."""
    age_h = p["age"] * 8766.0 + t
    LWt = p["LWadult"] - (p["LWadult"] - p["LW"]) * math.exp(
        -(age_h - p["age"] * 8766.0) / p["tauGrow"])
    dLW = (p["LWadult"] - LWt) / p["tauGrow"]
    return max(dLW / max(LWt, 1.0), 0.0)


# ---------------------------------------------------------------------------
# 3. RIGHT-HAND SIDE
# ---------------------------------------------------------------------------
def rhs(t, y, p):
    d = np.zeros(NST)
    Y = {s: y[i] for i, s in enumerate(STATES)}
    pos = lambda x: max(x, 0.0)

    AST, AGL, ACOL = pos(Y["AST"]), pos(Y["AGL"]), pos(Y["ACOL"])
    Gp = pos(Y["Gp"])
    G6P, Glyc = pos(Y["G6P"]), pos(Y["Glyc"])
    LV = pos(Y["LV"])
    Lac, Pyr, HCO3 = pos(Y["Lac"]), pos(Y["Pyr"]), pos(Y["HCO3"])
    MalCoA, TGliv, TGpl = pos(Y["MalCoA"]), pos(Y["TGliv"]), pos(Y["TGpl"])
    FFA, KB = pos(Y["FFA"]), pos(Y["KB"])
    Ins, Gcg, Epi, Cort, GH, IGF1 = (pos(Y["Ins"]), pos(Y["Gcg"]), pos(Y["Epi"]),
                                     pos(Y["Cort"]), pos(Y["GH"]), pos(Y["IGF1"]))
    PRPP, UA = pos(Y["PRPP"]), pos(Y["UA"])
    ALLOg, ALLOc, OXY = pos(Y["ALLOg"]), pos(Y["ALLOc"]), pos(Y["OXY"])
    GFRrel, Rdam, UACR, CitU, NCa = (pos(Y["GFRrel"]), pos(Y["Rdam"]),
                                     pos(Y["UACR"]), pos(Y["CitU"]), pos(Y["NCa"]))
    HCA, HCChaz = pos(Y["HCA"]), pos(Y["HCChaz"])
    TimeHypo, TimeNGC, AUCdef = Y["TimeHypo"], Y["TimeNGC"], Y["AUCdef"]
    AAVvg, MRNA, G6Pact, AntiAAV = (pos(Y["AAVvg"]), pos(Y["MRNA"]),
                                    pos(Y["G6Pact"]), pos(Y["AntiAAV"]))
    AG15, AG6Pn, ANC, EMPAc = (pos(Y["AG15"]), pos(Y["AG6Pn"]),
                               pos(Y["ANC"]), pos(Y["EMPAc"]))
    ACEIc = pos(Y["ACEIc"])

    LWg = max(LV * 1.05, 1.0)                 # g liver, from volume
    cs = p["chronic_scale"]

    # ---------------- hormones (fast, glucose-driven) ----------------------
    Iss = p["Imax"] * hill(Gp, p["Kins"], p["hins"])
    d[IX["Ins"]] = p["kIns"] * (Iss - Ins)
    Gss = p["Gcgmax"] * (1.0 - hill(Gp, p["Kgcg"], p["hgcg"])) + 4.0
    d[IX["Gcg"]] = p["kGcg"] * (Gss - Gcg)
    Ess = p["Epimax"] * (1.0 - hill(Gp, p["Kepi"], p["hepi"])) + 0.15 + 0.9 * p["illness"]
    d[IX["Epi"]] = p["kEpi"] * (Ess - Epi)
    Css = p["Cortmax"] * (1.0 - hill(Gp, p["Kcort"], 4.0)) + 180.0 + 220.0 * p["illness"]
    d[IX["Cort"]] = p["kCort"] * (Css - Cort)
    GHss = p["GHmax"] * (1.0 - hill(Gp, p["Kgh"], 4.0)) + 1.2
    d[IX["GH"]] = p["kGH"] * (GHss - GH)

    Insrel = Ins / 120.0
    Gcgrel = Gcg / 12.0
    Epirel = Epi / 0.5

    # ---------------- gut delivery -----------------------------------------
    Jdis = p["kdis"] * AST
    d[IX["AST"]] = -Jdis
    Ra_gut = p["ka"] * AGL
    Rd = p["Rdrip"] * (1.0 + p["RdripAmp"]
                       * math.sin(2.0 * math.pi * t / p["RdripPer"]))
    d[IX["AGL"]] = Jdis * p["Fabs"] + max(Rd, 0.0) - Ra_gut
    d[IX["ACOL"]] = Jdis * (1.0 - p["Fabs"]) - p["kcol"] * ACOL

    # ---------------- hepatic branch point ---------------------------------
    # (a) glucokinase: glucose IN, always forward
    Jgk = (p["Vgk"] * (1.0 + p["gkIns"] * min(Insrel, 2.5))
           * hill(Gp, p["Kgk"], p["ngk"]) * (LWg / p["LW"]))

    # (b) glycogen phosphorylase: activated by glucagon/epi, braked by glucose
    # Hill constants deliberately set ABOVE the resting hormone levels: with a
    # half-saturating K at the resting concentration the enzyme sits near its
    # ceiling and a glucagon challenge does almost nothing, which is not what
    # a glucagon test does to a normal liver.
    phos_act = (0.10 + 2.2 * hill(Gcgrel, 2.5, 2.0) + 0.90 * hill(Epirel, 2.5, 2.0)) \
        / (1.0 + 0.9 * hill(Gp, 6.0, 3.0)) / (1.0 + 0.6 * min(Insrel, 3.0))
    # umol/g x g / 1000 = mmol of glucosyl units in the store; x 1/h = mmol/h
    Jgp = p["kgp"] * Glyc * phos_act * LWg / 1000.0             # mmol/h
    # of which the alpha-1,6 branch-point residues emerge as FREE glucose --
    # this is the ONLY glycogen-derived glucose a GSD Ia liver can release,
    # and it is why residual EGP is small but not zero.
    Jdeb = p["fdeb"] * Jgp
    Jgp_g6p = Jgp - Jdeb
    Jlys = p["klys"] * Glyc * LWg / 1000.0                      # free glucose

    # (c) gluconeogenesis: still runs, still dumps into G6P
    gng_horm = 0.45 + 0.75 * hill(Gcgrel, 1.0, 1.4) + 0.45 * (Cort / 400.0)
    Jgng = p["Vgng"] * gng_horm * hill(Lac + 0.6 * FFA + p["Sub0"], p["Kgng"])

    # (d) G6Pase: the lesion.  a(t) = germline residual + expressed transgene,
    #     and in Ib the catalytic subunit is intact but the TRANSPORTER is not,
    #     so substrate never reaches it.
    a_enz = min(p["resid"] + G6Pact, 1.0)
    transporter = p["fG6PT"] if p["genotype"] != "control" else 1.0
    a_eff = a_enz * (transporter if p["genotype"] == "Ib" else 1.0)
    Jg6pase = p["VmaxG6P"] * a_eff * hill(G6P, p["KmG6P"]) * (LWg / p["LW"])
    Jextra = p["fEGPextra"] * p["VmaxG6P"] * a_eff * hill(G6P, p["KmG6P"])

    # (e) glycogen synthase: G6P is its ALLOSTERIC ACTIVATOR (Hill 2).  This is
    #     the single equation that turns a glucose-release defect into a
    #     storage disease.
    # Two arms, multiplied: covalent (insulin dephosphorylates GS, glucagon
    # phosphorylates it) and allosteric (G6P binds GS-b and switches it on
    # independently of phosphorylation state).  A normal liver stores glycogen
    # through the covalent arm; a GSD Ia liver stores it through the allosteric
    # arm even when the covalent arm says "stop".  That is the whole reason
    # this is a STORAGE disease and not merely a hypoglycaemia syndrome.
    gs_cov = (0.15 + 1.6 * hill(Ins, 110.0, 1.5)) / (1.0 + 0.8 * hill(Gcgrel, 1.2, 2.0))
    gs_allo = 0.22 + 0.78 * hill(G6P, p["Kgs"], 2.0)
    gs_act = gs_cov * gs_allo
    Jgs = (p["Vgs"] * gs_act * (1.0 - hill(Glyc, p["GlycMax"], 6.0))
           * (LWg / p["LW"]))

    # (f) glycolysis and PPP: the two overflow drains
    f26 = 1.0 + 0.55 * min(Insrel, 2.0) - 0.20 * hill(Gcgrel, 1.0, 1.5)
    # PFK-1 is strongly inhibited by intracellular acidosis.  This is the only
    # negative feedback that terminates a GSD Ia lactic crisis, and without it
    # the model has no brake at all -- lactate simply runs to the ceiling.
    acid_brake = hill(HCO3, p["KacidPFK"], 3.0) / hill(24.0, p["KacidPFK"], 3.0)
    Jgly = (p["Vpfk"] * hill(G6P, p["Kpfk"]) * max(f26, 0.25)
            * acid_brake * (LWg / p["LW"]))
    Jppp = p["Vppp"] * hill(G6P, p["Kppp"]) * (LWg / p["LW"])

    dG6P = (Jgk + Jgp_g6p + Jgng - Jg6pase - Jgs - Jgly - Jppp)
    d[IX["G6P"]] = dG6P * 1000.0 / LWg
    d[IX["Glyc"]] = (Jgs - Jgp - Jlys) * 1000.0 / LWg

    # ---------------- systemic glucose --------------------------------------
    sup_lac, need_glc, got_glc, FAI = cerebral_fuel(Gp, Lac, p)
    Uper = (p["Uper0"] + p["Si"] * Ins) * hill(Gp, 2.5) * (1.0 - 0.25 * p["illness"])
    # SGLT2 inhibition lowers the renal glucose threshold (Ib empagliflozin arm)
    Gthr = 10.0 - 3.4 * hill(EMPAc, p["EMPAec50"])
    Uren = 0.9 * p["BW"] / 10.0 * max(Gp - Gthr, 0.0) * GFRrel

    EGP = Jg6pase + Jextra + Jdeb + Jlys
    d[IX["Gp"]] = (Ra_gut + EGP - got_glc - Uper - Jgk - Uren) / p["Vg"]

    # ---------------- lactate / pyruvate / acid-base ------------------------
    Jlacox = p["Vlacox"] * hill(Lac, p["Klacox"]) * (1.0 + 0.25 * Epirel / 3.0)
    Jlacren = p["klacren"] * max(Lac - p["LacThrRen"], 0.0) * GFRrel
    LacProdPer = p["LacProdPer"] * (1.0 + 0.5 * hill(Epirel, 2.0)) \
        + 0.35 * p["LacProdPer"] * p["illness"]
    JlacCNS = sup_lac * 2.0                    # glucose-eq -> lactate mol
    d[IX["Lac"]] = (2.0 * Jgly * p["fLacOut"] + LacProdPer - Jlacox - Jlacren
                    - JlacCNS - 2.0 * Jgng * p["fGngLac"]) / p["Vlac"]
    d[IX["Pyr"]] = (Lac / p["LPratio"] - Pyr) * 6.0
    HCO3ss = p["HCO3n"] - p["cLac_HCO3"] * max(Lac - 1.0, 0.0) + 4.0 * p["citrate"]
    d[IX["HCO3"]] = p["kbuf"] * (max(HCO3ss, 6.0) - HCO3)

    # ---------------- lipid --------------------------------------------------
    chrebp = hill(G6P, 2.0, 1.8)
    # Carbon accounting: glycolytic pyruvate splits between LDH (-> lactate,
    # exported) and PDH (-> acetyl-CoA -> de novo lipogenesis / oxidation).
    # Counting the same pyruvate into both pools would inflate lactate and
    # triglyceride simultaneously, which is exactly the sort of error a second
    # implementation exists to catch.
    Jdnl = (p["kdnl"] * 2.0 * Jgly * (1.0 - p["fLacOut"])
            * (0.35 + 0.9 * chrebp) * (1.0 + 0.5 * min(Insrel, 2.0)))
    d[IX["MalCoA"]] = p["kmc"] * ((0.25 + 2.1 * chrebp
                                   + 0.5 * min(Insrel, 2.0)) - MalCoA)
    lipo = p["klip"] * (1.0 + 1.4 * hill(Epirel, 2.0) + 0.5 * hill(Gcgrel, 1.5)) \
        / (1.0 + 1.8 * min(Insrel, 3.0))
    d[IX["FFA"]] = lipo - p["kffau"] * FFA
    # ketogenesis is throttled by malonyl-CoA -> GSD Ia is HYPOketotic, unlike
    # GSD 0/III/VI which are ketotic.  The model must reproduce that unprompted.
    d[IX["KB"]] = p["kket"] * FFA / (1.0 + MalCoA / p["KmcKet"]) - p["kkbu"] * KB

    # VLDL-TG export saturates: ApoB-100 secretion capacity, not TG supply, is
    # rate-limiting once the hepatocyte is loaded.  That is why the liver
    # steatoses rather than exporting its way out of trouble.
    Jvldl = p["Vvldl"] * hill(TGliv, p["KvldlTG"], 1.0)          # mg/g/h
    Jbox = p["kbox"] * TGliv / (1.0 + MalCoA / p["Kmc"])
    d[IX["TGliv"]] = (Jdnl * 180.0 * 0.42 / max(LWg, 1.0)
                      + p["kffaTG"] * FFA * 12.0 - Jvldl - Jbox)
    lpl = p["fLPLdef"] * (1.0 + 0.55 * min(Insrel, 2.0)) \
        * (1.0 + 1.35 * p["fibrate"])
    Jlpl = (p["Vlpl"] * hill(TGpl, p["Klpl"]) + p["klinTG"] * TGpl) * lpl
    statin_eff = 1.0 - 0.22 * p["statin"]
    # mg/h delivered into a plasma pool of TGpl[mg/dL] x Vpl[L] x 10[dL/L]
    d[IX["TGpl"]] = Jvldl * LWg / (p["Vpl"] * 10.0) * statin_eff - Jlpl

    # ---------------- liver volume ------------------------------------------
    LVtar = (p["fLW"] * p["BW"] * 1000.0 / 1.05) * (
        1.0 + 0.62 * max(Glyc / p["Glyc0"] - 1.0, 0.0)
        + 0.30 * max(TGliv / 25.0 - 1.0, 0.0)) * (1.0 + 0.8 * HCA)
    d[IX["LV"]] = 0.0025 * (LVtar - LV)

    # ---------------- purine / urate ----------------------------------------
    # PRPP tracks PPP flux compressively (pool is buffered by consumption)
    prpp_tar = (Jppp / max(p["JpppRef"], 1e-9)) ** p["nPRPP"]
    d[IX["PRPP"]] = p["kprpp"] * prpp_tar - p["kprppUse"] * PRPP
    ATPdeg = p["kATPdeg"] * max(G6P - p["G6Pth_ATP"], 0.0)
    Eallo = p["EmaxAllo"] * hill(OXY, p["IC50oxy"])
    UAprod = (p["UAprod0"] * (0.45 + 0.75 * PRPP + ATPdeg)) * (1.0 - Eallo)
    UAcl = p["UAcl0"] * GFRrel / (1.0 + Lac / p["KiLacUA"] + KB / p["KiKBUA"])
    d[IX["UA"]] = (UAprod - UAcl * UA * 10.0) / p["Vua"] * 0.1
    d[IX["ALLOg"]] = -p["kaAllo"] * ALLOg
    d[IX["ALLOc"]] = p["kaAllo"] * ALLOg / p["VdAllo"] - p["kAllo"] * ALLOc
    d[IX["OXY"]] = p["fOxy"] * p["kAllo"] * ALLOc - p["kOxy"] * OXY

    # ---------------- kidney -------------------------------------------------
    Eacei = p["EmaxACEI"] * hill(ACEIc, p["IC50ACEI"])
    GFRtar = 1.0 + 0.30 * hill(max(Lac - 2.0, 0.0), 4.0) - 0.9 * Rdam
    d[IX["GFRrel"]] = cs * p["khf"] * (GFRtar - GFRrel)
    drive = (p["wLac"] * max(Lac - 2.0, 0.0)
             + p["wHF"] * max(GFRrel - 1.0, 0.0)
             + p["wUA"] * max(UA - 5.5, 0.0)
             + p["wTG"] * TGpl)
    d[IX["Rdam"]] = cs * (p["kdam"] * drive * (1.0 - Eacei) - p["krep"] * Rdam)
    d[IX["UACR"]] = cs * (p["kUACR"] * (Rdam * 2600.0 * (1.0 - Eacei) - UACR))
    d[IX["CitU"]] = cs * p["kCit"] * ((HCO3 / p["HCO3n"]) ** 2.2
                                      * (1.0 + 0.6 * p["citrate"]) - CitU)
    d[IX["NCa"]] = cs * p["kNCa"] * max(1.0 - CitU, 0.0) * (1.0 + 0.4 * max(UA - 6, 0))

    # ---------------- long-term liver ---------------------------------------
    # Metabolic control index: 0 = perfect, 1 = severely uncontrolled.
    MCI = min(1.0, (0.34 * max(Lac - 2.2, 0.0) / 5.0
                    + 0.26 * max(TGpl - 250.0, 0.0) / 900.0
                    + 0.22 * max(UA - 5.5, 0.0) / 5.0
                    + 0.18 * max(3.9 - Gp, 0.0) / 2.0) * 1.9)
    # logistic, not exponential: adenoma burden is a fraction of liver volume
    # and cannot exceed it
    d[IX["HCA"]] = cs * (p["kHCAini"]
                         + p["kHCAgro"] * HCA * max(1.0 - HCA / p["HCAmax"], 0.0)
                         ) * (MCI ** p["nMCI"]) * 4.0
    d[IX["HCChaz"]] = cs * p["kHCC"] * HCA * (1.0 + 2.0 * MCI)

    # ---------------- growth / bone -----------------------------------------
    d[IX["IGF1"]] = p["kIGF"] * (p["IGF1n"] * (1.0 - 0.55 * MCI)
                                 * (0.55 + 0.45 * hill(Gp, 3.5, 3.0)) - IGF1)
    acid = max(p["HCO3n"] - HCO3, 0.0) / 8.0
    age_now = p["age"] + t / 8766.0
    growing = 1.0 / (1.0 + math.exp((age_now - p["ageMature"]) / 1.5))
    HtTar = p["HtSDStar0"] - p["HtSDSmci"] * MCI
    d[IX["HtSDS"]] = cs * p["kgro"] * growing * (HtTar - Y["HtSDS"])
    BMDzTar = p["BMDztar0"] - p["BMDzAcid"] * acid - p["BMDzMci"] * MCI
    d[IX["BMDz"]] = cs * p["kbmd"] * (BMDzTar - Y["BMDz"])

    # ---------------- symptom integrators -----------------------------------
    # smoothed indicators -- a hard if/else here makes the RHS discontinuous
    # and LSODA will chatter across the switch
    d[IX["TimeHypo"]] = 1.0 / (1.0 + math.exp(6.0 * (Gp - 3.9)))
    d[IX["TimeNGC"]] = 1.0 / (1.0 + math.exp(120.0 * (FAI - p["FAIthr"])))
    d[IX["AUCdef"]] = max(p["Dcns"] * p["FAIthr"] - (sup_lac + got_glc), 0.0)

    # ---------------- gene / mRNA therapy ------------------------------------
    # THE durability term: AAV episomes do not replicate, so a growing liver
    # dilutes them.  In a 2-year-old this dominates every other loss process.
    gr = liver_growth_rate(t, p)
    d[IX["AAVvg"]] = -(gr + p["ksil"]) * AAVvg * cs
    d[IX["MRNA"]] = -p["kmrna"] * MRNA
    Ass = p["epsAAV"] * AAVvg + p["epsMRNA"] * MRNA
    d[IX["G6Pact"]] = p["kexpr"] * (Ass - G6Pact)
    d[IX["AntiAAV"]] = -p["kAb"] * AntiAAV * cs

    # ---------------- GSD Ib branch ------------------------------------------
    empa = 1.0 + (p["EMPAemax"] - 1.0) * hill(EMPAc, p["EMPAec50"])
    d[IX["AG15"]] = p["AG15in"] - p["kAG15"] * empa * AG15 * GFRrel
    # 1,5-AG is phosphorylated by hexokinase inside neutrophils; only G6PT can
    # get 1,5-AG6P back out.  In Ib it therefore accumulates and inhibits
    # hexokinase -> the neutrophil starves.
    d[IX["AG6Pn"]] = (p["kAGup"] * AG15 / 14.0
                      - p["kAG6Pout"] * p["fG6PT"] * AG6Pn
                      - p["kAG6Pdeg"] * AG6Pn)
    ancprod = p["ANCprod"] * (1.0 - hill(AG6Pn, p["KAG6P"], p["hAG6P"])) \
        * (1.0 + 1.8 * p["GCSF"])
    d[IX["ANC"]] = ancprod - p["kANC"] * ANC
    d[IX["EMPAc"]] = -p["kEMPA"] * EMPAc

    d[IX["ACEIc"]] = -p["kACEI"] * ACEIc
    return d


# ---------------------------------------------------------------------------
# 4. INITIAL CONDITIONS
# ---------------------------------------------------------------------------
def y0_for(p):
    y = np.zeros(NST)
    ctrl = p["genotype"] == "control"
    y[IX["Gp"]] = 4.8 if ctrl else 3.6
    y[IX["G6P"]] = 0.22 if ctrl else 3.6
    y[IX["Glyc"]] = p["Glyc0"] * (1.0 if ctrl else 2.1)
    y[IX["LV"]] = p["fLW"] * p["BW"] * 1000.0 / 1.05 * (1.0 if ctrl else 2.0)
    y[IX["Lac"]] = 1.0 if ctrl else 5.0
    y[IX["Pyr"]] = y[IX["Lac"]] / p["LPratio"]
    y[IX["HCO3"]] = 24.0 if ctrl else 19.0
    y[IX["MalCoA"]] = 1.0 if ctrl else 2.2
    y[IX["TGliv"]] = 22.0 if ctrl else 80.0
    y[IX["TGpl"]] = 70.0 if ctrl else 700.0
    y[IX["FFA"]] = 0.45
    y[IX["KB"]] = 0.15
    y[IX["Ins"]] = 60.0
    y[IX["Gcg"]] = 12.0
    y[IX["Epi"]] = 0.4
    y[IX["Cort"]] = 250.0
    y[IX["GH"]] = 2.0
    y[IX["IGF1"]] = p["IGF1n"] * (1.0 if ctrl else 0.6)
    y[IX["PRPP"]] = 1.0 if ctrl else 2.0
    y[IX["UA"]] = 3.4 if ctrl else 7.5
    y[IX["GFRrel"]] = 1.0
    y[IX["CitU"]] = 1.0 if ctrl else 0.7
    y[IX["HtSDS"]] = 0.0 if ctrl else -1.6
    y[IX["BMDz"]] = 0.0 if ctrl else -1.1
    y[IX["AG15"]] = p["AG15in"] / p["kAG15"]
    if p["genotype"] == "Ib":
        y[IX["AG6Pn"]] = 1.60
        y[IX["ANC"]] = 0.5
    else:
        y[IX["AG6Pn"]] = 0.30
        y[IX["ANC"]] = 3.9
    return y


# ---------------------------------------------------------------------------
# 5. SIMULATION DRIVER (discrete dosing events)
# ---------------------------------------------------------------------------
def simulate(p, y0, tmax, events=None, dt_out=0.05, rtol=1e-6, atol=1e-8):
    """
    events: list of (time_h, state_name, amount) bolus additions, or
            (time_h, "PARAM:name", value) parameter changes.
    """
    events = sorted(events or [], key=lambda e: e[0])
    times = [0.0] + [e[0] for e in events] + [tmax]
    times = sorted(set(t for t in times if 0.0 <= t <= tmax))
    ts_all, ys_all = [], []
    y = np.array(y0, dtype=float)
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        for e in events:
            if abs(e[0] - t0) < 1e-9:
                if str(e[1]).startswith("PARAM:"):
                    p[e[1].split(":", 1)[1]] = e[2]
                else:
                    y[IX[e[1]]] += e[2]
        if t1 <= t0:
            continue
        n = max(int((t1 - t0) / dt_out) + 1, 2)
        teval = np.linspace(t0, t1, n)
        sol = solve_ivp(rhs, (t0, t1), y, args=(p,), method="LSODA",
                        t_eval=teval, rtol=rtol, atol=atol, max_step=(t1 - t0))
        if not sol.success:
            raise RuntimeError(f"integration failed on [{t0},{t1}]: {sol.message}")
        ts_all.append(sol.t if i == 0 else sol.t[1:])
        ys_all.append(sol.y if i == 0 else sol.y[:, 1:])
        y = sol.y[:, -1].copy()
    for e in events:
        if abs(e[0] - tmax) < 1e-9 and not str(e[1]).startswith("PARAM:"):
            y[IX[e[1]]] += e[2]
    T = np.concatenate(ts_all)
    Y = np.concatenate(ys_all, axis=1)
    return T, Y


def col(Y, name):
    return Y[IX[name], :]


def fluxes(y, p):
    """
    Re-evaluate the hepatic branch point and the whole-body glucose balance at
    a given state, and return every flux by name (mmol/h).  Used for auditing
    the mass balance and for the flux-partition figure in the Shiny app --
    the whole argument of this model is which way carbon leaves G6P, so that
    partition has to be inspectable rather than implied.
    """
    Y = {s: max(y[i], 0.0) for i, s in enumerate(STATES)}
    Gp, G6P, Glyc, Lac = Y["Gp"], Y["G6P"], Y["Glyc"], Y["Lac"]
    Ins, Gcg, Epi, Cort = Y["Ins"], Y["Gcg"], Y["Epi"], Y["Cort"]
    LWg = max(Y["LV"] * 1.05, 1.0)
    Insrel, Gcgrel, Epirel = Ins / 120.0, Gcg / 12.0, Epi / 0.5

    Jgk = (p["Vgk"] * (1.0 + p["gkIns"] * min(Insrel, 2.5))
           * hill(Gp, p["Kgk"], p["ngk"]) * (LWg / p["LW"]))
    phos_act = (0.10 + 2.2 * hill(Gcgrel, 2.5, 2.0) + 0.90 * hill(Epirel, 2.5, 2.0)) \
        / (1.0 + 0.9 * hill(Gp, 6.0, 3.0)) / (1.0 + 0.6 * min(Insrel, 3.0))
    Jgp = p["kgp"] * Glyc * phos_act * LWg / 1000.0
    Jdeb = p["fdeb"] * Jgp
    Jlys = p["klys"] * Glyc * LWg / 1000.0
    gng_horm = 0.45 + 0.75 * hill(Gcgrel, 1.0, 1.4) + 0.45 * (Cort / 400.0)
    Jgng = p["Vgng"] * gng_horm * hill(Lac + 0.6 * Y["FFA"] + p["Sub0"], p["Kgng"])
    a_enz = min(p["resid"] + Y["G6Pact"], 1.0)
    a_eff = a_enz * (p["fG6PT"] if p["genotype"] == "Ib" else 1.0)
    Jg6pase = p["VmaxG6P"] * a_eff * hill(G6P, p["KmG6P"]) * (LWg / p["LW"])
    Jextra = p["fEGPextra"] * p["VmaxG6P"] * a_eff * hill(G6P, p["KmG6P"])
    gs_cov = (0.15 + 1.6 * hill(Ins, 110.0, 1.5)) / (1.0 + 0.8 * hill(Gcgrel, 1.2, 2.0))
    gs_allo = 0.22 + 0.78 * hill(G6P, p["Kgs"], 2.0)
    Jgs = (p["Vgs"] * gs_cov * gs_allo * (1.0 - hill(Glyc, p["GlycMax"], 6.0))
           * (LWg / p["LW"]))
    f26 = 1.0 + 0.55 * min(Insrel, 2.0) - 0.20 * hill(Gcgrel, 1.0, 1.5)
    acid_brake = hill(Y["HCO3"], p["KacidPFK"], 3.0) / hill(24.0, p["KacidPFK"], 3.0)
    Jgly = (p["Vpfk"] * hill(G6P, p["Kpfk"]) * max(f26, 0.25)
            * acid_brake * (LWg / p["LW"]))
    Jppp = p["Vppp"] * hill(G6P, p["Kppp"]) * (LWg / p["LW"])
    sup_lac, need_glc, got_glc, FAI = cerebral_fuel(Gp, Lac, p)
    Uper = (p["Uper0"] + p["Si"] * Ins) * hill(Gp, 2.5) * (1.0 - 0.25 * p["illness"])
    EGP = Jg6pase + Jextra + Jdeb + Jlys
    tot_in = Jgk + (Jgp - Jdeb) + Jgng
    tot_out = Jg6pase + Jgs + Jgly + Jppp
    return dict(
        Jgk=Jgk, Jgp=Jgp, Jdeb=Jdeb, Jlys=Jlys, Jgng=Jgng, Jg6pase=Jg6pase,
        Jextra=Jextra, Jgs=Jgs, Jgly=Jgly, Jppp=Jppp,
        G6P_in=tot_in, G6P_out=tot_out, G6P_net=tot_in - tot_out,
        EGP=EGP, EGP_mg_kg_min=EGP * 3.0 / p["BW"],
        CNS_glucose=got_glc, CNS_lactate_eq=sup_lac, CNS_FAI=FAI,
        U_peripheral=Uper,
        glucose_disposal=got_glc + Uper + Jgk,
        # branch fractions out of G6P -- the headline partition
        frac_to_glucose=Jg6pase / max(tot_out, 1e-9),
        frac_to_glycogen=Jgs / max(tot_out, 1e-9),
        frac_to_lactate=Jgly / max(tot_out, 1e-9),
        frac_to_PPP=Jppp / max(tot_out, 1e-9),
    )


def cornstarch_events(dose_g_per_kg, p, interval_h, tmax, t_first=0.0):
    """Cornstarch dose -> mmol glucose equivalents in the slow-starch depot."""
    per_dose = dose_g_per_kg * p["BW"] * 0.97 * 1000.0 / 180.0   # mmol glucose-eq
    ev, t = [], t_first
    while t <= tmax + 1e-9:
        ev.append((t, "AST", per_dose))
        t += interval_h
    return ev, per_dose


# ---------------------------------------------------------------------------
# 6. ANALYSES
# ---------------------------------------------------------------------------
RESULTS = {}


def steady_state(p, hours=400.0, drip_frac=None):
    """
    Bring a patient to a repeating-day quasi-steady state on a continuous
    glucose delivery equal to `drip_frac` of total demand, so that the slow
    pools (glycogen, TG, urate, liver volume) are settled before an acute
    experiment is run on top of them.
    """
    p = dict(p)
    if drip_frac is None:
        # a healthy liver supplies its own glucose between meals; a GSD Ia
        # liver cannot, so its "baseline" is a delivered-glucose fraction
        drip_frac = 1.05 if p["genotype"] == "control" else 0.72
    p["Rdrip"] = drip_frac * p["GUR_tot"]
    y = y0_for(p)
    T, Y = simulate(p, y, hours, dt_out=0.5, rtol=1e-6, atol=1e-8)
    return Y[:, -1].copy()


def control_state(age_y=1.0, bw=10.0, target_glucose=5.0, hours=300.0):
    """
    Reference healthy subject.  A normal liver makes its own glucose, so its
    "baseline" cannot be defined by a delivery fraction the way a GSD Ia
    patient's is.  Instead we bisect the enteral delivery until the model
    settles at a normal postabsorptive plasma glucose, which puts glycogen,
    G6P, insulin and lipids in their physiological ranges simultaneously.
    """
    p = set_genotype(base_params(age_y, bw), "control")
    lo, hi = 0.0, 1.2
    y = None
    for _ in range(16):
        mid = 0.5 * (lo + hi)
        y = steady_state(p, hours=hours, drip_frac=mid)
        if y[IX["Gp"]] > target_glucose:
            hi = mid
        else:
            lo = mid
    p["Rdrip"] = 0.5 * (lo + hi) * p["GUR_tot"]
    return p, y


def a1_fasting_tolerance():
    """
    RESULT 1 -- cornstarch coverage is limited by its RELEASE RATE, not by the
    dose.

    The obvious model of cornstarch is a reservoir: swallow D grams, spend them
    at the deficit rate, and the interval is D/deficit.  That prediction is
    made below as `t_reservoir_h` and it is WRONG -- the simulation always
    falls short of it, and the gap widens with body size.

    The reason is that a first-order starch depot delivers glucose at
    kdis*A(t), which decays.  Coverage ends when the DELIVERY RATE drops below
    the deficit, with starch still in the gut:

        release(t) = kdis * D * Fabs * exp(-kdis*t)  >  deficit
        =>  t_cover = (1/kdis) * ln( kdis * D * Fabs / deficit )

    which is LOGARITHMIC in the dose.  Doubling the dose therefore adds exactly
    ln2/kdis hours regardless of where you started -- 1.54 h for uncooked
    cornstarch (kdis 0.45/h) and 2.48 h for extended-release waxy maize
    (kdis 0.28/h).  That single fixed increment is why "just give more
    cornstarch" stops working, and why a slower starch and not a bigger one is
    what buys an unbroken night.
    """
    out = {}
    profiles = [("infant_6mo", 0.5, 7.0), ("child_1y", 1.0, 10.0),
                ("child_5y", 5.0, 18.0), ("adol_14y", 14.0, 50.0),
                ("adult_30y", 30.0, 70.0)]
    for name, age, bw in profiles:
        p = set_genotype(base_params(age, bw), "Ia")
        yss = steady_state(p)
        # net glucose deficit evaluated AT euglycaemia, mmol/h -- this is the
        # denominator of the ratio.  Everything in it is read straight out of
        # the audited flux decomposition, so the "prediction" below is not a
        # separate model, it is the same model solved by hand.
        pp = dict(p)
        pp["Rdrip"] = 0.0
        y = yss.copy()
        y[IX["Gp"]] = 4.5
        fx = fluxes(y, pp)
        deficit = fx["glucose_disposal"] - fx["EGP"]          # mmol/h
        EGPres = fx["EGP"]
        dose_gkg = 1.6 if age < 10 else 1.45

        def sim_tolerance(dose, kdis=None):
            q = dict(pp)
            if kdis is not None:
                q["kdis"] = kdis
            ev, _ = cornstarch_events(dose, q, 1e9, 0.0)
            y2 = yss.copy()
            y2[IX["Gp"]] = 4.2
            T, Y = simulate(q, y2, 30.0, events=ev, dt_out=0.02)
            below = np.where(col(Y, "Gp") < 3.3)[0]
            return float(T[below[0]]) if len(below) else 30.0

        def t_rate_limited(dose, kdis):
            D = dose * bw * 0.97 * 1000.0 / 180.0
            arg = kdis * D * pp["Fabs"] / max(deficit, 1e-9)
            return math.log(arg) / kdis if arg > 1 else 0.0

        delivered = dose_gkg * bw * 0.97 * 1000.0 / 180.0 * pp["Fabs"]   # mmol
        t_reservoir = delivered / max(deficit, 1e-6)
        t_meas = sim_tolerance(dose_gkg)
        out[name] = dict(
            age_y=age, bw_kg=bw,
            GUR_mg_kg_min=round(pp["GUR_mgkgmin"], 2),
            deficit_mmol_h=round(deficit, 2),
            deficit_g_h=round(deficit * 0.180, 2),
            deficit_mg_kg_min=round(deficit * 3.0 / bw, 2),
            residual_EGP_mg_kg_min=round(EGPres * 3.0 / bw, 2),
            residual_EGP_pct_of_demand=round(
                100.0 * EGPres * 3.0 / bw / pp["GUR_mgkgmin"], 1),
            cerebral_lactate_share=round(
                fx["CNS_lactate_eq"] / pp["Dcns"], 3),
            dose_g_per_kg=dose_gkg,
            delivered_mmol=round(delivered, 1),
            # the naive reservoir prediction -- always too optimistic
            t_reservoir_h=round(t_reservoir, 2),
            # the rate-limited prediction -- the one that matches
            t_rate_limited_h=round(t_rate_limited(dose_gkg, pp["kdis"]), 2),
            t_simulated_h=round(t_meas, 2),
            # doubling the dose
            t_simulated_double_dose_h=round(sim_tolerance(2 * dose_gkg), 2),
            gain_from_doubling_h=round(
                sim_tolerance(2 * dose_gkg) - t_meas, 2),
            # same dose, extended-release starch instead
            t_simulated_slow_release_h=round(
                sim_tolerance(dose_gkg, kdis=0.28), 2),
            gain_from_slower_starch_h=round(
                sim_tolerance(dose_gkg, kdis=0.28) - t_meas, 2),
        )
    RESULTS["R1_fasting_tolerance"] = dict(
        by_age=out,
        predicted_gain_from_doubling_h=round(math.log(2) / 0.45, 2),
        predicted_gain_switching_to_kdis_0p28_h=None,
        note=("gain from doubling the dose = ln2/kdis = 1.54 h at kdis 0.45/h, "
              "independent of dose and of body size -- compare the simulated "
              "gain_from_doubling_h column across the five ages"),
    )
    return RESULTS["R1_fasting_tolerance"]


def a2_cerebral_fuel():
    """
    RESULT 2 -- glucose and neuroglycopenia are separated by lactate.
    """
    p = set_genotype(base_params(1.0, 10.0), "Ia")
    tab = {}
    for lac in [0.8, 1.5, 2.5, 4.0, 6.0, 8.0, 10.0]:
        g = isofuel_glucose(lac, p)
        sup, need, got, fai = cerebral_fuel(g if math.isfinite(g) else 0.0, lac, p)
        tab[f"lac_{lac}"] = dict(
            lactate_mM=lac,
            isofuel_glucose_mM=round(g, 2) if math.isfinite(g) else None,
            lactate_share_of_brain_fuel=round(sup / p["Dcns"], 3),
        )
    # the clinical trap: normalise lactate without raising glucose
    g_low = 1.9
    _, _, _, fai_before = cerebral_fuel(g_low, 6.5, p)
    _, _, _, fai_after = cerebral_fuel(g_low, 1.6, p)
    RESULTS["R2_cerebral_fuel"] = dict(
        isofuel_curve=tab,
        threshold_FAI=p["FAIthr"],
        trap=dict(
            glucose_mM=g_low,
            FAI_at_lactate_6p5=round(fai_before, 3),
            FAI_at_lactate_1p6=round(fai_after, 3),
            symptomatic_before=bool(fai_before < p["FAIthr"]),
            symptomatic_after=bool(fai_after < p["FAIthr"]),
        ),
    )
    return RESULTS["R2_cerebral_fuel"]


def a3_counterregulation():
    """
    RESULT 3 -- in GSD Ia the counterregulatory response is not merely futile,
    it is harmful: glucagon raises lactate without raising glucose.
    """
    out = {}
    for gtype in ["control", "Ia"]:
        if gtype == "control":
            p, yss = control_state(1.0, 10.0)
        else:
            p = set_genotype(base_params(1.0, 10.0), "Ia")
            yss = steady_state(p, drip_frac=0.80)
        # glucagon challenge: 30 ug/kg IM equivalent -> clamp Gcg high 60 min
        pp = dict(p)
        pp["Rdrip"] = 0.0
        y = yss.copy()
        base_G, base_L = y[IX["Gp"]], y[IX["Lac"]]
        pp2 = dict(pp)
        pp2["Gcgmax"] = pp["Gcgmax"] * 20.0
        T, Y = simulate(pp2, y, 1.0, dt_out=0.01)
        y1 = Y[:, -1]
        T2, Y2 = simulate(pp, y1, 1.0, dt_out=0.01)
        Gtr = np.concatenate([col(Y, "Gp"), col(Y2, "Gp")])
        Ltr = np.concatenate([col(Y, "Lac"), col(Y2, "Lac")])
        out[gtype] = dict(
            baseline_glucose_mM=round(float(base_G), 2),
            peak_glucose_mM=round(float(Gtr.max()), 2),
            delta_glucose_mM=round(float(Gtr.max() - base_G), 2),
            baseline_lactate_mM=round(float(base_L), 2),
            peak_lactate_mM=round(float(Ltr.max()), 2),
            delta_lactate_mM=round(float(Ltr.max() - base_L), 2),
        )
    ia, ct = out["Ia"], out["control"]
    ratio = (round(ia["delta_lactate_mM"] / ia["delta_glucose_mM"], 2)
             if ia["delta_glucose_mM"] > 0.05 else None)
    out["interpretation"] = dict(
        glucose_response_ratio_Ia_vs_control=round(
            ia["delta_glucose_mM"] / max(ct["delta_glucose_mM"], 1e-9), 3),
        lactate_penalty_mM_per_mM_glucose_gained_Ia=ratio,
        lactate_penalty_is_unbounded=(ratio is None),
        reading=("in GSD Ia the glucagon response mobilises glycogen into a "
                 "blocked pathway: glucose does not move at all, so the "
                 "lactate cost per mmol/L of glucose gained is not merely "
                 "large, it is undefined -- there is no gain to divide by"),
    )
    RESULTS["R3_counterregulation"] = out
    return out


def a4_restored_activity():
    """
    RESULT 4 -- daily cornstarch requirement falls LINEARLY in restored
    G6Pase activity while fasting tolerance rises HYPERBOLICALLY, so there is
    a critical activity a* at which an 8-hour night becomes safe and below
    which nothing clinically changes.
    """
    p0 = set_genotype(base_params(14.0, 50.0), "Ia")
    rows = []
    for a in [0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.40]:
        p = dict(p0)
        p["resid"] = a
        yss = steady_state(p)
        pp = dict(p)
        pp["Rdrip"] = 0.0
        # (i) fasting tolerance: one 1.45 g/kg dose, time to Gp < 3.3
        ev, _ = cornstarch_events(1.45, pp, 1e9, 0.0)
        y = yss.copy()
        y[IX["Gp"]] = 4.2
        T, Y = simulate(pp, y, 30.0, events=ev, dt_out=0.02)
        G = col(Y, "Gp")
        below = np.where(G < 3.3)[0]
        t_tol = float(T[below[0]]) if len(below) else 30.0
        # (ii) daily requirement: continuous delivery needed to hold Gp >= 4.0
        lo, hi = 0.0, 1.4
        for _ in range(22):
            mid = 0.5 * (lo + hi)
            q = dict(pp)
            q["Rdrip"] = mid * q["GUR_tot"]
            _, Yq = simulate(q, yss.copy(), 60.0, dt_out=1.0)
            if col(Yq, "Gp")[-1] >= 4.0:
                hi = mid
            else:
                lo = mid
        req_frac = hi
        req_g_day = req_frac * p0["GUR_tot"] * 24.0 * 0.180 / 0.97 / p0["Fabs"]
        rows.append(dict(activity=a,
                         fasting_tolerance_h=round(t_tol, 2),
                         daily_starch_g=round(req_g_day, 1),
                         lactate_mM=round(float(yss[IX["Lac"]]), 2),
                         urate_mg_dL=round(float(yss[IX["UA"]]), 2),
                         TG_mg_dL=round(float(yss[IX["TGpl"]]), 0),
                         glycogen_mg_g=round(float(yss[IX["Glyc"]]) * 0.162, 1)))
    base = rows[0]
    a_star = None
    for r in rows:
        if r["fasting_tolerance_h"] >= 8.0:
            a_star = r["activity"]
            break
    RESULTS["R4_restored_activity"] = dict(
        rows=rows,
        a_star_for_8h_night=a_star,
        starch_reduction_at_a_star=(
            round(100 * (1 - [r for r in rows if r["activity"] == a_star][0]
                         ["daily_starch_g"] / base["daily_starch_g"]), 1)
            if a_star is not None else None),
        note=("cornstarch requirement is linear in a; fasting tolerance is "
              "1/(deficit - a*Vmax) and therefore diverges"),
    )
    return RESULTS["R4_restored_activity"]


def a5_aav_dilution():
    """
    RESULT 5 -- AAV durability in children is set by liver GROWTH, not by
    promoter silencing, and the decision is irreversible because anti-capsid
    antibody blocks redosing.
    """
    out = {}
    for age, bw in [(2.0, 12.0), (6.0, 20.0), (12.0, 40.0), (18.0, 62.0), (30.0, 70.0)]:
        p = set_genotype(base_params(age, bw), "Ia")
        LW_now = p["LW"]
        LW_adult = p["LWadult"]
        growth_mult = LW_adult / LW_now
        yss = steady_state(p)
        y = yss.copy()
        y[IX["AAVvg"]] = 0.22          # transduced load giving a0 = 0.22
        y[IX["G6Pact"]] = 0.0
        y[IX["AntiAAV"]] = 1.0
        q = dict(p)
        q["Rdrip"] = 0.72 * q["GUR_tot"]
        q["chronic_scale"] = 1.0
        horizon = min(30.0, max(1.0, 45.0 - age)) * 8766.0
        T, Y = simulate(q, y, horizon, dt_out=720.0, rtol=1e-5, atol=1e-7)
        act = col(Y, "G6Pact")
        peak = float(act.max())
        t_pk = float(T[int(act.argmax())])
        # activity at +5 y, +10 y, +20 y after dosing
        def at(yrs):
            tt = yrs * 8766.0
            if tt > T[-1]:
                return None
            return round(float(np.interp(tt, T, act)), 4)
        out[f"age_{int(age)}y"] = dict(
            age_at_dosing_y=age,
            liver_g_now=round(LW_now, 0),
            remaining_liver_growth_fold=round(growth_mult, 2),
            peak_activity=round(peak, 4),
            time_to_peak_d=round(t_pk / 24.0, 1),
            activity_5y=at(5), activity_10y=at(10), activity_20y=at(20),
            retained_fraction_10y=(round(at(10) / peak, 3)
                                   if at(10) is not None else None),
        )
    RESULTS["R5_aav_dilution"] = out
    return out


def a6_urate_lever():
    """
    RESULT 6 -- urate is a QUOTIENT, so its two levers MULTIPLY rather than add.

        urate_ss = production x (1 - E_allopurinol)
                   -------------------------------------------
                   clearance / (1 + lactate/Ki + ketones/Ki')

    The author's prior hypothesis going in was that lactate control would beat
    allopurinol outright, because lactate enters the denominator AND feeds the
    numerator through the pentose phosphate pathway while allopurinol touches
    only the numerator.  The simulation REFUTED that: the two single-agent
    effects come out nearly equal.  What it did confirm is the structural
    prediction that follows from the quotient -- the combination is markedly
    SUB-additive, so a clinician who adds allopurinol to a patient whose diet
    has just been intensified should expect considerably less than the sum of
    the two effects, and should not read the shortfall as non-adherence.
    """
    p = set_genotype(base_params(14.0, 50.0), "Ia")
    scenarios = {}

    def run(tag, frac=0.85, amp=0.95, allo_mg=None, hours=24 * 28):
        q = dict(p)
        q["RdripAmp"] = amp
        y = steady_state(q, hours=900.0, drip_frac=frac)
        q["Rdrip"] = frac * q["GUR_tot"]
        ev = []
        if allo_mg:
            t = 0.0
            while t < hours:
                ev.append((t, "ALLOg", allo_mg))
                t += 24.0
        T, Y = simulate(q, y, hours, events=ev, dt_out=2.0)
        last = T >= (hours - 24.0)          # average over the final day
        scenarios[tag] = dict(
            mean_lactate_mM=round(float(col(Y, "Lac")[last].mean()), 2),
            peak_lactate_mM=round(float(col(Y, "Lac")[last].max()), 2),
            urate_mg_dL=round(float(col(Y, "UA")[last].mean()), 2),
            oxypurinol_mg_L=round(float(col(Y, "OXY")[last].mean()), 2),
        )
        return scenarios[tag]

    base = run("baseline_poor_control")
    allo = run("allopurinol_300mg_only", allo_mg=300.0)
    diet = run("intensified_diet_only", frac=0.92, amp=0.20)
    both = run("diet_plus_allopurinol", frac=0.92, amp=0.20, allo_mg=300.0)
    d_allo = allo["urate_mg_dL"] - base["urate_mg_dL"]
    d_diet = diet["urate_mg_dL"] - base["urate_mg_dL"]
    d_both = both["urate_mg_dL"] - base["urate_mg_dL"]
    RESULTS["R6_urate"] = dict(
        scenarios=scenarios,
        delta_allopurinol=round(d_allo, 2),
        delta_diet=round(d_diet, 2),
        delta_both_observed=round(d_both, 2),
        delta_both_if_additive=round(d_allo + d_diet, 2),
        shortfall_vs_additive=round((d_allo + d_diet) - d_both, 2),
        multiplicative_prediction=round(
            base["urate_mg_dL"]
            * (allo["urate_mg_dL"] / base["urate_mg_dL"])
            * (diet["urate_mg_dL"] / base["urate_mg_dL"])
            - base["urate_mg_dL"], 2),
        prior_hypothesis="lactate control >> allopurinol",
        verdict_on_prior="REFUTED -- single-agent effects are nearly equal",
    )
    return RESULTS["R6_urate"]


def a7_ketone_signature():
    """
    CHECK -- the model must produce HYPOketotic hypoglycaemia in Ia (it is the
    single most useful bedside discriminator from GSD 0/III/VI) without having
    been told to.  Achieved via the malonyl-CoA/CPT1 brake.
    """
    out = {}
    for gtype in ["control", "Ia"]:
        if gtype == "control":
            p, yss = control_state(1.0, 10.0)
        else:
            p = set_genotype(base_params(1.0, 10.0), "Ia")
            yss = steady_state(p, drip_frac=0.80)
        q = dict(p)
        q["Rdrip"] = 0.0
        T, Y = simulate(q, yss, 8.0, dt_out=0.05)
        i = int(np.argmin(np.abs(T - 6.0)))
        out[gtype] = dict(
            glucose_6h_mM=round(float(col(Y, "Gp")[i]), 2),
            ketone_6h_mM=round(float(col(Y, "KB")[i]), 3),
            FFA_6h_mM=round(float(col(Y, "FFA")[i]), 2),
            lactate_6h_mM=round(float(col(Y, "Lac")[i]), 2),
            malonylCoA_6h_rel=round(float(col(Y, "MalCoA")[i]), 2),
        )
    ia = out["Ia"]
    out["ketone_to_FFA_ratio_Ia"] = round(ia["ketone_6h_mM"] / max(ia["FFA_6h_mM"], 1e-9), 3)
    RESULTS["R7_ketone_signature"] = out
    return out


def a8_gsd1b_empagliflozin():
    """
    RESULT 7 -- in GSD Ib the neutrophil count is a THRESHOLD function of
    1,5-AG6P, so empagliflozin looks like an all-or-nothing drug.
    """
    p = set_genotype(base_params(12.0, 38.0), "Ib")
    y = steady_state(p, hours=2000.0)
    q = dict(p)
    q["Rdrip"] = 0.72 * q["GUR_tot"]
    ev = []
    t = 24.0 * 30
    while t < 24.0 * 400:
        ev.append((t, "EMPAc", 780.0))     # ~0.4 mg/kg/day split, nmol/L bolus
        t += 24.0
    T, Y = simulate(q, y, 24.0 * 400, events=ev, dt_out=12.0)
    def at(day, name):
        return round(float(np.interp(day * 24.0, T, col(Y, name))), 3)
    RESULTS["R8_gsd1b"] = dict(
        baseline=dict(ANC=at(29, "ANC"), AG15=at(29, "AG15"), AG6Pn=at(29, "AG6Pn")),
        week4=dict(ANC=at(58, "ANC"), AG15=at(58, "AG15"), AG6Pn=at(58, "AG6Pn")),
        week12=dict(ANC=at(114, "ANC"), AG15=at(114, "AG15"), AG6Pn=at(114, "AG6Pn")),
        month12=dict(ANC=at(395, "ANC"), AG15=at(395, "AG15"), AG6Pn=at(395, "AG6Pn")),
    )
    return RESULTS["R8_gsd1b"]


def a9_regimens():
    """
    Head-to-head of the four real overnight strategies in a 5-year-old.
    """
    p = set_genotype(base_params(5.0, 18.0), "Ia")
    yss = steady_state(p)
    out = {}

    # 9 h is the question that is actually asked of an overnight regimen in a
    # school-age child; 12 h flatters continuous feeding and buries the rest.
    def score(tag, q, ev, tmax=9.0):
        y = yss.copy()
        y[IX["Gp"]] = 4.6
        T, Y = simulate(q, y, tmax, events=ev, dt_out=0.02)
        G, L = col(Y, "Gp"), col(Y, "Lac")
        FAI = np.array([cerebral_fuel(G[i], L[i], q)[3] for i in range(len(T))])
        out[tag] = dict(
            nadir_glucose_mM=round(float(G.min()), 2),
            time_below_3p9_h=round(float(np.trapezoid((G < 3.9).astype(float), T)), 2),
            time_below_3p0_h=round(float(np.trapezoid((G < 3.0).astype(float), T)), 2),
            time_FAI_below_thr_h=round(
                float(np.trapezoid((FAI < q["FAIthr"]).astype(float), T)), 2),
            mean_lactate_mM=round(float(L.mean()), 2),
            peak_lactate_mM=round(float(L.max()), 2),
            mean_insulin_pM=round(float(col(Y, "Ins").mean()), 1),
        )
        return out[tag]

    q = dict(p); q["Rdrip"] = 0.0
    ev, _ = cornstarch_events(1.6, q, 4.0, 9.0)
    score("UCCS_1.6gkg_q4h", q, ev)

    q = dict(p); q["Rdrip"] = 0.0
    ev, _ = cornstarch_events(1.6, q, 6.0, 9.0)
    score("UCCS_1.6gkg_q6h", q, ev)

    q = dict(p); q["Rdrip"] = 0.0
    ev, _ = cornstarch_events(2.4, q, 9.0, 0.0)   # single large UCCS dose
    score("UCCS_2.4gkg_single_dose", q, ev)

    q = dict(p); q["Rdrip"] = 0.0
    q["kdis"] = 0.28; q["Fabs"] = 0.78          # extended-release waxy maize
    ev, _ = cornstarch_events(2.0, q, 9.0, 0.0)
    score("Glycosade_2.0gkg_nightly", q, ev)

    q = dict(p)
    # continuous nasogastric drip at 7 mg/kg/min glucose polymer, F ~ 0.98
    q["Rdrip"] = 7.0 * q["BW"] / 3.0 * 0.98
    score("continuous_drip_7mgkgmin", q, [])

    # the classic disaster: the pump disconnects at 03:00.  Note what the
    # comparison shows -- the SAFEST regimen while it is running is the one
    # with the least reserve when it stops, because a continuously fed liver
    # has been kept insulinised and is not mobilising anything.
    q = dict(p)
    q["Rdrip"] = 7.0 * q["BW"] / 3.0 * 0.98
    ev = [(5.0, "PARAM:Rdrip", 0.0)]
    score("drip_then_pump_failure_at_5h", q, ev)

    RESULTS["R9_regimens"] = out
    return out


SLOW = ["GFRrel", "Rdam", "UACR", "CitU", "NCa", "HCA", "HCChaz",
        "HtSDS", "BMDz"]


def day_average(p, y0, hours=48.0):
    """
    Run one settled 24-hour cycle at fine resolution and return the daily-mean
    drivers.  Averaging over the day is not a shortcut for its own sake: the
    slow endpoints below respond to chronic exposure, and a 30-year integration
    that resolved 11 000 individual day-cycles would spend all of its effort on
    a waveform that is identical every time.
    """
    q = dict(p)
    T, Y = simulate(q, y0, hours, dt_out=0.05)
    m = T >= (hours - 24.0)
    out = {}
    for s in ["Gp", "Lac", "TGpl", "UA", "HCO3", "LV", "Glyc", "G6P", "KB"]:
        out[s] = float(col(Y, s)[m].mean())
    out["Gp_min"] = float(col(Y, "Gp")[m].min())
    out["Lac_max"] = float(col(Y, "Lac")[m].max())
    return out


def slow_rhs(t, z, drv, p, acei_eff):
    """
    The slow subsystem, driven by frozen daily-mean biochemistry.  Written out
    separately from rhs() rather than reusing it, precisely so that a
    discrepancy between the two would show up as a discrepancy in the answer.
    """
    GFRrel, Rdam, UACR, CitU, NCa, HCA, HCChaz, HtSDS, BMDz = z
    Lac, TGpl, UA, HCO3, Gp = (drv["Lac"], drv["TGpl"], drv["UA"],
                               drv["HCO3"], drv["Gp"])
    d = np.zeros(9)
    GFRtar = 1.0 + 0.30 * hill(max(Lac - 2.0, 0.0), 4.0) - 0.9 * Rdam
    d[0] = p["khf"] * (GFRtar - GFRrel)
    drive = (p["wLac"] * max(Lac - 2.0, 0.0)
             + p["wHF"] * max(GFRrel - 1.0, 0.0)
             + p["wUA"] * max(UA - 5.5, 0.0)
             + p["wTG"] * TGpl)
    d[1] = p["kdam"] * drive * (1.0 - acei_eff) - p["krep"] * Rdam
    d[2] = p["kUACR"] * (Rdam * 2600.0 * (1.0 - acei_eff) - UACR)
    d[3] = p["kCit"] * ((HCO3 / p["HCO3n"]) ** 2.2 * (1.0 + 0.6 * p["citrate"]) - CitU)
    d[4] = p["kNCa"] * max(1.0 - CitU, 0.0) * (1.0 + 0.4 * max(UA - 6.0, 0.0))
    MCI = min(1.0, (0.34 * max(Lac - 2.2, 0.0) / 5.0
                    + 0.26 * max(TGpl - 250.0, 0.0) / 900.0
                    + 0.22 * max(UA - 5.5, 0.0) / 5.0
                    + 0.18 * max(3.9 - drv["Gp_min"], 0.0) / 2.0) * 1.9)
    d[5] = (p["kHCAini"]
            + p["kHCAgro"] * HCA * max(1.0 - HCA / p["HCAmax"], 0.0)
            ) * (MCI ** p["nMCI"]) * 4.0
    d[6] = p["kHCC"] * HCA * (1.0 + 2.0 * MCI)
    acid = max(p["HCO3n"] - HCO3, 0.0) / 8.0
    age_now = p["age"] + t / 8766.0
    growing = 1.0 / (1.0 + math.exp((age_now - p["ageMature"]) / 1.5))
    d[7] = p["kgro"] * growing * ((p["HtSDStar0"] - p["HtSDSmci"] * MCI) - HtSDS)
    d[8] = p["kbmd"] * ((p["BMDztar0"] - p["BMDzAcid"] * acid
                         - p["BMDzMci"] * MCI) - BMDz)
    return d


def a10_lifetime():
    """
    RESULT 8 -- thirty years of good versus poor metabolic control, with and
    without an ACE inhibitor.

    Integrated as an explicit TWO-TIMESCALE problem: the fast system is solved
    for one settled day, its daily-mean biochemistry is frozen, and only the
    nine slow states are carried over 30 years.  The metabolic control index
    enters the adenoma term with exponent 2.6, so this is where the model makes
    its strongest long-horizon claim -- that the difference between good and
    poor control in adenoma burden is far larger than the difference in any
    single biomarker that was used to define it.
    """
    out = {}
    for tag, drip, amp, acei in [("poor_control", 0.85, 0.95, False),
                                 ("good_control", 0.92, 0.20, False),
                                 ("good_control_plus_ACEi", 0.92, 0.20, True)]:
        p = set_genotype(base_params(8.0, 26.0), "Ia")
        p["RdripAmp"] = amp
        y = steady_state(p, hours=900.0, drip_frac=drip)
        q = dict(p)
        q["Rdrip"] = drip * q["GUR_tot"]
        drv = day_average(q, y)
        acei_eff = (p["EmaxACEI"] * hill(55.0 / (p["kACEI"] * 24.0), p["IC50ACEI"])
                    if acei else 0.0)
        z0 = np.array([y[IX[s]] for s in SLOW])
        Ts = np.linspace(0.0, 30 * 8766.0, 3001)
        sol = solve_ivp(slow_rhs, (0.0, Ts[-1]), z0, args=(drv, p, acei_eff),
                        method="LSODA", t_eval=Ts, rtol=1e-8, atol=1e-11)
        Z = {s: sol.y[i] for i, s in enumerate(SLOW)}

        def at(yrs, name):
            return round(float(np.interp(yrs * 8766.0, sol.t, Z[name])), 4)

        out[tag] = dict(
            mean_glucose=round(drv["Gp"], 2), min_glucose=round(drv["Gp_min"], 2),
            mean_lactate=round(drv["Lac"], 2), peak_lactate=round(drv["Lac_max"], 2),
            urate=round(drv["UA"], 2), TG=round(drv["TGpl"], 0),
            bicarbonate=round(drv["HCO3"], 1),
            liver_volume_ratio=round(
                drv["LV"] / (p["fLW"] * p["BW"] * 1000.0 / 1.05), 2),
            HCA_10y=at(10, "HCA"), HCA_20y=at(20, "HCA"), HCA_30y=at(30, "HCA"),
            HCC_hazard_30y=at(30, "HCChaz"),
            UACR_10y=at(10, "UACR"), UACR_20y=at(20, "UACR"),
            UACR_30y=at(30, "UACR"),
            GFRrel_30y=at(30, "GFRrel"),
            HtSDS_final=at(30, "HtSDS"), BMDz_final=at(30, "BMDz"),
            nephrocalcinosis_30y=at(30, "NCa"),
        )
    g, b = out["good_control"], out["poor_control"]
    out["contrast"] = dict(
        lactate_ratio_poor_over_good=round(b["mean_lactate"] / g["mean_lactate"], 2),
        adenoma_ratio_poor_over_good_30y=round(
            b["HCA_30y"] / max(g["HCA_30y"], 1e-9), 1),
        UACR_ratio_poor_over_good_30y=round(
            b["UACR_30y"] / max(g["UACR_30y"], 1e-9), 1),
        reading=("the ratio in mean lactate is small and the ratio in adenoma "
                 "burden is large, because the control index enters the "
                 "adenoma term with exponent 2.6 AND that term is "
                 "self-amplifying -- so a biomarker difference a clinic would "
                 "call marginal is not marginal in the endpoint it predicts"),
        ACEi_UACR_reduction_30y_pct=round(
            100.0 * (1 - out["good_control_plus_ACEi"]["UACR_30y"]
                     / max(g["UACR_30y"], 1e-9)), 1),
    )
    RESULTS["R10_lifetime"] = out
    return out


def a11_illness():
    """
    Intercurrent illness with vomiting -- the commonest route to an ICU
    admission in GSD I.  Shows the coupling: catabolic stress raises lactate
    while removing the enteral glucose that was holding the system up.
    """
    p = set_genotype(base_params(2.0, 12.0), "Ia")
    yss = steady_state(p)
    q = dict(p); q["Rdrip"] = 0.0
    ev, _ = cornstarch_events(1.6, q, 4.0, 8.0)
    ev += [(8.0, "PARAM:illness", 1.0), (8.0, "PARAM:Fabs", 0.10)]
    y = yss.copy(); y[IX["Gp"]] = 4.5
    T, Y = simulate(q, y, 24.0, events=ev, dt_out=0.05)
    i12 = int(np.argmin(np.abs(T - 12.0)))
    i24 = -1
    untreated = dict(
        glucose_12h=round(float(col(Y, "Gp")[i12]), 2),
        glucose_24h=round(float(col(Y, "Gp")[i24]), 2),
        lactate_24h=round(float(col(Y, "Lac")[i24]), 2),
        HCO3_24h=round(float(col(Y, "HCO3")[i24]), 1),
        FAI_24h=round(cerebral_fuel(col(Y, "Gp")[i24], col(Y, "Lac")[i24], q)[3], 3),
    )
    # rescue: IV dextrose 8 mg/kg/min from h=12
    q2 = dict(p); q2["Rdrip"] = 0.0
    ev2 = list(ev) + [(12.0, "PARAM:Rdrip", 8.0 * q2["BW"] / 3.0)]
    y2 = yss.copy(); y2[IX["Gp"]] = 4.5
    T2, Y2 = simulate(q2, y2, 24.0, events=ev2, dt_out=0.05)
    rescued = dict(
        glucose_24h=round(float(col(Y2, "Gp")[-1]), 2),
        lactate_24h=round(float(col(Y2, "Lac")[-1]), 2),
        HCO3_24h=round(float(col(Y2, "HCO3")[-1]), 1),
    )
    RESULTS["R11_illness"] = dict(no_rescue=untreated, iv_dextrose_8mgkgmin=rescued)
    return RESULTS["R11_illness"]


def a12_overtreatment():
    """
    The mirror-image error: chasing normoglycaemia with excess carbohydrate.
    In GSD Ia every surplus glucose molecule is trapped as G6P and leaves as
    glycogen, lactate or fat, so overtreatment has a metabolic cost that a
    normal liver does not pay.
    """
    out = {}
    p = set_genotype(base_params(8.0, 26.0), "Ia")
    p["RdripAmp"] = 0.30
    for tag, frac in [("under_0.60", 0.60), ("under_0.75", 0.75),
                      ("target_0.90", 0.90), ("over_1.05", 1.05),
                      ("over_1.25", 1.25), ("over_1.45", 1.45)]:
        y = steady_state(p, hours=1200.0, drip_frac=frac)
        out[tag] = dict(
            delivery_frac_of_demand=frac,
            glucose_mM=round(float(y[IX["Gp"]]), 2),
            lactate_mM=round(float(y[IX["Lac"]]), 2),
            insulin_pM=round(float(y[IX["Ins"]]), 1),
            glycogen_mg_g=round(float(y[IX["Glyc"]]) * 0.162, 1),
            liver_TG_mg_g=round(float(y[IX["TGliv"]]), 1),
            plasma_TG_mg_dL=round(float(y[IX["TGpl"]]), 0),
            urate_mg_dL=round(float(y[IX["UA"]]), 2),
            liver_volume_ratio=round(
                float(y[IX["LV"]]) / (p["fLW"] * p["BW"] * 1000.0 / 1.05), 2),
        )
    RESULTS["R12_overtreatment"] = out
    return out


# ---------------------------------------------------------------------------
# 7. MAIN
# ---------------------------------------------------------------------------
def main():
    np.seterr(all="ignore")
    order = [
        ("R1  fasting tolerance is a ratio", a1_fasting_tolerance),
        ("R2  cerebral fuel / lactate substitution", a2_cerebral_fuel),
        ("R3  counterregulation is harmful", a3_counterregulation),
        ("R4  restored activity: linear vs hyperbolic", a4_restored_activity),
        ("R5  AAV dilution by liver growth", a5_aav_dilution),
        ("R6  urate is a quotient", a6_urate_lever),
        ("R7  hypoketotic signature (falsifiable check)", a7_ketone_signature),
        ("R8  GSD Ib / empagliflozin threshold", a8_gsd1b_empagliflozin),
        ("R9  overnight regimen head-to-head", a9_regimens),
        ("R10 30-year natural history", a10_lifetime),
        ("R11 intercurrent illness", a11_illness),
        ("R12 overtreatment has a cost", a12_overtreatment),
    ]
    for label, fn in order:
        sys.stderr.write(f"[run] {label} ...\n")
        sys.stderr.flush()
        res = fn()
        print("=" * 78)
        print(label)
        print("=" * 78)
        print(json.dumps(res, indent=2, ensure_ascii=False))
        print()

    with open(os.path.join(HERE, "gsd1a_scenario_results.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=2, ensure_ascii=False)
    sys.stderr.write("[done] wrote gsd1a_scenario_results.json\n")


if __name__ == "__main__":
    main()
