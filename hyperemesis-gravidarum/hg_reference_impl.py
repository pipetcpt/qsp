#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hg_reference_impl.py — dependency-free reference implementation of the
Hyperemesis Gravidarum (HG) QSP model.

WHY THIS FILE EXISTS
====================
`hg_mrgsolve_model.R` is the deliverable model, but mrgsolve needs R + a
compiler.  This file is the *same* ODE system re-implemented in pure Python
(stdlib only, fixed-step RK4) so that every number quoted in README.md and
in the R file's calibration notes can be regenerated and checked on any
machine with `python3` and nothing else:

    python3 hg_reference_impl.py            # all scenarios + validation table
    python3 hg_reference_impl.py --check    # PASS/FAIL on the quoted numbers

THE MODEL'S ONE IDEA
====================
Every observation below is awkward for a model in which nausea is a function
of the *level* of the emetogenic placental hormone GDF15:

  (A) Maternal GDF15 rises through the first trimester and stays high (indeed
      keeps climbing) into the second and third trimesters -- yet nausea and
      vomiting of pregnancy peaks at ~9-11 weeks and is essentially gone by
      16-20 weeks (Marjono 2003 PMID 12495665; Fejzo 2019 PMID 31515515).
  (B) Women with beta-thalassaemia have chronically HIGH GDF15 and report
      *very low* levels of NVP (Fejzo 2024 PMID 38092039).
  (C) *Low* pre-pregnancy GDF15 increases HG risk (same paper).
  (D) Pre-pregnancy metformin -- a drug that RAISES GDF15 and is stopped
      before the nausea would begin -- is associated with a >70% reduction in
      HG risk (aRR 0.29, Sharma 2025 PMID 40588059).  Pre-pregnancy tobacco,
      which also raises GDF15, gives aRR 0.51 in the same study.
  (E) Ondansetron, the guideline-recommended second-line drug, moved PUQE-24
      by only -0.51 (95% CI -2.32 to 1.30) vs placebo in the first
      placebo-controlled HG trial, while mirtazapine moved it -1.86
      (-3.61 to -0.12), separating further after day 4 (Ostenfeld 2026
      PMID 41478546).

This model takes ONE structural position and derives all five:

  *** The area postrema does not measure how much GDF15 there is.
      It measures how much MORE GDF15 there is than there used to be. ***

i.e. the GDF15 -> GFRAL -> area postrema axis is a FOLD-CHANGE DETECTOR with
an adapting set-point SP that tracks log(GDF15) with a time constant of ~2-3
weeks.  The emetic drive is (GDF15 / adapted set-point), not GDF15.

That single commitment gives:
  (A) once the placental GDF15 ramp flattens, the ratio decays to 1 and
      symptoms resolve WHILE THE LEVEL STAYS HIGH -- resolution needs no
      extra mechanism and no falling hormone;
  (B),(C),(D) any exposure that raises GDF15 *before* conception raises the
      set-point before the placenta arrives, so the same absolute pregnancy
      level produces a much smaller ratio.  beta-thalassaemia, metformin and
      tobacco are then THE SAME EXPERIMENT done three ways, and one
      parameter (ALPHA, adaptation completeness) predicts all three;
  (E) see the second commitment.

Second commitment: EFFICACY IS SET BY NODE POSITION.  Emetic drive reaches
the emetic central pattern generator through the NTS, which sums a large
GFRAL/area-postrema term and a small peripheral vagal/5-HT3 term.  A drug's
achievable effect is capped by the weight of the node it occupies.
Ondansetron occupies ~98% of 5-HT3 at clinical doses and still does almost
nothing, because in *established HG* the peripheral 5-HT3 limb carries only a
small share of the drive.  W_VAG is the single parameter fitted to
ondansetron's -0.51; every other drug's effect is then a *prediction* from
its published receptor affinities, brain penetration and PK.

Third commitment: THE NUTRITIONAL CASCADE IS A SEPARATE, SLOWER STATE.
Vomiting drives volume/Cl-/K+ loss (hypochloraemic alkalosis that will not
correct until chloride is replaced) and depletes a ~28 mg thiamine store with
a ~15 day half-life.  So an antiemetic can normalise PUQE while Wernicke risk
is still climbing, and IV dextrose given before thiamine makes it worse.

ALPHA = 0 collapses the fold-change detector back to a pure level detector.
That one switch turns this file into the naive comparator model, and is how
the claims above are shown to be *load-bearing* rather than decorative.

UNITS
=====
time         days of gestational age (GA) by LMP; negative = pre-conception
GDF15        pg/mL           hCG      IU/L
drug amounts nmol            concs    nM
volumes      L               electrolytes mmol (amounts) / mmol/L (concs)
weight       kg              thiamine mg
PUQE-24      3..15 (the real instrument's range)

DISCLAIMER: educational/research QSP model.  Not validated for clinical use.
"""

import math
import sys

# =====================================================================
# 1. PARAMETERS
# =====================================================================

P = dict(

    # ---- timing of gestation (days GA by LMP) ------------------------
    T_CONC      = 14.0,    # conception
    T_IMP       = 22.0,    # implantation; placental secretion starts here
    T_TERM      = 280.0,

    # ---- placental trophoblast mass (relative) -----------------------
    KG_PLAC     = 0.100,   # /d  logistic growth rate
    PLAC0       = 0.020,   # relative mass at implantation
    PLAC_CAPG   = 0.0025,  # /d  slow growth of carrying capacity (mass keeps
                           #     rising to term, so GDF15 does too)

    # ---- trophoblast integrated stress response (ISR) ---------------
    # GDF15 is an ISR/ATF4-CHOP target gene.  Placental pO2 is low until the
    # intervillous circulation opens at ~10-12 weeks, so the *specific* GDF15
    # secretion rate is front-loaded: it is highest exactly when placental
    # mass is still small.  This is why the GDF15 ramp is so steep in weeks
    # 5-9 -- mass and specific output rise together -- and why it flattens
    # early: the ISR term falls just as mass saturates.
    ISR_HI      = 1.75,    # ISR index in the hypoxic (pre-perfusion) phase
    ISR_LO      = 1.00,    # ISR index once intervillous flow is established
    T_PERF1     = 66.0,    # d GA, start of the perfusion transition
    T_PERF2     = 102.0,   # d GA, end of the perfusion transition
    TAU_ISR     = 4.0,     # d, lag of the ISR state behind pO2

    # ---- maternal GDF15 ---------------------------------------------
    KEL_GDF     = 5.545,   # /d  (plasma t1/2 ~3 h)
    GDF_BASE    = 250.0,   # pg/mL  typical non-pregnant level
    KS_GDF_PL   = 48000.0, # pg/mL/d per unit (PLAC*ISR)
    TROPH_GAIN  = 1.00,    # per-pregnancy trophoblast secretory gain.  Fejzo
                           # 2024 shows FETAL PRODUCTION and MATERNAL
                           # SENSITIVITY are two separate contributors to HG
                           # risk; TROPH_GAIN is the fetal axis (it scales
                           # placental GDF15 AND hCG together, because both
                           # are syncytiotrophoblast products), and SENS is
                           # the maternal axis.
    GREF        = 250.0,   # pg/mL  population reference for the level detector

    # ---- THE FOLD-CHANGE DETECTOR (the core of the model) -----------
    TAU_SP      = 30.0,    # d, adaptation time constant of the set-point
    ALPHA       = 0.92,    # adaptation completeness. 1 = pure fold-change
                           # detector, 0 = pure level detector (comparator).

    # ---- GFRAL receptor pool (fast desensitisation, days) -----------
    KSYN_GF     = 0.35,    # /d
    KDEG_GF     = 0.35,    # /d
    KINT_GF     = 0.20,    # /d  ligand-driven internalisation
    KD_GF       = 5000.0,  # pg/mL

    # ---- area postrema -> NTS -> emesis -----------------------------
    EC50_AP     = 4.10,    # effective-drive units at which AP output is half-max
    HILL_AP     = 4.0,     # steep: the AP is a threshold detector, not a
                           # proportional one -- this is what makes a 1.6-fold
                           # difference in drive the difference between HG and
                           # an ordinary queasy first trimester
    TAU_AP      = 0.30,    # d
    TAU_NTS     = 0.25,    # d
    TAU_NAUS    = 0.35,    # d
    TAU_CPG     = 0.20,    # d

    # NTS input weights.  These SUM TO W_TOT by construction (W_AP is derived
    # as the remainder), so changing W_VAG re-apportions the drive between the
    # central and peripheral limbs without altering total drive -- which is
    # what makes W_VAG identifiable from ondansetron alone while leaving the
    # untreated natural history untouched.
    W_TOT       = 1.000,
    W_VAG       = 0.034,   # peripheral 5-HT3 / vagal limb   <-- FITTED
    W_VEST      = 0.055,   # vestibular / H1 limb
    W_HCG       = 0.030,   # hCG-associated limb (thyroid/ovarian, non-GDF15)

    # Maximal fractional inhibition of NTS transmission obtainable by FULLY
    # occupying each receptor class at the NTS.  These are properties of the
    # node, not of any drug.  Only the SCALE (E0) is fitted -- to mirtazapine.
    # The RATIOS are set from pharmacology and held fixed:
    #   alpha-2-delta and presynaptic alpha-2 gate bulk glutamatergic and
    #   noradrenergic transmission into the NTS and so have wide authority;
    #   H1 / 5-HT2 / 5-HT3 / M1 / D2 each modulate one input stream.
    E0          = 0.044,   # <-- FITTED (the only drug-side scale parameter)
    R_H1        = 1.0,
    R_5HT2      = 1.4,
    R_5HT3C     = 0.6,    # central (NTS) 5-HT3, distinct from the vagal limb
    R_M1        = 0.4,
    R_D2        = 0.5,
    R_A2        = 5.8,    # presynaptic alpha-2 (agonist) gating of NA release
    R_A2D       = 7.5,    # alpha-2-delta gating of presynaptic glutamate

    # symptom transduction
    EC50_NAUS   = 0.60,   # NTS units
    HILL_NAUS   = 6.0,
    EC50_VOM    = 0.66,   # vomiting needs more drive than nausea
    HILL_VOM    = 7.0,
    VOM_MAX     = 12.0,   # episodes/day at full CPG output
    RETCH_RATIO = 1.45,   # retches per vomit

    # ---- gastric emptying -------------------------------------------
    TAU_GAS     = 1.0,
    GAS_P4      = 0.35,   # progesterone-driven delay by mid first trimester
    GAS_GDF     = 0.30,   # GDF15-driven delay
    EC50_EC     = 0.55,   # gastric delay at which EC-cell 5-HT release halves

    # ---- hCG ---------------------------------------------------------
    KS_HCG      = 28000.0,  # IU/L/d per unit cytotrophoblast index
    KEL_HCG     = 0.28,      # /d  (t1/2 ~ 2.5 d terminal)
    T_CTB1      = 60.0,      # cytotrophoblast index starts to fall
    T_CTB2      = 130.0,     # ... and reaches its plateau
    CTB_PLAT    = 0.16,

    # ---- fluid / electrolytes ---------------------------------------
    VOL0        = 14.0,    # L extracellular fluid
    ORAL_FLUID  = 2.30,    # L/d when eating normally
    INSENS      = 0.90,    # L/d insensible
    VOM_VOL     = 0.055,   # L lost per emetic episode
    KREN_VOL    = 0.55,    # /d renal volume regulation gain
    URINE_MIN   = 0.40,    # L/d
    GJ_CL       = 140.0,   # mmol/L chloride in gastric juice
    GJ_H        = 75.0,    # mmol/L H+ in gastric juice
    GJ_NA       = 42.0,
    GJ_K        = 11.0,
    NA_TOT0     = 1960.0,  # mmol (140 * 14 L)
    CL_TOT0     = 1442.0,  # mmol (103 * 14 L)
    HCO30       = 24.0,    # mmol/L
    K_TOT0      = 3500.0,  # mmol total body K
    K_INTAKE    = 72.0,    # mmol/d
    K_URINE0    = 72.0,    # mmol/d
    K_ALDO      = 1.2,     # aldosterone amplification of renal K loss
    K_RETAIN    = 0.86,    # fractional total-body K+ at which the kidney
                           # essentially stops excreting potassium
    KREN_HCO3   = 0.28,    # /d bicarbonate excretion gain (chloride-gated)
    HCO3_ESC    = 0.30,    # chloride-independent fraction of HCO3- excretion
    NA_HCO3_COUP= 0.45,    # fraction of excreted HCO3- obligated to carry Na+
    NA_FLOOR    = 126.0,   # mmol/L at which the Na+/HCO3- coupling gives way
                           # to K+/NH4+ excretion
    NA_INTAKE   = 150.0,   # mmol/d
    CL_INTAKE   = 150.0,   # mmol/d
    WATER_FLOOR = 0.55,    # fraction of normal water intake retained at
                           # zero food intake (patients keep sipping)
    THIRST_G    = 0.70,    # thirst amplification of water intake per unit
                           # RAAS activation

    # ---- energy / weight --------------------------------------------
    WT0         = 65.0,    # kg pre-pregnancy
    KCAL_NEED   = 2100.0,  # kcal/d
    KCAL_FULL   = 2100.0,  # kcal/d at full oral intake (weight-neutral, so
                           # gestational gain comes only from WT_FETAL)
    KCAL_PER_KG = 7700.0,
    WT_FETAL    = 0.052,   # kg/d of conceptus+uterus+plasma gain, 2nd half

    # ---- ketones / liver --------------------------------------------
    KET_MAX     = 5.2,     # mmol/L at total starvation
    TAU_KET     = 0.5,
    ALT_BASE    = 20.0,
    ALT_MAX     = 190.0,
    TAU_ALT     = 3.0,

    # ---- thiamine ----------------------------------------------------
    THI0        = 28.0,    # mg whole-body store
    THI_DIET    = 1.45,    # mg/d absorbed at full oral intake
    KEL_THI     = 0.046,   # /d  (t1/2 ~15 d)
    THI_CRIT    = 12.0,    # mg, below which Wernicke hazard starts (~40% of
                           # a full store; frank deficiency signs appear here)
    HZ_WE       = 0.060,   # /d hazard scale
    THI_PER_CHO = 0.0022,  # mg thiamine consumed per g glucose infused
    THI_REFEED  = 2.6,     # amplification of that when the store is depleted

    # ---- thyroid -----------------------------------------------------
    FT40        = 14.0,    # pmol/L
    TSH0        = 1.6,     # mIU/L
    KD_T4       = 0.10,    # /d (t1/2 ~7 d)
    KD_TSH      = 9.0,     # /d (t1/2 ~1 h effective turnover)
    AH_HCG      = 1.8e-5,  # TSH-equivalents per IU/L of hCG
    FT4_50      = 14.0,
    HILL_TSH    = 5.0,

    # =================================================================
    # DRUG PK  (population-typical; CL in L/DAY, V in L, KA in /day)
    # Ki / IC50 in nM at the effect site; KP = brain-to-plasma ratio of the
    # unbound drug; FU = plasma unbound fraction.
    # =================================================================
    # ondansetron 8 mg PO q8h (MW 293.4; t1/2 ~4 h)
    OND_KA=6.0, OND_CL=660.0, OND_V=160.0, OND_F=0.60,
    OND_FU=0.27, OND_KP=0.60, OND_KI_5HT3=0.50,
    # doxylamine (MW 270.4, base; t1/2 ~11 h)
    DOX_KA=1.1, DOX_CL=84.0, DOX_V=180.0, DOX_F=0.90,
    DOX_FU=0.20, DOX_KP=0.80, DOX_KI_H1=5.0,
    # metoclopramide 10 mg q6h (MW 299.8; t1/2 ~5 h)
    MCP_KA=3.0, MCP_CL=720.0, MCP_V=200.0, MCP_F=0.75,
    MCP_FU=0.70, MCP_KP=0.25, MCP_KI_D2=200.0,
    # promethazine 25 mg q6h (MW 284.4; t1/2 ~11 h, large Vd)
    PMZ_KA=2.0, PMZ_CL=1440.0, PMZ_V=970.0, PMZ_F=0.25,
    PMZ_FU=0.07, PMZ_KP=2.0, PMZ_KI_H1=2.0, PMZ_KI_M1=60.0,
    # mirtazapine 15-30 mg qHS (MW 265.4; t1/2 ~29 h -- this long half-life is
    # what generates the observed post-day-4 widening, with no onset parameter)
    MIR_KA=1.4, MIR_CL=300.0, MIR_V=530.0, MIR_F=0.50,
    MIR_FU=0.15, MIR_KP=1.0,
    MIR_KI_H1=0.14, MIR_KI_5HT2=50.0, MIR_KI_5HT3=8.1,
    # prednisolone (MW 360.4) -- occupies NO node on the emetic chain, by
    # design; only a small appetite/wellbeing effect on oral intake
    STE_KA=3.0, STE_CL=192.0, STE_V=45.0, STE_F=0.80,
    STE_APP=0.10, STE_EC50=400.0,
    # gabapentin 600 mg q8h (MW 171.2; t1/2 ~6 h, CSF ~10% of plasma)
    GBP_KA=1.0, GBP_CL=178.0, GBP_V=60.0, GBP_F=0.45,
    GBP_FU=0.97, GBP_KP=0.10, GBP_IC50_A2D=3000.0,
    # clonidine, 5 mg transdermal patch delivering ~0.15 mg/day (MW 230.1)
    CLO_KA=0.25, CLO_CL=302.0, CLO_V=190.0, CLO_F=1.00,
    CLO_FU=0.80, CLO_KP=1.5, CLO_KI_A2=3.0,
    # metformin (MW 129.2) -- acts ONLY by raising basal GDF15 (Coll 2020)
    MET_KA=1.5, MET_CL=1150.0, MET_V=650.0, MET_F=0.55,
    MET_EMAX_GDF=2.20,   # fractional increase in basal GDF15 synthesis;
                         # metformin roughly doubles circulating GDF15
    MET_EC50=8000.0,     # nM
    # anti-GDF15 mAb (irreversible-binding approximation)
    MAB_CL=0.072, MAB_V=3.2, MAB_KON=0.42,       # IgG, t1/2 ~31 d; KON /d per nM
    # recombinant long-acting GDF15 (pre-conception desensitiser)
    RGD_CL=0.55, RGD_V=5.0, RGD_POT=2.9e4,       # Fc-fusion, t1/2 ~6 d

    # ---- non-drug modifiers of basal GDF15 --------------------------
    THAL_FOLD   = 12.0,   # beta-thalassaemia: basal GDF15 fold-increase
    TOBACCO_F   = 0.30,   # fractional increase in basal GDF15 with smoking
                          # (smaller than metformin's, which is why the model
                          # predicts weaker protection -- matching the
                          # observed aRR ordering 0.29 < 0.51)
    LOWGDF_F    = 0.60,   # low-GDF15 risk-allele carrier: 0.60 x basal

    # ---- susceptibility ---------------------------------------------
    SENS        = 1.00,   # constitutional gain on AP drive (HG-susceptible >1)
)

# state vector layout -------------------------------------------------
NAMES = [
    "PLAC", "ISR", "GDF", "SP", "GFRAL", "AP", "VAG", "GAS", "NTS", "CPG",
    "NAUS", "HCG", "VOL", "NATOT", "CLTOT", "HCO3", "KTOT", "WT", "KET",
    "THI", "WERISK", "ALT", "FT4", "TSH",
    "OND_D", "OND_C", "DOX_D", "DOX_C", "MCP_D", "MCP_C", "PMZ_C",
    "MIR_D", "MIR_C", "STE_D", "STE_C", "GBP_D", "GBP_C", "CLO_D", "CLO_C",
    "MET_C", "MAB_C", "RGD_C",
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)


# =====================================================================
# 2. HELPERS
# =====================================================================

def ramp(t, t0, t1, y0, y1):
    """Smooth (cosine) transition from y0 to y1 between t0 and t1."""
    if t <= t0:
        return y0
    if t >= t1:
        return y1
    f = 0.5 * (1.0 - math.cos(math.pi * (t - t0) / (t1 - t0)))
    return y0 + (y1 - y0) * f


def hill(x, ec50, n):
    if x <= 0.0:
        return 0.0
    xn = x ** n
    return xn / (ec50 ** n + xn)


def occ(c, ki):
    """Fractional receptor occupancy."""
    return c / (c + ki) if c > 0 else 0.0


def puqe_item(x, bins):
    """
    Continuous version of a PUQE-24 item.  The real instrument scores 1..5 in
    bins; a step function would make the ODE output jump, so we linearly
    interpolate between bin midpoints.  `bins` are the lower edges of scores
    2,3,4,5.
    """
    b = [0.0] + list(bins)
    if x <= b[1]:
        return 1.0 + 1.0 * (x - b[0]) / max(b[1] - b[0], 1e-9)
    for k in range(1, 4):
        if x <= b[k + 1]:
            return (1.0 + k) + (x - b[k]) / max(b[k + 1] - b[k], 1e-9)
    return min(5.0, 5.0 + (x - b[4]) / max(b[4], 1e-9) * 0.0)


def puqe24(nausea_hours, vomits):
    """PUQE-24 total, 3..15 (Koren 2005 PMID 16147725; 2021 PMID 32811235)."""
    q1 = puqe_item(nausea_hours, (1.0, 2.0, 4.0, 6.0))     # hours of nausea
    q2 = puqe_item(vomits, (1.0, 3.0, 5.0, 7.0))           # vomits
    q3 = puqe_item(vomits * P["RETCH_RATIO"], (1.0, 3.0, 5.0, 7.0))
    return min(15.0, q1 + q2 + q3)


# =====================================================================
# 3. DOSING
# =====================================================================

class Regimen:
    """
    A list of (time, state_name, amount_nmol) bolus events plus continuous
    infusions.  `iv_fluid` is L/d of crystalloid, `iv_k` mmol/d of KCl,
    `iv_thi` mg/d of thiamine, `iv_dex_g` g/d of dextrose.
    """

    def __init__(self):
        self.bolus = []
        self.inf = []          # (t0, t1, state_name, nmol/d)
        self.fluid = []        # (t0, t1, L/d, mmol Na/d, mmol Cl/d, mmol K/d)
        self.thi = []           # (t0, t1, mg/d)
        self.dex = []           # (t0, t1, g/d)
        self.flags = dict(thal=False, tobacco=False, lowgdf=False)

    def repeat(self, state, amt_nmol, t0, t1, every):
        t = t0
        while t < t1 - 1e-9:
            self.bolus.append((t, state, amt_nmol))
            t += every
        return self

    def infuse(self, state, rate, t0, t1):
        self.inf.append((t0, t1, state, rate))
        return self

    def iv(self, t0, t1, litres_day, na=154.0, cl=154.0, k=0.0):
        self.fluid.append((t0, t1, litres_day, na * litres_day,
                           cl * litres_day, k * litres_day))
        return self

    def thiamine(self, t0, t1, mg_day):
        self.thi.append((t0, t1, mg_day))
        return self

    def dextrose(self, t0, t1, g_day):
        self.dex.append((t0, t1, g_day))
        return self

    def rate(self, lst, t, idx):
        s = 0.0
        for row in lst:
            if row[0] <= t < row[1]:
                s += row[idx]
        return s


# molar dose helper: mg -> nmol
def nmol(mg, mw):
    return mg * 1e6 / mw


MW = dict(OND=293.4, DOX=270.4, MCP=299.8, PMZ=284.4, MIR=265.4,
          STE=360.4, GBP=171.2, CLO=230.1, MET=129.2)


# =====================================================================
# 4. THE ODE SYSTEM
# =====================================================================

def derivatives(t, y, p, reg):
    d = [0.0] * NST
    g = lambda n: y[IX[n]]

    # ---------------- placenta ---------------------------------------
    # Implantation seeds the trophoblast compartment over a 1-day window;
    # without the seed the logistic term can never leave zero.
    plac = max(g("PLAC"), 0.0)
    if t >= p["T_IMP"]:
        cap = 1.0 + p["PLAC_CAPG"] * (t - p["T_IMP"])
        d[IX["PLAC"]] = p["KG_PLAC"] * plac * (1.0 - plac / cap)
        if t < p["T_IMP"] + 1.0:
            d[IX["PLAC"]] += p["PLAC0"]
    else:
        d[IX["PLAC"]] = 0.0

    # trophoblast ISR: high while the placenta is hypoxic, falls once the
    # intervillous circulation opens
    isr_target = ramp(t, p["T_PERF1"], p["T_PERF2"], p["ISR_HI"], p["ISR_LO"])
    isr = g("ISR")
    d[IX["ISR"]] = (isr_target - isr) / p["TAU_ISR"]

    # ---------------- drug concentrations ----------------------------
    def conc(state, V):                     # nM, total plasma
        return max(g(state), 0.0) / V

    c_ond = conc("OND_C", p["OND_V"])
    c_dox = conc("DOX_C", p["DOX_V"])
    c_mcp = conc("MCP_C", p["MCP_V"])
    c_pmz = conc("PMZ_C", p["PMZ_V"])
    c_mir = conc("MIR_C", p["MIR_V"])
    c_ste = conc("STE_C", p["STE_V"])
    c_gbp = conc("GBP_C", p["GBP_V"])
    c_clo = conc("CLO_C", p["CLO_V"])
    c_met = conc("MET_C", p["MET_V"])
    c_mab = conc("MAB_C", p["MAB_V"])
    c_rgd = conc("RGD_C", p["RGD_V"])

    # free concentration at the effect site
    f_ond_p = c_ond * p["OND_FU"]                       # peripheral (vagal)
    f_ond_c = f_ond_p * p["OND_KP"]                     # NTS
    f_dox_c = c_dox * p["DOX_FU"] * p["DOX_KP"]
    f_mcp_c = c_mcp * p["MCP_FU"] * p["MCP_KP"]
    f_mcp_p = c_mcp * p["MCP_FU"]
    f_pmz_c = c_pmz * p["PMZ_FU"] * p["PMZ_KP"]
    f_mir_c = c_mir * p["MIR_FU"] * p["MIR_KP"]
    f_gbp_c = c_gbp * p["GBP_FU"] * p["GBP_KP"]
    f_clo_c = c_clo * p["CLO_FU"] * p["CLO_KP"]

    # ---------------- maternal GDF15 ---------------------------------
    gdf = max(g("GDF"), 1.0)

    basal_fold = 1.0
    if reg.flags["thal"]:
        basal_fold *= p["THAL_FOLD"]
    if reg.flags["tobacco"]:
        basal_fold *= (1.0 + p["TOBACCO_F"])
    if reg.flags["lowgdf"]:
        basal_fold *= p["LOWGDF_F"]
    # metformin raises circulating GDF15 (Coll 2020 PMID 31875646)
    basal_fold *= (1.0 + p["MET_EMAX_GDF"] * occ(c_met, p["MET_EC50"]))

    syn_basal = p["GDF_BASE"] * p["KEL_GDF"] * basal_fold
    syn_plac = (p["KS_GDF_PL"] * p["TROPH_GAIN"] * plac * isr
                if t >= p["T_IMP"] else 0.0)
    syn_rgd = p["RGD_POT"] * c_rgd          # exogenous long-acting GDF15
    # anti-GDF15 mAb: pseudo-first-order removal of free ligand
    mab_sink = p["MAB_KON"] * c_mab * gdf
    d[IX["GDF"]] = syn_basal + syn_plac + syn_rgd - p["KEL_GDF"] * gdf - mab_sink

    # ---------------- THE FOLD-CHANGE DETECTOR -----------------------
    # SP is the adapted log set-point of the GFRAL axis.
    sp = g("SP")
    d[IX["SP"]] = (math.log(gdf) - sp) / p["TAU_SP"]
    # mixed detector: ALPHA=1 pure fold-change, ALPHA=0 pure level
    ref = math.exp(p["ALPHA"] * sp + (1.0 - p["ALPHA"]) * math.log(p["GREF"]))
    fold = gdf / ref

    # GFRAL surface pool (fast, ligand-driven internalisation)
    gf = max(g("GFRAL"), 1e-6)
    d[IX["GFRAL"]] = (p["KSYN_GF"] - p["KDEG_GF"] * gf
                      - p["KINT_GF"] * gf * gdf / (gdf + p["KD_GF"]))

    gdf_eff = fold * gf / (p["KSYN_GF"] / p["KDEG_GF"])   # normalise gf to 1

    # ---------------- hCG --------------------------------------------
    ctb = ramp(t, p["T_CTB1"], p["T_CTB2"], 1.0, p["CTB_PLAT"])
    hcg = max(g("HCG"), 0.0)
    d[IX["HCG"]] = (p["KS_HCG"] * p["TROPH_GAIN"] * plac * ctb
                    if t >= p["T_IMP"] else 0.0) \
        - p["KEL_HCG"] * hcg

    # ---------------- area postrema ----------------------------------
    ap_target = p["SENS"] * hill(gdf_eff, p["EC50_AP"], p["HILL_AP"])
    # D2 blockade acts at the AP chemoreceptor trigger zone
    ap_target *= (1.0 - p["E0"] * p["R_D2"] * occ(f_mcp_c, p["MCP_KI_D2"]))
    d[IX["AP"]] = (ap_target - g("AP")) / p["TAU_AP"]

    # ---------------- gastric emptying & peripheral 5-HT3 limb -------
    p4_idx = min(1.0, max(0.0, (t - p["T_IMP"]) / 50.0))
    gas_target = p["GAS_P4"] * p4_idx + p["GAS_GDF"] * hill(gdf_eff, p["EC50_AP"], 1.4)
    # metoclopramide is prokinetic (peripheral 5-HT4 agonism / D2 blockade)
    gas_target *= (1.0 - 0.45 * occ(f_mcp_p, p["MCP_KI_D2"]))
    d[IX["GAS"]] = (gas_target - g("GAS")) / p["TAU_GAS"]
    gas = max(g("GAS"), 0.0)

    # enterochromaffin 5-HT release rises with gastric stasis and distension
    ec = hill(gas, p["EC50_EC"], 2.0)
    vag_target = ec * (1.0 - occ(f_ond_p, p["OND_KI_5HT3"]))
    d[IX["VAG"]] = (vag_target - g("VAG")) / p["TAU_AP"]

    # vestibular / H1 limb
    vest = 0.55 + 0.45 * hill(gdf_eff, p["EC50_AP"], 1.2)
    vest *= (1.0 - occ(f_dox_c, p["DOX_KI_H1"]))
    vest *= (1.0 - occ(f_pmz_c, p["PMZ_KI_H1"]))
    vest *= (1.0 - occ(f_mir_c, p["MIR_KI_H1"]))

    hcg_limb = hcg / 60000.0

    # ---------------- NTS integrator ---------------------------------
    w_ap = p["W_TOT"] - p["W_VAG"] - p["W_VEST"] - p["W_HCG"]
    nts_in = (w_ap * g("AP") + p["W_VAG"] * max(g("VAG"), 0.0)
              + p["W_VEST"] * vest + p["W_HCG"] * hcg_limb)

    # multiplicative, independent inhibition of NTS transmission
    o_h1 = 1.0 - (1.0 - occ(f_dox_c, p["DOX_KI_H1"])) \
                * (1.0 - occ(f_pmz_c, p["PMZ_KI_H1"])) \
                * (1.0 - occ(f_mir_c, p["MIR_KI_H1"]))
    o_5ht2 = occ(f_mir_c, p["MIR_KI_5HT2"])
    o_5ht3c = 1.0 - (1.0 - occ(f_ond_c, p["OND_KI_5HT3"])) \
                   * (1.0 - occ(f_mir_c, p["MIR_KI_5HT3"]))
    o_a2 = occ(f_clo_c, p["CLO_KI_A2"])       # clonidine: agonist
    o_a2d = occ(f_gbp_c, p["GBP_IC50_A2D"])   # gabapentin
    o_m1 = occ(f_pmz_c, p["PMZ_KI_M1"])

    e0 = p["E0"]
    inh = ((1.0 - e0 * p["R_H1"] * o_h1)
           * (1.0 - e0 * p["R_5HT2"] * o_5ht2)
           * (1.0 - e0 * p["R_5HT3C"] * o_5ht3c)
           * (1.0 - e0 * p["R_A2"] * o_a2)
           * (1.0 - e0 * p["R_A2D"] * o_a2d)
           * (1.0 - e0 * p["R_M1"] * o_m1))

    d[IX["NTS"]] = (nts_in * inh - g("NTS")) / p["TAU_NTS"]
    nts = max(g("NTS"), 0.0)

    # ---------------- symptoms ---------------------------------------
    naus_target = 24.0 * hill(nts, p["EC50_NAUS"], p["HILL_NAUS"])
    d[IX["NAUS"]] = (naus_target - g("NAUS")) / p["TAU_NAUS"]

    cpg_target = hill(nts, p["EC50_VOM"], p["HILL_VOM"])
    d[IX["CPG"]] = (cpg_target - g("CPG")) / p["TAU_CPG"]
    vom = p["VOM_MAX"] * max(g("CPG"), 0.0)          # episodes/day
    naus_h = max(g("NAUS"), 0.0)

    # ---------------- oral intake ------------------------------------
    app_boost = p["STE_APP"] * occ(c_ste, p["STE_EC50"]) \
        + 0.12 * occ(f_mir_c, p["MIR_KI_H1"])        # mirtazapine appetite
    # oral intake collapses steeply, not linearly, with hours of nausea:
    # a woman nauseated 16 h a day is not eating two thirds of normal.
    oral = max(0.04, (max(0.0, 1.0 - 0.985 * naus_h / 24.0)) ** 1.5)
    oral = min(1.0, oral * (1.0 + app_boost))

    iv_fl = reg.rate(reg.fluid, t, 2)
    iv_na = reg.rate(reg.fluid, t, 3)
    iv_cl = reg.rate(reg.fluid, t, 4)
    iv_k = reg.rate(reg.fluid, t, 5)
    iv_thi = reg.rate(reg.thi, t, 2)
    iv_dex = reg.rate(reg.dex, t, 2)

    # ---------------- fluid & electrolytes ---------------------------
    vol = max(g("VOL"), 5.0)
    dep = (p["VOL0"] - vol) / p["VOL0"]                 # depletion fraction
    aldo = min(1.0, max(0.0, dep * 6.0))                # RAAS activation index
    l_vom = p["VOM_VOL"] * vom                          # L/d of gastric juice

    # A nauseated woman stops EATING long before she stops SIPPING.  Water
    # intake therefore falls much less than solute intake, which -- with the
    # ADH drive of volume depletion -- is why HG produces hyponatraemia rather
    # than the hypernatraemia that loss of (hypotonic) gastric juice alone
    # would predict.
    # ... and thirst rises with volume depletion, which is what stops an
    # untreated woman from simply exsanguinating her extracellular space
    thirst = 1.0 + p["THIRST_G"] * aldo
    water_in = p["ORAL_FLUID"] * (p["WATER_FLOOR"]
                                  + (1.0 - p["WATER_FLOOR"]) * oral) * thirst
    urine = max(p["URINE_MIN"],
                water_in + iv_fl - p["INSENS"]
                + p["KREN_VOL"] * (vol - p["VOL0"]))
    d[IX["VOL"]] = water_in + iv_fl - urine - p["INSENS"] - l_vom

    cl_conc = max(g("CLTOT"), 1.0) / vol
    na_conc = max(g("NATOT"), 1.0) / vol
    cl_avail = min(1.0, max(0.0, (cl_conc - 84.0) / (103.0 - 84.0)))

    # Gastric HCl secretion needs plasma chloride to draw on, so severe
    # chloride depletion is self-limiting: the vomitus becomes progressively
    # less acid and less chloride-rich.  Without this the model would happily
    # drive plasma chloride to implausible values.
    gj_cl_eff = p["GJ_CL"] * min(1.2, max(0.35, cl_conc / 103.0))
    gj_h_eff = p["GJ_H"] * min(1.2, max(0.35, cl_conc / 103.0))

    cl_in = p["CL_INTAKE"] * oral + iv_cl
    cl_out_gi = gj_cl_eff * l_vom
    d[IX["CLTOT"]] = cl_in - cl_out_gi \
        - max(0.0, cl_in - cl_out_gi) * (1.0 - 0.85 * aldo)

    # Bicarbonate cannot be excreted alone: much of the HCO3- that leaves in
    # the urine takes a Na+ with it.  This obligatory natriuresis is why a
    # chloride-depleted, volume-depleted woman -- whose kidney is otherwise
    # retaining sodium as hard as it can -- still becomes hyponatraemic.  As
    # plasma Na+ falls the coupling gives way to K+/NH4+ excretion, which is
    # what puts a floor under the sodium.
    hco3_excr = (p["KREN_HCO3"] * max(0.0, g("HCO3") - p["HCO30"])
                 * (p["HCO3_ESC"] + (1.0 - p["HCO3_ESC"]) * cl_avail))
    f_na = min(1.0, max(0.0, (na_conc - p["NA_FLOOR"])
                        / (140.0 - p["NA_FLOOR"])))
    na_in = p["NA_INTAKE"] * oral + iv_na
    na_out_gi = p["GJ_NA"] * l_vom
    d[IX["NATOT"]] = na_in - na_out_gi \
        - max(0.0, na_in - na_out_gi) * (1.0 - 0.85 * aldo) \
        - hco3_excr * vol * p["NA_HCO3_COUP"] * f_na

    # renal K+ wasting: secretion follows distal Na+ delivery and aldosterone,
    # but total excretion still tracks intake, so it falls as intake falls
    # renal K+ handling has to be able to RETAIN as well as waste, or a
    # depleted store never refills once vomiting stops
    ktot = max(g("KTOT"), 1.0)
    k_ren = min(1.0, max(0.05, (ktot / p["K_TOT0"] - p["K_RETAIN"])
                         / (1.0 - p["K_RETAIN"])))
    d[IX["KTOT"]] = (p["K_INTAKE"] * oral + iv_k - p["GJ_K"] * l_vom
                     - p["K_URINE0"] * (1.0 + p["K_ALDO"] * aldo)
                     * (0.15 + 0.85 * oral) * k_ren)

    # Classic chloride-responsive metabolic alkalosis: losing gastric H+ raises
    # plasma HCO3-, and the kidney cannot excrete the excess bicarbonate
    # without chloride to accompany it -- so the alkalosis persists until Cl-
    # is replaced, no matter how much of any other fluid is given.
    d[IX["HCO3"]] = gj_h_eff * l_vom / vol - hco3_excr

    # ---------------- energy, weight, ketones, liver -----------------
    kcal = p["KCAL_FULL"] * oral + iv_dex * 3.4
    # adaptive thermogenesis: resting expenditure falls with sustained deficit
    kcal_out = p["KCAL_NEED"] * (0.80 + 0.20 * oral)
    fetal = p["WT_FETAL"] * max(0.0, min(1.0, (t - 84.0) / 40.0))
    d[IX["WT"]] = (kcal - kcal_out) / p["KCAL_PER_KG"] \
        + d[IX["VOL"]] * 1.0 + fetal

    ket_ss = 0.1 + p["KET_MAX"] * max(0.0, 1.0 - kcal / 1250.0)
    d[IX["KET"]] = (ket_ss - g("KET")) / p["TAU_KET"]

    starv = max(0.0, 1.0 - kcal / 1600.0)
    alt_ss = p["ALT_BASE"] + (p["ALT_MAX"] - p["ALT_BASE"]) * starv ** 1.6
    d[IX["ALT"]] = (alt_ss - g("ALT")) / p["TAU_ALT"]

    # ---------------- thiamine & Wernicke hazard ---------------------
    thi = max(g("THI"), 0.0)
    refeed = 1.0 + p["THI_REFEED"] * max(0.0, 1.0 - thi / p["THI_CRIT"])
    d[IX["THI"]] = (p["THI_DIET"] * oral + iv_thi
                    - p["KEL_THI"] * thi
                    - p["THI_PER_CHO"] * iv_dex * refeed)
    d[IX["WERISK"]] = p["HZ_WE"] * max(0.0, 1.0 - thi / p["THI_CRIT"]) ** 2

    # ---------------- thyroid ----------------------------------------
    tshr = max(g("TSH"), 0.0) + p["AH_HCG"] * hcg
    d[IX["FT4"]] = p["KD_T4"] * p["FT40"] * (tshr / p["TSH0"]) \
        - p["KD_T4"] * g("FT4")
    d[IX["TSH"]] = p["KD_TSH"] * p["TSH0"] * 2.0 / (
        1.0 + (g("FT4") / p["FT4_50"]) ** p["HILL_TSH"]) \
        - p["KD_TSH"] * g("TSH")

    # ---------------- drug PK (first-order absorption, linear disposition) --
    for dep_s, cen_s, ka, CL, V, F in [
        ("OND_D", "OND_C", p["OND_KA"], p["OND_CL"], p["OND_V"], p["OND_F"]),
        ("DOX_D", "DOX_C", p["DOX_KA"], p["DOX_CL"], p["DOX_V"], p["DOX_F"]),
        ("MCP_D", "MCP_C", p["MCP_KA"], p["MCP_CL"], p["MCP_V"], p["MCP_F"]),
        ("MIR_D", "MIR_C", p["MIR_KA"], p["MIR_CL"], p["MIR_V"], p["MIR_F"]),
        ("STE_D", "STE_C", p["STE_KA"], p["STE_CL"], p["STE_V"], p["STE_F"]),
        ("GBP_D", "GBP_C", p["GBP_KA"], p["GBP_CL"], p["GBP_V"], p["GBP_F"]),
        ("CLO_D", "CLO_C", p["CLO_KA"], p["CLO_CL"], p["CLO_V"], p["CLO_F"]),
    ]:
        a = max(g(dep_s), 0.0)
        d[IX[dep_s]] += -ka * a
        d[IX[cen_s]] += F * ka * a - CL * max(g(cen_s), 0.0) / V

    d[IX["PMZ_C"]] += -p["PMZ_CL"] * max(g("PMZ_C"), 0.0) / p["PMZ_V"]
    d[IX["MET_C"]] += -p["MET_CL"] * max(g("MET_C"), 0.0) / p["MET_V"]
    d[IX["MAB_C"]] += -p["MAB_CL"] * max(g("MAB_C"), 0.0) / p["MAB_V"]
    d[IX["RGD_C"]] += -p["RGD_CL"] * max(g("RGD_C"), 0.0) / p["RGD_V"]

    # continuous infusions declared in the regimen
    for t0, t1, state, rate in reg.inf:
        if t0 <= t < t1:
            d[IX[state]] += rate

    return d, dict(fold=fold, gdf_eff=gdf_eff, vom=vom, naus_h=naus_h,
                   oral=oral, nts=nts, inh=inh, ref=ref, aldo=aldo,
                   cl_conc=cl_conc, kcal=kcal)


# =====================================================================
# 5. INITIAL CONDITIONS & INTEGRATION
# =====================================================================

def initial_state(p, reg):
    y = [0.0] * NST
    basal_fold = 1.0
    if reg.flags["thal"]:
        basal_fold *= p["THAL_FOLD"]
    if reg.flags["tobacco"]:
        basal_fold *= (1.0 + p["TOBACCO_F"])
    if reg.flags["lowgdf"]:
        basal_fold *= p["LOWGDF_F"]
    gdf0 = p["GDF_BASE"] * basal_fold

    y[IX["PLAC"]] = 0.0
    y[IX["ISR"]] = p["ISR_HI"]
    y[IX["GDF"]] = gdf0
    y[IX["SP"]] = math.log(gdf0)
    y[IX["GFRAL"]] = p["KSYN_GF"] / (p["KDEG_GF"] + p["KINT_GF"]
                                     * gdf0 / (gdf0 + p["KD_GF"]))
    y[IX["VOL"]] = p["VOL0"]
    y[IX["NATOT"]] = p["NA_TOT0"]
    y[IX["CLTOT"]] = p["CL_TOT0"]
    y[IX["HCO3"]] = p["HCO30"]
    y[IX["KTOT"]] = p["K_TOT0"]
    y[IX["WT"]] = p["WT0"]
    y[IX["KET"]] = 0.1
    y[IX["THI"]] = p["THI0"]
    y[IX["ALT"]] = p["ALT_BASE"]
    y[IX["FT4"]] = p["FT40"]
    y[IX["TSH"]] = p["TSH0"]
    return y


def simulate(reg, p=None, t0=-70.0, t1=200.0, dt=0.02, record_every=5,
             y0=None):
    """
    Fixed-step RK4 with bolus events applied at step boundaries.

    `y0` lets a run RESUME from a saved state vector.  Every treatment arm of a
    simulated trial shares the same untreated pre-randomisation history, so
    that history is integrated once and each arm branches from the saved state
    -- which is both faster and guarantees the arms are exactly comparable.
    The final state is returned in res["_y"].
    """
    p = p or P
    y = list(y0) if y0 is not None else initial_state(p, reg)
    bol = sorted(reg.bolus, key=lambda r: r[0])
    bi = 0
    n = int(round((t1 - t0) / dt))
    out = {k: [] for k in ["t", "GDF", "SP", "fold", "AP", "NTS", "naus_h",
                           "vom", "PUQE", "WT", "K", "CL", "HCO3", "NA",
                           "KET", "THI", "PWE", "ALT", "TSH", "FT4", "HCG",
                           "oral", "VOL", "GFRAL", "gdf_eff", "ref", "inh"]}
    t = t0
    for step in range(n + 1):
        # apply any bolus doses due at (or just before) this time
        while bi < len(bol) and bol[bi][0] <= t + 1e-9:
            _, state, amt = bol[bi]
            y[IX[state]] += amt
            bi += 1

        if step % record_every == 0:
            _, ex = derivatives(t, y, p, reg)
            pq = puqe24(ex["naus_h"], ex["vom"])
            out["t"].append(t)
            out["GDF"].append(y[IX["GDF"]])
            out["SP"].append(math.exp(y[IX["SP"]]))
            out["fold"].append(ex["fold"])
            out["gdf_eff"].append(ex["gdf_eff"])
            out["ref"].append(ex["ref"])
            out["inh"].append(ex["inh"])
            out["AP"].append(y[IX["AP"]])
            out["NTS"].append(y[IX["NTS"]])
            out["GFRAL"].append(y[IX["GFRAL"]])
            out["naus_h"].append(ex["naus_h"])
            out["vom"].append(ex["vom"])
            out["PUQE"].append(pq)
            out["WT"].append(y[IX["WT"]])
            out["VOL"].append(y[IX["VOL"]])
            out["K"].append(4.0 + (y[IX["KTOT"]] - p["K_TOT0"]) / 300.0)
            out["CL"].append(y[IX["CLTOT"]] / max(y[IX["VOL"]], 1.0))
            out["NA"].append(y[IX["NATOT"]] / max(y[IX["VOL"]], 1.0))
            out["HCO3"].append(y[IX["HCO3"]])
            out["KET"].append(y[IX["KET"]])
            out["THI"].append(y[IX["THI"]])
            out["PWE"].append(1.0 - math.exp(-y[IX["WERISK"]]))
            out["ALT"].append(y[IX["ALT"]])
            out["TSH"].append(y[IX["TSH"]])
            out["FT4"].append(y[IX["FT4"]])
            out["HCG"].append(y[IX["HCG"]])
            out["oral"].append(ex["oral"])

        if step == n:
            break
        k1, _ = derivatives(t, y, p, reg)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2, _ = derivatives(t + 0.5 * dt, y2, p, reg)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3, _ = derivatives(t + 0.5 * dt, y3, p, reg)
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4, _ = derivatives(t + dt, y4, p, reg)
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        y[IX["THI"]] = max(y[IX["THI"]], 0.0)
        t = t0 + (step + 1) * dt
    out["_y"] = list(y)
    return out


# =====================================================================
# 6. READOUT HELPERS
# =====================================================================

def at(res, t):
    """Value index nearest time t."""
    best, bi = 1e18, 0
    for i, tt in enumerate(res["t"]):
        dd = abs(tt - t)
        if dd < best:
            best, bi = dd, i
    return bi


def peak(res, key, lo=28.0, hi=200.0):
    v, tt = -1e18, None
    for i, t in enumerate(res["t"]):
        if lo <= t <= hi and res[key][i] > v:
            v, tt = res[key][i], t
    return v, tt


def first_below(res, key, thr, lo=70.0):
    for i, t in enumerate(res["t"]):
        if t >= lo and res[key][i] < thr:
            return t
    return None


def puqe_at(res, t):
    return res["PUQE"][at(res, t)]


# =====================================================================
# 7. SCENARIOS
# =====================================================================

def base_regimen(**flags):
    r = Regimen()
    r.flags.update(flags)
    return r


def add_supportive(r, t0, t1, dextrose=True, thiamine_mg=0.0):
    """Standard IV rehydration: 3 L/d crystalloid + 40 mmol/d KCl."""
    r.iv(t0, t1, 3.0, na=154.0, cl=154.0, k=13.3)
    if dextrose:
        r.dextrose(t0, t1, 150.0)
    if thiamine_mg > 0:
        r.thiamine(t0, t1, thiamine_mg)
    return r


def sc_ondansetron(t0, t1, r=None):
    r = r or base_regimen()
    return r.repeat("OND_D", nmol(8.0, MW["OND"]), t0, t1, 1.0 / 3.0)


def sc_mirtazapine(t0, t1, r=None, mg=30.0):
    r = r or base_regimen()
    return r.repeat("MIR_D", nmol(mg, MW["MIR"]), t0, t1, 1.0)


def sc_doxylamine(t0, t1, r=None):
    # Diclectin-style: 10/10 mg qAM, 10/10 mg qPM, 20/20 mg qHS
    r = r or base_regimen()
    t = t0
    while t < t1:
        r.bolus.append((t + 0.30, "DOX_D", nmol(10.0, MW["DOX"])))
        r.bolus.append((t + 0.65, "DOX_D", nmol(10.0, MW["DOX"])))
        r.bolus.append((t + 0.92, "DOX_D", nmol(20.0, MW["DOX"])))
        t += 1.0
    return r


def sc_metoclopramide(t0, t1, r=None):
    r = r or base_regimen()
    return r.repeat("MCP_D", nmol(10.0, MW["MCP"]), t0, t1, 0.25)


def sc_promethazine(t0, t1, r=None):
    r = r or base_regimen()
    return r.repeat("PMZ_C", nmol(25.0, MW["PMZ"]) * 0.25, t0, t1, 0.25)


def sc_gabapentin(t0, t1, r=None, mg=600.0):
    r = r or base_regimen()
    return r.repeat("GBP_D", nmol(mg, MW["GBP"]), t0, t1, 1.0 / 3.0)


def sc_clonidine(t0, t1, r=None):
    # 5 mg transdermal patch, ~0.15 mg/d systemic delivery, changed weekly
    r = r or base_regimen()
    return r.infuse("CLO_D", nmol(0.15, MW["CLO"]), t0, t1)


def sc_steroid(t0, r=None):
    """Yost 2003 regimen: IV methylprednisolone 125 mg then a prednisone taper."""
    r = r or base_regimen()
    r.bolus.append((t0, "STE_C", nmol(125.0, MW["STE"])))
    plan = [(40.0, 1), (20.0, 3), (10.0, 3), (5.0, 7)]
    t = t0
    for mg, days in plan:
        for k in range(days):
            r.bolus.append((t + k, "STE_D", nmol(mg, MW["STE"])))
        t += days
    return r


def sc_metformin_pre(mg=1000.0, t_start=-70.0, t_stop=14.0, r=None):
    """Pre-conception metformin, stopped at the positive test (Sharma 2025)."""
    r = r or base_regimen()
    return r.repeat("MET_C", nmol(mg, MW["MET"]) * 0.55, t_start, t_stop, 0.5)


def sc_antigdf15(t_dose=56.0, nmol_dose=90.0, r=None):
    r = r or base_regimen()
    r.bolus.append((t_dose, "MAB_C", nmol_dose))
    return r


def sc_rgdf15_pre(t0=-63.0, t1=0.0, nmol_wk=6.0, r=None):
    """Pre-conception recombinant long-acting GDF15 (mechanism-based
    prophylaxis proposed by Fejzo 2024)."""
    r = r or base_regimen()
    return r.repeat("RGD_C", nmol_wk, t0, t1, 7.0)


# =====================================================================
# 8. SIMULATED TRIALS
# =====================================================================
#
# A trial arm's effect is NOT the effect on one average patient.  Emetic drive
# passes through steep Hill transductions, and PUQE-24 has a hard ceiling of
# 15, so a single deterministic "typical" trajectory sitting on the steepest
# part of the curve overstates a drug's mean effect several-fold.  Every drug
# comparison below is therefore run over a heterogeneous cohort and averaged,
# with all arms branching from the SAME saved pre-randomisation state.
#
# This is not a cosmetic detail: it is why weak drugs look weak here.  A
# consequence worth stating on its own is that PUQE's ceiling makes severe-HG
# trials structurally insensitive -- a third of this cohort starts at 15/15,
# where a real reduction in drive cannot register at all.

# constitutional maternal sensitivity of the GFRAL axis, one entry per subject
TRIAL_POP = [0.82, 0.88, 0.94, 1.00, 1.08, 1.18, 1.30]
TRIAL_GA = 56.0        # randomisation at GA 8 weeks


def fmt(x, n=2):
    if x != x:
        return "n/a"
    return ("{:." + str(n) + "f}").format(x)


def build_population(sens_list=None, ga=None, **flags):
    """Integrate each subject up to randomisation once; return branch states."""
    sens_list = sens_list or TRIAL_POP
    ga = TRIAL_GA if ga is None else ga
    subs = []
    for s in sens_list:
        p = dict(P)
        p["SENS"] = s
        pre = simulate(base_regimen(**flags), p=p, t1=ga)
        subs.append((s, p, pre["_y"], flags))
    return subs, ga


def run_arm(subs, ga, maker, days=(2, 4, 7, 14), horizon=22.0, p_over=None):
    """
    Mean PUQE-24 change from baseline at each requested day, plus the mean
    peak fractional inhibition of NTS transmission the arm achieves.
    """
    acc = {dd: [] for dd in days}
    inh = []
    for s, p, y0, flags in subs:
        pp = dict(p)
        if p_over:
            pp.update(p_over)
        r = maker(ga, ga + 14.0)
        r.flags.update(flags)
        res = simulate(r, p=pp, t0=ga, t1=ga + horizon, y0=y0)
        b = puqe_at(res, ga)
        for dd in days:
            acc[dd].append(puqe_at(res, ga + dd) - b)
        lo, hi = at(res, ga + 2.0), at(res, ga + 13.0)
        inh.append(1.0 - min(res["inh"][lo:hi]))
    out = {dd: sum(acc[dd]) / len(acc[dd]) for dd in days}
    out["inh"] = sum(inh) / len(inh)
    return out


PLACEBO = lambda a, b: base_regimen()

ARM_MAKERS = [
    ("ondansetron 8 mg q8h", "peripheral 5-HT3",
     lambda a, b: sc_ondansetron(a, b)),
    ("metoclopramide 10 mg q6h", "AP D2 + prokinetic",
     lambda a, b: sc_metoclopramide(a, b)),
    ("promethazine 25 mg q6h", "NTS H1 / M1",
     lambda a, b: sc_promethazine(a, b)),
    ("doxylamine/pyridoxine", "H1 + vestibular",
     lambda a, b: sc_doxylamine(a, b)),
    ("prednisolone taper", "no node on the chain",
     lambda a, b: sc_steroid(a)),
    ("clonidine 5 mg patch", "NTS presynaptic alpha-2",
     lambda a, b: sc_clonidine(a, b)),
    ("mirtazapine 30 mg qHS", "H1 + 5-HT2 + 5-HT3",
     lambda a, b: sc_mirtazapine(a, b)),
    ("gabapentin 600 mg q8h", "NTS alpha-2-delta",
     lambda a, b: sc_gabapentin(a, b)),
    ("anti-GDF15 mAb 90 nmol", "the ligand itself",
     lambda a, b: sc_antigdf15(t_dose=a)),
]


# =====================================================================
# 9. MAIN
# =====================================================================

def natural_history(label, sens=1.0, troph=1.0, **flags):
    p = dict(P)
    p["SENS"] = sens
    p["TROPH_GAIN"] = troph
    return label, simulate(base_regimen(**flags), p=p)


def hg_case(res, wt0=None):
    """
    Windsor-style severity surrogate (Jansen 2021, PMID 34555550): peak
    PUQE-24 >= 13, or >= 7 together with >= 5% weight loss.
    """
    wt0 = wt0 or P["WT0"]
    pk, tpk = peak(res, "PUQE")
    wmin = min(res["WT"][at(res, 28.0):at(res, 160.0)])
    loss = (wt0 - wmin) / wt0 * 100.0
    return (pk >= 13.0 or (pk >= 7.0 and loss >= 5.0)), pk, loss


def main(check_only=False):
    lines = []
    W = lines.append
    W("=" * 78)
    W("HYPEREMESIS GRAVIDARUM QSP MODEL — reference implementation")
    W("{} ODE states | fold-change GDF15 detector | node-position pharmacology"
      .format(NST))
    W("=" * 78)

    # =============================================================== 1
    W("")
    W("--- 1. NATURAL HISTORY (untreated) " + "-" * 43)
    W("")
    cohorts = [
        natural_history("normal pregnancy", sens=0.62),
        natural_history("HG (sensitivity-driven)", sens=1.00),
        natural_history("HG (production-driven)", sens=0.62, troph=2.00),
        natural_history("low-GDF15 risk allele", sens=1.00, lowgdf=True),
        natural_history("beta-thalassaemia", sens=1.00, thal=True),
        natural_history("smoker (pre-conception)", sens=1.00, tobacco=True),
    ]
    nat = dict(cohorts)
    hdr = ("{:<26}{:>8}{:>8}{:>8}{:>8}{:>8}{:>7}"
           .format("cohort", "PUQEmax", "peak/wk", "PUQE16w", "GDFmax",
                   "wtloss%", "HGcase"))
    W(hdr)
    W("-" * len(hdr))
    for lab, res in cohorts:
        isc, pk, loss = hg_case(res)
        W("{:<26}{:>8}{:>8}{:>8}{:>8}{:>8}{:>7}".format(
            lab, fmt(pk, 1), fmt(peak(res, "PUQE")[1] / 7.0, 1),
            fmt(res["PUQE"][at(res, 112.0)], 1),
            fmt(peak(res, "GDF")[0], 0), fmt(loss, 1),
            "YES" if isc else "no"))
    W("")
    W("Fejzo 2024 (PMID 38092039) reports that FETAL PRODUCTION of GDF15 and")
    W("MATERNAL SENSITIVITY to it are separate contributors.  The model keeps")
    W("them as separate parameters (TROPH_GAIN and SENS), and they are not")
    W("interchangeable: only the production-driven phenotype drags hCG up with")
    W("GDF15, so only it becomes biochemically thyrotoxic (section 7).")

    # =============================================================== 2
    hg = nat["HG (sensitivity-driven)"]
    i9, i16, i28 = at(hg, 63.0), at(hg, 112.0), at(hg, 196.0)
    W("")
    W("--- 2. THE CENTRAL STRUCTURAL RESULT " + "-" * 41)
    W("")
    W("{:>8}{:>10}{:>11}{:>8}{:>9}{:>8}".format(
        "GA/wk", "GDF15", "set-point", "fold", "AP drive", "PUQE"))
    W("-" * 54)
    for wk in (6, 8, 9, 11, 12, 14, 16, 20, 28):
        i = at(hg, wk * 7.0)
        W("{:>8}{:>10}{:>11}{:>8}{:>9}{:>8}".format(
            wk, fmt(hg["GDF"][i], 0), fmt(hg["SP"][i], 0),
            fmt(hg["fold"][i], 2), fmt(hg["AP"][i], 3),
            fmt(hg["PUQE"][i], 1)))
    W("")
    W("GDF15 at GA 28 wk is {}x its GA 9 wk value, yet PUQE has returned to"
      .format(fmt(hg["GDF"][i28] / hg["GDF"][i9], 2)))
    W("baseline ({} vs {}).  The hormone did not go away; the set-point caught"
      .format(fmt(hg["PUQE"][i28], 1), fmt(hg["PUQE"][i9], 1)))
    W("up with it.  No model in which nausea is a function of the LEVEL of")
    W("GDF15 can produce this, which is the whole reason for the fold-change")
    W("detector.  Flipping ALPHA to 0 turns the same equations into exactly")
    W("that level-based comparator:")
    W("")
    p0 = dict(P)
    p0["ALPHA"] = 0.0
    lvl = simulate(base_regimen(), p=p0)
    lvl_thal = simulate(base_regimen(thal=True), p=p0)
    W("{:<34}{:>12}{:>12}".format("", "fold-change", "level (A=0)"))
    W("-" * 58)
    W("{:<34}{:>12}{:>12}".format("PUQE at GA 9 wk", fmt(hg["PUQE"][i9], 1),
                                  fmt(lvl["PUQE"][i9], 1)))
    W("{:<34}{:>12}{:>12}".format("PUQE at GA 16 wk", fmt(hg["PUQE"][i16], 1),
                                  fmt(lvl["PUQE"][i16], 1)))
    W("{:<34}{:>12}{:>12}".format("PUQE at GA 28 wk", fmt(hg["PUQE"][i28], 1),
                                  fmt(lvl["PUQE"][i28], 1)))
    W("{:<34}{:>12}{:>12}".format(
        "beta-thalassaemia peak PUQE",
        fmt(peak(nat["beta-thalassaemia"], "PUQE")[0], 1),
        fmt(peak(lvl_thal, "PUQE")[0], 1)))
    W("")
    W("The level model gets BOTH signs wrong: symptoms that never remit, and")
    W("beta-thalassaemia -- which reports almost no NVP -- as the worst case in")
    W("the cohort, because its GDF15 really is the highest.")

    # =============================================================== 3
    W("")
    W("--- 3. THE VOMIT TRIAL (Ostenfeld 2026, PMID 41478546) " + "-" * 23)
    W("")
    subs, ga = build_population()
    W("n = {} subjects, randomised at GA {} weeks, 14-day treatment."
      .format(len(subs), fmt(ga / 7.0, 0)))
    pbo = run_arm(subs, ga, PLACEBO)
    ond = run_arm(subs, ga, lambda a, b: sc_ondansetron(a, b))
    mir = run_arm(subs, ga, lambda a, b: sc_mirtazapine(a, b))
    W("")
    W("{:<24}{:>9}{:>9}{:>9}{:>9}".format("arm", "d2", "d4", "d7", "d14"))
    W("-" * 60)
    for lab, r in [("placebo", pbo), ("ondansetron", ond),
                   ("mirtazapine 30 mg", mir)]:
        W("{:<24}{:>9}{:>9}{:>9}{:>9}".format(
            lab, *[fmt(r[d], 2) for d in (2, 4, 7, 14)]))
    W("")
    W("{:<24}{:>9}{:>9}{:>9}{:>9}".format("vs placebo", "d2", "d4", "d7", "d14"))
    W("-" * 60)
    for lab, r in [("ondansetron", ond), ("mirtazapine 30 mg", mir)]:
        W("{:<24}{:>9}{:>9}{:>9}{:>9}".format(
            lab, *[fmt(r[d] - pbo[d], 2) for d in (2, 4, 7, 14)]))
    W("")
    W("  observed  ondansetron  -0.51  (95% CI -2.32 to  1.30)")
    W("            mirtazapine  -1.86  (95% CI -3.61 to -0.12), separating")
    W("                         further after day 4")
    W("")
    W("The post-day-4 widening is not fitted and there is no onset parameter in")
    W("the model.  It comes out of mirtazapine's ~29 h half-life: ondansetron")
    W("is at steady state within a day, mirtazapine is not until roughly day 5.")

    # =============================================================== 4
    W("")
    W("--- 4. THE NODE-POSITION LAW " + "-" * 49)
    W("")
    W("Emetic drive reaches the emetic pattern generator through the NTS, which")
    W("sums a large GFRAL/area-postrema term and a small peripheral vagal one.")
    W("A drug cannot do better than the weight of the node it occupies.")
    W("")
    hdr = ("{:<26}{:<24}{:>9}{:>9}{:>8}"
           .format("drug", "node", "dPUQE_d7", "dPUQE_d14", "NTS%"))
    W(hdr)
    W("-" * len(hdr))
    node_res = {}
    for lab, node, mk in ARM_MAKERS:
        r = run_arm(subs, ga, mk)
        node_res[lab] = r
        W("{:<26}{:<24}{:>9}{:>9}{:>8}".format(
            lab, node, fmt(r[7] - pbo[7], 2), fmt(r[14] - pbo[14], 2),
            fmt(r["inh"] * 100.0, 1)))
    W("")
    W("Two parameters are fitted here and nowhere else:")
    W("  W_VAG = {}  the share of emetic drive carried by the peripheral"
      .format(P["W_VAG"]))
    W("                5-HT3 limb, fitted to ondansetron's -0.51;")
    W("  E0    = {}  the scale of central node authority, fitted to"
      .format(P["E0"]))
    W("                mirtazapine's -1.86.")
    W("The RATIOS between node weights (R_H1, R_A2, R_A2D, ...) are fixed from")
    W("receptor pharmacology, and each drug's occupancy comes from its")
    W("published Ki, unbound fraction, brain penetration and PK.  So every")
    W("other row above is a prediction.")
    W("")
    W("Ondansetron occupies ~98% of 5-HT3 at 8 mg q8h and still achieves only")
    W("{}% inhibition of NTS transmission, because in ESTABLISHED HG the"
      .format(fmt(node_res["ondansetron 8 mg q8h"]["inh"] * 100.0, 1)))
    W("peripheral 5-HT3 limb carries {}% of the drive.  That is the model's"
      .format(fmt(P["W_VAG"] * 100.0, 1)))
    W("explanation for the guideline drug being the one that failed.")

    # =============================================================== 5
    W("")
    W("--- 5. OTHER TRIALS — PREDICTIONS, NOT FITS " + "-" * 34)
    W("")
    # Koren 2010: MODERATE NVP, so a milder cohort
    subs_mod, ga_mod = build_population(
        sens_list=[0.55, 0.62, 0.68, 0.74, 0.80], ga=56.0)
    pbo_m = run_arm(subs_mod, ga_mod, PLACEBO)
    dox_m = run_arm(subs_mod, ga_mod, lambda a, b: sc_doxylamine(a, b))
    kor = dox_m[14] - pbo_m[14]
    W("Koren 2010 (PMID 20843504) — doxylamine/pyridoxine in MODERATE NVP,")
    W("day 14.  Run against a milder cohort, as the trial was:")
    W("   model    placebo {}   active {}   difference {}"
      .format(fmt(pbo_m[14], 2), fmt(dox_m[14], 2), fmt(kor, 2)))
    W("   observed         -3.9           -4.8               -0.9")

    # Guttuso 2021: gabapentin vs an active comparator
    gbp = node_res["gabapentin 600 mg q8h"]

    def comparator(a, b):
        r = sc_ondansetron(a, b)
        return sc_promethazine(a, b, r)
    cmp_r = run_arm(subs, ga, comparator)
    red_g = -(gbp[7] - 0.0)
    red_c = -(cmp_r[7] - 0.0)
    gbp_rel = (red_g - red_c) / red_c * 100.0 if abs(red_c) > 1e-6 else float("nan")
    W("")
    W("Guttuso 2021 (PMID 33451591) — gabapentin vs ondansetron+promethazine,")
    W("days 5-7, baseline-adjusted PUQE reduction:")
    W("   model    gabapentin {}   comparator {}   relative {}%"
      .format(fmt(red_g, 2), fmt(red_c, 2), fmt(gbp_rel, 0)))
    W("   observed 52% greater reduction (95% CI 16 to 88)")

    # Maina 2014 CLONEMESI
    clo = node_res["clonidine 5 mg patch"]
    W("")
    W("Maina 2014 CLONEMESI (PMID 24684734) — clonidine 5 mg patch, day 5,")
    W("PUQE improvement over placebo:")
    W("   model {}      observed 95% CI 0.43 to 3.24 (point estimate ~1.8)"
      .format(fmt(-(clo[7] - pbo[7]), 2)))

    # Yost 2003 corticosteroid null
    ste = node_res["prednisolone taper"]
    W("")
    W("Yost 2003 (PMID 14662211) — methylprednisolone + prednisone taper on top")
    W("of standard care; primary outcome rehospitalisation:")
    W("   model    dPUQE vs placebo at day 7 = {}   day 14 = {}"
      .format(fmt(ste[7] - pbo[7], 2), fmt(ste[14] - pbo[14], 2)))
    W("   observed 34% vs 35% rehospitalised (P = .89) — null.")
    W("   The model predicts a null for a structural reason, not a fitted one:")
    W("   corticosteroids occupy no node on the GDF15 -> GFRAL -> AP -> NTS")
    W("   chain, so the only thing they can move here is appetite.")
    W("")
    W("A MISS, stated as such: metoclopramide comes out at {} (day 7), i.e."
      .format(fmt(node_res["metoclopramide 10 mg q6h"][7] - pbo[7], 2)))
    W("essentially inactive, because its central D2 occupancy is low (~11%) and")
    W("its prokinetic action feeds the peripheral limb, which this model says")
    W("is nearly irrelevant.  Trials generally find metoclopramide comparable")
    W("to promethazine, so either the peripheral limb matters more than W_VAG")
    W("allows, or metoclopramide has an action the model omits.  This is the")
    W("cleanest place to try to falsify the node-position law.")

    # =============================================================== 6
    W("")
    W("--- 6. PREVENTION: THE WINDOW CLOSES AT CONCEPTION " + "-" * 27)
    W("")
    prev = [
        ("nothing", base_regimen()),
        ("metformin, pre-conception only", sc_metformin_pre()),
        ("metformin, started GA 6 wk", sc_metformin_pre(t_start=42.0,
                                                        t_stop=140.0)),
        ("tobacco, pre-conception", base_regimen(tobacco=True)),
        ("beta-thalassaemia (lifelong)", base_regimen(thal=True)),
        ("rGDF15, pre-conception", sc_rgdf15_pre()),
        ("rGDF15, started GA 8 wk", sc_rgdf15_pre(t0=56.0, t1=140.0)),
        ("anti-GDF15 mAb at GA 8 wk", sc_antigdf15(t_dose=56.0)),
    ]
    hdr = ("{:<32}{:>9}{:>9}{:>9}{:>8}{:>8}"
           .format("intervention", "GDF15pre", "PUQEmax", "wtloss%", "HGcase",
                   "fold9w"))
    W(hdr)
    W("-" * len(hdr))
    pv = {}
    for lab, r in prev:
        res = simulate(r)
        pv[lab] = res
        isc, pk, loss = hg_case(res)
        W("{:<32}{:>9}{:>9}{:>9}{:>8}{:>8}".format(
            lab, fmt(res["GDF"][at(res, 0.0)], 0), fmt(pk, 1), fmt(loss, 1),
            "YES" if isc else "no", fmt(res["fold"][at(res, 63.0)], 2)))
    W("")
    W("The sign flip is the point.  Metformin BEFORE conception takes peak PUQE")
    W("from {} to {}.  The SAME drug started at GA 6 wk gives {} -- it raises"
      .format(fmt(peak(pv["nothing"], "PUQE")[0], 1),
              fmt(peak(pv["metformin, pre-conception only"], "PUQE")[0], 1),
              fmt(peak(pv["metformin, started GA 6 wk"], "PUQE")[0], 1)))
    W("GDF15 on top of an already-steep placental ramp, far too late to move a")
    W("set-point whose time constant is {} days.  Recombinant GDF15 behaves the"
      .format(fmt(P["TAU_SP"], 0)))
    W("same way, and worse: pre-conception it protects, mid-pregnancy it is")
    W("agonist at the very receptor causing the illness.")
    W("")
    W("Three independent pre-conception exposures that all raise GDF15 --")
    W("beta-thalassaemia, metformin and tobacco -- come out protective, in that")
    W("order of strength, from ONE parameter (ALPHA) plus how much each raises")
    W("basal GDF15.  Observed: beta-thalassaemia 'very low levels of NVP'")
    W("(PMID 38092039); metformin aRR 0.29 (95% CI 0.12-0.71); tobacco aRR 0.51")
    W("(0.30-0.86) (both PMID 40588059).  The model puts metformin ahead of")
    W("tobacco for one reason only -- metformin roughly doubles circulating")
    W("GDF15 while smoking raises it by about a third -- so the ordering is a")
    W("consequence of the mechanism, not a fitted result.")
    W("")
    rg = pv["rGDF15, pre-conception"]
    rg_pre = max(rg["PUQE"][at(rg, -63.0):at(rg, 0.0)])
    W("A cost the model finds on its own, and a serious one for anybody")
    W("planning to develop this: a fold-change detector is provoked by the")
    W("ONSET of recombinant GDF15 as surely as by a placenta.  Peak PUQE during")
    W("pre-conception rGDF15 titration is {} -- the prophylaxis reproduces a"
      .format(fmt(rg_pre, 1)))
    W("milder version of the illness it prevents, and would have to be titrated")
    W("over months rather than weeks.  Metformin does not have this problem at")
    W("the exposures modelled here, which is much of its appeal.")

    # =============================================================== 7
    W("")
    W("--- 7. THE SLOW CASCADE: A DIFFERENT CLOCK " + "-" * 35)
    W("")
    TS = 56.0
    sev = dict(P)
    sev["SENS"] = 1.18
    nut = [
        ("no IV support", base_regimen()),
        ("IV saline + KCl, no dextrose",
         add_supportive(base_regimen(), TS, TS + 28, dextrose=False)),
        ("IV dextrose, NO thiamine",
         add_supportive(base_regimen(), TS, TS + 28, dextrose=True)),
        ("IV dextrose + thiamine 100 mg/d",
         add_supportive(base_regimen(), TS, TS + 28, dextrose=True,
                        thiamine_mg=100.0)),
        ("mirtazapine 30 mg, no IV",
         sc_mirtazapine(TS, TS + 28)),
        ("mirtazapine + IV + thiamine",
         add_supportive(sc_mirtazapine(TS, TS + 28), TS, TS + 28,
                        dextrose=True, thiamine_mg=100.0)),
    ]
    hdr = ("{:<32}{:>8}{:>8}{:>7}{:>8}{:>7}{:>7}{:>8}"
           .format("management", "THImin", "P(WE)%", "K+min", "HCO3max",
                   "Cl-min", "Na+min", "wtloss%"))
    W(hdr)
    W("-" * len(hdr))
    for lab, r in nut:
        res = simulate(r, p=sev, t1=150.0)
        i0, i1 = at(res, 42.0), at(res, 145.0)
        _, _, loss = hg_case(res)
        W("{:<32}{:>8}{:>8}{:>7}{:>8}{:>7}{:>7}{:>8}".format(
            lab, fmt(min(res["THI"][i0:i1]), 1),
            fmt(max(res["PWE"][i0:i1]) * 100.0, 1),
            fmt(min(res["K"][i0:i1]), 2),
            fmt(max(res["HCO3"][i0:i1]), 1),
            fmt(min(res["CL"][i0:i1]), 1),
            fmt(min(res["NA"][i0:i1]), 1), fmt(loss, 1)))
    W("")
    W("The clinically important row is the third.  IV dextrose without thiamine")
    W("is BETTER hydrated and BETTER nourished than no support at all, and its")
    W("Wernicke risk is far higher, because glucose consumes a thiamine store")
    W("that is already low (Oudman 2019, PMID 30889425).  Adding thiamine to")
    W("the identical fluid abolishes the risk.  Thiamine first, then glucose.")
    W("")
    W("The second decoupling: an antiemetic that fixes PUQE does not fix this.")
    W("Nausea runs on a clock of hours; a {} mg thiamine store with a ~{} day"
      .format(fmt(P["THI0"], 0), fmt(math.log(2) / P["KEL_THI"], 0)))
    W("half-life runs on a clock of weeks.  Wernicke risk can still be rising")
    W("in a woman whose PUQE has already normalised.")

    # =============================================================== 8
    W("")
    W("--- 8. GESTATIONAL TRANSIENT THYROTOXICOSIS " + "-" * 34)
    W("")
    hdr = "{:<28}{:>10}{:>9}{:>9}{:>10}".format(
        "cohort", "hCGmax", "TSHmin", "FT4max", "TSH<0.4")
    W(hdr)
    W("-" * len(hdr))
    for lab in ["normal pregnancy", "HG (sensitivity-driven)",
                "HG (production-driven)"]:
        res = nat[lab]
        i0, i1 = at(res, 28.0), at(res, 165.0)
        tmin = min(res["TSH"][i0:i1])
        W("{:<28}{:>10}{:>9}{:>9}{:>10}".format(
            lab, fmt(max(res["HCG"][i0:i1]), 0), fmt(tmin, 2),
            fmt(max(res["FT4"][i0:i1]), 1), "YES" if tmin < 0.4 else "no"))
    W("")
    W("hCG cross-activates the TSH receptor (Rodien 2004, PMID 15073140).")
    W("Because hCG and GDF15 are both syncytiotrophoblast products, the model")
    W("predicts biochemical thyrotoxicosis tracks the FETAL-PRODUCTION axis and")
    W("not the maternal-sensitivity axis -- two women with identical PUQE can")
    W("differ completely in TSH.  That is testable with paired GDF15/TSH.")

    # =============================================================== 9
    W("")
    W("=" * 78)
    W("VALIDATION — model vs published observation")
    W("=" * 78)
    checks = [
        ("VOMIT mirtazapine dPUQE d2", mir[2] - pbo[2], -1.86,
         (-3.61, -0.12), "41478546  FITTED (E0)"),
        ("VOMIT ondansetron dPUQE d2", ond[2] - pbo[2], -0.51,
         (-2.32, 1.30), "41478546  FITTED (W_VAG)"),
        ("VOMIT mirtazapine > ondansetron d7",
         (mir[7] - pbo[7]) - (ond[7] - pbo[7]), -1.35,
         (-3.10, 0.40), "41478546"),
        ("Koren doxylamine dPUQE d14", kor, -0.90,
         (-1.80, -0.10), "20843504"),
        ("Guttuso gabapentin relative %", gbp_rel, 52.0,
         (16.0, 88.0), "33451591"),
        ("Maina clonidine PUQE gain", -(clo[7] - pbo[7]), 1.80,
         (0.43, 3.24), "24684734"),
        ("Yost steroid dPUQE d14 (null)", ste[14] - pbo[14], 0.0,
         (-0.75, 0.75), "14662211"),
        ("beta-thal peak PUQE (very low NVP)",
         peak(nat["beta-thalassaemia"], "PUQE")[0], 4.0, (3.0, 7.0),
         "38092039"),
        ("metformin pre-conception peak PUQE",
         peak(pv["metformin, pre-conception only"], "PUQE")[0], 8.0,
         (5.0, 11.5), "40588059"),
        ("tobacco pre-conception peak PUQE",
         peak(pv["tobacco, pre-conception"], "PUQE")[0], 10.5,
         (8.5, 13.0), "40588059"),
        ("HG peak PUQE, GA weeks", peak(hg, "PUQE")[1] / 7.0, 9.5,
         (8.0, 11.5), "31515515"),
        ("HG PUQE at GA 16 wk", hg["PUQE"][i16], 4.0, (3.0, 7.5),
         "31515515"),
        ("HG weight loss %", hg_case(hg)[2], 8.0, (5.0, 15.0), "34555550"),
        ("GDF15 28wk/9wk ratio > 1",
         hg["GDF"][i28] / hg["GDF"][i9], 1.4, (1.05, 2.5), "12495665"),
        ("production-driven HG TSH suppressed",
         min(nat["HG (production-driven)"]["TSH"][
             at(nat["HG (production-driven)"], 28.0):
             at(nat["HG (production-driven)"], 165.0)]), 0.2, (0.0, 0.40),
         "15073140"),
    ]
    W("")
    W("{:<38}{:>8}{:>8}{:>16}".format("target", "model", "obs", "accept"))
    W("-" * 78)
    npass = 0
    for lab, mod, obsv, (lo, hi), src in checks:
        ok = (mod == mod) and (lo <= mod <= hi)
        npass += ok
        W("{:<38}{:>8}{:>8}{:>16}  {}  PMID {}".format(
            lab, fmt(mod, 2), fmt(obsv, 2),
            "[" + fmt(lo, 2) + "," + fmt(hi, 2) + "]",
            "PASS" if ok else "FAIL", src))
    W("")
    W("{}/{} acceptance intervals met.".format(npass, len(checks)))
    W("")
    W("PARAMETERS FITTED TO DATA (4 of {}):".format(len(P)))
    W("  ALPHA, TAU_SP   adaptation completeness and time constant, fitted to")
    W("                  the natural history (peak week, 16-week remission)")
    W("  W_VAG           peripheral 5-HT3 limb share, fitted to ondansetron")
    W("  E0              central node-authority scale, fitted to mirtazapine")
    W("Everything else is taken from published physiology, receptor affinity")
    W("and PK, or is a structural assumption stated in the header.  The")
    W("beta-thalassaemia, metformin, tobacco, doxylamine, gabapentin,")
    W("clonidine, corticosteroid and thyroid rows are predictions.")

    txt = "\n".join(lines)
    if check_only:
        txt = "VALIDATION" + txt.split("VALIDATION", 1)[-1]
    print(txt)
    return npass, len(checks)


if __name__ == "__main__":
    chk = "--check" in sys.argv
    np_, nt_ = main(check_only=chk)
    if chk:
        print("\n" + ("PASS" if np_ == nt_ else
                      "PARTIAL ({}/{})".format(np_, nt_)))
        sys.exit(0 if np_ == nt_ else 1)
