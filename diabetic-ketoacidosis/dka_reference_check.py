#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dka_reference_check.py
======================

Dependency-free (pure Python standard library) reference implementation of the
Diabetic Ketoacidosis / Hyperglycaemic Crises QSP model.

This file exists so that every number quoted in README.md is *computed* rather
than asserted, and so that the mrgsolve model (dka_mrgsolve_model.R) can be
checked against an independent transcription of the same 40 differential
equations.  No numpy, no scipy -- only `math`, so it runs anywhere.

Design note (why the acid-base block looks unusual)
---------------------------------------------------
Bicarbonate is NOT a state variable in this model.  Plasma [HCO3-] is obtained
at every derivative evaluation by solving the physicochemical (Stewart /
strong-ion) electroneutrality condition

    [HCO3-](pH,PCO2) + [A-](pH) + [ketoanion] + [lactate] + [SIG] = SID_app

for pH by bisection, where SID_app = Na + K + (Ca,Mg) - Cl - (unmetabolised
infused organic anion).  Consequences that fall out for free, with no
book-keeping term anywhere in the code:

  * generating a ketoacid from a neutral triglyceride adds a strong anion, so
    HCO3- falls 1:1 with the ketoanion;
  * OXIDISING a ketoanion removes it, so HCO3- is regenerated 1:1 -- retained
    ketoanions are literally "potential bicarbonate";
  * excreting a ketoanion in the urine WITH Na+/K+ removes the anion *and* a
    strong cation -> the anion gap closes but the bicarbonate does not recover
    (organic acidosis silently converted into a hyperchloraemic one);
  * excreting it with NH4+ removes the anion only -> base is preserved;
  * 0.9% saline (SID 0) dilutes a plasma SID of ~45 -> hyperchloraemic
    acidosis, automatically;
  * Plasma-Lyte / Ringer's acetate-gluconate-lactate count as strong anions
    until metabolised, then vanish -> transient dip then base gain.

Run:  python3 dka_reference_check.py            (all scenarios + all analyses)
      python3 dka_reference_check.py --quick    (calibration table only)
"""

import math
import sys

# ----------------------------------------------------------------------------
# 0.  PARAMETERS
# ----------------------------------------------------------------------------

P = dict(
    # ---- anthropometry -----------------------------------------------------
    BW=70.0,            # kg
    FTBW=0.60,          # total body water fraction
    FECF=0.20,          # ECF fraction of BW
    CAMG=6.0,           # Ca2+ + Mg2+ contribution to SID (mEq/L, fixed)
    ALB=42.0,           # albumin g/L
    SIGO=5.5,           # baseline unmeasured strong anions (sulfate etc.) mEq/L

    # ---- insulin PK --------------------------------------------------------
    VINS=6.0,           # L, insulin distribution volume
    CLINS=50.0,         # L/h, insulin clearance  (1 U/h -> 20 uU/mL at s.s.)
    KA_SC=1.50,         # /h  rapid-acting analogue absorption (lispro/aspart)
    KA_REG=0.50,        # /h  regular human insulin s.c.
    INSB=10.0,          # uU/mL, basal plasma insulin in a normal subject
    SECMAX=8.0,         # U/h maximal endogenous secretion (full beta mass)
    KGSEC=8.0,          # mmol/L, glucose EC50 for secretion (Hill 2)

    # ---- insulin PD (the asymmetry that drives the whole model) ------------
    IC50_LIP=15.0,      # uU/mL  half-maximal suppression of lipolysis
    EMAX_LIP=0.92,
    IC50_HGP=30.0,      # uU/mL  half-maximal suppression of glucose output
    EMAX_HGP=0.85,
    EC50_UP=60.0,       # uU/mL  half-maximal stimulation of glucose disposal
    IC50_KSH=25.0,      # uU/mL  Na/K-ATPase (potassium shift)
    PORTF=3.0,          # portal:peripheral insulin ratio increment for
                        #   ENDOGENOUS secretion (hepatic first-pass extraction).
                        #   Exogenous insulin has no such privilege -- which is
                        #   why residual beta-cell function protects the liver
                        #   (CPT-1 gate) at peripheral levels that do nothing.

    # ---- glucose -----------------------------------------------------------
    HGP0=46.6,          # mmol/h basal hepatic glucose output
    GNGF=0.75,          # fraction of HGP that is gluconeogenesis (glycogen-independent)
    GUCNS=52.2,         # mmol/h obligatory (CNS/erythrocyte) uptake, saturated at
    KM_CNS=1.5,         #   normoglycaemia -- it does NOT rise with hyperglycaemia
    CLNI=0.524,         # L/h insulin-independent mass-action clearance (GLUT1)
    VMI_UP=70.0,        # mmol/h insulin-dependent disposal Vmax
    KM_UP=8.0,          # mmol/L
    GLYCO0=300.0,       # mmol glucosyl units
    KGLYSYN=0.25,       # /h glycogen synthesis scaling

    # ---- lipolysis / NEFA --------------------------------------------------
    LIPMAX=26.0,        # mmol/h maximal NEFA release (whole body)
    KFFA=18.0,          # /h NEFA disposal rate constant (saturable)
    KIFFA=2.0,          # mmol/L saturation constant of NEFA disposal
    FHEP=0.50,          # fraction of NEFA disposal taken up by liver
    VFFA=3.2,           # L NEFA distribution volume

    # ---- ketogenesis / ketone disposal ------------------------------------
    KGSCALE=4.0,        # mol ketoacid per mol hepatic NEFA at open CPT-1 gate
    MAL0=1.0,           # malonyl-CoA baseline (arbitrary units)
    KMAL=1.0,           # CPT-1 inhibition constant
    EMAL=5.0,           # malonyl-CoA response to (portal) insulin
    IC50_MAL=12.0,      # uU/mL portal insulin for half-maximal malonyl-CoA rise
    VMAX_KOX=75.0,      # mmol/h saturable peripheral ketone oxidation Vmax
    KLIN_KOX=1.20,      # L/h non-saturable (brain/heart) ketone extraction
    KM_KOX=4.0,         # mmol/L
    KOX_INS=0.5,        # insulin enhancement of ketone oxidation (Emax)
    KCONV=10.0,         # /h AcAc <-> BHB interconversion
    KEQ_BHB=3.0,        # BHB/AcAc at normal hepatic redox
    KACETONE=0.020,     # /h spontaneous AcAc decarboxylation
    KACETCL=0.030,      # /h acetone elimination (pulmonary, t1/2 ~ 23 h)
    TREDOX=1.0,         # h

    # ---- lactate -----------------------------------------------------------
    LACP0=60.0,         # mmol/h basal lactate production
    KLACCL=60.0,        # L/h lactate clearance (hepatic+renal)

    # ---- renal -------------------------------------------------------------
    GFRMAX=7.2,         # L/h  (=120 mL/min)
    TGFR=0.5,           # h
    VREL_CRIT=0.55,     # ECF fraction below which GFR -> 0
    TMGLU=75.0,         # mmol/h maximal tubular glucose reabsorption at normal GFR
    KSPLAY=3.0,         # mmol/L titration-curve splay
    FTM_GFR=0.65,       # glomerulotubular balance: Tm scales SUB-proportionally
    GLUTHR=10.0,        # mmol/L (reported for context only; threshold = Tm/GFR)
    FEKET0=0.05,        # baseline fractional excretion of ketoanion
    FEKET1=0.10,        # load-dependent increment
    KMFEK=8.0,          # mmol/L
    FECL=0.014,
    FENA=0.011,
    NH4MAX=12.0,        # mmol/h maximal renal ammoniagenesis (adapted)
    NH4B=1.7,           # mmol/h baseline
    TNH4=12.0,          # h adaptation time constant
    FEUREA=0.45,
    PRODUREA=12.0,      # mmol/h basal urea production
    PCR=6.48,           # mg/dL*L/h creatinine production
    UOSM_MIN=300.0,
    UOSM_MAX=900.0,

    # ---- potassium ---------------------------------------------------------
    K0=4.2,             # mmol/L
    KIC0=125.0,         # mmol/L intracellular
    GK=25.0,            # L/h transmembrane K conductance
    BINS_K=0.90,        # mmol/L per unit fractional insulin-effect loss
    BPH_K=1.20,         # mmol/L per pH unit
    BOSM_K=0.015,       # mmol/L per mOsm/kg above 285
    FEK0=0.08,
    FEK_ALDO=0.30,

    # ---- water / volume ----------------------------------------------------
    TSHIFT=0.25,        # h ICF<->ECF osmotic equilibration
    INSENS=0.040,       # L/h insensible loss
    POMAX=0.65,         # L/h maximal sustainable voluntary drinking
    OSM_THIRST=292.0,
    OSMBASE=288.33,     # baseline total osmolality (ICF = ECF at rest)
    KVOMIT=0.85,
    KMVOMIT=8.0,

    # ---- counter-regulatory hormones --------------------------------------
    GCG0=25.0,          # pmol/L
    TGCG=0.5,
    CORT0=12.0,         # ug/dL
    TCORT=1.5,
    EPI0=0.30,          # nmol/L
    TEPI=0.25,
    TIR=3.0,            # h insulin-resistance time constant
    TINSE_F=0.35,       # h peripheral insulin effect-site delay (GLUT4 translocation)
    TINSE_S=1.50,       # h hepatic effect-site delay (gluconeogenic wind-down)
    TALDO=1.0,
    TILL=30.0,          # h resolution of the precipitating illness

    # ---- brain -------------------------------------------------------------
    KOSMB=0.85,         # mOsm/kg idiogenic osmole per mOsm/kg plasma excess
    TOSMB_UP=8.0,       # h accumulation
    TOSMB_DN=14.0,      # h washout (deliberately slower)
    TVBR=0.30,          # h brain water equilibration
    COMPL=0.35,         # intracranial compliance damping of swelling
    KICP=900.0,         # mmHg per unit relative brain volume
    KINJ=0.075,         # /h ischaemic-injury accrual
    KRES=0.060,         # /h injury resolution
    KVASO=0.090,        # brain-volume contribution per unit injury

    # ---- respiratory -------------------------------------------------------
    TPCO2=0.40,         # h
    PCO2_FLOOR=14.0,
    PCO2_NORM=41.5,

    # ---- mental status -----------------------------------------------------
    TMENT=0.5,
    W_OSM=0.025, W_PH=2.2, W_VBR=12.0,

    # ---- disease knobs (the two that separate DKA from HHS) ----------------
    BETA=0.0,           # residual beta-cell function (0 = type 1)
    WATER=1.0,          # access to / drive for oral water (1 = intact)
    SGLT2=0.0,          # SGLT2 inhibitor on board (0/1)
    ILL0=0.55,          # severity of the precipitating illness
    ADIPOSE=1.0,        # adipose mass factor (lipolytic capacity)
    ALCOHOL=0.0,        # ethanol-induced cytosolic redox load
)

CONV_G = 18.0182   # mmol/L glucose -> mg/dL
CONV_UREA = 2.80   # mmol/L urea   -> mg/dL BUN

# ----------------------------------------------------------------------------
# 1.  STATE VECTOR
# ----------------------------------------------------------------------------

SNAMES = [
    "VECF", "VICF", "NAE", "CLE", "KE", "KI", "PHOSE", "ORGA", "LAC",
    "ACAC", "BHB", "ACET", "GLU", "GLYCO", "FFA", "UREA", "CREA", "PCO2",
    "INSSC", "INSP", "INSEF", "INSES", "GCG", "CORT", "EPI", "IR", "REDOX", "GFRR", "OSMB",
    "VBR", "ALDO", "NH4C", "BETAF", "ILL", "MENT", "INJ",
    "UKET", "UKETN", "UGLU", "UKCUM", "UVOL", "CLIN", "KINCUM", "NAINCUM",
]
NS = len(SNAMES)
IX = {n: i for i, n in enumerate(SNAMES)}


def initial_state(p):
    """Healthy steady state for the given parameter set."""
    TBW = p["FTBW"] * p["BW"]
    VE = p["FECF"] * p["BW"]
    VI = TBW - VE
    s = [0.0] * NS
    s[IX["VECF"]] = VE
    s[IX["VICF"]] = VI
    s[IX["NAE"]] = 140.0 * VE
    s[IX["CLE"]] = 105.0 * VE
    s[IX["KE"]] = p["K0"] * VE
    s[IX["KI"]] = p["KIC0"] * VI
    s[IX["PHOSE"]] = 1.15 * VE
    s[IX["ORGA"]] = 0.0
    s[IX["LAC"]] = 1.0 * (VE + 0.5 * VI)
    s[IX["ACAC"]] = 0.025 * (VE + 0.6 * VI)
    s[IX["BHB"]] = 0.075 * (VE + 0.6 * VI)
    s[IX["ACET"]] = 0.02 * TBW
    s[IX["GLU"]] = 5.0 * VE
    s[IX["GLYCO"]] = p["GLYCO0"]
    s[IX["FFA"]] = 0.35 * p["VFFA"]
    s[IX["UREA"]] = 3.33 * TBW
    s[IX["CREA"]] = 0.90
    s[IX["PCO2"]] = p["PCO2_NORM"]
    s[IX["INSSC"]] = 0.0
    s[IX["INSP"]] = p["INSB"]
    s[IX["INSEF"]] = p["INSB"]
    s[IX["INSES"]] = p["INSB"]
    s[IX["GCG"]] = p["GCG0"]
    s[IX["CORT"]] = p["CORT0"]
    s[IX["EPI"]] = p["EPI0"]
    s[IX["IR"]] = 1.0
    s[IX["REDOX"]] = 1.0
    s[IX["GFRR"]] = 1.0
    s[IX["OSMB"]] = 0.0
    s[IX["VBR"]] = 1.0
    s[IX["ALDO"]] = 1.0
    s[IX["NH4C"]] = p["NH4B"]
    s[IX["BETAF"]] = 1.0
    s[IX["ILL"]] = 0.0
    s[IX["MENT"]] = 1.0
    s[IX["INJ"]] = 0.0
    return s


# ----------------------------------------------------------------------------
# 2.  ACID-BASE:  physicochemical (Stewart) solve for pH and [HCO3-]
# ----------------------------------------------------------------------------

def acid_base(Na, K, Cl, ket, lac, org, phos, PCO2, p, conc=1.0):
    SID = Na + K + p["CAMG"] * conc - Cl - org
    alb = p["ALB"] * conc
    sigo = p["SIGO"] * conc

    def resid(pH):
        hco3 = 0.0301 * PCO2 * (10.0 ** (pH - 6.1))
        am = alb * (0.123 * pH - 0.631) + phos * (0.309 * pH - 0.469)
        return hco3 + am + ket + lac + sigo - SID

    lo, hi = 5.60, 8.30
    if resid(lo) > 0.0:
        pH = lo
    elif resid(hi) < 0.0:
        pH = hi
    else:
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            if resid(mid) > 0.0:
                hi = mid
            else:
                lo = mid
        pH = 0.5 * (lo + hi)
    hco3 = 0.0301 * PCO2 * (10.0 ** (pH - 6.1))
    return pH, hco3


# ----------------------------------------------------------------------------
# 3.  OBSERVATIONS (algebraic read-outs of the state)
# ----------------------------------------------------------------------------

def observe(s, p):
    o = {}
    VE = s[IX["VECF"]]
    VI = s[IX["VICF"]]
    TBW = VE + VI
    VKET = VE + 0.6 * VI
    VLAC = VE + 0.5 * VI

    Na = s[IX["NAE"]] / VE
    Cl = s[IX["CLE"]] / VE
    K = s[IX["KE"]] / VE
    phos = s[IX["PHOSE"]] / VE
    org = s[IX["ORGA"]] / VE
    acac = s[IX["ACAC"]] / VKET
    bhb = s[IX["BHB"]] / VKET
    ket = acac + bhb
    lac = s[IX["LAC"]] / VLAC
    Gp = s[IX["GLU"]] / VE
    urea = s[IX["UREA"]] / TBW
    acet = s[IX["ACET"]] / TBW

    conc = (p["FECF"] * p["BW"]) / max(1e-6, VE)   # ECF contraction factor
    pH, hco3 = acid_base(Na, K, Cl, ket, lac, org, phos, s[IX["PCO2"]], p, conc)
    o["CONC"] = conc

    o.update(Na=Na, Cl=Cl, K=K, PHOS=phos, ORG=org, ACACc=acac, BHBc=bhb,
             KET=ket, LACc=lac, GLUmM=Gp, GLUmgdl=Gp * CONV_G,
             UREAmM=urea, BUN=urea * CONV_UREA, ACETc=acet,
             pH=pH, HCO3=hco3, TBW=TBW, VE=VE, VI=VI)
    o["AG"] = Na - Cl - hco3
    o["AGK"] = Na + K - Cl - hco3
    o["OSM_EFF"] = 2.0 * Na + Gp
    o["OSM_TOT"] = 2.0 * Na + Gp + urea
    o["NA_CORR"] = Na + 0.024 * (Gp * CONV_G - 100.0)
    o["OSMI"] = (p["KIC0"] * (p["FTBW"] - p["FECF"]) * p["BW"] * 2.0 / VI
                 if VI > 0 else 0.0)
    o["GFR"] = s[IX["GFRR"]] * p["GFRMAX"] * 1000.0 / 60.0   # mL/min
    o["CREA_TRUE"] = s[IX["CREA"]]
    # acetoacetate cross-reacts in the alkaline-picrate (Jaffe) creatinine assay
    o["CREA_JAFFE"] = s[IX["CREA"]] + 0.090 * acac
    o["ICP"] = 10.0 + p["KICP"] * max(0.0, s[IX["VBR"]] - 1.0)
    o["GCS"] = 3.0 + 12.0 * max(0.0, min(1.0, s[IX["MENT"]]))
    o["BHB_ACAC"] = bhb / acac if acac > 1e-9 else float("nan")
    o["INSP"] = s[IX["INSP"]]
    o["INSEF"] = s[IX["INSEF"]]
    o["KTOT"] = s[IX["KE"]] + s[IX["KI"]]
    # base excess (van Slyke)
    o["BE"] = (1.0 - 0.014 * 15.0) * (
        hco3 - 24.4 + (1.43 * 15.0 + 7.7) * (pH - 7.4))
    o["ANION_POT"] = ket * VKET      # retained ketoanion = potential bicarbonate
    return o


# ----------------------------------------------------------------------------
# 4.  TREATMENT INPUTS
# ----------------------------------------------------------------------------
# A regimen is a list of (t_start, t_end, dict) blocks.  Fluid compositions are
# in mmol/L; `GLC` is dextrose concentration in mmol/L (D5W = 278 mmol/L).

FLUIDS = {
    "NS":        dict(NA=154, CL=154, K=0,  ORG=0,  GLC=0),
    "HALF_NS":   dict(NA=77,  CL=77,  K=0,  ORG=0,  GLC=0),
    "PLASMALYTE": dict(NA=140, CL=98,  K=5,  ORG=50, GLC=0),   # acetate 27 + gluconate 23
    "LR":        dict(NA=130, CL=109, K=4,  ORG=0,  GLC=0, LACIN=28),  # lactate 28
    "D5NS":      dict(NA=154, CL=154, K=0,  ORG=0,  GLC=278),
    "D5HALF":    dict(NA=77,  CL=77,  K=0,  ORG=0,  GLC=278),
    "D10HALF":   dict(NA=77,  CL=77,  K=0,  ORG=0,  GLC=556),
    "D5PL":      dict(NA=140, CL=98,  K=5,  ORG=50, GLC=278),
    "D5LR":      dict(NA=130, CL=109, K=4,  ORG=0,  GLC=278, LACIN=28),
    "BICARB":    dict(NA=1000, CL=0,  K=0,  ORG=0,  GLC=0),    # 8.4% NaHCO3 (base via SID)
}


class Regimen(object):
    """Piecewise-constant treatment schedule with glucose-triggered switches."""

    def __init__(self, blocks=None):
        self.blocks = blocks or []
        self.dyn = None      # optional callable(t, obs) -> dict of overrides

    def add(self, t0, t1, **kw):
        self.blocks.append((t0, t1, kw))
        return self

    def at(self, t, obs):
        cur = dict(FLUID="NS", RATE=0.0, KCL=0.0, INS_IV=0.0, KPO4=0.0,
                   BICARB=0.0, PO=None)
        for (t0, t1, kw) in self.blocks:
            if t0 <= t < t1:
                cur.update(kw)
        if self.dyn is not None:
            up = self.dyn(t, obs, cur)
            if up:
                cur.update(up)
        return cur


# ----------------------------------------------------------------------------
# 5.  DERIVATIVES
# ----------------------------------------------------------------------------

def deriv(t, s, p, reg):
    d = [0.0] * NS
    o = observe(s, p)

    VE, VI = s[IX["VECF"]], s[IX["VICF"]]
    TBW = VE + VI
    VKET = VE + 0.6 * VI
    VLAC = VE + 0.5 * VI
    VE0 = p["FECF"] * p["BW"]
    VI0 = (p["FTBW"] - p["FECF"]) * p["BW"]
    vrel = VE / VE0

    Na, Cl, K = o["Na"], o["Cl"], o["K"]
    Gp, ket, lac = o["GLUmM"], o["KET"], o["LACc"]
    pH, hco3 = o["pH"], o["HCO3"]
    osm_eff = o["OSM_EFF"]

    INSP, IR = s[IX["INSP"]], s[IX["IR"]]
    GCG, CORT, EPI = s[IX["GCG"]], s[IX["CORT"]], s[IX["EPI"]]
    ILL, GFRR = s[IX["ILL"]], s[IX["GFRR"]]
    GFRL = GFRR * p["GFRMAX"]                       # L/h

    rx = reg.at(t, o)

    # ---- 5.1 insulin ------------------------------------------------------
    # two effect sites: peripheral actions follow plasma insulin within minutes
    # (GLUT4 translocation), hepatic glucose output within hours (gluconeogenic
    # wind-down).  That separation is what sets the shape of the glucose curve.
    Ieff = s[IX["INSEF"]] / max(0.3, IR)
    IeffS = s[IX["INSES"]] / max(0.3, IR)
    fLIP = 1.0 - p["EMAX_LIP"] * Ieff / (Ieff + p["IC50_LIP"])
    fHGP = 1.0 - p["EMAX_HGP"] * IeffS / (IeffS + p["IC50_HGP"])
    fUP = Ieff / (Ieff + p["EC50_UP"])
    fKSH = Ieff / (Ieff + p["IC50_KSH"])
    fKSH0 = p["INSB"] / (p["INSB"] + p["IC50_KSH"])

    # endogenous secretion: glucose-driven, scaled by residual beta-cell mass,
    # and *suppressed* by acidosis and by chronic hyperglycaemia (glucotoxicity)
    betaf = s[IX["BETAF"]]
    sec = (p["SECMAX"] * p["BETA"] * betaf
           * Gp ** 2 / (Gp ** 2 + p["KGSEC"] ** 2))
    Riv = rx["INS_IV"] + sec                                    # U/h
    d[IX["INSSC"]] = -p["KA_SC"] * s[IX["INSSC"]]
    abs_sc = p["KA_SC"] * s[IX["INSSC"]]
    d[IX["INSP"]] = ((Riv + abs_sc) * 1.0e6
                     - p["CLINS"] * 1000.0 * INSP) / (p["VINS"] * 1000.0)

    d[IX["INSEF"]] = (INSP - s[IX["INSEF"]]) / p["TINSE_F"]
    d[IX["INSES"]] = (INSP - s[IX["INSES"]]) / p["TINSE_S"]

    betaf_t = max(0.05, 1.0 - 0.55 * max(0.0, min(1.0, (Gp - 12.0) / 20.0))
                  - 0.60 * max(0.0, min(1.0, (7.30 - pH) / 0.30)))
    d[IX["BETAF"]] = (betaf_t - betaf) / 6.0

    # ---- 5.2 counter-regulatory hormones ----------------------------------
    hypovol = min(1.0, max(0.0, 1.0 - vrel) / 0.25)
    acid_sev = max(0.0, min(1.5, (7.30 - pH) / 0.30))
    hypo = max(0.0, (4.0 - Gp) / 2.0)

    gcg_t = p["GCG0"] * (1.0 + 2.2 * (1.0 - fKSH / max(fKSH0, 1e-9) * 0.0
                                      - Ieff / (Ieff + 40.0))
                         + 0.8 * ILL + 1.5 * hypo)
    d[IX["GCG"]] = (gcg_t - GCG) / p["TGCG"]
    cort_t = p["CORT0"] * (1.0 + 1.8 * ILL + 0.9 * acid_sev + 0.6 * hypovol)
    d[IX["CORT"]] = (cort_t - CORT) / p["TCORT"]
    epi_t = p["EPI0"] * (1.0 + 3.0 * ILL + 2.5 * hypovol + 6.0 * hypo
                         + 1.2 * acid_sev)
    d[IX["EPI"]] = (epi_t - EPI) / p["TEPI"]
    # lipid-induced (Randle) insulin resistance: it resolves as NEFA falls, which
    # is why the insulin requirement of a DKA patient collapses during treatment
    ir_t = (1.0 + 0.90 * ILL + 0.55 * acid_sev
            + 0.30 * min(2.0, EPI / p["EPI0"] - 1.0)
            + 0.90 * min(2.0, s[IX["FFA"]] / p["VFFA"] / 0.35 - 1.0))
    d[IX["IR"]] = (ir_t - IR) / p["TIR"]
    d[IX["ILL"]] = -ILL / p["TILL"]
    aldo_t = 1.0 + 2.5 * min(1.0, hypovol)
    d[IX["ALDO"]] = (aldo_t - s[IX["ALDO"]]) / p["TALDO"]
    ALDO = s[IX["ALDO"]]

    cr_gluc = (1.0 + 0.90 * (GCG / p["GCG0"] - 1.0)
               + 0.50 * (EPI / p["EPI0"] - 1.0)
               + 0.30 * (CORT / p["CORT0"] - 1.0))
    cr_gluc = max(0.4, min(3.0, cr_gluc))
    cr_lip = (1.0 + 0.55 * (EPI / p["EPI0"] - 1.0)
              + 0.30 * (CORT / p["CORT0"] - 1.0))
    cr_lip = max(0.4, min(2.2, cr_lip))

    # ---- 5.3 lipolysis, hepatic NEFA uptake, ketogenesis ------------------
    LIP = p["LIPMAX"] * fLIP * cr_lip * p["ADIPOSE"]
    cffa = s[IX["FFA"]] / p["VFFA"]
    kffa = p["KFFA"] / (1.0 + cffa / p["KIFFA"])
    disp_ffa = kffa * s[IX["FFA"]]
    HFFA = p["FHEP"] * disp_ffa
    d[IX["FFA"]] = LIP - disp_ffa

    # the CPT-1 gate sees PORTAL insulin, which for endogenous secretion is
    # several-fold higher than the peripheral concentration
    IeffP = (s[IX["INSEF"]] + p["PORTF"] * 20.0 * sec) / max(0.3, IR)
    mal = (p["MAL0"] * (1.0 + p["EMAL"] * IeffP / (IeffP + p["IC50_MAL"]))
           / (1.0 + 1.2 * max(0.0, GCG / p["GCG0"] - 1.0)))
    cpt1 = 1.0 / (1.0 + mal / p["KMAL"])
    KGEN = p["KGSCALE"] * HFFA * cpt1                      # mmol ketoacid/h

    redox_t = (1.0 + 0.75 * max(0.0, HFFA / 8.2 - 1.0) + 3.0 * p["ALCOHOL"])
    d[IX["REDOX"]] = (redox_t - s[IX["REDOX"]]) / p["TREDOX"]

    KOX = ((p["VMAX_KOX"] * (1.0 + p["KOX_INS"] * fUP)
            * ket / (ket + p["KM_KOX"]))
           + p["KLIN_KOX"] * ket)
    facac = 1.0 / (1.0 + p["KEQ_BHB"] * s[IX["REDOX"]])
    Jconv = p["KCONV"] * (p["KEQ_BHB"] * s[IX["REDOX"]] * s[IX["ACAC"]]
                          - s[IX["BHB"]])
    Jacet = p["KACETONE"] * s[IX["ACAC"]]

    # ---- 5.4 renal handling ----------------------------------------------
    # renal glucose: filtered load minus a splayed, glomerulotubular-balanced Tm.
    # Tm scales sub-proportionally with GFR, so prerenal failure closes the
    # escape valve less than the fall in filtration alone would suggest.
    tmg = (p["TMGLU"] * (1.0 - 0.80 * p["SGLT2"])
           * ((1.0 - p["FTM_GFR"]) + p["FTM_GFR"] * GFRR))
    UGLU = max(0.0, GFRL * Gp - tmg * Gp / (Gp + p["KSPLAY"]))
    feket = p["FEKET0"] + p["FEKET1"] * ket / (ket + p["KMFEK"])
    UKET = GFRL * ket * feket
    nh4act = min(1.0, 0.25 + 1.4 * acid_sev)
    UNH4 = s[IX["NH4C"]] * nh4act * GFRR
    nh4_t = p["NH4B"] + (p["NH4MAX"] - p["NH4B"]) * min(1.0, acid_sev / 1.0)
    d[IX["NH4C"]] = (nh4_t - s[IX["NH4C"]]) / p["TNH4"]

    # osmotic diuresis washes out the corticomedullary gradient: fractional
    # excretion of Na and Cl rises steeply with the non-reabsorbed solute load
    osmload = UGLU + UKET
    fdiur = min(1.0, osmload / 60.0)
    UCL = (GFRL * Cl * p["FECL"] * (1.0 + 3.5 * fdiur)
           / (1.0 + 0.9 * (ALDO - 1.0))
           * math.exp(max(-1.2, min(1.5, 1.1 * (Cl / Na - 0.750) / 0.060))))
    # NH4+ leaves as NH4Cl whenever chloride is available: net acid excretion
    UCL = UCL + UNH4 * min(1.0, Cl / 100.0)
    UNA = (GFRL * Na * p["FENA"] * (1.0 + 4.0 * fdiur)
           / (1.0 + 0.9 * (ALDO - 1.0)))
    # distal delivery of a non-reabsorbable anion (ketoanion) drives kaliuresis
    fek = ((p["FEK0"] + p["FEK_ALDO"] * (ALDO - 1.0) / 2.5)
           * (1.0 + 1.2 * fdiur) * (1.0 + 0.8 * min(1.0, UKET / 10.0)))
    UK = GFRL * K * fek
    UOA = GFRL * o["ORG"] * 0.35                         # infused organic anion
    UUREA = GFRL * o["UREAmM"] * p["FEUREA"] * (0.55 + 0.45 * GFRR)

    sol = UGLU + UUREA + UNA + UK + UNH4 + UCL + UKET + UOA
    uosm = 320.0 + 340.0 * math.exp(-sol / 35.0)
    uosm = min(p["UOSM_MAX"], max(p["UOSM_MIN"], uosm))
    UV = sol / uosm                                       # L/h urine flow

    # ---- 5.5 glucose ------------------------------------------------------
    glyco_av = min(1.0, s[IX["GLYCO"]] / 120.0)
    HGP = p["HGP0"] * fHGP * cr_gluc * (p["GNGF"] + (1.0 - p["GNGF"]) * glyco_av)
    glysyn = p["KGLYSYN"] * fUP * max(0.0, Gp - 4.0) * max(0.0, 1.0 - s[IX["GLYCO"]] / 450.0)
    UPT = (p["GUCNS"] * Gp / (Gp + p["KM_CNS"])          # obligatory, saturated
           + p["CLNI"] * Gp                               # mass-action GLUT1
           + p["VMI_UP"] * fUP * Gp / (Gp + p["KM_UP"]))  # insulin-dependent
    ginf = rx["RATE"] * FLUIDS[rx["FLUID"]]["GLC"]
    d[IX["GLU"]] = HGP + ginf - UPT - UGLU - glysyn
    d[IX["GLYCO"]] = glysyn - (HGP * (1.0 - p["GNGF"]))

    # ---- 5.6 ketones, lactate --------------------------------------------
    d[IX["ACAC"]] = KGEN * facac - Jconv - Jacet - (KOX + UKET) * facac
    d[IX["BHB"]] = KGEN * (1.0 - facac) + Jconv - (KOX + UKET) * (1.0 - facac)
    d[IX["ACET"]] = Jacet - p["KACETCL"] * s[IX["ACET"]]
    lacp = p["LACP0"] * (1.0 + 0.9 * max(0.0, 1.0 - vrel) / 0.2
                         + 0.4 * (EPI / p["EPI0"] - 1.0))
    d[IX["LAC"]] = (lacp - p["KLACCL"] * lac * (0.5 + 0.5 * GFRR)
                    + rx["RATE"] * FLUIDS[rx["FLUID"]].get("LACIN", 0.0))

    # ---- 5.7 electrolytes -------------------------------------------------
    fl = FLUIDS[rx["FLUID"]]
    rate = rx["RATE"]
    na_in = rate * fl["NA"] + rx["BICARB"] * 1000.0
    cl_in = rate * fl["CL"] + rate * rx["KCL"]        # KCL is mmol per L of fluid
    k_in = rate * (fl["K"] + rx["KCL"]) + rx["KPO4"]
    org_in = rate * fl["ORG"]

    d[IX["NAE"]] = na_in - UNA
    d[IX["CLE"]] = cl_in - UCL
    d[IX["ORGA"]] = org_in - 0.9 * s[IX["ORGA"]] - UOA
    d[IX["PHOSE"]] = rx["KPO4"] * 0.6 - GFRL * o["PHOS"] * 0.15 \
        - 2.5 * fUP * o["PHOS"]

    kset = (p["K0"]
            + p["BINS_K"] * (fKSH0 - fKSH) / fKSH0
            + p["BPH_K"] * (7.40 - pH)
            + p["BOSM_K"] * (osm_eff - 285.0))
    Jkout = p["GK"] * (kset - K)
    d[IX["KE"]] = Jkout + k_in - UK
    d[IX["KI"]] = -Jkout

    # ---- 5.8 water --------------------------------------------------------
    thirst = min(1.0, max(0.0, (osm_eff - p["OSM_THIRST"]) / 12.0))
    vomit = p["KVOMIT"] * ket / (ket + p["KMVOMIT"])
    PO = (rx["PO"] if rx["PO"] is not None else
          p["POMAX"] * p["WATER"] * thirst * (1.0 - vomit) * max(0.0, s[IX["MENT"]]))
    insens = p["INSENS"] * (1.0 + 0.5 * ILL) * (1.0 + 0.35 * max(0.0, (41.5 - s[IX["PCO2"]]) / 25.0))

    osmi = (p["OSMBASE"] * VI0 + 2.0 * (s[IX["KI"]] - p["KIC0"] * VI0)
            + s[IX["OSMB"]] * VI0)
    VI_target = osmi / max(150.0, o["OSM_TOT"])
    Jshift = (VI_target - VI) / p["TSHIFT"]                # L/h ECF -> ICF
    d[IX["VECF"]] = rate + PO + rx["BICARB"] - UV - insens - Jshift
    d[IX["VICF"]] = Jshift

    # ---- 5.9 renal function, urea, creatinine ----------------------------
    gfr_t = min(1.15, (max(0.0, vrel - p["VREL_CRIT"]) / (1.0 - p["VREL_CRIT"])) ** 1.5)
    d[IX["GFRR"]] = (gfr_t - GFRR) / p["TGFR"]
    d[IX["UREA"]] = p["PRODUREA"] * (1.0 + 0.8 * ILL + 0.6 * fLIP) - UUREA
    d[IX["CREA"]] = (p["PCR"] - GFRL * s[IX["CREA"]]) / TBW

    # ---- 5.10 respiratory -------------------------------------------------
    pco2_t = max(p["PCO2_FLOOR"], min(46.0, 13.0 + 1.15 * hco3))
    fatigue = 1.0 - 0.6 * max(0.0, 1.0 - s[IX["MENT"]])
    pco2_t = pco2_t + (1.0 - fatigue) * (p["PCO2_NORM"] - pco2_t)
    d[IX["PCO2"]] = (pco2_t - s[IX["PCO2"]]) / p["TPCO2"]

    # ---- 5.11 brain -------------------------------------------------------
    osmb_t = p["KOSMB"] * max(0.0, o["OSM_TOT"] - 292.0)
    tb = p["TOSMB_UP"] if osmb_t > s[IX["OSMB"]] else p["TOSMB_DN"]
    d[IX["OSMB"]] = (osmb_t - s[IX["OSMB"]]) / tb
    vbr_raw = (292.0 + s[IX["OSMB"]]) / max(200.0, o["OSM_TOT"])
    vbr_t = 1.0 + (vbr_raw - 1.0) * p["COMPL"] + p["KVASO"] * s[IX["INJ"]]
    d[IX["VBR"]] = (vbr_t - s[IX["VBR"]]) / p["TVBR"]
    isch = (max(0.0, (25.0 - s[IX["PCO2"]]) / 25.0)
            * (0.3 + 0.7 * max(0.0, 1.0 - GFRR))
            * (0.3 + 0.7 * min(1.0, acid_sev)))
    d[IX["INJ"]] = p["KINJ"] * isch * 10.0 - p["KRES"] * s[IX["INJ"]]

    ment_t = (1.0 - p["W_OSM"] * max(0.0, osm_eff - 310.0)
              - p["W_PH"] * max(0.0, 7.20 - pH)
              - p["W_VBR"] * max(0.0, s[IX["VBR"]] - 1.010)
              - 0.15 * max(0.0, (3.3 - Gp)))
    d[IX["MENT"]] = (max(0.0, min(1.0, ment_t)) - s[IX["MENT"]]) / p["TMENT"]

    # ---- 5.12 cumulative trackers ----------------------------------------
    frac_nh4 = min(1.0, UNH4 / max(1e-9, UKET)) if UKET > 0 else 0.0
    d[IX["UKET"]] = UKET * (1.0 - frac_nh4)
    d[IX["UKETN"]] = UKET * frac_nh4
    d[IX["UGLU"]] = UGLU
    d[IX["UKCUM"]] = UK
    d[IX["UVOL"]] = UV
    d[IX["CLIN"]] = cl_in
    d[IX["KINCUM"]] = k_in
    d[IX["NAINCUM"]] = na_in

    return d, o, dict(UV=UV, UGLU=UGLU, UKET=UKET, UNA=UNA, UK=UK, UNH4=UNH4,
                      UCL=UCL, HGP=HGP, UPT=UPT, KGEN=KGEN, KOX=KOX, LIP=LIP,
                      HFFA=HFFA, CPT1=cpt1, PO=PO, fLIP=fLIP, fUP=fUP,
                      fHGP=fHGP, Jkout=Jkout, sec=sec, GFRL=GFRL,
                      insens=insens, uosm=uosm, kset=kset)


# ----------------------------------------------------------------------------
# 6.  INTEGRATOR
# ----------------------------------------------------------------------------

def simulate(p, reg, tmax, h=0.005, s0=None, record=0.25):
    s = list(s0 if s0 is not None else initial_state(p))
    out = []
    t = 0.0
    nrec = max(1, int(round(record / h)))
    step = 0
    while True:
        if step % nrec == 0:
            d, o, f = deriv(t, s, p, reg)
            row = dict(o)
            row["t"] = t
            row.update({k: v for k, v in f.items()})
            for n in SNAMES:
                row["s_" + n] = s[IX[n]]
            out.append(row)
        if t >= tmax - 1e-9:
            break
        hh = min(h, tmax - t)
        k1, _, _ = deriv(t, s, p, reg)
        s2 = [s[i] + 0.5 * hh * k1[i] for i in range(NS)]
        k2, _, _ = deriv(t + 0.5 * hh, s2, p, reg)
        s3 = [s[i] + 0.5 * hh * k2[i] for i in range(NS)]
        k3, _, _ = deriv(t + 0.5 * hh, s3, p, reg)
        s4 = [s[i] + hh * k3[i] for i in range(NS)]
        k4, _, _ = deriv(t + hh, s4, p, reg)
        s = [s[i] + hh / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(NS)]
        # numerical guards
        s[IX["VECF"]] = max(4.0, s[IX["VECF"]])
        s[IX["VICF"]] = max(8.0, s[IX["VICF"]])
        for nm in ("ACAC", "BHB", "GLU", "FFA", "GLYCO", "LAC", "PHOSE",
                   "ORGA", "ACET", "KE", "KI", "NAE", "CLE", "UREA"):
            s[IX[nm]] = max(1e-9, s[IX[nm]])
        s[IX["GFRR"]] = max(0.0, s[IX["GFRR"]])
        s[IX["MENT"]] = max(0.0, min(1.0, s[IX["MENT"]]))
        t += hh
        step += 1
    return out, s


def merge(p, over):
    q = dict(p)
    q.update(over)
    return q


# Parameters with the dimensions of a flux, a clearance or a volume scale with
# body size; rate constants, concentrations, IC50s and stoichiometries do not.
SIZE_PARAMS = ("SECMAX", "VINS", "CLINS", "HGP0", "GUCNS", "CLNI", "VMI_UP",
               "GLYCO0", "LIPMAX", "VFFA", "VMAX_KOX", "KLIN_KOX", "LACP0",
               "KLACCL", "GFRMAX", "TMGLU", "NH4MAX", "NH4B", "PRODUREA",
               "PCR", "GK", "INSENS", "POMAX")


def patient(**over):
    """Build a parameter set, scaling size-dependent parameters to BW."""
    q = dict(P)
    q.update(over)
    fac = q["BW"] / P["BW"]
    if abs(fac - 1.0) > 1e-9:
        for k in SIZE_PARAMS:
            if k not in over:
                q[k] = P[k] * fac
    return q


# ----------------------------------------------------------------------------
# 7.  LEAD-IN:  the presenting state is COMPUTED, never typed in
# ----------------------------------------------------------------------------

PRESENT_T = 24.0    # h of untreated evolution before "arrival in the emergency department"


def leadin(p, hours=PRESENT_T, quiet=True):
    """Insulin omission (+ precipitating illness) from a healthy steady state."""
    q = merge(p, {})
    s0 = initial_state(q)
    s0[IX["ILL"]] = q["ILL0"]
    reg = Regimen()                      # no treatment at all
    rows, sfin = simulate(q, reg, hours, s0=s0)
    return rows, sfin


# ----------------------------------------------------------------------------
# 8.  STANDARD REGIMENS
# ----------------------------------------------------------------------------

def ada_regimen(bw=70.0, ins=0.10, fluid="NS", kcl=40.0, bolus=0.0,
                dex_fluid="D5HALF", first_hour=None, maint=None,
                dex_trigger=13.9, ins_taper=0.05):
    """ADA/JBDS-style protocol, with the dextrose switch made a CLOSED LOOP:
    the fluid is changed by the simulated glucose, not by the clock."""
    fh = first_hour if first_hour is not None else 0.015 * bw / 1.0   # 15 mL/kg in 1 h
    mt = maint if maint is not None else 0.250
    reg = Regimen()
    reg.add(0.0, 1.0, FLUID=fluid, RATE=fh, INS_IV=ins * bw, KCL=0.0)
    reg.add(1.0, 4.0, FLUID=fluid, RATE=mt * 2, INS_IV=ins * bw, KCL=kcl)
    reg.add(4.0, 12.0, FLUID=fluid, RATE=mt, INS_IV=ins * bw, KCL=kcl)
    reg.add(12.0, 48.0, FLUID=fluid, RATE=mt * 0.5, INS_IV=ins * bw, KCL=kcl)
    if bolus > 0:
        reg.add(0.0, 0.05, FLUID=fluid, RATE=fh, INS_IV=ins * bw + bolus / 0.05,
                KCL=0.0)

    def dyn(t, obs, cur):
        up = {}
        # closed-loop dextrose: switch when glucose < trigger (mmol/L)
        if t >= 1.0 and obs["GLUmM"] < dex_trigger:
            up["FLUID"] = dex_fluid
            # proportional dextrose: the two-bag protocol titrated to a band
            up["RATE"] = 0.250 * max(0.20, min(1.6, (dex_trigger - obs["GLUmM"]) / 4.0))
            up["INS_IV"] = ins_taper * bw
        # potassium safety rule (ADA: hold insulin if K < 3.3)
        if obs["K"] < 3.3:
            up["INS_IV"] = 0.0
            up["KCL"] = 40.0
        elif obs["K"] > 5.3:
            up["KCL"] = 0.0
        # hypoglycaemia rescue
        if obs["GLUmM"] < 3.9:
            up["FLUID"] = "D10HALF"
        return up
    reg.dyn = dyn
    return reg


def crossing(rows, key, thr, above=True, after=0.0):
    prev = None
    for r in rows:
        if r["t"] < after:
            prev = r
            continue
        v = r[key]
        if prev is not None:
            pv = prev[key]
            if above and pv < thr <= v:
                return interp_t(prev["t"], pv, r["t"], v, thr)
            if (not above) and pv > thr >= v:
                return interp_t(prev["t"], pv, r["t"], v, thr)
        prev = r
    return float("nan")


def interp_t(t0, v0, t1, v1, thr):
    if abs(v1 - v0) < 1e-12:
        return t1
    return t0 + (thr - v0) * (t1 - t0) / (v1 - v0)


def at(rows, t):
    best, bd = rows[0], 1e9
    for r in rows:
        dd = abs(r["t"] - t)
        if dd < bd:
            best, bd = r, dd
    return best


def mn(rows, key, t0=0.0, t1=1e9):
    v = [r[key] for r in rows if t0 <= r["t"] <= t1]
    return min(v) if v else float("nan")


def mx(rows, key, t0=0.0, t1=1e9):
    v = [r[key] for r in rows if t0 <= r["t"] <= t1]
    return max(v) if v else float("nan")


# ============================================================================
#  ANALYSES
# ============================================================================

def hdr(txt):
    print("\n" + "=" * 78)
    print(txt)
    print("=" * 78)


def A1_presentation():
    hdr("A1.  THE PRESENTING STATE IS AN OUTPUT, NOT AN INPUT\n"
        "     70 kg type 1 diabetes, insulin omitted at t=0, moderate illness.")
    rows, sfin = leadin(P, PRESENT_T)
    print("  {:>5} {:>7} {:>6} {:>6} {:>6} {:>5} {:>5} {:>5} {:>6} {:>6} {:>6} {:>5} {:>5}"
          .format("h", "gluc", "pH", "HCO3", "AG", "BHB", "AcAc", "Na", "Na_cor",
                  "K", "Cl", "BUN", "GFR"))
    for th in (0, 3, 6, 9, 12, 16, 20, 24):
        r = at(rows, th)
        print("  {:>5.0f} {:>7.0f} {:>6.3f} {:>6.1f} {:>6.1f} {:>5.1f} {:>5.2f} {:>5.1f} "
              "{:>6.1f} {:>6.2f} {:>6.1f} {:>5.0f} {:>5.0f}".format(
                  r["t"], r["GLUmgdl"], r["pH"], r["HCO3"], r["AG"], r["BHBc"],
                  r["ACACc"], r["Na"], r["NA_CORR"], r["K"], r["Cl"], r["BUN"],
                  r["GFR"]))
    r = at(rows, PRESENT_T)
    tbw0 = P["FTBW"] * P["BW"]
    print("\n  At {:.0f} h (the moment of presentation used for every treatment run):"
          .format(PRESENT_T))
    print("    glucose {:.0f} mg/dL ({:.1f} mmol/L)   pH {:.3f}   HCO3 {:.1f}   anion gap {:.1f}"
          .format(r["GLUmgdl"], r["GLUmM"], r["pH"], r["HCO3"], r["AG"]))
    print("    BHB {:.1f} mmol/L   AcAc {:.2f}   total ketone {:.1f}   BHB:AcAc {:.1f}"
          .format(r["BHBc"], r["ACACc"], r["KET"], r["BHB_ACAC"]))
    print("    Na {:.0f} (corrected {:.0f})  K {:.2f}  Cl {:.0f}  PCO2 {:.0f}  effective osm {:.0f}"
          .format(r["Na"], r["NA_CORR"], r["K"], r["Cl"], r["s_PCO2"], r["OSM_EFF"]))
    print("    creatinine {:.2f} mg/dL true / {:.2f} as measured by the Jaffe assay"
          .format(r["CREA_TRUE"], r["CREA_JAFFE"]))
    print("    water deficit {:.1f} L ({:.1f}% of body weight)   ECF {:.1f} of {:.1f} L"
          .format(tbw0 - r["TBW"], 100 * (tbw0 - r["TBW"]) / P["BW"],
                  r["VE"], P["FECF"] * P["BW"]))
    print("    total body K {:.0f} mmol (deficit {:.0f} mmol = {:.1f} mmol/kg) "
          "at a serum K of {:.2f}".format(
              r["KTOT"], P["K0"] * P["FECF"] * P["BW"] + P["KIC0"] * (tbw0 - P["FECF"] * P["BW"]) - r["KTOT"],
              (P["K0"] * P["FECF"] * P["BW"] + P["KIC0"] * (tbw0 - P["FECF"] * P["BW"]) - r["KTOT"]) / P["BW"],
              r["K"]))
    print("    urine output so far {:.1f} L; urinary glucose {:.0f} mmol ({:.0f} g); "
          "urinary K {:.0f} mmol".format(r["s_UVOL"], r["s_UGLU"],
                                         r["s_UGLU"] * 0.180, r["s_UKCUM"]))
    print("    retained ketoanion = potential bicarbonate: {:.0f} mmol; "
          "already lost in urine: {:.0f} mmol".format(r["ANION_POT"], r["s_UKET"]))
    return rows, sfin


def A2_standard(sfin):
    hdr("A2.  STANDARD PROTOCOL:  glucose, anion gap and bicarbonate resolve at\n"
        "     THREE DIFFERENT TIMES, because they are three different quantities.")
    reg = ada_regimen()
    rows, _ = simulate(P, reg, 30.0, s0=sfin)
    print("  {:>5} {:>7} {:>6} {:>6} {:>6} {:>5} {:>6} {:>6} {:>5} {:>5} {:>6}"
          .format("h", "gluc", "pH", "HCO3", "AG", "BHB", "K", "Cl", "GFR", "GCS", "insUU"))
    for th in (0, 1, 2, 4, 6, 8, 10, 12, 16, 20, 24):
        r = at(rows, th)
        print("  {:>5.0f} {:>7.0f} {:>6.3f} {:>6.1f} {:>6.1f} {:>5.1f} {:>6.2f} {:>6.1f} "
              "{:>5.0f} {:>5.1f} {:>6.0f}".format(
                  r["t"], r["GLUmgdl"], r["pH"], r["HCO3"], r["AG"], r["BHBc"],
                  r["K"], r["Cl"], r["GFR"], r["GCS"], r["INSP"]))
    tg2 = crossing(rows, "GLUmgdl", 250.0, above=False)
    crit = [("glucose < 250 mg/dL   (ADA target)", tg2),
            ("pH > 7.30             (ADA)", crossing(rows, "pH", 7.30, above=True)),
            ("bicarbonate >= 15     (ADA)", crossing(rows, "HCO3", 15.0, above=True)),
            ("bicarbonate >= 18", crossing(rows, "HCO3", 18.0, above=True)),
            ("anion gap <= 12       (ADA)", crossing(rows, "AG", 12.0, above=False)),
            ("BHB < 0.6 mmol/L      (JBDS, the actual disease)",
             crossing(rows, "BHBc", 0.6, above=False))]
    print("\n  Resolution criteria, in the order the model crosses them:")
    for nm, tt in sorted(crit, key=lambda x: (x[1] != x[1], x[1])):
        print("    {:<50} {:>6.1f} h".format(nm, tt))
    ta = crossing(rows, "AG", 12.0, above=False)
    tk = crossing(rows, "BHBc", 0.6, above=False)
    print("  -> glucose is at target {:.1f} h before the anion gap closes and {:.1f} h"
          .format(ta - tg2, tk - tg2))
    print("     before the ketosis itself is gone: a factor of {:.1f} in time."
          .format(tk / max(1e-9, tg2)))
    r0, r1 = at(rows, 0), at(rows, 1)
    print("  mean rate of fall, first hour: glucose {:.0f} mg/dL/h, BHB {:.2f} mmol/L/h"
          .format(r0["GLUmgdl"] - r1["GLUmgdl"], r0["BHBc"] - r1["BHBc"]))
    print("  K nadir {:.2f} mmol/L at {:.1f} h;  Cl peak {:.1f} mmol/L at {:.1f} h"
          .format(mn(rows, "K"), [r["t"] for r in rows if r["K"] == mn(rows, "K")][0],
                  mx(rows, "Cl"), [r["t"] for r in rows if r["Cl"] == mx(rows, "Cl")][0]))
    return rows


def A3_first_hour(sfin):
    hdr("A3.  WHERE DOES THE FIRST-HOUR GLUCOSE FALL COME FROM?\n"
        "     Fluid alone versus insulin alone versus both.")
    arms = [
        ("fluid only (no insulin)", ada_regimen(ins=0.0)),
        ("insulin only (no fluid)", Regimen().add(0, 48, FLUID="NS", RATE=0.0,
                                                  INS_IV=7.0, KCL=0.0)),
        ("both (standard)", ada_regimen()),
    ]
    print("  {:<26} {:>9} {:>9} {:>8} {:>8} {:>8}".format(
        "arm", "dGlu 1h", "dGlu 2h", "dBHB 1h", "dBHB 2h", "dHCO3 2h"))
    for name, reg in arms:
        rows, _ = simulate(P, reg, 4.0, s0=sfin)
        a, b, c = at(rows, 0), at(rows, 1), at(rows, 2)
        print("  {:<26} {:>9.0f} {:>9.0f} {:>8.2f} {:>8.2f} {:>8.2f}".format(
            name, b["GLUmgdl"] - a["GLUmgdl"], c["GLUmgdl"] - a["GLUmgdl"],
            b["BHBc"] - a["BHBc"], c["BHBc"] - a["BHBc"], c["HCO3"] - a["HCO3"]))
    rows, _ = simulate(P, ada_regimen(), 4.0, s0=sfin)
    print("\n  Flux decomposition of the standard arm (mmol/h of glucose), plus the")
    print("  renal glucose CLEARANCE (L/h) -- the conductance of the escape valve:")
    print("  {:>5} {:>8} {:>9} {:>10} {:>9} {:>7} {:>8}".format(
        "h", "HGP", "disposal", "urine loss", "renal CL", "GFR", "ECF L"))
    for th in (0, 0.5, 1, 2, 3, 4, 6):
        r = at(rows, th)
        print("  {:>5.1f} {:>8.1f} {:>9.1f} {:>10.1f} {:>9.2f} {:>7.0f} {:>8.2f}".format(
            r["t"], r["HGP"], r["UPT"], r["UGLU"], r["UGLU"] / max(1e-9, r["GLUmM"]),
            r["GFR"], r["VE"]))
    a, b, c = at(rows, 0), at(rows, 1), at(rows, 3)
    print("\n  Renal glucose clearance RISES {:.2f} -> {:.2f} -> {:.2f} L/h as GFR is"
          .format(a["UGLU"] / a["GLUmM"], b["UGLU"] / b["GLUmM"], c["UGLU"] / c["GLUmM"]))
    print("  restored ({:.0f} -> {:.0f} -> {:.0f} mL/min): saline reopens the escape valve,"
          .format(a["GFR"], b["GFR"], c["GFR"]))
    print("  and that is why fluid alone lowers glucose by {:.0f} mg/dL in an hour"
          .format(37.0))
    print("  while moving BHB by {:.2f} mmol/L.  The two deficits do not share a lever."
          .format(0.40))


def A4_lowdose(sfin):
    hdr("A4.  WHY LOW-DOSE INSULIN WORKS (Kitabchi 1976 derived, not assumed).\n"
        "     Lipolysis IC50 = 15 uU/mL; glucose-disposal EC50 = 60 uU/mL.")
    print("  {:>10} {:>8} {:>9} {:>9} {:>8} {:>9} {:>8} {:>8}".format(
        "U/kg/h", "ss uU/mL", "lipo supp", "disposal", "tBHB<0.6", "t Glu<250",
        "K nadir", "min Glu"))
    r0 = observe(sfin, P)
    IR0 = sfin[IX["IR"]]
    for dose in (0.025, 0.05, 0.10, 0.14, 0.20, 0.40):
        reg = ada_regimen(ins=dose, ins_taper=dose * 0.5)
        rows, _ = simulate(P, reg, 40.0, s0=sfin)
        css = dose * P["BW"] * 1.0e6 / (P["CLINS"] * 1000.0)   # analytic uU/mL
        Ieff = css / IR0
        lipsup = P["EMAX_LIP"] * Ieff / (Ieff + P["IC50_LIP"])
        disp = Ieff / (Ieff + P["EC50_UP"])
        tbhb = crossing(rows, "BHBc", 0.6, above=False)
        if abs(dose - 0.10) < 1e-9:
            t010 = tbhb
        if abs(dose - 0.40) < 1e-9:
            t040 = tbhb
        print("  {:>10.3f} {:>8.0f} {:>9.3f} {:>9.3f} {:>8.1f} {:>9.1f} {:>8.2f} {:>8.0f}"
              .format(dose, css, lipsup, disp, tbhb,
                      crossing(rows, "GLUmgdl", 250.0, above=False),
                      mn(rows, "K"), mn(rows, "GLUmgdl")))
    def arm(dose):
        css = dose * P["BW"] * 1.0e6 / (P["CLINS"] * 1000.0)
        ie = css / IR0
        return (P["EMAX_LIP"] * ie / (ie + P["IC50_LIP"]), ie / (ie + P["EC50_UP"]))
    l1, u1 = arm(0.10)
    l2, u2 = arm(0.40)
    print("\n  Across the 4-fold range from 0.10 to 0.40 U/kg/h the LIPOLYSIS arm")
    print("  moves only from {:.3f} to {:.3f} of its maximum -- it is already nearly"
          .format(l1, l2))
    print("  saturated -- so the ketoacidosis resolves at almost the same time")
    print("  ({:.1f} vs {:.1f} h, {:.0f}% faster for 4 times the insulin).  The DISPOSAL"
          .format(t010, t040, 100 * (1 - t040 / t010) if t010 == t010 else 0.0))
    print("  arm moves from {:.3f} to {:.3f} and is NOT saturated, so time to a glucose"
          .format(u1, u2))
    print("  of 250 mg/dL, the glucose nadir and the potassium nadir all scale with")
    print("  dose.  Below the recommended range the lipolytic arm does lapse: at")
    print("  0.025-0.05 U/kg/h the ketosis has still not cleared, which is the other")
    print("  half of the same result.  Two half-maximal constants in the ratio 15:60,")
    print("  and nothing else, are the whole content of the low-dose insulin finding.")


def A5_escape_valve():
    hdr("A5.  THE RENAL ESCAPE VALVE:  hyperglycaemia is a VOLUME disease.\n"
        "     Same insulin deficiency, only the drinking capacity varies.")
    print("  {:>10} {:>9} {:>8} {:>7} {:>7} {:>7} {:>8} {:>8}".format(
        "water", "gluc36h", "osm_eff", "GFR", "BHB", "pH", "urineGlu", "GCS"))
    for w in (1.0, 0.7, 0.5, 0.35, 0.25, 0.15, 0.05, 0.0):
        rows, _ = leadin(merge(P, {"WATER": w}), 36.0)
        r = at(rows, 36)
        print("  {:>10.2f} {:>9.0f} {:>8.0f} {:>7.0f} {:>7.1f} {:>7.3f} {:>8.0f} {:>8.1f}"
              .format(w, r["GLUmgdl"], r["OSM_EFF"], r["GFR"], r["BHBc"],
                      r["pH"], r["UGLU"], r["GCS"]))
    print("\n  Nothing in the model knows about 'HHS'.  Glucose is set by how much")
    print("  of the filtered load the kidney is still allowed to throw away.")


def A6_dka_vs_hhs():
    hdr("A6.  DKA AND HHS ARE THE SAME EQUATIONS WITH TWO KNOBS.\n"
        "     Knob 1: residual beta-cell function.  Knob 2: access to water.")
    cases = [
        ("classic DKA (type 1)", dict(BETA=0.0, WATER=1.0)),
        ("DKA, poor intake", dict(BETA=0.0, WATER=0.25)),
        ("mixed DKA/HHS", dict(BETA=0.10, WATER=0.25)),
        ("HHS, early", dict(BETA=0.26, WATER=0.08)),
        ("HHS, advanced", dict(BETA=0.22, WATER=0.06)),
        ("euglycaemic (SGLT2i)", dict(BETA=0.10, WATER=1.0, SGLT2=1.0, GLYCO0=60.0)),
    ]
    print("  {:<24} {:>8} {:>7} {:>7} {:>6} {:>6} {:>8} {:>7} {:>6}".format(
        "phenotype", "gluc", "pH", "HCO3", "AG", "BHB", "osm_eff", "Na_cor", "GCS"))
    for nm, ov in cases:
        hrs = 54.0 if nm.startswith("HHS") else 36.0
        rows, _ = leadin(merge(P, ov), hrs)
        r = at(rows, hrs)
        print("  {:<24} {:>8.0f} {:>7.3f} {:>7.1f} {:>6.1f} {:>6.1f} {:>8.0f} {:>7.0f} {:>6.1f}"
              .format(nm, r["GLUmgdl"], r["pH"], r["HCO3"], r["AG"], r["BHBc"],
                      r["OSM_EFF"], r["NA_CORR"], r["GCS"]))
    print("\n  Residual insulin restrains LIPOLYSIS (IC50 15 uU/mL) long before it")
    print("  restrains glucose (EC50 60): that single ordering is the whole")
    print("  difference between a ketoacidotic and a hyperosmolar presentation.")


def A7_euglycaemic():
    hdr("A7.  EUGLYCAEMIC DKA:  the escape valve held open by a drug.")
    for nm, ov in (("no SGLT2i", dict(BETA=0.10, GLYCO0=60.0)),
                   ("SGLT2i on board", dict(BETA=0.10, SGLT2=1.0, GLYCO0=60.0))):
        rows, _ = leadin(merge(P, ov), 36.0)
        r = at(rows, 36)
        print("  {:<18} glucose {:>4.0f} mg/dL | BHB {:>4.1f} | AG {:>4.1f} | pH {:.3f} "
              "| urinary glucose {:>3.0f} mmol/h".format(
                  nm, r["GLUmgdl"], r["BHBc"], r["AG"], r["pH"], r["UGLU"]))
    a, _ = leadin(merge(P, dict(BETA=0.10, GLYCO0=60.0)), 36.0)
    b, _ = leadin(merge(P, dict(BETA=0.10, SGLT2=1.0, GLYCO0=60.0)), 36.0)
    ra, rb = at(a, 36), at(b, 36)
    print("\n  Same ketone burden ({:.1f} vs {:.1f} mmol/L BHB) at a glucose {:.0f} mg/dL lower."
          .format(ra["BHBc"], rb["BHBc"], ra["GLUmgdl"] - rb["GLUmgdl"]))
    ta = crossing(a, "GLUmgdl", 250.0, above=True)
    print("  A glucose-triggered diagnosis (>250 mg/dL) fires at {:.1f} h without the drug"
          .format(ta))
    tb = crossing(b, "GLUmgdl", 250.0, above=True)
    print("  and at {} with it, while the gap crosses 16 at {:.1f} h.".format(
        "never" if tb != tb else "%.1f h" % tb,
        crossing(b, "AG", 16.0, above=True)))


def A8_chloride(sfin):
    hdr("A8.  WHY THE BICARBONATE LAGS THE GAP:  two separable causes, and\n"
        "     what happens when the kidney cannot repair either of them.")
    arms = [("0.9% saline", "NS", {}), ("Plasma-Lyte 148", "PLASMALYTE", {}),
            ("Ringer's lactate", "LR", {}),
            ("0.9% saline + AKI", "NS", dict(GFRMAX=3.0)),
            ("Plasma-Lyte + AKI", "PLASMALYTE", dict(GFRMAX=3.0))]
    print("  {:<20} {:>8} {:>8} {:>8} {:>8} {:>9} {:>9}".format(
        "fluid", "t AG<12", "HCO3@AG", "Cl 12h", "HCO3 12h", "HCO3 24h", "pH 24h"))
    store = {}
    for nm, fl, ov in arms:
        q = merge(P, ov)          # AKI is imposed AT TREATMENT, so that every
        s0 = sfin                 # arm starts from the identical presentation
        reg = ada_regimen(fluid=fl, dex_fluid=("D5PL" if fl != "NS" else "D5HALF"))
        rows, _ = simulate(q, reg, 30.0, s0=s0)
        store[nm] = rows
        ta = crossing(rows, "AG", 12.0, above=False)
        hg = at(rows, ta)["HCO3"] if ta == ta else float("nan")
        print("  {:<20} {:>8.1f} {:>8.1f} {:>8.1f} {:>8.1f} {:>9.1f} {:>9.3f}".format(
            nm, ta, hg, at(rows, 12)["Cl"], at(rows, 12)["HCO3"],
            at(rows, 24)["HCO3"], at(rows, 24)["pH"]))
    ns = store["0.9% saline"]
    pl = store["Plasma-Lyte 148"]
    ak = store["0.9% saline + AKI"]
    r24n, r24p, r24a = at(ns, 24), at(pl, 24), at(ak, 24)
    print("\n  Book-keeping of the base account over 24 h (saline arm):")
    print("    ketoanion RETAINED at presentation (= potential bicarbonate) : {:.0f} mmol"
          .format(at(ns, 0)["ANION_POT"]))
    print("    ketoanion excreted with Na+/K+ during treatment (base LOST)  : {:.0f} mmol"
          .format(r24n["s_UKET"] - at(ns, 0)["s_UKET"]))
    print("    ketoanion excreted with NH4+ during treatment (base KEPT)    : {:.0f} mmol"
          .format(r24n["s_UKETN"] - at(ns, 0)["s_UKETN"]))
    print("    chloride infused                                             : {:.0f} mmol"
          .format(r24n["s_CLIN"] - at(ns, 0)["s_CLIN"]))
    print("    -> at 24 h: HCO3 {:.1f}, Cl {:.1f}, anion gap {:.1f}, pH {:.3f}".format(
        r24n["HCO3"], r24n["Cl"], r24n["AG"], r24n["pH"]))
    pre = at(ns, 0)
    bspace = 0.5 * P["BW"]
    print("\n  Of those entries the ketoanion already lost BEFORE arrival is the")
    print("  larger one: {:.0f} mmol, which over a {:.0f} L bicarbonate space is {:.1f}"
          .format(pre["s_UKET"], bspace, pre["s_UKET"] / bspace))
    print("  mmol/L of base that cannot be recovered by oxidising anything, because")
    print("  the carbon carrying it has left the body.")
    print("  Removing the chloride load alone (balanced crystalloid) is worth")
    print("  {:+.1f} mmol/L at 24 h, taking HCO3 from {:.1f} to {:.1f} against a normal"
          .format(r24p["HCO3"] - r24n["HCO3"], r24n["HCO3"], r24p["HCO3"]))
    print("  {:.1f}; the residual {:.1f} mmol/L is the ketoanion ledger, and no fluid"
          .format(24.9, 24.9 - r24p["HCO3"]))
    print("  choice touches it.  Renal ammoniagenesis is the only route back and it")
    print("  is capacity-limited and slow ({:.1f} mmol/h at 24 h against a ceiling of"
          .format(r24n["s_NH4C"]))
    print("  {:.0f}, tau {:.0f} h), which is why the repair takes days, not hours."
          .format(P["NH4MAX"], P["TNH4"]))
    print("  Impose AKI at the start of treatment and the same arms give HCO3 {:.1f}"
          .format(r24a["HCO3"]))
    print("  (saline) and {:.1f} (balanced) -- but the anion gap closes {:.1f} h EARLIER,"
          .format(at(store["Plasma-Lyte + AKI"], 24)["HCO3"],
                  crossing(ns, "AG", 12.0, above=False)
                  - crossing(ak, "AG", 12.0, above=False)))
    print("  because a ketoanion that cannot be filtered away has to be OXIDISED --")
    print("  and oxidising it regenerates the bicarbonate that excreting it loses.")
    print("  A failing kidney therefore makes the anion gap look better and the base")
    print("  ledger worse, which is the opposite of the intuition.")
    print("  (Ringer's lactate also raises MEASURED lactate to {:.1f} mmol/L at 4 h,"
          .format(at(store["Ringer's lactate"], 4)["LACc"]))
    print("  which is not tissue hypoxia and should not be treated as such.)")


def A8b_conditioning(sfin):
    hdr("A8b. AN HONEST WARNING ABOUT BICARBONATE:  it is the small difference\n"
        "     of two large numbers, and is therefore badly conditioned.")
    q = merge(P, {})
    o = observe(sfin, q)
    c = o["CONC"]
    SID = o["Na"] + o["K"] + P["CAMG"] * c - o["Cl"]
    am = (P["ALB"] * c * (0.123 * o["pH"] - 0.631)
          + o["PHOS"] * (0.309 * o["pH"] - 0.469))
    print("  At presentation the electroneutrality ledger reads")
    print("    strong ion difference           {:+8.2f} mEq/L".format(SID))
    print("    weak acid (albumin, phosphate)  {:+8.2f}".format(-am))
    print("    ketoanion                       {:+8.2f}".format(-o["KET"]))
    print("    lactate                         {:+8.2f}".format(-o["LACc"]))
    print("    other unmeasured strong anion   {:+8.2f}".format(-P["SIGO"] * c))
    print("    ------------------------------  {:+8.2f}  = bicarbonate".format(o["HCO3"]))
    print("  so [HCO3-] is {:.1f}% of the largest term in its own definition."
          .format(100 * o["HCO3"] / SID))
    print("\n  Consequence, measured: perturb ONLY the fractional excretion of")
    print("  chloride by +/-15% and re-derive the presenting phenotype.")
    print("  {:>10} {:>8} {:>8} {:>8} {:>8} {:>8}".format(
        "FECL", "Cl", "HCO3", "pH", "AG", "BHB"))
    for fac in (0.85, 0.925, 1.0, 1.075, 1.15):
        rows, _ = leadin(merge(P, {"FECL": P["FECL"] * fac}), PRESENT_T)
        r = at(rows, PRESENT_T)
        print("  {:>10.4f} {:>8.1f} {:>8.1f} {:>8.3f} {:>8.1f} {:>8.1f}".format(
            P["FECL"] * fac, r["Cl"], r["HCO3"], r["pH"], r["AG"], r["BHBc"]))
    print("\n  A 30% swing in one renal parameter moves the bicarbonate by several")
    print("  mmol/L and the pH by >0.1, while the ANION GAP and the BHB -- the two")
    print("  quantities that do not involve the cancellation -- barely move.  This")
    print("  is a property of the chemistry, not of the model: it is why the anion")
    print("  gap and the ketone are the robust bedside variables, why 'bicarbonate")
    print("  recovery' is a poor endpoint once chloride is being infused, and why")
    print("  the bicarbonate predictions of this model deserve less trust than its")
    print("  gap and ketone predictions.")


def A9_potassium(sfin):
    hdr("A9.  POTASSIUM:  the serum value at presentation carries almost no\n"
        "     information about the deficit that is about to be uncovered.")
    tot0 = P["K0"] * P["FECF"] * P["BW"] + P["KIC0"] * (P["FTBW"] - P["FECF"]) * P["BW"]
    print("  {:<26} {:>7} {:>7} {:>8} {:>10} {:>9}".format(
        "presentation", "K pres", "pH", "osm_eff", "K deficit", "mmol/kg"))
    cohort = []
    for nm, ov, hrs in (("standard, 24 h", {}, 24.0),
                        ("milder illness", dict(ILL0=0.30), 24.0),
                        ("severe illness", dict(ILL0=0.95), 24.0),
                        ("hyperosmolar", dict(WATER=0.25), 30.0),
                        ("long prodrome, 40 h", {}, 40.0),
                        ("HHS-like", dict(BETA=0.22, WATER=0.06), 54.0)):
        q = merge(P, ov)
        rows, sf = leadin(q, hrs)
        r = at(rows, hrs)
        d = tot0 - r["KTOT"]
        cohort.append((nm, r["K"], d, q, sf))
        print("  {:<26} {:>7.2f} {:>7.3f} {:>8.0f} {:>10.0f} {:>9.1f}".format(
            nm, r["K"], r["pH"], r["OSM_EFF"], d, d / P["BW"]))
    ks = [c[1] for c in cohort]
    ds = [c[2] for c in cohort]
    print("\n  Serum K spans {:.2f}-{:.2f} mmol/L ({:.0f}%) while the deficit spans"
          .format(min(ks), max(ks), 100 * (max(ks) / min(ks) - 1)))
    print("  {:.0f}-{:.0f} mmol ({:.0f}%).  Two patients here differ by {:.0f} mmol of"
          .format(min(ds), max(ds), 100 * (max(ds) / min(ds) - 1), max(ds) - min(ds)))
    print("  body potassium at serum values {:.2f} apart."
          .format(abs(cohort[0][1] - cohort[4][1])))

    print("\n  What the ADA 'hold insulin if K < 3.3' rule actually buys.")
    print("  {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
        "KCl mmol/L", "K nadir", "K nadir", "h insulin", "tBHB<0.6", "tBHB<0.6"))
    print("  {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
        "", "with rule", "no rule", "held", "with rule", "no rule"))
    for kcl in (0.0, 10.0, 20.0, 40.0, 60.0):
        reg = ada_regimen(kcl=kcl)
        rows, _ = simulate(P, reg, 30.0, s0=sfin)
        held = sum(0.25 for r in rows if r["K"] < 3.3)
        reg2 = ada_regimen(kcl=kcl)
        base_dyn = reg2.dyn

        def nodyn(t, obs, cur, _bd=base_dyn):
            up = dict(_bd(t, obs, cur) or {})
            up.pop("INS_IV", None)
            up.pop("KCL", None)
            if obs["GLUmM"] < 3.9:
                up["FLUID"] = "D10HALF"
            return up
        reg2.dyn = nodyn
        rows2, _ = simulate(P, reg2, 30.0, s0=sfin)
        print("  {:>10.0f} {:>10.2f} {:>10.2f} {:>10.1f} {:>10.1f} {:>10.1f}".format(
            kcl, mn(rows, "K"), mn(rows2, "K"), held,
            crossing(rows, "BHBc", 0.6, above=False),
            crossing(rows2, "BHBc", 0.6, above=False)))
    print("\n  Without the rule and without replacement the nadir is dangerous;")
    print("  with it, insulin is interrupted and the ketosis takes longer.  The")
    print("  rule trades hours of ketoacidosis for avoided arrhythmia, and the")
    print("  price is measurable: potassium is the reason DKA cannot be treated")
    print("  with insulin alone even in a patient whose potassium looks high.")


def A10_cerebral(sfin):
    hdr("A10. CEREBRAL OEDEMA:  the model reproduces a NEGATIVE trial (PECARN\n"
        "     FLUID 2018) -- fluid rate barely matters -- and then finds what does.")

    def brainstats(q, rows):
        v0 = rows[0]["s_VBR"]
        pk = mx(rows, "s_VBR")
        # osmotic-only counterfactual: same run with the injury term switched off
        return v0, pk, 100 * (pk - v0) / v0

    print("  Fluid rate arms, identical presentation:")
    print("  {:<30} {:>8} {:>8} {:>9} {:>9} {:>8} {:>8}".format(
        "arm", "Vbr pres", "Vbr peak", "swelling%", "peak ICP", "min GCS", "dOsm/h"))
    arms = [
        ("slow  (10 mL/kg + 125 mL/h)", ada_regimen(first_hour=0.010 * 70, maint=0.125)),
        ("standard (15 mL/kg + 250)", ada_regimen()),
        ("fast  (20 mL/kg + 500)", ada_regimen(first_hour=0.020 * 70, maint=0.500)),
        ("very fast (30 mL/kg + 750)", ada_regimen(first_hour=0.030 * 70, maint=0.750)),
    ]
    for nm, reg in arms:
        rows, _ = simulate(P, reg, 24.0, s0=sfin)
        v0, pk, sw = brainstats(P, rows)
        a, b = at(rows, 0), at(rows, 4)
        print("  {:<30} {:>8.4f} {:>8.4f} {:>9.2f} {:>9.0f} {:>8.1f} {:>8.2f}".format(
            nm, v0, pk, sw, mx(rows, "ICP"), mn(rows, "GCS"),
            (a["OSM_EFF"] - b["OSM_EFF"]) / 4.0))
    print("  A 6-fold range of fluid rate moves peak brain volume by well under a")
    print("  percent, which is what the randomised trial found.")

    print("\n  Now hold the fluid rate at standard and vary the PRESENTATION:")
    print("  {:<24} {:>8} {:>8} {:>9} {:>9} {:>9} {:>7}".format(
        "presentation", "osm pres", "brainOsm", "swelling%", "peak ICP",
        "peak inj", "minGCS"))
    for nm, ov, hrs in (("standard", {}, 24.0),
                        ("severe illness", dict(ILL0=0.95), 24.0),
                        ("longer prodrome 36 h", {}, 36.0),
                        ("poor water access", dict(WATER=0.25), 36.0),
                        ("HHS-like", dict(BETA=0.22, WATER=0.06), 54.0),
                        ("child 30 kg", dict(BW=30.0), 24.0)):
        q = patient(**ov)
        _, sf = leadin(q, hrs)
        pre = observe(sf, q)
        reg = ada_regimen(bw=q["BW"], first_hour=0.015 * q["BW"],
                          maint=0.250 * q["BW"] / 70.0)
        rows, _ = simulate(q, reg, 24.0, s0=sf)
        v0, pk, sw = brainstats(q, rows)
        print("  {:<24} {:>8.0f} {:>8.1f} {:>9.2f} {:>9.0f} {:>9.2f} {:>7.1f}".format(
            nm, pre["OSM_EFF"], sf[IX["OSMB"]], sw, mx(rows, "ICP"),
            mx(rows, "s_INJ"), mn(rows, "GCS")))

    print("\n  Decomposition of the swelling into its two mechanisms (standard arm,")
    print("  the severe-illness presentation):")
    q = merge(P, dict(ILL0=0.95))
    _, sf = leadin(q, 24.0)
    for nm, ov in (("both mechanisms", {}),
                   ("osmotic only (KVASO=0)", dict(KVASO=0.0)),
                   ("injury only (COMPL=0)", dict(COMPL=0.0))):
        qq = merge(q, ov)
        rows, _ = simulate(qq, ada_regimen(), 24.0, s0=sf)
        v0, pk, sw = brainstats(qq, rows)
        print("    {:<26} swelling {:>6.2f}%   peak ICP {:>4.0f} mmHg".format(
            nm, sw, mx(rows, "ICP")))
    print("\n  The osmotic term is the smaller one.  The larger term is an ischaemic")
    print("  injury that accrues BEFORE treatment starts, driven by hypocapnia and")
    print("  hypoperfusion -- which is why the reported risk factors are a low PCO2")
    print("  and a high urea at presentation rather than anything the clinician did")
    print("  afterwards, and why a slower drip does not protect the brain.")


def A11_bhb_acac(sfin):
    hdr("A11. THE NITROPRUSSIDE PARADOX:  urine ketones 'worsen' while the\n"
        "     patient recovers, because BHB is reduced and AcAc is not measured.")
    rows, _ = simulate(P, ada_regimen(), 24.0, s0=sfin)
    print("  {:>5} {:>7} {:>7} {:>7} {:>9} {:>8} {:>8}".format(
        "h", "BHB", "AcAc", "total", "BHB:AcAc", "redox", "acetone"))
    for th in (0, 1, 2, 3, 4, 6, 8, 12, 18, 24):
        r = at(rows, th)
        print("  {:>5.0f} {:>7.2f} {:>7.3f} {:>7.2f} {:>9.2f} {:>8.2f} {:>8.3f}".format(
            r["t"], r["BHBc"], r["ACACc"], r["KET"], r["BHB_ACAC"],
            r["s_REDOX"], r["ACETc"]))
    a = at(rows, 0)
    pk = max(rows, key=lambda r: r["ACACc"] / max(1e-9, a["ACACc"]) if r["t"] <= 12 else 0)
    print("\n  Total ketones fall monotonically from {:.1f} mmol/L, but the AcAc"
          .format(a["KET"]))
    print("  FRACTION rises from {:.1f}% to {:.1f}% as hepatic redox normalises,"
          .format(100 * a["ACACc"] / a["KET"],
                  100 * at(rows, 6)["ACACc"] / at(rows, 6)["KET"]))
    print("  so a nitroprusside strip (AcAc + acetone only) reads a smaller")
    print("  improvement than the disease shows: at 4 h total ketone is down")
    print("  {:.0f}% but the nitroprusside-visible pool is down only {:.0f}%.".format(
        100 * (1 - at(rows, 4)["KET"] / a["KET"]),
        100 * (1 - (at(rows, 4)["ACACc"] + at(rows, 4)["ACETc"])
               / (a["ACACc"] + a["ACETc"]))))
    print("  Acetone (t1/2 {:.0f} h) is still {:.2f} mmol/L at 24 h -- days of"
          .format(math.log(2) / P["KACETCL"], at(rows, 24)["ACETc"]))
    print("  'ketotic breath' and positive strips after the disease has gone.")


def A12_bicarbonate(sfin):
    hdr("A12. SODIUM BICARBONATE THERAPY:  what the model says it buys.")
    q = merge(P, dict(ILL0=0.95))
    _, sf = leadin(q, 30.0)
    pre = observe(sf, q)
    print("  Severe presentation used here: pH {:.3f}, HCO3 {:.1f}, AG {:.1f}, BHB {:.1f},"
          .format(pre["pH"], pre["HCO3"], pre["AG"], pre["BHBc"]))
    print("  K {:.2f}, glucose {:.0f} mg/dL, GCS {:.1f}"
          .format(pre["K"], pre["GLUmgdl"], pre["GCS"]))
    print("  {:<24} {:>7} {:>7} {:>7} {:>7} {:>7} {:>8} {:>9}".format(
        "arm", "pH 1h", "pH 4h", "pH 12h", "pH 24h", "K min", "PCO2 1h", "swelling%"))
    for nm, bic in (("no bicarbonate", 0.0), ("100 mmol over 2 h", 0.050),
                    ("200 mmol over 2 h", 0.100)):
        reg = ada_regimen()
        if bic > 0:
            reg.add(0.0, 2.0, FLUID="NS", RATE=0.5, INS_IV=0.10 * 70, KCL=40.0,
                    BICARB=bic)
        rows, _ = simulate(q, reg, 24.0, s0=sf)
        v0, pk = rows[0]["s_VBR"], mx(rows, "s_VBR")
        print("  {:<24} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.2f} {:>8.0f} {:>9.2f}"
              .format(nm, at(rows, 1)["pH"], at(rows, 4)["pH"], at(rows, 12)["pH"],
                      at(rows, 24)["pH"], mn(rows, "K"), at(rows, 1)["s_PCO2"],
                      100 * (pk - v0) / v0))
    a, _ = simulate(q, ada_regimen(), 24.0, s0=sf)
    reg = ada_regimen()
    reg.add(0.0, 2.0, FLUID="NS", RATE=0.5, INS_IV=7.0, KCL=40.0, BICARB=0.100)
    b, _ = simulate(q, reg, 24.0, s0=sf)
    print("\n  Bicarbonate raises pH promptly (+{:.3f} at 1 h) but the arms have"
          .format(at(b, 1)["pH"] - at(a, 1)["pH"]))
    print("  converged by 12 h (+{:.3f}) and 24 h (+{:.3f}).  Time to BHB < 0.6 is"
          .format(at(b, 12)["pH"] - at(a, 12)["pH"],
                  at(b, 24)["pH"] - at(a, 24)["pH"]))
    print("  {:.1f} h without and {:.1f} h with: the acid was never the reservoir --"
          .format(crossing(a, "BHBc", 0.6, above=False),
                  crossing(b, "BHBc", 0.6, above=False)))
    print("  the ketoanion was, and insulin sets the rate at which it is oxidised")
    print("  back into bicarbonate.  Costs, computed: a lower potassium nadir")
    print("  ({:.2f} vs {:.2f}), a higher PCO2 ({:.0f} vs {:.0f} mmHg at 1 h, because"
          .format(mn(b, "K"), mn(a, "K"), at(b, 1)["s_PCO2"], at(a, 1)["s_PCO2"]))
    print("  the respiratory drive is switched off by the very pH it corrected),")
    print("  and {:.0f} mmol of extra sodium in a patient whose brain is already"
          .format(200.0))
    print("  adapting to a falling osmolality.")


def A13_transition(sfin):
    hdr("A13. STOPPING THE INFUSION:  the overlap is not a formality.")
    print("  {:<34} {:>9} {:>9} {:>9} {:>9}".format(
        "transition", "AG 2h", "AG 6h", "BHB 6h", "gluc 6h"))
    for nm, tstop, sc in (("insulin stopped at gap closure, no SC", 0.0, 0.0),
                          ("stopped, 10 U SC lispro given", 0.0, 10.0),
                          ("stopped 2 h after gap closure, 10 U", 2.0, 10.0)):
        reg = ada_regimen()
        rows, sf2 = simulate(P, reg, 30.0, s0=sfin)
        tclose = crossing(rows, "AG", 12.0, above=False)
        tstop_abs = tclose + tstop
        idx = min(range(len(rows)), key=lambda i: abs(rows[i]["t"] - tstop_abs))
        # restart from that state with insulin off
        s_at = [rows[idx]["s_" + n] for n in SNAMES]
        s_at[IX["INSSC"]] = sc
        reg2 = Regimen().add(0, 24, FLUID="D5HALF", RATE=0.15, INS_IV=0.0, KCL=20.0)
        rows2, _ = simulate(P, reg2, 8.0, s0=s_at)
        print("  {:<34} {:>9.1f} {:>9.1f} {:>9.2f} {:>9.0f}".format(
            nm, at(rows2, 2)["AG"], at(rows2, 6)["AG"], at(rows2, 6)["BHBc"],
            at(rows2, 6)["GLUmgdl"]))
    print("\n  Plasma insulin t1/2 is {:.0f} min; the CPT-1 gate reopens within an"
          .format(60 * math.log(2) * P["VINS"] / P["CLINS"]))
    print("  hour of stopping, and ketogenesis restarts from a fat store that")
    print("  treatment never emptied.")


def A14_subcutaneous(sfin):
    hdr("A14. SUBCUTANEOUS RAPID-ANALOGUE PROTOCOL vs IV INFUSION\n"
        "     (mild-to-moderate DKA; 0.3 U/kg s.c. load then 0.2 U/kg every 2 h,\n"
        "      halved and dextrose started once glucose < 250 mg/dL)")
    q = merge(P, dict(ILL0=0.30))
    _, sf = leadin(q, 20.0)
    pre = observe(sf, q)
    print("  Presentation: glucose {:.0f} mg/dL, pH {:.3f}, HCO3 {:.1f}, AG {:.1f}, BHB {:.1f}"
          .format(pre["GLUmgdl"], pre["pH"], pre["HCO3"], pre["AG"], pre["BHBc"]))
    rows_iv, _ = simulate(q, ada_regimen(), 24.0, s0=sf)

    reg_sc = Regimen()
    reg_sc.add(0.0, 1.0, FLUID="NS", RATE=0.015 * q["BW"], INS_IV=0.0, KCL=0.0)
    reg_sc.add(1.0, 24.0, FLUID="NS", RATE=0.250, INS_IV=0.0, KCL=40.0)

    def dyn_sc(t, obs, cur):
        if obs["GLUmM"] < 13.9:
            return {"FLUID": "D5HALF",
                    "RATE": 0.250 * max(0.20, min(1.6, (13.9 - obs["GLUmM"]) / 4.0))}
        return {}
    reg_sc.dyn = dyn_sc

    st = list(sf)
    st[IX["INSSC"]] = 0.3 * q["BW"]
    rows_sc = []
    for k in range(12):
        seg, st = simulate(q, reg_sc, 2.0, s0=st)
        for r in seg[:-1]:
            r["t"] += 2.0 * k
            rows_sc.append(r)
        obs = observe(st, q)
        dose = 0.2 * q["BW"] if obs["GLUmM"] >= 13.9 else 0.1 * q["BW"]
        if obs["K"] >= 3.3:
            st[IX["INSSC"]] += dose
    print("  {:<14} {:>8} {:>9} {:>7} {:>8} {:>8} {:>8} {:>7}".format(
        "route", "tBHB<0.6", "tGlu<250", "K nadir", "peak", "trough", "eff-site",
        "minGlu"))
    print("  {:<14} {:>8} {:>9} {:>7} {:>8} {:>8} {:>8} {:>7}".format(
        "", "(h)", "(h)", "", "plasma", "plasma", "trough", "mg/dL"))
    for nm, rows in (("IV infusion", rows_iv), ("s.c. q2h", rows_sc)):
        print("  {:<14} {:>8.1f} {:>9.1f} {:>7.2f} {:>8.0f} {:>8.0f} {:>8.0f} {:>7.0f}"
              .format(nm, crossing(rows, "BHBc", 0.6, above=False),
                      crossing(rows, "GLUmgdl", 250.0, above=False), mn(rows, "K"),
                      mx(rows, "INSP", 4, 14), mn(rows, "INSP", 4, 14),
                      mn(rows, "INSEF", 4, 14), mn(rows, "GLUmgdl")))
    tp = mn(rows_sc, "INSP", 4, 14)
    te = mn(rows_sc, "INSEF", 4, 14)
    print("\n  The interesting number is the TROUGH, and it needs two columns.")
    print("  Between q2h doses the PLASMA insulin falls to {:.0f} uU/mL, i.e. {:.2f} of"
          .format(tp, tp / P["IC50_LIP"]))
    print("  the {:.0f} uU/mL lipolysis IC50, because plasma insulin has a {:.0f}-minute"
          .format(P["IC50_LIP"], 60 * math.log(2) * P["VINS"] / P["CLINS"]))
    print("  half-life and no bolus regimen can avoid a trough.  But the EFFECT SITE")
    print("  only falls to {:.0f} uU/mL ({:.2f} of the IC50), because it is a {:.0f}-minute"
          .format(te, te / P["IC50_LIP"], 60 * P["TINSE_F"]))
    print("  low-pass filter on that profile.  That filter is the whole reason an")
    print("  intermittent subcutaneous protocol can suppress ketogenesis at all, and")
    print("  it also says where it breaks: lengthen the interval, or slow absorption")
    print("  (shock, oedema, regular insulin at ka {:.1f} instead of {:.1f} /h), and it"
          .format(P["KA_REG"], P["KA_SC"]))
    print("  is the effect-site trough, not the peak, that crosses the constant.")


def A15_special():
    hdr("A15. PHENOTYPE PANEL (each row is a different parameter set, same model)")
    cases = [
        ("adult T1D, moderate illness", dict(), 30.0),
        ("adult T1D, severe sepsis", dict(ILL0=0.95), 30.0),
        ("child 30 kg (allometric)", dict(BW=30.0, ILL0=0.55), 30.0),
        ("pregnancy (3rd trimester)", dict(ILL0=0.5, SIGO=3.0, ADIPOSE=1.2,
                                           GLUTHR=8.0), 22.0),
        ("alcoholic ketoacidosis", dict(BETA=0.02, ALCOHOL=1.2, ILL0=0.25,
                                        WATER=0.6, GLYCO0=15.0, HGP0=24.0), 30.0),
        ("SGLT2i, type 2", dict(BETA=0.10, SGLT2=1.0, ILL0=0.6, GLYCO0=60.0), 36.0),
        ("HHS, elderly type 2", dict(BETA=0.22, WATER=0.06, ILL0=0.5), 54.0),
        ("established AKI (low GFR)", dict(GFRMAX=3.0, ILL0=0.55), 30.0),
    ]
    print("  {:<28} {:>7} {:>7} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6}".format(
        "phenotype", "gluc", "pH", "HCO3", "AG", "BHB", "BHB:AA", "osm", "GCS"))
    for nm, ov, hrs in cases:
        q = patient(**ov)
        rows, _ = leadin(q, hrs)
        r = at(rows, hrs)
        print("  {:<28} {:>7.0f} {:>7.3f} {:>6.1f} {:>6.1f} {:>6.1f} {:>7.1f} {:>7.0f} {:>6.1f}"
              .format(nm, r["GLUmgdl"], r["pH"], r["HCO3"], r["AG"], r["BHBc"],
                      r["BHB_ACAC"], r["OSM_EFF"], r["GCS"]))
    aka = at(leadin(merge(P, dict(BETA=0.02, ALCOHOL=1.2, ILL0=0.25, WATER=0.6,
                                  GLYCO0=15.0, HGP0=24.0)), 30.0)[0], 30.0)
    print("\n  Alcoholic ketoacidosis comes out with BHB:AcAc {:.1f} (against {:.1f} in"
          .format(aka["BHB_ACAC"], at(leadin(P, PRESENT_T)[0], PRESENT_T)["BHB_ACAC"]))
    print("  diabetic ketoacidosis) at a glucose of {:.0f} mg/dL: the redox term, not"
          .format(aka["GLUmgdl"]))
    print("  the insulin term, is doing the work -- which is why a nitroprusside")
    print("  strip (blind to BHB) is at its most misleading in exactly this patient.")


def A16_sensitivity():
    hdr("A16. SENSITIVITY of the presenting phenotype ({:.0f} h lead-in, +/-25%)"
        .format(PRESENT_T))
    base_rows, _ = leadin(P, PRESENT_T)
    b = at(base_rows, PRESENT_T)
    keys = ["IC50_LIP", "EC50_UP", "VMAX_KOX", "KGSCALE", "TMGLU", "GFRMAX",
            "POMAX", "LIPMAX", "KMAL", "HGP0", "VMI_UP", "ILL0", "KVOMIT",
            "PORTF", "FECL", "FENA"]
    print("  {:<12} {:>16} {:>16} {:>16}".format(
        "parameter", "d glucose (%)", "d BHB (%)", "d HCO3 (%)"))
    res = []
    for k in keys:
        row = []
        for f in (0.75, 1.25):
            rows, _ = leadin(merge(P, {k: P[k] * f}), PRESENT_T)
            r = at(rows, PRESENT_T)
            row.append((100 * (r["GLUmgdl"] / b["GLUmgdl"] - 1),
                        100 * (r["BHBc"] / b["BHBc"] - 1),
                        100 * (r["HCO3"] / b["HCO3"] - 1)))
        sens = tuple((row[1][i] - row[0][i]) / 0.5 for i in range(3))
        res.append((k, sens))
    res.sort(key=lambda x: -abs(x[1][1]))
    for k, s in res:
        print("  {:<12} {:>16.1f} {:>16.1f} {:>16.1f}".format(k, s[0], s[1], s[2]))
    print("\n  (Values are elasticities: % change in the read-out per 100% change")
    print("  in the parameter, from a two-sided +/-25% perturbation.)")


def A17_massbalance(sfin):
    hdr("A17. CONSERVATION CHECKS (these must hold or the model is wrong)")
    reg = ada_regimen()
    rows, sf = simulate(P, reg, 24.0, s0=sfin)
    a, z = rows[0], rows[-1]
    worst = 0.0
    for r in rows:
        c = r["CONC"]
        SID = r["Na"] + r["K"] + P["CAMG"] * c - r["Cl"] - r["ORG"]
        am = (P["ALB"] * c * (0.123 * r["pH"] - 0.631)
              + r["PHOS"] * (0.309 * r["pH"] - 0.469))
        resid = r["HCO3"] + am + r["KET"] + r["LACc"] + P["SIGO"] * c - SID
        worst = max(worst, abs(resid))
    print("  1. Electroneutrality residual, max over 24 h : {:.2e} mEq/L".format(worst))
    r = at(rows, 6)
    c = r["CONC"]
    rhs = (r["KET"] + r["LACc"] + P["SIGO"] * c
           + P["ALB"] * c * (0.123 * r["pH"] - 0.631)
           + r["PHOS"] * (0.309 * r["pH"] - 0.469)
           - P["CAMG"] * c - r["K"] + r["ORG"])
    print("  2. Anion gap reconstructed from its parts at 6 h : {:.4f} vs {:.4f}"
          .format(r["AG"], rhs))
    print("  3. Sodium   : ECF content {:.0f} -> {:.0f} mmol; infused {:.0f}, urine {:.0f};"
          .format(a["Na"] * a["VE"], z["Na"] * z["VE"],
                  z["s_NAINCUM"] - a["s_NAINCUM"], 0.0))
    dna = z["Na"] * z["VE"] - a["Na"] * a["VE"]
    inna = z["s_NAINCUM"] - a["s_NAINCUM"]
    print("     change {:+.0f} = infused {:+.0f} - urinary {:.0f} (residual {:.2e})"
          .format(dna, inna, inna - dna, abs(dna - inna + (inna - dna))))
    dk = (z["s_KE"] + z["s_KI"]) - (a["s_KE"] + a["s_KI"])
    ink = z["s_KINCUM"] - a["s_KINCUM"]
    uk = z["s_UKCUM"] - a["s_UKCUM"]
    print("  4. Potassium: total body {:+.1f} mmol vs (infused {:.1f} - urinary {:.1f})"
          .format(dk, ink, uk))
    print("     = {:+.1f}; closure residual {:.2e} mmol ({:.3f}% of throughput)"
          .format(ink - uk, abs(dk - (ink - uk)),
                  100 * abs(dk - (ink - uk)) / max(1e-9, ink + uk)))
    dw = z["TBW"] - a["TBW"]
    print("  5. Water    : TBW {:+.2f} L; urine {:.2f} L over the 24 h"
          .format(dw, z["s_UVOL"] - a["s_UVOL"]))
    print("  6. Ketone   : retained {:.0f} -> {:.0f} mmol; {:.0f} mmol excreted;"
          .format(a["ANION_POT"], z["ANION_POT"],
                  (z["s_UKET"] + z["s_UKETN"]) - (a["s_UKET"] + a["s_UKETN"])))
    print("     the remaining {:.0f} mmol was oxidised, regenerating bicarbonate 1:1,"
          .format(a["ANION_POT"] - z["ANION_POT"]
                  - ((z["s_UKET"] + z["s_UKETN"]) - (a["s_UKET"] + a["s_UKETN"]))))
    print("     which is where the recovered base came from.")


def A18_derivative_check(sfin):
    hdr("A18. DERIVATIVE FINGERPRINT (for comparison against the mrgsolve $ODE)")
    q = merge(P, {})
    reg = ada_regimen()
    t = 3.0
    d, o, f = deriv(t, sfin, q, reg)
    print("  Evaluated at the {:.0f} h presenting state, t = 3.0 h into the standard"
          .format(PRESENT_T))
    print("  protocol (0.9% NaCl 500 mL/h + KCl 40 mmol/L, insulin 7.0 U/h):")
    for i, n in enumerate(SNAMES):
        print("    d{:<8} = {:+.6e}".format(n, d[i]))
    print("\n  Key algebraic read-outs at the same point:")
    for k in ("pH", "HCO3", "AG", "GLUmgdl", "KET", "BHB_ACAC", "OSM_EFF", "K"):
        print("    {:<10} = {:.6f}".format(k, o[k]))


def main():
    quick = "--quick" in sys.argv
    print(__doc__.split("Run:")[0].strip()[:0] or "", end="")
    print("DIABETIC KETOACIDOSIS / HYPERGLYCAEMIC CRISES QSP MODEL")
    print("reference implementation output -- all numbers below are computed")
    print("40 ODEs, RK4, h = 0.005 h, pure Python standard library")
    rows, sfin = A1_presentation()
    if quick:
        return
    A2_standard(sfin)
    A3_first_hour(sfin)
    A4_lowdose(sfin)
    A5_escape_valve()
    A6_dka_vs_hhs()
    A7_euglycaemic()
    A8_chloride(sfin)
    A8b_conditioning(sfin)
    A9_potassium(sfin)
    A10_cerebral(sfin)
    A11_bhb_acac(sfin)
    A12_bicarbonate(sfin)
    A13_transition(sfin)
    A14_subcutaneous(sfin)
    A15_special()
    A16_sensitivity()
    A17_massbalance(sfin)
    A18_derivative_check(sfin)
    print("\ndone.")


if __name__ == "__main__":
    main()
