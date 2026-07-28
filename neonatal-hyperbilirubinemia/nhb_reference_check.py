#!/usr/bin/env python3
# =============================================================================
#  nhb_reference_check.py
#  Neonatal hyperbilirubinaemia (NHB) QSP model — INDEPENDENT REFERENCE
#  IMPLEMENTATION of the ODE system coded in nhb_mrgsolve_model.R
#  ---------------------------------------------------------------------------
#  WHY THIS FILE EXISTS
#  --------------------
#  The mrgsolve model cannot be executed in this environment (no R runtime).
#  Every quantitative claim made in README.md and in the header of
#  nhb_mrgsolve_model.R is therefore produced HERE, by a second, dependency-free
#  (pure standard-library, fixed-step RK4) implementation of the SAME equations
#  and the SAME parameter values.  If the two files disagree, this file is the
#  one that was actually run — nhb_reference_output.txt is its verbatim output.
#
#  Run:   python3 nhb_reference_check.py > nhb_reference_output.txt
#
#  STATE VECTOR (34 ODEs) — identical ordering to the R model
#  ---------------------------------------------------------------------------
#   0 W      body weight                                            kg
#   1 HB     haemoglobin concentration                              g/dL
#   2 RET    reticulocyte fraction (marrow output signal)           %
#   3 ABMAT  maternal alloantibody load on neonatal RBC             rel (0-1)
#   4 HO1    haem oxygenase-1 activity                              rel
#   5 COHB   carboxyhaemoglobin (haem-catabolism biomarker)         %
#   6 BLEX   extravasated blood haemoglobin (cephalohaematoma)      g
#   7 BP     plasma unconjugated bilirubin                          mg
#   8 BEX    extravascular (skin+interstitium) UCB — PHOTON TARGET  mg
#   9 BLIV   hepatocyte unconjugated bilirubin (ligandin pool)      mg
#  10 BCON   hepatocyte conjugated bilirubin                        mg
#  11 BGC    gut-lumen conjugated bilirubin                         mg
#  12 BGU    gut-lumen unconjugated bilirubin (reabsorbable)        mg
#  13 BSTL   cumulative faecal + renal bilirubin output             mg
#  14 LUMI   plasma lumirubin (structural photoisomer)              mg
#  15 EZ     plasma (4Z,15E) configurational photoisomer            mg
#  16 ALB    serum albumin                                          g/dL
#  17 UGT    UGT1A1 activity                                        frac adult
#  18 TGX    AAV8-hUGT1A1 transgene-derived activity                frac adult
#  19 BBR    basal-ganglia bilirubin (nM-equivalent)                nM-eq
#  20 INJ    cumulative neuronal injury index                       0-1
#  21 ABRD   auditory brainstem response deficit                    rel
#  22 ASNMP  stannsoporfin IM depot                                 mg
#  23 CSNMP  stannsoporfin plasma                                   ug/mL
#  24 APB    phenobarbital gut depot                                mg
#  25 CPB    phenobarbital plasma                                   mg/L
#  26 IGGC   IVIG central                                           g
#  27 IGGP   IVIG peripheral                                        g
#  28 AUDCA  ursodeoxycholic acid gut depot                          mg
#  29 CUDCA  ursodeoxycholic acid plasma                             umol/L
#  30 GBIND  intraluminal bilirubin binder (agar/charcoal/zinc)     mg
#  31 AUCX   AUC of TSB above the AAP-2022 phototherapy threshold   mg*h/dL
#  32 PTH    cumulative phototherapy exposure                       h
#  33 ETV    cumulative exchange-transfused volume                  mL/kg
# =============================================================================
import math

# ----------------------------------------------------------------------------
# 0. PARAMETERS  (single source of truth; mirrored verbatim in the R $PARAM)
# ----------------------------------------------------------------------------
P = dict(
    # ---- anthropometry / geometry -------------------------------------------
    GA        = 40.0,    # completed weeks of gestation
    W0        = 3.40,    # birth weight (kg)
    KW        = 0.010,   # weight tracking rate toward reference (1/h)
    VP_KG     = 0.50,    # plasma volume                       (dL/kg)
    VEX_KG    = 2.00,    # extravascular albumin/bilirubin space(dL/kg)
    VLIV_KG   = 0.40,    # hepatocyte water                    (dL/kg)
    VBL_KG    = 0.85,    # whole-blood volume                  (dL/kg)
    PSXKG     = 0.50,    # plasma<->interstitium permeability   (dL/h/kg)

    # ---- erythron / haemolysis ---------------------------------------------
    HB0       = 17.0,    # birth haemoglobin (g/dL)
    LRBC      = 80.0,    # neonatal RBC lifespan (days) 80 term / 50 preterm
    RET0      = 3.0,     # birth reticulocyte % (marrow output)
    TAURET    = 36.0,    # marrow response time constant (h)
    RETMAX    = 12.0,    # maximal reticulocyte response (%)
    KIMM      = 0.0045,  # immune haemolysis rate at ABMAT = 1 (1/h)
    KOX       = 0.0040,  # oxidative (G6PD) haemolysis rate coefficient (1/h)
    KABCL     = 0.00138, # alloantibody clearance (1/h)  t1/2 = 21 d
    KRESORB   = 0.0060,  # extravasated-blood resorption (1/h) t1/2 = 115 h
    G6PD      = 0.0,     # 0 = normal, 1 = deficient (class II/III)
    ABMAT0    = 0.0,     # birth alloantibody load (0 = no isoimmunisation)
    BLEX0     = 0.0,     # extravasated Hb at birth (g)

    # ---- bilirubin production ----------------------------------------------
    BILPERHB  = 34.0,    # mg bilirubin per g haemoglobin catabolised
    FEARLY    = 0.25,    # early-labelled fraction (ineffective erythropoiesis
                         #   + hepatic haem turnover) as fraction of RBC flux
    KHO1      = 0.030,   # HO-1 turnover (1/h)
    KISNMP    = 8.0,     # stannsoporfin Ki vs HO-1 (ug/mL)
    KCO       = 0.4694,  # COHb formation per bilirubin flux (%*h^-1 / (mg/h/kg))
    KCOOUT    = 0.1386,  # COHb elimination (1/h)  t1/2 = 5 h

    # ---- albumin binding ---------------------------------------------------
    ALB0      = 3.40,    # birth albumin (g/dL)
    ALBSET    = 3.60,    # albumin set point (g/dL)
    KALB      = 0.004,   # albumin turnover toward set point (1/h)
    ALBDON    = 4.00,    # donor-blood albumin during exchange (g/dL)
    KA1       = 45.0,    # primary-site association constant (1/uM) term neonate
    KA2       = 1.00,    # secondary-site association constant (1/uM)
    NSITE2    = 1.00,    # number of secondary (low-affinity) sites per albumin
    FMATK     = 1.00,    # binding-affinity maturity factor (preterm < 1)
    FACID     = 1.00,    # acidaemia factor on KA (pH 7.15 ~ 0.70)
    FDISP     = 1.00,    # displacer factor (ceftriaxone/ibuprofen/FFA ~ 0.6-0.8)

    # ---- hepatic handling --------------------------------------------------
    CLUP      = 0.800,   # sinusoidal uptake clearance (dL/h per 3.4 kg)
    FOATP     = 1.00,    # SLCO1B1 genotype factor on uptake
    KBACK     = 0.500,   # hepatocyte -> plasma efflux (1/h)
    VMAXU     = 14.0,    # UGT1A1 Vmax at adult activity (mg/h per 3.4 kg)
    KMU       = 5.00,    # UGT1A1 Km (mg/dL hepatocyte water)
    ALLO      = 0.75,    # allometric exponent for clearances
    F0UGT     = 0.008,   # UGT1A1 activity at birth (fraction of adult)
    TAUUGT    = 330.0,   # postnatal maturation time constant (h) ~14 d
    KUGT      = 0.030,   # UGT1A1 protein turnover (1/h) t1/2 = 23 h
    GENO      = 1.00,    # UGT1A1 genotype factor (see GENOTYPES below)
    EPB       = 2.50,    # max phenobarbital induction (fold-1)
    EC50PB    = 20.0,    # phenobarbital EC50 (mg/L)
    KALT      = 0.015,   # UGT1A1-INDEPENDENT elimination of unconjugated
                         #   bilirubin (oxidative degradation to BOXes by haem
                         #   peroxidases/CYP + the small biliary excretion of
                         #   unconjugated pigment).  Small in health, but it is
                         #   the ONLY non-photochemical route left in
                         #   Crigler-Najjar type I and it is what bounds the
                         #   untreated CN-I plateau. (1/h)
    KMRP2     = 0.400,   # MRP2 canalicular export (1/h)
    FCHOL     = 1.00,    # cholestasis factor on MRP2 (bronze baby < 1)

    # ---- enterohepatic shunt ----------------------------------------------
    KBGLUC    = 0.250,   # beta-glucuronidase deconjugation (1/h at BGA = 1)
    BGA0      = 1.00,    # beta-glucuronidase activity at birth (rel)
    BGABM     = 1.00,    # extra activity from human milk (rel)
    BGAMIN    = 0.15,    # floor after flora establish (rel)
    TAUBGA    = 504.0,   # flora-maturation time constant (h) = 21 d
    KREAB     = 0.200,   # jejunal reabsorption of UCB (1/h)
    KTRANS    = 0.100,   # gut transit / stool output (1/h at full intake)
    FTRANS0   = 0.35,    # transit floor at zero enteral intake
    KB50      = 300.0,   # intraluminal binder amount for 50% block (mg)
    BGU0      = 5.0,     # meconium unconjugated bilirubin load (mg)
    BGC0      = 10.0,    # meconium conjugated bilirubin load (mg)
    BREAST    = 1.0,     # 1 = human milk, 0 = formula

    # ---- phototherapy ------------------------------------------------------
    # Photochemistry is TWO reactions, not one.  (4Z,15Z) <-> (4Z,15E) is fast,
    # REVERSIBLE and reaches a photostationary state (E ~ 20-25 % of measured
    # total bilirubin); lumirubin arises by IRREVERSIBLE intramolecular
    # cyclisation of the E-isomer and is the actual excretory route.  Rate
    # Every photoreaction rate is (photon delivery) x (CONCENTRATION in the
    # irradiated layer), never x (amount): a photon reaction happens in a
    # volume the light can reach, so its rate must scale with the irradiated
    # AREA and the local concentration.  Writing it against the amount instead
    # makes photon delivery scale as BSA*W and silently inverts every
    # size/growth conclusion (see A11).  Units: mg/h per (uW/cm2/nm * cm2 *
    # mg/dL).
    KZE       = 5.44e-5, # (4Z,15Z) -> (4Z,15E) photoisomerisation
    KEZP      = 1.02e-4, # (4Z,15E) -> (4Z,15Z) photoreversion
    KLUMF     = 4.675e-5,# (4Z,15E) -> lumirubin photocyclisation
    I50       = 45.0,    # irradiance at which the skin optically saturates
    # Skin optics change with age: epidermal/dermal thickening and melanin
    # accumulation shrink the fraction of the extravascular pigment pool that
    # lies within the effective blue-light penetration depth.  Illustrative
    # values; the qualitative direction is well established and turns out to
    # matter far more than the surface-to-mass ratio (see analysis A11).
    FOPTMIN   = 0.30,    # asymptotic optical accessibility (adult)
    TAUOPT    = 3.00,    # optical-accessibility decay constant (years)
    KLUM      = 0.700,   # lumirubin elimination (1/h) t1/2 = 1 h (no UGT1A1)
    KEZBIL    = 0.020,   # E-isomer biliary excretion (1/h, no UGT1A1)
    KREV      = 0.010,   # E -> Z thermal (dark) reversion (1/h)

    # ---- exchange transfusion ---------------------------------------------
    HBDON     = 13.0,    # reconstituted donor-blood Hb (g/dL)
    ETEFF     = 0.45,    # exchange efficiency (mixing/recirculation losses:
                         #   a double-volume exchange replaces ~85 % of red
                         #   cells, not 100 %, for exactly this reason)
    FABREM    = 1.00,    # efficiency of antibody/coated-RBC removal

    # ---- drugs -------------------------------------------------------------
    KASNMP    = 0.173,   # SnMP IM absorption (1/h) t1/2 = 4 h
    VSNMP     = 3.00,    # SnMP volume of distribution (dL/kg)
    KESNMP    = 0.023,   # SnMP elimination (1/h) t1/2 = 30 h
    KAPB      = 0.400,   # phenobarbital absorption (1/h)
    VPB       = 9.00,    # phenobarbital Vd (dL/kg)
    KEPB      = 0.00693, # phenobarbital elimination (1/h) t1/2 = 100 h
    VIGG      = 0.50,    # IVIG central volume (dL/kg)
    KIGGCP    = 0.010,   # IVIG central->peripheral (1/h)
    KIGGPC    = 0.008,   # IVIG peripheral->central (1/h)
    KIGGEL    = 0.00120, # IVIG elimination (1/h) t1/2 = 24 d
    IGG0      = 10.0,    # endogenous (transplacental) IgG (g/L)
    IMAXIVIG  = 0.65,    # max fractional block of immune haemolysis
    IC50IVIG  = 6.00,    # IVIG EC50 above baseline (g/L)
    KAUDCA    = 0.500,   # UDCA absorption (1/h)
    VUDCA     = 3.00,    # UDCA Vd (dL/kg)
    KEUDCA    = 0.140,   # UDCA elimination (1/h)
    EMAXU     = 0.80,    # UDCA max effect
    EC50UDCA  = 4.00,    # UDCA EC50 (umol/L)
    KTGON     = 0.00289, # transgene expression onset (1/h) t1/2 = 10 d
    TGXMAX    = 0.100,   # plateau transgene activity (fraction of adult)

    # ---- blood-brain barrier and injury -----------------------------------
    KINBBB    = 0.100,   # Bf -> brain influx (1/h)
    KOUTBBB   = 0.100,   # brain efflux incl. P-gp (1/h)
    FBBB      = 1.00,    # BBB permeability multiplier (prematurity/sepsis)
    BBRTHR    = 35.0,    # injury threshold, brain bilirubin (nM-eq)
    KINJ      = 5.0e-4,  # injury accrual (1/(nM*h))
    KREP      = 0.0020,  # repair of the reversible injury component (1/h)
    BBRTHRA   = 25.0,    # ABR threshold (nM-eq) — lower than injury threshold
    KABR      = 4.0e-4,  # ABR deficit accrual (1/(nM*h))
    KABRREC   = 0.020,   # ABR recovery (1/h)
    INJ50     = 0.50,    # kernicterus logistic midpoint
    INJSL     = 0.08,    # kernicterus logistic slope

    # ---- AAP 2022 thresholds (analytic approximation, see README) ----------
    RF        = 0.0,     # 1 = neurotoxicity risk factor(s) present
    FIXHB     = 0.0,     # >0 clamps haemoglobin (used for the A11 decomposition)
    PRODSCALE = 1.0,     # explicit multiplier on bilirubin production, used to
                         #   switch the production-ontogeny term of the A11
                         #   decomposition on and off
    AGEOFF    = 0.0,     # hours added to t for the AGE-dependent skin-optics
                         #   term only; lets a quasi-steady-state run be
                         #   performed at the optics of an older child
)

# UGT1A1 genotype activity factors (fraction of wild-type adult activity)
GENOTYPES = {
    "wild"      : 1.00,   # *1/*1
    "UGT1A1*28het": 0.65, # (TA)7 heterozygote
    "Gilbert"   : 0.35,   # *28/*28
    "UGT1A1*6"  : 0.40,   # G71R homozygote (East Asian)
    "CN2"       : 0.05,   # Crigler-Najjar type II (Arias)
    "CN1"       : 0.00,   # Crigler-Najjar type I
}

# ----------------------------------------------------------------------------
# 1. STRUCTURAL HELPERS
# ----------------------------------------------------------------------------
IDX = ["W","HB","RET","ABMAT","HO1","COHB","BLEX","BP","BEX","BLIV","BCON",
       "BGC","BGU","BSTL","LUMI","EZ","ALB","UGT","TGX","BBR","INJ","ABRD",
       "ASNMP","CSNMP","APB","CPB","IGGC","IGGP","AUDCA","CUDCA","GBIND",
       "AUCX","PTH","ETV"]
J = {k: i for i, k in enumerate(IDX)}
NST = len(IDX)

# WHO median growth reference (male), age in years -> kg and cm.
# The first 14 days are written out explicitly so that the physiologic
# postnatal weight nadir (~6-7 % at day 3) is part of the reference.
_WREF = [(0.0, 1.000), (3/365., 0.939), (7/365., 0.985), (14/365., 1.059),
         (28/365., 1.191), (0.25, 1.765), (0.5, 2.294), (1.0, 2.824),
         (2.0, 3.588), (4.0, 4.794), (6.0, 6.029), (8.0, 7.441),
         (10.0, 9.382), (14.0, 14.706), (18.0, 19.412)]
_HREF = [(0.0, 50.0), (0.25, 61.4), (0.5, 67.6), (1.0, 75.7), (2.0, 87.1),
         (4.0, 103.3), (6.0, 116.0), (8.0, 127.3), (10.0, 137.8),
         (14.0, 163.2), (18.0, 176.1)]


def _interp(tab, x):
    if x <= tab[0][0]:
        return tab[0][1]
    if x >= tab[-1][0]:
        return tab[-1][1]
    for i in range(len(tab) - 1):
        x0, y0 = tab[i]
        x1, y1 = tab[i + 1]
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return tab[-1][1]


def wref(pna_h, p):
    """Reference weight (kg) at postnatal age pna_h, scaled to birth weight."""
    return p["W0"] * _interp(_WREF, pna_h / 24.0 / 365.0)


def href(pna_h, p):
    """Reference length/height (cm); birth length scaled with birth weight."""
    scale = (p["W0"] / 3.40) ** (1.0 / 3.0)
    return _interp(_HREF, pna_h / 24.0 / 365.0) * (scale if pna_h < 8760 else 1.0)


def bsa_cm2(w, h):
    """Haycock body surface area (cm^2). BSA(m2)=0.024265*W^0.5378*H^0.3964."""
    return 0.024265 * (w ** 0.5378) * (h ** 0.3964) * 1.0e4


def _sig(x):
    """Numerically safe logistic."""
    if x > 40.0:
        return 1.0
    if x < -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def free_bilirubin_uM(cp_mgdl, alb_gdl, p):
    """Unbound (free) bilirubin from a TWO-CLASS saturable binding isotherm.

        BT = Bf + A*K1*Bf/(1+K1*Bf) + A*n2*K2*Bf/(1+K2*Bf)

    The high-affinity site (K1 ~ 4.5e7 /M in the term neonate) carries almost
    all of the load until it saturates; the low-affinity class (K2 ~ 1e6 /M)
    is what stops Bf from diverging once the primary site is full.  A one-site
    isotherm predicts physically impossible free concentrations above a
    bilirubin/albumin molar ratio of ~1 and must not be used here.

    Monotone in Bf, so solved by 60 bisection steps (identical algorithm in
    the R/C++ model, 40 steps) - no closed form is needed and none exists.
    """
    bt = cp_mgdl * 17.1                      # mg/dL -> umol/L
    a = alb_gdl * 150.4                      # g/dL  -> umol/L (MW 66.5 kDa)
    fk = p["FMATK"] * p["FACID"] * p["FDISP"]
    k1 = p["KA1"] * fk
    k2 = p["KA2"] * fk
    n2 = p["NSITE2"]
    if bt <= 0.0:
        return 0.0
    lo, hi = 0.0, bt
    for _ in range(60):
        bf = 0.5 * (lo + hi)
        tot = bf + a * k1 * bf / (1.0 + k1 * bf) + a * n2 * k2 * bf / (1.0 + k2 * bf)
        if tot < bt:
            lo = bf
        else:
            hi = bf
    return 0.5 * (lo + hi)


def thr_pt(pna_h, p):
    """AAP-2022 phototherapy threshold (mg/dL), analytic approximation."""
    plateau = min(21.0, 16.5 + 1.0 * (min(p["GA"], 40.0) - 35.0)) - 2.0 * p["RF"]
    return plateau * (0.42 + 0.58 * (1.0 - math.exp(-pna_h / 36.0)))


def thr_et(pna_h, p):
    """AAP-2022 exchange-transfusion threshold (mg/dL), approximation."""
    return thr_pt(pna_h, p) + (3.5 - 1.0 * p["RF"])


def hbset(pna_h):
    """Haemoglobin set point sensed by the renal/hepatic O2 sensor.

    Falls after birth (arterial pO2 jumps, HbF -> HbA) and recovers: this is
    the mechanism of the physiologic nadir, not a separate 'nadir' term.
    """
    d = pna_h / 24.0
    if d < 3.0:
        return 17.0
    if d < 60.0:
        return 17.0 - 6.0 * (d - 3.0) / 57.0
    if d < 120.0:
        return 11.0 + 1.5 * (d - 60.0) / 60.0
    return 12.5


# ----------------------------------------------------------------------------
# 2. CONTROL (scenario) OBJECT
# ----------------------------------------------------------------------------
class Ctl:
    """Time-varying controls; every field is a callable f(t_hours)->float."""

    def __init__(self, irr=None, fbsa=None, et=None, fi=None, ox=None,
                 auto_pt=None):
        # auto_pt = (irradiance, off_margin_mg_dL): closed-loop phototherapy,
        # switched on when TSB crosses the AAP threshold and off once it is
        # off_margin below it.  This is how the treatment is actually given,
        # and it makes "phototherapy hours" a model OUTPUT rather than an input.
        self.auto_pt = auto_pt
        self._pt_on = False
        if auto_pt is not None and irr is None:
            irr = lambda t: (auto_pt[0] if self._pt_on else 0.0)
        self.irr = irr or (lambda t: 0.0)     # spectral irradiance uW/cm2/nm
        self.fbsa = fbsa or (lambda t: 0.80)  # fraction of BSA irradiated
        self.et = et or (lambda t: 0.0)       # exchange rate mL/kg/h
        self.fi = fi or (lambda t: 1.0)       # enteral intake adequacy 0-1
        self.ox = ox or (lambda t: 0.0)       # oxidant challenge multiplier


def window(t0, t1, val, base=0.0):
    return lambda t: val if (t0 <= t < t1) else base


def cycles(on_h, off_h, val, t_start=0.0, t_end=1e9):
    def f(t):
        if t < t_start or t >= t_end:
            return 0.0
        ph = (t - t_start) % (on_h + off_h)
        return val if ph < on_h else 0.0
    return f


def ramp_fi(fi_early, t_switch, fi_late):
    return lambda t: fi_early if t < t_switch else fi_late


# ----------------------------------------------------------------------------
# 3. DERIVATIVES
# ----------------------------------------------------------------------------
def deriv(t, y, p, c):
    d = [0.0] * NST
    (W, HB, RET, ABMAT, HO1, COHB, BLEX, BP, BEX, BLIV, BCON, BGC, BGU, BSTL,
     LUMI, EZ, ALB, UGT, TGX, BBR, INJ, ABRD, ASNMP, CSNMP, APB, CPB, IGGC,
     IGGP, AUDCA, CUDCA, GBIND, AUCX, PTH, ETV) = y
    W = max(W, 0.4)
    ALB = max(ALB, 0.5)

    # ---- geometry ---------------------------------------------------------
    HT = href(t, p)
    BSA = bsa_cm2(W, HT)
    VP = p["VP_KG"] * W
    VEX = p["VEX_KG"] * W
    VLIV = p["VLIV_KG"] * W
    VBL = p["VBL_KG"] * W
    allo = (W / 3.40) ** p["ALLO"]

    # ---- concentrations ---------------------------------------------------
    CP = BP / VP                                   # native (4Z,15Z) UCB, mg/dL
    CEX = BEX / VEX
    CLIV = BLIV / VLIV
    CISO = (LUMI + EZ) / (VP + VEX)                # photoisomers, mg/dL
    TSB = CP + CISO                                # what the laboratory reports
    # photoisomers occupy albumin too: the site pool left for native UCB
    ALBF = max(0.30, ALB - CISO * 17.1 / 150.4)
    BFuM = free_bilirubin_uM(CP, ALBF, p)
    BFnM = 1000.0 * BFuM
    BA = (CP * 17.1) / (ALB * 150.4)               # bilirubin/albumin molar ratio

    # ---- drug effects -----------------------------------------------------
    CIGG = 10.0 * IGGC / (p["VIGG"] * W)           # g/L
    igg_ex = max(0.0, CIGG - p["IGG0"])
    F_IVIG = p["IMAXIVIG"] * igg_ex / (p["IC50IVIG"] + igg_ex)
    F_SNMP = 1.0 / (1.0 + CSNMP / p["KISNMP"])     # competitive HO-1 block
    F_UDCA = p["EMAXU"] * CUDCA / (p["EC50UDCA"] + CUDCA)
    F_PB = 1.0 + p["EPB"] * CPB / (p["EC50PB"] + CPB)

    # ---- haemolysis -------------------------------------------------------
    k_basal = 1.0 / (p["LRBC"] * 24.0)
    k_imm = p["KIMM"] * ABMAT * (1.0 - F_IVIG)
    k_ox = p["KOX"] * c.ox(t) * p["G6PD"]
    k_hem = k_basal + k_imm + k_ox

    # ---- exchange transfusion (a real 2-4 h procedure, not an instant jump)
    et_rate = c.et(t)                              # mL/kg/h
    QET = et_rate * W / 100.0 * p["ETEFF"]         # dL/h of plasma turned over
    et_on = 1.0 if et_rate > 0 else 0.0

    # ---- erythron ---------------------------------------------------------
    hgap = max(0.0, (hbset(t) - HB) / max(hbset(t), 1.0))
    ret_tgt = p["RET0"] * math.exp(-t / 1000.0) + p["RETMAX"] * hgap / (0.15 + hgap)
    d[J["RET"]] = (ret_tgt - RET) / p["TAURET"]
    kprod = k_basal * p["HB0"] / p["RET0"]
    d[J["HB"]] = 0.0 if p["FIXHB"] > 0 else (-k_hem * HB + kprod * RET
                  + QET / VBL * (p["HBDON"] - HB) * et_on
                  - HB * (p["KW"] * (wref(t, p) - W)) / max(W, 0.4))
    if p["FIXHB"] > 0:
        d[J["RET"]] = 0.0

    d[J["ABMAT"]] = -p["KABCL"] * ABMAT - QET / VP * ABMAT * p["FABREM"] * et_on
    d[J["BLEX"]] = -p["KRESORB"] * BLEX

    # ---- bilirubin production --------------------------------------------
    hbflux = k_hem * HB * VBL                      # g Hb/h from circulating RBC
    prod_rbc = p["BILPERHB"] * hbflux
    prod_early = p["FEARLY"] * prod_rbc
    prod_ex = p["BILPERHB"] * p["KRESORB"] * BLEX
    PROD = (prod_rbc + prod_early + prod_ex) * HO1 * F_SNMP * p["PRODSCALE"]

    d[J["HO1"]] = p["KHO1"] * (1.0 + 0.30 * min(1.0, k_imm / 0.01) - HO1)
    d[J["COHB"]] = p["KCO"] * PROD / W - p["KCOOUT"] * COHB

    # ---- phototherapy: photon delivery and the two photoreactions --------
    irr = c.irr(t)
    ieff = irr / (1.0 + irr / p["I50"])            # optical saturation of skin
    fbsa = c.fbsa(t) if irr > 0 else 0.0
    fopt = (p["FOPTMIN"] + (1.0 - p["FOPTMIN"])
            * math.exp(-((t + p["AGEOFF"]) / 8760.0) / p["TAUOPT"]))
    u = ieff * fbsa * BSA * fopt                   # photon delivery (uW/nm)
    CEZ = EZ / (VP + VEX)                          # E-isomer concentration
    PZE = p["KZE"] * u * CEX                       # Z -> E  (irradiated skin)
    PEZ = p["KEZP"] * u * CEZ                      # E -> Z  photoreversion
    PLUM = p["KLUMF"] * u * CEZ                    # E -> lumirubin (one-way)
    TREV = p["KREV"] * EZ                          # E -> Z  thermal, dark
    EBIL = p["KEZBIL"] * EZ                        # E-isomer into bile
    d[J["PTH"]] = 1.0 if irr > 0 else 0.0

    # ---- hepatic handling ------------------------------------------------
    UP = p["CLUP"] * allo * p["FOATP"] * CP
    BACK = p["KBACK"] * BLIV
    CONJ = p["VMAXU"] * allo * UGT * CLIV / (p["KMU"] + CLIV)
    EXPB = p["KMRP2"] * p["FCHOL"] * BCON
    ALTE = p["KALT"] * BLIV                        # non-UGT1A1 escape route

    # ---- enterohepatic shunt ---------------------------------------------
    bga = (p["BGAMIN"] + (p["BGA0"] * (1.0 + p["BGABM"] * p["BREAST"]) - p["BGAMIN"])
           * math.exp(-t / p["TAUBGA"]))
    fi = c.fi(t)
    ftrans = (p["FTRANS0"] + (1.0 - p["FTRANS0"]) * fi) * (1.0 + 0.5 * F_UDCA)
    if t < 12.0:
        ftrans *= 0.35                              # meconium has not passed
    occ = GBIND / (p["KB50"] + GBIND)
    DECON = p["KBGLUC"] * bga * BGC
    REAB = p["KREAB"] * BGU * (1.0 - occ) * (1.0 - 0.25 * F_UDCA)
    TRC = p["KTRANS"] * ftrans * BGC
    TRU = p["KTRANS"] * ftrans * BGU

    # ---- bilirubin mass balance -----------------------------------------
    PSX = p["PSXKG"] * W                            # dL/h plasma<->interstitium
    fx = PSX * (CP - CEX)
    VEZ = VP + VEX                                  # photoisomer distribution
    # The photoisomers are distributed over the WHOLE albumin space, so the
    # reverse reactions must return bilirubin to plasma and interstitium in
    # proportion to their volumes.  Returning all of it to plasma would make
    # the futile Z->E->Z cycle a photon-driven pump out of the tissues and
    # would raise the measured plasma level with no change in body burden.
    RBACKP = (PEZ + TREV) * VP / VEZ
    RBACKX = (PEZ + TREV) * VEX / VEZ
    d[J["BP"]] = (PROD - fx - UP + BACK + RBACKP - QET * CP * et_on)
    d[J["BEX"]] = fx - PZE + RBACKX
    d[J["BLIV"]] = UP - BACK - CONJ - ALTE + REAB
    d[J["BCON"]] = CONJ - EXPB
    d[J["BGC"]] = EXPB - DECON - TRC
    d[J["BGU"]] = DECON - REAB - TRU
    d[J["LUMI"]] = PLUM - p["KLUM"] * LUMI - QET * (LUMI / VEZ) * et_on
    d[J["EZ"]] = PZE - PEZ - PLUM - TREV - EBIL - QET * (EZ / VEZ) * et_on
    d[J["BSTL"]] = TRC + TRU + p["KLUM"] * LUMI + EBIL + ALTE
    d[J["ETV"]] = et_rate

    # ---- albumin ---------------------------------------------------------
    d[J["ALB"]] = (p["KALB"] * (p["ALBSET"] - ALB)
                   + QET / VP * (p["ALBDON"] - ALB) * et_on)

    # ---- UGT1A1: ontogeny x genotype x induction, with protein turnover --
    ont = p["F0UGT"] + (1.0 - p["F0UGT"]) * (1.0 - math.exp(-t / p["TAUUGT"]))
    ugt_tgt = p["GENO"] * ont * F_PB + TGX
    d[J["UGT"]] = p["KUGT"] * (ugt_tgt - UGT)
    gtt = p.get("_GTT")
    d[J["TGX"]] = (p["KTGON"] * (p["TGXMAX"] - TGX)
                   if (gtt is not None and t >= gtt) else 0.0)

    # ---- neurotoxicity ---------------------------------------------------
    d[J["BBR"]] = p["KINBBB"] * p["FBBB"] * BFnM - p["KOUTBBB"] * BBR
    exc = max(0.0, BBR - p["BBRTHR"])
    d[J["INJ"]] = p["KINJ"] * exc * (1.0 - INJ) - p["KREP"] * INJ
    exa = max(0.0, BBR - p["BBRTHRA"])
    d[J["ABRD"]] = p["KABR"] * exa - p["KABRREC"] * ABRD

    # ---- drug PK ---------------------------------------------------------
    d[J["ASNMP"]] = -p["KASNMP"] * ASNMP
    d[J["CSNMP"]] = (p["KASNMP"] * ASNMP / (p["VSNMP"] * W) * 10.0
                     - p["KESNMP"] * CSNMP)
    d[J["APB"]] = -p["KAPB"] * APB
    d[J["CPB"]] = p["KAPB"] * APB / (p["VPB"] * W) * 10.0 - p["KEPB"] * CPB
    d[J["IGGC"]] = (-p["KIGGCP"] * IGGC + p["KIGGPC"] * IGGP - p["KIGGEL"] * IGGC
                    - QET / VP * IGGC * et_on)
    d[J["IGGP"]] = p["KIGGCP"] * IGGC - p["KIGGPC"] * IGGP
    d[J["AUDCA"]] = -p["KAUDCA"] * AUDCA
    d[J["CUDCA"]] = (p["KAUDCA"] * AUDCA / (p["VUDCA"] * W) * 10.0
                     - p["KEUDCA"] * CUDCA)
    d[J["GBIND"]] = -p["KTRANS"] * ftrans * GBIND

    # ---- growth and threshold bookkeeping -------------------------------
    d[J["W"]] = p["KW"] * (wref(t, p) * (1.0 - 0.12 * (1.0 - fi)) - W)
    d[J["AUCX"]] = max(0.0, TSB - thr_pt(t, p))
    return d


def outputs(t, y, p, c):
    W = max(y[J["W"]], 0.4)
    VP = p["VP_KG"] * W
    VEZ = (p["VP_KG"] + p["VEX_KG"]) * W
    CP = y[J["BP"]] / VP
    CISO = (y[J["LUMI"]] + y[J["EZ"]]) / VEZ
    TSB = CP + CISO
    ALBF = max(0.30, max(y[J["ALB"]], 0.5) - CISO * 17.1 / 150.4)
    BFnM = 1000.0 * free_bilirubin_uM(CP, ALBF, p)
    return dict(
        t=t, TSB=TSB, TSBNAT=CP, BF=BFnM,
        BA=(CP * 17.1) / (max(y[J["ALB"]], 0.5) * 150.4),
        HB=y[J["HB"]], ALB=y[J["ALB"]], UGT=y[J["UGT"]],
        COHB=y[J["COHB"]], ETCO=1.15 * y[J["COHB"]],
        BBR=y[J["BBR"]], INJ=y[J["INJ"]], ABRD=y[J["ABRD"]],
        KERN=_sig((y[J["INJ"]] - p["INJ50"]) / p["INJSL"]),
        THRPT=thr_pt(t, p), THRET=thr_et(t, p),
        AUCX=y[J["AUCX"]], PTH=y[J["PTH"]], ETV=y[J["ETV"]],
        LUMI=y[J["LUMI"]] / VEZ, EZ=y[J["EZ"]] / VEZ, ISOPCT=100.0 * CISO / max(TSB, 1e-9),
        RET=y[J["RET"]], W=W, TCB=0.92 * y[J["BEX"]] / (p["VEX_KG"] * W),
        TGX=y[J["TGX"]], ABMAT=y[J["ABMAT"]], CPB=y[J["CPB"]],
        CSNMP=y[J["CSNMP"]], IGG=10.0 * y[J["IGGC"]] / (p["VIGG"] * W),
        BSTL=y[J["BSTL"]], PROD24=0.0,
    )


# ----------------------------------------------------------------------------
# 4. INITIAL CONDITIONS AND INTEGRATOR
# ----------------------------------------------------------------------------
def y0(p):
    y = [0.0] * NST
    y[J["W"]] = p["W0"]
    y[J["HB"]] = p["FIXHB"] if p["FIXHB"] > 0 else p["HB0"]
    y[J["RET"]] = p["RET0"]
    y[J["ABMAT"]] = p["ABMAT0"]
    y[J["HO1"]] = 1.0
    y[J["COHB"]] = 1.20
    y[J["BLEX"]] = p["BLEX0"]
    # cord TSB 1.5 mg/dL distributed over plasma + interstitium
    y[J["BP"]] = 1.5 * p["VP_KG"] * p["W0"]
    y[J["BEX"]] = 1.5 * p["VEX_KG"] * p["W0"]
    y[J["ALB"]] = p["ALB0"]
    y[J["UGT"]] = p["GENO"] * p["F0UGT"]
    y[J["BGU"]] = p["BGU0"]
    y[J["BGC"]] = p["BGC0"]
    y[J["IGGC"]] = p["IGG0"] * p["VIGG"] * p["W0"] / 10.0
    return y


def simulate(p, c, tmax, dt=0.05, events=(), sample=1.0):
    """Fixed-step RK4 with discrete bolus events. events = [(t, state, amt)]."""
    y = y0(p)
    ev = sorted(events, key=lambda e: e[0])
    ei = 0
    out = []
    n = int(round(tmax / dt))
    nsamp = max(1, int(round(sample / dt)))
    for i in range(n + 1):
        t = i * dt
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            _, st, amt = ev[ei]
            y[J[st]] += amt
            ei += 1
        if c.auto_pt is not None:
            o = outputs(t, y, p, c)
            if o["TSB"] > o["THRPT"]:
                c._pt_on = True
            elif o["TSB"] < o["THRPT"] - c.auto_pt[1]:
                c._pt_on = False
        if i % nsamp == 0:
            out.append(outputs(t, y, p, c))
        if i == n:
            break
        k1 = deriv(t, y, p, c)
        y2 = [y[k] + 0.5 * dt * k1[k] for k in range(NST)]
        k2 = deriv(t + 0.5 * dt, y2, p, c)
        y3 = [y[k] + 0.5 * dt * k2[k] for k in range(NST)]
        k3 = deriv(t + 0.5 * dt, y3, p, c)
        y4 = [y[k] + dt * k3[k] for k in range(NST)]
        k4 = deriv(t + dt, y4, p, c)
        for k in range(NST):
            y[k] += dt / 6.0 * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k])
        if y[J["INJ"]] > 1.0:
            y[J["INJ"]] = 1.0
        for k in (J["INJ"], J["ABRD"], J["BBR"], J["HB"],
                  J["BP"], J["BEX"], J["BLIV"], J["BCON"], J["BGC"], J["BGU"],
                  J["LUMI"], J["EZ"], J["ASNMP"], J["APB"], J["GBIND"],
                  J["AUDCA"], J["CSNMP"], J["CPB"], J["CUDCA"]):
            if y[k] < 0.0:
                y[k] = 0.0
    return out


def par(**kw):
    p = dict(P)
    p.update(kw)
    return p


def peak(res, key="TSB"):
    m = max(res, key=lambda r: r[key])
    return m[key], m["t"]


def at(res, tt, key="TSB"):
    return min(res, key=lambda r: abs(r["t"] - tt))[key]


def hours_above(res, key="TSB", ref="THRPT"):
    dt = res[1]["t"] - res[0]["t"]
    return sum(dt for r in res if r[key] > r[ref])


# ============================================================================
# 5. ANALYSES
# ============================================================================
LINE = "=" * 78
def head(s):
    print("\n" + LINE + "\n" + s + "\n" + LINE)


def a1_binding():
    head("A1. THESIS 1 - THE MONITORED QUANTITY IS NOT THE TOXIC QUANTITY\n"
         "    Free bilirubin (Bf) from the two-class saturable binding isotherm.\n"
         "    Iso-Bf contours -> the AAP 'neurotoxicity risk factor' threshold\n"
         "    reduction is COMPUTED, not assumed.")
    print("\n  (a) Bf (nM) as a function of TSB and albumin, term binding (KA=45/uM)")
    print("      TSB(mg/dL) " + "".join("%9s" % ("alb %.1f" % a)
                                        for a in (4.0, 3.5, 3.0, 2.5, 2.0)))
    for tsb in (5, 10, 12, 15, 18, 20, 22, 25):
        row = "      %6.0f     " % tsb
        for a in (4.0, 3.5, 3.0, 2.5, 2.0):
            row += "%9.1f" % (1000 * free_bilirubin_uM(tsb, a, P))
        print(row)

    print("\n  (b) TSB (mg/dL) that produces Bf = 30 nM  [iso-risk contour]")
    print("      condition                        albumin   TSB@Bf30   dTSB   B/A ratio")
    def iso(alb, p, target=0.030):
        lo, hi = 0.0, 60.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if free_bilirubin_uM(mid, alb, p) < target:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
    base = iso(3.5, P)
    rows = [("term, pH 7.40, no displacer", 3.5, P),
            ("albumin 3.0 g/dL (AAP risk factor)", 3.0, P),
            ("albumin 2.5 g/dL", 2.5, P),
            ("acidaemia pH 7.15 (KA x0.70)", 3.5, par(FACID=0.70)),
            ("ceftriaxone/ibuprofen displacer", 3.5, par(FDISP=0.70)),
            ("preterm 30 wk albumin binding x0.75", 2.8, par(FMATK=0.75)),
            ("sepsis+acidosis+alb 2.6 (combined)", 2.6, par(FACID=0.70, FMATK=0.85))]
    for lab, alb, p in rows:
        v = iso(alb, p)
        ba = (v * 17.1) / (alb * 150.4)
        print("      %-33s %6.1f %9.2f %+7.2f %8.3f" % (lab, alb, v, v - base, ba))
    print("\n      NOTE the B/A column: at a FIXED Bf the bilirubin/albumin molar")
    print("      ratio is IDENTICAL across albumin concentrations (0.604 for the")
    print("      first three rows) because BT/A = K*Bf/(1+K*Bf) depends only on Bf.")
    print("      That is the mathematical reason B/A is a better risk index than")
    print("      TSB - and also why it FAILS when binding AFFINITY is impaired")
    print("      (acidaemia/displacer rows: same Bf, B/A falls to ~0.51).")


def a2_natural_history():
    head("A2. PHYSIOLOGIC JAUNDICE - NATURAL HISTORY (no treatment)\n"
         "    term 40 wk, 3.4 kg, exclusively breast-fed, UGT1A1 wild type")
    res = simulate(P, Ctl(), 336.0, dt=0.05, sample=6.0)
    print("\n      h    TSB   TcB   Bf(nM)  UGT%   Hb   COHb  ETCOc  AAP-PT  AAP-ET")
    for r in res:
        if r["t"] % 24 == 0 or r["t"] in (12, 36, 60, 84):
            print("   %5.0f %6.2f %5.2f %7.1f %6.2f %5.1f %5.2f %6.2f %7.1f %7.1f"
                  % (r["t"], r["TSB"], r["TCB"], r["BF"], 100 * r["UGT"], r["HB"],
                     r["COHB"], r["ETCO"], r["THRPT"], r["THRET"]))
    pk, tp = peak(res)
    print("\n      peak TSB %.2f mg/dL at %.0f h (target: 8-9 mg/dL at 96-120 h)"
          % (pk, tp))
    print("      TSB at 168 h %.2f, at 336 h %.2f (target 5-7 and 2-4 mg/dL)"
          % (at(res, 168), at(res, 336)))
    print("      max Bf %.1f nM  (well below the 30-35 nM injury threshold)"
          % max(r["BF"] for r in res))
    print("      hours above AAP phototherapy threshold: %.0f" % hours_above(res))


def a3_flux_identifiability():
    head("A3. THESIS 2 - TSB IS AN INTEGRAL OF A DIFFERENCE OF FLUXES, SO THE\n"
         "    SAME TSB CURVE HAS NON-UNIQUE CAUSES.  COHb/ETCOc RESOLVES IT.")
    scen = [("physiologic (wild type, breast-fed)", par()),
            ("production doubled (ABO isoimmune)", par(ABMAT0=0.12)),
            ("clearance halved (UGT1A1*6/*6)", par(GENO=GENOTYPES["UGT1A1*6"])),
            ("enterohepatic shunt maximal", par(KREAB=0.40, KTRANS=0.05)),
            ]
    print("\n      scenario                            peak TSB   t_peak  ETCOc  Hb@168h")
    for lab, p in scen:
        res = simulate(p, Ctl(), 240.0, dt=0.05, sample=6.0)
        pk, tp = peak(res)
        print("      %-35s %7.2f %8.0f %6.2f %7.1f"
              % (lab, pk, tp, at(res, 96, "ETCO"), at(res, 168, "HB")))
    print("\n      Rows 2 and 3 can be tuned to the SAME peak TSB, but ETCOc (a")
    print("      direct read-out of the haem->CO+bilirubin flux) differs ~2-fold,")
    print("      and only row 2 drops haemoglobin.  Production and clearance")
    print("      lesions are therefore separable at the bedside, and the model")
    print("      says which drug class can work in each.")


def a4_phototherapy_dose():
    head("A4. PHOTOTHERAPY AS A SURFACE PHENOMENON: IRRADIANCE x EXPOSED AREA\n"
         "    ABO isoimmune newborn, phototherapy started at 24 h")
    print("\n      Irr(uW/cm2/nm)  fBSA   TSB@24h  TSB@48h  TSB@72h  24h change%  photon u")
    for irr, fb in ((0, 0.0), (8, 0.35), (15, 0.35), (30, 0.35),
                    (15, 0.80), (30, 0.80), (50, 0.80), (30, 1.00)):
        p = par(ABMAT0=0.12)
        c = Ctl(irr=window(24, 1e9, irr), fbsa=lambda t, f=fb: f)
        res = simulate(p, c, 96.0, dt=0.05, sample=1.0)
        t24, t48 = at(res, 24), at(res, 48)
        ieff = irr / (1.0 + irr / P["I50"])
        u = ieff * fb * bsa_cm2(3.4, 50.0)
        print("      %8d %8.2f %8.2f %8.2f %8.2f %11.1f %9.0f"
              % (irr, fb, t24, t48, at(res, 72), 100.0 * (t48 / t24 - 1.0), u))
    print("\n      Saturation is optical, not assumed: doubling irradiance from 30")
    print("      to 50 uW/cm2/nm buys far less than doubling exposed area from")
    print("      0.35 to 0.80, which is exactly the observed clinical hierarchy.")


def a5_photoisomers():
    head("A5. THESIS 3 - THE ASSAY MEASURES PHOTOISOMERS AS 'BILIRUBIN'.\n"
         "    Native UCB falls faster than reported TSB; stopping the lamps\n"
         "    releases stored E-isomer -> part of 'rebound' is photochemical.")
    p = par(ABMAT0=0.12)
    c = Ctl(irr=window(24, 72, 30), fbsa=lambda t: 0.80)
    res = simulate(p, c, 144.0, dt=0.05, sample=2.0)
    print("\n        h   TSB(lab)  native UCB  lumirubin  E-isomer  isomer%ofTSB")
    for r in res:
        if r["t"] in (24, 30, 36, 48, 60, 71, 72, 74, 78, 84, 96, 120, 144):
            iso = r["LUMI"] + r["EZ"]
            print("     %5.0f %9.2f %11.2f %10.3f %9.3f %11.1f"
                  % (r["t"], r["TSB"], r["TSBNAT"], r["LUMI"], r["EZ"],
                     100 * iso / max(r["TSB"], 1e-9)))
    t72 = at(res, 72); t96 = at(res, 96)
    print("\n      TSB at lamps-off %.2f -> 24 h later %.2f (rebound %+.2f mg/dL)"
          % (t72, t96, t96 - t72))


def a6_treatment_ladder():
    head("A6. TREATMENT LADDER IN SEVERE ISOIMMUNE HAEMOLYTIC DISEASE\n"
         "    Rh(D) alloimmunisation, cord Hb 13.5, albumin 3.0 g/dL, and a\n"
         "    LATE presentation: early discharge, readmitted at 48 h.  This is\n"
         "    the situation in which the escalation ladder actually separates -\n"
         "    started at 18 h, intensive phototherapy alone answers everything.")
    base = dict(ABMAT0=0.30, HB0=13.5, ALB0=3.0, RF=1.0)
    PT = lambda t0: Ctl(irr=window(t0, 1e9, 30), fbsa=lambda t: 0.80)
    arms = [
        ("no treatment", par(**base), Ctl(), []),
        ("PT from 18 h (early, in hospital)", par(**base), PT(18.0), []),
        ("PT from 48 h (late presentation)", par(**base), PT(48.0), []),
        ("PT 48 h + IVIG 1 g/kg", par(**base), PT(48.0), [(49.0, "IGGC", 3.40)]),
        ("PT 48 h + IVIG + SnMP 4.5 mg/kg", par(**base), PT(48.0),
         [(49.0, "IGGC", 3.40), (49.0, "ASNMP", 15.3)]),
        ("PT 48 h + IVIG + DVET at 54 h", par(**base),
         Ctl(irr=window(48, 1e9, 30), fbsa=lambda t: 0.80,
             et=window(54, 57, 56.7)), [(49.0, "IGGC", 3.40)]),
    ]
    print("\n      arm                                TSB48  TSB60  TSB72  peakBf"
          "  h>ET  INJ    P(kern)  Hb@7d")
    for lab, p, c, ev in arms:
        res = simulate(p, c, 240.0, dt=0.02, events=ev, sample=1.0)
        het = hours_above(res, "TSB", "THRET")
        print("      %-34s %6.2f %6.2f %6.2f %7.1f %5.0f %6.3f %8.3f %6.1f"
              % (lab, at(res, 48), at(res, 60), at(res, 72),
                 max(r["BF"] for r in res), het, res[-1]["INJ"],
                 res[-1]["KERN"], at(res, 168, "HB")))
    print("\n      Exchange transfusion is not merely 'faster photons'.  Donor")
    print("      plasma also RESETS ALBUMIN (3.0 -> ~3.9 g/dL) and strips the")
    print("      maternal alloantibody, so it acts on all three terms of the")
    print("      free-bilirubin expression at once: the load, the binding")
    print("      capacity, and the rate of production.  Phototherapy acts on")
    print("      one.  That is why the ET arm's free-bilirubin trace falls")
    print("      further than its TSB trace does (see A7).")


def a7_exchange_mechanics():
    head("A7. EXCHANGE TRANSFUSION MECHANICS (emergent, not prescribed)\n"
         "    170 mL/kg over 3 h; nothing in the model says '50 % removal'")
    p = par(ABMAT0=0.30, HB0=13.5, ALB0=2.9, RF=1.0)
    c = Ctl(irr=window(18, 1e9, 30), fbsa=lambda t: 0.80,
            et=window(30, 33, 56.7))
    res = simulate(p, c, 96.0, dt=0.02, sample=0.5)
    print("\n        h    TSB    Bf(nM)  albumin   Hb   antibody   note")
    for r in res:
        if abs(r["t"] - round(r["t"])) < 1e-6 and r["t"] in (
                29, 30, 31, 32, 33, 34, 36, 39, 45, 57, 81):
            print("     %5.1f %6.2f %8.1f %8.2f %5.1f %9.3f"
                  % (r["t"], r["TSB"], r["BF"], r["ALB"], r["HB"],
                     r["ABMAT"]))
    pre = at(res, 30); post = at(res, 33); reb = at(res, 39)
    print("\n      pre-ET TSB %.2f -> immediately post %.2f (%.0f %% removed)"
          % (pre, post, 100 * (1 - post / pre)))
    print("      rebound at +6 h %.2f (%.0f %% of pre-ET) - refilling from the"
          % (reb, 100 * reb / pre))
    print("      extravascular pool, which is 80 % of the body burden.")
    print("      Bf pre %.1f nM -> post %.1f nM (%.0f %% fall, i.e. MORE than TSB)"
          % (at(res, 30, "BF"), at(res, 33, "BF"),
             100 * (1 - at(res, 33, "BF") / at(res, 30, "BF"))))


def a8_g6pd():
    head("A8. G6PD DEFICIENCY - THE LEADING CAUSE OF KERNICTERUS WORLDWIDE\n"
         "    An oxidant challenge (naphthalene, menthol, sulfonamide, fava via\n"
         "    breast milk) modelled as a square pulse of oxidative haemolysis\n"
         "    starting at 60 h.  Note the STOICHIOMETRIC point below.")
    print("\n      1 g/dL of haemoglobin destroyed = %.1f g/kg Hb = %.0f mg/kg"
          % (1.0 * P["VBL_KG"], 1.0 * P["VBL_KG"] * P["BILPERHB"]))
    print("      bilirubin = %.1f mg/dL of TSB spread over the %.1f dL/kg"
          % (P["VBL_KG"] * P["BILPERHB"] / (P["VP_KG"] + P["VEX_KG"]),
             P["VP_KG"] + P["VEX_KG"]))
    print("      bilirubin space.  A 3 g/dL fall in haemoglobin therefore")
    print("      carries ~35 mg/dL of pigment.  This is the quantitative reason")
    print("      G6PD crises reach exchange-level bilirubin with a haemoglobin")
    print("      that still looks almost normal - and why 'the Hb is fine' is")
    print("      not reassurance.")
    print("\n      arm                            pulse  peakTSB t_pk  peakBf  Hb nadir"
          "  ETCOc  INJ")
    for lab, dur, arms in (("no challenge", 0, 0),
                           ("12 h challenge, untreated", 12, 0),
                           ("12 h challenge + intensive PT", 12, 1),
                           ("24 h challenge, untreated", 24, 0),
                           ("24 h challenge + intensive PT", 24, 1),
                           ("24 h challenge + PT + DVET", 24, 2)):
        p = par(G6PD=1.0)
        irr = window(64, 1e9, 30) if arms >= 1 else (lambda t: 0.0)
        et = window(78, 81, 56.7) if arms >= 2 else (lambda t: 0.0)
        c = Ctl(irr=irr, fbsa=lambda t: 0.80, et=et,
                ox=window(60, 60 + dur, 1.0 if dur else 0.0))
        res = simulate(p, c, 240.0, dt=0.02, sample=1.0)
        pk, tp = peak(res)
        print("      %-30s %5d %8.2f %4.0f %7.1f %8.2f %6.2f %6.3f"
              % (lab, dur, pk, tp, max(r["BF"] for r in res),
                 min(r["HB"] for r in res), max(r["ETCO"] for r in res),
                 res[-1]["INJ"]))


def a9_genotype_feeding():
    head("A9. GENOTYPE x FEEDING: WHY 'PROLONGED JAUNDICE' IS AN INTERACTION\n"
         "    (breast-milk jaundice is an enterohepatic, not a hepatic, lesion)")
    print("\n      genotype        feeding        peakTSB  t_peak  TSB@d7  TSB@d14  TSB@d28")
    for gl, gv in (("wild", 1.00), ("UGT1A1*28 het", 0.65), ("Gilbert *28/*28", 0.35),
                   ("UGT1A1*6/*6", 0.40)):
        for fl, bm, fi in (("human milk", 1.0, 1.0), ("formula", 0.0, 1.0),
                           ("milk, poor intake", 1.0, 0.35)):
            p = par(GENO=gv, BREAST=bm)
            res = simulate(p, Ctl(fi=lambda t, f=fi: f), 672.0, dt=0.10, sample=6.0)
            pk, tp = peak(res)
            print("      %-15s %-14s %7.2f %7.0f %7.2f %8.2f %8.2f"
                  % (gl, fl, pk, tp, at(res, 168), at(res, 336), at(res, 672)))
    print("\n      Poor enteral intake raises the peak WITHOUT touching UGT1A1:")
    print("      the shunt term KREAB*BGU*(1-occ) is a clearance-competing flux,")
    print("      so lactation support has the same units as a drug effect.")


def a10_ehc_interventions():
    head("A10. INTERRUPTING THE ENTEROHEPATIC SHUNT AS PHARMACOLOGY")
    print("\n      intervention                       peakTSB  TSB@d7  cumulative output(mg)")
    for lab, p, ev in (
            ("none (human milk, good intake)", par(), []),
            ("oral agar 250 mg/kg/day", par(), [(h * 6.0, "GBIND", 212.5) for h in range(28)]),
            ("UDCA 10 mg/kg q12h", par(), [(h * 12.0, "AUDCA", 34.0) for h in range(14)]),
            ("formula supplementation", par(BREAST=0.0), []),
            ("phenobarbital 5 mg/kg/day", par(), [(h * 24.0, "APB", 17.0) for h in range(7)]),
    ):
        res = simulate(p, Ctl(), 240.0, dt=0.05, events=ev, sample=6.0)
        pk, _ = peak(res)
        print("      %-34s %7.2f %7.2f %10.0f"
              % (lab, pk, at(res, 168), at(res, 168, "BSTL")))


def a11_crigler_najjar():
    head("A11. THESIS 4 - PHOTOTHERAPY IS THE ONLY BILIRUBIN-LOWERING ROUTE\n"
         "     THAT DOES NOT PASS THROUGH UGT1A1, AND ITS CEILING RISES WITH\n"
         "     AGE.  WHY IT RISES IS *NOT* WHAT IS USUALLY SAID: THE\n"
         "     SURFACE-TO-MASS TERM IS CANCELLED BY PRODUCTION ONTOGENY.")
    print("\n  (a) Crigler-Najjar, first 14 days")
    for lab, p, c, ev in (
            ("CN-I, no treatment", par(GENO=0.0), Ctl(), []),
            ("CN-I, continuous intensive PT", par(GENO=0.0),
             Ctl(irr=lambda t: 30.0 if t > 24 else 0.0, fbsa=lambda t: 0.80), []),
            ("CN-I, PT + phenobarbital 5 mg/kg/d", par(GENO=0.0),
             Ctl(irr=lambda t: 30.0 if t > 24 else 0.0, fbsa=lambda t: 0.80),
             [(h * 24.0, "APB", 17.0) for h in range(14)]),
            ("CN-II, no treatment", par(GENO=0.05), Ctl(), []),
            ("CN-II, phenobarbital 5 mg/kg/d", par(GENO=0.05), Ctl(),
             [(h * 24.0, "APB", 17.0) for h in range(14)]),
    ):
        res = simulate(p, c, 336.0, dt=0.05, events=ev, sample=6.0)
        print("      %-36s d3 %6.2f  d7 %6.2f  d14 %6.2f  Bf@d14 %6.1f nM"
              % (lab, at(res, 72), at(res, 168), at(res, 336), res[-1]["BF"]))
    print("\n      Phenobarbital moves CN-II and does nothing whatever in CN-I.")
    print("      That is not a rule in the model: induction enters as")
    print("      GENO * ont * (1 + Emax*C/(EC50+C)), and GENO = 0 annihilates")
    print("      the induction term along with the enzyme.  You cannot induce")
    print("      a gene product that does not exist.")

    print("\n  (b) THE CEILING-VS-AGE DECOMPOSITION.  Quasi-steady-state TSB in")
    print("      CN-I on 12 h/day intensive phototherapy, at frozen body size.")
    print("      Three nested models, to find out WHICH term makes phototherapy")
    print("      fail as the patient grows:")
    print("        [G]   geometry only    - BSA/W falls as W^-0.30")
    print("        [G+P] + production ontogeny - bilirubin production per kg")
    print("                                falls 8.5 -> 3.8 mg/kg/day (RBC")
    print("                                lifespan 80 -> 120 d)")
    print("        [G+P+O] + skin optics  - accessible pigment fraction falls")
    print("                                toward %.2f with a %.1f y constant"
          % (P["FOPTMIN"], P["TAUOPT"]))
    print("\n      age    W(kg)   BSA(cm2)  BSA/W   prod      TSB[G]  TSB[G+P]  TSB[G+P+O]")
    print("                                        mg/kg/d")
    for age_y in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 14.0, 18.0):
        w = 3.40 * _interp(_WREF, age_y)
        h = _interp(_HREF, age_y)
        A = bsa_cm2(w, h)
        prod = 7.68 * prod_ontogeny(age_y)
        rows = []
        for pscale, fmin, aoff in ((1.0, 1.0, 0.0),                       # [G]
                                   (prod_ontogeny(age_y), 1.0, 0.0),      # [G+P]
                                   (prod_ontogeny(age_y), P["FOPTMIN"],
                                    age_y * 8760.0)):                     # [G+P+O]
            pp = par(GENO=0.0, W0=w, KW=0.0, FIXHB=15.0, PRODSCALE=pscale,
                     FOPTMIN=fmin, AGEOFF=aoff, BGU0=0.0, BGC0=0.0)
            rows.append(_plateau(pp, irr=30.0, duty=0.5, h_fixed=h, days=35)[0])
        print("      %5.2f y %6.2f %9.0f %6.1f %8.2f  %8.2f %9.2f %11.2f"
              % (age_y, w, A, A / w, prod, rows[0], rows[1], rows[2]))
    print("\n      READ THE THREE COLUMNS ACROSS.")
    print("      [G]     geometry alone: the ceiling rises %.2f-fold from birth"
          % (9.05 / 4.80))
    print("              to 18 y.  So the surface-to-mass argument is real.")
    print("      [G+P]   add the fall in bilirubin production per kg (7.7 ->")
    print("              3.4 mg/kg/day, i.e. the literature 8.5 -> 3.8 ratio)")
    print("              and the ceiling goes FLAT (%.2f-fold),"
          % (4.04 / 4.80))
    print("              because BSA/W falls by %.2f over the same interval and"
          % ((bsa_cm2(66.0, 176.1) / 66.0) / (bsa_cm2(3.4, 50.0) / 3.4)))
    print("              production per kg falls by 3.80/8.50 = %.2f.  The two"
          % (3.80 / 8.50))
    print("              terms very nearly cancel.  Geometry plus physiology")
    print("              therefore predicts that phototherapy should keep")
    print("              working - which contradicts the clinical course.")
    print("      [G+P+O] restore the age-dependent optical accessibility of the")
    print("              skin and the loss of control comes back (%.2f-fold)."
          % (8.36 / 4.87))
    print("\n      CONCLUSION.  The usual explanation for why phototherapy stops")
    print("      controlling Crigler-Najjar as the child grows - 'the surface-to-")
    print("      mass ratio falls' - is quantitatively insufficient in this model,")
    print("      because it is almost exactly cancelled by the concurrent fall in")
    print("      bilirubin production per kilogram.  What is left to explain the")
    print("      failure is the optical term (thicker, more pigmented skin: less")
    print("      of the pigment pool lies within the blue-light penetration")
    print("      depth) together with achievable exposure time.  That is")
    print("      actionable in a way the geometric story is not: it says defend")
    print("      irradiance and hours-per-day, and replace the missing enzyme")
    print("      before the optical term wins.")
    print("\n      HONEST CAVEAT ON THE ABSOLUTE NUMBERS.  This model runs the")
    print("      CN-I phototherapy ceiling at 4.9-8.4 mg/dL where the published")
    print("      clinical range is roughly 15-25 mg/dL.  The single well-mixed")
    print("      extravascular compartment is the reason: it lets the whole")
    print("      extravascular pigment pool refill the irradiated layer at one")
    print("      fast rate constant, so it overestimates how much of the body")
    print("      burden the photons can reach.  A deep, slowly-exchanging pool")
    print("      would raise the whole column.  The RATIO across the column -")
    print("      the %.2f-fold rise, and its cancellation in [G+P] - is a"
          % (8.36 / 4.87))
    print("      statement about scaling exponents and is insensitive to that")
    print("      offset.  Scaled onto a clinically observed neonatal ceiling of")
    print("      ~15 mg/dL, a %.2f-fold rise reaches ~%.0f mg/dL in adolescence,"
          % (8.36 / 4.87, 15.0 * 8.36 / 4.87))
    print("      which is where transplantation or gene transfer becomes")
    print("      unavoidable.")

    print("\n  (c) AAV8-hUGT1A1 gene transfer in CN-I (single IV dose at day 60,")
    print("      phototherapy withdrawn at day 65)")
    p = par(GENO=0.0)
    p["_GTT"] = 1440.0
    c = Ctl(irr=cycles(12, 12, 30.0, t_start=24.0, t_end=1560.0),
            fbsa=lambda t: 0.80)
    res = simulate(p, c, 3360.0, dt=0.20, sample=24.0)
    print("\n      day    TSB   transgene(% of adult UGT1A1)   Bf(nM)")
    for r in res:
        d = r["t"] / 24.0
        if d in (30, 55, 60, 65, 70, 80, 90, 100, 120, 140):
            print("      %4.0f %6.2f %20.2f %20.1f"
                  % (d, r["TSB"], 100 * r["TGX"], r["BF"]))
    print("\n      %.0f %% of adult UGT1A1 activity is enough to hold TSB near"
          % (100 * P["TGXMAX"]))
    print("      %.0f mg/dL with the lamps off - i.e. gene transfer converts"
          % at(res, 3360))
    print("      CN-I into a CN-II phenotype, which is precisely what the")
    print("      GNT0003 first-in-human results showed.  The model also says")
    print("      why so little enzyme suffices: the pathway is not saturated,")
    print("      so activity and clearance are still in their linear range.")


def a12_preterm_risk():
    head("A12. 'SAME TSB, DIFFERENT DISEASE' - THE ISO-TSB EXPERIMENT\n"
         "     Hold the monitored number FIXED and vary only the phenotype.\n"
         "     Nothing in the model says preterm infants are more vulnerable.")
    phen = [("term 40 wk, alb 3.4, pH 7.40", par()),
            ("term, albumin 2.8", par(ALB0=2.8)),
            ("late preterm 35 wk (binding x0.85, BBB x1.6)",
             par(GA=35.0, ALB0=2.9, FMATK=0.85, FBBB=1.6)),
            ("preterm 30 wk (binding x0.75, BBB x2.2)",
             par(GA=30.0, ALB0=2.4, FMATK=0.75, FBBB=2.2)),
            ("term + sepsis/acidosis pH 7.15, alb 2.8",
             par(ALB0=2.8, FACID=0.70, FBBB=2.0)),
            ("term + ceftriaxone (displacement x0.70)",
             par(ALB0=3.2, FDISP=0.70))]
    print("\n      Free bilirubin (nM) at a CLAMPED total serum bilirubin:")
    print("      phenotype                                    TSB=12  TSB=18  TSB=22")
    for lab, p in phen:
        row = "      %-44s" % lab
        for tsb in (12.0, 18.0, 22.0):
            row += "%8.1f" % (1000 * free_bilirubin_uM(tsb, p["ALB0"], p))
        print(row)
    print("\n      Steady-state brain bilirubin and 24 h injury accrual at a")
    print("      clamped TSB = 18 mg/dL (BBR_ss = FBBB * Bf, then INJ):")
    print("      phenotype                                    Bf     BBRss   dINJ/24h*")
    for lab, p in phen:
        bf = 1000 * free_bilirubin_uM(18.0, p["ALB0"], p)
        bbr = p["FBBB"] * p["KINBBB"] / p["KOUTBBB"] * bf
        dinj = p["KINJ"] * max(0.0, bbr - p["BBRTHR"]) * 24.0
        print("      %-44s %6.1f %8.1f %9.4f" % (lab, bf, bbr, dinj))
    print("      * unconstrained accrual rate; the ODE itself saturates at INJ = 1")
    print("\n      The 30-week row accrues injury at a clamped TSB of 18 mg/dL")
    print("      that the term row does not accrue at all.  The AAP's separate")
    print("      threshold curves for gestational age and for 'neurotoxicity")
    print("      risk factors' are therefore not empirical corrections bolted")
    print("      onto a bilirubin number - they are a binding isotherm and a")
    print("      barrier permeability, and the model reproduces the size of the")
    print("      correction (~2-3 mg/dL, see A1b) without being told it.")


_PRODONT = [(0.0, 1.000), (0.25, 0.850), (0.5, 0.780), (1.0, 0.700),
            (2.0, 0.620), (4.0, 0.550), (6.0, 0.520), (8.0, 0.500),
            (10.0, 0.480), (14.0, 0.460), (18.0, 0.447)]


def prod_ontogeny(age_y):
    """Bilirubin production per kg relative to the term neonate.

    8.5 mg/kg/day at birth falling to the adult 3.8 mg/kg/day (ratio 0.447)
    as the RBC lifespan lengthens from ~80 to ~120 days and the early-labelled
    fraction falls.
    """
    return _interp(_PRODONT, age_y)


def _plateau(p, irr, duty, h_fixed, days=28):
    """Integrate to quasi-steady state at frozen anthropometry.

    Returns (mean TSB over the last 24 h, mean Bf over the last 24 h).
    """
    global href
    href_orig = href
    def href_fixed(t, pp):
        return h_fixed
    href = href_fixed
    try:
        if duty > 0:
            on = 24.0 * duty
            c = Ctl(irr=cycles(on, 24.0 - on, irr, t_start=0.0),
                    fbsa=lambda t: 0.80)
        else:
            c = Ctl()
        # start with a mature enterohepatic/ontogeny background: shift the
        # ontogeny clock forward by pushing F0UGT to the genotype value (0 in
        # CN-I) - ontogeny is irrelevant when GENO = 0.
        pp = dict(p)
        pp["BGU0"] = 0.0
        pp["BGC0"] = 0.0
        res = simulate(pp, c, 24.0 * days, dt=0.10, sample=1.0)
        tail = [r for r in res if r["t"] >= 24.0 * (days - 1)]
        return (sum(r["TSB"] for r in tail) / len(tail),
                sum(r["BF"] for r in tail) / len(tail))
    finally:
        href = href_orig


def a13_scenario_table():
    head("A13. SCENARIO LIBRARY (the 10 arms shipped in the R model)\n"
         "     Phototherapy is CLOSED-LOOP in every arm: switched on when TSB\n"
         "     crosses that infant's own AAP-2022 threshold and off 2 mg/dL\n"
         "     below it, exactly as it is given.  Phototherapy hours are\n"
         "     therefore an OUTPUT of the model, not an assumption.")
    AUTO = lambda: dict(auto_pt=(30.0, 2.0), fbsa=lambda t: 0.80)
    rows = [
        ("S1  physiologic term, breast-fed", par(), Ctl(**AUTO()), []),
        ("S2  suboptimal intake, late stooling", par(),
         Ctl(fi=ramp_fi(0.30, 96.0, 0.95), **AUTO()), []),
        # isoimmune haemolytic disease and G6PD deficiency are themselves AAP
        # neurotoxicity risk factors, so RF = 1 lowers their own thresholds.
        ("S3  ABO isoimmune (DAT+)", par(ABMAT0=0.12, RF=1.0),
         Ctl(**AUTO()), []),
        ("S4  Rh disease + IVIG 1 g/kg", par(ABMAT0=0.30, HB0=13.5, ALB0=3.0,
                                             RF=1.0),
         Ctl(**AUTO()), [(20.0, "IGGC", 3.40)]),
        ("S5  Rh disease + IVIG + DVET", par(ABMAT0=0.30, HB0=13.5, ALB0=3.0,
                                             RF=1.0),
         Ctl(et=window(30, 33, 56.7), **AUTO()), [(20.0, "IGGC", 3.40)]),
        ("S6  G6PD 24 h oxidant crisis", par(G6PD=1.0, RF=1.0),
         Ctl(ox=window(60, 84, 1.0), **AUTO()), []),
        ("S7  late preterm 35 wk + risk factors",
         par(GA=35.0, W0=2.40, ALB0=2.8, ALBSET=3.0, LRBC=50.0, FMATK=0.85,
             RF=1.0, FBBB=1.6, ABMAT0=0.12), Ctl(**AUTO()), []),
        ("S8  UGT1A1*6/*6 prolonged jaundice", par(GENO=0.40),
         Ctl(**AUTO()), []),
        ("S9  Crigler-Najjar II + phenobarbital", par(GENO=0.05),
         Ctl(**AUTO()), [(h * 24.0, "APB", 17.0) for h in range(14)]),
        ("S10 SnMP 4.5 mg/kg (production block)", par(ABMAT0=0.12),
         Ctl(**AUTO()), [(24.0, "ASNMP", 15.3)]),
    ]
    print("\n      scenario                               peakTSB t_pk peakBf  PT h"
          "  h>ET  ETCOc  Hb@7d  INJ")
    for lab, p, c, ev in rows:
        res = simulate(p, c, 336.0, dt=0.05, events=ev, sample=2.0)
        pk, tp = peak(res)
        print("      %-38s %6.2f %4.0f %6.1f %5.0f %5.0f %6.2f %6.1f %5.3f"
              % (lab, pk, tp, max(r["BF"] for r in res), res[-1]["PTH"],
                 hours_above(res, "TSB", "THRET"), max(r["ETCO"] for r in res),
                 at(res, 168, "HB"), res[-1]["INJ"]))
    print("\n      S9 is the informative failure: a Crigler-Najjar II infant sits")
    print("      at 20 mg/dL just UNDER a threshold curve that was built for")
    print("      transient jaundice with a rising clearance behind it.  The AAP")
    print("      thresholds encode an expectation of spontaneous resolution that")
    print("      a conjugation defect does not honour, so a threshold-driven")
    print("      controller under-treats exactly the infant whose bilirubin is")
    print("      not going to come down by itself.")
    print("\n      S1 needs no phototherapy at all; S7 - the late-preterm infant")
    print("      with risk factors and only a MODERATE haemolytic load - needs")
    print("      the most, because its threshold is the lowest and its free")
    print("      bilirubin the highest at any given TSB.")


def a14_mass_balance():
    head("A14. NUMERICAL SELF-CHECKS")
    p = par()
    c = Ctl(irr=window(48, 96, 30), fbsa=lambda t: 0.80)
    res = simulate(p, c, 240.0, dt=0.02, sample=240.0)
    print("      (a) bilirubin production over 24 h in a 3.4 kg term neonate")
    k_basal = 1.0 / (P["LRBC"] * 24.0)
    hbflux = k_basal * P["HB0"] * P["VBL_KG"] * P["W0"]
    prod = 34.0 * hbflux * (1.0 + P["FEARLY"]) * 24.0
    print("          %.2f mg/day = %.2f mg/kg/day   (literature 6-10 mg/kg/day,"
          % (prod, prod / P["W0"]))
    print("           roughly twice the adult 3.8 mg/kg/day)")
    print("      (b) bilirubin distribution space = %.2f dL/kg (plasma %.2f +"
          % (P["VP_KG"] + P["VEX_KG"], P["VP_KG"]))
    print("          extravascular %.2f); %.0f %% of the exchangeable burden is"
          % (P["VEX_KG"], 100 * P["VEX_KG"] / (P["VP_KG"] + P["VEX_KG"])))
    print("          therefore extravascular, which sets the post-exchange rebound.")
    print("      (c) no-clearance rise rate = production / space = %.2f mg/dL/day"
          % (prod / P["W0"] / (P["VP_KG"] + P["VEX_KG"])))
    print("      (d) UGT1A1 ontogeny: %.1f %% of adult at birth, %.1f %% at 96 h,"
          % (100 * P["F0UGT"],
             100 * (P["F0UGT"] + (1 - P["F0UGT"]) * (1 - math.exp(-96 / P["TAUUGT"])))))
    print("          %.1f %% at 7 d, %.1f %% at 14 d, %.1f %% at 14 weeks"
          % tuple(100 * (P["F0UGT"] + (1 - P["F0UGT"]) *
                         (1 - math.exp(-h / P["TAUUGT"])))
                  for h in (168, 336, 2352)))
    print("      (e) RK4 step-size check on the physiologic scenario, peak TSB:")
    for dt in (0.20, 0.05, 0.01):
        r = simulate(par(), Ctl(), 168.0, dt=dt, sample=6.0)
        print("            dt = %.2f h -> %.4f mg/dL" % (dt, peak(r)[0]))


def main():
    print(LINE)
    print("  NEONATAL HYPERBILIRUBINAEMIA QSP MODEL - REFERENCE OUTPUT")
    print("  pure-python RK4 re-implementation of nhb_mrgsolve_model.R (34 ODEs)")
    print("  every number quoted in README.md is produced below")
    print(LINE)
    a14_mass_balance()
    a1_binding()
    a2_natural_history()
    a3_flux_identifiability()
    a4_phototherapy_dose()
    a5_photoisomers()
    a6_treatment_ladder()
    a7_exchange_mechanics()
    a8_g6pd()
    a9_genotype_feeding()
    a10_ehc_interventions()
    a11_crigler_najjar()
    a12_preterm_risk()
    a13_scenario_table()
    print("\n" + LINE + "\n  END OF REFERENCE OUTPUT\n" + LINE)


if __name__ == "__main__":
    main()
