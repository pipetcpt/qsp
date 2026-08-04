#!/usr/bin/env python3
# =====================================================================
#  Radiation-Induced Lung Injury (RILI) -- independent reference
#  implementation of the mrgsolve model in rili_mrgsolve_model.R
#
#  WHY THIS FILE EXISTS
#  --------------------
#  No R runtime is available in the build environment, so every equation
#  in the mrgsolve model was first written and executed here, in
#  dependency-free Python with a fixed-step RK4 integrator.  The R file
#  is a line-by-line port of what this file computes.  Numbers quoted in
#  README.md come from `python3 rili_reference_model.py`.
#
#  The defects this re-implementation exposed are documented at the point
#  of the fix, each tagged  # [DEFECT n].
#
#  STRUCTURE
#  ---------
#  The lung is NOT one compartment at a mean dose.  It is discretised
#  into NB=6 dose bins taken from a dose-volume histogram (DVH): bin b
#  holds fractional lung volume v[b] and receives total physical dose
#  D[b].  All 9 biological states exist per bin; organ-level endpoints
#  (DLCO, FVC, CTCAE grade, fibrosis score) are volume-weighted
#  integrals over bins.  A treatment plan is therefore a *vector*, and
#  two plans with identical mean lung dose can differ in every other
#  moment of the same histogram.
#
#      per-bin states (10): AT2 DOOM EC SURF PERM CYT TGFB MFB COL XCOL
#      global states:       CYTS, BEDL[6], BEDT, TUMLN,
#                           DEXd DEXc PIRd PIRc NINd NINc DURc DURp ACEd ACEc
#      total ODEs:          60 + 1 + 6 + 1 + 1 + 10 = 79
# =====================================================================

import math
import sys

# ---------------------------------------------------------------------
# 1.  DOSE-VOLUME GEOMETRY
# ---------------------------------------------------------------------
# Bin edges expressed as a FRACTION of the prescription dose, so the same
# histogram shape can be re-scaled to any prescription.  Representative
# dose of a bin is the mid-point of its edges.
NB = 6
EDGE_F = [0.00, 0.08, 0.24, 0.44, 0.66, 0.88, 1.05]
REP_F = [(EDGE_F[i] + EDGE_F[i + 1]) / 2.0 for i in range(NB)]
#        -> [0.04, 0.16, 0.34, 0.55, 0.77, 0.965]


def bin_doses(DRX):
    """Representative total physical dose (Gy) of each bin."""
    return [DRX * f for f in REP_F]


def bin_edges(DRX):
    return [DRX * f for f in EDGE_F]


def mean_lung_dose(v, DRX):
    D = bin_doses(DRX)
    return sum(v[b] * D[b] for b in range(NB))


def vx(v, DRX, x):
    """Vx = fractional lung volume receiving >= x Gy, with linear
    interpolation inside the bin that contains x."""
    e = bin_edges(DRX)
    tot = 0.0
    for b in range(NB):
        lo, hi = e[b], e[b + 1]
        if hi <= x:
            continue
        if lo >= x:
            tot += v[b]
        else:
            tot += v[b] * (hi - x) / (hi - lo)
    return tot


def gEUD(v, DRX, a):
    """Generalised equivalent uniform dose.  a = 1/n of the
    Lyman-Kutcher-Burman volume-effect exponent.  a -> 1 is pure mean
    dose (fully parallel organ); a -> inf is max dose (fully serial)."""
    D = bin_doses(DRX)
    s = sum(v[b] * (D[b] ** a) for b in range(NB))
    return s ** (1.0 / a)


# ---------------------------------------------------------------------
# 2.  PARAMETERS
# ---------------------------------------------------------------------
P = dict(
    # ---- radiobiology -------------------------------------------------
    ABL=3.0,        # alpha/beta of late-responding lung (Gy)
    ABT=10.0,       # alpha/beta of NSCLC (Gy)
    ALPHAT=0.31,    # tumour alpha (1/Gy)
    FXWK=5.0,       # fractions per week

    # ---- AT2 (type II pneumocyte) compartment ------------------------
    KDAM=0.0250,    # AT2 lethal hits per Gy3 of lung BED
    KMIT=0.0250,    # 1/40 d : mitosis-linked death of lethally hit
                    # cells.  Adult AT2 turnover is slow, and this rate
                    # constant -- not the dose -- is what sets the latency.
    KREP=0.0300,    # logistic AT2 repopulation rate (1/d)
    # [DEFECT 3] KDCYT was 0.050 with KMCYT 0.50, i.e. a bystander-kill
    # slope of 0.10.  Multiplied by KDAMP that put 0.30 of self-gain into
    # the cytokine equation on its own, so ANY perturbation in ANY bin
    # ignited the whole lung and pneumonitis never resolved -- every
    # scenario, SBRT included, ended grade 4.  The slope is now 0.015,
    # which keeps the acute loop sub-critical (see the gain budget below).
    KDCYT=0.0150,   # cytokine-mediated bystander AT2 loss (1/d)
    KMCYT=1.0000,   # half-max cytokine for bystander loss

    # ---- microvascular endothelium -----------------------------------
    KECDAM=0.0300,  # endothelium is ~1.2x more radiosensitive than AT2
    KECREP=0.0120,  # endothelial repair in irradiated lung is slow
    KECCYT=0.0300,
    KECIRR=0.0043,  # permanent loss of microvascular repair ceiling per
                    # Gy3 of cumulative BED (capillary rarefaction)

    # ---- surfactant ---------------------------------------------------
    KSYN=0.5000,
    KSDEG=0.5000,

    # ---- alveolar-capillary permeability / oedema --------------------
    # [DEFECT 4] KPRM 0.25 / KPCYT 0.30 / KPCL 0.10 put the oedema state
    # at a steady value of ~7 in an irradiated bin.  PERM/(1+PERM) then
    # sat at 0.88 in every bin above ~20 Gy, so the pneumonitis index
    # lost all dynamic range and the volume-effect question could not be
    # asked.  Rescaled so a severely injured bin reaches PERM ~2.
    # [DEFECT 10] one shared coefficient KPRM multiplied both leak
    # sources with WSURF re-weighting the epithelial one.  Endothelial
    # integrity recovers on KECREP while the epithelial pool keeps
    # falling until day ~84, and with the two arms weighted alike they
    # CANCELLED: oedema peaked 12 d after the last fraction and no
    # setting of KMIT moved it (13 d -> 8 d over a 3.3x change), because
    # KMIT was not the controlling constant at all.  Radiation
    # pneumonitis peaks 4-12 weeks after completion, so the two sources
    # are now separate coefficients and the epithelial arm carries it.
    KPRMEC=0.0250,  # leak from endothelial loss
    KPRMEP=0.0900,  # leak from surfactant/epithelial loss (dominant)
                    # KPRMEP 0.090 + KPCL 0.045 -> peak at 5.7 wk post-RT
    KPCYT=0.1000,
    KPCL=0.0450,    # lymphatic clearance, tau 22 d (impaired post-RT)

    # ---- fast loop: DAMP -> NF-kB -> cytokine ------------------------
    # GAIN BUDGET.  The acute loop must be sub-critical or pneumonitis
    # never resolves:
    #   GIMM*(KAMP + KDAMP*KDCYT/KMCYT) + KSYSIN*KSOUT/KCLS  <  KCE
    #   no drug     1.00*(0.025+0.030) + 0.0167 = 0.0717 < 0.15
    #   durvalumab  1.35*(0.025+0.030) + 0.0167 = 0.0910 < 0.15
    # Checkpoint blockade therefore lengthens the effective cytokine
    # time-constant from 12.8 d to 17.0 d rather than adding a term.
    KDAMP=2.0000,   # DAMP yield per unit dying-cell flux
    KDIR=0.0050,    # direct NF-kB activation per Gy3/d during delivery
    KAMP=0.0250,    # autocrine cytokine amplification
    CMAX=2.0000,    # saturation of the amplification term
    KCE=0.1500,     # cytokine elimination (1/d)
    KSYSIN=0.1000,  # systemic -> local cytokine spill-in
    KSOUT=0.0500,   # local -> systemic spill-out
    KCLS=0.3000,    # systemic cytokine clearance

    # ---- slow loop: TGF-beta1 ----------------------------------------
    # [DEFECT 1] With KTGFA 0.20 / KTGFEL 0.20 and a Hill-2 myofibroblast
    # response of half-max 0.60, the loop gain at the HEALTHY fixed point
    # was 2.74.  The unirradiated lung was therefore not a fixed point at
    # all: collagen climbed to the fibrotic attractor with the beam off,
    # which is why FIB730 came out at 3.02 for every plan from SBRT to
    # 74 Gy, matched to three decimals.  Gain at the healthy state is now
    # 0.43, rising to 10.8 near TGFB 0.5 -- stable AND bistable.
    # [DEFECT 2] the mechanotransduction limb was written km*CX/(1+CX),
    # which is very nearly LINEAR at small collagen excess.  That put
    # enough gain around the COL -> TGFb -> MFB -> COL loop at the
    # healthy state to destroy it: fixed-point analysis of the slow
    # subsystem returned ONE root, the fibrotic one at COL 3.70.  The
    # unirradiated lung had no stable state to sit in, which is what made
    # every plan -- SBRT at 4.3 Gy included -- end at the same collagen.
    # Latent TGF-b1 is unfolded by integrin-transmitted force only once
    # the matrix is stiff enough to resist myofibroblast contraction
    # (Wipff 2007, Klingberg 2014, Hinz 2015), so the limb is properly a
    # THRESHOLD, not a ratio.  With a Hill-3 stiffness term the subsystem
    # is genuinely bistable:
    #     healthy      TGFB 0.0401  MFB 0.00057  COL 1.0000   stable
    #     separatrix   TGFB 0.2115  MFB 0.3563   COL 1.5864   UNSTABLE
    #     fibrotic     TGFB 0.7767  MFB 0.9692   COL 2.9340   stable
    # and the critical sustained extra activation is 0.0207, which with
    # KTGFR 0.012 puts the radiation-alone latching dose at 44 Gy in
    # 2 Gy fractions -- the isodose at which radiographic fibrosis is
    # actually seen, not a fitted number.
    KTGF0=0.0100,   # constitutive latent-TGFb activation
    KTGFR=0.0120,   # radiation/ROS-driven activation per Gy3/d
    KTGFI=0.0500,   # inflammation-driven activation
    KTGFM=0.1600,   # mechanotransduction (stiffness) driven activation
    KMCOL=1.0000,   # collagen excess at half-max force transmission
    HTGFM=3.0000,   # Hill coefficient of the stiffness threshold
    KTGFA=0.0450,   # myofibroblast autocrine TGFb
    KTGFEL=0.2500,

    # ---- myofibroblast ------------------------------------------------
    KMFB=1.6000,
    KMTGF=0.5500,   # half-max TGFb for differentiation
    HMFB=4.0000,    # Hill coefficient; >2 is what buys the bistability
    KMFBS=0.5000,   # stiffness boost to differentiation
    KMFBD=0.0800,   # myofibroblast apoptosis (1/d)

    # ---- collagen -----------------------------------------------------
    # [DEFECT 6] KCOL0/KCOL/KCDEG were 0.0191/0.020/0.020, a collagen
    # time-constant of 1/KCDEG = 50 d.  Crossing the separatrix (COL 1.0
    # -> 1.586) at the myofibroblast density a 60 Gy plan sustains then
    # takes 82 d, which is LONGER THAN THE COURSE -- so 60 Gy/30 and
    # protons latched no fibrosis at all (Vfib 0.0%) while SBRT, whose
    # 4.2-day course delivers 8x the critical activation rate, latched
    # 3.3%.  The model said hypofractionation is fibrogenic and
    # conventional fractionation is not, which is the opposite of the
    # clinic.  All three coefficients are scaled by 2.5 together, so every
    # fixed point above is algebraically unchanged and only the time
    # constant moves, 50 d -> 20 d.
    KCOL0=0.0477625,  # constitutive synthesis, set so COL* = 1.0000
    KCOL=0.0500,    # myofibroblast-driven synthesis
    KCDEG=0.0500,   # MMP-mediated degradation
    KTIMP=1.2000,   # TGFb-driven TIMP inhibition of degradation
    # [DEFECT 11] with a single remodellable collagen pool, pirfenidone
    # started on day 365 -- a year after the switch had flipped --
    # abolished fibrosis completely (Vfib 8.0% -> 0.0%, FIB730 0.155 ->
    # 0.000), i.e. the model claimed an antifibrotic REVERSES established
    # radiation fibrosis.  It does not; pirfenidone slows progression.
    # What was missing is that LOX/LOXL2 cross-linked collagen is not an
    # MMP substrate: mature matrix is mechanical MEMORY.  Labile collagen
    # (COL) now matures irreversibly into a cross-linked pool (XCOL) that
    # turns over with a 1.8-year half-life, and stiffness/mechano-
    # transduction/DLCO all read the SUM.  The fibrotic state is then
    # held up by matrix the drug cannot reach, which is what converts the
    # antifibrotic from a dose-response drug into a race against a clock.
    # KXL had to be SLOW.  At 0.050/d cross-linking matured within weeks,
    # so the cross-linked pool was already carrying most of the stiffness
    # while the switch was still deciding.  That collapsed the critical
    # activation from 0.0207 to 0.0041 and dragged the fibrotic threshold
    # from 44 Gy down below 30 Gy -- every bin above 30 Gy latched, which
    # is not what the isodose-delimited clinical picture looks like.  LOX
    # maturation of newly deposited collagen takes months, and at 0.008/d
    # the labile loop alone decides whether the switch flips while the
    # cross-linked pool accumulates AFTERWARDS and locks it in.  The
    # threshold is therefore set by the same labile loop as before, and
    # the irreversibility is a consequence rather than a cause.
    KXL=0.0080,     # LOX-mediated maturation of labile -> cross-linked
    XMAX=1.0000,    # maximum cross-linkable collagen density
    KXDEG=0.0015,   # cross-linked collagen turnover, tau ~ 1.8 y

    # ---- function mapping --------------------------------------------
    PEC=0.7000,
    PSURF=0.3000,
    KPF=0.5000,     # permeability half-effect on diffusing capacity
    KCF=1.5000,     # collagen half-effect on diffusing capacity
    KCV=2.5000,     # collagen half-effect on vital capacity
    KPV=1.2000,

    # ---- clinical scoring --------------------------------------------
    # [DEFECT 5] the index used PERM/(1+PERM) and CYT/(1+CYT), i.e. it
    # imposed its own saturation on top of the states' own saturation.
    # That is an assumption about the volume effect smuggled into the
    # read-out: it flattens any high-dose bin to the same contribution
    # and would have produced a large LKB n by construction.  The index
    # is now LINEAR in the states, so whatever volume effect comes out
    # comes out of the biology.
    A1=1.0000,      # weight of oedema in the pneumonitis index
    A2=0.6000,      # weight of cytokine tone in the pneumonitis index
    DLREF=85.0,     # reference baseline DLCO (%pred) for reserve scaling
    PNI50=1.2900,   # fitted to QUANTEC (rili_calibration.py)
    PNISL=0.4100,   # logistic slope on ln(PNI)

    # ---- drug PK ------------------------------------------------------
    # prednisolone: 1-cmt, oral; a continuous daily input rate is used
    # because the PD read-out (NF-kB transcription) integrates over days
    KADEX=6.0000, CLDEX=8.0000, VDEX=50.0000,
    IMXDEX=0.8500, IC50DEX=0.0300,      # mg/L free prednisolone
    # [DEFECT 12] the steroid reached the oedema state only through the
    # cytokine arm, and because PERM is carried mostly by the epithelial
    # term the whole benefit was a 16% fall in the index -- far less than
    # the clinical response.  Glucocorticoids also up-regulate ENaC and
    # Na,K-ATPase and directly accelerate alveolar fluid clearance
    # (Folkesson & Matthay), which is a second, faster route to the same
    # state and does not disturb the latency structure.
    EMXCLR=1.2000,  # maximal fold-increase in alveolar fluid clearance

    # pirfenidone 801 mg tid
    KAPIR=4.0000, CLPIR=60.0000, VPIR=70.0000,
    IMXPIR=0.6000, IC50PIR=3.0000,      # mg/L

    # nintedanib 150 mg bid
    KANIN=1.5000, CLNIN=1390.0, VNIN=1050.0,
    IMXNIN=0.6500, IC50NIN=0.0150,      # mg/L

    # durvalumab 1500 mg q4w IV
    CLDUR=0.2320, VDUR=3.6000, QDUR=0.4800, VPDUR=3.0000,
    # [DEFECT 7] EMXDUR was 0.90, taking the fast-loop gain to 0.78
    # against an elimination of 0.15.  Adding durvalumab therefore
    # detonated the acute loop -- DLCO nadir 0.56 %pred in every arm --
    # instead of modulating it.  0.35 keeps the loop sub-critical and
    # reproduces the size of the PACIFIC pneumonitis excess.
    EMXDUR=0.3500, EC50DUR=8.0000,      # mg/L; gain on the fast loop
    KTUMIMM=0.0060,                     # immune tumour kill (1/d)

    # lisinopril 20 mg qd
    KAACE=1.0000, CLACE=6.0000, VACE=100.0,
    IMXACE=0.4500, IC50ACE=0.0200,      # mg/L

    # amifostine: handled algebraically (see protect_amifostine)
    PFAMI=0.6500,   # maximal radioprotection factor of lung
    THALFAMI=8.0,   # plasma half-life, minutes
    PFAVA=0.4500,   # avasopasem (SOD mimetic) maximal protection
    THALFAVA=27.0,  # minutes

    # ---- tumour -------------------------------------------------------
    KGROW=0.0347,   # 1/d  (NSCLC clonogen doubling ~ 20 d)
    TUMLN0=20.72,   # ln(clonogens) at start = ln(1e9)
)


# ---------------------------------------------------------------------
# 3.  SCENARIO / PLAN DEFINITIONS
# ---------------------------------------------------------------------
# A plan is (DRX, NFX, v[0..5]).  v must sum to 1 (whole-lung DVH).
PLANS = {
    # conventional NSCLC IMRT, 60 Gy / 30 fx
    "imrt60": dict(DRX=60.0, NFX=30, v=[0.40, 0.20, 0.12, 0.11, 0.09, 0.08]),
    # 3D conformal, few beams: less low-dose bath, bigger hot volume
    "crt60": dict(DRX=60.0, NFX=30, v=[0.52, 0.13, 0.08, 0.08, 0.09, 0.10]),
    # proton / highly conformal: spares low dose AND high dose
    "proton60": dict(DRX=60.0, NFX=30, v=[0.58, 0.16, 0.09, 0.07, 0.055, 0.045]),
    # peripheral SBRT 54 Gy / 3 fx
    "sbrt54": dict(DRX=54.0, NFX=3, v=[0.860, 0.075, 0.032, 0.018, 0.010, 0.005]),
    # central SBRT-like moderate hypofractionation 60 Gy / 8 fx
    "hypo60": dict(DRX=60.0, NFX=8, v=[0.700, 0.130, 0.070, 0.045, 0.033, 0.022]),
    # dose escalation arm of RTOG 0617
    "imrt74": dict(DRX=74.0, NFX=37, v=[0.34, 0.20, 0.13, 0.12, 0.11, 0.10]),
    # two plans matched on MEAN LUNG DOSE but not on histogram shape
    "bath": dict(DRX=60.0, NFX=30, v=[0.06, 0.55, 0.29, 0.05, 0.03, 0.02]),
    "hot": dict(DRX=60.0, NFX=30, v=[0.600, 0.110, 0.050, 0.050, 0.070, 0.120]),
}


def make_plan(name):
    p = PLANS[name]
    v = list(p["v"])
    s = sum(v)
    v = [x / s for x in v]            # renormalise defensively
    return p["DRX"], p["NFX"], v


# ---------------------------------------------------------------------
# 4.  STATE VECTOR LAYOUT
# ---------------------------------------------------------------------
SB = ["AT2", "DOOM", "EC", "SURF", "PERM", "CYT", "TGFB", "MFB", "COL",
      "XCOL"]
NSB = len(SB)
IX = {n: i for i, n in enumerate(SB)}


def sidx(b, name):
    return b * NSB + IX[name]


OFF = NB * NSB                     # 60
I_CYTS = OFF + 0
I_BEDL = OFF + 1                   # .. +6
I_BEDT = OFF + 7
I_TUMLN = OFF + 8
I_DEXD = OFF + 9
I_DEXC = OFF + 10
I_PIRD = OFF + 11
I_PIRC = OFF + 12
I_NIND = OFF + 13
I_NINC = OFF + 14
I_DURC = OFF + 15
I_DURP = OFF + 16
I_ACED = OFF + 17
I_ACEC = OFF + 18
NST = OFF + 19                     # 79


def initial_state(p=None):
    y = [0.0] * NST
    for b in range(NB):
        y[sidx(b, "AT2")] = 1.0
        y[sidx(b, "DOOM")] = 0.0
        y[sidx(b, "EC")] = 1.0
        y[sidx(b, "SURF")] = 1.0
        y[sidx(b, "PERM")] = 0.0
        y[sidx(b, "CYT")] = 0.0
        y[sidx(b, "TGFB")] = P["KTGF0"] / P["KTGFEL"]     # 0.05
        y[sidx(b, "MFB")] = 0.0
        y[sidx(b, "COL")] = 1.0
        y[sidx(b, "XCOL")] = 0.0
    y[I_TUMLN] = (p or P)["TUMLN0"]
    return y


# ---------------------------------------------------------------------
# 5.  RADIOPROTECTOR TIMING
# ---------------------------------------------------------------------
def protect_amifostine(on, delay_min, kind="ami"):
    """A thiol radioprotector only protects tissue that is irradiated
    WHILE the drug is present.  Amifostine's active metabolite WR-1065
    has a plasma half-life of ~8 min, so the protection factor decays
    exponentially with the drug-to-beam interval.  This is the single
    parameter that separates the trials which found benefit from those
    which did not."""
    if not on:
        return 0.0
    if kind == "ami":
        pf, th = P["PFAMI"], P["THALFAMI"]
    else:
        pf, th = P["PFAVA"], P["THALFAVA"]
    return pf * math.exp(-math.log(2.0) * delay_min / th)


# ---------------------------------------------------------------------
# 6.  DOSING SCHEDULES
# ---------------------------------------------------------------------
class Regimen(object):
    """Continuous-rate oral inputs plus discrete IV boluses.

    Oral drugs enter as a mg/day rate into the depot.  This is exact for
    the daily-averaged exposure and adequate here because every PD
    time-constant downstream is >= days.  Durvalumab is a true q28d
    bolus and is applied as a discrete event between integration steps.
    """

    def __init__(self, **kw):
        self.dex_mg = kw.get("dex_mg", 0.0)          # mg/day prednisolone
        self.dex_t0 = kw.get("dex_t0", 1e9)
        self.dex_dur = kw.get("dex_dur", 14.0)       # days at full dose
        self.dex_taper = kw.get("dex_taper", 42.0)   # days of linear taper
        self.pir_mg = kw.get("pir_mg", 0.0)
        self.pir_t0 = kw.get("pir_t0", 1e9)
        self.pir_t1 = kw.get("pir_t1", 1e9)
        self.nin_mg = kw.get("nin_mg", 0.0)
        self.nin_t0 = kw.get("nin_t0", 1e9)
        self.nin_t1 = kw.get("nin_t1", 1e9)
        self.ace_mg = kw.get("ace_mg", 0.0)
        self.ace_t0 = kw.get("ace_t0", 1e9)
        self.ace_t1 = kw.get("ace_t1", 1e9)
        self.dur_mg = kw.get("dur_mg", 0.0)
        self.dur_t0 = kw.get("dur_t0", 1e9)
        self.dur_n = kw.get("dur_n", 13)
        self.dur_q = kw.get("dur_q", 28.0)
        self.ami_on = kw.get("ami_on", False)
        self.ami_delay = kw.get("ami_delay", 20.0)   # min before beam
        self.ava_on = kw.get("ava_on", False)
        self.ava_delay = kw.get("ava_delay", 30.0)

    def dex_rate(self, t):
        if self.dex_mg <= 0 or t < self.dex_t0:
            return 0.0
        dt = t - self.dex_t0
        if dt <= self.dex_dur:
            return self.dex_mg
        if dt <= self.dex_dur + self.dex_taper:
            frac = 1.0 - (dt - self.dex_dur) / self.dex_taper
            return self.dex_mg * max(frac, 0.0)
        return 0.0

    @staticmethod
    def _win(t, mg, t0, t1):
        return mg if (t0 <= t <= t1) else 0.0

    def pir_rate(self, t):
        return self._win(t, self.pir_mg, self.pir_t0, self.pir_t1)

    def nin_rate(self, t):
        return self._win(t, self.nin_mg, self.nin_t0, self.nin_t1)

    def ace_rate(self, t):
        return self._win(t, self.ace_mg, self.ace_t0, self.ace_t1)

    def dur_events(self):
        if self.dur_mg <= 0:
            return []
        return [self.dur_t0 + i * self.dur_q for i in range(self.dur_n)]


# ---------------------------------------------------------------------
# 7.  RIGHT-HAND SIDE
# ---------------------------------------------------------------------
class Sim(object):
    def __init__(self, plan, reg, DLCO0=85.0, FVC0=95.0, perf=None,
                 gimm_base=1.0, params=None):
        self.p = dict(P)
        if params:
            self.p.update(params)
        self.DRX, self.NFX, self.v = plan
        self.reg = reg
        self.DLCO0 = DLCO0
        self.FVC0 = FVC0

        # course length: FXWK fractions per week
        self.TCOURSE = self.NFX / self.p["FXWK"] * 7.0

        # per-bin total lung BED (alpha/beta = ABL) and its delivery rate
        Db = bin_doses(self.DRX)
        self.Db = Db
        d_fx = [D / self.NFX for D in Db]
        self.BEDLtot = [Db[b] * (1.0 + d_fx[b] / self.p["ABL"])
                        for b in range(NB)]
        self.RB = [x / self.TCOURSE for x in self.BEDLtot]

        # tumour receives the full prescription
        dT = self.DRX / self.NFX
        self.BEDTtot = self.DRX * (1.0 + dT / self.p["ABT"])
        self.RBT = self.BEDTtot / self.TCOURSE

        # perfusion weighting of each bin (emphysema = low perfusion where
        # tissue is destroyed, so dose there costs less gas exchange)
        if perf is None:
            perf = [1.0] * NB
        w = [self.v[b] * perf[b] for b in range(NB)]
        sw = sum(w)
        self.w = [x / sw for x in w]

        # radioprotector factor, fixed for the course
        self.PROT = min(0.95,
                        protect_amifostine(reg.ami_on, reg.ami_delay, "ami") +
                        protect_amifostine(reg.ava_on, reg.ava_delay, "ava"))
        self.gimm_base = gimm_base

    # -- helpers -------------------------------------------------------
    def rt_on(self, t):
        return 1.0 if (0.0 <= t <= self.TCOURSE) else 0.0

    def rhs(self, t, y):
        p = self.p
        dy = [0.0] * NST
        on = self.rt_on(t)

        # ---------------- PK ------------------------------------------
        dexd, dexc = y[I_DEXD], y[I_DEXC]
        pird, pirc = y[I_PIRD], y[I_PIRC]
        nind, ninc = y[I_NIND], y[I_NINC]
        durc, durp = y[I_DURC], y[I_DURP]
        aced, acec = y[I_ACED], y[I_ACEC]

        dy[I_DEXD] = self.reg.dex_rate(t) - p["KADEX"] * dexd
        dy[I_DEXC] = p["KADEX"] * dexd - p["CLDEX"] / p["VDEX"] * dexc
        dy[I_PIRD] = self.reg.pir_rate(t) - p["KAPIR"] * pird
        dy[I_PIRC] = p["KAPIR"] * pird - p["CLPIR"] / p["VPIR"] * pirc
        dy[I_NIND] = self.reg.nin_rate(t) - p["KANIN"] * nind
        dy[I_NINC] = p["KANIN"] * nind - p["CLNIN"] / p["VNIN"] * ninc
        dy[I_ACED] = self.reg.ace_rate(t) - p["KAACE"] * aced
        dy[I_ACEC] = p["KAACE"] * aced - p["CLACE"] / p["VACE"] * acec
        dy[I_DURC] = (-p["CLDUR"] / p["VDUR"] * durc
                      - p["QDUR"] / p["VDUR"] * durc
                      + p["QDUR"] / p["VPDUR"] * durp)
        dy[I_DURP] = (p["QDUR"] / p["VDUR"] * durc
                      - p["QDUR"] / p["VPDUR"] * durp)

        CDEX = dexc / p["VDEX"]
        CPIR = pirc / p["VPIR"]
        CNIN = ninc / p["VNIN"]
        CACE = acec / p["VACE"]
        CDUR = durc / p["VDUR"]

        # ---------------- drug effects --------------------------------
        # steroid: transcriptional suppression of the FAST loop only
        SDEX = 1.0 - p["IMXDEX"] * CDEX / (p["IC50DEX"] + CDEX)
        # second, independent steroid route: ENaC / Na,K-ATPase driven
        # alveolar fluid clearance
        ECLR = p["EMXCLR"] * CDEX / (p["IC50DEX"] + CDEX)
        # pirfenidone: TGFb signalling + collagen synthesis (SLOW loop)
        IPIR = p["IMXPIR"] * CPIR / (p["IC50PIR"] + CPIR)
        # nintedanib: fibroblast proliferation (SLOW loop)
        ININ = p["IMXNIN"] * CNIN / (p["IC50NIN"] + CNIN)
        # ACE inhibitor: latent TGFb activation (SLOW loop)
        IACE = p["IMXACE"] * CACE / (p["IC50ACE"] + CACE)
        # durvalumab: raises the GAIN of the fast loop (checkpoint release)
        GIMM = self.gimm_base * (1.0 + p["EMXDUR"] * CDUR /
                                 (p["EC50DUR"] + CDUR))

        # ---------------- systemic cytokine ---------------------------
        CYTS = y[I_CYTS]
        spill = 0.0

        # ---------------- per-bin biology -----------------------------
        for b in range(NB):
            AT2 = y[sidx(b, "AT2")]
            DOOM = y[sidx(b, "DOOM")]
            EC = y[sidx(b, "EC")]
            SURF = y[sidx(b, "SURF")]
            PERM = y[sidx(b, "PERM")]
            CYT = y[sidx(b, "CYT")]
            TGFB = y[sidx(b, "TGFB")]
            MFB = y[sidx(b, "MFB")]
            COL = y[sidx(b, "COL")]
            XCOL = y[sidx(b, "XCOL")]
            COLX = max(COL - 1.0, 0.0)
            # stiffness, mechanotransduction and every functional
            # read-out see labile PLUS cross-linked collagen
            CXT = COLX + XCOL

            RB = self.RB[b] * on * (1.0 - self.PROT)

            # --- AT2: radiation puts cells in a doomed pool that only
            # dies when it next attempts mitosis.  This delay, not the
            # dose, is what sets the 4-12 week latency of pneumonitis.
            JRAD = p["KDAM"] * RB * AT2
            JCYT = p["KDCYT"] * (CYT / (p["KMCYT"] + CYT)) * AT2
            grow = p["KREP"] * AT2 * max(1.0 - (AT2 + DOOM), 0.0)
            dy[sidx(b, "AT2")] = grow - JRAD - JCYT
            dy[sidx(b, "DOOM")] = JRAD + JCYT - p["KMIT"] * DOOM
            DTH = p["KMIT"] * DOOM          # dying-cell flux -> DAMPs

            # --- microvascular endothelium.  Repair is toward a CEILING
            # that falls with cumulative dose: capillary dropout and loss
            # of endothelial progenitors are not repaired, so a bin that
            # has been irradiated can never return to EC = 1.  This is
            # the permanent, dose-graded component of DLCO loss and it
            # exists in bins that never become fibrotic.
            ECMAX = 1.0 / (1.0 + p["KECIRR"] * y[I_BEDL + b])
            dy[sidx(b, "EC")] = (p["KECREP"] * (ECMAX - EC)
                                 - p["KECDAM"] * RB * EC
                                 - p["KECCYT"] * CYT * EC)

            # --- surfactant.  [DEFECT 9] production was proportional to
            # VIABLE AT2 only, so the alveolus recovered the instant the
            # beam stopped and the whole injury peaked at t = end of RT
            # (tPk 0-3 d) instead of the 4-12 weeks every clinical series
            # reports.  A lethally hit pneumocyte keeps making surfactant
            # until it attempts mitosis and dies -- that is the entire
            # reason late-responding tissue injury is DELAYED -- so the
            # producing pool is AT2 + DOOM, and the alveolar lining keeps
            # thinning for weeks after the last fraction.
            dy[sidx(b, "SURF")] = (p["KSYN"] * (AT2 + DOOM)
                                   - p["KSDEG"] * SURF)

            # --- alveolar-capillary leak / oedema
            dy[sidx(b, "PERM")] = (p["KPRMEC"] * (1.0 - EC)
                                   + p["KPRMEP"] * (1.0 - SURF)
                                   + p["KPCYT"] * CYT
                                   - p["KPCL"] * (1.0 + ECLR) * PERM)

            # --- FAST loop.  GIMM multiplies the whole production side,
            # so checkpoint blockade is a GAIN change, not an added term.
            amp = p["KAMP"] * CYT / (1.0 + CYT / p["CMAX"])
            dy[sidx(b, "CYT")] = (SDEX * GIMM * (p["KDAMP"] * DTH
                                                 + p["KDIR"] * RB
                                                 + amp)
                                  + p["KSYSIN"] * CYTS
                                  - p["KCE"] * CYT)
            spill += self.v[b] * CYT

            # --- SLOW loop.  Four independent activators of latent
            # TGF-beta1; only the mechanotransduction and autocrine ones
            # close a positive feedback, which is what makes the
            # fibrotic state self-sustaining.
            mech = (p["KTGFM"] * CXT ** p["HTGFM"]
                    / (p["KMCOL"] ** p["HTGFM"] + CXT ** p["HTGFM"]))
            act = (p["KTGF0"]
                   + p["KTGFR"] * RB
                   + p["KTGFI"] * CYT
                   + mech
                   + p["KTGFA"] * MFB)
            dy[sidx(b, "TGFB")] = act * (1.0 - IACE) - p["KTGFEL"] * TGFB

            drive = TGFB ** p["HMFB"] / (p["KMTGF"] ** p["HMFB"]
                                         + TGFB ** p["HMFB"])
            dy[sidx(b, "MFB")] = (p["KMFB"] * drive * (1.0 - IPIR)
                                  * (1.0 + p["KMFBS"] * CXT)
                                  * (1.0 - ININ)
                                  * max(1.0 - MFB, 0.0)
                                  - p["KMFBD"] * MFB)

            deg = p["KCDEG"] * COL / (1.0 + p["KTIMP"] * TGFB / (1.0 + TGFB))
            xflux = p["KXL"] * COLX * max(p["XMAX"] - XCOL, 0.0)
            dy[sidx(b, "COL")] = (p["KCOL0"]
                                  + p["KCOL"] * MFB * (1.0 - IPIR)
                                  - deg - xflux)
            dy[sidx(b, "XCOL")] = xflux - p["KXDEG"] * XCOL

            # --- cumulative BED (report only)
            dy[I_BEDL + b] = self.RB[b] * on

        dy[I_CYTS] = p["KSOUT"] * spill - p["KCLS"] * CYTS

        # ---------------- tumour --------------------------------------
        dy[I_BEDT] = self.RBT * on
        rt_kill = p["ALPHAT"] * self.RBT * on
        imm_kill = p["KTUMIMM"] * CDUR / (p["EC50DUR"] + CDUR)
        # [DEFECT 8] repopulation was applied unconditionally, so a plan
        # that drove the clonogen number to 0.85 cells still "regrew" to
        # e^23 by day 730 and every TCP read 0%.  Below one surviving
        # clonogen there is nothing left to divide; growth is gated off
        # and TCP is read from the nadir, which is the only point at
        # which Poisson eradication is defined.
        alive = 1.0 if y[I_TUMLN] > 0.0 else 0.0
        dy[I_TUMLN] = p["KGROW"] * alive - rt_kill - imm_kill

        return dy

    # -- derived read-outs --------------------------------------------
    def readouts(self, t, y):
        p = self.p
        DL = 0.0
        FV = 0.0
        PNI = 0.0
        FIB = 0.0
        colmax = 0.0
        vfib = 0.0
        for b in range(NB):
            EC = max(y[sidx(b, "EC")], 1e-9)
            SURF = max(y[sidx(b, "SURF")], 1e-9)
            PERM = max(y[sidx(b, "PERM")], 0.0)
            CYT = max(y[sidx(b, "CYT")], 0.0)
            CXT = max(y[sidx(b, "COL")] - 1.0, 0.0) + y[sidx(b, "XCOL")]
            fd = (EC ** p["PEC"]) * (SURF ** p["PSURF"]) \
                / (1.0 + PERM / p["KPF"]) / (1.0 + CXT / p["KCF"])
            DL += self.w[b] * fd
            fv = 1.0 / (1.0 + CXT / p["KCV"]) / (1.0 + PERM / p["KPV"])
            FV += self.w[b] * fv
            PNI += self.v[b] * (p["A1"] * PERM + p["A2"] * CYT)
            FIB += self.v[b] * CXT
            colmax = max(colmax, CXT)
            # a bin is counted as CONVERTED once it is past the
            # separatrix (COL 1.586); below that it is still on its way
            # back to the healthy root
            if y[sidx(b, "COL")] + y[sidx(b, "XCOL")] > 1.586:
                vfib += self.v[b]
        # low baseline reserve amplifies the SAME absolute injury
        RF = math.sqrt(p["DLREF"] / max(self.DLCO0, 20.0))
        PNIe = PNI * RF
        return dict(t=t,
                    DLCO=self.DLCO0 * DL,
                    FVC=self.FVC0 * FV,
                    PNI=PNI, PNIe=PNIe,
                    FIB=FIB, COLMAX=colmax, VFIB=vfib,
                    NTCP=ntcp(PNIe, p),
                    TUMLN=y[I_TUMLN],
                    CYTS=y[I_CYTS])


def ntcp(pnie, p):
    """Probability of CTCAE grade >= 2 radiation pneumonitis.  The
    logistic is on ln(index), so the model produces a sigmoid NTCP curve
    in dose without one being assumed."""
    if pnie <= 1e-9:
        return 0.0
    z = (math.log(pnie) - math.log(p["PNI50"])) / p["PNISL"]
    return 1.0 / (1.0 + math.exp(-z))


def ctcae_grade(pnie, dlco_drop_pct, p):
    """CTCAE v5 radiation-pneumonitis grade, computed as an OUTPUT of the
    continuous state rather than supplied as an input."""
    if pnie < 0.45:
        return 0
    if pnie < 0.95:
        return 1                      # radiographic only
    if pnie < 1.55:
        return 2                      # symptomatic, steroids indicated
    if pnie < 2.30 or dlco_drop_pct < 45.0:
        return 3                      # oxygen indicated
    return 4                          # life-threatening


# ---------------------------------------------------------------------
# 8.  INTEGRATOR  (fixed-step RK4 with discrete IV events)
# ---------------------------------------------------------------------
def integrate(sim, tend, dt=0.05, tsave=1.0):
    y = initial_state(sim.p)
    t = 0.0
    events = sorted(sim.reg.dur_events())
    ei = 0
    out = []
    nextsave = 0.0
    dur_dose = sim.reg.dur_mg

    while t < tend - 1e-12:
        if t >= nextsave - 1e-9:
            out.append(sim.readouts(t, y))
            nextsave += tsave

        # apply any IV bolus that falls in [t, t+dt)
        while ei < len(events) and events[ei] < t + dt - 1e-12:
            if events[ei] >= t - 1e-12:
                y[I_DURC] += dur_dose
            ei += 1

        h = min(dt, tend - t)
        k1 = sim.rhs(t, y)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(NST)]
        k2 = sim.rhs(t + 0.5 * h, y2)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(NST)]
        k3 = sim.rhs(t + 0.5 * h, y3)
        y4 = [y[i] + h * k3[i] for i in range(NST)]
        k4 = sim.rhs(t + h, y4)
        y = [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(NST)]
        # states that are physically non-negative
        for i in range(NST):
            if i != I_TUMLN and y[i] < 0.0:
                y[i] = 0.0
        t += h

    out.append(sim.readouts(t, y))
    return out, y


# ---------------------------------------------------------------------
# 9.  SUMMARY OF A RUN
# ---------------------------------------------------------------------
def summarise(sim, out, label=""):
    peak = max(out, key=lambda r: r["PNIe"])
    dl_nadir = min(out, key=lambda r: r["DLCO"])
    tum_nadir = min(r["TUMLN"] for r in out)
    d365 = min(out, key=lambda r: abs(r["t"] - 365.0))
    d730 = out[-1]
    drop_peak = 100.0 * (sim.DLCO0 - dl_nadir["DLCO"]) / sim.DLCO0
    return dict(
        label=label,
        MLD=mean_lung_dose(sim.v, sim.DRX),
        V5=100 * vx(sim.v, sim.DRX, 5.0),
        V20=100 * vx(sim.v, sim.DRX, 20.0),
        V40=100 * vx(sim.v, sim.DRX, 40.0),
        gEUD1=gEUD(sim.v, sim.DRX, 1.0),
        PNIpk=peak["PNIe"],
        tPNIpk=peak["t"] - sim.TCOURSE,
        NTCP=100 * peak["NTCP"],
        GRADE=ctcae_grade(peak["PNIe"], drop_peak, sim.p),
        DLnadir=dl_nadir["DLCO"],
        tDLnadir=dl_nadir["t"] - sim.TCOURSE,
        DLCO365=d365["DLCO"],
        DLCO730=d730["DLCO"],
        dDLCO365=100 * (d365["DLCO"] - sim.DLCO0) / sim.DLCO0,
        FVC365=d365["FVC"],
        FIB365=d365["FIB"],
        FIB730=d730["FIB"],
        VFIB365=100 * d365["VFIB"],
        VFIB730=100 * d730["VFIB"],
        COLMAX730=d730["COLMAX"],
        TUMLN=tum_nadir,
        TCP=100 * math.exp(-math.exp(min(tum_nadir, 50.0))),
        BEDT=sim.BEDTtot,
        TCOURSE=sim.TCOURSE,
    )


HDR = ("{:<34s} {:>6s} {:>6s} {:>6s} {:>7s} {:>6s} {:>6s} {:>2s} "
       "{:>7s} {:>7s} {:>6s} {:>7s}")
ROW = ("{label:<34s} {MLD:6.2f} {V20:6.1f} {V40:6.1f} "
       "{PNIpk:7.4f} {tPNIpk:6.1f} {NTCP:6.1f} {GRADE:2d} {DLnadir:7.2f} "
       "{dDLCO365:7.2f} {VFIB730:6.1f} {FIB730:7.4f}")


def header():
    print(HDR.format("scenario", "MLD", "V20", "V40", "PNIpk", "tPk",
                     "NTCP%", "G", "DLnad", "dDL365", "Vfib%", "FIB730"))
    print("-" * 126)


def run(plan_name, reg, label, tend=730.0, dt=0.05, **kw):
    sim = Sim(make_plan(plan_name), reg, **kw)
    out, yend = integrate(sim, tend, dt=dt)
    s = summarise(sim, out, label)
    return sim, out, s


# =====================================================================
# 10.  SCENARIOS
# =====================================================================
def main():
    print(__doc__ or "")
    print("=" * 128)
    print("RADIATION-INDUCED LUNG INJURY -- reference implementation")
    print("79 ODEs = 6 dose bins x 10 states + 19 global")
    print("=" * 128)

    NONE = Regimen()
    results = {}

    # ---- A. natural history across plans -----------------------------
    print("\n[A] NATURAL HISTORY -- no drug, plan is the only variable")
    header()
    for nm, lab in [("sbrt54", "1  SBRT 54/3 peripheral"),
                    ("hypo60", "2  hypofx 60/8"),
                    ("proton60", "3  proton 60/30"),
                    ("imrt60", "4  IMRT 60/30 (reference)"),
                    ("crt60", "5  3D-CRT 60/30"),
                    ("imrt74", "6  IMRT 74/37 (RTOG 0617 high)")]:
        sim, out, s = run(nm, NONE, lab)
        results[nm] = (sim, out, s)
        print(ROW.format(**s))

    # ---- B. the matched-mean-dose experiment -------------------------
    print("\n[B] SAME MEAN LUNG DOSE, DIFFERENT HISTOGRAM SHAPE")
    header()
    for nm, lab in [("bath", "7  low-dose bath (V5 high)"),
                    ("hot", "8  focal hot spot (V40 high)")]:
        sim, out, s = run(nm, NONE, lab)
        results[nm] = (sim, out, s)
        print(ROW.format(**s))

    # ---- C. steroid --------------------------------------------------
    print("\n[C] CORTICOSTEROID -- acts on the FAST loop only")
    header()
    for t0, lab in [(70.0, "9  pred 60mg from d70 (onset)"),
                    (56.0, "10 pred 60mg from d56 (early)"),
                    (120.0, "11 pred 60mg from d120 (late)")]:
        reg = Regimen(dex_mg=60.0, dex_t0=t0)
        sim, out, s = run("imrt60", reg, lab)
        results["ster%d" % t0] = (sim, out, s)
        print(ROW.format(**s))

    # ---- D. antifibrotics -------------------------------------------
    print("\n[D] ANTIFIBROTICS -- act on the SLOW loop only")
    header()
    for kw, lab in [
        (dict(pir_mg=2403.0, pir_t0=0.0, pir_t1=365.0),
         "12 pirfenidone d0-365"),
        (dict(pir_mg=2403.0, pir_t0=180.0, pir_t1=545.0),
         "13 pirfenidone d180-545 (late)"),
        (dict(nin_mg=300.0, nin_t0=0.0, nin_t1=365.0),
         "14 nintedanib d0-365"),
        (dict(ace_mg=20.0, ace_t0=0.0, ace_t1=730.0),
         "15 lisinopril d0-730"),
    ]:
        sim, out, s = run("imrt60", Regimen(**kw), lab)
        results[lab[:2].strip()] = (sim, out, s)
        print(ROW.format(**s))

    # ---- E. radioprotector timing -----------------------------------
    print("\n[E] RADIOPROTECTOR -- protection requires TEMPORAL COINCIDENCE")
    header()
    for kw, lab in [
        (dict(ami_on=True, ami_delay=0.0), "16 amifostine 0 min pre-beam"),
        (dict(ami_on=True, ami_delay=15.0), "17 amifostine 15 min pre-beam"),
        (dict(ami_on=True, ami_delay=30.0), "18 amifostine 30 min pre-beam"),
        (dict(ami_on=True, ami_delay=60.0), "19 amifostine 60 min pre-beam"),
        (dict(ava_on=True, ava_delay=30.0), "20 avasopasem 30 min pre-beam"),
    ]:
        sim, out, s = run("imrt60", Regimen(**kw), lab)
        results[lab[:2].strip()] = (sim, out, s)
        print(ROW.format(**s))

    # ---- F. consolidation immunotherapy ------------------------------
    print("\n[F] CONSOLIDATION DURVALUMAB (PACIFIC) x mean lung dose")
    header()
    for nm, kw, lab in [
        ("imrt60", dict(), "21 chemoRT alone   MLD~18"),
        ("imrt60", dict(dur_mg=1500.0, dur_t0=56.0),
         "22 +durvalumab     MLD~18"),
        ("proton60", dict(), "23 chemoRT alone   MLD~13"),
        ("proton60", dict(dur_mg=1500.0, dur_t0=56.0),
         "24 +durvalumab     MLD~13"),
        ("imrt74", dict(), "25 chemoRT alone   MLD~22"),
        ("imrt74", dict(dur_mg=1500.0, dur_t0=56.0),
         "26 +durvalumab     MLD~22"),
    ]:
        sim, out, s = run(nm, Regimen(**kw), lab)
        results[lab[:2].strip()] = (sim, out, s)
        print(ROW.format(**s))

    # ---- G. baseline reserve ----------------------------------------
    print("\n[G] BASELINE RESERVE -- identical plan, different lung")
    header()
    for dl, lab in [(100.0, "27 DLCO0 100%pred"),
                    (85.0, "28 DLCO0  85%pred"),
                    (60.0, "29 DLCO0  60%pred (COPD)"),
                    (45.0, "30 DLCO0  45%pred (severe)")]:
        sim, out, s = run("imrt60", NONE, lab, DLCO0=dl)
        results[lab[:2].strip()] = (sim, out, s)
        print(ROW.format(**s))

    # ---- H. therapeutic ratio ---------------------------------------
    # The tumour and the lung read the SAME physical dose through two
    # different alpha/beta ratios, so fractionation moves them by
    # different amounts.  UCP = TCP x (1 - NTCP) is the only endpoint
    # that can see both at once.
    print("\n[H] THERAPEUTIC RATIO -- tumour BED10 vs lung BED3")
    print("%-30s %7s %7s %8s %8s %7s %7s %7s"
          % ("plan", "d/fx", "BEDT10", "MLD", "lungB3*", "TCP%", "NTCP%",
             "UCP%"))
    print("-" * 92)
    for nm, n0, lab in [("imrt60", 20.72, "31 IMRT 60/30   stage III"),
                        ("imrt74", 20.72, "32 IMRT 74/37   stage III"),
                        ("hypo60", 20.72, "33 hypofx 60/8  stage III"),
                        ("proton60", 20.72, "34 proton 60/30 stage III"),
                        ("sbrt54", 18.42, "35 SBRT 54/3    stage I"),
                        ("hypo60", 18.42, "36 hypofx 60/8  stage I")]:
        sim, out, s = run(nm, NONE, lab, params=dict(TUMLN0=n0))
        results[lab[:2].strip()] = (sim, out, s)
        ucp = s["TCP"] / 100.0 * (1.0 - s["NTCP"] / 100.0) * 100.0
        print("%-30s %7.2f %7.1f %8.2f %8.1f %7.1f %7.1f %7.1f"
              % (lab, sim.DRX / sim.NFX, sim.BEDTtot, s["MLD"],
                 sim.BEDLtot[5], s["TCP"], s["NTCP"], ucp))

    # ---- I. the antifibrotic therapeutic WINDOW ---------------------
    # The slow loop is bistable, so an antifibrotic is not a
    # dose-response drug but a race against the separatrix.
    print("\n[I] ANTIFIBROTIC TIME WINDOW -- pirfenidone start day swept")
    print("%-34s %8s %8s %8s %9s"
          % ("start day", "Vfib%", "FIB730", "dDL365", "dDL730"))
    print("-" * 74)
    for t0 in (0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 365.0, 1e9):
        lab = ("37 pirfenidone from d%-6.0f" % t0) if t0 < 1e8 \
            else "37 no antifibrotic        "
        reg = Regimen(pir_mg=2403.0, pir_t0=t0, pir_t1=t0 + 365.0)
        sim, out, s = run("imrt60", reg, lab)
        d730 = out[-1]
        print("%-34s %8.1f %8.4f %8.2f %9.2f"
              % (lab, s["VFIB730"], s["FIB730"], s["dDLCO365"],
                 100 * (d730["DLCO"] - sim.DLCO0) / sim.DLCO0))

    # ---- J. steroid symptom response --------------------------------
    # NTCP reads the PEAK of the index, which a steroid started at onset
    # can barely move.  What a steroid actually does is collapse the
    # index over the following weeks, which is what the patient feels
    # and what no fibrosis endpoint can see.
    print("\n[J] STEROID RESPONSE vs UNTREATED COUNTERFACTUAL (MLD 17.8)")
    print("%-28s %9s %9s %9s %9s %8s"
          % ("", "PNI d+0", "PNI d+14", "PNI d+28", "PNI d+56", "Vfib%"))
    print("-" * 78)
    T0 = 56.0
    for kw, lab in [(dict(), "38 no steroid"),
                    (dict(dex_mg=30.0, dex_t0=T0), "39 prednisolone 30 mg/d"),
                    (dict(dex_mg=60.0, dex_t0=T0), "40 prednisolone 60 mg/d")]:
        sim, out, s = run("imrt60", Regimen(**kw), lab)
        vals = []
        for dt in (0.0, 14.0, 28.0, 56.0):
            r = min(out, key=lambda r, d=dt: abs(r["t"] - (T0 + d)))
            vals.append(r["PNIe"])
        print("%-28s %9.4f %9.4f %9.4f %9.4f %8.1f"
              % (lab, vals[0], vals[1], vals[2], vals[3], s["VFIB730"]))

    return results


if __name__ == "__main__":
    main()
