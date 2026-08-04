#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cos_reference_check.py
======================
Independent re-implementation, in Python/scipy, of every ODE in
`cos_mrgsolve_model.R` (Controlled Ovarian Stimulation QSP model).

WHY THIS FILE EXISTS
--------------------
The R model cannot be executed in the build environment (no R/mrgsolve), so
every equation is re-typed here from the same specification and integrated
with scipy.  Anything the R file claims in a comment is computed here, and
the numbers printed by this script are the numbers quoted in README.md.
Writing the system twice is also the cheapest available bug detector: five
structural defects were found this way (see README, "검증에서 드러난 결함").

THE MODEL IN ONE PARAGRAPH
--------------------------
Ten follicle "slots" (equal-probability quantiles of a log-normal FSH-threshold
distribution, multiplicity AFC/10 each) grow if and only if the prevailing
total FSH concentration exceeds their own threshold.  Nothing in the parameter
list says how many follicles ovulate: mono-follicular selection in an
unstimulated cycle and 20-follicle multi-follicular growth under 150 IU/d of
recombinant FSH are the SAME equations with the FSH decline present or absent.
Granulosa mass is the single state variable that oestradiol, trigger-day
progesterone, oocyte yield and VEGF are all read from - which is why FSH dose
cannot separate yield from OHSS.  The LHCGR signal is read TWICE with different
kernels: oocyte maturation needs only the leading EDGE (a commitment switch plus
a ~34 h autonomous clock), while ovarian VEGF output is proportional to the
AREA under receptor occupancy.  hCG supplies 8 days of area; a GnRH-agonist
trigger supplies the edge and ~6 h of area.  Equal maturity, no OHSS - from one
asymmetry, not from a severity parameter.

Run:  python3 cos_reference_check.py            (full report)
      python3 cos_reference_check.py --brief    (headline table only)
"""

import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# 0.  index map. These are python state indices; the R model names its
#     compartments instead ($CMT in cos_mrgsolve_model.R), so the two files
#     share equations and parameter values, not slot numbers.
# ----------------------------------------------------------------------------
NF = 10
iD, iG, iV = 0, 10, 20
(FSHDEP, FSHC, CORID, CORIC, ANTD, ANTC, HCG, AGOD, AGOC, CAB, LET) = range(30, 41)
(SLH, LH, FSHE, E2, INHB, AMH, P4D, P4, CL, RS) = range(41, 51)
(VEGF, PERM, ASC, VP, PHCG, EXPLH, MCLK, OOC, MIIC, AUCE2, EXPV, ROV) = range(51, 63)
HCGD = 63          # hCG subcutaneous/IM depot
GRR = 64           # GnRH-receptor availability (0-1)
NST = 65

# equal-probability standard-normal quantile midpoints, i = 1..10
ZQ = np.array([-1.6449, -1.0364, -0.6745, -0.3853, -0.1257,
               0.1257, 0.3853, 0.6745, 1.0364, 1.6449])


def hill(x, k, n):
    x = max(x, 0.0)
    if x == 0.0:
        return 0.0
    xn = (x / k) ** n
    return xn / (1.0 + xn)


# ----------------------------------------------------------------------------
# 1.  parameters
# ----------------------------------------------------------------------------
P = dict(
    # --- follicle threshold distribution ------------------------------------
    T50=9.0,        # IU/L   median FSH threshold of the antral cohort
    SIGT=0.45,      # -      log-SD of the threshold distribution
    HT=4.0,         # -      Hill coefficient of the threshold (sharpness)
    KAMHT=0.15,     # -      AMH-driven elevation of the FSH threshold
    AMHREF=2.5,     # ng/mL  reference AMH for the threshold term
    # --- follicle growth ----------------------------------------------------
    KGR=2.90,       # mm/d   maximal diameter growth rate
    WGR=0.35,       # -      growth-rate floor for small antral follicles
    DGR=9.0,        # mm     diameter at which growth accelerates
    NGR=3.0,
    DMAX=24.0,      # mm     asymptotic follicle diameter
    KLHA=0.75,      # -      LH contribution to drive in large follicles
    DLH=11.0,       # mm     diameter at which LHCGR appears on granulosa
    KLHARR=11.0,    # IU/L   LH tone that arrests SMALL follicles (PCOS)
    NLHARR=4.0,
    # --- granulosa mass -----------------------------------------------------
    KG=0.85,        # /d     granulosa recruitment rate toward capacity
    KGD=0.55,       # /d     granulosa loss when starved of FSH
    DPRE=20.0,      # mm     diameter at which granulosa capacity = 1
    DARO=12.0,      # mm     diameter at which aromatase competence is half-maximal
    NARO=3.0,
    # --- atresia ------------------------------------------------------------
    KATR=0.12,      # /d     atresia rate at zero FSH drive
    # --- ovulation / aspiration ---------------------------------------------
    TMAT=0.85,      # d      maturation clock centre (post-commitment)
    NMAT=6.0,       # -      maturation clock steepness
    TRUP=1.55,      # d      post-commitment clock for wall rupture
    RUPX=0.90,      # d      cumulative LHCGR signal needed to COMPLETE rupture
    NRUP=60.0,      # -      rupture switch steepness
    NRUPX=4.0,      # -      steepness of the completion requirement
    KRUP=6.0,       # /d     follicle collapse rate at full rupture
    KCOMMIT=0.06,   # d      LHCGR exposure that commits the oocyte
    WASP=0.05,      # d      duration of the aspiration itself (~1.2 h)
    WCOL=0.30,      # d      duration of post-aspiration follicular collapse
    KCOL=3.0,       # -      collapse strength of the aspiration pulse
    FGRAN=0.55,     # -      granulosa fraction removed by aspiration
    ETARET=0.88,    # -      oocyte recovery efficiency per aspirated follicle
    DASP=11.0,      # mm     50% aspiration-yield diameter
    NASP=8.0,
    DMII=13.0,      # mm     50% competence diameter
    NMII=5.0,
    # --- gonadotropin secretion --------------------------------------------
    KFS=53.1,       # IU/L/d FSH secretion rate constant
    KELF=4.0,       # /d     endogenous FSH elimination
    KI=120.0,       # pg/mL  inhibin-B feedback constant
    KE2F=100.0,     # pg/mL  oestradiol feedback constant
    FANTF=0.45,     # -      maximal antagonist suppression of FSH
    FANTL=0.90,     # -      maximal antagonist suppression of LH
    KRELLH=0.1245,  # /d     fractional LH release per unit GnRH drive
    SLH0=8000.0,    # IU     releasable pituitary LH pool (basal steady state)
    SLHMAX=12000.0, # IU     storage capacity of the gonadotroph pool
    KSYNLH=2988.0,  # IU/d   LH synthesis at an empty pool
    VLH=10.0,       # L      LH distribution volume
    KELLH=16.6,     # /d     LH elimination (t1/2 = 1 h)
    AMPS=55.0,      # -      GnRH drive amplitude of the endogenous surge
    AMPA=44.0,      # -      GnRH drive amplitude of an agonist trigger
    KR=1.20,        # /d     surge-readiness accumulation
    DSURG=12.0,     # mm     follicle size that licenses the positive feedback
    NDSURG=8.0,
    KROFF=0.45,     # /d     surge-readiness decay
    E2S=400.0,      # pg/mL  E2 set-point for positive feedback
    NE2S=8.0,       # -      steepness of the E2 positive-feedback switch
    NRS=12.0,       # -      steepness of the readiness-to-surge switch
    KP4S=2.0,       # ng/mL  progesterone block of the surge
    # Progesterone slows the GnRH pulse generator. In an unstimulated luteal
    # phase (P4 ~13) this halves LH and the corpus luteum still lives 14 d; in
    # an IVF luteal phase (P4 ~35 from 15 corpora lutea) it shuts LH down to
    # 0.6 IU/L, so the luteal mass destroys its own support. Only hCG — which
    # bypasses the pituitary entirely — escapes the loop.
    KP4LH=12.0,     # ng/mL  P4 that halves GnRH pulse drive
    NP4LH=2.0,
    # --- steroids / peptides ------------------------------------------------
    KE2G=6000.0,    # pg/mL/d per unit granulosa mass
    KE2CL=900.0,    # pg/mL/d per unit corpus-luteum mass
    KTHECA=0.40,    # IU/L   LH requirement of theca androgen supply
    KELE2=12.0,     # /d
    E2BASE=18.0,    # pg/mL
    KIB=100.7,      # pg/mL/d per unit granulosa mass
    WSMIB=0.10,     # -      small-antral contribution to inhibin B
    KELIB=3.0,      # /d
    KAMH=0.0838,    # ng/mL/d per small antral follicle
    WAMHB=0.55,     # -      pre-antral pool that keeps secreting AMH during COS
    KELAMH=0.60,    # /d
    DAMH=8.0,       # mm     AMH switches off above this diameter
    KP4G=0.075,     # ng/mL/d per unit granulosa mass (pre-trigger leak)
    KP4CL=18.1,     # ng/mL/d per unit corpus-luteum mass
    KSATCL=1.46,    # units  luteal mass at which P4 output per cell halves
    KELP4=0.55,     # /d
    KAP4=6.0,       # /d     absorption of vaginal/IM progesterone
    FP4=0.025,      # (ng/mL)/mg bioavailability factor of exogenous P4
    # --- corpus luteum ------------------------------------------------------
    KFCL=1.0,       # -      granulosa-to-luteal conversion efficiency
    KLYS0=0.10,     # /d     intrinsic luteolysis
    KLYSM=0.50,     # /d     LH-withdrawal luteolysis
    # --- LHCGR: ONE ligand pool, THREE response thresholds -----------------
    # The receptor is the same; the DOWNSTREAM processes are not, and each has
    # its own signal requirement. This is the whole trigger argument.
    KSUP=2.0,       # IU/L   steroidogenic / luteal support (basal pulses suffice)
    KLHV=40.0,      # IU/L   VEGF-inducing luteinisation signal (needs a surge)
    NLHV=3.0,
    KLHM=25.0,      # IU/L   meiotic-commitment signal (needs a full surge)
    NLHM=6.0,
    KDTH=0.45,      # -      fall of the FSH threshold as a follicle enlarges
    DTH=10.0,       # mm     midpoint of that fall
    # --- OHSS ---------------------------------------------------------------
    KV=1.0,         # /d     VEGF production per unit mass per unit occupancy
    WCL=1.5,        # -      luteal weighting of VEGF production
    KELV=3.0,       # /d
    KVP=3.05,       # units  VEGF for half-maximal permeability
    NVP=5.3,
    KP=1.40,        # /d     permeability induction
    KPOFF=0.70,     # /d     permeability resolution
    EMAXCAB=0.50,   # -      cabergoline effect on VEGFR2 signalling
    EC50CAB=0.55,   # units  (1.0 = steady state on 0.5 mg/d)
    LP=0.60,        # L/d    trans-capillary leak at PERM = 1
    KREAB=0.25,     # /d     peritoneal reabsorption
    KVPC=0.90,      # /d     plasma-volume compensation
    VP0=2.60,       # L
    RBCV=1.60,      # L      red-cell volume
    # --- drug PK ------------------------------------------------------------
    KAF=2.50, KELFX=0.400, VF=26.0, FBIOF=0.80,          # rFSH
    KACO=0.60, KELCO=0.241, FBIOCO=0.60, POTCO=10.0,     # corifollitropin
    KAANT=100.0, KELANT=1.28, VANT=20.0, FBIOANT=0.91,   # GnRH antagonist
    IC50ANT=0.44,                                        # ng/mL
    KAHCG=1.40, KELHCG=0.50, VHCG=20.0, FBIOHCG=0.70,    # hCG
    KAAGO=40.0, KELAGO=5.55, VAGO=30.0, FBIOAGO=0.90,    # triptorelin
    EC50AGO=0.05,                                        # ng/mL
    # An agonist OCCUPIES and then DOWN-REGULATES the GnRH receptor; a
    # competitive antagonist only occupies it. That pharmacological
    # difference is what makes the agonist trigger cost the luteal phase.
    KGRDOWN=3.0,    # /d     receptor down-regulation at full agonist occupancy
    KGRREC=0.20,    # /d     receptor resynthesis (t1/2 3.5 d)
    KACAB=3.0, KELCAB=0.256,                             # cabergoline (units)
    KELLET=0.37, EMAXLET=0.80, EC50LET=0.50,             # letrozole (units)
    # --- pregnancy hCG (late OHSS) -----------------------------------------
    KELPHCG=0.55, PH0=8.0, GPH=0.35,
    # --- laboratory / embryology chain -------------------------------------
    FHMG=0.20,      # -      hCG-equivalent bioactivity per IU of hMG LH activity
    FERT=0.72,      # 2PN per MII (ICSI)
    LBREU=0.52,     # live birth per euploid blastocyst transferred
)


def patient(AFC=12.0, T50=9.0, TONE=1.0, AGE=32.0, SIGT=0.45, KATR=None):
    """Patient-level covariates. AFC is the only lever on cohort size."""
    q = dict(AFC=AFC, T50=T50, TONE=TONE, AGE=AGE, SIGT=SIGT)
    if KATR is not None:
        q['KATR'] = KATR
    return q


def blast_rate(age):
    """Blastulation per 2PN — mild age effect."""
    return 0.55 - 0.011 * max(0.0, age - 33.0)


def euploid_frac(age):
    """Euploid fraction of blastocysts (Franasiak 2014 shape)."""
    return 0.85 / (1.0 + math.exp((age - 38.2) / 4.0))


# ----------------------------------------------------------------------------
# 2.  protocol container
# ----------------------------------------------------------------------------
class Protocol:
    def __init__(self, **kw):
        self.fsh_dose = kw.get('fsh_dose', 150.0)      # IU/d
        self.fsh_start = kw.get('fsh_start', 2.0)      # cycle day
        self.fsh_stop_day = kw.get('fsh_stop_day', None)   # coasting
        self.coast_days = kw.get('coast_days', 3.0)
        self.hmg_lh = kw.get('hmg_lh', 0.0)            # IU LH-activity/d
        self.cori = kw.get('cori', 0.0)                # ug, single injection
        self.ant_start = kw.get('ant_start', 6.0)      # cycle day (0 = none)
        self.ant_dose = kw.get('ant_dose', 0.25)       # mg/d
        self.trigger = kw.get('trigger', 'hcg')        # hcg|ago|dual|none
        self.hcg_dose = kw.get('hcg_dose', 10000.0)    # IU
        self.ago_dose = kw.get('ago_dose', 0.2)        # mg triptorelin
        self.cab = kw.get('cab', False)
        self.letro = kw.get('letro', False)
        self.opu_delay = kw.get('opu_delay', 1.5)      # d after trigger
        self.luteal_p4 = kw.get('luteal_p4', 0.0)      # mg/d vaginal P4
        self.fresh_transfer = kw.get('fresh_transfer', True)
        self.ppos = kw.get('ppos', False)              # progestin priming
        self.trig_ncrit = kw.get('trig_ncrit', 3.0)    # follicles >= 17 mm
        self.trig_dcrit = kw.get('trig_dcrit', 17.0)
        self.max_stim = kw.get('max_stim', 16.0)
        self.iv_fluid = kw.get('iv_fluid', 0.0)        # L/d supportive care
        self.name = kw.get('name', 'protocol')


# ----------------------------------------------------------------------------
# 3.  the right-hand side
# ----------------------------------------------------------------------------
def rhs(t, y, p, ctx):
    d = np.zeros(NST)
    D = np.clip(y[iD:iD + NF], 1e-6, None)
    G = np.clip(y[iG:iG + NF], 0.0, None)
    V = np.clip(y[iV:iV + NF], 0.0, 1.0)
    m = ctx['m']
    Ti = ctx['T']

    # ---- drug concentrations ----------------------------------------------
    cfsh_ex = y[FSHC] / p['VF'] + p['POTCO'] * y[CORIC] / p['VF']
    cant = max(y[ANTC], 0.0) / p['VANT']
    cago = max(y[AGOC], 0.0) / p['VAGO']
    chcg = max(y[HCG], 0.0) / p['VHCG']
    ccab = max(y[CAB], 0.0)
    clet = max(y[LET], 0.0)

    fshe = max(y[FSHE], 0.0)
    lh = max(y[LH], 0.0)
    e2 = max(y[E2], 0.0)
    inhb = max(y[INHB], 0.0)
    amh = max(y[AMH], 0.0)
    p4 = max(y[P4], 0.0)
    cl = max(y[CL], 0.0)
    fsh_tot = fshe + cfsh_ex

    # ---- GnRH receptor: competitive occupancy ----------------------------
    xa = cago / p['EC50AGO']
    xn = cant / p['IC50ANT']
    occ_ago = xa / (1.0 + xa + xn)
    occ_ant = xn / (1.0 + xa + xn)

    # ---- LHCGR: one ligand pool read by three kernels ---------------------
    lheq = lh + chcg + max(y[PHCG], 0.0)
    sup = lheq / (lheq + p['KSUP'])                 # support   (K = 2 IU/L)
    occv = hill(lheq, p['KLHV'], p['NLHV'])         # VEGF      (K = 12 IU/L)
    occm = hill(lheq, p['KLHM'], p['NLHM'])         # meiosis   (K = 25 IU/L)

    # ---- follicle dynamics ------------------------------------------------
    mclk = max(y[MCLK], 0.0)
    # Rupture of the follicle wall is an INTEGRAL of the LHCGR signal, whereas
    # oocyte maturation (MCLK below) is an autonomous clock started by the
    # signal's leading EDGE. That difference — not any drug-specific parameter —
    # is why a GnRH-agonist trigger matures oocytes without rupturing follicles.
    rup = hill(mclk, p['TRUP'], p['NRUP']) \
        * hill(max(y[ROV], 0.0), p['RUPX'], p['NRUPX'])
    # Aspiration is written as TWO consecutive rectangular windows: the oocyte
    # is counted while the follicle is still intact, and the follicle collapses
    # afterwards. Overlapping the two would count follicles that are already
    # emptying, which is exactly the defect this fixed.
    pulse = 0.0     # normalised counting window (integrates to 1)
    coll = 0.0      # collapse window
    if ctx['t_opu'] is not None:
        t0 = ctx['t_opu']
        if t0 <= t < t0 + p['WASP']:
            pulse = 1.0 / p['WASP']
        if t0 + p['WASP'] <= t < t0 + p['WASP'] + p['WCOL']:
            coll = 1.0 / p['WCOL']

    Mg = 0.0            # total functional granulosa mass (follicles)
    MgA = 0.0           # aromatase-competent granulosa mass
    n_asp = 0.0
    n_mii = 0.0
    for i in range(NF):
        # the FSH threshold is not a constant: it falls as the follicle grows
        # (rising FSHR density + intrafollicular IGF amplification). This is
        # the positive feedback that makes ONE follicle dominant.
        Teff = Ti[i] * (1.0 + p['KAMHT'] * (amh / p['AMHREF'] - 1.0)) \
            * (1.0 - p['KDTH'] * hill(D[i], p['DTH'], 4.0))
        Teff = max(Teff, 0.5)
        A = hill(fsh_tot, Teff, p['HT'])
        lhr = hill(D[i], p['DLH'], 6.0)                 # LHCGR acquisition
        drive = min(1.0, A + p['KLHA'] * lhr * sup)
        # PCOS: LH tone arrests follicles that have not yet acquired LHCGR
        arrest = 1.0 / (1.0 + (lh / p['KLHARR']) ** p['NLHARR'] * (1.0 - lhr))
        drive *= arrest
        gcap = min(1.0, (D[i] / p['DPRE']) ** 2)
        rup_i = rup * lhr           # only LHCGR-competent follicles rupture
        # small antral follicles grow slowly and accelerate as they enlarge
        # (2->5 mm takes ~10 d, 10->20 mm takes ~5 d)
        gfac = p['WGR'] + (1.0 - p['WGR']) * hill(D[i], p['DGR'], p['NGR'])
        d[iD + i] = (p['KGR'] * drive * gfac * (1.0 - D[i] / p['DMAX'])
                     - p['KRUP'] * rup_i * (D[i] - 3.0)
                     - p['KCOL'] * coll * hill(D[i], p['DASP'], p['NASP'])
                     * (D[i] - 3.0))
        d[iG + i] = (p['KG'] * drive * (gcap - G[i]) - p['KGD'] * (1.0 - drive) * G[i]
                     - p['KRUP'] * rup_i * G[i]
                     - p['FGRAN'] * p['KCOL'] * coll
                     * hill(D[i], p['DASP'], p['NASP']) * G[i])
        d[iV + i] = -ctx['KATR'] * (1.0 - drive) ** 2 * V[i]
        Mg += m * G[i] * V[i]
        # aromatase competence is acquired during antral maturation, so E2 is
        # NOT proportional to granulosa mass alone — a cohort of 8 mm follicles
        # makes very little oestradiol while one 20 mm follicle makes a lot.
        MgA += m * G[i] * V[i] * hill(D[i], p['DARO'], p['NARO'])
        gate = hill(D[i], p['DASP'], p['NASP'])
        n_asp += m * V[i] * gate
        n_mii += m * V[i] * gate * hill(D[i], p['DMII'], p['NMII'])

    # ---- drug PK ----------------------------------------------------------
    d[FSHDEP] = -p['KAF'] * y[FSHDEP]
    d[FSHC] = p['KAF'] * y[FSHDEP] - p['KELFX'] * y[FSHC]
    d[CORID] = -p['KACO'] * y[CORID]
    d[CORIC] = p['KACO'] * y[CORID] - p['KELCO'] * y[CORIC]
    d[ANTD] = -p['KAANT'] * y[ANTD]
    d[ANTC] = p['KAANT'] * y[ANTD] - p['KELANT'] * y[ANTC]
    d[HCGD] = -p['KAHCG'] * y[HCGD]
    d[HCG] = p['KAHCG'] * y[HCGD] - p['KELHCG'] * y[HCG]
    d[AGOD] = -p['KAAGO'] * y[AGOD]
    d[AGOC] = p['KAAGO'] * y[AGOD] - p['KELAGO'] * y[AGOC]
    d[CAB] = -p['KELCAB'] * y[CAB]
    d[LET] = -p['KELLET'] * y[LET]

    # ---- pituitary / gonadotropins ---------------------------------------
    p4blk = 1.0 / (1.0 + (p4 / p['KP4S']) ** 2)
    # The positive feedback is licensed by a PREOVULATORY-SIZE follicle, not by
    # oestradiol alone: this is why a stimulated cycle whose E2 passes 300 pg/mL
    # on day 5 does not surge on day 5, and why the premature LH rise tracks
    # follicle diameter rather than E2.
    dmaxf = float(np.max(D))
    d[RS] = p['KR'] * hill(e2, p['E2S'], p['NE2S']) \
        * hill(dmaxf, p['DSURG'], p['NDSURG']) * p4blk \
        * (1.0 - p['FANTL'] * occ_ant) * (1.0 - y[RS]) \
        - p['KROFF'] * y[RS]
    surge = p['AMPS'] * hill(max(y[RS], 0.0), 0.5, p['NRS']) * p4blk
    basal = ctx['TONE']
    gnrh_lh = (basal + surge) * (1.0 - p['FANTL'] * occ_ant) + p['AMPA'] * occ_ago
    gnrh_f = (basal + 0.25 * surge) * (1.0 - p['FANTF'] * occ_ant) \
        + 0.15 * p['AMPA'] * occ_ago

    p4fb = 1.0 / (1.0 + (p4 / p['KP4LH']) ** p['NP4LH'])
    gnrh_lh *= p4fb
    gnrh_f *= (0.5 + 0.5 * p4fb)
    grr = min(max(y[GRR], 0.0), 1.0)
    d[GRR] = p['KGRREC'] * (1.0 - grr) - p['KGRDOWN'] * occ_ago * grr
    gnrh_lh *= grr
    gnrh_f *= grr
    rel = p['KRELLH'] * gnrh_lh * max(y[SLH], 0.0)
    d[SLH] = p['KSYNLH'] * (1.0 - max(y[SLH], 0.0) / p['SLHMAX']) - rel
    d[LH] = rel / p['VLH'] - p['KELLH'] * lh

    finh = 1.0 / (1.0 + inhb / p['KI'] + e2 / p['KE2F'])
    d[FSHE] = p['KFS'] * gnrh_f * finh - p['KELF'] * fshe

    # ---- steroids ---------------------------------------------------------
    theca = lheq / (lheq + p['KTHECA'])
    aro = 1.0 - p['EMAXLET'] * clet / (clet + p['EC50LET'])
    d[E2] = (p['KE2G'] * MgA * theca + p['KE2CL'] * cl * theca) * aro \
        + p['E2BASE'] * p['KELE2'] - p['KELE2'] * e2
    nsm = sum(m * V[i] * (1.0 - hill(D[i], p['DAMH'], 4.0)) for i in range(NF))
    d[INHB] = p['KIB'] * (Mg + p['WSMIB'] * nsm) - p['KELIB'] * inhb
    small = sum(m * (1.0 - hill(D[i], p['DAMH'], 4.0)) for i in range(NF))
    d[AMH] = p['KAMH'] * (p['WAMHB'] * ctx['AFC'] + small) - p['KELAMH'] * amh
    d[P4D] = -p['KAP4'] * y[P4D]
    d[P4] = p['KP4G'] * Mg \
        + p['KP4CL'] * cl / (1.0 + cl / p['KSATCL']) \
        + p['FP4'] * p['KAP4'] * y[P4D] + ctx['ppos_p4'] \
        - p['KELP4'] * p4
    lys = p['KLYS0'] + p['KLYSM'] * (1.0 - sup) ** 2
    d[CL] = p['KFCL'] * (p['KRUP'] * rup
                         * sum(m * G[i] * V[i] * hill(D[i], p['DLH'], 6.0)
                               for i in range(NF))
                         + (1.0 - p['FGRAN']) * p['KCOL'] * coll
                         * sum(m * G[i] * V[i] * hill(D[i], p['DASP'], p['NASP'])
                               for i in range(NF))) - lys * cl

    # ---- the two read-outs: EDGE (commitment clock) and AREA (VEGF) ------
    d[EXPLH] = occm
    d[MCLK] = hill(y[EXPLH], p['KCOMMIT'], 8.0)   # autonomous once committed
    d[ROV] = occm                                  # signal-dependent integral
    d[VEGF] = p['KV'] * (Mg + p['WCL'] * cl) * occv - p['KELV'] * y[VEGF]
    cabe = p['EMAXCAB'] * ccab / (ccab + p['EC50CAB'])
    d[PERM] = p['KP'] * hill(y[VEGF], p['KVP'], p['NVP']) * (1.0 - cabe) \
        - p['KPOFF'] * y[PERM]
    leak = p['LP'] * max(y[PERM], 0.0) * max(y[VP], 0.1) / p['VP0']
    d[ASC] = leak - p['KREAB'] * max(y[ASC], 0.0)
    d[VP] = p['KVPC'] * (p['VP0'] - y[VP]) - leak + ctx['iv'](t)
    d[PHCG] = ctx['phcg_in'](t) - p['KELPHCG'] * max(y[PHCG], 0.0)

    # ---- accumulators -----------------------------------------------------
    d[OOC] = pulse * n_asp * p['ETARET'] * hill(y[EXPLH], p['KCOMMIT'], 4.0)
    d[MIIC] = pulse * n_mii * p['ETARET'] * hill(mclk, p['TMAT'], p['NMAT'])
    d[AUCE2] = e2 / 1000.0
    d[EXPV] = occv
    return d


# ----------------------------------------------------------------------------
# 4.  simulation driver (adaptive trigger day)
# ----------------------------------------------------------------------------
def simulate(pt, proto, p=P, tend=32.0, dt=0.02):
    q = dict(p)
    q.update({k: v for k, v in pt.items() if k in q})
    AFC, T50, SIGT = pt['AFC'], pt['T50'], pt.get('SIGT', q['SIGT'])
    m = AFC / NF
    Ti = T50 * np.exp(SIGT * ZQ)
    KATR = pt.get('KATR', q['KATR'])

    y = np.zeros(NST)
    y[iD:iD + NF] = 5.0
    y[iG:iG + NF] = 0.03
    y[iV:iV + NF] = 1.0
    y[SLH] = q['SLH0']
    y[LH] = 5.0
    y[FSHE] = 6.0
    y[E2] = 40.0
    y[INHB] = 45.0
    y[AMH] = q['KAMH'] * (q['WAMHB'] * AFC
                          + AFC * (1.0 - hill(5.0, q['DAMH'], 4.0))) / q['KELAMH']
    y[P4] = 0.4
    y[VP] = q['VP0']
    y[GRR] = 1.0

    ctx = dict(m=m, T=Ti, AFC=AFC, TONE=pt['TONE'], KATR=KATR, t_opu=None,
               ppos_p4=0.0, phcg_in=lambda t: 0.0,
               iv=lambda t: 0.0)

    # daily bolus schedule ---------------------------------------------------
    trig_t = None
    opu_t = None
    cancelled = False
    pending = []
    t = 0.0
    ts, ys = [0.0], [y.copy()]
    day = 0
    trig_record = {}
    while t < tend - 1e-9:
        # ---- morning boluses ---------------------------------------------
        cd = t                                      # cycle day (t = 0 -> CD1)
        stim_on = (cd >= proto.fsh_start - 1e-9 and trig_t is None and
                   (proto.fsh_stop_day is None or cd < proto.fsh_stop_day))
        if stim_on and proto.fsh_dose > 0 and proto.cori == 0:
            y[FSHDEP] += q['FBIOF'] * proto.fsh_dose
        if stim_on and proto.hmg_lh > 0:
            y[HCGD] += q['FBIOHCG'] * proto.hmg_lh * q['FHMG']  # hMG LH activity
        if proto.cori > 0 and abs(cd - proto.fsh_start) < 1e-9:
            y[CORID] += q['FBIOCO'] * proto.cori * 100.0
        if (proto.cori > 0 and trig_t is None and
                cd >= proto.fsh_start + 7.0 and proto.fsh_dose > 0):
            y[FSHDEP] += q['FBIOF'] * proto.fsh_dose      # ENGAGE day-8 rescue
        if (proto.ant_start > 0 and cd >= proto.ant_start - 1e-9
                and (opu_t is None or cd <= opu_t)):
            y[ANTD] += q['FBIOANT'] * proto.ant_dose * 1000.0
        ctx['ppos_p4'] = 4.0 * q['KELP4'] if (
            proto.ppos and proto.fsh_start - 1e-9 <= cd and trig_t is None) else 0.0
        if proto.cab and trig_t is not None and cd >= trig_t:
            y[CAB] = min(y[CAB] + 0.256, 1.2)
        if proto.letro and opu_t is not None and cd >= opu_t:
            y[LET] = min(y[LET] + 0.37, 1.3)
        if proto.luteal_p4 > 0 and opu_t is not None and cd >= opu_t:
            y[P4D] += proto.luteal_p4

        # ---- trigger decision (evening scan) -----------------------------
        if trig_t is None and cd >= proto.fsh_start + 3.0:
            nbig = sum(m * y[iV + i] for i in range(NF)
                       if y[iD + i] >= proto.trig_dcrit)
            prem = (y[LH] > 15.0 and y[P4] > 1.8 and cd > proto.fsh_start + 3)
            if prem:
                cancelled = True
            coast_done = (proto.fsh_stop_day is not None and
                          cd >= proto.fsh_stop_day + proto.coast_days)
            if (nbig >= proto.trig_ncrit or coast_done
                    or cd >= proto.fsh_start + proto.max_stim
                    or cancelled) and proto.trigger != 'none':
                trig_t = cd + 0.4
                opu_t = trig_t + proto.opu_delay
                ctx['t_opu'] = opu_t
                trig_record = dict(
                    trig_day=trig_t, stim_days=trig_t - proto.fsh_start,
                    nbig=nbig, cancelled=cancelled)
                hd = proto.hcg_dose if proto.trigger in ('hcg',) else (
                    1500.0 if proto.trigger == 'dual' else 0.0)
                ad = proto.ago_dose if proto.trigger in ('ago', 'dual') else 0.0
                # trigger injections are given at trig_t EXACTLY, mid-segment:
                # queue them as timed boluses rather than folding them into the
                # morning of the trigger day (which advanced the agonist by 9.6 h
                # and silently ovulated the dual-trigger arm before retrieval).
                if hd > 0:
                    pending.append((trig_t, HCGD, q['FBIOHCG'] * hd))
                if ad > 0:
                    pending.append((trig_t, AGOD, q['FBIOAGO'] * ad * 1000.0))
                if proto.fresh_transfer and proto.trigger != 'none':
                    def phcg_in(tt, t0=opu_t + 9.0):
                        if tt < t0:
                            return 0.0
                        return q['PH0'] * math.exp(min(q['GPH'] * (tt - t0), 8.0))
                    ctx['phcg_in'] = phcg_in
        if proto.iv_fluid > 0:
            ctx['iv'] = lambda tt: (proto.iv_fluid
                                    if (opu_t is not None and tt > opu_t + 3)
                                    else 0.0)

        # ---- integrate one day, breaking at any mid-day dose --------------
        breaks = sorted(set(e[0] for e in pending if t < e[0] < t + 1.0))
        for tb in breaks + [t + 1.0]:
            if tb - t > 1e-6:
                tev = np.clip(np.arange(t, tb + 1e-9, dt), t, tb)
                if tev.size < 2 or tev[-1] < tb - 1e-9:
                    tev = np.append(tev, tb)
                seg = solve_ivp(rhs, (t, tb), y, args=(q, ctx), method='LSODA',
                                rtol=1e-7, atol=1e-9, max_step=0.01,
                                dense_output=False, t_eval=tev)
                if not seg.success:
                    raise RuntimeError(seg.message)
                for k in range(1, seg.t.size):
                    ts.append(seg.t[k])
                    ys.append(seg.y[:, k].copy())
                y = seg.y[:, -1].copy()
                y[iV:iV + NF] = np.clip(y[iV:iV + NF], 0.0, 1.0)
                t = tb
            for (te, cmt, amt) in [e for e in pending if abs(e[0] - tb) < 1e-9]:
                y[cmt] += amt
            pending = [e for e in pending if e[0] > tb + 1e-9]
        day += 1

    ts = np.array(ts)
    Y = np.array(ys).T
    out = dict(t=ts, Y=Y, trig=trig_record, opu=opu_t, trig_t=trig_t,
               m=m, T=Ti, cancelled=cancelled, proto=proto, pt=pt)
    out.update(endpoints(out, q))
    return out


def at(out, tt, idx):
    return float(np.interp(tt, out['t'], out['Y'][idx]))


def endpoints(out, q):
    Y, ts = out['Y'], out['t']
    e = {}
    trig = out['trig_t']
    pre = ts < (trig if trig is not None else ts[-1])
    e['LH_pre_max'] = float(np.max(Y[LH][pre])) if pre.any() else 0.0
    e['prem_surge'] = bool(e['LH_pre_max'] > 12.0 and
                           (trig is None or trig > 4.0))
    if trig is not None:
        e['E2_trig'] = at(out, trig, E2)
        e['P4_trig'] = at(out, trig, P4)
        e['LH_trig'] = at(out, trig, LH)
        e['FSH_trig'] = at(out, trig, FSHE) + at(out, trig, FSHC) / q['VF'] \
            + q['POTCO'] * at(out, trig, CORIC) / q['VF']
        e['nfoll_11'] = sum(out['m'] * at(out, trig, iV + i)
                            for i in range(NF) if at(out, trig, iD + i) >= 11.0)
        e['nfoll_17'] = sum(out['m'] * at(out, trig, iV + i)
                            for i in range(NF) if at(out, trig, iD + i) >= 17.0)
        e['stim_days'] = out['trig'].get('stim_days', float('nan'))
    e['oocytes'] = float(Y[OOC, -1])
    e['MII'] = float(Y[MIIC, -1])
    e['LHexp_trig_opu'] = (at(out, out['opu'], EXPLH) - at(out, trig, EXPLH)) \
        if trig is not None else 0.0
    e['LHexp_total'] = float(Y[EXPLH, -1] - (at(out, trig, EXPLH) if trig else 0))
    e['AUCV'] = float(Y[EXPV, -1] - (at(out, trig, EXPV) if trig else 0))
    e['AUCV_all'] = float(Y[EXPV, -1])
    e['VEGF_max'] = float(np.max(Y[VEGF]))
    e['PERM_max'] = float(np.max(Y[PERM]))
    e['ASC_max'] = float(np.max(Y[ASC]))
    hct = q['RBCV'] / (q['RBCV'] + np.clip(Y[VP], 0.2, None)) * 100.0
    e['HCT_max'] = float(np.max(hct))
    e['E2_max'] = float(np.max(Y[E2]))
    e['P4_max'] = float(np.max(Y[P4]))
    if out['opu'] is not None:
        e['P4_lut7'] = at(out, out['opu'] + 7.0, P4)
        e['P4_lut5'] = at(out, out['opu'] + 5.0, P4)
    age = out['pt']['AGE']
    e['blast'] = e['MII'] * q['FERT'] * blast_rate(age)
    e['euploid'] = e['blast'] * euploid_frac(age)
    e['CLBR'] = 1.0 - (1.0 - q['LBREU']) ** e['euploid']
    e['OHSS'] = ohss_grade(e['HCT_max'], e['ASC_max'])
    e['MII_frac'] = e['MII'] / e['oocytes'] if e['oocytes'] > 0.05 else 0.0
    return e


def ohss_grade(hct, asc):
    """Golan-style grading read out of the fluid states, not asserted."""
    if hct >= 55.0 or asc >= 4.0:
        return 'critical'
    if hct >= 45.0 or asc >= 1.5:
        return 'severe'
    if hct >= 43.0 or asc >= 0.35:
        return 'moderate'
    if asc >= 0.10:
        return 'mild'
    return 'none'


# ----------------------------------------------------------------------------
# 5.  scenarios
# ----------------------------------------------------------------------------
NORMAL = patient(AFC=12.0, T50=9.0, TONE=1.0, AGE=32.0)
PCOS = patient(AFC=25.0, T50=9.0, TONE=1.8, AGE=31.0)
POOR = patient(AFC=5.0, T50=12.0, TONE=1.0, AGE=38.0)
OLDER = patient(AFC=8.0, T50=10.0, TONE=1.0, AGE=41.0)
EXTREME = patient(AFC=32.0, T50=9.0, TONE=1.7, AGE=29.0)

SCEN = [
    ('S0  natural cycle (no drugs)', NORMAL,
     Protocol(name='natural', fsh_dose=0, ant_start=0, trigger='none',
              fresh_transfer=False)),
    ('S1  antagonist + rFSH 150 + hCG 10000', NORMAL,
     Protocol(name='standard', luteal_p4=600)),
    ('S2  individualised low dose (rFSH 100)', NORMAL,
     Protocol(name='indiv', fsh_dose=100, luteal_p4=600)),
    ('S3  high responder, rFSH 150 + hCG', PCOS,
     Protocol(name='hi-hcg', luteal_p4=600)),
    ('S4  high responder, agonist trigger + freeze-all', PCOS,
     Protocol(name='hi-ago', trigger='ago', fresh_transfer=False)),
    ('S5  high responder, rFSH 75 + hCG (dose reduction)', PCOS,
     Protocol(name='hi-lowdose', fsh_dose=75, luteal_p4=600)),
    ('S6  high responder, hCG + cabergoline 0.5 mg', PCOS,
     Protocol(name='hi-cab', cab=True, luteal_p4=600)),
    ('S7  high responder, coasting (FSH stop CD8)', PCOS,
     Protocol(name='hi-coast', fsh_stop_day=8.0, luteal_p4=600)),
    ('S8  poor responder, rFSH 300', POOR,
     Protocol(name='poor', fsh_dose=300, luteal_p4=600)),
    ('S9  poor responder, rFSH 225 + hMG LH 75', POOR,
     Protocol(name='poor-lh', fsh_dose=225, hmg_lh=75, luteal_p4=600)),
    ('S10 corifollitropin alfa 150 ug single', NORMAL,
     Protocol(name='cori', cori=150.0, fsh_dose=150, luteal_p4=600)),
    ('S11 NO antagonist (premature LH surge)', NORMAL,
     Protocol(name='no-ant', ant_start=0.0, luteal_p4=600)),
    ('S12 dual trigger (agonist + hCG 1500)', PCOS,
     Protocol(name='dual', trigger='dual', luteal_p4=600)),
    ('S13 agonist trigger + FRESH transfer', PCOS,
     Protocol(name='ago-fresh', trigger='ago', fresh_transfer=True,
              luteal_p4=600)),
    ('S14 hCG trigger, retrieval delayed to 40 h', NORMAL,
     Protocol(name='late-opu', opu_delay=1.667, luteal_p4=600)),
    ('S15 age 41, AFC 8, rFSH 300', OLDER,
     Protocol(name='age41', fsh_dose=300, luteal_p4=600)),
    ('S16 high responder, hCG + freeze-all', PCOS,
     Protocol(name='hi-hcg-fa', fresh_transfer=False)),
    ('S17 high responder, PPOS + agonist trigger', PCOS,
     Protocol(name='ppos', ppos=True, ant_start=0.0, trigger='ago',
              fresh_transfer=False)),
    ('S18 AFC 32, hCG 10000 + fresh transfer', EXTREME,
     Protocol(name='xs-hcg', luteal_p4=600)),
    ('S19 AFC 32, agonist trigger + freeze-all', EXTREME,
     Protocol(name='xs-ago', trigger='ago', fresh_transfer=False)),
    ('S20 AFC 32, hCG + cabergoline + freeze-all', EXTREME,
     Protocol(name='xs-cab', cab=True, fresh_transfer=False)),
    ('S21 AFC 32, hCG 10000 + freeze-all', EXTREME,
     Protocol(name='xs-fa', fresh_transfer=False)),
]


def run_all(verbose=True):
    res = {}
    for label, pt, pr in SCEN:
        out = simulate(pt, pr, tend=34.0)
        res[label] = out
    return res


def fmt(x, n=1):
    return ('%.' + str(n) + 'f') % x if isinstance(x, float) else str(x)


def report(res, brief=False):
    print('=' * 108)
    print('CONTROLLED OVARIAN STIMULATION — QSP reference implementation '
          '(python/scipy, %d ODEs)' % NST)
    print('=' * 108)
    hdr = ('%-42s %5s %6s %6s %5s %5s %5s %6s %5s %5s %8s' %
           ('scenario', 'dstim', 'E2trg', 'P4trg', 'foll', 'OPU', 'MII',
            'AUCveg', 'Hct', 'asc', 'OHSS'))
    print(hdr)
    print('-' * 108)
    for label, _, _ in SCEN:
        o = res[label]
        print('%-42s %5s %6s %6s %5s %5s %5s %6s %5s %5s %8s' % (
            label,
            fmt(o.get('stim_days', float('nan')), 1),
            fmt(o.get('E2_trig', 0.0), 0),
            fmt(o.get('P4_trig', 0.0), 2),
            fmt(o.get('nfoll_11', 0.0), 1),
            fmt(o['oocytes'], 1),
            fmt(o['MII'], 1),
            fmt(o['AUCV'], 2),
            fmt(o['HCT_max'], 1),
            fmt(o['ASC_max'], 2),
            o['OHSS']))
    print('-' * 108)
    if brief:
        return
    print()
    print('embryology chain')
    print('%-42s %6s %6s %6s %6s %7s' %
          ('scenario', 'MII', '2PN', 'blast', 'eupl', 'CLBR%'))
    for label, pt, _ in SCEN:
        o = res[label]
        print('%-42s %6.1f %6.1f %6.1f %6.2f %7.1f' % (
            label, o['MII'], o['MII'] * P['FERT'], o['blast'],
            o['euploid'], 100 * o['CLBR']))


# ----------------------------------------------------------------------------
# 6.  sweeps — the quantitative content of the argument
# ----------------------------------------------------------------------------
def sweep_dose():
    print()
    print('A. FSH DOSE SWEEP — does dose separate yield from OHSS?')
    print('   (same patient, same antagonist protocol, same hCG 10000 trigger)')
    print('%-22s %5s %6s %6s %6s %6s %6s %6s %9s' % (
        'patient / dose', 'dstim', 'E2trg', 'OPU', 'MII', 'eupl', 'Hct',
        'ascL', 'OHSS'))
    for label, pt in (('AFC 12', NORMAL), ('AFC 25', PCOS)):
        for dose in (75, 112, 150, 225, 300, 450):
            o = simulate(pt, Protocol(fsh_dose=dose, luteal_p4=600), tend=30.0)
            print('%-22s %5.1f %6.0f %6.1f %6.1f %6.2f %6.1f %6.2f %9s' % (
                '%s / %d IU' % (label, dose), o.get('stim_days', 0),
                o.get('E2_trig', 0), o['oocytes'], o['MII'], o['euploid'],
                o['HCT_max'], o['ASC_max'], o['OHSS']))


def sweep_trigger():
    print()
    print('B. TRIGGER SWEEP in the AFC-32 patient — ONE ligand, THREE kernels')
    print('   AUCsup = luteal support (K 2 IU/L) | AUCveg = VEGF (K 40, n 3)')
    print('   AUCmat = meiosis (K 25, n 6) | edge = h to commitment')
    print('%-26s %7s %7s %7s %6s %6s %6s %6s %6s %8s' % (
        'trigger', 'peak', 'AUCveg', 'AUCmat', 'OPU', 'MII', 'MII%',
        'Hct', 'ascL', 'P4d7'))
    trigs = [('hCG 10000 IU', dict(trigger='hcg', hcg_dose=10000)),
             ('hCG 5000 IU', dict(trigger='hcg', hcg_dose=5000)),
             ('hCG 2500 IU', dict(trigger='hcg', hcg_dose=2500)),
             ('hCG 1500 IU', dict(trigger='hcg', hcg_dose=1500)),
             ('dual (ago + hCG 1500)', dict(trigger='dual')),
             ('agonist 0.2 mg alone', dict(trigger='ago'))]
    for label, kw in trigs:
        o = simulate(EXTREME, Protocol(fresh_transfer=False, **kw), tend=30.0)
        Y, t = o['Y'], o['t']
        tr = o['trig_t']
        post = t >= tr
        peak = float(np.max(Y[LH][post] + Y[HCG][post] / P['VHCG']))
        print('%-26s %7.0f %7.2f %7.2f %6.1f %6.1f %6.2f %6.1f %6.2f %8.1f' % (
            label, peak, o['AUCV'], o['LHexp_total'], o['oocytes'], o['MII'],
            o['MII_frac'], o['HCT_max'], o['ASC_max'],
            o.get('P4_lut7', float('nan'))))


def sweep_opu():
    print()
    print('C. RETRIEVAL-TIMING SWEEP (hCG trigger, AFC 12) — the 34-38 h window')
    print('%-14s %8s %8s %8s %8s' % ('trigger->OPU', 'oocytes', 'MII',
                                     'MII frac', 'ruptured'))
    for h in (28, 32, 34, 36, 38, 40, 44):
        o = simulate(NORMAL, Protocol(opu_delay=h / 24.0, luteal_p4=600),
                     tend=30.0)
        ref = 8.0
        print('%-14s %8.1f %8.1f %8.2f %8s' % (
            '%d h' % h, o['oocytes'], o['MII'], o['MII_frac'],
            '%.0f%%' % (100 * max(0.0, 1 - o['oocytes'] / ref))))


def sweep_afc():
    print()
    print('D. OVARIAN-RESERVE SWEEP (rFSH 150, hCG, age 32) — yield and CLBR')
    print('%-8s %6s %7s %7s %7s %7s %7s %8s' % (
        'AFC', 'AMH', 'OPU', 'MII', 'blast', 'eupl', 'CLBR%', 'OHSS'))
    for afc in (4, 6, 9, 12, 16, 20, 25, 32):
        pt = patient(AFC=float(afc), T50=9.0, TONE=1.0, AGE=32.0)
        o = simulate(pt, Protocol(luteal_p4=600), tend=30.0)
        amh = float(o['Y'][AMH, 0])
        print('%-8d %6.2f %7.1f %7.1f %7.1f %7.2f %7.1f %8s' % (
            afc, amh, o['oocytes'], o['MII'], o['blast'], o['euploid'],
            100 * o['CLBR'], o['OHSS']))


def sweep_antagonist():
    print()
    print('E. ANTAGONIST START DAY (AFC 12, rFSH 150) — surge protection')
    print('%-14s %7s %8s %8s %8s %8s' % ('antagonist', 'LHmax', 'surge?',
                                         'dstim', 'oocytes', 'MII'))
    for st in (0.0, 4.0, 5.0, 6.0, 7.0, 8.0):
        o = simulate(NORMAL, Protocol(ant_start=st, luteal_p4=600), tend=30.0)
        print('%-14s %7.1f %8s %8.1f %8.1f %8.1f' % (
            'none' if st == 0 else 'CD%d' % st, o['LH_pre_max'],
            'YES' if o['prem_surge'] else 'no', o.get('stim_days', 0),
            o['oocytes'], o['MII']))


def sweep_age():
    print()
    print('F. AGE SWEEP at fixed AFC 12 — oocytes are not the limiting factor')
    print('%-8s %7s %7s %7s %7s %8s' % ('age', 'OPU', 'MII', 'blast', 'eupl',
                                        'CLBR%'))
    for age in (28, 32, 36, 39, 42, 44):
        pt = patient(AFC=12.0, T50=9.0, TONE=1.0, AGE=float(age))
        o = simulate(pt, Protocol(luteal_p4=600), tend=30.0)
        print('%-8d %7.1f %7.1f %7.1f %7.2f %8.1f' % (
            age, o['oocytes'], o['MII'], o['blast'], o['euploid'],
            100 * o['CLBR']))


if __name__ == '__main__':
    brief = '--brief' in sys.argv
    r = run_all()
    report(r, brief=brief)
    if '--no-sweeps' not in sys.argv:
        sweep_dose()
        sweep_trigger()
        sweep_opu()
        sweep_afc()
        sweep_antagonist()
        sweep_age()
    print()
    print('done.')
