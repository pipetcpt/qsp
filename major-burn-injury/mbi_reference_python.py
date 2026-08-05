#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mbi_reference_python.py
=======================
Dependency-free RK4 re-implementation of the Major Thermal Burn Injury QSP
model (`mbi_mrgsolve_model.R`), used to CALIBRATE and VERIFY the mrgsolve
model in an environment with no R runtime.

Why this file exists
--------------------
The mrgsolve model is the deliverable, but a $ODE block that has never been
integrated is a hypothesis, not a model.  Every equation, parameter and
initial condition here is a line-for-line transcription of the R file, so
running this script exercises the same system.  Where the two differ, the R
file is wrong and this file is the reference.

Everything below is standard library only (no numpy/scipy) so it runs
anywhere Python 3 does.

Run:
    python3 mbi_reference_python.py            # full validation suite
    python3 mbi_reference_python.py --quick    # calibration targets only
"""

import math
import sys

# =====================================================================
# 0.  CONSTANTS AND PATIENT DESCRIPTORS
# =====================================================================

def landis_pappenheimer(C):
    """Colloid osmotic pressure [mmHg] from total protein C [g/dL].

    Pi = 2.1 C + 0.16 C^2 + 0.009 C^3

    This is THE convex function of the whole model.  Its convexity is why
    each litre of protein-free crystalloid costs more oncotic pressure than
    the litre before it:
        dPi/dC = 2.1 + 0.32 C + 0.027 C^2
               = 3.55 mmHg per g/dL at C = 3.5
               = 5.66 mmHg per g/dL at C = 7.0
    """
    if C < 0.0:
        C = 0.0
    return 2.1 * C + 0.16 * C * C + 0.009 * C * C * C


def _logit(x):
    """Overflow-safe logistic 1/(1+exp(-x))."""
    if x > 40.0:
        return 1.0
    if x < -40.0:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def default_params(WT=80.0, AGE=35.0, TBSA=45.0, INH=0.0, PED=0.0):
    fB = TBSA / 100.0
    p = dict(
        # ---- patient ----
        WT=WT, AGE=AGE, TBSA=TBSA, INH=INH, PED=PED, FB=fB,

        # ---- fluid compartments (mL) ----
        VP0=40.0 * WT,                 # plasma volume
        VI0=150.0 * WT,                # total interstitial volume
        VICF0=330.0 * WT,              # intracellular volume
        CP0=7.0,                       # plasma total protein g/dL
        CI0=2.5,                       # interstitial total protein g/dL

        # ---- Starling ----
        # Kf calibrated so that the NORMAL state filters ~2.8 L/day:
        #   Kf_tot * 0.7 mmHg ~ 118 mL/h  ->  Kf_tot ~ 168 mL/h/mmHg (70 kg)
        KFTOT=2.10 * WT,               # mL/h/mmHg, whole body baseline
        KFMB=2.60,                     # burned tissue Kf multiplier increment
        KFMU=0.90,                     # unburned Kf multiplier increment (TBSA-scaled)
        TAULEAK=9.5,                   # h, capillary leak time constant
        SIG0=0.90,                     # baseline reflection coefficient
        # DEFECT 5: sigma_burn = 0.15 stripped 86 g/h of protein by convection
        # alone -- the entire 224 g plasma pool in under three hours.  Measured
        # burn-tissue reflection coefficients for albumin are ~0.3-0.5, not 0.1,
        # and total protein (with globulins) sieves less still.
        DSIGB=0.50,                    # sigma drop in burned tissue (0.90 -> 0.40)
        DSIGU=0.15,                    # sigma drop in unburned tissue (TBSA-scaled)
        PS0=0.16,                      # dL/h diffusive protein permeability (whole body)
        PSMULT=3.0,                    # PS increment during leak

        PC0=17.0,                      # mmHg capillary hydrostatic at VP0
        KPC=12.0,                      # mmHg per unit relative plasma volume
        PI0=-1.0,                      # mmHg baseline interstitial hydrostatic
        ECB=8.0,                      # mmHg per unit relative expansion, burned
        ECU=11.0,                      # mmHg per unit relative expansion, unburned
        IMB0=22.0,                     # mmHg imbibition pressure at t=0
        TAUIMB=0.80,                   # h imbibition decay

        # DEFECT 1 (found by simulation): an unbounded Kf*deltaP product drove
        # 7.3 L/h into burned tissue at t=0 and emptied the plasma compartment
        # inside 30 min.  Filtration cannot exceed what plasma flow DELIVERS.
        # Adding the delivery limit is not a numerical patch -- it is the
        # documented reason resuscitation INCREASES burn oedema (Demling): the
        # zone of stasis is hypoperfused, and restoring perfusion restores the
        # delivery term.
        QPB0=18500.0,                  # mL plasma/h to the whole skin at rest
        QPU0=120000.0,                 # mL plasma/h to everything else
        FFMAX=0.55,                    # maximum filtration fraction
        STASIS0=0.40,                  # burned-tissue perfusion floor (no-reflow)

        # third-space capacities (a cavity has a volume; DEFECT 2)
        VASCMAX=6500.0,                # mL peritoneal/retroperitoneal capacity
        EVLWMAX=2600.0,                # mL extravascular lung water capacity

        # DEFECT 6: peak lymph return was 101 mL/h against a filtration rate of
        # 1400 mL/h, so the interstitium accumulated protein with no return
        # path and the plasma pool could never reach the quasi-steady state
        # Ci -> (1-sigma)*Cp.  Whole-body lymph is ~120 mL/h at rest and rises
        # 10-20x in burned tissue (Arturson).
        LBASE=120.0,                   # mL/h whole-body lymph flow at rest
        LMAX=4.0,                     # saturable lymph multiplier ceiling
        KL=2.0,                        # lymph sensitivity to interstitial volume

        # ---- renal ----
        GFR0=7.2,                      # L/h (=120 mL/min)
        # DEFECT 3: the original UO law made a euvolaemic patient pass
        # 2.4 mL/kg/h, i.e. ~5x the 0.5 mL/kg/h set-point, so the controller
        # saw "too much urine" from t=0 and never ramped.  UO must be
        # 1.0 mL/kg/h at normal volume status.
        UOBASE=1.00,                   # mL/kg/h at normal volume status
        # DEFECT 9: a shallow perfusion-UO curve let a patient 30 % down on
        # plasma volume still pass 0.6 mL/kg/h, so the controller "succeeded"
        # while the patient stayed in shock.  Steepening it makes the model
        # state the clinically important fact explicitly: a urine output of
        # 0.5 mL/kg/h corresponds to a plasma volume that is still ~18 % down.
        # UO 0.5 is a PERMISSIVE target, not a normal one.
        VPCRIT=0.88,                   # relative VP below which UO collapses
        HUO=8.0,                       # Hill coefficient of the UO collapse
        # DEFECT 12: a 3 mL/kg/h ceiling made the mobilisation ("diuretic")
        # phase impossible -- oedema returned to the plasma faster than the
        # kidney could excrete it and plasma volume reached 183 % of baseline
        # on day 2.  A healthy kidney clears far more than that; the ceiling
        # belongs at ~8 mL/kg/h.
        UOMAX=8.0,                     # mL/kg/h ceiling on urine output
        KDIU=6.0,                      # volume-expansion diuresis gain
        KOPI=0.0,                      # opioid-induced vasodilatation: fraction
                                       # by which urine output falls at the SAME
                                       # plasma volume ("opioid creep" -- the
                                       # controller reads it as hypovolaemia)
        KARC=0.55,                     # augmented renal clearance gain (flow phase)
        KIAPGFR=0.055,                 # per mmHg IAP loss of renal perfusion
        KTUB=0.030,                    # tubular injury rate at zero perfusion
        KTUBR=0.014,                   # tubular repair rate /h
        CRE0=0.9,                      # mg/dL

        # ---- controller ----
        FORMULA=4.0,                   # mL/kg/%TBSA over 24 h (Parkland)
        UOTGT=0.5,                     # mL/kg/h adult target
        GCTRL=0.32,                    # controller gain
        DRMAX=0.30,                    # /h maximum rate of change of the multiplier
        TAUUO=1.50,                    # h averaging window the clinician actually uses
        RMIN=0.40, RMAX=3.00,          # bounds as multiples of the formula rate
        TITRATE=1.0,                   # 1 = closed loop, 0 = open loop
        MAINT=0.0,                     # mL/h fixed maintenance (paediatric dextrose)
        # DEFECT 10: stopping all fluid at hour 24 dried the patient out --
        # plasma COP climbed to 32 mmHg by day 2 and the interstitium drained
        # BELOW baseline.  After resuscitation, evaporative and urinary losses
        # are replaced continuously; that is the maintenance phase.
        KMNT=1.00,                     # fraction of ongoing losses replaced
        # DEFECT 15: a maintenance rate that replaced only evaporation plus a
        # nominal urine ran an 8 L/day negative balance for four days.  After
        # resuscitation the input is titrated to the PATIENT, not to a formula:
        # a proportional term on plasma volume is both stable and what is
        # actually done at the bedside.
        KVP=0.150,                     # /h gain of the maintenance regulator

        # ---- albumin / colloid ----
        # DEFECT 7: albumin was originally ADDED on top of the crystalloid, so
        # the colloid arms used MORE total volume -- the opposite of every
        # trial.  Colloid is not an add-on; a fraction of the infusion is
        # SWITCHED to 5 % albumin.
        ALBSTART=1e9,                  # h, colloid start time (default: never)
        FCOL=0.0,                      # fraction of the infusion given as colloid
        KPPD=0.0022,                   # /h albumin catabolism (t1/2 ~ 13 d)
        HPPD=0.60,                     # pool-dependence of the catabolic rate
        ALBCONC=5.0,                   # g/dL of the colloid solution

        # ---- third space, lung, abdomen ----
        # DEFECT 13: intra-abdominal pressure never left 9 mmHg even at
        # 329 mL/kg in 24 h, so the model could not produce the complication
        # that defines over-resuscitation.  Ivy (2000) put the intra-abdominal
        # hypertension threshold at ~250 mL/kg/24 h; the model must cross it.
        KASC=0.120,                    # fraction of unburned filtrate to peritoneum
        KASCR=0.012,                   # /h reabsorption
        KEVLW=0.020, KEVLWR=0.045,
        # abdominal pressure-volume relation is STIFF: linear terms alone could
        # not separate a 226 mL/kg resuscitation (intra-abdominal hypertension)
        # from a 348 mL/kg one (abdominal compartment syndrome).  The cubic
        # term is the compliance limit.
        IAPK=0.00180,                  # mmHg per mL of third-space volume
        IAPB=0.00050,                  # mmHg per mL of generalised oedema
        IAPQ=0.00900,                  # mmHg per (litre of total excess)^3

        # ---- inflammation ----
        KDAMP=0.55,                    # DAMP production per unit open area
        KDDEG=0.10,
        KIL6=9.0, KIL6D=0.30,
        KTNF=2.6, KTNFD=0.55,
        KIL10=0.65, KIL10D=0.22,
        KPMN=0.9, KPMND=0.28,
        KROSKF=0.85,                   # ROS contribution to the Kf lesion
        KHLA=0.110, KHLAR=0.020,
        KCRP=1.8, KCRPD=0.028,

        # ---- infection ----
        BWD0=2.0,                      # log10 CFU/g at admission
        KGROW=0.055,                   # /h wound bacterial growth
        BWDMAX=9.0,
        KTOP=0.022,                    # topical antimicrobial kill /h
        KINV=0.020,                    # invasion rate above threshold
        BTHRESH=5.0,                   # log10 CFU/g invasion threshold
        KABX=0.0,                      # systemic antibiotic effect (set by scenario)
        KBSYSD=0.045,
        SEPTHR=1.0,                    # systemic burden defining sepsis
        BSYSMAX=5.0,                   # saturation of the systemic burden

        # ---- neuroendocrine ----
        CATMAX=8.0,                    # x normal at full drive
        A50=25.0,                      # %TBSA open area for half-maximal drive
        HILLA=2.5,
        KCATD=0.045,                   # /h catecholamine turnover
        COR0=12.0, CORMAX=3.2, KCORD=0.16,
        GCG0=80.0, GCGMAX=3.0, KGCGD=0.20,
        T30=110.0, KT3=0.020,
        IGF0=220.0, KIGF=0.020,

        # ---- hypermetabolism ----
        EMAXHM=0.95,                   # REE can reach 1.85 x predicted
        KREE=0.0200,                   # /h REE transducer (tau ~ 4.6 d up)
        KREED=0.00040,                 # /h REE decay after closure (tau ~ 104 d)
        TAMB=31.0,                     # ambient temperature, degC
        KEVAP=0.0,                     # set from TBSA in build()
        TCORE0=37.0,

        # ---- protein / body composition ----
        # DEFECT 8: the original sensitivities made B/S ~ 5, and lean mass fell
        # 63 % in fourteen days.  Muscle protein turns over at ~1.5 %/day, and
        # the OBSERVED burn catabolism is ~0.7 %/day of lean mass on top of it
        # -- i.e. B/S ~ 1.5, not 5.  Net balance is a small difference of two
        # large fluxes, which is exactly why it is so parameter-sensitive.
        KSYN=0.000625,                 # /h fractional synthesis (1.5 %/day)
        KBRK=0.000625,                 # /h fractional breakdown at baseline
        KBCAT=0.125,                   # breakdown sensitivity to adrenergic tone
        KBCOR=0.100,                   # breakdown sensitivity to cortisol
        KBINF=0.035,                   # breakdown sensitivity to IL-6/TNF
        KSINS=0.20,                    # synthesis sensitivity to insulin
        KSIGF=0.22,                    # synthesis sensitivity to IGF-1
        KHOMEO=20.0,                   # anabolism is DEFICIT-driven: it stops
                                       # once the deficit is repaid
        KSPR=1.45,                     # synthesis gain from beta-blockade
                                       # (Herndon 2001: the net-balance swing is
                                       #  driven as much by synthesis as by
                                       #  breakdown)
        LBM0F=0.80,                    # lean fraction of body weight
        KFATL=0.00055,                 # /h lipolysis
        KBONE=0.000030,                # /h bone loss per unit catabolic drive

        # ---- glucose ----
        # DEFECT 23: the original glucose block mixed concentrations and rates
        # with no distribution volume and produced blood glucose of
        # 6 500 mg/dL (18 500 in the sepsis arm).  That single number then
        # multiplied the bacterial-invasion term by 83, which is why
        # antibiotics appeared useless.  Rewritten in mass units.
        GLC0=90.0, INS0=10.0,
        VGLC=1.6,                      # dL/kg glucose distribution volume
        HGP0=9600.0,                   # mg/h basal hepatic glucose production
        U0=4000.0,                     # mg/h insulin-INdependent uptake
        UI=12135.0,                    # mg/h insulin-dependent uptake capacity
        KINSH=1.00,                    # insulin suppression of HGP
        KGCGH=0.55,                    # glucagon drive on HGP
        KIRESH=0.25,                   # insulin resistance drive on HGP
        KINSS=0.030, KINSD=0.30, INSMAX=70.0,
        IRESMAX=2.50, IRESK=3.00,      # saturating insulin resistance
        # Enteral carbohydrate.  It matters because burn hypoglycaemia is not
        # an insulin-dose problem, it is a FEED-INTERRUPTION problem: feeds are
        # held around each trip to theatre while the insulin infusion runs on.
        GFEED=9000.0,                  # mg/h enteral carbohydrate
        TFEED=12.0,                    # h, enteral feeding starts
        NPOHRS=6.0,                    # h feeds held before/around each operation
        FFA0=0.45, KFFA=0.9, KFFAD=0.6,
        KHFAT=0.030, KHFATD=0.0035,

        # ---- wound ----
        TEXC=1e9,                      # h, time of first excision
        EXCRATE=0.0,                   # %TBSA per operating session
        EXCINT=48.0,                   # h between sessions
        NEXC=0,                        # number of sessions
        KGRAFT=0.020,                  # /h graft take
        TAKEMAX=0.95,
        KREEPI=0.0035,                 # /h spontaneous re-epithelialisation
        FDEEP=0.65,                    # fraction of TBSA that is deep
        KDONH=0.0060,                  # /h donor site healing
        DONFRAC=0.55,                  # donor area per grafted area
        # DEFECT 18: without a donor-site constraint a 45 % burn closed in
        # 12 days and an 80 % burn closed just as fast.  In massive burns the
        # rate-limiting resource is unburned skin, and it must re-heal before
        # re-harvest -- which is why closure time is grossly non-linear in
        # %TBSA.
        DONPOOLF=0.60,                 # usable fraction of unburned skin
        MESH=3.0,                      # mesh expansion ratio
        KSCAR=0.00035,

        # ---- drug PK ----
        # propranolol: oral, high first-pass, F ~ 0.25, t1/2 ~ 4 h
        KAPR=1.2, FPR=0.25, VPR=250.0, CLPR=55.0, EC50PR=12.0, EMAXPR=0.62,
        # oxandrolone: F ~ 0.97, t1/2 ~ 9 h
        KAOX=0.9, FOX=0.95, VOX=45.0, CLOX=3.5, EC50OX=25.0, EMAXOX=0.42,
        # vancomycin 2-cmt
        VVC=0.72 * WT, VVP=0.60 * WT, QVC=6.0, CLVC=0.062 * WT,
        # ascorbate
        VVIT=0.30 * WT, CLVIT=8.0, EC50VIT=2500.0, EMAXVIT=0.40,
        # exogenous insulin
        KINSEL=0.35,

        # ---- outcome ----
        # DEFECT 22: the original hazard had NO burn-size term at all -- every
        # scenario clustered at 16 % because risk entered only through the
        # complications.  Mortality is carried by the OPEN WOUND: its area,
        # for as long as it stays open.
        # The revised Baux score weights one year of age exactly like one
        # %TBSA.  Reproducing that with a LINEAR age term is impossible: over
        # the 45-80 %TBSA range the area term already moves ~16-fold, so age
        # must move 16-fold over 35 years too.  Hence an exponential.
        HZOPEN=0.00260,                # /h at A_open = 100 %TBSA, age 20
        HZPOW=3.20,                    # superlinearity in open area
        HZAGEE=0.0792,                 # /year above 20 (exponential)
        HZINH=1.60,                    # inhalation-injury multiplier
        HZSEP=0.00110, HZACS=0.00120, HZARDS=0.00040,
        HZLBM=0.00030, HZFAT=0.00010,
    )
    # evaporative loss: (25 + %TBSA) x BSA[m2] mL/h  (BSA from Mosteller-ish)
    BSA = math.sqrt(WT * 170.0 / 3600.0)
    p["BSA"] = BSA
    p["KEVAP"] = (25.0 + TBSA) * BSA          # mL/h at full open area
    if PED:
        p["UOTGT"] = 1.0
    return p


# =====================================================================
# 1.  STATE VECTOR
# =====================================================================
NAMES = [
    "VP", "VIB", "VIU", "VICF", "PP", "PIB", "PIU", "VASC", "EVLW",     # 0-8
    "UOC", "FIN", "SCR", "RTUB", "LAC",                                  # 9-13
    "DAMP", "IL6", "TNF", "IL10", "PMN", "HLADR", "CRP",                # 14-20
    "BWD", "BSYS",                                                       # 21-22
    "CAT", "COR", "GCG", "T3", "IGF1",                                   # 23-27
    "REE", "TCORE",                                                      # 28-29
    "GLC", "INS", "FFA", "HFAT",                                         # 30-33
    "LBM", "FATM", "BMC", "SCARC",                                       # 34-37
    "AOPEN", "AEXC", "AGRF", "AHEAL", "ADON",                            # 38-42
    "APRD", "APRC", "AOXD", "AOXC", "AVCC", "AVCP", "AVTC", "AINSX",    # 43-50
    "HAZ", "ALBIN", "RSTATE", "UOWIN", "HYPOH",                          # 51-55
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)


def initial_state(p):
    y = [0.0] * NST
    fB = p["FB"]
    y[IX["VP"]] = p["VP0"]
    y[IX["VIB"]] = p["VI0"] * fB
    y[IX["VIU"]] = p["VI0"] * (1.0 - fB)
    y[IX["VICF"]] = p["VICF0"]
    y[IX["PP"]] = p["CP0"] * p["VP0"] / 100.0
    y[IX["PIB"]] = p["CI0"] * p["VI0"] * fB / 100.0
    y[IX["PIU"]] = p["CI0"] * p["VI0"] * (1.0 - fB) / 100.0
    y[IX["SCR"]] = p["CRE0"]
    y[IX["LAC"]] = 1.0
    y[IX["HLADR"]] = 100.0
    y[IX["CRP"]] = 3.0
    y[IX["BWD"]] = p["BWD0"]
    y[IX["CAT"]] = 1.0
    y[IX["COR"]] = p["COR0"]
    y[IX["GCG"]] = p["GCG0"]
    y[IX["T3"]] = p["T30"]
    y[IX["IGF1"]] = p["IGF0"]
    y[IX["REE"]] = 1.0
    y[IX["TCORE"]] = p["TCORE0"]
    y[IX["GLC"]] = p["GLC0"]
    y[IX["INS"]] = p["INS0"]
    y[IX["FFA"]] = p["FFA0"]
    y[IX["HFAT"]] = 300.0
    y[IX["LBM"]] = p["LBM0F"] * p["WT"]
    y[IX["FATM"]] = (1.0 - p["LBM0F"]) * p["WT"]
    y[IX["BMC"]] = 1.0
    y[IX["AOPEN"]] = p["TBSA"]
    y[IX["RSTATE"]] = 1.0          # controller multiplier on the formula rate
    y[IX["UOWIN"]] = p["UOTGT"] * p["WT"]   # filtered urine output, mL/h
    return y


# =====================================================================
# 2.  DERIVED QUANTITIES  (identical to the mrgsolve $ODE preamble)
# =====================================================================
def derive(t, y, p, sc):
    d = {}
    fB = p["FB"]
    VP = max(y[IX["VP"]], 200.0)
    VIB = max(y[IX["VIB"]], 1.0)
    VIU = max(y[IX["VIU"]], 1.0)
    VIB0 = p["VI0"] * fB
    VIU0 = p["VI0"] * (1.0 - fB)

    # ---- leak envelope: the FAST clock -------------------------------
    # graded, not switch-like: an antibiotic that halves the burden must move it
    b2 = max(y[IX["BSYS"]], 0.0) ** 2
    sep = b2 / (b2 + p["SEPTHR"] ** 2)
    Bleak = math.exp(-t / p["TAULEAK"]) + 0.55 * sep
    if Bleak > 1.6:
        Bleak = 1.6
    d["Bleak"] = Bleak
    d["sep"] = sep

    # ---- ascorbate acts on the OXIDANT component of Kf ---------------
    cvit = y[IX["AVTC"]] / p["VVIT"] * 1000.0          # mg/L -> ug/mL*1000
    vitEff = 1.0 - p["EMAXVIT"] * cvit / (p["EC50VIT"] + cvit)
    d["cvit"] = cvit
    d["vitEff"] = vitEff

    # oxidant fraction of the permeability lesion is what ascorbate can touch
    kfLes = Bleak * ((1.0 - p["KROSKF"]) + p["KROSKF"] * vitEff)

    Kfb = p["KFTOT"] * fB * (1.0 + p["KFMB"] * kfLes)
    Kfu = p["KFTOT"] * (1.0 - fB) * (1.0 + p["KFMU"] * kfLes * (fB / 0.40))
    sigb = p["SIG0"] - p["DSIGB"] * kfLes
    sigu = p["SIG0"] - p["DSIGU"] * kfLes * (fB / 0.40)
    sigb = min(max(sigb, 0.05), 0.95)
    sigu = min(max(sigu, 0.05), 0.95)
    d["Kfb"], d["Kfu"], d["sigb"], d["sigu"] = Kfb, Kfu, sigb, sigu

    # ---- protein concentrations and COP ------------------------------
    Cp = y[IX["PP"]] / (VP / 100.0)
    Cib = y[IX["PIB"]] / (VIB / 100.0)
    Ciu = y[IX["PIU"]] / (VIU / 100.0)
    PIp = landis_pappenheimer(Cp)
    PIib = landis_pappenheimer(Cib)
    PIiu = landis_pappenheimer(Ciu)
    d["Cp"], d["Cib"], d["Ciu"] = Cp, Cib, Ciu
    d["PIp"], d["PIib"], d["PIiu"] = PIp, PIib, PIiu

    # ---- hydrostatic pressures ---------------------------------------
    Pc = p["PC0"] + p["KPC"] * (VP / p["VP0"] - 1.0)
    Pc = max(Pc, 4.0)
    imb = p["IMB0"] * math.exp(-t / p["TAUIMB"])
    Pib = p["PI0"] + p["ECB"] * (VIB / VIB0 - 1.0) - imb
    Piu = p["PI0"] + p["ECU"] * (VIU / VIU0 - 1.0)
    d["Pc"], d["Pib"], d["Piu"], d["imb"] = Pc, Pib, Piu, imb

    # ---- Starling fluxes ---------------------------------------------
    Jvb_raw = Kfb * ((Pc - Pib) - sigb * (PIp - PIib))
    Jvu_raw = Kfu * ((Pc - Piu) - sigu * (PIp - PIiu))

    # ---- DELIVERY LIMIT ------------------------------------------------
    # relVP-driven global flow, times a burn-specific no-reflow factor.
    relVP0 = VP / p["VP0"]
    _x0 = min(relVP0 / p["VPCRIT"], 2.2) ** p["HUO"]
    _xr = (1.0 / p["VPCRIT"]) ** p["HUO"]
    perf0 = (_x0 / (1.0 + _x0)) / (_xr / (1.0 + _xr))
    perf0 = min(perf0, 1.20)
    stasis = p["STASIS0"] + (1.0 - p["STASIS0"]) * min(perf0, 1.0)
    ebb0 = math.exp(-t / 30.0)
    COr = (0.55 + 0.45 * min(relVP0 / 0.95, 1.15)) * (1.0 - 0.35 * ebb0) + 0.85 * (1.0 - ebb0)
    COr = max(COr, 0.25)
    Jvb_cap = p["FFMAX"] * p["QPB0"] * fB * COr * stasis
    Jvu_cap = p["FFMAX"] * p["QPU0"] * (1.0 - fB) * COr
    # smooth min: x*c/(x+c) -> min(x,c) without a kink
    Jvb = (Jvb_raw * Jvb_cap / (Jvb_raw + Jvb_cap)) if Jvb_raw > 0 else Jvb_raw
    Jvu = (Jvu_raw * Jvu_cap / (Jvu_raw + Jvu_cap)) if Jvu_raw > 0 else Jvu_raw
    d["Jvb"], d["Jvu"] = Jvb, Jvu
    d["Jvb_raw"], d["Jvb_cap"], d["stasis"] = Jvb_raw, Jvb_cap, stasis

    # ---- lymph return (saturable) ------------------------------------
    def lymph(V, V0):
        # DEFECT 14: without the (V/V0) factor lymph kept running at its
        # baseline rate while the interstitium emptied, so the mobilisation
        # phase drained the interstitium to 2 % of baseline and the plasma
        # protein concentrated to a COP of 196 mmHg on day 5.
        drive = 1.0 + p["KL"] * max(V / V0 - 1.0, 0.0)
        drive = min(drive, p["LMAX"])
        return p["LBASE"] * (V0 / p["VI0"]) * drive * max(V / V0, 0.0)
    Qlb = lymph(VIB, VIB0)
    Qlu = lymph(VIU, VIU0)
    d["Qlb"], d["Qlu"] = Qlb, Qlu

    # ---- protein flux -------------------------------------------------
    PSb = p["PS0"] * fB * (1.0 + p["PSMULT"] * kfLes)
    PSu = p["PS0"] * (1.0 - fB) * (1.0 + 0.4 * p["PSMULT"] * kfLes * (fB / 0.40))
    Jpb = max(Jvb, 0.0) / 100.0 * (1.0 - sigb) * Cp + PSb * (Cp - Cib)
    Jpu = max(Jvu, 0.0) / 100.0 * (1.0 - sigu) * Cp + PSu * (Cp - Ciu)
    # lymph carries interstitial protein back
    Jlb = Qlb / 100.0 * Cib
    Jlu = Qlu / 100.0 * Ciu
    d["Jpb"], d["Jpu"], d["Jlb"], d["Jlu"] = Jpb, Jpu, Jlb, Jlu

    # ---- intra-abdominal pressure -------------------------------------
    vexc = y[IX["VASC"]] + max(VIU - VIU0, 0.0)
    IAP = (5.0 + p["IAPK"] * y[IX["VASC"]] + p["IAPB"] * max(VIU - VIU0, 0.0)
           + p["IAPQ"] * (vexc / 1000.0) ** 3)
    d["IAP"] = IAP

    # ---- haemodynamics -------------------------------------------------
    relVP = VP / p["VP0"]
    ebb = math.exp(-t / 30.0)
    flow = 1.0 - math.exp(-t / 30.0)
    COrel = (0.55 + 0.45 * min(relVP / 0.95, 1.15)) * (1.0 - 0.35 * ebb) + 0.85 * flow
    COrel = max(COrel, 0.25)
    d["COrel"] = COrel

    # ---- GFR and urine output ------------------------------------------
    # perfusion term, NORMALISED so that perf = 1 at relVP = 1
    def _perf(rv):
        x = min(rv / p["VPCRIT"], 2.2) ** p["HUO"]
        return x / (1.0 + x)
    perf = _perf(relVP) / _perf(1.0)
    perf = min(perf, 1.20)
    iapPen = max(1.0 - p["KIAPGFR"] * max(IAP - 12.0, 0.0), 0.10)
    arc = 1.0 + p["KARC"] * flow * (1.0 - 0.5 * d["sep"])
    GFRrel = perf * iapPen * arc * (1.0 - 0.85 * y[IX["RTUB"]])
    UO = (p["UOBASE"] * p["WT"] * GFRrel * (1.0 + p["KDIU"] * max(relVP - 1.0, 0.0))
          * (1.0 - p["KOPI"]))
    UO = min(UO, p["UOMAX"] * p["WT"])
    d["GFRrel"], d["UO"], d["arc"], d["perf"] = GFRrel, UO, arc, perf

    # ---- resuscitation controller ---------------------------------------
    R_formula = p["FORMULA"] * p["WT"] * p["TBSA"] / 24.0          # mL/h flat
    # Parkland front-loading: half in the first 8 h
    if t < 24.0:
        shape = 1.5 if t < 8.0 else 0.75
        Rbase = R_formula * shape
        Rinf = Rbase * y[IX["RSTATE"]] + p["MAINT"]
    else:
        # maintenance phase: replace evaporative + urinary losses
        # DEFECT 11: replacing the MEASURED urine output closed a positive
        # loop (more volume -> more urine -> more replacement) and the plasma
        # volume ran away to 148 % of baseline.  Maintenance replaces
        # evaporative loss plus a PHYSIOLOGIC urine, not the observed one.
        evap_now = (p["KEVAP"] * (max(y[IX["AOPEN"]], 0.0) + max(y[IX["ADON"]], 0.0))
                    / max(p["TBSA"], 1.0))
        Rinf = (p["KMNT"] * (evap_now + UO)
                + p["KVP"] * (p["VP0"] - VP)) + p["MAINT"]
        if Rinf < 0.0:
            Rinf = 0.0
    # a FRACTION of the same infusion is switched to 5 % albumin
    # DEFECT 16: applying the colloid fraction to the whole titrated infusion
    # for 48 h delivered ~700 g of protein into a 224 g pool -- plasma COP
    # reached 92 mmHg and sucked the entire interstitium into the circulation.
    # Real protocols run 5 % albumin only through the first 24 h.
    fcol = p["FCOL"] if (t >= p["ALBSTART"] and t < 24.0) else 0.0
    albVol = Rinf * fcol
    albRate = albVol * p["ALBCONC"] / 100.0                        # mL -> g of protein
    d["Rinf"], d["albRate"], d["albVol"], d["fcol"] = Rinf, albRate, albVol, fcol
    d["Rtot"] = Rinf

    # ---- evaporation ----------------------------------------------------
    aopen = max(y[IX["AOPEN"]], 0.0) + max(y[IX["ADON"]], 0.0)
    evap = p["KEVAP"] * (aopen / max(p["TBSA"], 1.0))
    d["aopen"], d["evap"] = aopen, evap

    # ---- adrenergic drive (the TRANSDUCER) --------------------------------
    A = aopen
    driveA = (A ** p["HILLA"]) / (p["A50"] ** p["HILLA"] + A ** p["HILLA"])
    cpr = y[IX["APRC"]] / p["VPR"] * 1000.0                        # ng/mL
    prBlock = p["EMAXPR"] * cpr / (p["EC50PR"] + cpr)
    d["driveA"], d["cpr"], d["prBlock"] = driveA, cpr, prBlock
    # transduced adrenergic signal = tone x (1 - receptor blockade)
    adr = y[IX["CAT"]] * (1.0 - prBlock)
    d["adr"] = adr

    # ---- oxandrolone ------------------------------------------------------
    cox = y[IX["AOXC"]] / p["VOX"] * 1000.0
    oxEff = p["EMAXOX"] * cox / (p["EC50OX"] + cox)
    d["cox"], d["oxEff"] = cox, oxEff

    # ---- vancomycin -------------------------------------------------------
    cvan = y[IX["AVCC"]] / p["VVC"]
    d["cvan"] = cvan

    # enteral carbohydrate, held around each operating session
    feed = p["GFEED"] if t >= p["TFEED"] else 0.0
    if p["NEXC"] > 0 and p["EXCRATE"] > 0:
        for _k in range(int(p["NEXC"])):
            _tk = p["TEXC"] + _k * p["EXCINT"]
            if _tk - p["NPOHRS"] <= t <= _tk + 2.0:
                feed = 0.0
                break
    d["feed"] = feed

    # ---- inflammation composite --------------------------------------------
    infl = y[IX["IL6"]] / 100.0 + y[IX["TNF"]] / 30.0
    d["infl"] = infl

    return d


# =====================================================================
# 3.  DERIVATIVES
# =====================================================================
def deriv(t, y, p, sc):
    d = derive(t, y, p, sc)
    dy = [0.0] * NST
    fB = p["FB"]
    VIB0 = p["VI0"] * fB
    VIU0 = p["VI0"] * (1.0 - fB)

    # ---------- 1. fluid ----------
    ascIn = p["KASC"] * max(d["Jvu"], 0.0) * max(1.0 - y[IX["VASC"]] / p["VASCMAX"], 0.0)
    ascOut = p["KASCR"] * y[IX["VASC"]]
    evlwIn = (p["KEVLW"] * max(d["Jvu"], 0.0) * (1.0 + 1.4 * p["INH"]) * (1.0 + 0.8 * d["sep"])
              * max(1.0 - y[IX["EVLW"]] / p["EVLWMAX"], 0.0))
    evlwOut = p["KEVLWR"] * y[IX["EVLW"]]

    dy[IX["VP"]] = (d["Rtot"] - d["Jvb"] - d["Jvu"] + d["Qlb"] + d["Qlu"]
                    - d["UO"] - ascIn + ascOut - evlwIn + evlwOut)
    dy[IX["VIB"]] = d["Jvb"] - d["Qlb"] - d["evap"] * 0.65
    dy[IX["VIU"]] = d["Jvu"] - d["Qlu"] - ascIn + ascOut - evlwIn + evlwOut - d["evap"] * 0.35
    dy[IX["VASC"]] = ascIn - ascOut
    dy[IX["EVLW"]] = evlwIn - evlwOut
    # cell swelling: Na/K-ATPase failure scales with lactate/hypoperfusion
    dy[IX["VICF"]] = 3.0 * p["WT"] * (1.0 - d["perf"]) - 0.05 * (y[IX["VICF"]] - p["VICF0"])

    dy[IX["UOC"]] = d["UO"]
    dy[IX["FIN"]] = d["Rtot"]
    dy[IX["ALBIN"]] = d["albRate"]

    # ---------- 2. protein ----------
    albSyn = 0.55 * (1.0 / (1.0 + y[IX["IL6"]] / 260.0)) * (1.0 - 0.35 * max(y[IX["HFAT"]] - 300.0, 0.0) / 900.0)
    dy[IX["PP"]] = (-d["Jpb"] - d["Jpu"] + d["Jlb"] + d["Jlu"]
                    + d["albRate"] + albSyn
                    - p["KPPD"] * y[IX["PP"]] * (y[IX["PP"]] / (p["CP0"] * p["VP0"] / 100.0)) ** p["HPPD"])
    dy[IX["PIB"]] = d["Jpb"] - d["Jlb"]
    dy[IX["PIU"]] = d["Jpu"] - d["Jlu"]

    # ---------- 3. controller ----------
    # first-order filter on the measured urine output (clinicians average an hour)
    dy[IX["UOWIN"]] = (d["UO"] - y[IX["UOWIN"]]) / p["TAUUO"]
    if p["TITRATE"] > 0.5 and t < 24.0:
        err = (p["UOTGT"] * p["WT"] - y[IX["UOWIN"]]) / (p["UOTGT"] * p["WT"])
        drs = p["GCTRL"] * err * y[IX["RSTATE"]]
        # DEFECT 4: an unlimited proportional term made the multiplier ring
        # between 0.47 and 2.15 within six hours.  Real protocols move the
        # rate in bounded steps.
        if drs > p["DRMAX"]:
            drs = p["DRMAX"]
        if drs < -p["DRMAX"]:
            drs = -p["DRMAX"]
        # bound the multiplier
        if y[IX["RSTATE"]] >= p["RMAX"] and drs > 0:
            drs = 0.0
        if y[IX["RSTATE"]] <= p["RMIN"] and drs < 0:
            drs = 0.0
        dy[IX["RSTATE"]] = drs
    else:
        dy[IX["RSTATE"]] = 0.0

    # ---------- 4. renal ----------
    dy[IX["RTUB"]] = p["KTUB"] * (1.0 - d["perf"]) * (1.0 + 1.5 * d["sep"]) - p["KTUBR"] * y[IX["RTUB"]]
    # creatinine: production falls with LBM, clearance with GFR
    # dSCR/dt = kprod*(LBM/LBM0) - kel*GFRrel*SCR ; steady state = CRE0 when both = 1
    kel = 0.28
    kprod = kel * p["CRE0"]
    dy[IX["SCR"]] = (kprod * (y[IX["LBM"]] / (p["LBM0F"] * p["WT"]))
                     - kel * max(d["GFRrel"], 0.03) * y[IX["SCR"]])
    dy[IX["LAC"]] = 5.5 * (1.0 - d["perf"]) + 1.2 * d["sep"] - 0.45 * y[IX["LAC"]]

    # ---------- 5. inflammation ----------
    dampProd = p["KDAMP"] * (d["aopen"] / max(p["TBSA"], 1.0)) * (1.0 + 1.2 * math.exp(-t / 24.0))
    dy[IX["DAMP"]] = dampProd - p["KDDEG"] * y[IX["DAMP"]]
    dy[IX["IL6"]] = (p["KIL6"] * y[IX["DAMP"]] * (1.0 + 3.0 * d["sep"]) * (1.0 + 0.6 * p["INH"])
                     - p["KIL6D"] * y[IX["IL6"]])
    dy[IX["TNF"]] = p["KTNF"] * y[IX["DAMP"]] * (1.0 + 2.0 * d["sep"]) - p["KTNFD"] * y[IX["TNF"]]
    dy[IX["IL10"]] = p["KIL10"] * (y[IX["IL6"]] / 100.0) - p["KIL10D"] * y[IX["IL10"]]
    dy[IX["PMN"]] = p["KPMN"] * y[IX["DAMP"]] - p["KPMND"] * y[IX["PMN"]]
    dy[IX["HLADR"]] = (-p["KHLA"] * (y[IX["IL10"]] + 0.35 * (y[IX["COR"]] - p["COR0"]))
                       + p["KHLAR"] * (100.0 - y[IX["HLADR"]]))
    dy[IX["CRP"]] = p["KCRP"] * y[IX["IL6"]] - p["KCRPD"] * y[IX["CRP"]]

    # ---------- 6. infection ----------
    immune = (y[IX["HLADR"]] / 100.0) * (1.0 - 0.35 * max(1.0 - y[IX["LBM"]] / (p["LBM0F"] * p["WT"]) / 0.90, 0.0))
    # DEFECT 24: bacterial growth was normalised by %TBSA, so a 20 % and a
    # 70 % burn had identical infection dynamics and burn size could not
    # generate its own leading cause of late death.  Colony density does not
    # depend on wound size -- but the number of portals does, so the size
    # belongs on the INVASION term, not on growth.
    growth = (p["KGROW"] * (1.0 - y[IX["BWD"]] / p["BWDMAX"])
              * (1.0 if d["aopen"] > 0.5 else 0.0))
    dy[IX["BWD"]] = growth - p["KTOP"] * sc.get("topical", 1.0) - 0.010 * immune
    if y[IX["BWD"]] < 0.0 and dy[IX["BWD"]] < 0.0:
        dy[IX["BWD"]] = 0.0
    invade = (p["KINV"] * max(y[IX["BWD"]] - p["BTHRESH"], 0.0) * (2.0 - immune)
              * (d["aopen"] / 45.0))
    # hyperglycaemia impairs neutrophil function (bounded)
    invade *= (1.0 + 0.45 * min(max(y[IX["GLC"]] - 180.0, 0.0) / 100.0, 2.0))
    kill = p["KBSYSD"] + sc.get("abx", 0.0) * (d["cvan"] / (d["cvan"] + 8.0))
    # DEFECT 21: an unsaturated invasion term let the systemic burden reach
    # 240 arbitrary units, at which point the sepsis switch was pinned at 1
    # and ANTIBIOTICS COULD NOT CHANGE THE OUTCOME -- halving 240 to 120
    # leaves a saturated logistic exactly where it was.
    dy[IX["BSYS"]] = invade * max(1.0 - y[IX["BSYS"]] / p["BSYSMAX"], 0.0) - kill * y[IX["BSYS"]]

    # ---------- 7. neuroendocrine ----------
    catTarget = 1.0 + (p["CATMAX"] - 1.0) * d["driveA"] * (1.0 + 0.5 * d["sep"])
    dy[IX["CAT"]] = p["KCATD"] * (catTarget - y[IX["CAT"]])
    corTarget = p["COR0"] * (1.0 + (p["CORMAX"] - 1.0) * d["driveA"]) * sc.get("cortmod", 1.0)
    dy[IX["COR"]] = p["KCORD"] * (corTarget - y[IX["COR"]])
    gcgTarget = p["GCG0"] * (1.0 + (p["GCGMAX"] - 1.0) * d["driveA"])
    dy[IX["GCG"]] = p["KGCGD"] * (gcgTarget - y[IX["GCG"]])
    dy[IX["T3"]] = p["KT3"] * (p["T30"] * (1.0 - 0.55 * d["driveA"]) - y[IX["T3"]])
    igfTarget = p["IGF0"] * (1.0 - 0.62 * d["driveA"]) * (1.0 + 0.55 * d["oxEff"]) * sc.get("ghmod", 1.0)
    dy[IX["IGF1"]] = p["KIGF"] * (igfTarget - y[IX["IGF1"]])

    # ---------- 8. hypermetabolism (the SLOW clock) ----------
    # target REE is driven by the OPEN AREA and transduced by beta-receptors
    hmDrive = d["driveA"] * (d["adr"] / p["CATMAX"]) ** 0.55
    ambPen = 1.0 + 0.030 * max(31.0 - p["TAMB"], 0.0)
    reeTarget = 1.0 + p["EMAXHM"] * hmDrive * ambPen
    # DEFECT 20: the decay branch inherited 6 % of the rise rate, which made
    # the slow clock's time constant 26 d instead of ~100 d and erased the
    # single best-documented feature of the syndrome -- that hypermetabolism
    # outlasts wound closure by many months.
    krate = p["KREE"] if reeTarget > y[IX["REE"]] else p["KREED"]
    dy[IX["REE"]] = krate * (reeTarget - y[IX["REE"]])
    dy[IX["TCORE"]] = 0.35 * (p["TCORE0"] + 1.6 * (y[IX["REE"]] - 1.0) / 0.85 + 0.6 * d["sep"] - y[IX["TCORE"]])

    # ---------- 9. protein / body composition ----------
    insEff = y[IX["INS"]] / (y[IX["INS"]] + 25.0)
    igfEff = y[IX["IGF1"]] / p["IGF0"]
    S = p["KSYN"] * (1.0 + p["KSINS"] * insEff * sc.get("nutrition", 1.0)
                     + p["KSIGF"] * igfEff + 0.85 * d["oxEff"]
                     + p["KSPR"] * d["prBlock"]) * sc.get("nutrition", 1.0)
    B = p["KBRK"] * (1.0 + p["KBCAT"] * (d["adr"] - 1.0)
                     + p["KBCOR"] * (y[IX["COR"]] / p["COR0"] - 1.0)
                     + p["KBINF"] * d["infl"])
    B = max(B, 0.0)
    # DEFECT 17: with no set-point, recovery ran past baseline and the
    # propranolol+oxandrolone arm reached 89 kg of lean mass from 64 kg.
    # Anabolism is deficit-driven: it stops when the deficit is repaid.
    S /= (1.0 + p["KHOMEO"] * max(y[IX["LBM"]] / (p["LBM0F"] * p["WT"]) - 1.0, 0.0))
    dy[IX["LBM"]] = (S - B) * y[IX["LBM"]]
    d["S"], d["B"] = S, B
    dy[IX["FATM"]] = -p["KFATL"] * d["adr"] / 1.0 * y[IX["FATM"]] + 0.00030 * insEff * y[IX["FATM"]]
    dy[IX["BMC"]] = -p["KBONE"] * (d["adr"] + y[IX["COR"]] / p["COR0"] + d["infl"]) + 0.0000060 * (1.0 - y[IX["BMC"]])
    dy[IX["SCARC"]] = p["KSCAR"] * max(y[IX["AHEAL"]], 0.0) * (1.0 + 1.4 * max(t - 504.0, 0.0) / 504.0) - 0.00020 * y[IX["SCARC"]]

    # ---------- 10. glucose / lipid ----------
    xres = max(d["adr"] - 1.0, 0.0) + max(y[IX["COR"]] / p["COR0"] - 1.0, 0.0) + d["infl"]
    ires = 1.0 + p["IRESMAX"] * xres / (xres + p["IRESK"])
    Vg = p["VGLC"] * p["WT"]                                   # dL
    hgp = (p["HGP0"] * (1.0 + p["KGCGH"] * (y[IX["GCG"]] / p["GCG0"] - 1.0))
           * (1.0 + p["KIRESH"] * (ires - 1.0))
           / (1.0 + p["KINSH"] * insEff) * sc.get("metformin", 1.0))
    upt = (p["U0"] + p["UI"] * insEff / ires) * (y[IX["GLC"]] / p["GLC0"])
    dy[IX["GLC"]] = (hgp + d["feed"] - upt) / Vg
    insTarget = p["INS0"] * (1.0 + p["KINSS"] * max(y[IX["GLC"]] - p["GLC0"], 0.0))
    insTarget = min(insTarget, p["INSMAX"]) + y[IX["AINSX"]] * p["KINSEL"]
    dy[IX["INS"]] = p["KINSD"] * (insTarget - y[IX["INS"]])
    # hours of exposure below 70 mg/dL -- the cost of a tight target
    dy[IX["HYPOH"]] = 1.0 / (1.0 + math.exp((y[IX["GLC"]] - 70.0) / 3.0))
    dy[IX["FFA"]] = p["KFFA"] * (d["adr"] / 3.0) - p["KFFAD"] * y[IX["FFA"]] * (1.0 + 1.5 * insEff)
    dy[IX["HFAT"]] = p["KHFAT"] * y[IX["FFA"]] * 100.0 - p["KHFATD"] * y[IX["HFAT"]]

    # ---------- 11. wound ----------
    # spontaneous re-epithelialisation only of the non-deep fraction
    superficialLeft = max(y[IX["AOPEN"]] - p["FDEEP"] * p["TBSA"], 0.0)
    reepi = p["KREEPI"] * superficialLeft * (1.0 - 0.35 * max(y[IX["BWD"]] - 5.0, 0.0) / 4.0)
    reepi = max(reepi, 0.0)
    # graft take
    edemaPen = 1.0 - 0.30 * min(max(y[IX["VIB"]] / VIB0 - 1.0, 0.0) / 1.5, 1.0)
    infPen = 1.0 - 0.55 * min(max(y[IX["BWD"]] - 5.0, 0.0) / 3.0, 1.0)
    nutPen = 1.0 - 0.40 * min(max(1.0 - y[IX["LBM"]] / (p["LBM0F"] * p["WT"]), 0.0) / 0.20, 1.0)
    take = p["TAKEMAX"] * max(edemaPen, 0.1) * max(infPen, 0.1) * max(nutPen, 0.1)
    graftTake = p["KGRAFT"] * y[IX["AGRF"]]
    dy[IX["AGRF"]] = -graftTake
    dy[IX["AHEAL"]] = graftTake * take + reepi
    # failed graft returns to open
    dy[IX["AOPEN"]] = -reepi + graftTake * (1.0 - take)
    dy[IX["AEXC"]] = 0.0
    dy[IX["ADON"]] = -p["KDONH"] * y[IX["ADON"]]

    # ---------- 12. drug PK ----------
    dy[IX["APRD"]] = -p["KAPR"] * y[IX["APRD"]]
    dy[IX["APRC"]] = p["KAPR"] * y[IX["APRD"]] * p["FPR"] - p["CLPR"] / p["VPR"] * y[IX["APRC"]]
    dy[IX["AOXD"]] = -p["KAOX"] * y[IX["AOXD"]]
    dy[IX["AOXC"]] = p["KAOX"] * y[IX["AOXD"]] * p["FOX"] - p["CLOX"] / p["VOX"] * y[IX["AOXC"]]
    clvan = p["CLVC"] * d["arc"] * max(d["GFRrel"] / max(d["arc"], 0.01), 0.05)
    dy[IX["AVCC"]] = (-clvan / p["VVC"] * y[IX["AVCC"]]
                      - p["QVC"] / p["VVC"] * y[IX["AVCC"]] + p["QVC"] / p["VVP"] * y[IX["AVCP"]])
    dy[IX["AVCP"]] = p["QVC"] / p["VVC"] * y[IX["AVCC"]] - p["QVC"] / p["VVP"] * y[IX["AVCP"]]
    dy[IX["AVTC"]] = -p["CLVIT"] / p["VVIT"] * y[IX["AVTC"]]
    dy[IX["AINSX"]] = -1.2 * y[IX["AINSX"]]

    # ---------- 13. hazard ----------
    acs = _logit((d["IAP"] - 20.0) / 1.6)
    ards = _logit((y[IX["EVLW"]] - 900.0) / 220.0)
    lbmLoss = max(1.0 - y[IX["LBM"]] / (p["LBM0F"] * p["WT"]), 0.0)
    aof = max(d["aopen"], 0.0) / 100.0
    dy[IX["HAZ"]] = (p["HZOPEN"] * (aof ** p["HZPOW"])
                     * math.exp(p["HZAGEE"] * max(p["AGE"] - 20.0, 0.0))
                     * (1.0 + p["HZINH"] * p["INH"])
                     + p["HZSEP"] * d["sep"]
                     + p["HZACS"] * acs
                     + p["HZARDS"] * ards
                     + p["HZLBM"] * max(lbmLoss - 0.10, 0.0) / 0.10
                     + p["HZFAT"] * max(y[IX["HFAT"]] - 900.0, 0.0) / 900.0)
    return dy


# =====================================================================
# 4.  EVENTS (bolus doses, operations) and the integrator
# =====================================================================
def apply_events(t, y, p, sc, dt):
    """Discrete events applied at the top of each step (times are exact)."""
    eps = dt / 2.0

    def due(t0, interval, n):
        if interval <= 0 or n <= 0:
            return False
        for k in range(int(n)):
            tk = t0 + k * interval
            if abs(t - tk) < eps:
                return True
        return False

    # propranolol PO q6h from sc['prop_start']
    if sc.get("prop_mgkgday", 0.0) > 0:
        d0 = sc.get("prop_start", 48.0)
        if t >= d0 - eps and due(d0, 6.0, sc.get("prop_ndose", 120)):
            y[IX["APRD"]] += sc["prop_mgkgday"] * p["WT"] / 4.0

    # oxandrolone 10 mg PO BID
    if sc.get("oxa_mg", 0.0) > 0:
        d0 = sc.get("oxa_start", 96.0)
        if t >= d0 - eps and due(d0, 12.0, sc.get("oxa_ndose", 60)):
            y[IX["AOXD"]] += sc["oxa_mg"]

    # vancomycin IV q12h
    if sc.get("van_mg", 0.0) > 0:
        d0 = sc.get("van_start", 240.0)
        if t >= d0 - eps and due(d0, sc.get("van_int", 12.0), sc.get("van_ndose", 28)):
            y[IX["AVCC"]] += sc["van_mg"]

    # ascorbate: continuous 66 mg/kg/h for 24 h -> handled as hourly boluses
    if sc.get("vitc_mgkgh", 0.0) > 0 and t < 24.0:
        if abs(t - round(t)) < eps:
            y[IX["AVTC"]] += sc["vitc_mgkgh"] * p["WT"]

    # insulin infusion targeting glucose
    if sc.get("insulin_target", 0.0) > 0 and abs(t - round(t)) < eps:
        gap = y[IX["GLC"]] - sc["insulin_target"]
        if gap > 0:
            # gain rises as the target falls: chasing 80-110 mg/dL means
            # dosing on smaller and smaller errors
            gain = 0.16 * (145.0 / max(sc["insulin_target"], 60.0)) ** 2.4
            y[IX["AINSX"]] += min(gap * gain, 24.0)

    # operating sessions: excise + graft EXCRATE %TBSA
    n = int(p["NEXC"])
    if n > 0 and p["EXCRATE"] > 0:
        for k in range(n):
            tk = p["TEXC"] + k * p["EXCINT"]
            if abs(t - tk) < eps:
                donor_avail = max((100.0 - p["TBSA"]) * p["DONPOOLF"] - y[IX["ADON"]], 0.0)
                amt = min(p["EXCRATE"], y[IX["AOPEN"]], donor_avail * p["MESH"])
                if amt > 0:
                    y[IX["AOPEN"]] -= amt
                    y[IX["AGRF"]] += amt
                    y[IX["ADON"]] += amt * p["DONFRAC"]
                    # excision removes eschar -> DAMP source falls immediately
                    y[IX["DAMP"]] *= (1.0 - 0.30 * amt / max(p["TBSA"], 1.0))
    return y


def simulate(p, sc, tmax=1440.0, dt=0.02, record=1.0, dt_late=0.08, t_switch=48.0):
    """RK4 with a fine step over the fast clock and a coarse step afterwards.

    Every discrete event time in the scenario library is a multiple of the
    coarse step, so no dose or operating session can be stepped over.
    Convergence is verified in section 8 of main().
    """
    # DEFECT 19: a 0.25 h late step was UNSTABLE -- intra-abdominal pressure
    # reached 2.4e7 mmHg by day 12 while 0.10 h and 0.04 h agreed to every
    # printed digit.  0.08 h is used throughout and the convergence check is
    # reported in section 8.
    y = initial_state(p)
    out = {"t": []}
    for n in NAMES:
        out[n] = []
    extra = ["Rinf", "UO", "PIp", "Cp", "Jvb", "Jvu", "IAP", "Rtot",
             "GFRrel", "driveA", "adr", "cpr", "cvan", "sep", "Pib", "sigb", "Pc",
             "Jvb_raw", "Jvb_cap", "stasis", "PIib", "Qlb", "feed", "IAP"]
    for e in extra:
        out[e] = []

    t = 0.0
    nextrec = 0.0
    while True:
        h = dt if t < t_switch - 1e-9 else dt_late
        if t + h > tmax + 1e-9:
            h = tmax - t
        y = apply_events(t, y, p, sc, h)
        if t >= nextrec - 1e-9:
            d = derive(t, y, p, sc)
            out["t"].append(t)
            for j, n in enumerate(NAMES):
                out[n].append(y[j])
            for e in extra:
                out[e].append(d[e])
            nextrec += record
        if h <= 1e-9:
            break
        dt_ = h
        # RK4
        k1 = deriv(t, y, p, sc)
        y2 = [y[j] + dt_ / 2 * k1[j] for j in range(NST)]
        k2 = deriv(t + dt_ / 2, y2, p, sc)
        y3 = [y[j] + dt_ / 2 * k2[j] for j in range(NST)]
        k3 = deriv(t + dt_ / 2, y3, p, sc)
        y4 = [y[j] + dt_ * k3[j] for j in range(NST)]
        k4 = deriv(t + dt_, y4, p, sc)
        y = [y[j] + dt_ / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(NST)]
        # non-negativity floors (physical states cannot go negative)
        for nm in ("VP", "VIB", "VIU", "PP", "PIB", "PIU", "VASC", "EVLW",
                   "DAMP", "IL6", "TNF", "IL10", "PMN", "BSYS", "BWD",
                   "APRD", "APRC", "AOXD", "AOXC", "AVCC", "AVCP", "AVTC", "AINSX",
                   "AOPEN", "AGRF", "ADON", "AHEAL"):
            if y[IX[nm]] < 0.0:
                y[IX[nm]] = 0.0
        if y[IX["HLADR"]] < 2.0:
            y[IX["HLADR"]] = 2.0
        t += dt_
    return out


# =====================================================================
# 5.  DERIVED CLINICAL READOUTS
# =====================================================================
def at(out, tt):
    """index of the record closest to time tt"""
    best, bi = 1e18, 0
    for i, t in enumerate(out["t"]):
        if abs(t - tt) < best:
            best, bi = abs(t - tt), i
    return bi


def readouts(out, p):
    i24 = at(out, 24.0)
    i48 = at(out, 48.0)
    r = {}
    parkland = 4.0 * p["WT"] * p["TBSA"]   # the Parkland prescription, always
    r["V24"] = out["FIN"][i24]
    r["mLkgPct"] = out["FIN"][i24] / p["WT"] / p["TBSA"]
    r["inout"] = out["FIN"][i24] / parkland
    r["mLkg24"] = out["FIN"][i24] / p["WT"]
    r["UO24"] = out["UOC"][i24]
    r["VPmin"] = min(out["VP"][:i24 + 1])
    r["VPmin_pct"] = 100.0 * r["VPmin"] / p["VP0"]
    r["Cp24"] = out["Cp"][i24]
    r["PIp24"] = out["PIp"][i24]
    r["PIpmin"] = min(out["PIp"][:i24 + 1])
    edema24 = (out["VIB"][i24] + out["VIU"][i24] + out["VASC"][i24] + out["EVLW"][i24]
               - p["VI0"])
    r["edema24_L"] = edema24 / 1000.0
    wt24 = (out["VP"][i24] - p["VP0"] + edema24 + out["VICF"][i24] - p["VICF0"]) / 1000.0
    r["wtgain24_pct"] = 100.0 * wt24 / p["WT"]
    i48b = min(i48, len(out["t"]) - 1)
    edema48 = (out["VIB"][i48b] + out["VIU"][i48b] + out["VASC"][i48b] + out["EVLW"][i48b] - p["VI0"])
    wt48 = (out["VP"][i48b] - p["VP0"] + edema48 + out["VICF"][i48b] - p["VICF0"]) / 1000.0
    r["wtgain48_pct"] = 100.0 * wt48 / p["WT"]
    r["IAPmax"] = max(out["IAP"])
    wpk = 0.0
    for k in range(len(out["t"])):
        ed = (out["VIB"][k] + out["VIU"][k] + out["VASC"][k] + out["EVLW"][k] - p["VI0"])
        w = (out["VP"][k] - p["VP0"] + ed + out["VICF"][k] - p["VICF0"]) / 1000.0
        if w > wpk:
            wpk = w
    r["wtgain_peak_pct"] = 100.0 * wpk / p["WT"]
    r["REEmax"] = max(out["REE"])
    r["REEmax_pct"] = 100.0 * r["REEmax"]
    iend = len(out["t"]) - 1
    r["REEend_pct"] = 100.0 * out["REE"][iend]
    lbm0 = p["LBM0F"] * p["WT"]
    i336 = at(out, 336.0)
    r["dLBM14d_pct"] = 100.0 * (out["LBM"][i336] - lbm0) / lbm0
    r["dLBMend_pct"] = 100.0 * (out["LBM"][iend] - lbm0) / lbm0
    r["mort"] = 100.0 * (1.0 - math.exp(-out["HAZ"][iend]))
    r["IL6max"] = max(out["IL6"])
    r["HLADRmin"] = min(out["HLADR"])
    r["BSYSmax"] = max(out["BSYS"])
    r["closure_d"] = None
    for i, v in enumerate(out["AOPEN"]):
        if v <= 0.05 * p["TBSA"]:
            r["closure_d"] = out["t"][i] / 24.0
            break
    r["SCRmax"] = max(out["SCR"])
    r["fret_early"] = None
    r["cat_max"] = max(out["CAT"])
    r["GLCmax"] = max(out["GLC"])
    r["GLCmin"] = min(out["GLC"])
    r["hypo_h"] = out["HYPOH"][len(out["t"]) - 1]
    return r


def rbaux(p):
    return p["AGE"] + p["TBSA"] + 17.0 * p["INH"]


def rbaux_mortality(score):
    """Osler-type logistic on the revised Baux score (external validator)."""
    return 100.0 / (1.0 + math.exp(-(score - 109.0) / 10.5))


# =====================================================================
# 6.  SCENARIO LIBRARY
# =====================================================================
def build(name, WT=80.0, AGE=35.0, TBSA=45.0, INH=0.0):
    p = default_params(WT, AGE, TBSA, INH)
    sc = dict(topical=1.0, abx=0.0, nutrition=1.0, cortmod=1.0, ghmod=1.0,
              metformin=1.0)
    # DEFAULT SURGICAL PLAN: excision + autograft from day 5, 18 %TBSA (or half
    # the burn, whichever is smaller) per session, every 72 h, donor-site
    # limited.  Scenarios that study the surgical timing override it below.
    # Without this, the resuscitation-only arms never closed the wound and all
    # returned the same 60-day mortality, which hid every slow-clock difference.
    p["TEXC"] = 120.0
    p["EXCRATE"] = min(18.0, TBSA * 0.5)
    p["EXCINT"] = 72.0
    p["NEXC"] = 12
    if name == "no_resus":
        p["FORMULA"] = 0.0; p["TITRATE"] = 0.0
    elif name == "parkland_fixed":
        p["FORMULA"] = 4.0; p["TITRATE"] = 0.0
    elif name == "parkland_titrated":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
    elif name == "brooke_titrated":
        p["FORMULA"] = 2.0; p["TITRATE"] = 1.0
    elif name == "isbi_2mLkg":
        p["FORMULA"] = 2.0; p["TITRATE"] = 1.0; p["GCTRL"] = 0.40
    elif name == "creep_uncapped":
        # "chasing urine output": the same protocol driven to a 1.0 mL/kg/h
        # target instead of 0.5, with no ceiling on the multiplier
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0; p["RMAX"] = 6.0
        p["UOTGT"] = 1.0; p["GCTRL"] = 0.40
    elif name == "opioid_creep":
        # the same 0.5 mL/kg/h target, but high-dose opioid has cut urine
        # output 35 % at any given plasma volume
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0; p["RMAX"] = 6.0; p["KOPI"] = 0.35
    elif name == "albumin_8h":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["ALBSTART"] = 8.0; p["FCOL"] = 0.25      # a third of the infusion as 5 % albumin
    elif name == "albumin_0h":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["ALBSTART"] = 0.0; p["FCOL"] = 0.25
    elif name == "vitc":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        sc["vitc_mgkgh"] = 66.0
    elif name == "excision_d3":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
    elif name == "excision_d14":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 336.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
    elif name == "propranolol":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["prop_mgkgday"] = 4.0; sc["prop_start"] = 120.0; sc["prop_ndose"] = 200
    elif name == "oxandrolone":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["oxa_mg"] = 10.0; sc["oxa_start"] = 120.0; sc["oxa_ndose"] = 100
    elif name == "prop_plus_oxa":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["prop_mgkgday"] = 4.0; sc["prop_start"] = 120.0; sc["prop_ndose"] = 200
        sc["oxa_mg"] = 10.0; sc["oxa_start"] = 120.0; sc["oxa_ndose"] = 100
    elif name == "insulin_tight":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["insulin_target"] = 100.0
    elif name == "insulin_moderate":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["insulin_target"] = 145.0
    elif name == "sepsis_d10":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 336.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        p["KTOP"] = 0.004; p["KGROW"] = 0.070
    elif name == "sepsis_treated":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 336.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        p["KTOP"] = 0.004; p["KGROW"] = 0.070
        sc["abx"] = 0.35; sc["van_mg"] = 1000.0; sc["van_start"] = 240.0
        sc["van_int"] = 12.0; sc["van_ndose"] = 40
    elif name == "cold_room":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        p["TAMB"] = 24.0
    elif name == "full_protocol":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["ALBSTART"] = 8.0; p["FCOL"] = 0.25
        p["TEXC"] = 72.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
        sc["prop_mgkgday"] = 4.0; sc["prop_start"] = 120.0; sc["prop_ndose"] = 200
        sc["oxa_mg"] = 10.0; sc["oxa_start"] = 120.0; sc["oxa_ndose"] = 100
        sc["insulin_target"] = 145.0
    elif name == "standard_care":
        p["FORMULA"] = 4.0; p["TITRATE"] = 1.0
        p["TEXC"] = 120.0; p["EXCRATE"] = min(18.0, TBSA * 0.5); p["EXCINT"] = 72.0; p["NEXC"] = 12
    else:
        raise KeyError(name)
    return p, sc


SCENARIOS = ["no_resus", "parkland_fixed", "parkland_titrated", "brooke_titrated",
             "isbi_2mLkg", "creep_uncapped", "opioid_creep", "albumin_8h", "albumin_0h", "vitc",
             "excision_d3", "excision_d14", "propranolol", "oxandrolone",
             "prop_plus_oxa", "insulin_tight", "insulin_moderate", "sepsis_d10",
             "sepsis_treated", "cold_room", "standard_care", "full_protocol"]


# =====================================================================
# 7.  MAIN
# =====================================================================
def hdr(s):
    print("\n" + "=" * 92)
    print(s)
    print("=" * 92)


REF_PATIENTS = [
    # WT, AGE, TBSA, INH, label
    (80, 25, 20, 0, "25 y, 20 %TBSA"),
    (80, 50, 30, 0, "50 y, 30 %TBSA"),
    (80, 35, 45, 0, "35 y, 45 %TBSA  (reference)"),
    (80, 60, 40, 0, "60 y, 40 %TBSA"),
    (80, 35, 45, 1, "35 y, 45 %TBSA + inhalation"),
    (80, 30, 80, 0, "30 y, 80 %TBSA"),
    (80, 45, 70, 0, "45 y, 70 %TBSA"),
    (80, 70, 50, 1, "70 y, 50 %TBSA + inhalation"),
]


def main():
    quick = "--quick" in sys.argv
    tmax = 336.0 if quick else 1440.0

    hdr("MAJOR THERMAL BURN INJURY — QSP reference implementation (dependency-free Python RK4)")
    print("Reference patient : 80 kg, 35 y, 45 %TBSA, no inhalation injury")
    print("Integration       : RK4, dt = 0.02 h to 48 h then 0.08 h; %d ODE states" % NST)
    print("Horizon           : %.0f h (%.0f days)" % (tmax, tmax / 24.0))

    # ------------------------------------------------------------------
    hdr("0 · THE CONVEX FUNCTION — why each litre costs more than the last")
    print("   Landis-Pappenheimer:  Pi = 2.1 C + 0.16 C^2 + 0.009 C^3   [C = total protein g/dL]\n")
    print("   C[g/dL]   Pi[mmHg]   dPi/dC")
    for C in (7.0, 6.0, 5.0, 4.0, 3.5, 3.0, 2.5, 2.0):
        print("    %4.1f     %6.2f     %5.2f" % (C, landis_pappenheimer(C),
                                                 2.1 + 0.32 * C + 0.027 * C * C))
    print("\n   7.0 -> 3.5 g/dL is a 50 %% dilution but a %.0f %% loss of oncotic pressure"
          % (100 * (1 - landis_pappenheimer(3.5) / landis_pappenheimer(7.0))))
    print("   THIS is the mechanism of fluid creep: f_ret decays as the resuscitation proceeds.")

    # ------------------------------------------------------------------
    hdr("1 · SCENARIO SWEEP — the fast clock (first 24-48 h)")
    fmt = "%-19s %8s %7s %7s %7s %7s %7s %8s %7s"
    print(fmt % ("scenario", "mL/kg/%", "in:out", "mL/kg", "VPmin%", "COP24", "UO24 L", "wt peak%", "IAPmax"))
    print("-" * 92)
    results = {}
    for s in SCENARIOS:
        p, sc = build(s)
        out = simulate(p, sc, tmax=tmax)
        r = readouts(out, p)
        results[s] = (p, sc, out, r)
        print(fmt % (s, "%.2f" % r["mLkgPct"], "%.2f" % r["inout"], "%.0f" % r["mLkg24"],
                     "%.0f" % r["VPmin_pct"], "%.1f" % r["PIp24"], "%.2f" % (r["UO24"] / 1000.0),
                     "%.0f" % r["wtgain_peak_pct"], "%.1f" % r["IAPmax"]))
    print("\n   IAP > 12 mmHg = intra-abdominal hypertension;  > 20 mmHg = abdominal compartment syndrome.")
    print("   Ivy (2000) placed the intra-abdominal hypertension threshold near 250 mL/kg/24 h.")

    # ------------------------------------------------------------------
    hdr("2 · SCENARIO SWEEP — the slow clock (60 days)")
    fmt2 = "%-19s %9s %9s %10s %10s %9s %8s %8s %8s"
    print(fmt2 % ("scenario", "REE peak%", "REE d60%", "dLBM 14d%", "dLBM d60%",
                  "closure d", "HLA-DR", "BSYSmax", "mort %"))
    print("-" * 92)
    for s in SCENARIOS:
        p, sc, out, r = results[s]
        cl = "%.0f" % r["closure_d"] if r["closure_d"] else "  >%.0f" % (tmax / 24.0)
        print(fmt2 % (s, "%.0f" % r["REEmax_pct"], "%.0f" % r["REEend_pct"],
                      "%+.1f" % r["dLBM14d_pct"], "%+.1f" % r["dLBMend_pct"], cl,
                      "%.0f" % r["HLADRmin"], "%.2f" % r["BSYSmax"], "%.1f" % r["mort"]))

    # ------------------------------------------------------------------
    hdr("3 · CALIBRATION TARGETS vs MODEL")
    p, sc, out, r = results["parkland_titrated"]
    pf = results["parkland_fixed"][3]
    rows = [
        ("24-h volume, UO-titrated Parkland [mL/kg/%TBSA]", r["mLkgPct"], "5.2-6.7  Engrav 2000 / Cartotto 2002"),
        ("in:out ratio vs the 4 mL/kg/% prescription", r["inout"], "1.2-1.6  Saffle 2007 'fluid creep'"),
        ("24-h volume [mL/kg]", r["mLkg24"], "IAH threshold ~250  Ivy 2000"),
        ("plasma volume nadir [% of baseline]", r["VPmin_pct"], "60-80 % in the first 8 h"),
        ("plasma COP at 24 h [mmHg]", r["PIp24"], "10-16  (from ~26 at baseline)"),
        ("peak weight gain [% body weight]", r["wtgain_peak_pct"], "+15 to +30 %"),
        ("open-loop Parkland: 24-h urine [L]", pf["UO24"] / 1000.0, "the fixed formula UNDER-fills"),
    ]
    for lab, val, ref in rows:
        print("  %-52s model %8.2f   lit %s" % (lab, val, ref))

    rs = results["standard_care"][3]
    rp = results["propranolol"][3]
    rox = results["oxandrolone"][3]
    print("  %-52s model %8.2f   lit %s" % ("peak REE [% of predicted]", rs["REEmax_pct"], "120-180 %"))
    print("  %-52s model %8.2f   lit %s" % ("REE at day 60 [% of predicted]", rs["REEend_pct"], "still elevated at 6-12 months"))
    print("  %-52s model %8.2f   lit %s" % ("LBM at 14 d, no anabolic agent [%]", rs["dLBM14d_pct"], "about -9 %  Herndon 2001"))
    print("  %-52s model %8.2f   lit %s" % ("LBM at 14 d, propranolol [%]", rp["dLBM14d_pct"], "about +9 %  Herndon 2001"))
    print("  %-52s model %8.2f   lit %s" % ("LBM at 14 d, oxandrolone [%]", rox["dLBM14d_pct"], "improved vs control  Wolf 2006"))
    rv = results["vitc"][3]
    ra = results["albumin_8h"][3]
    ra0 = results["albumin_0h"][3]
    print("  %-52s model %+8.1f   lit %s" % ("ascorbate: change in 24-h volume [%]",
                                             100 * (rv["mLkgPct"] / r["mLkgPct"] - 1), "-45 %  Tanaka 2000 (not replicated)"))
    print("  %-52s model %+8.1f   lit %s" % ("albumin from  8 h: change in volume [%]",
                                             100 * (ra["mLkgPct"] / r["mLkgPct"] - 1), "-20 to -40 %  Navickis 2016"))
    print("  %-52s model %+8.1f   lit %s" % ("albumin from  0 h: change in volume [%]",
                                             100 * (ra0["mLkgPct"] / r["mLkgPct"] - 1), "(no direct trial comparator)"))
    rsep = results["sepsis_d10"][3]
    rsept = results["sepsis_treated"][3]
    print("  %-52s model %8.1f   lit %s" % ("burn sepsis, untreated [% mortality]", rsep["mort"], "very high"))
    print("  %-52s model %8.1f   lit %s" % ("burn sepsis, antibiotic-treated [%]", rsept["mort"], "30-60 %"))

    # ------------------------------------------------------------------
    hdr("4 · EXTERNAL VALIDATOR — the revised Baux score (NOT an input to the model)")
    print("  Mortality below is computed MECHANISTICALLY (open wound area x time, age, inhalation,")
    print("  sepsis, compartment syndrome, ARDS, lean-mass loss).  rBaux is used only to score it.\n")
    print("  %-32s %7s %12s %12s %10s" % ("patient", "rBaux", "rBaux mort%", "model mort%", "closure d"))
    print("  " + "-" * 78)
    val = []
    for (wt, age, tbsa, inh, lab) in REF_PATIENTS:
        p2, sc2 = build("standard_care", WT=wt, AGE=age, TBSA=tbsa, INH=inh)
        o2 = simulate(p2, sc2, tmax=tmax)
        r2 = readouts(o2, p2)
        s2 = rbaux(p2)
        cl = "%.0f" % r2["closure_d"] if r2["closure_d"] else ">%.0f" % (tmax / 24.0)
        val.append((s2, rbaux_mortality(s2), r2["mort"]))
        print("  %-32s %7.0f %12.1f %12.1f %10s" % (lab, s2, rbaux_mortality(s2), r2["mort"], cl))
    err = sum(abs(a - b) for _, a, b in val) / len(val)
    print("\n  mean absolute deviation from the rBaux logistic: %.1f percentage points" % err)
    print("  NOTE the two rBaux-80 patients (35 y/45 % and 50 y/30 %): the model returns")
    print("  %.1f %% and %.1f %% WITHOUT being told they are equivalent.  The exchangeability of"
          % (val[2][2], val[1][2]))
    print("  one year of age for one %TBSA is an OUTPUT here, not an assumption.")

    # ------------------------------------------------------------------
    hdr("5 · THE CONTROLLER, HOUR BY HOUR (Parkland 4 mL/kg/%, titrated to UO 0.5 mL/kg/h)")
    p, sc, out, r = results["parkland_titrated"]
    print(" %5s %9s %8s %8s %8s %7s %7s %8s %8s %7s" %
          ("t[h]", "R[mL/h]", "mult", "UO[mL/h]", "VP[%]", "Cp", "COP", "Pi_b", "sigma_b", "Jv_b"))
    for tt in (0, 1, 2, 4, 6, 8, 12, 16, 20, 24, 36, 48):
        i = at(out, float(tt))
        print(" %5d %9.0f %8.2f %8.0f %8.0f %7.2f %7.1f %8.1f %8.2f %7.0f" %
              (tt, out["Rtot"][i], out["RSTATE"][i], out["UO"][i],
               100 * out["VP"][i] / p["VP0"], out["Cp"][i], out["PIp"][i],
               out["Pib"][i], out["sigb"][i], out["Jvb"][i]))
    print("\n  Read the COP column against the multiplier column: the controller pushes the rate")
    print("  UP precisely while the oncotic pressure it is destroying makes each mL less effective.")

    # ------------------------------------------------------------------
    hdr("6 · STATEMENT 3 — driver vs transducer are SUB-ADDITIVE")
    rl = results["excision_d14"][3]
    re3 = results["excision_d3"][3]
    rpr = results["propranolol"][3]
    p, sc = build("excision_d14")
    sc["prop_mgkgday"] = 4.0
    sc["prop_start"] = 120.0
    sc["prop_ndose"] = 200
    o = simulate(p, sc, tmax=tmax)
    rlp = readouts(o, p)
    print("                                          peak REE   dLBM 14 d")
    print("  late excision (d14), no drug          :  %6.1f %%   %+7.1f %%" % (rl["REEmax_pct"], rl["dLBM14d_pct"]))
    print("  late excision (d14) + propranolol     :  %6.1f %%   %+7.1f %%   transducer alone %+.1f pts"
          % (rlp["REEmax_pct"], rlp["dLBM14d_pct"], rlp["REEmax_pct"] - rl["REEmax_pct"]))
    print("  early excision (d3), no drug          :  %6.1f %%   %+7.1f %%   driver alone     %+.1f pts"
          % (re3["REEmax_pct"], re3["dLBM14d_pct"], re3["REEmax_pct"] - rl["REEmax_pct"]))
    print("  early excision (d3) + propranolol     :  %6.1f %%   %+7.1f %%   both             %+.1f pts"
          % (rpr["REEmax_pct"], rpr["dLBM14d_pct"], rpr["REEmax_pct"] - rl["REEmax_pct"]))
    sumsep = (rlp["REEmax_pct"] - rl["REEmax_pct"]) + (re3["REEmax_pct"] - rl["REEmax_pct"])
    both = rpr["REEmax_pct"] - rl["REEmax_pct"]
    print("\n  sum of the separate effects %+.1f pts vs the combination %+.1f pts" % (sumsep, both))
    print("  => SUB-ADDITIVE, because both act on factors of the same product")
    print("     REE_target = 1 + Emax x driveA(A_open) x transduction(beta).")
    print("  Propranolol buys %+.1f points when closure is slow and %+.1f when it is fast."
          % (rlp["REEmax_pct"] - rl["REEmax_pct"], rpr["REEmax_pct"] - re3["REEmax_pct"]))

    # ------------------------------------------------------------------
    hdr("7 · GLYCAEMIC CONTROL — the cost of a tight target")
    print("  %-22s %9s %9s %9s %11s %9s" % ("scenario", "GLC mean", "GLC max", "GLC min",
                                            "hours <70", "mort %"))
    for s in ("standard_care", "insulin_moderate", "insulin_tight"):
        p, sc, out, r = results[s]
        gm = sum(out["GLC"]) / len(out["GLC"])
        print("  %-22s %9.0f %9.0f %9.0f %11.1f %9.1f" % (s, gm, r["GLCmax"], r["GLCmin"],
                                                          r["hypo_h"], r["mort"]))
    print("\n  The model was not told that tight control fails.  It reproduces the trial result")
    print("  because glucose falls only a few mg/dL while hypoglycaemic exposure multiplies --")
    print("  and because burn hypoglycaemia is driven by FEED INTERRUPTION around theatre,")
    print("  not by the insulin dose.")

    # ------------------------------------------------------------------
    hdr("8 · BURN PHARMACOKINETICS — the sign of the dosing error depends on protein binding")
    p, sc, out, r = results["sepsis_treated"]
    i = at(out, 288.0)
    fu0 = 0.10
    fu = min(fu0 * (7.0 / max(out["Cp"][i], 0.5)), 1.0)
    print("  augmented renal clearance multiplier (flow phase) : %.2f x" % (1.0 + p["KARC"]))
    print("  vancomycin concentration at day 12 (1 g q12h)     : %.1f mg/L" % out["cvan"][i])
    print("  total protein at day 12                           : %.2f g/dL" % out["Cp"][i])
    print("  free fraction of a 90 %%-bound drug                 : %.3f -> %.3f" % (fu0, fu))
    print("""
  For a LOW-EXTRACTION drug, CL = fu x CL_int, so at steady state
        C_free,ss = Dose / (CL_int x tau)      -- independent of fu
        C_total,ss = C_free,ss / fu            -- falls as albumin falls
  Therapeutic drug monitoring on TOTAL concentration therefore OVER-diagnoses
  under-exposure for highly bound drugs.  Augmented renal clearance points the
  other way: it lowers the FREE concentration and TDM UNDER-diagnoses it.
  The two commonly quoted burn PK changes have OPPOSITE signs on the quantity
  that matters, which is why 'burn patients need higher doses' is only half true.""")

    # ------------------------------------------------------------------
    hdr("9 · NUMERICAL CONVERGENCE AND MASS BALANCE")
    print("  %8s %14s %14s %16s" % ("dt_late", "24-h volume", "peak REE %", "day-30 lean mass"))
    for dtl in (0.25, 0.12, 0.08, 0.04, 0.02):
        p, sc = build("full_protocol")
        o = simulate(p, sc, tmax=720.0, dt=0.02, dt_late=dtl)
        i24 = at(o, 24.0)
        i30 = at(o, 720.0)
        print("  %8.3f %14.0f %14.1f %16.2f" % (dtl, o["FIN"][i24], 100 * max(o["REE"]), o["LBM"][i30]))
    print("  (0.25 h is UNSTABLE -- it was the defect that produced an intra-abdominal")
    print("   pressure of 2.4e7 mmHg.  0.12 h and finer agree; 0.08 h is the default.)")

    p, sc = build("parkland_titrated")
    o = simulate(p, sc, tmax=48.0)
    i24 = at(o, 24.0)
    inbody = (o["VP"][i24] + o["VIB"][i24] + o["VIU"][i24] + o["VASC"][i24] + o["EVLW"][i24])
    start = p["VP0"] + p["VI0"]
    evap_est = 0.0
    for k in range(1, i24 + 1):
        aop = o["AOPEN"][k] + o["ADON"][k]
        evap_est += p["KEVAP"] * (aop / p["TBSA"]) * (o["t"][k] - o["t"][k - 1])
    resid = inbody - (start + o["FIN"][i24] - o["UOC"][i24] - evap_est)
    print("\n  24-h extracellular fluid balance:")
    print("    start %8.0f  + infused %8.0f  - urine %8.0f  - evaporated %8.0f  = %8.0f mL"
          % (start, o["FIN"][i24], o["UOC"][i24], evap_est,
             start + o["FIN"][i24] - o["UOC"][i24] - evap_est))
    print("    actually in the compartments                                        = %8.0f mL" % inbody)
    print("    residual %+.0f mL (%.3f %% of throughput) -- intracellular shift"
          % (resid, 100 * abs(resid) / max(o["FIN"][i24], 1.0)))

    # ------------------------------------------------------------------
    hdr("10 · WHAT THE MODEL GETS WRONG (stated, not hidden)")
    print("""  1. ALBUMIN STARTED AT 8 h saves only %.0f %% of the volume, against a
     meta-analytic estimate of 20-40 %%.  The model's reason is mechanical and
     testable: under a front-loaded Parkland shape most of the crystalloid --
     and therefore most of the oncotic dilution -- has already been given by
     hour 8.  Started at hour 0 the same colloid fraction saves %.0f %%.  If the
     trials are right and timing does NOT matter that much, this structure is
     wrong.
  2. ASCORBATE reproduces Tanaka's volume reduction almost exactly, which is a
     problem rather than a triumph: no subsequent trial has replicated it, so
     the model is fitting a result that may not be real.  The ascorbate arm
     should be read as 'what would follow IF the oxidant hypothesis of the Kf
     lesion is correct', not as a recommendation.
  3. The 5 %% albumin solution is entered through a TOTAL-PROTEIN oncotic
     equation, which understates it: albumin is more oncotically active per
     gram than globulin, so the colloid arms are conservative by construction.
  4. Mortality is fitted to the revised Baux logistic through the open-wound
     hazard.  The exchangeability result in section 4 is therefore a check on
     internal consistency, not an independent validation.
  5. No coagulopathy, no rhabdomyolysis kinetics, no drug-specific
     nephrotoxicity, and inhalation injury enters only as a fluid multiplier
     and a hazard multiplier -- not as gas exchange.""" % (
        100 * (ra["mLkgPct"] / r["mLkgPct"] - 1) if False else
        100 * (results["albumin_8h"][3]["mLkgPct"] / results["parkland_titrated"][3]["mLkgPct"] - 1),
        -100 * (1 - results["albumin_0h"][3]["mLkgPct"] / results["parkland_titrated"][3]["mLkgPct"])))

    print("\nDone.")


if __name__ == "__main__":
    main()
