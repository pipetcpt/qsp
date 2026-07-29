#!/usr/bin/env python3
# =============================================================================
#  arvc_reference_impl.py
#  Arrhythmogenic Right Ventricular Cardiomyopathy (ARVC / arrhythmogenic
#  cardiomyopathy) — dependency-free reference implementation of the QSP model
#  specified in arvc_qsp_model.dot and arvc_mrgsolve_model.R
# =============================================================================
#
#  WHY THIS FILE EXISTS
#  --------------------
#  The R/mrgsolve file is the deliverable model, but it needs R + mrgsolve + a
#  C++ toolchain to run.  This file is a numerical twin of the SAME equations
#  that needs nothing but python3, so every number quoted in the README, in the
#  .dot header and in the .R header can be regenerated and checked by anyone.
#
#      python3 arvc_reference_impl.py            # run all scenarios, print tables
#      python3 arvc_reference_impl.py --check    # PASS/FAIL assertions
#      python3 arvc_reference_impl.py --falsify  # run the load-independent twin
#
#  ---------------------------------------------------------------------------
#  THE FIVE OBSERVATIONS THE MODEL EXISTS TO EXPLAIN
#  ---------------------------------------------------------------------------
#  (A) PKP2 truncating variants are far commoner in the population than ARVC
#      is.  Penetrance in genotype-positive relatives is roughly 35-50% by
#      mid-life, and a large biobank screen found most carriers of ARVC-gene
#      loss-of-function variants had no ARVC phenotype at all.
#      [Carruth 2019 PMID 31638835; Groeneweg 2015 PMID 25820315]
#      => the genotype cannot be the clock.
#
#  (B) Exercise dose is graded and large.  Desmosomal-variant carriers in the
#      highest exercise tertile reached the same phenotype far younger and had
#      about a 3-fold hazard of VT/death versus the lowest.
#      [James 2013 PMID 23871885; Ruwald 2015 PMID 25896080]
#
#  (C) The same phenotype occurs in athletes with NO desmosomal variant
#      ("gene-elusive" ARVC), and those patients had done MORE exercise than
#      the gene-positive ones.  [Sawant 2014 PMID 25516436]
#      => load alone can reach the phenotype; genotype only changes the price
#      per unit load.
#
#  (D) The right ventricle goes first, although PKP2 is expressed equally in
#      both ventricles.  Left-dominant arrhythmogenic cardiomyopathy is
#      overwhelmingly a DSP / FLNC phenomenon, not a PKP2 one.
#      [Corrado 2024 PMID 37844667; Brandao 2023 PMID 37048743]
#
#  (E) The therapy hierarchy is strange if you think of ARVC as one disease
#      with one arrhythmia.  Amiodarone clearly reduces arrhythmia; sotalol,
#      once first-line, came out no better than nothing in the North American
#      registry; flecainide -- a class IC drug, i.e. the CAST drug, in a
#      structural cardiomyopathy -- reduces VT as an add-on; endocardial-only
#      ablation recurs in most patients while combined endo-epicardial ablation
#      does not; and the ICD, which does nothing whatever to the myocardium,
#      is the only therapy with an unambiguous mortality benefit.
#      [Marcus 2009 PMID 19660690; Ermakov 2017 PMID 27939893;
#       Santangeli 2015 PMID 26546346; Berruezo 2012 PMID 22205683]
#
#  ---------------------------------------------------------------------------
#  THREE STRUCTURAL COMMITMENTS
#  ---------------------------------------------------------------------------
#  1. THE CLOCK IS CUMULATIVE MECHANICAL WORK, NOT TIME.
#     A desmosomal variant does not destroy myocardium; it lowers the load at
#     which the intercalated disc fails.  Myocyte loss is therefore a
#     FATIGUE-FAILURE process: rate = beats x (wall stress / reserve)^NMECH,
#     with NMECH fixed a priori at 4 (the Basquin exponent range for
#     load-bearing biological tissue), NOT fitted.  Genotype enters only
#     through the reserve term.  Because wall stress depends on wall thickness,
#     and wall thickness depends on how much myocardium is left, the process is
#     autocatalytic -- which is what produces a decades-long concealed phase
#     followed by rapid overt progression, with no phase-transition parameter.
#
#     The falsifier is one line: PHI_EX = 0 makes damage load-INDEPENDENT and
#     turns this file into a time-based degeneration model.  Then (A), (B), (C)
#     and (D) all invert simultaneously.  Run --falsify.
#
#     This commitment has a direct experimental test in the literature that the
#     model is NOT fitted to: in plakoglobin-deficient mice, LOAD-REDUCING
#     THERAPY (furosemide + nitrate) prevented RV enlargement, prevented the
#     arrhythmic phenotype and preserved conduction -- a drug with no
#     desmosomal, ion-channel or anti-fibrotic action at all.
#     [Fabritz 2011 PMID 21292134; Kirchhof 2006 PMID 17030684]
#
#  2. RV SELECTIVITY IS LAPLACE, NOT BIOLOGY.
#     No chamber-specific gene expression is used for PKP2 disease.  The RV
#     free wall is thin (~4 mm vs ~9 mm), so baseline wall stress per unit
#     pressure is higher; and during exercise RV end-systolic wall stress rises
#     about 125% while LV rises about 14%, because pulmonary vascular
#     resistance falls much less than systemic.  [La Gerche 2011 PMID 21085033;
#     La Gerche 2012 PMID 22160404; La Gerche 2017 PMID 28957535]
#     Raised to the power NMECH = 4 that 125%-vs-14% asymmetry becomes a
#     ~15-fold difference in fatigue dose per exercising beat.  So the RV goes
#     first, and the RV/LV gap WIDENS with training -- both predictions, not
#     inputs.  Left dominance requires an extra, genotype-specific LV term
#     (KAPPA_LV), which is switched on only for DSP and FLNC.
#
#  3. THERE ARE TWO ARRHYTHMIA GENERATORS ON DIFFERENT CLOCKS, AND A DRUG
#     CANNOT OUTPERFORM THE GENERATOR IT OCCUPIES.
#       GENERATOR I  (early, catecholaminergic, structure-independent):
#         PKP2 loss strips Nav1.5 and Cx43 from the intercalated disc and
#         destabilises RyR2, producing diastolic Ca leak and triggered
#         activity in a structurally normal heart.
#         [Sato 2009 PMID 19661460; Cerrone 2014 PMID 24352520;
#          van Opbergen 2019 PMID 31438494]
#         Occupied by: beta-blockade, flecainide (RyR2 open-state block),
#         exercise restriction (removes the catecholamine surges).
#       GENERATOR II (late, scar-dependent, re-entrant):
#         patchy fibrofatty replacement, maximal when viable and replaced
#         tissue are interdigitated ~50/50, with slow conduction shortening
#         the re-entrant wavelength.  Subepicardial before endocardial.
#         [Basso 2009; Tschabrunn 2022 PMID 34883271]
#         Occupied by: epicardial ablation, amiodarone; NOT by beta-blockade.
#     This is what makes flecainide's sign depend on where the patient is:
#     the same drug removes Generator I trigger and slows conduction, which
#     FEEDS Generator II.  And it is why the ICD -- attached to neither
#     generator but to the outcome -- is the only unambiguous mortality
#     therapy.
#
#  ---------------------------------------------------------------------------
#  CALIBRATION: FITTED vs PREDICTED
#  ---------------------------------------------------------------------------
#  FITTED (4 parameters, all listed in FITTED_PARAMS below):
#    K_INJ    fatigue-failure scale     -> median age at definite diagnosis
#                                          ~36 y in PKP2 carriers on ordinary
#                                          recreational activity [PMID 25820315]
#    H0_VA    arrhythmia hazard scale   -> ~10%/yr sustained VA in definite ARVC
#    LAM2     Generator II weight       -> the amiodarone > sotalol gap
#                                          [PMID 19660690]
#    K_DIL    dilatation gain           -> RVEDVi trajectory to overt disease
#
#  PREDICTED (not fitted -- these are the model's exposed neck):
#    * the ~3x exercise hazard ratio and its dose-response (James 2013)
#    * incomplete penetrance in sedentary carriers (Carruth 2019)
#    * gene-elusive ARVC in extreme-dose athletes only (Sawant 2014)
#    * RV before LV with no chamber-specific biology, and RV/LV widening with
#      training (La Gerche 2011)
#    * load-reducing therapy preventing the phenotype (Fabritz 2011)
#    * the sotalol null (Marcus 2009)
#    * flecainide helping early and hurting late (same drug, opposite sign)
#    * endocardial-only vs endo-epicardial ablation recurrence gap
#      (Santangeli 2015, Berruezo 2012)
#    * ICD: full mortality benefit, exactly zero substrate benefit
#    * AAV-PKP2 gene therapy: large benefit given in the concealed phase,
#      small benefit once fibrofatty replacement is established
#      (Bradford 2023 PMID 39196150)
#
#  A STATED MISS: the model puts beta-blocker monotherapy at a larger
#  arrhythmia reduction than the North American registry observed (which found
#  no significant benefit for beta-blockers on VT endpoints).  See section 9.
#
#  A SECOND STATED MISS -- WHICH TURNS OUT TO BE THE SAME MISS.
#  Ermakov 2017 (PMID 27939893) found flecainide ADDED to a beta-blocker
#  reduced VT.  This model makes that combination slightly WORSE at every age
#  (+3% to +7%), because beta-blockade has already driven Generator I to zero,
#  so all that is left of flecainide is its conduction cost feeding Generator
#  II.  Note that this is the FIRST miss seen from the other side: both say
#  that part of the RyR2 leak must be adrenergic-INDEPENDENT.  Add such a
#  component and beta-blockade stops abolishing Generator I (repairing miss 1)
#  while flecainide-on-top-of-beta-blockade regains something to act on
#  (repairing miss 2).  Two independent discrepancies converging on one
#  missing term, with one testable repair, is the most useful thing this model
#  produced about its own structure.
#
#  DISCLAIMER: educational / research QSP model.  Semi-quantitative, not
#  independently validated, NOT for clinical decision-making.
# =============================================================================

import math
import sys

# =============================================================================
#  SECTION 1 -- PARAMETERS
# =============================================================================

FITTED_PARAMS = ("K_INJ", "H0_VA", "LAM2", "K_DIL")

P = dict(
    # ---- time base ----------------------------------------------------------
    AGE0=12.0,          # y      simulation starts in early adolescence
    DT=1.0,             # d      RK4 step (convergence checked against 0.25 d)

    # ---- mechanical / fatigue-failure core (COMMITMENT 1) -------------------
    NMECH=4.0,          # -      fatigue exponent, fixed a priori (Basquin range)
    PHI_EX=1.0,         # -      1 = load-dependent damage; 0 = FALSIFIER
    HR_REST=62.0,       # bpm
    HR_EX=158.0,        # bpm    vigorous-intensity heart rate
    MET_VIG=8.0,        # MET    intensity credited as "vigorous"
    K_EX_RV=1.25,       # -      fractional rise in RV ESWS at max exercise
    K_EX_LV=0.14,       # -      fractional rise in LV ESWS at max exercise
    GAMMA_G=2.0,        # -      how steeply reserve loss raises failure rate
    K_INJ=2.80e-6,      # /d     FITTED fatigue-failure scale
    FRAILTY=1.0,        # -      between-subject multiplier on K_INJ (see
                        #        cohort_penetrance: reserve, load history and
                        #        hormonal milieu are not identical between
                        #        carriers, and penetrance is a POPULATION
                        #        quantity, not a single trajectory)
    SIGMA_FRAILTY=0.55, # -      log-normal SD of FRAILTY across carriers
    SEX_K_FEMALE=0.72,  # -      female multiplier on the fatigue rate; a
                        #        PREDICTION of the observed male predominance
                        #        [Akdis 2017 PMID 28329361]
    K_REG=1.0e-6,       # /d     adult myocyte regeneration (essentially nil)
    MYO_FLOOR=0.38,     # -      REGIONAL SPARING.  This is a lumped free-wall
                        #        model, but the disease is regional: the
                        #        subtricuspid basal wall, apex and infundibulum
                        #        carry the highest local stress and go first,
                        #        while septum-adjacent and outflow myocardium
                        #        stays subcritical.  A lumped average therefore
                        #        cannot go to zero, and it does not in
                        #        pathological series.  Injury acts on
                        #        (MYO - MYO_FLOOR), which also caps the
                        #        autocatalytic runaway without introducing a
                        #        phase-transition parameter.
                        #        [Corrado 1997 PMID 9362410]
    KAPPA_LV=1.0,       # -      genotype-specific LV vulnerability (DSP/FLNC>1)
    XI_FIB=0.70,        # -      load-bearing contribution of fibrous tissue
    XI_FAT=0.35,        # -      bulk (not contractile) contribution of fat

    # ---- desmosome / genotype ----------------------------------------------
    PKP2_SET=1.0,       # -      set by genotype()
    DSP_SET=1.0,
    TAU_PKP2=30.0,      # d
    K_PG=0.06, K_PG_OFF=0.05,      # nuclear plakoglobin turnover
    KM_PG=0.55,
    K_W=0.09, K_W_OFF=0.09,        # canonical Wnt/beta-catenin
    K_H=0.05, K_H2=0.010, K_H_OFF=0.06,   # Hippo activation
    ETA_NAV=0.55,       # -      exponent PKP2 -> Nav1.5 at the ID
    ETA_CX=0.75,        # -      exponent PKP2 -> Cx43 at the ID
    ECV=0.30,           # -      exponent (gNa x coupling) -> conduction velocity
    TAU_ID=14.0,        # d

    # ---- Generator I: Ca handling / triggered activity ----------------------
    K_RY=0.10, K_RY_OFF=0.10,
    BASE_RY=0.10,       # -      baseline RyR2 lability in WT
    RY_GAIN=0.30,       # -      extra lability per unit PKP2 lost.  Kept
                        #        deliberately modest: PKP2 haploinsufficiency
                        #        raises diastolic leak by a factor of about
                        #        2-3, not by an order of magnitude
                        #        [van Opbergen 2019 PMID 31438494]
    K_UP=1.0, K_REL=0.22, K_EXT=1.4,
    PVC_MAX=4000.0,     # /24h   maximal PVC burden
    PVC_BASE=40.0,      # /24h   normal adult PVC count
    K50_TRIG=0.13, HILL_TRIG=4.0,
    TAU_PVC=5.0,        # d
    W_SOURCE_SINK=6.0,  # -      myocyte loss removes the electrotonic
                        #        (source-sink) suppression of ectopy.  This is
                        #        the dominant reason PVC burden explodes in
                        #        overt disease while a concealed carrier stays
                        #        near normal.

    # ---- Generator II: scar re-entry ---------------------------------------
    CV0=55.0,           # cm/s   healthy RV free-wall conduction velocity
    KM_CVFIB=0.55,      # -      fibrofatty tortuosity half-effect
    FF_THRESH=0.05,     # -      MINIMUM excess fibrofatty replacement that can
                        #        host a re-entrant circuit.  Re-entry needs a
                        #        path longer than the wavelength; a few percent
                        #        of diffuse interstitial fibrosis cannot supply
                        #        one, which is why ordinary age-related
                        #        fibrosis is not arrhythmogenic and why
                        #        Generator II has a genuine onset rather than
                        #        rising from zero with the first lost myocyte.
    TAU_CV=20.0,        # d
    N_CV=1.5,           # -      how strongly slow conduction feeds re-entry
    TAU_HET=60.0,       # d
    LAM2=1.91,          # -      FITTED Generator II weight in the hazard

    # ---- inflammation / fibro-fatty replacement -----------------------------
    K_NEC=0.14,         # /d     clearance of the acutely injured pool
    K_INF_ON=2.6, K_INF_OFF=0.05,
    IMM_SENS=0.55,      # -      NF-kB priming from desmosome loss
    PHI_INF=1.30,       # -      inflammation amplifies myocyte injury
    K_FAP=0.020, K_FAP_OFF=0.030,
    W_HIPPO=0.55, W_WNT=0.85,
    K_FIB=0.075, K_FIB_DEG=0.00035,
    K_FAT=0.115, K_FAT_DEG=0.00025,

    # ---- chamber mechanics --------------------------------------------------
    RVEDVI0=88.0,       # mL/m2  normal indexed RV end-diastolic volume
    LVEDVI0=68.0,       # mL/m2
    EF0_RV=0.55, EF0_LV=0.62,
    K_DIL=1.55,         # -      FITTED dilatation gain per unit function lost
    K_ATH=2.2,          # -      physiological athlete's RV dilatation
    TAU_CONT=45.0,      # d
    TAU_VOL=60.0,       # d
    K_PRV=0.55,         # -      RV systolic pressure rise as the RV fails
    W_INF_CONT=0.30,    # -      reversible stunning during a hot phase

    # ---- neurohormonal -----------------------------------------------------
    K_SNS=1.25, K_SNS2=0.85,
    TAU_SNS=20.0,
    KM_B1=1.10, TAU_B1=25.0,
    NT0=45.0,           # pg/mL
    K_BNP_RV=520.0, K_BNP_LV=900.0, TAU_BNP=10.0,

    # ---- event hazards -----------------------------------------------------
    H0_VA=0.2010,       # /y     FITTED sustained-VA hazard scale
    HMAX_VA=0.35,       # /y     saturation ceiling on the annual hazard.  A
                        #        hazard has to saturate: a patient cannot have
                        #        an unbounded number of sustained episodes per
                        #        year before an intervention or a death
                        #        intervenes.  h = HMAX*h_raw/(h_raw+HMAX),
                        #        which is indistinguishable from h_raw while
                        #        h_raw is small.
    P1=1.15, P2=1.30,   # -      generator exponents in the hazard
    P_SCD=0.30,         # -      P(sudden death | sustained VA, no ICD)
    ICD_EFF=0.88,       # -      P(VA terminated | ICD in situ)
    H_HF0=0.0022,       # /y     baseline HF-death hazard
    K_HF=0.115,         # /y     HF-death hazard gain on biventricular failure
    H_ICD_CX=0.055,     # /y     ICD lead/pocket complication hazard
    R_ICD_INAPPROP=0.042,  # /y  inappropriate shock rate

    # ---- drug PK (steady-state, chronic oral) --------------------------------
    #   nadolol
    BB_DOSE=80.0, BB_F=0.30, BB_CL=200.0, BB_MW=309.4, BB_FU=0.70, BB_KD=4.0,
    #   flecainide
    FL_DOSE=200.0, FL_F=0.90, FL_CL=600.0, FL_MW=414.3, FL_FU=0.60,
    FL_IC50_NA=6.0e3, FL_IC50_RYR=2.5e3,      # nM (use-dependent Na; RyR2)
    #   sotalol
    SO_DOSE=320.0, SO_F=0.95, SO_CL=150.0, SO_MW=272.4, SO_FU=1.00,
    SO_KD_B1=320.0, SO_IC50_IKR=8.0e3,
    #   amiodarone (2 compartments, slow enough to integrate directly)
    AM_DOSE=200.0, AM_F=0.50, AM_MW=645.3,
    AM_VC=60.0, AM_VP=4600.0, AM_CLD=48.0, AM_CL=110.0,
    AM_IC50_G1=1.2e3, AM_IC50_G2=1.6e3,
    AM_KTOX=1.0/3000.0,
    #   load-reducing / neurohormonal arms (fractional effects at target dose)
    E_LOADRED_MAX=0.26,     # diuretic + nitrate, or ARNI: RVEDV target
    E_MRA_FIB=0.30,         # spironolactone on fibrogenesis
    E_IL1_MAX=0.70,         # anakinra / canakinumab on INF decay
    E_GC_MAX=0.45,          # corticosteroid on INF decay
    E_GSK_MAX=0.60,         # GSK-3beta inhibition on the adipogenic switch
    #   AAV9-PKP2 gene therapy
    VEC_TD_MAX=0.62,        # -   transduced cardiomyocyte fraction at high dose
    TAU_TD=9.0,             # d   transduction/uncoating
    TAU_TG=21.0,            # d   transgene expression rise
    TG_GAIN=0.85,           # -   PKP2 restored per unit transduced fraction
    TG_HALF=3650.0,         # d   episomal transgene loss half-life

    # ---- ablation ----------------------------------------------------------
    ABL_ENDO=0.35,      # -      Generator II substrate removed, endo only
    ABL_EPI=0.78,       # -      combined endo-epicardial
    TAU_ABL=2600.0,     # d      substrate re-accumulation after ablation

    # ---- exposures set per scenario ----------------------------------------
    EX=15.0,            # MET-h/wk
    BSA=1.90,           # m2
    MALE=1.0,
    ON_BB=0.0, ON_FL=0.0, ON_SO=0.0, ON_AM=0.0,
    ON_LOADRED=0.0, ON_MRA=0.0, ON_IL1=0.0, ON_GC=0.0, ON_GSK=0.0,
    ON_ICD=0.0,
    BIOPSY=0.0,
)

# genotype presets: (PKP2 set-point, DSP set-point, LV vulnerability, variant?)
GENOTYPES = {
    "WT":     dict(PKP2_SET=1.00, DSP_SET=1.00, KAPPA_LV=1.00, VARIANT=0),
    "PKP2tv": dict(PKP2_SET=0.50, DSP_SET=1.00, KAPPA_LV=1.00, VARIANT=1),
    #   KAPPA_LV has to exceed LOAD_RV/LOAD_LV (about 3.4 in a competitive
    #   athlete) for the LV to actually go first.  Anything less makes DSP a
    #   slightly-less-RV-dominant disease, which is not what DSP is.
    "DSPtv":  dict(PKP2_SET=0.88, DSP_SET=0.50, KAPPA_LV=5.00, VARIANT=1),
    "FLNCtv": dict(PKP2_SET=0.95, DSP_SET=0.95, KAPPA_LV=4.60, VARIANT=1),
}

# exercise presets (MET-hours per week of vigorous-equivalent activity)
EXERCISE = {
    "sedentary":    6.0,    # below any guideline target
    "guideline":   15.0,    # ~150 min/wk moderate-to-vigorous
    "competitive": 60.0,    # ~8 h/wk endurance training
    "elite":      100.0,    # professional endurance load
}

STATES = [
    "PKP2_ID", "DSP_ID", "PG_NUC", "WNT_ACT", "HIPPO",
    "NAV_ID", "CX43_ID",
    "RYR_LEAK", "CA_SR", "CA_DIA",
    "D_RV", "MYO_RV", "NEC_RV", "INF_RV", "FAP_RV", "FIB_RV", "FAT_RV",
    "D_LV", "MYO_LV", "NEC_LV", "INF_LV", "FAP_LV", "FIB_LV", "FAT_LV",
    "RVEDV", "RV_CONT", "LVEDV", "LV_CONT",
    "CV_RV", "CV_LV", "SCAR_HET_RV", "SCAR_HET_LV",
    "PVC24", "SNS", "BETA1_D", "NTBNP", "ABL_HOM",
    "H_VT", "H_DEATH", "H_HF",
    "A_AM_C", "A_AM_P", "AMIO_TOX",
    "VEC_TD", "PKP2_TG",
    "ICD_SHK_A", "ICD_SHK_I",
]
IDX = {s: i for i, s in enumerate(STATES)}
NSTATE = len(STATES)

# normal myocardial composition of the RV free wall (fractions summing to 1)
MYO0, FIB0, FAT0 = 0.92, 0.03, 0.05
H_EFF0 = MYO0 + P["XI_FIB"] * FIB0 + P["XI_FAT"] * FAT0


def initial_state(p):
    y = [0.0] * NSTATE
    y[IDX["PKP2_ID"]] = p["PKP2_SET"]
    y[IDX["DSP_ID"]] = p["DSP_SET"]
    y[IDX["PG_NUC"]] = 0.0
    y[IDX["WNT_ACT"]] = 1.0
    y[IDX["HIPPO"]] = 0.0
    y[IDX["NAV_ID"]] = p["PKP2_SET"] ** p["ETA_NAV"]
    y[IDX["CX43_ID"]] = p["PKP2_SET"] ** p["ETA_CX"]
    y[IDX["RYR_LEAK"]] = 0.0
    y[IDX["CA_SR"]] = 1.0
    y[IDX["CA_DIA"]] = 0.10
    y[IDX["MYO_RV"]] = MYO0
    y[IDX["FIB_RV"]] = FIB0
    y[IDX["FAT_RV"]] = FAT0
    y[IDX["MYO_LV"]] = MYO0
    y[IDX["FIB_LV"]] = FIB0
    y[IDX["FAT_LV"]] = FAT0
    y[IDX["RVEDV"]] = p["RVEDVI0"] * p["BSA"]
    y[IDX["RV_CONT"]] = 1.0
    y[IDX["LVEDV"]] = p["LVEDVI0"] * p["BSA"]
    y[IDX["LV_CONT"]] = 1.0
    y[IDX["CV_RV"]] = p["CV0"] * ((p["PKP2_SET"] ** p["ETA_NAV"])
                                  * (p["PKP2_SET"] ** p["ETA_CX"])) ** p["ECV"]
    y[IDX["CV_LV"]] = p["CV0"] * 1.05
    y[IDX["PVC24"]] = p["PVC_BASE"]
    y[IDX["SNS"]] = 1.0
    y[IDX["BETA1_D"]] = 1.0
    y[IDX["NTBNP"]] = p["NT0"]
    # give the fast ID/Ca compartments their algebraic steady state so the
    # trajectory does not start with a spurious transient
    y[IDX["RYR_LEAK"]] = p["RY_GAIN"] * (1.0 - p["PKP2_SET"]) + p["BASE_RY"]
    y[IDX["CA_DIA"]] = p["K_REL"] * y[IDX["RYR_LEAK"]] * 1.0 / p["K_EXT"]
    return y


# =============================================================================
#  SECTION 2 -- STEADY-STATE PK FOR THE CHRONIC ORAL DRUGS
# =============================================================================
#  The disease clock in this model runs for decades.  The fast oral drugs
#  (nadolol t1/2 ~20 h, flecainide ~14 h, sotalol ~12 h) reach steady state in
#  days, so integrating their absorption/distribution ODEs at a 0.5-day step
#  would be both unstable and pointless.  Here they are represented by the
#  EXACT average steady-state concentration for the same linear model that
#  arvc_mrgsolve_model.R integrates:
#
#      Css,avg = F * Dose / (CL * tau)
#
#  --check verifies this against a fine-step numerical integration of the same
#  one-compartment oral ODEs (section 10, check "pk_css_matches_ode"), so the
#  simplification is demonstrated rather than asserted.
#
#  Amiodarone is different -- a ~58-day terminal half-life and a huge
#  peripheral volume mean its loading kinetics matter on the disease timescale,
#  so it keeps real ODE compartments.
# =============================================================================

def css_free_nM(dose_mg_per_day, F, CL_L_per_day, MW, fu):
    """Average free steady-state concentration in nM for a chronic oral drug."""
    css_mg_per_L = F * dose_mg_per_day / CL_L_per_day
    return fu * css_mg_per_L / MW * 1e6


def pk_effects(y, p):
    """All drug/intervention effects, as dimensionless fractions in [0,1]."""
    e = {}

    # --- beta-blockade (nadolol), Generator I via adrenergic drive ----------
    c_bb = p["ON_BB"] * css_free_nM(p["BB_DOSE"], p["BB_F"], p["BB_CL"],
                                    p["BB_MW"], p["BB_FU"])
    occ_bb = c_bb / (c_bb + p["BB_KD"]) if c_bb > 0 else 0.0

    # --- sotalol: weak beta-block PLUS IKr block ---------------------------
    c_so = p["ON_SO"] * css_free_nM(p["SO_DOSE"], p["SO_F"], p["SO_CL"],
                                    p["SO_MW"], p["SO_FU"])
    occ_so_b1 = c_so / (c_so + p["SO_KD_B1"]) if c_so > 0 else 0.0
    occ_so_ikr = c_so / (c_so + p["SO_IC50_IKR"]) if c_so > 0 else 0.0
    e["OCC_B1"] = 1.0 - (1.0 - occ_bb) * (1.0 - occ_so_b1)

    # --- flecainide: RyR2 open-state block (helps I) AND Na block (feeds II)
    c_fl = p["ON_FL"] * css_free_nM(p["FL_DOSE"], p["FL_F"], p["FL_CL"],
                                    p["FL_MW"], p["FL_FU"])
    e["E_FL_RYR"] = c_fl / (c_fl + p["FL_IC50_RYR"]) if c_fl > 0 else 0.0
    e["E_NA_BLOCK"] = c_fl / (c_fl + p["FL_IC50_NA"]) if c_fl > 0 else 0.0

    # --- amiodarone from its own compartments ------------------------------
    c_am = y[IDX["A_AM_C"]] / p["AM_VC"] / p["AM_MW"] * 1e6   # nM, total
    e["E_AM_G1"] = c_am / (c_am + p["AM_IC50_G1"])
    e["E_AM_G2"] = c_am / (c_am + p["AM_IC50_G2"])
    e["C_AM_nM"] = c_am

    # --- repolarisation lengthening: helps re-entry, hurts triggering -------
    #   IKr block widens the re-entrant wavelength (anti-Generator-II) but
    #   promotes early afterdepolarisations in diseased tissue (pro-Gen-I).
    e["E_WL"] = 1.0 - (1.0 - 0.62 * occ_so_ikr) * (1.0 - 0.55 * e["E_AM_G2"])
    #   The proarrhythmic arm is deliberately SMALL relative to the
    #   wavelength benefit, and much smaller for amiodarone than for sotalol:
    #   amiodarone prolongs QT markedly yet carries a low torsades rate, which
    #   is one of the better-established asymmetries in antiarrhythmic
    #   pharmacology and is not something this model gets to choose.
    e["PROARR_REPOL"] = 0.18 * occ_so_ikr + 0.08 * e["E_AM_G2"]
    #   And the more important arm: IKr block in a HETEROGENEOUS substrate
    #   increases dispersion of repolarisation, which promotes re-entry
    #   INITIATION.  Without this term a lengthened wavelength makes sotalol
    #   look good in exactly the patients in whom it demonstrably is not, so
    #   the term is what turns the drug's two arms into the observed null.
    #   The magnitude is set by a harder piece of evidence than the ARVC
    #   registry: d-sotalol INCREASED mortality in patients with LV dysfunction
    #   after myocardial infarction (SWORD, Waldo 1996 PMID 8691967).  A model
    #   in which IKr block is net-beneficial in a scarred ventricle is
    #   contradicted by a randomised trial, so the coefficient is chosen large
    #   enough that the two arms at least cancel.
    e["DISP_REPOL"] = 0.80 * occ_so_ikr + 0.10 * e["E_AM_G2"]

    # --- structural / upstream arms ----------------------------------------
    e["E_LOADRED"] = p["ON_LOADRED"] * p["E_LOADRED_MAX"]
    e["E_MRA"] = p["ON_MRA"] * p["E_MRA_FIB"]
    e["E_ANTIINF"] = 1.0 - (1.0 - p["ON_IL1"] * p["E_IL1_MAX"]) \
        * (1.0 - p["ON_GC"] * p["E_GC_MAX"])
    e["E_GSK"] = p["ON_GSK"] * p["E_GSK_MAX"]
    return e


# =============================================================================
#  SECTION 3 -- MECHANICS: THE LAPLACE ENGINE (COMMITMENT 2)
# =============================================================================

def mechanics(y, p):
    """Relative wall stress and the fatigue load index for each chamber."""
    m = {}

    # ---- effective (load-bearing) wall thickness --------------------------
    #   Myocardium bears load; fibrous replacement bears some; fat bears none.
    #   This is the autocatalytic link: losing myocytes thins the wall, which
    #   raises stress, which accelerates losing myocytes.
    h_rv = (y[IDX["MYO_RV"]] + p["XI_FIB"] * y[IDX["FIB_RV"]]
            + p["XI_FAT"] * y[IDX["FAT_RV"]]) / H_EFF0
    h_lv = (y[IDX["MYO_LV"]] + p["XI_FIB"] * y[IDX["FIB_LV"]]
            + p["XI_FAT"] * y[IDX["FAT_LV"]]) / H_EFF0
    h_rv = max(h_rv, 0.22)
    h_lv = max(h_lv, 0.22)

    # ---- chamber radius from volume (sphere-equivalent) -------------------
    r_rv = max(1e-6, y[IDX["RVEDV"]] / (p["RVEDVI0"] * p["BSA"])) ** (1.0 / 3.0)
    r_lv = max(1e-6, y[IDX["LVEDV"]] / (p["LVEDVI0"] * p["BSA"])) ** (1.0 / 3.0)

    # ---- chamber pressure: rises as the failing RV meets its afterload ----
    p_rv = 1.0 + p["K_PRV"] * max(0.0, r_rv ** 3 - 1.0)
    p_lv = 1.0

    # ---- Laplace, normalised to 1.0 in a healthy resting chamber ----------
    m["SIG_RV"] = p_rv * r_rv / h_rv
    m["SIG_LV"] = p_lv * r_lv / h_lv

    # ---- duty cycle of vigorous exercise ---------------------------------
    f_ex = min(0.35, p["EX"] / (p["MET_VIG"] * 168.0))
    m["F_EX"] = f_ex

    n = p["NMECH"]
    phi = p["PHI_EX"]

    def load(sig, k_ex):
        # mean-field fatigue dose rate: beats x (stress)^n, averaged over the
        # week, normalised so a perfectly rested healthy chamber scores 1.0.
        # PHI_EX = 0 removes the stress dependence entirely (the falsifier).
        s_rest = sig ** n if phi else 1.0
        s_ex = (sig * (1.0 + phi * k_ex)) ** n if phi else 1.0
        return ((1.0 - f_ex) * p["HR_REST"] * s_rest
                + f_ex * p["HR_EX"] * s_ex) / p["HR_REST"]

    m["LOAD_RV"] = load(m["SIG_RV"], p["K_EX_RV"])
    m["LOAD_LV"] = load(m["SIG_LV"], p["K_EX_LV"])
    return m


# =============================================================================
#  SECTION 4 -- RIGHT-HAND SIDE
# =============================================================================

def rhs(t, y, p, ev):
    d = [0.0] * NSTATE
    e = pk_effects(y, p)
    m = mechanics(y, p)

    # ---------------------------------------------------------------- desmosome
    # PKP2 protein at the ID relaxes to its genotype set point; AAV-delivered
    # transgene adds to it.
    pkp2_eff = min(1.0, y[IDX["PKP2_ID"]] + y[IDX["PKP2_TG"]])
    dsp_eff = y[IDX["DSP_ID"]]
    desmo = max(0.12, (pkp2_eff ** 0.7) * (dsp_eff ** 0.3))

    d[IDX["PKP2_ID"]] = (p["PKP2_SET"] - y[IDX["PKP2_ID"]]) / p["TAU_PKP2"]
    d[IDX["DSP_ID"]] = (p["DSP_SET"] - y[IDX["DSP_ID"]]) / p["TAU_PKP2"]

    # nuclear plakoglobin rises as the junction fails (the ARVC "signal-loss"
    # immunohistochemistry finding)
    d[IDX["PG_NUC"]] = p["K_PG"] * (1.0 - desmo) - p["K_PG_OFF"] * y[IDX["PG_NUC"]]

    # nuclear plakoglobin competes off beta-catenin -> canonical Wnt falls.
    # GSK-3beta inhibition restores it.
    wnt_target = 1.0 / (1.0 + (y[IDX["PG_NUC"]] / p["KM_PG"]) ** 2)
    wnt_target = wnt_target + e["E_GSK"] * (1.0 - wnt_target)
    d[IDX["WNT_ACT"]] = p["K_W"] * wnt_target - p["K_W_OFF"] * y[IDX["WNT_ACT"]]

    # Hippo activation: junction loss plus mechanical load
    d[IDX["HIPPO"]] = (p["K_H"] * (1.0 - desmo)
                       + p["K_H2"] * max(0.0, m["LOAD_RV"] - 1.0)
                       - p["K_H_OFF"] * y[IDX["HIPPO"]] * (1.0 + e["E_GSK"]))

    # ------------------------------------------------- ID electrical remodelling
    nav_t = pkp2_eff ** p["ETA_NAV"]
    cx_t = pkp2_eff ** p["ETA_CX"]
    d[IDX["NAV_ID"]] = (nav_t - y[IDX["NAV_ID"]]) / p["TAU_ID"]
    d[IDX["CX43_ID"]] = (cx_t - y[IDX["CX43_ID"]]) / p["TAU_ID"]

    # ------------------------------------------- GENERATOR I: Ca leak / triggers
    beta_act = y[IDX["SNS"]] * (1.0 - e["OCC_B1"]) * y[IDX["BETA1_D"]]
    ry_drive = (p["RY_GAIN"] * (1.0 - pkp2_eff) + p["BASE_RY"]) * beta_act \
        * (1.0 - e["E_FL_RYR"]) * (1.0 - 0.45 * e["E_AM_G1"])
    d[IDX["RYR_LEAK"]] = p["K_RY"] * ry_drive - p["K_RY_OFF"] * y[IDX["RYR_LEAK"]]

    d[IDX["CA_SR"]] = p["K_UP"] * (1.0 - y[IDX["CA_SR"]]) \
        - p["K_REL"] * y[IDX["RYR_LEAK"]] * y[IDX["CA_SR"]]
    d[IDX["CA_DIA"]] = p["K_REL"] * y[IDX["RYR_LEAK"]] * y[IDX["CA_SR"]] \
        - p["K_EXT"] * y[IDX["CA_DIA"]]

    # triggered-activity index: Ca leak, amplified as myocyte loss removes the
    # source-sink protection that normally suppresses ectopy
    myo_loss_rv = max(0.0, min(MYO0, MYO0 - y[IDX["MYO_RV"]])) / MYO0
    trig = y[IDX["CA_DIA"]] * (1.0 + p["W_SOURCE_SINK"] * myo_loss_rv) \
        * (1.0 + e["PROARR_REPOL"])
    hill = p["HILL_TRIG"]
    pvc_t = p["PVC_BASE"] + p["PVC_MAX"] * trig ** hill \
        / (p["K50_TRIG"] ** hill + trig ** hill)
    d[IDX["PVC24"]] = (pvc_t - y[IDX["PVC24"]]) / p["TAU_PVC"]

    # ------------------------------------------ COMMITMENT 1: the fatigue clock
    vuln = (1.0 / desmo) ** p["GAMMA_G"]
    d[IDX["D_RV"]] = m["LOAD_RV"] * vuln / 365.25
    d[IDX["D_LV"]] = m["LOAD_LV"] * vuln * p["KAPPA_LV"] / 365.25

    def chamber(tag, load, vuln_ch, kappa):
        iM, iN, iI, iP, iF, iA = (IDX["MYO_" + tag], IDX["NEC_" + tag],
                                  IDX["INF_" + tag], IDX["FAP_" + tag],
                                  IDX["FIB_" + tag], IDX["FAT_" + tag])
        myo = max(0.0, y[iM])
        nec, inf, fap = max(0.0, y[iN]), max(0.0, y[iI]), max(0.0, y[iP])
        fib, fat = max(0.0, y[iF]), max(0.0, y[iA])

        # fatigue-failure flux: beats x stress^n x (1/reserve)^gamma
        k_inj = p["K_INJ"] * p["FRAILTY"] \
            * (1.0 if p["MALE"] else p["SEX_K_FEMALE"])
        # injury acts on the SUSCEPTIBLE myocardium only -- see MYO_FLOOR
        suscept = max(0.0, myo - p["MYO_FLOOR"])
        inj = k_inj * load * vuln_ch * kappa * suscept \
            * (1.0 + p["PHI_INF"] * inf)
        space = max(0.0, (MYO0 - myo) - (fib - FIB0) - (fat - FAT0))

        d[iM] = -inj + p["K_REG"] * space
        d[iN] = inj - p["K_NEC"] * nec
        # myocyte death is immunogenic: NF-kB-driven inflammation both follows
        # injury and feeds back on it (the "hot phase")
        d[iI] = (p["K_INF_ON"] * nec * (1.0 + p["IMM_SENS"] * (1.0 - desmo))
                 - p["K_INF_OFF"] * inf * (1.0 + 2.2 * e["E_ANTIINF"]))
        # fibro-adipogenic progenitors: inflammation + Hippo + Wnt-suppression
        fap_drive = (inf + p["W_HIPPO"] * y[IDX["HIPPO"]]
                     + p["W_WNT"] * max(0.0, 1.0 - y[IDX["WNT_ACT"]]))
        d[iP] = p["K_FAP"] * fap_drive * (1.0 - fap) - p["K_FAP_OFF"] * fap
        # replacement only fills space the myocytes have vacated
        adipo = max(0.0, 1.0 - y[IDX["WNT_ACT"]]) * y[IDX["HIPPO"]]
        d[iF] = p["K_FIB"] * fap * space * (1.0 - e["E_MRA"]) - p["K_FIB_DEG"] * (fib - FIB0)
        d[iA] = p["K_FAT"] * fap * space * adipo - p["K_FAT_DEG"] * (fat - FAT0)
        return inj

    chamber("RV", m["LOAD_RV"], vuln, 1.0)
    chamber("LV", m["LOAD_LV"], vuln, p["KAPPA_LV"])

    # ------------------------------------------------------- chamber mechanics
    rv_cont_t = (max(0.0, y[IDX["MYO_RV"]]) / MYO0) ** 1.30 \
        * (1.0 - p["W_INF_CONT"] * min(1.0, y[IDX["INF_RV"]]))
    lv_cont_t = (max(0.0, y[IDX["MYO_LV"]]) / MYO0) ** 1.30 \
        * (1.0 - p["W_INF_CONT"] * min(1.0, y[IDX["INF_LV"]]))
    d[IDX["RV_CONT"]] = (rv_cont_t - y[IDX["RV_CONT"]]) / p["TAU_CONT"]
    d[IDX["LV_CONT"]] = (lv_cont_t - y[IDX["LV_CONT"]]) / p["TAU_CONT"]

    # dilatation: pathological (function loss) plus the physiological
    # athlete's RV.  Load-reducing therapy acts HERE -- and therefore, through
    # Laplace, on the clock itself.
    rvedv_t = p["RVEDVI0"] * p["BSA"] \
        * (1.0 + p["K_DIL"] * max(0.0, 1.0 - y[IDX["RV_CONT"]])) \
        * (1.0 + p["K_ATH"] * m["F_EX"]) * (1.0 - e["E_LOADRED"])
    lvedv_t = p["LVEDVI0"] * p["BSA"] \
        * (1.0 + 0.75 * p["K_DIL"] * max(0.0, 1.0 - y[IDX["LV_CONT"]])) \
        * (1.0 + 0.55 * p["K_ATH"] * m["F_EX"]) * (1.0 - 0.6 * e["E_LOADRED"])
    d[IDX["RVEDV"]] = (rvedv_t - y[IDX["RVEDV"]]) / p["TAU_VOL"]
    d[IDX["LVEDV"]] = (lvedv_t - y[IDX["LVEDV"]]) / p["TAU_VOL"]

    # --------------------------------------------- conduction and Generator II
    fibfat_rv = (y[IDX["FIB_RV"]] - FIB0) + (y[IDX["FAT_RV"]] - FAT0)
    fibfat_lv = (y[IDX["FIB_LV"]] - FIB0) + (y[IDX["FAT_LV"]] - FAT0)
    cv_rv_t = p["CV0"] * (max(1e-4, y[IDX["NAV_ID"]])
                          * max(1e-4, y[IDX["CX43_ID"]])) ** p["ECV"] \
        / (1.0 + max(0.0, fibfat_rv) / p["KM_CVFIB"]) \
        * (1.0 - 0.55 * e["E_NA_BLOCK"])
    cv_lv_t = p["CV0"] * 1.05 * (max(1e-4, y[IDX["NAV_ID"]])) ** p["ECV"] \
        / (1.0 + max(0.0, fibfat_lv) / p["KM_CVFIB"]) \
        * (1.0 - 0.55 * e["E_NA_BLOCK"])
    d[IDX["CV_RV"]] = (cv_rv_t - y[IDX["CV_RV"]]) / p["TAU_CV"]
    d[IDX["CV_LV"]] = (cv_lv_t - y[IDX["CV_LV"]]) / p["TAU_CV"]

    # scar heterogeneity: maximal when viable and replaced tissue interdigitate
    # ~50/50.  Pure myocardium has no channels; a homogeneous scar has none
    # either -- which is exactly what ablation exploits.
    het_rv = 4.0 * max(0.0, fibfat_rv) * (y[IDX["MYO_RV"]] / MYO0)
    het_lv = 4.0 * max(0.0, fibfat_lv) * (y[IDX["MYO_LV"]] / MYO0)
    d[IDX["SCAR_HET_RV"]] = (het_rv * (1.0 - y[IDX["ABL_HOM"]])
                             - y[IDX["SCAR_HET_RV"]]) / p["TAU_HET"]
    d[IDX["SCAR_HET_LV"]] = (het_lv - y[IDX["SCAR_HET_LV"]]) / p["TAU_HET"]
    d[IDX["ABL_HOM"]] = -y[IDX["ABL_HOM"]] / p["TAU_ABL"]

    # ------------------------------------------------------------ neurohormonal
    sns_t = 1.0 + p["K_SNS"] * max(0.0, 1.0 - y[IDX["RV_CONT"]]) \
        + p["K_SNS2"] * max(0.0, 1.0 - y[IDX["LV_CONT"]])
    d[IDX["SNS"]] = (sns_t - y[IDX["SNS"]]) / p["TAU_SNS"]
    b1_t = 1.0 / (1.0 + max(0.0, y[IDX["SNS"]] - 1.0) / p["KM_B1"])
    d[IDX["BETA1_D"]] = (b1_t - y[IDX["BETA1_D"]]) / p["TAU_B1"]

    nt_t = p["NT0"] * (1.0 + 0.6 * (y[IDX["SNS"]] - 1.0)) \
        + p["K_BNP_RV"] * max(0.0, y[IDX["RVEDV"]] / (p["RVEDVI0"] * p["BSA"]) - 1.0) \
        + p["K_BNP_LV"] * max(0.0, 1.0 - y[IDX["LV_CONT"]])
    d[IDX["NTBNP"]] = (nt_t - y[IDX["NTBNP"]]) / p["TAU_BNP"]

    # -------------------------------------------------------------- amiodarone
    am_in = p["ON_AM"] * p["AM_F"] * p["AM_DOSE"]           # mg/d
    d[IDX["A_AM_C"]] = (am_in
                        - p["AM_CL"] * y[IDX["A_AM_C"]] / p["AM_VC"]
                        - p["AM_CLD"] * (y[IDX["A_AM_C"]] / p["AM_VC"]
                                         - y[IDX["A_AM_P"]] / p["AM_VP"]))
    d[IDX["A_AM_P"]] = p["AM_CLD"] * (y[IDX["A_AM_C"]] / p["AM_VC"]
                                      - y[IDX["A_AM_P"]] / p["AM_VP"])
    d[IDX["AMIO_TOX"]] = p["AM_KTOX"] * e["E_AM_G1"]

    # ------------------------------------------------------- AAV gene transfer
    d[IDX["VEC_TD"]] = -y[IDX["VEC_TD"]] * 0.0        # set by dose events
    d[IDX["PKP2_TG"]] = (p["TG_GAIN"] * y[IDX["VEC_TD"]] - y[IDX["PKP2_TG"]]) \
        / p["TAU_TG"] - math.log(2.0) / p["TG_HALF"] * y[IDX["PKP2_TG"]]

    # -------------------------------------------------------------- the hazards
    gen1, gen2 = generators(y, p, e)
    h_raw = p["H0_VA"] * (gen1 ** p["P1"] + p["LAM2"] * gen2 ** p["P2"])
    h_va = p["HMAX_VA"] * h_raw / (h_raw + p["HMAX_VA"])
    d[IDX["H_VT"]] = h_va / 365.25

    biv = max(0.0, 1.0 - y[IDX["RV_CONT"]]) + 1.4 * max(0.0, 1.0 - y[IDX["LV_CONT"]])
    h_hf = p["H_HF0"] + p["K_HF"] * biv ** 2
    icd = p["ON_ICD"]
    h_death = (1.0 - icd * p["ICD_EFF"]) * p["P_SCD"] * h_va + h_hf \
        + icd * p["H_ICD_CX"] * 0.10
    d[IDX["H_HF"]] = h_hf / 365.25
    d[IDX["H_DEATH"]] = h_death / 365.25
    d[IDX["ICD_SHK_A"]] = icd * h_va / 365.25
    d[IDX["ICD_SHK_I"]] = icd * p["R_ICD_INAPPROP"] / 365.25
    return d


def generators(y, p, e):
    """The two arrhythmia generators (COMMITMENT 3), dimensionless."""
    # Generator I -- triggered activity, read out as EXCESS PVC burden over the
    # normal background.  A normal heart scores zero, so the model does not
    # attribute a hazard to ordinary ectopy.
    gen1 = max(0.0, y[IDX["PVC24"]] - p["PVC_BASE"]) / 1000.0
    # Generator II -- re-entry: patchy substrate above the circuit-size
    # threshold, times slow conduction, minus whatever the ablation
    # homogenised, and widened by IKr block
    cv = max(6.0, y[IDX["CV_RV"]])
    het = max(0.0, y[IDX["SCAR_HET_RV"]] - 4.0 * p["FF_THRESH"]) \
        + 0.45 * max(0.0, y[IDX["SCAR_HET_LV"]] - 4.0 * p["FF_THRESH"])
    gen2 = het * (p["CV0"] / cv) ** p["N_CV"] * (1.0 - 0.45 * e["E_WL"]) \
        * (1.0 + e["DISP_REPOL"]) * (1.0 - 0.30 * e["E_AM_G2"])
    return max(0.0, gen1), max(0.0, gen2)


# =============================================================================
#  SECTION 5 -- OBSERVABLES AND THE 2010 TASK FORCE CRITERIA
# =============================================================================

def observe(t, y, p):
    o = {}
    e = pk_effects(y, p)
    m = mechanics(y, p)
    o["age"] = p["AGE0"] + t / 365.25
    o["MYO_RV"] = y[IDX["MYO_RV"]]
    o["FIBFAT_RV"] = (y[IDX["FIB_RV"]] - FIB0) + (y[IDX["FAT_RV"]] - FAT0)
    o["FAT_RV"] = y[IDX["FAT_RV"]]
    o["FIB_RV"] = y[IDX["FIB_RV"]]
    o["MYO_LV"] = y[IDX["MYO_LV"]]
    o["FIBFAT_LV"] = (y[IDX["FIB_LV"]] - FIB0) + (y[IDX["FAT_LV"]] - FAT0)
    o["FIB_LV"] = y[IDX["FIB_LV"]]
    o["FAT_LV"] = y[IDX["FAT_LV"]]
    o["RVEF"] = p["EF0_RV"] * y[IDX["RV_CONT"]]
    o["LVEF"] = p["EF0_LV"] * y[IDX["LV_CONT"]]
    o["RVEDVI"] = y[IDX["RVEDV"]] / p["BSA"]
    o["LVEDVI"] = y[IDX["LVEDV"]] / p["BSA"]
    o["CV_RV"] = y[IDX["CV_RV"]]
    o["PVC24"] = y[IDX["PVC24"]]
    o["NTBNP"] = y[IDX["NTBNP"]]
    o["SIG_RV"] = m["SIG_RV"]
    o["SIG_LV"] = m["SIG_LV"]
    o["LOAD_RV"] = m["LOAD_RV"]
    o["LOAD_LV"] = m["LOAD_LV"]
    o["D_RV"] = y[IDX["D_RV"]]
    o["D_LV"] = y[IDX["D_LV"]]
    o["PKP2_EFF"] = min(1.0, y[IDX["PKP2_ID"]] + y[IDX["PKP2_TG"]])
    o["INF_RV"] = y[IDX["INF_RV"]]
    o["AMIO_TOX"] = y[IDX["AMIO_TOX"]]

    # residual myocytes as the pathologist would score the biopsy
    tot = max(0.0, y[IDX["MYO_RV"]]) + y[IDX["FIB_RV"]] + y[IDX["FAT_RV"]]
    o["RESID_MYO_PCT"] = 100.0 * max(0.0, y[IDX["MYO_RV"]]) / max(1e-9, tot)

    # terminal activation duration scales inversely with conduction velocity
    o["TAD_ms"] = 32.0 * p["CV0"] / max(6.0, y[IDX["CV_RV"]])
    # right precordial T-wave inversion extent tracks substrate extent
    o["N_TWI"] = min(4, int(3.0 * min(1.0, max(0.0, (MYO0 - y[IDX["MYO_RV"]]) / 0.30))))
    o["EPSILON"] = 1 if y[IDX["CV_RV"]] < 22.0 else 0
    o["LATE_POT"] = 1 if o["TAD_ms"] >= 50.0 else 0

    g1, g2 = generators(y, p, e)
    o["GEN1"] = g1
    o["GEN2"] = g2
    h_raw = p["H0_VA"] * (g1 ** p["P1"] + p["LAM2"] * g2 ** p["P2"])
    o["H_VA_RAW"] = h_raw
    o["H_VA_PER_YR"] = p["HMAX_VA"] * h_raw / (h_raw + p["HMAX_VA"])
    o["S_VA_5Y"] = math.exp(-o["H_VA_PER_YR"] * 5.0)
    o["H_VT"] = y[IDX["H_VT"]]
    o["H_DEATH"] = y[IDX["H_DEATH"]]
    o["SURV"] = math.exp(-y[IDX["H_DEATH"]])
    o["VA_FREE"] = math.exp(-y[IDX["H_VT"]])
    o["ICD_SHK_A"] = y[IDX["ICD_SHK_A"]]
    o["ICD_SHK_I"] = y[IDX["ICD_SHK_I"]]

    # ---- 2010 Task Force Criteria (Marcus 2010 PMID 20172912) --------------
    major, minor = [], []
    regional = y[IDX["MYO_RV"]] < 0.85 * MYO0          # wall-motion abnormality
    edvi_major = 110.0 if p["MALE"] else 100.0
    edvi_minor = 100.0 if p["MALE"] else 90.0
    if regional and (o["RVEF"] <= 0.40 or o["RVEDVI"] >= edvi_major):
        major.append("I:structure")
    elif regional and (0.40 < o["RVEF"] <= 0.45 or o["RVEDVI"] >= edvi_minor):
        minor.append("I:structure")
    if p["BIOPSY"]:
        if o["RESID_MYO_PCT"] < 60.0:
            major.append("II:tissue")
        elif o["RESID_MYO_PCT"] < 75.0:
            minor.append("II:tissue")
    if o["N_TWI"] >= 3:
        major.append("III:repol")
    elif o["N_TWI"] == 2:
        minor.append("III:repol")
    if o["EPSILON"]:
        major.append("IV:depol")
    elif o["TAD_ms"] >= 55.0 or o["LATE_POT"]:
        minor.append("IV:depol")
    # arrhythmia category: LV-morphology NSVT is major, RV-morphology minor
    nsvt = o["H_VA_PER_YR"] > 0.045
    if nsvt and y[IDX["MYO_LV"]] < 0.85 * MYO0:
        major.append("V:arrhythmia")
    elif nsvt or o["PVC24"] > 500.0:
        minor.append("V:arrhythmia")
    if p.get("VARIANT", 0):
        major.append("VI:genetics")

    nM, nm = len(major), len(minor)
    o["TFC_MAJOR"], o["TFC_MINOR"] = nM, nm
    if nM >= 2 or (nM == 1 and nm >= 2) or nm >= 4:
        o["TFC"] = "definite"
    elif (nM == 1 and nm == 1) or nm == 3:
        o["TFC"] = "borderline"
    elif nM == 1 or nm == 2:
        o["TFC"] = "possible"
    else:
        o["TFC"] = "none"
    o["TFC_LIST"] = "+".join(major + [x + "(m)" for x in minor])
    return o


# =============================================================================
#  SECTION 6 -- INTEGRATOR WITH SCHEDULED EVENTS
# =============================================================================
#  Events are (age_years, kind, value):
#     ("EX", v)        change exercise dose to v MET-h/wk
#     ("DRUG", name)   turn a drug on ("BB","FL","SO","AM","LOADRED","MRA",
#                      "IL1","GC","GSK") -- prefix "-" turns it off
#     ("ICD", 1/0)     implant / explant
#     ("ABL", "endo"|"epi")
#     ("AAV", frac)    single AAV9-PKP2 infusion, transduced fraction target
# =============================================================================

def simulate(params=None, genotype="PKP2tv", exercise="guideline",
             events=(), age_end=70.0, record_every=30.0, stop_when=None):
    p = dict(P)
    p.update(GENOTYPES[genotype])
    p["EX"] = EXERCISE[exercise] if isinstance(exercise, str) else float(exercise)
    if params:
        p.update(params)

    y = initial_state(p)
    t = 0.0
    dt = p["DT"]
    t_end = (age_end - p["AGE0"]) * 365.25
    ev = sorted([( (a - p["AGE0"]) * 365.25, k, v) for (a, k, v) in events])
    ei = 0
    out = []
    next_rec = 0.0

    def apply_event(kind, val):
        if kind == "EX":
            p["EX"] = EXERCISE[val] if isinstance(val, str) else float(val)
        elif kind == "DRUG":
            off = val.startswith("-")
            key = "ON_" + val.lstrip("-")
            p[key] = 0.0 if off else 1.0
        elif kind == "ICD":
            p["ON_ICD"] = float(val)
        elif kind == "BIOPSY":
            p["BIOPSY"] = float(val)
        elif kind == "ABL":
            frac = p["ABL_EPI"] if val == "epi" else p["ABL_ENDO"]
            y[IDX["ABL_HOM"]] = min(0.95, y[IDX["ABL_HOM"]] + frac)
            y[IDX["SCAR_HET_RV"]] *= (1.0 - frac)
        elif kind == "AAV":
            y[IDX["VEC_TD"]] = min(p["VEC_TD_MAX"], y[IDX["VEC_TD"]] + float(val))
        else:
            raise ValueError("unknown event " + kind)

    while t <= t_end + 1e-9:
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            apply_event(ev[ei][1], ev[ei][2])
            ei += 1
        if t >= next_rec - 1e-9:
            o_now = observe(t, y, p)
            out.append(o_now)
            next_rec += record_every
            if stop_when is not None and stop_when(o_now):
                return out, p
        # RK4
        k1 = rhs(t, y, p, ev)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NSTATE)]
        k2 = rhs(t + 0.5 * dt, y2, p, ev)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NSTATE)]
        k3 = rhs(t + 0.5 * dt, y3, p, ev)
        y4 = [y[i] + dt * k3[i] for i in range(NSTATE)]
        k4 = rhs(t + dt, y4, p, ev)
        for i in range(NSTATE):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        # hard physical bounds
        for nm2 in ("MYO_RV", "MYO_LV", "FIB_RV", "FIB_LV", "FAT_RV", "FAT_LV",
                    "WNT_ACT", "PKP2_ID", "DSP_ID", "PKP2_TG", "BETA1_D",
                    "RV_CONT", "LV_CONT", "CA_SR"):
            y[IDX[nm2]] = min(1.0, max(0.0, y[IDX[nm2]]))
        for nm2 in ("PG_NUC", "HIPPO", "INF_RV", "INF_LV", "FAP_RV", "FAP_LV",
                    "NEC_RV", "NEC_LV", "RYR_LEAK", "CA_DIA", "SCAR_HET_RV",
                    "SCAR_HET_LV", "PVC24", "A_AM_C", "A_AM_P"):
            y[IDX[nm2]] = max(0.0, y[IDX[nm2]])
        t += dt
    out.append(observe(t_end, y, p))
    return out, p


def age_at(traj, pred):
    for o in traj:
        if pred(o):
            return o["age"]
    return float("inf")


def at_age(traj, age):
    best = traj[0]
    for o in traj:
        if o["age"] <= age + 1e-6:
            best = o
        else:
            break
    return best


# =============================================================================
#  SECTION 6b -- PENETRANCE AS A POPULATION QUANTITY
# =============================================================================
#  A single trajectory cannot answer "is this variant penetrant?", because
#  penetrance is a statement about a cohort.  Carriers differ in desmosomal
#  reserve, in lifetime load history and in hormonal milieu, so the fatigue
#  scale K_INJ carries a between-subject multiplier FRAILTY.
#
#  The cohort is built from FIXED standard-normal quantiles rather than random
#  draws, so the penetrance curves are exactly reproducible and the file stays
#  dependency-free (no RNG, no seed to remember).
#
#  Nothing here is fitted.  SIGMA_FRAILTY = 0.55 is a plausible log-normal
#  spread chosen a priori; the female multiplier 0.72 is a single number
#  standing in for the observed male predominance.  What comes OUT -- roughly a
#  third of sedentary carriers penetrant by 60 against essentially all
#  competitive athletes penetrant by 40 -- is the prediction.
# =============================================================================

# standard-normal quantiles at (i+0.5)/n for n = 15, hard-coded so the cohort
# needs no statistics library
Z15 = (-1.834, -1.282, -0.967, -0.728, -0.524, -0.341, -0.168, 0.0,
       0.168, 0.341, 0.524, 0.728, 0.967, 1.282, 1.834)
Z9 = (-1.593, -1.036, -0.674, -0.319, 0.0, 0.319, 0.674, 1.036, 1.593)


def cohort_penetrance(genotype="PKP2tv", exercise="guideline", male=True,
                      ages=(30.0, 40.0, 50.0, 60.0), n=15, age_end=75.0,
                      params=None):
    """Fraction of a carrier cohort meeting definite TFC by each age."""
    zs = Z15 if n == 15 else Z9
    sigma = P["SIGMA_FRAILTY"]
    dx = []
    for z in zs:
        pr = dict(FRAILTY=math.exp(sigma * z), MALE=1.0 if male else 0.0,
                  DT=1.0)
        if params:
            pr.update(params)
        traj, _ = simulate(genotype=genotype, exercise=exercise, params=pr,
                           events=[(12.0, "BIOPSY", 1)], age_end=age_end,
                           record_every=180.0,
                           stop_when=lambda o: o["TFC"] == "definite")
        dx.append(age_at(traj, lambda o: o["TFC"] == "definite"))
    return {a: sum(1 for d in dx if d <= a) / float(len(dx)) for a in ages}, dx


def scen_penetrance():
    print("\n" + "=" * 78)
    print(" SCENARIO 1b -- PENETRANCE IS A POPULATION QUANTITY")
    print("=" * 78)
    print(" 15-subject cohort per stratum, log-normal spread on the fatigue")
    print(" scale (SD 0.55, fixed quantiles -- no RNG).  Percentage meeting")
    print(" definite Task Force Criteria by each age.\n")
    print("  %-9s %-6s %-13s %8s %8s %8s %8s" %
          ("genotype", "sex", "exercise", "by 30", "by 40", "by 50", "by 60"))
    for gt in ("PKP2tv", "WT"):
        for male in (True, False):
            for ex in ("sedentary", "guideline", "competitive"):
                pen, _ = cohort_penetrance(gt, ex, male)
                print("  %-9s %-6s %-13s %7.0f%% %7.0f%% %7.0f%% %7.0f%%" %
                      (gt, "male" if male else "female", ex,
                       100 * pen[30.0], 100 * pen[40.0],
                       100 * pen[50.0], 100 * pen[60.0]))
    print("\n  Observed anchors, none of them fitted: penetrance in")
    print("  genotype-positive relatives is incomplete and roughly 35-50% by")
    print("  mid-life [PMID 25820315], most carriers of ARVC-gene")
    print("  loss-of-function variants in an unselected biobank have no")
    print("  phenotype [PMID 31638835], male carriers do worse than female")
    print("  [PMID 28329361], and competitive athletes are the group in which")
    print("  penetrance approaches complete [PMID 23871885, 25896080].")


# =============================================================================
#  SECTION 7 -- FATIGUE-DOSE ARITHMETIC (shown explicitly, no ODE needed)
# =============================================================================

def load_index(ex_metwk, k_ex, nmech=None, phi=1.0, sig=1.0):
    p = P
    n = P["NMECH"] if nmech is None else nmech
    f = min(0.35, ex_metwk / (p["MET_VIG"] * 168.0))
    if not phi:
        return ((1 - f) * p["HR_REST"] + f * p["HR_EX"]) / p["HR_REST"]
    return ((1 - f) * p["HR_REST"] * sig ** n
            + f * p["HR_EX"] * (sig * (1 + k_ex)) ** n) / p["HR_REST"]


def print_load_table():
    print("\n" + "=" * 78)
    print(" SECTION 7 -- THE FATIGUE CLOCK, BEFORE ANY ODE IS SOLVED")
    print("=" * 78)
    print(" Mean-field fatigue dose rate = beats x (wall stress)^NMECH, weekly")
    print(" average, normalised to a perfectly rested healthy chamber = 1.00.")
    print(" NMECH = %.1f fixed a priori.  Exercise stress rise: RV +%.0f%%, "
          "LV +%.0f%% [PMID 21085033]" % (P["NMECH"], 100 * P["K_EX_RV"],
                                          100 * P["K_EX_LV"]))
    print("")
    print("  %-13s %8s   %8s %8s   %8s" % ("exercise", "MET-h/wk", "LOAD_RV",
                                            "LOAD_LV", "RV/LV"))
    base = None
    for nmv in ("sedentary", "guideline", "competitive", "elite"):
        ex = EXERCISE[nmv]
        lr = load_index(ex, P["K_EX_RV"])
        ll = load_index(ex, P["K_EX_LV"])
        if base is None:
            base = lr
        print("  %-13s %8.0f   %8.3f %8.3f   %8.2f" % (nmv, ex, lr, ll, lr / ll))
    lr_sed = load_index(EXERCISE["sedentary"], P["K_EX_RV"])
    lr_comp = load_index(EXERCISE["competitive"], P["K_EX_RV"])
    print("")
    print("  PREDICTED RV fatigue-rate ratio, competitive vs sedentary : %.2f"
          % (lr_comp / lr_sed))
    print("  OBSERVED  hazard ratio for VT/death, highest vs lowest")
    print("            exercise tertile in desmosomal carriers          : 3.16")
    print("            [James 2013 PMID 23871885]")
    print("  -> not fitted: NMECH was fixed at 4 and the two stress")
    print("     coefficients come straight from exercise CMR.")


# =============================================================================
#  SECTION 8 -- SCENARIOS
# =============================================================================

def scen_natural_history():
    print("\n" + "=" * 78)
    print(" SCENARIO 1 -- PENETRANCE IS SET BY LOAD, NOT BY GENOTYPE")
    print("=" * 78)
    print(" Age at which the 2010 Task Force Criteria first reach 'definite'.")
    print(" Same equations, same genotype; only the exercise dose differs.\n")
    print("  %-9s %-13s %8s %8s %8s %8s %8s" %
          ("genotype", "exercise", "age_dx", "RVEF@40", "RVEDVi", "PVC/24h", "%myo"))
    rows = {}
    for gt in ("WT", "PKP2tv"):
        for ex in ("sedentary", "guideline", "competitive", "elite"):
            traj, p = simulate(genotype=gt, exercise=ex,
                               events=[(12.0, "BIOPSY", 1)], age_end=80.0)
            adx = age_at(traj, lambda o: o["TFC"] == "definite")
            o40 = at_age(traj, 40.0)
            rows[(gt, ex)] = (adx, traj)
            print("  %-9s %-13s %8s %8.2f %8.0f %8.0f %8.1f" %
                  (gt, ex, ("%.1f" % adx) if adx < 1e6 else "never",
                   o40["RVEF"], o40["RVEDVI"], o40["PVC24"], o40["RESID_MYO_PCT"]))
    print("\n  Reading: a PKP2 truncating variant in a sedentary life is")
    print("  frequently a non-event (incomplete penetrance, PMID 31638835);")
    print("  the SAME variant in a competitive athlete is a disease of the")
    print("  twenties (PMID 23871885); and a variant-negative athlete can")
    print("  reach the phenotype on load alone if the load is extreme enough")
    print("  (gene-elusive ARVC, PMID 25516436).")
    return rows


def scen_chamber_selectivity():
    print("\n" + "=" * 78)
    print(" SCENARIO 2 -- WHY THE RIGHT VENTRICLE, WITH NO RV-SPECIFIC BIOLOGY")
    print("=" * 78)
    print("  %-9s %-13s %9s %9s %9s %9s %8s" %
          ("genotype", "exercise", "D_RV", "D_LV", "fibfatRV", "fibfatLV",
           "RV/LV"))
    for gt in ("PKP2tv", "DSPtv"):
        for ex in ("sedentary", "competitive"):
            traj, p = simulate(genotype=gt, exercise=ex, age_end=50.0)
            o = traj[-1]
            r = (o["FIBFAT_RV"] / o["FIBFAT_LV"]) if o["FIBFAT_LV"] > 1e-9 else float("inf")
            print("  %-9s %-13s %9.0f %9.0f %9.4f %9.4f %8.2f" %
                  (gt, ex, o["D_RV"], o["D_LV"], o["FIBFAT_RV"],
                   o["FIBFAT_LV"], r))
    print("\n  PKP2 is expressed equally in both ventricles and the model")
    print("  applies no chamber-specific biology to it.  The RV:LV ratio is")
    print("  pure Laplace + the exercise stress asymmetry, and it WIDENS with")
    print("  training.  Left dominance needs the genotype-specific KAPPA_LV")
    print("  term, which is on only for DSP/FLNC -- matching where")
    print("  left-dominant arrhythmogenic cardiomyopathy actually comes from.")


def scen_load_reduction():
    print("\n" + "=" * 78)
    print(" SCENARIO 3 -- LOAD-REDUCING THERAPY (the external validation)")
    print("=" * 78)
    print(" Fabritz 2011 (PMID 21292134) gave plakoglobin-deficient mice")
    print(" furosemide + nitrate -- no desmosomal, ion-channel or anti-fibrotic")
    print(" action whatever -- and prevented RV enlargement and the arrhythmic")
    print(" phenotype.  The model was NOT fitted to this.  Here the same")
    print(" intervention is a 26% reduction in the RV volume set point from")
    print(" age 16, in a competitive athlete carrying PKP2tv.\n")
    print("  %-34s %8s %8s %8s %8s" %
          ("arm", "age_dx", "RVEDVi", "RVEF", "VA/yr"))
    arms = [
        ("no therapy", []),
        ("exercise restriction only", [(16.0, "EX", "guideline")]),
        ("load-reducing therapy only", [(16.0, "DRUG", "LOADRED")]),
        ("both", [(16.0, "EX", "guideline"), (16.0, "DRUG", "LOADRED")]),
        ("beta-blocker only", [(16.0, "DRUG", "BB")]),
    ]
    res = {}
    for nm2, ev in arms:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=[(12.0, "BIOPSY", 1)] + ev, age_end=60.0)
        adx = age_at(traj, lambda o: o["TFC"] == "definite")
        o = at_age(traj, 40.0)
        res[nm2] = (adx, o)
        print("  %-34s %8s %8.0f %8.2f %8.3f" %
              (nm2, ("%.1f" % adx) if adx < 1e6 else "never",
               o["RVEDVI"], o["RVEF"], o["H_VA_PER_YR"]))
    print("\n  Note what beta-blockade does NOT do: it lowers the arrhythmia")
    print("  hazard without moving the structural trajectory, because it")
    print("  occupies Generator I and not the clock.")
    return res


def scen_antiarrhythmics():
    print("\n" + "=" * 78)
    print(" SCENARIO 4 -- A DRUG CANNOT OUTPERFORM THE GENERATOR IT OCCUPIES")
    print("=" * 78)
    print(" A PKP2tv patient with definite ARVC at age 38 (overt substrate).")
    print(" Each drug is started at 38 and the sustained-VA hazard is read at")
    print(" age 39.  Nothing here is fitted except the Generator II weight.\n")
    base_ev = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline")]
    ref = None
    print("  %-30s %9s %9s %9s %9s" %
          ("therapy", "GEN1", "GEN2", "VA/yr", "vs none"))
    out = {}
    for nm2, ev in [
        ("none", []),
        ("beta-blocker (nadolol)", [(38.0, "DRUG", "BB")]),
        ("sotalol", [(38.0, "DRUG", "SO")]),
        ("flecainide alone", [(38.0, "DRUG", "FL")]),
        ("beta-blocker + flecainide", [(38.0, "DRUG", "BB"), (38.0, "DRUG", "FL")]),
        ("amiodarone", [(38.0, "DRUG", "AM")]),
        ("ablation, endocardial only", [(38.0, "ABL", "endo")]),
        ("ablation, endo + epicardial", [(38.0, "ABL", "epi")]),
    ]:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=base_ev + ev, age_end=42.0, record_every=15.0)
        o = at_age(traj, 39.0)
        if ref is None:
            ref = o["H_VA_PER_YR"]
        out[nm2] = o
        print("  %-30s %9.2f %9.3f %9.3f %8.0f%%" %
              (nm2, o["GEN1"], o["GEN2"], o["H_VA_PER_YR"],
               100.0 * (o["H_VA_PER_YR"] / ref - 1.0)))
    print("\n  Observed anchors: amiodarone was the only drug with clear")
    print("  efficacy and sotalol was not better than nothing in the North")
    print("  American ARVC registry [PMID 19660690]; flecainide added to a")
    print("  beta-blocker reduced VT [PMID 27939893]; endocardial-only")
    print("  ablation recurs in most patients whereas combined")
    print("  endo-epicardial ablation does not [PMID 26546346, 22205683].")
    return out


def scen_flecainide_sign():
    print("\n" + "=" * 78)
    print(" SCENARIO 5 -- THE SAME DRUG WITH TWO SIGNS")
    print("=" * 78)
    print(" Flecainide blocks the RyR2 open state (removes Generator I) and")
    print(" blocks Nav1.5 (slows conduction, which FEEDS Generator II).  So")
    print(" its sign must depend on where the patient is on the trajectory.")
    print(" This is a prediction, not a fit.\n")
    print("  %6s %9s %9s %11s %11s %9s %11s" %
          ("age", "fibfat", "GEN1", "VA/yr base", "VA/yr flec", "change",
           "+BB change"))
    for age in (22.0, 28.0, 34.0, 40.0, 46.0):
        common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline")]
        t0, _ = simulate(genotype="PKP2tv", exercise="competitive",
                         events=common, age_end=age + 1.0, record_every=15.0)
        t1, _ = simulate(genotype="PKP2tv", exercise="competitive",
                         events=common + [(age, "DRUG", "FL")],
                         age_end=age + 1.0, record_every=15.0)
        t2, _ = simulate(genotype="PKP2tv", exercise="competitive",
                         events=common + [(16.0, "DRUG", "BB")],
                         age_end=age + 1.0, record_every=15.0)
        t3, _ = simulate(genotype="PKP2tv", exercise="competitive",
                         events=common + [(16.0, "DRUG", "BB"),
                                          (age, "DRUG", "FL")],
                         age_end=age + 1.0, record_every=15.0)
        a = at_age(t0, age + 0.9)
        b = at_age(t1, age + 0.9)
        c = at_age(t2, age + 0.9)
        d2 = at_age(t3, age + 0.9)
        print("  %6.0f %9.4f %9.2f %11.3f %11.3f %8.0f%% %10.0f%%" %
              (age, a["FIBFAT_RV"], a["GEN1"], a["H_VA_PER_YR"],
               b["H_VA_PER_YR"],
               100.0 * (b["H_VA_PER_YR"] / a["H_VA_PER_YR"] - 1.0),
               100.0 * (d2["H_VA_PER_YR"] / c["H_VA_PER_YR"] - 1.0)))
    print("\n  Early (trigger-dominant) flecainide helps; late")
    print("  (scar-dominant) the conduction cost overtakes the trigger")
    print("  benefit and the sign flips.  That is the CAST lesson arriving")
    print("  from the geometry rather than from a warning label.")


def scen_icd():
    print("\n" + "=" * 78)
    print(" SCENARIO 6 -- THE ICD: ALL OF THE MORTALITY, NONE OF THE DISEASE")
    print("=" * 78)
    print("  %-26s %10s %10s %10s %10s %9s" %
          ("arm", "RVEF@50", "fibfat@50", "VA/yr@50", "surv@60", "shocks"))
    for nm2, ev in [("no ICD", []), ("ICD at 38", [(38.0, "ICD", 1)])]:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=[(12.0, "BIOPSY", 1), (16.0, "EX", "guideline"),
                                   (16.0, "DRUG", "BB")] + ev, age_end=60.0)
        o50 = at_age(traj, 50.0)
        o60 = traj[-1]
        print("  %-26s %10.2f %10.4f %10.3f %10.3f %9.1f" %
              (nm2, o50["RVEF"], o50["FIBFAT_RV"], o50["H_VA_PER_YR"],
               o60["SURV"], o60["ICD_SHK_A"] + o60["ICD_SHK_I"]))
    print("\n  Identical myocardium, identical arrhythmia hazard, different")
    print("  survival.  The ICD is attached to the outcome, not to either")
    print("  generator -- which is exactly why it is the only therapy in ARVC")
    print("  with an unambiguous mortality benefit, and why it cannot be a")
    print("  substitute for the clock-directed interventions.")


def scen_gene_therapy():
    print("\n" + "=" * 78)
    print(" SCENARIO 7 -- AAV9-PKP2 GENE THERAPY: TIMING IS THE WHOLE STORY")
    print("=" * 78)
    print(" A single infusion restoring PKP2 to ~0.5 + 0.53 of normal at the")
    print(" intercalated disc, given at different points on the same")
    print(" trajectory (PKP2tv, competitive athlete, beta-blocked).")
    print(" Gene therapy lowers the price per unit load; it does not remove")
    print(" scar that has already formed.  [PMID 39196150]\n")
    print("  %-22s %9s %9s %10s %10s %9s" %
          ("infusion", "PKP2eff", "fibfat@55", "RVEF@55", "VA/yr@55", "age_dx"))
    common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline"), (16.0, "DRUG", "BB")]
    for label, ev in [("none", []),
                      ("age 18 (concealed)", [(18.0, "AAV", 0.62)]),
                      ("age 26 (early)", [(26.0, "AAV", 0.62)]),
                      ("age 34 (overt)", [(34.0, "AAV", 0.62)]),
                      ("age 45 (late)", [(45.0, "AAV", 0.62)])]:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=common + ev, age_end=60.0)
        o = at_age(traj, 55.0)
        adx = age_at(traj, lambda q: q["TFC"] == "definite")
        print("  %-22s %9.2f %9.4f %10.2f %10.3f %9s" %
              (label, o["PKP2_EFF"], o["FIBFAT_RV"], o["RVEF"],
               o["H_VA_PER_YR"], ("%.1f" % adx) if adx < 1e6 else "never"))


def scen_hot_phase():
    print("\n" + "=" * 78)
    print(" SCENARIO 8 -- THE HOT PHASE AND THE UPSTREAM ARMS")
    print("=" * 78)
    print(" Myocyte death in arrhythmogenic cardiomyopathy is immunogenic;")
    print(" NF-kB-driven inflammation both follows injury and feeds back on")
    print(" it [PMID 31533459, 27170944].  Anti-inflammatory and")
    print(" anti-adipogenic arms therefore act on the amplifier, not the")
    print(" clock -- so they slow the disease without stopping it.\n")
    print("  %-34s %9s %9s %9s %9s" %
          ("arm", "fibfat@50", "fat@50", "RVEF@50", "age_dx"))
    common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline"), (16.0, "DRUG", "BB")]
    for label, ev in [
        ("none", []),
        ("IL-1 blockade from 24", [(24.0, "DRUG", "IL1")]),
        ("corticosteroid from 24", [(24.0, "DRUG", "GC")]),
        ("MRA from 24", [(24.0, "DRUG", "MRA")]),
        ("GSK-3beta inhibition from 24", [(24.0, "DRUG", "GSK")]),
        ("all four from 24", [(24.0, "DRUG", "IL1"), (24.0, "DRUG", "MRA"),
                              (24.0, "DRUG", "GSK")]),
        ("exercise restriction from 24", [(24.0, "EX", "sedentary")]),
    ]:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=common + ev, age_end=60.0)
        o = at_age(traj, 50.0)
        adx = age_at(traj, lambda q: q["TFC"] == "definite")
        print("  %-34s %9.4f %9.3f %9.2f %9s" %
              (label, o["FIBFAT_RV"], o["FAT_RV"], o["RVEF"],
               ("%.1f" % adx) if adx < 1e6 else "never"))


def scen_amiodarone_cost():
    print("\n" + "=" * 78)
    print(" SCENARIO 9 -- THE MOST EFFECTIVE DRUG IS THE WORST LIFETIME BET")
    print("=" * 78)
    print(" Amiodarone's 58-day terminal half-life and enormous peripheral")
    print(" volume are in the model, so cumulative exposure is a state")
    print(" variable.  ARVC patients are typically diagnosed in their")
    print(" thirties, so 'most effective per year' and 'right choice for the")
    print(" next forty years' are different questions.\n")
    print("  %-24s %10s %10s %10s %10s" %
          ("arm", "VA/yr@40", "VA/yr@60", "cumtox@60", "VA-free@60"))
    common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline"), (16.0, "DRUG", "BB"),
              (38.0, "ICD", 1)]
    for label, ev in [("beta-blocker only", []),
                      ("+ amiodarone from 38", [(38.0, "DRUG", "AM")]),
                      ("+ amiodarone 38-45 only", [(38.0, "DRUG", "AM"),
                                                   (45.0, "DRUG", "-AM")]),
                      ("+ endo-epi ablation 38", [(38.0, "ABL", "epi")])]:
        traj, p = simulate(genotype="PKP2tv", exercise="competitive",
                           events=common + ev, age_end=60.0)
        print("  %-24s %10.3f %10.3f %10.2f %10.3f" %
              (label, at_age(traj, 40.0)["H_VA_PER_YR"],
               traj[-1]["H_VA_PER_YR"], traj[-1]["AMIO_TOX"],
               traj[-1]["VA_FREE"]))
    print("\n  STATED MISS: this model gives beta-blocker monotherapy a larger")
    print("  arrhythmia reduction than the North American registry observed")
    print("  (which found no significant beta-blocker benefit on VT")
    print("  endpoints, PMID 19660690).  The model routes almost all of")
    print("  Generator I through beta-adrenergic drive; if the registry is")
    print("  right, part of the RyR2 leak must be adrenergic-INDEPENDENT.")
    print("  That is the cleanest place to try to falsify commitment 3.")
    print("")
    print("  A SECOND MISS, WHICH IS THE SAME MISS: scenario 5 shows")
    print("  flecainide ADDED to a beta-blocker coming out slightly harmful at")
    print("  every age, whereas Ermakov 2017 (PMID 27939893) found it reduced")
    print("  VT.  The reason is that beta-blockade has already driven")
    print("  Generator I to zero here, leaving only flecainide's conduction")
    print("  cost.  Both misses therefore say the same thing -- some of the")
    print("  RyR2 leak has to be adrenergic-INDEPENDENT.  One added term")
    print("  repairs both, and that is a testable prediction about the model")
    print("  rather than an excuse for it.")


def fit_kinj(phi_ex, target_age=39.4, genotype="PKP2tv", exercise="guideline",
             lo=1.0e-7, hi=1.0e-3, iters=22):
    """Refit the fatigue scale so a model variant hits the SAME anchor.

    This is what makes the falsifier a fair test rather than a rigged one.
    Setting PHI_EX = 0 does not merely delete the exercise effect, it also
    makes the whole disease slower, so a naive comparison would just show the
    calendar-clock model failing to produce any disease at all.  The honest
    comparison is: fit the calendar-clock model to the SAME calibration anchor
    the fatigue model was fitted to (median age at definite diagnosis on
    ordinary recreational activity), and then ask whether it still reproduces
    the predictions the fatigue model got for free.
    """
    def age_dx(k):
        tr, _ = simulate(genotype=genotype, exercise=exercise,
                         params=dict(K_INJ=k, PHI_EX=phi_ex),
                         events=[(12.0, "BIOPSY", 1)], age_end=95.0,
                         record_every=180.0,
                         stop_when=lambda o: o["TFC"] == "definite")
        return age_at(tr, lambda o: o["TFC"] == "definite")
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        a = age_dx(mid)
        if a > target_age:          # too slow -> need a larger scale
            lo = mid
        else:
            hi = mid
        if hi / lo < 1.02:
            break
    return math.sqrt(lo * hi)


def scen_falsifier():
    print("\n" + "=" * 78)
    print(" SCENARIO 10 -- THE FALSIFIER: PHI_EX = 0 (LOAD-INDEPENDENT DAMAGE)")
    print("=" * 78)
    print(" One parameter turns the fatigue clock into a calendar clock.")
    print(" Everything else -- every rate constant, every Ki, the Task Force")
    print(" scoring -- is untouched.\n")
    k0 = fit_kinj(0.0)
    print(" Both variants are FITTED TO THE SAME ANCHOR first: the median age")
    print(" at definite diagnosis on ordinary recreational activity.  The")
    print(" calendar-clock model needs K_INJ = %.3g to hit it (the fatigue" % k0)
    print(" model uses %.3g).  Only then is it asked to reproduce the" % P["K_INJ"])
    print(" predictions the fatigue model got for free.\n")
    print("  %-9s %-13s %14s %14s" %
          ("genotype", "exercise", "age_dx fatigue", "age_dx calendar"))
    for gt in ("WT", "PKP2tv"):
        for ex in ("sedentary", "guideline", "competitive"):
            a1 = age_at(simulate(genotype=gt, exercise=ex,
                                 events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                        lambda o: o["TFC"] == "definite")
            a0 = age_at(simulate(genotype=gt, exercise=ex,
                                 params=dict(PHI_EX=0.0, K_INJ=k0),
                                 events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                        lambda o: o["TFC"] == "definite")
            print("  %-9s %-13s %14s %14s" %
                  (gt, ex, ("%.1f" % a1) if a1 < 1e6 else "never",
                   ("%.1f" % a0) if a0 < 1e6 else "never"))
    tr1, _ = simulate(genotype="PKP2tv", exercise="competitive", age_end=50.0)
    tr0, _ = simulate(genotype="PKP2tv", exercise="competitive",
                      params=dict(PHI_EX=0.0, K_INJ=k0), age_end=50.0)
    r1 = tr1[-1]["FIBFAT_RV"] / max(1e-9, tr1[-1]["FIBFAT_LV"])
    r0 = tr0[-1]["FIBFAT_RV"] / max(1e-9, tr0[-1]["FIBFAT_LV"])
    print("\n  RV:LV fibrofatty ratio at 50 y, PHI_EX=1 : %.2f  (RV disease)" % r1)
    print("  RV:LV fibrofatty ratio at 50 y, PHI_EX=0 : %.2f  (no chamber"
          " preference)" % r0)
    lr = simulate(genotype="PKP2tv", exercise="competitive",
                  params=dict(PHI_EX=0.0, K_INJ=k0),
                  events=[(12.0, "BIOPSY", 1), (16.0, "DRUG", "LOADRED")],
                  age_end=60.0)[0]
    ln = simulate(genotype="PKP2tv", exercise="competitive",
                  params=dict(PHI_EX=0.0, K_INJ=k0), events=[(12.0, "BIOPSY", 1)],
                  age_end=60.0)[0]
    print("  Load-reducing therapy with PHI_EX=0: fibfat@50 %.4f vs %.4f"
          " (no effect)" % (at_age(lr, 50.0)["FIBFAT_RV"],
                            at_age(ln, 50.0)["FIBFAT_RV"]))
    print("\n  With the load term removed the model predicts: exercise is")
    print("  irrelevant, penetrance is complete and age-fixed, gene-elusive")
    print("  ARVC is impossible, the ventricles are affected equally, and")
    print("  load-reducing therapy does nothing.  All five are contradicted")
    print("  by data the model was not fitted to.  That is what makes")
    print("  commitment 1 load-bearing rather than decorative.")


# =============================================================================
#  SECTION 9 -- PK CONSISTENCY CHECK (the Css shortcut, demonstrated)
# =============================================================================

def pk_ode_css(dose, F, CL, V, tau=1.0, ka=12.0, n_days=60, dt=1.0 / 2880.0):
    """Fine-step integration of the same 1-cmt oral model the R file uses."""
    Ag, Ac = 0.0, 0.0
    k = CL / V
    t = 0.0
    nxt = 0.0
    auc, tlast = 0.0, 0.0
    while t < n_days:
        if t >= nxt - 1e-12:
            Ag += F * dose
            nxt += tau
        dAg = -ka * Ag
        dAc = ka * Ag - k * Ac
        Ag += dt * dAg
        Ac += dt * dAc
        if t >= n_days - 7.0:                      # average over the last week
            auc += Ac / V * dt
            tlast += dt
        t += dt
    return auc / tlast


# =============================================================================
#  SECTION 10 -- SELF-CHECKS
# =============================================================================

def run_checks():
    checks = []

    def chk(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    # --- 1. the fatigue arithmetic reproduces the observed exercise HR ------
    lr_sed = load_index(EXERCISE["sedentary"], P["K_EX_RV"])
    lr_comp = load_index(EXERCISE["competitive"], P["K_EX_RV"])
    ratio = lr_comp / lr_sed
    chk("exercise_hazard_ratio_predicted", 2.4 <= ratio <= 3.9,
        "RV fatigue rate ratio competitive/sedentary = %.2f (observed HR 3.16, "
        "PMID 23871885)" % ratio)

    # --- 2. RV vs LV asymmetry emerges from Laplace alone ------------------
    rv_lv_comp = (load_index(EXERCISE["competitive"], P["K_EX_RV"])
                  / load_index(EXERCISE["competitive"], P["K_EX_LV"]))
    rv_lv_sed = (load_index(EXERCISE["sedentary"], P["K_EX_RV"])
                 / load_index(EXERCISE["sedentary"], P["K_EX_LV"]))
    chk("rv_lv_gap_widens_with_training", rv_lv_comp > rv_lv_sed > 1.0,
        "RV/LV fatigue ratio %.2f sedentary -> %.2f competitive" %
        (rv_lv_sed, rv_lv_comp))

    # --- 3. calibration target: median age at definite diagnosis -----------
    tr, _ = simulate(genotype="PKP2tv", exercise="guideline",
                     events=[(12.0, "BIOPSY", 1)], age_end=80.0)
    adx = age_at(tr, lambda o: o["TFC"] == "definite")
    chk("age_at_diagnosis_calibrated", 32.0 <= adx <= 46.0,
        "PKP2tv + guideline activity, median frailty: definite TFC at age %.1f "
        "(observed median 31-36 y, PMID 25820315)" % adx)

    # --- 4. incomplete penetrance in sedentary carriers -------------------
    tr, _ = simulate(genotype="PKP2tv", exercise="sedentary",
                     events=[(12.0, "BIOPSY", 1)], age_end=80.0)
    adx_sed = age_at(tr, lambda o: o["TFC"] == "definite")
    chk("incomplete_penetrance_when_sedentary", adx_sed > adx + 8.0,
        "sedentary PKP2tv reaches definite at %s vs %.1f on guideline activity"
        % (("%.1f" % adx_sed) if adx_sed < 1e6 else "never", adx))

    # --- 5. gene-elusive ARVC requires extreme load -----------------------
    a_wt_elite = age_at(simulate(genotype="WT", exercise="elite",
                                 events=[(12.0, "BIOPSY", 1)], age_end=90.0)[0],
                        lambda o: o["TFC"] == "definite")
    a_wt_sed = age_at(simulate(genotype="WT", exercise="sedentary",
                               events=[(12.0, "BIOPSY", 1)], age_end=90.0)[0],
                      lambda o: o["TFC"] == "definite")
    chk("gene_elusive_needs_extreme_load", a_wt_elite < 75.0 and a_wt_sed > 85.0,
        "variant-negative: elite athlete %s, sedentary %s (PMID 25516436)" %
        (("%.1f" % a_wt_elite) if a_wt_elite < 1e6 else "never",
         ("%.1f" % a_wt_sed) if a_wt_sed < 1e6 else "never"))

    # --- 6. RV before LV in PKP2, with no chamber-specific biology ---------
    tr, _ = simulate(genotype="PKP2tv", exercise="competitive", age_end=50.0)
    o = tr[-1]
    chk("rv_dominant_in_pkp2", o["FIBFAT_RV"] > 3.0 * o["FIBFAT_LV"],
        "fibrofatty RV %.4f vs LV %.4f (ratio %.1f) with KAPPA_LV = 1" %
        (o["FIBFAT_RV"], o["FIBFAT_LV"], o["FIBFAT_RV"] / max(1e-9, o["FIBFAT_LV"])))

    # --- 7. left dominance only with the DSP/FLNC term --------------------
    trd, _ = simulate(genotype="DSPtv", exercise="competitive", age_end=50.0)
    od = trd[-1]
    chk("left_dominance_requires_genotype_term",
        (od["FIBFAT_LV"] / max(1e-9, od["FIBFAT_RV"]))
        > 2.5 * (o["FIBFAT_LV"] / max(1e-9, o["FIBFAT_RV"])),
        "LV:RV ratio DSPtv %.2f vs PKP2tv %.2f" %
        (od["FIBFAT_LV"] / max(1e-9, od["FIBFAT_RV"]),
         o["FIBFAT_LV"] / max(1e-9, o["FIBFAT_RV"])))

    # --- 8. arrhythmia hazard scale in definite disease -------------------
    #   read on the SAME trajectory the scale was fitted to: a carrier on
    #   ordinary recreational activity, just after definite diagnosis
    tr, _ = simulate(genotype="PKP2tv", exercise="guideline",
                     events=[(12.0, "BIOPSY", 1)], age_end=44.0)
    h_def = at_age(tr, 42.0)["H_VA_PER_YR"]
    chk("va_hazard_definite_disease", 0.06 <= h_def <= 0.16,
        "sustained-VA hazard just after definite diagnosis = %.1f%%/yr "
        "(observed ~10%%/yr)" % (100 * h_def))

    # --- 9. and it is far lower before the phenotype ---------------------
    #   measured in a genuinely concealed carrier: PKP2tv, sedentary, age 25,
    #   still phenotype-negative
    trc, _ = simulate(genotype="PKP2tv", exercise="sedentary",
                      events=[(12.0, "BIOPSY", 1)], age_end=26.0)
    oc = at_age(trc, 25.0)
    h_conceal = oc["H_VA_PER_YR"]
    chk("va_hazard_concealed_phase",
        h_conceal < 0.015 and oc["TFC"] != "definite",
        "concealed-phase hazard = %.2f%%/yr with TFC '%s' (observed ~0.5-1%%/yr "
        "in phenotype-negative carriers)" % (100 * h_conceal, oc["TFC"]))

    # --- 10. load-reducing therapy prevents the phenotype ----------------
    base = simulate(genotype="PKP2tv", exercise="competitive",
                    events=[(12.0, "BIOPSY", 1)], age_end=60.0)[0]
    lr = simulate(genotype="PKP2tv", exercise="competitive",
                  events=[(12.0, "BIOPSY", 1), (16.0, "DRUG", "LOADRED")],
                  age_end=60.0)[0]
    f_b = at_age(base, 45.0)["FIBFAT_RV"]
    f_l = at_age(lr, 45.0)["FIBFAT_RV"]
    chk("load_reduction_prevents_phenotype", f_l < 0.72 * f_b,
        "fibrofatty at 45 y: %.4f -> %.4f with load reduction alone "
        "(PMID 21292134)" % (f_b, f_l))

    # --- 11. beta-blockade moves arrhythmia but not structure ------------
    #   The claim is NOT "beta-blockade is powerful".  It is that beta-blockade
    #   occupies Generator I and nothing else, so it must (a) leave structure
    #   untouched and (b) lose most of its effect once Generator II dominates.
    #   Testing it as a fixed percentage would be testing the wrong thing.
    base_g = simulate(genotype="PKP2tv", exercise="guideline",
                      events=[(12.0, "BIOPSY", 1)], age_end=50.0)[0]
    bb = simulate(genotype="PKP2tv", exercise="guideline",
                  events=[(12.0, "BIOPSY", 1), (16.0, "DRUG", "BB")],
                  age_end=50.0)[0]
    f_b = at_age(base_g, 45.0)["FIBFAT_RV"]
    f_bb = at_age(bb, 45.0)["FIBFAT_RV"]
    r_early = (at_age(bb, 28.0)["H_VA_PER_YR"]
               / at_age(base_g, 28.0)["H_VA_PER_YR"])
    r_late = (at_age(bb, 48.0)["H_VA_PER_YR"]
              / at_age(base_g, 48.0)["H_VA_PER_YR"])
    chk("beta_blocker_is_generator_I_only",
        abs(f_bb - f_b) / f_b < 0.02 and r_early < 0.55 and r_late > r_early + 0.20,
        "structure %.4f vs %.4f (unchanged to %.1f%%); hazard ratio vs no drug "
        "%.2f at 28 y (Generator I dominant) but only %.2f at 48 y "
        "(Generator II dominant)" %
        (f_b, f_bb, 100 * abs(f_bb - f_b) / f_b, r_early, r_late))

    # --- 12. the sotalol null ------------------------------------------
    ev = [(12.0, "BIOPSY", 1)]
    h_none = at_age(simulate(genotype="PKP2tv", exercise="guideline",
                             events=ev, age_end=44.0, record_every=15.0)[0],
                    42.0)["H_VA_PER_YR"]
    h_sot = at_age(simulate(genotype="PKP2tv", exercise="guideline",
                            events=ev + [(40.0, "DRUG", "SO")], age_end=44.0,
                            record_every=15.0)[0], 42.0)["H_VA_PER_YR"]
    h_amio = at_age(simulate(genotype="PKP2tv", exercise="guideline",
                             events=ev + [(40.0, "DRUG", "AM")], age_end=44.0,
                             record_every=15.0)[0], 42.0)["H_VA_PER_YR"]
    chk("sotalol_null_amiodarone_works",
        0.92 * h_none < h_sot < 1.35 * h_none and h_amio < 0.80 * h_none,
        "vs no drug: sotalol %.0f%%, amiodarone %.0f%% (PMID 19660690 found "
        "amiodarone effective and sotalol no better than nothing)" %
        (100 * h_sot / h_none, 100 * h_amio / h_none))

    # --- 13. flecainide changes sign along the trajectory ---------------
    def flec_delta(age):
        #   measured WITHOUT a beta-blocker, so the RyR2 arm has something left
        #   to do; the add-on case is printed in scenario 5
        common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline")]
        a = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                            events=common, age_end=age + 1.0,
                            record_every=15.0)[0], age + 0.9)["H_VA_PER_YR"]
        b = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                            events=common + [(age, "DRUG", "FL")],
                            age_end=age + 1.0, record_every=15.0)[0],
                   age + 0.9)["H_VA_PER_YR"]
        return b / a - 1.0
    d_early, d_late = flec_delta(24.0), flec_delta(46.0)
    chk("flecainide_sign_flips", d_early < -0.02 and d_late > d_early + 0.05,
        "hazard change: %+.0f%% at 24 y, %+.0f%% at 46 y (PMID 27939893 + the "
        "CAST caution, from geometry)" % (100 * d_early, 100 * d_late))

    # --- 14. endo-only vs endo-epi ablation ----------------------------
    h0 = at_age(simulate(genotype="PKP2tv", exercise="competitive", events=ev,
                         age_end=42.0, record_every=15.0)[0], 39.0)["GEN2"]
    h_en = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                           events=ev + [(38.0, "ABL", "endo")], age_end=42.0,
                           record_every=15.0)[0], 39.0)["GEN2"]
    h_ep = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                           events=ev + [(38.0, "ABL", "epi")], age_end=42.0,
                           record_every=15.0)[0], 39.0)["GEN2"]
    chk("epicardial_ablation_beats_endocardial",
        h_ep < 0.55 * h_en < h0,
        "Generator II: none %.3f, endo %.3f, endo+epi %.3f "
        "(PMID 26546346, 22205683)" % (h0, h_en, h_ep))

    # --- 15. ICD: survival benefit, zero substrate benefit -------------
    #   measured on the guideline trajectory: in the end-stage athlete both
    #   arms are dominated by heart-failure death, so the device has little
    #   left to save -- itself a prediction, and the reason the model puts
    #   exercise restriction upstream of the device rather than beside it.
    noicd = simulate(genotype="PKP2tv", exercise="guideline",
                     events=ev + [(16.0, "DRUG", "BB")], age_end=60.0)[0]
    icd = simulate(genotype="PKP2tv", exercise="guideline",
                   events=ev + [(16.0, "DRUG", "BB"), (40.0, "ICD", 1)],
                   age_end=60.0)[0]
    chk("icd_outcome_not_substrate",
        icd[-1]["SURV"] > noicd[-1]["SURV"] + 0.05
        and abs(icd[-1]["FIBFAT_RV"] - noicd[-1]["FIBFAT_RV"]) < 1e-9,
        "survival at 60 y %.3f -> %.3f, fibrofatty identical to %.0e "
        "(all of the mortality benefit, none of the disease)" %
        (noicd[-1]["SURV"], icd[-1]["SURV"],
         abs(icd[-1]["FIBFAT_RV"] - noicd[-1]["FIBFAT_RV"])))

    # --- 16. gene therapy: timing dominates dose -----------------------
    common = [(12.0, "BIOPSY", 1), (16.0, "EX", "guideline"), (16.0, "DRUG", "BB")]
    f_none = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                             events=common, age_end=60.0)[0], 55.0)["FIBFAT_RV"]
    f_early = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                              events=common + [(18.0, "AAV", 0.62)],
                              age_end=60.0)[0], 55.0)["FIBFAT_RV"]
    f_late = at_age(simulate(genotype="PKP2tv", exercise="competitive",
                             events=common + [(45.0, "AAV", 0.62)],
                             age_end=60.0)[0], 55.0)["FIBFAT_RV"]
    ben_early = 1.0 - f_early / f_none
    ben_late = 1.0 - f_late / f_none
    chk("gene_therapy_timing_dominates",
        ben_early >= 2.0 * ben_late,
        "substrate spared at 55 y: %.0f%% if infused at 18 y vs %.0f%% at 45 y "
        "(none %.3f, early %.3f, late %.3f) -- timing worth >2x, PMID 39196150"
        % (100 * ben_early, 100 * ben_late, f_none, f_early, f_late)),

    # --- 17. the falsifier really does invert the predictions ----------
    #   refit the calendar-clock variant to the SAME anchor before comparing,
    #   so the falsifier is a fair test and not a straw man
    k0 = fit_kinj(0.0)
    a1 = age_at(simulate(genotype="PKP2tv", exercise="sedentary",
                         events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                lambda o: o["TFC"] == "definite")
    a2 = age_at(simulate(genotype="PKP2tv", exercise="competitive",
                         events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                lambda o: o["TFC"] == "definite")
    b1 = age_at(simulate(genotype="PKP2tv", exercise="sedentary",
                         params=dict(PHI_EX=0.0, K_INJ=k0),
                         events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                lambda o: o["TFC"] == "definite")
    b2 = age_at(simulate(genotype="PKP2tv", exercise="competitive",
                         params=dict(PHI_EX=0.0, K_INJ=k0),
                         events=[(12.0, "BIOPSY", 1)], age_end=95.0)[0],
                lambda o: o["TFC"] == "definite")
    chk("falsifier_removes_exercise_effect",
        (a1 - a2) > 8.0 and abs(b1 - b2) < 2.5,
        "age-at-dx gap sedentary-vs-competitive: %.1f y with the fatigue "
        "clock, %.1f y with a calendar clock refitted to the same anchor "
        "(K_INJ %.3g)" % (a1 - a2, abs(b1 - b2), k0))

    trf, _ = simulate(genotype="PKP2tv", exercise="competitive",
                      params=dict(PHI_EX=0.0, K_INJ=k0), age_end=50.0)
    of = trf[-1]
    chk("falsifier_removes_chamber_preference",
        (of["FIBFAT_RV"] / max(1e-9, of["FIBFAT_LV"])) < 1.6,
        "RV:LV fibrofatty ratio without the load term = %.2f (vs %.1f with it)"
        % (of["FIBFAT_RV"] / max(1e-9, of["FIBFAT_LV"]),
           o["FIBFAT_RV"] / max(1e-9, o["FIBFAT_LV"])))

    # --- 17b. PENETRANCE AS A POPULATION QUANTITY (all predictions) -----
    pen_sed, _ = cohort_penetrance("PKP2tv", "sedentary", True, n=9)
    pen_gl, _ = cohort_penetrance("PKP2tv", "guideline", True, n=9)
    pen_gl_f, _ = cohort_penetrance("PKP2tv", "guideline", False, n=9)
    pen_comp, _ = cohort_penetrance("PKP2tv", "competitive", True, n=9)
    pen_wt_sed, _ = cohort_penetrance("WT", "sedentary", True, n=9)
    pen_wt_elite, _ = cohort_penetrance("WT", "elite", True, n=9)

    chk("penetrance_incomplete_by_midlife",
        0.30 <= pen_gl[40.0] <= 0.60,
        "PKP2tv carriers on guideline activity: %.0f%% definite by 40 y "
        "(observed roughly 35-50%% by mid-life, PMID 25820315)"
        % (100 * pen_gl[40.0]))
    chk("penetrance_low_when_sedentary",
        pen_sed[40.0] < pen_gl[40.0] and pen_sed[60.0] < 0.72,
        "sedentary carriers: %.0f%% by 40 y, %.0f%% by 60 y (PMID 31638835)"
        % (100 * pen_sed[40.0], 100 * pen_sed[60.0]))
    chk("penetrance_near_complete_in_athletes",
        pen_comp[40.0] >= 0.88,
        "competitive athletes carrying PKP2tv: %.0f%% definite by 40 y "
        "(PMID 23871885, 25896080)" % (100 * pen_comp[40.0]))
    chk("male_predominance_predicted",
        pen_gl[40.0] > pen_gl_f[40.0],
        "by 40 y: male %.0f%% vs female %.0f%% (PMID 28329361; one parameter, "
        "SEX_K_FEMALE, and the ordering is the output)"
        % (100 * pen_gl[40.0], 100 * pen_gl_f[40.0]))
    chk("gene_elusive_cohort_needs_extreme_load",
        pen_wt_sed[60.0] == 0.0 and pen_wt_elite[60.0] >= 0.5,
        "variant-negative cohort by 60 y: sedentary %.0f%%, elite athlete "
        "%.0f%% (PMID 25516436)"
        % (100 * pen_wt_sed[60.0], 100 * pen_wt_elite[60.0]))

    # --- 18. the Css shortcut equals the ODE it replaces ---------------
    css_alg = P["FL_F"] * P["FL_DOSE"] / P["FL_CL"]
    css_ode = pk_ode_css(P["FL_DOSE"] / 2.0, P["FL_F"], P["FL_CL"],
                         8.0 * 75.0, tau=0.5)
    chk("pk_css_matches_ode", abs(css_ode / css_alg - 1.0) < 0.02,
        "flecainide Css,avg analytic %.4f mg/L vs fine-step ODE %.4f mg/L "
        "(%.1f%%)" % (css_alg, css_ode, 100 * (css_ode / css_alg - 1)))

    # --- 19. mass balance of myocardial composition -------------------
    tr, _ = simulate(genotype="PKP2tv", exercise="competitive", age_end=60.0)
    worst = max(abs(o["MYO_RV"] + o["FIB_RV"] + o["FAT_RV"] - 1.0) for o in tr)
    chk("myocardial_composition_bounded", worst < 0.12,
        "max |MYO+FIB+FAT - 1| over 48 y = %.4f (replacement fills vacated "
        "space, it is not forced to)" % worst)

    # --- 20. integrator convergence ----------------------------------
    c1 = simulate(genotype="PKP2tv", exercise="competitive",
                  params=dict(DT=1.0), age_end=50.0)[0][-1]["FIBFAT_RV"]
    c2 = simulate(genotype="PKP2tv", exercise="competitive",
                  params=dict(DT=0.25), age_end=50.0)[0][-1]["FIBFAT_RV"]
    chk("integrator_converged", abs(c1 / c2 - 1.0) < 0.015,
        "fibrofatty at 50 y: dt=1.0 d %.5f vs dt=0.25 d %.5f (%.2f%%)" %
        (c1, c2, 100 * abs(c1 / c2 - 1)))

    print("\n" + "=" * 78)
    print(" SELF-CHECKS")
    print("=" * 78)
    npass = 0
    for name, ok, detail in checks:
        print(" [%s] %-42s %s" % ("PASS" if ok else "FAIL", name, detail))
        npass += ok
    print("-" * 78)
    print(" %d/%d checks passed" % (npass, len(checks)))
    print(" fitted parameters: %s" % ", ".join(
        "%s=%.4g" % (k, P[k]) for k in FITTED_PARAMS))
    print("=" * 78)
    return npass == len(checks)


# =============================================================================
#  MAIN
# =============================================================================

def main():
    args = sys.argv[1:]
    if "--check" in args:
        sys.exit(0 if run_checks() else 1)
    if "--falsify" in args:
        scen_falsifier()
        return
    print("=" * 78)
    print(" ARVC QSP MODEL -- REFERENCE IMPLEMENTATION")
    print(" %d ODE states | %d fitted parameters | dependency-free"
          % (NSTATE, len(FITTED_PARAMS)))
    print("=" * 78)
    print_load_table()
    scen_natural_history()
    scen_penetrance()
    scen_chamber_selectivity()
    scen_load_reduction()
    scen_antiarrhythmics()
    scen_flecainide_sign()
    scen_icd()
    scen_gene_therapy()
    scen_hot_phase()
    scen_amiodarone_cost()
    scen_falsifier()
    run_checks()


if __name__ == "__main__":
    main()
