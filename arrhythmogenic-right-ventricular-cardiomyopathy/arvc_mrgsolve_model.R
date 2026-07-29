# =============================================================================
#  arvc_mrgsolve_model.R
#  Arrhythmogenic Right Ventricular Cardiomyopathy (ARVC / arrhythmogenic
#  cardiomyopathy) — QSP model for mrgsolve
# =============================================================================
#
#  54 ODE compartments.  TIME IS DAYS SINCE AGE 12 YEARS — the model runs over
#  a lifetime, because its central claim is that the disease has no calendar
#  clock: age enters only as the integral of mechanical load.
#
#  ---------------------------------------------------------------------------
#  WHAT THIS MODEL IS FOR
#  ---------------------------------------------------------------------------
#  ARVC is normally written as a linear cascade — desmosomal variant, junction
#  failure, myocyte death, fibrofatty replacement, arrhythmia.  That cascade
#  cannot survive five observations:
#
#   (A) PKP2 truncating variants are far commoner in the population than ARVC
#       is; most carriers of ARVC-gene loss-of-function variants in an
#       unselected biobank have no phenotype.  A cascade starting at the
#       genotype predicts complete penetrance.
#       [Carruth 2019 PMID 31638835; Groeneweg 2015 PMID 25820315]
#   (B) Exercise dose is graded and large: desmosomal carriers in the highest
#       exercise tertile reached the phenotype far younger, hazard ratio about
#       3.16 for VT/death.  [James 2013 PMID 23871885; Ruwald 2015 PMID 25896080]
#   (C) The same phenotype occurs in athletes with NO desmosomal variant, and
#       those patients had done MORE exercise.  [Sawant 2014 PMID 25516436]
#   (D) The right ventricle goes first, although PKP2 is expressed equally in
#       both ventricles.  Left dominance is a DSP/FLNC phenomenon.
#       [Corrado 2024 PMID 37844667; Brandao 2023 PMID 37048743]
#   (E) The therapy hierarchy is upside down for a single-cascade disease:
#       sotalol was no better than nothing in the North American registry,
#       flecainide (a class IC drug, in a structural cardiomyopathy) reduces VT
#       as an add-on, endocardial-only ablation recurs where endo-epicardial
#       does not, and the ICD — which touches no part of the cascade — is the
#       only therapy with an unambiguous mortality benefit.
#       [Marcus 2009 PMID 19660690; Ermakov 2017 PMID 27939893;
#        Santangeli 2015 PMID 26546346; Berruezo 2012 PMID 22205683]
#
#  ---------------------------------------------------------------------------
#  THREE STRUCTURAL COMMITMENTS
#  ---------------------------------------------------------------------------
#  1. THE CLOCK IS CUMULATIVE MECHANICAL WORK, NOT TIME.
#     A desmosomal variant does not destroy myocardium; it lowers the load at
#     which the intercalated disc fails.  Myocyte loss is a FATIGUE-FAILURE
#     process:  rate = beats x (wall stress / reserve)^NMECH.  NMECH is fixed
#     a priori at 4 (the Basquin exponent range for load-bearing biological
#     tissue) and is NOT fitted.  Genotype enters at exactly one place, the
#     reserve term.
#     Because wall stress depends on wall thickness and wall thickness depends
#     on how much myocardium is left, the process is autocatalytic — which is
#     what produces a decades-long concealed phase followed by rapid overt
#     progression, with no phase-transition parameter anywhere in the file.
#
#     THE FALSIFIER IS ONE PARAMETER.  PHI_EX = 0 makes damage
#     load-independent and turns this model into a calendar-clock model.  Then
#     (A), (B), (C) and (D) invert simultaneously: exercise stops mattering,
#     penetrance becomes complete and age-fixed, gene-elusive ARVC becomes
#     impossible, and the two ventricles are affected equally.
#
#     The commitment has an experimental test it was NOT fitted to: in
#     plakoglobin-deficient mice, LOAD-REDUCING THERAPY (furosemide + nitrate)
#     — no desmosomal, ion-channel or anti-fibrotic action whatever —
#     prevented RV enlargement, preserved conduction and prevented the
#     arrhythmic phenotype.  [Fabritz 2011 PMID 21292134;
#     Kirchhof 2006 PMID 17030684].  It is parameter ON_LOADRED here.
#
#  2. RV SELECTIVITY IS LAPLACE, NOT BIOLOGY.
#     No chamber-specific gene expression is used for PKP2 disease.  Two
#     measured asymmetries do the work: the RV free wall is ~4 mm against the
#     LV's ~9 mm, and at maximal exercise RV end-systolic wall stress rises
#     about 125% while the LV rises about 14%, because pulmonary vascular
#     resistance falls far less than systemic.
#     [La Gerche 2011 PMID 21085033; La Gerche 2012 PMID 22160404]
#     Raised to the fourth power, 125%-vs-14% becomes roughly fifteen-fold per
#     exercising beat.  So the RV goes first, and the RV/LV gap WIDENS with
#     training — both predictions.  Left dominance requires KAPPA_LV, the only
#     chamber-specific term in the model, switched on only for DSP and FLNC.
#
#  3. TWO ARRHYTHMIA GENERATORS, AND A DRUG CANNOT OUTPERFORM THE GENERATOR IT
#     OCCUPIES.
#       GENERATOR I  — early, catecholaminergic, structure-independent.  PKP2
#         loss strips Nav1.5 and Cx43 from the intercalated disc and
#         destabilises RyR2, so diastolic Ca leak and triggered activity exist
#         in a structurally normal heart.  This is why a young carrier can die
#         during exertion with a normal echocardiogram.
#         [Sato 2009 PMID 19661460; Cerrone 2014 PMID 24352520;
#          van Opbergen 2019 PMID 31438494]
#       GENERATOR II — late, re-entrant, scar-dependent.  Maximal when viable
#         myocardium and fibrofatty replacement interdigitate about 50/50: a
#         homogeneous scar has no conducting channels and neither does normal
#         myocardium.  Subepicardial before endocardial.
#         [Tschabrunn 2022 PMID 34883271; Berruezo 2012 PMID 22205683]
#     Drugs are attached to specific NODES with occupancy from published
#     Ki/IC50 and PK, not to "arrhythmia".  Two consequences are the model's
#     most falsifiable output:
#       * FLECAINIDE HAS TWO SIGNS — RyR2 block removes Generator I while
#         Nav1.5 block slows conduction and FEEDS Generator II, so its sign
#         depends on where the patient is on the trajectory.
#       * THE ICD IS ATTACHED TO THE OUTCOME, not to either generator: all of
#         the mortality benefit, exactly none of the disease.
#
#  ---------------------------------------------------------------------------
#  CALIBRATION — WHAT IS FITTED AND WHAT IS NOT
#  ---------------------------------------------------------------------------
#  FITTED (4 parameters out of 130):
#    K_INJ   fatigue-failure scale   -> median age at definite Task Force
#                                       diagnosis in PKP2 carriers on ordinary
#                                       recreational activity [PMID 25820315]
#    H0_VA   arrhythmia hazard scale -> ~10%/yr sustained VA in definite ARVC
#    LAM2    Generator II weight     -> the amiodarone-versus-sotalol gap
#                                       [PMID 19660690]
#    K_DIL   dilatation gain         -> the RVEDVi trajectory into overt disease
#
#  PREDICTED (not fitted — this is where the model is exposed):
#    * the ~3x exercise hazard ratio and its dose-response [PMID 23871885]
#    * incomplete penetrance: ~44% of carriers definite by 40 y on guideline
#      activity, ~22% if sedentary, ~100% if competitive [PMID 31638835]
#    * gene-elusive ARVC in extreme-dose athletes only [PMID 25516436]
#    * male predominance, from one multiplier [PMID 28329361]
#    * RV before LV with no chamber-specific biology, and the RV/LV gap
#      widening with training [PMID 21085033]
#    * load-reducing therapy preventing the phenotype [PMID 21292134]
#    * the sotalol null — its IKr and proarrhythmic arms cancel [PMID 19660690]
#    * flecainide helping early and hurting late: one drug, two signs
#    * the endocardial-only versus endo-epicardial ablation gap
#      [PMID 26546346, 22205683]
#    * the ICD's outcome-only benefit
#    * AAV-PKP2 gene therapy in which TIMING dominates DOSE [PMID 39196150]
#
#  A STATED MISS: this model gives beta-blocker monotherapy a larger arrhythmia
#  reduction than the North American registry observed (which found no
#  significant beta-blocker benefit on VT endpoints, PMID 19660690).  Almost
#  all of Generator I is routed through beta-adrenergic drive here; if the
#  registry is right, part of the RyR2 leak must be adrenergic-INDEPENDENT.
#  That is the cleanest place to try to falsify commitment 3.
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
#  ---------------------------------------------------------------------------
#  VERIFICATION
#  ---------------------------------------------------------------------------
#  Every number quoted above is regenerated by the dependency-free Python twin
#  of this same system (47 states there against 54 here -- see the two stated
#  differences below):
#
#      python3 arvc_reference_impl.py            # all scenarios
#      python3 arvc_reference_impl.py --check    # PASS/FAIL
#      python3 arvc_reference_impl.py --falsify  # PHI_EX = 0
#
#  The twin needs only python3, so the calibration can be checked without R,
#  mrgsolve or a compiler.  It differs from this file in exactly two respects,
#  both stated here rather than buried:
#    (i)  the twin represents the fast oral drugs (nadolol, flecainide,
#         sotalol) by their exact average steady-state concentration
#         Css = F*Dose/(CL*tau) instead of integrating absorption at a 1-day
#         step, and it verifies that shortcut numerically against the same
#         ODEs (check "pk_css_matches_ode").  THIS file carries the full
#         absorption/distribution compartments, so it can additionally
#         simulate loading, missed doses and non-adherence.
#    (ii) the twin deposits an AAV dose straight into the transduced-cell
#         state, whereas this file routes it through a vector compartment with
#         a 9-day transduction delay.  Over a decades-long trajectory the
#         difference is immaterial, but it is a difference.
#
#  ---------------------------------------------------------------------------
#  UNITS
#    time            days since age 12 y      (AGE = 12 + TIME/365.25)
#    drug amounts    mg                       concentrations mg/L and nM
#    volumes         L (drug) / mL (cardiac)
#    tissue          dimensionless fractions of original RV free-wall mass
#    conduction      cm/s                     PVC burden per 24 h
#    NT-proBNP       pg/mL                    exercise dose MET-hours/week
#    hazards         per year (integrated states are cumulative hazards)
#
#  DISCLAIMER: educational / research QSP model.  Semi-quantitative, not
#  independently validated, and NOT for clinical decision-making.  The
#  AAV-PKP2, GSK-3beta-inhibition, IL-1-blockade and load-reduction arms are
#  investigational or extrapolated from animal data.
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

arvc_code <- '
$PROB
# ARVC / arrhythmogenic cardiomyopathy — 54-compartment QSP model
# The clock is cumulative mechanical work, not time.

$PARAM @annotated
// ---------------------------------------------------------------- time base
AGE0      :  12.0  : simulation start age (y)

// ------------------------------- COMMITMENT 1: the fatigue-failure clock ---
NMECH     :   4.0  : fatigue exponent, FIXED A PRIORI (Basquin range) (-)
PHI_EX    :   1.0  : 1 = load-dependent damage; 0 = THE FALSIFIER (-)
HR_REST   :  62.0  : resting heart rate (bpm)
HR_EX     : 158.0  : vigorous-intensity heart rate (bpm)
MET_VIG   :   8.0  : intensity credited as vigorous (MET)
EX        :  15.0  : exercise dose (MET-hours/week)
EX2       :  15.0  : exercise dose AFTER restriction (MET-hours/week)
T_RESTR   : 1.0E9  : day on which the exercise dose changes to EX2 (d)
K_EX_RV   :  1.25  : fractional rise in RV ESWS at max exercise, PMID 21085033
K_EX_LV   :  0.14  : fractional rise in LV ESWS at max exercise, PMID 21085033
GAMMA_G   :   2.0  : steepness of reserve loss into failure rate (-)
K_INJ     : 2.80E-6: FITTED fatigue-failure scale (1/d)
FRAILTY   :   1.0  : between-subject multiplier on K_INJ (-)
SEX_K     :   1.0  : 1.0 male, 0.72 female (PMID 28329361) (-)
K_REG     : 1.0E-6 : adult myocyte regeneration, essentially nil (1/d)
MYO_FLOOR :  0.38  : regionally spared myocardium; injury acts on (MYO-floor)
KAPPA_LV  :   1.0  : LV-intrinsic vulnerability; >1 ONLY for DSP/FLNC (-)
XI_FIB    :  0.70  : load-bearing contribution of fibrous tissue (-)
XI_FAT    :  0.35  : bulk (non-load-bearing) contribution of fat (-)

// -------------------------------------------------- genotype / desmosome ---
PKP2_SET  :  0.50  : PKP2 set point (1.0 WT, 0.50 truncating het) (-)
DSP_SET   :  1.00  : desmoplakin set point (0.50 for DSP truncating) (-)
VARIANT   :   1.0  : pathogenic variant present -> TFC category VI major (-)
TAU_PKP2  :  30.0  : protein turnover at the disc (d)
K_PG      :  0.06  : nuclear plakoglobin appearance (1/d)
K_PG_OFF  :  0.05  : nuclear plakoglobin clearance (1/d)
KM_PG     :  0.55  : plakoglobin half-effect on Wnt (-)
K_W       :  0.09  : Wnt activity turnover in (1/d)
K_W_OFF   :  0.09  : Wnt activity turnover out (1/d)
K_H       :  0.05  : Hippo activation by junction loss (1/d)
K_H2      : 0.010  : Hippo activation by mechanical load (1/d)
K_H_OFF   :  0.06  : Hippo deactivation (1/d)
ETA_NAV   :  0.55  : exponent PKP2 -> Nav1.5 at the disc (-)
ETA_CX    :  0.75  : exponent PKP2 -> Cx43 at the disc (-)
ECV       :  0.30  : exponent (gNa x coupling) -> conduction velocity (-)
TAU_ID    :  14.0  : disc protein re-equilibration (d)

// ------------------------------------- GENERATOR I: Ca / triggered activity
K_RY      :  0.10  : RyR2 leak appearance (1/d)
K_RY_OFF  :  0.10  : RyR2 leak resolution (1/d)
BASE_RY   :  0.10  : baseline RyR2 lability in WT (-)
RY_GAIN   :  0.30  : extra lability per unit PKP2 lost, PMID 31438494 (-)
K_UP      :   1.0  : SR Ca uptake (1/d, scaled)
K_REL     :  0.22  : SR Ca release through the leak (1/d)
K_EXT     :   1.4  : cytosolic Ca extrusion (1/d)
PVC_MAX   : 4000   : maximal PVC burden (per 24 h)
PVC_BASE  :  40.0  : normal adult PVC burden (per 24 h)
K50_TRIG  :  0.13  : half-effect of the trigger index (-)
HILL_TRIG :   4.0  : Hill coefficient, trigger -> ectopy (-)
TAU_PVC   :   5.0  : PVC burden re-equilibration (d)
W_SRCSNK  :   6.0  : loss of source-sink suppression per unit myocyte loss (-)

// ---------------------------------------- GENERATOR II: scar re-entry ------
CV0       :  55.0  : healthy RV free-wall conduction velocity (cm/s)
KM_CVFIB  :  0.55  : fibrofatty tortuosity half-effect (-)
FF_THRESH :  0.05  : minimum excess fibrofatty able to host a re-entry circuit (-)
TAU_CV    :  20.0  : conduction re-equilibration (d)
N_CV      :   1.5  : how strongly slow conduction feeds re-entry (-)
TAU_HET   :  60.0  : substrate heterogeneity re-equilibration (d)
LAM2      :  1.91  : FITTED Generator II weight in the hazard (-)

// ------------------------------ inflammation / fibro-fatty replacement -----
K_NEC     :  0.14  : clearance of the acutely injured pool (1/d)
K_INF_ON  :   2.6  : inflammation driven by injured myocytes (1/d)
K_INF_OFF :  0.05  : inflammation resolution (1/d)
IMM_SENS  :  0.55  : NF-kB priming by desmosome loss, PMID 31533459 (-)
PHI_INF   :  1.30  : inflammation AMPLIFIES injury (the amplifier, not clock)
K_FAP     : 0.020  : fibro-adipogenic progenitor activation (1/d)
K_FAP_OFF : 0.030  : progenitor deactivation (1/d)
W_HIPPO   :  0.55  : Hippo weight on progenitor activation (-)
W_WNT     :  0.85  : Wnt-suppression weight on progenitor activation (-)
K_FIB     : 0.075  : fibrous replacement of vacated space (1/d)
K_FIB_DEG : 3.5E-4 : fibrous regression (1/d)
K_FAT     : 0.115  : fatty replacement of vacated space (1/d)
K_FAT_DEG : 2.5E-4 : fatty regression (1/d)

// ------------------------------------------------- chamber mechanics -------
BSA       :  1.90  : body surface area (m2)
RVEDVI0   :  88.0  : normal indexed RVEDV (mL/m2)
LVEDVI0   :  68.0  : normal indexed LVEDV (mL/m2)
EF0_RV    :  0.55  : normal RVEF (-)
EF0_LV    :  0.62  : normal LVEF (-)
K_DIL     :  1.55  : FITTED dilatation gain per unit function lost (-)
K_ATH     :   2.2  : physiological athlete RV dilatation (-)
TAU_CONT  :  45.0  : contractile re-equilibration (d)
TAU_VOL   :  60.0  : volume re-equilibration (d)
K_PRV     :  0.55  : RV systolic pressure rise as the RV fails (-)
W_INFCONT :  0.30  : reversible stunning during a hot phase (-)
MALE      :   1.0  : sex, for the TFC RVEDVi thresholds (-)

// ------------------------------------------------------ neurohormonal -----
K_SNS     :  1.25  : sympathetic gain on RV dysfunction (-)
K_SNS2    :  0.85  : sympathetic gain on LV dysfunction (-)
TAU_SNS   :  20.0  : sympathetic re-equilibration (d)
KM_B1     :  1.10  : beta1 downregulation half-effect (-)
TAU_B1    :  25.0  : beta1 density re-equilibration (d)
NT0       :  45.0  : baseline NT-proBNP (pg/mL)
K_BNP_RV  : 520.0  : NT-proBNP gain on RV dilatation (pg/mL)
K_BNP_LV  : 900.0  : NT-proBNP gain on LV dysfunction (pg/mL)
TAU_BNP   :  10.0  : NT-proBNP kinetics (d)

// ------------------------------------------------------- event hazards ----
H0_VA     : 0.2010 : FITTED sustained-VA hazard scale (1/y)
HMAX_VA   :  0.35  : saturation ceiling on the annual VA hazard (1/y)
P1        :  1.15  : Generator I exponent in the hazard (-)
P2        :  1.30  : Generator II exponent in the hazard (-)
P_SCD     :  0.30  : P(sudden death | sustained VA, no ICD) (-)
ICD_EFF   :  0.88  : P(episode terminated | ICD in situ) (-)
H_HF0     : 0.0022 : baseline heart-failure death hazard (1/y)
K_HF      : 0.115  : HF-death hazard gain on biventricular failure (1/y)
H_ICD_CX  : 0.055  : ICD lead/pocket complication hazard (1/y)
R_ICD_INA : 0.042  : inappropriate shock rate (1/y)
ON_ICD    :   0.0  : ICD implanted (0/1)

// ------------------------------------------------------------- drug PK ----
// nadolol
BB_KA     :  12.0  : absorption (1/d)
BB_F      :  0.30  : bioavailability (-)
BB_CL     : 200.0  : clearance (L/d)
BB_V      : 140.0  : volume (L)
BB_MW     : 309.4  : molecular weight (g/mol)
BB_FU     :  0.70  : unbound fraction (-)
BB_KD     :   4.0  : beta1 dissociation constant, free (nM)
// flecainide
FL_KA     :  20.0  : absorption (1/d)
FL_F      :  0.90  : bioavailability (-)
FL_CL     : 600.0  : clearance (L/d)
FL_V      : 600.0  : volume (L)
FL_MW     : 414.3  : molecular weight (g/mol)
FL_FU     :  0.60  : unbound fraction (-)
FL_IC50NA : 6000.0 : use-dependent Nav1.5 block IC50, free (nM)
FL_IC50RY : 2500.0 : RyR2 open-state block IC50, free (nM)
// sotalol
SO_KA     :  14.0  : absorption (1/d)
SO_F      :  0.95  : bioavailability (-)
SO_CL     : 150.0  : clearance (L/d)
SO_V      : 100.0  : volume (L)
SO_MW     : 272.4  : molecular weight (g/mol)
SO_KD_B1  : 320.0  : beta1 dissociation constant (nM)
SO_IC50IK : 8000.0 : IKr block IC50 (nM)
// amiodarone
AM_KA     :   4.0  : absorption (1/d)
AM_F      :  0.50  : bioavailability (-)
AM_VC     :  60.0  : central volume (L)
AM_VP     : 4600.0 : peripheral volume (L)
AM_CLD    :  48.0  : intercompartmental clearance (L/d)
AM_CL     : 110.0  : clearance (L/d)
AM_MW     : 645.3  : molecular weight (g/mol)
AM_IC50G1 : 1200.0 : Generator I IC50 (nM)
AM_IC50G2 : 1600.0 : Generator II IC50 (nM)
AM_KTOX   : 3.3E-4 : cumulative extracardiac toxicity accrual (1/d)

// ---------------------------------------- upstream / structural arms ------
ON_LOADRED:   0.0  : load-reducing therapy on (0/1), PMID 21292134
E_LOADRED :  0.26  : fractional reduction in the RV volume set point (-)
ON_MRA    :   0.0  : mineralocorticoid antagonist on (0/1)
E_MRA_FIB :  0.30  : fractional reduction in fibrogenesis (-)
ON_IL1    :   0.0  : IL-1 blockade on (0/1)
E_IL1     :  0.70  : fractional increase in inflammation resolution (-)
ON_GC     :   0.0  : corticosteroid on (0/1)
E_GC      :  0.45  : fractional increase in inflammation resolution (-)
ON_GSK    :   0.0  : GSK-3beta inhibition on (0/1), PMID 27170944
E_GSK     :  0.60  : fractional restoration of Wnt / suppression of Hippo (-)

// ------------------------------------------ AAV9-PKP2 gene therapy --------
TAU_TD    :   9.0  : vector uncoating / transduction (d)
TAU_TG    :  21.0  : transgene expression rise (d)
TG_GAIN   :  0.85  : PKP2 restored per unit transduced fraction (-)
TG_HALF   : 3650.0 : episomal transgene loss half-life (d)

// ----------------------------------------------------------- ablation -----
TAU_ABL   : 2600.0 : substrate re-accumulation after ablation (d)

// ----------------------------------------------- reference composition ----
MYO0      :  0.92  : normal myocyte fraction of the RV free wall (-)
FIB0      :  0.03  : normal interstitial fibrous fraction (-)
FAT0      :  0.05  : normal interstitial/epicardial fat fraction (-)

$CMT @annotated
// ---- drug PK (13) -------------------------------------------------------
A_BB_G    : nadolol, gut (mg)
A_BB_C    : nadolol, central (mg)
A_FL_G    : flecainide, gut (mg)
A_FL_C    : flecainide, central (mg)
A_SO_G    : sotalol, gut (mg)
A_SO_C    : sotalol, central (mg)
A_AM_G    : amiodarone, gut (mg)
A_AM_C    : amiodarone, central (mg)
A_AM_P    : amiodarone, peripheral (mg)
AMIO_TOX  : cumulative amiodarone extracardiac toxicity index (-)
A_VEC     : AAV9-PKP2 vector, dosing compartment (vg-equivalents, scaled)
VEC_TD    : transduced cardiomyocyte fraction (-)
PKP2_TG   : transgene-derived PKP2 at the disc (-)
// ---- desmosome and disc (6) --------------------------------------------
PKP2_ID   : PKP2 protein at the intercalated disc (fraction of WT)
DSP_ID    : desmoplakin at the intercalated disc (fraction of WT)
PG_NUC    : nuclear plakoglobin (-)
WNT_ACT   : canonical Wnt/beta-catenin activity (-)
HIPPO     : Hippo pathway activity (-)
NAV_ID    : Nav1.5 at the disc (fraction of WT)
// ---- Generator I (4) ---------------------------------------------------
CX43_ID   : Cx43 at the disc (fraction of WT)
RYR_LEAK  : diastolic SR Ca leak (-)
CA_SR     : SR calcium load (-)
CA_DIA    : diastolic cytosolic calcium (-)
// ---- RV tissue (7) ----------------------------------------------------
D_RV      : cumulative RV mechanical dose (load-years)
MYO_RV    : viable RV myocyte fraction (-)
NEC_RV    : acutely injured RV myocyte pool (-)
INF_RV    : RV inflammatory state (-)
FAP_RV    : RV fibro-adipogenic progenitor activation (-)
FIB_RV    : RV fibrous fraction (-)
FAT_RV    : RV fatty fraction (-)
// ---- LV tissue (7) ----------------------------------------------------
D_LV      : cumulative LV mechanical dose (load-years)
MYO_LV    : viable LV myocyte fraction (-)
NEC_LV    : acutely injured LV myocyte pool (-)
INF_LV    : LV inflammatory state (-)
FAP_LV    : LV fibro-adipogenic progenitor activation (-)
FIB_LV    : LV fibrous fraction (-)
FAT_LV    : LV fatty fraction (-)
// ---- chamber mechanics (4) -------------------------------------------
RVEDV     : RV end-diastolic volume (mL)
RV_CONT   : RV contractile function (-)
LVEDV     : LV end-diastolic volume (mL)
LV_CONT   : LV contractile function (-)
// ---- electrophysiology (5) -------------------------------------------
CV_RV     : RV free-wall conduction velocity (cm/s)
CV_LV     : LV conduction velocity (cm/s)
SCARHETRV : RV scar heterogeneity / conducting-channel index (-)
SCARHETLV : LV scar heterogeneity (-)
PVC24     : PVC burden (per 24 h)
// ---- neurohormonal and events (8) ------------------------------------
SNS       : sympathetic tone (-)
BETA1_D   : beta1-adrenoceptor density (-)
NTBNP     : NT-proBNP (pg/mL)
ABL_HOM   : ablation-induced substrate homogenisation (-)
H_VT      : cumulative sustained-VA hazard (-)
H_DEATH   : cumulative all-cause death hazard (-)
H_HF      : cumulative heart-failure death hazard (-)
ICD_SHK   : cumulative ICD shocks, appropriate + inappropriate (-)

$MAIN
// initial conditions: a healthy heart carrying whatever genotype it carries
PKP2_ID_0  = PKP2_SET;
DSP_ID_0   = DSP_SET;
WNT_ACT_0  = 1.0;
NAV_ID_0   = pow(PKP2_SET, ETA_NAV);
CX43_ID_0  = pow(PKP2_SET, ETA_CX);
RYR_LEAK_0 = RY_GAIN * (1.0 - PKP2_SET) + BASE_RY;
CA_SR_0    = 1.0;
CA_DIA_0   = K_REL * (RY_GAIN * (1.0 - PKP2_SET) + BASE_RY) / K_EXT;
MYO_RV_0   = MYO0;
FIB_RV_0   = FIB0;
FAT_RV_0   = FAT0;
MYO_LV_0   = MYO0;
FIB_LV_0   = FIB0;
FAT_LV_0   = FAT0;
RVEDV_0    = RVEDVI0 * BSA;
LVEDV_0    = LVEDVI0 * BSA;
RV_CONT_0  = 1.0;
LV_CONT_0  = 1.0;
CV_RV_0    = CV0 * pow(pow(PKP2_SET, ETA_NAV) * pow(PKP2_SET, ETA_CX), ECV);
CV_LV_0    = CV0 * 1.05;
PVC24_0    = PVC_BASE;
SNS_0      = 1.0;
BETA1_D_0  = 1.0;
NTBNP_0    = NT0;

$ODE
// =========================================================================
//  0. helpers
// =========================================================================
double H_EFF0 = MYO0 + XI_FIB * FIB0 + XI_FAT * FAT0;

double myo_rv = fmax(0.0, MYO_RV);
double myo_lv = fmax(0.0, MYO_LV);
double fib_rv = fmax(0.0, FIB_RV);
double fib_lv = fmax(0.0, FIB_LV);
double fat_rv = fmax(0.0, FAT_RV);
double fat_lv = fmax(0.0, FAT_LV);
double inf_rv = fmax(0.0, INF_RV);
double inf_lv = fmax(0.0, INF_LV);

// =========================================================================
//  1. DRUG PK AND RECEPTOR OCCUPANCY
//     Each drug is attached to the NODE it occupies, with occupancy computed
//     from published Ki/IC50, unbound fraction and PK — not to a generic
//     antiarrhythmic effect.  This is COMMITMENT 3 made numerical.
// =========================================================================
double c_bb  = BB_FU * (A_BB_C / BB_V) / BB_MW * 1.0E6;          // nM, free
double occ_bb = c_bb / (c_bb + BB_KD);

double c_so  = (A_SO_C / SO_V) / SO_MW * 1.0E6;                  // nM
double occ_so_b1  = c_so / (c_so + SO_KD_B1);
double occ_so_ikr = c_so / (c_so + SO_IC50IK);
double OCC_B1 = 1.0 - (1.0 - occ_bb) * (1.0 - occ_so_b1);

double c_fl  = FL_FU * (A_FL_C / FL_V) / FL_MW * 1.0E6;          // nM, free
double E_FL_RYR   = c_fl / (c_fl + FL_IC50RY);                   // helps Gen I
double E_NA_BLOCK = c_fl / (c_fl + FL_IC50NA);                   // feeds Gen II

double c_am  = (A_AM_C / AM_VC) / AM_MW * 1.0E6;                 // nM
double E_AM_G1 = c_am / (c_am + AM_IC50G1);
double E_AM_G2 = c_am / (c_am + AM_IC50G2);

// repolarisation lengthening: widens the re-entrant wavelength (anti-Gen-II)
// but promotes afterdepolarisations in diseased tissue (pro-Gen-I).  These two
// arms are why sotalol comes out at approximately zero.
double E_WL         = 1.0 - (1.0 - 0.62 * occ_so_ikr) * (1.0 - 0.55 * E_AM_G2);
// The proarrhythmic arm is small relative to the wavelength benefit, and much
// smaller for amiodarone than for sotalol: amiodarone prolongs QT markedly yet
// carries a low torsades rate.  That asymmetry is data, not a modelling choice.
double PROARR_REPOL = 0.18 * occ_so_ikr + 0.08 * E_AM_G2;
// The more important arm: IKr block in a HETEROGENEOUS substrate increases
// dispersion of repolarisation and so promotes re-entry INITIATION.  Without
// this, a lengthened wavelength makes sotalol look good in exactly the
// patients in whom it demonstrably is not.  This term is what turns the drug's
// two arms into the observed null.
// Magnitude set by a harder piece of evidence than the ARVC registry:
// d-sotalol INCREASED mortality in patients with LV dysfunction after
// myocardial infarction (SWORD, PMID 8691967).  A model in which IKr block is
// net-beneficial in a scarred ventricle is contradicted by a randomised trial.
double DISP_REPOL = 0.80 * occ_so_ikr + 0.10 * E_AM_G2;

double EL_LOADRED = ON_LOADRED * E_LOADRED;
double EL_MRA     = ON_MRA * E_MRA_FIB;
double EL_ANTIINF = 1.0 - (1.0 - ON_IL1 * E_IL1) * (1.0 - ON_GC * E_GC);
double EL_GSK     = ON_GSK * E_GSK;

// =========================================================================
//  2. COMMITMENT 2 — THE LAPLACE ENGINE
//     No chamber-specific biology.  The RV goes first because it is thin and
//     because its afterload rises disproportionately on exercise.
// =========================================================================
double h_rv = fmax(0.22, (myo_rv + XI_FIB * fib_rv + XI_FAT * fat_rv) / H_EFF0);
double h_lv = fmax(0.22, (myo_lv + XI_FIB * fib_lv + XI_FAT * fat_lv) / H_EFF0);

double r_rv = pow(fmax(1.0E-6, RVEDV / (RVEDVI0 * BSA)), 1.0 / 3.0);
double r_lv = pow(fmax(1.0E-6, LVEDV / (LVEDVI0 * BSA)), 1.0 / 3.0);

// RV systolic pressure rises as the failing, dilating RV meets its afterload
double p_rv = 1.0 + K_PRV * fmax(0.0, r_rv * r_rv * r_rv - 1.0);

double SIG_RV = p_rv * r_rv / h_rv;      // Laplace, normalised to healthy = 1
double SIG_LV = 1.0   * r_lv / h_lv;

// duty cycle of vigorous exercise over the week.  The dose is allowed to
// change once (T_RESTR), because "the same patient with a different load
// history" is the comparison the model exists to make.
double EXnow = (SOLVERTIME >= T_RESTR) ? EX2 : EX;
double f_ex = fmin(0.35, EXnow / (MET_VIG * 168.0));

// mean-field fatigue dose rate: beats x stress^NMECH, normalised so that a
// perfectly rested healthy chamber scores 1.0.  PHI_EX = 0 removes the stress
// dependence entirely and this becomes a calendar clock (THE FALSIFIER).
double s_rest_rv = (PHI_EX > 0.5) ? pow(SIG_RV, NMECH) : 1.0;
double s_ex_rv   = (PHI_EX > 0.5) ? pow(SIG_RV * (1.0 + K_EX_RV), NMECH) : 1.0;
double s_rest_lv = (PHI_EX > 0.5) ? pow(SIG_LV, NMECH) : 1.0;
double s_ex_lv   = (PHI_EX > 0.5) ? pow(SIG_LV * (1.0 + K_EX_LV), NMECH) : 1.0;

double LOAD_RV = ((1.0 - f_ex) * HR_REST * s_rest_rv
                  + f_ex * HR_EX * s_ex_rv) / HR_REST;
double LOAD_LV = ((1.0 - f_ex) * HR_REST * s_rest_lv
                  + f_ex * HR_EX * s_ex_lv) / HR_REST;

// =========================================================================
//  3. DESMOSOME — the ONLY route from genotype into the model
// =========================================================================
double pkp2_eff = fmin(1.0, fmax(0.0, PKP2_ID) + fmax(0.0, PKP2_TG));
double dsp_eff  = fmax(0.0, DSP_ID);
double desmo    = fmax(0.12, pow(pkp2_eff, 0.7) * pow(dsp_eff, 0.3));
double VULN     = pow(1.0 / desmo, GAMMA_G);

dxdt_PKP2_ID = (PKP2_SET - PKP2_ID) / TAU_PKP2;
dxdt_DSP_ID  = (DSP_SET  - DSP_ID)  / TAU_PKP2;

// nuclear plakoglobin rises as the junction fails (the diagnostic signal loss)
dxdt_PG_NUC = K_PG * (1.0 - desmo) - K_PG_OFF * PG_NUC;

// nuclear plakoglobin competes beta-catenin off TCF/LEF; GSK-3beta inhibition
// restores it.  Wnt is the BRAKE on adipogenesis.
double wnt_t0 = 1.0 / (1.0 + pow(PG_NUC / KM_PG, 2.0));
double wnt_t  = wnt_t0 + EL_GSK * (1.0 - wnt_t0);
dxdt_WNT_ACT = K_W * wnt_t - K_W_OFF * WNT_ACT;

// Hippo is activated by junction loss AND by mechanical load
dxdt_HIPPO = K_H * (1.0 - desmo) + K_H2 * fmax(0.0, LOAD_RV - 1.0)
             - K_H_OFF * HIPPO * (1.0 + EL_GSK);

// PKP2 scaffolds Nav1.5 and the Cx43 plaque: electrical remodelling precedes
// any structural change, which is Generator I existing in a normal heart
dxdt_NAV_ID  = (pow(pkp2_eff, ETA_NAV) - NAV_ID)  / TAU_ID;
dxdt_CX43_ID = (pow(pkp2_eff, ETA_CX)  - CX43_ID) / TAU_ID;

// =========================================================================
//  4. GENERATOR I — Ca leak and triggered activity
// =========================================================================
double beta_act = SNS * (1.0 - OCC_B1) * fmax(0.0, BETA1_D);
double ry_drive = (RY_GAIN * (1.0 - pkp2_eff) + BASE_RY) * beta_act
                  * (1.0 - E_FL_RYR) * (1.0 - 0.45 * E_AM_G1);
dxdt_RYR_LEAK = K_RY * ry_drive - K_RY_OFF * RYR_LEAK;

dxdt_CA_SR  = K_UP * (1.0 - CA_SR) - K_REL * fmax(0.0, RYR_LEAK) * CA_SR;
dxdt_CA_DIA = K_REL * fmax(0.0, RYR_LEAK) * fmax(0.0, CA_SR) - K_EXT * CA_DIA;

// myocyte loss removes the electrotonic (source-sink) suppression of ectopy
double myo_loss_rv = fmax(0.0, fmin(MYO0, MYO0 - myo_rv)) / MYO0;
double trig = fmax(0.0, CA_DIA) * (1.0 + W_SRCSNK * myo_loss_rv)
              * (1.0 + PROARR_REPOL);
double pvc_t = PVC_BASE + PVC_MAX * pow(trig, HILL_TRIG)
               / (pow(K50_TRIG, HILL_TRIG) + pow(trig, HILL_TRIG));
dxdt_PVC24 = (pvc_t - PVC24) / TAU_PVC;

// =========================================================================
//  5. COMMITMENT 1 — THE FATIGUE CLOCK IN THE TISSUE
// =========================================================================
dxdt_D_RV = LOAD_RV * VULN / 365.25;
dxdt_D_LV = LOAD_LV * VULN * KAPPA_LV / 365.25;

double k_inj = K_INJ * FRAILTY * SEX_K;

// ---- right ventricle ----------------------------------------------------
// injury acts on the SUSCEPTIBLE myocardium only.  This is a lumped free-wall
// model but the disease is regional: subtricuspid wall, apex and infundibulum
// carry the highest local stress and go first, while septum-adjacent and
// outflow myocardium stays subcritical, so a lumped average cannot reach zero.
// It also caps the autocatalytic runaway with no phase-transition parameter.
double sus_rv = fmax(0.0, myo_rv - MYO_FLOOR);
double inj_rv = k_inj * LOAD_RV * VULN * 1.0 * sus_rv * (1.0 + PHI_INF * inf_rv);
double space_rv = fmax(0.0, (MYO0 - myo_rv) - (fib_rv - FIB0) - (fat_rv - FAT0));
double adipo = fmax(0.0, 1.0 - WNT_ACT) * fmax(0.0, HIPPO);

dxdt_MYO_RV = -inj_rv + K_REG * space_rv;
dxdt_NEC_RV = inj_rv - K_NEC * fmax(0.0, NEC_RV);
dxdt_INF_RV = K_INF_ON * fmax(0.0, NEC_RV) * (1.0 + IMM_SENS * (1.0 - desmo))
              - K_INF_OFF * inf_rv * (1.0 + 2.2 * EL_ANTIINF);
double fapdrv_rv = inf_rv + W_HIPPO * fmax(0.0, HIPPO)
                   + W_WNT * fmax(0.0, 1.0 - WNT_ACT);
dxdt_FAP_RV = K_FAP * fapdrv_rv * (1.0 - fmax(0.0, FAP_RV))
              - K_FAP_OFF * fmax(0.0, FAP_RV);
dxdt_FIB_RV = K_FIB * fmax(0.0, FAP_RV) * space_rv * (1.0 - EL_MRA)
              - K_FIB_DEG * (fib_rv - FIB0);
dxdt_FAT_RV = K_FAT * fmax(0.0, FAP_RV) * space_rv * adipo
              - K_FAT_DEG * (fat_rv - FAT0);

// ---- left ventricle: identical equations, KAPPA_LV = 1 unless DSP/FLNC ---
double sus_lv = fmax(0.0, myo_lv - MYO_FLOOR);
double inj_lv = k_inj * LOAD_LV * VULN * KAPPA_LV * sus_lv
                * (1.0 + PHI_INF * inf_lv);
double space_lv = fmax(0.0, (MYO0 - myo_lv) - (fib_lv - FIB0) - (fat_lv - FAT0));

dxdt_MYO_LV = -inj_lv + K_REG * space_lv;
dxdt_NEC_LV = inj_lv - K_NEC * fmax(0.0, NEC_LV);
dxdt_INF_LV = K_INF_ON * fmax(0.0, NEC_LV) * (1.0 + IMM_SENS * (1.0 - desmo))
              - K_INF_OFF * inf_lv * (1.0 + 2.2 * EL_ANTIINF);
double fapdrv_lv = inf_lv + W_HIPPO * fmax(0.0, HIPPO)
                   + W_WNT * fmax(0.0, 1.0 - WNT_ACT);
dxdt_FAP_LV = K_FAP * fapdrv_lv * (1.0 - fmax(0.0, FAP_LV))
              - K_FAP_OFF * fmax(0.0, FAP_LV);
dxdt_FIB_LV = K_FIB * fmax(0.0, FAP_LV) * space_lv * (1.0 - EL_MRA)
              - K_FIB_DEG * (fib_lv - FIB0);
dxdt_FAT_LV = K_FAT * fmax(0.0, FAP_LV) * space_lv * adipo
              - K_FAT_DEG * (fat_lv - FAT0);

// =========================================================================
//  6. CHAMBER MECHANICS — where load-reducing therapy acts on the CLOCK
// =========================================================================
double rvct = pow(myo_rv / MYO0, 1.30) * (1.0 - W_INFCONT * fmin(1.0, inf_rv));
double lvct = pow(myo_lv / MYO0, 1.30) * (1.0 - W_INFCONT * fmin(1.0, inf_lv));
dxdt_RV_CONT = (rvct - RV_CONT) / TAU_CONT;
dxdt_LV_CONT = (lvct - LV_CONT) / TAU_CONT;

double rvedv_t = RVEDVI0 * BSA
                 * (1.0 + K_DIL * fmax(0.0, 1.0 - RV_CONT))
                 * (1.0 + K_ATH * f_ex)
                 * (1.0 - EL_LOADRED);
double lvedv_t = LVEDVI0 * BSA
                 * (1.0 + 0.75 * K_DIL * fmax(0.0, 1.0 - LV_CONT))
                 * (1.0 + 0.55 * K_ATH * f_ex)
                 * (1.0 - 0.6 * EL_LOADRED);
dxdt_RVEDV = (rvedv_t - RVEDV) / TAU_VOL;
dxdt_LVEDV = (lvedv_t - LVEDV) / TAU_VOL;

// =========================================================================
//  7. GENERATOR II — conduction and scar re-entry
// =========================================================================
double ff_rv = (fib_rv - FIB0) + (fat_rv - FAT0);
double ff_lv = (fib_lv - FIB0) + (fat_lv - FAT0);

double cv_rv_t = CV0 * pow(fmax(1.0E-4, NAV_ID) * fmax(1.0E-4, CX43_ID), ECV)
                 / (1.0 + fmax(0.0, ff_rv) / KM_CVFIB)
                 * (1.0 - 0.55 * E_NA_BLOCK);
double cv_lv_t = CV0 * 1.05 * pow(fmax(1.0E-4, NAV_ID), ECV)
                 / (1.0 + fmax(0.0, ff_lv) / KM_CVFIB)
                 * (1.0 - 0.55 * E_NA_BLOCK);
dxdt_CV_RV = (cv_rv_t - CV_RV) / TAU_CV;
dxdt_CV_LV = (cv_lv_t - CV_LV) / TAU_CV;

// Scar heterogeneity peaks at ~50/50 interdigitation: normal myocardium has no
// conducting channels and neither does a homogeneous scar.  That is exactly
// the principle substrate ablation exploits.
double het_rv = 4.0 * fmax(0.0, ff_rv) * (myo_rv / MYO0);
double het_lv = 4.0 * fmax(0.0, ff_lv) * (myo_lv / MYO0);
dxdt_SCARHETRV = (het_rv * (1.0 - fmin(1.0, fmax(0.0, ABL_HOM))) - SCARHETRV)
                 / TAU_HET;
dxdt_SCARHETLV = (het_lv - SCARHETLV) / TAU_HET;
dxdt_ABL_HOM   = -ABL_HOM / TAU_ABL;

// =========================================================================
//  8. NEUROHORMONAL
// =========================================================================
double sns_t = 1.0 + K_SNS * fmax(0.0, 1.0 - RV_CONT)
               + K_SNS2 * fmax(0.0, 1.0 - LV_CONT);
dxdt_SNS = (sns_t - SNS) / TAU_SNS;
double b1_t = 1.0 / (1.0 + fmax(0.0, SNS - 1.0) / KM_B1);
dxdt_BETA1_D = (b1_t - BETA1_D) / TAU_B1;

double nt_t = NT0 * (1.0 + 0.6 * (SNS - 1.0))
              + K_BNP_RV * fmax(0.0, RVEDV / (RVEDVI0 * BSA) - 1.0)
              + K_BNP_LV * fmax(0.0, 1.0 - LV_CONT);
dxdt_NTBNP = (nt_t - NTBNP) / TAU_BNP;

// =========================================================================
//  9. DRUG PK COMPARTMENTS
// =========================================================================
dxdt_A_BB_G = -BB_KA * A_BB_G;
dxdt_A_BB_C =  BB_KA * A_BB_G - BB_CL / BB_V * A_BB_C;
dxdt_A_FL_G = -FL_KA * A_FL_G;
dxdt_A_FL_C =  FL_KA * A_FL_G - FL_CL / FL_V * A_FL_C;
dxdt_A_SO_G = -SO_KA * A_SO_G;
dxdt_A_SO_C =  SO_KA * A_SO_G - SO_CL / SO_V * A_SO_C;
dxdt_A_AM_G = -AM_KA * A_AM_G;
dxdt_A_AM_C =  AM_KA * A_AM_G - AM_CL / AM_VC * A_AM_C
               - AM_CLD * (A_AM_C / AM_VC - A_AM_P / AM_VP);
dxdt_A_AM_P =  AM_CLD * (A_AM_C / AM_VC - A_AM_P / AM_VP);
dxdt_AMIO_TOX = AM_KTOX * E_AM_G1;

// AAV9-PKP2: vector clears fast, expression rises over weeks, and the
// transgene restores RESERVE — i.e. it lowers the price per beat.  It removes
// no scar, which is why timing dominates dose.
dxdt_A_VEC   = -A_VEC / TAU_TD;
dxdt_VEC_TD  =  A_VEC / TAU_TD - 0.0 * VEC_TD;
dxdt_PKP2_TG = (TG_GAIN * fmax(0.0, VEC_TD) - PKP2_TG) / TAU_TG
               - log(2.0) / TG_HALF * PKP2_TG;

// =========================================================================
// 10. THE TWO GENERATORS AND THE HAZARDS
// =========================================================================
// Generator I is EXCESS ectopy over the normal background, so a normal heart
// scores zero.  Generator II requires a substrate big enough to host a circuit
// (FF_THRESH), so ordinary age-related interstitial fibrosis is not counted.
double GEN1 = fmax(0.0, PVC24 - PVC_BASE) / 1000.0;
double cvr  = fmax(6.0, CV_RV);
double hetx = fmax(0.0, SCARHETRV - 4.0 * FF_THRESH)
              + 0.45 * fmax(0.0, SCARHETLV - 4.0 * FF_THRESH);
double GEN2 = fmax(0.0, hetx * pow(CV0 / cvr, N_CV)
                        * (1.0 - 0.45 * E_WL) * (1.0 + DISP_REPOL)
                        * (1.0 - 0.30 * E_AM_G2));

// The hazard saturates: a patient cannot have an unbounded number of
// sustained episodes per year before an intervention or a death intervenes.
double h_raw = H0_VA * (pow(GEN1, P1) + LAM2 * pow(GEN2, P2));
double h_va  = HMAX_VA * h_raw / (h_raw + HMAX_VA);
dxdt_H_VT = h_va / 365.25;

double biv  = fmax(0.0, 1.0 - RV_CONT) + 1.4 * fmax(0.0, 1.0 - LV_CONT);
double h_hf = H_HF0 + K_HF * biv * biv;

// The ICD is attached to the OUTCOME, not to either generator.
double h_death = (1.0 - ON_ICD * ICD_EFF) * P_SCD * h_va + h_hf
                 + ON_ICD * H_ICD_CX * 0.10;
dxdt_H_HF    = h_hf / 365.25;
dxdt_H_DEATH = h_death / 365.25;
dxdt_ICD_SHK = ON_ICD * (h_va + R_ICD_INA) / 365.25;

$TABLE
double AGE      = AGE0 + TIME / 365.25;
double RVEF     = EF0_RV * RV_CONT;
double LVEF     = EF0_LV * LV_CONT;
double RVEDVI   = RVEDV / BSA;
double LVEDVI   = LVEDV / BSA;
double FIBFATRV = (FIB_RV - FIB0) + (FAT_RV - FAT0);
double FIBFATLV = (FIB_LV - FIB0) + (FAT_LV - FAT0);
double TOTRV    = fmax(1.0E-9, fmax(0.0, MYO_RV) + FIB_RV + FAT_RV);
double RESIDMYO = 100.0 * fmax(0.0, MYO_RV) / TOTRV;

// ---- ECG phenotype: OUTPUTS of the state variables, never inputs ---------
double TAD_MS   = 32.0 * CV0 / fmax(6.0, CV_RV);
double NTWI     = floor(3.0 * fmin(1.0, fmax(0.0, (MYO0 - MYO_RV) / 0.30)));
double EPSILONW = (CV_RV < 22.0) ? 1.0 : 0.0;
double LATEPOT  = (TAD_MS >= 50.0) ? 1.0 : 0.0;

// ---- the two generators and the instantaneous hazard --------------------
double G1OUT    = fmax(0.0, PVC24 - PVC_BASE) / 1000.0;
double CVR2     = fmax(6.0, CV_RV);
double C_SO2    = (A_SO_C / SO_V) / SO_MW * 1.0E6;
double OSOIKR   = C_SO2 / (C_SO2 + SO_IC50IK);
double C_AM2    = (A_AM_C / AM_VC) / AM_MW * 1.0E6;
double EAMG2    = C_AM2 / (C_AM2 + AM_IC50G2);
double EWL2     = 1.0 - (1.0 - 0.62 * OSOIKR) * (1.0 - 0.55 * EAMG2);
double HETX2    = fmax(0.0, SCARHETRV - 4.0 * FF_THRESH)
                  + 0.45 * fmax(0.0, SCARHETLV - 4.0 * FF_THRESH);
double DISP2    = 0.80 * OSOIKR + 0.10 * EAMG2;
double G2OUT    = fmax(0.0, HETX2 * pow(CV0 / CVR2, N_CV)
                            * (1.0 - 0.45 * EWL2) * (1.0 + DISP2)
                            * (1.0 - 0.30 * EAMG2));
double HRAWOUT  = H0_VA * (pow(G1OUT, P1) + LAM2 * pow(G2OUT, P2));
double HVA_YR   = HMAX_VA * HRAWOUT / (HRAWOUT + HMAX_VA);
double VAFREE   = exp(-H_VT);
double SURV     = exp(-H_DEATH);

// ---- 2010 Task Force Criteria (Marcus 2010 PMID 20172912) ---------------
//  Implemented, not summarised: category counts then the definite/borderline/
//  possible rule.  Category VI is a FREE major for any variant carrier, which
//  is why genotyped relatives are diagnosed earlier than probands.
double regional  = (MYO_RV < 0.85 * MYO0) ? 1.0 : 0.0;
double edvi_maj  = (MALE > 0.5) ? 110.0 : 100.0;
double edvi_min  = (MALE > 0.5) ? 100.0 :  90.0;
double NMAJ = 0.0;
double NMIN = 0.0;
if (regional > 0.5 && (RVEF <= 0.40 || RVEDVI >= edvi_maj)) {
  NMAJ = NMAJ + 1.0;
} else if (regional > 0.5 && ((RVEF > 0.40 && RVEF <= 0.45) || RVEDVI >= edvi_min)) {
  NMIN = NMIN + 1.0;
}
if (RESIDMYO < 60.0) {
  NMAJ = NMAJ + 1.0;                       // biopsy assumed available
} else if (RESIDMYO < 75.0) {
  NMIN = NMIN + 1.0;
}
if (NTWI >= 3.0)      { NMAJ = NMAJ + 1.0; }
else if (NTWI >= 2.0) { NMIN = NMIN + 1.0; }
if (EPSILONW > 0.5)                          { NMAJ = NMAJ + 1.0; }
else if (TAD_MS >= 55.0 || LATEPOT > 0.5)    { NMIN = NMIN + 1.0; }
double nsvt = (HVA_YR > 0.045) ? 1.0 : 0.0;
if (nsvt > 0.5 && MYO_LV < 0.85 * MYO0) { NMAJ = NMAJ + 1.0; }
else if (nsvt > 0.5 || PVC24 > 500.0)   { NMIN = NMIN + 1.0; }
if (VARIANT > 0.5) { NMAJ = NMAJ + 1.0; }

// 3 = definite, 2 = borderline, 1 = possible, 0 = none
double TFC = 0.0;
if (NMAJ >= 2.0 || (NMAJ >= 1.0 && NMIN >= 2.0) || NMIN >= 4.0)      { TFC = 3.0; }
else if ((NMAJ >= 1.0 && NMIN >= 1.0) || NMIN >= 3.0)                { TFC = 2.0; }
else if (NMAJ >= 1.0 || NMIN >= 2.0)                                 { TFC = 1.0; }

// ---- mechanics, exposed so the clock can be inspected directly ----------
double EXNOWO  = (TIME >= T_RESTR) ? EX2 : EX;
double FEX     = fmin(0.35, EXNOWO / (MET_VIG * 168.0));
double HEFF0   = MYO0 + XI_FIB * FIB0 + XI_FAT * FAT0;
double HRVOUT  = fmax(0.22, (MYO_RV + XI_FIB * FIB_RV + XI_FAT * FAT_RV) / HEFF0);
double RRVOUT  = pow(fmax(1.0E-6, RVEDV / (RVEDVI0 * BSA)), 1.0 / 3.0);
double PRVOUT  = 1.0 + K_PRV * fmax(0.0, RRVOUT * RRVOUT * RRVOUT - 1.0);
double SIGRV   = PRVOUT * RRVOUT / HRVOUT;
double HLVOUT  = fmax(0.22, (MYO_LV + XI_FIB * FIB_LV + XI_FAT * FAT_LV) / HEFF0);
double SIGLV   = pow(fmax(1.0E-6, LVEDV / (LVEDVI0 * BSA)), 1.0 / 3.0) / HLVOUT;
double SRRV    = (PHI_EX > 0.5) ? pow(SIGRV, NMECH) : 1.0;
double SXRV    = (PHI_EX > 0.5) ? pow(SIGRV * (1.0 + K_EX_RV), NMECH) : 1.0;
double SRLV    = (PHI_EX > 0.5) ? pow(SIGLV, NMECH) : 1.0;
double SXLV    = (PHI_EX > 0.5) ? pow(SIGLV * (1.0 + K_EX_LV), NMECH) : 1.0;
double LOADRVO = ((1.0 - FEX) * HR_REST * SRRV + FEX * HR_EX * SXRV) / HR_REST;
double LOADLVO = ((1.0 - FEX) * HR_REST * SRLV + FEX * HR_EX * SXLV) / HR_REST;

$CAPTURE
AGE EXNOWO RVEF LVEF RVEDVI LVEDVI FIBFATRV FIBFATLV RESIDMYO
TAD_MS NTWI EPSILONW LATEPOT
G1OUT G2OUT HVA_YR VAFREE SURV
NMAJ NMIN TFC
SIGRV SIGLV LOADRVO LOADLVO FEX
'

# =============================================================================
#  BUILD
# =============================================================================
mod <- mcode("arvc", arvc_code)

# Convenience: named dose compartments for the event tables below
CMT_BB  <- which(cmt(mod) == "A_BB_G")
CMT_FL  <- which(cmt(mod) == "A_FL_G")
CMT_SO  <- which(cmt(mod) == "A_SO_G")
CMT_AM  <- which(cmt(mod) == "A_AM_G")
CMT_VEC <- which(cmt(mod) == "A_VEC")

# genotype presets -----------------------------------------------------------
#  Note that only PKP2_SET / DSP_SET / KAPPA_LV / VARIANT change.  There is no
#  chamber-specific biology anywhere except KAPPA_LV, and that is ON only for
#  DSP and FLNC.
GENOTYPE <- list(
  WT     = list(PKP2_SET = 1.00, DSP_SET = 1.00, KAPPA_LV = 1.00, VARIANT = 0),
  PKP2tv = list(PKP2_SET = 0.50, DSP_SET = 1.00, KAPPA_LV = 1.00, VARIANT = 1),
  # KAPPA_LV has to exceed LOAD_RV/LOAD_LV (about 3.4 in a competitive
  # athlete) for the LV to actually go first.  Anything less makes DSP a
  # slightly-less-RV-dominant disease, which is not what DSP is.
  DSPtv  = list(PKP2_SET = 0.88, DSP_SET = 0.50, KAPPA_LV = 5.00, VARIANT = 1),
  FLNCtv = list(PKP2_SET = 0.95, DSP_SET = 0.95, KAPPA_LV = 4.60, VARIANT = 1)
)

# exercise presets, MET-hours per week of vigorous-equivalent activity -------
EXERCISE <- c(sedentary = 6, guideline = 15, competitive = 60, elite = 100)

YR <- 365.25

age_to_day <- function(age, age0 = 12) (age - age0) * YR

# ---------------------------------------------------------------------------
#  Chronic oral dosing helper.  Doses run from `start_age` to `stop_age`.
# ---------------------------------------------------------------------------
chronic <- function(cmt, amt, ii, start_age, stop_age, age0 = 12) {
  n <- max(1, floor((stop_age - start_age) * YR / ii))
  ev(time = age_to_day(start_age, age0), amt = amt, cmt = cmt,
     ii = ii, addl = n - 1)
}

# =============================================================================
#  SCENARIO 1 — PENETRANCE IS SET BY LOAD, NOT BY GENOTYPE
#  Same genotype, same equations; only the exercise dose differs.
# =============================================================================
run_natural_history <- function(genotype = "PKP2tv", age_end = 80) {
  bind_rows(lapply(names(EXERCISE), function(ex) {
    mod %>%
      param(GENOTYPE[[genotype]]) %>%
      param(EX = EXERCISE[[ex]]) %>%
      mrgsim(end = age_to_day(age_end), delta = 30) %>%
      as_tibble() %>%
      mutate(genotype = genotype, exercise = ex)
  }))
}

# =============================================================================
#  SCENARIO 2 — THE FALSIFIER
#  PHI_EX = 0 makes damage load-independent.  Every exercise-related
#  prediction should collapse.
# =============================================================================
run_falsifier <- function(genotype = "PKP2tv", age_end = 80) {
  bind_rows(lapply(c(1, 0), function(phi) {
    bind_rows(lapply(c("sedentary", "competitive"), function(ex) {
      mod %>%
        param(GENOTYPE[[genotype]]) %>%
        param(EX = EXERCISE[[ex]], PHI_EX = phi) %>%
        mrgsim(end = age_to_day(age_end), delta = 90) %>%
        as_tibble() %>%
        mutate(PHI_EX = phi, exercise = ex)
    }))
  }))
}

# =============================================================================
#  SCENARIO 3 — LOAD-REDUCING THERAPY (the external validation)
#  Fabritz 2011 PMID 21292134: furosemide + nitrate prevented the phenotype in
#  plakoglobin-deficient mice.  The model was not fitted to this.
# =============================================================================
run_load_reduction <- function(age_end = 60) {
  base <- mod %>% param(GENOTYPE$PKP2tv) %>% param(EX = EXERCISE[["competitive"]])
  arms <- list(
    "no therapy"                 = function(m) m,
    "exercise restriction"       = function(m) param(m, EX2 = EXERCISE[["guideline"]],
                                                     T_RESTR = age_to_day(16)),
    "load-reducing therapy"      = function(m) param(m, ON_LOADRED = 1),
    "both"                       = function(m) param(m, EX2 = EXERCISE[["guideline"]],
                                                     T_RESTR = age_to_day(16),
                                                     ON_LOADRED = 1),
    "beta-blocker only"          = function(m) m
  )
  bind_rows(lapply(names(arms), function(nm) {
    m <- arms[[nm]](base)
    d <- if (nm == "beta-blocker only")
      chronic(CMT_BB, 80, 1, 16, age_end) else ev(time = 0, amt = 0, cmt = CMT_BB)
    m %>% mrgsim(events = d, end = age_to_day(age_end), delta = 30) %>%
      as_tibble() %>% mutate(arm = nm)
  }))
}

# =============================================================================
#  SCENARIO 4 — A DRUG CANNOT OUTPERFORM THE GENERATOR IT OCCUPIES
#  Started at age 38 in overt disease; read the hazard at 39.
# =============================================================================
run_antiarrhythmics <- function(age_start = 38, age_end = 42) {
  base <- mod %>% param(GENOTYPE$PKP2tv) %>% param(EX = EXERCISE[["competitive"]])
  none <- ev(time = 0, amt = 0, cmt = CMT_BB)
  arms <- list(
    "none"                      = list(m = base, d = none),
    "beta-blocker (nadolol)"    = list(m = base,
                                       d = chronic(CMT_BB, 80, 1, age_start, age_end)),
    "sotalol"                   = list(m = base,
                                       d = chronic(CMT_SO, 160, 0.5, age_start, age_end)),
    "flecainide alone"          = list(m = base,
                                       d = chronic(CMT_FL, 100, 0.5, age_start, age_end)),
    "beta-blocker + flecainide" = list(m = base,
                                       d = c(chronic(CMT_BB, 80, 1, age_start, age_end),
                                             chronic(CMT_FL, 100, 0.5, age_start, age_end))),
    "amiodarone"                = list(m = base,
                                       d = chronic(CMT_AM, 200, 1, age_start, age_end)),
    "ablation, endo only"       = list(m = base, d = none, abl = 0.35),
    "ablation, endo + epi"      = list(m = base, d = none, abl = 0.78)
  )
  bind_rows(lapply(names(arms), function(nm) {
    a <- arms[[nm]]
    d <- a$d
    if (!is.null(a$abl)) {
      # ablation is a step change in the substrate state, not a dose
      d <- c(d, ev(time = age_to_day(age_start), amt = a$abl, cmt =
                     which(cmt(mod) == "ABL_HOM")))
    }
    a$m %>% mrgsim(events = d, end = age_to_day(age_end), delta = 15) %>%
      as_tibble() %>% mutate(arm = nm)
  }))
}

# =============================================================================
#  SCENARIO 5 — FLECAINIDE HAS TWO SIGNS
#  The same drug, started at different points on the same trajectory.
# =============================================================================
run_flecainide_sign <- function(ages = c(22, 28, 34, 40, 46)) {
  bind_rows(lapply(ages, function(a) {
    bind_rows(lapply(c(FALSE, TRUE), function(flec) {
      d <- chronic(CMT_BB, 80, 1, 16, a + 1)
      if (flec) d <- c(d, chronic(CMT_FL, 100, 0.5, a, a + 1))
      mod %>% param(GENOTYPE$PKP2tv) %>%
        param(EX = EXERCISE[["competitive"]]) %>%
        mrgsim(events = d, end = age_to_day(a + 1), delta = 15) %>%
        as_tibble() %>% mutate(start_age = a, flecainide = flec)
    }))
  }))
}

# =============================================================================
#  SCENARIO 6 — THE ICD: ALL OF THE MORTALITY, NONE OF THE DISEASE
# =============================================================================
run_icd <- function(age_end = 60) {
  bind_rows(lapply(c(0, 1), function(icd) {
    mod %>% param(GENOTYPE$PKP2tv) %>%
      param(EX = EXERCISE[["competitive"]], ON_ICD = icd) %>%
      mrgsim(events = chronic(CMT_BB, 80, 1, 16, age_end),
             end = age_to_day(age_end), delta = 30) %>%
      as_tibble() %>% mutate(ICD = ifelse(icd > 0, "ICD", "no ICD"))
  }))
}

# =============================================================================
#  SCENARIO 7 — AAV9-PKP2 GENE THERAPY: TIMING DOMINATES DOSE
#  The transgene restores RESERVE, i.e. the rate.  It removes no existing scar.
# =============================================================================
run_gene_therapy <- function(infusion_ages = c(18, 26, 34, 45), age_end = 60,
                             td_target = 0.62) {
  bind_rows(lapply(c(NA, infusion_ages), function(a) {
    d <- chronic(CMT_BB, 80, 1, 16, age_end)
    if (!is.na(a)) d <- c(d, ev(time = age_to_day(a), amt = td_target,
                                cmt = CMT_VEC))
    mod %>% param(GENOTYPE$PKP2tv) %>%
      param(EX = EXERCISE[["competitive"]]) %>%
      mrgsim(events = d, end = age_to_day(age_end), delta = 30) %>%
      as_tibble() %>%
      mutate(infusion = ifelse(is.na(a), "none", paste0("age ", a)))
  }))
}

# =============================================================================
#  SCENARIO 8 — THE UPSTREAM / AMPLIFIER ARMS
#  These act on inflammation and the adipogenic switch, i.e. on the amplifier
#  and not on the clock, so they slow the disease without stopping it.
# =============================================================================
run_upstream <- function(age_end = 60) {
  arms <- list(
    "none"                    = list(),
    "IL-1 blockade"           = list(ON_IL1 = 1),
    "corticosteroid"          = list(ON_GC = 1),
    "MRA"                     = list(ON_MRA = 1),
    "GSK-3beta inhibition"    = list(ON_GSK = 1),
    "all four"                = list(ON_IL1 = 1, ON_MRA = 1, ON_GSK = 1),
    "exercise restriction"    = list(EX2 = EXERCISE[["sedentary"]],
                                     T_RESTR = age_to_day(24))
  )
  bind_rows(lapply(names(arms), function(nm) {
    m <- mod %>% param(GENOTYPE$PKP2tv) %>% param(EX = EXERCISE[["competitive"]])
    if (length(arms[[nm]])) m <- do.call(param, c(list(m), arms[[nm]]))
    m %>% mrgsim(events = chronic(CMT_BB, 80, 1, 16, age_end),
                 end = age_to_day(age_end), delta = 30) %>%
      as_tibble() %>% mutate(arm = nm)
  }))
}

# =============================================================================
#  SCENARIO 9 — CHAMBER SELECTIVITY WITH NO CHAMBER-SPECIFIC BIOLOGY
# =============================================================================
run_chamber <- function(age_end = 50) {
  bind_rows(lapply(c("PKP2tv", "DSPtv", "FLNCtv"), function(g) {
    bind_rows(lapply(c("sedentary", "competitive"), function(ex) {
      mod %>% param(GENOTYPE[[g]]) %>% param(EX = EXERCISE[[ex]]) %>%
        mrgsim(end = age_to_day(age_end), delta = 60) %>%
        as_tibble() %>% mutate(genotype = g, exercise = ex)
    }))
  }))
}

# =============================================================================
#  SCENARIO 10 — PENETRANCE AS A POPULATION QUANTITY
#  Penetrance is a statement about a cohort, so FRAILTY carries the
#  between-subject spread.  Nothing here is fitted; the ~44%-by-40 on guideline
#  activity against ~100%-by-40 in competitive athletes is the prediction.
# =============================================================================
run_penetrance <- function(genotype = "PKP2tv", n = 15, age_end = 75) {
  z <- qnorm((seq_len(n) - 0.5) / n)
  bind_rows(lapply(names(EXERCISE), function(ex) {
    bind_rows(lapply(c(1, 0.72), function(sk) {
      idata <- data.frame(ID = seq_len(n),
                          FRAILTY = exp(0.55 * z),
                          SEX_K = sk,
                          MALE = ifelse(sk > 0.9, 1, 0))
      mod %>% param(GENOTYPE[[genotype]]) %>% param(EX = EXERCISE[[ex]]) %>%
        idata_set(idata) %>%
        mrgsim(end = age_to_day(age_end), delta = 180) %>%
        as_tibble() %>%
        mutate(exercise = ex, sex = ifelse(sk > 0.9, "male", "female"))
    }))
  }))
}

penetrance_table <- function(sim, ages = c(30, 40, 50, 60)) {
  dx <- sim %>%
    group_by(exercise, sex, ID) %>%
    summarise(age_dx = suppressWarnings(min(AGE[TFC >= 3])), .groups = "drop")
  strata <- dx %>% distinct(exercise, sex)
  bind_cols(
    strata,
    as.data.frame(setNames(lapply(ages, function(a) {
      vapply(seq_len(nrow(strata)), function(i) {
        d <- dx$age_dx[dx$exercise == strata$exercise[i] & dx$sex == strata$sex[i]]
        mean(d <= a)
      }, numeric(1))
    }), paste0("by_", ages)))
  )
}

# =============================================================================
#  PLOTS
# =============================================================================
plot_natural_history <- function(sim) {
  sim %>%
    tidyr::pivot_longer(c(RVEF, RVEDVI, FIBFATRV, HVA_YR)) %>%
    ggplot(aes(AGE, value, colour = exercise)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "age (years)", y = NULL,
         title = "ARVC: same genotype, four exercise doses",
         subtitle = "the clock is cumulative mechanical work, not time") +
    theme_bw()
}

plot_generators <- function(sim) {
  sim %>%
    tidyr::pivot_longer(c(G1OUT, G2OUT)) %>%
    ggplot(aes(AGE, value, colour = name)) +
    geom_line(linewidth = 0.8) +
    labs(x = "age (years)", y = "generator drive",
         colour = NULL,
         title = "Two arrhythmia generators on different clocks",
         subtitle = "Generator I is trigger-driven and early; Generator II is scar re-entry") +
    theme_bw()
}

# =============================================================================
#  EXAMPLE SESSION
# =============================================================================
if (interactive()) {
  nh <- run_natural_history()
  print(plot_natural_history(nh))

  nh %>% group_by(exercise) %>%
    summarise(age_definite = suppressWarnings(min(AGE[TFC >= 3])))

  fx <- run_falsifier()
  fx %>% group_by(PHI_EX, exercise) %>%
    summarise(age_definite = suppressWarnings(min(AGE[TFC >= 3])), .groups = "drop")
  # With PHI_EX = 0 the two exercise arms should become indistinguishable.

  lr <- run_load_reduction()
  lr %>% filter(abs(AGE - 45) < 0.2) %>% select(arm, RVEDVI, RVEF, FIBFATRV)

  aa <- run_antiarrhythmics()
  aa %>% filter(abs(AGE - 39) < 0.05) %>% select(arm, G1OUT, G2OUT, HVA_YR)

  gt <- run_gene_therapy()
  gt %>% filter(abs(AGE - 55) < 0.2) %>% select(infusion, FIBFATRV, RVEF, HVA_YR)

  pen <- run_penetrance()
  penetrance_table(pen)
}
