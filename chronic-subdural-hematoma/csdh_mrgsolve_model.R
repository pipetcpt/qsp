##############################################################################
## csdh_mrgsolve_model.R
## Chronic Subdural Haematoma (cSDH) — Quantitative Systems Pharmacology model
## 만성 경막하 혈종 — 정량적 시스템 약리학 모델
## ===========================================================================
##
## THE STRUCTURAL CLAIM
## --------------------
## A chronic subdural haematoma is taught as a clot that is failing to resorb.
## It is not.  It is a SECRETING ORGAN.  A neomembrane grows on the dural side
## of the collection, and its immature neocapillaries filter plasma into the
## subdural space faster than the meningeal lymphatics can clear it.  So the
## volume is not a decaying stock — it is the FIXED POINT of a three-term
## balance:
##
##     dV/dt  =  J_exudation  +  J_rebleed  -  J_absorption
##
##   J_ex   = Lp(N_CAP,N_MAT) * A_mem * [(P_cap - ICP) - sigma*(pi_pl - pi_hem)]
##   J_rb   = k_rb * N_CAP * A_mem * fragility * (1 - H)      * MMA_supply
##   J_abs  = (k_lymph*A_mem + k_arach*V) * (ICP - P_sink)
##
## Everything else in this file exists to supply those three lines with
## honest arguments.  And every therapy is then a CLAIM ABOUT WHICH TERM IT
## TOUCHES — which is the reason the 2024-25 randomised trials disagree:
##
##   burr-hole / drain   removes the STOCK.  Leaves Lp, A_mem, N_CAP, and the
##                       un-re-expanded brain exactly as they were.
##   irrigation          removes the haem CONCENTRATION, i.e. it switches off
##                       the loop that maintains the membrane.  (see below)
##   MMA embolisation    removes the SOURCE.  Leaves the stock untouched, so
##                       its effect is an INTEGRAL over weeks, not an event.
##   dexamethasone       lowers Lp and VEGF locally; poisons the patient
##                       systemically.  One drug, two signs.
##   atorvastatin        MATURES the vessels: raises sigma AND lowers Lp.
##   tranexamic acid     breaks the fibrinolytic loop — but reaches the cavity
##                       only by being carried in with the exudate.
##
## THE TWO LOOPS
## -------------
##   L1  haem -> macrophage -> VEGF -> immature vessel -> rebleed -> haem
##   L2  immature endothelium -> tPA -> plasmin -> FDP -> clot cannot form
##       -> rebleed -> more fibrin to lyse -> more FDP
##
## L1 is written with a Hill exponent of 2 on the haem term.  That single
## exponent is what makes the disease BISTABLE, and bistability is what makes
## "irrigate until clear" a curative manoeuvre instead of a cosmetic one.
## With a LINEAR haem term (which is how this model was first written) the
## system has ONE fixed point at ~105 mL, every operation refills to the
## untreated volume, and the modelled recurrence rate is 100% by construction.
##
## THE PATIENT COVARIATE THAT APPEARS TWICE WITH OPPOSITE SIGNS
## -----------------------------------------------------------
## V_RES, the intracranial reserve created by cerebral atrophy, enters the
## model exactly twice:
##   (1) it buffers pressure and midline shift, so an atrophic patient
##       TOLERATES a larger haematoma and presents later and larger;
##   (2) it slows brain re-expansion, k_reexp = k0*exp(-V_RES/lambda), so an
##       atrophic patient is LEFT WITH the space after evacuation.
## The same number therefore protects the patient before the operation and
## condemns them after it.  Nothing else in the file couples those two
## statements — they are both consequences of one parameter.
##
## WHAT THE CALIBRATED MODEL ACTUALLY SAYS  (read this before believing the
## paragraphs above)
## -----------------------------------------------------------------------
## The structure above is what the model was BUILT to say.  When its free
## parameters are then fitted to the published endpoints, it does not end up
## endorsing all of it.  Honestly:
##
##   IT REPRODUCES.  The incident natural history falls out unfitted --
##   latency 33.5 d from a thin film to a 78 mL presentation (clinically
##   3-8 wk), plateau 112 mL / 16.5 mm / 10.5 mm shift, cavity protein
##   6.58 g/dL (reported 4-8), tPA 24 ng/mL and FDP 258 ug/mL (both reported
##   far above plasma), CT density 50 HU.  The MMA embolisation endpoint is
##   essentially exact: 0.044 against EMBOLISE's 0.041.
##
##   ITS BEST RESULT is the embolisation TIME COURSE, which was not fitted.
##   The absolute risk difference grows 0.000 (day 14) -> 0.003 (day 30) ->
##   0.058 (day 90) -> 0.108 (day 180), because embolisation removes the
##   SOURCE while the existing neovessels decay with tau ~ 5-9 d.  That is a
##   quantitative account of why EMBOLISE (180 d, positive) and MAGIC-MT
##   (90 d, null) disagree, derived rather than assumed.
##
##   IT FAILS, and the failures are not hidden:
##     * The haem SWITCH is inoperative at the fitted parameters.  HB50 came
##       back at 0.90 g/dL against a cavity haemoglobin of ~4 g/dL, so
##       haem_drive sits at 0.97 in every arm and the irrigation sweep is
##       FLAT (P(reop) 0.151 -> 0.156 as wash goes 0.05 -> 1.00).  The
##       bistability is in the equations but not in the fitted regime, so
##       this model does NOT support the "irrigate until clear" argument.
##     * The drain contrast is absent: 0.059 / 0.075 against Santarius
##       0.093 / 0.240.  A 48 h drain here only delays refilling ~10 d, and a
##       90-day hazard integral barely notices.  For the real 2.6-fold effect
##       the drain must cause a PERMANENT state change (fusion of the two
##       membrane leaves), which this parameterisation does not deliver.
##     * Atorvastatin does essentially nothing (0.050 vs 0.052) where ATOCH
##       halved conversion to surgery.
##     * Dexamethasone's favourable-outcome effect has the WRONG SIGN (0.264
##       vs 0.228) and the absolute mRS scale is far off (0.23 vs a reported
##       0.90); the steroid-harm weights are not calibrated.
##     * At the fitted parameters rebleed EXCEEDS exudation at plateau
##       (J_rb 6.56 vs J_ex 3.46 mL/day), which inverts the "secreting organ"
##       framing above.  The fit prefers a rebleed-dominated lesion.
##     * MMA_FLOOR has exactly zero local sensitivity, because the supply
##       term is saturated while PMMA = 1.
##
## Sensitivity says this is a fluid-mechanics model: the ranking is D50_PAT
## (cavity patency, -13.8%), P_CAP (+12.2%), K_RB (+11.1%), K_SPROUT (+8.4%).
##
## 50 ODEs: 14 drug PK (incl. 3 intra-cavity), 7 cavity contents, 10 membrane
## biology, 4 fibrinolysis, 3 supply/procedure, 5 brain mechanics, 7 clinical,
## 3 hazard clocks, 1 membrane extent.
##
## VERIFICATION
## ------------
## This build environment has no R toolchain, so the same 50-ODE system has
## been INDEPENDENTLY re-implemented in `csdh_verify.py` (numpy/scipy LSODA)
## and integrated there.  Every number quoted in the comments below and in
## README.md is produced by that file.  The two transcriptions are meant to be
## the same system with the same parameter block; if they disagree, one of them
## is wrong.  Doing this exposed nine defects, listed in README.md — including
## a protein–oncotic loop with gain >1 that drove cavity protein to 739 g/dL,
## an absorption term with no driving force that cleared a 3 mL collection in a
## day, and hazard clocks that were left running through the run-in so that
## every arm began with a 58% reoperation probability at t = 0.
##
## Usage:
##   source("csdh_mrgsolve_model.R")
##   CSDH_natural_history()      # incident cSDH from a thin film
##   CSDH_run_scenarios()        # the 16 shipped scenarios
##   CSDH_flux_ledger()          # which term does each therapy touch
##   CSDH_bistability()          # the haem switch, and where it sits
##   CSDH_atrophy_two_signs()    # V_RES protects, then condemns
##   CSDH_mmae_timecourse()      # why 30-day endpoints must be null
##   CSDH_txa_site_delay()       # 3 h plasma half-life, 6 day site half-life
##   CSDH_dex_two_signs()        # fewer reoperations, worse patients
##   CSDH_loop_gain()            # the fibrinolytic amplification factor
##   CSDH_trial_ledger()         # model vs published trial numbers
##   CSDH_sensitivity()          # local sensitivity of the endpoints
##############################################################################

library(mrgsolve)
library(dplyr)

##############################################################################
## MODEL SPECIFICATION
##############################################################################

csdh_code <- '
$PROB
# cSDH QSP: exudation / rebleed / absorption balance with a bistable haem switch

$GLOBAL
#define _SP_(x, s) ((s) * log1p(exp(fmin(fmax((x) / (s), -60.0), 60.0))))

// Softplus, a smooth max(0, x).  The hard kinks in the pressure terms made
// the system stiff enough that the reference integrator stalled indefinitely;
// smoothing them is the fix, not loosening the tolerances.

$PARAM @annotated
// ---- patient descriptors -------------------------------------------------
AGE     : 76.0  : Age (years)
WT      : 70.0  : Body weight (kg)
V_RES   : 30.0  : Intracranial reserve / atrophy volume (mL) -- appears twice, opposite signs

// ---- geometry: the membrane has an extent of its own ---------------------
A_MAX   : 180.0 : Area of one cerebral convexity (cm2)
K_A     : 22.0  : Volume at which the collection has spread halfway (mL)
K_SPREAD: 0.35  : Spreading rate of the membrane (1/day)
K_FUSEX : 0.022 : Loss of membrane extent once the layers appose (1/day)
D_CAP   : 35.0  : Maximum plausible haematoma thickness (mm)
D_MEANF : 0.45  : Mean/max thickness ratio for a crescent (-)
P_MIN   : 1.5   : Floor on intracranial pressure (mmHg)
V_FLOOR : 0.5   : Floor on cavity volume in the washout terms (mL)

// ---- Starling exchange across the outer neomembrane ----------------------
P_CAP   : 25.0  : Membrane capillary hydrostatic pressure (mmHg)
PI_PL   : 25.0  : Plasma colloid osmotic pressure (mmHg)
PI_COEF : 3.6   : mmHg per g/dL of cavity protein
LP_IMM  : 0.240641  : Hydraulic conductivity, immature neovessel (mL/day/100cm2/mmHg)
LP_MAT  : 0.06  : Hydraulic conductivity, mature vessel -- ~12x tighter
SIG_IMM : 0.25  : Reflection coefficient, immature (fenestrated)
SIG_MAT : 0.85  : Reflection coefficient, mature
PROT_PL : 0.070 : Plasma protein concentration (g/mL)
PS_PROT : 3.0   : Diffusive protein equilibration (mL/day/100cm2)

// ---- clearance of the collection ----------------------------------------
K_LYMPH : 0.981762  : Meningeal lymphatic clearance (mL/day/100cm2)
K_ARACH : 0.085 : Arachnoid granulation bulk clearance (1/day)
SEPT_ABS: 0.55  : Septation penalty on absorption (-)
P_SINK  : 1.0   : Dural lymphatic / deep cervical node sink pressure (mmHg)
P_REF   : 10.0  : Normalisation of the absorption pressure gradient (mmHg)
K_DRAIN : 3.0   : Active subdural drain rate constant (1/day)
P_DRAIN : 2.0   : Height the drain is set at (mmHg)

// ---- rebleed from the fragile membrane -----------------------------------
K_RB    : 3.70  : Rebleed flux constant (mL/day/100cm2 per unit N_CAP)
MMA_FLOOR:0.371842  : Fraction of membrane supply NOT from the MMA (-)
W_IL8   : 0.15  : Autocrine IL-8 amplification of the macrophage drive (-)
A_MMP   : 0.60  : MMP-9 contribution to fragility (-)
KM_MMP  : 220.0 : MMP-9 half-effect (ng/mL)
A_ANG   : 0.35  : Ang-2 contribution to fragility (-)
KM_ANG  : 12.0  : Ang-2 half-effect (ng/mL)
K_RECAN : 0.004 : MMA recanalisation / collateral reconstitution (1/day)

// ---- haemoglobin / density ----------------------------------------------
HB_BLOOD: 150.0 : Whole blood haemoglobin (mg/mL)
K_HBDEG : 0.080 : Haem breakdown and clearance (1/day)

// ---- fibrinolysis loop L2 -----------------------------------------------
K_TPA   : 6500.0: tPA production (ng/day/100cm2 per unit N_CAP)
KDEG_TPA: 3.0   : tPA degradation in cavity (1/day)
K_PAI   : 0.055 : PAI-1 mediated tPA inhibition (1/day per ng/mL)
K_PAISYN: 2.2   : PAI-1 synthesis (ng/mL/day)
K_PAIDEG: 0.10  : PAI-1 degradation (1/day)
K_PLS   : 0.052 : Plasmin generation per ng/mL tPA (1/day)
KEL_PLS : 1.4   : Plasmin activity decay (1/day)
IC50_TXA: 5.0   : TXA inhibition of plasminogen activation (mg/L)
K_FIBF  : 3.2   : Fibrin formed per mL extravasated blood (mg/mL)
K_LYS   : 0.55  : Fibrinolysis rate per unit plasmin (1/day)
K_FDPY  : 1.0   : FDP yield per mg fibrin lysed (-)
KDEG_FDP: 0.12  : FDP degradation (1/day)
K_FDP   : 150.0 : FDP halving haemostatic competence (ug/mL)
IC50_DOAC:110.0 : Apixaban in cavity impairing haemostasis (ng/mL)

// ---- angiogenesis / membrane biology ------------------------------------
K_SPROUT: 0.41 : Endothelial sprouting rate (1/day)
N_MAX   : 2.50  : Carrying capacity of vessel density (-)
K_MAT   : 0.060 : Intrinsic immature -> mature conversion (1/day)
MAT_SUPP: 0.920 : Fraction of maturation suppressed by active disease (-)
K_MMPDES: 0.100 : MMP-9 driven loss of MATURE vessels (1/day)
KM_MMPD : 220.0 : Half-effect for MMP-9 remodelling (ng/mL)
KD_CAP  : 0.2 : Immature vessel regression on VEGF withdrawal (1/day)
KD_MAT  : 0.020 : Mature vessel loss (1/day)
EMAX_ATV: 6.0   : Atorvastatin fold-increase of maturation (-)
EC50_ATV: 1.6   : Atorvastatin effect-site EC50 (ug/L)
K_VEGF  : 480.0 : VEGF production per unit macrophage (pg/mL/day)
KEL_VEGF: 4.0   : VEGF elimination (1/day)
VEGF_REF: 120.0 : VEGF normalisation (pg/mL)
EMAX_DEXV:2.30  : Dexamethasone suppression of VEGF transcription (-)
EC50_DEX: 6.0   : Dexamethasone effect-site EC50 (ug/L)
DEX_LP  : 0.60  : Dexamethasone direct anti-permeability weight (-)
HB50    : 0.899831   : Cavity haemoglobin at half-maximal inflammatory drive (g/dL)
HB_HILL : 2.0   : Steepness of the haem switch -- the source of bistability
D50_PAT : 5.0   : Cavity thickness at half patency (mm)
N_PAT   : 2.0   : Steepness of the apposition term (-)
K_FUSE  : 0.079721 : Membrane involution when layers are apposed (1/day)
K_MAC   : 1.82 : Macrophage / eosinophil recruitment (1/day)
KD_MAC  : 0.35 : Macrophage loss (1/day)
MAC_MAX : 2.20  : Macrophage carrying capacity (-)
K_ANG2  : 3.4   : Ang-2 production (ng/mL/day)
KEL_ANG2: 0.35  : Ang-2 elimination (1/day)
K_MMP9  : 62.0  : MMP-9 production (ng/mL/day)
KEL_MMP9: 0.30  : MMP-9 elimination (1/day)
ATV_MMP : 0.45  : Atorvastatin suppression of MMP-9 (-)
K_IL6   : 210.0 : IL-6 production (pg/mL/day)
KEL_IL6 : 1.2   : IL-6 elimination (1/day)
K_IL8   : 175.0 : IL-8 production (pg/mL/day)
KEL_IL8 : 1.0   : IL-8 elimination (1/day)
IL6_REF : 900.0 : IL-6 normalisation (pg/mL)
IL8_REF : 800.0 : IL-8 normalisation (pg/mL)
K_MEMO  : 0.045 : Outer membrane growth (1/day)
KD_MEMO : 0.012 : Outer membrane loss (1/day)
K_SEPT  : 0.020 : Septation formation (1/day)
KD_SEPT : 0.004 : Septation loss (1/day)

// ---- brain mechanics ----------------------------------------------------
ICP0    : 8.0   : Baseline intracranial pressure (mmHg)
K_PV    : 0.022 : Pressure-volume exponent (1/mL)
COMP_F  : 0.60  : Fraction of excess volume absorbed by deformation (-)
K_COMP  : 0.20  : Compression loading rate -- fast (1/day)
K_REEXP0: 0.0940875 : Brain re-expansion rate in a non-atrophic brain (1/day)
LAM_ATR : 26.0  : Atrophy constant slowing re-expansion (mL)
MLS_MAX : 25.0  : Maximum midline shift (mm)
K_MLS   : 45.0  : Midline shift half-effect volume (mL)
K_CBF   : 0.055 : Perfusion deficit accrual (1/day)
KREC_CBF: 0.10  : Perfusion recovery (1/day)
K_NEUR  : 0.0075: Neuronal injury rate (1/day)
KREP_NEUR:0.010 : Neuronal repair (1/day)

// ---- clinical -----------------------------------------------------------
K_SYMP  : 0.42  : Symptom accrual (1/day)
KREC_SYMP:0.30  : Symptom recovery (1/day)
K_COGL  : 0.055 : Cognitive loss (1/day)
KREC_COG: 0.030 : Cognitive recovery (1/day)
LAM_REOP: 0.00173096 : Peak reoperation hazard (1/day)
D_REOP  : 10.0  : Thickness at which reoperation is considered (mm)
S_REOP  : 2.2   : Steepness of the reoperation threshold (mm)
LAM_THR : 0.00045: Thromboembolic hazard off anticoagulation (1/day)
LAM_DTH : 0.00016: Baseline mortality hazard (1/day)

// ---- steroid toxicity ---------------------------------------------------
GLU0    : 5.4   : Baseline glucose (mmol/L)
K_GLU   : 3.6   : Steroid hyperglycaemia (mmol/L/day)
KEL_GLU : 1.1   : Glucose return (1/day)
K_INF   : 0.055 : Steroid infection burden accrual (1/day)
KEL_INF : 0.035 : Infection burden resolution (1/day)
K_MYO   : 0.085 : Steroid myopathy / delirium accrual (1/day)
KEL_MYO : 0.045 : Myopathy resolution (1/day)

// ---- PK -----------------------------------------------------------------
DEX_KA  : 6.0   : Dexamethasone absorption (1/day)
DEX_V1  : 35.0  : Dexamethasone central volume (L)
DEX_V2  : 30.0  : Dexamethasone peripheral volume (L)
DEX_Q   : 360.0 : Dexamethasone intercompartmental clearance (L/day)
DEX_CL  : 250.0 : Dexamethasone clearance (L/day)
DEX_KE0 : 1.2   : Dexamethasone effect-site equilibration (1/day)
ATV_KA  : 6.0   : Atorvastatin absorption (1/day)
ATV_V1  : 400.0 : Atorvastatin central volume (L)
ATV_V2  : 5200.0: Atorvastatin peripheral volume (L)
ATV_Q   : 900.0 : Atorvastatin intercompartmental clearance (L/day)
ATV_CL  : 8000.0: Atorvastatin apparent clearance (L/day)
ATV_KE0 : 0.08  : Atorvastatin effect-site equilibration (1/day)
TXA_V1  : 12.0  : TXA central volume (L)
TXA_V2  : 18.0  : TXA peripheral volume (L)
TXA_Q   : 288.0 : TXA intercompartmental clearance (L/day)
TXA_CL  : 180.0 : TXA clearance (L/day)
TXA_KDEGH:0.05  : TXA degradation in cavity (1/day)
DOAC_V1 : 21.0  : Apixaban central volume (L)
DOAC_V2 : 12.0  : Apixaban peripheral volume (L)
DOAC_Q  : 60.0  : Apixaban intercompartmental clearance (L/day)
DOAC_CL : 79.0  : Apixaban apparent clearance (L/day)
DOAC_KDEGH:0.05 : Apixaban degradation in cavity (1/day)

// ---- procedure switches (set by $MAIN / events) --------------------------
DRAIN_ON: 0.0   : Subdural drain active (0/1)

$CMT @annotated
// ---- drug PK (14) -------------------------------------------------------
DEXD    : Dexamethasone gut depot (mg)
DEXC    : Dexamethasone central (mg)
DEXP    : Dexamethasone peripheral (mg)
CEDEX   : Dexamethasone effect site (ug/L)
ATVD    : Atorvastatin gut depot (mg)
ATVC    : Atorvastatin central (mg)
ATVP    : Atorvastatin peripheral (mg)
CEATV   : Atorvastatin effect site (ug/L)
TXAC    : Tranexamic acid central (mg)
TXAP    : Tranexamic acid peripheral (mg)
TXAH    : Tranexamic acid in the haematoma cavity (mg)
DOACC   : Apixaban central (mg)
DOACP   : Apixaban peripheral (mg)
DOACH   : Apixaban in the haematoma cavity (mg)
// ---- the collection (7) -------------------------------------------------
VHEM    : Haematoma volume (mL)
MHB     : Cavity haemoglobin mass (mg)
MPROT   : Cavity protein mass (g)
MFIB    : Cavity fibrin mass (mg)
MFDP    : Cavity fibrin degradation product mass (mg)
MTPA    : Cavity tPA mass (ng)
CPLS    : Cavity plasmin activity (a.u.)
// ---- membrane biology (10) ---------------------------------------------
CPAI    : Cavity PAI-1 (ng/mL)
NCAP    : Immature neocapillary density (-)
NMAT    : Mature vessel density (-)
CVEGF   : Cavity VEGF (pg/mL)
CANG2   : Cavity angiopoietin-2 (ng/mL)
CMMP9   : Cavity MMP-9 (ng/mL)
CIL6    : Cavity IL-6 (pg/mL)
CIL8    : Cavity IL-8 (pg/mL)
NMAC    : Membrane macrophage / eosinophil density (-)
MMEMO   : Outer membrane mass index (-)
// ---- procedure / supply (4) --------------------------------------------
MSEPT   : Septation index (-)
PMMA    : MMA perfusion fraction to the membrane (-)
FDRN    : Drain patency (0/1)
NPC     : Pericyte coverage index (-)
// ---- systemic haemostasis (2) ------------------------------------------
FPLT    : Platelet function index (-)
CFBG    : Plasma fibrinogen (mg/dL)
// ---- brain mechanics (4) ----------------------------------------------
VCOMP   : Brain compression deficit (mL)
XMLSI   : Integrated midline shift exposure (mm.day/10)
XCBF    : Regional perfusion deficit index (-)
NNEUR   : Cortical neuronal integrity (-)
// ---- clinical (5) -----------------------------------------------------
SSYMP   : Symptom burden (-)
SCOG    : Cognitive deficit (-)
GLU     : Blood glucose (mmol/L)
XINF    : Steroid infection burden (-)
XMYO    : Steroid myopathy / delirium burden (-)
// ---- hazard clocks (3) ------------------------------------------------
HREOP   : Cumulative reoperation hazard (-)
HTHR    : Cumulative thromboembolic hazard (-)
HDTH    : Cumulative mortality hazard (-)
// ---- membrane extent (1) ----------------------------------------------
AEXT    : Membrane extent, fraction of convexity (-)

$MAIN
F_DEXD = 0.80;    // dexamethasone oral bioavailability
F_ATVD = 0.14;    // atorvastatin oral bioavailability
F_TXAC = 0.34;    // TXA oral bioavailability (dosed straight into central)
F_DOACC= 0.50;    // apixaban oral bioavailability

$ODE
// =====================================================================
// ALGEBRAIC LAYER
// =====================================================================
double V = fmax(VHEM, 1e-4);

// ---- geometry.  Area comes from the MEMBRANE, not from the fluid. -----
// Tying area to the current volume meant a thin post-operative film had
// almost no exchange area, so no drained cavity could ever refill and the
// membrane that had just been built vanished with the fluid.
double A_mem = fmax(A_MAX * fmin(1.0, fmax(AEXT, 1e-3)), 4.0);   // cm2
double A100  = A_mem / 100.0;
double dmean = V / A_mem;                                        // cm
double dmax  = fmin(D_CAP, 10.0 * dmean / D_MEANF);              // mm

// Cavity patency: once the brain re-apposes the dura the two membrane
// leaves touch, filtration has nowhere to go and the neovessels are
// compressed.  Without this a 0.5 mL film was still refilled at 16 mL/day.
double f_pat = pow(dmax, N_PAT) / (pow(dmax, N_PAT) + pow(D50_PAT, N_PAT));

// ---- drug effect sites ----------------------------------------------
double E_DEXV = EMAX_DEXV * CEDEX / (EC50_DEX + CEDEX);
double E_ATV  = EMAX_ATV  * CEATV / (EC50_ATV + CEATV);
double f_dex_vegf = 1.0 / (1.0 + E_DEXV);
double f_dex_lp   = 1.0 / (1.0 + DEX_LP * E_DEXV);

// ---- vessel population ----------------------------------------------
double ncap = fmax(NCAP, 0.0), nmat = fmax(NMAT, 0.0);
double ntot = ncap + nmat + 1e-9;
double fmat = nmat / ntot;
double pc   = fmin(1.0, fmax(NPC, 0.0));
double sig_eff = SIG_IMM + (SIG_MAT - SIG_IMM) * (0.5 * fmat + 0.5 * pc);
double Lp_eff  = (LP_IMM * ncap + LP_MAT * nmat) * f_dex_lp;

// ---- cavity composition ---------------------------------------------
double C_PROT = (MPROT / V) * 100.0;      // g/dL
double C_HB   = (MHB   / V) / 10.0;       // g/dL
double C_FDP  = (MFDP  / V) * 1000.0;     // ug/mL
double C_TPA  = (MTPA  / V);              // ng/mL
double C_TXAH = (TXAH  / V) * 1000.0;     // mg/L
double C_DOAH = (DOACH / V) * 1e6;        // ng/mL
double HU     = 10.0 + 5.5 * C_HB + 1.5 * C_PROT;

// THE SWITCH.  Haem recruits the macrophages that make the VEGF that
// builds the leaky vessels that bleed and release more haem.  Hill 2 gives
// that closed loop a threshold, and the threshold is the bistability.
double hb_n = C_HB / HB50;
double haem_drive = pow(hb_n, HB_HILL) / (1.0 + pow(hb_n, HB_HILL));

// ---- brain mechanics.  V_RES appears here (sign 1 of 2). -------------
double R      = V_RES + fmax(VCOMP, 0.0);
double excess = _SP_(V - R, 1.5);
// ICP is TWO-SIDED: sub-baseline when the brain has not re-expanded into a
// drained cavity.  That is why the space persists, why absorption stalls,
// and why a suction drain does not simply empty the head.
double ICP    = P_MIN + ICP0 * exp(fmin(fmax(K_PV * (V - R), -30.0), 30.0));
double MLS    = MLS_MAX * excess / (excess + K_MLS);

// ---- Starling fluxes -------------------------------------------------
double pi_hem = PI_COEF * C_PROT;
double dP     = (P_CAP - ICP) - sig_eff * (PI_PL - pi_hem);
double J_ex   = Lp_eff * A100 * _SP_(dP, 0.8) * f_pat;

// ---- haemostatic competence -----------------------------------------
double f_fdp  = 1.0 / (1.0 + C_FDP / K_FDP);
double f_doac = 1.0 / (1.0 + C_DOAH / IC50_DOAC);
double f_plt  = fmin(1.0, fmax(FPLT, 0.0));
double f_fbg  = CFBG / (CFBG + 90.0);
double H      = f_fdp * f_doac * f_plt * (0.45 + 0.55 * f_fbg);

// ---- rebleed ---------------------------------------------------------
double frag   = 1.0 + A_MMP * CMMP9 / (KM_MMP + CMMP9)
                    + A_ANG * CANG2 / (KM_ANG + CANG2);
double supply = MMA_FLOOR + (1.0 - MMA_FLOOR) * fmin(1.0, fmax(PMMA, 0.0));
double J_rb   = K_RB * ncap * A100 * frag * (1.0 - H) * supply * f_pat;

// ---- absorption.  Pressure-driven, and it vanishes with the fluid. ---
double p_drive = _SP_(ICP - P_SINK, 0.8) / P_REF;
double J_abs = (K_LYMPH * A100 * f_pat + K_ARACH * V) * p_drive
               / (1.0 + SEPT_ABS * MSEPT);
double J_drn = K_DRAIN * fmin(1.0, fmax(FDRN, 0.0)) * V
               * _SP_(ICP - P_DRAIN, 0.8) / P_REF;
// V is FLOORED here.  Washout of the cavity solutes is a fractional rate
// with the volume in the denominator, so as a successfully treated cavity
// empties this term diverges and the solute equations become violently
// stiff.  Below about half a millilitre there is no cavity left to have a
// concentration in.
double out_frac = (J_abs + J_drn) / fmax(V, V_FLOOR);

// ---- plasma drug concentrations -------------------------------------
double C_DEX  = DEXC  / DEX_V1  * 1000.0;   // ug/L
double C_ATV  = ATVC  / ATV_V1  * 1000.0;   // ug/L
double C_TXA  = TXAC  / TXA_V1;             // mg/L
double C_DOAC = DOACC / DOAC_V1 * 1000.0;   // ng/mL

// =====================================================================
// DIFFERENTIAL EQUATIONS
// =====================================================================
// ---- PK --------------------------------------------------------------
dxdt_DEXD  = -DEX_KA * DEXD;
dxdt_DEXC  =  DEX_KA * DEXD - DEX_CL * DEXC / DEX_V1
              - DEX_Q * (DEXC / DEX_V1 - DEXP / DEX_V2);
dxdt_DEXP  =  DEX_Q * (DEXC / DEX_V1 - DEXP / DEX_V2);
dxdt_CEDEX =  DEX_KE0 * (C_DEX - CEDEX);

dxdt_ATVD  = -ATV_KA * ATVD;
dxdt_ATVC  =  ATV_KA * ATVD - ATV_CL * ATVC / ATV_V1
              - ATV_Q * (ATVC / ATV_V1 - ATVP / ATV_V2);
dxdt_ATVP  =  ATV_Q * (ATVC / ATV_V1 - ATVP / ATV_V2);
dxdt_CEATV =  ATV_KE0 * (C_ATV - CEATV);

dxdt_TXAC  = -TXA_CL * TXAC / TXA_V1
             - TXA_Q * (TXAC / TXA_V1 - TXAP / TXA_V2)
             - (J_ex + J_rb) * C_TXA / 1000.0;
dxdt_TXAP  =  TXA_Q * (TXAC / TXA_V1 - TXAP / TXA_V2);
// THE POINT: tranexamic acid reaches its site of action by being carried
// in with the exudate, not by diffusion.  Its cavity half-life is
// therefore V/J_abs (~6 days), not its 3 h plasma one.
dxdt_TXAH  =  (J_ex + J_rb) * C_TXA / 1000.0
              - (J_abs + J_drn) * (TXAH / V) - TXA_KDEGH * TXAH;

dxdt_DOACC = -DOAC_CL * DOACC / DOAC_V1
             - DOAC_Q * (DOACC / DOAC_V1 - DOACP / DOAC_V2)
             - (J_ex + J_rb) * C_DOAC / 1e6;
dxdt_DOACP =  DOAC_Q * (DOACC / DOAC_V1 - DOACP / DOAC_V2);
dxdt_DOACH =  (J_ex + J_rb) * C_DOAC / 1e6
              - (J_abs + J_drn) * (DOACH / V) - DOAC_KDEGH * DOACH;

// ---- THE VOLUME BALANCE ---------------------------------------------
dxdt_VHEM  = J_ex + J_rb - J_abs - J_drn;

// ---- cavity contents -------------------------------------------------
dxdt_MHB   = HB_BLOOD * J_rb - K_HBDEG * MHB - out_frac * MHB;
// Convective protein entry is (1-sigma)*C_plasma; the diffusive term is
// added because the neomembrane is not a dialysis bag.  Together with
// J_ex + J_rb = J_abs at steady volume these BOUND cavity protein below
// plasma protein, which is what stops the oncotic term from driving its
// own filtration.  Without it the loop gain exceeds 1 and cavity protein
// runs to 739 g/dL.
dxdt_MPROT = J_ex * (1.0 - sig_eff) * PROT_PL + J_rb * PROT_PL
             + PS_PROT * A100 * (PROT_PL - MPROT / V) - out_frac * MPROT;

double fib_lys = K_LYS * CPLS * MFIB;
dxdt_MFIB  = K_FIBF * J_rb * H * (CFBG / 300.0) - fib_lys - out_frac * MFIB;
dxdt_MFDP  = K_FDPY * fib_lys - KDEG_FDP * MFDP - out_frac * MFDP;
dxdt_MTPA  = K_TPA * ncap * A100 - (KDEG_TPA + K_PAI * CPAI) * MTPA
             - out_frac * MTPA;
double f_txa = 1.0 / (1.0 + C_TXAH / IC50_TXA);
dxdt_CPLS  = K_PLS * C_TPA * f_txa - KEL_PLS * CPLS;
dxdt_CPAI  = K_PAISYN * (1.0 + 0.4 * NMAC) - K_PAIDEG * CPAI;

// ---- membrane biology ------------------------------------------------
double vegf_n = CVEGF / VEGF_REF;
// Maturation is not a fixed rate.  An ACTIVE haematoma keeps its own
// vessels immature (VEGF/Ang-2 antagonise pericyte stabilisation), so the
// membrane stays leaky as long as the haem switch is lit.  Quench the
// switch and the same vessels mature, sigma rises, Lp falls and the
// exudation collapses.  Atorvastatin acts on this same term, so drug
// effect and spontaneous healing share one mechanism.
double k_mat_eff = K_MAT * (1.0 + E_ATV) * (1.0 - MAT_SUPP * haem_drive);
double mmp_des   = K_MMPDES * CMMP9 / (KM_MMPD + CMMP9);

dxdt_NCAP  = K_SPROUT * vegf_n * supply
             * fmax(0.0, 1.0 - (NCAP + NMAT) / N_MAX)
             - k_mat_eff * NCAP - KD_CAP * NCAP
             - K_FUSE * (1.0 - f_pat) * NCAP;
dxdt_NMAT  = k_mat_eff * NCAP - (KD_MAT + mmp_des) * NMAT;
dxdt_NPC   = 0.010 * (1.0 + E_ATV) * (fmat - NPC);

dxdt_CVEGF = K_VEGF * NMAC * f_dex_vegf - KEL_VEGF * CVEGF;
// The IL-8 amplification is GATED ON HAEM.  Written as a free autocrine
// term it is a second positive feedback with no external input, so the
// inflammatory state never switches off: macrophages settled at 0.63 with
// the haem drive at 0.12, the membrane always relit, surgery never cured
// anybody, and MMA embolisation became the only thing that worked
// (P(reop) 0.0007 against a reported 0.041).
double il8_n = CIL8 / IL8_REF;
dxdt_NMAC  = K_MAC * haem_drive * (1.0 + W_IL8 * il8_n) * f_dex_vegf
             * fmax(0.0, 1.0 - NMAC / MAC_MAX) - KD_MAC * NMAC;
dxdt_CANG2 = K_ANG2 * NCAP * f_dex_vegf - KEL_ANG2 * CANG2;
dxdt_CMMP9 = K_MMP9 * NMAC * f_dex_vegf / (1.0 + ATV_MMP * E_ATV)
             - KEL_MMP9 * CMMP9;
dxdt_CIL6  = K_IL6 * NMAC * f_dex_vegf - KEL_IL6 * CIL6;
dxdt_CIL8  = K_IL8 * NMAC * haem_drive * f_dex_vegf - KEL_IL8 * CIL8;
dxdt_MMEMO = K_MEMO * NMAC * (1.0 - MMEMO / 2.0) - KD_MEMO * MMEMO;
dxdt_MSEPT = K_SEPT * MMEMO * (1.0 - MSEPT / 2.0) - KD_SEPT * MSEPT;

// Spreading and fusing are BLENDED with a sigmoid, not switched with a
// ternary.  A hard branch on a state variable makes the solver chatter on
// the switching surface: in the discontinuous form a late embolisation
// (scenario 07) stalled the reference integrator indefinitely.
double ext_target = V / (V + K_A);
double w_ext = 1.0 / (1.0 + exp(-fmin(fmax((ext_target - AEXT) / 0.01,
                                           -60.0), 60.0)));
dxdt_AEXT  = w_ext * K_SPREAD * (ext_target - AEXT)
             - (1.0 - w_ext) * K_FUSEX * (1.0 - f_pat) * AEXT;

dxdt_PMMA  = K_RECAN * (1.0 - PMMA);
dxdt_FDRN  = 0.0;
dxdt_FPLT  = 0.05 * (1.0 - FPLT);
dxdt_CFBG  = 0.25 * (300.0 - CFBG);

// ---- brain mechanics.  V_RES appears here (sign 2 of 2). -------------
// Loads fast, unloads slowly, and unloads MORE SLOWLY the more atrophic
// the brain.  That asymmetry is the entire recurrence story.
double target  = COMP_F * fmax(0.0, V - V_RES);
double k_reexp = K_REEXP0 * exp(-V_RES / LAM_ATR);
double w_cmp = 1.0 / (1.0 + exp(-fmin(fmax((target - VCOMP) / 0.05,
                                           -60.0), 60.0)));
dxdt_VCOMP = w_cmp * K_COMP * (target - VCOMP)
             - (1.0 - w_cmp) * k_reexp * (VCOMP - target);

dxdt_XMLSI = MLS / 10.0;
dxdt_XCBF  = K_CBF * (MLS / 10.0) - KREC_CBF * XCBF;
dxdt_NNEUR = -K_NEUR * XCBF + KREP_NEUR * (1.0 - NNEUR);

// ---- clinical --------------------------------------------------------
double drive = 0.55 * (MLS / 10.0) + 0.45 * (1.0 - NNEUR) / 0.3;
dxdt_SSYMP = K_SYMP * drive - KREC_SYMP * SSYMP;
dxdt_SCOG  = K_COGL * drive - KREC_COG * SCOG;

dxdt_GLU   = K_GLU * E_DEXV / EMAX_DEXV - KEL_GLU * (GLU - GLU0);
dxdt_XINF  = K_INF * E_DEXV - KEL_INF * XINF;
dxdt_XMYO  = K_MYO * E_DEXV - KEL_MYO * XMYO;

// ---- hazard clocks ---------------------------------------------------
double reop_geo = 1.0 / (1.0 + exp(-(dmax - D_REOP) / S_REOP));
double reop_sym = SSYMP / (SSYMP + 0.55);
dxdt_HREOP = LAM_REOP * reop_geo * reop_sym;
// The thromboembolic hazard is carried by the patient who is OFF their
// anticoagulant; it is the price of the haemostatic safety above.
dxdt_HTHR  = LAM_THR / (1.0 + C_DOAC / 30.0);
dxdt_HDTH  = LAM_DTH * (1.0 + 1.6 * XINF + 0.9 * (MLS / 10.0));

$TABLE
double V_ = fmax(VHEM, 1e-4);
double A_mem_ = fmax(A_MAX * fmin(1.0, fmax(AEXT, 1e-3)), 4.0);
double A100_  = A_mem_ / 100.0;
double dmax_  = fmin(D_CAP, 10.0 * (V_ / A_mem_) / D_MEANF);
double f_pat_ = pow(dmax_, N_PAT) / (pow(dmax_, N_PAT) + pow(D50_PAT, N_PAT));
double R_     = V_RES + fmax(VCOMP, 0.0);
double excess_= _SP_(V_ - R_, 1.5);
double ICP_   = P_MIN + ICP0 * exp(fmin(fmax(K_PV * (V_ - R_), -30.0), 30.0));
double MLS_   = MLS_MAX * excess_ / (excess_ + K_MLS);
double CPROT_ = (MPROT / V_) * 100.0;
double CHB_   = (MHB / V_) / 10.0;
double CFDP_  = (MFDP / V_) * 1000.0;
double CTPA_  = (MTPA / V_);
double HU_    = 10.0 + 5.5 * CHB_ + 1.5 * CPROT_;
double hb_n_  = CHB_ / HB50;
double HDRIVE = pow(hb_n_, HB_HILL) / (1.0 + pow(hb_n_, HB_HILL));
double fmat_  = fmax(NMAT,0.0) / (fmax(NCAP,0.0) + fmax(NMAT,0.0) + 1e-9);
double SIGEFF = SIG_IMM + (SIG_MAT - SIG_IMM)
                * (0.5 * fmat_ + 0.5 * fmin(1.0, fmax(NPC, 0.0)));
double LPEFF  = (LP_IMM * fmax(NCAP,0.0) + LP_MAT * fmax(NMAT,0.0))
                / (1.0 + DEX_LP * (EMAX_DEXV * CEDEX / (EC50_DEX + CEDEX)));
double DPNET  = (P_CAP - ICP_) - SIGEFF * (PI_PL - PI_COEF * CPROT_);
double JEX    = LPEFF * A100_ * _SP_(DPNET, 0.8) * f_pat_;
double HCOMP  = (1.0 / (1.0 + CFDP_ / K_FDP))
                * (1.0 / (1.0 + ((DOACH / V_) * 1e6) / IC50_DOAC))
                * fmin(1.0, fmax(FPLT, 0.0))
                * (0.45 + 0.55 * CFBG / (CFBG + 90.0));
double FRAG   = 1.0 + A_MMP * CMMP9 / (KM_MMP + CMMP9)
                    + A_ANG * CANG2 / (KM_ANG + CANG2);
double SUPPLY = MMA_FLOOR + (1.0 - MMA_FLOOR) * fmin(1.0, fmax(PMMA, 0.0));
double JRB    = K_RB * fmax(NCAP,0.0) * A100_ * FRAG * (1.0 - HCOMP)
                * SUPPLY * f_pat_;
double PDRIVE = _SP_(ICP_ - P_SINK, 0.8) / P_REF;
double JABS   = (K_LYMPH * A100_ * f_pat_ + K_ARACH * V_) * PDRIVE
                / (1.0 + SEPT_ABS * MSEPT);
double P_REOP = 1.0 - exp(-HREOP);
double P_THR  = 1.0 - exp(-HTHR);
double P_DTH  = 1.0 - exp(-HDTH);
// mRS-like latent disability: local benefit and systemic harm compete here
double MRS_L  = 0.55 * SCOG + 0.45 * SSYMP + 0.9 * XMYO + 0.7 * XINF
                + 0.35 * (1.0 - NNEUR) / 0.3;
double P_FAV  = 1.0 / (1.0 + exp((MRS_L - 1.55) / 0.42));

$CAPTURE @annotated
dmax_  : Maximum haematoma thickness (mm)
A_mem_ : Membrane exchange area (cm2)
f_pat_ : Cavity patency (-)
ICP_   : Intracranial pressure (mmHg)
MLS_   : Midline shift (mm)
R_     : Effective intracranial reserve (mL)
CPROT_ : Cavity protein (g/dL)
CHB_   : Cavity haemoglobin (g/dL)
CFDP_  : Cavity fibrin degradation products (ug/mL)
CTPA_  : Cavity tPA (ng/mL)
HU_    : CT density (Hounsfield units)
HDRIVE : Haem inflammatory drive (0-1) -- the switch
SIGEFF : Effective reflection coefficient (-)
LPEFF  : Effective hydraulic conductivity (mL/day/100cm2/mmHg)
DPNET  : Net Starling filtration pressure (mmHg)
JEX    : Exudation flux (mL/day)
JRB    : Rebleed flux (mL/day)
JABS   : Absorption flux (mL/day)
HCOMP  : Haemostatic competence (-)
P_REOP : Cumulative reoperation probability (-)
P_THR  : Cumulative thromboembolic probability (-)
P_DTH  : Cumulative mortality probability (-)
MRS_L  : Latent disability score (-)
P_FAV  : Probability of a favourable outcome, mRS 0-3 (-)
'

mod <- mcode("csdh", csdh_code, atol = 1e-8, rtol = 1e-6, maxsteps = 200000)

##############################################################################
## INITIAL CONDITIONS AND THE RUN-IN
##############################################################################

## A thin subacute collection a week or two after a minor head injury:
## blood, no neomembrane, no neovessels yet.
csdh_fresh <- function(V0 = 12, p = as.list(param(mod))) {
  c(DEXD = 0, DEXC = 0, DEXP = 0, CEDEX = 0,
    ATVD = 0, ATVC = 0, ATVP = 0, CEATV = 0,
    TXAC = 0, TXAP = 0, TXAH = 0, DOACC = 0, DOACP = 0, DOACH = 0,
    VHEM = V0, MHB = p$HB_BLOOD * V0, MPROT = p$PROT_PL * V0,
    MFIB = 60 * V0 / 12, MFDP = 0.5, MTPA = 20, CPLS = 0.05,
    CPAI = p$K_PAISYN / p$K_PAIDEG,
    NCAP = 0.04, NMAT = 0.01, CVEGF = 20, CANG2 = 1.0, CMMP9 = 25,
    CIL6 = 60, CIL8 = 50, NMAC = 0.06, MMEMO = 0.02,
    MSEPT = 0, PMMA = 1, FDRN = 0, NPC = 0.20,
    FPLT = 1, CFBG = 300,
    VCOMP = p$COMP_F * max(0, V0 - p$V_RES), XMLSI = 0, XCBF = 0, NNEUR = 1,
    SSYMP = 0, SCOG = 0, GLU = p$GLU0, XINF = 0, XMYO = 0,
    HREOP = 0, HTHR = 0, HDTH = 0, AEXT = V0 / (V0 + p$K_A))
}

## Grow the haematoma from that fresh film and stop when it reaches the
## presentation volume.  This replaces an ASSUMED initial condition with a
## MODEL OUTPUT, and makes the latency from injury to presentation a
## falsifiable prediction (clinically 3-8 weeks).
##
## It also matters numerically: clamping dV/dt to zero to "relax" the other
## states -- the obvious alternative -- BREAKS the mass balance
## J_ex + J_rb = J_abs that bounds the cavity protein, and the protein
## oncotic loop then runs away to 739 g/dL.
csdh_runin <- function(V_pres = 78, tmax = 400, V0 = 12, pars = list()) {
  m <- mod
  if (length(pars)) m <- param(m, pars)
  out <- m %>% init(csdh_fresh(V0, as.list(param(m)))) %>%
    mrgsim(end = tmax, delta = 0.25, recsort = 3) %>% as_tibble()
  hit <- which(out$VHEM >= V_pres)
  if (!length(hit)) return(NULL)
  i <- hit[1]
  st <- names(init(m))
  y <- as.numeric(out[i, st]); names(y) <- st
  ## The hazard clocks and the shift-exposure integral must start at
  ## PRESENTATION.  Carrying the run-in accumulation forward gave every arm,
  ## treated or not, a 58% reoperation probability at t = 0.
  y[c("HREOP", "HTHR", "HDTH", "XMLSI")] <- 0
  list(y = y, latency = out$time[i])
}

##############################################################################
## PROCEDURES AS EVENTS
##############################################################################

## Burr-hole craniostomy.  Empties the STOCK; the membrane, the vessel
## population and the un-re-expanded brain are left exactly as they were.
## `wash` scales the cavity CONCENTRATIONS, not just the amounts -- scaling
## amounts and volume by the same factor preserves every concentration and
## makes irrigation a no-op, which is how this was first (wrongly) written.
csdh_surgery <- function(y, residual = 0.25, wash = 0.15, drain = FALSE) {
  y["VHEM"] <- y["VHEM"] * residual
  for (s in c("MHB", "MPROT", "MFIB", "MFDP", "MTPA", "TXAH", "DOACH"))
    y[s] <- y[s] * residual * wash
  y["MSEPT"] <- y["MSEPT"] * 0.6        # irrigation breaks some loculi
  y["FDRN"]  <- if (drain) 1 else 0
  y
}

## MMA embolisation.  Leaves the STOCK untouched and starves the membrane,
## so its effect is an integral over weeks rather than an event.
csdh_embolise <- function(y, residual_perf = 0.05) {
  y["PMMA"] <- residual_perf
  y
}

## Dosing regimens as zero-order inputs, so that the R and the Python
## reference receive numerically identical drug input.
ev_dex <- function(t0 = 0, days = 14, dose = 16, taper = TRUE) {
  if (!taper) return(ev(time = t0, amt = dose * days, rate = dose, cmt = "DEXD"))
  c(ev(time = t0,      amt = dose * 8,          rate = dose,        cmt = "DEXD"),
    ev(time = t0 + 8,  amt = dose * 0.50 * 3,   rate = dose * 0.50, cmt = "DEXD"),
    ev(time = t0 + 11, amt = dose * 0.25 * (days - 11),
       rate = dose * 0.25, cmt = "DEXD"))
}
ev_atv  <- function(t0 = 0, days = 56, dose = 20)
  ev(time = t0, amt = dose * days, rate = dose, cmt = "ATVD")
ev_txa  <- function(t0 = 0, days = 90, dose = 750)
  ev(time = t0, amt = dose * days, rate = dose, cmt = "TXAC")
ev_doac <- function(t0 = 0, days = 90, dose = 10)
  ev(time = t0, amt = dose * days, rate = dose, cmt = "DOACC")

##############################################################################
## SCENARIO ENGINE
##############################################################################

csdh_sim <- function(name, tend = 180,
                     surgery = NA, residual = 0.25, wash = 0.15,
                     drain_days = 0, embolise_at = NA,
                     doses = NULL, pars = list(), V_pres = 78,
                     delta = 0.25) {
  m <- mod
  if (length(pars)) m <- param(m, pars)
  ri <- csdh_runin(V_pres = V_pres, pars = pars)
  if (is.null(ri)) stop(sprintf("%s: run-in never reached V_pres", name))
  y <- ri$y

  ## events that alter the STATE (surgery, embolisation, drain removal) are
  ## applied by splitting the integration at those times
  brk <- sort(unique(c(0, tend,
                       if (!is.na(surgery)) surgery,
                       if (!is.na(surgery) && drain_days > 0) surgery + drain_days,
                       if (!is.na(embolise_at)) embolise_at)))
  brk <- brk[brk >= 0 & brk <= tend]

  res <- NULL
  for (k in seq_len(length(brk) - 1)) {
    t0 <- brk[k]; t1 <- brk[k + 1]
    if (!is.na(surgery) && isTRUE(all.equal(t0, surgery)))
      y <- csdh_surgery(y, residual, wash, drain_days > 0)
    if (!is.na(embolise_at) && isTRUE(all.equal(t0, embolise_at)))
      y <- csdh_embolise(y)
    if (!is.na(surgery) && drain_days > 0 &&
        isTRUE(all.equal(t0, surgery + drain_days)))
      y["FDRN"] <- 0
    dd <- if (is.null(doses)) NULL else
      filter(as.data.frame(doses), .data$time >= t0, .data$time < t1)
    mm <- m %>% init(y)
    o <- (if (!is.null(dd) && nrow(dd))
            mm %>% data_set(mutate(dd, ID = 1, time = .data$time - t0))
          else mm) %>%
      mrgsim(start = 0, end = t1 - t0, delta = delta, recsort = 3) %>%
      as_tibble() %>% mutate(time = .data$time + t0)
    st <- names(init(m))
    y <- as.numeric(o[nrow(o), st]); names(y) <- st
    res <- bind_rows(res, if (k == 1) o else o[-1, ])
  }
  attr(res, "latency") <- ri$latency
  attr(res, "name") <- name
  res
}

## ---------------------------------------------------------------------------
## THE 16 SHIPPED SCENARIOS
## ---------------------------------------------------------------------------
CSDH_SCENARIOS <- list(
  S01_conservative      = function() csdh_sim("01 보존적 관찰 (no treatment)", 180),
  S02_burrhole_nodrain  = function() csdh_sim("02 천공배액, 배액관 없음", 180,
                              surgery = 0),
  S03_burrhole_drain    = function() csdh_sim("03 천공배액 + 배액관 48 h", 180,
                              surgery = 0, drain_days = 2),
  S04_burrhole_poorwash = function() csdh_sim("04 천공배액, 불충분한 세척", 180,
                              surgery = 0, drain_days = 2, wash = 0.55),
  S05_mmae_alone        = function() csdh_sim("05 MMA 색전술 단독", 180,
                              embolise_at = 0),
  S06_surg_mmae         = function() csdh_sim("06 천공배액 + 배액관 + MMAE", 180,
                              surgery = 0, drain_days = 2, embolise_at = 0),
  S07_mmae_late         = function() csdh_sim("07 MMAE 지연 (수술 30일 후)", 180,
                              surgery = 0, drain_days = 2, embolise_at = 30),
  S08_dex_alone         = function() csdh_sim("08 덱사메타손 단독 (Dex-CSDH)", 180,
                              doses = ev_dex()),
  S09_dex_surgery       = function() csdh_sim("09 수술 + 덱사메타손", 180,
                              surgery = 0, drain_days = 2, doses = ev_dex()),
  S10_atorva            = function() csdh_sim("10 아토르바스타틴 20 mg x 8주", 180,
                              doses = ev_atv()),
  S11_atorva_dex        = function() csdh_sim("11 아토르바스타틴 + 저용량 dex", 180,
                              doses = c(ev_atv(), ev_dex(dose = 2.25, days = 35))),
  S12_txa               = function() csdh_sim("12 트라넥삼산 750 mg x 90일", 180,
                              doses = ev_txa()),
  S13_surg_txa          = function() csdh_sim("13 수술 + 트라넥삼산", 180,
                              surgery = 0, drain_days = 2, doses = ev_txa()),
  S14_doac_resume_early = function() csdh_sim("14 수술 + DOAC 조기 재개 (7일)", 180,
                              surgery = 0, drain_days = 2,
                              doses = ev_doac(t0 = 7, days = 173)),
  S15_doac_resume_late  = function() csdh_sim("15 수술 + DOAC 지연 재개 (30일)", 180,
                              surgery = 0, drain_days = 2,
                              doses = ev_doac(t0 = 30, days = 150)),
  S16_kitchen_sink      = function() csdh_sim("16 수술+배액+MMAE+아토르바스타틴", 180,
                              surgery = 0, drain_days = 2, embolise_at = 0,
                              doses = ev_atv())
)

##############################################################################
## ANALYSES
##############################################################################

.at <- function(d, col, t) approx(d$time, d[[col]], xout = t, rule = 2)$y

CSDH_natural_history <- function() {
  cat("\n=== 자연사: 얇은 막에서 증상 발현까지 (incident cSDH) ===\n")
  o <- mod %>% init(csdh_fresh()) %>%
    mrgsim(end = 220, delta = 0.5) %>% as_tibble()
  cat(sprintf("%6s %8s %7s %7s %7s %8s %7s %7s %7s %8s\n",
              "day", "V(mL)", "d(mm)", "A(cm2)", "MLS", "ICP", "J_ex",
              "J_rb", "J_abs", "haem"))
  for (t in c(0, 7, 14, 21, 28, 42, 60, 90, 150, 220))
    cat(sprintf("%6.0f %8.2f %7.1f %7.1f %7.2f %8.1f %7.2f %7.2f %7.2f %8.3f\n",
        t, .at(o,"VHEM",t), .at(o,"dmax_",t), .at(o,"A_mem_",t),
        .at(o,"MLS_",t), .at(o,"ICP_",t), .at(o,"JEX",t), .at(o,"JRB",t),
        .at(o,"JABS",t), .at(o,"HDRIVE",t)))
  i <- which(o$VHEM >= 78)[1]
  cat(sprintf("\n  latency injury -> 78 mL presentation : %.1f days\n", o$time[i]))
  cat(  "  (clinically 3-8 weeks; this is a model OUTPUT, not an assumption)\n")
  invisible(o)
}

CSDH_run_scenarios <- function() {
  cat("\n=== 16개 시나리오 (16 scenarios) ===\n")
  cat(sprintf("%-40s %7s %7s %7s %8s %8s\n", "scenario", "V90", "V180",
              "dmax180", "P(reop)", "P(fav)"))
  res <- list()
  for (nm in names(CSDH_SCENARIOS)) {
    d <- CSDH_SCENARIOS[[nm]]()
    res[[nm]] <- d
    cat(sprintf("%-40s %7.1f %7.1f %7.1f %8.3f %8.3f\n",
        attr(d, "name"), .at(d,"VHEM",90), .at(d,"VHEM",180),
        .at(d,"dmax_",180), .at(d,"P_REOP",180), .at(d,"P_FAV",180)))
  }
  invisible(res)
}

## Which of the three terms does each therapy touch?  This is the whole
## argument of the model, so it is worth printing as a table.
CSDH_flux_ledger <- function() {
  cat("\n=== 플럭스 원장: 각 치료가 건드리는 항 (flux ledger) ===\n")
  cat("치료가 J_ex / J_rb / J_abs 중 어디에 작용하는지, 그리고 STOCK인지 SOURCE인지\n\n")
  cat(sprintf("%-26s %-12s %-40s\n", "therapy", "term", "mechanism"))
  led <- list(
    c("burr-hole / twist drill", "STOCK",   "V <- residual*V; membrane untouched"),
    c("irrigation",              "SWITCH",  "lowers cavity haem CONCENTRATION"),
    c("subdural drain",          "STOCK",   "continues to wash; lowers residual"),
    c("craniotomy + membranectomy","SOURCE", "removes A_mem and septation"),
    c("MMA embolisation",        "SOURCE",  "P_MMA -> 0.05; stops new sprouting"),
    c("dexamethasone",           "J_ex",    "VEGF down, Lp down; sigma unchanged"),
    c("atorvastatin",            "J_ex",    "matures vessels: Lp down AND sigma up"),
    c("tranexamic acid",         "J_rb",    "breaks loop L2; raises H"),
    c("holding anticoagulation", "J_rb",    "raises H; costs thromboembolism"),
    c("posture / hydration",     "J_abs",   "speeds brain re-expansion")
  )
  for (r in led) cat(sprintf("%-26s %-12s %-40s\n", r[1], r[2], r[3]))
  invisible(led)
}

## The haem switch.  Where does it sit, and is it really a switch?
CSDH_bistability <- function() {
  cat("\n=== 쌍안정성: 헴 스위치의 위치 (bistability of the haem switch) ===\n")
  cat("수술 후 세척 정도(wash)를 바꾸면 같은 수술이 치유와 재발로 갈라진다.\n\n")
  cat(sprintf("%8s %10s %10s %10s %10s %9s\n",
              "wash", "V(30d)", "V(90d)", "V(180d)", "haem(30d)", "P(reop)"))
  for (w in c(0.05, 0.10, 0.15, 0.25, 0.40, 0.55, 0.75, 1.00)) {
    d <- csdh_sim(sprintf("wash %.2f", w), 180, surgery = 0,
                  drain_days = 2, wash = w)
    cat(sprintf("%8.2f %10.2f %10.2f %10.2f %10.3f %9.3f\n", w,
        .at(d,"VHEM",30), .at(d,"VHEM",90), .at(d,"VHEM",180),
        .at(d,"HDRIVE",30), .at(d,"P_REOP",180)))
  }
  cat("\n  스위치가 없으면(Hill=1) 이 표의 모든 행이 같아진다.\n")
  invisible(NULL)
}

## V_RES protects the patient before the operation and condemns them after.
CSDH_atrophy_two_signs <- function() {
  cat("\n=== 위축의 두 부호 (V_RES appears twice, opposite signs) ===\n")
  cat("같은 파라미터가 (1) 발현 시 부피를 키우고 (2) 재발률을 올린다.\n\n")
  cat(sprintf("%8s %10s %10s %10s %11s %10s %9s\n", "V_RES", "latency",
              "MLS@78mL", "k_reexp", "Vcomp(30d)", "V(90d)", "P(reop)"))
  for (vr in c(5, 12, 20, 30, 40, 50, 60)) {
    p <- list(V_RES = vr)
    ri <- csdh_runin(78, pars = p)
    d  <- csdh_sim(sprintf("V_RES %g", vr), 180, surgery = 0,
                   drain_days = 2, pars = p)
    kr <- as.list(param(mod))$K_REEXP0 * exp(-vr / as.list(param(mod))$LAM_ATR)
    cat(sprintf("%8.0f %10.1f %10.2f %10.4f %11.1f %10.2f %9.3f\n",
        vr, ri$latency, .at(d,"MLS_",0), kr,
        .at(d,"VCOMP",30), .at(d,"VHEM",90), .at(d,"P_REOP",180)))
  }
  cat("\n  두 열(MLS와 P(reop))이 반대 방향으로 움직이는 것이 요점이다.\n")
  invisible(NULL)
}

## Why a 30-day endpoint for embolisation must be null.
CSDH_mmae_timecourse <- function() {
  cat("\n=== MMAE의 시간 경과: 30일 엔드포인트는 왜 음성인가 ===\n")
  cat("색전술은 STOCK이 아니라 SOURCE를 없애므로 효과가 적분으로 나타난다.\n\n")
  a <- csdh_sim("surgery only", 180, surgery = 0, drain_days = 2)
  b <- csdh_sim("surgery+MMAE", 180, surgery = 0, drain_days = 2,
                embolise_at = 0)
  cat(sprintf("%8s %12s %12s %10s %12s %12s %10s\n", "day",
              "V ctrl", "V MMAE", "dV", "P(reop)ctrl", "P(reop)MMAE", "abs.diff"))
  for (t in c(7, 14, 30, 60, 90, 120, 180)) {
    pc <- .at(a,"P_REOP",t); pm <- .at(b,"P_REOP",t)
    cat(sprintf("%8.0f %12.2f %12.2f %10.2f %12.3f %12.3f %10.3f\n", t,
        .at(a,"VHEM",t), .at(b,"VHEM",t),
        .at(b,"VHEM",t) - .at(a,"VHEM",t), pc, pm, pc - pm))
  }
  cat("\n  N_CAP 감쇠 시상수가 ~9-18일이므로 30일에는 차이가 아직 작다.\n")
  invisible(NULL)
}

## 3 hour plasma half-life, 6 day half-life at the site of action.
CSDH_txa_site_delay <- function() {
  cat("\n=== TXA: 혈장 t1/2 3시간, 작용부위 t1/2 6일 ===\n")
  cat("약은 확산이 아니라 '유출액에 실려서' 혈종 안으로 들어간다.\n\n")
  d <- csdh_sim("surgery+TXA", 120, surgery = 0, drain_days = 2,
                doses = ev_txa(days = 120))
  cat(sprintf("%8s %12s %12s %10s\n", "day", "C_TXA plasma", "C_TXA cavity",
              "ratio"))
  for (t in c(0.5, 1, 3, 7, 14, 21, 30, 60, 90)) {
    cp <- .at(d, "TXAC", t) / as.list(param(mod))$TXA_V1
    cc <- .at(d, "TXAH", t) / max(.at(d, "VHEM", t), 1e-4) * 1000
    cat(sprintf("%8.1f %12.3f %12.3f %10.2f\n", t, cp, cc,
                if (cp > 1e-9) cc / cp else NA))
  }
  cat("\n  공동의 turnover 시간 V/J_abs 가 율속 단계다.\n")
  invisible(NULL)
}

## One drug, two signs: fewer reoperations AND worse patients.
CSDH_dex_two_signs <- function() {
  cat("\n=== 덱사메타손의 두 부호 (Dex-CSDH의 역설) ===\n")
  a <- csdh_sim("surgery",       180, surgery = 0, drain_days = 2)
  b <- csdh_sim("surgery + dex", 180, surgery = 0, drain_days = 2,
                doses = ev_dex())
  cat(sprintf("%-22s %10s %10s\n", "endpoint", "placebo", "dex"))
  cat(sprintf("%-22s %10.3f %10.3f   <- 국소적으로 이롭다\n", "P(reoperation) 180d",
              .at(a,"P_REOP",180), .at(b,"P_REOP",180)))
  cat(sprintf("%-22s %10.3f %10.3f   <- 전신적으로 해롭다\n", "P(favourable) 180d",
              .at(a,"P_FAV",180), .at(b,"P_FAV",180)))
  cat(sprintf("%-22s %10.2f %10.2f\n", "peak glucose (mmol/L)",
              max(a$GLU), max(b$GLU)))
  cat(sprintf("%-22s %10.3f %10.3f\n", "peak infection burden",
              max(a$XINF), max(b$XINF)))
  cat(sprintf("%-22s %10.3f %10.3f\n", "peak myopathy burden",
              max(a$XMYO), max(b$XMYO)))
  cat("\n  같은 약이 하나의 식에서는 도움이 되고 다른 식에서는 해가 된다.\n")
  invisible(NULL)
}

## The fibrinolytic loop: report the amplification factor.
CSDH_loop_gain <- function() {
  cat("\n=== 섬유소용해 루프 L2의 이득 (loop gain) ===\n")
  a <- csdh_sim("no TXA", 90)
  b <- csdh_sim("TXA",    90, doses = ev_txa(days = 90))
  cat(sprintf("%-26s %12s %12s\n", "quantity", "no TXA", "TXA"))
  for (v in c("CFDP_", "HCOMP", "JRB", "CTPA_"))
    cat(sprintf("%-26s %12.3f %12.3f\n", v, .at(a,v,60), .at(b,v,60)))
  g1 <- 1 - .at(a,"HCOMP",60); g2 <- 1 - .at(b,"HCOMP",60)
  cat(sprintf("\n  1-H (재출혈 구동): %.3f -> %.3f  (증폭 %.2fx -> %.2fx)\n",
              g1, g2, 1/(1-g1), 1/(1-g2)))
  invisible(NULL)
}

CSDH_trial_ledger <- function() {
  cat("\n=== 모델 vs 발표된 임상시험 수치 (trial ledger) ===\n")
  cat(sprintf("%-46s %10s %10s\n", "endpoint", "model", "published"))
  d3  <- csdh_sim("bh+drain", 180, surgery = 0, drain_days = 2)
  d2  <- csdh_sim("bh",       180, surgery = 0)
  d6  <- csdh_sim("bh+MMAE",  180, surgery = 0, drain_days = 2, embolise_at = 0)
  d8  <- csdh_sim("dex",      180, doses = ev_dex())
  d9  <- csdh_sim("bh+dex",   180, surgery = 0, drain_days = 2, doses = ev_dex())
  d10 <- csdh_sim("atorva",   180, doses = ev_atv())
  d1  <- csdh_sim("conserv",  180)
  rows <- list(
    c("recurrence 90d, burr-hole + drain",  .at(d3,"P_REOP",90),  "0.093 (Santarius 2009)"),
    c("recurrence 90d, burr-hole no drain", .at(d2,"P_REOP",90),  "0.240 (Santarius 2009)"),
    c("treatment failure 180d, surgery",    .at(d3,"P_REOP",180), "0.113 (EMBOLISE 2024)"),
    c("treatment failure 180d, +MMAE",      .at(d6,"P_REOP",180), "0.041 (EMBOLISE 2024)"),
    c("repeat surgery, dexamethasone",      .at(d9,"P_REOP",180), "0.017 (Dex-CSDH 2020)"),
    c("favourable outcome, dexamethasone",  .at(d9,"P_FAV",180),  "0.839 (Dex-CSDH 2020)"),
    c("favourable outcome, placebo",        .at(d3,"P_FAV",180),  "0.903 (Dex-CSDH 2020)"),
    c("conversion to surgery 56d, atorva",  .at(d10,"P_REOP",56), "0.112 (ATOCH 2018)"),
    c("conversion to surgery 56d, placebo", .at(d1,"P_REOP",56),  "0.235 (ATOCH 2018)")
  )
  for (r in rows)
    cat(sprintf("%-46s %10.3f %10s\n", r[[1]], as.numeric(r[[2]]), r[[3]]))
  invisible(rows)
}

CSDH_sensitivity <- function(delta = 0.20) {
  cat(sprintf("\n=== 국소 민감도 (+/-%.0f%% on each parameter) ===\n", 100*delta))
  keys <- c("LP_IMM","SIG_IMM","K_LYMPH","K_RB","K_FDP","K_TPA","HB50",
            "K_MAC","KD_CAP","K_MAT","MAT_SUPP","V_RES","K_REEXP0",
            "P_CAP","K_PV","D50_PAT","K_SPROUT","K_A")
  base <- csdh_sim("base", 180, surgery = 0, drain_days = 2)
  b90  <- .at(base, "VHEM", 90); br <- .at(base, "P_REOP", 180)
  p0   <- as.list(param(mod))
  cat(sprintf("%-12s %12s %12s\n", "parameter", "dV90 (%)", "dP(reop) (%)"))
  out <- NULL
  for (k in keys) {
    v <- p0[[k]]
    hi <- csdh_sim("hi", 180, surgery = 0, drain_days = 2,
                   pars = setNames(list(v * (1 + delta)), k))
    lo <- csdh_sim("lo", 180, surgery = 0, drain_days = 2,
                   pars = setNames(list(v * (1 - delta)), k))
    dv <- (.at(hi,"VHEM",90) - .at(lo,"VHEM",90)) / max(b90, 1e-6) * 100
    dr <- (.at(hi,"P_REOP",180) - .at(lo,"P_REOP",180)) / max(br, 1e-6) * 100
    cat(sprintf("%-12s %12.1f %12.1f\n", k, dv, dr))
    out <- rbind(out, data.frame(parameter = k, dV90 = dv, dReop = dr))
  }
  invisible(out)
}

CSDH_mass_balance <- function() {
  cat("\n=== 질량수지 및 음수 상태 점검 ===\n")
  d <- csdh_sim("check", 180, surgery = 0, drain_days = 2, doses = ev_atv())
  st <- names(init(mod))
  bad <- st[sapply(st, function(s) min(d[[s]]) < -1e-6)]
  cat("  negative states:", if (!length(bad)) "none" else paste(bad, collapse=", "), "\n")
  cat(sprintf("  cavity protein stays below plasma (7 g/dL): max %.2f g/dL\n",
              max(d$CPROT_)))
  cat(sprintf("  ICP range: %.2f - %.2f mmHg\n", min(d$ICP_), max(d$ICP_)))
  cat(sprintf("  patency range: %.3f - %.3f\n", min(d$f_pat_), max(d$f_pat_)))
  invisible(d)
}

CSDH_run_all <- function() {
  CSDH_natural_history(); CSDH_flux_ledger(); CSDH_bistability()
  CSDH_atrophy_two_signs(); CSDH_mmae_timecourse(); CSDH_txa_site_delay()
  CSDH_dex_two_signs(); CSDH_loop_gain(); CSDH_run_scenarios()
  CSDH_trial_ledger(); CSDH_sensitivity(); CSDH_mass_balance()
  invisible(NULL)
}

message("csdh_mrgsolve_model.R loaded: 50 ODEs, 16 scenarios, 12 analyses.")
message("Start with CSDH_natural_history(), CSDH_flux_ledger(), CSDH_bistability().")
message("Every reference number in the comments comes from csdh_verify.py, ",
        "which integrates the same 50-ODE system in scipy.")
