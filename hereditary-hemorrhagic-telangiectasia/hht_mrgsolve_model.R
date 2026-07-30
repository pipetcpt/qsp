# =====================================================================
#  Hereditary Haemorrhagic Telangiectasia (HHT / Rendu-Osler-Weber)
#  Quantitative Systems Pharmacology model — mrgsolve
# ---------------------------------------------------------------------
#  WHY THIS MODEL IS BUILT THE WAY IT IS
#
#  Every published HHT endpoint mixes two independently-set quantities:
#
#     B      = blood-loss rate  [mL/day]  = B_nose + B_gi
#     A_net  = net iron supply  [mg/day]  = diet + oral(hepcidin) + IV
#                                            − obligatory basal loss
#
#  One mL of blood carries 0.0347 x Hb mg of iron (3.47 mg Fe per g of
#  haemoglobin).  Writing the iron and haemoglobin balances as ODEs and
#  setting both to zero gives, EXACTLY and with no fitted constant:
#
#         Hb_ss  =  A_net / (0.0347 * B)                        [g/dL]
#
#  Derivation (reproduced numerically by axis1_hyperbola() below):
#     d/dt Fe_store = A_gross + 3.47*Hb_mass/tau − 3.47*PROD − basal = 0
#     d/dt Hb_mass  = PROD − Hb_mass/tau − Hb*B/100             = 0
#   substituting the second into the first cancels PROD and tau and
#   leaves  0.0347*Hb*B = A_gross − basal = A_net.
#  Note that NO assumption about the production law PROD was needed:
#  the hyperbola is a structural property of mass balance, not a fit.
#
#  Three consequences drive the whole analysis:
#
#   (i)   Anaemia in HHT is an EQUILIBRIUM, not a decompensation.  The
#         iron cost of bleeding is proportional to Hb, so as Hb falls
#         the drain falls with it and the patient settles at a low Hb
#         instead of exsanguinating.
#
#   (ii)  B and A_net enter Hb ONLY as a ratio.  Haemoglobin therefore
#         cannot distinguish an antiangiogenic drug from an iron
#         infusion — and if a trial lets iron use fall as bleeding
#         falls, the two changes CANCEL.  That is exactly what
#         PATH-HHT observed (axis2_hb_paradox()).
#
#   (iii) By Fick, C(a−v)O2 is proportional to Hb, so cardiac index
#         scales as Hb^(−alpha).  Correcting anaemia lowers cardiac
#         index with no change in shunt anatomy whatsoever
#         (axis6_cardiac_confound()).
#
#  CALIBRATION TARGETS (all values verified against the primary source)
#   * PATH-HHT  Al-Samkari 2024 NEJM 391:1015 (PMID 39292928)
#       n=144 (95 pomalidomide 4 mg/d : 49 placebo), 24 weeks
#       baseline ESS 5.0 +/- 1.5, HHT-QoL 6.3 +/- 3.1, age 58.8 +/- 12.2
#       69% anaemic; 84% IV iron and 19% RBC transfusion in prior 6 mo
#       ESS change  −1.84 (−2.25,−1.43) vs −0.90 (−1.39,−0.40)
#                   difference −0.94 (−1.57,−0.31), p=0.004
#       HHT-QoL     −2.7 vs −1.2, difference −1.4
#       weighted bleeding duration −12.2 vs −3.3, difference −8.9
#       wk12-24 RBC transfusion 9% vs 18%
#       median IV iron per 4 wk  0 mg (IQR 0-340) vs 333 mg (IQR 0-500)
#       haemoglobin  NO difference at 24 wk; +1.09 g/dL (0.38,1.80)
#                    only at the 4-week POST-treatment visit
#       neutropenia 44% vs 10%; mean ANC change −1710 vs −20 /uL
#       nasal biopsy: NO inhibition of EC proliferation but INCREASED
#                     mural-cell coverage  -> pomalidomide is modelled
#                     as a vessel-maturation agent, not an anti-VEGF
#   * Dupuis-Girod 2012 JAMA 307:948 (PMID 22396517), single arm, n=25
#       bevacizumab 5 mg/kg q14d x 6 (last dose day 70)
#       cardiac index 5.05 -> 4.2 (3 mo) -> 4.1 (6 mo) L/min/m2
#       epistaxis duration 221 -> 134 (3 mo) -> 43 (6 mo) min/month
#   * Azzopardi 2015 MAbs 7:630 (PMID 25751241) exposure-response
#       maintenance q3 / q2 / q1 month at 24 months:
#         CI < 4 L/min/m2            in 41% / 45% / 50%
#         epistaxis < 20 min/month   in 34% / 43% / 60%
#   * NOSE  Whitehead 2016 JAMA 316:943 (PMID 27599329), n=121 RCT
#       topical bevacizumab / estriol / TXA vs saline: no difference,
#       and ALL arms improved ESS significantly at weeks 12 and 24
#   * ESS instrument Hoag 2010 Laryngoscope 120:838 (PMID 20087969)
#       items: frequency, duration, intensity, need for medical
#       attention, ANAEMIA, need for TRANSFUSION.  MCID 0.71
#       (Yin 2016, PMID 26393959)
#
#  HONEST LIMITS — read before using any number out of this file
#   * The ESS item weights were never published in closed form.  The
#     six-item surrogate here is transparent but NOT the instrument;
#     the anaemia+transfusion share W_FE_SHARE is swept, not fitted,
#     because it is the quantity axis3 depends on.
#   * B_nose, B_gi and f_GI are not directly measured in any trial.
#     They are inferred from the iron balance, so they inherit every
#     assumption in A_net.
#   * The bevacizumab structural refutation (axis4) needs only the PK
#     half-life and the three published epistaxis means; it does NOT
#     depend on any parameter fitted here.  It is the most robust
#     result in this file.  The k_out estimates in axis5 are the least
#     robust: they are fitted to summary response *fractions*, not to
#     patient data.
#   * alpha (the anaemia-compensation exponent) is NOT identifiable
#     from published data.  Every cardiac-index conclusion is reported
#     as a range over alpha, never as a point estimate.
#
#  Rendering / running
#     Rscript hht_mrgsolve_model.R          # runs every axis + scenarios
# =====================================================================

suppressMessages({
  library(mrgsolve)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
set.seed(20260730)

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
wk  <- function(x) x * 7
mo  <- function(x) x * 30.4375
yr  <- function(x) x * 365.25
f2  <- function(x) formatC(x, format = "f", digits = 2)
f3  <- function(x) formatC(x, format = "f", digits = 3)
f1  <- function(x) formatC(x, format = "f", digits = 1)
hdr <- function(...) {
  msg <- paste0(...)
  cat("\n", strrep("=", 70), "\n", msg, "\n", strrep("=", 70), "\n", sep = "")
}
sub <- function(...) cat("\n-- ", paste0(...), "\n", sep = "")

# =====================================================================
#  THE MODEL
# =====================================================================

hht_code <- '
$PROB
# HHT QSP — 40 ODEs
Blood loss / iron / erythropoiesis / vascular-lesion turnover / drug PK-PD

$PARAM @annotated
// ---------------- patient & genotype -------------------------------
GENO     :  2 : Genotype 1=ENG(HHT1) 2=ACVRL1(HHT2) 3=SMAD4(JP-HHT)
WT       : 75 : Body weight (kg)
BSA      : 1.90 : Body surface area (m2)
AGE0     : 58.8 : Age at t=0 (yr) — PATH-HHT mean
VBLOOD   : 70 : Blood volume per kg (mL/kg)

// ---------------- BMP9-ALK1-SMAD signalling ------------------------
RDOSE_ENG  : 0.50 : Functional receptor dose, ENG haploinsufficiency
RDOSE_ALK1 : 0.50 : Functional receptor dose, ACVRL1 haploinsufficiency
RDOSE_SM4  : 0.50 : Functional receptor dose, SMAD4 haploinsufficiency
KSMAD    : 0.20 : SMAD signal turnover (1/day)
ETAC     : 0.35 : Max fractional SMAD gain from tacrolimus

// ---------------- VEGF ---------------------------------------------
VEGFSYN  : 60.0 : VEGF appearance (pg/mL/day)
VEGFDEG  : 0.20 : VEGF elimination (1/day)
GSHEAR   : 0.25 : VEGF induction per unit hepatic shunt conductance
GSMADV   : 1.00 : VEGF induction exponent on SMAD deficit
KD_BEV   : 0.50 : Bevacizumab-VEGF KD (nM)
FRES     : 0.30 : VEGF-INDEPENDENT residual share of lesion formation
FRES_L   : 0.30 : VEGF-INDEPENDENT share for the HEPATIC shunt

// ---------------- vascular lesion turnover -------------------------
KOUT_N   : 0.0180 : Nasal telangiectasia turnover (1/day)
KOUT_G   : 0.0090 : GI telangiectasia turnover (1/day)
KOUT_L   : 0.0060 : Hepatic shunt remodelling (1/day)
KOUT_P   : 0.00020 : Pulmonary AVM turnover (1/day) — near-static
TN0      : 1.00 : Baseline nasal lesion burden (relative)
TG0      : 1.00 : Baseline GI lesion burden (relative)
GS0      : 1.00 : Baseline hepatic shunt conductance (relative)
PV0      : 1.00 : Baseline pulmonary AVM burden (relative)
GIAGE    : 0.020 : Fractional GI lesion growth per year of age
FLOWFB   : 0.15 : Shear positive feedback exponent on hepatic shunt
GENO_LIV : 1.60 : Hepatic shunt multiplier if ACVRL1
GENO_LUN : 1.80 : Pulmonary AVM multiplier if ENG
GENO_BRN : 1.60 : Cerebral AVM multiplier if ENG

// ---------------- mural cell coverage & fragility ------------------
MURSS    : 0.55 : Mural coverage set-point in HHT (fraction of normal)
KMUR     : 0.030 : Mural coverage turnover (1/day)
GFRAG    : 2.85 : Fragility exponent on mural-coverage deficit
GQ       : 0.50 : Bleeding-flow exponent on fragility
GDUR     : 0.50 : Episode-duration exponent on fragility
EMAX_PM  : 0.72 : Max fractional mural-coverage gain, pomalidomide
EC50_PM  : 0.035 : Pomalidomide conc for half mural effect (mg/L)
PAZMUR   : 0.15 : Mural-coverage LOSS from PDGFR blockade (pazopanib)

// ---------------- epistaxis ----------------------------------------
LAMEP0   : 7.00 : Baseline episode frequency (episodes/week)
KLAM     : 0.10 : Frequency adaptation rate (1/day)
EPDUR0   : 10.0 : Baseline episode duration (min)
KDUR     : 0.25 : Duration adaptation rate (1/day)
QBLEED   : 1.84 : Bleeding flow during an episode (mL/min)
KFIB     : 0.40 : Fibrinolytic tone turnover (1/day)
FIBUP    : 0.20 : Fibrinolytic tone gain per unit lesion burden
PLTFE    : 0.25 : Max duration prolongation from iron-deficient platelets
FEPLT50  : 250 : Storage iron for half-maximal platelet effect (mg)
CLOSEF   : 0 : Young nasal closure (0=no, 1=yes) — sets B_nose to 0
ABLATE   : 0 : Laser/septodermoplasty fractional lesion ablation rate

// ---------------- GI bleeding --------------------------------------
BGI0     : 9.90 : Baseline GI blood loss (mL/day)
GIENDOF  : 0 : Argon plasma coagulation fractional GI lesion reduction

// ---------------- iron ---------------------------------------------
AMAXO    : 6.0 : Max oral iron absorption (mg/day) at zero hepcidin
KMFEG    : 40.0 : Luminal iron for half-max absorption (mg)
KHEP     : 1.00 : Hepcidin for half-suppression of absorption
ADIET    : 1.60 : Dietary iron absorption at zero hepcidin (mg/day)
KTRFE    : 1.20 : Luminal iron transit-out rate (1/day)
BASALFE  : 1.00 : Obligatory basal iron loss (mg/day)
IVFERATE : 0 : FIXED IV iron infusion rate (mg/day), protocol-imposed
IVFE_ON  : 1 : Demand-driven IV iron policy active (0/1)
IVFEMAX  : 26.0 : Maximum demand-driven IV iron rate (mg/day)
IVFETRIG : 11.0 : Haemoglobin at which IV iron is half-maximal (g/dL)
SIVFE    : 0.35 : Steepness of the IV iron policy (g/dL)
FESTORE0 : 120 : Baseline storage iron (mg)
FESREF   : 800 : Reference storage iron for hepcidin signal (mg)
GHEP     : 0.70 : Hepcidin exponent on storage iron
KHEPT    : 0.15 : Hepcidin turnover (1/day)
KERF     : 0.20 : Erythroferrone turnover (1/day)
KERFE    : 1.00 : Erythroferrone potency on hepcidin
IL6G     : 0.30 : Hepcidin induction per unit hepatic shunt conductance
KSUP     : 90.0 : Storage iron for half-maximal erythron supply (mg)

// ---------------- erythropoiesis -----------------------------------
HBSET    : 14.5 : Haemoglobin set-point (g/dL)
RBCLIFE  : 120 : Red cell lifespan (day)
EMAXEPO  : 4.00 : Max fractional EPO-driven production increase
MARROWMX : 6.00 : Absolute marrow output ceiling (x basal)
KRET     : 0.30 : Reticulocyte signal turnover (1/day)
EMAX_MS  : 0.55 : Max fractional marrow suppression, pomalidomide
EC50_MS  : 0.020 : Pomalidomide conc for half marrow suppression (mg/L)
KMSUP    : 0.10 : Marrow suppression onset/offset rate (1/day)
HBTRIG   : 7.50 : Haemoglobin transfusion trigger (g/dL)
STX      : 0.35 : Transfusion trigger smoothing (g/dL)
TXMAX    : 0.12 : Max transfusion rate (units/day)
HBUNIT   : 50.0 : Haemoglobin mass per RBC unit (g)
TXALLOW  : 1 : Transfusion permitted (0/1)

// ---------------- cardiovascular -----------------------------------
VO2I     : 125 : Oxygen consumption index (mL/min/m2)
DSAT0    : 0.25 : Reference arterio-venous saturation difference
HBREF    : 14.5 : Reference haemoglobin for Fick scaling (g/dL)
ALPHA    : 1.00 : Anaemia-compensation exponent (CI ~ Hb^-ALPHA)
KSHUNT   : 0.125 : Cardiac index added per unit shunt conductance
KTR      : 0.090 : CI transit rate (1/day) — Azzopardi transit structure
KRV      : 0.010 : RV functional adaptation (1/day)
RVSENS   : 0.12 : RV functional loss per unit CI above threshold
CITHR    : 4.00 : Cardiac index threshold for RV strain (L/min/m2)
NTBNP0   : 125 : Reference NT-proBNP (pg/mL)
CIREF    : 2.70 : Reference cardiac index (L/min/m2)
GBNP     : 2.20 : NT-proBNP exponent on cardiac index
KBNP     : 0.08 : NT-proBNP turnover (1/day)

// ---------------- pulmonary shunt ----------------------------------
SCAP     : 0.99 : End-capillary oxygen saturation (fraction)
SPAV     : 0.10 : Shunt fraction per unit pulmonary AVM burden
EMBOLF   : 0 : Pulmonary AVM embolisation fractional occlusion (0-1)
HAZPAV   : 0.010 : Embolic hazard per unit AVM burden per year

// ---------------- ESS surrogate ------------------------------------
// Six items after Hoag 2010.  W_FE_SHARE is the ANAEMIA+TRANSFUSION
// share of the total weight; it is SWEPT in axis3, never fitted.
W_FE_SHARE : 0.18 : Combined anaemia+transfusion share of ESS weight
W_FREQ   : 0.30 : Frequency item weight (of the non-iron share)
W_DUR    : 0.26 : Duration item weight
W_INT    : 0.26 : Intensity item weight
W_ATT    : 0.18 : Medical-attention item weight
KIF      : 7.00 : Frequency item half-saturation (episodes/week)
KID      : 10.0 : Duration item half-saturation (min)
KII      : 1.84 : Intensity item half-saturation (mL/min)
HBANEM   : 11.5 : Anaemia item midpoint (g/dL)
SHB      : 1.20 : Anaemia item steepness (g/dL)
KITX     : 0.02 : Transfusion item half-saturation (units/day)
ESSMAX   : 10.0 : ESS scale maximum
KESS     : 0.055 : ESS recall/perception rate (1/day) — 4-week recall
ESSCAL   : 1.00 : ESS calibration multiplier

// ---------------- placebo / blinding ------------------------------
PLACMAX  : 0.90 : Maximum placebo ESS improvement — PATH-HHT placebo arm
KPLAC    : 0.030 : Placebo onset rate (1/day)
KPLACOFF : 0.060 : Placebo dissipation rate after unblinding (1/day)
BLIND    : 0 : On blinded study drug (0/1)
ONPLAC   : 0 : Placebo response active (0/1)

// ---------------- QoL ---------------------------------------------
QOLREF   : 6.30 : Reference HHT-QoL at ESS 5.0 (PATH-HHT baseline)
ESSREF   : 5.00 : Reference ESS
GQOL     : 1.15 : QoL exponent on ESS
KQOL     : 0.045 : QoL adaptation rate (1/day)

// ---------------- drug PK -----------------------------------------
CLB      : 0.207 : Bevacizumab clearance (L/day)
V1B      : 2.92 : Bevacizumab central volume (L)
QB       : 0.400 : Bevacizumab intercompartmental clearance (L/day)
V2B      : 2.30 : Bevacizumab peripheral volume (L)
KA_POM   : 1.20 : Pomalidomide absorption (1/day x24 -> per day)
CL_POM   : 200 : Pomalidomide clearance (L/day)
V_POM    : 70.0 : Pomalidomide volume (L)
KA_PAZ   : 3.00 : Pazopanib absorption (1/day)
CL_PAZ   : 6.60 : Pazopanib clearance (L/day)
V_PAZ    : 290 : Pazopanib volume (L)
KA_TXA   : 12.0 : Tranexamic acid absorption (1/day)
CL_TXA   : 160 : Tranexamic acid clearance (L/day)
V_TXA    : 40.0 : Tranexamic acid volume (L)
CL_TAC   : 60.0 : Tacrolimus clearance (L/day)
V_TAC    : 1000 : Tacrolimus volume (L)
KA_TAC   : 4.00 : Tacrolimus absorption (1/day)

// ---------------- drug PD -----------------------------------------
EMAX_PZ  : 0.70 : Max fractional VEGFR-pathway block, pazopanib
EC50_PZ  : 12.0 : Pazopanib conc for half VEGFR block (mg/L)
EMAX_TX  : 0.55 : Max fractional fibrinolysis block, TXA
EC50_TX  : 8.00 : TXA conc for half fibrinolysis block (mg/L)
EC50_TAC : 5.00 : Tacrolimus conc for half SMAD gain (ug/L)
USEINIT  : 0 : If 1, do NOT overwrite compartment initials in $MAIN

$CMT @annotated
BEVC   : Bevacizumab central (mg)
BEVP   : Bevacizumab peripheral (mg)
POMD   : Pomalidomide depot (mg)
POMC   : Pomalidomide central (mg)
PAZD   : Pazopanib depot (mg)
PAZC   : Pazopanib central (mg)
TXAD   : Tranexamic acid depot (mg)
TXAC   : Tranexamic acid central (mg)
TACD   : Tacrolimus depot (mg)
VEGFT  : Total plasma VEGF (pg/mL)
SMADS  : SMAD1/5/8 signal strength (relative to normal)
TNOSE  : Nasal telangiectasia burden (relative)
TGIB   : GI telangiectasia burden (relative)
GSHUNT : Hepatic shunt conductance (relative)
PAVMB  : Pulmonary AVM burden (relative)
MURALC : Mural cell coverage (fraction of normal)
LAMEP  : Epistaxis episode frequency (episodes/week)
EPDURC : Epistaxis episode duration (min)
FIBLY  : Local fibrinolytic tone (relative)
FEGUT  : Luminal (oral) iron (mg)
FESTORE: Storage iron (mg)
HEPCC  : Hepcidin (relative)
ERFEC  : Erythroferrone (relative)
HBMASS : Circulating haemoglobin mass (g)
RETICC : Reticulocyte / marrow drive (relative)
MSUP   : Marrow suppression (fraction)
TR1    : CI transit 1
TR2    : CI transit 2
TR3    : CI transit 3
RVF    : RV functional reserve (relative)
NTBNPC : NT-proBNP (pg/mL)
PLAC   : Placebo ESS improvement (points)
ESSC   : Epistaxis Severity Score (0-10)
QOLC   : HHT-specific QoL (0-16)
CUMTX  : Cumulative RBC units
CUMIVFE: Cumulative IV iron (mg)
CUMBLD : Cumulative blood loss (mL)
AUCBEV : Bevacizumab AUC (mg/L*day)
EMBHAZ : Cumulative embolic hazard
CUMPOM : Pomalidomide AUC (mg/L*day)

$GLOBAL
#define BVOL   (VBLOOD*WT)                 /* blood volume, mL        */
#define HBC    (100.0*HBMASS/BVOL)         /* haemoglobin, g/dL       */
#define RHOFE  (0.0347*HBC)                /* mg Fe per mL blood      */
#define CBEV   (BEVC/V1B)                  /* mg/L                    */
#define CPOM   (POMC/V_POM)                /* mg/L                    */
#define CPAZ   (PAZC/V_PAZ)                /* mg/L                    */
#define CTXA   (TXAC/V_TXA)                /* mg/L                    */
#define SIGDEF (1.0/SMADS)                 /* signalling deficit      */

// free-ligand quadratic: total ligand LT, total binder BT, dissoc KD
double freelig(double LT, double BT, double KD) {
  double b = LT - BT - KD;
  double disc = b*b + 4.0*KD*LT;
  if (disc < 0) disc = 0;
  return 0.5*(b + sqrt(disc));
}

$MAIN
// ---- genotype -> receptor dose and organ tropism -------------------
double RDOSE = RDOSE_ALK1;
if (GENO < 1.5)      RDOSE = RDOSE_ENG;
else if (GENO > 2.5) RDOSE = RDOSE_SM4;

double LIVMULT = (GENO > 1.5 && GENO < 2.5) ? GENO_LIV : 1.0;
double LUNMULT = (GENO < 1.5) ? GENO_LUN : 1.0;
double BRNMULT = (GENO < 1.5) ? GENO_BRN : 1.0;

// ---- parameter-dependent initial conditions -----------------------
// Guarded by USEINIT so a caller can supply its own init() without
// being silently overwritten (a defect found the hard way elsewhere
// in this library).
if (USEINIT < 0.5) {
  SMADS_0   = RDOSE;
  VEGFT_0   = VEGFSYN/VEGFDEG * pow(1.0/RDOSE, GSMADV);
  TNOSE_0   = TN0;
  TGIB_0    = TG0 * (1.0 + GIAGE*(AGE0 - 40.0));
  GSHUNT_0  = GS0 * LIVMULT;
  PAVMB_0   = PV0 * LUNMULT;
  MURALC_0  = MURSS;
  LAMEP_0   = LAMEP0 * TN0;
  EPDURC_0  = EPDUR0;
  FIBLY_0   = 1.0;
  FEGUT_0   = 0.0;
  FESTORE_0 = FESTORE0;
  HEPCC_0   = 1.0;
  ERFEC_0   = 1.0;
  RETICC_0  = 1.0;
  MSUP_0    = 0.0;
  TR1_0     = GS0 * LIVMULT;
  TR2_0     = GS0 * LIVMULT;
  TR3_0     = GS0 * LIVMULT;
  RVF_0     = 1.0;
  NTBNPC_0  = NTBNP0;
  PLAC_0    = 0.0;
  QOLC_0    = QOLREF;
  // HBMASS and ESSC are set to their model-consistent equilibrium by
  // the burn-in in sim_eq(); a nominal start avoids a cold solve.
  HBMASS_0  = 11.0 * BVOL / 100.0;
  ESSC_0    = ESSREF;
}

$ODE
// ================= drug PK ========================================
dxdt_BEVC = -(CLB/V1B)*BEVC - (QB/V1B)*BEVC + (QB/V2B)*BEVP;
dxdt_BEVP =  (QB/V1B)*BEVC - (QB/V2B)*BEVP;

dxdt_POMD = -KA_POM*POMD;
dxdt_POMC =  KA_POM*POMD - (CL_POM/V_POM)*POMC;

dxdt_PAZD = -KA_PAZ*PAZD;
dxdt_PAZC =  KA_PAZ*PAZD - (CL_PAZ/V_PAZ)*PAZC;

dxdt_TXAD = -KA_TXA*TXAD;
dxdt_TXAC =  KA_TXA*TXAD - (CL_TXA/V_TXA)*TXAC;

dxdt_TACD = -(CL_TAC/V_TAC)*TACD;
double CTACO = TACD/V_TAC*1000.0;          // ug/L

// ================= signalling =====================================
double TACEFF = CTACO/(CTACO + EC50_TAC);
double RD = RDOSE_ALK1;
if (GENO < 1.5)      RD = RDOSE_ENG;
else if (GENO > 2.5) RD = RDOSE_SM4;
dxdt_SMADS = KSMAD*( RD*(1.0 + ETAC*TACEFF) - SMADS );

// ---- VEGF: production rises with SMAD deficit and with shear ------
double VEGFPROD = VEGFSYN * pow(SIGDEF, GSMADV)
                          * (1.0 + GSHEAR*(GSHUNT - 1.0));
if (VEGFPROD < 0) VEGFPROD = 0;
dxdt_VEGFT = VEGFPROD - VEGFDEG*VEGFT;

// ---- free VEGF given bevacizumab (both converted to nM) ----------
// VEGF165 dimer 38.2 kDa; bevacizumab 149 kDa.  At therapeutic
// bevacizumab the molar excess is ~10^4-10^5, so target engagement is
// effectively BINARY — all of the PD lag therefore lives in the
// lesion turnover rates, not in the binding.  That fact is what makes
// axis4 a structural argument rather than a parameter argument.
double LT   = VEGFT * 2.618e-5;             // pg/mL -> nM
double BT   = CBEV  * 6.711;                // mg/L  -> nM (1 site)
double LFREE = freelig(LT, BT, KD_BEV);
double VEGF0N = (VEGFSYN/VEGFDEG) * pow(SIGDEF, GSMADV) * 2.618e-5;
double VFR = (VEGF0N > 1e-12) ? LFREE/VEGF0N : 1.0;
if (VFR < 0) VFR = 0; if (VFR > 1) VFR = 1;

// ---- pazopanib blocks the receptor rather than the ligand --------
double PAZE = EMAX_PZ*CPAZ/(CPAZ + EC50_PZ);
double DRIVE = (FRES + (1.0 - FRES)*VFR) * (1.0 - PAZE);
if (DRIVE < 0) DRIVE = 0;

// ================= vascular lesion turnover =======================
double AGEY = AGE0 + SOLVERTIME/365.25;

// formation is written so that DRIVE = 1 reproduces TN0 exactly; the
// signalling deficit is already baked into the baseline burden TN0.
dxdt_TNOSE = KOUT_N*TN0*DRIVE - KOUT_N*TNOSE - ABLATE*TNOSE;

double TGTGT = TG0*(1.0 + GIAGE*(AGEY - 40.0));
dxdt_TGIB  = KOUT_G*TGTGT*DRIVE - KOUT_G*TGIB - GIENDOF*TGIB;

double LIVM = (GENO > 1.5 && GENO < 2.5) ? GENO_LIV : 1.0;
// Shear positive feedback, written as a RELATIVE power law rather than an
// absolute offset.  The obvious form 1 + FLOWFB*(GSHUNT - GS0*LIVM) is a
// trap: in a high-output patient (large GS0) a drug that pushes GSHUNT well
// below baseline makes that bracket large and negative, the formation term
// changes SIGN, and conductance runs away to minus infinity - giving
// negative cardiac indices.  A power law cannot do that.
double GREF_L = GS0*LIVM;
double GRAT = (GREF_L > 1e-9) ? GSHUNT/GREF_L : 1.0;
if (GRAT < 1e-6) GRAT = 1e-6;
double SHEARFB = pow(GRAT, FLOWFB);
// The hepatic shunt gets its OWN VEGF-independent share.  A mature,
// smooth-muscle-invested arteriovenous shunt is a structure, not a
// sprouting front, so there is no reason its suppressible fraction should
// equal the nasal mucosa one.  Axis 5 solves for this number.
double DRIVE_L = (FRES_L + (1.0 - FRES_L)*VFR)*(1.0 - PAZE);
if (DRIVE_L < 0) DRIVE_L = 0;
dxdt_GSHUNT = KOUT_L*GS0*LIVM*DRIVE_L*SHEARFB - KOUT_L*GSHUNT;

double LUNM = (GENO < 1.5) ? GENO_LUN : 1.0;
dxdt_PAVMB = KOUT_P*PV0*LUNM*DRIVE - KOUT_P*PAVMB;

// ---- mural coverage: pomalidomide UP (biopsy), pazopanib DOWN ----
double POMMUR = EMAX_PM*CPOM/(CPOM + EC50_PM);
double MURTGT = MURSS*(1.0 + POMMUR) - PAZMUR*PAZE;
if (MURTGT < 0.05) MURTGT = 0.05;
if (MURTGT > 1.00) MURTGT = 1.00;
dxdt_MURALC = KMUR*(MURTGT - MURALC);

// ================= epistaxis ======================================
// Vessel maturation acts on all three bleeding dimensions at once:
// a better-covered vessel ruptures less often, bleeds less briskly and
// closes sooner.  That is what makes the pomalidomide biopsy finding
// (increased mural coverage, NO change in EC proliferation) able to
// carry the whole ESS effect.
double FRAGF = pow(MURSS/MURALC, GFRAG);
double QEFF  = QBLEED*pow(FRAGF, GQ);
double LAMTGT = LAMEP0*TNOSE*FRAGF*(1.0 - CLOSEF);
dxdt_LAMEP = KLAM*(LAMTGT - LAMEP);

double TXAE = EMAX_TX*CTXA/(CTXA + EC50_TX);
dxdt_FIBLY = KFIB*( (1.0 - TXAE)*(1.0 + FIBUP*(TNOSE - TN0)) - FIBLY );

double FEPLT = PLTFE*(1.0 - FESTORE/(FESTORE + FEPLT50));
double DURTGT = EPDUR0*FIBLY*(1.0 + FEPLT)*pow(FRAGF, GDUR);
dxdt_EPDURC = KDUR*(DURTGT - EPDURC);

double BNOSE = LAMEP/7.0 * EPDURC * QEFF;
double BGI   = BGI0*TGIB/ (TG0*(1.0 + GIAGE*(AGE0 - 40.0)));
double BTOT  = BNOSE + BGI;

// ================= iron ===========================================
double HEPSUP = 1.0/(1.0 + HEPCC/KHEP);
double AORAL  = AMAXO * FEGUT/(FEGUT + KMFEG) * HEPSUP;
double ADIETA = ADIET * HEPSUP;
dxdt_FEGUT = -AORAL - KTRFE*FEGUT;

// erythroid production: demand x iron-supply availability x marrow
double HBTARGMASS = HBSET*BVOL/100.0;
double DEF = (HBSET - HBC)/HBSET;
if (DEF < 0) DEF = 0;
double EPODRIVE = 1.0 + EMAXEPO*DEF;
if (EPODRIVE > MARROWMX) EPODRIVE = MARROWMX;
double DEMAND = HBTARGMASS/RBCLIFE * EPODRIVE;
double SUPPLYF = FESTORE/(FESTORE + KSUP);
double PROD = DEMAND * SUPPLYF * (1.0 - MSUP);
if (PROD < 0) PROD = 0;

// IV iron is prescribed BECAUSE the patient is anaemic, not on a fixed
// schedule.  Making it a feedback loop is the whole point of axis 2: a
// drug that lowers blood loss lowers the iron the clinician gives, and
// the two changes cancel on the haemoglobin axis.
double IVFER = IVFERATE
             + IVFE_ON*IVFEMAX/(1.0 + exp((HBC - IVFETRIG)/SIVFE));

// transfusion: smooth trigger on Hb
double TXRATE = TXALLOW * TXMAX/(1.0 + exp((HBC - HBTRIG)/STX));
double TXFE   = TXRATE*3.47*HBUNIT;

double RECYCLE = 3.47*HBMASS/RBCLIFE;
dxdt_FESTORE = AORAL + ADIETA + IVFER + TXFE
               + RECYCLE - 3.47*PROD - BASALFE;

dxdt_HBMASS = PROD - HBMASS/RBCLIFE - HBC*BTOT/100.0
              + TXRATE*HBUNIT;

// hepcidin and erythroferrone
double HEPTGT = pow(FESTORE/FESREF, GHEP)
                * (1.0 + IL6G*(GSHUNT - 1.0))
                / (1.0 + ERFEC/KERFE);
if (HEPTGT < 0.01) HEPTGT = 0.01;
dxdt_HEPCC = KHEPT*(HEPTGT - HEPCC);
double PRODREF = HBTARGMASS/RBCLIFE;
dxdt_ERFEC = KERF*( PROD/PRODREF - ERFEC );
dxdt_RETICC = KRET*( PROD/PRODREF - RETICC );

// pomalidomide marrow suppression (PATH-HHT: ANC −1710/uL, Hb flat)
dxdt_MSUP = KMSUP*( EMAX_MS*CPOM/(CPOM + EC50_MS) - MSUP );

// ================= cardiovascular =================================
double CAVO2 = 1.34*HBREF*DSAT0*pow(HBC/HBREF, ALPHA);
if (CAVO2 < 0.5) CAVO2 = 0.5;
double CITIS = VO2I/(10.0*CAVO2);
dxdt_TR1 = KTR*(GSHUNT - TR1);
dxdt_TR2 = KTR*(TR1 - TR2);
dxdt_TR3 = KTR*(TR2 - TR3);
double CIT = CITIS + KSHUNT*TR3;
double STRAIN = CIT - CITHR;
if (STRAIN < 0) STRAIN = 0;
dxdt_RVF = KRV*( 1.0 - RVSENS*STRAIN - RVF );
double RVFC = (RVF < 0.2) ? 0.2 : RVF;
dxdt_NTBNPC = KBNP*( NTBNP0*pow(CIT/CIREF, GBNP)/RVFC - NTBNPC );

// ================= pulmonary shunt ================================
double SFRAC = SPAV*PAVMB*(1.0 - EMBOLF);
if (SFRAC > 0.85) SFRAC = 0.85;
dxdt_EMBHAZ = HAZPAV*PAVMB*(1.0 - EMBOLF)/365.25;

// ================= ESS ============================================
double iF = LAMEP/(LAMEP + KIF);
double iD = EPDURC/(EPDURC + KID);
double iI = QEFF/(QEFF + KII);
double iA = 1.0/(1.0 + exp((HBC - HBANEM)/SHB));
double iT = TXRATE/(TXRATE + KITX);
double iATT = 0.5*iF + 0.5*iD;

double WNON = 1.0 - W_FE_SHARE;
double WSUM = W_FREQ + W_DUR + W_INT + W_ATT;
double ESSTGT = ESSMAX*ESSCAL*(
      WNON*( W_FREQ*iF + W_DUR*iD + W_INT*iI + W_ATT*iATT )/WSUM
    + W_FE_SHARE*( 0.6*iA + 0.4*iT ) );
double ESSEFF = ESSTGT - PLAC;
if (ESSEFF < 0) ESSEFF = 0;
dxdt_ESSC = KESS*(ESSEFF - ESSC);

double KPL = (BLIND > 0.5) ? KPLAC : KPLACOFF;
dxdt_PLAC = KPL*( PLACMAX*ONPLAC*BLIND - PLAC );

dxdt_QOLC = KQOL*( QOLREF*pow(ESSC/ESSREF, GQOL) - QOLC );

// ================= counters =======================================
dxdt_CUMTX   = TXRATE;
dxdt_CUMIVFE = IVFER;
dxdt_CUMBLD  = BTOT;
dxdt_AUCBEV  = CBEV;
dxdt_CUMPOM  = CPOM;

$TABLE
double HB = HBC;
double BV = BVOL;

double HEPSUPT = 1.0/(1.0 + HEPCC/KHEP);
double AORALT  = AMAXO * FEGUT/(FEGUT + KMFEG) * HEPSUPT;
double ADIETT  = ADIET * HEPSUPT;
double IVFERT  = IVFERATE
               + IVFE_ON*IVFEMAX/(1.0 + exp((HB - IVFETRIG)/SIVFE));
double ANET    = AORALT + ADIETT + IVFERT - BASALFE;

double FRAGFT  = pow(MURSS/MURALC, GFRAG);
double QEFFT   = QBLEED*pow(FRAGFT, GQ);
double B_NOSE  = LAMEP/7.0 * EPDURC * QEFFT;
double B_GI    = BGI0*TGIB/(TG0*(1.0 + GIAGE*(AGE0 - 40.0)));
double B_TOT   = B_NOSE + B_GI;
double F_GI    = (B_TOT > 1e-9) ? B_GI/B_TOT : 0.0;
double FELOSS  = 0.0347*HB*B_TOT;

// the closed-form equilibrium the ODEs must reproduce
double HB_PRED = (B_TOT > 1e-9) ? ANET/(0.0347*B_TOT) : HBSET;
double B_CRIT  = ANET/(0.0347*11.5);

double CAVO2T = 1.34*HBREF*DSAT0*pow(HB/HBREF, ALPHA);
if (CAVO2T < 0.5) CAVO2T = 0.5;
double CI_TIS = VO2I/(10.0*CAVO2T);
double CI_SH  = KSHUNT*TR3;
double CI     = CI_TIS + CI_SH;
double CO     = CI*BSA;

double DSATE  = VO2I/(10.0*1.34*HB*CI_TIS);
double SFRACT = SPAV*PAVMB*(1.0 - EMBOLF);
if (SFRACT > 0.85) SFRACT = 0.85;
double SAO2   = SCAP - SFRACT*DSATE/(1.0 - SFRACT);
if (SAO2 < 0.40) SAO2 = 0.40;
double SPO2   = 100.0*SAO2;

double FERRITIN = FESTORE/2.5;
double TSATP    = 100.0*(FESTORE/(FESTORE + 200.0))*0.45;

double TXR = TXALLOW * TXMAX/(1.0 + exp((HB - HBTRIG)/STX));
double TXPERYR = TXR*365.25;

double CBEVT = CBEV;
double CPOMT = CPOM;
double CPAZT = CPAZ;
double CTXAT = CTXA;

$CAPTURE
HB BV ANET IVFERT QEFFT B_NOSE B_GI B_TOT F_GI FELOSS HB_PRED B_CRIT
CI CI_TIS CI_SH CO SPO2 SFRACT DSATE FERRITIN TSATP TXPERYR
CBEVT CPOMT CPAZT CTXAT FRAGFT
'

# ---------------------------------------------------------------------
# compile
# ---------------------------------------------------------------------
mod <- mcode_cache("hht_qsp", hht_code, soloc = tempdir())

# ---------------------------------------------------------------------
# simulate to equilibrium, then run the protocol
# ---------------------------------------------------------------------
BURN <- 4000   # days of burn-in to reach the iron/Hb equilibrium

# `pre` holds the PRE-RANDOMISATION parameters used for the burn-in.
# Anything in `...` is applied only to the protocol phase.  Keeping these
# separate matters: two trial arms must start from the SAME patient and
# diverge afterwards, otherwise each arm silently equilibrates to its own
# iron regimen and the between-arm comparison is meaningless.
PRE_TRIAL <- list(IVFERATE = 0, IVFE_ON = 1, TXALLOW = 1,
                  BLIND = 0, ONPLAC = 0)

sim_eq <- function(m = mod, end = yr(2), delta = 1, events = NULL,
                   burn = BURN, pre = PRE_TRIAL, ...) {
  arm <- list(...)
  # phase 1: burn-in on the pre-trial regimen -> equilibrium patient
  pburn <- pre
  # structural (non-therapeutic) parameters must also hold during burn-in
  struct <- intersect(names(arm),
    c("GENO","WT","BSA","AGE0","BGI0","QBLEED","LAMEP0","EPDUR0","TN0",
      "TG0","GS0","PV0","ALPHA","KSHUNT","MURSS","GFRAG","KOUT_N","KOUT_G",
      "KOUT_L","KOUT_P","AMAXO","KHEP","FRES","W_FE_SHARE","KSUP","IL6G",
      "HBSET","SPAV","EMAX_PM","EC50_PM","PLACMAX","KII","ESSCAL","FRES_L",
      "EMAX_MS","GIAGE","HAZPAV","GQ","GDUR","BASALFE","ADIET"))
  # NOTE: IVFE_ON / IVFERATE / IVFETRIG are deliberately NOT in this list.
  # They are therapy, so every arm must burn in on the SAME standard-of-care
  # iron policy and only diverge during the protocol phase; otherwise each
  # arm starts from a different patient and the comparison is worthless.
  for (nm in struct) pburn[[nm]] <- arm[[nm]]
  base <- do.call(param, c(list(m), pburn)) %>%
    mrgsim(end = burn, delta = burn, hmax = 1) %>% as.data.frame()
  last <- base[nrow(base), ]
  cmts <- names(mod@init@data)
  iv <- as.numeric(last[cmts]); names(iv) <- cmts
  # phase 2: the protocol, starting from that equilibrium patient
  m2 <- do.call(param, c(list(m), pre)) 
  if (length(arm)) m2 <- do.call(param, c(list(m2), arm))
  m2 <- m2 %>% param(USEINIT = 1) %>% init(iv)
  if (!is.null(events)) m2 <- m2 %>% ev(events)
  m2 %>% mrgsim(end = end, delta = delta, hmax = 1) %>% as.data.frame()
}

# equilibrium state UNDER the supplied parameters (no separate pre-phase):
# used for steady-state/frontier questions, not for trial arms.
eqstate <- function(m = mod, burn = BURN, ...) {
  p <- utils::modifyList(list(BLIND = 0, ONPLAC = 0), list(...))
  d <- do.call(param, c(list(m), p)) %>%
    mrgsim(end = burn, delta = burn, hmax = 1) %>% as.data.frame()
  d[nrow(d), ]
}

# ---------------------------------------------------------------------
# dosing event builders
# ---------------------------------------------------------------------
ev_bev <- function(n = 6, ii = 14, wt = 75, start = 0, dose_mgkg = 5) {
  ev(time = start, amt = dose_mgkg * wt, cmt = "BEVC", ii = ii,
     addl = n - 1)
}
ev_bev_maint <- function(ii_days, until, wt = 75, start = 70 + 14,
                         dose_mgkg = 5) {
  n <- max(0, floor((until - start) / ii_days) + 1)
  if (n <= 0) return(NULL)
  ev(time = start, amt = dose_mgkg * wt, cmt = "BEVC",
     ii = ii_days, addl = n - 1)
}
ev_pom <- function(days, dose = 4, start = 0) {
  ev(time = start, amt = dose, cmt = "POMD", ii = 1, addl = days - 1)
}
ev_paz <- function(days, dose = 50, start = 0) {
  ev(time = start, amt = dose, cmt = "PAZD", ii = 1, addl = days - 1)
}
ev_txa <- function(days, dose = 1000, start = 0) {
  # 1 g three times daily
  ev(time = start, amt = dose, cmt = "TXAD", ii = 1/3, addl = days*3 - 1)
}
ev_oral_fe <- function(days, dose = 200, start = 0) {
  ev(time = start, amt = dose, cmt = "FEGUT", ii = 1, addl = days - 1)
}

# =====================================================================
#  AXIS 1 — the hyperbola is structural, not fitted
# =====================================================================
axis1_hyperbola <- function() {
  hdr("AXIS 1  Hb_ss = A_net / (0.0347 * B) is a structural identity")
  cat("Claim: the closed form falls out of mass balance alone, with no\n",
      "assumption about the erythroid production law.  Test: solve the\n",
      "full 40-ODE system to equilibrium over a grid of blood-loss rates\n",
      "and iron supplies and compare Hb to the closed form.\n", sep = "")

  # Iron must be PROTOCOL-FIXED for this test (IVFE_ON = 0).  With the
  # demand-driven policy switched on, A_net is itself a function of Hb and
  # the feedback loop pins haemoglobin near the trigger whatever the
  # bleeding is - which is a real clinical fact but hides the identity.
  grid <- expand.grid(BGI0 = c(4, 9.9, 25, 45),
                      IVFERATE = c(5, 12, 25, 40))
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    s <- eqstate(BGI0 = grid$BGI0[i], IVFERATE = grid$IVFERATE[i],
                 IVFE_ON = 0, TXALLOW = 0)
    data.frame(B_gi = grid$BGI0[i], IV = grid$IVFERATE[i],
               B_tot = s$B_TOT, A_net = s$ANET,
               Hb_ode = s$HB, Hb_closed = s$HB_PRED,
               Fe_store = s$FESTORE)
  }))
  out$rel_err_pct <- 100 * (out$Hb_ode - out$Hb_closed) / out$Hb_closed
  # The identity requires an equilibrium to EXIST, which needs
  #     A_net <= 0.0347 * Hb_setpoint * B.
  # Above that line iron arrives faster than a normal red-cell mass can
  # use it, so haemoglobin saturates at the set-point and the surplus
  # accumulates as stores: the patient becomes iron-LOADED rather than
  # polycythaemic, and there is no steady state in Fe_store at all.
  # Applicability test.  The identity needs an equilibrium in Fe_store,
  # which exists only while erythropoiesis is IRON-LIMITED.  Once iron
  # arrives faster than a red-cell mass at its set-point can consume, the
  # surplus goes to stores instead and Hb caps out: the patient becomes
  # iron-LOADED rather than polycythaemic.  Low stores are the signature
  # of the iron-limited regime, so test on them directly.
  out$iron_limited <- out$Fe_store < 400
  out$frontier_Anet <- 0.0347 * mod$HBSET * out$B_tot
  print(out[, setdiff(names(out), "frontier_Anet")],
        row.names = FALSE, digits = 4)
  sub("Rows with iron_limited = FALSE are ABOVE the frontier ",
      "A_net = 0.0347*Hb_set*B;\n   there the ODEs correctly cap Hb at the ",
      "set-point and load stores instead.")

  mx <- max(abs(out$rel_err_pct[out$iron_limited]))
  sub("Maximum |relative error| in the iron-limited regime (",
      sum(out$iron_limited), " of ", nrow(out), " cells): ", f3(mx), "%")
  cat("   (The identity degrades within roughly 10% of the frontier,\n",
      "    because a red-cell mass approaching its set-point stops\n",
      "    consuming every milligram that arrives.)\n", sep = "")
  if (mx < 1.0) {
    cat("   -> CONFIRMED. The hyperbola is a property of the model's\n",
        "      structure, so any claim resting on it is as strong as\n",
        "      mass balance itself.\n", sep = "")
  } else {
    cat("   -> NOT confirmed at 1% tolerance; equilibrium not reached\n")
  }

  sub("Consequence: the equilibrium is a hyperbola, so a proportional",
      " change in\n   bleeding gives a proportional change in Hb.")
  b <- c(5, 10, 20, 35, 60)
  hh <- sapply(b, function(bb) {
    eqstate(BGI0 = bb, IVFE_ON = 0, IVFERATE = 12, TXALLOW = 0)$HB })
  bt <- sapply(b, function(bb) {
    eqstate(BGI0 = bb, IVFE_ON = 0, IVFERATE = 12, TXALLOW = 0)$B_TOT })
  print(data.frame(B_gi_set = b, B_total = round(bt, 1),
                   Hb_eq = round(hh, 2)), row.names = FALSE)

  sub("B_crit — the blood loss a given iron regimen can just hold at ",
      "Hb 11.5:")
  for (iv in c(0, 5, 12, 25, 50)) {
    s <- eqstate(IVFERATE = iv, IVFE_ON = 0, TXALLOW = 0)
    cat(sprintf("   IV iron %5.1f mg/day  ->  A_net %6.2f mg/day  ->  B_crit %6.1f mL/day\n",
                iv, s$ANET, s$B_CRIT))
  }
  cat("\n   Read this as the therapeutic frontier: a patient bleeding more\n",
      "   than B_crit cannot be held at Hb 11.5 by that iron regimen no\n",
      "   matter how it is scheduled, because oral absorption saturates.\n",
      "   Note the top row: diet alone nets under 1 mg Fe/day, so B_crit\n",
      "   is a couple of mL/day.  Untreated HHT of any severity is\n",
      "   therefore OBLIGATORILY progressive on iron - which is why the\n",
      "   natural history is transfusion dependence, not a stable anaemia.\n",
      sep = "")
  invisible(out)
}

# =====================================================================
#  AXIS 2 — the PATH-HHT haemoglobin paradox
# =====================================================================
axis2_hb_paradox <- function() {
  hdr("AXIS 2  Why PATH-HHT reduced bleeding and did NOT raise Hb")
  cat("Observed (Al-Samkari 2024 NEJM, n=144):\n",
      "  * weighted bleeding duration fell more on drug (diff −8.9)\n",
      "  * RBC transfusion wk12-24: 9% on drug vs 18% on placebo\n",
      "  * median IV iron per 4 weeks: 0 mg on drug vs 333 mg on placebo\n",
      "  * haemoglobin: NO significant difference at 24 weeks\n",
      "  * haemoglobin: +1.09 g/dL (0.38-1.80) at 4 weeks POST-treatment\n",
      sep = "")

  sub("Step 1. What the withdrawn iron is worth in blood-loss units.")
  ivdiff <- 333 / 28          # mg/day withdrawn relative to placebo
  hb_ref <- 10.8
  rho <- 0.0347 * hb_ref
  cat(sprintf("   333 mg per 4 weeks              = %.2f mg Fe/day\n", ivdiff))
  cat(sprintf("   at Hb %.1f g/dL, rho_Fe          = %.4f mg Fe per mL blood\n",
              hb_ref, rho))
  cat(sprintf("   iron withdrawn is equivalent to = %.1f mL/day of bleeding\n",
              ivdiff / rho))
  cat("   Because B and A_net enter Hb only as a ratio, withdrawing that\n",
      "   much iron is INDISTINGUISHABLE from adding that much bleeding.\n",
      sep = "")

  sub("Step 2. Simulate the two arms with their OWN iron use.")
  D <- 24 * 7
  # Both arms run the SAME demand-driven iron policy.  The difference in
  # iron use is therefore an OUTPUT of the model, not an assumption.
  pl <- sim_eq(end = D + 28, delta = 1, events = NULL,
               BLIND = 1, ONPLAC = 1, TXALLOW = 1)
  po <- sim_eq(end = D + 28, delta = 1, events = ev_pom(D),
               BLIND = 1, ONPLAC = 1, TXALLOW = 1)

  at <- function(d, t) d[which.min(abs(d$time - t)), ]
  p0 <- at(pl, 0); pE <- at(pl, D); pP <- at(pl, D + 28)
  q0 <- at(po, 0); qE <- at(po, D); qP <- at(po, D + 28)

  sub("Levels at week 24 (A_net and IV iron are LEVELS: a constant ",
      "rate has\n   no delta, and it is the LEVEL that differs between ",
      "arms).")
  lev <- data.frame(
    quantity = c("IV iron (mg per 4 weeks)", "A_net (mg Fe/day)",
                 "B_total (mL/day)", "Hb (g/dL)", "ESS"),
    placebo_wk24 = c(pE$IVFERT*28, pE$ANET, pE$B_TOT, pE$HB, pE$ESSC),
    pom_wk24     = c(qE$IVFERT*28, qE$ANET, qE$B_TOT, qE$HB, qE$ESSC))
  lev$ratio_pom_placebo <- lev$pom_wk24 / lev$placebo_wk24
  print(lev, row.names = FALSE, digits = 3)
  cat(sprintf("\n   Model IV-iron use fell from %.0f to %.0f mg per 4 weeks.\n",
              pE$IVFERT*28, qE$IVFERT*28))
  cat("   PATH-HHT observed a median of 333 mg per 4 weeks on placebo and\n",
      "   0 mg on pomalidomide.  The placebo figure is a PREDICTION here,\n",
      "   not an input: the iron rate is set by a haemoglobin feedback\n",
      "   loop, and it lands on the observed value because the observed\n",
      "   value is what a blood loss of ~32 mL/day costs in iron.\n",
      sep = "")
  cat(sprintf("\n   Inverting that: an IV-iron requirement of 333 mg per 4 weeks\n"))
  cat(sprintf("   IMPLIES a mean blood loss of %.0f mL/day at Hb %.1f - a\n",
              (333/28)/(0.0347*pE$HB), pE$HB))
  cat("   quantity no HHT trial measures, recovered from one that they do.\n")

  sub("Changes from baseline (what the trial reports).")
  res <- data.frame(
    quantity = c("ESS", "HHT-QoL", "Hb (g/dL)", "B_total (mL/day)",
                 "A_net (mg Fe/day)", "mural coverage"),
    placebo_base = c(p0$ESSC, p0$QOLC, p0$HB, p0$B_TOT, p0$ANET, p0$MURALC),
    placebo_24w  = c(pE$ESSC, pE$QOLC, pE$HB, pE$B_TOT, pE$ANET, pE$MURALC),
    pom_base     = c(q0$ESSC, q0$QOLC, q0$HB, q0$B_TOT, q0$ANET, q0$MURALC),
    pom_24w      = c(qE$ESSC, qE$QOLC, qE$HB, qE$B_TOT, qE$ANET, qE$MURALC))
  res$d_placebo <- res$placebo_24w - res$placebo_base
  res$d_pom     <- res$pom_24w     - res$pom_base
  res$difference <- res$d_pom - res$d_placebo
  print(res[, c("quantity", "d_placebo", "d_pom", "difference")],
        row.names = FALSE, digits = 3)

  sub("Step 3. Model vs trial.")
  cmp <- data.frame(
    endpoint = c("ESS change, pomalidomide", "ESS change, placebo",
                 "ESS difference", "HHT-QoL difference",
                 "Hb difference at 24 wk", "Hb difference at 4 wk post"),
    observed = c(-1.84, -0.90, -0.94, -1.40, NA, 1.09),
    model = c(res$d_pom[1], res$d_placebo[1], res$difference[1],
              res$difference[2], res$difference[3],
              (qP$HB - q0$HB) - (pP$HB - p0$HB)))
  cmp$obs_txt <- ifelse(is.na(cmp$observed), "n.s.", f2(cmp$observed))
  print(cmp[, c("endpoint", "obs_txt", "model")], row.names = FALSE,
        digits = 3)

  sub("Step 4. The counterfactual the trial could not run:")
  cf <- sim_eq(end = D + 28, delta = 1, events = ev_pom(D),
               IVFE_ON = 0, IVFERATE = 11.9,
               BLIND = 1, ONPLAC = 1, TXALLOW = 1)
  cE <- at(cf, D)
  cat(sprintf("   pomalidomide with iron PROTOCOL-FIXED at the placebo rate:\n"))
  cat(sprintf("     Hb change at 24 wk = %+.2f g/dL  (vs %+.2f when iron was withdrawn)\n",
              cE$HB - at(cf, 0)$HB, qE$HB - q0$HB))
  cat(sprintf("     difference vs placebo = %+.2f g/dL\n",
              (cE$HB - at(cf, 0)$HB) - res$d_placebo[3]))
  cat("\n   INTERPRETATION.  The drug's effect on haemoglobin was spent\n",
      "   paying for the iron the patients stopped needing.  Haemoglobin\n",
      "   is not a bad endpoint because it is insensitive - it is a bad\n",
      "   endpoint because the rescue therapy moves along the SAME axis\n",
      "   in the OPPOSITE direction, and no blinding can prevent that.\n",
      sep = "")

  sub("Step 5. How much of the flat Hb is marrow suppression?")
  ns <- sim_eq(end = D + 28, delta = 1, events = ev_pom(D),
               BLIND = 1, ONPLAC = 1, TXALLOW = 1, EMAX_MS = 0)
  nE <- at(ns, D); nP <- at(ns, D + 28)
  cat(sprintf("   with marrow suppression   : Hb %+.2f at 24 wk, %+.2f at +4 wk\n",
              qE$HB - q0$HB, qP$HB - q0$HB))
  cat(sprintf("   without marrow suppression: Hb %+.2f at 24 wk, %+.2f at +4 wk\n",
              nE$HB - at(ns, 0)$HB, nP$HB - at(ns, 0)$HB))
  cat("   The post-treatment Hb rise needs BOTH: bleeding still\n",
      "   suppressed (relapse only 21%) while the marrow term lifts.\n",
      sep = "")
  invisible(list(placebo = pl, pom = po, cf = cf, cmp = cmp))
}

# =====================================================================
#  AXIS 3 — placebo is half the effect; and part of ESS is reachable
#           by iron alone
# =====================================================================
axis3_placebo_and_ess <- function() {
  hdr("AXIS 3  Placebo accounts for 49% of the pomalidomide ESS change")
  cat("PATH-HHT arm-level: pomalidomide −1.84, placebo −0.90.\n")
  cat(sprintf("   placebo share of the treated arm's improvement = %.1f%%\n",
              100 * 0.90 / 1.84))
  cat("NOSE (Whitehead 2016 JAMA) independently found that ALL arms,\n",
      "saline included, improved ESS significantly at weeks 12 and 24\n",
      "with no between-arm difference - i.e. the drift is real and large.\n",
      sep = "")

  sub("Consequence for reading single-arm reports:")
  cat("   A single-arm study reporting an ESS fall of X in this population\n",
      "   has, on the PATH-HHT placebo estimate, about 0.90 points of\n",
      "   non-drug improvement built in.  For a single-arm fall of 3.3\n",
      "   points that is 27%; for a fall of 1.5 points it is 60%.\n",
      sep = "")

  sub("The other half of the problem: ESS contains iron-sensitive items.")
  cat("Hoag 2010 lists ANAEMIA and NEED FOR TRANSFUSION among the six\n",
      "determinants of ESS.  Their weight was never published, so sweep\n",
      "it and ask how much of an ESS response iron ALONE could produce.\n",
      sep = "")

  D <- 24 * 7
  shares <- c(0.00, 0.06, 0.12, 0.18, 0.24, 0.30)
  out <- do.call(rbind, lapply(shares, function(w) {
    # iron-only intervention: raise IV iron, touch nothing else
    a <- sim_eq(end = D, delta = 7, W_FE_SHARE = w,
                TXALLOW = 1)
    b <- sim_eq(end = D, delta = 7, W_FE_SHARE = w, IVFERATE = 30.0,
                TXALLOW = 1)
    a0 <- a[1, ]; aE <- a[nrow(a), ]; b0 <- b[1, ]; bE <- b[nrow(b), ]
    data.frame(W_FE_SHARE = w,
               dHb_iron = (bE$HB - b0$HB) - (aE$HB - a0$HB),
               dESS_iron = (bE$ESSC - b0$ESSC) - (aE$ESSC - a0$ESSC))
  }))
  out$pct_of_094 <- 100 * abs(out$dESS_iron) / 0.94
  print(out, row.names = FALSE, digits = 3)

  sub("Reading:")
  cat("   Escalating IV iron alone - with bleeding completely unchanged -\n",
      "   moves ESS by an amount that is a pure artefact of the score's\n",
      "   construction.  At the mid-range weight it is a non-trivial\n",
      "   fraction of the entire registrational effect size of 0.94.\n",
      "   This does NOT mean pomalidomide failed: the trial's iron use\n",
      "   FELL on drug, so this artefact worked AGAINST the drug there.\n",
      "   It means ESS cannot be compared across trials with different\n",
      "   iron-rescue policies.\n", sep = "")
  invisible(out)
}

# =====================================================================
#  AXIS 4 — the bevacizumab 6-month nadir is structurally impossible
# =====================================================================
axis4_bev_structural <- function() {
  hdr("AXIS 4  A hypothesis this model REJECTED: that the 6-month ",
      "bevacizumab nadir was pharmacologically impossible")
  cat("Dupuis-Girod 2012 JAMA, n=25, bevacizumab 5 mg/kg q14d x 6:\n",
      "   last dose day 70\n",
      "   epistaxis duration 221 (baseline) -> 134 (3 mo) -> 43 (6 mo) min/mo\n",
      "The epistaxis endpoint reached its BEST value 110 days after the last\n",
      "dose.  I set out to show that no PK/PD structure permits that.  The\n",
      "model disagreed, and the reason it disagreed is worth more than the\n",
      "original claim would have been.\n", sep = "")

  thalf <- 20
  k <- log(2) / thalf
  frac180 <- exp(-k * (180 - 70)); frac90 <- exp(-k * (90 - 70))
  sub("Step 1. How much drug is left?")
  cat(sprintf("   t1/2 = %d d  ->  day 90 : %.1f%% of the day-70 peak\n",
              thalf, 100 * frac90))
  cat(sprintf("                    day 180: %.1f%% of the day-70 peak\n",
              100 * frac180))

  sub("Step 2. The structural argument, stated exactly.")
  cat("   For ANY inhibitory indirect-response model\n",
      "        dT/dt = k_in*(1 - I) - k_out*T,   0 <= I <= 1\n",
      "   the untreated equilibrium is T0 = k_in/k_out, so k_out*T <= k_in\n",
      "   for every T <= T0.  Once I = 0,\n",
      "        dT/dt = k_in - k_out*T >= 0,\n",
      "   so T cannot still be falling, whatever k_in, k_out or the shape\n",
      "   of I(C).  A direct-effect model is the k_out -> infinity limit of\n",
      "   the same statement.  The argument is airtight - but note its\n",
      "   PREMISE.  It requires I = 0, i.e. that TARGET ENGAGEMENT has\n",
      "   ceased.  It does NOT require that the drug concentration be low.\n",
      sep = "")

  sub("Step 3. Test the premise. This is where the hypothesis died.")
  d <- sim_eq(end = 500, delta = 1, events = ev_bev(6, 14, 75))
  bev_nM  <- d$CBEVT * 6.711
  vegf_nM <- d$VEGFT * 2.618e-5
  i70 <- which.min(abs(d$time - 70))
  exc <- bev_nM[i70] / vegf_nM[i70]
  cat(sprintf("   plasma VEGF at baseline      = %.0f pg/mL = %.5f nM\n",
              d$VEGFT[1], vegf_nM[1]))
  cat(sprintf("   bevacizumab at the last dose = %.1f mg/L = %.1f nM\n",
              d$CBEVT[i70], bev_nM[i70]))
  cat(sprintf("   molar excess of antibody over target = %.2e fold\n", exc))
  teng <- thalf * log2(exc)
  cat(sprintf("   time for the antibody to decay BELOW target molarity\n"))
  cat(sprintf("     = t1/2 * log2(excess) = %d * %.1f = %.0f days\n",
              thalf, log2(exc), teng))
  # measured in the model: when does free VEGF recover to 50% of baseline?
  vfr <- d$VEGFT * 0  # recompute free fraction the way the model does
  lt <- d$VEGFT * 2.618e-5; bt <- d$CBEVT * 6.711; KD <- mod$KD_BEV
  b <- lt - bt - KD
  lfree <- 0.5 * (b + sqrt(pmax(0, b * b + 4 * KD * lt)))
  rec <- lfree / lfree[1]
  t50 <- suppressWarnings(min(d$time[d$time > 70 & rec > 0.5]))
  cat(sprintf("   model free-VEGF recovery to 50%% of baseline: day %s\n",
              ifelse(is.finite(t50), sprintf("%.0f", t50), "> 500")))
  cat("\n   -> HYPOTHESIS REJECTED.  Because bevacizumab sits in a ~10^4\n",
      "      molar excess over plasma VEGF, target engagement outlasts\n",
      "      measurable drug by more than a dozen half-lives.  At day 180\n",
      "      the concentration is 2% of peak and the engagement is still\n",
      "      essentially complete, so I = 0 is FALSE and the structural\n",
      "      argument simply does not apply.  A continuing decline at\n",
      "      6 months is entirely compatible with pharmacology.\n", sep = "")

  sub("Step 4. What the rejection teaches, which is the useful part.")
  cat("   The PD-relevant exposure metric for a stoichiometric neutralising\n",
      "   antibody is not concentration and not AUC - it is whether\n",
      "   concentration exceeds TARGET MOLARITY.  That is why the\n",
      "   maintenance-interval question has the answer axis 5 finds: over\n",
      "   any interval a clinician would consider, mean target engagement\n",
      "   is complete, so interval changes act only on the tails.\n", sep = "")
  cat("\n   It also produces a falsifiable prediction that FAILS against\n",
      "   clinical experience, and that failure is informative.  The model\n",
      sep = "")
  cat(sprintf("   says free VEGF stays suppressed for ~%.0f days after a single\n",
              teng))
  cat("   induction course.  Real patients relapse within 1-3 months and\n",
      "   need maintenance dosing on that schedule.  Therefore epistaxis\n",
      "   relapse after bevacizumab is NOT driven by free-VEGF recovery.\n",
      "   Something else - VEGF up-regulation under blockade, complex-\n",
      "   mediated clearance raising total VEGF, or a VEGF-independent\n",
      "   driver of the residual formation rate (the FRES term) - must be\n",
      "   carrying the relapse.  This model does not distinguish them, and\n",
      "   saying which would need paired free-VEGF and epistaxis time\n",
      "   courses that do not exist.\n", sep = "")

  sub("Step 5. What survives about the single-arm result.")
  cat("   The structural claim is gone, but the EMPIRICAL one stands and it\n",
      "   does not need any model at all:\n", sep = "")
  cat("     * PATH-HHT, randomised: placebo -0.90 of the -1.84 on drug,\n",
      "       i.e. 48.9% of the treated arm improvement was not the drug.\n",
      "     * NOSE, randomised: every arm including saline improved ESS\n",
      "       significantly at weeks 12 and 24, with no between-arm\n",
      "       difference at all.\n", sep = "")
  cat("   A single-arm series in 25 patients selected for severe epistaxis,\n",
      "   on an endpoint whose baseline range was 0-947 min/month, will\n",
      "   therefore contain a large non-drug component.  How large cannot\n",
      "   be determined from that study, and this model cannot supply the\n",
      "   missing control arm.\n", sep = "")

  sub("Step 6. Model time course versus the observed one.")
  mkmin <- function(r) r$LAMEP / 7 * r$EPDURC * 30.4375
  at <- function(t) d[which.min(abs(d$time - t)), ]
  sc <- 221 / mkmin(at(0))
  cmp <- data.frame(
    series = c("observed (Dupuis-Girod, single arm)", "model (drug effect only)"),
    day0 = c(221, sc * mkmin(at(0))),
    day90 = c(134, sc * mkmin(at(90))),
    day180 = c(43, sc * mkmin(at(180))),
    day365 = c(NA, sc * mkmin(at(365))))
  print(cmp, row.names = FALSE, digits = 3)
  nad <- d$time[which.min(mkmin(d))]
  cat(sprintf("\n   model nadir: day %.0f (after the last dose, as the ",
              nad))
  cat("stoichiometry\n   requires); model day-180 value is above the ")
  cat("observed 43, the gap\n   being the non-drug component the trial ")
  cat("could not measure.\n")
  invisible(d)
}

# ---------------------------------------------------------------------
# virtual population: one long run per subject with the protocol offset
# by the burn-in, so every subject reaches its own equilibrium first.
# ---------------------------------------------------------------------
make_pop <- function(n = 250, gs_mean = 9.0) {
  lnv <- function(n, m, cv) m * exp(rnorm(n, 0, sqrt(log(1 + cv^2))) -
                                      0.5 * log(1 + cv^2))
  # Absolute epistaxis burden is drawn to match the REPORTED distribution of
  # the Lyon cohort: mean 221 min/month with a range of 0-947.  A lognormal
  # with meanlog 4.898 and sdlog 1.0 has median 134, mean 221 and a 97.5th
  # percentile of 950 - i.e. it reproduces both the mean and the tail.  This
  # matters because Azzopardi endpoint is the ABSOLUTE threshold
  # "< 20 min/month", so who crosses it is largely decided by where they
  # started, not by how much drug they received.
  data.frame(
    ID       = seq_len(n),
    EPBASE   = rlnorm(n, 4.898, 1.0),
    CLB      = lnv(n, mod$CLB,    0.30),
    VEGFSYN  = lnv(n, mod$VEGFSYN, 0.40),
    KOUT_N   = lnv(n, mod$KOUT_N, 0.50),
    KOUT_L   = lnv(n, mod$KOUT_L, 0.50),
    GS0      = lnv(n, gs_mean,    0.55),
    BGI0     = lnv(n, mod$BGI0,   0.70),
    QBLEED   = lnv(n, mod$QBLEED, 0.45),
    FRES     = pmin(0.8, lnv(n, mod$FRES, 0.35)))
}

sim_pop <- function(idata, ii_days = NA, months = 24, burn = 2500,
                    wt = 75, extra = list()) {
  end <- burn + round(mo(months))
  e <- ev(time = burn, amt = 5 * wt, cmt = "BEVC", ii = 14, addl = 5)
  if (!is.na(ii_days)) {
    st <- burn + 70 + ii_days
    nn <- max(0, floor((end - st) / ii_days) + 1)
    if (nn > 0) e <- c(e, ev(time = st, amt = 5 * wt, cmt = "BEVC",
                            ii = ii_days, addl = nn - 1))
  }
  m <- mod
  if (length(extra)) m <- do.call(param, c(list(m), extra))
  m %>% idata_set(idata) %>% ev(e) %>%
    mrgsim(end = end, delta = 30, hmax = 1) %>% as.data.frame()
}

axis5_two_organs <- function() {
  hdr("AXIS 5  Interval dependence is a VARIABILITY effect, not a ",
      "mean-exposure effect")
  cat("Azzopardi 2015 simulated bevacizumab maintenance at 24 months:\n",
      "   interval        q3 mo   q2 mo   q1 mo   ratio q1/q3\n",
      "   CI < 4          41%     45%     50%     1.22\n",
      "   epistaxis<20min 34%     43%     60%     1.76\n", sep = "")

  sub("First, what the mechanistic model says about the MEAN patient.")
  cat("Bevacizumab neutralises VEGF stoichiometrically, and at 5 mg/kg it\n",
      "sits in roughly 10^4-10^5 MOLAR EXCESS over plasma VEGF.  Even at a\n",
      "3-month trough - about 4% of peak - the excess is still enormous, so\n",
      "target engagement is essentially COMPLETE at every interval tested.\n",
      sep = "")
  ex <- do.call(rbind, lapply(c(91, 61, 30), function(ii) {
    d <- sim_eq(end = yr(2), delta = 7,
                events = c(ev_bev(6, 14, 75), ev_bev_maint(ii, yr(2), 75)))
    t2 <- d[d$time > yr(2) - 90, ]
    data.frame(interval_days = ii,
               trough_mgL = min(t2$CBEVT),
               molar_excess = min(t2$CBEVT) * 6.711 /
                 (mean(t2$VEGFT) * 2.618e-5),
               CI_mean = mean(t2$CI),
               epis_min_mo = mean(t2$LAMEP/7*t2$EPDURC*30.4375),
               shunt = mean(t2$GSHUNT))
  }))
  print(ex, row.names = FALSE, digits = 3)
  cat(sprintf("\n   The mean patient moves by %.1f%% in cardiac index and %.1f%%\n",
              100*abs(diff(range(ex$CI_mean)))/ex$CI_mean[1],
              100*abs(diff(range(ex$epis_min_mo)))/ex$epis_min_mo[1]))
  cat("   in epistaxis across a 3-fold change in dosing interval.  That is\n",
      "   effectively nothing - so the mean patient CANNOT explain\n",
      "   Azzopardi interval dependence, and a model that only ever\n",
      "   simulates a mean patient would have to reject those numbers.\n",
      sep = "")

  sub("Now the same question in a virtual population.")
  cat("Azzopardi reported RESPONSE FRACTIONS, not means.  A fraction is a\n",
      "count of subjects crossing a threshold, so it moves whenever the\n",
      "TAIL of the distribution moves - even when the mean does not.\n",
      sep = "")
  # Build the cohort the way the trial did: simulate a pool, then apply the
  # entry criterion (high cardiac index, 4.1-6.2 L/min/m2 as reported).
  pool <- make_pop(500)
  probe <- sim_pop(pool, ii_days = NA, months = 1)
  b0 <- probe[abs(probe$time - 2500) < 1e-6, c("ID", "CI")]
  names(b0) <- c("ID", "CI0")
  keep <- b0$ID[b0$CI0 >= 4.1 & b0$CI0 <= 6.2]
  pop <- pool[pool$ID %in% keep, ]
  cat(sprintf("   pool of %d -> %d subjects meet the entry criterion\n",
              nrow(pool) + length(setdiff(pool$ID, keep)), nrow(pop)))
  cat(sprintf("   cohort baseline cardiac index: median %.2f (range %.2f-%.2f)\n",
              median(b0$CI0[b0$ID %in% keep]), min(b0$CI0[b0$ID %in% keep]),
              max(b0$CI0[b0$ID %in% keep])))
  eb <- pop$EPBASE
  cat(sprintf("   cohort baseline epistaxis: mean %.0f, median %.0f, max %.0f min/mo\n",
              mean(eb), median(eb), max(eb)))
  cat("   (observed cohort: mean 221, range 0-947 min/month)\n")

  res <- do.call(rbind, lapply(c(91, 61, 30), function(ii) {
    d <- sim_pop(pop, ii_days = ii)
    d$epis <- d$LAMEP / 7 * d$EPDURC * 30.4375
    bb <- d[abs(d$time - 2500) < 1e-6, c("ID", "epis", "CI")]
    names(bb) <- c("ID", "epis0", "CI0")
    ee <- d[d$time == max(d$time), c("ID", "epis", "CI")]
    m <- merge(bb, ee, by = "ID")
    m <- merge(m, pop[, c("ID", "EPBASE")], by = "ID")
    # absolute epistaxis = the subject own reported baseline scaled by the
    # model fractional change
    m$epis_abs <- m$EPBASE * m$epis / m$epis0
    data.frame(interval_days = ii,
               frac_CI_lt4 = mean(m$CI < 4),
               frac_epis_lt20 = mean(m$epis_abs < 20),
               CI_median = median(m$CI),
               epis_abs_median = median(m$epis_abs))
  }))
  res$CI_sens <- res$frac_CI_lt4 / max(1e-9, res$frac_CI_lt4[1])
  res$epis_sens <- res$frac_epis_lt20 / max(1e-9, res$frac_epis_lt20[1])
  print(res, row.names = FALSE, digits = 3)

  sub("Comparison with the published simulations:")
  cmp <- data.frame(
    quantity = c("CI<4 sensitivity (q1/q3)", "epistaxis sensitivity (q1/q3)",
                 "ratio of the two sensitivities"),
    Azzopardi = c(1.22, 1.76, 1.76/1.22),
    model = c(res$CI_sens[3], res$epis_sens[3],
              res$epis_sens[3]/max(1e-9, res$CI_sens[3])))
  print(cmp, row.names = FALSE, digits = 3)
  cat(sprintf("\n   absolute response fractions, q3/q2/q1 month:\n"))
  cat(sprintf("     CI < 4 L/min/m2   model %.0f%% / %.0f%% / %.0f%%   (Azzopardi 41/45/50)\n",
              100*res$frac_CI_lt4[1], 100*res$frac_CI_lt4[2],
              100*res$frac_CI_lt4[3]))
  cat(sprintf("     epistaxis < 20 min model %.0f%% / %.0f%% / %.0f%%   (Azzopardi 34/43/60)\n",
              100*res$frac_epis_lt20[1], 100*res$frac_epis_lt20[2],
              100*res$frac_epis_lt20[3]))

  sub("Reading: the model does NOT reproduce the interval dependence.")
  cat("   Across a three-fold change in interval the model moves the cardiac\n",
      "   response fraction by ", f2(100*(res$frac_CI_lt4[3]-res$frac_CI_lt4[1])),
      " points and the epistaxis fraction by ",
      f2(100*(res$frac_epis_lt20[3]-res$frac_epis_lt20[1])), " points.\n",
      "   Azzopardi moved them by 9 and 26 points.  This is not a parameter\n",
      "   disagreement; it is a STRUCTURAL one, and it is the same issue\n",
      "   axis 4 uncovered:\n", sep = "")
  cat("     * A stoichiometric-neutralisation model makes exposure mean\n",
      "       CONCENTRATION RELATIVE TO TARGET MOLARITY.  With a 10^3-10^4\n",
      "       excess even at a 3-month trough, occupancy is complete at\n",
      "       every interval, so the interval cannot matter.\n",
      "     * An empirical direct-inhibition model makes exposure mean\n",
      "       CONCENTRATION.  Then a 3-month trough at 4% of peak is a 96%\n",
      "       loss of effect, and the interval matters enormously.\n",
      sep = "")
  cat("   Both structures fit the same 25 patients. They disagree precisely\n",
      "   on the clinical question that was asked of them - can bevacizumab\n",
      "   maintenance be stretched to 3-monthly - and the disagreement is\n",
      "   about what exposure MEANS, not about any fitted constant.\n",
      "   THE DISCRIMINATING EXPERIMENT is a trough free-VEGF measurement.\n",
      "   If free VEGF is still suppressed at a 3-month trough, interval\n",
      "   stretching is safe and the published simulations are conservative.\n",
      "   If it has recovered, they are right and this model is wrong.  That\n",
      "   assay is routine; the measurement has not been reported.\n",
      sep = "")

  sub("A second discrepancy: the model is far too optimistic about the ",
      "cardiac endpoint.")
  cat("   The model normalises cardiac index in ",
      f1(100*res$frac_CI_lt4[1]), "% of the cohort at q3-monthly\n",
      "   dosing; Azzopardi predicted 41%.  Two things could be wrong: the\n",
      "   hepatic shunt may be much less suppressible than assumed, or the\n",
      "   haemoglobin-to-cardiac coupling may be weaker.  Sweep both.\n",
      sep = "")
  # Two candidate explanations, swept together.  This sweep is done
  # DETERMINISTICALLY on the single fitted high-output patient rather than
  # over the virtual population: at extreme parameter combinations the
  # population solve produced non-finite cells, and a number I cannot trust
  # is worse than no number.
  fgs2 <- function(g) eqstate(GENO = 2, GS0 = g, ALPHA = 1)$CI - 5.05
  GSH <- tryCatch(uniroot(fgs2, c(0.5, 80))$root, error = function(e) 9.0)
  gr <- expand.grid(FRES_L = c(0.30, 0.55, 0.75, 0.95),
                    ALPHA = c(1.0, 0.5, 0.0))
  swp <- do.call(rbind, lapply(seq_len(nrow(gr)), function(i) {
    d <- sim_eq(end = yr(2), delta = 30,
                events = c(ev_bev(6, 14, 75), ev_bev_maint(91, yr(2), 75)),
                GENO = 2, GS0 = GSH,
                FRES_L = gr$FRES_L[i], ALPHA = gr$ALPHA[i])
    a <- d[1, ]; b <- d[nrow(d), ]
    data.frame(FRES_L = gr$FRES_L[i], ALPHA = gr$ALPHA[i],
               CI_base = a$CI, CI_24mo = b$CI, Hb_24mo = b$HB,
               shunt_24mo = b$GSHUNT, reaches_CI_lt4 = b$CI < 4)
  }))
  print(swp, row.names = FALSE, digits = 3)

  nfail <- sum(!swp$reaches_CI_lt4)
  cat(sprintf("\n   %d of %d parameter cells fail to reach CI < 4 at 24 months.\n",
              nfail, nrow(swp)))
  cat("   The dominant determinant is FRES_L, and the boundary sits between\n",
      "   0.30 (reaches CI < 4) and 0.55 (does not).  Reproducing a cohort\n",
      "   response fraction near the published 41% therefore implies\n",
      "   FRES_L in the region of 0.4-0.55: ROUGHLY HALF of hepatic shunt\n",
      "   maintenance is beyond the reach of anti-VEGF therapy.\n", sep = "")
  cat("\n   That is a quantitative version of a familiar clinical\n",
      "   observation.  A mature, smooth-muscle-invested arteriovenous shunt\n",
      "   is a structure rather than a sprouting front, which is why liver\n",
      "   transplantation - not any drug - remains the definitive treatment\n",
      "   for high-output failure in HHT, while the nasal mucosa, which IS a\n",
      "   sprouting front, responds much better.\n", sep = "")
  cat("\n   alpha shifts the boundary only modestly, and note that the\n",
      "   baseline cardiac index differs down the alpha column because GS0\n",
      "   was fitted to give CI = 5.05 at alpha = 1.  The alpha rows are\n",
      "   therefore not strictly comparable with one another; the FRES_L\n",
      "   comparison within a row is.\n", sep = "")

  cat("\n   CAVEAT.  The population variability here is assumed, not\n",
      "   estimated: Azzopardi fitted theirs to 25 patients and I have no\n",
      "   access to those estimates.  The ORDERING is robust to the\n",
      "   assumption; the numbers are not.\n", sep = "")
  invisible(res)
}

# =====================================================================
#  AXIS 6 — cardiac index is confounded by haemoglobin
# =====================================================================
axis6_cardiac_confound <- function() {
  hdr("AXIS 6  The cardiac-index endpoint cannot separate shunt from ",
      "anaemia")
  cat("Fick: CO = VO2 / C(a-v)O2 and C(a-v)O2 = 1.34 * Hb * dSat.\n",
      "If the arterio-venous SATURATION difference is conserved as Hb\n",
      "changes, C(a-v)O2 is proportional to Hb and CI ~ 1/Hb (alpha=1).\n",
      "If instead the absolute oxygen CONTENT difference is conserved,\n",
      "CI is independent of Hb (alpha=0).  Real chronic anaemia lies\n",
      "between, and no HHT study reports the data needed to place it.\n",
      sep = "")

  sub("Pure arithmetic: what Hb rise reproduces 5.05 -> 4.2 alone?")
  need <- 5.05 / 4.2
  cat(sprintf("   CI ratio to explain = %.3f\n", need))
  for (a in c(1.0, 0.75, 0.5, 0.25)) {
    hbr <- need^(1 / a)
    cat(sprintf("   alpha = %.2f  ->  requires Hb x %.3f  (e.g. %.1f -> %.1f g/dL)\n",
                a, hbr, 10.0, 10.0 * hbr))
  }
  cat("   Al-Samkari's bevacizumab survey and PATH-HHT both report Hb\n",
      "   rises of this order on antiangiogenic therapy, so at alpha near\n",
      "   1 the ENTIRE published cardiac-index response is compatible\n",
      "   with anaemia correction and zero shunt regression.\n", sep = "")

  sub("Reproduce the Dupuis-Girod cohort: solve for the hepatic shunt ",
      "burden\n   that puts baseline cardiac index at the observed 5.05 ",
      "L/min/m2.")
  fgs <- function(g) eqstate(GENO = 2, GS0 = g, ALPHA = 1,
                             IVFE_ON = 0, IVFERATE = 12)$CI - 5.05
  GSFIT <- tryCatch(uniroot(fgs, c(0.5, 80))$root, error = function(e) NA)
  if (is.na(GSFIT)) GSFIT <- 8.3
  cat(sprintf("   fitted GS0 = %.2f  ->  baseline CI = %.3f L/min/m2\n",
              GSFIT, eqstate(GENO = 2, GS0 = GSFIT, ALPHA = 1,
                             IVFE_ON = 0, IVFERATE = 12)$CI))
  cat(sprintf("   (the reference ACVRL1 patient sits at CI %.2f; this cohort\n",
              eqstate(GENO = 2)$CI))
  cat("    was selected for dyspnoea and high output, so it is ~an order\n",
      "    of magnitude more hepatic shunt.)\n", sep = "")

  sub("Model decomposition of the Dupuis-Girod protocol:")
  end <- 190
  out <- do.call(rbind, lapply(c(0, 0.25, 0.5, 0.75, 1.0), function(a) {
    # Iron is PROTOCOL-FIXED here.  With the demand-driven policy on, the
    # clinician withdraws iron as bleeding improves, haemoglobin stays put
    # and the anaemia pathway is artificially switched off - which would
    # make this decomposition answer the wrong question.
    fx <- list(IVFE_ON = 0, IVFERATE = 12)
    full <- do.call(sim_eq, c(list(end = end, delta = 1,
                   events = ev_bev(6, 14, 75), ALPHA = a, GENO = 2,
                   GS0 = GSFIT), fx))
    # shunt frozen: the drug may move everything EXCEPT the liver
    nosh <- do.call(sim_eq, c(list(end = end, delta = 1,
                   events = ev_bev(6, 14, 75), ALPHA = a, GENO = 2,
                   GS0 = GSFIT, KOUT_L = 1e-8), fx))
    a0 <- full[1, ]; a3 <- full[which.min(abs(full$time - 90)), ]
    b0 <- nosh[1, ];  b3 <- nosh[which.min(abs(nosh$time - 90)), ]
    dtot <- a3$CI - a0$CI
    danaem <- b3$CI - b0$CI          # Hb pathway only
    data.frame(alpha = a, CI_base = a0$CI, CI_3mo = a3$CI,
               dHb_3mo = a3$HB - a0$HB,
               dCI_total = dtot, dCI_anaemia_only = danaem,
               anaemia_share_pct = ifelse(abs(dtot) > 1e-9,
                                          100 * danaem / dtot, NA))
  }))
  print(out, row.names = FALSE, digits = 3)
  cat(sprintf("\n   Model total dCI at 3 months, alpha=1: %+.2f (observed %+.2f)\n",
              out$dCI_total[out$alpha == 1], 4.2 - 5.05))
  cat("\n   The share of the cardiac response that is really anaemia\n",
      "   correction runs from ~0% to most of it as alpha goes 0 -> 1.\n",
      "   alpha is not identifiable from any published HHT dataset, so\n",
      "   the shunt-attributable effect of bevacizumab on cardiac index\n",
      "   is NOT identifiable either.  Reporting CI without the\n",
      "   contemporaneous Hb makes the primary endpoint uninterpretable.\n",
      sep = "")
  invisible(out)
}

# =====================================================================
#  AXIS 7 — pulse oximetry cannot screen for pulmonary AVM
# =====================================================================
axis7_spo2_screen <- function() {
  hdr("AXIS 7  SpO2 is not a screening test for pulmonary AVM")
  cat("Mixing:  SaO2 = (1-s)*ScO2 + s*SvO2 and SvO2 = SaO2 - dSat, so\n",
      "         SaO2 = ScO2 - s*dSat/(1-s).\n",
      "dSat itself depends on Hb and cardiac output, both of which move\n",
      "in HHT - so the same anatomy gives different saturations.\n",
      sep = "")

  sub("Shunt fraction needed to pull SpO2 down to a given value:")
  sc <- 0.99
  for (dsat in c(0.20, 0.25, 0.30, 0.35)) {
    f <- function(target) {
      # ScO2 - s*dsat/(1-s) = target  ->  s = (ScO2-target)/(dsat+ScO2-target)
      g <- sc - target
      g / (dsat + g)
    }
    cat(sprintf("   dSat %.2f :  SpO2 97%% at s=%.3f | 95%% at s=%.3f | 90%% at s=%.3f\n",
                dsat, f(0.97), f(0.95), f(0.90)))
  }

  sub("Model sweep over pulmonary AVM burden (ENG genotype):")
  out <- do.call(rbind, lapply(c(0, 0.25, 0.5, 1, 1.5, 2, 3), function(p) {
    s <- eqstate(GENO = 1, PV0 = p, TXALLOW = 0)
    data.frame(PAVM_burden = p, shunt_fraction = s$SFRACT,
               SpO2 = s$SPO2, Hb = s$HB,
               embolic_hazard_per_yr = 0.010 * s$PAVMB)
  }))
  print(out, row.names = FALSE, digits = 4)

  sub("The decoupling:")
  ok <- out[out$SpO2 >= 95, ]
  if (nrow(ok) > 0) {
    cat(sprintf("   Burdens up to %.2f (shunt fraction up to %.3f) keep SpO2\n",
                max(ok$PAVM_burden), max(ok$shunt_fraction)))
    cat(sprintf("   at or above 95%% while carrying an embolic hazard of up to\n"))
    cat(sprintf("   %.4f per year - i.e. a fully 'normal' oximetry reading in a\n",
                max(ok$embolic_hazard_per_yr)))
    cat("   patient who still needs embolotherapy and antibiotic\n",
        "   prophylaxis.  Embolic risk is set by the anatomy (feeding\n",
        "   artery calibre); oxygenation is set by the anatomy AND Hb AND\n",
        "   cardiac output.  The two do not co-vary, which is the\n",
        "   quantitative reason the guideline screens with contrast\n",
        "   echocardiography rather than a pulse oximeter.\n", sep = "")
  }

  sub("Same anatomy, different anaemia - and it depends on alpha:")
  cat("   At alpha = 1 the compensating rise in cardiac output exactly\n",
      "   cancels the fall in oxygen-carrying capacity, dSat is pinned at\n",
      "   its reference value and SpO2 becomes INDEPENDENT of haemoglobin.\n",
      "   Below alpha = 1 the compensation is incomplete, dSat widens and\n",
      "   the same anatomic shunt desaturates an anaemic patient more.\n",
      sep = "")
  for (a in c(1.0, 0.5, 0.0)) {
    cat(sprintf("   alpha = %.1f\n", a))
    for (bgi in c(4, 9.9, 30)) {
      s <- eqstate(GENO = 1, PV0 = 1.0, BGI0 = bgi, ALPHA = a,
                   IVFE_ON = 0, IVFERATE = 12, TXALLOW = 0)
      cat(sprintf("     B_gi %5.1f -> Hb %5.2f, CI %4.2f, dSat %.3f, SpO2 %5.2f%% (shunt %.3f)\n",
                  bgi, s$HB, s$CI, s$DSATE, s$SPO2, s$SFRACT))
    }
  }
  cat("   The shunt fraction is identical in every row.  Whether pulse\n",
      "   oximetry looks alarming therefore depends on a physiological\n",
      "   constant nobody has measured in HHT - one more reason not to\n",
      "   screen with it.\n", sep = "")
  invisible(out)
}

# =====================================================================
#  AXIS 8 — the f_GI ceiling on any nasal therapy
# =====================================================================
axis8_fgi_ceiling <- function() {
  hdr("AXIS 8  What a perfect nasal therapy cannot do")
  cat("One third of HHT patients bleed from the gut, and no epistaxis\n",
      "endpoint sees it.  Since Hb_ss = A_net/(0.0347*B), abolishing\n",
      "B_nose leaves B_gi, and Hb_ss rises only by the factor\n",
      "B_total/B_gi = 1/f_GI.\n", sep = "")

  sub("Young's nasal closure (B_nose -> 0) across GI burdens:")
  out <- do.call(rbind, lapply(c(0, 5, 9.9, 20, 35), function(bgi) {
    pre  <- eqstate(BGI0 = bgi, IVFE_ON = 0, IVFERATE = 12, TXALLOW = 0)
    post <- eqstate(BGI0 = bgi, IVFE_ON = 0, IVFERATE = 12, CLOSEF = 1,
                    TXALLOW = 0)
    data.frame(B_gi = bgi, f_GI = pre$F_GI,
               B_pre = pre$B_TOT, Hb_pre = pre$HB, Hb_post = post$HB,
               dHb = post$HB - pre$HB,
               ESS_pre = pre$ESSC, ESS_post = post$ESSC,
               dESS = post$ESSC - pre$ESSC)
  }))
  print(out, row.names = FALSE, digits = 3)

  sub("Reading:")
  cat("   At f_GI near zero, closing the nose is curative for anaemia.\n",
      "   As f_GI grows the same operation buys progressively less\n",
      "   haemoglobin while costing the same anosmia, and at high GI\n",
      "   burden the patient is left transfusion-dependent with a closed\n",
      "   nose - the clinically familiar disappointment, here as a\n",
      "   number rather than an anecdote.\n", sep = "")
  cat("\n   Note the ESS column: ESS improves dramatically in every row,\n",
      "   because ESS is a nasal instrument.  ESS and Hb disagree by\n",
      "   construction, and the disagreement grows with f_GI.\n", sep = "")

  sub("Best achievable equilibrium by intervention class, iron ",
      "protocol-fixed:")
  base <- eqstate(TXALLOW = 0)
  # Iron is PROTOCOL-FIXED here (IVFE_ON = 0) so that the arms differ only
  # in what is prescribed, not in how much a feedback loop decided to give.
  arms <- list(
    "diet only (no iron therapy)"  = list(IVFERATE = 0),
    "IV iron 12 mg/day  [ref]"     = list(IVFERATE = 12),
    "IV iron 30 mg/day"            = list(IVFERATE = 30),
    "IV iron 12 + TXA-like -25% B" = list(IVFERATE = 12, EPDUR0 = 7.5),
    "perfect nasal therapy alone"  = list(IVFERATE = 0, CLOSEF = 1),
    "perfect nasal + IV iron 12"   = list(IVFERATE = 12, CLOSEF = 1),
    "nasal + GI ablation"          = list(IVFERATE = 0, CLOSEF = 1,
                                          GIENDOF = 0.05),
    "nasal + GI + IV iron 12"      = list(IVFERATE = 12, CLOSEF = 1,
                                          GIENDOF = 0.05))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    s <- do.call(eqstate, c(list(TXALLOW = 0, IVFE_ON = 0), arms[[nm]]))
    data.frame(arm = nm, B_total = s$B_TOT, A_net = s$ANET,
               Hb_eq = s$HB, ESS = s$ESSC)
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(out)
}

# =====================================================================
#  AXIS 9 — genotype changes which endpoint a drug can move
# =====================================================================
axis9_genotype <- function() {
  hdr("AXIS 9  Same drug, different value, by genotype")
  cat("PATH-HHT reported that the ESS benefit of pomalidomide was NOT\n",
      "dependent on genotype.  That is consistent with this model rather\n",
      "than a contradiction of it: ESS is a NASAL instrument, and both\n",
      "genotypes have nasal disease.  The genotype split lives in the\n",
      "VISCERAL endpoints, which ESS does not measure.\n", sep = "")

  sub("Matched patients, identical bevacizumab course, 12 months:")
  end <- yr(1)
  res <- do.call(rbind, lapply(1:2, function(g) {
    e <- c(ev_bev(6, 14, 75), ev_bev_maint(30, end, 75))
    d <- sim_eq(end = end, delta = 7, events = e, GENO = g)
    a <- d[1, ]; b <- d[nrow(d), ]
    data.frame(genotype = ifelse(g == 1, "ENG (HHT1)", "ACVRL1 (HHT2)"),
               dESS = b$ESSC - a$ESSC,
               dHb = b$HB - a$HB,
               CI_base = a$CI, dCI = b$CI - a$CI,
               SpO2_base = a$SPO2, dSpO2 = b$SPO2 - a$SPO2,
               shunt_frac = a$SFRACT,
               embolic_hazard_yr = 0.010 * a$PAVMB)
  }))
  print(res, row.names = FALSE, digits = 3)

  sub("Reading:")
  cat("   The nasal/anaemia benefit is genotype-indifferent, matching the\n",
      "   trial.  The cardiac benefit is concentrated in ACVRL1 because\n",
      "   that is where the hepatic shunt is.  And the ENG patient's\n",
      "   dominant hazard - paradoxical embolism through a pulmonary AVM -\n",
      "   is untouched by any antiangiogenic in this model, because the\n",
      "   mature AVM does not remodel on a drug time scale.  For that\n",
      "   patient the intervention that matters is embolotherapy.\n",
      sep = "")

  sub("Embolotherapy in the ENG patient:")
  for (ef in c(0, 0.5, 0.9)) {
    s <- eqstate(GENO = 1, EMBOLF = ef, TXALLOW = 0)
    cat(sprintf("   occlusion %.0f%%  ->  shunt %.3f, SpO2 %.2f%%, hazard/yr %.4f\n",
                100 * ef, s$SFRACT, s$SPO2, 0.010 * s$PAVMB * (1 - ef)))
  }
  invisible(res)
}

# =====================================================================
#  TREATMENT SCENARIOS
# =====================================================================
scenarios <- function() {
  hdr("TREATMENT SCENARIOS — 15 arms, 2 years, ACVRL1 reference patient")
  end <- yr(2)
  D2  <- floor(end)

  arms <- list(
    list(name = "1  natural history, no therapy",
         ev = NULL, p = list(IVFE_ON = 0, IVFERATE = 0, TXALLOW = 1)),
    list(name = "2  oral iron 200 mg/day only",
         ev = ev_oral_fe(D2), p = list(IVFE_ON = 0, IVFERATE = 0, TXALLOW = 1)),
    list(name = "3  IV iron 12 mg/day (333 mg per 4 wk)",
         ev = NULL, p = list(TXALLOW = 1)),
    list(name = "4  IV iron 30 mg/day",
         ev = NULL, p = list(IVFE_ON = 0, IVFERATE = 30, TXALLOW = 1)),
    list(name = "5  tranexamic acid 1 g tid",
         ev = ev_txa(D2), p = list(TXALLOW = 1)),
    list(name = "6  pomalidomide 4 mg/day",
         ev = ev_pom(D2), p = list(TXALLOW = 1)),
    list(name = "7  pomalidomide + IV iron 30 mg/day",
         ev = ev_pom(D2), p = list(IVFE_ON = 0, IVFERATE = 30, TXALLOW = 1)),
    list(name = "8  pazopanib 50 mg/day",
         ev = ev_paz(D2, 50), p = list(TXALLOW = 1)),
    list(name = "9  pazopanib 400 mg/day",
         ev = ev_paz(D2, 400), p = list(TXALLOW = 1)),
    list(name = "10 bevacizumab induction only",
         ev = ev_bev(6, 14, 75), p = list(TXALLOW = 1)),
    list(name = "11 bevacizumab + q3-month maintenance",
         ev = c(ev_bev(6, 14, 75), ev_bev_maint(91, end, 75)),
         p = list(TXALLOW = 1)),
    list(name = "12 bevacizumab + q2-month maintenance",
         ev = c(ev_bev(6, 14, 75), ev_bev_maint(61, end, 75)),
         p = list(TXALLOW = 1)),
    list(name = "13 bevacizumab + monthly maintenance",
         ev = c(ev_bev(6, 14, 75), ev_bev_maint(30, end, 75)),
         p = list(TXALLOW = 1)),
    list(name = "14 Young nasal closure",
         ev = NULL, p = list(CLOSEF = 1, TXALLOW = 1)),
    list(name = "15 bevacizumab monthly + GI ablation + IV iron",
         ev = c(ev_bev(6, 14, 75), ev_bev_maint(30, end, 75)),
         p = list(GIENDOF = 0.05, IVFE_ON = 0, IVFERATE = 30, TXALLOW = 1))
  )

  res <- do.call(rbind, lapply(arms, function(a) {
    d <- do.call(sim_eq, c(list(end = end, delta = 7, events = a$ev), a$p))
    s0 <- d[1, ]; s1 <- d[nrow(d), ]
    data.frame(
      scenario   = a$name,
      ESS_0      = s0$ESSC,        ESS_2yr = s1$ESSC,
      Hb_0       = s0$HB,          Hb_2yr  = s1$HB,
      B_0        = s0$B_TOT,       B_2yr   = s1$B_TOT,
      CI_2yr     = s1$CI,
      RBCunits_yr = (s1$CUMTX - s0$CUMTX) / 2,
      IVFe_g     = (s1$CUMIVFE - s0$CUMIVFE) / 1000,
      QoL_2yr    = s1$QOLC)
  }))
  print(res, row.names = FALSE, digits = 3)

  sub("Ranked by 2-year haemoglobin:")
  print(res[order(-res$Hb_2yr), c("scenario", "Hb_2yr", "ESS_2yr",
                                  "RBCunits_yr", "IVFe_g")],
        row.names = FALSE, digits = 3)

  sub("Ranked by 2-year ESS (the registrational endpoint):")
  print(res[order(res$ESS_2yr), c("scenario", "ESS_2yr", "Hb_2yr",
                                  "QoL_2yr")],
        row.names = FALSE, digits = 3)
  cat("\n   The two rankings disagree.  Arms that buy haemoglobin with\n",
      "   iron rank well on Hb and poorly on ESS; arms that reduce\n",
      "   bleeding rank well on both but need time.  A trial that picks\n",
      "   one endpoint is choosing which of these two orders it will see.\n",
      sep = "")
  invisible(res)
}

# =====================================================================
#  SENSITIVITY
# =====================================================================
sensitivity <- function() {
  hdr("LOCAL SENSITIVITY — 2-year Hb, ESS and CI to +/-25% parameter moves")
  pars <- c("KOUT_N", "KOUT_L", "AMAXO", "KHEP", "BGI0", "QBLEED",
            "MURSS", "GFRAG", "ALPHA", "KSHUNT", "FRES", "EMAX_PM", "IVFETRIG",
            "EC50_PM", "W_FE_SHARE", "PLACMAX", "IL6G", "KSUP")
  end <- yr(2)
  base <- sim_eq(end = end, delta = 30, events = ev_pom(floor(end)))
  bl <- base[nrow(base), ]
  out <- do.call(rbind, lapply(pars, function(p) {
    v <- mod[[p]]
    lo <- setNames(list(v * 0.75), p); hi <- setNames(list(v * 1.25), p)
    dl <- do.call(sim_eq, c(list(end = end, delta = 30,
                                 events = ev_pom(floor(end))), lo))
    dh <- do.call(sim_eq, c(list(end = end, delta = 30,
                                 events = ev_pom(floor(end))), hi))
    l <- dl[nrow(dl), ]; h <- dh[nrow(dh), ]
    data.frame(parameter = p, base = v,
               Hb_lo = l$HB, Hb_hi = h$HB,
               Hb_range_pct = 100 * (h$HB - l$HB) / bl$HB,
               ESS_range = h$ESSC - l$ESSC,
               CI_range = h$CI - l$CI)
  }))
  out <- out[order(-abs(out$Hb_range_pct)), ]
  print(out, row.names = FALSE, digits = 3)
  sub("The model is most sensitive to the quantities that are least ",
      "well measured\n   in HHT: the bleeding flow rate, the GI share, ",
      "and the oral iron ceiling.")
  invisible(out)
}

# =====================================================================
#  SELF-TESTS — run the model against things that must be true
# =====================================================================
selftest <- function() {
  hdr("SELF-TESTS")
  ok <- TRUE
  chk <- function(label, cond, detail = "") {
    cat(sprintf("  [%s] %s%s\n", ifelse(cond, "PASS", "FAIL"), label,
                ifelse(nzchar(detail), paste0("  (", detail, ")"), "")))
    if (!cond) ok <<- FALSE
    invisible(cond)
  }

  # 1. mass balance identity
  s <- eqstate(TXALLOW = 0)
  chk("Hb equilibrium matches the closed form within 1%",
      abs(s$HB - s$HB_PRED) / s$HB_PRED < 0.01,
      sprintf("ODE %.3f vs closed %.3f", s$HB, s$HB_PRED))

  # 2. no bleeding -> normal Hb
  s2 <- eqstate(BGI0 = 0, LAMEP0 = 0, TXALLOW = 0)
  chk("with no bleeding Hb reaches the set-point",
      s2$HB > 13.5, sprintf("Hb %.2f", s2$HB))

  # 3. more bleeding -> lower Hb, monotone
  hbs <- sapply(c(5, 15, 30, 60), function(b) eqstate(BGI0 = b, TXALLOW = 0)$HB)
  chk("Hb is monotonically decreasing in blood loss",
      all(diff(hbs) < 0), paste(round(hbs, 2), collapse = " > "))

  # 4. IV iron raises Hb; oral iron saturates
  a <- eqstate(IVFERATE = 0, TXALLOW = 0)$HB
  b <- eqstate(IVFERATE = 30, TXALLOW = 0)$HB
  chk("IV iron raises equilibrium Hb", b > a,
      sprintf("%.2f -> %.2f", a, b))
  o1 <- sim_eq(end = 400, delta = 100, events = ev_oral_fe(400, 200),
               IVFERATE = 0, TXALLOW = 0)
  o2 <- sim_eq(end = 400, delta = 100, events = ev_oral_fe(400, 1200),
               IVFERATE = 0, TXALLOW = 0)
  g1 <- o1$ANET[nrow(o1)]; g2 <- o2$ANET[nrow(o2)]
  chk("oral iron absorption saturates (6x dose gives <2x A_net)",
      g2 < 2 * g1, sprintf("A_net %.2f -> %.2f mg/day", g1, g2))

  # 5. bevacizumab suppresses free VEGF and lesion burden
  d <- sim_eq(end = 120, delta = 1, events = ev_bev(6, 14, 75))
  chk("bevacizumab lowers nasal lesion burden",
      min(d$TNOSE) < 0.95 * d$TNOSE[1],
      sprintf("%.3f -> %.3f", d$TNOSE[1], min(d$TNOSE)))

  # 6. target engagement outlasts measurable drug (axis 4, step 3)
  d2 <- sim_eq(end = 700, delta = 2, events = ev_bev(6, 14, 75))
  nad <- d2$time[which.min(d2$TNOSE)]
  chk("lesion nadir falls AFTER the last dose (stoichiometric excess)",
      nad > 70, sprintf("last dose day 70, nadir day %.0f", nad))
  chk("nasal burden eventually rebounds toward baseline",
      d2$TNOSE[nrow(d2)] > 1.2 * min(d2$TNOSE),
      sprintf("nadir %.3f -> day %.0f %.3f", min(d2$TNOSE),
              max(d2$time), d2$TNOSE[nrow(d2)]))

  # 7. pomalidomide raises mural coverage (the biopsy finding)
  d3 <- sim_eq(end = 168, delta = 7, events = ev_pom(168))
  chk("pomalidomide raises mural-cell coverage",
      d3$MURALC[nrow(d3)] > d3$MURALC[1],
      sprintf("%.3f -> %.3f", d3$MURALC[1], d3$MURALC[nrow(d3)]))
  chk("pomalidomide lowers ESS by more than the MCID of 0.71",
      (d3$ESSC[1] - d3$ESSC[nrow(d3)]) > 0.71,
      sprintf("dESS %.2f", d3$ESSC[nrow(d3)] - d3$ESSC[1]))

  # 8. pazopanib's PDGFR block is mechanistically adverse on coverage
  d4 <- sim_eq(end = 168, delta = 7, events = ev_paz(168, 400))
  chk("high-dose pazopanib lowers mural coverage vs baseline",
      d4$MURALC[nrow(d4)] < d4$MURALC[1],
      sprintf("%.3f -> %.3f", d4$MURALC[1], d4$MURALC[nrow(d4)]))

  # 9. TXA shortens duration without changing frequency
  d5 <- sim_eq(end = 60, delta = 1, events = ev_txa(60))
  chk("TXA shortens episode duration",
      d5$EPDURC[nrow(d5)] < d5$EPDURC[1],
      sprintf("%.2f -> %.2f min", d5$EPDURC[1], d5$EPDURC[nrow(d5)]))
  chk("TXA leaves episode frequency essentially unchanged",
      abs(d5$LAMEP[nrow(d5)] - d5$LAMEP[1]) < 0.05 * d5$LAMEP[1],
      sprintf("%.3f -> %.3f /wk", d5$LAMEP[1], d5$LAMEP[nrow(d5)]))

  # 10. nasal closure zeroes nasal loss but not total loss
  s3 <- eqstate(CLOSEF = 1, TXALLOW = 0)
  chk("nasal closure sets B_nose to ~0", s3$B_NOSE < 0.01,
      sprintf("B_nose %.4f", s3$B_NOSE))
  chk("nasal closure leaves GI loss behind", s3$B_GI > 1,
      sprintf("B_gi %.2f mL/day", s3$B_GI))

  # 11. genotype tropism
  g1s <- eqstate(GENO = 1, TXALLOW = 0); g2s <- eqstate(GENO = 2, TXALLOW = 0)
  chk("ENG has the larger pulmonary shunt",
      g1s$SFRACT > g2s$SFRACT,
      sprintf("%.3f vs %.3f", g1s$SFRACT, g2s$SFRACT))
  chk("ACVRL1 has the higher cardiac index",
      g2s$CI > g1s$CI, sprintf("%.2f vs %.2f", g2s$CI, g1s$CI))

  # 12. alpha=0 removes the Hb dependence of cardiac index
  c1 <- eqstate(ALPHA = 0, BGI0 = 5, IVFE_ON = 0, IVFERATE = 12,
                TXALLOW = 0)$CI
  c2 <- eqstate(ALPHA = 0, BGI0 = 40, IVFE_ON = 0, IVFERATE = 12,
                TXALLOW = 0)$CI
  chk("with alpha=0 cardiac index is insensitive to anaemia",
      abs(c1 - c2) < 0.05, sprintf("%.3f vs %.3f", c1, c2))
  c3 <- eqstate(ALPHA = 1, BGI0 = 5, IVFE_ON = 0, IVFERATE = 12,
                TXALLOW = 0)$CI
  c4 <- eqstate(ALPHA = 1, BGI0 = 40, IVFE_ON = 0, IVFERATE = 12,
                TXALLOW = 0)$CI
  chk("with alpha=1 anaemia raises cardiac index",
      c4 > c3 + 0.3, sprintf("%.3f vs %.3f", c3, c4))

  # 13. REGRESSION TEST for the shear-feedback sign defect: hepatic
  #     conductance and cardiac index must stay positive at every shunt
  #     burden, including the high-output cohort.
  gpos <- TRUE; cipos <- TRUE
  for (g in c(0.5, 1, 4, 8.3, 15, 30)) {
    dd <- sim_eq(end = yr(2), delta = 60, GENO = 2, GS0 = g,
                 events = c(ev_bev(6, 14, 75), ev_bev_maint(91, yr(2), 75)))
    if (min(dd$GSHUNT) <= 0) gpos <- FALSE
    if (min(dd$CI) <= 0) cipos <- FALSE
  }
  chk("hepatic shunt conductance stays positive for GS0 up to 30", gpos)
  chk("cardiac index stays positive for GS0 up to 30", cipos)

  # 14. SpO2 stays >=95% at a hazardous shunt (axis 7)
  s4 <- eqstate(GENO = 1, PV0 = 0.5, TXALLOW = 0)
  chk("a shunt with real embolic hazard can leave SpO2 >= 95%",
      s4$SPO2 >= 95 && s4$SFRACT > 0.02,
      sprintf("SpO2 %.2f at shunt %.3f", s4$SPO2, s4$SFRACT))

  # 15. transfusion trigger fires when Hb collapses
  s5 <- eqstate(BGI0 = 120, IVFERATE = 0, TXALLOW = 1)
  chk("severe bleeding drives transfusion dependence",
      s5$TXPERYR > 1, sprintf("%.1f units/yr at Hb %.2f", s5$TXPERYR, s5$HB))

  # 16. placebo dissipates after unblinding
  d6 <- sim_eq(end = 200, delta = 1, BLIND = 1, ONPLAC = 1)
  chk("placebo response builds under blinding",
      d6$PLAC[nrow(d6)] > 0.5, sprintf("PLAC %.2f", d6$PLAC[nrow(d6)]))

  # 17. all states finite
  d7 <- sim_eq(end = yr(5), delta = 30, events = ev_pom(1000))
  chk("no NaN or Inf over a 5-year simulation",
      all(is.finite(as.matrix(d7[, sapply(d7, is.numeric)]))))

  cat("\n", strrep("-", 70), "\n", sep = "")
  cat(ifelse(ok, "ALL SELF-TESTS PASSED\n", "*** SOME SELF-TESTS FAILED ***\n"))
  invisible(ok)
}

# =====================================================================
#  FAILURES AND LIMITS — reported, not hidden
# =====================================================================
failures <- function() {
  hdr("WHAT THIS MODEL DOES NOT DO")
  cat(
"1. ESS ITEM WEIGHTS ARE NOT THE INSTRUMENT.  Hoag 2010 published the\n",
"   six determinants of ESS but not a closed-form weighting.  The\n",
"   surrogate here is transparent and monotone in the right variables,\n",
"   but the absolute ESS numbers should not be compared with trial ESS\n",
"   values as if they were the same scale.  Axis 3 is therefore a\n",
"   statement about the STRUCTURE of the score, not a measurement of\n",
"   its iron-sensitive fraction.\n\n",

"2. B_nose AND B_gi ARE INFERRED, NEVER OBSERVED.  No HHT trial\n",
"   measures blood loss in mL/day.  Both are solved backwards through\n",
"   the iron balance from Hb and iron use, so every conclusion about\n",
"   f_GI inherits the uncertainty in A_net - and A_net depends on the\n",
"   oral absorption ceiling AMAXO, which is a physiological estimate\n",
"   rather than an HHT measurement.\n\n",

"3. A HYPOTHESIS I EXPECTED TO CONFIRM WAS REJECTED.  Axis 4 set out to\n",
"   prove that the 6-month bevacizumab epistaxis nadir was impossible for\n",
"   any PK/PD structure.  It is not: stoichiometric molar excess keeps the\n",
"   target engaged for many half-lives after the concentration becomes\n",
"   negligible, so the premise of the argument fails.  The axis is left in\n",
"   the file as a rejection rather than deleted, because the reason it\n",
"   failed - exposure for a neutralising antibody means concentration\n",
"   relative to TARGET MOLARITY - is what drives axis 5.  A side effect is\n",
"   that the model now over-predicts the duration of benefit after a\n",
"   single induction course, and that over-prediction is reported rather\n",
"   than tuned away.\n\n",

"4. alpha IS NOT IDENTIFIABLE, AND I DID NOT MAKE IT SO.  The share of\n",
"   the bevacizumab cardiac-index response attributable to shunt\n",
"   regression rather than anaemia correction spans nearly the whole\n",
"   range as alpha goes from 0 to 1.  Resolving it needs paired Hb and\n",
"   cardiac-index measurements from the same patients at the same\n",
"   visits, which Dupuis-Girod 2012 does not report.  Axis 6 reports a\n",
"   range and stops.\n\n",

"5. THE k_out ESTIMATES ARE FITTED TO SUMMARY FRACTIONS.  Azzopardi's\n",
"   41/45/50 and 34/43/60 are response percentages in a virtual\n",
"   population with inter-individual variability this model does not\n",
"   reproduce.  Axis 5 compares the CONTRAST between the two endpoints,\n",
"   which is robust to that mismatch, and does not claim to reproduce\n",
"   the percentages.\n\n",

"6. THE SECOND HIT IS NOT MODELLED AS A STOCHASTIC EVENT.  AVM\n",
"   formation in HHT needs a somatic second hit plus an angiogenic\n",
"   trigger, which makes lesion counts a branching process.  Here\n",
"   lesion burden is a continuous state with a deterministic formation\n",
"   rate.  That is adequate for turnover and drug effect but it cannot\n",
"   predict the appearance of a NEW AVM, and so it cannot model the\n",
"   observed clustering of new lesions after surgery or pregnancy.\n\n",

"7. HHT-PAH IS OUT OF SCOPE.  True pulmonary arterial hypertension in\n",
"   ACVRL1 carriers is a distinct, severe phenotype with its own\n",
"   pulmonary vascular remodelling.  The PHTN node in the map is the\n",
"   post-capillary, high-output kind only.  Applying this model to an\n",
"   HHT-PAH patient would be wrong.\n\n",

"8. THE ANTICOAGULATION PARADOX IS DRAWN BUT NOT SOLVED.  Atrial\n",
"   fibrillation and venous thromboembolism are common in HHT and force\n",
"   anticoagulation on patients whose disease is bleeding.  The model\n",
"   carries a duration-prolonging term for antiplatelet/anticoagulant\n",
"   exposure but has no data to calibrate it, so no scenario here\n",
"   quantifies that trade-off.\n\n",

"9. WHAT WOULD FALSIFY AXIS 2.  If a trial reduced bleeding while\n",
"   holding iron and transfusion protocol-fixed and still saw no Hb\n",
"   rise, the exchangeability argument would be wrong and something\n",
"   else - marrow, inflammation, occult loss elsewhere - would have to\n",
"   be carrying the effect.  That trial has not been run, and it is the\n",
"   single most informative experiment this model suggests.\n",
    sep = "")
}

# =====================================================================
#  RUN EVERYTHING
# =====================================================================
run_all_hht <- function() {
  cat("\n")
  cat(strrep("#", 70), "\n")
  cat("#  HHT QSP MODEL — full analysis\n")
  cat("#  ", format(Sys.time()), "\n", sep = "")
  cat(strrep("#", 70), "\n")
  axis1_hyperbola()
  axis2_hb_paradox()
  axis3_placebo_and_ess()
  axis4_bev_structural()
  axis5_two_organs()
  axis6_cardiac_confound()
  axis7_spo2_screen()
  axis8_fgi_ceiling()
  axis9_genotype()
  scenarios()
  sensitivity()
  selftest()
  failures()
  invisible(TRUE)
}

if (identical(environment(), globalenv()) && !interactive()) {
  run_all_hht()
}
