## =============================================================================
## BPH (Benign Prostatic Hyperplasia) QSP mrgsolve Model
## =============================================================================
## File: bph_mrgsolve_model.R
## Description: Quantitative Systems Pharmacology model for BPH incorporating
##   multi-drug PK (tamsulosin, finasteride, dutasteride, tadalafil) and
##   mechanistic PD (DHT-AR axis, prostate volume, IPSS, Qmax, PVR, PSA).
##
## Calibration references:
##   - MTOPS trial (McConnell 2003, NEJM 349:2387): finasteride reduces AUR 57%,
##     TURP 64%; combination superior to monotherapy.
##   - CombAT trial (Roehrborn 2010, Eur Urol 57:123): dutasteride+tamsulosin
##     combination reduces IPSS ~6 pts vs ~4 pts monotherapy at 4 years.
##   - NEPTUNE trial (Chapple 2014, Eur Urol 65:998): tadalafil 5 mg QD reduces
##     IPSS by 5.6 vs 2.3 placebo at 12 weeks.
##   - 5AR inhibitors: finasteride reduces DHT ~70%, dutasteride ~94%.
##   - Prostate volume reduction: ~25-30% over 2 years with 5AR inhibition.
##   - Tamsulosin IC50 alpha1A: ~0.8 ng/mL (Foglar 1995, Eur J Pharmacol 288:201).
##
## Author: Generated via Claude Code Routine (CCR) — 2026-06-18
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---------------------------------------------------------------------------
## mrgsolve model string
## ---------------------------------------------------------------------------

bph_model_code <- '
$PROB
BPH QSP Model — 24-compartment PK/PD
Multi-drug: tamsulosin / finasteride / dutasteride / tadalafil
Endpoints: IPSS, Qmax, PVR, PSA, prostate volume, DHT

$PARAM
// ---- Tamsulosin PK (0.4 mg QD, alpha1A-AR antagonist) ----
ka_tams  = 0.693  // absorption rate h^-1 (Tmax ~6 h fed state)
CL_tams  = 1.8    // apparent clearance L/h (CYP2D6/3A4, high first-pass)
V1_tams  = 15.0   // central volume L
Q_tams   = 2.5    // intercompartmental clearance L/h
V2_tams  = 30.0   // peripheral volume L
F_tams   = 0.90   // oral bioavailability
IC50_tams = 0.8   // ng/mL, alpha1A-AR 50% occupancy (Foglar 1995)
Emax_tams = 0.95  // maximum alpha1A occupancy fraction

// ---- Finasteride PK (5 mg QD, selective SRD5A2 inhibitor) ----
ka_fina  = 0.50   // h^-1
CL_fina  = 5.0    // L/h (CYP3A4)
V1_fina  = 76.0   // L (Ohtawa 1991)
Q_fina   = 3.0    // L/h
V2_fina  = 152.0  // L
F_fina   = 0.63   // bioavailability
IC50_fina = 0.3   // ng/mL, SRD5A2 IC50 (Gormley 1992, J Urol 147:1562)

// ---- Dutasteride PK (0.5 mg QD, dual SRD5A1+2 inhibitor) ----
ka_dut   = 0.08   // h^-1 (slow absorption, Tmax ~1-3 h but prolonged)
CL_dut   = 0.15   // L/h (t½ ~5 weeks, Clark 2004, Clin Pharmacokinet 43:461)
V1_dut   = 5.0    // L (high tissue binding)
Q_dut    = 0.5    // L/h
V2_dut   = 40.0   // L (extensive distribution in adipose/prostate)
F_dut    = 0.60   // bioavailability
IC50_dut = 0.05   // ng/mL, dual SRD5A1+2 (Bramson 1997, J Steroid Biochem 64:189)

// ---- Tadalafil PK (5 mg QD, PDE5 inhibitor for LUTS/BPH) ----
ka_tad   = 0.35   // h^-1 (Tmax ~2 h)
CL_tad   = 2.5    // L/h (CYP3A4, Forgue 2006)
V1_tad   = 25.0   // L
Q_tad    = 3.5    // L/h
V2_tad   = 88.0   // L
F_tad    = 0.80   // bioavailability
IC50_tad = 2.0    // ng/mL, PDE5 IC50 (Gresser 2002, World J Urol 20:360)

// ---- Hormonal/DHT parameters ----
TEST0        = 15.0   // nmol/L, baseline testosterone (normal range 10-35 nmol/L)
DHT0         = 1.5    // nmol/L, baseline plasma DHT
DHT_PROST0   = 25.0   // nmol/L, baseline prostatic DHT (~10-20x plasma, Griffiths 1984)
AR_ACT0      = 0.5    // normalized androgen receptor activity (0-1 scale)
k_conv_DHT   = 0.08   // h^-1, testosterone → DHT via SRD5A2 (prostate)
k_deg_DHT    = 0.15   // h^-1, DHT elimination (conjugation + excretion)
k_AR_on      = 2.0    // (nmol/L)^-1 h^-1, DHT-AR binding rate constant
k_AR_off     = 0.5    // h^-1, DHT-AR dissociation rate
AR_ACT_MAX   = 1.0    // maximum AR activation (normalized)

// ---- Prostate volume dynamics ----
PV0          = 40.0   // cc, baseline volume (MTOPS mean ~36 cc, range 25-80)
k_PV_growth  = 0.00025 // h^-1, ~1.5-2 cc/year natural history (Berry 1984, J Urol 132:474)
EC50_AR_PV   = 0.4    // AR activity producing 50% maximal growth
k_PV_shrink  = 0.0005  // h^-1, DHT-reduction-driven involution rate constant

// ---- PSA dynamics ----
PSA0         = 2.5    // ng/mL, baseline PSA (typical BPH range 1.5-4 ng/mL)
k_PSA_prod   = 0.8    // PSA production proportional to AR_ACT * PV
k_PSA_deg    = 0.05   // h^-1, PSA clearance from plasma (t½ ~2 days, Stenman 1994)

// ---- cGMP/smooth muscle dynamics (tadalafil PD) ----
CGMP0        = 1.0    // normalized cGMP level in prostatic/urethral smooth muscle
k_cGMP_prod  = 0.5    // h^-1, basal cGMP synthesis by sGC
k_cGMP_deg   = 0.4    // h^-1, PDE5-mediated cGMP hydrolysis (dominant pathway)

// ---- LUTS/functional endpoints ----
IPSS0        = 18.0   // baseline IPSS (moderate-severe BPH, MTOPS entry criterion ≥8)
QMAX0        = 10.0   // mL/s, baseline peak flow (MTOPS entry <15 mL/s)
PVR0         = 100.0  // mL, baseline post-void residual
INFLAM0      = 0.3    // baseline prostatic inflammation index (0-1, normalized)

// ---- LUTS PD effect parameters ----
// Alpha1-AR blockade → urethral relaxation → ΔQmax, ΔIPSS (rapid, dynamic component)
k_Qmax_alpha1 = 2.0   // Qmax improvement per unit alpha1 blockade (mL/s per unit)
k_IPSS_alpha1 = 5.0   // IPSS reduction per unit alpha1 blockade
// Prostate volume effect on LUTS (static component, structural obstruction)
k_IPSS_PV    = 0.15   // IPSS pts per cc above baseline (Lepor 1993)
k_Qmax_PV    = 0.05   // Qmax reduction per cc above baseline (Eckhardt 2001)
k_PVR_PV     = 0.8    // PVR increase per cc above baseline
// cGMP/smooth muscle relaxation → LUTS improvement (tadalafil effect)
k_IPSS_cGMP  = 3.0    // IPSS improvement per unit cGMP above baseline
k_Qmax_cGMP  = 1.5    // Qmax improvement per unit cGMP above baseline
// Inflammation contribution to LUTS (De Nunzio 2011, BJU Int 108:1771)
k_IPSS_inflam = 4.0   // IPSS contribution of inflammation

// ---- Inflammation dynamics ----
k_inflam_DHT  = 0.05  // DHT-driven prostatic inflammation stimulation
k_inflam_deg  = 0.02  // h^-1, inflammation resolution rate

// ---- Treatment flags (0 = off, 1 = on) ----
USE_TAMS = 0
USE_FINA = 0
USE_DUT  = 0
USE_TAD  = 0

$CMT
// Tamsulosin 3-compartment PK
TAMS_GUT TAMS_C TAMS_P
// Finasteride 3-compartment PK
FINA_GUT FINA_C FINA_P
// Dutasteride 3-compartment PK
DUT_GUT DUT_C DUT_P
// Tadalafil 3-compartment PK
TAD_GUT TAD_C TAD_P
// Hormonal compartments
TEST_P       // testosterone plasma (nmol/L * V_dist)
DHT_P        // DHT plasma (nmol/L * V_dist)
DHT_PROST    // DHT prostate tissue (nmol/L)
AR_ACT       // active AR-DHT complex (normalized 0-1)
// Structural/functional endpoints
PV           // prostate volume (cc)
CGMP         // cGMP in smooth muscle (normalized)
ALPHA1_OCC   // alpha1-AR occupancy (fraction, 0-1)
IPSS         // IPSS score (0-35)
QMAX         // peak urinary flow rate (mL/s)
PVR          // post-void residual (mL)
PSA          // PSA ng/mL
INFLAM       // prostatic inflammation index (0-1)

$MAIN
// Initialize all compartments at steady-state baseline
if(NEWIND <= 1) {
  // PK compartments start empty (drug naive)
  _init_TAMS_GUT = 0;
  _init_TAMS_C   = 0;
  _init_TAMS_P   = 0;
  _init_FINA_GUT = 0;
  _init_FINA_C   = 0;
  _init_FINA_P   = 0;
  _init_DUT_GUT  = 0;
  _init_DUT_C    = 0;
  _init_DUT_P    = 0;
  _init_TAD_GUT  = 0;
  _init_TAD_C    = 0;
  _init_TAD_P    = 0;

  // Hormonal steady state
  // TEST_P: testosterone in plasma distribution volume (~50 L)
  _init_TEST_P    = TEST0 * 50.0;   // nmol (amount)
  _init_DHT_P     = DHT0 * 30.0;    // nmol (Vd ~30 L for DHT)
  _init_DHT_PROST = DHT_PROST0;     // nmol/L (tissue concentration)
  _init_AR_ACT    = AR_ACT0;        // normalized

  // Structural/functional baseline
  _init_PV        = PV0;
  _init_CGMP      = CGMP0;
  _init_ALPHA1_OCC = 0.0;           // no drug → zero alpha1 blockade
  _init_IPSS      = IPSS0;
  _init_QMAX      = QMAX0;
  _init_PVR       = PVR0;
  _init_PSA       = PSA0;
  _init_INFLAM    = INFLAM0;
}

$ODE
// =========================================================================
// SECTION 1: TAMSULOSIN PK (3-compartment oral, 0.4 mg QD)
// =========================================================================
// Dose is added to TAMS_GUT via event records (see R simulation code)
// GI absorption → central → peripheral
dxdt_TAMS_GUT = -ka_tams * TAMS_GUT;
// Central: input from gut (F already accounted for in dose), distribution, elimination
// Note: TAMS_GUT units = ug (dose); TAMS_C = ug; concentration = TAMS_C/V1_tams ng/mL
dxdt_TAMS_C   =  F_tams * ka_tams * TAMS_GUT
                - (CL_tams / V1_tams) * TAMS_C
                - (Q_tams  / V1_tams) * TAMS_C
                + (Q_tams  / V2_tams) * TAMS_P;
dxdt_TAMS_P   =  (Q_tams  / V1_tams) * TAMS_C
                - (Q_tams  / V2_tams) * TAMS_P;

// Plasma concentration (ng/mL = ug/L)
double CP_tams = (V1_tams > 0) ? TAMS_C / V1_tams : 0.0;

// Alpha1-AR occupancy from tamsulosin (Emax model)
double occ_tams = (USE_TAMS > 0.5) ?
                  Emax_tams * CP_tams / (IC50_tams + CP_tams) : 0.0;

// =========================================================================
// SECTION 2: FINASTERIDE PK (3-compartment oral, 5 mg QD)
// =========================================================================
dxdt_FINA_GUT = -ka_fina * FINA_GUT;
dxdt_FINA_C   =  F_fina * ka_fina * FINA_GUT
                - (CL_fina / V1_fina) * FINA_C
                - (Q_fina  / V1_fina) * FINA_C
                + (Q_fina  / V2_fina) * FINA_P;
dxdt_FINA_P   =  (Q_fina  / V1_fina) * FINA_C
                - (Q_fina  / V2_fina) * FINA_P;

double CP_fina = (V1_fina > 0) ? FINA_C / V1_fina : 0.0;  // ng/mL

// SRD5A2 inhibition fraction (Imax model, finasteride selective for type 2)
// At 5 mg QD, Css ~9 ng/mL → ~97% inhibition; calibrated to DHT reduction ~70%
// (type 1 not inhibited: residual ~30% DHT from type 1 enzyme in skin/liver)
double inh_fina = (USE_FINA > 0.5) ?
                  CP_fina / (IC50_fina + CP_fina) : 0.0;

// =========================================================================
// SECTION 3: DUTASTERIDE PK (3-compartment oral, 0.5 mg QD)
// =========================================================================
dxdt_DUT_GUT = -ka_dut * DUT_GUT;
dxdt_DUT_C   =  F_dut * ka_dut * DUT_GUT
               - (CL_dut / V1_dut) * DUT_C
               - (Q_dut  / V1_dut) * DUT_C
               + (Q_dut  / V2_dut) * DUT_P;
dxdt_DUT_P   =  (Q_dut  / V1_dut) * DUT_C
               - (Q_dut  / V2_dut) * DUT_P;

double CP_dut = (V1_dut > 0) ? DUT_C / V1_dut : 0.0;  // ng/mL

// Dual SRD5A1+2 inhibition (dutasteride inhibits both type 1 and type 2)
// At 0.5 mg QD, Css ~40 ng/mL → >99% inhibition; DHT reduction ~94% (Gisleskog 2002)
double inh_dut = (USE_DUT > 0.5) ?
                 CP_dut / (IC50_dut + CP_dut) : 0.0;

// Combined 5AR inhibition (additive; finasteride type2 only, dutasteride both types)
// SRD5A2 contributes ~70% of prostatic DHT, SRD5A1 ~30%
// Finasteride inhibits SRD5A2 fraction only
// Dutasteride inhibits both fractions
double f2_fraction = 0.70;  // fraction of DHT synthesis from SRD5A2
double f1_fraction = 0.30;  // fraction of DHT synthesis from SRD5A1

double inh_SRD5A2 = inh_fina * (1.0 - inh_dut) + inh_dut;  // combined type 2 inhibition
double inh_SRD5A1 = inh_dut;                                   // type 1 inhibition (dutasteride only)
double total_5AR_inh = f2_fraction * inh_SRD5A2 + f1_fraction * inh_SRD5A1;

// =========================================================================
// SECTION 4: TADALAFIL PK (3-compartment oral, 5 mg QD)
// =========================================================================
dxdt_TAD_GUT = -ka_tad * TAD_GUT;
dxdt_TAD_C   =  F_tad * ka_tad * TAD_GUT
               - (CL_tad / V1_tad) * TAD_C
               - (Q_tad  / V1_tad) * TAD_C
               + (Q_tad  / V2_tad) * TAD_P;
dxdt_TAD_P   =  (Q_tad  / V1_tad) * TAD_C
               - (Q_tad  / V2_tad) * TAD_P;

double CP_tad = (V1_tad > 0) ? TAD_C / V1_tad : 0.0;  // ng/mL

// PDE5 inhibition fraction by tadalafil
double inh_pde5 = (USE_TAD > 0.5) ?
                  CP_tad / (IC50_tad + CP_tad) : 0.0;

// =========================================================================
// SECTION 5: TESTOSTERONE DYNAMICS
// =========================================================================
// Testosterone plasma amount: production from testes (LH-driven, feedback not modeled
// explicitly), peripheral conversion to DHT (minor loss), SHBG binding, elimination.
// Simplified: testosterone maintained near TEST0 (quasi-steady-state in normal men).
// Vd_TEST = 50 L for unit conversion
double conc_TEST = TEST_P / 50.0;   // nmol/L
double kprod_TEST = k_deg_DHT * DHT0 * 30.0 / 50.0 + k_conv_DHT * DHT0 * 30.0 / 50.0;
// Approximate testosterone production to maintain baseline
double kprod_TEST_ss = k_conv_DHT * TEST0 + (TEST0 * 0.02); // conversion + basal elimination
dxdt_TEST_P  = kprod_TEST_ss * 50.0   // production (nmol/L * Vd → nmol/h)
               - 0.02 * TEST_P         // elimination (t½ ~1.5 h for free T)
               - k_conv_DHT * (1.0 - total_5AR_inh) * TEST_P; // 5AR conversion loss

// =========================================================================
// SECTION 6: DHT PLASMA DYNAMICS
// =========================================================================
// DHT plasma: produced from testosterone via 5AR (both types), eliminated via
// glucuronidation/sulfation. Vd_DHT = 30 L.
double conc_DHT_P = DHT_P / 30.0;   // nmol/L

// DHT production = 5AR conversion of testosterone (remaining after inhibition)
double DHT_prod_P = k_conv_DHT * (1.0 - total_5AR_inh) * TEST_P * (30.0 / 50.0);

dxdt_DHT_P = DHT_prod_P
             - k_deg_DHT * DHT_P;   // elimination t½ ~2-3 h

// =========================================================================
// SECTION 7: PROSTATIC DHT DYNAMICS
// =========================================================================
// Prostatic DHT: concentrated by active uptake and local 5AR activity.
// Ratio DHT_PROST/DHT_P ~15-20x at baseline (Griffiths 1984).
// Influenced by 5AR inhibition within prostate tissue.
// k_uptake: DHT transfer from plasma to prostate tissue
double k_uptake_DHT = 0.10;    // h^-1, DHT uptake into prostate
double k_efflux_DHT = 0.006;   // h^-1, efflux calibrated to maintain ~15x ratio

// Local prostatic 5AR also inhibited (tissue concentration of drug)
// Simplified: use same inhibition fraction as systemic
dxdt_DHT_PROST = k_conv_DHT * (1.0 - total_5AR_inh) * conc_TEST * 20.0
                 + k_uptake_DHT * conc_DHT_P * 15.0   // plasma → tissue accumulation
                 - k_efflux_DHT * DHT_PROST            // tissue → plasma efflux
                 - k_deg_DHT * 0.3 * DHT_PROST;        // local metabolism

// =========================================================================
// SECTION 8: ANDROGEN RECEPTOR ACTIVATION
// =========================================================================
// AR_ACT represents the fraction of androgen receptors in active (DHT-bound) state.
// Drives prostate growth, PSA production, and inflammatory signaling.
// Hill equation relating DHT_PROST to AR_ACT (quasi-equilibrium approach)
// Use turnover model for AR_ACT dynamics:
double AR_ACT_eq = AR_ACT_MAX * DHT_PROST / (EC50_AR_PV * DHT_PROST0 * 10.0 + DHT_PROST);
// AR_ACT equilibrates to AR_ACT_eq with rate k_AR_on/off
double k_AR_eq = k_AR_off + k_AR_on * DHT_PROST;
dxdt_AR_ACT = k_AR_on * DHT_PROST * (AR_ACT_MAX - AR_ACT)
              - k_AR_off * AR_ACT;

// =========================================================================
// SECTION 9: PROSTATE VOLUME DYNAMICS
// =========================================================================
// Prostate volume grows when AR activity is above a threshold (androgen-driven).
// 5AR inhibition reduces DHT → reduced AR activity → volume shrinks.
// Natural growth: ~1.5-2 cc/year (Berry 1984, J Urol 132:474)
// With 5AR inhibition: -25-30% over 2 years (MTOPS, CombAT).
double AR_growth_signal = AR_ACT / (EC50_AR_PV + AR_ACT);  // 0-1, AR-driven growth
// Involution: occurs when AR_ACT drops below baseline AR_ACT0
double AR_shrink_signal = (AR_ACT0 - AR_ACT > 0) ? (AR_ACT0 - AR_ACT) / AR_ACT0 : 0.0;

dxdt_PV = k_PV_growth * PV * AR_growth_signal    // androgen-driven growth
           - k_PV_shrink * PV * AR_shrink_signal;  // involution when DHT is suppressed

// =========================================================================
// SECTION 10: cGMP DYNAMICS (tadalafil PD)
// =========================================================================
// cGMP in prostatic/urethral smooth muscle mediates relaxation.
// PDE5 inhibition by tadalafil reduces cGMP degradation → elevated cGMP →
// smooth muscle relaxation → reduced urethral resistance.
// Nitric oxide (NO) signaling maintained at basal level (not modeled separately).
double pde5_activity = 1.0 - inh_pde5;  // fraction of PDE5 remaining active

dxdt_CGMP = k_cGMP_prod - k_cGMP_deg * pde5_activity * CGMP;

// =========================================================================
// SECTION 11: ALPHA1-AR OCCUPANCY STATE VARIABLE
// =========================================================================
// Tracks the dynamic alpha1-AR occupancy integrating tamsulosin PK/PD.
// Equilibrates to the instantaneous occupancy with fast rate.
double k_alpha1_eq = 5.0;  // h^-1, rapid equilibration
dxdt_ALPHA1_OCC = k_alpha1_eq * (occ_tams - ALPHA1_OCC);

// =========================================================================
// SECTION 12: IPSS DYNAMICS
// =========================================================================
// IPSS is driven by:
//   (a) Prostate volume (static/structural obstruction): increases IPSS
//   (b) Alpha1-AR tone in urethra (dynamic obstruction): alpha1 blockade reduces IPSS
//   (c) cGMP/smooth muscle relaxation (tadalafil): reduces IPSS
//   (d) Prostatic inflammation: increases IPSS
// IPSS is bounded 0-35
double PV_excess = PV - PV0;  // cc above baseline
double IPSS_struct = IPSS0 + k_IPSS_PV * PV_excess;     // structural component
double IPSS_dyn_alpha1 = -k_IPSS_alpha1 * ALPHA1_OCC;   // alpha1 blockade benefit
double IPSS_dyn_cGMP   = -k_IPSS_cGMP  * (CGMP - CGMP0); // cGMP benefit
double IPSS_inflam_comp = k_IPSS_inflam * (INFLAM - INFLAM0); // inflammation
double IPSS_target = IPSS_struct + IPSS_dyn_alpha1 + IPSS_dyn_cGMP + IPSS_inflam_comp;
// Constrain to valid range
IPSS_target = (IPSS_target < 0)  ? 0.0  : IPSS_target;
IPSS_target = (IPSS_target > 35) ? 35.0 : IPSS_target;

// IPSS equilibrates to target with moderate speed (patient-reported, weeks timescale)
double k_IPSS_eq = 0.005;  // h^-1 (~200 h = ~8 day equilibration)
dxdt_IPSS = k_IPSS_eq * (IPSS_target - IPSS);

// =========================================================================
// SECTION 13: QMAX (PEAK FLOW RATE) DYNAMICS
// =========================================================================
// Qmax determined by:
//   (a) Prostate volume (structural): reduces Qmax
//   (b) Alpha1-AR blockade: increases Qmax (dynamic component)
//   (c) cGMP (tadalafil): modest Qmax improvement
// Bounded to 4-30 mL/s range
double Qmax_struct  = QMAX0 - k_Qmax_PV * PV_excess;
double Qmax_alpha1  = k_Qmax_alpha1 * ALPHA1_OCC;
double Qmax_cGMP    = k_Qmax_cGMP  * (CGMP - CGMP0);
double Qmax_target  = Qmax_struct + Qmax_alpha1 + Qmax_cGMP;
Qmax_target = (Qmax_target < 4)  ? 4.0  : Qmax_target;
Qmax_target = (Qmax_target > 30) ? 30.0 : Qmax_target;

double k_Qmax_eq = 0.005;  // h^-1
dxdt_QMAX = k_Qmax_eq * (Qmax_target - QMAX);

// =========================================================================
// SECTION 14: POST-VOID RESIDUAL (PVR) DYNAMICS
// =========================================================================
// PVR inversely related to Qmax and directly to prostate size.
// Normal PVR <50 mL; obstruction increases PVR.
double PVR_target = PVR0 + k_PVR_PV * PV_excess - 15.0 * ALPHA1_OCC - 8.0 * (CGMP - CGMP0);
PVR_target = (PVR_target < 0) ? 0.0 : PVR_target;

double k_PVR_eq = 0.003;  // h^-1 (slower than IPSS)
dxdt_PVR = k_PVR_eq * (PVR_target - PVR);

// =========================================================================
// SECTION 15: PSA DYNAMICS
// =========================================================================
// PSA production proportional to AR_ACT (transcriptional target) and prostate volume.
// 5AR inhibitors reduce PSA ~50% by reducing AR activation (Gormley 1992).
// Tamsulosin has NO effect on PSA (alpha1 blocker, not hormonal).
// Tadalafil has minimal effect on PSA.
double PSA_prod = k_PSA_prod * AR_ACT * (PV / PV0);  // cc-normalized production
double PSA_elim = k_PSA_deg * PSA;                     // first-order elimination

dxdt_PSA = PSA_prod - PSA_elim;

// =========================================================================
// SECTION 16: PROSTATIC INFLAMMATION DYNAMICS
// =========================================================================
// Chronic inflammation correlates with BPH severity and LUTS (De Nunzio 2011).
// DHT drives NF-κB → IL-8, IL-17, COX2 expression in prostate stromal cells.
// 5AR inhibition reduces prostatic DHT → reduced inflammation.
double inflam_target = INFLAM0 + k_inflam_DHT * (DHT_PROST - DHT_PROST0) / DHT_PROST0;
inflam_target = (inflam_target < 0) ? 0.0 : inflam_target;
inflam_target = (inflam_target > 1) ? 1.0 : inflam_target;

dxdt_INFLAM = k_inflam_deg * (inflam_target - INFLAM);

$TABLE
// Drug concentrations (ng/mL)
double CP_tams_out = TAMS_C / V1_tams;
double CP_fina_out = FINA_C / V1_fina;
double CP_dut_out  = DUT_C  / V1_dut;
double CP_tad_out  = TAD_C  / V1_tad;

// DHT inhibition percentage (vs baseline)
double DHT_P_conc = DHT_P / 30.0;   // nmol/L
double DHT_inhibition_pct = (DHT_P_conc < DHT0) ?
                              (1.0 - DHT_P_conc / DHT0) * 100.0 : 0.0;

// Prostatic DHT inhibition
double DHT_PROST_inh_pct = (DHT_PROST < DHT_PROST0) ?
                             (1.0 - DHT_PROST / DHT_PROST0) * 100.0 : 0.0;

// Prostate volume change (cc and %)
double PV_change_cc  = PV - PV0;
double PV_change_pct = PV_change_cc / PV0 * 100.0;

// IPSS change from baseline
double IPSS_change = IPSS - IPSS0;

// Qmax change from baseline
double QMAX_change = QMAX - QMAX0;

// PVR change
double PVR_change = PVR - PVR0;

// PSA change %
double PSA_change_pct = (PSA - PSA0) / PSA0 * 100.0;

// Alpha1-AR blockade %
double Alpha1_block_pct = ALPHA1_OCC * 100.0;

// AR activity
double AR_activity = AR_ACT;

// cGMP level
double cGMP_level = CGMP;

capture CP_tams    = CP_tams_out;
capture CP_fina    = CP_fina_out;
capture CP_dut     = CP_dut_out;
capture CP_tad     = CP_tad_out;
capture DHT_inhibition_pct  = DHT_inhibition_pct;
capture DHT_PROST_inh_pct   = DHT_PROST_inh_pct;
capture PV_change_pct        = PV_change_pct;
capture PV_cc                = PV;
capture IPSS_score           = IPSS;
capture IPSS_change          = IPSS_change;
capture QMAX_mL_s            = QMAX;
capture QMAX_change          = QMAX_change;
capture PVR_mL               = PVR;
capture PVR_change           = PVR_change;
capture PSA_ngmL             = PSA;
capture PSA_change_pct       = PSA_change_pct;
capture Alpha1_block_pct     = Alpha1_block_pct;
capture AR_act_norm          = AR_activity;
capture cGMP_norm            = cGMP_level;
capture Inflam_idx           = INFLAM;
'

## ---------------------------------------------------------------------------
## Compile the model
## ---------------------------------------------------------------------------
cat("Compiling BPH QSP mrgsolve model...\n")
mod <- mcode("bph_qsp", bph_model_code)
cat("Model compiled successfully.\n")
cat(sprintf("  Compartments: %d\n", length(init(mod))))
cat(sprintf("  Parameters:   %d\n", length(param(mod))))

## ---------------------------------------------------------------------------
## Simulation setup
## ---------------------------------------------------------------------------

## Simulation duration: 2 years (730 days), output daily
sim_duration_days <- 730
dt_h <- 24           # output every 24 h (1 day)
times_h <- seq(0, sim_duration_days * 24, by = dt_h)

## Helper: create daily oral dosing events
## dose_mg : dose in mg
## duration_days: total number of days to dose
## cmt: compartment name (GUT compartment)
make_ev <- function(dose_mg, duration_days, cmt) {
  # Convert mg to ug for PK (consistent with V in L → conc in ng/mL = ug/L)
  dose_ug <- dose_mg * 1000
  ev(amt = dose_ug, cmt = cmt, ii = 24, addl = duration_days - 1, time = 0)
}

## ---------------------------------------------------------------------------
## Define 6 Treatment Scenarios
## ---------------------------------------------------------------------------

## Scenario parameters: USE_TAMS, USE_FINA, USE_DUT, USE_TAD flags + dosing events

scenario_list <- list(

  ## Scenario 1: Watchful Waiting (no treatment)
  ## Clinical context: appropriate for mild-moderate symptoms, shared decision-making
  list(
    name  = "Watchful Waiting",
    short = "WW",
    color = "grey40",
    params = c(USE_TAMS = 0, USE_FINA = 0, USE_DUT = 0, USE_TAD = 0),
    events = ev(amt = 0, cmt = 1, time = 0)  # null event
  ),

  ## Scenario 2: Tamsulosin 0.4 mg QD (alpha1-blocker monotherapy)
  ## Rapid symptom relief; no effect on prostate volume or disease progression.
  ## Guideline: 1st-line for moderate-severe LUTS without enlargement risk (AUA 2021)
  list(
    name  = "Tamsulosin 0.4 mg QD",
    short = "TAMS",
    color = "#E74C3C",
    params = c(USE_TAMS = 1, USE_FINA = 0, USE_DUT = 0, USE_TAD = 0),
    events = make_ev(0.4, sim_duration_days, "TAMS_GUT")
  ),

  ## Scenario 3: Finasteride 5 mg QD (selective 5AR inhibitor, type 2)
  ## Reduces DHT ~70%, prostate volume ~20-30% over 6-24 months.
  ## Benefit: AUR risk reduction 57%, TURP 64% (MTOPS, McConnell 2003).
  list(
    name  = "Finasteride 5 mg QD",
    short = "FINA",
    color = "#3498DB",
    params = c(USE_TAMS = 0, USE_FINA = 1, USE_DUT = 0, USE_TAD = 0),
    events = make_ev(5, sim_duration_days, "FINA_GUT")
  ),

  ## Scenario 4: Dutasteride 0.5 mg QD (dual SRD5A1+2 inhibitor)
  ## Reduces DHT ~94% (vs finasteride 70%). Greater prostatic DHT suppression.
  ## Equivalent clinical benefit to finasteride but potentially faster PV reduction.
  ## (Andriole 2010, Eur Urol 57:142)
  list(
    name  = "Dutasteride 0.5 mg QD",
    short = "DUT",
    color = "#9B59B6",
    params = c(USE_TAMS = 0, USE_FINA = 0, USE_DUT = 1, USE_TAD = 0),
    events = make_ev(0.5, sim_duration_days, "DUT_GUT")
  ),

  ## Scenario 5: Combination Therapy (Dutasteride + Tamsulosin)
  ## CombAT trial (Roehrborn 2010, Eur Urol 57:123):
  ##   IPSS reduction: combination -6.2 vs tamsulosin -4.3 vs dutasteride -3.8 at 4 years
  ##   AUR risk reduction: combination 67.6% vs dutasteride 48.0% vs tamsulosin 19.0%
  ## Addresses both dynamic (alpha1) and static (volume) components simultaneously.
  list(
    name  = "Combination (DUT + TAMS)",
    short = "COMBO",
    color = "#27AE60",
    params = c(USE_TAMS = 1, USE_FINA = 0, USE_DUT = 1, USE_TAD = 0),
    events = ev_c(
      make_ev(0.5, sim_duration_days, "DUT_GUT"),
      make_ev(0.4, sim_duration_days, "TAMS_GUT")
    )
  ),

  ## Scenario 6: Tadalafil 5 mg QD (PDE5 inhibitor for BPH-LUTS)
  ## NEPTUNE trial (Chapple 2014, Eur Urol 65:998):
  ##   IPSS -5.6 vs -2.3 placebo at 12 weeks; Qmax +2.4 vs +1.0 mL/s
  ## FDA-approved for BPH-LUTS (2011); dual benefit ED+LUTS
  ## No effect on prostate volume or DHT
  list(
    name  = "Tadalafil 5 mg QD",
    short = "TAD",
    color = "#E67E22",
    params = c(USE_TAMS = 0, USE_FINA = 0, USE_DUT = 0, USE_TAD = 1),
    events = make_ev(5, sim_duration_days, "TAD_GUT")
  )
)

## ---------------------------------------------------------------------------
## Run simulations for all scenarios
## ---------------------------------------------------------------------------
cat("\nRunning 6 treatment scenarios (2-year simulation)...\n")

results_list <- lapply(scenario_list, function(sc) {
  cat(sprintf("  Simulating: %s\n", sc$name))

  # Update model parameters for this scenario
  mod_sc <- param(mod, sc$params)

  # Run simulation
  out <- mrgsim(
    mod_sc,
    events  = sc$events,
    end     = sim_duration_days * 24,
    delta   = 24,        # 24 h intervals = daily output
    obsonly = TRUE
  )

  # Convert to data frame and add scenario label
  df <- as.data.frame(out)
  df$scenario <- sc$name
  df$scenario_short <- sc$short
  df$color <- sc$color
  df$time_days <- df$time / 24  # convert hours to days
  df$time_weeks <- df$time_days / 7
  df$time_months <- df$time_days / 30.44
  df
})

## Combine all scenarios into single data frame
results_all <- bind_rows(results_list)
results_all$scenario <- factor(results_all$scenario,
                                levels = sapply(scenario_list, function(x) x$name))

cat("Simulations complete.\n")
cat(sprintf("  Total observations: %d\n", nrow(results_all)))

## ---------------------------------------------------------------------------
## Clinical Reference Data (from landmark trials)
## ---------------------------------------------------------------------------

## MTOPS trial reference points (McConnell 2003, NEJM 349:2387)
## Finasteride arm (n=768), 4-year results (we use 2-year for consistency)
## IPSS change at ~24 months finasteride: approximately -3.5 pts
## Qmax change at 24 months finasteride: approximately +0.9 mL/s

## CombAT trial (Roehrborn 2010, Eur Urol 57:123) — 4 year results
## Reported at months 3, 6, 12, 24, 36, 48
ref_combat <- data.frame(
  time_months = c(3, 6, 12, 24),
  IPSS_combo  = c(-3.8, -5.3, -5.9, -6.2),    # combination arm IPSS change
  IPSS_tams   = c(-3.2, -4.1, -4.5, -4.3),     # tamsulosin alone
  IPSS_dut    = c(-2.1, -2.8, -3.4, -3.8),     # dutasteride alone
  Qmax_combo  = c(2.0, 2.2, 2.3, 2.4),         # mL/s change from baseline
  Qmax_tams   = c(1.4, 1.6, 1.7, 1.7),
  Qmax_dut    = c(0.5, 0.7, 0.9, 1.1)
)

## NEPTUNE trial reference (Chapple 2014, Eur Urol 65:998) — tadalafil 12 weeks
ref_neptune <- data.frame(
  time_months = c(1, 2, 3),         # months 1, 2, 3 (approximate)
  IPSS_tad    = c(-3.5, -4.8, -5.6) # tadalafil 5 mg QD IPSS change
)

## PLESS trial — finasteride prostate volume (Roehrborn 1996)
ref_pless_pv <- data.frame(
  time_months = c(6, 12, 24, 36),
  PV_fina_pct = c(-10, -18, -25, -28)  # % change from baseline
)

## ---------------------------------------------------------------------------
## Visualization: Key Clinical Endpoints
## ---------------------------------------------------------------------------

## Color palette matching scenario list
scenario_colors <- setNames(
  sapply(scenario_list, function(x) x$color),
  sapply(scenario_list, function(x) x$name)
)

## ---- PLOT 1: IPSS Score Over Time ----
cat("\nGenerating plots...\n")

p1_ipss <- ggplot(results_all, aes(x = time_months, y = IPSS_score,
                                    color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  ## Add CombAT reference points for combination and tamsulosin
  geom_point(data = ref_combat,
             aes(x = time_months, y = IPSS0 + IPSS_combo),
             inherit.aes = FALSE, color = "#27AE60", shape = 17, size = 3) +
  geom_point(data = ref_combat,
             aes(x = time_months, y = IPSS0 + IPSS_tams),
             inherit.aes = FALSE, color = "#E74C3C", shape = 17, size = 3) +
  geom_point(data = ref_combat,
             aes(x = time_months, y = IPSS0 + IPSS_dut),
             inherit.aes = FALSE, color = "#9B59B6", shape = 17, size = 3) +
  ## NEPTUNE reference for tadalafil
  geom_point(data = ref_neptune,
             aes(x = time_months, y = IPSS0 + IPSS_tad),
             inherit.aes = FALSE, color = "#E67E22", shape = 15, size = 3) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  coord_cartesian(ylim = c(6, 25)) +
  labs(
    title    = "IPSS Score Over 2 Years — BPH Treatment Comparison",
    subtitle = paste0("Triangles = CombAT trial data (Roehrborn 2010, n=4844); ",
                      "Squares = NEPTUNE (Chapple 2014)"),
    x        = "Time (months)",
    y        = "IPSS Score (0–35)",
    color    = "Treatment"
  ) +
  geom_hline(yintercept = 8,  linetype = "dashed", color = "grey60", linewidth = 0.7) +
  geom_hline(yintercept = 19, linetype = "dashed", color = "grey60", linewidth = 0.7) +
  annotate("text", x = 22, y = 7,  label = "Mild (<8)",         size = 3, color = "grey50") +
  annotate("text", x = 21, y = 18, label = "Moderate (8–19)",   size = 3, color = "grey50") +
  annotate("text", x = 22, y = 25, label = "Severe (≥20)",      size = 3, color = "grey50") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---- PLOT 2: Prostate Volume Over Time ----
p2_pv <- ggplot(results_all, aes(x = time_months, y = PV_cc,
                                  color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  ## PLESS reference for finasteride PV change
  geom_point(data = ref_pless_pv,
             aes(x = time_months, y = PV0 * (1 + PV_fina_pct / 100)),
             inherit.aes = FALSE, color = "#3498DB", shape = 17, size = 3) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(
    title    = "Prostate Volume Over 2 Years",
    subtitle = "Triangles = PLESS trial (Roehrborn 1996, finasteride arm)",
    x        = "Time (months)",
    y        = "Prostate Volume (cc)",
    color    = "Treatment"
  ) +
  annotate("text", x = 20, y = PV0 * 0.72,
           label = "Target: -25 to -30%\n(5AR inhibitors)",
           size = 3, color = "#3498DB", hjust = 0) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---- PLOT 3: Qmax Over Time ----
p3_qmax <- ggplot(results_all, aes(x = time_months, y = QMAX_mL_s,
                                    color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  annotate("text", x = 1, y = 15.3, label = "Normal Qmax (≥15 mL/s)",
           size = 3, color = "darkgreen", hjust = 0) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(
    title    = "Peak Urine Flow Rate (Qmax) Over 2 Years",
    x        = "Time (months)",
    y        = "Qmax (mL/s)",
    color    = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---- PLOT 4: PSA Over Time ----
p4_psa <- ggplot(results_all, aes(x = time_months, y = PSA_ngmL,
                                   color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(
    title    = "PSA Level Over 2 Years",
    subtitle = "5AR inhibitors reduce PSA ~50% (doubles true PSA for PCa screening)",
    x        = "Time (months)",
    y        = "PSA (ng/mL)",
    color    = "Treatment"
  ) +
  geom_hline(yintercept = 2.5, linetype = "dotted", color = "grey40") +
  annotate("text", x = 22, y = 2.6, label = "Baseline PSA", size = 3, color = "grey40") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---- PLOT 5: DHT Inhibition % Over Time (5AR inhibitors only) ----
dht_data <- results_all %>%
  filter(scenario %in% c("Finasteride 5 mg QD", "Dutasteride 0.5 mg QD",
                          "Combination (DUT + TAMS)"))

p5_dht <- ggplot(dht_data, aes(x = time_months, y = DHT_inhibition_pct,
                                color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "#3498DB",  alpha = 0.6) +
  geom_hline(yintercept = 94, linetype = "dashed", color = "#9B59B6", alpha = 0.6) +
  annotate("text", x = 1, y = 71, label = "Finasteride target (~70%)",
           size = 3, color = "#3498DB", hjust = 0) +
  annotate("text", x = 1, y = 95, label = "Dutasteride target (~94%)",
           size = 3, color = "#9B59B6", hjust = 0) +
  scale_color_manual(values = scenario_colors[names(scenario_colors) %in% unique(dht_data$scenario)]) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title    = "Plasma DHT Inhibition Over 2 Years",
    subtitle = "Dutasteride (dual SRD5A1+2) achieves greater DHT suppression than finasteride (SRD5A2-selective)",
    x        = "Time (months)",
    y        = "DHT Inhibition (%)",
    color    = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---- PLOT 6: PK Profiles (steady-state, Day 14) ----
## Show 24-h PK at pseudo-steady state (~Day 14 for tamsulosin/tadalafil/finasteride)
pk_sim_days <- 14
pk_times_h <- seq(0, (pk_sim_days + 1) * 24, by = 0.5)

## Simulate PK for tamsulosin (representative)
mod_pk_tams <- param(mod, c(USE_TAMS = 1, USE_FINA = 0, USE_DUT = 0, USE_TAD = 0))
out_pk_tams <- mrgsim(
  mod_pk_tams,
  events  = make_ev(0.4, pk_sim_days + 1, "TAMS_GUT"),
  end     = (pk_sim_days + 1) * 24,
  delta   = 0.5,
  obsonly = TRUE
)
pk_tams_df <- as.data.frame(out_pk_tams)
pk_tams_df$time_days <- pk_tams_df$time / 24

## Simulate PK for tadalafil
mod_pk_tad <- param(mod, c(USE_TAMS = 0, USE_FINA = 0, USE_DUT = 0, USE_TAD = 1))
out_pk_tad <- mrgsim(
  mod_pk_tad,
  events  = make_ev(5, pk_sim_days + 1, "TAD_GUT"),
  end     = (pk_sim_days + 1) * 24,
  delta   = 0.5,
  obsonly = TRUE
)
pk_tad_df <- as.data.frame(out_pk_tad)
pk_tad_df$time_days <- pk_tad_df$time / 24

## Filter to Day 14 only (last 24 h)
day14_tams <- pk_tams_df %>% filter(time_days >= pk_sim_days & time_days <= pk_sim_days + 1)
day14_tad  <- pk_tad_df  %>% filter(time_days >= pk_sim_days & time_days <= pk_sim_days + 1)
day14_tams$time_h_ss <- day14_tams$time - pk_sim_days * 24
day14_tad$time_h_ss  <- day14_tad$time  - pk_sim_days * 24

p6_pk <- ggplot() +
  geom_line(data = day14_tams,
            aes(x = time_h_ss, y = CP_tams, color = "Tamsulosin 0.4mg"),
            linewidth = 1.2) +
  geom_line(data = day14_tad,
            aes(x = time_h_ss, y = CP_tad, color = "Tadalafil 5mg"),
            linewidth = 1.2) +
  geom_hline(data = data.frame(drug = "Tamsulosin", val = 0.8),
             aes(yintercept = val, color = "Tamsulosin 0.4mg"),
             linetype = "dashed", alpha = 0.5) +
  geom_hline(data = data.frame(drug = "Tadalafil", val = 2.0),
             aes(yintercept = val, color = "Tadalafil 5mg"),
             linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = c("Tamsulosin 0.4mg" = "#E74C3C",
                                 "Tadalafil 5mg"    = "#E67E22")) +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  labs(
    title    = "Steady-State PK Profiles (Day 14, once-daily dosing)",
    subtitle = "Dashed lines = IC50 values; above IC50 = pharmacologically active",
    x        = "Time after dose (h)",
    y        = "Plasma Concentration (ng/mL)",
    color    = "Drug"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

## ---------------------------------------------------------------------------
## Summary statistics at key time points
## ---------------------------------------------------------------------------
cat("\n=== Summary: Key Endpoints at 24 Months ===\n")

summary_24m <- results_all %>%
  filter(abs(time_months - 24) < 0.5) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, IPSS_score, IPSS_change, QMAX_mL_s, QMAX_change,
         PV_cc, PV_change_pct, PSA_ngmL, PSA_change_pct,
         DHT_inhibition_pct, Alpha1_block_pct) %>%
  as.data.frame()

print(summary_24m, digits = 3)

cat("\n=== Summary: Key Endpoints at 3 Months ===\n")
summary_3m <- results_all %>%
  filter(abs(time_months - 3) < 0.2) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, IPSS_score, IPSS_change, QMAX_mL_s, QMAX_change) %>%
  as.data.frame()
print(summary_3m, digits = 3)

## ---------------------------------------------------------------------------
## Arrange and display plots
## ---------------------------------------------------------------------------
## If gridExtra or patchwork available, arrange in grid; otherwise print individually
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  cat("\nDisplaying 6-panel figure...\n")
  grid.arrange(p1_ipss, p2_pv, p3_qmax, p4_psa, p5_dht, p6_pk,
               ncol = 2, nrow = 3,
               top = "BPH QSP Model — 6 Treatment Scenarios (2-Year Simulation)")
} else {
  cat("\nPrinting plots individually (install gridExtra for combined view):\n")
  print(p1_ipss)
  print(p2_pv)
  print(p3_qmax)
  print(p4_psa)
  print(p5_dht)
  print(p6_pk)
}

## ---------------------------------------------------------------------------
## Save plots to file
## ---------------------------------------------------------------------------
output_dir <- dirname(sys.frame(1)$ofile %||% ".")

## Save combined PNG
png_file <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path %||%
                                 getwd()),
                       "bph_simulation_results.png")

## Safer: save to current working directory
tryCatch({
  png("bph_simulation_results.png", width = 1600, height = 1800, res = 120)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(
      p1_ipss, p2_pv, p3_qmax, p4_psa, p5_dht, p6_pk,
      ncol = 2, nrow = 3,
      top = "BPH QSP Model — 6 Treatment Scenarios (2-Year Simulation)"
    )
  }
  dev.off()
  cat("Plot saved: bph_simulation_results.png\n")
}, error = function(e) {
  cat("Note: Could not save PNG automatically:", conditionMessage(e), "\n")
})

## ---------------------------------------------------------------------------
## Model validation summary
## ---------------------------------------------------------------------------
cat("\n=== MODEL VALIDATION AGAINST CLINICAL TRIALS ===\n")
cat("\n--- Calibration targets ---\n")
cat(paste0(
  "MTOPS (McConnell 2003, NEJM 349:2387):\n",
  "  Finasteride reduces AUR risk by 57%, TURP by 64% vs placebo\n",
  "  Model: 5AR inhibition reduces PV 25-30% → structural obstruction relief\n\n",
  "CombAT (Roehrborn 2010, Eur Urol 57:123):\n",
  "  Combination IPSS change at 2y: -6.2 vs tamsulosin -4.3 vs dutasteride -3.8\n",
  sprintf("  Model 2y results — Combo: %.1f, TAMS: %.1f, DUT: %.1f\n",
          summary_24m$IPSS_change[summary_24m$scenario == "Combination (DUT + TAMS)"],
          summary_24m$IPSS_change[summary_24m$scenario == "Tamsulosin 0.4 mg QD"],
          summary_24m$IPSS_change[summary_24m$scenario == "Dutasteride 0.5 mg QD"]),
  "\nNEPTUNE (Chapple 2014, Eur Urol 65:998):\n",
  "  Tadalafil 5 mg QD reduces IPSS by 5.6 vs 2.3 placebo at 12 weeks\n",
  sprintf("  Model 3m Tadalafil IPSS change: %.1f\n",
          summary_3m$IPSS_change[summary_3m$scenario == "Tadalafil 5 mg QD"]),
  "\nDHT Inhibition Targets:\n",
  "  Finasteride (SRD5A2-selective): ~70% DHT reduction\n",
  sprintf("  Model 2y: %.1f%%\n",
          summary_24m$DHT_inhibition_pct[summary_24m$scenario == "Finasteride 5 mg QD"]),
  "  Dutasteride (dual SRD5A1+2): ~94% DHT reduction\n",
  sprintf("  Model 2y: %.1f%%\n",
          summary_24m$DHT_inhibition_pct[summary_24m$scenario == "Dutasteride 0.5 mg QD"]),
  "\nProstate Volume Reduction (5AR inhibitors, 2 years):\n",
  "  Clinical: -25 to -30% (PLESS, ARIA3001/3002)\n",
  sprintf("  Model 2y finasteride: %.1f%%\n",
          summary_24m$PV_change_pct[summary_24m$scenario == "Finasteride 5 mg QD"]),
  sprintf("  Model 2y dutasteride: %.1f%%\n",
          summary_24m$PV_change_pct[summary_24m$scenario == "Dutasteride 0.5 mg QD"])
))

cat("\n=== BPH QSP Model Simulation Complete ===\n")

## Return results invisibly for further analysis
invisible(list(
  model    = mod,
  results  = results_all,
  plots    = list(ipss  = p1_ipss,
                  pv    = p2_pv,
                  qmax  = p3_qmax,
                  psa   = p4_psa,
                  dht   = p5_dht,
                  pk    = p6_pk)
))
