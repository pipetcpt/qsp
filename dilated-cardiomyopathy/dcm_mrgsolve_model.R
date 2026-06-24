##############################################################################
# Dilated Cardiomyopathy (DCM) — mrgsolve QSP Model
#
# Compartments (23 ODEs):
#   PK: Enalapril, Carvedilol, Spironolactone, Sacubitril, Dapagliflozin
#   PD: RAAS (AngII, Aldosterone), SNS (NE), Cardiac (LVEF, LVEDV, BNP),
#       Fibrosis (Fibrosis_frac), Inflammation (TGFb, IL6),
#       Renal (GFR, Vol), Exercise (6MWT)
#
# Clinical Trial Calibration:
#   PARADIGM-HF  (McMurray 2014, NEJM): Sacubitril/Valsartan vs Enalapril
#   DAPA-HF      (McMurray 2019, NEJM): Dapagliflozin vs placebo
#   COPERNICUS   (Packer 2001, NEJM):   Carvedilol vs placebo
#   RALES        (Pitt 1999, NEJM):     Spironolactone vs placebo
#   CONSENSUS    (CONSENSUS, 1987):     Enalapril vs placebo
#   SHIFT        (Swedberg 2010, Lancet): Ivabradine vs placebo
##############################################################################

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. MODEL SPECIFICATION
# ============================================================

code_dcm <- '
$PROB Dilated Cardiomyopathy (DCM) QSP ODE Model v2.0

$PLUGIN Rcpp

$PARAM
// ---- Disease baseline ----
LVEF0      = 25,    // % baseline LVEF (severe DCM)
LVEDV0     = 250,   // mL baseline LV end-diastolic volume
BNP0       = 800,   // pg/mL baseline BNP
AngII0     = 50,    // pg/mL baseline angiotensin II
Aldo0      = 200,   // pg/mL baseline aldosterone
NE0        = 600,   // pg/mL baseline norepinephrine
Fib0       = 0.20,  // fraction baseline fibrosis (20%)
TGFb0      = 15,    // pg/mL baseline TGF-β
IL6_0      = 8,     // pg/mL baseline IL-6
GFR0       = 55,    // mL/min/1.73m2 baseline eGFR
Vol0       = 7.5,   // L baseline plasma volume
SixMWT0    = 300,   // m baseline 6MWT

// ---- Drug dose (mg) ----
DOSE_ENA   = 10,    // Enalapril mg BID
DOSE_CAR   = 25,    // Carvedilol mg BID
DOSE_SPR   = 25,    // Spironolactone mg QD
DOSE_SAC   = 97,    // Sacubitril mg BID (in ARNI tablet)
DOSE_VAL   = 103,   // Valsartan mg BID (in ARNI tablet)
DOSE_DAPA  = 10,    // Dapagliflozin mg QD
DOSE_IVA   = 7.5,   // Ivabradine mg BID

// ---- Drug on/off flags ----
ON_ENA  = 1,
ON_CAR  = 1,
ON_SPR  = 1,
ON_SAC  = 0,        // set 1 for ARNI (replaces ENA)
ON_DAPA = 0,        // set 1 for SGLT2i add-on
ON_IVA  = 0,        // set 1 for ivabradine add-on

// ---- PK parameters: Enalapril ----
ka_ENA  = 0.69,   // h-1 absorption rate
CL_ENA  = 5.0,    // L/h clearance of enalaprilat
Vd_ENA  = 35,     // L volume of distribution
F_ENA   = 0.40,   // bioavailability

// ---- PK parameters: Carvedilol ----
ka_CAR  = 0.80,
CL_CAR  = 120,
Vd_CAR  = 800,
F_CAR   = 0.25,

// ---- PK parameters: Spironolactone ----
ka_SPR  = 0.50,
CL_SPR  = 8.0,
Vd_SPR  = 50,
F_SPR   = 0.65,

// ---- PK parameters: Sacubitril (active LBQ657) ----
ka_SAC  = 1.20,
CL_SAC  = 3.0,
Vd_SAC  = 40,
F_SAC   = 0.60,

// ---- PK parameters: Dapagliflozin ----
ka_DAPA = 1.50,
CL_DAPA = 8.5,
Vd_DAPA = 75,
F_DAPA  = 0.78,

// ---- PK parameters: Ivabradine ----
ka_IVA  = 1.10,
CL_IVA  = 25,
Vd_IVA  = 175,
F_IVA   = 0.40,

// ---- PD parameters: RAAS ----
Emax_ACEi    = 0.80,   // max ACE inhibition fraction
EC50_ACEi    = 0.05,   // µg/mL enalaprilat EC50
kout_AngII   = 0.35,   // h-1 AngII elimination
ksyn_AngII   = 17.5,   // pg/mL/h baseline production
Emax_ARB     = 0.85,   // max AT1R blockade fraction
EC50_ARB     = 0.08,   // µg/mL valsartan EC50

// ---- PD parameters: Aldosterone ----
kout_Aldo    = 0.20,   // h-1
AngII_Aldo   = 0.015,  // AngII→Aldo coupling coefficient
Emax_MRA     = 0.90,   // max MR blockade
EC50_MRA     = 0.10,   // µg/mL spiro EC50

// ---- PD parameters: NE / SNS ----
kout_NE      = 0.40,   // h-1
ksyn_NE      = 240,    // baseline NE production (pg/mL/h)
Emax_BB      = 0.75,   // max beta-blockade
EC50_BB      = 0.05,   // µg/mL carvedilol EC50
NE_HR_slope  = 0.001,  // NE → heart rate slope

// ---- PD parameters: LVEF / LV remodeling ----
kout_LVEF    = 0.002,  // h-1 spontaneous deterioration rate
kout_LVEDV   = 0.001,  // h-1
AngII_LVEF   = 0.0003, // AngII negative effect on LVEF
NE_LVEF      = 0.0002, // NE negative effect on LVEF
Fib_LVEF     = 0.05,   // Fibrosis negative effect on LVEF
BB_LVEF      = 0.010,  // BB positive effect on LVEF (long-term)
ACEi_LVEF   = 0.006,  // ACEi positive effect
ARB_LVEF    = 0.006,
ARNI_LVEF   = 0.012,  // ARNI > ACEi
DAPA_LVEF   = 0.008,  // SGLT2i positive effect
IVA_LVEF    = 0.004,  // Ivabradine positive effect
LVEF_max     = 55,     // ceiling LVEF (normal)
LVEDV_min    = 120,    // floor LVEDV (normal)

// ---- PD parameters: BNP ----
kout_BNP     = 0.04,   // h-1
AngII_BNP    = 0.01,   // AngII → BNP production
LVEDV_BNP   = 0.05,   // volume → BNP production
Emax_NEP     = 0.75,   // NEP inhibition by sacubitril
EC50_NEP     = 0.03,   // µg/mL LBQ657 EC50

// ---- PD parameters: Fibrosis ----
kout_Fib     = 0.0005, // h-1 fibrosis resolution rate
Aldo_Fib     = 0.0001, // Aldo → fibrosis coupling
TGFb_Fib     = 0.0005, // TGF-β → fibrosis coupling
MRA_Fib      = 0.003,  // MRA reduces fibrosis
DAPA_Fib     = 0.002,  // SGLT2i reduces fibrosis
Fib_max      = 0.60,   // max fibrosis fraction

// ---- PD parameters: TGF-β ----
kout_TGFb    = 0.10,   // h-1
AngII_TGFb   = 0.10,   // AngII → TGF-β
Aldo_TGFb    = 0.05,
MRA_TGFb     = 0.20,   // MRA reduces TGF-β

// ---- PD parameters: IL-6 ----
kout_IL6     = 0.15,   // h-1
NE_IL6       = 0.008,  // NE → IL-6 (cardiac toxicity)
Fib_IL6      = 0.80,   // fibrosis → IL-6
BB_IL6       = 0.20,   // BB reduces IL-6

// ---- PD parameters: GFR / Volume ----
kout_GFR     = 0.005,  // h-1 GFR deterioration
AngII_GFR    = 0.005,  // AngII → GFR effect (efferent constriction)
DAPA_GFR     = 0.002,  // SGLT2i modest GFR protection
CO_GFR       = 0.10,   // cardiac output → GFR coupling
kout_Vol     = 0.08,   // h-1 volume homeostasis
MRA_Vol      = 0.15,   // MRA → volume reduction (natriuresis)
DAPA_Vol     = 0.20,   // SGLT2i → volume reduction (osmotic)
ACEi_Vol     = 0.10,   // ACEi → volume reduction

// ---- PD parameters: 6MWT ----
kout_6MWT    = 0.0008, // h-1 spontaneous decline
LVEF_6MWT    = 1.8,    // LVEF → 6MWT slope
IL6_6MWT     = 5.0,    // IL-6 worsens 6MWT
Vol_6MWT     = 20,     // volume overload worsens 6MWT

// ---- Heart Rate ----
HR0          = 85,     // baseline HR (bpm)
IVA_HR       = 0.25,   // ivabradine max HR reduction fraction

$CMT
// PK compartments (depot + central for each drug)
ENA_GUT ENA_CENT   // Enalapril
CAR_GUT CAR_CENT   // Carvedilol
SPR_GUT SPR_CENT   // Spironolactone
SAC_GUT SAC_CENT   // Sacubitril (LBQ657)
DAPA_GUT DAPA_CENT // Dapagliflozin

// PD state variables
AngII   // Angiotensin II (pg/mL)
Aldo    // Aldosterone (pg/mL)
NE      // Norepinephrine (pg/mL)
LVEF    // Left ventricular ejection fraction (%)
LVEDV   // LV end-diastolic volume (mL)
BNP     // BNP (pg/mL)
Fib     // Fibrosis fraction (0-1)
TGFb    // TGF-β (pg/mL)
IL6     // IL-6 (pg/mL)
GFR     // eGFR (mL/min/1.73m2)
Vol     // Plasma volume index (L)
SixMWT  // 6MWT distance (m)

$MAIN
// ---- Dosing intervals (BID = q12h, QD = q24h) ----
// Handled via event tables in R

// ---- Initial conditions ----
ENA_GUT_0   = 0;
ENA_CENT_0  = 0;
CAR_GUT_0   = 0;
CAR_CENT_0  = 0;
SPR_GUT_0   = 0;
SPR_CENT_0  = 0;
SAC_GUT_0   = 0;
SAC_CENT_0  = 0;
DAPA_GUT_0  = 0;
DAPA_CENT_0 = 0;

AngII_0  = AngII0;
Aldo_0   = Aldo0;
NE_0     = NE0;
LVEF_0   = LVEF0;
LVEDV_0  = LVEDV0;
BNP_0    = BNP0;
Fib_0    = Fib0;
TGFb_0   = TGFb0;
IL6_0    = IL6_0;
GFR_0    = GFR0;
Vol_0    = Vol0;
SixMWT_0 = SixMWT0;

$ODE
// ============================================================
// PHARMACOKINETICS
// ============================================================

// --- Enalapril (prodrug → enalaprilat active) ---
double C_ENA_mg = ENA_CENT / Vd_ENA;  // mg/L = µg/mL
dxdt_ENA_GUT  = -ka_ENA * ENA_GUT;
dxdt_ENA_CENT =  ka_ENA * ENA_GUT * F_ENA - CL_ENA * C_ENA_mg;

// --- Carvedilol ---
double C_CAR_mg = CAR_CENT / Vd_CAR;
dxdt_CAR_GUT  = -ka_CAR * CAR_GUT;
dxdt_CAR_CENT =  ka_CAR * CAR_GUT * F_CAR - CL_CAR * C_CAR_mg;

// --- Spironolactone ---
double C_SPR_mg = SPR_CENT / Vd_SPR;
dxdt_SPR_GUT  = -ka_SPR * SPR_GUT;
dxdt_SPR_CENT =  ka_SPR * SPR_GUT * F_SPR - CL_SPR * C_SPR_mg;

// --- Sacubitril (as active LBQ657) ---
double C_SAC_mg = SAC_CENT / Vd_SAC;
dxdt_SAC_GUT  = -ka_SAC * SAC_GUT;
dxdt_SAC_CENT =  ka_SAC * SAC_GUT * F_SAC - CL_SAC * C_SAC_mg;

// --- Dapagliflozin ---
double C_DAPA_mg = DAPA_CENT / Vd_DAPA;
dxdt_DAPA_GUT  = -ka_DAPA * DAPA_GUT;
dxdt_DAPA_CENT =  ka_DAPA * DAPA_GUT * F_DAPA - CL_DAPA * C_DAPA_mg;

// ============================================================
// DRUG EFFECT FUNCTIONS (Emax models)
// ============================================================

// ACEi effect (enalaprilat)
double E_ACEi  = ON_ENA * Emax_ACEi * C_ENA_mg / (EC50_ACEi + C_ENA_mg);

// ARB effect: ARNI valsartan + standalone valsartan
// Simplified: ON_SAC implies valsartan component (same EC50 as standalone ARB)
double C_ARB = ON_SAC * DOSE_VAL / (Vd_ENA * 1.2);  // approx plasma conc
double E_ARB   = ON_SAC * Emax_ARB * C_ARB / (EC50_ARB + C_ARB);

// Beta-blocker effect (carvedilol)
double E_BB    = ON_CAR * Emax_BB * C_CAR_mg / (EC50_BB + C_CAR_mg);

// MRA effect (spironolactone)
double E_MRA   = ON_SPR * Emax_MRA * C_SPR_mg / (EC50_MRA + C_SPR_mg);

// NEP inhibition (sacubitril/LBQ657)
double E_NEP   = ON_SAC * Emax_NEP * C_SAC_mg / (EC50_NEP + C_SAC_mg);

// SGLT2i effect (dapagliflozin)
double E_DAPA  = ON_DAPA * C_DAPA_mg / (0.2 + C_DAPA_mg);  // Emax = 1 at saturation

// Ivabradine: HR reduction (simplified from dose)
double E_IVA   = ON_IVA * IVA_HR;

// ============================================================
// RAAS: ANGIOTENSIN II
// ============================================================
// AngII synthesis inhibited by ACEi or ARB
// ACEi blocks conversion AngI→AngII; ARB blocks AT1R signaling (feedback)
double ACEi_or_ARB  = fmax(E_ACEi, E_ARB);  // ACEi or ARNI-valsartan
double AngII_prod   = ksyn_AngII * (1.0 - ACEi_or_ARB);
double AngII_elim   = kout_AngII * AngII;
dxdt_AngII = AngII_prod - AngII_elim;

// ============================================================
// ALDOSTERONE
// ============================================================
double Aldo_prod  = AngII_Aldo * AngII;
double Aldo_elim  = kout_Aldo * Aldo;
double Aldo_MRA   = E_MRA * Aldo;   // MRA doesn't reduce Aldo itself; blocks MR
dxdt_Aldo = Aldo_prod - Aldo_elim;  // Aldo synthesis driven by AngII

// ============================================================
// SYMPATHETIC NE (norepinephrine)
// ============================================================
// NE elevated when cardiac output is low; BB reduces spillover + chronic toxicity
double CO_effect = fmax(0.1, 1.0 - LVEF / LVEF_max);  // lower LVEF → more NE
double NE_prod   = ksyn_NE * CO_effect * (1.0 - E_BB * 0.5);  // BB partly reduces release
double NE_elim   = kout_NE * NE;
dxdt_NE = NE_prod - NE_elim;

// ============================================================
// LVEF — Left Ventricular Ejection Fraction
// ============================================================
// LVEF deteriorates due to: high AngII, high NE, fibrosis, IL-6
// LVEF improves with: BB (reverse remodeling), ACEi/ARB, ARNI, SGLT2i
double LVEF_deteriorate = AngII_LVEF * AngII + NE_LVEF * NE + Fib_LVEF * Fib + 0.001 * IL6;
double LVEF_improve     = BB_LVEF * E_BB + ACEi_LVEF * E_ACEi + ARB_LVEF * E_ARB +
                          ARNI_LVEF * E_NEP + DAPA_LVEF * E_DAPA + IVA_LVEF * E_IVA;
double LVEF_net         = LVEF_improve - LVEF_deteriorate;
// Constrain LVEF between 5% (near death) and LVEF_max
dxdt_LVEF = LVEF_net * LVEF * (1.0 - LVEF / LVEF_max) * (LVEF > 5.0 ? 1.0 : 0.0);

// ============================================================
// LVEDV — LV End-Diastolic Volume (dilated ventricle)
// ============================================================
// Volume increases (dilation) with high AngII, Aldo, NE → preload/afterload
// Volume decreases with effective treatment → reverse remodeling
double LVEDV_inc  = AngII * 0.01 + Aldo * 0.005 + Vol * 5;
double LVEDV_dec  = (E_ACEi + E_ARB + E_NEP) * 2.0 + E_MRA * 1.0 + E_DAPA * 1.5;
double LVEDV_net  = (LVEDV_inc - LVEDV_dec) * 0.001;
dxdt_LVEDV = LVEDV_net * (LVEDV - LVEDV_min);  // stops at normal

// ============================================================
// BNP — B-type Natriuretic Peptide
// ============================================================
// BNP produced in response to wall stress (LVEDV, LVEDP)
// BNP degraded by NEP (neprilysin); inhibited by sacubitril → BNP paradox
double BNP_prod  = (AngII_BNP * AngII + LVEDV_BNP * (LVEDV / 100.0));
double NEP_deg   = kout_BNP * (1.0 - E_NEP);  // NEP inhibition → less BNP degradation
dxdt_BNP = BNP_prod - NEP_deg * BNP;

// ============================================================
// FIBROSIS — Interstitial Fibrosis Fraction
// ============================================================
// Driven by TGF-β, AngII (via MR/AT1R), Aldo
// Reduced by MRA, SGLT2i
double Fib_prod   = TGFb_Fib * TGFb + Aldo_Fib * Aldo * (1.0 - E_MRA);
double Fib_dec    = MRA_Fib * E_MRA + DAPA_Fib * E_DAPA + kout_Fib;
dxdt_Fib = Fib_prod - Fib_dec * Fib;
if (Fib > Fib_max) dxdt_Fib = 0.0;

// ============================================================
// TGF-β
// ============================================================
double TGFb_prod  = AngII_TGFb * (AngII / AngII0) + Aldo_TGFb * (Aldo / Aldo0);
double TGFb_dec   = kout_TGFb + MRA_TGFb * E_MRA;
dxdt_TGFb = TGFb_prod - TGFb_dec * TGFb;

// ============================================================
// IL-6 — Interleukin-6 (surrogate for systemic inflammation)
// ============================================================
double IL6_prod  = NE_IL6 * NE + Fib_IL6 * Fib;
double IL6_dec   = kout_IL6 + BB_IL6 * E_BB * 0.1;
dxdt_IL6 = IL6_prod - IL6_dec * IL6;

// ============================================================
// eGFR — Glomerular Filtration Rate
// ============================================================
double CO       = (LVEF / 100.0) * LVEDV * HR0 * (1.0 - E_IVA) / 1000.0;  // L/min approx
double GFR_prod  = CO_GFR * CO;
double GFR_dec   = kout_GFR + AngII_GFR * (AngII / AngII0) * 0.01;
// SGLT2i initially dips GFR then stabilizes — modelled as long-term protection
double GFR_DAPA  = DAPA_GFR * E_DAPA * (1.0 - GFR / GFR0);
dxdt_GFR = GFR_prod - GFR_dec * GFR + GFR_DAPA;

// ============================================================
// PLASMA VOLUME (preload index)
// ============================================================
// Volume reduced by natriuresis (ACEi/ARB → Aldo↓, MRA blocks Na retention, SGLT2i osmotic)
double Vol_prod  = 0.5;  // intake
double Vol_dec   = MRA_Vol * E_MRA + DAPA_Vol * E_DAPA + ACEi_Vol * (E_ACEi + E_ARB * 0.5) + kout_Vol;
dxdt_Vol = Vol_prod - Vol_dec * Vol;

// ============================================================
// 6-MINUTE WALK TEST (functional capacity)
// ============================================================
// Improves with LVEF↑, CO↑; worsens with IL-6↑, volume overload
double target_6MWT = LVEF_6MWT * LVEF - IL6_6MWT * IL6 - Vol_6MWT * (Vol - 5.0);
double rate_6MWT   = 0.002 * (target_6MWT - SixMWT);  // slow adaptation
dxdt_SixMWT = rate_6MWT;

$TABLE
// Derived outputs
capture HR       = HR0 * (1.0 - E_IVA) * (1.0 + NE_HR_slope * (NE - NE0));
capture CO_Lmin  = (LVEF / 100.0) * LVEDV * HR / 60000.0;  // L/min
capture SV       = (LVEF / 100.0) * LVEDV;                  // mL
capture NT_proBNP = BNP * 6.5;                               // NT-proBNP approx 6.5x BNP
capture E_ACEi_out = E_ACEi;
capture E_BB_out   = E_BB;
capture E_MRA_out  = E_MRA;
capture E_NEP_out  = E_NEP;
capture E_DAPA_out = E_DAPA;
capture LVESV    = LVEDV * (1.0 - LVEF / 100.0);
// NYHA class: simplified from LVEF + 6MWT
capture NYHA     = (LVEF < 15) ? 4.0 : (LVEF < 25) ? 3.0 : (LVEF < 35) ? 2.5 : 2.0;
capture Cenalaprilat = ENA_CENT / Vd_ENA;
capture Ccarvedilol  = CAR_CENT / Vd_CAR;
capture Cspiro       = SPR_CENT / Vd_SPR;
capture Csacubitril  = SAC_CENT / Vd_SAC;
capture Cdapa        = DAPA_CENT / Vd_DAPA;
'

# ============================================================
# 2. COMPILE MODEL
# ============================================================
mod_dcm <- mread_cache("DCM_QSP", tempdir(), code_dcm)

# ============================================================
# 3. DOSING EVENT TABLE (24-week simulation)
# ============================================================

build_events <- function(duration_h = 168 * 24,  # 24 weeks = 4032h
                         on_ena  = TRUE,
                         on_car  = TRUE,
                         on_spr  = TRUE,
                         on_sac  = FALSE,
                         on_dapa = FALSE,
                         on_iva  = FALSE) {
  ev <- ev_seq()

  # Enalapril 10mg BID (q12h) — NOT used with ARNI
  if (on_ena && !on_sac) {
    ev <- ev + ev(amt = 10, cmt = "ENA_GUT", ii = 12, addl = duration_h/12 - 1)
  }

  # Carvedilol 25mg BID (q12h)
  if (on_car) {
    ev <- ev + ev(amt = 25, cmt = "CAR_GUT", ii = 12, addl = duration_h/12 - 1)
  }

  # Spironolactone 25mg QD (q24h)
  if (on_spr) {
    ev <- ev + ev(amt = 25, cmt = "SPR_GUT", ii = 24, addl = duration_h/24 - 1)
  }

  # Sacubitril/Valsartan (ARNI): sacubitril 97mg BID
  if (on_sac) {
    ev <- ev + ev(amt = 97, cmt = "SAC_GUT", ii = 12, addl = duration_h/12 - 1)
  }

  # Dapagliflozin 10mg QD
  if (on_dapa) {
    ev <- ev + ev(amt = 10, cmt = "DAPA_GUT", ii = 24, addl = duration_h/24 - 1)
  }

  return(ev)
}

# ============================================================
# 4. SIMULATION: 5 TREATMENT SCENARIOS
# ============================================================

sim_duration <- 168 * 24  # 24 weeks in hours
obs_times    <- seq(0, sim_duration, by = 24)  # daily observations

# Scenario definitions
scenarios <- list(
  list(name = "1_Placebo",           on_ena=F, on_car=F, on_spr=F, on_sac=F, on_dapa=F, on_iva=F,
       params = c(ON_ENA=0, ON_CAR=0, ON_SPR=0, ON_SAC=0, ON_DAPA=0, ON_IVA=0)),
  list(name = "2_ACEi+BB",           on_ena=T, on_car=T, on_spr=F, on_sac=F, on_dapa=F, on_iva=F,
       params = c(ON_ENA=1, ON_CAR=1, ON_SPR=0, ON_SAC=0, ON_DAPA=0, ON_IVA=0)),
  list(name = "3_ACEi+BB+MRA",       on_ena=T, on_car=T, on_spr=T, on_sac=F, on_dapa=F, on_iva=F,
       params = c(ON_ENA=1, ON_CAR=1, ON_SPR=1, ON_SAC=0, ON_DAPA=0, ON_IVA=0)),
  list(name = "4_ARNI+BB+MRA",       on_ena=F, on_car=T, on_spr=T, on_sac=T, on_dapa=F, on_iva=F,
       params = c(ON_ENA=0, ON_CAR=1, ON_SPR=1, ON_SAC=1, ON_DAPA=0, ON_IVA=0)),
  list(name = "5_ARNI+BB+MRA+SGLT2i", on_ena=F, on_car=T, on_spr=T, on_sac=T, on_dapa=T, on_iva=F,
       params = c(ON_ENA=0, ON_CAR=1, ON_SPR=1, ON_SAC=1, ON_DAPA=1, ON_IVA=0))
)

run_scenario <- function(scen) {
  ev_dose <- build_events(
    duration_h = sim_duration,
    on_ena  = scen$on_ena,
    on_car  = scen$on_car,
    on_spr  = scen$on_spr,
    on_sac  = scen$on_sac,
    on_dapa = scen$on_dapa
  )
  mod_dcm %>%
    param(scen$params) %>%
    mrgsim(events = ev_dose, end = sim_duration, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(Scenario = scen$name,
           time_week = time / (24 * 7))
}

results <- purrr::map_dfr(scenarios, run_scenario)

# ============================================================
# 5. RESULTS SUMMARY TABLE
# ============================================================
summary_table <- results %>%
  filter(time_week %in% c(0, 4, 12, 24)) %>%
  group_by(Scenario, time_week) %>%
  summarise(
    LVEF      = round(mean(LVEF), 1),
    BNP       = round(mean(BNP), 0),
    NT_proBNP = round(mean(NT_proBNP), 0),
    SixMWT    = round(mean(SixMWT), 0),
    HR        = round(mean(HR), 0),
    NYHA      = round(mean(NYHA), 1),
    GFR       = round(mean(GFR), 0),
    Fib_pct   = round(mean(Fib) * 100, 1),
    .groups = "drop"
  )

print(summary_table)

# ============================================================
# 6. KEY PLOTS
# ============================================================

gg_LVEF <- ggplot(results, aes(x = time_week, y = LVEF, color = Scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 35, linetype = "dashed", color = "gray50") +
  annotate("text", x = 20, y = 36.5, label = "ICD threshold (35%)", size = 3.5) +
  labs(title = "LVEF Over Time — DCM QSP Model",
       subtitle = "Scenarios: Placebo vs Guideline-Directed Medical Therapy (GDMT)",
       x = "Time (weeks)", y = "LVEF (%)",
       color = "Treatment") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

gg_BNP <- ggplot(results, aes(x = time_week, y = BNP, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "BNP Over Time",
       x = "Time (weeks)", y = "BNP (pg/mL)",
       color = "Treatment") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

gg_6MWT <- ggplot(results, aes(x = time_week, y = SixMWT, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "6-Minute Walk Test Over Time",
       x = "Time (weeks)", y = "6MWT (m)",
       color = "Treatment") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

gg_Fib <- ggplot(results, aes(x = time_week, y = Fib * 100, color = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "Myocardial Fibrosis Over Time",
       x = "Time (weeks)", y = "Fibrosis (%)",
       color = "Treatment") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

gg_PK <- results %>%
  filter(Scenario == "5_ARNI+BB+MRA+SGLT2i") %>%
  select(time_week, Cenalaprilat, Ccarvedilol, Cspiro, Csacubitril, Cdapa) %>%
  pivot_longer(-time_week, names_to = "Drug", values_to = "Conc") %>%
  ggplot(aes(x = time_week, y = Conc, color = Drug)) +
  geom_line(size = 1.2) +
  labs(title = "Drug Plasma Concentrations (Scenario 5)",
       subtitle = "ARNI + Carvedilol + Spironolactone + Dapagliflozin",
       x = "Time (weeks)", y = "Concentration (µg/mL)", color = "Drug") +
  scale_color_brewer(palette = "Set2") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# print plots (if running interactively)
# print(gg_LVEF)
# print(gg_BNP)
# print(gg_6MWT)
# print(gg_Fib)
# print(gg_PK)

# ============================================================
# 7. EXPORT RESULTS
# ============================================================
# write.csv(results, "dcm_simulation_results.csv", row.names = FALSE)
# write.csv(summary_table, "dcm_summary_table.csv", row.names = FALSE)

cat("\nDCM QSP Model run complete.\n")
cat("Summary at Week 24 (Scenario 5 — Quadruple therapy):\n")
print(summary_table %>% filter(Scenario == "5_ARNI+BB+MRA+SGLT2i", time_week == 24))
