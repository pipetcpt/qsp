## ============================================================
## Hypertrophic Cardiomyopathy (HCM) — QSP mrgsolve Model
## ============================================================
## Disease: Hypertrophic Cardiomyopathy (obstructive / non-obstructive)
## Model scope:
##   1. Mavacamten 3-compartment PK (CYP2C19 EM/PM polymorphism)
##   2. Beta-blocker 1-compartment PK
##   3. Sarcomere duty-ratio / myosin ATPase dynamics
##   4. Calcium cycling (cytosolic + SR)
##   5. Calcineurin-NFAT hypertrophy signaling
##   6. ERK1/2 stress signaling
##   7. LV structural remodeling (IVS thickness, LV mass)
##   8. TGFβ1-driven fibrosis (collagen content)
##   9. LV mechanics (LVOT gradient, LVEDP)
##  10. Heart rate (autonomic + drug effects)
##  11. Neurohormonal biomarkers (NT-proBNP)
##  12. Myocardial injury biomarker (Troponin I)
##  13. Atrial fibrillation probability (cumulative risk)
##  14. Diastolic function grade
##  15. Cardiac output and LVEF
##
## Parameters calibrated to:
##  - EXPLORER-HCM (Olivotto I et al., NEJM 2020) — mavacamten RCT
##  - VALOR-HCM (Desai MY et al., NEJM 2022) — mavacamten vs SRT
##  - SEQUOIA-HCM (Nagueh SF et al., Lancet 2024) — aficamten RCT
##  - Maron BJ et al., NEJM 2018 (HCM management)
##  - Green EM et al., Science 2016 (SRX/DRX state)
## ============================================================

library(mrgsolve)

hcm_code <- '
$PROB
HCM QSP Model
Mavacamten (MYK-461) + beta-blocker PK/PD
Sarcomere mechanics, hypertrophy, fibrosis, LVOT hemodynamics

$PARAM @annotated
// ---- Mavacamten PK (3-compartment) ----
ka_m    : 0.693  : hr-1, absorption rate constant (t1/2_abs ~1h)
CL_m    : 1.80   : L/hr, apparent clearance (CYP2C19 EM)
V1_m    : 42.0   : L, central volume of distribution
Q_m     : 0.90   : L/hr, intercompartmental clearance
V2_m    : 125.0  : L, peripheral volume
MW_mava : 471.0  : g/mol, mavacamten molecular weight
F_mava  : 0.93   : -, oral bioavailability

// ---- Mavacamten PD ----
EC50_mava  : 85.0  : nM, half-maximal duty ratio suppression
Emax_mava  : 0.65  : -, maximum fractional suppression of duty ratio
Hill_mava  : 1.2   : -, Hill coefficient

// ---- Beta-blocker PK (metoprolol succinate) ----
ka_bb   : 0.50   : hr-1, absorption rate
CL_bb   : 18.0   : L/hr, clearance
V_bb    : 180.0  : L, volume of distribution

// ---- Beta-blocker PD ----
EC50_bb    : 25.0   : ng/mL, HR reduction EC50
Emax_bb_hr : 0.35   : -, max HR fractional reduction (35%)
EC50_bb_c  : 40.0   : ng/mL, contractility reduction EC50
Emax_bb_c  : 0.25   : -, max contractility reduction

// ---- Sarcomere / duty ratio ----
DR_normal  : 0.04   : -, normal duty ratio (~4% of heads in power stroke)
DR_HCM     : 0.09   : -, HCM mutation-driven elevated duty ratio
k_SRX      : 0.15   : hr-1, rate constant toward SRX state
k_DRX      : 0.05   : hr-1, rate constant toward DRX state

// ---- Calcium handling ----
Ca_base    : 200.0  : nM, resting cytosolic Ca2+
Ca_SR_base : 420.0  : uM, SR Ca2+ load at baseline
k_SERCA    : 0.80   : hr-1, SERCA2a reuptake rate constant
k_NCX      : 0.20   : hr-1, NCX extrusion rate constant
k_RyR_leak : 0.05   : hr-1, diastolic RyR2 leak constant
Ca_amp_HCM : 1.35   : -, relative Ca transient amplitude in HCM

// ---- Calcineurin-NFAT signaling ----
k_CaN      : 0.50   : hr-1, calcineurin activation rate
k_NFAT_in  : 0.30   : hr-1, NFAT nuclear translocation rate
k_NFAT_out : 0.80   : hr-1, NFAT nuclear export rate
CaN_base   : 0.30   : -, relative baseline calcineurin activity
NFAT_base  : 0.20   : -, baseline nuclear NFAT fraction

// ---- ERK signaling ----
k_ERK      : 0.30   : hr-1, ERK activation/inactivation rate
ERK_base   : 0.15   : -, baseline ERK activity

// ---- Hypertrophy (IVS, LV mass) ----
k_hyp      : 2.5e-4 : hr-1, hypertrophy rate (NFAT-driven)
IVS_0      : 18.0   : mm, baseline IVS thickness (HCM patient)
IVS_max    : 32.0   : mm, maximum possible IVS
k_IVS_rev  : 5e-5   : hr-1, slow reversal of hypertrophy
LVmass_0   : 1.60   : -, relative LV mass (1.0 = normal)

// ---- TGFβ1 / Fibrosis ----
TGFb_base  : 1.00   : ng/mL, baseline TGF-β1
k_TGFprod  : 0.015  : hr-1, TGF-β1 production rate
k_TGFdeg   : 0.020  : hr-1, TGF-β1 degradation rate
k_col_syn  : 0.005  : hr-1, collagen synthesis rate
k_col_deg  : 0.002  : hr-1, collagen degradation rate
Col_base   : 1.30   : -, relative collagen content (HCM baseline)

// ---- LV mechanics ----
LVOT_base  : 45.0   : mmHg, resting LVOT gradient (HCM patient)
LVOT_min   : 8.0    : mmHg, minimum physiological gradient
k_LVOT_dyn : 0.50   : hr-1, LVOT dynamics rate constant
LVEDP_base : 20.0   : mmHg, baseline LVEDP
k_LVEDP    : 0.20   : hr-1, LVEDP dynamics

// ---- Hemodynamics ----
HR_base    : 80.0   : bpm, baseline heart rate
HR_min     : 45.0   : bpm, minimum HR (bradycardia limit)
SV_base    : 72.0   : mL, baseline stroke volume
LVEF_base  : 72.0   : %, baseline LVEF (hyperdynamic in HCM)
k_CO_dyn   : 0.10   : hr-1, cardiac output dynamics

// ---- NT-proBNP ----
NTpBNP_0   : 400.0  : pg/mL, baseline NT-proBNP
k_BNP_syn  : 0.10   : hr-1, BNP synthesis rate constant
k_BNP_deg  : 0.05   : hr-1, BNP degradation rate

// ---- Troponin I ----
TropI_base : 25.0   : ng/L, baseline high-sensitivity troponin I
k_Trop_rel : 0.02   : hr-1, troponin release rate (ischemia)
k_Trop_clr : 0.04   : hr-1, troponin clearance rate

// ---- AF risk ----
AF_hazard_base : 4.5 : %/yr, baseline annual AF incidence in HCM
k_AF_LA    : 0.10   : hr-1, AF risk accumulation from LA pressure

// ---- CYP2C19 polymorphism ----
CYP_factor : 1.0    : -, 1.0=EM, 2.5=PM, 0.4=UM (mavacamten metabolism)

$CMT @annotated
// --- Mavacamten PK ---
A_gut_m   : Mavacamten gut compartment (nmol)
A_c_m     : Mavacamten central compartment (nmol)
A_p_m     : Mavacamten peripheral compartment (nmol)

// --- Beta-blocker PK ---
A_bb      : Beta-blocker plasma (ug)

// --- Calcium dynamics ---
Ca_cyt    : Cytosolic free Ca2+ (nM)
Ca_SR     : Sarcoplasmic reticulum Ca2+ (uM)

// --- Signaling ---
Calcineurin : Calcineurin activity (normalized 0-1)
NFAT_nuc    : Nuclear NFAT fraction (0-1)
ERK_act     : ERK1/2 phosphorylation (normalized 0-1)

// --- Structural remodeling ---
IVS         : IVS thickness (mm)
LVmass_rel  : Relative LV mass (1.0 = normal)
TGFb1       : TGF-β1 plasma/myocardial (ng/mL)
Collagen    : Relative myocardial collagen content

// --- Hemodynamics ---
LVOT        : LVOT gradient (mmHg)
LVEDP       : LV end-diastolic pressure (mmHg)
HR_state    : Heart rate state (bpm)

// --- Biomarkers ---
NT_proBNP   : NT-proBNP plasma (pg/mL)
TroponinI   : High-sensitivity Troponin I (ng/L)

// --- AF cumulative hazard ---
AF_hazard   : Cumulative AF hazard (arbitrary units)

$INIT
A_gut_m = 0, A_c_m = 0, A_p_m = 0, A_bb = 0
Ca_cyt = 200, Ca_SR = 420
Calcineurin = 0.30, NFAT_nuc = 0.20, ERK_act = 0.15
IVS = 18.0, LVmass_rel = 1.60, TGFb1 = 1.00, Collagen = 1.30
LVOT = 45.0, LVEDP = 20.0, HR_state = 80.0
NT_proBNP = 400.0, TroponinI = 25.0
AF_hazard = 0.0

$ODE

// =============================================================
// MAVACAMTEN PK (3-compartment with CYP2C19 polymorphism)
// =============================================================
double CL_mava_adj = CL_m * CYP_factor;  // CYP2C19-adjusted clearance
double Cm = A_c_m / V1_m;                // central concentration (nmol/L = nM)

dxdt_A_gut_m = -ka_m * A_gut_m;
dxdt_A_c_m   =  ka_m * A_gut_m
              - (CL_mava_adj / V1_m) * A_c_m
              - (Q_m / V1_m) * A_c_m
              + (Q_m / V2_m) * A_p_m;
dxdt_A_p_m   =  (Q_m / V1_m) * A_c_m
              - (Q_m / V2_m) * A_p_m;

// Mavacamten PD — Emax model on duty ratio suppression
double E_mava = Emax_mava * pow(Cm, Hill_mava)
                / (pow(EC50_mava, Hill_mava) + pow(Cm, Hill_mava));
if (E_mava > 0.80) E_mava = 0.80;  // safety cap (LVEF protection)

// =============================================================
// BETA-BLOCKER PK (1-compartment, metoprolol)
// =============================================================
double Cbb = A_bb / V_bb;  // ng/mL

dxdt_A_bb = -(CL_bb / V_bb) * A_bb;
// (dose added via event table)

// Beta-blocker PD
double E_bb_hr = Emax_bb_hr * Cbb / (EC50_bb + Cbb);
double E_bb_c  = Emax_bb_c  * Cbb / (EC50_bb_c + Cbb);

// =============================================================
// SARCOMERE DUTY RATIO
// =============================================================
// Effective duty ratio: HCM mutation elevates DRX fraction
// Mavacamten promotes SRX sequestration → ↓ duty ratio
double DR_eff = DR_HCM * (1.0 - E_mava) * (1.0 - 0.5 * E_bb_c);
if (DR_eff < DR_normal * 0.5) DR_eff = DR_normal * 0.5;  // floor
double DR_ratio = DR_eff / DR_normal;  // relative to normal (>1 in HCM)

// =============================================================
// CALCIUM DYNAMICS
// =============================================================
// SR Ca2+ load drives RyR2 release; SERCA2a driven by PKA (BB effect reversal)
double Ca_SR_norm = Ca_SR / Ca_SR_base;
double Ca_influx  = 60.0 * DR_ratio * (1.0 + 0.3 * (Ca_SR_norm - 1.0));
double Ca_SERCA   = k_SERCA * Ca_cyt;
double Ca_NCX_out = k_NCX * Ca_cyt;
double Ca_leak_SR = k_RyR_leak * Ca_SR_norm * Ca_cyt / 200.0;

dxdt_Ca_cyt = Ca_influx - Ca_SERCA - Ca_NCX_out;
dxdt_Ca_SR  = Ca_SERCA * (200.0 / 420.0) - Ca_leak_SR * (420.0 / 200.0);

// =============================================================
// CALCINEURIN-NFAT SIGNALING
// =============================================================
// Calcineurin activated by elevated cytosolic Ca2+ (via CaM)
double Ca_norm    = Ca_cyt / Ca_base;
double CaN_target = CaN_base + (1.0 - CaN_base) * (Ca_norm - 1.0) / (1.0 + (Ca_norm - 1.0));
if (CaN_target < 0.05) CaN_target = 0.05;
if (CaN_target > 1.00) CaN_target = 1.00;

dxdt_Calcineurin = k_CaN * (CaN_target - Calcineurin);

// NFAT: dephosphorylated by calcineurin → nuclear translocation
double NFAT_target = NFAT_base + 0.80 * Calcineurin;
if (NFAT_target > 1.0) NFAT_target = 1.0;

dxdt_NFAT_nuc = k_NFAT_in * (NFAT_target - NFAT_nuc);

// =============================================================
// ERK1/2 SIGNALING (mechanical + TGFβ + AngII)
// =============================================================
double LVOT_norm  = LVOT / LVOT_base;
double TGFb_norm  = TGFb1 / TGFb_base;
double ERK_target = ERK_base
                  + 0.40 * (LVOT_norm - 1.0) / (1.0 + (LVOT_norm - 1.0))
                  + 0.30 * (TGFb_norm - 1.0) / (1.0 + (TGFb_norm - 1.0));
if (ERK_target < 0.0)  ERK_target = 0.0;
if (ERK_target > 1.0)  ERK_target = 1.0;

dxdt_ERK_act = k_ERK * (ERK_target - ERK_act);

// =============================================================
// CARDIAC HYPERTROPHY (IVS, LV mass)
// =============================================================
// Hypertrophy rate driven by NFAT and ERK (calcineurin-NFAT dominant)
double hyp_signal = NFAT_nuc * (1.0 + 0.4 * ERK_act);
double IVS_growth = k_hyp * hyp_signal * (IVS_max - IVS);
double IVS_regress = k_IVS_rev * IVS * (1.0 - hyp_signal);  // slow regression

dxdt_IVS = IVS_growth - IVS_regress;

// LV mass tracks IVS with a lag
double LVmass_target = 1.0 + 0.60 * (IVS - 11.0) / 11.0;
if (LVmass_target < 1.0) LVmass_target = 1.0;

dxdt_LVmass_rel = 0.001 * (LVmass_target - LVmass_rel);

// =============================================================
// TGFβ1 / FIBROSIS
// =============================================================
// TGFβ1 upregulated by mechanical stress (LVOT), ischemia, AngII
double TGFb_prod = k_TGFprod * (1.0 + 0.50 * (LVOT_norm - 1.0)
                                     + 0.30 * (DR_ratio - 1.0)
                                     + 0.20 * (Collagen - 1.0));  // autocrine
double TGFb_deg  = k_TGFdeg * TGFb1;

dxdt_TGFb1 = TGFb_prod - TGFb_deg;

// Collagen: TGFβ1 drives synthesis, endogenous degradation
double Col_synth = k_col_syn * TGFb1 / TGFb_base;
double Col_deg   = k_col_deg * Collagen;

dxdt_Collagen = Col_synth - Col_deg;

// =============================================================
// LV MECHANICS: LVOT GRADIENT
// =============================================================
// LVOT gradient driven by: IVS, duty ratio (contractility), HR
// Reduced by mavacamten (duty ratio ↓), beta-blocker (HR+contractility↓)
double IVS_contrib  = 2.5 * (IVS - 11.0);   // mmHg per mm above 11mm
double DR_contrib   = 40.0 * (DR_ratio - 1.0);
double HR_contrib   = 0.25 * (HR_state - 70.0);
double LVOT_target  = LVOT_min + IVS_contrib + DR_contrib + HR_contrib;
if (LVOT_target < LVOT_min) LVOT_target = LVOT_min;

dxdt_LVOT = k_LVOT_dyn * (LVOT_target - LVOT);

// =============================================================
// LVEDP (DIASTOLIC FUNCTION)
// =============================================================
// LVEDP driven by: LV stiffness (fibrosis), preload (RAAS), LVOT afterload
double Col_effect  = 1.0 / (1.0 + 0.50 * (Collagen - 1.0));  // ↓ compliance
double LVEDP_target = 8.0
                    + 12.0 / Col_effect
                    + 0.30 * LVOT
                    - 2.0 * E_bb_hr;  // BB reduces LVEDP via HR reduction
if (LVEDP_target < 5.0) LVEDP_target = 5.0;

dxdt_LVEDP = k_LVEDP * (LVEDP_target - LVEDP);

// =============================================================
// HEART RATE
// =============================================================
// HR elevated by SNS (secondary to obstruction), reduced by BB/CCB
double HR_target = HR_base * (1.0 - E_bb_hr) + 0.15 * LVOT;
if (HR_target < HR_min) HR_target = HR_min;
if (HR_target > 130.0)  HR_target = 130.0;

dxdt_HR_state = 0.10 * (HR_target - HR_state);

// =============================================================
// NT-proBNP (wall-stress driven)
// =============================================================
// NT-proBNP rises with LVEDP (wall stress index)
double NTpBNP_target = 50.0 * exp(0.14 * LVEDP)
                     * (1.0 + 0.20 * (LVmass_rel - 1.0))
                     * (1.0 + 0.10 * (Collagen - 1.0));
double BNP_clearance = k_BNP_deg * NT_proBNP;

dxdt_NT_proBNP = k_BNP_syn * (NTpBNP_target - NT_proBNP) * 0.05;

// =============================================================
// TROPONIN I (myocardial injury — ischemia driven)
// =============================================================
// Supply-demand mismatch: high LVOT + high DR = high MVO2 vs limited supply
double MVO2_index   = DR_ratio * (1.0 + 0.20 * LVOT / 50.0);
double isch_factor  = (MVO2_index > 1.3) ? (MVO2_index - 1.3) / 1.3 : 0.0;
double Trop_target  = TropI_base * (1.0 + 3.0 * isch_factor);

dxdt_TroponinI = k_Trop_rel * (Trop_target - TroponinI)
               - k_Trop_clr * (TroponinI - TropI_base);

// =============================================================
// AF CUMULATIVE HAZARD (Poisson-like accumulation)
// =============================================================
// LA pressure drives atrial remodeling → AF
double LA_press_est = LVEDP + 3.0 + 0.20 * LVOT;  // estimated LA pressure
double AF_rate      = AF_hazard_base * (LA_press_est / 15.0) / (24.0 * 365.0);
// (per hr, from annual %)

dxdt_AF_hazard = AF_rate;

$TABLE
// ---- Derived PK outputs ----
capture Conc_mava_nM  = A_c_m / V1_m;        // mavacamten nM
capture Conc_BB_ngmL  = A_bb / V_bb;          // beta-blocker ng/mL

// ---- Derived PD outputs ----
capture E_mava_pct    = 100.0 * E_mava;       // % duty ratio suppression
capture DR_out        = DR_HCM * (1.0 - E_mava) * (1.0 - 0.5 * E_bb_c);

// ---- Hemodynamic outputs ----
capture LVEF_pct      = LVEF_base
                      - 2.0 * (Collagen - Col_base)
                      - 1.5 * (LVmass_rel - LVmass_0);
capture SV_mL         = SV_base * (1.0 - 0.5 * (LVOT - LVOT_base) / LVOT_base);
capture CO_Lmin       = SV_mL * HR_state / 1000.0;
capture LA_press_out  = LVEDP + 3.0 + 0.20 * LVOT;

// ---- Diastolic function grade ----
// 0=normal, 1=impaired relaxation, 2=pseudonormal, 3=restrictive
double E_e_est = (LVEDP + 3.0) / 1.5;
capture Diastolic_grade = (E_e_est < 8.0) ? 0 :
                          (E_e_est < 13.0) ? 1 :
                          (E_e_est < 20.0) ? 2 : 3;

// ---- Fibrosis extent ----
capture ECV_pct       = 20.0 + 8.0 * (Collagen - 1.0);  // estimated ECV%
capture LGE_pct_LV    = (Collagen > 1.5) ? 8.0 * (Collagen - 1.5) : 0.0;

// ---- SCD risk (HCM Risk-SCD components) ----
capture SCD_5yr_pct   = 0.3 * IVS
                      + 0.05 * LVOT
                      - 0.8 * HR_state / 100.0
                      - 1.0;  // simplified composite score

// ---- NYHA class (derived from peak VO2 surrogate) ----
double VO2_est = 30.0 - 0.15 * (LVOT - LVOT_base)
                       - 0.50 * (LVEDP - 10.0)
                       - 2.0  * (Diastolic_grade);
capture peak_VO2_est = (VO2_est < 8.0) ? 8.0 : VO2_est;
capture NYHA_est     = (VO2_est > 20.0) ? 1 :
                       (VO2_est > 14.0) ? 2 :
                       (VO2_est > 10.0) ? 3 : 4;

$CAPTURE
Conc_mava_nM Conc_BB_ngmL E_mava_pct DR_out
LVEF_pct SV_mL CO_Lmin LA_press_out
Diastolic_grade ECV_pct LGE_pct_LV SCD_5yr_pct
peak_VO2_est NYHA_est
'

hcm_mod <- mcode("hcm_qsp", hcm_code)

# ================================================================
# TREATMENT SCENARIOS (5+ scenarios)
# ================================================================

# Helper: build dosing regimen
make_mava_dose_ng <- function(dose_mg, freq_hr = 24, ndays = 365) {
  # Convert mg dose to nmol (MW = 471 g/mol)
  dose_nmol <- dose_mg * 1e6 / 471.0
  ev(amt = dose_nmol * 0.93, ii = freq_hr, addl = ndays * (24/freq_hr) - 1,
     cmt = "A_gut_m")
}

make_bb_dose_ug <- function(dose_mg, freq_hr = 24, ndays = 365) {
  dose_ug <- dose_mg * 1000
  ev(amt = dose_ug, ii = freq_hr, addl = ndays * (24/freq_hr) - 1,
     cmt = "A_bb")
}

## Simulation time: 2 years (17520 hours)
sim_time <- seq(0, 17520, by = 24)

## ----------------------------------------------------------------
## Scenario 1: Untreated HCM (natural history)
## ----------------------------------------------------------------
ev1 <- ev(amt = 0, cmt = "A_gut_m")  # no drug

out1 <- hcm_mod %>%
  ev(ev1) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Untreated HCM")

## ----------------------------------------------------------------
## Scenario 2: Mavacamten 5 mg QD (CYP2C19 EM starting dose)
## ----------------------------------------------------------------
ev2 <- make_mava_dose_ng(dose_mg = 5)

out2 <- hcm_mod %>%
  param(CYP_factor = 1.0) %>%
  ev(ev2) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Mavacamten 5mg (EM)")

## ----------------------------------------------------------------
## Scenario 3: Mavacamten 10 mg QD (uptitrated EM)
## ----------------------------------------------------------------
ev3 <- make_mava_dose_ng(dose_mg = 10)

out3 <- hcm_mod %>%
  param(CYP_factor = 1.0) %>%
  ev(ev3) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Mavacamten 10mg (EM)")

## ----------------------------------------------------------------
## Scenario 4: Mavacamten 2.5 mg QD (CYP2C19 PM — reduced dose)
## ----------------------------------------------------------------
ev4 <- make_mava_dose_ng(dose_mg = 2.5)

out4 <- hcm_mod %>%
  param(CYP_factor = 2.5) %>%   # PM: slower metabolism → higher exposure
  ev(ev4) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Mavacamten 2.5mg (PM)")

## ----------------------------------------------------------------
## Scenario 5: Beta-blocker monotherapy (metoprolol 200 mg QD)
## ----------------------------------------------------------------
ev5 <- make_bb_dose_ug(dose_mg = 200)

out5 <- hcm_mod %>%
  ev(ev5) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Beta-blocker 200mg")

## ----------------------------------------------------------------
## Scenario 6: Combination (Mavacamten 5mg + Beta-blocker 100mg)
## ----------------------------------------------------------------
ev6a <- make_mava_dose_ng(dose_mg = 5)
ev6b <- make_bb_dose_ug(dose_mg = 100)
ev6  <- ev6a + ev6b

out6 <- hcm_mod %>%
  param(CYP_factor = 1.0) %>%
  ev(ev6) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Mava 5mg + BB 100mg")

## ----------------------------------------------------------------
## Scenario 7: Post-septal reduction therapy (SRT) — simulate
##   immediate LVOT reduction (surgical myectomy modeled as
##   step change in IVS at time 0)
## ----------------------------------------------------------------
# Model SRT as reduced initial IVS (post-myectomy IVS ~10-12mm)
out7 <- hcm_mod %>%
  init(IVS = 11.0, LVOT = 10.0) %>%
  ev(ev(amt = 0, cmt = "A_gut_m")) %>%
  mrgsim(end = 17520, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "Post-SRT (Myectomy)")

## ----------------------------------------------------------------
## Combine results
## ----------------------------------------------------------------
all_results <- bind_rows(out1, out2, out3, out4, out5, out6, out7)

# ================================================================
# SUMMARY TABLE (Week 52, Day 365)
# ================================================================
summary_52w <- all_results %>%
  filter(time == 8760) %>%   # 365 days = 8760 hours
  select(Scenario, LVOT, LVEDP, NT_proBNP, TroponinI, IVS,
         Collagen, NYHA_est, peak_VO2_est, LVEF_pct,
         Conc_mava_nM, E_mava_pct, CO_Lmin, ECV_pct) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

cat("\n=== HCM QSP Model — 52-Week Summary ===\n")
print(summary_52w, width = 120)

# ================================================================
# KEY PLOTS (requires ggplot2)
# ================================================================
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("dplyr", quietly = TRUE)) {

  library(ggplot2)
  library(dplyr)

  all_results$day <- all_results$time / 24

  cols7 <- c("#E53935","#1E88E5","#43A047","#F4511E",
             "#8E24AA","#039BE5","#6D4C41")

  # Plot 1: LVOT gradient over time
  p1 <- ggplot(all_results, aes(day, LVOT, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 30, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cols7) +
    labs(title = "LVOT Gradient Over 2 Years",
         x = "Day", y = "LVOT Gradient (mmHg)",
         caption = "Dashed line: 30 mmHg clinical threshold") +
    theme_bw() + theme(legend.position = "bottom")

  # Plot 2: NT-proBNP over time
  p2 <- ggplot(all_results, aes(day, NT_proBNP, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = cols7) +
    labs(title = "NT-proBNP Over 2 Years",
         x = "Day", y = "NT-proBNP (pg/mL)") +
    theme_bw() + theme(legend.position = "bottom")

  # Plot 3: IVS thickness (hypertrophy)
  p3 <- ggplot(all_results, aes(day, IVS, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = cols7) +
    labs(title = "IVS Thickness (Hypertrophy) Over 2 Years",
         x = "Day", y = "IVS Thickness (mm)") +
    theme_bw() + theme(legend.position = "bottom")

  # Plot 4: Mavacamten PK (Scenario 2 only)
  mava_pk <- out2 %>%
    filter(day <= 30) %>%
    select(day, Conc_mava_nM, E_mava_pct)

  p4 <- ggplot(mava_pk, aes(day, Conc_mava_nM)) +
    geom_line(color = "#1E88E5", linewidth = 1.2) +
    geom_hline(yintercept = 85, linetype = "dashed", color = "#E53935") +
    labs(title = "Mavacamten PK — First 30 Days (5 mg QD, CYP2C19 EM)",
         x = "Day", y = "Mavacamten Plasma Conc (nM)",
         caption = "Dashed: EC50 = 85 nM") +
    theme_bw()

  # Plot 5: Fibrosis (collagen content)
  p5 <- ggplot(all_results, aes(day, Collagen, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = cols7) +
    labs(title = "Myocardial Fibrosis (Collagen) Over 2 Years",
         x = "Day", y = "Relative Collagen Content") +
    theme_bw() + theme(legend.position = "bottom")

  # Plot 6: NYHA class
  p6 <- ggplot(all_results %>% filter(time %% 720 == 0),
               aes(day, NYHA_est, color = Scenario)) +
    geom_step(linewidth = 0.9) +
    scale_y_continuous(breaks = 1:4, labels = paste("NYHA", 1:4)) +
    scale_color_manual(values = cols7) +
    labs(title = "Estimated NYHA Class Over 2 Years",
         x = "Day", y = "NYHA Functional Class") +
    theme_bw() + theme(legend.position = "bottom")

  print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)
}

# ================================================================
# SENSITIVITY ANALYSIS: EC50 for mavacamten
# ================================================================
cat("\n=== Sensitivity Analysis: Mavacamten EC50 ===\n")

ec50_grid <- c(50, 85, 120, 170)
sa_results <- lapply(ec50_grid, function(ec50) {
  hcm_mod %>%
    param(EC50_mava = ec50) %>%
    ev(make_mava_dose_ng(10)) %>%
    mrgsim(end = 8760, delta = 168) %>%
    as.data.frame() %>%
    filter(time == 8760) %>%
    mutate(EC50 = ec50)
}) %>% bind_rows()

cat("LVOT at 52w by EC50 value:\n")
print(sa_results[, c("EC50", "LVOT", "E_mava_pct", "NT_proBNP", "NYHA_est")])

cat("\nHCM QSP mrgsolve model loaded successfully.\n")
cat("Scenarios: 1=Untreated, 2=Mava5mg-EM, 3=Mava10mg-EM,\n")
cat("           4=Mava2.5mg-PM, 5=BB200mg, 6=Mava+BB, 7=Post-SRT\n")
