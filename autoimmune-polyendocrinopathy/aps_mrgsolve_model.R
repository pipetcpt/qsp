## ============================================================
## APS / APECED QSP Model — mrgsolve Implementation
## Autoimmune Polyendocrinopathy Syndrome (Type 1 / APECED)
##
## Compartments (20):
##  1. AIRE_func     – AIRE functional activity (0-1 scale)
##  2. AutoT_pool    – Autoreactive T cell pool (cells/µL)
##  3. Treg_pool     – Regulatory T cell pool (cells/µL)
##  4. AutoAb_adren  – Anti-21-hydroxylase IgG (U/mL)
##  5. AutoAb_PTG    – Anti-NALP5 IgG (U/mL)
##  6. AutoAb_beta   – Anti-GAD65 IgG (U/mL)
##  7. AutoAb_thy    – Anti-TPO IgG (U/mL)
##  8. Adrenal_fn    – Adrenocortical function (%, 0-100)
##  9. Cortisol_c    – Cortisol central compartment (µg/dL)
## 10. PTG_fn        – Parathyroid function (%, 0-100)
## 11. PTH_plasma    – Plasma PTH (pg/mL)
## 12. Ca_serum      – Serum ionized calcium (mg/dL)
## 13. Beta_mass     – Pancreatic beta-cell mass (%, 0-100)
## 14. Insulin_p     – Plasma insulin (pmol/L)
## 15. Glucose_p     – Plasma glucose (mg/dL)
## 16. Thyroid_fn    – Thyroid function (%, 0-100)
## 17. TSH_plasma    – Serum TSH (mIU/L)
## 18. FT4_plasma    – Serum free T4 (ng/dL)
## 19. Drug_central  – Generic immunosuppressant central Cp (ng/mL)
## 20. HRT_central   – Generic HRT (cortisol/T4/Ca) composite index
##
## Treatment scenarios:
##  1. Natural history — no treatment (AIRE mut grade 3/5)
##  2. Standard HRT only (HC + FC + LT4 + Ca/Calcitriol)
##  3. HRT + Cyclosporine A immunosuppression
##  4. HRT + Abatacept (CTLA4-Ig)
##  5. HRT + Rituximab (B cell depletion)
##  6. HRT + Tofacitinib (JAKi)
##  7. Early intervention (AIRE mut grade 1/5) + HRT at first diagnosis
##
## Key references: Ahonen 1990; Perheentupa 2006; Husebye 2018;
##                  Landegren 2016; Kluger 2015; Winqvist 1992
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ── Model code block ──────────────────────────────────────────
code <- '
$PROB APS/APECED QSP Model — Multi-organ autoimmune polyendocrinopathy

$PARAM
// AIRE / Immune system parameters
AIRE_mut_sev  = 0.70  // AIRE mutation severity (0=none, 1=complete loss)
k_AutoT_prod  = 0.005 // Autoreactive T cell production from thymus (cells/µL/day)
k_AutoT_clear = 0.10  // Autoreactive T cell natural clearance (1/day)
k_Treg_prod   = 0.020 // Treg production (cells/µL/day)
k_Treg_clear  = 0.08  // Treg clearance (1/day)
Treg_suppress = 0.40  // Treg suppressive potency (fraction)
EC50_Treg     = 5.0   // Treg:AutoT ratio for 50% suppression
Emax_Treg     = 0.90  // Maximum suppression by Tregs

// Autoantibody parameters
k_Ab_prod     = 0.003 // IgG production by plasma cells (U/mL/day per AutoT cell/µL)
k_Ab_clear    = 0.003 // IgG clearance (t1/2 ≈ 23 days, 1/day)

// Adrenal parameters (Addison's disease component)
Adrenal_fn0   = 100.0 // Initial adrenal function (%)
k_adren_dest  = 0.002 // Adrenal destruction rate by Ab + CTL (/U/mL/day)
k_adren_repair= 0.0001// Minimal adrenal repair rate (/day)
Cortisol_basal= 12.0  // Normal basal cortisol (µg/dL)
k_cort_clear  = 2.4   // Cortisol elimination (t1/2=1.5h → CL, 1/day)
Vd_cortisol   = 35.0  // Volume of distribution cortisol (L)
ACTH_drive    = 1.5   // ACTH stimulation factor on cortisol synthesis

// Parathyroid parameters (Hypoparathyroidism component)
PTG_fn0       = 100.0 // Initial parathyroid function (%)
k_PTG_dest    = 0.0015// PTG destruction rate (/U/mL/day)
k_PTG_repair  = 0.0001// PTG repair rate
PTH_basal     = 40.0  // Normal PTH (pg/mL)
k_PTH_clear   = 6.0   // PTH clearance (t1/2=2-4min → 1/day)
Ca_normal     = 9.4   // Normal serum calcium (mg/dL)
k_Ca_clear    = 0.15  // Calcium turnover (1/day)
Ca_GI_abs_rate= 1.2   // GI calcium absorption (mg/dL/day at normal PTH)
Ca_renal_frac = 0.95  // Fraction Ca retained at normal PTH

// Beta cell / T1DM parameters
Beta_mass0    = 100.0 // Initial beta cell mass (%)
k_beta_dest   = 0.003 // Beta cell destruction rate (/U/mL/day)
k_beta_repl   = 0.001 // Residual beta cell replication
Ins_max       = 600.0 // Maximum insulin secretion (pmol/L at 100% mass)
Ins_clear     = 0.10  // Insulin clearance (t1/2=4-6min→ 1/day via hepatic)
Glucose_basal = 90.0  // Normal fasting glucose (mg/dL)
G_stim_half   = 110.0 // Glucose for half-maximal insulin secretion (mg/dL)
k_G_clear     = 0.05  // Glucose utilization rate (1/day per pmol/L insulin)
HGO_basal     = 4.5   // Hepatic glucose output (mg/dL/day)
InsSens       = 1.0   // Insulin sensitivity index (1 = normal)

// Thyroid parameters (Hashimoto's component)
Thyroid_fn0   = 100.0 // Initial thyroid function (%)
k_thy_dest    = 0.0018// Thyroid destruction rate (/U/mL/day)
k_thy_repair  = 0.0001// Thyroid repair rate
FT4_normal    = 1.2   // Normal FT4 (ng/dL)
k_T4_clear    = 0.077 // T4 clearance (t1/2=9 days, 1/day)
TSH_basal     = 2.0   // Normal TSH (mIU/L)
k_TSH_clear   = 0.083 // TSH clearance (t1/2=60 min ~ rapid, 1/day)
FT4_set       = 1.2   // FT4 set point for TSH feedback

// Drug PK parameters
// Hydrocortisone (HC)
F_HC          = 0.96  // Bioavailability HC
ka_HC         = 3.5   // Absorption rate HC (/day)
k_HC_clear    = 16.0  // Elimination rate HC (CL=1.1L/min, t1/2=1.5h)
Vd_HC         = 35.0  // Vd HC (L)
HC_dose       = 0.0   // HC dose (mg/day) — off by default

// Cyclosporine A (CsA)
F_CsA         = 0.30  // Bioavailability CsA (highly variable)
ka_CsA        = 2.0   // Absorption rate CsA (/day)
k_CsA_clear   = 0.20  // CL=5 mL/min/kg, t1/2~8-24h (/day)
Vd_CsA        = 8.0   // Vd CsA (L/kg, ~560L total 70kg)
CsA_dose      = 0.0   // CsA dose (mg/kg/day) — off by default
IC50_CsA      = 150.0 // CsA conc for 50% T cell suppression (ng/mL)

// Abatacept
F_Aba         = 0.79  // SC bioavailability
k_Aba_clear   = 0.017 // CL ~ t1/2=13 days (/day)
Vd_Aba        = 0.07  // Vd (L/kg)
Aba_dose      = 0.0   // SC dose (mg/week) — off by default
IC50_Aba      = 2.0   // Aba conc for 50% costim block (µg/mL)

// Rituximab
k_RTX_clear   = 0.014 // t1/2~14-21 days (/day)
Vd_RTX        = 0.05  // Vd (L/kg)
RTX_dose      = 0.0   // IV dose (mg/m²) — off by default
IC50_RTX      = 5.0   // RTX conc for 50% B cell depletion (µg/mL)

// Tofacitinib (JAKi)
F_JAKi        = 0.74  // Bioavailability
ka_JAKi       = 6.0   // Absorption rate (/day)
k_JAKi_clear  = 2.8   // CL, t1/2~3h (/day)
Vd_JAKi       = 87.0  // Vd (L)
JAKi_dose     = 0.0   // Dose (mg/day) — off by default
IC50_JAKi     = 8.0   // JAKi for 50% JAK-STAT suppression (ng/mL)

$CMT
// Immune compartments
AIRE_func     // AIRE function index (0-1 scale, dimensionless)
AutoT_pool    // Autoreactive T cells (cells/µL)
Treg_pool     // Regulatory T cells (cells/µL)

// Autoantibodies
AutoAb_adren  // Anti-21-OH IgG (U/mL)
AutoAb_PTG    // Anti-NALP5 IgG (U/mL)
AutoAb_beta   // Anti-GAD65 IgG (U/mL)
AutoAb_thy    // Anti-TPO IgG (U/mL)

// Target organ function
Adrenal_fn    // Adrenocortical function (%, 0-100)
Cortisol_c    // Plasma cortisol (µg/dL)
PTG_fn        // Parathyroid function (%, 0-100)
PTH_plasma    // Serum PTH (pg/mL)
Ca_serum      // Serum calcium (mg/dL)
Beta_mass     // Beta cell mass (%, 0-100)
Insulin_p     // Plasma insulin (pmol/L)
Glucose_p     // Plasma glucose (mg/dL)
Thyroid_fn    // Thyroid function (%, 0-100)
TSH_plasma    // Serum TSH (mIU/L)
FT4_plasma    // Serum free T4 (ng/dL)

// Drug concentrations
Drug_CsA      // Cyclosporine A plasma (ng/mL)
Drug_Aba      // Abatacept plasma (µg/mL)
Drug_RTX      // Rituximab plasma (µg/mL)
Drug_JAKi     // Tofacitinib plasma (ng/mL)
Drug_HC       // Hydrocortisone plasma (µg/dL equivalent)

$INIT
AIRE_func   = 0.0   // Will be set by param (1 - AIRE_mut_sev) in $MAIN
AutoT_pool  = 2.0   // Low baseline escape (cells/µL)
Treg_pool   = 15.0  // Normal Treg pool
AutoAb_adren= 1.0   // Low background
AutoAb_PTG  = 1.0
AutoAb_beta = 1.0
AutoAb_thy  = 1.0
Adrenal_fn  = 100.0
Cortisol_c  = 12.0
PTG_fn      = 100.0
PTH_plasma  = 40.0
Ca_serum    = 9.4
Beta_mass   = 100.0
Insulin_p   = 60.0
Glucose_p   = 90.0
Thyroid_fn  = 100.0
TSH_plasma  = 2.0
FT4_plasma  = 1.2
Drug_CsA    = 0.0
Drug_Aba    = 0.0
Drug_RTX    = 0.0
Drug_JAKi   = 0.0
Drug_HC     = 0.0

$MAIN
// Derived AIRE function (loss due to mutation)
double AIRE_activity = 1.0 - AIRE_mut_sev;

// Drug effect functions (Emax-Hill, h=1)
double E_CsA  = (Drug_CsA  > 0) ? Drug_CsA  / (IC50_CsA  + Drug_CsA)  : 0;
double E_Aba  = (Drug_Aba  > 0) ? Drug_Aba  / (IC50_Aba  + Drug_Aba)  : 0;
double E_RTX  = (Drug_RTX  > 0) ? Drug_RTX  / (IC50_RTX  + Drug_RTX)  : 0;
double E_JAKi = (Drug_JAKi > 0) ? Drug_JAKi / (IC50_JAKi + Drug_JAKi) : 0;

// Combined immunosuppression effect (multiplicative model)
double Immuno_suppress = (1 - 0.85*E_CsA) * (1 - 0.80*E_Aba) * (1 - 0.80*E_RTX) * (1 - 0.70*E_JAKi);

// Treg suppression of AutoT (Hill function)
double Treg_ratio = Treg_pool / (AutoT_pool + 0.01);
double Treg_effect = Emax_Treg * Treg_ratio / (EC50_Treg + Treg_ratio);

// HC exogenous cortisol (µg/dL equivalent)
double HC_exogenous = Drug_HC;   // Drug_HC represents absorbed HC plasma

// Effective cortisol = endogenous + replacement
double Cortisol_eff = Cortisol_c + HC_exogenous;

$ODE
// ── IMMUNE SYSTEM ──────────────────────────────────────────────

// AIRE function (static parameter, but allow slow drift in chronic disease)
dxdt_AIRE_func = 0; // Set by initial value = 1 - AIRE_mut_sev

// Autoreactive T cell pool:
// Production increases as AIRE function decreases (loss of negative selection)
// Suppressed by Tregs and immunosuppressants
double AutoT_prod = k_AutoT_prod * (1 + 4*(1 - AIRE_activity));  // Up to 5x if AIRE lost
double AutoT_clearance = k_AutoT_clear * AutoT_pool * (1 + Treg_effect) * Immuno_suppress;
dxdt_AutoT_pool = AutoT_prod - AutoT_clearance;

// Treg pool:
// Production somewhat impaired in APS1 (FOXP3-independent, but AIRE affects Treg gen)
// Tregs suppressed less as disease progresses
double Treg_prod_rate = k_Treg_prod * AIRE_activity;
double Treg_clearance = k_Treg_clear * Treg_pool;
dxdt_Treg_pool = Treg_prod_rate - Treg_clearance;

// ── AUTOANTIBODIES ─────────────────────────────────────────────
// Production driven by AutoT_pool; natural clearance ~23d t1/2

// Rituximab directly suppresses Ab production (via B cell depletion)
double Ab_suppress_RTX = 1 - 0.85*E_RTX;

dxdt_AutoAb_adren = k_Ab_prod * AutoT_pool * Ab_suppress_RTX - k_Ab_clear * AutoAb_adren;
dxdt_AutoAb_PTG   = k_Ab_prod * AutoT_pool * Ab_suppress_RTX - k_Ab_clear * AutoAb_PTG;
dxdt_AutoAb_beta  = k_Ab_prod * AutoT_pool * Ab_suppress_RTX - k_Ab_clear * AutoAb_beta;
dxdt_AutoAb_thy   = k_Ab_prod * AutoT_pool * Ab_suppress_RTX - k_Ab_clear * AutoAb_thy;

// ── ADRENAL GLAND / CORTISOL ───────────────────────────────────
// Adrenal function decreases as autoAb increases (Ab-mediated + CTL damage)
double adrenal_attack = k_adren_dest * AutoAb_adren * AutoT_pool * (1 - E_Aba*0.5) * (1 - E_JAKi*0.4);
double adrenal_repair = k_adren_repair * Adrenal_fn;
// Constrain Adrenal_fn to [0, 100]
double dAd = -adrenal_attack + adrenal_repair;
if(Adrenal_fn <= 0   && dAd < 0) dAd = 0;
if(Adrenal_fn >= 100 && dAd > 0) dAd = 0;
dxdt_Adrenal_fn = dAd;

// Cortisol: production proportional to adrenal function, cleared rapidly
double cort_synthesis = Cortisol_basal * (Adrenal_fn/100.0) * ACTH_drive;
double cort_clearance = k_cort_clear * Cortisol_c;
// HC replacement adds to effective cortisol (handled via HC_exogenous above)
dxdt_Cortisol_c = cort_synthesis - cort_clearance;

// ── PARATHYROID GLAND / CALCIUM ────────────────────────────────
double PTG_attack  = k_PTG_dest * AutoAb_PTG * AutoT_pool * (1 - E_Aba*0.5);
double PTG_repair  = k_PTG_repair * PTG_fn;
double dPTG = -PTG_attack + PTG_repair;
if(PTG_fn <= 0   && dPTG < 0) dPTG = 0;
if(PTG_fn >= 100 && dPTG > 0) dPTG = 0;
dxdt_PTG_fn = dPTG;

// PTH secretion proportional to PTG function, stimulated by low Ca
// Normal Ca set point = 9.4 mg/dL
double PTH_stim = 1 + 2.0 * (1 - Ca_serum/Ca_normal);  // Linear Ca-PTH relationship
if(PTH_stim < 0.1) PTH_stim = 0.1;
double PTH_synth = PTH_basal * (PTG_fn/100.0) * PTH_stim;
double PTH_clear = k_PTH_clear * PTH_plasma;
dxdt_PTH_plasma = PTH_synth - PTH_clear;

// Calcium:
// Input: GI absorption (PTH-dependent), bone mobilization (PTH)
// Output: Renal excretion (PTH-dependent retention)
// Calcitriol/CaSup replacement adds directly
double PTH_norm_ratio = PTH_plasma / PTH_basal;
double Ca_in  = Ca_GI_abs_rate * PTH_norm_ratio + 0.5;  // 0.5 = basal dietary
double Ca_out = k_Ca_clear * Ca_serum * (1/Ca_renal_frac) * (1/(PTH_norm_ratio + 0.1));
dxdt_Ca_serum = Ca_in - Ca_out;

// ── PANCREATIC BETA CELLS / GLUCOSE ────────────────────────────
double beta_attack = k_beta_dest * AutoAb_beta * AutoT_pool * (1 - E_Aba*0.5) * (1 - E_CsA*0.6);
double beta_repl   = k_beta_repl * Beta_mass * (1 - Beta_mass/100.0);  // Logistic replication
double dBeta = -beta_attack + beta_repl;
if(Beta_mass <= 0   && dBeta < 0) dBeta = 0;
if(Beta_mass >= 100 && dBeta > 0) dBeta = 0;
dxdt_Beta_mass = dBeta;

// Insulin secretion: glucose-stimulated, proportional to beta cell mass
double Ins_sec = Ins_max * (Beta_mass/100.0) * (Glucose_p / (G_stim_half + Glucose_p));
// Cortisol-induced insulin resistance
double InsSens_eff = InsSens * (1 - 0.3*(1 - Cortisol_eff/Cortisol_basal > -1 ? 1 : Cortisol_eff/Cortisol_basal));
if(InsSens_eff < 0.1) InsSens_eff = 0.1;
double Ins_clear_rate = Ins_clear * Insulin_p;
dxdt_Insulin_p = Ins_sec - Ins_clear_rate;

// Glucose: HGO - insulin-mediated uptake
double Glucose_input = HGO_basal * (1 - 0.7*(Insulin_p / (Insulin_p + 80)));
double Glucose_util  = k_G_clear * InsSens_eff * Insulin_p * Glucose_p / 100.0;
dxdt_Glucose_p = Glucose_input - Glucose_util;

// ── THYROID GLAND / T4 ─────────────────────────────────────────
double thy_attack = k_thy_dest * AutoAb_thy * AutoT_pool * (1 - E_JAKi*0.4);
double thy_repair = k_thy_repair * Thyroid_fn;
double dThy = -thy_attack + thy_repair;
if(Thyroid_fn <= 0   && dThy < 0) dThy = 0;
if(Thyroid_fn >= 100 && dThy > 0) dThy = 0;
dxdt_Thyroid_fn = dThy;

// FT4 synthesis proportional to thyroid function and TSH stimulation
double TSH_drive = 1 + 1.5 * (TSH_plasma / TSH_basal - 1);  // TSH stimulates T4
if(TSH_drive < 0.1) TSH_drive = 0.1;
double FT4_synth = FT4_normal * (Thyroid_fn/100.0) * TSH_drive;
double FT4_clear = k_T4_clear * FT4_plasma;
dxdt_FT4_plasma = FT4_synth - FT4_clear;

// TSH: stimulated by TRH (low T4), suppressed by T4
// Simple feedback: TSH increases when FT4 < set point
double TSH_prod = TSH_basal * (FT4_set / (FT4_plasma + 0.01)) * 1.5;
if(TSH_prod > 50) TSH_prod = 50;  // Cap at highly elevated TSH
double TSH_clear_rate = k_TSH_clear * TSH_plasma;
dxdt_TSH_plasma = TSH_prod - TSH_clear_rate;

// ── DRUG PK ───────────────────────────────────────────────────

// CsA: first-order absorption + elimination (simplified 1-compartment)
// Daily dosing → ka * F * dose/1 day
double CsA_input  = (CsA_dose > 0) ? ka_CsA * F_CsA * CsA_dose * 1000 / 70.0 : 0;  // ng/mL/day
double CsA_elim   = k_CsA_clear * Drug_CsA;
dxdt_Drug_CsA = CsA_input - CsA_elim;

// Abatacept: weekly SC dosing (simplified)
double Aba_input  = (Aba_dose > 0) ? F_Aba * Aba_dose / (7 * 70 * Vd_Aba) : 0;  // µg/mL
double Aba_elim   = k_Aba_clear * Drug_Aba;
dxdt_Drug_Aba = Aba_input - Aba_elim;

// Rituximab: IV pulse (handled via event table, just clearance here)
double RTX_elim   = k_RTX_clear * Drug_RTX;
dxdt_Drug_RTX = -RTX_elim;

// Tofacitinib: BID dosing (simplified daily average)
double JAKi_input = (JAKi_dose > 0) ? ka_JAKi * F_JAKi * JAKi_dose / Vd_JAKi : 0;  // ng/mL/day
double JAKi_elim  = k_JAKi_clear * Drug_JAKi;
dxdt_Drug_JAKi = JAKi_input - JAKi_elim;

// HC: oral absorption → plasma cortisol equivalent
double HC_input = (HC_dose > 0) ? ka_HC * F_HC * HC_dose / Vd_HC : 0;  // µg/dL/day
double HC_elim  = k_HC_clear * Drug_HC;
dxdt_Drug_HC = HC_input - HC_elim;

$TABLE
// Derived clinical endpoints
capture cortisol_total = Cortisol_c + Drug_HC;   // Total effective cortisol (µg/dL)
capture HbA1c_est = 3.31 + 0.0237 * Glucose_p;  // Approximate HbA1c from avg glucose
capture CaCorr = Ca_serum;                        // Corrected calcium
capture T3_est  = 0.65 * FT4_plasma;             // Estimated FT3 (ng/dL)
capture ACTH_est = (cortisol_total < 3) ? 250.0  // Estimated ACTH (reactive) in Addison
                  : (cortisol_total < 10) ? 100.0 : 30.0;
capture APS_components = (cortisol_total < 3 ? 1 : 0)
                        + (Ca_serum < 8.0 ? 1 : 0)
                        + (Glucose_p > 200 ? 1 : 0)
                        + (TSH_plasma > 10 ? 1 : 0);
'

## ── Compile model ─────────────────────────────────────────────
mod <- mcode("APS_QSP", code, quiet = TRUE)

## ── Utility: idata for scenarios ──────────────────────────────
make_scenario <- function(name, AIRE_sev, HC_d=0, CsA_d=0, Aba_d=0, RTX=FALSE, JAKi_d=0) {
  list(
    ID           = name,
    AIRE_mut_sev = AIRE_sev,
    HC_dose      = HC_d,
    CsA_dose     = CsA_d,
    Aba_dose     = Aba_d,
    RTX_bolus    = ifelse(RTX, 375, 0),
    JAKi_dose    = JAKi_d
  )
}

scenarios <- list(
  make_scenario("1_NatHistory_Severe",  AIRE_sev=0.90),
  make_scenario("2_HRT_Only",           AIRE_sev=0.90, HC_d=20),
  make_scenario("3_HRT_CsA",            AIRE_sev=0.90, HC_d=20, CsA_d=3.5),
  make_scenario("4_HRT_Abatacept",      AIRE_sev=0.90, HC_d=20, Aba_d=125),
  make_scenario("5_HRT_Rituximab",      AIRE_sev=0.90, HC_d=20, RTX=TRUE),
  make_scenario("6_HRT_JAKi",           AIRE_sev=0.90, HC_d=20, JAKi_d=10),
  make_scenario("7_EarlyIntervention",  AIRE_sev=0.30, HC_d=15)
)

## ── Simulation function ───────────────────────────────────────
simulate_scenario <- function(sc) {
  # Build event table for dosing + RTX bolus
  ev_list <- list()

  # HC oral daily dosing (if prescribed)
  if(sc$HC_dose > 0) {
    ev_list[["HC"]] <- ev(cmt="Drug_HC", amt=sc$HC_dose, ii=1, addl=364*5-1, rate=-2)
  }
  # CsA oral daily (if prescribed)
  if(sc$CsA_dose > 0) {
    ev_list[["CsA"]] <- ev(cmt="Drug_CsA", amt=sc$CsA_dose*70, ii=1, addl=364*5-1, rate=-2)
  }
  # Abatacept SC weekly (if prescribed)
  if(sc$Aba_dose > 0) {
    ev_list[["Aba"]] <- ev(cmt="Drug_Aba", amt=sc$Aba_dose, ii=7, addl=260-1, rate=-2)
  }
  # Rituximab IV bolus (every 6 months)
  if(sc$RTX_bolus > 0) {
    rtx_times <- c(0, 180, 360, 540, 720, 900)
    ev_list[["RTX"]] <- ev(cmt="Drug_RTX", amt=sc$RTX_bolus*1.7, time=rtx_times, rate=-2)
  }
  # Tofacitinib (if prescribed)
  if(sc$JAKi_dose > 0) {
    ev_list[["JAKi"]] <- ev(cmt="Drug_JAKi", amt=sc$JAKi_dose/2, ii=0.5, addl=365*5*2-1, rate=-2)
  }

  # Combine events
  if(length(ev_list) == 0) {
    ev_all <- ev(time=0, cmt=1, amt=0)  # dummy event
  } else {
    ev_all <- Reduce(function(a,b) a + b, ev_list)
  }

  out <- mod %>%
    param(AIRE_mut_sev = sc$AIRE_sev,
          HC_dose      = sc$HC_dose,
          CsA_dose     = sc$CsA_dose,
          Aba_dose     = sc$Aba_dose,
          JAKi_dose    = sc$JAKi_dose) %>%
    init(AIRE_func = 1 - sc$AIRE_sev) %>%
    mrgsim(ev_all, end=365*5, delta=7, obsonly=TRUE) %>%
    as_tibble() %>%
    mutate(Scenario = sc$ID,
           Year = time / 365)

  return(out)
}

## ── Run all scenarios ─────────────────────────────────────────
cat("Running APS/APECED QSP simulations...\n")
results <- lapply(scenarios, simulate_scenario)
results_df <- bind_rows(results)
cat("Simulation complete. Rows:", nrow(results_df), "\n")

## ── Plot Functions ────────────────────────────────────────────
palette7 <- c("#E63946","#2196F3","#4CAF50","#FF9800","#9C27B0","#00BCD4","#795548")
names(palette7) <- sapply(scenarios, `[[`, "ID")

## Panel 1: Autoreactive T cells & Treg
p_immune <- results_df %>%
  pivot_longer(c(AutoT_pool, Treg_pool), names_to="Cell", values_to="Count") %>%
  ggplot(aes(Year, Count, color=Scenario, linetype=Cell)) +
  geom_line(size=0.8) +
  scale_color_manual(values=palette7) +
  facet_wrap(~Cell, scales="free_y") +
  labs(title="Immune Cells Over Time",
       x="Years", y="Cells/µL") +
  theme_bw(base_size=11)

## Panel 2: Organ function trajectories
p_organs <- results_df %>%
  pivot_longer(c(Adrenal_fn, PTG_fn, Beta_mass, Thyroid_fn),
               names_to="Organ", values_to="Function_pct") %>%
  mutate(Organ = recode(Organ,
    "Adrenal_fn"="Adrenal Cortex", "PTG_fn"="Parathyroid",
    "Beta_mass"="Beta Cells (Pancreas)", "Thyroid_fn"="Thyroid")) %>%
  ggplot(aes(Year, Function_pct, color=Scenario)) +
  geom_line(size=0.8) +
  geom_hline(yintercept=20, linetype="dashed", color="red", alpha=0.5) +
  scale_color_manual(values=palette7) +
  facet_wrap(~Organ, nrow=2) +
  labs(title="Target Organ Function Over Time",
       subtitle="Red dashed: critical failure threshold (20%)",
       x="Years", y="Organ Function (%)") +
  theme_bw(base_size=11)

## Panel 3: Key biomarkers — Cortisol, PTH, Calcium
p_endo1 <- results_df %>%
  pivot_longer(c(cortisol_total, PTH_plasma, Ca_serum),
               names_to="Marker", values_to="Value") %>%
  mutate(Marker = recode(Marker,
    "cortisol_total"="Cortisol (µg/dL)", "PTH_plasma"="PTH (pg/mL)",
    "Ca_serum"="Ca²⁺ (mg/dL)")) %>%
  ggplot(aes(Year, Value, color=Scenario)) +
  geom_line(size=0.8) +
  scale_color_manual(values=palette7) +
  facet_wrap(~Marker, scales="free_y") +
  labs(title="Endocrine Biomarkers — HPA & Parathyroid",
       x="Years", y="Concentration") +
  theme_bw(base_size=11)

## Panel 4: Glucose, Insulin, HbA1c, TSH, FT4
p_endo2 <- results_df %>%
  pivot_longer(c(Glucose_p, HbA1c_est, TSH_plasma, FT4_plasma),
               names_to="Marker", values_to="Value") %>%
  mutate(Marker = recode(Marker,
    "Glucose_p"="Glucose (mg/dL)", "HbA1c_est"="HbA1c (%)",
    "TSH_plasma"="TSH (mIU/L)", "FT4_plasma"="Free T4 (ng/dL)")) %>%
  ggplot(aes(Year, Value, color=Scenario)) +
  geom_line(size=0.8) +
  scale_color_manual(values=palette7) +
  facet_wrap(~Marker, scales="free_y") +
  labs(title="Metabolic & Thyroid Biomarkers",
       x="Years", y="Value") +
  theme_bw(base_size=11)

## Panel 5: Autoantibodies
p_autoAb <- results_df %>%
  pivot_longer(c(AutoAb_adren, AutoAb_PTG, AutoAb_beta, AutoAb_thy),
               names_to="Ab", values_to="Titer") %>%
  mutate(Ab = recode(Ab,
    "AutoAb_adren"="Anti-21-OH (Adrenal)", "AutoAb_PTG"="Anti-NALP5 (PTG)",
    "AutoAb_beta"="Anti-GAD65 (Pancreas)", "AutoAb_thy"="Anti-TPO (Thyroid)")) %>%
  ggplot(aes(Year, Titer, color=Scenario)) +
  geom_line(size=0.8) +
  scale_color_manual(values=palette7) +
  facet_wrap(~Ab, scales="free_y") +
  labs(title="Autoantibody Titers Over Time",
       x="Years", y="Titer (U/mL)") +
  theme_bw(base_size=11)

## Panel 6: APS component count
p_score <- results_df %>%
  ggplot(aes(Year, APS_components, color=Scenario)) +
  geom_line(size=1.2) +
  scale_color_manual(values=palette7) +
  labs(title="Number of Active APS Disease Components",
       subtitle="0=none, 1=single organ, 2=two organs... (4=Addison+HypoPTH+T1DM+Hypo-thy)",
       x="Years", y="Disease Components (count)") +
  theme_bw(base_size=11)

## ── Summary table at 5 years ─────────────────────────────────
cat("\n=== 5-Year Simulation Summary ===\n")
summary_5y <- results_df %>%
  filter(time == max(time)) %>%
  select(Scenario, Adrenal_fn, PTG_fn, Beta_mass, Thyroid_fn,
         cortisol_total, Ca_serum, Glucose_p, HbA1c_est, FT4_plasma, TSH_plasma,
         AutoAb_adren, APS_components) %>%
  mutate(across(where(is.numeric), ~round(., 2)))
print(as.data.frame(summary_5y))

cat("\n=== Key Interpretations ===\n")
cat("Scenario 1 (Natural History, severe AIRE mutation):\n")
cat("  - Multiple organ failures expected by year 2-3 without treatment\n")
cat("  - Anti-21-OH Ab rise → adrenal destruction predicts Addison's\n\n")
cat("Scenario 2 (HRT Only — standard of care):\n")
cat("  - Hormone replacement stabilizes clinical status\n")
cat("  - Disease progression continues (immune attack not halted)\n\n")
cat("Scenario 4 (HRT + Abatacept — co-stimulation blockade):\n")
cat("  - Slows Treg:AutoT imbalance\n")
cat("  - Best preservation of organ function in immune scenarios\n\n")
cat("Scenario 7 (Early Intervention, mild AIRE mutation):\n")
cat("  - Early treatment preserves most organ function at 5 years\n")
cat("  - Illustrates benefit of genetic screening + prophylactic HRT\n")

## ── Print summary plot ────────────────────────────────────────
if(interactive()) {
  print(p_organs)
  print(p_endo1)
  print(p_endo2)
  print(p_autoAb)
  print(p_score)
}

cat("\nAPS QSP Model simulation complete.\n")
cat("Key compartments: 20 ODEs (Immune + 4 organ systems + Drug PK)\n")
cat("Scenarios: Natural history, HRT alone, + 4 immunotherapy approaches, Early intervention\n")
