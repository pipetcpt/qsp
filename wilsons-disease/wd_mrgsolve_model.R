## ============================================================
## Wilson's Disease (WD) QSP Model — mrgsolve Implementation
## File: wd_mrgsolve_model.R
## Author: QSP Library (CCR)
## Date: 2026-06-24
##
## Model Overview:
##   24-compartment ODE system covering:
##   - Copper kinetics (GI, hepatocyte, systemic NCBC,
##     brain, kidney, cornea)
##   - Ceruloplasmin synthesis & secretion
##   - Liver pathology (ROS, ALT, fibrosis)
##   - Neurological degeneration
##   - Drug PK: D-Penicillamine (DPA), Zinc, Trientine, ALXN1840
##
## Key References:
##   - Bandmann et al. Lancet Neurol 2015
##   - Członkowska et al. Nat Rev Dis Primers 2018
##   - Schilsky et al. AASLD Practice Guidance 2023
##   - ATLAS trial (ALXN1840): Schilsky NEJM 2022
##   - Weiss et al. J Hepatol 2011 (pharmacokinetics)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## Model Code
## ============================================================
wd_code <- '
$PROB
Wilson Disease QSP Model — Copper Metabolism & Pharmacological Interventions

$PARAM @annotated
// ---- Copper Dietary Input ----
Cu_intake  : 1.2  : Dietary copper intake (mg/day)
f_abs_Cu   : 0.60 : Fractional GI absorption of copper (normal)

// ---- ATP7B Function ----
ATP7B_func : 0.05 : Residual ATP7B function (0=null; 1=WT; WD≈0.05)
k_bil_WT   : 3.0  : Biliary Cu excretion rate constant (/day, WT)
k_Cp_synth : 0.15 : Ceruloplasmin synthesis rate (mg/dL/day per Cu unit)
k_Cp_deg   : 0.03 : Ceruloplasmin degradation rate (/day)
Cp_baseline: 30.0 : Baseline ceruloplasmin (mg/dL, WT = 30)

// ---- Copper Redistribution ----
k_hep_NCBC : 0.08 : Rate of hepatic Cu overflow to NCBC (/day)
k_NCBC_uri : 0.15 : NCBC renal clearance to urine (L/day, CL/V)
k_NCBC_brain: 0.005 : NCBC -> brain transfer (/day)
k_NCBC_kid : 0.003 : NCBC -> kidney transfer (/day)
k_NCBC_corn: 0.0005: NCBC -> cornea transfer (/day)
k_brain_out: 0.001 : Brain Cu slow clearance (/day)
k_kidney_out: 0.002: Kidney Cu clearance (/day)
k_cornea_out: 0.0002: Corneal Cu clearance (/day, very slow)

// ---- Liver Pathology ----
k_ROS_gen  : 0.4  : ROS generation rate from hepatic Cu
k_ROS_scav : 0.2  : ROS scavenging capacity (/day, SOD1/catalase)
k_ALT_prod : 5.0  : ALT release rate from hepatocyte damage (U/L/day)
k_ALT_elim : 0.1  : ALT elimination from serum (/day, t½≈5d)
ALT_base   : 20.0 : Baseline serum ALT (U/L)
k_fib_prog : 0.002: Fibrosis progression rate from ROS (/day)
k_fib_reg  : 0.001: Fibrosis regression rate (treatment) (/day)
Fib_max    : 4.0  : Maximum Metavir fibrosis score

// ---- Neurological ----
k_neuro_prog: 0.003: Neurodegeneration progression from brain Cu (/day)
k_neuro_reg : 0.001: Neurological recovery rate (/day)
UWDRS_base : 0.0  : Baseline UWDRS (0 = pre-symptomatic)

// ---- Metallothionein (Hepatic) ----
MT_max     : 50.0 : Maximum MT binding capacity (μg Cu)
MT_Km      : 25.0 : MT Cu binding half-saturation (μg Cu)
k_MT_degrad: 0.05 : MT degradation rate (/day)

// ---- D-Penicillamine (DPA) PK ----
ka_DPA     : 2.0  : DPA absorption rate constant (/hr)
F_DPA      : 0.55 : DPA oral bioavailability
Vd_DPA     : 0.5  : DPA volume of distribution (L/kg × 70 kg = 35 L → /L)
CL_DPA     : 3.5  : DPA total clearance (L/hr)
Kchel_DPA  : 0.3  : DPA chelation potency constant (/hr per μM)
DPA_dose_mg: 0.0  : DPA dose per occasion (mg, 0=off)

// ---- Zinc PK ----
ka_Zn      : 1.5  : Zinc absorption rate constant (/hr)
F_Zn       : 0.15 : Zinc oral bioavailability (variable 10-40%)
Vd_Zn      : 0.25 : Zinc Vd (L/kg)
CL_Zn      : 1.0  : Zinc clearance (L/hr)
IC50_Zn    : 2.0  : Zinc plasma IC50 for Cu absorption inhibition (μg/mL)
Zn_dose_mg : 0.0  : Zinc dose per occasion (mg, 0=off)

// ---- Trientine PK ----
ka_TRI     : 1.8  : Trientine absorption rate constant (/hr)
F_TRI      : 0.45 : Trientine oral bioavailability
Vd_TRI     : 0.4  : Trientine Vd (L/kg)
CL_TRI     : 4.0  : Trientine clearance (L/hr)
Kchel_TRI  : 0.15 : Trientine chelation potency (/hr per μg/mL)
TRI_dose_mg: 0.0  : Trientine dose per occasion (mg, 0=off)

// ---- ALXN1840 (bis-TTM) PK ----
ka_TTM     : 0.8  : ALXN1840 absorption rate constant (/hr)
F_TTM      : 0.30 : ALXN1840 oral bioavailability
Vd_TTM     : 1.0  : ALXN1840 Vd (L/kg)
CL_TTM     : 0.7  : ALXN1840 clearance (L/hr, t½≈19h)
Emax_TTM   : 0.98 : Maximum NCBC reduction by TTM (98%)
EC50_TTM   : 0.5  : TTM EC50 for NCBC reduction (μg/mL)
TTM_dose_mg: 0.0  : ALXN1840 dose per occasion (mg, 0=off)

$INIT @annotated
// Drug PK compartments
GUT_DPA  : 0   : DPA gut depot (mg)
CENT_DPA : 0   : DPA central compartment (mg)
GUT_ZN   : 0   : Zinc gut depot (mg)
CENT_ZN  : 0   : Zinc central (mg)
GUT_TRI  : 0   : Trientine gut depot (mg)
CENT_TRI : 0   : Trientine central (mg)
GUT_TTM  : 0   : ALXN1840 gut depot (mg)
CENT_TTM : 0   : ALXN1840 central (mg)

// Copper kinetics
CU_GI    : 0.5 : GI lumen copper (mg)
CU_HEP   : 80.0 : Hepatic Cu pool (μg/g dw equivalent, WD baseline elevated)
MT_HEP   : 20.0 : Metallothionein-bound hepatic Cu (μg)
CU_NCBC  : 25.0 : Non-ceruloplasmin bound Cu in serum (μg/dL × L)
CP_SERUM : 15.0 : Serum ceruloplasmin (mg/dL; WD baseline low)
CU_URINE : 0.0 : Urinary copper excretion (cumulative, μg/day)
CU_BRAIN : 5.0 : Brain Cu pool (μg; elevated in neuro-WD)
CU_KIDNEY: 2.0 : Kidney Cu pool (μg)
CU_CORNEA: 1.0 : Corneal Cu pool (arbitrary units)

// Liver pathology
ROS_HEP  : 5.0 : Hepatic ROS index (0=normal, 100=max)
ALT_SERUM: 45.0 : Serum ALT (U/L; WD baseline elevated)
FIBROSIS : 0.5 : Hepatic fibrosis score (Metavir F0-F4)

// Neurological
NEURODEGENERATION : 2.0 : Neurodegeneration index (0=none)

$ODE
// ============================================================
// DRUG PK — D-Penicillamine
// ============================================================
double Vd_DPA_tot = Vd_DPA * 70;            // Total Vd (L)
double Cc_DPA = CENT_DPA / Vd_DPA_tot;      // Plasma conc (mg/L ≈ μg/mL)

dxdt_GUT_DPA  = -ka_DPA * GUT_DPA;
dxdt_CENT_DPA = ka_DPA * GUT_DPA - (CL_DPA / Vd_DPA_tot) * CENT_DPA;

// ============================================================
// DRUG PK — Zinc
// ============================================================
double Vd_ZN_tot = Vd_Zn * 70;
double Cc_ZN = CENT_ZN / Vd_ZN_tot;         // plasma Zn (mg/L)

dxdt_GUT_ZN  = -ka_Zn * GUT_ZN;
dxdt_CENT_ZN = ka_Zn * GUT_ZN - (CL_Zn / Vd_ZN_tot) * CENT_ZN;

// Zinc effect: reduces Cu absorption via metallothionein induction
double Zn_eff = Cc_ZN / (IC50_Zn + Cc_ZN);  // 0..1

// ============================================================
// DRUG PK — Trientine
// ============================================================
double Vd_TRI_tot = Vd_TRI * 70;
double Cc_TRI = CENT_TRI / Vd_TRI_tot;

dxdt_GUT_TRI  = -ka_TRI * GUT_TRI;
dxdt_CENT_TRI = ka_TRI * GUT_TRI - (CL_TRI / Vd_TRI_tot) * CENT_TRI;

// ============================================================
// DRUG PK — ALXN1840 (TTM)
// ============================================================
double Vd_TTM_tot = Vd_TTM * 70;
double Cc_TTM = CENT_TTM / Vd_TTM_tot;      // plasma conc (mg/L)

dxdt_GUT_TTM  = -ka_TTM * GUT_TTM;
dxdt_CENT_TTM = ka_TTM * GUT_TTM - (CL_TTM / Vd_TTM_tot) * CENT_TTM;

// TTM effect on NCBC: Emax model
double TTM_eff = Emax_TTM * Cc_TTM / (EC50_TTM + Cc_TTM);  // 0..Emax

// ============================================================
// COPPER GI ABSORPTION
// ============================================================
// Daily Cu intake as zero-order input to CU_GI
// Zinc reduces absorption; DPA/Trientine do NOT block absorption
double Cu_absorb_rate = Cu_intake / 24.0;    // mg/hr continuous
double f_abs_eff = f_abs_Cu * (1.0 - Zn_eff * 0.8); // Zinc reduces absorption 80%

dxdt_CU_GI = Cu_absorb_rate - f_abs_eff * CU_GI;

// Absorbed Cu into hepatic pool
double Cu_in_hep = f_abs_eff * CU_GI;       // mg/hr

// ============================================================
// HEPATIC COPPER — ATP7B Dysfunction Core
// ============================================================
// Biliary excretion: reduced by ATP7B dysfunction
double k_bil_eff = k_bil_WT * ATP7B_func;    // /day effective biliary clearance

// MT saturation (sigmoidal buffer)
double MT_saturation = CU_HEP / (MT_Km + CU_HEP);  // 0..1
double MT_binding = k_MT_degrad * MT_HEP;           // MT turnover
double Cu_to_MT = (1.0 - MT_saturation) * CU_HEP * 0.1;

// DPA chelation effect (mobilizes hepatic Cu)
double DPA_chel = Kchel_DPA * Cc_DPA;
double TRI_chel = Kchel_TRI * Cc_TRI;
double total_chelation = (DPA_chel + TRI_chel) * CU_HEP / (CU_HEP + 10.0);

// Hepatic Cu ODE: input - biliary excretion - MT buffer - chelation
dxdt_CU_HEP = Cu_in_hep * 24.0                 // daily flux (μg → scale)
             - k_bil_eff * CU_HEP               // biliary excretion
             - Cu_to_MT                          // MT binding
             + MT_binding * 0.5                  // MT Cu release
             - total_chelation;                  // drug chelation

dxdt_MT_HEP = Cu_to_MT - MT_binding - MT_HEP * 0.01;  // MT dynamics

// ============================================================
// CERULOPLASMIN SYNTHESIS
// ============================================================
// Cp synthesis requires Cu via ATP7B; low in WD
double Cp_synth = k_Cp_synth * CU_HEP * ATP7B_func;
dxdt_CP_SERUM = Cp_synth - k_Cp_deg * CP_SERUM;

// ============================================================
// NON-CERULOPLASMIN BOUND COPPER (NCBC)
// ============================================================
// NCBC = overflow from saturated hepatic MT + reduced biliary excretion
double k_NCBC_gen = k_hep_NCBC * MT_saturation * CU_HEP;

// TTM and DPA/Trientine reduce NCBC
double NCBC_chel_out = (DPA_chel * 0.6 + TRI_chel * 0.4) * CU_NCBC / (CU_NCBC + 5.0);
double NCBC_TTM_out  = TTM_eff * CU_NCBC;  // TTM very effective

dxdt_CU_NCBC = k_NCBC_gen
              - k_NCBC_uri   * CU_NCBC
              - k_NCBC_brain * CU_NCBC
              - k_NCBC_kid   * CU_NCBC
              - k_NCBC_corn  * CU_NCBC
              - NCBC_chel_out
              - NCBC_TTM_out;

// ============================================================
// URINARY COPPER (cumulative daily marker)
// ============================================================
double Cu_urine_rate = k_NCBC_uri * CU_NCBC                // baseline filtration
                     + DPA_chel * CU_HEP * 0.01            // DPA-chelated Cu
                     + TRI_chel * CU_HEP * 0.005;          // TRI-chelated Cu
dxdt_CU_URINE = Cu_urine_rate;   // tracks rate; reset daily in table

// ============================================================
// BRAIN COPPER
// ============================================================
dxdt_CU_BRAIN = k_NCBC_brain * CU_NCBC
               - k_brain_out  * CU_BRAIN
               - TTM_eff * 0.01 * CU_BRAIN;  // TTM slight CNS effect

// ============================================================
// KIDNEY COPPER
// ============================================================
dxdt_CU_KIDNEY = k_NCBC_kid * CU_NCBC
                - k_kidney_out * CU_KIDNEY;

// ============================================================
// CORNEAL COPPER (Kayser-Fleischer proxy)
// ============================================================
dxdt_CU_CORNEA = k_NCBC_corn * CU_NCBC
                - k_cornea_out * CU_CORNEA;

// ============================================================
// HEPATIC REACTIVE OXYGEN SPECIES (ROS)
// ============================================================
// Fenton-like: Cu1+ reacts with H2O2 → OH•
double ROS_gen  = k_ROS_gen * (CU_HEP / 100.0) * (CU_HEP / 100.0);  // quadratic
double ROS_scav = k_ROS_scav * ROS_HEP;                               // linear scavenging
dxdt_ROS_HEP = ROS_gen - ROS_scav;
double ROS_norm = ROS_HEP;  // alias for readability

// ============================================================
// SERUM ALT (Liver Damage Marker)
// ============================================================
// ALT generated from hepatocyte damage (ROS-driven)
double ALT_gen = k_ALT_prod * (ROS_norm / 20.0);  // normalized to ROS index
dxdt_ALT_SERUM = ALT_gen - k_ALT_elim * (ALT_SERUM - ALT_base);

// ============================================================
// HEPATIC FIBROSIS
// ============================================================
// Driven by cumulative ROS; slow process
double fib_prog  = k_fib_prog * ROS_norm * FIBROSIS * (1.0 - FIBROSIS / Fib_max);
double fib_reg   = k_fib_reg  * (DPA_chel + TRI_chel + TTM_eff) * FIBROSIS;
dxdt_FIBROSIS = fib_prog - fib_reg;
if(FIBROSIS > Fib_max) dxdt_FIBROSIS = 0.0;
if(FIBROSIS < 0.0)     dxdt_FIBROSIS = 0.0;

// ============================================================
// NEURODEGENERATION
// ============================================================
double neuro_prog = k_neuro_prog * CU_BRAIN;
double neuro_reg  = k_neuro_reg * (TTM_eff * 0.5 + DPA_chel * 0.1);
dxdt_NEURODEGENERATION = neuro_prog - neuro_reg;
if(NEURODEGENERATION < 0.0) dxdt_NEURODEGENERATION = 0.0;

$TABLE
// ============================================================
// Derived PK/PD Outputs
// ============================================================
capture Cc_DPA_out   = CENT_DPA / (Vd_DPA * 70);   // DPA plasma (mg/L ≈ μg/mL)
capture Cc_ZN_out    = CENT_ZN  / (Vd_Zn  * 70);   // Zn plasma (mg/L)
capture Cc_TRI_out   = CENT_TRI / (Vd_TRI * 70);   // Trientine plasma (mg/L)
capture Cc_TTM_out   = CENT_TTM / (Vd_TTM * 70);   // ALXN1840 plasma (mg/L)

capture Cu_hep_out   = CU_HEP;               // Hepatic Cu (μg/g dw)
capture NCBC_out     = CU_NCBC;              // NCBC (μg/dL equivalent)
capture Cp_out       = CP_SERUM;             // Ceruloplasmin (mg/dL)
capture Cu_urine_out = CU_URINE;             // Urinary Cu (μg/day rate)
capture Cu_brain_out = CU_BRAIN;             // Brain Cu index
capture Cu_corn_out  = CU_CORNEA;            // Corneal Cu index (KF proxy)

capture ROS_out      = ROS_HEP;             // Hepatic ROS index
capture ALT_out      = ALT_SERUM;            // Serum ALT (U/L)
capture Fib_out      = FIBROSIS;             // Fibrosis score (F0-F4)
capture Neuro_out    = NEURODEGENERATION;    // Neurodegeneration index

// Clinical derived indices
capture Cu_serum_total = Cp_out * 3.15 + NCBC_out;  // Total serum Cu (approx μg/dL)
capture KF_index       = Cu_corn_out;                // KF ring severity proxy
capture UWDRS_approx   = NEURODEGENERATION * 5.0;    // Approximate UWDRS

'

## ============================================================
## Build & Compile Model
## ============================================================
wd_mod <- mcode("WilsonDisease_QSP", wd_code)

## ============================================================
## Helper: Build dosing regimen
## ============================================================
build_dose <- function(drug = "DPA",
                       dose_mg = 500,
                       interval_hr = 8,
                       n_doses = 365 * 3,
                       cmt_gut = "GUT_DPA",
                       bioavail = 0.55) {
  doses <- ev(
    amt    = dose_mg * bioavail,
    cmt    = cmt_gut,
    ii     = interval_hr,
    addl   = n_doses - 1,
    time   = 0
  )
  return(doses)
}

## ============================================================
## Simulation Scenarios
## ============================================================

# Common settings
sim_years <- 5
sim_hrs   <- sim_years * 365 * 24
tsamp     <- seq(0, sim_hrs, by = 24)   # daily samples

## ---- Scenario 1: Untreated WD (Hepatic) ----
params_WD_base <- list(
  ATP7B_func = 0.05,
  Cu_intake  = 1.2,
  DPA_dose_mg = 0, Zn_dose_mg = 0, TRI_dose_mg = 0, TTM_dose_mg = 0
)

dose_none <- ev(time = 0, amt = 0, cmt = "GUT_DPA")

sim_S1 <- wd_mod %>%
  param(params_WD_base) %>%
  mrgsim(ev = dose_none, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S1: Untreated WD", Day = time / 24)

## ---- Scenario 2: D-Penicillamine 500 mg TID ----
dose_DPA <- ev(
  amt  = 500 * 0.55,   # bioavailable amount
  cmt  = "GUT_DPA",
  ii   = 8,
  addl = sim_years * 365 * 3 - 1,
  time = 0
)

params_DPA <- modifyList(params_WD_base, list(
  Kchel_DPA = 0.3,
  DPA_dose_mg = 500
))

sim_S2 <- wd_mod %>%
  param(params_DPA) %>%
  mrgsim(ev = dose_DPA, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S2: DPA 500mg TID", Day = time / 24)

## ---- Scenario 3: Zinc Acetate 50 mg TID (Maintenance) ----
dose_ZN <- ev(
  amt  = 50 * 0.15,
  cmt  = "GUT_ZN",
  ii   = 8,
  addl = sim_years * 365 * 3 - 1,
  time = 0
)

params_ZN <- modifyList(params_WD_base, list(
  IC50_Zn    = 2.0,
  Zn_dose_mg = 50
))

sim_S3 <- wd_mod %>%
  param(params_ZN) %>%
  mrgsim(ev = dose_ZN, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S3: Zinc 50mg TID", Day = time / 24)

## ---- Scenario 4: Trientine 500 mg TID ----
dose_TRI <- ev(
  amt  = 500 * 0.45,
  cmt  = "GUT_TRI",
  ii   = 8,
  addl = sim_years * 365 * 3 - 1,
  time = 0
)

params_TRI <- modifyList(params_WD_base, list(
  Kchel_TRI  = 0.15,
  TRI_dose_mg = 500
))

sim_S4 <- wd_mod %>%
  param(params_TRI) %>%
  mrgsim(ev = dose_TRI, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S4: Trientine 500mg TID", Day = time / 24)

## ---- Scenario 5: ALXN1840 15 mg QD ----
dose_TTM <- ev(
  amt  = 15 * 0.30,
  cmt  = "GUT_TTM",
  ii   = 24,
  addl = sim_years * 365 - 1,
  time = 0
)

params_TTM <- modifyList(params_WD_base, list(
  Emax_TTM   = 0.98,
  EC50_TTM   = 0.5,
  TTM_dose_mg = 15
))

sim_S5 <- wd_mod %>%
  param(params_TTM) %>%
  mrgsim(ev = dose_TTM, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S5: ALXN1840 15mg QD", Day = time / 24)

## ---- Scenario 6: DPA -> Zinc Maintenance Switch (1 yr DPA then switch) ----
dose_DPA_lead <- ev(
  amt  = 500 * 0.55, cmt = "GUT_DPA",
  ii   = 8, addl = 365 * 3 - 1, time = 0
)
dose_ZN_maint <- ev(
  amt  = 50 * 0.15, cmt = "GUT_ZN",
  ii   = 8, addl = 4 * 365 * 3 - 1, time = 365 * 24  # start at 1 year
)
dose_S6 <- c(dose_DPA_lead, dose_ZN_maint)

sim_S6 <- wd_mod %>%
  param(modifyList(params_DPA, list(IC50_Zn = 2.0, Zn_dose_mg = 50))) %>%
  mrgsim(ev = dose_S6, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S6: DPA→Zinc Switch (yr 1)", Day = time / 24)

## ---- Scenario 7: ALXN1840 + Trientine Combination ----
dose_combo <- c(
  ev(amt = 15 * 0.30,  cmt = "GUT_TTM", ii = 24, addl = sim_years * 365 - 1, time = 0),
  ev(amt = 250 * 0.45, cmt = "GUT_TRI", ii = 8,  addl = sim_years * 365 * 3 - 1, time = 0)
)

params_combo <- modifyList(params_TTM, list(
  Kchel_TRI = 0.15, TRI_dose_mg = 250
))

sim_S7 <- wd_mod %>%
  param(params_combo) %>%
  mrgsim(ev = dose_combo, end = sim_hrs, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "S7: ALXN1840+Trientine", Day = time / 24)

## ---- Scenario 8: Healthy (WT ATP7B Control) ----
params_WT <- list(
  ATP7B_func = 1.0,
  Cu_intake  = 1.2
)

sim_S8 <- wd_mod %>%
  param(params_WT) %>%
  mrgsim(ev = dose_none, end = sim_hrs, delta = 24,
         init = list(CU_HEP=15, CP_SERUM=30, CU_NCBC=5, ALT_SERUM=20,
                     FIBROSIS=0, NEURODEGENERATION=0, CU_BRAIN=1, ROS_HEP=1)) %>%
  as.data.frame() %>%
  mutate(Scenario = "S8: Healthy WT Control", Day = time / 24)

## ============================================================
## Combine Results
## ============================================================
all_sims <- bind_rows(sim_S1, sim_S2, sim_S3, sim_S4,
                      sim_S5, sim_S6, sim_S7, sim_S8)

## ============================================================
## Plotting Functions
## ============================================================

scenario_colors <- c(
  "S1: Untreated WD"          = "#e74c3c",
  "S2: DPA 500mg TID"         = "#e67e22",
  "S3: Zinc 50mg TID"         = "#2ecc71",
  "S4: Trientine 500mg TID"   = "#3498db",
  "S5: ALXN1840 15mg QD"      = "#9b59b6",
  "S6: DPA→Zinc Switch (yr 1)"= "#1abc9c",
  "S7: ALXN1840+Trientine"    = "#f39c12",
  "S8: Healthy WT Control"    = "#95a5a6"
)

plot_copper_kinetics <- function(df) {
  df %>%
    pivot_longer(c(NCBC_out, Cp_out, Cu_serum_total),
                 names_to = "Biomarker", values_to = "Value") %>%
    mutate(Biomarker = recode(Biomarker,
      NCBC_out = "NCBC (free Cu, μg/dL)",
      Cp_out = "Ceruloplasmin (mg/dL)",
      Cu_serum_total = "Total Serum Cu (μg/dL)"
    )) %>%
    ggplot(aes(x = Day / 365, y = Value, color = Scenario)) +
    geom_line(size = 0.8) +
    facet_wrap(~Biomarker, scales = "free_y", ncol = 1) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Copper Kinetics by Treatment Scenario",
         x = "Years", y = "Concentration") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          strip.background = element_rect(fill = "#2c3e50"),
          strip.text = element_text(color = "white"))
}

plot_liver_outcomes <- function(df) {
  df %>%
    pivot_longer(c(ALT_out, Fib_out, ROS_out),
                 names_to = "Marker", values_to = "Value") %>%
    mutate(Marker = recode(Marker,
      ALT_out = "Serum ALT (U/L)",
      Fib_out = "Fibrosis Score (Metavir)",
      ROS_out = "Hepatic ROS Index"
    )) %>%
    ggplot(aes(x = Day / 365, y = Value, color = Scenario)) +
    geom_line(size = 0.8) +
    facet_wrap(~Marker, scales = "free_y", ncol = 1) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Hepatic Outcomes by Treatment Scenario",
         x = "Years", y = "Value") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          strip.background = element_rect(fill = "#1a5276"),
          strip.text = element_text(color = "white"))
}

plot_neuro_outcomes <- function(df) {
  df %>%
    pivot_longer(c(Cu_brain_out, Neuro_out, UWDRS_approx, KF_index),
                 names_to = "Marker", values_to = "Value") %>%
    mutate(Marker = recode(Marker,
      Cu_brain_out = "Brain Cu Pool (μg)",
      Neuro_out    = "Neurodegeneration Index",
      UWDRS_approx = "Approx. UWDRS Score",
      KF_index     = "KF Ring Index (Corneal Cu)"
    )) %>%
    ggplot(aes(x = Day / 365, y = Value, color = Scenario)) +
    geom_line(size = 0.8) +
    facet_wrap(~Marker, scales = "free_y", ncol = 2) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Neurological & Corneal Outcomes",
         x = "Years", y = "Value") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          strip.background = element_rect(fill = "#154360"),
          strip.text = element_text(color = "white"))
}

plot_drug_pk <- function(df, scenario_filter = c("S2: DPA 500mg TID",
                                                   "S3: Zinc 50mg TID",
                                                   "S4: Trientine 500mg TID",
                                                   "S5: ALXN1840 15mg QD")) {
  df %>%
    filter(Scenario %in% scenario_filter, Day <= 30) %>%
    pivot_longer(c(Cc_DPA_out, Cc_ZN_out, Cc_TRI_out, Cc_TTM_out),
                 names_to = "Drug", values_to = "Conc") %>%
    mutate(Drug = recode(Drug,
      Cc_DPA_out = "D-Penicillamine",
      Cc_ZN_out  = "Zinc",
      Cc_TRI_out = "Trientine",
      Cc_TTM_out = "ALXN1840"
    )) %>%
    filter(Conc > 0.001) %>%
    ggplot(aes(x = Day, y = Conc, color = Scenario)) +
    geom_line(size = 0.8) +
    facet_wrap(~Drug, scales = "free_y", ncol = 2) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Drug Plasma Concentrations (First 30 Days)",
         x = "Day", y = "Plasma Concentration (mg/L)") +
    theme_bw(base_size = 12)
}

## ============================================================
## Run and Display All Plots (when run interactively)
## ============================================================
if (interactive()) {
  print(plot_copper_kinetics(all_sims))
  print(plot_liver_outcomes(all_sims))
  print(plot_neuro_outcomes(all_sims))
  print(plot_drug_pk(all_sims))

  # Summary table at 1-year, 3-year, 5-year
  cat("\n=== KEY OUTCOMES SUMMARY ===\n")
  summary_tbl <- all_sims %>%
    filter(Day %in% c(365, 1095, 1825)) %>%
    group_by(Scenario, Day) %>%
    summarise(
      NCBC_mean      = mean(NCBC_out, na.rm=TRUE),
      ALT_mean       = mean(ALT_out, na.rm=TRUE),
      Fibrosis_mean  = mean(Fib_out, na.rm=TRUE),
      UWDRS_mean     = mean(UWDRS_approx, na.rm=TRUE),
      KF_mean        = mean(KF_index, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(Year = Day / 365)
  print(summary_tbl)
}

cat("Wilson Disease mrgsolve model loaded. Run plot functions to visualize.\n")
cat("Scenarios:\n")
cat("  S1: Untreated WD (hepatic presentation)\n")
cat("  S2: D-Penicillamine 500 mg TID\n")
cat("  S3: Zinc Acetate 50 mg TID\n")
cat("  S4: Trientine 500 mg TID\n")
cat("  S5: ALXN1840 (bis-TTM) 15 mg QD [ATLAS regimen]\n")
cat("  S6: DPA induction → Zinc maintenance switch at 1 year\n")
cat("  S7: ALXN1840 + Trientine combination\n")
cat("  S8: Healthy WT control (reference)\n")
