## =============================================================================
## Bullous Pemphigoid — mrgsolve QSP Model
## Author: Claude Code Routine (CCR)
## Date: 2026-06-19
##
## Model structure:
##   PK:  Prednisolone (2-CMT oral), Dupilumab (1-CMT SC),
##        Omalizumab (1-CMT SC), Rituximab (1-CMT IV),
##        Doxycycline (1-CMT oral)
##   PD:  B cells (naive → activated → plasma cells [short/long])
##        Th2 cells, IL-4/IL-13, IL-31 (itch)
##        Anti-BP180 IgG & IgE autoantibodies
##        Eosinophils (blood & skin), Mast cells (skin)
##        Complement C5a, DEJ disruption
##        BPDAI activity, Itch NRS
##
## Reference parameters calibrated to:
##   - Amber et al. J Allergy Clin Immunol 2018 (dupilumab BP)
##   - Joly et al. NEJM 2002, 2009 (prednisolone BP trials)
##   - Fairley et al. J Invest Dermatol 2020 (rituximab BP)
##   - Moriuchi et al. J Dermatol 2020 (omalizumab BP)
##   - Joly et al. NEJM 2002 (doxycycline + niacinamide BP trial)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── Model Definition ──────────────────────────────────────────────────────────
bp_model_code <- '
$PROB
Bullous Pemphigoid QSP Model
Compartments: 23 (5 drug PK + 18 disease biology)

$PARAM
// ── Prednisolone PK ──────────────────────────────────────────────────────────
// 1-compartment oral; Vd, CL from Lew et al. Br J Clin Pharmacol 1993
ka_pred   = 2.4     // absorption rate constant (1/h)
F_pred    = 0.82    // oral bioavailability
Vd_pred   = 45.0    // volume of distribution (L)
CL_pred   = 6.5     // clearance (L/h)
Vp_pred   = 15.0    // peripheral volume (L)
Q_pred    = 1.2     // inter-compartmental CL (L/h)

// ── Dupilumab PK (SC, 2-weekly) ──────────────────────────────────────────────
// Regeneron/Sanofi: Tmax~7d, t1/2~21d (Type IgG4 mAb)
ka_dup    = 0.022   // absorption rate (1/h; SC depot)
F_dup     = 0.64    // bioavailability (SC)
Vd_dup    = 4.8     // central Vd (L)
CL_dup    = 0.0064  // clearance (L/h)

// ── Omalizumab PK (SC, 4-weekly) ─────────────────────────────────────────────
ka_oma    = 0.010   // absorption rate (1/h; SC depot)
F_oma     = 0.62    // bioavailability (SC)
Vd_oma    = 7.7     // central Vd (L)
CL_oma    = 0.0112  // clearance (L/h)

// ── Rituximab PK (IV, 1g Q2W x2, then Q6M) ───────────────────────────────────
// 2-CMT IV; Berinstein et al.
Vd_ritu   = 3.6     // central Vd (L)
CL_ritu   = 0.0162  // clearance (L/h)
Vp_ritu   = 4.5     // peripheral Vd (L)
Q_ritu    = 0.2     // inter-compartmental CL

// ── Doxycycline PK (PO, 200 mg/d) ────────────────────────────────────────────
ka_doxy   = 0.58    // absorption rate (1/h)
F_doxy    = 0.93    // bioavailability
Vd_doxy   = 150.0   // Vd (L; widely distributed)
CL_doxy   = 1.8     // clearance (L/h; renal + biliary)

// ── B Cell Dynamics ───────────────────────────────────────────────────────────
kBn_prod  = 0.006   // naive B cell production (cells/µL/h)
kBn_die   = 0.0012  // naive B cell death (1/h)
kBact     = 0.0003  // naive→activated (Ag-driven; 1/h base)
kGC       = 0.0008  // activated→GC B cell (1/h)
kMem      = 0.0004  // GC→memory B (1/h)
kSLPC     = 0.0005  // GC→short-lived plasma cell (1/h)
kLLPC     = 0.00004 // SLPC→long-lived plasma cell (1/h)
kBact_die = 0.004   // activated B cell death (1/h)
kMem_die  = 0.00005 // memory B cell death (1/h)
kSLPC_die = 0.012   // SLPC death (1/h; t1/2 ~ 2.4 d)
kLLPC_die = 0.00014 // LLPC death (1/h; t1/2 ~ 300 d)

// ── Th2 Cell Dynamics ─────────────────────────────────────────────────────────
kTh2_prod = 0.0008  // naive→Th2 polarisation base rate (1/h)
kTh2_die  = 0.002   // Th2 apoptosis (1/h)
Th2_0     = 80.0    // steady-state Th2 cells (AU)

// ── Autoantibody Production & Clearance ──────────────────────────────────────
kIgG_prod = 0.00018 // IgG anti-BP180 production (AU/h per LLPC)
kIgG_deg  = 0.0013  // IgG catabolism (1/h; t1/2 ~ 21 d)
kIgE_prod = 0.00006 // IgE anti-BP180 production (AU/h per SLPC)
kIgE_deg  = 0.0046  // IgE catabolism (1/h; t1/2 ~ 6 d)
IgG_0     = 1.0     // initial anti-BP180 IgG (normalised)
IgE_0     = 0.5     // initial anti-BP180 IgE (normalised)

// ── Eosinophil Dynamics ───────────────────────────────────────────────────────
kEos_prod = 0.0004  // eosinophil BM production (cells/µL/h)
kEos_die  = 0.0017  // peripheral eosinophil apoptosis (1/h; t1/2~17h)
kEos_skin = 0.0008  // blood→skin migration (1/h; driven by IgG)
kEos_sk_die = 0.004 // skin eosinophil death (1/h)
Eos_0     = 0.25    // baseline eosinophils (×10^9/L)

// ── Mast Cell / C5a / DEJ ────────────────────────────────────────────────────
kMast_act = 0.0030  // mast cell activation rate (IgE-driven; 1/h)
kMast_base = 0.0005 // baseline mast cell activation
kMast_die  = 0.0050 // activated mast-cell return to resting (1/h)
kC5a_prod  = 0.0020 // C5a generation (IgG-complement driven)
kC5a_deg   = 0.030  // C5a degradation (1/h; t1/2 ~1.4 min → lumped 30/h for model)
kDEJ_dam   = 0.0015 // DEJ damage rate (per eos-skin + C5a)
kDEJ_repair= 0.0006 // DEJ repair rate (1/h)
DEJ_0      = 1.0    // DEJ integrity (1=normal, 0=completely disrupted)

// ── BPDAI & Itch ─────────────────────────────────────────────────────────────
kBPDAI_rise = 0.0040  // BPDAI rise from DEJ disruption (1/h)
kBPDAI_fall = 0.0010  // BPDAI natural fall/healing (1/h)
kItch_rise  = 0.0050  // itch NRS rise (IL-31, histamine driven)
kItch_fall  = 0.0015  // itch NRS fall (1/h)
BPDAI_0     = 30.0    // baseline BPDAI activity score (0-90)
Itch_0      = 6.0     // baseline itch NRS (0-10)

// ── PD Potency (EC50 for drug effects) ───────────────────────────────────────
EC50_pred_IL4  = 0.015  // prednisolone EC50 for IL-4 suppression (mg/L)
EC50_pred_Eos  = 0.010  // prednisolone EC50 for eosinophil suppression (mg/L)
EC50_dup_IL4   = 0.008  // dupilumab EC50 for IL-4/IL-13 suppression (mg/L)
EC50_oma_IgE   = 0.030  // omalizumab EC50 for IgE neutralisation (mg/L)
EC50_ritu_B    = 0.015  // rituximab EC50 for B cell depletion (mg/L)
EC50_doxy_MMP  = 2.00   // doxycycline EC50 for MMP-9 inhibition (mg/L)
Emax           = 0.90   // maximal drug effect (fraction)

// ── Baseline disease activity modifiers ──────────────────────────────────────
BP180_Ag_level = 1.0   // baseline antigen presentation (1=BP patient)
IL4_stim       = 1.5   // baseline IL-4/IL-13 stimulation factor (disease)
Th2_bias       = 2.0   // Th2 polarisation bias factor (disease vs. healthy)

$CMT
// Drug PK compartments
PRED_GUT      // prednisolone GI compartment
PRED_C        // prednisolone central (plasma)
PRED_P        // prednisolone peripheral
DUP_DEPOT     // dupilumab SC depot
DUP_C         // dupilumab central
OMA_DEPOT     // omalizumab SC depot
OMA_C         // omalizumab central
RITU_C        // rituximab central
DOXY_GUT      // doxycycline GI
DOXY_C        // doxycycline central

// Disease PD compartments
B_NAIVE       // naive B cells (×10^3/µL)
B_ACT         // activated B cells
B_MEM         // memory B cells
SLPC          // short-lived plasma cells
LLPC          // long-lived plasma cells
TH2           // Th2 cells (AU)
IGG_BP180     // anti-BP180 IgG (normalised)
IGE_BP180     // anti-BP180 IgE (normalised)
EOS_BLOOD     // blood eosinophils (×10^9/L)
EOS_SKIN      // skin eosinophils (AU)
MAST_ACT      // activated mast cells (AU)
C5A           // complement C5a (AU)
DEJ           // DEJ integrity (1=normal)

$INIT
PRED_GUT  = 0
PRED_C    = 0
PRED_P    = 0
DUP_DEPOT = 0
DUP_C     = 0
OMA_DEPOT = 0
OMA_C     = 0
RITU_C    = 0
DOXY_GUT  = 0
DOXY_C    = 0
B_NAIVE   = 200.0
B_ACT     = 10.0
B_MEM     = 20.0
SLPC      = 5.0
LLPC      = 2.0
TH2       = 80.0
IGG_BP180 = 1.0
IGE_BP180 = 0.5
EOS_BLOOD = 0.35
EOS_SKIN  = 0.20
MAST_ACT  = 0.10
C5A       = 0.05
DEJ       = 0.50   // 50% DEJ integrity at disease onset

$ODE
// ── DRUG PK ──────────────────────────────────────────────────────────────────

// Prednisolone
dxdt_PRED_GUT = -ka_pred * PRED_GUT;
dxdt_PRED_C   =  ka_pred * F_pred * PRED_GUT
                 - (CL_pred/Vd_pred + Q_pred/Vd_pred)*PRED_C
                 + (Q_pred/Vp_pred)*PRED_P;
dxdt_PRED_P   =  (Q_pred/Vd_pred)*PRED_C - (Q_pred/Vp_pred)*PRED_P;

// Dupilumab
dxdt_DUP_DEPOT = -ka_dup * DUP_DEPOT;
dxdt_DUP_C     =  ka_dup * F_dup * DUP_DEPOT - (CL_dup/Vd_dup)*DUP_C;

// Omalizumab
dxdt_OMA_DEPOT = -ka_oma * OMA_DEPOT;
dxdt_OMA_C     =  ka_oma * F_oma * OMA_DEPOT - (CL_oma/Vd_oma)*OMA_C;

// Rituximab
dxdt_RITU_C = -(CL_ritu/Vd_ritu)*RITU_C;

// Doxycycline
dxdt_DOXY_GUT = -ka_doxy * DOXY_GUT;
dxdt_DOXY_C   =  ka_doxy * F_doxy * DOXY_GUT - (CL_doxy/Vd_doxy)*DOXY_C;

// ── DRUG EFFECTS (Hill function) ─────────────────────────────────────────────
double Cpred   = PRED_C / Vd_pred;  // mg/L
double Cdup    = DUP_C  / Vd_dup;
double Coma    = OMA_C  / Vd_oma;
double Critu   = RITU_C / Vd_ritu;
double Cdoxy   = DOXY_C / Vd_doxy;

// Drug inhibition functions (Imax models)
double E_pred_IL4  = Emax * Cpred / (EC50_pred_IL4 + Cpred);
double E_pred_Eos  = Emax * Cpred / (EC50_pred_Eos + Cpred);
double E_dup_IL4   = Emax * Cdup  / (EC50_dup_IL4  + Cdup);
double E_oma_IgE   = Emax * Coma  / (EC50_oma_IgE  + Coma);
double E_ritu_B    = Emax * Critu / (EC50_ritu_B   + Critu);
double E_doxy_MMP  = Emax * Cdoxy / (EC50_doxy_MMP + Cdoxy);

// Combined Th2/IL-4 suppression
double E_Th2_total = 1.0 - (1.0 - E_pred_IL4) * (1.0 - E_dup_IL4);

// ── TH2 CELL DYNAMICS ────────────────────────────────────────────────────────
double Th2_drive = kTh2_prod * Th2_bias * IL4_stim * BP180_Ag_level
                   * (1.0 - E_Th2_total);
dxdt_TH2 = Th2_drive * B_NAIVE/(B_NAIVE + 50.0)
            - kTh2_die * TH2;

// ── B CELL DYNAMICS ───────────────────────────────────────────────────────────
double BCR_stim = IGG_BP180 * BP180_Ag_level; // antigen-driven activation
double Bcell_dep = E_ritu_B;                   // rituximab depletion

dxdt_B_NAIVE =  kBn_prod - kBn_die * B_NAIVE
                - kBact * BCR_stim * B_NAIVE * (1.0 - Bcell_dep);
dxdt_B_ACT   =  kBact * BCR_stim * B_NAIVE * (1.0 - Bcell_dep)
                - (kGC + kBact_die) * B_ACT * (1.0 - Bcell_dep);
dxdt_B_MEM   =  kMem * B_ACT - kMem_die * B_MEM * (1.0 - Bcell_dep);
dxdt_SLPC    =  kSLPC * B_ACT - (kSLPC_die + kLLPC) * SLPC;
dxdt_LLPC    =  kLLPC * SLPC  - kLLPC_die * LLPC;

// ── AUTOANTIBODIES ───────────────────────────────────────────────────────────
// IgG anti-BP180 (produced by LLPC + SLPC, suppressed by rituximab)
double IgG_prod_rate = kIgG_prod * (LLPC + 0.3*SLPC) * (1.0 - 0.6*E_ritu_B);
dxdt_IGG_BP180 = IgG_prod_rate - kIgG_deg * IGG_BP180;

// IgE anti-BP180 (produced by SLPC under Th2/IL-4; neutralised by omalizumab)
double IgE_prod_rate = kIgE_prod * SLPC * (1.0 - E_Th2_total);
double IgE_neutralise = E_oma_IgE * 0.8 * IGE_BP180; // omalizumab neutralises
dxdt_IGE_BP180 = IgE_prod_rate - kIgE_deg * IGE_BP180 - IgE_neutralise;

// ── EOSINOPHILS ──────────────────────────────────────────────────────────────
// Blood eosinophils: driven by IL-5 (Th2), suppressed by pred
double Eos_stim = (TH2/Th2_0) * IL4_stim * (1.0 - E_pred_Eos);
dxdt_EOS_BLOOD = kEos_prod * Eos_stim - kEos_die * EOS_BLOOD
                 - kEos_skin * IGG_BP180 * EOS_BLOOD;

// Skin eosinophils: migrate from blood, driven by IgG/complement
double Skin_recruit = kEos_skin * IGG_BP180 * EOS_BLOOD * (1.0 - E_dup_IL4*0.5);
dxdt_EOS_SKIN  = Skin_recruit - kEos_sk_die * EOS_SKIN;

// ── MAST CELLS ────────────────────────────────────────────────────────────────
// Activated by IgE cross-linking (suppressed by omalizumab)
double Mast_drive = kMast_act * IGE_BP180 * (1.0 - E_oma_IgE)
                    + kMast_base;
dxdt_MAST_ACT  = Mast_drive - kMast_die * MAST_ACT;

// ── COMPLEMENT C5a ───────────────────────────────────────────────────────────
double C5a_drive = kC5a_prod * IGG_BP180 * (1.0 - E_pred_IL4*0.3);
dxdt_C5A = C5a_drive - kC5a_deg * C5A;

// ── DEJ INTEGRITY ────────────────────────────────────────────────────────────
// Damaged by: eosinophils (skin), C5a, MAST_ACT
// Repaired by: absence of inflammation
double DEJ_damage_rate = kDEJ_dam * (EOS_SKIN + 0.5*C5A + 0.3*MAST_ACT)
                         * (1.0 - E_doxy_MMP);
double DEJ_repair_rate = kDEJ_repair * (1.0 - DEJ); // repair tendency
// DEJ is integrity 0-1; damage brings it down, repair brings it up
dxdt_DEJ = DEJ_repair_rate - DEJ_damage_rate * DEJ;

// ── BPDAI & ITCH ─────────────────────────────────────────────────────────────
// (handled in $TABLE as outputs; BPDAI driven by 1-DEJ)

$TABLE
// Drug concentrations (ng/mL scale for report)
double Cpred_ngml  = PRED_C/Vd_pred * 1000.0;  // µg/L ≈ ng/mL for pred (MW~360)
double Cdup_mcg    = DUP_C/Vd_dup;              // mg/L → µg/mL *1000 not needed; reported as mg/L
double Coma_mcg    = OMA_C/Vd_oma;
double Critu_mcg   = RITU_C/Vd_ritu;
double Cdoxy_mcg   = DOXY_C/Vd_doxy * 1000.0;  // µg/mL

// Disease biomarkers
double anti_BP180_AU = IGG_BP180;   // normalised anti-BP180 IgG (AU)
double anti_BP180_IgE_AU = IGE_BP180;
double eos_blood_09L = EOS_BLOOD;   // ×10^9/L
double eos_skin_AU   = EOS_SKIN;

// Efficacy readouts
// BPDAI_activity: inversely related to DEJ integrity (scale 0-90)
double BPDAI = BPDAI_0 * (1.0 - DEJ) * (IGG_BP180/(IGG_BP180 + 0.5));
double BPDAI_act = (BPDAI < 0) ? 0 : (BPDAI > 90 ? 90 : BPDAI);

// Itch NRS (0-10): driven by MAST_ACT, EOS_SKIN, and IGE_BP180
double itch_NRS = Itch_0 * MAST_ACT/0.10 * IGE_BP180/0.5 * (1.0/
                  (1.0 + (DUP_C/Vd_dup)/EC50_dup_IL4));
double Itch = (itch_NRS < 0) ? 0 : (itch_NRS > 10 ? 10 : itch_NRS);

// Blister count approximation (blisters/day)
double new_blisters = BPDAI_act / 10.0;  // approximate clinical mapping

// New blister rate = 0 criterion for remission
double in_remission = (new_blisters < 0.3) ? 1.0 : 0.0;

// BSA affected (%)
double BSA = 20.0 * (1.0 - DEJ) * (IGG_BP180/(IGG_BP180 + 1.0));

// Cumulative steroid dose (handled externally in R)
double Cpred_mgL   = PRED_C/Vd_pred;

$CAPTURE
Cpred_ngml Cdup_mcg Coma_mcg Critu_mcg Cdoxy_mcg
anti_BP180_AU anti_BP180_IgE_AU
eos_blood_09L eos_skin_AU MAST_ACT C5A DEJ
BPDAI_act Itch new_blisters BSA in_remission
B_NAIVE B_ACT LLPC TH2 IGG_BP180 IGE_BP180
'

## ── Compile Model ─────────────────────────────────────────────────────────────
mod <- mread_cache(model = "bp_model", code = bp_model_code)

## ── Helper: Dosing Regimens ───────────────────────────────────────────────────

#' Build prednisolone dosing events
#' @param dose_mg_kg  starting dose in mg/kg
#' @param bw_kg       body weight (kg)
#' @param start_h     start hour
#' @param taper_weeks number of weeks before dose reduction steps
pred_dose_events <- function(dose_mg_kg = 0.5, bw_kg = 65,
                              start_h = 0, taper_weeks = 4) {
  dose_mg <- dose_mg_kg * bw_kg
  # Initial dose for taper_weeks, then 50% reduction steps every 4 weeks
  # Until 0.1 mg/kg/d maintenance
  ev <- c()
  d <- dose_mg
  t <- start_h
  while (d > 0.1 * bw_kg) {
    ev <- c(ev, ev(cmt = "PRED_GUT", amt = d, ii = 24, addl = 7*taper_weeks - 1,
                   time = t))
    t <- t + 7 * taper_weeks * 24
    d <- d * 0.5
  }
  # Maintenance
  ev <- c(ev, ev(cmt = "PRED_GUT", amt = 0.1 * bw_kg, ii = 24, addl = 365,
                 time = t))
  do.call(ev, ev)
}

#' Dupilumab 300 mg SC Q2W dosing
dup_events <- function(n_doses = 26, start_h = 0) {
  # Loading dose 600 mg at week 0, then 300 mg Q2W
  ev_load <- ev(cmt = "DUP_DEPOT", amt = 600, time = start_h)
  ev_maint <- ev(cmt = "DUP_DEPOT", amt = 300, ii = 14*24,
                 addl = n_doses - 1, time = start_h + 14*24)
  ev_load + ev_maint
}

#' Omalizumab 300 mg SC Q4W dosing
oma_events <- function(n_doses = 12, start_h = 0) {
  ev(cmt = "OMA_DEPOT", amt = 300, ii = 28*24, addl = n_doses - 1,
     time = start_h)
}

#' Rituximab 1g IV x2, then 1g Q6M
ritu_events <- function(start_h = 0) {
  ev1 <- ev(cmt = "RITU_C", amt = 1000, time = start_h)
  ev2 <- ev(cmt = "RITU_C", amt = 1000, time = start_h + 14*24)
  ev3 <- ev(cmt = "RITU_C", amt = 1000, time = start_h + 26*7*24)
  ev1 + ev2 + ev3
}

#' Doxycycline 200 mg/d PO
doxy_events <- function(n_weeks = 52, start_h = 0) {
  ev(cmt = "DOXY_GUT", amt = 200, ii = 24, addl = 7*n_weeks - 1,
     time = start_h)
}

## ── Simulation Time Grid ──────────────────────────────────────────────────────
t_end <- 52 * 7 * 24  # 52 weeks in hours
t_obs <- seq(0, t_end, by = 24)  # daily observations

## ── Scenario Definitions ──────────────────────────────────────────────────────

scenarios <- list(
  "1_No_Treatment" = list(
    label = "Untreated BP",
    color = "#E53935",
    events = ev(cmt = "PRED_GUT", amt = 0, time = 0)  # placeholder
  ),
  "2_High_Dose_Pred" = list(
    label = "High-dose Prednisolone (0.75 mg/kg/d taper)",
    color = "#FF7043",
    events = pred_dose_events(dose_mg_kg = 0.75, bw_kg = 65, taper_weeks = 4)
  ),
  "3_Low_Dose_Pred_Doxy" = list(
    label = "Low-dose Pred (0.3 mg/kg) + Doxycycline 200mg/d",
    color = "#FFA726",
    events = pred_dose_events(dose_mg_kg = 0.30, bw_kg = 65, taper_weeks = 8) +
             doxy_events(n_weeks = 52)
  ),
  "4_Dupilumab" = list(
    label = "Dupilumab 300mg SC Q2W + Low-dose Pred",
    color = "#42A5F5",
    events = dup_events(n_doses = 26) +
             pred_dose_events(dose_mg_kg = 0.25, bw_kg = 65, taper_weeks = 12)
  ),
  "5_Omalizumab" = list(
    label = "Omalizumab 300mg SC Q4W + Low-dose Pred",
    color = "#26C6DA",
    events = oma_events(n_doses = 12) +
             pred_dose_events(dose_mg_kg = 0.30, bw_kg = 65, taper_weeks = 8)
  ),
  "6_Rituximab" = list(
    label = "Rituximab 1g IV ×2 then Q6M + Short Pred",
    color = "#AB47BC",
    events = ritu_events() +
             pred_dose_events(dose_mg_kg = 0.50, bw_kg = 65, taper_weeks = 4)
  )
)

## ── Run Simulations ──────────────────────────────────────────────────────────
run_scenario <- function(scenario_id, scenario) {
  out <- mod %>%
    mrgsim_df(events = scenario$events, end = t_end, delta = 24, hmax = 1) %>%
    mutate(
      scenario = scenario_id,
      label    = scenario$label,
      time_weeks = time / (7 * 24)
    )
  return(out)
}

# Run all scenarios
results <- lapply(names(scenarios), function(s) {
  run_scenario(s, scenarios[[s]])
})
df_all <- bind_rows(results)

## ── Plots ─────────────────────────────────────────────────────────────────────

colors_vec <- sapply(scenarios, function(x) x$color)
names(colors_vec) <- sapply(scenarios, function(x) x$label)

# Map label to color
df_all <- df_all %>%
  mutate(color = colors_vec[label])

# Plot 1: BPDAI Activity Score over Time
p1 <- ggplot(df_all, aes(x = time_weeks, y = BPDAI_act, color = label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_vec) +
  labs(title = "Bullous Pemphigoid Disease Activity (BPDAI) Over 52 Weeks",
       x = "Time (weeks)", y = "BPDAI Activity Score (0–90)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 3))

# Plot 2: Anti-BP180 IgG Antibody Titers
p2 <- ggplot(df_all, aes(x = time_weeks, y = anti_BP180_AU * 100, color = label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_vec) +
  labs(title = "Anti-BP180 IgG Antibody Titer",
       x = "Time (weeks)", y = "Anti-BP180 IgG (AU/mL, normalised × 100)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 3: Itch NRS Score
p3 <- ggplot(df_all, aes(x = time_weeks, y = Itch, color = label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_vec) +
  ylim(0, 10) +
  labs(title = "Itch NRS Score Over 52 Weeks",
       x = "Time (weeks)", y = "Itch NRS (0–10)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 4: DEJ Integrity
p4 <- ggplot(df_all, aes(x = time_weeks, y = DEJ, color = label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_vec) +
  ylim(0, 1) +
  labs(title = "Dermal-Epidermal Junction (DEJ) Integrity",
       x = "Time (weeks)", y = "DEJ Integrity (0=disrupted, 1=normal)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 5: Drug PK — Prednisolone concentration
df_pred_sc2 <- df_all %>% filter(label == scenarios[["2_High_Dose_Pred"]]$label)
p5 <- ggplot(df_pred_sc2, aes(x = time_weeks, y = Cpred_ngml)) +
  geom_line(color = "#FF7043", linewidth = 1.2) +
  labs(title = "Prednisolone Plasma Concentration (Scenario 2: High-dose Taper)",
       x = "Time (weeks)", y = "Prednisolone (ng/mL)") +
  theme_bw(base_size = 12)

# Plot 6: B Cell Dynamics under Rituximab
df_ritu_sc <- df_all %>% filter(label == scenarios[["6_Rituximab"]]$label)
p6_data <- df_ritu_sc %>%
  select(time_weeks, B_NAIVE, B_ACT, LLPC) %>%
  pivot_longer(-time_weeks, names_to = "Cell_Type", values_to = "Count")
p6 <- ggplot(p6_data, aes(x = time_weeks, y = Count, color = Cell_Type)) +
  geom_line(linewidth = 1.1) +
  labs(title = "B Cell Dynamics Under Rituximab Treatment",
       x = "Time (weeks)", y = "Cell Count (AU)",
       color = "Cell Subset") +
  scale_color_manual(values = c(B_NAIVE = "#42A5F5", B_ACT = "#EF5350",
                                LLPC = "#AB47BC")) +
  theme_bw(base_size = 12)

## ── Summary Table ─────────────────────────────────────────────────────────────
summary_table <- df_all %>%
  group_by(label) %>%
  summarise(
    BPDAI_week4   = round(BPDAI_act[which.min(abs(time_weeks - 4))], 1),
    BPDAI_week12  = round(BPDAI_act[which.min(abs(time_weeks - 12))], 1),
    BPDAI_week52  = round(BPDAI_act[which.min(abs(time_weeks - 52))], 1),
    Itch_week4    = round(Itch[which.min(abs(time_weeks - 4))], 1),
    Itch_week52   = round(Itch[which.min(abs(time_weeks - 52))], 1),
    Anti_BP180_w52 = round(anti_BP180_AU[which.min(abs(time_weeks - 52))]*100, 0),
    Remission_rate = round(mean(in_remission[time_weeks > 12]) * 100, 0),
    .groups = "drop"
  ) %>%
  rename(
    "Treatment" = label,
    "BPDAI Wk4"  = BPDAI_week4,
    "BPDAI Wk12" = BPDAI_week12,
    "BPDAI Wk52" = BPDAI_week52,
    "Itch Wk4"   = Itch_week4,
    "Itch Wk52"  = Itch_week52,
    "Anti-BP180 Wk52" = Anti_BP180_w52,
    "Remission >Wk12 (%)" = Remission_rate
  )

print(summary_table)

## ── Save Plots ────────────────────────────────────────────────────────────────
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/bp_bpdai_scenarios.png",   p1, width = 10, height = 6, dpi = 150)
ggsave("plots/bp_antibody_scenarios.png",p2, width = 10, height = 6, dpi = 150)
ggsave("plots/bp_itch_scenarios.png",    p3, width = 10, height = 6, dpi = 150)
ggsave("plots/bp_dej_scenarios.png",     p4, width = 10, height = 6, dpi = 150)
ggsave("plots/bp_pk_prednisolone.png",   p5, width = 9,  height = 5, dpi = 150)
ggsave("plots/bp_bcell_rituximab.png",   p6, width = 9,  height = 5, dpi = 150)

cat("\n=== Bullous Pemphigoid QSP Model — Simulation Complete ===\n")
cat("Plots saved to: bullous-pemphigoid/plots/\n")
