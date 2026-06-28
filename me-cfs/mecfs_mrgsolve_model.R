## ============================================================
## ME/CFS QSP Model — mrgsolve ODE System
## Myalgic Encephalomyelitis / Chronic Fatigue Syndrome
## ============================================================
## Reference: Naviaux et al. (2016) PNAS; Fluge et al. (2017);
##            Tomas et al. (2017) PNAS; Montoya et al. (2017) PNAS
##            Davis et al. (2023); Sweetman et al. (2023)
##
## Compartments (25 ODE state variables):
##  Immune:  V, IFN, NK, Tex, AutoAb
##  Cytokines: IL6, TNFa, IFNg, NLRP3state
##  MCAS:    MC_act, Histamine
##  HPA:     CRH, Cortisol
##  ANS:     NE_plasma, HRV_index
##  Mito:    PDH_act, ATP_state, ROS_state
##  CNS:     Neuro_inf, Cog_func
##  PEM:     PEM_sens, Fatigue
##  Drug PK: LDN_cp, Pyr_cp, Rit_cp, NADpool
## ============================================================

library(mrgsolve)

mecfs_model_code <- '
$PROB ME/CFS QSP Model — Multi-system ODE Model
Mechanistic QSP model for Myalgic Encephalomyelitis/Chronic Fatigue Syndrome

$PARAM
// --- Disease Parameters ---
kV_prod     = 0.0005  // Viral reactivation rate (latent pool)
kV_clear    = 0.15    // NK-mediated viral clearance rate
kV_IFNclear = 0.10    // IFN-γ mediated viral clearance
V_max       = 10.0    // Max viral load (normalized)

// --- NK Cell Dynamics ---
kNK_prod     = 0.05   // NK cell baseline production
kNK_decay    = 0.05   // NK cell natural decay
kNK_exhaust  = 0.02   // NK exhaustion by viral load
kNK_IFN      = 0.008  // NK exhaustion by chronic IFN

// --- T Cell Exhaustion ---
kTex_induce  = 0.03   // Rate of CD8+ T cell exhaustion induction
kTex_recover = 0.005  // Recovery from exhaustion (slow)

// --- Autoantibody Dynamics ---
kAutoAb_prod = 0.002  // Autoantibody production (B cell mediated)
kAutoAb_decay = 0.015 // Autoantibody half-life decay
kAutoAb_mol  = 0.05   // Molecular mimicry induction coefficient

// --- Cytokine Dynamics ---
kIL6_prod    = 0.10   // IL-6 production rate (inflammation)
kIL6_decay   = 0.25   // IL-6 degradation
kIL6_V       = 0.08   // IL-6 production by viral antigen
kTNFa_prod   = 0.08   // TNF-α production
kTNFa_decay  = 0.30   // TNF-α degradation
kIFNg_prod   = 0.06   // IFN-γ production (Th1/NK)
kIFNg_decay  = 0.20   // IFN-γ degradation
kIFNg_V      = 0.10   // IFN-γ induction by viral load
kNLRP3_act   = 0.05   // NLRP3 activation rate
kNLRP3_decay = 0.15   // NLRP3 return to baseline

// --- IFN (Type I) ---
kIFN_prod    = 0.15   // Type I IFN production rate
kIFN_decay   = 0.30   // Type I IFN degradation

// --- Mast Cell Activation ---
kMC_act      = 0.04   // Mast cell activation rate
kMC_decay    = 0.08   // MC activation decay
kHist_prod   = 0.20   // Histamine release per MC activation
kHist_decay  = 0.35   // Histamine degradation

// --- HPA Axis ---
kCRH_prod    = 0.10   // CRH production (stress)
kCRH_decay   = 0.25   // CRH degradation
kCort_prod   = 0.15   // Cortisol production (ACTH driven)
kCort_decay  = 0.12   // Cortisol degradation
Cort_base    = 1.0    // Normal baseline cortisol (normalized)
IL6_GR_IC50  = 1.5    // IL-6 concentration for 50% GR resistance

// --- ANS ---
kNE_prod     = 0.08   // NE production (sympathetic)
kNE_decay    = 0.20   // NE reuptake/degradation
kHRV_restore = 0.05   // HRV restoration rate
kHRV_suppress = 0.15  // HRV suppression by NE
NE_base      = 1.0    // Normal NE baseline
HRV_base     = 1.0    // Normal HRV baseline

// --- Mitochondrial Dynamics ---
kPDH_base    = 1.0    // PDH baseline activity
kPDK1_IFNg   = 0.30   // PDK1 induction by IFN-γ (inhibits PDH)
kPDK1_TNFa   = 0.20   // PDK1 induction by TNF-α
kROS_prod    = 0.15   // ROS production rate (inverse of PDH activity)
kROS_decay   = 0.10   // ROS scavenging (SOD/GSH)
kATP_prod    = 1.2    // ATP production rate (normal)
kATP_decay   = 0.30   // ATP consumption rate
kATP_ROS     = 0.25   // ATP reduction by ROS (ETC damage)

// --- Neuroinflammation ---
kNI_prod     = 0.08   // Neuroinflammation induction rate
kNI_decay    = 0.06   // Neuroinflammation resolution rate
kNI_IL6      = 0.12   // NI driven by IL-6
kNI_TNFa     = 0.10   // NI driven by TNF-α
kCog_impair  = 0.10   // Cognitive impairment by NI
kCog_restore = 0.05   // Cognitive function restoration

// --- PEM & Fatigue ---
kPEM_ATP     = 0.20   // PEM sensitivity from ATP deficit
kPEM_NI      = 0.15   // PEM sensitivity from neuroinflammation
kFat_prod    = 0.15   // Fatigue accumulation rate
kFat_decay   = 0.08   // Fatigue recovery rate
PEM_threshold = 2.0   // PEM threshold (normalized exertion)

// --- Drug PK Parameters ---
// LDN (Low-Dose Naltrexone)
LDN_dose     = 0.0    // LDN daily dose (mg) — set in scenarios
LDN_F        = 0.96   // Bioavailability
LDN_ka       = 1.2    // Absorption rate constant (1/h)
LDN_CL       = 9.7    // Clearance (L/h)
LDN_V1       = 28.0   // Central volume (L)
LDN_Q        = 5.0    // Intercompartmental clearance
LDN_V2       = 45.0   // Peripheral volume (L)
LDN_Emax_TLR4 = 0.7   // Max TLR4 inhibition by LDN
LDN_EC50_TLR4 = 0.05  // LDN Cp for 50% TLR4 inhibition (ng/mL)

// Pyridostigmine
Pyr_dose     = 0.0    // Pyridostigmine dose (mg) — set in scenarios
Pyr_F        = 0.20   // Oral bioavailability
Pyr_ka       = 0.8    // Absorption rate (1/h)
Pyr_CL       = 10.0   // Clearance (L/h)
Pyr_V        = 50.0   // Volume of distribution (L)
Pyr_Emax_ANS = 0.6    // Max ANS restoration
Pyr_EC50     = 30.0   // Cp50 for AChE inhibition (ng/mL)

// Rituximab (simplified 1-compartment)
Rit_dose     = 0.0    // Rituximab dose (mg) — IV
Rit_CL       = 0.23   // Clearance (L/day → converted)
Rit_V        = 4.4    // Central volume (L)
Rit_Emax_B   = 0.95   // Max B cell depletion
Rit_EC50     = 0.1    // EC50 for B cell depletion (µg/mL)

// NAD+ Supplementation
NAD_dose     = 0.0    // NAD precursor dose (mg) — set in scenarios
NAD_ka       = 0.5    // Absorption rate (1/h)
NAD_F        = 0.35   // Bioavailability (NMN/NR)
NAD_CL       = 2.0    // Clearance from pool
NAD_V        = 10.0   // Volume (cellular pool equivalent)
NAD_Emax_mito = 0.5   // Max mitochondrial restoration
NAD_EC50     = 200.0  // EC50 for mito restoration

$INIT
// Immune compartments (normalized to 1.0 = healthy baseline)
V        = 0.01  // Viral latency (low baseline)
IFN      = 0.1   // Type I IFN (slight elevation in ME/CFS)
NK       = 0.6   // NK cell count (reduced at baseline in ME/CFS)
Tex      = 0.4   // T cell exhaustion index (elevated)
AutoAb   = 0.5   // Autoantibody level (elevated)

// Cytokines (normalized, 1.0 = healthy)
IL6      = 1.3   // IL-6 (mildly elevated)
TNFa     = 1.2   // TNF-α (mildly elevated)
IFNg     = 1.4   // IFN-γ (elevated)
NLRP3state = 0.8 // NLRP3 activation state

// MCAS
MC_act   = 0.6   // Mast cell activation (elevated)
Histamine = 0.8  // Histamine level (elevated)

// HPA Axis
CRH      = 0.8   // CRH (variable, can be ↑ or ↓)
Cortisol = 0.6   // Cortisol (hypocortisolism)

// ANS
NE_plasma = 1.4  // Norepinephrine (elevated)
HRV_index = 0.5  // HRV index (reduced)

// Mitochondrial
PDH_act  = 0.4   // PDH activity (↓ in ME/CFS)
ATP_state = 0.4  // ATP level (reduced)
ROS_state = 1.6  // ROS (elevated)

// CNS
Neuro_inf = 0.8  // Neuroinflammation index
Cog_func  = 0.5  // Cognitive function (impaired)

// PEM & Fatigue
PEM_sens = 1.8   // PEM sensitivity (elevated)
Fatigue  = 0.7   // Fatigue score (high, normalized to 0-1)

// Drug PK
LDN_cp   = 0.0   // LDN central plasma (ng/mL)
Pyr_cp   = 0.0   // Pyridostigmine central plasma (ng/mL)
Rit_cp   = 0.0   // Rituximab central plasma (µg/mL)
NADpool  = 0.3   // NAD+ pool (↓ in ME/CFS, normalized)

$ODE
// ============================================================
// Drug PK ODEs
// ============================================================
// LDN 2-compartment PK
double LDN_dose_rate = LDN_dose * LDN_F * LDN_ka / LDN_V1;
dxdt_LDN_cp = LDN_dose_rate - (LDN_CL/LDN_V1)*LDN_cp - (LDN_Q/LDN_V1)*LDN_cp;
// (simplified: peripheral not tracked, absorbed from depot externally)

// Pyridostigmine 1-compartment
dxdt_Pyr_cp = (Pyr_dose * Pyr_F * Pyr_ka / Pyr_V) - (Pyr_CL/Pyr_V)*Pyr_cp;

// Rituximab 1-compartment
dxdt_Rit_cp = (Rit_dose / Rit_V) - (Rit_CL/Rit_V)*Rit_cp;

// NAD+ Pool
dxdt_NADpool = (NAD_dose * NAD_F * NAD_ka / NAD_V) - (NAD_CL/NAD_V)*NADpool;

// ============================================================
// Drug Effect Calculations (PD)
// ============================================================
// LDN TLR4 inhibition (sigmoid Emax)
double LDN_TLR4_inh = LDN_Emax_TLR4 * LDN_cp / (LDN_EC50_TLR4 + LDN_cp + 1e-9);

// Pyridostigmine ANS restoration
double Pyr_ANS_rest = Pyr_Emax_ANS * Pyr_cp / (Pyr_EC50 + Pyr_cp + 1e-9);

// Rituximab B cell depletion → AutoAb reduction
double Rit_B_depl = Rit_Emax_B * Rit_cp / (Rit_EC50 + Rit_cp + 1e-9);

// NAD+ mitochondrial restoration
double NAD_mito_rest = NAD_Emax_mito * NADpool / (NAD_EC50/1000.0 + NADpool + 1e-9);

// ============================================================
// IMMUNE COMPARTMENTS
// ============================================================
// Viral Load (latent + reactivation)
dxdt_V = kV_prod * V * (1.0 - V/V_max)
         - kV_clear * NK * V
         - kV_IFNclear * IFNg * V;

// Type I IFN (driven by viral load and latent reservoir)
dxdt_IFN = kIFN_prod * (V + 0.1) * (1.0 - LDN_TLR4_inh)
           - kIFN_decay * IFN;

// NK Cell Count (depleted/exhausted in ME/CFS)
dxdt_NK = kNK_prod
          - kNK_decay * NK
          - kNK_exhaust * V * NK
          - kNK_IFN * IFN * NK;

// CD8+ T Cell Exhaustion Index (higher = more exhausted)
dxdt_Tex = kTex_induce * IFN * V
           - kTex_recover * Tex * Cortisol;

// Autoantibody Level (β2AR + M1R autoantibodies)
double AutoAb_prod_mol_mimicry = kAutoAb_mol * V * (1.0 - Rit_B_depl);
dxdt_AutoAb = kAutoAb_prod * (1.0 - Rit_B_depl) + AutoAb_prod_mol_mimicry
              - kAutoAb_decay * AutoAb;

// ============================================================
// CYTOKINE DYNAMICS
// ============================================================
// IL-6 (driven by viral Ag, NK dysfunction, LDN-sensitive TLR4)
dxdt_IL6 = kIL6_prod * (V + 0.2) * (1.0 - LDN_TLR4_inh * 0.5)
           + kIL6_V * (1.0 - NK)   // NK dysfunction → less IL-6 suppression
           - kIL6_decay * IL6
           - kIL6_decay * 0.5 * Cortisol * IL6; // GR suppression

// TNF-α
dxdt_TNFa = kTNFa_prod * (V + 0.15) * (1.0 - LDN_TLR4_inh * 0.4)
            - kTNFa_decay * TNFa
            - kTNFa_decay * 0.3 * Cortisol * TNFa;

// IFN-γ (from Th1 cells and NK cells)
dxdt_IFNg = kIFNg_prod * (1.0 - NK) * V  // impaired NK → less NK-IFN, but Th1 compensates
            + kIFNg_V * V
            - kIFNg_decay * IFNg;

// NLRP3 Inflammasome State
dxdt_NLRP3state = kNLRP3_act * (ROS_state - 1.0) * (IL6 - 1.0 + 1e-9)
                  - kNLRP3_decay * NLRP3state;

// ============================================================
// MAST CELL ACTIVATION & MCAS
// ============================================================
// MC Activation State (driven by IL-33, complement, autoantigens)
dxdt_MC_act = kMC_act * (IL6 * 0.5 + V * 0.5 + AutoAb * 0.3)
              - kMC_decay * MC_act
              - kMC_decay * 0.5 * Cortisol * MC_act; // Cortisol suppresses MC

// Histamine (released by activated MC)
dxdt_Histamine = kHist_prod * MC_act
                 - kHist_decay * Histamine;

// ============================================================
// HPA AXIS
// ============================================================
// CRH (driven by stress, cytokines; suppressed by cortisol feedback)
dxdt_CRH = kCRH_prod * (1.0 + Neuro_inf * 0.3)
            - kCRH_decay * CRH
            - kCRH_decay * Cortisol * CRH * 0.5; // Cortisol feedback

// Cortisol (driven by CRH/ACTH; suppressed by GR resistance)
double GR_resistance = IL6 / (IL6_GR_IC50 + IL6); // GR resistance by IL-6
dxdt_Cortisol = kCort_prod * CRH * (1.0 - GR_resistance * 0.5)
                - kCort_decay * Cortisol;

// ============================================================
// AUTONOMIC NERVOUS SYSTEM
// ============================================================
// Norepinephrine (sympathetic activation from AutoAb β2AR, POTS)
dxdt_NE_plasma = kNE_prod * (AutoAb * 0.5 + 0.5)
                 - kNE_decay * NE_plasma
                 + Pyr_ANS_rest * (NE_base - NE_plasma) * 0.3; // Pyr normalizes

// HRV Index (reduced by NE, reduced by low vagal tone; Pyr restores)
dxdt_HRV_index = kHRV_restore * (1.0 - Pyr_ANS_rest * 0.0) + kHRV_restore * Pyr_ANS_rest
                 - kHRV_suppress * NE_plasma * HRV_index
                 + kHRV_restore * Pyr_ANS_rest;

// ============================================================
// MITOCHONDRIAL FUNCTION & ENERGY METABOLISM
// ============================================================
// PDH Activity (inhibited by PDK1, which is induced by IFN-γ and TNF-α)
double PDK1_activity = kPDK1_IFNg * IFNg + kPDK1_TNFa * TNFa;
dxdt_PDH_act = -PDK1_activity * PDH_act        // PDK1 inhibits PDH
               + 0.02 * (1.0 - PDH_act)         // Basal recovery
               + NAD_mito_rest * (1.0 - PDH_act) * 0.3; // NAD+ restores PDH

// ATP State (produced by OxPhos, reduced by ROS damage)
dxdt_ATP_state = kATP_prod * PDH_act
                 - kATP_decay * ATP_state
                 - kATP_ROS * ROS_state * ATP_state
                 + Mito_support_effect * 0.3
                 + NAD_mito_rest * 0.2;

// ROS State (elevated when ETC is impaired)
dxdt_ROS_state = kROS_prod * (1.0 - PDH_act)  // More ROS when PDH impaired
                 + kROS_prod * 0.5 * (IL6 - 1.0 + 1e-9)
                 - kROS_decay * ROS_state
                 - kROS_decay * NAD_mito_rest * 0.3;  // NAD reduces ROS

// ============================================================
// NEUROINFLAMMATION & COGNITIVE FUNCTION
// ============================================================
// Neuroinflammation Index (microglia activation, TSPO signal)
dxdt_Neuro_inf = kNI_prod * (kNI_IL6 * IL6 + kNI_TNFa * TNFa)
                 - kNI_decay * Neuro_inf
                 - kNI_decay * LDN_TLR4_inh * 0.5 * Neuro_inf // LDN reduces NI
                 + kNI_prod * 0.3 * Histamine; // Histamine promotes NI

// Cognitive Function (impaired by NI, ROS; improved by cortisol, treatment)
dxdt_Cog_func = kCog_restore * Cortisol * 0.3
                - kCog_impair * Neuro_inf
                - kCog_impair * 0.5 * ROS_state
                + kCog_restore * ATP_state * 0.2
                - kCog_impair * 0.3 * (1.0 - ATP_state);

// ============================================================
// PEM SENSITIVITY & FATIGUE
// ============================================================
// PEM Sensitivity (lower threshold = easier to trigger PEM)
dxdt_PEM_sens = kPEM_ATP * (1.0 - ATP_state)    // Low ATP → high PEM
                + kPEM_NI * Neuro_inf
                - 0.05 * PEM_sens * Cortisol     // Cortisol reduces PEM slightly
                + 0.03 * AutoAb * PEM_sens;      // AutoAb amplifies PEM

// Fatigue Score (accumulated from multiple pathways)
dxdt_Fatigue = kFat_prod * (1.0 - ATP_state) * 0.4
               + kFat_prod * Neuro_inf * 0.3
               + kFat_prod * (IL6 - 1.0) * 0.2
               + kFat_prod * (1.0 - HRV_index) * 0.2
               - kFat_decay * Fatigue * Cortisol
               - kFat_decay * 0.5 * Fatigue;

$CAPTURE
// Key model outputs
V IFN NK Tex AutoAb
IL6 TNFa IFNg NLRP3state
MC_act Histamine
CRH Cortisol
NE_plasma HRV_index
PDH_act ATP_state ROS_state
Neuro_inf Cog_func
PEM_sens Fatigue
LDN_cp Pyr_cp Rit_cp NADpool
GR_resistance PDK1_activity LDN_TLR4_inh Pyr_ANS_rest

$GLOBAL
double Mito_support_effect = 0.0; // placeholder for CoQ10/L-carnitine effect

$MAIN
// Initialize mito support effect from dose (simplified)
Mito_support_effect = (NAD_dose > 0) ? 0.3 : 0.0;

$SET delta = 0.1, end = 365, hmax = 0.5
'

# ============================================================
# Build & Compile Model
# ============================================================

model <- mrgsolve::mcode("mecfs_qsp", mecfs_model_code)

# ============================================================
# SCENARIO DEFINITIONS
# ============================================================

# Shared event table: daily dosing for 12 months
sim_duration <- 365  # days
dose_start   <- 7    # start treatment at day 7

scenarios <- list(

  # ---- Scenario 1: No Treatment (Natural History) ----
  s1_no_treatment = list(
    name   = "No Treatment (Natural History)",
    params = list(
      LDN_dose = 0, Pyr_dose = 0, Rit_dose = 0,
      NAD_dose = 0
    ),
    events = mrgsolve::ev(time = 0, amt = 0, cmt = 0)
  ),

  # ---- Scenario 2: Low-Dose Naltrexone (LDN 4.5 mg/day) ----
  s2_LDN = list(
    name   = "LDN 4.5 mg/day",
    params = list(
      LDN_dose = 4.5, Pyr_dose = 0, Rit_dose = 0,
      NAD_dose = 0
    ),
    events = mrgsolve::ev(
      ID   = 1,
      time = dose_start + seq(0, sim_duration - dose_start, by = 1),
      amt  = 4.5,
      cmt  = "LDN_cp",
      rate = -2
    )
  ),

  # ---- Scenario 3: Pyridostigmine (30 mg TID) ----
  s3_pyridostigmine = list(
    name   = "Pyridostigmine 30 mg TID (ANS support)",
    params = list(
      LDN_dose = 0, Pyr_dose = 30, Rit_dose = 0,
      NAD_dose = 0
    ),
    events = mrgsolve::ev(
      ID   = 1,
      time = rep(dose_start + seq(0, sim_duration - dose_start, by = 1), each = 3) +
             rep(c(0, 8, 16)/24, times = sim_duration - dose_start + 1),
      amt  = 30,
      cmt  = "Pyr_cp",
      rate = -2
    )
  ),

  # ---- Scenario 4: Rituximab (1000 mg IV x 2, 2-week apart) ----
  s4_rituximab = list(
    name   = "Rituximab 1g IV (Weeks 1 & 3)",
    params = list(
      LDN_dose = 0, Pyr_dose = 0, Rit_dose = 1000,
      NAD_dose = 0
    ),
    events = mrgsolve::ev(
      ID   = 1,
      time = c(dose_start, dose_start + 14),
      amt  = c(1000, 1000),
      cmt  = "Rit_cp",
      rate = c(10, 10)  # infuse over ~100h
    )
  ),

  # ---- Scenario 5: NAD+ Supplementation (500 mg/day NMN) ----
  s5_NAD = list(
    name   = "NAD+ Precursor (NMN 500 mg/day)",
    params = list(
      LDN_dose = 0, Pyr_dose = 0, Rit_dose = 0,
      NAD_dose = 500
    ),
    events = mrgsolve::ev(
      ID   = 1,
      time = dose_start + seq(0, sim_duration - dose_start, by = 1),
      amt  = 500,
      cmt  = "NADpool",
      rate = -2
    )
  ),

  # ---- Scenario 6: Combination (LDN + NAD+ + Pyridostigmine) ----
  s6_combination = list(
    name   = "Combination: LDN + NAD+ + Pyridostigmine",
    params = list(
      LDN_dose = 4.5, Pyr_dose = 30, Rit_dose = 0,
      NAD_dose = 500
    ),
    events = mrgsolve::ev(
      ID   = 1,
      time = c(
        dose_start + seq(0, sim_duration - dose_start, by = 1),     # LDN
        rep(dose_start + seq(0, sim_duration - dose_start, by = 1), each = 3) + # Pyr TID
          rep(c(0, 8, 16)/24, times = sim_duration - dose_start + 1),
        dose_start + seq(0, sim_duration - dose_start, by = 1)      # NAD
      ),
      amt  = c(
        rep(4.5, sim_duration - dose_start + 1),    # LDN
        rep(30, (sim_duration - dose_start + 1)*3), # Pyr
        rep(500, sim_duration - dose_start + 1)     # NAD
      ),
      cmt  = c(
        rep("LDN_cp", sim_duration - dose_start + 1),
        rep("Pyr_cp", (sim_duration - dose_start + 1)*3),
        rep("NADpool", sim_duration - dose_start + 1)
      ),
      rate = -2
    )
  )
)

# ============================================================
# RUN SIMULATIONS
# ============================================================

run_scenario <- function(model, scenario, n_id = 1) {
  model %>%
    mrgsolve::param(scenario$params) %>%
    mrgsolve::ev(scenario$events) %>%
    mrgsolve::mrgsim(nid = n_id, end = sim_duration, delta = 0.5) %>%
    as.data.frame() %>%
    dplyr::mutate(Scenario = scenario$name)
}

# Run all scenarios (requires dplyr + mrgsolve)
if (interactive()) {
  library(dplyr)
  library(ggplot2)

  sim_results <- lapply(scenarios, function(s) {
    tryCatch(
      run_scenario(model, s),
      error = function(e) { message("Error in ", s$name, ": ", e$message); NULL }
    )
  })

  all_results <- do.call(rbind, Filter(Negate(is.null), sim_results))

  # ============================================================
  # SUMMARY PLOTS
  # ============================================================

  # 1. Fatigue Score over time
  p_fatigue <- ggplot(all_results, aes(x = time, y = Fatigue, color = Scenario)) +
    geom_line(size = 0.8) +
    labs(title = "ME/CFS Fatigue Score Over Time",
         subtitle = "Normalized fatigue (0 = healthy, 1 = maximal fatigue)",
         x = "Time (days)", y = "Fatigue Score") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 2. ATP State (Mitochondrial Energy)
  p_atp <- ggplot(all_results, aes(x = time, y = ATP_state, color = Scenario)) +
    geom_line(size = 0.8) +
    labs(title = "Cellular ATP State Over Time",
         subtitle = "Normalized ATP (1.0 = healthy; <0.5 = energy crisis)",
         x = "Time (days)", y = "ATP State (normalized)") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 3. PEM Sensitivity
  p_pem <- ggplot(all_results, aes(x = time, y = PEM_sens, color = Scenario)) +
    geom_line(size = 0.8) +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray50") +
    labs(title = "Post-Exertional Malaise (PEM) Sensitivity",
         subtitle = "Higher = more sensitive; dashed = baseline threshold",
         x = "Time (days)", y = "PEM Sensitivity Index") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 4. Neuroinflammation
  p_ni <- ggplot(all_results, aes(x = time, y = Neuro_inf, color = Scenario)) +
    geom_line(size = 0.8) +
    labs(title = "Neuroinflammation Index Over Time",
         x = "Time (days)", y = "Neuroinflammation Index") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 5. Autoantibody Levels
  p_autoab <- ggplot(all_results, aes(x = time, y = AutoAb, color = Scenario)) +
    geom_line(size = 0.8) +
    labs(title = "Autoantibody Levels (β2AR + M1R AutoAb)",
         x = "Time (days)", y = "Autoantibody Level (normalized)") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 6. HRV Index (ANS)
  p_hrv <- ggplot(all_results, aes(x = time, y = HRV_index, color = Scenario)) +
    geom_line(size = 0.8) +
    labs(title = "Heart Rate Variability (HRV) Index",
         subtitle = "1.0 = healthy; <0.5 = autonomic dysfunction",
         x = "Time (days)", y = "HRV Index (normalized)") +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # 7. Cytokine Profile (IL-6, TNF-α, IFN-γ)
  p_cytokines <- all_results %>%
    dplyr::select(time, Scenario, IL6, TNFa, IFNg) %>%
    tidyr::pivot_longer(cols = c(IL6, TNFa, IFNg), names_to = "Cytokine", values_to = "Level") %>%
    ggplot(aes(x = time, y = Level, color = Scenario, linetype = Cytokine)) +
    geom_line(size = 0.7) +
    labs(title = "Cytokine Dynamics Over Time",
         x = "Time (days)", y = "Cytokine Level (normalized)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  # Combine plots
  library(gridExtra)
  combined_plot <- gridExtra::arrangeGrob(
    p_fatigue, p_atp, p_pem, p_ni, p_autoab, p_hrv,
    ncol = 2, nrow = 3,
    top = "ME/CFS QSP Model — 6-Scenario Simulation Results"
  )

  # Print summary table at 6 months
  summary_6mo <- all_results %>%
    dplyr::filter(abs(time - 180) < 1) %>%
    dplyr::group_by(Scenario) %>%
    dplyr::summarise(
      Fatigue      = mean(Fatigue),
      ATP_state    = mean(ATP_state),
      PEM_sens     = mean(PEM_sens),
      Neuro_inf    = mean(Neuro_inf),
      AutoAb       = mean(AutoAb),
      HRV_index    = mean(HRV_index),
      Cortisol     = mean(Cortisol),
      .groups = "drop"
    )

  cat("\n=== ME/CFS QSP Model: 6-Month Outcome Summary ===\n")
  print(summary_6mo, digits = 3)

  invisible(list(
    model   = model,
    results = all_results,
    plots   = list(fatigue = p_fatigue, atp = p_atp, pem = p_pem,
                   ni = p_ni, autoab = p_autoab, hrv = p_hrv,
                   cytokines = p_cytokines),
    summary = summary_6mo
  ))
}
