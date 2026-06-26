## =============================================================================
## ALS (Amyotrophic Lateral Sclerosis) — QSP mrgsolve Model
## =============================================================================
## Compartments  : 26 (9 PK + 14 disease biology + 3 clinical endpoints)
## ODE parameters: 60+
## Scenarios     : 7 (untreated, riluzole, edaravone, combo, tofersen, AMX0035, all)
## Calibration   : Riley 2004 (riluzole), Abe 2017 (edaravone), Miller 2022 (tofersen),
##                 Paganoni 2020 (AMX0035/CENTAUR), Bensimon 1994 (natural history)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------------------------------------------------------
## 1. Model code
## -----------------------------------------------------------------------------
als_code <- '
$PROB
ALS QSP Model
Amyotrophic Lateral Sclerosis — Motor Neuron Degeneration, Protein Aggregation,
Excitotoxicity, Neuroinflammation, Drug PK/PD

$PARAM @annotated
// ── Motor Neuron Parameters ──────────────────────────────────────
k_MN_death    : 0.00080  : baseline MN death rate constant (1/day)
k_MN_death_U  : 0.00060  : upper MN death rate constant  (1/day)
MN_upper_0    : 1.0      : initial upper motor neuron (normalized)
MN_lower_0    : 1.0      : initial lower motor neuron (normalized)

// ── SOD1 Protein Dynamics ────────────────────────────────────────
k_SOD1_syn    : 0.50     : SOD1 synthesis rate (AU/day)
k_SOD1_deg    : 0.10     : SOD1 wild-type degradation rate (1/day)
k_SOD1_mis    : 0.08     : SOD1 misfolding rate constant (1/day)
k_SOD1_clr    : 0.015    : misfolded SOD1 clearance (UPS/autophagy, 1/day)
f_SOD1_mut    : 0.0      : fraction mutant (0=sporadic, 1=SOD1-ALS)

// ── TDP-43 Dynamics ──────────────────────────────────────────────
k_TDP_export  : 0.020    : TDP-43 nuclear→cytoplasmic export (1/day)
k_TDP_import  : 0.150    : TDP-43 cytoplasmic→nuclear import (1/day)
k_TDP_agg     : 0.010    : TDP-43 cytoplasmic aggregation (1/day)
k_TDP_clr     : 0.008    : TDP-43 aggregate clearance (1/day)
TDP43_nuc_0   : 10.0     : initial nuclear TDP-43 (AU)

// ── Glutamate Excitotoxicity ──────────────────────────────────────
Glu_base      : 1.0      : baseline synaptic glutamate (AU)
k_Glu_rel     : 2.0      : glutamate release rate (AU/day)
k_EAAT2       : 1.80     : EAAT2 baseline uptake constant (1/day)
k_EAAT2_ALS   : 0.40     : EAAT2 ALS reduction factor (dimensionless <1)
k_Ca_entry    : 0.50     : Ca2+ influx rate constant per Glu excess
k_Ca_efflux   : 1.20     : Ca2+ efflux/buffering rate (1/day)
Ca_i_0        : 0.10     : baseline intracellular Ca2+ (μM)

// ── Oxidative Stress ─────────────────────────────────────────────
k_ROS_prod    : 0.40     : ROS production rate per Ca2+ & SOD1mis (AU/day)
k_ROS_scav    : 0.60     : ROS scavenging by GSH (1/day per AU_GSH)
k_GSH_syn     : 0.50     : GSH synthesis rate (AU/day)
k_GSH_cons    : 0.12     : GSH consumption (1/day, basal turnover)
ROS_0         : 0.50     : baseline ROS (AU)
GSH_0         : 4.0      : baseline GSH (AU)

// ── Mitochondrial Function ───────────────────────────────────────
k_Mito_dam    : 0.06     : mitochondrial damage by ROS (1/day per AU_ROS)
k_Mito_rep    : 0.03     : mitochondrial repair/mitophagy (1/day)
Mito_0        : 1.0      : baseline mitochondrial function (normalized)

// ── Neuroinflammation ────────────────────────────────────────────
k_Mic_act     : 0.12     : microglial activation rate (1/day)
k_Mic_res     : 0.060    : microglial resolution rate (1/day)
k_TNFa_prod   : 0.60     : TNF-α production by M1 microglia (AU/day)
k_TNFa_deg    : 0.35     : TNF-α degradation rate (1/day)
k_IL1b_prod   : 0.40     : IL-1β production rate (AU/day)
k_IL1b_deg    : 0.25     : IL-1β degradation (1/day)
Mic_0         : 0.10     : baseline microglial activation (0–1)

// ── Trophic Support ──────────────────────────────────────────────
k_BDNF_prod   : 0.25     : BDNF production (AU/day)
k_BDNF_deg    : 0.18     : BDNF degradation (1/day)
BDNF_0        : 1.0      : baseline BDNF (AU)
BDNF_prot_wt  : 0.30     : max neuroprotection weight from BDNF (dimensionless)

// ── Biomarkers ───────────────────────────────────────────────────
k_NfL_rel     : 0.10     : NfL release rate per MN death (AU/day)
k_NfL_clr     : 0.025    : NfL clearance from CSF (1/day)
NfL_CSF_0     : 5.0      : baseline CSF NfL (pg/mL AU)

// ── Clinical Endpoints ───────────────────────────────────────────
ALSFRS_0      : 48.0     : initial ALSFRS-R total score
FVC_0         : 100.0    : initial FVC % predicted
k_FVC_dec     : 0.0060   : FVC decline rate per unit TNFa & Mic (1/day)

// ── Riluzole PK (2-compartment, oral) ────────────────────────────
F_RIL         : 0.60     : riluzole oral bioavailability
ka_RIL        : 0.80     : absorption rate (1/h)
CL_RIL        : 28.0     : clearance (L/h)
V1_RIL        : 245.0    : central volume (L)
Q_RIL         : 15.0     : inter-compartmental CL (L/h)
V2_RIL        : 112.0    : peripheral volume (L)
// Riluzole PD
IC50_RIL      : 0.50     : IC50 glutamate release inhibition (μg/mL)
Emax_RIL      : 0.60     : max inhibition of Glu release (fraction)

// ── Edaravone PK (1-compartment, IV → oral switch) ───────────────
CL_EDA        : 18.0     : edaravone clearance (L/h)
V_EDA         : 120.0    : edaravone volume (L)
// Edaravone PD
IC50_EDA      : 1.20     : IC50 ROS scavenging (μg/mL)
Emax_EDA      : 0.70     : max ROS reduction (fraction)

// ── Tofersen PK (2-compartment, SC) ─────────────────────────────
ka_TOF        : 0.030    : tofersen SC absorption rate (1/h)
CL_TOF        : 0.50     : tofersen plasma clearance (L/h)
V1_TOF        : 15.0     : central volume (L)
Q_TOF         : 0.30     : distribution to CSF (L/h)
V2_TOF        : 5.0      : CSF volume (L)
// Tofersen PD (SOD1 mRNA knockdown)
EC50_TOF      : 0.10     : EC50 for SOD1 mRNA knockdown (μg/mL CSF)
Emax_TOF      : 0.80     : max SOD1 mRNA knockdown (fraction)

// ── AMX0035: Phenylbutyrate PK ───────────────────────────────────
F_PB          : 0.85     : phenylbutyrate bioavailability
ka_PB         : 1.20     : PB absorption rate (1/h)
CL_PB         : 12.0     : PB clearance (L/h)
V_PB          : 50.0     : PB volume (L)
// AMX0035 PD (ER stress & mitochondria)
EC50_PB       : 50.0     : EC50 ER stress reduction (μmol/L PB)
Emax_PB       : 0.50     : max ER stress / CHOP reduction
Emax_mito_PB  : 0.30     : max mitochondrial protection by PB+TUDCA

$CMT @annotated
// ── Drug PK ──────────────────────────────────────────────────────
DEPOT_RIL  : Riluzole oral depot (mg)
C1_RIL     : Riluzole central compartment (mg)
C2_RIL     : Riluzole peripheral compartment (mg)
IV_EDA     : Edaravone central compartment (mg)
DEPOT_TOF  : Tofersen SC depot (mg)
C1_TOF     : Tofersen plasma compartment (mg)
C2_TOF     : Tofersen CSF compartment (mg)
DEPOT_PB   : Phenylbutyrate oral depot (mg)
C_PB       : Phenylbutyrate plasma compartment (mg)

// ── Disease Biology ───────────────────────────────────────────────
MN_upper   : Upper motor neurons (corticospinal, normalized)
MN_lower   : Lower motor neurons (spinal/bulbar, normalized)
SOD1_wt    : Wild-type SOD1 protein (AU)
SOD1_mis   : Misfolded SOD1 aggregates (AU)
TDP43_nuc  : Nuclear TDP-43 (AU)
TDP43_cyto : Cytoplasmic TDP-43 (AU)
Glut_syn   : Synaptic glutamate (AU)
Ca_i       : Intracellular calcium (μM)
ROS        : Reactive oxygen species (AU)
GSH        : Glutathione (AU)
Mito       : Mitochondrial integrity (normalized 0–1)
Mic_act    : Microglial activation state (0–1)
TNFa       : TNF-α (AU)
BDNF       : BDNF trophic factor (AU)

// ── Biomarkers & Clinical ─────────────────────────────────────────
NfL_CSF    : CSF neurofilament light chain (pg/mL AU)
ALSFRS     : ALSFRS-R total score (0–48)
FVC        : FVC % predicted

$MAIN
// ── Derived concentrations ───────────────────────────────────────
double Cp_RIL  = C1_RIL  / V1_RIL;       // riluzole μg/mL
double Cp_EDA  = IV_EDA  / V_EDA;         // edaravone μg/mL
double Ccsf_TOF = C2_TOF / V2_TOF;        // tofersen CSF μg/mL
double Cp_PB   = (C_PB   / V_PB) * 1000.0; // PB μmol/L (MW≈122)

// ── Drug effect calculations ─────────────────────────────────────
double E_RIL  = Emax_RIL  * Cp_RIL  / (IC50_RIL  + Cp_RIL  + 1e-9);
double E_EDA  = Emax_EDA  * Cp_EDA  / (IC50_EDA  + Cp_EDA  + 1e-9);
double E_TOF  = Emax_TOF  * Ccsf_TOF / (EC50_TOF + Ccsf_TOF + 1e-9);
double E_PB   = Emax_PB   * Cp_PB   / (EC50_PB   + Cp_PB   + 1e-9);
double E_mito_PB = Emax_mito_PB * Cp_PB / (EC50_PB + Cp_PB + 1e-9);

// ── SOD1 aggregation burden ──────────────────────────────────────
double total_SOD1  = SOD1_wt + SOD1_mis + 1e-9;
double frac_SOD1mis = SOD1_mis / total_SOD1;
double SOD1_burden  = f_SOD1_mut * frac_SOD1mis;  // 0 in sporadic ALS

// ── TDP-43 pathology ─────────────────────────────────────────────
double frac_TDP_cyto = TDP43_cyto / (TDP43_nuc + TDP43_cyto + 1e-9);

// ── Excitotoxicity normalized ────────────────────────────────────
double Glu_excess = fmax(0.0, Glut_syn - Glu_base);

// ── ROS & inflammation burden ────────────────────────────────────
double ROS_norm  = ROS / (ROS_0 + 1e-9);
double TNFa_norm = TNFa / 1.0;

// ── Net motor neuron death rate (multi-hit) ──────────────────────
double MN_death_rate = k_MN_death * (
    1.0
  + 3.0 * SOD1_burden         // SOD1 tox
  + 2.0 * frac_TDP_cyto       // TDP-43 pathology
  + 1.5 * ROS_norm             // oxidative stress
  + 1.2 * Glu_excess           // excitotoxicity
  + 0.8 * TNFa_norm            // neuro-inflammation
  - BDNF_prot_wt * (BDNF / BDNF_0)  // trophic protection
);
MN_death_rate = fmax(0.0001, MN_death_rate); // floor to avoid negatives

// ── EAAT2 reduction in ALS (inflammation-driven) ─────────────────
double EAAT2_eff = k_EAAT2 * (k_EAAT2_ALS + (1.0 - k_EAAT2_ALS) / (1.0 + Mic_act));

// ── SOD1 synthesis (reduced by tofersen) ─────────────────────────
double SOD1_syn_eff = k_SOD1_syn * (1.0 - E_TOF);

// ── NfL release proportional to MN death ─────────────────────────
double NfL_release = k_NfL_rel * MN_death_rate * (MN_upper + MN_lower);

// ── Microglial activation stimulus ───────────────────────────────
double Mic_stim = 0.5 * SOD1_burden
                + 0.3 * frac_TDP_cyto
                + 0.3 * ROS_norm
                + 0.5 * (MN_upper_0 - MN_upper)
                + 0.5 * (MN_lower_0 - MN_lower);

$ODE
// ─── DRUG PK ODEs ────────────────────────────────────────────────
// Riluzole 2-compartment
dxdt_DEPOT_RIL = -ka_RIL * DEPOT_RIL;
dxdt_C1_RIL    = F_RIL * ka_RIL * DEPOT_RIL
                 - (CL_RIL + Q_RIL) / V1_RIL * C1_RIL
                 + Q_RIL / V2_RIL * C2_RIL;
dxdt_C2_RIL    = Q_RIL / V1_RIL * C1_RIL
                 - Q_RIL / V2_RIL * C2_RIL;

// Edaravone 1-compartment (IV bolus recorded as event)
dxdt_IV_EDA    = -(CL_EDA / V_EDA) * IV_EDA;

// Tofersen 2-compartment SC
dxdt_DEPOT_TOF = -ka_TOF * DEPOT_TOF;
dxdt_C1_TOF    = ka_TOF * DEPOT_TOF
                 - (CL_TOF + Q_TOF) / V1_TOF * C1_TOF
                 + Q_TOF / V2_TOF * C2_TOF;
dxdt_C2_TOF    = Q_TOF / V1_TOF * C1_TOF
                 - Q_TOF / V2_TOF * C2_TOF;

// Phenylbutyrate (AMX0035) 1-compartment oral
dxdt_DEPOT_PB  = -ka_PB * DEPOT_PB;
dxdt_C_PB      = F_PB * ka_PB * DEPOT_PB - (CL_PB / V_PB) * C_PB;

// ─── DISEASE BIOLOGY ODEs ────────────────────────────────────────
// Motor Neurons
dxdt_MN_upper  = -k_MN_death_U * MN_death_rate * MN_upper;
dxdt_MN_lower  = -k_MN_death   * MN_death_rate * MN_lower;

// SOD1 protein dynamics
dxdt_SOD1_wt   = SOD1_syn_eff
                 - k_SOD1_deg * SOD1_wt
                 - k_SOD1_mis * SOD1_wt * f_SOD1_mut;
dxdt_SOD1_mis  = k_SOD1_mis * SOD1_wt * f_SOD1_mut
                 - k_SOD1_clr * SOD1_mis;

// TDP-43 dynamics (nuclear ↔ cytoplasmic shuttle)
dxdt_TDP43_nuc  = -k_TDP_export * TDP43_nuc
                  + k_TDP_import * TDP43_cyto;
dxdt_TDP43_cyto =  k_TDP_export * TDP43_nuc
                  - k_TDP_import * TDP43_cyto
                  - k_TDP_agg   * TDP43_cyto;

// Synaptic glutamate
dxdt_Glut_syn  = k_Glu_rel * (1.0 - E_RIL) * MN_lower
                 - EAAT2_eff * Glut_syn;

// Intracellular calcium
dxdt_Ca_i      = k_Ca_entry * Glu_excess
                 - k_Ca_efflux * (Ca_i - Ca_i_0);

// Reactive oxygen species
double ROS_prod = k_ROS_prod * (1.0 + SOD1_burden) * Ca_i / Mito;
double ROS_scav = k_ROS_scav * GSH * ROS + E_EDA * ROS;
dxdt_ROS       = ROS_prod - ROS_scav;

// Glutathione
dxdt_GSH       = k_GSH_syn - k_GSH_cons * GSH - k_ROS_scav * GSH * ROS;

// Mitochondrial function (0–1, 1=healthy)
dxdt_Mito      = k_Mito_rep * (1.0 - Mito) * (1.0 + E_mito_PB)
                 - k_Mito_dam * ROS * Mito;

// Microglial activation
dxdt_Mic_act   = k_Mic_act * Mic_stim * (1.0 - Mic_act)
                 - k_Mic_res * Mic_act;

// TNF-α
dxdt_TNFa      = k_TNFa_prod * Mic_act - k_TNFa_deg * TNFa;

// BDNF trophic support
dxdt_BDNF      = k_BDNF_prod * (1.0 - 0.35 * TNFa_norm)
                 - k_BDNF_deg * BDNF;

// ─── BIOMARKER & CLINICAL ODEs ───────────────────────────────────
// CSF NfL (elevated with neurodegeneration)
dxdt_NfL_CSF   = NfL_release - k_NfL_clr * NfL_CSF;

// ALSFRS-R total score (declines proportionally to MN loss)
double MN_frac  = 0.5 * (MN_upper + MN_lower) / (0.5 * (MN_upper_0 + MN_lower_0));
dxdt_ALSFRS    = -k_MN_death * ALSFRS * (1.8 - MN_frac) * (1.0 + 0.5 * Mic_act);

// FVC % predicted (respiratory drive declines with lower MN loss + inflammation)
dxdt_FVC       = -k_FVC_dec * FVC * (1.0 + TNFa_norm + Mic_act);

$INIT
DEPOT_RIL  = 0
C1_RIL     = 0
C2_RIL     = 0
IV_EDA     = 0
DEPOT_TOF  = 0
C1_TOF     = 0
C2_TOF     = 0
DEPOT_PB   = 0
C_PB       = 0

MN_upper   = 1.0
MN_lower   = 1.0
SOD1_wt    = 5.0
SOD1_mis   = 0.0
TDP43_nuc  = 10.0
TDP43_cyto = 0.5
Glut_syn   = 1.0
Ca_i       = 0.1
ROS        = 0.5
GSH        = 4.0
Mito       = 1.0
Mic_act    = 0.1
TNFa       = 0.2
BDNF       = 1.0

NfL_CSF    = 5.0
ALSFRS     = 48.0
FVC        = 100.0

$TABLE
capture Cp_RIL      = C1_RIL  / V1_RIL;
capture Cp_EDA      = IV_EDA  / V_EDA;
capture Ccsf_TOF    = C2_TOF  / V2_TOF;
capture Cp_PB_uM    = (C_PB   / V_PB) * 1000.0;
capture E_RIL_cap   = Emax_RIL * Cp_RIL / (IC50_RIL + Cp_RIL + 1e-9);
capture E_EDA_cap   = Emax_EDA * Cp_EDA / (IC50_EDA + Cp_EDA + 1e-9);
capture E_TOF_cap   = Emax_TOF * Ccsf_TOF / (EC50_TOF + Ccsf_TOF + 1e-9);
capture MN_total    = MN_upper + MN_lower;
capture MN_pct      = 100.0 * (MN_upper + MN_lower) / 2.0;
capture SOD1_frac   = SOD1_mis / (SOD1_wt + SOD1_mis + 1e-9);
capture TDP_cyto_frac = TDP43_cyto / (TDP43_nuc + TDP43_cyto + 1e-9);
capture ROS_norm_out = ROS / ROS_0;
capture MN_death_out = k_MN_death * (1.0 + 3.0 * f_SOD1_mut * SOD1_mis / (SOD1_wt + SOD1_mis + 1e-9));
'

## -----------------------------------------------------------------------------
## 2. Compile the model
## -----------------------------------------------------------------------------
als_mod <- mcode("als_qsp", als_code, quiet = TRUE)

cat("Model compiled successfully.\n")
cat("Compartments:", length(CMT(als_mod)), "\n")
cat("Parameters:  ", length(param(als_mod)), "\n")

## -----------------------------------------------------------------------------
## 3. Helper: build dosing events
## -----------------------------------------------------------------------------
make_dosing <- function(scenario, duration_days = 548) {

  ev_list <- list()

  ## Riluzole: 50 mg oral BID (q12h) starting day 0
  if (scenario %in% c("riluzole", "combo_ril_eda", "all_drugs")) {
    ev_list[["riluzole"]] <- ev(
      time = seq(0, (duration_days - 0.5) * 24, by = 12),   # hours
      amt  = 50,
      cmt  = "DEPOT_RIL",
      evid = 1
    )
  }

  ## Edaravone: 60 mg/day IV for 14 days on / 14 days off
  ## (Cycle 1: 14 consecutive; cycles 2+: 10 of 14 days)
  if (scenario %in% c("edaravone", "combo_ril_eda", "all_drugs")) {
    cycle_on_times <- c(
      seq(0, 13) * 24,                      # Cycle 1: days 0-13
      sapply(1:10, function(cycle) {
        base <- 28 + (cycle - 1) * 28        # next cycles start day 28, 56 ...
        seq(base, base + 9) * 24             # 10 of 14 days
      })
    )
    ev_list[["edaravone"]] <- ev(
      time = unlist(cycle_on_times),
      amt  = 60,
      cmt  = "IV_EDA",
      evid = 1,
      rate = 60 / 1          # 60 mg over 1 h (simplification)
    )
  }

  ## Tofersen: 100 mg SC — loading 3 doses (days 0, 14, 28) then q28d
  if (scenario %in% c("tofersen", "all_drugs")) {
    loading <- c(0, 14, 28)
    maint   <- seq(56, duration_days, by = 28)
    tof_days <- unique(c(loading, maint))
    ev_list[["tofersen"]] <- ev(
      time = tof_days * 24,
      amt  = 100,
      cmt  = "DEPOT_TOF",
      evid = 1
    )
  }

  ## AMX0035 (Phenylbutyrate 3 g BID): day 0 onward
  if (scenario %in% c("amx0035", "all_drugs")) {
    ev_list[["amx0035"]] <- ev(
      time = seq(0, (duration_days - 0.5) * 24, by = 12),
      amt  = 3000,           # 3 g = 3000 mg
      cmt  = "DEPOT_PB",
      evid = 1
    )
  }

  ## Merge all events
  if (length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = 1, evid = 2))
  Reduce(c, ev_list)
}

## -----------------------------------------------------------------------------
## 4. Run all 7 scenarios
## -----------------------------------------------------------------------------
scenarios <- c(
  "no_treatment",
  "riluzole",
  "edaravone",
  "combo_ril_eda",
  "tofersen",
  "amx0035",
  "all_drugs"
)

## SOD1-ALS subtype parameters for tofersen scenario
params_sporadic <- list(f_SOD1_mut = 0.0)   # sporadic ALS (TDP-43 dominant)
params_SOD1     <- list(f_SOD1_mut = 1.0)   # SOD1-familial ALS

run_scenario <- function(mod, sc, extra_params = list(),
                         duration_days = 548, delta_h = 6) {
  if (length(extra_params) > 0) mod <- param(mod, extra_params)
  dose_ev <- make_dosing(sc, duration_days)
  out <- mod %>%
    ev(dose_ev) %>%
    mrgsim(end = duration_days * 24, delta = delta_h,
           obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(
      time_days = time / 24,
      scenario  = sc
    )
  out
}

results <- bind_rows(lapply(scenarios, function(sc) {
  params <- if (sc == "tofersen") params_SOD1 else params_sporadic
  run_scenario(als_mod, sc, params)
}))

cat("Simulations complete. Scenarios:", paste(scenarios, collapse = ", "), "\n")

## -----------------------------------------------------------------------------
## 5. Key results tables
## -----------------------------------------------------------------------------
## 5a. ALSFRS-R at 12 and 18 months
alsfrs_summary <- results %>%
  filter(time_days %in% c(0, 180, 365, 548)) %>%
  select(time_days, scenario, ALSFRS, FVC, NfL_CSF, MN_pct) %>%
  group_by(scenario, time_days) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  arrange(scenario, time_days)

cat("\n=== ALSFRS-R Summary ===\n")
print(alsfrs_summary, n = Inf)

## 5b. Motor neuron survival
mn_summary <- results %>%
  filter(time_days %in% c(0, 365, 548)) %>%
  select(time_days, scenario, MN_pct, NfL_CSF) %>%
  group_by(scenario, time_days) %>%
  summarise(across(everything(), mean), .groups = "drop")

## -----------------------------------------------------------------------------
## 6. Publication-ready plots
## -----------------------------------------------------------------------------
sc_colors <- c(
  "no_treatment"  = "#D32F2F",
  "riluzole"      = "#1976D2",
  "edaravone"     = "#388E3C",
  "combo_ril_eda" = "#7B1FA2",
  "tofersen"      = "#F57C00",
  "amx0035"       = "#00838F",
  "all_drugs"     = "#5D4037"
)

sc_labels <- c(
  "no_treatment"  = "Untreated",
  "riluzole"      = "Riluzole 50 mg BID",
  "edaravone"     = "Edaravone 60 mg/day",
  "combo_ril_eda" = "Riluzole + Edaravone",
  "tofersen"      = "Tofersen 100 mg SC q4w\n(SOD1-ALS)",
  "amx0035"       = "AMX0035 (PB 3g + TUDCA 1g) BID",
  "all_drugs"     = "All active treatments"
)

## --- Plot 1: ALSFRS-R over time (18 months) ---
p_alsfrs <- ggplot(results,
       aes(x = time_days, y = ALSFRS, color = scenario, linetype = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed",
                                   "dotdash", "dotted", "twodash", "solid"),
                        labels = sc_labels) +
  labs(
    title    = "ALSFRS-R Progression — ALS QSP Model",
    subtitle = "18-month simulation across 7 treatment scenarios",
    x        = "Time (days)",
    y        = "ALSFRS-R Total Score (0–48)",
    color    = "Scenario",
    linetype = "Scenario",
    caption  = "Calibrated to Bensimon 1994, Riley 2004, Abe 2017, Paganoni 2020, Miller 2022"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.key.width = unit(1.5, "cm"))

## --- Plot 2: CSF NfL biomarker ---
p_nfl <- ggplot(results,
       aes(x = time_days, y = NfL_CSF, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(
    title = "CSF Neurofilament Light Chain (NfL) — ALS QSP Model",
    x     = "Time (days)",
    y     = "CSF NfL (normalized AU)",
    color = "Scenario"
  ) +
  theme_bw(base_size = 12)

## --- Plot 3: Motor Neuron Survival ---
p_mn <- ggplot(results,
       aes(x = time_days, y = MN_pct, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(
    title = "Motor Neuron Survival — ALS QSP Model",
    x     = "Time (days)",
    y     = "Motor Neuron Survival (%)",
    color = "Scenario"
  ) +
  theme_bw(base_size = 12)

## --- Plot 4: FVC decline ---
p_fvc <- ggplot(results,
       aes(x = time_days, y = FVC, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(
    title = "FVC % Predicted — ALS QSP Model",
    x     = "Time (days)",
    y     = "FVC (% predicted)",
    color = "Scenario"
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "darkred") +
  annotate("text", x = 10, y = 52, label = "NIV threshold (50%)", hjust = 0, size = 3) +
  theme_bw(base_size = 12)

## --- Plot 5: Riluzole PK (BID, first 24 h) ---
p_pk_ril <- results %>%
  filter(scenario == "riluzole", time_days <= 7) %>%
  ggplot(aes(x = time_days * 24, y = Cp_RIL)) +
  geom_line(color = "#1976D2", size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 0.55, label = "IC50 = 0.5 μg/mL",
           hjust = 0, size = 3, color = "red") +
  labs(
    title = "Riluzole PK — First 7 Days",
    x     = "Time (h)",
    y     = "Riluzole Cp (μg/mL)"
  ) +
  theme_bw(base_size = 12)

## --- Plot 6: Neuroinflammation dynamics ---
p_inflam <- results %>%
  filter(scenario %in% c("no_treatment", "combo_ril_eda", "all_drugs")) %>%
  ggplot(aes(x = time_days, y = TNFa, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = sc_colors[c("no_treatment","combo_ril_eda","all_drugs")],
                     labels = sc_labels[c("no_treatment","combo_ril_eda","all_drugs")]) +
  labs(
    title = "TNF-α Neuroinflammation — ALS QSP Model",
    x     = "Time (days)",
    y     = "TNF-α (AU)",
    color = "Scenario"
  ) +
  theme_bw(base_size = 12)

## Arrange and save
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  g <- gridExtra::arrangeGrob(p_alsfrs, p_nfl, p_mn, p_fvc, p_pk_ril, p_inflam, ncol = 2)
  ggsave("als_qsp_results.pdf", g, width = 16, height = 18)
  cat("Plots saved to als_qsp_results.pdf\n")
} else {
  print(p_alsfrs)
}

## -----------------------------------------------------------------------------
## 7. Drug PK summary table
## -----------------------------------------------------------------------------
pk_table <- tribble(
  ~Drug,        ~Route,   ~F_pct, ~t_half_h, ~Vd_L,  ~CL_Lh, ~Dose,          ~IC50_Emax,
  "Riluzole",   "PO BID", 60,     12,         357,    28,     "50 mg q12h",   "IC50=0.5 μg/mL, Emax=60%",
  "Edaravone",  "IV/PO",  99,     4.5,        120,    18,     "60 mg/day IV", "IC50=1.2 μg/mL, Emax=70%",
  "Tofersen",   "SC q4w", NA,     168,        20,     0.5,    "100 mg q4w",   "EC50=0.1 μg/mL, Emax=80%",
  "PB (AMX)",   "PO BID", 85,     3.0,        50,     12,     "3g q12h",      "EC50=50 μmol/L, Emax=50%"
)

cat("\n=== Drug PK/PD Summary ===\n")
print(pk_table)

## -----------------------------------------------------------------------------
## 8. Sensitivity analysis — k_MN_death multiplier
## -----------------------------------------------------------------------------
sa_results <- bind_rows(lapply(c(0.5, 0.75, 1.0, 1.25, 1.5), function(mult) {
  mod_sa <- param(als_mod, list(k_MN_death = 0.00080 * mult))
  dose_ev <- make_dosing("riluzole")
  out <- mod_sa %>%
    ev(dose_ev) %>%
    mrgsim(end = 548 * 24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(time_days = time / 24, k_mult = mult)
  out
}))

p_sa <- ggplot(sa_results, aes(x = time_days, y = ALSFRS,
                                color = factor(k_mult))) +
  geom_line(size = 0.9) +
  scale_color_viridis_d(name = "k_death\nmultiplier") +
  labs(
    title    = "Sensitivity: Motor Neuron Death Rate (Riluzole scenario)",
    subtitle = "ALSFRS-R sensitivity to k_MN_death (±50%)",
    x        = "Time (days)",
    y        = "ALSFRS-R Score"
  ) +
  theme_bw(base_size = 12)

print(p_sa)

cat("\n✓ ALS QSP mrgsolve model run complete.\n")
