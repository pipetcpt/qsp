## ============================================================
## Primary Open-Angle Glaucoma (POAG)
## Quantitative Systems Pharmacology (QSP) Model
## mrgsolve ODE implementation
##
## Pathways modelled:
##   1. Topical drug PK (5 drug classes, first-order aqueous compartment)
##   2. Aqueous humor dynamics → IOP (Goldmann equation, dynamic)
##   3. Trabecular meshwork (TM) biology (ECM accumulation, C_tm decline)
##   4. Optic nerve head (ONH) BDNF dynamics
##   5. RGC apoptosis (Caspase-3 → RGC loss)
##   6. Structural & functional endpoints (RNFL, VF-MD)
##
## Treatment scenarios:
##   1.  Untreated POAG (natural history)
##   2.  Latanoprost 0.005% QD (prostaglandin analog, PGA)
##   3.  Timolol 0.5% BID (beta-blocker, BB)
##   4.  Dorzolamide 2% TID (carbonic anhydrase inhibitor, CAI)
##   5.  Brimonidine 0.2% BID (alpha-2 agonist, A2A + neuroprotection)
##   6.  Netarsudil 0.02% QD (ROCK inhibitor)
##   7.  Fixed combo: Latanoprost + Timolol BID (DuoTrav)
##   8.  Triple therapy: PGA + BB + CAI
##
## Key references:
##   Goldmann (1951) Ann Ophthalmol; AGIS (2000) Am J Ophthalmol;
##   OHTS (2002) JAMA; CNTGS; Leske et al. (2003) Arch Ophthalmol;
##   Quigley & Broman (2006) Br J Ophthalmol;
##   Crowston & Weinreb (2008) Exp Eye Res; Almasieh et al. (2012) PRSB
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---- Model definition ------------------------------------------------
code <- '
$PROB POAG QSP Model - mrgsolve

$PARAM
// ---- Drug PK parameters ----
// Prostaglandin analog (PGA, e.g. Latanoprost acid)
ke_PGA   = 0.693    // Elimination from aqueous (1/h); t1/2 ~ 1 h
F_abs_PGA = 0.10    // Bioavailability into aqueous (fraction of dose)
Dose_PGA  = 0.0     // ug dose per administration (50 ug Latanoprost = activated acid ~ 2 ug)

// Beta-blocker (BB, e.g. Timolol)
ke_BB    = 0.231    // Elimination (1/h); t1/2 ~ 3 h in aqueous
F_abs_BB = 0.025    // Bioavailability into aqueous
Dose_BB  = 0.0      // ug

// Carbonic anhydrase inhibitor (CAI, e.g. Dorzolamide)
ke_CAI   = 0.139    // Elimination (1/h); t1/2 ~ 5 h in aqueous
F_abs_CAI = 0.04    // Bioavailability into aqueous
Dose_CAI  = 0.0     // ug

// Alpha-2 agonist (A2A, e.g. Brimonidine)
ke_A2A   = 0.462    // Elimination (1/h); t1/2 ~ 1.5 h
F_abs_A2A = 0.03    // Bioavailability into aqueous
Dose_A2A  = 0.0     // ug

// ROCK inhibitor (ROCKi, e.g. Netarsudil)
ke_ROCK  = 0.347    // Elimination (1/h); t1/2 ~ 2 h
F_abs_ROCK = 0.05   // Bioavailability into aqueous
Dose_ROCK = 0.0     // ug

// ---- Aqueous humor dynamics ----
F_prod_base  = 2.5  // Baseline aqueous production (uL/min)
F_uv_base    = 0.40 // Baseline uveoscleral outflow (uL/min)
C_tm_base    = 0.30 // Baseline TM outflow facility (uL/min/mmHg)
P_ep         = 8.0  // Episcleral venous pressure (mmHg)
V_aq         = 250  // Aqueous volume (uL)
// IOP dynamics time constant
tau_IOP      = 0.1  // Hours to equilibrate IOP (fast ~ 6 min)

// ---- Drug PD: Maximum effects (Emax) ----
// PGA → Increase F_uv
Emax_PGA  = 1.0     // Emax fractional increase in F_uv (100%)
EC50_PGA  = 2.0     // EC50 in ng/mL
hill_PGA  = 1.5

// BB → Decrease F_prod
Emax_BB   = 0.30    // Emax fractional decrease (30%)
EC50_BB   = 150.0   // EC50 in ng/mL
hill_BB   = 1.2

// CAI → Decrease F_prod
Emax_CAI  = 0.25    // Emax fractional decrease (25%)
EC50_CAI  = 500.0   // EC50 in ng/mL
hill_CAI  = 1.0

// A2A → Decrease F_prod + increase F_uv
Emax_A2A_prod = 0.25 // Emax ↓ F_prod (25%)
EC50_A2A      = 1.5  // EC50 in ng/mL
hill_A2A      = 1.2
Emax_A2A_uv   = 0.10 // Emax ↑ F_uv (10%)

// ROCK inhibitor → Increase C_tm
Emax_ROCK = 0.35    // Emax fractional increase in C_tm (35%)
EC50_ROCK = 30.0    // EC50 in ng/mL
hill_ROCK = 1.5
Emax_ROCK_ep = 3.0  // Emax decrease in P_ep (mmHg)

// ---- TM biology ----
k_ECM_prod  = 0.001  // ECM accumulation rate (1/h per unit TM stress)
k_ECM_clear = 0.0005 // Baseline ECM clearance (1/h)
k_TM_IOP    = 0.01   // IOP-driven ECM accumulation rate
k_C_tm_decline = 2.0e-4 // ECM → C_tm decline (1/h per ECM unit)
C_tm_min    = 0.10   // Minimum C_tm (severe TM damage)

// ---- BDNF dynamics ----
BDNF0       = 50.0   // Baseline BDNF (pg/mL in ONH)
k_BDNF_prod = 0.05   // BDNF production rate (pg/mL/h)
k_BDNF_deg  = 0.001  // BDNF degradation (1/h)
k_BDNF_IOP  = 0.002  // IOP-dependent BDNF depletion (/mmHg/h above threshold)
IOP_thresh_BDNF = 18.0 // IOP threshold above which BDNF falls
k_A2A_BDNF  = 0.20   // A2A agonism increases BDNF synthesis (fraction)

// ---- RGC apoptosis ----
RGC0        = 1.2    // Baseline RGC count (millions)
k_RGC_base  = 5.7e-6 // Normal age-related loss (1/h; ~0.5%/yr)
EC50_Casp3  = 0.3    // Casp3 EC50 for RGC loss
hill_Casp3  = 2.0
k_Casp3_act = 0.005  // Caspase-3 activation rate (proportional to IOP stress)
k_Casp3_inact = 0.02 // Caspase-3 inactivation (1/h)
IOP_Casp3_thresh = 15.0 // IOP threshold for Casp3 activation (mmHg)
hill_IOP_Casp3 = 2.5

// ---- Structural/Functional endpoints ----
RNFL0       = 100.0  // Baseline RNFL (um)
k_RNFL_RGC  = 1.0    // RNFL ∝ RGC (um per million RGC)
VF_MD_floor = -30.0  // Minimum VF-MD (dB; total blindness)
k_VF_RGC    = 25.0   // Slope of VF-MD vs RGC fraction (dB)
n_VF        = 2.5    // Hill coefficient for VF-MD vs RGC

// ---- Diurnal IOP variation ----
diurnal_amp  = 0.10  // Amplitude as fraction of baseline IOP (10%)

$CMT
// Drug concentrations in aqueous (ng/mL)
C_PGA   // Prostaglandin analog (ng/mL)
C_BB    // Beta-blocker (ng/mL)
C_CAI   // CAI (ng/mL)
C_A2A   // Alpha-2 agonist (ng/mL)
C_ROCK  // ROCK inhibitor (ng/mL)

// Aqueous humor dynamics
F_aq    // Aqueous production rate (uL/min) - dynamic state
F_uv    // Uveoscleral outflow rate (uL/min) - dynamic state
C_tm    // Trabecular outflow facility (uL/min/mmHg) - dynamic state
IOP     // Intraocular pressure (mmHg)

// TM biology
ECM_TM  // ECM accumulation index (0 = normal, dimensionless)

// ONH/neuroprotection
BDNF    // BDNF in ONH (pg/mL)
Casp3   // Caspase-3 activity (0-1 index)

// Cell count and structural
RGC     // Retinal ganglion cells (millions)
RNFL    // RNFL thickness (um)
VF_MD   // VF mean deviation (dB)

$MAIN
// Initial conditions
F_aq_0  = F_prod_base;
F_uv_0  = F_uv_base;
C_tm_0  = C_tm_base;
IOP_0   = (F_prod_base - F_uv_base) / C_tm_base + P_ep;  // ~14.7 mmHg

ECM_TM_0 = 0.0;
BDNF_0   = BDNF0;
Casp3_0  = 0.01;  // Low baseline
RGC_0    = RGC0;
RNFL_0   = RNFL0;
VF_MD_0  = 0.0;   // Normal (dB)

$ODE
// ============================================================
// SECTION 1: Drug PK in aqueous humor
// ============================================================
// Drug appears in aqueous via zero-order input events (bolus → 1-cpt)
dxdt_C_PGA  = -ke_PGA  * C_PGA;
dxdt_C_BB   = -ke_BB   * C_BB;
dxdt_C_CAI  = -ke_CAI  * C_CAI;
dxdt_C_A2A  = -ke_A2A  * C_A2A;
dxdt_C_ROCK = -ke_ROCK * C_ROCK;

// ============================================================
// SECTION 2: Drug PD — Emax models
// ============================================================
// Hill equations for each drug effect
double E_PGA  = Emax_PGA  * pow(C_PGA,  hill_PGA)  / (pow(EC50_PGA,  hill_PGA)  + pow(C_PGA,  hill_PGA));
double E_BB   = Emax_BB   * pow(C_BB,   hill_BB)   / (pow(EC50_BB,   hill_BB)   + pow(C_BB,   hill_BB));
double E_CAI  = Emax_CAI  * pow(C_CAI,  hill_CAI)  / (pow(EC50_CAI,  hill_CAI)  + pow(C_CAI,  hill_CAI));
double E_A2A  = Emax_A2A_prod * pow(C_A2A, hill_A2A) / (pow(EC50_A2A, hill_A2A) + pow(C_A2A, hill_A2A));
double E_A2A_uv = Emax_A2A_uv * pow(C_A2A, hill_A2A) / (pow(EC50_A2A, hill_A2A) + pow(C_A2A, hill_A2A));
double E_ROCK = Emax_ROCK * pow(C_ROCK, hill_ROCK) / (pow(EC50_ROCK, hill_ROCK) + pow(C_ROCK, hill_ROCK));
double E_ROCK_ep = Emax_ROCK_ep * pow(C_ROCK, hill_ROCK) / (pow(EC50_ROCK, hill_ROCK) + pow(C_ROCK, hill_ROCK));

// Diurnal variation (cosine, peak at 09:00 = t mod 24 == 9)
double t_mod24 = fmod(SOLVERTIME / 1.0, 24.0);  // t is in hours
double diurnal_factor = 1.0 + diurnal_amp * cos(2 * 3.14159 * (t_mod24 - 9.0) / 24.0);

// ============================================================
// SECTION 3: Aqueous humor dynamics
// ============================================================
// Target values modified by drugs
double F_aq_target = F_prod_base * diurnal_factor * (1.0 - E_BB - E_CAI - E_A2A);
double F_uv_target = F_uv_base * (1.0 + E_PGA + E_A2A_uv);
// C_tm: disease reduces it via ECM, drugs can increase
double C_tm_disease = C_tm_base * exp(-k_C_tm_decline * ECM_TM);
if(C_tm_disease < C_tm_min) C_tm_disease = C_tm_min;
double C_tm_target = C_tm_disease * (1.0 + E_ROCK);

// P_ep modified by ROCK inhibitor
double P_ep_eff = P_ep - E_ROCK_ep;

// Differential equations for F_aq, F_uv, C_tm (slow adaptation)
dxdt_F_aq  = 5.0 * (F_aq_target  - F_aq);   // fast adjustment (h scale)
dxdt_F_uv  = 2.0 * (F_uv_target  - F_uv);
dxdt_C_tm  = 0.5 * (C_tm_target  - C_tm);

// Goldmann equation: equilibrium IOP
double IOP_eq = (F_aq - F_uv) / C_tm + P_ep_eff;
if(IOP_eq < 5.0) IOP_eq = 5.0;   // physiological minimum

dxdt_IOP = (1.0 / tau_IOP) * (IOP_eq - IOP);

// ============================================================
// SECTION 4: TM biology — ECM accumulation
// ============================================================
// ECM accumulates due to IOP-driven stress and aging
double IOP_stress = (IOP > 18.0) ? (IOP - 18.0) : 0.0;
double ECM_prod_rate = k_ECM_prod + k_TM_IOP * IOP_stress;
double ECM_clear_rate = k_ECM_clear;  // ROCKi can enhance; simplified here
double ECM_ROCK_clear = 0.5 * E_ROCK * ECM_TM;  // ROCK-I reduces ECM

dxdt_ECM_TM = ECM_prod_rate - ECM_clear_rate * ECM_TM - ECM_ROCK_clear;
if(ECM_TM < 0.0) ECM_TM = 0.0;

// ============================================================
// SECTION 5: BDNF dynamics in ONH
// ============================================================
// BDNF depleted by IOP (>threshold) and replenished by A2A neuroprotection
double IOP_above = (IOP > IOP_thresh_BDNF) ? (IOP - IOP_thresh_BDNF) : 0.0;
double A2A_BDNF_effect = 1.0 + k_A2A_BDNF * E_A2A;
double BDNF_prod = k_BDNF_prod * A2A_BDNF_effect;
double BDNF_degrad = k_BDNF_deg + k_BDNF_IOP * IOP_above;

dxdt_BDNF = BDNF_prod - BDNF_degrad * BDNF;

// ============================================================
// SECTION 6: Caspase-3 activation → RGC apoptosis
// ============================================================
// Caspase-3 driven by IOP stress and BDNF deprivation
double IOP_Casp_stress = pow(IOP / IOP_Casp3_thresh, hill_IOP_Casp3);
double BDNF_survival = BDNF / (EC50_Casp3 * 100 + BDNF);   // BDNF protective
double Casp3_act_rate = k_Casp3_act * IOP_Casp_stress * (1.0 - BDNF_survival);

dxdt_Casp3 = Casp3_act_rate - k_Casp3_inact * Casp3;
if(Casp3 < 0.0) Casp3 = 0.0;
if(Casp3 > 1.0) Casp3 = 1.0;

// ============================================================
// SECTION 7: RGC loss
// ============================================================
double Casp3_eff = pow(Casp3, hill_Casp3) / (pow(EC50_Casp3, hill_Casp3) + pow(Casp3, hill_Casp3));
double k_RGC_loss = k_RGC_base + k_Casp3_act * 2.0 * Casp3_eff;

dxdt_RGC = -k_RGC_loss * RGC;
if(RGC < 0.0) RGC = 0.0;

// ============================================================
// SECTION 8: Structural/Functional endpoints
// ============================================================
// RNFL ∝ RGC
dxdt_RNFL = k_RNFL_RGC * (RGC - RNFL / k_RNFL_RGC) * 0.001;  // slow structural loss

// VF-MD: nonlinear relationship to RGC (dB)
double RGC_frac = RGC / RGC0;
if(RGC_frac < 0.01) RGC_frac = 0.01;
double VF_MD_eq = -k_VF_RGC * pow(1.0 - RGC_frac, n_VF);
if(VF_MD_eq < VF_MD_floor) VF_MD_eq = VF_MD_floor;
dxdt_VF_MD = 0.05 * (VF_MD_eq - VF_MD);

$CAPTURE
// IOP and aqueous dynamics
IOP F_aq F_uv C_tm P_ep_eff IOP_eq
// Drug concentrations
C_PGA C_BB C_CAI C_A2A C_ROCK
// Drug effects
E_PGA E_BB E_CAI E_A2A E_ROCK
// TM & BDNF
ECM_TM BDNF Casp3
// Endpoints
RGC RNFL VF_MD
// Derived
RGC_frac Casp3_eff
'

# Compile model
mod <- mread_cache("poag_qsp", tempdir(), code)

# ============================================================
# SIMULATION SCENARIOS
# ============================================================

# Daily dosing events: e.g., Latanoprost QD at 21:00 h (dose = 50 ug × 0.10 = 5 ng/mL effective)
# For simplicity we use bolus events that directly add to C_PGA etc.
# Scaling: 50 ug × F_abs × 1000/250uL ~ effective aqueous Cmax in ng/mL
# Latanoprost acid Cmax in aqueous ~ 15-25 nM = 5-9 ng/mL → use 8 ng/mL per dose
# Timolol Cmax ~ 800 nM = ~220 ng/mL → use 200 ng/mL per dose (BID)
# Dorzolamide Cmax ~ 5 uM = 1600 ng/mL → use 1500 ng/mL
# Brimonidine Cmax ~ 1.4 nM = 0.38 ng/mL → use 0.5 ng/mL
# Netarsudil Cmax ~ 0.1 uM = 53 ng/mL → use 50 ng/mL

sim_years <- 10   # Simulate 10 years
sim_hours <- sim_years * 365 * 24

# Helper: create dosing events for each scenario
make_events <- function(scenario = "untreated") {
  # Base time vector (h)
  end_h <- sim_hours

  ev_list <- list()

  if(scenario %in% c("PGA", "PGA_BB", "triple")) {
    # Latanoprost QD at 21h (evening)
    n_doses <- floor(end_h / 24)
    ev_PGA <- ev(cmt = "C_PGA", amt = 8, time = seq(21, n_doses * 24 + 21, by = 24),
                 rate = -2)  # instantaneous bolus
    ev_list[[length(ev_list)+1]] <- ev_PGA
  }
  if(scenario %in% c("BB", "PGA_BB", "triple")) {
    # Timolol BID at 8h and 20h
    n_doses <- floor(end_h / 24)
    times_BB <- c(sapply(0:(n_doses-1), function(d) c(d*24+8, d*24+20)))
    ev_BB <- ev(cmt = "C_BB", amt = 200, time = times_BB, rate = -2)
    ev_list[[length(ev_list)+1]] <- ev_BB
  }
  if(scenario %in% c("CAI", "triple")) {
    # Dorzolamide TID at 8, 14, 20h
    n_doses <- floor(end_h / 24)
    times_CAI <- c(sapply(0:(n_doses-1), function(d) c(d*24+8, d*24+14, d*24+20)))
    ev_CAI <- ev(cmt = "C_CAI", amt = 1500, time = times_CAI, rate = -2)
    ev_list[[length(ev_list)+1]] <- ev_CAI
  }
  if(scenario == "A2A") {
    # Brimonidine BID
    n_doses <- floor(end_h / 24)
    times_A2A <- c(sapply(0:(n_doses-1), function(d) c(d*24+8, d*24+20)))
    ev_A2A <- ev(cmt = "C_A2A", amt = 0.5, time = times_A2A, rate = -2)
    ev_list[[length(ev_list)+1]] <- ev_A2A
  }
  if(scenario == "ROCK") {
    # Netarsudil QD at 21h
    n_doses <- floor(end_h / 24)
    ev_ROCK <- ev(cmt = "C_ROCK", amt = 50, time = seq(21, n_doses*24+21, by=24), rate=-2)
    ev_list[[length(ev_list)+1]] <- ev_ROCK
  }

  if(length(ev_list) == 0) return(NULL)
  ev_combined <- Reduce(c, ev_list)
  return(ev_combined)
}

# Initial condition parameters for POAG patient
# Baseline elevated IOP = 24 mmHg → adjust initial states accordingly
init_poag <- c(
  C_PGA = 0, C_BB = 0, C_CAI = 0, C_A2A = 0, C_ROCK = 0,
  F_aq  = 2.5, F_uv  = 0.40, C_tm  = 0.25,  # Diseased C_tm (normal 0.30 → POAG 0.25)
  IOP   = 24.0,   # Elevated baseline
  ECM_TM = 0.15,  # Some existing ECM in POAG
  BDNF  = 35.0,   # Mildly reduced BDNF
  Casp3 = 0.05,   # Low-level activation
  RGC   = 1.1,    # Mild early RGC loss
  RNFL  = 88.0,   # Thinned RNFL (early-moderate POAG)
  VF_MD = -3.5    # Mild VF defect (early)
)

# Observation times: sample every 24 h (daily average), but run for 10 years
out_times <- seq(0, sim_hours, by = 24)

# Run all scenarios
scenarios <- c("untreated", "PGA", "BB", "CAI", "A2A", "ROCK", "PGA_BB", "triple")
scenario_labels <- c(
  "Untreated (Natural History)",
  "Latanoprost QD (PGA)",
  "Timolol BID (Beta-Blocker)",
  "Dorzolamide TID (CAI)",
  "Brimonidine BID (A2A + Neuroprot.)",
  "Netarsudil QD (ROCK-I)",
  "Latanoprost + Timolol BID (DuoTrav)",
  "Triple: PGA + BB + CAI"
)

run_scenario <- function(scen_name) {
  events <- make_events(scen_name)
  if(is.null(events)) {
    out <- mod %>%
      init(init_poag) %>%
      mrgsim(end = sim_hours, delta = 24)
  } else {
    out <- mod %>%
      init(init_poag) %>%
      mrgsim(events = events, end = sim_hours, delta = 24)
  }
  out_df <- as.data.frame(out)
  out_df$scenario <- scen_name
  out_df$time_yr  <- out_df$time / (365 * 24)
  return(out_df)
}

results_list <- lapply(scenarios, run_scenario)
results <- do.call(rbind, results_list)

# Map scenario names to labels
scen_map <- setNames(scenario_labels, scenarios)
results$scenario_label <- scen_map[results$scenario]
results$scenario_label <- factor(results$scenario_label, levels = scenario_labels)

# ============================================================
# KEY PLOTS
# ============================================================

# Color palette for 8 scenarios
cols8 <- c("#D32F2F", "#2196F3", "#4CAF50", "#FF9800", "#9C27B0",
           "#00BCD4", "#795548", "#607D8B")

# --- Plot 1: IOP over time
p1 <- ggplot(results, aes(x = time_yr, y = IOP, color = scenario_label)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 21, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 18, linetype = "dotted", color = "blue", alpha = 0.5) +
  scale_color_manual(values = cols8) +
  labs(title = "IOP Dynamics Over 10 Years",
       x = "Time (Years)", y = "IOP (mmHg)",
       color = "Treatment Scenario",
       caption = "Red dashed: OHT threshold (21 mmHg); Blue dotted: target IOP (18 mmHg)") +
  theme_bw(base_size = 11) +
  coord_cartesian(ylim = c(8, 30))

# --- Plot 2: VF-MD progression
p2 <- ggplot(results, aes(x = time_yr, y = VF_MD, color = scenario_label)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = -6,  linetype = "dashed", color = "#FF8C00", alpha=0.6) +
  geom_hline(yintercept = -12, linetype = "dashed", color = "#D32F2F", alpha=0.6) +
  scale_color_manual(values = cols8) +
  labs(title = "Visual Field Progression (VF-MD)",
       x = "Time (Years)", y = "VF Mean Deviation (dB)",
       color = "Treatment Scenario",
       caption = "Dashed: Stage boundaries (−6 moderate, −12 advanced)") +
  theme_bw(base_size = 11)

# --- Plot 3: RNFL thinning
p3 <- ggplot(results, aes(x = time_yr, y = RNFL, color = scenario_label)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = cols8) +
  labs(title = "RNFL Thickness Over Time",
       x = "Time (Years)", y = "RNFL Thickness (μm)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11)

# --- Plot 4: RGC count
p4 <- ggplot(results, aes(x = time_yr, y = RGC, color = scenario_label)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = cols8) +
  labs(title = "Retinal Ganglion Cell (RGC) Count",
       x = "Time (Years)", y = "RGC Count (millions)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11)

# --- Plot 5: BDNF in ONH
p5 <- ggplot(results, aes(x = time_yr, y = BDNF, color = scenario_label)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = cols8) +
  labs(title = "BDNF in Optic Nerve Head",
       x = "Time (Years)", y = "BDNF (pg/mL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11)

# --- Plot 6: Caspase-3 activity
p6 <- ggplot(results, aes(x = time_yr, y = Casp3, color = scenario_label)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = cols8) +
  labs(title = "Caspase-3 Apoptotic Activity",
       x = "Time (Years)", y = "Caspase-3 Index (0-1)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11)

# --- Summary table at Year 10
summary_y10 <- results %>%
  filter(time_yr >= 9.9) %>%
  group_by(scenario_label) %>%
  summarise(
    IOP_mmHg    = round(mean(IOP, na.rm = TRUE), 1),
    VF_MD_dB    = round(mean(VF_MD, na.rm = TRUE), 2),
    RNFL_um     = round(mean(RNFL, na.rm = TRUE), 1),
    RGC_M       = round(mean(RGC,  na.rm = TRUE), 3),
    BDNF_pg     = round(mean(BDNF, na.rm = TRUE), 1),
    Casp3_idx   = round(mean(Casp3, na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\n=== POAG QSP Model — 10-Year Summary ===\n")
print(summary_y10, n = Inf)

# --- IOP reduction summary
iop_reduction <- results %>%
  filter(time_yr >= 0.49 & time_yr <= 0.51) %>%   # 6-month IOP
  group_by(scenario_label) %>%
  summarise(IOP_6mo = round(mean(IOP, na.rm = TRUE), 1), .groups = "drop") %>%
  mutate(IOP_baseline = 24.0,
         IOP_pct_reduction = round((IOP_baseline - IOP_6mo) / IOP_baseline * 100, 1))

cat("\n=== IOP Reduction at 6 Months ===\n")
print(iop_reduction)

# ============================================================
# DOSE-RESPONSE SENSITIVITY: IOP vs PGA dose
# ============================================================
dose_levels <- c(0, 1, 2, 4, 8, 12, 16) # ng/mL effective Cmax
dose_iop <- sapply(dose_levels, function(d) {
  params_d <- param(mod, Dose_PGA = d)
  # Steady-state IOP with constant PGA
  # Approximate analytically: at Css, C_PGA = dose; E_PGA = Emax*d^h/(EC50^h+d^h)
  E <- 1.0 * d^1.5 / (2.0^1.5 + d^1.5)
  F_uv_d <- 0.40 * (1 + E)
  IOP_d <- (2.5 - F_uv_d) / 0.25 + 8.0
  return(IOP_d)
})

dr_df <- data.frame(PGA_Cmax = dose_levels, IOP_ss = dose_iop)
p_dose <- ggplot(dr_df, aes(x = PGA_Cmax, y = IOP_ss)) +
  geom_line(color = "#2196F3", size = 1.2) +
  geom_point(size = 3, color = "#2196F3") +
  labs(title = "Dose-Response: PGA Concentration vs Steady-State IOP",
       x = "PGA Aqueous Concentration (ng/mL)",
       y = "Steady-State IOP (mmHg)") +
  theme_bw(base_size = 11)

# Print all plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p_dose)

message("POAG QSP model simulation complete.")
