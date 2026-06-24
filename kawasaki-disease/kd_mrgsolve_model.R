################################################################################
# Kawasaki Disease (KD) — QSP Model (mrgsolve)
#
# Mechanistic coverage:
#   • IVIG PK (2-compartment, FcRn recycling)
#   • Aspirin PK (1-compartment + salicylate metabolite)
#   • Methylprednisolone PK (2-compartment)
#   • Infliximab PK (2-compartment)
#   • Anakinra PK (1-compartment)
#   • Cytokine dynamics: IL-1β, IL-6, TNF-α
#   • Macrophage activation
#   • Endothelial activation
#   • Fever (body temperature)
#   • Platelet dynamics (thrombocytosis)
#   • C-Reactive Protein (CRP)
#   • Coronary Artery Z-score (CAL)
#   • IVIG response probability
#
# Parameter calibration references:
#   Burns 2020 (IVIG PD), Uehara 2019 (cytokine),
#   Son 2011 (IVIG resistance), Kobayashi 2006 (risk score)
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------
# Model code block
# ------------------------------------------------------------------
kd_code <- '
$PROB Kawasaki Disease QSP Model — PK/PD with 5 Treatment Scenarios

$PARAM
// ---- IVIG PK parameters ----
CL_IVIG    = 0.0033  // L/h/kg — IgG clearance
Vc_IVIG    = 0.050   // L/kg  — central volume
Vp_IVIG    = 0.048   // L/kg  — peripheral volume
Qp_IVIG    = 0.010   // L/h/kg — intercompartmental CL
ka_FcRn    = 0.015   // 1/h   — FcRn recycling rate
F_FcRn     = 0.60    // fraction recycled by FcRn

// ---- Aspirin PK parameters ----
ka_ASA     = 0.80    // 1/h   — aspirin absorption
Vc_ASA     = 0.18    // L/kg  — aspirin distribution
CL_ASA     = 0.60    // L/h/kg — aspirin clearance
f_SA       = 0.85    // fraction metabolised to salicylate
CL_SA      = 0.010   // L/h/kg — salicylate clearance (saturable at hi dose)
Vc_SA      = 0.16    // L/kg  — salicylate distribution

// ---- Methylprednisolone PK ----
CL_MP      = 0.48    // L/h/kg
Vc_MP      = 0.70    // L/kg
Vp_MP      = 0.55    // L/kg
Qp_MP      = 0.060   // L/h/kg

// ---- Infliximab PK ----
CL_IFX     = 0.0062  // L/h/kg
Vc_IFX     = 0.052   // L/kg
Vp_IFX     = 0.046   // L/kg
Qp_IFX     = 0.0038  // L/h/kg

// ---- Anakinra PK ----
CL_ANK     = 0.28    // L/h/kg  — renal clearance dominant
Vc_ANK     = 0.17    // L/kg
ka_ANK     = 0.30    // 1/h     — SC absorption

// ---- Disease dynamics (calibrated to Uehara 2019, Burns 2020) ----
// IL-1beta
kprod_IL1  = 0.12    // ng/mL/h — baseline production
kdeg_IL1   = 0.15    // 1/h     — natural degradation
ksyn_IL1   = 0.25    // amplification by macrophage activation
IL1_base   = 0.8     // ng/mL   — baseline IL-1β

// IL-6
kprod_IL6  = 0.28    // ng/mL/h
kdeg_IL6   = 0.10    // 1/h
ksyn_IL6   = 0.40    // from IL-1 and macrophage
IL6_base   = 2.5     // ng/mL

// TNF-alpha
kprod_TNF  = 0.15    // ng/mL/h
kdeg_TNF   = 0.20    // 1/h
ksyn_TNF   = 0.30
TNF_base   = 1.2     // ng/mL

// Macrophage activation
Mac0       = 1.0     // baseline (normalized)
kact_Mac   = 0.08    // 1/h — activation rate by trigger
kinact_Mac = 0.05    // 1/h — spontaneous deactivation
hill_n     = 2.0     // Hill coefficient

// Endothelial activation
kact_EC    = 0.10    // 1/h
kinact_EC  = 0.04    // 1/h
EC0        = 1.0     // baseline

// Fever dynamics
T_base     = 36.8    // °C baseline
krise_T    = 0.35    // °C/h — pyrogenic rise rate
kfall_T    = 0.10    // °C/h — spontaneous fall
T_max      = 40.5    // °C maximum

// CRP dynamics
kprod_CRP  = 3.5     // mg/L/h
kdeg_CRP   = 0.040   // 1/h
CRP_base   = 1.0     // mg/L baseline

// Platelet dynamics (thrombocytosis peak at ~2 weeks)
PLT_base   = 250     // ×10³/μL baseline
kprod_PLT  = 0.12    // ×10³/μL/h
kdeg_PLT   = 0.0060  // 1/h

// Coronary Z-score
Z0         = 0.0     // initial z-score
krise_Z    = 0.0015  // 1/h — z-score rise driven by inflammation
kfall_Z    = 0.0008  // 1/h — remodeling
Z_max      = 12.0    // maximum z-score (giant aneurysm upper bound)

// ---- PD: drug effect parameters ----
// IVIG — immunomodulation
EC50_IVIG  = 8.0     // g/L — EC50 for IVIG anti-inflammatory
Emax_IVIG  = 0.80    // max fractional inhibition of cytokines
n_IVIG     = 1.5     // Hill exponent

// Aspirin — COX inhibition (anti-pyretic)
EC50_ASA   = 15.0    // mg/L salicylate for 50% fever reduction
Emax_ASA   = 0.60    // max antipyretic effect
// Aspirin anti-platelet (COX-1 irreversible, low-dose)
IC50_COX1  = 3.0     // mg/L ASA for 50% COX-1 inhibition
Emax_COX1  = 0.95    // near-complete COX-1 inhibition

// Methylprednisolone PD
EC50_MP    = 0.25    // mg/L for 50% cytokine suppression
Emax_MP    = 0.75    // max suppression

// Infliximab PD — TNF neutralization
EC50_IFX   = 2.5     // μg/mL
Emax_IFX   = 0.90    // nearly complete TNF neutralization
n_IFX      = 1.8

// Anakinra PD — IL-1 receptor blockade
EC50_ANK   = 1.0     // μg/mL
Emax_ANK   = 0.85

// ---- Scaling ----
WT         = 15      // kg — typical KD child (~2–3yr)
IVIG_RES   = 0       // 0=responder, 1=IVIG-resistant

$CMT
// IVIG
A_IVIG_c A_IVIG_p
// Aspirin gut → plasma → salicylate
A_ASA_gut A_ASA_c A_SA_c
// Methylprednisolone
A_MP_c A_MP_p
// Infliximab
A_IFX_c A_IFX_p
// Anakinra
A_ANK_gut A_ANK_c
// Disease state compartments
IL1b IL6 TNFa Mac_act EC_act Fever CRP PLT_c CAL_Z

$INIT
A_IVIG_c = 0
A_IVIG_p = 0
A_ASA_gut = 0
A_ASA_c  = 0
A_SA_c   = 0
A_MP_c   = 0
A_MP_p   = 0
A_IFX_c  = 0
A_IFX_p  = 0
A_ANK_gut = 0
A_ANK_c  = 0
// Disease — initialise at diseased state at t=0 (symptom onset)
IL1b     = IL1_base * 6.0
IL6      = IL6_base * 8.0
TNFa     = TNF_base * 5.0
Mac_act  = Mac0 * 3.5
EC_act   = EC0 * 2.5
Fever    = 39.5
CRP      = 45.0
PLT_c    = PLT_base
CAL_Z    = 0.0

$ODE
// ---- Concentrations (per kg) ----
double C_IVIG  = A_IVIG_c / (Vc_IVIG * WT);   // g/L
double C_ASA   = A_ASA_c  / (Vc_ASA  * WT);   // mg/L
double C_SA    = A_SA_c   / (Vc_SA   * WT);   // mg/L
double C_MP    = A_MP_c   / (Vc_MP   * WT);   // mg/L
double C_IFX   = A_IFX_c  / (Vc_IFX  * WT);   // mg/L (=μg/mL × 1000 adjust below)
double C_ANK   = A_ANK_c  / (Vc_ANK  * WT);   // mg/L

// ---- Drug effect calculation (Emax models) ----
// IVIG immunomodulation
double E_IVIG  = Emax_IVIG * pow(C_IVIG, n_IVIG) /
                 (pow(EC50_IVIG, n_IVIG) + pow(C_IVIG, n_IVIG));
// IVIG additional effect in resistant patients (reduced)
if(IVIG_RES == 1) E_IVIG = E_IVIG * 0.30;

// Aspirin antipyretic
double E_ASA_pyr = Emax_ASA * C_SA / (EC50_ASA + C_SA);
// Aspirin anti-platelet (low dose, irreversible COX-1)
double E_ASA_plt = Emax_COX1 * C_ASA / (IC50_COX1 + C_ASA);

// Methylprednisolone
double E_MP    = Emax_MP * C_MP / (EC50_MP + C_MP);

// Infliximab (μg/mL ≈ mg/L for antibody concentrations)
double C_IFX_ugmL = C_IFX * 1000.0;   // convert to μg/mL
double E_IFX   = Emax_IFX * pow(C_IFX_ugmL, n_IFX) /
                 (pow(EC50_IFX, n_IFX) + pow(C_IFX_ugmL, n_IFX));

// Anakinra
double C_ANK_ugmL = C_ANK * 1000.0;
double E_ANK   = Emax_ANK * C_ANK_ugmL / (EC50_ANK + C_ANK_ugmL);

// Combined drug suppression of cytokines (maximum of individual effects)
double E_total_cyt = 1.0 - (1.0 - E_IVIG) * (1.0 - E_MP) * (1.0 - E_ANK);
double E_TNF_drug  = 1.0 - (1.0 - E_IFX)  * (1.0 - E_MP);

// ---- IVIG PK ODEs ----
dxdt_A_IVIG_c = -(CL_IVIG * WT) * C_IVIG - (Qp_IVIG * WT) * C_IVIG
                + (Qp_IVIG * WT) * (A_IVIG_p / (Vp_IVIG * WT))
                + F_FcRn * ka_FcRn * A_IVIG_c;    // FcRn recycling
dxdt_A_IVIG_p = (Qp_IVIG * WT) * C_IVIG
                - (Qp_IVIG * WT) * (A_IVIG_p / (Vp_IVIG * WT));

// ---- Aspirin PK ODEs ----
dxdt_A_ASA_gut = -ka_ASA * A_ASA_gut;
dxdt_A_ASA_c   = ka_ASA * A_ASA_gut
                 - (CL_ASA * WT) * C_ASA;
dxdt_A_SA_c    = f_SA * (CL_ASA * WT) * C_ASA
                 - (CL_SA * WT) * C_SA;

// ---- Methylprednisolone PK ODEs ----
dxdt_A_MP_c = -(CL_MP * WT) * C_MP
              - (Qp_MP * WT) * C_MP
              + (Qp_MP * WT) * (A_MP_p / (Vp_MP * WT));
dxdt_A_MP_p =  (Qp_MP * WT) * C_MP
              - (Qp_MP * WT) * (A_MP_p / (Vp_MP * WT));

// ---- Infliximab PK ODEs ----
dxdt_A_IFX_c = -(CL_IFX * WT) * C_IFX
               - (Qp_IFX * WT) * C_IFX
               + (Qp_IFX * WT) * (A_IFX_p / (Vp_IFX * WT));
dxdt_A_IFX_p =  (Qp_IFX * WT) * C_IFX
               - (Qp_IFX * WT) * (A_IFX_p / (Vp_IFX * WT));

// ---- Anakinra PK ODEs ----
dxdt_A_ANK_gut = -ka_ANK * A_ANK_gut;
dxdt_A_ANK_c   =  ka_ANK * A_ANK_gut - (CL_ANK * WT) * C_ANK;

// ---- Disease dynamics ODEs ----
// Macrophage activation
double kact_eff = kact_Mac * (1.0 - E_total_cyt);
dxdt_Mac_act = kact_eff * Mac_act * (1.0 - Mac_act / 10.0)
               - kinact_Mac * Mac_act;

// IL-1β (driven by macrophage, suppressed by IVIG/steroids/anakinra)
double IL1_drive = ksyn_IL1 * Mac_act * (1.0 - E_total_cyt) * (1.0 - E_ANK);
dxdt_IL1b = kprod_IL1 + IL1_drive * IL1_base
            - kdeg_IL1 * IL1b;

// IL-6 (driven by IL-1β and macrophage)
double IL6_drive = ksyn_IL6 * (IL1b / IL1_base) * Mac_act * (1.0 - E_total_cyt);
dxdt_IL6 = kprod_IL6 + IL6_drive * IL6_base
           - kdeg_IL6 * IL6;

// TNF-α
double TNF_drive = ksyn_TNF * Mac_act * (1.0 - E_TNF_drug);
dxdt_TNFa = kprod_TNF + TNF_drive * TNF_base
            - kdeg_TNF * TNFa;

// Endothelial activation (driven by TNF + IL1)
double EC_stim = (TNFa / TNF_base + IL1b / IL1_base) / 2.0;
dxdt_EC_act = kact_EC * EC_stim * (1.0 - E_total_cyt)
              - kinact_EC * EC_act;

// Fever (body temperature)
double fever_drive = krise_T * (IL1b/IL1_base + TNFa/TNF_base + IL6/IL6_base) / 3.0;
double antipyr = E_ASA_pyr + E_MP * 0.5;  // aspirin + steroid antipyresis
dxdt_Fever = fever_drive * (1.0 - antipyr)
             - kfall_T * (Fever - T_base) * (1.0 + E_IVIG);

// CRP (IL-6 driven acute phase)
dxdt_CRP = kprod_CRP * (IL6 / IL6_base) * (1.0 - E_total_cyt)
           - kdeg_CRP * CRP;

// Platelets (thrombocytosis: peaks at ~2 weeks)
double PLT_stim = IL6 / IL6_base;  // IL-6 drives thrombopoiesis
dxdt_PLT_c = kprod_PLT * PLT_stim - kdeg_PLT * PLT_c;

// Coronary Z-score (inflammatory endothelial damage → aneurysm)
double inflam_index = (EC_act - EC0) * (TNFa / TNF_base);
double drug_protect = E_IVIG + E_MP * 0.5 + E_IFX * 0.3;
drug_protect = (drug_protect > 1.0) ? 1.0 : drug_protect;
dxdt_CAL_Z = krise_Z * inflam_index * (1.0 - drug_protect) * (Z_max - CAL_Z)
             - kfall_Z * CAL_Z * E_IVIG;  // IVIG aids regression

$TABLE
capture C_IVIG_gL  = A_IVIG_c / (Vc_IVIG * WT);
capture C_ASA_mgL  = A_ASA_c  / (Vc_ASA  * WT);
capture C_SA_mgL   = A_SA_c   / (Vc_SA   * WT);
capture C_MP_mgL   = A_MP_c   / (Vc_MP   * WT);
capture C_IFX_ugmL = A_IFX_c  / (Vc_IFX  * WT) * 1000.0;
capture C_ANK_ugmL = A_ANK_c  / (Vc_ANK  * WT) * 1000.0;
capture Fever_C    = Fever;
capture CRP_mgL    = CRP;
capture PLT_k      = PLT_c;
capture Z_score    = CAL_Z;
capture IL1_ngmL   = IL1b;
capture IL6_ngmL   = IL6;
capture TNF_ngmL   = TNFa;
capture Mac_norm   = Mac_act;
capture EC_norm    = EC_act;
// Derived: IVIG resistance probability (logistic, Kobayashi score)
capture p_resist   = 1.0 / (1.0 + exp(-(0.8 * (Mac_act - 3.0))));
// Giant aneurysm risk flag
capture giant_CAL  = (CAL_Z >= 10.0) ? 1.0 : 0.0;
'

# Build the model
kd_mod <- mcode("kawasaki_disease", kd_code)

# ==================================================================
# Helper function: schedule IVIG dose
# ==================================================================
make_ivig_event <- function(dose_gkg = 2.0, wt_kg = 15, start_h = 0, dur_h = 12) {
  rate <- dose_gkg * wt_kg * 1000 / dur_h  # mg/h
  # IVIG goes directly to central compartment (IV infusion → cmt = A_IVIG_c)
  ev(amt = dose_gkg * wt_kg * 1000, cmt = "A_IVIG_c",
     time = start_h, rate = rate)
}

make_asa_event <- function(dose_mgkg_day = 80, wt_kg = 15, times_day = 4, n_days = 14) {
  dose_each <- dose_mgkg_day * wt_kg / times_day
  ev_times  <- seq(0, by = 24 / times_day, length.out = times_day * n_days)
  ev(amt = dose_each, cmt = "A_ASA_gut", time = ev_times)
}

make_asa_lo_event <- function(dose_mgkg_day = 5, wt_kg = 15, start_day = 14, n_days = 42) {
  dose_each <- dose_mgkg_day * wt_kg
  ev_times  <- seq(start_day * 24, by = 24, length.out = n_days)
  ev(amt = dose_each, cmt = "A_ASA_gut", time = ev_times)
}

make_mp_event <- function(dose_mgkg_day = 2, wt_kg = 15, n_days = 5) {
  dose_each <- dose_mgkg_day * wt_kg
  ev_times  <- seq(0, by = 24, length.out = n_days)
  ev(amt = dose_each, cmt = "A_MP_c", time = ev_times)
}

make_ifx_event <- function(dose_mgkg = 5, wt_kg = 15, start_h = 48) {
  ev(amt = dose_mgkg * wt_kg, cmt = "A_IFX_c", time = start_h,
     rate = dose_mgkg * wt_kg / 2)  # 2h infusion
}

make_ank_event <- function(dose_mgkg_day = 4, wt_kg = 15, n_days = 14) {
  dose_each <- dose_mgkg_day * wt_kg
  ev_times  <- seq(0, by = 24, length.out = n_days)
  ev(amt = dose_each, cmt = "A_ANK_gut", time = ev_times)
}

# ==================================================================
# Scenario 1: Standard IVIG (2 g/kg) + High-dose Aspirin → Low-dose
# ==================================================================
run_scenario1 <- function(wt = 15) {
  e_ivig <- make_ivig_event(2.0, wt, start_h = 0)
  e_asa_hi <- make_asa_event(80, wt, times_day = 4, n_days = 14)
  e_asa_lo <- make_asa_lo_event(5, wt, start_day = 14, n_days = 42)
  events <- e_ivig + e_asa_hi + e_asa_lo
  out <- mrgsim(kd_mod, events = events, end = 56 * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = 0))
  as.data.frame(out) %>% mutate(Scenario = "S1: IVIG + Aspirin (Standard)")
}

# ==================================================================
# Scenario 2: IVIG + Aspirin + Methylprednisolone (Primary Adjunct,
#             Kobayashi high-risk score ≥ 4)
# ==================================================================
run_scenario2 <- function(wt = 15) {
  e_ivig <- make_ivig_event(2.0, wt, start_h = 0)
  e_asa_hi <- make_asa_event(80, wt, 4, 14)
  e_asa_lo <- make_asa_lo_event(5, wt, 14, 42)
  e_mp   <- make_mp_event(2, wt, 5)
  events <- e_ivig + e_asa_hi + e_asa_lo + e_mp
  out <- mrgsim(kd_mod, events = events, end = 56 * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = 0))
  as.data.frame(out) %>% mutate(Scenario = "S2: IVIG + Aspirin + Steroids (High-risk)")
}

# ==================================================================
# Scenario 3: IVIG-Resistant — 2nd IVIG dose at Day 2
# ==================================================================
run_scenario3 <- function(wt = 15) {
  e_ivig1 <- make_ivig_event(2.0, wt, start_h = 0)
  e_ivig2 <- make_ivig_event(2.0, wt, start_h = 48)  # 2nd dose
  e_asa_hi <- make_asa_event(80, wt, 4, 14)
  e_asa_lo <- make_asa_lo_event(5, wt, 14, 42)
  events <- e_ivig1 + e_ivig2 + e_asa_hi + e_asa_lo
  out <- mrgsim(kd_mod, events = events, end = 56 * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = 1))
  as.data.frame(out) %>% mutate(Scenario = "S3: IVIG-Resistant → 2nd IVIG")
}

# ==================================================================
# Scenario 4: IVIG-Resistant — Infliximab rescue
# ==================================================================
run_scenario4 <- function(wt = 15) {
  e_ivig <- make_ivig_event(2.0, wt, start_h = 0)
  e_ifx  <- make_ifx_event(5, wt, start_h = 48)
  e_asa_hi <- make_asa_event(80, wt, 4, 14)
  e_asa_lo <- make_asa_lo_event(5, wt, 14, 42)
  events <- e_ivig + e_ifx + e_asa_hi + e_asa_lo
  out <- mrgsim(kd_mod, events = events, end = 56 * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = 1))
  as.data.frame(out) %>% mutate(Scenario = "S4: IVIG-Resistant → Infliximab Rescue")
}

# ==================================================================
# Scenario 5: IVIG-Resistant — Anakinra (IL-1 blockade)
# ==================================================================
run_scenario5 <- function(wt = 15) {
  e_ivig <- make_ivig_event(2.0, wt, start_h = 0)
  e_ank  <- make_ank_event(4, wt, 14)
  e_asa_hi <- make_asa_event(80, wt, 4, 14)
  e_asa_lo <- make_asa_lo_event(5, wt, 14, 42)
  events <- e_ivig + e_ank + e_asa_hi + e_asa_lo
  out <- mrgsim(kd_mod, events = events, end = 56 * 24, delta = 1,
                param = list(WT = wt, IVIG_RES = 1))
  as.data.frame(out) %>% mutate(Scenario = "S5: IVIG-Resistant → Anakinra Rescue")
}

# ==================================================================
# Run all scenarios and combine
# ==================================================================
run_all_scenarios <- function(wt = 15) {
  bind_rows(
    run_scenario1(wt),
    run_scenario2(wt),
    run_scenario3(wt),
    run_scenario4(wt),
    run_scenario5(wt)
  ) %>%
    mutate(Day = time / 24)
}

# ------------------------------------------------------------------
# Quick demo run
# ------------------------------------------------------------------
if(interactive()) {
  cat("Running KD QSP model — 5 treatment scenarios...\n")
  results <- run_all_scenarios(wt = 15)

  # --- Plot 1: Fever trajectory
  p1 <- ggplot(results, aes(Day, Fever_C, color = Scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 38.5, linetype = "dashed", color = "red") +
    labs(title = "Body Temperature (Fever)", x = "Day", y = "Temperature (°C)") +
    theme_bw() + theme(legend.position = "bottom")

  # --- Plot 2: CRP
  p2 <- ggplot(results, aes(Day, CRP_mgL, color = Scenario)) +
    geom_line(linewidth = 1) +
    labs(title = "C-Reactive Protein", x = "Day", y = "CRP (mg/L)") +
    theme_bw() + theme(legend.position = "bottom")

  # --- Plot 3: Coronary Z-score
  p3 <- ggplot(results, aes(Day, Z_score, color = Scenario)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 2.5, linetype = "dashed") +
    geom_hline(yintercept = 10,  linetype = "dashed", color = "red") +
    labs(title = "Coronary Artery Z-score", x = "Day",
         y = "Z-score", caption = "Dashed: z=2.5 (CAL threshold), z=10 (giant)") +
    theme_bw() + theme(legend.position = "bottom")

  # --- Plot 4: Cytokines
  cyt <- results %>%
    select(Day, Scenario, IL1_ngmL, IL6_ngmL, TNF_ngmL) %>%
    pivot_longer(c(IL1_ngmL, IL6_ngmL, TNF_ngmL),
                 names_to = "Cytokine", values_to = "Conc")
  p4 <- ggplot(cyt, aes(Day, Conc, color = Scenario)) +
    geom_line() +
    facet_wrap(~Cytokine, scales = "free_y") +
    labs(title = "Cytokine Dynamics", x = "Day", y = "ng/mL") +
    theme_bw() + theme(legend.position = "bottom")

  # --- Plot 5: Platelet count
  p5 <- ggplot(results, aes(Day, PLT_k, color = Scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 500, linetype = "dashed", color = "purple") +
    labs(title = "Platelet Count (Thrombocytosis)", x = "Day",
         y = "Platelets (×10³/μL)", caption = "Dashed: 500×10³ thrombocytosis threshold") +
    theme_bw() + theme(legend.position = "bottom")

  # --- Plot 6: IVIG concentration
  p6 <- results %>%
    filter(Scenario %in% c("S1: IVIG + Aspirin (Standard)",
                           "S4: IVIG-Resistant → Infliximab Rescue")) %>%
    ggplot(aes(Day, C_IVIG_gL, color = Scenario)) +
    geom_line(linewidth = 1) +
    labs(title = "IVIG Plasma Concentration", x = "Day", y = "IVIG (g/L)") +
    theme_bw() + theme(legend.position = "bottom")

  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)

  # --- Summary table
  summary_tbl <- results %>%
    group_by(Scenario) %>%
    summarise(
      Fever_Day3   = round(Fever_C[which.min(abs(Day - 3))], 1),
      CRP_Day7     = round(CRP_mgL[which.min(abs(Day - 7))], 1),
      Z_score_peak = round(max(Z_score), 2),
      PLT_peak     = round(max(PLT_k), 0),
      .groups = "drop"
    )
  print(summary_tbl)

  # --- Sensitivity analysis: weight effect on IVIG PK
  cat("\nSensitivity: body weight effect on IVIG exposure (AUC 0–14d)\n")
  wt_range <- c(8, 12, 15, 20, 25)
  sens <- lapply(wt_range, function(w) {
    r <- run_scenario1(w)
    r %>% filter(Day <= 14) %>%
      summarise(AUC_IVIG = sum(C_IVIG_gL) * (14 * 24 / nrow(.)),
                WT = w, .groups = "drop")
  }) %>% bind_rows()
  print(sens)
}
