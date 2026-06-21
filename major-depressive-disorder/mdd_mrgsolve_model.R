# =============================================================================
# Major Depressive Disorder (MDD) — QSP mrgsolve Model
# =============================================================================
# Description : Quantitative Systems Pharmacology ODE model for MDD
#               covering drug PK (escitalopram, venlafaxine, ketamine),
#               serotonin/NE/DA dynamics, HPA axis, neuroinflammation,
#               BDNF/neuroplasticity, and clinical endpoint mapping.
#
# Author      : QSP Disease Model Library (CCR automated session)
# Date        : 2026-06-20
# References  : See mdd_references.md for full citation list
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# =============================================================================
# MODEL DEFINITION
# =============================================================================

mdd_model_code <- '
$PROB
MDD QSP Model — Escitalopram SSRI prototype + multi-drug scenarios
Compartments: Oral depot, Plasma PK, SERT occupancy, 5-HT synapse,
              NE synapse, DA synapse, HPA axis (CRH/ACTH/Cortisol),
              Neuroinflammation (IL-6), BDNF, Neurogenesis, HDRS score

$PARAM @annotated
// -----------------------------------------------------------------------
// Drug PK Parameters — Escitalopram (SSRI prototype)
// Ref: Rao (2007) Clin Pharmacokinet; Herrmann (2007) Eur Neuropsychopharm
// -----------------------------------------------------------------------
DOSE_ESC   : 10    : Escitalopram oral dose [mg/day]
KA_ESC     : 0.46  : Absorption rate constant escitalopram [1/h]
CL_ESC     : 37    : Clearance escitalopram [L/h]
V_ESC      : 1090  : Volume of distribution escitalopram [L]
F_ESC      : 0.80  : Bioavailability escitalopram [-]
KI_SERT_ESC: 1.1   : Ki SERT inhibition escitalopram [nM]

// Drug PK Parameters — Venlafaxine (SNRI)
// Ref: Klamerus (1992) J Clin Pharmacol
DOSE_VEN   : 150   : Venlafaxine oral dose [mg/day]
KA_VEN     : 1.2   : Absorption rate constant venlafaxine [1/h]
CL_VEN     : 90    : Clearance venlafaxine [L/h]
V_VEN      : 500   : Volume of distribution venlafaxine [L]
F_VEN      : 0.45  : Bioavailability venlafaxine [-]
KI_SERT_VEN: 7.5   : Ki SERT inhibition venlafaxine [nM]
KI_NET_VEN : 2.7   : Ki NET inhibition venlafaxine [nM]

// Drug PK Parameters — Ketamine IV
// Ref: Clements (1982) Br J Anaesth; Zarate (2006) ARGH Psychiatry
DOSE_KET   : 0     : Ketamine IV dose [mg/kg]
KE_KET     : 1.73  : Elimination rate ketamine [1/h]  (t1/2~0.4h)
V_KET      : 3.0   : Volume of distribution ketamine [L/kg]
KI_NMDA_KET: 3.0   : Ki NMDA block ketamine [uM]
BWT        : 70    : Body weight [kg]

// -----------------------------------------------------------------------
// Serotonin Dynamics
// Ref: Dayan & Huys (2008) PLoS Comp Biol; Qi (2008) Pharmacol Biochem
// -----------------------------------------------------------------------
KSYN_5HT   : 0.12  : 5-HT synthesis rate [nM/h]
KDEG_5HT   : 0.08  : 5-HT degradation rate [1/h]
KREUP_5HT  : 0.15  : 5-HT reuptake rate (SERT) [1/h]
KREL_5HT   : 0.10  : 5-HT release rate [1/h]
SS_5HT     : 1.5   : Steady-state synaptic 5-HT [nM, baseline]
EC50_SERT  : 50    : EC50 for SERT occupancy [nM plasma]
EMAX_SERT  : 0.95  : Maximum SERT occupancy [-]
K_AUTO5HT  : 0.3   : 5-HT1A autoreceptor feedback strength [-]
DESENS_RATE: 0.005 : 5-HT1A autoreceptor desensitization rate [1/h]
RESENS_RATE: 0.002 : 5-HT1A autoreceptor resensitization rate [1/h]

// -----------------------------------------------------------------------
// Norepinephrine Dynamics
// Ref: Moret & Briley (2011) Neuropsychiatr Dis Treat
// -----------------------------------------------------------------------
KSYN_NE    : 0.08  : NE synthesis rate [nM/h]
KDEG_NE    : 0.06  : NE degradation rate [1/h]
KREUP_NE   : 0.12  : NE reuptake rate (NET) [1/h]
SS_NE      : 1.0   : Steady-state synaptic NE [nM, baseline]
EC50_NET   : 80    : EC50 for NET occupancy [nM plasma]

// -----------------------------------------------------------------------
// Dopamine Dynamics
// Ref: Stahl (2013) Stahl Psychopharmacology
// -----------------------------------------------------------------------
KSYN_DA    : 0.06  : DA synthesis rate [nM/h]
KDEG_DA    : 0.04  : DA degradation rate [1/h]
KREUP_DA   : 0.10  : DA reuptake rate (DAT) [1/h]
SS_DA      : 0.8   : Steady-state synaptic DA [nM, baseline]

// -----------------------------------------------------------------------
// HPA Axis Parameters
// Ref: Sriram (2012) J Theor Biol; Gupta (2007) IET Syst Biol
// -----------------------------------------------------------------------
KCRH_SYNTH : 0.10  : CRH synthesis rate [pg/mL/h]
KCRH_DEG   : 0.15  : CRH degradation rate [1/h]
KACTH_SYN  : 0.20  : ACTH synthesis by CRH [1/h]
KACTH_DEG  : 0.12  : ACTH degradation rate [1/h]
KCORT_SYN  : 0.30  : Cortisol synthesis by ACTH [nmol/L/h]
KCORT_DEG  : 0.05  : Cortisol degradation rate [1/h]
CORT_FB    : 0.04  : Cortisol negative feedback on CRH [1/(nmol/L)]
SS_CRH     : 0.67  : Baseline CRH [pg/mL]
SS_ACTH    : 1.67  : Baseline ACTH [pg/mL]
SS_CORT    : 15.0  : Baseline morning cortisol [nmol/L]
STRESS_INPUT: 0.0  : External stress forcing [pg/mL/h]

// -----------------------------------------------------------------------
// Neuroinflammation — IL-6
// Ref: Dowlati (2010) Biol Psychiatry; Haapakoski (2015) Brain Behav Immun
// -----------------------------------------------------------------------
KIL6_SYNTH : 0.05  : IL-6 basal synthesis rate [pg/mL/h]
KIL6_DEG   : 0.08  : IL-6 degradation [1/h]
KIL6_CORT  : 0.03  : Cortisol stimulation of IL-6 (chronic) [-]
SS_IL6     : 0.625 : Baseline IL-6 [pg/mL]
IL6_IDO_K  : 0.1   : IL-6 effect on IDO/kynurenine [pg/mL^-1]

// -----------------------------------------------------------------------
// BDNF / Neuroplasticity
// Ref: Castrén (2014) Nat Rev Neurosci; Autry & Bhattacharya (2012) Cell
// -----------------------------------------------------------------------
KBDNF_SYNTH: 0.04  : BDNF synthesis rate [ng/mL/h]
KBDNF_DEG  : 0.02  : BDNF degradation rate [1/h]
SS_BDNF    : 2.0   : Baseline serum BDNF [ng/mL]
BDNF_5HT_K : 0.3   : 5-HT stimulation of BDNF [-]
BDNF_CORT_K: 0.15  : Cortisol suppression of BDNF [-]
BDNF_IL6_K : 0.10  : IL-6 suppression of BDNF [-]
BDNF_KET_K : 2.0   : Ketamine rapid BDNF stimulation [-]

// MTOR/Neurogenesis
KMTOR_ACT  : 0.05  : mTOR activation rate by BDNF [1/h]
KMTOR_DEG  : 0.03  : mTOR deactivation [1/h]
KNEURO_SYN : 0.01  : Neurogenesis rate [1/h]
KNEURO_DEG : 0.005 : Neurogenesis decay [1/h]
SS_NEURO   : 1.0   : Baseline neurogenesis index [-]

// -----------------------------------------------------------------------
// Clinical Score Mapping
// Ref: Rush (2006) STAR*D; Oliva (2021) J Affect Disord
// -----------------------------------------------------------------------
HDRS_BASE  : 22.0  : Baseline HDRS-17 score (moderate-severe MDD)
HDRS_MAX_5HT: 6.0  : Maximum HDRS reduction from 5-HT normalization
HDRS_MAX_NE : 4.0  : Maximum HDRS reduction from NE normalization
HDRS_MAX_BDNF: 5.0 : Maximum HDRS reduction from BDNF normalization
HDRS_MAX_NEURO: 3.0: Maximum HDRS reduction from neurogenesis
HDRS_MAX_CORT: 2.0 : Maximum HDRS reduction from cortisol normalization
EC50_HDRS_5HT: 0.5 : EC50 for 5-HT effect on HDRS [nM above baseline]
WEEK_SCALE : 168   : Hours per week (for time axis)

// -----------------------------------------------------------------------
// MDD Disease Parameters (Baseline perturbations)
// -----------------------------------------------------------------------
MDD_5HT_DEF : 0.60 : MDD 5-HT deficit (fraction of normal) [-]
MDD_NE_DEF  : 0.70 : MDD NE deficit (fraction of normal) [-]
MDD_DA_DEF  : 0.65 : MDD DA deficit (fraction of normal) [-]
MDD_BDNF_DEF: 0.65 : MDD BDNF deficit (fraction of normal) [-]
MDD_CORT_EX : 1.40 : MDD cortisol excess (fraction of normal) [-]
MDD_IL6_EX  : 1.80 : MDD IL-6 elevation (fraction of normal) [-]

$CMT @annotated
// Compartments (18 total)
DEPOT_ESC  : Escitalopram oral absorption depot [mg]
CENTRAL_ESC: Escitalopram plasma concentration [mg/L]
DEPOT_VEN  : Venlafaxine oral absorption depot [mg]
CENTRAL_VEN: Venlafaxine plasma concentration [mg/L]
CENTRAL_KET: Ketamine plasma concentration [mg/L]
SERT_OCC   : SERT occupancy [fraction 0-1]
NET_OCC    : NET occupancy [fraction 0-1]
X5HT_SYN   : Synaptic serotonin concentration [nM]
NE_SYN     : Synaptic norepinephrine concentration [nM]
DA_SYN     : Synaptic dopamine concentration [nM]
CRH        : CRH (hypothalamus) [pg/mL]
ACTH       : ACTH (pituitary) [pg/mL]
CORTISOL   : Cortisol [nmol/L]
IL6        : Interleukin-6 [pg/mL]
BDNF       : Brain-derived neurotrophic factor [ng/mL]
MTOR       : mTOR activity index [-]
NEURO      : Neurogenesis index [-]
HDRS       : HDRS-17 clinical score [-]

$MAIN
// Initialize MDD disease state (perturbed from healthy baseline)
if(NEWIND <= 1) {
  // Starting values reflect MDD pathophysiology
  X5HT_SYN_0  = SS_5HT  * MDD_5HT_DEF;
  NE_SYN_0    = SS_NE   * MDD_NE_DEF;
  DA_SYN_0    = SS_DA   * MDD_DA_DEF;
  CRH_0       = SS_CRH  * MDD_CORT_EX;
  ACTH_0      = SS_ACTH * MDD_CORT_EX;
  CORTISOL_0  = SS_CORT * MDD_CORT_EX;
  IL6_0       = SS_IL6  * MDD_IL6_EX;
  BDNF_0      = SS_BDNF * MDD_BDNF_DEF;
  NEURO_0     = SS_NEURO * 0.75;
  HDRS_0      = HDRS_BASE;
  MTOR_0      = 0.5;
}

// Derived PK quantities
double Cp_ESC_nM = CENTRAL_ESC / 324.39 * 1e6;   // mg/L -> nM (MW=324.39)
double Cp_VEN_nM = CENTRAL_VEN / 277.40 * 1e6;   // mg/L -> nM (MW=277.40)
double Cp_KET_nM = CENTRAL_KET / 237.73 * 1e6;   // mg/L -> nM (MW=237.73)

// SERT occupancy (combined SSRI + SNRI; Hill equation)
double SERT_ESC  = EMAX_SERT * pow(Cp_ESC_nM, 1.5) / (pow(KI_SERT_ESC, 1.5) + pow(Cp_ESC_nM, 1.5));
double SERT_VEN  = EMAX_SERT * pow(Cp_VEN_nM, 1.0) / (pow(KI_SERT_VEN, 1.0) + pow(Cp_VEN_nM, 1.0));
double SERT_TOT  = 1.0 - (1.0 - SERT_ESC) * (1.0 - SERT_VEN);  // combined
SERT_OCC_0 = SERT_TOT;  // track occupancy

// NET occupancy
double NET_VEN   = EMAX_SERT * Cp_VEN_nM / (KI_NET_VEN + Cp_VEN_nM);
NET_OCC_0 = NET_VEN;

// Ketamine NMDA block fraction
double Cp_KET_uM = Cp_KET_nM / 1000.0;
double NMDA_block = Cp_KET_uM / (KI_NMDA_KET + Cp_KET_uM);

$ODE
// -----------------------------------------------------------------------
// PK: Escitalopram
// -----------------------------------------------------------------------
double dose_rate_ESC = DOSE_ESC / 24.0;  // daily dose -> hourly rate
dxdt_DEPOT_ESC  = F_ESC * dose_rate_ESC - KA_ESC * DEPOT_ESC;
dxdt_CENTRAL_ESC = KA_ESC * DEPOT_ESC - (CL_ESC / V_ESC) * CENTRAL_ESC;

// PK: Venlafaxine
double dose_rate_VEN = DOSE_VEN / 24.0;
dxdt_DEPOT_VEN  = F_VEN * dose_rate_VEN - KA_VEN * DEPOT_VEN;
dxdt_CENTRAL_VEN = KA_VEN * DEPOT_VEN - (CL_VEN / V_VEN) * CENTRAL_VEN;

// PK: Ketamine IV (single-dose bolus handled via event; first-order elimination)
dxdt_CENTRAL_KET = -KE_KET * CENTRAL_KET;

// -----------------------------------------------------------------------
// SERT / NET occupancy dynamics
// -----------------------------------------------------------------------
double Cp_ESC_nM_ = CENTRAL_ESC / 324.39 * 1e6;
double Cp_VEN_nM_ = CENTRAL_VEN / 277.40 * 1e6;
double Cp_KET_uM_ = CENTRAL_KET / 237.73 * 1e3;   // mg/L -> uM

double SERT_ESC_ = EMAX_SERT * pow(Cp_ESC_nM_, 1.5) / (pow(KI_SERT_ESC, 1.5) + pow(Cp_ESC_nM_, 1.5));
double SERT_VEN_ = EMAX_SERT * Cp_VEN_nM_ / (KI_SERT_VEN + Cp_VEN_nM_);
double SERT_target = 1.0 - (1.0 - SERT_ESC_) * (1.0 - SERT_VEN_);
dxdt_SERT_OCC = 0.1 * (SERT_target - SERT_OCC);  // fast equilibration

double NET_VEN_ = EMAX_SERT * Cp_VEN_nM_ / (KI_NET_VEN + Cp_VEN_nM_);
dxdt_NET_OCC = 0.1 * (NET_VEN_ - NET_OCC);

// Ketamine NMDA block
double NMDA_block_ = Cp_KET_uM_ / (KI_NMDA_KET + Cp_KET_uM_);

// -----------------------------------------------------------------------
// 5-HT Synaptic Dynamics
// Reuptake blocked by SERT occupancy
// Autoreceptor feedback included (desensitizes with chronic SSRI)
// -----------------------------------------------------------------------
double reuptake_5HT = KREUP_5HT * (1.0 - SERT_OCC) * X5HT_SYN;
double feedback_5HT = K_AUTO5HT * X5HT_SYN / (X5HT_SYN + SS_5HT);  // autoreceptor
double synth_5HT = KSYN_5HT * (1.0 - feedback_5HT * 0.5);  // feedback reduces synth
dxdt_X5HT_SYN = synth_5HT - KDEG_5HT * X5HT_SYN - reuptake_5HT;

// -----------------------------------------------------------------------
// NE Synaptic Dynamics
// -----------------------------------------------------------------------
double reuptake_NE = KREUP_NE * (1.0 - NET_OCC) * NE_SYN;
dxdt_NE_SYN = KSYN_NE - KDEG_NE * NE_SYN - reuptake_NE;

// -----------------------------------------------------------------------
// DA Synaptic Dynamics
// -----------------------------------------------------------------------
dxdt_DA_SYN = KSYN_DA - KDEG_DA * DA_SYN - KREUP_DA * DA_SYN;

// -----------------------------------------------------------------------
// HPA Axis: CRH → ACTH → Cortisol (with neg feedback)
// -----------------------------------------------------------------------
double stress_t = STRESS_INPUT;
double fb_cort = 1.0 / (1.0 + CORT_FB * CORTISOL);  // neg feedback
dxdt_CRH     = (KCRH_SYNTH + stress_t) * fb_cort - KCRH_DEG * CRH;
dxdt_ACTH    = KACTH_SYN * CRH - KACTH_DEG * ACTH;
dxdt_CORTISOL = KCORT_SYN * ACTH - KCORT_DEG * CORTISOL;

// -----------------------------------------------------------------------
// Neuroinflammation: IL-6
// Elevated by chronic stress/cortisol; antidepressants reduce modestly
// -----------------------------------------------------------------------
double il6_stress = KIL6_CORT * (CORTISOL - SS_CORT) / SS_CORT;  // cortisol-driven
double il6_antidep = 0.1 * SERT_OCC;  // SSRI mild anti-inflammatory
dxdt_IL6 = KIL6_SYNTH * (1.0 + il6_stress) - KIL6_DEG * IL6 - il6_antidep * IL6;

// -----------------------------------------------------------------------
// BDNF Dynamics
// Increased by 5-HT, NE, ketamine; suppressed by cortisol, IL-6
// -----------------------------------------------------------------------
double bdnf_5ht_stim = BDNF_5HT_K * (X5HT_SYN - SS_5HT * MDD_5HT_DEF) / SS_5HT;
double bdnf_cort_sup = BDNF_CORT_K * (CORTISOL - SS_CORT) / SS_CORT;
double bdnf_il6_sup  = BDNF_IL6_K  * (IL6 - SS_IL6) / SS_IL6;
double bdnf_ket_stim = BDNF_KET_K  * NMDA_block_;  // rapid ketamine effect
double bdnf_ne_stim  = 0.15 * (NE_SYN - SS_NE * MDD_NE_DEF) / SS_NE;
double bdnf_total_stim = bdnf_5ht_stim + bdnf_ne_stim + bdnf_ket_stim
                         - bdnf_cort_sup - bdnf_il6_sup;

dxdt_BDNF = KBDNF_SYNTH * (1.0 + bdnf_total_stim) - KBDNF_DEG * BDNF;

// -----------------------------------------------------------------------
// mTOR Activation (by BDNF and ketamine)
// -----------------------------------------------------------------------
double mtor_bdnf = KMTOR_ACT * (BDNF / SS_BDNF);
double mtor_ket  = 0.5 * NMDA_block_;  // ketamine activates mTOR rapidly
dxdt_MTOR = mtor_bdnf + mtor_ket - KMTOR_DEG * MTOR;

// -----------------------------------------------------------------------
// Neurogenesis (hippocampal) — slow process (weeks)
// -----------------------------------------------------------------------
double neuro_bdnf = KNEURO_SYN * (BDNF / SS_BDNF) * MTOR;
double neuro_cort_sup = 0.02 * (CORTISOL - SS_CORT) / SS_CORT;
dxdt_NEURO = neuro_bdnf - KNEURO_DEG * NEURO - neuro_cort_sup * NEURO;

// -----------------------------------------------------------------------
// HDRS-17 Clinical Score Mapping
// Nonlinear saturation of each component
// Score decreases (improves) with: ↑5-HT, ↑NE, ↑BDNF, ↑neurogenesis, ↓cortisol
// -----------------------------------------------------------------------
double d5HT   = X5HT_SYN - SS_5HT * MDD_5HT_DEF;
double dNE    = NE_SYN   - SS_NE  * MDD_NE_DEF;
double dBDNF  = BDNF     - SS_BDNF * MDD_BDNF_DEF;
double dNEURO = NEURO    - 0.75;
double dCORT  = CORTISOL - SS_CORT * MDD_CORT_EX;

double hdrs_5ht   = HDRS_MAX_5HT  * d5HT   / (EC50_HDRS_5HT + fabs(d5HT));
double hdrs_ne    = HDRS_MAX_NE   * dNE    / (0.3 + fabs(dNE));
double hdrs_bdnf  = HDRS_MAX_BDNF * dBDNF  / (0.5 + fabs(dBDNF));
double hdrs_neuro = HDRS_MAX_NEURO * dNEURO / (0.2 + fabs(dNEURO));
double hdrs_cort  = HDRS_MAX_CORT * (-dCORT) / (5.0 + fabs(dCORT));  // lower cortisol -> improve

double hdrs_target = HDRS_BASE - (hdrs_5ht + hdrs_ne + hdrs_bdnf + hdrs_neuro + hdrs_cort);
hdrs_target = fmax(0.0, fmin(52.0, hdrs_target));  // clamp 0-52
dxdt_HDRS = 0.005 * (hdrs_target - HDRS);  // slow convergence (score changes over weeks)

$TABLE
// Derived outputs for plotting
double PHQ9    = HDRS * 0.55;        // approximate PHQ-9 from HDRS
double MADRS_S = HDRS * 1.7;        // approximate MADRS from HDRS
double SERT_PCT = SERT_OCC * 100;
double NET_PCT  = NET_OCC  * 100;
double RESPONSE_FLAG = (HDRS <= HDRS_BASE * 0.5) ? 1.0 : 0.0;  // ≥50% reduction
double REMISSION_FLAG = (PHQ9 < 5.0) ? 1.0 : 0.0;

// Convert concentrations to clinical units
double ESC_nM   = CENTRAL_ESC / 324.39 * 1e6;
double VEN_nM   = CENTRAL_VEN / 277.40 * 1e6;
double KET_uM   = CENTRAL_KET / 237.73 * 1e3;

// Kynurenine proxy (inversely related to 5-HT via tryptophan competition)
double KYN_RATIO = (SS_IL6 * MDD_IL6_EX / IL6) * (SS_5HT / X5HT_SYN);

$CAPTURE
ESC_nM VEN_nM KET_uM SERT_PCT NET_PCT
X5HT_SYN NE_SYN DA_SYN
CORTISOL IL6 BDNF MTOR NEURO
HDRS PHQ9 MADRS_S
RESPONSE_FLAG REMISSION_FLAG KYN_RATIO
'

# Compile the model
mdd_mod <- mcode("mdd_qsp", mdd_model_code)

cat("Model compiled successfully.\n")
cat("Compartments:", nrow(init(mdd_mod)), "\n")
cat("Parameters:", length(param(mdd_mod)), "\n")

# =============================================================================
# TREATMENT SCENARIOS
# =============================================================================

# Simulation time: 0 to 8 weeks (56 days = 1344 hours)
sim_time <- seq(0, 1344, by = 2)  # every 2 hours
weeks    <- sim_time / 168

# Helper: run simulation and add scenario label
run_scenario <- function(mod, params = list(), name = "Unnamed") {
  mod %>%
    param(params) %>%
    mrgsim(end = 1344, delta = 2) %>%
    as.data.frame() %>%
    mutate(scenario = name, week = time / 168)
}

# -----------------------------------------------------------------------
# Scenario 1: No Treatment Baseline (MDD natural history)
# -----------------------------------------------------------------------
scen1 <- mdd_mod %>%
  param(DOSE_ESC = 0, DOSE_VEN = 0, DOSE_KET = 0, STRESS_INPUT = 0.02) %>%
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "1. No Treatment (MDD Baseline)", week = time / 168)

cat("Scenario 1: No Treatment — Final HDRS:", round(tail(scen1$HDRS, 1), 1), "\n")

# -----------------------------------------------------------------------
# Scenario 2: Escitalopram 10mg/day (SSRI)
# -----------------------------------------------------------------------
scen2 <- mdd_mod %>%
  param(DOSE_ESC = 10, DOSE_VEN = 0, DOSE_KET = 0, STRESS_INPUT = 0.01) %>%
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "2. Escitalopram 10mg/day (SSRI)", week = time / 168)

cat("Scenario 2: Escitalopram 10mg — Final HDRS:", round(tail(scen2$HDRS, 1), 1), "\n")

# -----------------------------------------------------------------------
# Scenario 3: Venlafaxine 150mg/day (SNRI)
# -----------------------------------------------------------------------
scen3 <- mdd_mod %>%
  param(DOSE_ESC = 0, DOSE_VEN = 150, DOSE_KET = 0, STRESS_INPUT = 0.01) %>%
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "3. Venlafaxine 150mg/day (SNRI)", week = time / 168)

cat("Scenario 3: Venlafaxine 150mg — Final HDRS:", round(tail(scen3$HDRS, 1), 1), "\n")

# -----------------------------------------------------------------------
# Scenario 4: Ketamine IV acute (0.5 mg/kg bolus at t=0)
# Using event object for single IV bolus
# -----------------------------------------------------------------------
ket_dose_mg <- 0.5 * 70  # 0.5 mg/kg × 70 kg = 35 mg
ket_V_L     <- 3.0 * 70  # L (3 L/kg × 70 kg)
ket_C0      <- ket_dose_mg / ket_V_L  # mg/L initial concentration

scen4 <- mdd_mod %>%
  param(DOSE_ESC = 0, DOSE_VEN = 0, DOSE_KET = 0, STRESS_INPUT = 0.01) %>%
  init(CENTRAL_KET = ket_C0) %>%  # IV bolus as initial condition
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "4. Ketamine IV 0.5mg/kg (single dose)", week = time / 168)

cat("Scenario 4: Ketamine IV — Final HDRS:", round(tail(scen4$HDRS, 1), 1), "\n")
cat("Scenario 4: HDRS at 24h:", round(scen4$HDRS[scen4$time == 24][1], 1), "\n")

# -----------------------------------------------------------------------
# Scenario 5: Escitalopram + Aripiprazole Augmentation
# Aripiprazole effect modeled as additional DA/5-HT modulation
# -----------------------------------------------------------------------
scen5 <- mdd_mod %>%
  param(DOSE_ESC = 10, DOSE_VEN = 0, DOSE_KET = 0, STRESS_INPUT = 0.01,
        # Aripiprazole augmentation: partial D2 agonist → ↑ DA/NE via presynaptic D2
        # modeled as 15% boost to NE and DA synthesis
        KSYN_DA = 0.06 * 1.15, KSYN_NE = 0.08 * 1.15,
        # 5-HT2A antagonism → ↑ 5-HT indirectly
        KSYN_5HT = 0.12 * 1.10) %>%
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "5. Escitalopram + Aripiprazole (Augmentation)", week = time / 168)

cat("Scenario 5: ESC+Aripiprazole — Final HDRS:", round(tail(scen5$HDRS, 1), 1), "\n")

# -----------------------------------------------------------------------
# Scenario 6: Treatment-Resistant Depression (TRD)
# High stress, blunted SERT response, elevated neuroinflammation
# -----------------------------------------------------------------------
scen6 <- mdd_mod %>%
  param(DOSE_ESC = 20,  # maximal dose
        DOSE_VEN = 0,
        DOSE_KET = 0,
        STRESS_INPUT = 0.05,  # ongoing high stress
        EMAX_SERT = 0.90,     # normal SERT pharmacology
        MDD_5HT_DEF = 0.50,   # more severe deficit
        MDD_BDNF_DEF = 0.50,  # worse neuroplasticity
        MDD_CORT_EX = 1.60,   # more cortisol dysregulation
        MDD_IL6_EX = 2.50,    # high inflammation
        HDRS_BASE = 26.0,     # severe MDD
        KIL6_SYNTH = 0.08     # higher IL-6 production
        ) %>%
  mrgsim(end = 1344, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "6. Treatment-Resistant MDD (TRD)", week = time / 168)

cat("Scenario 6: TRD — Final HDRS:", round(tail(scen6$HDRS, 1), 1), "\n")

# -----------------------------------------------------------------------
# Combine all scenarios
# -----------------------------------------------------------------------
all_scenarios <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6)

# =============================================================================
# SUMMARY TABLE
# =============================================================================

summary_table <- all_scenarios %>%
  group_by(scenario) %>%
  summarise(
    HDRS_baseline = first(HDRS),
    HDRS_wk2  = HDRS[which.min(abs(week - 2))],
    HDRS_wk4  = HDRS[which.min(abs(week - 4))],
    HDRS_wk8  = HDRS[which.min(abs(week - 8))],
    SERT_max  = max(SERT_PCT, na.rm = TRUE),
    BDNF_wk8  = BDNF[which.min(abs(week - 8))],
    .groups = "drop"
  ) %>%
  mutate(
    pct_reduction_wk8 = round((HDRS_baseline - HDRS_wk8) / HDRS_baseline * 100, 1),
    response_wk8 = pct_reduction_wk8 >= 50
  )

cat("\n=== 8-Week Summary ===\n")
print(summary_table)

# =============================================================================
# VISUALIZATION
# =============================================================================

# Color palette for scenarios
scenario_colors <- c(
  "1. No Treatment (MDD Baseline)"               = "#E74C3C",
  "2. Escitalopram 10mg/day (SSRI)"              = "#3498DB",
  "3. Venlafaxine 150mg/day (SNRI)"              = "#2ECC71",
  "4. Ketamine IV 0.5mg/kg (single dose)"        = "#9B59B6",
  "5. Escitalopram + Aripiprazole (Augmentation)"= "#F39C12",
  "6. Treatment-Resistant MDD (TRD)"             = "#1ABC9C"
)

# Plot 1: HDRS-17 over time
p1 <- ggplot(all_scenarios, aes(x = week, y = HDRS, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "gray40") +
  annotate("text", x = 7.5, y = 6.5, label = "Remission threshold (HDRS≤7)",
           size = 3, color = "gray40") +
  geom_hline(yintercept = 22 * 0.5, linetype = "dotted", color = "gray60") +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "MDD QSP Model: HDRS-17 Over 8 Weeks",
       subtitle = "Treatment scenarios comparison",
       x = "Time (weeks)", y = "HDRS-17 Score",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 2: Serotonin dynamics
p2 <- ggplot(all_scenarios, aes(x = week, y = X5HT_SYN, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 7, y = 1.55, label = "Healthy baseline", size = 3) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "Synaptic Serotonin (5-HT) Dynamics",
       x = "Time (weeks)", y = "[5-HT] Synaptic (nM)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 3: SERT occupancy (SSRI/SNRI)
p3 <- ggplot(filter(all_scenarios, scenario %in% c(
  "2. Escitalopram 10mg/day (SSRI)",
  "3. Venlafaxine 150mg/day (SNRI)",
  "5. Escitalopram + Aripiprazole (Augmentation)")),
  aes(x = week, y = SERT_PCT, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "SERT Occupancy Over Time",
       x = "Time (weeks)", y = "SERT Occupancy (%)",
       color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Plot 4: BDNF dynamics
p4 <- ggplot(all_scenarios, aes(x = week, y = BDNF, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "BDNF Dynamics Over 8 Weeks",
       x = "Time (weeks)", y = "BDNF (ng/mL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Plot 5: Cortisol
p5 <- ggplot(all_scenarios, aes(x = week, y = CORTISOL, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 15.0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "Cortisol (HPA Axis) Over 8 Weeks",
       x = "Time (weeks)", y = "Cortisol (nmol/L)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Plot 6: IL-6 Neuroinflammation
p6 <- ggplot(all_scenarios, aes(x = week, y = IL6, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.625, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = "IL-6 (Neuroinflammation) Over 8 Weeks",
       x = "Time (weeks)", y = "IL-6 (pg/mL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

cat("\nAll plots generated.\n")
cat("Summary: MDD QSP Model with", nrow(init(mdd_mod)), "ODE compartments.\n")
cat("Treatment scenarios:", length(unique(all_scenarios$scenario)), "\n")
