##############################################################################
# Migraine QSP — mrgsolve ODE Model
# Disease: Migraine (Episodic & Chronic)
# Author:  QSP Disease Model Library (CCR)
# Date:    2026-06-20
#
# Key references for parameterisation:
#  - Edvinsson et al. Nat Rev Neurol 2018 (CGRP mechanism)
#  - Olesen et al. Neurology 2009 (NTG model)
#  - Voss et al. Cephalalgia 2016 (erenumab PK/PD)
#  - Chan et al. Neurology 2021 (rimegepant acute)
#  - SAMURAI trial (lasmiditan), ARTISAN-EM (rimegepant prevention)
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ── Model code ───────────────────────────────────────────────────────────────
code <- '
$PROB Migraine QSP Model — Trigeminovascular/CGRP/CSD/Drug PK-PD

$PARAM @annotated
// ── CGRP synthesis & release ─────────────────────────────────────────────
ksyn_CGRP   : 0.05  : CGRP synthesis rate (pmol/min)
kdeg_CGRP   : 0.15  : CGRP basal degradation (1/min)
krel_CGRP   : 0.30  : CGRP release rate constant triggered by CSD/trigger (1/min)
krel_base   : 0.01  : Basal CGRP release (pmol/min)
EC50_CGRP_R : 0.50  : EC50 CGRP at CLR/RAMP1 receptor (pmol/L)
Emax_CGRP_R : 1.0   : Emax vasodilation (dimensionless)

// ── CSD propagation ───────────────────────────────────────────────────────
CSD_thresh  : 3.0   : K+ threshold for CSD initiation (mM above baseline)
kCSD_prop   : 0.10  : CSD propagation rate (1/min)
kCSD_decay  : 0.05  : CSD decay rate (1/min)
CSD_CGRP_amp: 2.0   : Amplitude factor: CSD → CGRP release

// ── Trigeminovascular sensitization ──────────────────────────────────────
kTG_act     : 0.20  : TG neuron activation by CGRP/PGE2 (1/min per unit)
kTG_decay   : 0.10  : TG activation decay (1/min)
kCS_init    : 0.05  : Central sensitization initiation rate (1/min)
kCS_decay   : 0.02  : Central sensitization decay (1/min)
CS_thresh   : 0.5   : TG activation threshold for central sensitization

// ── Pain signal ───────────────────────────────────────────────────────────
kpain_onset : 0.30  : Pain signal build-up rate (1/min)
kpain_decay : 0.05  : Spontaneous pain decay rate (1/min)
pain_max    : 10.0  : Maximum VAS pain score

// ── Nitric oxide & PGE2 ──────────────────────────────────────────────────
kNO_synth   : 0.08  : NO synthesis (baseline, pmol/min)
kNO_deg     : 0.50  : NO degradation (1/min)
kPGE2_synth : 0.06  : PGE2 synthesis (pg/min)
kPGE2_deg   : 0.20  : PGE2 degradation (1/min)
PGE2_sens   : 0.80  : PGE2 sensitization factor on TG

// ── Serotonin system ─────────────────────────────────────────────────────
kSER_base   : 0.50  : Baseline platelet 5-HT release (ng/min)
kSER_deg    : 0.30  : 5-HT degradation/reuptake (1/min)

// ── Sumatriptan SC PK (2-compartment, SC dosing) ─────────────────────────
ka_SUM      : 3.00  : SC absorption rate (1/h)
F_SUM       : 0.97  : SC bioavailability
Vc_SUM      : 2.4   : Central Vd (L/kg, 70 kg → 168 L)
Vp_SUM      : 1.8   : Peripheral Vd (L/kg, 70 kg → 126 L)
CL_SUM      : 72.0  : Total clearance (L/h, 1.2 L/min)
Q_SUM       : 20.0  : Inter-compartment CL (L/h)
MW_SUM      : 413.5 : Molecular weight sumatriptan succinate (g/mol)
EC50_1B_SUM : 3.2   : EC50 5-HT1B agonism (nM)
EC50_1D_SUM : 3.2   : EC50 5-HT1D agonism (nM)

// ── Erenumab SC PK (2-compartment, TMDD simplified) ──────────────────────
ka_ERE      : 0.012 : SC absorption rate erenumab (1/h, Tmax≈6d)
F_ERE       : 0.82  : SC bioavailability erenumab
Vc_ERE      : 3.86  : Central Vd erenumab (L)
Vp_ERE      : 2.5   : Peripheral Vd (L)
CL_ERE      : 0.14  : Linear clearance erenumab (L/h, t½≈28d)
Q_ERE       : 0.10  : Inter-compartment CL (L/h)
KD_ERE      : 0.01  : KD for CGRP-R binding (nM, ~10 pM)
CGRPR_tot   : 5.0   : Total CGRP receptor (nmol/L, tissue)

// ── Rimegepant Oral PK (gepant, 1-compartment) ───────────────────────────
ka_RIM      : 0.80  : Oral absorption rate rimegepant (1/h, Tmax≈1.5h)
F_RIM       : 0.64  : Oral bioavailability rimegepant
Vc_RIM      : 113.0 : Vd rimegepant (L)
CL_RIM      : 7.3   : Total clearance rimegepant (L/h, t½≈11h)
Ki_RIM      : 0.027 : Ki CGRP-R antagonism (nM)

// ── Topiramate PK/PD (preventive) ────────────────────────────────────────
ka_TOP      : 0.35  : Oral absorption topiramate (1/h, Tmax≈2h)
F_TOP       : 0.80  : Bioavailability topiramate
Vc_TOP      : 0.65  : Vd topiramate (L/kg, 70 kg → 45.5 L)
CL_TOP      : 0.030 : CL topiramate (L/h/kg → 2.1 L/h)
EC50_TOP    : 4000  : EC50 topiramate CSD reduction (ng/mL)

// ── Monthly migraine days (population PD) ────────────────────────────────
MMD_base    : 12.0  : Baseline monthly migraine days
MMD_min     : 0.5   : Minimum achievable MMD
kMMD_sens   : 0.5   : Sensitivity parameter MMD → pain signal

// ── Body weight ───────────────────────────────────────────────────────────
BW          : 70.0  : Body weight (kg)

$CMT @annotated
// Sumatriptan
DEPOT_SUM  : SC depot sumatriptan (mg)
CENT_SUM   : Central compartment sumatriptan (mg)
PERI_SUM   : Peripheral compartment sumatriptan (mg)

// Erenumab
DEPOT_ERE  : SC depot erenumab (mg)
CENT_ERE   : Central compartment erenumab (mg)
PERI_ERE   : Peripheral compartment erenumab (mg)
CGRPR_FREE : Free CGRP receptor (nmol/L)

// Rimegepant
CENT_RIM   : Central compartment rimegepant (mg)

// Topiramate
DEPOT_TOP  : Oral depot topiramate (mg)
CENT_TOP   : Central compartment topiramate (mg)

// Disease compartments
CGRP_TG    : CGRP in trigeminal ganglion/plasma (pmol/L)
CSD_ACT    : CSD activity state (dimensionless, 0–1)
TG_ACT     : TG neuron activation (dimensionless, 0–1)
CS_STATE   : Central sensitization state (dimensionless, 0–1)
PGE2_COMP  : PGE2 tissue level (pg/mL)
NO_COMP    : Nitric oxide level (pmol/L)
SEROTONIN  : Platelet 5-HT (ng/mL)
PAIN_SCORE : VAS pain score (0–10)

$MAIN
// Derived PK parameters (normalised to BW)
double Vc_SUM_L = Vc_SUM * BW;
double Vp_SUM_L = Vp_SUM * BW;
double Vc_TOP_L = Vc_TOP * BW;

// ── Sumatriptan concentrations (μg/L = μg/mL approximation)
double C_SUM = CENT_SUM / Vc_SUM_L;           // mg/L = μg/mL
double C_SUM_nM = C_SUM * 1000.0 / MW_SUM;    // convert to nM

// ── Erenumab concentration (mg/L)
double C_ERE  = CENT_ERE / Vc_ERE;            // mg/L

// ── Rimegepant concentration (mg/L → nM)
double C_RIM = CENT_RIM / Vc_RIM;
double MW_RIM = 534.6;
double C_RIM_nM = C_RIM * 1e6 / MW_RIM;       // nM

// ── Topiramate concentration (μg/mL)
double C_TOP = CENT_TOP / Vc_TOP_L * 1000.0;  // μg/mL

// ── 5-HT1B/1D agonism (sumatriptan) — Hill equation
double Emax_5HT1B = 1.0;
double occ_SUM_1B = Emax_5HT1B * C_SUM_nM / (EC50_1B_SUM + C_SUM_nM);
double occ_SUM_1D = Emax_5HT1B * C_SUM_nM / (EC50_1D_SUM + C_SUM_nM);

// ── CGRP receptor occupancy by erenumab (competitive binding approximation)
double KD_ERE_mgL = KD_ERE * MW_SUM / 1e6;   // rough unit conversion
double occ_ERE = C_ERE / (C_ERE + KD_ERE_mgL + 1e-12);

// ── CGRP-R occupancy by rimegepant (competitive inhibition)
// fraction blocked = C_RIM_nM / (C_RIM_nM + Ki_RIM)
double block_RIM = C_RIM_nM / (C_RIM_nM + Ki_RIM + 1e-12);

// ── Combined CGRP receptor blockade (erenumab OR rimegepant)
double CGRPR_block = 1.0 - (1.0 - occ_ERE) * (1.0 - block_RIM);

// ── Topiramate effect on CSD threshold (E_max model, inhibitory)
double Imax_TOP = 0.60;   // max 60% CSD frequency reduction
double top_inh = Imax_TOP * C_TOP / (EC50_TOP + C_TOP + 1e-12);

// ── CGRP-driven vasodilation / TG sensitization
double CGRP_eff = CGRP_TG / (EC50_CGRP_R + CGRP_TG + 1e-12);

// ── PGE2-driven TG sensitization amplification
double PGE2_amp = 1.0 + PGE2_sens * PGE2_COMP / 50.0;

// ── CSD-triggered CGRP release (sigmoid)
double CSD_frac = CSD_ACT / (CSD_ACT + 0.5 + 1e-12);

// ── Sumatriptan 5-HT1D → inhibits CGRP release
double CGRP_rel_inh = 1.0 - 0.70 * occ_SUM_1D;   // max 70% inhibition

// ── 5-HT1B vasoconstriction correction factor
double vasoconstr = 1.0 + 0.80 * occ_SUM_1B;      // counters vasodilation

// ── Pain generation (driven by TG activation and central sensitization)
double pain_drive = TG_ACT * (1.0 + 2.0 * CS_STATE) * PGE2_amp;
double pain_inh   = 0.5 * occ_SUM_1B + 0.5 * occ_SUM_1D + CGRPR_block;

$ODE
// ── Sumatriptan PK
dxdt_DEPOT_SUM = -ka_SUM * DEPOT_SUM;
dxdt_CENT_SUM  =  ka_SUM * DEPOT_SUM * F_SUM
                 - (CL_SUM + Q_SUM) / Vc_SUM_L * CENT_SUM
                 + Q_SUM / Vp_SUM_L * PERI_SUM;
dxdt_PERI_SUM  =  Q_SUM / Vc_SUM_L * CENT_SUM
                 - Q_SUM / Vp_SUM_L * PERI_SUM;

// ── Erenumab PK
dxdt_DEPOT_ERE = -ka_ERE * DEPOT_ERE;
dxdt_CENT_ERE  =  ka_ERE * DEPOT_ERE * F_ERE
                 - (CL_ERE + Q_ERE) / Vc_ERE * CENT_ERE
                 + Q_ERE / Vp_ERE * PERI_ERE;
dxdt_PERI_ERE  =  Q_ERE / Vc_ERE * CENT_ERE
                 - Q_ERE / Vp_ERE * PERI_ERE;
// Free CGRP receptor (simplified, not full TMDD)
dxdt_CGRPR_FREE = -occ_ERE * CGRPR_tot * 0.01 + 0.01 * (CGRPR_tot - CGRPR_FREE);

// ── Rimegepant PK
dxdt_CENT_RIM  =  ka_RIM * F_RIM  // dose administered separately as event
                 - CL_RIM / Vc_RIM * CENT_RIM;

// ── Topiramate PK
dxdt_DEPOT_TOP = -ka_TOP * DEPOT_TOP;
dxdt_CENT_TOP  =  ka_TOP * DEPOT_TOP * F_TOP
                 - CL_TOP / Vc_TOP_L * CENT_TOP;

// ── CGRP dynamics (trigeminal release model)
// Sources: basal synthesis, CSD-triggered release (modulated by sumatriptan 1D)
// Sink:    degradation
dxdt_CGRP_TG   =  ksyn_CGRP
                 + krel_base
                 + krel_CGRP * CSD_frac * CSD_CGRP_amp * CGRP_rel_inh
                 - kdeg_CGRP * CGRP_TG * (1.0 + 3.0 * CGRPR_block);

// ── CSD dynamics (logistic-like propagation, inhibited by topiramate)
double CSD_baseline = 0.05;   // small basal cortical activity
dxdt_CSD_ACT   =  kCSD_prop * CSD_ACT * (1.0 - CSD_ACT) * (1.0 - top_inh)
                 + CSD_baseline
                 - kCSD_decay * CSD_ACT;

// ── TG neuron activation
// Driven by CGRP, PGE2, NO, CSD
// Inhibited by sumatriptan 5-HT1D and CGRP-R blockade
dxdt_TG_ACT    =  kTG_act * (CGRP_eff * PGE2_amp + NO_COMP / 10.0) * (1.0 - CGRPR_block)
                 - kTG_decay * TG_ACT * (1.0 + occ_SUM_1D);

// ── Central sensitization
// Initiated once TG activation exceeds threshold
double CS_drive = (TG_ACT > CS_thresh) ? kCS_init * TG_ACT : 0.0;
dxdt_CS_STATE  =  CS_drive - kCS_decay * CS_STATE;

// ── PGE2 dynamics (COX pathway)
dxdt_PGE2_COMP =  kPGE2_synth * (1.0 + 2.0 * TG_ACT)   // inflammation amplifies
                 - kPGE2_deg * PGE2_COMP;

// ── NO dynamics
// NTG-triggered NO, CSD-triggered NOS, basal
dxdt_NO_COMP   =  kNO_synth * (1.0 + CSD_frac)
                 - kNO_deg * NO_COMP;

// ── Platelet serotonin
dxdt_SEROTONIN =  kSER_base - kSER_deg * SEROTONIN * (1.0 - 0.3 * occ_SUM_1B);

// ── Pain score (VAS 0–10)
// Driven by TG activation + CS; damped by drug effects
dxdt_PAIN_SCORE =  kpain_onset * pain_drive * (1.0 - pain_inh) * pain_max
                 - kpain_decay * PAIN_SCORE;

// Clamp pain score to [0, 10]
if (PAIN_SCORE < 0.0) PAIN_SCORE = 0.0;
if (PAIN_SCORE > 10.0) PAIN_SCORE = 10.0;

$TABLE
double CGRP_plasma   = CGRP_TG;
double Pain_VAS      = PAIN_SCORE;
double TG_activation = TG_ACT;
double Central_sens  = CS_STATE;
double SUM_Cp_nM     = C_SUM_nM;
double ERE_Cp_mgL    = C_ERE;
double RIM_Cp_nM     = C_RIM_nM;
double TOP_Cp_ugmL   = C_TOP;
double CGRPR_occ_ERE = occ_ERE;
double CGRPR_blk_RIM = block_RIM;
double HT5_1B_occ    = occ_SUM_1B;
double CSD_activity  = CSD_ACT;
double PGE2_level    = PGE2_COMP;
double Pain_free_2h  = (PAIN_SCORE < 0.5) ? 1.0 : 0.0;

$CAPTURE CGRP_plasma Pain_VAS TG_activation Central_sens SUM_Cp_nM
         ERE_Cp_mgL RIM_Cp_nM TOP_Cp_ugmL CGRPR_occ_ERE CGRPR_blk_RIM
         HT5_1B_occ CSD_activity PGE2_level Pain_free_2h
'

mod <- mread_cache("mgr_qsp", tempdir(), code, quiet = TRUE)

##############################################################################
# INITIAL CONDITIONS
##############################################################################
init_cond <- init(mod,
  CGRP_TG   = 0.3,    # pmol/L baseline CGRP
  CSD_ACT   = 0.1,    # low baseline cortical activity
  TG_ACT    = 0.05,   # minimal resting TG activation
  CS_STATE  = 0.0,    # no central sensitization at rest
  PGE2_COMP = 20.0,   # pg/mL baseline PGE2
  NO_COMP   = 0.5,    # pmol/L baseline NO
  SEROTONIN = 1.0,    # ng/mL baseline platelet 5-HT
  PAIN_SCORE = 0.0,   # pain-free at start
  CGRPR_FREE = 5.0    # all receptors free at baseline
)
mod <- mod %>% init(init_cond)

##############################################################################
# SCENARIO 1 — UNTREATED ACUTE MIGRAINE ATTACK
# CSD initiated at t=0; pain peaks ~2–4 h; resolves ~12–24 h
##############################################################################
cat("\n=== Scenario 1: Untreated Acute Migraine Attack ===\n")

# Trigger CSD with initial condition boost
init_attack <- init_cond
init_attack$CSD_ACT <- 0.80   # CSD fully triggered

ev_no_tx <- ev(time = 0, cmt = "CSD_ACT", amt = 0, rate = 0)

sim_untreated <- mod %>%
  init(init_attack) %>%
  mrgsim(end = 24, delta = 0.1, add = c(2, 6, 12, 24)) %>%
  as.data.frame()

cat("Peak pain VAS:", round(max(sim_untreated$Pain_VAS, na.rm = TRUE), 2), "\n")
cat("Pain at 2h:", round(sim_untreated$Pain_VAS[sim_untreated$time == 2][1], 2), "\n")
cat("Pain at 6h:", round(sim_untreated$Pain_VAS[sim_untreated$time == 6][1], 2), "\n")

##############################################################################
# SCENARIO 2 — SUMATRIPTAN SC 6 mg AT ATTACK ONSET (T=0)
##############################################################################
cat("\n=== Scenario 2: Sumatriptan SC 6 mg at Onset ===\n")

ev_sumat <- ev(time = 0, cmt = "DEPOT_SUM", amt = 6)

sim_sumat <- mod %>%
  init(init_attack) %>%
  ev(ev_sumat) %>%
  mrgsim(end = 24, delta = 0.1, add = c(2, 6, 12, 24)) %>%
  as.data.frame()

cat("Pain at 2h:", round(sim_sumat$Pain_VAS[sim_sumat$time == 2][1], 2), "\n")
cat("Pain freedom 2h:", sim_sumat$Pain_free_2h[sim_sumat$time == 2][1], "\n")

##############################################################################
# SCENARIO 3 — LASMIDITAN 200 mg ORAL AT ONSET (5-HT1F agonist)
# Modelled via: R_5HT1F → direct TNC inhibition (separate PK)
##############################################################################
cat("\n=== Scenario 3: Lasmiditan 200 mg Oral at Onset ===\n")

# Lasmiditan PK params (1-cmpt, ka=0.5/h, F=0.38, Vd=2.0 L/kg, CL=2.2 L/h/kg)
lam_param <- list(
  ka_LAS = 0.50, F_LAS = 0.38, Vc_LAS = 140.0, CL_LAS = 154.0
)
# Using surrogate: lasmiditan reduces TG_ACT by 5-HT1F mechanism
# Simplified: overlay as reduction factor on top of scenario 1
sim_lasm <- sim_untreated %>%
  mutate(
    LAS_Cp  = 200 * lam_param$F_LAS / lam_param$Vc_LAS *
               exp(-lam_param$CL_LAS / lam_param$Vc_LAS * time),
    MW_LAS  = 439.4,
    LAS_nM  = LAS_Cp * 1e6 / MW_LAS,
    EC50_LAS = 2.21,
    occ_LAS  = LAS_nM / (LAS_nM + EC50_LAS),
    # 5-HT1F effect: reduce TG activation → reduce pain
    Pain_VAS_LAS = pmax(Pain_VAS * (1 - 0.55 * occ_LAS), 0),
    Pain_free_2h_LAS = (Pain_VAS_LAS < 0.5) * 1.0
  )

cat("Pain at 2h (lasmiditan):", round(
  sim_lasm$Pain_VAS_LAS[sim_lasm$time == 2][1], 2), "\n")
cat("Pain freedom 2h:", sim_lasm$Pain_free_2h_LAS[sim_lasm$time == 2][1], "\n")

##############################################################################
# SCENARIO 4 — RIMEGEPANT 75 mg ORAL (ACUTE + PREVENTIVE)
##############################################################################
cat("\n=== Scenario 4: Rimegepant 75 mg Oral ===\n")

ev_rim <- ev(time = 0, cmt = "CENT_RIM", amt = 75 * 0.64)  # absorbed fraction

sim_rim <- mod %>%
  init(init_attack) %>%
  ev(ev_rim) %>%
  mrgsim(end = 24, delta = 0.1, add = c(2, 6, 12, 24)) %>%
  as.data.frame()

cat("Pain at 2h:", round(sim_rim$Pain_VAS[sim_rim$time == 2][1], 2), "\n")
cat("CGRP-R blockade at 1h:", round(sim_rim$CGRPR_blk_RIM[sim_rim$time == 1][1], 4), "\n")
cat("CGRP-R blockade at 6h:", round(sim_rim$CGRPR_blk_RIM[sim_rim$time == 6][1], 4), "\n")

##############################################################################
# SCENARIO 5 — ERENUMAB 140 mg SC MONTHLY PREVENTION (3-month simulation)
##############################################################################
cat("\n=== Scenario 5: Erenumab 140 mg SC Monthly Prevention (3 months) ===\n")

# Three monthly doses: day 0, 28, 56
ev_ere3 <- ev(time = 0,  cmt = "DEPOT_ERE", amt = 140) +
           ev(time = 672, cmt = "DEPOT_ERE", amt = 140) +   # 28 days
           ev(time = 1344, cmt = "DEPOT_ERE", amt = 140)    # 56 days

# Reset to healthy baseline for prevention scenario
init_healthy <- init_cond
sim_ere <- mod %>%
  init(init_healthy) %>%
  ev(ev_ere3) %>%
  mrgsim(end = 2016, delta = 1.0) %>%   # 84 days in hours
  as.data.frame()

cat("Erenumab Cmax (peak day ~6):",
    round(max(sim_ere$ERE_Cp_mgL, na.rm = TRUE), 3), "mg/L\n")
cat("CGRP-R occupancy at 1 month:",
    round(sim_ere$CGRPR_occ_ERE[sim_ere$time == 672][1], 3), "\n")
cat("CGRP-R occupancy at 3 months:",
    round(sim_ere$CGRPR_occ_ERE[sim_ere$time == 2016][1], 3), "\n")

##############################################################################
# SCENARIO 6 — TOPIRAMATE 100 mg/day PREVENTION (steady-state)
##############################################################################
cat("\n=== Scenario 6: Topiramate 100 mg/day Prevention (30-day) ===\n")

# Twice-daily dosing: 50 mg every 12 h
ev_top <- ev_ss(
  amt = 50, ii = 12, cmt = "DEPOT_TOP"
)

sim_top <- mod %>%
  init(init_healthy) %>%
  ev(ev_top) %>%
  mrgsim(end = 720, delta = 0.5) %>%   # 30 days
  as.data.frame()

cat("Topiramate Css (steady-state approx):",
    round(mean(tail(sim_top$TOP_Cp_ugmL, 48), na.rm = TRUE), 2), "μg/mL\n")
top_inh_ss <- with(tail(sim_top, 1),
  0.60 * TOP_Cp_ugmL / (4000 + TOP_Cp_ugmL))
cat("CSD inhibition at SS:", round(top_inh_ss, 3), "\n")

##############################################################################
# SCENARIO 7 — CHRONIC MIGRAINE: PROGRESSION & CGRP mAb (1-year simulation)
# Model: MMD increases when CS_STATE accumulates; erenumab halts progression
##############################################################################
cat("\n=== Scenario 7: Chronic Migraine Progression vs Erenumab (1 year) ===\n")

# Simulate uncontrolled chronic migraine: MMD ≥ 15
# Proxy: high CSD_ACT baseline representing medication overuse / kindling
init_chronic <- init_cond
init_chronic$CSD_ACT   <- 0.60
init_chronic$CS_STATE  <- 0.40
init_chronic$TG_ACT    <- 0.35
init_chronic$PAIN_SCORE <- 6.0

# Untreated chronic migraine (1 year)
sim_chr_untx <- mod %>%
  init(init_chronic) %>%
  mrgsim(end = 8760, delta = 6.0) %>%
  as.data.frame() %>%
  mutate(day = time / 24, scenario = "Untreated Chronic")

# With erenumab 140 mg monthly
ev_ere_1yr <- lapply(0:12, function(i)
  ev(time = i * 730, cmt = "DEPOT_ERE", amt = 140)) %>%
  do.call(c, .)

sim_chr_ere <- mod %>%
  init(init_chronic) %>%
  ev(ev_ere_1yr) %>%
  mrgsim(end = 8760, delta = 6.0) %>%
  as.data.frame() %>%
  mutate(day = time / 24, scenario = "Erenumab 140 mg Monthly")

sim_scen7 <- bind_rows(sim_chr_untx, sim_chr_ere)
cat("Mean pain (untreated):",
    round(mean(sim_chr_untx$Pain_VAS, na.rm = TRUE), 2), "\n")
cat("Mean pain (erenumab):",
    round(mean(sim_chr_ere$Pain_VAS, na.rm = TRUE), 2), "\n")
cat("Pain reduction %:",
    round((1 - mean(sim_chr_ere$Pain_VAS) / mean(sim_chr_untx$Pain_VAS)) * 100, 1), "%\n")

##############################################################################
# PLOTS
##############################################################################

# Combine acute scenarios
sim_acute <- bind_rows(
  sim_untreated %>% mutate(scenario = "Untreated"),
  sim_sumat     %>% mutate(scenario = "Sumatriptan SC 6 mg"),
  sim_rim       %>% mutate(scenario = "Rimegepant 75 mg"),
  sim_lasm      %>%
    select(time, Pain_VAS = Pain_VAS_LAS) %>%
    mutate(scenario = "Lasmiditan 200 mg",
           TG_activation = NA, Central_sens = NA, CGRP_plasma = NA,
           CSD_activity = NA)
)

p_acute <- ggplot(sim_acute %>% filter(time <= 24),
       aes(time, Pain_VAS, colour = scenario, linetype = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 23, y = 0.7, label = "Pain freedom threshold", size = 3) +
  labs(title = "Acute Migraine Attack — Pain VAS Over Time",
       subtitle = "Sumatriptan SC vs Lasmiditan vs Rimegepant vs Untreated",
       x = "Time (hours)", y = "VAS Pain Score (0–10)",
       colour = "Treatment", linetype = "Treatment") +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

p_cgrp <- ggplot(sim_acute %>% filter(time <= 24, !is.na(CGRP_plasma)),
       aes(time, CGRP_plasma, colour = scenario)) +
  geom_line(size = 1.0) +
  labs(title = "Plasma CGRP Dynamics",
       x = "Time (hours)", y = "CGRP (pmol/L)",
       colour = "Treatment") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

p_ere_pk <- ggplot(sim_ere, aes(time / 24, ERE_Cp_mgL)) +
  geom_line(colour = "#2c7bb6", size = 1.1) +
  geom_line(aes(y = CGRPR_occ_ERE * 5), colour = "#d7191c",
            linetype = "dashed", size = 1.0) +
  scale_y_continuous(
    name = "Erenumab (mg/L)",
    sec.axis = sec_axis(~./5, name = "CGRP-R Occupancy (fraction)")
  ) +
  labs(title = "Erenumab PK/PD — 3-Month Monthly Dosing",
       x = "Time (days)") +
  theme_bw(base_size = 13)

p_chronic <- ggplot(sim_scen7,
       aes(day, Pain_VAS, colour = scenario)) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, size = 1.1) +
  labs(title = "Chronic Migraine: Untreated vs Erenumab Prevention (1 Year)",
       x = "Day", y = "VAS Pain Score",
       colour = "Scenario") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# Save plots
plot_list <- list(
  acute    = p_acute,
  cgrp     = p_cgrp,
  ere_pk   = p_ere_pk,
  chronic  = p_chronic
)

for (nm in names(plot_list)) {
  ggsave(
    filename = file.path("plots", paste0("mgr_", nm, ".png")),
    plot     = plot_list[[nm]],
    width = 10, height = 6, dpi = 150
  )
}

cat("\n========================================================\n")
cat("All 7 scenarios complete. Plots saved to plots/ directory.\n")
cat("========================================================\n")
