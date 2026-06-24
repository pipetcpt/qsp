## ============================================================
##  Fibromyalgia (FM) — mrgsolve QSP Model
##  Compartments:
##    Drug PK  : Duloxetine (DUL), Pregabalin (PRE),
##               Milnacipran (MIL), Amitriptyline (TCA)
##    Peripheral: DRG sensitization, NGF, PGE2
##    Spinal   : Dorsal-horn WDR, Substance P (CSF),
##               NMDA-receptor state, Wind-up, LTP (central sensitization)
##    Brain    : Synaptic NE, 5-HT, descending inhibition (DPMS)
##    HPA Axis : CRH, ACTH, Cortisol, negative feedback
##    ANS      : SNS tone, HRV
##    Sleep    : Sleep pressure (adenosine), SWS depth
##    Immune   : Microglia activation, IL-1β (spinal)
##    Outcomes : Pain score, FIQ, fatigue, depression
##  Total: 30 ODE compartments
## ============================================================

library(mrgsolve)

fm_code <- '
$PROB
Fibromyalgia QSP — Duloxetine / Pregabalin / Milnacipran / Amitriptyline PK-PD
30-compartment ODE model including:
  Drug PK (4 drugs), Peripheral sensitization, Central sensitization (spinal/brain),
  HPA axis, ANS, Sleep dynamics, Neuroinflammation, Clinical outcomes

$PARAM
// ---- Drug PK ----
ka_DUL  = 0.80   // h^-1   oral absorption, duloxetine
CL_DUL  = 54.0   // L/h    total clearance
V1_DUL  = 1640   // L      central volume
Q_DUL   = 18.0   // L/h    inter-compartmental clearance
V2_DUL  = 820    // L      peripheral volume
F_DUL   = 0.50   // bioavailability

ka_PRE  = 1.30   // h^-1   pregabalin (rapid absorption)
CL_PRE  = 6.8    // L/h    renal clearance (CLcr-based)
V_PRE   = 42.0   // L      1-compartment
F_PRE   = 0.90

ka_MIL  = 0.90   // h^-1
CL_MIL  = 50.0   // L/h
V_MIL   = 300    // L
F_MIL   = 0.85

ka_TCA  = 0.70   // h^-1   amitriptyline
CL_TCA  = 40.0   // L/h
V_TCA   = 1500   // L
F_TCA   = 0.48

// ---- Drug PD: SNRI (Duloxetine, Milnacipran) ----
IC50_SERT_DUL  = 0.003  // mg/L  SERT IC50
IC50_NET_DUL   = 0.012  // mg/L  NET IC50
IC50_SERT_MIL  = 0.015  // mg/L
IC50_NET_MIL   = 0.018  // mg/L
IC50_SERT_TCA  = 0.010  // mg/L
IC50_NET_TCA   = 0.025  // mg/L
Emax_SNRI      = 1.0    // dimensionless

// ---- Drug PD: Pregabalin (alpha2delta) ----
IC50_PRE_alpha2d = 0.05 // mg/L
Emax_PRE         = 0.70 // max fractional block of Ca-channel

// ---- Peripheral sensitization ----
kprod_NGF   = 0.05   // baseline NGF synthesis h^-1
kdeg_NGF    = 0.10   // NGF degradation h^-1
kact_TRPV1  = 0.20   // PGE2/NGF → TRPV1 sensitization
kdeg_PGE2   = 0.30   // h^-1
kprod_PGE2  = 0.08   // h^-1 baseline
kstim_DRG   = 0.40   // DRG firing rate constant
kdeg_DRG    = 0.50   // DRG decay h^-1

// ---- Spinal / Central sensitization ----
kprod_SP    = 0.15   // Substance P production h^-1
kdeg_SP     = 0.20   // h^-1
kWU         = 0.05   // wind-up accumulation h^-1
kWU_decay   = 0.02   // wind-up decay h^-1
kLTP        = 0.08   // LTP induction rate h^-1
kLTP_decay  = 0.005  // LTP spontaneous decay h^-1
kNMDA_act   = 0.12   // NMDA receptor activation by Glu/SP
kNMDA_decay = 0.08   // NMDA inactivation h^-1
Emax_inhib  = 0.70   // max inhibitory effect of descending DPMS

// ---- Synaptic transmitters (brain) ----
ksyn_NE     = 0.30   // NE synthesis/pool h^-1
ksyn_5HT    = 0.25   // 5-HT synthesis h^-1
kdeg_NE     = 0.40   // NE reuptake/degradation h^-1
kdeg_5HT    = 0.35   // 5-HT reuptake/degradation h^-1
kdesc_NE    = 0.20   // descending NE → spinal inhibition
kdesc_5HT   = 0.15   // descending 5-HT → spinal inhibition

// ---- HPA axis ----
kprod_CRH   = 0.50   // h^-1
kdeg_CRH    = 0.80   // h^-1
kprod_ACTH  = 0.40   // h^-1
kdeg_ACTH   = 0.60   // h^-1
kprod_CORT  = 0.30   // h^-1
kdeg_CORT   = 0.20   // h^-1
kfb_CORT    = 0.60   // cortisol negative feedback strength

// ---- ANS ----
kSNS_base   = 0.60   // baseline SNS tone
kSNS_stress = 0.20   // stress-to-SNS coupling
kSNS_decay  = 0.40   // h^-1
kHRV_base   = 0.80   // baseline HRV (normalized)

// ---- Sleep ----
kaden_prod  = 0.15   // adenosine production h^-1
kaden_clear = 0.12   // adenosine clearance h^-1
kSWS_drive  = 0.20   // adenosine → SWS drive
kSWS_decay  = 0.25   // SWS decay h^-1
kSWS_pain_inh = 0.30 // pain → SWS inhibition coefficient
kSWS_TCA    = 0.15   // TCA sedation → SWS support

// ---- Neuroinflammation ----
kprod_MG    = 0.10   // microglia activation h^-1
kdeg_MG     = 0.12   // h^-1
kprod_IL1b  = 0.18   // IL-1β production h^-1
kdeg_IL1b   = 0.25   // h^-1
kMG_cortisol = 0.15  // cortisol → microglia suppression

// ---- Clinical outcomes ----
k_pain_LTP  = 0.40   // LTP → pain score scaling
k_pain_SP   = 0.20   // CSF SP → pain score
k_FIQ_pain  = 0.35   // pain → FIQ
k_FIQ_sleep = 0.25   // sleep deprivation → FIQ
k_FIQ_dep   = 0.20   // depression → FIQ
k_fatigue   = 0.30   // sleep loss → fatigue
k_dep_LTP   = 0.15   // central sensitization → depression
pain_base   = 5.0    // baseline NRS pain score
FIQ_base    = 55.0   // baseline FIQR score
fatigue_base = 60.0  // baseline fatigue VAS

// ---- Dosing flags ----
use_DUL = 0     // 0=off, 1=on
use_PRE = 0
use_MIL = 0
use_TCA = 0

$CMT
// Drug PK compartments (10)
DUL_gut DUL_cent DUL_peri   // duloxetine: gut, central, peripheral
PRE_gut PRE_cent             // pregabalin: gut, central
MIL_gut MIL_cent             // milnacipran
TCA_gut TCA_cent             // amitriptyline

// Peripheral sensitization (3)
NGF PGE2 DRG_act

// Spinal (6)
SP_csf NMDA_state WindUp LTP_cs NE_syn SHT_syn

// HPA axis (3)
CRH ACTH CORT

// ANS + Sleep (4)
SNS_tone SWS_depth Adenosine DPMS

// Neuroinflammation (2)
MG_act IL1b_sp

// Clinical outcomes (4) -- tracked as ODEs for smoothing
Pain_score FIQ_score Fatigue_VAS Depression_score

$MAIN
// ----- PK initial steady-state approximations -----
// Pain begins at FM steady state (untreated)
if(NEWIND <= 1) {
  // start at baseline disease state
}

$ODE
// ============================================================
// Drug PK
// ============================================================
// -- Duloxetine --
double ka_DUL_eff = use_DUL * ka_DUL;
dxdt_DUL_gut  = -ka_DUL_eff * DUL_gut;
dxdt_DUL_cent =  ka_DUL_eff * DUL_gut * F_DUL
                 - (CL_DUL/V1_DUL) * DUL_cent
                 - (Q_DUL/V1_DUL)  * DUL_cent
                 + (Q_DUL/V2_DUL)  * DUL_peri;
dxdt_DUL_peri =  (Q_DUL/V1_DUL)   * DUL_cent
                 - (Q_DUL/V2_DUL)  * DUL_peri;

// Plasma concentrations (mg/L)
double Cp_DUL  = DUL_cent / V1_DUL;

// -- Pregabalin --
double ka_PRE_eff = use_PRE * ka_PRE;
dxdt_PRE_gut  = -ka_PRE_eff * PRE_gut;
dxdt_PRE_cent =  ka_PRE_eff * PRE_gut * F_PRE - (CL_PRE/V_PRE) * PRE_cent;
double Cp_PRE  = PRE_cent / V_PRE;

// -- Milnacipran --
double ka_MIL_eff = use_MIL * ka_MIL;
dxdt_MIL_gut  = -ka_MIL_eff * MIL_gut;
dxdt_MIL_cent =  ka_MIL_eff * MIL_gut * F_MIL - (CL_MIL/V_MIL) * MIL_cent;
double Cp_MIL  = MIL_cent / V_MIL;

// -- Amitriptyline --
double ka_TCA_eff = use_TCA * ka_TCA;
dxdt_TCA_gut  = -ka_TCA_eff * TCA_gut;
dxdt_TCA_cent =  ka_TCA_eff * TCA_gut * F_TCA - (CL_TCA/V_TCA) * TCA_cent;
double Cp_TCA  = TCA_cent / V_TCA;

// ============================================================
// Drug PD — fractional inhibition
// ============================================================
// SERT inhibition (combined)
double inh_SERT = Emax_SNRI * (
    Cp_DUL / (IC50_SERT_DUL + Cp_DUL) +
    Cp_MIL / (IC50_SERT_MIL + Cp_MIL) +
    Cp_TCA / (IC50_SERT_TCA + Cp_TCA)
);
inh_SERT = (inh_SERT > 1.0) ? 1.0 : inh_SERT;

// NET inhibition
double inh_NET  = Emax_SNRI * (
    Cp_DUL / (IC50_NET_DUL + Cp_DUL) +
    Cp_MIL / (IC50_NET_MIL + Cp_MIL) +
    Cp_TCA / (IC50_NET_TCA + Cp_TCA)
);
inh_NET = (inh_NET > 1.0) ? 1.0 : inh_NET;

// Alpha2-delta channel block (pregabalin)
double inh_Ca   = Emax_PRE * Cp_PRE / (IC50_PRE_alpha2d + Cp_PRE);

// TCA sedation → SWS support
double eff_TCA_sleep = Cp_TCA / (0.05 + Cp_TCA);

// ============================================================
// Peripheral Sensitization
// ============================================================
// NGF dynamics (NGF elevated in FM muscle/skin)
dxdt_NGF  = kprod_NGF * 1.5          // elevated production in FM
             - kdeg_NGF * NGF;

// PGE2 dynamics (mast cell / COX pathway)
dxdt_PGE2 = kprod_PGE2 * (1.0 + IL1b_sp * 0.5)  // IL-1β amplifies PGE2
             - kdeg_PGE2 * PGE2;

// DRG afferent activity (A-delta/C fibers)
double TRPV1_sens = kact_TRPV1 * (PGE2 + NGF * 0.5);
dxdt_DRG_act = kstim_DRG * TRPV1_sens
               - kdeg_DRG * DRG_act;

// ============================================================
// Spinal Dorsal Horn — Central Sensitization
// ============================================================
// Substance P (CSF proxy)
dxdt_SP_csf = kprod_SP * DRG_act * (1.0 - inh_Ca * 0.6)
              - kdeg_SP * SP_csf;

// NMDA receptor state (0–1 scale, 0=resting)
dxdt_NMDA_state = kNMDA_act * (SP_csf + DRG_act * 0.5) * (1.0 - NMDA_state)
                  - kNMDA_decay * NMDA_state;

// Wind-up (WU): repetitive C-fiber activity → temporal summation
dxdt_WindUp = kWU * DRG_act * NMDA_state
              - kWU_decay * WindUp;

// Long-term potentiation / central sensitization index (0–1)
dxdt_LTP_cs = kLTP * WindUp * (1.0 + IL1b_sp * 0.4)
              - kLTP_decay * LTP_cs
              - LTP_cs * DPMS * Emax_inhib; // descending inhibition

// ============================================================
// Supraspinal — Synaptic Monoamines & Descending Inhibition
// ============================================================
// Synaptic NE (descending)
dxdt_NE_syn  = ksyn_NE  * (1.0 + inh_NET  * 2.0)  // NET block → ↑NE
               - kdeg_NE * (1.0 - inh_NET) * NE_syn;

// Synaptic 5-HT (descending)
dxdt_SHT_syn = ksyn_5HT * (1.0 + inh_SERT * 2.0)  // SERT block → ↑5-HT
               - kdeg_5HT * (1.0 - inh_SERT) * SHT_syn;

// Descending Pain Modulating System (DPMS, 0–1 scale)
dxdt_DPMS    = kdesc_NE * NE_syn + kdesc_5HT * SHT_syn
               - 0.30 * DPMS      // intrinsic turnover
               - DPMS * SNS_tone * 0.10;  // SNS stress partially offsets

// ============================================================
// HPA Axis
// ============================================================
// CRH (hypothalamus) — stress + pain driven
dxdt_CRH  = kprod_CRH * (1.0 + LTP_cs * 0.5 + SNS_tone * 0.3)
             - kdeg_CRH  * CRH
             - kfb_CORT  * CORT * CRH;  // negative feedback

// ACTH (pituitary)
dxdt_ACTH = kprod_ACTH * CRH - kdeg_ACTH * ACTH;

// Cortisol (adrenal)
dxdt_CORT  = kprod_CORT * ACTH - kdeg_CORT * CORT;

// ============================================================
// ANS — Sympathetic tone
// ============================================================
dxdt_SNS_tone = kSNS_base + kSNS_stress * (LTP_cs + 0.5 * (1.0 - CORT * 0.5))
                - kSNS_decay * SNS_tone;

// ============================================================
// Sleep — SWS depth and adenosine
// ============================================================
// Adenosine (sleep pressure)
dxdt_Adenosine = kaden_prod - kaden_clear * Adenosine
                 + SNS_tone * 0.05;  // stress slows clearance

// Slow-wave sleep depth (0–1)
dxdt_SWS_depth = kSWS_drive * Adenosine
                 - kSWS_decay * SWS_depth
                 - kSWS_pain_inh * LTP_cs * SWS_depth   // pain disrupts SWS
                 + kSWS_TCA * eff_TCA_sleep;              // TCA sedation

// ============================================================
// Neuroinflammation
// ============================================================
// Microglia activation (spinal)
dxdt_MG_act  = kprod_MG * DRG_act * (1.0 + WindUp * 0.5)
               - kdeg_MG * MG_act
               - kMG_cortisol * CORT * MG_act;

// IL-1β (spinal)
dxdt_IL1b_sp = kprod_IL1b * MG_act - kdeg_IL1b * IL1b_sp;

// ============================================================
// Clinical Outcomes — smoothed ODE representations
// ============================================================
// Pain NRS (0–10): driven by LTP, SP, offset by DPMS + drug effect
double pain_target = pain_base
                     + k_pain_LTP * LTP_cs * 6.0        // max +6
                     + k_pain_SP  * SP_csf * 2.0
                     - DPMS * 4.0                         // descending inhibition
                     - (inh_NET + inh_SERT) * 2.0;       // SNRI analgesic
pain_target = (pain_target < 0) ? 0 : (pain_target > 10 ? 10 : pain_target);
dxdt_Pain_score    = 0.15 * (pain_target - Pain_score);

// FIQR score (0–100): composite
double FIQ_target  = FIQ_base
                     + k_FIQ_pain  * (Pain_score - 5.0) * 8.0
                     + k_FIQ_sleep * (1.0 - SWS_depth) * 20.0
                     + k_FIQ_dep   * Depression_score * 0.5;
FIQ_target = (FIQ_target < 0) ? 0 : (FIQ_target > 100 ? 100 : FIQ_target);
dxdt_FIQ_score     = 0.10 * (FIQ_target - FIQ_score);

// Fatigue VAS (0–100)
double fatigue_tgt = fatigue_base
                     + k_fatigue * (1.0 - SWS_depth) * 30.0
                     - inh_NET * 20.0;    // NE↑ reduces fatigue
fatigue_tgt = (fatigue_tgt < 0) ? 0 : (fatigue_tgt > 100 ? 100 : fatigue_tgt);
dxdt_Fatigue_VAS   = 0.12 * (fatigue_tgt - Fatigue_VAS);

// Depression (PHQ-9 scaled 0–27)
double dep_target  = 8.0
                     + k_dep_LTP * LTP_cs * 12.0
                     + (1.0 - SWS_depth) * 4.0
                     - (inh_SERT + inh_NET) * 5.0;
dep_target = (dep_target < 0) ? 0 : (dep_target > 27 ? 27 : dep_target);
dxdt_Depression_score = 0.08 * (dep_target - Depression_score);

$TABLE
double Cp_DUL_out = DUL_cent / V1_DUL;
double Cp_PRE_out = PRE_cent / V_PRE;
double Cp_MIL_out = MIL_cent / V_MIL;
double Cp_TCA_out = TCA_cent / V_TCA;
double inh_SERT_pct = 100 * (Cp_DUL_out/(IC50_SERT_DUL+Cp_DUL_out) +
                              Cp_MIL_out/(IC50_SERT_MIL+Cp_MIL_out));
double inh_NET_pct  = 100 * (Cp_DUL_out/(IC50_NET_DUL+Cp_DUL_out)  +
                              Cp_MIL_out/(IC50_NET_MIL+Cp_MIL_out));
double Ca_block_pct = 100 * Emax_PRE * Cp_PRE_out/(IC50_PRE_alpha2d+Cp_PRE_out);

$CAPTURE
Cp_DUL_out Cp_PRE_out Cp_MIL_out Cp_TCA_out
inh_SERT_pct inh_NET_pct Ca_block_pct
SP_csf NMDA_state WindUp LTP_cs DPMS
NE_syn SHT_syn
CRH ACTH CORT SNS_tone
SWS_depth Adenosine
MG_act IL1b_sp
Pain_score FIQ_score Fatigue_VAS Depression_score
'

## ============================================================
## Build and compile model
## ============================================================
fm_mod <- mcode("fibromyalgia_qsp", fm_code)

## ============================================================
## Initial conditions — FM steady state (untreated)
## ============================================================
FM_init <- init(fm_mod,
  NGF      = 1.5,    # elevated in FM
  PGE2     = 1.2,
  DRG_act  = 0.8,
  SP_csf   = 1.8,    # high CSF SP
  NMDA_state = 0.4,
  WindUp   = 0.3,
  LTP_cs   = 0.55,   # significant central sensitization
  NE_syn   = 0.7,
  SHT_syn  = 0.6,
  DPMS     = 0.35,   # blunted descending inhibition
  CRH      = 0.9,
  ACTH     = 0.85,
  CORT     = 0.8,
  SNS_tone = 0.75,   # elevated SNS
  Adenosine = 1.2,
  SWS_depth = 0.35,  # poor SWS
  MG_act   = 0.65,
  IL1b_sp  = 0.55,
  Pain_score   = 6.5,
  FIQ_score    = 68,
  Fatigue_VAS  = 72,
  Depression_score = 12
)

## ============================================================
## Dosing Events
## ============================================================
# Duloxetine 60 mg QD
dose_DUL_60 <- ev(amt = 60, ii = 24, addl = 83, cmt = "DUL_gut",
                  param = list(use_DUL = 1))

# Pregabalin 150 mg BID
dose_PRE_150bid <- ev(amt = 150, ii = 12, addl = 167, cmt = "PRE_gut",
                      param = list(use_PRE = 1))

# Milnacipran 50 mg BID
dose_MIL_50bid <- ev(amt = 50, ii = 12, addl = 167, cmt = "MIL_gut",
                     param = list(use_MIL = 1))

# Amitriptyline 25 mg QHS
dose_TCA_25 <- ev(amt = 25, ii = 24, addl = 83, cmt = "TCA_gut", time = 22,
                  param = list(use_TCA = 1))

# Combination: DUL + PRE
dose_combo <- ev_seq(dose_DUL_60, dose_PRE_150bid)

## ============================================================
## Simulation Scenarios
## ============================================================
sim_time <- seq(0, 84 * 24, by = 1)  # 12 weeks in hours

# Scenario 1: Untreated FM baseline
sim_base <- fm_mod %>%
  init(FM_init) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Untreated FM", time_days = time / 24)

# Scenario 2: Duloxetine 60 mg QD
sim_DUL <- fm_mod %>%
  init(FM_init) %>%
  param(use_DUL = 1) %>%
  ev(dose_DUL_60) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Duloxetine 60 mg QD", time_days = time / 24)

# Scenario 3: Pregabalin 150 mg BID
sim_PRE <- fm_mod %>%
  init(FM_init) %>%
  param(use_PRE = 1) %>%
  ev(dose_PRE_150bid) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Pregabalin 150mg BID", time_days = time / 24)

# Scenario 4: Milnacipran 50 mg BID
sim_MIL <- fm_mod %>%
  init(FM_init) %>%
  param(use_MIL = 1) %>%
  ev(dose_MIL_50bid) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Milnacipran 50mg BID", time_days = time / 24)

# Scenario 5: Duloxetine + Pregabalin combination
sim_COMBO <- fm_mod %>%
  init(FM_init) %>%
  param(use_DUL = 1, use_PRE = 1) %>%
  ev(dose_combo) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "DUL + PRE Combo", time_days = time / 24)

# Scenario 6: Low-dose Amitriptyline (sleep-targeted)
sim_TCA <- fm_mod %>%
  init(FM_init) %>%
  param(use_TCA = 1) %>%
  ev(dose_TCA_25) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Amitriptyline 25mg QHS", time_days = time / 24)

# Combine all scenarios
sim_all <- bind_rows(sim_base, sim_DUL, sim_PRE, sim_MIL, sim_COMBO, sim_TCA)

## ============================================================
## Visualization
## ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)

colors_scen <- c(
  "Untreated FM"         = "#E74C3C",
  "Duloxetine 60 mg QD"  = "#3498DB",
  "Pregabalin 150mg BID" = "#2ECC71",
  "Milnacipran 50mg BID" = "#9B59B6",
  "DUL + PRE Combo"      = "#E67E22",
  "Amitriptyline 25mg QHS" = "#1ABC9C"
)

# Plot 1: Pain score over time
p1 <- ggplot(sim_all, aes(time_days, Pain_score, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "NRS Pain Score (0-10) — 12-Week Treatment",
       x = "Time (days)", y = "Pain NRS", color = "Scenario") +
  theme_bw(14) + ylim(0, 10) +
  geom_hline(yintercept = c(3, 5), linetype = "dashed", color = "grey50")

# Plot 2: FIQ score
p2 <- ggplot(sim_all, aes(time_days, FIQ_score, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "FIQR Score (0-100) — Functional Impact",
       x = "Time (days)", y = "FIQR", color = "Scenario") +
  theme_bw(14)

# Plot 3: SWS depth and fatigue
p3 <- ggplot(sim_all, aes(time_days, SWS_depth, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Slow-Wave Sleep Depth (normalized)",
       x = "Time (days)", y = "SWS Depth (a.u.)", color = "Scenario") +
  theme_bw(14)

# Plot 4: Central sensitization (LTP)
p4 <- ggplot(sim_all, aes(time_days, LTP_cs, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Central Sensitization Index (Spinal LTP)",
       x = "Time (days)", y = "LTP_cs (0-1)", color = "Scenario") +
  theme_bw(14)

# Plot 5: PK — duloxetine plasma concentration
p5 <- ggplot(filter(sim_DUL, time_days <= 14),
             aes(time_days, Cp_DUL_out)) +
  geom_line(color = "#3498DB", linewidth = 1) +
  labs(title = "Duloxetine PK — Plasma Concentration (first 14 days)",
       x = "Time (days)", y = "Cp (mg/L)") +
  theme_bw(14)

# Plot 6: PD biomarker — CSF Substance P
p6 <- ggplot(sim_all, aes(time_days, SP_csf, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "CSF Substance P (biomarker of central sensitization)",
       x = "Time (days)", y = "SP_csf (a.u.)", color = "Scenario") +
  theme_bw(14)

# Print summary at week 12
summary_wk12 <- sim_all %>%
  filter(time_days >= 83 & time_days <= 84) %>%
  group_by(scenario) %>%
  summarise(
    Pain_NRS     = round(mean(Pain_score), 2),
    FIQR         = round(mean(FIQ_score), 1),
    Fatigue_VAS  = round(mean(Fatigue_VAS), 1),
    Depression   = round(mean(Depression_score), 1),
    SWS_depth    = round(mean(SWS_depth), 3),
    LTP_cs       = round(mean(LTP_cs), 3),
    SP_csf       = round(mean(SP_csf), 3),
    DPMS         = round(mean(DPMS), 3),
    .groups = "drop"
  ) %>%
  arrange(Pain_NRS)

cat("\n=== FM QSP Model — 12-week Treatment Outcome Summary ===\n")
print(summary_wk12)

## ============================================================
## Responder analysis: ≥30% pain reduction
## ============================================================
baseline_pain <- mean(filter(sim_base, time_days <= 1)$Pain_score)

resp_30 <- sim_all %>%
  filter(time_days >= 83) %>%
  group_by(scenario) %>%
  summarise(
    pain_wk12 = mean(Pain_score),
    pct_change = 100 * (pain_wk12 - baseline_pain) / baseline_pain,
    responder_30 = pct_change <= -30,
    responder_50 = pct_change <= -50,
    .groups = "drop"
  )

cat("\n=== Responder Analysis ===\n")
print(resp_30)

## ============================================================
## Output plots
## ============================================================
library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange(p5, p6, nrow = 1)

message("FM QSP model simulation complete.")
