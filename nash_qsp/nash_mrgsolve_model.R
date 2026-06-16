## ============================================================
## NASH / NAFLD QSP Model — mrgsolve Implementation
## Disease: 비알코올성 지방간질환 (Non-Alcoholic Steatohepatitis)
## Version: 1.0  |  Date: 2026-06-16
## Author:  QSP Model Generator
##
## Compartments & Modules:
##   1. Drug PK  (FXR agonist / pan-PPAR agonist)
##   2. Hepatic lipid metabolism (FFA, TG, DNL, FAO)
##   3. Inflammation (TNF-α, IL-1β, IL-6, NF-κB)
##   4. Oxidative stress (ROS, JNK, ASK1)
##   5. Fibrosis (TGF-β1, HSC, Collagen)
##   6. Clinical endpoints (NAS, Fibrosis, ALT, Steatosis%)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ------------------------------------------------------------------
## 1. mrgsolve MODEL CODE BLOCK
## ------------------------------------------------------------------
model_code <- '

$PROB
NASH / NAFLD QSP Model
- FXR agonist (OCA-like) and pan-PPAR agonist (Lanifibranor-like) PK/PD
- Key modules: lipid, inflammation, oxidative stress, fibrosis
- Reference: NASH QSP framework (Vasilyeva et al. 2021, Diehl & Day 2017)

$PARAM @annotated
// ---- Drug PK Parameters (FXR agonist, OCA-like) ----
DOSE_FXR  : 10     : mg   // Daily oral dose FXR agonist
ka_FXR    : 1.2    : 1/h  // Absorption rate FXR agonist
CL_FXR    : 8.5    : L/h  // Clearance
Vc_FXR    : 30.0   : L    // Central volume
// ---- Drug PK Parameters (pan-PPAR agonist, Lanifibranor-like) ----
DOSE_PPAR : 800    : mg   // Daily oral dose pan-PPAR agonist (mg)
ka_PPAR   : 0.9    : 1/h  // Absorption rate
CL_PPAR   : 12.0   : L/h  // Clearance
Vc_PPAR   : 45.0   : L    // Central volume
// ---- FXR agonist PD ----
EC50_FXR  : 0.5    : mg/L // EC50 for FXR activation (OCA ~0.099 uM x MW 429)
Emax_FXR  : 1.0    : dml  // Max fractional FXR activation
// ---- pan-PPAR agonist PD ----
EC50_PPAR : 2.0    : mg/L // EC50 pan-PPAR activation
Emax_PPAR : 1.0    : dml  // Max fractional PPAR activation
// ---- Hepatic lipid parameters ----
kFFA_in   : 0.15   : 1/h  // Peripheral FFA flux into liver
kFFA_ox   : 0.08   : 1/h  // Basal FAO rate constant
kDNL      : 0.04   : 1/h  // Basal de novo lipogenesis rate
kTG_ester : 0.25   : 1/h  // FFA -> TG esterification
kTG_export: 0.06   : 1/h  // TG -> VLDL export
kTG_lipo  : 0.03   : 1/h  // Basal intrahepatic lipolysis
TG_base   : 2.5    : dml  // Baseline hepatic TG (relative units)
FFA_base  : 1.0    : dml  // Baseline hepatic FFA
// ---- Insulin resistance effect on lipid ----
IR_base   : 0.5    : dml  // Baseline insulin resistance index (0=none, 1=max)
kIR_DNL   : 1.8    : dml  // Fold-increase in DNL from IR
kIR_FAO   : 0.5    : dml  // Fold-decrease in FAO from IR
// ---- Inflammation parameters ----
k_NFkB    : 0.12   : 1/h  // NF-κB activation rate (from LPS/lipotoxicity)
k_NFkB_deg: 0.18   : 1/h  // NF-κB degradation
TNFa_base : 1.0    : dml  // Baseline TNF-α production (rel. units)
kTNFa_pr  : 0.25   : 1/h  // TNF-α production from NFkB
kTNFa_deg : 0.20   : 1/h  // TNF-α degradation
IL1b_base : 1.0    : dml  // IL-1β baseline
kIL1b_pr  : 0.20   : 1/h  // IL-1β production
kIL1b_deg : 0.22   : 1/h  // IL-1β degradation
// ---- Oxidative stress ----
kROS_gen  : 0.10   : 1/h  // ROS generation from mitochondrial FAO excess
kROS_scav : 0.15   : 1/h  // Antioxidant scavenging rate (Nrf2, GSH)
ROS_base  : 1.0    : dml  // Baseline ROS (rel. units)
kASK1_act : 0.08   : 1/h  // ASK1 activation from ROS
kASK1_deg : 0.12   : 1/h  // ASK1 inactivation
kJNK_act  : 0.15   : 1/h  // JNK activation from ASK1
kJNK_deg  : 0.18   : 1/h  // JNK inactivation
// ---- Fibrosis ----
kHSC_act  : 0.05   : 1/h  // HSC activation from TGFb+TNFa+JNK
kHSC_res  : 0.02   : 1/h  // HSC reversion to quiescence
kTGFb_pr  : 0.08   : 1/h  // TGF-β1 production from activated HSC
kTGFb_deg : 0.10   : 1/h  // TGF-β1 degradation
kCOL_pr   : 0.04   : 1/h  // Collagen production (from Smad2/3)
kCOL_deg  : 0.005  : 1/h  // Collagen degradation (MMP activity)
COL_base  : 1.0    : dml  // Baseline collagen index
// ---- Clinical endpoints ----
ALT_base  : 30     : U/L  // Baseline serum ALT
kALT_inj  : 0.80   : dml  // Scaling: hepatocyte injury -> ALT release
kALT_cl   : 0.04   : 1/h  // ALT clearance from plasma
// ---- PPAR effect on FAO ----
Emax_FAO_PPAR : 1.5 : dml // Max fold increase in FAO from pan-PPAR
EC50_FAO_PPAR : 1.2 : mg/L // EC50
// ---- FXR effect on inflammation ----
Emax_FXR_antiinf : 0.70 : dml // Max fractional suppression of NFkB by FXR
EC50_FXR_antiinf : 0.4  : mg/L // EC50 for anti-inflammatory effect

$CMT @annotated
// Drug PK
GUT_FXR  : Gut lumen FXR agonist (mg)
CENT_FXR : Central plasma FXR agonist (mg)
GUT_PPAR : Gut lumen pan-PPAR agonist (mg)
CENT_PPAR: Central plasma pan-PPAR agonist (mg)
// Lipid metabolism (relative to baseline = 1.0)
FFA_L    : Hepatic free fatty acid pool
TG_L     : Hepatic triglyceride pool
// Inflammation (relative)
NFkB_A   : Active NF-kB (nuclear)
TNFa     : TNF-alpha (liver)
IL1b     : IL-1 beta (liver)
// Oxidative stress
ROS      : Reactive oxygen species (relative)
ASK1_A   : Activated ASK1
JNK_A    : Activated JNK
// Fibrosis
HSC_A    : Activated hepatic stellate cells (fraction 0-1)
TGFb1    : TGF-beta 1 (liver)
COL      : Hepatic collagen content
// Clinical biomarkers
ALT_p    : Plasma ALT (U/L)

$MAIN
// ---- Drug concentrations (mg/L) ----
double Cp_FXR  = CENT_FXR  / Vc_FXR;
double Cp_PPAR = CENT_PPAR / Vc_PPAR;

// ---- PD occupancy (0 to 1) ----
double occ_FXR  = Emax_FXR  * Cp_FXR  / (EC50_FXR  + Cp_FXR);
double occ_PPAR = Emax_PPAR * Cp_PPAR / (EC50_PPAR + Cp_PPAR);

// ---- Lipotoxicity signal (drives inflammation & ROS) ----
double lip_tox = (FFA_L / FFA_base) * (TG_L / TG_base);

// ---- FAO rate (baseline + PPAR up-regulation) ----
double FAO_fold = 1.0 + Emax_FAO_PPAR * Cp_PPAR / (EC50_FAO_PPAR + Cp_PPAR);
double FAO_rate = kFFA_ox * FAO_fold;

// ---- DNL rate (reduced by FXR via SREBP1c suppression) ----
double FXR_DNL_supp = 1.0 - 0.60 * occ_FXR;
double DNL_rate = kDNL * (1.0 + kIR_DNL * IR_base) * FXR_DNL_supp;

// ---- NFkB production driven by LPS/lipotoxicity, suppressed by FXR ----
double FXR_NFkB_supp = 1.0 - Emax_FXR_antiinf * Cp_FXR / (EC50_FXR_antiinf + Cp_FXR);
double NFkB_prod = k_NFkB * lip_tox * FXR_NFkB_supp;

// ---- HSC activation driven by TGFb + TNFa + JNK ----
double HSC_drive = TGFb1 * (1.0 + 0.5 * TNFa) * (1.0 + 0.3 * JNK_A);

// ---- Collagen production, suppressed by PPARγ component of pan-PPAR ----
double PPAR_antifib = 1.0 - 0.50 * occ_PPAR;
double COL_prod = kCOL_pr * HSC_A * PPAR_antifib;

// ---- ALT injury signal ----
double inj_signal = lip_tox + JNK_A + 0.5 * TNFa;

$ODE
// ---- FXR agonist PK ----
dxdt_GUT_FXR   = -ka_FXR  * GUT_FXR;
dxdt_CENT_FXR  =  ka_FXR  * GUT_FXR  - (CL_FXR  / Vc_FXR)  * CENT_FXR;

// ---- pan-PPAR agonist PK ----
dxdt_GUT_PPAR  = -ka_PPAR * GUT_PPAR;
dxdt_CENT_PPAR =  ka_PPAR * GUT_PPAR - (CL_PPAR / Vc_PPAR) * CENT_PPAR;

// ---- Hepatic FFA ----
// Input: peripheral FFA influx + DNL; Output: FAO + esterification to TG
dxdt_FFA_L = kFFA_in * FFA_base
           + DNL_rate
           - FAO_rate  * FFA_L
           - kTG_ester * FFA_L;

// ---- Hepatic TG ----
// Input: esterification from FFA; Output: VLDL export + intrahepatic lipolysis
dxdt_TG_L  = kTG_ester * FFA_L
           - kTG_export * TG_L
           - kTG_lipo   * TG_L;

// ---- NF-κB ----
dxdt_NFkB_A = NFkB_prod - k_NFkB_deg * NFkB_A;

// ---- TNF-α ----
dxdt_TNFa   = kTNFa_pr * NFkB_A - kTNFa_deg * TNFa;

// ---- IL-1β (also driven by NLRP3 / Caspase-1 from lipotoxicity) ----
dxdt_IL1b   = kIL1b_pr * NFkB_A * (1.0 + 0.4 * lip_tox) - kIL1b_deg * IL1b;

// ---- ROS: generated by lipotoxicity and FAO overload; scavenged by Nrf2/GSH ----
// NB: Nrf2 upregulated by PPAR agonist (mild)
double Nrf2_fold = 1.0 + 0.3 * occ_PPAR;
dxdt_ROS    = kROS_gen * lip_tox - kROS_scav * Nrf2_fold * ROS;

// ---- ASK1 activation by ROS ----
dxdt_ASK1_A = kASK1_act * ROS - kASK1_deg * ASK1_A;

// ---- JNK activation by ASK1 ----
dxdt_JNK_A  = kJNK_act  * ASK1_A - kJNK_deg * JNK_A;

// ---- HSC activation ----
dxdt_HSC_A  = kHSC_act * HSC_drive * (1.0 - HSC_A) - kHSC_res * HSC_A;

// ---- TGF-β1: autocrine from activated HSC ----
dxdt_TGFb1  = kTGFb_pr * HSC_A - kTGFb_deg * TGFb1;

// ---- Collagen ----
dxdt_COL    = COL_prod - kCOL_deg * COL;

// ---- Plasma ALT ----
dxdt_ALT_p  = kALT_inj * inj_signal - kALT_cl * ALT_p;

$TABLE
// ---- Derived clinical endpoints ----
// Steatosis % (0-100): sigmoid function of hepatic TG
double Steatosis_pct = 100.0 * (TG_L / TG_base) / (1.0 + TG_L / TG_base) * 0.66;

// NAS component scores (0-3 steatosis, 0-3 inflammation, 0-2 ballooning)
double NAS_steat  = 3.0 * Steatosis_pct / 100.0;
double NAS_inflam = 3.0 * (TNFa + IL1b) / 2.0 / (1.0 + (TNFa + IL1b) / 2.0);
double NAS_bal    = 2.0 * JNK_A / (1.0 + JNK_A);
double NAS_total  = NAS_steat + NAS_inflam + NAS_bal;

// Fibrosis stage proxy (0-4 Metavir): sigmoid of collagen index
double Fibrosis_stage = 4.0 * (COL - 1.0) / (1.0 + (COL - 1.0)) ;
double Fibrosis_stage2 = (Fibrosis_stage < 0.0) ? 0.0 : Fibrosis_stage;

// FIB-4 proxy (simplified)
double FIB4_proxy = ALT_p / 30.0 * (COL / COL_base);

// Drug concentrations
double Cp_FXR_out  = CENT_FXR  / Vc_FXR;
double Cp_PPAR_out = CENT_PPAR / Vc_PPAR;

// PD occupancy
double Eff_FXR  = Emax_FXR  * Cp_FXR_out  / (EC50_FXR  + Cp_FXR_out);
double Eff_PPAR = Emax_PPAR * Cp_PPAR_out / (EC50_PPAR + Cp_PPAR_out);

$CAPTURE
Cp_FXR_out Cp_PPAR_out Eff_FXR Eff_PPAR
Steatosis_pct NAS_total NAS_steat NAS_inflam NAS_bal
Fibrosis_stage2 FIB4_proxy
TG_L FFA_L NFkB_A TNFa IL1b ROS ASK1_A JNK_A HSC_A TGFb1 COL ALT_p
'

## ------------------------------------------------------------------
## 2. COMPILE MODEL
## ------------------------------------------------------------------
nash_mod <- mcode("nash_qsp", model_code)

## ------------------------------------------------------------------
## 3. INITIAL CONDITIONS (homeostatic baseline)
## ------------------------------------------------------------------
init_baseline <- c(
  GUT_FXR   = 0,
  CENT_FXR  = 0,
  GUT_PPAR  = 0,
  CENT_PPAR = 0,
  FFA_L     = 1.0,  # normalized baseline
  TG_L      = 1.0,  # normalized baseline (mild NAFL = ~2.5)
  NFkB_A    = 0.5,
  TNFa      = 0.5,
  IL1b      = 0.5,
  ROS       = 1.0,
  ASK1_A    = 0.5,
  JNK_A     = 0.5,
  HSC_A     = 0.1,
  TGFb1     = 0.2,
  COL       = 1.0,
  ALT_p     = 30.0
)

## ------------------------------------------------------------------
## 4. DISEASE PROGRESSION SCENARIO (no treatment, 52 weeks)
## ------------------------------------------------------------------
run_disease_progression <- function(
  IR_level = 0.7,   # 0 = no IR, 1 = max IR
  NASH_init = TRUE  # start with established NAFL
) {
  params <- param(nash_mod, IR_base = IR_level)

  init_vals <- init_baseline
  if (NASH_init) {
    init_vals["TG_L"]  <- 2.5  # moderate steatosis
    init_vals["NFkB_A"] <- 1.2
    init_vals["TNFa"]  <- 1.0
    init_vals["IL1b"]  <- 1.0
    init_vals["ROS"]   <- 1.5
    init_vals["HSC_A"] <- 0.25
    init_vals["COL"]   <- 1.5
  }

  ev   <- ev(time = 0, amt = 0, cmt = 1)  # no treatment
  sims <- params %>%
    init(init_vals) %>%
    mrgsim(events = ev, end = 52 * 7 * 24, delta = 24) %>%
    as_tibble()

  sims
}

## ------------------------------------------------------------------
## 5. TREATMENT SCENARIO: FXR agonist (OCA-like), daily oral
## ------------------------------------------------------------------
run_FXR_treatment <- function(
  dose_fxr  = 10,   # mg once daily
  duration_wk = 52,
  IR_level  = 0.7
) {
  params <- param(nash_mod, IR_base = IR_level, DOSE_FXR = dose_fxr)

  init_vals <- init_baseline
  init_vals["TG_L"]  <- 2.5
  init_vals["NFkB_A"] <- 1.2
  init_vals["TNFa"]  <- 1.0
  init_vals["IL1b"]  <- 1.0
  init_vals["ROS"]   <- 1.5
  init_vals["ASK1_A"] <- 0.7
  init_vals["JNK_A"] <- 0.6
  init_vals["HSC_A"] <- 0.25
  init_vals["COL"]   <- 1.5
  init_vals["ALT_p"] <- 60

  # Once-daily dosing events
  n_doses <- duration_wk * 7
  dose_ev <- ev(amt  = dose_fxr,
                cmt  = 1,          # GUT_FXR
                ii   = 24,         # 24h interval
                addl = n_doses - 1)

  sims <- params %>%
    init(init_vals) %>%
    mrgsim(events = dose_ev, end = n_doses * 24, delta = 12) %>%
    as_tibble()

  sims
}

## ------------------------------------------------------------------
## 6. TREATMENT SCENARIO: pan-PPAR agonist (Lanifibranor-like)
## ------------------------------------------------------------------
run_PPAR_treatment <- function(
  dose_ppar = 800,   # mg once daily
  duration_wk = 52,
  IR_level  = 0.7
) {
  params <- param(nash_mod, IR_base = IR_level, DOSE_PPAR = dose_ppar)

  init_vals <- init_baseline
  init_vals["TG_L"]  <- 2.5
  init_vals["NFkB_A"] <- 1.2
  init_vals["TNFa"]  <- 1.0
  init_vals["IL1b"]  <- 1.0
  init_vals["ROS"]   <- 1.5
  init_vals["ASK1_A"] <- 0.7
  init_vals["JNK_A"] <- 0.6
  init_vals["HSC_A"] <- 0.25
  init_vals["COL"]   <- 1.5
  init_vals["ALT_p"] <- 60

  n_doses <- duration_wk * 7
  dose_ev <- ev(amt  = dose_ppar,
                cmt  = 3,          # GUT_PPAR
                ii   = 24,
                addl = n_doses - 1)

  sims <- params %>%
    init(init_vals) %>%
    mrgsim(events = dose_ev, end = n_doses * 24, delta = 12) %>%
    as_tibble()

  sims
}

## ------------------------------------------------------------------
## 7. COMBINATION THERAPY: FXR + pan-PPAR
## ------------------------------------------------------------------
run_combo_treatment <- function(
  dose_fxr  = 10,
  dose_ppar = 800,
  duration_wk = 52,
  IR_level  = 0.7
) {
  params <- param(nash_mod, IR_base = IR_level,
                  DOSE_FXR = dose_fxr, DOSE_PPAR = dose_ppar)

  init_vals <- init_baseline
  init_vals["TG_L"]  <- 2.5
  init_vals["NFkB_A"] <- 1.2
  init_vals["TNFa"]  <- 1.0
  init_vals["IL1b"]  <- 1.0
  init_vals["ROS"]   <- 1.5
  init_vals["ASK1_A"] <- 0.7
  init_vals["JNK_A"] <- 0.6
  init_vals["HSC_A"] <- 0.25
  init_vals["COL"]   <- 1.5
  init_vals["ALT_p"] <- 60

  n_doses <- duration_wk * 7
  fxr_ev  <- ev(amt  = dose_fxr,  cmt = 1, ii = 24, addl = n_doses - 1)
  ppar_ev <- ev(amt  = dose_ppar, cmt = 3, ii = 24, addl = n_doses - 1)
  combo_ev <- c(fxr_ev, ppar_ev)

  sims <- params %>%
    init(init_vals) %>%
    mrgsim(events = combo_ev, end = n_doses * 24, delta = 12) %>%
    as_tibble()

  sims
}

## ------------------------------------------------------------------
## 8. PARAMETER SENSITIVITY ANALYSIS (one-at-a-time)
## ------------------------------------------------------------------
run_sensitivity <- function(param_name, fold_range = c(0.5, 2.0), n_steps = 5) {
  base_val  <- param(nash_mod)[[param_name]]
  fold_vals <- seq(fold_range[1], fold_range[2], length.out = n_steps)

  results <- lapply(fold_vals, function(f) {
    new_val <- base_val * f
    p <- param(nash_mod)
    p[[param_name]] <- new_val

    init_vals <- init_baseline
    init_vals["TG_L"] <- 2.5; init_vals["COL"] <- 1.5

    sim <- p %>%
      init(init_vals) %>%
      mrgsim(end = 52 * 7 * 24, delta = 24) %>%
      as_tibble() %>%
      filter(time == max(time)) %>%
      mutate(fold = f, param = param_name, param_val = new_val)
    sim
  })

  bind_rows(results)
}

## ------------------------------------------------------------------
## 9. RUN ALL SCENARIOS
## ------------------------------------------------------------------
message("Running NASH QSP simulations...")

sim_notreat <- run_disease_progression()  %>% mutate(arm = "No Treatment")
sim_fxr     <- run_FXR_treatment()        %>% mutate(arm = "FXR Agonist (OCA-like)")
sim_ppar    <- run_PPAR_treatment()       %>% mutate(arm = "pan-PPAR (Lanifibranor-like)")
sim_combo   <- run_combo_treatment()      %>% mutate(arm = "Combination")

all_sims <- bind_rows(sim_notreat, sim_fxr, sim_ppar, sim_combo) %>%
  mutate(time_wk = time / (7 * 24))

## ------------------------------------------------------------------
## 10. PLOTS
## ------------------------------------------------------------------
theme_qsp <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#E3F2FD"),
    legend.position  = "bottom",
    legend.title     = element_blank()
  )

arm_colors <- c(
  "No Treatment"               = "#B71C1C",
  "FXR Agonist (OCA-like)"     = "#1565C0",
  "pan-PPAR (Lanifibranor-like)" = "#2E7D32",
  "Combination"                = "#6A1B9A"
)

## -- NAS score over time
p1 <- ggplot(all_sims, aes(time_wk, NAS_total, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "NAS Score (0–8)", x = "Time (weeks)", y = "NAS Total") +
  theme_qsp

## -- Fibrosis stage over time
p2 <- ggplot(all_sims, aes(time_wk, Fibrosis_stage2, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "Fibrosis Stage (0–4 Metavir proxy)", x = "Time (weeks)", y = "Fibrosis Stage") +
  theme_qsp

## -- Hepatic TG (steatosis)
p3 <- ggplot(all_sims, aes(time_wk, TG_L, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "Hepatic TG (normalized)", x = "Time (weeks)", y = "TG_L (rel. to baseline)") +
  theme_qsp

## -- ALT
p4 <- ggplot(all_sims, aes(time_wk, ALT_p, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "Plasma ALT (U/L)", x = "Time (weeks)", y = "ALT (U/L)") +
  theme_qsp

## -- Inflammation: TNF-α
p5 <- ggplot(all_sims, aes(time_wk, TNFa, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "Hepatic TNF-α", x = "Time (weeks)", y = "TNF-α (relative)") +
  theme_qsp

## -- Collagen (fibrosis marker)
p6 <- ggplot(all_sims, aes(time_wk, COL, color = arm)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = arm_colors) +
  labs(title = "Hepatic Collagen Index", x = "Time (weeks)", y = "Collagen (relative)") +
  theme_qsp

combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title   = "NASH QSP Model — Treatment Comparison (52 weeks)",
    subtitle = "비알코올성 지방간염 QSP 시뮬레이션",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(combined_plot)
ggsave("nash_simulation_results.pdf", combined_plot, width = 14, height = 12)

## ------------------------------------------------------------------
## 11. SENSITIVITY ANALYSIS PLOT
## ------------------------------------------------------------------
sens_params <- c("kDNL", "kFFA_ox", "kHSC_act", "kTGFb_pr", "IR_base")
sens_results <- lapply(sens_params, run_sensitivity) %>% bind_rows()

p_sens <- ggplot(sens_results, aes(x = fold, y = Fibrosis_stage2, color = param)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title    = "One-at-a-Time Sensitivity: Final Fibrosis Stage (52 wk)",
    subtitle = "Parameter fold-change vs. No-Treatment baseline",
    x        = "Parameter fold-change from baseline",
    y        = "Fibrosis Stage at 52 wk (proxy)"
  ) +
  theme_qsp

print(p_sens)
ggsave("nash_sensitivity_analysis.pdf", p_sens, width = 9, height = 6)

message("Done. Outputs: nash_simulation_results.pdf, nash_sensitivity_analysis.pdf")
