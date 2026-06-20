## =============================================================================
## Atopic Dermatitis (AD) — QSP Model (mrgsolve)
## =============================================================================
## Mechanistic scope:
##   • Dupilumab 2-compartment PK (SC Q2W) — anti-IL-4Rα mAb
##   • Upadacitinib 1-compartment PK (oral QD) — JAK1 inhibitor
##   • Nemolizumab 1-compartment PK (SC Q4W) — anti-IL-31Ra mAb
##   • IL-4Rα receptor occupancy (RO) with target-mediated drug disposition (TMDD)
##   • pSTAT6 inhibition → downstream cytokine/biomarker dynamics
##   • Skin compartment: Th2 inflammation, ILC2, eosinophils, barrier function
##   • IL-31 / pruritus NRS dynamics
##   • Clinical endpoints: EASI, IGA, NRS itch
## References:
##   Beck et al. NEJM 2014 (Dupilumab Phase 2)
##   Simpson et al. NEJM 2016 (SOLO-1/2)
##   Reich et al. NEJM 2020 (ECZTRA-1/2)
##   Guttman-Yassky et al. Lancet 2020 (Rising Up — Upadacitinib)
##   Silverberg et al. NEJM 2021 (JADE MONO-1/2 — Abrocitinib)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## =============================================================================
## Model Code
## =============================================================================

ad_model_code <- '
$PROB Atopic Dermatitis QSP Model
  Dupilumab + Upadacitinib + Nemolizumab PK/PD
  Skin inflammation, barrier, pruritus, EASI dynamics

$PARAM @annotated
  // --- Dupilumab PK (2-compartment, SC) ---
  KA_DUP    : 0.12    : /day    : SC absorption rate constant
  CL_DUP    : 0.21    : L/day   : Dupilumab clearance
  V1_DUP    : 3.5     : L       : Central volume
  V2_DUP    : 4.2     : L       : Peripheral volume
  Q_DUP     : 0.8     : L/day   : Inter-compartment clearance
  F_SC_DUP  : 0.64    :         : SC bioavailability
  MW_DUP    : 148000  : g/mol   : Dupilumab MW (144kDa mAb)

  // --- Dupilumab PD: IL-4Rα Receptor Occupancy (RO) ---
  Kd_IL4Ra  : 0.027   : nM      : Dupilumab-IL4Rα KD
  Rtot_IL4Ra: 2.5     : nM      : Total IL-4Rα (skin compartment equiv.)
  kon_IL4Ra : 0.26    : /nM/day : Association rate
  koff_IL4Ra: 0.007   : /day    : Dissociation rate

  // --- STAT6 phosphorylation dynamics ---
  kSTAT6_in  : 0.5    : /day    : Basal pSTAT6 induction rate
  kSTAT6_out : 0.8    : /day    : pSTAT6 dephosphorylation
  EC50_STAT6 : 0.45   :         : Dupilumab RO for 50% STAT6 inhibition
  Emax_STAT6 : 0.95   :         : Max STAT6 inhibition fraction

  // --- JAK inhibitor (Upadacitinib) PK ---
  KA_UPA    : 3.2     : /day    : Oral absorption rate
  CL_UPA    : 38.0    : L/day   : Apparent clearance
  Vd_UPA    : 166     : L       : Apparent volume
  F_UPA     : 0.79    :         : Oral bioavailability

  // --- Upadacitinib PD: JAK1 inhibition ---
  IC50_JAK1  : 0.045  : ng/mL   : Upadacitinib JAK1 IC50
  Imax_JAK1  : 0.97   :         : Max JAK1 inhibition

  // --- Nemolizumab PK (anti-IL-31Ra, SC Q4W) ---
  KA_NEMO   : 0.09    : /day    : SC absorption rate
  CL_NEMO   : 0.15    : L/day   : Clearance
  V1_NEMO   : 3.8     : L       : Volume
  F_SC_NEMO : 0.72    :         : SC bioavailability

  // --- IL-31 / Pruritus dynamics ---
  IL31_base  : 45.0   : pg/mL   : Baseline IL-31 (AD patients elevated)
  kIL31_prod : 0.18   : /day    : IL-31 production by Th2/ILC2
  kIL31_deg  : 0.22   : /day    : IL-31 degradation
  EC50_NEMO  : 12.0   : pg/mL   : IL-31 EC50 for pruritus
  Emax_itch  : 0.85   :         : Max itch reduction from IL-31 block
  IC50_itch_JAK : 0.08: ng/mL   : Upadacitinib IC50 for itch (JAK1)

  // --- Th2 inflammation dynamics (skin) ---
  Th2_base   : 100.0  :         : Baseline Th2 cell index (AU)
  kTh2_in    : 8.5    : /day    : Th2 recruitment rate
  kTh2_out   : 0.09   : /day    : Th2 elimination rate
  TARC_base  : 3500   : pg/mL   : Baseline TARC/CCL17 in AD
  kTARC_in   : 2.0    : /day    : TARC production
  kTARC_out  : 0.35   : /day    : TARC clearance

  // --- ILC2 activation ---
  ILC2_base  : 50.0   :         : Baseline ILC2 activity index
  kILC2_in   : 3.0    : /day    : TSLP/IL33 driven ILC2 activation
  kILC2_out  : 0.12   : /day    : ILC2 deactivation

  // --- Eosinophil dynamics ---
  Eos_base   : 350.0  : cells/uL: Baseline blood eosinophil count
  kEos_in    : 25.0   : cells/uL/day : IL-5 driven eosinophil production
  kEos_out   : 0.07   : /day    : Eosinophil clearance
  kEos_tiss  : 0.03   : /day    : Tissue recruitment rate

  // --- IgE dynamics ---
  IgE_base   : 2000   : IU/mL   : Baseline total IgE (AD patients)
  kIgE_in    : 0.05   : /day    : IL-4 driven IgE production
  kIgE_out   : 0.005  : /day    : IgE catabolism (slow, t1/2~70d)

  // --- Epidermal barrier (Filaggrin / TEWL) ---
  FLG_base   : 0.40   :         : Baseline relative FLG expression (normal=1.0)
  kFLG_rest  : 0.02   : /day    : FLG restoration rate (treatment effect)
  kFLG_sup   : 0.04   : /day    : STAT6-driven FLG suppression rate
  TEWL_base  : 25.0   : g/m2/h  : Baseline TEWL in AD
  TEWL_norm  : 8.0    : g/m2/h  : Normal skin TEWL

  // --- Topical corticosteroid (TCS) additive effect ---
  TCS_effect : 0.0    :         : TCS effect on STAT6 (0=no TCS, 0.4=TCS)

  // --- Skin inflammation composite (maps to EASI) ---
  kInfl_prod : 0.8    : /day    : Inflammatory signal production
  kInfl_elim : 0.15   : /day    : Inflammatory signal resolution

  // --- EASI score dynamics ---
  EASI_base  : 29.0   :         : Baseline EASI (severe=29-72, mod=7-28)
  kEASI_resp : 0.10   : /day    : EASI response rate constant
  EASI_min   : 0.5    :         : Minimum possible EASI

  // --- IGA score ---
  IGA_thresh : 4.0    :         : IGA 0/1 threshold EASI

  // --- Simulation flags ---
  DUP_on     : 0      :         : Dupilumab dosing flag
  UPA_on     : 0      :         : Upadacitinib dosing flag
  NEMO_on    : 0      :         : Nemolizumab dosing flag

$CMT @annotated
  // Dupilumab PK
  DUP_SC     : Dupilumab SC depot (mg)
  DUP_C      : Dupilumab central compartment (mg)
  DUP_P      : Dupilumab peripheral compartment (mg)

  // IL-4Rα receptor
  Rfree      : Free IL-4Rα receptor (nM)
  RC         : Receptor-drug complex (nM)

  // pSTAT6
  pSTAT6     : Phosphorylated STAT6 (AU)

  // Upadacitinib PK
  UPA_GI     : Upadacitinib GI absorption (mg)
  UPA_C      : Upadacitinib central (mg)

  // Nemolizumab PK
  NEMO_SC    : Nemolizumab SC depot (mg)
  NEMO_C     : Nemolizumab central (mg)

  // IL-31
  IL31       : IL-31 plasma/skin level (pg/mL)

  // Th2 inflammation
  Th2        : Th2 cell index (AU)
  TARC       : TARC/CCL17 level (pg/mL)
  ILC2       : ILC2 activation index (AU)

  // Eosinophil
  Eos_blood  : Blood eosinophil count (cells/uL)

  // IgE
  IgE        : Total IgE level (IU/mL)

  // Barrier
  FLG        : Filaggrin expression (relative)
  TEWL       : Trans-epidermal water loss (g/m2/h)

  // Skin inflammation composite
  SkinInfl   : Skin inflammation index (AU)

  // EASI score
  EASI       : Eczema Area & Severity Index

$INIT
  DUP_SC  = 0
  DUP_C   = 0
  DUP_P   = 0
  Rfree   = 2.5
  RC      = 0
  pSTAT6  = 100
  UPA_GI  = 0
  UPA_C   = 0
  NEMO_SC = 0
  NEMO_C  = 0
  IL31    = 45.0
  Th2     = 100.0
  TARC    = 3500.0
  ILC2    = 50.0
  Eos_blood = 350.0
  IgE     = 2000.0
  FLG     = 0.40
  TEWL    = 25.0
  SkinInfl = 100.0
  EASI    = 29.0

$ODE
  // ---------------------------------------------------------------
  // Dupilumab 2-compartment PK (SC administration)
  // ---------------------------------------------------------------
  double dose_rate_dup = KA_DUP * DUP_SC * F_SC_DUP;
  dxdt_DUP_SC  = -KA_DUP * DUP_SC;
  double Cp_dup = DUP_C / V1_DUP;   // mg/L = ug/mL
  double Cp2_dup = DUP_P / V2_DUP;
  dxdt_DUP_C   = dose_rate_dup
                 - (CL_DUP/V1_DUP) * DUP_C
                 - (Q_DUP/V1_DUP)  * DUP_C
                 + (Q_DUP/V2_DUP)  * DUP_P;
  dxdt_DUP_P   = (Q_DUP/V1_DUP) * DUP_C - (Q_DUP/V2_DUP) * DUP_P;

  // Convert Cp to nM: MW_DUP=148000 g/mol; 1 ug/mL = 1e-3 g/L
  double Cp_dup_nM = (Cp_dup * 1000.0) / MW_DUP * 1e6;

  // ---------------------------------------------------------------
  // IL-4Rα Receptor Occupancy (simplified quasi-steady-state TMDD)
  // ---------------------------------------------------------------
  double Rfree_eq = Rtot_IL4Ra * Kd_IL4Ra / (Kd_IL4Ra + Cp_dup_nM);
  // Dynamic RO model
  dxdt_Rfree = koff_IL4Ra * RC - kon_IL4Ra * Rfree * Cp_dup_nM;
  dxdt_RC    = kon_IL4Ra * Rfree * Cp_dup_nM - koff_IL4Ra * RC;
  double RO  = RC / (Rtot_IL4Ra + 1e-10);  // fraction occupied 0-1

  // ---------------------------------------------------------------
  // pSTAT6 dynamics (suppressed by RO and JAK1 inhibition)
  // ---------------------------------------------------------------
  // JAK1 inhibition by upadacitinib
  double Cp_upa  = UPA_C / Vd_UPA;  // mg/L = ug/mL -> ng/mL * 1000
  double Cp_upa_ngmL = Cp_upa * 1000.0;
  double JAK1_inh = Imax_JAK1 * Cp_upa_ngmL / (IC50_JAK1 + Cp_upa_ngmL + 1e-10);
  // Combined STAT6 inhibition
  double STAT6_inh_dup = Emax_STAT6 * RO / (EC50_STAT6 + RO + 1e-10);
  double STAT6_inh_total = 1.0 - (1.0 - STAT6_inh_dup) * (1.0 - JAK1_inh) * (1.0 - TCS_effect);
  if (STAT6_inh_total > 0.99) STAT6_inh_total = 0.99;
  dxdt_pSTAT6 = kSTAT6_in * (1.0 - STAT6_inh_total) * 100.0 - kSTAT6_out * pSTAT6;

  // ---------------------------------------------------------------
  // Upadacitinib 1-compartment PK
  // ---------------------------------------------------------------
  dxdt_UPA_GI = -KA_UPA * UPA_GI;
  dxdt_UPA_C  = KA_UPA * F_UPA * UPA_GI - (CL_UPA / Vd_UPA) * UPA_C;

  // ---------------------------------------------------------------
  // Nemolizumab PK (anti-IL-31Ra)
  // ---------------------------------------------------------------
  dxdt_NEMO_SC = -KA_NEMO * NEMO_SC;
  double Cp_nemo = NEMO_C / V1_NEMO;  // mg/L
  dxdt_NEMO_C  = KA_NEMO * F_SC_NEMO * NEMO_SC
                 - (CL_NEMO / V1_NEMO) * NEMO_C;
  // Nemolizumab IL-31 receptor occupancy (simplified)
  double RO_NEMO = Cp_nemo / (Cp_nemo + 0.015);  // KD~0.015 mg/L for nemolizumab

  // ---------------------------------------------------------------
  // IL-31 dynamics
  // ---------------------------------------------------------------
  double Th2_drive = Th2 / 100.0;
  double ILC2_drive = ILC2 / 50.0;
  double IL31_prod = kIL31_prod * (Th2_drive * 0.7 + ILC2_drive * 0.3) * 45.0;
  double IL31_block_frac = RO_NEMO * 0.85 + JAK1_inh * 0.60;  // JAKi also reduces IL-31 downstream
  if (IL31_block_frac > 0.95) IL31_block_frac = 0.95;
  dxdt_IL31 = IL31_prod * (1.0 - IL31_block_frac) - kIL31_deg * IL31;

  // ---------------------------------------------------------------
  // ILC2 activation (driven by TSLP/IL-33; suppressed by JAKi)
  // ---------------------------------------------------------------
  double ILC2_drive_in = kILC2_in * (1.0 - JAK1_inh * 0.5);
  dxdt_ILC2 = ILC2_drive_in * 50.0 / (1.0 + ILC2/50.0) - kILC2_out * ILC2;

  // ---------------------------------------------------------------
  // Th2 dynamics (skin)
  // ---------------------------------------------------------------
  double TARC_drive = TARC / 3500.0;
  double Th2_in = kTh2_in * TARC_drive * (1.0 - STAT6_inh_total * 0.7);
  dxdt_Th2 = Th2_in - kTh2_out * Th2;

  // ---------------------------------------------------------------
  // TARC/CCL17 dynamics
  // ---------------------------------------------------------------
  double pSTAT6_norm = pSTAT6 / 100.0;
  dxdt_TARC = kTARC_in * pSTAT6_norm * 3500.0 - kTARC_out * TARC;

  // ---------------------------------------------------------------
  // Blood eosinophil count
  // ---------------------------------------------------------------
  double IL5_equiv = Th2 / 100.0 * ILC2 / 50.0;  // proxy for IL-5
  double Eos_in = kEos_in * IL5_equiv * (1.0 - STAT6_inh_total * 0.8);
  dxdt_Eos_blood = Eos_in - kEos_out * Eos_blood;

  // ---------------------------------------------------------------
  // IgE dynamics (slow — months to years)
  // ---------------------------------------------------------------
  double IgE_prod = kIgE_in * (Th2/100.0) * (1.0 - STAT6_inh_total * 0.6) * 2000.0;
  dxdt_IgE = IgE_prod - kIgE_out * IgE;

  // ---------------------------------------------------------------
  // Filaggrin expression and TEWL
  // ---------------------------------------------------------------
  double FLG_sup_STAT6 = kFLG_sup * pSTAT6_norm;  // STAT6 suppresses FLG
  double FLG_rest = kFLG_rest * (1.0 - FLG);       // Restoration toward 1.0
  dxdt_FLG  = FLG_rest - FLG_sup_STAT6 * FLG;
  // TEWL inverse of FLG
  double TEWL_target = TEWL_norm + (TEWL_base - TEWL_norm) * (1.0 - FLG) / 0.60;
  dxdt_TEWL = 0.15 * (TEWL_target - TEWL);

  // ---------------------------------------------------------------
  // Skin inflammation composite index
  // ---------------------------------------------------------------
  double Infl_driver = (Th2/100.0) * (pSTAT6/100.0) * (Eos_blood/350.0);
  dxdt_SkinInfl = kInfl_prod * Infl_driver * 100.0 - kInfl_elim * SkinInfl;

  // ---------------------------------------------------------------
  // EASI score dynamics
  // ---------------------------------------------------------------
  double EASI_target = EASI_min + (EASI_base - EASI_min) * (SkinInfl / 100.0);
  dxdt_EASI = kEASI_resp * (EASI_target - EASI);

$TABLE
  // Derived PK/PD outputs
  double Cp_dup_ugmL = DUP_C / V1_DUP;
  double RO_pct      = RC / (Rtot_IL4Ra + 1e-10) * 100.0;
  double pSTAT6_pct  = pSTAT6;
  double JAK1inh_pct = (Imax_JAK1 * (UPA_C/Vd_UPA*1000.0) /
                        (IC50_JAK1 + UPA_C/Vd_UPA*1000.0 + 1e-10)) * 100.0;

  // Pruritus NRS (0-10): driven by IL-31 and JAK1 inhibition
  double itch_IL31   = 10.0 * (IL31 / IL31_base) * (1.0 - RO_NEMO * 0.85);
  double itch_JAKblock = JAK1inh_pct / 100.0 * 0.7;
  double NRS_itch    = itch_IL31 * (1.0 - itch_JAKblock);
  if (NRS_itch < 0.2) NRS_itch = 0.2;
  if (NRS_itch > 10.0) NRS_itch = 10.0;

  // IGA: 0=clear, 1=almost clear, 2=mild, 3=moderate, 4=severe
  double IGA = 0.0;
  if (EASI > 21.0)      IGA = 4.0;
  else if (EASI > 14.0) IGA = 3.0;
  else if (EASI > 7.0)  IGA = 2.0;
  else if (EASI > 2.0)  IGA = 1.0;
  else                  IGA = 0.0;

  // EASI response endpoints
  double EASI75_achieved = (EASI <= EASI_base * 0.25) ? 1.0 : 0.0;
  double EASI90_achieved = (EASI <= EASI_base * 0.10) ? 1.0 : 0.0;

  // SCORAD (proxy): linear map from EASI
  double SCORAD_est = EASI * 2.1 + NRS_itch * 1.5;

  // Blood biomarkers
  double TARC_pct_change = (TARC - TARC_base) / TARC_base * 100.0;
  double Eos_pct_change  = (Eos_blood - Eos_base) / Eos_base * 100.0;

  // Capture variables
  capture Cp_dup_ugmL, RO_pct, pSTAT6_pct, JAK1inh_pct,
          NRS_itch, IGA, EASI75_achieved, EASI90_achieved,
          TARC, SCORAD_est, Eos_blood,
          TARC_pct_change, Eos_pct_change,
          FLG, TEWL, IgE, IL31

$CAPTURE EASI NRS_itch IGA RO_pct pSTAT6_pct Cp_dup_ugmL JAK1inh_pct
'

ad_mod <- mcode("atopic_dermatitis_qsp", ad_model_code)

## =============================================================================
## Dosing Functions
## =============================================================================

# Dupilumab: 600mg loading dose, then 300mg Q2W SC
make_dup_dosing <- function(n_doses = 28, start_day = 0) {
  ev_load <- ev(amt = 600, cmt = "DUP_SC", time = start_day)
  maint_times <- seq(start_day + 14, by = 14, length.out = n_doses - 1)
  ev_maint <- ev(amt = 300, cmt = "DUP_SC", time = maint_times)
  c(ev_load, ev_maint)
}

# Upadacitinib: 30mg QD oral
make_upa_dosing <- function(n_days = 360, dose_mg = 30, start_day = 0) {
  ev(amt = dose_mg, cmt = "UPA_GI", time = seq(start_day, by = 1, length.out = n_days))
}

# Nemolizumab: 60mg SC Q4W (loading 30mg at day 0)
make_nemo_dosing <- function(n_doses = 13, start_day = 0) {
  ev(amt = 60, cmt = "NEMO_SC", time = seq(start_day, by = 28, length.out = n_doses))
}

## =============================================================================
## Simulation Setup: 5 Treatment Scenarios (52 weeks)
## =============================================================================

sim_time <- seq(0, 364, by = 1)   # 52 weeks daily output

run_scenario <- function(scenario_name, dup_on, upa_on, nemo_on,
                         tcs_effect = 0.0, upa_dose = 30) {
  params <- c(DUP_on    = dup_on,
              UPA_on    = upa_on,
              NEMO_on   = nemo_on,
              TCS_effect = tcs_effect)

  events <- NULL
  if (dup_on  == 1) events <- c(events, make_dup_dosing())
  if (upa_on  == 1) events <- c(events, make_upa_dosing(dose_mg = upa_dose))
  if (nemo_on == 1) events <- c(events, make_nemo_dosing())

  if (is.null(events)) {
    out <- ad_mod %>%
      param(params) %>%
      mrgsim(delta = 1, end = 364)
  } else {
    out <- ad_mod %>%
      param(params) %>%
      ev(events) %>%
      mrgsim(delta = 1, end = 364)
  }

  as_tibble(out) %>%
    mutate(scenario = scenario_name)
}

# Run all 5 scenarios
scenarios <- list(
  # 1. No treatment (disease progression)
  list(name = "1. No Treatment",
       dup = 0, upa = 0, nemo = 0, tcs = 0.0),

  # 2. Topical corticosteroid (TCS) only — represented by TCS_effect
  list(name = "2. TCS Only",
       dup = 0, upa = 0, nemo = 0, tcs = 0.35),

  # 3. Dupilumab 300mg Q2W (SOLO trials)
  list(name = "3. Dupilumab Q2W",
       dup = 1, upa = 0, nemo = 0, tcs = 0.0),

  # 4. Upadacitinib 30mg QD (Rising Up trial)
  list(name = "4. Upadacitinib 30mg QD",
       dup = 0, upa = 1, nemo = 0, tcs = 0.0),

  # 5. Nemolizumab 60mg Q4W (itch-focused)
  list(name = "5. Nemolizumab Q4W",
       dup = 0, upa = 0, nemo = 1, tcs = 0.0),

  # 6. Dupilumab + TCS combination
  list(name = "6. Dupilumab + TCS",
       dup = 1, upa = 0, nemo = 0, tcs = 0.35)
)

results <- map_dfr(scenarios, function(s) {
  message("Simulating: ", s$name)
  run_scenario(s$name, s$dup, s$upa, s$nemo, s$tcs)
})

## =============================================================================
## Key Clinical Response Endpoints at Week 16 & 52
## =============================================================================

week_summary <- results %>%
  filter(time %in% c(0, 112, 224, 364)) %>%
  group_by(scenario, time) %>%
  summarise(
    EASI_mean    = round(mean(EASI), 2),
    NRS_mean     = round(mean(NRS_itch), 2),
    IGA_01_pct   = round(mean(IGA <= 1) * 100, 1),
    EASI75_pct   = round(mean(EASI75_achieved) * 100, 1),
    EASI90_pct   = round(mean(EASI90_achieved) * 100, 1),
    RO_pct       = round(mean(RO_pct), 1),
    TARC_pct_chg = round(mean(TARC_pct_change), 1),
    Eos_pct_chg  = round(mean(Eos_pct_change), 1),
    FLG          = round(mean(FLG), 3),
    TEWL         = round(mean(TEWL), 1),
    .groups = "drop"
  ) %>%
  mutate(week = time / 7)

print(week_summary)

## =============================================================================
## Plots
## =============================================================================

# Color palette
scenario_colors <- c(
  "1. No Treatment"        = "#e74c3c",
  "2. TCS Only"            = "#f39c12",
  "3. Dupilumab Q2W"       = "#2980b9",
  "4. Upadacitinib 30mg QD"= "#8e44ad",
  "5. Nemolizumab Q4W"     = "#27ae60",
  "6. Dupilumab + TCS"     = "#1a5276"
)

# --- Plot 1: EASI over time ---
p_easi <- ggplot(results, aes(x = time / 7, y = EASI, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  geom_hline(yintercept = c(7, 14, 21), linetype = "dashed",
             color = c("#27ae60", "#f39c12", "#e74c3c"), alpha = 0.6) +
  annotate("text", x = 1, y = c(5.5, 12.5, 19.5, 27.5),
           label = c("Mild (≤7)", "Mod (≤14)", "Severe (≤21)", "V-Severe (>21)"),
           hjust = 0, size = 3, color = "gray40") +
  labs(title = "Atopic Dermatitis QSP — EASI Score Over 52 Weeks",
       subtitle = "6 treatment scenarios including no treatment, TCS, biologics, JAK inhibitor",
       x = "Week", y = "EASI Score", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

# --- Plot 2: NRS Itch over time ---
p_nrs <- ggplot(results, aes(x = time / 7, y = NRS_itch, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Pruritus NRS (0-10) Over 52 Weeks",
       subtitle = "Nemolizumab shows fastest itch reduction; JAK inhibitor broadband",
       x = "Week", y = "NRS Itch Score (0-10)", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  scale_y_continuous(limits = c(0, 10.5)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# --- Plot 3: TARC/CCL17 biomarker ---
p_tarc <- ggplot(results, aes(x = time / 7, y = TARC, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  geom_hline(yintercept = 450, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = 1, y = 300, label = "Normal range (<450 pg/mL)",
           hjust = 0, size = 3, color = "darkgreen") +
  labs(title = "TARC/CCL17 Biomarker Over 52 Weeks",
       subtitle = "Key Th2 biomarker; IL-4Rα blockade (dupilumab) strongly suppresses TARC",
       x = "Week", y = "TARC (pg/mL)", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# --- Plot 4: Blood eosinophil count ---
p_eos <- ggplot(results, aes(x = time / 7, y = Eos_blood, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = scenario_colors) +
  geom_hline(yintercept = c(500, 1500), linetype = "dashed",
             color = c("#f39c12", "#e74c3c"), alpha = 0.7) +
  labs(title = "Blood Eosinophil Count Over 52 Weeks",
       subtitle = "Dupilumab transiently increases eosinophils before long-term normalization",
       x = "Week", y = "Eosinophils (cells/μL)", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# --- Plot 5: Receptor occupancy (Dupilumab only) ---
p_ro <- results %>%
  filter(grepl("Dupilumab", scenario)) %>%
  ggplot(aes(x = time / 7, y = RO_pct, color = scenario)) +
  geom_line(linewidth = 1.3) +
  geom_hline(yintercept = c(70, 90), linetype = "dashed",
             color = c("#f39c12", "#27ae60")) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IL-4Rα Receptor Occupancy (RO%) — Dupilumab Scenarios",
       subtitle = "Target RO >70% for clinical efficacy; Q2W dosing maintains >85%",
       x = "Week", y = "IL-4Rα Receptor Occupancy (%)", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 12)

# --- Plot 6: Filaggrin and TEWL ---
p_barrier <- results %>%
  select(time, scenario, FLG, TEWL) %>%
  pivot_longer(c(FLG, TEWL), names_to = "marker", values_to = "value") %>%
  ggplot(aes(x = time / 7, y = value, color = scenario)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~marker, scales = "free_y",
             labeller = labeller(marker = c(FLG = "Filaggrin (rel. expression)",
                                            TEWL = "TEWL (g/m²/h)"))) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Epidermal Barrier Recovery",
       subtitle = "Dupilumab restores FLG via STAT6 inhibition and reduces TEWL",
       x = "Week", y = "Value", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 52, 8)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Print all plots
print(p_easi)
print(p_nrs)
print(p_tarc)
print(p_eos)
print(p_ro)
print(p_barrier)

## =============================================================================
## Sensitivity Analysis: Dupilumab Dose (150 vs 300 vs 600 mg Q2W)
## =============================================================================

sens_results <- map_dfr(c(150, 300, 600), function(dose_load) {
  dose_maint <- dose_load / 2
  ev_sc <- ev(amt = dose_load, cmt = "DUP_SC", time = 0)
  maint <- ev(amt = dose_maint, cmt = "DUP_SC",
              time = seq(14, by = 14, length.out = 27))
  out <- ad_mod %>%
    param(c(DUP_on = 1, UPA_on = 0, NEMO_on = 0, TCS_effect = 0)) %>%
    ev(c(ev_sc, maint)) %>%
    mrgsim(delta = 1, end = 364)
  as_tibble(out) %>%
    mutate(dose_label = paste0("Load: ", dose_load, "mg / Maint: ", dose_maint, "mg Q2W"))
})

p_sens <- ggplot(sens_results, aes(x = time / 7, y = EASI, color = dose_label)) +
  geom_line(linewidth = 1.3) +
  scale_color_brewer(palette = "Blues", direction = 1) +
  labs(title = "Sensitivity Analysis: Dupilumab Dose-Response (EASI)",
       subtitle = "Standard 600/300mg regimen outperforms lower doses",
       x = "Week", y = "EASI Score", color = "Dupilumab Dosing") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_sens)

## =============================================================================
## PK/PD Diagnostics
## =============================================================================

# Dupilumab PK profile (first 12 weeks)
pk_diag <- results %>%
  filter(scenario == "3. Dupilumab Q2W", time <= 84) %>%
  ggplot(aes(x = time / 7, y = Cp_dup_ugmL)) +
  geom_line(color = "#2980b9", linewidth = 1.3) +
  geom_hline(yintercept = c(0.027, 2.0), linetype = c("dashed", "dotted"),
             color = c("red", "darkgreen")) +
  annotate("text", x = 0.5, y = 0.05, label = "Kd IL-4Rα ~27 ng/mL",
           hjust = 0, size = 3, color = "red") +
  annotate("text", x = 0.5, y = 2.3, label = "Trough target >2 μg/mL",
           hjust = 0, size = 3, color = "darkgreen") +
  labs(title = "Dupilumab PK Profile — First 12 Weeks",
       subtitle = "2-compartment SC model; Q2W 600mg load → 300mg maintenance",
       x = "Week", y = "Dupilumab Cp (μg/mL)") +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  theme_bw(base_size = 12)

print(pk_diag)

message("Simulation complete. Key results at Week 16:")
print(week_summary %>% filter(week == 16) %>%
      select(scenario, EASI_mean, NRS_mean, EASI75_pct, IGA_01_pct))
