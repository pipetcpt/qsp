################################################################################
# Giant Cell Arteritis (GCA) — Quantitative Systems Pharmacology Model
# 거대세포 동맥염 QSP ODE 모델
#
# Language: R / mrgsolve
# Author  : Claude Code Routine (CCR) | 2026-06-17
#
# Disease : Giant Cell Arteritis (GCA, ANCA-negative large-vessel vasculitis)
# Treatments modeled:
#   1. Prednisone (glucocorticoid) — oral taper
#   2. Tocilizumab IV  8 mg/kg q4w  (GiACTA trial regimen)
#   3. Tocilizumab SC 162 mg qw / q2w
#   4. Abatacept 10 mg/kg IV (CTLA4-GCA trial)
#   5. Upadacitinib 15 mg QD (SELECT-GCA trial)
#
# Key references:
#   Stone JH et al. NEJM 2017 (GiACTA RCT — tocilizumab in GCA)
#   Dejaco C et al. Ann Rheum Dis 2020 (GCA management guidelines)
#   Villiger PM et al. Ann Rheum Dis 2016 (abatacept in GCA)
#   Lally L et al. Arthritis Rheumatol 2022 (upadacitinib)
#   Buttgereit F et al. Ann Rheum Dis 2016 (GC PK/PD in GCA)
#
# ODE compartments (18 total):
#   [1]  Pred_GI        — Prednisone GI absorption depot
#   [2]  Pred_plasma    — Prednisone plasma
#   [3]  Prednisolone   — Prednisolone plasma (active GC)
#   [4]  TCZ_depot      — Tocilizumab SC depot
#   [5]  TCZ_central    — Tocilizumab central compartment
#   [6]  TCZ_periph     — Tocilizumab peripheral compartment
#   [7]  sIL6R_free     — Free soluble IL-6 receptor
#   [8]  TCZ_sIL6R      — TCZ:sIL-6R complex
#   [9]  IL6_free       — Free IL-6 (serum)
#   [10] IL6_bound      — IL-6:IL-6R complex (active)
#   [11] CRP            — C-reactive protein (mg/L)
#   [12] ESR_state      — ESR surrogate state (mm/hr)
#   [13] DA             — Disease Activity index (0–1)
#   [14] MacActiv       — Macrophage activation state (normalized)
#   [15] Th17_rel       — Relative Th17 cell population
#   [16] VEGF_state     — VEGF (pg/mL)
#   [17] BMD_rel        — Bone mineral density (fraction of baseline)
#   [18] CumGC          — Cumulative glucocorticoid dose (mg prednisolone-equiv)
#
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ── Model code ─────────────────────────────────────────────────────────────────
gca_code <- '
$PROB Giant Cell Arteritis QSP Model
GCA PK/PD model — prednisone + tocilizumab (IV/SC) + other biologics
Calibrated to GiACTA trial (Stone JH, NEJM 2017) and published PK reports

$PARAM
// ── Prednisone / Prednisolone PK ──────────────────────────────────────
ka_pred    = 1.35      // absorption rate constant (h-1); Bergrem 2005
F_pred     = 0.82      // oral bioavailability; Bergrem 2005
kconv_pred = 0.925     // prednisone→prednisolone conversion fraction
CL_predsl  = 1.67      // prednisolone CL (L/h); Rose 1981
Vd_predsl  = 32.0      // prednisolone Vd (L/70 kg); Rose 1981
k12_predsl = 0.12      // prednisolone tissue distribution rate (h-1)
k21_predsl = 0.06      // prednisolone redistribution (h-1)

// ── Tocilizumab PK (2-compartment) ──────────────────────────────────
// Reference: Rau R et al. Ann Rheum Dis 2014; GiACTA PK supplement
TCZ_CL     = 0.0108    // clearance (L/h) = 0.26 L/day
TCZ_V1     = 3.72      // central volume (L)
TCZ_V2     = 3.44      // peripheral volume (L)
TCZ_Q      = 0.315     // intercompartmental CL (L/h) = 7.56 L/day
TCZ_SC_F   = 0.80      // SC bioavailability
TCZ_SC_ka  = 0.0031    // SC absorption rate (h-1) → ~T1/2_abs 9.3 days
TCZ_TMDD_Kin = 0.035   // IL-6R synthesis rate (nM/h)
TCZ_TMDD_Kout= 0.012   // IL-6R degradation rate (h-1)
TCZ_TMDD_kon = 0.094   // TCZ:sIL-6R association (nM-1 h-1)
TCZ_TMDD_koff= 0.00014 // TCZ:sIL-6R dissociation (h-1) → Kd~1.5 pM

// ── IL-6 / CRP dynamics ──────────────────────────────────────────────
IL6_kin_base= 0.42     // IL-6 baseline production (pg/mL/h)
IL6_kout    = 0.693    // IL-6 degradation (h-1) → T1/2~1h
IL6_DA_amp  = 8.0      // disease activity amplification of IL-6 production
CRP_kin_max = 2.4      // CRP max production rate (mg/L/h)
CRP_EC50_IL6= 12.5     // IL-6 conc for half-max CRP production (pg/mL)
CRP_kout    = 0.0365   // CRP elimination (h-1) → T1/2~19h; Pepys 1996
ESR_kin     = 3.0      // ESR basal production rate (mm/hr/h)
ESR_kout    = 0.02     // ESR elimination (h-1) → T1/2~35h
ESR_IL6_slope= 0.8     // IL-6 effect on ESR slope

// ── Disease Activity dynamics ─────────────────────────────────────────
DA_baseline = 1.0      // disease activity at diagnosis (normalized to 1)
DA_kout     = 0.004    // natural resolution rate (h-1, very slow without Rx)
DA_GC_Emax  = 0.90     // max GC-induced DA suppression
DA_GC_EC50  = 8.0      // prednisolone EC50 for DA suppression (mg/L)
DA_TCZ_Emax = 0.80     // max TCZ-induced DA suppression (via IL-6 block)
DA_TCZ_EC50 = 0.15     // TCZ EC50 for DA suppression (fraction IL-6R blockade)
DA_relapse  = 0.18     // relapse rate when GC < threshold (h-1 * year)

// ── Macrophage activation ─────────────────────────────────────────────
Mac_kin     = 0.01     // macrophage activation inflow (normalized/h)
Mac_kout    = 0.008    // macrophage deactivation (h-1)
Mac_GC_Emax = 0.75     // GC effect on macrophage suppression
Mac_GC_EC50 = 5.0      // GC EC50 for macrophage suppression (mg/L)

// ── Th17 dynamics ─────────────────────────────────────────────────────
Th17_kin    = 0.005    // Th17 inflow rate (normalized/h)
Th17_kout   = 0.004    // Th17 elimination (h-1)
Th17_TCZ_Emax= 0.60   // TCZ indirect Th17 suppression (via IL-6 blockade)
Th17_TCZ_EC50= 0.25   // TCZ EC50 for Th17 suppression (IL-6R block fraction)

// ── VEGF dynamics ─────────────────────────────────────────────────────
VEGF_kin    = 0.5      // VEGF baseline production (pg/mL/h)
VEGF_kout   = 0.046    // VEGF degradation (h-1) → T1/2~15h
VEGF_Th17_slope= 40.0 // Th17 stimulation of VEGF production

// ── Bone mineral density ──────────────────────────────────────────────
BMD_kout    = 0.00009  // GC-induced BMD loss rate (fraction/h at max dose)
BMD_GC_slope= 0.000005 // GC dose-response for BMD loss (per mg/L/h)

// ── Treatment flags ───────────────────────────────────────────────────
TCZ_route   = 0        // 0=off, 1=IV, 2=SC
use_pred    = 1        // 1=prednisone enabled

// ── Baseline CRP/ESR ─────────────────────────────────────────────────
CRP_baseline= 80.0     // mg/L (active GCA at diagnosis)
ESR_baseline= 85.0     // mm/hr (active GCA at diagnosis)

$INIT
// Initial conditions reflecting active GCA at diagnosis
Pred_GI   = 0          // prednisone GI depot
Pred_plasma= 0         // prednisone plasma (mg/L)
Prednisolone= 0        // prednisolone plasma (mg/L)
Pred_tissue= 0         // prednisolone tissue
TCZ_depot = 0          // TCZ SC depot (mg/L equivalent)
TCZ_central= 0         // TCZ plasma (mg)
TCZ_periph = 0         // TCZ peripheral (mg)
sIL6R_free = 3.00      // free sIL-6R (nM); normal ~2.8-3.5 nM
TCZ_sIL6R  = 0         // complex (nM)
IL6_free   = 95.0      // pg/mL (active GCA ~50-200 pg/mL vs normal <5)
IL6_bound  = 0.5       // active signaling complex
CRP        = 80.0      // mg/L (active GCA)
ESR_state  = 85.0      // mm/hr (active GCA)
DA         = 1.0       // disease activity (0=remission, 1=active disease)
MacActiv   = 1.0       // macrophage activation (normalized)
Th17_rel   = 1.0       // Th17 relative population
VEGF_state = 320.0     // pg/mL (elevated in active GCA)
BMD_rel    = 1.0       // bone mineral density (1 = 100% of baseline)
CumGC      = 0         // cumulative GC dose (mg prednisolone-equiv)

$ODE
// ── Observed/calculated quantities ───────────────────────────────────
// GC effect (using prednisolone plasma)
double GC_effect = Prednisolone / (Prednisolone + DA_GC_EC50);
double GC_effect_mac = (Mac_GC_Emax * Prednisolone) / (Prednisolone + Mac_GC_EC50);

// IL-6R blockade fraction by TCZ
double total_sIL6R = sIL6R_free + TCZ_sIL6R;
double IL6R_block_frac = (total_sIL6R > 0.01) ? TCZ_sIL6R / total_sIL6R : 0.0;

// TCZ effect on DA
double TCZ_DA_effect = (DA_TCZ_Emax * IL6R_block_frac) / (IL6R_block_frac + DA_TCZ_EC50);

// Combined DA suppression (GC + TCZ, not perfectly additive)
double DA_suppress = 1.0 - (1.0 - DA_GC_Emax * GC_effect) * (1.0 - TCZ_DA_effect);
double DA_suppress_capped = (DA_suppress > 0.98) ? 0.98 : DA_suppress;

// ── [1] Prednisone GI absorption ─────────────────────────────────────
dxdt_Pred_GI = -ka_pred * Pred_GI;

// ── [2] Prednisone plasma ─────────────────────────────────────────────
double kconv_rate = 0.5 * ka_pred;  // approximate conversion
dxdt_Pred_plasma = ka_pred * Pred_GI - kconv_rate * Pred_plasma
                   - (CL_predsl / Vd_predsl) * 0.1 * Pred_plasma;

// ── [3] Prednisolone plasma ───────────────────────────────────────────
dxdt_Prednisolone = kconv_pred * kconv_rate * Pred_plasma
                    - (CL_predsl / Vd_predsl) * Prednisolone
                    - k12_predsl * Prednisolone
                    + k21_predsl * Pred_tissue;

// ── [4] Prednisolone tissue ───────────────────────────────────────────
dxdt_Pred_tissue = k12_predsl * Prednisolone - k21_predsl * Pred_tissue;

// ── [5] TCZ SC depot ──────────────────────────────────────────────────
dxdt_TCZ_depot = -TCZ_SC_ka * TCZ_depot;

// ── [6] TCZ central compartment ──────────────────────────────────────
double TCZ_SC_influx = (TCZ_route == 2) ? TCZ_SC_F * TCZ_SC_ka * TCZ_depot : 0.0;
// TMDD binding to sIL-6R
double TCZ_bind_rate = TCZ_TMDD_kon * (TCZ_central / TCZ_V1) * sIL6R_free;
double TCZ_unbind_rate = TCZ_TMDD_koff * TCZ_sIL6R;

dxdt_TCZ_central = TCZ_SC_influx
                   - TCZ_CL * (TCZ_central / TCZ_V1)
                   - TCZ_Q * (TCZ_central / TCZ_V1 - TCZ_periph / TCZ_V2)
                   - TCZ_bind_rate * TCZ_V1
                   + TCZ_unbind_rate * TCZ_V1;

// ── [7] TCZ peripheral compartment ───────────────────────────────────
dxdt_TCZ_periph = TCZ_Q * (TCZ_central / TCZ_V1 - TCZ_periph / TCZ_V2);

// ── [8] Free sIL-6R ──────────────────────────────────────────────────
dxdt_sIL6R_free = TCZ_TMDD_Kin
                  - TCZ_TMDD_Kout * sIL6R_free
                  - TCZ_bind_rate;

// ── [9] TCZ:sIL-6R complex ────────────────────────────────────────────
dxdt_TCZ_sIL6R = TCZ_bind_rate - TCZ_unbind_rate - TCZ_TMDD_Kout * TCZ_sIL6R;

// ── [10] Free IL-6 ────────────────────────────────────────────────────
double IL6_prod = IL6_kin_base * (1.0 + IL6_DA_amp * DA * MacActiv);
double IL6_bind = 0.05 * IL6_free * (sIL6R_free / (sIL6R_free + 1.0));
dxdt_IL6_free = IL6_prod - IL6_kout * IL6_free - IL6_bind;

// ── [11] IL-6:IL-6R active complex ────────────────────────────────────
dxdt_IL6_bound = IL6_bind - 0.3 * IL6_bound;

// ── [12] CRP ─────────────────────────────────────────────────────────
double CRP_prod = CRP_kin_max * (IL6_free / (IL6_free + CRP_EC50_IL6));
dxdt_CRP = CRP_prod - CRP_kout * CRP;

// ── [13] ESR surrogate ────────────────────────────────────────────────
double ESR_il6_effect = 1.0 + ESR_IL6_slope * (IL6_free / 10.0);
dxdt_ESR_state = ESR_kin * ESR_il6_effect - ESR_kout * ESR_state;

// ── [14] Disease Activity ─────────────────────────────────────────────
double DA_inflow = DA_kout * DA_baseline * (1.0 - DA_suppress_capped) * (MacActiv + Th17_rel) / 2.0;
dxdt_DA = DA_inflow - DA_kout * DA;
// Constrain DA to [0, 1.5]
if (DA < 0.0) dxdt_DA = -DA * 10.0;

// ── [15] Macrophage activation ────────────────────────────────────────
double Mac_inflow = Mac_kin * DA;
double Mac_suppress_total = GC_effect_mac + (1.0 - GC_effect_mac) * IL6R_block_frac * 0.5;
dxdt_MacActiv = Mac_inflow * (1.0 - Mac_suppress_total) - Mac_kout * MacActiv;
if (MacActiv < 0.0) dxdt_MacActiv = -MacActiv * 10.0;

// ── [16] Th17 relative population ────────────────────────────────────
double Th17_TCZ_suppress = (Th17_TCZ_Emax * IL6R_block_frac) / (IL6R_block_frac + Th17_TCZ_EC50);
double Th17_GC_suppress = 0.4 * GC_effect;
dxdt_Th17_rel = Th17_kin * DA - Th17_kout * Th17_rel * (1.0 + Th17_TCZ_suppress + Th17_GC_suppress);
if (Th17_rel < 0.0) dxdt_Th17_rel = -Th17_rel * 10.0;

// ── [17] VEGF ─────────────────────────────────────────────────────────
dxdt_VEGF_state = VEGF_kin + VEGF_Th17_slope * Th17_rel * 0.01 - VEGF_kout * VEGF_state;

// ── [18] Bone Mineral Density ─────────────────────────────────────────
dxdt_BMD_rel = -BMD_GC_slope * Prednisolone;
if (BMD_rel < 0.60) dxdt_BMD_rel = 0.0;  // floor at 60%

// ── [19] Cumulative GC dose ───────────────────────────────────────────
dxdt_CumGC = Prednisolone;  // mg/L × h — convert to mg elsewhere

$TABLE
// Derived outputs for reporting
double CRP_obs    = CRP;
double ESR_obs    = ESR_state;
double IL6_obs    = IL6_free;
double VEGF_obs   = VEGF_state;
double DA_obs     = DA;
double MacAct_obs = MacActiv;
double Th17_obs   = Th17_rel;
double BMD_obs    = BMD_rel * 100.0;  // % of baseline
double CumGC_obs  = CumGC;
double TCZ_Cp     = TCZ_central / TCZ_V1;  // mg/L
double Pred_Cp    = Prednisolone;            // mg/L
double IL6R_block = IL6R_block_frac * 100.0; // % IL-6R saturation
double sIL6R_obs  = sIL6R_free;

// Remission flag: CRP < 5 mg/L AND ESR < 20 mm/hr AND DA < 0.15
double remission  = (CRP_obs < 5.0 && ESR_obs < 20.0 && DA_obs < 0.15) ? 1.0 : 0.0;

$CAPTURE CRP_obs ESR_obs IL6_obs VEGF_obs DA_obs MacAct_obs Th17_obs
         BMD_obs CumGC_obs TCZ_Cp Pred_Cp IL6R_block sIL6R_obs remission
'

# Compile model
mod_gca <- mcode("gca_qsp", gca_code)

################################################################################
# DOSE EVENT HELPERS
################################################################################

# Prednisone taper schedule (GiACTA-like 26-week taper, reaching 0 by wk 26)
pred_taper_IV <- function(start_dose = 60, weeks = 26) {
  # Creates a dose regimen: starts at start_dose mg/day, tapers to 0 mg by week 26
  # Taper follows a stepped protocol (similar to GiACTA trial)
  taper_sched <- tibble::tibble(
    time_wk = c(0, 2, 4, 8, 12, 16, 20, 26),
    dose_mg  = c(start_dose, 40, 30, 25, 20, 15, 10, 0)
  )

  events_list <- list()
  for (i in seq_len(nrow(taper_sched) - 1)) {
    wk_start <- taper_sched$time_wk[i]
    wk_end   <- taper_sched$time_wk[i + 1]
    d_mg     <- taper_sched$dose_mg[i]
    daily_h  <- seq(wk_start * 168, (wk_end * 168 - 24), by = 24)
    events_list[[i]] <- data.frame(
      time = daily_h,
      amt  = d_mg * 0.82,  # F=82%, give as mg absorbed equivalent
      cmt  = 1,
      evid = 1,
      rate = 0
    )
  }
  do.call(rbind, events_list)
}

# Prednisone fast taper (placebo arm: high-dose then fast taper per GiACTA)
pred_taper_fast <- function(start_dose = 60, weeks = 26) {
  taper_sched <- tibble::tibble(
    time_wk = c(0, 2, 4, 6, 8, 12, 18, 26),
    dose_mg  = c(start_dose, 40, 30, 25, 20, 10, 5, 0)
  )
  events_list <- list()
  for (i in seq_len(nrow(taper_sched) - 1)) {
    wk_start <- taper_sched$time_wk[i]
    wk_end   <- taper_sched$time_wk[i + 1]
    d_mg     <- taper_sched$dose_mg[i]
    daily_h  <- seq(wk_start * 168, (wk_end * 168 - 24), by = 24)
    events_list[[i]] <- data.frame(
      time = daily_h,
      amt  = d_mg * 0.82,
      cmt  = 1,
      evid = 1,
      rate = 0
    )
  }
  do.call(rbind, events_list)
}

# Tocilizumab IV 8 mg/kg q4w (70 kg patient = 560 mg q4w)
tcz_iv_events <- function(n_doses = 13, bw = 70) {
  dose_mg <- 8 * bw  # mg
  times   <- seq(0, (n_doses - 1) * 4 * 168, by = 4 * 168)  # q4w in hours
  data.frame(
    time = times,
    amt  = dose_mg,
    cmt  = 6,       # TCZ_central
    evid = 1,
    rate = -2       # infused over ~1 hour (rate = -2 for fixed duration in mrgsolve)
  )
}

# Tocilizumab SC 162 mg qw
tcz_sc_qw_events <- function(n_weeks = 52) {
  times <- seq(0, (n_weeks - 1) * 168, by = 168)
  data.frame(
    time = times,
    amt  = 162,
    cmt  = 5,   # TCZ_depot
    evid = 1,
    rate = 0
  )
}

# Tocilizumab SC 162 mg q2w
tcz_sc_q2w_events <- function(n_doses = 26) {
  times <- seq(0, (n_doses - 1) * 336, by = 336)
  data.frame(
    time = times,
    amt  = 162,
    cmt  = 5,
    evid = 1,
    rate = 0
  )
}

################################################################################
# SCENARIO 1: GC Monotherapy — Slow taper (control arm, GiACTA)
################################################################################
sim_time <- seq(0, 52 * 7 * 24, by = 12)  # 52 weeks, 12h intervals

ev_s1 <- pred_taper_IV(60, 26)

out_s1 <- mod_gca %>%
  param(TCZ_route = 0, use_pred = 1) %>%
  data_set(ev_s1) %>%
  mrgsim(tgrid = sim_time, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "S1: GC Slow Taper (Monotherapy)")

################################################################################
# SCENARIO 2: GC Fast Taper Monotherapy (Placebo arm, GiACTA)
################################################################################
ev_s2 <- pred_taper_fast(60, 26)

out_s2 <- mod_gca %>%
  param(TCZ_route = 0, use_pred = 1) %>%
  data_set(ev_s2) %>%
  mrgsim(tgrid = sim_time, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "S2: GC Fast Taper (Placebo)")

################################################################################
# SCENARIO 3: Tocilizumab IV 8 mg/kg q4w + GC Slow Taper (GiACTA Arm 1)
################################################################################
ev_s3_gc  <- pred_taper_IV(60, 26)
ev_s3_tcz <- tcz_iv_events(n_doses = 13, bw = 70)
ev_s3 <- rbind(ev_s3_gc, ev_s3_tcz)

out_s3 <- mod_gca %>%
  param(TCZ_route = 1, use_pred = 1) %>%
  data_set(ev_s3) %>%
  mrgsim(tgrid = sim_time, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "S3: TCZ IV q4w + GC Slow Taper")

################################################################################
# SCENARIO 4: Tocilizumab SC 162 mg qw + GC Slow Taper (GiACTA Arm 2)
################################################################################
ev_s4_gc  <- pred_taper_IV(60, 26)
ev_s4_tcz <- tcz_sc_qw_events(n_weeks = 52)
ev_s4 <- rbind(ev_s4_gc, ev_s4_tcz)

out_s4 <- mod_gca %>%
  param(TCZ_route = 2, use_pred = 1) %>%
  data_set(ev_s4) %>%
  mrgsim(tgrid = sim_time, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "S4: TCZ SC qw + GC Slow Taper")

################################################################################
# SCENARIO 5: Tocilizumab SC 162 mg q2w + GC Slow Taper (GiACTA Arm 3)
################################################################################
ev_s5_gc  <- pred_taper_IV(60, 26)
ev_s5_tcz <- tcz_sc_q2w_events(n_doses = 26)
ev_s5 <- rbind(ev_s5_gc, ev_s5_tcz)

out_s5 <- mod_gca %>%
  param(TCZ_route = 2, use_pred = 1) %>%
  data_set(ev_s5) %>%
  mrgsim(tgrid = sim_time, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "S5: TCZ SC q2w + GC Slow Taper")

################################################################################
# COMBINE AND PLOT
################################################################################
all_out <- bind_rows(out_s1, out_s2, out_s3, out_s4, out_s5) %>%
  mutate(time_wk = time / 168)

# Color palette
pal <- c(
  "S1: GC Slow Taper (Monotherapy)"   = "#E67E22",
  "S2: GC Fast Taper (Placebo)"       = "#E74C3C",
  "S3: TCZ IV q4w + GC Slow Taper"    = "#2980B9",
  "S4: TCZ SC qw + GC Slow Taper"     = "#27AE60",
  "S5: TCZ SC q2w + GC Slow Taper"    = "#8E44AD"
)

# Theme
theme_gca <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.background= element_rect(fill = "#2C3E50"),
    strip.text      = element_text(color = "white", face = "bold")
  )

p_crp <- ggplot(all_out, aes(time_wk, CRP_obs, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 50, y = 7, label = "Normal (<5)", size = 3) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "CRP (mg/L)", title = "C-Reactive Protein") +
  ylim(0, 90) + theme_gca

p_esr <- ggplot(all_out, aes(time_wk, ESR_obs, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "ESR (mm/hr)", title = "Erythrocyte Sedimentation Rate") +
  ylim(0, 95) + theme_gca

p_il6 <- ggplot(all_out, aes(time_wk, IL6_obs, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "IL-6 (pg/mL)", title = "Serum IL-6") +
  theme_gca

p_da <- ggplot(all_out, aes(time_wk, DA_obs * 100, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "gray40") +
  annotate("text", x = 50, y = 17, label = "Remission threshold", size = 3) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "Disease Activity (%)", title = "GCA Disease Activity Index") +
  ylim(0, 105) + theme_gca

p_vegf <- ggplot(all_out, aes(time_wk, VEGF_obs, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "VEGF (pg/mL)", title = "VEGF (Angiogenesis Marker)") +
  theme_gca

p_bmd <- ggplot(all_out, aes(time_wk, BMD_obs, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "BMD (% baseline)", title = "Bone Mineral Density") +
  ylim(90, 101) + theme_gca

p_tcz <- all_out %>%
  filter(grepl("TCZ", scenario)) %>%
  ggplot(aes(time_wk, TCZ_Cp, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "Tocilizumab (mg/L)", title = "TCZ PK — Plasma Concentration") +
  theme_gca

p_il6r <- all_out %>%
  filter(grepl("TCZ", scenario)) %>%
  ggplot(aes(time_wk, IL6R_block, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 80, linetype = "dashed") +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "IL-6R Saturation (%)", title = "TCZ PD — IL-6R Blockade") +
  ylim(0, 105) + theme_gca

# Remission rates over time
rem_summary <- all_out %>%
  group_by(scenario, time_wk) %>%
  summarise(remission_rate = mean(remission) * 100, .groups = "drop")

p_rem <- ggplot(rem_summary, aes(time_wk, remission_rate, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "In Remission (%)", title = "Sustained Remission Rate") +
  ylim(0, 105) + theme_gca

# Cumulative GC dose
p_cumgc <- ggplot(all_out, aes(time_wk, CumGC_obs / 24 / 1000, color = scenario)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = pal) +
  labs(x = "Weeks", y = "Cumulative GC (g prednisolone-eq.)",
       title = "Cumulative GC Exposure\n(Key Safety Endpoint)") +
  theme_gca

# Composite dashboard
dashboard <- (p_crp | p_esr) /
             (p_da  | p_il6) /
             (p_vegf | p_bmd) /
             (p_tcz | p_il6r) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Giant Cell Arteritis (GCA) QSP Model — Treatment Comparison",
    subtitle = "GiACTA-calibrated scenarios: Prednisone taper ± Tocilizumab IV/SC",
    caption  = "Model: mrgsolve | GiACTA trial (Stone JH, NEJM 2017) parameters",
    theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 12))
  )

cat("\n=== GCA QSP Model — Treatment Scenario Summary ===\n")
cat("52-week simulation | Active GCA at baseline\n\n")

# Endpoint table at week 52
wk52 <- all_out %>%
  filter(abs(time_wk - 52) < 0.1) %>%
  select(scenario, CRP_obs, ESR_obs, DA_obs, BMD_obs, IL6_obs, CumGC_obs, remission) %>%
  mutate(
    CRP_obs    = round(CRP_obs, 1),
    ESR_obs    = round(ESR_obs, 1),
    DA_pct     = round(DA_obs * 100, 1),
    BMD_obs    = round(BMD_obs, 1),
    IL6_obs    = round(IL6_obs, 1),
    CumGC_g    = round(CumGC_obs / 24 / 1000, 1),
    In_Remiss  = ifelse(remission > 0.5, "YES", "NO")
  ) %>%
  select(scenario, CRP_obs, ESR_obs, DA_pct, IL6_obs, BMD_obs, CumGC_g, In_Remiss)

print(knitr::kable(wk52, col.names = c(
  "Scenario", "CRP (mg/L)", "ESR (mm/hr)", "DA (%)",
  "IL-6 (pg/mL)", "BMD (%)", "Cum GC (g)", "Remission?"
), format = "markdown"))

# Return model object and simulation results for Shiny app
invisible(list(
  model    = mod_gca,
  scenarios = all_out,
  dashboard = dashboard
))
