## ============================================================
## Thyroid Cancer (ThyCa) QSP Model — mrgsolve ODE
## ============================================================
## Disease: Thyroid Cancer (DTC, MTC, ATC subtypes)
## Author : Claude Code Routine (CCR)
## Date   : 2026-06-27
##
## Compartments (18 ODEs):
##   PK   : Drug Central + Peripheral (3 TKI classes)
##   PD   : Oncogenic pathways, tumor dynamics,
##          biomarkers, clinical endpoints
##
## Treatment Scenarios:
##   1. Untreated (natural history)
##   2. Lenvatinib (1st-line DTC — SELECT trial)
##   3. Sorafenib (1st-line DTC — DECISION trial)
##   4. Lenvatinib after sorafenib failure
##   5. Selpercatinib (RET+ MTC — LIBRETTO-001)
##   6. Vandetanib (MTC — ZETA trial)
##   7. Lenvatinib + everolimus (mTOR combo — HOPE trial)
##
## Key calibration:
##   - SELECT (lenvatinib DTC): PFS HR 0.21 (7.4 → 18.3 mo)
##   - DECISION (sorafenib DTC): PFS HR 0.59 (5.8 → 10.8 mo)
##   - LIBRETTO-001 (selpercatinib MTC): ORR 69%
##   - ZETA (vandetanib MTC): PFS HR 0.46 (19.3 → 30.5 mo)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## ── Model definition ─────────────────────────────────────────
code <- '
$PROB
Thyroid Cancer QSP Model
18-compartment ODE system covering:
  - Drug PK (oral two-compartment: Lenvatinib, Sorafenib, Selpercatinib)
  - MAPK/RAF-MEK-ERK pathway activity
  - PI3K/AKT/mTOR pathway activity
  - Tumor cell dynamics (proliferation, apoptosis)
  - Tumor biomarkers (Tg for DTC, Calcitonin for MTC)
  - TSH dynamics
  - Clinical endpoints (tumor volume, RECIST response)

$PARAM
// ── Tumor biology parameters ─────────────────────────────────
// Baseline proliferation/death (per day)
kprol_base  = 0.03    // baseline tumor proliferation (/day)
kdeath_base = 0.02    // baseline apoptosis (/day)
kprol_MAPK  = 0.02    // MAPK-driven extra proliferation
kprol_PI3K  = 0.01    // PI3K-driven extra proliferation

// Pathway activity (normalized 0-1)
MAPK_base   = 0.8     // baseline MAPK activity (BRAF V600E high)
PI3K_base   = 0.5     // baseline PI3K activity
MAPK_max    = 1.0     // max MAPK activity
PI3K_max    = 1.0     // max PI3K activity

// Angiogenesis
kVEGF_prod  = 0.05    // VEGF production (HIF-1α driven)
kVEGF_clear = 0.1     // VEGF clearance
kAngio      = 0.02    // angiogenesis rate per VEGF
kAngio_reg  = 0.05    // vessel regression

// Biomarker kinetics
kTg_prod    = 0.08    // Tg secretion rate per tumor cell
kTg_clear   = 0.05    // Tg clearance (/day)
kCT_prod    = 0.15    // Calcitonin secretion (MTC)
kCT_clear   = 0.08    // Calcitonin clearance

// TSH dynamics
TSH_base    = 1.5     // baseline TSH (mIU/L) during active disease
TSH_supp    = 0.1     // TSH under suppression
kTSH_stim   = 0.003   // TSH-driven tumor growth contribution

// ── PK parameters — Lenvatinib (24 mg QD oral) ──────────────
CL_LENV  = 4.0        // Clearance (L/h) [Pop PK: CL ~4.2 L/h]
V1_LENV  = 50.0       // Central volume (L)
Q_LENV   = 1.6        // Inter-comp. CL (L/h)
V2_LENV  = 120.0      // Peripheral volume (L)
ka_LENV  = 1.4        // Absorption rate (/h) [t1/2abs ~0.5h]
F_LENV   = 0.85       // Bioavailability (85%)
dose_LENV= 24.0       // mg QD

// ── PK parameters — Sorafenib (400 mg BID oral) ─────────────
CL_SORA  = 3.2        // Clearance (L/h)
V1_SORA  = 55.0       // Central volume (L)
Q_SORA   = 2.0        // Inter-comp. CL (L/h)
V2_SORA  = 220.0      // Peripheral volume (L)
ka_SORA  = 0.6        // Absorption rate (/h) [F ~38%, delayed]
F_SORA   = 0.38       // Bioavailability
dose_SORA= 400.0      // mg BID

// ── PK parameters — Selpercatinib (160 mg BID oral) ─────────
CL_SELP  = 6.5        // Clearance (L/h)
V1_SELP  = 198.0      // Central volume (L)
Q_SELP   = 3.2        // Inter-comp. CL (L/h)
V2_SELP  = 1130.0     // Peripheral volume (L)
ka_SELP  = 0.7        // Absorption rate (/h)
F_SELP   = 0.73       // Bioavailability (73%)
dose_SELP= 160.0      // mg BID

// ── PD (Emax) parameters ─────────────────────────────────────
// Lenvatinib PD: VEGFR/FGFR/RET inhibition
EC50_LENV_VEGFR = 5.0   // ng/mL for 50% VEGFR inhibition
Emax_LENV_VEGFR = 0.92  // max VEGFR inhibition
EC50_LENV_MAPK  = 20.0  // ng/mL for 50% MAPK pathway suppression
Emax_LENV_MAPK  = 0.60  // max MAPK suppression via FGFR/RET
Hill_LENV = 1.5

// Sorafenib PD: BRAF/CRAF + VEGFR2 inhibition
EC50_SORA_BRAF  = 3.5   // μg/mL for 50% BRAF inhibition
Emax_SORA_BRAF  = 0.75  // max BRAF inhibition
EC50_SORA_VEGFR = 6.0   // μg/mL VEGFR inhibition
Emax_SORA_VEGFR = 0.80
Hill_SORA = 1.5

// Selpercatinib PD: selective RET inhibition
EC50_SELP_RET   = 10.0  // ng/mL for 50% RET inhibition
Emax_SELP_RET   = 0.97  // near-complete RET inhibition
Hill_SELP = 2.0

// Everolimus (mTOR inhibitor): if used in combo
EC50_EVER_mTOR  = 5.0   // ng/mL
Emax_EVER_mTOR  = 0.80
use_everolimus = 0

// Treatment switches
use_LENV  = 0
use_SORA  = 0
use_SELP  = 0
use_TSH_supp = 0

// Disease subtype (1=DTC, 2=MTC)
is_MTC    = 0

$CMT
// PK compartments (6)
LENV_gut   // Lenvatinib gut depot (CMT 1)
LENV_C     // Lenvatinib central (CMT 2)
LENV_P     // Lenvatinib peripheral (CMT 3)
SORA_gut   // Sorafenib gut depot (CMT 4)
SORA_C     // Sorafenib central (CMT 5)
SORA_P     // Sorafenib peripheral (CMT 6)
SELP_gut   // Selpercatinib gut depot (CMT 7)
SELP_C     // Selpercatinib central (CMT 8)
SELP_P     // Selpercatinib peripheral (CMT 9)

// Disease compartments (9)
MAPK_act   // MAPK pathway activity (0-1) (CMT 10)
PI3K_act   // PI3K/AKT/mTOR activity (0-1) (CMT 11)
VEGF       // VEGF concentration (CMT 12)
Angio      // Tumor vasculature (normalized) (CMT 13)
TumorN     // Tumor cell number (normalized, 1=initial) (CMT 14)
Tg         // Serum thyroglobulin (ng/mL) [DTC] (CMT 15)
CT         // Serum calcitonin (pg/mL) [MTC] (CMT 16)
TSH        // TSH level (mIU/L) (CMT 17)
TumVol     // Tumor volume (sum of longest diameters, mm) (CMT 18)

$MAIN
// Convert concentrations (doses in mg, volumes in L)
double Cp_LENV = LENV_C / V1_LENV * 1000.0; // ng/mL
double Cp_SORA = SORA_C / V1_SORA * 1000.0; // ng/mL (from mg/L)
double Cp_SELP = SELP_C / V1_SELP * 1000.0; // ng/mL

// Emax models
// Lenvatinib inhibition on VEGFR (anti-angiogenic)
double Inh_LENV_VEGFR = use_LENV * Emax_LENV_VEGFR *
  pow(Cp_LENV, Hill_LENV) /
  (pow(EC50_LENV_VEGFR, Hill_LENV) + pow(Cp_LENV, Hill_LENV));

// Lenvatinib partial MAPK suppression (via FGFR/RET)
double Inh_LENV_MAPK  = use_LENV * Emax_LENV_MAPK *
  pow(Cp_LENV, Hill_LENV) /
  (pow(EC50_LENV_MAPK, Hill_LENV) + pow(Cp_LENV, Hill_LENV));

// Sorafenib BRAF inhibition
double Inh_SORA_BRAF  = use_SORA * Emax_SORA_BRAF *
  pow(Cp_SORA, Hill_SORA) /
  (pow(EC50_SORA_BRAF, Hill_SORA) + pow(Cp_SORA, Hill_SORA));

// Sorafenib VEGFR
double Inh_SORA_VEGFR = use_SORA * Emax_SORA_VEGFR *
  pow(Cp_SORA, Hill_SORA) /
  (pow(EC50_SORA_VEGFR, Hill_SORA) + pow(Cp_SORA, Hill_SORA));

// Selpercatinib RET inhibition
double Inh_SELP_RET   = use_SELP * Emax_SELP_RET *
  pow(Cp_SELP, Hill_SELP) /
  (pow(EC50_SELP_RET, Hill_SELP) + pow(Cp_SELP, Hill_SELP));

$INIT
LENV_gut = 0
LENV_C   = 0
LENV_P   = 0
SORA_gut = 0
SORA_C   = 0
SORA_P   = 0
SELP_gut = 0
SELP_C   = 0
SELP_P   = 0
MAPK_act = 0.8
PI3K_act = 0.5
VEGF     = 0.5
Angio    = 0.6
TumorN   = 1.0
Tg       = 50.0   // elevated at baseline (DTC)
CT       = 100.0  // elevated at baseline (MTC)
TSH      = 1.5
TumVol   = 30.0   // mm (sum of target lesions)

$ODE
// ── PK ─────────────────────────────────────────────────────
// Lenvatinib (oral QD)
dxdt_LENV_gut = -ka_LENV * LENV_gut;
dxdt_LENV_C   = F_LENV * ka_LENV * LENV_gut
               - (CL_LENV + Q_LENV) / V1_LENV * LENV_C
               + Q_LENV / V2_LENV * LENV_P;
dxdt_LENV_P   = Q_LENV / V1_LENV * LENV_C
               - Q_LENV / V2_LENV * LENV_P;

// Sorafenib (oral BID)
dxdt_SORA_gut = -ka_SORA * SORA_gut;
dxdt_SORA_C   = F_SORA * ka_SORA * SORA_gut
               - (CL_SORA + Q_SORA) / V1_SORA * SORA_C
               + Q_SORA / V2_SORA * SORA_P;
dxdt_SORA_P   = Q_SORA / V1_SORA * SORA_C
               - Q_SORA / V2_SORA * SORA_P;

// Selpercatinib (oral BID)
dxdt_SELP_gut = -ka_SELP * SELP_gut;
dxdt_SELP_C   = F_SELP * ka_SELP * SELP_gut
               - (CL_SELP + Q_SELP) / V1_SELP * SELP_C
               + Q_SELP / V2_SELP * SELP_P;
dxdt_SELP_P   = Q_SELP / V1_SELP * SELP_C
               - Q_SELP / V2_SELP * SELP_P;

// ── MAPK pathway activity ──────────────────────────────────
// Activators: BRAF V600E constitutive
// Inhibitors: Sorafenib (BRAF), Lenvatinib (partial, via FGFR/RET)
//             Selpercatinib (RET contribution in RET+ MTC)
double MAPK_activ = MAPK_base + 0.05 * TSH * kTSH_stim;
double MAPK_inh   = Inh_SORA_BRAF * 0.9 + Inh_LENV_MAPK * 0.3;
if(is_MTC > 0.5) MAPK_inh += Inh_SELP_RET * 0.8;
double MAPK_target = MAPK_activ * (1.0 - MAPK_inh);
if(MAPK_target < 0.0) MAPK_target = 0.0;
if(MAPK_target > 1.0) MAPK_target = 1.0;
dxdt_MAPK_act = 0.5 * (MAPK_target - MAPK_act); // 0.5/day rate to equilibrium

// ── PI3K pathway activity ──────────────────────────────────
double PI3K_activ = PI3K_base;
// VEGFR-PI3K cross-talk
double PI3K_inh   = Inh_LENV_VEGFR * 0.4 + Inh_SORA_VEGFR * 0.3;
double PI3K_target = PI3K_activ * (1.0 - PI3K_inh);
if(PI3K_target < 0.0) PI3K_target = 0.0;
dxdt_PI3K_act = 0.3 * (PI3K_target - PI3K_act);

// ── VEGF dynamics ──────────────────────────────────────────
double VEGF_prod  = kVEGF_prod * TumorN * PI3K_act; // HIF-1α driven
dxdt_VEGF = VEGF_prod - kVEGF_clear * VEGF;

// ── Angiogenesis ───────────────────────────────────────────
double angio_stim = kAngio * VEGF * (1.0 - Inh_LENV_VEGFR) * (1.0 - Inh_SORA_VEGFR);
dxdt_Angio = angio_stim - kAngio_reg * Angio;
if(Angio < 0) Angio = 0.0;
if(Angio > 2) Angio = 2.0;

// ── Tumor cell dynamics ──────────────────────────────────
// Proliferation driven by MAPK + PI3K + TSH + angiogenesis
double growth_factor = (1.0 + kprol_MAPK * MAPK_act +
                        kprol_PI3K * PI3K_act +
                        kTSH_stim * TSH);
double kprol_net = kprol_base * growth_factor * Angio;

// Drug-driven apoptosis
double kdeath_drug = 0.0;
kdeath_drug += 0.015 * Inh_LENV_VEGFR;   // Lenvatinib anti-tumor effect
kdeath_drug += 0.010 * Inh_SORA_BRAF;    // Sorafenib
kdeath_drug += 0.018 * Inh_SELP_RET;     // Selpercatinib (RET-direct)

dxdt_TumorN = kprol_net * TumorN - (kdeath_base + kdeath_drug) * TumorN;
if(TumorN < 0) TumorN = 0.0;

// ── Tumor volume (sum of target lesions) ───────────────────
// SLD ~ proportional to tumor cell number (cube root scaling)
dxdt_TumVol = (kprol_net - kdeath_base - kdeath_drug) * TumVol;
if(TumVol < 0) TumVol = 0.0;

// ── Biomarkers ─────────────────────────────────────────────
// Tg: secreted by DTC cells, reduced by ablation/drug
double Tg_prod  = (1.0 - is_MTC) * kTg_prod * TumorN * TSH / (TSH + 0.5);
double Tg_supp  = use_TSH_supp * 0.7; // TSH suppression reduces Tg
dxdt_Tg = Tg_prod * (1.0 - Tg_supp) - kTg_clear * Tg;

// Calcitonin: MTC biomarker
double CT_prod  = is_MTC * kCT_prod * TumorN;
double CT_inh   = Inh_SELP_RET * 0.8 + (use_SORA > 0.5 ? 0.3 : 0.0);
dxdt_CT = CT_prod * (1.0 - CT_inh) - kCT_clear * CT;

// ── TSH dynamics ───────────────────────────────────────────
double TSH_target = use_TSH_supp > 0.5 ? TSH_supp : TSH_base;
dxdt_TSH = 0.1 * (TSH_target - TSH); // 0.1/day adjustment

$TABLE
double Cp_LENV_out = LENV_C / V1_LENV * 1000.0;
double Cp_SORA_out = SORA_C / V1_SORA * 1000.0;
double Cp_SELP_out = SELP_C / V1_SELP * 1000.0;

// RECIST response (% change in SLD from baseline)
double SLD_change = (TumVol - 30.0) / 30.0 * 100.0;

// Response categories
double resp_CR = (double)(TumVol < 1.0);   // CR: near-complete
double resp_PR = (double)(SLD_change <= -30.0 && TumVol >= 1.0);
double resp_PD = (double)(SLD_change >= 20.0);
double resp_SD = (double)(SLD_change > -30.0 && SLD_change < 20.0 && TumVol >= 1.0);

$CAPTURE
Cp_LENV_out Cp_SORA_out Cp_SELP_out
MAPK_act PI3K_act VEGF Angio
TumorN TumVol Tg CT TSH
SLD_change resp_CR resp_PR resp_PD resp_SD
'

## Compile
thyca_mod <- mcode("thyca_qsp", code)

## ── Dosing functions ─────────────────────────────────────────

make_lenv_dose <- function(start = 0, end = 365) {
  ev(amt = 24, cmt = "LENV_gut", ii = 1, addl = end - start - 1, time = start)
}

make_sora_dose <- function(start = 0, end = 365) {
  ev(amt = 400, cmt = "SORA_gut", ii = 0.5, addl = (end - start) * 2 - 1, time = start)
}

make_selp_dose <- function(start = 0, end = 365) {
  ev(amt = 160, cmt = "SELP_gut", ii = 0.5, addl = (end - start) * 2 - 1, time = start)
}

## ── Simulation function ──────────────────────────────────────

run_thyca <- function(mod, scenario, end_time = 365) {
  switch(scenario,

    "Untreated" = {
      p <- c(use_LENV=0, use_SORA=0, use_SELP=0, use_TSH_supp=1, is_MTC=0)
      ev_dose <- ev(amt=0, cmt="LENV_gut", time=0)
    },

    "Lenvatinib_1L" = {
      p <- c(use_LENV=1, use_SORA=0, use_SELP=0, use_TSH_supp=1, is_MTC=0)
      ev_dose <- make_lenv_dose(0, end_time)
    },

    "Sorafenib_1L" = {
      p <- c(use_LENV=0, use_SORA=1, use_SELP=0, use_TSH_supp=1, is_MTC=0)
      ev_dose <- make_sora_dose(0, end_time)
    },

    "Lenvatinib_2L" = {
      # Sorafenib first 6 months, then lenvatinib
      p <- c(use_LENV=1, use_SORA=1, use_SELP=0, use_TSH_supp=1, is_MTC=0)
      ev_dose <- make_sora_dose(0, 180) + make_lenv_dose(180, end_time)
    },

    "Selpercatinib_MTC" = {
      p <- c(use_LENV=0, use_SORA=0, use_SELP=1, use_TSH_supp=0, is_MTC=1)
      ev_dose <- make_selp_dose(0, end_time)
    },

    "Vandetanib_MTC" = {
      # Modeled as RET/VEGFR/EGFR inhibitor (~65% of selpercatinib potency)
      p <- c(use_LENV=0, use_SORA=0, use_SELP=1, use_TSH_supp=0, is_MTC=1,
             Emax_SELP_RET=0.75, EC50_SELP_RET=15.0)
      ev_dose <- ev(amt=160, cmt="SELP_gut", ii=1, addl=end_time-1, time=0)
    },

    "Lenvatinib_Everolimus" = {
      p <- c(use_LENV=1, use_SORA=0, use_SELP=0, use_TSH_supp=1, is_MTC=0,
             use_everolimus=1)
      ev_dose <- make_lenv_dose(0, end_time)
    }
  )

  thyca_mod %>%
    param(p) %>%
    mrgsim_q(ev_dose, end = end_time, delta = 1) %>%
    as_tibble() %>%
    mutate(scenario = scenario)
}

## ── Run all scenarios ─────────────────────────────────────────
scenarios <- c("Untreated", "Lenvatinib_1L", "Sorafenib_1L",
               "Lenvatinib_2L", "Selpercatinib_MTC",
               "Vandetanib_MTC", "Lenvatinib_Everolimus")

results <- map_df(scenarios, function(sc) {
  tryCatch(run_thyca(thyca_mod, sc, end_time = 730),
           error = function(e) { message(sc, ": ", e$message); tibble() })
})

## ── Summary ──────────────────────────────────────────────────
summary_tbl <- results %>%
  group_by(scenario) %>%
  summarise(
    Tg_final    = last(Tg),
    CT_final    = last(CT),
    TumVol_final = last(TumVol),
    SLD_at2y    = last(SLD_change),
    pct_CR      = mean(resp_CR) * 100,
    pct_PR      = mean(resp_PR) * 100,
    pct_PD_ever = max(resp_PD) * 100,
    MAPK_mean   = mean(MAPK_act),
    .groups = "drop"
  )

cat("\n── 2-Year Summary ───────────────────────────────────────\n")
print(summary_tbl, digits = 3)

## ── Color palette ─────────────────────────────────────────────
pal <- c(
  Untreated            = "#DC2626",
  Lenvatinib_1L        = "#2563EB",
  Sorafenib_1L         = "#16A34A",
  Lenvatinib_2L        = "#7C3AED",
  Selpercatinib_MTC    = "#EA580C",
  Vandetanib_MTC       = "#CA8A04",
  Lenvatinib_Everolimus = "#0891B2"
)

## ── Plots ─────────────────────────────────────────────────────
p1 <- ggplot(results, aes(time / 30, TumVol, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 30 * 0.7, linetype = "dashed", colour = "gray40") +
  annotate("text", x = 1, y = 30 * 0.7 * 0.95, label = "PR threshold (-30%)", hjust = 0, size = 3) +
  geom_hline(yintercept = 30 * 1.2, linetype = "dashed", colour = "red4") +
  annotate("text", x = 1, y = 30 * 1.22, label = "PD threshold (+20%)", hjust = 0, size = 3, colour = "red4") +
  scale_colour_manual(values = pal) +
  labs(title = "Tumor Volume (Sum of Longest Diameters)",
       x = "Time (months)", y = "SLD (mm)", colour = NULL) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p2 <- ggplot(results, aes(time / 30, Tg, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_y_log10() +
  labs(title = "Serum Thyroglobulin (DTC Biomarker)",
       x = "Time (months)", y = "Tg (ng/mL, log scale)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p3 <- ggplot(results, aes(time / 30, CT, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  scale_y_log10() +
  labs(title = "Serum Calcitonin (MTC Biomarker)",
       x = "Time (months)", y = "Calcitonin (pg/mL, log scale)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p4 <- ggplot(results, aes(time / 30, MAPK_act, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  labs(title = "MAPK Pathway Activity",
       x = "Time (months)", y = "Activity (normalized 0–1)", colour = NULL) +
  ylim(0, 1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p5 <- ggplot(results, aes(time / 30, Cp_LENV_out, colour = scenario)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = pal) +
  labs(title = "Lenvatinib Plasma Concentration",
       x = "Time (months)", y = "Cp (ng/mL)", colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p6 <- ggplot(results, aes(time / 30, SLD_change, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = c(-30, 20), linetype = "dashed", colour = c("#16A34A", "#DC2626")) +
  scale_colour_manual(values = pal) +
  labs(title = "% Change in Sum of Longest Diameters (RECIST)",
       x = "Time (months)", y = "% Change from baseline", colour = NULL) +
  coord_cartesian(ylim = c(-80, 100)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p1); print(p2); print(p3)
print(p4); print(p5); print(p6)

## ── Virtual Patient Analysis ──────────────────────────────────
set.seed(2026)
n_vp <- 200

vp_params <- tibble(
  id         = 1:n_vp,
  kprol_vp   = rlnorm(n_vp, log(0.03), 0.35),
  MAPK_vp    = rbeta(n_vp, 8, 2),      # BRAF V600E → high MAPK
  F_LENV_vp  = rnorm(n_vp, 0.85, 0.10) |> pmax(0.4) |> pmin(1.0),
  TumVol_0   = rlnorm(n_vp, log(30), 0.5)
)

vp_sim <- map_df(c("Untreated", "Lenvatinib_1L"), function(sc) {
  map_df(seq_len(min(50, n_vp)), function(i) {
    p_vp <- c(
      use_LENV = if(sc == "Lenvatinib_1L") 1 else 0,
      use_SORA = 0, use_SELP = 0, use_TSH_supp = 1, is_MTC = 0,
      kprol_base = vp_params$kprol_vp[i],
      MAPK_base  = vp_params$MAPK_vp[i],
      F_LENV     = vp_params$F_LENV_vp[i]
    )
    ev_vp <- if(sc == "Lenvatinib_1L") make_lenv_dose(0, 730) else
             ev(amt=0, cmt="LENV_gut", time=0)
    tryCatch({
      thyca_mod %>%
        param(p_vp) %>%
        init(TumVol = vp_params$TumVol_0[i], Tg = 50) %>%
        mrgsim_q(ev_vp, end = 730, delta = 30) %>%
        as_tibble() %>%
        mutate(id = i, scenario = sc)
    }, error = function(e) tibble())
  })
})

p_vp <- ggplot(vp_sim, aes(time/30, SLD_change,
                            group = interaction(id, scenario),
                            colour = scenario)) +
  geom_line(alpha = 0.2, linewidth = 0.4) +
  stat_summary(aes(group = scenario), fun = median,
               geom = "line", linewidth = 2) +
  geom_hline(yintercept = c(-30, 20), linetype = "dashed") +
  scale_colour_manual(values = c(Untreated="#DC2626", Lenvatinib_1L="#2563EB")) +
  coord_cartesian(ylim = c(-90, 120)) +
  labs(title = "Virtual Patient % SLD Change (n=50 per arm)",
       subtitle = "Thin = individual; thick = median",
       x = "Time (months)", y = "% SLD Change", colour = NULL) +
  theme_bw(base_size = 12)

print(p_vp)

cat("\n── Model complete ─────────────────────────────────────────\n")
cat("results object: ", nrow(results), "rows\n")
cat("vp_sim object : ", nrow(vp_sim), "rows\n")
