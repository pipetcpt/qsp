##############################################################################
# Chronic Hypothyroidism QSP Model — mrgsolve ODE Implementation
# 만성 갑상선 기능 저하증 정량적 시스템 약리학 모델
#
# Author  : Claude Code Routine (CCR)
# Date    : 2026-06-20
# Disease : Chronic Hypothyroidism (원발성 / 중추성 / 약물 유발)
#
# Model Scope
# -----------
# 1. Hypothalamic-Pituitary-Thyroid (HPT) axis with Hill-function feedback
# 2. T4/T3/rT3 distribution and peripheral deiodinase metabolism
# 3. Levothyroxine (LT4) 2-compartment oral PK (MW = 776.87 g/mol)
# 4. Liothyronine (LT3) 1-compartment oral PK  (MW = 650.97 g/mol)
# 5. PD biomarkers: HR, LDL, BMR, symptom score, BMD change
#
# References
# ----------
# Berberich AJ et al. (2017) Thyroid 27(12):1454-1462
# Jonklaas J et al. (2014) Thyroid 24(12):1670-1751
# Dietrich JW et al. (2016) Front Endocrinol 7:29
# Larsen PR et al. (2012) Williams Textbook of Endocrinology
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ============================================================
# MODEL CODE
# ============================================================
code <- '
$PLUGIN Rcpp

$PARAM @annotated
// --- HPT Axis Physiological Parameters ---
TRH0        : 5.0    : TRH baseline concentration (pmol/L)
TSH0        : 2.0    : TSH baseline (mIU/L)
TT40        : 100.0  : Total T4 baseline (nmol/L)
TT30        : 1.8    : Total T3 baseline (nmol/L)
rT30        : 0.45   : Reverse T3 baseline (nmol/L)

// --- TRH Kinetics (Hypothalamus) ---
ksyn_TRH    : 5.0    : TRH synthesis rate (pmol/L/h)
kdeg_TRH    : 1.04   : TRH degradation rate (h-1, t1/2 ≈ 40 min)
EC50_TRH    : 5.4    : fT3 EC50 for TRH inhibition (pmol/L)
n_TRH       : 2.0    : Hill coefficient, TRH-T3 feedback

// --- TSH Kinetics (Anterior Pituitary) ---
ksyn_TSH    : 2.08   : TSH synthesis rate (mIU/L/h)
kdeg_TSH    : 0.693  : TSH degradation rate (h-1, t1/2 ≈ 1 h)
EC50_T4_fb  : 20.0   : fT4 EC50 for TSH inhibitory feedback (pmol/L)
EC50_T3_fb  : 5.4    : fT3 EC50 for TSH inhibitory feedback (pmol/L)
n_fb        : 2.0    : Hill coefficient, T4/T3 → TSH feedback
n_TRH_stim  : 1.5    : Hill coefficient, TRH stimulation of TSH

// --- Thyroid Hormone Synthesis & Secretion ---
ksec_T4     : 1.4    : T4 thyroid secretion rate at baseline TSH (nmol/L/h)
ksec_T3     : 0.052  : T3 thyroid secretion rate at baseline TSH (nmol/L/h)
thyroid_cap : 1.0    : Thyroid functional capacity (0=ablated, 1=normal)

// --- T4/T3/rT3 Kinetics ---
kconv_T4T3  : 0.0043 : T4 to T3 peripheral conversion (h-1) [D1/D2]
kconv_T4rT3 : 0.0019 : T4 to rT3 conversion (h-1)          [D3]
kel_T4      : 0.0041 : T4 elimination rate (h-1, t1/2 = 168h ≈ 7d)
kel_T3      : 0.029  : T3 elimination rate (h-1, t1/2 = 24h)
kel_rT3     : 0.145  : rT3 elimination rate (h-1, t1/2 ≈ 4.8h)
ff_T4       : 0.0002 : Free fraction of total T4 (unitless, 0.02%)
ff_T3       : 0.003  : Free fraction of total T3 (unitless, 0.30%)

// --- Levothyroxine (LT4) PK Parameters ---
Ka_LT4      : 0.35   : LT4 absorption rate constant (h-1, Tmax ≈ 3 h)
F_LT4       : 0.70   : LT4 bioavailability (fasted; ranges 0.65-0.80)
V1_LT4      : 15.0   : LT4 central volume (L) [effective plasma]
Vss_LT4     : 700.0  : LT4 Vss (L, ~10 L/kg, 70-kg patient)
k12_LT4     : 0.020  : LT4 central-to-peripheral rate (h-1)
k21_LT4     : 0.004  : LT4 peripheral-to-central rate (h-1)
CL_LT4      : 0.50   : LT4 systemic clearance (L/h, t1/2 ≈ 168h in Vss)
MW_T4       : 776.87 : Molecular weight of T4 (g/mol)

// --- Liothyronine (LT3) PK Parameters ---
Ka_LT3      : 0.80   : LT3 absorption rate (h-1, Tmax ≈ 1 h)
F_LT3       : 0.95   : LT3 bioavailability
V1_LT3      : 15.0   : LT3 central volume (L)
CL_LT3      : 2.5    : LT3 clearance (L/h, t1/2 ≈ 24h)
MW_T3       : 650.97 : Molecular weight of T3 (g/mol)

// --- PD Parameters ---
// Heart Rate
HR0         : 70.0   : Baseline heart rate (bpm)
Emax_HR     : 25.0   : Max HR change attributable to T3 (bpm)
EC50_HR     : 5.4    : fT3 EC50 for HR (pmol/L)
kout_HR     : 0.10   : HR response rate constant (h-1)
// LDL Cholesterol
LDL0        : 3.2    : Baseline LDL (mmol/L)
Emax_LDL    : 0.60   : Max LDL fold-increase when T3 → 0
EC50_LDL    : 5.4    : fT3 EC50 for LDL effect (pmol/L)
kout_LDL    : 0.005  : LDL response rate constant (h-1)
// BMR
BMR0        : 100.0  : Baseline BMR (% of normal)
Emax_BMR    : 0.35   : Max fractional BMR decrease when T3 → 0
EC50_BMR    : 5.4    : fT3 EC50 for BMR effect (pmol/L)
kout_BMR    : 0.008  : BMR response rate constant (h-1)
// Symptom Score (0 = asymptomatic, 100 = maximally symptomatic)
Emax_Sym    : 95.0   : Max hypothyroid symptom score
EC50_Sym    : 5.4    : fT3 EC50 for symptom (pmol/L)
kout_Sym    : 0.02   : Symptom response rate constant (h-1)
// BMD (% change from baseline, tracked over months-years)
kout_BMD    : 0.0001 : BMD change rate constant (h-1)

$CMT @annotated
A_TRH       : Hypothalamic TRH (pmol/L)
A_TSH       : Pituitary TSH (mIU/L)
A_TT4       : Total T4 (nmol/L)
A_TT3       : Total T3 (nmol/L)
A_rT3       : Reverse T3 (nmol/L)
A_LT4_gut   : LT4 gut absorption compartment (ug)
A_LT4_c     : LT4 central compartment (ug)
A_LT4_p     : LT4 peripheral compartment (ug)
A_LT3_gut   : LT3 gut absorption compartment (ug)
A_LT3_c     : LT3 central compartment (ug)
Eff_HR      : Heart rate (bpm)
Eff_LDL     : LDL cholesterol (mmol/L)
Eff_BMR     : Basal metabolic rate (% normal)
Eff_Sym     : Hypothyroid symptom score (0-100)
Eff_BMD     : Bone mineral density change from baseline (%)

$INIT
A_TRH    = 5.0
A_TSH    = 2.0
A_TT4    = 100.0
A_TT3    = 1.8
A_rT3    = 0.45
A_LT4_gut = 0
A_LT4_c  = 0
A_LT4_p  = 0
A_LT3_gut = 0
A_LT3_c  = 0
Eff_HR   = 70.0
Eff_LDL  = 3.2
Eff_BMR  = 100.0
Eff_Sym  = 5.0
Eff_BMD  = 0.0

$MAIN
// ---- Derived concentrations (free fractions) ----
double fT4 = A_TT4 * ff_T4 * 1000.0;   // pmol/L
double fT3 = A_TT3 * ff_T3 * 1000.0;   // pmol/L
double fT3_base = TT30 * ff_T3 * 1000.0; // baseline fT3 in pmol/L

// ---- Thyroid secretion (linear with TSH relative to baseline) ----
double TSH_rel  = A_TSH / TSH0;  // TSH ratio (1 = normal)
double T4_sec   = ksec_T4 * TSH_rel * thyroid_cap;
double T3_sec   = ksec_T3 * TSH_rel * thyroid_cap;

// ---- LT4 contribution to TT4 pool (ug → nmol/L in V1) ----
// [A_LT4_c ug] / [V1_LT4 L] / [MW g/mol] * 1e6 = nmol/L
// But absorption gives ug absorbed; central compartment feeds TT4 pool via
// conversion to endogenous-equivalent T4 in full distribution volume.
// Simplified: contribution rate = absorption_rate_into_central / Vss_LT4 / MW * 1e6
// In ODE below, we add Ka*A_LT4_gut*F_LT4 to a separate track

// ---- LT3 contribution to TT3 (similar scaling) ----

// ---- PD equilibrium targets ----
// Heart rate: T3-mediated, increases with T3
double HR_target  = HR0  * (fT3 / fT3_base);
if (HR_target < 35.0)  HR_target = 35.0;
if (HR_target > 130.0) HR_target = 130.0;

// LDL: inversely related to T3 (LDL rises when T3 low)
double LDL_target = LDL0 * (1.0 + Emax_LDL * (1.0 - fT3 / (EC50_LDL + fT3)) /
                              (1.0 - fT3_base / (EC50_LDL + fT3_base)));
if (LDL_target < 0.5) LDL_target = 0.5;

// BMR: proportional to T3 relative to baseline
double BMR_target = BMR0 * (fT3 / fT3_base);
if (BMR_target < 50.0)  BMR_target = 50.0;
if (BMR_target > 130.0) BMR_target = 130.0;

// Symptom score: high when T3 low
double fT3_safe = (fT3 < 0.01) ? 0.01 : fT3;
double Sym_target = Emax_Sym * (1.0 - fT3_safe / (EC50_Sym + fT3_safe)) /
                    (1.0 - fT3_base / (EC50_Sym + fT3_base));
if (Sym_target < 0.0)  Sym_target = 0.0;
if (Sym_target > 100.0) Sym_target = 100.0;

// BMD: slight loss with T3 excess (over-treatment); slight stagnation with deficiency
// Net effect: over-treatment accelerates bone loss; under-treatment slows turnover
double BMD_target = -0.5 * (fT3 / fT3_base - 1.0) * 12.0;  // %/year projection
if (BMD_target < -8.0) BMD_target = -8.0;
if (BMD_target > 2.0)  BMD_target = 2.0;

$ODE
// ---- TRH (hypothalamus) ----
// Inhibited by fT3 (main feedback) via Hill function
double Inh_TRH = pow(EC50_TRH, n_TRH) /
                 (pow(EC50_TRH, n_TRH) + pow(fT3, n_TRH));
dxdt_A_TRH = ksyn_TRH * Inh_TRH - kdeg_TRH * A_TRH;

// ---- TSH (anterior pituitary) ----
// Stimulated by TRH (Hill); inhibited by fT4 and fT3
double Stim_TRH = pow(A_TRH / TRH0, n_TRH_stim);
double Inh_T4_pit = pow(EC50_T4_fb, n_fb) /
                    (pow(EC50_T4_fb, n_fb) + pow(fT4, n_fb));
double Inh_T3_pit = pow(EC50_T3_fb, n_fb) /
                    (pow(EC50_T3_fb, n_fb) + pow(fT3, n_fb));
dxdt_A_TSH = ksyn_TSH * Stim_TRH * Inh_T4_pit * Inh_T3_pit - kdeg_TSH * A_TSH;

// ---- Total T4 ----
// Exogenous LT4: absorbed from gut at Ka*F; scaled to TT4 in Vss_LT4 (nmol/L)
double LT4_abs_rate = Ka_LT4 * A_LT4_gut * F_LT4;              // ug/h
double LT4_to_TT4   = LT4_abs_rate / Vss_LT4 / MW_T4 * 1.0e6; // nmol/L/h

dxdt_A_TT4 = T4_sec
             + LT4_to_TT4
             - kconv_T4T3  * A_TT4
             - kconv_T4rT3 * A_TT4
             - kel_T4       * A_TT4;

// ---- Total T3 ----
double LT3_abs_rate = Ka_LT3 * A_LT3_gut * F_LT3;              // ug/h
double LT3_to_TT3   = LT3_abs_rate / V1_LT3 / MW_T3 * 1.0e6;  // nmol/L/h

dxdt_A_TT3 = T3_sec
             + kconv_T4T3 * A_TT4
             + LT3_to_TT3
             - kel_T3 * A_TT3;

// ---- Reverse T3 ----
dxdt_A_rT3 = kconv_T4rT3 * A_TT4 - kel_rT3 * A_rT3;

// ---- LT4 gut compartment ----
dxdt_A_LT4_gut = -Ka_LT4 * A_LT4_gut;

// ---- LT4 central compartment (ug) ----
double k10_LT4 = CL_LT4 / V1_LT4;
dxdt_A_LT4_c = Ka_LT4 * A_LT4_gut * F_LT4
               - (k10_LT4 + k12_LT4) * A_LT4_c
               + k21_LT4 * A_LT4_p;

// ---- LT4 peripheral compartment (ug) ----
dxdt_A_LT4_p = k12_LT4 * A_LT4_c - k21_LT4 * A_LT4_p;

// ---- LT3 gut compartment ----
dxdt_A_LT3_gut = -Ka_LT3 * A_LT3_gut;

// ---- LT3 central compartment (ug) ----
double k10_LT3 = CL_LT3 / V1_LT3;
dxdt_A_LT3_c = Ka_LT3 * A_LT3_gut * F_LT3 - k10_LT3 * A_LT3_c;

// ---- PD effect compartments (indirect response) ----
dxdt_Eff_HR  = kout_HR  * (HR_target  - Eff_HR);
dxdt_Eff_LDL = kout_LDL * (LDL_target - Eff_LDL);
dxdt_Eff_BMR = kout_BMR * (BMR_target - Eff_BMR);
dxdt_Eff_Sym = kout_Sym * (Sym_target - Eff_Sym);
dxdt_Eff_BMD = kout_BMD * (BMD_target - Eff_BMD);

$TABLE
// Output derived quantities
double fT4_pmol  = A_TT4 * ff_T4 * 1000.0;   // pmol/L
double fT3_pmol  = A_TT3 * ff_T3 * 1000.0;   // pmol/L
double fT4_ngdL  = fT4_pmol * 0.0777;         // ng/dL (1 pmol/L T4 = 0.0777 ng/dL)
double fT3_pgmL  = fT3_pmol * 0.651;          // pg/mL (1 pmol/L T3 = 0.651 pg/mL)
double TT4_ugdL  = A_TT4 * 0.0777;            // μg/dL
double TT3_ngdL  = A_TT3 * 65.1;              // ng/dL
double LT4_conc  = A_LT4_c / V1_LT4;          // ug/L in central compartment
double WeightEst = (100.0 - Eff_BMR) * 0.15;  // rough body weight gain estimate (kg)

$CAPTURE fT4_pmol fT3_pmol fT4_ngdL fT3_pgmL TT4_ugdL TT3_ngdL
         A_TSH A_TT4 A_TT3 A_rT3 LT4_conc
         Eff_HR Eff_LDL Eff_BMR Eff_Sym Eff_BMD WeightEst
'

# Compile model
mod <- mcode("hypo_qsp", code)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Dosing regimen: daily oral LT4 (midnight) and optional LT3 (morning)
make_ev <- function(LT4_dose_ug = 0, LT3_dose_ug = 0,
                    n_days = 365, LT3_times = c(8, 20)) {
  evs <- list()
  if (LT4_dose_ug > 0) {
    evs[[1]] <- ev(amt = LT4_dose_ug, cmt = "A_LT4_gut",
                   ii = 24, addl = n_days - 1, time = 0)
  }
  if (LT3_dose_ug > 0) {
    for (tt in LT3_times) {
      evs[[length(evs) + 1]] <- ev(amt = LT3_dose_ug / length(LT3_times),
                                   cmt = "A_LT3_gut",
                                   ii = 24, addl = n_days - 1, time = tt)
    }
  }
  if (length(evs) == 0) return(ev(amt = 0, cmt = 1, time = 1e6))
  do.call(c, evs)
}

run_sim <- function(mod, ev_obj, thyroid_cap = 1.0, t_max = 8760, delta = 4) {
  mrgsim(mod, ev = ev_obj, end = t_max, delta = delta,
         param = list(thyroid_cap = thyroid_cap)) %>%
    as.data.frame()
}

# ============================================================
# SCENARIO DEFINITIONS (7 SCENARIOS)
# ============================================================

# Total simulation: 2 years (17520 h); display weeks in output
T_MAX <- 8760    # 1 year (hours)
DELTA <- 6       # output every 6 hours

# 1. Normal euthyroid (no disease, no treatment)
scenario_labels <- c(
  "1. 정상 갑상선 (Euthyroid)",
  "2. 미치료 갑저 (Untreated Hypothyroidism)",
  "3. LT4 표준 치료 (LT4 100 μg/day)",
  "4. LT4 저용량 (준임상적 갑저, LT4 50 μg/day)",
  "5. LT4 과치료 (TSH 억제, LT4 175 μg/day)",
  "6. T4+T3 병용 (Combination LT4 100+LT3 10 μg/day)",
  "7. 갑상선 전절제 후 LT4 (Post-Thyroidectomy LT4 125 μg/day)"
)

scenarios <- list(
  list(LT4 = 0,   LT3 = 0,  cap = 1.0, label = scenario_labels[1]),
  list(LT4 = 0,   LT3 = 0,  cap = 0.0, label = scenario_labels[2]),
  list(LT4 = 100, LT3 = 0,  cap = 0.0, label = scenario_labels[3]),
  list(LT4 = 50,  LT3 = 0,  cap = 0.3, label = scenario_labels[4]),
  list(LT4 = 175, LT3 = 0,  cap = 0.0, label = scenario_labels[5]),
  list(LT4 = 100, LT3 = 10, cap = 0.0, label = scenario_labels[6]),
  list(LT4 = 125, LT3 = 0,  cap = 0.0, label = scenario_labels[7])
)

cat("Running 7 QSP scenarios for Chronic Hypothyroidism...\n")

results <- lapply(seq_along(scenarios), function(i) {
  sc <- scenarios[[i]]
  ev_obj <- make_ev(LT4_dose_ug = sc$LT4, LT3_dose_ug = sc$LT3,
                    n_days = T_MAX / 24)
  df <- run_sim(mod, ev_obj, thyroid_cap = sc$cap,
                t_max = T_MAX, delta = DELTA)
  df$scenario  <- sc$label
  df$scenario_id <- i
  df
})

df_all <- bind_rows(results)
df_all$time_weeks <- df_all$time / 168  # convert hours to weeks

# ============================================================
# STEADY-STATE SUMMARY TABLE
# ============================================================
ss_summary <- df_all %>%
  filter(time > T_MAX * 0.85) %>%   # last 15% = approximate SS
  group_by(scenario) %>%
  summarise(
    TSH_mIU_L    = round(mean(A_TSH),   2),
    FT4_pmol_L   = round(mean(fT4_pmol), 1),
    FT3_pmol_L   = round(mean(fT3_pmol), 1),
    TT4_nmol_L   = round(mean(A_TT4),   1),
    TT3_nmol_L   = round(mean(A_TT3),   2),
    HR_bpm       = round(mean(Eff_HR),  0),
    LDL_mmol_L   = round(mean(Eff_LDL), 2),
    BMR_pct      = round(mean(Eff_BMR), 1),
    Symptom_score = round(mean(Eff_Sym), 1),
    .groups = "drop"
  )

cat("\n========== STEADY-STATE SUMMARY (Year 1) ==========\n")
print(ss_summary, n = 20, width = 120)

# Normal reference ranges
cat("\n--- Reference Ranges ---\n")
cat("TSH:  0.5-4.5 mIU/L   |  FT4: 12-22 pmol/L  |  FT3: 3.5-6.5 pmol/L\n")
cat("TT4: 58-161 nmol/L    |  TT3: 1.2-2.7 nmol/L\n")
cat("HR:  60-100 bpm       |  LDL: <3.4 mmol/L   |  BMR: ~100%\n")

# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

cols7 <- c("#1b9e77","#d95f02","#7570b3","#e7298a",
           "#66a61e","#e6ab02","#a6761d")

plot_var <- function(df, var, ylab, title, ref_lo = NA, ref_hi = NA) {
  p <- ggplot(df, aes(x = time_weeks, y = .data[[var]],
                      color = scenario, group = scenario)) +
    geom_line(size = 0.8, alpha = 0.85) +
    scale_color_manual(values = cols7) +
    labs(x = "시간 (주, Weeks)", y = ylab, title = title,
         color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7),
          plot.title = element_text(face = "bold", size = 10))
  if (!is.na(ref_lo))
    p <- p + geom_hline(yintercept = ref_lo, linetype = "dashed",
                        color = "grey40", alpha = 0.6)
  if (!is.na(ref_hi))
    p <- p + geom_hline(yintercept = ref_hi, linetype = "dashed",
                        color = "grey40", alpha = 0.6)
  p
}

# HPT axis plots
p_tsh <- plot_var(df_all, "A_TSH",   "TSH (mIU/L)",    "혈청 TSH", 0.5, 4.5)
p_ft4 <- plot_var(df_all, "fT4_pmol","Free T4 (pmol/L)","유리 T4", 12, 22)
p_ft3 <- plot_var(df_all, "fT3_pmol","Free T3 (pmol/L)","유리 T3", 3.5, 6.5)
p_rt3 <- plot_var(df_all, "A_rT3",   "rT3 (nmol/L)",   "역-T3 (Reverse T3)")

# PD effect plots
p_hr  <- plot_var(df_all, "Eff_HR",  "Heart Rate (bpm)","심박수", 60, 100)
p_ldl <- plot_var(df_all, "Eff_LDL", "LDL (mmol/L)",   "LDL 콜레스테롤", NA, 3.4)
p_bmr <- plot_var(df_all, "Eff_BMR", "BMR (%)",        "기초대사율 (BMR)", 85, 115)
p_sym <- plot_var(df_all, "Eff_Sym", "증상 점수 (0-100)","갑저 증상 점수")

# Combine plots
hpt_panel <- (p_tsh | p_ft4) / (p_ft3 | p_rt3) +
  plot_annotation(title = "HPT 축 역학 — 7가지 치료 시나리오",
                  theme = theme(plot.title = element_text(face="bold")))

pd_panel  <- (p_hr | p_ldl) / (p_bmr | p_sym) +
  plot_annotation(title = "PD 효과 — 7가지 치료 시나리오",
                  theme = theme(plot.title = element_text(face="bold")))

# Save plots
ggsave("hypo_hpt_panel.png", hpt_panel, width = 14, height = 9, dpi = 150)
ggsave("hypo_pd_panel.png",  pd_panel,  width = 14, height = 9, dpi = 150)

cat("\n[Saved] hypo_hpt_panel.png, hypo_pd_panel.png\n")

# ============================================================
# DOSE-RESPONSE ANALYSIS
# ============================================================
cat("\n--- Dose-response: LT4 dose vs. steady-state TSH (thyroid_cap = 0) ---\n")
lt4_doses <- c(25, 50, 75, 100, 125, 150, 175, 200)
dr_results <- lapply(lt4_doses, function(dose) {
  ev_obj <- make_ev(LT4_dose_ug = dose, n_days = 365)
  df_ss  <- mrgsim(mod, ev = ev_obj, end = 8760, delta = 24,
                   param = list(thyroid_cap = 0)) %>%
    as.data.frame() %>%
    tail(14) %>%   # last 2 weeks
    summarise(
      LT4_dose     = dose,
      TSH_ss       = round(mean(A_TSH),    2),
      FT4_ss       = round(mean(fT4_pmol), 1),
      FT3_ss       = round(mean(fT3_pmol), 1),
      HR_ss        = round(mean(Eff_HR),   0),
      LDL_ss       = round(mean(Eff_LDL),  2),
      .groups      = "drop"
    )
  df_ss
})
dr_df <- bind_rows(dr_results)
cat("\nDose-Response Table (Post-Thyroidectomy, thyroid_cap = 0):\n")
print(dr_df)

# ============================================================
# SENSITIVITY ANALYSIS: bioavailability effect
# ============================================================
cat("\n--- Sensitivity: F_LT4 effect on TSH at 100 ug/day ---\n")
F_vals <- seq(0.50, 0.90, by = 0.05)
sa_F <- lapply(F_vals, function(f) {
  ev_obj <- make_ev(LT4_dose_ug = 100, n_days = 365)
  df_ss  <- mrgsim(mod, ev = ev_obj, end = 8760, delta = 24,
                   param = list(thyroid_cap = 0, F_LT4 = f)) %>%
    as.data.frame() %>% tail(14) %>%
    summarise(F_LT4 = f,
              TSH   = round(mean(A_TSH), 2),
              FT4   = round(mean(fT4_pmol), 1))
})
print(bind_rows(sa_F))

cat("\n[Done] Chronic Hypothyroidism QSP model simulation complete.\n")
cat("Files created: hypo_hpt_panel.png, hypo_pd_panel.png\n")
