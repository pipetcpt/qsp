##############################################################################
# Pernicious Anemia (악성 빈혈) — QSP Model
# mrgsolve ODE-based PK/PD Simulation
#
# Disease: Pernicious anemia (autoimmune gastric atrophy → IF deficiency →
#          cobalamin malabsorption → megaloblastic anemia + SCD)
#
# Drug targets modeled:
#   - Cyanocobalamin IM (loading + maintenance)
#   - Hydroxocobalamin IM (higher binding affinity)
#   - High-dose oral cyanocobalamin (passive absorption)
#   - Methylcobalamin oral
#
# Key compartments (15 ODEs):
#   1.  DEPOT       — IM injection depot (µg)
#   2.  ORAL        — GI-tract oral B12 (µg)
#   3.  IF_POOL     — Functional intrinsic factor (relative, 0-1)
#   4.  CBA_PORTAL  — Portal-absorbed cobalamin (µg) [IF-mediated + passive]
#   5.  PLASMA      — Plasma total cobalamin (pg/mL)
#   6.  HOLOTC      — HoloTranscobalamin II (active fraction, pmol/L)
#   7.  LIVER       — Hepatic cobalamin stores (µg)
#   8.  BONE_MARROW — Bone marrow cobalamin (µg)
#   9.  NERVE       — Nervous system cobalamin (µg)
#  10.  AUTO_AB     — Autoimmune antibody burden (anti-IF Ab, relative 0-1)
#  11.  PARIETAL    — Parietal cell function index (0-1)
#  12.  HGB         — Hemoglobin (g/dL)
#  13.  MCV         — Mean corpuscular volume (fL)
#  14.  RETIC       — Reticulocyte count (× 10⁹/L)
#  15.  NEURO       — Neurological disability score (0-10; higher = worse)
#
# Treatment scenarios:
#   S1: Untreated PA (disease progression only)
#   S2: IM loading (1000 µg/day × 7d) then monthly maintenance
#   S3: IM maintenance only (1000 µg/month)
#   S4: High-dose oral (1000 µg/day, passive 1%)
#   S5: Aggressive loading (2000 µg/day × 14d) then biweekly
#
# Calibration references:
#   - Carmel 2008 (Am J Clin Nutr): B12 absorption physiology
#   - Stabler 2013 (NEJM): B12 deficiency review
#   - Andres et al. 2004 (QJM): oral vs IM efficacy
#   - Oh & Brown 2003 (Am Fam Physician): PA clinical management
#   - Wolffenbuttel et al. 2019 (BMJ): holotranscobalamin cutoffs
#
# Author: Claude Code Routine — QSP Disease Library
# Date:   2026-06-20
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ============================================================
# MODEL SPECIFICATION
# ============================================================
pa_model_code <- '
$PROB Pernicious Anemia QSP Model — Cobalamin PK/PD

$PARAM @annotated
// --- PK Parameters ---
KA_IM     : 0.40  : IM absorption rate constant (h⁻¹)
KA_ORAL   : 0.08  : Oral gastric emptying rate (h⁻¹)
F_PASSIVE  : 0.01  : Passive oral absorption fraction (IF-independent, ~1%)
KIF_ABS   : 2.50  : IF-mediated absorption rate constant (h⁻¹, ileum)
KIF_MAX   : 1.50  : Maximum IF-mediated absorption capacity scaling
CL_PLASMA  : 28.0  : Plasma cobalamin clearance to tissues (L/h)
V_PLASMA   : 6.0   : Central plasma volume (L) — for free B12
K12_LIVER  : 0.30  : Rate of liver uptake from plasma (h⁻¹)
K21_LIVER  : 0.015 : Rate of liver release to plasma (h⁻¹)
K12_BM     : 0.08  : Rate of bone marrow uptake from plasma (h⁻¹)
K21_BM     : 0.025 : Rate of bone marrow release (h⁻¹)
K12_NERVE  : 0.04  : Rate of nerve uptake from plasma (h⁻¹)
K21_NERVE  : 0.002 : Rate of nerve release (h⁻¹)
KREN       : 0.25  : Renal excretion rate of free plasma B12 (h⁻¹)
ENTERO_F   : 0.005 : Enterohepatic recycling fraction (/h)

// --- Immune/Autoimmune Parameters ---
KAB_PROD  : 0.0005 : Auto-antibody production rate (/h)
KAB_DECAY : 0.0002 : Auto-antibody natural decay (/h)
AB0        : 0.80  : Baseline auto-antibody burden (PA established, 0-1)
KIF_DECAY : 0.003  : Parietal cell destruction rate by autoimmunity (/h)
KIF_REGEN : 0.0001 : Residual parietal cell regeneration (/h)

// --- PD Parameters: Hematopoiesis ---
HGB0       : 7.5   : Baseline hemoglobin in PA (g/dL)
HGB_MAX    : 14.5  : Maximum hemoglobin (g/dL)
EC50_HGB   : 300.0 : B12 plasma conc for half-max Hb response (pg/mL)
EMAX_HGB   : 1.0   : Maximum Emax for Hb restoration
KOUT_HGB   : 0.004 : Hb equilibration rate (/h, ~120-day RBC half-life)
MCV0       : 115.0 : Baseline MCV in PA (fL)
MCV_NL     : 88.0  : Normal MCV target (fL)
K_MCV      : 0.005 : MCV normalization rate (/h)
EC50_MCV   : 200.0 : B12 for MCV normalization (pg/mL)
RETIC_BASE : 30.0  : Baseline reticulocyte count (×10⁹/L, PA)
RETIC_MAX  : 400.0 : Peak reticulocyte count (crisis response)
RETIC_EC50 : 250.0 : B12 for half-max retic response (pg/mL)
K_RETIC    : 0.05  : Retic response rate constant (/h)

// --- PD Parameters: Neurological ---
NEURO0     : 4.0   : Baseline neurological score (0-10, PA)
NEURO_MAX  : 8.0   : Maximum neurological disability (untreated)
NEURO_MIN  : 1.0   : Minimum achievable score with treatment
KNEURO_PROG: 0.0008: Neurological progression rate (without B12, /h)
KNEURO_RECOV: 0.003: Neurological recovery rate (with B12, /h)
EC50_NEURO : 400.0 : Plasma B12 for neuro protection (pg/mL)

// --- Biomarker Parameters ---
MMA_BASE   : 3.5   : Baseline MMA (µmol/L) in untreated PA
MMA_NORM   : 0.15  : Normal MMA (µmol/L)
K_MMA      : 0.05  : MMA normalization rate (/h)
HCY_BASE   : 28.0  : Baseline homocysteine (µmol/L) in PA
HCY_NORM   : 9.0   : Normal homocysteine (µmol/L)
K_HCY      : 0.04  : Homocysteine normalization rate (/h)

// --- Initial States ---
PARIETAL0  : 0.20  : Residual parietal cell function in PA (20%)
LIVER0     : 1000  : Hepatic B12 stores at PA diagnosis (µg, depleted)
PLASMA0    : 80.0  : Plasma B12 at PA diagnosis (pg/mL, deficient)

$CMT @annotated
DEPOT     : IM injection depot (µg)
ORAL_GI   : Oral B12 in GI tract (µg)
IF_POOL   : Intrinsic factor functional capacity (0-1)
CBA_PORT  : Portal-absorbed cobalamin (µg)
PLASMA    : Plasma total cobalamin (pg/mL × V_PLASMA)
HOLOTC    : HoloTranscobalamin II (active B12, pmol/L)
LIVER     : Hepatic cobalamin stores (µg)
BONE_MRW  : Bone marrow cobalamin pool (µg)
NERVE     : Nervous system cobalamin pool (µg)
AUTO_AB   : Autoimmune antibody burden (0-1)
PARIETAL  : Parietal cell function index (0-1)
HGB       : Hemoglobin (g/dL)
MCV       : Mean corpuscular volume (fL)
RETIC     : Reticulocyte count (×10⁹/L)
NEURO     : Neurological disability score (0-10)

$INIT
DEPOT    = 0
ORAL_GI  = 0
IF_POOL  = 0.20
CBA_PORT = 0
PLASMA   = 80.0 * 6.0    // pg/mL × volume = pg (treat as amount)
HOLOTC   = 12.0           // pmol/L — severely depleted
LIVER    = 1000.0         // µg (depleted stores)
BONE_MRW = 0.8            // µg (very low)
NERVE    = 0.5            // µg (depleted)
AUTO_AB  = 0.80           // high autoimmune burden
PARIETAL = 0.20           // 20% residual function
HGB      = 7.5            // g/dL (severe PA anemia)
MCV      = 115.0          // fL (macrocytic)
RETIC    = 30.0           // ×10⁹/L (low in PA)
NEURO    = 4.0            // moderate neurological disability

$ODE
// ------------------------------------------------------------
// Derived plasma concentration (pg/mL)
double Cp = PLASMA / V_PLASMA;  // pg/mL

// IF-mediated absorption efficiency (modulated by auto-Ab)
double IF_eff = IF_POOL * (1.0 - AUTO_AB * 0.9);  // Ab blocks IF ~90% at full Ab burden
double IF_abs_rate = KIF_ABS * IF_eff * KIF_MAX;

// ------------------------------------------------------------
// 1. DEPOT — IM injection depot
dxdt_DEPOT = -KA_IM * DEPOT;

// 2. ORAL_GI — GI tract B12
double IF_mediated = (ORAL_GI > 0) ? IF_abs_rate * ORAL_GI / (ORAL_GI + 0.1) : 0.0;
double passive_abs  = KA_ORAL * F_PASSIVE * ORAL_GI;
dxdt_ORAL_GI = -KA_ORAL * ORAL_GI;  // gastric emptying

// 3. IF_POOL — Parietal cell function (drives IF production)
dxdt_IF_POOL = KIF_REGEN * (PARIETAL - IF_POOL) - KIF_DECAY * AUTO_AB * IF_POOL;

// 4. CBA_PORT — Portal blood
double IM_to_plasma  = KA_IM * DEPOT;
double oral_to_portal = IF_mediated + passive_abs;
double entero_return  = ENTERO_F * LIVER;
dxdt_CBA_PORT = oral_to_portal + entero_return - 0.5 * CBA_PORT;

// 5. PLASMA (amount = Cp × V_PLASMA)
double influx = IM_to_plasma + 0.5 * CBA_PORT;
double liver_uptake   = K12_LIVER * PLASMA;
double liver_release  = K21_LIVER * LIVER;
double bm_uptake      = K12_BM * PLASMA;
double bm_release     = K21_BM * BONE_MRW;
double nerve_uptake   = K12_NERVE * PLASMA;
double nerve_release  = K21_NERVE * NERVE;
double renal_cl       = (Cp > 300.0) ? KREN * (Cp - 300.0) * V_PLASMA : 0.0;  // excrete excess
dxdt_PLASMA = influx + liver_release + bm_release + nerve_release
              - liver_uptake - bm_uptake - nerve_uptake - renal_cl;

// 6. HOLOTC — active B12 (pmol/L), tracks ~20% of plasma B12
// HoloTC ~ 0.2 × Cp (pg/mL) × 0.738 (pmol conversion factor)
double holoTC_target = 0.20 * Cp * 0.738;
dxdt_HOLOTC = 0.10 * (holoTC_target - HOLOTC);

// 7. LIVER
dxdt_LIVER = liver_uptake - liver_release - entero_return;

// 8. BONE_MRW
dxdt_BONE_MRW = bm_uptake - bm_release;

// 9. NERVE
dxdt_NERVE = nerve_uptake - nerve_release;

// 10. AUTO_AB — autoimmune antibody dynamics
// Slow natural decay; treatment does not significantly reduce Ab
double ab_prod = KAB_PROD * (1.0 - AUTO_AB);
double ab_decay = KAB_DECAY * AUTO_AB;
dxdt_AUTO_AB = ab_prod - ab_decay;

// 11. PARIETAL — ongoing destruction by autoimmunity
dxdt_PARIETAL = -KIF_DECAY * AUTO_AB * PARIETAL + KIF_REGEN * (0.05 - PARIETAL);
// minimal regeneration; floor at ~5% in severe PA

// 12. HGB — driven by bone marrow B12 availability
double HGB_SS = HGB0 + (HGB_MAX - HGB0) * BONE_MRW / (BONE_MRW + EC50_HGB / 1000.0);
// Ramp: use plasma B12 for fast response signal
double HGB_target = HGB0 + (HGB_MAX - HGB0) * (Cp / (Cp + EC50_HGB)) * EMAX_HGB;
dxdt_HGB = KOUT_HGB * (HGB_target - HGB);

// 13. MCV — normalization tracks B12 correction
double MCV_target = MCV_NL + (MCV0 - MCV_NL) * (1.0 - Cp / (Cp + EC50_MCV));
dxdt_MCV = K_MCV * (MCV_target - MCV);

// 14. RETIC — reticulocyte crisis response (peaks 4-10 days after treatment)
double RETIC_target = RETIC_BASE + (RETIC_MAX - RETIC_BASE) * (Cp / (Cp + RETIC_EC50));
dxdt_RETIC = K_RETIC * (RETIC_target - RETIC);

// 15. NEURO — progressive damage without B12, partial recovery with
double neuro_protect = Cp / (Cp + EC50_NEURO);
double neuro_progress = KNEURO_PROG * (NEURO_MAX - NEURO) * (1.0 - neuro_protect);
double neuro_recovery = KNEURO_RECOV * (NEURO - NEURO_MIN) * neuro_protect;
dxdt_NEURO = neuro_progress - neuro_recovery;

$TABLE
double Cp_pgmL   = PLASMA / V_PLASMA;          // pg/mL
double MMA       = MMA_NORM + (MMA_BASE - MMA_NORM) * exp(-K_MMA * (Cp_pgmL / 100.0));
double HCY       = HCY_NORM + (HCY_BASE - HCY_NORM) * exp(-K_HCY * (Cp_pgmL / 200.0));
double IF_func   = IF_POOL * 100.0;             // % parietal function
double B12_STORE = LIVER;                       // µg
capture Cp_pgmL MMA HCY IF_func B12_STORE

$CAPTURE
Cp_pgmL HGB MCV RETIC NEURO HOLOTC MMA HCY IF_func PARIETAL AUTO_AB B12_STORE LIVER BONE_MRW NERVE
'

pa_model <- mrgsolve::mcode("PerniciousAnemia_QSP", pa_model_code)

# ============================================================
# TREATMENT EVENT SCHEDULES
# ============================================================

make_events <- function(scenario) {
  switch(scenario,
    # S1: Untreated PA
    "S1_untreated" = ev(),

    # S2: IM loading 1000 µg/day × 7 days, then monthly
    "S2_IM_standard" = ev(
      data.frame(
        ID   = 1,
        time = c(0, 24, 48, 72, 96, 120, 144,          # loading 7 days
                 seq(720, 720*24, by=720)),              # monthly (720h = 30 days)
        amt  = c(rep(1000, 7), rep(1000, 24)),
        cmt  = 1,  # DEPOT (IM)
        evid = 1
      )
    ),

    # S3: IM maintenance only (monthly, no loading)
    "S3_IM_maintenance" = ev(
      data.frame(
        ID   = 1,
        time = seq(0, 720*24, by=720),
        amt  = 1000,
        cmt  = 1,
        evid = 1
      )
    ),

    # S4: High-dose oral 1000 µg/day (1% passive absorption)
    "S4_oral_HD" = ev(
      data.frame(
        ID   = 1,
        time = seq(0, 24*730, by=24),   # daily for 2 years
        amt  = 1000,
        cmt  = 2,  # ORAL_GI
        evid = 1
      )
    ),

    # S5: Aggressive — IM 2000 µg/day × 14 days, then biweekly
    "S5_aggressive" = ev(
      data.frame(
        ID   = 1,
        time = c(seq(0, 13*24, by=24),             # loading 14 days
                 seq(336, 336+360*24*2, by=336)),   # biweekly (336h = 14 days)
        amt  = c(rep(2000, 14), rep(1000, ceiling(360*24*2/336))),
        cmt  = 1,
        evid = 1
      )
    )
  )
}

# ============================================================
# SIMULATION PARAMETERS
# ============================================================
TSIM <- 365 * 24  # 1 year in hours
DT   <- 12        # output every 12 hours

run_scenario <- function(scen_name, scenario_id) {
  evt <- make_events(scen_name)
  if (inherits(evt, "ev") && nrow(as.data.frame(evt)) == 0) {
    # No events — just simulate disease course
    out <- pa_model %>%
      mrgsim(end=TSIM, delta=DT) %>%
      as.data.frame()
  } else {
    df_ev <- as.data.frame(evt)
    df_ev <- df_ev[df_ev$time <= TSIM, ]
    out <- pa_model %>%
      ev(df_ev) %>%
      mrgsim(end=TSIM, delta=DT) %>%
      as.data.frame()
  }
  out$Scenario <- scenario_id
  out$time_days <- out$time / 24
  out
}

# ============================================================
# RUN ALL SCENARIOS
# ============================================================
cat("Running Pernicious Anemia QSP simulations...\n")

scenarios <- list(
  list(name="S1_untreated",    id="S1: Untreated PA"),
  list(name="S2_IM_standard",  id="S2: IM Loading + Monthly Maintenance"),
  list(name="S3_IM_maintenance",id="S3: IM Monthly Maintenance Only"),
  list(name="S4_oral_HD",      id="S4: High-dose Oral (1000 µg/day)"),
  list(name="S5_aggressive",   id="S5: Aggressive IM (2000 µg × 14d → Biweekly)")
)

results_list <- lapply(scenarios, function(s) {
  cat("  Scenario:", s$id, "\n")
  tryCatch(
    run_scenario(s$name, s$id),
    error = function(e) {
      cat("    ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
})
results_list <- Filter(Negate(is.null), results_list)
results <- bind_rows(results_list)

cat("Simulation complete. Rows:", nrow(results), "\n\n")

# ============================================================
# SUMMARY TABLE AT 6 MONTHS AND 12 MONTHS
# ============================================================
summary_table <- results %>%
  filter(time_days %in% c(0, 30, 90, 180, 365)) %>%
  select(Scenario, time_days, Cp_pgmL, HGB, MCV, RETIC, NEURO,
         HOLOTC, MMA, HCY, B12_STORE) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("=== Summary at Key Timepoints ===\n")
print(summary_table)

# ============================================================
# VISUALIZATION
# ============================================================
colors_scen <- c(
  "S1: Untreated PA"                              = "#E74C3C",
  "S2: IM Loading + Monthly Maintenance"          = "#2ECC71",
  "S3: IM Monthly Maintenance Only"               = "#3498DB",
  "S4: High-dose Oral (1000 µg/day)"             = "#F39C12",
  "S5: Aggressive IM (2000 µg × 14d → Biweekly)" = "#9B59B6"
)

theme_qsp <- theme_bw(base_size=11) +
  theme(
    legend.position="bottom",
    legend.title=element_blank(),
    legend.text=element_text(size=8),
    plot.title=element_text(face="bold", size=12),
    strip.background=element_rect(fill="#F0F3F4")
  )

# --- Plot 1: Plasma B12 ---
p1 <- ggplot(results, aes(x=time_days, y=Cp_pgmL, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=200, linetype="dashed", color="gray40", linewidth=0.7) +
  geom_hline(yintercept=900, linetype="dotted", color="gray40", linewidth=0.7) +
  annotate("text", x=350, y=215, label="Deficiency threshold (200 pg/mL)",
           color="gray40", size=3, hjust=1) +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(ylim=c(0, 2000)) +
  labs(title="Plasma Cobalamin (B12)", x="Time (days)", y="Plasma B12 (pg/mL)") +
  theme_qsp

# --- Plot 2: Hemoglobin ---
p2 <- ggplot(results, aes(x=time_days, y=HGB, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=12, linetype="dashed", color="gray40", linewidth=0.7) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=12, alpha=0.05, fill="red") +
  annotate("text", x=350, y=12.3, label="Anemia threshold (12 g/dL)", color="gray40", size=3, hjust=1) +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(ylim=c(5, 16)) +
  labs(title="Hemoglobin", x="Time (days)", y="Hb (g/dL)") +
  theme_qsp

# --- Plot 3: MCV ---
p3 <- ggplot(results, aes(x=time_days, y=MCV, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=100, linetype="dashed", color="gray40", linewidth=0.7) +
  annotate("text", x=350, y=101.5, label="Macrocytosis threshold (100 fL)", color="gray40", size=3, hjust=1) +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(ylim=c(80, 125)) +
  labs(title="Mean Corpuscular Volume (MCV)", x="Time (days)", y="MCV (fL)") +
  theme_qsp

# --- Plot 4: Reticulocyte Response ---
p4 <- ggplot(results, aes(x=time_days, y=RETIC, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_vline(xintercept=7, linetype="dotted", color="darkred", linewidth=0.6) +
  annotate("text", x=8, y=350, label="Retic crisis\n(Day 4-10)", color="darkred", size=3) +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(xlim=c(0, 30)) +  # zoom first 30 days
  labs(title="Reticulocyte Crisis (First 30 Days)", x="Time (days)", y="Reticulocytes (×10⁹/L)") +
  theme_qsp

# --- Plot 5: Neurological Score ---
p5 <- ggplot(results, aes(x=time_days, y=NEURO, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=2, linetype="dashed", color="gray40") +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=5, ymax=10, alpha=0.05, fill="red") +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(ylim=c(0, 9)) +
  labs(title="Neurological Disability Score (0-10)", x="Time (days)", y="Neuro Score") +
  theme_qsp

# --- Plot 6: MMA and Homocysteine ---
p6a <- ggplot(results, aes(x=time_days, y=MMA, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=0.4, linetype="dashed", color="gray40") +
  scale_color_manual(values=colors_scen) +
  labs(title="Methylmalonic Acid (MMA)", x="Time (days)", y="MMA (µmol/L)") +
  theme_qsp + theme(legend.position="none")

p6b <- ggplot(results, aes(x=time_days, y=HCY, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=15, linetype="dashed", color="gray40") +
  scale_color_manual(values=colors_scen) +
  labs(title="Homocysteine", x="Time (days)", y="Hcy (µmol/L)") +
  theme_qsp + theme(legend.position="none")

# --- Plot 7: Liver B12 Stores ---
p7 <- ggplot(results, aes(x=time_days, y=B12_STORE, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=2000, linetype="dashed", color="gray40") +
  annotate("text", x=350, y=2100, label="Normal stores (2000 µg)", color="gray40", size=3, hjust=1) +
  scale_color_manual(values=colors_scen) +
  labs(title="Hepatic B12 Stores", x="Time (days)", y="Liver B12 (µg)") +
  theme_qsp

# --- Plot 8: HoloTC ---
p8 <- ggplot(results, aes(x=time_days, y=HOLOTC, color=Scenario)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=35, linetype="dashed", color="gray40") +
  annotate("text", x=350, y=37, label="Deficiency <35 pmol/L", color="gray40", size=3, hjust=1) +
  scale_color_manual(values=colors_scen) +
  coord_cartesian(ylim=c(0, 200)) +
  labs(title="Holo-Transcobalamin II (Active B12)", x="Time (days)", y="HoloTC (pmol/L)") +
  theme_qsp

# --- Combined dashboard ---
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p7) +
  plot_annotation(
    title    = "Pernicious Anemia QSP Model — Treatment Scenario Comparison",
    subtitle = "악성 빈혈 정량적 시스템 약리학 모델 | 5가지 치료 시나리오 비교",
    theme    = theme(
      plot.title    = element_text(face="bold", size=14),
      plot.subtitle = element_text(size=11, color="gray40")
    )
  )

biomarker_plot <- (p6a + p6b) / (p8 + p7) +
  plot_annotation(title="Functional B12 Biomarkers & Stores")

# ============================================================
# SAVE PLOTS
# ============================================================
tryCatch({
  ggsave("pna_treatment_comparison.png",   combined_plot,  width=14, height=18, dpi=150)
  ggsave("pna_biomarker_dynamics.png",     biomarker_plot, width=14, height=10, dpi=150)
  cat("Plots saved.\n")
}, error=function(e) cat("Plot save error:", e$message, "\n"))

# ============================================================
# INDIVIDUAL PARAMETER SENSITIVITY
# ============================================================
cat("\n=== Sensitivity Analysis: Effect of F_PASSIVE on Oral B12 ==\n")
F_vals <- c(0.005, 0.01, 0.015, 0.02)  # 0.5%, 1%, 1.5%, 2% passive absorption

sens_results <- lapply(F_vals, function(f) {
  mod_f <- param(pa_model, F_PASSIVE=f)
  oral_ev <- ev(
    data.frame(ID=1, time=seq(0, 24*365, by=24), amt=1000, cmt=2, evid=1)
  )
  out <- mod_f %>% ev(oral_ev) %>% mrgsim(end=TSIM, delta=DT) %>% as.data.frame()
  out$F_passive_pct <- paste0(f*100, "% passive")
  out$time_days <- out$time / 24
  out
})
sens_df <- bind_rows(sens_results)

p_sens <- ggplot(sens_df, aes(x=time_days, y=HGB, color=F_passive_pct)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=12, linetype="dashed") +
  scale_color_brewer(type="qual", palette="Set2") +
  labs(
    title    = "Sensitivity: Passive Oral B12 Absorption Fraction",
    subtitle = "Effect of F_PASSIVE on Hemoglobin response (oral 1000 µg/day)",
    x = "Time (days)", y = "Hb (g/dL)", color = "Absorption %"
  ) +
  theme_qsp

print(p_sens)

# ============================================================
# KEY SIMULATION OUTPUTS
# ============================================================
cat("\n=== Key Clinical Outcomes at 1 Year ===\n")
year1 <- results %>%
  filter(abs(time_days - 365) < 1) %>%
  select(Scenario, Cp_pgmL, HGB, MCV, NEURO, MMA, HCY, HOLOTC, B12_STORE) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
print(year1)

cat("\n=== Reticulocyte Peak (Crisis Response) ===\n")
retic_peak <- results %>%
  group_by(Scenario) %>%
  filter(time_days <= 30) %>%
  slice_max(RETIC, n=1) %>%
  select(Scenario, time_days, RETIC) %>%
  mutate(RETIC=round(RETIC, 1), time_days=round(time_days, 1))
print(retic_peak)

cat("\n=== Time to Hb >12 g/dL ===\n")
time_to_remission <- results %>%
  filter(HGB >= 12) %>%
  group_by(Scenario) %>%
  slice_min(time_days, n=1) %>%
  select(Scenario, time_days) %>%
  mutate(time_days=round(time_days, 1))
print(time_to_remission)
