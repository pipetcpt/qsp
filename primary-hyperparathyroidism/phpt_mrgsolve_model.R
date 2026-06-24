## ============================================================
## Primary Hyperparathyroidism (PHPT) QSP Model — mrgsolve
## ============================================================
## Compartments (20 ODEs):
##  1.  A_cin       Cinacalcet gut (absorption)
##  2.  C_cin1      Cinacalcet central compartment
##  3.  C_cin2      Cinacalcet peripheral compartment
##  4.  C_deno      Denosumab plasma (mAb)
##  5.  C_rankl     Free circulating RANKL
##  6.  C_complex   Denosumab-RANKL complex (TMDD)
##  7.  C_alen      Alendronate plasma
##  8.  C_alen_bone Alendronate bone-bound
##  9.  PTH         Intact PTH [plasma] (pmol/L)
## 10.  Ca          Ionized Ca2+ [plasma] (mM)
## 11.  PO4         Phosphate [plasma] (mM)
## 12.  VitD25      25(OH)D [plasma] (nmol/L)
## 13.  VitD125     1,25(OH)2D (Calcitriol) [plasma] (pmol/L)
## 14.  OB          Active osteoblast mass (relative)
## 15.  OC          Active osteoclast mass (relative)
## 16.  RANKL_s     Soluble RANKL (relative, bone origin)
## 17.  BMD_LS      Lumbar spine BMD (g/cm2)
## 18.  BMD_FN      Femoral neck BMD (g/cm2)
## 19.  Ca_urine    Urinary calcium excretion (mg/day)
## 20.  GFR         eGFR (mL/min/1.73m2)
##
## Clinical Scenarios:
##  0. Normal healthy subject (baseline validation)
##  1. Untreated PHPT (mild) — Ca 2.7 mM, PTH 150 pmol/L
##  2. Untreated PHPT (severe) — Ca 3.0 mM, PTH 300 pmol/L
##  3. Cinacalcet 60 mg/day (calcimimetic, medical management)
##  4. Denosumab 60 mg SC q6mo (skeletal protection)
##  5. Parathyroidectomy — PTH normalized (simulated resection)
##  6. Cinacalcet + Denosumab (combination: PHPT + skeletal)
##  7. CKD-PHPT: secondary HPT with cinacalcet + calcitriol
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- mrgsolve model code ----

code <- '
$PARAM
// ---- Disease severity parameters ----
PHPT_severity = 1    // 0=normal, 1=mild PHPT, 2=severe PHPT
PTx_time      = 9999 // Parathyroidectomy time (days; 9999=no surgery)
PTx_success   = 1    // 1=cure (PTH→normal), 0=persistent PHPT

// ---- PTH kinetics & secretion ----
kel_PTH   = 10.4     // PTH elimination rate (/day; t1/2 ~ 4 min)
PTH_ss0   = 5.5      // Baseline intact PTH (pmol/L; normal ~1-5 pmol/L)
PTH_max   = 300      // Maximum PTH secretion (PHPT severe)
kPTH_Ca   = 2.8      // Ca-mediated PTH suppression gain
Ca_set    = 1.15     // CaSR calcium set-point (mM; ionized)
PTH_base  = 1.0      // Baseline PTH secretion (relative)
nH_CaSR   = 3.5      // Hill coefficient CaSR-Ca response
CaSR_Cin_shift = 0   // Cinacalcet shift in Ca set-point (mM; PD effect)

// ---- Calcium homeostasis ----
Ca0       = 1.20     // Baseline ionized Ca2+ (mM)
kel_Ca    = 0.8      // Ca buffering rate constant (/day)
PTH_Ca_bone = 0.025  // PTH effect on bone Ca release (mM/pmol/L/day)
PTH_Ca_TmCa = 0.003  // PTH effect on renal Ca reabsorption (mM/pmol/L/day)
VitD125_Ca_int = 0.0015 // 1,25D effect on intestinal Ca absorption
Ca_loss_base = 0.15  // Baseline urinary Ca loss rate (mM/day)
Ca_GFR_factor = 0.003 // GFR impact on Ca filtration

// ---- Phosphate homeostasis ----
PO4_0     = 1.0      // Baseline phosphate (mM)
kel_PO4   = 1.2      // PO4 equilibration rate (/day)
PTH_PO4   = -0.08    // PTH phosphaturic effect (mM/pmol/L/day)

// ---- Vitamin D metabolism ----
VitD25_0  = 75       // Baseline 25(OH)D (nmol/L)
VitD125_0 = 90       // Baseline 1,25(OH)2D (pmol/L)
k_CYP27B1 = 0.05     // PTH stimulation of 1α-hydroxylase
k_CYP24A1 = 0.008    // VitD125 autoinduction of 24-hydroxylase
k_FGF23_VD = 0.003   // FGF23 suppression of CYP27B1
kel_VD25  = 0.007    // 25-OH-D clearance (/day)
kel_VD125 = 0.15     // 1,25-OH2-D clearance (/day)
GFR_VD_factor = 0.015 // GFR impact on CYP27B1 activity

// ---- Bone remodeling (OB/OC coupling) ----
OB0       = 1.0      // Baseline OB mass (normalized)
OC0       = 1.0      // Baseline OC mass (normalized)
kOB_form  = 0.05     // OB formation rate (/day)
kOC_form  = 0.08     // OC formation rate (/day)
kOB_death = 0.05     // OB death rate (/day)
kOC_death = 0.08     // OC death rate (/day)
PTH_OC    = 0.008    // PTH → RANKL↑ → OC increase gain
PTH_OB    = 0.003    // PTH → OB anabolic gain (intermittent; chronic net catabolic)
RANKL_OPG = 1.5      // RANKL/OPG ratio (PHPT elevated)
RANKL_s0  = 1.0      // Baseline soluble RANKL (normalized)
kRANKL    = 0.05     // RANKL production rate
kRANKL_deg= 0.05     // RANKL degradation rate
BMD_LS0   = 0.960    // Baseline lumbar spine BMD (g/cm2; adult female)
BMD_FN0   = 0.730    // Baseline femoral neck BMD (g/cm2)
kBMD_form = 0.0003   // BMD increase rate (OB contribution, /day)
kBMD_resorb = 0.0004 // BMD decrease rate (OC contribution, /day)

// ---- Renal function & urinary Ca ----
GFR0      = 85       // Baseline eGFR (mL/min/1.73m2)
kGFR_loss = 0.00005  // GFR decline rate (/day) — nephrocalcinosis
Ca_urine0 = 180      // Baseline urinary Ca (mg/day)
PTH_UCa   = 1.8      // PTH-driven Ca filtered load effect (mg/day/pmol/L)
Ca_UCa    = 95       // Ca plasma-driven calciuria (mg/day/mM)

// ---- Cinacalcet PK/PD parameters ----
Cin_dose  = 0        // Cinacalcet dose (mg/day; 0=not used)
ka_cin    = 6.0      // Absorption rate constant (/day; Tmax ~3-4h)
F_cin     = 0.22     // Bioavailability (with food ~0.40; fasted ~0.22)
Vc_cin    = 55       // Central Vd (L)
Vp_cin    = 180      // Peripheral Vd (L)
CL_cin    = 125      // Clearance (L/day)
Q_cin     = 167      // Intercompartment CL (L/day)
EC50_cin  = 3.0      // Cinacalcet EC50 for CaSR shift (ng/mL)
Emax_cin  = 0.35     // Maximum CaSR set-point shift (mM)

// ---- Denosumab PK/PD parameters (TMDD) ----
Deno_dose = 0        // Denosumab dose (mg per administration; 0=not used)
Deno_interval = 180  // Dosing interval (days)
ka_deno   = 0.069    // SC absorption rate (/day; Tmax ~10 days)
F_deno    = 0.62     // Bioavailability (SC)
Vd_deno   = 3.1      // Volume of distribution (L)
kel_deno  = 0.025    // Linear clearance rate (/day; t1/2 ~28 days)
kon_deno  = 0.001    // RANKL-Deno association rate (1/nM/day)
koff_deno = 0.003    // RANKL-Deno dissociation rate (/day)
kdeg_complex = 0.05  // Complex degradation rate (/day)
EC50_deno = 1000     // Denosumab EC50 for OC suppression (ng/mL)
Emax_deno = 0.85     // Maximum OC suppression (fraction)

// ---- Alendronate PK/PD parameters ----
Alen_dose = 0        // Alendronate dose (mg/week; 0=not used)
Alen_interval = 7    // Dosing interval (days)
F_alen    = 0.006    // Oral bioavailability (0.6%)
ka_alen   = 12.0     // GI absorption rate (/day)
kel_alen  = 2.88     // Plasma elimination (/day; t1/2 ~6h)
k_bone_bind = 1.44   // Bone binding rate (/day)
k_bone_rel  = 0.001  // Bone release rate (/day; long retention)
EC50_alen = 0.05     // Bone concentration EC50 for OC inhibition (µg/g bone)
Emax_alen = 0.75     // Maximum OC apoptosis induction

$CMT
A_cin C_cin1 C_cin2
C_deno C_rankl C_complex
C_alen C_alen_bone
PTH Ca PO4 VitD25 VitD125
OB OC RANKL_s
BMD_LS BMD_FN
Ca_urine GFR

$INIT
A_cin    = 0
C_cin1   = 0
C_cin2   = 0
C_deno   = 0
C_rankl  = 1.0
C_complex = 0
C_alen   = 0
C_alen_bone = 0
PTH      = 5.5
Ca       = 1.20
PO4      = 1.0
VitD25   = 75
VitD125  = 90
OB       = 1.0
OC       = 1.0
RANKL_s  = 1.0
BMD_LS   = 0.960
BMD_FN   = 0.730
Ca_urine = 180
GFR      = 85

$MAIN
// PHPT disease state: elevate PTH secretion rate
double PTH_secretion_rate;
if (SOLVERTIME < PTx_time) {
  if (PHPT_severity == 0) {
    PTH_secretion_rate = PTH_ss0;             // Normal
  } else if (PHPT_severity == 1) {
    PTH_secretion_rate = PTH_ss0 * 5.0;      // Mild PHPT ~25 pmol/L
  } else {
    PTH_secretion_rate = PTH_ss0 * 22.0;     // Severe PHPT ~120 pmol/L
  }
} else {
  // Post-PTx: PTH returns to near-normal (if successful)
  if (PTx_success == 1) {
    PTH_secretion_rate = PTH_ss0 * 1.1;
  } else {
    PTH_secretion_rate = PTH_ss0 * 4.0;     // Persistent PHPT
  }
}

// Cinacalcet PK dosing (continuous oral daily dose)
double Cin_conc = C_cin1 / Vc_cin;           // ng/mL
double CaSR_shift = (Emax_cin * Cin_conc) / (EC50_cin + Cin_conc);
double Ca_set_eff = Ca_set - CaSR_shift;     // Effective Ca set-point

// CaSR-mediated PTH suppression (Hill function)
double CaSR_effect = pow(Ca / Ca_set_eff, nH_CaSR) /
                     (1.0 + pow(Ca / Ca_set_eff, nH_CaSR));

// Denosumab concentration (ng/mL in Vd=3.1L → µg/mL for EC50)
double Deno_conc = C_deno / Vd_deno;         // ng/mL
double OC_inh_deno = (Emax_deno * Deno_conc) / (EC50_deno + Deno_conc);

// Alendronate bone concentration effect on OC
double Alen_bone_conc = C_alen_bone;          // relative units
double OC_inh_alen = (Emax_alen * Alen_bone_conc) / (EC50_alen + Alen_bone_conc);
double OC_total_inh = OC_inh_deno + OC_inh_alen - OC_inh_deno * OC_inh_alen;

// RANKL-mediated OC stimulation
double RANKL_effect = RANKL_s / RANKL_s0;

// VitD125 effect on intestinal Ca absorption (proportion of total)
double VitD_Ca_absorb = VitD125_Ca_int * VitD125;

// GFR-adjusted CYP27B1 activity
double GFR_VD_effect = GFR_VD_factor * (GFR / GFR0);

// PTH-driven Ca mobilization from bone
double PTH_bone_Ca = PTH_Ca_bone * PTH * OC / OB;

// PTH-driven renal Ca reabsorption (TmCa increase)
double PTH_TmCa = PTH_Ca_TmCa * PTH;

// Net Ca change from inputs
double Ca_in = VitD_Ca_absorb + PTH_bone_Ca + PTH_TmCa;
double Ca_out = Ca_loss_base + Ca_GFR_factor * GFR * Ca;

$ODE
// ---- Cinacalcet PK (2-compartment oral) ----
dxdt_A_cin = -ka_cin * A_cin;
dxdt_C_cin1 = ka_cin * F_cin * A_cin
              - (CL_cin / Vc_cin + Q_cin / Vc_cin) * C_cin1
              + (Q_cin / Vp_cin) * C_cin2;
dxdt_C_cin2 = (Q_cin / Vc_cin) * C_cin1 - (Q_cin / Vp_cin) * C_cin2;

// ---- Denosumab PK (1-compartment SC + TMDD) ----
dxdt_C_deno = -kel_deno * C_deno
              - kon_deno * C_deno * C_rankl
              + koff_deno * C_complex;
dxdt_C_rankl = kRANKL - kRANKL_deg * C_rankl
               - kon_deno * C_deno * C_rankl
               + koff_deno * C_complex;
dxdt_C_complex = kon_deno * C_deno * C_rankl
                 - koff_deno * C_complex
                 - kdeg_complex * C_complex;

// ---- Alendronate PK (1-compartment oral + bone binding) ----
dxdt_C_alen = ka_alen * F_alen * A_cin * 0  // placeholder (dose via event)
              - kel_alen * C_alen
              - k_bone_bind * C_alen;
dxdt_C_alen_bone = k_bone_bind * C_alen - k_bone_rel * C_alen_bone;

// ---- PTH dynamics ----
double PTH_release = PTH_secretion_rate * (1.0 - CaSR_effect) * (1.0 - 0.7 * CaSR_shift);
dxdt_PTH = PTH_release - kel_PTH * PTH;

// ---- Calcium dynamics ----
dxdt_Ca = Ca_in - Ca_out;

// ---- Phosphate dynamics ----
double PO4_in  = 0.3 + 0.1 * VitD125 / 90.0;
double PO4_out = kel_PO4 * PO4 + fabs(PTH_PO4) * PTH;  // PTH phosphaturic
dxdt_PO4 = PO4_in - PO4_out;

// ---- Vitamin D metabolism ----
double CYP27B1_activity = k_CYP27B1 * PTH + GFR_VD_effect;
double CYP24A1_activity = k_CYP24A1 * VitD125;
double VD25_prod = 0.52;  // skin + diet synthesis/input (nmol/L/day)
dxdt_VitD25  = VD25_prod - kel_VD25 * VitD25 - CYP27B1_activity * VitD25 * 0.1;
dxdt_VitD125 = CYP27B1_activity * VitD25 * 0.1
               - kel_VD125 * VitD125
               - CYP24A1_activity * VitD125;

// ---- Bone remodeling (OB/OC) ----
double RANKL_net = RANKL_s * (1.0 + PTH_OC * PTH) * (1.0 - OC_total_inh);
dxdt_OB = kOB_form * (1.0 + PTH_OB * PTH * 0.5) - kOB_death * OB;
dxdt_OC = kOC_form * RANKL_net - kOC_death * OC;
dxdt_RANKL_s = PTH_OC * PTH * kRANKL - kRANKL_deg * RANKL_s;

// ---- BMD ----
double BMD_formation = kBMD_form * OB;
double BMD_resorption = kBMD_resorb * OC;
dxdt_BMD_LS = BMD_formation - BMD_resorption;
dxdt_BMD_FN = (BMD_formation - BMD_resorption) * 0.85;  // cortical more affected

// ---- Urinary calcium excretion ----
double UCa_target = Ca_urine0 + PTH_UCa * (PTH - PTH_ss0) + Ca_UCa * (Ca - Ca0);
dxdt_Ca_urine = 0.5 * (UCa_target - Ca_urine);

// ---- eGFR decline ----
double GFR_decline = kGFR_loss * (Ca_urine / Ca_urine0 - 1.0) * GFR;
dxdt_GFR = -GFR_decline;

$TABLE
double Cin_ng_mL = C_cin1 / Vc_cin;
double Deno_ug_mL = C_deno / (Vd_deno * 1000.0);  // convert ng/mL → µg/mL
double Ca_total_mM = Ca * 2.0;   // approximate: ionized ≈ total/2 (simplified)
double PTH_pg_mL  = PTH * 9.425; // pmol/L → pg/mL (MW 9425 Da)
double OC_OB_ratio = OC / OB;
double CTX_index  = OC * 1.5;    // proxy for serum CTX
double P1NP_index = OB * 1.2;    // proxy for serum P1NP
double T_score_LS = (BMD_LS - 0.96) / 0.12;  // T-score (reference mean ± SD)
double T_score_FN = (BMD_FN - 0.73) / 0.10;

$CAPTURE
PTH PTH_pg_mL Ca PO4 VitD125 VitD25
OB OC RANKL_s BMD_LS BMD_FN Ca_urine GFR
Cin_ng_mL C_deno OC_OB_ratio CTX_index P1NP_index
T_score_LS T_score_FN Ca_total_mM
'

## Compile model
mod <- mread_cache("phpt_qsp", tempdir(), code)

## ============================================================
## Dosing event builders
## ============================================================

make_cin_events <- function(dose_mg, days, interval_days = 1) {
  ## Continuous daily oral dosing (split into twice daily for Tmax accuracy)
  ev(amt = dose_mg, cmt = "A_cin", time = seq(0, days, by = interval_days),
     rate = 0, addl = 0)
}

make_deno_events <- function(dose_mg = 60, times_days = c(0, 180, 360)) {
  ev(amt = dose_mg * 1000, cmt = "C_deno", time = times_days,
     rate = 0.069 * dose_mg * 1000)  ## SC absorption modeled via rate
}

make_alen_events <- function(dose_mg = 70, days = 730, interval = 7) {
  ev(amt = dose_mg, cmt = "C_alen", time = seq(0, days, by = interval))
}

## ============================================================
## Scenario Definitions
## ============================================================

scenarios <- list(
  sc0 = list(
    label    = "0. Normal (Healthy)",
    params   = list(PHPT_severity = 0),
    events   = NULL,
    duration = 365
  ),
  sc1 = list(
    label    = "1. Untreated PHPT (Mild)",
    params   = list(PHPT_severity = 1),
    events   = NULL,
    duration = 1825  # 5 years
  ),
  sc2 = list(
    label    = "2. Untreated PHPT (Severe)",
    params   = list(PHPT_severity = 2),
    events   = NULL,
    duration = 1825
  ),
  sc3 = list(
    label    = "3. PHPT + Cinacalcet 60 mg/day",
    params   = list(PHPT_severity = 1, Cin_dose = 60),
    events   = make_cin_events(60, 1825),
    duration = 1825
  ),
  sc4 = list(
    label    = "4. PHPT + Denosumab 60 mg q6mo",
    params   = list(PHPT_severity = 1, Deno_dose = 60),
    events   = make_deno_events(60, c(0, 180, 360, 540, 720, 900, 1080,
                                      1260, 1440, 1620, 1800)),
    duration = 1825
  ),
  sc5 = list(
    label    = "5. PHPT → Parathyroidectomy (day 90)",
    params   = list(PHPT_severity = 1, PTx_time = 90, PTx_success = 1),
    events   = NULL,
    duration = 1825
  ),
  sc6 = list(
    label    = "6. PHPT + Cinacalcet + Denosumab",
    params   = list(PHPT_severity = 1, Cin_dose = 60, Deno_dose = 60),
    events   = bind_rows(
      as_data_frame(make_cin_events(60, 1825)),
      as_data_frame(make_deno_events(60, c(0, 180, 360, 540, 720, 900,
                                           1080, 1260, 1440, 1620, 1800)))
    ),
    duration = 1825
  ),
  sc7 = list(
    label    = "7. CKD-PHPT + Cinacalcet 90 mg",
    params   = list(PHPT_severity = 1, Cin_dose = 90,
                    GFR0 = 30, GFR_VD_factor = 0.004,
                    kGFR_loss = 0.0002),
    events   = make_cin_events(90, 1825),
    duration = 1825
  )
)

## ============================================================
## Run all scenarios
## ============================================================

run_scenario <- function(sc) {
  p <- sc$params
  m <- mod %>% param(p)
  out <- if (is.null(sc$events)) {
    m %>% mrgsim(end = sc$duration, delta = 1, obsonly = TRUE)
  } else {
    m %>% ev(sc$events) %>% mrgsim(end = sc$duration, delta = 1, obsonly = TRUE)
  }
  as_tibble(out) %>% mutate(scenario = sc$label)
}

results <- bind_rows(lapply(scenarios, run_scenario))

## ============================================================
## Summary Table — Week 52 Values
## ============================================================

summary_tbl <- results %>%
  filter(time %in% c(0, 90, 182, 365, 730, 1095, 1825)) %>%
  group_by(scenario, time) %>%
  summarise(
    PTH_pmol_L = round(mean(PTH), 1),
    PTH_pg_mL  = round(mean(PTH_pg_mL), 0),
    Ca_ionized  = round(mean(Ca), 3),
    Ca_total    = round(mean(Ca_total_mM), 2),
    PO4_mM      = round(mean(PO4), 2),
    VitD125     = round(mean(VitD125), 1),
    BMD_LS      = round(mean(BMD_LS), 3),
    T_score_LS  = round(mean(T_score_LS), 2),
    Ca_urine_mg = round(mean(Ca_urine), 0),
    GFR_ml_min  = round(mean(GFR), 1),
    OC_OB_ratio = round(mean(OC_OB_ratio), 2),
    CTX_index   = round(mean(CTX_index), 2),
    .groups = "drop"
  )

print(summary_tbl)

## ============================================================
## Plots
## ============================================================

theme_qsp <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "#E3F2FD"),
        panel.grid.minor = element_blank())

sc_colors <- c(
  "0. Normal (Healthy)"              = "#9E9E9E",
  "1. Untreated PHPT (Mild)"         = "#F44336",
  "2. Untreated PHPT (Severe)"       = "#B71C1C",
  "3. PHPT + Cinacalcet 60 mg/day"   = "#2196F3",
  "4. PHPT + Denosumab 60 mg q6mo"   = "#4CAF50",
  "5. PHPT → Parathyroidectomy (day 90)" = "#FF9800",
  "6. PHPT + Cinacalcet + Denosumab" = "#9C27B0",
  "7. CKD-PHPT + Cinacalcet 90 mg"   = "#795548"
)

p1 <- ggplot(results, aes(time / 365, PTH_pg_mL, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "PTH [Plasma]", x = "Time (years)", y = "iPTH (pg/mL)",
       color = NULL) +
  geom_hline(yintercept = c(15, 65), linetype = "dashed", color = "gray60", size = 0.4) +
  annotate("text", x = 0.1, y = 65, label = "ULN 65 pg/mL", size = 2.5, hjust = 0) +
  theme_qsp

p2 <- ggplot(results, aes(time / 365, Ca_total_mM, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "Total Serum Calcium", x = "Time (years)", y = "Ca (mM)",
       color = NULL) +
  geom_hline(yintercept = c(2.2, 2.6), linetype = "dashed", color = "gray60", size = 0.4) +
  annotate("text", x = 0.1, y = 2.65, label = "ULN 2.6 mM", size = 2.5, hjust = 0) +
  theme_qsp

p3 <- ggplot(results, aes(time / 365, BMD_LS, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "Lumbar Spine BMD", x = "Time (years)", y = "BMD (g/cm²)",
       color = NULL) +
  geom_hline(yintercept = 0.96, linetype = "dashed", color = "gray40", size = 0.4) +
  theme_qsp

p4 <- ggplot(results, aes(time / 365, Ca_urine, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "Urinary Calcium Excretion", x = "Time (years)", y = "CaU (mg/day)",
       color = NULL) +
  geom_hline(yintercept = 300, linetype = "dashed", color = "red3", size = 0.4) +
  annotate("text", x = 0.1, y = 310, label = "Hypercalciuria threshold", size = 2.5, hjust = 0) +
  theme_qsp

p5 <- ggplot(results, aes(time / 365, GFR, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "eGFR Trend", x = "Time (years)", y = "eGFR (mL/min/1.73m²)",
       color = NULL) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "orange3", size = 0.4) +
  theme_qsp

p6 <- ggplot(results, aes(time / 365, OC_OB_ratio, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "OC/OB Ratio (Bone Remodeling Balance)", x = "Time (years)",
       y = "OC/OB ratio", color = NULL) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 0.4) +
  theme_qsp

p7 <- ggplot(results, aes(time / 365, VitD125, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "1,25(OH)₂D (Calcitriol)", x = "Time (years)",
       y = "Calcitriol (pmol/L)", color = NULL) +
  theme_qsp

p8 <- ggplot(results, aes(time / 365, T_score_LS, color = scenario)) +
  geom_line(size = 0.7, alpha = 0.85) +
  scale_color_manual(values = sc_colors) +
  labs(title = "Lumbar Spine T-score", x = "Time (years)",
       y = "T-score", color = NULL) +
  geom_hline(yintercept = c(-1, -2.5), linetype = "dashed",
             color = c("orange3", "red3"), size = 0.4) +
  annotate("text", x = 3, y = -1.1, label = "Osteopenia", size = 2.5) +
  annotate("text", x = 3, y = -2.6, label = "Osteoporosis", size = 2.5) +
  theme_qsp

## Combine
panel <- (p1 | p2) / (p3 | p4) / (p5 | p6)
panel_biomarkers <- (p7 | p8)

print(panel)
print(panel_biomarkers)

## Cinacalcet PK profile (scenario 3, first 7 days)
cin_pk <- results %>%
  filter(scenario == "3. PHPT + Cinacalcet 60 mg/day", time <= 30)

p_cin_pk <- ggplot(cin_pk, aes(time, Cin_ng_mL)) +
  geom_line(color = "#2196F3", size = 0.8) +
  labs(title = "Cinacalcet PK — First 30 Days (60 mg/day)",
       x = "Time (days)", y = "Cinacalcet (ng/mL)") +
  geom_hline(yintercept = 3.0, linetype = "dashed", color = "red3", size = 0.4) +
  annotate("text", x = 2, y = 3.3, label = "EC50 = 3 ng/mL", size = 2.5) +
  theme_qsp

print(p_cin_pk)

## Week-52 Bar comparison
w52 <- summary_tbl %>%
  filter(time == 365) %>%
  select(scenario, PTH_pg_mL, Ca_total, BMD_LS, Ca_urine_mg, GFR_ml_min)

cat("\n=== Year-1 Summary Table ===\n")
print(as.data.frame(w52))

## ============================================================
## Model documentation: PK parameter table
## ============================================================

pk_params <- data.frame(
  Drug = c("Cinacalcet", "Cinacalcet", "Cinacalcet", "Cinacalcet",
           "Denosumab", "Denosumab", "Denosumab", "Denosumab",
           "Alendronate", "Alendronate", "Alendronate"),
  Parameter = c("F (bioavailability)", "Vc (central Vd, L)",
                "CL (L/day)", "t½ (h)",
                "F (SC)", "Vd (L)", "kel (linear, /day)", "t½ (days)",
                "F (oral)", "k_bone_bind (/day)", "t½ bone (years)"),
  Value = c("22% (fasted), 40% (fed)", "55", "125", "~6-8",
            "62%", "3.1", "0.025", "~28",
            "0.6%", "1.44", "~10"),
  Source = c("Sensipar PI", "Sensipar PI", "Pop-PK", "Clinical data",
             "Prolia PI", "Pop-PK (mAb)", "Prolia PI", "Clinical data",
             "Fosamax PI", "Pop-PK model", "Bone kinetics"),
  stringsAsFactors = FALSE
)

cat("\n=== Drug PK Parameters ===\n")
print(pk_params)
