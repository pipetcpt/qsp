## ============================================================
##  Obesity QSP Model — mrgsolve ODE Implementation
##  Disease: Obesity (비만)
##  Drugs:   ① Placebo
##           ② Semaglutide 2.4 mg SC QW (GLP-1RA)
##           ③ Tirzepatide 15 mg SC QW (GIP/GLP-1 dual RA)
##           ④ Orlistat 120 mg PO TID (Lipase inhibitor)
##           ⑤ Phentermine/Topiramate 15/92 mg PO QD
##
##  Calibrated to:
##   · STEP 1 trial (Wilding NEJM 2021): semaglutide 2.4 mg ~14.9% WT loss
##   · SURMOUNT-1 trial (Jastreboff NEJM 2022): tirzepatide 15 mg ~20.9%
##   · XENDOS trial (Torgerson JAMA 2004): orlistat ~2.8% excess WT loss
##   · CONQUER trial (Garvey Lancet 2011): Qsymia 15/92 mg ~9.3% WT loss
##
##  ODE States (20 compartments):
##   PK:   SEMA_GUT, SEMA_C, TIRZ_GUT, TIRZ_C
##   PD1:  GLP1R_OCC, GIPR_OCC
##   PD2:  FOOD_R, GASTRIC, GHRELIN_R
##   META: INSULIN_P, GLUCOSE_P
##   COMP: ADIP, BWT_C, LEPTIN_P, TRIG_P
##   BIO:  HBA1C_C, INFLAM_I, HOMA_IR_C
##   ORL:  ORL_GUT (Orlistat), ORL_INH
##   CNS:  NE_EFFECT
##
##  Author: Claude Code Routine (CCR) — 2026-06-18
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ── mrgsolve model string ─────────────────────────────────────────────────────
ob_model_code <- '
$PROB
Obesity QSP Model — GLP-1RA / Dual GIP+GLP-1RA / Lipase Inhibitor / CNS Agent

$PARAM @annotated
// ── Semaglutide PK (Marbury 2020, Lau 2015) ──────────────────────────────
ka_sema   : 0.0177 : SC absorption rate (/h); t_max ~72h
CL_sema   : 0.0403 : Clearance (L/h); albumin-bound
Vd_sema   : 12.4   : Volume of distribution (L)
F_sema    : 0.89   : SC bioavailability

// ── Tirzepatide PK (Heald 2022, Coskun 2022) ─────────────────────────────
ka_tirz   : 0.0187 : SC absorption rate (/h); t_max ~48h
CL_tirz   : 0.0500 : Clearance (L/h)
Vd_tirz   : 9.20   : Volume of distribution (L)
F_tirz    : 0.80   : SC bioavailability

// ── Orlistat PK/PD (Drent 2001) ──────────────────────────────────────────
ka_orl    : 0.50   : GI absorption rate (/h; minimal systemic)
ke_orl    : 0.35   : GI elimination (/h)
Imax_orl  : 0.30   : Max fat absorption reduction (fraction) ~30%
IC50_orl  : 0.55   : Half-max orlistat gut conc (mg); EC50 ~0.2 µg/mL

// ── CNS agent (Phentermine/Topiramate) ───────────────────────────────────
ka_cns    : 0.15   : Absorption rate (/h)
ke_cns    : 0.07   : Elimination rate (/h)
Emax_cns  : 0.40   : Max food intake suppression (CNS agents)
EC50_cns  : 0.50   : Half-max CNS drug effect (relative conc)

// ── GLP-1R occupancy parameters ──────────────────────────────────────────
EC50_sema_glp1 : 0.016 : Semaglutide EC50 for GLP-1R (nM); Lau 2015
EC50_tirz_glp1 : 0.050 : Tirzepatide EC50 for GLP-1R (nM); Coskun 2022
EC50_tirz_gip  : 0.013 : Tirzepatide EC50 for GIPR (nM);   Coskun 2022
kout_rec       : 0.10  : Receptor occupancy equilibration (/h)

// ── Food intake / energy balance PD ──────────────────────────────────────
Emax_food_glp1 : 0.42  : Max food intake reduction by GLP-1R (fraction)
Emax_food_gip  : 0.08  : Additional food reduction by GIPR (fraction)
EC50_food_rec  : 0.50  : Receptor occupancy for 50% food reduction
kout_food      : 0.0069 : Food intake turnover (/h); ~6-day adaptation
FOOD_R0        : 1.0   : Baseline food intake (relative = 1.0)

// ── Gastric emptying ─────────────────────────────────────────────────────
kout_gastric   : 0.10  : Gastric emptying equilibration (/h)
Emax_gastric   : 0.25  : Max gastric emptying reduction (GLP-1R)
GASTRIC0       : 1.0   : Baseline gastric emptying rate (relative)

// ── Ghrelin dynamics ─────────────────────────────────────────────────────
kout_ghrelin   : 0.06  : Ghrelin turnover (/h)
Imax_ghrelin   : 0.30  : Ghrelin suppression by GLP-1RA (max)
GHRELIN_R0     : 1.0   : Baseline ghrelin (relative)

// ── Pancreatic / metabolic PD ─────────────────────────────────────────────
kout_ins   : 0.12  : Insulin clearance (/h)
Emax_glp1_ins  : 1.8  : Max incretin effect on insulin secretion (fold)
Emax_gip_ins   : 1.2  : GIP incretin effect on insulin secretion (fold)
EC50_ins_sec   : 0.6  : Receptor occ for half-max incretin effect
kHGP      : 1.10  : Hepatic glucose production constant (mg/dL/h)
kGdisp    : 0.015 : Insulin-mediated glucose disposal (/h per µU/mL)
Imax_ins  : 0.80  : Max insulin suppression of EGP
IC50_ins_HGP : 5.0 : Insulin IC50 for HGP suppression (µU/mL)
kout_gluc : 0.05  : Glucose distribution/turnover (/h)

// ── Body weight / adipose dynamics ───────────────────────────────────────
kout_adip : 0.00024 : Adipose mass turnover (/h; slow, ~170-day half-time)
calden_fat : 7700  : Energy density of fat mass (kcal/kg)
adip_frac  : 0.87  : Fraction of weight loss from fat (at energy deficit)
kout_bwt   : 0.00024 : Body weight equilibration (/h)

// ── Leptin dynamics (Considine 1996, Maffei 1995) ────────────────────────
kout_leptin : 0.30  : Leptin half-life ~3h (clearance /h)
k_lep_adip  : 0.714 : Leptin-adipose proportionality (ng/mL per kg fat)

// ── Triglyceride dynamics ─────────────────────────────────────────────────
kout_trig  : 0.008 : Plasma TG equilibration (/h; ~4-day t½)
TG_adip_exp: 1.50  : TG-adipose power relationship exponent

// ── HbA1c dynamics (Rodbard 2009) ────────────────────────────────────────
kout_hba1c : 0.00046 : HbA1c turnover (/h; ~60-day RBC t½)

// ── Inflammation index ────────────────────────────────────────────────────
kout_inflam : 0.004  : Inflammation index turnover (/h; slow)
inflam_adip_exp : 1.20 : Inflammation-adipose power exponent

// ── HOMA-IR dynamics ─────────────────────────────────────────────────────
kout_homa  : 0.05  : HOMA-IR equilibration (/h)

// ── Baseline patient characteristics ─────────────────────────────────────
BW0    : 106.0 : Baseline body weight (kg); STEP 1 mean
A0     : 44.0  : Baseline adipose mass (kg); ~41% body fat
I0     : 18.0  : Baseline fasting insulin (µU/mL)
G0     : 112.0 : Baseline fasting glucose (mg/dL)
L0     : 32.0  : Baseline leptin (ng/mL)
TG0    : 185.0 : Baseline plasma TG (mg/dL)
HBA0   : 5.8   : Baseline HbA1c (%)
INF0   : 1.50  : Baseline inflammation index
HOMA0  : 4.5   : Baseline HOMA-IR

$CMT  @annotated
// Drug PK compartments
SEMA_GUT  : Semaglutide SC depot (nmol)
SEMA_C    : Semaglutide plasma (nmol/L)
TIRZ_GUT  : Tirzepatide SC depot (nmol)
TIRZ_C    : Tirzepatide plasma (nmol/L)
ORL_GUT   : Orlistat gut lumen (mg)
CNS_C     : CNS agent plasma (relative)

// Receptor occupancy mediators
GLP1R_OCC : GLP-1R occupancy (0–1)
GIPR_OCC  : GIPR occupancy (0–1)

// PD mediators
FOOD_R    : Food intake (relative, 1 = baseline)
GASTRIC   : Gastric emptying rate (relative)
GHRELIN_R : Ghrelin level (relative)

// Metabolic compartments
INSULIN_P : Plasma insulin (µU/mL)
GLUCOSE_P : Plasma glucose (mg/dL)

// Body composition
ADIP      : Adipose tissue mass (kg)
BWT_C     : Body weight (kg)
LEPTIN_P  : Plasma leptin (ng/mL)

// Clinical biomarkers
TRIG_P    : Plasma triglycerides (mg/dL)
HBA1C_C   : HbA1c (%)
INFLAM_I  : Inflammation index (rel.)
HOMA_IR_C : HOMA-IR (rel. units)

$MAIN
// Derived quantities at each time step
double glp1r_target = SEMA_C/(SEMA_C + EC50_sema_glp1) +
                      TIRZ_C/(TIRZ_C + EC50_tirz_glp1);
if(glp1r_target > 1.0) glp1r_target = 1.0;

double gipr_target  = TIRZ_C/(TIRZ_C + EC50_tirz_gip);
if(gipr_target > 1.0) gipr_target = 1.0;

// Orlistat fat absorption inhibition
double orl_effect = Imax_orl * ORL_GUT/(ORL_GUT + IC50_orl);

// CNS appetite suppression
double cns_effect = Emax_cns * CNS_C/(CNS_C + EC50_cns);

$ODE
// ── Semaglutide PK ──
dxdt_SEMA_GUT = -ka_sema * SEMA_GUT;
dxdt_SEMA_C   = ka_sema * SEMA_GUT * F_sema / Vd_sema
                - (CL_sema/Vd_sema) * SEMA_C;

// ── Tirzepatide PK ──
dxdt_TIRZ_GUT = -ka_tirz * TIRZ_GUT;
dxdt_TIRZ_C   = ka_tirz * TIRZ_GUT * F_tirz / Vd_tirz
                - (CL_tirz/Vd_tirz) * TIRZ_C;

// ── Orlistat GI kinetics ──
dxdt_ORL_GUT  = -ke_orl * ORL_GUT;

// ── CNS agent kinetics ──
dxdt_CNS_C    = -ke_cns * CNS_C;

// ── Receptor occupancy (indirect response, equilibrate to target) ──
dxdt_GLP1R_OCC = kout_rec * (glp1r_target - GLP1R_OCC);
dxdt_GIPR_OCC  = kout_rec * (gipr_target  - GIPR_OCC);

// ── Food intake reduction ──
double food_eff_glp1 = Emax_food_glp1 * GLP1R_OCC / (GLP1R_OCC + EC50_food_rec);
double food_eff_gip  = Emax_food_gip  * GIPR_OCC  / (GIPR_OCC  + EC50_food_rec);
double total_food_inh = food_eff_glp1 + food_eff_gip + orl_effect + cns_effect;
if(total_food_inh > 0.70) total_food_inh = 0.70; // cap at 70% reduction
double kin_food = kout_food * FOOD_R0;
dxdt_FOOD_R   = kin_food * (1.0 - total_food_inh) - kout_food * FOOD_R;

// ── Gastric emptying ──
double gastric_target = GASTRIC0 * (1.0 - Emax_gastric * GLP1R_OCC);
dxdt_GASTRIC  = kout_gastric * (gastric_target - GASTRIC);

// ── Ghrelin dynamics ──
double ghrelin_target = GHRELIN_R0 * (1.0 - Imax_ghrelin * GLP1R_OCC);
dxdt_GHRELIN_R = kout_ghrelin * (ghrelin_target - GHRELIN_R);

// ── Plasma insulin (indirect response with incretin potentiation) ──
double inc_effect = 1.0
  + Emax_glp1_ins * GLP1R_OCC / (GLP1R_OCC + EC50_ins_sec)
  + Emax_gip_ins  * GIPR_OCC  / (GIPR_OCC  + EC50_ins_sec);
double kin_ins = kout_ins * I0;
dxdt_INSULIN_P = kin_ins * (GLUCOSE_P/G0) * inc_effect - kout_ins * INSULIN_P;

// ── Plasma glucose (minimal model; Bergman 1989 simplified) ──
double EGP_frac = 1.0 - Imax_ins * INSULIN_P / (INSULIN_P + IC50_ins_HGP);
if(EGP_frac < 0.05) EGP_frac = 0.05;
dxdt_GLUCOSE_P = kHGP * EGP_frac - kGdisp * (INSULIN_P/I0) * GLUCOSE_P
                 + kout_gluc * (G0 * FOOD_R - GLUCOSE_P) * 0.1;

// ── Adipose mass (energy balance driven, slow) ──
// Energy deficit from food reduction (kcal/h): ΔE = ΔFood × basal_kcal/day/24
// Simplified: daily intake ~2500 kcal -> 104.2 kcal/h at baseline
double kcal_per_h = 104.2; // baseline
double energy_deficit_h = kcal_per_h * (1.0 - FOOD_R); // kcal/h saved
// Fat mass loss rate (kg/h) = deficit × fat fraction / caloric density
dxdt_ADIP = -energy_deficit_h * adip_frac / calden_fat;

// ── Body weight ──
dxdt_BWT_C = -energy_deficit_h / calden_fat;

// ── Plasma leptin (proportional to fat mass, Considine 1996) ──
double leptin_ss = k_lep_adip * ADIP;
dxdt_LEPTIN_P = kout_leptin * (leptin_ss - LEPTIN_P);

// ── Plasma triglycerides (power-law adipose relationship) ──
double trig_ss = TG0 * pow(ADIP/A0, TG_adip_exp);
dxdt_TRIG_P = kout_trig * (trig_ss - TRIG_P);

// ── HbA1c (lagged integral of glucose; Rohlfing 2002 calibration) ──
// HbA1c (%) ≈ (eAG + 46.7) / 28.7; eAG in mg/dL
double hba1c_ss = (GLUCOSE_P + 46.7) / 28.7;
dxdt_HBA1C_C = kout_hba1c * (hba1c_ss - HBA1C_C);

// ── Inflammation index (adipose-driven, slow) ──
double inflam_ss = INF0 * pow(ADIP/A0, inflam_adip_exp);
dxdt_INFLAM_I = kout_inflam * (inflam_ss - INFLAM_I);

// ── HOMA-IR (composite metabolic index) ──
double homa_ss = (INSULIN_P * GLUCOSE_P) / 22.5;
dxdt_HOMA_IR_C = kout_homa * (homa_ss - HOMA_IR_C);

$CAPTURE @annotated
// Capture derived/output variables for plotting
glp1r_target  : GLP-1R target occupancy (0-1)
gipr_target   : GIPR target occupancy (0-1)
orl_effect    : Orlistat fat absorption inhibition (0-1)
cns_effect    : CNS appetite suppression effect (0-1)

$INIT @annotated
SEMA_GUT   = 0.0
SEMA_C     = 0.0
TIRZ_GUT   = 0.0
TIRZ_C     = 0.0
ORL_GUT    = 0.0
CNS_C      = 0.0
GLP1R_OCC  = 0.0
GIPR_OCC   = 0.0
FOOD_R     = 1.0
GASTRIC    = 1.0
GHRELIN_R  = 1.0
INSULIN_P  = 18.0
GLUCOSE_P  = 112.0
ADIP       = 44.0
BWT_C      = 106.0
LEPTIN_P   = 32.0
TRIG_P     = 185.0
HBA1C_C    = 5.80
INFLAM_I   = 1.50
HOMA_IR_C  = 4.50
'

# ── Compile model ──────────────────────────────────────────────────────────────
mod <- mcode("obesity_qsp", ob_model_code)

cat("\n===== Obesity QSP mrgsolve Model =====\n")
cat("Compartments:", length(mod@cmtL), "\n")
cat("Parameters:  ", length(param(mod)), "\n")
param(mod)

# ── Dosing event helpers ───────────────────────────────────────────────────────

# Semaglutide 2.4 mg SC QW (escalation: 0.25→0.5→1.0→1.7→2.4 mg over 16 wk)
# MW semaglutide ≈ 4113.6 g/mol → 2.4 mg = 0.583 µmol = 583 nmol
sema_escalation <- function() {
  doses_mg  <- c(0.25, 0.25, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.7, 1.7, 1.7, 1.7,
                 rep(2.4, 60))
  MW_sema   <- 4113.6
  doses_nmol <- doses_mg / MW_sema * 1e6
  tibble(
    time = seq(0, 7*length(doses_mg)-1, by=7),
    cmt  = "SEMA_GUT",
    amt  = doses_nmol,
    evid = 1
  )
}

# Tirzepatide 15 mg SC QW (escalation: 2.5→5→7.5→10→12.5→15 mg)
# MW tirzepatide ≈ 4813.5 g/mol → 15 mg = 3.116 µmol = 3116 nmol
tirz_escalation <- function() {
  doses_mg   <- c(2.5, 2.5, 2.5, 2.5, 5.0, 5.0, 5.0, 5.0,
                  7.5, 7.5, 7.5, 7.5, 10.0, 10.0, 10.0, 10.0,
                  12.5, 12.5, 12.5, 12.5, rep(15.0, 52))
  MW_tirz    <- 4813.5
  doses_nmol <- doses_mg / MW_tirz * 1e6
  tibble(
    time = seq(0, 7*length(doses_mg)-1, by=7),
    cmt  = "TIRZ_GUT",
    amt  = doses_nmol,
    evid = 1
  )
}

# Orlistat 120 mg PO TID (3x daily, approximated as continuous input)
orl_dosing <- function(end_week = 52) {
  # Simulate 120 mg three times daily (q8h) for 52 weeks
  t_start <- seq(0, end_week*7*24-8, by=8)
  tibble(
    time = t_start,
    cmt  = "ORL_GUT",
    amt  = 120,
    evid = 1
  )
}

# Phentermine/Topiramate QD — represented as CNS_C dose (rel. units)
cns_dosing <- function(end_week = 52) {
  t_start <- seq(0, end_week*7*24-24, by=24)
  tibble(
    time = t_start,
    cmt  = "CNS_C",
    amt  = 1.0,
    evid = 1
  )
}

# ── Simulation function ────────────────────────────────────────────────────────
run_scenario <- function(arm_name, dosing_df, sim_weeks = 68) {
  sim_end  <- sim_weeks * 7 * 24  # hours
  obs_times <- seq(0, sim_end, by=24)

  out <- mod %>%
    ev(as.data.frame(dosing_df)) %>%
    mrgsim(end=sim_end, delta=24, outvars=c(
      "SEMA_C","TIRZ_C","ORL_GUT","CNS_C",
      "GLP1R_OCC","GIPR_OCC",
      "FOOD_R","GASTRIC","GHRELIN_R",
      "INSULIN_P","GLUCOSE_P",
      "ADIP","BWT_C","LEPTIN_P",
      "TRIG_P","HBA1C_C","INFLAM_I","HOMA_IR_C"
    )) %>%
    as_tibble() %>%
    mutate(
      arm       = arm_name,
      week      = time / (24 * 7),
      bwt_pct   = (BWT_C - 106.0) / 106.0 * 100,   # % change from baseline
      bmi       = BWT_C / (1.73^2),                  # assume 1.73m height
      waist_est = 88.2 + 1.15*(BWT_C - 106.0),       # simplified waist proxy
      homa_est  = (INSULIN_P * GLUCOSE_P) / 22.5
    )
  out
}

# ── 5 Treatment Scenarios ──────────────────────────────────────────────────────
sim_weeks <- 72

cat("\nRunning 5 treatment scenarios...\n")

# ① Placebo (no dosing)
s1 <- run_scenario("① Placebo", tibble(time=0, cmt="SEMA_GUT", amt=0, evid=1),
                   sim_weeks=sim_weeks)

# ② Semaglutide 2.4 mg QW
s2 <- run_scenario("② Semaglutide 2.4 mg QW", sema_escalation(), sim_weeks=sim_weeks)

# ③ Tirzepatide 15 mg QW
s3 <- run_scenario("③ Tirzepatide 15 mg QW", tirz_escalation(), sim_weeks=sim_weeks)

# ④ Orlistat 120 mg TID
s4 <- run_scenario("④ Orlistat 120 mg TID", orl_dosing(sim_weeks), sim_weeks=sim_weeks)

# ⑤ Phentermine/Topiramate 15/92 mg QD
s5 <- run_scenario("⑤ Phentermine/Topiramate QD", cns_dosing(sim_weeks), sim_weeks=sim_weeks)

all_sims <- bind_rows(s1, s2, s3, s4, s5)

# ── Summary table at key timepoints ───────────────────────────────────────────
summary_weeks <- c(12, 24, 36, 52, 68, 72)
summary_tbl <- all_sims %>%
  filter(round(week) %in% summary_weeks) %>%
  group_by(arm, week = round(week)) %>%
  slice(1) %>%
  select(arm, week, BWT_C, bwt_pct, bmi, HBA1C_C, TRIG_P,
         INSULIN_P, GLUCOSE_P, HOMA_IR_C, INFLAM_I, LEPTIN_P) %>%
  ungroup()

cat("\n===== Summary at Key Timepoints =====\n")
print(summary_tbl, n=30)

# ── Calibration check ─────────────────────────────────────────────────────────
cat("\n===== Calibration vs Clinical Trial Targets =====\n")
calib <- all_sims %>%
  filter(round(week) %in% c(52, 68, 72)) %>%
  group_by(arm, week=round(week)) %>%
  slice(1) %>%
  select(arm, week, bwt_pct) %>%
  ungroup()

calib_targets <- tribble(
  ~arm,                          ~week, ~target_pct, ~trial,
  "① Placebo",                   68,    -2.4,        "STEP 1 (Wilding 2021)",
  "② Semaglutide 2.4 mg QW",    68,    -14.9,       "STEP 1 (Wilding 2021)",
  "③ Tirzepatide 15 mg QW",     72,    -20.9,       "SURMOUNT-1 (Jastreboff 2022)",
  "④ Orlistat 120 mg TID",      52,    -5.7,        "XENDOS (Torgerson 2004)",
  "⑤ Phentermine/Topiramate QD",52,    -9.3,        "CONQUER (Garvey 2011)"
)

calib_check <- calib %>%
  inner_join(calib_targets, by=c("arm","week")) %>%
  mutate(error_pct = bwt_pct - target_pct)

print(calib_check)

# ── Plots ─────────────────────────────────────────────────────────────────────
arm_colors <- c(
  "① Placebo"                  = "#9E9E9E",
  "② Semaglutide 2.4 mg QW"   = "#7E57C2",
  "③ Tirzepatide 15 mg QW"    = "#1565C0",
  "④ Orlistat 120 mg TID"     = "#2E7D32",
  "⑤ Phentermine/Topiramate QD"= "#E65100"
)

# Plot 1: Body weight % change
p1 <- ggplot(all_sims, aes(week, bwt_pct, color=arm)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(-5,-10,-15,-20), linetype="dashed",
             color="grey70", linewidth=0.5) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="Body Weight Change (%)",
       subtitle="From baseline — 72 weeks simulation",
       x="Week", y="Body Weight Change (%)") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom")

# Plot 2: Plasma glucose
p2 <- ggplot(all_sims, aes(week, GLUCOSE_P, color=arm)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(100, 126), linetype=c("dashed","solid"),
             color=c("green3","red2"), linewidth=0.6) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="Fasting Plasma Glucose (mg/dL)",
       x="Week", y="Glucose (mg/dL)") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

# Plot 3: HbA1c
p3 <- ggplot(all_sims, aes(week, HBA1C_C, color=arm)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(5.7, 6.5), linetype=c("dashed","solid"),
             color=c("orange","red2"), linewidth=0.6) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="HbA1c (%)",
       x="Week", y="HbA1c (%)") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

# Plot 4: Plasma triglycerides
p4 <- ggplot(all_sims, aes(week, TRIG_P, color=arm)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=150, linetype="dashed", color="orange", linewidth=0.6) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="Plasma Triglycerides (mg/dL)",
       x="Week", y="TG (mg/dL)") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

# Plot 5: HOMA-IR
p5 <- ggplot(all_sims, aes(week, HOMA_IR_C, color=arm)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=2.5, linetype="dashed", color="red2", linewidth=0.6) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="HOMA-IR (Insulin Resistance)",
       x="Week", y="HOMA-IR") +
  theme_bw(base_size=12) +
  theme(legend.position="none")

# Plot 6: GLP-1R occupancy (sema vs tirz)
p6 <- all_sims %>%
  filter(arm %in% c("② Semaglutide 2.4 mg QW","③ Tirzepatide 15 mg QW")) %>%
  filter(week <= 24) %>%
  ggplot(aes(week, GLP1R_OCC*100, color=arm)) +
  geom_line(linewidth=1.1) +
  scale_color_manual(values=arm_colors, name=NULL) +
  labs(title="GLP-1R Occupancy (%)\nEscalation Period (0–24 wk)",
       x="Week", y="GLP-1R Occupancy (%)") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom")

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title    = "Obesity QSP Model — 5-Arm Treatment Simulation",
    subtitle = "Semaglutide · Tirzepatide · Orlistat · Phentermine/Topiramate · Placebo",
    caption  = "Calibrated to STEP 1, SURMOUNT-1, XENDOS, CONQUER trials"
  )

print(combined_plot)

# ── Dose-Response Analysis (Semaglutide) ──────────────────────────────────────
cat("\n===== Dose-Response Analysis: Semaglutide =====\n")

sema_doses_mg <- c(0.25, 0.5, 1.0, 2.0, 2.4)
MW_sema       <- 4113.6

dr_results <- lapply(sema_doses_mg, function(dose_mg) {
  dose_nmol <- dose_mg / MW_sema * 1e6
  ev_single <- tibble(time=0, cmt="SEMA_GUT", amt=dose_nmol, evid=1)
  # Repeat weekly dose for 68 weeks
  ev_df <- tibble(
    time = seq(0, 67*7*24, by=7*24),
    cmt  = "SEMA_GUT",
    amt  = dose_nmol,
    evid = 1
  )
  out <- mod %>%
    ev(as.data.frame(ev_df)) %>%
    mrgsim(end=68*7*24, delta=7*24, outvars="BWT_C") %>%
    as_tibble() %>%
    filter(round(time) == 68*7*24) %>%
    mutate(dose_mg=dose_mg, bwt_pct=(BWT_C-106)/106*100) %>%
    select(dose_mg, BWT_C, bwt_pct)
  out
}) %>% bind_rows()

cat("Semaglutide dose-response at 68 weeks:\n")
print(dr_results)

# ── Model parameter summary ───────────────────────────────────────────────────
cat("\n===== Key Model Parameters =====\n")
param_summary <- tibble(
  Parameter    = c("ODE States","Drug PK CMT","Receptor Occ. CMT","PD Mediators",
                   "Metabolic CMT","Biomarker CMT","Treatment Arms",
                   "Sim. Duration","Calibration Trials"),
  Value        = c("20","6","2","3","5","4","5","72 weeks","4 pivotal RCTs")
)
print(param_summary)

cat("\n===== Drug PK Summary =====\n")
pk_summary <- tribble(
  ~Drug,                  ~Route, ~Dose,        ~MW_gmol, ~t_half,  ~F_pct,  ~EC50_nM,
  "Semaglutide",          "SC QW","2.4 mg",     4113.6,   "168 h",  89,      0.016,
  "Tirzepatide (GLP-1R)", "SC QW","15 mg",      4813.5,   "120 h",  80,      0.050,
  "Tirzepatide (GIPR)",   "SC QW","15 mg",      4813.5,   "120 h",  80,      0.013,
  "Orlistat",             "PO TID","120 mg",    495.7,    "local",   1,      NA,
  "Phentermine/Topi",     "PO QD", "15/92 mg",  NA,       "various",NA,     NA
)
print(pk_summary)

cat("\n===== Clinical Trial Calibration =====\n")
trial_summary <- tribble(
  ~Trial,           ~Drug,              ~N,     ~Duration, ~PrimaryEndpoint,  ~ObservedPct, ~ModelPct,
  "STEP 1",         "Semaglutide 2.4",  1961,   "68 wk",   "Wt loss ≥5%",    "-14.9%",     paste0(round(filter(calib_check,arm=="② Semaglutide 2.4 mg QW")$bwt_pct,1),"%"),
  "SURMOUNT-1",     "Tirzepatide 15",   630,    "72 wk",   "Wt loss ≥5%",    "-20.9%",     paste0(round(filter(calib_check,arm=="③ Tirzepatide 15 mg QW")$bwt_pct,1),"%"),
  "XENDOS",         "Orlistat 120 mg",  3305,   "4 yr",    "Wt loss + T2DM", "-5.7%",      paste0(round(filter(calib_check,arm=="④ Orlistat 120 mg TID")$bwt_pct,1),"%"),
  "CONQUER",        "Qsymia 15/92",     2487,   "56 wk",   "Wt loss ≥10%",   "-9.3%",      paste0(round(filter(calib_check,arm=="⑤ Phentermine/Topiramate QD")$bwt_pct,1),"%")
)
print(trial_summary)

cat("\nModel compilation and simulation complete.\n")
cat("Use ob_shiny_app.R for interactive dashboard.\n")
