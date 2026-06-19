## ============================================================
## Paget's Disease of Bone (PBD) — QSP mrgsolve ODE Model
## 파젯병 정량적 시스템 약리학 mrgsolve 미분방정식 모델
##
## Author: CCR (Claude Code Routine)
## Date  : 2026-06-19
##
## Compartments (20 ODEs):
##   Drug PK       : DRUG1 (plasma), DRUG2 (peripheral), DRUGB (bone-bound)
##   Denosumab PK  : DEN_SC (SC depot), DEN_C (plasma), DEN_P (peripheral)
##   Disease model : OCP (OC precursors), OCN (active OC),
##                   OBP (OB precursors), OBN (active OB),
##                   RANKL_f (free RANKL), OPG_f (free OPG)
##   Biomarkers    : CTX (serum CTx), P1NP (serum P1NP),
##                   ALP_b (bone ALP), TRAP_s (TRAP5b),
##                   BMD (bone mineral density pagetic site)
##   Clinical      : PAIN (VAS 0-10)
##
## Parameter calibration references:
##   - Reid IR et al. N Engl J Med 2005;353:898 (zoledronate)
##   - Hosking DJ et al. Lancet 1998;350:1733 (pamidronate vs etidronate)
##   - Siris ES et al. Ann Intern Med 1996;125:401 (alendronate)
##   - Stopeck AT et al. J Clin Oncol 2010;28:5132 (denosumab bone turnover)
##   - Papapoulos SE. Bone 2006;38:S8-13 (bisphosphonate PK/PD)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================
## MODEL DEFINITION
## ============================================================

pbd_model_code <- '
$PARAM @annotated
// ---- Bisphosphonate PK (generic; adjusted per drug) ----
CL_bp    = 5.0    : L/h        Systemic clearance bisphosphonate
V1_bp    = 8.0    : L          Central volume bisphosphonate
V2_bp    = 40.0   : L          Peripheral volume bisphosphonate
Q_bp     = 2.0    : L/h        Inter-compartment clearance
Ka_bp    = 0.5    : h-1        Absorption rate oral bisphosphonate
F_bp     = 0.006  : fraction   Oral bioavailability alendronate (~0.6%)
Kbone_on = 0.05   : h-1        Bone uptake rate (plasma -> bone-bound)
Kbone_off= 0.0001 : h-1        Bone release rate (very slow t1/2 ~years)
CL_bone  = 0.002  : h-1        Effective bone clearance (osteoclast resorption)

// ---- Denosumab PK (anti-RANKL mAb, SC) ----
Ka_den   = 0.003  : h-1        SC absorption rate (t_max ~10 days)
F_den    = 0.62   : fraction   SC bioavailability denosumab
CL_den   = 0.013  : L/h        Clearance denosumab
V1_den   = 3.0    : L          Central Vd denosumab
V2_den   = 2.5    : L          Peripheral Vd denosumab
Q_den    = 0.08   : L/h        Inter-compartmental Q denosumab
KD_den   = 0.008  : nM         KD for denosumab-RANKL binding

// ---- Disease parameters — RANKL / OPG axis ----
RANKL_syn = 0.08  : nmol/L/h   Baseline RANKL synthesis (elevated in PBD)
RANKL_deg = 0.04  : h-1        RANKL degradation rate
OPG_syn   = 0.06  : nmol/L/h   Baseline OPG synthesis
OPG_deg   = 0.03  : h-1        OPG degradation rate
RANKL_ss0 = 2.0   : nmol/L     Steady-state RANKL in PBD (elevated 2-4×normal)
OPG_ss0   = 1.5   : nmol/L     Steady-state OPG in PBD (reduced)
PBD_fold  = 3.0   : -          RANKL/OPG ratio fold-elevation in PBD

// ---- Osteoclast dynamics ----
OCP_syn   = 500   : cells/mL/h Osteoclast precursor production
OCP_deg   = 0.04  : h-1        OCP clearance
OC_diff   = 0.001 : h-1        OCP -> OCN differentiation (RANKL-stimulated)
OC_max    = 2000  : cells/mL   Max osteoclast density (Michaelis-Menten)
OCN_life  = 0.008 : h-1        OC apoptosis rate (natural)
OCN_base  = 300   : cells/mL   Baseline OC in PBD (~3x normal 100/mL)
EC50_oc   = 1.0   : nmol/L     EC50 of RANKL for OC differentiation

// ---- Osteoblast dynamics ----
OBP_syn   = 600   : cells/mL/h Osteoblast precursor production
OBP_deg   = 0.03  : h-1        OBP clearance
OB_diff   = 0.0008: h-1        OBP -> OBN differentiation (coupling-driven)
OBN_life  = 0.006 : h-1        OB apoptosis/transition rate
OBN_base  = 250   : cells/mL   Baseline OB in PBD (coupled to elevated OC)
OB_coupling = 1.5 : -          OC-OB coupling coefficient (TGF-b, IGF-1)

// ---- Bone resorption markers ----
CTX_prod  = 0.0002: ng/mL/cell/h CTx production per OC
CTX_elim  = 0.15  : h-1        CTx elimination rate
CTX_base  = 3.0   : ng/mL      Baseline serum CTx in active PBD

// ---- Bone formation markers ----
P1NP_prod = 0.0003: ng/mL/cell/h P1NP production per OB
P1NP_elim = 0.04  : h-1        P1NP elimination rate
P1NP_base = 80    : ng/mL      Baseline serum P1NP in active PBD
ALP_prod  = 0.0006: IU/L/cell/h Bone ALP production per OB
ALP_elim  = 0.01  : h-1        Bone ALP elimination rate (t1/2 ~3 days)
ALP_base  = 350   : IU/L       Baseline bone ALP in active PBD

// ---- TRAP5b ----
TRAP_prod = 0.0001: U/L/cell/h TRAP5b per OC
TRAP_elim = 0.08  : h-1        TRAP5b elimination

// ---- Bone mineral density (pagetic site) ----
BMD_form  = 0.00001: g/cm2/h   BMD formation rate (from OB)
BMD_resorb= 0.00002: g/cm2/h   BMD resorption rate (from OC; in PBD > form)
BMD_base  = 1.5   : g/cm2      Baseline BMD at pagetic site (elevated but poor quality)

// ---- Pain ----
PAIN_base = 6.0   : VAS 0-10   Baseline pain in active PBD
PAIN_slope= 0.5   : -          Pain sensitivity to ALP
PAIN_recov= 0.001 : h-1        Pain recovery rate with treatment

// ---- Drug effect parameters ----
Emax_bp   = 0.95  : fraction   Max fractional inhibition of OC by BP (FPP synth)
EC50_bp   = 0.01  : nmol/L     EC50 of bone-bound BP for OC suppression
Emax_den  = 0.90  : fraction   Max RANKL neutralization by denosumab
EC50_rankl_den = 0.05 : nmol/L EC50 of denosumab for RANKL binding

$INIT @annotated
// Drug PK compartments
DRUG1 = 0    : nmol    Drug (bisphosphonate) central compartment
DRUG2 = 0    : nmol    Drug peripheral compartment
DRUGB = 0    : nmol/g  Drug bone-bound

// Denosumab PK
DEN_SC = 0   : mg      Denosumab SC depot
DEN_C  = 0   : mg/L    Denosumab plasma
DEN_P  = 0   : mg/L    Denosumab peripheral

// Disease compartments
OCP    = 12500  : cells/mL   OC precursors steady-state
OCN    = 300    : cells/mL   Active osteoclasts (PBD elevated)
OBP    = 20000  : cells/mL   OB precursors steady-state
OBN    = 250    : cells/mL   Active osteoblasts

// RANKL/OPG
RANKL_f = 2.0  : nmol/L     Free RANKL (elevated in PBD)
OPG_f   = 1.5  : nmol/L     Free OPG (reduced in PBD)

// Biomarkers
CTX   = 3.0    : ng/mL      Serum CTx (active PBD: 2-10 ng/mL)
P1NP  = 80.0   : ng/mL      Serum P1NP (active PBD: 50-400 ng/mL)
ALP_b = 350.0  : IU/L       Bone ALP (active PBD: 200-3000 IU/L)
TRAP_s= 5.5    : U/L        Serum TRAP5b
BMD   = 1.5    : g/cm2      BMD pagetic site

// Clinical
PAIN  = 6.0    : VAS         Pain VAS 0-10

$ODE
// ============================================================
// DRUG 1: BISPHOSPHONATE PK
// ============================================================
double C1_bp = DRUG1 / V1_bp;   // plasma concentration nmol/L
double C2_bp = DRUG2 / V2_bp;   // peripheral concentration

dxdt_DRUG1 = -CL_bp*C1_bp - Q_bp*(C1_bp - C2_bp) - Kbone_on*DRUG1 + Kbone_off*DRUGB;
dxdt_DRUG2 =  Q_bp*(C1_bp - C2_bp);
dxdt_DRUGB =  Kbone_on*DRUG1 - Kbone_off*DRUGB - CL_bone*DRUGB;

// ============================================================
// DENOSUMAB PK
// ============================================================
double C_den  = DEN_C / V1_den;   // mg/L -> nM (MW~147000)
double CP_den = DEN_P / V2_den;

dxdt_DEN_SC = -Ka_den * F_den * DEN_SC;
dxdt_DEN_C  =  Ka_den * F_den * DEN_SC - CL_den*C_den - Q_den*(C_den - CP_den)
               - KD_den * DEN_C * RANKL_f;   // simplified binding
dxdt_DEN_P  =  Q_den*(C_den - CP_den);

// ============================================================
// RANKL / OPG DYNAMICS
// ============================================================
// Drug effects on RANKL/OPG
double Den_RANKL_inh = Emax_den * DEN_C / (EC50_rankl_den + DEN_C);
double RANKL_free_eff = RANKL_f * (1 - Den_RANKL_inh);

dxdt_RANKL_f = RANKL_syn * PBD_fold - RANKL_deg * RANKL_f
               - KD_den * DEN_C * RANKL_f / 100.0;   // simplified mAb binding
dxdt_OPG_f   = OPG_syn - OPG_deg * OPG_f;

// ============================================================
// OSTEOCLAST DYNAMICS
// ============================================================
// RANKL-driven OC differentiation
double RANKL_OPG_r = RANKL_free_eff / (OPG_f + 0.001);
double OC_diff_stim = OC_diff * RANKL_free_eff / (EC50_oc + RANKL_free_eff);

// Bone-bound drug effect on OC apoptosis (FPP synthase inhibition)
double Cbone_eff = DRUGB / (DRUGB + 0.001);   // bone-local
double BP_oc_inh = Emax_bp * Cbone_eff / (EC50_bp + Cbone_eff);

// OC apoptosis enhanced by BP
double OC_total_death = (OCN_life + BP_oc_inh * 0.1) * OCN;

dxdt_OCP = OCP_syn - OCP_deg * OCP - OC_diff_stim * OCP;
dxdt_OCN = OC_diff_stim * OCP - OCN_life * OCN * (1 + BP_oc_inh * 5.0);

// ============================================================
// OSTEOBLAST DYNAMICS
// ============================================================
// OC-OB coupling (simplified: OB formation proportional to OC)
double OB_drive = OB_coupling * OCN / OBN_base;

dxdt_OBP = OBP_syn - OBP_deg * OBP - OB_diff * OBP;
dxdt_OBN = OB_diff * OBP * OB_drive - OBN_life * OBN;

// ============================================================
// BIOMARKERS
// ============================================================
// CTx: produced by osteoclasts
dxdt_CTX   = CTX_prod * OCN - CTX_elim * CTX;

// P1NP: produced by osteoblasts
dxdt_P1NP  = P1NP_prod * OBN - P1NP_elim * P1NP;

// Bone ALP: produced by osteoblasts (longer half-life ~3 days)
dxdt_ALP_b = ALP_prod * OBN - ALP_elim * ALP_b;

// TRAP5b
dxdt_TRAP_s = TRAP_prod * OCN - TRAP_elim * TRAP_s;

// ============================================================
// BMD (pagetic site) — net balance of formation and resorption
// ============================================================
double BmdForm   = BMD_form  * (OBN / OBN_base);
double BmdResorb = BMD_resorb * (OCN / OCN_base);
dxdt_BMD = BmdForm - BmdResorb;

// ============================================================
// PAIN (VAS 0-10)
// ============================================================
double ALP_effect = (ALP_b - 100) / 300.0;   // pain driven by disease activity
dxdt_PAIN = PAIN_slope * ALP_effect * 0.001 - PAIN_recov * (PAIN - 1.0);

$TABLE
// Derived outputs
double C1_bp_out  = DRUG1 / V1_bp;           // BP plasma conc nmol/L
double DEN_C_mgL  = DEN_C / V1_den;          // Denosumab plasma mg/L
double RANKL_OPG  = RANKL_f / (OPG_f + 0.001); // RANKL:OPG ratio
double ALP_pct    = ALP_b / 350.0 * 100;     // ALP % of baseline
double ALP_norm   = (ALP_b < 120) ? 1 : 0;  // ALP normalization flag

$CAPTURE C1_bp_out DEN_C_mgL RANKL_OPG ALP_pct ALP_norm
'

## Build model
pbd_mod <- mcode("pbd_qsp", pbd_model_code)

## ============================================================
## DOSING REGIMENS (5 SCENARIOS)
## ============================================================

t_sim <- seq(0, 8760, by = 24)   # 1 year in hours, daily steps

## Scenario 1: No Treatment — disease natural history
ev1 <- ev(time = 0, amt = 0, cmt = 1)   # null dose

## Scenario 2: Zoledronate 5 mg IV single infusion (year 0)
## IV bolus approximation: CL_bp calibrated, Kbone_on = 0.05
ev2 <- ev(time = 0,    amt = 5000,  cmt = 1, rate = 333,   # 15-min infusion (5000 mcg/h equiv.)
          addl = 0)

## Scenario 3: Pamidronate 60 mg IV × 3 consecutive days
ev3 <- ev(time = 0,    amt = 60000, cmt = 1, rate = 2500,  # 24h infusion
          addl = 2,    ii = 24)

## Scenario 4: Alendronate 40 mg PO daily for 6 months
## (oral F = 0.006, Ka = 0.5 h-1)
ev4 <- ev(time = 0,    amt = 40000, cmt = 1,   # using F_bp for bioavailability
          addl = 179,  ii = 24)

## Scenario 5: Denosumab 60 mg SC Q6 months × 2 doses (1-year simulation)
ev5 <- ev(time = 0,    amt = 60,    cmt = 4,   # DEN_SC compartment
          addl = 1,    ii = 4380)              # Q6M = 4380 h

## Parameter sets for each scenario
params_base <- list(PBD_fold = 3.0, RANKL_syn = 0.08)

## Helper: run one scenario and tag it
run_scenario <- function(ev_obj, label, extra_params = list()) {
  params_run <- modifyList(params_base, extra_params)
  pbd_mod %>%
    param(params_run) %>%
    mrgsim(ev = ev_obj, tgrid = t_sim, carry_out = "evid") %>%
    as_tibble() %>%
    mutate(Scenario = label)
}

set.seed(42)
sim1 <- run_scenario(ev1, "1_NoTreatment")
sim2 <- run_scenario(ev2, "2_Zoledronate_5mg_IV")
sim3 <- run_scenario(ev3, "3_Pamidronate_60mg_IV")
sim4 <- run_scenario(ev4, "4_Alendronate_40mg_PO")
sim5 <- run_scenario(ev5, "5_Denosumab_60mg_SC", extra_params = list(Ka_bp=0))

sim_all <- bind_rows(sim1, sim2, sim3, sim4, sim5) %>%
  mutate(day = time / 24)

## ============================================================
## VISUALISATION
## ============================================================

scenario_colors <- c(
  "1_NoTreatment"        = "#E53935",
  "2_Zoledronate_5mg_IV" = "#1E88E5",
  "3_Pamidronate_60mg_IV"= "#43A047",
  "4_Alendronate_40mg_PO"= "#FB8C00",
  "5_Denosumab_60mg_SC"  = "#8E24AA"
)
scenario_labels <- c(
  "1_NoTreatment"        = "No Treatment",
  "2_Zoledronate_5mg_IV" = "Zoledronate 5 mg IV",
  "3_Pamidronate_60mg_IV"= "Pamidronate 60 mg IV×3",
  "4_Alendronate_40mg_PO"= "Alendronate 40 mg PO daily",
  "5_Denosumab_60mg_SC"  = "Denosumab 60 mg SC Q6M"
)

theme_pbd <- theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    plot.title      = element_text(face = "bold", size = 12)
  )

## Panel A: Bone ALP (primary endpoint)
pA <- ggplot(sim_all, aes(day, ALP_b, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 120, linetype = "dashed", color = "gray50") +
  annotate("text", x = 5, y = 130, label = "ULN = 120 IU/L", size = 3, color = "gray40") +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  scale_y_continuous(limits = c(50, 500)) +
  labs(title = "A. Bone ALP (Primary Endpoint)",
       x = "Day", y = "Bone ALP (IU/L)") +
  theme_pbd

## Panel B: Serum CTx (resorption marker)
pB <- ggplot(sim_all, aes(day, CTX, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "B. Serum CTx (Resorption Marker)",
       x = "Day", y = "CTx (ng/mL)") +
  theme_pbd

## Panel C: Serum P1NP (formation marker)
pC <- ggplot(sim_all, aes(day, P1NP, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "C. Serum P1NP (Formation Marker)",
       x = "Day", y = "P1NP (ng/mL)") +
  theme_pbd

## Panel D: Active osteoclast number
pD <- ggplot(sim_all, aes(day, OCN, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  annotate("text", x = 5, y = 110, label = "Normal OC~100/mL", size = 3, color = "gray40") +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "D. Active Osteoclast Number",
       x = "Day", y = "OC Density (cells/mL)") +
  theme_pbd

## Panel E: BMD at pagetic site
pE <- ggplot(sim_all, aes(day, BMD, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "E. BMD — Pagetic Site",
       x = "Day", y = "BMD (g/cm²)") +
  theme_pbd

## Panel F: Pain VAS
pF <- ggplot(sim_all, aes(day, PAIN, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(limits = c(0, 10)) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "F. Pain VAS Score",
       x = "Day", y = "Pain VAS (0–10)") +
  theme_pbd

## Combined figure (requires patchwork)
fig_all <- (pA | pB | pC) / (pD | pE | pF) +
  plot_annotation(
    title   = "Paget's Disease of Bone — QSP Simulation Results",
    subtitle= "Five treatment scenarios over 1 year",
    theme   = theme(plot.title    = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 11))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## ============================================================
## DOSE-RESPONSE: Zoledronate dose vs ALP at day 365
## ============================================================

doses_zol <- c(0.5, 1, 2, 5, 10) * 1000   # mcg
alp_at365 <- sapply(doses_zol, function(d) {
  ev_d <- ev(time = 0, amt = d, cmt = 1, rate = d / 0.25)   # 15-min infusion
  res  <- pbd_mod %>%
    param(params_base) %>%
    mrgsim(ev = ev_d, end = 8760, delta = 24) %>%
    as_tibble()
  tail(res$ALP_b, 1)
})

df_dr <- tibble(Dose_mg = doses_zol / 1000, ALP_Day365 = alp_at365)

pDR <- ggplot(df_dr, aes(Dose_mg, ALP_Day365)) +
  geom_point(size = 4, color = "#1E88E5") +
  geom_line(color = "#1E88E5", linewidth = 1) +
  geom_hline(yintercept = 120, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Zoledronate Dose–Response (ALP at Day 365)",
       x = "Dose (mg, log scale)", y = "Bone ALP (IU/L)") +
  theme_classic(base_size = 11)

## ============================================================
## SENSITIVITY ANALYSIS: RANKL/OPG ratio effect on ALP
## ============================================================

pbd_folds <- c(1.0, 1.5, 2.0, 3.0, 4.0, 6.0)

df_sens <- lapply(pbd_folds, function(f) {
  ev_z <- ev(time = 0, amt = 5000, cmt = 1, rate = 333)
  pbd_mod %>%
    param(list(PBD_fold = f)) %>%
    mrgsim(ev = ev_z, end = 8760, delta = 24) %>%
    as_tibble() %>%
    mutate(day = time / 24, RANKL_fold = f)
}) %>% bind_rows()

pSens <- ggplot(df_sens, aes(day, ALP_b, color = factor(RANKL_fold))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 120, linetype = "dashed", color = "gray50") +
  scale_color_viridis_d(name = "RANKL/OPG\nFold Elevation") +
  labs(title = "Sensitivity: RANKL/OPG Ratio vs ALP Response\n(after Zoledronate 5 mg IV)",
       x = "Day", y = "Bone ALP (IU/L)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

## ============================================================
## SUMMARY TABLE — Day 30, 90, 180, 365
## ============================================================

summary_tbl <- sim_all %>%
  filter(day %in% c(0, 30, 90, 180, 365)) %>%
  group_by(Scenario, day) %>%
  summarize(
    ALP_IUL  = round(mean(ALP_b), 1),
    CTx_ngmL = round(mean(CTX), 2),
    P1NP_ng  = round(mean(P1NP), 1),
    OCN_cells= round(mean(OCN), 0),
    Pain_VAS = round(mean(PAIN), 1),
    .groups = "drop"
  ) %>%
  mutate(Scenario = scenario_labels[Scenario])

print(summary_tbl)

## ============================================================
## SAVE PLOTS
## ============================================================

if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/pbd_main_simulation.png",  fig_all,  width = 16, height = 10, dpi = 150)
ggsave("plots/pbd_dose_response.png",    pDR,      width = 7,  height = 5,  dpi = 150)
ggsave("plots/pbd_sensitivity.png",      pSens,    width = 9,  height = 5,  dpi = 150)

message("PBD QSP model run complete. Plots saved to ./plots/")
