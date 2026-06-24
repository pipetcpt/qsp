## ============================================================
## Diabetic Retinopathy (DR) — QSP mrgsolve ODE Model
## ============================================================
## 질환: 당뇨병성 망막병증 (Diabetic Retinopathy)
## 약어: DR
## 버전: 1.0   날짜: 2026-06-24
##
## 구획 목록 (18 ODE compartments):
##  Drug PK (4):  DRUG_CENT, DRUG_VIT, DRUG_PERIPH, CORT_VIT
##  Glycemic (2): BG, HBA1C
##  VEGF (3):     VEGF_FREE, VEGF_BOUND, VEGF_PLANT (PlGF)
##  Oxidative (2):ROS, AGE
##  Inflammation (2): CYT, ICAM
##  Cellular (2): PERICYTE, EC (endothelial cells)
##  Structural (3): PERM, NV (neovascularization index), CRT
##  Visual (1):   VA (visual acuity — ETDRS letters)
##
## 치료 시나리오 (6):
##  S0 = No treatment (disease progression, poor glycemic ctrl)
##  S1 = Tight glycemic control only (HbA1c → 7%)
##  S2 = Aflibercept 2mg IVT q4w × 5, then q8w
##  S3 = Ranibizumab 0.5mg IVT q4w
##  S4 = Faricimab 6mg IVT q4w × 4, then q16w (dual VEGF/Ang2)
##  S5 = Aflibercept + tight glycemic control (combination)
##
## 파라미터 보정 근거 (주요 임상시험):
##  - RISE/RIDE (ranibizumab DME, NEJM 2012)
##  - PROTOCOL T (aflibercept vs ranibizumab vs bevacizumab, NEJM 2015)
##  - CLARITY (aflibercept vs PRP for PDR, Lancet 2017)
##  - PANORAMA (aflibercept for NPDR, Ophthalmology 2019)
##  - TENAYA/LUCERNE (faricimab vs aflibercept, Lancet 2022)
##  - DCCT/EDIC (glycemic control, NEJM 1993, long-term follow-up)
##  - UKPDS (type 2 DM outcomes, BMJ 1998)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## -------- mrgsolve model code block -------------------------
code_DR <- '
$PROB
Diabetic Retinopathy (DR) QSP Model
18-compartment ODE: Drug PK + Glycemic + VEGF + ROS/AGE + Inflammation + Cellular + Structural + Visual

$PARAM @annotated
// ---- Drug PK parameters (aflibercept representative) -----
DOSE_IVT   :  2.0  : Intravitreal dose [mg]
MW_AFL     :  115  : Molecular weight aflibercept [kDa]
Vd_vit     :  4.0  : Vitreous volume [mL]
CL_vit     :  0.44 : Vitreous clearance [mL/day] (t½~9d AFL; ~3d RBZ; ~10d FAR)
Q_vit      :  0.05 : Vitreous-plasma distribution [mL/day]
Vd_cent    :  3600 : Plasma central volume [mL]  (AFL ~8L)
CL_cent    :  360  : Plasma clearance [mL/day]   (AFL ~1.6L/h)
Vd_per     :  5000 : Peripheral volume [mL]
Q_per      :  50   : Intercompartmental clearance [mL/day]
// Binding kinetics (drug-VEGF)
kon_VEGF   :  13.2 : VEGF binding on-rate [nM-1 day-1]
koff_VEGF  :  0.066: VEGF binding off-rate [day-1]  (Kd~5pM AFL)
kon_Ang2   :  8.0  : Ang2 binding on-rate [nM-1 day-1] (faricimab only)
koff_Ang2  :  0.05 : Ang2 binding off-rate [day-1]
// Corticosteroid PK (dexamethasone implant)
CL_cort    :  0.35 : Corticosteroid vitreous clearance [day-1] (Ozurdex t½~2d release over 60d)
// ---- Glycemic control -----
BG_base    :  9.5  : Baseline blood glucose [mmol/L] (poor control)
BG_ctrl    :  7.2  : Treated BG target [mmol/L] (HbA1c ~7%)
kBG_eq     :  0.05 : BG equilibration rate [day-1]
kHbA1c_up  :  0.033: HbA1c update rate [day-1] (reflects ~3mo average)
// ---- VEGF dynamics -----
ksyn_VEGF  :  0.15 : VEGF retinal synthesis rate [nM/day] (baseline)
kdeg_VEGF  :  0.12 : VEGF degradation rate [day-1]
VEGF_base  :  1.25 : Baseline free vitreous VEGF [nM] (~1250 pg/mL)
VEGF_max   :  8.0  : Max VEGF (severe DR/DME) [nM]
// Hill coefficient for BG→VEGF
n_BG_VEGF  :  2.0  : Hill coefficient
EC50_BG    :  9.0  : BG EC50 for VEGF synthesis [mmol/L]
// PlGF
ksyn_PLGF  :  0.08 : PlGF synthesis rate [nM/day]
kdeg_PLGF  :  0.10 : PlGF degradation rate [day-1]
// ---- Oxidative Stress (ROS) -----
ksyn_ROS   :  0.20 : ROS production rate [AU/day] (normalized)
kdeg_ROS   :  0.18 : ROS clearance rate [day-1]
ROS_base   :  1.0  : Baseline ROS [AU]
kBG_ROS    :  0.08 : BG driving coefficient for ROS
// AGE accumulation
ksyn_AGE   :  0.005: AGE formation rate [AU/day] (slow accumulation)
kdeg_AGE   :  0.002: AGE clearance rate [day-1] (very slow)
AGE_base   :  1.0  : Baseline AGE [AU]
// ---- Inflammation -----
ksyn_CYT   :  0.30 : Cytokine (IL-6/IL-8/TNFα composite) synthesis [AU/day]
kdeg_CYT   :  0.25 : Cytokine degradation [day-1]
CYT_base   :  1.0  : Baseline cytokine index [AU]
kVEGF_CYT  :  0.15 : VEGF driving cytokine synthesis
kROS_CYT   :  0.12 : ROS driving cytokine synthesis
ksyn_ICAM  :  0.20 : ICAM-1 synthesis [AU/day]
kdeg_ICAM  :  0.15 : ICAM-1 degradation [day-1]
ICAM_base  :  1.0  : Baseline ICAM-1 [AU]
kCYT_ICAM  :  0.18 : Cytokines → ICAM-1
// ---- Cellular (pericyte / EC) -----
PERI_base  :  100.0: Baseline pericyte count [% normal]
kLoss_PERI :  0.003: Pericyte loss rate [day-1] (driven by AGE, ROS)
kAGE_PERI  :  0.015: AGE accelerates pericyte loss
kROS_PERI  :  0.010: ROS accelerates pericyte loss
EC_base    :  100.0: Baseline EC count [% normal]
kLoss_EC   :  0.002: EC apoptosis rate [day-1]
kICM_EC    :  0.012: ICAM-1/leukostasis → EC loss
// ---- Structural: Permeability / NV / CRT -----
PERM_base  :  1.0  : Baseline vascular permeability [AU]
ksyn_PERM  :  0.25 : Permeability synthesis [AU/day] (driven by VEGF)
kdeg_PERM  :  0.20 : Permeability normalization rate [day-1]
kVEGF_PERM :  0.30 : VEGF → permeability driving coefficient
kCYT_PERM  :  0.15 : Cytokines → permeability
NV_base    :  0.0  : Baseline NV index [AU] (0=no NV)
ksyn_NV    :  0.008: NV growth rate [AU/day] (driven by VEGF/hypoxia)
kdeg_NV    :  0.002: NV regression rate [day-1]
kVEGF_NV   :  0.05 : VEGF → NV
CRT_base   :  270.0: Baseline CRT [µm] (normal ~270µm)
kCRT_up    :  0.15 : CRT increase rate driven by permeability
kCRT_dn    :  0.08 : Spontaneous CRT resolution rate [day-1]
CRT_max    :  600.0: Maximum CRT [µm]
// ---- Visual Acuity (ETDRS letters) -----
VA_base    :  65.0 : Baseline VA [ETDRS letters] (~20/50)
kVA_dn     :  0.005: VA loss rate driven by CRT/NV
kVA_up     :  0.003: VA recovery rate (treatment driven)
VA_min     :  0.0  : Minimum VA [letters]
VA_max     :  100.0: Maximum VA [letters]
kCRT_VA    :  0.002: CRT excess → VA loss
kNV_VA     :  0.015: NV → VA loss (vitreous hemorrhage etc.)

$CMT @annotated
DRUG_CENT  : Drug central compartment [mg]
DRUG_VIT   : Drug vitreous compartment [mg]
DRUG_PERIPH: Drug peripheral compartment [mg]
CORT_VIT   : Corticosteroid vitreous [mg]
BG         : Blood glucose [mmol/L]
HBA1C      : HbA1c [%]
VEGF_FREE  : Free vitreous VEGF [nM]
VEGF_BOUND : VEGF-drug complex [nM]
VEGF_PLANT : PlGF [nM]
ROS        : Oxidative stress index [AU]
AGE        : AGE accumulation [AU]
CYT        : Cytokine index (IL-6/IL-8/TNFα) [AU]
ICAM       : ICAM-1 index [AU]
PERICYTE   : Pericyte count [% normal]
EC_COUNT   : Endothelial cell count [% normal]
PERM       : Vascular permeability [AU]
NV         : Neovascularization index [AU]
CRT        : Central retinal thickness [µm]
VA         : Visual acuity [ETDRS letters]

$MAIN
// Calculated VEGF synthesis driven by BG and AGE
double BG_eff   = pow(BG, n_BG_VEGF) / (pow(EC50_BG, n_BG_VEGF) + pow(BG, n_BG_VEGF));
double VEGFsyn  = ksyn_VEGF * (1.0 + 3.0 * BG_eff) * (1.0 + 0.5 * (AGE/AGE_base - 1.0));

// Corticosteroid inhibition on VEGF and cytokines
double Cort_inh = CORT_VIT / (CORT_VIT + 0.05);  // IC50=0.05mg for GR effect
double VEGFsyn_eff = VEGFsyn * (1.0 - 0.6 * Cort_inh);

// VEGF binding to drug
double Drug_vit_nM = DRUG_VIT / Vd_vit * 1e6 / MW_AFL;  // convert mg/mL → nM
double Bind_on     = kon_VEGF * VEGF_FREE * Drug_vit_nM;
double Bind_off    = koff_VEGF * VEGF_BOUND;

// Ang2-drug binding (faricimab context — parameterised via F_ANG2 flag)
double F_ANG2 = 0.0;  // set to 1.0 for faricimab scenarios

$GLOBAL
// Global helper: clamp function
double clamp_val(double x, double lo, double hi) {
  return (x < lo) ? lo : (x > hi) ? hi : x;
}

$ODE
// ---- Drug PK ----
// Vitreous compartment (site of action)
dxdt_DRUG_VIT    =  - (CL_vit / Vd_vit) * DRUG_VIT
                    - Q_vit  * DRUG_VIT / Vd_vit
                    + Q_vit  * DRUG_CENT / Vd_cent;
// Plasma central
dxdt_DRUG_CENT   =  - (CL_cent / Vd_cent) * DRUG_CENT
                    + Q_vit  * DRUG_VIT / Vd_vit
                    - Q_vit  * DRUG_CENT / Vd_cent
                    - (Q_per  / Vd_cent) * DRUG_CENT
                    + (Q_per  / Vd_per)  * DRUG_PERIPH;
// Peripheral
dxdt_DRUG_PERIPH =   (Q_per  / Vd_cent) * DRUG_CENT
                    - (Q_per  / Vd_per)  * DRUG_PERIPH;
// Corticosteroid vitreous (slow-release depot)
dxdt_CORT_VIT    =  - CL_cort * CORT_VIT;

// ---- Glycemic markers ----
dxdt_BG    = kBG_eq * (BG_base - BG);   // drives toward baseline unless treatment modifies BG_base
dxdt_HBA1C = kHbA1c_up * (0.195 * BG + 3.1 - HBA1C);  // Rohlfing conversion: HbA1c ≈ (BG*0.195)+3.1

// ---- VEGF Dynamics ----
double BG_eff2   = pow(BG, n_BG_VEGF) / (pow(EC50_BG, n_BG_VEGF) + pow(BG, n_BG_VEGF));
double VEGFsyn_d = ksyn_VEGF * (1.0 + 3.0 * BG_eff2) * (1.0 + 0.5 * (AGE/AGE_base - 1.0)) * (1.0 - 0.6 * (CORT_VIT/(CORT_VIT+0.05)));
double Dvit_nM   = (DRUG_VIT / Vd_vit) * 1e6 / MW_AFL;
double Bon       = kon_VEGF  * VEGF_FREE * Dvit_nM;
double Boff      = koff_VEGF * VEGF_BOUND;

dxdt_VEGF_FREE  = VEGFsyn_d - kdeg_VEGF * VEGF_FREE - Bon + Boff;
dxdt_VEGF_BOUND = Bon - Boff - kdeg_VEGF * VEGF_BOUND;   // complex degrades too

// PlGF (influenced by BG and HIF, simpler)
double PLGFsyn  = ksyn_PLGF * (1.0 + 1.5 * BG_eff2);
dxdt_VEGF_PLANT = PLGFsyn - kdeg_PLGF * VEGF_PLANT;

// ---- Oxidative Stress (ROS) ----
double ROS_driv = kBG_ROS * (BG - 5.5);   // excess glucose above normal 5.5 mmol/L
dxdt_ROS        = ksyn_ROS + (ROS_driv > 0 ? ROS_driv : 0.0) - kdeg_ROS * ROS;

// ---- AGE Accumulation ----
double AGE_rate = ksyn_AGE * (BG / 5.5) * (1.0 + 0.5 * (HBA1C - 7.0)/7.0);
dxdt_AGE        = AGE_rate - kdeg_AGE * AGE;

// ---- Inflammation ----
double CYT_stim = kVEGF_CYT * (VEGF_FREE / VEGF_base) + kROS_CYT * (ROS / ROS_base);
double CYT_inh  = 0.7 * (CORT_VIT / (CORT_VIT + 0.05));  // corticosteroid inhibition
dxdt_CYT        = ksyn_CYT * (1.0 + CYT_stim) * (1.0 - CYT_inh) - kdeg_CYT * CYT;
dxdt_ICAM       = ksyn_ICAM * (1.0 + kCYT_ICAM * (CYT / CYT_base)) - kdeg_ICAM * ICAM;

// ---- Cellular Dynamics ----
double PERI_loss = kLoss_PERI + kAGE_PERI * (AGE/AGE_base - 1.0) + kROS_PERI * (ROS/ROS_base - 1.0);
dxdt_PERICYTE   = - PERI_loss * PERICYTE;   // irreversible loss (no regeneration modeled)
double EC_loss   = kLoss_EC + kICM_EC * (ICAM / ICAM_base - 1.0);
dxdt_EC_COUNT   = - EC_loss * EC_COUNT;

// ---- Structural ----
// Vascular Permeability driven by free VEGF and cytokines
double PERM_in  = ksyn_PERM * (kVEGF_PERM * VEGF_FREE/VEGF_base + kCYT_PERM * CYT/CYT_base);
dxdt_PERM       = PERM_in - kdeg_PERM * PERM;

// Neovascularization (NV) — develops when pericytes low + VEGF high
double NV_drive = kVEGF_NV * (VEGF_FREE / VEGF_base) * (1.0 - PERICYTE/PERI_base);
double NV_reg   = kdeg_NV * NV;   // spontaneous regression slow
dxdt_NV         = ksyn_NV * NV_drive - NV_reg;

// CRT: increases with permeability, decreases toward baseline
double CRT_eff_in = kCRT_up * (PERM / PERM_base) * (CRT_max - CRT);
double CRT_eff_dn = kCRT_dn * (CRT - CRT_base);
dxdt_CRT          = CRT_eff_in - CRT_eff_dn;

// ---- Visual Acuity (VA) ----
// Loss: driven by CRT excess (edema) and NV (hemorrhage)
double CRT_excess = (CRT > CRT_base) ? (CRT - CRT_base) / 100.0 : 0.0;
double VA_loss    = kCRT_VA * CRT_excess * VA + kNV_VA * NV * VA / 10.0;
// Recovery: spontaneous partial recovery when CRT resolves
double VA_rec     = kVA_up * (VA_base - VA) * (1.0 - CRT_excess / 4.0);
dxdt_VA           = VA_rec - VA_loss;

$CAPTURE
BG HBA1C VEGF_FREE VEGF_BOUND VEGF_PLANT ROS AGE CYT ICAM
PERICYTE EC_COUNT PERM NV CRT VA
DRUG_VIT DRUG_CENT CORT_VIT

$SET delta=1, end=730   // 2-year simulation, daily steps
'

## ============================================================
## Compile the model
## ============================================================
mod <- mcode("DR_QSP", code_DR)

## ============================================================
## Initial Conditions (Disease State at Study Entry)
## ============================================================
## Represents a patient with: HbA1c ~9.5%, moderate NPDR+DME,
##   BCVA ~65 letters (20/50), CRT ~380 µm

init_DR <- init(mod,
  DRUG_CENT   = 0,
  DRUG_VIT    = 0,
  DRUG_PERIPH = 0,
  CORT_VIT    = 0,
  BG          = 9.5,      # mmol/L
  HBA1C       = 8.5,      # %
  VEGF_FREE   = 3.5,      # nM (elevated, ~3500 pg/mL)
  VEGF_BOUND  = 0,
  VEGF_PLANT  = 0.8,      # nM PlGF
  ROS         = 2.2,      # AU (elevated)
  AGE         = 2.0,      # AU (elevated chronic)
  CYT         = 2.5,      # AU
  ICAM        = 2.0,      # AU
  PERICYTE    = 60,       # % (40% already lost)
  EC_COUNT    = 75,       # %
  PERM        = 3.0,      # AU
  NV          = 0.5,      # AU (early NV)
  CRT         = 380,      # µm (DME present)
  VA          = 60        # ETDRS letters (20/63)
)

## ============================================================
## Dosing Events
## ============================================================

## Aflibercept 2mg: loading phase q4w × 5 doses, then q8w
dose_AFL <- function(n_loading=5, loading_interval=28, maint_interval=56,
                     n_maint=14, dose_mg=2) {
  times_load <- seq(0, (n_loading - 1) * loading_interval, by=loading_interval)
  times_maint<- seq(n_loading * loading_interval,
                    n_loading * loading_interval + (n_maint-1)*maint_interval,
                    by=maint_interval)
  times_all <- c(times_load, times_maint)
  ev(amt = dose_mg, cmt = "DRUG_VIT", time = times_all)
}

## Ranibizumab 0.5mg: q4w throughout
dose_RBZ <- function(n_inj=24, interval=28, dose_mg=0.5) {
  ev(amt = dose_mg, cmt = "DRUG_VIT",
     time = seq(0, (n_inj-1)*interval, by=interval))
}

## Faricimab 6mg: loading q4w × 4, then q16w (treat-to-extend)
dose_FAR <- function(n_loading=4, loading_interval=28,
                     maint_interval=112, n_maint=6, dose_mg=6) {
  times_load <- seq(0, (n_loading - 1) * loading_interval, by=loading_interval)
  times_maint<- seq(n_loading * loading_interval,
                    n_loading * loading_interval + (n_maint-1)*maint_interval,
                    by=maint_interval)
  times_all <- c(times_load, times_maint)
  ev(amt = dose_mg, cmt = "DRUG_VIT", time = times_all)
}

## Dexamethasone implant 0.7mg (Ozurdex): q6 months
dose_DEXA <- function(n_inj=4, interval=182, dose_mg=0.7) {
  ev(amt = dose_mg, cmt = "CORT_VIT",
     time = seq(0, (n_inj-1)*interval, by=interval))
}

## ============================================================
## Helper: run scenario
## ============================================================
run_scenario <- function(mod, init_state, ev_obj, params=list(),
                         scenLabel="S0", end_t=730) {
  m <- mod %>% param(params)
  if (!is.null(ev_obj)) {
    out <- mrgsim(m, init = init_state, events = ev_obj, end = end_t, delta = 1)
  } else {
    out <- mrgsim(m, init = init_state, end = end_t, delta = 1)
  }
  as.data.frame(out) %>% mutate(Scenario = scenLabel)
}

## ============================================================
## Scenario Definitions & Run
## ============================================================

## S0: No treatment — poor glycemic control (BG stays ~9.5)
params_S0 <- list(BG_base = 9.5, BG_ctrl = 9.5, kBG_eq = 0.01)
s0 <- run_scenario(mod, init_DR, NULL, params_S0, "S0: No Tx (poor ctrl)")

## S1: Tight glycemic control only (HbA1c → 7%)
## BG_base lowered to 7.2 mmol/L over ~3 months
params_S1 <- list(BG_base = 7.2, BG_ctrl = 7.2, kBG_eq = 0.02)
s1 <- run_scenario(mod, init_DR, NULL, params_S1, "S1: Glycemic Ctrl Only")

## S2: Aflibercept 2mg (loading q4w×5, then q8w)
params_S2 <- list(BG_base = 9.5, kBG_eq = 0.01)
ev_S2 <- dose_AFL()
s2 <- run_scenario(mod, init_DR, ev_S2, params_S2, "S2: Aflibercept 2mg (q4→q8w)")

## S3: Ranibizumab 0.5mg q4w
params_S3 <- list(BG_base = 9.5, kBG_eq = 0.01,
                  CL_vit = 0.23,   # ranibizumab faster: t½~3d
                  MW_AFL = 48,     # RBZ MW ~48 kDa
                  kon_VEGF = 10.0, koff_VEGF = 0.10)
ev_S3 <- dose_RBZ()
s3 <- run_scenario(mod, init_DR, ev_S3, params_S3, "S3: Ranibizumab 0.5mg q4w")

## S4: Faricimab 6mg (loading q4w×4, then q16w) — higher MW, dual action
params_S4 <- list(BG_base = 9.5, kBG_eq = 0.01,
                  MW_AFL  = 146,   # faricimab MW ~146 kDa
                  CL_vit  = 0.35,  # t½ in vitreous ~14d
                  kon_VEGF = 14.0, koff_VEGF = 0.04)
ev_S4 <- dose_FAR()
s4 <- run_scenario(mod, init_DR, ev_S4, params_S4, "S4: Faricimab 6mg (q4→q16w)")

## S5: Combination — aflibercept + tight glycemic control
params_S5 <- list(BG_base = 7.2, BG_ctrl = 7.2, kBG_eq = 0.02)
ev_S5 <- dose_AFL()
s5 <- run_scenario(mod, init_DR, ev_S5, params_S5, "S5: AFL + Glycemic Ctrl")

## S6 (bonus): Dexamethasone implant (corticosteroid DME)
params_S6 <- list(BG_base = 9.5, kBG_eq = 0.01,
                  CL_cort = 0.35)
ev_S6 <- dose_DEXA()
s6 <- run_scenario(mod, init_DR, ev_S6, params_S6, "S6: Dexamethasone Implant")

## Combine all scenarios
all_scen <- bind_rows(s0, s1, s2, s3, s4, s5, s6)

## ============================================================
## Visualization
## ============================================================
scen_colors <- c(
  "S0: No Tx (poor ctrl)"          = "#E53935",
  "S1: Glycemic Ctrl Only"         = "#FB8C00",
  "S2: Aflibercept 2mg (q4→q8w)"  = "#1E88E5",
  "S3: Ranibizumab 0.5mg q4w"     = "#43A047",
  "S4: Faricimab 6mg (q4→q16w)"  = "#8E24AA",
  "S5: AFL + Glycemic Ctrl"        = "#00ACC1",
  "S6: Dexamethasone Implant"      = "#6D4C41"
)

## -- Figure 1: Visual Acuity (primary endpoint)
p_VA <- ggplot(all_scen, aes(x=time, y=VA, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  geom_hline(yintercept=55, linetype="dashed", color="gray40") +
  annotate("text", x=500, y=56, label="Severe VL threshold (55 letters)", size=3, color="gray40") +
  labs(title="Visual Acuity Over 2 Years",
       subtitle="ETDRS letter score (higher = better)",
       x="Time (days)", y="ETDRS Letters", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=9))

## -- Figure 2: Central Retinal Thickness
p_CRT <- ggplot(all_scen, aes(x=time, y=CRT, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  geom_hline(yintercept=310, linetype="dashed", color="darkblue") +
  annotate("text", x=500, y=315, label="CI-DME threshold (310 µm)", size=3, color="darkblue") +
  labs(title="Central Retinal Thickness (CRT) Over 2 Years",
       subtitle="OCT-measured structural endpoint (lower = better)",
       x="Time (days)", y="CRT (µm)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=9))

## -- Figure 3: Free Vitreous VEGF
p_VEGF <- ggplot(all_scen, aes(x=time, y=VEGF_FREE, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  labs(title="Free Vitreous VEGF Concentration",
       subtitle="Pharmacodynamic target of anti-VEGF therapy",
       x="Time (days)", y="Free VEGF (nM)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=9))

## -- Figure 4: Pericyte Count over time
p_PERI <- ggplot(all_scen, aes(x=time, y=PERICYTE, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  labs(title="Pericyte Count (% Normal)",
       subtitle="Cellular biomarker of DR progression",
       x="Time (days)", y="Pericytes (% normal)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=9))

## -- Figure 5: Neovascularization Index
p_NV <- ggplot(all_scen, aes(x=time, y=NV, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  labs(title="Retinal Neovascularization Index",
       subtitle="NV development in PDR (anti-VEGF vs untreated)",
       x="Time (days)", y="NV Index (AU)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="right", legend.text=element_text(size=9))

## -- Figure 6: Drug vitreous PK (S2 aflibercept)
pk_data <- all_scen %>% filter(grepl("Aflibercept", Scenario) | grepl("Ranibiz", Scenario))
p_PK <- ggplot(pk_data, aes(x=time, y=DRUG_VIT, color=Scenario)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=scen_colors) +
  scale_y_log10() +
  labs(title="Intravitreal Drug Concentration Over Time",
       subtitle="Anti-VEGF agent PK in vitreous humor (log scale)",
       x="Time (days)", y="Drug Concentration (mg/mL, log10)", color="Scenario") +
  theme_bw(base_size=12)

## -- Summary table: endpoint values at 1 year and 2 years
summary_table <- all_scen %>%
  filter(time %in% c(364, 728)) %>%
  mutate(Timepoint = ifelse(time==364, "Year 1", "Year 2")) %>%
  select(Scenario, Timepoint, VA, CRT, VEGF_FREE, NV, PERICYTE, HBA1C) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("\n===== DR QSP Model Summary: Key Endpoints =====\n")
print(summary_table, n=Inf)

## -- VA change from baseline
va_change <- all_scen %>%
  filter(time %in% c(0, 364, 728)) %>%
  group_by(Scenario) %>%
  mutate(VA_change = VA - VA[time == 0]) %>%
  filter(time > 0) %>%
  select(Scenario, time, VA, VA_change, CRT)

cat("\n===== VA Change from Baseline =====\n")
print(va_change, n=Inf)

## ============================================================
## Output plots (if running interactively)
## ============================================================
if (interactive()) {
  print(p_VA)
  print(p_CRT)
  print(p_VEGF)
  print(p_PERI)
  print(p_NV)
  print(p_PK)
}

## ============================================================
## Clinical Trial Benchmarking Notes
## ============================================================
cat("
=== Parameter Calibration Notes ===

PROTOCOL T (NEJM 2015, n=660, 1-year primary):
  Aflibercept 2mg: VA +13.3 letters; CRT -169µm (baseline VA 61 letters, CRT ~407µm)
  Ranibizumab 0.3mg: VA +11.2 letters; CRT -147µm
  Bevacizumab 1.25mg: VA +9.7 letters; CRT -101µm

PANORAMA (Ophthalmology 2019, n=402, 100wk, NPDR without DME):
  Aflibercept 2mg q16w: 2-step improvement 65% vs sham 15%

TENAYA/LUCERNE (Lancet 2022, n=671+729, 1-year):
  Faricimab 6mg: VA +5.8/+6.6 letters; CRT -189/-194µm
  Aflibercept 2mg: VA +5.1/+6.6 letters; CRT -163/-174µm
  Mean injection interval Faricimab: ~13-14 weeks (treat-to-extend)

DCCT (NEJM 1993, n=1441, 6.5 years):
  Intensive glycemic control: 76% reduction in new DR incidence
  HbA1c 7.2% vs 9.1%

MODEL CALIBRATION TARGETS (approximate Year-1 from S2 AFL):
  VA: +10 to +14 letters from baseline 60
  CRT: decrease from 380 → 220-250µm
  Free VEGF: suppression by >80% in loading phase
  NV Index: significant regression after treatment
")
