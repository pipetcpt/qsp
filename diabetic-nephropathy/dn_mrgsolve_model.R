## =============================================================================
## Diabetic Nephropathy (DN) QSP Model — mrgsolve Implementation
## =============================================================================
## Disease: Diabetic Nephropathy (ICD-10 N08/E11.21)
## Purpose: Simulate progressive CKD in type 2 diabetes with multi-drug PK/PD
##
## KEY PATHWAYS MODELLED:
##  1. Glucose metabolism / hyperglycemia
##  2. RAAS (renin-angiotensin-aldosterone system)
##  3. TGF-β / ECM fibrosis
##  4. Oxidative stress
##  5. Podocyte damage → proteinuria
##  6. Tubular injury → interstitial fibrosis → GFR decline
##
## DRUG PK/PD (4 agents):
##  Drug 1 — ACE Inhibitor   (enalapril prototype)
##  Drug 2 — ARB             (losartan prototype)
##  Drug 3 — SGLT2 Inhibitor (empagliflozin prototype)
##  Drug 4 — Finerenone      (non-steroidal MRA)
##
## TREATMENT SCENARIOS (5 defined below):
##  S0: No treatment
##  S1: ACEi monotherapy
##  S2: ARB monotherapy
##  S3: SGLT2i monotherapy
##  S4: ACEi + SGLT2i dual
##  S5: SGLT2i + Finerenone (FIDELIO-DKD inspired)
##  S6: ACEi + SGLT2i + Finerenone (triple)
##
## CALIBRATION NOTES:
##  • GFR decline rate ~3-4 mL/min/yr in untreated DKD (Perkins 2003 NEJM)
##  • ACEi reduces UACR ~30-35% (Lewis 1993 NEJM; BENEDICT trial)
##  • ARB reduces UACR ~25-30% (RENAAL, IDNT 2001 NEJM)
##  • SGLT2i reduces eGFR slope ~2 mL/min/yr (EMPA-REG OUTCOME, CREDENCE)
##  • Finerenone reduces UACR ~31%, eGFR slope -0.9 mL/min/yr (FIDELIO-DKD NEJM 2020)
##  • Combined SGLT2i + MRA may provide additive renoprotection (CONFIDENCE trial)
##
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## -----------------------------------------------------------------------------
## mrgsolve MODEL CODE
## -----------------------------------------------------------------------------

code <- '
$PROB
Diabetic Nephropathy QSP Model v1.0
19 ODE compartments | 4 drug PK | Disease PD | Clinical Endpoints

$PARAM
// ---- Drug PK Parameters ----
// ACEi (enalapril-like)
ka_acei  = 1.2    // h-1 oral absorption rate
Vd_acei  = 25     // L   central volume
CL_acei  = 5.5    // L/h clearance
// ARB (losartan-like)
ka_arb   = 0.9
Vd_arb   = 32
CL_arb   = 7.8
// SGLT2i (empagliflozin-like)
ka_sglt2 = 1.5
Vd_sglt2 = 74
CL_sglt2 = 9.2
// Finerenone
ka_fine  = 1.1
Vd_fine  = 52
CL_fine  = 6.3

// ---- Drug PD (Emax) ----
Emax_acei   = 0.90  // max ACE inhibition
EC50_acei   = 2.5   // ng/mL half-max conc.
Emax_arb    = 0.85
EC50_arb    = 80    // ng/mL
Emax_sglt2  = 0.85  // max SGLT2 blockade
EC50_sglt2  = 15    // ng/mL
Emax_fine   = 0.88  // max MR blockade
EC50_fine   = 120   // ng/mL

// ---- Glucose / Metabolic ----
BG_base     = 8.5   // mmol/L baseline BG (poorly controlled DM2)
BG_min      = 4.5   // mmol/L normal fasting
kBG_in      = 0.05  // h-1  glucose production rate const.
kBG_out     = 0.008 // h-1  insulin-mediated clearance
HbA1c_base  = 8.2   // % baseline
kHbA1c      = 0.003 // day-1  HbA1c equilibration rate
AGE_base    = 1.0   // AU  baseline AGE level
kAGE_in     = 0.01  // AU/day  AGE formation rate (glucose-driven)
kAGE_out    = 0.002 // day-1   AGE clearance
AGE_thresh  = 0.5   // AU  threshold above baseline

// ---- RAAS ----
AngII_base  = 1.0   // AU  (normalized)
kAngII_in   = 0.15  // AU/day  AngII production
kAngII_out  = 0.12  // day-1   AngII clearance
Pglo_base   = 48    // mmHg  normal intraglomerular pressure
Pglo_max    = 65    // mmHg  max with uncontrolled AngII
kPglo_AngII = 8.0   // mmHg/(AU) AngII -> Pglo sensitivity

// ---- TGF-beta / Fibrosis ----
TGF_base    = 1.0   // AU
kTGF_in     = 0.08  // day-1
kTGF_out    = 0.06  // day-1
TGF_AngII   = 0.35  // AU/AU AngII-driven TGF induction
TGF_AGE     = 0.20  // AU/AU AGE-driven
TGF_ROS     = 0.25  // AU/AU ROS-driven
kECM_in     = 0.04  // day-1  ECM accumulation
kECM_out    = 0.01  // day-1  ECM turnover
ECM_base    = 1.0   // AU

// ---- Oxidative Stress ----
ROS_base    = 1.0   // AU
kROS_in     = 0.12
kROS_out    = 0.10
ROS_BG      = 0.15  // glucose contribution to ROS
ROS_AngII   = 0.20  // AngII contribution
ROS_TGF     = 0.10  // TGF contribution (Nox4)

// ---- Podocyte Damage ----
Pod_base    = 1.0   // AU  (1=intact)
kPod_loss   = 0.003 // day-1 podocyte loss rate
Pod_TGF     = 0.8   // sensitivity to TGF
Pod_AngII   = 0.6   // sensitivity to AngII
Pod_ROS     = 0.5   // sensitivity to ROS
Pod_Pglo    = 0.004 // mmHg-1 day-1  mechanical stress
Pod_min     = 0.1   // AU  minimum podocyte density

// ---- Proteinuria (UACR) ----
UACR_base   = 300   // mg/g  baseline (overt proteinuria)
kUACR_Pod   = 200   // mg/g per AU pod. damage
kUACR_Pglo  = 3.0   // mg/g per mmHg excess pressure
kUACR_out   = 0.05  // day-1  equilibration

// ---- Tubular Compartment ----
Tub_base    = 1.0   // AU (1=healthy)
kTub_loss   = 0.002 // day-1  tubular injury accumulation
Tub_UACR    = 0.001 // AU per mg/g  protein toxicity
Tub_hypoxia = 0.3   // hypoxia contribution to tubular loss
kTub_out    = 0.005 // day-1  tubular recovery (limited)

// ---- Interstitial Fibrosis ----
Fib_base    = 0.1   // AU (10% fibrosis at baseline)
kFib_in     = 0.012 // day-1
kFib_out    = 0.002 // day-1
Fib_TGF     = 0.6   // TGF contribution
Fib_Tub     = 0.4   // tubular injury contribution

// ---- GFR / Renal Function ----
GFR_0       = 55    // mL/min/1.73m2  baseline eGFR (CKD G3a)
GFR_min     = 5     // mL/min/1.73m2  ESKD threshold
kGFR_Fib    = 0.06  // mL/min per AU fibrosis decline
kGFR_ECM    = 0.04  // mL/min per AU ECM decline
kGFR_Pglo   = 0.003 // mL/min per mmHg mechanical
// SGLT2i acute GFR dip (hemodynamic; see CREDENCE)
SGLT2_GFR_dip = 3.0 // mL/min acute dip on SGLT2i start

// ---- BP ----
SBP_base    = 145   // mmHg
DBP_base    = 88    // mmHg
kSBP_AngII  = 10    // mmHg per AU excess AngII
kSBP_Na     = 5     // mmHg per AU Na retention (Aldosterone)

$CMT
// Drug PK compartments (8)
GI_acei  CENT_acei   // [1-2] ACEi GI + Central
GI_arb   CENT_arb    // [3-4] ARB
GI_sglt2 CENT_sglt2  // [5-6] SGLT2i
GI_fine  CENT_fine   // [7-8] Finerenone

// Disease state compartments (11)
BG        // [9]  Blood glucose (mmol/L)
AGE_cmpt  // [10] AGE accumulation (AU)
AngII_cmpt // [11] AngII (AU)
TGF_cmpt  // [12] TGF-β (AU)
ROS_cmpt  // [13] Oxidative stress (AU)
ECM_cmpt  // [14] ECM / glomerulosclerosis (AU)
Pod_cmpt  // [15] Podocyte integrity (AU)
UACR_cmpt // [16] UACR (mg/g)
Tub_cmpt  // [17] Tubular integrity (AU)
Fib_cmpt  // [18] Interstitial fibrosis (AU)
GFR_cmpt  // [19] eGFR (mL/min/1.73m2)

$INIT
GI_acei=0 CENT_acei=0
GI_arb=0  CENT_arb=0
GI_sglt2=0 CENT_sglt2=0
GI_fine=0 CENT_fine=0
BG=8.5
AGE_cmpt=1.0
AngII_cmpt=1.0
TGF_cmpt=1.0
ROS_cmpt=1.0
ECM_cmpt=1.0
Pod_cmpt=1.0
UACR_cmpt=300
Tub_cmpt=1.0
Fib_cmpt=0.1
GFR_cmpt=55

$ODE
// ============================================================
// DRUG PK — first-order 2-cpt oral (1 cpt GI + 1 cpt central)
// ============================================================
// ACEi
dxdt_GI_acei  = -ka_acei  * GI_acei;
dxdt_CENT_acei = ka_acei  * GI_acei - (CL_acei/Vd_acei) * CENT_acei;

// ARB
dxdt_GI_arb  = -ka_arb   * GI_arb;
dxdt_CENT_arb = ka_arb   * GI_arb  - (CL_arb/Vd_arb)   * CENT_arb;

// SGLT2i
dxdt_GI_sglt2  = -ka_sglt2 * GI_sglt2;
dxdt_CENT_sglt2 = ka_sglt2 * GI_sglt2 - (CL_sglt2/Vd_sglt2) * CENT_sglt2;

// Finerenone
dxdt_GI_fine  = -ka_fine  * GI_fine;
dxdt_CENT_fine = ka_fine  * GI_fine  - (CL_fine/Vd_fine)  * CENT_fine;

// ============================================================
// DRUG EFFECT (Emax, inhibitory)
// ============================================================
// ACEi: inhibits ACE → reduces AngII production
double Cp_acei   = CENT_acei / Vd_acei * 1000; // ng/mL (rough unit)
double E_acei    = Emax_acei  * Cp_acei  / (EC50_acei  + Cp_acei);

// ARB: blocks AT1R → reduces AngII-mediated effects
double Cp_arb    = CENT_arb  / Vd_arb  * 1000;
double E_arb     = Emax_arb   * Cp_arb   / (EC50_arb   + Cp_arb);

// SGLT2i: blocks SGLT2 → glucosuria → BG reduction, tubular O2 demand
double Cp_sglt2  = CENT_sglt2/ Vd_sglt2*1000;
double E_sglt2   = Emax_sglt2 * Cp_sglt2 / (EC50_sglt2 + Cp_sglt2);

// Finerenone: blocks MR → reduces aldosterone-driven fibrosis/inflammation
double Cp_fine   = CENT_fine / Vd_fine *1000;
double E_fine    = Emax_fine  * Cp_fine  / (EC50_fine  + Cp_fine);

// Combined RAAS blockade
double RAAS_block = 1 - (1-E_acei)*(1-E_arb);  // combined (no double count)

// ============================================================
// BLOOD GLUCOSE  (units: mmol/L)
// ============================================================
double BG_excess  = (BG > BG_min) ? (BG - BG_min) : 0;
double BG_in      = kBG_in * BG_base;        // constant glucose load
double BG_out     = kBG_out * BG;
double BG_sglt2   = E_sglt2 * 0.015 * BG;    // SGLT2i-mediated glucosuria
dxdt_BG = BG_in - BG_out - BG_sglt2;

// ============================================================
// AGE (Advanced Glycation End-products)
// ============================================================
double AGE_drive  = kAGE_in * (BG_excess / (BG_base - BG_min)); // scaled
dxdt_AGE_cmpt = AGE_drive - kAGE_out * AGE_cmpt;

// ============================================================
// AngII — RAAS (AU)
// ============================================================
double AngII_in  = kAngII_in * (1 - RAAS_block);
double AngII_out = kAngII_out * AngII_cmpt;
dxdt_AngII_cmpt = AngII_in - AngII_out;

// ============================================================
// TGF-β (AU)
// ============================================================
double TGF_drive_AngII = TGF_AngII * (AngII_cmpt - 1) * (1 - RAAS_block) * (1 - E_fine);
double TGF_drive_AGE   = TGF_AGE   * (AGE_cmpt   - 1);
double TGF_drive_ROS   = TGF_ROS   * (ROS_cmpt   - 1);
double TGF_sglt2_red   = E_sglt2 * 0.3 * (TGF_cmpt - TGF_base); // SGLT2i direct effect
double TGF_glp1_red    = 0.0;   // placeholder if GLP-1 RA added

double TGF_in  = kTGF_in * TGF_base
               + (TGF_drive_AngII > 0 ? TGF_drive_AngII : 0)
               + (TGF_drive_AGE   > 0 ? TGF_drive_AGE   : 0)
               + (TGF_drive_ROS   > 0 ? TGF_drive_ROS   : 0);
double TGF_out = kTGF_out * TGF_cmpt + TGF_sglt2_red;
dxdt_TGF_cmpt = TGF_in - TGF_out;

// ============================================================
// Oxidative Stress (AU)
// ============================================================
double ROS_in  = kROS_in  + ROS_BG  * BG_excess/3.0
               + ROS_AngII * (AngII_cmpt - 1)
               + ROS_TGF  * (TGF_cmpt - 1);
double ROS_out = kROS_out * ROS_cmpt;
dxdt_ROS_cmpt = ROS_in - ROS_out;

// ============================================================
// ECM / Glomerulosclerosis (AU)
// ============================================================
double ECM_in  = kECM_in * TGF_cmpt * AngII_cmpt * (1 - E_fine * 0.5);
double ECM_out = kECM_out * ECM_cmpt;
dxdt_ECM_cmpt = ECM_in - ECM_out;

// ============================================================
// Podocyte Integrity (AU, 1=intact, 0=lost)
// ============================================================
double Pod_loss_rate =
  kPod_loss
  + Pod_TGF   * (TGF_cmpt  - 1) * 0.01
  + Pod_AngII * (AngII_cmpt - 1) * 0.005 * (1 - RAAS_block)
  + Pod_ROS   * (ROS_cmpt   - 1) * 0.004
  + Pod_Pglo  * (AngII_cmpt * kPglo_AngII - Pglo_base) * 0.0005;
// Prevent below minimum
double Pod_eff = (Pod_cmpt > Pod_min) ? Pod_cmpt : Pod_min;
dxdt_Pod_cmpt = -Pod_loss_rate * Pod_eff;

// ============================================================
// UACR (mg/g)
// ============================================================
double Pod_damage = 1.0 - Pod_cmpt;   // 0=no damage → 1=full loss
double Pglo_curr  = Pglo_base + kPglo_AngII * (AngII_cmpt - 1) * (1 - RAAS_block);
double Pglo_excess = (Pglo_curr > Pglo_base) ? (Pglo_curr - Pglo_base) : 0;
double UACR_ss    = UACR_base * (1 + kUACR_Pod * Pod_damage / UACR_base
                                  + kUACR_Pglo * Pglo_excess / UACR_base);
dxdt_UACR_cmpt = kUACR_out * (UACR_ss - UACR_cmpt);

// ============================================================
// Tubular Integrity (AU, 1=healthy)
// ============================================================
double Tub_loss  = kTub_loss
                 + Tub_UACR    * (UACR_cmpt - 30) / 1000
                 + Tub_hypoxia * E_sglt2 * (-0.15)   // SGLT2i relieves hypoxia
                 + Tub_hypoxia * (BG_excess / 5.0) * 0.002;
double Tub_eff = (Tub_cmpt > 0.05) ? Tub_cmpt : 0.05;
dxdt_Tub_cmpt = -Tub_loss * Tub_eff + kTub_out * (1.0 - Tub_cmpt); // partial repair

// ============================================================
// Interstitial Fibrosis (AU)
// ============================================================
double Fib_in  = kFib_in * (Fib_TGF * TGF_cmpt + Fib_Tub * (1.0 - Tub_cmpt))
               * (1 - E_fine * 0.55);  // finerenone reduces fibrosis
double Fib_out = kFib_out * Fib_cmpt;
dxdt_Fib_cmpt = Fib_in - Fib_out;

// ============================================================
// eGFR (mL/min/1.73m2)
// ============================================================
// Driving forces for GFR decline
double GFR_decline = kGFR_Fib  * (Fib_cmpt - Fib_base)
                   + kGFR_ECM  * (ECM_cmpt  - 1.0)
                   + kGFR_Pglo * Pglo_excess;
// SGLT2i hemodynamic effect: acute dip then stabilization
double SGLT2_hemodynamic = E_sglt2 * SGLT2_GFR_dip * 0.05;  // small ongoing
double GFR_in  = 0;
double GFR_out = GFR_decline + SGLT2_hemodynamic;
double GFR_floor = (GFR_cmpt > GFR_min) ? 1.0 : 0.0;
dxdt_GFR_cmpt = -(GFR_out * GFR_floor);

$TABLE
// Derived PK
double Cp_acei_out  = CENT_acei  / Vd_acei  * 1e3;  // ng/mL
double Cp_arb_out   = CENT_arb   / Vd_arb   * 1e3;
double Cp_sglt2_out = CENT_sglt2 / Vd_sglt2 * 1e3;
double Cp_fine_out  = CENT_fine  / Vd_fine  * 1e3;

// Drug effect magnitudes
double Inh_ACE   = E_acei  * 100;  // %
double Inh_AT1R  = E_arb   * 100;
double Inh_SGLT2 = E_sglt2 * 100;
double Inh_MR    = E_fine  * 100;

// Clinical outputs
double eGFR   = GFR_cmpt;
double UACR_o = UACR_cmpt;
double HbA1c  = HbA1c_base - 0.5 * E_sglt2;   // simplified
double SBP    = SBP_base - 10*(E_acei + E_arb*0.8) - 5*E_sglt2;
double DBP    = DBP_base - 5*(E_acei + E_arb*0.8);

// CKD Stage classification
double CKD_Stage_num;
if      (eGFR >= 90)              CKD_Stage_num = 1;
else if (eGFR >= 60)              CKD_Stage_num = 2;
else if (eGFR >= 45)              CKD_Stage_num = 3.1;
else if (eGFR >= 30)              CKD_Stage_num = 3.2;
else if (eGFR >= 15)              CKD_Stage_num = 4;
else                              CKD_Stage_num = 5;

// % change from baseline
double pct_UACR  = (UACR_cmpt - 300) / 300 * 100;
double pct_GFR   = (GFR_cmpt  -  55) /  55 * 100;

$CAPTURE
Cp_acei_out Cp_arb_out Cp_sglt2_out Cp_fine_out
Inh_ACE Inh_AT1R Inh_SGLT2 Inh_MR
eGFR UACR_o HbA1c SBP DBP
TGF_cmpt ROS_cmpt ECM_cmpt Pod_cmpt Fib_cmpt AngII_cmpt
CKD_Stage_num pct_UACR pct_GFR BG AGE_cmpt
'

## -----------------------------------------------------------------------------
## Compile Model
## -----------------------------------------------------------------------------
mod <- mcode("DN_QSP", code)
cat("Model compiled successfully.\n")
cat("Compartments:", length(mod@cmtL), "\n")

## -----------------------------------------------------------------------------
## Dosing Regimens (helper)
## -----------------------------------------------------------------------------
build_events <- function(scenario = "S0",
                         end_days = 730,
                         acei_dose = 10,    # mg enalapril BID equiv
                         arb_dose  = 100,   # mg losartan QD
                         sglt2_dose = 25,   # mg empagliflozin QD
                         fine_dose  = 20) { # mg finerenone QD

  ev_null <- ev(time=0, amt=0, cmt=1)   # no treatment

  acei_ev  <- ev(amt = acei_dose,  cmt = 1, ii = 12, addl = end_days*2 - 1, time=0)
  arb_ev   <- ev(amt = arb_dose,   cmt = 3, ii = 24, addl = end_days   - 1, time=0)
  sglt2_ev <- ev(amt = sglt2_dose, cmt = 5, ii = 24, addl = end_days   - 1, time=0)
  fine_ev  <- ev(amt = fine_dose,  cmt = 7, ii = 24, addl = end_days   - 1, time=0)

  evts <- switch(scenario,
    "S0" = ev_null,
    "S1" = acei_ev,
    "S2" = arb_ev,
    "S3" = sglt2_ev,
    "S4" = ev_c(acei_ev, sglt2_ev),
    "S5" = ev_c(sglt2_ev, fine_ev),
    "S6" = ev_c(acei_ev, sglt2_ev, fine_ev),
    ev_null
  )
  evts
}

## -----------------------------------------------------------------------------
## Run All Treatment Scenarios — 2 years (730 days)
## -----------------------------------------------------------------------------
scenarios <- c("S0","S1","S2","S3","S4","S5","S6")
labels    <- c(
  "S0: No Treatment",
  "S1: ACEi (Enalapril)",
  "S2: ARB (Losartan)",
  "S3: SGLT2i (Empa)",
  "S4: ACEi + SGLT2i",
  "S5: SGLT2i + Finerenone",
  "S6: ACEi + SGLT2i + Finerenone"
)

END  <- 730   # days
STEP <- 7     # weekly output
results_list <- list()

for (i in seq_along(scenarios)) {
  sc  <- scenarios[i]
  evts <- build_events(sc, END)
  out  <- mrgsim(mod, ev = evts, end = END, delta = STEP, obsonly = TRUE)
  df   <- as.data.frame(out)
  df$Scenario  <- sc
  df$ScenLabel <- labels[i]
  df$time_yr   <- df$time / 365
  results_list[[sc]] <- df
}

results <- bind_rows(results_list)
results$ScenLabel <- factor(results$ScenLabel, levels = labels)

cat("Simulation complete:", nrow(results), "rows across", length(scenarios), "scenarios\n")

## -----------------------------------------------------------------------------
## Summary Table at Year 2
## -----------------------------------------------------------------------------
summary_tab <- results %>%
  filter(abs(time - 730) < 8) %>%
  group_by(Scenario, ScenLabel) %>%
  summarise(
    eGFR_yr2     = round(mean(eGFR), 1),
    dGFR         = round(mean(eGFR) - 55, 1),
    UACR_yr2     = round(mean(UACR_o), 0),
    pct_UACR_chg = round(mean(pct_UACR), 1),
    SBP_yr2      = round(mean(SBP), 1),
    TGF_yr2      = round(mean(TGF_cmpt), 2),
    Fib_yr2      = round(mean(Fib_cmpt), 3),
    Pod_yr2      = round(mean(Pod_cmpt), 3),
    CKD_Stage    = round(mean(CKD_Stage_num), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(eGFR_yr2))

cat("\n=== 2-Year Outcome Summary ===\n")
print(summary_tab, n = 20)

## -----------------------------------------------------------------------------
## Plotting
## -----------------------------------------------------------------------------

col_pal <- c(
  "S0: No Treatment"             = "#E53935",
  "S1: ACEi (Enalapril)"         = "#1E88E5",
  "S2: ARB (Losartan)"           = "#43A047",
  "S3: SGLT2i (Empa)"            = "#FB8C00",
  "S4: ACEi + SGLT2i"            = "#8E24AA",
  "S5: SGLT2i + Finerenone"      = "#00897B",
  "S6: ACEi + SGLT2i + Finerenone" = "#F4511E"
)

## Plot 1: eGFR over time
p1 <- ggplot(results, aes(time_yr, eGFR, color = ScenLabel)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(15, 30, 45, 60), linetype = "dashed",
             color = "grey60", alpha = 0.7) +
  annotate("text", x = 0.05, y = c(15,30,45,60)+1.5,
           label = c("G5","G4","G3b","G3a"), size = 3, color = "grey40") +
  scale_color_manual(values = col_pal) +
  labs(title = "eGFR Trajectory (Diabetic Nephropathy)",
       x = "Time (years)", y = "eGFR (mL/min/1.73m²)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Plot 2: UACR over time
p2 <- ggplot(results, aes(time_yr, UACR_o, color = ScenLabel)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(30, 300), linetype = "dashed", color = "grey60") +
  annotate("text", x = 0.05, y = c(50, 320),
           label = c("Microalbuminuria\n(30)", "Macroalbuminuria\n(300)"),
           size = 2.8, color = "grey40") +
  scale_color_manual(values = col_pal) +
  scale_y_log10() +
  labs(title = "UACR over Time",
       x = "Time (years)", y = "UACR (mg/g, log scale)",
       color = "Treatment") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

## Plot 3: TGF-β
p3 <- ggplot(results, aes(time_yr, TGF_cmpt, color = ScenLabel)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = col_pal) +
  labs(title = "TGF-β (Fibrosis Driver)",
       x = "Time (years)", y = "TGF-β (AU)",
       color = "Treatment") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

## Plot 4: Interstitial Fibrosis
p4 <- ggplot(results, aes(time_yr, Fib_cmpt, color = ScenLabel)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = col_pal) +
  labs(title = "Interstitial Fibrosis",
       x = "Time (years)", y = "Fibrosis (AU)",
       color = "Treatment") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

## Plot 5: Podocyte integrity
p5 <- ggplot(results, aes(time_yr, Pod_cmpt, color = ScenLabel)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = col_pal) +
  labs(title = "Podocyte Integrity",
       x = "Time (years)", y = "Podocyte Density (AU)",
       color = "Treatment") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

## Plot 6: Blood Pressure
p6 <- ggplot(results, aes(time_yr, SBP, color = ScenLabel)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "green4") +
  scale_color_manual(values = col_pal) +
  labs(title = "Systolic Blood Pressure",
       x = "Time (years)", y = "SBP (mmHg)",
       color = "Treatment") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

## Combine plots
combined <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title = "Diabetic Nephropathy QSP — Multi-Treatment Simulation",
    subtitle = "2-year trajectories across 7 treatment strategies",
    theme = theme(plot.title = element_text(size=16, face="bold"))
  )

print(combined)

## -----------------------------------------------------------------------------
## Drug PK profiles (first 48 hours)
## -----------------------------------------------------------------------------
run_pk <- function(drug = "sglt2", dose = 25, cmt = 5, ka = 1.5, vd = 74, cl = 9.2) {
  ev_pk <- ev(amt = dose, cmt = cmt, ii = 24, addl = 1, time = 0)
  out <- mrgsim(mod, ev = ev_pk, end = 48, delta = 0.25, obsonly = TRUE)
  as.data.frame(out)
}

pk_sglt2 <- run_pk("sglt2", 25, 5)
ggplot(pk_sglt2, aes(time, Cp_sglt2_out)) +
  geom_line(color = "#FB8C00", size = 1.2) +
  labs(title = "SGLT2i (Empagliflozin 25mg QD) — PK Profile",
       x = "Time (h)", y = "Plasma concentration (ng/mL)") +
  theme_bw(base_size = 12)

## -----------------------------------------------------------------------------
## GFR Slope Analysis (annualized mL/min/yr)
## -----------------------------------------------------------------------------
gfr_slopes <- results %>%
  group_by(Scenario, ScenLabel) %>%
  summarise(
    slope_yr1 = {
      d <- . %>% filter(time_yr <= 1)
      lm(eGFR ~ time_yr, data = d)$coefficients[2]
    },
    slope_yr2 = {
      lm(eGFR ~ time_yr, data = .)$coefficients[2]
    },
    .groups = "drop"
  )

cat("\n=== eGFR Slope (mL/min/year) ===\n")
print(gfr_slopes %>% select(ScenLabel, slope_yr1, slope_yr2) %>%
      mutate(across(c(slope_yr1,slope_yr2), round, 1)), n=10)

cat("\n\n")
cat("===========================================================\n")
cat("Diabetic Nephropathy QSP Model — Run complete\n")
cat("===========================================================\n")
cat("Model compartments: 19 ODEs\n")
cat("Treatment scenarios: 7 (S0–S6)\n")
cat("Simulation horizon: 2 years (730 days)\n")
cat("Key calibration targets:\n")
cat("  - Natural history GFR decline: ~3-4 mL/min/yr (Perkins 2003)\n")
cat("  - ACEi UACR reduction: ~30-35% (Lewis 1993 NEJM)\n")
cat("  - SGLT2i eGFR preservation: ~2 mL/min/yr (CREDENCE)\n")
cat("  - Finerenone UACR reduction: ~31% (FIDELIO-DKD NEJM 2020)\n")
cat("===========================================================\n")
